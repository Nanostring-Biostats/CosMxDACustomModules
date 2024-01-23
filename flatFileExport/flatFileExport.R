# Name: Flat File Export Custom Module
message("Custom Script Version: 1.1.1")

# Copyright 2023 NanoString Technologies, Inc.
# This software and any associated files are distributed pursuant to the NanoString AtoMx Spatial 
# Information Platform Software as a Service Agreement, available on the NanoString Technologies, Inc 
# website at www.nanostring.com, as updated.  All rights reserved.  No permission is granted to modify, 
# publish, distribute, sublicense or sell copies of this software.

# Description: Export Flat files used in Seurat's LoadNanostring function

# To get around v1.1 2-hour limit, suggested pipeline contains 3 flat file 
# export modules in any order:
# 1. export counts, metadata, & fovPositions
# 2. export transcripts
# 3. export polygons

# All default=TRUE files are required for Seurat LoadNanostring function

# User Defined Variables
# countMatrix - Generate count matrix file (default=TRUE)
# cellMetadata - Generate cell metadata file (default=TRUE)
# transcripts - Generate aligned transcripts file (default=TRUE)
# polygons - Generate cell boundaries (polygon) file (default=TRUE)
# fovPositions - Generate fov position file (default=FALSE)
# gzip - gzip files, does not affect LoadNanostring function
# studyName - Output Folder Name
# outPath - Destination S3 file path

# AWS Credentials:
# access_key    - Destination AWS access key
# secret_key    - Destination AWS secret key
# s3Region      - Destination AWS region
# session_token - Destination AWS session token, if needed

# Load Packages
library(nanopipeline)
library(stringr)
library(data.table)
library(tiledb)
library(reticulate)
library(rminiconda)
library(Matrix)

# Test User Variables
variableTest <- function(varName, varType, msg, required = TRUE){
  if(!varType %in% c("character", "logical", "numeric")){
    stop("varType must be \"character\", \"logical\", or \"numeric\"")
  }
  
  varMsg <- NULL
  
  typeClass <- switch(varType, 
                      character={"STRING"},
                      logical={"BOOL"},
                      numeric={"NUM"})
  
  if(required == TRUE){
    if(!exists(varName)){
      varMsg <- paste0("\n\"", varName, "\" ", "was not set in custom module creation")
    }
  }
  
  if(exists(varName)){
    if(varType == "logical"){
      if(is.null(get(varName))){
        assign(varName, FALSE, envir = .GlobalEnv)
      }
    }

    if(!is(get(varName), varType)){
      if(!is.null(get(varName)) | required == TRUE){
        varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" varType was not set as ", typeClass, 
                                        " in custom module creation"))
      }
    }
    
    if(required == TRUE){
      if(is.null(get(varName))){
        varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" was not set and is required"))
      }else{
        if(get(varName) == ""){
          varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" was not set and is required")) 
        }
      }
    }
  }
  
  if(!is.null(varMsg)){
    varMsg <- paste0(msg, varMsg)
  }else{
    varMsg <- msg
  }
  
  return(varMsg)
}

variableMsg <- NULL
variableMsg <- variableTest(varName = "countMatrix", varType = "logical", msg = variableMsg)
variableMsg <- variableTest(varName = "cellMetadata", varType = "logical", msg = variableMsg)
variableMsg <- variableTest(varName = "transcripts", varType = "logical", msg = variableMsg)
variableMsg <- variableTest(varName = "polygons", varType = "logical", msg = variableMsg)
variableMsg <- variableTest(varName = "fovPositions", varType = "logical", msg = variableMsg)
variableMsg <- variableTest(varName = "gzip", varType = "logical", msg = variableMsg)
variableMsg <- variableTest(varName = "studyName", varType = "character", msg = variableMsg)
variableMsg <- variableTest(varName = "outPath", varType = "character", msg = variableMsg)
variableMsg <- variableTest(varName = "access_key", varType = "character", msg = variableMsg)
variableMsg <- variableTest(varName = "secret_key", varType = "character", msg = variableMsg)
variableMsg <- variableTest(varName = "s3Region", varType = "character", msg = variableMsg)
variableMsg <- variableTest(varName = "session_token", varType = "character", required = FALSE, msg = variableMsg)

# stop on variable errors if applicable
if(!is.null(variableMsg)){
  stop(variableMsg)
}

# Module Code
dir.create("output")

fileType <- ifelse(gzip == TRUE, yes = ".csv.gz", no = ".csv")

# change AWS credentials to user environment
if(exists("session_token")){
  session_token <- paste("aws_session_token =", session_token)
}else{
  session_token <- NULL
}

dir.create("~/.aws")

fileConn<-file("~/.aws/credentials", )
writeLines(c("[export]",paste("aws_access_key_id =", access_key), 
             paste("aws_secret_access_key =", secret_key), session_token), fileConn)
close(fileConn)

############################### Helper functions ###############################
getAnnots <- function(){
  # read in annotation file
  # get input folder location and slide names
  all_files <- system2(command = "aws",arg = c("s3","ls", paste0(studyDirectory, "logs/")), stdout = TRUE)
  studyCreation <- all_files[grep("study_creation_", all_files)]
  studyCreation <- strsplit(studyCreation, " ")[[1]]
  studyCreation <- studyCreation[length(studyCreation)]

  studyCreation <- paste0(paste0(studyDirectory, "logs/"), studyCreation)

  annots <- system2(command = "aws",arg = c("s3","ls", studyCreation), stdout = TRUE)
  annots <- annots[grep(".csv", annots)]
  annots <- annots[grep("Resources", annots, invert = TRUE)]
  annots <- strsplit(annots, " ")[[1]]
  annots <- annots[length(annots)]
  annots <- paste0(studyCreation, annots)

  bucket <- strsplit(annots,"/")[[1]][3]
  file_key <- strsplit(annots, paste0(bucket,"/"))[[1]][2]
  system2(command = "aws", arg = c("s3api", "get-object","--bucket",bucket,"--key", file_key,paste0("/tmp/tmp.csv")), stdout = TRUE)
  annots <- read.delim("/tmp/tmp.csv", header = TRUE, sep = ",")
  file.remove("/tmp/tmp.csv")
  rm(bucket)
  rm(file_key)
  
  return(annots)
}

getPxToUmValue <- function(){
  annots <- getAnnots()
  
  px_conversion_file <- fs::path(annots$folder_path,
                                 annots$slidefolders,
                                 "RunSummary")[1L]
  
  px_conversion_file <- paste0(gsub(pattern = "s3:/", replacement = "s3://", px_conversion_file), "/")
  all_files <- system2(command = "aws",arg = c("s3","ls", px_conversion_file), stdout = TRUE)
  
  temp_file_vec <- all_files[grep("ExptConfig.txt", all_files)]
  temp_file_vec <- strsplit(temp_file_vec, " ")[[1L]]
  temp_file_vec <- temp_file_vec[length(temp_file_vec)]
  
  if(grepl("alpha|dash", tolower(temp_file_vec))){
    stop("Instrument version must be BETA to use this custom module")
  }
  
  px_conversion_file <- paste0(px_conversion_file, temp_file_vec)
  
  bucket <- strsplit(px_conversion_file,"/")[[1L]][3L]
  file_key <- strsplit(px_conversion_file, paste0(bucket,"/"))[[1L]][2L]
  system2(command = "aws", arg = c("s3api", "get-object","--bucket",bucket,"--key", file_key,paste0("/tmp/tmp.csv")), stdout = TRUE)
  px_conversion_file <- read.delim("/tmp/tmp.csv")
  file.remove("/tmp/tmp.csv")
  rm(bucket)
  rm(file_key)
  
  px_line <- grep(pattern = "ImPixel_nm", px_conversion_file[[1L]])
  pixel_conversion_value <- px_conversion_file[[1L]][px_line]
  pixel_conversion_value <- as.numeric(stringr::str_split(string = pixel_conversion_value, 
                                                          pattern = ": ", simplify = TRUE)[,2L])
  pixel_conversion_value <- (pixel_conversion_value/1e3) # nm to um
  
  return(pixel_conversion_value)
}

align_latest.fovs <- function(fov_position,
                              spatial_columns = c("x", "y"), 
                              fov_column = "FOV", 
                              prefixFOV = "FOV_", 
                              slide_column = NULL){
  # data.frame
  locs <- as.data.frame(fov_position)
  
  fov_spacing_value <- 1

  cols <- colnames(locs)[!colnames(locs) %in% c("locfile", "Run_Tissue_name")]
  for(i in cols){
    locs[[i]] <- as.numeric(locs[[i]])
  }
  
  # change the address to be relative by shifting xy in each slot/slide
  # so we can set the arg relative_to_slot to be `TRUE` later.
  loc_annot_list <- lapply(seq_along(unique(locs$Slide)), function(i){
    # grab all x y coord for the ith slot/slide
    loc_annot_tmp <- locs[which(locs$Slide==i), ]
    # shift the x coord to the left
    loc_annot_tmp$X_mm <- (loc_annot_tmp$X_mm - min(loc_annot_tmp$X_mm))*fov_spacing_value
    # shift the y coord to the bottom
    loc_annot_tmp$Y_mm <- (loc_annot_tmp$Y_mm - min(loc_annot_tmp$Y_mm))*fov_spacing_value
    return(loc_annot_tmp)
  })
  
  locs <- do.call(rbind, loc_annot_list)
  
  # check whether spatial_columns are correctly specified 
  if ( any(!(c(fov_column, spatial_columns) %in% colnames(locs))) ){
    stop(sprintf("%s, %s are not found in provided fov_position.", fov_column,
                 paste0(spatial_columns, collapse = ', ')))
  }
  
  # get the informative columns
  # strip off prefix for FOV name
  locs[[paste0('ori_', fov_column)]] <- locs[[fov_column]]
  locs[[fov_column]] <- gsub(prefixFOV, "", locs[[fov_column]])
  
  # create numeric values for fov and slide
  fov_converter <- seq_len(length(unique(locs[[fov_column]])))
  names(fov_converter) <- unique(locs[[fov_column]])
  locs[[fov_column]] <- fov_converter[locs[[fov_column]]]
  if(!is.null(slide_column)){
    locs[[paste0('ori_', slide_column)]] <- locs[[slide_column]]
    slide_converter <- seq_len(length(unique(locs[[slide_column]])))
    names(slide_converter) <- unique(locs[[slide_column]])
    locs[[slide_column]] <- slide_converter[locs[[slide_column]]]
  }
  
  locs_info <- data.frame(slide= locs[[slide_column]],
                          fov = locs[[fov_column]],
                          x = locs[[spatial_columns[1L]]],
                          y = locs[[spatial_columns[2L]]])
  
  locs <- locs_info[c("x","y","fov","slide")]
  
  locs$x <- -locs$x + max(locs$x)
  locs$y <- -locs$y + max(locs$y)
  
  return(locs)
}

usableCores <- function(percentCores = 0.9, minNotUsedCores = 2){
  
  if(percentCores > 1 & percentCores <= 0){
    stop("percentCores is not a valid number, must be between 0-1")
  }
  
  nCores <- parallel::detectCores()
  
  if(nCores <= minNotUsedCores){
    stop("minNotUsedCores must be fewer than available cores")
  }
  
  usableCores <- min(floor(nCores*percentCores), nCores-minNotUsedCores)
  
  return(usableCores)
}

################################ Cell Metadata #################################
message("Generating Cell Metadata file")

pxToUm <- getPxToUmValue()
pixel_size <- um_to_mm(pxToUm)
annots <- getAnnots()

obs <- study$somas$RNA$obs$to_dataframe()
colnames(obs)[grep("slide_ID_numeric", colnames(obs))] <- "slide_ID"
colnames(obs)[grep("orig", colnames(obs))] <- "cell"
obs$cell <- obs$cell_ID 
obs$cell_ID <- as.numeric(str_split(obs$cell_ID, "_", simplify = TRUE)[,4])
obs <- obs[order(obs$slide_ID, obs$fov, obs$cell_ID),]

obs$x_slide_mm <- mm_to_pixel(obs$x_slide_mm, pxToUm)
obs$y_slide_mm <- mm_to_pixel(obs$y_slide_mm, pxToUm)

colnames(obs)[grep("x_FOV_px", colnames(obs))] <- "CenterX_local_px"
colnames(obs)[grep("y_FOV_px", colnames(obs))] <- "CenterY_local_px"
colnames(obs)[grep("x_slide_mm", colnames(obs))] <- "CenterX_global_px"
colnames(obs)[grep("y_slide_mm", colnames(obs))] <- "CenterY_global_px"

for(slide in unique(obs$Run_Tissue_name)){
  slideName <- gsub("\\W", "", slide)
  outDir <- paste0("output/",slideName, "/")
  dir.create(outDir)
  
  if(cellMetadata == TRUE){
    data.table::fwrite(file = paste0(outDir, slideName,"_metadata_file", fileType), 
                       x = obs[obs$Run_Tissue_name == slide,], quote = FALSE, 
                       row.names = FALSE, sep = ",")
  }
}

slideFOVs <- unique(obs[,c("slide_ID", "fov", "Run_Tissue_name")])

################################# Cell Counts ##################################
if(countMatrix == TRUE){
  message("Generating Cell Counts file")
  
  RNAcounts <- Matrix::t(study$somas$RNA$X$members$counts$to_matrix(batch_mode = TRUE))
  RNAcounts <- nanopipeline:::impute_missing_cells(RNAcounts, obs$cell, 0)
  
  NEGcounts <- Matrix::t(study$somas$negprobes$X$members$counts$to_matrix(batch_mode = TRUE))
  NEGcounts <- nanopipeline:::impute_missing_cells(NEGcounts, obs$cell, 0)
  
  counts <- cbind(RNAcounts[obs$cell,], NEGcounts[obs$cell,])
  
  rm(RNAcounts, NEGcounts)
  gc()
  
  if("falsecode" %in% names(study$somas)){
    FALSEcounts <- Matrix::t(study$somas$falsecode$X$members$counts$to_matrix(batch_mode = TRUE))
    FALSEcounts <- nanopipeline:::impute_missing_cells(FALSEcounts, obs$cell, 0)
    
    counts <- cbind(counts, FALSEcounts[obs$cell,])
    
    rm(FALSEcounts)
    gc()
  }
  
  for(idx in seq_len(nrow(slideFOVs))){
    slideName <- slideFOVs$Run_Tissue_name[idx]
    cells <- which(obs$Run_Tissue_name == slideName & 
                     obs$fov == slideFOVs$fov[idx])
    
    slideName <- gsub("\\W", "", slideName)
    
    countsFOV <- cbind(obs[cells,c("fov","cell_ID")], counts[cells,,drop=FALSE])
    
    data.table::fwrite(file = paste0("output/",slideName, "/", slideName, "_exprMat_file", fileType), 
                       x = countsFOV, quote = FALSE, 
                       append = ifelse(slideFOVs$fov[idx] == 1, FALSE, TRUE),
                       row.names = FALSE, sep = ",")
  }
  
  rm(countsFOV,counts, cells)
  gc()
}

################################ FOV Positions #################################
message("Generating FOV Positions file")

latest.fovs <- study$somas$RNA$obsm$members$latest.fovs$to_matrix(batch_mode = TRUE)
latest.fovs <- as.data.frame(latest.fovs)
latest.fovs <- latest.fovs[order(latest.fovs$Slide, as.integer(latest.fovs$FOV)),]

locs <- align_latest.fovs(fov_position = latest.fovs,
                          slide_column = "Slide",
                          spatial_columns = c("X_mm", "Y_mm"))

colnames(locs) <- c("X_mm", "Y_mm", "FOV", "Slide")

locs <- locs[order(locs$Slide, as.integer(locs$FOV)),
             c("Slide", "FOV", "X_mm", "Y_mm")]

if(fovPositions == TRUE){
  for(slide in unique(obs$slide_ID)){
    slideName <- unique(obs$Run_Tissue_name[obs$slide_ID == slide])
    slideName <- gsub("\\W", "", slideName)
    
    data.table::fwrite(file = paste0("output/",slideName, "/", slideName, "_fov_positions_file", fileType), 
                       x = locs[locs$Slide == slide,], quote = FALSE, 
                       row.names = FALSE, sep = ",")
  }
}

rm(latest.fovs)
gc()

############################### Cell Transcripts ###############################
if(transcripts == TRUE){
  if(any(annots$assay_type == "RNA")){
    message("Generating Cell Transcripts file")
    
    coords <- c("x_FOV_px", "y_FOV_px","x_slide_mm", "y_slide_mm")
  
    slideFOVsAll <- unique(obs[,c("Run_Tissue_name", "slide_ID", "fov")])
    binSize <- 40
    
    for(slideName in unique(slideFOVs$Run_Tissue_name)){
      slideFOVs <- slideFOVsAll[slideFOVsAll$Run_Tissue_name == slideName,]
      numFOVs <- nrow(slideFOVs)
      nbins <- ceiling(numFOVs / binSize)
      
      startOfBin <- sequence(from = 1, by = binSize, nvec = nbins)
      endOfBin <- sequence(from = binSize, by = binSize, nvec = nbins)
      
      if(nbins == 1){
        endOfBin[nbins] <- numFOVs
      }else if(numFOVs - endOfBin[nbins - 1] < binSize/4){
        endOfBin[nbins - 1] <- numFOVs
        endOfBin <- endOfBin[-nbins]
        startOfBin <- startOfBin[-nbins]
      }else{
        endOfBin[nbins] <- numFOVs
      }
      
      bins <- cbind(slideFOVs[startOfBin,c("Run_Tissue_name", "slide_ID")], 
                    lower.fov=slideFOVs[startOfBin, "fov"], 
                    upper.fov=slideFOVs[endOfBin,"fov"])
      
      NULLList <- lapply(seq_len(nrow(bins)), function(idx){
        qc <- tiledb::tiledb_query_condition_combine(
          tiledb::tiledb_query_condition_init("slideID", bins$slide_ID[idx], "INT32", "EQ"),
          tiledb::tiledb_query_condition_combine(
            tiledb::tiledb_query_condition_init("fov", bins$lower.fov[idx], "INT32", "GE"), 
            tiledb::tiledb_query_condition_init("fov", bins$upper.fov[idx], "INT32", "LE"),
            "AND"), "AND")
        
        tCoords <- rbind(tiledb_array(study$somas$RNA$obsm$members$transcriptCoords$uri, 
                                      query_condition=qc, return_as="data.frame")[], 
                         tiledb_array(study$somas$negprobes$obsm$members$transcriptCoords$uri, 
                                      query_condition=qc, return_as="data.frame")[])
        
        
        if("falsecode" %in% names(study$somas)){
          tCoords <- rbind(tCoords, tiledb_array(study$somas$falsecode$obsm$members$transcriptCoords$uri, 
                                                 query_condition=qc, return_as="data.frame")[])
        }
        
        if(!"y_slide_mm" %in% study$somas$RNA$obsm$members$transcriptCoords$attrnames()){
          tCoords$x_slide_mm <- 0
          tCoords$y_slide_mm <- 0
          for(fov in bins$lower.fov[idx]:bins$upper.fov[idx]){
            fovTx <- which(tCoords$fov == fov)
            
            abs_locs <- locs[locs$Slide == bins$slide_ID[1L] & 
                               locs$FOV == fov,]
            
            tCoords[fovTx, "x_slide_mm"] <- (tCoords[fovTx, "x_FOV_px"])*pixel_size + abs_locs[["X_mm"]]
            tCoords[fovTx, "y_slide_mm"] <- (-tCoords[fovTx, "y_FOV_px"])*pixel_size + abs_locs[["Y_mm"]]
          }
        }
        
        tCoords <- tCoords[order(tCoords$slideID, tCoords$fov, tCoords$CellId), 
                           c("fov", "CellId", "cell_id", coords, "z_FOV_slice", 
                             "target", "CellComp")]
        
        # Seurat export needs to be in pixels not mm will rename column before exporting
        tCoords$x_slide_mm <- mm_to_pixel(tCoords$x_slide_mm, pxToUm)
        tCoords$y_slide_mm <- mm_to_pixel(tCoords$y_slide_mm, pxToUm)
        
        colnames(tCoords)[grep("x_FOV_px", colnames(tCoords))] <- "x_local_px"
        colnames(tCoords)[grep("y_FOV_px", colnames(tCoords))] <- "y_local_px"
        colnames(tCoords)[grep("x_slide_mm", colnames(tCoords))] <- "x_global_px"
        colnames(tCoords)[grep("y_slide_mm", colnames(tCoords))] <- "y_global_px"
        colnames(tCoords)[grep("z_FOV_slice", colnames(tCoords))] <- "z"
        colnames(tCoords)[grep("cell_id", colnames(tCoords))] <- "cell"
        colnames(tCoords)[grep("CellId", colnames(tCoords))] <- "cell_ID"
        
        slideName <- gsub("\\W", "",slideFOVs$Run_Tissue_name[idx])
        
        data.table::fwrite(file = paste0("output/",slideName, "/", slideName,"_tx_file", fileType), 
                           x = tCoords, quote = FALSE, sep = ",",
                           row.names = FALSE, append = ifelse(idx == 1, FALSE, TRUE), 
                           col.names = ifelse(idx == 1, TRUE, FALSE))
      })
    }
    
    rm(NULLList)
    gc()
  }else{
    message("Protein data does not contain transcripts, no transcript file will be produced.")
  }
}

################################## Polygons ####################################
if(polygons == TRUE){
  message("Generating Cell Polygon file")
  
  if(grepl("cn", s3Region)){
    rminiconda::rminiconda_pip_install(pkg_name = c("opencv-python==4.8.0.74"), name = "nanopipeline_python",
                                       "-i https://mirrors.aliyun.com/pypi/simple/")
  }else{
    rminiconda::rminiconda_pip_install(pkg_name = c("opencv-python==4.8.0.74"), name = "nanopipeline_python")
  }
  
  # Module Code
  pyscript <- 
    "import cv2 as cv
import numpy as np

def GetPolygon(image, cell_ids):

  im = cv.imread(image, cv.IMREAD_UNCHANGED)
  assert im is not None, 'file could not be read, check with os.path.exists()'
  
  boundaries = []
  
  for cellID in cell_ids:
    cell = cv.inRange(im, cellID,cellID)
  
    contours, heirarchy = cv.findContours(cell, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
  
    if len(contours) > 0:
      hull = []
  
      # calculate points for each contour
      for i in range(len(contours)):
        # creating convex hull object for each contour
        hull.append(cv.convexHull(contours[i], False))
  
      if len(hull) > 1:
        lengths = []
        for i in range(len(hull)):
          lengths.append(len(hull[i]))
        
        a = np.asarray(hull[lengths.index(max(lengths))]) 
        if(a.shape[0] > 3):
          a = a.reshape(a.shape[0],2)
          boundaries.append(a)
        else:
          boundaries.append(np.nan)
      else:
        a = np.asarray(hull)
        if(a.shape[1] > 3):
          a = a.reshape(a.shape[1],2)
          boundaries.append(a)
        else:
          boundaries.append(np.nan)
    else:
      boundaries.append(np.nan)
  
  return(boundaries)"
  writeLines(pyscript, "GetPolygon.py")
  
  py <- reticulate::py_run_file("GetPolygon.py")
  
  cellCounts <- aggregate(obs$cell_ID, list(obs$slide_ID, obs$fov), max)
  colnames(cellCounts) <- c("slide_ID_numeric", "fov", "maxCellID")
  
  for(slideID in seq_len(nrow(annots))){
    CellStatsDir <- paste0(annots$folder_path[slideID],
                           annots$slidefolders[slideID], 
                           "/CellStatsDir/")
    FOVs <- system2(command = "aws",arg = c("s3","ls", CellStatsDir), stdout = TRUE)
    FOVs <- FOVs[grep("FOV", FOVs)]
    if(length(FOVs) == 0){
      message("Study has a different file structure than expected, no polygons will be generated")
    }else{
      boundaryList <- parallel::mclapply(FOVs, function(FOV){
        FOV <- strsplit(FOV, " ")[[1]]
        FOV <- FOV[length(FOV)]
        FOV_ID <- as.numeric(gsub("/", "", gsub("FOV", "", FOV)))
        FOV <- paste0(CellStatsDir, FOV)
        
        CellLabel <- system2(command = "aws",arg = c("s3","ls", FOV), stdout = TRUE)
        CellLabel <- CellLabel[grep("CellLabel", CellLabel)]
        CellLabel <- strsplit(CellLabel, " ")[[1]]
        CellLabel <- paste0(FOV, CellLabel[length(CellLabel)])
        
        numCells <- cellCounts$maxCellID[cellCounts$slide_ID_numeric == slideID &
                                           cellCounts$fov == FOV_ID]
        
        abs_locs <- locs[locs$Slide == slideID & locs$FOV == FOV_ID,]
        
        if(nrow(abs_locs) == 0){
          slide <- annots$Run_Tissue_name[slideID]
          abs_locs <- locs[locs$Slide == slide & locs$FOV == FOV_ID,]
        }
        
        if(length(numCells) > 0){
          if(numCells > 1){
            system2(command = "aws",arg = c("s3","cp", CellLabel, "."), stdout = FALSE)
            
            boundList <- py$GetPolygon(image = basename(CellLabel), cell_ids= seq_len(numCells))
            names(boundList) <- 1:length(boundList)
            
            boundaries <- sapply(seq_len(length(boundList)), function(i){
              if(!any(is.nan(boundList[[i]]))){
                cbind(slideID = slideID, fov = FOV_ID, 
                      CellID=i, vertexID=seq_len(nrow(boundList[[i]])), boundList[[i]],
                      x_slide_mm = 0, y_slide_mm = 0)
              }
            },simplify = FALSE)
            
            unlink(basename(CellLabel))
            
            boundaries <- do.call("rbind.data.frame", boundaries)
            
            colnames(boundaries)[5:6] <- c("x_FOV_px", "y_FOV_px")
            
            boundaries[, "x_slide_mm"] <- (boundaries[, "x_FOV_px"])*pixel_size + abs_locs[["X_mm"]]
            boundaries[, "y_slide_mm"] <- (-boundaries[, "y_FOV_px"])*pixel_size + abs_locs[["Y_mm"]]
            
            boundaries$cell_id <- paste("c", boundaries$slideID, boundaries$fov, boundaries$CellID, sep = "_")
            
            # Seurat export needs to be in pixels not mm will rename column before exporting
            boundaries$x_slide_mm <- mm_to_pixel(boundaries$x_slide_mm, pxToUm)
            boundaries$y_slide_mm <- mm_to_pixel(boundaries$y_slide_mm, pxToUm)
            
            return(boundaries)
          }else{
            return(NULL)
          }
        }else{
          return(NULL)
        }
        
      }, mc.cores = usableCores())
      
      boundaries <- do.call(rbind.data.frame, boundaryList[which(unlist(lapply(boundaryList,function(x){!is.null(x)})))])
      rm(boundaryList)
      
      boundaries <- boundaries[,c("fov", "CellID", "cell_id",
                                  "x_FOV_px", "y_FOV_px",
                                  "x_slide_mm", "y_slide_mm")]
      
      colnames(boundaries)[grep("CellID", colnames(boundaries))] <- "cellID"
      colnames(boundaries)[grep("cell_id", colnames(boundaries))] <- "cell"
      colnames(boundaries)[grep("x_FOV_px", colnames(boundaries))] <- "x_local_px"
      colnames(boundaries)[grep("y_FOV_px", colnames(boundaries))] <- "y_local_px"
      colnames(boundaries)[grep("x_slide_mm", colnames(boundaries))] <- "x_global_px"
      colnames(boundaries)[grep("y_slide_mm", colnames(boundaries))] <- "y_global_px"
      
      slideName <- gsub("\\W", "", annots$Run_Tissue_name[slideID])
      
      data.table::fwrite(file = paste0("output/",slideName, "/", slideName,"-polygons", fileType), 
                         x = boundaries, quote = FALSE, row.names = FALSE)
    }
  }
}

################################## Move to S3 ##################################
system2(command = "aws", args = c("s3", "mv", "output/", paste0(outPath, studyName, "/"),
                                  "--recursive", "--profile", "export"), stdout = FALSE)

unlink("output/", recursive = TRUE)
