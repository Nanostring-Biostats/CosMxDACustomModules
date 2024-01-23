# CosMxDA Export Custom Module
message("Custom Script Version: 1.2.3")

# Copyright 2023 NanoString Technologies, Inc.
# This software and any associated files are distributed pursuant to the NanoString AtoMx Spatial 
# Information Platform Software as a Service Agreement, available on the NanoString Technologies, Inc 
# website at www.nanostring.com, as updated.  All rights reserved.  No permission is granted to modify, 
# publish, distribute, sublicense or sell copies of this software.

# Export all files associated with current study: Seurat, TileDB array, and/or raw files. 
# For booleans: rawFiles, tiledbArray, & SeuratObject determine what type of data is exported. 
# exportFOVImages, spotFiles, FullSeuratObject, & transcripts determine what is included in the data type export. 
# To export spotFiles, rawFiles & spotFiles MUST be checked.

# User defined variables
# studyName     - Output Folder Name
# outPath       - Destination S3 file path
# access_key    - Destination AWS access key
# secret_key    - Destination AWS secret key
# s3Region      - Destination AWS region
# session_token - Destination AWS session token, if needed

# Advanced Boolean Options
# rawFiles         - Export Raw Files
# SeuratObject     - Export a Seurat Object
# FullSeuratObject - Seurat object contains all previous module output
# transcripts      - Export Seurat contains transcript coordinates (large data)
# tiledbArray      - Export TileDB array
# spotFiles        - Export spot files to redo Target Decoding (large data)
# exportFOVImages  - Export all FOV Images (large data)


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
variableMsg <- variableTest(varName = "studyName", varType = "character", msg = variableMsg)
variableMsg <- variableTest(varName = "outPath", varType = "character", msg = variableMsg)
variableMsg <- variableTest(varName = "access_key", varType = "character", msg = variableMsg)
variableMsg <- variableTest(varName = "secret_key", varType = "character", msg = variableMsg)
variableMsg <- variableTest(varName = "s3Region", varType = "character", msg = variableMsg)
variableMsg <- variableTest(varName = "session_token", varType = "character", required = FALSE, msg = variableMsg)
variableMsg <- variableTest(varName = "rawFiles", varType = "logical", msg = variableMsg)
variableMsg <- variableTest(varName = "SeuratObject", varType = "logical", msg = variableMsg)
variableMsg <- variableTest(varName = "FullSeuratObject", varType = "logical", msg = variableMsg)
variableMsg <- variableTest(varName = "transcripts", varType = "logical", msg = variableMsg)
variableMsg <- variableTest(varName = "tiledbArray", varType = "logical", msg = variableMsg)
variableMsg <- variableTest(varName = "spotFiles", varType = "logical", msg = variableMsg)
variableMsg <- variableTest(varName = "exportFOVImages", varType = "logical", msg = variableMsg)

# stop on variable errors if applicable
if(!is.null(variableMsg)){
  stop(variableMsg)
}

library(nanopipeline)
library(parallel)

# nanopipeline v1.2 functions
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

readInTranscriptCoords <- function(tCoords){
  if(!is(tCoords, "AnnotationMatrix")){
    stop("tCoords must be a AnnotationMatrix")
  }
  
  if(!grepl("transcriptCoords", tCoords$uri)){
    stop("Function is only meant for transcript coordinates")
  }
  
  message(paste("Reading AnnotationMatrix into memory from", tCoords$uri))
  tCoords <- tiledb::tiledb_array(tCoords$uri, return_as="data.frame")[]
  
  tCoords <- tCoords[,which(colnames(tCoords) != "__tiledb_rows")]
  
  tCoords$x_FOV_px <- as.integer(tCoords$x_FOV_px)
  tCoords$y_FOV_px <- as.integer(tCoords$y_FOV_px)
  
  if(any(c("x_slide_mm", "y_slide_mm") %in% colnames(tCoords))){
    tCoords$x_slide_mm <- as.numeric(tCoords$x_slide_mm)
    tCoords$y_slide_mm <- as.numeric(tCoords$y_slide_mm)
  }else{
    tCoords$z_FOV_slice <- as.integer(tCoords$z_FOV_slice)
  }
  
  tCoords$fov <- as.integer(tCoords$fov)
  tCoords$CellId <- as.integer(tCoords$CellId)
  
  return(tCoords)
}

toFullSeuratObject <- function(tiledbsc_dataset, transcripts = TRUE){
  if(!is(tiledbsc_dataset, "SOMACollection")){
    stop("tiledbsc_dataset must be a SOMACollection")
  }
  
  sobject <- suppressWarnings(tiledbsc_dataset$to_seurat(batch_mode = TRUE))
  Seurat::DefaultAssay(sobject) <- "RNA"
  
  obs <- tiledbsc_dataset$somas$RNA$obs$to_dataframe()
  sobject@meta.data <- obs[match(colnames(sobject), rownames(obs)),]
  
  if(any(grepl("RNA_normalized", names(sobject@assays)))){
    assayNames <- names(sobject@assays)[grepl("RNA_normalized", names(sobject@assays))]
    
    # removing raw data from normalized assays to avoid duplications
    for(i in assayNames){
      sobject[[i]]@counts <- matrix(data = 0, nrow = 0, ncol = 0)
    }
  }
  
  if("transcriptCoords" %in% names(tiledbsc_dataset$somas$RNA$obsm$members) & transcripts == TRUE){
    sobject@misc$transcriptCoords <- NULL
    for(i in names(tiledbsc_dataset$somas)){
      if("transcriptCoords" %in% names(tiledbsc_dataset$somas[[i]]$obsm$members)){
        sobject@misc$transcriptCoords <- rbind(sobject@misc$transcriptCoords,
                                               readInTranscriptCoords(tiledbsc_dataset$somas[[i]]$obsm$members$transcriptCoords))
      }
    }
  }
  
  if("latest.fovs" %in% names(tiledbsc_dataset$somas$RNA$obsm$members)){
    sobject@misc$latest.fovs <- as.data.frame(tiledbsc_dataset$somas$RNA$obsm$members$latest.fovs$to_matrix())
  }
  
  ### OBSM ###
  
  slots <- names(tiledbsc_dataset$somas$RNA$obsm$members)
  slots <- slots[!grepl("dimreduction|latest.fovs|transcriptCoords", slots)]
  
  for(i in slots){
    sobject@misc[[i]] <- tryCatch(tiledbsc_dataset$somas$RNA$obsm$members[[i]]$to_dataframe(), 
                                  error = function(e){tiledbsc_dataset$somas$RNA$obsm$members[[i]]$to_matrix(batch_mode = TRUE)})
    
    if(rownames(sobject@misc[[i]])[1L] %in% colnames(sobject)){
      sobject@misc[[i]] <- sobject@misc[[i]][colnames(sobject),]
    }
  }
  
  for(i in names(sobject@reductions)){
    sobject@reductions[[i]]@assay.used <- "RNA"
  }
  
  ### OBSP ###
  
  for(i in names(sobject@graphs)){
    sobject@graphs[[i]] <- NULL
  }
  
  for (i in tiledbsc_dataset$somas$RNA$obsp$members) {
    if(i$exists()){
      name <- i$uri
      name <- strsplit(name, "/")[[1L]]
      
      slot <- name[length(name)]
      
      mtx <- tiledbsc::AnnotationPairwiseMatrix$new(i$uri)$to_sparse_matrix()
      missing <- colnames(sobject)[which(!colnames(sobject) %in% colnames(mtx))]
      
      if(length(missing) > 0){
        missingCols <- matrix(nrow = nrow(mtx), ncol = length(missing), 
                              dimnames = list(rownames(mtx), missing), data = 0)
        mtx <- cbind(mtx, missingCols)
        
        missingRows <- matrix(nrow =length(missing), ncol = ncol(mtx), 
                              dimnames = list(missing, colnames(mtx)), data = 0)
        mtx <- rbind(mtx, missingRows)
      }
      
      mtx <- mtx[colnames(sobject), ]
      mtx <- mtx[, colnames(sobject)]
      
      sobject@graphs[[slot]] <- mtx
    }
  }
  
  ### VARM ### 
  slots <- names(tiledbsc_dataset$somas$RNA$varm$members)
  slots <- slots[!grepl("dimreduction", slots)]
  
  for(i in slots){
    sobject[["RNA"]]@misc[[i]] <- tiledbsc_dataset$somas$RNA$varm$members[[i]]$to_matrix(batch_mode = TRUE)
  }
  
  ### VARP ###
  
  slots <- names(tiledbsc_dataset$somas$RNA$varp$members)
  
  for(i in slots){
    if(grepl("log_likelyhood", i)){
      sobject[["RNA"]]@misc[[i]] <- tiledbsc_dataset$somas$RNA$varp$members[[i]]$to_seurat_graph()
    }else{
      sobject[["RNA"]]@misc[[i]] <- tiledbsc_dataset$somas$RNA$varp$members[[i]]$to_matrix()
      
      # sorts target names, ensuring that any target that starts with a number is not in the first slot.
      if(sort(colnames(sobject[["RNA"]]@misc[[i]]), decreasing = TRUE)[1L] %in% rownames(sobject)){
        editedColNames <- gsub(pattern = "\\W" , ".", rownames(sobject))
        if(any(grepl("^[[:digit:]]+", editedColNames))){
          num <- grep("^[[:digit:]]+", editedColNames)
          editedColNames[num] <- paste0("X", editedColNames[num])
        }
        
        sobject[["RNA"]]@misc[[i]] <- sobject[["RNA"]]@misc[[i]][, editedColNames]
        colnames(sobject[["RNA"]]@misc[[i]]) <- rownames(sobject)
      }
    }
  }
  
  ### UNS ###
  
  for(i in names(tiledbsc_dataset$somas$RNA$uns$members)){
    if(tiledbsc_dataset$somas$RNA$uns$members[[i]]$class() == "TileDBGroup"){
      if(grepl("cellproximity", i)){
        if(any(grepl(tiledbsc_dataset$somas$RNA$uns$members[[i]]$list_objects()$URI, 
                     pattern = "which"))){
          # v1.1
          objectURIs <- tiledbsc_dataset$somas$RNA$uns$members[[i]]$list_object_uris()
        }else{
          # v1.2
          objectURIs <- unlist(lapply(tiledbsc_dataset$somas$RNA$uns$members[[i]]$members, 
                                      function(self) self$list_object_uris()))
        }
      }else{
        objectURIs <- tiledbsc_dataset$somas$RNA$uns$members[[i]]$list_object_uris()
      }
      
      for(j in 1:length(objectURIs)){
        name <- paste0(i, "_", names(objectURIs)[j])
        
        sobject@misc[[name]] <- nanopipeline:::from_tdb(objectURIs[j])
      }
    }else{
      sobject@misc[[i]] <- nanopipeline:::from_tdb(tiledbsc_dataset$somas$RNA$uns$members[[i]]$uri)
    }
  }
  
  return(sobject)
}

# Variables for future filtering capability, not currently used 
FOVs <- NULL
flowcell <- NULL 

# Check the nested booleans
if(FullSeuratObject == TRUE){SeuratObject <- TRUE}
if(transcripts == TRUE){SeuratObject <- TRUE}
if(spotFiles == TRUE){rawFiles <- TRUE}
if(exportFOVImages == TRUE){rawFiles <- TRUE}

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

# read in annotation file 
# get input folder location and slide names
all_files <- system2(command = "aws",arg = c("s3","ls", paste0(studyDirectory, "logs/")), stdout = TRUE)
studyCreation <- all_files[grep("study_creation_", all_files)]
studyCreation <- strsplit(studyCreation, " ")[[1L]]
studyCreation <- studyCreation[length(studyCreation)]

studyCreation <- paste0(paste0(studyDirectory, "logs/"), studyCreation)

annots <- system2(command = "aws",arg = c("s3","ls", studyCreation), stdout = TRUE)
annots <- annots[grep(".csv", annots)]
annots <- annots[grep("Resources", annots, invert = TRUE)]
annots <- strsplit(annots, " ")[[1L]]
annots <- annots[length(annots)]
annots <- paste0(studyCreation, annots)

bucket <- strsplit(annots,"/")[[1L]][3L]
file_key <- strsplit(annots, paste0(bucket,"/"))[[1L]][2L]
system2(command = "aws", arg = c("s3api", "get-object","--bucket",bucket,"--key", file_key,paste0("/tmp/tmp.csv")), stdout = TRUE)
annots <- read.delim("/tmp/tmp.csv", header = TRUE, sep = ",")
file.remove("/tmp/tmp.csv")
rm(bucket)
rm(file_key)

# get slide names
slideName <- annots$Run_Tissue_name
slideName <- gsub("\\W", "", slideName)

studyName <- gsub("\\W", "", studyName)

# create local output folder
if(!dir.exists(studyName)){
  dir.create(studyName)
}

# get tiledb object name
tiledbName <- strsplit(study$uri, "/")[[1L]]
tiledbName <- tiledbName[length(tiledbName)]

if(!endsWith(outPath, "/")){
  outPath <- paste0(outPath, "/")
}

if(!is.null(FOVs) | !is.null(flowcell)){
  if(!is.null(FOVs)){
    FOVs <- eval(parse(text = FOVs))
    if(!is.numeric(FOVs)){
      stop("resulting FOVs must be numeric")
    }
  }
  
  if(!is.null(flowcell)){
    if(!is.numeric(flowcell) | length(flowcell) > 1){
      stop("flowcell must be 1 numeric value")
    }
    
    flowcellName <- annots$Run_Tissue_name[flowcell]
  }
  else{
    flowcellName <- NULL
  }
  
  
  config_subset <- list(cells = NULL, 
                        fovs = FOVs, 
                        targets = NULL, 
                        cellTypes = NULL,
                        cellTypeColumn = "RNA_nbclust_clusters",
                        slides = flowcellName, 
                        remove = FALSE,
                        and = TRUE,
                        removeFlagged = FALSE)
  
}else{
  config_subset <- NULL
}

if(tiledbArray == TRUE | (SeuratObject == TRUE & !is.null(config_subset))){
  message("Exporting TileDB")
  
  localTileDB <- paste0(studyName, "/", tiledbName, "_TileDB/")
  
  # Copy Tiledb Object
  system2(command = "aws", args = c("s3", "cp", study$uri, localTileDB, "--recursive"), stdout = FALSE)
  
  if(!is.null(config_subset)){
    study <- tiledbsc::SOMACollection$new(uri = localTileDB, verbose = FALSE)
    nanopipeline::runTileDBSubset(config = config_subset, scdataset = study)
  }
  
  if(SeuratObject == FALSE | is.null(config_subset)){
    # move folder to user S3 bucket
    system2(command = "aws", args = c("s3", "mv", studyName, paste0(outPath, studyName, "/"), "--recursive", "--profile", "export"), stdout = FALSE)
  }
}else{
  localTileDB <- NULL
}

# read in Seurat Object
if(SeuratObject == TRUE){
  message("Exporting Seurat Object")
  if(FullSeuratObject == TRUE){
    sem <- suppressWarnings(toFullSeuratObject(tiledbsc_dataset = study,transcripts = transcripts))
  }else{
    sem <- suppressWarnings(study$to_seurat(batch_mode = TRUE))
    Seurat::DefaultAssay(sem) <- "RNA"
    
    obs <- study$somas$RNA$obs$to_dataframe()
    sem@meta.data <- obs[match(colnames(sem), rownames(obs)),]
    
    for(i in names(sem@reductions)){
      sem@reductions[[i]]@assay.used <- "RNA"
    }
    
    if(transcripts == TRUE){
      sem@misc$transcriptCoords <- NULL
      for(i in names(study$somas)){
        if("transcriptCoords" %in% names(study$somas[[i]]$obsm$members)){
          sem@misc$transcriptCoords <- rbind(sem@misc$transcriptCoords,
                                             readInTranscriptCoords(study$somas[[i]]$obsm$members$transcriptCoords))
        }
      }
    }
  }
  
  if(tiledbArray == FALSE & !is.null(localTileDB)){
    unlink(localTileDB, recursive = TRUE)
  }
  
  # Save Seurat Object
  saveRDS(sem, paste0(studyName, "/", tiledbName, "_seuratObject.RDS"))
  
  rm(sem)
  gc()
  
  # move folder to user S3 bucket
  system2(command = "aws", args = c("s3", "mv", studyName, paste0(outPath, studyName, "/"), "--recursive", "--profile", "export"), stdout = FALSE)
}

if(rawFiles == TRUE){
  message("Exporting Raw Files")
  # Path to highest data folder
  paths <- strsplit(annots$slidefolders, "/")
  for(i in 1:nrow(annots)){
    paths[i] <- paste0(annots$folder_path[i], paste(paths[[i]][-length(paths[[i]])], collapse = "/"), "/")
  }
  
  paths <- unlist(paths)
  paths <- gsub("//", "/", paths)
  paths <- gsub("s3:/", "s3://", paths)
  paths <- unique(paths)
  
  if(!is.null(flowcell)){
    paths <- paths[flowcell]
    slideName <- slideName[flowcell]
    annots <- annots[flowcell,]
  }
  
  for(i in 1:length(paths)){
    print(paste(slideName[[i]], "-", i, "/", length(paths)))
    
    all_files <- system2(command = "aws", args = c("s3", "ls", paths[[i]], "--recursive"), stdout = TRUE)
    
    all_files <- unlist(lapply(all_files, function(x){
      temp_file_vec <- strsplit(x, " ")[[1L]]
      
      return(paste0(annots$folder_path[i], temp_file_vec[length(temp_file_vec)]))
    }))
    
    if(spotFiles == FALSE){
      all_files <- all_files[grep("SpotFiles", all_files, invert = TRUE)]
    }
    
    if(exportFOVImages == FALSE){
      all_files <- all_files[grep("png$|tif$|tiff$|jpg$", tolower(all_files), invert = TRUE)]
    }
    
    if(!is.null(FOVs)){
      FOVstrings <- stringr::str_pad(FOVs, width = 3, pad = 0)
      
      flowcellFiles <- grep(pattern = "FOV[0-9][0-9][0-9]|F[0-9][0-9][0-9]", all_files, invert = TRUE)
      
      FOVfiles <- NULL
      for(f in FOVstrings){
        FOVfiles <- c(FOVfiles, grep(pattern = paste0(paste0("FOV", f), "|", paste0("F", f)), all_files))
      }
      
      all_files <- all_files[c(flowcellFiles, FOVfiles)]
    }
    
    NULLList <- parallel::mclapply(all_files, function(currentFile){
      resultsFile <- gsub(paste0(annots$folder_path[i], dirname(annots$slidefolders[i]), "/"), "", currentFile)
      if(startsWith(prefix = basename(annots$slidefolders[i]), resultsFile) | startsWith(prefix = "Logs", resultsFile)){
        resultsFile <- paste0(studyName, "/", slideName[[i]], "/", resultsFile)
        
        # copy to working environment using AtoMx credentials
        system2(command = "aws", args = c("s3", "cp", currentFile, resultsFile), stdout = FALSE)
        # move folder to user S3 bucket
        system2(command = "aws", args = c("s3", "mv", resultsFile, paste0(outPath, resultsFile), "--profile", "export"), stdout = FALSE)
      }
    }, mc.cores = usableCores())
    
    rm(NULLList)
    gc()
  }
}

unlink(studyName, recursive = TRUE)
