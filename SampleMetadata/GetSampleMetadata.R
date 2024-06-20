# Name: Get Sample Metadata Custom Module
message("Custom Script Version: 1.0")
# Supported AtoMx Versions: v1.3.2

# Copyright 2024 Bruker Spatial Biology, Inc.
# This software and any associated files are distributed pursuant to the NanoString AtoMx Spatial 
# Information Platform Software as a Service Agreement, available on the NanoString Technologies, Inc 
# website at www.nanostring.com, as updated.  All rights reserved.  No permission is granted to modify, 
# publish, distribute, sublicense or sell copies of this software.

# Description: Download sample metadata file

# User Defined Variables

# fullDataset - return the full metadata array
# qcFlagsOnly - only return the qcFlags and index columns 
# indexOnly - only return the index columns
# cellLevel - return cell level metadata 
# fovLevel - return FOV level metadata
# flowcellLevel - return flowcell level metadata
# metaSummary - return summary of metadata columns

# Load Packages
library(nanopipeline)

# Test User Variables
variableTest <- function(varName, varType, msg, required = TRUE){
  if(!varType %in% c("character", "logical", "numeric", "file")){
    stop("varType must be \"character\", \"logical\", \"numeric\", or \"file\"")
  }
  
  varMsg <- NULL
  
  typeClass <- switch(varType, 
                      character={"STRING"},
                      logical={"BOOL"},
                      numeric={"NUM"},
                      file={"FILE"})
  
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
    
    if(varType == "file"){
      if(!is(get(varName), "character")){
        if(!is.null(get(varName)) | required == TRUE){
          varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" varType was not set as ", typeClass, 
                                          " in custom module creation"))
        }
      }
      if(!file.exists(get(varName))){
        varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" was not uploaded properly, please reupload"))
      }
    }else{
      if(!is(get(varName), varType)){
        if(!is.null(get(varName)) | required == TRUE){
          varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" varType was not set as ", typeClass, 
                                          " in custom module creation"))
        }
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
variableMsg <- variableTest(varName = "fullDataset", varType = "logical", required = TRUE, msg = variableMsg)
variableMsg <- variableTest(varName = "qcFlagsOnly", varType = "logical", required = TRUE, msg = variableMsg)
variableMsg <- variableTest(varName = "indexOnly", varType = "logical", required = TRUE, msg = variableMsg)
variableMsg <- variableTest(varName = "cellLevel", varType = "logical", required = TRUE, msg = variableMsg)
variableMsg <- variableTest(varName = "fovLevel", varType = "logical", required = TRUE, msg = variableMsg)
variableMsg <- variableTest(varName = "flowcellLevel", varType = "logical", required = TRUE, msg = variableMsg)
variableMsg <- variableTest(varName = "metaSummary", varType = "logical", required = TRUE, msg = variableMsg)

# stop on variable errors if applicable
if(!is.null(variableMsg)){
  stop(variableMsg)
}

if(fullDataset == FALSE & qcFlagsOnly == FALSE){
  indexOnly <- TRUE
}

# Module Code:
outDir <- "/output"
dir.create(outDir)
setwd(outDir)
zipFile <- paste0("sampleMetadata.zip")

#' @title Get column names for unique columns
#' 
#' @description From cell metadata, get columns that are unique for a subset of columns.
#' This allows subsetting of cell metadata to FOV or flowcell level data. 
#'
#' @param columns Column to subset to for comparison: FOV and/or slide
#' @param obs Cell metadata
#' 
levelData <- function(columns, obs){
  comp <- unique(obs[,columns])
  if(is(comp,"numeric")){
    comp <- length(comp)
  }else{
    comp <- nrow(comp)
  }
  cols <- unlist(parallel::mclapply(colnames(obs), function(i){
    temp <- unique(obs[,c(columns, i)])
    if(nrow(temp) == comp){
      return(i)
    }
  }, mc.cores = nanopipeline::usableCores()))

  cols <- unique(c(columns,cols))

  return(unique(obs[,cols]))
}

obs <- study$somas$RNA$obs$to_dataframe()

outFiles <- NULL

# subset metadata down to FOV and/or flowcell level data
if(cellLevel == TRUE){
  index <- c("cell_ID", "slide_ID_numeric", "fov")
  
  if(fullDataset){
    outFile <- "cellMetadata.csv"
    outFiles <- c(outFiles, outFile)
    
    data.table::fwrite(file = outFile, x = obs, quote = FALSE, 
                       row.names = FALSE, sep = ",")
  }
  
  if(qcFlagsOnly){
    outFile <- "cellMetadata_qc.csv"
    outFiles <- c(outFiles, outFile)
    
    cols <- colnames(obs)[grep("qcFlags|\\.qc", colnames(obs))]
    obsCells <- obs[,c(index, cols), drop = FALSE]
    data.table::fwrite(file = outFile, x = obsCells, quote = FALSE, 
                       row.names = FALSE, sep = ",")
  }
  
  if(indexOnly){
    outFile <- "cellMetadata_index.csv"
    outFiles <- c(outFiles, outFile)
    
    obsCells <- obs[,index, drop = FALSE]
    data.table::fwrite(file = outFile, x = obsCells, quote = FALSE, 
                       row.names = FALSE, sep = ",")
  }
}

if(fovLevel == TRUE){
  index <- c("slide_ID_numeric", "fov")
  obsFOV <- levelData(columns = index, obs)
  
  if(fullDataset){
    outFile <- "fovMetadata.csv"
    outFiles <- c(outFiles, outFile)
    
    data.table::fwrite(file = outFile, x = obsFOV, quote = FALSE, 
                       row.names = FALSE, sep = ",")
  }
  
  if(qcFlagsOnly){
    outFile <- "fovMetadata_qc.csv"
    outFiles <- c(outFiles, outFile)
    
    cols <- colnames(obsFOV)[grep("qcFlags|\\.qc", colnames(obsFOV))]
    obsFOV <- obsFOV[,c(index, cols), drop = FALSE]
    data.table::fwrite(file = outFile, x = obsFOV, quote = FALSE, 
                       row.names = FALSE, sep = ",")
  }
  
  if(indexOnly){
    outFile <- "fovMetadata_index.csv"
    outFiles <- c(outFiles, outFile)
    
    obsFOV <- obsFOV[,index, drop = FALSE]
    data.table::fwrite(file = outFile, x = obsFOV, quote = FALSE, 
                       row.names = FALSE, sep = ",")
  }
}

if(flowcellLevel == TRUE){
  index <- "slide_ID_numeric"
  obsFC <- levelData(index, obs)
  obsFC <- obsFC[,which(!colnames(obsFC) %in% c("qcFlagsFOV"))]
  
  if(fullDataset){
    outFile <- "flowcellMetadata.csv"
    outFiles <- c(outFiles, outFile)
    
    data.table::fwrite(file = outFile, x = obsFC, quote = FALSE, 
                       row.names = FALSE, sep = ",")
  }
  
  if(qcFlagsOnly){
    outFile <- "flowcellMetadata_qc.csv"
    outFiles <- c(outFiles, outFile)
    
    cols <- colnames(obsFC)[grep("qcFlags|\\.qc", colnames(obsFC))]
    obsFC <- obsFC[,c(index, cols), drop = FALSE]
    data.table::fwrite(file = outFile, x = obsFC, quote = FALSE, 
                       row.names = FALSE, sep = ",")
  }
  
  if(indexOnly){
    outFile <- "flowcellMetadata_index.csv"
    outFiles <- c(outFiles, outFile)
    
    obsFC <- obsFC[,index, drop = FALSE]
    data.table::fwrite(file = outFile, x = obsFC, quote = FALSE, 
                       row.names = FALSE, sep = ",")
  }
}

if(metaSummary == TRUE){
  vals <- NULL
  
  for(i in colnames(obs)){
    colClass <- class(obs[[i]])
    
    if(inherits(obs[[i]], "numeric") |
       inherits(obs[[i]], "integer")){
      n <- ""
      values <- round(summary(obs[[i]]), 4)
    }else{
      n <- length(unique(obs[[i]]))
      values <- table(obs[[i]])[1:(min(40,n))]
    }
    
    names(values) <- gsub("\\W", "", names(values))
    
    if(n == nrow(obs)){
      values <- paste0(names(values)[1:20], collapse = ",")
    }else{
      values <- paste0(paste0(names(values), ":", values), collapse = ",")
    }
    
    
    vals <- rbind.data.frame(vals,cbind(columnName=i, class=colClass, 
                                        nUnique=n, 
                                        "values:counts"=values))
  }
  
  outFile <- paste0("sampleMetadataSummary.tsv")
  outFiles <- c(outFiles, outFile)
  data.table::fwrite(vals, outFile, quote = FALSE, sep = "\t", row.names = FALSE)
}

zip(zipfile = zipFile, files = outFiles)

fileSize <- as.numeric(strsplit(system2("ls", args = c("-al", zipFile), 
                                    stdout = TRUE), split = " ")[[1L]][5])

if(fileSize > 6e+8){
  if("cellMetadata.csv" %in% outFiles){
    unlink("cellMetadata.csv")
  }
  warning("Zip folder is too large for local download, please download cell metadata using the Flat File Export custom module or button")
}

unlink(zipFile)
