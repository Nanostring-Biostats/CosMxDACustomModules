# Name: Update Sample Metadata Custom Module
message("Custom Script Version: 1.0")
# Supported AtoMx Versions: v1.3.2

# Copyright 2024 Bruker Spatial Biology, Inc.
# This software and any associated files are distributed pursuant to the NanoString AtoMx Spatial 
# Information Platform Software as a Service Agreement, available on the NanoString Technologies, Inc 
# website at www.nanostring.com, as updated.  All rights reserved.  No permission is granted to modify, 
# publish, distribute, sublicense or sell copies of this software.

# Description: Update sample metadata from uploaded file

# User Defined Variables
# metaFile - File containing new metadata

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
variableMsg <- variableTest(varName = "metaFile", varType = "file", msg = variableMsg)

# stop on variable errors if applicable
if(!is.null(variableMsg)){
  stop(variableMsg)
}

# Module Code

updateCellMetadata <- function(uri, obs, indexCol = "cell_ID"){
  if(grepl("obs$", uri)){
    obs$obs_id <- obs[[indexCol]] # keep index column in data frame
    
    tiledb::tiledb_vfs_remove_dir(uri = uri)
    tiledb::fromDataFrame(obj = as.data.frame(obs), uri = uri, col_index = "obs_id")
  }else{
    stop("This function is only meant for cell metadata in the obs slot")
  }
}

# get metadata for comparison
obs <- study$somas$RNA$obs$to_dataframe()

# unzip if necessary
if(endsWith(metaFile, ".zip")){
  metaFile <- unzip(metaFile)
}

# read in new metadata file
newObs <- data.table::fread(metaFile, data.table = FALSE)

newCols <- colnames(newObs)[which(!colnames(newObs) %in% colnames(obs))]

if(length(newCols) == 0){
  warning("No new columns added to given file, skipping")
}

cellCols <- c("cell_ID")
fovCols <- c("slide_ID_numeric","fov")
fcCols <-c("slide_ID_numeric")

# add FOV or flowcell metadata if given
if(nrow(newObs) == length(unique(newObs[,fcCols]))){
  rowMatch <- match(paste(obs[,fcCols]),paste(newObs[,fcCols]))
  
  for(i in newCols){
    obs[[i]] <- newObs[rowMatch, i]
  }
}else if(nrow(newObs) == nrow(unique(newObs[,fovCols]))){
  rowMatch <- match(paste(obs[,fovCols[1]], obs[,fovCols[2]]),
                    paste(newObs[,fovCols[1]], newObs[,fovCols[2]]))
  
  for(i in newCols){
    obs[[i]] <- newObs[rowMatch, i]
  }
}else if(nrow(obs) != nrow(newObs)){
  stop("All rows must be present")
}else{
  rowMatch <- match(obs[,cellCols], newObs[,cellCols])
  for(i in newCols){
    obs[[i]] <- newObs[rowMatch, i]
  }
}

updateCellMetadata(uri = study$somas$RNA$obs$uri, obs = obs)


