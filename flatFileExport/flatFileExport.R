# Name: Flat File Export Custom Module
message("Custom Script Version: 1.2.0")

# Copyright 2023-2025 Bruker Spatial Biology, Inc.
# This software and any associated files are distributed pursuant to the Bruker Spatial Biology AtoMx Spatial
# Informatics Platform Software as a Service Agreement, available on the Bruker Spatial Biology, Inc
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
outDir <- "output/"
dir.create(outDir)

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
             paste("aws_secret_access_key =", secret_key), session_token,
             paste("ignore_configured_endpoint_urls =", "true")), fileConn)
close(fileConn)

############################### Helper functions ###############################
getAnnots <- function(){
  # read in annotation file
  # get input folder location and slide names
  all_files <- s3_ls(paste0(studyDirectory, "logs/"), arry = study, fileType = "folder")
  studyCreation <- all_files[grep("study_creation_", all_files)]
  
  annots <- s3_ls(studyCreation, arry = study, fileType = "folder")
  annots <- annots[grep(".csv", annots)]
  annots <- annots[grep("Resources", annots, invert = TRUE)]
  
  return(annots)
}

annots <- getAnnots()


# DO NOT CHANGE THESE VALUES ERRORS OR INCORRECT DATA COULD OCCUR 
config_loading <- list(annotfile = annots,
                       folderpathColumn = "folder_path", 
                       slidefoldersColumn = "slidefolders", 
                       slidenameColumn = "Run_Tissue_name",
                       storage_type = "s3")

config_align <- list(locfile = "../Logs/SpatialProfiling_*.fovs.csv", 
                     millimeter_conversion_toggle = TRUE,
                     instrument_version = NULL, 
                     pixel_conversion_value = NULL,
                     zstep_mm = NULL)

s3_get(S3Path = annots, arry = study, localPath = "annots.csv", fileType = "file")
annots <- read.delim("annots.csv", header = TRUE, sep = ",")
file.remove("annots.csv")

flatFileExport(outPath = outDir, study = study, annots = annots,
               countMatrix = countMatrix, cellMetadata = cellMetadata, 
               transcriptsFile = transcripts, polygons = polygons, 
               FOVpositions = fovPositions, config_align = config_align, 
               config_loading = config_loading)

################################## Move to S3 ##################################
system2(command = "aws", args = c("s3", "mv", outDir, paste0(outPath, studyName, "/"),
                                  "--recursive", "--profile", "export"), stdout = FALSE)

unlink(outDir, recursive = TRUE)
