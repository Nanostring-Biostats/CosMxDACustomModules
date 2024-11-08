# Name: Combine Sample Metadata Custom Module
message("Custom Script Version: 1.0.0")
# Supported AtoMx Versions: v1.3.2

# Copyright 2024 Bruker Spatial Biology, Inc.
# This software and any associated files are distributed pursuant to the 
# NanoString AtoMx Spatial Information Platform Software as a Service Agreement, 
# available on the NanoString Technologies, Inc website at www.nanostring.com, 
# as updated.  All rights reserved.  No permission is granted to modify, publish, 
# distribute, sublicense or sell copies of this software.

# Description: Combine logical queries on sample metadata to generate a binary 
#              column with results that can be used for filtering in later steps of analysis.

# columnName - Name of new metadata column
# column - Column(s) to query on
# logic - Logic for query (ex: '>', '<=', '=')
# value - Query comparison (ex: 5, 'B-cell')
# condition - If combining queries, AND or OR. If one query, leave blank
# fullCondition - Advanced feature: Full query condition. No error checking on this. 
#                 Overrules other variables if available.
# trueResult - Advanced feature: value for true query result, default 1
# falseResult - Advanced feature: value for false query result, default 0

# Load Packages
library(tiledb)

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
variableMsg <- variableTest(varName = "columnName", varType = "character", 
                            msg = variableMsg, required = TRUE)
variableMsg <- variableTest(varName = "column", varType = "character", 
                            msg = variableMsg, required = FALSE)
variableMsg <- variableTest(varName = "logic", varType = "character", 
                            msg = variableMsg, required = FALSE)
variableMsg <- variableTest(varName = "value", varType = "character", 
                            msg = variableMsg, required = FALSE)
variableMsg <- variableTest(varName = "condition", varType = "character", 
                            msg = variableMsg, required = FALSE)
variableMsg <- variableTest(varName = "fullCondition", varType = "character", 
                            msg = variableMsg, required = FALSE)
variableMsg <- variableTest(varName = "trueResult", varType = "character", 
                            msg = variableMsg, required = TRUE)
variableMsg <- variableTest(varName = "falseResult", varType = "character", 
                            msg = variableMsg, required = TRUE)

# stop on variable errors if applicable
if(!is.null(variableMsg)){
  stop(variableMsg)
}

updateCellMetadata <- function(uri, obs, indexCol = "cell_ID"){
  if(grepl("obs$", uri)){
    obs$obs_id <- obs[[indexCol]] # keep index column in data frame
    
    tiledb::tiledb_vfs_remove_dir(uri = uri)
    tiledb::fromDataFrame(obj = as.data.frame(obs), uri = uri, 
                          col_index = "obs_id")
  }else{
    stop("This function is only meant for cell metadata in the obs slot")
  }
}

# Module Code
obs <- study$somas$RNA$obs$to_dataframe()
cols <- study$somas$RNA$obs$attrnames()

# Ensure columnName is unique
if(columnName %in% cols){
  num <- 0
  columnNameRaw <- columnName
  while(columnName %in% cols){
    num <- num + 1
    columnName <- paste0(columnNameRaw, num)
  }
}

if(is.null(fullCondition) | fullCondition == ""){
  if(grepl(pattern = "^c\\(",x = column)){
    column <- eval(parse(text = column))
  }
  if(grepl(pattern = "^c\\(",x = logic)){
    logic <- eval(parse(text = logic))
  }
  if(grepl(pattern = "^c\\(",x = value)){
    value <- eval(parse(text = value))
  }
  
  if(grepl(pattern = "^c\\(",x = condition)){
    condition <- eval(parse(text = condition))
    if(length(condition) > 1){
      warning("Only first condition used")
      condition <- condition[1L]
    }
  }
  condition <- toupper(condition)
  
  # Check for valid condition
  if(condition == "AND"){
    condition <- "&" 
  }else if(condition == "OR"){
    condition <- "|" 
  }else if(condition == ""){
    condition <- "&" 
  }else if(is.null(condition)){
    condition <- "&" 
  }else{
    stop("Condition must be either AND or OR")
  }
  
  # if one logic is provided for multiple queries, replicate logic for all
  if(length(logic) == 1 & length(column) > 1){
    logic <- rep(logic, length(column))
    warning("Same logic is used for all queries")
  }
  
  # replicate conditions for all queries
  if(length(condition) == 1 & length(column) > 2){
    condition <- rep(condition, length(column)-1)
  }
  
  # ensure all equals logic matches R format
  if(any(logic == "=")){
    logic[which(logic == "=")] <- "==" 
  }
  
  # ensure length of column, value, and logic all match
  if(length(column) != length(value) | length(value) != length(logic)){
    stop("column, value, and logic lengths do not match, please review your inputs")
  }
  
  # ensure valid column names
  if(!any(column %in% cols)){
    stop("Invalid column name, please reference the sampleMetadataSummary.tsv file from GetSampleMetadata for column names")
  }
  
  # ensure valid column values
  err <- NULL
  for(i in seq_len(length(value))){
    if(is.na(suppressWarnings(as.numeric(value[i])))){
      if(!value[i] %in% unique(obs[[column[i]]])){
        err <- paste(err, paste0(value[i], " is not in column ", column[i]), sep = "\n")
      }
    }else{
      val <- as.numeric(value[i])
      if(val < min(obs[[column[i]]]) | val > max(obs[[column[i]]])){
        err <- paste(err, paste0(val, " is not in valid range for column ", column[i]), sep = "\n")
      }
    }
  }
  if(!is.null(err)){
    stop(paste0(err, "\nPlease refer to sampleMetadataSummary.tsv file from GetSampleMetadata for valid column values"))
  }
  
  # write full condition
  fullCondition <- NULL
  for(i in seq_len(length(column))){
    if(suppressWarnings(!is.na(as.numeric(value[i])))){
      val <- as.numeric(value[i])
    }else{
      val <- paste0("\"",value[i],"\"")
    }
    fullCondition <- paste0(fullCondition, "obs$", column[i], " ", logic[i], " ", val)
    if(i < length(column)){
      fullCondition <- paste0(fullCondition, " ", condition[i], " ")
    }
  }
}

# run query on condition
rows <- eval(parse(text = fullCondition))

if(sum(rows == TRUE) == 0){
  stop("Zero cells matched the query, please check spellings or update values")
}

if(!is.na(suppressWarnings(as.numeric(falseResult)))){
  falseResult <- as.numeric(falseResult)
}

if(!is.na(suppressWarnings(as.numeric(trueResult)))){
  trueResult <- as.numeric(trueResult)
}

if(class(trueResult) != class(falseResult)){
  falseResult <- as.character(falseResult)
  trueResult <- as.character(trueResult)
}

# create column from query results
obs[[columnName]] <- falseResult
obs[[columnName]][rows] <- trueResult

# update cell metadata
updateCellMetadata(uri = study$somas$RNA$obs$uri, obs = obs)

# print results to console
print(fullCondition)
print(columnName)
print(table(obs[[columnName]]))
