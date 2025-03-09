library(shiny)
library(shinydashboard)
library(markdown)
library(shinyjs)
library(tools)
library(Seurat)
library(dplyr)
library(DT)
library(shinydashboardPlus)
library(ggplot2)
library(shinybusy)
library(glue)
library(ggthemes)
library("readxl")
library(tidyr)
library(hash) 

# Read in file and perform validation.
load_seurat_obj <- function(path){
  errors <- c()
  # check file extension
  if (!tolower(tools::file_ext(path)) == "xlsx") { # ignores case
    errors <- c(errors, "Invalid xlsx file.")
    return(errors)
  }
  
  # try to reead in fil
  tryCatch(
    {
      obj <-read_excel(path)
    },
    error = function(e) {
      errors <- c(errors, "Invalid xlsx file.")
      return(errors)
    }
  )
}




