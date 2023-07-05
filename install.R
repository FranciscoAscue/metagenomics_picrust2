# Package R dependencies
if( !is.element("devtools",rownames(installed.packages() ) ) ){
  install.packages("devtools")
  install.packages("BiocManager")
  
}

library(devtools)

################################################################################

# Install missing packages

missingPackages <- function(pkg){
  if( !is.element(pkg,rownames(installed.packages() ) ) ){
    message(pkg, "-----> Package is not installed ")
    BiocManager::install(pkg)
  }
}

dependencies <- c("dada2","phyloseq","Biostrings",
                  "ape","microbiome","DT","vegan","phyloseq",
                  "tidyverse","magrittr","readxl","multcomp","mvtnorm",
                  "survival","TH.data","MASS")

for(i in dependencies){
  missingPackages(i)
  library(i, character.only = TRUE)
}

if( !is.element("microbiomeutilities",rownames(installed.packages() ) ) ){
  devtools::install_github("microsud/microbiomeutilities")
}
library(microbiomeutilities)
################################################################################
