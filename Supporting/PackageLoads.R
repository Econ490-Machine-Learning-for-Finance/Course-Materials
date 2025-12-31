# Purpose: Load packages
# Notes: 1) checks if installed. 2) if not, then installs, then 3) calls library
# Usage: 
## Add to the pkgs list as needed
## source("PackageLoads.R") within other files as needed. 
# Authors
## Mike Aguilar (mike.aguilar@duke.edu) & Ziming Huang

# Define the packages you need
pkgs <- c("tidyverse", "data.table", "lubridate", 
          "dplyr", "tidyr", "ggplot2",
          "zoo", "PerformanceAnalytics", "gridExtra",
          "reshape2", "corrplot", "scales", 
          "moments", "lubridate", "quantmod",
          "tidyquant")

# Install any missing packages
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}
