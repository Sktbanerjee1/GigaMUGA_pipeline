source("src/func_room/dependency_manager.R")

deps <- c(
  "optparse",
  "readr",
  "dplyr",
  "ggplot2",
  "logger",
  "plyr"
)

# install dependencies
dependency_manager(deps)

# load dependencies
sapply(deps, function(package) {
  library(package, character.only = TRUE)
})