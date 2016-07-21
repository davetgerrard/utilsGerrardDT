
# returns the current script if called from command line. Empty character string if not.
getCurrentScript <- function() {
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  return(script.name)
}