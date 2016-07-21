addCaption <- function(filename, caption, suffix=".txt", verbose=T) {
  # writes a text file 
  if(is.null(filename) ) {
    warning("No filename given to addCaption!")
  } else {
  if(is.null(caption) ) {
    warning("No caption given to addCaption!")
  } else {
    captionFile <- paste0(filename, suffix)
   write(caption, file=captionFile)
   if(verbose)  print(paste("Wrote caption to ", captionFile)) 
  }
  
}
}