library(PTXQC)
QCReport <- function(mq.path){
  file.path <- paste(getwd(), mq.path, sep = '/')
  createReport(file.path)
}