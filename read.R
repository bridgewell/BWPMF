suppressPackageStartupMessages({
  library(RcppHDFS)
})
# library(magrittr)
# library(data.table)
# library(dplyr)
path.fmt <- "/wush/userchannelcount/2015-10-01.txt/part-%05d"
# # for(i in 0:2343) {
out <- file("2015-10-01.txt", "wb")
for(i in 0:312) {
  writeBin(tmp <- hdfs.read(sprintf(path.fmt, i)), out)
}