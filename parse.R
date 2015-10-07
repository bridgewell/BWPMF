Sys.setenv("PKG_CXXFLAGS" = "-std=c++11", "PKG_LIBS" = "-lboost_serialization -lboost_iostreams")
Rcpp::sourceCpp("parse.cpp", verbose = TRUE)
src <- serialize_cookie()
deserialize_cookie(src)

function() {
chunk_size <- 1000
con <- file("2015-10-01.txt", "r")
counter.max <- counter <- 100
gc_times <- 0
pb <- txtProgressBar(max = 98272059, style = 3)
while(length(src <- readLines(con, n = chunk_size)) > 0) {
  encode(src)
  if (counter == 0) {
    counter <- 100
    gc()
    gc_times <- gc_times + 1
    setTxtProgressBar(pb, gc_times * chunk_size * counter.max)
  } else counter <- counter - 1
}
close(pb)
close(con)
}