argv <- commandArgs(TRUE)
library(BWPMF)

encode(argv[1])

saveRDS(list(cookie = serialize_cookie(), hostname = serialize_hostname()), 
        sprintf("%s.encode", argv))
