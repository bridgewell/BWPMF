argv <- commandArgs(TRUE)
library(BWPMF)

local({
  encode <- readRDS(sprintf("%s.encode", argv))
  deserialize_cookie(encode$cookie)
  deserialize_hostname(encode$hostname)
})
gc()
user_data <- encode_data(argv[1])
check_history(user_data)
saveRDS(serialize_history(user_data), sprintf("%s.history", argv))
