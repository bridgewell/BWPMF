#'@export
deserialize_cookie <- function(src) {
  switch(class(src),
         "raw" = deserialize_cookie_raw(src),
         "character" = deserialize_cookie_path(src)
  )
}

#'@export
deserialize_hostname <- function(src) {
  switch(class(src),
         "raw" = deserialize_hostname_raw(src),
         "character" = deserialize_hostname_path(src)
  )
}

#'@export
deserialize_history <- function(src) {
  switch(class(src),
         "raw" = deserialize_history_raw(src),
         "character" = deserialize_history_path(src))
}