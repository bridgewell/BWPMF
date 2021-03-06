library(BWPMF)
src.path <- system.file("2015-10-01-100.txt", package = "BWPMF")
encode(src.path)
history <- encode_data(src.path)
count_cookie()
count_hostname()
history_size <- check_history(history)
stopifnot(history_size == 3970)
history_non_zero_size <- count_non_zero_of_history(history)
stopifnot(history_non_zero_size == 885)
testing_id <- c(154, 397, 513, 818, 273, 3, 862, 635)
testing_history <- extract_history(training_history <- history, testing_id)

stopifnot(check_history(training_history) + check_history(testing_history) == history_size)
stopifnot(count_non_zero_of_history(training_history) + count_non_zero_of_history(testing_history) == history_non_zero_size)


m1 <- init_model(.1, .1, .1, .1, .1, .1, 10, training_history)
phi1 <- init_phi(m1, training_history)
m2 <- new(BWPMF::Model, m1)
phi2 <- init_phi(m1, training_history, .tmp_path <- tempfile(), 10)

train_once(m1, training_history, phi1, function(msg) {})
stopifnot(!isTRUE(all.equal(m1$export_user(), m2$export_user())))

train_once(m2, training_history, phi2, function(msg) {})
stopifnot(max(abs(m1$export_user() - m2$export_user())) < 1e-5)
stopifnot(isTRUE(all.equal(dim(dphi1 <- dump_phi(phi1)), dim(dphi2 <- dump_phi(phi2, training_history, 10)))))
stopifnot(max(abs(dphi1 - dphi2)) < 1e-5)
stopifnot(max(abs(m1$export_item() - m2$export_item())) < 1e-5)
