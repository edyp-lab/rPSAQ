## code to prepare `10ng_r4_sample` dataset goes here

sample_10ng_r4 = openxlsx::read.xlsx("inst/extdata/7843_1ng_4.xlsx", sheet = 1, colNames = TRUE)

usethis::use_data(sample_10ng_r4, overwrite = TRUE)
