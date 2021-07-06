test_that("probility is uniform when variables are equal", {
  d <- data.frame(
    aoa = c(1, 1),
    a = c(1, 1),
    b = c(1, 1),
    learned = c(TRUE, TRUE),
    unknown = c(TRUE, TRUE)
  )
  expect_equal(netgrowr:::ratio_of_strengths(d, aoa ~ 0 + a + b, beta = c(1,1)), matrix(c(0.5, 0.5), nrow = 2, ncol = 1, dimnames = list(1:2, NULL)))
  expect_equal(netgrowr:::ratio_of_strengths(d, aoa ~ 0 + a + b, beta = c(1,2)), matrix(c(0.5, 0.5), nrow = 2, ncol = 1, dimnames = list(1:2, NULL)))
})
test_that("unlearned words contribute to computation but are not reported in output", {
  d <- data.frame(
    aoa = c(1, 1),
    a = c(1, 1),
    b = c(1, 1),
    learned = c(FALSE, TRUE),
    unknown = c(TRUE, TRUE)
  )
  expect_equal(netgrowr:::ratio_of_strengths(d, aoa ~ 0 + a + b, beta = c(1,1)), matrix(0.5, dimnames = list("2", NULL)))
  expect_equal(netgrowr:::ratio_of_strengths(d, aoa ~ 0 + a + b, beta = c(1,2)), matrix(0.5, dimnames = list("2", NULL)))
})

test_that("weights (and ratio of strengths) applied correctly", {
  d <- data.frame(
    aoa = c(1, 1),
    a = c(1, 1),
    b = c(1, 0),
    learned = c(TRUE, TRUE),
    unknown = c(TRUE, TRUE)
  )
  actual <- netgrowr:::ratio_of_strengths(d, aoa ~ 0 + a + b, beta = c(1,0))
  expect_equal(actual, matrix(c(0.5, 0.5), nrow = 2, dimnames = list(1:2, NULL)), ignore_attr = TRUE)

  actual <- netgrowr:::ratio_of_strengths(d, learned ~ 0 + a + b, beta = c(0,1))
  x <- exp(c(1, 0))
  expect_equal(actual, matrix(x / sum(x), nrow = 2, dimnames = list(1:2, NULL)))
})
