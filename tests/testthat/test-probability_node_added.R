test_that("probility is uniform when variables are equal", {
  d <- data.frame(
    aoa   = c(1, 1, 2),
    gval  = c(1, 1, 0),
    tstep = c(1, 1, 2),
    id    = c(1, 2, 3)
  )
  actual <- netgrowr:::probability_node_added(
    beta = c(1, 1),
    formula = aoa ~ gval,
    data = d,
    split_by = "tstep",
    label_with = "id"
  )
  expect_equal(actual, c("1" = 0.5, "2" = 0.5, "3" = 1.0))
})
test_that("unlearned words, but not previously learned words, contribute to computation", {
  d <- data.frame(
    aoa   = c(1, 1, 2, 2,  1,  1, 2, 2),
    tstep = c(1, 1, 1, 1,  2,  2, 2, 2),
    gval  = c(3, 2, 1, 1, NA, NA, 2, 1),
    id    = c(1, 2, 3, 4,  1,  2, 3, 4)
  )

  actual <- netgrowr:::probability_node_added(
    beta = c(1),
    formula = aoa ~ 1,
    d,
    split_by = "tstep",
    label_with = "id"
  )
  expect_equal(actual, c('1' = 0.25, '2' = 0.25, '3' = 0.5, '4' = 0.5))

  actual <- netgrowr:::probability_node_added(
    beta = c(1, 2),
    formula = aoa ~ gval,
    d,
    split_by = "tstep",
    label_with = "id"
  )
  expect_equal(as.vector(actual[1] > actual[2]), TRUE)
  expect_equal(as.vector(actual[3] + actual[4]), 1.0)
})
