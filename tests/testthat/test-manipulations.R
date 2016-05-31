context("manipulations")

test_that("Subset Non-connected Regions: get_elements", {

  m = matrix(1:12, nrow = 4)

  row_index = c(1, 2, 1)
  col_index = c(2, 2, 3)

  expect_equal(get_elements(m, row_index - 1, col_index - 1),
               as.matrix(m[cbind(row_index, col_index)]),
               check.attributes = F)

})

