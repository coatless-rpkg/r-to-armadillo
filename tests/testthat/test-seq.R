context("seq Set of Functions in Armadillo")

test_that("Test seq_default: Sequence Integer", {

  expect_equal(seq_int(-1, 2), as.matrix(seq(from = -1, to = 2, by = 1)))

  expect_equal(seq_int(2, 1), as.matrix(seq(from = 2, to = 1, by = -1)))
})


test_that("Test seq_default: Sequence Default", {

  expect_equal(seq_default(-1, 2, 10), as.matrix(seq(from = -1, to = 2, length.out = 10)))

  expect_equal(seq_default(2, 1, 10), as.matrix(seq(from = 2, to = 1, length.out = 10)))
})

test_that("Test seq_default_a: Sequence A Default", {

  expect_equal(seq_default_a(-1, 12), as.matrix(seq(from = 1, to = -1, length.out = 12)))

  expect_equal(seq_default_a(2, 4), as.matrix(seq(from = -2, to = 2, length.out = 4)))
})

test_that("Test seq_along: Sequence Along", {

  expect_equal(seq_along_cpp(1:10), as.matrix(seq_along(1:10) - 1))

  expect_equal(seq_along_cpp(5:10), as.matrix(seq_along(5:10) - 1))
})


test_that("Test seq_Length: Sequence Length", {

  expect_equal(seq_len_cpp(2), as.matrix(seq_len(2) - 1))

  expect_equal(seq_len_cpp(5), as.matrix(seq_len(5) - 1))

})
