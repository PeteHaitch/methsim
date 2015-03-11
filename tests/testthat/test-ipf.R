### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ipf
###
context("ipf function")

test_that("ipf gives identical output to mipfp::ipf", {
  or <- 1
  seed <- matrix(c(or, 1, 1, 1), ncol = 2)
  row_margins <- c(0.2, 0.8)
  col_margins <- c(0.73, 0.27)
  expect_identical(
    mipfp::Ipfp(seed = seed,
                target.list = list(1, 2),
                target.data = list(row_margins, col_margins),
                print = FALSE,
                iter = 1000,
                tol = 1e-10)$xi.hat,
    ipf(seed = seed,
        row_margins = row_margins,
        col_margins = col_margins,
        iter = 1000,
        tol = 1e-10)
  )
  or <- 3.1236
  seed <- matrix(c(or, 1, 1, 1), ncol = 2)
  row_margins <- c(0.44, 0.56)
  col_margins <- c(0.91, 0.09)
  expect_identical(
    mipfp::Ipfp(seed = seed,
                target.list = list(1, 2),
                target.data = list(row_margins, col_margins),
                print = FALSE,
                iter = 1000,
                tol = 1e-10)$xi.hat,
    ipf(seed =  seed,
          row_margins = row_margins,
          col_margins = col_margins,
          iter = 1000,
          tol = 1e-10)
  )
})
