# Copyright (C) 2015 Peter Hickey
#
# This file is part of methsim.
#
# methsim is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# methsim is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with methsim  If not, see <http://www.gnu.org/licenses/>.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### findMostOverlapping
###

context("findMostOverlapping")

test_that("findMostOverlapping works for simple case", {
  query <- GRanges(paste0("chr", c(1, 2)),
                   IRanges(start = c(1, 1),
                           end = c(20, 20)))
  subject <- GRanges(paste0("chr", c(1, 1, 2, 2)),
                     IRanges(start = c(1, 10, 1, 5),
                             end = c(10, 20, 10, 15)))
  expect_identical(findMostOverlapping(query, subject), c(2L, 4L))
})

test_that("findMostOverlapping breaks ties at random", {
  query <- GRanges(paste0("chr", c(1, 2)),
                   IRanges(start = c(1, 1),
                           end = c(19, 19)))
  subject <- GRanges(paste0("chr", c(1, 1, 2, 2)),
                     IRanges(start = c(1, 10, 1, 5),
                             end = c(10, 20, 10, 15)))
  set.seed(1)
  expect_identical(findMostOverlapping(query, subject), c(1L, 4L))
  set.seed(2)
  expect_identical(findMostOverlapping(query, subject), c(1L, 4L))
})

test_that("findMostOverlapping returns NA if there is no hit", {
  seqinfo <- Seqinfo(paste0("chr", 1:3))
  query <- GRanges(paste0("chr", c(1, 2)),
                   IRanges(start = c(1, 1),
                           end = c(20, 20)),
                   seqinfo = seqinfo)
  subject <- GRanges(paste0("chr", c(1, 1, 3, 3)),
                     IRanges(start = c(1, 10, 1, 5),
                             end = c(10, 20, 10, 15)),
                     seqinfo = seqinfo)
  expect_identical(findMostOverlapping(query, subject), c(2L, NA_integer_))
  expect_identical(findMostOverlapping(query[2], subject), NA_integer_)
})

test_that("findMostOverlapping ignore.strand works", {
  query <- GRanges(paste0("chr1"),
                   IRanges(start = 1,
                           end = 20),
                   strand = "+")
  subject <-  GRanges(paste0("chr", c(1, 1)),
                      IRanges(start = c(1, 1),
                              end = c(5, 20)),
                      strand = c("+", "-"))
  expect_identical(findMostOverlapping(query, subject), 1L)
  expect_identical(findMostOverlapping(query, subject, ignore.strand = TRUE),
                   2L)
})

test_that("findMostOverlapping returns error for GTuples/GTuplesList", {
  expect_error(findMostOverlapping(GTuples(), GRanges()),
               paste0("'findMostOverlapping' not yet implemented for ",
                      "'GTuples'- and GTuplesList'-based objects."))
  expect_error(findMostOverlapping(GRanges(), GTuples()),
               paste0("'findMostOverlapping' not yet implemented for ",
                      "'GTuples'- and GTuplesList'-based objects."))
  expect_error(findMostOverlapping(GTuplesList(), GRanges()),
               paste0("'findMostOverlapping' not yet implemented for ",
                      "'GTuples'- and GTuplesList'-based objects."))
  expect_error(findMostOverlapping(GRanges(), GTuplesList()),
               paste0("'findMostOverlapping' not yet implemented for ",
                      "'GTuples'- and GTuplesList'-based objects."))

})
