v1 = (c(1,1,1,0))
v2 = c(0,0,0,1)
v3 = c(1,1,1,0)
v4 = c(0,0,0,1)

w = matrix(data = c(v1,v2,v3,v4), ncol = 4, nrow = 4)
w = Matrix::Matrix(w, sparse = TRUE)
rownames = c("A", "B", "C", "D")
seed_proteins = c("A", "C")
seed = 1
p <- sparseRWR(seed, w, norm = TRUE)

test_that("sparseRWR returns a list of 3 ouputs of the correct types and lengths", {
  expect_true(is.double(p[[1]]))
  expect_equal(length(p[[1]]), 4)
  expect_true(is.double(p[[2]]))
  expect_equal(length(p[[2]]), 1)
  expect_equal(length(p[[3]]), 1)
  expect_true(is.integer(p[[3]]))
})

test_that("sparseRWR completes in < 1000 iterations", {
  expect_lt(p[[3]], 1000)
})

test_that("sparseRWR output equals expected", {
  expect_equal(sparseRWR(1,w, norm = TRUE)[[1]], c(0.709, 0, 0.109, 0), tolerance = 0.001)
  expect_equal(sparseRWR(2,w, norm = TRUE)[[1]], c(0.1090909, 0.6, 0.1090909, 0), tolerance = 0.001)
  expect_equal(sparseRWR(3,w, norm = TRUE)[[1]], c(0.1090909, 0, 0.7090909, 0), tolerance = 0.001)
  expect_equal(sparseRWR(4,w, norm = TRUE)[[1]], c(0.07272, 0.4, 0.07272, 1), tolerance = 0.001)
  expect_equal(sparseRWR(c(1,2), w, norm = TRUE)[[1]],  c(0.40909, 0.3, 0.10909, 0), tolerance = 0.001)
  expect_lte(sum(sparseRWR(c(1,3), w, norm = TRUE)[[1]]), 1)
})

test_that("sparseRWR doesn't break if you give seeds that don't exist in the data", {
  expect_equal(sparseRWR(c(1,6),w, norm = TRUE)[[1]], c(0.709, 0, 0.109, 0), tolerance = 0.001)
})

g <- igraph::sample_gnp(n = 1000, p = 0/1000)
w <- igraph::as_adjacency_matrix(g)

test_that("sparseRWR breaks for really sparse graphs", {
  expect_error(sparseRWR(c(1,2,3,4), w = w))
})


