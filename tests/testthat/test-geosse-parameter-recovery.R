test_that("identical epochs reproduce the ordinary ClaSSE likelihood", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  states <- c(a = 1, b = 2, c = 3)

  geosse.parameters <- c(
    sA = 0.20,
    sB = 0.20,
    sAB = 0.20,
    xA = 0.05,
    xB = 0.05,
    dA = 0.10,
    dB = 0.10
  )

  classe.parameters <- pars.geosse2classe(geosse.parameters)
  parameter.names <- names(classe.parameters)

  ordinary.likelihood <- diversitree::make.classe(
    tree,
    states,
    k = 3
  )

  episodic.likelihood <- make.classe.td(
    tree,
    states,
    k = 3,
    n.epoch = 2
  )

  episodic.parameters <- c(
    t.1 = 0.5,
    setNames(classe.parameters, paste0(parameter.names, ".1")),
    setNames(classe.parameters, paste0(parameter.names, ".2"))
  )

  ordinary.value <- ordinary.likelihood(
    classe.parameters,
    condition.surv = FALSE,
    root = diversitree::ROOT.GIVEN,
    root.p = c(1, 0, 0)
  )

  episodic.value <- episodic.likelihood(
    episodic.parameters,
    condition.surv = FALSE,
    root = diversitree::ROOT.GIVEN,
    root.p = c(1, 0, 0)
  )

  expect_equal(
    episodic.value,
    ordinary.value,
    tolerance = 1e-7
  )
})

test_that("GeoSSE conversion places rates in the expected ClaSSE entries", {
  geosse.parameters <- c(
    sA = 0.11,
    sB = 0.22,
    sAB = 0.33,
    xA = 0.04,
    xB = 0.05,
    dA = 0.16,
    dB = 0.17
  )

  classe.parameters <- pars.geosse2classe(geosse.parameters)

  expect_equal(unname(classe.parameters["lambda111"]), 0.11)
  expect_equal(unname(classe.parameters["lambda313"]), 0.11)
  expect_equal(unname(classe.parameters["lambda222"]), 0.22)
  expect_equal(unname(classe.parameters["lambda323"]), 0.22)
  expect_equal(unname(classe.parameters["lambda312"]), 0.33)
  expect_equal(unname(classe.parameters["mu1"]), 0.04)
  expect_equal(unname(classe.parameters["q32"]), 0.04)
  expect_equal(unname(classe.parameters["mu2"]), 0.05)
  expect_equal(unname(classe.parameters["q31"]), 0.05)
  expect_equal(unname(classe.parameters["q13"]), 0.16)
  expect_equal(unname(classe.parameters["q23"]), 0.17)
})
