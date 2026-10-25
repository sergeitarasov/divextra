test_that("rounded probes at nonzero epoch boundaries use actual weights", {
  schedule <- .classe_td_schedule(c(t.1=9, lambda111.1=0, mu1.1=0,
    lambda111.2=0, mu1.2=.0155), 1, "A")
  ages <- c(0, 8, 9, 9+10^seq(-12,-2,length.out=60), 10, 22.73)
  expected <- -expm1(-.0155*pmax(ages-9,0))
  for (control in list(list(),list(rtol=1e-12,atol=1e-15))) {
    cache <- .classe_td_extinction_cache(schedule,22.73,control)
    actual <- vapply(ages,function(t)cache$E(t)[1],0)
    expect_equal(actual,expected,tolerance=1e-6)
    near <- ages>9 & ages<9.01
    expect_lt(max(abs(log(actual[near])-log(expected[near]))),1e-6)
    for (d in cache$diagnostics) {
      expect_lte(d$max.absolute.error,1e-6)
      expect_lte(d$max.log.error,1e-6)
    }
  }
})
