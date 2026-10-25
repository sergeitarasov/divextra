test_that("epoch targets stop exactly and timed epochs continue unsplit stems", {
  p <- setNames(rep(0,10),diversitree:::default.argnames.classe(2))
  p["lambda111"] <- 1
  set.seed(310)
  tb <- make.tree.classe.episodic(rbind(p,p,p),2,ntips=c(4,9,15))
  e <- attr(tb,"epochs")
  expect_equal(e$n.end,c(4,9,15))
  expect_true(all(e$reason=="ntips"))
  expect_true(all(diff(e$end)>0))
  expect_equal(sum(!tb$split & !tb$extinct),15)
  expect_true(attr(tb,"completed"))
  expect_equal(e$end.age,attr(tb,"t.total")-e$end)
  zero <- p*0
  tb <- make.tree.classe.episodic(rbind(zero,zero,zero),2,duration=c(1,2,3))
  expect_equal(attr(tb,"epochs")$end,c(1,3,6))
  expect_equal(tb$len,6)
  expect_true(attr(tb,"completed"))
})

test_that("both limits, downward targets and explicit early exits work", {
  p <- setNames(rep(0,10),diversitree:::default.argnames.classe(2))
  p["lambda111"] <- 1
  dead <- p*0
  dead["mu1"] <- 2
  set.seed(311)
  tb <- make.tree.classe.episodic(rbind(p,dead),2,ntips=c(5,2))
  expect_equal(attr(tb,"epochs")$n.end,c(5,2))
  tb <- make.tree.classe.episodic(rbind(p*0,p),2,ntips=c(5,4),duration=c(2,Inf))
  expect_equal(attr(tb,"epochs")$reason,c("duration","ntips"))
  tb <- make.tree.classe.episodic(rbind(dead,p),2,duration=c(100,1))
  expect_false(attr(tb,"completed"))
  expect_equal(attr(tb,"stop.reason"),"extinct")
  expect_equal(nrow(attr(tb,"epochs")),1)
  tb <- make.tree.classe.episodic(rbind(p*0,p),2,ntips=c(5,9))
  expect_equal(attr(tb,"stop.reason"),"stalled")
  tb <- make.tree.classe.episodic(rbind(p,p),2,ntips=c(20,30),max.taxa=4)
  expect_equal(sum(!tb$split & !tb$extinct),4)
  expect_equal(attr(tb,"stop.reason"),"max.taxa")
  tb <- make.tree.classe.episodic(rbind(p,p),2,ntips=c(20,30),max.events=2)
  expect_equal(attr(tb,"stop.reason"),"max.events")
  expect_equal(sum(attr(tb,"epochs")$events),2)
  tb <- make.tree.classe.episodic(rbind(p,p),2,ntips=c(1,3))
  expect_equal(attr(tb,"epochs")$duration[1],0)
})

test_that("histories, naming and input validation are preserved", {
  p <- pars.geosse2classe(c(1,1,1,.1,.1,1,1))
  set.seed(312)
  tb <- make.tree.classe.episodic(rbind(p,p),3,ntips=c(10,20),duration=c(10,10))
  h <- attr(tb,"hist")
  expect_true(nrow(h)>0)
  expect_true(all(h$t>=tb$start[h$idx]))
  expect_true(all(h$t<=tb$start[h$idx]+tb$len[h$idx]+1e-10))
  tree <- table2simmap(tb,c("A","B","AB"))
  expect_s3_class(tree,"simmap")
  expect_equal(unname(vapply(tree$maps,sum,numeric(1))),tree$edge.length)
  set.seed(313)
  a <- make.tree.classe.episodic(rbind(p,p),3,duration=c(2,3))
  set.seed(313)
  b <- make.tree.classe.episodic(rbind(p,p)[,27:1],3,duration=c(2,3))
  expect_identical(a,b)
  expect_error(make.tree.classe.episodic(rbind(p,p),3),"finite")
  expect_error(make.tree.classe.episodic(rbind(p,p),3,ntips=c(2,3,4)),"ntips")
  expect_error(make.tree.classe.episodic(rbind(p,p),3,ntips=1.5),"ntips")
  expect_error(make.tree.classe.episodic(rbind(p,p),3,duration=-1),"duration")
  expect_error(make.tree.classe.episodic(rbind(p,p),3,duration=1,x0=1.5),"x0")
})

test_that("pure-birth hitting times agree with the analytic expectation", {
  p <- setNames(rep(0,10),diversitree:::default.argnames.classe(2))
  p["lambda111"] <- 1
  recent <- p*2
  set.seed(314)
  times <- numeric(400)
  for(i in seq_along(times)) {
    tb <- make.tree.classe.episodic(rbind(p,recent),2,ntips=c(5,10))
    times[i] <- attr(tb,"t.total")
  }
  expectation <- sum(1/(1:4))+sum(1/(5:9))/2
  standard.error <- sqrt((sum(1/(1:4)^2)+sum(1/(5:9)^2)/4)/length(times))
  expect_lt(abs(mean(times)-expectation),5*standard.error)
})

test_that("epoch-specific cladogenetic and transition rates are used", {
  p <- setNames(rep(0,10),diversitree:::default.argnames.classe(2))
  p["lambda112"] <- 1
  q <- p*0
  q["q12"] <- 100
  set.seed(315)
  tb <- make.tree.classe.episodic(rbind(p,q),2,ntips=c(2,Inf),duration=c(Inf,1))
  first <- attr(tb,"info.epochs")[[1]]
  expect_equal(sort(first$state[!first$split]),c(1,2))
  expect_equal(tb$state[!tb$split],c(2,2))
  h <- attr(tb,"hist")
  expect_equal(nrow(h),1)
  expect_gt(h$t[1],attr(tb,"epochs")$end[1])
})
