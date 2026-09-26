tree_cond_args <- function() {
  list(tree=ape::read.tree(text="(a:2,b:2);"),
    pars=c(lambda111.1=.5,mu1.1=.05),tip.weights=c(a=1,b=1),
    k=1,state.labels="A",backbone.method="tree-conditioned",root.prior=1)
}

test_that("tree-conditioned hidden births follow the analytical Poisson law", {
  args<-tree_cond_args()
  sims<-do.call(make.simmap.classe.td,c(args,list(nsim=1000,seed=1401)))
  E<-function(t).05*(-expm1(-.45*t))/(.5-.05*exp(-.45*t))
  expected<-2*integrate(function(t)2*.5*E(t),0,2)$value
  counts<-vapply(sims,function(s)sum(s$events$event=="hidden_speciation"),0)
  expect_lt(abs(mean(counts)-expected),5*sqrt(expected/length(sims)))
  expect_lt(abs(mean(counts==0)-exp(-expected)),.05)
  expect_true(all(vapply(sims,function(s)all(s$branch.attempts==1L),TRUE)))
  expect_true(all(vapply(sims,function(s)all(s$tips$fate[s$tips$generated]=="extinct"),TRUE)))
  for(s in sims[1:10]) {
    expect_equal(vapply(s$tree$maps,sum,0),s$tree$edge.length)
    expect_identical(s$method,"tree-conditioned-extinction-only")
  }
})

tree_cond_two_state <- function() {
  c(t.1=.7,lambda111.1=.3,lambda112.1=.15,lambda122.1=.08,
    lambda211.1=.05,lambda212.1=.12,lambda222.1=.4,
    mu1.1=.2,mu2.1=.6,q12.1=.25,q21.1=.1,
    lambda111.2=.4,lambda112.2=.2,lambda122.2=.03,
    lambda211.2=.1,lambda212.2=.15,lambda222.2=.2,
    mu1.2=.5,mu2.2=.1,q12.2=.1,q21.2=.3)
}

test_that("likelihood and sampled ancestral states match diversitree across epochs", {
  tree<-ape::read.tree(text="((a:.5,b:.5):.5,c:1);")
  pars<-tree_cond_two_state()
  states<-c(a=1,b=2,c=1)
  prior<-c(.4,.6)
  lik<-make.classe.td(tree,states,k=2,n.epoch=2,sampling.f=c(1,1),
    control=list(tol=1e-10))
  exact<-lik(pars,root=diversitree::ROOT.GIVEN,root.p=prior,condition.surv=FALSE)
  asr<-asr.marginal.classe(lik,pars,root=diversitree::ROOT.GIVEN,
    root.p=prior,condition.surv=FALSE)
  sims<-make.simmap.classe.td(tree,pars,tip.weights=states,k=2,
    state.labels=c("A","B"),backbone.method="tree-conditioned",root.prior=prior,
    nsim=800,seed=2841)
  expect_equal(sims[[1]]$root.loglik,as.numeric(exact),tolerance=2e-4)
  expect_equal(sims[[1]]$root.p,as.numeric(asr[,1]),tolerance=1e-4)
  for(node in c(4,5)) {
    p<-mean(vapply(sims,function(s) {
      e<-s$events[s$events$event=="observed_speciation" & s$events$node==paste0("b",node),]
      e$from==1
    },TRUE))
    expect_equal(p,as.numeric(asr[1,node-length(tree$tip.label)]),tolerance=.07)
  }
  expect_true(all(vapply(sims,function(s)
    all(s$tips$state[!s$tips$generated]==states[s$tips$label[!s$tips$generated]]),TRUE)))
})

test_that("zero extinction gives CTMC bridges and no hidden side births", {
  pars<-c(lambda111.1=.4,lambda112.1=0,lambda122.1=0,
    lambda211.1=0,lambda212.1=0,lambda222.1=.4,
    mu1.1=0,mu2.1=0,q12.1=.3,q21.1=.1)
  P<-expm::expm(matrix(c(-.3,.3,.1,-.1),2,byrow=TRUE))
  expected<-c(.5,.5)*P[,1]*P[,2];expected<-expected/sum(expected)
  sims<-make.simmap.classe.td(ape::read.tree(text="(a:1,b:1);"),pars,
    tip.weights=c(a=1,b=2),k=2,state.labels=c("A","B"),
    backbone.method="tree-conditioned",root.prior=c(.5,.5),nsim=500,seed=54)
  expect_equal(sims[[1]]$root.p,unname(expected),tolerance=1e-5)
  expect_true(all(vapply(sims,function(s)sum(s$tips$generated)==0,TRUE)))
  expect_equal(mean(vapply(sims,function(s)s$root.state==1,TRUE)),
               expected[1],tolerance=.07)
})

test_that("new mode refuses ambiguous conditioning inputs", {
  args<-tree_cond_args()
  expect_error(do.call(make.simmap.classe.td,c(args,list(node.probs=matrix(1,1,1)))),"node.probs")
  expect_error(do.call(make.simmap.classe.td,c(args,list(side.lineage.mode="all"))),"extinct_clades")
  args$root.prior<-NULL
  expect_error(do.call(make.simmap.classe.td,args),"root.prior")
})

test_that("selected side expansion leaves the conditioned backbone unchanged", {
  args<-tree_cond_args()
  args$pars<-c(lambda111.1=.5,mu1.1=1)
  full<-do.call(make.simmap.classe.td,c(args,list(seed=25)))
  focus<-do.call(make.simmap.classe.td,c(args,list(seed=25,side.lineage.tips="a")))
  expect_equal(prune.unobserved.simmap(full)$tree,prune.unobserved.simmap(focus)$tree)
  expect_equal(full$root.p,focus$root.p)
  expect_setequal(focus$tips$label[!focus$tips$generated],c("a","b"))
})
