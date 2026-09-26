transition_fixture <- function(clado=FALSE, zero=FALSE) {
  edge <- function(parent,child,start,end,map,generated=FALSE,b=1L)
    list(parent=parent,child=child,start.age=start,end.age=end,length=start-end,
         map=map,generated=generated,backbone.edge=if(generated) NA_integer_ else b)
  edges <- list(
    edge("r","h",5,3,if(clado || zero) c(A=2) else c(A=1,A.S=1)),
    edge("h","a",3,0,if(zero) c(A=3) else if(clado) c(A.S=3) else c(A=1.5,A.S=1.5)),
    edge("r","b",5,0,if(clado) c(A.S=5) else c(A=5),b=2L),
    edge("h","g",3,2.5,c(A=.25,A.S=.25),TRUE),
    edge("g","e",2.5,1,c(A.S=1.5),TRUE),
    edge("g","f",2.5,0,c(A.S=2.5),TRUE))
  tip <- function(node,age,generated,state) data.frame(node=node,label=node,age=age,
    state=if(state=="A") 1L else 2L,generated=generated,state.label=state,
    fate=if(!generated) "observed" else if(age>0) "extinct" else "extra_extant")
  sim <- .classe_td_graph_to_simmap(edges,list(a=tip("a",0,FALSE,if(zero) "A" else "A.S"),
    b=tip("b",0,FALSE,if(clado) "A.S" else "A"),e=tip("e",1,TRUE,"A.S"),
    f=tip("f",0,TRUE,"A.S")),c("A","A.S"))
  sim$events <- data.frame(event.id=1:4,
    event=c("observed_speciation","hidden_speciation","speciation","transition"),
    age=c(5,3,2.5,2.75),node=c("r","h","g","h"),
    from.label=c("A",if(clado || zero) "A" else "A.S","A.S","A"),
    backbone.edge=c(NA,1L,NA,NA))
  sim
}

test_that("pruning removes generated extant and extinct tips without changing maps", {
  sim <- transition_fixture()
  p <- prune.unobserved.simmap(sim)
  expect_setequal(unname(p$tree$tip.label),c("a","b"))
  expect_false(any(p$tips$generated))
  expect_false(any(p$lineages$generated))
  expect_equal(p$lineages$map,sim$lineages$map[!sim$lineages$generated])
  expect_equal(vapply(p$tree$maps,sum,0),p$tree$edge.length)
  expect_false(any(p$events$event=="transition")) # side event at surviving h
  expect_equal(prune.unobserved.simmap(p),p)
  focal <- prune.unobserved.simmap(keep.lineage.simmap(sim,"a"))
  expect_equal(unname(focal$tree$tip.label),"a")
  expect_equal(max(focal$lineages$start.age),5)
})

test_that("anagenetic arrival ages and zero maps have the correct denominators", {
  sims <- lapply(list(transition_fixture(),transition_fixture(zero=TRUE)),prune.unobserved.simmap)
  x <- transitions.through.time.simmap(sims,0:5,to="A.S")
  expect_equal(x$events$age,c(4,1.5))
  expect_true(all(x$events$event.type=="anagenetic"))
  expect_equal(x$per.map$n.transitions,c(2L,0L))
  expect_equal(unname(x$counts[,1]),c(0,1,0,0,1))
  expect_equal(x$overall$probability,.5)
  expect_equal(x$summary$p.any,c(0,.5,0,0,.5))
  y <- transitions.through.time.simmap(sims,0:5,to="A.S",quantile.type=7)
  expect_equal(y$summary$lower,c(0,.025,0,0,.025))
  expect_equal(y$summary$upper,c(0,.975,0,0,.975))
  expect_equal(y$age.summary$median,c(4,1.5,2.75))
  expect_equal(y$overall,x$overall)
  expect_equal(y$summary$prob.lower,x$summary$prob.lower)
  expect_error(transitions.through.time.simmap(sims,0:5,to="A.S",quantile.type=10),"quantile.type")
  expect_equal(x$age.summary$median,c(4,1.5,1.5))
  expect_equal(x$event.age.distribution$probability,c(0,.5,0,0,.5))
  expect_equal(x$first.age.distribution$probability,c(0,0,0,0,1))
  expect_true(is.na(x$per.map$first.age[2]))
  expect_equal(x$count.distribution$probability,c(.5,0,.5))
  ci <- stats::binom.test(1,2)$conf.int
  expect_equal(c(x$overall$lower,x$overall$upper),as.numeric(ci))
  expect_warning(transitions.through.time.simmap(sims,0:2,to="A.S"),"outside")
  expect_error(transitions.through.time.simmap(transition_fixture(),0:5,to="A.S"),"prune.unobserved")
  expect_error(transitions.through.time.simmap(sims,0:5,to="typo"),"Unknown")
})

test_that("cladogenetic arrivals use only retained daughters including root changes", {
  sim <- transition_fixture(clado=TRUE)
  p <- prune.unobserved.simmap(sim)
  x <- transitions.through.time.simmap(p,0:5,to="A.S")
  expect_equal(x$events$age,c(5,3))
  expect_true(all(x$events$event.type=="cladogenetic"))
  expect_equal(unname(x$counts[,1]),c(0,0,0,1,1))
  focal <- prune.unobserved.simmap(keep.lineage.simmap(sim,"a"))
  f <- transitions.through.time.simmap(focal,0:5,to="A.S")
  expect_equal(f$events$age,3)
  expect_equal(transitions.through.time.simmap(p,0:5,to="A.S",
    event.types="anagenetic")$overall$probability,0)
})

test_that("no-event age distributions are undefined and plotting works", {
  p <- prune.unobserved.simmap(transition_fixture(zero=TRUE))
  x <- transitions.through.time.simmap(p,0:5,to="A.S")
  expect_equal(x$overall$probability,0)
  expect_true(x$overall$upper>0)
  expect_true(all(is.na(x$age.summary$median)))
  expect_true(all(is.na(x$first.age.distribution$probability)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(x))
  expect_invisible(plot(x,type="counts"))
  expect_invisible(plot(x,type="first.age"))
})

test_that("existing root occupancy is flagged but is not a migration", {
  sim <- transition_fixture(zero=TRUE)
  sim$events$from.label[sim$events$node=="r"] <- "A.S"
  p <- prune.unobserved.simmap(sim)
  x <- transitions.through.time.simmap(p,0:5,to="A.S")
  expect_true(x$per.map$root.in.destination)
  expect_equal(x$overall$n.root.in.destination,1)
  expect_equal(x$overall$probability,0)
  expect_true(is.na(x$per.map$first.age))
})
