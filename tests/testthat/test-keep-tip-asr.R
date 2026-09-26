test_that("pruning preserves node identity, ASR and distances", {
  tree <- ape::read.tree(text = "(((a:1,b:1):1,c:2):1,(d:2,e:2):1);")
  tree$node.label <- paste0("original", 6:9)
  asr <- rbind(A = (1:4)/5, B = 1-(1:4)/5)
  tips <- setNames(c(1,2,1,2,1), tree$tip.label)
  x <- extract.clade.asr(tree, 6, asr, tips)
  y <- keep.tip.asr(x, c("a", "c", "e"))
  expect_equal(ape::cophenetic.phylo(y$tree),
               ape::cophenetic.phylo(tree)[y$tree$tip.label, y$tree$tip.label])
  for (i in seq_len(y$tree$Nnode)) {
    descendants <- ape::extract.clade(y$tree, ape::Ntip(y$tree)+i)$tip.label
    original <- ape::getMRCA(tree, descendants)
    expect_equal(y$node.map$original.node[i], original)
    expect_equal(unname(y$node.probs[, i]), unname(asr[, original-5]))
    expect_equal(y$tree$node.label[i], tree$node.label[original-5])
  }
  expect_equal(y$tip.states, tips[y$tree$tip.label])
  expect_equal(colnames(y$node.probs), as.character(4:5))
  expect_equal(keep.tip.asr(y, c("a", "c")), keep.tip.asr(x, c("a", "c")))
  expect_equal(keep.tip.asr(x, tree$tip.label), x)
  expect_equal(keep.tip.asr(x, c("a", "b"))$node.map$original.node, 8L)
  x$tip.states <- rbind(A = as.numeric(tips==1), B = as.numeric(tips==2))
  colnames(x$tip.states) <- names(tips)
  expect_equal(keep.tip.asr(x, c("a", "c"))$tip.states,
               x$tip.states[, c("a", "c")])
  expect_error(keep.tip.asr(x, c("a", "typo")), "Unknown tip")
  expect_error(keep.tip.asr(x, "a"), "at least two")
  expect_error(keep.tip.asr(x, c("a", "a")), "unique")
})
