# Run after devtools::load_all('/Users/taravser/GitHub/divextra').
# Three-state, three-tip illustration; forward simulation uses the same rates.
library(divextra)
pars.ge <- matrix(c(
  .1, .1, .1, 0, 0, .1, .1,
  .3, .3, .3, 0, 0, .1, .1
), 2, 7, byrow = TRUE)
colnames(pars.ge) <- diversitree:::default.argnames.geosse()
pars.cl <- t(apply(pars.ge, 1, pars.geosse2classe))
set.seed(123)
tb <- make.tree.classe.td(pars.cl, k = 3, max.t1 = 10,
                         max.t2 = 15, x0 = 2, single.lineage = TRUE)

# Forward rows run old -> young. Inference parameters run young -> old.
to.inference <- function(x) c(
  t.1 = 5,
  stats::setNames(x[2, ], paste0(colnames(x), '.1')),
  stats::setNames(x[1, ], paste0(colnames(x), '.2'))
)
tree <- ape::read.tree(text = '((a:3,b:3):4,c:7);')
states <- c(a = 1L, b = 2L, c = 2L)
labels <- c('A', 'B', 'AB')
root.p <- c(A = 0, B = 1, AB = 0)
tip.weights <- matrix(0, 3, 3, dimnames = list(labels, names(states)))
tip.weights[cbind(states, seq_along(states))] <- 1

run.example <- function(forward.pars, seed) {
  pars <- to.inference(forward.pars)
  lik <- make.classe.td(tree, states, k = 3, n.epoch = 2, strict = FALSE)
  node.probs <- asr.marginal.classe(
    lik, pars, root = diversitree::ROOT.GIVEN, root.p = unname(root.p)
  )
  dimnames(node.probs) <- list(labels, as.character(4:5))
  make.simmap.classe.td(
    tree, pars, tip.weights, node.probs, k = 3,
    state.labels = labels, nsim = 3, seed = seed
  )
}
zero.extinction <- run.example(pars.cl, 123)
positive.ge <- pars.ge
positive.ge[, c('xA', 'xB')] <- .15
positive.cl <- t(apply(positive.ge, 1, pars.geosse2classe))
positive.extinction <- run.example(positive.cl, 124)
counts <- function(sims) t(vapply(sims, function(sim) c(
  observed = sum(sim$tips$fate == 'observed'),
  extinct = sum(sim$tips$fate == 'extinct'),
  extra_extant = sum(sim$tips$fate == 'extra_extant')
), numeric(3)))
print(counts(zero.extinction))
print(counts(positive.extinction))

colors <- c(A = '#D55E00', B = '#0072B2', AB = '#009E73')
if (interactive() && requireNamespace('phytools', quietly = TRUE)) {
  phytools::plotSimmap(positive.extinction[[1]]$tree,
                      colors = colors, split.vertical = TRUE, ftype = 'i')
}
