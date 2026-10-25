# GeoSSE: backward likelihoods and an explicit forward history ----
# Run from the divextra package folder, or source this file in RStudio.
# No fitted parameters, marginal ASRs, fossils or incomplete sampling.
# This is our event-time conditional sampler, NOT the Freyman time-slice algorithm.
# The forward section below unfolds ONE deliberately short seeded realization.
# It is a teaching example, not a replacement for make.simmap.classe.td().

library(ape)
library(diversitree)
library(deSolve)
library(phytools)

# Tree, observed states, known rates and root prior ----
tree <- read.tree(text = "(tip_A:5,tip_B:5);")
plot(tree)
tip.states <- c(tip_A = 1, tip_B = 2)
# diversitree codes: 0 = AB, 1 = A, 2 = B.
# Its internal probability vectors are ordered AB, A, B (NOT A, B, AB).
states <- c("AB", "A", "B")
colors <- c(AB = "#CC79A7", A = "#0072B2", B = "#E69F00")
pars <- c(sA = .21, sB = .16, sAB = .09,
          xA = .13, xB = .11, dA = .07, dB = .05)
root.prior <- c(AB = .40, A = .35, B = .25)
ages <- seq(0, 5, by = .001)
print(pars)
plot(tree, main = "Observed tree: two living species")
tiplabels(c("A", "B"), bg = colors[c("A", "B")])

# Backward equations: use diversitree's own GeoSSE derivative ----
# deSolve needs a list, whereas this internal diversitree function returns a vector.
# ::: is intentional: this teaching example uses a non-exported function.
geosse.ode <- function(age, values, parameters) {
  change <- diversitree:::derivs.geosse(age, values, parameters)
  list(change)
}

# E_i(0) = 0: there are no unsampled living species.
# D for tip_A starts at (0,1,0); D for tip_B starts at (0,0,1).
initial.A <- c(E_AB = 0, E_A = 0, E_B = 0, D_AB = 0, D_A = 1, D_B = 0)
initial.B <- c(E_AB = 0, E_A = 0, E_B = 0, D_AB = 0, D_A = 0, D_B = 1)
backward.A <- ode(initial.A, ages, geosse.ode, pars, rtol = 1e-11, atol = 1e-13)
backward.B <- ode(initial.B, ages, geosse.ode, pars, rtol = 1e-11, atol = 1e-13)
backward.A <- as.data.frame(backward.A)
backward.B <- as.data.frame(backward.B)

# Inspect exact rows every 1 Ma. D is NOT a normalized state probability.
print(backward.A[match(0:5, backward.A$time), ])
print(backward.B[match(0:5, backward.B$time), ])
stopifnot(max(abs(backward.A[, 2:4] - backward.B[, 2:4])) < 1e-9)

# See how E and the two branch-specific D vectors change with age ----
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
matplot(ages, backward.A[, 2:4], type = "l", lty = 1, lwd = 2,
        col = colors[states], xlab = "Age (Ma)", ylab = "E", main = "Extinction probability")
legend("topleft", states, col = colors[states], lty = 1, bty = "n")
matplot(ages, backward.A[, 5:7], type = "l", lty = 1, lwd = 2,
        col = colors[states], xlab = "Age (Ma)", ylab = "D", main = "Branch ending in A")
matplot(ages, backward.B[, 5:7], type = "l", lty = 1, lwd = 2,
        col = colors[states], xlab = "Age (Ma)", ylab = "D", main = "Branch ending in B")
par(mfrow = c(1, 1))

# Combine daughters at the observed root speciation ----
left <- unlist(backward.A[nrow(backward.A), 5:7], use.names = FALSE)
right <- unlist(backward.B[nrow(backward.B), 5:7], use.names = FALSE)
names(left) <- names(right) <- states

# An unordered mixed-state split has two assignments to the labelled branches.
# Each assignment gets half the rate, matching diversitree's node convention.
root.events <- data.frame(
  parent = c("A", "B", "AB", "AB", "AB", "AB", "AB", "AB"),
  left = c("A", "B", "A", "AB", "B", "AB", "A", "B"),
  right = c("A", "B", "AB", "A", "AB", "B", "B", "A"),
  rate = c(pars["sA"], pars["sB"], pars["sA"] / 2, pars["sA"] / 2,
           pars["sB"] / 2, pars["sB"] / 2, pars["sAB"] / 2, pars["sAB"] / 2))
root.events$likelihood <- root.events$rate * left[root.events$left] * right[root.events$right]
root.D <- c(AB = sum(root.events$likelihood[root.events$parent == "AB"]),
            A = sum(root.events$likelihood[root.events$parent == "A"]),
            B = sum(root.events$likelihood[root.events$parent == "B"]))
loglik.manual <- log(sum(root.prior * root.D))

# Cross-check both the root combination and the full log likelihood ----
tip.root.values <- cbind(unlist(backward.A[nrow(backward.A), 2:7]),
                         unlist(backward.B[nrow(backward.B), 2:7]))
root.builtin <- diversitree:::initial.conditions.geosse(tip.root.values, pars, 5, NULL)
stopifnot(max(abs(root.D - root.builtin[4:6])) < 1e-10)
lik <- make.geosse(tree, tip.states, sampling.f = c(1, 1, 1), strict = FALSE,
                  control = list(rtol = 1e-10, atol = 1e-12))
# FALSE avoids an extra crown-survival normalization; the tree/tips still condition maps.
loglik.diversitree <- lik(pars, condition.surv = FALSE, root = ROOT.GIVEN, root.p = root.prior)
print(c(manual = loglik.manual, diversitree = loglik.diversitree))
stopifnot(abs(loglik.manual - loglik.diversitree) < 1e-7)
root.posterior <- root.prior * root.D / sum(root.prior * root.D)
print(root.posterior)

# Interpolation: query the backward curves at any event age ----
E.AB <- approxfun(ages, backward.A$E_AB)
E.A <- approxfun(ages, backward.A$E_A)
E.B <- approxfun(ages, backward.A$E_B)
D.left.AB <- approxfun(ages, backward.A$D_AB)
D.left.A <- approxfun(ages, backward.A$D_A)
D.left.B <- approxfun(ages, backward.A$D_B)
D.right.AB <- approxfun(ages, backward.B$D_AB)
D.right.A <- approxfun(ages, backward.B$D_A)
D.right.B <- approxfun(ages, backward.B$D_B)
E <- function(age) c(AB = E.AB(age), A = E.A(age), B = E.B(age))
D.left <- function(age) c(AB = D.left.AB(age), A = D.left.A(age), B = D.left.B(age))
D.right <- function(age) c(AB = D.right.AB(age), A = D.right.A(age), B = D.right.B(age))

# All original GeoSSE events: inspect this table before going forward ----
events <- data.frame(
  parent = c("A", "A", "A", "B", "B", "B", "AB", "AB", "AB", "AB", "AB"),
  type = c("split", "change", "death", "split", "change", "death",
           "split", "split", "split", "change", "change"),
  daughter1 = c("A", "AB", NA, "B", "AB", NA, "A", "B", "A", "B", "A"),
  daughter2 = c("A", NA, NA, "B", NA, NA, "AB", "AB", "B", NA, NA),
  rate = unname(pars[c("sA", "dA", "xA", "sB", "dB", "xB",
                      "sA", "sB", "sAB", "xA", "xB")]))
print(events)

# One forward step: the only reusable sampling helper ----
# h = D.left or D.right for an observed branch; h = E for an extinct side lineage.
# This function samples just ONE event, never a whole branch or tree.
# Its returned candidates table exposes every weight and normalized probability.
next.event <- function(age, state, h, extinct = FALSE) {
  possible <- events[events$parent == state, ]
  total.rate <- sum(possible$rate)
  h.start <- h(age)[state]
  stopifnot(age > 0, h.start > 0)

  # P(no event until age u | required descendants) = exp(-r*(age-u))*h(u)/h(age).
  # Draw an Exp(1) threshold in cumulative-hazard space, not Exp(r) waiting time.
  threshold <- rexp(1)
  hazard <- function(u) total.rate * (age - u) + log(h.start) - log(h(u)[state])
  p.no.event <- exp(-hazard(0))
  if (threshold >= hazard(0)) {
    return(list(type = "end", age = 0, state = state, p.no.event = p.no.event,
                threshold = threshold, candidates = NULL))
  }
  # uniroot can handle the +Inf hazard at u=0 when h(0)=0.
  equation <- function(u) hazard(u) - threshold
  event.age <- uniroot(equation, c(0, age), tol = 1e-10)$root
  ee <- E(event.age)
  hh <- h(event.age)
  splits <- possible[possible$type == "split", ]
  changes <- possible[possible$type == "change", ]

  if (!extinct) {
    # Either daughter can carry the unique observed descendant.
    # Equal-state splits also have two possibilities: total weight 2*lambda*D*E.
    first <- data.frame(type = "split", state = splits$daughter1,
                        side = splits$daughter2,
                        weight = splits$rate * hh[splits$daughter1] * ee[splits$daughter2])
    second <- data.frame(type = "split", state = splits$daughter2,
                         side = splits$daughter1,
                         weight = splits$rate * hh[splits$daughter2] * ee[splits$daughter1])
    change <- data.frame(type = "change", state = changes$daughter1, side = NA_character_,
                         weight = changes$rate * hh[changes$daughter1])
    candidates <- rbind(first, second, change)
  } else {
    # Both daughters must go extinct, so there is no factor of two here.
    first <- data.frame(type = "split", state = splits$daughter1, side = splits$daughter2,
                        weight = splits$rate * ee[splits$daughter1] * ee[splits$daughter2])
    change <- data.frame(type = "change", state = changes$daughter1, side = NA_character_,
                         weight = changes$rate * ee[changes$daughter1])
    death <- data.frame(type = "death", state = NA_character_, side = NA_character_,
                        weight = sum(possible$rate[possible$type == "death"]))
    candidates <- rbind(first, change, death)
  }
  candidates$probability <- candidates$weight / sum(candidates$weight)
  selected <- sample.int(nrow(candidates), 1, prob = candidates$probability)
  list(type = candidates$type[selected], age = event.age, state = candidates$state[selected],
       side = candidates$side[selected], p.no.event = p.no.event, threshold = threshold,
       selected = selected, candidates = candidates)
}

# Forward: sample the root and its two daughter states ----
# Seed is selected for a short, illustrative history, NOT for a statistical analysis.
set.seed(70)
root.state <- sample(states, 1, prob = root.posterior)
root.choices <- root.events[root.events$parent == root.state, ]
root.choices$probability <- root.choices$likelihood / sum(root.choices$likelihood)
root.choice <- sample.int(nrow(root.choices), 1, prob = root.choices$probability)
left.state <- root.choices$left[root.choice]
right.state <- root.choices$right[root.choice]
print(root.state)
print(root.choices)
print(c(left = left.state, right = right.state))

# Left branch, step 1: draw the next age and event ----
# Initially AB, but its observed endpoint is A. No-event probability must be zero.
left.1 <- next.event(age = 5, state = left.state, h = D.left)
print(left.1)
# With seed 70: AB -> A, by local extinction in B, at about 4.535 Ma.
stopifnot(root.state == "AB", left.state == "AB", right.state == "B",
          left.1$type == "change", left.1$state == "A")

# Left branch, step 2: now start from that newly sampled state and age ----
left.2 <- next.event(age = left.1$age, state = left.1$state, h = D.left)
print(left.2)
# A -> A + A at about 2.938 Ma. One A continues to tip_A; store the other A.
# The two identical split rows represent the two possible surviving daughters.
stopifnot(left.2$type == "split", left.2$state == "A", left.2$side == "A")
side.birth.age <- left.2$age
side.birth.state <- left.2$side

# Left branch, step 3: does anything else happen before the present? ----
left.3 <- next.event(age = left.2$age, state = left.2$state, h = D.left)
print(left.3)
# This time the Exp(1) threshold exceeds the cumulative hazard to the endpoint.
# The draw is "end": retain A all the way to the present. No forced endpoint change.
stopifnot(left.3$type == "end", left.3$state == "A")

# Right branch: independently simulate conditional on its observed B endpoint ----
right.1 <- next.event(age = 5, state = right.state, h = D.right)
print(right.1)
# In this realization, B remains B for the entire branch.
stopifnot(right.1$type == "end", right.1$state == "B")

# Expand the stored side lineage, using E rather than D ----
side.1 <- next.event(age = side.birth.age, state = side.birth.state, h = E, extinct = TRUE)
print(side.1)
# This A daughter dies at about 1.432 Ma. Other seeds can produce whole side clades.
# No artificial killing and no rejection of surviving proposals were used.
stopifnot(side.1$type == "death", side.1$age > 0)

# Assemble this particular complete history into a simmap tree ----
# Tips: 1 = observed A; 2 = observed B; 3 = extinct A.
# Nodes: 4 = observed root; 5 = the hidden speciation.
mapped.tree <- list(
  edge = matrix(c(4, 5, 5, 1, 5, 3, 4, 2), ncol = 2, byrow = TRUE),
  edge.length = c(5 - side.birth.age, side.birth.age,
                  side.birth.age - side.1$age, 5),
  tip.label = c("tip_A", "tip_B", "extinct_A"),
  Nnode = 2L)
mapped.tree$maps <- list(
  c(AB = 5 - left.1$age, A = left.1$age - side.birth.age),
  c(A = side.birth.age),
  c(A = side.birth.age - side.1$age),
  c(B = 5))
mapped.tree$mapped.edge <- matrix(0, nrow = 4, ncol = 3,
                                  dimnames = list(NULL, states))
mapped.tree$mapped.edge[1, ] <- c(5 - left.1$age, left.1$age - side.birth.age, 0)
mapped.tree$mapped.edge[2, "A"] <- side.birth.age
mapped.tree$mapped.edge[3, "A"] <- side.birth.age - side.1$age
mapped.tree$mapped.edge[4, "B"] <- 5
class(mapped.tree) <- c("simmap", "phylo")
attr(mapped.tree, "order") <- "cladewise"

# Separate event record preserves the ancestral state at the observed root.
# This manually built simmap is for plotting; it is NOT a full DivExtra result object.
history <- data.frame(
  age = c(5, left.1$age, side.birth.age, side.1$age, 0),
  event = c("Root: AB -> AB + B", "Left: AB -> A", "Left: A -> A + A",
            "Side: A dies", "Observed tips: A and B"))
print(history)

# Structural checks: map lengths, extinction age and original backbone ----
stopifnot(max(abs(rowSums(mapped.tree$mapped.edge) - mapped.tree$edge.length)) < 1e-10)
stopifnot(abs(sum(mapped.tree$maps[[1]]) - mapped.tree$edge.length[1]) < 1e-10,
          abs(sum(mapped.tree$maps[[2]]) - mapped.tree$edge.length[2]) < 1e-10,
          abs(sum(mapped.tree$maps[[3]]) - mapped.tree$edge.length[3]) < 1e-10,
          abs(sum(mapped.tree$maps[[4]]) - mapped.tree$edge.length[4]) < 1e-10)
depths <- node.depth.edgelength(mapped.tree)
stopifnot(abs(depths[1] - 5) < 1e-10, abs(depths[2] - 5) < 1e-10,
          abs(depths[3] - (5 - side.1$age)) < 1e-10)
observed.tree <- mapped.tree
class(observed.tree) <- "phylo"
observed.tree$maps <- observed.tree$mapped.edge <- NULL
observed.tree <- drop.tip(observed.tree, "extinct_A")
stopifnot(isTRUE(all.equal.phylo(tree, observed.tree, use.edge.length = TRUE)))

# Plot the actual sampled history ----
plotSimmap(mapped.tree, colors = colors, split.vertical = TRUE, ftype = "i",
           fsize = 1, lwd = 3, xlim = c(0, 6.5), ylim = c(.7, 3.7),
           mar = c(4, 3, 3, 1))
axis(1, at = 0:5, labels = 5:0)
mtext("Age (Ma)", side = 1, line = 2.5)
nodelabels("AB", node = 4, bg = colors["AB"], cex = .8)
legend("topright", legend = states, fill = colors[states], bty = "n", cex = .8)
title("GeoSSE: two living tips and one extinct tip")

# Interpretation and limits ----
# The single extinction is one selected realization, NOT an expected number.
# All rates differ, but are constant through time here; there are no epochs to cross.
# The .001 Ma grid interpolates probabilities; event times are continuous, not snapped to it.
# Refine ages to check interpolation error. Production DivExtra uses adaptive curves.
# This code uses D/E-weighted event times, not Freyman's short-interval forward ODE pass.
# It demonstrates complete sampling: E(0)=0. It does not infer missing living species.
# To change parameters/seed freely, use the package simulator: a new history can need
# arbitrarily many steps/side descendants, so it cannot generally be unfolded without loops.
