# Simulating episodic ClaSSE histories on a fixed reconstructed tree.
#
# The event engine follows the Gillespie construction used in diversitree's
# make.tree.classe(), but keeps absolute ages so that a lineage can cross any
# number of make.classe.td() epochs.  The implementation is owned here because
# the diversitree conversion helpers discard the founder lineage and are not
# suitable for grafting complete side trees onto a supplied backbone.

.classe_td_state_vector <- function(x, ids, k, state.labels, what) {
  if (is.matrix(x)) {
    if (nrow(x) != k)
      stop(what, " probability matrix must have k rows")
    x <- apply(x, 2L, which.max)
  }
  if (is.null(names(x))) {
    if (length(x) != length(ids))
      stop(what, " must be named or have one value per requested node/tip")
    names(x) <- as.character(ids)
  }
  key <- as.character(ids)
  if (!all(key %in% names(x)))
    stop(what, " is missing: ", paste(setdiff(key, names(x)), collapse = ", "))
  x <- x[key]
  if (is.character(x)) {
    x <- match(x, state.labels)
  }
  x <- as.integer(x)
  if (anyNA(x) || any(x < 1L | x > k))
    stop(what, " must contain states in 1:k or matching state labels")
  stats::setNames(x, key)
}

.classe_td_state_data <- function(x, ids, k, state.labels, what) {
  key <- as.character(ids)
  if (!is.matrix(x)) {
    states <- .classe_td_state_vector(x, ids, k, state.labels, what)
    probabilities <- matrix(
      0, nrow = k, ncol = length(key),
      dimnames = list(state.labels, key)
    )
    probabilities[cbind(unname(states), seq_along(states))] <- 1
    return(list(states = states, probabilities = probabilities, soft = FALSE))
  }

  if (nrow(x) != k)
    stop(what, " probability matrix must have k rows")
  storage.mode(x) <- "double"
  if (any(!is.finite(x)) || any(x < 0))
    stop(what, " probabilities must be finite and non-negative")
  if (!is.null(rownames(x)) && all(state.labels %in% rownames(x)))
    x <- x[state.labels, , drop = FALSE]
  if (is.null(colnames(x))) {
    if (ncol(x) != length(key))
      stop(what, " probability matrix must be named or have one column per requested node/tip")
    colnames(x) <- key
  }
  if (!all(key %in% colnames(x)))
    stop(what, " is missing: ", paste(setdiff(key, colnames(x)), collapse = ", "))
  x <- x[, key, drop = FALSE]
  totals <- colSums(x)
  if (any(totals <= 0))
    stop(what, " probability columns must have positive sums")
  x <- sweep(x, 2L, totals, "/")
  rownames(x) <- state.labels
  states <- stats::setNames(max.col(t(x), ties.method = "first"), key)
  list(states = states, probabilities = x, soft = TRUE)
}

.classe_td_schedule <- function(pars, k, state.labels = as.character(seq_len(k))) {
  if (is.null(names(pars)))
    stop("pars must be a named full parameter vector from make.classe.td()")
  time.idx <- grep("^t\\.[0-9]+$", names(pars))
  n.epoch <- length(time.idx) + 1L
  boundaries <- if (length(time.idx)) as.numeric(pars[time.idx]) else numeric()
  if (any(!is.finite(boundaries)) || is.unsorted(boundaries, strictly = TRUE))
    stop("Epoch boundaries t.1, t.2, ... must be finite and strictly increasing")

  values <- if (length(time.idx)) pars[-time.idx] else pars
  if (length(values) %% n.epoch != 0L)
    stop("Parameter count is not divisible by the inferred number of epochs")
  n.par <- length(values) %/% n.epoch
  blocks <- split(values, rep(seq_len(n.epoch), each = n.par))

  arrays <- lapply(seq_along(blocks), function(i) {
    block <- blocks[[i]]
    suffix <- paste0("\\.", i, "$")
    if (!all(grepl(suffix, names(block))))
      stop("Epoch ", i, " parameter names do not all end in .", i)
    names(block) <- sub(suffix, "", names(block))
    if (any(!is.finite(block)) || any(block < 0))
      stop("All ClaSSE rates must be finite and non-negative")
    pars_to_arrays(block, k, state.labels)
  })

  list(k = k, n.epoch = n.epoch, boundaries = boundaries,
       arrays = arrays, state.labels = state.labels, pars = pars)
}

.classe_td_epoch <- function(age, boundaries) {
  1L + sum(age > boundaries)
}

.classe_td_lower_boundary <- function(age, boundaries) {
  lower <- boundaries[boundaries < age]
  if (length(lower)) max(lower) else 0
}

.classe_td_transition_matrix <- function(start.age, stop.age, schedule) {
  k <- schedule$k
  ans <- diag(k)
  age <- start.age
  while (age > stop.age + 1e-12) {
    epoch <- .classe_td_epoch(age, schedule$boundaries)
    boundary <- max(
      stop.age,
      .classe_td_lower_boundary(age, schedule$boundaries)
    )
    generator <- schedule$arrays[[epoch]]$Q
    diag(generator) <- -rowSums(generator)
    ans <- ans %*% expm::expm(generator * (age - boundary))
    age <- boundary
  }
  ans
}

.classe_td_observed_pair <- function(array, parent.state, child.scores) {
  lambda <- array$lam.tensor[[parent.state]]
  ij <- which(lambda > 0, arr.ind = TRUE)
  if (!nrow(ij))
    return(NULL)
  candidates <- vector("list", nrow(ij))
  weights <- vector("list", nrow(ij))
  for (z in seq_len(nrow(ij))) {
    i <- ij[z, 1]
    j <- ij[z, 2]
    rate <- lambda[i, j]
    if (i == j) {
      candidates[[z]] <- matrix(c(i, j), nrow = 1L)
      weights[[z]] <- rate * child.scores[[1]][i] * child.scores[[2]][j]
    } else {
      candidates[[z]] <- rbind(c(i, j), c(j, i))
      weights[[z]] <- c(
        rate * 0.5 * child.scores[[1]][i] * child.scores[[2]][j],
        rate * 0.5 * child.scores[[1]][j] * child.scores[[2]][i]
      )
    }
  }
  candidates <- do.call(rbind, candidates)
  weights <- unlist(weights, use.names = FALSE)
  if (!length(weights) || sum(weights) <= 0)
    return(NULL)
  candidates[sample.int(nrow(candidates), 1L, prob = weights), , drop = TRUE]
}

.classe_td_observed_weight <- function(array, parent.state, child.scores) {
  lambda <- array$lam.tensor[[parent.state]]
  ij <- which(lambda > 0, arr.ind = TRUE)
  if (!nrow(ij))
    return(0)
  sum(vapply(seq_len(nrow(ij)), function(z) {
    i <- ij[z, 1]
    j <- ij[z, 2]
    rate <- lambda[i, j]
    if (i == j)
      return(rate * child.scores[[1]][i] * child.scores[[2]][j])
    rate * 0.5 * (
      child.scores[[1]][i] * child.scores[[2]][j] +
        child.scores[[1]][j] * child.scores[[2]][i]
    )
  }, numeric(1)))
}

.append_map <- function(map, state, duration, state.labels) {
  if (duration <= 0)
    return(map)
  nm <- state.labels[state]
  if (length(map) && identical(names(map)[length(map)], nm)) {
    map[length(map)] <- map[length(map)] + duration
  } else {
    map <- c(map, stats::setNames(duration, nm))
  }
  map
}

.sample_lambda_pair <- function(array, parent.state) {
  lambda <- array$lam.tensor[[parent.state]]
  ij <- which(lambda > 0, arr.ind = TRUE)
  if (!nrow(ij))
    return(NULL)
  w <- lambda[ij]
  ij[sample.int(nrow(ij), 1L, prob = w), , drop = TRUE]
}

.classe_td_event_rates <- function(array, state) {
  lambda <- array$lam.tensor[[state]]
  ij <- which(lambda > 0, arr.ind = TRUE)
  q <- array$Q[state, ]
  q[state] <- 0
  list(
    lambda.ij = ij,
    lambda = if (nrow(ij)) lambda[ij] else numeric(),
    mu = unname(array$mu[state]),
    q.to = which(q > 0),
    q = q[q > 0]
  )
}

#' Mark a reconstructed tree as an episodic ClaSSE simulation backbone
#'
#' Creates a stable representation of a supplied extant tree, its tip states,
#' and reconstructed ancestral-state probabilities. Internal-node states are
#' interpreted as states immediately before the observed cladogenetic event.
#'
#' @param tree A rooted, bifurcating, ultrametric `phylo` object.
#' @param node.states Named states for all internal nodes, or a `k` by `Nnode`
#'   matrix of reconstructed state probabilities. Matrix columns are retained as
#'   soft endpoint constraints rather than reduced to maximum-probability states.
#' @param tip.states Named states for all tips.
#' @param k Number of ClaSSE states.
#' @param state.labels Optional state labels in model order.
#' @param stem.length Length of the stem lineage. By default, uses
#'   `tree$root.edge` when present, otherwise zero.
#' @param stem.state State at the old end of the stem. By default, uses the
#'   reconstructed crown-root state.
#'
#' @return An object of class `classe_td_backbone`.
#'
#' @details A vector of reconstructed node states is represented internally as
#'   one-hot probabilities. If a probability matrix is supplied, each column is
#'   normalized and retained. During simulation, arrival in state `i` at that
#'   node is accepted with the reconstructed probability for state `i`; a
#'   rejected arrival triggers a local retry. Thus uncertain reconstructions do
#'   not impose a single maximum-probability endpoint on every branch.
#' @export
mark.classe.td.backbone <- function(tree, node.states, tip.states, k,
                                    state.labels = as.character(seq_len(k)),
                                    stem.length = NULL, stem.state = NULL) {
  if (!inherits(tree, "phylo"))
    stop("tree must inherit from 'phylo'")
  if (!ape::is.rooted(tree) || !ape::is.binary(tree))
    stop("tree must be rooted and bifurcating")
  if (!ape::is.ultrametric(tree))
    stop("the reconstructed backbone tree must be ultrametric")
  if (length(state.labels) != k || anyDuplicated(state.labels))
    stop("state.labels must contain k unique labels")

  ntip <- ape::Ntip(tree)
  internal <- ntip + seq_len(tree$Nnode)
  root <- setdiff(unique(tree$edge[, 1]), unique(tree$edge[, 2]))
  if (length(root) != 1L)
    stop("tree must have exactly one root")

  node.data <- .classe_td_state_data(
    node.states, internal, k, state.labels, "node.states"
  )
  tip.data <- .classe_td_state_data(
    tip.states, tree$tip.label, k, state.labels, "tip.states"
  )
  if (is.null(stem.length))
    stem.length <- if (is.null(tree$root.edge)) 0 else tree$root.edge
  if (length(stem.length) != 1L || !is.finite(stem.length) || stem.length < 0)
    stop("stem.length must be one finite non-negative number")
  if (is.null(stem.state))
    stem.state <- node.data$states[as.character(root)]
  stem.state <- .classe_td_state_vector(
    unname(stem.state), "stem", k, state.labels, "stem.state"
  )[[1]]

  depth <- ape::node.depth.edgelength(tree)
  height <- max(depth[seq_len(ntip)])
  age <- height - depth
  edge.table <- data.frame(
    edge.id = seq_len(nrow(tree$edge)),
    parent = tree$edge[, 1], child = tree$edge[, 2],
    parent.age = age[tree$edge[, 1]], child.age = age[tree$edge[, 2]],
    length = tree$edge.length,
    stringsAsFactors = FALSE
  )

  out <- list(
    tree = tree, root = root, height = height, age = age,
    node.states = node.data$states, tip.states = tip.data$states,
    node.probs = node.data$probabilities, tip.probs = tip.data$probabilities,
    probabilistic.nodes = node.data$soft,
    edge.table = edge.table, stem.length = stem.length, stem.state = stem.state,
    k = k, state.labels = state.labels
  )
  class(out) <- "classe_td_backbone"
  out
}

.simulate_classe_td_once <- function(backbone, schedule, max.events, max.taxa,
                                     max.branch.tries = 100L,
                                     reject.extra.extant = FALSE,
                                     guided = FALSE, initial.root.state = NULL,
                                     side.lineage.edges = NULL,
                                     side.lineage.mode = "all",
                                     max.side.tries = 10000L,
                                     extinction.cache = NULL,
                                     descendant.cache = NULL) {
  tree <- backbone$tree
  state.labels <- schedule$state.labels
  edges <- list()
  tips <- list()
  events <- list()
  pending.sides <- list()
  branch.attempts <- integer()
  side.attempts <- integer()
  condition.sides <- identical(side.lineage.mode, "extinct_clades")
  node.counter <- 0L
  event.counter <- 0L
  extinct.counter <- 0L
  extant.counter <- 0L
  n.events <- 0L
  n.terminals <- 0L
  failed <- NULL
  root.state.used <- NA_integer_

  new.node <- function(prefix = "g") {
    node.counter <<- node.counter + 1L
    paste0(prefix, node.counter)
  }
  bnode <- function(node) paste0("b", node)
  add.event <- function(type, age, epoch, node, from = NA_integer_,
                        to1 = NA_integer_, to2 = NA_integer_, backbone.edge = NA_integer_) {
    event.counter <<- event.counter + 1L
    events[[event.counter]] <<- data.frame(
      event.id = event.counter, event = type, age = age, epoch = epoch,
      node = node, from = from, to1 = to1, to2 = to2,
      backbone.edge = backbone.edge, stringsAsFactors = FALSE
    )
  }
  add.edge <- function(parent, child, start.age, end.age, map, backbone.edge = NA_integer_,
                       generated = TRUE) {
    if (!length(map) || abs(sum(map) - (start.age - end.age)) > 1e-7)
      stop("Internal error: state map does not sum to edge length")
    edges[[length(edges) + 1L]] <<- list(
      parent = parent, child = child, start.age = start.age, end.age = end.age,
      length = start.age - end.age, map = map,
      backbone.edge = backbone.edge, generated = generated
    )
  }
  add.tip <- function(node, label, age, state, fate, generated) {
    n.terminals <<- n.terminals + 1L
    if (n.terminals > max.taxa)
      stop(structure(list(message = "max.taxa exceeded"), class = c("classe_limit", "error", "condition")))
    tips[[node]] <<- data.frame(
      node = node, label = label, age = age, state = state,
      fate = fate, generated = generated, stringsAsFactors = FALSE
    )
  }
  count.event <- function() {
    n.events <<- n.events + 1L
    if (n.events > max.events)
      stop(structure(list(message = "max.events exceeded"), class = c("classe_limit", "error", "condition")))
  }
  snapshot <- function() {
    list(
      n.edges = length(edges), n.tips = length(tips), n.event.rows = length(events),
      n.pending = length(pending.sides),
      node.counter = node.counter, event.counter = event.counter,
      extinct.counter = extinct.counter, extant.counter = extant.counter,
      n.events = n.events, n.terminals = n.terminals
    )
  }
  restore <- function(x) {
    pending.sides <<- pending.sides[seq_len(x$n.pending)]
    if (length(edges) > x$n.edges)
      edges <<- edges[seq_len(x$n.edges)]
    if (x$n.edges == 0L)
      edges <<- list()
    if (length(tips) > x$n.tips)
      tips <<- tips[seq_len(x$n.tips)]
    if (x$n.tips == 0L)
      tips <<- list()
    if (length(events) > x$n.event.rows)
      events <<- events[seq_len(x$n.event.rows)]
    if (x$n.event.rows == 0L)
      events <<- list()
    node.counter <<- x$node.counter
    event.counter <<- x$event.counter
    extinct.counter <<- x$extinct.counter
    extant.counter <<- x$extant.counter
    n.events <<- x$n.events
    n.terminals <<- x$n.terminals
    failed <<- NULL
  }

  draw.event <- function(age, state, stop.age, map,
                         allow.extinction = TRUE,
                         retain.spine.state = FALSE) {
    repeat {
      epoch <- .classe_td_epoch(age, schedule$boundaries)
      boundary <- max(stop.age, .classe_td_lower_boundary(age, schedule$boundaries))
      rates <- .classe_td_event_rates(schedule$arrays[[epoch]], state)
      if (retain.spine.state && length(rates$lambda)) {
        compatible <- rates$lambda.ij[, 1] == state |
          rates$lambda.ij[, 2] == state
        rates$lambda.ij <- rates$lambda.ij[compatible, , drop = FALSE]
        rates$lambda <- rates$lambda[compatible]
      }
      mu.rate <- if (allow.extinction) rates$mu else numeric()
      all.rates <- c(rates$lambda, mu.rate, rates$q)
      total <- sum(all.rates)
      span <- age - boundary
      if (total <= 0) {
        map <- .append_map(map, state, span, state.labels)
        age <- boundary
      } else {
        dt <- stats::rexp(1L, total)
        if (dt >= span) {
          map <- .append_map(map, state, span, state.labels)
          age <- boundary
        } else {
          map <- .append_map(map, state, dt, state.labels)
          age <- age - dt
          z <- sample.int(length(all.rates), 1L, prob = all.rates)
          nl <- length(rates$lambda)
          if (z <= nl) {
            return(list(type = "speciation", age = age, epoch = epoch,
                        pair = rates$lambda.ij[z, ], state = state, map = map))
          }
          nmu <- length(mu.rate)
          if (nmu && z == nl + 1L) {
            return(list(type = "extinction", age = age, epoch = epoch,
                        state = state, map = map))
          }
          to <- rates$q.to[z - nl - nmu]
          return(list(type = "transition", age = age, epoch = epoch,
                      state = state, to = to, map = map))
        }
      }
      if (age <= stop.age + 1e-12)
        return(list(type = "end", age = stop.age,
                    epoch = .classe_td_epoch(stop.age, schedule$boundaries),
                    state = state, map = map))
    }
  }

  simulate.side <- function(parent, state, start.age) {
    edge.start <- start.age
    map <- numeric()
    age <- start.age
    repeat {
      ev <- if (condition.sides && !is.null(extinction.cache)) {
        .classe_td_draw_extinct_event(age, state, map, schedule, extinction.cache)
      } else {
        draw.event(age, state, 0, map, allow.extinction = TRUE)
      }
      age <- ev$age
      map <- ev$map
      # Stop at the first survivor; do not generate the rest of a rejected tree.
      if (condition.sides && ev$type == "end")
        stop(structure(list(message = "side clade reached extinction deadline"),
                       class = c("classe_reject", "error", "condition")))
      if (ev$type == "transition") {
        count.event()
        add.event("transition", age, ev$epoch, parent, state, ev$to)
        state <- ev$to
        next
      }
      if (ev$type == "speciation") {
        count.event()
        nd <- new.node("s")
        add.edge(parent, nd, edge.start, age, map)
        add.event("speciation", age, ev$epoch, nd, state, ev$pair[1], ev$pair[2])
        simulate.side(nd, ev$pair[1], age)
        simulate.side(nd, ev$pair[2], age)
        return(invisible(NULL))
      }
      terminal <- new.node("t")
      add.edge(parent, terminal, edge.start, age, map)
      if (ev$type == "extinction") {
        count.event()
        extinct.counter <<- extinct.counter + 1L
        label <- sprintf("extinct_%05d", extinct.counter)
        add.tip(terminal, label, age, state, "extinct", TRUE)
        add.event("extinction", age, ev$epoch, terminal, state)
      } else {
        if (reject.extra.extant)
          stop(structure(
            list(message = "generated lineage survived to the present"),
            class = c("classe_reject", "error", "condition")
          ))
        extant.counter <<- extant.counter + 1L
        label <- sprintf("extant_%05d", extant.counter)
        add.tip(terminal, label, 0, state, "extra_extant", TRUE)
        add.event("present", 0, 1L, terminal, state)
      }
      return(invisible(NULL))
    }
  }

  simulate.spine <- function(parent, state, start.age, stop.age, target.probs,
                             child, backbone.edge) {
    edge.start <- start.age
    map <- numeric()
    age <- start.age
    repeat {
      # This is an observed lineage in the supplied reconstructed tree. Its
      # survival is conditioned on the observation; extinction is simulated
      # only for the unobserved daughter created at a hidden speciation event.
      ev <- if (!is.null(descendant.cache)) {
        .classe_td_draw_backbone_event(age, state, stop.age, map,
          backbone.edge, schedule, extinction.cache, descendant.cache)
      } else draw.event(
        age, state, stop.age, map, allow.extinction = FALSE,
        retain.spine.state = !guided
      )
      age <- ev$age
      map <- ev$map
      if (ev$type == "transition") {
        count.event()
        add.event("transition", age, ev$epoch, parent, state, ev$to,
                  backbone.edge = backbone.edge)
        state <- ev$to
        next
      }
      if (ev$type == "extinction") {
        failed <<- "backbone_extinct"
        return(list(ok = FALSE, state = NA_integer_))
      }
      if (ev$type == "speciation") {
        count.event()
        nd <- new.node("h")
        add.edge(parent, nd, edge.start, age, map, backbone.edge, generated = FALSE)
        # The legacy sampler retains the parental spine state. The guided
        # sampler instead randomly assigns both daughters from the full tensor.
        spine.state <- state
        side.state <- if (ev$pair[1] == state) ev$pair[2] else ev$pair[1]
        if (guided) {
          oriented <- ev$pair
          if (is.null(descendant.cache) && oriented[1] != oriented[2] && stats::runif(1) < 0.5)
            oriented <- rev(oriented)
          spine.state <- oriented[1]
          side.state <- oriented[2]
        }
        add.event("hidden_speciation", age, ev$epoch, nd, state,
                  spine.state, side.state, backbone.edge)
        side.ok <- tryCatch({
          if (guided) {
            pending.sides[[length(pending.sides) + 1L]] <<-
              list(parent = nd, state = side.state, age = age,
                   backbone.edge = backbone.edge)
          } else {
            simulate.side(nd, side.state, age)
          }
          TRUE
        }, classe_reject = function(e) FALSE)
        if (!side.ok) {
          failed <<- "extra_extant"
          return(list(ok = FALSE, state = NA_integer_))
        }
        parent <- nd
        state <- spine.state
        edge.start <- age
        map <- numeric()
        next
      }
      if (is.null(descendant.cache) && stats::runif(1) > target.probs[state]) {
        failed <<- "endpoint_state_rejected"
        return(list(ok = FALSE, state = state))
      }
      add.edge(parent, child, edge.start, stop.age, map,
               backbone.edge, generated = FALSE)
      return(list(ok = TRUE, state = state))
    }
  }

  walk.backbone <- function(node, parent.state) {
    if (!is.null(failed))
      return(FALSE)
    rows <- which(tree$edge[, 1] == node)
    if (length(rows) != 2L) {
      failed <<- "backbone_not_binary"
      return(FALSE)
    }
    age <- backbone$age[node]
    epoch <- .classe_td_epoch(age, schedule$boundaries)
    children <- tree$edge[rows, 2]
    target.probs <- lapply(children, function(child.node) {
      if (child.node <= ape::Ntip(tree)) {
        backbone$tip.probs[, tree$tip.label[child.node]]
      } else {
        backbone$node.probs[, as.character(child.node)]
      }
    })
    child.scores <- if (!is.null(descendant.cache)) descendant.cache$child.scores(node) else if (guided) rep(list(rep(1, backbone$k)), 2L) else lapply(seq_along(children), function(z) {
      transition <- .classe_td_transition_matrix(
        age, backbone$age[children[z]], schedule
      )
      as.numeric(transition %*% target.probs[[z]])
    })

    node.ok <- FALSE
    last.failure <- NA_character_
    for (attempt in seq_len(max.branch.tries)) {
      key <- as.character(node)
      if (is.na(branch.attempts[key])) branch.attempts[key] <<- 0L
      branch.attempts[key] <<- branch.attempts[key] + 1L
      before <- snapshot()
      pair <- .classe_td_observed_pair(
        schedule$arrays[[epoch]], parent.state, child.scores
      )
      if (is.null(pair)) {
        failed <<- "incompatible_observed_node"
        return(FALSE)
      }
      add.event("observed_speciation", age, epoch, bnode(node), parent.state,
                pair[1], pair[2])
      node.ok <- TRUE
      child.states <- integer(length(children))
      for (z in seq_along(rows)) {
        edge.row <- rows[z]
        child.node <- children[z]
        branch <- simulate.spine(
          bnode(node), pair[z], age, backbone$age[child.node], target.probs[[z]],
          bnode(child.node), edge.row
        )
        if (!branch$ok) {
          node.ok <- FALSE
          break
        }
        child.states[z] <- branch$state
      }
      # After accepting this pair, process descendant observed nodes. In the
      # guided sampler their retry exhaustion fails the draw, preserving the
      # rule that accepted ancestral branches are never reweighted by retries.
      if (node.ok) {
        for (z in seq_along(children)) {
          child.node <- children[z]
          if (child.node <= ape::Ntip(tree)) {
            add.tip(
              bnode(child.node), tree$tip.label[child.node], 0,
              child.states[z], "observed", FALSE
            )
          } else if (!walk.backbone(child.node, child.states[z])) {
            if (guided) return(FALSE)
            node.ok <- FALSE
            break
          }
        }
      }
      if (node.ok)
        break
      last.failure <- failed
      restore(before)
    }
    if (!node.ok) {
      failed <<- paste0(
        "max_branch_tries_node_", node,
        if (!is.na(last.failure)) paste0("_", last.failure) else ""
      )
      return(FALSE)
    }

    TRUE
  }

  run.backbone <- function() {
    if (guided) {
      root.state <- initial.root.state
    } else if (backbone$stem.length > 0) {
      stem <- simulate.spine(
        "stem_origin", backbone$stem.state,
        backbone$height + backbone$stem.length, backbone$height,
        backbone$node.probs[, as.character(backbone$root)],
        bnode(backbone$root), 0L
      )
      if (!stem$ok)
        return(FALSE)
      root.state <- stem$state
    } else {
      root.probs <- backbone$node.probs[, as.character(backbone$root)]
      rows <- which(tree$edge[, 1] == backbone$root)
      children <- tree$edge[rows, 2]
      child.probs <- lapply(children, function(child.node) {
        if (child.node <= ape::Ntip(tree))
          backbone$tip.probs[, tree$tip.label[child.node]]
        else
          backbone$node.probs[, as.character(child.node)]
      })
      child.scores <- lapply(seq_along(children), function(z) {
        transition <- .classe_td_transition_matrix(
          backbone$age[backbone$root], backbone$age[children[z]], schedule
        )
        as.numeric(transition %*% child.probs[[z]])
      })
      epoch <- .classe_td_epoch(
        backbone$age[backbone$root], schedule$boundaries
      )
      compatibility <- vapply(seq_len(backbone$k), function(state) {
        .classe_td_observed_weight(
          schedule$arrays[[epoch]], state, child.scores
        )
      }, numeric(1))
      root.weights <- root.probs * compatibility
      if (sum(root.weights) <= 0) {
        failed <<- "incompatible_root_state"
        return(FALSE)
      }
      root.state <- sample.int(backbone$k, 1L, prob = root.weights)
    }
    root.state.used <<- root.state
    walk.backbone(backbone$root, root.state)
  }
  ok <- tryCatch(
    {
      accepted <- run.backbone()
      if (isTRUE(accepted) && guided) {
        for (side in pending.sides) {
          if (is.null(side.lineage.edges) || side$backbone.edge %in% side.lineage.edges) {
            before <- snapshot()
            side.ok <- FALSE
            for (attempt in seq_len(if (condition.sides && is.null(extinction.cache)) max.side.tries else 1L)) {
              side.ok <- tryCatch({
                simulate.side(side$parent, side$state, side$age)
                TRUE
              }, classe_reject = function(e) FALSE)
              if (side.ok) break
              restore(before)
            }
            side.attempts[side$parent] <- attempt
            if (!side.ok) {
              failed <- paste0("max_side_tries_", side$parent,
                               "_extinction_before_present",
                               " (founding state ", state.labels[side$state],
                               ", age ", signif(side$age, 7), ")")
              return(list(success = FALSE, reason = failed))
            }
          }
        }
      }
      accepted
    },
    classe_reject = function(e) {
      failed <<- "extra_extant"
      FALSE
    },
    classe_limit = function(e) {
      failed <<- if (grepl("max.taxa", e$message)) "max_taxa" else "max_events"
      FALSE
    }
  )
  if (!isTRUE(ok))
    return(list(success = FALSE, reason = failed, events = events,
                n.events = n.events, extra.extant = extant.counter,
                root.state = root.state.used))

  event.table <- if (length(events)) do.call(rbind, events) else data.frame()
  graph <- .classe_td_graph_to_simmap(
    edges, tips, state.labels,
    root.start = if (backbone$stem.length > 0) "stem_origin" else NULL
  )
  event.table$phylo.node <- unname(graph$node.map$phylo.node[
    match(event.table$node, graph$node.map$graph.node)
  ])
  event.table$from.label <- ifelse(
    is.na(event.table$from), NA_character_, state.labels[event.table$from]
  )
  event.table$to1.label <- ifelse(
    is.na(event.table$to1), NA_character_, state.labels[event.table$to1]
  )
  event.table$to2.label <- ifelse(
    is.na(event.table$to2), NA_character_, state.labels[event.table$to2]
  )
  graph$tips$state.label <- state.labels[graph$tips$state]
  graph$tree$tip.status <- graph$tips
  graph$tree$node.map <- graph$node.map
  graph$tree$cladogenetic.events <- event.table[
    event.table$event %in% c("observed_speciation", "hidden_speciation", "speciation"),
    , drop = FALSE
  ]
  graph$tree$event.history <- event.table
  graph$tree$lineage.history <- graph$lineages
  class(graph$tree) <- c("classe_td_simmap", "simmap", "phylo")
  list(success = TRUE, reason = if (extant.counter) "extra_extant" else NA_character_,
       tree = graph$tree, tips = graph$tips, lineages = graph$lineages,
       events = event.table, node.map = graph$node.map,
       extra.extant = extant.counter,
       n.events = n.events, root.state = root.state.used,
       branch.attempts = branch.attempts, side.attempts = side.attempts)
}

.classe_td_graph_to_simmap <- function(edges, tips, state.labels,
                                       root.start = NULL) {
  all.edges <- edges
  root.edge <- NULL
  root.map <- NULL
  if (!is.null(root.start)) {
    root.idx <- which(vapply(edges, `[[`, character(1), "parent") == root.start)
    if (length(root.idx) != 1L)
      stop("Internal error: mapped stem must have one initial edge")
    root.edge <- edges[[root.idx]]$length
    root.map <- edges[[root.idx]]$map
    edges <- edges[-root.idx]
  }
  terminals <- names(tips)
  original <- terminals[!vapply(tips, function(x) x$generated, logical(1))]
  generated <- setdiff(terminals, original)
  terminal.order <- c(original, generated)
  all.nodes <- unique(c(
    vapply(edges, `[[`, character(1), "parent"),
    vapply(edges, `[[`, character(1), "child")
  ))
  internal <- setdiff(all.nodes, terminal.order)
  ids <- stats::setNames(seq_along(c(terminal.order, internal)), c(terminal.order, internal))
  n.tip <- length(terminal.order)
  ids[internal] <- n.tip + seq_along(internal)

  edge <- do.call(rbind, lapply(edges, function(x) c(ids[x$parent], ids[x$child])))
  storage.mode(edge) <- "integer"
  edge.length <- vapply(edges, `[[`, numeric(1), "length")
  maps <- lapply(edges, `[[`, "map")
  tip.label <- vapply(tips[terminal.order], function(x) x$label, character(1))
  phy <- list(edge = edge, edge.length = edge.length,
              tip.label = tip.label, Nnode = length(internal), maps = maps)
  if (!is.null(root.edge)) {
    phy$root.edge <- root.edge
    phy$root.map <- root.map
  }
  class(phy) <- c("classe_td_simmap", "simmap", "phylo")

  old.keys <- paste(edge[, 1], edge[, 2], sep = ",")
  phy <- ape::reorder.phylo(phy, "cladewise")
  new.keys <- paste(phy$edge[, 1], phy$edge[, 2], sep = ",")
  ord <- match(new.keys, old.keys)
  phy$maps <- maps[ord]
  phy$edge.length <- edge.length[ord]
  mapped <- matrix(0, nrow(phy$edge), length(state.labels),
                   dimnames = list(new.keys, state.labels))
  for (i in seq_along(phy$maps)) {
    totals <- tapply(as.numeric(phy$maps[[i]]), names(phy$maps[[i]]), sum)
    mapped[i, names(totals)] <- totals
  }
  phy$mapped.edge <- mapped
  attr(phy, "map.order") <- "right-to-left"

  tip.table <- do.call(rbind, tips[terminal.order])
  rownames(tip.table) <- NULL
  lineage.table <- do.call(rbind, lapply(seq_along(all.edges), function(i) {
    x <- all.edges[[i]]
    data.frame(
      lineage.id = i, parent = x$parent, child = x$child,
      start.age = x$start.age, end.age = x$end.age, length = x$length,
      backbone.edge = x$backbone.edge, generated = x$generated,
      root.edge = !is.null(root.start) && x$parent == root.start,
      map = paste(paste(names(x$map), signif(x$map, 10), sep = ":"), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  node.map <- data.frame(
    graph.node = names(ids), phylo.node = as.integer(ids),
    stringsAsFactors = FALSE
  )
  list(tree = phy, tips = tip.table, lineages = lineage.table,
       node.map = node.map)
}

.classe_td_nsim <- function(nsim) {
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1 || nsim != as.integer(nsim))
    stop("nsim must be one positive integer")
  as.integer(nsim)
}

.classe_td_simulation_set <- function(x) {
  if (length(x) == 1L)
    return(x[[1L]])
  structure(x, class = c("classe_td_tree_simulations", "list"))
}

#' Simulate an episodic ClaSSE candidate on a reconstructed tree
#'
#' Simulates geographic transitions and hidden speciation events along a fixed
#' reconstructed backbone. Side trees are simulated without conditioning, so
#' the returned candidate may contain additional extant tips.
#'
#' @param backbone A `classe_td_backbone` object, or a tree when `node.states`
#'   and `tip.states` are supplied.
#' @param pars Named full parameter vector from `make.classe.td()`.
#' @param node.states,tip.states,k,state.labels Passed to
#'   [mark.classe.td.backbone()] when `backbone` is a tree.
#' @param max.tries Maximum attempts to obtain a candidate that preserves all
#'   supplied backbone lineages and endpoint states.
#' @param max.events,max.taxa Safety limits for each attempt.
#' @param max.branch.tries Maximum local rejection attempts for a pair of
#'   observed daughter branches to satisfy their endpoint constraints.
#' @param nsim Number of independent stochastic maps to generate. The default
#'   returns one `classe_td_tree_sim`; values above one return a list of such
#'   objects with class `classe_td_tree_simulations`.
#' @param seed Optional random seed.
#'
#' @return A `classe_td_tree_sim` object whose `tree` element is a stochastic
#'   mapping tree of class `simmap` and `phylo`, or for `nsim > 1`, a
#'   `classe_td_tree_simulations` list containing those objects.
#'
#' @details The supplied backbone consists of observed lineages, so their
#'   survival is conditioned upon. Along each backbone branch, anagenetic and
#'   hidden speciation events are simulated while unconditional extinction is
#'   omitted for the continuing spine. At every hidden speciation, the second,
#'   unobserved daughter is simulated with the full process and may either go
#'   extinct or survive to the present.
#'
#' @examples
#' tree <- ape::read.tree(text = "(a:0.5,b:0.5);")
#' pars <- c(
#'   t.1 = 0.25,
#'   lambda111.1 = 0.8, mu1.1 = 2,
#'   lambda111.2 = 0.8, mu1.2 = 2
#' )
#' sim <- simulate.classe.td.on.tree(
#'   tree, pars, node.states = c(`3` = 1),
#'   tip.states = c(a = 1, b = 1), k = 1,
#'   state.labels = "A", seed = 3
#' )
#' sim$tree$maps
#' @export simulate.classe.td.on.tree
simulate.classe.td.on.tree <- function(
    backbone, pars, node.states = NULL, tip.states = NULL, k = NULL,
    state.labels = NULL, max.tries = 1000L, max.events = 100000L,
    max.taxa = 100000L, max.branch.tries = 100L, nsim = 1L,
    seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  nsim <- .classe_td_nsim(nsim)
  if (!inherits(backbone, "classe_td_backbone")) {
    if (is.null(k))
      stop("k is required when backbone is a phylo object")
    if (is.null(state.labels))
      state.labels <- as.character(seq_len(k))
    backbone <- mark.classe.td.backbone(
      backbone, node.states, tip.states, k, state.labels
    )
  }
  if (is.null(state.labels))
    state.labels <- backbone$state.labels
  schedule <- .classe_td_schedule(pars, backbone$k, state.labels)

  simulations <- vector("list", nsim)
  for (simulation.index in seq_len(nsim)) {
    failures <- integer()
    accepted <- FALSE
    for (attempt in seq_len(max.tries)) {
      ans <- .simulate_classe_td_once(
        backbone, schedule, max.events, max.taxa, max.branch.tries
      )
      if (isTRUE(ans$success)) {
        ans$backbone <- backbone
        ans$parameters <- pars
        ans$attempts <- attempt
        ans$failures <- failures
        ans$conditioned <- FALSE
        ans$endpoint.constraint <- if (backbone$probabilistic.nodes) "probability" else "fixed"
        ans$simulation <- simulation.index
        class(ans) <- "classe_td_tree_sim"
        simulations[[simulation.index]] <- ans
        accepted <- TRUE
        break
      }
      reason <- ans$reason
      if (is.null(reason) || is.na(reason))
        reason <- "unknown"
      failures[reason] <- if (is.na(failures[reason])) 1L else failures[reason] + 1L
    }
    if (!accepted)
      stop("No backbone-compatible simulation ", simulation.index, " of ", nsim,
           " in ", max.tries, " attempts. Failure counts: ",
           paste(names(failures), failures, sep = "=", collapse = ", "))
  }
  .classe_td_simulation_set(simulations)
}

#' Reconstruct the extant tree from a simulated ClaSSE history
#'
#' Prunes genuinely extinct generated tips, then checks whether any additional
#' extant lineage remains. Extra extant tips are never silently removed.
#'
#' @param simulation A `classe_td_tree_sim` object.
#'
#' @return A list containing the reconstructed `phylo`, a backbone-match flag,
#'   and the labels of extra extant tips.
#' @export
reconstruct.classe.td <- function(simulation) {
  if (inherits(simulation, "classe_td_tree_simulations"))
    return(lapply(simulation, reconstruct.classe.td))
  if (!inherits(simulation, "classe_td_tree_sim"))
    stop("simulation must be a classe_td_tree_sim object")
  extinct <- simulation$tips$label[simulation$tips$fate == "extinct"]
  extra <- simulation$tips$label[simulation$tips$fate == "extra_extant"]
  phy <- simulation$tree
  class(phy) <- "phylo"
  phy$maps <- phy$mapped.edge <- NULL
  if (length(extinct))
    phy <- ape::drop.tip(phy, extinct)
  target <- simulation$backbone$tree
  # ape::drop.tip() discards the surviving continuation of a pruned split on
  # the stem instead of adding it to root.edge. The fixed backbone supplies
  # the exact collapsed stem length, so restore it before comparison.
  phy$root.edge <- target$root.edge
  root.same <- (is.null(phy$root.edge) && is.null(target$root.edge)) ||
    (!is.null(phy$root.edge) && !is.null(target$root.edge) &&
       isTRUE(all.equal(phy$root.edge, target$root.edge, tolerance = 1e-8)))
  same <- !length(extra) && root.same && isTRUE(all.equal(
    ape::reorder.phylo(phy, "cladewise"),
    ape::reorder.phylo(target, "cladewise"),
    use.edge.length = TRUE, use.tip.label = TRUE, use.node.label = FALSE
  ))
  list(tree = phy, matches.backbone = same, extra.extant = extra,
       rejected = !same)
}

#' Simulate extinct ClaSSE side lineages conditional on an observed tree
#'
#' Uses local rejection at observed nodes until all generated side lineages are
#' extinct and pruning those lineages recovers the supplied backbone.
#'
#' @inheritParams simulate.classe.td.on.tree
#' @param max.tries Maximum number of whole-history rejection attempts.
#' @param min.extinct Minimum number of generated extinct terminal species in
#'   an accepted history. Values above zero add an explicit ascertainment
#'   condition and are mainly useful for examples and visualization.
#'
#' @return A conditioned `classe_td_tree_sim` object with a `simmap` tree.
#'   For `nsim > 1`, returns a `classe_td_tree_simulations` list.
#'
#' @details Generated side lineages are rejected as soon as one reaches the
#'   present. Rejection is local to the observed node proposal, so already
#'   accepted parts of the backbone are not needlessly resimulated. Because
#'   backbone lineages are observed, their survival is conditioned upon:
#'   extinction events are simulated for generated daughters, not as
#'   unconditional termination events on the observed spine. Setting
#'   `min.extinct > 0` adds an ascertainment condition and can still require
#'   many attempts when the fitted extinction rates are low.
#'
#' @examples
#' tree <- ape::read.tree(text = "(a:0.5,b:0.5);")
#' pars <- c(
#'   t.1 = 0.25,
#'   lambda111.1 = 0.8, mu1.1 = 2,
#'   lambda111.2 = 0.8, mu1.2 = 2
#' )
#' sim <- simulate.extinct.classe.td(
#'   tree, pars, node.states = c(`3` = 1),
#'   tip.states = c(a = 1, b = 1), k = 1,
#'   state.labels = "A", seed = 3
#' )
#' if (requireNamespace("phytools", quietly = TRUE)) {
#'   phytools::plotSimmap(sim$tree, ftype = "off")
#' }
#' @export simulate.extinct.classe.td
simulate.extinct.classe.td <- function(
    backbone, pars, node.states = NULL, tip.states = NULL, k = NULL,
    state.labels = NULL, max.tries = 10000L, max.events = 100000L,
    max.taxa = 100000L, max.branch.tries = 100L, min.extinct = 0L,
    nsim = 1L, seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  nsim <- .classe_td_nsim(nsim)
  if (!inherits(backbone, "classe_td_backbone")) {
    if (is.null(k))
      stop("k is required when backbone is a phylo object")
    if (is.null(state.labels))
      state.labels <- as.character(seq_len(k))
    backbone <- mark.classe.td.backbone(
      backbone, node.states, tip.states, k, state.labels
    )
  }
  if (is.null(state.labels))
    state.labels <- backbone$state.labels
  schedule <- .classe_td_schedule(pars, backbone$k, state.labels)
  if (length(min.extinct) != 1L || is.na(min.extinct) || min.extinct < 0 ||
      min.extinct != as.integer(min.extinct))
    stop("min.extinct must be one non-negative integer")
  min.extinct <- as.integer(min.extinct)
  simulations <- vector("list", nsim)
  for (simulation.index in seq_len(nsim)) {
    failures <- integer()
    accepted <- FALSE
    for (attempt in seq_len(max.tries)) {
      ans <- .simulate_classe_td_once(
        backbone, schedule, max.events, max.taxa, max.branch.tries,
        reject.extra.extant = TRUE
      )
      reason <- ans$reason
      n.extinct <- if (isTRUE(ans$success)) sum(ans$tips$fate == "extinct") else 0L
      if (isTRUE(ans$success) && n.extinct >= min.extinct) {
        ans$backbone <- backbone
        ans$parameters <- pars
        ans$attempts <- attempt
        ans$failures <- failures
        ans$conditioned <- TRUE
        ans$endpoint.constraint <- if (backbone$probabilistic.nodes) "probability" else "fixed"
        ans$simulation <- simulation.index
        class(ans) <- "classe_td_tree_sim"
        check <- reconstruct.classe.td(ans)
        if (isTRUE(check$matches.backbone)) {
          simulations[[simulation.index]] <- ans
          accepted <- TRUE
          break
        }
        reason <- "reconstruction_mismatch"
      } else if (isTRUE(ans$success) && n.extinct < min.extinct) {
        reason <- "too_few_extinct"
      }
      if (is.null(reason) || is.na(reason))
        reason <- "unknown"
      failures[reason] <- if (is.na(failures[reason])) 1L else failures[reason] + 1L
    }
    if (!accepted)
      stop("No all-extra-extinct simulation ", simulation.index, " of ", nsim,
           " in ", max.tries, " attempts. Failure counts: ",
           paste(names(failures), failures, sep = "=", collapse = ", "))
  }
  .classe_td_simulation_set(simulations)
}
