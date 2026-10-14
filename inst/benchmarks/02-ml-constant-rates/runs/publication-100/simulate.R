# Run locally from benchmarks/publication/constant-rates-100 ----
library(divextra)
dir.create("datasets", showWarnings = FALSE)
dir.create("preparation", showWarnings = FALSE)
stopifnot(!file.exists("manifest.csv"))

# Local-only fix for a living unsplit lineage at the bookkeeping boundary.
simulate <- divextra::make.tree.classe.td
parts <- as.list(body(simulate))
found <- 0L
for (i in seq_along(parts)) {
  if (identical(parts[[i]], quote(if (is.null(tb.e1)) return(NULL)))) {
    parts[[i]] <- quote(if (n.taxa == 0) return(NULL))
    found <- found + 1L
  }
}
stopifnot(found == 1L)
body(simulate) <- as.call(parts)
dir.create("provenance", showWarnings = FALSE)
dump("simulate", "provenance/simulator-used.R")

# Fixed proposal pool, rather than stopping a tree at its target size ----
design <- expand.grid(scenario = c("Slow", "Medium", "Fast"),
  medium.duration = c(30, 40, 50, 60), stringsAsFactors = FALSE)
rates <- c(Slow = .01, Medium = .1, Fast = .3)
design$s <- unname(rates[design$scenario])
design$d <- design$s
design$x <- design$s / 2
design$duration <- design$medium.duration * .1 / design$s
write.csv(design, "design.csv", row.names = FALSE)
records <- data.frame()
pool <- list()
started <- proc.time()[3]
for (cell in seq_len(nrow(design))) {
  row <- design[cell, ]
  pars <- pars.geosse2classe(c(row$s, row$s, row$s, row$x, row$x, row$d, row$d))
  for (replicate in 1:800) {
    seed <- 210000000 + cell * 10000 + replicate
    set.seed(seed)
    root <- sample(1:3, 1)
    tb <- simulate(rbind(pars, pars), k = 3, max.taxa = 2000,
      max.t1 = row$duration / 2, max.t2 = row$duration,
      x0 = root, single.lineage = TRUE)
    status <- "extinct"
    states <- integer()
    if (!is.null(tb)) {
      states <- tb$state[!tb$split & !tb$extinct]
      if (length(states) > 0) {
        status <- if (abs(attr(tb, "t.total") - row$duration) < 1e-8) "complete" else "capped"
      }
    }
    N <- length(states)
    counts <- tabulate(states, nbins = 3)
    eligible <- status == "complete" && N >= 50 && N <= 2000
    id <- paste0("c", cell, "-r", replicate)
    records <- rbind(records, data.frame(id = id, cell = cell, scenario = row$scenario,
      duration = row$duration, replicate = replicate, seed = seed, root = root,
      status = status, N = N, A = counts[1], B = counts[2], AB = counts[3], eligible = eligible))
    if (eligible) pool[[id]] <- tb
  }
  write.csv(records, "preparation/proposals.csv", row.names = FALSE)
  message("Simulation cell ", cell, "/12; elapsed ", round(proc.time()[3] - started), "s")
}
records$tip.group <- cut(records$N, c(50, 100, 250, 500, 1000, 2001), right = FALSE,
  labels = c("50-99", "100-249", "250-499", "500-999", "1000-2000"))
availability <- table(records$scenario[records$eligible], records$tip.group[records$eligible])
write.csv(as.data.frame(availability), "preparation/availability.csv", row.names = FALSE)
stopifnot(all(availability >= 100))

# Select 100 distinct trees per rate scenario and size group ----
set.seed(2026101401)
selected <- integer()
for (scenario in c("Slow", "Medium", "Fast")) {
  for (group in levels(records$tip.group)) {
    candidates <- which(records$scenario == scenario & records$eligible & records$tip.group == group)
    selected <- c(selected, candidates[sample.int(length(candidates), 100)])
  }
}
manifest <- records[selected, ]
manifest$tree.index <- seq_len(nrow(manifest))
manifest$task <- (manifest$tree.index - 1L) %/% 10L + 1L
stopifnot(nrow(manifest) == 1500, !anyDuplicated(manifest$seed))
diagnostics <- data.frame()
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  tb <- pool[[row$id]]
  tree <- diversitree:::me.to.ape.bisse(tb[-1, ], tb$state[1])
  tree <- diversitree::prune(tree)
  stopifnot(ape::Ntip(tree) == row$N, length(tree$tip.state) == row$N)
  stopifnot(all(tabulate(tree$tip.state, nbins = 3) == unlist(row[c("A", "B", "AB")])) )
  rate <- rates[row$scenario]
  z <- list(id = row$id, scenario = row$scenario, tip.group = as.character(row$tip.group),
    tree = tree, table = tb, seed = row$seed, root = row$root, duration = row$duration,
    recent.time = row$duration / 2,
    links = c(s = TRUE, d = TRUE, x = TRUE),
    truth = c(s_recent = unname(rate), d_recent = unname(rate), x_recent = unname(rate / 2),
      s_old = unname(rate), d_old = unname(rate), x_old = unname(rate / 2)))
  saveRDS(z, sprintf("datasets/tree-%03d.rds", i))
  depth <- ape::node.depth.edgelength(tree)
  for (fraction in .5) {
    cut <- max(depth) - row$duration * fraction
    crossing <- if (cut < 0) 1L else sum(depth[tree$edge[, 1]] < cut & depth[tree$edge[, 2]] >= cut)
    diagnostics <- rbind(diagnostics, data.frame(tree.index = i, fraction = fraction,
      boundary.age = row$duration * fraction, crown.age = max(depth), crossing = crossing,
      boundary.before.crown = cut < 0))
  }
}
write.csv(manifest, "manifest.csv", row.names = FALSE)
write.csv(diagnostics, "tree-diagnostics.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "preparation/session.txt")
print(table(manifest$scenario, manifest$tip.group))
