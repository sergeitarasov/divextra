# Ten trees per array task; three starts at the known generating boundary ----
source("model.R")
manifest <- read.csv("manifest.csv")
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
stopifnot(task >= 1, task <= max(manifest$task))
indices <- manifest$tree.index[manifest$task == task]
out <- sprintf("results/task-%03d", task)
dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))
jobs <- list()
for (i in indices) {
  z <- readRDS(sprintf("datasets/tree-%03d.rds", i))
  for (start in 1:3) {
    set.seed(z$seed + start * 10000000)
    # Identical absolute starting strategy across scenarios; no truth input.
    initial <- c(s_recent=.08, d_recent=.08, x_recent=.03) *
      c(.5, 1, 2)[start] * exp(runif(3, -log(2), log(2)))
    jobs[[length(jobs) + 1L]] <- list(tree.index=i, variant="linked50",
      model="linked_classe", fraction=.5, start=start, initial=initial)
  }
}
saveRDS(jobs, file.path(out, "jobs.rds"))
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
pending <- integer()
for (j in seq_along(jobs)) {
  if (!file.exists(file.path(out, sprintf("fit-%03d.rds", j)))) pending <- c(pending, j)
}
active <- list()
while (length(pending) > 0 || length(active) > 0) {
  while (length(pending) > 0 && length(active) < workers) {
    j <- pending[1]
    pending <- pending[-1]
    job <- jobs[[j]]
    z <- readRDS(sprintf("datasets/tree-%03d.rds", job$tree.index))
    z$recent.time <- z$duration * job$fraction
    p <- callr::r_bg(function(z, job) {
      source("model.R")
      ans <- fit.one(z, job)
      ans$tree.index <- job$tree.index
      ans$variant <- job$variant
      ans$fraction <- job$fraction
      ans
    }, args = list(z = z, job = job), libpath = .libPaths(), supervise = TRUE,
      stdout = file.path(out, "logs", sprintf("fit-%03d.out", j)),
      stderr = file.path(out, "logs", sprintf("fit-%03d.err", j)))
    active[[as.character(j)]] <- list(process = p, begin = Sys.time())
  }
  for (key in names(active)) {
    item <- active[[key]]
    expired <- as.numeric(difftime(Sys.time(), item$begin, units = "secs")) > 1200
    if (expired && item$process$is_alive()) item$process$kill()
    if (!item$process$is_alive()) {
      result <- tryCatch(item$process$get_result(), error = function(e)
        list(error = conditionMessage(e), timeout = expired))
      if (is.null(result)) result <- list(error = "Worker returned NULL", timeout = expired)
      path <- file.path(out, sprintf("fit-%03d.rds", as.integer(key)))
      saveRDS(result, paste0(path, ".tmp"))
      stopifnot(file.rename(paste0(path, ".tmp"), path))
      active[[key]] <- NULL
      message("Saved fit ", key, "/", length(jobs))
    }
  }
  Sys.sleep(.2)
}

# Smoke tasks cover all three rate scales; fail the dependency gate on invalid fits.
if (task %in% c(1L, 51L, 101L)) {
  for (j in seq_along(jobs)) {
    fit <- readRDS(file.path(out, sprintf("fit-%03d.rds", j)))
    stopifnot(is.null(fit$error), length(fit$par) == 3,
      is.finite(fit$lnLik), fit$lnLik > -1e49)
  }
}
