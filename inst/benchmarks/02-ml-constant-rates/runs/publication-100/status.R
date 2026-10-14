# Update completeness without presenting partial recovery as final ----
m <- read.csv("manifest.csv")
present <- valid <- converged <- timeouts <- 0L
for(i in m$tree.index) for(start in 1:3) {
  j <- ((i-1) %% 10)*3+start
  path <- sprintf("results/task-%03d/fit-%03d.rds",m$task[i],j)
  if(!file.exists(path)) next
  present <- present+1L
  f <- readRDS(path)
  ok <- is.null(f$error) && length(f$lnLik)==1 && is.finite(f$lnLik) && f$lnLik > -1e49 &&
    length(f$par)==3 && all(is.finite(f$par))
  valid <- valid+ok
  converged <- converged+isTRUE(ok && f$convergence==0)
  timeouts <- timeouts+isTRUE(f$timeout)
}
terminal <- FALSE
if(file.exists("accounting-inference.txt")) {
  a <- read.table("accounting-inference.txt",sep="|",fill=TRUE,stringsAsFactors=FALSE)
  terminal <- nrow(a)>0 && all(grepl("^(COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE)",a$V2))
}
finished <- present==4500 && terminal
write.csv(data.frame(expected=4500,present=present,valid=valid,converged=converged,
  timeouts=timeouts,remote.terminal=terminal,ready.for.analysis=finished),
  "retrieval-status.csv",row.names=FALSE)
writeLines(c("# Publication run status","",
  paste("Checked:",format(Sys.time(),tz="Europe/Helsinki",usetz=TRUE)),
  sprintf("Retrieved %d/4500 expected fits; valid: %d; convergence code 0: %d; timeouts: %d.",
    present,valid,converged,timeouts),
  if(finished) "All inference tasks are terminal and all fit records are available. Final local analysis follows."
    else "INCOMPLETE: do not interpret these partial results as final recovery. Remote jobs may still be running; see accounting-inference.txt.",
  "", "Remote jobs: 2300270 (smoke), 2300275 (main), 2300286 (summary).",
  "Original project and remote run paths were preserved during migration."),"run-status.md")
if(finished) source("analyse.R")
cat(present,"/4500 fits;",if(finished) "analysis complete" else "final analysis deferred","\n")
