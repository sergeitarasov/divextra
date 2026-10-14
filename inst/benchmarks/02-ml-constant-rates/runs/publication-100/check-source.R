# Verify that archived critical functions match the installed runtime ----
library(divextra)
folder <- tempfile("divextra-source-")
dir.create(folder)
untar("provenance/divextra-source.tar.gz",exdir=folder)
definitions <- readRDS("provenance/runtime-functions.rds")
found <- character()
files <- c("R/classe_td.R","R/sim-classe-td.R","R/utils-yaml.R")
for(file in files) {
  expressions <- parse(file.path(folder,"divextra",file),keep.source=FALSE)
  for(expression in expressions) {
    if(!is.call(expression) || length(expression)!=3) next
    if(!as.character(expression[[1]]) %in% c("<-","=")) next
    if(!is.symbol(expression[[2]])) next
    name <- as.character(expression[[2]])
    if(!name %in% names(definitions)) next
    f <- eval(expression[[3]],envir=asNamespace("divextra"))
    stopifnot(identical(definitions[[name]],c(deparse(formals(f)),deparse(body(f)))))
    found <- c(found,name)
  }
}
stopifnot(setequal(found,names(definitions)))
writeLines(c("Source commit: 4aa1513ca3c1f9aa493f0ffae779cf67505f2e8e",
  "Clean source snapshot; all five critical source definitions match local installed runtime.",
  "Simulation additionally applies the documented local unsplit-stem safeguard."),
  "provenance/source-verification.txt")
