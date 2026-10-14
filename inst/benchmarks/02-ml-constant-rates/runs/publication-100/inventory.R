# Checksums for the actual local inputs, scripts and retrieved fit records ----
files <- c(list.files("datasets",full.names=TRUE),
  list.files("results",pattern="^fit-[0-9]+[.]rds$",recursive=TRUE,full.names=TRUE),
  list.files(pattern="[.](R|yml|sh)$"), "manifest.csv", "design.csv",
  "provenance/divextra-source.tar.gz", "provenance/runtime-functions.rds")
files <- sort(unique(files[file.exists(files)]))
checksum <- tools::md5sum(files)
write.csv(data.frame(path=files,md5=unname(checksum)),"provenance/files-md5.csv",row.names=FALSE)
cat("Checksummed",length(files),"local files; missing remote fits are not included.\n")
