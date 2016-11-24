generator <- function(scenario){

  dirout <- paste0("scenario", scenario)
  stopifnot(!dir.exists(dirout))

  dir.create(dirout)

  boot <- readLines("Rboot.R")
  writeLines(sprintf(boot, scenario), con = paste0(dirout, "/Rboot.R"))

  file <- readLines("Rfile.R")
  writeLines(sprintf(file, scenario), con = paste0(dirout, "/Rfile.R"))

  csh <- readLines("biowulf-run.csh")
  writeLines(sprintf(csh, scenario), con = paste0(dirout, "/biowulf-run.csh"))

  file.copy(c("seedfile.txt", "seedfile2.txt"), paste0(dirout, "/", c("seedfile.txt", "seedfile2.txt")))


}
