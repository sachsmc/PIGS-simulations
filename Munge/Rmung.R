source("simulate-data.R")
source("assemble-results.R")
scenario <- 2

files <- c(gsub("Rmung", "Rfile-bytrial", "DUMY1"),
           gsub("Rmung", "Rfile-overall", "DUMY1"),
           gsub("Rmung", "Rfile-leftout", "DUMY1"))

for(ffs in files) load(ffs)

thisS <- unlist(map(bytrial$data, function(d) d$S.obs[!is.na(d$S.obs)]))
thisS <- sort(unique(thisS))

oot <- c(run_tests(bytrial), D.pee = calc_DjD0_pee(bytrial, leftout) < .05,
  CEP.coverage(scenario, bytrial, leftout, overall$pseval, Sobs = thisS))

write.table(t(oot), file = gsub(".RData", ".dat", "DUMY1", fixed = TRUE), 
            row.names = FALSE, col.names = FALSE)
