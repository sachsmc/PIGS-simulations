write.table(round(10000 * runif(1000), 0), file = "seedfile.txt", row.names = FALSE, col.names = FALSE)

Rswarm --rfile=Rfile.R --sfile=seedfile.txt --path=//data//sachsmc//pigs-sims//scenario%d/ --reps=300 --sims=1 --start=0 --ext1=.RData

swarm -f Rfile.sw --module R --sbatch "--mail-type=END"

Rswarm --rfile=Rboot.R --sfile=seedfile2.txt --path=//data//sachsmc//pigs-sims//scenario%d/ --reps=300 --sims=1 --start=0 --ext1=.RData

swarm -f Rboot.sw --module R --sbatch "--mail-type=END"
