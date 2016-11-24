PIGS Paper Code
==============

Code to run and analyze simulation study for the PIGS papers. 

### Organization

1. Run contains the code to generate the data, analyze it, and run that on Biowulf. 
2. Munge contains the code to analyze the resulting data. 

### Running on Biowulf

1. Run `generator.R` for the deisred scenario (1:6)
2. Run the following code, taken from `biowulf-run.csh`, replacing sachsmc with your username 

        Rswarm --rfile=Rfile.R --sfile=seedfile.txt --path=//data//sachsmc//pigs-sims//scenario1/ --reps=300 --sims=1 --start=0 --ext1=.RData
        swarm -f Rfile.sw --module R --sbatch "--mail-type=END"
   
3. Once step 2 has run, add additional bootstrap runs with 

        Rswarm --rfile=Rboot.R --sfile=seedfile2.txt --path=//data//sachsmc//pigs-sims//scenario%d/ --reps=300 --sims=1 --start=0 --ext1=.RData
        swarm -f Rboot.sw --module R --sbatch "--mail-type=END"

## License

The MIT License (MIT)
Copyright (c) 2016 Michael C Sachs

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
