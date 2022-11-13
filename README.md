# Coverage_Analysis

usage: coverageAnalysisV4.py [-h] -r <input.bam> -t <input.bam> [-reg <roi.txt>] -o <out.txt> [-rreads <input.bam>] [--index]

Check coverages in given region

optional arguments: \
  -h, --help           show this help message and exit 
  
  -r <input.bam>       reference alignment file 
  
  -t <input.bam>       target alignment file
  
  -reg <roi.txt>       Copy number variation region
  
  -o <out.txt>         outfile name
  
  -rreads <input.bam>  reference reads
  
  --index              will index bam files if index files have not been created
