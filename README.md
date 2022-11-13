# Coverage_Analysis

## Description

The correct identification of gene copy number variations is of great interest for determining the biological significance of the variations in assemblies of the same species. A naive way to identify gene copy number variation is to compare the gene annotations of two different gene assemblies. However any gene copy number variation between genome assemblies is not necessarily a true gene copy number variation; Genome assemblies often have assembly errors, so a region of suspected CNV might be an incorrectly assembled region, instead of true CNV. Assemblers are known to have difficulties with separating large repeats in areas such as the telomeres and centromeres, especially when using short read data (Kelley and Salzberg. 2010). Additionally, gene annotations are often incomplete, so a gene copy may be missing from the annotation despite being in the assembly. One possible way to provide additional evidence to support true gene copy number between two assemblies is to map the DNA-sequencing reads from one genome to each assembly and compare the read coverages in the regions of interest. Therefore, the purpose of this project is to create a Python tool to evaluate read coverage and mapping quality over a user-specified region to distinguish true Copy Number Variation from misassembled regions in an assembly or annotation errors. 

The user inputs to the tool are two Binary Alignment Map (BAM) files, index files of the BAM files, and a text file for the regions of interest. One BAM file contains alignments of the sequencing reads of one genome aligned to itself. The other BAM file contains the alignments of the reads to the other genome. The text file has the genomic coordinates of the original region in the comparison genome followed by the location of the region of suspected Copy Number Variation. 

To calculate coverage across the user-specified regions, The analysis tool uses pybedtools   to run bedtools genomecov with the -ibam and -d features for the read alignments to Genome 1 and Genome 2 (Dale et al., 2011). This creates a data structure of the read coverage at every nucleotide location, By iterating through this structure, the tool calculates the whole genome coverage as well as the coverages of the user specified regions for Genome 1 and Genome 2 . In the output file, the coverage of each region is normalized by dividing the region’s coverage by the total average coverage of its respective genome. 

To calculate mapping quality across the user-specified regions, pySam (Li et al., 2009) was implemented to process each read in the BAM files. We extracted the mapping quality for each read and used these values to calculate the average mapping quality of the genome. For every region, we found the reads fully or partially contained in the region, from which the average of the region was calculated. The mapping qualities across the regions were again normalized by dividing them by the average mapping quality of their respective genomes.
The output is a text file that includes the locations of the input files as well as the  normalized coverage and normalized mapping quality for the regions in the reference and target file. An example of the output format is shown in Table 1


## Usage 
coverageAnalysisV4.py [-h] -r <input.bam> -t <input.bam> [-reg <roi.txt>] -o <out.txt> [-rreads <input.bam>] [--index]

Check coverages in given region

optional arguments: \
  -h, --help           show this help message and exit 
  
  -r <input.bam>       reference alignment file 
  
  -t <input.bam>       target alignment file
  
  -reg <roi.txt>       Copy number variation region
  
  -o <out.txt>         outfile name
  
  -rreads <input.bam>  reference reads
  
  --index              will index bam files if index files have not been created
