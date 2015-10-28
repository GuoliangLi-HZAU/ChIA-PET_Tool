# ChIA-PET Tool

## Introduction
**Ch**romatin **I**nteraction **A**nalysis with **P**aired-**E**nd **T**ag (**ChIA-PET**) sequencing is a genome-wide high-throughput technology that detects chromatin interactions associated with a specific protein of interest (Fullwood *et al.,* 2009). **ChIA-PET Tool** (Li *et al.,* 2010) is a computational package to process the next-generation sequence data generated from ChIA-PET wet-lab experiments, 
<br /><br />
which contains 7 steps: 
> 1) &nbsp; linker filtering<br />
> 2) &nbsp; mapping the paired-end reads to a reference genome<br />
> 3) &nbsp; purifying the mapped reads<br />
> 4) &nbsp; dividing the reads into different categories<br />
> 5) &nbsp; peak calling<br />
> 6) &nbsp; interaction calling<br />
> 7) &nbsp; visualizing the results<br />

**ChIA-PET Tool** was originally published in the journal *Genome Biology* in 2010. After that, the package and its modifications were used in many research projects for publications in high-profile journals. The modifications include revising the linker filtering scripts, adopting the state-of-the-art mapping tools (such as BWA and Bowtie), generating the statistics of the data, and evaluating the quality of the data. In this updated package, we demonstrate how to apply the latest ChIA-PET Tool to the publicly available ChIA-PET data and illustrate the details and interpretation of the results to facilitate the usage of ChIA-PET Tool.
<br />
### ChIA-PET Tool framework

The current ChIA-PET Tool is a command-line program whose execution requires a terminal program. ChIA-PET Tool is mainly coded in Java. Shell scripts are used to glue the different steps in ChIA-PET Tool as a single pipeline. R scripts are used to calculate p-values and generate figures. The package can be downloaded from https://github.com/GuoliangLi-HZAU/ChIA-PET_Tool/archive/master.zip, which includes ChIA-PET Tool source codes in Java, a precompiled JAR file, shell scripts, R scripts and some example files.

> program
>> LGL.jar	&nbsp;&nbsp;&nbsp; // the kernel program written in JAVA<br />
>> LGL  
>>> src	
>>>> LGL	&nbsp;&nbsp;&nbsp;// source codes in JAVA<br />
>>>> path.txt	&nbsp;&nbsp;&nbsp;// necessary file used to compiled source codes<br />

>> hypergeometric.r	&nbsp;&nbsp;&nbsp;// assessing the statistical significance of the interaction from hypergeometric model in R<br />
>> pois.r		&nbsp;&nbsp;&nbsp;// assessing the statistical significance of the peaks from Poisson distribution model in R<br />
>> cutoff_hist_binsize_10bp.r   &nbsp;&nbsp;&nbsp;// assessing the border between self-ligation and inter-ligation PETs <br />
>> peakHeader.txt <br />
>> hg19.chromSize.txt  	&nbsp;&nbsp;&nbsp;// a file contains the length of each chromosome<br />
>> mm10.chromSize.txt <br />
>> ChIA-PET_Tool_Report  	&nbsp;&nbsp;&nbsp;// an empty template for generating report<br />
>> Rscript_and_genome_data  	&nbsp;&nbsp;&nbsp;// R scripts for generating report<br />
>>> ChIA-PET_Tool_Report.r<br />
>>> Plotting_functions.R<br />
>>> hg19_cytoBandIdeo.txt<br />
>>> mm10_cytoBandIdeo.txt


> linker_set_1.with-barcode-info.txt	&nbsp;&nbsp;&nbsp;// linker file 1

> linker_set_2.with-barcode-info.txt  &nbsp;&nbsp;&nbsp;// linker file 2

> MCF7.input.information.txt	&nbsp;&nbsp;&nbsp;// input files and core parameters of linker filtering

> run.MCF7.sh	&nbsp;&nbsp;&nbsp;// shell scripts gluing the whole steps

> deletion.sh   &nbsp;&nbsp;&nbsp;// delete some temperary files and move files for generating visualization report to a new folder

#### If you need to modify the source codes, using the following commands to pack the files. (change your working directory to: ChIA-PET_Tool/program/LGL/src/)
```shell
mkdir ../classes
javac -d ../classes @path.txt
cd ../classes/
jar -cvf LGL.jar LGL
rm ../../LGL.jar 
cp LGL.jar ../../
```
## Installation and download

### Supporting software
ChIA-PET Tool is a pipeline based primarily on JAVA (http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html). At the same time, it also depends on the following softwares. <br />
BWA (http://bio-bwa.sourceforge.net/) is used to map ChIA-PET sequencing reads to a reference genome. BWA can be replaced by other mapping tools, such as Bowtie (http://sourceforge.net/projects/bowtie-bio/files/bowtie). The corresponding mapping tools in the scripts and the genome index should be modified for this purpose. <br /> 
SAMtools (http://samtools.sourceforge.net/) is used to convert the alignment output from SAM format to BAM format. <br />
Bedtools (https://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz) is required to convert the files from BAM format to bedpe format. <br />
R (http://www.r-project.org/) environment is used to compute the p-values and R packages xtable (http://cran.r-project.org/web/packages/xtable/index.html) and RCircos (http://cran.r-project.org/web/packages/RCircos/index.html) are used to generate the graphs for visualization. <br />
Install each software package according to the corresponding instructions and test each software to be run properly.<br />

### Required data
To run ChIA-PET Tool, the genome sequence, chromosome sizes, and cytoband data of the interested genome are required. The genome index needs to be built with BWA (if BWA is used for mapping) in advance.
<br /><br />
In our test, human hg19 reference genome (ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz), chromosome sizes (ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes) and cytoband data (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBandIdeo.txt.gz) were all downloaded from UCSC. So did the mouse data. 

### Example ChIA-PET data
In our test, we used published ChIA-PET data associated with RNA polymerase II (RNAPII) from human breast cell line MCF7 and leukemia cell line K562 (Li *et al.,* 2012), which could be downloaded from ([GEO with accession number GSE33664](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33664)). <br />

## Execution
ChIA-PET Tool is an easy-to-use pipeline and you can simply run it with one command line after you setup all the required tools, data and parameters:<br />
```shell
sh run.MCF7.sh
```
**Before you run the pipeline, you need to modify the variables in the shell scripts, especially for the required tools.** The details of parameters and their meanings can be found in ***user_manual.pdf***. There are different output files from ChIA-PET Tool. The format of the result files and the interpretations of the results are in ***Outputs.docx***, and the running information and summary statistics are shown in ***ChIA-PET_Tool_Report*** generated with ChIA-PET Tool.

## References

* Guoliang Li, Xiaoan Ruan, Raymond K. Auerbach, Michael Snyder, Yijun Ruan, et al. [Extensive Promoter-Centered Chromatin Interactions Provide a Topological Basis for Transcription Regulation](http://www.cell.com/abstract/S0092-8674%2811%2901517-0). **Cell 148(1), 84-98 (2012)** ([Date Sets](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33664))
*  Li G, Fullwood MJ, Xu H et al. [ChIA-PET tool for comprehensive chromatin interaction analysis with paired-end tag sequencing](http://genomebiology.com/2010/11/2/R22). **Genome Biology 11(2):R22 (2010)**
* Fullwood, M. J. et al. [An oestrogen-receptor-alpha-bound human chromatin interactome](http://dx.doi.org/10.1038/nature08497). **Nature 462, 58-64 (2009)**

## Contact us
If you have any problems or suggestions, you could send email to Dr. Guoliang Li (guoliang.li@mail.hzau.edu.cn).
 

### ---END---

##### &nbsp;CopyRight &#169; 2015 Guoliang's Lab, All Rights Reserved





