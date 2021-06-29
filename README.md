# nanopore_genomic_tutorial



## Preparation
Make sure you have > 10GB of space available

### Download Bandage, the assembly graph viewer

https://rrwick.github.io/Bandage/


### Add conda channels

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

### Make conda environment with all needed packages

    conda create -n nanopore_tutorial flye unicycler minimap2 samtools pilon racon medaka seqkit kraken2 nanoplot
Select yes to install. This will take a while and download about 700MB.

This can be done to detect contamination. 

    kraken2 --use-names --threads 50 --db ~/valentin/kraken2_gtdb/ --report kraken2_gtdb.report knoellia_reads.fastq

### Install pavian in Rstudio

Follow instructions on https://github.com/fbreitwieser/pavian

### Activate conda environment
    conda activate nanopore_tutorial

### Put all files in folder nanopore_tutorial and navigate to folder nanopore_tutorial




## Read QC
Long reads are already basecalled and adapters trimmed (done with guppy basecaller). Short reads are already trimmed (received from sequencing centre). Let's check out the quality of the long reads before assembly.

To check for possible contamination in our sample, you can use kraken2 to classfy all the reads, then check report. The databases required are quite big though (30-140GB), which means that they need a lot of RAM, so we won't do it in this tutorial. The small minikraken databases can be used as an alternative, but only cover the species most abundant in databases. Kraken can then be run like this:

    kraken2 --use-names --threads 50 --db ~/valentin/kraken2_gtdb/ --report kraken2_gtdb.report knoellia_reads.fastq 

I added the resulting report in the repository. Let's check it in pavian.

In RStudio type in

    pavian::runApp(port=5000) 
Then load the report file and visualise the plot in "Samples".

Use seqkit stats to check stats on N50, number of reads, max length, quality scores, etc

    seqkit stats -a knoellia__reads.fastq

Use nanoplot to plot quality vs length (optional)

    NanoPlot --fastq knoellia_reads.fastq -t 8 --N50 --plots -o knoellia_nanoplot


## Hybrid Assembly
Assemble long reads and short reads with unicycler. Should take 30-90 mins on 8 cores.

    unicycler -1 30885_3I2SR_1_trimmed.fastq.gz -2 30885_3ISR_2_trimmed.fastq.gz -l knoellia_reads.fastq -t 8 -o unicycler


## Long Read Assembly

### Long Read Assembly Polishing with Racon (4x) and Medaka (1x)

### Long Read Assembly Polishing with Pilon

### Short Read assembly

## Compare Assemblies
We can compare assemblies using different tools. E.g. seqkit. We can also use Bandage.app to check out the .gfa assembly graph.


