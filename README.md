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
Select yes to install. This will take a while and download about 500MB.

### Activate conda environment
    conda activate nanopore_tutorial

### Put all files in folder nanopore_tutorial

## Read QC
Long reads are already basecalled and adapters trimmed (done with guppy basecaller). Short reads are already trimmed (received from sequencing centre). Let's check out the quality of the long reads before assembly.

Use seqkit stats to check stats on N50, number of reads, max length, etc.

    seqkit stats -a knoellia_long_reads.fastq

Use nanoplot to plot quality vs length. You will notice that the reads are all above 7000, because I removed everything smaller than that.



## Hybrid Assembly
Assemble long reads and short reads with unicycler.

## Long Read Assembly

### Long Read Assembly Polishing with Racon (4x) and Medaka (1x)

### Long Read Assembly Polishing with Pilon

### Short Read assembly

## Compare Assemblies
We can compare assemblies using different tools. E.g. seqkit. We can also use Bandage.app to check out the .gfa assembly graph.


