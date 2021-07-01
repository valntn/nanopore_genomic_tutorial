# nanopore_genomic_tutorial



## Preparation
Make sure you have > 10GB of space available
General note: I have used all commands with `--threads 8 / -t 8`. This means that the command will use 8 of your CPU threads. If your computer only has 8 threads available (e.g. Macbooks), this might make everything else VERY SLOW. Modify the multithreading by using e.g. 6 or 7 instead of 8.

If you're unsure about the file formats fasta, fastq, sam & bam, check out: https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/file-formats-tutorial/#
### Download Bandage, the assembly graph viewer

https://rrwick.github.io/Bandage/


### Add conda channels

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

### Make conda environment with all needed packages

    conda create -n nanopore_tutorial flye unicycler minimap2 samtools pilon racon medaka seqkit kraken2 nanoplot bowtie2
Select yes to install. This will take a while and download about 700MB.


### Install pavian in Rstudio

Follow instructions on https://github.com/fbreitwieser/pavian

### Activate conda environment
    conda activate nanopore_tutorial

### Put all files in folder nanopore_tutorial and navigate to folder nanopore_tutorial




## Read QC
Long reads are already basecalled and adapters trimmed (done with guppy basecaller). To put all reads of one barcode in one file, navigate to the barcode folder and use `cat ./*.fastq > all_reads.fastq`. For this example, all reads are already in one file.
Short reads (forward and reverse) are already trimmed (received from sequencing centre). 

Let's check out the quality of the long reads before assembly.

To check for possible contamination in our sample, you can use kraken2 to classify all the reads, then check report. The databases required are quite big though (30-140GB), which means that they need a lot of RAM, so we won't do it in this tutorial. The small minikraken databases can be used as an alternative, but only cover the species most abundant in databases. Kraken can then be run like this:

    kraken2 --use-names --threads 50 --db ~/valentin/kraken2_gtdb/ --report kraken2_gtdb.report knoellia_reads.fastq 

I added the resulting report in the repository. Let's check it in pavian.

In RStudio type in

    pavian::runApp(port=5000) 
Then load the report file and visualise the plot in "Samples".

Use seqkit stats to check stats on N50, number of reads, max length, quality scores, etc

    seqkit stats -a knoellia__reads.fastq

Use nanoplot to plot quality vs length (optional)

    NanoPlot --fastq knoellia_reads.fastq -t 8 --N50 --plots -o knoellia_nanoplot

## Assemblies
We're going to do short-read assembly, hybrid assembly and long-read assembly followed by polishing with short reads (so kind of hybrid).

## Short Read assembly with SPAdes
SPAdes just uses the short reads to assemble. Very accurate, maybe not so contiguous - let's see for ourselves. Main output if scaffolds.fasta.

    spades.py --isolate -1 ./30885_3I2SR_1_trimmed.fastq.gz -2 ./30885_3I2SR_2_trimmed.fastq.gz -t 8 -o spades_assembly


## Hybrid Assembly with Unicycler
Unicycler uses spades to assemble the short read contigs, but then also uses the long reads to bridge the gaps in between, and then polishes off the gaps with short reads. Result should be much more contiguous than short read assembly. Should take 30-90 mins on 8 cores, main output is assembly.fasta

    unicycler -1 30885_3I2SR_1_trimmed.fastq.gz -2 30885_3ISR_2_trimmed.fastq.gz -l knoellia_reads.fastq -t 8 -o unicycler


## Long Read Assembly with Flye
Assemble long reads with flye. This usually produces the most contiguous assembly, but it still needs polishing both with long-read specific tools and with the short reads of the same sample. This should take 20 minutes or so.

    flye --threads 8 --genome-size 4g --nano-raw knoellia_reads.fastq --out-dir flye_assembly

You can inspect the assembly outcome by checking assembly_info.txt. 

    less flye_assembly/assembly_info.txt
The slightly contaminated long reads combined with flye's ability to assemble contigs at low coverage have lead to the presence several other fragments apart from the circular contig.

### Long Read Assembly Polishing with Racon (4x) and Medaka (1x)
The common workflow is to use racon in 4 iterations before feeding the resulting file into medaka. Racon uses the basecalled reads to find a "most likely" correct base. Medaka uses a machine learning model trained on genomic data to correct assemblies.

#### Racon
First, we need to map the reads to the assembly, creating a .paf (similar to sam and bam) file. Then we polish with racon. This takes minutes.

    mkdir racon
    minimap2 -x map-ont -t 8 ./flye_assembly/assembly.fasta knoellia_reads.fastq > ./racon/knoellia.racon.1.paf
    
Then we polish with racon.

    racon -u -t 8 -m 8 -x -6 -g -8 -w 500 knoellia_reads.fastq ./racon/knoellia.racon.1.paf ./flye_assembly/assembly.fasta > knoellia.racon.1.fasta

Then we repeat the process three more times. You can automate this with the script provided in the files.

#### Medaka
Medaka is a machine-learning model based polisher, so it is super important to correctly specific the flow cell (R9.41) and the basecaller used (guppy 3.2.2 high accuracy which has the same model as guppy 3.03), because that will influence the way that medaka corrects the assembly. Main output is consensus.fasta in the medaka folder.
    
    medaka_consensus -i ./racon/assembly.racon.4.fasta -d knoellia_reads.fastq -m r941_min_high_g303 -t 8 -o medaka


### Long Read Assembly Polishing with Pilon
After medaka, the assembly should be pretty accurate already. Now polish with short reads it until it's done. This could be anything between 1 and 4 rounds of pilon polishing.

First, map the short to the medaka output with bowtie2.

    mkdir pilon
    cd pilon
    bowtie2-build --threads 8 ../medaka/consensus.fasta knoellia_pilon_1
    bowtie2 -x knoellia_pilon_1 -1 ../30885_3I2SR_1_trimmed.fastq.gz -2 ../30885_3I2SR_2_trimmed.fastq.gz -S knoellia_pilon_1.sam --local --very-sensitive-local -I 0 -X 1000 --threads 8


Then, convert the .sam format file to .bam format file, while only keeping mapped reads (-F 4)

Polish with pilon.

## Compare Assemblies
We can compare assemblies using different tools. E.g. seqkit:

    seqkit stats -a spades_assembly/scaffolds.fasta
    seqkit stats -a unicycler/assembly.fasta
    
    
We can also use Bandage.app to check out the .gfa assembly graphs that spades, unicycler and flye produce.


