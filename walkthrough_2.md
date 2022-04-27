# nanopore_genomic_tutorial

If you're unsure about the file formats fasta, fastq, sam & bam, check out: https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/file-formats-tutorial/#
### Download Bandage, the assembly graph viewer

https://rrwick.github.io/Bandage/


### Add conda channels

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

### Make conda environment with all needed packages


## Read QC
Long reads are *already* basecalled and adapters trimmed (done with guppy basecaller). The read files were put into one file using the cat command.
To put all reads of one barcode in one file, navigate to the barcode folder and use `cat ./*.fastq > barcode_xy_all_reads.fastq`. For this example, all reads are already in one file.
Short reads (forward and reverse) are already trimmed (received from sequencing centre). 

Let's check out the quality of the long reads before assembly.

To check for possible contamination in our sample, you can use kraken2 to classify all the reads, then check report - how many % are E. coli?

    kraken2 --db ~/software/kraken_db/ --threads 6 --output barcode_18.output --report barcode_18.report --use-names ~/nanopore_data/fastq_pass/barcode18.fastq

You can also visualise the reports using pavian in Rstudio on your laptop (if you have it).

<details>
<summary>Click to expand</summary>

whatever

In RStudio type in

    pavian::runApp(port=5000) 
Then load the report file and visualise the plot in "Samples". We can see that there is some contamination with *Streptomyces* in the long reads file, but the short reads seem fine.
</details>

Use seqkit stats to check stats on N50, number of reads, max length, quality scores, etc

    seqkit stats -a knoellia__reads.fastq


## Assemblies
We're going to do short-read assembly, hybrid assembly and long-read assembly followed by polishing with short reads (i.e. kind of hybrid).

## Short Read assembly with SPAdes
SPAdes just uses the short reads to assemble. Very accurate, maybe not so contiguous - let's see for ourselves. Main output is scaffolds.fasta.

    spades.py --isolate -1 ./30885_3I2SR_1_trimmed.fastq.gz -2 ./30885_3I2SR_2_trimmed.fastq.gz -t 8 -o spades_assembly


## Hybrid Assembly with Unicycler
Unicycler uses spades to assemble the short read contigs, but then also uses the long reads to bridge the gaps in between, and then polishes off the gaps with short reads. Result should be much more contiguous than short read assembly. Should take 30-90 mins on 8 cores, main output is assembly.fasta

    unicycler -1 30885_3I2SR_1_trimmed.fastq.gz -2 30885_3ISR_2_trimmed.fastq.gz -l knoellia_reads.fastq -t 8 -o unicycler


## Long Read Assembly with Flye
Assemble long reads with flye. This usually produces the most contiguous assembly, but it still needs polishing both with long-read specific tools and with the short reads of the same sample. This should take 20 minutes or so.

    flye --threads 8 --genome-size 4m --nano-raw knoellia_reads.fastq --out-dir flye_assembly

You can inspect the assembly outcome by checking assembly_info.txt. 

    less flye_assembly/assembly_info.txt
The slightly contaminated long reads combined with flye's ability to assemble contigs at low coverage have lead to the presence several other fragments apart from the circular contig.

### Long Read Assembly Polishing with Racon (4x) and Medaka (1x)

~~The common workflow is to use racon in 4 iterations before feeding the resulting file into medaka. Racon uses the basecalled reads to find a "most likely" correct base. Medaka uses a machine learning model trained on genomic data to correct assemblies.~~

This has been superseded - _new versions of medaka now work with the output of flye directly._

~~#### Racon
First, we need to map the reads to the assembly, creating a .paf (similar to sam and bam) file. Then we polish with racon. This takes minutes.

    mkdir racon
    minimap2 -x map-ont -t 8 ./flye_assembly/assembly.fasta knoellia_reads.fastq > ./racon/knoellia.racon.1.paf
    
~~Then we polish with racon.~~

    racon -u -t 8 -m 8 -x -6 -g -8 -w 500 knoellia_reads.fastq ./racon/knoellia.racon.1.paf ./flye_assembly/assembly.fasta > knoellia.racon.1.fasta

~~Then we repeat the process three more times. You can automate this with the script provided in the files.~~

#### Medaka
Medaka is a machine-learning model based polisher, so it is super important to correctly specific the flow cell (R9.41) and the basecaller used (guppy 3.2.2 high accuracy which has the same model as guppy 3.03), because that will influence the way that medaka corrects the assembly. Main output is consensus.fasta in the medaka folder.
    
    medaka_consensus -d ./racon/assembly.racon.4.fasta -i knoellia_reads.fastq -m r941_min_high_g303 -t 8 -o medaka


### Long Read Assembly Polishing with Pilon
After medaka, the assembly should be pretty accurate already. Now polish with short reads it until it's done. This could be anything between 1 and 4 rounds of pilon polishing.

First, map the short reads to the medaka output with bowtie2. 

    mkdir pilon
    cd pilon
    bowtie2-build --threads 8 ../medaka/consensus.fasta knoellia_pilon_1
    bowtie2 -x knoellia_pilon_1 -1 ../30885_3I2SR_1_trimmed.fastq.gz -2 ../30885_3I2SR_2_trimmed.fastq.gz -S knoellia_pilon_1.sam --local --threads 8


Then, convert the .sam format file to .bam format file, while only keeping mapped reads (-F 4); then sort & index.

    samtools view -@ 8 -b -S -F 4 knoellia_pilon_1.sam -o knoellia_pilon_1.bam
    samtools sort -@ 8 knoellia_pilon_1.bam -o knoellia_pilon_1.sorted.bam
    samtools index -@ 8 knoellia_pilon_1.sorted.bam

Polish with pilon using the sorted & indexed bam. This might also take a while. If CPU/RAM usage is a problem, try adding the `--fix bases` flag. It limits the amount of errors is corrects, but reduces computational cost.

    pilon --genome ../medaka/consensus.fasta --frags knoellia_pilon_1.sorted.bam --output knoellia_pilon_1 --outdir . --threads 8 --changes 
    
ONLY IF a `java.lang.OutOfMemoryError: Java heap space` comes up, do the following: Type `which pilon`. This will point you to the place of the actual executable, and give you something like `/home/valentin/miniconda3/envs/nanopore_tutorial/bin/pilon`. Edit the file by typing `nano /home/valentin/miniconda3/envs/nanopore_tutorial/bin/pilon (or whatever your path is)`. Then find the line that specifies the memory limits: `default_jvm_mem_opts = ['-Xms512m', '-Xmx1g']` and change the `-Xmx1g` into `-Xmx8g`. Save by pressing ctrl+X then type 'y'. Then try running the pilon command again. 

The output should be a polished genome. Thanks to the `--changes` flag, pilon also generates a file detailing the changes it made in each round. When this file is empty, it means that pilon has finished and no further rounds are necessary.

There is a script in the files that automates six rounds of mapping & polishing, but check the script for correct flags etc before running. No guarantees it will work.

#### Select the contigs
Since we saw from the kraken report on the long reads, there was some contamination of the reads. In the assembly_info.txt generated by flye, we could see that there are several short, non-circular fragments with low coverage in addition to the large, circular contig that is the chromosome. The short contigs are likely from the contaminated reads. Therefore, in order to only carry forward the genome that we want, we can remove the others with seqkit. In this case, since the closed genome is the only large contig, a simple selection by size (i.e. select all contigs > 1MB) works. If you want to be sure that it's not an actual real plasmid or something we are removing, you can also copy and blast some of the short sequences to check where they come from.

    seqkit seq -m 1000000 consensus_6_pilon.fasta -o knoellia_flye_pilon6_filtered.fasta

## Compare Assemblies
We can compare assemblies using different tools. E.g. a very obvious one using seqkit to check min & max sequence lengths, N50s, etc.

    seqkit stats -a spades_assembly/scaffolds.fasta
    seqkit stats -a unicycler/assembly.fastaknoellia_flye_pilon6_filtered.fasta
    seqkit stats -a pilon/pilon_polished/knoellia_flye_pilon6_filtered.fasta
    
We can also use Bandage.app to check out the .gfa assembly graphs that spades, unicycler and flye produce.

### BUSCO  
Busco uses single copy core genes (SCCGs) to assess completeness and assembly quality. We can use BUSCO to assess the different assemblies we created, as well as to show the improvement of the polishing steps in the course of the nanopore assembly.

First, select the appropriate dataset listed by BUSCO that fits with the phylogeny of our *Knoellia* genome. We should use the most specific dataset we can, so we can list all available datasets and choose the most specific one:

    busco --list-datasets
    
Another option is to use the --auto-lineage flag, then BUSCO will try selecting the best lineage itself. Alternatives to BUSCO, especially for MAGs, include checkm (more computing power required).

The result will show how many of the expected SCCGs are present & complete (C), if they are fragmented (F), or if they are missing completely (M). This indicates whether a genome is complete or a part is missing; it also gives a measure of sequence quality, since lots of fragmented SCCGs mean high rates of indels.
    
    
    mkdir busco
    cd busco
    busco -m genome -i ../flye_assembly/assembly.fasta -o flye_only -l micrococcales_odb10 -c 8
    busco -m genome -i ../pilon/pilon_polished/knoellia_flye_pilon6_filtered.fasta -o pilon_6 -l micrococcales_odb10 -c 8
    busco -m genome -i ../unicycler/assembly.fasta -o unicycler -l micrococcales_odb10 -c 8
    busco -m genome -i ../spades_assembly/scaffolds.fasta -o spades -l micrococcales_odb10 -c 8
    busco -m genome -i ../medaka/consensus.fasta -o medaka -l micrococcales_odb10 -c 8
    
Results:

    Flye only:
    C:67.1%[S:53.3%,D:13.8%],F:19.2%,M:13.7%,n:537 
    
    After medaka:
    C:76.9%[S:56.8%,D:20.1%],F:12.8%,M:10.3%,n:537 
    
    After pilon:
    C:98.5%[S:98.3%,D:0.2%],F:0.9%,M:0.6%,n:537 
    
    Unicycler:
    C:98.0%[S:97.8%,D:0.2%],F:1.1%,M:0.9%,n:537  
    
    SPAdes:
    C:99.1%[S:98.9%,D:0.2%],F:0.4%,M:0.5%,n:537  

From these results, it is  evident how much short read polishing improves the long-read assembly, even after medaka polishing. If you have read any release notes on these tools, you might think that medaka should produce much better results. This is likely due to the fact that the medaka model provided by nanopore is trained only on E. coli, human genome and yeast. Therefore it will only produce the best results for these organisms.

The contamination can also be assessed from these numbers: the % duplicated (D) shows how many of the SCCGs are present twice. in the flye & medaka assembly, we can see that there are a lot of duplicated SCCGs, while that is not the case afterwards. *This is because after the pilon polishing we actively selected to keep only the one circular large circular sequence, and the duplicated genes came from the contaminations*. 

We can also see that SPAdes produces the most accurate assembly (highest completeness), closely followed by the pilon-polished flye assembly and unicycler. So while it seems that SPAdes wins in terms of accuracy, it is important to keep in mind the contiguity of the assembly - the largest fragment in the SPAdes assembly is 132kb long, while both unicycler and flye assemble the complete chromosome. So it seems like a trade-off between a slight decrease in accuracy and a high increase in assembly contiguity.

