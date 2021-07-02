#!/bin/bash
filein=$1 	#assembly in
longreads=$2	#basecalled reads in
fwd=$3		#illumina fwd
rev=$4		#illunmina rev
counter=1

filename=$(basename -- $filein .fasta)


if [ ! -d ./pilon_polished ]
then
mkdir pilon_polished
fi

if [ ! -d ./sams ]
then
mkdir sams
fi

if [ ! -d ./bams ]
then
mkdir bams
fi

if [ ! -f sams/$filename.sam ]
then

echo "Bowtie mapback of $filename"
bowtie2-build --threads 8 $filein $filename
bowtie2 -x $filename -1 $3 -2 $4 -S sams/$filename.sam --local --threads 8
fi

if [ ! -f ./bams/$filename.sorted.bam ]
then

echo "samtools conversion & sorting of $filename"
samtools view -b -S -F 4 ./sams/$filename.sam -o ./bams/$filename.bam
samtools sort ./bams/$filename.bam -o ./bams/$filename.sorted.bam
samtools index ./bams/$filename.sorted.bam
fi

if [ ! -f ./pilon_polished/"$filename"_"$counter"_pilon.fasta ]
then

echo 'Pilon polishing of $filename'
pilon --genome $filein --frags ./bams/$filename.sorted.bam --output "$filename"_"$counter"_pilon --outdir ./pilon_polished --threads 8 --changes --fix bases
fi

while [ $counter -le 5 ]
do

filein=./pilon_polished/"$filename"_"$counter"_pilon.fasta
let counter=counter+1

echo "Bowtie mapback $counter of $filename"
bowtie2-build --threads 8 $filein "$filename"_"$counter"_pilon
bowtie2 -x "$filename"_"$counter"_pilon -1 $3 -2 $4 -S sams/"$filename"_"$counter"_pilon.sam --local --threads 8

echo "samtools conversion & sorting of "$filename"_"$counter"_pilon"
samtools view -b -S -F 4 sams/"$filename"_"$counter"_pilon.sam -o ./bams/"$filename"_"$counter"_pilon.bam
samtools sort ./bams/"$filename"_"$counter"_pilon.bam -o ./bams/"$filename"_"$counter"_pilon.sorted.bam
samtools index ./bams/"$filename"_"$counter"_pilon.sorted.bam

echo "Pilon polishing of "$filename"_"$counter"_pilon"
pilon --genome $filein --frags ./bams/"$filename"_"$counter"_pilon.sorted.bam --output "$filename"_"$counter"_pilon --outdir ./pilon_polished --threads 8 --changes --fix bases

done
