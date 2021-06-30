#!/bin/bash

assembly=$1
reads=$2

#$1 is the assembly.
#$2 is the reads fastq file.
shortname=$(basename -- "$assembly" .fasta)
counter=1

while [ $counter -le 4 ]
do
	if [ ! -f ./racon"$shortname".racon."$counter".paf ];
	then
		echo "Minimap mapback $counter"
		minimap2 -x map-ont -t 8 $assembly $reads > ./racon/"$shortname".racon."$counter".paf

	else
		echo "./racon/"$shortname".racon."$counter".paf already exists. Skipping mapping step"
	fi

	if [ ! -f ./racon/"$shortname".racon."$counter".fasta ];
	then
	       echo "Racon iteration $counter"
	      	racon -u -t 8 -m 8 -x -6 -g -8 -w 500 $reads ./racon/"$shortname".racon."$counter".paf $assembly > ./racon/"$shortname".racon."$counter".fasta
		assembly="./racon/$shortname".racon."$counter".fasta
		rm ./racon/"$shortname".racon."$counter".paf
	else
		echo "./racon/"$shortname".racon."$counter".fasta already exists. Skipping correction step"
	fi
		let counter=counter+1
done
