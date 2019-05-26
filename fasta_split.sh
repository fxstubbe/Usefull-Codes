#!/bin/bash

#Genome to split
file="trial.fa"

#Splitting seq by seq
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.final.contigs.fasta
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < $file
