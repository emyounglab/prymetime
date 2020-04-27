#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -J eng_sig_blast

cd "$1"

makeblastdb -in "$1"_final.fasta -dbtype nucl

blastn -task blastn -query ../"$2" -db "$1"_final.fasta -perc_identity 98 -qcov_hsp_perc 98 -out "$1"_blastn.txt

perl ../blastn_parse.pl "$1"_blastn.txt > "$1"_blastn.parsed.txt

bp_search2gff.pl --input "$1"_blastn.txt --addid --version 3 --type hit -o "$1"_blastn.gff -f blast --method engineered_region

samtools faidx "$1"_final.fasta

awk 'BEGIN {FS="\t"}; {print $1 FS "1" FS $2}' "$1"_final.fasta.fai > "$1"_final.bed
