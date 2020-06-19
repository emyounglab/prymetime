#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -J eng_sig_blast

EXECDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"

cd "$1"
PREFIX=$(basename "$1")

makeblastdb -in "$PREFIX"_final.fasta -dbtype nucl

blastn -task blastn -query "$2" -db "$PREFIX"_final.fasta -perc_identity 98 -qcov_hsp_perc 98 -out "$PREFIX"_blastn.txt

perl "$EXECDIR"/blastn_parse.pl "$PREFIX"_blastn.txt > "$PREFIX"_blastn.parsed.txt

bp_search2gff --input "$PREFIX"_blastn.txt --addid --version 3 --type hit -o "$PREFIX"_blastn.gff -f blast --method engineered_region

samtools faidx "$PREFIX"_final.fasta

awk 'BEGIN {FS="\t"}; {print $1 FS "1" FS $2}' "$PREFIX"_final.fasta.fai > "$PREFIX"_final.bed
