#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -J eng_sig_blast

EXECDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"

cd "$1"

makeblastdb -in "$1"_final.fasta -dbtype nucl

#eng_sig blast
blastn -task blastn -query ../"$2" -db "$1"_final.fasta -perc_identity 98 -qcov_hsp_perc 98 -out "$1"_blastn_eng_sig.txt

perl ../blastn_parse.pl "$1"_blastn_eng_sig.txt > "$1"_blastn_eng_sig.parsed.txt

bp_search2gff.pl --input "$1"_blastn_eng_sig.txt --addid --version 3 --type hit -o "$1"_blastn_eng_sig.gff -f blast --method eng_sig

#telomere L blast
blastn -task blastn -query ../TELO_L.fasta -db "$1"_final.fasta -out "$1"_blastn_telo_L.txt -max_target_seqs 1 -max_hsps 1

perl ../blastn_parse.pl "$1"_blastn_telo_L.txt > "$1"_blastn_telo_L.parsed.txt

bp_search2gff --input "$1"_blastn_telo_L.txt --addid --version 3 --type hit -o "$1"_blastn_telo_L.gff -f blast --method telomere

#telomere R blast
blastn -task blastn -query ../TELO_R.fasta -db "$1"_final.fasta -out "$1"_blastn_telo_R.txt -max_target_seqs 1 -max_hsps 1

perl ../blastn_parse.pl "$1"_blastn_telo_R.txt > "$1"_blastn_telo_R.parsed.txt

bp_search2gff --input "$1"_blastn_telo_R.txt --addid --version 3 --type hit -o "$1"_blastn_telo_R.gff -f blast --method telomere

#centromere blast
blastn -task blastn -query ../CEN.fasta -db "$1"_final.fasta -out "$1"_blastn_cent.txt -max_target_seqs 1 -max_hsps 1

perl ../blastn_parse.pl "$1"_blastn_cent.txt > "$1"_blastn_cent.parsed.txt

bp_search2gff --input "$1"_blastn_cent.txt --addid --version 3 --type hit -o "$1"_blastn_cent.gff -f blast --method centromere

#mitochondrion blast
blastn -task blastn -query ../MITO.fasta -db "$1"_final.fasta -out "$1"_blastn_mito.txt -max_target_seqs 1 -max_hsps 1

perl ../blastn_parse.pl "$1"_blastn_mito.txt > "$1"_blastn_mito.parsed.txt

bp_search2gff --input "$1"_blastn_mito.txt --addid --version 3 --type hit -o "$1"_blastn_mito.gff -f blast --method mito

#genome_bed
samtools faidx "$1"_final.fasta

awk 'BEGIN {FS="\t"}; {print $1 FS "1" FS $2}' "$1"_final.fasta.fai > "$1"_final.bed

#alitv

#eng_sig alitv
gff2bed < "$1"_blastn_eng_sig.gff > "$1"_blastn_eng_sig.bed

awk -F'\t' -v OFS="\t" '{ print $1, $2, $3, -1, $4 }' "$1"_blastn_eng_sig.bed > eng_sig.tsv

#telomere L alitv
gff2bed < "$1"_blastn_telo_L.gff > "$1"_blastn_telo_L.bed

awk -F'\t' -v OFS="\t" '{ print $1, $2, $3, -1, $4 }' "$1"_blastn_telo_L.bed > telomere_L.tsv

#telomere R alitv
gff2bed < "$1"_blastn_telo_R.gff > "$1"_blastn_telo_R.bed

awk -F'\t' -v OFS="\t" '{ print $1, $2, $3, -1, $4 }' "$1"_blastn_telo_R.bed > telomere_R.tsv

cat telomere_L.tsv telomere_R.tsv > telomere.tsv

#centromere alitv
gff2bed < "$1"_blastn_cent.gff > "$1"_blastn_cent.bed

awk -F'\t' -v OFS="\t" '{ print $1, $2, $3, -1, $4 }' "$1"_blastn_cent.bed > centromere.tsv

#mitochondrion alitv
gff2bed < "$1"_blastn_mito.gff > "$1"_blastn_mito.bed

awk -F'\t' -v OFS="\t" '{ print $1, $2, $3, -1, $4 }' "$1"_blastn_mito.bed > mitochondrion.tsv

#chromomap

#eng_sig chromomap
gff2bed < "$1"_blastn_eng_sig.gff > "$1"_blastn_eng_sig.bed

awk -F'\t' -v OFS="\t" '{ print $4, $1, $2, $3, $8 }' "$1"_blastn_eng_sig.bed > "$1"_blastn_eng_sig_cmap.txt

#telomere L chromomap
gff2bed < "$1"_blastn_telo_L.gff > "$1"_blastn_telo_L.bed

awk -F'\t' -v OFS="\t" '{ print $4, $1, $2, $3, $8 }' "$1"_blastn_telo_L.bed > "$1"_blastn_telo_L_cmap.txt

#telomere R chromomap
gff2bed < "$1"_blastn_telo_R.gff > "$1"_blastn_telo_R.bed

awk -F'\t' -v OFS="\t" '{ print $4, $1, $2, $3, $8 }' "$1"_blastn_telo_R.bed > "$1"_blastn_telo_R_cmap.txt

#centromere chromomap
gff2bed < "$1"_blastn_cent.gff > "$1"_blastn_cent.bed

awk -F'\t' -v OFS="\t" '{ print $4, $1, $2, $3, $8 }' "$1"_blastn_cent.bed > "$1"_blastn_cent_cmap.txt

#mitochondrion chromomap
gff2bed < "$1"_blastn_mito.gff > "$1"_blastn_mito.bed

awk -F'\t' -v OFS="\t" '{ print $4, $1, $2, $3, $8 }' "$1"_blastn_mito.bed > "$1"_blastn_mito_cmap.txt

#combine all chromomap
cat "$1"_blastn_eng_sig_cmap.txt "$1"_blastn_telo_L_cmap.txt "$1"_blastn_telo_R_cmap.txt "$1"_blastn_cent_cmap.txt "$1"_blastn_mito_cmap.txt > "$1"_cmap.txt

