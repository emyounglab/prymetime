#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p long
#SBATCH --mem-per-cpu=10000MB
#SBATCH -t 7-00:00:00
#SBATCH -o unicycler_%j.out
#SBATCH -e unicycler_%j.err
#SBATCH -J unicycler

#fail if there's a typo in variable names
set -u
#fail if any command fails
set -e

PREFIX=$(basename "$4")
cd "$4"
cd unicycler

for f in *.fasta
do
minimap2 -ax map-ont "$f" "$1" | samtools fastq -n -F 4 - > "$f"_nano_map.fastq

minimap2 -ax sr "$f" "$2" | samtools fastq -n -F 4 - > "$f"_ill_map_1.fastq

minimap2 -ax sr "$f" "$3" | samtools fastq -n -F 4 - > "$f"_ill_map_2.fastq

fastq_pair "$f"_ill_map_1.fastq "$f"_ill_map_2.fastq

unicycler -1 "$f"_ill_map_1.fastq.paired.fq -2 "$f"_ill_map_2.fastq.paired.fq -l "$f"_nano_map.fastq -o "$f"_unicycler
done

cat *_unicycler/assembly.fasta > ../unicycler_contigs.fasta

cd ../

cat unicycler_contigs.fasta polished_contigs.fasta > "$PREFIX"_comb.fasta

seqkit rename "$PREFIX"_comb.fasta | seqkit seq -m 1000 | seqkit sort --by-length --reverse | seqkit replace -p '.+' -r 'scaffold_{nr}' > "$PREFIX"_final.fasta
