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

PREFIX=$(basename "$3")

cd "$3"
mkdir -p unicycler

cd unicycler

echo "WARNING: only short contigs found, performing Unicycler only"

unicycler --threads ${N_THREADS} -1 "$1" -2 "$2" -o "$PREFIX"_unicycler
cat *_unicycler/assembly.fasta > ../unicycler_contigs.fasta

cd ../

# seqkit has threads on by default
seqkit seq unicycler_contigs.fasta -m 1000 > unicycler_contigs_filtered.fasta
seqkit rename unicycler_contigs_filtered.fasta | seqkit sort --by-length --reverse \
	| awk '/^>/ {if (/circular=true/) \
	{printf(">scaffold_%d_circ\n", ++i)} \
        else {printf(">scaffold_%d\n", ++i)} next} \
	{ print }' unicycler_contigs_filtered.fasta > "$PREFIX"_final.fasta

