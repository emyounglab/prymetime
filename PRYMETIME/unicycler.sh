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

if [[ -s "$OUTDIR/cir_rep_contigs.fasta" ]]; then

  cd unicycler

  for f in *.fasta; do
    minimap2 -t ${N_THREADS} -ax map-ont "$f" "$1" | samtools fastq --threads ${N_THREADS} -n -F 4 - > "$f"_nano_map.fastq

    bowtie2-build --threads ${N_THREADS} "$f" refgenome
    bowtie2 -x refgenome --threads ${N_THREADS} --no-unal -1 "$2" -2 "$3" | samtools fastq --threads ${N_THREADS} -n -f 2 -1 "$f"_ill_map_1.fastq -2 "$f"_ill_map_2.fastq -
    rm -f refgenome*.bt2

    unicycler --threads ${N_THREADS} -1 "$f"_ill_map_1.fastq -2 "$f"_ill_map_2.fastq -l "$f"_nano_map.fastq -o "$f"_unicycler
  done

  cat *_unicycler/assembly.fasta > ../unicycler_contigs.fasta

  cd ../

  cat unicycler_contigs.fasta polished_contigs.fasta > "$PREFIX"_comb.fasta

  # seqkit has threads on by default
  seqkit rename "$PREFIX"_comb.fasta | seqkit seq -m 1000 | seqkit sort --by-length --reverse | seqkit replace -p '.+' -r 'scaffold_{nr}' > "$PREFIX"_final.fasta
else
  # if there are no cir_rep_contigs, treat only linear
  seqkit rename polished_contigs.fasta | seqkit seq -m 1000 | seqkit sort --by-length --reverse | seqkit replace -p '.+' -r 'scaffold_{nr}' > "$PREFIX"_final.fasta
fi
