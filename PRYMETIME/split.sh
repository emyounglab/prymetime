#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o split_%j.out
#SBATCH -e split_%j.err
#SBATCH -J split

#fail if there's a typo in variable names
set -u
#fail if any command fails
set -e

# arg1: output directory
# arg2: circular contigs file

if [[ -s "$2" ]]; then

  # if there are cir_rep_contigs
  cd "$1"
  mkdir -p unicycler
  cp "$2" unicycler/
  cd unicycler
  awk '/^>/{s=++d".fasta"} {print > s}' "$2"
  rm "$2"

else

  echo "WARNING: no circular contigs found; continuing with only linear contigs" >&2

fi
