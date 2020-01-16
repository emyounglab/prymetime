#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem-per-cpu=10000MB
#SBATCH -o medaka_%j.out
#SBATCH -e medaka_%j.err
#SBATCH -J medaka

source $HOME/medaka/bin/activate
#medaka_consensus -i ~/nanopore_6/fastq/demulti/BC08.fastq -d ~/nanopore_6/FEY_48/scaffolds.fasta -o ~/nanopore_6/FEY_48/FEY_48_Flye_All_med
medaka_consensus -i "$1" -d "$2" -o "$3"
