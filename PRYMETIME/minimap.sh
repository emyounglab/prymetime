#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=160G
#SBATCH -o minimap_%j.out
#SBATCH -e minimap_%j.err
#SBATCH -J minimap


#minimap2 -ax sr ~/nanopore_6/FEY_48/FEY_48_Flye_All_med/consensus.fasta ~/illumina_2/FEY48_ill_com.fastq > ~/nanopore_6/FEY_48/FEY_48_Flye_All_med.sam
minimap2 -ax sr "$1"  "$2" > "$3"

