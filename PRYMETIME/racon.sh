#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=160G
#SBATCH -o racon_%j.out
#SBATCH -e racon_%j.err
#SBATCH -J racon

#racon ~/illumina_2/FEY48_ill_com.fastq ~/nanopore_6/FEY_48/FEY_48_Flye_All_med.sam ~/nanopore_6/FEY_48/FEY_48_Flye_All_med/consensus.fasta > ~/nanopore_6/FEY_48/FEY_48_Flye_All_med_rac.fasta
racon --threads ${N_THREADS} "$1" "$2" "$3" > "$4"

bwa index "$4"
