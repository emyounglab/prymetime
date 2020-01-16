#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem-per-cpu=20000MB
#SBATCH -o flye_%j.out
#SBATCH -e flye_%j.err
#SBATCH -J flye

# source activate flye_24

#flye --nano-raw ~/nanopore_6/fastq/demulti/BC05.fastq --genome-size 12100000 --plasmids -o ~/nanopore_6/FEY_27
flye --nano-raw "$1" --genome-size "$2" --meta --plasmids -o "$3"



