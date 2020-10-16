#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem-per-cpu=20000MB
#SBATCH -o flye_%j.out
#SBATCH -e flye_%j.err
#SBATCH -J flye

flye --nano-raw "$1" --genome-size "$2" --meta --plasmids -o "$3" --threads 128
