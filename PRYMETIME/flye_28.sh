#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16

flye --threads ${N_THREADS} --nano-raw "$1" --genome-size "$2" --meta -o "$3"
