#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16

flye --nano-raw "$1" --genome-size "$2" --meta --plasmids -o "$3"
