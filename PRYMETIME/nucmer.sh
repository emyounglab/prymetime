#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o nucmer_%j.out
#SBATCH -e nucmer_%j.err
#SBATCH -J nucmer

cp nucmer4.py "$1"
cd "$1"
python nucmer4.py
