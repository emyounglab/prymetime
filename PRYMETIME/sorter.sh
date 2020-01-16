#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o sorter_%j.out
#SBATCH -e sorter_%j.err
#SBATCH -J sorter_merge

cp sorter.py "$1"
cd "$1"
python sorter.py
