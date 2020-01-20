#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o nucmer_%j.out
#SBATCH -e nucmer_%j.err
#SBATCH -J nucmer

EXECDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"
cp "${EXECDIR}/nucmer4.py" "$1"
cd "$1"
python3 nucmer4.py
