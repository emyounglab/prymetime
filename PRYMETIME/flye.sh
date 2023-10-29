#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16


if [[ "$LONG_READ_TYPE" == "pacbio" ]]; then

	echo "Long reads identified as pacbio hifi"
	flye --threads ${N_THREADS} --min-overlap 5000 --pacbio-hifi "$1" --meta -o "$2"

else

	echo "Long reads identified as nanopore"
	flye --threads ${N_THREADS} --min-overlap 5000 --nano-raw "$1" --meta -o "$2"

fi
