#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p long
#SBATCH --mem-per-cpu=10000MB
#SBATCH -t 7-00:00:00
#SBATCH -o unicycler_%j.out
#SBATCH -e unicycler_%j.err
#SBATCH -J unicycler

#fail if there's a typo in variable names
set -u
#fail if any command fails
set -e

PREFIX=$(basename "$4")
cd "$4"

if [[ -e "./cir_rep_contigs.fasta" ]]; then

	echo "cir_rep_contigs.fasta exists"

	if [[ -s "./cir_rep_contigs.fasta" ]]; then

	  echo "circular repetitve contigs found, performing Unicycler"

	  cd unicycler

	  for f in *.fasta; do
	    minimap2 -t ${N_THREADS} -ax map-ont "$f" "$1" | samtools fastq --threads ${N_THREADS} -n -f 4 - > "$f"_nano_map.fastq

	    bowtie2-build --threads ${N_THREADS} "$f" "${f}-idx"
	    bowtie2 -x "${f}-idx" --threads ${N_THREADS} --no-unal -1 "$2" -2 "$3" | samtools fastq --threads ${N_THREADS} -n -f 2 -1 "$f"_ill_map_1.fastq -2 "$f"_ill_map_2.fastq -
	    rm -f "${f}-idx"*.bt2

	    if [[ -s "${f}"_ill_map_1.fastq ]]; then

	    unicycler --threads ${N_THREADS} -1 "$f"_ill_map_1.fastq -2 "$f"_ill_map_2.fastq -l "$f"_nano_map.fastq -o "$f"_unicycler

	    fi

	  done

	  cat *_unicycler/assembly.fasta > ../unicycler_contigs.fasta

	  cd ../

	  cat unicycler_contigs.fasta polished_contigs.fasta > "$PREFIX"_comb.fasta

	  # seqkit has threads on by default
	  seqkit seq "$PREFIX"_comb.fasta -m 5000 > "$PREFIX"_comb_filtered.fasta
	  seqkit rename "$PREFIX"_comb_filtered.fasta | seqkit sort --by-length --reverse \
		| awk '/^>/ {if (/circular=true/) \
		{printf(">scaffold_%d_circ\n", ++i)} \
		else {printf(">scaffold_%d\n", ++i)} next} \
		{ print }' "$PREFIX"_comb_filtered.fasta > "$PREFIX"_final.fasta
	else

	  echo "WARNING: no circular contigs, treat only linear"
          seqkit seq polished_contigs.fasta -m 5000  > polished_contigs_filtered.fasta
	  seqkit rename polished_contigs_filtered.fasta | seqkit sort --by-length --reverse | seqkit replace -p '.+' -r 'scaffold_{nr}' > "$PREFIX"_final.fasta

	fi

else
	cd "$4"
	mkdir -p unicycler

	cd unicycler

	echo "WARNING: only linear contigs found, performing Unicycler only"

        unicycler --threads ${N_THREADS} -1 "$2" -2 "$3" -l "$1" -o "$PREFIX"_unicycler

        cat *_unicycler/assembly.fasta > ../unicycler_contigs.fasta

        cd ../

        # seqkit has threads on by default
        seqkit seq unicycler_contigs.fasta -m 5000 > unicycler_contigs_filtered.fasta
        seqkit rename unicycler_contigs_filtered.fasta | seqkit sort --by-length --reverse \
              | awk '/^>/ {if (/circular=true/) \
	      {printf(">scaffold_%d_circ\n", ++i)} \
              else {printf(">scaffold_%d\n", ++i)} next} \
	      { print }' unicycler_contigs_filtered.fasta > "$PREFIX"_final.fasta

fi
