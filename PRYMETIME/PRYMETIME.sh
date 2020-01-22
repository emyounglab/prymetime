#!/bin/bash


EXECDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"

function usage {
        cat <<EOF
Usage: $0 [-help] [-verbose] [-nanopore <file>] [-illumina_1 <file>] [-illumina_2 <file>] [-outdir <dir>] [-genome-size <size>]
Processes fastq nanopore plus illumina files

	-help		Print Help
	-verbose	Be Verbose
	-nanopore	Specify nanopore fastq file
	-illumina_1	Paired-end read 1
	-illumina_2	Paired-end read 2
	-outdir		Specify an output directory
	-genome-size	Specify genome size
	-v	Verbose
EOF
}

VERBOSE=no

while [ $# -gt 0 ]; do
    case "$1" in
        -help)		usage; exit 0;;
        -verbose)	VERBOSE=yes;;
	-nanopore)	shift;IN_FASTQ_NANOPORE="$1";;
	-illumina_1)	shift;IN_FASTQ_ILLUMINA_1="$1";;
	-illumina_2)	shift;IN_FASTQ_ILLUMINA_2="$1";;
	-outdir)	shift;OUTDIR="$1";;
	-genome-size)	shift;GENOME_SIZE="$1";;
        -)     shift; break;;
        -*)
		usage;
                exit 1;;
        *)      break;;         # terminate while loop
    esac
    shift
done

if [[ -z "$IN_FASTQ_NANOPORE" ]]; then
	echo "ERROR: Must specify a nanopore fastq file with -nanopore <file>" >&2
	exit 2
fi
if [[ -z "$IN_FASTQ_ILLUMINA_1" ]]; then
	echo "ERROR: Must specify a merged illumina fastq file with -illumina_1 <file>" >&2
	exit 2
fi
if [[ -z "$IN_FASTQ_ILLUMINA_2" ]]; then
	echo "ERROR: Must specify a merged illumina fastq file with -illumina_2 <file>" >&2
	exit 2
fi
if [[ ! -f "$IN_FASTQ_NANOPORE" ]]; then
	echo "ERROR: Nanopore file '$IN_FASTQ_NANOPORE' does not exist" >&2
	exit 2
fi
if [[ ! -f "$IN_FASTQ_ILLUMINA_1" ]]; then
	echo "ERROR: Illumina file '$IN_FASTQ_ILLUMINA_1' does not exist" >&2
	exit 2
fi
if [[ ! -f "$IN_FASTQ_ILLUMINA_2" ]]; then
	echo "ERROR: Illumina file '$IN_FASTQ_ILLUMINA_2' does not exist" >&2
	exit 2
fi

if [[ -z "$OUTDIR" ]]; then
	OUTDIR=/output
	echo "WARNING: -outdir not specified; outputting to $OUTDIR" >&2
else
	OUTDIR="$OUTDIR"
fi

if [[ -z "$GENOME_SIZE" ]]; then
	GENOME_SIZE=12100000
	echo "WARNING: -genome-size not specified; using $GENOME_SIZE" >&2
fi

#fail if there's a typo in variable names
set -u
#fail if any command fails
set -e

#REAL PART OF SCRIPT BEGINS HERE

if [[ "$VERBOSE" = "yes" ]]; then
	echo "Submitting jobs" >&2
    set -v
    set -x
fi

$EXECDIR/flye.sh "$IN_FASTQ_NANOPORE" $GENOME_SIZE "$OUTDIR"

$EXECDIR/sorter.sh "$OUTDIR"

$EXECDIR/medaka.sh "$IN_FASTQ_NANOPORE" "$OUTDIR/lin_contigs.fasta" "$OUTDIR/medaka"

$EXECDIR/illumina_merge.sh "$IN_FASTQ_ILLUMINA_1" "$IN_FASTQ_ILLUMINA_2" "$OUTDIR/illumina_merge.fastq"

$EXECDIR/minimap.sh "$OUTDIR/medaka/consensus.fasta" "$OUTDIR/illumina_merge.fastq" "$OUTDIR/minimap.sam"

$EXECDIR/racon.sh "$OUTDIR/illumina_merge.fastq" "$OUTDIR/minimap.sam" \
    "$OUTDIR/medaka/consensus.fasta" "$OUTDIR/racon.fasta"

$EXECDIR/pilon.sh "$OUTDIR/racon.fasta" "$OUTDIR/illumina_merge.fastq" \
    "$OUTDIR/pilon.bam" "$OUTDIR/pilon" "$OUTDIR"

$EXECDIR/nucmer.sh "$OUTDIR"

$EXECDIR/split.sh "$OUTDIR"

$EXECDIR/unicycler.sh "$IN_FASTQ_NANOPORE" "$IN_FASTQ_ILLUMINA_1" \
    "$IN_FASTQ_ILLUMINA_2" "$OUTDIR"
