#!/bin/bash


EXECDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"

function usage {
        cat <<EOF
Usage: $0 [-help] [-verbose] [-long <file>] [-illumina_1 <file>] [-illumina_2 <file>]\
 [-outdir <dir>] [-eng_sig <file>] [-ref_genome <file>] [-preferred_assembly]\
 [-read_type]

Processes fastq long reads plus illumina files

	-help		Print Help
	-verbose	Be Verbose
	-long		Specify long read fastq file
	-illumina_1	Paired-end read 1
	-illumina_2	Paired-end read 2
	-outdir		Specify an output directory
	-eng_sig	Fasta file with engineering signatures
	-ref_genome	Reference genome fasta file for comparison
	-preferred_assembly	Tags if organism should be assembled with short-read preference
	-read_type	Tags if providing only short or long reads or assembly to visualize edits
	-long_read_type	Tags if long reads are Nanopore or HiFi PacBio
	-v	Verbose
EOF
}

VERBOSE=no

while [ $# -gt 0 ]; do
    case "$1" in
        -help)		usage; exit 0;;
        -verbose)	VERBOSE=yes;;
	-threads)	shift;N_THREADS="$1";;
	-long)		shift;IN_FASTQ_LONG="$1";;
	-illumina_1)	shift;IN_FASTQ_ILLUMINA_1="$1";;
	-illumina_2)	shift;IN_FASTQ_ILLUMINA_2="$1";;
	-outdir)	shift;OUTDIR="$1";;
	-eng_sig)	shift;ENG_SIG="$1";;
	-ref_genome)	shift;REF_GENOME="$1";;
	-preferred_assembly)	shift;SHORT="$1";;
	-read_type)	shift;READ_TYPE="$1";;
	-long_read_type)	shift;LONG_READ_TYPE="$1";;
        -)     shift; break;;
        -*)
		usage;
                exit 1;;
        *)      break;;         # terminate while loop
    esac
    shift
done


if [[ -z "$IN_FASTQ_LONG" ]]; then
	echo "WARNING: Must specify a nanopore or hifi pacbio fastq file with -long <file>" >&2
fi
if [[ -z "$IN_FASTQ_ILLUMINA_1" ]]; then
	echo "WARNING: Must specify a merged illumina fastq file with -illumina_1 <file>" >&2
fi
if [[ -z "$IN_FASTQ_ILLUMINA_2" ]]; then
	echo "WARNING: Must specify a merged illumina fastq file with -illumina_2 <file>" >&2
fi
if [[ ! -f "$IN_FASTQ_LONG" ]]; then
        echo "WARNING: Long read file '$IN_FASTQ_LONG' does not exist" >&2
fi
if [[ ! -f "$IN_FASTQ_ILLUMINA_1" ]]; then
	echo "WARNING: Illumina file '$IN_FASTQ_ILLUMINA_1' does not exist" >&2
fi
if [[ ! -f "$IN_FASTQ_ILLUMINA_2" ]]; then
	echo "WARNING: Illumina file '$IN_FASTQ_ILLUMINA_2' does not exist" >&2
fi

if [[ -z "$OUTDIR" ]]; then
	OUTDIR=/output
	echo "WARNING: -outdir not specified; outputting to $OUTDIR" >&2
else
	OUTDIR="$OUTDIR"
fi

if [[ -z "$N_THREADS" ]]; then
    N_THREADS=8
    echo "WARNING: -threads not specified; using $N_THREADS" >&2
fi
export N_THREADS

        if [[ "$VERBOSE" = "yes" ]]; then
                echo "Submitting jobs" >&2
                set -v
                set -x
        fi

#fail if there's a typo in variable names
#set -u
#fail if any command fails
set -e

### ENG SIGS ONLY

if [[ "$READ_TYPE" == "assembly" ]]; then

        # if ENG_SIG & REF_GENOME argument was provided, do some more work
        if [[ ! -z "$ENG_SIG" ]]; then

                echo "engineered signatures found"

                $EXECDIR/eng_sig_genome_3.sh "$OUTDIR" "$ENG_SIG"
                $EXECDIR/eng_sig_cmap.sh "$OUTDIR"

                else echo "no engineered signatures found"
        fi
exit 2
fi


### SHORT READS ONLY

if [[ "$READ_TYPE" == "short" ]]; then
        echo "WARNING: only short reads detected, continuing without long reads" >&2

	$EXECDIR/unicycler_short.sh "$IN_FASTQ_ILLUMINA_1" \
        	"$IN_FASTQ_ILLUMINA_2" "$OUTDIR"

	# if ENG_SIG & REF_GENOME argument was provided, do some more work
	if [[ ! -z "$ENG_SIG" ]]; then

		echo "engineered signatures found"

	        $EXECDIR/eng_sig_genome_3.sh "$OUTDIR" "$ENG_SIG"
        	$EXECDIR/eng_sig_cmap.sh "$OUTDIR"

		else echo "no engineered signatures found"
	fi
exit 2
fi

### LONG READS ONLY

if [[ "$READ_TYPE" == "long" ]]; then
        echo "WARNING: only long reads detected, continuing without short reads" >&2

$EXECDIR/flye.sh "$IN_FASTQ_LONG" "$OUTDIR"

$EXECDIR/sorter2.sh "$OUTDIR"

if [[ -s "$OUTDIR/lin_contigs.fasta" ]]; then

        $EXECDIR/medaka.sh "$IN_FASTQ_LONG" "$OUTDIR/lin_contigs.fasta" "$OUTDIR/medaka"

        $EXECDIR/racon_long.sh "$OUTDIR/medaka/consensus.fasta" "$OUTDIR/racon.fasta"

        $EXECDIR/nucmer_long.sh "$OUTDIR"

        $EXECDIR/split.sh "$OUTDIR" "$OUTDIR/cir_rep_contigs.fasta"

        $EXECDIR/unicycler_long.sh "$IN_FASTQ_LONG" "$OUTDIR"

        # if ENG_SIG & REF_GENOME argument was provided, do some more work
        if [[ ! -z "$ENG_SIG" ]]; then

                echo "engineered signatures found"

                $EXECDIR/eng_sig_genome_3.sh "$OUTDIR" "$ENG_SIG"

                $EXECDIR/eng_sig_cmap.sh "$OUTDIR"

        else echo "no engineered signatures found"
        fi
exit 2
else
        #skip medaka, racon, and pilon if no linear contigs
        echo "WARNING: no linear contigs found; continuing with only circular contigs" >&2

        $EXECDIR/nucmer.sh "$OUTDIR"

        $EXECDIR/split.sh "$OUTDIR" "$OUTDIR/cir_rep_contigs.fasta"

        $EXECDIR/unicycler_long.sh "$IN_FASTQ_LONG" "$OUTDIR"

        # if ENG_SIG & REF_GENOME argument was provided, do some more work
        if [[ ! -z "$ENG_SIG" ]]; then

                echo "engineered signatures found"

                $EXECDIR/eng_sig_genome_3.sh "$OUTDIR" "$ENG_SIG"

                $EXECDIR/eng_sig_cmap.sh "$OUTDIR"

        else echo "no engineered signatures found"
        fi
exit 2
fi
fi



### BOTH READS, EITHER PREFERENCE

if [[ -z "$SHORT" ]]; then
	echo "Not tagged as short preferred, continuting to assemble with Flye"
else
	echo "Tagged as short preferred assembly, skipping to Unicycler"

        $EXECDIR/unicycler.sh "$IN_FASTQ_LONG" "$IN_FASTQ_ILLUMINA_1" \
                "$IN_FASTQ_ILLUMINA_2" "$OUTDIR"

        # if ENG_SIG & REF_GENOME argument was provided, do some more work
        if [[ ! -z "$ENG_SIG" ]]; then

                echo "engineered signatures found"

                $EXECDIR/eng_sig_genome_3.sh "$OUTDIR" "$ENG_SIG"

                $EXECDIR/eng_sig_cmap.sh "$OUTDIR"

        fi

exit 2

fi


$EXECDIR/flye.sh "$IN_FASTQ_LONG" "$OUTDIR"

$EXECDIR/sorter2.sh "$OUTDIR"

if [[ -s "$OUTDIR/lin_contigs.fasta" ]]; then

	$EXECDIR/medaka.sh "$IN_FASTQ_LONG" "$OUTDIR/lin_contigs.fasta" "$OUTDIR/medaka"

	$EXECDIR/illumina_merge.sh "$IN_FASTQ_ILLUMINA_1" "$IN_FASTQ_ILLUMINA_2" "$OUTDIR/illumina_merge.fastq"

	$EXECDIR/bowtie2.sh "$OUTDIR/lin_contigs.fasta" "$OUTDIR/illumina_merge.fastq" "$OUTDIR/bowtie2.sam" "$OUTDIR"

	$EXECDIR/racon.sh "$OUTDIR/illumina_merge.fastq" "$OUTDIR/bowtie2.sam" \
		"$OUTDIR/medaka/consensus.fasta" "$OUTDIR/racon.fasta"

	$EXECDIR/pilon.sh "$OUTDIR/racon.fasta" "$OUTDIR/illumina_merge.fastq" \
		"$OUTDIR/pilon.bam" "$OUTDIR/pilon" "$OUTDIR"

	$EXECDIR/nucmer.sh "$OUTDIR"

	$EXECDIR/split.sh "$OUTDIR" "$OUTDIR/cir_rep_contigs.fasta"

	$EXECDIR/unicycler.sh "$IN_FASTQ_LONG" "$IN_FASTQ_ILLUMINA_1" \
		"$IN_FASTQ_ILLUMINA_2" "$OUTDIR"

        # if ENG_SIG & REF_GENOME argument was provided, do some more work
        if [[ ! -z "$ENG_SIG" ]]; then

                echo "engineered signatures found"

                $EXECDIR/eng_sig_genome_3.sh "$OUTDIR" "$ENG_SIG"

                $EXECDIR/eng_sig_cmap.sh "$OUTDIR"

	else echo "no engineered signatures found"
        fi

else
	#skip medaka, racon, and pilon if no linear contigs
	echo "WARNING: no linear contigs found; continuing with only circular contigs" >&2

	$EXECDIR/nucmer.sh "$OUTDIR"

	$EXECDIR/split.sh "$OUTDIR" "$OUTDIR/cir_rep_contigs.fasta"

	$EXECDIR/unicycler.sh "$IN_FASTQ_LONG" "$IN_FASTQ_ILLUMINA_1" \
		"$IN_FASTQ_ILLUMINA_2" "$OUTDIR"

        # if ENG_SIG & REF_GENOME argument was provided, do some more work
        if [[ ! -z "$ENG_SIG" ]]; then

                echo "engineered signatures found"

                $EXECDIR/eng_sig_genome_3.sh "$OUTDIR" "$ENG_SIG"

                $EXECDIR/eng_sig_cmap.sh "$OUTDIR"

        else echo "no engineered signatures found"

        fi
fi
