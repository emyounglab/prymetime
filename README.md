# Prymetime

Prymetime is a de novo genome assembly pipeline that uses long reads from Oxford Nanopore Technologies and short reads from Illumina. It was designed to produce high-quality genome assemblies from engineered yeast strains. Prymetime relies on the long read de novo assembler Flye for linear contigs and the hybrid assembler Unicycler for circular contigs.

All software requirements for Prymetime have been packaged together into a Docker image. Docker is available freely here: https://hub.docker.com/search?offering=community&type=edition

Although it is possible to run the Prymetime Docker image on a desktop computer, we strongly recommend running the pipeline on a server. The memory requirements of Flye and Unicycler at the recommended 40X genome coverage for nanopore and Illumina reads are likely not possible on a "normal" desktop computer.

## Build Docker image

Download Docker image
```shell
git clone https://github.com/emyounglab/prymetime.git
```
Build Docker image
```shell
docker build --tag prymetime prymetime
```

Install time is around one hour on a desktop computer.

## Run Docker image with data

Mount a directory with the `-v` flag. The directory before the `:`
must be an absolute path to a file or directory, and the directory
after the `:` is where it will be mounted inside the container.

Run Prymetime assembly pipeline
```shell
docker run -it --rm \
    -v /path/to/input_dir:/input \
    -v /path/to/output_dir:/output \
    prymetime \
    -nanopore /input/<my_nanopore.fastq> \
    -illumina_1 /input/<my_illumina_1.fastq> \
    -illumina_2 /input/<my_illumina_2.fastq>
    -genome-size <approximate genome size in bp>
    -outdir /output/<my_directory>
```
The final genome assembly will be the <my_directory>_final.fasta file.

Run Prymetime with engineering signatures search
```shell
docker run -it --rm \
    -v /path/to/input_dir:/input \
    -v /path/to/output_dir:/output \
    prymetime \
    -nanopore /input/<my_nanopore.fastq> \
    -illumina_1 /input/<my_illumina_1.fastq> \
    -illumina_2 /input/<my_illumina_2.fastq>
    -genome-size <approximate genome size in bp>
    -outdir /output/<my_directory>
    -eng_sig /input/<my_eng_sig.fasta>
```

The -eng_sig option will also produce a PDF displaying engineering signatures that were found in the genome assembly, shown below:

![FEY_15 engineering signatures](https://github.com/emyounglab/prymetime/blob/master/docs/eng_sig_figure_FEY_15.jpg)

The eng_sig_felix.fasta file (provided in the PRYMETIME folder) contains all engineering signatures used in this study.

## Run Docker image interactively

The entrypoint script can be overridden for debugging using the
`--entrypoint` argument to docker run. Using `/bin/bash` as the
entrypoint starts an interactive shell when the docker image is
run. Here is an example:

```shell
docker run -it --rm \
    -v $(realpath ../data):/input \
    -v $(realpath output):/output \
    --entrypoint /bin/bash \
    prymetime
```

The run time of Prymetime will depend highly on the computer or server used, and the size of the read libraries. On a desktop computer with a small 10X genome coverage read library, Prymetime took approximately 7 hours.

### Detailed PRYMETIME genome assembly pipeline

![PRYMETIME_pipeline](https://github.com/emyounglab/prymetime/blob/master/docs/PRYMETIME_pipeline_description.png)

### Supporting software
Prymetime utilizes the following software packages:
* [Flye](https://github.com/fenderglass/Flye)
* [Medaka](https://github.com/nanoporetech/medaka)
* [Racon](https://github.com/lbcb-sci/racon)
* [Pilon](https://github.com/broadinstitute/pilon)
* [Unicycler](https://github.com/rrwick/Unicycler)
* [Minimap2](https://github.com/lh3/minimap2)
* [BWA](https://github.com/lh3/bwa)
* [Samtools](https://github.com/samtools/samtools)
* [Fastq-pair](https://github.com/linsalrob/fastq-pair)
* [Mummer](https://github.com/mummer4/mummer)
* [karyoploteR](https://https://github.com/bernatgel/karyoploteR)

### Publication

Joseph H. Collins, Kevin W. Keating, Trent R. Jones, Shravani Balaji, Celeste B. Marsan, Marina Como, Zachary J. Newlon, Tom Mitchell, Bryan Bartley, Aaron Adler, Nicholas Roehner, and Eric M. Young. Verification of genetic engineering in yeasts with nanopore whole genome sequencing. In preparation. (2020).
