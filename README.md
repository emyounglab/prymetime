# Prymetime

Prymetime is a de novo genome assembly pipeline that uses long reads from Oxford Nanopore Technologies and short reads from Illumina. It was designed to produce high-quality genome assemblies from engineered yeast strains. Prymetime relies on the long read de novo assembler Flye for linear contigs and the hybrid assembler Unicycler for circular contigs.

## Docker

### Build Image

```shell
docker build --tag prymetime:0.2 .
```

### Run container with data

Mount a directory with the `-v` flag. The directory before the `:`
must be an absolute path to a file or directory, and the directory
after the `:` is where it will be mounted inside the container.

```shell
docker run -it --rm \
    -v /path/to/input_dir:/input \
    -v /path/to/output_dir:/output \
    prymetime:0.2 \
    -nanopore /input/my_nanopore.fastq \
    -illumina_1 /input/my_illumina_1.fastq \
    -illumina_2 /input/my_illumina_2.fastq
```



# Run the docker container interactively

The entrypoint script can be overridden for debugging using the
`--entrypoint` argument to docker run. Using `/bin/bash` as the
entrypoint starts an interactive shell when the docker image is
run. Here is an example:

```shell
docker run -it --rm \
    -v $(realpath ../wpi-data):/input \
    -v $(realpath output):/output \
    --entrypoint /bin/bash \
    prymetime:0.2
```

### Detailed PRYMETIME genome assembly pipeline

![PRYMETIME_pipeline](https://github.com/jcollins1x/PRYMETIME/blob/master/docs/PRYMETIME_pipeline_description.png)

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

### Publication

Joseph H. Collins, Kevin W. Keating, Trent R. Jones, Shravani Balaji, Celeste B. Marsan, Marina Como, Zachary J. Newlon, Tom Mitchell, Bryan Bartley, Aaron Adler, Nicholas Roehner, and Eric M. Young. Verification of genetic engineering in yeasts with nanopore whole genome sequencing. In preparation. (2020).
