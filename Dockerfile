FROM ubuntu:bionic as build

RUN apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -y install \
      autoconf \
      automake \
      cmake \
      curl \
      gcc \
      git \
      libbz2-dev \
      libcurl4-gnutls-dev \
      liblzma-dev \
      libncurses5-dev \
      libssl-dev \
      make \
      perl \
      python3-pip \
      wget \
      zlib1g-dev \
      && \
    apt-get clean

# Install Flye 2.8
RUN pip3 install git+https://github.com/fenderglass/Flye@2.8

# Install Medaka (https://github.com/nanoporetech/medaka)

# Medaka depends on bgzip, minimap2, samtools, bcftools, and tabix. Install
# those dependencies.

# Install minimap2 (https://github.com/lh3/minimap2)
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 \
    | tar -jxvf - -C /usr/local \
    && ln -s /usr/local/minimap2-2.17_x64-linux/minimap2 /usr/local/bin

# Install packages needed to build HTSlib and samtools
RUN apt-get -y install \
    autoconf \
    automake \
    make \
    gcc \
    perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    && \
    apt-get clean

# Install HTSlib for bgzip, tabix
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
    && bunzip2 -c htslib-1.9.tar.bz2 | tar xf - \
    && cd htslib-1.9 \
    && ./configure \
    && make \
    && make install

# Install samtools 1.9 (requires HTSlib?)
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && bunzip2 -c samtools-1.9.tar.bz2 | tar xf - \
    && cd samtools-1.9 \
    && ./configure \
    && make \
    && make install
    
# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 \
    && bunzip2 -c bcftools-1.10.2.tar.bz2 | tar xf - \
    && cd bcftools-1.10.2 \
    && ./configure \
    && make \
    && make install

# Install idna
RUN pip3 install idna

# Medaka only runs on python 3.5 and python 3.6 as of January,
# 2020. This ties us to Ubuntu 18.04.
RUN pip3 install medaka==1.0.3

# Install racon
RUN git clone --recursive https://github.com/isovic/racon.git racon \
    && cd racon \
    && mkdir build \
    && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release .. \
    && make \
    && make install

# Install Unicycler 0.4.8
RUN pip3 install git+https://github.com/rrwick/Unicycler.git@v0.4.8

# Install fastq_pair
RUN curl -s -L https://github.com/linsalrob/fastq-pair/archive/v1.0.tar.gz | tar xzf - \
    && cd fastq-pair-1.0 \
    && mkdir build \
    && cd build \
    && cmake ../ \
    && make \
    && make install

# Install SPAdes for unicycler
RUN curl -s http://cab.spbu.ru/files/release3.14.0/SPAdes-3.14.0-Linux.tar.gz | tar xzf - \
    && cd SPAdes-3.14.0-Linux \
    && cp bin/* /usr/local/bin/ \
    && cp -r share/* /usr/local/share/

# Install seqkit v0.12.0 in /usr/local/bin
RUN curl -s -L \
    https://github.com/shenwei356/seqkit/releases/download/v0.12.0/seqkit_linux_amd64.tar.gz \
    | tar -xzf - -C /usr/local/bin

# Install Blast+ v2.10.0 in /usr/local
RUN curl -s -L \
    https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/ncbi-blast-2.10.0+-x64-linux.tar.gz \
    | tar -xzf - -C /usr/local \
    && ln -s /usr/local/ncbi-blast-2.10.0+/bin/* /usr/local/bin

# ----------------------------------------------------------------------
#
# Now construct the final docker image without all of the development
# gunk so that it is a leaner docker image.
#
# ----------------------------------------------------------------------
FROM ubuntu:bionic

# Add R package repository
RUN apt -y update && \
    DEBIAN_FRONTEND=noninteractive \
    apt -y --no-install-recommends install \
      gnupg \
      software-properties-common \
    && \
    apt-key adv --keyserver keyserver.ubuntu.com \
      --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

# sorter needs pandas and biopython
RUN apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -y --no-install-recommends install \
      bioperl \
      bowtie2 \
      build-essential \
      bwa \
      default-jdk-headless \
      gfortran \
      libbz2-dev \
      libcurl4-gnutls-dev \
      libjpeg9-dev \
      liblzma-dev \
      libmariadbclient-dev \
      libpng-dev \
      libssl-dev \
      libxml2-dev \
      python3 \
      python3-biopython \
      python3-idna \
      python3-pandas \
      python3-pymummer \
      r-base \
      zlib1g \
      && \
    apt-get clean

COPY --from=build /usr/local /usr/local

#install bedtools
RUN apt-get -y install \
    bedtools \
    && \
    apt-get clean

# Install AliTV

# AliTV needs libyaml-perl libhash-merge-perl bioperl perl cpanm lastz
RUN apt-get -y install \
    libyaml-perl \
    libhash-merge-perl \
    bioperl \
    perl \
    cpanminus \
    && \
    apt-get clean

# Download, compile and install lastz
RUN apt-get -y install \
    wget \
    && \
    apt-get clean

RUN wget https://github.com/lastz/lastz/archive/1.04.00.tar.gz && \
tar -xf 1.04.00.tar.gz && \
cd lastz-1.04.00 && \
# remove -Werror from Makefile to fix compile errors
sed -i 's/-Werror //' src/Makefile && \
make && \
install -m 0755 src/lastz /usr/local/bin/ && \
install -m 0755 src/lastz_D /usr/local/bin/ && \
cd .. && rm -rf lastz-*


# AliTV v1.0.6 install
RUN wget https://github.com/AliTVTeam/AliTV/archive/v1.0.6.tar.gz && \
tar -xf AliTV-1.0.6.tar.gz

# Install chromoMap
RUN R -e "install.packages('chromoMap', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('htmltools', repos = 'http://cran.us.r-project.org')"

# pilon
ADD https://github.com/broadinstitute/pilon/releases/download/v1.23/pilon-1.23.jar /usr/local/bin
ADD pilon /usr/local/bin

# Add the entrypoint script
COPY PRYMETIME /usr/local/bin/prymetime

# When the docker image is launched the resulting container runs the entrypoint script
ENTRYPOINT ["/usr/local/bin/prymetime/PRYMETIME.sh"]
