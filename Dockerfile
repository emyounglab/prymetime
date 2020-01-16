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
      mummer \
      perl \
      python-pip \
      python3-pip \
      python3-pymummer \
      wget \
      zlib1g-dev \
      && \
    apt-get clean

# Install Flye 2.4.2
#
# Use Flye version 2.4.2 per Joe @ WPI. Flye 2.4.2 requires Python
# 2. Can we update to use Flye 2.6 which is a new release that
# supports Python 3?
RUN pip install git+https://github.com/fenderglass/Flye@2.4.2

# Install Medaka (https://github.com/nanoporetech/medaka)

# Medaka depends on bgzip, minimap2, samtools, and tabix. Install
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

# Medaka only runs on python 3.5 and python 3.6 as of January,
# 2020. This ties us to Ubuntu 18.04.
RUN pip3 install medaka


# ----------------------------------------------------------------------
#
# Now construct the final docker image without all of the development
# gunk so that it is a leaner docker image.
#
# ----------------------------------------------------------------------
FROM ubuntu:bionic

# sorter needs pandas and biopython

RUN apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -y --no-install-recommends install \
      libcurl4-gnutls-dev \
      python \
      python3 \
      python3-biopython \
      python3-pandas \
      && \
    apt-get clean

COPY --from=build /usr/local /usr/local

# Add the entrypoint script
COPY PRYMETIME /usr/local/bin/prymetime

# When the docker image is launched the resulting container runs the entrypoint script
ENTRYPOINT ["/usr/local/bin/prymetime/PRYMETIME.sh"]
