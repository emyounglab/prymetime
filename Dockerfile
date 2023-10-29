FROM ubuntu:focal as build

#ENV CUDO_VISIBLE_DEVICES=-1

# Set locale settings
ENV LANGUAGE=en_US.en
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
RUN apt-get update && apt-get install -y locales && \
        sed -i -e "s/# $LANG.*/$LANG UTF-8/" /etc/locale.gen && \
        dpkg-reconfigure --frontend=noninteractive locales && \
        update-locale LANG=$LANG

# Generate locale settings
RUN locale-gen en_US.UTF-8

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
      libtool \
      make \
      wget \
      apt-utils \
      zip \
      zlib1g-dev \
      yasm \
      build-essential \
      python3-dev \
      python3-pip \
      zlib1g-dev \
      libncursesw5-dev \
      gfortran \
      libreadline8 \
      libreadline-dev \
      libx11-dev \
      libxt6 \
      xorg-dev \
      libpcre2-posix2 \
      libpcre2-dev \
      && \
    apt-get clean

# downgrade numpy for deprecated np.bool
RUN pip3 install numpy==1.23.1

# Install Flye 2.9
RUN pip3 install git+https://github.com/fenderglass/Flye@2.9.2

# Install Medaka (https://github.com/nanoporetech/medaka)

# Medaka depends on bgzip, minimap2, samtools, bcftools, and tabix. Install
# those dependencies.

# Install minimap2 (https://github.com/lh3/minimap2)
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 \
    | tar -jxvf - -C /usr/local \
    && ln -s /usr/local/minimap2-2.24_x64-linux/minimap2 /usr/local/bin

# Install bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
RUN curl -L -o bowtie2-2.4.5-linux-x86_64.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.5/bowtie2-2.4.5-linux-x86_64.zip/download \
    && unzip bowtie2-2.4.5-linux-x86_64.zip -d /usr/local \
    && ln -s /usr/local/bowtie2-2.4.5-linux-x86_64/bowtie2 /usr/local/bin \
    && ln -s /usr/local/bowtie2-2.4.5-linux-x86_64/bowtie2-build /usr/local/bin

# Install packages needed to build HTSlib and samtools
RUN apt-get -y install \
    autoconf \
    automake \
    make \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    && \
    apt-get clean

# Install HTSlib for bgzip, tabix
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 \
    && bunzip2 -c htslib-1.17.tar.bz2 | tar xf - \
    && cd htslib-1.17 \
    && ./configure \
    && make \
    && make install

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 \
    && bunzip2 -c samtools-1.17.tar.bz2 | tar xf - \
    && cd samtools-1.17 \
    && ./configure \
    && make \
    && make install

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 \
    && bunzip2 -c bcftools-1.17.tar.bz2 | tar xf - \
    && cd bcftools-1.17 \
    && ./configure \
    && make \
    && make install

# Install idna
RUN pip3 install idna

# Install Cython
RUN pip3 install --upgrade cython

# Install Setuptools
RUN pip3 install setuptools-scm==6.4.2

# Install protobuf
RUN pip3 install protobuf==4.22.1

# Install Medaka
RUN pip3 install medaka==1.7.3

# Install racon
RUN git clone --recursive https://github.com/isovic/racon.git racon \
    && cd racon \
    && mkdir build \
    && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release .. \
    && make \
    && make install

# Install Unicycler 0.4.8
RUN pip3 install git+https://github.com/rrwick/Unicycler.git@v0.5.0

# Install fastq_pair
RUN curl -s -L https://github.com/linsalrob/fastq-pair/archive/v1.0.tar.gz | tar xzf - \
    && cd fastq-pair-1.0 \
    && mkdir build \
    && cd build \
    && cmake ../ \
    && make \
    && make install

# Install SPAdes for unicycler
RUN wget http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Linux.tar.gz \
    && tar -xzf SPAdes-3.15.5-Linux.tar.gz \
    && cd SPAdes-3.15.5-Linux \
    && cp bin/* /usr/local/bin/ \
    && cp -r share/* /usr/local/share/

# Install seqkit v0.12.0 in /usr/local/bin
RUN curl -s -L \
    https://github.com/shenwei356/seqkit/releases/download/v0.12.0/seqkit_linux_amd64.tar.gz \
    | tar -xzf - -C /usr/local/bin

# Install Blast+ v2.10.0 in /usr/local
RUN curl -s -L \
    https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz \
    | tar -xzf - -C /usr/local \
    && ln -s /usr/local/ncbi-blast-2.13.0+/bin/* /usr/local/bin

# R
RUN wget -c https://cloud.r-project.org/src/base/R-4/R-4.2.3.tar.gz \
        && tar -zxf R-4.2.3.tar.gz \
        && cd R-4.2.3 \
        && ./configure \
        && make -j9 \
        && make install


# ----------------------------------------------------------------------
#
# Now construct the final docker image without all of the development
# gunk so that it is a leaner docker image.
#
# ----------------------------------------------------------------------
FROM ubuntu:focal

# Install locales package

# Set locale settings
ENV LANGUAGE=en_US.en
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
RUN apt-get update && apt-get install -y locales && \
	sed -i -e "s/# $LANG.*/$LANG UTF-8/" /etc/locale.gen && \
	dpkg-reconfigure --frontend=noninteractive locales && \
	update-locale LANG=$LANG

# sorter needs pandas and biopython
RUN apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -y --no-install-recommends install \
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
      bedops \
      python3 \
      python3-biopython \
      python3-idna \
      python3-pandas \
      python3-pymummer \
      python3-pip \
      bedtools \
      zlib1g \
      && \
    apt-get clean

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

# downgrade numpy for deprecated np.bool
RUN pip3 install numpy==1.23.1

COPY --from=build /usr/local /usr/local

# Install AliTV

# AliTV needs libyaml-perl libhash-merge-perl bioperl perl cpanm lastz
RUN apt-get -y update && \
  DEBIAN_FRONTEND=noninteractive \
  apt-get -y --no-install-recommends install \
    libyaml-perl \
    libhash-merge-perl \
    bioperl \
    perl \
    wget \
    git \	
    cpanminus && \
  apt-get clean

#RUN cpan FindBin::Real
#RUN cpan Log::Log4perl
#RUN cpanm JSON
#RUN cpanm Bio::Perl Bio::FeatureIO

# Download, compile and install lastz
#RUN wget https://github.com/lastz/lastz/archive/1.04.22.tar.gz && \
#	tar -xf 1.04.22.tar.gz && \
#	cd lastz-1.04.22 && \
#	make && \
#	make install
	# remove -Werror from Makefile to fix compile errors
#	sed -i 's/-Werror //' src/Makefile && \
#	make && \
#	install -m 0755 src/lastz /usr/local/bin/ && \
#	install -m 0755 src/lastz_D /usr/local/bin/ && \
#	cd .. && rm -rf lastz-*

#WORKDIR /app
#COPY . /app
#ENV PERL5LIB="/app/lib:${PERL5LIB}"

# AliTV v1.0.6 install
#RUN git clone https://github.com/AliTVTeam/AlitTV-perl-interface && \
#	cd AliTV-perl-interface \
#	cpanm --installdeps .

#RUN chmod 755 AliTV-perl-interface-1.0.6/bin/alitv.pl

# Install chromoMap
RUN R -e "install.packages('chromoMap', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('htmltools', repos = 'http://cran.us.r-project.org')"

# pilon
ADD https://github.com/broadinstitute/pilon/releases/download/v1.23/pilon-1.23.jar /usr/local/bin/
RUN chmod 755 /usr/local/bin/pilon-1.23.jar
ADD pilon /usr/local/bin/

# Add the entrypoint script
COPY PRYMETIME /usr/local/bin/prymetime

# When the docker image is launched the resulting container runs the entrypoint script
ENTRYPOINT ["/usr/local/bin/prymetime/PRYMETIME.sh"]
