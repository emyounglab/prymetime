FROM ubuntu:focal as build

RUN apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -y install \
      cmake \
      curl \
      git \
      mummer \
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


# ----------------------------------------------------------------------
#
# Now construct the final docker image without all of the development
# gunk so that it is a leaner docker image.
#
# ----------------------------------------------------------------------
FROM ubuntu:focal

RUN apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -y --no-install-recommends install \
      python \
      && \
    apt-get clean

COPY --from=build /usr/local /usr/local

# Add the entrypoint script
COPY PRYMETIME /usr/local/bin/prymetime

# When the docker image is launched the resulting container runs the entrypoint script
ENTRYPOINT ["/usr/local/bin/prymetime/PRYMETIME.sh"]
