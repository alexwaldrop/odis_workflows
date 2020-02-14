# Base Image
FROM python:3.7.6-slim-stretch

# Metadata
LABEL base.image="python:3.7.6-slim-stretch"
LABEL version="1"
LABEL software="CAT"
LABEL software.version="latest"
LABEL description="Bioinformatics utility for concatenating large numbers of reference contigs into a smaller number of artificial contigs"

# Maintainer
MAINTAINER Alex Waldrop <awaldrop@rti.org>


ENV DEBIAN_FRONTEND noninteractive

RUN mkdir prodigal

COPY diamond /opt/
COPY Prodigal prodigal
COPY CAT /opt/

RUN ls -l
RUN ls -l /opt/
RUN ls -l /prodigal

# Install prodigal
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libglib2.0 \
        wget \
        make \
        gcc \
        g++ \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
    && cd /prodigal \
    && make install INSTALLDIR=/opt/ \
    && rm -rf prodigal

# Add files to folder
RUN chmod -R 755 /opt/

ENV PATH /opt/:$PATH
ENV PATH /opt/CAT_pack/:$PATH

CMD ["CAT", "--help"]