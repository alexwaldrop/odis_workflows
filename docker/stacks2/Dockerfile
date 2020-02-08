# Base image
FROM ubuntu:18.04

# Maintainer and author
LABEL maintainer="Alex Waldrop <awaldrop@rti.org>"

#### Basic Ubuntu setup ####

ENV DEBIAN_FRONTEND noninteractive

RUN echo 'deb http://mirror.math.princeton.edu/pub/ubuntu/ bionic main' >> /etc/apt/sources.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        zlib1g-dev \
        gcc-7 \
        wget \
        locales \
    && rm -rf /var/lib/apt/lists/*

# Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN apt-get update \
    && wget http://catchenlab.life.illinois.edu/stacks/source/stacks-2.5.tar.gz \
    && tar xzf stacks-2.5.tar.gz \
    && cd stacks-2.5 \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -rf stacks-2.5

WORKDIR /

CMD ["process_radtags", "-h"]