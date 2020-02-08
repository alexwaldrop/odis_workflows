# Base Image
FROM python:3.7.6-slim-stretch

# Metadata
LABEL base.image="python:3.7.6-slim-stretch"
LABEL version="1"
LABEL software="merge_ref_contigs"
LABEL software.version="latest"
LABEL description="Bioinformatics utility for concatenating large numbers of reference contigs into a smaller number of artificial contigs"

# Maintainer
MAINTAINER Alex Waldrop <awaldrop@rti.org>

# update the OS related packages
RUN apt-get update

# install required dependencies for QCParser
ADD requirements.txt .
RUN pip install -r requirements.txt && rm requirements.txt

# Add files to folder
ADD merge_ref_contigs.py /opt/
RUN chmod 755 /opt/merge_ref_contigs.py

ENV PATH /opt/:$PATH

CMD ["merge_ref_contigs.py", "--help"]