FROM debian:testing
MAINTAINER ap13@sanger.ac.uk

RUN apt-get update -qq && apt-get install -y kmc git python3 python3-setuptools python3-biopython python3-pip spades parallel

RUN pip3 install git+git://github.com/sanger-pathogens/plasmidtron.git

WORKDIR /data
