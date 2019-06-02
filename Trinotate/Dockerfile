FROM ubuntu:18.04
ENV PATH=$PATH:/opt/Trinotate-v3.1.1/:/opt/TransDecoder-v5.5.0/ \
    TRINOTATE_HOME=/opt/Trinotate-v3.1.1
# ENV ENV="/etc/profile"
# Installing wget
RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*
# Installing some Perl packages
RUN apt-get update \
  && apt-get install -y libfindbin-libs-perl liburi-perl \
  && rm -rf /var/lib/apt/lists/*
# Installing Blast, SQLite and HMMER
RUN apt-get update \
  && apt-get install -y ncbi-blast+ sqlite3 hmmer\
  && rm -rf /var/lib/apt/lists/*
# Adding packages to /opt
WORKDIR /opt/
# Download and extract Trinotate
RUN wget https://github.com/Trinotate/Trinotate/archive/Trinotate-v3.1.1.tar.gz
RUN mkdir Trinotate-v3.1.1
RUN tar -xzf Trinotate-v3.1.1.tar.gz -C Trinotate-v3.1.1 --strip-components=1
# Download and extract TransDecoder
RUN wget https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v5.5.0.tar.gz
RUN mkdir TransDecoder-v5.5.0
RUN tar -xzf TransDecoder-v5.5.0.tar.gz -C TransDecoder-v5.5.0 --strip-components=1
# Prepare databases for Trinotate
WORKDIR /home/
RUN ${TRINOTATE_HOME}/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
RUN makeblastdb -in uniprot_sprot.pep -dbtype prot
RUN gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm