FROM ubuntu:18.04
ENV PATH=$PATH:/opt/Trinotate-v3.1.1/:/opt/TransDecoder-v5.5.0/:/opt/Blast+-v2.9.0/bin \
    TRINOTATE_HOME=/opt/Trinotate-v3.1.1
# ENV ENV="/etc/profile"
# Installing wget, make, gcc, and other stuff
RUN  apt-get update \
  && apt-get install -y wget make gcc \
  && rm -rf /var/lib/apt/lists/*
# Installing some Perl packages
RUN apt-get update \
  && apt-get install -y libfindbin-libs-perl liburi-perl \
  && rm -rf /var/lib/apt/lists/*
# Installing SQLite
RUN apt-get update \
  && apt-get install -y sqlite3 \
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
# Dowload latest version of BLAST+
RUN wget -O ncbi-blast-v2.9.0.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz
RUN mkdir Blast+-v2.9.0
RUN tar -xzf ncbi-blast-v2.9.0.tar.gz -C Blast+-v2.9.0 --strip-components=1
# Dowload latest version of HMMER
RUN wget http://eddylab.org/software/hmmer/hmmer-3.2.1.tar.gz
RUN tar -xzf hmmer-3.2.1.tar.gz
WORKDIR /opt/hmmer-3.2.1
RUN ./configure \
    && make \
    && make check \
    && make install
# Prepare databases for Trinotate
WORKDIR /home/
RUN ${TRINOTATE_HOME}/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
RUN makeblastdb -in uniprot_sprot.pep -dbtype prot
RUN gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm