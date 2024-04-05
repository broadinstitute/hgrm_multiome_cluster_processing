FROM mambaorg/micromamba:1.5.7

ENV PATH=$PATH:/app

RUN micromamba create -y -n tools2 snapatac2 "scanpy==1.10.0" gcsfs fastparquet -c bioconda -c conda-forge

WORKDIR /app