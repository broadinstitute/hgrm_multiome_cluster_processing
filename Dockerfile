FROM mambaorg/micromamba:1.5.7

ENV PATH=$PATH:/app

RUN micromamba create -y -n tools2 snapatac2 "scanpy==1.10.0" gcsfs fastparquet -c bioconda -c conda-forge
RUN micromamba run -n tools2 python3 -m pip install argparse

USER root
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git

WORKDIR /app

COPY monitor_script.sh .
RUN chmod a+rx /app/monitor_script.sh