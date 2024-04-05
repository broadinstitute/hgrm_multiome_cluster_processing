FROM us.gcr.io/landerlab-atacseq-200218/gcloud-samtools:0.1

ENV PATH=$PATH:/app
COPY requirements.txt .
RUN pip3 install --break-system-packages -r requirements.txt
RUN rm requirements.txt  # cleanup

RUN curl -LO https://github.com/kaizhang/SnapATAC2/releases/download/nightly/snapatac2-2.6.1.dev0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
RUN pip3 install snapatac2-2.6.1.dev0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
RUN rm snapatac2-2.6.1.dev0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl # cleanup

WORKDIR /app