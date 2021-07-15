#FROM ubuntu:21.04
FROM continuumio/miniconda3:latest

ENV DEBIAN_FRONTEND=noninteractive

#RUN apt-get update && \
#    apt-get -y install python3-pip python3-rdkit python3-pandas python3-ipython && \
#    apt-get -y clean && \
#    rm -rf \
#      /var/lib/apt/lists/* \
#      /usr/share/doc \
#      /usr/share/doc-base \
#      /usr/share/man \
#      /usr/share/locale \
#      /usr/share/zoneinfo

#ADD requirements.txt /
#RUN pip install --trusted-host pypi.python.org -r requirements.txt

RUN conda create -c conda-forge -n my-rdkit-env rdkit

RUN conda activate my-rdkit-env ; conda install -n my-rdkit-env ipython
RUN conda activate my-rdkit-env ; conda install -n my-rdkit-env pandas

WORKDIR /tmp
COPY classify_homols.py /tmp/
COPY input /tmp/input

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "my-rdkit-env", "python3", "classify_homols.py"]
