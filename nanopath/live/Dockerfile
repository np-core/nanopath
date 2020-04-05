FROM continuumio/miniconda3

LABEL name="nanopath-live"
LABEL version="0.1.0"
LABEL author="esteinig"

RUN apt-get update && apt-get install curl build-essential -y

ENV CONDA_DIR=/opt/conda

RUN conda install -c conda-forge -c bioconda -c esteinig --yes \
    nodejs flask flask-socketio flask-cors\
    && conda clean -a \
    && find $CONDA_DIR -follow -type f -name '*.a' -delete \
    && find $CONDA_DIR -follow -type f -name '*.pyc' -delete

COPY . /live

WORKDIR /live/client
RUN npm install

CMD ["npm", "run", "serve"]


