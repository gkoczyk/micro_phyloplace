FROM mambaorg/micromamba:1.3.1-bionic
COPY --chown=$MAMBA_USER:$MAMBA_USER data/phyloplace.yml /tmp/phyloplace.yaml
ARG conda_env=phyloplace
USER root
RUN apt-get update && apt-get install libgl1-mesa-glx libxcb-util1 libxcb-render-util0 libxcb-keysyms1  libxcb-icccm4 libxcb-image0 xvfb -y
USER mambauser
RUN micromamba install -y -n base -f /tmp/phyloplace.yaml && \
    micromamba clean --all --yes
## RUN conda env create -f data/phyloplace.yml -n phyloplace
## RUN echo "conda activate phyloplace" >> ~/.bashrc
## ENV PATH /opt/conda/envs/phyloplace/bin:$PATH
## ENV CONDA_DEFAULT_ENV $conda_env
ENV WISECONFIGDIR=/opt/conda/share/wise2/wisecfg/
ENV XDG_RUNTIME_DIR=/tmp/runtime-mambauser
##CMD [ "python", "test.py" ]
ADD scripts/ /scripts/

## RUN chown -R root:root /workdir
## RUN chmod 755 /workdir
ADD data /data/
ADD scripts/runner.sh /usr/local/bin/run
USER root
RUN chmod ugo+rwx /data
RUN chmod 755 /usr/local/bin/run
USER mambauser
VOLUME /current
WORKDIR /current
ENTRYPOINT ["/usr/local/bin/run"]