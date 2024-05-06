FROM mrcieu/genehackman:1.0.0

RUN apt update && apt install -y --no-install-recommends libgmp-dev libmpfr-dev

RUN wget http://christianbenner.com/finemap_v1.4.2_x86_64.tgz && \
    tar -xf finemap_v1.4.2_x86_64.tgz && \
    mv finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 /bin/finemap

COPY docker/requirements.* docker/
RUN Rscript docker/requirements.R
