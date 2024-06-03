FROM mrcieu/genehackman:1.0.0

RUN wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list && \
    apt update && apt install -y --no-install-recommends libgmp-dev libmpfr-dev jq csvtool libgsl-dev intel-oneapi-mkl

RUN wget http://christianbenner.com/finemap_v1.4.2_x86_64.tgz && \
    tar -xf finemap_v1.4.2_x86_64.tgz && \
    mv finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 /bin/finemap

ENV LD_PRELOAD=/opt/intel/oneapi/mkl/2024.1/lib/libmkl_rt.so
COPY docker/requirements.* docker/
RUN Rscript docker/requirements.R
COPY common_cis_variants /home/common_cis_variants/

