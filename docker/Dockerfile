FROM mrcieu/genehackman:1.0.0

RUN wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list && \
    apt update && apt install -y --no-install-recommends libhiredis-dev less libgmp-dev libmpfr-dev jq csvtool libgsl-dev intel-oneapi-mkl default-jre rsync
    # export PERL_MM_USE_DEFAULT=1 && cpan DBI DBD::mysql Archize::Zip
    #uncomment this if you need to start using GPU / cuda again
    # wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-keyring_1.1-1_all.deb && \
    # dpkg -i cuda-keyring_1.1-1_all.deb && \
    # apt-get update && \
    # apt-get -y install cuda-toolkit-12-6

#Installing ssimp, finemap, smr, gcta, DENTIST, picard, and gwcoco here...
RUN cd /home && \
    wget http://christianbenner.com/finemap_v1.4.2_x86_64.tgz && \
    tar -xf finemap_v1.4.2_x86_64.tgz && \
    mv finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 /bin/finemap && \
    wget https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.1-linux-x86_64.zip && \
    unzip smr-1.3.1-linux-x86_64.zip && \
    mv smr-1.3.1-linux-x86_64/smr /bin/smr && \
    # wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip && \
    # unzip gcta-1.94.1-linux-kernel-3-x86_64.zip && \
    # cp gcta-1.94.1-linux-kernel-3-x86_64/gcta64 /bin/gcta64 && \
    wget https://github.com/Kitware/CMake/releases/download/v3.30.0-rc1/cmake-3.30.0-rc1-linux-x86_64.tar.gz && \
    tar xf cmake-3.30.0-rc1-linux-x86_64.tar.gz && \
    git clone https://github.com/jwr-git/pwcoco.git && \
    cd pwcoco && \
    mkdir build && cd build && /home/cmake-3.30.0-rc1-linux-x86_64/bin/cmake .. && make && \
    mv pwcoco /bin/pwcoco && \
    cd /home && \
    wget https://github.com/zkutalik/ssimp_software/raw/master/compiled/ssimp-linux-0.5.6 && \
    chmod 775 ssimp-linux-0.5.6 && \
    mv ssimp-linux-0.5.6 /bin/ssimp && \
    wget -O DENTIST.gz 'https://www.dropbox.com/scl/fi/1yna9rizfy6s0ixek0cca/DENTIST.1.3.0.0.gz?rlkey=0ak0cdsev152jteibceit4mje&dl=0' && \
    gunzip DENTIST.gz && \
    chmod 775 DENTIST && \
    mv DENTIST /bin/DENTIST && \
    wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240820.zip && \
    unzip plink2_linux_x86_64_20240820.zip && \
    mv plink2 /usr/bin/ && \
    wget https://github.com/broadinstitute/picard/releases/download/3.2.0/picard.jar && \
    chmod 775 picard.jar && \
    mv picard.jar /usr/bin && \
    rm -r bcftools htslib || true && \
    git clone --recurse-submodules https://github.com/samtools/htslib.git && \
    git clone https://github.com/samtools/bcftools.git && \
    cd bcftools && \
    make && \
    wget https://software.broadinstitute.org/software/score/score_1.20-20240505.zip && \
    unzip score_1.20-20240505.zip && \
    cd /home

ENV BCFTOOLS_PLUGINS=/home/bcftools
ENV LD_PRELOAD=/opt/intel/oneapi/mkl/latest/lib/libmkl_rt.so

COPY docker/requirements.* docker/
RUN Rscript docker/requirements.R
RUN pip install -r docker/requirements.txt

COPY . /home/pipeline
