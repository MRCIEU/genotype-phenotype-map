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
COPY common_cis_trans_variants /home/common_cis_trans_variants/

#delete me after
#root@53d5340a7d9c:/home# history
 #    1  wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
 #    2  ls
 #    3  wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg
 #    4  less /etc/apt/sources.list.d/oneAPI.list
 #    5  cat /etc/apt/sources.list.d/oneAPI.list
 #    6  wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
 #    7  wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg -o --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
 #    8  ls
 #    9  echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
 #   10  echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
 #   11  sudo apt update
 #   12  apt update
 #   13  sudo apt install intel-oneapi-mkl
 #   14  apt install intel-oneapi-mkl
 #   15  wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg -o --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
 #   16  wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg -o --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
 #   17  less /usr/share/keyrings/oneapi-archive-keyring.gpg
 #   18  cat /usr/share/keyrings/oneapi-archive-keyring.gpg
 #   19  echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
 #   20  apt update
 #   21  apt install intel-oneapi-mkl
 #   22  wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg -o --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg
 #   23  wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg
 #   24  cat /usr/share/keyrings/oneapi-archive-keyring.gpg
 #   25  head /usr/share/keyrings/oneapi-archive-keyring.gpg
 #   26  echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
 #   27  sudo apt update
 #   28  apt update
 #   29  sudo apt install intel-oneapi-mkl
 #   30  apt install intel-oneapi-mkl
 #   31  apt install intel-oneapi-mkl-devel
 #   32  apt install intel-oneapi-mkl-devel
 #   33  R
 #   34  cd /usr/lib/
 #   35  ls libmkl_rt.so
 #   36  find . | grep mkl
 #   37  cd ../local/lib/
 #   38  ls
 #   39  cd /lib
 #   40  ls
 #   41  find . | grep mkl
 #   42  find . | grep lib
 #   43  apt-file search libmkl_rt.so
 #   44  dpkg -L intel-oneapi-mkl-devel
 #   45  cd /home/
 #   46  dpkg -L intel-oneapi-mkl-devel
 #   47  dpkg -contents intel-oneapi-mkl-devel
 #   48  dpkg --contents intel-oneapi-mkl-devel
 #   49  dpkg -S libmkl_rt.so
 #   50  LD_PRELOAD=/opt/intel/oneapi/mkl/2024.1/lib/libmkl_rt.so
 #   51  LD_PRELOAD=/opt/intel/oneapi/mkl/2024.1/lib/libmkl_rt.so R
 #   52  history