#!/bin/bash
set -e
set -x

start_dir=$(pwd)

KMC_VERSION=${KMC_VERSION:-"3.0.0"}
SPADES_VERSION=3.10.1

KMC3_DOWNLOAD_URL="https://github.com/refresh-bio/KMC/releases/download/v${KMC_VERSION}/KMC3.linux.tar.gz"
KMC2_DOWNLOAD_URL_BASE="http://sun.aei.polsl.pl/REFRESH/kmc/downloads/${KMC_VERSION}/linux/"

SPADES_URL="http://cab.spbu.ru/files/release${SPADES_VERSION}/SPAdes-${SPADES_VERSION}-Linux.tar.gz"

PARALLEL_VERSION="20170122"
PARALLEL_DOWNLOAD_FILENAME="parallel-${PARALLEL_VERSION}.tar.bz2" 
PARALLEL_URL="http://ftp.gnu.org/gnu/parallel/${PARALLEL_DOWNLOAD_FILENAME}"

# Make an install location
if [ ! -d 'build' ]; then
  mkdir build
fi
cd build
build_dir=$(pwd)

# DOWNLOAD ALL THE THINGS
download () {
  url=$1
  download_location=$2

  if [ -e $download_location ]; then
    echo "Skipping download of $url, $download_location already exists"
  else
    echo "Downloading $url to $download_location"
    wget $url -O $download_location
  fi
}

# --------------- KMC ------------------
cd $build_dir

if [ "${KMC_VERSION}" == "2.3.0" ]; then 
  download "${KMC2_DOWNLOAD_URL_BASE}kmc" "kmc"
  download "${KMC2_DOWNLOAD_URL_BASE}kmc_tools" "kmc_tools"
  download "${KMC2_DOWNLOAD_URL_BASE}kmc_dump" "kmc_dump"
else
  download "${KMC3_DOWNLOAD_URL}" "KMC3.linux.tar.gz"
  tar xzf KMC3.linux.tar.gz
fi
chmod +x kmc
chmod +x kmc_tools
chmod +x kmc_dump

# --------------- SPAdes ------------------
cd $build_dir
download $SPADES_URL "SPAdes-${SPADES_VERSION}-Linux.tar.gz"
tar zxf "SPAdes-${SPADES_VERSION}-Linux.tar.gz"
spades_dir="$build_dir/SPAdes-${SPADES_VERSION}-Linux"
cd $spades_dir

# --------------- parallel ------------------
cd $build_dir
PARALLEL_BUILD_DIR="${build_dir}/parallel-${PARALLEL_VERSION}"
download $PARALLEL_URL $PARALLEL_DOWNLOAD_FILENAME
tar xjvf $PARALLEL_DOWNLOAD_FILENAME
cd $PARALLEL_BUILD_DIR
./configure
make

cd $start_dir

update_path () {
  new_dir=$1
  if [[ ! "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]; then
    export PATH=${new_dir}:${PATH}
  fi
}

update_path ${build_dir}
update_path "${spades_dir}/bin"
update_path "${PARALLEL_BUILD_DIR}/src"

pip install pyfastaq biopython matplotlib

