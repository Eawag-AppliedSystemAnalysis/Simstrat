FROM ubuntu:bionic
LABEL org.opencontainers.imag.authors="SURF Team" email="davide.vanzo@eawag.ch"

# get the tools we need
RUN DEBIAN_FRONTEND=noninteractive apt-get update && apt-get install -y --no-install-recommends \
apt-utils \
gfortran gfortran-mingw-w64 \
make \
nano \
less \
tree \
python3 python3-pip \
doxygen \
git \
graphviz \
locales \
wget \
ca-certificates && \
# Install cmake 4.0.2 manually
wget https://github.com/Kitware/CMake/releases/download/v4.0.2/cmake-4.0.2-linux-x86_64.sh && \
bash cmake-4.0.2-linux-x86_64.sh --skip-license --prefix=/usr/local && \
rm cmake-4.0.2-linux-x86_64.sh && \
rm -rf /var/lib/apt/lists/*

# install some needed python modules
RUN pip3 install pip --upgrade
RUN pip3 install setuptools
RUN pip3 install wheel
RUN pip3 install future
RUN pip3 install FoBiS.py ford pygooglechart

# set alternatives for the compiler
RUN update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/x86_64-linux-gnu-gfortran 100
RUN update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/x86_64-w64-mingw32-gfortran 50

# to change manually the alternatives with (requires root user)
# update-alternatives --config gfortran

# Set the locale
RUN locale-gen en_US.UTF-8
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'

# root dir
RUN mkdir /home/Simstrat
WORKDIR /home/Simstrat

# calls that are needed to build and start the container with the build environment:
# docker build -t simstrat:ubuntu "$(pwd)/$(dirname $0)"
# docker create --name simstrat -it -v "$(pwd)/$(dirname $0)/../../":/home/Simstrat simstrat:ubuntu
# docker start simstrat

# enter the terminal
# docker exec -it simstrat bash

# change the compiler
# docker exec -it simstrat bash -c "./misc/docker_scripts/set_MinGW_compiler.sh"
# docker exec -it simstrat bash -c "./misc/docker_scripts/set_GNU_compiler.sh"