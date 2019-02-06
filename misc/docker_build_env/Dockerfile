FROM ubuntu:bionic
MAINTAINER SURF Team "davide.vanzo@eawag.ch"

# get the tools we need
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils \
gfortran nano less tree python3 python3-pip doxygen git graphviz locales

RUN pip3 install setuptools wheel
RUN pip3 install FoBiS.py ford pygooglechart

# Set the locale
RUN locale-gen en_US.UTF-8
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'

# add local user
ARG developer_user_id
RUN useradd -u $developer_user_id -m developer

# root dir
RUN mkdir /home/Simstrat
WORKDIR /home/Simstrat

# calls that are needed to build and start the container with the build environment:
# docker build --build-arg developer_user_id=$usr_id -t simstrat:ubuntu "$(pwd)/$(dirname $0)"
# docker create --name simstrat -it -v "$(pwd)/$(dirname $0)/../../":/home/Simstrat simstrat:ubuntu
# docker start simstrat
# docker exec -it -u developer simstrat bash
# docker exec -it -u root simstrat bash