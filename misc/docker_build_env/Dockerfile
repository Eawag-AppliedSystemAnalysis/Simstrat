FROM ubuntu:bionic
MAINTAINER SURF Team "davide.vanzo@eawag.ch"

# get the tools we need
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils \
gfortran nano less tree python python-pip doxygen git graphviz

RUN pip install setuptools wheel
RUN pip install FoBiS.py ford pygooglechart

# root dir
RUN mkdir /home/Simstrat
WORKDIR /home/Simstrat

# calls that are needed to build and start the container with the build environment:
# docker build -t simstrat:ubuntu .
# docker create --name simstrat -it -v <pathToLocalGitRepoDirectory>:/home/Simstrat simstrat:ubuntu
# docker start simstrat
# docker exec -it simstrat bash
