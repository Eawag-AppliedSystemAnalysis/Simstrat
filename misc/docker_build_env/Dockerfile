FROM alpine:latest
MAINTAINER SURF Team "davide.vanzo@eawag.ch"

# get the tools we need
RUN apk update && apk add gfortran \
	musl-dev bash python py-pip doxygen git

RUN pip install FoBiS.py pygooglechart

# root dir
RUN mkdir /home/Simstrat
WORKDIR /home/Simstrat

# calls that are needed to build and start the container with the build environment:
# docker build -t simstrat:alpine .
# docker create --name simstrat -it -v <pathToLocalGitRepoDirectory>:/home/Simstrat simstrat:alpine
# docker start simstrat
# docker exec -it simstrat bash
