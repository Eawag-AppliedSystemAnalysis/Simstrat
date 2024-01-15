FROM ubuntu
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils \
gfortran gfortran-mingw-w64 make nano less tree python3 python3-pip doxygen git graphviz locales
RUN pip3 install pip --upgrade
RUN pip3 install setuptools
RUN pip3 install wheel
RUN pip3 install future
RUN pip3 install FoBiS.py ford pygooglechart
RUN ln -sf /usr/bin/x86_64-linux-gnu-gfortran /usr/bin/gfortran
RUN locale-gen en_US.UTF-8
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'
RUN mkdir -p /simstrat/run
COPY build /simstrat/build
COPY lib /simstrat/lib
COPY src /simstrat/src
ENV F90=gfortran
WORKDIR /simstrat/lib/libaed2
RUN make
WORKDIR /simstrat/build
RUN FoBiS.py build
COPY entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]