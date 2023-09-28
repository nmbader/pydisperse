FROM ubuntu:18.04 as builder
MAINTAINER nmbader@sep.stanford.edu
RUN apt-get -y update
RUN apt-get -y install build-essential
RUN apt-get -y install gcc g++ gfortran vim lsof
RUN apt-get -y install pkg-config

RUN apt-get update
RUN apt-get -y  install python3-pip
RUN python3 -m pip install --no-cache-dir --upgrade pip

RUN python3 -m pip install --no-cache-dir numpy &&\
    python3 -m pip install --no-cache-dir jupyter &&\
    python3 -m pip install --no-cache-dir scipy &&\
    python3 -m pip install --no-cache-dir pandas &&\
    python3 -m pip install --no-cache-dir wheel &&\
    python3 -m pip install --no-cache-dir scikit-build &&\
    python3 -m pip install --no-cache-dir matplotlib &&\
    python3 -m pip install --no-cache-dir pybind11

RUN mkdir -p /home
WORKDIR /home

ADD cpp /home/cpp
ADD notebooks /home/notebooks
ADD python /home/python

RUN cd /home &&\
    mkdir -p local &&\
    c++ -std=c++11 -ffast-math -shared -Ofast -fPIC -fopenmp -Icpp -c cpp/functions.cpp -o local/libdispersion.so &&\
    c++ -std=c++11 -ffast-math -shared -Ofast -fPIC -fopenmp `python3 -m pybind11 --includes` -L./local -ldispersion cpp/pybind.cpp  -o local/dispersion`python3-config --extension-suffix`

RUN apt-get -y clean

ENV HOME=/home 
ENV PATH="/home/python:${PATH}"
ENV PYTHONPATH="/home/python:/home/local:${PYTHONPATH}"
RUN echo 'alias python=python3' >> ~/.bashrc