FROM ubuntu:14.04

RUN apt-get update -y && apt-get install -y build-essential git cmake libgsl0-dev python-dev python-numpy libboost-dev libboost-serialization-dev libboost-iostreams-dev libboost-program-options-dev libboost-filesystem-dev libboost-system-dev swig

RUN git clone https://github.com/jgasthaus/libplump

RUN cd libplump && mkdir build && cd build && cmake .. && make
