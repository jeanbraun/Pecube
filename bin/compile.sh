#!/bin/sh
cd bin
make clean
cd ..

cd src
make clean
make all
make MPI
cd ..

