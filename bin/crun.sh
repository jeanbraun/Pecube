#!/bin/sh

cd src
make all
cd ..
bin/Test $1
bin/Pecube $1
bin/Vtk $1
