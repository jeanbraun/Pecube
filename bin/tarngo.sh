#!/bin/sh

mkdir Pecube

echo 'This is version '$1 > Pecube/Version.txt

mkdir Pecube/src
cp src/*.f90 Pecube/src
cp src/*.f Pecube/src
cp src/*.c Pecube/src
cp src/*.inc Pecube/src
cp src/*.h Pecube/src
cp src/Makefile Pecube/src

cp -r bin Pecube

cp -r docs Pecube

cp -r old_input Pecube

cp -r tools Pecube

mkdir Pecube/EXMP1
mkdir Pecube/EXMP1/input
cp EXMP1/input/Pecube.in Pecube/EXMP1/input/
cp -r EXMP1/data Pecube/EXMP1

mkdir Pecube/EXMP2
mkdir Pecube/EXMP2/input
cp EXMP2/input/Pecube.in Pecube/EXMP2/input/
cp -r EXMP2/data Pecube/EXMP2

mkdir Pecube/EXMP3
mkdir Pecube/EXMP3/input
cp EXMP3/input/Pecube.in Pecube/EXMP3/input/
cp -r EXMP3/data Pecube/EXMP3

mkdir Pecube/EXMP4
mkdir Pecube/EXMP4/input
cp EXMP4/input/Pecube.in Pecube/EXMP4/input/
cp -r EXMP4/data Pecube/EXMP4

mkdir Pecube/EXMP5
mkdir Pecube/EXMP5/input
cp EXMP5/input/Pecube.in Pecube/EXMP5/input/
cp -r EXMP5/data Pecube/EXMP5

mkdir Pecube/EXMP6
mkdir Pecube/EXMP6/input
cp EXMP6/input/Pecube.in Pecube/EXMP6/input/
cp -r EXMP6/data Pecube/EXMP6

mkdir Pecube/EXMP7
mkdir Pecube/EXMP7/input
cp EXMP7/input/Pecube.in Pecube/EXMP7/input/
cp -r EXMP7/data Pecube/EXMP7

mkdir Pecube/EXMP8
mkdir Pecube/EXMP8/input
cp EXMP8/input/Pecube.in Pecube/EXMP8/input/
cp -r EXMP8/data Pecube/EXMP8

mkdir Pecube/EXMP9
mkdir Pecube/EXMP9/input
cp EXMP9/input/Pecube.in Pecube/EXMP9/input/
cp -r EXMP9/data Pecube/EXMP9

mkdir Pecube/EXMPA
mkdir Pecube/EXMPA/input
cp EXMPA/input/Pecube.in Pecube/EXMPA/input/
cp -r EXMPA/data Pecube/EXMPA

mkdir Pecube/EXMPB
mkdir Pecube/EXMPB/input
cp EXMPB/input/Pecube.in Pecube/EXMPB/input/
cp -r EXMPB/data Pecube/EXMPB

mkdir Pecube/topo

mkdir Pecube/tmp

tar -cvzf Pecube_$1.tar.gz Pecube

rm -r Pecube
