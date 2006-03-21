#!/bin/sh

tmpdir="tmp.$$"
mkdir $tmpdir
cp ZMAT $tmpdir/ZMAT.sym
test -f ECPDATA && cp ECPDATA $tmpdir
cd $tmpdir
sed 's/SYMMETRY=ON/SYMMETRY=OFF/' ZMAT.sym >ZMAT
cp ../GENBAS .
jaces2 --xjoda
cp MOL ../mol
cd ..
rm -rf $tmpdir
