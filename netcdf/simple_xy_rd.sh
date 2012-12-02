#!/bin/bash
#
cp ~/include/netcdf.mod .
gfortran -c -g simple_xy_rd.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simple_xy_rd.f90"
  exit
fi
rm compiler.txt
rm netcdf.mod
#
gfortran simple_xy_rd.o -L$HOME/lib/$ARCH -lnetcdf
if [ $? -ne 0 ]; then
  echo "Errors linking and loading simple_xy_rd.o"
  exit
fi
rm simple_xy_rd.o
#
mv a.out simple_xy_rd
./simple_xy_rd > simple_xy_rd_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running simple_xy_rd"
  exit
fi
rm simple_xy_rd
#
echo "Test results written to simple_xy_rd_output.txt."
