#!/bin/bash
#
gfortran -c mgs.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling mgs.f90"
  exit
fi
#
echo "The mgs.f90 file was compiled."
