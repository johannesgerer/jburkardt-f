#!/bin/bash
#
mkdir temp
cd temp
rm *
#
#  Need these first...
#
gfortran -c -g ../buckley.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling buckley.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../internal_types.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling internal_types.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../ds_routines.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ds_routines.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../error_handling.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling error_handling.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../rule_1.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rule_1.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../rule_c2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rule_c2.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../rule_c3.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rule_c3.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../rule_cn.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rule_cn.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../rule_t2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rule_t2.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../rule_t3.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rule_t3.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../rule_tn.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rule_tn.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../rule_general.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rule_general.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../volume.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling volume.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../divide.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling divide.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../region_processor.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling region_processor.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../global_all.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling global_all.f90."
  exit
fi
rm compiler.txt
#
#  Now the rest.
#
gfortran -c -g ../check.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling check.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g ../cui.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cui.f90."
  exit
fi
rm compiler.txt
#
echo "Create the archive."
ar qc libcubpack.a *.o *.mod
rm *.o
rm *.mod
#
echo "Store the archive."
mv libcubpack.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "A new version of cubpack has been created."
