#!/bin/bash
#
gfortran -c -g discrete_pdf_sample.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling discrete_pdf_sample.f90"
  exit
fi
rm compiler.txt
#
gfortran discrete_pdf_sample.o
if [ $? -ne 0 ]; then
  echo "Errors while loading discrete_pdf_sample.o"
  exit
fi
rm discrete_pdf_sample.o
#
mv a.out ~/bin/$ARCH/discrete_pdf_sample
#
echo "Executable installed as ~/bin/$ARCH/discrete_pdf_sample"
