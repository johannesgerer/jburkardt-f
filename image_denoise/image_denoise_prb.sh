#!/bin/bash
#
gfortran -c -g image_denoise_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling image_denoise_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran image_denoise_prb.o -L$HOME/lib/$ARCH -limage_denoise
if [ $? -ne 0 ]; then
  echo "Errors linking and loading image_denoise_prb.o"
  exit
fi
rm image_denoise_prb.o
#
mv a.out image_denoise_prb
./image_denoise_prb > image_denoise_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running image_denoise_prb"
  exit
fi
rm image_denoise_prb
#
echo "Program output written to image_denoise_prb_output.txt"
