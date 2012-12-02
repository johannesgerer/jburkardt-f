#!/bin/bash
#
gfortran -c -g sint_test.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sint_test.f90"
  exit
fi
rm compiler.txt
#
gfortran sint_test.o -L$HOME/lib/$ARCH -lfftpack5
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sint_test.o"
  exit
fi
rm sint_test.o
#
mv a.out sint_test
./sint_test > sint_test_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sint_test"
  exit
fi
rm sint_test
#
echo "Test program output written to sint_test_output.txt."
