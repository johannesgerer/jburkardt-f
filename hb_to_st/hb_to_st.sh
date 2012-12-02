#!/bin/bash
#
gfortran -c -g hb_to_st.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hb_to_st.f90"
  exit
fi
rm compiler.txt
#
gfortran hb_to_st.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hb_to_st.o"
  exit
fi
rm hb_to_st.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/hb_to_st
#
echo "Executable installed as ~/bin/$ARCH/hb_to_st"
