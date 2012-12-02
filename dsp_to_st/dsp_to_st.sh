#!/bin/bash
#
gfortran -c -g dsp_to_st.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dsp_to_st.f90"
  exit
fi
rm compiler.txt
#
gfortran dsp_to_st.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dsp_to_st.o"
  exit
fi
rm dsp_to_st.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/dsp_to_st
#
echo "Executable program installed as ~/bin/$ARCH/dsp_to_st"
