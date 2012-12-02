#!/bin/bash
#
gfortran -c -g calendar_nyt_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling calendar_nyt_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran calendar_nyt_prb.o -L$HOME/lib/$ARCH -lcalendar_nyt
if [ $? -ne 0 ]; then
  echo "Errors linking and loading calendar_nyt_prb.o"
  exit
fi
rm calendar_nyt_prb.o
#
mv a.out calendar_nyt_prb
./calendar_nyt_prb > calendar_nyt_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running calendar_nyt_prb"
  exit
fi
rm calendar_nyt_prb
#
echo "Test program output written to calendar_nyt_prb_output.txt."
