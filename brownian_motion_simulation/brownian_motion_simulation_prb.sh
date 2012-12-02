#!/bin/bash
#
gfortran -c -g brownian_motion_simulation_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling brownian_motion_simulation_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran brownian_motion_simulation_prb.o -L$HOME/lib/$ARCH -lbrownian_motion_simulation
if [ $? -ne 0 ]; then
  echo "Errors linking and loading brownian_motion_simulation_prb.o"
  exit
fi
rm brownian_motion_simulation_prb.o
#
mv a.out brownian_motion_simulation_prb
./brownian_motion_simulation_prb > brownian_motion_simulation_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running brownian_motion_simulation_prb"
  exit
fi
rm brownian_motion_simulation_prb
#
echo "Test program output written to brownian_motion_simulation_prb_output.txt."
