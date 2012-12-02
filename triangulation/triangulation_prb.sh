#!/bin/bash
#
gfortran -c -g triangulation_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran triangulation_prb.o -L$HOME/lib/$ARCH -ltriangulation
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_prb.o"
  exit
fi
rm triangulation_prb.o
#
mv a.out triangulation_prb
./triangulation_prb > triangulation_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running triangulation_prb"
  exit
fi
rm triangulation_prb
#
echo "Test program output written to triangulation_prb_output.txt."
#
if [ -e ns_triangulation.eps ]; then
  convert ns_triangulation.eps ns_triangulation.png
  rm ns_triangulation.eps
fi
#
if [ -e triangulation_order3_plot.eps ]; then
  convert triangulation_order3_plot.eps triangulation_order3_plot.png
  rm triangulation_order3_plot.eps
fi
#
if [ -e triangulation_order3_plot2.eps ]; then
  convert triangulation_order3_plot2.eps triangulation_order3_plot2.png
  rm triangulation_order3_plot2.eps
fi
#
if [ -e triangulation_order6_plot.eps ]; then
  convert triangulation_order6_plot.eps triangulation_order6_plot.png
  rm triangulation_order6_plot.eps
fi
#
echo "EPS figures converted to PNG format."
