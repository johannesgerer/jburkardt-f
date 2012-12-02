#!/bin/bash
#
gfortran -c -g graph_paper_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling graph_paper_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran graph_paper_prb.o -L$HOME/lib/$ARCH -lgraph_paper
if [ $? -ne 0 ]; then
  echo "Errors linking and loading graph_paper_prb.o"
  exit
fi
rm graph_paper_prb.o
#
mv a.out graph_paper_prb
./graph_paper_prb > graph_paper_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running graph_paper_prb"
  exit
fi
rm graph_paper_prb
#
echo "Test program output written to graph_paper_prb_output.txt."
#
convert hexagonal_1.eps hexagonal_1.png
rm hexagonal_1.eps
#
convert hexagonal_2.eps hexagonal_2.png
rm hexagonal_2.eps
#
convert hexagonal_3.eps hexagonal_3.png
rm hexagonal_3.eps
#
convert hexagonal_4.eps hexagonal_4.png
rm hexagonal_4.eps
#
convert hexagonal_5.eps hexagonal_5.png
rm hexagonal_5.eps
#
convert polar_1.eps polar_1.png
rm polar_1.eps
#
convert sudoku_blank.eps sudoku_blank.png
rm sudoku_blank.eps
#
convert sudoku_filled.eps sudoku_filled.png
rm sudoku_filled.eps
#
convert staggered_2.eps staggered_2.png
rm staggered_2.eps
#
convert triangular_1.eps triangular_1.png
rm triangular_1.eps
#
convert triangular_2.eps triangular_2.png
rm triangular_2.eps
#
convert uniform_1.eps uniform_1.png
rm uniform_1.eps
#
convert uniform_2.eps uniform_2.png
rm uniform_2.eps
#
convert uniform_3.eps uniform_3.png
rm uniform_3.eps
#
convert uniform_4.eps uniform_4.png
rm uniform_4.eps
#
echo "EPS graphics files converted to PNG format."

