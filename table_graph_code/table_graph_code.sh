#!/bin/bash
#
gfortran -c -g table_graph_code.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_graph_code.f90"
  exit
fi
rm compiler.txt
#
gfortran table_graph_code.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_graph_code.o"
  exit
fi
rm table_graph_code.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_graph_code
#
echo "Program installed as ~/bin/$ARCH/table_graph_code"
