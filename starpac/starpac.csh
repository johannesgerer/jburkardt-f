#!/bin/csh
#
mkdir temp
cd temp
rm *
f90split ../starpac.f90
#
foreach FILE (`ls -1 *.f90`)
  F90 -c -g $FILE >& compiler.txt
  if ( $status != 0 ) then
    echo "Errors compiling " $FILE
    exit
  endif
  rm compiler.txt
end
rm *.f90
#
ar qc libstarpac.a *.o
rm *.o
#
mv libstarpac.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libstarpac.a"
