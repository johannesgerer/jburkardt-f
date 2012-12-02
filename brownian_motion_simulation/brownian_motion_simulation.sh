#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../brownian_motion_simulation.f90
#
for FILE in `ls -1 *.f90`;
do
  gfortran -c -g $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
rm *.f90
#
ar qc libbrownian_motion_simulation.a *.o
rm *.o
#
mv libbrownian_motion_simulation.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libbrownian_motion_simulation.a"
