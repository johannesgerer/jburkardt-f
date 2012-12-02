#!/bin/bash
#
gfortran fft_serial.f90
mv a.out ~/bin/OSX/fft_serial
echo "Executable installed as ~/bin/OSX/fft_serial"
