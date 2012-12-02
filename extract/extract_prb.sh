#!/bin/bash
#
extract beta extract_prb.f90
#
cat beta.f90
rm beta.f90
#
extract theta extract_prb.f90
cat theta.f90
rm theta.f90
#
echo "Normal end of execution of EXTRACT tests."
