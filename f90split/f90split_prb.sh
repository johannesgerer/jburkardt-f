#!/bin/bash
#
./f90split f90split_prb.f90
#
#  Discard the files we made.
#
rm alpha.f90
rm beta.f90
rm blockdata.f90
rm delta.f90
rm enid.f90
rm epsilon.f90
rm eta.f90
rm gamma.f90
rm iota.f90
rm theta.f90
rm zeta.f90
#
echo "The f90split test problem has been executed."
