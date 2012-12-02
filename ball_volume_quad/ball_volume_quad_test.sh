#!/bin/bash
#
echo "Estimate volume of 6D sphere using sparse grid rules."
#
cd ../../datasets/sparse_grid_cc
#
~/bin/$ARCH/ball_volume_quad cc_d6_level0 >  ../../f_src/ball_volume_quad/ball_volume_quad_output.txt
~/bin/$ARCH/ball_volume_quad cc_d6_level1 >> ../../f_src/ball_volume_quad/ball_volume_quad_output.txt
~/bin/$ARCH/ball_volume_quad cc_d6_level2 >> ../../f_src/ball_volume_quad/ball_volume_quad_output.txt
~/bin/$ARCH/ball_volume_quad cc_d6_level3 >> ../../f_src/ball_volume_quad/ball_volume_quad_output.txt
~/bin/$ARCH/ball_volume_quad cc_d6_level4 >> ../../f_src/ball_volume_quad/ball_volume_quad_output.txt
~/bin/$ARCH/ball_volume_quad cc_d6_level5 >> ../../f_src/ball_volume_quad/ball_volume_quad_output.txt
#
echo "Estimate volume of 6D sphere using equivalent Monte Carlo rules."
#
cd ../../datasets/quadrature_rules_uniform
#
~/bin/$ARCH/ball_volume_quad uniform_d6_00001 >> ../../f_src/ball_volume_quad/ball_volume_quad_output.txt
~/bin/$ARCH/ball_volume_quad uniform_d6_00013 >> ../../f_src/ball_volume_quad/ball_volume_quad_output.txt
~/bin/$ARCH/ball_volume_quad uniform_d6_00085 >> ../../f_src/ball_volume_quad/ball_volume_quad_output.txt
~/bin/$ARCH/ball_volume_quad uniform_d6_00389 >> ../../f_src/ball_volume_quad/ball_volume_quad_output.txt
~/bin/$ARCH/ball_volume_quad uniform_d6_01457 >> ../../f_src/ball_volume_quad/ball_volume_quad_output.txt
~/bin/$ARCH/ball_volume_quad uniform_d6_04865 >> ../../f_src/ball_volume_quad/ball_volume_quad_output.txt
#
echo " "
echo "(Moral: do not use interpolatory rules on discontinuous data.)"
echo " "
echo "The output was written to ball_volume_quad_test_output.txt"


