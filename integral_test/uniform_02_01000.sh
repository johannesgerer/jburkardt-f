#/bin/bash
#
cp ../../datasets/uniform/uniform_02_01000.txt .
cp ../voronoi_weight/uniform_02_01000_weight_*.txt .
#
integral_test uniform_02_01000.txt                              > uniform_02_01000_output.txt
integral_test uniform_02_01000.txt uniform_02_01000_weight_a.txt > uniform_02_01000_a_output.txt
integral_test uniform_02_01000.txt uniform_02_01000_weight_b.txt > uniform_02_01000_b_output.txt
integral_test uniform_02_01000.txt uniform_02_01000_weight_c.txt > uniform_02_01000_c_output.txt
integral_test uniform_02_01000.txt uniform_02_01000_weight_d.txt > uniform_02_01000_d_output.txt
integral_test uniform_02_01000.txt uniform_02_01000_weight_e.txt > uniform_02_01000_e_output.txt
integral_test uniform_02_01000.txt uniform_02_01000_weight_f.txt > uniform_02_01000_f_output.txt
integral_test uniform_02_01000.txt uniform_02_01000_weight_g.txt > uniform_02_01000_g_output.txt
#
rm uniform_02_01000.txt
rm uniform_02_01000_weight_*.txt

