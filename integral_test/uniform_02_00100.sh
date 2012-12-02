#/bin/bash
#
cp ../../datasets/uniform/uniform_02_00100.txt .
cp ../voronoi_weight/uniform_02_00100_weight_*.txt .
#
integral_test uniform_02_00100.txt                              > uniform_02_00100_output.txt
integral_test uniform_02_00100.txt uniform_02_00100_weight_a.txt > uniform_02_00100_a_output.txt
integral_test uniform_02_00100.txt uniform_02_00100_weight_b.txt > uniform_02_00100_b_output.txt
integral_test uniform_02_00100.txt uniform_02_00100_weight_c.txt > uniform_02_00100_c_output.txt
integral_test uniform_02_00100.txt uniform_02_00100_weight_d.txt > uniform_02_00100_d_output.txt
integral_test uniform_02_00100.txt uniform_02_00100_weight_e.txt > uniform_02_00100_e_output.txt
integral_test uniform_02_00100.txt uniform_02_00100_weight_f.txt > uniform_02_00100_f_output.txt
integral_test uniform_02_00100.txt uniform_02_00100_weight_g.txt > uniform_02_00100_g_output.txt
#
rm uniform_02_00100.txt
rm uniform_02_00100_weight_*.txt

