#/bin/bash
#
cp ../../datasets/uniform/uniform_02_00010.txt .
cp ../voronoi_weight/uniform_02_00010_weight_*.txt .
#
integral_test uniform_02_00010.txt                              > uniform_02_00010_output.txt
integral_test uniform_02_00010.txt uniform_02_00010_weight_a.txt > uniform_02_00010_a_output.txt
integral_test uniform_02_00010.txt uniform_02_00010_weight_b.txt > uniform_02_00010_b_output.txt
integral_test uniform_02_00010.txt uniform_02_00010_weight_c.txt > uniform_02_00010_c_output.txt
integral_test uniform_02_00010.txt uniform_02_00010_weight_d.txt > uniform_02_00010_d_output.txt
integral_test uniform_02_00010.txt uniform_02_00010_weight_e.txt > uniform_02_00010_e_output.txt
integral_test uniform_02_00010.txt uniform_02_00010_weight_f.txt > uniform_02_00010_f_output.txt
integral_test uniform_02_00010.txt uniform_02_00010_weight_g.txt > uniform_02_00010_g_output.txt
#
rm uniform_02_00010.txt
rm uniform_02_00010_weight_*.txt

