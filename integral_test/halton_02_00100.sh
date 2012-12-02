#/bin/bash
#
cp ../../datasets/halton/halton_02_00100.txt .
cp ../voronoi_weight/halton_02_00100_weight_*.txt .
#
integral_test halton_02_00100.txt                              > halton_02_00100_output.txt
integral_test halton_02_00100.txt halton_02_00100_weight_a.txt > halton_02_00100_a_output.txt
integral_test halton_02_00100.txt halton_02_00100_weight_b.txt > halton_02_00100_b_output.txt
integral_test halton_02_00100.txt halton_02_00100_weight_c.txt > halton_02_00100_c_output.txt
integral_test halton_02_00100.txt halton_02_00100_weight_d.txt > halton_02_00100_d_output.txt
integral_test halton_02_00100.txt halton_02_00100_weight_e.txt > halton_02_00100_e_output.txt
integral_test halton_02_00100.txt halton_02_00100_weight_f.txt > halton_02_00100_f_output.txt
integral_test halton_02_00100.txt halton_02_00100_weight_g.txt > halton_02_00100_g_output.txt
#
rm halton_02_00100.txt
rm halton_02_00100_weight_*.txt

