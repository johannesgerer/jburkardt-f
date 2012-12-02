#/bin/bash
#
cp ../../datasets/halton/halton_02_01000.txt .
cp ../voronoi_weight/halton_02_01000_weight_*.txt .
#
integral_test halton_02_01000.txt                              > halton_02_01000_output.txt
integral_test halton_02_01000.txt halton_02_01000_weight_a.txt > halton_02_01000_a_output.txt
integral_test halton_02_01000.txt halton_02_01000_weight_b.txt > halton_02_01000_b_output.txt
integral_test halton_02_01000.txt halton_02_01000_weight_c.txt > halton_02_01000_c_output.txt
integral_test halton_02_01000.txt halton_02_01000_weight_d.txt > halton_02_01000_d_output.txt
integral_test halton_02_01000.txt halton_02_01000_weight_e.txt > halton_02_01000_e_output.txt
integral_test halton_02_01000.txt halton_02_01000_weight_f.txt > halton_02_01000_f_output.txt
integral_test halton_02_01000.txt halton_02_01000_weight_g.txt > halton_02_01000_g_output.txt
#
rm halton_02_01000.txt
rm halton_02_01000_weight_*.txt

