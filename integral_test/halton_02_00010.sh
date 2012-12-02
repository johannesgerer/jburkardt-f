#/bin/bash
#
cp ../../datasets/halton/halton_02_00010.txt .
cp ../voronoi_weight/halton_02_00010_weight_*.txt .
#
integral_test halton_02_00010.txt                              > halton_02_00010_output.txt
integral_test halton_02_00010.txt halton_02_00010_weight_a.txt > halton_02_00010_a_output.txt
integral_test halton_02_00010.txt halton_02_00010_weight_b.txt > halton_02_00010_b_output.txt
integral_test halton_02_00010.txt halton_02_00010_weight_c.txt > halton_02_00010_c_output.txt
integral_test halton_02_00010.txt halton_02_00010_weight_d.txt > halton_02_00010_d_output.txt
integral_test halton_02_00010.txt halton_02_00010_weight_e.txt > halton_02_00010_e_output.txt
integral_test halton_02_00010.txt halton_02_00010_weight_f.txt > halton_02_00010_f_output.txt
integral_test halton_02_00010.txt halton_02_00010_weight_g.txt > halton_02_00010_g_output.txt
#
rm halton_02_00010.txt
rm halton_02_00010_weight_*.txt

