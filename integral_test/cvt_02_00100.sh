#/bin/bash
#
cp ../../datasets/cvt/cvt_02_00100.txt .
cp ../voronoi_weight/cvt_02_00100_weight_*.txt .
#
integral_test cvt_02_00100.txt                              > cvt_02_00100_output.txt
integral_test cvt_02_00100.txt cvt_02_00100_weight_a.txt > cvt_02_00100_a_output.txt
integral_test cvt_02_00100.txt cvt_02_00100_weight_b.txt > cvt_02_00100_b_output.txt
integral_test cvt_02_00100.txt cvt_02_00100_weight_c.txt > cvt_02_00100_c_output.txt
integral_test cvt_02_00100.txt cvt_02_00100_weight_d.txt > cvt_02_00100_d_output.txt
integral_test cvt_02_00100.txt cvt_02_00100_weight_e.txt > cvt_02_00100_e_output.txt
integral_test cvt_02_00100.txt cvt_02_00100_weight_f.txt > cvt_02_00100_f_output.txt
integral_test cvt_02_00100.txt cvt_02_00100_weight_g.txt > cvt_02_00100_g_output.txt
#
rm cvt_02_00100.txt
rm cvt_02_00100_weight_*.txt

