#/bin/bash
#
cp ../../datasets/cvt/cvt_02_01000.txt .
cp ../voronoi_weight/cvt_02_01000_weight_*.txt .
#
integral_test cvt_02_01000.txt                              > cvt_02_01000_output.txt
integral_test cvt_02_01000.txt cvt_02_01000_weight_a.txt > cvt_02_01000_a_output.txt
integral_test cvt_02_01000.txt cvt_02_01000_weight_b.txt > cvt_02_01000_b_output.txt
integral_test cvt_02_01000.txt cvt_02_01000_weight_c.txt > cvt_02_01000_c_output.txt
integral_test cvt_02_01000.txt cvt_02_01000_weight_d.txt > cvt_02_01000_d_output.txt
integral_test cvt_02_01000.txt cvt_02_01000_weight_e.txt > cvt_02_01000_e_output.txt
integral_test cvt_02_01000.txt cvt_02_01000_weight_f.txt > cvt_02_01000_f_output.txt
integral_test cvt_02_01000.txt cvt_02_01000_weight_g.txt > cvt_02_01000_g_output.txt
#
rm cvt_02_01000.txt
rm cvt_02_01000_weight_*.txt

