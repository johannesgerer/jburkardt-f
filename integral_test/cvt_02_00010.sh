#/bin/bash
#
cp ../../datasets/cvt/cvt_02_00010.txt .
cp ../voronoi_weight/cvt_02_00010_weight_*.txt .
#
integral_test cvt_02_00010.txt                           > cvt_02_00010_output.txt
integral_test cvt_02_00010.txt cvt_02_00010_weight_a.txt > cvt_02_00010_a_output.txt
integral_test cvt_02_00010.txt cvt_02_00010_weight_b.txt > cvt_02_00010_b_output.txt
integral_test cvt_02_00010.txt cvt_02_00010_weight_c.txt > cvt_02_00010_c_output.txt
integral_test cvt_02_00010.txt cvt_02_00010_weight_d.txt > cvt_02_00010_d_output.txt
integral_test cvt_02_00010.txt cvt_02_00010_weight_e.txt > cvt_02_00010_e_output.txt
integral_test cvt_02_00010.txt cvt_02_00010_weight_f.txt > cvt_02_00010_f_output.txt
integral_test cvt_02_00010.txt cvt_02_00010_weight_g.txt > cvt_02_00010_g_output.txt
#
rm cvt_02_00010.txt
rm cvt_02_00010_weight_*.txt
