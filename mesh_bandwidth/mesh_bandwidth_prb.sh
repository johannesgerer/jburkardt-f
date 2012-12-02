#!/bin/bash
#
rm mesh_bandwidth_prb_output.txt
#
cp ../../data/polygonal_surface/sphere_q4_elements.txt sphere_q4_elements.txt
mesh_bandwidth sphere_q4_elements.txt > mesh_bandwidth_prb_output.txt
rm sphere_q4_elements.txt
#
cp ../../data/polygonal_surface/sphere_t3_elements.txt sphere_t3_elements.txt
mesh_bandwidth sphere_t3_elements.txt >> mesh_bandwidth_prb_output.txt
rm sphere_t3_elements.txt
#
cp ../../data/tet_mesh_order4/cube_order4_tetras.txt cube_order4_tetras.txt
mesh_bandwidth cube_order4_tetras.txt >> mesh_bandwidth_prb_output.txt
rm cube_order4_tetras.txt
#
cp ../../data/tet_mesh_order4/twenty_order4_tetras.txt twenty_order4_tetras.txt
mesh_bandwidth twenty_order4_tetras.txt >> mesh_bandwidth_prb_output.txt
rm twenty_order4_tetras.txt
#
cp ../../data/tet_mesh_order10/cube_order10_tetras.txt cube_order10_tetras.txt
mesh_bandwidth cube_order10_tetras.txt >> mesh_bandwidth_prb_output.txt
rm cube_order10_tetras.txt
#
cp ../../data/tet_mesh_order10/oneoneeight_order10_tetras.txt oneoneeight_order10_tetras.txt
mesh_bandwidth oneoneeight_order10_tetras.txt >> mesh_bandwidth_prb_output.txt
rm oneoneeight_order10_tetras.txt
#
cp ../../data/triangulation_order3/ell_tri3.txt ell_tri3.txt
mesh_bandwidth ell_tri3.txt >> mesh_bandwidth_prb_output.txt
rm ell_tri3.txt
#
cp ../../data/triangulation_order3/hex_cvt_tri3.txt hex_cvt_tri3.txt
mesh_bandwidth hex_cvt_tri3.txt  >> mesh_bandwidth_prb_output.txt
rm hex_cvt_tri3.txt
#
cp ../../data/triangulation_order3/hex_triangle_tri3.txt hex_triangle_tri3.txt
mesh_bandwidth hex_triangle_tri3.txt >> mesh_bandwidth_prb_output.txt
rm hex_triangle_tri3.txt
#
cp ../../data/triangulation_order3/hot_pipe_tri3.txt hot_pipe_tri3.txt
mesh_bandwidth hot_pipe_tri3.txt >> mesh_bandwidth_prb_output.txt
rm hot_pipe_tri3.txt
#
cp ../../data/triangulation_order6/ell_tri6.txt ell_tri6.txt
mesh_bandwidth ell_tri6.txt >> mesh_bandwidth_prb_output.txt
rm ell_tri6.txt
#
cp ../../data/triangulation_order6/hex_jeff_tri6.txt hex_jeff_tri6.txt
mesh_bandwidth hex_jeff_tri6.txt >> mesh_bandwidth_prb_output.txt
rm hex_jeff_tri6.txt
#
echo "Program output written to mesh_bandwidth_prb_output.txt"
