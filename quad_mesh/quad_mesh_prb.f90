program main

!*****************************************************************************80
!
!! MAIN is the main program for QUAD_MESH_PRB.
!
!  Discussion:
!
!    QUAD_MESH_PRB tests routines from the QUAD_MESH library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUAD_MESH_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the QUAD_MESH library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test062 ( )
  call test07 ( )
  call test072 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )
  call test105 ( )
  call test11 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUAD_MESH_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests ADJ_SIZE_Q4_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: adj_col
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) k
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  ADJ_SIZE counts the node adjacencies.'
!
!  Get the sizes of the example.
!
  call example1_q4_mesh_size ( node_num, element_num, hole_num )

  allocate ( element_neighbor(4,element_num) )
  allocate ( element_node(4,element_num) )
  allocate ( node_xy(2,1:node_num) )

  call example1_q4_mesh ( node_num, element_num, node_xy, element_node, &
    element_neighbor )
!
!  Get the count of the node adjacencies.
!
  allocate ( adj_col(node_num+1) )

  call adj_size_q4_mesh ( node_num, element_num, element_node, &
    element_neighbor, adj_num, adj_col )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of adjacency entries is ', adj_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Adjacency pointers:'
  write ( *, '(a)' ) ' '
  do node = 1, node_num
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) node, adj_col(node), adj_col(node+1)-1
  end do

  deallocate ( adj_col )
  deallocate ( node_xy )
  deallocate ( element_neighbor )
  deallocate ( element_node )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests AREA_Q4_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: element_area(:)
  integer ( kind = 4 ), allocatable :: element_neighbor(:,:)
  integer ( kind = 4 ), allocatable :: element_node(:,:)
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) hole_num
  real ( kind = 8 ) mesh_area
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable :: node_xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  AREA_Q4_MESH computes the area of each element'
  write ( *, '(a)' ) '  in a Q4 mesh.'

  call example2_q4_mesh_size ( node_num, element_num, hole_num )

  allocate ( element_neighbor(4,element_num) )
  allocate ( element_node(4,element_num) )
  allocate ( node_xy(2,1:node_num) )

  call example2_q4_mesh ( node_num, element_num, node_xy, element_node, &
    element_neighbor )

  allocate ( element_area(element_num) )

  call area_q4_mesh ( node_num, element_num, node_xy, element_node, &
    element_area, mesh_area )

  call r8vec_print ( element_num, element_area, '  Element areas:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g16.8)' ) '   Mesh =   ', mesh_area

  deallocate ( element_area )
  deallocate ( element_neighbor )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 demonstrates AREA_QUAD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) :: quad_xy(2,4) = reshape ( (/ &
    1.0D+00, 2.0D+00, &
    5.0D+00, 2.0D+00, &
    5.0D+00, 3.0D+00, &
    4.0D+00, 4.0D+00 /), (/ 2, 4 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  AREA_QUAD computes the area of a quadrilateral.'

  call area_quad ( quad_xy, area )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g16.8)' ) '  Area =    ', area

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests BOUNDARY_EDGE_COUNT_Q4_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) boundary_edge_num
  integer ( kind = 4 ), allocatable :: element_neighbor(:,:)
  integer ( kind = 4 ), allocatable :: element_node(:,:)
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable :: node_xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  BOUNDARY_EDGE_COUNT_Q4_MESH counts the'
  write ( *, '(a)' ) '  boundary edges by looking at the mesh'
  write ( *, '(a)' ) '  and counting unpaired edges.'
!
!  The value of HOLE_NUM, which is supplied by this call, is not
!  needed for this calculation.
!
  call example1_q4_mesh_size ( node_num, element_num, hole_num )

  allocate ( element_neighbor(4,element_num) )
  allocate ( element_node(4,element_num) )
  allocate ( node_xy(2,1:node_num) )

  call example1_q4_mesh ( node_num, element_num, node_xy, element_node, &
    element_neighbor )

  call boundary_edge_count_q4_mesh ( element_num, element_node, &
    boundary_edge_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary edges = ', boundary_edge_num
  write ( *, '(a,i8)' ) '  Correct number =           ', 22

  deallocate ( element_neighbor )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests BOUNDARY_EDGE_COUNT_EULER_Q4_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) boundary_edge_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  BOUNDARY_EDGE_COUNT_EULER_Q4_MESH counts the'
  write ( *, '(a)' ) '    boundary edges using Euler''s formula.'

  call example1_q4_mesh_size ( node_num, element_num, hole_num )

  call boundary_edge_count_euler_q4_mesh ( node_num, element_num, hole_num, &
    boundary_edge_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary edges = ', boundary_edge_num
  write ( *, '(a,i8)' ) '  Correct number =           ', 22

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests EXAMPLE1_Q4_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: element_area(:)
  integer ( kind = 4 ), allocatable :: element_neighbor(:,:)
  integer ( kind = 4 ), allocatable :: element_node(:,:)
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) hole_num
  real ( kind = 8 ) mesh_area
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable :: node_xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  EXAMPLE1_Q4_MESH sets up example #1 Q4 mesh.'

  call example1_q4_mesh_size ( node_num, element_num, hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes =    ', node_num
  write ( *, '(a,i8)' ) '  Number of elements = ', element_num
  write ( *, '(a,i8)' ) '  Number of holes =    ', hole_num

  allocate ( element_neighbor(4,element_num) )
  allocate ( element_node(4,element_num) )
  allocate ( node_xy(2,1:node_num) )

  call example1_q4_mesh ( node_num, element_num, node_xy, element_node, &
    element_neighbor )

  call r8mat_transpose_print ( 2, node_num, node_xy, '  Node coordinates:' )
  call i4mat_transpose_print ( 4, element_num, element_node, '  Elements:' )
  call i4mat_transpose_print ( 4, element_num, element_neighbor, &
    '  Element neighbors' )

  deallocate ( element_neighbor )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test062 ( )

!*****************************************************************************80
!
!! TEST062 tests EXAMPLE1_Q4_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: element_area(:)
  integer   ( kind = 4 ), allocatable :: element_neighbor(:,:)
  integer   ( kind = 4 ), allocatable :: element_node(:,:)
  integer   ( kind = 4 ) element_num
  integer   ( kind = 4 ) :: element_show = 2
  integer   ( kind = 4 ) hole_num
  real ( kind = 8 ) mesh_area
  integer   ( kind = 4 ) node_num
  integer   ( kind = 4 ) :: node_show = 2
  real ( kind = 8 ), allocatable :: node_xy(:,:)
  character ( len = 100 ) :: output_filename = 'q4_mesh_ex1.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST062'
  write ( *, '(a)' ) '  EXAMPLE1_Q4_MESH sets up example #1 Q4 mesh.'
  write ( *, '(a)' ) '  PLOT_Q4_MESH plots it.'

  call example1_q4_mesh_size ( node_num, element_num, hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes =    ', node_num
  write ( *, '(a,i8)' ) '  Number of elements = ', element_num
  write ( *, '(a,i8)' ) '  Number of holes =    ', hole_num

  allocate ( element_neighbor(4,element_num) )
  allocate ( element_node(4,element_num) )
  allocate ( node_xy(2,1:node_num) )

  call example1_q4_mesh ( node_num, element_num, node_xy, element_node, &
    element_neighbor )

  call plot_q4_mesh ( node_num, element_num, node_xy, element_node, &
    node_show, element_show, output_filename )

  deallocate ( element_neighbor )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests EXAMPLE2_Q4_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: element_area(:)
  integer ( kind = 4 ), allocatable :: element_neighbor(:,:)
  integer ( kind = 4 ), allocatable :: element_node(:,:)
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) hole_num
  real ( kind = 8 ) mesh_area
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable :: node_xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  EXAMPLE2_Q4_MESH sets up example #2 Q4 mesh.'

  call example2_q4_mesh_size ( node_num, element_num, hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes =    ', node_num
  write ( *, '(a,i8)' ) '  Number of elements = ', element_num
  write ( *, '(a,i8)' ) '  Number of holes =    ', hole_num

  allocate ( element_neighbor(4,element_num) )
  allocate ( element_node(4,element_num) )
  allocate ( node_xy(2,1:node_num) )

  call example2_q4_mesh ( node_num, element_num, node_xy, element_node, &
    element_neighbor )

  call r8mat_transpose_print ( 2, node_num, node_xy, '  Node coordinates:' )
  call i4mat_transpose_print ( 4, element_num, element_node, '  Elements:' )
  call i4mat_transpose_print ( 4, element_num, element_neighbor, &
    '  Element neighbors' )

  deallocate ( element_neighbor )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test072 ( )

!*****************************************************************************80
!
!! TEST072 tests EXAMPLE2_Q4_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: element_area(:)
  integer   ( kind = 4 ), allocatable :: element_neighbor(:,:)
  integer   ( kind = 4 ), allocatable :: element_node(:,:)
  integer   ( kind = 4 ) element_num
  integer   ( kind = 4 ) :: element_show = 2
  integer   ( kind = 4 ) hole_num
  real ( kind = 8 ) mesh_area
  integer   ( kind = 4 ) node_num
  integer   ( kind = 4 ) :: node_show = 2
  real ( kind = 8 ), allocatable :: node_xy(:,:)
  character ( len = 100 ) :: output_filename = 'q4_mesh_ex2.eps'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST072'
  write ( *, '(a)' ) '  EXAMPLE2_Q4_MESH sets up example #2 Q4 mesh.'
  write ( *, '(a)' ) '  PLOT_Q4_MESH plots it.'

  call example2_q4_mesh_size ( node_num, element_num, hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes =    ', node_num
  write ( *, '(a,i8)' ) '  Number of elements = ', element_num
  write ( *, '(a,i8)' ) '  Number of holes =    ', hole_num

  allocate ( element_neighbor(4,element_num) )
  allocate ( element_node(4,element_num) )
  allocate ( node_xy(2,1:node_num) )

  call example2_q4_mesh ( node_num, element_num, node_xy, element_node, &
    element_neighbor )

  call plot_q4_mesh ( node_num, element_num, node_xy, element_node, &
    node_show, element_show, output_filename )

  deallocate ( element_neighbor )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests NEIGHBOR_ELEMENTS_Q4_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), allocatable :: element_neighbor(:,:)
  integer   ( kind = 4 ), allocatable :: element_neighbor2(:,:)
  integer   ( kind = 4 ), allocatable :: element_node(:,:)
  integer   ( kind = 4 )  element_num
  integer   ( kind = 4 )  hole_num
  integer   ( kind = 4 )  node_num
  real ( kind = 8 ), allocatable :: node_xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  NEIGHBOR_ELEMENTS_Q4_MESH determines the'
  write ( *, '(a)' ) '  adjacency relationships between elements.'

  call example1_q4_mesh_size ( node_num, element_num, hole_num )

  allocate ( element_neighbor(4,element_num) )
  allocate ( element_node(4,element_num) )
  allocate ( node_xy(2,1:node_num) )

  call example1_q4_mesh ( node_num, element_num, node_xy, element_node, &
    element_neighbor )

  call i4mat_transpose_print ( 4, element_num, element_neighbor, &
    '  Element neighbors as reported by EXAMPLE1_Q4_MESH:' )

  allocate ( element_neighbor2(4,element_num) )

  call neighbor_elements_q4_mesh ( element_num, element_node, &
    element_neighbor2 )

  call i4mat_transpose_print ( 4, element_num, element_neighbor2, &
    '  Element neighbors computed by NEIGHBOR_ELEMENTS_Q4_MESH:' )

  deallocate ( element_neighbor )
  deallocate ( element_neighbor2 )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 writes data to files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), allocatable :: element_neighbor(:,:)
  integer   ( kind = 4 ), allocatable :: element_node(:,:)
  integer   ( kind = 4 )  element_num
  integer   ( kind = 4 )  hole_num
  integer   ( kind = 4 )  node_num
  real ( kind = 8 ), allocatable :: node_xy(:,:)
  character ( len = 100 ) output_filename

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  Write Q4 Mesh Example #1 to files.'

  call example2_q4_mesh_size ( node_num, element_num, hole_num )

  allocate ( element_neighbor(4,element_num) )
  allocate ( element_node(4,element_num) )
  allocate ( node_xy(2,1:node_num) )

  call example2_q4_mesh ( node_num, element_num, node_xy, element_node, &
    element_neighbor )

  output_filename = 'q4_mesh_ex2_element_neighbors.txt'
  call i4mat_write ( output_filename, 4, element_num, element_neighbor )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Element neighbors written to "' &
    // trim ( output_filename ) // '".'

  output_filename = 'q4_mesh_ex2_elements.txt'
  call i4mat_write ( output_filename, 4, element_num, element_node )
  write ( *, '(a)' ) '  Elements written to "' &
    // trim ( output_filename ) // '".'

  output_filename = 'q4_mesh_ex2_xy.txt'
  call r8mat_write ( output_filename, 2, node_num, node_xy )
  write ( *, '(a)' ) '  Node coordinates written to "' &
    // trim ( output_filename ) // '".'

  deallocate ( element_neighbor )
  deallocate ( element_node )
  deallocate ( node_xy )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests NODE_ORDER_Q4_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: element_neighbor(:,:)
  integer ( kind = 4 ), allocatable :: element_node(:,:)
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ), allocatable :: node_order(:)
  real ( kind = 8 ), allocatable :: node_xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  NODE_ORDER_4_MESH computes the order'
  write ( *, '(a)' ) '  of the nodes in a Q4 mesh.'

  call example1_q4_mesh_size ( node_num, element_num, hole_num )

  allocate ( element_neighbor(4,element_num) )
  allocate ( element_node(4,element_num) )
  allocate ( node_order(1:node_num) )
  allocate ( node_xy(2,1:node_num) )

  call example1_q4_mesh ( node_num, element_num, node_xy, element_node, &
    element_neighbor )

  call node_order_q4_mesh ( element_num, element_node, node_num, node_order )

  call i4vec_print ( node_num, node_order, '      NODE         ORDER' )

  deallocate ( element_neighbor )
  deallocate ( element_node )
  deallocate ( node_order )
  deallocate ( node_xy )

  return
end
subroutine test105 ( )

!*****************************************************************************80
!
!! TEST105 tests SAMPLE_Q4_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: element_neighbor(:,:)
  integer ( kind = 4 ), allocatable :: element_node(:,:)
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable :: node_xy(:,:)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ), allocatable :: sample_element(:)
  real ( kind = 8 ), allocatable :: sample_xy(:,:)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST105'
  write ( *, '(a)' ) '  SAMPLE_Q4_MESH returns uniform sample points from'
  write ( *, '(a)' ) '  a Q4 mesh.'

  call example1_q4_mesh_size ( node_num, element_num, hole_num )

  allocate ( element_neighbor(4,element_num) )
  allocate ( element_node(4,element_num) )
  allocate ( node_xy(2,1:node_num) )

  call example1_q4_mesh ( node_num, element_num, node_xy, element_node, &
    element_neighbor )

  sample_num = 20

  allocate ( sample_xy(1:2,1:sample_num) )
  allocate ( sample_element(1:sample_num) )

  seed = 123456789

  call sample_q4_mesh ( node_num, node_xy, element_num, element_node, &
    sample_num, seed, sample_xy, sample_element )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X        Y     Element'
  write ( *, '(a)' ) ' '

  do sample = 1, sample_num

    write ( *, '(2x,i8,2x,f8.4,2x,f8.4,2x,i8)' ) &
      sample, sample_xy(1:2,sample), sample_element(sample)

  end do

  deallocate ( element_neighbor )
  deallocate ( element_node )
  deallocate ( node_xy )
  deallocate ( sample_element )
  deallocate ( sample_xy )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 demonstrates SAMPLE_QUAD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 )  n
  character ( len = 100 ) output_filename
  real ( kind = 8 ) :: quad_xy(2,4) = reshape ( (/ &
    1.0D+00, 2.0D+00, &
    5.0D+00, 2.0D+00, &
    5.0D+00, 3.0D+00, &
    4.0D+00, 4.0D+00 /), (/ 2, 4 /) )
  integer   ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  SAMPLE_QUAD computes N random points in a quadrilateral.'
  write ( *, '(a)' ) '  Write them to a file.'

  n = 5000

  allocate ( xy(2,n) )

  seed = 123456789

  call sample_quad ( quad_xy, n, seed, xy )

  output_filename = 'sample_quad.txt'
  call r8mat_write ( output_filename, 2, n, xy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Point coordinates written to "' &
    // trim ( output_filename ) // '".'

  deallocate ( xy )

  return
end
