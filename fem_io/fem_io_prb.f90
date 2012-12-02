program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM_IO_PRB.
!
!  Discussion:
!
!    FEM_IO_PRB runs the FEM_IO tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM_IO_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FEM_IO library.'

  call test01 ( )
  call test02 ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM_IO_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests FEM_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  character ( len = 80 ) :: element_file_name = 'ell_elements.txt'
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_coord
  character ( len = 80 ) :: node_coord_file_name = 'ell_nodes.txt'
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_data
  character ( len = 80 ) :: node_data_file_name = 'ell_values.txt'  
  integer ( kind = 4 ) node_data_num
  integer ( kind = 4 ) node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  FEM_READ reads finite element data from files.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The node coordinate file name is "' // &
    trim ( node_coord_file_name ) // '".'
  write ( *, '(a)' ) '  The element file name is "' // &
    trim ( element_file_name ) // '".'
  write ( *, '(a)' ) '  The node data file name is "' // &
    trim ( node_data_file_name ) // '".'

  call fem_header_read ( node_coord_file_name, element_file_name, &
    node_data_file_name, dim_num, node_num, element_num, &
    element_order, node_data_num )

  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( node_coord(1:dim_num,1:node_num) )
  allocate ( node_data(1:node_data_num,1:node_num) )

  call fem_data_read ( node_coord_file_name, element_file_name, &
    node_data_file_name, dim_num, node_num, element_num, &
    element_order, node_data_num, node_coord, element_node, node_data )

  call fem_header_print ( dim_num, node_num, element_order, element_num, &
    node_data_num )

  call r8mat_transpose_print_some ( dim_num, node_num, node_coord, 1, 1, &
    dim_num,  10, '  First 10 node coordindates:' )

  call i4mat_transpose_print_some ( element_order, element_num, &
    element_node, 1, 1, element_order, 10, '  First 10 elements' )

  call r8mat_transpose_print_some ( node_data_num, node_num, node_data, &
    1, 1, node_data_num, 10, '  First 10 node data sets:' )

  deallocate ( element_node )
  deallocate ( node_coord )
  deallocate ( node_data )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! FEM_IO_TEST02 tests FEM_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: node_num = 5
  integer ( kind = 4 ), parameter :: element_num = 3
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ), parameter :: node_data_num = 2

  character ( len = 80 ) :: element_file_name = 'tiny_elements.txt'
  integer ( kind = 4 ), dimension (element_order,element_num) :: element_node = &
    reshape ( (/ &
    1, 2, 4, &
    5, 4, 2, &
    2, 3, 5 /), (/ element_order, element_num /) )
  real ( kind = 8 ), dimension (dim_num,node_num) :: node_coord = &
    reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    2.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00 /), (/ dim_num, node_num /) )
  character ( len = 80 ) :: node_coord_file_name = 'tiny_nodes.txt'
  real ( kind = 8 ), dimension (node_data_num,node_num) :: node_data = &
    reshape ( (/ &
    1.0D+00, 0.0D+00, &
    0.8D+00, 0.2D+00, &
    0.6D+00, 0.4D+00, &
    0.9D+00, 0.1D+00, &
    0.5D+00, 0.5D+00 /), (/ node_data_num, node_num /) )
  character ( len = 80 ) :: node_data_file_name = 'tiny_values.txt'  
    
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM_TEST02'
  write ( *, '(a)' ) '  Demonstrate the use of FEM_WRITE to write finite'
  write ( *, '(a)' ) '  element data to files.'

  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) '  The node coordinate file name is "' // &
    trim ( node_coord_file_name ) // '".'
  write ( *, '(a)' ) '  The element file name is "' // &
    trim ( element_file_name ) // '".'
  write ( *, '(a)' ) '  The node data file name is "' // &
    trim ( node_data_file_name ) // '".'

  call fem_header_print ( dim_num, node_num, element_order, element_num, &
    node_data_num )

  call r8mat_transpose_print ( dim_num, node_num, node_coord, &
    '  Node coordindates:' )

  call i4mat_transpose_print ( element_order, element_num, &
    element_node, '  Elements' )

  call r8mat_transpose_print ( node_data_num, node_num, node_data, &
    '  Node data sets:' )

  call fem_write ( node_coord_file_name, element_file_name, &
    node_data_file_name, dim_num, node_num, element_num, &
    element_order, node_data_num, node_coord, element_node, node_data )

  return
end
