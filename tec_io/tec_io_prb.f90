program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_IO_PRB.
!
!  Discussion:
!
!    TEC_IO_PRB tests the TEC_IO routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEC_IO_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEC_IO library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEC_IO_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests TEC_DATA_READ and TEC_HEADER_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: node_coord
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: node_data
  integer ( kind = 4 ) node_data_num
  integer ( kind = 4 ) node_num
  character ( len = 80 ) :: tec_file_name = 'ell.dat';
  integer ( kind = 4 ) tec_file_unit

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEC_IO_TEST01'
  write ( *, '(a)' ) '  TEC_HEADER_READ can read the header of a TEC file.'
  write ( *, '(a)' ) '  TEC_DATA_READ can read the data of a TEC file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we will read data from "' &
    // trim ( tec_file_name ) // '".'

  call tec_open_read ( tec_file_name, tec_file_unit )

  call tec_header_read ( tec_file_name, tec_file_unit, dim_num, node_num, &
    element_num, element_order, node_data_num )

  call tec_header_print ( dim_num, node_num, element_num, &
    element_order, node_data_num )

  allocate ( node_coord(1:dim_num,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( node_data(1:node_data_num,1:node_num) )

  call tec_data_read ( tec_file_name, tec_file_unit, dim_num, &
    node_num, element_num, element_order, node_data_num, node_coord, &
    element_node, node_data )

  close ( unit = tec_file_unit )

  call r8mat_transpose_print_some ( dim_num, node_num, node_coord, 1, 1, &
    dim_num, 10, '  Coordinates of first 10 nodes:' )

  call i4mat_transpose_print_some ( element_order, element_num, &
    element_node, 1, 1, element_order, 10, '  Nodes of first 10 elements:' )

  call r8mat_transpose_print_some ( node_data_num, node_num, node_data, 1, 1, &
    node_data_num, 10, '  Node data for first 10 nodes:' )

  deallocate ( node_coord )
  deallocate ( element_node )
  deallocate ( node_data )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests TEC_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: element_num = 3
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ), parameter :: node_data_num = 2
  integer ( kind = 4 ), parameter :: node_num = 5

  integer ( kind = 4 ), dimension ( element_order, element_num ) :: element_node = &
    reshape ( (/ &
    1, 2, 4, &
    5, 4, 2, &
    2, 3, 5 /), (/ element_order, element_num /) )
  real    ( kind = 8 ), dimension ( dim_num, node_num ) :: node_coord = &
    reshape ( (/ &
    0.0, 0.0, &
    1.0, 0.0, &
    2.0, 0.0, &
    0.0, 1.0, &
    1.0, 1.0 /), (/ dim_num, node_num /) )
  real    ( kind = 8 ), dimension ( node_data_num, node_num ) :: node_data = &
    reshape ( (/ &
    1.0, 0.0, &
    0.8, 0.2, &
    0.6, 0.4, &
    0.9, 0.1, &
    0.5, 0.5 /), (/ node_data_num, node_num /) )
  character ( len = 80 ) :: tec_file_name = 'tiny.dat'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEC_IO_TEST02'
  write ( *, '(a)' ) '  TEC_WRITE can write finite element data to a TEC file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we will write data to "' &
    // trim ( tec_file_name ) // '".'

  call tec_header_print ( dim_num, node_num, element_num, &
    element_order, node_data_num )

  call r8mat_transpose_print ( dim_num, node_num, node_coord, &
    '  Coordinates of nodes:' )

  call i4mat_transpose_print ( element_order, element_num, element_node, &
    '  Nodes of elements:' )

  call r8mat_transpose_print ( node_data_num, node_num, node_data, &
    '  Node data for nodes:' )

  call tec_write ( tec_file_name, dim_num, node_num, element_num, &
    element_order, node_data_num, node_coord, element_node, node_data )

  return
end
