program main

!*****************************************************************************80
!
!! MAIN is the main program for POLY_IO_PRB.
!
!  Discussion:
!
!    POLY_IO_PRB runs the tests of the POLY_IO routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLY_IO_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the POLY_IO library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLY_IO_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests POLY_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: edge_num = 16
  integer ( kind = 4 ), parameter :: hole_num = 2
  integer ( kind = 4 ), parameter :: node_num = 16

  integer ( kind = 4 ), dimension(2,edge_num) :: edge_nodes = reshape ( (/ &
     1,  2, &
     2,  3, &
     3,  4, &
     4,  1, &
     5,  6, &
     6,  7, &
     7,  8, &
     8,  9, &
     9, 10, &
    10,  5, &
    11, 12, &
    12, 13, &
    13, 14, &
    14, 15, &
    15, 16, &
    16, 11 /), (/ 2, edge_num /) )
  real ( kind = 8 ), dimension(2,hole_num) :: hole_point = reshape ( (/ &
    0.25D+00, 0.75D+00, &
    0.60D+00, 0.40D+00 /), (/ 2, hole_num /) )
  character ( len = 80 ) :: poly_file_name = 'shape.poly'
  real ( kind = 8 ), dimension(2,node_num) :: segment = reshape ( (/ &
    0.000D+00, 0.000D+00, &
    1.000D+00, 0.000D+00, &
    1.000D+00, 1.000D+00, &
    0.000D+00, 1.000D+00, &
    0.350D+00, 0.750D+00, &
    0.300D+00, 0.836D+00, &
    0.200D+00, 0.836D+00, &
    0.150D+00, 0.750D+00, &
    0.200D+00, 0.663D+00, &
    0.300D+00, 0.663D+00, &
    0.700D+00, 0.400D+00, &
    0.650D+00, 0.486D+00, &
    0.550D+00, 0.486D+00, &
    0.500D+00, 0.400D+00, &
    0.550D+00, 0.313D+00, &
    0.650D+00, 0.313D+00 /), (/ 2, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  POLY_WRITE writes data to a POLY file.'

  call poly_write ( poly_file_name, node_num, segment, edge_num, &
    edge_nodes, hole_num, hole_point )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POLY_WRITE created the file "' &
    // trim ( poly_file_name ) // '".'

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests POLY_HEADER_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  character ( len = 80 ) :: poly_file_name = 'shape.poly'
  integer ( kind = 4 ) region_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  POLY_HEADER_READ reads header data from a POLY file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We are going to read from the file "' &
    // trim ( poly_file_name ) // '".'

  call poly_header_read ( poly_file_name, node_num, edge_num, hole_num, &
    region_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes =   ', node_num
  write ( *, '(a,i8)' ) '  Number of edges =   ', edge_num
  write ( *, '(a,i8)' ) '  Number of holes =   ', hole_num
  write ( *, '(a,i8)' ) '  Number of regions = ', region_num

  return
end
