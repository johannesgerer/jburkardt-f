program main

!*****************************************************************************80
!
!! MAIN is the main program for ST_IO_PRB.
!
!  Discussion:
!
!    ST_IO_PRB is the main program for the ST_IO tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ST_IO_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ST_IO library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ST_IO_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests ST_HEADER_WRITE and ST_DATA_WRITE.
!
!  Discussion:
!
!    The matrix is:
!
!      11  12   0   0  15
!      21  22   0   0   0
!       0   0  33   0  35
!       0   0   0  44   0
!      51   0  53   0  55
!
!    The index vectors are 1 based, and so have to be converted to
!    0-base before being written.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nrow = 5
  integer ( kind = 4 ), parameter :: ncol = 5
  integer ( kind = 4 ), parameter :: nnzero = 11

  real ( kind = 8 ), dimension ( nnzero ) :: a = (/ &
    51.0D+00, 12.0D+00, 11.0D+00, 33.0D+00, 15.0D+00, &
    53.0D+00, 55.0D+00, 22.0D+00, 35.0D+00, 44.0D+00, &
    21.0D+00 /)
  integer ( kind = 4 ) base0
  integer ( kind = 4 ) base1
  integer ( kind = 4 ), dimension ( nnzero ) :: col= (/ &
     1, 2, 1, 3, 5, 3, 5, 2, 5, 4, 1 /)
  integer ( kind = 4 ), dimension ( nnzero ) :: row = (/ &
     5, 1, 1, 3, 1, 5, 5, 2, 3, 4, 2 /)
  character ( len = 80 ) :: output_filename = 'a5by5.st'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  ST_HEADER_WRITE writes the header of an ST file;'
  write ( *, '(a)' ) '  ST_DATA_WRITE writes the data of an ST file.'

  base1 = 1
  base0 = 0
  call st_rebase ( base1, base0, nnzero, row )
  call st_rebase ( base1, base0, nnzero, col )

  call st_header_print ( nrow, ncol, nnzero )

  call st_data_print ( nrow, ncol, nnzero, row, col, a, &
    '  TEST01 matrix data to be written to a file:' )

  call st_header_write ( output_filename, nrow, ncol, nnzero )

  call st_data_write ( output_filename, nrow, ncol, nnzero, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Wrote the matrix data to "' &
    // trim ( output_filename ) // '".'

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests ST_HEADER_READ, ST_DATA_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( : ) :: a
  integer ( kind = 4 ), allocatable, dimension ( : ) :: col
  character ( len = 80 ) :: input_filename = 'kershaw.st'
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ), allocatable, dimension ( : ) :: row

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  ST_HEADER_READ reads the header from an ST file.'
  write ( *, '(a)' ) '  ST_DATA_READ reads the data from an ST file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data from "' // input_filename // '".'

  call st_header_read ( input_filename, nrow, ncol, nnzero )

  call st_header_print ( nrow, ncol, nnzero )

  allocate ( a(1:nnzero) )
  allocate ( col(1:nnzero) )
  allocate ( row(1:nnzero) )

  call st_data_read ( input_filename, nrow, ncol, nnzero, row, col, a )

  call st_data_print ( nrow, ncol, nnzero, row, col, a, &
    '  TEST02 matrix data read from file' )

  deallocate ( a )
  deallocate ( col )
  deallocate ( row )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests ST_SORT_A.
!
!  Discussion:
!
!    The matrix is:
!
!      11  12   0   0  15
!      21  22   0   0   0
!       0   0  33   0  35
!       0   0   0  44   0
!      51   0  53   0  55
!
!    The index vectors are 1 based, and so have to be converted to
!    0-base before being written.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nrow = 5
  integer ( kind = 4 ), parameter :: ncol = 5
  integer ( kind = 4 ), parameter :: nnzero = 11

  real ( kind = 8 ), dimension ( nnzero ) :: a = (/ &
    51.0D+00, 12.0D+00, 11.0D+00, 33.0D+00, 15.0D+00, &
    53.0D+00, 55.0D+00, 22.0D+00, 35.0D+00, 44.0D+00, &
    21.0D+00 /)
  integer ( kind = 4 ) base0
  integer ( kind = 4 ) base1
  integer ( kind = 4 ), dimension ( nnzero ) :: col= (/ &
     1, 2, 1, 3, 5, 3, 5, 2, 5, 4, 1 /)
  integer ( kind = 4 ), dimension ( nnzero ) :: row = (/ &
     5, 1, 1, 3, 1, 5, 5, 2, 3, 4, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  ST_SORT_A sorts an ST matrix by columns.'

  base1 = 1
  base0 = 0
  call st_rebase ( base1, base0, nnzero, row )
  call st_rebase ( base1, base0, nnzero, col )

  call st_header_print ( nrow, ncol, nnzero )

  call st_data_print ( nrow, ncol, nnzero, row, col, a, &
    '  Matrix data before sorting:' )

  call st_sort_a ( nrow, ncol, nnzero, row, col, a )

  call st_data_print ( nrow, ncol, nnzero, row, col, a, &
    '  Matrix data sorted by column:' )

  return
end
