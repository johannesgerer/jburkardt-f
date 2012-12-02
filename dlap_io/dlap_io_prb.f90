program main

!*****************************************************************************80
!
!! MAIN is the main program for DLAP_IO_PRB.
!
!  Discussion:
!
!    DLAP_IO_PRB is the main program for the DLAP_IO tests.
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

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DLAP_IO_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the DLAP_IO library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DLAP_IO_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DLAP_FILE_WRITE.
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

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nelt = 11

  real ( kind = 8 ), dimension ( nelt ) :: a = (/ &
    51.0D+00, 12.0D+00, 11.0D+00, 33.0D+00, 15.0D+00, &
    53.0D+00, 55.0D+00, 22.0D+00, 35.0D+00, 44.0D+00, &
    21.0D+00 /)
  integer ( kind = 4 ), dimension ( nelt ) :: ia = (/ &
     5,  1,  1,  3,  1,  5,  5,  2,  3,  4,  2 /)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) :: irhs = 1
  integer ( kind = 4 ) :: isoln = 1
  integer ( kind = 4 ) :: isym = 0
  character ( len = 80 ) :: output_file_name = 'a5by5.dlap'
  integer ( kind = 4 ) output_file_unit
  integer ( kind = 4 ), dimension ( nelt ) :: ja = (/ &
     1,  2,  1,  3,  5,  3,  5,  2,  5,  4,  1 /)
  real ( kind = 8 ), dimension ( n ) :: rhs = (/ &
    110.0D+00, 65.0D+00, 274.0D+00, 176.0D+00, 485.0D+00 /)
  real ( kind = 8 ), dimension ( n ) :: soln = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  DLAP_FILE_WRITE writes a matrix in SLAP Triad format'
  write ( *, '(a)' ) '  to a DLAP sparse matrix file.'

  call dlap_file_print ( n, nelt, isym, irhs, isoln, ia, ja, a, rhs, &
    soln, '  The DLAP data to be written to the file.' )

  call get_unit ( output_file_unit )

  open ( unit = output_file_unit, file = output_file_name, &
    status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' &
      // trim ( output_file_name ) // '".'
    return
  end if

  call dlap_file_write ( n, nelt, isym, irhs, isoln, ia, ja, a, rhs, &
    soln, output_file_unit )

  close ( unit = output_file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Wrote the matrix data to "' &
    // trim ( output_file_name ) // '".'

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests DLAP_FILE_READ.
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

  integer ( kind = 4 ), parameter :: n_max = 5
  integer ( kind = 4 ), parameter :: nelt_max = 11

  real ( kind = 8 ), dimension ( nelt_max ) :: a
  integer ( kind = 4 ), dimension ( nelt_max ) :: ia
  character ( len = 80 ) :: input_file_name = 'a5by5.dlap'
  integer ( kind = 4 ) input_file_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) irhs
  integer ( kind = 4 ) isoln
  integer ( kind = 4 ) :: isym = 0
  integer ( kind = 4 ), dimension ( nelt_max ) :: ja
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt
  real ( kind = 8 ), dimension ( n_max ) :: rhs
  real ( kind = 8 ), dimension ( n_max ) :: soln

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  DLAP_FILE_READ reads a matrix from a SLAP sparse'
  write ( *, '(a)' ) '  matrix file into DLAP Triad format.'

  call get_unit ( input_file_unit )

  open ( unit = input_file_unit, file = input_file_name, &
    status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' &
      // trim ( input_file_name ) // '".'
    return
  end if

  call dlap_file_read ( n_max, nelt_max, n, nelt, isym, irhs, isoln, &
    ia, ja, a, rhs, soln, input_file_unit )

  close ( unit = input_file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Read the matrix data from "' &
    // trim ( input_file_name ) // '".'

  call dlap_file_print ( n, nelt, isym, irhs, isoln, ia, ja, a, rhs, &
    soln, '  The DLAP data read from the file.' )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests DLAP_FILE_WRITE.
!
!  Discussion:
!
!    The symmetric matrix is:
!
!      11   0  31   0  51
!       0  22   0   0   0
!      31   0  33  43   0
!       0   0  43  44  54
!      51   0   0  54  55
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

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nelt = 9

  real ( kind = 8 ), dimension ( nelt ) :: a = (/ &
    11.0D+00, 31.0D+00, 51.0D+00, 22.0D+00, 33.0D+00, &
    43.0D+00, 44.0D+00, 54.0D+00, 55.0D+00 /)
  integer ( kind = 4 ), dimension ( nelt ) :: ia = (/ &
    1, 3, 5, 2, 3, 4, 4, 5, 5  /)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) :: irhs = 1
  integer ( kind = 4 ) :: isoln = 1
  integer ( kind = 4 ) :: isym = 1
  character ( len = 80 ) :: output_file_name = 'a5by5_sym.dlap'
  integer ( kind = 4 ) output_file_unit
  integer ( kind = 4 ), dimension ( nelt ) :: ja = (/ &
     1, 1, 1, 2, 3, 3, 4, 4, 5 /)
  real ( kind = 8 ), dimension ( n ) :: rhs = (/ &
    359.0D+00, 44.0D+00, 302.0D+00, 575.0D+00, 542.0D+00 /)
  real ( kind = 8 ), dimension ( n ) :: soln = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  DLAP_FILE_WRITE writes a matrix in SLAP Triad format'
  write ( *, '(a)' ) '  to a DLAP sparse matrix file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, symmetric storage is used.'

  call dlap_file_print ( n, nelt, isym, irhs, isoln, ia, ja, a, rhs, &
    soln, '  The DLAP data to be written to the file.' )

  call get_unit ( output_file_unit )

  open ( unit = output_file_unit, file = output_file_name, &
    status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' &
      // trim ( output_file_name ) // '".'
    return
  end if

  call dlap_file_write ( n, nelt, isym, irhs, isoln, ia, ja, a, rhs, &
    soln, output_file_unit )

  close ( unit = output_file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Wrote the matrix data to "' &
    // trim ( output_file_name ) // '".'

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests DLAP_FILE_READ.
!
!  Discussion:
!
!    In this example, a matrix in symmetric storage will be read.
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

  integer ( kind = 4 ), parameter :: n_max = 5
  integer ( kind = 4 ), parameter :: nelt_max = 11

  real ( kind = 8 ), dimension ( nelt_max ) :: a
  integer ( kind = 4 ), dimension ( nelt_max ) :: ia
  character ( len = 80 ) :: input_file_name = 'a5by5_sym.dlap'
  integer ( kind = 4 ) input_file_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) irhs
  integer ( kind = 4 ) isoln
  integer ( kind = 4 ) :: isym = 0
  integer ( kind = 4 ), dimension ( nelt_max ) :: ja
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt
  real ( kind = 8 ), dimension ( n_max ) :: rhs
  real ( kind = 8 ), dimension ( n_max ) :: soln

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  DLAP_FILE_READ reads a matrix from a SLAP sparse'
  write ( *, '(a)' ) '  matrix file into DLAP Triad format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, symmetric storage is used.'

  call get_unit ( input_file_unit )

  open ( unit = input_file_unit, file = input_file_name, &
    status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' &
      // trim ( input_file_name ) // '".'
    return
  end if

  call dlap_file_read ( n_max, nelt_max, n, nelt, isym, irhs, isoln, &
    ia, ja, a, rhs, soln, input_file_unit )

  close ( unit = input_file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Read the matrix data from "' &
    // trim ( input_file_name ) // '".'

  call dlap_file_print ( n, nelt, isym, irhs, isoln, ia, ja, a, rhs, &
    soln, '  The DLAP data read from the file.' )

  return
end

