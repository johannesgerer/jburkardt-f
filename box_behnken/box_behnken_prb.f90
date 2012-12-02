program main

!*****************************************************************************80
!
!! MAIN calls a set of problems for BOX_BEHNKEN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BOX_BEHNKEN_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BOX_BEHNKEN library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BOX_BEHNKEN_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BOX_BEHNKEN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ), dimension (dim_num,2) :: range = reshape ( (/ &
    0.0D+00, 10.0D+00,  5.0D+00, &
    1.0D+00, 11.0D+00, 15.0D+00 /), (/ dim_num, 2 /) )
  integer ( kind = 4 ) x_num
  real ( kind = 8 ), allocatable, dimension(:,:) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  BOX_BEHNKEN computes a Box-Behnken dataset.'

  call r8mat_transpose_print ( dim_num, 2, range, &
    '  The ranges:' )

  call box_behnken_size ( dim_num, x_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a,i6)' ) '  For dimension DIM_NUM = ', dim_num, &
    ' the Box-Behnken design is of size ', x_num

  allocate ( x(1:dim_num,1:x_num) )

  call box_behnken ( dim_num, x_num, range, x )

  call r8mat_transpose_print ( dim_num, x_num, x, &
    '  The Box-Behnken design:' )

  deallocate ( x )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 writes a BOX_BEHNKEN dataset to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 4

  character ( len = 80 ) file_out_name
  real ( kind = 8 ), dimension (dim_num,2) :: range = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /), (/ dim_num, 2 /) )
  integer ( kind = 4 ) x_num
  real ( kind = 8 ), allocatable, dimension(:,:) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  R8MAT_WRITE writes a Box-Behnken dataset'
  write ( *, '(a)' ) '  to a file.'

  call r8mat_transpose_print ( dim_num, 2, range, &
    '  The ranges:' )

  call box_behnken_size ( dim_num, x_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a,i6)' ) '  For dimension DIM_NUM = ', dim_num, &
    ' the Box-Behnken design is of size ', x_num

  allocate ( x(1:dim_num,1:x_num) )

  call box_behnken ( dim_num, x_num, range, x )

  file_out_name = 'box_behnken_04_33.txt'

  call r8mat_write ( file_out_name, dim_num, x_num, x )

  deallocate ( x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data was written to the file "' // &
    trim ( file_out_name ) // '".'

  return
end

