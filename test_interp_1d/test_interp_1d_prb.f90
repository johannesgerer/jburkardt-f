program main

!*****************************************************************************80
!
!! TEST_INTERP_1D_TEST tests the TEST_INTERP_1D library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) nd

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INTERP_1D_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_INTERP_1D library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'

  call test01 ( )

  nd = 11
  call test02 ( nd )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INTERP_1D_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return;
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 simply prints the title of each function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Print the title of each function.'

  call p00_prob_num ( prob_num )
  
  write ( *, '(a)' ) ' '
  write ( *, '(a,i2,a)' ) '  There are ', prob_num, ' functions available:'
  write ( *, '(a)' ) '  Index  Title'
  write ( *, '(a)' ) ' '

  do prob = 1, prob_num

    call p00_title ( prob, title )

    write ( *, '(2x,i2,2x,a)' ) prob, trim ( title )

  end do

  return
end
subroutine test02 ( nd )

!*****************************************************************************80
!
!! TEST_INTERP_1D_TEST02 evaluates each test function at ND sample points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ND, the number of sample points.
!
  implicit none

  integer ( kind = 4 ) nd

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  character ( len = 80 ) title
  real ( kind = 8 ) x(nd)
  real ( kind = 8 ) y(nd)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INTERP_1D_TEST02'
  write ( *, '(a)' ) '  Use P00_F to sample each function.'

  call p00_prob_num ( prob_num )

  a = 0.0D+00
  b = 1.0D+00
  call r8vec_linspace ( nd, a, b, x )

  write ( *, '(a)' ) ' '

  do prob = 1, prob_num

    call p00_f ( prob, nd, x, y )
    write ( title, '(a,i2)' ) 'X, Y for problem ', prob
    call r8vec2_print ( nd, x, y, trim ( title ) )

  end do

  return
end
