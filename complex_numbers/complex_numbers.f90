program main

!*****************************************************************************80
!
!! COMPLEX_NUMBERS is a program which demonstrates the use of complex numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMPLEX_NUMBERS:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Demonstrate complex number usage.'
!
!  Single precision complex.
!
  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Double precision complex.
!
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMPLEX_NUMBERS:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 declaration and assignment.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none
!
!  Declare a complex number A.
!  Declare a complex vector B.
!  Declare a complex array C.
!
  complex a
  complex b(3)
  complex c(2,2)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Declare a COMPLEX variable.'
  write ( *, '(a)' ) '  Assign value with an = statement.'
!
!  Assign values to A, B, and C.
!
  a = ( 1.0, 2.0 )

  b(1) = ( 1.0, 2.0 )
  b(2) = ( 3.0, 4.0 )
  b(3) = ( 5.0, 6.0 )

  c(1,1) = ( 1.0, 0.1 )
  c(2,1) = ( 2.0, 0.1 )
  c(1,2) = ( 1.0, 0.2 )
  c(2,2) = ( 2.0, 0.2 )
!
!  Print them.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Scalar A:'
  write ( *, '(a)' ) ' '

  print *, a
  write ( *, * ) a
  write ( *, '(2g14.6)' ) a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vector B:'
  write ( *, '(a)' ) ' '

  do i = 1, 3
    write ( *, '(2g14.6)' ) b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Array C:'
  write ( *, '(a)' ) ' '

  do i = 1, 2
    write ( *, '(2(2g14.6))' ) c(i,1:2)
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02: declaration with initialization.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none
!
!  Declare and initialize a complex number A.
!  Declare and initialize a complex vector B.
!  Declare and initialize a complex array C.
!
  complex :: a = ( 1.0, 2.0 )
  complex :: b(3) = (/ &
    ( 1.0, 2.0 ), ( 3.0, 4.0 ), ( 5.0, 6.0 ) /)
  complex :: c(2,2) = reshape ( (/ &
    ( 1.0, 0.1 ), ( 2.0, 0.1 ), ( 1.0, 0.2 ), ( 2.0, 0.2 ) /), (/ 2, 2 /) )
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Declare a COMPLEX variable.'
  write ( *, '(a)' ) '  Initialize as part of the declaration.'
!
!  Print them.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Scalar A:'
  write ( *, '(a)' ) ' '

  print *, a
  write ( *, * ) a
  write ( *, '(2g14.6)' ) a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vector B:'
  write ( *, '(a)' ) ' '

  do i = 1, 3
    write ( *, '(2g14.6)' ) b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Array C:'
  write ( *, '(a)' ) ' '

  do i = 1, 2
    write ( *, '(2(2g14.6))' ) c(i,1:2)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03: intrinsic functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex a
  complex b(3)
  complex c(3,3)
  complex d(3)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Apply intrinsic functions to COMPLEX variables'

  a = ( 1.0, 2.0 )

  b = (/ ( 1.0, 2.0 ), ( 3.0, 4.0 ), ( 5.0, 6.0 ) /)

  c = reshape ( (/ &
    ( 1.0, 0.1 ), ( 2.0, 0.1 ), ( 3.0, 0.1 ), &
    ( 1.0, 0.2 ), ( 2.0, 0.2 ), ( 3.0, 0.2 ), &
    ( 1.0, 0.3 ), ( 2.0, 0.3 ), ( 3.0, 0.3 ) /), (/ 3, 3 /) )
!
!  Print them.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  a =              ', a
  write ( *, '(a,2g14.6)' ) '  - a =            ', - a
  write ( *, '(a,2g14.6)' ) '  a + 3 =          ', a + 3
  write ( *, '(a,2g14.6)' ) '  a + (0,5) =      ', a + ( 0, 5 )
  write ( *, '(a,2g14.6)' ) '  4 * a =          ', 4 * a
  write ( *, '(a,2g14.6)' ) '  a / 8 =          ', a / 8
  write ( *, '(a,2g14.6)' ) '  a * a =          ', a * a
  write ( *, '(a,2g14.6)' ) '  a**2 =           ', a**2
  write ( *, '(a,2g14.6)' ) '  1/a =            ', 1.0 / a
  write ( *, '(a)' ) ' '
  write ( *, '(a, g14.6)' ) '  abs(a) =         ', abs ( a )
  write ( *, '(a,2g14.6)' ) '  acos(a) =        ', acos ( a )
  write ( *, '(a,2g14.6)' ) '  asin(a) =        ', asin ( a )
  write ( *, '(a,2g14.6)' ) '  atan(a) =        ', atan ( a )
  write ( *, '(a,2g14.6)' ) '  cmplx(1) =       ', cmplx ( 1 )
  write ( *, '(a,2g14.6)' ) '  cmplx(2,3) =     ', cmplx ( 2, 3 )
  write ( *, '(a,2g14.6)' ) '  cmplx(4.0) =     ', cmplx ( 4.0 )
  write ( *, '(a,2g14.6)' ) '  cmplx(5.0,6.0) = ', cmplx ( 5.0, 6.0 )
  write ( *, '(a,2g14.6)' ) '  conjg(a) =       ', conjg ( a )
  write ( *, '(a,2g14.6)' ) '  cos(a) =         ', cos ( a )
  write ( *, '(a,2g14.6)' ) '  cosh(a) =        ', cosh ( a )
  write ( *, '(a,2g14.6)' ) '  exp(a) =         ', exp ( a )
  write ( *, '(a, g14.6)' ) '  imag(a) =        ', imag ( a )
  write ( *, '(a, i8)'    ) '  int(a) =         ', int ( a )
  write ( *, '(a,2g14.6)' ) '  log(a) =         ', log ( a )
  write ( *, '(a, g14.6)' ) '  real(a) =        ', real ( a )
  write ( *, '(a,2g14.6)' ) '  sin(a) =         ', sin ( a )
  write ( *, '(a,2g14.6)' ) '  sinh(a) =        ', sinh ( a )
  write ( *, '(a,2g14.6)' ) '  sqrt(a) =        ', sqrt ( a )
  write ( *, '(a,2g14.6)' ) '  tan(a) =         ', tan ( a )
  write ( *, '(a,2g14.6)' ) '  tanh(a) =        ', tanh ( a )
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  sum(b) =         ', sum ( b )
  write ( *, '(a,2g14.6)' ) '  product(b) =     ', product ( b )
  write ( *, '(a,2g14.6)' ) '  dot_product(b,b)=', dot_product ( b, b )

  d = matmul ( c, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  D = C * B = matmul ( c, b ):'
  write ( *, '(a)' ) ' '

  do i = 1, 3
    write ( *, '(2g14.6)' ) d(i)
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 declaration and assignment for double precision complex variables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none
!
!  Declare a double precision complex number A.
!  Declare a double precision complex vector B.
!  Declare a double precision complex array C.
!
  double complex a
  double complex b(3)
  double complex c(2,2)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Declare a DOUBLE COMPLEX variable.'
  write ( *, '(a)' ) '  Assign value with an = statement.'
!
!  Assign values to A, B, and C.
!
  a = ( 1.0D+00, 2.0D+00 )

  b(1) = ( 1.0D+00, 2.0D+00 )
  b(2) = ( 3.0D+00, 4.0D+00 )
  b(3) = ( 5.0D+00, 6.0D+00 )

  c(1,1) = ( 1.0D+00, 0.1D+00 )
  c(2,1) = ( 2.0D+00, 0.1D+00 )
  c(1,2) = ( 1.0D+00, 0.2D+00 )
  c(2,2) = ( 2.0D+00, 0.2D+00 )
!
!  Print them.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Scalar A:'
  write ( *, '(a)' ) ' '

  print *, a
  write ( *, * ) a
  write ( *, '(2g14.6)' ) a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vector B:'
  write ( *, '(a)' ) ' '

  do i = 1, 3
    write ( *, '(2g14.6)' ) b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Array C:'
  write ( *, '(a)' ) ' '

  do i = 1, 2
    write ( *, '(2(2g14.6))' ) c(i,1:2)
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05: declaration with initialization for double precision complex variables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none
!
!  Declare and initialize a double complex number A.
!  Declare and initialize a double complex vector B.
!  Declare and initialize a double complex array C.
!
  double complex :: a = ( 1.0D+00, 2.0D+00 )
  double complex :: b(3) = (/ &
    ( 1.0D+00, 2.0D+00 ), &
    ( 3.0D+00, 4.0D+00 ), &
    ( 5.0D+00, 6.0D+00 ) /)
  double complex :: c(2,2) = reshape ( (/ &
    ( 1.0D+00, 0.1D+00 ), &
    ( 2.0D+00, 0.1D+00 ), &
    ( 1.0D+00, 0.2D+00 ), &
    ( 2.0D+00, 0.2D+00 ) /), (/ 2, 2 /) )
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Declare a DOUBLE COMPLEX variable.'
  write ( *, '(a)' ) '  Initialize as part of the declaration.'
!
!  Print them.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Scalar A:'
  write ( *, '(a)' ) ' '

  print *, a
  write ( *, * ) a
  write ( *, '(2g14.6)' ) a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vector B:'
  write ( *, '(a)' ) ' '

  do i = 1, 3
    write ( *, '(2g14.6)' ) b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Array C:'
  write ( *, '(a)' ) ' '

  do i = 1, 2
    write ( *, '(2(2g14.6))' ) c(i,1:2)
  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST03: intrinsic functions for double precision complex variables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  double complex a
  double complex b(3)
  double complex c(3,3)
  double complex d(3)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Apply intrinsic functions to DOUBLE COMPLEX variable'

  a = ( 1.0D+00, 2.0D+00 )

  b = (/ ( 1.0D+00, 2.0D+00 ), ( 3.0D+00, 4.0D+00 ), ( 5.0D+00, 6.0D+00 ) /)

  c = reshape ( (/ &
    ( 1.0D+00, 0.1D+00 ), ( 2.0D+00, 0.1D+00 ), ( 3.0D+00, 0.1D+00 ), &
    ( 1.0D+00, 0.2D+00 ), ( 2.0D+00, 0.2D+00 ), ( 3.0D+00, 0.2D+00 ), &
    ( 1.0D+00, 0.3D+00 ), ( 2.0D+00, 0.3D+00 ), ( 3.0D+00, 0.3D+00 ) /), &
    (/ 3, 3 /) )
!
!  Print them.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  a =              ', a
  write ( *, '(a,2g14.6)' ) '  - a =            ', - a
  write ( *, '(a,2g14.6)' ) '  a + 3 =          ', a + 3
  write ( *, '(a,2g14.6)' ) '  a + (0,5) =      ', a + ( 0, 5 )
  write ( *, '(a,2g14.6)' ) '  4 * a =          ', 4 * a
  write ( *, '(a,2g14.6)' ) '  a / 3 =          ', a / 3
  write ( *, '(a,2g14.6)' ) '  a * a =          ', a * a
  write ( *, '(a,2g14.6)' ) '  a**2 =           ', a**2
  write ( *, '(a,2g14.6)' ) '  1/a =            ', 1.0 / a
  write ( *, '(a)' ) ' '
  write ( *, '(a, g14.6)' ) '  abs(a) =         ', abs ( a )
  write ( *, '(a,2g14.6)' ) '  acos(a) =        ', acos ( a )
  write ( *, '(a,2g14.6)' ) '  asin(a) =        ', asin ( a )
  write ( *, '(a,2g14.6)' ) '  atan(a) =        ', atan ( a )
  write ( *, '(a,2g14.6)' ) '  cmplx(1) =       ', cmplx ( 1 )
  write ( *, '(a,2g14.6)' ) '  cmplx(2,3) =     ', cmplx ( 2, 3 )
  write ( *, '(a,2g14.6)' ) '  cmplx(4.0) =     ', cmplx ( 4.0 )
  write ( *, '(a,2g14.6)' ) '  cmplx(5.0,6.0) = ', cmplx ( 5.0, 6.0 )
  write ( *, '(a,2g14.6)' ) '  conjg(a) =       ', conjg ( a )
  write ( *, '(a,2g14.6)' ) '  cos(a) =         ', cos ( a )
  write ( *, '(a,2g14.6)' ) '  cosh(a) =        ', cosh ( a )
  write ( *, '(a,2g14.6)' ) '  exp(a) =         ', exp ( a )
  write ( *, '(a, g14.6)' ) '  imag(a) =        ', imag ( a )
  write ( *, '(a, i8)'    ) '  int(a) =         ', int ( a )
  write ( *, '(a,2g14.6)' ) '  log(a) =         ', log ( a )
  write ( *, '(a, g14.6)' ) '  real(a) =        ', real ( a )
  write ( *, '(a,2g14.6)' ) '  sin(a) =         ', sin ( a )
  write ( *, '(a,2g14.6)' ) '  sinh(a) =        ', sinh ( a )
  write ( *, '(a,2g14.6)' ) '  sqrt(a) =        ', sqrt ( a )
  write ( *, '(a,2g14.6)' ) '  tan(a) =         ', tan ( a )
  write ( *, '(a,2g14.6)' ) '  tanh(a) =        ', tanh ( a )
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  sum(b) =         ', sum ( b )
  write ( *, '(a,2g14.6)' ) '  product(b) =     ', product ( b )
  write ( *, '(a,2g14.6)' ) '  dot_product(b,b)=', dot_product ( b, b )

  d = matmul ( c, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  D = C * B = matmul ( c, b ):'
  write ( *, '(a)' ) ' '

  do i = 1, 3
    write ( *, '(2g14.6)' ) d(i)
  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
