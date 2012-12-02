program main

!*****************************************************************************80
!
!! MAIN is the main program for PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PIECEWISE_LINEAR_PRODUCT_INTEGRAL library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.
!
!  Discussion:
!
!    For the first test, we use the same single "piece" for both F and G.
!    Hence, we are actually integrating X^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 April 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 2
  integer ( kind = 4 ), parameter :: g_num = 2

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) exact
  real ( kind = 8 ) :: f_v(f_num) = (/ 0.0D+00, 5.0D+00 /)
  real ( kind = 8 ) :: f_x(f_num) = (/ 0.0D+00, 5.0D+00 /)
  real ( kind = 8 ) :: g_v(g_num) = (/ 0.0D+00, 5.0D+00 /)
  real ( kind = 8 ) :: g_x(g_num) = (/ 0.0D+00, 5.0D+00 /)
  integer ( kind = 4 ) i
  real ( kind = 8 ) integral

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) &
    '  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL on a very simple problem.'
  write ( *, '(a)' ) '  F and G are both defined over a single common'
  write ( *, '(a)' ) '  interval, so that F(X) = G(X) = X.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           A           B      Integral        Exact'
  write ( *, '(a)' ) ' '

  a = 1.0D+00
  do i = 1, 5
    b = real ( i, kind = 8 )
    call piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, g_num, &
      g_x, g_v, integral )
    exact = ( b * b * b - a * a * a ) / 3.0D+00
    write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' ) a, b, integral, exact
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.
!
!  Discussion:
!
!    For this test, we use multiple "pieces" for both F and G,
!    but we define the values so that we are still actually integrating X^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 April 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 3
  integer ( kind = 4 ), parameter :: g_num = 4

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) exact
  real ( kind = 8 ) :: f_v(f_num) = (/ 0.0D+00, 2.0D+00, 5.0D+00 /)
  real ( kind = 8 ) :: f_x(f_num) = (/ 0.0D+00, 2.0D+00, 5.0D+00 /)
  real ( kind = 8 ) :: g_v(g_num) = (/ 0.0D+00, 1.5D+00, 3.0D+00, 5.0D+00 /)
  real ( kind = 8 ) :: g_x(g_num) = (/ 0.0D+00, 1.5D+00, 3.0D+00, 5.0D+00 /)
  integer ( kind = 4 ) i
  real ( kind = 8 ) integral

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) &
    '  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL on a simple problem.'
  write ( *, '(a)' ) '  F and G are both defined over separate, multiple'
  write ( *, '(a)' ) '  intervals, but still true that F(X) = G(X) = X.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           A           B      Integral        Exact'
  write ( *, '(a)' ) ' '

  a = 1.0D+00
  do i = 1, 5
    b = real ( i, kind = 8 )
    call piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, g_num, &
      g_x, g_v, integral )
    exact = ( b * b * b - a * a * a ) / 3.0D+00
    write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' ) a, b, integral, exact
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.
!
!  Discussion:
!
!    For this test, F(X) and G(X) are piecewise linear interpolants to
!    SIN(X) and 2 * COS(X), so we know the exact value of the integral
!    of the product of the original functions, but this is only an estimate
!    of the exact value of the integral of the product of the interpolants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 11
  integer ( kind = 4 ), parameter :: g_num = 31

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) exact
  real ( kind = 8 ) :: f_v(f_num)
  real ( kind = 8 ) :: f_x(f_num)
  real ( kind = 8 ) :: g_v(g_num)
  real ( kind = 8 ) :: g_x(g_num)
  integer ( kind = 4 ) i
  real ( kind = 8 ) integral
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  integer ( kind = 4 ) quad_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) &
    '  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL on a simple problem.'
  write ( *, '(a)' ) '  F and G are defined over separate, multiple'
  write ( *, '(a)' ) '  intervals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X) interpolates SIN(X),'
  write ( *, '(a)' ) '  G(X) interpolates 2*COS(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We compare:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  INTEGRAL, our value for the integral,'
  write ( *, '(a)' ) '  QUAD, a quadrature estimate for the integral, and'
  write ( *, '(a)' ) '  CLOSE, the value of the integral of 2*COS(X)*SIN(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '           A           B      Integral        Quad            Close'
  write ( *, '(a)' ) ' '

  do i = 1, f_num
    f_x(i) = ( real ( f_num - i,     kind = 8     ) * 0.0D+00 &
             + real (         i - 1, kind = 8 ) * pi ) &
             / real ( f_num     - 1, kind = 8 )
    f_v(i) = sin ( f_x(i) )
  end do

  do i = 1, g_num
    g_x(i) = ( real ( g_num - i,     kind = 8     ) * 0.0D+00 &
             + real (         i - 1, kind = 8 ) * pi ) &
             / real ( g_num     - 1, kind = 8 )
    g_v(i) = 2.0D+00 * cos ( g_x(i) )
  end do

  a = 0.0D+00
  do i = 0, 6
    b = real ( i, kind = 8 ) * pi / 6.0D+00
    call piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, g_num, &
      g_x, g_v, integral )
    exact = - ( cos ( 2.0D+00 * b ) - cos ( 2.0D+00 * a ) ) / 2.0D+00
    quad_num = 2000
    call piecewise_linear_product_quad ( a, b, f_num, f_x, f_v, g_num, &
      g_x, g_v, quad_num, quad )
    write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      a, b, integral, quad, exact
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.
!
!  Discussion:
!
!    For this test, we compute the integrals of a hat function with itself,
!    and a hat function with its neighbor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 3
  integer ( kind = 4 ), parameter :: g_num = 3

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) exact
  real ( kind = 8 ) :: f_v(f_num) = (/ 0.0D+00, 1.0D+00, 0.0D+00 /)
  real ( kind = 8 ) :: f_x(f_num) = (/ 0.0D+00, 1.0D+00, 2.0D+00 /)
  real ( kind = 8 ) :: g_v(g_num) = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) :: g_x(g_num) = (/ 0.0D+00, 1.0D+00, 2.0D+00 /)
  integer ( kind = 4 ) i
  real ( kind = 8 ) integral

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL.'
  write ( *, '(a)' ) '  The nodes are at 0, 1, and 2.'
  write ( *, '(a)' ) '  F(X) = ( 0, 1, 0 ).'
  write ( *, '(a)' ) '  G(X) = ( 1, 0, 0 ).'
  write ( *, '(a)' ) ' '

  a = 0.0D+00
  b = 2.0D+00

  call piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, f_num, &
    f_x, f_v, integral )

  write ( *, '(a,g14.6)' ) '  Integral F(X) * F(X) dx = ', integral

  call piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, g_num, &
    g_x, g_v, integral )

  write ( *, '(a,g14.6)' ) '  Integral F(X) * G(X) dx = ', integral

  call piecewise_linear_product_integral ( a, b, g_num, g_x, g_v, g_num, &
    g_x, g_v, integral )

  write ( *, '(a,g14.6)' ) '  Integral G(X) * G(X) dx = ', integral

  return
end
