program main

!*****************************************************************************80
!
!! CHEBYSHEV_TEST tests CHEBYSHEV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBYSHEV_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the CHEBYSHEV library.'

  call chebyshev_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBYSHEV_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine chebyshev_test01 ( )

!*****************************************************************************80
!
!! CHEBYSHEV_TEST01 tests CHEBYSHEV_COEFFICIENTS and CHEBYSHEV_INTERPOLANT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable :: c(:)
  real ( kind = 8 ), external :: f1
  real ( kind = 8 ), external :: f2
  real ( kind = 8 ), external :: f3
  real ( kind = 8 ), allocatable :: fc(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: x(:)
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBYSHEV_TEST01'
  write ( *, '(a)' ) '  CHEBYSHEV_COEFFICIENTS computes the coefficients of the'
  write ( *, '(a)' ) '  Chebyshev interpolant.'
  write ( *, '(a)' ) '  CHEBYSHEV_INTERPOLANT evaluates the interpolant.'

  n = 5
  a = -1.0D+00
  b = +1.0D+00

  allocate ( c(1:n) )
  allocate ( fc(1:n) )
  allocate ( x(1:n) )

  call chebyshev_coefficients ( a, b, n, f1, c )

  call chebyshev_zeros ( n, x )
  x(1:n) = 0.5D+00 * ( a + b ) + x(1:n) * 0.5D+00 * ( b - a )

  m = n
  call chebyshev_interpolant ( a, b, n, c, m, x, fc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X) is a trig function:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X           C(I)        F(X)       C(F)(X)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)' ) x(i), c(i), f1( x(i) ), fc(i)
  end do

  deallocate ( c )
  deallocate ( fc )
  deallocate ( x )
!
!  Try a variant interval.
!
  n = 5
  a = 0.0D+00
  b = +3.0D+00

  allocate ( c(1:n) )
  allocate ( fc(1:n) )
  allocate ( x(1:n) )

  call chebyshev_coefficients ( a, b, n, f1, c )

  call chebyshev_zeros ( n, x )
  x(1:n) = 0.5D+00 * ( a + b ) + x(1:n) * 0.5D+00 * ( b - a )

  m = n
  call chebyshev_interpolant ( a, b, n, c, m, x, fc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Consider the same F(X), but now over [0,3]:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X           C(I)        F(X)       C(F)(X)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)' ) x(i), c(i), f1( x(i) ), fc(i)
  end do

  deallocate ( c )
  deallocate ( fc )
  deallocate ( x )
!
!  Try a higher order.
!
  n = 10
  a = -1.0D+00
  b = +1.0D+00

  allocate ( c(1:n) )
  allocate ( fc(1:n) )
  allocate ( x(1:n) )

  call chebyshev_coefficients ( a, b, n, f1, c )

  call chebyshev_zeros ( n, x )
  x(1:n) = 0.5D+00 * ( a + b ) + x(1:n) * 0.5D+00 * ( b - a )

  m = n
  call chebyshev_interpolant ( a, b, n, c, m, x, fc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Consider the same F(X), but now with higher order:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X           C(I)        F(X)       C(F)(X)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)' ) x(i), c(i), f1( x(i) ), fc(i)
  end do

  deallocate ( c )
  deallocate ( fc )
  deallocate ( x )
!
!  Try a polynomial.
!
  n = 10
  a = -1.0D+00
  b = +1.0D+00

  allocate ( c(1:n) )
  allocate ( fc(1:n) )
  allocate ( x(1:n) )

  call chebyshev_coefficients ( a, b, n, f3, c )

  call chebyshev_zeros ( n, x )
  x(1:n) = 0.5D+00 * ( a + b ) + x(1:n) * 0.5D+00 * ( b - a )

  m = n
  call chebyshev_interpolant ( a, b, n, c, m, x, fc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X) is a degree 4 polynomial:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X           C(I)        F(X)       C(F)(X)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)' ) x(i), c(i), f3( x(i) ), fc(i)
  end do

  deallocate ( c )
  deallocate ( fc )
  deallocate ( x )
!
!  Try a function with decaying behavior.
!
  n = 10
  a = -1.0D+00
  b = +1.0D+00

  allocate ( c(1:n) )
  allocate ( fc(1:n) )
  allocate ( x(1:n) )

  call chebyshev_coefficients ( a, b, n, f2, c )

  call chebyshev_zeros ( n, x )
  x(1:n) = 0.5D+00 * ( a + b ) + x(1:n) * 0.5D+00 * ( b - a )

  m = n
  call chebyshev_interpolant ( a, b, n, c, m, x, fc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The polynomial approximation to F(X) decays:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X           C(I)        F(X)       C(F)(X)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)' ) x(i), c(i), f2( x(i) ), fc(i)
  end do

  deallocate ( c )
  deallocate ( fc )
  deallocate ( x )

  return
end
function f1 ( x )

!*****************************************************************************80
!
!! F1 evaluates a function that can be used for Chebyshev interpolation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, a point where the function is to be evaluated.
!
!    Output, real ( kind = 8 ) F1, the function value.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  f1 = cos ( 2.0D+00 * pi * x ) * sin ( 3.0D+00 * pi * x )

  return
end
function f2 ( x )

!*****************************************************************************80
!
!! F2 evaluates a function that can be used for Chebyshev interpolation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input,  real ( kind = 8 ) X, a point where the function is to be evaluated.
!
!    Output, real ( kind = 8 ) F2, the function value.
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) x

  f2 = exp ( x )

  return
end
function f3 ( x )

!*****************************************************************************80
!
!! F3 evaluates a function that can be used for Chebyshev interpolation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, a point where the function is to be evaluated.
!
!    Output, real ( kind = 8 ) F3, the function values.
!
  implicit none

  real ( kind = 8 ) f3
  real ( kind = 8 ) x

  f3 = ( x - 3.0D+00 ) * ( x - 1.0D+00 ) * ( x - 1.0D+00 ) * ( x + 2.0D+00 )

  return
end
