program main

!*****************************************************************************80
!
!! FEM1D_HEAT_STEADY_PRB tests the routines in FEM1D_HEAT_STEADY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_HEAT_STEADY_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FEM1D_HEAT_STEADY library.'

  call fem1d_heat_steady_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_HEAT_STEADY_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine fem1d_heat_steady_test01 ( )

!*****************************************************************************80
!
!! FEM1D_HEAT_STEADY_TEST01 carries out test case #1.
!
!  Discussion:
!
!    Use K1, F1, EXACT1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) exact1
  real ( kind = 8 ), external :: f1
  integer ( kind = 4 ) i
  real ( kind = 8 ), external :: k1
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ua
  real ( kind = 8 ) ub
  real ( kind = 8 ) uexact
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_HEAT_STEADY_TEST01'
  write ( *, '(a)' ) '  K1(X)  = 1.0'
  write ( *, '(a)' ) '  F1(X)  = X * ( X + 3 ) * exp ( X )'
  write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
!
!  Geometry definitions.
!
  a = 0.0D+00
  b = 1.0D+00
  ua = 0.0D+00
  ub = 0.0D+00
  call r8vec_even ( n, a, b, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
  write ( *, '(a,f10.4)' ) '  Left endpoint A = ', a
  write ( *, '(a,f10.4)' ) '  Right endpoint B = ', b
  write ( *, '(a,f10.4)' ) '  Prescribed U(A) = ', ua
  write ( *, '(a,f10.4)' ) '  Prescribed U(B) = ', ub

  call fem1d_heat_steady ( n, a, b, ua, ub, k1, f1, x, u )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    X         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    uexact = exact1 ( x(i) )
    write ( *, '(2x,i4,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, x(i), u(i), uexact, abs ( u(i) - uexact )
  end do

  return
end
function k1 ( x )

!*****************************************************************************80
!
!! K1 evaluates K function #1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) K1, the value of K(X).
!
  implicit none

  real ( kind = 8 ) k1
  real ( kind = 8 ) x

  k1 = 1.0D+00

  return
end
function f1 ( x )

!*****************************************************************************80
!
!! F1 evaluates right hand side function #1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) F1, the value of F(X).
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  f1 = x * ( x + 3.0D+00 ) * exp ( x )

  return
end
function exact1 ( x )

!*****************************************************************************80
!
!! EXACT1 evaluates exact solution #1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) EXACT1, the value of U(X).
!
  implicit none

  real ( kind = 8 ) exact1
  real ( kind = 8 ) x

  exact1 = x * ( 1.0D+00 - x ) * exp ( x )

  return
end

