program main

!*****************************************************************************80
!
!! MAIN is the main program for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real      ( kind = 8 ) a
  real      ( kind = 8 ) b
  real      ( kind = 8 ), external :: f1
  logical                header
  real      ( kind = 8 ), external :: k1
  real      ( kind = 8 ) u(n)
  character ( len = 80 ) u_file
  real      ( kind = 8 ) ua
  real      ( kind = 8 ) ub
  real      ( kind = 8 ) x(n)
  character ( len = 80 ) x_file

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROBLEM1:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A test problem for FD1D_HEAT_STEADY.'

  a = 0.0D+00
  b = 1.0D+00

  ua = 0.0D+00
  ub = 1.0D+00

  call fd1d_heat_steady ( n, a, b, ua, ub, k1, f1, x, u )

  x_file = 'problem1_nodes.txt'
  call r8mat_write ( x_file, 1, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X data written to "' // trim ( x_file ) // '".'

  u_file = 'problem1_values.txt'
  call r8mat_write ( u_file, 1, n, u )

  write ( *, '(a)' ) '  U data written to "' // trim ( u_file ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROBLEM1:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function k1 ( x )

!*****************************************************************************80
!
!! K1 evaluates the heat transfer coefficient K(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the position.
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
!! F1 evaluates the right hand side of the steady state heat equation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the position.
!
!    Output, real ( kind = 8 ) F1, the value of F(X).
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  f1 = 0.0D+00

  return
end


