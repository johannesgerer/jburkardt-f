program main

!*****************************************************************************80
!
!! MAIN is the main program for problem 4.
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

  integer ( kind = 4 ), parameter :: n = 21

  real      ( kind = 8 ) a
  real      ( kind = 8 ) b
  real      ( kind = 8 ), external :: f4
  logical                header
  real      ( kind = 8 ), external :: k4
  real      ( kind = 8 ) u(n)
  character ( len = 80 ) u_file
  real      ( kind = 8 ) ua
  real      ( kind = 8 ) ub
  real      ( kind = 8 ) x(n)
  character ( len = 80 ) x_file

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROBLEM4:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A test problem for FD1D_HEAT_STEADY.'
  write ( *, '(a)' ) '  A heat source and a heat sink.'

  a = 0.0D+00
  b = 1.0D+00

  ua = 0.0D+00
  ub = 0.0D+00

  call fd1d_heat_steady ( n, a, b, ua, ub, k4, f4, x, u )

  x_file = 'problem4_nodes.txt'
  call r8mat_write ( x_file, 1, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X data written to "' // trim ( x_file ) // '".'

  u_file = 'problem4_values.txt'
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
function k4 ( x )

!*****************************************************************************80
!
!! K4 evaluates the heat transfer coefficient K(X).
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
!    Output, real ( kind = 8 ) K4, the value of K(X).
!
  implicit none

  real ( kind = 8 ) k4
  real ( kind = 8 ) x

  k4 = 1.0D+00

  return
end
function f4 ( x )

!*****************************************************************************80
!
!! F4 evaluates the right hand side of the steady state heat equation.
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
!    Output, real ( kind = 8 ) F4, the value of F(X).
!
  implicit none

  real ( kind = 8 ) f4
  real ( kind = 8 ) x

  if ( x < 0.15D+00 ) then
    f4 = 0.0D+00
  else if ( x < 0.35D+00 ) then
    f4 = 1.0D+00
  else if ( x < 0.75D+00 ) then
    f4 = 0.0D+00
  else if ( x < 0.85D+00 ) then
    f4 = - 2.0D+00
  else
    f4 = 0.0D+00
  end if

  return
end


