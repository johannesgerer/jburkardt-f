program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS178_PRB.
!
!  Discussion:
!
!    TOMS178_PRB calls sample problems for the TOMS178 library.
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS178_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS178 library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS178_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests HOOKE with the Rosenbrock function.
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nvars = 2

  real    ( kind = 8 ) endpt(nvars)
  real    ( kind = 8 ) eps
  integer ( kind = 4 ) hooke
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itermax
  real    ( kind = 8 ) rho
  real    ( kind = 8 ), external :: rosenbrock
  real    ( kind = 8 ) startpt(nvars)
  real    ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  HOOKE seeks a minimizer of F(X).'
  write ( *, '(a)' ) '  Here we use the Rosenbrock function.'
!
!  Starting guess for Rosenbrock.
!
  startpt(1) = -1.2D+00
  startpt(2) = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Initial estimate X = '
  write ( *, '(a)' ) ' '
  do i = 1, nvars
    write ( *, '(2x,i8,2x,g14.6)' ) i, startpt(i)
  end do

  value = rosenbrock ( startpt, nvars )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X) = ', value
!
!  Call HOOKE.
!
  itermax = 5000
  rho = 0.5D+00
  eps = 1.0D-06

  it = hooke ( nvars, startpt, endpt, rho, eps, itermax, rosenbrock )
!
!  Results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X* = '
  write ( *, '(a)' ) ' '
  do i = 1, nvars
    write ( *, '(2x,i8,2x,g14.6)' ) i, endpt(i)
  end do

  value = rosenbrock ( endpt, nvars )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', value

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests HOOKE with the WOODS function.
!
!  Discussion:
!
!    The Hooke and Jeeves algorithm works well when RHO = 0.5, but
!    does poorly when RHO = 0.6, and better when RHO = 0.8
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nvars = 4

  real    ( kind = 8 ) endpt(nvars)
  real    ( kind = 8 ) eps
  integer ( kind = 4 ) hooke
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itermax
  real    ( kind = 8 ) rho
  real    ( kind = 8 ) startpt(nvars)
  real    ( kind = 8 ) value
  real    ( kind = 8 ), external :: woods

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  HOOKE seeks a minimizer of F(X).'
  write ( *, '(a)' ) '  Here we use the Woods function.'
!
!  Starting guess.
!
  startpt(1) = -3.0D+00
  startpt(2) = -1.0D+00
  startpt(3) = -3.0D+00
  startpt(4) = -1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Initial estimate X = '
  write ( *, '(a)' ) ' '
  do i = 1, nvars
    write ( *, '(2x,i8,2x,g14.6)' ) i, startpt(i)
  end do

  value = woods ( startpt, nvars )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X) = ', value
!
!  Call HOOKE.
!
  itermax = 5000
  rho = 0.5D+00
  eps = 1.0D-06

  it = hooke ( nvars, startpt, endpt, rho, eps, itermax, woods )
!
!  Results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X* = '
  write ( *, '(a)' ) ' '
  do i = 1, nvars
    write ( *, '(2x,i8,2x,g14.6)' ) i, endpt(i)
  end do

  value = woods ( endpt, nvars )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', value

  return
end
function rosenbrock ( x, n )

!*****************************************************************************80
!
!! ROSENBROCK evaluates the Rosenbrock function.
!
!  Discussion:
!
!    The Hooke and Jeeves algorithm works reasonably well on
!    Rosenbrock's test function, depending on the value of RHO chosen.
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the function.
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, real ( kind = 8 ) ROSENBROCK, the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) rosenbrock
  real    ( kind = 8 ) x(n)

  rosenbrock = 100.0 * ( x(2) - x(1) * x(1) )**2 &
             +         ( 1.0D+00 - x(1) )**2

  return
end
function woods ( x, n )

!*****************************************************************************80
!
!! WOODS evaluates the Woods function.
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the function.
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, real ( kind = 8 ) WOODS, the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) s1
  real    ( kind = 8 ) s2
  real    ( kind = 8 ) s3
  real    ( kind = 8 ) t1
  real    ( kind = 8 ) t2
  real    ( kind = 8 ) t3
  real    ( kind = 8 ) t4
  real    ( kind = 8 ) t5
  real    ( kind = 8 ) woods
  real    ( kind = 8 ) x(n)

  s1 = x(2) - x(1) * x(1)
  s2 = 1.0D+00 - x(1)
  s3 = x(2) - 1.0D+00
  t1 = x(4) - x(3) * x(3)
  t2 = 1.0D+00 - x(3)
  t3 = x(4) - 1.0D+00
  t4 = s3 + t3
  t5 = s3 - t3

  woods = 100.0D+00 * s1**2 &
        +             s2**2 &
        +  90.0D+00 * t1**2 &
        +             t2**2 &
        +  10.0D+00 * t4**2 &
        +   0.1D+00 * t5**2

  return
end

