program main

!*****************************************************************************80
!
!! RK4_PRB demonstrates the use of the RK4 one-step ODE solver.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Local, real ( kind = 8 ) DT, the time step.
!
!    Local, real ( kind = 8 ) T0, the time at which the solution is known.
!
!    Local, real ( kind = 8 ) TMAX, the maximum time at which a solution is desired.
!
!    Local, real ( kind = 8 ) U0, the estimated solution at time T0.
!
  implicit none

  real ( kind = 8 ), parameter :: dt = 0.1D+00
  external f3
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  real ( kind = 8 ), parameter :: tmax = 12.0D+00 * pi
  real ( kind = 8 ) u0
  real ( kind = 8 ) u1

  t0 = 0.0D+00
  u0 = 0.5D+00

  do
!
!  Print (T0,U0).
!
    write ( *, '(2x,g14.6,2x,g14.6)' ) t0, u0
!
!  Stop if we've exceeded TMAX.
!
    if ( tmax <= t0 ) then
      exit
    end if
!
!  Otherwise, advance to time T1, and have RK4 estimate 
!  the solution U1 there.
!
    t1 = t0 + dt
    call rk4 ( t0, u0, dt, f3, u1 )
!
!  Shift the data to prepare for another step.
!
    t0 = t1
    u0 = u1

  end do

  stop
end
subroutine f3 ( t, u, uprime )

!*****************************************************************************80
!
!! F3 evaluates the right hand side of a particular ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) U, the current solution value.
!
!    Output, real ( kind = 8 ) UPRIME, the value of the derivative, dU/dT.
!
  implicit none

  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) uprime
  
  uprime = u * cos ( t )
 
  return
end
