program main

!*****************************************************************************80
!
!! MAIN is the main program for STOCHASTIC_RK_PRB.
!
!  Discussion:
!
!    STOCHASTIC_RK_PRB calls a set of problems for STOCHASTIC_RK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STOCHASTIC_RK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the STOCHASTIC_RK library.'
 
  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STOCHASTIC_RK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
   
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests RK1_TI_STEP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ), external :: fi
  real    ( kind = 8 ), external :: gi
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  real    ( kind = 8 ) q
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) t
  real    ( kind = 8 ), parameter :: t0 = 0.0D+00
  real    ( kind = 8 ), parameter :: tn = 1.0D+00
  real    ( kind = 8 ) x(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  RK1_TI_STEP uses a first order RK method'
  write ( *, '(a)' ) '  for a problem whose right hand side does not'
  write ( *, '(a)' ) '  depend explicitly on time.'

  h = ( tn - t0 ) / real ( n, kind = 8 )
  q = 1.0D+00
  seed = 123456789

  i = 0
  t = t0
  x(i) = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I           T             X'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8,2x,f14.6,2x,g14.6)' ) i, t, x(i)

  do i = 1, n

    t = ( real ( n - i, kind = 8 ) * t0   &
        + real (     i, kind = 8 ) * tn ) &
        / real ( n,     kind = 8 )

    call rk1_ti_step ( x(i-1), t, h, q, fi, gi, seed, x(i) )

    write ( *, '(2x,i8,2x,f14.6,2x,g14.6)' ) i, t, x(i)

  end do

  return
end
function fi ( x )

!*****************************************************************************80
!
!! FI is a time invariant deterministic right hand side.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) FI, the value.
!
  implicit none

  real    ( kind = 8 ) fi
  real    ( kind = 8 ) x

  fi = 1.0D+00

  return
end
function gi ( x )

!*****************************************************************************80
!
!! GI is a time invariant stochastic right hand side.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) GI, the value.
!
  implicit none

  real    ( kind = 8 ) gi
  real    ( kind = 8 ) x

  gi = 1.0D+00

  return
end
