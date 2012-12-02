program main

!*****************************************************************************80
!
!! FD1D_BURGERS_LEAP solves the nonviscous Burgers equation using leapfrogging.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2010
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

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) dt
  real ( kind = 8 ) dx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) n
  integer ( kind = 4 ) step
  integer ( kind = 4 ) step_num
  real ( kind = 8 ) t
  real ( kind = 8 ) t_init
  real ( kind = 8 ) t_last
  real ( kind = 8 ), allocatable :: uc(:)
  real ( kind = 8 ), allocatable :: un(:)
  real ( kind = 8 ), allocatable :: uo(:)
  real ( kind = 8 ), allocatable :: x(:)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_BURGERS_LEAP:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Solve the non-viscous time-dependent Burgers equation,'
  write ( *, '(a)' ) '  using the leap-frog method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Equation to be solved:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    du/dt + u * du/dx = 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  for x in [ a, b ], for t in [t_init, t_last]'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  with initial conditions:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    u(x,o) = u_init'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  and boundary conditions:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    u(a,t) = u_a(t), u(b,t) = u_b(t)'
!
!  Set and report the problem parameters.
!
  n = 21
  a = -1.0D+00
  b = +1.0D+00
  dx = ( b - a ) / real ( n - 1, kind = 8 )
  step_num = 30
  t_init = 0.0D+00
  t_last = 3.0D+00
  dt = ( t_last - t_init ) / real ( step_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6,a,g14.6)' ) a, ' <= X <= ', b
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
  write ( *, '(a,g14.6)' ) '  DX = ', dx
  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6,a,g14.6)' ) t_init, ' <= T <= ', t_last
  write ( *, '(a,i8)' ) '  Number of time steps = ', step_num
  write ( *, '(a,g14.6)' ) '  DT = ', dt

  allocate ( uc(1:n) )
  allocate ( un(1:n) )
  allocate ( uo(1:n) )
  allocate ( x(1:n) )

  call r8vec_even ( n, a, b, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X:'
  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 4, n )
    write ( *, '(2x,5g14.6)' ) x(ilo:ihi)
  end do
!
!  Set the initial condition,
!  and apply boundary conditions to first and last entries.
!
  step = 0
  t = t_init
  call u_init ( n, x, t, un )
  call u_a ( x(1), t, un(1) )
  call u_b ( x(n), t, un(n) )

  call report ( step, step_num, n, x, t, un )
!
!  Use Euler's method to get the first step.
!
  step = 1
  t = ( real ( step_num - step, kind = 8 ) * t_init   &
      + real (            step, kind = 8 ) * t_last ) &
      / real ( step_num,        kind = 8 )

  uc(1:n) = un(1:n)

  un(2:n-1) = uc(2:n-1) - dt * uc(2:n-1) * ( uc(3:n) - uc(1:n-2) ) / 2.0D+00 / dx

  call u_a ( x(1), t, un(1) )
  call u_b ( x(n), t, un(n) )

  call report ( step, step_num, n, x, t, un )
!
!  Subsequent steps use the leapfrog method.
!
  do step = 2, step_num
 
    t = ( real ( step_num - step, kind = 8 ) * t_init   &
        + real (            step, kind = 8 ) * t_last ) &
        / real ( step_num,        kind = 8 )

    uo(1:n) = uc(1:n)
    uc(1:n) = un(1:n)

    un(2:n-1) = uo(2:n-1) - dt * uc(2:n-1) * ( uc(3:n) - uc(1:n-2) ) / dx

    call u_a ( x(1), t, un(1) )
    call u_b ( x(n), t, un(n) )

    call report ( step, step_num, n, x, t, un )

  end do
!
!  Free memory.
!
  deallocate ( uc )
  deallocate ( un )
  deallocate ( uo )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_BURGERS_LEAP:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns an R8VEC of evenly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If N is 1, then the midpoint is returned.
!
!    Otherwise, the two endpoints are returned, and N-2 evenly
!    spaced points between them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) ALO, AHI, the low and high values.
!
!    Output, real ( kind = 8 ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) ahi
  real ( kind = 8 ) alo
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * alo   &
             + real (     i - 1, kind = 8 ) * ahi ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine report ( step, step_num, n, x, t, u )

!*****************************************************************************80
!
!! REPORT prints or plots or saves the data at the current time step.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STEP, the index of the current step,
!    between 0 and STEP_NUM.
!
!    Input, integer ( kind = 4 ) STEP_NUM, the number of steps to take.
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, real ( kind = 8 ) X(N), the coordinates of the nodes.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) U(N), the initial values U(X,T).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) step
  integer ( kind = 4 ) step_num
  real ( kind = 8 ) t
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  STEP = ', step
  write ( *, '(a,g14.6)' ) '  TIME = ', t
  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 4, n )
    write ( *, '(2x,5g14.6)' ) u(ilo:ihi)
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
subroutine u_a ( x, t, ua )

!*****************************************************************************80
!
!! U_A sets the boundary condition for U at A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, T, the position and time.
!
!    Output, real ( kind = 8 ) UA, the prescribed value of U(X,T).
!
  implicit none

  real ( kind = 8 ) t
  real ( kind = 8 ) ua
  real ( kind = 8 ) x

  ua = + 0.5D+00

  return
end
subroutine u_b ( x, t, ub )

!*****************************************************************************80
!
!! U_B sets the boundary condition for U at B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, T, the position and time.
!
!    Output, real ( kind = 8 ) UB, the prescribed value of U(X,T).
!
  implicit none

  real ( kind = 8 ) t
  real ( kind = 8 ) ub
  real ( kind = 8 ) x

  ub = - 0.5D+00

  return
end
subroutine u_init ( n, x, t, u )

!*****************************************************************************80
!
!! U_INIT sets the initial condition for U.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, real ( kind = 8 ) X(N), the coordinates of the nodes.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Output, real ( kind = 8 ) U(N), the initial values U(X,T).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ua
  real ( kind = 8 ) ub
  real ( kind = 8 ) x(n)

  call u_a ( x(1), t, ua )
  call u_b ( x(n), t, ub )

  q = 2.0D+00 * ( ua - ub ) / pi
  r = ( ua + ub ) / 2.0D+00
!
!  S can be varied.  It is the slope of the initial condition at the midpoint.
!
  s = 1.0D+00

  u(1:n) = - q * atan ( s * ( 2.0D+00 * x(1:n) - x(1) - x(n) ) &
    / ( x(n) - x(1) ) ) + r

  return
end

