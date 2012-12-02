program main

!*****************************************************************************80
!
!! MAIN is the main program for SPRING_ODE.
!
!  Discussion:
!
!    This is a simple example of how to plot when you don't have a plotter.
!    This is a particular kind of "ASCII graphics", or "typewriter graphics"
!    or "lineprinter graphics", and shows you how valuable an illustration 
!    can be, even when it's as crude as this example.
!
!    Hooke's law for a spring observes that the restoring force is
!    proportional to the displacement: F = - k x
!
!    Newton's law relates the force to acceleration: F = m a
!
!    Putting these together, we have
!
!      m * d^2 x/dt^2 = - k * x
!
!    We can add a damping force with coefficient c:
!
!      m * d^2 x/dt^2 = - k * x - c * dx/dt
!
!    If we write this as a pair of first order equations for (x,v), we have
!
!          dx/dt = v
!      m * dv/dt = - k * x - c * v
!
!    and now we can approximate these values for small time steps.
!
!    Note that the plotting assumes that the value of X will always be
!    between -1 and +1.  If the initial condition uses V = 0, and X starts
!    between -1 and +1, then this will be OK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2012
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

  real ( kind = 8 ) c
  real ( kind = 8 ) dt
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) k
  real ( kind = 8 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  real ( kind = 8 ) t
  real ( kind = 8 ) t_final
  real ( kind = 8 ) v
  real ( kind = 8 ) v_old
  real ( kind = 8 ) x
  real ( kind = 8 ) x_old
  character z(21)

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPRING_ODE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Approximate the solution of a spring equation.'
  write ( *, '(a)' ) '  Display the solution with line printer graphics.'
  write ( *, '(a)' ) ' '
!
!  Data
!
  m = 1.0
  k = 1.0
  c = 0.3
  t_final = 20.0
  n = 100
  dt = t_final / real ( n )
!
!  Initial conditions.
!
  x = 1.0
  v = 0.0
!
!  Compute the approximate solution at equally spaced times.
!
  do i = 0, n

    x_old = x
    v_old = v

    t = real ( i ) * t_final / real ( n )
    x = x_old + dt * v_old
    v = v_old + ( dt / m ) * ( - k * x_old - c * v_old )
!
!  Approximate the position of X in [-1,+1] to within 1/10.
!
    p = int ( ( (  1 * ( 1.0 - x ) ) + 21 * ( 1.0 + x ) ) / 2.0 ) 
    p = max ( p, 1 )
    p = min ( p, 21 )
!
!  Fill in the next line of the plot, placing 'x' in the p position.
!
    if ( mod ( i, 10 ) == 0 ) then
      z(1:21) = '-'
    else
      z(1:21) = ' '
    end if
    z(1) = '|'
    z(6) = '.'
    z(11) = '+'
    z(16) = '.'
    z(21) = '|'

    z(p) = 'x'

    write ( *, '(21a)' ) z(1:21)

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPRING_ODE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
