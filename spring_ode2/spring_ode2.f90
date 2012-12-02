program main

!*****************************************************************************80
!
!! MAIN is the main program for SPRING_ODE2.
!
!  Discussion:
!
!    This is a revision of the SPRING_ODE code.
!
!    In this revision of the program, we want to use vectors (C arrays) to 
!    store the data, and we want to write the data out to a file in a form 
!    that Gnuplot (or other plotting programs) can use.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2012
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

  integer ( kind = 4 ), parameter :: n = 100

  real ( kind = 8 ) c
  real ( kind = 8 ) dt
  integer ( kind = 4 ) i
  real ( kind = 8 ) k
  real ( kind = 8 ) m
  real ( kind = 8 ) t(0:n)
  real ( kind = 8 ) t_final
  real ( kind = 8 ) v(0:n)
  real ( kind = 8 ) x(0:n)

  call timestamp ( )
  write ( *, '(a)' ) '#'
  write ( *, '(a)' ) '#SPRING_ODE2'
  write ( *, '(a)' ) '#  FORTRAN90 version'
  write ( *, '(a)' ) '#  Approximate the solution of a spring equation.'
  write ( *, '(a)' ) '#  Write data to a file for use by gnuplot.'
  write ( *, '(a)' ) '#'
!
!  Data
!
  m = 1.0D+00
  k = 1.0D+00
  c = 0.3D+00
  t_final = 20.0D+00
  dt = t_final / real ( n, kind = 8 )
!
!  Initial conditions.
!
  t(0) = 0.0D+00
  x(0) = 1.0D+00
  v(0) = 0.0D+00
!
!  Compute the approximate solution at equally spaced times.
!
  do i = 1, n

    t(i) = real ( i, kind = 8 ) * t_final / real ( n, kind = 8 )
    x(i) = x(i-1) + dt * v(i-1)
    v(i) = v(i-1) + ( dt / m ) * ( - k * x(i-1) - c * v(i-1) )

  end do
!
!  Write the data to a file for plotting, possibly by gnuplot.
!  gnuplot expects T, X, and V to be columns of output.
!
  do i = 0, n
    write ( *, '(g14.6,2x,g14.6,2x,g14.6)' ) t(i), x(i), v(i)
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) '#'
  write ( *, '(a)' ) '#SPRING_ODE2:'
  write ( *, '(a)' ) '#  Normal end of execution.'
  write ( *, '(a)' ) '#'
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

  write ( *, '(a1,1x,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    '#', d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
