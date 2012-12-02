program main

!*****************************************************************************80
!
!! MAIN is the main program for QUAD_SERIAL.
!
!  Modified:
!
!    24 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  double precision a
  double precision b
  double precision error
  double precision exact
  external f
  double precision f
  integer i
  integer n
  double precision total
  double precision wtime
  double precision wtime1
  double precision wtime2
  double precision x

  a =  0.0
  b = 10.0
  n = 10000000
  exact = 0.49936338107645674464D+00

  call timestamp ( )
  write ( *, * ) ' '
  write ( *, * ) 'QUAD_SERIAL:'
  write ( *, * ) '  FORTRAN90 version'
  write ( *, * ) '  Estimate the integral of f(x) from A to B.'
  write ( *, * ) '  f(x) = 50 / ( pi * ( 2500 * x * x + 1 ) ).'
  write ( *, * ) ' '
  write ( *, * ) '  A        = ', a
  write ( *, * ) '  B        = ', b
  write ( *, * ) '  N        = ', n
  write ( *, * ) '  Exact    = ', exact

  call cpu_time ( wtime1 )

  total = 0.0D+00
  do i = 1, n
    x = ( ( n - i ) * a + ( i - 1 ) * b ) / ( n - 1 )
    total = total + f ( x )
  end do

  call cpu_time ( wtime2 )

  total = ( b - a ) * total / dble ( n )
  error = abs ( total - exact )
  wtime = wtime2 - wtime1
 
  write ( *, * ) ' '
  write ( *, * ) '  Estimate = ', total
  write ( *, * ) '  Error    = ', error
  write ( *, * ) '  Time     = ', wtime
!
!  Terminate.
!
  write ( *, * ) ' '
  write ( *, * ) 'QUAD_SERIAL:'
  write ( *, * ) '  Normal end of execution.'
  write ( *, * ) ' '
  call timestamp ( )

  stop
end
function f ( x )

!*****************************************************************************80
!
!! F evaluates the function.
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
!    Input, double precision X, the evaluation point.
!
!    Output, double precision F, the function value.
!
  double precision f
  double precision pi
  double precision x

  pi = 3.141592653589793D+00
  f = 50.0D+00 / ( pi * ( 2500.0D+00 * x * x + 1.0D+00 ) )

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

  character ( len = 8  ) ampm
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
