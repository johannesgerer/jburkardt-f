program main

!*****************************************************************************80
!
!! MAIN is the main program for QUAD2D_SERIAL.
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
  double precision a
  double precision b
  double precision error
  double precision exact
  external f
  double precision f
  integer i
  integer j
  integer n
  integer nx
  integer ny
  double precision pi
  double precision total
  double precision wtime
  double precision wtime1
  double precision wtime2
  double precision x
  double precision y

  a = 0.0
  b = 1.0
  nx = 32768
  ny = 32768
  n = nx * ny
  pi = 3.141592653589793D+00
  exact = pi * pi / 6.0D+00

  call timestamp ( )
  write ( *, * ) ' '
  write ( *, * ) 'QUAD2D_SERIAL:'
  write ( *, * ) '  FORTRAN90 version'
  write ( *, * ) '  Estimate the integral of f(x,y) over [0,1]x[0,1].'
  write ( *, * ) '  f(x,y) = 1 / ( 1 - x * y ).'
  write ( *, * ) ' '
  write ( *, * ) '  A        = ', a
  write ( *, * ) '  B        = ', b
  write ( *, * ) '  NX       = ', nx
  write ( *, * ) '  NY       = ', ny
  write ( *, * ) '  N        = ', n
  write ( *, * ) '  Exact    = ', exact

  call cpu_time ( wtime1 )

  total = 0.0D+00
  do i = 1, nx
    x = ( ( 2 * nx - 2 * i + 1 ) * a + ( 2 * i - 1 ) * b ) / ( 2 * nx )
    do j = 1, ny
      y = ( ( 2 * ny - 2 * j + 1 ) * a + ( 2 * j - 1 ) * b ) / ( 2 * ny )
      total = total + f ( x, y )
    end do
  end do

  call cpu_time ( wtime2 )

  total = ( b - a ) * ( b - a ) * total / dble ( nx ) / dble ( ny )
  error = abs ( total - exact )
  wtime = wtime2 - wtime1
 
  write ( *, * ) ' '
  write ( *, * ) '  Estimate = ', total
  write ( *, * ) '  Error    = ', error
  write ( *, * ) '  Time     = ', wtime
!
!  Termiante.
!
  write ( *, * ) ' '
  write ( *, * ) 'QUAD2D_SERIAL:'
  write ( *, * ) '  Normal end of execution.'
  write ( *, * ) ' '
  call timestamp ( )

  stop
end
function f ( x, y )

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
!    Input, double precision X, Y, the coordinates of a point.
!
!    Output, double precision F, the function value at (X,Y).
!
  double precision f
  double precision x
  double precision y

  f = 1.0D+00 / ( 1.0D+00 - x * y )

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
