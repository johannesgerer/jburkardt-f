subroutine llsq ( n, x, y, a, b )

!*****************************************************************************80
!
!! LLSQ solves a linear least squares problem matching a line to data.
!
!  Discussion:
!
!    A formula for a line of the form Y = A * X + B is sought, which
!    will minimize the root-mean-square error to N data points ( X(I), Y(I) );
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the data points.
!
!    Output, real ( kind = 8 ) A, B, the slope and Y-intercept of the 
!    least-squares approximant to the data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bot
  real ( kind = 8 ) top
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xbar
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ybar
!
!  Special case.
!
  if ( n == 1 ) then
    a = 0.0D+00
    b = y(1)
    return
  end if
!
!  Average X and Y.
!
  xbar = sum ( x(1:n) ) / real ( n, kind = 8 )
  ybar = sum ( y(1:n) ) / real ( n, kind = 8 )
!
!  Compute Beta.
!
  top = dot_product ( x(1:n) - xbar, y(1:n) - ybar )
  bot = dot_product ( x(1:n) - xbar, x(1:n) - xbar )

  a = top / bot

  b = ybar - a * xbar

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
