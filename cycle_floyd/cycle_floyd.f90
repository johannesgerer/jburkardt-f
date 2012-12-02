subroutine cycle_floyd ( f, x0, lam, mu )

!*****************************************************************************80
!
!! CYCLE_FLOYD finds a cycle in an iterated mapping using Floyd's method.
!
!  Discussion:
!
!    Suppose we a repeatedly apply a function f(), starting with the argument
!    x0, then f(x0), f(f(x0)) and so on.  Suppose that the range of f is finite.
!    Then eventually the iteration must reach a cycle.  Once the cycle is reached,
!    succeeding values stay within that cycle.
!
!    Starting at x0, there is a "nearest element" of the cycle, which is
!    reached after MU applications of f.
!
!    Once the cycle is entered, the cycle has a length LAM, which is the number
!    of steps required to first return to a given value.
!
!    This function uses Floyd's method to determine the values of MU and LAM,
!    given F and X0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 2, Seminumerical Algorithms,
!    Third Edition,
!    Addison Wesley, 1997,
!    ISBN: 0201896842,
!    LC: QA76.6.K64.
!
!  Parameters:
!
!    Input, external integer ( kind = 4 ) F(), the name of the function 
!    to be analyzed.
!
!    Input, integer ( kind = 4 ) X0, the starting point.
!
!    Output, integer ( kind = 4 ) LAM, the length of the cycle.
!
!    Output, integer ( kind = 4 ) MU, the index in the sequence starting
!    at X0, of the first appearance of an element of the cycle.
!
  implicit none

  external f
  integer ( kind = 4 ) f
  integer ( kind = 4 ) hare
  integer ( kind = 4 ) lam
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) tortoise
  integer ( kind = 4 ) x0

  tortoise = f ( x0 )
  hare = f ( tortoise )

  do while ( tortoise /= hare )
    tortoise = f ( tortoise )
    hare = f ( f ( hare ) )
  end do

  mu = 0
  tortoise = x0

  do while ( tortoise /= hare )
    tortoise = f ( tortoise )
    hare = f ( hare )
    mu = mu + 1
  end do

  lam = 1
  hare = f ( tortoise )
  do while ( tortoise /= hare )
    hare = f ( hare )
    lam = lam + 1
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
