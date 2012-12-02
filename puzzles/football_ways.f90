program main

!*****************************************************************************80
!
!! MAIN is the main program for FOOTBALL_WAYS.
!
!  Discussion:
!
!    FOOTBALL_WAYS counts the ways of getting a football score.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 June 2005
!
!  Author:
!
!    John Burkardt
!
  integer, parameter :: n_max = 100

  integer n
  integer ( kind = 8 ) s(0:n_max)

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FOOTBALL_WAYS'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program determines the number of ways'
  write ( *, '(a)' ) '  of achieving a particular score in football.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We assume scoring options are:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  +2 for a safety;'
  write ( *, '(a)' ) '  +3 for a field goal;'
  write ( *, '(a)' ) '  +6 for a touchdown with no followup;'
  write ( *, '(a)' ) '  +7 for a touchdown with a point bonus;'
  write ( *, '(a)' ) '  +8 for a touchdown with two point conversion;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If the other team forfeits, a token score of 1'
  write ( *, '(a)' ) '  is assigned.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N                  S(N)'
  write ( *, '(a)' ) ' '

  s(0) = 1
  s(1) = 1
  s(2) = 1
  s(3) = 1
  s(4) = 1
  s(5) = 2
  s(6) = 3
  s(7) = 4
  s(8) = 7

  do n = 0, 8
    write ( *, '(2x,i6,2x,i20)' ) n, s(n)
  end do

  do n = 9, n_max
    s(n) = s(n-2) + s(n-3) + s(n-6) + s(n-7) + s(n-8)
    write ( *, '(2x,i6,2x,i20)' ) n, s(n)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FOOTBALL_WAYS'
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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2001
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
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

