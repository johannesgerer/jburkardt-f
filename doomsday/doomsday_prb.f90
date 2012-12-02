program main

!*****************************************************************************80
!
!! DOOMSDAY_PRB tests DOOMSDAY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DOOMSDAY_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the DOOMSDAY library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DOOMSDAY_PRB:'
  write ( *, '(a)' ) '  Test the DOOMSDAY library.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DOOMSDAY against a couple of test dates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) w
  character ( len = 10 ) s1
  character ( len = 10 ) s2
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Try a couple selected dates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YYYY  MM  DD  Weekday    Weekday'
  write ( *, '(a)' ) '                Tabulated  Computed'
  write ( *, '(a)' ) ' '

  y = 1989
  m = 7
  d = 13
  call doomsday_gregorian ( y, m, d, w )
  call weekday_to_name_common ( w, s1 )
  s2 = 'Thursday'
  write ( *, '(2x,i4,2x,i2,2x,i2,2x,a10,2x,a10)' ) y, m, d, s1, s2

  y = 2012
  m = 5
  d = 26
  call doomsday_gregorian ( y, m, d, w )
  call weekday_to_name_common ( w, s1 )
  s2 = 'Saturday'
  write ( *, '(2x,i4,2x,i2,2x,i2,2x,a10,2x,a10)' ) y, m, d, s1, s2

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests DOOMSDAY against a number of known values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) w1
  integer ( kind = 4 ) w2
  character ( len = 10 ) s1
  character ( len = 10 ) s2
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  WEEKDAY_VALUES supplies a list of dates and weekdays.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YYYY  MM  DD  Weekday    Weekday'
  write ( *, '(a)' ) '                Tabulated  Computed'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call weekday_values ( n_data, y, m, d, w1 )

    if ( n_data <= 0 ) then
      exit
    end if
!
!  The transition from Julian to Gregorian calendars occurred in 1582
!  (for some people).  The data in "WEEKDAY_VALUES" before the transition
!  is stored in Julian format, which DOOMSDAY_GREGORIAN can't handle.
!  So let's just refuse to handle 1582 or earlier!
!
    if ( y <= 1582 ) then
      cycle
    end if

    call doomsday_gregorian ( y, m, d, w2 )

    call weekday_to_name_common ( w1, s1 )
    call weekday_to_name_common ( w2, s2 )

    write ( *, '(2x,i4,2x,i2,2x,i2,2x,a10,2x,a10)' ) y, m, d, s1, s2

  end do

  return
end
