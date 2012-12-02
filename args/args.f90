program main

!*****************************************************************************80
!
!! MAIN is the main program for ARGS.
!
!  Discussion:
!
!    ARGS demonstrates the use of the (semi-standard) GETARGS utility.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 September 2002
!
!  Author:
!
!    John Burkardt
!
!  Usage:
!
!    args arg1 arg2 arg3 ...
!
  implicit none

  character ( len = 80 ) arg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) numarg

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ARGS'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate the use of command line argument'
  write ( *, '(a)' ) '  routines in a FORTRAN program.'
  write ( *, '(a)' ) '  These include GETARG and IARGC.'

  numarg = iargc ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) &
    '  ARGS was called with IARGC() = ', numarg, ' arguments.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CALL GETARG(I,ARG) returns the arguments:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I     ARG '
  write ( *, '(a)' ) ' '

  do i = 0, numarg
    call getarg ( i, arg )
    write ( *, '(2x,i3,2x,a20)' ) i, arg
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ARGS:'
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
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer values(8)
  integer ( kind = 4 ) y
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
