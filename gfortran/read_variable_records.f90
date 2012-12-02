program main

!*****************************************************************************80
!
!! MAIN is the main program for READ_VARIABLE_RECORDS.
!
!  Discussion:
!
!    READ_VARIABLE_RECORDS demonstrates how to read a variable amount of data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: k_max = 10

  integer i
  integer ios
  integer j
  integer k
  integer values(k_max)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'READ_VARIABLE_RECORDS:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate how to read an unknown number of values.'
  write ( *, '(a)' ) '  Although the number of values is unknown, we will'
  write ( *, '(a)' ) '  assume that the formatting of the values is fixed.'
  write ( *, '(a)' ) '  In particular, we will here assume that each'
  write ( *, '(a)' ) '  record of the file contains integers in the I4 format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  It turns out that the first integer will indicate'
  write ( *, '(a)' ) '  the number of remaining data items on the line,'
  write ( *, '(a)' ) '  but we won''t actually use that fact.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   K   1   2   3   4   5   6   7   8   9  10'
  write ( *, '(a)' ) '  --  --  --  --  --  --  --  --  --  --  --'
  write ( *, '(a)' ) ' '

  open ( unit = 1, file = 'read_variable_records.txt', status = 'old' )

  do
!
!  The "advance = 'no'" means that the next READ statement will pick
!  up exactly where this one left off.
!
    read ( 1, '(i4)', iostat = ios, advance = 'no' ) k

    if ( ios /= 0 ) then
      exit
    end if

    i = 0

    do

      i = i + 1
!
!  The "IOSTAT=IOS" argument will be nonzero when we run out of stuff to read.
!
      read ( 1, '(i4)', advance = 'no', iostat = ios ) j

      if ( ios /= 0 ) then
        exit
      end if

      if ( i < k_max ) then
        values(i) = j
      end if

    end do

    write ( *, '(10i4)' ) k, values(1:min(k,k_max))

  end do

  close ( unit = 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'READ_VARIABLE_RECORDS:'
  write ( *, '(a)' ) '  Normal end of execution.'
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
!    26 February 2005
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
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
