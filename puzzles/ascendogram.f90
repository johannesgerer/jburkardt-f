program main

!*****************************************************************************80
!
!! MAIN is the main program for ASCENDOGRAM.
!
!  Discussion:
!
!    ASCENDOGRAM searches for "ascendograms".
!
!    An ascendogram is a word in which the letters are in alphabetical order.
!    In a "loose" ascendogram, letters may be repeated, as in "abbot".
!    In a "strict" ascendogram, no letter may be repeated, as in "sty".
!
!    This program reads words from a file, and writes each loose or strict
!    ascendogram to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer ios
  integer, parameter :: length_min = 5
  logical loose
  character ( len = 80 ) :: loose_file = 'loose.txt'
  integer loose_number
  integer, parameter :: loose_unit = 2
  integer n
  logical strict
  character ( len = 80 ) :: strict_file = 'strict.txt'
  integer strict_number
  integer, parameter :: strict_unit = 3
  character ( len = 80 ) string
  character ( len = 80 ) :: word_file = 'wordlist.txt'
  integer word_long_number
  integer word_number
  integer, parameter :: word_unit = 1

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASCENDOGRAM:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine a file "' // trim ( word_file ) // '"'
  write ( *, '(a)' ) '  full of words, and look for "ascendograms" '
  write ( *, '(a)' ) '  (words whose letters are in alphabetical order).'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  We only consider words of length at least ', length_min
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Copy loose ascendograms to "' // trim ( loose_file ) // '".'
  write ( *, '(a)' ) &
    '  Copy strict ascendograms to "' // trim ( strict_file ) // '".'

  word_long_number = 0
  word_number = 0
  loose_number = 0
  strict_number = 0

  open ( unit = word_unit, file = word_file, status = 'old' )
  open ( unit = loose_unit, file = loose_file, status = 'replace' )
  open ( unit = strict_unit, file = strict_file, status = 'replace' )
!
!  Read a word from the word file, and judge it.
!
  do

    read ( word_unit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then 
      exit
    end if

    word_number = word_number + 1

    call s_low ( string )

    n = len_trim ( string )

    if ( length_min <= n ) then

      word_long_number = word_long_number + 1

      loose = .true.
      strict = .true.

      do i = 1, n-1

        if ( lgt ( string(i:i), string(i+1:i+1) ) ) then
          loose = .false.
        end if

        if ( lge ( string(i:i), string(i+1:i+1) ) ) then
          strict = .false.
        end if

      end do

      if ( loose ) then
        write ( loose_unit, '(a)' ) trim ( string )
        loose_number = loose_number + 1
      end if

      if ( strict ) then
        write ( strict_unit, '(a)' ) trim ( string )
        strict_number = strict_number + 1
      end if

    end if

  end do

  close ( unit = word_unit )
  close ( unit = loose_unit )
  close ( unit = strict_unit )

  write ( *, * ) ' '
  write ( *, * ) '  Number of words checked                  ', word_number
  write ( *, * ) '  Number of words that were "long enough"  ', word_long_number
  write ( *, * ) '  Number of loose ascendogram words found  ', loose_number
  write ( *, * ) '  Number of strict ascendogram words found ', strict_number

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASCENDOGRAM:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine s_low ( s )

!*****************************************************************************80
!
!! S_LOW replaces all uppercase letters by lowercase ones.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be
!    transformed.  On output, the string is all lowercase.
!
  implicit none

  integer i
  integer nchar
  character ( len = * ) s

  nchar = len_trim ( s )

  do i = 1, nchar
    call c_low ( s(i:i) )
  end do

  return
end
subroutine c_low ( c )

!*****************************************************************************80
!
!! C_LOW lowercases a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to be lowercased.
!
  implicit none

  character c
  integer i

  i = ichar ( c )

  if ( 65 <= i .and. i <= 90 ) then
    c = char ( i + 32 )
  end if

  return
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
