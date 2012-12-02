program main

!*****************************************************************************80
!
!! MAIN is the main program for F77_CLEANUP.
!
!  Discussion:
!
!    F77_CLEANUP makes a copy of a FORTRAN77 file in which:
!
!    * TAB characters are replaced by six spaces;
!    * all comments begin with a lowercase "c";
!    * all continuation characters are "&";
!    * all (noncomment) lines are no longer than 72 characters.
!
!  Usage:
!
!    f77_cleanup input.f77 output.f77
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  character ( len = 100 ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) num_arg
  character ( len = 100 ) output_file
  integer ( kind = 4 ) output_unit

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'F77_CLEANUP'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read a FORTRAN77 file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Make a copy in which:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' TAB characters are replaced by six spaces;'
  write ( *, '(a)' ) ' * comments begin with a lowercase "c";'
  write ( *, '(a)' ) ' * continuation characters are "&";'
  write ( *, '(a)' ) ' * (noncomment) lines are no longer than 72 characters.'
!
!  Get the number of command line arguments.
!
  num_arg = iargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( num_arg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the input file name:'
    read ( *, '(a)', iostat = ios ) input_file

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'F77_CLEANUP - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1
    call getarg ( iarg, input_file )

  end if
!
!  If two command line arguments, the second one is the output file name.
!
  if ( num_arg < 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the output file name:'
    read ( *, '(a)', iostat = ios ) output_file

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'F77_CLEANUP - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 2
    call getarg ( iarg, output_file )

  end if
!
!  Now we know what to do.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'F77_CLEANUP'
  write ( *, '(a)' ) '  Read F77 file:  "' // trim ( input_file ) // '".'
  write ( *, '(a)' ) '  Write F77 file: "' // trim ( output_file ) // '".'

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F77_CLEANUP - Fatal error!'
    write ( *, '(a)' ) '  There was an error while trying to open the'
    write ( *, '(a)' ) '  input file "' // trim ( input_file ) // '".'
    stop
  end if

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace' )

  call cleanup_f77 ( input_unit, output_unit )

  close ( unit = input_unit )
  close ( unit = output_unit )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'F77_CLEANUP'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine cleanup_f77 ( input, output )

!*****************************************************************************80
!
!! CLEANUP_F77 copies an F77 file, using "&" for continuation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT, OUTPUT, the I/O units associated
!    with the input and output files respectively.
!
  implicit none

  logical, save :: found_tab = .false.
  logical have_line
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  integer ( kind = 4 ) output
  character, parameter :: TAB = char ( 9 )

  have_line = .false.

  do
!
!  Read the next line.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  If the line contains any TAB characters, then
!
    if ( index ( line, TAB ) /= 0 ) then

      if ( .not. found_tab ) then
        found_tab = .true.
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'F77_TO_F77_CONTINUATION - Warning!'
        write ( *, '(a)' ) '  This file contains loathsome TAB characters.'
        write ( *, '(a)' ) '  They will be replaced by 6 blanks.'
        write ( *, '(a)' ) '  The results may not be what you want, but they'
        write ( *, '(a)' ) '  are better than you deserve.'
      end if

      call s_tab_blanks ( line )

    end if
!
!  If it's an F77 comment line, write it out immediately, with a "c".
!
    if ( &
      line(1:1) == '!' .or. &
      line(1:1) == '*' .or. &
      line(1:1) == 'c' .or. &
      line(1:1) == 'C' .or. &
      line(1:1) == 'd' .or. &
      line(1:1) == 'D' ) then

      line(1:1) = 'c'
      write ( output, '(a)' ) trim ( line )
      cycle

    end if
!
!  Now is the time to truncate the input line to 72 columns.
!
    line = line(1:72)
!
!  If the new line is a continuation of the old line, then...
!
    if ( line(6:6) /= ' ' ) then
      line(6:6) = '&'
    end if

    write ( output, '(a)' ) trim ( line )

  end do

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine s_tab_blanks ( s )

!*****************************************************************************80
!
!! S_TAB_BLANKS replaces TAB characters by 6 spaces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be modified.  On
!    output, some significant characters at the end of S may have
!    been lost.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) ntab
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
!
!  If no TAB's occur in the line, there is nothing to do.
!
  if ( index ( s, TAB ) == 0 ) then
    return
  end if
!
!  Otherwise, find out how long the string is.
!
  lenc = len_trim ( s )
  lens = len ( s )
!
!  Count the TAB's.
!
  ntab = 0
  do i = 1, lenc
    if ( s(i:i) == TAB ) then
      ntab = ntab + 1
    end if
  end do
!
!  Now copy the string onto itself, going backwards.
!  As soon as we've processed the first TAB, we're done.
!
  iput = lenc + 5 * ntab

  do iget = lenc, 1, - 1

    if ( s(iget:iget) /= TAB ) then

      if ( iput <= lens ) then
        s(iput:iput) = s(iget:iget)
      end if

      iput = iput - 1

    else

      do i = iput, iput - 5, -1
        if ( i <= lens ) then
          s(i:i) = ' '
        end if
      end do

      iput = iput - 6
      ntab = ntab - 1

      if ( ntab == 0 ) then
        return
      end if

    end if

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
  integer ( kind = 4 ) values(8)
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
