program main

!*****************************************************************************80
!
!! MAIN is the main program for FILE_ROW_REVERSE.
!
!  Usage:
!
!    file_row_reverse input_filename output_filename
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  character ( len = 255 ) input_filename
  character ( len = 255 ) output_filename

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_ROW_REVERSE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Copy a file with the lines in reverse order.'
!
!  Count the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the input filename.
!
  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, input_filename )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the input filename:'
    read ( *, '(a)' ) input_filename
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input file is  "' // trim ( input_filename ) // '".'
!
!  Get the output filename.
!
  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, output_filename )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the output filename:'
    read ( *, '(a)' ) output_filename
  end if
!
!  Copy the input file to the output.
!
  call file_reverse_rows ( input_filename, output_filename )

  write ( *, '(a)' ) '  Output file is "' // trim ( output_filename ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_ROW_REVERSE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine file_reverse_rows ( input_file_name, output_file_name )

!*****************************************************************************80
!
!! FILE_REVERSE_ROWS makes a copy of a file with the lines in reverse order.
!
!  Example:
!
!    Input file:
!
!      This is the tale
!      of three little pigs
!      and their tails.
!
!    Output file:
!
!      and their tails.
!      of three little pigs
!      This is the tale
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the file to
!    be reversed.
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the name of the file to be
!    created, contained a reversed copy of the input file.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) input_count
  character ( len = * ) input_file_name
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  character ( len = * ) output_file_name
  integer ( kind = 4 ) output_unit
!
!  Open the input file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_REVERSE_ROWS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file:'
    write ( *, '(4x,a)' ) '"' // trim ( input_file_name ) // '".'
    return
  end if
!
!  Move to the end of the input file.
!
  input_count = 0

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    input_count = input_count + 1

  end do
!
!  Open the output file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file_name, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_REVERSE_ROWS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file:'
    write ( *, '(4x,a)' ) '"' // trim ( output_file_name ) // '".'
    return
  end if
!
!  Read backwards.
!
  backspace ( unit = input_unit, iostat = ios )

  do i = input_count, 1, -1

    backspace ( unit = input_unit, iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_REVERSE_ROWS - Fatal error!'
      write ( *, '(a)' ) '  IOS nonzero on backspace.'
      exit
    end if

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)') 'FILE_REVERSE_ROWS - Fatal error!'
      write ( *, '(a)' ) '  IOS nonzero in read'
      exit
    end if

    write ( output_unit, '(a)' ) trim ( line )

    backspace ( unit = input_unit )

  end do

  close ( unit = input_unit )

  endfile ( unit = output_unit )
  close ( unit = output_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  The input file had ', input_count, ' lines.'

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
!    26 October 2008
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
