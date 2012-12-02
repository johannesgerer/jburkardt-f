program main

!*****************************************************************************80
!
!! MAIN is the main program for READ_ALIGN.
!
!  Discussion:
!
!    READ_ALIGN reads a multiple alignment file and rewrites it in PIR format.
!
!    The program can optionally read in pairs of numbers in free format 
!    from yet another file to specify fragments.
!
!    The program is primarily intended to create a PIR format file
!    that can be "fed" to the DISTANCES program.
!
!  Modified:
!
!    29 June 2000
!
!  Author:
!
!    Des Higgins, 
!    EMBL Data Library, 
!    April 1992
!
  implicit none

  integer ( kind = 4 ), parameter :: line_len = 200
  integer ( kind = 4 ), parameter :: maxl = 5000
  integer ( kind = 4 ), parameter :: sequence_max = 500
  integer ( kind = 4 ), parameter :: name_len = 16

  logical blank_line
  integer ( kind = 4 ) block_end
  character ch
  integer ( kind = 4 ) found
  integer ( kind = 4 ) i
  character ( len = line_len ) input_line
  character ( len = 100 ) input_name
  integer ( kind = 4 ), parameter :: input_unit = 44
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  character ( len = line_len ) linbuff
  integer ( kind = 4 ) name_end
  integer ( kind = 4 ) name_start
  character ( len = 100 ) output_name
  integer ( kind = 4 ), parameter :: output_unit = 45
  integer ( kind = 4 ) segment_from
  character ( len = 100 ) segment_name
  integer ( kind = 4 ) segment_to
  integer ( kind = 4 ), parameter :: segment_unit = 46
  logical segment_use
  integer ( kind = 4 ) sequence_length
  character ( len = name_len ) sequence_name(sequence_max)
  integer ( kind = 4 ) sequence_number
  character sequence(sequence_max,maxl)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'READ_ALIGN'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read in a multiple alignment file,'
  write ( *, '(a)' ) '  Write out a version using the PIR format.'

  found = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the input file name:'
  read ( *, '(a)' ) input_name

  open ( unit = input_unit, file = input_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_ALIGN - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the segment file name,'
  write ( *, '(a)' ) '  or <RETURN> if none:'
  read ( *, '(a)' ) segment_name

  if ( segment_name /= ' ' ) then
    segment_use = .true.
  else
    segment_use = .false.
  end if

  if ( segment_use ) then

    open ( unit = segment_unit, file = segment_name, status = 'old', &
      iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'READ_ALIGN - Fatal error!'
      write ( *, '(a)' ) '  Could not open the segment file.'
      stop
    end if

  end if
!
!  Skip the first few lines of text.
!  Read until completely blank line.
!
  do

    read ( input_unit, '(a)', iostat = ios ) input_line

    if ( ios /= 0 ) then
      go to 100
    end if

    if ( blank_line ( input_line ) ) then
      exit
    end if

  end do

  sequence_length = 0
!
!  Skip any more blank lines.
!
50  continue

  do

    read ( input_unit, '(a)', iostat = ios ) input_line

    if ( ios /= 0 ) then
      go to 100
    end if

    if ( .not. blank_line ( input_line ) ) then
      exit
    end if

  end do

  backspace ( input_unit )
  sequence_number  = 0
!
!  Now read a block of alignment.
!
  do

    read ( input_unit, '(a)', iostat = ios ) input_line

    if ( ios /= 0 ) then
      exit
    end if

    if ( blank_line ( input_line ) ) then
      sequence_length = sequence_length + found
      go to 50
    end if

    sequence_number = sequence_number + 1

    if ( sequence_number > sequence_max ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'READ_ALIGN - Fatal error!'
      write ( *, '(a)' ) '  Maximum number of sequences exceeded.'
      write ( *, '(a,i6)' ) '  Internal limit is ', sequence_max
      stop
    end if

    do i = line_len, 1, -1
      if ( input_line(i:i) /= ' ' ) then
        block_end = i
        exit
      end if
    end do
!
!  Find first non blank character (start of name).
!
    do i = 1, line_len
      if ( input_line(i:i) /= ' ' ) then
        name_start = i
        exit
      end if
    end do
!
!  Find first blank character after the name.
!
    do i = name_start, line_len
      if ( input_line(i:i) == ' ' ) then
        name_end = i
        exit
      end if
    end do

    found = 0
    do i = block_end, name_end, -1
      ch = input_line(i:i)
      if ( ch /= ' ' ) then
        found = found + 1
        linbuff(found:found) = ch
      end if
    end do

    do i = 1, found
      sequence(sequence_number,sequence_length+i) = linbuff(found-i+1:found-i+1)
    end do

    if ( sequence_length == 0 ) then
      sequence_name(sequence_number) = input_line(1:min(name_len,name_end))
    end if

  end do
!
!  Finished reading input file.
!
100   continue

  close ( unit = input_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of sequences read = ', sequence_number
  write ( *, '(a,i6)' ) '  Alignment length  = ', sequence_length
!
!  Write the output file.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the output file name:'
  read ( *, '(a)' ) output_name

  open ( unit = output_unit, file = output_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_ALIGN - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  do i = 1, sequence_number

    write ( output_unit, '(a,a)' ) '>P1;', sequence_name(i)
    write ( output_unit, '(a)' ) ' '

    if ( .not. segment_use ) then

      write ( output_unit, '(1x,60a1)' ) sequence(i,1:sequence_length)

    else

      rewind ( segment_unit )

      do

        read ( segment_unit, *, iostat = ios ) segment_from, segment_to

        if ( ios /= 0 ) then
          exit
        end if

        write ( output_unit, '(1x,60a1)' ) sequence(i,segment_from:segment_to)

      end do

    end if

    write ( output_unit, '(a)' ) '*'

  end do
!
!  Normal end of execution.
!
  close ( unit = output_unit )

  if ( segment_use ) then
    close ( unit = segment_unit )
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'READ_ALIGN:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function blank_line ( line )

!*****************************************************************************80
!
!! BLANK_LINE looks for a "blank" line.
!
!  Discussion:
!
!    A "blank" line is one that contains only spaces, dots,
!    asterisks or digits.  This only works for ASCII files.
!
!  Modified:
!
!    06 April 2000
!
!  Author:
!
!    Des Higgins, 
!    EMBL Data Library, 
!    April 1992
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, the line of text to be examined.
!
!    Output, logical BLANK_LINE, is TRUE if the line is "blank", and
!    FALSE otherwise.
!
  implicit none

  logical blank_line
  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ich
  character ( len = * ) line

  do i = 1, len_trim ( line )

    c = line(i:i)

    if ( c /= ' '  .and. c /= '.' .and. c /= '*' ) then

      ich = ichar ( c )
      if ( ich < 48 .or. 57 < ich ) then
        blank_line = .false.
        return
      end if

    end if

  end do

  blank_line = .true.

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
