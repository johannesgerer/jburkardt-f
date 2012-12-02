program main

!*****************************************************************************80
!
!! MAIN is the main program for FILE_TRANSPOSE.
!
!  Discussion:
!
!    This program creates a "transposed" copy of a file.  That is, roughly
!    speaking, the I-th word of the J-th line becomes the J-th word of the 
!    I-th line.
!
!    Here a word is a sequence of characters delimited by line breaks,
!    spaces, or tabs.
!
!    This is only a clean operation if all the lines have the same number of
!    words.  
!    
!  Usage:
!
!    file_transpose input_filename output_filename
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2012
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
  logical :: verbose = .false.

  if ( verbose ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_TRANSPOSE'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) '  Make a "transposed" copy of a file.'
  end if
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
  call file_transpose ( input_filename, output_filename )
!
!  Terminate.
!
  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_TRANSPOSE:'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

  return
end
subroutine file_transpose ( input_filename, output_filename )

!*****************************************************************************80
!
!! FILE_TRANSPOSE makes a copy of a file with the words tranposed.
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
!      This of and
!      is three their
!      the little tails.
!      tale pigs
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the input filename.
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output filename.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) input_column
  character ( len = * ) input_filename
  character ( len = 255 ) input_line
  integer ( kind = 4 ) input_row
  integer ( kind = 4 ) input_unit
  character ( len = 255 ) input_word
  integer ( kind = 4 ) input_word_length
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) output_column
  character ( len = * ) output_filename
  character ( len = 255 ) output_line
  integer ( kind = 4 ) output_unit
!
!  Open the original file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_TRANSPOSE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file:'
    write ( *, '(4x,a)' ) '"' // trim ( input_filename ) // '".'
    return
  end if
!
!  Open the transpose unit.
!
  call get_unit ( output_unit )
  open ( unit = output_unit, file = output_filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )
!
!  Line by line, construct the transpose file.
!
  input_column = 0

  do

    input_column = input_column + 1
    input_row = 0
    output_column = 0
    output_line = ''

    do 

      read ( input_unit, '(a)', iostat = ios ) input_line

      if ( ios /= 0 ) then
        exit
      end if

      call s_word_find ( input_line, input_column, input_word, input_word_length )

      if ( 0 < input_word_length ) then
        output_column = output_column + 1
        if ( output_column == 1 ) then
          output_line = input_word
        else
          output_line = trim ( output_line ) // ' ' // trim ( input_word )
        end if
      end if

    end do

    if ( output_column == 0 ) then
      exit
    end if

    write ( output_unit, '(a)' ) trim ( output_line )

    rewind ( input_unit )

  end do

  close ( unit = input_unit )
  close ( unit = output_unit )
 
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
subroutine s_word_find ( s, iword, word, nchar )

!*****************************************************************************80
!
!! S_WORD_FIND finds the word of a given index in a string.
!
!  Discussion:
!
!    A "word" is any string of nonblank characters, separated from other
!    words by one or more blanks or TABS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, integer ( kind = 4 ) IWORD, the index of the word to be
!    searched for.  If IWORD is positive, then the IWORD-th
!    word is sought.  If IWORD is zero or negative, then
!    assuming that the string has N words in it, the
!    N+IWORD-th word will be sought.
!
!    Output, character ( len = * ) WORD, the IWORD-th word of the
!    string, or ' ' if the WORD could not be found.
!
!    Output, integer ( kind = 4 ) NCHAR, the number of characters in WORD,
!    or 0 if the word could not be found.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iblank
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) iword
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jword
  integer ( kind = 4 ) kword
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  integer ( kind = 4 ) s_len
  character, parameter :: TAB = achar ( 9 )
  character ( len = * ) word

  ilo = 0
  ihi = 0
  s_len = len_trim ( s )
  word = ''
  nchar = 0

  if ( s_len <= 0 ) then
    return
  end if

  if ( 0 < iword ) then

    if ( s(1:1) == ' ' .or. s(1:1) == TAB ) then
      iblank = 1
      jword = 0
      jlo = 0
      jhi = 0
    else
      iblank = 0
      jword = 1
      jlo = 1
      jhi = 1
    end if

    i = 1

    do

      i = i + 1

      if ( s_len < i ) then

        if ( jword == iword ) then
          ilo = jlo
          ihi = s_len
          nchar = s_len + 1 - jlo
          word = s(ilo:ihi)
        else
          ilo = 0
          ihi = 0
          nchar = 0
          word = ' '
        end if

        return

      end if

      if ( ( s(i:i) == ' ' .or. s(i:i) == TAB ) .and. iblank == 0 ) then

        jhi = i - 1
        iblank = 1
        if ( jword == iword ) then
          ilo = jlo
          ihi = jhi
          nchar = jhi + 1 - jlo
          word = s(ilo:ihi)
          return
        end if

      else if ( s(i:i) /= ' ' .and. s(i:i) /= TAB .and. iblank == 1 ) then

        jlo = i
        jword = jword + 1
        iblank = 0

      end if

    end do

  else

    iblank = 0
    kword = 1 - iword
    jword = 1
    jlo = s_len
    jhi = s_len
    i = s_len

    do

      i = i - 1

      if ( i <= 0 ) then

        if ( jword == kword ) then
          ilo = 1
          ihi = jhi
          nchar = jhi
          word = s(ilo:ihi)
        else
          ilo = 0
          ihi = 0
          nchar = 0
          word = ' '
        end if

        return

      end if

      if ( ( s(i:i) == ' ' .or. s == TAB ) .and. iblank == 0 ) then

        jlo = i + 1
        iblank = 1

        if ( jword == kword ) then
          ilo = jlo
          ihi = jhi
          nchar = jhi + 1 - jlo
          word = s(ilo:ihi)
          return
        end if

      else if ( s(i:i) /= ' ' .and. s(i:i) /= TAB .and. iblank == 1 ) then

        jhi = i
        jword = jword + 1
        iblank = 0

      end if

    end do

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
