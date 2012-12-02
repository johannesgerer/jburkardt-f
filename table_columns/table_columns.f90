program main

!*****************************************************************************80
!
!! MAIN is the main routine for TABLE_COLUMNS.
!
!  Discuaaion:
!
!    TABLE_COLUMNS extracts specific columns from a file.
!
!    The file is presumed to consist of lines of data.
!    Each line of data contains several words, separated by one or more
!    spaces or TAB characters.
!    The first "word" in each line is the first "column", and so on.
!
!  Example:
!
!    The input file contains:
!
!      #  myfile.txt
!      #
!         1   2    3    4    5
!        18  19   20   21   22
!      1440 778 9909 2333 0999
!
!    The command:
!
!      table_columns 3 input_file output_file
!
!    will create an output file containing column 3:
!
!      #  myfile.txt
!      #
!         3
!        20
!      9909
!
!    More than one column may be requested.  Simply separate them
!    by spaces.  To specify a range of columns, use a colon.
!
!    The command:
!
!      table_columns 2:4 1 3 input_file output_file
!
!    will create an output file containing columns 2, 3, 4, 1 3:
!
!      #  myfile.txt
!      #
!         2    3    4    1    3
!        19   20   21   18   20
!       778 9909 2333 1440 9909
!
!    If a particular line is missing a column of data, then the
!    program will behave as though that value was "***".
!
!  Usage:
!
!    table_columns column_spec input_file output_file
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Local parameters:
!
!    Local, integer LIST_MAX, limits the total number of columns that
!    can be copied.
!
!    Local, integer VAL_MAX, limits the total number of columns that can
!    be specified in a single "colon" command.  Thus 2:20 will be legal,
!    but 2:50 is asking for 49 columns, which is too many.
!
  implicit none

  integer ( kind = 4 ), parameter :: list_max = 25
  integer ( kind = 4 ), parameter :: val_max = 25

  logical, parameter :: debug = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  character ( len = 255 ) input_filename
  integer ( kind = 4 ) list(list_max)
  integer ( kind = 4 ) list_num
  integer ( kind = 4 ) numarg
  character ( len = 255 ) output_filename
  integer ( kind = 4 ) val(val_max)
  integer ( kind = 4 ) val_num
  character ( len = 255 ) word

  if ( debug ) then

    call timestamp ( )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TABLE_COLUMNS'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) '  Extract a column of data from a file.'

  end if

  ierror = 0
  numarg = iargc ( )
!
!  Get the column numbers.
!
  if ( numarg <= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Usage:'
    write ( *, '(a)' ) '  table_columns column_list input_file output_file'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  column_list is a list of columns, '
    write ( *, '(a)' ) '  separated by spaces.  You may include a column range,'
    write ( *, '(a)' ) '  listing an initial and final column with a colon.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Example:'
    write ( *, '(a)' ) '  table_columns 2:4 1 3 input output'
    stop
  end if
!
!  Read list of columns.
!
  list_num = 0

  do iarg = 1, numarg - 2

    call getarg ( iarg, word )

    call ranger ( word, val_max, val_num, val )

    do i = 1, val_num
      list_num = list_num + 1
      list(list_num) = val(i)
    end do

  end do
!
!  Argument N-1 is presumed to be input filename.
!
  iarg = numarg-1
  call getarg ( iarg, input_filename )
!
!  Argument N is presumed to be output filename.
!
  iarg = numarg
  call getarg ( iarg, output_filename )

  call copy_columns ( input_filename, output_filename, list_num, list )
!
!  Terminate.
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TABLE_COLUMNS'
    write ( *, '(a)' ) '  Normal end of execution.'
    call timestamp ( )
  end if

  stop
end
subroutine copy_columns ( input_filename, output_filename, list_num, list )

!*****************************************************************************80
!
!! COPY_COLUMNS copies specified columns from one file to another.
!
!  Discussion:
!
!    Each line of the input file is either
!    a) BLANK
!    b) a COMMENT (beginning with the '#' character in column 1;
!    c) a STRING of symbols separated by spaces.
!
!    We think of each line that is a "string" as being composed of
!    words or columns of data.  Our intent is to make a new file,
!    in which certain words of the input file are copied to the
!    output file, in any order, and with any amount of repetition.
!
!    Thus, if one line of the input file reads:
!
!      ant  botfly  coati  duykerbok
!
!    and the list of columns is given as
!
!      3 3 2 4 3
!
!    then the corresponding line of the output file will read
!
!      coati coati botfly duykerbok coati
!
!    Meanwhile, blanks and comments are copied directly from the
!    input to the output file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the
!    file which is to be read in.
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the
!    file to be created.
!
!    Input, integer ( kind = 4 ) LIST_NUM, the number of column indices.
!
!    Input, integer ( kind = 4 ) LIST(LIST_NUM), the column indices.
!
  implicit none

  integer ( kind = 4 ) list_num
  integer ( kind = 4 ), parameter :: word_max = 100

  integer ( kind = 4 ) i
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 1000 ) line
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) list(list_num)
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
  character ( len = 255 ) string
  integer ( kind = 4 ) w
  integer ( kind = 4 ) word_end(word_max)
  integer ( kind = 4 ) word_num
  integer ( kind = 4 ) word_start(word_max)
!
!  Open the input file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COPY_COLUMNS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '".'
    return
  end if
!
!  Open the output file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COPY_COLUMNS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '".'
    return
  end if
!
!  Read a line of input.
!  If it's blank, or a comment, skip it.
!  Otherwise, extract the columns, and write them out.
!
  line_num = 0

  do

    line_num = line_num + 1

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      write ( output_unit, '(a)' ) ' '
      cycle
    end if

    if ( line(1:1) == '#' ) then
      write ( output_unit, '(a)' ) trim ( line )
      cycle
    end if
!
!  Now we have the line.
!
    call s_word_count ( line, word_num )
!
!  Make a word index for this line.
!
    call word_bounds ( line, word_num, word_start, word_end )
!
!  Print the selected words into a string.
!
    string = ' '

    do i = 1, list_num
      w = list(i)
      if ( 0 < w .and. w <= word_num ) then
        j = word_start(w)
        k = word_end(w)
        call s_cat2 ( string, line(j:k), string )
      else
        call s_cat2 ( string, '***', string )
      end if
    end do
!
!  Print the string.
!
    write ( output_unit, '(a)' ) trim ( string )

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
subroutine ranger ( s, maxval, nval, ival )

!*****************************************************************************80
!
!! RANGER "understands" a range defined by a string like '4:8'.
!
!  Discussion:
!
!    The range can be much more complicated, as in
!
!      '4:8 10 2 14:20'
!
!    or (commas are optional)
!
!      '4:8,10, 2 , 14:20'
!
!    RANGER will return the values
!
!      4, 5, 6, 7, 8, 10, 2, 14, 15, 16, 17, 18, 19, 20
!
!    0 and negative integers are acceptable.  So are pairs
!    of values that are equal, as in '4:4', which just represents
!    4, and pairs that represent descending sequences, as
!    in '4:-2' which represents 4, 3, 2, 1, 0, -1, -2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, contains a string of integers,
!    representing themselves, and pairs of integers representing
!    themselves and all integers between them.
!
!    Input, integer ( kind = 4 ) MAXVAL, the dimension of the IVAL vector,
!    which represents the maximum number of integers that may
!    be read from the string.
!
!    Output, integer ( kind = 4 ) NVAL, the number of integers read
!    from the string.
!
!    Output, integer ( kind = 4 ) IVAL(MAXVAL).  The first NVAL entries of
!    IVAL contain the integers read from the string.
!
  implicit none

  integer ( kind = 4 ) maxval

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ival(maxval)
  integer ( kind = 4 ) length
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nval
  character ( len = * ) s

  nval = 0
!
!  Replace all commas by blanks.
!
  call s_ch_blank ( s, ',' )
!
!  Replace multiple consecutive blanks by one blank.
!
  call s_blanks_delete ( s )
!
!  Get the length of the string to the last nonblank.
!
  lens = len_trim ( s )
!
!  Set a pointer to the next location to be examined.
!
  next = 1

  do while ( next <= lens )
!
!  Find the next integer in the string.
!
    call s_to_i4 ( s(next:), intval, ierror, length )

    if ( ierror /= 0 ) then
      return
    end if
!
!  Move the pointer.
!
    next = next + length
!
!  If there's room, add the value to the list.
!
    if ( maxval <= nval ) then
      return
    end if

    nval = nval + 1
    ival(nval) = intval
!
!  Have we reached the end of the string?
!
    if ( lens < next ) then
      return
    end if
!
!  Skip past the next character if it is a space.
!
    if ( s(next:next) == ' ' ) then
      next = next + 1
      if ( lens < next ) then
        return
      end if
    end if
!
!  Is the next character a colon?
!
    if ( s(next:next) /= ':' ) then
      cycle
    end if
!
!  Increase the pointer past the colon.
!
    next = next + 1

    if ( lens < next ) then
      return
    end if
!
!  Find the next integer in the string.
!
    call s_to_i4 ( s(next:), intval, ierror, length )

    if ( ierror /= 0 ) then
      return
    end if
!
!  Move the pointer.
!
    next = next + length
!
!  Generate integers between the two values.
!
    ilo = ival(nval)

    if ( ilo <= intval ) then
      inc = + 1
    else
      inc = - 1
    end if

    do i = ilo+inc, intval, inc

      if ( maxval <= nval ) then
        return
      end if

      nval = nval + 1
      ival(nval) = i

    end do

  end do

  return
end
subroutine s_blanks_delete ( s )

!*****************************************************************************80
!
!! S_BLANKS_DELETE replaces consecutive blanks by one blank.
!
!  Discussion:
!
!    The remaining characters are left justified and right padded with blanks.
!    TAB characters are converted to spaces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character newchr
  character oldchr
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  j = 0
  newchr = ' '

  do i = 1, len ( s )

    oldchr = newchr
    newchr = s(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    s(i:i) = ' '

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      s(j:j) = newchr
    end if

  end do

  return
end
subroutine s_cat2 ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT2 concatenates two strings to make a third string.
!
!  Discussion:
!
!    This variation of S_CAT removes any initial blanks from S2,
!    and then inserts two blanks between S1 and S2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  if ( s1 == ' ' .and. s2 == ' ' ) then
    s3 = ' '
  else if ( s1 == ' ' ) then
    s3 = s2
  else if ( s2 == ' ' ) then
    s3 = s1
  else
    s3 = trim ( s1 ) // '  ' // trim ( adjustl ( s2 ) )
  end if

  return
end
subroutine s_ch_blank ( s, c )

!*****************************************************************************80
!
!! S_CH_BLANK replaces each occurrence of a particular character by a blank.
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
!    Input/output, character ( len = * ) S, the string to be transformed.
!
!    Input, character C, the character to be removed.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) s

  nchar = len_trim ( s )

  do i = 1, nchar

    if ( s(i:i) == c ) then
      s(i:i) = ' '
    end if

  end do

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
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
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If STRING is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character of S used.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine s_word_count ( s, word_num )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) WORD_NUM, the number of "words"
!    in the string.  Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  character ( len = * ) s
  integer ( kind = 4 ) s_len
  integer ( kind = 4 ) word_num

  word_num = 0
  s_len = len ( s )

  if ( s_len <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, s_len

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      word_num = word_num + 1
      blank = .false.
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
subroutine word_bounds ( line, word_num, word_start, word_end )

!*****************************************************************************80
!
!! WORD_BOUNDS returns the start and end of each word in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string containing words
!    separated by spaces.
!
!    Input, integer ( kind = 4 ) WORD_NUM, the number of words in the line.
!
!    Output, integer ( kind = 4 ) WORD_START(WORD_NUM), WORD_END(WORD_NUM),
!    the locations in LINE of the beginning and end of each word.
!
  implicit none

  integer ( kind = 4 ) word_num

  logical blank
  character c
  logical, parameter :: debug = .true.
  integer ( kind = 4 ) i
  character ( len = * ) line
  integer ( kind = 4 ) line_len
  integer ( kind = 4 ) w
  integer ( kind = 4 ) word_end(word_num)
  integer ( kind = 4 ) word_start(word_num)

  i = 0
  w = 0
  blank = .true.

  line_len = len_trim ( line )

  do i = 1, line_len + 1

    if ( i <= line_len ) then
      c = line(i:i)
    else
      c = ' '
    end if

    if ( c == ' ' ) then

      if ( .not. blank ) then
        word_end(w) = i-1
        if ( w == word_num ) then
          exit
        end if
      end if

      blank = .true.

    else

      if ( blank ) then
        w = w + 1
        word_start(w) = i
      end if

      blank = .false.

    end if

  end do

  if ( w /= word_num ) then
    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WORD_BOUNDS - Warning:'
      write ( *, '(a)' ) '  Found fewer words than requested.'
    end if
  end if

  return
end
subroutine word_find ( s, iword, word, nchar )

!*****************************************************************************80
!
!! WORD_FIND finds the word of a given index in a string.
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
!    01 April 2001
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
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
  character ( len = * ) word

  ilo = 0
  ihi = 0
  length = len_trim ( s )

  if ( length <= 0 ) then
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

      if ( length < i ) then

        if ( jword == iword ) then
          ilo = jlo
          ihi = length
          nchar = length + 1 - jlo
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
    jlo = length
    jhi = length
    i = length

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
