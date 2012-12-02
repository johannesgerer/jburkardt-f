subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  character ch
  integer ( kind = 4 ) itemp

  itemp = ichar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( ch1, ch2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
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
!    Input, character CH1, CH2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character ch1
  character ch1_cap
  character ch2
  character ch2_cap

  ch1_cap = ch1
  ch2_cap = ch2

  call ch_cap ( ch1_cap )
  call ch_cap ( ch2_cap )

  if ( ch1_cap == ch2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
function ch_index_last ( s, ch )

!*****************************************************************************80
!
!! CH_INDEX_LAST is the last occurrence of a character in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character CH, the character to be searched for.
!
!    Output, integer ( kind = 4 ) CH_INDEX_LAST, the location of the last 
!    occurrence of the character in the string, or -1 if it does not occur.
!
  implicit none

  character ch
  integer ( kind = 4 ) ch_index_last
  integer ( kind = 4 ) i
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  ch_index_last = -1
  s_length = len_trim ( s )

  do i = s_length, 1, -1

    if ( s(i:i) == ch ) then
      ch_index_last = i
      return
    end if

  end do
 
  return
end
function ch_is_digit ( ch )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character CH, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if the character is a digit, 
!    FALSE otherwise.
!
  implicit none

  character ch
  logical ch_is_digit

  if ( lge ( ch, '0' ) .and. lle ( ch, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
  end if

  return
end
subroutine ch_low ( ch )

!*****************************************************************************80
!
!! CH_LOW lowercases a single character.
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
!    Input/output, character CH, the character to be lowercased.
!
  implicit none

  character ch
  integer ( kind = 4 ) itemp

  itemp = ichar ( ch )

  if ( 65 <= itemp .and. itemp <= 90 ) then
    ch = char ( itemp + 32 )
  end if

  return
end
subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     CH  DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character CH, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  
!    If the character was 'illegal', then DIGIT is -1.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( lge ( ch, '0' ) .and. lle ( ch, '9' ) ) then

    digit = ichar ( ch ) - 48

  else if ( ch == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
function ch_to_rot13 ( ch )

!*****************************************************************************80
!
!! CH_TO_ROT13 converts a character to its ROT13 equivalent.
!
!  Discussion:
!
!    Two applications of CH_TO_ROT13 to a character will return the original.
!
!    As a further scrambling, digits are similarly rotated using
!    a "ROT5" scheme.
!
!  Example:
!
!    Input:  Output:
!
!    a       n
!    C       P
!    J       W
!    1       6
!    5       0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character CH, the character to be converted.
!
!    Output, character CH_TO_ROT13, the ROT13 equivalent of the character.
!
  implicit none

  character ch
  character ch_to_rot13
  integer ( kind = 4 ) itemp

  itemp = ichar ( ch )
!
!  [0:4] -> [5:9]
!
  if ( 48 <= itemp .and. itemp <= 52 ) then
    itemp = itemp + 5
!
!  [5:9] -> [0:4]
!
  else if ( 53 <= itemp .and. itemp <= 57 ) then
    itemp = itemp - 5
!
!  [A:M] -> [N:Z]
!
  else if ( 65 <= itemp .and. itemp <= 77 ) then
    itemp = itemp + 13
!
!  [N:Z] -> [A:M]
!
  else if ( 78 <= itemp .and. itemp <= 90 ) then
    itemp = itemp - 13
!
!  [a:m] -> [n:z]
!
  else if ( 97 <= itemp .and. itemp <= 109 ) then
    itemp = itemp + 13
!
!  [n:z] -> [a:m]
!
  else if ( 110 <= itemp .and. itemp <= 122 ) then
    itemp = itemp - 13
  end if

  ch_to_rot13 = char ( itemp )

  return
end
subroutine digit_inc ( c )

!*****************************************************************************80
!
!! DIGIT_INC increments a decimal digit.
!
!  Example:
!
!    Input  Output
!    -----  ------
!    '0'    '1'
!    '1'    '2'
!    ...
!    '8'    '9'
!    '9'    '0'
!    'A'    'A'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, a digit to be incremented.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  call ch_to_digit ( c, digit )

  if ( digit == -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, c )

  return
end
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine file_advance_to_string ( iunit, s, line, ierror )

!*****************************************************************************80
!
!! FILE_ADVANCE_TO_STRING searches ahead in a text file for a string.
!
!  Discussion:
!
!    The file should already have been opened.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the unit number associated 
!    with the open file.
!
!    Input, character ( len = * ) S, a string to search for.
!
!    Output, character ( len = * ) LINE:
!    If IERROR = 0, the line of the file that was just read,
!      and which contains the string.
!    If IERROR = 1, then LINE is blank.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error, the string was found.
!    1, error, the end of the file was reached.
!
  implicit none

  logical, parameter :: DEBUG = .true.
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = * ) line
  integer ( kind = 4 ) line_num
  character ( len = * ) s

  ierror = 0
  line_num = 0

  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

    if ( index ( line, trim ( s ) ) /= 0 ) then
      return
    end if

  end do

  line = ' '
  ierror = 1

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ADVANCE_TO_STRING - Warning!'
    write ( *, '(a)' ) '  Did not find the string:'
    write ( *, '(a)' ) '    ' // trim ( s )
    write ( *, '(a,i8)' ) '  Number of lines read was ', line_num
  end if

  return
end
subroutine file_append ( filename, ierror, iunit, rec_num )

!*****************************************************************************80
!
!! FILE_APPEND allows a user to append new information to an old file.
!
!  Discussion:
!
!    This routine was created to address the fact that ANSI FORTRAN
!    does not let one easily append information to a sequential
!    access file once it has been closed.  In order to allow a user
!    to append new information, we create a new, writeable copy
!    of the file by means of a temporary copy.
!
!    On input, the file should not be open.  On output, the file is
!    open, the file is writeable, and the file I/O pointer is
!    ready to write a new line following the last line of the
!    original contents of the file.
!
!    It is assumed that each line of the file is no longer than
!    256 characters.
!
!    The copied lines will not have any trailing blanks.
!
!    A temporary file will be created with the name "append.tmp".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file to
!    which more information is to be appended.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred while trying to open the new file.
!
!    Output, integer ( kind = 4 ) IUNIT, the FORTRAN unit number associated 
!    with the file.
!
!    Output, integer ( kind = 4 ) REC_NUM, the number of records (that is, 
!    lines) that are already in the file.
!
  implicit none

  character ( len = * ) filename
  character ( len = 255 ) file_temp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) iunit2
  character ( len = 255 ) line
  integer ( kind = 4 ) rec_num

  ierror = 0

  file_temp = 'append.tmp'
!
!  Open old file as readable.  If it doesn't exist, we can
!  skip ahead.  Otherwise, also open new file as writeable.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_APPEND - Note:'
    write ( *, '(a)' ) '  This is a new file.'

    open ( unit = iunit, file = filename, status = 'new', iostat = ios )

    if ( ios /= 0 ) then
      ierror = 4
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_APPEND - Fatal error!'
      write ( *, '(a)' ) '  Unexpected error while opening the file:'
      write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
      return
    end if

    rec_num = 0
 
    return
  end if

  rewind iunit

  call get_unit ( iunit2 )

  open ( unit = iunit2, file = file_temp, status = 'new', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_APPEND - Fatal error!'
    write ( *, '(a)' ) '  Unexpected error while opening the file:'
    write ( *, '(4x,a)' ) '"' // trim ( file_temp ) // '".'
    return
  end if
!
!  Copy data from old file into temporary file.
!
  rec_num = 0
 
  do
 
    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    rec_num = rec_num + 1
    write ( iunit2, '(a)' ) trim ( line )
 
  end do
 
  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_APPEND - Note:'
  write ( *, '(a,i8)' ) '  The number of records in the file is ', rec_num
!
!  Delete the old copy of the file.
!
  call file_delete ( filename )
!
!  Mark the end of the temporary file, close it, 
!  then reopen it with "OLD" status.
!
  endfile ( unit = iunit2 )

  close ( unit = iunit2 )

  open ( unit = iunit2, file = file_temp, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_APPEND - Fatal error!'
    write ( *, '(a)' ) '  Unexpected error while opening the file:'
    write ( *, '(4x,a)' ) '"' // trim ( file_temp ) // '".'
    return
  end if

  rewind iunit2
!
!  Create a new version of the old file, opening it with
!  "STATUS = 'NEW'" so that it is writeable.
!
  open ( unit = iunit, file = filename, status = 'new', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_APPEND - Fatal error!'
    write ( *, '(a)' ) '  Unexpected error while opening the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    return
  end if
!
!  Read each line from the temporary file, and copy it
!  back into the old file.
!
  do
 
    read ( iunit2, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    write ( iunit, '(a)' ) trim ( line )
 
  end do

  close ( unit = iunit2 )
!
!  Delete the temporary file, and return.
!
  call file_delete ( file_temp )
 
  return
end
subroutine file_char_count ( filename, char_num )

!*****************************************************************************80
!
!! FILE_CHAR_COUNT counts the number of characters in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file.
!
!    Output, integer ( kind = 4 ) CHAR_NUM, the number of characters in 
!    the file.
!
  implicit none

  integer ( kind = 4 ) char_num
  character ( len = * ) filename
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 255 ) line

  char_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    char_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_CHAR_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    return
  end if
!
!  Read the lines, count the characters.
!
!  Right now, I don't know how to determine which characters in LINE
!  were read from the file and which represent padding.
!
!  I tried going to the last nonnull, but it looks like, at least on
!  my system, LINE gets filled with blanks.  So I'll just count 
!  til the last non-blank, but that doesn't distinguish between a
!  blank that really was in the original file and should be counted,
!  and a blank that is just padding.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    char_num = char_num + len_trim ( line )

  end do

  close ( unit = iunit )

  return
end
subroutine file_column_count ( filename, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file are presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some 
!    comment lines, which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer ( kind = 4 ) column_num
  character ( len = * ) filename
  logical got_one
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 255 ) line
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    got_one = .true.
    exit

  end do

  if ( .not. got_one ) then

    rewind ( iunit )

    do 

      read ( iunit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = iunit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

  return
end
subroutine file_column_range ( filename, column_num, col_min, col_max )

!*****************************************************************************80
!
!! FILE_COLUMN_RANGE determines the minimum and maximum ranges of each column.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Each line of the file is presumed to consist of COLUMN_NUM real numbers,
!    separated by spaces.
!
!    The routine computes the range of each column.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
!    Input, integer ( kind = 4 ) COLUMN_NUM, the number of columns assumed 
!    to be in the file.
!
!    Output, real ( kind = 8 ) COL_MIN(COLUMN_NUM), COL_MAX(COLUMN_NUM), 
!    the minimum and maximum for each column.
!
  implicit none

  integer ( kind = 4 ) column_num

  real ( kind = 8 ) col_max(column_num)
  real ( kind = 8 ) col_min(column_num)
  character ( len = * ) filename
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  integer ( kind = 4 ) nrow
  real ( kind = 8 ) x(column_num)
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_RANGE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    return
  end if

  nrow = 0
  col_min(1:column_num) = 0.0D+00
  col_max(1:column_num) = 0.0D+00

  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    call s_to_r8vec ( line, column_num, x, ierror )

    if ( ierror /= 0 ) then
      exit
    end if

    nrow = nrow + 1

    if ( nrow == 1 ) then
      col_min(1:column_num) = x(1:column_num)
      col_max(1:column_num) = x(1:column_num)
    else
      do j = 1, column_num
        col_min(j) = min ( col_min(j), x(j) )
        col_max(j) = max ( col_max(j), x(j) )
      end do
    end if

  end do

  close ( unit = iunit )

  return
end
subroutine file_copy ( old_filename, new_filename, ierror )

!*****************************************************************************80
!
!! FILE_COPY makes a copy of a file.
!
!  Discussion:
!
!    The file is assumed to be sequential access, with variable 
!    length records.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OLD_filename, the name of the file 
!    to be copied.
!
!    Input, character ( len = * ) NEW_filename, the name of the copy of 
!    the file.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    1, the file names are the same.
!    2, a free unit number could not be found for the old file.
!    3, the routine could not open the old file.
!    4, a free unit number could not be found for the new file.
!    5, the routine could not open the new file.
!
  implicit none

  logical file_exist
  logical file_is_open
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  character ( len = * ) new_filename
  integer ( kind = 4 ) new_unit
  character ( len = * ) old_filename
  integer ( kind = 4 ) old_unit

  ierror = 0
!
!  Does the original file exist?
!
  if ( .not. file_exist ( old_filename ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  The old file does not exist.'
    return
  end if
!
!  Is the original file open?
!
  if ( file_is_open ( old_filename ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  The old file is open.'
    write ( *, '(a)' ) '  It must be closed before it can be copied.'
    return
  end if
!
!  Make sure the file names aren't the same.
!
  if ( new_filename == old_filename ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  The old and new file names are identical.'
    return
  end if
!
!  Does the new file exist?
!
  if ( file_exist ( new_filename ) ) then

    if ( file_is_open ( new_filename ) ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
      write ( *, '(a)' ) '  A file is already open with the new name.'
      write ( *, '(a)' ) '  It must be closed before it can be overwritten.'
      return
    end if
  
    call file_delete ( new_filename )

  end if
!
!  At this point:
!
!    The old file exists, and is not open.
!    The new file does not exist, and has a different name.
!
!  Open the old file.
!
  call get_unit ( old_unit )

  if ( old_unit == 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  Could not get a unit number for the old file.'
    return
  end if

  open ( unit = old_unit, file = old_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  Could not open the old file:'
    write ( *, '(4x,a)' ) '"' // trim ( old_filename ) // '".'
    return
  end if
!
!  Open the new file.
!
  call get_unit ( new_unit )

  if ( new_unit == 0 ) then
    ierror = 4
    close ( unit = old_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free unit for the copy file.'
    return
  end if

  open ( unit = new_unit, file = new_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 5
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COPY - Fatal error!'
    write ( *, '(a)' ) '  Could not open the new file:'
    write ( *, '(4x,a)' ) '"' // trim ( new_filename ) // '".'
    close ( unit = old_unit )
    return
  end if
!
!  Read an old line, write a new line.
!
  do
 
    read ( old_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    write ( new_unit, '(a)' ) trim ( line )

  end do
 
  close ( unit = old_unit )

  endfile ( unit = new_unit )
  close ( unit = new_unit )

  return
end
subroutine file_delete ( filename )

!*****************************************************************************80
!
!! FILE_DELETE deletes a named file if it exists.
!
!  Discussion:
!
!    You might want to call this routine to get rid of any old copy
!    of a file, before trying to open a new copy with the OPEN argument:
!      status = 'new'.
!
!    It's not always safe to open a file with " STATUS = 'UNKNOWN' ".
!    For instance, on the SGI, the most recent version of the FORTRAN
!    compiler seems to go crazy when I open an unformatted direct
!    access file this way.  It creates an enormous file (of somewhat
!    random size).  The problem goes away if I delete any old copy
!    using this routine, and then open a fresh copy with
!    " STATUS = 'NEW' ".  It's a scary world.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
  implicit none

  logical file_exist
  logical file_is_open
  character ( len = * ) filename
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical, parameter :: verbose = .false.
!
!  Does the file exist?
!
  if ( .not. file_exist ( filename ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE - Warning!'
    write ( *, '(a)' ) '  There is no file of the given name.'
    return
  end if
!
!  Is the file open?
!
  if ( file_is_open ( filename ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE - Warning!'
    write ( *, '(a)' ) '  The file is currently open.'
    write ( *, '(a)' ) '  It must be closed before it can be deleted.'
    return
  end if
!
!  Get a free unit number.
!
  call get_unit ( iunit )

  if ( iunit == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE: Warning!'
    write ( *, '(a)' ) '  A free FORTRAN unit could not be found.'
    return
  end if

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE:'
    write ( *, '(a)' ) '  Deleting "' // trim ( filename ) // '".'
  end if

  open ( unit = iunit, file = filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE: Warning!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    return
  end if

  close ( unit = iunit, status = 'delete' )

  return
end
function file_exist ( filename )

!*****************************************************************************80
!
!! FILE_EXIST reports whether a file exists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
!    Output, logical FILE_EXIST, is TRUE if the file exists.
!
  implicit none

  logical file_exist
  character ( len = * ) filename

  inquire ( file = filename, exist = file_exist )

  return
end
subroutine file_get_next_integer ( file_unit, more, value )

!*****************************************************************************80
!
!! FILE_GET_NEXT_INTEGER returns the next integer from a file.
!
!  Discussion:
!
!    The file should have been opened before calling this routine.
!
!    The routine will read ANY string of characters that looks like an integer,
!    and does not require a specific separate.  Garbage characters are ignored.
!
!    If a very long string of digits is encountered, the routine will trigger
!    integer overflow as it tries to pack them all into one integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) FILE_UNIT, the unit number associated with 
!    the file.
!
!    Input/output, logical MORE.  On first call, the user should set MORE to 
!    be FALSE, which signals the routine to initialize.  Each time the 
!    routine reads another integer from the file, it returns it in VALUE, 
!    and returns MORE as TRUE.  When no more integers can be read, MORE 
!    is returned as FALSE.
!  
!    Output, integer ( kind = 4 ) VALUE, the next integer in the file.  However,
!    if MORE is FALSE on output, then VALUE is set to 0.
!
  implicit none

  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) file_status
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  character ( len = 255 ) line
  logical more
  integer ( kind = 4 ) value

  save line

  if ( .not. more ) then
    line = ' '
    more = .true.
  end if
!
!  Try to extract an integer from the current line.
!
  do
!
!  If LINE is blank, read another line from the file.
!
    do while ( len_trim ( line ) == 0 )

      read ( file_unit, '(a)', iostat = file_status ) line

      if ( file_status /= 0 ) then
        more = .false.
        value = 0
        line = ' '
        return
      end if

    end do
!
!  Try to extract the next integer from LINE.
!
    call s_to_i4 ( line, value, ierror, length )
!
!  If we got a value, then chop out the used bits, and return.
!
    if ( ierror == 0 ) then
      line(1:length) = ' '
      line = adjustl ( line )
      exit
    end if
!
!  If we could not read an integer value, this line is useless or exhausted.
!  Loop again.
!
    line = ' '

  end do

  return
end
subroutine file_get_next_word ( iunit, word, line, line_num, ierror )

!*****************************************************************************80
!
!! FILE_GET_NEXT_WORD returns the next word and trailing context from a file.
!
!  Discussion:
!
!    The file should have been opened before calling this routine.
!    The file should contain ASCII text, which can be thought of as
!    words separated by one or more blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) IUNIT, the unit number associated with the file.
!
!    Output, character ( len = * ) WORD, the next word in the file.  If the
!    current line of the file is blank, or if the file has been exhausted, 
!    WORD will be set to ' '.
!
!    Input/output, character ( len = * ) LINE, the remaining text of the line 
!    that contains the information in WORD.  On each call, the next word
!    in LINE is extracted until LINE is empty, when it is refilled by
!    reading another line from the file.  Because LINE contains information
!    needed by this routine, it should not be altered by the user
!    between calls.
!
!    Input/output, integer ( kind = 4 ) LINE_NUM, the number of lines read from 
!    the file.  Before the first call to this routine, the user should set 
!    LINE_NUM to 0.
!
!    Output, integer ( kind = 4 ) IERROR, error flag. 
!    0, no error, another word was read, and returned in WORD.
!    1, end of file.  WORD and LINE were set to ' '.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) lenc
  character ( len = * ) line
  integer ( kind = 4 ) line_num
  character ( len = * ) word

  ierror = 0
!
!  If LINE_NUM is zero, then initialize LINE.
!
  if ( line_num <= 0 ) then
    line_num = 0
    line = ' '
  end if
!
!  If LINE is blank, try to read a new line from the file.
!
  if ( line == ' ' ) then

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = 1
      word = ' '
      line = ' '
      return
    end if

    line_num = line_num + 1

    if ( line == ' ' ) then
      word = ' '
      return
    end if

  end if
!
!  Extract the next word from LINE into WORD and return.
!
  lenc = len_trim ( line )
!
!  Find ILO, the index of the first nonblank in LINE.
!
  ilo = 1

  do while ( line(ilo:ilo) == ' ' ) 
    ilo = ilo + 1
  end do
!
!  Find IHI, the index of the last consecutive nonblank after the one at ILO.
!
  ihi = ilo

  do while ( ihi + 1 <= lenc )
    if ( line(ihi+1:ihi+1) == ' ' ) then
      exit
    end if
    ihi = ihi + 1
  end do
!
!  Set WORD.
!
  word = line(ilo:ihi)
!
!  Slide TEXT to the left.
!
  if ( ihi + 1 <= lenc ) then
    line = line(ihi+1:)
  else
    line = ' '
  end if

  return
end
subroutine file_insert ( input_unit, output_unit, line_num )

!*****************************************************************************80
!
!! FILE_INSERT copies the contents of an input file into an output file.
!
!  Discussion:
!
!    Both the input and output files should already be opened by the 
!    user.  The routine simply reads a line from the input file and 
!    writes it to the output file.  The input file is assumed to be 
!    a simple, sequential access text file, with records no longer 
!    than 256 characters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT_UNIT, the unit number associated 
!    with the input file.
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the unit number associated 
!    with the output file.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of records copied 
!    from the input file to the output file.
!
  implicit none

  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) output_unit

  line_num = 0
!
!  Make sure the file names aren't the same.
!
  if ( input_unit == output_unit ) then
    return
  end if

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      return
    end if

    line_num = line_num + 1
    write ( output_unit, '(a)' ) trim ( line )

  end do

  return
end
function file_is_open ( filename )

!*****************************************************************************80
!
!! FILE_IS_OPEN reports whether a file (specified by filename) is open.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
!    Output, logical FILE_IS_OPEN, is TRUE if the file is open.
!
  implicit none

  character ( len = * ) filename
  logical file_is_open

  inquire ( file = filename, opened = file_is_open )

  return
end
subroutine file_line_get ( filename, line_index, line )

!*****************************************************************************80
!
!! FILE_LINE_GET gets a particular line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
!    Input, integer ( kind = 4 ) LINE_INDEX, the index of the line to be read.
!
!    Output, character ( len = * ) LINE, the text of the line.
!
  implicit none

  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = * ) line
  integer ( kind = 4 ) line_index
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    line = ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_GET - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    return
  end if
!
!  Count the lines.
!
  do i = 1, line_index

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      line = ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_LINE_GET - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file.'
      return
    end if

  end do

  close ( unit = iunit )

  return
end
subroutine file_line_uniform ( filename, seed, line, line_index, line_num )

!*****************************************************************************80
!
!! FILE_LINE_UNIFORM returns a random line from a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    The algorithm used is interesting because it does not require
!    the number of lines in the file to be known in advance, and it
!    only reads the file once.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Tom Christiansen and Nathan Torkington,
!    "8.6: Picking a Random Line from a File",
!    Perl Cookbook, pages 284-285,
!    O'Reilly, 1999.
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, character ( len = * ) LINE, a random line from the file.
!
!    Output, integer ( kind = 4 ) LINE_INDEX, the index of the chosen line.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of lines in the file.
!
  implicit none

  character ( len = * ) filename
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = * ) line
  integer ( kind = 4 ) line_index
  integer ( kind = 4 ) line_num
  character ( len = 255 ) line_read
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  line_num = 0
  line_index = -1
  line = ' '
  line_read = ' '
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    return
  end if
!
!  Read the lines.
!
  do

    read ( iunit, '(a)', iostat = ios ) line_read

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

    r = r8_uniform_01 ( seed )

    if ( r * real ( line_num, kind = 8 ) <= 1.0D+00 ) then

      line = line_read
      line_index = line_num

    end if

  end do

  close ( unit = iunit )

  return
end
function file_line_width ( filename )

!*****************************************************************************80
!
!! FILE_LINE_WIDTH returns the length of the longest line in a file.
!
!  Discussion:
!
!    Whether or not this routine works properly depends on some 
!    system dependent matters.  For instance, the routine tries to
!    access the file a character at a time.  To do this, it seems
!    necessary to treat the file, presumably a text file, as though
!    it were an unformatted file, with direct access.  Direct access
!    requires a record length whose units are left undecided by the
!    standard.  We here assume that the unit is a byte.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file to be read.
!
!    Output, integer ( kind = 4 ) FILE_LINE_WIDTH, the length of the 
!    longest line.
!
  implicit none

  character ch
  integer ( kind = 4 ) ch_num
  integer ( kind = 4 ) file_line_width
  character ( len = * ) filename
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) line_width
  integer ( kind = 4 ) record
  integer ( kind = 4 ) value

  value = -1
!
!  Open the file.
!
!  The smallest amount of information we can write at a time is
!  1 word = 4 bytes = 32 bits.
!
  call get_unit ( file_unit )

  open ( unit = file_unit, file = filename, status = 'old', &
    form = 'unformatted', access = 'direct', recl = 1 )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_WIDTH - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    file_line_width = value
    close (  unit = file_unit )
    return
  end if

  ch_num = 0
  record = 0
  line_width = 0
  value = 0

  do

    record = record + 1
    read ( file_unit, rec = record, iostat = ios ) ch

    if ( ios /= 0 ) then
      exit
    end if

    ch_num = ch_num + 1
!
!  Character 10 is the LF character.
!  Character 13 is the CR character.
!
    if ( iachar ( ch ) == 10 .or. iachar ( ch ) == 13 ) then
      line_width = 0
    else
      line_width = line_width + 1
      value = max ( value, line_width )
    end if

  end do

  close ( unit = file_unit )

  file_line_width = value

  return
end
subroutine file_lines_uniform ( filename, seed, n, line, line_index, line_num )

!*****************************************************************************80
!
!! FILE_LINES_UNIFORM selects N random lines from a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    The algorithm used is interesting because it does not require
!    the number of lines in the file to be known in advance, and it
!    only reads the file once.
!
!    If the number of lines requested is more than the number of lines
!    in the file, then all the lines in the file will be selected, and
!    extra blank lines will be appended, with index -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Tom Christiansen and Nathan Torkington,
!    "8.6: Picking a Random Line from a File",
!    Perl Cookbook, pages 284-285,
!    O'Reilly, 1999.
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Input, integer ( kind = 4 ) N, the number of lines to be extracted.
!
!    Output, character ( len = * ) LINE(N), N random lines from the file.
!
!    Output,  integer ( kind = 4 ) LINE_INDEX(N), the indices of the lines.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of lines in the file.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = * ) line(n)
  integer ( kind = 4 ) line_index(n)
  integer ( kind = 4 ) line_num
  character ( len = 255 ) line_read
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  line_num = 0
  line_index(1:n) = -1
  do i = 1, n
    line(i) = ' '
  end do
  line_read = ' '
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINES_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    return
  end if
!
!  Read the lines.
!
  do

    read ( iunit, '(a)', iostat = ios ) line_read

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

    if ( line_num <= n ) then

      i = line_num
      line(i) = line_read
      line_index(i) = line_num

    else

      r = r8_uniform_01 ( seed )

      if ( r * real ( line_num, kind = 8 ) <= real ( n ) ) then

        i = i4_uniform ( 1, n, seed )
        line(i) = line_read
        line_index(i) = line_num

      end if

    end if

  end do

  close ( unit = iunit )

  return
end
subroutine file_merge ( filename_1, filename_2, filename_3, n1, n2, n3 )

!*****************************************************************************80
!
!! FILE_MERGE merges two sorted files into a third.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename_1, filename_2, the names of the
!    two input files to be merged.
!
!    Input, character ( len = * ) filename_3, the name of the output file to
!    be created.
!
!    Output, integer ( kind = 4 ) N1, N2, N3, the number of lines of text in the
!    two input files and the output file.
!
  implicit none

  character ( len = * ) filename_1
  character ( len = * ) filename_2
  character ( len = * ) filename_3
  integer ( kind = 4 ) file_unit_1
  integer ( kind = 4 ) file_unit_2
  integer ( kind = 4 ) file_unit_3
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  character ( len = 255 ) word1
  character ( len = 255 ) word2
  character ( len = 255 ) word3

  n1 = 0
  n2 = 0
  n3 = 0

  call get_unit ( file_unit_1 )

  open ( unit = file_unit_1, file = filename_1, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_MERGE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file #1:'
    write ( *, '(4x,a)' ) '"' // trim ( filename_1 ) // '".'
    return
  end if

  call get_unit ( file_unit_2 )

  open ( unit = file_unit_2, file = filename_2, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_MERGE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file #2:'
    write ( *, '(4x,a)' ) '"' // trim ( filename_2 ) // '".'
    return
  end if

  call get_unit ( file_unit_3 )

  open ( unit = file_unit_3, file = filename_3, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_MERGE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename_3 ) // '".'
    return
  end if

  word1 = ' '
  word2 = ' '

  do

    if ( word1 == ' ' ) then

      read ( file_unit_1, '(a)', iostat = ios ) word1

      if ( ios == 0 ) then
        call s_low ( word1 )
        n1 = n1 + 1
      end if

      if ( word1 == ' ' ) then
        word1 = '_END_'
        close ( unit = file_unit_1 )
      end if

    end if

    if ( word2 == ' ' ) then
 
      read ( file_unit_2, '(a)', iostat = ios ) word2

      if ( ios == 0 ) then
        call s_low ( word2 )
        n2 = n2 + 1
      end if

      if ( word2 == ' ' ) then
        word2 = '_END_'
        close ( unit = file_unit_2 )
      end if

    end if

    if ( word1 == '_END_' .and. word2 == '_END_' ) then

      exit

    else if ( word1 /= '_END_' .and. word2 == '_END_' ) then

      word3 = word1
      word1 = ' '

    else if ( word1 == '_END_' .and. word2 /= '_END_' ) then

      word3 = word2
      word2 = ' '

    else

      if ( llt ( word1, word2 ) ) then

        word3 = word1
        word1 = ' '

      else if ( word1 == word2 ) then

        word3 = word1
        word1 = ' '
        word2 = ' '

      else

        word3 = word2
        word2 = ' '

      end if

    end if

    write ( file_unit_3, '(a)' ) trim ( word3 )
    n3 = n3 + 1

  end do

  endfile ( unit = file_unit_3 )

  close ( unit = file_unit_3 )

  return
end
subroutine filename_append ( filename, append )

!*****************************************************************************80
!
!! FILENAME_APPEND appends a string to a filename, before the extension.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  
!
!    A file with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!    The idea is that the string in APPEND is to be appended to
!    the part of the filename that precedes the extension.
!
!    The intended purpose of this routine is to be able to easily
!    generate a filename that indicates its relation to another file.
!
!  Example:
!
!          Input             Output
!    ===================     =========
!    filename    APPEND     filename
!
!    bob.for      6          bob6.for
!    bob.bob.bob  JOB        bob.bobJOB.bob
!    bob          yak        bobyak
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILENAME, a file name.
!    On output, the file name has been modified.
!
!    Input, character ( len = * ) APPEND, the string to be appended.
!
  implicit none

  character ( len = * ) append
  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) len_append
  integer ( kind = 4 ) len_name

  call filename_ext_get ( filename, i, j )
!
!  If there is no extension, then simply slap APPEND on the end.
!
  if ( i == -1 ) then

    len_name = len_trim ( filename )
    filename(len_name+1:) = append
!
!  If there is an extension, then insert APPEND.
!
  else

    len_append = len_trim ( append )
    filename(i:) = append(1:len_append) // filename(i:j)

  end if

  return
end
subroutine filename_dec ( filename )

!*****************************************************************************80
!
!! FILENAME_DEC decrements a partially numeric file name.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be decreased by 1 on
!    each call.  If this number is all 0's on input, the output number
!    is all 9's.  Non-numeric letters of the name are unaffected.
!
!    If the name is empty, then the routine stops.
!
!    If the name contains no digits, the empty string is returned.
!
!  Example:
!
!      Input            Output
!      -----            ------
!      'a7to12.txt'     'a7to11.txt' (typical case.  Last digit decremented)
!      'a7to00.txt'     'a8to99.txt' (last digit decremented, with carry.)
!      'a0to00.txt'     'a9to99.txt' (wrap around)
!      'cat.txt'        ' '          (no digits in input name.)
!      ' '              STOP!        (error.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILENAME.
!    On input, a character string to be decremented.
!    On output, the decremented string.
!
  implicit none

  character c
  integer ( kind = 4 ) change
  integer ( kind = 4 ) digit
  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( filename )

  if ( lens <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILENAME_DEC - Fatal error!'
    write ( *, '(a)' ) '  The input filename is empty.'
    stop
  end if

  change = 0

  do i = lens, 1, -1

    c = filename(i:i)

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      change = change + 1

      digit = ichar ( c ) - 48
      digit = digit - 1

      if ( digit == -1 ) then
        digit = 9
      end if

      c = char ( digit + 48 )

      filename(i:i) = c

      if ( c /= '9' ) then
        return
      end if

    end if

  end do

  if ( change == 0 ) then
    filename = ' '
  end if

  return
end
subroutine filename_ext_get ( filename, i, j )

!*****************************************************************************80
!
!! FILENAME_EXT_GET determines the "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!    Blanks are unusual in filenames.  This routine ignores all
!    trailing blanks, but will treat initial or internal blanks
!    as regular characters acceptable in a file name.
!
!  Example:
!
!    filename   I  J
!
!    bob.for     4  7
!    N.B.C.D     6  7
!    Naomi.      6  6
!    Arthur     -1 -1 
!    .com        1  1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, a file name to be examined.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the first and 
!    last characters in the file extension.  
!    If no period occurs in filename, then
!      I = J = -1;
!    Otherwise,
!      I is the position of the LAST period in filename, and J is the
!      position of the last nonblank character following the period.
!
  implicit none

  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ch_index_last

  i = ch_index_last ( filename, '.' )

  if ( i == -1 ) then

    j = -1

  else

    j = len_trim ( filename )

  end if

  return
end
subroutine filename_ext_swap ( filename, ext )

!*****************************************************************************80
!
!! FILENAME_EXT_SWAP replaces the current "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!  Example:
!
!          Input           Output
!    ================     =========
!    filename    EXT     filename
!
!    bob.for      obj     bob.obj
!    bob.bob.bob  txt     bob.bob.txt
!    bob          yak     bob.yak
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) filename, a file name.
!    On output, the extension of the file has been changed.
!
!    Input, character ( len = * ) EXT, the extension to be used on the output
!    copy of filename, replacing the current extension if any.
!
  implicit none

  character ( len = * ) ext
  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) len_max
  integer ( kind = 4 ) len_name

  len_max = len ( filename )
  len_name = len_trim ( filename )

  call filename_ext_get ( filename, i, j )

  if ( i == -1 ) then

    if ( len_max < len_name + 1 ) then
      return
    end if

    len_name = len_name + 1
    filename(len_name:len_name) = '.'
    i = len_name + 1

  else

    i = i + 1
    filename(i:j) = ' '

  end if

  filename(i:) = ext

  return
end
subroutine filename_inc ( filename )

!*****************************************************************************80
!
!! FILENAME_INC increments a partially numeric filename.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected.
!
!    If the name is empty, then the routine stops.
!
!    If the name contains no digits, the empty string is returned.
!
!  Example:
!
!      Input            Output
!      -----            ------
!      'a7to11.txt'     'a7to12.txt'
!      'a7to99.txt'     'a8to00.txt'
!      'a9to99.txt'     'a0to00.txt'
!      'cat.txt'        ' '
!      ' '              STOP!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILENAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  integer ( kind = 4 ) change
  integer ( kind = 4 ) digit
  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( filename )

  if ( lens <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILENAME_INC - Fatal error!'
    write ( *, '(a)' ) '  The input string is empty.'
    stop
  end if

  change = 0

  do i = lens, 1, -1

    c = filename(i:i)

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      change = change + 1

      digit = ichar ( c ) - 48
      digit = digit + 1

      if ( digit == 10 ) then
        digit = 0
      end if

      c = char ( digit + 48 )

      filename(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do
!
!  No digits were found.  Return blank.
!
  if ( change == 0 ) then
    filename = ' '
    return
  end if

  return
end
subroutine filename_inc_nowrap ( filename )

!*****************************************************************************80
!
!! FILENAME_INC_NOWRAP increments a partially numeric filename.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  Non-numeric letters of the name are unaffected.
!
!    If the (nonempty) name contains no digits, or all the digits are
!    9, then the empty string is returned.
!
!    If the empty string is input, the routine stops.
!
!  Example:
!
!      Input            Output
!      -----            ------
!      'a7to11.txt'     'a7to12.txt'
!      'a7to99.txt'     'a8to00.txt'
!      'a8to99.txt'     'a9to00.txt'
!      'a9to99.txt'     ' '
!      'cat.txt'        ' '
!      ' '              STOP!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILENAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  integer ( kind = 4 ) carry
  integer ( kind = 4 ) change
  integer ( kind = 4 ) digit
  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( filename )

  if ( lens <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILENAME_INC_NOWRAP - Fatal error!'
    write ( *, '(a)' ) '  The input string is empty.'
    stop
  end if

  change = 0
  carry = 0

  do i = lens, 1, -1

    c = filename(i:i)

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      change = change + 1
      carry = 0

      digit = ichar ( c ) - 48
      digit = digit + 1

      if ( digit == 10 ) then
        digit = 0
        carry = 1
      end if

      c = char ( digit + 48 )

      filename(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do
!
!  Unsatisfied carry.  The input digits were all 9.  Return blank.
!
  if ( carry == 1 ) then
    filename = ' '
    return
  end if
!
!  No digits were found.  Return blank.
!
  if ( change == 0 ) then
    filename = ' '
    return
  end if

  return
end
subroutine file_para_count ( filename, para_num )

!*****************************************************************************80
!
!! FILE_PARA_COUNT counts the number of paragraphs in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.  A paragraph is
!    a sequence of nonblank lines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
!    Output, integer ( kind = 4 ) PARA_NUM, the number of paragraphs found in 
!    the file.
!
  implicit none

  character ( len = * ) filename
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lenc_old
  character ( len = 255 ) line
  integer ( kind = 4 ) para_num

  para_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    para_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_PARA_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    return
  end if
!
!  Count the paragraphs.
!
  lenc = 0

  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    lenc_old = lenc
    lenc = len_trim ( line )

    if ( 0 < lenc .and. lenc_old <= 0 ) then
      para_num = para_num + 1
    end if

  end do

  close ( unit = iunit )

  return
end
subroutine file_paren_check ( filename )

!*****************************************************************************80
!
!! FILE_PAREN_CHECK checks a file for generalized parenthesis errors.
!
!  Discussion:
!
!    The check made is that the current number of left parentheses read must
!    always be at least as great as the number of right parentheses read.
!    Moreover, when we reach the end of the file, the numbers must be equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
  implicit none

  character ( len = * ) filename
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  integer ( kind = 4 ) line_len
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) sum_p

  sum_p = 0
!
!  Open the file.
!
  call get_unit ( file_unit )

  open ( unit = file_unit, file = filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_PAREN_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    return
  end if

  line_num = 0

  do

    read ( file_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

    line_len = len_trim ( line )

    do i = 1, line_len

      if ( line(i:i) == '(' ) then

        sum_p = sum_p + 1

      else if ( line(i:i) == ')' ) then

        sum_p = sum_p - 1

        if ( sum_p < 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FILE_PAREN_CHECK - Warning!'
          write ( *, '(a)' ) '  Parenthesis error in the file:'
          write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
          write ( *, '(a,i8)' ) &
            '  An illegal right parenthesis occurs on line', line_num
          write ( *, '(a)' ) '    ' // trim ( line )
          close ( unit = file_unit )
          return
        end if

      end if

    end do

  end do

  close ( unit = file_unit )

  if ( sum_p /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_PAREN_CHECK - Warning!'
    write ( *, '(a)' ) '  Parenthesis error in the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    write ( *, '(a,i8)' ) '  Number of missing right parentheses: ', sum_p
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_PAREN_CHECK - Note:'
    write ( *, '(a)' ) '  Parenthesis checks passed for file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
  end if

  return
end
subroutine file_print ( filename )

!*****************************************************************************80
!
!! FILE_PRINT prints the contents of a text file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
  implicit none

  character ( len = * ) filename
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 255 ) line
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_PRINT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    return
  end if

  do

    read ( iunit, '(a)', iostat = ios ) line 

    if ( ios /= 0 ) then
      exit
    end if

    write ( *, '(a)' ) trim ( line )
 
  end do

  close ( unit = iunit )

  return
end
subroutine file_rename ( filename_old, filename_new )

!*****************************************************************************80
!
!! FILE_RENAME renames a file.
!
!  Discussion:
!
!    Actually, this routine copies the file, and deletes the original.
!    But to the user, it should look like a rename, just a little slower.
!
!    If a file already exists with the new name, it is deleted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename_OLD, the name of the original file.
!
!    Output, character ( len = * ) filename_NEW, the name of the new file.
!
  implicit none

  logical file_exist
  logical file_is_open
  character ( len = * ) filename_new
  character ( len = * ) filename_old
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical s_eqi
!
!  Does the old file exist?
!
  if ( .not. file_exist ( filename_old ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RENAME - Error!'
    write ( *, '(a)' ) '  The original file to be renamed does not exist.'
    return
  end if
!
!  Is the old file open?
!
  if ( file_is_open ( filename_old ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RENAME - Error!'
    write ( *, '(a)' ) '  The original file is open.'
    write ( *, '(a)' ) '  It must be closed before it can be renamed.'
    return
  end if
!
!  Does old file name = new file name?
!
  if ( s_eqi ( filename_new, filename_old ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RENAME: Warning!'
    write ( *, '(a)' ) '  The old and new file names are the same.'
    write ( *, '(a)' ) '  I suppose this means there is nothing to do.'
    return
  end if
!
!  Does the new file exist?
!
  if ( file_exist ( filename_new ) ) then
!
!  Is the new file open?
!
    if ( file_is_open ( filename_new ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_RENAME - Error!'
      write ( *, '(a)' ) '  The new file is already open.'
      write ( *, '(a)' ) '  It must be closed before it can be overwritten.'
      return
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RENAME:'
    write ( *, '(a)' ) '  Deleting pre-existing file with new name.'
    call file_delete ( filename_new )

  end if
!
!  Copy old into new.
!
  call file_copy ( filename_old, filename_new, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RENAME: Warning!'
    write ( *, '(a)' ) '  Could not copy the old file!'
    return
  end if
!
!  Delete the old file.
!
  call file_delete ( filename_old )

  return
end
subroutine file_reverse_columns ( input_filename, output_filename )

!*****************************************************************************80
!
!! FILE_REVERSE_COLUMNS makes a copy of a file with each lines reversed.
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
!      elat eht si sihT
!      sgip elttil eerht fo
!      .sliat rieht dna
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
!    Input, character ( len = * ) INPUT_filename, the name of the file to 
!    be reversed.
!
!    Input, character ( len = * ) OUTPUT_filename, the name of the file to be
!    created, contained a copy of the input file with the columns reversed.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
!
!  Open the input file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_REVERSE_COLUMNS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file:'
    write ( *, '(4x,a)' ) '"' // trim ( input_filename ) // '".'
    return
  end if
!
!  Open the output file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_REVERSE_COLUMNS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file:'
    write ( *, '(4x,a)' ) '"' // trim ( output_filename ) // '".'
    return
  end if

  do 

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    call s_reverse ( line )

    write ( output_unit, '(a)' ) trim ( line )

  end do

  close ( unit = input_unit )

  endfile ( unit = output_unit )
  close ( unit = output_unit )

  return
end
subroutine file_reverse_rows ( input_filename, output_filename )

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
!    Input, character ( len = * ) INPUT_filename, the name of the file to 
!    be reversed.
!
!    Input, character ( len = * ) OUTPUT_filename, the name of the file to be
!    created, contained a reversed copy of the input file.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) input_count
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
!
!  Open the input file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_REVERSE_ROWS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file:'
    write ( *, '(4x,a)' ) '"' // trim ( input_filename ) // '".'
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

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_REVERSE_ROWS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file:'
    write ( *, '(4x,a)' ) '"' // trim ( output_filename ) // '".'
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

  return
end
subroutine file_rot13 ( input_filename, output_filename )

!*****************************************************************************80
!
!! FILE_ROT13 makes a ROT13-encoded copy of a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_filename, the name of the file to 
!    be reversed.
!
!    Input, character ( len = * ) OUTPUT_filename, the name of the file to be
!    created, contained a ROT13-encoded copy of the input file.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
!
!  Open the input file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROT13 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file:'
    write ( *, '(4x,a)' ) '"' // trim ( input_filename ) // '".'
    return
  end if
!
!  Open the output file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROT13 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file:'
    write ( *, '(4x,a)' ) '"' // trim ( output_filename ) // '".'
    return
  end if

  do 

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    call s_to_rot13 ( line )

    write ( output_unit, '(a)' ) trim ( line )

  end do

  close ( unit = input_unit )

  endfile ( unit = output_unit )
  close ( unit = output_unit )

  return
end
subroutine file_row_count ( filename, line_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of rows in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Blank lines and comment lines, which begin with '#', are not counted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of lines found in the 
!    file.  If the file could not be opened, then LINE_NUM is returned as -1.
!
  implicit none

  character ( len = * ) filename
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 255 ) line
  integer ( kind = 4 ) line_num
  logical, parameter :: verbose = .false.

  line_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then

    line_num = -1

    if ( verbose ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
      write ( *, '(a)' ) '  Could not open the file:'
      write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    end if

    return

  end if
!
!  Count the lines.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    line_num = line_num + 1

  end do

  close ( unit = iunit )

  return
end
subroutine file_sequence_delete ( filename, delete_count )

!*****************************************************************************80
!
!! FILE_SEQUENCE_DELETE deletes a file sequence.
!
!  Discussion:
!
!    We suppose the user has a set of files whose names differ only
!    in some numeric tag that is sequentially increasing, as, perhaps,
!    "file001.txt", "file002.txt" through "file137.txt", say.
!
!    The user specifies filename as the name of the first file in the
!    sequence.  This function deletes that file, generates the next
!    name in the sequence, and, if a file with that name exists, it
!    deletes it as well.  The process continues until a file name is
!    reached for which there is no existing file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the first file
!    in the sequence.
!
!    Output, integer ( kind = 4 ) DELETE_COUNT, the number of files deleted.
!
  implicit none

  integer ( kind = 4 ) delete_count
  logical file_exist
  character ( len = * ) filename
  character ( len = 255 ) filename2

  delete_count = 0
  filename2 = filename

  do while ( file_exist ( filename2 ) )

    call file_delete ( filename2 )

    delete_count = delete_count + 1

    call filename_inc ( filename2 )
     
  end do

  return
end
subroutine file_sequence_size ( filename, file_dim, file_num )

!*****************************************************************************80
!
!! FILE_SEQUENCE_SIZE sizes a file sequence.
!
!  Discussion:
!
!    We suppose the user has a set of files whose names differ only
!    in some numeric tag that is sequentially increasing, as, perhaps,
!    "file001.txt", "file002.txt" through "file137.txt", say.
!
!    The user specifies the name of the first file in the sequence.
!    This function determines the number of files in the sequence, 
!    and makes a guess for the "dimension" of the files, that is, the number
!    of numeric data items.
!
!    Note that the function only checks the dimension of the data in
!    the first file.  It is up to the user to determine whether this
!    dimension is used for every file in the sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the first file in 
!    the sequence.
!
!    Output, integer ( kind = 4 ) FILE_DIM, the dimension of the data 
!    in one file.
!
!    Output, integer ( kind = 4 ) FILE_NUM, the number of files.
!
  implicit none

  integer ( kind = 4 ) file_dim
  logical file_exist
  character ( len = * ) filename
  character ( len = 255 ) filename2
  integer ( kind = 4 ) file_num
  integer ( kind = 4 ) file_unit

  file_num = 0
  file_dim = 0
  filename2 = filename

  do

    if ( .not. file_exist ( filename2 ) ) then
      exit
    end if

    file_num = file_num + 1

    if ( file_num == 1 ) then

      call file_word_count ( filename2, file_dim )

    end if

    call filename_inc ( filename2 )

  end do

  return
end
function file_tag_check ( filename, left, right )

!*****************************************************************************80
!
!! FILE_TAG_CHECK checks a file for generalized parenthesis errors.
!
!  Discussion:
!
!    The check made is that the current number of left "parentheses" read must
!    always be at least as great as the number of right "parentheses" read.
!    Moreover, when we reach the end of the file, the numbers must be equal.
!
!    Typical examples of left and right parentheses might be:
!
!    (), [], {}, <>, <P> </P>, 'do' 'end do'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
!    Input, character ( len = * ) LEFT, RIGHT, the left and right
!    parentheses marks.
!
!    Output, logical FILE_TAG_CHECK, is true if the file passed the check.
!
  implicit none

  character ( len = * ) filename
  logical file_tag_check
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) ios
  character ( len = * ) left
  integer ( kind = 4 ) left_len
  integer ( kind = 4 ) left_pos
  character ( len = 255 ) line
  character ( len = 255 ) line_copy
  integer ( kind = 4 ) line_len
  integer ( kind = 4 ) line_num
  character ( len = * ) right
  integer ( kind = 4 ) right_len
  integer ( kind = 4 ) right_pos
  integer ( kind = 4 ) sum_p

  sum_p = 0
  left_len = len ( left )
  right_len = len ( right )
!
!  Open the file.
!
  call get_unit ( file_unit )

  open ( unit = file_unit, file = filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_TAG_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    file_tag_check = .false.
    return
  end if

  line_num = 0

  do

    read ( file_unit, '(a)', iostat = ios ) line
    line_copy = line

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

    line_len = len_trim ( line )

    do

      left_pos = index ( line, left )
      if ( left_pos == 0 ) then
        left_pos = line_len + 1
      end if

      right_pos = index ( line, right )
      if ( right_pos == 0 ) then
        right_pos = line_len + 1
      end if

      if ( left_pos < right_pos ) then

        sum_p = sum_p + 1
        line = adjustl ( line(left_pos+left_len:) )
        line_len = len_trim ( line )

      else if ( right_pos < left_pos ) then

        sum_p = sum_p - 1
        line = adjustl ( line(right_pos+right_len:) )
        line_len = len_trim ( line )

        if ( sum_p < 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FILE_TAG_CHECK - Warning!'
          write ( *, '(a)' ) '  Tag error in the file:'
          write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
          write ( *, '(a,i8)' ) '  An illegal right tag occurs on line', &
            line_num
          write ( *, '(a)' ) '    ' // trim ( line_copy )
          close ( unit = file_unit )
          file_tag_check = .false.
          return
        end if

      else

        exit

      end if

    end do

  end do

  close ( unit = file_unit )

  if ( sum_p /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_TAG_CHECK - Warning!'
    write ( *, '(a)' ) '  Tag error in the file:'
    write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
    write ( *, '(a,i8)' ) '  Number of missing right tags: ', sum_p
    file_tag_check = .false.
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_TAG_CHECK - Note:'
    write ( *, '(a)' ) '  Tag checks passed for file:'
    write ( *, '(a)' ) '    "' // trim ( filename ) // '".'
    file_tag_check = .true.
  end if

  return
end
subroutine file_unique_lines ( input_filename, output_filename, &
  input_line_num, output_line_num )

!*****************************************************************************80
!
!! FILE_UNIQUE_LINES makes a copy of the unique lines of a sorted file.
!
!  Discussion:
!
!    Actually, the input file doesn't have to be sorted.  The routine
!    simply reads each line of the input file, and writes it to the
!    output file if it is distinct from the previous input line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_filename, the name of the input file.
!
!    Input, character ( len = * ) OUTPUT_filename, the name of the output file.
!
!    Output, integer ( kind = 4 ) INPUT_LINE_NUM, the number of lines in the 
!    input file.
!
!    Output, integer ( kind = 4 ) OUTPUT_LINE_NUM, the number of lines in the 
!    output file.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) input_line_num
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  character ( len = 255 ) line_old
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_line_num
  integer ( kind = 4 ) output_unit

  input_line_num = 0
  output_line_num = 0
  line_old = ' '

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_UNIQUE_LINES - Error!'
    write ( *, '(a)' ) '  Could not open the input file:'
    write ( *, '(4x,a)' ) '"' // trim ( input_filename ) // '".'
    return
  end if

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_UNIQUE_LINES - Error!'
    write ( *, '(a)' ) '  Could not open the output file:'
    write ( *, '(4x,a)' ) '"' // trim ( output_filename ) // '".'
    close ( unit = input_unit )
    return
  end if

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    input_line_num = input_line_num + 1

    if ( input_line_num == 1 .or. line /= line_old ) then

      write ( output_unit, '(a)' ) trim ( line )
      output_line_num = output_line_num + 1
      line_old = line

    end if

  end do

  close ( unit = input_unit )

  endfile ( unit = output_unit )
  close ( unit = output_unit )

  return
end
subroutine file_word_count ( filename, word_num )

!*****************************************************************************80
!
!! FILE_WORD_COUNT counts the number of words in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) filename, the name of the file.
!
!    Output, integer ( kind = 4 ) WORD_NUM, the number of words found in the file.
!
  implicit none

  character ( len = * ) filename
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 255 ) line
  integer ( kind = 4 ) nplus
  integer ( kind = 4 ) word_num

  word_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filename, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    word_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_WORD_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( filename ) // '".'
    return
  end if
!
!  Read the lines.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    call s_word_count ( line, nplus )

    word_num = word_num + nplus

  end do

  close ( unit = iunit )

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
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
function len_nonnull ( s )

!*****************************************************************************80
!
!! LEN_NONNULL returns the length of a string up to the last non-null character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to measure.
!
!    Output, integer ( kind = 4 ) LEN_NONNULL, the length of the string, up to 
!    the last non-null character.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len_nonnull
  integer ( kind = 4 ) len_s
  character, parameter :: NULL = char ( 0 )
  character ( len = * ) s

  len_s = len ( s )
 
  do i = len_s, 1, -1
    if ( s(i:i) /= NULL ) then
      len_nonnull = i
      return
    end if
  end do

  len_nonnull = 0

  return
end
subroutine number_inc ( s )

!*****************************************************************************80
!
!! NUMBER_INC increments the integer represented by a string.
!
!  Example:
!
!    Input      Output
!    -----      ------
!    '17'       '18'
!    'cat3'     'cat4'
!    '2for9'    '3for0'
!    '99thump'  '00thump'
!
!  Discussion:
!
!    If the string contains characters that are not digits, they will
!    simply be ignored.  If the integer is all 9's on input, then
!    the output will be all 0's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, a string representing an integer.
!
  implicit none

  logical ch_is_digit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  character ( len = * ) s

  lens = len_trim ( s )

  do i = lens, 1, -1

    if ( ch_is_digit ( s(i:i) ) ) then

      call digit_inc ( s(i:i) )

      if ( s(i:i) /= '0' ) then
        return
      end if

    end if

  end do

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
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
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
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

  integer ( kind = 4 ) i
  character ( len = * ) s

  do i = 1, len_trim ( s )
    call ch_low ( s(i:i) )
  end do

  return
end
subroutine s_reverse ( s )

!*****************************************************************************80
!
!! S_REVERSE reverses the characters in a string.
!
!  Example:
!
!    Input        Output
!
!    ' Cat'       'taC '
!    'Goo gol  '  'log ooG  '
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to reverse.
!    Trailing blanks are ignored.
!
  implicit none

  character ch
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len_trim ( s )
 
  do i = 1, s_length / 2
    j = s_length + 1 - i
    ch     = s(i:i)
    s(i:i) = s(j:j)
    s(j:j) = ch
  end do
 
  return
end
subroutine s_to_i4 ( s, value, ierror, length )

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
!    13 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) VALUE, the integer value read from the string.
!    If the string is blank, then VALUE will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters used. 
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * ) s
  integer ( kind = 4 ) state
  integer ( kind = 4 ) value

  value = 0
  ierror = 0
  length = 0

  state = 0
  isgn = 1

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  STATE = 0, haven't read anything.
!
    if ( state == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        state = 1
        isgn = -1
      else if ( c == '+' ) then
        state = 1
        isgn = +1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = ichar ( c ) - ichar ( '0' )
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 2, have read at least one digit, expecting more.
!
    else if ( state == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + ichar ( c ) - ichar ( '0' )

      else

        value = isgn * value
        ierror = 0
        length = i - 1
        return

      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( state == 2 ) then

    value = isgn * value
    ierror = 0
    length = len_trim ( s )

  else

    value = 0
    ierror = 1
    length = 0

  end if

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the real value that was read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    * 0, no errors occurred.
!    * 1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 .or. nchar <= lchar+1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lchar
  real ( kind = 8 ) rvec(n)
  character ( len = * ) s

  i = 0

  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

  end do

  return
end
subroutine s_to_rot13 ( s )

!*****************************************************************************80
!
!! S_TO_ROT13 "rotates" the alphabetical characters in a string by 13 positions.
!
!  Discussion:
!
!    Two applications of the routine will return the original string.
!
!  Example:
!
!    Input:                      Output:
!
!    abcdefghijklmnopqrstuvwxyz  nopqrstuvwxyzabcdefghijklm
!    Cher                        Pure
!    James Thurston Howell       Wnzrf Guhefgba Ubjryy
!    0123456789                  5678901234
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!   07 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, a string to be "rotated".
!
  implicit none

  character ch_to_rot13
  integer ( kind = 4 ) i
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len_trim ( s )

  do i = 1, s_length
    s(i:i) = ch_to_rot13 ( s(i:i) )
  end do

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
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) WORD_NUM, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  character ( len = * ) s
  integer ( kind = 4 ) word_num

  word_num = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

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
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
