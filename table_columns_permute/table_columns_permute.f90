program main

!*****************************************************************************80
!
!! MAIN is the main program for TABLE_COLUMNS_PERMUTE.
!
!  Discussion:
!
!    TABLE_COLUMNS_PERMUTE reads a table file, which it regards as
!    a table of M rows and N columns of "words".
!
!    A word begins at the first position in a line of text, and
!    at the first blank following a non-blank character.
!
!    A word terminates at the last non-blank character.
!
!    Thus, words will include the blank spaces in front of them.
!    (We need this, in order to be able to hope that the permuted
!    columns will still properly line up.)
!
!    The user supplies a permutation to be applied to the columns.
!    That is, the user types in the new column for the data that
!    is currently in column 1, followed by the new column for the
!    data in column 2, and so on.
!
!    The program permutes the columns and writes out a copy of
!    the file.
!
!    Lines of the table that begin with the comment character "#"
!    are not affected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Usage:
!
!    table_columns_permute input_filename
!
  implicit none

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) column_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) :: ierror = 0
  character ( len = 255 ) :: input_filename = ' '
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  integer ( kind = 4 ) line_num
  character ( len = 255 ) line2
  character ( len = 255 ) :: output_filename = ' '
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ), allocatable, dimension ( : ) :: perm
  integer ( kind = 4 ) row_num

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TABLE_COLUMNS_PERMUTE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read a dataset X(1:M,1:N) of M rows and N columns,'
  write ( *, '(a)' ) '  get a column permutation vector PERM(1:N),'
  write ( *, '(a)' ) '  write the permuted data X(1:M,PERM(1:N)) to a file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data is regarded as "words", and an effort'
  write ( *, '(a)' ) '  is made to process the data in a way that preserves'
  write ( *, '(a)' ) '  column alignment.  This works best if the data in'
  write ( *, '(a)' ) '  column 1 always begins with at least one blank!'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, input_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TABLE_COLUMNS_PERMUTE:'
    write ( *, '(a)' ) '  Please enter the name of the input file.'

    read ( *, '(a)' ) input_filename

  end if
!
!  Get the number of rows and columns of data.
!
  call file_column_count ( input_filename, column_num )

  call file_row_count ( input_filename, row_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header of "' // trim ( input_filename ) //'".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of columns = ', column_num
  write ( *, '(a,i8)' ) '  Number of rows    = ', row_num
!
!  Get the permutation vector.
!
  allocate ( perm(1:column_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the permutation vector:'

  read ( *, * ) perm(1:column_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The permutation vector:'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,30i3)' ) perm(1:column_num)
!
!  Read lines, permute them, write them out.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TABLE_COLUMNS_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
    stop
  end if

  output_filename = input_filename
  call file_name_ext_swap ( output_filename, 'permuted.txt' )

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TABLE_COLUMNS_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file: ' // &
      trim ( output_filename )
    stop
  end if

  line_num = 0

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then

      write ( output_unit, '(a)', iostat = ios ) trim ( line )

    else

      call s_word_permute ( line, column_num, perm, line2 )
      write ( output_unit, '(a)', iostat = ios ) trim ( line2 )

    end if

  end do

  close ( unit = input_unit )
  close ( unit = output_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data from "' // trim ( input_filename ) //'".'
  write ( *, '(a)' ) '  Wrote the data to  "' // trim ( output_filename ) //'".'
!
!  Terminate
!
  deallocate ( perm )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TABLE_COLUMNS_PERMUTE'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ch_cap ( c )

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
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine file_column_count ( input_filename, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some
!    comment lines,
!    which have a "#" in column 1.
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
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer ( kind = 4 ) column_num
  logical got_one
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 256 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( input_filename )
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = ios ) line

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

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = ios ) line

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

  close ( unit = input_unit )

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
subroutine file_name_ext_get ( file_name, i, j )

!*****************************************************************************80
!
!! FILE_NAME_EXT_GET determines the "extension" of a file name.
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
!    FILE_NAME   I  J
!
!    bob.for     4  7
!    N.B.C.D     6  7
!    Naomi.      6  6
!    Arthur      0  0
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
!    Input, character ( len = * ) FILE_NAME, a file name to be examined.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the first and last characters
!    in the file extension.
!
!    If no period occurs in FILE_NAME, then
!      I = J = 0;
!    Otherwise,
!      I is the position of the LAST period in FILE_NAME, and J is the
!      position of the last nonblank character following the period.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s_index_last

  i = s_index_last ( file_name, '.' )

  if ( i /= 0 ) then

    j = len_trim ( file_name )

  else

    j = 0

  end if

  return
end
subroutine file_name_ext_swap ( file_name, ext )

!*****************************************************************************80
!
!! FILE_NAME_EXT_SWAP replaces the current "extension" of a file name.
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
!    FILE_NAME    EXT     FILE_NAME
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
!    Input/output, character ( len = * ) FILE_NAME, a file name.
!    On output, the extension of the file has been changed.
!
!    Input, character ( len = * ) EXT, the extension to be used on the output
!    copy of FILE_NAME, replacing the current extension if any.
!
  implicit none

  character ( len = * ) ext
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) len_max
  integer ( kind = 4 ) len_name

  len_max = len ( file_name )
  len_name = len_trim ( file_name )

  call file_name_ext_get ( file_name, i, j )

  if ( i == 0 ) then

    if ( len_max < len_name + 1 ) then
      return
    end if

    len_name = len_name + 1
    file_name(len_name:len_name) = '.'
    i = len_name + 1

  else

    i = i + 1
    file_name(i:j) = ' '

  end if

  file_name(i:) = ext

  return
end
subroutine file_row_count ( input_filename, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 100 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    row_num = row_num + 1

  end do

  close ( unit = input_unit )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
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
subroutine perm_inverse3 ( n, perm, perm_inv )

!*****************************************************************************80
!
!! PERM_INVERSE3 produces the inverse of a given permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items permuted.
!
!    Input, integer ( kind = 4 ) PERM(N), a permutation.
!
!    Output, integer ( kind = 4 ) PERM_INV(N), the inverse permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)

  do i = 1, n
    perm_inv(perm(i)) = i
  end do

  return
end
subroutine s_copy ( s1, s2 )

!*****************************************************************************80
!
!! S_COPY copies one string into another.
!
!  Discussion:
!
!    If S1 is shorter than S2, the rest of S2 is blank.
!    If S1 is longer than S2, then the excess information is lost.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the string to be copied.
!
!    Output, character ( len = * ) S2, the copy.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2

  s2(1:min(len(s1),len(s2))) = s1(1:min(len(s1),len(s2)))
  s2(len(s1)+1:len(s2)) = ' '

  return
end
function s_index_last ( s, sub )

!*****************************************************************************80
!
!! S_INDEX_LAST finds the LAST occurrence of a given substring.
!
!  Discussion:
!
!    It returns the location in the string at which the substring SUB is
!    first found, or 0 if the substring does not occur at all.
!
!    The routine is also trailing blank insensitive.  This is very
!    important for those cases where you have stored information in
!    larger variables.  If S is of length 80, and SUB is of
!    length 80, then if S = 'FRED' and SUB = 'RED', a match would
!    not be reported by the standard FORTRAN INDEX, because it treats
!    both variables as being 80 characters long!  This routine assumes that
!    trailing blanks represent garbage!
!
!    This means that this routine cannot be used to find, say, the last
!    occurrence of a substring 'A ', since it assumes the blank space
!    was not specified by the user, but is, rather, padding by the
!    system.  However, as a special case, this routine can properly handle
!    the case where either S or SUB is all blanks.
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
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character ( len = * ) SUB, the substring to search for.
!
!    Output, integer ( kind = 4 ) S_INDEX_LAST.  0 if SUB does not occur in
!    the string.  Otherwise S_INDEX_LAST = I, where S(I:I+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the last place
!    this happens.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) llen1
  integer ( kind = 4 ) llen2
  character ( len = * ) s
  integer ( kind = 4 ) s_index_last
  character ( len = * ) sub

  s_index_last = 0

  llen1 = len_trim ( s )
  llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN
!
  if ( llen1 == 0 ) then
    llen1 = len ( s )
  end if

  if ( llen2 == 0 ) then
    llen2 = len ( sub )
  end if

  if ( llen1 < llen2 ) then
    return
  end if

  do j = 1, llen1+1-llen2

    i = llen1 + 2 - llen2 - j

    if ( s(i:i+llen2-1) == sub ) then
      s_index_last = i
      return
    end if

  end do

  return
end
subroutine s_word_count ( s, nword )

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
!    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) nword
  character ( len = * ) s

  nword = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      nword = nword + 1
      blank = .false.
    end if

  end do

  return
end
subroutine s_word_permute ( s1, n, perm, s2 )

!*****************************************************************************80
!
!! S_WORD_PERMUTE permutes the words in a string.
!
!  Discussion:
!
!    A word is a blank-delimited sequence of characters.
!
!    The string is assumed to contain N "words".  If more words are
!    in the string, their position is not affected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, a line of text.
!
!    Input, integer ( kind = 4 ) N, the number of words to permute.
!
!    Input, integer ( kind = 4 ) PERM(N), the permutation.  PERM(1) is the new
!    location of the item whose original location was 1.
!
!    Output, character ( len = * ) S2, a copy of S1 with the
!    first N words permuted.
!
  implicit none

  integer ( kind = 4 ) n

  character c1
  character c2
  integer ( kind = 4 ) index1
  integer ( kind = 4 ) index2
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)
  character ( len = * ) s1
  integer ( kind = 4 ) s1_length
  integer ( kind = 4 ) s1_pos
  integer ( kind = 4 ) s1_word_index(n)
  integer ( kind = 4 ) s1_word_length(n)
  character ( len = * ) s2
  integer ( kind = 4 ) s2_pos
  integer ( kind = 4 ) word_length
!
!  Set up word position and length vectors.
!
  s1_length = len ( s1 )

  s1_word_length(1:n) = 0
  s1_word_index(1:n) = 0

  index1 = 0
  c2 = ' '

  do s1_pos = 1, s1_length

    c1 = c2
    c2 = s1(s1_pos:s1_pos)

    if ( s1_pos == 1 .or. ( c1 /= ' ' .and. c2 == ' ' ) ) then

      if ( n <= index1 ) then
        exit
      end if

      index1 = index1 + 1

      s1_word_index(index1) = s1_pos

    end if

    s1_word_length(index1) = s1_word_length(index1) + 1

  end do
!
!  Invert the permutation.
!
  call perm_inverse3 ( n, perm, perm_inv )
!
!  Copy S1 into S2, so we get any trailing information.
!
  call s_copy ( s1, s2 )
!
!  Copy the first N words of S1 into S2 in permuted order.
!
  s2_pos = 1

  do index2 = 1, n

    index1 = perm_inv(index2)

    s1_pos = s1_word_index(index1)

    word_length = s1_word_length(index1)

    s2(s2_pos:s2_pos+word_length-1) = s1(s1_pos:s1_pos+word_length-1)

    s2_pos = s2_pos + word_length

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
