program main

!*****************************************************************************80
!
!! MAIN is the main program for the SUBANAGRAM program.
!
!  Discussion:
!
!    An ANAGRAM of a master word is a new word formed by using all the letters
!    of the master word.
!
!    A SUBANAGRAM of a master word is a new word formed by using some of the
!    letters of the master word.
!
!    TEXTROPOLIS is an iPhone game that seeks as many subanagrams of
!    certain city names as possible, although they are required to
!    have at least 4 letters.
!
!  Dedication:
!
!    This program is dedicated to Vivian Benton.
!
!  Usage:
!
!    subanagram length word
!
!    where
!
!    * length is the minimum length of a subanagram;
!    * word is the master word (no spaces).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 )  arg
  integer   ( kind = 4 )  arg_num
  integer   ( kind = 4 )  i
  integer   ( kind = 4 )  ierror
  integer   ( kind = 4 )  ios
  integer   ( kind = 4 )  lchar
  character ( len = 255 ) master_word
  integer   ( kind = 4 )  master_word_length
  character ( len = 255 ) master_word2
  logical                 s_s_subanagram_sorted
  character ( len = 255 ) string
  logical                 subanagram
  character ( len = 255 ) subanagram_filename
  integer   ( kind = 4 )  subanagram_number
  integer   ( kind = 4 )  subanagram_unit
  character ( len = 255 ) word
  character ( len = 255 ) word_filename
  integer   ( kind = 4 )  word_length
  integer   ( kind = 4 )  word_length_max
  integer   ( kind = 4 )  word_length_min
  integer   ( kind = 4 )  word_long_number
  integer   ( kind = 4 )  word_number
  integer   ( kind = 4 )  word_unit
  character ( len = 255 ) word2

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SUBANAGRAM'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Find subanagrams of a given word.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Given an master word (or phrase), this'
  write ( *, '(a)' ) '  program examines a file full of words, and looks '
  write ( *, '(a)' ) '  for "subanagrams", that is, words formed by using '
  write ( *, '(a)' ) '  some (possibly all) of the letters of the master word.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )

  if ( arg_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the minimum word length:'
    read ( *, * ) word_length_min
  else
    arg = 1
    call getarg ( arg, string )
    call s_to_i4 ( string, word_length_min, ierror, lchar )
  end if

  if ( arg_num < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the master word:'
    read ( *, '(a)' ) master_word
  else
    arg = 2
    call getarg ( arg, master_word )
  end if
!
!  Preprocess the master string.
!
  master_word2 = master_word
  call s_nonalpha_delete ( master_word2 )
  master_word_length = len_trim ( master_word2 )
  call s_low ( master_word2(1:master_word_length) )
  call s_cat ( master_word2, '.txt', subanagram_filename )
  call s_sort_a ( master_word2(1:master_word_length) )
!
!  Report the command line arguments.
!
  word_length_max = master_word_length

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Master word: "' // trim ( master_word ) // '".'
  write ( *, '(a)' ) '  Sorted:      "' // trim ( master_word2 ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Minimum subanagrams length: ', word_length_min
  write ( *, '(a,i8)' ) '  Maximum subanagrams length: ', word_length_max
!
!  Open the files.
!
  word_filename = '/Users/burkardt/public_html/datasets/words/wordlist.txt'

  call get_unit ( word_unit )
  open ( unit = word_unit, file = word_filename, status = 'old', iostat = ios )
  if ( ios .ne. 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUBANAGRAM - Fatal error!'
    write ( *, '(a)' ) '  Could not open the word list file "' &
      // trim ( word_filename ) // '".'
    stop
  end if

  call get_unit ( subanagram_unit )
  open ( unit = subanagram_unit, file = subanagram_filename, status = 'replace' )
!
!  Read a word from the word file, and judge it.
!
  word_number = 0
  word_long_number = 0
  subanagram_number = 0

  do

    read ( word_unit, '(a)', iostat = ios ) word

    if ( ios /= 0 ) then
      exit
    end if

    word_number = word_number + 1

    word_length = len_trim ( word )

    if ( word_length_min <= word_length .and. &
                            word_length <= word_length_max ) then

      word_long_number = word_long_number + 1

      word2 = word

      call s_low ( word2(1:word_length) )

      call s_sort_a ( word2(1:word_length) )

      subanagram = s_s_subanagram_sorted ( trim ( master_word2 ), trim ( word2 ) )

      if ( subanagram ) then
        write ( subanagram_unit, '(a)' ) trim ( word )
        subanagram_number = subanagram_number + 1
      end if

    end if

  end do

  close ( unit = word_unit )
  close ( unit = subanagram_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8,a)' ) word_number,       '  words checked;'
  write ( *, '(2x,i8,a)' ) word_long_number,  &
    '  words checked that were the right length;'
  write ( *, '(2x,i8,a)' ) subanagram_number, &
    '  words checked that were the right length AND subanagrams.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The subanagrams were written to the file "' &
    // trim ( subanagram_filename ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SUBANAGRAM:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ch_low ( ch )

!*****************************************************************************80
!
!! CH_LOW lowercases a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
!    which guarantee the ASCII collating sequence.
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

  character              ch
  integer   ( kind = 4 ) i

  i = iachar ( ch )

  if ( 65 <= i .and. i <= 90 ) then
    ch = achar ( i + 32 )
  end if

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
  logical              lopen

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
subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT concatenates two strings to make a third string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2000
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
    s3 = trim ( s1 ) // trim ( s2 )
  end if

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

  integer   ( kind = 4 ) i
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len_trim ( s )

  do i = 1, s_length
    call ch_low ( s(i:i) )
  end do

  return
end
subroutine s_nonalpha_delete ( s )

!*****************************************************************************80
!
!! S_NONALPHA_DELETE removes nonalphabetic characters from a string.
!
!  Discussion:
!
!    The remaining characters are left justified and blank padded.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2009
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

  character              ch
  integer   ( kind = 4 ) get
  integer   ( kind = 4 ) put
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  put = 0
  s_length = len_trim ( s )

  do get = 1, s_length

    ch = s(get:get)

    if ( ( lle ( 'A', ch ) .and. lle ( ch, 'Z' ) ) .or. &
         ( lle ( 'a', ch ) .and. lle ( ch, 'z' ) ) ) then
      put = put + 1
      s(put:put) = ch
    end if

  end do

  s(put+1:s_length) = ' '

  return
end
function s_s_subanagram_sorted ( s1, s2 )

!*****************************************************************************80
!
!! S_S_SUBANAGRAM_SORTED determines if S2 is a "subanagram" of S1.
!
!  Discussion:
!
!    This routine assumes that S1 and S2 have already been sorted.
!
!    S2 is an anagram of S1 if S2 can be formed by permuting the letters
!    of S1
!
!    S2 is an subanagram of S1 if S2 can be formed by selecting SOME of
!    the letters of S1 and permuting them.
!
!    Blanks (trailing or otherwise), punctuation, and capitalization
!    are all significant, so be sure to input exactly the information
!    you want to check.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the master string.
!
!    Input, character ( len = * ) S2, the second string.
!
!    Output, logical S_S_SUBANAGRAM_SORTED is TRUE if S2 is a subanagram of S1.
!
  implicit none

  integer   ( kind = 4 ) i1
  integer   ( kind = 4 ) i2
  logical                s_s_subanagram_sorted
  character ( len = * )  s1
  integer   ( kind = 4 ) s1_length
  character ( len = * )  s2
  integer   ( kind = 4 ) s2_length

  s_s_subanagram_sorted = .false.

  s1_length = len ( s1 )
  s2_length = len ( s2 )

  i1 = 0

  do i2 = 1, s2_length

    do

      i1 = i1 + 1
!
!  Ran out of S1 before finishing.  No match is possible.
!
      if ( s1_length < i1 ) then
        return
      end if
!
!  The current character in S1 is already greater than the character in S2.
!  No match is possible.
!
      if ( llt ( s2(i2:i2), s1(i1:i1) ) ) then
        return
      end if
!
!  Found an exact match for current character.  Keep going.
!
      if ( s1(i1:i1) == s2(i2:i2) ) then
        exit
      end if
!
!  Didn't find a match, but one might be possible if we increase I1.
!
    end do

  end do
!
!  We matched every character of S2 with something in S1.
!
  s_s_subanagram_sorted = .true.

  return
end
subroutine s_sort_a ( s )

!*****************************************************************************80
!
!! S_SORT_A sorts a string into ascending order.
!
!  Discussion:
!
!    The string is assumed to be short, and so a simple bubble sort is used.
!
!    ALL the characters are sorted, including blanks and punctuation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be sorted.
!
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) k
  character ( len = * )  s
  integer   ( kind = 4 ) s_length

  s_length = len ( s )

  do i = 1, s_length - 1

    c = s(i:i)
    j = i

    do k = i + 1, s_length
      if ( iachar ( s(k:k) ) < iachar ( s(j:j) ) ) then
        j = k
      end if
    end do

    if ( i /= j ) then
      s(i:i) = s(j:j)
      s(j:j) = c
    end if

  end do

  return
end
subroutine s_to_i4 ( s, value, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2009
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
!    Output, integer ( kind = 4 ) LENGTH, the number of characters
!    of S used to make the integer.
!
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) length
  character ( len = * )  s
  integer   ( kind = 4 ) state
  character              :: TAB = achar ( 9 )
  integer   ( kind = 4 ) value

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

      if ( c == ' ' .or. c == TAB ) then

      else if ( c == '-' ) then
        state = 1
        isgn = -1
      else if ( c == '+' ) then
        state = 1
        isgn = +1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' .or. c == TAB ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 2, have read at least one digit, expecting more.
!
    else if ( state == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + iachar ( c ) - iachar ( '0' )

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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
