program main

!*****************************************************************************80
!
!! AMSTERDAM works on the "AMSTERDAM" problem.
!
!  Discussion:
!
!    TEXTROPOLIS is an iPhone game that seeks as many subanagrams of
!    certain city names as possible, although they are required to
!    have at least 4 letters.
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
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  logical amsterdam
  character ( len = 80 ) :: amsterdam_file = 'amsterdam.txt'
  integer ( kind = 4 ) amsterdam_number
  character ( len = 13 ) :: amsterdam_string = 'amsterdam'
  integer ( kind = 4 ), parameter :: amsterdam_unit = 3
  logical s_s_subanagram_custom
  character ( len = 80 ) string
  integer ( kind = 4 ) string_length
  integer ( kind = 4 ), parameter :: string_length_min = 4
  character ( len = 80 ) string2
  character ( len = 80 ) :: word_file = 'wordlist.txt'
  integer ( kind = 4 ) word_long_number
  integer ( kind = 4 ) word_number
  integer ( kind = 4 ), parameter :: word_unit = 1

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AMSTERDAM'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Work on the "AMSTERDAM" problem.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  As in the iPhone application Textropolis,'
  write ( *, '(a)' ) '  find as many words as possible using the letters '
  write ( *, '(a)' ) '  in "AMSTERDAM".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine a file "' // trim ( word_file ) // '"'
  write ( *, '(a)' ) '  full of words, and look for "sub-anagrams" '
  write ( *, '(a)' ) '  (words whose letters are a subset of the target ones).'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  We only consider words of length at least ', string_length_min
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Copy sub-anagrams to "' // trim ( amsterdam_file ) // '".'

  word_long_number = 0
  word_number = 0
  amsterdam_number = 0
!
!  Ascending sort the master string.
!
  call s_sort_a ( amsterdam_string )

  open ( unit = word_unit, file = word_file, status = 'old' )
  open ( unit = amsterdam_unit, file = amsterdam_file, status = 'replace' )
!
!  Read a word from the word file, and judge it.
!
  do

    read ( word_unit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then 
      exit
    end if

    word_number = word_number + 1

    string_length = len_trim ( string )

    if ( string_length_min <= string_length ) then

      word_long_number = word_long_number + 1

      string2 = string

      call s_low ( string2 )

      call s_sort_a ( string2(1:string_length) )

      amsterdam = s_s_subanagram_custom ( amsterdam_string, trim ( string2 ) )

      if ( amsterdam ) then
        write ( amsterdam_unit, '(a)' ) trim ( string )
        amsterdam_number = amsterdam_number + 1
      end if

    end if

  end do

  close ( unit = word_unit )
  close ( unit = amsterdam_unit )

  write ( *, * ) ' '
  write ( *, * ) '  Number of words checked                  ', word_number
  write ( *, * ) '  Number of words that were "long enough"  ', word_long_number
  write ( *, * ) '  Number of sub-anagrams found             ', amsterdam_number

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AMSTERDAM:'
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

  character ch
  integer ( kind = 4 ) i

  i = iachar ( ch )

  if ( 65 <= i .and. i <= 90 ) then
    ch = achar ( i + 32 )
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

  integer ( kind = 4 ) i
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len_trim ( s )

  do i = 1, s_length
    call ch_low ( s(i:i) )
  end do

  return
end
function s_s_subanagram_custom ( s1, s2 )

!*****************************************************************************80
!
!! S_S_SUBANAGRAM_CUSTOM determines if S2 is a "subanagram" of S1.
!
!  Discussion:
!
!    This routine has been "customized" because we can assume that
!    S1 and S2 have already been sorted.
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
!    Output, logical S_S_SUBANAGRAM_CUSTOM is TRUE if S2 is a subanagram of S1.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  logical s_s_subanagram_custom
  character ( len = * ) s1
  integer ( kind = 4 ) s1_length
  character ( len = * ) s2
  integer ( kind = 4 ) s2_length

  s_s_subanagram_custom = .false.
!
!  Sort both.
!
! call s_sort_a ( s1 )
! call s_sort_a ( s2 )

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
  s_s_subanagram_custom = .true.

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

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = * ) s
  integer ( kind = 4 ) s_length

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
