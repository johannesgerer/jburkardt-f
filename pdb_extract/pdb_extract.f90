program main

!*****************************************************************************80
!
!! MAIN is the main program for PDB_EXTRACT.
!
!  Description:
!
!    PDB_EXTRACT extracts data from a PDB file.
!
!    The current version of this program searches a PDB file for lines
!    whose fourth word is 'CA'.
!
!    From those lines, it reads the 9th and 10th words.
!
!    It interprets the 10th word as a real number, multiplying it
!      by 3/(8*pi**2).
!
!    It then writes the 9th word and the scaled 10th word to a new file.
!
!  Usage:
!
!    pdb_extract input_name output_name
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2000
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilen
  integer ( kind = 4 ) input_line
  character ( len = 80 ) input_name
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ipxfargc
  integer ( kind = 4 ) lchar
  character ( len = 256 ) line
  integer ( kind = 4 ) numarg
  integer ( kind = 4 ) nword
  integer ( kind = 4 ) output_line
  character ( len = 80 ) output_name
  integer ( kind = 4 ) output_unit
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  character ( len = 80 ) word1
  character ( len = 80 ) word2
  character ( len = 80 ) word3

  call timestamp ( )

  ierror = 0
!
!  Old style:
!
  numarg = iargc ( )
!
!  New style:
!
! numarg = ipxfargc ( )
!
!  Get the input filename.
!
  if ( numarg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_EXTRACT'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Input required!'
    write ( *, '(a)' ) '  Enter the input filename.'
    read ( *, '(a)' ) input_name

  else

    iarg = 1
!
!  Old style:
!
    call getarg ( iarg, input_name )
!
!  New style:
!
!   call pxfgetarg ( iarg, input_name, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'PDB_EXTRACT - Error!'
!     write ( *, '(a)' ) '  Could not read the argument.'
!     stop
!   end if

  end if
!
!  Get the output filename.
!
  if ( numarg < 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_EXTRACT - Input required!'
    write ( *, '(a)' ) '  Enter the output filename.'
    read ( *, '(a)' ) output_name

  else

    iarg = 2
!
!  Old style:
!
    call getarg ( iarg, output_name )
!
!  New style:
!
!   call pxfgetarg ( iarg, output_name, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'PDB_EXTRACT - Error!'
!     write ( *, '(a)' ) '  Could not read the argument.'
!     stop
!   end if

  end if
!
!  Open the input file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_EXTRACT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file:' // trim ( input_name )
    stop
  end if

  input_line = 0
!
!  Open the output file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_EXTRACT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file: ' &
      // trim ( output_name )
    close ( unit = input_unit )
    stop
  end if

  output_line = 0
!
!  Read a line of input.
!
  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    input_line = input_line + 1
!
!  Get the word in column 4.
!
    call word_find ( line, 4, word1, nword )
!
!  If that word is 'CA', then
!  * get columns 9 and 10 of that line,
!  * interpret column 10 as a real number and multiply it by 8 / pi**2.
!  * write the two items to the output file.
!
    if ( trim ( word1 ) == 'CA' ) then

      call word_find ( line, 9, word2, nword )

      call word_find ( line, 10, word3, nword )

      call s_to_r8 ( word3, r, ierror, lchar )

      r = 3.0D+00 * r / ( 8.0D+00 * pi**2 )

      write ( output_unit, '(a,2x,g14.6)' ) trim ( word2 ), r

      output_line = output_line + 1

    end if

  end do

  close ( unit = input_unit )

  close ( unit = output_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PDB_EXTRACT - Note:'
  write ( *, '(a,i6)' ) '  Input file lines:  ', input_line
  write ( *, '(a,i6)' ) '  Output file lines: ', output_line
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PDB_EXTRACT'
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
function ch_eqi ( c1, c2 )

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
!    14 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c2
  character cc1
  character cc2

  cc1 = c1
  cc2 = c2

  call ch_cap ( cc1 )
  call ch_cap ( cc2 )

  if ( cc1 == cc2 ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
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
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

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
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT.
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

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
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!    Output, integer LCHAR, the number of characters read from
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
      else if ( ihave > 1 ) then
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
      else if ( ihave >= 6 .and. ihave <= 8 ) then
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
    if ( iterm == 1 .or. lchar+1 >= nchar ) then
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
!
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
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
  character ( len = * ) word

  ilo = 0
  ihi = 0
  lchar = len_trim ( s )

  if ( lchar <= 0 ) then
    return
  end if

  if ( iword > 0 ) then

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

      if ( i > lchar ) then

        if ( jword == iword ) then
          ilo = jlo
          ihi = lchar
          nchar = lchar + 1 - jlo
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
    jlo = lchar
    jhi = lchar
    i = lchar

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
