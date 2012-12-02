subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
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
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  character ch
  integer ( kind = 4 ) itemp

  itemp = iachar ( ch )
 
  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
  end if
 
  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.  
!
!  Discussion:
!
!    CH_EQI ( 'A', 'a' ) is TRUE.
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
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
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
!    If CH was 'illegal', then DIGIT is -1.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then
 
    digit = iachar ( ch ) - 48
 
  else if ( ch == ' ' ) then
 
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

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * )  s
  integer ( kind = 4 ) state
  character :: TAB = achar ( 9 )
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
subroutine s_to_i4vec ( s, n, i4vec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an integer vector from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2003
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
!    Output, integer ( kind = 4 ) I4VEC(N), the values read from the string.
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
  integer ( kind = 4 ) i4vec(n)
  integer ( kind = 4 ) length
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_i4 ( s(ilo:), i4vec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

  end do

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 value from a string.
!
!  Discussion:
!
!    An "R8" value is simply a real number to be stored as a
!    variable of type "real ( kind = 8 )".
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
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
!    12 January 2009
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
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  integer ( kind = 4 ) s_length
  character :: TAB = achar ( 9 )

  s_length = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( s_length < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
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
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
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
!  Scientific notation exponent marker.
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
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

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
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
  if ( iterm /= 1 .and. length+1 == s_length ) then
    length = s_length
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, r8vec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Discussion:
!
!    An R8VEC is a vector of real values, of type "real ( kind = 8 )".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2005
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
!    Output, real ( kind = 8 ) R8VEC(N), the values read from the string.
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
  real ( kind = 8 ) r8vec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), r8vec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

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
!    Output, integer ( kind = 4 ) WORD_NUM, the number of "words" in the 
!    string.  Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  character ( len = * ) s
  integer ( kind = 4 ) s_length
  integer ( kind = 4 ) word_num

  word_num = 0
  s_length = len ( s )

  if ( s_length <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, s_length

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
subroutine xy_data_print ( point_num, xy )

!*****************************************************************************80
!
!! XY_DATA_PRINT prints the data of an XY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) XY(2,POINT_NUM), the point coordinates.
!
  implicit none

  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) j
  real ( kind = 8 ) xy(2,point_num)

  do j = 1, point_num
    write ( *, '(2x,f14.6,2x,f14.6)' ) xy(1:2,j)
  end do

  return
end
subroutine xy_data_read ( input_filename, point_num, xy )

!*****************************************************************************80
!
!! XY_DATA_READ reads data from an XY file.
!
!  Discussion:
!
!    The number of points in the file can be determined by calling
!    XY_HEADER_READ first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.  The program
!    will stop reading data once POINT_NUM values have been read.
!
!    Output, real ( kind = 8 ) XY(2,POINT_NUM), the point coordinates.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 )  point_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  real ( kind = 8 ) temp(dim_num)
  character ( len = 255 ) text
  real ( kind = 8 ) xy(2,point_num)

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XY_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
    stop
  end if

  i = 0

  do while ( i < point_num )

    read ( input_unit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XY_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected I/O error.'
      stop
    end if

    if ( text(1:1) == '#' .or. len_trim ( text ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( text, dim_num, temp, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XY_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  S_TO_R8VEC returned an error flag.'
      stop
    end if

    i = i + 1

    xy(1,i) = temp(1)
    xy(2,i) = temp(2)

  end do

  close ( unit = input_unit )

  return
end
subroutine xy_data_write ( output_unit, point_num, xy )

!*****************************************************************************80
!
!! XY_DATA_WRITE writes the data of an XY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) XY(2,POINT_NUM), the point coordinates.
!
  implicit none

  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) j
  integer ( kind = 4 ) output_unit
  real ( kind = 8 ) xy(2,point_num)

  do j = 1, point_num
    write ( output_unit, '(2x,f14.6,2x,f14.6)' ) xy(1:2,j)
  end do

  return
end
subroutine xy_example ( point_num, xy )

!*****************************************************************************80
!
!! XY_EXAMPLE sets up sample XY data suitable for an XY file.
!
!  Discussion:
!
!    The points lie on a spiral in the unit square.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Output, real ( kind = 8 ) XY(2,POINT_NUM), the point coordinates.
!
  implicit none

  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  integer ( kind = 4 ), parameter :: turns = 5
  real ( kind = 8 ) xy(2,point_num)

  do j = 1, point_num
    r = real ( j - 1, kind = 8 ) / real ( point_num - 1, kind = 8 )
    theta = real ( turns, kind = 8 ) * r * ( 2.0D+00 * r8_pi )
    xy(1,j) = r * cos ( theta )
    xy(2,j) = r * sin ( theta )
  end do

  xy(1,1:point_num) = 0.5D+00 * ( 1.0D+00 + xy(1,1:point_num) )
  xy(2,1:point_num) = 0.5D+00 * ( 1.0D+00 + xy(2,1:point_num) )

  return
end
subroutine xy_header_print ( point_num )

!*****************************************************************************80
!
!! XY_HEADER_PRINT prints the header of an XY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
  implicit none

  integer ( kind = 4 ) point_num

  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  return
end
subroutine xy_header_read ( input_filename, point_num )

!*****************************************************************************80
!
!! XY_HEADER_READ determines the number of pairs of data in an XY file.
!
!  Discussion:
!
!    This routine assumes that the file contains exactly three kinds of
!    records:
!
!    COMMENTS which begin with a '#' character in column 1;
!    BLANKS which contain nothing but 'whitespace';
!    XY coordinates, which each contain one pair of real values.
!
!    The routine ignores comments and blanks and returns
!    the number of records containing XY coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the number of points in the file.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 )  point_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  real ( kind = 8 ) temp(dim_num)
  character ( len = 255 ) text

  point_num = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XY_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '".'
    stop
  end if

  do

    read ( input_unit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    if ( text(1:1) == '#' .or. len_trim ( text ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( text, dim_num, temp, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XY_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  S_TO_R8VEC returned an error flag.'
      stop
    end if

    point_num = point_num + 1

  end do

  close ( unit = input_unit )

  return
end
subroutine xy_header_write ( output_filename, output_unit, point_num )

!*****************************************************************************80
!
!! XY_HEADER_WRITE writes the header of an XY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
  implicit none

  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) point_num
!
!  Write the header.
!
  write ( output_unit, '(a)'    ) '#  ' // trim ( output_filename )
  write ( output_unit, '(a)'    ) '#  created by xy_io::xy_header_write.f90'
  write ( output_unit, '(a)'    ) '#'
  write ( output_unit, '(a,i8)' ) '#  Number of points = ', point_num
  write ( output_unit, '(a)'    ) '#'

  return
end
subroutine xy_read_test ( input_filename )

!*****************************************************************************80
!
!! XY_READ_TEST tests the XY file read routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file
!    containing the ASCII portable pixel map data.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: xy(:,:)
!
!  Read the header.
!
  call xy_header_read ( input_filename, point_num )
!
!  Allocate the data.
!
  allocate ( xy(2,point_num) )
!
!  Read the data.
!
  call xy_data_read ( input_filename, point_num, xy )

  deallocate ( xy )

  return
end
subroutine xy_write ( output_filename, point_num, xy )

!*****************************************************************************80
!
!! XY_WRITE writes an XY file.
!
!  Example:
!
!    # my_file.xy
!    # created by XY_IO::XY_WRITE.
!    #
!    #  Number of points = 5
!    #
!    0.0  0.0
!    1.0  2.0
!    3.0  5.0
!    2.0  1.0
!    8.0  7.5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the file
!    to which the data should be written.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) XY(2,POINT_NUM), the point coordinates.
!
  implicit none

  integer ( kind = 4 ) point_num

  logical, parameter :: debug = .false.
  integer ( kind = 4 ) ios
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
  real ( kind = 8 ) xy(2,point_num)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XY_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '".'
    stop
  end if
!
!  Write the header.
!
  call xy_header_write ( output_filename, output_unit, point_num )
!
!  Write the data.
!
  call xy_data_write ( output_unit, point_num, xy )
!
!  Close the file.
!
  close ( unit = output_unit )
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XY_WRITE - Note:'
    write ( *, '(a)' ) '  The data was written.'
    write ( *, '(a,i8)' ) '  Number of points =    ', point_num
  end if

  return
end
subroutine xy_write_test ( output_filename )

!*****************************************************************************80
!
!! XY_WRITE_TEST tests the XY write routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the file
!    to contain the data.
!
  implicit none

  character ( len = * ) output_filename
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable :: xy(:,:)

  point_num = 100
!
!  Allocate memory.
!
  allocate ( xy(2,point_num) )
!
!  Set the data.
!
  call xy_example ( point_num, xy )
!
!  Write the data to the file.
!
  call xy_write ( output_filename, point_num, xy )

  deallocate ( xy )

  return
end
subroutine xyf_data_print ( point_num, face_num, face_data_num, &
  face_pointer, face_data )

!*****************************************************************************80
!
!! XYF_DATA_PRINT prints the data of an XYF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_DATA_NUM, the number of face items.
!
!    Input, integer ( kind = 4 ) FACE_POINTER(LINE_NUM+1), pointers to the
!    first line item for each face.
!
!    Input, integer ( kind = 4 ) FACE_DATA(LINE_DATA_NUM), indices
!    of points that form faces.
!
  implicit none

  integer ( kind = 4 ) face_data_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_data(face_data_num)
  integer ( kind = 4 ) face_pointer(face_num+1)

  do face = 1, face_num
    write ( *, '(2x,i4,2x,i8,2x,i8)' ) &
      face, face_pointer(face), face_pointer(face+1)-1
  end do

  write ( *, '(a)' ) ' '

  do face = 1, face_num
    do i = face_pointer(face), face_pointer(face+1) - 1
      write ( *, '(2x,i8)', advance = 'NO' ) face_data(i)
    end do
    write ( *, '(1x)', advance = 'YES' )
  end do

  return
end
subroutine xyf_data_read ( input_filename, face_num, face_data_num, &
  face_pointer, face_data )

!*****************************************************************************80
!
!! XYF_DATA_READ reads the data in an XYF file.
!
!  Discussion:
!
!    This routine assumes that the file contains exactly three kinds of
!    records:
!
!    COMMENTS which begin with a '#' character in column 1;
!    BLANKS which contain nothing but 'whitespace';
!    FACE ITEMS, which are indices of points on a face.
!
!    The routine ignores comments and blanks and returns
!    the number of face items.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_DATA_NUM, the number of face items.
!
!    Output, integer ( kind = 4 ) FACE_POINTER(FACE_NUM+1), pointers to the
!    first face item for each face.
!
!    Output, integer ( kind = 4 ) FACE_DATA(FACE_DATA_NUM), the face items.
!
  implicit none

  integer ( kind = 4 ) face_data_num
  integer ( kind = 4 ) face_num

  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_data(face_data_num)
  integer ( kind = 4 ) face_pointer(face_num+1)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) n
  character ( len = 255 ) text

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYF_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '".'
    stop
  end if

  face = 0
  face_pointer(1) = 1

  do while ( face < face_num )

    read ( input_unit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XYF_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of information.'
      stop
    end if

    if ( text(1:1) == '#' .or. len_trim ( text ) == 0 ) then
      cycle
    end if

    face = face + 1

    call s_word_count ( text, n )
    face_pointer(face+1) = face_pointer(face) + n
 
    ilo = face_pointer(face)
    ihi = face_pointer(face+1) - 1

    call s_to_i4vec ( text, n, face_data(ilo:ihi), ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XYF_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected error from S_TO_I4VEC.'
      stop
    end if

  end do

  close ( unit = input_unit )

  return
end
subroutine xyf_data_write ( output_unit, point_num, face_num, face_data_num, &
  face_pointer, face_data )

!*****************************************************************************80
!
!! XYF_DATA_WRITE writes the data of an XYF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_DATA_NUM, the number of face items.
!
!    Input, integer ( kind = 4 ) FACE_POINTER(LINE_NUM+1), pointers to the
!    first item for each face
!
!    Input, integer ( kind = 4 ) FACE_DATA(LINE_DATA_NUM), indices
!    of points that form faces.
!
  implicit none

  integer ( kind = 4 ) face_data_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_data(face_data_num)
  integer ( kind = 4 ) face_pointer(face_num+1)
  integer ( kind = 4 ) output_unit

  do face = 1, face_num
    do i = face_pointer(face), face_pointer(face+1) - 1
      write ( output_unit, '(2x,i8)', advance = 'NO' ) face_data(i)
    end do
    write ( output_unit, '(1x)', advance = 'YES' )
  end do

  return
end
subroutine xyf_example ( point_num, face_num, face_data_num, xy, &
  face_pointer, face_data )

!*****************************************************************************80
!
!! XYF_EXAMPLE sets data suitable for a pair of XY and XYF files.
!
!  Discussion:
!
!    There are 65 points.
!    There are 48 faces.
!    There are 48*4=192 face items.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) FACE_DATA_NUM, the number of line items.
!
!    Output, real ( kind = 8 ) XY(2,POINT_NUM), the point coordinates.
!
!    Output, integer ( kind = 4 ) FACE_POINTER(FACE_NUM+1), pointers to the
!    first face item for each face.
!
!    Output, integer ( kind = 4 ) FACE_DATA(FACE_DATA_NUM), indices
!    of points that form faces.
!
  implicit none

  integer ( kind = 4 ) face_data_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), parameter :: n_t = 13
  integer ( kind = 4 ), parameter :: n_r = 5
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_data(face_data_num)
  integer ( kind = 4 ) face_pointer(face_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: r_min = 1.0D+00
  real ( kind = 8 ), parameter :: r_max = 3.0D+00
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: t_min = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: t_max = 0.0D+00
  real ( kind = 8 ) xy(2,point_num)

  k = 0
  do j = 1, n_r

    r = ( real ( n_r - j ,    kind = 8 ) * r_min &
        + real (       j - 1, kind = 8 ) * r_max ) &
        / real ( n_r     - 1, kind = 8 )

    do i = 1, n_t

      t = ( real ( n_t - i,     kind = 8 ) * t_min &
          + real (       i - 1, kind = 8 ) * t_max ) &
          / real ( n_t     - 1, kind = 8 )

      k = k + 1
      xy(1,k) = r * cos ( t )
      xy(2,k) = r * sin ( t )

    end do
  end do

  face = 0
  k = 0

  do j = 1, n_r - 1
    do i = 1, n_t - 1

      face = face + 1
      face_pointer(face) = k + 1

      k = k + 1
      face_data(k) = i     + ( j - 1 ) * n_t
      k = k + 1
      face_data(k) = i + 1 + ( j - 1 ) * n_t
      k = k + 1
      face_data(k) = i + 1 + ( j     ) * n_t
      k = k + 1
      face_data(k) = i     + ( j     ) * n_t

    end do
  end do

  face = face + 1
  face_pointer(face) = k + 1

  return
end
subroutine xyf_example_size ( point_num, face_num, face_data_num )

!*****************************************************************************80
!
!! XYF_EXAMPLE_SIZE sizes the data to be created by XYF_EXAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of face.
!
!    Output, integer ( kind = 4 ) FACE_DATA_NUM, the number of face items.
!
  implicit none

  integer ( kind = 4 ) face_data_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), parameter :: n_t = 13
  integer ( kind = 4 ), parameter :: n_r = 5
  integer ( kind = 4 ) point_num

  face_data_num = 4 * ( n_t - 1 ) * ( n_r - 1 )
  face_num = ( n_t - 1 ) * ( n_r - 1 )
  point_num = n_t * n_r

  return
end
subroutine xyf_header_print ( point_num, face_num, face_data_num )

!*****************************************************************************80
!
!! XYF_HEADER_PRINT prints the header of an XYF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of face.
!
!    Input, integer ( kind = 4 ) FACE_DATA_NUM, the number of face items.
!
  implicit none

  integer ( kind = 4 ) face_data_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) point_num

  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of points =     ', point_num
  write ( *, '(a,i8)' ) '  Number of faces =      ', face_num
  write ( *, '(a,i8)' ) '  Number of face items = ', face_data_num

  return
end
subroutine xyf_header_read ( input_filename, face_num, face_data_num )

!*****************************************************************************80
!
!! XYF_HEADER_READ determines the number of face items in an XZF file.
!
!  Discussion:
!
!    This routine assumes that the file contains exactly three kinds of
!    records:
!
!    COMMENTS which begin with a '#' character in column 1;
!    BLANKS which contain nothing but 'whitespace';
!    FACE ITEMS, which are indices of points on a face;
!
!    The routine ignores comments and blanks and returns
!    the number of face items.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) FACE_DATA_NUM, the number of face items.
!
  implicit none

  integer ( kind = 4 ) face_data_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_val
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) length
  integer ( kind = 4 ) n
  character ( len = 255 ) text

  face_data_num = 0
  face_num = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYF_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '".'
    stop
  end if

  do

    read ( input_unit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    if ( text(1:1) == '#' .or. len_trim ( text ) == 0 ) then
      cycle
    end if

    call s_word_count ( text, n )

    face_data_num = face_data_num + n

    face_num = face_num + 1

  end do

  close ( unit = input_unit )

  return
end
subroutine xyf_header_write ( output_filename, output_unit, point_num, &
  face_num, face_data_num )

!*****************************************************************************80
!
!! XYF_HEADER_WRITE writes the header of an XYF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_DATA_NUM, the number of face items.
!
  implicit none

  integer ( kind = 4 ) face_data_num
  integer ( kind = 4 ) face_num
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) point_num
!
!  Write the header.
!
  write ( output_unit, '(a)'    ) '#  ' // trim ( output_filename )
  write ( output_unit, '(a)'    ) '#  created by xy_io::xyf_header_write.f90'
  write ( output_unit, '(a)'    ) '#'
  write ( output_unit, '(a,i8)' ) '#  Number of points =     ', point_num
  write ( output_unit, '(a,i8)' ) '#  Number of faces =      ', face_num
  write ( output_unit, '(a,i8)' ) '#  Number of face items = ', face_data_num
  write ( output_unit, '(a)'    ) '#'

  return
end
subroutine xyf_write ( output_filename, point_num, face_num, face_data_num, &
  face_pointer, face_data )

!*****************************************************************************80
!
!! XYF_WRITE writes an XYF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the file
!    to which the data should be written.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of face.
!
!    Input, integer ( kind = 4 ) FACE_DATA_NUM, the number of face items.
!
!    Input, integer ( kind = 4 ) FACE_POINTER(FACE_NUM+1), pointers to the
!    first face item for each line.
!
!    Input, integer ( kind = 4 ) FACE_DATA(FACE_DATA_NUM), indices
!    of points that form face.
!
  implicit none

  integer ( kind = 4 ) face_data_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) ios
  integer ( kind = 4 ) face_data(face_data_num)
  integer ( kind = 4 ) face_pointer(face_num+1)
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYF_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '".'
    stop
  end if
!
!  Write the header.
!
  call xyf_header_write ( output_filename, output_unit, point_num, face_num, &
    face_data_num )
!
!  Write the data.
!
  call xyf_data_write ( output_unit, point_num, face_num, face_data_num, &
    face_pointer, face_data )
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine xyl_data_print ( point_num, line_num, line_data_num, &
  line_pointer, line_data )

!*****************************************************************************80
!
!! XYL_DATA_PRINT prints the data of an XYL file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
!    Input, integer ( kind = 4 ) LINE_POINTER(LINE_NUM+1), pointers to the
!    first line item for each line.
!
!    Input, integer ( kind = 4 ) LINE_DATA(LINE_DATA_NUM), indices
!    of points that form lines.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) line
  integer ( kind = 4 ) line_data(line_data_num)
  integer ( kind = 4 ) line_pointer(line_num+1)

  do line = 1, line_num
    write ( *, '(2x,i4,2x,i8,2x,i8)' ) &
      line, line_pointer(line), line_pointer(line+1)-1
  end do

  write ( *, '(a)' ) ' '

  do line = 1, line_num
    do i = line_pointer(line), line_pointer(line+1) - 1
      write ( *, '(2x,i8)', advance = 'NO' ) line_data(i)
    end do
    write ( *, '(1x)', advance = 'YES' )
  end do

  return
end
subroutine xyl_data_read ( input_filename, line_num, line_data_num, &
  line_pointer, line_data )

!*****************************************************************************80
!
!! XYL_DATA_READ reads the data in an XYL file.
!
!  Discussion:
!
!    This routine assumes that the file contains exactly three kinds of
!    records:
!
!    COMMENTS which begin with a '#' character in column 1;
!    BLANKS which contain nothing but 'whitespace';
!    LINE ITEMS, which are indices of points on a line, or -1 to terminate a line.
!
!    The routine ignores comments and blanks and returns
!    the number of line items.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
!    Output, integer ( kind = 4 ) LINE_POINTER(LINE_NUM+1), pointers to the
!    first line item for each line.
!
!    Output, integer ( kind = 4 ) LINE_DATA(LINE_DATA_NUM), the line items.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) line
  integer ( kind = 4 ) line_data(line_data_num)
  integer ( kind = 4 ) line_pointer(line_num+1)
  integer ( kind = 4 ) n
  character ( len = 255 ) text

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYL_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '".'
    stop
  end if

  line = 0
  line_pointer(1) = 1

  do while ( line < line_num )

    read ( input_unit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XYL_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of information.'
      stop
    end if

    if ( text(1:1) == '#' .or. len_trim ( text ) == 0 ) then
      cycle
    end if

    line = line + 1

    call s_word_count ( text, n )
    line_pointer(line+1) = line_pointer(line) + n
 
    ilo = line_pointer(line)
    ihi = line_pointer(line+1) - 1

    call s_to_i4vec ( text, n, line_data(ilo:ihi), ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XYL_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected error from S_TO_I4VEC.'
      stop
    end if

  end do

  close ( unit = input_unit )

  return
end
subroutine xyl_data_write ( output_unit, point_num, line_num, line_data_num, &
  line_pointer, line_data )

!*****************************************************************************80
!
!! XYL_DATA_WRITE writes the data of an XYL file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
!    Input, integer ( kind = 4 ) LINE_POINTER(LINE_NUM+1), pointers to the
!    first line item for each line.
!
!    Input, integer ( kind = 4 ) LINE_DATA(LINE_DATA_NUM), indices
!    of points that form lines.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) line
  integer ( kind = 4 ) line_data(line_data_num)
  integer ( kind = 4 ) line_pointer(line_num+1)
  integer ( kind = 4 ) output_unit

  do line = 1, line_num
    do i = line_pointer(line), line_pointer(line+1) - 1
      write ( output_unit, '(2x,i8)', advance = 'NO' ) line_data(i)
    end do
    write ( output_unit, '(1x)', advance = 'YES' )
  end do

  return
end
subroutine xyl_example ( point_num, line_num, line_data_num, xy, &
  line_pointer, line_data )

!*****************************************************************************80
!
!! XYL_EXAMPLE sets data suitable for a pair of XY and XYL files.
!
!  Discussion:
!
!    There are 13 points.
!    There are 3 lines.
!    There are 15 line data items.
!
!         4 12-11
!         /\ | |
!        /  \| |
!       /   13 |
!      /      \10
!     /        \
!    5          3
!    |          |
!    |     9--8 |
!    |     |  | |
!    |     |  | |
!    |     6--7 |
!    |          |
!    1----------2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
!    Output, real ( kind = 8 ) XY(2,POINT_NUM), the point coordinates.
!
!    Output, integer ( kind = 4 ) LINE_POINTER(LINE_NUM+1), pointers to the
!    first line item for each line.
!
!    Output, integer ( kind = 4 ) LINE_DATA(LINE_DATA_NUM), indices
!    of points that form lines.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) line_data(line_data_num)
  integer ( kind = 4 ) line_pointer(line_num+1)
  real ( kind = 8 ) xy(2,point_num)

  xy(1:2,1:point_num) = reshape ( (/ &
     0.0D+00,   0.0D+00, & 
     6.0D+00,   0.0D+00, &
     6.0D+00,   7.0D+00, &
     3.0D+00,  10.0D+00, &
     0.0D+00,   7.0D+00, &
     4.0D+00,   1.0D+00, &
     5.0D+00,   1.0D+00, &
     5.0D+00,   4.0D+00, &
     4.0D+00,   4.0D+00, &
     5.0D+00,   8.0D+00, &
     5.0D+00,  11.0D+00, &
     4.0D+00,  11.0D+00, &
     4.0D+00,   9.0D+00 /), (/ 2, point_num /) )

  line_pointer(1:line_num+1) = (/ 1, 7, 12, 16 /)

  line_data(1:line_data_num) = (/ &
     1,  2,  3,  4,  5,  1, &
     6,  7,  8,  9,  6, &
    10, 11, 12, 13 /)

  return
end
subroutine xyl_example_size ( point_num, line_num, line_data_num )

!*****************************************************************************80
!
!! XYL_EXAMPLE_SIZE sizes the data to be created by XYL_EXAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Output, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num

  line_data_num = 15
  line_num = 3
  point_num = 13

  return
end
subroutine xyl_header_print ( point_num, line_num, line_data_num )

!*****************************************************************************80
!
!! XYL_HEADER_PRINT prints the header of an XYL file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num

  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of points =     ', point_num
  write ( *, '(a,i8)' ) '  Number of lines =      ', line_num
  write ( *, '(a,i8)' ) '  Number of line items = ', line_data_num

  return
end
subroutine xyl_header_read ( input_filename, line_num, line_data_num )

!*****************************************************************************80
!
!! XYL_HEADER_READ determines the number of line items in an XYL file.
!
!  Discussion:
!
!    This routine assumes that the file contains exactly three kinds of
!    records:
!
!    COMMENTS which begin with a '#' character in column 1;
!    BLANKS which contain nothing but 'whitespace';
!    LINE ITEMS, which are indices of points on a line, or -1 to terminate a line.
!
!    The routine ignores comments and blanks and returns
!    the number of line items.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Output, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_val
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) length
  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) n
  character ( len = 255 ) text

  line_data_num = 0
  line_num = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYL_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '".'
    stop
  end if

  do

    read ( input_unit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    if ( text(1:1) == '#' .or. len_trim ( text ) == 0 ) then
      cycle
    end if

    call s_word_count ( text, n )

    line_data_num = line_data_num + n

    line_num = line_num + 1

  end do

  close ( unit = input_unit )

  return
end
subroutine xyl_header_write ( output_filename, output_unit, point_num, &
  line_num, line_data_num )

!*****************************************************************************80
!
!! XYL_HEADER_WRITE writes the header of an XYL file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) point_num
!
!  Write the header.
!
  write ( output_unit, '(a)'    ) '#  ' // trim ( output_filename )
  write ( output_unit, '(a)'    ) '#  created by xy_io::xyl_header_write.f90'
  write ( output_unit, '(a)'    ) '#'
  write ( output_unit, '(a,i8)' ) '#  Number of points =     ', point_num
  write ( output_unit, '(a,i8)' ) '#  Number of lines =      ', line_num
  write ( output_unit, '(a,i8)' ) '#  Number of line items = ', line_data_num
  write ( output_unit, '(a)'    ) '#'

  return
end
subroutine xyl_write ( output_filename, point_num, line_num, line_data_num, &
  line_pointer, line_data )

!*****************************************************************************80
!
!! XYL_WRITE writes an XYL file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the file
!    to which the data should be written.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
!    Input, integer ( kind = 4 ) LINE_DATA_NUM, the number of line items.
!
!    Input, integer ( kind = 4 ) LINE_POINTER(LINE_NUM+1), pointers to the
!    first line item for each line.
!
!    Input, integer ( kind = 4 ) LINE_DATA(LINE_DATA_NUM), indices
!    of points that form lines.
!
  implicit none

  integer ( kind = 4 ) line_data_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) ios
  integer ( kind = 4 ) line_data(line_data_num)
  integer ( kind = 4 ) line_pointer(line_num+1)
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_unit
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYL_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '".'
    stop
  end if
!
!  Write the header.
!
  call xyl_header_write ( output_filename, output_unit, point_num, line_num, &
    line_data_num )
!
!  Write the data.
!
  call xyl_data_write ( output_unit, point_num, line_num, line_data_num, &
    line_pointer, line_data )
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
