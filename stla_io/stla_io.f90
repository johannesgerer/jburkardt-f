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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

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
subroutine r8vec_cross_3d ( v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_CROSS_3D computes the cross product of two vectors in 3D.
!
!  Discussion:
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), the two vectors.
!
!    Output, real ( kind = 8 ) V3(3), the cross product vector.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v2(dim_num)
  real ( kind = 8 ) v3(dim_num)

  v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

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
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
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
!    07 September 2004
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

  logical ch_eqi
  character c
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

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

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

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
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length+1 == nchar ) then
    length = nchar
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
function stla_check ( input_file_name )

!*****************************************************************************80
!
!! STLA_CHECK checks an ASCII StereoLithography file.
!
!  Example:
!
!    solid MYSOLID
!      facet normal 0.4 0.4 0.2
!        outerloop
!          vertex  1.0 2.1 3.2
!          vertex  2.1 3.7 4.5
!          vertex  3.1 4.5 6.7
!        end loop
!      end facet
!      ...
!      facet normal 0.2 0.2 0.4
!        outerloop
!          vertex  2.0 2.3 3.4
!          vertex  3.1 3.2 6.5
!          vertex  4.1 5.5 9.0
!        end loop
!      end facet
!    end solid MYSOLID
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
!  Reference:
!
!    3D Systems, Inc,
!    Stereolithography Interface Specification,
!    October 1989.
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, logical STLA_CHECK, is TRUE if the file is legal.
!
  implicit none

  logical check
  logical done
  real ( kind = 8 ) dval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  character ( len = * ) input_file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) lchar
  logical s_eqi
  integer ( kind = 4 ) state
  logical stla_check
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num
  integer ( kind = 4 ) vertex
  character ( len = 256 ) word1
  character ( len = 256 ) word2

  state = 0
  text_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = input_file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' &
      // trim ( input_file_name ) // '".'
    stla_check = .false.
    return
  end if
!
!  Read the next line of text.
!
  do

    read ( iunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      if ( state /= 0 .and. state /= 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        write ( *, '(a)' ) '  End-of-file, but model not finished.'
        stla_check = .false.
        return
      end if
      exit
    end if

    text_num = text_num + 1

    done = .true.
!
!  Read the first word in the line.
!
    call word_next_read ( text, word1, done )

    if ( done ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
      write ( *, '(a,i8)' ) '  File line number = ', text_num
      write ( *, '(a)' ) '  No information on line.'
      stla_check = .false.
      return
    end if
!
!  "Doctor" the text, changing a beginning occurrence of:
!
!      END FACET to ENDFACET
!      END LOOP to ENDLOOP
!      END SOLID to ENDSOLID
!      FACET NORMAL to FACETNORMAL
!      OUTER LOOP to OUTERLOOP
!
    if ( s_eqi ( word1, 'END' ) ) then

      call word_next_read ( text, word2, done )

      if ( .not. s_eqi ( word2, 'FACET' ) .and. &
           .not. s_eqi ( word2, 'LOOP' ) .and. &
           .not. s_eqi ( word2, 'SOLID' ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        write ( *, '(a)' ) '  The tag END was followed by an illegal '
        write ( *, '(a)' ) '  word: "' // trim ( word2 ) // '", when expecting'
        write ( *, '(a)' ) '  "FACET", "LOOP", or "SOLID".'
        stla_check = .false.
        return
      end if

      call s_cat ( word1, word2, word1 )

    else if ( s_eqi ( word1, 'FACET' ) ) then

      call word_next_read ( text, word2, done )

      if ( .not. s_eqi ( word2, 'NORMAL' ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        write ( *, '(a)' ) '  The tag FACET was followed by an illegal '
        write ( *, '(a)' ) '  word: "' // trim ( word2 ) // '", when expecting'
        write ( *, '(a)' ) '  "NORMAL".'
        stla_check = .false.
        return
      end if

      call s_cat ( word1, word2, word1 )

    else if ( s_eqi ( word1, 'OUTER' ) ) then

      call word_next_read ( text, word2, done )

      if ( .not. s_eqi ( word2, 'LOOP' ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        write ( *, '(a)' ) '  The tag OUTER was followed by an illegal '
        write ( *, '(a)' ) '  word: "' // trim ( word2 ) // '", when expecting'
        write ( *, '(a)' ) '  "LOOP".'
        stla_check = .false.
        return
      end if

      call s_cat ( word1, word2, word1 )

    end if
!
!  This first word tells us what to do.
!
!  SOLID - begin a new solid.
!    Valid in state 0, moves to state 1.
!  ENDSOLID - end current solid.
!    Valid in state 1, moves to state 0.
!
!  FACETNORMAL - begin a new facet.
!    Valid in state 0 or 1, moves to state 2.
!  ENDFACET - end current facet.
!    Valid in state 2, moves to state 1.
!
!  OUTERLOOP - begin a list of vertices.
!    Valid in state 2, moves to state 3.
!  ENDLOOP - end vertex list.
!    Valid in state 3, moves to state 2.
!
!  VERTEX - give coordinates of next vertex.
!    Valid in state 3 if current vertex count is 0, 1 or 2.
!
!  End of file -
!    Valid in state 0 or 1.
!
    if ( s_eqi ( word1, 'SOLID' ) ) then

      if ( state /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        write ( *, '(a)' ) '  A new SOLID statement was encountered, but we'
        write ( *, '(a)' ) '  have not finished processing the current solid.'
        stla_check = .false.
        return
      end if

      state = 1

    else if ( s_eqi ( word1, 'ENDSOLID' ) ) then

      if ( state /= 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        write ( *, '(a)' ) '  An END SOLID statement was encountered, but'
        write ( *, '(a)' ) '  either we have not begun a solid at all, or we'
        write ( *, '(a)' ) '  are not at an appropriate point to finish the'
        write ( *, '(a)' ) '  current solid.'
        stla_check = .false.
        return
      end if

      state = 0

    else if ( s_eqi ( word1, 'FACETNORMAL' ) ) then

      if ( state /= 0 .and. state /= 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        write ( *, '(a)' ) '  Model not in right state for FACET.'
        stla_check = .false.
        return
      end if

      state = 2

      do i = 1, 3

        call word_next_read ( text, word2, done )

        if ( done ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
          write ( *, '(a,i8)' ) '  File line number = ', text_num
          write ( *, '(a)' ) '  End of information while reading a component'
          write ( *, '(a)' ) '  of the normal vector.'
          stla_check = .false.
          return
        end if

        call s_to_r8 ( word2, dval, ierror, lchar )

        if ( ierror /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
          write ( *, '(a,i8)' ) '  File line number = ', text_num
          write ( *, '(a)' ) &
            '  Error while reading a component of the normal vector.'
          stla_check = .false.
          return
        end if

      end do

    else if ( s_eqi ( word1, 'ENDFACET' ) ) then

      if ( state /= 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        write ( *, '(a)' ) '  Model not in right state for ENDFACET.'
        stla_check = .false.
        return
      end if

      state = 1

    else if ( s_eqi ( word1, 'OUTERLOOP' ) ) then

      if ( state /= 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        write ( *, '(a)' ) '  Model not in right state for OUTERLOOP.'
        stla_check = .false.
        return
      end if

      state = 3
      vertex = 0

    else if ( s_eqi ( word1, 'ENDLOOP' ) ) then

      if ( state /= 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        write ( *, '(a)' ) '  Model not in right state for ENDLOOP.'
        stla_check = .false.
        return
      end if

      state = 2

    else if ( s_eqi ( word1, 'VERTEX' ) ) then

      if ( state /= 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        write ( *, '(a)' ) '  Model not in right state for VERTEX.'
        stla_check = .false.
        return
      end if

      if ( 3 <= vertex ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        write ( *, '(a)' ) '  More than 3 vertices specified for a face.'
        stla_check = .false.
        return
      end if

      do i = 1, 3

        call word_next_read ( text, word2, done )

        if ( done ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
          write ( *, '(a,i8)' ) '  File line number = ', text_num
          write ( *, '(a)' ) '  The value of a vertex coordinate is missing.'
          stla_check = .false.
          return
        end if

        call s_to_r8 ( word2, dval, ierror, lchar )

        if ( ierror /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
          write ( *, '(a,i8)' ) '  File line number = ', text_num
          write ( *, '(a)' ) '  The value of a vertex coordinate makes'
          write ( *, '(a)' ) '  no sense.'
          stla_check = .false.
          return
        end if

      end do

      vertex = vertex + 1

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
      write ( *, '(a,i8)' ) '  File line number = ', text_num
      write ( *, '(a)' ) '  Unrecognized line in file.'
      stla_check = .false.
      return

    end if

  end do
!
!  Close the file.
!
  close ( unit = iunit )

  stla_check = .true.

  return
end
subroutine stla_face_node_print ( face_num, face_node )

!*****************************************************************************80
!
!! STLA_FACE_NODE_PRINT prints the node indices for each face.
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
!  Reference:
!
!    3D Systems, Inc,
!    Stereolithography Interface Specification,
!    October 1989.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NODE(3,FACE_NUM), the nodes that 
!    make up each face.
!
  implicit none

  integer ( kind = 4 ) face_num

  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_node(3,face_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Face         Nodes'
  write ( *, '(a)' ) ' '

  do face = 1, face_num

    write ( *, '(2x,i8,3(2x,i8))' ) face, face_node(1:3,face)

  end do

  return
end
subroutine stla_face_normal_compute ( node_num, face_num, node_xyz, &
  face_node, face_normal )

!*****************************************************************************80
!
!! STLA_FACE_NORMAL_COMPUTE computes normal vectors for an STLA file.
!
!  Discussion:
!
!    This routine computes the normal vector to each triangular face
!    in the STLA solid.  If the nodes of each triangular face are
!    listed in counterclockwise order (as seen from outside the solid),
!    then the normal vectors will be properly outward facing.
!
!    The normal vectors will have unit Euclidean norm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    3D Systems, Inc,
!    Stereolithography Interface Specification,
!    October 1989.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the node coordinates.
!
!    Input, integer ( kind = 4 ) FACE_NODE(3,FACE_NUM), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Output, real ( kind = 8 ) FACE_NORMAL(3,FACE_NUM), the normal vector
!    at each face.
!
  implicit none

  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_node(3,face_num)
  real ( kind = 8 ) face_normal(3,face_num)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  real ( kind = 8 ) node_xyz(3,node_num)
  real ( kind = 8 ) norm
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)

  do face = 1, face_num

    n1 = face_node(1,face)
    n2 = face_node(2,face)
    n3 = face_node(3,face)

    v1(1:3) = node_xyz(1:3,n2) - node_xyz(1:3,n1)
    v2(1:3) = node_xyz(1:3,n3) - node_xyz(1:3,n1)

    call r8vec_cross_3d ( v1, v2, face_normal(1:3,face ) )

    norm = sqrt ( sum ( ( face_normal(1:3,face) )**2 ) )

    if ( norm /= 0.0D+00 ) then
      face_normal(1:3,face) = face_normal(1:3,face) / norm
    end if

  end do

  return
end
subroutine stla_face_normal_print ( face_num, face_normal )

!*****************************************************************************80
!
!! STLA_FACE_NORMAL_PRINT prints the normal vectors.
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
!  Reference:
!
!    3D Systems, Inc,
!    Stereolithography Interface Specification,
!    October 1989.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, real ( kind = 8 ) FACE_NORMAL(3,FACE_NUM), the normal vector
!    at each face.
!
  implicit none

  integer ( kind = 4 ) face_num

  integer ( kind = 4 ) face
  real ( kind = 8 ) face_normal(3,face_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Face         Normal Vectors'
  write ( *, '(a)' ) ' '

  do face = 1, face_num

    write ( *, '(2x,i8,3(2x,g14.6))' ) face, face_normal(1:3,face)

  end do

  return
end
subroutine stla_node_xyz_print ( node_num, node_xyz )

!*****************************************************************************80
!
!! STLA_NODE_XYZ_PRINT prints the node coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    3D Systems, Inc,
!    Stereolithography Interface Specification,
!    October 1989.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the coordinates
!    of the nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xyz(3,node_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Node         Coordinates'
  write ( *, '(a)' ) ' '

  do node = 1, node_num

    write ( *, '(2x,i8,3(2x,g14.6))' ) node, node_xyz(1:3,node)

  end do

  return
end
subroutine stla_read ( input_file_name, node_num, face_num, node_xyz, &
  face_node, face_normal, ierror )

!*****************************************************************************80
!
!! STLA_READ reads graphics information from an ASCII StereoLithography file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    3D Systems, Inc,
!    Stereolithography Interface Specification,
!    October 1989.
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of vertices defined.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces defined.
!
!    Output, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the coordinates of points.
!
!    Output, integer ( kind = 4 ) FACE_NODE(3,FACE_NUM), the nodes that
!    make up each face.
!
!    Output, real ( kind = 8 ) FACE_NORMAL(3,FACE_NUM), the normal vector
!    at each face.
!
!    Output, integer ( kind = 4 ) IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) node_num

  logical done
  real ( kind = 8 ) dval
  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_node(3,face_num)
  real ( kind = 8 ) face_normal(3,face_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  character ( len = * ) input_file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xyz(3,node_num)
  logical s_eqi
  integer ( kind = 4 ) state
  real ( kind = 8 ) temp(3)
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num
  integer ( kind = 4 ) vertex
  character ( len = 256 ) word1
  character ( len = 256 ) word2

  ierror = 0
  state = 0
  text_num = 0
  face = 0
  node = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = input_file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STLA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' // &
      trim ( input_file_name ) // '".'
    ierror = 1
    return
  end if
!
!  Read the next line of text.
!
  do

    read ( iunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      if ( state /= 0 .and. state /= 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning.'
        write ( *, '(a)' ) '  End-of-file, but model not finished.'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        ierror = 1
        return
      end if
      exit
    end if

    text_num = text_num + 1
    done = .true.
!
!  Read the first word in the line.
!
    call word_next_read ( text, word1, done )
!
!  "Doctor" the text, changing a beginning occurrence of:
!
!      END FACET to ENDFACET
!      END LOOP to ENDLOOP
!      END SOLID to ENDSOLID
!      FACET NORMAL to FACETNORMAL
!      OUTER LOOP to OUTERLOOP
!
    if ( s_eqi ( word1, 'END' ) .or. &
         s_eqi ( word1, 'FACET' ) .or. &
         s_eqi ( word1, 'OUTER' ) ) then

      call word_next_read ( text, word2, done )
      call s_cat ( word1, word2, word1 )

    end if
!
!  This first word tells us what to do.
!
!  SOLID - begin a new solid.
!    Valid in state 0, moves to state 1.
!  ENDSOLID - end current solid.
!    Valid in state 1, moves to state 0.
!
!  FACETNORMAL - begin a new facet.
!    Valid in state 0 or 1, moves to state 2.
!  ENDFACET - end current facet.
!    Valid in state 2, moves to state 1.
!
!  OUTERLOOP - begin a list of vertices.
!    Valid in state 2, moves to state 3, sets vertex count to 0.
!  ENDLOOP - end vertex list.
!    Valid in state 3, moves to state 2.
!
!  VERTEX - give coordinates of next vertex.
!    Valid in state 3.
!
!  End of file -
!    Valid in state 0 or 1.
!
    if ( s_eqi ( word1, 'SOLID' ) ) then

      if ( state /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for SOLID.'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        ierror = 1
        return
      end if

      state = 1

    else if ( s_eqi ( word1, 'ENDSOLID' ) ) then

      if ( state /= 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for ENDSOLID.'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        ierror = 1
        return
      end if

      state = 0

    else if ( s_eqi ( word1, 'FACETNORMAL' ) ) then

      if ( state /= 0 .and. state /= 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for FACET.'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        ierror = 1
        return
      end if

      state = 2
      face = face + 1

      if ( face_num < face ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  More faces being read than expected.'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        ierror = 1
        return
      end if

      do i = 1, 3
        face_normal(i,face) = 0.0D+00
        call word_next_read ( text, word2, done )
        if ( .not. done ) then
          call s_to_r8 ( word2, dval, ierror, lchar )
          if ( ierror == 0 ) then
            face_normal(i,face) = dval
          end if
        end if
      end do

    else if ( s_eqi ( word1, 'ENDFACET' ) ) then

      if ( state /= 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for ENDFACET.'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        ierror = 1
        return
      end if

      state = 1

    else if ( s_eqi ( word1, 'OUTERLOOP' ) ) then

      if ( state /= 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for OUTERLOOP.'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        ierror = 1
        return
      end if

      state = 3
      vertex = 0

    else if ( s_eqi ( word1, 'ENDLOOP' ) ) then

      if ( state /= 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for ENDLOOP.'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        ierror = 1
        return
      end if

      state = 2

    else if ( s_eqi ( word1, 'VERTEX' ) ) then

      if ( state /= 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for VERTEX.'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        ierror = 1
        return
      end if

      if ( 3 <= vertex ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Too many vertices for face.'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        ierror = 1
        return
      end if

      do i = 1, 3
        call word_next_read ( text, word2, done )
        call s_to_r8 ( word2, dval, ierror, lchar )
        temp(i) = dval
      end do

      if ( node_num <= node ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  More nodes being read than expected.'
        write ( *, '(a,i8)' ) '  File line number = ', text_num
        ierror = 1
        return
      end if

      node = node + 1
      node_xyz(1:3,node) = temp(1:3)

      vertex = vertex + 1
      face_node(vertex,face) = node

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STLA_READ - Warning!'
      write ( *, '(a)' ) '  Unrecognized line in file.'
      write ( *, '(a,i8)' ) '  File line number = ', text_num
      ierror = 1
      return

    end if

  end do
!
!  Close the file.
!
  close ( unit = iunit )

  return
end
subroutine stla_size ( input_file_name, solid_num, node_num, face_num, &
  text_num )

!*****************************************************************************80
!
!! STLA_SIZE determines sizes associated with an STLA file.
!
!  Discussion:
!
!    This routine assumes that the file is a legal STLA file.
!
!    To perform checks on the file, call STLA_CHECK first.
!
!    Note that the counts for the number of nodes and edges are
!    overestimates, since presumably, most nodes will be defined several
!    times, once for each face they are part of, and most edges will
!    be defined twice.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    3D Systems, Inc,
!    Stereolithography Interface Specification,
!    October 1989.
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) SOLID_NUM, the number of solids defined.
!    Presumably, this is 1.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of vertices defined.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of faces defined.
!
!    Output, integer ( kind = 4 ) TEXT_NUM, the number of lines of text in the file.
!
  implicit none

  logical check
  logical done
  real ( kind = 8 ) dval
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  character ( len = * ) input_file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) node_num
  logical s_eqi
  integer ( kind = 4 ) solid_num
  integer ( kind = 4 ) state
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num
  integer ( kind = 4 ) vertex
  character ( len = 256 ) word1
  character ( len = 256 ) word2

  ierror = 0

  state = 0
  text_num = 0

  solid_num = 0
  node_num = 0
  face_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = input_file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STLA_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' &
      // trim ( input_file_name ) // '".'
    ierror = 1
    return
  end if
!
!  Read the next line of text.
!
  do

    read ( iunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      if ( state /= 0 .and. state /= 1 ) then
        return
      end if
      exit
    end if

    text_num = text_num + 1

    done = .true.
!
!  Read the first word in the line.
!
    call word_next_read ( text, word1, done )

    if ( done ) then
      return
    end if
!
!  "Doctor" the text, changing a beginning occurrence of:
!
!      END FACET to ENDFACET
!      END LOOP to ENDLOOP
!      END SOLID to ENDSOLID
!      FACET NORMAL to FACETNORMAL
!      OUTER LOOP to OUTERLOOP
!
    if ( s_eqi ( word1, 'END' ) ) then

      call word_next_read ( text, word2, done )

      if ( .not. s_eqi ( word2, 'FACET' ) .and. &
           .not. s_eqi ( word2, 'LOOP' ) .and. &
           .not. s_eqi ( word2, 'SOLID' ) ) then
        return
      end if

      call s_cat ( word1, word2, word1 )

    else if ( s_eqi ( word1, 'FACET' ) ) then

      call word_next_read ( text, word2, done )

      if ( .not. s_eqi ( word2, 'NORMAL' ) ) then
        return
      end if

      call s_cat ( word1, word2, word1 )

    else if ( s_eqi ( word1, 'OUTER' ) ) then

      call word_next_read ( text, word2, done )

      if ( .not. s_eqi ( word2, 'LOOP' ) ) then
        return
      end if

      call s_cat ( word1, word2, word1 )

    end if
!
!  This first word tells us what to do.
!
!  SOLID - begin a new solid.
!    Valid in state 0, moves to state 1.
!  ENDSOLID - end current solid.
!    Valid in state 1, moves to state 0.
!
!  FACETNORMAL - begin a new facet.
!    Valid in state 0 or 1, moves to state 2.
!  ENDFACET - end current facet.
!    Valid in state 2, moves to state 1.
!
!  OUTERLOOP - begin a list of vertices.
!    Valid in state 2, moves to state 3.
!  ENDLOOP - end vertex list.
!    Valid in state 3, moves to state 2.
!
!  VERTEX - give coordinates of next vertex.
!    Valid in state 3.
!
!  End of file -
!    Valid in state 0 or 1.
!
    if ( s_eqi ( word1, 'SOLID' ) ) then

      if ( state /= 0 ) then
        return
      end if

      state = 1

    else if ( s_eqi ( word1, 'ENDSOLID' ) ) then

      if ( state /= 1 ) then
        return
      end if

      state = 0

      solid_num = solid_num + 1

    else if ( s_eqi ( word1, 'FACETNORMAL' ) ) then

      if ( state /= 0 .and. state /= 1 ) then
        return
      end if

      state = 2

      do i = 1, 3

        call word_next_read ( text, word2, done )

        if ( done ) then
          return
        end if

        call s_to_r8 ( word2, dval, ierror, lchar )

        if ( ierror /= 0 ) then
          return
        end if

      end do

    else if ( s_eqi ( word1, 'ENDFACET' ) ) then

      if ( state /= 2 ) then
        return
      end if

      state = 1

      face_num = face_num + 1

    else if ( s_eqi ( word1, 'OUTERLOOP' ) ) then

      if ( state /= 2 ) then
        return
      end if

      state = 3
      vertex = 0

    else if ( s_eqi ( word1, 'ENDLOOP' ) ) then

      if ( state /= 3 ) then
        return
      end if

      state = 2

    else if ( s_eqi ( word1, 'VERTEX' ) ) then

      if ( state /= 3 ) then
        return
      end if

      if ( 3 <= vertex ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Too many vertices for a face.'
        ierror = 1
        return
      end if

      do i = 1, 3

        call word_next_read ( text, word2, done )

        if ( done ) then
          return
        end if

        call s_to_r8 ( word2, dval, ierror, lchar )

        if ( ierror /= 0 ) then
          return
        end if

      end do

      vertex = vertex + 1
      node_num = node_num + 1

    else

      return

    end if

  end do
!
!  Close the file.
!
  close ( unit = iunit )

  return
end
subroutine stla_size_print ( input_file_name, solid_num, node_num, face_num, &
  text_num )

!*****************************************************************************80
!
!! STLA_SIZE_PRINT prints sizes associated with an STLA file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    3D Systems, Inc,
!    Stereolithography Interface Specification,
!    October 1989.
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) SOLID_NUM, the number of solids defined.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of vertices defined.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces defined.
!
!    Input, integer ( kind = 4 ) TEXT_NUM, the number of lines of text in the file.
!
  implicit none

  integer ( kind = 4 ) face_num
  character ( len = * ) input_file_name
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) solid_num
  integer ( kind = 4 ) text_num

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  Object sizes for STLA file "' // &
    trim ( input_file_name ) // '":'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Solids =                   ', solid_num
  write ( *, '(a,i8)' ) '  Nodes (may be repeated) =  ', node_num
  write ( *, '(a,i8)' ) '  Faces (triangular only) =  ', face_num
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of lines of text =  ', text_num

  return
end
subroutine stla_write ( output_file_name, node_num, face_num, node_xyz, &
  face_node, face_normal )

!*****************************************************************************80
!
!! STLA_WRITE writes graphics information to an ASCII StereoLithography file.
!
!  Example:
!
!    solid MYSOLID
!      facet normal 0.4 0.4 0.2
!        outerloop
!          vertex  1.0 2.1 3.2
!          vertex  2.1 3.7 4.5
!          vertex  3.1 4.5 6.7
!        end loop
!      end facet
!      ...
!      facet normal 0.2 0.2 0.4
!        outerloop
!          vertex  2.0 2.3 3.4
!          vertex  3.1 3.2 6.5
!          vertex  4.1 5.5 9.0
!        end loop
!      end facet
!    end solid MYSOLID
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    3D Systems, Inc,
!    Stereolithography Interface Specification,
!    October 1989.
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the node coordinates.
!
!    Input, integer ( kind = 4 ) FACE_NODE(3,FACE_NUM), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, real ( kind = 8 ) FACE_NORMAL(3,FACE_NUM), the normal vector
!    at each face.
!
  implicit none

  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_node(3,face_num)
  real ( kind = 8 ) face_normal(3,face_num)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xyz(3,node_num)
  character ( len = * ) output_file_name
  integer ( kind = 4 ) text_num
  integer ( kind = 4 ) vertex

  text_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = output_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STLA_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' &
      // trim ( output_file_name ) // '".'
    stop
  end if

  write ( iunit, '(a)' ) 'solid MYSOLID'
  text_num = text_num + 1

  do face = 1, face_num

    write ( iunit, '(a,3(2x,g14.6))' ) '  facet normal', face_normal(1:3,face)
    text_num = text_num + 1

    write ( iunit, '(a)' ) '    outer loop'
    text_num = text_num + 1

    do vertex = 1, 3

      node = face_node(vertex,face)

      write ( iunit, '(a,2x,3(2x,g14.6))' ) '      vertex', node_xyz(1:3,node)
      text_num = text_num + 1

    end do

    write ( iunit, '(a)' ) '    end loop'
    text_num = text_num + 1
    write ( iunit, '(a)' ) '  end facet'
    text_num = text_num + 1

  end do

  write ( iunit, '(a)' ) 'end solid MYSOLID'
  text_num = text_num + 1

  close ( unit = iunit )

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
subroutine word_next_read ( s, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_READ "reads" words from a string, one at a time.
!
!  Special cases:
!
!    The following characters are considered to be a single word,
!    whether surrounded by spaces or not:
!
!      " ( ) { } [ ]
!
!    Also, if there is a trailing comma on the word, it is stripped off.
!    This is to facilitate the reading of lists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string, presumably containing words
!    separated by spaces.
!
!    Output, character ( len = * ) WORD.
!    If DONE is FALSE, then WORD contains the "next" word read.
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!    On input with a fresh string, set DONE to TRUE.
!    On output, the routine sets DONE:
!      FALSE if another word was read,
!      TRUE if no more words could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), save :: lenc = 0
  integer ( kind = 4 ), save :: next = 1
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
  character ( len = * ) word
!
!  We "remember" LENC and NEXT from the previous call.
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
  if ( done ) then

    next = 1
    done = .false.
    lenc = len_trim ( s )

    if ( lenc <= 0 ) then
      done = .true.
      word = ' '
      return
    end if

  end if
!
!  Beginning at index NEXT, search the string for the next nonblank,
!  which signals the beginning of a word.
!
  ilo = next
!
!  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  do

    if ( lenc < ilo ) then
      word = ' '
      done = .true.
      next = lenc + 1
      return
    end if
!
!  If the current character is blank, skip to the next one.
!
    if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
  if ( s(ilo:ilo) == '"' .or. &
       s(ilo:ilo) == '(' .or. &
       s(ilo:ilo) == ')' .or. &
       s(ilo:ilo) == '{' .or. &
       s(ilo:ilo) == '}' .or. &
       s(ilo:ilo) == '[' .or. &
       s(ilo:ilo) == ']' ) then

    word = s(ilo:ilo)
    next = ilo + 1
    return

  end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
  next = ilo + 1

  do while ( next <= lenc )

    if ( s(next:next) == ' ' ) then
      exit
    else if ( s(next:next) == TAB ) then
      exit
    else if ( s(next:next) == '"' ) then
      exit
    else if ( s(next:next) == '(' ) then
      exit
    else if ( s(next:next) == ')' ) then
      exit
    else if ( s(next:next) == '{' ) then
      exit
    else if ( s(next:next) == '}' ) then
      exit
    else if ( s(next:next) == '[' ) then
      exit
    else if ( s(next:next) == ']' ) then
      exit
    end if

    next = next + 1

  end do

  if ( s(next-1:next-1) == ',' ) then
    word = s(ilo:next-2)
  else
    word = s(ilo:next-1)
  end if

  return
end
