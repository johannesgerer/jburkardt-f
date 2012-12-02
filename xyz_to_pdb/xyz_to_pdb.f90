program main

!*****************************************************************************80
!
!! MAIN is the main program for XYZ_TO_PDB.
!
!  Discussion:
!
!    XYZ_TO_PDB writes the data from an XYZ file as spatial coordinates
!    of ATOM records in a PDB file.
!
!    The program may be invoked with both files specified on the command
!    line.
!
!      xyz_to_pdb file.xyz file.pdb
!
!    If either or both files are missing, the program will ask for them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) numarg
  character ( len = 255 ) pdb_file_name

  integer ( kind = 4 ) pdb_file_unit
  character ( len = 40 ) string
  integer ( kind = 4 ) xyz_file_line_num
  character ( len = 255 ) xyz_file_name
  integer ( kind = 4 ) xyz_file_unit
  integer ( kind = 4 ) xyz_num

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XYZ_TO_PDB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read XYZ data from a file.'
  write ( *, '(a)' ) '  Rewrite it as atomic coordinates in'
  write ( *, '(a)' ) '  ATOM records of a PDB file.'
  write ( *, '(a)' ) ' '
!
!  Get the number of command line arguments.
!
  numarg = iargc ( )

  if ( 1 <= numarg ) then

    iarg = 1
    call getarg ( iarg, xyz_file_name )

  else

    write ( *, '(a)' ) '  Enter the XYZ file name:'
    read ( *, '(a)', iostat = ios ) xyz_file_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XYZ_TO_PDB - Fatal error!'
      write ( *, '(a)' ) '  Could not read the XYZ file name.'
      stop
    end if

  end if

  if ( 2 <= numarg ) then

    iarg = 2
    call getarg ( iarg, pdb_file_name )

  else

    write ( *, '(a)' ) '  Enter the PDB file name:'
    read ( *, '(a)', iostat = ios ) pdb_file_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XYZ_TO_PDB - Fatal error!'
      write ( *, '(a)' ) '  Could not read the PDB file name.'
      stop
    end if

  end if

  write ( *, '(a)') 'Reading XYZ file "' // trim ( xyz_file_name ) // '".'
!
!  Open the XYZ file.
!
  call get_unit ( xyz_file_unit )

  open ( unit = xyz_file_unit, file = xyz_file_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYZ_TO_PDB - Fatal error!'
    write ( *, '(a)' ) '  Could not open the XYZ file.'
    stop
  end if
!
!  Open the PDB file.
!
  call get_unit ( pdb_file_unit )

  open ( unit = pdb_file_unit, file = pdb_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYZ_TO_PDB - Fatal error!'
    write ( *, '(a)' ) '  Could not open the PDB file.'
    stop
  end if

  call timestring ( string )
!
!  TITLE section.
!
  write ( pdb_file_unit, '(a)' ) &
    'HEADER    MUSCLE PROTEIN                          02-JUN-93   1MYS'

  write ( pdb_file_unit, '(a)' ) 'TITLE     XYZ data.'
  write ( pdb_file_unit, '(a)' ) &
    'CAVEAT     1MYS    Artificial data.  Only the ATOM XYZ is important.'
  write ( pdb_file_unit, '(a)' ) 'COMPND    MOL_ID: 1'
  write ( pdb_file_unit, '(a)' ) 'SOURCE    MOL_ID: 1;'
  write ( pdb_file_unit, '(a)' ) 'KEYWDS    NONE'
  write ( pdb_file_unit, '(a)' ) 'EXPDTA    X-RAY DIFFRACTION'
  write ( pdb_file_unit, '(a)' ) 'AUTHOR    JOHN BURKARDT'
  write ( pdb_file_unit, '(a)' ) 'REVDAT       09-JAN-06'
!
!  REMARKS section.
!
  write ( pdb_file_unit, '(a)' ) 'REMARK   1'
  write ( pdb_file_unit, '(a)' ) 'REMARK   1 REFERENCE 1'
  write ( pdb_file_unit, '(a)' ) 'REMARK   2'
  write ( pdb_file_unit, '(a)' ) 'REMARK   2 RESOLUTION. NOT APPLICABLE'
  write ( pdb_file_unit, '(a)' ) 'REMARK   3'
  write ( pdb_file_unit, '(a)' ) 'REMARK   3 REFINEMENT.'
  write ( pdb_file_unit, '(a)' ) &
    'REMARK   4 XXXX COMPLIES WITH FORMAT V. 2.1, 25-OCT-1996'
  write ( pdb_file_unit, '(a)' ) 'REMARK   5'
  write ( pdb_file_unit, '(a)' ) 'REMARK   5 WARNING.'
  write ( pdb_file_unit, '(a)' ) 'REMARK   6'
  write ( pdb_file_unit, '(a)' ) 'REMARK   6 ' // trim ( pdb_file_name )
  write ( pdb_file_unit, '(a)' ) 'REMARK   6'
  write ( pdb_file_unit, '(a)' ) 'REMARK   6 Created by XYZ_TO_PDB on'
  write ( pdb_file_unit, '(a)' ) 'REMARK   6 ' // trim ( string )
  write ( pdb_file_unit, '(a)' ) 'REMARK   6'
!
!  PRIMARY STRUCTURE section.
!
!  HETEROGEN section.
!
!  SECONDARY STRUCTURE section.
!
!  CONNECTIVITY ANNOTATION section.
!
!  MISCELLANEOUS FEATURES section.
!
!  CRYSTALLOGRAPHIC section.
!
  write ( pdb_file_unit, '(a)' ) &
    'CRYST1      1.0      1.0      1.0     90     90     90 P 1           1'
!
!  COORDINATE TRANSFORMATION section.
!
  write ( pdb_file_unit, '(a)' ) &
    'ORIGX1      1.000000  0.000000  0.000000       0.00000'
  write ( pdb_file_unit, '(a)' ) &
    'ORIGX2      0.000000  1.000000  0.000000       0.00000'
  write ( pdb_file_unit, '(a)' ) &
    'ORIGX3      0.000000  0.000000  1.000000       0.00000'
  write ( pdb_file_unit, '(a)' ) &
    'SCALE1      1.000000  0.000000  0.000000       0.00000'
  write ( pdb_file_unit, '(a)' ) &
    'SCALE2      0.000000  1.000000  0.000000       0.00000'
  write ( pdb_file_unit, '(a)' ) &
    'SCALE3      0.000000  0.000000  1.000000       0.00000'
!
!  ATOMIC COORDINATE DATA section.
!
  write ( pdb_file_unit, '(a)' ) 'MODEL        1'
!
!  Transfer the XYZ coordinate data to the PDB file as ATOM records.
!
  call xyz_to_pdb_atom ( xyz_file_unit, pdb_file_unit, xyz_file_line_num, &
    xyz_num )

  write ( pdb_file_unit, '(a)' ) 'ENMDL'
!
!  CHEMICAL CONNECTIVITY section.
!
!  BOOKKEEPING section.
!
  write ( pdb_file_unit, '(a,i5,a)' ) &
    'MASTER       15    0    0    0    0    0    0    6', &
    xyz_num, '    1    0    0' 

  write ( pdb_file_unit, '(a)' ) 'END'

  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8,a)' ) xyz_file_line_num, ' lines read from XYZ file.'
  write ( *, '(2x,i8,a)' ) xyz_num, ' coordinate records found in XYZ file.'
!
!  Close the XYZ file.
!
  close ( unit = xyz_file_unit )
!
!  Close the PDB file.
!
  close ( unit = pdb_file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XYZ_TO_PDB'
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

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

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
function s_eqi ( strng1, strng2 )

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
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
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
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LENGTH, the number of characters read
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
!    07 September 2004
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
  ierror = 0
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
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = '31 May 2001   9:45:54.872 AM'
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
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
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
  character ( len = * ) string
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

  write ( string, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine xyz_to_pdb_atom ( xyz_file_unit, pdb_file_unit, xyz_file_line_num, &
  xyz_num )

!*****************************************************************************80
!
!! XYZ_TO_PDB_ATOM_TO_XYZ writes XYZ data as ATOM records in a PDB file.
!
!  Discussion:
!
!    The XYZ and PDB files are presumed to have been opened by the user.
!
!  Format:
!
!    COLUMNS  DATA TYPE     FIELD       DEFINITION
!    --------------------------------------------------------------------------
!     1 -  6  Record name   "ATOM  "
!     7 - 11  Integer       serial      Atom serial number.
!    13 - 16  Atom          name        Atom name.
!    17       Character     altLoc      Alternate location indicator.
!    18 - 20  Residue name  resName     Residue name.
!    22       Character     chainID     Chain identifier.
!    23 - 26  Integer       resSeq      Residue sequence number.
!    27       AChar         iCode       Code for insertion of residues.
!    31 - 38  Real(8.3)     x           Orthogonal coordinates for X, Angstroms.
!    39 - 46  Real(8.3)     y           Orthogonal coordinates for Y, Angstroms.
!    47 - 54  Real(8.3)     z           Orthogonal coordinates for Z, Angstroms.
!    55 - 60  Real(6.2)     occupancy   Occupancy.
!    61 - 66  Real(6.2)     tempFactor  Temperature factor.
!    73 - 76  LString(4)    segID       Segment identifier, left-justified.
!    77 - 78  LString(2)    element     Element symbol, right-justified.
!    79 - 80  LString(2)    charge      Charge on the atom.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PDB_FILE_UNIT, the FORTRAN unit number associated with
!    the PDB file.
!
!    Input, integer ( kind = 4 ) XYZ_FILE_UNIT, the FORTRAN unit number associated with
!    the XYZ file.
!
!    Output, integer ( kind = 4 ) XYZ_FILE_LINE_NUM, the number of lines read from the
!    XYZ file.
!
!    Output, integer ( kind = 4 ) XYZ_NUM, the number of XYZ coordinates read.
!
  implicit none

  character altloc
  character ( len = 4 ) atom_name
  character chains
  character ( len = 2 ) charge
  character ( len = 2 ) element
  character ( kind = 4 ) icode
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  real ( kind = 8 ) occ
  integer ( kind = 4 ) pdb_file_unit
  character ( len = 3 ) resname
  integer ( kind = 4 ) resno
  character ( len = 4 ) segid
  character ( len = 80 ) string
  real ( kind = 8 ) temp
  real ( kind = 8 ) xyz(3)
  integer ( kind = 4 ) xyz_file_line_num
  integer ( kind = 4 ) xyz_file_unit
  integer ( kind = 4 ) xyz_num

  atom_name = 'H   '
  altloc = ' '
  resname = '   '
  chains = ' '
  resno = 0
  icode = ' '
  occ = 1.0D+00
  temp = 1.0D+00
  segid = '    '
  element = '  '
  charge = '  '

  xyz_file_line_num = 0
  xyz_num = 0

  do

    read ( xyz_file_unit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      exit
    end if

    xyz_file_line_num = xyz_file_line_num + 1
 
    if ( string(1:1) == '#' .or. len_trim ( string ) == 0 ) then

      cycle

    end if

    call s_to_r8vec ( string, 3, xyz, ierror )

    xyz_num = xyz_num + 1

    write ( pdb_file_unit, &
      '(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
      'ATOM  ', xyz_num, atom_name, altloc, resname, chains, resno, &
      icode, xyz(1:3), occ, temp, segid, element, charge

  end do
 
  write ( pdb_file_unit, &
    '(a6,i5,1x,4x,a1,a3,1x,a1,i4,a1)' ) &
    'TER   ', xyz_num+1, altloc, resname, chains, resno, icode

  return
end
