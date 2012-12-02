program main

!*****************************************************************************80
!
!! MAIN is the main program for SPACER_DATA_CONVERT.
!
!  Discussion:
!
!    SPACER_DATA_CONVERT converts a data file to SPACER format.
!
!    The input file has a fairly simple form.
!
!  Usage:
!
!    spacer_data_convert input_name output_name
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: maxn = 100

  integer iarg
  integer iargc
  integer ierror
  integer ilen
  character ( len = 255 ) input_file_name
  character ( len = 10 ) name(maxn)
  integer ncol
  integer nrow
  integer numarg
  character ( len = 255 ) output_file_name
  real xx(maxn,maxn)

  ierror = 0
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPACER_DATA_CONVERT'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Convert input data to SPACER format.'
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
    write ( *, '(a)' ) 'SPACER_DATA_CONVERT - Input required!'
    write ( *, '(a)' ) '  Enter the input filename.'
    read ( *, '(a)' ) input_file_name

  else

    iarg = 1
!
!  Old style:
!
    call getarg ( iarg, input_file_name )
!
!  New style:
!
!   call pxfgetarg ( iarg, input_file_name, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'SPACER_DATA_CONVERT - Error!'
!     write ( *, '(a)' ) '  Could not read the argument.'
!     stop
!   end if

  end if
!
!  Get the output filename.
!
  if ( numarg < 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPACER_DATA_CONVERT - Input required!'
    write ( *, '(a)' ) '  Enter the output filename.'
    read ( *, '(a)' ) output_file_name

  else

    iarg = 2
!
!  Old style:
!
    call getarg ( iarg, output_file_name )
!
!  New style:
!
!   call pxfgetarg ( iarg, output_file_name, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'SPACER_DATA_CONVERT - Error!'
!     write ( *, '(a)' ) '  Could not read the argument.'
!     stop
!   end if

  end if

  call input_read ( input_file_name, maxn, nrow, ncol, xx, name, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPACER_DATA_CONVERT'
    write ( *, '(a)' ) '  Failure while reading the input data.'
    stop
  end if
!
!  Print the data.
!
  call input_print ( maxn, nrow, ncol, xx, name )
!
!  Write the data to a file.
!
  call output_write ( output_file_name, maxn, nrow, ncol, xx, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPACER_DATA_CONVERT'
    write ( *, '(a)' ) '  Failure while writing the output data.'
    stop
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPACER_DATA_CONVERT'
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
  integer itemp

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
!  Examples:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
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
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer digit

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
!    Output, integer IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer i
  integer ios
  integer iunit
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
subroutine i4_extract ( s, i, ierror )

!*****************************************************************************80
!
!! I4_EXTRACT "extracts" an I4 from the beginning of a string.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S; on input, a string from
!    whose beginning an integer is to be extracted.  On output,
!    the integer, if found, has been removed.
!
!    Output, integer I.  If IERROR is 0, then I contains the
!    "next" integer read from S; otherwise I is 0.
!
!    Output, integer IERROR.
!    0, no error.
!    nonzero, an integer could not be extracted from the beginning of the
!    string.  I is 0 and S is unchanged.
!
  implicit none

  integer i
  integer ierror
  integer lchar
  character ( len = * ) s

  i = 0

  call s_to_i4 ( s, i, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    ierror = 1
    i = 0
  else
    call s_shift_left ( s, lchar )
  end if

  return
end
subroutine input_print ( maxn, nrow, ncol, xx, name )

!*****************************************************************************80
!
!! INPUT_PRINT prints the input data.
!
!  Modified:
!
!    03 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXN, the maximum values for NROW and NCOL.
!
!    Input, integer NROW, the number of rows of data.
!
!    Input, integer NCOL, the number of columns of data.
!
!    Input, real XX(MAXN,MAXN), the NROW by NCOL distance data.
!
!    Input, character ( len = 10 ) NAME(MAXN), the names of the objects.
!
  implicit none

  integer maxn

  integer i
  integer j
  integer jhi
  integer jlo
  character ( len = 10 ) name(maxn)
  integer ncol
  integer nrow
  real xx(maxn,maxn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Object names:'
  write ( *, '(a)' ) ' '
  do i = 1, nrow
    write ( *, '(i6,2x,a10)' ) i, name(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Object distances:'
  write ( *, '(a)' ) ' '
  do i = 1, nrow
    do jlo = 1, ncol, 10
      jhi = min ( jlo + 9, ncol )
      if ( jlo == 1 ) then
        write ( *, '(i6,2x,10f8.3)' ) i, xx(i,jlo:jhi)
      else
        write ( *, '(6x,2x,10f8.3)' ) xx(i,jlo:jhi)
      end if
    end do
  end do

  return
end
subroutine input_read ( input_filename, maxn, nrow, ncol, xx, name, ierror )

!*****************************************************************************80
!
!! INPUT_READ reads the distance data from a file.
!
!  Modified:
!
!    05 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 80 ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer MAXN, the maximum values for NROW and NCOL.
!
!    Output, integer NROW, the number of rows of data.
!
!    Output, integer NCOL, the number of columns of data.
!
!    Output, real XX(MAXN,MAXN), the NROW by NCOL distance data.
!
!    Output, character ( len = 10 ) NAME(MAXN), the names of the objects.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, an error occurred.
!
  implicit none

  integer maxn

  integer i
  integer ierror
  character ( len = 80 ) input_filename
  integer input_unit
  integer ios
  integer j
  character ( len = 10 ) name(maxn)
  integer ncol
  integer nrow
  character ( len = 200 ) s
  real xx(maxn,maxn)

  ierror = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INPUT_READ'
  write ( *, '(a)' ) '  Reading data from ' // trim ( input_filename )

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INPUT_READ'
    write ( *, '(a)' ) '  An error occurred!'
    return
  end if

  read ( input_unit, '(a)', iostat = ios ) s

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INPUT_READ'
    write ( *, '(a)' ) '  An error occurred!'
    return
  end if

  call i4_extract ( s, nrow, ierror )
  ncol = nrow

  do i = 1, nrow

    read ( input_unit, '(a)', iostat = ios ) s

    if ( ios /= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INPUT_READ'
      write ( *, '(a)' ) '  An error occurred!'
      return
    end if

    call word_extract ( s, name(i) )

    do j = 1, ncol

      do

        call r4_extract ( s, xx(i,j), ierror )

        if ( ierror == 0 ) then
          exit
        end if

        read ( input_unit, '(a)', iostat = ios ) s

        if ( ios /= 0 ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'INPUT_READ'
          write ( *, '(a)' ) '  An error occurred!'
          return
        end if

      end do

    end do

  end do

  close ( unit = input_unit )

  return
end
subroutine output_write ( output_filename, maxn, nrow, ncol, xx, ierror )

!*****************************************************************************80
!
!! OUTPUT_WRITE writes the distance data to a file in SPACER format.
!
!  Modified:
!
!    03 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 80 ) OUTPUT_FILENAME, the name of the output file.
!
!    Input, integer MAXN, the maximum values for NROW and NCOL.
!
!    Input, integer NROW, the number of rows of data.
!
!    Input, integer NCOL, the number of columns of data.
!
!    Input, real XX(MAXN,MAXN), the NROW by NCOL distance data.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, an error occurred.
!
  implicit none

  integer maxn

  integer i
  integer ierror
  integer ios
  integer j
  integer ncol
  integer nrow
  character ( len = 80 ) output_filename
  integer output_unit
  real xx(maxn,maxn)

  ierror = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OUTPUT_WRITE'
  write ( *, '(a)' ) '  Writing data to ' // trim ( output_filename )

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OUTPUT_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( output_unit, '(i6)' ) nrow

  do i = 1, nrow
    write ( output_unit, '(10f8.3)' ) xx(i,1:i)
  end do

  close ( unit = output_unit )

  return
end
subroutine r4_extract ( s, r, ierror )

!*****************************************************************************80
!
!! R4_EXTRACT "extracts" an R4 from the beginning of a string.
!
!  Modified:
!
!    02 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S; on input, a string from
!    whose beginning a real is to be extracted.  On output,
!    the real, if found, has been removed.
!
!    Output, real R.  If IERROR is 0, then R contains the
!    next real read from the string; otherwise R is 0.
!
!    Output, integer IERROR.
!    0, no error.
!    nonzero, a real could not be extracted from the beginning of the
!    string.  R is 0.0 and S is unchanged.
!
  implicit none

  integer ierror
  integer lchar
  real r
  character ( len = * ) s

  r = 0.0E+00

  call s_to_r4 ( s, r, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    ierror = 1
    r = 0.0E+00
  else
    call s_shift_left ( s, lchar )
  end if

  return
end
subroutine s_shift_left ( s, ishft )

!*****************************************************************************80
!
!! S_SHIFT_LEFT shifts the characters in a string to the left and blank pads.
!
!  Discussion:
!
!    A shift of 2 would change "Violin" to "olin  ".
!    A shift of -2 would change "Violin" to "  Violin".
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be shifted.
!
!    Input, integer ISHFT, the number of positions to the
!    left to shift the characters.
!
  implicit none

  integer i
  integer ishft
  integer nchar
  character ( len = * ) s

  nchar = len ( s )

  if ( ishft > 0 ) then

    do i = 1, nchar - ishft
      s(i:i) = s(i+ishft:i+ishft)
    end do

    do i = nchar - ishft + 1, nchar
      s(i:i) = ' '
    end do

  else if ( ishft < 0 ) then

    do i = nchar, - ishft + 1, - 1
      s(i:i) = s(i+ishft:i+ishft)
    end do

    do i = - ishft, 1, -1
      s(i:i) = ' '
    end do

  end if

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
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
!    Output, integer IVAL, the integer value read from the string.
!    If STRING is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer last
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
subroutine s_to_r4 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R4 reads an R4 from a string.
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
!  Examples:
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
!    Output, real R, the real value that was read from the string.
!
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer lchar
  integer nchar
  integer ndig
  real r
  real rbot
  real rexp
  real rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0E+00
  lchar = - 1
  isgn = 1
  rtop = 0.0E+00
  rbot = 1.0E+00
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
        rtop = 10.0E+00 * rtop + real ( ndig )
      else if ( ihave == 5 ) then
        rtop = 10.0E+00 * rtop + real ( ndig )
        rbot = 10.0E+00 * rbot
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
    rexp = 1.0E+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0E+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0E+00**rexp
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
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
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
subroutine word_extract ( s, w )

!*****************************************************************************80
!
!! WORD_EXTRACT extracts the next word from a string.
!
!  Discussion:
!
!    A "word" is a string of characters terminated by a blank or
!    the end of the string.
!
!  Modified:
!
!    31 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string.  On output, the first
!    word has been removed, and the remaining string has been shifted left.
!
!    Output, character ( len = * ) W, the leading word of the string.
!
  implicit none

  integer iget1
  integer iget2
  integer lchar
  character ( len = * ) s
  character ( len = * ) w

  w = ' '

  lchar = len_trim ( s )
!
!  Find the first nonblank.
!
  iget1 = 0

  do

    iget1 = iget1 + 1

    if ( iget1 > lchar ) then
      return
    end if

    if ( s(iget1:iget1) /= ' ' ) then
      exit
    end if

  end do
!
!  Look for the last contiguous nonblank.
!
  iget2 = iget1

  do

    if ( iget2 >= lchar ) then
      exit
    end if

    if ( s(iget2+1:iget2+1) == ' ' ) then
      exit
    end if

    iget2 = iget2 + 1

  end do
!
!  Copy the word.
!
  w = s(iget1:iget2)
!
!  Shift the string.
!
  s(1:iget2) = ' '
  s = adjustl ( s(iget2+1:) )

  return
end
