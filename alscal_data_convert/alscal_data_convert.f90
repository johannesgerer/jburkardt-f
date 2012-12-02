program main

!*****************************************************************************80
!
!! MAIN is the main program for ALSCAL_DATA_CONVERT.
!
!  Discussion:
!
!    ALSCAL_DATA_CONVERT converts a data file to ALSCAL format.
!
!    The input data file has a fairly simple form.  The output
!    data file tries to write an appropriate value for every ALSCAL option.
!
!  Modified:
!
!    05 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: maxn = 100

  integer   ( kind = 4 ) ierror
  character ( len = 80 ) input_filename
  character ( len = 10 ) name(maxn)
  integer   ( kind = 4 ) ncol
  integer   ( kind = 4 ) nrow
  character ( len = 80 ) output_filename
  real      ( kind = 8 ) xx(maxn,maxn)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ALSCAL_DATA_CONVERT'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read Hugh''s data, '
  write ( *, '(a)' ) '  Print it out,'
  write ( *, '(a)' ) '  Write it to a file for ALSCAL.'
  write ( *, '(a)' ) ' '
!
!  Read Hugh's data.
!
  input_filename = 'input_data.txt'

  call input_read ( input_filename, maxn, nrow, ncol, xx, name, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ALSCAL_DATA_CONVERT'
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
  output_filename = 'alscal_data.txt'

  call output_write ( output_filename, maxn, nrow, ncol, xx, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ALSCAL_DATA_CONVERT'
    write ( *, '(a)' ) '  Failure while writing the output data.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ALSCAL_DATA_CONVERT'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
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
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is TRUE.
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
subroutine i4_extract ( s, i, ierror )

!*****************************************************************************80
!
!! I4_EXTRACT "extracts" an I4 from the beginning of a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Output, integer ( kind = 4 ) I.  If IERROR is 0, then I contains the
!    next integer read from S; otherwise I is 0.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, no error.
!    nonzero, an integer could not be extracted from the beginning of the
!    string.  I is 0 and S is unchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  character ( len = * ) s

  i = 0

  call s_to_i4 ( s, i, ierror, length )

  if ( ierror /= 0 .or. length == 0 ) then
    ierror = 1
    i = 0
  else
    call s_shift_left ( s, length )
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
!    05 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXN, the maximum values for NROW and NCOL.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows of data.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns of data.
!
!    Input, real ( kind = 8 ) XX(MAXN,MAXN), the NROW by NCOL distance data.
!
!    Input, character ( len = 10 ) NAME(MAXN), the names of the objects.
!
  implicit none

  integer   ( kind = 4 ) maxn

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = 10 ) name(maxn)
  integer   ( kind = 4 ) ncol
  integer   ( kind = 4 ) nrow
  real      ( kind = 8 ) xx(maxn,maxn)

  write ( *, * ) ' '
  write ( *, * ) 'Object names:'
  write ( *, * ) ' '
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
!
!*****************************************************************************80
!
!! INPUT_READ reads the distance data from a file.
!
!
!  Modified:
!
!    29 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 80 ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) MAXN, the maximum values for NROW and NCOL.
!
!    Output, integer ( kind = 4 ) NROW, the number of rows of data.
!
!    Output, integer ( kind = 4 ) NCOL, the number of columns of data.
!
!    Output, real ( kind = 8 ) XX(MAXN,MAXN), the NROW by NCOL distance data.
!
!    Output, character ( len = 10 ) NAME(MAXN), the names of the objects.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    nonzero, an error occurred.
!
  implicit none

  integer   ( kind = 4 ) maxn

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  character ( len = 80 ) input_filename
  integer   ( kind = 4 ) input_unit
  integer   ( kind = 4 ) ios
  integer   ( kind = 4 ) j
  character ( len = 10 ) name(maxn)
  integer   ( kind = 4 ) ncol
  integer   ( kind = 4 ) nrow
  character ( len = 200 ) s
  real      ( kind = 8 ) xx(maxn,maxn)

  ierror = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INPUT_READ'
  write ( *, '(a)' ) '  Reading data from "' // trim ( input_filename ) // '".'

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INPUT_READ - Fatal error'
    write ( *, '(a)' ) '  Could not open the file!'
    return
  end if

  read ( input_unit, '(a)', iostat = ios ) s

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INPUT_READ - Fatal error'
    write ( *, '(a)' ) '  A read error occurred!'
    return
  end if

  call i4_extract ( s, nrow, ierror )
  ncol = nrow

  do i = 1, nrow

    read ( input_unit, '(a)', iostat = ios ) s

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INPUT_READ - Fatal error'
      write ( *, '(a)' ) '  A read error occurred!'
      return
    end if

    call s_word_extract_first ( s, name(i) )

    do j = 1, ncol

      do

        call r8_extract ( s, xx(i,j), ierror )

        if ( ierror == 0 ) then
          exit
        end if

        read ( input_unit, '(a)', iostat = ios ) s

        if ( ios /= 0 ) then
          ierror = ios
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'INPUT_READ - Fatal error'
          write ( *, '(a)' ) '  A read error occurred!'
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
!! OUTPUT_WRITE writes the distance data to a file in ALSCAL format.
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
!    Input, integer ( kind = 4 ) MAXN, the maximum values for NROW and NCOL.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows of data.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns of data.
!
!    Input, real ( kind = 8 ) XX(MAXN,MAXN), the NROW by NCOL distance data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    nonzero, an error occurred.
!
  implicit none

  integer   ( kind = 4 ) maxn

  real      ( kind = 8 ) cut
  integer   ( kind = 4 ) debug
  real      ( kind = 8 ) epsi
  character ( len = 80 ) fmt
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) icnstr
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) indata
  integer   ( kind = 4 ) initw
  integer   ( kind = 4 ) initws
  integer   ( kind = 4 ) initx
  integer   ( kind = 4 ) initxc
  character ( len = 80 ) input_filename
  integer   ( kind = 4 ) ios
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) maxit
  character ( len = 10 ) name(maxn)
  integer   ( kind = 4 ) ncc
  integer   ( kind = 4 ) ncol
  integer   ( kind = 4 ) ndeg
  integer   ( kind = 4 ) ndim
  integer   ( kind = 4 ) ndir
  integer   ( kind = 4 ) ndmn
  integer   ( kind = 4 ) ndmx
  integer   ( kind = 4 ) ndt
  integer   ( kind = 4 ) ndtyp
  integer   ( kind = 4 ) noulb
  integer   ( kind = 4 ) nph
  integer   ( kind = 4 ) nps
  integer   ( kind = 4 ) npt
  integer   ( kind = 4 ) nrow
  integer   ( kind = 4 ) ns
  integer   ( kind = 4 ) nsim
  integer   ( kind = 4 ) nwc
  integer   ( kind = 4 ) nwe
  character ( len = 80 ) output_filename
  integer   ( kind = 4 ) output_unit
  real      ( kind = 8 ) stmin
  real      ( kind = 8 ) xx(maxn,maxn)

  ierror = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OUTPUT_WRITE'
  write ( *, '(a)' ) '  Writing data to "' // trim ( output_filename ) // '".'

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

  write ( output_unit, '(a)' ) 'Hugh''s Data'

  ns = 1
  ndtyp = 2
  nsim = 0
  nps = 1
  nwc = 1
  ndeg = 1
  ndmx = 999
  cut = 0.0D+00

  write ( output_unit, '(9i4,f8.4)' ) nrow, ncol, ns, ndtyp, nsim, nps, &
    nwc, ndeg, ndmx, cut

  nwe = 0
  ndim = 2
  ndmn = 1
  ncc = 0
  maxit = 30
  epsi = 0.001D+00
  stmin = 0.005D+00
  ndir = 1

  write ( output_unit, '(5i4,2f8.4,i4)' ) nwe, ndim, ndmn, ncc, maxit, epsi, &
    stmin, ndir

  ndt = 1
  npt = 1
  nph = 0
  indata = 0
  initx = 1
  initxc = 1
  initw = 1
  initws = 1
  noulb = 0
  icnstr = 0
  debug = 0

  write ( output_unit, '(12i4)' ) ndt, npt, nph, indata, initx, initxc, &
    initw, initws, noulb, icnstr, debug

  fmt = '(10f8.3)'
  write ( output_unit, '(a)' ) fmt

  do i = 1, nrow
    write ( output_unit, fmt ) ( xx(i,j), j = 1, ncol )
  end do

  write ( output_unit, '(a)' ) 'END'

  close ( unit = output_unit )

  return
end
subroutine r8_extract ( s, r8, ierror )

!*****************************************************************************80
!
!! R8_EXTRACT "extracts" an R8 from the beginning of a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Output, real ( kind = 8 ) R8.  If IERROR is 0, then R4 contains the
!    next real read from the string; otherwise R4 is 0.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, no error.
!    nonzero, a real could not be extracted from the beginning of the
!    string.  R4 is 0.0 and S is unchanged.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  real ( kind = 8 ) r8
  character ( len = * ) s

  r8 = 0.0D+00

  call s_to_r8 ( s, r8, ierror, length )

  if ( ierror /= 0 .or. length == 0 ) then
    ierror = 1
    r8 = 0.0D+00
  else
    call s_shift_left ( s, length )
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, integer ( kind = 4 ) ISHFT, the number of positions to the
!    left to shift the characters.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ishft
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len ( s )

  if ( 0 < ishft ) then

    do i = 1, s_length - ishft
      s(i:i) = s(i+ishft:i+ishft)
    end do

    do i = s_length - ishft + 1, s_length
      s(i:i) = ' '
    end do

  else if ( ishft < 0 ) then

    do i = s_length, - ishft + 1, - 1
      s(i:i) = s(i+ishft:i+ishft)
    end do

    do i = -ishft, 1, -1
      s(i:i) = ' '
    end do

  end if

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
!    Output, integer ( kind = 4 ) LENGTH, the number of characters
!    of S used to make the integer.
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
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' ) then

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
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  integer ( kind = 4 ) s_length

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
    if ( .false. ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
      write ( *, '(a)' ) '  Illegal or nonnumeric input:'
      write ( *, '(a)' ) '    "' // trim ( s ) // '".'
    end if
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
subroutine s_word_extract_first ( s, w )

!*****************************************************************************80
!
!! S_WORD_EXTRACT_FIRST extracts the first word from a string.
!
!  Discussion:
!
!    A "word" is a string of characters terminated by a blank or
!    the end of the string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 2006
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

  integer ( kind = 4 ) get1
  integer ( kind = 4 ) get2
  character ( len = * ) s
  integer ( kind = 4 ) s_length
  character ( len = * ) w

  w = ' '

  s_length = len_trim ( s )

  if ( s_length < 1 ) then
    return
  end if
!
!  Find the first nonblank.
!
  get1 = 0

  do

    get1 = get1 + 1

    if ( s_length < get1 ) then
      return
    end if

    if ( s(get1:get1) /= ' ' ) then
      exit
    end if

  end do
!
!  Look for the last contiguous nonblank.
!
  get2 = get1

  do

    if ( s_length <= get2 ) then
      exit
    end if

    if ( s(get2+1:get2+1) == ' ' ) then
      exit
    end if

    get2 = get2 + 1

  end do
!
!  Copy the word.
!
  w = s(get1:get2)
!
!  Shift the string.
!
  s(1:get2) = ' '
  s = adjustl ( s(get2+1:) )

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
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
