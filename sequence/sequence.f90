program main

!*****************************************************************************80
!
!! MAIN is the main program for SEQUENCE.
!
!  Discussion:
!
!    SEQUENCE reads in a sequence and calls SEQUENCE_FILLER to fill it in.
!
!  Discussion:
!
!    The input to the program is a string representing a sequence,
!    with missing values.  For instance,
!
!      '? 3 6 ? 15'
!
!    The program assumes that the question marks represent missing
!    values that are to be filled in once a rule has been deduced
!    for the sequence.
!
!    The rule deduced for the sequence is found by constructing the
!    polynomial that interpolates the known data.  Abscissas for the
!    data are assigned by position in the sequence.  Thus for the
!    above sequence, the three known values have abscissas of 2, 3
!    and 5.
!
!    Once the interpolating polynomial is found, it is evaluated at
!    the points where data was not given, and the results are reported
!    back to the user.  In the example case, we would get the
!    filled in sequence:
!
!      '1 3 6 10 15'
!
!    Note that this procedure will always be able to produce a result,
!    but it may not be the expected result.  This is particularly so
!    when the sequence is most easily represented in terms of a geometric
!    calculation, as in
!
!      '1 2 4 8 16 32'
!
!    or as in the Fibonacci sequence:
!
!      '1 1 2 3 5 8 13 21 34 55'
!
!    both of which this program will be unable to recognize.
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
  implicit none

  integer ( kind = 4 ) ios
  integer ( kind = 4 ) nseq
  character ( len = 80 ) string

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SEQUENCE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Fill in missing entries in a numeric sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Specify a sequence problem by listing the'
  write ( *, '(a)' ) '  numbers, using "?" for missing values, as in:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    1 2 ? 4 5 ?'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program will print back its best guess for'
  write ( *, '(a)' ) '  the sequence; in this case, of course:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    1 2 3 4 5 6'

  nseq = 0

  do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter a sequence (or hit RETURN to quit)'

    read ( *, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( string ) == 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your sequence is ' // trim ( string )

    call sequence_filler ( string )
    nseq = nseq + 1

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of sequences processed: ', nseq
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SEQUENCE:'
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
!    Output, integer DIGIT, the corresponding integer value.  If C was
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
subroutine data_to_dif ( diftab, ntab, xtab, ytab )

!*****************************************************************************80
!
!! DATA_TO_DIF sets up a divided difference table from raw data.
!
!  Discussion:
!
!    Space can be saved by using a single array for both the DIFTAB and
!    YTAB dummy parameters.  In that case, the divided difference table will
!    overwrite the Y data without interfering with the computation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real DIFTAB(NTAB), the divided difference coefficients
!    corresponding to the input (XTAB,YTAB).
!
!    Input, integer NTAB, the number of pairs of points
!    (XTAB(I),YTAB(I)) which are to be used as data.  The
!    number of entries to be used in DIFTAB, XTAB and YTAB.
!
!    Input, real XTAB(NTAB), the X values at which data was taken.
!    These values must be distinct.
!
!    Input, real YTAB(NTAB), the corresponding Y values.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 4 ) diftab(ntab)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical r4vec_distinct
  real ( kind = 4 ) xtab(ntab)
  real ( kind = 4 ) ytab(ntab)

  if ( .not. r4vec_distinct ( ntab, xtab ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_TO_DIF - Fatal error!'
    write ( *, '(a)' ) '  Two entries of XTAB are equal!'
    return
  end if
!
!  Copy the data values into DIFTAB.
!
  diftab(1:ntab) = ytab(1:ntab)
!
!  Compute the divided differences.
!
  do i = 2, ntab
    do j = ntab, i, -1

      diftab(j) = ( diftab(j) - diftab(j-1) ) / ( xtab(j) - xtab(j+1-i) )

    end do
  end do

  return
end
subroutine dif_to_poly ( diftab, ntab, xtab )

!*****************************************************************************80
!
!! DIF_TO_POLY converts a divided difference polynomial to standard form.
!
!  Discussion:
!
!    The vector DIFTAB, containing the divided difference polynomial
!    coefficients is overwritten with the standard form polynomial
!    coefficients, but the abscissa vector XTAB is unchanged.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real DIFTAB(NTAB).
!
!    On input, DIFTAB contains the coefficients of the divided difference
!    polynomials, corresponding to the XTAB array.
!
!    On output, DIFTAB contains the standard form polyomial coefficients.
!    DIFTAB(1) is the constant term, and DIFTAB(NTAB) is the coefficient
!    of X**(NTAB-1).
!
!    Input, integer NTAB, the number of coefficients, and abscissas.
!
!    Input, real XTAB(NTAB), the X values used in the divided difference
!    representation of the polynomial.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 4 ) diftab(ntab)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 4 ) xtab(ntab)
!
!  Recompute the divided difference coefficients.
!
  do j = 1, ntab-1
    do i = 1, ntab-j
      diftab(ntab-i) = diftab(ntab-i) - xtab(ntab-i-j+1) * diftab(ntab-i+1)
    end do
  end do

  return
end
subroutine dif_val ( diftab, ntab, xtab, xval, yval )

!*****************************************************************************80
!
!! DIF_VAL evaluates a divided difference polynomial at a point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real DIFTAB(NTAB), the divided difference polynomial coefficients.
!
!    Input, integer NTAB, the number of divided difference
!    coefficients in DIFTAB, and the number of points XTAB.
!
!    Input, real XTAB(NTAB), the X values upon which the
!    divided difference polynomial is based.
!
!    Input, real XVAL, a value of X at which the polynomial
!    is to be evaluated.
!
!    Output, real YVAL, the value of the polynomial at XVAL.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 4 ) diftab(ntab)
  integer ( kind = 4 ) i
  real ( kind = 4 ) xtab(ntab)
  real ( kind = 4 ) xval
  real ( kind = 4 ) yval

  yval = diftab(ntab)
  do i = 1, ntab-1
    yval = diftab(ntab-i) + ( xval - xtab(ntab-i) ) * yval
  end do

  return
end
subroutine poly_print ( n, poly_cof )

!*****************************************************************************80
!
!! POLY_PRINT prints out a polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the degree of the polynomial.
!
!    Input, real POLY_COF(0:N), the polynomial coefficients.
!    POLY_COF(0) is the constant term and
!    POLY_COF(N) is the coefficient of X**N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 4 ) mag
  integer ( kind = 4 ) n2
  character plus_minus
  real ( kind = 4 ) poly_cof(0:n)
!
!  Search for the first nonzero coefficient.  If none, we'll take
!  the constant term.
!
  do i = n, 1, -1

    n2 = i
    if ( poly_cof(i) /= 0 ) then
      go to 10
    end if

  end do

  n2 = 0

10    continue

  write ( *, '(a)' ) ' '

  if ( poly_cof(n2) < 0.0 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( poly_cof(n2) )

  if ( n2 >= 2 ) then
    write ( *, '( '' p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2
  else if ( n2 == 1 ) then
    write ( *, '( '' p(x) = '', a1, g14.6, '' * x'' )' ) plus_minus, mag
  else if ( n2 == 0 ) then
    write ( *, '( '' p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n2-1, 0, -1

    if ( poly_cof(i) < 0.0 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( poly_cof(i) )

    if ( mag /= 0.0 ) then

      if ( i >= 2 ) then
        write ( *, ' ( ''        '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''        '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''        '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
function r4vec_distinct ( n, a )

!*****************************************************************************80
!
!! R4VEC_DISTINCT is true if the entries in an R4VEC are distinct.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real A(N), the vector to be checked.
!
!    Output, logical R4VEC_DISTINCT is .TRUE. if all N elements of A
!    are distinct.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical r4vec_distinct

  r4vec_distinct = .false.

  do i = 2, n
    do j = 1, i - 1
      if ( a(i) == a(j) ) then
        return
      end if
    end do
  end do

  r4vec_distinct = .true.

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
subroutine s_is_r ( s, r, lval )

!*****************************************************************************80
!
!! S_IS_R is TRUE if a string represents a real number.
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
!    Input, character ( len = * ) S, the string to be checked.
!
!    Output, real R.  If the string represents a real number, then R
!    is the real number represented.  Otherwise R is 0.
!
!    Output, logical LVAL, is TRUE if the string represents a real number.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) lenc
  logical lval
  real ( kind = 4 ) r
  character ( len = * ) s

  lenc = len_trim ( s )

  call s_to_r ( s, r, ierror, lchar )

  if ( ierror == 0 .and. lchar >= lenc ) then
    lval = .true.
  else
    lval = .false.
    r = 0.0
  end if

  return
end
subroutine s_to_r ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R reads a real number from a string.
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
  real ( kind = 4 ) r
  real ( kind = 4 ) rbot
  real ( kind = 4 ) rexp
  real ( kind = 4 ) rtop
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
subroutine sequence_filler ( string )

!*****************************************************************************80
!
!! SEQUENCE_FILLER "fills in" a sequence with missing values.
!
!  Example:
!
!    Input:
!
!      STRING = '4 5 ? 14 24 ?'
!
!    Output (printed):
!
!      4.0
!      5.0
!      8.0
!     14.0
!     24.0
!     39.0
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
!    Input, character ( len = * ) STRING, contains a sequence with missing
!    values.  The missing values are represented by question marks.
!
  implicit none

  integer ( kind = 4 ), parameter :: MAXSEQ = 20

  real ( kind = 4 ) diftab(MAXSEQ)
  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iseq(MAXSEQ)
  logical lval
  integer ( kind = 4 ) ndata
  integer ( kind = 4 ) nseq
  character ( len = * ) string
  character ( len = 10 ) word
  real ( kind = 4 ) xdata(MAXSEQ)
  real ( kind = 4 ) xseq(MAXSEQ)
  real ( kind = 4 ) xval
  real ( kind = 4 ) yseq(MAXSEQ)
  real ( kind = 4 ) ydata(MAXSEQ)
  real ( kind = 4 ) yval
!
!  Read the tokens in the string.
!
  ndata = 0
  nseq = 0
  done = .true.

  do

    call word_next_read ( string, word, done )

    if ( done ) then
      exit
    end if

    if ( word == '?' ) then

      nseq = nseq + 1
      xseq(nseq) = nseq
      yseq(nseq) = 0.0E+00
      iseq(nseq) = 0

    else

      call s_is_r ( word, yval, lval )

      if ( .not. lval ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SEQUENCE_FILLER - Fatal error!'
        write ( *, '(a)' ) '  Could not recognize the input.'
        stop
      end if

      nseq = nseq + 1
      xseq(nseq) = nseq
      yseq(nseq) = yval
      iseq(nseq) = 1

      ndata = ndata + 1
      xdata(ndata) = nseq
      ydata(ndata) = yval

    end if

  end do
!
!  There should be at least one piece of data in the string.
!
  if ( ndata <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SEQUENCE_FILLER - Note'
    write ( *, '(a)' ) '  There wasn''t any numerical data in the'
    write ( *, '(a)' ) '  sequence, so no filing-in can be done.'
    return
  end if
!
!  If there were no question marks, then behave as though there was
!  one at the end.
!
  if ( ndata == nseq ) then
    nseq = nseq + 1
    xseq(nseq) = nseq
    yseq(nseq) = 0.0E+00
    iseq(nseq) = 0
  end if
!
!  Set up the divided difference table.
!
  call data_to_dif ( diftab, ndata, xdata, ydata )
!
!  Evaluate the divided difference polynomial at each unknown point.
!
  do i = 1, nseq
    if ( iseq(i) == 0 ) then
      xval = xseq(i)
      call dif_val ( diftab, ndata, xdata, xval, yval )
      yseq(i) = yval
    end if
  end do
!
!  Print the filled in sequence.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Sequence as completed by SEQUENCE_FILLER'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X         P(X)'
  write ( *, '(a)' ) ' '
  do i = 1, nseq
    write ( *, '(i6,g14.6)' ) i, yseq(i)
  end do
!
!  Convert the difference polynomial to standard form.
!
  call dif_to_poly ( diftab, ndata, xdata )
!
!  Print the polynomial.
!
  call poly_print ( ndata-1, diftab )

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
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
!
!    If DONE is FALSE, then WORD contains the "next" word read.
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!
!    On input with a fresh string, set DONE to TRUE.
!
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

10    continue
!
!  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  if ( ilo > lenc ) then
    word = ' '
    done = .true.
    next = lenc + 1
    return
  end if
!
!  If the current character is blank, skip to the next one.
!
  if ( s(ilo:ilo) == ' ' .or. s(ilo:ilo) == TAB ) then
    ilo = ilo + 1
    go to 10
  end if
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

20    continue

  if ( next > lenc ) then
    word = s(ilo:next-1)
    return
  end if

  if ( s(next:next) /= ' ' .and. &
       s(next:next) /= TAB .and. &
       s(next:next) /= '"' .and. &
       s(next:next) /= '(' .and. &
       s(next:next) /= ')' .and. &
       s(next:next) /= '{' .and. &
       s(next:next) /= '}' .and. &
       s(next:next) /= '[' .and. &
       s(next:next) /= ']' ) then

    next = next + 1
    go to 20

  end if

  word = s(ilo:next-1)

  return
end
