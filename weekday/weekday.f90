subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
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
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the positive remainder when I is divided by J.
!
!  Discussion:
!
!    NREM = I4_MODP ( I, J )
!    NMULT = ( I - NREM ) / J
!
!    I = J * NMULT + NREM
!
!  Example:
!
!        I         J   NMULT  NREM    Factorization
!
!      107        50       2     7    107 =  2 *  50 + 7
!      107       -50      -2     7    107 = -2 * -50 + 7
!     -107        50      -3    43   -107 = -3 *  50 + 43
!     -107       -50       3    43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the positive remainder
!    when I is divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) i4_modp

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i6)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
subroutine i4_to_s_left ( intval, s )

!*****************************************************************************80
!
!! I4_TO_S_LEFT converts an integer to a left-justified string.
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  1
!        -1  -1
!         0  0
!      1952  1952
!    123456  123456
!   1234567  ******  <-- Not enough room!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be left-justified.  If there is not enough space,
!    the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * ) s

  if ( intval == 0 ) then
    s = '0'
    return
  end if

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Strip off the last digit of IVAL and stick it into the string.
!
  do while ( ival /= 0 )

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

  end do
!
!  Shift the string to the left.
!
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '

  return
end
subroutine i4_to_s_zero ( intval, s )

!*****************************************************************************80
!
!! I4_TO_S_ZERO converts an integer to a string, with zero padding.
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  000001
!        -1  -00001
!         0  000000
!      1952  001952
!    123456  123456
!   1234567  ******  <-- Not enough room!
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
!    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be right justified, and zero padded.
!    If there is not enough space, the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * ) s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  Working from right to left, strip off the digits of the integer
!  and place them into S(ILO:IHI).
!
  ipos = ihi

  do while ( ival /= 0 .or. ipos == ihi )

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

  end do
!
!  Fill the empties with zeroes.
!
  do i = ilo, ipos
    s(i:i) = '0'
  end do

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an integer to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the
!    integer value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    i4_wrap = jlo
  else
    i4_wrap = jlo + i4_modp ( ival - jlo, wide )
  end if

  return
end
subroutine jed_to_weekday ( jed, w, f )

!*****************************************************************************80
!
!! JED_TO_WEEKDAY computes the day of the week from a JED.
!
!  Discussion:
!
!    BC 4713/01/01 => JED = 0.0 was noon on a Monday.
!
!    jedmod = mod ( 0.0D+00, 7.0D+00 ) = 0.0D+00
!    j = mod ( nint ( 0 ), 7 ) = 0
!    f = ( 0.0D+00 + 0.5D+00 ) - real ( j ) = 0.5D+00
!    w = i4_wrap ( 0 + 2, 1, 7 ) = 2 = MONDAY
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) W, the day of the week of the date.
!    The days are numbered from Sunday through Saturday, 1 through 7.
!
!    Output, real ( kind = 8 ) F, the fractional part of the day.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  real ( kind = 8 ) jedmod
  integer ( kind = 4 ) w

  jedmod = mod ( jed, 7.0D+00 )

  j = mod ( nint ( jedmod ), 7 )

  f = ( jedmod + 0.5D+00 ) - real ( j, kind = 8 )

  w = i4_wrap ( j + 2, 1, 7 )

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
!    11 May 2000
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

  s3 = trim ( s1 ) // trim ( s2 )

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
subroutine weekday_to_name_common ( w, s )

!*****************************************************************************80
!
!! WEEKDAY_TO_NAME_COMMON returns the name of a Common weekday.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) W, the weekday index.
!
!    Output, character ( len = * ) S, the weekday name.
!
  implicit none

  character ( len = 9 ), parameter, dimension(7) :: name = (/ &
    'Sunday   ', 'Monday   ', 'Tuesday  ', 'Wednesday', &
    'Thursday ', 'Friday   ', 'Saturday ' /)
  character ( len = * ) s
  integer ( kind = 4 ) w
  integer ( kind = 4 ) w2
!
!  Make a local copy of the weekday number.
!
  w2 = w
!
!  Return the weekday name.
!
  s = name ( w2 )

  return
end
subroutine weekday_values ( n_data, y, m, d, w )

!*****************************************************************************80
!
!! WEEKDAY_VALUES returns the day of the week for various dates.
!
!  Discussion:
!
!    The CE or Common Era calendar is used, under the
!    hybrid Julian/Gregorian Calendar, with a transition from Julian
!    to Gregorian.  The day after 04 October 1582 was 15 October 1582.
!
!    The year before 1 AD or CE is 1 BC or BCE.  In this data set,
!    years BC/BCE are indicated by a negative year value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Reingold, Nachum Dershowitz,
!    Calendrical Calculations: The Millennium Edition,
!    Cambridge University Press, 2001,
!    ISBN: 0 521 77752 6
!    LC: CE12.R45.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) Y, M, D, the Common Era date.
!
!    Output, integer ( kind = 8 ) W, the day of the week.  Sunday = 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 34

  integer ( kind = 4 ) d
  integer ( kind = 4 ), save, dimension ( n_max ) :: d_vec = (/ &
    30, &
     8, &
    26, &
     3, &
     7, &
    18, &
     7, &
    19, &
    14, &
    18, &
    16, &
     3, &
    26, &
    20, &
     4, &
    25, &
    31, &
     9, &
    24, &
    10, &
    30, &
    24, &
    19, &
     2, &
    27, &
    19, &
    25, &
    29, &
    19, &
     7, &
    17, &
    25, &
    10, &
    18 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
     7, &
    12, &
     9, &
    10, &
     1, &
     5, &
    11, &
     4, &
    10, &
     5, &
     3, &
     3, &
     3, &
     4, &
     6, &
     1, &
     3, &
     9, &
     2, &
     6, &
     6, &
     7, &
     6, &
     8, &
     3, &
     4, &
     8, &
     9, &
     4, &
    10, &
     3, &
     2, &
    11, &
     7 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) w
  integer ( kind = 4 ), save, dimension ( n_max ) :: w_vec = (/ &
    1, &
    4, &
    4, &
    1, &
    4, &
    2, &
    7, &
    1, &
    7, &
    1, &
    6, &
    7, &
    6, &
    1, &
    1, &
    4, &
    7, &
    7, &
    7, &
    4, &
    1, &
    6, &
    1, &
    2, &
    4, &
    1, &
    1, &
    2, &
    2, &
    5, &
    3, &
    1, &
    4, &
    1 /)
  integer ( kind = 4 ) y
  integer ( kind = 4 ), save, dimension ( n_max ) :: y_vec = (/ &
    - 587, &
    - 169, &
       70, &
      135, &
      470, &
      576, &
      694, &
     1013, &
     1066, &
     1096, &
     1190, &
     1240, &
     1288, &
     1298, &
     1391, &
     1436, &
     1492, &
     1553, &
     1560, &
     1648, &
     1680, &
     1716, &
     1768, &
     1819, &
     1839, &
     1903, &
     1929, &
     1941, &
     1943, &
     1943, &
     1992, &
     1996, &
     2038, &
     2094 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    y = 0
    m = 0
    d = 0
    w = 0
  else
    y = y_vec(n_data)
    m = m_vec(n_data)
    d = d_vec(n_data)
    w = w_vec(n_data)
  end if

  return
end
subroutine y_common_to_astronomical ( y, y2 )

!*****************************************************************************80
!
!! Y_COMMON_TO_ASTRONOMICAL converts a Common year to an Astronomical year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the Common year.
!
!    Output, integer ( kind = 4 ) Y2, the Astronomical year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2

  if ( y < 0 ) then
    y2 = y + 1
  else if ( y == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Y_COMMON_TO_ASTRONOMICAL - Fatal error!'
    write ( *, '(a)' ) '  COMMON calendar does not have a year 0.'
    stop
  else
    y2 = y
  end if

  return
end
subroutine ymd_to_s_common ( y, m, d, s )

!*****************************************************************************80
!
!! YMD_TO_S_COMMON writes a Common YMD date into a string.
!
!  Format:
!
!    CE YYYY/MM/DD
!    BCE YYYY/MM/DD
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, the YMD date.
!
!    Output, character ( len = * ) S, a representation of the date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  character ( len = 20 ) s1
  character ( len = 2 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y
  m2 = m
  d2 = d

  if ( 0 <= y2 ) then
    s1 = 'CE '
    call i4_to_s_left ( y2, s1(4:) )
  else
    s1 = 'BCE '
    call i4_to_s_left (  - y2, s1(5:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( m2, s2 )

  call s_cat ( s1, s2, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( d2, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine ymd_to_weekday_common ( y, m, d, w )

!*****************************************************************************80
!
!! YMD_TO_WEEKDAY_COMMON returns the weekday of a Common YMD date.
!
!  Discussion:
!
!    The "common" calendar is meant to be the calendar which is Julian up to
!    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, the YMD date.
!
!    Output, integer ( kind = 4 ) W, is the week day number of the date, with
!    1 for Sunday, through 7 for Saturday.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y

  f = 0.5D+00

  call ymdf_to_jed_common ( y, m, d, f, jed )

  call jed_to_weekday ( jed, w, f2 )

  return
end
subroutine ymdf_compare ( y1, m1, d1, f1, y2, m2, d2, f2, cmp )

!*****************************************************************************80
!
!! YMDF_COMPARE compares two YMDF dates.
!
!  Discussion:
!
!    The comparison should work for a pair of dates in any calendar.
!
!    No check is made that the dates are actually legitimate.  It is
!    assumed that the calling routine has already ensured that the
!    dates are properly "normalized".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1, the
!    first YMDF date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2, the
!    second YMDF date.
!
!    Output, character CMP:
!    '<' if date 1 precedes date 2;
!    '=' if date 1 equals date 2;
!    '>' if date 1 follows date 2;
!
  implicit none

  character cmp
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  cmp = '?'
!
!  Compare years...
!
  if ( y1 < y2 ) then
    cmp = '<'
  else if ( y1 > y2 ) then
    cmp = '>'
  else
!
!  ...if necessary, compare months in equal years...
!
    if ( m1 < m2 ) then
      cmp = '<'
    else if ( m1 > m2 ) then
      cmp = '>'
    else
!
!  ...if necessary, compare days in equal months...
!
      if ( d1 < d2 ) then
        cmp = '<'
      else if ( d1 > d2 ) then
        cmp = '>'
      else
!
!  ...if necessary, compare fractional parts.
!
        if ( f1 < f2 ) then
          cmp = '<'
        else if ( f1 > f2 ) then
          cmp = '>'
        else
          cmp = '='
        end if

      end if

    end if

  end if

  return
end
subroutine ymdf_to_jed_common ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_COMMON converts a Common YMDF date to a JED.
!
!  Discussion:
!
!    The "common" calendar is meant to be the calendar which is Julian up to
!    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
!
!    The Julian Ephemeris Date is essentially a count of the number
!    of days that have elapsed since noon, 1 January 4713 BC, at
!    Greenwich, England.  Strictly speaking, the Julian Ephemeris Date
!    is counted from noon, and thus day "0" began at noon on 1 January 4713 BC,
!    and ended at noon on 2 January 4713 BC.
!
!    The Julian Ephemeris Date was devised by Joseph Scaliger in 1583.
!
!    The Julian Ephemeris Date has been adopted by astronomers as
!    a convenient reference for dates.
!
!  Example:
!
!       Y   M     D         JED
!    --------------     -------
!    BC 4713 Jan  1           0
!    AD 1968 May 23     2440000
!    AD 1984 Dec 31     2446065
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  character cmp
  integer ( kind = 4 ) d
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the month and year.
!
  y1 = y
  m1 = m
  d1 = d
  f1 = f

  y2 = 1582
  m2 = 10
  d2 = 4+1
  f2 = 0.0D+00

  call ymdf_compare ( y1, m1, d1, f1, y2, m2, d2, f2, cmp )

  if ( cmp == '<' ) then
    call ymdf_to_jed_julian ( y1, m1, d1, f1, jed )
    return
  end if
!
!  Use the Gregorian calendar for dates strictly after 1752/9/13.
!
  y2 = 1582
  m2 = 10
  d2 = 15-1
  f2 = 0.0D+00

  call ymdf_compare ( y1, m1, d1, f1, y2, m2, d2, f2, cmp )

  if ( cmp == '>' ) then
    call ymdf_to_jed_gregorian ( y1, m1, d1, f1, jed )
    return
  end if

  jed = -1.0D+00
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'YMDF_TO_JED_COMMON - Error!'
  write ( *, '(a)' ) '  Illegal date!'

  return
end
subroutine ymdf_to_jed_gregorian ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_GREGORIAN converts a Gregorian YMDF date to a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Algorithm E,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, pages 323-324.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, real ( kind = 8 ) JED, the corresponding JED.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  integer ( kind = 4 ) g
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) jed
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y_prime
!
!  Account for the missing year 0 by moving negative years up one.
!
  call y_common_to_astronomical ( y, y2 )
!
!  Convert the calendar date to a computational date.
!
  y_prime = y2 + 4716 - ( 14 - m ) / 12
  m_prime = mod ( m + 9, 12 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  j1 = ( 1461 * y_prime ) / 4

  j2 = ( 153 * m_prime + 2 ) / 5

  g = ( 3 * ( ( y_prime + 184 ) / 100 ) / 4 ) - 38

  jed = real ( j1 + j2 + d_prime - 1401 - g, kind = 8 ) - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_julian ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_JULIAN converts a Julian YMDF date to a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Algorithm E,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, pages 323-324.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) jed
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y_prime
!
!  Account for the missing year 0 by moving negative years up one.
!
  call y_common_to_astronomical ( y, y2 )
!
!  Convert the calendar date to a computational date.
!
  y_prime = y2 + 4716 - ( 14 - m ) / 12
  m_prime = mod ( m + 9, 12 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  j1 = ( 1461 * y_prime ) / 4

  j2 = ( 153 * m_prime + 2 ) / 5

  jed = real ( j1 + j2 + d_prime - 1401, kind = 8 ) - 0.5D+00
  jed = jed + f

  return
end
