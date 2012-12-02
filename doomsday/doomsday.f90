subroutine doomsday_gregorian ( y, m, d, w )

!*****************************************************************************80
!
!! DOOMSDAY_GREGORIAN: weekday given any date in Gregorian calendar.
!
!  Discussion:
!
!    This procedure does not include any procedure to switch to the Julian
!    calendar for dates early enough that that calendar was used instead.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Conway,
!    Tomorrow is the Day After Doomsday,
!    Eureka,
!    Volume 36, October 1973, pages 28-31.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, the year, month and day of the date.
!    Note that the year must be positive.
!
!    Output, integer ( kind = 4  ) W, the weekday of the date.
!
  implicit none

  integer ( kind = 4 ), dimension ( 4 ) :: anchor = (/ 1, 6, 4, 3 /)
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) drd
  integer ( kind = 4 ) drdr
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ), dimension ( 12 ) :: mdoom = (/ &
    3, 28, 0, 4, 9, 6, 11, 8, 5, 10, 7, 12 /)
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) ydoom
  logical year_is_leap_gregorian
  integer ( kind = 4 ) yy
  integer ( kind = 4 ) yy12d
  integer ( kind = 4 ) yy12r
  integer ( kind = 4 ) yy12r4d
!
!  Refuse to handle Y <= 0.
!
  if ( y <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DOOMSDAY_GREGORIAN - Fatal error!'
    write ( *, '(a)' ) '  Y <= 0.'
    stop
  end if
!
!  Determine the century C.
!
  c = y / 100
!
!  Determine the last two digits of the year, YY
!
  yy = mod ( y, 100 )
!
!  Divide the last two digits of the year by 12.
!
  yy = mod ( y, 100 )
  yy12d = yy / 12
  yy12r = mod ( yy, 12 ) 
  yy12r4d = yy12r / 4
  drd = yy12d + yy12r + yy12r4d
  drdr = mod ( drd, 7 )
  ydoom = anchor( mod ( c-1, 4 ) + 1 ) + drdr
  ydoom = i4_wrap ( ydoom, 1, 7 )
!
!  If M = 1 or 2, and leap year, add 1.
!
  if ( ( m == 1 .or. m == 2 ) .and. year_is_leap_gregorian ( y ) ) then
    l = 1
  else
    l = 0
  end if

  w = ydoom + ( d -  mdoom(m) - l )
  w = i4_wrap ( w, 1, 7 )

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
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
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
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
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, a value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) value
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

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
function year_is_leap_gregorian ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_GREGORIAN returns TRUE if the Gregorian year was a leap year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_GREGORIAN, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_gregorian

  if ( y <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YEAR_IS_LEAP_GREGORIAN - Fatal error!'
    write ( *, '(a)' ) '  This function will not accept nonpositive years.'
    stop
  end if

  if ( mod ( y, 400 ) == 0 ) then
    year_is_leap_gregorian = .true.
  else if ( mod ( y, 100 ) == 0 ) then
    year_is_leap_gregorian = .false.
  else if ( mod ( y, 4 ) == 0 ) then
    year_is_leap_gregorian = .true.
  else
    year_is_leap_gregorian = .false.
  end if

  return
end
