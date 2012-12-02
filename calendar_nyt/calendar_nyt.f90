subroutine day_borrow_common ( y, m, d )

!*****************************************************************************80
!
!! DAY_BORROW_COMMON borrows days from months in a Common date.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, a year, month, and day
!    representing a date.  On input, D might be negative.  On output,
!    M should have decreased by one month, and D gone up by the
!    number of days in the month we "cashed in".  Y may be affected
!    if the input value of M was 1.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) y

  do while ( d <= 0 )

    m = m - 1

    call month_borrow_common ( y, m )

    days = month_length_common ( y, m )

    d = d + days

  end do

  return
end
subroutine day_borrow_gregorian ( y, m, d )

!*****************************************************************************80
!
!! DAY_BORROW_GREGORIAN borrows days from months in a Gregorian date.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, a year, month, and day
!    representing a date.  On input, D might be negative.  On output,
!    M should have decreased by one month, and D gone up by the
!    number of days in the month we "cashed in".  Y may be affected
!    if the input value of M was 1.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_gregorian
  integer ( kind = 4 ) y

  do while ( d <= 0 )

    m = m - 1

    call month_borrow_gregorian ( y, m )

    days = month_length_gregorian ( y, m )

    d = d + days

  end do

  return
end
subroutine day_borrow_julian ( y, m, d )

!*****************************************************************************80
!
!! DAY_BORROW_JULIAN borrows days from months in a Julian date.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, a year, month, and day
!    representing a date.  On input, D might be negative.  On output,
!    M should have decreased by one month, and D gone up by the
!    number of days in the month we "cashed in".  Y may be affected
!    if the input value of M was 1.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_julian
  integer ( kind = 4 ) y

  do while ( d <= 0 )

    m = m - 1

    call month_borrow_julian ( y, m )

    days = month_length_julian ( y, m )

    d = d + days

  end do

  return
end
subroutine day_carry_common ( y, m, d )

!*****************************************************************************80
!
!! DAY_CARRY_COMMON carries days to months in a Common date.
!
!  Discussion:
!
!    While ( number of days in M ) < D:
!      decrease the day D by the number of days in the month M;
!      increase M by 1;
!      if necessary, adjust Y.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, the YMD date.
!    On output, D is between 1 and the number of days in M.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) y
!
!  If the date is in the transition month, deflate it,
!  so we can perform ordinary arithmetic.
!
  call deflate_common ( y, m, d )

  days = month_length_common ( y, m )

  do while ( days < d )

    d = d - days
    m = m + 1
    days = month_length_common ( y, m )
!
!  Make sure the month isn't too big.
!
    call month_carry_common ( y, m )

  end do
!
!  If the date is in the transition month, inflate it.
!
  call inflate_common ( y, m, d )

  return
end
subroutine day_carry_gregorian ( y, m, d )

!*****************************************************************************80
!
!! DAY_CARRY_GREGORIAN carries days to months in a Gregorian date.
!
!  Discussion:
!
!    While ( number of days in M ) < D:
!      decrease the day D by the number of days in the month M;
!      increase M by 1;
!      if necessary, adjust Y.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, the YMD date.
!    On output, D is between 1 and the number of days in M.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_gregorian
  integer ( kind = 4 ) y

  days = month_length_gregorian ( y, m )

  do while ( days < d )

    d = d - days
    m = m + 1
    days = month_length_gregorian ( y, m )
!
!  Make sure the month isn't too big.
!
    call month_carry_gregorian ( y, m )

  end do

  return
end
subroutine day_carry_julian ( y, m, d )

!*****************************************************************************80
!
!! DAY_CARRY_JULIAN carries days to months in a Julian date.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, the YMD date.
!    On output, D is between 1 and the number of days in M.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_julian
  integer ( kind = 4 ) y

  days = month_length_julian ( y, m )

  do while ( days < d )

    d = d - days
    m = m + 1
    days = month_length_julian ( y, m )
!
!  Make sure the month isn't too big.
!
    call month_carry_julian ( y, m )

  end do

  return
end
subroutine deflate_common ( y, m, d )

!*****************************************************************************80
!
!! DEFLATE_COMMON "deflates" dates in the Common Calendar transition month.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, the YMD date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  if ( y == 1582 ) then
    if ( m == 10 ) then
      if ( 15 <= d ) then
        d = d - 10
      end if
    end if
  end if

  return
end
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
subroutine epoch_to_jed_nyt ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_NYT returns the epoch of the NYT calendar as a JED.
!
!  Discussion:
!
!    The "epoch" of the NYT calendar is the mythical date when issue "0"
!    would have been printed, namely, a tad past midnight, 17 September 1851.
!
!    Volume #1, Issue #1 was printed on 18 September 1851.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date of the epoch.
!
  implicit none

  real ( kind = 8 ) jed

  jed = 2397382.5D+00
!
!  The following value is effectively the JED we are using for an
!  epoch set to the nominal issue number 50,000.
!
! jed = 2449790.5D+00

  return
end
subroutine frac_borrow_common ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_BORROW_COMMON borrows fractions from days in a Common YMDF date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F,
!    a YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  do while ( f < 0.0D+00 )

    f = f + 1.0D+00

    d = d - 1

  end do

  call day_borrow_common ( y, m, d )

  return
end
subroutine frac_borrow_julian ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_BORROW_JULIAN borrows fractions from days in a Julian YMDF date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F,
!    a YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  do while ( f < 0.0D+00 )

    f = f + 1.0D+00

    d = d - 1

  end do

  call day_borrow_julian ( y, m, d )

  return
end
subroutine frac_carry_common ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_CARRY_COMMON carries fractions to days in a Common YMDF date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F, the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  real ( kind = 8 ) f
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  if ( f < 1.0D+00 ) then
    return
  end if

  days = int ( f )

  f = f - real ( days, kind = 8 )
  d = d + days

  call day_carry_common ( y, m, d )

  return
end
subroutine frac_carry_julian ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_CARRY_JULIAN carries fractions from days in a Julian YMDF date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F, a YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  real ( kind = 8 ) f
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  if ( f < 1.0D+00 ) then
    return
  end if

  days = int ( f )

  f = f - real ( days, kind = 8 )
  d = d + days

  call day_carry_julian ( y, m, d )

  return
end
subroutine frac_to_s ( f, s )

!*****************************************************************************80
!
!! FRAC_TO_S writes a positive fraction into a left justified character string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) F, the number to be written into the string.
!    F should be between 0.0 and 1.0.
!
!    Output, character ( len = * ) S, a representation of F.
!
  implicit none

  real ( kind = 8 ) f
  character ( len = * ) s
  character ( len = 14 ) s2

  if ( f < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FRAC_TO_S - Fatal error!'
    write ( *, '(a)' ) '  The input fraction was negative:'
    write ( *, '(g14.6)' ) f
    stop
  else if ( 1.0D+00 <= f ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FRAC_TO_S - Fatal error!'
    write ( *, '(a)' ) '  The input fraction was 1 or more:'
    write ( *, '(g14.6)' ) f
    stop
  end if

  write ( s2, '(f11.10)' ) f

  s = s2

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
subroutine inflate_common ( y, m, d )

!*****************************************************************************80
!
!! INFLATE_COMMON "inflates" dates in the Common Calendar transition month.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, the YMD date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  if ( y == 1582 ) then
    if ( m == 10 ) then
      if ( 5 <= d ) then
        d = d + 10
      end if
    end if
  end if

  return
end
subroutine jed_to_nyt ( jed, volume, issue )

!*****************************************************************************80
!
!! JED_TO_NYT converts a JED to an NYT date.
!
!  Discussion:
!
!    The New York Times began publication with Volume 1, Issue 1 on
!    Thursday, 18 September 1851.
!
!    The Volume number is incremented annually, on 18 September.
!
!    It seemed an initriguing idea, then, to devise a formula that would
!    produce the New York Times issue number for a given date, or that
!    could start with the issue number and return the date on which that
!    issue appeared.
!
!    In a simple world, this would have been essentially a translation
!    of the JED, that is, the first approximation would be
!
!      NYT(today) - NYT(18 September 1851) = JED(today) - JED(18 September 1851)
!
!    so
!
!      NYT(today) = NYT(18 September 1851) + JED(today) - JED(18 September 1851)
!
!    and we're done.
!
!    However, the first problem involved Sunday issues.  No Sunday paper was
!    printed at all, until 21 April 1861.  Moreover, that paper was given
!    the issue number 2990, which was the same issue number as the Saturday
!    paper.  This persisted until some time around April 1905, when Sunday
!    papers were assigned their own issue number.   Once this was done, the
!    New York Times issue number began to "track" the Julian Ephemeris Date
!    in a simple way.
!
!    The second obvious problem occurred because there was an 88 day strike
!    in 1978.  The issue for August 9 was 44027, and the issue for November 6
!    was 44028 (I THINK, I AM NOT COMPLETELY SURE HERE).  During other strikes,
!    the New York Times had increased the issue number each day, even if no
!    paper was printed.  This was the first time that a strike caused the
!    issue number sequence to halt.
!
!    The third problem was more subtle.  An article printed on 14 March 1995
!    heralded the printing of issue 50,000 of the New York Times.  It also
!    mentioned issues and corresponding dates for several points in the past,
!    explained the 88 day strike lacuna, and the fact that there were no
!    Sunday papers at all until 21 April 1861.  This information seemed enough
!    to define a new formula that would work for the present era, at least,
!    after Sunday papers were printed and given their own issue number.
!    We could do this by basing the formula on the JED for issue 50,000, which
!    turned out to have the value 2449790.5.
!
!    For days on or after 6 November 1978,
!
!      NYT(today) = NYT(14 March 1995) + JED(today) - JED(14 March 1995)
!
!    For days on or before 9 August 1978,
!
!      NYT(today) = NYT(14 March 1995) + JED(today) - JED(14 March 1995) + 88
!
!    I set up this formula, and it worked pretty well for my list of known
!    dates and issue numbers between 1909 and 1995.
!
!    Then I pulled out the New York Times that I had bought that day
!    (22 November 2007), and tried out the formula.  To my dismay, the value
!    returned by the formula was exactly 500 higher than the value printed
!    on my paper.  This was very disturbing!
!
!    Going online, I tried to find more examples of issues and dates in the
!    interval between 14 March 1995 and 22 November 2007.  This was harder
!    than you might think.  Almost no one refers to the issue number.  Even
!    the article indexes for the New York Times, whether printed or online,
!    refer only to the date.  I ended up having to scan for images of the
!    front page.  There were surprisingly many, but most were of such poor
!    quality that the issue number could not be read.  Patience rewarded
!    me, though, with data for 1997, then for 2005, then for 2003, then
!    2002.  Gradually, I began to jokingly assume that the dreaded Year 2000
!    catastrophe had somehow hit the New York Times!
!
!    Imagine my shock when a colleague whom I had dragged into the search
!    with me discovered that this was true in a way.  On the front page of
!    the New York Times for 1 January 2000 was the statement that a mistake
!    in issue numbering, made in 1898 and never noticed until recently,
!    was being corrected.  The issue numbers had been accidentally "inflated"
!    by 500 back then, so they were now being "corrected" by dropping 500.
!
!    The ghastly details were:
!
!      Date              Issue
!      ----------------  ------
!      06 February 1898  14,499
!      07 February 1898  15,000
!      31 December 1999  51,753
!      01 January  2000  51,254
!
!    With this information, it becomes possible to adjust the formula to
!    be correct for current issues, and back over the "hiccup" in 1898.
!    The formula, however, obviously becomes more complicated, and, what's
!    worse, the issue number itself no longer can be used to deduce the
!    date, since there is now a sequence of 500 issue numbers that were used
!    twice.  Luckily, if we require the Volume number as well, we have
!    enough information to go back and forth.
!
!    The formula for the New York Times Volume Number is not as simple
!    as it could be.  The Volume started at 1 on 18 September 1851, and
!    increases by 1 each successive 18 September.  To determine the
!    volume number for a given date, you need to go "backwards" to the
!    most recent elapsed 18 September, note the year in which that occurred,
!    subtract 1851 from that, and add 1!
!
!      NYT_VOLUME = Year(Most-recently-elapsed-18-September) - 1851 + 1.
!
!    Now I have to work out the details of the formula to allow for the
!    two hiccups and the strike, and I should have a start on a usable pair
!    of formulas.
!
!    This excruciating (and unfinished) effort demonstrates, I hope, that
!    calendars are human creations, which aspire to mathematical regularity,
!    but which inevitably acquire the irregularities and ambiguities of all
!    human creations!
!
!    Surprisingly, computing the correct issue number from the date
!    is complicated.  Here are a few of the misadventures:
!
!      Fri,  2 Jan 1852, no issue.
!            6 Jul 1852, no issue
!      Sat,  2 Jul 1853, no issue, would have been 559.
!      Mon,  4 Jul 1853, INCORRECT issue number 560 (559 not used).
!      Tue,  5 Jul 1853, correct issue number 560.
!            6 Jul 1854, issue, but same issue number as 5 Jul 1854.
!      Thu,  5 Jul 1855, issue, but same issue number as 4 Jul 1855 (#1184)
!      Tue, 25 Sep 1855, issue jumps by 2, from 1253 to 1255!
!      Sat, 29 Sep 1856, issue, but same issue number as 28 Sep 1855 (#1258).
!      Fri,  4 Jan 1856, issue, but same issue number as 3 Jan 1856, (#1340).
!      Mon,  7 Jul 1856, issue, but same issue number as 5 Jul 1856, (#1497).
!      Sat,  3 Jan 1857, issue, but same issue number as 2 Jan 1857, (#1651).
!      Sat,  2 Jan 1858, issue, but same issue number as 1 Jan 1858, (#1962).
!      Tue,  6 Jul 1858, issue, but same issue number as 5 Jul 1858, (#2119).
!      Tue,  5 Jul 1859, no issue.
!      Tue,  3 Jan 1860, no issue.
!      Thu,  5 Jul 1860, no issue.
!      Wed,  2 Jan 1861, no issue
!      Sun, 21 Apr 1861, first Sunday issue.  First two Sundays get distinct
!                        issue numbers.  Thereafter, a "correction" is made.
!      Fri,  5 Jul 1861, no issue.
!      Thu,  2 Jan 1862, no issue.
!      Sat,  5 Jul 1862, no issue.
!      Fri,  2 Jan 1863, no issue.
!      Sat,  2 Jan 1864, no issue.
!      Tue,  5 Jul 1864, no issue.
!      Wed,  5 Jul 1865, no issue.
!      Tue,  2 Jan 1866, no issue.
!      Wed,  2 Jan 1867, no issue.
!      Sat,  5 Feb 1898, issue 14499.
!      Mon,  7 Feb 1898, issue 15000 (incremented by 501 instead of by 1)
!      Sun, 23 Apr 1905, Sunday paper gets distinct issue number.
!      Wed,  9 Aug 1978, last prestrike issue.  Issue numbers halt.
!      Mon,  6 Nov 1978, first poststrike issue, issue numbers resume.
!      Sat,  1 Jan 2000, issue numbers "corrected" downwards by 500.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Anonymous,
!    A Correction Welcome to 51,254,
!    The New York Times,
!    01 January 2000, Volume 149, Issue 51254.
!
!    James Barron,
!    What's in a Number? 143 Years of News,
!    The New York Times,
!    14 March 1995, Volume 144, Issue 50000.
!
!    The New York Times,
!    Page One, 1896-1996, A Special Commemorative Edition Celebrating the
!    100th Anniversary of the Purchase of the New York Times by Adolph S Ochs,
!    Galahad Books, 1996,
!    ISBN: 0-88365-961-1,
!    LC: D411.P25.
!
!    The New York Times,
!    The Complete First Pages, 1851-2008,
!    Black Dog & Leventhal Publishers, 2008,
!    ISBN13: 978-1-57912-749-7,
!    LC: D351.N53.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) VOLUME, ISSUE, the New York Times
!   volume and issue.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) days
  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  integer ( kind = 4 ) issue
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_02_01_1852
  real ( kind = 8 ) jed_07_06_1852
  real ( kind = 8 ) jed_07_02_1853
  real ( kind = 8 ) jed_07_06_1854
  real ( kind = 8 ) jed_07_05_1855
  real ( kind = 8 ) jed_25_09_1855
  real ( kind = 8 ) jed_29_09_1855
  real ( kind = 8 ) jed_04_01_1856
  real ( kind = 8 ) jed_07_07_1856
  real ( kind = 8 ) jed_03_01_1857
  real ( kind = 8 ) jed_02_01_1858
  real ( kind = 8 ) jed_06_07_1858
  real ( kind = 8 ) jed_05_07_1859
  real ( kind = 8 ) jed_03_01_1860
  real ( kind = 8 ) jed_05_07_1860
  real ( kind = 8 ) jed_02_01_1861
  real ( kind = 8 ) jed_04_21_1861
  real ( kind = 8 ) jed_04_28_1861
  real ( kind = 8 ) jed_05_05_1861
  real ( kind = 8 ) jed_05_07_1861
  real ( kind = 8 ) jed_02_01_1862
  real ( kind = 8 ) jed_05_07_1862
  real ( kind = 8 ) jed_02_01_1863
  real ( kind = 8 ) jed_28_09_1863
  real ( kind = 8 ) jed_30_09_1863
  real ( kind = 8 ) jed_02_01_1864
  real ( kind = 8 ) jed_05_07_1864
  real ( kind = 8 ) jed_03_01_1865
  real ( kind = 8 ) jed_05_07_1865
  real ( kind = 8 ) jed_02_01_1866
  real ( kind = 8 ) jed_02_01_1867
  real ( kind = 8 ) jed_07_02_1898
  real ( kind = 8 ) jed_22_04_1905
  real ( kind = 8 ) jed_10_08_1978
  real ( kind = 8 ) jed_05_11_1978
  real ( kind = 8 ) jed_01_01_2000
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) sundays
  integer ( kind = 4 ) volume
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!    The "epoch" of the NYT calendar is the mythical date when issue "0"
!    would have been printed, namely, a tad past midnight, 17 September 1851.
!
!    Volume #1, Issue #1 was printed on 18 September 1851.
!
  y2 = 1851
  m2 = 9
  d2 = 17
  f2 = 0.0D+00
  call ymdf_to_jed_common ( y2, m2, d2, f2, jed_epoch )
!
!  We start out by computing the number of elapsed days, which is
!  our initial estimate of the issue number.
!
  issue = jed - jed_epoch
!
!  If the user has given a JED before the epoch, return now.
!
  if ( issue <= 0 ) then
    volume = -1
    return
  end if
!
!  For dates on or after issue #1, the volume computation is easy.
!
  call jed_to_ymdf_common ( jed, y, m, d, f )

  volume = y - 1851 + 1

  if ( ( m == 9 .and. d < 18 ) .or. m < 9 ) then
    volume = volume - 1
  end if

  f = 0.0D+00
!
!  CORRECTION #1
!  Deal with nonissue on Friday, 2 January 1852
!
  call ymdf_to_jed_common ( 1852, 1, 2, f, jed_02_01_1852 )

  if ( jed_02_01_1852 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #2
!  Deal with nonissue on Tuesday, 6 July 1852
!
  call ymdf_to_jed_common ( 1852, 7, 6, f, jed_07_06_1852 )

  if ( jed_07_06_1852 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #3
!  Deal with nonissue on Saturday, 2 July 1853
!
  call ymdf_to_jed_common ( 1853, 7, 2, f, jed_07_02_1853 )

  if ( jed_07_02_1853 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #4
!  Deal with use of single issue number 873 for both
!  Wednesday, 5 July 1854 and Thursday, 6 July 1854.
!
  call ymdf_to_jed_common ( 1854, 7, 6, f, jed_07_06_1854 )

  if ( jed_07_06_1854 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #5
!  Deal with use of single issue number 1184 for both
!  Wednesday, 4 July 1855 and Thursday, 5 July 1855.
!
  call ymdf_to_jed_common ( 1855, 7, 5, f, jed_07_05_1855 )

  if ( jed_07_05_1855 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #6
!  They skipped one!
!  Issue 1253 = 24 September 1855
!  Issue 1255 = 25 September 1855
!
  call ymdf_to_jed_common ( 1855, 9, 25, f, jed_25_09_1855 )

  if ( jed_25_09_1855 <= jed ) then
    issue = issue + 1
  end if
!
!  CORRECTION #7
!  They "fixedÓ it.
!  Issue 1258 = 28 September 1855
!  Issue 1258 = 29 September 1855
!
  call ymdf_to_jed_common ( 1855, 9, 29, f, jed_29_09_1855 )

  if ( jed_29_09_1855 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #8
!  Deal with use of single issue number 1340 for both
!  Thursday, 3 January 1856 and Friday, 4 January 1856.
!
  call ymdf_to_jed_common ( 1856, 1, 4, f, jed_04_01_1856 )

  if ( jed_04_01_1856 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #8
!  Deal with use of single issue number 1497 for both
!  Saturday, 5 July 1856 and Monday, 7 July 1856.
!
  call ymdf_to_jed_common ( 1856, 7, 7, f, jed_07_07_1856 )

  if ( jed_07_07_1856 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #9
!  Deal with use of single issue number 1651 for both
!  Friday, 2 January 1857 and Saturday, 3 January 1857.
!
  call ymdf_to_jed_common ( 1857, 1, 3, f, jed_03_01_1857 )

  if ( jed_03_01_1857 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #10
!  Deal with use of single issue number 1962 for both
!  Friday, 1 January 1858 and Saturday, 2 January 1858.
!
  call ymdf_to_jed_common ( 1858, 1, 2, f, jed_02_01_1858 )

  if ( jed_02_01_1858 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #11
!  Deal with use of single issue number 2119 for both
!  Monday, 5 July 1858 and Tuesday, 6 July 1858.
!
  call ymdf_to_jed_common ( 1858, 7, 6, f, jed_06_07_1858 )

  if ( jed_06_07_1858 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #12
!  Deal with nonissue on Tuesday, 5 July 1859:
!
  call ymdf_to_jed_common ( 1859, 7, 5, f, jed_05_07_1859 )

  if ( jed_05_07_1859 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #13
!  Deal with nonissue on Tuesday, 3 January 1860:
!
  call ymdf_to_jed_common ( 1860, 1, 3, f, jed_03_01_1860 )

  if ( jed_03_01_1860 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #14
!  Deal with nonissue on Thursday, 5 July 1860:
!
  call ymdf_to_jed_common ( 1860, 7, 5, f, jed_05_07_1860 )

  if ( jed_05_07_1860 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #15
!  Deal with nonissue on Wednesday, 2 January 1861:
!
  call ymdf_to_jed_common ( 1861, 1, 2, f, jed_02_01_1861 )

  if ( jed_02_01_1861 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #16
!  Sunday issue got its own issue number, 21 April 1861.
!
  call ymdf_to_jed_common ( 1861, 4, 21, f, jed_04_21_1861 )

  if ( jed_04_21_1861 <= jed ) then
    issue = issue + 1
  end if
!
!  CORRECTION #17
!  Sunday issue got its own issue number, 28 April 1861.
!
  call ymdf_to_jed_common ( 1861, 4, 28, f, jed_04_28_1861 )

  if ( jed_04_28_1861 <= jed ) then
    issue = issue + 1
  end if
!
!  CORRECTION #18
!  Two Sunday issues retroactively "corrected" back down, 5 May 1861.
!
  call ymdf_to_jed_common ( 1861, 5, 5, f, jed_05_05_1861 )

  if ( jed_05_05_1861 <= jed ) then
    issue = issue - 2
  end if
!
!  CORRECTION #19
!  Deal with nonissue on Friday, 5 July 1861:
!
  call ymdf_to_jed_common ( 1861, 7, 5, f, jed_05_07_1861 )

  if ( jed_05_07_1861 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #20
!  Deal with nonissue on Thursday, 2 January 1862:
!
  call ymdf_to_jed_common ( 1862, 1, 2, f, jed_02_01_1862 )

  if ( jed_02_01_1862 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #21
!  Deal with nonissue on Saturday, 5 July 1862:
!
  call ymdf_to_jed_common ( 1862, 7, 5, f, jed_05_07_1862 )

  if ( jed_05_07_1862 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #22
!  Deal with nonissue on Friday, 2 January 1863:
!
  call ymdf_to_jed_common ( 1863, 1, 2, f, jed_02_01_1863 )

  if ( jed_02_01_1863 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #23
!  Deal with failure of issue increment on Monday, 28 September 1863:
!
  call ymdf_to_jed_common ( 1863, 9, 28, f, jed_28_09_1863 )

  if ( jed_28_09_1863 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #23
!  Deal with double issue increment on Wednesday, 30 September 1863:
!
  call ymdf_to_jed_common ( 1863, 9, 30, f, jed_30_09_1863 )

  if ( jed_30_09_1863 <= jed ) then
    issue = issue + 1
  end if
!
!  CORRECTION #24
!  Deal with nonissue on Saturday, 2 January 1864:
!
  call ymdf_to_jed_common ( 1864, 1, 2, f, jed_02_01_1864 )

  if ( jed_02_01_1864 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #25
!  Deal with nonissue on Tuesday, 5 July 1864:
!
  call ymdf_to_jed_common ( 1864, 7, 5, f, jed_05_07_1864 )

  if ( jed_05_07_1864 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #26
!  Deal with nonissue on Monday, 3 January 1865:
!
  call ymdf_to_jed_common ( 1865, 1, 3, f, jed_03_01_1865 )

  if ( jed_03_01_1865 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #27
!  Deal with nonissue on Wednesday, 5 July 1865:
!
  call ymdf_to_jed_common ( 1865, 7, 5, f, jed_05_07_1865 )

  if ( jed_05_07_1865 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #28
!  Deal with nonissue on Tuesday, 2 January 1866:
!
  call ymdf_to_jed_common ( 1866, 1, 2, f, jed_02_01_1866 )

  if ( jed_02_01_1866 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #29
!  Deal with nonissue on Wednesday, 2 January 1867:
!
  call ymdf_to_jed_common ( 1867, 1, 2, f, jed_02_01_1867 )

  if ( jed_02_01_1867 <= jed ) then
    issue = issue - 1
  end if
!
!  CORRECTION #30
!  Deal with the interval from Thursday, 18 September 1851
!  to Saturday, 22 April 1905.
!
!  During this period, there were no Sunday issues.
!
  call ymdf_to_jed_common ( 1905, 4, 22, f, jed_22_04_1905 )
  days = min ( jed, jed_22_04_1905 ) - jed_epoch
  sundays =  ( days + 3 ) / 7
  issue = issue - sundays
!
!  CORRECTION #31
!  Issues jumped by 500 because of mistake on 7 February 1898.
!
  call ymdf_to_jed_common ( 1898, 2, 7, f, jed_07_02_1898 )

  if ( jed_07_02_1898 <= jed ) then
    issue = issue + 500
  end if
!
!  CORRECTION #32
!  No issues from 10 August 1978 through 5 November 1978.
!
  call ymdf_to_jed_common ( 1978, 8, 10, f, jed_10_08_1978 )
  call ymdf_to_jed_common ( 1978, 11, 5, f, jed_05_11_1978 )

  if ( jed_10_08_1978 <= jed ) then
    issue = issue - ( min ( jed_05_11_1978, jed ) - jed_10_08_1978 ) - 1
  end if
!
!  CORRECTION #33
!  Issues decreased by 500 to correct previous mistake, 1 January 2000.
!
  call ymdf_to_jed_common ( 2000, 1, 1, f, jed_01_01_2000 )

  if ( jed_01_01_2000 <= jed ) then
    issue = issue - 500
  end if

  return
end
subroutine jed_to_nyt_issue_ideal ( jed, issue_ideal )

!*****************************************************************************80
!
!! JED_TO_NYT_ISSUE_IDEAL returns an ideal issue number for the New York Times.
!
!  Discussion:
!
!    The New York Times began publication with Volume 1, Issue 1 on
!    Thursday, 18 September 1851.
!
!    Each issue was assigned a volume number and an issue number.
!    The issue number did not restart at 1 with each new volume, but
!    continued to increment.
!
!    A new of interruptions, mistakes, and unusual conventions meant that
!    it is actually not easy to determine the issue number from the date;
!    moreover, it's not easy to determine the "real" issue number from
!    the nominal issue number.
!
!    For instance, in the early years, the paper would not publish
!    on a day shortly after New Year's day (usually 02 January, but
!    not always.)  It would similarly not publish on a day shortly after
!    (or before!) July 4th.
!
!    There was no paper published on Sundays until 1861; the first two
!    Sunday papers came out with their own issue numbers; this was
!    "corrected" by the time the third Sunday issue came out, and thereafter
!    the Sunday paper had the same nominal issuse number as the Saturday
!    paper, until 1905.
!
!    There was a truly bizarre accident in 1898 where the issue number
!    was incremented 501 instead of by 1.
!
!    This was matched by a misguided "correction" in which the 01 January 2000
!    issue was decremented by 499 instead of being incremented by 1.
!
!    There was an anomaly in 1978 when a strike shut the paper down;
!    while the paper was shut down, time moved on but the issue number did not.
!
!    Finally, there were a number of accidents in which the issue number
!    was incorrectly assigned for a day or two and then corrected.
!
!    The attempt to compute the nominal issue number from the JED is
!    handled by the function JED_TO_NYT.
!
!    In contrast, what this function does is to try to compute an
!    ideal issue number, in which a perfect librarian has stacked a single
!    copy of each day's New York Times, and numbered them consecutively.
!    Our task is to determine, given a date, what the corresponding
!    ideal issue number is for the paper that was printed on that date,
!    or most recently and before.
!
!    Thus, if JED represents today, then JED_TO_NYT_ISSUE_IDEAL ( JED )
!    would be exactly the number of distinc issues of the New York Times that
!    have been printed.
!
!    The cover of "The Complete First Pages" seems to claim that
!    JED_TO_NYT_ISSUE_IDEAL ( 01 April 2008 ) is 54,267.
!    However, since that is actually the nominal issue number for that
!    date, this is only "nominally" true and factually false!
!
!
!    Here are a few of the misadventures in the issue numbering scheme.
!    For the purpose of the ideal issue system, the only thing that matters
!    is days on which there was no issue.  We don't care what nominal issue
!    number was printed each day, or how badly it was garbled:
!
!      Fri,  2 Jan 1852, no issue.
!            6 Jul 1852, no issue
!      Sat,  2 Jul 1853, no issue, would have been 559.
!      Mon,  4 Jul 1853, INCORRECT issue number 560 (559 not used).
!      Tue,  5 Jul 1853, correct issue number 560.
!            6 Jul 1854, issue, but same issue number as 5 Jul 1854.
!      Thu,  5 Jul 1855, issue, but same issue number as 4 Jul 1855 (#1184)
!      Tue, 25 Sep 1855, issue jumps by 2, from 1253 to 1255!
!      Sat, 29 Sep 1856, issue, but same issue number as 28 Sep 1855 (#1258).
!      Fri,  4 Jan 1856, issue, but same issue number as 3 Jan 1856, (#1340).
!      Mon,  7 Jul 1856, issue, but same issue number as 5 Jul 1856, (#1497).
!      Sat,  3 Jan 1857, issue, but same issue number as 2 Jan 1857, (#1651).
!      Sat,  2 Jan 1858, issue, but same issue number as 1 Jan 1858, (#1962).
!      Tue,  6 Jul 1858, issue, but same issue number as 5 Jul 1858, (#2119).
!      Tue,  5 Jul 1859, no issue.
!      Tue,  3 Jan 1860, no issue.
!      Thu,  5 Jul 1860, no issue.
!      Wed,  2 Jan 1861, no issue
!      Sun, 21 Apr 1861, first Sunday issue.  First two Sundays get distinct
!                        issue numbers.  Thereafter, a "correction" is made.
!      Fri,  5 Jul 1861, no issue.
!      Thu,  2 Jan 1862, no issue.
!      Sat,  5 Jul 1862, no issue.
!      Fri,  2 Jan 1863, no issue.
!      Sat,  2 Jan 1864, no issue.
!      Tue,  5 Jul 1864, no issue.
!      Wed,  5 Jul 1865, no issue.
!      Tue,  2 Jan 1866, no issue.
!      Wed,  2 Jan 1867, no issue.
!      Sat,  5 Feb 1898, issue 14499.
!      Mon,  7 Feb 1898, issue 15000 (incremented by 501 instead of by 1)
!      Sun, 23 Apr 1905, Sunday paper gets distinct issue number.
!      Wed,  9 Aug 1978, last prestrike issue.  Issue numbers halt.
!      Mon,  6 Nov 1978, first poststrike issue, issue numbers resume.
!      Sat,  1 Jan 2000, issue numbers "corrected" downwards by 499.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Anonymous,
!    A Correction Welcome to 51,254,
!    The New York Times,
!    01 January 2000, Volume 149, Issue 51254.
!
!    James Barron,
!    What's in a Number? 143 Years of News,
!    The New York Times,
!    14 March 1995, Volume 144, Issue 50000.
!
!    The New York Times,
!    Page One, 1896-1996, A Special Commemorative Edition Celebrating the
!    100th Anniversary of the Purchase of the New York Times by Adolph S Ochs,
!    Galahad Books, 1996,
!    ISBN: 0-88365-961-1,
!    LC: D411.P25.
!
!    The New York Times,
!    The Complete First Pages, 1851-2008,
!    Black Dog & Leventhal Publishers, 2008,
!    ISBN13: 978-1-57912-749-7,
!    LC: D351.N53.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) ISSUE_IDEAL, the ideal NYT issue number.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) days
  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  integer ( kind = 4 ) issue_ideal
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_02_01_1852
  real ( kind = 8 ) jed_07_06_1852
  real ( kind = 8 ) jed_07_02_1853
  real ( kind = 8 ) jed_07_06_1854
  real ( kind = 8 ) jed_07_05_1855
  real ( kind = 8 ) jed_25_09_1855
  real ( kind = 8 ) jed_29_09_1855
  real ( kind = 8 ) jed_04_01_1856
  real ( kind = 8 ) jed_07_07_1856
  real ( kind = 8 ) jed_03_01_1857
  real ( kind = 8 ) jed_02_01_1858
  real ( kind = 8 ) jed_06_07_1858
  real ( kind = 8 ) jed_05_07_1859
  real ( kind = 8 ) jed_03_01_1860
  real ( kind = 8 ) jed_05_07_1860
  real ( kind = 8 ) jed_02_01_1861
  real ( kind = 8 ) jed_20_04_1861
  real ( kind = 8 ) jed_04_21_1861
  real ( kind = 8 ) jed_04_28_1861
  real ( kind = 8 ) jed_05_05_1861
  real ( kind = 8 ) jed_05_07_1861
  real ( kind = 8 ) jed_02_01_1862
  real ( kind = 8 ) jed_05_07_1862
  real ( kind = 8 ) jed_02_01_1863
  real ( kind = 8 ) jed_28_09_1863
  real ( kind = 8 ) jed_30_09_1863
  real ( kind = 8 ) jed_02_01_1864
  real ( kind = 8 ) jed_05_07_1864
  real ( kind = 8 ) jed_03_01_1865
  real ( kind = 8 ) jed_05_07_1865
  real ( kind = 8 ) jed_02_01_1866
  real ( kind = 8 ) jed_02_01_1867
  real ( kind = 8 ) jed_07_02_1898
  real ( kind = 8 ) jed_22_04_1905
  real ( kind = 8 ) jed_10_08_1978
  real ( kind = 8 ) jed_05_11_1978
  real ( kind = 8 ) jed_01_01_2000
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) sundays
  integer ( kind = 4 ) volume
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!    The "epoch" of the NYT calendar is the mythical date when issue "0"
!    would have been printed, namely, a tad past midnight, 17 September 1851.
!
!    Volume #1, Issue #1 was printed on 18 September 1851.
!
  y2 = 1851
  m2 = 9
  d2 = 17
  f2 = 0.0D+00
  call ymdf_to_jed_common ( y2, m2, d2, f2, jed_epoch )

  if ( jed <= jed_epoch ) then
    issue_ideal = 0
    return
  end if

  issue_ideal = jed - jed_epoch

  f = 0.0D+00
!
!  CORRECTION #1
!  Deal with nonissue on Friday, 2 January 1852
!
  call ymdf_to_jed_common ( 1852, 1, 2, f, jed_02_01_1852 )

  if ( jed_02_01_1852 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #2
!  Deal with nonissue on Tuesday, 6 July 1852
!
  call ymdf_to_jed_common ( 1852, 7, 6, f, jed_07_06_1852 )

  if ( jed_07_06_1852 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #3
!  Deal with nonissue on Saturday, 2 July 1853
!
  call ymdf_to_jed_common ( 1853, 7, 2, f, jed_07_02_1853 )

  if ( jed_07_02_1853 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #12
!  Deal with nonissue on Tuesday, 5 July 1859:
!
  call ymdf_to_jed_common ( 1859, 7, 5, f, jed_05_07_1859 )

  if ( jed_05_07_1859 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #13
!  Deal with nonissue on Tuesday, 3 January 1860:
!
  call ymdf_to_jed_common ( 1860, 1, 3, f, jed_03_01_1860 )

  if ( jed_03_01_1860 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #14
!  Deal with nonissue on Thursday, 5 July 1860:
!
  call ymdf_to_jed_common ( 1860, 7, 5, f, jed_05_07_1860 )

  if ( jed_05_07_1860 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #15
!  Deal with nonissue on Wednesday, 2 January 1861:
!
  call ymdf_to_jed_common ( 1861, 1, 2, f, jed_02_01_1861 )

  if ( jed_02_01_1861 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #19
!  Deal with nonissue on Friday, 5 July 1861:
!
  call ymdf_to_jed_common ( 1861, 7, 5, f, jed_05_07_1861 )

  if ( jed_05_07_1861 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #20
!  Deal with nonissue on Thursday, 2 January 1862:
!
  call ymdf_to_jed_common ( 1862, 1, 2, f, jed_02_01_1862 )

  if ( jed_02_01_1862 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #21
!  Deal with nonissue on Saturday, 5 July 1862:
!
  call ymdf_to_jed_common ( 1862, 7, 5, f, jed_05_07_1862 )

  if ( jed_05_07_1862 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #22
!  Deal with nonissue on Friday, 2 January 1863:
!
  call ymdf_to_jed_common ( 1863, 1, 2, f, jed_02_01_1863 )

  if ( jed_02_01_1863 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #24
!  Deal with nonissue on Saturday, 2 January 1864:
!
  call ymdf_to_jed_common ( 1864, 1, 2, f, jed_02_01_1864 )

  if ( jed_02_01_1864 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #25
!  Deal with nonissue on Tuesday, 5 July 1864:
!
  call ymdf_to_jed_common ( 1864, 7, 5, f, jed_05_07_1864 )

  if ( jed_05_07_1864 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #26
!  Deal with nonissue on Monday, 3 January 1865:
!
  call ymdf_to_jed_common ( 1865, 1, 3, f, jed_03_01_1865 )

  if ( jed_03_01_1865 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #27
!  Deal with nonissue on Wednesday, 5 July 1865:
!
  call ymdf_to_jed_common ( 1865, 7, 5, f, jed_05_07_1865 )

  if ( jed_05_07_1865 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #28
!  Deal with nonissue on Tuesday, 2 January 1866:
!
  call ymdf_to_jed_common ( 1866, 1, 2, f, jed_02_01_1866 )

  if ( jed_02_01_1866 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #29
!  Deal with nonissue on Wednesday, 2 January 1867:
!
  call ymdf_to_jed_common ( 1867, 1, 2, f, jed_02_01_1867 )

  if ( jed_02_01_1867 <= jed ) then
    issue_ideal = issue_ideal - 1
  end if
!
!  CORRECTION #30
!  Deal with the interval from Thursday, 18 September 1851
!  to Saturday, 20 April 1861.
!
!  During this period, there were no Sunday issues.
!
  call ymdf_to_jed_common ( 1861, 4, 20, f, jed_20_04_1861 )
  days = min ( jed, jed_20_04_1861 ) - jed_epoch
  sundays =  ( days + 3 ) / 7
  issue_ideal = issue_ideal - sundays
!
!  CORRECTION #32
!  No issues from 10 August 1978 through 5 November 1978.
!
  call ymdf_to_jed_common ( 1978, 8, 10, f, jed_10_08_1978 )
  call ymdf_to_jed_common ( 1978, 11, 5, f, jed_05_11_1978 )

  if ( jed_10_08_1978 <= jed ) then
    issue_ideal = issue_ideal &
      - ( min ( jed_05_11_1978, jed ) - jed_10_08_1978 ) - 1
  end if

  return
end
subroutine jed_to_ymdf_common ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_COMMON converts a JED to a Common YMDF date.
!
!  Discussion:
!
!    The "common" calendar is meant to be the calendar which is Julian up to
!    JED = 2299160.5, and Gregorian thereafter.
!
!    There is no year 0.  BC years are specified using a negative value.
!
!  Example:
!
!        JED            Y    M   D
!    -------    ------------------
!          0    BCE  4713  Jan   1
!    2440000    CE   1968  May  23
!    2446065    CE   1984  Dec  31
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, integr ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F, the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_transition
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  call transition_to_jed_common ( jed_transition )

  if ( jed <= jed_transition ) then
    call jed_to_ymdf_julian ( jed, y, m, d, f )
  else
    call jed_to_ymdf_gregorian ( jed, y, m, d, f )
  end if

  return
end
subroutine jed_to_ymdf_gregorian ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_GREGORIAN converts a JED to a Gregorian YMDF date.
!
!  Discussion:
!
!    This Gregorian calendar is extended backwards in time before
!    its actual adoption.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Algorithm F,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, pages 324-325.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  integer ( kind = 4 ) g
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) j_prime
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) t_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
!
!  Determine the computational date (Y'/M'/D').
!
  j = int ( jed + 0.5D+00 )
  f = ( jed + 0.5D+00 ) - real ( j, kind = 8 )

  g = ( 4 * j + 274277 ) / 146097
  g = ( 3 * g ) / 4 - 38

  j_prime = j + 1401 + g

  y_prime =     ( 4 * j_prime + 3 ) / 1461
  t_prime = mod ( 4 * j_prime + 3, 1461 ) / 4
  m_prime =     ( 5 * t_prime + 2 ) / 153
  d_prime = mod ( 5 * t_prime + 2, 153 ) / 5
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime + 2, 12 ) + 1
  y = y_prime - 4716 + ( 14 - m ) / 12
!
!  Any year before 1 AD must be moved one year further back, since
!  this calendar does not include a year 0.
!
  call y_astronomical_to_common ( y, y )

  return
end
subroutine jed_to_ymdf_julian ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_JULIAN converts a JED to a Julian YMDF date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Algorithm F,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, pages 324-325.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F, the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) j_prime
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) t_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
!
!  Determine the computational date (Y'/M'/D').
!
  j = int ( jed + 0.5D+00 )
  f = ( jed + 0.5D+00 ) - real ( j, kind = 8 )

  j_prime = j + 1401

  y_prime =     ( 4 * j_prime + 3 ) / 1461
  t_prime = mod ( 4 * j_prime + 3, 1461 ) / 4
  m_prime =     ( 5 * t_prime + 2 ) / 153
  d_prime = mod ( 5 * t_prime + 2, 153 ) / 5
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime + 2, 12 ) + 1
  y = y_prime - 4716 + ( 14 - m ) / 12
!
!  Any year before 1 AD must be moved one year further back, since
!  this calendar does not include a year 0.
!
  call y_astronomical_to_common ( y, y )

  return
end
subroutine month_borrow_common ( y, m )

!*****************************************************************************80
!
!! MONTH_BORROW_COMMON borrows a year of months on the Common calendar.
!
!  Discussion:
!
!    If the month index is legal, nothing is done.  If the month index
!    is too small, then one or more years are "cashed in" to bring the
!    month index up to a legal value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M, the YM date.
!    On input, M might be negative.  On output, Y should have decreased by
!    one, and M gone up by the number of months in the year that we
!    "cashed in".  The routine knows there was no year 0.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_common

  do while ( m <= 0 )

    months = year_length_months_common ( y )

    m = m + months
    y = y - 1

    if ( y == 0 ) then
      y = - 1
    end if

  end do

  return
end
subroutine month_borrow_gregorian ( y, m )

!*****************************************************************************80
!
!! MONTH_BORROW_GREGORIAN borrows a year of months on the Gregorian calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M, the YM date.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_gregorian

  do while ( m <= 0 )

    months = year_length_months_gregorian ( y )

    m = m + months
    y = y - 1

    if ( y == 0 ) then
      y = - 1
    end if

  end do

  return
end
subroutine month_borrow_julian ( y, m )

!*****************************************************************************80
!
!! MONTH_BORROW_JULIAN borrows a year of months on the Julian calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M, the YM date.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_julian

  do while ( m <= 0 )

    months = year_length_months_julian ( y )

    m = m + months
    y = y - 1

    if ( y == 0 ) then
      y = - 1
    end if

  end do

  return
end
subroutine month_carry_common ( y, m )

!*****************************************************************************80
!
!! MONTH_CARRY_COMMON carries a year of months on the Common calendar.
!
!  Algorithm:
!
!    While 12 < M:
!
!      decrease M by 12;
!      increase Y by 1;
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
!    Input/output, integer ( kind = 4 ) Y, M, the year and month.
!    On output, M is no greater than 12.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_common

  do

    months = year_length_months_common ( y )

    if ( m <= months ) then
      exit
    end if

    m = m - months
    y = y + 1

  end do

  return
end
subroutine month_carry_gregorian ( y, m )

!*****************************************************************************80
!
!! MONTH_CARRY_GREGORIAN carries a year of months on the Gregorian calendar.
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
!    Input/output, integer ( kind = 4 ) Y, M, the year and month.
!    On output, M is no greater than 12.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_gregorian

  do

    months = year_length_months_gregorian ( y )

    if ( m <= months ) then
      exit
    end if

    m = m - months
    y = y + 1

  end do

  return
end
subroutine month_carry_julian ( y, m )

!*****************************************************************************80
!
!! MONTH_CARRY_JULIAN carries a year of months on the Julian calendar.
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
!    Input/output, integer ( kind = 4 ) Y, M, the year and month.
!    On output, M is no greater than 12.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_julian

  do

    months = year_length_months_julian ( y )

    if ( m <= months ) then
      exit
    end if

    m = m - months
    y = y + 1

  end do

  return
end
function month_length_common ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_COMMON returns the number of days in a Common month.
!
!  Discussion:
!
!    The "common" calendar is meant to be the calendar which is Julian up to
!    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
!
!    The routine knows that February has 28 days, except in leap years,
!    when it has 29.
!
!    In the Common calendar, October 1582 had only 21 days
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
!    Input, integer ( kind = 4 ) Y, the year in which the month occurred.
!
!    Input, integer ( kind = 4 ) M, the number of the month.
!
!    Output, integer ( kind = 4 ) MONTH_LENGTH_COMMON, the number of
!    days in the month.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(12) :: mdays = (/ &
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_common
!
!  Copy the input.
!
  m2 = m
  y2 = y
!
!  Check the input.
!
  call ym_check_common ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    month_length_common = - 1
    return
  end if
!
!  Take care of the special case.
!
  if ( y2 == 1582 ) then
    if ( m2 == 10 ) then
      month_length_common = 21
      return
    end if
  end if
!
!  Get the number of days in the month.
!
  month_length_common = mdays ( m2 )
!
!  If necessary, add 1 day for February 29.
!
  if ( m2 == 2 .and. year_is_leap_common ( y2 ) ) then
    month_length_common = month_length_common + 1
  end if

  return
end
function month_length_gregorian ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_GREGORIAN returns the number of days in a Gregorian month.
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
!    Input, integer ( kind = 4 ) Y, the year in which the month occurred.
!
!    Input, integer ( kind = 4 ) M, the number of the month.
!
!    Output, integer ( kind = 4 ) MONTH_LENGTH_GREGORIAN,
!    the number of days in the month.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(12) :: mdays = (/ &
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  integer ( kind = 4 ) month_length_gregorian
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_gregorian
!
!  Copy the input.
!
  m2 = m
  y2 = y
!
!  Check the input.
!
  call ym_check_gregorian ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    month_length_gregorian = - 1
    return
  end if
!
!  Get the number of days in the month.
!
  month_length_gregorian = mdays ( m2 )
!
!  If necessary, add 1 day for February 29.
!
  if ( m2 == 2 ) then
    if ( year_is_leap_gregorian ( y2 ) ) then
      month_length_gregorian = month_length_gregorian + 1
    end if
  end if

  return
end
function month_length_julian ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_JULIAN returns the number of days in a Julian month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year in which the month occurred.
!
!    Input, integer ( kind = 4 ) M, the number of the month.
!
!    Output, integer ( kind = 4 ) MONTH_LENGTH_JULIAN, the number of days
!    in the month.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(12) :: mdays = (/ &
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  integer ( kind = 4 ) month_length_julian
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_julian
!
!  Copy the input.
!
  m2 = m
  y2 = y
!
!  Check the input.
!
  call ym_check_julian ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    month_length_julian = - 1
    return
  end if
!
!  Get the number of days in the month.
!
  month_length_julian = mdays ( m2 )
!
!  If necessary, add 1 day for February 29.
!
  if ( m2 == 2 .and. year_is_leap_julian ( y2 ) ) then
    month_length_julian = month_length_julian + 1
  end if

  return
end
subroutine now_to_jed ( jed )

!*****************************************************************************80
!
!! NOW_TO_JED expresses the current date as JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  character ( len = 8 ) date
  real ( kind = 8 ) f
  integer ( kind = 4 ) h
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 ) time
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
  mu = values(8)

  f = real ( mu, kind = 8 )
  f = real ( s, kind = 8 ) + f / 1000.0D+00
  f = real ( n, kind = 8 ) + f / 60.0D+00
  f = real ( h, kind = 8 ) + f / 60.0D+00
  f = f / 24.0D+00

  call ymdf_to_jed_common ( y, m, d, f, jed )

  return
end
subroutine nyt_to_jed ( volume, issue, jed )

!*****************************************************************************80
!
!! NYT_TO_JED converts an NYT date to a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) VOLUME, ISSUE, the New York Times
!    volume and issue.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) issue
  real ( kind = 8 ) jed
  real ( kind = 8 ), parameter :: jed_epoch_50000 = 2449790.5D+00
  integer ( kind = 4 ) volume

  if ( 149 < volume ) then
    jed = jed_epoch_50000 + real ( issue - 50000 + 500, kind = 8 )
!
!  Take care of the bizarre case of the second half of Volume 149,
!  Jan 1 2000 to Sep 17 2000, issues 51254 through ?, which were also
!  lowered by 500.
!
  else if ( volume == 149 .and. issue < 51600 ) then
    jed = jed_epoch_50000 + real ( issue - 50000 + 500, kind = 8 )
  else if ( 44028 <= issue ) then
    jed = jed_epoch_50000 + real ( issue - 50000, kind = 8 )
!
!  Factor in the strike of 1978.
!
  else
    jed = jed_epoch_50000 + real ( issue - 50000 - 88, kind = 8 )
  end if

  return
end
subroutine nyt_to_ymd ( volume, issue, y, m, d )

!*****************************************************************************80
!
!! NYT_TO_YMD converts an NYT date to a YMD date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) VOLUME, ISSUE, the New York Times
!    volume and issue.
!
!    Output, integer ( kind = 4 ) Y, M, D, the year, month and day.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) issue
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) volume
  integer ( kind = 4 ) y

  call nyt_to_jed ( volume, issue, jed )

  call jed_to_ymdf_common ( jed, y, m, d, f )

  return
end
function r8_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R8_UNIFORM returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
subroutine transition_to_jed_common ( jed )

!*****************************************************************************80
!
!! TRANSITION_TO_JED_COMMON returns the Common calendar transition as a JED.
!
!  Discussion:
!
!    In the Common calendar, the last moment of the Julian calendar was
!      11:59 pm, 4 October 1582 Julian/CE,
!      11:59 pm, 14 October 1582 Gregorian.
!    The first minute of the Gregorian calendar ended at
!      12:01 am, 5 October 1582 Julian,
!      12:01 am, 15 October 1582 Gregorian/CE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date of the date.
!
  implicit none

  real ( kind = 8 ) jed

  jed = 2299160.5D+00

  return
end
subroutine y_astronomical_to_common ( y, y2 )

!*****************************************************************************80
!
!! Y_ASTRONOMICAL_TO_COMMON converts an Astronomical year to a Common year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the astronomical year.
!
!    Output, integer ( kind = 4 ) Y2, the Common year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2

  if ( y <= 0 ) then
    y2 = y - 1
  else
    y2 = y
  end if

  return
end
subroutine y_check_common ( y, ierror )

!*****************************************************************************80
!
!! Y_CHECK_COMMON checks a Common year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year, which must not be 0.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if Y is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y

  if ( y /= 0 ) then
    ierror = 0
  else
    ierror = 1
  end if

  return
end
subroutine y_check_gregorian ( y, ierror )

!*****************************************************************************80
!
!! Y_CHECK_GREGORIAN checks a Gregorian year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year, which must not be 0.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if Y is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y

  if ( y /= 0 ) then
    ierror = 0
  else
    ierror = 1
  end if

  return
end
subroutine y_check_julian ( y, ierror )

!*****************************************************************************80
!
!! Y_CHECK_JULIAN checks a Julian year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year, which must not be 0.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if Y is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y

  if ( y /= 0 ) then
    ierror = 0
  else
    ierror = 1
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
!    23 July 2000
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
    y2 = - huge ( 1 )
  else
    y2 = y
  end if

  return
end
function year_is_leap_common ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_COMMON returns TRUE if the Common year was a leap year.
!
!  Discussion:
!
!    The "common" calendar is meant to be the calendar which is Julian up to
!    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
!
!  Algorithm:
!
!    If ( the year is less than 0 ) then
!
!      if the year+1 is divisible by 4 then
!        the year is a leap year.
!
!    else if ( the year is 0 ) then
!
!      the year is not a leap year ( in fact, it's illegal )
!
!    else if ( the year is no greater than 1582 ) then
!
!      if the year is divisible by 4 then
!        the year is a leap year.
!
!    else if (
!      the year is divisible by 4 and
!      ( the year is not divisible by 100
!      or
!      the year is divisible by 400 )
!      ) then
!        the year is a leap year.
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
!    Output, logical YEAR_IS_LEAP_COMMON, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_common

  if ( y == 0 ) then
    year_is_leap_common = .false.
    return
  end if
!
!  BC years have to have 1 added to them to make a proper leap year evaluation.
!
  call y_common_to_astronomical ( y, y2 )

  if ( y2 <= 1582 ) then

    if ( i4_modp ( y2, 4 ) == 0 ) then
      year_is_leap_common = .true.
    else
      year_is_leap_common = .false.
    end if

  else

    if ( i4_modp ( y2, 400 ) == 0 ) then
      year_is_leap_common = .true.
    else if ( i4_modp ( y2, 100 ) == 0 ) then
      year_is_leap_common = .false.
    else if ( i4_modp ( y2, 4 ) == 0 ) then
      year_is_leap_common = .true.
    else
      year_is_leap_common = .false.
    end if

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
  integer ( kind = 4 ) y2
  logical year_is_leap_gregorian

  if ( y == 0 ) then
    year_is_leap_gregorian = .false.
    return
  end if
!
!  BC years have to have 1 added to them to make a proper leap year evaluation.
!
  call y_common_to_astronomical ( y, y2 )

  if ( mod ( y2, 400 ) == 0 ) then
    year_is_leap_gregorian = .true.
  else if ( mod ( y2, 100 ) == 0 ) then
    year_is_leap_gregorian = .false.
  else if ( mod ( y2, 4 ) == 0 ) then
    year_is_leap_gregorian = .true.
  else
    year_is_leap_gregorian = .false.
  end if

  return
end
function year_is_leap_julian ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_JULIAN returns TRUE if the Julian year was a leap year.
!
!  Algorithm:
!
!    If ( Y < 0 and Y+1 is divisible by 4 ) then
!      the year is a leap year.
!    else if ( Y == 0 ) then
!      the year is illegal
!    else if ( 0 < Y and Y is divisible by 4 ) then
!      the year is a leap year.
!    else
!      the year is NOT a leap year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_JULIAN, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_julian

  if ( y == 0 ) then
    year_is_leap_julian = .false.
    return
  end if

  call y_common_to_astronomical ( y, y2 )

  if ( i4_modp ( y2, 4 ) == 0 ) then
    year_is_leap_julian = .true.
  else
    year_is_leap_julian = .false.
  end if

  return
end
function year_length_months_common ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_COMMON returns the number of months in a Common year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_COMMON, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_common

  year_length_months_common = 12

  return
end
function year_length_months_gregorian ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_GREGORIAN returns the number of months in a Gregorian year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer Y, the year to be checked.
!
!    Output, integer YEAR_LENGTH_MONTHS_GREGORIAN, the number of months
!    in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_gregorian

  year_length_months_gregorian = 12

  return
end
function year_length_months_julian ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_JULIAN returns the number of months in a Julian year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_JULIAN, the number of months
!    in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_julian

  year_length_months_julian = 12

  return
end
subroutine ym_check_common ( y, m, ierror )

!*****************************************************************************80
!
!! YM_CHECK_COMMON checks a Common YM date.
!
!  Discussion:
!
!    If the month is less than 1, then the month is incremented
!    by 12, and the year decremented by 1, repeatedly, until
!    the month is greater than or equal to 1.
!
!    If the month is greater than 12, then the month is decremented
!    by 12, and the year incremented by 1, repeatedly, until the
!    month is less than or equal to 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found,
!    and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_common ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure the month isn't too small or too big.
!
  call month_borrow_common ( y, m )

  call month_carry_common ( y, m )

  return
end
subroutine ym_check_gregorian ( y, m, ierror )

!*****************************************************************************80
!
!! YM_CHECK_GREGORIAN checks a Gregorian YM date.
!
!  Discussion:
!
!    If the month is less than 1, then the month is incremented
!    by 12, and the year decremented by 1, repeatedly, until
!    the month is greater than or equal to 1.
!
!    If the month is greater than 12, then the month is decremented
!    by 12, and the year incremented by 1, repeatedly, until the
!    month is less than or equal to 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found,
!    and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_gregorian ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure the month isn't too small or too big.
!
  call month_borrow_gregorian ( y, m )

  call month_carry_gregorian ( y, m )

  return
end
subroutine ym_check_julian ( y, m, ierror )

!*****************************************************************************80
!
!! YM_CHECK_JULIAN checks a Julian YM date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the date,
!    and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_julian ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure the month isn't too small or too big.
!
  call month_borrow_julian ( y, m )

  call month_carry_julian ( y, m )

  return
end
subroutine ymd_check_common ( y, m, d, ierror )

!*****************************************************************************80
!
!! YMD_CHECK_COMMON checks a Common YMD date.
!
!  Discussion:
!
!    Certain simple errors in dates will be corrected, such as
!      "31 September 1996"
!    which will become
!      "1 October 1996".
!
!    The routine also knows that in the Common calendar, the dates
!    5 October 1582 through 14 October 1582 are illegal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date,
!    which may be corrected if necessary and possible.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if the date is legal.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  character ( len = 30 ) s
  integer ( kind = 4 ) y

  ierror = 0
!
!  Check the year.
!
  if ( y == 0 ) then
    ierror = 1
    return
  end if
!
!  Check the month.
!
  call month_borrow_common ( y, m )

  call month_carry_common ( y, m )
!
!  Check the day.
!
  call day_borrow_common ( y, m, d )

  call day_carry_common ( y, m, d )
!
!  Now make sure that the date does not fall in the
!  Julian-to-Gregorian calendar switchover limbo.
!
  if ( y == 1582 ) then
    if ( m == 10 ) then
      if ( 5 <= d .and. d <= 14 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'YMD_CHECK_COMMON - Warning!'
        write ( *, '(a)' ) '  Illegal date:'
        call ymd_to_s_numeric ( y, m, d, s )
        write ( *, '(4x,a)' ) trim ( s )
      end if
    end if
  end if

  return
end
subroutine ymd_check_gregorian ( y, m, d, ierror )

!*****************************************************************************80
!
!! YMD_CHECK_GREGORIAN checks a Gregorian YMD date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date, which may
!    be corrected if necessary and possible.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if the date is legal.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the month.
!
  call ym_check_gregorian ( y, m, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Check the day.
!
  call day_borrow_gregorian ( y, m, d )

  call day_carry_gregorian ( y, m, d )

  return
end
subroutine ymd_check_julian ( y, m, d, ierror )

!*****************************************************************************80
!
!! YMD_CHECK_JULIAN checks a Julian YMD date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date, which may
!    be corrected if necessary and possible.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if the date is legal.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the month.
!
  call ym_check_julian ( y, m, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Check the day.
!
  call day_borrow_julian ( y, m, d )

  call day_carry_julian ( y, m, d )

  return
end
subroutine ymd_to_nyt ( y, m, d, volume, issue )

!*****************************************************************************80
!
!! YMD_TO_NYT converts a YMD date to an NYT date.
!
!  Discussion:
!
!    The New York Times began publication with Volume 1, Issue 1 on
!    Thursday, 18 September 1851.
!
!    The Volume number is incremented annually, on 18 September.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Anonymous,
!    A Correction; Welcome to 51,254,
!    The New York Times,
!    01 January 2000, Volume 149, Issue 51254.
!
!    James Barron,
!    What's in a Number? 143 Years of News,
!    The New York Times,
!    14 March 1995, Volume 144, Issue 50000.
!
!    The New York Times,
!    Page One, 1896-1996, A Special Commemorative Edition Celebrating the
!    100th Anniversary of the Purchase of the New York Times by Adolph S Ochs,
!    Galahad Books, 1996,
!    ISBN: 0-88365-961-1,
!    LC: D411.P25.
!
!    The New York Times,
!    The Complete First Pages, 1851-2008,
!    Black Dog & Leventhal Publishers, 2008,
!    ISBN13: 978-1-57912-749-7,
!    LC: D351.N53.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, the YMD date.
!
!    Output, integer ( kind = 4 ) VOLUME, ISSUE, the New York Times
!   volume and issue.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) issue
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) volume
  integer ( kind = 4 ) y

  f = 0.0D+00

  call ymdf_to_jed_common ( y, m, d, f, jed )

  call jed_to_nyt ( jed, volume, issue )

  return
end
subroutine ymd_to_s_numeric ( y, m, d, s )

!*****************************************************************************80
!
!! YMD_TO_S_NUMERIC writes a YMD date into a string.
!
!  Format:
!
!       YYYY/MM/DD
!    or
!      -YYYY/MM/DD
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2000
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
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  character ( len = 20 ) s1
  character ( len = 2 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Make local copies of the input.
!
  y2 = y
  m2 = m
  d2 = d

  call i4_to_s_left ( y2, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( m2, s2 )

  call s_cat ( s1, s2, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( d2, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine ymdf_check_common ( y, m, d, f, ierror )

!*****************************************************************************80
!
!! YMDF_CHECK_COMMON checks a Common YMDF date.
!
!  Discussion:
!
!    Certain simple errors in dates will be corrected, such as
!      "31 September 1996"
!    which will become
!      "1 October 1996".
!
!    The routine also knows that in the Common calendar, the dates
!    5 October 1582 through 14 October 1582 are illegal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F, the
!    YMDF date, which may be corrected if necessary and possible.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if the date is legal.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  ierror = 0

  call ymd_check_common ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call frac_borrow_common ( y, m, d, f )

  call frac_carry_common ( y, m, d, f )

  return
end
subroutine ymdf_check_julian ( y, m, d, f, ierror )

!*****************************************************************************80
!
!! YMDF_CHECK_JULIAN checks a Julian YMDF date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F, the
!    YMDF date, which may be corrected if necessary and possible.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if the date is legal.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  ierror = 0

  call ymd_check_julian ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call frac_borrow_julian ( y, m, d, f )

  call frac_carry_julian ( y, m, d, f )

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
!    Input, integer ( kind = 4 ) Y1, integer ( kind = 4 ) M1,
!    integer ( kind = 4 ) D1, real ( kind = 8 ) F1, the
!    first YMDF date.
!
!    Input, integer ( kind = 4 ) Y2, integer ( kind = 4 ) M2,
!    integer ( kind = 4 ) D2, real ( kind = 8 ) F2, the
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
!    Input, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F, the YMDF date.
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

  call ymdf_check_common ( y1, m1, d1, f1, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if

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
!    Input, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F, the YMDF date.
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
!  Check the date.
!
  call ymd_check_gregorian ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
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

  g = ( y_prime + 184 ) / 100
  g = ( 3 * g ) / 4 - 38

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
!    Input, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F, the YMDF date.
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
!  Check the date.
!
  call ymdf_check_julian ( y, m, d, f, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
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
subroutine ymdf_to_s_common ( y, m, d, f, s )

!*****************************************************************************80
!
!! YMDF_TO_S_COMMON writes a Common YMDF date into a string.
!
!  Format:
!
!    CE YYYY/MM/DD.FF
!    BCE YYYY/MM/DD.FF
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
!    Input, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, character ( len = * ) S, a representation of the date.
!
  implicit none

  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) d2
  real      ( kind = 8 ) f
  real      ( kind = 8 ) f2
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) m2
  character ( len = 20 ) s1
  character ( len = 2 )  s2
  character ( len = * )  s
  integer   ( kind = 4 ) y
  integer   ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y
  m2 = m
  d2 = d
  f2 = f
!
!  Check the input.
!
  call ymdf_check_common ( y2, m2, d2, f2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

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

  call frac_to_s ( f2, s1 )

  call s_cat ( s, s1(1:3), s )

  return
end
