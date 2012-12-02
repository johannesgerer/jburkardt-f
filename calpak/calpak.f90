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
subroutine cws_to_jed_gps ( c, w, s, jed )

!*****************************************************************************80
!
!! CWS_TO_JED_GPS converts a GPS CWS date to a JED.
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
!    Input, integer ( kind = 4 ) C, integer ( kind = 4 ) W, real ( kind = 8 ) S,
!    the GPS cycle/week/second date.
!
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) s
  integer ( kind = 4 ) w

  call epoch_to_jed_gps ( jed_epoch )

  d = real ( 7 * ( 1024 * c + w ), kind = 8 ) &
    + s / ( 24.0D+00 * 60.0D+00 * 60.0D+00 )

  jed = jed_epoch + d

  return
end
subroutine cws_to_s_gps ( c, w, sec, s )

!*****************************************************************************80
!
!! CWS_TO_S_GPS writes a GPS CWS date into a string.
!
!  Format:
!
!    CC/WWWW/SSSSSS.SS GPS
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
!    Input, integer ( kind = 4 ) C,  integer ( kind = 4 ) W,
!    real ( kind = 8 ) SEC, the GPS cycle/week/second date.
!
!    Output, character ( len = * ) S, a representation of the date.
!
  implicit none

  integer ( kind = 4 ) c
  character ( len = 25 ) s1
  character ( len = 4 ) s2
  character ( len = 9 ) s3
  character ( len = * ) s
  real ( kind = 8 ) sec
  integer ( kind = 4 ) w

  call i4_to_s_left ( c, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( w, s2 )

  call s_cat ( s1, s2, s1 )

  call s_cat ( s1, '/', s1 )

  write ( s3, '(f9.2)' ) sec
  s3 = adjustl ( s3 )

  call s_cat ( s1, s3, s )

  call s_cat ( s, ' GPS', s )

  return
end
subroutine datenum_to_jed ( dn, jed )

!*****************************************************************************80
!
!! DATENUM_TO_JED converts a MATLAB date number to a JED.
!
!  Discussion:
!
!    The MATLAB "datenum" function accepts a string defining
!    a date and returns a MATLAB date number:
!
!      dn = datenum ( 'Aug 17 1939' )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DN, a MATLAB date number.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  real ( kind = 8 ) dn
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch

  call epoch_to_jed_matlab ( jed_epoch )

  jed = dn + jed_epoch

  return
end
subroutine day_borrow_alexandrian ( y, m, d )

!*****************************************************************************80
!
!! DAY_BORROW_ALEXANDRIAN borrows days from months in an Alexandrian date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
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
  integer ( kind = 4 ) month_length_alexandrian
  integer ( kind = 4 ) y

  do while ( d <= 0 )

    m = m - 1

    call month_borrow_alexandrian ( y, m )

    days = month_length_alexandrian ( y, m )

    d = d + days

  end do

  return
end
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
subroutine day_borrow_eg_civil ( y, m, d )

!*****************************************************************************80
!
!! DAY_BORROW_EG_CIVIL borrows days from months in an Egyptian Civil date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, a year, month, and day
!    representing a date.  On input, D might be negative.  On output,
!    M should have decreased by one month, and D gone up by the
!    number of days in the month we "cashed in".  Y may be affected
!    if the input value of M was 1.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_eg_civil
  integer ( kind = 4 ) y

  do while ( d <= 0 )

    m = m - 1

    call month_borrow_eg_civil ( y, m )

    days = month_length_eg_civil ( y, m )

    d = d + days

  end do

  return
end
subroutine day_borrow_english ( y, m, d )

!*****************************************************************************80
!
!! DAY_BORROW_ENGLISH borrows days from months in an English date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, a year, month, and day
!    representing a date.  On input, D might be negative.  On output,
!    M should have decreased by one month, and D gone up by the
!    number of days in the month we "cashed in".  Y may be affected
!    if the input value of M was 1.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_english
  integer ( kind = 4 ) y

  do while ( d <= 0 )

    m = m - 1

    call month_borrow_english ( y, m )

    days = month_length_english ( y, m )

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
!    Input/output, integer ( kind = 4 ) Y, M, D, a year, month, and day
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
subroutine day_borrow_hebrew ( y, m, d )

!*****************************************************************************80
!
!! DAY_BORROW_HEBREW borrows days from months in a Hebrew date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, a year, month, and day
!    representing a date.  On input, D might be negative.  On output,
!    M should have decreased by one month, and D gone up by the
!    number of days in the month we "cashed in".  Y may be affected
!    if the input value of M was 1.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_hebrew
  integer ( kind = 4 ) y

  do while ( d <= 0 )

    m = m - 1

    call month_borrow_hebrew ( y, m )

    days = month_length_hebrew ( y, m )

    d = d + days

  end do

  return
end
subroutine day_borrow_islamic ( y, m, d )

!*****************************************************************************80
!
!! DAY_BORROW_ISLAMIC borrows days from months in an Islamic date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, a year, month, and day
!    representing a date.  On input, D might be negative.  On output,
!    M should have decreased by one month, and D gone up by the
!    number of days in the month we "cashed in".  Y may be affected
!    if the input value of M was 1.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_islamic
  integer ( kind = 4 ) y

  do while ( d <= 0 )

    m = m - 1

    call month_borrow_islamic ( y, m )

    days = month_length_islamic ( y, m )

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
!    Input/output, integer ( kind = 4 ) Y, M, D, a year, month, and day
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
subroutine day_borrow_republican ( y, m, d )

!*****************************************************************************80
!
!! DAY_BORROW_REPUBLICAN borrows days from months in a Republican date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, a year, month, and day
!    representing a date.  On input, D might be negative.  On output,
!    M should have decreased by one month, and D gone up by the
!    number of days in the month we "cashed in".  Y may be affected
!    if the input value of M was 1.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_republican
  integer ( kind = 4 ) y

  do while ( d <= 0 )

    m = m - 1

    call month_borrow_republican ( y, m )

    days = month_length_republican ( y, m )

    d = d + days

  end do

  return
end
subroutine day_borrow_roman ( y, m, d )

!*****************************************************************************80
!
!! DAY_BORROW_ROMAN borrows days from months in a Roman date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, a year, month, and day
!    representing a date.  On input, D might be negative.  On output,
!    M should have decreased by one month, and D gone up by the
!    number of days in the month we "cashed in".  Y may be affected
!    if the input value of M was 1.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_roman
  integer ( kind = 4 ) y

  do while ( d <= 0 )

    m = m - 1

    call month_borrow_roman ( y, m )

    days = month_length_roman ( y, m )

    d = d + days

  end do

  return
end
subroutine day_carry_alexandrian ( y, m, d )

!*****************************************************************************80
!
!! DAY_CARRY_ALEXANDRIAN carries days to months in an Alexandrian date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
!    On output, D is between 1 and the number of days in M.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_alexandrian
  integer ( kind = 4 ) y

  days = month_length_alexandrian ( y, m )

  do while ( days < d )

    d = d - days
    m = m + 1
    days = month_length_alexandrian ( y, m )
!
!  Make sure the month isn't too big.
!
    call month_carry_alexandrian ( y, m )

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
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
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
subroutine day_carry_eg_civil ( y, m, d )

!*****************************************************************************80
!
!! DAY_CARRY_EG_CIVIL carries days to months in an Egyptian Civil date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
!    On output, D is between 1 and the number of days in M.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_eg_civil
  integer ( kind = 4 ) y

  days = month_length_eg_civil ( y, m )

  do while ( days < d )

    d = d - days
    m = m + 1
    days = month_length_eg_civil ( y, m )
!
!  Make sure the month isn't too big.
!
    call month_carry_eg_civil ( y, m )

  end do

  return
end
subroutine day_carry_english ( y, m, d )

!*****************************************************************************80
!
!! DAY_CARRY_ENGLISH carries days to months in an English date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
!    On output, D is between 1 and the number of days in M.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_english
  integer ( kind = 4 ) y
!
!  If the date is in the transition month, deflate it,
!  so we can perform ordinary arithmetic.
!
  call deflate_english ( y, m, d )

  days = month_length_english ( y, m )

  do while ( days < d )

    d = d - days
    m = m + 1
    days = month_length_english ( y, m )
!
!  Make sure the month isn't too big.
!
    call month_carry_english ( y, m )

  end do
!
!  If the date is in the transition month, inflate it.
!
  call inflate_english ( y, m, d )

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
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
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
subroutine day_carry_hebrew ( y, m, d )

!*****************************************************************************80
!
!! DAY_CARRY_HEBREW carries days to months in a Hebrew date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
!    On output, D is between 1 and the number of days in M.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_hebrew
  integer ( kind = 4 ) y

  days = month_length_hebrew ( y, m )

  do while ( days < d )

    d = d - days
    m = m + 1
    days = month_length_hebrew ( y, m )
!
!  Make sure the month isn't too big.
!
    call month_carry_hebrew ( y, m )

  end do

  return
end
subroutine day_carry_islamic ( y, m, d )

!*****************************************************************************80
!
!! DAY_CARRY_ISLAMIC carries days to months in an Islamic date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
!    On output, D is between 1 and the number of days in M.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_islamic
  integer ( kind = 4 ) y

  days = month_length_islamic ( y, m )

  do while ( days < d )

    d = d - days
    m = m + 1
    days = month_length_islamic ( y, m )
!
!  Make sure the month isn't too big.
!
    call month_carry_islamic ( y, m )

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
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
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
subroutine day_carry_republican ( y, m, d )

!*****************************************************************************80
!
!! DAY_CARRY_REPUBLICAN carries days to months in a Republican date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
!    On output, D is between 1 and the number of days in M.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_republican
  integer ( kind = 4 ) y

  days = month_length_republican ( y, m )

  do while ( days < d )

    d = d - days
    m = m + 1
    days = month_length_republican ( y, m )
!
!  Make sure the month isn't too big.
!
    call month_carry_republican ( y, m )

  end do

  return
end
subroutine day_carry_roman ( y, m, d )

!*****************************************************************************80
!
!! DAY_CARRY_ROMAN carries days to months in a Roman date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
!    On output, D is between 1 and the number of days in M.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_roman
  integer ( kind = 4 ) y

  days = month_length_roman ( y, m )

  do while ( days < d )

    d = d - days
    m = m + 1
    days = month_length_roman ( y, m )
!
!  Make sure the month isn't too big.
!
    call month_carry_roman ( y, m )

  end do

  return
end
subroutine day_list_common ( y1, m1, d1, y2, m2, d2 )

!*****************************************************************************80
!
!! DAY_LIST_COMMON prints a list of days between two dates.
!
!  Discussion:
!
!    Given the dates of September 25, 2005 and October 2, 2005,
!    the routine should print out:
!
!    Sun, Sep 25 2005 -
!    Mon, Sep 26 2005 -
!    Tue, Sep 27 2005 -
!    Wed, Sep 28 2005 -
!    Thu, Sep 29 2005 -
!    Fri, Sep 30 2005 -
!    Sat, Oct 01 2005 -
!    Sun, Oct 02 2005 -
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, the first date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2, the second date.
!
  implicit none

  character cmp
  integer ( kind = 4 ) d
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  integer ( kind = 4 ) m
  character ( len = 3 ) m_name
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) w
  character ( len = 3 ) w_name
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y = y1
  m = m1
  d = d1
  f = 0.0D+00

  cmp = '<'

  do while ( cmp /= '>' )

    call ymdf_to_weekday_common ( y, m, d, f, w )

    call weekday_to_name_common3 ( w, w_name )

    call month_to_month_name_common3 ( m, m_name )

    write ( *, '(a3,'','',1x,a3,1x,i2,1x,i4,'' - '')' ) w_name, m_name, d, y

    call ymdf_next_common ( y, m, d, f, y, m, d, f )

    call ymdf_compare ( y, m, d, f, y2, m2, d2, f, cmp )

  end do

  return
end
function days_before_month_common ( y, m )

!*****************************************************************************80
!
!! DAYS_BEFORE_MONTH_COMMON returns the number of days before a Common month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
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
!    Output, integer ( kind = 4 ) DAYS_BEFORE_MONTH_COMMON, the number of
!    days in the year before the first day of the given month.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(12) :: mdays = (/ &
     0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
  integer ( kind = 4 ) days_before_month_common
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
    days_before_month_common = 0
    return
  end if

  days_before_month_common = mdays ( m2 )

  if ( 2 < m2 .and. year_is_leap_common ( y2 ) ) then
    days_before_month_common = days_before_month_common + 1
  end if

  return
end
function days_before_month_gregorian ( y, m )

!*****************************************************************************80
!
!! DAYS_BEFORE_MONTH_GREGORIAN: number of days before a Gregorian month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
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
!    Output, integer ( kind = 4 ) DAYS_BEFORE_MONTH_GREGORIAN, the number of
!    days in the year before the first day of the given month.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(12) :: mdays = (/ &
     0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
  integer ( kind = 4 ) days_before_month_gregorian
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
    days_before_month_gregorian = 0
    return
  end if

  days_before_month_gregorian = mdays ( m2 )

  if ( 2 < m2 .and. year_is_leap_gregorian ( y2 ) ) then
    days_before_month_gregorian = days_before_month_gregorian + 1
  end if

  return
end
function days_before_month_julian ( y, m )

!*****************************************************************************80
!
!! DAYS_BEFORE_MONTH_JULIAN returns the number of days before a Julian month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
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
!    Output, integer ( kind = 4 ) DAYS_BEFORE_MONTH_JULIAN, the number of
!    days in the year before the first day of the given month.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(12) :: mdays = (/ &
     0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
  integer ( kind = 4 ) days_before_month_julian
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
    days_before_month_julian = 0
    return
  end if

  days_before_month_julian = mdays ( m2 )

  if ( 2 < m2 .and. year_is_leap_julian ( y2 ) ) then
    days_before_month_julian = days_before_month_julian + 1
  end if

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
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
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
subroutine deflate_english ( y, m, d )

!*****************************************************************************80
!
!! DEFLATE_ENGLISH "deflates" dates in the English Calendar transition month.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  if ( y == 1752 ) then
    if ( m == 9 ) then
      if ( 14 <= d ) then
        d = d - 11
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
subroutine easter_ds ( y, m, d )

!*****************************************************************************80
!
!! EASTER_DS computes the month and day of Easter for a Gregorian year.
!
!  Example:
!
!    Input:
!
!      Y = 2000
!
!    Output:
!
!      M = 4
!      D = 23
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
!  Reference:
!
!    Peter Duffett-Smith,
!    Practical Astronomy With Your Calculator,
!    Third Edition,
!    Cambridge University Press, 1996,
!    ISBN: 0-521-35699-7,
!    LC: QB62.5.D83.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year, which must be 1583 or greater.
!    (The formula is only valid for years after the Gregorian calendar
!    was adopted.)
!
!    Output, integer ( kind = 4 ) M, D, the month and day of Easter.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) dd
  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  integer ( kind = 4 ) g
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) y

  if ( y <= 0 ) then
    m = -1
    d = -1
    return
  end if

  call year_to_golden_number ( y, a )

  a = a - 1

  b = y / 100
  c = mod ( y, 100 )

  dd = b / 4
  e = mod ( b, 4 )

  f = ( b + 8 ) / 25
  g = ( b - f + 1 ) / 3
  h = mod ( 19 * a + b - dd - g + 15, 30 )

  i = c / 4
  k = mod ( c, 4 )

  l = mod ( 32 + 2 * e + 2 * i - h - k, 7 )
  mm = ( a + 11 * h + 22 * l ) / 451

  m = ( h + l - 7 * mm + 114 ) / 31
  d = mod ( h + l - 7 * mm + 114, 31 ) + 1

  return
end
subroutine easter_egr ( y, m, d )

!*****************************************************************************80
!
!! EASTER_EGR computes the month and day of Easter for a Common year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Algorithm O,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 375.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, integer ( kind = 4 ) M, D, the month and day of Easter.
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) g
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) u
  integer ( kind = 4 ) vp
  integer ( kind = 4 ) y

  if ( y <= 0 ) then
    m = -1
    d = -1
    return
  end if

  p = y + ( y / 4 ) - ( y / 100 ) + ( y / 400 ) - 1
  n = 7 - mod ( p, 7 )
  h = y / 100
  q = h - h / 4
  g = 1 + mod ( y, 19 )
  e = mod ( 57 + 11 * g - q + ( h - ( h - 17 ) / 25 ) / 3, 30 )
  u = mod ( 53 - e, 30 )
  vp = ( g - 1 + 11 * u ) / 319
  r = 22 + u - vp
  c = i4_wrap ( r + 3, 1, 7 )
  s = r + mod ( 7 + n - c, 7 )

  m = 3 + ( s / 32 )
  d = i4_wrap ( s, 1, 31 )

  return
end
subroutine easter_egr2 ( y, m, d )

!*****************************************************************************80
!
!! EASTER_EGR2 computes the month and day of Easter for a Common year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Algorithm P,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 376.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, integer ( kind = 4 ) M, D, the month and day of Easter.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) m
  integer ( kind = 4 ) s
  integer ( kind = 4 ) y

  if ( y <= 0 ) then
    m = -1
    d = -1
    return
  end if

  a = y / 100
  b = a - ( a / 4 )
  c = mod ( y, 19 )
  d = mod ( 15 + 19 * c + b - ( a - ( a - 17 ) / 25 ) / 3, 30 )
  e = d - ( c + 11 * d ) / 319
  s = 22 + e + mod ( 140004 - y - ( y / 4 ) + b - e, 7 )

  m = 3 + ( s / 32 )
  d = i4_wrap ( s, 1, 31 )

  return
end
subroutine easter_julian ( y, m, d )

!*****************************************************************************80
!
!! EASTER_JULIAN computes the date of Easter in the Julian calendar.
!
!  Discussion:
!
!    This computation for the date of Easter uses the Dionysian
!    canon that applied to the Julian calendar.  The determination
!    of the date of Easter changed at the same time that the calendar
!    was modified to use the Gregorian system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Algorithm M,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 365.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, integer ( kind = 4 ) M, D, the month and day of the Julian
!    calendar on which Easter occurs.
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) y

  if ( y <= 0 ) then
    m = -1
    d = -1
    return
  end if

  p = y + ( y / 4 ) + 4
  n = 7 - mod ( p, 7 )

  call year_to_epact_julian ( y, e )

  r = 22 + mod ( 53 - e, 30 )

  c = i4_wrap ( r + 3, 1, 7 )

  s = r + mod ( 7 + n - c, 7 )

  m = 3 + ( s / 32 )
!
!  Use wrapping so that 1 <= D <= 31.
!
  d = i4_wrap ( s, 1, 31 )

  return
end
subroutine easter_julian2 ( y, m, d )

!*****************************************************************************80
!
!! EASTER_JULIAN2 computes the date of Easter in the Julian calendar.
!
!  Discussion:
!
!    This computation for the date of Easter uses the Dionysian
!    canon that applied to the Julian calendar.  The determination
!    of the date of Easter changed at the same time that the calendar
!    was modified to use the Gregorian system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Algorithm N,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 365.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, integer ( kind = 4 ) M, D, the month and day of the Julian calendar
!    on which Easter occurs.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) m
  integer ( kind = 4 ) s
  integer ( kind = 4 ) y

  if ( y <= 0 ) then
    m = -1
    d = -1
    return
  end if

  call year_to_golden_number ( y, a )
  a = a - 1

  b = 22 + mod ( 225 - 11 * a, 30 )
  s = b + mod ( 56 + 6 * y - ( y / 4 ) - b, 7 )

  m = 3 + ( s / 32 )
!
!  Use wrapping to ensure that 1 <= D <= 31.
!
  d = i4_wrap ( s, 1, 31 )

  return
end
subroutine easter_knuth ( y, m, d )

!*****************************************************************************80
!
!! EASTER_KNUTH computes the month and day of Easter for a Gregorian year.
!
!  Discussion:
!
!    Knuth attributes the algorithm to Aloysius Lilius and Christopher Clavius
!    in the late 16th century.  The algorithm is for use with the Gregorian
!    calendar.
!
!  Example:
!
!    Input:
!
!      Y = 2000
!
!    Output:
!
!      M = 4
!      D = 23
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1: Fundamental Algorithms,
!    Addison Wesley, 1968, pages 155-156.
!
!    Donald Knuth,
!    The Calculation of Easter,
!    Communications of the ACM,
!    Volume 5, Number 4, April 1962, pages 209-210.
!
!    Thomas O'Beirne,
!    Puzzles and Paradoxes,
!    Oxford University Press, 1965, chapter 10.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year, which must be 1583 or greater.
!    (The formula is only valid for years after the Gregorian calendar
!    was adopted.)
!
!    Output, integer ( kind = 4 ) M, D, the month and day of Easter.
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) dd
  integer ( kind = 4 ) e
  integer ( kind = 4 ) g
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
  integer ( kind = 4 ) z

  if ( y <= 0 ) then
    m = -1
    d = -1
    return
  end if
!
!  E1: Set the golden number of the year in the 19-year Metonic cycle.
!
  call year_to_golden_number ( y, g )
!
!  E2: Set the century.
!
  c = ( y / 100 ) + 1
!
!  E3: Corrections.
!  X is the number of years divisible by 100 in which leap year was dropped.
!  Z is a special correction to synchronize Easter with the moon's orbit.
!
  x = ( 3 * c / 4 ) - 12
  z = ( 8 * c + 5 ) / 25 - 5
!
!  E4: Find Sunday.
!
  dd = ( 5 * y / 4 ) - x - 10
!
!  E5: Epact
!
  e = i4_modp ( 11 * g + 20 + z - x, 30 )

  if ( ( e == 25 .and. 11 < g ) .or. ( e == 24 ) ) then
    e = e + 1
  end if
!
!  E6: Find the full moon.
!
  n = 44 - e
  if ( n < 21 ) then
    n = n + 30
  end if
!
!  E7: Advance to Sunday.
!
  n = n + 7 - mod ( dd + n, 7 )
!
!  E8: Get month.
!
  if ( 31 < n ) then
    d = n - 31
    m = 4
  else
    d = n
    m = 3
  end if

  return
end
subroutine easter_stewart ( y, m, d )

!*****************************************************************************80
!
!! EASTER_STEWART computes the month and day of Easter for a Gregorian year.
!
!  Example:
!
!    Y = 2001
!
!    A = 6
!    B = 20
!    C = 1
!    DD = 5
!    E = 0
!    G = 6
!    H = 18
!    MM = 0
!    J = 0
!    K = 1
!    L = 6
!    M = 4
!    D = 15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas O'Beirne,
!    Puzzles and Paradoxes,
!    Oxford University Press, 1965.
!
!    Ian Stewart,
!    Easter is a Quasicrystal,
!    Scientific American,
!    March 2001, pages 80-83.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, integer ( kind = 4 ) M, D, the month and day of Easter.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) dd
  integer ( kind = 4 ) e
  integer ( kind = 4 ) g
  integer ( kind = 4 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) y

  a = mod ( y, 19 )
  b = y / 100
  c = mod ( y, 100 )
  dd = b / 4
  e = mod ( b, 4 )
  g = ( 8 * b + 13 ) / 25
  h = mod ( 19 * a + b - dd - g + 15, 30 )
  mm = ( a + 11 * h ) / 319
  j = c / 4
  k = mod ( c, 4 )
  l = mod ( 2 * e + 2 * j - k - h + mm + 32, 7 )

  m = ( h - mm + l + 90 ) / 25
  d = mod ( h - mm + l + m + 19 , 32 )

  return
end
subroutine epoch_to_jed_akbar ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_AKBAR returns the epoch of the Akbar calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 2002
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

  jed = 2289425.5D+00

  return
end
subroutine epoch_to_jed_alexandrian ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_ALEXANDRIAN: epoch of the Alexandrian calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
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

  jed = 1713262.5D+00

  return
end
subroutine epoch_to_jed_armenian ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_ARMENIAN returns the epoch of the Armenian calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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

  jed = 1922867.5D+00

  return
end
subroutine epoch_to_jed_bahai ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_BAHAI returns the epoch of the Bahai calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2002
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

  jed = 2394646.5D+00

  return
end
subroutine epoch_to_jed_bessel ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_BESSEL returns the epoch of the Bessel calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2001
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

  jed = 2415020.31352D+00

  return
end
subroutine epoch_to_jed_chinese ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_CHINESE returns the epoch of the Chinese calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2000
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

  jed = 758325.5D+00

  return
end
subroutine epoch_to_jed_common ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_COMMON returns the epoch of the Common calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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

  jed = 1721423.5D+00

  return
end
subroutine epoch_to_jed_coptic ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_COPTIC returns the epoch of the Coptic calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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

  jed = 1825029.5D+00

  return
end
subroutine epoch_to_jed_deccan ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_DECCAN returns the epoch of the Fasli Deccan calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 2001
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

  jed = 1936747.5D+00

  return
end
subroutine epoch_to_jed_eg_civil ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_EG_CIVIL: epoch of the Egyptian Civil calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 2002
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

  jed = 1448637.5D+00

  return
end
subroutine epoch_to_jed_eg_lunar ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_EG_LUNAR: epoch of the Egyptian Lunar calendar as a JED.
!
!  Discussion:
!
!    This is just a fake value, making the Egyptian Lunar calendar start
!    at the same data as the Egyptian Civil calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2000
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

  jed = 1448637.5D+00

  return
end
subroutine epoch_to_jed_english ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_ENGLISH returns the epoch of the English calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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

  jed = 1721423.5D+00

  return
end
subroutine epoch_to_jed_ethiopian ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_ETHIOPIAN returns the epoch of the Ethiopian calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2002
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

  jed = 1724220.5D+00

  return
end
subroutine epoch_to_jed_gps ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_GPS returns the epoch of the GPS calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2001
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

  jed = 2444244.5D+00

  return
end
subroutine epoch_to_jed_greek ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_GREEK returns the epoch of the Greek calendar as a JED.
!
!  Discussion:
!
!    The Greek Olympiad calendar began on 9 July 776 BC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2000
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

  jed = 1438178.5D+00

  return
end
subroutine epoch_to_jed_gregorian ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_GREGORIAN returns the epoch of the Gregorian calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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

  jed = 1721425.5D+00

  return
end
subroutine epoch_to_jed_hebrew ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_HEBREW returns the epoch of the Hebrew calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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

  jed = 347998.5D+00

  return
end
subroutine epoch_to_jed_hindu_lunar ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_HINDU_LUNAR: epoch of the Hindu lunar calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 2002
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

  jed = 1741959.5D+00

  return
end
subroutine epoch_to_jed_hindu_solar ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_HINDU_SOLAR: epoch of the Hindu solar calendar as a JED.
!
!  Discussion:
!
!    This is the beginning of the Kali Yuga era.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2002
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

  jed = 588465.75D+00

  return
end
subroutine epoch_to_jed_islamic_a ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_ISLAMIC_A returns the epoch of the Islamic A calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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

  jed = 1948438.5D+00

  return
end
subroutine epoch_to_jed_islamic_b ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_ISLAMIC_B returns the epoch of the Islamic B calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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

  jed = 1948439.5D+00

  return
end
subroutine epoch_to_jed_jed ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_JED returns the epoch of the JED as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 2001
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

  jed = 0.0D+00

  return
end
subroutine epoch_to_jed_jelali ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_JELALI returns the epoch of the Jelali calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 November 2000
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

  jed = 2114872.5D+00

  return
end
subroutine epoch_to_jed_julian ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_JULIAN returns the epoch of the Julian calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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

  jed = 1721423.5D+00

  return
end
subroutine epoch_to_jed_khwarizmian ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_KHWARIZMIAN: epoch of the Khwarizmian calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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

  jed = 1952067.5D+00

  return
end
subroutine epoch_to_jed_macedonian ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_MACEDONIAN: epoch of the Macedonian calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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

  jed = 1607708.5D+00

  return
end
subroutine epoch_to_jed_matlab ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_MATLAB: epoch of the "MATLAB calendar" as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2010
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

  jed = 1721058.5D+00

  return
end
subroutine epoch_to_jed_mayan_long ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_MAYAN_LONG: epoch of the Mayan long count calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2002
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

  jed = 584282.5D+00

  return
end
subroutine epoch_to_jed_mjd ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_MJD returns the epoch of the MJD calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2001
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

  jed = 2400000.5D+00

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
subroutine epoch_to_jed_persian ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_PERSIAN returns the epoch of the Persian calendar as a JED.
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
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date of the epoch.
!
  implicit none

  real ( kind = 8 ) jed

  jed = 1952062.5D+00

  return
end
subroutine epoch_to_jed_persian_solar ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_PERSIAN_SOLAR: epoch of the Persian solar calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 2002
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

  jed = 1948320.5D+00

  return
end
subroutine epoch_to_jed_rd ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_RD returns the epoch of the RD calendar as a JED.
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
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date of the epoch.
!
  implicit none

  real ( kind = 8 ) jed

  jed = 1721425.5D+00

  return
end
subroutine epoch_to_jed_republican ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_REPUBLICAN: epoch of the Republican calendar as a JED.
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
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date of the epoch.
!
  implicit none

  real ( kind = 8 ) jed

  jed = 2375839.5D+00

  return
end
subroutine epoch_to_jed_roman ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_ROMAN returns the epoch of the Roman calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!     16 March 2001
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

  jed = 1446389.5D+00

  return
end
subroutine epoch_to_jed_saka ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_SAKA returns the epoch of the Saka calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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

  jed = 1749994.5D+00

  return
end
subroutine epoch_to_jed_soor_san ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_SOOR_SAN: epoch of the Fasli Soor San calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 November 2000
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

  jed = 1940351.5D+00

  return
end
subroutine epoch_to_jed_syrian ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_SYRIAN returns the epoch of the Syrian calendar as a JED.
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
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date of the epoch.
!
  implicit none

  real ( kind = 8 ) jed

  jed = 1607738.5D+00

  return
end
subroutine epoch_to_jed_unix ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_UNIX returns the epoch of the UNIX calendar as a JED.
!
!  Discussion:
!
!    The UNIX Epoch is taken to be the first second of 1 January 1970.
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
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date of the epoch.
!
  implicit none

  real ( kind = 8 ) jed

  jed = 2440587.50D+00

  return
end
subroutine epoch_to_jed_y2k ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_Y2K returns the epoch of the Y2K calendar as a JED.
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
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date of the epoch.
!
  implicit none

  real ( kind = 8 ) jed

  jed = 2451544.5D+00

  return
end
subroutine epoch_to_jed_zoroastrian ( jed )

!*****************************************************************************80
!
!! EPOCH_TO_JED_ZOROASTRIAN: epoch of the Zoroastrian calendar as a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2002
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

  jed = 1862836.5D+00

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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
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
subroutine frac_borrow_english ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_BORROW_ENGLISH borrows fractions from days in an English YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
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

  call day_borrow_english ( y, m, d )

  return
end
subroutine frac_borrow_gregorian ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_BORROW_GREGORIAN borrows fractions from days in a Gregorian YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
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

  call day_borrow_gregorian ( y, m, d )

  return
end
subroutine frac_borrow_hebrew ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_BORROW_HEBREW borrows fractions from days in a Hebrew YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
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

  call day_borrow_hebrew ( y, m, d )

  return
end
subroutine frac_borrow_islamic ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_BORROW_ISLAMIC borrows fractions from days in an Islamic YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
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

  call day_borrow_islamic ( y, m, d )

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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
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
subroutine frac_borrow_republican ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_BORROW_REPUBLICAN borrows fractions from days in a Republican YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
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

  call day_borrow_republican ( y, m, d )

  return
end
subroutine frac_borrow_roman ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_BORROW_ROMAN borrows fractions from days in a Roman YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
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

  call day_borrow_roman ( y, m, d )

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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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
subroutine frac_carry_english ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_CARRY_ENGLISH carries fractions to days in an English YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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

  call day_carry_english ( y, m, d )

  return
end
subroutine frac_carry_gregorian ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_CARRY_GREGORIAN carries fractions from days in a Gregorian YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    a YMDF date.
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

  call day_carry_gregorian ( y, m, d )

  return
end
subroutine frac_carry_hebrew ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_CARRY_HEBREW carries fractions from days in a Hebrew YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    a YMDF date.
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

  call day_carry_hebrew ( y, m, d )

  return
end
subroutine frac_carry_islamic ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_CARRY_ISLAMIC carries fractions from days in an Islamic YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    a YMDF date.
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

  call day_carry_islamic ( y, m, d )

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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    a YMDF date.
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
subroutine frac_carry_republican ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_CARRY_REPUBLICAN carries fractions from days in a Republican YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    a YMDF date.
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

  call day_carry_republican ( y, m, d )

  return
end
subroutine frac_carry_roman ( y, m, d, f )

!*****************************************************************************80
!
!! FRAC_CARRY_ROMAN carries fractions to days in a Roman YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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

  call day_carry_roman ( y, m, d )

  return
end
subroutine frac_to_hms ( f, h, m, s )

!*****************************************************************************80
!
!! FRAC_TO_HMS converts a fractional day into hours, minutes, seconds.
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
!    Input, real ( kind = 8 ) F, a day fraction between 0.0 and 1.0.
!
!    Output, integer ( kind = 4 ) H, M, S, the equivalent hours, minutes
!    and seconds.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) s

  f2 = f

  f2 = 24.0D+00 * f2
  h = int ( f2 )
  f2 = f2 - real ( h, kind = 8 )

  f2 = 60.0D+00 * f2
  m = int ( f2 )
  f2 = f2 - real ( m, kind = 8 )

  f2 = 60.0D+00 * f2
  s = int ( f2 )
  f2 = f2 - real ( s, kind = 8 )

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
subroutine hms_to_s ( h, n, second, s )

!*****************************************************************************80
!
!! HMS_TO_S "prints" an HMS date into a string.
!
!  Format:
!
!    HH:MM:SS
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
!    Input, integer ( kind = 4 ) H, N, SECOND, the HMS date.
!
!    Output, character ( len = * ) S, contains a representation of the date.
!
  implicit none

  integer ( kind = 4 ) h
  integer ( kind = 4 ) n
  character ( len = * ) s
  character ( len = 8 ) s1
  integer ( kind = 4 ) second

  call i4_to_s_zero ( h, s1(1:2) )
  s1(3:3) = ':'
  call i4_to_s_zero ( n, s1(4:5) )
  s1(6:6) = ':'
  call i4_to_s_zero ( second, s1(7:8) )

  s = s1

  return
end
subroutine hour_borrow_common ( y, m, d, h )

!*****************************************************************************80
!
!! HOUR_BORROW_COMMON "borrows" a day of hours.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, H, the year, month, day
!    and hour of the date.  The value of H is presumably negative, and
!    so hours will be "borrowed" to make H positive.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  do while ( h < 0 )

    h = h + 24
    d = d - 1

    call day_borrow_common ( y, m, d )

  end do

  return
end
subroutine hour_carry_common ( y, m, d, h )

!*****************************************************************************80
!
!! HOUR_CARRY_COMMON is given a YMDH date, and carries hours to days.
!
!  Algorithm:
!
!    While 24 <= H:
!
!      decrease H by the number of hours in a day;
!      increase D by 1;
!      if necessary, adjust M and Y.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, H, the year, month, day
!    and hour of the date.  On input, H is presumably 24 or greater.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  do while ( 24 <= h )

    h = h - 24
    d = d + 1

    call day_carry_common ( y, m, d )

  end do

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
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two integers.
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
!    Input/output, integer ( kind = 4 ) I, J, the two integers to be swapped.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
function i4_to_a ( i )

!*****************************************************************************80
!
!! I4_TO_A returns the I-th alphabetic character.
!
!  Example:
!
!    I  I4_TO_A
!
!   -8  ' '
!    0  ' '
!    1  'A'
!    2  'B'
!   ..
!   26  'Z'
!   27  'a'
!   52  'z'
!   53  ' '
!   99  ' '
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the letter to be returned.
!    0 is a space;
!    1 through 26 requests 'A' through 'Z', (ASCII 65:90);
!    27 through 52 requests 'a' through 'z', (ASCII 97:122);
!
!    Output, character I4_TO_A, the requested alphabetic letter.
!
  implicit none

  integer ( kind = 4 ), parameter :: cap_shift = 64
  integer ( kind = 4 ) i
  character i4_to_a
  integer ( kind = 4 ), parameter :: low_shift = 96

  if ( i <= 0 ) then
    i4_to_a = ' '
  else if ( 1 <= i .and. i <= 26 ) then
    i4_to_a = char ( cap_shift + i )
  else if ( 27 <= i .and. i <= 52 ) then
    i4_to_a = char ( low_shift + i - 26 )
  else if ( 53 <= i ) then
    i4_to_a = ' '
  end if

  return
end
subroutine i4_to_roman ( intval, s )

!*****************************************************************************80
!
!! I4_TO_ROMAN converts an integer to a string of Roman numerals.
!
!  Example:
!
!    INTVAL  S
!
!        -2  -II <-- Not a Roman numeral
!        -1  -I  <-- Not a Roman numeral
!         0   0  <-- Not a Roman numeral
!         1   I
!         2   II
!         3   III
!         4   IV
!         5   V
!        10   X
!        20   XX
!        30   XXX
!        40   XL
!        50   L
!        60   LX
!        70   LXX
!        80   LXXX
!        90   XC
!       100   C
!       500   D
!      1000   M
!      4999   MMMMCMLXLIX
!
!  Discussion:
!
!    To generate numbers greater than 4999, the numeral 'V' had a bar
!    above it, representing a value of 5000, a barred 'X' represented
!    10,000 and so on.
!
!    In the subtractive representation of 4 by 'IV', 9 by 'IX' and so on,
!    'I' can only subtract from 'V' or 'X',
!    'X' can only subtract from 'L' or 'C',
!    'C' can only subtract from 'D' or 'M'.
!    Under these rules, 1999 cannot be written IMM!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.  If the
!    integer has absolute value greater than 4999, the string '?' will be
!    returned.  If the integer is 0, then the string '0' will be returned.  If
!    the integer is negative, then a minus sign will precede it, even
!    though this has nothing to do with Roman numerals.
!
!    Output, character ( len = * ) S, the representation of the integer
!    as a Roman numeral.
!
  implicit none

  integer ( kind = 4 ) icopy
  integer ( kind = 4 ) intval
  character ( len = * ) s

  s = ' '
  icopy = intval

  if ( 4999 < abs ( icopy ) ) then
    s = '?'
    return
  end if

  if ( icopy == 0 ) then
    s = '0'
    return
  end if

  if ( icopy <= 0 ) then
    s = '-'
    icopy = - icopy
  end if

  do while ( 0 < icopy )

    if ( 1000 <= icopy ) then
      call s_cat ( s, 'M', s )
      icopy = icopy - 1000
    else if ( 900 <= icopy ) then
      call s_cat ( s, 'CM', s )
      icopy = icopy - 900
    else if ( 500 <= icopy ) then
      call s_cat ( s, 'D', s )
      icopy = icopy - 500
    else if ( 400 <= icopy ) then
      call s_cat ( s, 'CD', s )
      icopy = icopy - 400
    else if ( 100 <= icopy ) then
      call s_cat ( s, 'C', s )
      icopy = icopy - 100
    else if ( 90 <= icopy ) then
      call s_cat ( s, 'XC', s )
      icopy = icopy - 90
    else if ( 50 <= icopy ) then
      call s_cat ( s, 'L', s )
      icopy = icopy - 50
    else if ( 40 <= icopy ) then
      call s_cat ( s, 'XL', s )
      icopy = icopy - 40
    else if ( 10 <= icopy ) then
      call s_cat ( s, 'X', s )
      icopy = icopy - 10
    else if ( 9 <= icopy ) then
      call s_cat ( s, 'IX', s )
      icopy = icopy - 9
    else if ( 5 <= icopy ) then
      call s_cat ( s, 'V', s )
      icopy = icopy - 5
    else if ( 4 <= icopy ) then
      call s_cat ( s, 'IV', s )
      icopy = icopy - 4
    else
      call s_cat ( s, 'I', s )
      icopy = icopy - 1
    end if

  end do

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
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, the YMD date.
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
subroutine inflate_english ( y, m, d )

!*****************************************************************************80
!
!! INFLATE_ENGLISH "inflates" dates in the English Calendar transition month.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) M,
!    integer ( kind = 4 ) D, the YMD date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  if ( y == 1752 ) then
    if ( m == 9 ) then
      if ( 3 <= d ) then
        d = d + 11
      end if
    end if
  end if

  return
end
subroutine j_borrow_common ( y, j )

!*****************************************************************************80
!
!! J_BORROW_COMMON borrows year-days from years in a Common date.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_common

  do while ( j <= 0 )

    y = y - 1

    days = year_length_common ( y )

    j = j + days

  end do

  return
end
subroutine j_borrow_english ( y, j )

!*****************************************************************************80
!
!! J_BORROW_ENGLISH borrows year-days from years in an English date.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_english

  do while ( j <= 0 )

    y = y - 1

    days = year_length_english ( y )

    j = j + days

  end do

  return
end
subroutine j_borrow_gregorian ( y, j )

!*****************************************************************************80
!
!! J_BORROW_GREGORIAN borrows year-days from years in a Gregorian date.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_gregorian

  do while ( j <= 0 )

    y = y - 1

    days = year_length_gregorian ( y )

    j = j + days

  end do

  return
end
subroutine j_borrow_hebrew ( y, j )

!*****************************************************************************80
!
!! J_BORROW_HEBREW borrows year-days from years in a Hebrew date.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_hebrew

  do while ( j <= 0 )

    y = y - 1

    days = year_length_hebrew ( y )

    j = j + days

  end do

  return
end
subroutine j_borrow_islamic ( y, j )

!*****************************************************************************80
!
!! J_BORROW_ISLAMIC borrows year-days from years in an Islamic date.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_islamic

  do while ( j <= 0 )

    y = y - 1

    days = year_length_islamic ( y )

    j = j + days

  end do

  return
end
subroutine j_borrow_julian ( y, j )

!*****************************************************************************80
!
!! J_BORROW_JULIAN borrows year-days from years in a Julian date.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_julian

  do while ( j <= 0 )

    y = y - 1

    days = year_length_julian ( y )

    j = j + days

  end do

  return
end
subroutine j_borrow_republican ( y, j )

!*****************************************************************************80
!
!! J_BORROW_REPUBLICAN borrows year-days from years in a Republican date.
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
!    Input/output, integer ( kind = 4 ) Y, integer ( kind = 4 ) J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_republican

  do while ( j <= 0 )

    y = y - 1

    days = year_length_republican ( y )

    j = j + days

  end do

  return
end
subroutine j_borrow_roman ( y, j )

!*****************************************************************************80
!
!! J_BORROW_ROMAN borrows year-days from years in a Roman date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_roman

  do while ( j <= 0 )

    y = y - 1

    days = year_length_roman ( y )

    j = j + days

  end do

  return
end
subroutine j_carry_common ( y, j )

!*****************************************************************************80
!
!! J_CARRY_COMMON carries year-days to years in a Common date.
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
!    Input/output, integer ( kind = 4 ) Y, J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_common

  do

    days = year_length_common ( y )

    if ( j < days ) then
      exit
    end if

    j = j - days
    y = y + 1

  end do

  return
end
subroutine j_carry_english ( y, j )

!*****************************************************************************80
!
!! J_CARRY_ENGLISH carries year-days to years in an English date.
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
!    Input/output, integer ( kind = 4 ) Y, J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_english

  do

    days = year_length_english ( y )

    if ( j < days ) then
      exit
    end if

    j = j - days
    y = y + 1

  end do

  return
end
subroutine j_carry_gregorian ( y, j )

!*****************************************************************************80
!
!! J_CARRY_GREGORIAN carries year-days to years in a Gregorian date.
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
!    Input/output, integer ( kind = 4 ) Y, J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_gregorian

  do

    days = year_length_gregorian ( y )

    if ( j < days ) then
      exit
    end if

    j = j - days
    y = y + 1

  end do

  return
end
subroutine j_carry_hebrew ( y, j )

!*****************************************************************************80
!
!! J_CARRY_HEBREW carries year-days to years in a Hebrew date.
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
!    Input/output, integer ( kind = 4 ) Y, J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_hebrew

  do

    days = year_length_hebrew ( y )

    if ( j < days ) then
      exit
    end if

    j = j - days
    y = y + 1

  end do

  return
end
subroutine j_carry_islamic ( y, j )

!*****************************************************************************80
!
!! J_CARRY_ISLAMIC carries year-days to years in an Islamic date.
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
!    Input/output, integer ( kind = 4 ) Y, J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_islamic

  do

    days = year_length_islamic ( y )

    if ( j < days ) then
      exit
    end if

    j = j - days
    y = y + 1

  end do

  return
end
subroutine j_carry_julian ( y, j )

!*****************************************************************************80
!
!! J_CARRY_JULIAN carries year-days to years in a Julian date.
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
!    Input/output, integer ( kind = 4 ) Y, J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_julian

  do

    days = year_length_julian ( y )

    if ( j < days ) then
      exit
    end if

    j = j - days
    y = y + 1

  end do

  return
end
subroutine j_carry_republican ( y, j )

!*****************************************************************************80
!
!! J_CARRY_REPUBLICAN carries year-days to years in a Republican date.
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
!    Input/output, integer ( kind = 4 ) Y, J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_republican

  do

    days = year_length_republican ( y )

    if ( j < days ) then
      exit
    end if

    j = j - days
    y = y + 1

  end do

  return
end
subroutine j_carry_roman ( y, j )

!*****************************************************************************80
!
!! J_CARRY_ROMAN carries year-days to years in a Roman date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, J, a YJ date.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_roman

  do

    days = year_length_roman ( y )

    if ( j < days ) then
      exit
    end if

    j = j - days
    y = y + 1

  end do

  return
end
subroutine jed_check ( jed, ierror )

!*****************************************************************************80
!
!! JED_CHECK checks a Julian Ephemeris Date.
!
!  Discussion:
!
!    The routine returns an error if JED < 0, although there is no
!    reason why such dates are invalid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if JED is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  real ( kind = 8 ) jed

  if ( 0.0D+00 <= jed ) then
    ierror = 0
  else
    ierror = 1
  end if

  return
end
subroutine jed_test ( i, jed )

!*****************************************************************************80
!
!! JED_TEST returns some "interesting" JED's.
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
!    Bonnie Blackburn, Leofranc Holford-Stevens,
!    The Oxford Companion to the Year,
!    Oxford, 1999.
!
!    Frank Parise, editor,
!    The Book of Calendars,
!    Facts on File, Inc, 1982,
!    CE11.K4 / 529.3.
!
!    Edward Reingold, Nachum Dershowitz,
!    Calendrical Calculations, the Millennium Edition,
!    Cambridge, 2002,
!    CE12.R45 / 529.3-dc21
!
!    Edward Richards,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the test date requested.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!    If I is less than 1, or greater than the number of test dates
!    available, JED is returned as -1.0.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch_50000
!
!  JED Epoch:
!  Beginning of current Scaliger cycle.
!  Monday, Noon, 1 January 4713 BCE/Julian
!
  if ( i == 1 ) then

    jed = 0.0D+00
!
!  The day after the JED Epoch.
!  Tuesday, Noon, 2 January 4713 BCE/Julian
!
  else if ( i == 2 ) then

    jed = 1.0D+00
!
!  Archbishop James Ussher's estimate of the date of Creation,
!  (Noon), 23 October 4004 BCE/Julian
!
  else if ( i == 3 ) then

    jed = 259258.000D+00
!
!  Hebrew Epoch.
!  7 October 3761 BCE/Julian
!
  else if ( i == 4 ) then

    jed = 347998.5D+00
!
!  Mayan Long Count Epoch.
!  6 September 3114 BCE/Julian
!  (Reingold and Dershowitz)
!
  else if ( i == 5 ) then

    jed = 584282.5D+00
!
!  Hindu Solar Epoch.
!  Beginning of the Kali Yuga age.
!  18 February 3102 BCE/Julian
!
  else if ( i == 6 ) then

    jed = 588465.75D+00
!
!  Chinese Epoch.
!  8 March 2637 BCE/Julian
!
  else if ( i == 7 ) then

    jed = 758325.5D+00
!
!  Greek Olympiad Epoch
!  9 July 776 BCE/Julian
!
  else if ( i == 8 ) then

    jed = 1438178.5D+00
!
!  Roman Epoch
!  Ab Urbe Condita
!  1 January 753 BCE/Julian
!
  else if ( i == 9 ) then

    jed = 1446389.5D+00
!
!  Egyptian Civil Calendar Epoch.
!  Ascension of Nabonassar to throne of Babylon.
!  26 February 747 BCE/Julian
!
  else if ( i == 10 ) then

    jed = 1448637.5D+00
!
!  Egyptian Lunar Calendar Epoch.
!  (Don't really know where to set this...)
!  Ascension of Nabonassar to throne of Babylon.
!  26 February 747 BCE/Julian
!
  else if ( i == 11 ) then

    jed = 1448637.5D+00
!
!  Macedonian Epoch
!  1 September 312 BCE/Julian
!
  else if ( i == 12 ) then

    jed = 1607708.5D+00
!
!  Syrian Epoch
!  1 October 312 BCE/Julian
!
  else if ( i == 13 ) then

    jed = 1607738.5D+00
!
!  Alexandrian Epoch
!  29 August 23 BCE/Julian
!
  else if ( i == 14 ) then

    jed = 1713262.5D+00
!
!  "1 January, 0 BC"?  MATLAB epoch?
!
  else if ( i == 15 ) then

    jed = 1721058.5D+00
!
!  Julian Epoch MINUS ONE DAY
!  Friday, 31 December 1 BCE/Julian
!
  else if ( i == 16 ) then

    jed = 1721423.5D+00
    jed = jed - 1.0D+00
!
!  Julian Epoch
!  Saturday, 1 January 1 CE/Julian
!
  else if ( i == 17 ) then

    jed = 1721423.5D+00
!
!  Gregorian Epoch
!  Monday, 3 January 1 CE/Julian
!  Monday, 1 January 1 Gregorian
!
  else if ( i == 18 ) then

    jed = 1721425.5D+00
!
!  RD: Reingold and Dershowitz Epoch
!  Monday, 3 January 1 CE/Julian
!  Monday, 1 January 1 Gregorian
!
  else if ( i == 19 ) then

    jed = 1721425.5D+00
!
!  Ethiopian Epoch
!  29 August 8 CE/Julian
!  (Reingold and Dershowitz)
!
  else if ( i == 20 ) then

    jed = 1724220.5D+00
!
!  Hindu Lunar Epoch, the Vikrama
!  24 March 57 CE/Julian
!  (The actual day and month are not specified by RD)
!  (Reingold and Dershowitz)
!
  else if ( i == 21 ) then

    jed = 1741959.5D+00
!
!  Saka Epoch
!  4 March 79 CE/Julian
!
  else if ( i == 22 ) then

    jed = 1749994.5D+00
!
!  Coptic Epoch
!  29 August 284 CE/Julian
!
  else if ( i == 23 ) then

    jed = 1825029.5D+00
!
!  Zoroastrian Epoch.
!  3 March 388 CE/Julian
!
   else if ( i == 24 ) then

     jed = 1862836.5D+00
!
!  Armenian Epoch
!  11 July 552 CE/Julian
!
  else if ( i == 25 ) then

    jed = 1922867.5D+00
!
!  Fasli Deccan Epoch
!  12 July 590 CE/Julian
!
   else if ( i == 26 ) then

     jed = 1936747.5D+00
!
!  Fasli Soor San Epoch
!  24 May 600 CE/Julian
!
   else if ( i == 27 ) then

     jed = 1940351.5D+00
!
!  Persian Solar Epoch
!  19 March 622 CE/Julian
!
  else if ( i == 28 ) then

    jed = 1948320.5D+00
!
!  Islamic A Epoch
!  Thursday, 15 July 622 CE/Julian
!
  else if ( i == 29 ) then

    jed = 1948438.5D+00
!
!  Islamic B Epoch
!  Friday, 16 July 622 CE/Julian
!
  else if ( i == 30 ) then

    jed = 1948439.5D+00
!
!  Yazdegerd Epoch
!  16 June 632 CE
!
  else if ( i == 31 ) then

    jed = 1952062.5D+00
!
!  Khwarizmian Epoch
!  21 June 632 CE/Julian
!
  else if ( i == 32 ) then

    jed = 1952067.5D+00
!
!  Battle of Hastings.
!  Saturday, 14 October 1066 CE/Julian.
!           (20 October 1066 Gregorian.)
!
  else if ( i == 33 ) then

    jed = 2110700.5D+00
!
!  Jelali Epoch
!  17 March 1078 CE/Julian
!
  else if ( i == 34 ) then

    jed = 2114872.5D+00
!
!  Akbar Epoch
!  9 February 1556 CE/Julian
!  19 February 1556 Gregorian
!
  else if ( i == 35 ) then

    jed = 2289425.5D+00
!
!  Common Era calendar transition:
!  Noon of the last day of Julian calendar usage.
!  Thursday, 04 October 1582 CE/English/Julian
!  Thursday, 14 October 1582 Gregorian
!
  else if ( i == 36 ) then

    jed = 2299160.5D+00
    jed = jed - 0.5D+00
!
!  Common Era calendar transition:
!  Noon of the first day of Gregorian calendar usage.
!  Friday, 05 October 1582 English/Julian
!  Friday, 15 October 1582 CE/Gregorian
!
  else if ( i == 37 ) then

    jed = 2299160.5D+00
    jed = jed + 0.5D+00
!
!  A day chosen by Lewis Carroll to test his day-of-the-week algorithm,
!  Wednesday, 4 March 1676 CE/Gregorian
!  Wednesday, 23 February 1676 English/Julian
!
  else if ( i == 38 ) then

    jed = 2333269.5D+00
!
!  English calendar
!  noon of the last day of Julian calendar usage.
!  02 September 1752 English/Julian
!  13 September 1752 CE/Gregorian
!
  else if ( i == 39 ) then

    jed = 2361221.5D+00
    jed = jed - 0.5D+00
!
!  English calendar,
!  noon of the first day of Gregorian calendar usage.
!  03 September 1752 Julian
!  14 September 1752 CE/English/Gregorian
!
  else if ( i == 40 ) then

    jed = 2361221.5D+00
    jed = jed + 0.5D+00
!
!  A day chosen by Lewis Carroll to test his day-of-the-week algorithm,
!  Thursday, 18 September 1783 CE/Gregorian
!
  else if ( i == 41 ) then

    jed = 2372547.5D+00
!
!  French Republican Epoch
!  Saturday, 11 September 1792 Julian
!  Saturday, 22 September 1792 CE/Gregorian
!
  else if ( i == 42 ) then

    jed = 2375839.5D+00
!
!  Bahai Epoch.
!  9 March 1844 Julian
!  21 March 1844 CE/Gregorian
!
  else if ( i == 43 ) then

    jed = 2394646.5D+00
!
!  Clive James Lucas test date.
!
  else if ( i == 44 ) then

    jed = 2394710.50D+00
!
!  New York Times "epoch" date,
!  fictitious Volume 1, issue #0,
!  17 September 1851
!  (issue #1 was on 18 September 1851):
!
  else if ( i == 45 ) then

    jed = 2397383.50D+00
!
!  Modified Julian Date Epoch.
!  17 November 1858 CE/Gregorian
!
  else if ( i == 46 ) then

    jed = 2400000.5D+00
!
!  NYT issue 10,000
!  24 September 1883
!
  else if ( i == 47 ) then

    jed_epoch_50000 = 2449790.5D+00
    jed = jed_epoch_50000 - 40000.0D+00 - 88.0D+00
!
!  Bessel Year Count Epoch.
!  1 January 1900 CE/Gregorian
!
  else if ( i == 48 ) then

    jed = 2415020.31352D+00
!
!  NYT issue 30,000
!  14 March 1940
!
  else if ( i == 49 ) then

    jed_epoch_50000 = 2449790.5D+00
    jed = jed_epoch_50000 - 20000.0D+00 - 88.0D+00
!
!  NYT issue 40,000
!  ???
!
  else if ( i == 50 ) then

    jed_epoch_50000 = 2449790.5D+00
    jed = jed_epoch_50000 - 10000.0D+00 - 88.0D+00
!
!  UNIX epoch.
!  1 January 1970 CE/Gregorian.
!
  else if ( i == 51 ) then

    jed = 2440587.50D+00
!
!  NYT issue 44027
!  ???
!
  else if ( i == 52 ) then

    jed_epoch_50000 = 2449790.5D+00
    jed = jed_epoch_50000 - 5973
!
!  NYT issue 44028
!  ???
!
  else if ( i == 53 ) then

    jed_epoch_50000 = 2449790.5D+00
    jed = jed_epoch_50000 - 5972
!
!  GPS epoch.
!  6 January 1980 CE/Gregorian
!
  else if ( i == 54 ) then

    jed = 2444244.5D+00
!
!  NYT issue 50,000
!  14 March 1995
!
  else if ( i == 55 ) then

    jed_epoch_50000 = 2449790.5D+00
    jed = jed_epoch_50000
!
!  25 February 1996
!  A Reingold/Dershowitz test date.
!
  else if ( i == 56 ) then

    jed = 2450138.5D+00
!
!  Y2K day
!  1 January 2000 CE/Gregorian
!
  else if ( i == 57 ) then

    jed = 2451544.5D+00
!
!  Today
!
  else if ( i == 58 ) then

    call now_to_jed ( jed )
!
!  End of Current Mayan Great Cycle
!  23 December 2012 CE/Gregorian
!
  else if ( i == 59 ) then

    jed = 2456284.5D+00
!
!  Scaliger cycle repeats.
!  1 January 3266 CE/Gregorian
!
  else if ( i == 60 ) then

    jed = 2913943.0D+00

  else

    jed = -1.0D+00

  end if

  return
end
subroutine jed_to_cws_gps ( jed, c, w, s )

!*****************************************************************************80
!
!! JED_TO_CWS_GPS converts a JED to a GPS CWS date.
!
!  Discussion:
!
!    The GPS time keeping is in terms of seconds, weeks, and cycles
!    of 1024 weeks.  The weeks and cycles begin numbering at 0.
!
!    The computation is only valid for dates after the GPS epoch,
!    that is, after 6 January 1980.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) C, W, real ( kind = 8 ) S, the
!    corresponding GPS cycles/weeks/seconds date.
!
  implicit none

  integer ( kind = 4 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) w
  real ( kind = 8 ) s

  call epoch_to_jed_gps ( jed_epoch )

  d = jed - jed_epoch

  if ( d < 0.0D+00 ) then
    s = -1.0
    w = -1
    c = -1
    return
  end if

  w = int ( d ) / 7
  d = d - real ( 7 * w, kind = 8 )

  c = w / 1024
  w = w - 1024 * c

  s = d * real ( 24.0D+00 * 60.0D+00 * 60.0D+00, kind = 8 )

  return
end
subroutine jed_to_datenum ( jed, dn )

!*****************************************************************************80
!
!! JED_TO_DATENUM converts a JED to a MATLAB date number.
!
!  Discussion:
!
!    The MATLAB "datenum" function accepts a string defining
!    a date and returns a datenumber:
!
!      dn = datenum ( 'Aug 17 1939' )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, real ( kind = 8 ) DN, a MATLAB date number.
!
  implicit none

  real ( kind = 8 ) dn
  real ( kind = 8 ) jed

  dn = jed - 1721058.5D+00

  return
end
subroutine jed_to_mayan_long ( jed, pictun, baktun, katun, tun, uinal, kin, f )

!*****************************************************************************80
!
!! JED_TO_MAYAN_LONG converts a JED to a Mayan long count date.
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
!  Reference:
!
!    Edward Richards,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, chapter 27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) PICTUN, BAKTUN, KATUN, TUN, UINAL, KIN, values
!    defining the Mayan long date.
!
!    Output, real ( kind = 8 ) F, the fractional part of the date.
!
  implicit none

  integer ( kind = 4 ) baktun
  integer ( kind = 4 ) days
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) katun
  integer ( kind = 4 ) kin
  integer ( kind = 4 ) pictun
  integer ( kind = 4 ) tun
  integer ( kind = 4 ) uinal

  call epoch_to_jed_mayan_long ( jed_epoch )

  j = int ( jed - jed_epoch )
  f = ( jed - jed_epoch ) - real ( j, kind = 8 )

  days = j

  if ( 0 <= days ) then

    pictun = days / 2880000
    days = days - pictun * 2880000

  else

    pictun = 0
    do while ( days < 0 )
      pictun = pictun - 1
      days = days + 2880000
    end do

  end if

  baktun = days / 144000
  days = days - baktun * 144000
  katun = days / 7200
  days = days - katun * 7200
  tun = days / 360
  days = days - tun * 360
  uinal = days / 20
  days = days - uinal * 20
  kin = days

  return
end
subroutine jed_to_mayan_round ( jed, y, a, b, c, d, f )

!*****************************************************************************80
!
!! JED_TO_MAYAN_ROUND converts a JED to a Mayan round date.
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
!  Reference:
!
!    Edward Richards,
!    Algorithm K,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 340.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, A, B, C, D, values defining the Mayan
!    round date.
!
!    Output, real ( kind = 8 ) F, the fractional part of the date.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  real ( kind = 8 ) f
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) n
  integer ( kind = 4 ) y

  call epoch_to_jed_mayan_long ( jed_epoch )

  j = int ( jed - jed_epoch )
  f = ( jed - jed_epoch ) - real ( j, kind = 8 )

  days = j

  y = 0

  do while ( days < 0 )
    days = days + 18980
    y = y - 1
  end do

  y = y + days / 18980

  days = mod ( days, 18980 )

  a = i4_wrap ( days + 4, 1, 13 )
  b = i4_wrap ( days, 1, 20 )

  n = mod ( days + 348, 365 )
  c = mod ( n, 20 )
  d = n / 20

  return
end
subroutine jed_to_mjd ( jed, mjd )

!*****************************************************************************80
!
!! JED_TO_MJD converts a JED to a modified JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, real ( kind = 8 ) MJD, the modified Julian Ephemeris Date.
!
  implicit none

  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) mjd

  call epoch_to_jed_mjd ( jed_epoch )

  mjd = jed - jed_epoch

  return
end
subroutine jed_to_nearest_noon ( jed1, jed2 )

!*****************************************************************************80
!
!! JED_TO_NEAREST_NOON converts a JED to the JED of the nearest noon.
!
!  Discussion:
!
!    This is primarily to make a fair test of the weekday routines,
!    which have trouble when the JED is at midnight.
!
!    Note that noon corresponds to an integral JED value.
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
!  Parameters:
!
!    Input, real ( kind = 8 ) JED1, the Julian Ephemeris Date.
!
!    Output, real ( kind = 8 ) JED2, the Julian Ephemeris Date
!    of the nearest noon.  If JED1 was at midnight, JED2 is
!    advanced to the NEXT noon, not the previous one.
!
  implicit none

  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2

  jed2 = anint ( jed1 )

  return
end
subroutine jed_to_next_noon ( jed1, jed2 )

!*****************************************************************************80
!
!! JED_TO_NEXT_NOON converts a JED to the JED of the next noon.
!
!  Discussion:
!
!    This is primarily to make a fair test of the weekday routines,
!    which have trouble when the JED is at midnight.
!
!    Note that noon corresponds to an integral JED value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED1, the Julian Ephemeris Date.
!
!    Output, real ( kind = 8 ) JED2, the Julian Ephemeris Date
!    of the next noon.
!
  implicit none

  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2

  jed2 = aint ( jed1 )
!
!  The integer part of JED1 is one of the two integers that
!  bracket JED1.  If it's the smaller one (which it should
!  be as long as JED1 is positive), make it the bigger one.
!
!  This correctly leaves undisturbed cases where JED1 is
!  already an integer, and where JED1 is negative (which
!  is not a case we expect to occur often).
!
  if ( jed2 < jed1 ) then
    jed2 = jed2 + 1.0D+00
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
!    paper was printed.  This was the first time that a strike caused the issue
!    number sequence to halt.
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
  sundays = ( days + 3 ) / 7
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
    issue = issue - int ( min ( jed_05_11_1978, jed ) - jed_10_08_1978 ) - 1
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
subroutine jed_to_rd ( jed, rd )

!*****************************************************************************80
!
!! JED_TO_RD converts a JED to an RD.
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
!  Reference:
!
!    Edward Reingold, Nachum Dershowitz,
!    Calendrical Calculations, the Millennium Edition,
!    Cambridge, 2002.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, real ( kind = 8 ) RD, the RD date.
!
  implicit none

  real ( kind = 8 ) jed
  real ( kind = 8 ) rd
  real ( kind = 8 ) rd_epoch

  call epoch_to_jed_rd ( rd_epoch )

  rd = jed - rd_epoch

  return
end
subroutine jed_to_ss_unix ( jed, s )

!*****************************************************************************80
!
!! JED_TO_SS_UNIX converts a JED to a UNIX SS date.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, real ( kind = 8 ) S, the corresponding UNIX SS date.
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) s

  call epoch_to_jed_unix ( jed_epoch )

  d = jed - jed_epoch

  s = d * 24.0D+00 * 60.0D+00 * 60.0D+00

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
subroutine jed_to_year_hebrew ( jed, y )

!*****************************************************************************80
!
!! JED_TO_YEAR_HEBREW: the year in the Hebrew calendar when a JED occurred.
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
!    Algorithm H,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 331.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, the year in the Hebrew calendar that
!    included the JED.  If the input JED is less than the epoch of the Hebrew
!    calendar, then Y is always returned as -1.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed2
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  call epoch_to_jed_hebrew ( jed_epoch )

  if ( jed < jed_epoch ) then
    y = -1
    return
  end if
!
!  Using integer arithmetic in this computation may cause overflow.
!
!  Compute the number of months elapsed up to the date.
!
  m = 1 + int ( ( 25920.0D+00 * ( jed - jed_epoch + 2.5D+00 ) ) / 765433.0D+00 )
!
!  Estimate the number of years represented by these months.
!
  y = 19 * ( m / 235 ) + ( 19 * ( i4_modp ( m, 235 ) - 2 ) ) / 235 + 1
  y = max ( y, 1 )
!
!  Determine the JED of the first day of that year.
!
  call new_year_to_jed_hebrew ( y, jed2 )
!
!  We might have been off by 1 year.
!
  if ( jed < jed2 ) then
    y = y - 1
  end if

  return
end
subroutine jed_to_yearcount_bessel ( jed, bessel )

!*****************************************************************************80
!
!! JED_TO_YEARCOUNT_BESSEL converts a JED to Bessel year count.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, real ( kind = 8 ) BESSEL, the Bessel year.
!
  implicit none

  real ( kind = 8 ) bessel
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ), parameter :: year_length = 365.242198781D+00

  call epoch_to_jed_bessel ( jed_epoch )
  bessel = 1900.0D+00 + ( jed - jed_epoch ) / year_length

  return
end
subroutine jed_to_yearcount_julian ( jed, julian )

!*****************************************************************************80
!
!! JED_TO_YEARCOUNT_JULIAN converts a JED to a Julian year count.
!
!  Discussion:
!
!    An average year in the Julian calendar is exactly 365.25 days long.
!    This calculation counts the number of average Julian years from
!    the beginning of the year 2000.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, real ( kind = 8 ) JULIAN, the Julian year.
!
  implicit none

  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) julian
  real ( kind = 8 ), parameter :: year_length = 365.25D+00

  call epoch_to_jed_y2k ( jed_epoch )

  julian = 2000.0D+00 + ( jed - jed_epoch ) / year_length

  return
end
subroutine jed_to_yjf_common ( jed, y, j, f )

!*****************************************************************************80
!
!! JED_TO_YJF_COMMON converts a JED to a Common YJF date.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call jed_to_ymdf_common ( jed, y1, m1, d1, f1 )

  call ymdf_to_yjf_common ( y1, m1, d1, f1, y, j, f )

  return
end
subroutine jed_to_yjf_english ( jed, y, j, f )

!*****************************************************************************80
!
!! JED_TO_YJF_ENGLISH converts a JED to an English YJF date.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call jed_to_ymdf_english ( jed, y1, m1, d1, f1 )

  call ymdf_to_yjf_english ( y1, m1, d1, f1, y, j, f )

  return
end
subroutine jed_to_yjf_gregorian ( jed, y, j, f )

!*****************************************************************************80
!
!! JED_TO_YJF_GREGORIAN converts a JED to a Gregorian YJF date.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call jed_to_ymdf_gregorian ( jed, y1, m1, d1, f1 )

  call ymdf_to_yjf_gregorian ( y1, m1, d1, f1, y, j, f )

  return
end
subroutine jed_to_yjf_hebrew ( jed, y, j, f )

!*****************************************************************************80
!
!! JED_TO_YJF_HEBREW converts a JED to a Hebrew YJF date.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call jed_to_ymdf_hebrew ( jed, y1, m1, d1, f1 )

  call ymdf_to_yjf_hebrew ( y1, m1, d1, f1, y, j, f )

  return
end
subroutine jed_to_yjf_islamic_a ( jed, y, j, f )

!*****************************************************************************80
!
!! JED_TO_YJF_ISLAMIC_A converts a JED to an Islamic-A YJF date.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call jed_to_ymdf_islamic_a ( jed, y1, m1, d1, f1 )

  call ymdf_to_yjf_islamic ( y1, m1, d1, f1, y, j, f )

  return
end
subroutine jed_to_yjf_islamic_b ( jed, y, j, f )

!*****************************************************************************80
!
!! JED_TO_YJF_ISLAMIC_B converts a JED to an Islamic-B YJF date.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call jed_to_ymdf_islamic_b ( jed, y1, m1, d1, f1 )

  call ymdf_to_yjf_islamic ( y1, m1, d1, f1, y, j, f )

  return
end
subroutine jed_to_yjf_julian ( jed, y, j, f )

!*****************************************************************************80
!
!! JED_TO_YJF_JULIAN converts a JED to a Julian YJF date.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call jed_to_ymdf_julian ( jed, y1, m1, d1, f1 )

  call ymdf_to_yjf_julian ( y1, m1, d1, f1, y, j, f )

  return
end
subroutine jed_to_yjf_republican ( jed, y, j, f )

!*****************************************************************************80
!
!! JED_TO_YJF_REPUBLICAN converts a JED to a Republican YJF date.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call jed_to_ymdf_republican ( jed, y1, m1, d1, f1 )

  call ymdf_to_yjf_republican ( y1, m1, d1, f1, y, j, f )

  return
end
subroutine jed_to_yjf_roman ( jed, y, j, f )

!*****************************************************************************80
!
!! JED_TO_YJF_ROMAN converts a JED to a Roman YJF date.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call jed_to_ymdf_roman ( jed, y1, m1, d1, f1 )

  call ymdf_to_yjf_roman ( y1, m1, d1, f1, y, j, f )

  return
end
subroutine jed_to_ymdf_alexandrian ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_ALEXANDRIAN converts a JED to an Alexandrian YMDF date.
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
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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

  j_prime = j + 124

  y_prime =     ( 4 * j_prime + 3 ) / 1461
  t_prime = mod ( 4 * j_prime + 3, 1461 ) / 4
  m_prime = t_prime / 30
  d_prime = mod ( t_prime, 30 )
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime, 13 ) + 1
  y = y_prime - 4690 + ( 13 - m ) / 13

  return
end
subroutine jed_to_ymdf_armenian ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_ARMENIAN converts a JED to an Armenian YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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

  j_prime = j + 317

  y_prime =       j_prime / 365
  t_prime = mod ( j_prime, 365 )
  m_prime =       t_prime / 30
  d_prime = mod ( t_prime, 30 )
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime, 13 ) + 1
  y = y_prime - 5268 + ( 13 - m ) / 13

  return
end
subroutine jed_to_ymdf_bahai ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_BAHAI converts a JED to a Bahai YMDF date.
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
!    Algorithm F,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, pages 324-325.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  integer ( kind = 4 ) g
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

  g = ( 3 * ( ( 4 * j + 274273 ) / 146097 ) / 4 ) - 50
  j_prime = j + 1412 + g

  y_prime =     ( 4 * j_prime + 3 ) / 1461
  t_prime = mod ( 4 * j_prime + 3, 1461 ) / 4
  m_prime =       t_prime / 19
  d_prime = mod ( t_prime, 19 )
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime + 19, 20 ) + 1
  y = y_prime - 6560 + ( 39 - m ) / 20

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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
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
subroutine jed_to_ymdf_coptic ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_COPTIC converts a JED to a Coptic YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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

  j_prime = j + 124

  y_prime =     ( 4 * j_prime + 3 ) / 1461
  t_prime = mod ( 4 * j_prime + 3, 1461 ) / 4
  m_prime = t_prime / 30
  d_prime = mod ( t_prime, 30 )
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime, 13 ) + 1
  y = y_prime - 4996 + ( 13 - m ) / 13

  return
end
subroutine jed_to_ymdf_eg_civil ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_EG_CIVIL converts a JED to an Egyptian Civil YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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

  j_prime = j + 47 - 1

  y_prime =       j_prime / 365
  t_prime = mod ( j_prime, 365 )
  m_prime =       t_prime / 30
  d_prime = mod ( t_prime, 30 )
!
!  Convert the computational date to calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime, 13 ) + 1
  y = y_prime - 3968 + ( 13 - m ) / 13

  return
end
subroutine jed_to_ymdf_eg_lunar ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_EG_LUNAR converts a JED to an Egyptian Lunar YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_eg_lunar
  integer ( kind = 4 ) ncycle
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_eg_lunar

  call epoch_to_jed_eg_lunar ( jed_epoch )

  j = int ( jed - jed_epoch )
  f = ( jed - jed_epoch ) - real ( j, kind = 8 )

  d = 1 + j
  m = 1
  y = 1
!
!  Account for the number of 25 year cycles of 9125 days.
!
  ncycle = d / 9125
  y = y + 25 * ncycle
  d = d - ncycle * 9125

  do while ( year_length_eg_lunar ( y ) < d )
    d = d - year_length_eg_lunar ( y )
    y = y + 1
  end do

  do while ( month_length_eg_lunar ( y, m ) < d )
    d = d - month_length_eg_lunar ( y, m )
    m = m + 1
  end do

  return
end
subroutine jed_to_ymdf_english ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_ENGLISH converts a JED to an English YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_transition
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  call transition_to_jed_english ( jed_transition )

  if ( jed <= jed_transition ) then
    call jed_to_ymdf_julian ( jed, y, m, d, f )
  else
    call jed_to_ymdf_gregorian ( jed, y, m, d, f )
  end if

  return
end
subroutine jed_to_ymdf_ethiopian ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_ETHIOPIAN converts a JED to an Ethiopian YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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

  j_prime = j + 124

  y_prime =     ( 4 * j_prime + 3 ) / 1461
  t_prime = mod ( 4 * j_prime + 3, 1461 ) / 4
  m_prime =       t_prime / 30
  d_prime = mod ( t_prime, 30 )
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime, 13 ) + 1
  y = y_prime - 4720 + ( 13 - m ) / 13

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
!    18 January 2009
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
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

  g = ( 3 * ( ( 4 * j + 274277 ) / 146097 ) / 4 ) - 38
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
subroutine jed_to_ymdf_gregorian2 ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_GREGORIAN2 converts a JED to a Gregorian YMDF date.
!
!  Discussion:
!
!    The theory behind this routine is very clean.  The Gregorian
!    calendar has cycles of 1, 4, 100 and 400 years, and we can
!    analyze a date by determining where it lies within these cycles.
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
!    Edward Reingold, Nachum Dershowitz, Stewart Clamen,
!    Calendrical Calculations, II: Three Historical Calendars,
!    Software - Practice and Experience,
!    Volume 23, Number 4, pages 383-404, April 1993.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d0
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) d3
  integer ( kind = 4 ) d4
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ), parameter :: g1 = 365
  integer ( kind = 4 ), parameter :: g4 = 1461
  integer ( kind = 4 ), parameter :: g100 = 36524
  integer ( kind = 4 ), parameter :: g400 = 146097
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n100
  integer ( kind = 4 ) n400
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call epoch_to_jed_gregorian ( jed_epoch )

  j = int ( jed - jed_epoch )
  f1 = ( jed - jed_epoch ) - real ( j, kind = 8 )

  d0 = j
  n400 = 0

  do while ( d0 < 0 )
    d0 = d0 + g400
    n400 = n400 - 1
  end do

  n400 = n400 + d0 / g400
  d1 = i4_modp ( d0, g400 )

  n100 = d1 / g100
  d2 = i4_modp ( d1, g100 )

  n4 = d2 / g4
  d3 = i4_modp ( d2, g4 )

  n1 = d3 / g1
  d4 = i4_modp ( d3, g1 )

  if ( n100 == 4 .or. n1 == 4 ) then
    j1 = 366
    y1 = 400 * n400 + 100 * n100 + 4 * n4 + n1
  else
    j1 = d4 + 1
    y1 = 400 * n400 + 100 * n100 + 4 * n4 + n1 + 1
  end if
!
!  Any year before 1 AD must be moved one year further back, since
!  this calendar does not include a year 0.
!
  call y_astronomical_to_common ( y1, y1 )

  call yjf_to_ymdf_gregorian ( y1, j1, f1, y, m, d, f )

  return
end
subroutine jed_to_ymdf_hebrew ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_HEBREW converts a JED to a Hebrew YMDF date.
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
!    Algorithm I,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 334.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) j1
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) type
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call jed_to_year_hebrew ( jed, y1 )

  call new_year_to_jed_hebrew ( y1, jed2 )

  call year_to_type_hebrew ( y1, type )

  j1 = int ( jed - jed2 )
  f1 = ( jed - jed2 ) - real ( j1, kind = 8 )

  j1 = j1 + 1

  call yjf_to_ymdf_hebrew ( y1, j1, f1, y, m, d, f )

  return
end
subroutine jed_to_ymdf_hindu_solar ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_HINDU_SOLAR converts a JED to a Hindu solar YMDF date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Reingold, Nachum Dershowitz, Stewart Clamen,
!    Calendrical Calculations, II: Three Historical Calendars,
!    Software - Practice and Experience,
!    Volume 23, Number 4, pages 383-404, April 1993.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jf
  integer ( kind = 4 ) m
  real ( kind = 8 ) month_length_hindu_solar
  integer ( kind = 4 ) y
  real ( kind = 8 ) year_length_hindu_solar

!  Date = JF.
!
  jf = jed
!
!  Date = JED_EPOCH + JF
!
  call epoch_to_jed_hindu_solar ( jed_epoch )
  jf = jf - jed_epoch
!
!  Date = JED_EPOCH + Y years + JF
!
  y =  int ( jf / year_length_hindu_solar ( ) )
  jf = jf - y * year_length_hindu_solar ( )
!
!  Date = JED_EPOCH + Y years + ( M - 1 ) months + JF
!
  m = 1 + int ( jf / month_length_hindu_solar ( ) )
  jf = jf - ( real ( m - 1, kind = 8 ) * month_length_hindu_solar ( ) )
!
!  Date = JED_EPOCH + Y years + ( M - 1 ) months + ( D - 1 ) days + f
!
  d = int ( jf ) + 1
  f = jf - ( d - 1 )

  return
end
subroutine jed_to_ymdf_islamic_a ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_ISLAMIC_A converts a JED to an Islamic A YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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

  j_prime = j + 7665

  y_prime =     ( 30 * j_prime + 15 ) / 10631
  t_prime = mod ( 30 * j_prime + 15, 10631 ) / 30
  m_prime =     ( 100 * t_prime + 10 ) / 2951
  d_prime = mod ( 100 * t_prime + 10, 2951 ) / 100
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime, 12 ) + 1
  y = y_prime - 5519 + ( 12 - m ) / 12

  return
end
subroutine jed_to_ymdf_islamic_b ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_ISLAMIC_B converts a JED to an Islamic B YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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

  j_prime = j + 7664

  y_prime =     ( 30 * j_prime + 15 ) / 10631
  t_prime = mod ( 30 * j_prime + 15, 10631 ) / 30
  m_prime =     ( 100 * t_prime + 10 ) / 2951
  d_prime = mod ( 100 * t_prime + 10, 2951 ) / 100
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime, 12 ) + 1
  y = y_prime - 5519 + ( 12 - m ) / 12

  return
end
subroutine jed_to_ymdf_jelali ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_JELALI converts a JED to a Jelali YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) y

  call epoch_to_jed_jelali ( jed_epoch )

  j = int ( jed - jed_epoch )
  f = ( jed - jed_epoch ) - real ( j, kind = 8 )

  d = 1 + j
  m = 1
  y = 1
!
!  Account for the number of completed 4 year cycles of 1461 days.
!
  n = ( d - 1 ) / 1461
  y = y + 4 * n
  d = d - n * 1461
!
!  Account for the number of completed 365 day years.
!
  n = ( d - 1 ) / 365
  y = y + n
  d = d - n * 365
!
!  Account for the number of completed 30 day months.
!
  n = ( d - 1 ) / 30
  m = m + n
  d = d - n * 30

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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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
subroutine jed_to_ymdf_julian2 ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_JULIAN2 converts a JED to a Julian YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jd
  integer ( kind = 4 ) je
  real ( kind = 8 ) jed
  integer ( kind = 4 ) jg
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the input.
!
  call jed_check ( jed, ierror )

  if ( ierror /= 0 ) then
    y = -1
    m = -1
    d = -1
    f = -1.0D+00
    return
  end if

  j = int ( jed + 0.5D+00 )
  f = ( jed + 0.5D+00 ) - real ( j, kind = 8 )

  jd = int ( ( real ( j + 1524, kind = 8 ) - 122.1D+00 ) / 365.25D+00 )
  je = int ( 365.25D+00 * real ( jd, kind = 8 ) )
  jg = int ( real ( j + 1524 - je, kind = 8 ) / 30.6001D+00 )
!
!  Now compute D, M and Y.
!
  d = j + 1524 - je - int ( 30.6001D+00 * jg )

  if ( jg <= 13 ) then
    m = jg - 1
  else
    m = jg - 13
  end if

  if ( 2 < m ) then
    y = jd - 4716
  else
    y = jd - 4715
  end if
!
!  Any year before 1 AD must be moved one year further back, since
!  this calendar does not include a year 0.
!
  call y_astronomical_to_common ( y, y )

  return
end
subroutine jed_to_ymdf_julian3 ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_JULIAN3 converts a JED to a Julian YMDF date.
!
!  Discussion:
!
!    The theory behind this routine is very clean.  The Julian
!    calendar has cycles of 1 and 4 years, and we can analyze a date
!    by determining where it lies within these cycles.
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
!  Reference:
!
!    Edward Reingold, Nachum Dershowitz, Stewart Clamen,
!    Calendrical Calculations, II: Three Historical Calendars,
!    Software - Practice and Experience,
!    Volume 23, Number 4, pages 383-404, April 1993.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ), parameter :: cycle_length = 1461
  integer ( kind = 4 ) d
  integer ( kind = 4 ) d0
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ), parameter :: year_length = 365

  call epoch_to_jed_julian ( jed_epoch )

  j = int ( jed - jed_epoch )
  f1 = ( jed - jed_epoch ) - real ( j, kind = 8 )

  if ( f1 < 0.0D+00 ) then
    f1 = f1 + 1.0D+00
    j = j - 1
  end if

  d0 = j
  n4 = 0

  do while ( d0 <= 0 )
    d0 = d0 + cycle_length
    n4 = n4 - 1
  end do

  n4 = n4 + d0 / cycle_length
  d1 = i4_modp ( d0, cycle_length )
  n1 = d1 / year_length
  d2 = i4_modp ( d1, year_length )

  if ( n1 == 4 ) then
    j1 = 366
    y1 = 4 * n4 + n1
  else
    j1 = d2 + 1
    y1 = 4 * n4 + n1 + 1
  end if
!
!  Any year before 1 AD must be moved one year further back, since
!  this calendar does not include a year 0.
!
  call y_astronomical_to_common ( y1, y1 )

  call yjf_to_ymdf_julian ( y1, j1, f1, y, m, d, f )

  return
end
subroutine jed_to_ymdf_khwarizmian ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_KHWARIZMIAN converts a JED to a Khwarizmian YMDF date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2010
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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

  j_prime = j + 317

  y_prime =       j_prime / 365
  t_prime = mod ( j_prime, 365 )
  m_prime =       t_prime / 30
  d_prime = mod ( t_prime, 30 )
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime, 13 ) + 1
  y = y_prime - 5348 + ( 13 - m ) / 13

  return
end
subroutine jed_to_ymdf_macedonian ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_MACEDONIAN converts a JED to a Macedonian YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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
  m = mod ( m_prime + 6, 12 ) + 1
  y = y_prime - 4405 + ( 18 - m ) / 12

  return
end
subroutine jed_to_ymdf_persian ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_PERSIAN converts a JED to a Persian YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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

  j_prime = j + 77

  y_prime =       j_prime / 365
  t_prime = mod ( j_prime, 365 )
  m_prime =       t_prime / 30
  d_prime = mod ( t_prime, 30 )
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime + 9, 13 ) + 1
  y = y_prime - 5348 + ( 22 - m ) / 13

  return
end
subroutine jed_to_ymdf_republican ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_REPUBLICAN converts a JED to a Republican YMDF date.
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
!    Algorithm F,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, pages 324-325.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
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

  g = ( 3 * ( ( 4 * j + 578797 ) / 146097 ) / 4 ) - 51
  j_prime = j + 111 + g

  y_prime =     ( 4 * j_prime + 3 ) / 1461
  t_prime = mod ( 4 * j_prime + 3, 1461 ) / 4
  m_prime =     t_prime / 30
  d_prime = mod ( t_prime, 30 )
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime, 13 ) + 1
  y = y_prime - 6504 + ( 13 - m ) / 13

  return
end
subroutine jed_to_ymdf_roman ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_ROMAN converts a JED to a Roman YMDF date.
!
!  Discussion:
!
!    The Roman calendar used here is artificial.  It is assumed to begin
!    on the Julian calendar date 1 January 753 BC, and to be simply a
!    copy of the Julian calendar, shifted by 753 years.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
  integer ( kind = 4 ) yj

  call jed_to_ymdf_julian ( jed, yj, m, d, f )

  call y_julian_to_roman ( yj, y )

  return
end
subroutine jed_to_ymdf_saka ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_SAKA converts a JED to a Saka YMDF date.
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
!    Algorithm F,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, pages 324-325.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
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
  integer ( kind = 4 ) s
  integer ( kind = 4 ) t_prime
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
  integer ( kind = 4 ) z
!
!  Determine the computational date (Y'/M'/D').
!
  j = int ( jed + 0.5D+00 )
  f = ( jed + 0.5D+00 ) - real ( j, kind = 8 )

  g = ( 3 * ( ( 4 * j + 274073 ) / 146097 ) / 4 ) - 36

  j_prime = j + 1348 + g

  y_prime =     ( 4 * j_prime + 3 ) / 1461
  t_prime = mod ( 4 * j_prime + 3, 1461 ) / 4

  x = t_prime / 365
  z = t_prime / 185 - x
  s = 31 - z

  m_prime =             ( t_prime - 5 * z ) / s
  d_prime = 6 * x + mod ( t_prime - 5 * z, s )
!
!  Convert the computational date to a calendar date.
!
  d = d_prime + 1
  m = mod ( m_prime + 1, 12 ) + 1
  y = y_prime - 4794 + ( 13 - m ) / 12

  return
end
subroutine jed_to_ymdf_soor_san ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_SOOR_SAN converts a JED to a Soor San YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) y

  call epoch_to_jed_soor_san ( jed_epoch )

  j = int ( jed - jed_epoch )
  f = ( jed - jed_epoch ) - real ( j, kind = 8 )

  d = 1 + j
  m = 1
  y = 1
!
!  Account for the number of completed 4 year cycles of 1461 days.
!
  n = ( d - 1 ) / 1461
  y = y + 4 * n
  d = d - n * 1461
!
!  Account for the number of completed 365 day years.
!
  n = ( d - 1 ) / 365
  y = y + n
  d = d - n * 365
!
!  Account for the number of completed 30 day months.
!
  n = ( d - 1 ) / 30
  m = m + n
  d = d - n * 30

  return
end
subroutine jed_to_ymdf_syrian ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_SYRIAN converts a JED to a Syrian YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
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
  m = mod ( m_prime + 5, 12 ) + 1
  y = y_prime - 4405 + ( 17 - m ) / 12

  return
end
subroutine jed_to_ymdf_zoroastrian ( jed, y, m, d, f )

!*****************************************************************************80
!
!! JED_TO_YMDF_ZOROASTRIAN converts a JED to a Zoroastrian YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) years

  call epoch_to_jed_zoroastrian ( jed_epoch )

  j = int ( jed - jed_epoch )
  f = ( jed - jed_epoch ) - real ( j, kind = 8 )

  d = 1 + j
  m = 1
  y = 1

  years = ( d - 1 ) / 365
  y = y + years
  d = d - years * 365

  months = ( d - 1 ) / 30
  m = m + months
  d = d - months * 30

  return
end
subroutine jed_to_ymdhms_common ( jed, y, m, d, h, n, s )

!*****************************************************************************80
!
!! JED_TO_YMDHMS_COMMON converts a JED to a Common YMDHMS date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
!    Output, integer ( kind = 4 ) Y, M, D, H, N, S, the YMDHMS date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  integer ( kind = 4 ) h
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_transition
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call transition_to_jed_common ( jed_transition )

  if ( jed <= jed_transition ) then
    call jed_to_ymdf_julian ( jed, y1, m1, d1, f1 )
  else
    call jed_to_ymdf_gregorian ( jed, y1, m1, d1, f1 )
  end if

  call ymdf_to_ymdhms_common ( y1, m1, d1, f1, y, m, d, h, n, s )

  return
end
function lower ( s )

!*****************************************************************************80
!
!! LOWER returns a lowercase version of a string.
!
!  Discussion:
!
!    LOWER is a string function of undeclared length.  The length
!    of the argument returned is determined by the declaration of
!    LOWER in the calling routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string.
!
!    Output, character ( len = * ) LOWER, a lowercase copy of the string.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character ( len = * ) lower
  integer ( kind = 4 ) n
  character ( len = * ) s

  lower = s

  n = len_trim ( lower )

  do i = 1, n

    j = ichar ( lower(i:i) )

    if ( 65 <= j .and. j <= 90 ) then
      lower(i:i) = char ( j + 32 )
    end if

  end do

  return
end
subroutine mayan_long_to_jed ( pictun, baktun, katun, tun, uinal, kin, f, jed )

!*****************************************************************************80
!
!! MAYAN_LONG_TO_JED converts a Mayan long date to a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Algorithm L,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 341.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PICTUN, BAKTUN, KATUN, TUN, UINAL, KIN, values
!    defining the Mayan long date.
!
!    Input, real ( kind = 8 ) F, the fractional part of the date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) baktun
  integer ( kind = 4 ) days
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) katun
  integer ( kind = 4 ) kin
  integer ( kind = 4 ) pictun
  integer ( kind = 4 ) tun
  integer ( kind = 4 ) uinal

  days = (((((   pictun   * 20 + baktun ) * 20 + katun  ) * 20 &
    + tun    ) * 18 + uinal  ) * 20 + kin )

  call epoch_to_jed_mayan_long ( jed_epoch )

  jed = jed_epoch + real ( days, kind = 8 ) + f

  return
end
subroutine mayan_round_to_jed ( y, a, b, c, d, f, jed )

!*****************************************************************************80
!
!! MAYAN_ROUND_TO_JED converts a Mayan round date to a JED.
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
!    Algorithm L,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 341.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, A, B, C, D, values defining the Mayan
!    round date.
!
!    Input, real ( kind = 8 ) F, the fractional part of the date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) i4_modp
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) r
  integer ( kind = 4 ) y

  m = 13 * i4_modp ( 60 + 3 * ( a - b ), 20 ) + a - 1
  m = i4_modp ( m + 101, 260 )
  n = 20 * d + c
  n = i4_modp ( n + 17, 365 )
  r = 365 * i4_modp ( 364 + m - n, 52 ) + n

  call epoch_to_jed_mayan_long ( jed_epoch )
  jed = jed_epoch + real ( 18980 * y + r, kind = 8 ) + f

  return
end
subroutine minute_borrow_common ( y, m, d, h, n )

!*****************************************************************************80
!
!! MINUTE_BORROW_COMMON "borrows" an hour of minutes in a Common date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, H, N, the year,
!    month, day, hour and minute representing a date.  On input, N
!    might be negative.
!    On output, H should have decreased by one, and N gone up by 60.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) y

  do while ( n < 0 )

    n = n + 60
    h = h - 1

    call hour_borrow_common ( y, m, d, h )

  end do

  return
end
subroutine minute_carry_common ( y, m, d, h, n )

!*****************************************************************************80
!
!! MINUTE_CARRY_COMMON: given a Common YMDHMS date, carries minutes to hours.
!
!  Algorithm:
!
!    While 60 <= N:
!
!      decrease N by the number of minutes in an hour;
!      increase H by 1;
!      if necessary, adjust Y, M and D.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, H, N, the date.
!    On output, N is between 0 and 59.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) y

  do while ( 60 <= n )

    n = n - 60
    h = h + 1

    call hour_carry_common ( y, m, d, h )

  end do

  return
end
subroutine mjd_to_jed ( mjd, jed )

!*****************************************************************************80
!
!! MJD_TO_JED converts a modified JED to a JED.
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
!    Input, real ( kind = 8 ) MJD, the modified Julian Ephemeris Date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) mjd

  call epoch_to_jed_mjd ( jed_epoch )
  jed = mjd + jed_epoch

  return
end
subroutine month_borrow_alexandrian ( y, m )

!*****************************************************************************80
!
!! MONTH_BORROW_ALEXANDRIAN borrows a year of months on Alexandrian calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!    On input, M might be negative.  On output, Y should have decreased by
!    one, and M gone up by the number of months in the year that we
!    "cashed in".  The routine knows there was no year 0.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_alexandrian

  do while ( m <= 0 )

    months = year_length_months_alexandrian ( y )

    m = m + months
    y = y - 1

  end do

  return
end
subroutine month_borrow_bahai ( y, m )

!*****************************************************************************80
!
!! MONTH_BORROW_BAHAI borrows a year of months on the Bahai calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!    On input, M might be negative.  On output, Y should have decreased by
!    one, and M gone up by the number of months in the year that we
!    "cashed in".  The routine knows there was no year 0.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_bahai

  do while ( m <= 0 )

    months = year_length_months_bahai ( y )

    m = m + months
    y = y - 1

  end do

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
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
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
subroutine month_borrow_eg_civil ( y, m )

!*****************************************************************************80
!
!! MONTH_BORROW_EG_CIVIL borrows a year of months on Egyptian Civil calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_eg_civil

  do while ( m <= 0 )

    months = year_length_months_eg_civil ( y )

    m = m + months
    y = y - 1

  end do

  return
end
subroutine month_borrow_english ( y, m )

!*****************************************************************************80
!
!! MONTH_BORROW_ENGLISH borrows a year of months on the English calendar.
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
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_english

  do while ( m <= 0 )

    months = year_length_months_english ( y )

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
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
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
subroutine month_borrow_hebrew ( y, m )

!*****************************************************************************80
!
!! MONTH_BORROW_HEBREW borrows a year of months on the Hebrew calendar.
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
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_hebrew

  do while ( m <= 0 )

    months = year_length_months_hebrew ( y )

    m = m + months
    y = y - 1

  end do

  return
end
subroutine month_borrow_islamic ( y, m )

!*****************************************************************************80
!
!! MONTH_BORROW_ISLAMIC borrows a year of months on the Islamic calendar.
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
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_islamic

  do while ( m <= 0 )

    months = year_length_months_islamic ( y )

    m = m + months
    y = y - 1

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
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
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
subroutine month_borrow_republican ( y, m )

!*****************************************************************************80
!
!! MONTH_BORROW_REPUBLICAN borrows a year of months on the Republican calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_republican

  do while ( m <= 0 )

    months = year_length_months_republican ( y )

    m = m + months
    y = y - 1

  end do

  return
end
subroutine month_borrow_roman ( y, m )

!*****************************************************************************80
!
!! MONTH_BORROW_ROMAN borrows a year of months on the Roman calendar.
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
!    27 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_roman

  do while ( m <= 0 )

    months = year_length_months_roman ( y )

    m = m + months
    y = y - 1

    if ( y == 0 ) then
      y = - 1
    end if

  end do

  return
end
subroutine month_cal_common ( y, m )

!*****************************************************************************80
!
!! MONTH_CAL_COMMON prints a Common month calendar.
!
!  Discussion:
!
!    The "common" calendar is meant to be the calendar which is Julian up to
!    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
!
!  Format:
!
!    COMMON CALENDAR
!    APRIL 1997 CE
!
!    Su  M Tu  W Th  F Sa
!           1  2  3  4  5
!     6  7  8  9 10 11 12
!    13 14 15 16 17 18 19
!    20 21 22 23 24 25 26
!    27 28 29 30
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_max
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  integer ( kind = 4 ) iday
  integer ( kind = 4 ) ierror
  character ( len = 2 ) label(7)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) month_length_common
  character ( len = 9 ) s1
  character ( len = 10 ) s2
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Make local copies of the input.
!
  m2 = m
  y2 = y
!
!  Check the month and year.  After this call, month is
!  guaranteed to be between 1 and 12.
!
  call ym_check_common ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Find the day of the week for Y M 1.
!
  d = 1
  f = 0.0D+00

  call ymdf_to_weekday_common ( y2, m2, d, f, w )
!
!  Find the appropriate label for the first box in the calendar.
!
  iday = 2 - w
!
!  Print out a heading.
!
  call month_to_month_name_common ( m2, s1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMMON CALENDAR'
  call y_to_s_common ( y2, s2 )
  write ( *, '(a,1x,a)' ) trim ( s1 ), trim ( s2 )
  write ( *, '(a)' ) ' '
!
!  Get the days of the week.
!
  do w = 1, 7
    call weekday_to_name_common2 ( w, label(w) )
  end do

  write ( *, '(1x,7a3)' ) label(1:7)
!
!  Print out a line of day numbers.
!  IDAY keeps track of the numerical day,
!  D2 keeps track of the label for the day, which  only differed
!  from IDAY in October 1582.
!
  d2 = iday
  f2 = 0.0D+00
  d_max = month_length_common ( y2, m2 )

  do while ( iday <= d_max )

    do w = 1, 7

      if ( iday < 1 ) then
        label(w) = ' '
      else if ( d_max < iday ) then
        label(w) = ' '
      else
        write ( label(w), '(i2)' ) d2
      end if

      iday = iday + 1

      call ymdf_next_common ( y2, m2, d2, f2, y2, m2, d2, f2 )

    end do

    write ( *, '(1x,7a3)' ) label(1:7)

  end do

  return
end
subroutine month_cal_english ( y, m )

!*****************************************************************************80
!
!! MONTH_CAL_ENGLISH prints an English month calendar.
!
!  Format:
!
!    ENGLISH CALENDAR
!        APRIL 1997 AD NS
!
!    Su  M Tu  W Th  F Sa
!           1  2  3  4  5
!     6  7  8  9 10 11 12
!    13 14 15 16 17 18 19
!    20 21 22 23 24 25 26
!    27 28 29 30
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
!    Input, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_max
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  integer ( kind = 4 ) iday
  integer ( kind = 4 ) ierror
  character ( len = 2 ) label(7)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) month_length_english
  character ( len = 9 ) s1
  character ( len = 10 ) s2
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Make local copies of the input.
!
  m2 = m
  y2 = y
!
!  Check the month and year.  After this call, month is
!  guaranteed to be between 1 and 12.
!
  call ym_check_english ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Find the day of the week for Y M 1.
!
  d = 1
  f = 0.0D+00

  call ymdf_to_weekday_english ( y2, m2, d, f, w )
!
!  Find the appropriate label for the first box in the calendar.
!
  iday = 2 - w
!
!  Print out a heading.
!
  call month_to_month_name_common ( m2, s1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ENGLISH CALENDAR'
  call y_to_s_english ( y2, s2 )
  write ( *, '(a,1x,a)' ) trim ( s1 ), trim ( s2 )
  write ( *, '(a)' ) ' '
!
!  Get the days of the week.
!
  do w = 1, 7
    call weekday_to_name_common2 ( w, label(w) )
  end do

  write ( *, '(1x,7a3)' ) label(1:7)
!
!  Print out a line of day numbers.
!  IDAY keeps track of the numerical day,
!  D2 keeps track of the label for the day, which
!  only differed from IDAY in September 1752.
!
  d2 = iday
  f2 = 0.0D+00
  d_max = month_length_english ( y2, m2 )

  do while ( iday <= d_max )

    do w = 1, 7

      if ( iday < 1 ) then
        label(w) = ' '
      else if ( d_max < iday ) then
        label(w) = ' '
      else
        write ( label(w), '(i2)' ) d2
      end if

      iday = iday + 1

      call ymdf_next_english ( y2, m2, d2, f2, y2, m2, d2, f2 )

    end do

    write ( *, '(1x,7a3)' ) label(1:7)

  end do

  return
end
subroutine month_cal_gregorian ( y, m )

!*****************************************************************************80
!
!! MONTH_CAL_GREGORIAN prints a Gregorian month calendar.
!
!  Format:
!
!    GREGORIAN CALENDAR
!        APRIL 1997 AD
!
!    Su  M Tu  W Th  F Sa
!           1  2  3  4  5
!     6  7  8  9 10 11 12
!    13 14 15 16 17 18 19
!    20 21 22 23 24 25 26
!    27 28 29 30
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) iday
  integer ( kind = 4 ) ierror
  character ( len = 2 ) label(7)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) month_length_gregorian
  character ( len = 9 ) s1
  character ( len = 10 ) s2
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Make local copies of the input.
!
  m2 = m
  y2 = y
!
!  Check the month and year.  After this call, month is
!  guaranteed to be between 1 and 12.
!
  call ym_check_gregorian ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Find the day of the week for Y M 1.
!
  d = 1
  f = 0.0D+00

  call ymdf_to_weekday_gregorian ( y2, m2, d, f, w )
!
!  Find the appropriate label for the first box in the calendar.
!
  iday = 2 - w
!
!  Print out a heading.
!
  call month_to_month_name_common ( m2, s1 )
  call y_to_s_gregorian ( y, s2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GREGORIAN CALENDAR'
  write ( *, '(a,1x,a)' ) trim ( s1 ), trim ( s2 )
  write ( *, '(a)' ) ' '
!
!  Get the days of the week.
!
  do w = 1, 7
    call weekday_to_name_common2 ( w, label(w) )
  end do

  write ( *, '(1x,7a3)' ) label(1:7)
!
!  Print out a line of day numbers.
!
  do while ( iday <= month_length_gregorian ( y2, m2 ) )

    do w = 1, 7

      if ( iday < 1 ) then
        label(w) = ' '
      else if ( month_length_gregorian ( y2, m2 ) < iday ) then
        label(w) = ' '
      else
        write ( label(w), '(i2)' ) iday
      end if

      iday = iday + 1

    end do

    write ( *, '(1x,7a3)' ) label(1:7)

  end do

  return
end
subroutine month_cal_hebrew ( y, m )

!*****************************************************************************80
!
!! MONTH_CAL_HEBREW prints a Hebrew month calendar.
!
!  Format:
!
!    HEBREW CALENDAR
!    month year AM
!
!    Su  M Tu  W Th  F Sa
!           1  2  3  4  5
!     6  7  8  9 10 11 12
!    13 14 15 16 17 18 19
!    20 21 22 23 24 25 26
!    27 28 29 30
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
!    Input, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) iday
  integer ( kind = 4 ) ierror
  character ( len = 2 ) label(7)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) month_length_hebrew
  character ( len = 9 ) s1
  character ( len = 10 ) s2
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Make local copies of the input.
!
  m2 = m
  y2 = y
!
!  Check the month and year.  After this call, month is
!  guaranteed to be between 1 and 12.
!
  call ym_check_hebrew ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Find the day of the week for Y M 1.
!
  d = 1
  f = 0.0D+00

  call ymdf_to_weekday_hebrew ( y2, m2, d, f, w )
!
!  Find the appropriate label for the first box in the calendar.
!
  iday = 2 - w
!
!  Print out a heading.
!
  call month_to_month_name_hebrew ( y2, m2, s1 )
  call y_to_s_hebrew ( y, s2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEBREW CALENDAR'
  write ( *, '(a,1x,a)' ) trim ( s1 ), trim ( s2 )
  write ( *, '(a)' ) ' '
!
!  Get the days of the week.
!
  do w = 1, 7
    call weekday_to_name_common2 ( w, label(w) )
  end do

  write ( *, '(1x,7a3)' ) label(1:7)
!
!  Print out a line of day numbers.
!
  do while ( iday <= month_length_hebrew ( y2, m2 ) )

    do w = 1, 7

      if ( iday < 1 ) then
        label(w) = ' '
      else if ( month_length_hebrew ( y2, m2 ) < iday ) then
        label(w) = ' '
      else
        write ( label(w), '(i2)' ) iday
      end if

      iday = iday + 1

    end do

    write ( *, '(1x,7a3)' ) label(1:7)

  end do

  return
end
subroutine month_cal_islamic_a ( y, m )

!*****************************************************************************80
!
!! MONTH_CAL_ISLAMIC_A prints an Islamic-A month calendar.
!
!  Format:
!
!    ISLAMIC-A CALENDAR
!    month year AH
!
!    Su  M Tu  W Th  F Sa
!           1  2  3  4  5
!     6  7  8  9 10 11 12
!    13 14 15 16 17 18 19
!    20 21 22 23 24 25 26
!    27 28 29 30
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
!    Input, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) iday
  integer ( kind = 4 ) ierror
  character ( len = 2 ) label(7)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) month_length_islamic
  character ( len = 9 ) s1
  character ( len = 10 ) s2
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Make local copies of the input.
!
  m2 = m
  y2 = y
!
!  Check the month and year.  After this call, month is
!  guaranteed to be between 1 and 12.
!
  call ym_check_islamic ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Find the day of the week for Y M 1.
!
  d = 1
  f = 0.0D+00

  call ymdf_to_weekday_islamic_a ( y2, m2, d, f, w )
!
!  Find the appropriate label for the first box in the calendar.
!
  iday = 2 - w
!
!  Print out a heading.
!
  call month_to_month_name_islamic ( m2, s1 )
  call y_to_s_islamic ( y, s2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ISLAMIC-A CALENDAR'
  write ( *, '(a,1x,a)' ) trim ( s1 ), trim ( s2 )
  write ( *, '(a)' ) ' '
!
!  Get the days of the week.
!
  do w = 1, 7
    call weekday_to_name_common2 ( w, label(w) )
  end do

  write ( *, '(1x,7a3)' ) label(1:7)
!
!  Print out a line of day numbers.
!
  do while ( iday <= month_length_islamic ( y2, m2 ) )

    do w = 1, 7

      if ( iday < 1 ) then
        label(w) = ' '
      else if ( month_length_islamic ( y2, m2 ) < iday ) then
        label(w) = ' '
      else
        write ( label(w), '(i2)' ) iday
      end if

      iday = iday + 1

    end do

    write ( *, '(1x,7a3)' ) label(1:7)

  end do

  return
end
subroutine month_cal_julian ( y, m )

!*****************************************************************************80
!
!! MONTH_CAL_JULIAN prints a Julian month calendar.
!
!  Format:
!
!    JULIAN CALENDAR
!    APRIL 1997 AD
!
!    Su  M Tu  W Th  F Sa
!           1  2  3  4  5
!     6  7  8  9 10 11 12
!    13 14 15 16 17 18 19
!    20 21 22 23 24 25 26
!    27 28 29 30
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
!    Input, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) iday
  integer ( kind = 4 ) ierror
  character ( len = 2 ) label(7)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) month_length_julian
  character ( len = 9 ) s1
  character ( len = 10 ) s2
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Make local copies of the input.
!
  m2 = m
  y2 = y
!
!  Check the month and year.  After this call, month is
!  guaranteed to be between 1 and 12.
!
  call ym_check_julian ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Find the day of the week for Y/M/1.
!
  d = 1
  f = 0.0D+00

  call ymdf_to_weekday_julian ( y2, m2, d, f, w )
!
!  Find the appropriate label for the first box in the calendar.
!
  iday = 2 - w
!
!  Print out a heading.
!
  call month_to_month_name_common ( m2, s1 )
  call y_to_s_julian ( y, s2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'JULIAN CALENDAR'
  write ( *, '(a,1x,a)' ) trim ( s1 ), trim ( s2 )
  write ( *, '(a)' ) ' '
!
!  Get the days of the week.
!
  do w = 1, 7
    call weekday_to_name_common2 ( w, label(w) )
  end do

  write ( *, '(1x,7a3)' ) label(1:7)
!
!  Print out a line of day numbers.
!
  do while ( iday <= month_length_julian ( y2, m2 ) )

    do w = 1, 7

      if ( iday < 1 ) then
        label(w) = ' '
      else if ( month_length_julian ( y2, m2 ) < iday ) then
        label(w) = ' '
      else
        write ( label(w), '(i2)' ) iday
      end if

      iday = iday + 1

    end do

    write ( *, '(1x,7a3)' ) label(1:7)

  end do

  return
end
subroutine month_cal_republican ( y, m )

!*****************************************************************************80
!
!! MONTH_CAL_REPUBLICAN prints a Republican month calendar.
!
!  Format:
!
!    REPUBLICAN CALENDAR
!    Brumaire 3 ER
!
!     1  2  3  4  5  6  7  8  9 10
!    11 12 13 14 15 16 17 18 19 20
!    21 22 23 24 25 26 27 28 29 30
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) iday
  integer ( kind = 4 ) ierror
  character ( len = 2 ) label(10)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) month_length_republican
  character ( len = 15 ) s1
  character ( len = 10 ) s2
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Make local copies of the input.
!
  m2 = m
  y2 = y
!
!  Check the month and year.
!
  call ym_check_republican ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Print out a heading.
!
  call month_to_month_name_republican ( m2, s1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REPUBLICAN CALENDAR'
  call y_to_s_republican ( y2, s2 )
  write ( *, '(a,1x,a)' ) trim ( s1 ), trim ( s2 )
  write ( *, '(a)' ) ' '
!
!  Print out a line of day numbers.
!
  iday = 1

  do while ( iday <= month_length_republican ( y2, m2 ) )

    do w = 1, 10

      if ( month_length_republican ( y2, m2 ) < iday ) then
        label(w) = ' '
      else
        write ( label(w), '(i2)' ) iday
      end if

      iday = iday + 1

    end do

    write ( *, '(1x,10a3)' ) label(1:10)

  end do

  return
end
subroutine month_cal_roman ( y, m )

!*****************************************************************************80
!
!! MONTH_CAL_ROMAN prints a Roman month calendar.
!
!  Format:
!
!    ROMAN CALENDAR
!    Aprilis DVI AUC
!
!    Kalends Aprilis
!    Ante diem iv Nones Aprilis
!    Ante diem iii Nones Aprilis
!    Pridie Nones Aprilis
!    Nones Aprilis
!    Ante diem viii Ides Aprilis
!    Ante diem vii Ides Aprilis
!    Ante diem vi Ides Aprilis
!    Ante diem v Ides Aprilis
!    Ante diem iv Ides Aprilis
!    Ante diem iii Ides Aprilis
!    Pridie Ides Aprilis
!    Ides Aprilis
!    Ante diem xvii Kalends Maius
!    Ante diem xvi Kalends Maius
!    ...
!    Ante diem iv Kalends Maius
!    Ante diem iii Kalends Maius
!    Pridie Kalends Maius
!
!  Discussion:
!
!    "AUC" means "ab urbe condita", or "from the founding of the city".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) iday
  integer ( kind = 4 ) ides
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) jday
  integer ( kind = 4 ) last
  character ( len = 10 ) lower
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ) month_length_roman
  integer ( kind = 4 ) nones
  character ( len = 40 ) s_date
  character ( len = 10 ) s_day
  character ( len = 15 ) s_month
  character ( len = 15 ) s_month_next
  character ( len = 10 ) s_year
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical   year_is_leap_roman
!
!  Make local copies of the input.
!
  m2 = m
  y2 = y
!
!  Check the month and year.
!
  call ym_check_roman ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Print out a heading.
!
  call month_to_month_name_roman ( m2, s_month )

  m3 = i4_wrap ( m2+1, 1, 12 )

  call month_to_month_name_roman ( m3, s_month_next )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ROMAN CALENDAR'
  call i4_to_roman ( y2, s_year )
  call s_cat ( s_year, ' AUC', s_year )
  write ( *, '(a,1x,a)' ) trim ( s_month ), trim ( s_year )
  write ( *, '(a)' ) ' '

  call month_to_nones_roman ( m2, nones )
  call month_to_ides_roman ( m2, ides )
  last = month_length_roman ( y2, m2 )

  do iday = 1, last

    if ( iday == 1 ) then
      s_date = 'Kalends ' // s_month
    else if ( iday < nones - 1 ) then
      jday = nones + 1 - iday
      call i4_to_roman ( jday, s_day )
      s_day = lower ( s_day )
      s_date = 'Ante diem ' // trim ( s_day ) // ' Nones ' // s_month
    else if ( iday == nones - 1 ) then
      s_date = 'Pridie Nones ' // s_month
    else if ( iday == nones ) then
      s_date = 'Nones ' // s_month
    else if ( iday < ides - 1 ) then
      jday = ides + 1 - iday
      call i4_to_roman ( jday, s_day )
      s_day = lower ( s_day )
      s_date = 'Ante diem ' // trim ( s_day ) // ' Ides ' // s_month
    else if ( iday == ides - 1 ) then
      s_date = 'Pridie Ides ' // s_month
    else if ( iday == ides ) then
      s_date = 'Ides ' // s_month

    else if ( m2 == 2 .and. year_is_leap_roman ( y2 ) ) then

      if ( iday < 25 ) then
        jday = last + 1 - iday
        call i4_to_roman ( jday, s_day )
        s_day = lower ( s_day )
        s_date = 'Ante diem ' // trim ( s_day ) // ' Kalends ' // s_month_next
      else if ( iday == 25 ) then
        jday = last + 2 - iday
        call i4_to_roman ( jday, s_day )
        s_day = lower ( s_day )
        s_date = 'Ante diem Bis ' // trim ( s_day ) // ' Kalends ' &
          // s_month_next
      else if ( iday < last ) then
        jday = last + 2 - iday
        call i4_to_roman ( jday, s_day )
        s_day = lower ( s_day )
        s_date = 'Ante diem ' // trim ( s_day ) // ' Kalends ' // s_month_next
      else
        s_date = 'Pridie Kalends ' // s_month_next
      end if

    else if ( iday < last ) then
      jday = last + 2 - iday
      call i4_to_roman ( jday, s_day )
      s_day = lower ( s_day )
      s_date = 'Ante diem ' // trim ( s_day ) // ' Kalends ' // s_month_next
    else
      s_date = 'Pridie Kalends ' // s_month_next
    end if

    write ( *, '(a)' ) trim ( s_date )

  end do

  return
end
subroutine month_cal_store_common ( y, m, lines )

!*****************************************************************************80
!
!! MONTH_CAL_STORE_COMMON stores a Common month calendar.
!
!  Discussion:
!
!    The "common" calendar is meant to be the calendar which is Julian before
!    the transition date, and Gregorian afterwards, with the transition date
!    best specified as as JED = 2299160.
!
!  Format:
!
!           1  2  3  4  5
!     6  7  8  9 10 11 12
!    13 14 15 16 17 18 19
!    20 21 22 23 24 25 26
!    27 28 29 30
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
!    Input, integer ( kind = 4 ) Y,  M, the YM date.
!
!    Output, character ( len = 20 ) LINES(6), the six lines of the calendar.
!
  implicit none

  integer ( kind = 4 ) days
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) iday
  integer ( kind = 4 ) ierror
  character ( len = 20 ) lines(6)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) n_line
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2

  lines(1:6) = ' '
!
!  Make local copies of the input.
!
  m2 = m
  y2 = y
!
!  Check the month and year.  After this call, month is
!  guaranteed to be between 1 and 12.
!
  call ym_check_common ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Find the day of the week for Y M 1.
!
  d2 = 1
  f2 = 0.0D+00

  call ymdf_to_weekday_common ( y2, m2, d2, f2, w )

  days = month_length_common ( y2, m2 )
!
!  Find the appropriate label for the first box in the calendar.
!
  iday = 2 - w
!
!  Print out a line of day numbers.
!  IDAY keeps track of the numerical day,
!  JDAY keeps track of the label for the day, which differed in October 1582.
!
  d2 = iday
  f2 = 0.0D+00

  n_line = 0

  do while ( iday <= days )

    n_line = n_line + 1

    do w = 1, 7

      i1 = 3 * ( w - 1 ) + 1
      i2 = i1 + 1

      if ( 1 <= iday .and. iday <= days ) then
        write ( lines(n_line)(i1:i2), '(i2)' ) d2
      end if

      iday = iday + 1

      call ymdf_next_common ( y2, m2, d2, f2, y2, m2, d2, f2 )

    end do

  end do

  return
end
subroutine month_carry_alexandrian ( y, m )

!*****************************************************************************80
!
!! MONTH_CARRY_ALEXANDRIAN carries a year of months on the Alexandrian calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_alexandrian

  months = year_length_months_alexandrian ( y )

  do

    if ( m <= months ) then
      exit
    end if

    m = m - months
    y = y + 1

  end do

  return
end
subroutine month_carry_bahai ( y, m )

!*****************************************************************************80
!
!! MONTH_CARRY_BAHAI carries a year of months on the Bahai calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_bahai

  months = year_length_months_bahai ( y )

  do

    if ( m <= months ) then
      exit
    end if

    m = m - months
    y = y + 1

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
subroutine month_carry_eg_civil ( y, m )

!*****************************************************************************80
!
!! MONTH_CARRY_EG_CIVIL carries a year of months on the Egyptian Civil calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the year and month.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_eg_civil

  months = year_length_months_eg_civil ( y )

  do

    if ( m <= months ) then
      exit
    end if

    m = m - months
    y = y + 1

  end do

  return
end
subroutine month_carry_english ( y, m )

!*****************************************************************************80
!
!! MONTH_CARRY_ENGLISH carries a year of months on the English calendar.
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
!    Input/output, integer ( kind = 4 ) Y, M, the year and month.
!    On output, M is no greater than 12.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_english

  do

    months = year_length_months_english ( y )

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
subroutine month_carry_hebrew ( y, m )

!*****************************************************************************80
!
!! MONTH_CARRY_HEBREW carries a year of months on the Hebrew calendar.
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
!    Input/output, integer ( kind = 4 ) Y, M, the year and month.
!    On output, M is no greater than 12.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_hebrew

  do

    months = year_length_months_hebrew ( y )

    if ( m <= months ) then
      exit
    end if

    m = m - months
    y = y + 1

  end do

  return
end
subroutine month_carry_islamic ( y, m )

!*****************************************************************************80
!
!! MONTH_CARRY_ISLAMIC carries a year of months on the Islamic calendar.
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
!    Input/output, integer ( kind = 4 ) Y, M, the year and month.
!    On output, M is no greater than 12.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_islamic

  do

    months = year_length_months_islamic ( y )

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
subroutine month_carry_republican ( y, m )

!*****************************************************************************80
!
!! MONTH_CARRY_REPUBLICAN carries a year of months on the Republican calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
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
  integer ( kind = 4 ) year_length_months_republican

  do

    months = year_length_months_republican ( y )

    if ( m <= months ) then
      exit
    end if

    m = m - months
    y = y + 1

  end do

  return
end
subroutine month_carry_roman ( y, m )

!*****************************************************************************80
!
!! MONTH_CARRY_ROMAN carries a year of months on the Roman calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2000
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
  integer ( kind = 4 ) year_length_months_roman

  do

    months = year_length_months_roman ( y )

    if ( m <= months ) then
      exit
    end if

    m = m - months
    y = y + 1

  end do

  return
end
function month_length_alexandrian ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_ALEXANDRIAN returns the number of days in an Alexandrian month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_ALEXANDRIAN, the number of
!    days in the month.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(13) :: mdays = (/ &
    30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 5 /)
  integer ( kind = 4 ) month_length_alexandrian
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_alexandrian
!
!  Copy the input.
!
  m2 = m
  y2 = y

  if ( m2 < 1 .or. 13 < m2 ) then
    month_length_alexandrian = 0
  else
    month_length_alexandrian = mdays(m2)
  end if

  if ( m2 == 13 .and. year_is_leap_alexandrian ( y2 ) ) then
    month_length_alexandrian = month_length_alexandrian + 1
  end if

  return
end
function month_length_bahai ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_BAHAI returns the number of days in a Bahai month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2000
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_BAHAI, the number of
!    days in the month.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) month_length_bahai
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_bahai
!
!  Copy the input.
!
  m2 = m
  y2 = y
!
!  Check the input.
!
  call ym_check_bahai ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    month_length_bahai = 0
    return
  end if

  if ( m2 <= 18 .or. m2 == 20 ) then
    month_length_bahai = 19
  else if ( year_is_leap_bahai ( y2 ) ) then
    month_length_bahai = 5
  else
    month_length_bahai = 4
  end if

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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_COMMON, the number of days
!    in the month.
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
    month_length_common = 0
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
  month_length_common = mdays(m2)
!
!  If necessary, add 1 day for February 29.
!
  if ( m2 == 2 .and. year_is_leap_common ( y2 ) ) then
    month_length_common = month_length_common + 1
  end if

  return
end
function month_length_coptic ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_COPTIC returns the number of days in a Coptic month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_COPTIC, the number of days
!    in the month.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(13) :: mdays = (/ &
    30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 5 /)
  integer ( kind = 4 ) month_length_coptic
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_coptic
!
!  Copy the input.
!
  m2 = m
  y2 = y

  if ( m2 < 1 .or. 13 < m2 ) then
    month_length_coptic = 0
  else
    month_length_coptic = mdays(m2)
  end if

  if ( m2 == 13 .and. year_is_leap_coptic ( y2 ) ) then
    month_length_coptic = month_length_coptic + 1
  end if

  return
end
function month_length_eg_civil ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_EG_CIVIL returns the number of days in an Egyptian Civil month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2000
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_EG_CIVIL, the number of days
!    in the month.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_eg_civil
  integer ( kind = 4 ) y

  if ( m < 1 ) then
    month_length_eg_civil = 0
  else if ( m <= 12 ) then
    month_length_eg_civil = 30
  else if ( m == 13 ) then
    month_length_eg_civil = 5
  else
    month_length_eg_civil = 0
  end if

  return
end
function month_length_eg_lunar ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_EG_LUNAR returns the number of days in an Egyptian Lunar month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2000
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_EG_LUNAR, the number of days
!    in the month.
!
  implicit none

  integer ( kind = 4 ) last
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(13) :: mdays = (/ &
    29, 30, 29, 30, 29, 30, 29, 30, 29, 30, 29, 30, 30 /)
  integer ( kind = 4 ) month_length_eg_lunar
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_eg_lunar
  integer ( kind = 4 ) year_length_months_eg_lunar
!
!  Copy the input.
!
  m2 = m
  y2 = y

  last = year_length_months_eg_lunar ( y2 )

  if ( m2 < 1 ) then
    month_length_eg_lunar = 0
  else if ( m2 <= last ) then
    month_length_eg_lunar = mdays(m2)
  else
    month_length_eg_lunar = 0
  end if

  if ( m2 == last ) then
    if ( year_is_leap_eg_lunar ( y ) ) then
      month_length_eg_lunar = month_length_eg_lunar + 1
    end if
  end if

  return
end
function month_length_english ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_ENGLISH returns the number of days in an English month.
!
!  Discussion:
!
!    In the English calendar, September 1752 had only 19 days.
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_ENGLISH, the number of days
!    in the month.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(12) :: mdays = (/ &
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  integer ( kind = 4 ) month_length_english
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_english
!
!  Copy the input.
!
  m2 = m
  y2 = y
!
!  Check the input.
!
  call ym_check_english ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    month_length_english = 0
    return
  end if
!
!  Take care of special cases:
!
  if ( y2 == 1752 ) then
    if ( m2 == 9 ) then
      month_length_english = 19
      return
    end if
  end if
!
!  Get the number of days in the month.
!
  month_length_english = mdays(m2)
!
!  If necessary, add 1 day for February 29.
!
  if ( m2 == 2 .and. year_is_leap_english ( y2 ) ) then
    month_length_english = month_length_english + 1
  end if

  return
end
function month_length_ethiopian ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_ETHIOPIAN returns the number of days in an Ethiopian month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_ETHIOPIAN, the number of days
!    in the month.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(13) :: mdays = (/ &
    30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 5 /)
  integer ( kind = 4 ) month_length_ethiopian
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_ethiopian
!
!  Copy the input.
!
  m2 = m
  y2 = y

  if ( m2 < 1 .or. 13 < m2 ) then
    month_length_ethiopian = 0
  else
    month_length_ethiopian = mdays(m2)
  end if

  if ( m2 == 13 .and. year_is_leap_ethiopian ( y2 ) ) then
    month_length_ethiopian = month_length_ethiopian + 1
  end if

  return
end
function month_length_greek ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_GREEK returns the number of days in a Greek month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2000
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_GREEK, the number of days
!    in the month.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(13) :: mdays = (/ &
    30, 29, 30, 29, 30, 29, 29, 30, 29, 30, 29, 30, 29 /)
  integer ( kind = 4 ) month_length_greek
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_greek
  logical year_is_embolismic_greek
!
!  Copy the input.
!
  m2 = m
  y2 = y

  if ( m2 < 1 ) then
    month_length_greek = 0
    return
  end if
!
!  A 13-month year.
!
  if ( year_is_embolismic_greek ( y2 ) ) then

    if ( 13 < m2 ) then
      month_length_greek = 0
      return
    end if

    month_length_greek = mdays(m2)

    if ( m2 == 7 .and. year_is_leap_greek ( y2 ) ) then
      month_length_greek = month_length_greek + 1
    end if
!
!  A 12 month year.
!
  else

    if ( m2 <= 6 ) then
      month_length_greek = mdays(m2)
    else if ( m2 <= 12 ) then
      month_length_greek = mdays(m2+1)
    else
      month_length_greek = 0
    end if

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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_GREGORIAN, the number of days
!    in the month.
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
    month_length_gregorian = 0
    return
  end if
!
!  Get the number of days in the month.
!
  month_length_gregorian = mdays(m2)
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
function month_length_hebrew ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_HEBREW returns the number of days in a Hebrew month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 333.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, the year and month number.  Note that
!    some years only had 12 months in them, while others have 13.  If
!    Y only has 12 months, then the length of the 13th month is
!    returned as 0 days.
!
!    Output, integer ( kind = 4 ) MONTH_LENGTH_HEBREW, the number of days
!    in the month.
!
  implicit none

  integer ( kind = 4 ), dimension (6,13) :: a = reshape ( source = (/ &
    30,  30,  30,  30,  30,  30, &
    29,  29,  30,  29,  29,  30, &
    29,  30,  30,  29,  30,  30, &
    29,  29,  29,  29,  29,  29, &
    30,  30,  30,  30,  30,  30, &
    29,  29,  29,  30,  30,  30, &
    30,  30,  30,  29,  29,  30, &
    29,  29,  29,  30,  30,  29, &
    30,  30,  30,  29,  29,  29, &
    29,  29,  29,  30,  30,  30, &
    30,  30,  30,  29,  29,  29, &
    29,  29,  29,  30,  30,  30, &
     0,   0,   0,  29,  29,  29 /), shape = (/ 6, 13 /) )

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) month_length_hebrew
  integer ( kind = 4 ) type
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Copy the input
!
  y2 = y
  m2 = m
!
!  Check the input.
!
  call ym_check_hebrew ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    month_length_hebrew = 0
    return
  end if

  call year_to_type_hebrew ( y2, type )

  if ( type < 1 .or. 6 < type ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MONTH_LENGTH_HEBREW - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal year TYPE = ', type
    write ( *, '(a,i6)' ) '  Y = ', y2
    stop
  end if

  if ( m2 < 1 .or. 13 < m2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MONTH_LENGTH_HEBREW - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal MONTH = ', m2
    stop
  end if

  month_length_hebrew = a(type,m2)

  return
end
function month_length_hindu_solar ( )

!*****************************************************************************80
!
!! MONTH_LENGTH_HINDU_SOLAR returns the number of days in a Hindu solar month.
!
!  Discussion:
!
!    Warning: this is a DOUBLE PRECISION quantity, with a fractional part!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MONTH_LENGTH_HINDU_SOLAR, the number of
!    days in the month.
!
  implicit none

  real ( kind = 8 ) month_length_hindu_solar
  real ( kind = 8 ) year_length_hindu_solar

  month_length_hindu_solar = year_length_hindu_solar ( ) / 12.0D+00

  return
end
function month_length_iranian ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_IRANIAN returns the number of days in an Iranian month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 November 2000
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_IRANIAN, the number of
!    days in the month.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(12) :: mdays = (/ &
    31, 31, 31, 31, 31, 31, 30, 30, 30, 30, 30, 29 /)
  integer ( kind = 4 ) month_length_iranian
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_iranian
!
!  Copy the input.
!
  m2 = m
  y2 = y
!
!  Get the number of days in the month.
!
  month_length_iranian = mdays(m2)
!
!  If necessary, add 1 day for a leap year..
!
  if ( m2 == 12 .and. year_is_leap_iranian ( y2 ) ) then
    month_length_iranian = month_length_iranian + 1
  end if

  return
end
function month_length_islamic ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_ISLAMIC returns the number of days in an Islamic month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_ISLAMIC, the number of days
!    in the month.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(12) :: mdays = (/ &
    30, 29, 30, 29, 30, 29, 30, 29, 30, 29, 30, 29 /)
  integer ( kind = 4 ) month_length_islamic
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_islamic
!
!  Copy the input.
!
  m2 = m
  y2 = y
!
!  Check the input.
!
  call ym_check_islamic ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    month_length_islamic = 0
    return
  end if
!
!  Get the number of days in the month.
!
  month_length_islamic = mdays(m2)
!
!  If necessary, add 1 day for a leap year.
!
  if ( m2 == 12 .and. year_is_leap_islamic ( y2 ) ) then
    month_length_islamic = month_length_islamic + 1
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
    month_length_julian = 0
    return
  end if
!
!  Get the number of days in the month.
!
  month_length_julian = mdays(m2)
!
!  If necessary, add 1 day for February 29.
!
  if ( m2 == 2 .and. year_is_leap_julian ( y2 ) ) then
    month_length_julian = month_length_julian + 1
  end if

  return
end
function month_length_lunar ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_LUNAR returns the number of days in a lunar month.
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
!    Input, integer ( kind = 4 ) Y, the year in which the month occurred.
!
!    Input, integer ( kind = 4 ) M, the number of the month.
!
!    Output, real ( kind = 8 ) MONTH_LENGTH_LUNAR, the number of days in
!    the month.
!
  implicit none

  integer ( kind = 4 ) m
  real ( kind = 8 ) month_length_lunar
  integer ( kind = 4 ) y

  if ( m < 1 .or. 12 < m ) then
    month_length_lunar = 0
  else
    month_length_lunar = 29.53058D+00
  end if

  return
end
function month_length_persian ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_PERSIAN returns the number of days in a Persian month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2000
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_PERSIAN, the number of days
!    in the month.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(12) :: mdays = (/ &
    31, 31, 31, 31, 31, 31, 30, 30, 30, 30, 30, 29 /)
  integer ( kind = 4 ) month_length_persian
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_persian
!
!  Copy the input.
!
  m2 = m
  y2 = y
!
!  Get the number of days in the month.
!
  month_length_persian = mdays(m2)
!
!  If necessary, add 1 day for a leap year.
!
  if ( m2 == 12 .and. year_is_leap_persian ( y2 ) ) then
    month_length_persian = month_length_persian + 1
  end if

  return
end
function month_length_republican ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_REPUBLICAN returns the number of days in a Republican month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_REPUBLICAN, the number of days
!    in the month.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) month_length_republican
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_republican
!
!  Copy the input.
!
  m2 = m
  y2 = y
!
!  Check the input.
!
  call ym_check_republican ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    month_length_republican = 0
    return
  end if
!
!  Get the number of days in the month.
!
  if ( 1 <= m2 .and. m2 <= 12 ) then
    month_length_republican = 30
  else if ( m2 == 13 ) then
    if ( year_is_leap_republican ( y2 ) ) then
      month_length_republican = 6
    else
      month_length_republican = 5
    end if
  end if

  return
end
function month_length_roman ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_ROMAN returns the number of days in a Roman month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2000
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
!    Output, integer ( kind = 4 ) MONTH_LENGTH_ROMAN, the number of days
!    in the month.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(12) :: mdays = (/ &
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  integer ( kind = 4 ) month_length_roman
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_roman
!
!  Copy the input.
!
  m2 = m
  y2 = y
!
!  Check the input.
!
  call ym_check_roman ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    month_length_roman = 0
    return
  end if

  month_length_roman = mdays(m2)

  if ( m2 == 2 .and. year_is_leap_roman ( y2 ) ) then
    month_length_roman = month_length_roman + 1
  end if

  return
end
function month_length_synodic ( )

!*****************************************************************************80
!
!! MONTH_LENGTH_SYNODIC returns the mean synodic month length.
!
!  Discussion:
!
!    The synodic month is the time from one new moon to the next, that is,
!    when the moon and Sun are in conjunction.
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
!    Output, real ( kind = 8 ) MONTH_LENGTH_SYNODIC, the length of the
!    mean synodic month,
!    in days.
!
  implicit none

  real ( kind = 8 ) month_length_synodic

  month_length_synodic = 29.530588853D+00

  return
end
subroutine month_name_to_month_common ( month_name, m )

!*****************************************************************************80
!
!! MONTH_NAME_TO_MONTH_COMMON returns the month number of a Common month
!
!  Discussion:
!
!    Capitalization is ignored.  The month name has to match up to
!    the unique beginning of a month name, and the rest is ignored.
!    Here are the limits:
!
!      JAnuary
!      February
!      MARch
!      APril
!      MAY
!      JUNe
!      JULy
!      AUgust
!      September
!      October
!      November
!      December
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) MONTH_NAME, a string containing a month
!    name or abbreviation.
!
!    Output, integer ( kind = 4 ) M, the number of the month, or -1 if the name
!    could not be recognized.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 3 ) string

  string = month_name
  call s_cap ( string )

       if ( string(1:2) == 'JA' ) then
    m = 1
  else if ( string(1:1) == 'F' ) then
    m = 2
  else if ( string(1:3) == 'MAR' ) then
    m = 3
  else if ( string(1:2) == 'AP' ) then
    m = 4
  else if ( string(1:3) == 'MAY' ) then
    m = 5
  else if ( string(1:3) == 'JUN' ) then
    m = 6
  else if ( string(1:3) == 'JUL' ) then
    m = 7
  else if ( string(1:2) == 'AU' ) then
    m = 8
  else if ( string(1:1) == 'S' ) then
    m = 9
  else if ( string(1:1) == 'O' ) then
    m = 10
  else if ( string(1:1) == 'N' ) then
    m = 11
  else if ( string(1:1) == 'D' ) then
    m = 12
  else
    m = - 1
  end if

  return
end
subroutine month_to_ides_roman ( m, d )

!*****************************************************************************80
!
!! MONTH_TO_IDES_ROMAN returns the day of the ides of a Roman month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, integer ( kind = 4 ) D, the day of the ides of the month.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter, dimension(12) :: ides = (/ &
    13, 13, 15, 13, 15, 13, 15, 13, 13, 15, 13, 13 /)

  if ( m < 1 .or. 12 < m ) then
    d = -1
  else
    d = ides(m)
  end if

  return
end
subroutine month_to_month_name_bahai ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_BAHAI returns the name of a Bahai month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 10 ), parameter, dimension(20) :: name = (/ &
    'Baha      ', 'Jalal     ', 'Jamal     ', 'Azamat    ', &
    'Nur       ', 'Rahmat    ', 'Kalimat   ', 'Kamal     ', &
    'Asma      ', 'Izzat     ', 'Mashiyyat ', 'Ilm       ', &
    'Qudrat    ', 'Qawl      ', 'Masail    ', 'Sharaf    ', &
    'Sultan    ', 'Mulk      ', 'Ayyam-i-Ha', 'Ala       ' /)

  if ( m < 1 .or. 20 < m ) then

    do i = 1, len ( month_name )
      month_name(i:i) = '?'
    end do

  else

    month_name = name(m)

  end if

  return
end
subroutine month_to_month_name_common ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_COMMON returns the name of a Common month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 9 ), parameter, dimension(12) :: name = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)

  if ( m < 1 .or. 12 < m ) then

    do i = 1, len ( month_name )
      month_name(i:i) = '?'
    end do

  else

    month_name = name(m)

  end if

  return
end
subroutine month_to_month_name_common3 ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_COMMON3 returns an abbreviation of a Common month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, character ( len = 3 ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  character ( len = 3 ) month_name
  character ( len = 3 ), parameter, dimension(12) :: name = (/ &
    'Jan', 'Feb', 'Mar', 'Apr', &
    'May', 'Jun', 'Jul', 'Aug', &
    'Sep', 'Oct', 'Nov', 'Dec' /)

  if ( m < 1 .or. 12 < m ) then

    do i = 1, len ( month_name )
      month_name(i:i) = '?'
    end do

  else

    month_name = name(m)

  end if

  return
end
subroutine month_to_month_name_coptic ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_COPTIC returns the name of a Coptic month.
!
!  Discussion:
!
!    The names are closely related to the month names of the Egyptian
!    Civil calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 10 ), parameter, dimension(13) :: name = (/ &
    'Tut       ', 'Babah     ',  'Hatur     ', 'Kiyahk    ', &
    'Tubah     ', 'Amshir    ',  'Baramhat  ', 'Baramundah', &
    'Bashans   ', 'Ba''unah   ', 'Abib      ', 'Misra     ', &
    'al-Nasi   ' /)

  if ( m < 1 .or. 13 < m ) then
    month_name = '?????'
  else
    month_name = name(m)
  end if

  return
end
subroutine month_to_month_name_eg_civil ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_EG_CIVIL returns the name of an Egyptian Civil month.
!
!  Discussion:
!
!    The 13th month had only 5 days, which were treated as the birthdays
!    of Osiris, Horus, Set, Isis and Nephthys.
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
!    Input, integer ( kind = 4 ) M, the month.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 10 ), parameter, dimension(13) :: name = (/ &
    'Thoth     ', 'Phaophi   ', 'Hathyr    ', 'Choiak    ', &
    'Tybi      ', 'Mecheir   ', 'Phamenoth ', 'Pharmouthi', &
    'Pachon    ', 'Payni     ', 'Epeiph    ', 'Mesore    ', &
    'Epagomenai' /)

  if ( m < 1 .or. 13 < m ) then
    month_name = '?????'
  else
    month_name = name(m)
  end if

  return
end
subroutine month_to_month_name_eg_lunar ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_EG_LUNAR returns the name of an Egyptian Lunar month.
!
!  Discussion:
!
!    "Akhet" means "flood",
!    "Peret" means "going forth" (for planting),
!    "Shomu" means "deficiency" (the dry season).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 10 ), parameter, dimension(13) :: name = (/ &
    'Akhet I   ', 'Akhet II  ', 'Akhet III ', 'Akhet IV  ', &
    'Peret I   ', 'Peret II  ', 'Peret III ', 'Peret IV  ', &
    'Shomu I   ', 'Shomu II  ', 'Shomu III ', 'Shomu IV  ', &
    'Shomu V   ' /)

  if ( m < 1 .or. 13 < m ) then
    month_name = '?????'
  else
    month_name = name(m)
  end if

  return
end
subroutine month_to_month_name_ethiopian ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_ETHIOPIAN returns the name of an Ethiopian month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 10 ), parameter, dimension(13) :: name = (/ &
    'Maskaram  ', 'Teqemt    ', 'Khedar    ', 'Takhsas   ', &
    'Ter       ', 'Yakatit   ', 'Magabit   ', 'Miyazya   ', &
    'Genbot    ', 'Sane      ', 'Hamle     ', 'Nahase    ', &
    'Paguemen  ' /)

  if ( m < 1 .or. 13 < m ) then
    month_name = '?????'
  else
    month_name = name(m)
  end if

  return
end
subroutine month_to_month_name_greek ( y, m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_GREEK returns the name of a Greek month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, the year and month.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 12 ), parameter, dimension(13) :: name = (/ &
    'Hecatombaeon', 'Metageitnion', 'Boedromion  ', 'Pyanepsion  ', &
    'Maemacterion', 'Poseidon    ', 'Poseidon II ', 'Gamelion    ', &
    'Anthesterion', 'Elaphebolion', 'Munychion   ', 'Thargelion  ', &
    'Scirophorion' /)
  integer ( kind = 4 ) y
  logical   year_is_embolismic_greek
!
!  13 month year.
!
  if ( year_is_embolismic_greek ( y ) ) then

    if ( m < 1 .or. 13 < m ) then

      month_name = '?????'

    else

      month_name = name(m)

    end if
!
!  12 month year.
!
  else

    if ( m < 1 .or. 12 < m ) then

      month_name = '?????'

    else if ( m <= 6 ) then

      month_name = name(m)

    else

      month_name = name(m+1)

    end if

  end if

  return
end
subroutine month_to_month_name_hebrew ( y, m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_HEBREW returns the name of a Hebrew month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, the year and month.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 9 ), parameter, dimension(13) :: name = (/ &
    'Tishri   ', 'Heshvan  ', 'Kislev   ', 'Tebet    ', 'Shebat   ', &
    'Adar     ', 'Veadar   ', 'Nisan    ', 'Iyar     ', 'Sivan    ', &
    'Tammuz   ', 'Ab       ', 'Elul     ' /)
  integer ( kind = 4 ) y
  logical   year_is_embolismic_hebrew

  if ( year_is_embolismic_hebrew ( y ) ) then

    if ( m < 1 .or. 13 < m ) then

      month_name = '?????'

    else

      month_name = name(m)

    end if

  else

    if ( m < 1 .or. 12 < m ) then

      month_name = '?????'

    else if ( m <= 6 ) then

      month_name = name(m)

    else

      month_name = name(m+1)

    end if

  end if

  return
end
subroutine month_to_month_name_hindu_lunar ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_HINDU_LUNAR returns the name of a Hindu lunar month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 10 ), parameter, dimension(12) :: name = (/ &
    'Chaitra   ', 'Vaisakha  ', 'Jyaishtha ', 'Ashadha   ', &
    'Sravana   ', 'Bhadrapada', 'Asvina    ', 'Karttika  ', &
    'Margasira ', 'Pausha    ', 'Magha     ', 'Phalguna  ' /)

  if ( m < 1 .or. 12 < m ) then

    month_name = '?????'

  else

    month_name = name(m)

  end if

  return
end
subroutine month_to_month_name_hindu_solar ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_HINDU_SOLAR returns the name of a Hindu solar month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 10 ), parameter, dimension(12) :: name = (/ &
    'Mesha     ', 'Vrshabha  ', 'Mithuna   ', 'Karka     ', &
    'Simha     ', 'Kanya     ', 'Tula      ', 'Vrischika ', &
    'Dhanus    ', 'Makara    ', 'Kumbha    ', 'Mina      ' /)

  if ( m < 1 .or. 12 < m ) then

    month_name = '?????'

  else

    month_name = name(m)

  end if

  return
end
subroutine month_to_month_name_iranian ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_IRANIAN returns the name of an Iranian month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 11 ), parameter, dimension(12) :: name = (/ &
    'Farvardin  ', 'Ordibehesht', 'Xordad     ', 'Tir        ', &
    'Mordad     ', 'Shahrivar  ', 'Mehr       ', 'Aban       ', &
    'Azar       ', 'Dey        ', 'Bahman     ', 'Esfand     ' /)

  if ( m < 1 .or. 12 < m ) then

    month_name = '?????'

  else

    month_name = name(m)

  end if

  return
end
subroutine month_to_month_name_islamic ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_ISLAMIC returns the name of an Islamic month.
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
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 10 ), parameter, dimension(12) :: name = (/ &
    'Muharram  ', 'Safar     ', 'Rabi I    ', 'Rabi II   ', &
    'Jumada I  ', 'Jumada II ', 'Rajab     ', 'Shaban    ', &
    'Ramadan   ', 'Shawwal   ', 'Dhul-quda ', 'Dhul-hejji' /)

  if ( m < 1 .or. 12 < m ) then

    month_name = '?????'

  else

    month_name = name(m)

  end if

  return
end
subroutine month_to_month_name_persian ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_PERSIAN returns the name of a Persian month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 11 ), parameter, dimension(12) :: name = (/ &
    'Farvardin  ', 'Ordibehesht', 'Khordad    ', 'Tir        ', &
    'Mordad     ', 'Shahrivar  ', 'Mehr       ', 'Aban       ', &
    'Azar       ', 'Dey        ', 'Bahman     ', 'Esfand     ' /)

  if ( m < 1 .or. 12 < m ) then

    month_name = '?????'

  else

    month_name = name(m)

  end if

  return
end
subroutine month_to_month_name_republican ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_REPUBLICAN returns the name of a Republican month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 14 ), parameter, dimension(13) :: name = (/ &
    'Vendemaire    ', 'Brumaire      ', 'Frimaire      ', 'Nivose        ', &
    'Pluviose      ', 'Ventose       ', 'Germinal      ', 'Floreal       ', &
    'Prairial      ', 'Messidor      ', 'Thermidor     ', 'Fructidor     ', &
    'Sansculottides' /)

  if ( m < 1 .or. 13 < m ) then

    month_name = '?????'

  else

    month_name = name(m)

  end if

  return
end
subroutine month_to_month_name_roman ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_ROMAN returns the name of a Roman month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 10 ), parameter, dimension(12) :: name = (/ &
    'Januarius ', 'Februarius', 'Martius   ', 'Aprilis   ', &
    'Maius     ', 'Junius    ', 'Julius    ', 'Augustus  ', &
    'September ', 'October   ', 'November  ', 'December  ' /)

  if ( m < 1 .or. 12 < m ) then

    do i = 1, len ( month_name )
      month_name(i:i) = '?'
    end do

  else

    month_name = name(m)

  end if

  return
end
subroutine month_to_month_name_soor_san ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_SOOR_SAN returns the name of a Soor San month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 7 ), parameter, dimension(12) :: name = (/ &
    'Baune  ', 'Abib   ', 'Meshri ', 'Tot    ', &
    'Babe   ', 'Hatur  ', 'Kyak   ', 'Tabe   ', &
    'Mashir ', 'Buramat', 'Barsude', 'Bashans' /)

  if ( m < 1 .or. 12 < m ) then

    month_name = '?????'

  else

    month_name = name(m)

  end if

  return
end
subroutine month_to_month_name_zoroastrian ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME_ZOROASTRIAN returns the name of a Zoroastrian month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, character ( len = * ) MONTH_NAME, the month name.
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 11 ), parameter, dimension(12) :: name = (/ &
    'Furvurdeen ', 'Ardibehest ', 'Khordad    ', 'Tir        ', &
    'Amerdad    ', 'Sherever   ', 'Moher      ', 'Aban       ', &
    'Adur       ', 'Deh        ', 'Bahman     ', 'Aspendadmad' /)

  if ( m < 1 .or. 12 < m ) then

    month_name = '?????'

  else

    month_name = name(m)

  end if

  return
end
subroutine month_to_nones_roman ( m, d )

!*****************************************************************************80
!
!! MONTH_TO_NONES_ROMAN returns the day of the nones of a Roman month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the month index.
!
!    Output, integer ( kind = 4 ) D, the day of the nones of the month.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter, dimension(12) :: nones = (/ &
    5, 5, 7, 5, 7, 5, 7, 5, 5, 7, 5, 5 /)

  if ( m < 1 .or. 12 < m ) then
    d = -1
  else
    d = nones(m)
  end if

  return
end
subroutine moon_phase_to_jed ( n, phase, jed )

!*****************************************************************************80
!
!! MOON_PHASE_TO_JED calculates the JED of a moon phase.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 June 2001
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes: The Art of Scientific Computing,
!    Cambridge University Press.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, specifies that the N-th such phase
!    of the moon since January 1900 is to be computed.
!
!    Input, integer ( kind = 4 ) PHASE, specifies which phase is to be computed.
!    0=new moon,
!    1=first quarter,
!    2=full,
!    3=last quarter.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date on which the
!    requested phase occurs.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: degrees_to_radians = pi / 180.0D+00

  real ( kind = 8 ) am
  real ( kind = 8 ) as
  real ( kind = 8 ) c
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) n
  integer ( kind = 4 ) phase
  real ( kind = 8 ) t
  real ( kind = 8 ) xtra
!
!  First estimate.
!
  j = 2415020 + 28 * n + 7 * phase
!
!  Compute a correction term.
!
  c = n + phase / 4.0D+00

  t = c / 1236.85D+00

  xtra = 0.75933D+00 + 1.53058868D+00 * c &
    + ( 0.0001178D+00 - 0.000000155D+00 * t ) * t * t

  as = degrees_to_radians * ( 359.2242D+00 + 29.105356D+00 * c )

  am = degrees_to_radians * &
    ( 306.0253D+00 + 385.816918D+00 * c + 0.010730D+00 * t * t )

  if ( phase == 0 .or. phase == 2 ) then

    xtra = xtra + ( 0.1734D+00 - 0.000393D+00 * t ) * sin ( as ) &
      - 0.4068D+00 * sin ( am )

  else if ( phase == 1 .or. phase == 3 ) then

    xtra = xtra + ( 0.1721D+00 - 0.0004D+00 * t ) * sin ( as ) &
      - 0.6280D+00 * sin ( am )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MOON_PHASE_TO_JED - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal PHASE option = ', phase
    stop

  end if

  jed = real ( j, kind = 8 ) + xtra

  return
end
subroutine mothers_day ( y, m, d )

!*****************************************************************************80
!
!! MOTHERS_DAY computes the date of Mother's Day (US) for a Common year.
!
!  Discussion:
!
!    Mother's Day occurs on the second Sunday in May.
!
!  Example:
!
!    Input:
!
!      Y = 2003
!
!    Output:
!
!      M = 5
!      D = 11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, integer ( kind = 4 ) M, D, the month and day of Mother's Day.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) m
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
!
!  Determine the day of the week for 8 May, the earliest day
!  that Mother's day can occur.
!
  m = 5
  d = 8
  f = 0.0D+00

  call ymdf_to_weekday_common ( y, m, d, f, w )
!
!  W = 1 means this day is Sunday, and day D is Mother's day.
!  Otherwise, figure out how to increment W to 8 (Sunday again);
!  The same increment makes D the correct day number.
!
  if ( w /= 1 ) then
    d = d + 8 - w
  end if

  return
end
subroutine new_year_to_jed_hebrew ( y, jed )

!*****************************************************************************80
!
!! NEW_YEAR_TO_JED_HEBREW returns the JED of the beginning of a Hebrew year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Algorithm G,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, pages 330.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the Hebrew year.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) e_prime
  integer ( kind = 4 ) i4_wrap
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) t_prime
  integer ( kind = 4 ) tc
  integer ( kind = 4 ) th
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
!
!  integer ( kind = 4 ) option.
!
!     t = 31524 + 765433 * ( ( 235 * y - 234 ) / 19 )
!     d = t / 25920
!     t_prime = mod ( t, 25920 )
!
!  Real option.
!
  mu = ( 235 * y - 234 ) / 19
  tc = 204 + 793 * mu
  th = 5 + 12 * mu + tc / 1080
  d = 1 + 29 * mu + th / 24
  t_prime = mod ( tc, 1080 ) + 1080 * mod ( th, 24 )

  w = i4_wrap ( d+1, 1, 7 )

  e = mod ( 7 * y + 13, 19 ) / 12
  e_prime = mod ( 7 * y + 6, 19 ) / 12

  if ( 19940 <= t_prime .or. &
     (  9924 <= t_prime .and. w == 3 .and. e == 0 ) .or. &
     ( 16788 <= t_prime .and. w == 2 .and. e == 0 .and. e_prime == 1 ) ) then
    d = d + 1
  end if

  call epoch_to_jed_hebrew ( jed_epoch )

  jed = jed_epoch - 1 + real ( d + mod ( mod ( d + 5, 7 ), 2 ), kind = 8 )

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
subroutine now_to_yjf_common ( y, j, f )

!*****************************************************************************80
!
!! NOW_TO_YJF_COMMON expresses the current date as a Common YJF date.
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
!    Output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
  implicit none

  character ( len = 8 ) date
  integer ( kind = 4 ) d1
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) h1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) mu1
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) s1
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y1 = values(1)
  m1 = values(2)
  d1 = values(3)
  h1 = values(5)
  n1 = values(6)
  s1 = values(7)
  mu1 = values(8)

  f1 = real ( mu1, kind = 8 )
  f1 = real ( s1, kind = 8 ) + f1 / 1000.0D+00
  f1 = real ( n1, kind = 8 ) + f1 / 60.0D+00
  f1 = real ( h1, kind = 8 ) + f1 / 60.0D+00
  f1 = f1 / 24.0D+00

  call ymdf_to_yjf_common ( y1, m1, d1, f1, y, j, f )

  return
end
subroutine now_to_ymdf_common ( y, m, d, f )

!*****************************************************************************80
!
!! NOW_TO_YMDF_COMMON expresses the current date as a Common YMDF date.
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
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  character ( len = 8 ) date
  real ( kind = 8 ) f
  integer ( kind = 4 ) h
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

  return
end
subroutine now_to_ymdhms_common ( y, m, d, h, n, s )

!*****************************************************************************80
!
!! NOW_TO_YMDHMS_COMMON expresses the current date as a Common YMDHMS date.
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
!    Output, integer ( kind = 4 ) Y, M, D, H, N, S, the YMDHMS date.
!
  implicit none

  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
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
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP switches two R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r8_to_s_left ( d, s )

!*****************************************************************************80
!
!! R8_TO_S_LEFT writes an R8 into a left justified string.
!
!  Method:
!
!    A 'G20.12' format is used with a WRITE statement.
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
!    Input, real ( kind = 8 ) D, the number to be written into the string.
!
!    Output, character ( len = * ) S, the string into which
!    the real number is to be written.  If the string is less than 20
!    characters long, it will will be returned as a series of
!    asterisks.
!
  implicit none

  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character ( len = 20 ) s2

  nchar = len ( s )

  if ( nchar < 20 ) then

    do i = 1, nchar
      s(i:i) = '*'
    end do

  else if ( d == 0.0D+00 ) then
    s(1:20) = '     0.0      '
  else
    write ( s2, '(g20.12)' ) d
    s(1:20) = s2
  end if
!
!  Shift the string left.
!
  s = adjustl ( s )

  return
end
function r8_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! R8_UNIFORM_AB returns a scaled pseudorandom R8.
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
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_AB, a number strictly between A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed <= 0 ) then
    seed = seed + 2147483647
  end if

  r8_uniform_ab = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine random_initialize ( seed )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN 90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed ( )
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!    And this is the FORTRAN 90 people's idea of convenience?
!
!    And I still get poorly randomized values, somehow, having to do
!    with a bad seed, or something.  I am about ready to go back to
!    using my own damn routine!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a seed value.
!
  implicit none

  integer ( kind = 4 ) date_time(8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
  real t
!
!  Initialize the random number seed.
!
  call random_seed ( )
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )
!
!  Get the current date and time.
!
  call date_and_time ( values = date_time )
!
!  Construct a slightly random value.
!
  seed = 0
  do i = 1, 8
    seed = ieor ( seed, date_time(i) )
  end do
!
!  Make slightly random assignments to SEED_VECTOR.
!
  do i = 1, seed_size
    seed_vector(i) = ieor ( seed, i )
  end do
!
!  Set the random number seed value.
!
  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Because EVEN THIS DOESN'T SEEM TO PROPERLY MIX UP THE RANDOM
!  NUMBERS, call the random number routine a bunch of times.
!
  do i = 1, 100
    call random_number ( harvest = t )
  end do
!
!  I STILL GET LOUSY RESULTS.  THE HELL WITH IT!
!
  return
end
subroutine rd_to_jed ( rd, jed )

!*****************************************************************************80
!
!! RD_TO_JED converts an RD to a JED.
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
!    Input, real ( kind = 8 ) RD, the RD Date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  real ( kind = 8 ) jed
  real ( kind = 8 ) rd
  real ( kind = 8 ) rd_epoch

  call epoch_to_jed_rd ( rd_epoch )
  jed = rd_epoch + rd

  return
end
subroutine s_cap ( s )

!*****************************************************************************80
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) s

  nchar = len ( s )

  do i = 1, nchar

    c = s(i:i)
    call ch_cap ( c )
    s(i:i) = c

  end do

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
subroutine s_cat1 ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT1 concatenates two strings, with a single blank separator.
!
!  Example:
!
!    S1       S2       S
!
!    'cat'    'dog'    'cat dog'
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
!    Output, character ( len = * ) S3, the string made by concatenating
!    S1, a blank, and S2, ignoring any trailing blanks.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  s3 = trim ( s1 ) // ' ' // trim ( s2 )

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
subroutine s_to_hms ( s, pat, h, n, second )

!*****************************************************************************80
!
!! S_TO_HMS converts a string into a H:M:S date.
!
!  Discussion:
!
!    The characters in PAT indicate where the data is stored.  A particular
!    letter, such as "H", indicates, an hour field, while the number of "H"
!    characters indicates the width of the field.
!
!    The codes are:
!
!    'H' or 'h' an hour field
!    'M' or 'm' a minute field
!    'S' or 's' a second field
!
!  Example:
!
!    S                    PAT
!    ------------         ------------
!    '230859'             'hhmmss'
!    '10:30'              'hh:mm'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string containing the data.
!
!    Input, character ( len = * ) PAT, describes how the data is stored.
!    PAT must be the same length as S.
!
!    Output, integer ( kind = 4 ) H, N, SECOND, the hour, minute and second
!    represented by the string.  Any item not read from the string will
!    have a value of -1.
!
  implicit none

  integer ( kind = 4 ) h
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  integer ( kind = 4 ) n
  character ( len = * ) pat
  character ( len = * ) s
  integer ( kind = 4 ) second

  h = 0
  n = 0
  second = 0

  length = min ( len ( s ), len ( pat ) )

  ihi = 0

  do while ( ihi < length )

    ilo = ihi + 1
    ihi = ilo

    do while ( ihi + 1 <= length )
      if ( pat(ihi+1:ihi+1) /= pat(ilo:ilo) ) then
        exit
      end if
      ihi = ihi + 1
    end do

         if ( pat(ilo:ilo) == 'H' .or. pat(ilo:ilo) == 'h' ) then
      call s_to_i4 ( s(ilo:ihi), h, ierror, last )
    else if ( pat(ilo:ilo) == 'M' .or. pat(ilo:ilo) == 'm' ) then
      call s_to_i4 ( s(ilo:ihi), n, ierror, last )
    else if ( pat(ilo:ilo) == 'S' .or. pat(ilo:ilo) == 's' ) then
      call s_to_i4 ( s(ilo:ihi), second, ierror, last )
    end if

  end do

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If STRING is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
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
subroutine s_to_ymd_common ( s, pat, y, m, d )

!*****************************************************************************80
!
!! S_TO_YMD_COMMON converts a string into a Common YMD date.
!
!  Discussion:
!
!    The characters in PAT indicate where the day, month and year data
!    is stored.  A particular letter, such as "Y", indicates, a year
!    field, while the number of "Y" characters indicates the width of
!    the field.
!
!    The codes are:
!
!    'Y' or 'y', a year field
!    'M' or 'm', a numeric month field
!    'N' or 'n', a literal month field
!    'D' or 'd', a day field
!
!  Example:
!
!    S              PAT
!    ------------   ------------
!    '19991031'     'YYYYMMDD'
!    '10-31-99'     'MM-DD-YY'
!    '10-31-99'     'MM/DD/YY'
!    'Oct 31 1999'  'NNN DD YYYY'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string containing the data.
!
!    Input, character ( len = * ) PAT, describes how the data is stored.
!    PAT must be the same length as S.
!
!    Output, integer ( kind = 4 ) Y, M, D, the YMD date
!    represented by the string.  Any item not read from the string will
!    have a value of -1.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  integer ( kind = 4 ) m
  character ( len = * ) pat
  character ( len = * ) s
  integer ( kind = 4 ) y

  d = 0
  m = 0
  y = 0

  length = min ( len ( s ), len ( pat ) )

  ihi = 0

  do while ( ihi < length )

    ilo = ihi + 1
    ihi = ilo

    do while ( ihi + 1 <= length )
      if ( pat(ihi+1:ihi+1) /= pat(ilo:ilo) ) then
        exit
      end if
      ihi = ihi + 1
    end do

         if ( pat(ilo:ilo) == 'Y' .or. pat(ilo:ilo) == 'y' ) then
      call s_to_i4 ( s(ilo:ihi), y, ierror, last )
    else if ( pat(ilo:ilo) == 'M' .or. pat(ilo:ilo) == 'm' ) then
      call s_to_i4 ( s(ilo:ihi), m, ierror, last )
    else if ( pat(ilo:ilo) == 'N' .or. pat(ilo:ilo) == 'n' ) then
      call month_name_to_month_common ( s(ilo:ihi), m )
    else if ( pat(ilo:ilo) == 'D' .or. pat(ilo:ilo) == 'd' ) then
      call s_to_i4 ( s(ilo:ihi), d, ierror, last )
    end if

  end do

  return
end
subroutine s_to_ymdhms_common ( s, pat, y, m, d, h, n, second )

!*****************************************************************************80
!
!! S_TO_YMDHMS_COMMON converts a string into a Common YMD H:M:S date.
!
!  Discussion:
!
!    The characters in PAT indicate where the data is stored.  A particular
!    letter, such as "Y", indicates, a year field, while the number of "Y"
!    characters indicates the width of the field.
!
!    The codes are:
!
!    'Y' a year field
!    'M' a numeric month field
!    'N' a literal month field
!    'D' a day field
!    'h' an hour field
!    'm' a minute field
!    's' a second field
!
!  Example:
!
!    S                    PAT
!    ------------         ------------
!    '19991031230859'     'YYYYMMDDhhmmss'
!    '10-31-99'           'MM-DD-YY'
!    '10-31-99'           'MM/DD/YY'
!    'Oct 31 1999'        'NNN DD YYYY'
!    '10:30'              'hh:mm'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string containing the data.
!
!    Input, character ( len = * ) PAT, describes how the data is stored.
!    PAT must be the same length as S.
!
!    Output, integer ( kind = 4 ) Y, M, D, HOUR, N, SECOND, the YMDHMS
!    date represented by the string.  Any item not
!    read from the string will have a value of -1.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = * ) pat
  character ( len = * ) s
  integer ( kind = 4 ) second
  integer ( kind = 4 ) y

  y = 0
  m = 0
  d = 0
  h = 0
  n = 0
  second = 0

  length = min ( len ( s ), len ( pat ) )

  ihi = 0

  do while ( ihi < length )

    ilo = ihi + 1
    ihi = ilo

    do while ( ihi + 1 <= length )
      if ( pat(ihi+1:ihi+1) /= pat(ilo:ilo) ) then
        exit
      end if
      ihi = ihi + 1
    end do

         if ( pat(ilo:ilo) == 'Y' ) then
      call s_to_i4 ( s(ilo:ihi), y, ierror, last )
    else if ( pat(ilo:ilo) == 'M' ) then
      call s_to_i4 ( s(ilo:ihi), m, ierror, last )
    else if ( pat(ilo:ilo) == 'N' ) then
      call month_name_to_month_common ( s(ilo:ihi), m )
    else if ( pat(ilo:ilo) == 'D' ) then
      call s_to_i4 ( s(ilo:ihi), d, ierror, last )
    else if ( pat(ilo:ilo) == 'h' ) then
      call s_to_i4 ( s(ilo:ihi), h, ierror, last )
    else if ( pat(ilo:ilo) == 'm' ) then
      call s_to_i4 ( s(ilo:ihi), n, ierror, last )
    else if ( pat(ilo:ilo) == 's' ) then
      call s_to_i4 ( s(ilo:ihi), second, ierror, last )
    end if

  end do

  return
end
subroutine second_borrow_common ( y, m, d, h, n, s )

!*****************************************************************************80
!
!! SECOND_BORROW_COMMON "borrows" a minute of seconds in a common date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, H, N, S, the YMDHMS date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) y

  do while ( s < 0 )

    s = s + 60
    n = n - 1

    call minute_borrow_common ( y, m, d, h, n )

  end do

  return
end
subroutine second_carry_common ( y, m, d, h, n, s )

!*****************************************************************************80
!
!! SECOND_CARRY_COMMON: given a Common YMDHMS date, carries seconds to minutes.
!
!  Algorithm:
!
!    While 60 <= S:
!
!      decrease S by 60;
!      increase N by 1;
!      if necessary, adjust H, D, M and Y.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, H, N, S,
!    the year, month, day, hours, minutes, seconds,
!    On output, S is between 0 and 59.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) y

  do while ( 60 <= s )

    s = s - 60
    n = n + 1

    call minute_carry_common ( y, m, d, h, n )

  end do

  return
end
subroutine ss_to_jed_unix ( s, jed )

!*****************************************************************************80
!
!! SS_TO_JED_UNIX converts a UNIX SS date to a JED.
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
!    Input, real ( kind = 8 ) S, the UNIX date.
!
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) s

  call epoch_to_jed_unix ( jed_epoch )

  d = s / ( 24.0D+00 * 60.0D+00 * 60.0D+00 )

  jed = jed_epoch + d

  return
end
subroutine thanksgiving_canada ( y, m, d )

!*****************************************************************************80
!
!! THANKSGIVING_CANADA computes Canadian Thanksgiving for a Common year.
!
!  Discussion:
!
!    Canadian Thanksgiving occurs on the second Monday in October.
!
!  Example:
!
!    Input:
!
!      Y = 2002
!
!    Output:
!
!      M = 11
!      D = 28
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, integer ( kind = 4 ) M, D, the month and day of Thanksgiving.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) m
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
!
!  Determine the day of the week for 8 October, the earliest day
!  that Thanksgiving can occur.
!
  m = 10
  d = 8
  f = 0.0E+00

  call ymdf_to_weekday_common ( y, m, d, f, w )
!
!  If W = 2 means this day is Monday, and day D is Thanksgiving.
!  Otherwise, figure out how to increment W to 2;
!  The same increment makes D the correct day number.
!
  if ( w < 2 ) then
    d = d + 2 - w
  else if ( 2 < w ) then
    d = d + 2 + 7 - w
  end if

  return
end
subroutine thanksgiving_us ( y, m, d )

!*****************************************************************************80
!
!! THANKSGIVING_US computes the date of Thanksgiving (US) for a Common year.
!
!  Discussion:
!
!    Thanksgiving (US) occurs on the fourth Thursday in November.
!
!  Example:
!
!    Input:
!
!      Y = 2002
!
!    Output:
!
!      M = 11
!      D = 28
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, integer ( kind = 4 ) M, D, the month and day of Thanksgiving.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) m
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
!
!  Determine the day of the week for 22 November, the earliest day
!  that Thanksgiving can occur.
!
  m = 11
  d = 22
  f = 0.0E+00

  call ymdf_to_weekday_common ( y, m, d, f, w )
!
!  W = 5 means this day is Thursday, and day D is Thanksgiving.
!  Otherwise, figure out how to increment W to 5;
!  The same increment makes D the correct day number.
!
  if ( w < 5 ) then
    d = d + 5 - w
  else if ( 5 < w ) then
    d = d + 12 - w
  end if

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
subroutine transition_to_jed_english ( jed )

!*****************************************************************************80
!
!! TRANSITION_TO_JED_ENGLISH returns the English calendar transition as a JED.
!
!  Discussion:
!
!    In the English calendar, the last moment of the Julian calendar was
!      11:59 pm, 2 September 1752 Julian/English,
!      11:59 pm, 13 September 1752 Gregorian/CE.
!    The first minute of the Gregorian calendar ended at
!      12:01 am, 3 September 1752 Julian,
!      12:01 am, 15 September 1752 Gregorian/CE/English.
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

  jed = 2361221.5D+00

  return
end
subroutine transition_to_jed_jed ( jed )

!*****************************************************************************80
!
!! TRANSITION_TO_JED_JED returns the JED calendar transition as a JED.
!
!  Discussion:
!
!    In Scaliger's design of the JED, three cycles with different periods
!    began on JED = 0.  These three cycles coincide once more on the
!    transition day.
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

  jed = 2913943.0D+00

  return
end
subroutine transition_to_jed_mayan_long ( jed )

!*****************************************************************************80
!
!! TRANSITION_TO_JED_MAYAN_LONG: Mayan long count calendar transition as a JED.
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

  jed = 2456284.5D+00

  return
end
subroutine weekday_check_common ( w )

!*****************************************************************************80
!
!! WEEKDAY_CHECK_COMMON makes sure the Common weekday number is between 1 and 7.
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
!    Input/output, integer ( kind = 4 ) W, the weekday index.
!
  implicit none

  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) w

  w = i4_wrap ( w, 1, 7 )

  return
end
subroutine weekday_to_name_bahai ( w, s )

!*****************************************************************************80
!
!! WEEKDAY_TO_NAME_BAHAI returns the name of a Bahai weekday.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2000
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

  character ( len = 8 ), parameter, dimension(7) :: weekday_name = (/ &
    'Jalal   ', 'Jamal   ', 'Kamal   ', 'Fidal   ', &
    'Idal    ', 'Istijlal', 'Istiqlal' /)
  character ( len = * ) s
  integer ( kind = 4 ) w
  integer ( kind = 4 ) w2
!
!  Make a local copy of the weekday number.
!
  w2 = w
!
!  Check the weekday number.
!
  call weekday_check_common ( w2 )
!
!  Return the weekday name.
!
  s = weekday_name(w2)

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

  character ( len = 9 ), parameter, dimension(7) :: weekday_name = (/ &
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
!  Check the weekday number.
!
  call weekday_check_common ( w2 )
!
!  Return the weekday name.
!
  s = weekday_name(w2)

  return
end
subroutine weekday_to_name_common2 ( w, s )

!*****************************************************************************80
!
!! WEEKDAY_TO_NAME_COMMON2 returns the abbreviated name of a Common weekday.
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
!    Output, character ( len = 2 ) S, the abbreviated weekday name.
!
  implicit none

  character ( len = 2 ), parameter, dimension(7) :: weekday_name = (/ &
    'Su', ' M', 'Tu', ' W', 'Th', ' F', 'Sa' /)
  character ( len = * ) s
  integer ( kind = 4 ) w
  integer ( kind = 4 ) w2
!
!  Make a local copy of the weekday number.
!
  w2 = w
!
!  Check the weekday number.
!
  call weekday_check_common ( w2 )
!
!  Return the weekday name.
!
  s = weekday_name(w2)

  return
end
subroutine weekday_to_name_common3 ( w, s )

!*****************************************************************************80
!
!! WEEKDAY_TO_NAME_COMMON3 returns the abbreviated name of a Common weekday.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) W, the weekday index.
!
!    Output, character ( len = 3 ) S, the abbreviated weekday name.
!
  implicit none

  character ( len = 3 ), parameter, dimension(7) :: weekday_name = (/ &
    'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat' /)
  character ( len = * ) s
  integer ( kind = 4 ) w
  integer ( kind = 4 ) w2
!
!  Make a local copy of the weekday number.
!
  w2 = w
!
!  Check the weekday number.
!
  call weekday_check_common ( w2 )
!
!  Return the weekday name.
!
  s = weekday_name(w2)

  return
end
subroutine weekday_to_name_french ( w, s )

!*****************************************************************************80
!
!! WEEKDAY_TO_NAME_FRENCH returns the name of a French weekday.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 2012
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

  character ( len = 8 ), parameter, dimension(7) :: weekday_name = (/ &
    'Dimanche', 'Lundi   ', 'Mardi   ', 'Mercredi', 'Jeudi   ', &
    'Vendredi', 'Samedi  ' /)
  character ( len = * ) s
  integer ( kind = 4 ) w
  integer ( kind = 4 ) w2
!
!  Make a local copy of the weekday number.
!
  w2 = w
!
!  Check the weekday number.
!
  call weekday_check_common ( w2 )
!
!  Return the weekday name.
!
  s = weekday_name(w2)

  return
end
subroutine weekday_to_name_german ( w, s )

!*****************************************************************************80
!
!! WEEKDAY_TO_NAME_GERMAN returns the name of a German weekday.
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

  character ( len = 10 ), parameter, dimension(7) :: weekday_name = (/ &
    'Sonntag   ', 'Montag    ', 'Dienstag  ', 'Mittwoch  ', &
    'Donnerstag', 'Freitag   ', 'Samstag   ' /)
  character ( len = * ) s
  integer ( kind = 4 ) w
  integer ( kind = 4 ) w2
!
!  Make a local copy of the weekday number.
!
  w2 = w
!
!  Check the weekday number.
!
  call weekday_check_common ( w2 )
!
!  Return the weekday name.
!
  s = weekday_name(w2)

  return
end
subroutine weekday_to_name_hebrew ( w, s )

!*****************************************************************************80
!
!! WEEKDAY_TO_NAME_HEBREW returns the name of a Hebrew weekday.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 2000
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

  character ( len = 12 ), parameter, dimension(7) :: weekday_name = (/ &
    'Yom rishon  ', 'Yom sheni   ', 'Yom shelishi', 'Yom revii   ', &
    'Yom hamishi ', 'Yom shishi  ', 'Sabbath     ' /)
  character ( len = * ) s
  integer ( kind = 4 ) w
  integer ( kind = 4 ) w2
!
!  Make a local copy of the weekday number.
!
  w2 = w
!
!  Check the weekday number.
!
  call weekday_check_common ( w2 )
!
!  Return the weekday name.
!
  s = weekday_name(w2)

  return
end
subroutine weekday_to_name_islamic ( w, s )

!*****************************************************************************80
!
!! WEEKDAY_TO_NAME_ISLAMIC returns the name of an Islamic weekday.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
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

  character ( len = 13 ), parameter, dimension(7) :: weekday_name = (/ &
    'Yom ilHadd   ', 'Yom litneen  ', 'Yom ittalat  ', 'Yom larba    ', &
    'Yom ilkhamiis', 'Yom ilguma   ', 'Yom issabt   ' /)
  character ( len = * ) s
  integer ( kind = 4 ) w
  integer ( kind = 4 ) w2
!
!  Make a local copy of the weekday number.
!
  w2 = w
!
!  Check the weekday number.
!
  call weekday_check_common ( w2 )
!
!  Return the weekday name.
!
  s = weekday_name(w2)

  return
end
subroutine weekday_to_name_italian ( w, s )

!*****************************************************************************80
!
!! WEEKDAY_TO_NAME_ITALIAN returns the name of an Italian weekday.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 2012
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

  character ( len = 9 ), parameter, dimension(7) :: weekday_name = (/ &
    'domenica ', 'lunedi   ', 'martedi  ', 'mercoledi', 'giovedi  ', &
    'venerdi  ', 'sabato   ' /)
  character ( len = * ) s
  integer ( kind = 4 ) w
  integer ( kind = 4 ) w2
!
!  Make a local copy of the weekday number.
!
  w2 = w
!
!  Check the weekday number.
!
  call weekday_check_common ( w2 )
!
!  Return the weekday name.
!
  s = weekday_name(w2)

  return
end
subroutine weekday_to_name_republican ( w, s )

!*****************************************************************************80
!
!! WEEKDAY_TO_NAME_REPUBLICAN returns the name of a Republican weekday.
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
!    Input, integer ( kind = 4 ) W, the weekday index.
!
!    Output, character ( len = * ) S, the weekday name.
!
  implicit none

  character ( len = 9 ), parameter, dimension(10) :: weekday_name = (/ &
    'Primedi ', 'Duodi   ', 'Tridi   ', 'Quartidi', 'Quintidi', &
    'Sextidi ', 'Septidi ', 'Octidi  ', 'Nonidi  ', 'Decadi  ' /)
  character ( len = * ) s
  integer ( kind = 4 ) w

  if ( w < 1 .or. 10 < w ) then
    s = '?????'
  else
    s = weekday_name(w)
  end if

  return
end
subroutine weekday_to_name_roman ( w, s )

!*****************************************************************************80
!
!! WEEKDAY_TO_NAME_ROMAN returns the name of a Roman weekday.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
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

  character ( len = 13 ), parameter, dimension(7) :: weekday_name = (/ &
    'Dies Solis   ', 'Dies Lunae   ', 'Dies Martis  ', 'Dies Mercurii', &
    'Dies Iovis   ', 'Dies Veneris ', 'Dies Saturni ' /)
  character ( len = * ) s
  integer ( kind = 4 ) w
  integer ( kind = 4 ) w2
!
!  Make a local copy of the weekday number.
!
  w2 = w
!
!  Check the weekday number.
!
  call weekday_check_common ( w2 )
!
!  Return the weekday name.
!
  s = weekday_name(w2)

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
subroutine y_check_alexandrian ( y, ierror )

!*****************************************************************************80
!
!! Y_CHECK_ALEXANDRIAN checks an Alexandrian year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year, which must be positive.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if Y is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y

  if ( 0 < y ) then
    ierror = 0
  else
    ierror = 1
  end if

  return
end
subroutine y_check_bahai ( y, ierror )

!*****************************************************************************80
!
!! Y_CHECK_BAHAI checks a Bahai year.
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
!    Input, integer ( kind = 4 ) Y, the year, which must not be 0.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if Y is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y

  if ( y <= 0 ) then
    ierror = 1
  else
    ierror = 0
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
subroutine y_check_eg_civil ( y, ierror )

!*****************************************************************************80
!
!! Y_CHECK_EG_CIVIL checks an Egyptian Civil year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year, which must be positive.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if Y is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y

  if ( 0 < y ) then
    ierror = 0
  else
    ierror = 1
  end if

  return
end
subroutine y_check_english ( y, ierror )

!*****************************************************************************80
!
!! Y_CHECK_ENGLISH checks an English year.
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
subroutine y_check_greek ( y, ierror )

!*****************************************************************************80
!
!! Y_CHECK_GREEK checks a Greek year.
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
!    Input, integer ( kind = 4 ) Y, the year, which must not be 0.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if Y is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y

  if ( y <= 0 ) then
    ierror = 1
  else
    ierror = 0
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
subroutine y_check_hebrew ( y, ierror )

!*****************************************************************************80
!
!! Y_CHECK_HEBREW checks a Hebrew year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year, which must be positive.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if Y is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y

  if ( 0 < y ) then
    ierror = 0
  else
    ierror = 1
  end if

  return
end
subroutine y_check_islamic ( y, ierror )

!*****************************************************************************80
!
!! Y_CHECK_ISLAMIC checks an Islamic year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year, which must be positive.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if Y is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y

  if ( 0 < y ) then
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
subroutine y_check_republican ( y, ierror )

!*****************************************************************************80
!
!! Y_CHECK_REPUBLICAN checks a Republican year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year, which must be positive.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if Y is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y

  if ( 0 < y ) then
    ierror = 0
  else
    ierror = 1
  end if

  return
end
subroutine y_check_roman ( y, ierror )

!*****************************************************************************80
!
!! Y_CHECK_ROMAN checks a Roman year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year, which must be positive.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if Y is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y

  if ( 0 < y ) then
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
subroutine y_julian_to_roman ( y, y2 )

!*****************************************************************************80
!
!! Y_JULIAN_TO_ROMAN converts a Julian year to a Roman year.
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
!    Input, integer ( kind = 4 ) Y, the Julian year.
!
!    Output, integer ( kind = 4 ) Y2, the corresponding Roman year.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2

  call y_check_julian ( y, ierror )

  if ( ierror /= 0 ) then
    y2 = -1
    return
  end if

  if ( y < 0 ) then
    y = y + 1
  end if

  y2 = y + 753

  return
end
subroutine y_roman_to_julian ( y, y2 )

!*****************************************************************************80
!
!! Y_ROMAN_TO_JULIAN converts a Roman year to a Julian year.
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
!    Input, integer ( kind = 4 ) Y, the Roman year.
!
!    Output, integer ( kind = 4 ) Y2, the corresponding Julian year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2

  y2 = y - 753

  if ( y2 <= 0 ) then
    y2 = y2 - 1
  end if

  return
end
subroutine y_to_s ( y, s )

!*****************************************************************************80
!
!! Y_TO_S writes a year into a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  character ( len = * ) s
  integer ( kind = 4 ) y

  call i4_to_s_left ( y, s )

  return
end
subroutine y_to_s_alexandrian ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_ALEXANDRIAN writes an Alexandrian year into a string.
!
!  Format:
!
!    AX YearNumber
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  character ( len = * ) s
  integer ( kind = 4 ) y

  s(1:3) = 'AX '
  call i4_to_s_left ( y, s(4:) )

  return
end
subroutine y_to_s_bahai ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_BAHAI writes a Bahai year into a string.
!
!  Format:
!
!    Bahai YearNumber
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
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = * ) s
  integer ( kind = 4 ) y

  call y_check_bahai ( y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Y_TO_S_BAHAI - Fatal error!'
    write ( *, '(a)' ) '  Illegal year.'
    stop
  end if

  s(1:6) = 'Bahai '
  call i4_to_s_left ( y, s(7:) )

  return
end
subroutine y_to_s_common ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_COMMON writes a Common year into a string.
!
!  Format:
!
!    YearNumber BCE
!    YearNumber CE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = * ) s
  integer ( kind = 4 ) y

  call y_check_common ( y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Y_TO_S_COMMON - Fatal error!'
    write ( *, '(a)' ) '  Illegal year.'
    stop
  end if

  if ( y < 0 ) then
    s(1:4) = 'BCE '
    call i4_to_s_left ( -y, s(5:) )
  else if ( 0 < y ) then
    s(1:3) = 'CE '
    call i4_to_s_left ( y, s(4:) )
  end if

  return
end
subroutine y_to_s_coptic ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_COPTIC writes a Coptic year into a string.
!
!  Format:
!
!    Coptic YearNumber
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  character ( len = * ) s
  integer ( kind = 4 ) y

  s(1:7) = 'Coptic '
  call i4_to_s_left ( y, s(8:) )

  return
end
subroutine y_to_s_eg_civil ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_EG_CIVIL writes an Egyptian Civil year into a string.
!
!  Format:
!
!    EN YearNumber
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Discussion:
!
!    "EN" stands for the Era of Nabonassar, a Babylonian king who
!    acceded in 747 BC, used by the astronomer Ptolemy to assign
!    an artificial starting year for the Egyptian calendar.
!
!  Modified:
!
!    14 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  character ( len = * ) s
  integer ( kind = 4 ) y

  s(1:3) = 'EN '
  call i4_to_s_left ( y, s(4:) )

  return
end
subroutine y_to_s_eg_lunar ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_EG_LUNAR writes an Egyptian Lunar year into a string.
!
!  Format:
!
!    EL YearNumber
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  character ( len = * ) s
  integer ( kind = 4 ) y

  s(1:3) = 'EL '
  call i4_to_s_left ( y, s(4:) )

  return
end
subroutine y_to_s_english ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_ENGLISH writes an English year into a string.
!
!  Format:
!
!    YearNumber BC OS
!    YearNumber AD OS
!    YearNumber AD NS
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
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = * ) s
  integer ( kind = 4 ) y

  call y_check_english ( y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Y_TO_S_ENGLISH - Fatal error!'
    write ( *, '(a)' ) '  Illegal year.'
    stop
  end if

  if ( y < 0 ) then
    s(1:6) = 'BC OS '
    call i4_to_s_left ( -y, s(7:) )
  else if ( y <= 1752 ) then
    s(1:6) = 'AD OS'
    call i4_to_s_left ( y, s(7:) )
  else if ( 1752 < y ) then
    s(1:6) = 'AD NS'
    call i4_to_s_left ( y, s(7:) )
  end if

  return
end
subroutine y_to_s_ethiopian ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_ETHIOPIAN writes an Ethiopian year into a string.
!
!  Format:
!
!    Ethiopian YearNumber
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  character ( len = * ) s
  integer ( kind = 4 ) y

  s(1:10) = 'Ethiopian '
  call i4_to_s_left ( y, s(11:) )

  return
end
subroutine y_to_s_greek ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_GREEK writes a Greek year into a string.
!
!  Format:
!
!    OL 87.1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) o
  character ( len = * ) s
  character ( len = 5 ) so
  character ( len = 1 ) syy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) yy

  call y_check_greek ( y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Y_TO_S_GREEK - Fatal error!'
    write ( *, '(a)' ) '  Illegal year.'
    stop
  end if

  o = 1 + ( ( y - 1 ) / 4 )
  yy = i4_wrap ( y, 1, 4 )

  call i4_to_s_left ( o, so )
  call i4_to_s_left ( yy, syy )

  s = 'OL ' // trim ( so ) // '.' // trim ( syy )

  return
end
subroutine y_to_s_gregorian ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_GREGORIAN writes a Gregorian year into a string.
!
!  Format:
!
!    YearNumber BC
!    YearNumber AD
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = * ) s
  integer ( kind = 4 ) y

  call y_check_gregorian ( y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Y_TO_S_GREGORIAN - Fatal error!'
    write ( *, '(a)' ) '  Illegal year.'
    stop
  end if

  if ( y < 0 ) then
    s(1:3) = 'BC '
    call i4_to_s_left ( -y, s(4:) )
  else if ( 0 < y ) then
    s(1:3) = 'AD'
    call i4_to_s_left ( y, s(4:) )
  end if

  return
end
subroutine y_to_s_hebrew ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_HEBREW writes a Hebrew year into a string.
!
!  Format:
!
!    YearNumber AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = 10 ) sy
  character ( len = * ) s
  integer ( kind = 4 ) y

  call y_check_hebrew ( y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Y_TO_S_HEBREW - Fatal error!'
    write ( *, '(a)' ) '  Illegal year.'
    stop
  end if

  call i4_to_s_left ( y, sy )

  call s_cat1 ( sy, 'AM', s )

  return
end
subroutine y_to_s_islamic ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_ISLAMIC writes an Islamic year into a string.
!
!  Format:
!
!    YearNumber AH
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
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = 10 ) sy
  character ( len = * ) s
  integer ( kind = 4 ) y

  call y_check_islamic ( y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Y_TO_S_ISLAMIC - Fatal error!'
    write ( *, '(a)' ) '  Illegal year.'
    stop
  end if

  call i4_to_s_left ( y, sy )

  call s_cat1 ( sy, 'AH', s )

  return
end
subroutine y_to_s_julian ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_JULIAN writes a Julian year into a string.
!
!  Format:
!
!    YearNumber BC
!    YearNumber AD
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = * ) s
  integer ( kind = 4 ) y

  call y_check_julian ( y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Y_TO_S_JULIAN - Fatal error!'
    write ( *, '(a)' ) '  Illegal year.'
    stop
  end if

  if ( y < 0 ) then
    s(1:3) = 'BC '
    call i4_to_s_left ( -y, s(4:) )
  else if ( 0 < y ) then
    s(1:3) = 'AD'
    call i4_to_s_left ( y, s(4:) )
  end if

  return
end
subroutine y_to_s_persian ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_PERSIAN writes a Persian year into a string.
!
!  Format:
!
!    AP YearNumber
!
!  Discussion:
!
!    "AP" stands for "anno Persico".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  character ( len = * ) s
  integer ( kind = 4 ) y

  s(1:3) = 'AP '
  call i4_to_s_left ( y, s(4:) )

  return
end
subroutine y_to_s_republican ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_REPUBLICAN writes a Republican year into a string.
!
!  Format:
!
!    YearNumber ER
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = 10 ) sy
  character ( len = * ) s
  integer ( kind = 4 ) y

  call y_check_republican ( y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Y_TO_S_REPUBLICAN - Fatal error!'
    write ( *, '(a)' ) '  Illegal year.'
    stop
  end if

  call i4_to_s_left ( y, sy )

  call s_cat1 ( sy, 'ER', s )

  return
end
subroutine y_to_s_roman ( y, s )

!*****************************************************************************80
!
!! Y_TO_S_ROMAN writes a Roman year into a string.
!
!  Format:
!
!    YearNumber AUC
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, character ( len = * ) S, a representation of the year.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = 30 ) sy
  character ( len = * ) s
  integer ( kind = 4 ) y

  call y_check_roman ( y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Y_TO_S_ROMAN - Fatal error!'
    write ( *, '(a)' ) '  Illegal year.'
    stop
  end if

  call i4_to_roman ( y, sy )

  call s_cat1 ( sy, 'AUC', s )

  return
end
subroutine year_cal_common ( y )

!*****************************************************************************80
!
!! YEAR_CAL_COMMON prints out a calendar for a Common year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  character ( len = 20 ) lines1(6)
  character ( len = 20 ) lines2(6)
  character ( len = 20 ) lines3(6)
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  character ( len = 10 ) s1
  character ( len = 10 ) s2
  character ( len = 10 ) s3
  character ( len = 10 ) s4
  integer ( kind = 4 ) w
  character ( len = 2 ) weekdays(7)
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Make local copies of the input.
!
  y2 = y
!
!  Check the year.
!
  call y_check_common ( y2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make the year heading.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMMON CALENDAR'
  call y_to_s_common ( y2, s4 )
  write ( *, '(6x,a)' ) trim ( s4 )
!
!  Print out the month headings.
!
  do m1 = 1, 12, 3

    m2 = m1 + 1
    m3 = m2 + 1

    call month_to_month_name_common ( m1, s1 )
    call month_to_month_name_common ( m2, s2 )
    call month_to_month_name_common ( m3, s3 )

    write ( *, '(a)' ) ' '
    write ( *, '(5x,a,5x,2x,5x,a,5x,2x,5x,a)' ) s1, s2, s3
    write ( *, '(a)' ) ' '
!
!  Get the days of the week.
!
    do w = 1, 7
      call weekday_to_name_common2 ( w, weekdays(w) )
    end do

    write ( *, '(7(a2,1x),1x,7(a2,1x),1x,7(a2,1x))' ) &
      weekdays(1:7), weekdays(1:7), weekdays(1:7)
    write ( *, '(a)' ) ' '

    call month_cal_store_common ( y, m1, lines1 )
    call month_cal_store_common ( y, m2, lines2 )
    call month_cal_store_common ( y, m3, lines3 )

    do i = 1, 6
      write ( *, '(a,2x,a,2x,a)' ) lines1(i), lines2(i), lines3(i)
    end do

  end do

  return
end
function year_is_embolismic_eg_lunar ( y )

!*****************************************************************************80
!
!! YEAR_IS_EMBOLISMIC_EG_LUNAR: TRUE if the Egyptian Lunar year was embolismic.
!
!  Discussion:
!
!    This is just a "fake" function, which does repeat every 25 years,
!    and has 9 embolismic and 16 common years in that cycle, but with
!    a pattern I just made up for now.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_EMBOLISMIC_EG_LUNAR, TRUE if the year
!    was embolismic.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_embolismic_eg_lunar

  y2 = mod ( y - 1, 25 )

  if ( mod ( y2, 3 ) == 0 ) then
    year_is_embolismic_eg_lunar = .true.
  else
    year_is_embolismic_eg_lunar = .false.
  end if

  return
end
function year_is_embolismic_greek ( y )

!*****************************************************************************80
!
!! YEAR_IS_EMBOLISMIC_GREEK returns TRUE if the Greek year was embolismic.
!
!  Discussion:
!
!    Apparently, the Greek calendar was emended haphazardly.  This
!    routine does not attempt to follow that historical pattern, and
!    just uses the Hebrew calendar pattern for now.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_EMBOLISMIC_GREEK, TRUE if the year was embolismic.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) y
  logical year_is_embolismic_greek

  if ( 12 <= i4_modp ( 7 * y + 13, 19 ) ) then
    year_is_embolismic_greek = .true.
  else
    year_is_embolismic_greek = .false.
  end if

  return
end
function year_is_embolismic_hebrew ( y )

!*****************************************************************************80
!
!! YEAR_IS_EMBOLISMIC_HEBREW returns TRUE if the Hebrew year was embolismic.
!
!  Discussion:
!
!    In a 19 year cycle, there are 7 embolismic years.  During these years,
!    an extra month, "Adar II", (sometimes called "Veadar") is inserted after
!    the month of Adar.  Nonembolismic years are called "common" years.
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
!    Output, logical YEAR_IS_EMBOLISMIC_HEBREW, TRUE if the year was embolismic.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) y
  logical year_is_embolismic_hebrew

  if ( 12 <= i4_modp ( 7 * y + 13, 19 ) ) then
    year_is_embolismic_hebrew = .true.
  else
    year_is_embolismic_hebrew = .false.
  end if

  return
end
function year_is_leap_alexandrian ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_ALEXANDRIAN: TRUE if the Alexandrian year was a leap year.
!
!  Discussion:
!
!    The Alexandrian year, which started on the 29th of August of the Julian
!    year, was a leap year if it included the bissextile day of the Julian
!    calendar.  In other words, if the Alexandrian year BEGAN in a Julian year
!    that preceded a Julian leap year, then the Alexandrian year was a leap
!    year.
!
!    We deem year AX 1 to have begun in Julian 23 BC.  Julian 21 BC was
!    theoretically a leap year, so AX 2 was a leap year, as was AX 6, AX 10,
!    and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_ALEXANDRIAN, TRUE if the year was a leap year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_alexandrian

  if ( mod ( y, 4 ) == 2 ) then
    year_is_leap_alexandrian = .true.
  else
    year_is_leap_alexandrian = .false.
  end if

  return
end
function year_is_leap_bahai ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_BAHAI returns TRUE if the Bahai year was a leap year.
!
!  Discussion:
!
!    The leap year rules are the same as those used in the Gregorian
!    calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_BAHAI, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_bahai

  if ( y <= 0 ) then
    year_is_leap_bahai = .false.
    return
  end if

  if ( mod ( y, 400 ) == 0 ) then
    year_is_leap_bahai = .true.
  else if ( mod ( y, 100 ) == 0 ) then
    year_is_leap_bahai = .false.
  else if ( mod ( y, 4 ) == 0 ) then
    year_is_leap_bahai = .true.
  else
    year_is_leap_bahai = .false.
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
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YEAR_IS_LEAP_COMMON - Fatal error!'
    write ( *, '(a)' ) '  Year 0 is illegal.'
    stop
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
function year_is_leap_coptic ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_COPTIC returns TRUE if the Coptic year was a leap year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Nachum Dershowitz, Edward Reingold,
!    Calendrical Calculations,
!    Cambridge, 1997, page 58.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_COPTIC, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_coptic

  if ( y <= 0 ) then
    year_is_leap_coptic = .false.
    return
  end if

  if ( mod ( y, 4 ) == 3 ) then
    year_is_leap_coptic = .true.
  else
    year_is_leap_coptic = .false.
  end if

  return
end
function year_is_leap_eg_lunar ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_EG_LUNAR: TRUE if the Egyptian Lunar year was a leap year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_EG_LUNAR, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_eg_lunar

  if ( y <= 0 ) then
    year_is_leap_eg_lunar = .false.
    return
  end if

  if ( mod ( y, 5 ) == 0 ) then
    year_is_leap_eg_lunar = .true.
  else
    year_is_leap_eg_lunar = .false.
  end if

  return
end
function year_is_leap_english ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_ENGLISH returns TRUE if the English year was a leap year.
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
!    else if ( the year is no greater than 1752 ) then
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
!    Output, logical YEAR_IS_LEAP_ENGLISH, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_english

  if ( y == 0 ) then
    year_is_leap_english = .false.
    return
  end if
!
!  BC years have to have 1 added to them to make a proper leap year evaluation.
!
  call y_common_to_astronomical ( y, y2 )

  if ( y2 <= 1752 ) then

    if ( i4_modp ( y2, 4 ) == 0 ) then
      year_is_leap_english = .true.
    else
      year_is_leap_english = .false.
    end if

  else

    if ( i4_modp ( y2, 400 ) == 0 ) then
      year_is_leap_english = .true.
    else if ( i4_modp ( y2, 100 ) == 0 ) then
      year_is_leap_english = .false.
    else if ( i4_modp ( y2, 4 ) == 0 ) then
      year_is_leap_english = .true.
    else
      year_is_leap_english = .false.
    end if

  end if

  return
end
function year_is_leap_ethiopian ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_ETHIOPIAN returns TRUE if the Ethiopian year was a leap year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Nachum Dershowitz, Edward Reingold,
!    Calendrical Calculations,
!    Cambridge, 1997, page 58.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_ETHIOPIAN, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_ethiopian

  if ( y <= 0 ) then
    year_is_leap_ethiopian = .false.
    return
  end if

  if ( mod ( y, 4 ) == 3 ) then
    year_is_leap_ethiopian = .true.
  else
    year_is_leap_ethiopian = .false.
  end if

  return
end
function year_is_leap_greek ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_GREEK returns TRUE if the Greek year was a leap year.
!
!  Discussion:
!
!    The actual practice of adding the extra day to the Greek calendar
!    seems to have been unmethodical.  Here, we simply make up a rule
!    as a placeholder for now.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_GREEK, TRUE if the year was a leap year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_embolismic_greek
  logical year_is_leap_greek

  if ( year_is_embolismic_greek ( y ) .and. ( mod ( y, 3 ) == 0 ) ) then
    year_is_leap_greek = .true.
  else
    year_is_leap_greek = .false.
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
function year_is_leap_iranian ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_IRANIAN returns TRUE if the Iranian year was a leap year.
!
!  Discussion:
!
!    I don't know the rule for this, so I'm just setting it FALSE for now.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_IRANIAN, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_iranian

  year_is_leap_iranian = .false.

  return
end
function year_is_leap_islamic ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_ISLAMIC returns TRUE if the Islamic year was a leap year.
!
!  Discussion:
!
!    In a 30 year cycle, there are 11 leap years, years 2, 5, 7, 10, 13,
!    16, 18, 21, 24, 26 and 29.  During these years, the 12th month has
!    30 days instead of 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_ISLAMIC, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) y
  logical year_is_leap_islamic

  if ( i4_modp ( 11 * y + 14, 30 ) < 11 ) then
    year_is_leap_islamic = .true.
  else
    year_is_leap_islamic = .false.
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
function year_is_leap_persian ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_PERSIAN returns TRUE if the Persian year was a leap year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Nachum Dershowitz, Edward Reingold,
!    Calendrical Calculations,
!    Cambridge, 1997, page 58.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_PERSIAN, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3
  logical year_is_leap_persian

  if ( y <= 0 ) then
    y2 = y - 473
  else
    y2 = y - 474
  end if

  y3 = 474 + mod ( y2, 2820 )

  if ( mod ( 682 * ( y3 + 38 ), 2816 ) < 682 ) then
    year_is_leap_persian = .true.
  else
    year_is_leap_persian = .false.
  end if

  return
end
function year_is_leap_republican ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_REPUBLICAN returns TRUE if the Republican year was a leap year.
!
!  Discussion:
!
!    The French Republican calendar was in use for 14 years.
!    In that time, years 3, 7 and 11 were designated as leap years.
!    The easiest way to harmonize the rules and history is to apply
!    the leap year rules to Y+1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_REPUBLICAN, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_republican

  y2 = y

  call y_check_republican ( y2, ierror )

  if ( ierror /= 0 ) then
    year_is_leap_republican = .false.
    return
  end if

  year_is_leap_republican = .false.

  if ( mod ( y2+1, 4 ) == 0 ) then
    year_is_leap_republican = .true.
    if ( mod ( y2+1, 100 ) == 0 ) then
      year_is_leap_republican = .false.
      if ( mod ( y2+1, 400 ) == 0 ) then
        year_is_leap_republican = .true.
        if ( mod ( y2+1, 4000 ) == 0 ) then
          year_is_leap_republican = .false.
        end if
      end if
    end if
  end if

  return
end
function year_is_leap_roman ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_ROMAN returns TRUE if the Roman year was a leap year.
!
!  Discussion:
!
!    For our unrealistic and idealized Roman calendar, we are going to
!    take a year to have been a leap year if the corresponding year in
!    the idealized Julian calendar was a leap year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_ROMAN, TRUE if the year was a leap year.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_julian
  logical year_is_leap_roman

  call y_check_roman ( y, ierror )

  if ( ierror /= 0 ) then
    year_is_leap_roman = .false.
    return
  end if

  call y_roman_to_julian ( y, y2 )

  year_is_leap_roman = year_is_leap_julian ( y2 )

  return
end
function year_length_alexandrian ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_ALEXANDRIAN returns the number of days in an Alexandrian year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_ALEXANDRIAN, the number of
!    days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_alexandrian
  integer ( kind = 4 ) year_length_alexandrian

  if ( year_is_leap_alexandrian ( y ) ) then
    year_length_alexandrian = 366
  else
    year_length_alexandrian = 365
  end if

  return
end
function year_length_bahai ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_BAHAI returns the number of days in a Bahai year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_BAHAI, the number of
!    days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_bahai
  integer ( kind = 4 ) year_length_bahai

  if ( year_is_leap_bahai ( y ) ) then
    year_length_bahai = 366
  else
    year_length_bahai = 365
  end if

  return
end
function year_length_common ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_COMMON returns the number of days in a Common year.
!
!  Discussion:
!
!    The "common" calendar is meant to be the calendar which is Julian up to
!    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
!
!    If Y is 0, then the routine returns 0, reflecting the fact that
!    there was officially no year 0.
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
!    Output, integer ( kind = 4 ) YEAR_LENGTH_COMMON, the number of
!    days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_common
  integer ( kind = 4 ) year_length_common

  if ( y == 0 ) then
    year_length_common = 0
  else if ( y == 1582 ) then
    year_length_common = 355
  else if ( year_is_leap_common ( y ) ) then
    year_length_common = 366
  else
    year_length_common = 365
  end if

  return
end
function year_length_coptic ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_COPTIC returns the number of days in a Coptic year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_COPTIC, the number of
!    days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_coptic
  integer ( kind = 4 ) year_length_coptic

  if ( year_is_leap_coptic ( y ) ) then
    year_length_coptic = 366
  else
    year_length_coptic = 365
  end if

  return
end
function year_length_eg_civil ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_EG_CIVIL returns the number of days in an Egyptian Civil year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_EG_CIVIL, the number of
!    days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_eg_civil

  year_length_eg_civil = 365

  return
end
function year_length_eg_lunar ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_EG_LUNAR returns the number of days in an Egyptian Lunar year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_EG_LUNAR, the number of
!    days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_embolismic_eg_lunar
  logical year_is_leap_eg_lunar
  integer ( kind = 4 ) year_length_eg_lunar

  if ( .not. year_is_embolismic_eg_lunar ( y ) ) then
    year_length_eg_lunar = 354
  else
    year_length_eg_lunar = 384
  end if

  if ( year_is_leap_eg_lunar ( y ) ) then
    year_length_eg_lunar = year_length_eg_lunar + 1
  end if

  return
end
function year_length_english ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_ENGLISH returns the number of days in an English year.
!
!  Discussion:
!
!    The "English" calendar is meant to be the calendar which is Julian before
!    the transition date, and Gregorian afterwards.
!
!    1752 was a special year with only 355 days instead of 366.
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
!    Output, integer ( kind = 4 ) YEAR_LENGTH_ENGLISH, the number of
!    days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_english
  integer ( kind = 4 ) year_length_english

  if ( y == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YEAR_LENGTH_ENGLISH - Fatal error!'
    write ( *, '(a)' ) '  Illegal Y = 0.'
    stop
  end if

  if ( y == 1752 ) then
    year_length_english = 355
  else if ( year_is_leap_english ( y ) ) then
    year_length_english = 366
  else
    year_length_english = 365
  end if

  return
end
function year_length_ethiopian ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_ETHIOPIAN returns the number of days in an Ethiopian year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_ETHIOPIAN, the number of
!    days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_ethiopian
  integer ( kind = 4 ) year_length_ethiopian

  if ( year_is_leap_ethiopian ( y ) ) then
    year_length_ethiopian = 366
  else
    year_length_ethiopian = 365
  end if

  return
end
function year_length_greek ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_GREEK returns the number of days in a Greek year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_GREEK, the number of
!    days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_embolismic_greek
  logical year_is_leap_greek
  integer ( kind = 4 ) year_length_greek

  if ( year_is_embolismic_greek ( y ) ) then
    year_length_greek = 386
    if ( year_is_leap_greek ( y ) ) then
      year_length_greek = year_length_greek + 1
    end if
  else
    year_length_greek = 357
  end if

  return
end
function year_length_gregorian ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_GREGORIAN returns the number of days in a Gregorian year.
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
!    Output, integer ( kind = 4 ) YEAR_LENGTH_GREGORIAN, the number of
!    days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_gregorian
  integer ( kind = 4 ) year_length_gregorian

  if ( y == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'YEAR_LENGTH_GREGORIAN - Fatal error!'
    write ( *, * ) '  Illegal Y = 0.'
    stop
  end if

  if ( year_is_leap_gregorian ( y ) ) then
    year_length_gregorian = 366
  else
    year_length_gregorian = 365
  end if

  return
end
function year_length_hebrew ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_HEBREW returns the number of days in a Hebrew year.
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
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_HEBREW, the number of
!    days in the year.
!
  implicit none

  real ( kind = 8 ) jed
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) year_length_hebrew

  call new_year_to_jed_hebrew ( y, jed )

  y2 = y + 1
  call new_year_to_jed_hebrew ( y2, jed2 )

  year_length_hebrew = nint ( jed2 - jed )

  return
end
function year_length_hindu_solar ( )

!*****************************************************************************80
!
!! YEAR_LENGTH_HINDU_SOLAR returns the number of days in a Hindu solar year.
!
!  Discussion:
!
!    Warning: This is a DOUBLE PRECISION quantity.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) YEAR_LENGTH_HINDU_SOLAR, the number of days
!    in the year.
!
  implicit none

  real ( kind = 8 ) year_length_hindu_solar

  year_length_hindu_solar = real ( 1577917828, kind = 8 ) &
    / real ( 4320000, kind = 8 )

  return
end
function year_length_islamic ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_ISLAMIC returns the number of days in an Islamic year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_ISLAMIC, the number of
!    days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_islamic
  integer ( kind = 4 ) year_length_islamic

  if ( year_is_leap_islamic ( y ) ) then
    year_length_islamic = 355
  else
    year_length_islamic = 354
  end if

  return
end
function year_length_julian ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_JULIAN returns the number of days in a Julian year.
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
!    Output, integer ( kind = 4 ) YEAR_LENGTH_JULIAN, the number of
!    days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_julian
  integer ( kind = 4 ) year_length_julian

  if ( y == 0 ) then
    year_length_julian = 0
  else if ( year_is_leap_julian ( y ) ) then
    year_length_julian = 366
  else
    year_length_julian = 365
  end if

  return
end
function year_length_lunar ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_LUNAR returns the number of days in a "lunar year".
!
!  Discussion:
!
!    The "lunar year" is taken to be the length of 12 mean lunations.
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
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, real ( kind = 8 ) YEAR_LENGTH_LUNAR, the number of days
!    in the year.
!
  implicit none

  integer ( kind = 4 ) y
  real ( kind = 8 ) year_length_lunar

  year_length_lunar = 354.3671E+00

  return
end
function year_length_months_alexandrian ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_ALEXANDRIAN: number of months in an Alexandrian year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_ALEXANDRIAN, the
!    number of months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_alexandrian

  year_length_months_alexandrian = 13

  return
end
function year_length_months_bahai ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_BAHAI returns the number of months in a Bahai year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_BAHAI, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_bahai

  year_length_months_bahai = 20

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
function year_length_months_coptic ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_COPTIC returns the number of months in a Coptic year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_COPTIC, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_coptic

  year_length_months_coptic = 13

  return
end
function year_length_months_eg_civil ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_EG_CIVIL: number of months in an Egyptian Civil year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_EG_CIVIL, the number
!    of months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_eg_civil

  year_length_months_eg_civil = 13

  return
end
function year_length_months_eg_lunar ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_EG_LUNAR: number of months in an Egyptian Lunar year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_EG_LUNAR, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_embolismic_eg_lunar
  integer ( kind = 4 ) year_length_months_eg_lunar

  if ( year_is_embolismic_eg_lunar ( y ) ) then
    year_length_months_eg_lunar = 13
  else
    year_length_months_eg_lunar = 12
  end if

  return
end
function year_length_months_english ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_ENGLISH returns the number of months in an English year.
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
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_ENGLISH, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_english

  year_length_months_english = 12

  return
end
function year_length_months_ethiopian ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_ETHIOPIAN: number of months in an Ethiopian year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_ETHIOPIAN, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_ethiopian

  year_length_months_ethiopian = 13

  return
end
function year_length_months_greek ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_GREEK returns the number of months in a Greek year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_GREEK, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_embolismic_greek
  integer ( kind = 4 ) year_length_months_greek

  if ( year_is_embolismic_greek ( y ) ) then
    year_length_months_greek = 13
  else
    year_length_months_greek = 12
  end if

  return
end
function year_length_months_gregorian ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_GREGORIAN: number of months in a Gregorian year.
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
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_GREGORIAN, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_gregorian

  year_length_months_gregorian = 12

  return
end
function year_length_months_hebrew ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_HEBREW returns the number of months in a Hebrew year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_HEBREW, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_embolismic_hebrew
  integer ( kind = 4 ) year_length_months_hebrew

  if ( year_is_embolismic_hebrew ( y ) ) then
    year_length_months_hebrew = 13
  else
    year_length_months_hebrew = 12
  end if

  return
end
function year_length_months_hindu_lunar ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_HINDU_LUNAR: number of months in a Hindu lunar year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_HINDU_LUNAR, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_hindu_lunar

  year_length_months_hindu_lunar = 12

  return
end
function year_length_months_hindu_solar ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_HINDU_SOLAR: number of months in a Hindu solar year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_HINDU_SOLAR, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_hindu_solar

  year_length_months_hindu_solar = 12

  return
end
function year_length_months_islamic ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_ISLAMIC returns the number of months in an Islamic year.
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
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_ISLAMIC, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_islamic

  year_length_months_islamic = 12

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
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_JULIAN, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_julian

  year_length_months_julian = 12

  return
end
function year_length_months_persian ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_PERSIAN returns the number of months in a Persian year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_PERSIAN, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_persian

  year_length_months_persian = 12

  return
end
function year_length_months_republican ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_REPUBLICAN: number of months in a Republican year.
!
!  Discussion:
!
!    The routine always returns the value 13.  The 13-th month was not
!    regarded as a month, but as a special 5 or 6 day period known as
!    the "Sansculottides".  For our purposes, it's the 13th month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_REPUBLICAN, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_republican

  year_length_months_republican = 13

  return
end
function year_length_months_roman ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_MONTHS_ROMAN returns the number of months in a Roman year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_MONTHS_ROMAN, the number of
!    months in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_roman

  year_length_months_roman = 12

  return
end
function year_length_persian ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_PERSIAN returns the number of days in a Persian year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_PERSIAN, the number of days
!    in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_persian
  integer ( kind = 4 ) year_length_persian

  if ( year_is_leap_persian ( y ) ) then
    year_length_persian = 366
  else
    year_length_persian = 365
  end if

  return
end
function year_length_republican ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_REPUBLICAN returns the number of days in a Republican year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_ENGLISH, the number of days
!    in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_republican
  integer ( kind = 4 ) year_length_republican

  if ( year_is_leap_republican ( y ) ) then
    year_length_republican = 366
  else
    year_length_republican = 365
  end if

  return
end
function year_length_roman ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_ROMAN returns the number of days in a Roman year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_ROMAN, the number of days
!    in the year.
!
  implicit none

  integer ( kind = 4 ) y
  logical year_is_leap_roman
  integer ( kind = 4 ) year_length_roman

  if ( year_is_leap_roman ( y ) ) then
    year_length_roman = 366
  else
    year_length_roman = 365
  end if

  return
end
function year_length_solar ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_SOLAR returns the number of days in a "solar" year.
!
!  Discussion:
!
!    The "solar" year is taken to be the mean tropical year.
!    The number of days in a mean tropical year has varied from
!    365.2424992 in 4000 BC to 365.2421897 in 2000 AD.
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
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, real ( kind = 8 ) YEAR_LENGTH_SOLAR, the number of days
!    in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  real ( kind = 8 ) year_length_solar

  if ( y < 1 ) then
    y2 = y + 1
  else
    y2 = y
  end if

  if ( y2 < - 4000 ) then
    year_length_solar = 365.2424992D+00
  else if ( y2 <= 2000 ) then
    year_length_solar = &
      ( real ( 2000 - y2, kind = 8 ) * 365.2424992D+00 &
      + real ( 4000 + y2, kind = 8 ) * 365.2421897D+00 ) &
      /        6000.0D+00
  else
    year_length_solar = 365.2421897D+00
  end if

  return
end
subroutine year_to_dominical_common ( y, n1, n2 )

!*****************************************************************************80
!
!! YEAR_TO_DOMINICAL_COMMON: dominical numbers, Common calendar.
!
!  Discussion:
!
!    The Julian calendar calculations are used through the year 1582,
!    and the Gregorian thereafter.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, integer ( kind = 4 ) N1, N2, the dominical numbers for the year.
!    If Y is a leap year, then N1 applies before March 1, and N2 after.
!    If Y is not a leap year, then N1 applies throughout the year,
!    and N2 is returned as N1.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) y

  if ( y <= 1582 ) then
    call year_to_dominical_julian ( y, n1, n2 )
  else
    call year_to_dominical_gregorian ( y, n1, n2 )
  end if

  return
end
subroutine year_to_dominical_gregorian ( y, n1, n2 )

!*****************************************************************************80
!
!! YEAR_TO_DOMINICAL_GREGORIAN: dominical numbers, Gregorian calendar.
!
!  Discussion:
!
!    The days of each year are numbered with "calendar letters", with
!    January 1 having letter 'A', January 7 having letter 'G', and
!    the cycle then repeating with January 8 having letter 'A'.
!
!    This cycle is independent of the weekday cycle.  If a year is
!    not a leap year, then all Sundays have the same calendar letter.
!    This is called the dominical letter of the year.  If a year is
!    a leap year, then all Sundays before March 1 have one calendar
!    letter, and all Sundays after have another (namely, the calendar
!    letter one position earlier in the cycle).
!
!    Using the correspondence A = 1, B = 2, ..., we may speak of
!    the dominical number of a year, or dominical numbers for a leap year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, integer ( kind = 4 ) N1, N2, the dominical numbers for the year.
!    If Y is a leap year, then N1 applies before March 1, and N2 after.
!    If Y is not a leap year, then N1 applies throughout the year,
!    and N2 is returned as N1.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) p1
  integer ( kind = 4 ) y
  logical year_is_leap_gregorian
  integer ( kind = 4 ) y2

  if ( y == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YEAR_TO_DOMINICAL_GREGORIAN - Fatal error!'
    write ( *, '(a)' ) '  Illegal input Y = 0.'
    stop
  end if

  call y_common_to_astronomical ( y, y2 )

  p1 = y2 + ( y2 / 4 ) - ( y2 / 100 ) + ( y2 / 400 ) - 1
  n1 = 7 - i4_modp ( p1, 7 )

  if ( year_is_leap_gregorian ( y2 ) ) then
    n2 = n1
    p1 = p1 - 1
    n1 = 7 - i4_modp ( p1, 7 )
  else
    n2 = n1
  end if

  return
end
subroutine year_to_dominical_julian ( y, n1, n2 )

!*****************************************************************************80
!
!! YEAR_TO_DOMINICAL_JULIAN: dominical numbers, Julian calendar.
!
!  Discussion:
!
!    The days of each year are numbered with "calendar letters", with
!    January 1 having letter 'A', January 7 having letter 'G', and
!    the cycle then repeating with January 8 having letter 'A'.
!
!    This cycle is independent of the weekday cycle.  If a year is
!    not a leap year, then all Sundays have the same calendar letter.
!    This is called the dominical letter of the year.  If a year is
!    a leap year, then all Sundays before March 1 have one calendar
!    letter, and all Sundays after have another (namely, the calendar
!    letter one position earlier in the cycle).
!
!    Using the correspondence A = 1, B = 2, ..., we may speak of
!    the dominical number of a year, or dominical numbers for a leap year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.  The year 0 is illegal input.
!
!    Output, integer ( kind = 4 ) N1, N2, the dominical numbers for the year.
!    If Y is a leap year, then N1 applies before March 1, and N2 after.
!    If Y is not a leap year, then N1 applies throughout the year,
!    and N2 is returned as N1.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) p1
  integer ( kind = 4 ) y
  logical year_is_leap_julian
  integer ( kind = 4 ) y2

  if ( y == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YEAR_TO_DOMINICAL_JULIAN - Fatal error!'
    write ( *, '(a)' ) '  Illegal input Y = 0.'
    stop
  end if

  call y_common_to_astronomical ( y, y2 )

  p1 = y2 + ( y2 / 4 ) + 4
  n1 = 7 - i4_modp ( p1, 7 )

  if ( year_is_leap_julian ( y2 ) ) then
    n2 = n1
    p1 = p1 - 1
    n1 = 7 - i4_modp ( p1, 7 )
  else
    n2 = n1
  end if

  return
end
subroutine year_to_epact_gregorian ( y, e )

!*****************************************************************************80
!
!! YEAR_TO_EPACT_GREGORIAN returns the epact of a Gregorian year.
!
!  Discussion:
!
!    The epact of a year is the age in days of the notional moon on
!    the first day of the year.  If the year begins with a new moon,
!    the epact is zero.  If the new moon occurred the day before,
!    the epact is 1.  There is a unique epact for every golden number.
!
!    The Gregorian epact calculation is an adjustment to the Julian
!    calculation that takes into account the shift of the calendar
!    to restore the vernal equinox to March 21, and the adjustment to
!    the average length of the year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2000
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
!    Input, integer ( kind = 4 ) Y, the year.  The year 0 is illegal input.
!
!    Output, integer ( kind = 4 ) E, the epact, between 0 and 28.
!
  implicit none

  integer ( kind = 4 ) e
  integer ( kind = 4 ) g
  integer ( kind = 4 ) h
  integer ( kind = 4 ) q
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2

  if ( y == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YEAR_TO_EPACT_GREGORIAN - Fatal error!'
    write ( *, '(a)' ) '  Illegal input Y = 0.'
    stop
  end if

  call y_common_to_astronomical ( y, y2 )

  call year_to_golden_number ( y, g )

  h = ( y2 / 100 )

  q = h - ( h / 4 )

  e = mod ( 57 + 11 * g - q + ( h - ( h - 17 ) / 25 ) / 3, 30 )

  if ( e == 24 .or. ( e == 25 .and. 12 <= g ) ) then
    e = e + 1
  end if

  return
end
subroutine year_to_epact_julian ( y, e )

!*****************************************************************************80
!
!! YEAR_TO_EPACT_JULIAN returns the epact of a Julian year.
!
!  Discussion:
!
!    The epact of a year is the age in days of the notional moon on
!    the first day of the year.  If the year begins with a new moon,
!    the epact is zero.  If the new moon occurred the day before,
!    the epact is 1.  There is a unique epact for every golden number.
!
!    Bear in mind that the notional moon is not the one in the sky,
!    but a theoretical one that satisfactorily approximates the behavior
!    of the real one, but which is tame enough to be described by a formula.
!
!  Example:
!
!    Year  Golden Number  Epact
!
!      1 BC     1           8
!      1 AD     2          19
!      2 AD     3           0
!      3 AD     4          11
!      4 AD     5          22
!      5 AD     6           3
!      6 AD     7          14
!      7 AD     8          25
!      8 AD     9           6
!      9 AD    10          17
!     10 AD    11          28
!     11 AD    12           9
!     12 AD    13          20
!     13 AD    14           1
!     14 AD    15          12
!     15 AD    16          23
!     16 AD    17           4
!     17 AD    18          15
!     18 AD    19          26
!     19 AD     1           8
!     20 AD     2          19
!   1066 AD     3           0
!   1900 AD     1           8
!   1919 AD     1           8
!   1938 AD     1           8
!   1957 AD     1           8
!   1976 AD     1           8
!   1995 AD     1           8
!   2014 AD     1           8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2000
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
!    Input, integer ( kind = 4 ) Y, the year.  The year 0 is illegal input.
!
!    Output, integer ( kind = 4 ) E, the epact, between 0 and 28.
!
  implicit none

  integer ( kind = 4 ) e
  integer ( kind = 4 ) g
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) y

  if ( y == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YEAR_TO_EPACT_JULIAN - Fatal error!'
    write ( *, '(a)' ) '  Illegal input Y = 0.'
    stop
  end if

  call year_to_golden_number ( y, g )

  e = i4_wrap ( 11 * g - 3, 0, 29 )

  return
end
subroutine year_to_golden_number ( y, g )

!*****************************************************************************80
!
!! YEAR_TO_GOLDEN_NUMBER returns the golden number of a Common year.
!
!  Discussion:
!
!    Nineteen solar years are very close to 235 lunations.  Calendars
!    that try to keep track of both the sun and moon often make use of
!    this fact, ascribed to the Greek astronomer Meton.
!
!    While trying to determine a formula for Easter, Dionysus Exiguus
!    symbolized the place of each year in its Metonic cycle by a
!    "golden number" between 1 and 19.  The numbering began with the
!    year 1 BC, assigned the golden number of 1.  The following year,
!    1 AD, got the golden number of 2, and after that it gets easier.
!
!    The same golden year calculation is done for years in the Julian
!    or Gregorian calendar.
!
!  Example:
!
!    Year  Golden Number
!
!      1 BC     1
!      1 AD     2
!      2 AD     3
!     18 AD    19
!     19 AD     1
!     20 AD     2
!   1066 AD     3
!   1900 AD     1
!   1919 AD     1
!   1938 AD     1
!   1957 AD     1
!   1976 AD     1
!   1995 AD     1
!   2014 AD     1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, integer ( kind = 4 ) G, the golden number, between 1 and 19.  This
!    records the position of the year in the 19 year Metonic cycle.
!
  implicit none

  integer ( kind = 4 ) g
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2

  if ( y == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YEAR_TO_GOLDEN_NUMBER - Fatal error!'
    write ( *, '(a)' ) '  Illegal input Y = 0.'
    stop
  end if
!
!  We assume that BC years come in as negative numbers, and that
!  the year before 1 AD is 1 BC.  So add 1 to any negative value
!  so that the arithmetic works.
!
  call y_common_to_astronomical ( y, y2 )

  g = i4_wrap ( y2 + 1, 1, 19 )

  return
end
subroutine year_to_indiction_common ( y, i )

!*****************************************************************************80
!
!! YEAR_TO_INDICTION_COMMON returns the indiction number of a Common year.
!
!  Discussion:
!
!    The Roman empire had a taxation cycle that, at one time, comprised
!    15 years.  As is typical in calendrical matters, the actual length
!    of this cycle and the time that the cycle began varied from place
!    to place and time to time, and historians even disagree about the
!    indiction cycle given a specific place and time.  Nonetheless,
!    it is customary to retrospectively impose a uniform and regular
!    indiction cycle on the ancient world.  (The 15 year indiction cycle,
!    in fact, was factored into Scaliger's determination of an appropriate
!    starting point for the Julian Ephemeris Date.)
!
!  Example:
!
!    Year  Indiction Number
!
!      3 BC     1
!      2 BC     2
!      1 BC     3
!      1 AD     4
!     10 AD    13
!     11 AD    14
!     12 AD    15
!     13 AD     1
!     14 AD     2
!     15 AD     3
!     26 AD    14
!     27 AD    15
!     28 AD     1
!   1900 AD    13
!   2000 AD     8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year.
!
!    Output, integer ( kind = 4 ) I, the indiction number, between 1 and 15.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2

  if ( y == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YEAR_TO_INDICTION_COMMON - Fatal error!'
    write ( *, '(a)' ) '  Illegal input Y = 0.'
    stop
  end if
!
!  We assume that BC years come in as negative numbers, and that
!  the year before 1 AD is 1 BC.  So add 1 to any negative value
!  so that the arithmetic works.
!
  call y_common_to_astronomical ( y, y2 )

  i = i4_wrap ( y2 + 3, 1, 15 )

  return
end
subroutine year_to_scaliger_common ( y, c1, c2, c3, r1, r2, r3 )

!*****************************************************************************80
!
!! YEAR_TO_SCALIGER_COMMON converts a Common year to its Scaliger indices.
!
!  Discussion:
!
!    The year 4713 BCE was chosen by Joseph Scaliger for the start of
!    his Julian Ephemeris Date system, because three cycles coincided
!    in that year, the 28 year Julian calendar cycle, the 19 year Metonic
!    cycle, and the 15 year Roman Indiction cycle.  Thus, the year
!    4713 BCE has Scaliger index (1,1,1).  Each subsequent year has a distinct
!    set of Scaliger indices until 7980 years later, when the year
!    3266 CE will again have the Scaliger index (1,1,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the Common year.
!
!    Output, integer ( kind = 4 ) C1, C2, C3, the number of completed
!    Julian, Metonic and Indiction cycles.
!
!    Output, integer ( kind = 4 ) R1, R2, R3, the Julian, Metonic and
!    Indiction cycle numbers that make up the Scaliger index.
!
  implicit none

  integer ( kind = 4 ) c1
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) c3
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) r1
  integer ( kind = 4 ) r2
  integer ( kind = 4 ) r3
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2

  if ( y == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YEAR_TO_SCALIGER_COMMON - Fatal error!'
    write ( *, '(a)' ) '  Illegal input Y = 0.'
    stop
  end if
!
!  Adjust for missing year 0.
!
  if ( y < 0 ) then
    y2 = y + 1
  else
    y2 = y
  end if
!
!  Now shift so 4713 BC becomes the year 1.
!
  y2 = y2 + 4713

  c1 = ( y2 - 1 ) / 28
  c2 = ( y2 - 1 ) / 19
  c3 = ( y2 - 1 ) / 15

  r1 = i4_wrap ( y2, 1, 28 )
  r2 = i4_wrap ( y2, 1, 19 )
  r3 = i4_wrap ( y2, 1, 15 )

  return
end
subroutine year_to_type_hebrew ( y, type )

!*****************************************************************************80
!
!! YEAR_TO_TYPE_HEBREW returns the type of a Hebrew year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Richards,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, pages 332.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the Hebrew year.
!    Nonpositive years are illegal input.
!
!    Output, integer ( kind = 4 ) TYPE, the year type.
!    1, Common, Deficient, 12 months, 353 days;
!    2, Common, Regular, 12 months, 354 days;
!    3, Common, Abundant, 12 months, 355 days;
!    4, Embolismic, Deficient, 13 months, 383 days;
!    5, Embolismic, Regular, 13 months, 384 days;
!    6, Embolismic, Abundant, 13 months, 385 days.
!
  implicit none

  real ( kind = 8 ) jed
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) type
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length

  if ( y <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YEAR_TO_TYPE_HEBREW - Fatal error!'
    write ( *, '(a)' ) '  Illegal input Y = 0.'
    stop
  end if

  call new_year_to_jed_hebrew ( y, jed )

  call new_year_to_jed_hebrew ( y+1, jed2 )

  year_length = nint ( jed2 - jed )

       if ( year_length == 353 ) then
    type = 1
  else if ( year_length == 354 ) then
    type = 2
  else if ( year_length == 355 ) then
    type = 3
  else if ( year_length == 383 ) then
    type = 4
  else if ( year_length == 384 ) then
    type = 5
  else if ( year_length == 385 ) then
    type = 6
  else
    type = 0
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YEAR_TO_TYPE_HEBREW - Fatal error!'
    write ( *, '(a,i6)' ) '  Computed an illegal type = ', type
    stop
  end if

  return
end
subroutine yj_check_common ( y, j, ierror )

!*****************************************************************************80
!
!! YJ_CHECK_COMMON checks a Common YJ date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no serious error was found in
!    the date, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_common ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure J is not too small or too big.
!
  call j_borrow_common ( y, j )

  call j_carry_common ( y, j )

  return
end
subroutine yj_check_english ( y, j, ierror )

!*****************************************************************************80
!
!! YJ_CHECK_ENGLISH checks an English YJ date.
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
!    Input/output, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no serious error was found in
!    the date, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_english ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure J is not too small or too big.
!
  call j_borrow_english ( y, j )

  call j_carry_english ( y, j )

  return
end
subroutine yj_check_gregorian ( y, j, ierror )

!*****************************************************************************80
!
!! YJ_CHECK_GREGORIAN checks a Gregorian YJ date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no serious error was found in
!    the date, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_gregorian ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure J is not too small or too big.
!
  call j_borrow_gregorian ( y, j )

  call j_carry_gregorian ( y, j )

  return
end
subroutine yj_check_hebrew ( y, j, ierror )

!*****************************************************************************80
!
!! YJ_CHECK_HEBREW checks a Hebrew YJ date.
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
!    Input/output, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no serious error was found in
!    the date, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_hebrew ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure J is not too small or too big.
!
  call j_borrow_hebrew ( y, j )

  call j_carry_hebrew ( y, j )

  return
end
subroutine yj_check_islamic ( y, j, ierror )

!*****************************************************************************80
!
!! YJ_CHECK_ISLAMIC checks an Islamic YJ date.
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
!    Input/output, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no serious error was found in
!    the date, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_islamic ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure J is not too small or too big.
!
  call j_borrow_islamic ( y, j )

  call j_carry_islamic ( y, j )

  return
end
subroutine yj_check_julian ( y, j, ierror )

!*****************************************************************************80
!
!! YJ_CHECK_JULIAN checks a Julian YJ date.
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
!    Input/output, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no serious error was found in
!    the date, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_julian ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure J is not too small or too big.
!
  call j_borrow_julian ( y, j )

  call j_carry_julian ( y, j )

  return
end
subroutine yj_check_republican ( y, j, ierror )

!*****************************************************************************80
!
!! YJ_CHECK_REPUBLICAN checks a Republican YJ date.
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
!    Input/output, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no serious error was found in
!    the date, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_republican ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure J is not too small or too big.
!
  call j_borrow_republican ( y, j )

  call j_carry_republican ( y, j )

  return
end
subroutine yj_check_roman ( y, j, ierror )

!*****************************************************************************80
!
!! YJ_CHECK_ROMAN checks a Roman YJ date.
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
!    Input/output, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no serious error was found in
!    the date, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_roman ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure J is not too small or too big.
!
  call j_borrow_roman ( y, j )

  call j_carry_roman ( y, j )

  return
end
subroutine yj_to_s_common ( y, j, s )

!*****************************************************************************80
!
!! YJ_TO_S_COMMON writes a Common YJ date into a string.
!
!  Format:
!
!    CE YYYY/JJJ
!    BCE YYYY/JJJ
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_common ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( 0 <= y ) then
    s1 = 'CE '
    call i4_to_s_left ( y, s1(4:) )
  else
    s1 = 'BCE '
    call i4_to_s_left (  - y, s1(5:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine yj_to_s_english ( y, j, s )

!*****************************************************************************80
!
!! YJ_TO_S_ENGLISH writes an English YJ date into a string.
!
!  Format:
!
!    BC OS YYYY/JJJ
!    AD OS YYYY/JJJ
!    AD NS YYYY/JJJ
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
!    Input, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_english ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( y < 0 ) then
    s1 = 'BC OS '
    call i4_to_s_left (  - y, s1(7:) )
  else if ( y < 1752 .or. ( y == 1752 .and. j < 278 ) ) then
    s1 = 'AD OS'
    call i4_to_s_left ( y, s1(7:) )
  else
    s1 = 'AD NS'
    call i4_to_s_left ( y, s1(7:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine yj_to_s_gregorian ( y, j, s )

!*****************************************************************************80
!
!! YJ_TO_S_GREGORIAN writes a Gregorian YJ date into a string.
!
!  Format:
!
!    AD YYYY/JJJ
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_gregorian ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( 0 <= y ) then
    s1 = 'AD '
    call i4_to_s_left ( y, s1(4:) )
  else
    s1 = 'BC '
    call i4_to_s_left (  - y, s1(4:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine yj_to_s_hebrew ( y, j, s )

!*****************************************************************************80
!
!! YJ_TO_S_HEBREW writes a Hebrew YJ date into a string.
!
!  Format:
!
!    AM YYYY/JJJ
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
!    Input, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_hebrew ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  s1 = 'AM '
  call i4_to_s_left ( y, s1(4:) )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine yj_to_s_islamic ( y, j, s )

!*****************************************************************************80
!
!! YJ_TO_S_ISLAMIC writes an Islamic YJ date into a string.
!
!  Format:
!
!    AH YYYY/JJJ
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_islamic ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  s1 = 'AH '
  call i4_to_s_left ( y, s1(4:) )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine yj_to_s_julian ( y, j, s )

!*****************************************************************************80
!
!! YJ_TO_S_JULIAN writes a Julian YJ date into a string.
!
!  Format:
!
!    AD YYYY/JJJ
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
!    Input, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_julian ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( 0 <= y ) then
    s1 = 'AD '
    call i4_to_s_left ( y, s1(4:) )
  else
    s1 = 'BC '
    call i4_to_s_left (  - y, s1(4:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine yj_to_s_numeric ( y, j, s )

!*****************************************************************************80
!
!! YJ_TO_S_NUMERIC "prints" a YJ date into a string.
!
!  Format:
!
!    YYYY/JJJ
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_common ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  call i4_to_s_left ( y, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine yj_to_s_republican ( y, j, s )

!*****************************************************************************80
!
!! YJ_TO_S_REPUBLICAN writes a Republican YJ date into a string.
!
!  Format:
!
!    ER YYYY/JJJ
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_republican ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  s1 = 'ER '
  call i4_to_s_left ( y, s1(4:) )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine yj_to_s_roman ( y, j, s )

!*****************************************************************************80
!
!! YJ_TO_S_ROMAN writes a Roman YJ date into a string.
!
!  Format:
!
!    AUC YYYY/JJJ
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, J, the YJ date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_roman ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  s1 = 'AUC '
  call i4_to_s_left ( y, s1(5:) )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine yjf_check_common ( y, j, f, ierror )

!*****************************************************************************80
!
!! YJF_CHECK_COMMON normalizes a Common YJF date.
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
!    Input/output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if there was an error.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y

  call yj_check_common ( y, j, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Force the fraction to lie between 0 and 1.
!
  do while ( f < 0.0D+00 )
    f = f + 1.0D+00
    j = j - 1
  end do

  do while ( 1.0D+00 <= f )
    f = f - 1.0D+00
    j = j + 1
  end do

  return
end
subroutine yjf_check_english ( y, j, f, ierror )

!*****************************************************************************80
!
!! YJF_CHECK_ENGLISH normalizes an English YJF date.
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
!    Input/output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if there was an error.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) y

  call yj_check_english ( y, j, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Force the fraction to lie between 0 and 1.
!
  do while ( f < 0.0D+00 )
    f = f + 1.0D+00
    j = j - 1
  end do

  do while ( 1.0D+00 <= f )
    f = f - 1.0D+00
    j = j + 1
  end do

  return
end
subroutine yjf_compare ( y1, j1, f1, y2, j2, f2, cmp )

!*****************************************************************************80
!
!! YJF_COMPARE compares two YJF dates.
!
!  Discussion:
!
!    The routine is "generic" and does not assume a particular calendar.
!    However, it does assume that the calendar dates are "normalized".
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
!    Input, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1,
!    the first YJF date.
!
!    Input, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2,
!    the second YJF date.
!
!    Output, character CMP:
!    '<' if date 1 precedes date 2;
!    '=' if date 1 equals date 2;
!    '>' if date 1 follows date 2;
!
  implicit none

  character cmp
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  if ( y1 < y2 ) then
    cmp = '<'
  else if ( y1 > y2 ) then
    cmp = '>'
  else

    if ( j1 < j2 ) then
      cmp = '<'
    else if ( j1 > j2 ) then
      cmp = '>'
    else

      if ( f1 < f2 ) then
        cmp = '<'
      else if ( f1 > f2 ) then
        cmp = '>'
      else
        cmp = '='
      end if

    end if

  end if

  return
end
subroutine yjf_dif_common ( y1, j1, f1, y2, j2, f2, days, ierror )

!*****************************************************************************80
!
!! YJF_DIF_COMMON computes day difference between two Common YJF dates.
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
!    Input, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1,
!    the first YJF date.
!
!    Input, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2,
!    the second YJF date.
!
!    Output, real ( kind = 8 ) DAYS, the day difference between the two dates.
!
!    Output, integer ( kind = 4 ) IERROR, is 1 if either date is illegal,
!    0 otherwise.
!
  implicit none

  real ( kind = 8 ) days
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Check the dates.
!
  call yjf_check_common ( y1, j1, f1, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call yjf_check_common ( y2, j2, f2, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call yjf_to_jed_common ( y1, j1, f1, jed1 )

  call yjf_to_jed_common ( y2, j2, f2, jed2 )

  days = jed2 - jed1

  return
end
subroutine yjf_swap ( y1, j1, f1, y2, j2, f2 )

!*****************************************************************************80
!
!! YJF_SWAP swaps the data defining two YJF dates.
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
!    Input/output, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1,
!    the first date.
!
!    Input/output, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2,
!    the second date.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  call i4_swap ( y1, y2 )

  call i4_swap ( j1, j2 )

  call r8_swap ( f1, f2 )

  return
end
subroutine yjf_to_jed_common ( y, j, f, jed )

!*****************************************************************************80
!
!! YJF_TO_JED_COMMON converts a Common YJF date to a JED.
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y1 = y
  j1 = j
  f1 = f
!
!  Check the input.
!
  call yjf_check_common ( y1, j1, f1, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Convert the input.
!
  call yjf_to_ymdf_common ( y1, j1, f1, y2, m2, d2, f2 )

  call ymdf_to_jed_common ( y2, m2, d2, f2, jed )

  return
end
subroutine yjf_to_jed_english ( y, j, f, jed )

!*****************************************************************************80
!
!! YJF_TO_JED_ENGLISH converts an English YJF date to a JED.
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y1 = y
  j1 = j
  f1 = f
!
!  Check the input.
!
  call yj_check_english ( y1, j1, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Convert the input.
!
  call yjf_to_ymdf_english ( y1, j1, f1, y2, m2, d2, f2 )

  call ymdf_to_jed_english ( y2, m2, d2, f2, jed )

  return
end
subroutine yjf_to_jed_gregorian ( y, j, f, jed )

!*****************************************************************************80
!
!! YJF_TO_JED_GREGORIAN converts a Gregorian YJF date to a JED.
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y1 = y
  j1 = j
  f1 = f
!
!  Check the input.
!
  call yj_check_gregorian ( y1, j1, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Convert the input.
!
  call yjf_to_ymdf_gregorian ( y1, j1, f1, y2, m2, d2, f2 )

  call ymdf_to_jed_gregorian ( y2, m2, d2, f2, jed )

  return
end
subroutine yjf_to_jed_hebrew ( y, j, f, jed )

!*****************************************************************************80
!
!! YJF_TO_JED_HEBREW converts a Hebrew YJF date to a JED.
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y1 = y
  j1 = j
  f1 = f
!
!  Check the input.
!
  call yj_check_hebrew ( y1, j1, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Convert the input.
!
  call yjf_to_ymdf_hebrew ( y1, j1, f1, y2, m2, d2, f2 )

  call ymdf_to_jed_hebrew ( y2, m2, d2, f2, jed )

  return
end
subroutine yjf_to_jed_islamic_a ( y, j, f, jed )

!*****************************************************************************80
!
!! YJF_TO_JED_ISLAMIC_A converts an Islamic-A YJF date to a JED.
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y1 = y
  j1 = j
  f1 = f
!
!  Check the input.
!
  call yj_check_islamic ( y1, j1, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Convert the input.
!
  call yjf_to_ymdf_islamic ( y1, j1, f1, y2, m2, d2, f2 )

  call ymdf_to_jed_islamic_a ( y2, m2, d2, f2, jed )

  return
end
subroutine yjf_to_jed_islamic_b ( y, j, f, jed )

!*****************************************************************************80
!
!! YJF_TO_JED_ISLAMIC_B converts an Islamic-B YJF date to a JED.
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y1 = y
  j1 = j
  f1 = f
!
!  Check the input.
!
  call yj_check_islamic ( y1, j1, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Convert the input.
!
  call yjf_to_ymdf_islamic ( y1, j1, f1, y2, m2, d2, f2 )

  call ymdf_to_jed_islamic_b ( y2, m2, d2, f2, jed )

  return
end
subroutine yjf_to_jed_julian ( y, j, f, jed )

!******************************************************************************
!
!! YJF_TO_JED_JULIAN converts a Julian YJF date to a JED.
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y1 = y
  j1 = j
  f1 = f
!
!  Check the input.
!
  call yj_check_julian ( y1, j1, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Convert the input.
!
  call yjf_to_ymdf_julian ( y1, j1, f1, y2, m2, d2, f2 )

  call ymdf_to_jed_julian ( y2, m2, d2, f2, jed )

  return
end
subroutine yjf_to_jed_republican ( y, j, f, jed )

!*****************************************************************************80
!
!! YJF_TO_JED_REPUBLICAN converts a Republican YJF date to a JED.
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y1 = y
  j1 = j
  f1 = f
!
!  Check the input.
!
  call yj_check_republican ( y1, j1, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Convert the input.
!
  call yjf_to_ymdf_republican ( y1, j1, f1, y2, m2, d2, f2 )

  call ymdf_to_jed_republican ( y2, m2, d2, f2, jed )

  return
end
subroutine yjf_to_jed_roman ( y, j, f, jed )

!*****************************************************************************80
!
!! YJF_TO_JED_ROMAN converts a Roman YJF date to a JED.
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y1 = y
  j1 = j
  f1 = f
!
!  Check the input.
!
  call yj_check_roman ( y1, j1, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Convert the input.
!
  call yjf_to_ymdf_roman ( y1, j1, f1, y2, m2, d2, f2 )

  call ymdf_to_jed_roman ( y2, m2, d2, f2, jed )

  return
end
subroutine yjf_to_s_common ( y, j, f, s )

!*****************************************************************************80
!
!! YJF_TO_S_COMMON "prints" a Common YJF date into a string.
!
!  Format:
!
!    CE Y/J.F
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, character ( len = * ) S, contains a representation of the date.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 20 ) s1
  character ( len = 3 ) s2
  character ( len = 8 ) s3
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yjf_check_common ( y, j, f, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( 0 <= y ) then
    s1 = 'CE '
    call i4_to_s_left ( y, s1(4:) )
  else
    s1 = 'BCE '
    call i4_to_s_left (  - y, s1(5:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s1 )

  call frac_to_s ( f, s3 )

  call s_cat ( s1, s3(1:3), s )

  return
end
subroutine yjf_to_s_english ( y, j, f, s )

!*****************************************************************************80
!
!! YJF_TO_S_ENGLISH writes an English YJF date into a string.
!
!  Format:
!
!    BC OS YYYY/JJJ.FF
!    AD OS YYYY/JJJ.FF
!    AD NS YYYY/JJJ.FF
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 20 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yjf_check_english ( y, j, f, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( y < 0 ) then
    s1 = 'BC OS '
    call i4_to_s_left (  -y, s1(7:) )
  else if ( y < 1752 .or. ( y == 1752 .and. j < 278 ) ) then
    s1 = 'AD OS '
    call i4_to_s_left ( y, s1(7:) )
  else
    s1 = 'AD NS '
    call i4_to_s_left ( y, s1(7:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2(1:3) )

  call s_cat ( s1, s2(1:3), s1 )

  call frac_to_s ( f, s2(1:3) )

  call s_cat ( s1, s2(1:3), s1 )

  s = s1

  return
end
subroutine yjf_to_s_gregorian ( y, j, f, s )

!*****************************************************************************80
!
!! YJF_TO_S_GREGORIAN writes a Gregorian YJF date into a string.
!
!  Format:
!
!    AD YYYY/JJJ.FF
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_gregorian ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( 0 <= y ) then
    s1 = 'AD '
    call i4_to_s_left ( y, s1(4:) )
  else
    s1 = 'BC '
    call i4_to_s_left (  - y, s1(4:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  call frac_to_s ( f, s1 )

  call s_cat ( s, s1(1:3), s )

  return
end
subroutine yjf_to_s_hebrew ( y, j, f, s )

!*****************************************************************************80
!
!! YJF_TO_S_HEBREW writes a Hebrew YJF date into a string.
!
!  Format:
!
!    AM YYYY/JJJ.FF
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_hebrew ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  s1 = 'AM '
  call i4_to_s_left ( y, s1(4:) )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  call frac_to_s ( f, s1 )

  call s_cat ( s, s1(1:3), s )

  return
end
subroutine yjf_to_s_islamic ( y, j, f, s )

!*****************************************************************************80
!
!! YJF_TO_S_ISLAMIC writes an Islamic YJF date into a string.
!
!  Format:
!
!    AH YYYY/JJJ.FF
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_islamic ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  s1 = 'AH '
  call i4_to_s_left ( y, s1(4:) )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  call frac_to_s ( f, s1 )

  call s_cat ( s, s1(1:3), s )

  return
end
subroutine yjf_to_s_julian ( y, j, f, s )

!*****************************************************************************80
!
!! YJF_TO_S_JULIAN writes a Julian YJF date into a string.
!
!  Format:
!
!    AD YYYY/JJJ.FF
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_julian ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  s1 = 'AD '
  call i4_to_s_left ( y, s1(4:) )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  call frac_to_s ( f, s1 )

  call s_cat ( s, s1(1:3), s )

  return
end
subroutine yjf_to_s_republican ( y, j, f, s )

!*****************************************************************************80
!
!! YJF_TO_S_REPUBLICAN writes a Republican YJF date into a string.
!
!  Format:
!
!    ER YYYY/JJJ.FF
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_republican ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  s1 = 'ER '
  call i4_to_s_left ( y, s1(4:) )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  call frac_to_s ( f, s1 )

  call s_cat ( s, s1(1:3), s )

  return
end
subroutine yjf_to_s_roman ( y, j, f, s )

!*****************************************************************************80
!
!! YJF_TO_S_ROMAN writes a Roman YJF date into a string.
!
!  Format:
!
!    AUC YYYY/JJJ.FF
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, character ( len = * ) S, the representation of the date.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 11 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y

  call yj_check_roman ( y, j, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  s1 = 'AUC '
  call i4_to_s_left ( y, s1(5:) )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( j, s2 )

  call s_cat ( s1, s2, s )

  call frac_to_s ( f, s1 )

  call s_cat ( s, s1(1:3), s )

  return
end
subroutine yjf_to_weekday_common ( y, j, f, w )

!*****************************************************************************80
!
!! YJF_TO_WEEKDAY_COMMON returns the weekday of a Common YJF date.
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
!    Input, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the YJF date.
!
!    Output, integer ( kind = 4 ) W, is the week day number of the date, with
!    1 for Sunday, through 7 for Saturday.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  integer ( kind = 4 ) j
  real ( kind = 8 ) jed
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y

  call yjf_to_jed_common ( y, j, f, jed )

  call jed_to_weekday ( jed, w, f2 )

  return
end
subroutine yjf_to_ymdf_common ( y1, j1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YJF_TO_YMDF_COMMON converts a Common date from YJF to YMDF format.
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
!    Input, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1, the YJF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y1
  j2 = j1
  f2 = f1
!
!  Check the input.
!
  call yjf_check_common ( y2, j2, f2, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    m2 = 0
    d2 = 0
    f2 = 0.0D+00
    return
  end if
!
!  Convert the input.
!
  d2 = j2
  m2 = 1

  call day_borrow_common ( y2, m2, d2 )

  call day_carry_common ( y2, m2, d2 )

  return
end
subroutine yjf_to_ymdf_english ( y1, j1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YJF_TO_YMDF_ENGLISH converts an English date from YJF to YMDF format.
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
!    Input, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1, the YJF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y1
  j2 = j1
  f2 = f1
!
!  Check the input.
!
  call yj_check_english ( y2, j2, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    m2 = 0
    d2 = 0
    f2 = 0.0D+00
    return
  end if
!
!  Convert the input.
!
  m2 = 1
  d2 = j2

  call day_borrow_english ( y2, m2, d2 )

  call day_carry_english ( y2, m2, d2 )

  return
end
subroutine yjf_to_ymdf_gregorian ( y1, j1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YJF_TO_YMDF_GREGORIAN converts a Gregorian date from YJF to YMDF format.
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
!    Input, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1, the YJF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y1
  j2 = j1
  f2 = f1
!
!  Check the input.
!
  call yj_check_gregorian ( y2, j2, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    m2 = 0
    d2 = 0
    f2 = 0.0D+00
    return
  end if
!
!  Convert the input.
!
  m2 = 1
  d2 = j2

  call day_borrow_gregorian ( y2, m2, d2 )

  call day_carry_gregorian ( y2, m2, d2 )

  return
end
subroutine yjf_to_ymdf_hebrew ( y1, j1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YJF_TO_YMDF_HEBREW converts a YJF to YMDF date, both in the Hebrew calendar.
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
!    Input, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1, the YJF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y1
  j2 = j1
  f2 = f1
!
!  Check the input.
!
  call yj_check_hebrew ( y2, j2, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    m2 = 0
    d2 = 0
    f2 = 0.0D+00
    return
  end if
!
!  Convert the input.
!
  m2 = 1
  d2 = j2

  call day_borrow_hebrew ( y2, m2, d2 )

  call day_carry_hebrew ( y2, m2, d2 )

  return
end
subroutine yjf_to_ymdf_islamic ( y1, j1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YJF_TO_YMDF_ISLAMIC: YJF to YMDF date, both in the Islamic calendar.
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
!    Input, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1, the YJF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y1
  j2 = j1
  f2 = f1
!
!  Check the input.
!
  call yj_check_islamic ( y2, j2, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    m2 = 0
    d2 = 0
    f2 = 0.0D+00
    return
  end if
!
!  Convert the input.
!
  m2 = 1
  d2 = j2

  call day_borrow_islamic ( y2, m2, d2 )

  call day_carry_islamic ( y2, m2, d2 )

  return
end
subroutine yjf_to_ymdf_julian ( y1, j1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YJF_TO_YMDF_JULIAN converts a YJF to YMDF date, both in the Julian calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   16 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1, the YJF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y1
  j2 = j1
  f2 = f1
!
!  Check the input.
!
  call yj_check_julian ( y2, j2, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    m2 = 0
    d2 = 0
    f2 = 0.0D+00
    return
  end if
!
!  Convert the input.
!
  m2 = 1
  d2 = j2

  call day_borrow_julian ( y2, m2, d2 )

  call day_carry_julian ( y2, m2, d2 )

  return
end
subroutine yjf_to_ymdf_republican ( y1, j1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YJF_TO_YMDF_REPUBLICAN: YJF to YMDF date in the Republican calendar.
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
!    Input, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1, the YJF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y1
  j2 = j1
  f2 = f1
!
!  Check the input.
!
  call yj_check_republican ( y2, j2, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    m2 = 0
    d2 = 0
    f2 = 0.0D+00
    return
  end if
!
!  Convert the input.
!
  m2 = 1
  d2 = j2

  call day_borrow_republican ( y2, m2, d2 )

  call day_carry_republican ( y2, m2, d2 )

  return
end
subroutine yjf_to_ymdf_roman ( y1, j1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YJF_TO_YMDF_ROMAN converts a YJF to YMDF date in the Roman calendar.
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
!    Input, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1, the YJF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2, the
!    YMDF date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y1
  j2 = j1
  f2 = f1
!
!  Check the input.
!
  call yj_check_roman ( y2, j2, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    m2 = 0
    d2 = 0
    f2 = 0.0D+00
    return
  end if
!
!  Convert the input.
!
  m2 = 1
  d2 = j2

  call day_borrow_roman ( y2, m2, d2 )

  call day_carry_roman ( y2, m2, d2 )

  return
end
subroutine yjf_to_ymdhms_common ( y1, j1, f1, y2, m2, d2, h2, n2, s2 )

!*****************************************************************************80
!
!! YJF_TO_YMDHMS_COMMON converts a Common YJF date to a YMDHMS date.
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
!    Input, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1, the YJF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, H2, N2, S2, the YMDHMS date.
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) h2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  call yjf_check_common ( y1, j1, f1, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    m2 = 0
    d2 = 0
    h2 = 0
    n2 = 0
    s2 = 0
  end if

  call yjf_to_ymdf_common ( y1, j1, f1, y2, m2, d2, f2 )

  call frac_to_hms ( f2, h2, n2, s2 )

  return
end
subroutine yjf_uniform_common ( y1, j1, f1, y2, j2, f2, seed, y, j, f )

!*****************************************************************************80
!
!! YJF_UNIFORM_COMMON picks a random Common YJF date between two given dates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, J1, real ( kind = 8 ) F1,
!    the first YJF date.
!
!    Input, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2,
!    the second YJF date.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) Y, J, real ( kind = 8 ) F, the random
!    YJF date.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  call yjf_to_jed_common ( y1, j1, f1, jed1 )
  call yjf_to_jed_common ( y2, j2, f2, jed2 )

  jed = r8_uniform_ab ( jed1, jed2, seed )

  call jed_to_yjf_common ( jed, y, j, f )

  return
end
subroutine ym_check_alexandrian ( y, m, ierror )

!*****************************************************************************80
!
!! YM_CHECK_ALEXANDRIAN checks an Alexandrian YM date.
!
!  Discussion:
!
!    If the month is less than 1, then the month is incremented
!    by the number of months in the PREVIOUS year, and the year is
!    decremented by 1.
!
!    If the month is greater than the number of months in the CURRENT year,
!    then the month is decremented by the number of months in the CURRENT year,
!    and the year incremented by 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the
!    date, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_alexandrian ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure the month isn't too small or too big.
!
  call month_borrow_alexandrian ( y, m )

  call month_carry_alexandrian ( y, m )

  return
end
subroutine ym_check_bahai ( y, m, ierror )

!*****************************************************************************80
!
!! YM_CHECK_BAHAI checks a Bahai YM date.
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
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the
!    date, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_bahai
!
!  Check the year.
!
  call y_check_bahai ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  if ( m < 1 .or. year_length_months_bahai ( y ) < m ) then
    ierror = 1
    return
  end if

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
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the date
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
subroutine ym_check_eg_civil ( y, m, ierror )

!*****************************************************************************80
!
!! YM_CHECK_EG_CIVIL checks an Egyptian Civil YM date.
!
!  Discussion:
!
!    If the month is less than 1, then the month is incremented
!    by the number of months in the PREVIOUS year, and the year is
!    decremented by 1.
!
!    If the month is greater than the number of months in the CURRENT year,
!    then the month is decremented by the number of months in the CURRENT year,
!    and the year incremented by 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the date
!    and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_eg_civil ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure the month isn't too small or too big.
!
  call month_borrow_eg_civil ( y, m )

  call month_carry_eg_civil ( y, m )

  return
end
subroutine ym_check_english ( y, m, ierror )

!*****************************************************************************80
!
!! YM_CHECK_ENGLISH checks an English YM date.
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
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the date
!    and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_english ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure the month isn't too small or too big.
!
  call month_borrow_english ( y, m )

  call month_carry_english ( y, m )

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
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the date
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
subroutine ym_check_hebrew ( y, m, ierror )

!*****************************************************************************80
!
!! YM_CHECK_HEBREW checks a Hebrew YM date.
!
!  Discussion:
!
!    If the month is less than 1, then the month is incremented
!    by the number of months in the PREVIOUS year, and the year is
!    decremented by 1.
!
!    If the month is greater than the number of months in the CURRENT year,
!    then the month is decremented by the number of months in the CURRENT year,
!    and the year incremented by 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the date
!    and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_hebrew ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure the month isn't too small or too big.
!
  call month_borrow_hebrew ( y, m )

  call month_carry_hebrew ( y, m )

  return
end
subroutine ym_check_islamic ( y, m, ierror )

!*****************************************************************************80
!
!! YM_CHECK_ISLAMIC checks an Islamic YM date.
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
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the date
!    and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_islamic ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure the month isn't too small or too big.
!
  call month_borrow_islamic ( y, m )

  call month_carry_islamic ( y, m )

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
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the date
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
subroutine ym_check_republican ( y, m, ierror )

!*****************************************************************************80
!
!! YM_CHECK_REPUBLICAN checks a Republican YM date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the date
!    and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_republican ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure the month isn't too small or too big.
!
  call month_borrow_republican ( y, m )

  call month_carry_republican ( y, m )

  return
end
subroutine ym_check_roman ( y, m, ierror )

!*****************************************************************************80
!
!! YM_CHECK_ROMAN checks a Roman YM date.
!
!  Discussion:
!
!    If the month is less than 1, then the month is incremented
!    by the number of months in the PREVIOUS year, and the year is
!    decremented by 1.
!
!    If the month is greater than the number of months in the CURRENT year,
!    then the month is decremented by the number of months in the CURRENT year,
!    and the year incremented by 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the date
!    and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call y_check_roman ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Make sure the month isn't too small or too big.
!
  call month_borrow_roman ( y, m )

  call month_carry_roman ( y, m )

  return
end
subroutine ym_to_decimal ( y, m, yf )

!*****************************************************************************80
!
!! YM_TO_DECIMAL converts a Y/M date to a Decimal Y.F date.
!
!  Discussion:
!
!    Each month is take to be 1/12 of a year long, and the decimal value
!    is returned for the middle of the month.
!
!    1980 January  => 1980.04
!    1980 February => 1980.12
!    1980 March    => 1980.21
!    1980 December => 1980.96
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, the YM date.
!
!    Output, real ( kind = 8 ) YF, the Decimal date.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
  real ( kind = 8 ) yf

  yf = real ( y, kind = 8 ) + real ( 2 * m - 1, kind = 8 ) / 24.0D+00

  return
end
subroutine ymd_check_alexandrian ( y, m, d, ierror )

!*****************************************************************************80
!
!! YMD_CHECK_ALEXANDRIAN checks an Alexandrian YMD date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
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

  ierror = 0
!
!  Check the year.
!
  if ( y <= 0 ) then
    ierror = 1
    return
  end if
!
!  Check the month.
!
  call ym_check_alexandrian ( y, m, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Check the day.
!
  call day_borrow_alexandrian ( y, m, d )

  call day_carry_alexandrian ( y, m, d )

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
subroutine ymd_check_eg_civil ( y, m, d, ierror )

!*****************************************************************************80
!
!! YMD_CHECK_EG_CIVIL checks an Egyptian Civil YMD date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
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

  ierror = 0
!
!  Check the year.
!
  if ( y <= 0 ) then
    ierror = 1
    return
  end if
!
!  Check the month.
!
  call ym_check_eg_civil ( y, m, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Check the day.
!
  call day_borrow_eg_civil ( y, m, d )

  call day_carry_eg_civil ( y, m, d )

  return
end
subroutine ymd_check_english ( y, m, d, ierror )

!*****************************************************************************80
!
!! YMD_CHECK_ENGLISH checks an English YMD date.
!
!  Discussion:
!
!    Certain simple errors in dates will be corrected, such as
!      "31 September 1996"
!    which will become
!      "1 October 1996".
!
!    The routine also knows that in the English calendar, the dates
!    3 September 1752 through 13 September 1752 are illegal.
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
  character ( len = 30 ) s
  integer ( kind = 4 ) y
!
!  Check the month.
!
  call ym_check_english ( y, m, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Check the day.
!
  call day_borrow_english ( y, m, d )

  call day_carry_english ( y, m, d )
!
!  Now make sure that the date does not fall in the
!  Julian-to-Gregorian calendar switchover limbo.
!
  if ( y == 1752 ) then
    if ( m == 9 ) then
      if ( 3 <= d .and. d <= 13 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'YMD_CHECK_ENGLISH - Warning!'
        write ( *, '(a)' ) '  Illegal date!'
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
subroutine ymd_check_hebrew ( y, m, d, ierror )

!*****************************************************************************80
!
!! YMD_CHECK_HEBREW checks a Hebrew YMD date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2000
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
  ierror = 0
!
!  Check the year.
!
  if ( y <= 0 ) then
    ierror = 1
    return
  end if
!
!  Check the month.
!
  call ym_check_hebrew ( y, m, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Check the day.
!
  call day_borrow_hebrew ( y, m, d )

  call day_carry_hebrew ( y, m, d )

  return
end
subroutine ymd_check_islamic ( y, m, d, ierror )

!*****************************************************************************80
!
!! YMD_CHECK_ISLAMIC checks an Islamic YMD date.
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

  ierror = 0
!
!  Check the year.
!
  if ( y <= 0 ) then
    ierror = 1
    return
  end if
!
!  Check the month.
!
  call ym_check_islamic ( y, m, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Check the day.
!
  call day_borrow_islamic ( y, m, d )

  call day_carry_islamic ( y, m, d )

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
subroutine ymd_check_republican ( y, m, d, ierror )

!*****************************************************************************80
!
!! YMD_CHECK_REPUBLICAN checks a Republican YMD date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
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

  ierror = 0
!
!  Check the year.
!
  if ( y <= 0 ) then
    ierror = 1
    return
  end if
!
!  Check the month.
!
  call ym_check_republican ( y, m, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Check the day.
!
  call day_borrow_republican ( y, m, d )

  call day_carry_republican ( y, m, d )

  return
end
subroutine ymd_check_roman ( y, m, d, ierror )

!*****************************************************************************80
!
!! YMD_CHECK_ROMAN checks a Roman YMD date.
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

  ierror = 0
!
!  Check the year.
!
  if ( y <= 0 ) then
    ierror = 1
    return
  end if
!
!  Check the month.
!
  call month_borrow_roman ( y, m )

  call month_carry_roman ( y, m )
!
!  Check the day.
!
  call day_borrow_roman ( y, m, d )

  call day_carry_roman ( y, m, d )

  return
end
subroutine ymd_compare ( y1, m1, d1, y2, m2, d2, cmp )

!*****************************************************************************80
!
!! YMD_COMPARE compares two YMD dates.
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
!    13 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, the first YMD date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2, the second YMD date.
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
        cmp = '='
      end if

    end if

  end if

  return
end
subroutine ymd_dif_common ( y1, m1, d1, y2, m2, d2, days, ierror )

!*****************************************************************************80
!
!! YMD_DIF_COMMON gets the day difference between two Common YMD dates.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, the first YMD date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2, the second YMD date.
!
!    Output, integer ( kind = 4 ) DAYS, the number of days between the dates.
!
!    Output, integer ( kind = 4 ) IERROR, is 1 if either date is illegal,
!    0 otherwise.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) days
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  days = 0
!
!  Check the dates.
!
  call ymd_check_common ( y1, m1, d1, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call ymd_check_common ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call ymd_to_jed_common ( y1, m1, d1, jed1 )

  call ymd_to_jed_common ( y2, m2, d2, jed2 )

  days = nint ( jed2 - jed1 )

  return
end
subroutine ymd_inc_ymd_common ( y1, m1, d1, yn, mn, dn, y2, m2, d2 )

!*****************************************************************************80
!
!! YMD_INC_YMD_COMMON increments a Common YMD date by a YMD increment.
!
!  Discussion:
!
!    You often see on old gravestones statements like
!
!      "Joe Blow died on May 8 1784 aged 38 Years, 7 Months and 5 Days."
!
!    It's not exactly clear how to interpret such a statement, since
!    we can't actually convert 38 Years, 7 Months and 5 Days to a number
!    of days.  (Years and months vary in their day length).  However,
!    we can assume that what was meant was, if you take the year, month
!    and day of Joe Blow's birthday, and you:
!
!      add 38 to the year,
!      add 7 to the month, and if you go past December, subtract 12 and
!        increment the year,
!      add 5 to the day, and if you go past the length of the month,
!        increment the month and decrement the day appropriately.
!
!    Notice, in particular, that if you do the operations in the reverse
!    order, you may get a different answer, since what you do with a large
!    day value depends on the month you assume you are working in.
!
!    Just warning you that this is a poorly posed problem.
!
!    Thanks to Charlie Cullen for pointing out this little problem to me.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, the YMD date.
!
!    Input, integer ( kind = 4 ) YN, MN, DN, the increment to the YMD date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, the incremented YMD date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) dn
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) fn
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) yn
!
!  TEMPORARY
!
  f1 = 0.0D+00
  fn = 0.0D+00

  y2 = y1 + yn
  m2 = m1 + mn
  d2 = d1 + dn
  f2 = f1 + fn

  call ymdf_check_common ( y2, m2, d2, f2, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    m2 = 0
    d2 = 0
    f2 = 0.0D+00
    return
  end if

  return
end
subroutine ymd_to_decimal ( y, m, d, yf )

!*****************************************************************************80
!
!! YMD_TO_DECIMAL converts a Y/M/D date to a Decimal Y.F date.
!
!  Discussion:
!
!    The day is assumed to be at noon.  In other words, 1983 January 1st has
!    a decimal value of 1983 + 0.5 / 365.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, the YMD date.
!
!    Output, real ( kind = 8 ) YF, the Decimal date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) day_max
  integer ( kind = 4 ) days
  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_common
  real ( kind = 8 ) yf
!
!  How many days between Y/1/1 and Y/M/D?
!
  m2 = 1
  d2 = 1
  call ymd_dif_common ( y, m2, d2, y, m, d, days, ierror )
!
!  How many days in this year total?
!
  day_max = year_length_common ( y )
!
!  The decimal part of the year is ( D + 0.5 ) / DMAX.
!
  f = ( real ( days, kind = 8 ) + 0.5D+00 ) &
      / real ( day_max, kind = 8 )

  yf = real ( y, kind = 8 ) + f

  return
end
subroutine ymd_to_jed_common ( y, m, d, jed )

!*****************************************************************************80
!
!! YMD_TO_JED_COMMON converts a Common YMD date to a JED.
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
!    Input, integer ( kind = 4 ) Y, M, D, the YMD date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  character cmp
  integer ( kind = 4 ) d
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
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

  call ymd_check_common ( y1, m1, d1, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if

  y2 = 1582
  m2 = 10
  d2 = 4+1

  call ymd_compare ( y1, m1, d1, y2, m2, d2, cmp )

  if ( cmp == '<' ) then
    call ymd_to_jed_julian ( y1, m1, d1, jed )
    return
  end if
!
!  Use the Gregorian calendar for dates strictly after 1752/9/13.
!
  y2 = 1582
  m2 = 10
  d2 = 15-1

  call ymd_compare ( y1, m1, d1, y2, m2, d2, cmp )

  if ( cmp == '>' ) then
    call ymd_to_jed_gregorian ( y1, m1, d1, jed )
    return
  end if

  jed = -1.0D+00
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'YMD_TO_JED_COMMON - Error!'
  write ( *, '(a)' ) '  Illegal date!'

  return
end
subroutine ymd_to_jed_gregorian ( y, m, d, jed )

!*****************************************************************************80
!
!! YMD_TO_JED_GREGORIAN converts a Gregorian YMD date to a JED.
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
!    Input, integer ( kind = 4 ) Y, M, D, the YMD date.
!
!    Output, real ( kind = 8 ) JED, the corresponding JED.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
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

  g = ( 3 * ( ( y_prime + 184 ) / 100 ) / 4 ) - 38

  jed = real ( j1 + j2 + d_prime - 1401 - g, kind = 8 ) - 0.5D+00

  return
end
subroutine ymd_to_jed_julian ( y, m, d, jed )

!*****************************************************************************80
!
!! YMD_TO_JED_JULIAN converts a Julian YMD date to a JED.
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
!    Input, integer ( kind = 4 ) Y, M, D, the YMD date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
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
  call ymd_check_julian ( y, m, d, ierror )

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
!    volume and issue.
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
subroutine ymd_to_s_alexandrian ( y, m, d, s )

!*****************************************************************************80
!
!! YMD_TO_S_ALEXANDRIAN "prints" an Alexandrian YMD date into a string.
!
!  Format:
!
!    DayNumber MonthName YearNumber AX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2000
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
  character ( len = 2 ) sd
  character ( len = 10 ) sm
  character ( len = 10 ) sy
  character ( len = * ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y
  m2 = m
  d2 = d
!
!  Check the input.
!
  call ymd_check_alexandrian ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  call i4_to_s_left ( y2, sy )

  call month_to_month_name_eg_civil ( m2, sm )

  call i4_to_s_left ( d2, sd )

  call s_cat1 ( sd, sm, s )
  call s_cat1 ( s, sy, s )
  call s_cat1 ( s, 'AX', s )

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
!
!  Check the input.
!
  call ymd_check_common ( y2, m2, d2, ierror )

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

  return
end
subroutine ymd_to_s_eg_civil ( y, m, d, s )

!*****************************************************************************80
!
!! YMD_TO_S_EG_CIVIL "prints" an Egyptian Civil YMD date into a string.
!
!  Format:
!
!    DayNumber MonthName YearNumber EN
!
!  Discussion:
!
!    "EN" stands for the Era of Nabonassar, a Babylonian king who
!    acceded in 747 BC, used by the astronomer Ptolemy to assign
!    an artificial starting year for the Egyptian calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2000
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
  character ( len = 2 ) sd
  character ( len = 10 ) sm
  character ( len = 10 ) sy
  character ( len = * ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y
  m2 = m
  d2 = d
!
!  Check the input.
!
  call ymd_check_eg_civil ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  call i4_to_s_left ( y2, sy )

  call month_to_month_name_eg_civil ( m2, sm )

  call i4_to_s_left ( d2, sd )

  call s_cat1 ( sd, sm, s )
  call s_cat1 ( s, sy, s )
  call s_cat1 ( s, 'EN', s )

  return
end
subroutine ymd_to_s_eg_lunar ( y, m, d, s )

!*****************************************************************************80
!
!! YMD_TO_S_EG_LUNAR "prints" an Egyptian Lunar YMD date into a string.
!
!  Format:
!
!    DayNumber MonthName YearNumber EL
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2000
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
!
  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  character ( len = 2 ) sd
  character ( len = 10 ) sm
  character ( len = 10 ) sy
  character ( len = * ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y
  m2 = m
  d2 = d

  call i4_to_s_left ( y2, sy )

  call month_to_month_name_eg_lunar ( m2, sm )

  call i4_to_s_left ( d2, sd )

  call s_cat1 ( sd, sm, s )
  call s_cat1 ( s, sy, s )
  call s_cat1 ( s, 'EL', s )

  return
end
subroutine ymd_to_s_english ( y, m, d, s )

!*****************************************************************************80
!
!! YMD_TO_S_ENGLISH writes an English YMD date into a string.
!
!  Format:
!
!    BC OS YYYY/MM/DD
!    AD OS YYYY/MM/DD
!    AD NS YYYY/MM/DD
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
!
!  Check the input.
!
  call ymd_check_english ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( y2 < 0 ) then
    s1 = 'BC OS '
    call i4_to_s_left (  - y2, s1(7:) )
  else if ( y2 < 1752 .or. ( y2 == 1752 .and. m2 < 9 ) .or. &
    ( y2 == 1752 .and. m2 == 9 .and. d2 < 3 ) ) then
    s1 = 'AD OS '
    call i4_to_s_left ( y2, s1(7:) )
  else
    s1 = 'AD NS '
    call i4_to_s_left ( y2, s1(7:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( m2, s2 )

  call s_cat ( s1, s2, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( d2, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine ymd_to_s_gregorian ( y, m, d, s )

!*****************************************************************************80
!
!! YMD_TO_S_GREGORIAN writes a Gregorian YMD date into a string.
!
!  Format:
!
!    AD YYYY/MM/DD
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
!
!  Check the input.
!
  call ymd_check_gregorian ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( 0 <= y2 ) then
    s1 = 'AD '
    call i4_to_s_left ( y2, s1(4:) )
  else
    s1 = 'BC '
    call i4_to_s_left (  - y2, s1(4:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( m2, s2 )

  call s_cat ( s1, s2, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( d2, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine ymd_to_s_hebrew ( y, m, d, s )

!*****************************************************************************80
!
!! YMD_TO_S_HEBREW "prints" a Hebrew YMD date into a string.
!
!  Format:
!
!    DayNumber MonthName YearNumber AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2000
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
  character ( len = 2 ) sd
  character ( len = 10 ) sm
  character ( len = 10 ) sy
  character ( len = * ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y
  m2 = m
  d2 = d
!
!  Check the input.
!
  call ymd_check_hebrew ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  call i4_to_s_left ( y2, sy )

  call month_to_month_name_hebrew ( y2, m2, sm )

  call i4_to_s_left ( d2, sd )

  call s_cat1 ( sd, sm, s )
  call s_cat1 ( s, sy, s )
  call s_cat1 ( s, 'AM', s )

  return
end
subroutine ymd_to_s_islamic ( y, m, d, s )

!*****************************************************************************80
!
!! YMD_TO_S_ISLAMIC writes an Islamic YMD date into a string.
!
!  Format:
!
!    DayNumber MonthName YearNumber AH
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2000
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
  character ( len = 2 ) sd
  character ( len = 10 ) sm
  character ( len = 10 ) sy
  character ( len = * ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Copy the input.
!
  y2 = y
  m2 = m
  d2 = d
!
!  Check the input.
!
  call ymd_check_islamic ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  call i4_to_s_left ( y2, sy )

  call month_to_month_name_islamic ( m2, sm )

  call i4_to_s_left ( d2, sd )

  call s_cat1 ( sd, sm, s )
  call s_cat1 ( s, sy, s )
  call s_cat1 ( s, 'AH', s )

  return
end
subroutine ymd_to_s_julian ( y, m, d, s )

!*****************************************************************************80
!
!! YMD_TO_S_JULIAN writes a Julian YMD date into a string.
!
!  Format:
!
!    AD YYYY/MM/DD
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
!
!  Check the input.
!
  call ymd_check_julian ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( 0 <= y2 ) then
    s1 = 'AD '
    call i4_to_s_left ( y2, s1(4:) )
  else
    s1 = 'BC '
    call i4_to_s_left (  - y2, s1(4:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( m2, s2 )

  call s_cat ( s1, s2, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( d2, s2 )

  call s_cat ( s1, s2, s )

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
subroutine ymd_to_s_republican ( y, m, d, s )

!*****************************************************************************80
!
!! YMD_TO_S_REPUBLICAN writes a Republican YMD date into a string.
!
!  Format:
!
!    ER YYYY/MM/DD
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
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
!
!  Check the input.
!
  call ymd_check_republican ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( 0 <= y2 ) then
    s1 = 'ER '
    call i4_to_s_left ( y2, s1(4:) )
  else
    s1 = '-ER '
    call i4_to_s_left (  - y2, s1(4:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( m2, s2 )

  call s_cat ( s1, s2, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( d2, s2 )

  call s_cat ( s1, s2, s )

  return
end
subroutine ymd_to_s_roman ( y, m, d, s )

!*****************************************************************************80
!
!! YMD_TO_S_ROMAN writes a Roman YMD date into a string.
!
!  Example:
!
!     Y  M   D  S
!    --  -  --  -----------------------------------
!    56  4   1  Kalends Aprilis DVI AUC
!    56  4   2  Ante diem iv Nones Aprilis DVI AUC
!    56  4   3  Ante diem iii Nones Aprilis DVI AUC
!    56  4   4  Pridie Nones Aprilis DVI AUC
!    56  4   5  Nones Aprilis DVI AUC
!    56  4   6  Ante diem viii Ides Aprilis DVI AUC
!    56  4   7  Ante diem vii Ides Aprilis DVI AUC
!    56  4   8  Ante diem vi Ides Aprilis DVI AUC
!    56  4   9  Ante diem v Ides Aprilis DVI AUC
!    56  4  10  Ante diem iv Ides Aprilis DVI AUC
!    56  4  11  Ante diem iii Ides Aprilis DVI AUC
!    56  4  12  Pridie Ides Aprilis DVI AUC
!    56  4  13  Ides Aprilis DVI AUC
!    56  4  14  Ante diem xvii Kalends Maius DVI AUC
!    56  4  15  Ante diem xvi Kalends Maius DVI AUC
!    ...
!    56  4  28  Ante diem iv Kalends Maius DVI AUC
!    56  4  29  Ante diem iii Kalends Maius DVI AUC
!    56  4  30  Pridie Kalends Maius DVI AUC
!
!  Discussion:
!
!    "AUC" means "ab urbe condita", or "from the founding of the city".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, the YMD date.
!
!    Output, character ( len = * ) S, a string representing the date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ides
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) jday
  integer ( kind = 4 ) last
  character ( len = 10 ) lower
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ) month_length_roman
  integer ( kind = 4 ) nones
  character ( len = * ) s
  character ( len = 10 ) s_day
  character ( len = 15 ) s_month
  character ( len = 15 ) s_month_next
  character ( len = 10 ) s_year
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical   year_is_leap_roman
!
!  Copy the input.
!
  y2 = y
  m2 = m
  d2 = d
!
!  Check the input.
!
  call ymd_check_roman ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  call month_to_month_name_roman ( m2, s_month )
!
!  Get the next month's name.
!
  m3 = i4_wrap ( m2+1, 1, 12 )
  call month_to_month_name_roman ( m3, s_month_next )

  call month_to_nones_roman ( m2, nones )
  call month_to_ides_roman ( m2, ides )
  last = month_length_roman ( y2, m2 )

  if ( d2 == 1 ) then
    s = 'Kalends ' // s_month
  else if ( d2 < nones - 1 ) then
    jday = nones + 1 - d2
    call i4_to_roman ( jday, s_day )
    s_day = lower ( s_day )
    s = 'Ante diem ' // trim ( s_day ) // ' Nones ' // s_month
  else if ( d2 == nones - 1 ) then
    s = 'Pridie Nones ' // s_month
  else if ( d2 == nones ) then
    s = 'Nones ' // s_month
  else if ( d2 < ides - 1 ) then
    jday = ides + 1 - d2
    call i4_to_roman ( jday, s_day )
    s_day = lower ( s_day )
    s = 'Ante diem ' // trim ( s_day ) // ' Ides ' // s_month
  else if ( d2 == ides - 1 ) then
    s = 'Pridie Ides ' // s_month
  else if ( d2 == ides ) then
    s = 'Ides ' // s_month

  else if ( m2 == 2 .and. year_is_leap_roman ( y2 ) ) then

    if ( d2 < 25 ) then
      jday = last + 1 - d2
      call i4_to_roman ( jday, s_day )
      s_day = lower ( s_day )
      s = 'Ante diem ' // trim ( s_day ) // ' Kalends ' // s_month_next
    else if ( d2 == 25 ) then
      jday = last + 2 - d2
      call i4_to_roman ( jday, s_day )
      s_day = lower ( s_day )
      s = 'Ante diem Bis ' // trim ( s_day ) // ' Kalends ' // s_month_next
    else if ( d2 < last ) then
      jday = last + 2 - d2
      call i4_to_roman ( jday, s_day )
      s_day = lower ( s_day )
      s = 'Ante diem ' // trim ( s_day ) // ' Kalends ' // s_month_next
    else
      s = 'Pridie Kalends ' // s_month_next
    end if

  else if ( d2 < last ) then
    jday = last + 2 - d2
    call i4_to_roman ( jday, s_day )
    s_day = lower ( s_day )
    s = 'Ante diem ' // trim ( s_day ) // ' Kalends ' // s_month_next
  else
    s = 'Pridie Kalends ' // s_month_next
  end if

  call i4_to_roman ( y2, s_year )

  call s_cat1 ( s, s_year, s )
  call s_cat ( s, ' AUC', s )

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
!    09 June 2012
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

  f = 0.0D+00
  call ymdf_to_jed_common ( y, m, d, f, jed )

  call jed_to_weekday ( jed, w, f2 )

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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the
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
subroutine ymdf_check_english ( y, m, d, f, ierror )

!*****************************************************************************80
!
!! YMDF_CHECK_ENGLISH checks an English YMDF date.
!
!  Discussion:
!
!    Certain simple errors in dates will be corrected, such as
!      "31 September 1996"
!    which will become
!      "1 October 1996".
!
!    The routine also knows that in the English calendar, the dates
!    3 September 1752 through 13 September 1752 are illegal.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the
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

  call ymd_check_english ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call frac_borrow_english ( y, m, d, f )

  call frac_carry_english ( y, m, d, f )

  return
end
subroutine ymdf_check_gregorian ( y, m, d, f, ierror )

!*****************************************************************************80
!
!! YMDF_CHECK_GREGORIAN checks a Gregorian YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the
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

  call ymd_check_gregorian ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call frac_borrow_gregorian ( y, m, d, f )

  call frac_carry_gregorian ( y, m, d, f )

  return
end
subroutine ymdf_check_hebrew ( y, m, d, f, ierror )

!*****************************************************************************80
!
!! YMDF_CHECK_HEBREW checks a Hebrew YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the
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

  call ymd_check_hebrew ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call frac_borrow_hebrew ( y, m, d, f )

  call frac_carry_hebrew ( y, m, d, f )

  return
end
subroutine ymdf_check_islamic ( y, m, d, f, ierror )

!*****************************************************************************80
!
!! YMDF_CHECK_ISLAMIC checks an Islamic YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the
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

  call ymd_check_islamic ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call frac_borrow_islamic ( y, m, d, f )

  call frac_carry_islamic ( y, m, d, f )

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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the
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
subroutine ymdf_check_republican ( y, m, d, f, ierror )

!*****************************************************************************80
!
!! YMDF_CHECK_REPUBLICAN checks a Republican YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the
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

  call ymd_check_republican ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call frac_borrow_republican ( y, m, d, f )

  call frac_carry_republican ( y, m, d, f )

  return
end
subroutine ymdf_check_roman ( y, m, d, f, ierror )

!*****************************************************************************80
!
!! YMDF_CHECK_ROMAN checks a Roman YMDF date.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the
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

  call ymd_check_roman ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call frac_borrow_roman ( y, m, d, f )

  call frac_carry_roman ( y, m, d, f )

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
subroutine ymdf_dif_common ( y1, m1, d1, f1, y2, m2, d2, f2, days, ierror )

!*****************************************************************************80
!
!! YMDF_DIF_COMMON gets the day difference between two Common YMDF dates.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the first YMDF date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the second YMDF date.
!
!    Output, real ( kind = 8 ) DAYS, the number of days between the two dates.
!
!    Output, integer ( kind = 4 ) IERROR, is 1 if either date is illegal,
!    0 otherwise.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  days = 0.0D+00
!
!  Check the dates.
!
  call ymdf_check_common ( y1, m1, d1, f1, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call ymdf_check_common ( y2, m2, d2, f2, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call ymdf_to_jed_common ( y1, m1, d1, f1, jed1 )

  call ymdf_to_jed_common ( y2, m2, d2, f2, jed2 )

  days = jed2 - jed1

  return
end
subroutine ymdf_dif_english ( y1, m1, d1, f1, y2, m2, d2, f2, days, ierror )

!*****************************************************************************80
!
!! YMDF_DIF_ENGLISH gets the day difference between two English YMDF dates.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, the first YMDF date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2, the second YMDF date.
!
!    Output, real ( kind = 8 ) DAYS, the number of days between the two dates.
!
!    Output, integer ( kind = 4 ) IERROR, is 1 if either date is illegal,
!    0 otherwise.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  days = 0.0D+00
!
!  Check the dates.
!
  call ymdf_check_english ( y1, m1, d1, f1, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call ymdf_check_english ( y2, m2, d2, f2, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call ymdf_to_jed_english ( y1, m1, d1, f1, jed1 )

  call ymdf_to_jed_english ( y2, m2, d2, f2, jed2 )

  days = jed2 - jed1

  return
end
subroutine ymdf_dif_ymdf_common ( y1, m1, d1, f1, y2, m2, d2, f2, yn, mn, &
  dn, fn, ierror )

!*****************************************************************************80
!
!! YMDF_DIF_YMDF_COMMON gets the YMDF difference between two Common YMDF dates.
!
!  Discussion:
!
!    This difference is not well defined.  A reasonable way to define this
!    difference is:
!
!      Use Y1, M1, D1, F1 as a base,
!
!      Increment Y1 by 1 repeatedly, until your date is less than
!      a year before Y2/M2/D2/F2.
!
!      Increment M1 by 1 repeatedly, until your date is less than a
!      month before Y2/M2/D2/F2.
!
!      Increment D1 by 1 repeatedly, until your data is less than a
!      day before Y2/M2/D2/F2.
!
!      Measure the fractional day difference.
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
!    Output, integer ( kind = 4 ) YN, MN, DN, real ( kind = 8 ) FN,
!    the difference in years, months, days, and fractional days from the
!    first date to the second.
!
!    Output, integer ( kind = 4 ) IERROR, is 1 if either date is illegal,
!    0 otherwise.
!
  implicit none

  character cmp
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) dn
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) fn
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3
  integer ( kind = 4 ) yn
!
!  Check the dates.
!
  call ymdf_check_common ( y1, m1, d1, f1, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call ymdf_check_common ( y2, m2, d2, f2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Compare the dates.
!
  call ymdf_compare ( y1, m1, d1, f1, y2, m2, d2, f2, cmp )
!
!  We swap dates, if necessary, so that date 1 is never greater
!  than date 2.
!
  if ( cmp == '=' ) then
    yn = 0
    mn = 0
    dn = 0
    fn = 0.0D+00
    return
  else if ( cmp == '>' ) then
    call ymdf_swap ( y1, m1, d1, f1, y2, m2, d2, f2 )
  end if
!
!  Year difference.
!
  yn = y2 - y1
!
!  Month difference.
!
  y3 = y2

  if ( m2 < m1 ) then
    yn = yn - 1
    mn = m2 - m1 + 12
    y3 = y2 - 1
  else
    mn = m2 - m1
  end if

  m3 = m2
!
!  Day difference.
!
  if ( d2 < d1 ) then
    mn = mn - 1
    m3 = m2 - 1
    if ( m3 == 0 ) then
      m3 = 12
      y3 = y3 - 1
    end if
    dn = d2 - d1 + month_length_common ( y3, m3 )
  else
    dn = d2 - d1
  end if
!
!  Fractional difference.
!  THERE'S MORE TO THIS CODE, AFTER ALL, WHAT IF DN < 0?
!
  if ( f2 < f1 ) then
    dn = dn - 1
  end if

  fn = f2 - f1

  return
end
subroutine ymdf_inc_common ( y1, m1, d1, f1, days, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_INC_COMMON increments a Common YMDF date by DAYS days.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1, the
!    YMDF date.
!
!    Input, real ( kind = 8 ) DAYS, the number of days to advance the date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2, the
!    incremented YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the parameters.
!
  y2 = y1
  m2 = m1
  d2 = d1
  f2 = f1 + days
!
!  Check the parameters.
!
  call ymdf_check_common ( y2, m2, d2, f2, ierror )

  return
end
subroutine ymdf_inc_english ( y1, m1, d1, f1, days, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_INC_ENGLISH increments an English YMDF date by DAYS days.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Input, real ( kind = 8 ) DAYS, the number of days to advance the date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2, the
!    incremented YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the parameters.
!
  y2 = y1
  m2 = m1
  d2 = d1
  f2 = f1 + days
!
!  Check the parameters.
!
  call ymdf_check_english ( y2, m2, d2, f2, ierror )

  return
end
subroutine ymdf_inc_gregorian ( y1, m1, d1, f1, days, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_INC_GREGORIAN increments a Gregorian YMDF date by DAYS days.
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
!    Input, integer ( kind = 4 ) Y1,  M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Input, real ( kind = 8 ) DAYS, the number of days to advance the date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the incremented YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the parameters.
!
  y2 = y1
  m2 = m1
  d2 = d1
  f2 = f1 + days
!
!  Check the parameters.
!
  call ymdf_check_gregorian ( y2, m2, d2, f2, ierror )

  return
end
subroutine ymdf_inc_hebrew ( y1, m1, d1, f1, days, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_INC_HEBREW increments a Hebrew YMDF date by DAYS days.
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
!    Input, integer ( kind = 4 ) Y1,  M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Input, real ( kind = 8 ) DAYS, the number of days to advance the date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the incremented YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the parameters.
!
  y2 = y1
  m2 = m1
  d2 = d1
  f2 = f1 + days
!
!  Check the parameters.
!
  call ymdf_check_hebrew ( y2, m2, d2, f2, ierror )

  return
end
subroutine ymdf_inc_islamic ( y1, m1, d1, f1, days, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_INC_ISLAMIC increments an Islamic YMDF date by DAYS days.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Input, real ( kind = 8 ) DAYS, the number of days to advance the date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the incremented YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the parameters.
!
  y2 = y1
  m2 = m1
  d2 = d1
  f2 = f1 + days
!
!  Check the parameters.
!
  call ymdf_check_islamic ( y2, m2, d2, f2, ierror )

  return
end
subroutine ymdf_inc_julian ( y1, m1, d1, f1, days, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_INC_JULIAN increments a Julian YMDF date by DAYS days.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Input, real ( kind = 8 ) DAYS, the number of days to advance the date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the incremented YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the parameters.
!
  y2 = y1
  m2 = m1
  d2 = d1
  f2 = f1 + days
!
!  Check the parameters.
!
  call ymdf_check_julian ( y2, m2, d2, f2, ierror )

  return
end
subroutine ymdf_inc_republican ( y1, m1, d1, f1, days, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_INC_REPUBLICAN increments a Republican YMDF date by DAYS days.
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
!    Input, integer ( kind = 4 ) Y1, M1,  D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Input, real ( kind = 8 ) DAYS, the number of days to advance the date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the incremented YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the parameters.
!
  y2 = y1
  m2 = m1
  d2 = d1
  f2 = f1 + days
!
!  Check the parameters.
!
  call ymdf_check_republican ( y2, m2, d2, f2, ierror )

  return
end
subroutine ymdf_inc_roman ( y1, m1, d1, f1, days, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_INC_ROMAN increments a Roman YMDF date by a DAYS days.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Input, real ( kind = 8 ) DAYS, the number of days to advance the date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the incremented YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Copy the parameters.
!
  y2 = y1
  m2 = m1
  d2 = d1
  f2 = f1 + days
!
!  Check the parameters.
!
  call ymdf_check_roman ( y2, m2, d2, f2, ierror )

  return
end
subroutine ymdf_next_common ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_NEXT_COMMON returns the Common YMDF date of the next day.
!
!  Discussion:
!
!    The routine knows that in the Common calendar, the day after
!    4 October 1582 was 15 October 1582.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    tomorrow's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 + 1
  f2 = f1

  call day_carry_common ( y2, m2, d2 )

  return
end
subroutine ymdf_next_english ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_NEXT_ENGLISH returns the English YMD date of the next day.
!
!  Discussion:
!
!    The routine knows that in the English calendar,
!    the day after 2 September 1752 was 14 September 1752.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    tomorrow's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 + 1
  f2 = f1

  call day_carry_english ( y2, m2, d2 )

  return
end
subroutine ymdf_next_gregorian ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_NEXT_GREGORIAN returns the Gregorian YMDF date of the next day.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    tomorrow's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 + 1
  f2 = f1

  call day_carry_gregorian ( y2, m2, d2 )

  return
end
subroutine ymdf_next_hebrew ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_NEXT_HEBREW returns the Hebrew YMDF date of the next day.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    tomorrow's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 + 1
  f2 = f1

  call day_carry_hebrew ( y2, m2, d2 )

  return
end
subroutine ymdf_next_islamic ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_NEXT_ISLAMIC returns the Islamic YMDF date of the next day.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    tomorrow's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 + 1
  f2 = f1

  call day_carry_islamic ( y2, m2, d2 )

  return
end
subroutine ymdf_next_julian ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_NEXT_JULIAN returns the Julian YMDF date of the next day.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    tomorrow's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 + 1
  f2 = f1

  call day_carry_julian ( y2, m2, d2 )

  return
end
subroutine ymdf_next_republican ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_NEXT_REPUBLICAN returns the Republican YMD date of the next day.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    tomorrow's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 + 1
  f2 = f1

  call day_carry_republican ( y2, m2, d2 )

  return
end
subroutine ymdf_next_roman ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_NEXT_ROMAN returns the Roman YMDF date of the next day.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    tomorrow's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 + 1
  f2 = f1

  call day_carry_roman ( y2, m2, d2 )

  return
end
subroutine ymdf_prev_common ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_PREV_COMMON returns the Common YMDF date of the previous day.
!
!  Discussion:
!
!    The routine knows that in the Common calendar, the day before
!    15 October 1582 was 4 October 1582.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   16 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    yesterday's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 - 1
  f2 = f1

  call day_borrow_common ( y2, m2, d2 )

  return
end
subroutine ymdf_prev_english ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_PREV_ENGLISH returns the English YMDF date of the previous day.
!
!  Discussion:
!
!    The routine knows that in the English calendar,
!    the day before 14 September 1752 was 2 September 1752.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    yesterday's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 - 1
  f2 = f1

  call day_borrow_english ( y2, m2, d2 )

  return
end
subroutine ymdf_prev_gregorian ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_PREV_GREGORIAN returns the Gregorian YMDF date of the previous day.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    yesterday's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 - 1
  f2 = f1

  call day_borrow_gregorian ( y2, m2, d2 )

  return
end
subroutine ymdf_prev_hebrew ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_PREV_HEBREW returns the Hebrew YMDF date of the previous day.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    yesterday's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 - 1
  f2 = f1

  call day_borrow_hebrew ( y2, m2, d2 )

  return
end
subroutine ymdf_prev_islamic ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_PREV_ISLAMIC returns the Islamic YMDF date of the previous day.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    yesterday's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 - 1
  f2 = f1

  call day_borrow_islamic ( y2, m2, d2 )

  return
end
subroutine ymdf_prev_julian ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_PREV_JULIAN returns the Julian YMDF date of the previous day.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    yesterday's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 - 1
  f2 = f1

  call day_borrow_julian ( y2, m2, d2 )

  return
end
subroutine ymdf_prev_republican ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_PREV_REPUBLICAN returns the Republican YMDF date of the previous day.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    yesterday's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 - 1
  f2 = f1

  call day_borrow_republican ( y2, m2, d2 )

  return
end
subroutine ymdf_prev_roman ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_PREV_ROMAN returns the Roman YMDF date of the previous day.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    yesterday's YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  m2 = m1
  d2 = d1 - 1
  f2 = f1

  call day_borrow_roman ( y2, m2, d2 )

  return
end
subroutine ymdf_swap ( y1, m1, d1, f1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDF_SWAP swaps two YMDF dates.
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
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  call i4_swap ( y1, y2 )

  call i4_swap ( m1, m2 )

  call i4_swap ( d1, d2 )

  call r8_swap ( f1, f2 )

  return
end
subroutine ymdf_to_jed_alexandrian ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_ALEXANDRIAN converts an Alexandrian YMDF date to a JED.
!
!  Discussion:
!
!    This code needs to be adjusted to fit the Alexandrian model.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 4690 - ( 13 - m ) / 13
  m_prime = mod ( m + 12, 13 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  jed = real ( &
    ( 1461 * y_prime ) / 4 + 30 * m_prime + d_prime - 124, kind = 8 ) - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_armenian ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_ARMENIAN converts an Armenian YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none
!
  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 5268 - ( 13 - m ) / 13
  m_prime = mod ( m + 12, 13 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  jed = real ( 365 * y_prime + 30 * m_prime + d_prime - 317, kind = 8 ) &
    - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_bahai ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_BAHAI converts a Bahai YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none
!
  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  integer ( kind = 4 ) g
  real ( kind = 8 ) jed
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 6560 - ( 39 - m ) / 20
  m_prime = mod ( m, 20 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  j1 = ( 1461 * y_prime ) / 4

  j2 = 19 * m_prime

  g = ( 3 * ( ( y_prime + 184 ) / 100 ) / 4 ) - 50
  jed = real ( j1 + j2 + d_prime - 1412 - g, kind = 8 ) - 0.5D+00
  jed = jed + f

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
subroutine ymdf_to_jed_coptic ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_COPTIC converts a Coptic YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime

!  Convert the calendar date to a computational date.
!
  y_prime = y + 4996 - ( 13 - m ) / 13
  m_prime = mod ( m + 12, 13 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  jed = real ( &
    ( 1461 * y_prime ) / 4 + 30 * m_prime + d_prime - 124, kind = 8 ) - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_eg_civil ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_EG_CIVIL converts an Egyptian Civil YMDF date to a JED.
!
!  Discussion:
!
!    The Egyptian Civil calendar used a year of 365 days.  The year comprised
!    12 months of 30 days, with 5 epagomenal days occurring at the end of
!    the year.  Since the observed year is about 365.25 days long, and no
!    attempt was made to adjust the Egyptian Civil year to the observed year,
!    the calendar dates gradually drifted with respect to the observed dates.
!
!    The epoch or first day of the Egyptian Civil calendar is taken as
!    JED = 1448638.5.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 3968 - ( 13 - m ) / 13
  m_prime = mod ( m + 12, 13 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  jed = real ( 365 * y_prime + 30 * m_prime + d_prime - 47 + 1, kind = 8 ) &
    - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_eg_lunar ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_EG_LUNAR converts an Egyptian Lunar YMDF date to a JED.
!
!  Discussion:
!
!    Count
!      the days up to the day before the start of the calendar,
!      the days in the current month,
!      the 29 days guaranteed in the previous months of this year,
!      the (months/2) 30th days in the previous months of this year,
!      the 354 days guaranteed in each of the previous years,
!      the extra leap days in the preceding years,
!      the extra 30 days in the leap months in the preceding years.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  call epoch_to_jed_eg_lunar ( jed_epoch )

  jed = jed_epoch + real ( - 1 + d + 29 * ( m - 1 ) + ( m - 1 ) / 2 &
    + 354 * ( y - 1 ) + ( y - 1 ) / 5 &
    + 30 * ( ( ( y - 1 ) / 25 ) * 9 + ( mod ( ( y - 1 ), 25 ) + 2 ) / 3 ), &
    kind = 8 )

  jed = jed + f

  return
end
subroutine ymdf_to_jed_english ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_ENGLISH converts an English YMDF date to a JED.
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

  character ( len = 1 ) cmp
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
!  Check the date.
!
  call ymd_check_english ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Use the Julian Calendar for dates strictly before 1752/9/3.
!
  y1 = 1752
  m1 = 9
  d1 = 3
  f1 = 0.0D+00

  call ymdf_compare ( y, m, d, f, y1, m1, d1, f1, cmp )

  if ( cmp == '<' ) then
    call ymdf_to_jed_julian ( y, m, d, f, jed )
    return
  end if
!
!  Use the Gregorian calendar for dates strictly after 1752/9/13.
!
  y2 = 1752
  m2 = 9
  d2 = 13
  f2 = 0.0D+00

  call ymdf_compare ( y, m, d, f, y2, m2, d2, f2, cmp )

  if ( cmp == '>' ) then
    call ymdf_to_jed_gregorian ( y, m, d, f, jed )
    return
  end if
!
!  Error return if the date falls between the transition dates.
!
  jed = -1.0D+00

  return
end
subroutine ymdf_to_jed_ethiopian ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_ETHIOPIAN converts an Ethiopian YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 4720 - ( 13 - m ) / 13
  m_prime = mod ( m + 12, 13 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  j1 = ( 1461 * y_prime ) / 4

  j2 = 30 * m_prime

  jed = real ( j1 + j2 + d_prime - 124, kind = 8 ) - 0.5D+00
  jed = jed + f

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

  g = ( 3 * ( ( y_prime + 184 ) / 100 ) / 4 ) - 38

  jed = real ( j1 + j2 + d_prime - 1401 - g, kind = 8 ) - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_hebrew ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_HEBREW converts a Hebrew YMDF date to a JED.
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
!    Algorithm J,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 334.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, real ( kind = 8 ) JED, the corresponding JED.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) month_length_hebrew
  integer ( kind = 4 ) y
!
!  Check the date.
!
  call ymd_check_hebrew ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YMDF_TO_JED_HEBREW - Fatal error!'
    write ( *, '(a)' ) '  Illegal date!'
    write ( *, '(a,3i6)' ) '  Y/M/D = ', y, m, d
    jed = -1.0D+00
    return
  end if
!
!  Determine the JED of the beginning of the year.
!
  call new_year_to_jed_hebrew ( y, jed )
!
!  Work through the preceding months.
!
  do m2 = 1, m - 1
    jed = jed + real ( month_length_hebrew ( y, m2 ), kind = 8 )
  end do
!
!  Add on the days.
!
  jed = jed + real ( d - 1, kind = 8 )
  jed = jed + f

  return
end
subroutine ymdf_to_jed_hindu_solar ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_HINDU_SOLAR converts a Hindu solar YMDF date to a JED.
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
!    Edward Reingold, Nachum Dershowitz, Stewart Clamen,
!    Calendrical Calculations, II: Three Historical Calendars,
!    Software - Practice and Experience,
!    Volume 23, Number 4, pages 383-404, April 1993.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  real ( kind = 8 ) month_length_hindu_solar
  integer ( kind = 4 ) y
  real ( kind = 8 ) year_length_hindu_solar

  call epoch_to_jed_hindu_solar ( jed_epoch )

  jed = jed_epoch + &
      real ( d - 1, kind = 8 ) &
    + real ( m - 1, kind = 8 ) * month_length_hindu_solar ( ) &
    + real ( y, kind = 8 ) * year_length_hindu_solar ( )

  jed = jed + f

  return
end
subroutine ymdf_to_jed_islamic_a ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_ISLAMIC_A converts an Islamic A YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
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
  integer ( kind = 4 ) y_prime
!
!  Check the date.
!
  call ymd_check_islamic ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 5519 - ( 12 - m ) / 12
  m_prime = mod ( m + 11, 12 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  j1 = ( 10631 * y_prime + 14 ) / 30

  j2 = ( 2951 * m_prime + 51 ) / 100

  jed = real ( j1 + j2 + d_prime - 7665, kind = 8 ) - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_islamic_a2 ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_ISLAMIC_A2 converts an Islamic A YMDF date to a JED.
!
!  Discussion:
!
!    The algorithm has the beauty of being comprehensible!
!
!    Count the days up to the day before the start of the calendar,
!    plus the days in the current month, the 29 days guaranteed
!    in the previous months of this year, the (months/2) 30th days,
!    the 354 days in each of the previous years, plus the total number
!    of leap days in the preceding years.
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
!    Edward Reingold, Nachum Dershowitz,
!    Calendrical Calculations I,
!    Software - Practice and Experience,
!    Volume 20, Number 9, September 1990, pages 899-928.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the date.
!
  call ymd_check_islamic ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if

  call epoch_to_jed_islamic_a ( jed_epoch )

  jed = jed_epoch + real ( - 1 + d + 29 * ( m - 1 ) + ( m / 2 ) &
    + 354 * ( y - 1 ) + ( 11 * y + 3 ) / 30, kind = 8 )
  jed = jed + f

  return
end
subroutine ymdf_to_jed_islamic_b ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_ISLAMIC_B converts an Islamic B YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
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
  integer ( kind = 4 ) y_prime
!
!  Check the date.
!
  call ymd_check_islamic ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 5519 - ( 12 - m ) / 12
  m_prime = mod ( m + 11, 12 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  j1 = ( 10631 * y_prime + 14 ) / 30

  j2 = ( 2951 * m_prime + 51 ) / 100

  jed = real ( j1 + j2 + d_prime - 7664, kind = 8 ) - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_jelali ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_JELALI converts a Jelali YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  call epoch_to_jed_jelali ( jed_epoch )

  jed = jed_epoch + real ( ( d - 1 ) + 30 * ( m - 1 ) + 365 * ( y - 1 ) &
    + ( y - 1 ) / 4, kind = 8 )
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
subroutine ymdf_to_jed_julian2 ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_JULIAN2 converts a Julian YMDF date to a JED.
!
!  Example:
!
!          Y  M  D          JED
!    --------------     -------
!    BC 4713  1  1            0
!    AD    1  1  1      1721424
!    AD 1844  5 11      2394710
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 2008
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

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) days_before_month_julian
  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3

  y2 = y
  m2 = m
  d2 = d

  call ymd_check_julian ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Account for the missing year 0 by moving negative years up one.
!
  call y_common_to_astronomical ( y2, y3 )
!
!  The JED is the number of days in years past, plus the number of days in
!  the previous months this year, plus the number of days.
!
  jed = real ( ( ( 1461 * ( y3 + 4715 ) ) / 4 ) - 1095 &
    + days_before_month_julian ( y2, m2 ) + d2 - 1, kind = 8 ) - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_khwarizmian ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_KHWARIZMIAN converts a Khwarizmian YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 5348 - ( 13 - m ) / 13
  m_prime = mod ( m + 12, 13 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  jed = real ( 365 * y_prime + 30 * m_prime + d_prime - 317, kind = 8 ) &
    - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_macedonian ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_MACEDONIAN converts a Macedonian YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 4405 - ( 18 - m ) / 12
  m_prime = mod ( m + 5, 12 )
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
subroutine ymdf_to_jed_persian ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_PERSIAN converts a Persian YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 5348 - ( 22 - m ) / 13
  m_prime = mod ( m + 3, 13 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  jed = real ( 365 * y_prime + 30 * m_prime + d_prime - 77, kind = 8 ) - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_republican ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_REPUBLICAN converts a Republican YMDF date to a JED.
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
  integer ( kind = 4 ) y_prime
!
!  Check the date.
!
  call ymd_check_republican ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 6504 - ( 13 - m ) / 13
  m_prime = mod ( m + 12, 13 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  j1 = ( 1461 * y_prime ) / 4

  j2 = 30 * m_prime

  g = ( 3 * ( ( y_prime + 396 ) / 100 ) / 4 ) - 51
  jed = real ( j1 + j2 + d_prime - 111 - g, kind = 8 ) - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_roman ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_ROMAN converts a Roman YMDF date to a JED.
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

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
!
!  Check the date.
!
  call ymd_check_roman ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    jed = -1.0D+00
    return
  end if

  call y_roman_to_julian ( y, y2 )

  call ymdf_to_jed_julian ( y2, m, d, f, jed )

  return
end
subroutine ymdf_to_jed_saka ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_SAKA converts a Saka YMDF date to a JED.
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
  real ( kind = 8 ) jed
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
  integer ( kind = 4 ) z
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 4794 - ( 13 - m ) / 12
  m_prime = mod ( m + 10, 12 )
  d_prime = d - 1
!
!  Convert the computational date to a JED.
!
  j1 = ( 1461 * y_prime ) / 4

  z = m_prime / 6

  j2 = ( 31 - z ) * m_prime + 5 * z

  g = ( 3 * ( ( y_prime + 184 ) / 100 ) / 4 ) - 36

  jed = real ( j1 + j2 + d_prime - 1348 - g, kind = 8 ) - 0.5D+00
  jed = jed + f

  return
end
subroutine ymdf_to_jed_soor_san ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_SOOR_SAN converts a Soor San YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  call epoch_to_jed_soor_san ( jed_epoch )

  jed = jed_epoch + real ( ( d - 1 ) + 30 * ( m - 1 ) + 365 * ( y - 1 ) &
    + ( y - 1 ) / 4, kind = 8 )
  jed = jed + f

  return
end
subroutine ymdf_to_jed_syrian ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_SYRIAN converts a Syrian YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d_prime
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
!
!  Convert the calendar date to a computational date.
!
  y_prime = y + 4405 - ( 17 - m ) / 12
  m_prime = mod ( m + 6, 12 )
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
subroutine ymdf_to_jed_zoroastrian ( y, m, d, f, jed )

!*****************************************************************************80
!
!! YMDF_TO_JED_ZOROASTRIAN converts a Zoroastrian YMDF date to a JED.
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
!    Output, real ( kind = 8 ) JED, the corresponding Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  call epoch_to_jed_zoroastrian ( jed_epoch )

  jed = jed_epoch + real ( &
    ( d - 1 ) + 30 * ( m - 1 ) + 365 * ( y - 1 ), kind = 8 )
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
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, character ( len = * ) S, a representation of the date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f2
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
subroutine ymdf_to_s_english ( y, m, d, f, s )

!*****************************************************************************80
!
!! YMDF_TO_S_ENGLISH writes an English YMDF date into a string.
!
!  Format:
!
!    BC OS YYYY/MM/DD.FF
!    AD OS YYYY/MM/DD.FF
!    AD NS YYYY/MM/DD.FF
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
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, character ( len = * ) S, a representation of the date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  character ( len = 20 ) s1
  character ( len = 3 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
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
  call ymdf_check_english ( y2, m2, d2, f2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( y2 < 0 ) then
    s1 = 'BC OS '
    call i4_to_s_left (  - y2, s1(7:) )
  else if ( y2 < 1752 .or. ( y2 == 1752 .and. m2 < 9 ) .or. &
    ( y2 == 1752 .and. m2 == 9 .and. d2 < 3 ) ) then
    s1 = 'AD OS '
    call i4_to_s_left ( y2, s1(7:) )
  else
    s1 = 'AD NS '
    call i4_to_s_left ( y2, s1(7:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( m2, s2(1:2) )

  call s_cat ( s1, s2(1:2), s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( d2, s2(1:2) )

  call s_cat ( s1, s2(1:2), s1 )

  call frac_to_s ( f2, s2(1:3) )

  call s_cat ( s1, s2(1:3), s1 )

  s = s1

  return
end
subroutine ymdf_to_s_gregorian ( y, m, d, f, s )

!*****************************************************************************80
!
!! YMDF_TO_S_GREGORIAN writes a Gregorian YMDF date into a string.
!
!  Format:
!
!    AD YYYY/MM/DD.FF
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
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, character ( len = * ) S, a representation of the date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  character ( len = 25 ) s1
  character ( len = 5 ) s2
  character ( len = * ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
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
  call ymdf_check_gregorian ( y2, m2, d2, f2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( y2 < 0 ) then
    s1 = 'BC '
    call i4_to_s_left ( -y2, s1(4:) )
  else
    s1 = 'AD '
    call i4_to_s_left ( y2, s1(4:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( m2, s2(1:2) )

  call s_cat ( s1, s2(1:2), s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( d2, s2(1:3) )

  call s_cat ( s1, s2(1:3), s1 )

  call frac_to_s ( f2, s2 )

  call s_cat ( s1, s2(1:3), s1 )

  s = s1

  return
end
subroutine ymdf_to_s_hebrew ( y, m, d, f, s )

!*****************************************************************************80
!
!! YMDF_TO_S_HEBREW "prints" a Hebrew YMDF date into a string.
!
!  Format:
!
!    DayNumber.Fraction MonthName YearNumber AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, character ( len = * ) S, a representation of the date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  character ( len = 2 ) sd
  character ( len = 3 ) sf
  character ( len = 10 ) sm
  character ( len = 10 ) sy
  character ( len = * ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
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
  call ymd_check_hebrew ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  call i4_to_s_left ( y2, sy )

  call month_to_month_name_hebrew ( y2, m2, sm )

  call i4_to_s_left ( d2, sd )

  call frac_to_s ( f2, sf )

  call s_cat ( sd, sf, s )
  call s_cat1 ( s, sm, s )
  call s_cat1 ( s, sy, s )
  call s_cat1 ( s, 'AM', s )

  return
end
subroutine ymdf_to_s_islamic ( y, m, d, f, s )

!*****************************************************************************80
!
!! YMDF_TO_S_ISLAMIC writes an Islamic YMDF date into a string.
!
!  Format:
!
!    AH YYYY/MM/DD.FF
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, character ( len = * ) S, a representation of the date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
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
!
!  Check the input.
!
  call ymd_check_islamic ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  s1 = 'AH '
  call i4_to_s_left ( y2, s1(4:) )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( m2, s2 )

  call s_cat ( s1, s2, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( d2, s2 )

  call s_cat ( s1, s2, s )

  call frac_to_s ( f, s1 )

  call s_cat ( s, s1(1:3), s )

  return
end
subroutine ymdf_to_s_julian ( y, m, d, f, s )

!*****************************************************************************80
!
!! YMDF_TO_S_JULIAN writes a Julian YMDF date into a string.
!
!  Format:
!
!    AD YYYY/MM/DD.FF
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
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, character ( len = * ) S, a representation of the date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
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
!
!  Check the input.
!
  call ymd_check_julian ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( 0 <= y2 ) then
    s1 = 'AD '
    call i4_to_s_left ( y2, s1(4:) )
  else
    s1 = 'BC '
    call i4_to_s_left (  - y2, s1(4:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( m2, s2 )

  call s_cat ( s1, s2, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( d2, s2 )

  call s_cat ( s1, s2, s )

  call frac_to_s ( f, s1 )

  call s_cat ( s, s1(1:3), s )

  return
end
subroutine ymdf_to_s_numeric ( y, m, d, f, s )

!*****************************************************************************80
!
!! YMDF_TO_S_NUMERIC writes a YMDF date into a string.
!
!  Format:
!
!       YYYY/MM/DD.FF
!    or
!      -YYYY/MM/DD.FF
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
!    Output, character ( len = * ) S, a representation of the date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
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

  call frac_to_s ( f, s1 )

  call s_cat ( s, s1(1:3), s )

  return
end
subroutine ymdf_to_s_republican ( y, m, d, f, s )

!*****************************************************************************80
!
!! YMDF_TO_S_REPUBLICAN writes a Republican YMDF date into a string.
!
!  Format:
!
!    ER YYYY/MM/DD.FF
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
!    Output, character ( len = * ) S, a representation of the date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
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
!
!  Check the input.
!
  call ymd_check_republican ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  if ( 0 <= y2 ) then
    s1 = 'ER '
    call i4_to_s_left ( y2, s1(4:) )
  else
    s1 = '-ER '
    call i4_to_s_left (  - y2, s1(4:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( m2, s2 )

  call s_cat ( s1, s2, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( d2, s2 )

  call s_cat ( s1, s2, s )

  call frac_to_s ( f, s1 )

  call s_cat ( s, s1(1:3), s )

  return
end
subroutine ymdf_to_s_roman ( y, m, d, f, s )

!*****************************************************************************80
!
!! YMDF_TO_S_ROMAN writes a Roman YMDF date into a string.
!
!  Example:
!
!     Y  M   D  F    S
!    --  -  --  ---  -----------------------------------
!    56  4   1  0.1  Kalends Aprilis DVI AUC
!    56  4   2  0.2  Ante diem iv Nones Aprilis DVI AUC
!    56  4   3  0.3  Ante diem iii Nones Aprilis DVI AUC
!    56  4   4  0.4  Pridie Nones Aprilis DVI AUC
!    56  4   5  0.5  Nones Aprilis DVI AUC
!    56  4   6  0.6  Ante diem viii Ides Aprilis DVI AUC
!    56  4   7  0.7  Ante diem vii Ides Aprilis DVI AUC
!    56  4   8  0.8  Ante diem vi Ides Aprilis DVI AUC
!    56  4   9  0.9  Ante diem v Ides Aprilis DVI AUC
!    56  4  10  0.0  Ante diem iv Ides Aprilis DVI AUC
!    56  4  11  0.0  Ante diem iii Ides Aprilis DVI AUC
!    56  4  12  0.0  Pridie Ides Aprilis DVI AUC
!    56  4  13  0.0  Ides Aprilis DVI AUC
!    56  4  14  0.0  Ante diem xvii Kalends Maius DVI AUC
!    56  4  15  0.0  Ante diem xvi Kalends Maius DVI AUC
!    ...
!    56  4  28  0.0  Ante diem iv Kalends Maius DVI AUC
!    56  4  29  0.0  Ante diem iii Kalends Maius DVI AUC
!    56  4  30  0.0  Pridie Kalends Maius DVI AUC
!
!  Discussion:
!
!    "AUC" means "ab urbe condita", or "from the founding of the city".
!
!    At the moment, we ignore F.
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
!    Output, character ( len = * ) S, a string representing the date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ides
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) jday
  integer ( kind = 4 ) last
  character ( len = 10 ) lower
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ) month_length_roman
  integer ( kind = 4 ) nones
  character ( len = * ) s
  character ( len = 10 ) s_day
  character ( len = 15 ) s_month
  character ( len = 15 ) s_month_next
  character ( len = 10 ) s_year
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical   year_is_leap_roman
!
!  Copy the input.
!
  y2 = y
  m2 = m
  d2 = d
!
!  Check the input.
!
  call ymd_check_roman ( y2, m2, d2, ierror )

  if ( ierror /= 0 ) then
    s = '?'
    return
  end if

  call month_to_month_name_roman ( m2, s_month )
!
!  Get the next month's name.
!
  m3 = i4_wrap ( m2+1, 1, 12 )
  call month_to_month_name_roman ( m3, s_month_next )

  call month_to_nones_roman ( m2, nones )
  call month_to_ides_roman ( m2, ides )
  last = month_length_roman ( y2, m2 )

  if ( d2 == 1 ) then
    s = 'Kalends ' // s_month
  else if ( d2 < nones - 1 ) then
    jday = nones + 1 - d2
    call i4_to_roman ( jday, s_day )
    s_day = lower ( s_day )
    s = 'Ante diem ' // trim ( s_day ) // ' Nones ' // s_month
  else if ( d2 == nones - 1 ) then
    s = 'Pridie Nones ' // s_month
  else if ( d2 == nones ) then
    s = 'Nones ' // s_month
  else if ( d2 < ides - 1 ) then
    jday = ides + 1 - d2
    call i4_to_roman ( jday, s_day )
    s_day = lower ( s_day )
    s = 'Ante diem ' // trim ( s_day ) // ' Ides ' // s_month
  else if ( d2 == ides - 1 ) then
    s = 'Pridie Ides ' // s_month
  else if ( d2 == ides ) then
    s = 'Ides ' // s_month

  else if ( m2 == 2 .and. year_is_leap_roman ( y2 ) ) then

    if ( d2 < 25 ) then
      jday = last + 1 - d2
      call i4_to_roman ( jday, s_day )
      s_day = lower ( s_day )
      s = 'Ante diem ' // trim ( s_day ) // ' Kalends ' // s_month_next
    else if ( d2 == 25 ) then
      jday = last + 2 - d2
      call i4_to_roman ( jday, s_day )
      s_day = lower ( s_day )
      s = 'Ante diem Bis ' // trim ( s_day ) // ' Kalends ' // s_month_next
    else if ( d2 < last ) then
      jday = last + 2 - d2
      call i4_to_roman ( jday, s_day )
      s_day = lower ( s_day )
      s = 'Ante diem ' // trim ( s_day ) // ' Kalends ' // s_month_next
    else
      s = 'Pridie Kalends ' // s_month_next
    end if

  else if ( d2 < last ) then
    jday = last + 2 - d2
    call i4_to_roman ( jday, s_day )
    s_day = lower ( s_day )
    s = 'Ante diem ' // trim ( s_day ) // ' Kalends ' // s_month_next
  else
    s = 'Pridie Kalends ' // s_month_next
  end if

  call i4_to_roman ( y2, s_year )

  call s_cat1 ( s, s_year, s )
  call s_cat ( s, ' AUC', s )

  return
end
subroutine ymdf_to_week_common ( y, m, d, f, iweek )

!*****************************************************************************80
!
!! YMDF_TO_WEEK_COMMON returns the week number for a Common YMDF date.
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
!    Input, integer ( kind = 4 ) Y, M, D, is the date.
!
!    Output, integer ( kind = 4 ) IWEEK, is the week number of the date, with
!    1 for the first week, and 53 or 54 for the last.  A week begins
!    on Sunday.  The first and last weeks may be partial weeks.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) iday1
  integer ( kind = 4 ) iday2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iweek
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) ndays
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Make a local copy of the input date.
!
  d2 = d
  m2 = m
  y2 = y
  f2 = f
!
!  Check the input date.
!
  call ymdf_check_common ( y2, m2, d2, f2, ierror )

  if ( ierror /= 0 ) then
    iweek = 0
    return
  end if
!
!  Find the number of days between Y/1/1 and Y/M1/D1.
!
  d1 = 1
  m1 = 1
  y1 = y
  f1 = 0.0D+00

  call ymdf_dif_common ( y1, m1, d1, f1, y2, m2, d2, f2, days, ierror )
!
!  Find the day of the week of Y/1/1.
!
  call ymdf_to_weekday_common ( y1, m1, d1, f1, iday1 )
!
!  Find the day of the week of Y/M2/D2.
!
  call ymdf_to_weekday_common ( y2, m2, d2, f2, iday2 )
!
!  Expand the week containing Y/1/1 to begin on Sunday.
!
  ndays = int ( days ) + ( iday1 - 1 )
!
!  Expand the week containing Y/M2/D2 to end on Saturday.
!
  ndays = ndays + ( 8 - iday2 )
!
!  Now NDAYS should be an exact multiple of 7, and IWEEK is easy.
!
  iweek = ndays / 7

  return
end
subroutine ymdf_to_weekday_common ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_COMMON returns the weekday of a Common YMDF date.
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
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
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

  call ymdf_to_jed_common ( y, m, d, f, jed )

  call jed_to_weekday ( jed, w, f2 )

  return
end
subroutine ymdf_to_weekday_english ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_ENGLISH returns the weekday of an English YMDF date.
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

  call ymdf_to_jed_english ( y, m, d, f, jed )

  call jed_to_weekday ( jed, w, f2 )

  return
end
subroutine ymdf_to_weekday_english2 ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_ENGLISH2 returns the weekday of an English YMDF date.
!
!  Discussion:
!
!    Lewis Carroll's method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gary Meisters,
!    Lewis Carroll's Day-of-the-Week Algorithm,
!    Math Horizons, November 2002, pages 24-25.
!
!    Lewis Carroll (Charles Dodgson),
!    To Find the Day of the Week for Any Given Date,
!    Nature, 31 March 1887.
!
!    John Conway,
!    Tomorrow is the Day After Doomsday,
!    Eureka, 1973.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, integer ( kind = 4 ) W, is the week day number of the date, with
!    1 for Sunday, through 7 for Saturday.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) cc
  integer ( kind = 4 ) cval
  character cmp1
  character cmp2
  integer ( kind = 4 ) d
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) dval
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( 12 ) :: m_table = &
    (/ 0, 3, 3, 6, 1, 4, 6, 2, 5, 0, 3, 12 /)
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) mval
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  logical   year_is_leap_english
  integer ( kind = 4 ) yval
  integer ( kind = 4 ) yy
!
!  Dates of 2 September 1752 and earlier fall under one scheme.
!
  y1 = 1752
  m1 = 9
  d1 = 3
  f1 = 0.0D+00

  call ymdf_compare ( y, m, d, f, y1, m1, d1, f1, cmp1 )
!
!  Dates of 14 September 1752 and later fall under another.
!
  y1 = 1752
  m1 = 9
  d1 = 13
  f1 = 0.0D+00

  call ymdf_compare ( y, m, d, f, y1, m1, d1, f1, cmp2 )

  if ( cmp1 == '<' ) then

    cval = 18 - y / 100

  else if ( cmp2 == '>' ) then

    cval = 2 * ( 3 - mod ( ( y / 100 ), 4 ) )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YMDF_TO_WEEKDAY_ENGLISH2 - Fatal error!'
    write ( *, '(a)' ) '  The given date is illegal.'
    stop

  end if

  cc = y / 100
  yy = y - cc * 100
  a = yy / 12
  b = yy - 12 * a
  yval = a + b + ( b / 4 )

  mval = m_table(m)

  dval = d

  if ( ( m == 1 .or. m == 2 ) .and. year_is_leap_english ( y ) ) then
    mval = mval + 6
  end if

  w = i4_modp ( cval + yval + mval + dval, 7 ) + 1

  return
end
subroutine ymdf_to_weekday_gregorian ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_GREGORIAN returns the weekday of a Gregorian YMDF date.
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

  call ymdf_to_jed_gregorian ( y, m, d, f, jed )

  call jed_to_weekday ( jed, w, f2 )

  return
end
subroutine ymdf_to_weekday_gregorian2 ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_GREGORIAN2 returns the weekday of a Gregorian YMDF date.
!
!  Discussion:
!
!    This routine computes the day of the week from the date in
!    the Gregorian calendar.
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
!    Algorithm B,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 308.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, integer ( kind = 4 ) W, the day of the week of the date.
!    The days are numbered from Sunday through Saturday, 1 through 7.
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) t
  integer ( kind = 4 ) v
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
  integer ( kind = 4 ) z
!
!  Check the input.
!
  call ymd_check_gregorian ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    w = 0
    return
  end if

  m_prime = mod ( 9 + m, 12 )
  q = m_prime / 10
  z = ( 13 * m_prime + 2 ) / 5
  t = 28 * m_prime + z + d - 365 * q + 59

  c = i4_wrap ( t, 1, 7 )

  y_prime = y - q
  v = ( y / 4 - y_prime / 4 ) - ( y / 100 - y_prime / 100 ) &
    + ( y / 400 - y_prime / 400 )
  p = y + y / 4 - y / 100 + y / 400 - 1 - v
  n = 7 - i4_modp ( p, 7 )
  w = 1 + i4_modp ( 7 + c - n, 7 )

  return
end
subroutine ymdf_to_weekday_gregorian3 ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_GREGORIAN3 returns the weekday of a Gregorian YMDF date.
!
!  Discussion:
!
!    The algorithm is also valid for BC years.
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
!    Algorithm D,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 309.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, integer ( kind = 4 ) W, the day of the week of the date.
!    The days are numbered from Sunday through Saturday, 1 through 7.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
  integer ( kind = 4 ) y2

  call ymd_check_gregorian ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    w = 0
    return
  end if

  if ( y < 0 ) then
    y2 = 1 - y
  else
    y2 = y
  end if

  m_prime = mod ( 9 + m, 12 )
  y_prime = y2 - m_prime / 10

  do while ( y_prime < 0 )
    y_prime = y_prime + 400
  end do

  w = 1 + i4_modp ( 2 + d + ( 13 * m_prime + 2 ) / 5 &
    + y_prime + y_prime / 4 - y_prime / 100 + y_prime / 400, 7 )

  return
end
subroutine ymdf_to_weekday_gregorian4 ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_GREGORIAN4 returns the weekday of a Gregorian YMDF date.
!
!  Discussion:
!
!    This routine uses "Zeller's congruence"
!
!    W = 1 + ( [ 2.6*M-0.2 ] + K + Y + [ Y/4 ] + [ C/4 ] - 2 * C ) mod 7
!
!    where
!
!    M = month, but counting so that March = 1, April = 2, ..., Feb = 12.
!    K = the day of the month;
!    Y = the century year (that is, the last two digits of the year);
!    C = the century (the year divided by 100)
!    W = the day of the week, with Sunday = 1,
!    [X] = integer part of X.
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
!    Output, integer ( kind = 4 ) W, is the week day number of the date, with
!    1 for Sunday, through 7 for Saturday.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  real ( kind = 8 ) f
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) u
  integer ( kind = 4 ) v
  integer ( kind = 4 ) w
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
  integer ( kind = 4 ) z

  call ymd_check_gregorian ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    w = 0
    return
  end if

  if ( m <= 2 ) then
    a = m + 10
  else
    a = m - 2
  end if
!
!  What do you want to happen when Y is negative?
!
  c = i4_modp ( y, 100 )
  e = ( y - c ) / 100

  v = ( 13 * a - 1 ) / 5
  x = c / 4
  u = e / 4
  z = v + x + u + d + c - 2 * e

  w = i4_modp ( z, 7 ) + 1

  return
end
subroutine ymdf_to_weekday_gregorian5 ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_GREGORIAN5 returns the weekday of a Gregorian YMDF date.
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
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    CRC Press, 2000, page 738.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, integer ( kind = 4 ) W, is the week day number of the date, with
!    1 for Sunday, through 7 for Saturday.
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  real ( kind = 8 ) f
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) yy

  call ymd_check_gregorian ( y, m, d, ierror )

  if ( ierror /= 0 ) then
    w = 0
    return
  end if

  c = ( y / 100 )
  yy = y - c * 100

  if ( m < 3 ) then
    mm = m + 10
    yy = yy - 1
  else
    mm = m - 2
  end if

  days = d + int ( 2.6D+00 * real ( mm, kind = 8 ) - 0.2D+00 ) &
    - 2 * c + yy + ( yy / 4 ) + ( c / 4 )

  w = i4_modp ( days, 7 ) + 1

  return
end
subroutine ymdf_to_weekday_hebrew ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_HEBREW returns the weekday of a Hebrew YMDF date.
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

  call ymdf_to_jed_hebrew ( y, m, d, f, jed )

  call jed_to_weekday ( jed, w, f2 )

  return
end
subroutine ymdf_to_weekday_islamic_a ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_ISLAMIC_A returns the weekday of an Islamic A YMDF date.
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

  call ymdf_to_jed_islamic_a ( y, m, d, f, jed )

  call jed_to_weekday ( jed, w, f2 )

  return
end
subroutine ymdf_to_weekday_julian ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_JULIAN computes the weekday of a Julian YMDF date.
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
!    Output, integer ( kind = 4 ) W, the day of the week of the date.
!    The days are numbered from Sunday through Saturday, 1 through 7.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y

  call ymdf_to_jed_julian ( y, m, d, f, jed )

  call jed_to_weekday ( jed, w, f2 )

  return
end
subroutine ymdf_to_weekday_julian2 ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_JULIAN2 returns the weekday of a Julian YMDF date.
!
!  Discussion:
!
!    This routine computes the day of the week from the date in
!    the Julian calendar, that is, the calendar in force before the
!    Gregorian calendar, in which every fourth year was a leap year.
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
!    Algorithm A,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, pages 307-308.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, integer ( kind = 4 ) W, the day of the week of the date.
!    The days are numbered from Sunday through Saturday, 1 through 7.
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) t
  integer ( kind = 4 ) v
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_prime
  integer ( kind = 4 ) z

  m_prime = mod ( 9 + m, 12 )
  q = m_prime / 10
  z = ( 13 * m_prime + 2 ) / 5
  t = 28 * m_prime + z + d - 365 * q + 59

  c = i4_wrap ( t, 1, 7 )

  y_prime = y - q
  v = y / 4 - y_prime / 4
  p = y + y / 4 + 4 - v
  n = 7 - mod ( p, 7 )

  w = i4_wrap ( c + 1 - n, 1, 7 )

  return
end
subroutine ymdf_to_weekday_julian3 ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_JULIAN3 returns the week day of a Julian YMD date.
!
!  Discussion:
!
!    The algorithm is also valid for BC years.
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
!    Algorithm C,
!    Mapping Time, The Calendar and Its History,
!    Oxford, 1999, page 309.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the YMDF date.
!
!    Output, integer ( kind = 4 ) W, the day of the week of the date.
!    The days are numbered from Sunday through Saturday, 1 through 7.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_prime
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y_prime

  if ( y < 0 ) then
    y2 = 1 - y
  else
    y2 = y
  end if

  m_prime = mod ( 9 + m, 12 )
  y_prime = y2 - m_prime / 10

  do while ( y_prime < 0 )
    y_prime = y_prime + 28
  end do

  w = 1 + mod ( d + ( 13 * m_prime + 2 ) / 5 + y_prime + y_prime / 4, 7 )

  return
end
subroutine ymdf_to_weekday_republican ( y, m, d, f, w )

!*****************************************************************************80
!
!! YMDF_TO_WEEKDAY_REPUBLICAN returns the weekday of a Republican YMDF date.
!
!  Discussion:
!
!    The Republican calendar used a 10 day week.
!    There was a final "month" of 5 or 6 days.
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
!    Output, integer ( kind = 4 ) W, the day of the week of the date.
!    The days are numbered from Sunday through Saturday, 1 through 7.
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) m
  integer ( kind = 4 ) w
  integer ( kind = 4 ) y

  w = i4_wrap ( d, 1, 10 )

  return
end
subroutine ymdf_to_yjf_common ( y1, m1, d1, f1, y2, j2, f2 )

!*****************************************************************************80
!
!! YMDF_TO_YJF_COMMON converts from YMDF to YJF form in the Common calendar.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2, YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Check the date.
!
  call ymdf_check_common ( y1, m1, d1, f1, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    j2 = 0
    f2 = 0.0D+00
    return
  end if

  y2 = y1
  j2 = d1
  f2 = f1
!
!  Add in the days of the elapsed months.
!
  do m = 1, m1 - 1
    j2 = j2 + month_length_common ( y2, m )
  end do

  return
end
subroutine ymdf_to_yjf_english ( y1, m1, d1, f1, y2, j2, f2 )

!*****************************************************************************80
!
!! YMDF_TO_YJF_ENGLISH converts from YMDF to YJF form in the English calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) month_length_english
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Check the date.
!
  call ymdf_check_english ( y1, m1, d1, f1, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    j2 = 0
    f2 = 0.0D+00
    return
  end if

  y2 = y1
  j2 = d1
  f2 = f1
!
!  Add in the days of the elapsed months.
!
  do m = 1, m1 - 1
    j2 = j2 + month_length_english ( y2, m )
  end do

  return
end
subroutine ymdf_to_yjf_gregorian ( y1, m1, d1, f1, y2, j2, f2 )

!*****************************************************************************80
!
!! YMDF_TO_YJF_GREGORIAN: YMDF to YJF form in the Gregorian calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) month_length_gregorian
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Check the date.
!
  call ymd_check_gregorian ( y1, m1, d1, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    j2 = 0
    f2 = 0.0D+00
    return
  end if

  y2 = y1
  j2 = d1
  f2 = f1
!
!  Add in the days of the elapsed months.
!
  do m = 1, m1 - 1
    j2 = j2 + month_length_gregorian ( y2, m )
  end do

  return
end
subroutine ymdf_to_yjf_hebrew ( y1, m1, d1, f1, y2, j2, f2 )

!*****************************************************************************80
!
!! YMDF_TO_YJF_HEBREW converts from YMDF to YJF form in the Hebrew calendar.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) month_length_hebrew
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y2 = y1
  j2 = d1
  f2 = f1

  do m = 1, m1 - 1
    j2 = j2 + month_length_hebrew ( y1, m )
  end do

  return
end
subroutine ymdf_to_yjf_islamic ( y1, m1, d1, f1, y2, j2, f2 )

!*****************************************************************************80
!
!! YMDF_TO_YJF_ISLAMIC converts from YMDF to YJF form in the Islamic calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) month_length_islamic
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Check the date.
!
  call ymd_check_islamic ( y1, m1, d1, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    j2 = 0
    f2 = 0.0D+00
    return
  end if

  y2 = y1
  j2 = d1
  f2 = f1
!
!  Add in the days of the elapsed months.
!
  do m = 1, m1 - 1
    j2 = j2 + month_length_islamic ( y2, m )
  end do

  return
end
subroutine ymdf_to_yjf_julian ( y1, m1, d1, f1, y2, j2, f2 )

!*****************************************************************************80
!
!! YMDF_TO_YJF_JULIAN converts from YMDF to YJF form in the Julian calendar.
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
!    YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) month_length_julian
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Check the date.
!
  call ymd_check_julian ( y1, m1, d1, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    j2 = 0
    f2 = 0.0D+00
    return
  end if

  y2 = y1
  j2 = d1
  f2 = f1
!
!  Add in the days of the elapsed months.
!
  do m = 1, m1 - 1
    j2 = j2 + month_length_julian ( y2, m )
  end do

  return
end
subroutine ymdf_to_yjf_republican ( y1, m1, d1, f1, y2, j2, f2 )

!*****************************************************************************80
!
!! YMDF_TO_YJF_REPUBLICAN: YMDF to YJF form in the Republican calendar.
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
!    YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) month_length_republican
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Check the date.
!
  call ymd_check_republican ( y1, m1, d1, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    j2 = 0
    f2 = 0.0D+00
    return
  end if

  y2 = y1
  j2 = d1
  f2 = f1
!
!  Add in the days of the elapsed months.
!
  do m = 1, m1 - 1
    j2 = j2 + month_length_republican ( y2, m )
  end do

  return
end
subroutine ymdf_to_yjf_roman ( y1, m1, d1, f1, y2, j2, f2 )

!*****************************************************************************80
!
!! YMDF_TO_YJF_ROMAN converts from YMDF to YJF form in the Roman calendar.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1, the
!    YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) month_length_roman
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Check the date.
!
  call ymdf_check_roman ( y1, m1, d1, f1, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    j2 = 0
    f2 = 0.0D+00
    return
  end if

  y2 = y1
  j2 = d1
  f2 = f1
!
!  Add in the days of the elapsed months.
!
  do m = 1, m1 - 1
    j2 = j2 + month_length_roman ( y2, m )
  end do

  return
end
subroutine ymdf_to_ymdhms_common ( y1, m1, d1, f1, y2, m2, d2, h2, n2, s2 )

!*****************************************************************************80
!
!! YMDF_TO_YMDHMS_COMMON converts from YMDF to YMDHMS in the Common calendar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the YMDF date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, H2, N2, S2, the YMDHMS date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) h2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Check the date.
!
  call ymdf_check_common ( y1, m1, d1, f1, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YMDF_TO_YMDHMS_COMMON - Fatal error!'
    stop
  end if

  y2 = y1
  m2 = m1
  d2 = d1

  f2 = f1
  f2 = f2 * 24.0D+00
  h2 = int ( f2 )
  f2 = f2 - h2
  f2 = f2 * 60.0D+00
  n2 = int ( f2 )
  f2 = f2 - n2
  f2 = f2 * 60.0D+00
  s2 = nint ( f2 )

  return
end
subroutine ymdf_uniform_common ( y1, m1, d1, f1, y2, m2, d2, f2, seed, &
  y, m, d, f )

!*****************************************************************************80
!
!! YMDF_UNIFORM_COMMON picks a random Common YMDF date between two given dates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the first YMDF date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the second YMDF date.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F,
!    the random YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  call ymdf_to_jed_common ( y1, m1, d1, f1, jed1 )
  call ymdf_to_jed_common ( y2, m2, d2, f2, jed2 )

  jed = r8_uniform_ab ( jed1, jed2, seed )

  call jed_to_ymdf_common ( jed, y, m, d, f )

  return
end
subroutine ymdf_uniform_english ( y1, m1, d1, f1, y2, m2, d2, f2, seed, &
  y, m, d, f )

!*****************************************************************************80
!
!! YMDF_UNIFORM_ENGLISH: random English YMDF date between two given dates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, real ( kind = 8 ) F1,
!    the first YMDF date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2,
!    the second YMDF date.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) Y, M, D, real ( kind = 8 ) F, the
!    random YMDF date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  call ymdf_to_jed_english ( y1, m1, d1, f1, jed1 )
  call ymdf_to_jed_english ( y2, m2, d2, f2, jed2 )

  jed = r8_uniform_ab ( jed1, jed2, seed )

  call jed_to_ymdf_english ( jed, y, m, d, f )

  return
end
subroutine ymdhms_check_common ( y, m, d, h, n, s, ierror )

!*****************************************************************************80
!
!! YMDHMS_CHECK_COMMON checks a Common YMDHMS date.
!
!  Discussion:
!
!    The routine will correct certain simple errors in dates, such as
!      "11:03:42 31 September 1996"
!    which will become
!      "11:03:42 1 October 1996".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, H, N, S.
!    These items have the obvious meanings.
!    The routine may change any of these values to more reasonable values.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was detected in the
!    date, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) y
!
!  Check that the second is between 0 and 59.
!  N may get bumped up or down.
!
  call second_borrow_common ( y, m, d, h, n, s )

  call second_carry_common ( y, m, d, h, n, s )
!
!  Check that the minute is between 0 and 59.
!  H may get bumped up or down.
!
  call minute_borrow_common ( y, m, d, h, n )

  call minute_carry_common ( y, m, d, h, n )
!
!  Check that the hour is between 0 and 23.
!  D may get bumped up or down.
!
  call hour_borrow_common ( y, m, d, h )

  call hour_carry_common ( y, m, d, h )
!
!  Now make adjustments to D, M, and Y.
!
  call ymd_check_common ( y, m, d, ierror )

  return
end
subroutine ymdhms_compare ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, n2, s2, &
  cmp, ierror )

!*****************************************************************************80
!
!! YMDHMS_COMPARE compares two YMDHMS dates.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, H1, N1, S1, the first date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2,  H2, N2, S2, the second date.
!
!    Output, character CMP:
!    '<' if date 1 precedes date 2;
!    '=' if date 1 equals date 2;
!    '>' if date 1 follows date 2;
!    '?' if one of the dates was illegal.
!
!    Output, integer ( kind = 4 ) IERROR, is 1 if either date was illegal, and
!    0 otherwise.
!
  implicit none

  character cmp
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) h1
  integer ( kind = 4 ) h2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  cmp = '?'
!
!  Make local copies of the input, and then check them.
!  We need local copies because the checking routine can
!  change the input values.
!
  call ymdhms_check_common ( y1, m1, d1, h1, n1, s1, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call ymdhms_check_common ( y2, m2, d2, h2, n2, s2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
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
!  ...if necessary, compare hours in equal days...
!
        if ( h1 < h2 ) then
          cmp = '<'
        else if ( h1 > h2 ) then
          cmp = '>'
        else
!
!  ...if necessary, compare minutes in equal hours...
!
          if ( n1 < n2 ) then
            cmp = '<'
          else if ( n1 > n2 ) then
            cmp = '>'
          else
!
!  ...if necessary, compare seconds in equal minutes...
!
            if ( s1 < s2 ) then
              cmp = '<'
            else if ( s1 > s2 ) then
              cmp = '>'
            else
              cmp = '='
            end if

          end if
        end if
      end if
    end if
  end if

  return
end
subroutine ymdhms_dif_dhms ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, n2, s2, &
  nday, nhour, nminute, nsecond, ierror )

!*****************************************************************************80
!
!! YMDHMS_DIF_DHMS computes the DHMS difference between two YMDHMS dates.
!
!  Discussion:
!
!    The result is reported in days, minutes, hours and seconds.
!    The result is POSITIVE if the second date is later than the first.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, H1, N1, S1, the first date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2, H2, N2, S2, the second date.
!
!    Output, integer ( kind = 4 ) NDAY, NHOUR, NMINUTE, NSECOND,
!    the number of days, hours, minutes, seconds between the two dates.
!
!    Output, integer ( kind = 4 ) IERROR, is 1 if either date is illegal,
!    0 otherwise.
!
  implicit none

  character cmp
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) h1
  integer ( kind = 4 ) h2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) nday
  integer ( kind = 4 ) nhour
  integer ( kind = 4 ) nminute
  integer ( kind = 4 ) nsecond
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) year_length_common
!
!  Compare the dates.
!
  call ymdhms_compare ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, n2, s2, &
    cmp, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  We swap dates, if necessary, so that date 1 is never greater
!  than date 2.
!
  if ( cmp == '>' ) then
    call ymdhms_swap ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, n2, s2 )
  end if
!
!  Set beginning to 00:00:00 1 Jan Y1.
!    end       to 00:00:00 1 Jan Y1.
!
  nhour   = h2 - h1
  nminute = n2 - n1
  nsecond = s2 - s1
  nday    = 0
!
!  Subtract seconds, minutes, hours
!
  do while ( nsecond < 0 )
    nsecond = nsecond + 60
    nminute = nminute - 1
  end do

  do while ( nminute < 0 )
    nminute = nminute + 60
    nhour = nhour - 1
  end do

  do while ( nhour < 0 )
    nhour = nhour + 24
    nday = nday - 1
  end do
!
!  Set beginning to H1:M1:S1 1 Jan Y1
!    end          H2:M2:S2 1 Jan Y2.
!
  do y = y1, y2-1
    nday = nday + year_length_common ( y )
  end do
!
!  Set beginning to H1:M1:S1 1 Jan Y1
!    end          H2:M2:S2 1 M2  Y2.
!
  do m = 1, m2-1
    nday = nday + month_length_common ( y2, m )
  end do
!
!  Set beginning to H1:M1:S1 1  Jan Y1
!    end          H2:M2:S2 D2 M2  Y2.
!
  nday = nday + ( d2 - 1 )
!
!  Set beginning to H1:M1:S1 1  M1 Y1
!    end          H2:M2:S2 D2 M2 Y2.
!
  do m = 1, m1-1
    nday = nday - month_length_common ( y1, m )
  end do
!
!  Set beginning to H1:M1:S1 D1 M1 Y1
!    end          H2:M2:S2 D2 M2 Y2.
!
  nday = nday - ( d1 - 1 )
!
!  If we swapped dates, then the differences are
!  correct, but the signs should be negated.
!
!  Set beginning to H2:M2:S2 D2 M2 Y2,
!    end       to H1:M1:S1 D1 M1 Y1.
!
  if ( cmp == '>' ) then

    nday = - nday
    nhour = - nhour
    nminute = - nminute
    nsecond = - nsecond

    call ymdhms_swap ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, n2, s2 )

  end if

  return
end
subroutine ymdhms_swap ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, n2, s2 )

!*****************************************************************************80
!
!! YMDHMS_SWAP swaps the data defining two YMDHMS dates.
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
!    Input/output, integer ( kind = 4 ) Y1, M1, D1, H1, M1, S1, the first date.
!
!    Input/output, integer ( kind = 4 ) Y2, M2, D2, H2, M2, S2, the second date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) h1
  integer ( kind = 4 ) h2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  call i4_swap ( y1, y2 )

  call i4_swap ( m1, m2 )

  call i4_swap ( d1, d2 )

  call i4_swap ( h1, h2 )

  call i4_swap ( n1, n2 )

  call i4_swap ( s1, s2 )

  return
end
subroutine ymdhms_to_decimal ( y, m, d, h, n, s, yf )

!*****************************************************************************80
!
!! YMDHMS_TO_DECIMAL converts a YMDHMS date to a Decimal Y.F date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, H, N, S, the YMDHMS date.
!
!    Output, real ( kind = 8 ) YF, the Decimal date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) day_max
  integer ( kind = 4 ) days
  real ( kind = 8 ) f
  integer ( kind = 4 ) h
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_common
  real ( kind = 8 ) yf
!
!  How many days between January 1st and day D?
!
  d1 = 1
  m1 = 1
  call ymd_dif_common ( y, m1, d1, y, m, d, days, ierror )
!
!  How many days in this year total?
!
  day_max = year_length_common ( y )
!
!  The decimal part of the year is ( D + H/24 + N/24*60 + S/24*60*60 ) / DMAX.
!
  f =       real ( s, kind = 8 )      / 60.0D+00
  f = ( f + real ( n, kind = 8 )    ) / 60.0D+00
  f = ( f + real ( h, kind = 8 )    ) / 24.0D+00
  f = ( f + real ( days, kind = 8 ) ) / real ( day_max, kind = 8 )

  yf = real ( y, kind = 8 ) + f

  return
end
subroutine ymdhms_to_jed_common ( y, m, d, h, n, s, jed )

!*****************************************************************************80
!
!! YMDHMS_TO_JED_COMMON converts a Common YMDHMS date to a JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, H, N, S, the YMDHMS date.
!
!    Output, real ( kind = 8 ) JED, the Julian Ephemeris Date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  integer ( kind = 4 ) h
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1

  call ymdhms_to_ymdf_common ( y, m, d, h, n, s, y1, m1, d1, f1 )

  call ymdf_to_jed_common ( y1, m1, d1, f1, jed )

  return
end
subroutine ymdhms_to_s_common ( y, m, d, h, n, second, s )

!*****************************************************************************80
!
!! YMDHMS_TO_S_COMMON "prints" a Common YMDHMS date into a string.
!
!  Format:
!
!    CE YYYY/MM/DD HH:MM:SS
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, M, D, H, N, SECOND, the YMDHMS date.
!
!    Output, character ( len = * ) S, contains a representation of the date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) second
  character ( len = * ) s
  character ( len = 30 ) s1
  character ( len = 2 ) s2
  integer ( kind = 4 ) y

  if ( 0 <= y ) then
    s1 = 'CE '
    call i4_to_s_left ( y, s1(4:) )
  else
    s1 = 'BCE '
    call i4_to_s_left (  - y, s1(5:) )
  end if

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( m, s2 )

  call s_cat ( s1, s2, s1 )

  call s_cat ( s1, '/', s1 )

  call i4_to_s_zero ( d, s2 )

  call s_cat ( s1, s2, s1 )
!
!  Insert a pseudoblank.
!
  i = len_trim ( s1 )
  i = i + 1
  s1(i:i) = '*'

  call i4_to_s_zero ( h, s2 )

  call s_cat ( s1, s2, s1 )
  call s_cat ( s1, ':', s1 )

  call i4_to_s_zero ( n, s2 )

  call s_cat ( s1, s2, s1 )
  call s_cat ( s1, ':', s1 )

  call i4_to_s_zero ( second, s2 )

  call s_cat ( s1, s2, s1 )
!
!  Replace pseudoblank by the true thing.
!
  s1(i:i) = ' '

  s = s1

  return
end
subroutine ymdhms_to_yjf_common ( y1, m1, d1, h1, n1, s1, y2, j2, f2 )

!*****************************************************************************80
!
!! YMDHMS_TO_YJF_COMMON converts a YMDHMS date to a YJF date.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, H1, N1, S1,
!    the year, month, day, hour, minute and second of the date.
!
!    Output, integer ( kind = 4 ) Y2, J2, real ( kind = 8 ) F2, the YJF date.
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) h1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Check the date.
!
  call ymdhms_check_common ( y1, m1, d1, h1, n1, s1, ierror )

  if ( ierror /= 0 ) then
    y2 = 0
    j2 = 0
    f2 = 0.0D+00
    return
  end if
!
!  Start the day number at 0.
!
  y2 = y1
  j2 = 0
!
!  Add in the days of the elapsed months.
!
  do i = 1, m1 - 1
    j2 = j2 + month_length_common ( y2, i )
  end do
!
!  Add in the elapsed days of the current month.
!
  j2 = j2 + d1
!
!  Now compute the day fraction.
!
  f2 = real ( ( ( h1 * 60 + n1 ) * 60 + s1 ), kind = 8 ) &
    / real ( 24 * 60 * 60, kind = 8 )

  return
end
subroutine ymdhms_to_ymdf_common ( y1, m1, d1, h1, n1, s1, y2, m2, d2, f2 )

!*****************************************************************************80
!
!! YMDHMS_TO_YMDF_COMMON converts a YMDHMS date to a YMDF date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, H1, N1, S1,
!    the year, month, day, hour, minute and second of the date.
!
!    Output, integer ( kind = 4 ) Y2, M2, D2, real ( kind = 8 ) F2, 
!    the YMDF date.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) h1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
!
!  Check the date.
!
  call ymdhms_check_common ( y1, m1, d1, h1, n1, s1, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'YMDHMS_TO_YMDF_COMMON - Fatal error!'
    stop
  end if

  y2 = y1
  m2 = m1
  d2 = d1
!
!  Now compute the day fraction.
!
  f2 = real ( ( ( h1 * 60 + n1 ) * 60 + s1 ), kind = 8 ) &
    / real ( 24 * 60 * 60, kind = 8 )

  return
end
subroutine ymdhms_uniform_common ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, &
  n2, s2, seed, y, m, d, h, n, s )

!*****************************************************************************80
!
!! YMDHMS_UNIFORM_COMMON: random Common YMDHMS date between two given dates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, H1, N1, S1,
!    the first YMDHMS date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2, H2, N2, S2,
!    the second YMDHMS date.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) Y, M, D, H, N, S,
!    the random YMDHMS date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) h
  integer ( kind = 4 ) h1
  integer ( kind = 4 ) h2
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) s
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  call ymdhms_to_jed_common ( y1, m1, d1, h1, n1, s1, jed1 )

  call ymdhms_to_jed_common ( y2, m2, d2, h2, n2, s2, jed2 )

  jed = r8_uniform_ab ( jed1, jed2, seed )

  call jed_to_ymdhms_common ( jed, y, m, d, h, n, s )

  return
end
