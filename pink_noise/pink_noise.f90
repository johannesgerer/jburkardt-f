subroutine cdelay2 ( m, q )

!*****************************************************************************80
!
!! CDELAY2 is a circular buffer implementation of M-fold delay.
!
!  Example:
!
!    Suppose we call CDELAY2 12 times, always with M = 3, and with
!    Q having the input value 3 on the first call.  Q will go through
!    the following sequence of values over the 12 calls:
!
!    I   M  Qin  Qout
!
!    1   3   3   2
!    2   3   2   1
!    3   3   1   0
!    4   3   0   3
!    5   3   3   2
!    6   3   2   1
!    7   3   1   0
!    8   3   0   3
!    9   3   3   2
!   10   3   2   1
!   11   3   1   0
!   12   3   0   3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2010
!
!  Author:
!
!    Original C version by Sophocles Orfanidis.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sophocles Orfanidis,
!    Introduction to Signal Processing,
!    Prentice-Hall, 1995,
!    ISBN: 0-13-209172-0,
!    LC: TK5102.5.O246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the maximum value that Q can have.
!
!    Input/output, integer ( kind = 4 ) Q, a counter which is decremented 
!    on every call.
!    However, the value "after" 0 is M.  
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
!
!  Decrement the offset.
!
  q = q - 1
!
!  Q = - 1 wraps to Q = M.
!
  call wrap2 ( m, q )

  return
end
subroutine corr ( n, x, m, r )

!*****************************************************************************80
!
!! CORR computes the sample correlation of a signal sample.
!
!  Discussion:
!
!    The sample correlation is defined, for 0 <= i < N, as
!
!      R(i) = 1/N * sum ( 0 <= j <= N - 1 - i ) X(i+j) * X(j)
!
!    The sample correlation is an estimate of the correlation function.
!
!    It is usually the case that the signal X is assumed to
!    have zero mean.  Here, we compute the mean and adjust the
!    calculation accordingly:
!
!      R(i) = 1/N * sum ( 0 <= j <= N - 1 - i ) 
!        ( X(i+j) - Xbar ) * ( X(j) - Xbar )
!
!    Experience suggests that only the first 5 or 10 percent of
!    the lags are statistically reliable, so that one might choose
!    M = N / 20 or M = N / 10, for instance.
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
!    Sophocles Orfanidis,
!    Introduction to Signal Processing,
!    Prentice-Hall, 1995,
!    ISBN: 0-13-209172-0,
!    LC: TK5102.5.O246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equally spaced signal
!    samples.
!
!    Input, real ( kind = 8 ) X(0:N-1), the signal samples.
!
!    Input, integer ( kind = 4 ) M, the maximum lag to consider.
!    0 <= M < N.
!
!    Output, real ( kind = 8 ) R(0:M), the sample correlations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(0:m)
  real ( kind = 8 ) x(0:n-1)
  real ( kind = 8 ) xbar

  r(0:m) = 0.0D+00

  xbar = sum ( x(0:n-1) ) / real ( n, kind = 8 )

  do i = 0, m
    do j = 0, n - i - 1
      r(i) = r(i) + ( x(i+j) - xbar ) * ( x(j) - xbar )
    end do
  end do

  r(0:m) = r(0:m) / real ( n, kind = 8 )

  return
end
subroutine cross_corr ( n, x, y, m, r )

!*****************************************************************************80
!
!! CROSS_CORR computes the sample cross correlation between two signal samples.
!
!  Discussion:
!
!    The sample cross correlation is defined, for 0 <= i < N, as
!
!      R(i) = 1/N * sum ( 0 <= j <= N - 1 - i ) X(i+j) * Y(j)
!
!    The sample cross correlation is an estimate of the cross 
!    correlation function.
!
!    It is usually the case that the signals X and Y are assumed to
!    have zero mean.  Here, we compute the means and adjust the
!    calculation accordingly:
!
!      R(i) = 1/N * sum ( 0 <= j <= N - 1 - i ) 
!        ( X(i+j) - Xbar ) * ( Y(j) - Ybar )
!
!    Experience suggests that only the first 5 or 10 percent of
!    the lags are statistically reliable, so that one might choose
!    M = N / 20 or M = N / 10, for instance.
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
!    Sophocles Orfanidis,
!    Introduction to Signal Processing,
!    Prentice-Hall, 1995,
!    ISBN: 0-13-209172-0,
!    LC: TK5102.5.O246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equally spaced signal
!    samples.
!
!    Input, real ( kind = 8 ) X(0:N-1), Y(0:N-1), the signal samples.
!
!    Input, integer ( kind = 4 ) M, the maximum lag to consider.
!    0 <= M < N.
!
!    Output, real ( kind = 8 ) R(0:M), the sample correlations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(0:m)
  real ( kind = 8 ) x(0:n-1)
  real ( kind = 8 ) xbar
  real ( kind = 8 ) y(0:n-1)
  real ( kind = 8 ) ybar

  r(0:m) = 0.0D+00

  xbar = sum ( x(0:n-1) ) / real ( n, kind = 8 )
  ybar = sum ( y(0:n-1) ) / real ( n, kind = 8 )

  do i = 0, m
    do j = 0, n - i - 1
      r(i) = r(i) + ( x(i+j) - xbar ) * ( y(j) - ybar )
    end do
  end do

  r(0:m) = r(0:m) / real ( n, kind = 8 )

  return
end
function ran1f ( b, u, q )

!*****************************************************************************80
!
!! RAN1F is a 1/F random number generator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2010
!
!  Author:
!
!    Original C version by Sophocles Orfanidis.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sophocles Orfanidis,
!    Introduction to Signal Processing,
!    Prentice-Hall, 1995,
!    ISBN: 0-13-209172-0,
!    LC: TK5102.5.O246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) B, the number of signals to combine.
!    For this algorithm, B cannot be more than 31!
!
!    Input/output, real ( kind = 8 ) U(B), the signals to combine.  It is 
!    expected that each of the initial values of U will be drawn from a 
!    distribution with zero mean.
!
!    Input/output, integer ( kind = 4 ) Q(B), a set of counters that determine 
!    when each entry of U is to be updated.
!
!    Output, real ( kind = 8 ) RAN1F, the value.
!
  implicit none

  integer ( kind = 4 ) b

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) q(b)
  real ( kind = 8 ) ran1f
  real ( kind = 8 ) ranh
  real ( kind = 8 ) u(b)
  real ( kind = 8 ) y

  if ( 31 < b ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAN1F - Fatal error!'
    write ( *, '(a)' ) '  32 <= B, too many signals.'
    stop
  end if

  y = 0.0D+00

  j = 1
  do i = 1, b
    y = y + ranh ( j, u(i), q(i) )
    j = j * 2
  end do

  if ( 0 < b ) then
    y = y / real ( b, kind = 8 )
  end if

  ran1f = y

  return
end
function ranh ( d, u, q )

!*****************************************************************************80
!
!  Purpose:
!
!    RANH is a hold random number generator of period D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2010
!
!  Author:
!
!    Original C version by Sophocles Orfanidis.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sophocles Orfanidis,
!    Introduction to Signal Processing,
!    Prentice-Hall, 1995,
!    ISBN: 0-13-209172-0,
!    LC: TK5102.5.O246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) D, the hold period.  D must be at least 1.
!
!    Input/output, real ( kind = 8 ) U, a value to be held until Q has 
!    decremented to 0, when Q will be reset to D, and U will be randomly reset.
!
!    Input/output, integer ( kind = 4 ) Q, a counter which is decremented by 1 
!    on each call until reaching 0.
!
!    Output, real ( kind = 8 ) RANH, the input value of U.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) q
  real ( kind = 8 ) ranh
  real ( kind = 8 ) u
  real ( kind = 8 ) y

  if ( d < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANH - Fatal error!'
    write ( *, '(a)' ) '  D < 1.'
    stop
  end if
!
!  Hold this sample for D calls.
!
  y = u
!
!  Decrement Q and wrap mod D.
!
  call cdelay2 ( d - 1, q )
!
!  Every D calls, get a new U with zero mean.
!
  if ( q == 0 ) then
    call random_number ( harvest = u )
    u = 2.0D+00 * u - 1.0D+00
  end if

  ranh = y

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
subroutine wrap2 ( m, q )

!*****************************************************************************80
!
!! WRAP2 is a circular wrap of the pointer offset Q.
!
!  Discussion:
!
!    Input values of Q between 0 and M are "legal".
!    Values of Q below 0 are incremented by M + 1 until they are legal.
!    Values of Q above M are decremented by M + 1 until they become legal.
!    The legal value is the output value of the function.
!
!  Example:
!
!    M  Qin  Qout
!
!    3  -5   3
!    3  -4   0
!    3  -3   1
!    3  -2   2
!    3  -1   3
!    3   0   0
!    3   1   1
!    3   2   2
!    3   3   3
!    3   4   0
!    3   5   1
!    3   6   2
!    3   7   3
!    3   8   0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2010
!
!  Author:
!
!    Original C version by Sophocles Orfanidis.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sophocles Orfanidis,
!    Introduction to Signal Processing,
!    Prentice-Hall, 1995,
!    ISBN: 0-13-209172-0,
!    LC: TK5102.5.O246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the maximum acceptable value for outputs.
!    M must be at least 0.
!
!    Input/output, integer ( kind = 4 ) Q, the value to be wrapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) q

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRAP2 - Fatal error!'
    write ( *, '(a)' ) '  M < 0.'
    stop
  end if
!
!  When Q = M + 1, it wraps to Q = 0.
!
  do while ( m < q )
    q = q - m - 1
  end do
!
!  When Q = - 1, it wraps to Q = M.
!
  do while ( q < 0 )
    q = q + m + 1
  end do

  return
end
