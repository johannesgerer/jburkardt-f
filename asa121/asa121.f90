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
function trigamma ( x, ifault )

!*****************************************************************************80
!
!! TRIGAMMA calculates trigamma(x) = d**2 log(gamma(x)) / dx**2
!
!  Modified:
!
!    19 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by BE Schneider.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    BE Schneider,
!    Algorithm AS 121:
!    Trigamma Function,
!    Applied Statistics, 
!    Volume 27, Number 1, pages 97-99, 1978.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the trigamma function.
!    0 < X.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error.
!    1, X <= 0.
!
!    Output, real ( kind = 8 ) TRIGAMMA, the value of the trigamma function.
!
  implicit none

  real ( kind = 8 ), parameter :: a = 0.0001D+00
  real ( kind = 8 ), parameter :: b = 5.0D+00
  real ( kind = 8 ), parameter :: b2 =  0.1666666667D+00
  real ( kind = 8 ), parameter :: b4 = -0.03333333333D+00
  real ( kind = 8 ), parameter :: b6 =  0.02380952381D+00
  real ( kind = 8 ), parameter :: b8 = -0.03333333333D+00
  integer ( kind = 4 ) ifault
  real ( kind = 8 ) trigamma
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!
!  Check the input.
!
  if ( x <= 0.0D+00 ) then
    ifault = 1
    trigamma = 0.0D+00
    return
  end if

  ifault = 0
  z = x
!
!  Use small value approximation if X <= A.
!
  if ( x <= a ) then
    trigamma = 1.0D+00 / x / x
    return
  end if
!
!  Increase argument to ( X + I ) >= B.
!
  trigamma = 0.0D+00

  do while ( z < b )
    trigamma = trigamma + 1.0D+00 / z / z
    z = z + 1.0D+00
  end do
!
!  Apply asymptotic formula if argument is B or greater.
!
  y = 1.0D+00 / z / z

  trigamma = trigamma + 0.5D+00 * &
      y + ( 1.0D+00 &
    + y * ( b2  &
    + y * ( b4  &
    + y * ( b6  &
    + y *   b8 )))) / z

  return
end
subroutine trigamma_values ( n_data, x, fx )

!*****************************************************************************80
!
!! TRIGAMMA_VALUES returns some values of the TriGamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      PolyGamma[1,x]
!
!    TriGamma(X) = d^2 ln ( Gamma ( X ) ) / d X^2
!
!  Modified:
!
!    16 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1644934066848226D+01, &
     0.1433299150792759D+01, &
     0.1267377205423779D+01, &
     0.1134253434996619D+01, &
     0.1025356590529597D+01, &
     0.9348022005446793D+00, &
     0.8584318931245799D+00, &
     0.7932328301639984D+00, &
     0.7369741375017002D+00, &
     0.6879720582426356D+00, &
     0.6449340668482264D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
