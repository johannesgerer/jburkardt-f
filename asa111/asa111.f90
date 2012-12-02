subroutine normal_01_cdf_values ( n_data, x, fx )

!*****************************************************************************80
!
!! NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NormalDistribution [ 0, 1 ]
!      CDF [ dist, x ]
!
!  Modified:
!
!    28 August 2004
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

  integer ( kind = 4 ), parameter :: n_max = 17

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.5398278372770290D+00, &
    0.5792597094391030D+00, &
    0.6179114221889526D+00, &
    0.6554217416103242D+00, &
    0.6914624612740131D+00, &
    0.7257468822499270D+00, &
    0.7580363477769270D+00, &
    0.7881446014166033D+00, &
    0.8159398746532405D+00, &
    0.8413447460685429D+00, &
    0.9331927987311419D+00, &
    0.9772498680518208D+00, &
    0.9937903346742239D+00, &
    0.9986501019683699D+00, &
    0.9997673709209645D+00, &
    0.9999683287581669D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000D+00, &
    0.1000000000000000D+00, &
    0.2000000000000000D+00, &
    0.3000000000000000D+00, &
    0.4000000000000000D+00, &
    0.5000000000000000D+00, &
    0.6000000000000000D+00, &
    0.7000000000000000D+00, &
    0.8000000000000000D+00, &
    0.9000000000000000D+00, &
    0.1000000000000000D+01, &
    0.1500000000000000D+01, &
    0.2000000000000000D+01, &
    0.2500000000000000D+01, &
    0.3000000000000000D+01, &
    0.3500000000000000D+01, &
    0.4000000000000000D+01 /)

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
function ppnd ( p, ifault )

!*****************************************************************************80
!
!! PPND produces the normal deviate value corresponding to lower tail area = P.
!
!  Modified:
!
!    21 January 2008
!
!  Author:
!
!    J Beasley, S Springer
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    J Beasley, S Springer,
!    Algorithm AS 111:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 26, Number 1, 1977, pages 118-121.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the value of the cumulative probability
!    densitity function.  0 < P < 1.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error.
!    1, P <= 0 or P >= 1.  PPND is returned as 0.
!
!    Output, real ( kind = 8 ) PPND, the normal deviate value with the property
!    that the probability of a standard normal deviate being less than or
!    equal to PPND is P.
!
  implicit none

  real ( kind = 8 ), parameter :: a0 = 2.50662823884D+00
  real ( kind = 8 ), parameter :: a1 = -18.61500062529D+00
  real ( kind = 8 ), parameter :: a2 = 41.39119773534D+00
  real ( kind = 8 ), parameter :: a3 = -25.44106049637D+00
  real ( kind = 8 ), parameter :: b1 = -8.47351093090D+00
  real ( kind = 8 ), parameter :: b2 = 23.08336743743D+00
  real ( kind = 8 ), parameter :: b3 = -21.06224101826D+00
  real ( kind = 8 ), parameter :: b4 = 3.13082909833D+00
  real ( kind = 8 ), parameter :: c0 = -2.78718931138D+00
  real ( kind = 8 ), parameter :: c1 = -2.29796479134D+00
  real ( kind = 8 ), parameter :: c2 = 4.85014127135D+00
  real ( kind = 8 ), parameter :: c3 = 2.32121276858D+00
  real ( kind = 8 ), parameter :: d1 = 3.54388924762D+00
  real ( kind = 8 ), parameter :: d2 = 1.63706781897D+00
  integer ( kind = 4 ) ifault
  real ( kind = 8 ) p
  real ( kind = 8 ) ppnd
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: split = 0.42D+00
  real ( kind = 8 ) value

  ifault = 0
!
!  0.08 < P < 0.92
!
  if ( abs ( p - 0.5D+00 ) <= split ) then

    r = ( p - 0.5D+00 ) * ( p - 0.5D+00 )

    value = ( p - 0.5D+00 ) * ( ( ( &
        a3   * r &
      + a2 ) * r &
      + a1 ) * r &
      + a0 ) / ( ( ( ( &
        b4   * r &
      + b3 ) * r &
      + b2 ) * r &
      + b1 ) * r &
      + 1.0D+00 )
!
!  P < 0.08 or P > 0.92,
!  R = min ( P, 1-P )
!
  else if ( 0.0D+00 < p .and. p < 1.0D+00 ) then

    if ( 0.5D+00 < p ) then
      r = sqrt ( - log ( 1.0D+00 - p ) )
    else
      r = sqrt ( - log ( p ) )
    end if

    value = ( ( ( &
        c3   * r &
      + c2 ) * r &
      + c1 ) * r &
      + c0 ) / ( ( &
        d2   * r &
      + d1 ) * r &
      + 1.0D+00 )

    if ( p < 0.5D+00 ) then
      value = - value
    end if
!
!  P <= 0.0 or 1.0 <= P
!
  else

    ifault = 1
    value = 0.0D+00

  end if

  ppnd = value

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
