function alogam ( x, ifault )

!*****************************************************************************80
!
!! ALOGAM computes the logarithm of the Gamma function.
!
!  Modified:
!
!    28 March 1999
!
!  Author:
!
!    Malcolm Pike, David Hill
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Malcolm Pike, David Hill,
!    Algorithm 291:
!    Logarithm of Gamma Function,
!    Communications of the ACM,
!    Volume 9, Number 9, September 1966, page 684.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function.
!    X should be greater than 0.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error.
!    1, X <= 0.
!
!    Output, real ( kind = 8 ) ALOGAM, the logarithm of the Gamma
!    function of X.
!
  implicit none

  real    ( kind = 8 ) alogam
  real    ( kind = 8 ) f
  integer ( kind = 4 ) ifault
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) z

  if ( x <= 0.0D+00 ) then
    ifault = 1
    alogam = 0.0D+00
    return
  end if

  ifault = 0
  y = x

  if ( x < 7.0D+00 ) then

    f = 1.0D+00
    z = y

    do while ( z < 7.0D+00 )
      f = f * z
      z = z + 1.0D+00
    end do

    y = z
    f = - log ( f )

  else

    f = 0.0D+00

  end if

  z = 1.0D+00 / y / y

  alogam = f + ( y - 0.5D+00 ) * log ( y ) - y &
    + 0.918938533204673D+00 + &
    ((( &
    - 0.000595238095238D+00   * z &
    + 0.000793650793651D+00 ) * z &
    - 0.002777777777778D+00 ) * z &
    + 0.083333333333333D+00 ) / y

  return
end
subroutine beta_cdf_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! BETA_CDF_VALUES returns some values of the Beta CDF.
!
!  Discussion:
!
!    The incomplete Beta function may be written
!
!      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
!                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
!
!    Thus,
!
!      BETA_INC(A,B,0.0) = 0.0
!      BETA_INC(A,B,1.0) = 1.0
!
!    The incomplete Beta function is also sometimes called the
!    "modified" Beta function, or the "normalized" Beta function
!    or the Beta CDF (cumulative density function.
!
!    In Mathematica, the function can be evaluated by:
!
!      BETA[X,A,B] / BETA[A,B]
!
!    The function can also be evaluated by using the Statistics package:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = BetaDistribution [ a, b ]
!      CDF [ dist, x ]
!
!  Modified:
!
!    04 January 2006
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
!    Karl Pearson,
!    Tables of the Incomplete Beta Function,
!    Cambridge University Press, 1968,
!    LC: QA351.P38.
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
!    Output, real ( kind = 8 ) A, B, the parameters of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 42

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.5D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    30.0D+00, &
    30.0D+00, &
    40.0D+00, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.2D+01, &
     0.3D+01, &
     0.4D+01, &
     0.5D+01 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.0D+00, &
     0.5D+00, &
     5.0D+00, &
     5.0D+00, &
    10.0D+00, &
     5.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.2D+01, &
     0.3D+01, &
     0.4D+01, &
     0.5D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6376856085851985D-01, &
    0.2048327646991335D+00, &
    0.1000000000000000D+01, &
    0.0000000000000000D+00, &
    0.5012562893380045D-02, &
    0.5131670194948620D-01, &
    0.2928932188134525D+00, &
    0.5000000000000000D+00, &
    0.2800000000000000D-01, &
    0.1040000000000000D+00, &
    0.2160000000000000D+00, &
    0.3520000000000000D+00, &
    0.5000000000000000D+00, &
    0.6480000000000000D+00, &
    0.7840000000000000D+00, &
    0.8960000000000000D+00, &
    0.9720000000000000D+00, &
    0.4361908850559777D+00, &
    0.1516409096347099D+00, &
    0.8978271484375000D-01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.4598773297575791D+00, &
    0.2146816102371739D+00, &
    0.9507364826957875D+00, &
    0.5000000000000000D+00, &
    0.8979413687105918D+00, &
    0.2241297491808366D+00, &
    0.7586405487192086D+00, &
    0.7001783247477069D+00, &
    0.5131670194948620D-01, &
    0.1055728090000841D+00, &
    0.1633399734659245D+00, &
    0.2254033307585166D+00, &
    0.3600000000000000D+00, &
    0.4880000000000000D+00, &
    0.5904000000000000D+00, &
    0.6723200000000000D+00, &
    0.2160000000000000D+00, &
    0.8370000000000000D-01, &
    0.3078000000000000D-01, &
    0.1093500000000000D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.01D+00, &
    0.10D+00, &
    1.00D+00, &
    0.00D+00, &
    0.01D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    0.50D+00, &
    0.90D+00, &
    0.50D+00, &
    1.00D+00, &
    0.50D+00, &
    0.80D+00, &
    0.60D+00, &
    0.80D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.70D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_log_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_LOG_VALUES returns some values of the Log Gamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Log[Gamma[x]]
!
!  Modified:
!
!    14 August 2004
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

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1524063822430784D+01, &
     0.7966778177017837D+00, &
     0.3982338580692348D+00, &
     0.1520596783998375D+00, &
     0.0000000000000000D+00, &
    -0.4987244125983972D-01, &
    -0.8537409000331584D-01, &
    -0.1081748095078604D+00, &
    -0.1196129141723712D+00, &
    -0.1207822376352452D+00, &
    -0.1125917656967557D+00, &
    -0.9580769740706586D-01, &
    -0.7108387291437216D-01, &
    -0.3898427592308333D-01, &
    0.00000000000000000D+00, &
    0.69314718055994530D+00, &
    0.17917594692280550D+01, &
    0.12801827480081469D+02, &
    0.39339884187199494D+02, &
    0.71257038967168009D+02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.20D+00, &
     0.40D+00, &
     0.60D+00, &
     0.80D+00, &
     1.00D+00, &
     1.10D+00, &
     1.20D+00, &
     1.30D+00, &
     1.40D+00, &
     1.50D+00, &
     1.60D+00, &
     1.70D+00, &
     1.80D+00, &
     1.90D+00, &
     2.00D+00, &
     3.00D+00, &
     4.00D+00, &
    10.00D+00, &
    20.00D+00, &
    30.00D+00 /)

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
subroutine mdbeta ( x, p, q, prob, ier )

!*****************************************************************************80
!
!! MDBETA evaluates the incomplete beta function.
!
!  Modified:
!
!    30 January 2008
!
!  Author:
!
!    Oliver Ludwig
!    Modifications by John Burkardt
!
!  Reference:
!
!    Oliver Ludwig,
!    Algorithm 179:
!    Incomplete Beta Ratio,
!    Communications of the ACM,
!    Volume 6, Number 6, June 1963, page 314.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value to which function is to be
!    integrated.  X must be in the range [0,1] inclusive.
!
!    Input, real ( kind = 8 ) P, the first parameter.  P must be greater
!    than 0.0.
!
!    Input, real ( kind = 8 ) Q, the second parameter.  Q must be greater
!    than 0.0.
!
!    Output, real ( kind = 8 ) PROB.  The probability that a random variable
!    from a Beta distribution having parameters P and Q will be less than
!    or equal to X.
!
!    Output, integer IER, error parameter.
!    0, normal exit.
!    1, X is not in the range [0,1] inclusive.
!    2, P or Q is less than or equal to 0.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) ALEPS, the logarithm of EPS1.
!
!    Local, real ( kind = 8 ) EPS, the machine precision.
!
!    Local, real ( kind = 8 ) EPS1, the smallest representable number.
!
  implicit none

  real    ( kind = 8 ), parameter :: aleps = - 179.6016D+00
  real    ( kind = 8 ) alogam
  real    ( kind = 8 ) c
  real    ( kind = 8 ) cnt
  real    ( kind = 8 ) d4
  real    ( kind = 8 ) dp
  real    ( kind = 8 ) dq
  real    ( kind = 8 ), parameter :: eps = 2.2D-16
  real    ( kind = 8 ), parameter :: eps1 = 1.0D-78
  real    ( kind = 8 ) finsum
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ifault
  real    ( kind = 8 ) infsum
  integer ( kind = 4 ) interval
  real    ( kind = 8 ) p
  real    ( kind = 8 ) p1
  real    ( kind = 8 ) pq
  real    ( kind = 8 ) prob
  real    ( kind = 8 ) ps
  real    ( kind = 8 ) px
  real    ( kind = 8 ) q
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) wh
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xb
  real    ( kind = 8 ) y
!
!  Check ranges of the arguments.
!
  prob = 0.0D+00
  y = x

  if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
    ier = 1
    return
  end if

  if ( p <= 0.0D+00 .or. q <= 0.0D+00 ) then
    ier = 2
    return
  end if

  ier = 0

  if ( x <= 0.5D+00 ) then
    interval = 0
  else
    interval = 1
    temp = p
    p = q
    q = temp
    y = 1.0D+00 - y
  end if

  if ( x == 0.0D+00 .or. x == 1.0D+00 ) then

    prob = 0.0D+00

    if ( interval /= 0 ) then
      prob = 1.0D+00 - prob
      temp = p
      p = q
      q = temp
    end if

    return
  end if

  ib = q
  temp = ib
  ps = q - real ( ib, kind = 8 )

  if ( q == temp ) then
    ps = 1.0D+00
  end if

  dp = p
  dq = q
  px = dp * log ( y )
  pq = alogam ( dp + dq, ifault )
  p1 = alogam ( dp, ifault )
  c = alogam ( dq, ifault )
  d4 = log ( dp )
  xb = px + alogam ( ps + dp, ifault ) - alogam ( ps, ifault ) - d4 - p1
!
!  Scaling
!
  ib = int ( xb / aleps )
  infsum = 0.0D+00
!
!  First term of a decreasing series will underflow.
!
  if ( ib == 0 ) then

    infsum = exp ( xb )
    cnt = infsum * dp
!
!  CNT will equal exp ( temp ) * ( 1.d0 - ps ) * i * p * y**i / factorial ( i ).
!
    wh = 0.0D+00

    do

      wh = wh + 1.0D+00
      cnt = cnt * ( wh - ps ) * y / wh
      xb = cnt / ( dp + wh )
      infsum = infsum + xb

      if ( xb / eps < infsum ) then
        exit
      end if

    end do

  end if

  finsum = 0.0D+00

  if ( dq <= 1.0D+00 ) then

    prob = finsum + infsum

    if ( interval /= 0 ) then
      prob = 1.0D+00 - prob
      temp = p
      p = q
      q = temp
    end if

    return
  end if

  xb = px + dq * log ( 1.0D+00 - y ) + pq - p1 - log ( dq ) - c
!
!  Scaling.
!
  ib = int ( xb / aleps )

  if ( ib < 0 ) then
    ib = 0
  end if

  c = 1.0D+00 / ( 1.0D+00 - y )
  cnt = exp ( xb - real ( ib, kind = 8 ) * aleps )
  ps = dq
  wh = dq

  do

    wh = wh - 1.0D+00

    if ( wh <= 0.0D+00 ) then

      prob = finsum + infsum

      if ( interval /= 0 ) then
        prob = 1.0D+00 - prob
        temp = p
        p = q
        q = temp
      end if

      exit

    end if

    px = ( ps * c ) / ( dp + wh )

    if ( px <= 1.0D+00 ) then

      if ( cnt / eps <= finsum .or. cnt <= eps1 / px ) then

        prob = finsum + infsum

        if ( interval /= 0 ) then
          prob = 1.0D+00 - prob
          temp = p
          p = q
          q = temp
        end if

        exit

      end if

    end if

    cnt = cnt * px
!
!  Rescale.
!
    if ( 1.0D+00 < cnt ) then
      ib = ib - 1
      cnt = cnt * eps1
    end if

    ps = wh

    if ( ib == 0 ) then
      finsum = finsum + cnt
    end if

  end do

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
