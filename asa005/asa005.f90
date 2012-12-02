function alngam ( xvalue, ifault )

!*****************************************************************************80
!
!! ALNGAM computes the logarithm of the gamma function.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Allan Macleod.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Allan Macleod,
!    Algorithm AS 245,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the Gamma function.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!    1, XVALUE is less than or equal to 0.
!    2, XVALUE is too big.
!
!    Output, real ( kind = 8 ) ALNGAM, the logarithm of the gamma function of X.
!
  implicit none

  real ( kind = 8 ) alngam
  real ( kind = 8 ), parameter :: alr2pi = 0.918938533204673D+00
  integer ( kind = 4 ) ifault
  real ( kind = 8 ), dimension ( 9 ) :: r1 = (/ &
    -2.66685511495D+00, &
    -24.4387534237D+00, &
    -21.9698958928D+00, &
     11.1667541262D+00, &
     3.13060547623D+00, &
     0.607771387771D+00, &
     11.9400905721D+00, &
     31.4690115749D+00, &
     15.2346874070D+00 /)
  real ( kind = 8 ), dimension ( 9 ) :: r2 = (/ &
    -78.3359299449D+00, &
    -142.046296688D+00, &
     137.519416416D+00, &
     78.6994924154D+00, &
     4.16438922228D+00, &
     47.0668766060D+00, &
     313.399215894D+00, &
     263.505074721D+00, &
     43.3400022514D+00 /)
  real ( kind = 8 ), dimension ( 9 ) :: r3 = (/ &
    -2.12159572323D+05, &
     2.30661510616D+05, &
     2.74647644705D+04, &
    -4.02621119975D+04, &
    -2.29660729780D+03, &
    -1.16328495004D+05, &
    -1.46025937511D+05, &
    -2.42357409629D+04, &
    -5.70691009324D+02 /)
  real ( kind = 8 ), dimension ( 5 ) :: r4 = (/ &
     0.279195317918525D+00, &
     0.4917317610505968D+00, &
     0.0692910599291889D+00, &
     3.350343815022304D+00, &
     6.012459259764103D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ), parameter :: xlge = 5.10D+05
  real ( kind = 8 ), parameter :: xlgst = 1.0D+30
  real ( kind = 8 ) xvalue
  real ( kind = 8 ) y

  x = xvalue
  alngam = 0.0D+00
!
!  Check the input.
!
  if ( xlgst <= x ) then
    ifault = 2
    return
  end if

  if ( x <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  ifault = 0
!
!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
!
  if ( x < 1.5D+00 ) then

    if ( x < 0.5D+00 ) then

      alngam = - log ( x )
      y = x + 1.0D+00
!
!  Test whether X < machine epsilon.
!
      if ( y == 1.0D+00 ) then
        return
      end if

    else

      alngam = 0.0D+00
      y = x
      x = ( x - 0.5D+00 ) - 0.5D+00

    end if

    alngam = alngam + x * (((( &
        r1(5)   * y &
      + r1(4) ) * y &
      + r1(3) ) * y &
      + r1(2) ) * y &
      + r1(1) ) / (((( &
                  y &
      + r1(9) ) * y &
      + r1(8) ) * y &
      + r1(7) ) * y &
      + r1(6) )

    return

  end if
!
!  Calculation for 1.5 <= X < 4.0.
!
  if ( x < 4.0D+00 ) then

    y = ( x - 1.0D+00 ) - 1.0D+00

    alngam = y * (((( &
        r2(5)   * x &
      + r2(4) ) * x &
      + r2(3) ) * x &
      + r2(2) ) * x &
      + r2(1) ) / (((( &
                  x &
      + r2(9) ) * x &
      + r2(8) ) * x &
      + r2(7) ) * x &
      + r2(6) )
!
!  Calculation for 4.0 <= X < 12.0.
!
  else if ( x < 12.0D+00 ) then

    alngam = (((( &
        r3(5)   * x &
      + r3(4) ) * x &
      + r3(3) ) * x &
      + r3(2) ) * x &
      + r3(1) ) / (((( &
                  x &
      + r3(9) ) * x &
      + r3(8) ) * x &
      + r3(7) ) * x &
      + r3(6) )
!
!  Calculation for 12.0 <= X.
!
  else

    y = log ( x )
    alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

    if ( x <= xlge ) then

      x1 = 1.0D+00 / x
      x2 = x1 * x1

      alngam = alngam + x1 * ( ( &
             r4(3)   * &
        x2 + r4(2) ) * &
        x2 + r4(1) ) / ( ( &
        x2 + r4(5) ) * &
        x2 + r4(4) )

    end if

  end if

  return
end
function alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    David Hill
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a2 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a3 = 5.92885724438D+00
  real ( kind = 8 ) alnorm
  real ( kind = 8 ), parameter :: b1 = -29.8213557807D+00
  real ( kind = 8 ), parameter :: b2 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: c1 = -0.000000038052D+00
  real ( kind = 8 ), parameter :: c2 = 0.000398064794D+00
  real ( kind = 8 ), parameter :: c3 = -0.151679116635D+00
  real ( kind = 8 ), parameter :: c4 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: c5 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: c6 = 3.99019417011D+00
  real ( kind = 8 ), parameter :: con = 1.28D+00
  real ( kind = 8 ), parameter :: d1 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: d2 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: d3 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: d4 = -15.1508972451D+00
  real ( kind = 8 ), parameter :: d5 = 30.789933034D+00
  real ( kind = 8 ), parameter :: ltone = 7.0D+00
  real ( kind = 8 ), parameter :: p = 0.398942280444D+00
  real ( kind = 8 ), parameter :: q = 0.39990348504D+00
  real ( kind = 8 ), parameter :: r = 0.398942280385D+00
  logical up
  logical upper
  real ( kind = 8 ), parameter :: utzero = 18.66D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  up = upper
  z = x

  if ( z < 0.0D+00 ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0D+00
    else
      alnorm = 1.0D+00
    end if

    return

  end if

  y = 0.5D+00 * z * z

  if ( z <= con ) then

    alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 &
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0D+00 - alnorm
  end if

  return
end
function prncst ( st, idf, d, ifault )

!*****************************************************************************80
!
!! PRNCST computes the lower tail of noncentral T distribution.
!
!  Modified:
!
!    26 January 2008
!
!  Author:
!
!    BE Cooper
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    BE Cooper,
!    Algorithm AS 5:
!    The Integral of the Non-Central T-Distribution,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, page 193.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ST, the argument.
!
!    Input, integer ( kind = 4 ) IDF, the number of degrees of freedom.
!
!    Input, real ( kind = 8 ) D, the noncentrality parameter.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Output, real ( kind = 8 ) PRNCST, the value of the lower tail of
!    the noncentral T distribution.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) G1, 1.0 / sqrt(2.0 * pi)
!
!    Local, real ( kind = 8 ) G2, 1.0 / (2.0 * pi)
!
!    Local, real ( kind = 8 ) G3, sqrt(2.0 * pi)
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ak
  real ( kind = 8 ) alngam
  real ( kind = 8 ) alnorm
  real ( kind = 8 ) b
  real ( kind = 8 ) d
  real ( kind = 8 ) da
  real ( kind = 8 ) drb
  real ( kind = 8 ), parameter :: emin = 12.5D+00
  real ( kind = 8 ) f
  real ( kind = 8 ) fk
  real ( kind = 8 ) fkm1
  real ( kind = 8 ) fmkm1
  real ( kind = 8 ) fmkm2
  real ( kind = 8 ), parameter :: g1 = 0.3989422804D+00
  real ( kind = 8 ), parameter :: g2 = 0.1591549431D+00
  real ( kind = 8 ), parameter :: g3 = 2.5066282746D+00
  integer ( kind = 4 ) idf
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ioe
  integer ( kind = 4 ) k
  real ( kind = 8 ) prncst
  real ( kind = 8 ) rb
  real ( kind = 8 ) st
  real ( kind = 8 ) sum
  real ( kind = 8 ) tfn

  f = real ( idf, kind = 8 )
!
!  For very large IDF, use the normal approximation.
!
  if ( 100 < idf ) then

    ifault = 1

    a = sqrt ( 0.5D+00 * f ) &
    * exp ( alngam ( 0.5D+00 * ( f - 1.0D+00 ), k ) &
    - alngam ( 0.5D+00 * f, k ) ) * d

    prncst = alnorm ( ( st - a ) / sqrt ( f * ( 1.0D+00 + d * d ) &
    / ( f - 2.0D+00 ) - a * a ), .false. )

    return
  end if

  ifault = 0
  ioe = mod ( idf, 2 )
  a = st / sqrt ( f )
  b = f / ( f + st * st )
  rb = sqrt ( b )
  da = d * a
  drb = d * rb

  if ( idf == 1 ) then
    prncst = alnorm ( drb, .true. ) + 2.0D+00 * tfn ( drb, a )
    return
  end if

  sum = 0.0D+00

  if ( abs ( drb ) < emin ) then
    fmkm2 = a * rb * exp ( - 0.5D+00 * drb * drb ) &
    * alnorm ( a * drb, .false. ) * g1
  else
    fmkm2 = 0.0D+00
  end if

  fmkm1 = b * da * fmkm2

  if ( abs ( d ) < emin ) then
    fmkm1 = fmkm1 + b * a * g2 * exp ( - 0.5D+00 * d * d )
  end if

  if ( ioe == 0 ) then
    sum = fmkm2
  else
    sum = fmkm1
  end if

  ak = 1.0D+00
  fk = 2.0D+00

  do k = 2, idf - 2, 2

    fkm1 = fk - 1.0D+00
    fmkm2 = b * ( da * ak * fmkm1 + fmkm2 ) * fkm1 / fk
    ak = 1.0D+00 / ( ak * fkm1 )
    fmkm1 = b * ( da * ak * fmkm2 + fmkm1 ) * fk / ( fk + 1.0D+00 )

    if ( ioe == 0 ) then
      sum = sum + fmkm2
    else
      sum = sum + fmkm1
    end if

    ak = 1.0D+00 / ( ak * fk )
    fk = fk + 2.0D+00

  end do

  if ( ioe == 0 ) then
    prncst = alnorm ( d, .true. ) + sum * g3
  else
    prncst = alnorm ( drb, .true. ) + 2.0D+00 * ( sum + tfn ( drb, a ) )
  end if

  return
end
subroutine student_noncentral_cdf_values ( n_data, df, lambda, x, fx )

!*****************************************************************************80
!
!! STUDENT_NONCENTRAL_CDF_VALUES returns values of the noncentral Student CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NoncentralStudentTDistribution [ df, lambda ]
!      CDF [ dist, x ]
!
!    Mathematica seems to have some difficulty computing this function
!    to the desired number of digits.
!
!  Modified:
!
!    01 September 2004
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
!    Output, integer ( kind = 4 ) DF, real ( kind = 8 ) LAMBDA, the parameters
!    of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30

  integer ( kind = 4 ) df
  integer ( kind = 4 ), save, dimension ( n_max ) :: df_vec = (/ &
     1,  2,  3, &
     1,  2,  3, &
     1,  2,  3, &
     1,  2,  3, &
     1,  2,  3, &
    15, 20, 25, &
     1,  2,  3, &
    10, 10, 10, &
    10, 10, 10, &
    10, 10, 10 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.8975836176504333D+00, &
    0.9522670169D+00, &
    0.9711655571887813D+00, &
    0.8231218864D+00, &
    0.9049021510D+00, &
    0.9363471834D+00, &
    0.7301025986D+00, &
    0.8335594263D+00, &
    0.8774010255D+00, &
    0.5248571617D+00, &
    0.6293856597D+00, &
    0.6800271741D+00, &
    0.20590131975D+00, &
    0.2112148916D+00, &
    0.2074730718D+00, &
    0.9981130072D+00, &
    0.9994873850D+00, &
    0.9998391562D+00, &
    0.168610566972D+00, &
    0.16967950985D+00, &
    0.1701041003D+00, &
    0.9247683363D+00, &
    0.7483139269D+00, &
    0.4659802096D+00, &
    0.9761872541D+00, &
    0.8979689357D+00, &
    0.7181904627D+00, &
    0.9923658945D+00, &
    0.9610341649D+00, &
    0.8688007350D+00 /)
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_max ) :: lambda_vec = (/ &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00, &
    4.0D+00, &
    4.0D+00, &
    4.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
    15.00D+00, &
    15.00D+00, &
    15.00D+00, &
     0.05D+00, &
     0.05D+00, &
     0.05D+00, &
     4.00D+00, &
     4.00D+00, &
     4.00D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00, &
     6.00D+00, &
     6.00D+00, &
     6.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    df = 0
    lambda = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    df = df_vec(n_data)
    lambda = lambda_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function tfn ( x, fx )

!*****************************************************************************80
!
!! TFN calculates the T-function of Owen.
!
!  Modified:
!
!    16 January 2008
!
!  Author:
!
!    JC Young, Christoph Minder
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    MA Porter, DJ Winstanley,
!    Remark AS R30:
!    A Remark on Algorithm AS76:
!    An Integral Useful in Calculating Noncentral T and Bivariate
!    Normal Probabilities,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, page 113.
!
!    JC Young, Christoph Minder,
!    Algorithm AS 76:
!    An Algorithm Useful in Calculating Non-Central T and
!    Bivariate Normal Distributions,
!    Applied Statistics,
!    Volume 23, Number 3, 1974, pages 455-457.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, FX, the parameters of the function.
!
!    Output, real ( kind = 8 ) TFN, the value of the T-function.
!
  implicit none

  integer ( kind = 4 ), parameter :: ng = 5

  real ( kind = 8 ) fx
  real ( kind = 8 ) fxs
  integer ( kind = 4 ) i
  real ( kind = 8 ), dimension ( ng ) :: r = (/ &
    0.1477621D+00, &
    0.1346334D+00, &
    0.1095432D+00, &
    0.0747257D+00, &
    0.0333357D+00 /)
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) rt
  real ( kind = 8 ) tfn
  real ( kind = 8 ), parameter :: tp = 0.159155D+00
  real ( kind = 8 ), parameter :: tv1 = 1.0D-35
  real ( kind = 8 ), parameter :: tv2 = 15.0D+00
  real ( kind = 8 ), parameter :: tv3 = 15.0D+00
  real ( kind = 8 ), parameter :: tv4 = 1.0D-06
  real ( kind = 8 ), dimension ( ng ) :: u = (/ &
    0.0744372D+00, &
    0.2166977D+00, &
    0.3397048D+00, &
    0.4325317D+00, &
    0.4869533D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xs
!
!  Test for X near zero.
!
  if ( abs ( x ) < tv1 ) then
    tfn = tp * atan ( fx )
    return
  end if
!
!  Test for large values of abs(X).
!
  if ( tv2 < abs ( x ) ) then
    tfn = 0.0D+00
    return
  end if
!
!  Test for FX near zero.
!
  if ( abs ( fx ) < tv1 ) then
    tfn = 0.0D+00
    return
  end if
!
!  Test whether abs ( FX ) is so large that it must be truncated.
!
  xs = - 0.5D+00 * x * x
  x2 = fx
  fxs = fx * fx
!
!  Computation of truncation point by Newton iteration.
!
  if ( tv3 <= log ( 1.0D+00 + fxs ) - xs * fxs ) then

    x1 = 0.5D+00 * fx
    fxs = 0.25D+00 * fxs

    do

      rt = fxs + 1.0D+00

      x2 = x1 + ( xs * fxs + tv3 - log ( rt ) ) &
      / ( 2.0D+00 * x1 * ( 1.0D+00 / rt - xs ) )

      fxs = x2 * x2

      if ( abs ( x2 - x1 ) < tv4 ) then
        exit
      end if

      x1 = x2

    end do

  end if
!
!  Gaussian quadrature.
!
  rt = 0.0D+00

  do i = 1, ng

    r1 = 1.0D+00 + fxs * ( 0.5D+00 + u(i) )**2
    r2 = 1.0D+00 + fxs * ( 0.5D+00 - u(i) )**2

    rt = rt + r(i) * ( exp ( xs * r1 ) / r1 + exp ( xs * r2 ) / r2 )

  end do

  tfn = rt * x2 * tp

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
