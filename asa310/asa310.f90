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
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
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
subroutine beta_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

!*****************************************************************************80
!
!! BETA_NONCENTRAL_CDF_VALUES returns some values of the noncentral Beta CDF.
!
!  Discussion:
!
!    The values presented here are taken from the reference, where they
!    were given to a limited number of decimal places.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R Chattamvelli, R Shanmugam,
!    Algorithm AS 310:
!    Computing the Non-central Beta Distribution Function,
!    Applied Statistics,
!    Volume 46, Number 1, 1997, pages 146-156.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the shape parameters.
!
!    Output, real ( kind = 8 ) LAMBDA, the noncentrality parameter.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 25

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
        5.0D+00, &
        5.0D+00, &
        5.0D+00, &
       10.0D+00, &
       10.0D+00, &
       10.0D+00, &
       20.0D+00, &
       20.0D+00, &
       20.0D+00, &
       10.0D+00, &
       10.0D+00, &
       15.0D+00, &
       20.0D+00, &
       20.0D+00, &
       20.0D+00, &
       30.0D+00, &
       30.0D+00, &
       10.0D+00, &
       10.0D+00, &
       10.0D+00, &
       15.0D+00, &
       10.0D+00, &
       12.0D+00, &
       30.0D+00, &
       35.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
        5.0D+00, &
        5.0D+00, &
        5.0D+00, &
       10.0D+00, &
       10.0D+00, &
       10.0D+00, &
       20.0D+00, &
       20.0D+00, &
       20.0D+00, &
       20.0D+00, &
       10.0D+00, &
        5.0D+00, &
       10.0D+00, &
       30.0D+00, &
       50.0D+00, &
       20.0D+00, &
       40.0D+00, &
        5.0D+00, &
       10.0D+00, &
       30.0D+00, &
       20.0D+00, &
        5.0D+00, &
       17.0D+00, &
       30.0D+00, &
       30.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
       0.4563021D+00, &
       0.1041337D+00, &
       0.6022353D+00, &
       0.9187770D+00, &
       0.6008106D+00, &
       0.0902850D+00, &
       0.9998655D+00, &
       0.9925997D+00, &
       0.9641112D+00, &
       0.9376626573D+00, &
       0.7306817858D+00, &
       0.1604256918D+00, &
       0.1867485313D+00, &
       0.6559386874D+00, &
       0.9796881486D+00, &
       0.1162386423D+00, &
       0.9930430054D+00, &
       0.0506899273D+00, &
       0.1030959706D+00, &
       0.9978417832D+00, &
       0.2555552369D+00, &
       0.0668307064D+00, &
       0.0113601067D+00, &
       0.7813366615D+00, &
       0.8867126477D+00 /)
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_max ) :: lambda_vec = (/ &
        54.0D+00, &
       140.0D+00, &
       170.0D+00, &
        54.0D+00, &
       140.0D+00, &
       250.0D+00, &
        54.0D+00, &
       140.0D+00, &
       250.0D+00, &
       150.0D+00, &
       120.0D+00, &
        80.0D+00, &
       110.0D+00, &
        65.0D+00, &
       130.0D+00, &
        80.0D+00, &
       130.0D+00, &
        20.0D+00, &
        54.0D+00, &
        80.0D+00, &
       120.0D+00, &
        55.0D+00, &
        64.0D+00, &
       140.0D+00, &
        20.0D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
       0.8640D+00, &
       0.9000D+00, &
       0.9560D+00, &
       0.8686D+00, &
       0.9000D+00, &
       0.9000D+00, &
       0.8787D+00, &
       0.9000D+00, &
       0.9220D+00, &
       0.868D+00, &
       0.900D+00, &
       0.880D+00, &
       0.850D+00, &
       0.660D+00, &
       0.720D+00, &
       0.720D+00, &
       0.800D+00, &
       0.644D+00, &
       0.700D+00, &
       0.780D+00, &
       0.760D+00, &
       0.795D+00, &
       0.560D+00, &
       0.800D+00, &
       0.670D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    lambda = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    lambda = lambda_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function betain ( x, p, q, beta, ifault )

!*****************************************************************************80
!
!! BETAIN computes the incomplete Beta function ratio.
!
!  Modified:
!
!    12 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    KL Majumder, GP Bhattacharjee,
!    Algorithm AS 63:
!    The incomplete Beta Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 409-411.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument, between 0 and 1.
!
!    Input, real ( kind = 8 ) P, Q, the parameters, which
!    must be positive.
!
!    Input, real ( kind = 8 ) BETA, the logarithm of the complete
!    beta function.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error.
!    nonzero, an error occurred.
!
!    Output, real ( kind = 8 ) BETAIN, the value of the incomplete
!    Beta function ratio.
!
  implicit none

  real ( kind = 8 ), parameter :: acu = 0.1D-14
  real ( kind = 8 ) ai
  real ( kind = 8 ) beta
  real ( kind = 8 ) betain
  real ( kind = 8 ) cx
  integer ( kind = 4 ) ifault
  logical indx
  integer ( kind = 4 ) ns
  real ( kind = 8 ) p
  real ( kind = 8 ) pp
  real ( kind = 8 ) psq
  real ( kind = 8 ) q
  real ( kind = 8 ) qq
  real ( kind = 8 ) rx
  real ( kind = 8 ) temp
  real ( kind = 8 ) term
  real ( kind = 8 ) x
  real ( kind = 8 ) xx

  betain = x
  ifault = 0
!
!  Check the input arguments.
!
  if ( p <= 0.0D+00 .or. q <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
    ifault = 2
    return
  end if
!
!  Special cases.
!
  if ( x == 0.0D+00 .or. x == 1.0D+00 ) then
    return
  end if
!
!  Change tail if necessary and determine S.
!
  psq = p + q
  cx = 1.0D+00 - x

  if ( p < psq * x ) then
    xx = cx
    cx = x
    pp = q
    qq = p
    indx = .true.
  else
    xx = x
    pp = p
    qq = q
    indx = .false.
  end if

  term = 1.0D+00
  ai = 1.0D+00
  betain = 1.0D+00
  ns = int ( qq + cx * psq )
!
!  Use Soper's reduction formula.
!
  rx = xx / cx
  temp = qq - ai
  if ( ns == 0 ) then
    rx = xx
  end if

  do

    term = term * temp * rx / ( pp + ai )
    betain = betain + term
    temp = abs ( term )

    if ( temp <= acu .and. temp <= acu * betain ) then

      betain = betain * exp ( pp * log ( xx ) &
      + ( qq - 1.0D+00 ) * log ( cx ) - beta ) / pp

      if ( indx ) then
        betain = 1.0D+00 - betain
      end if

      exit

    end if

    ai = ai + 1.0D+00
    ns = ns - 1

    if ( 0 <= ns ) then
      temp = qq - ai
      if ( ns == 0 ) then
        rx = xx
      end if
    else
      temp = psq
      psq = psq + 1.0D+00
    end if

  end do

  return
end
function betanc ( x, a, b, lambda, ifault )

!*****************************************************************************80
!
!! BETANC computes the tail of the noncentral Beta distribution.
!
!  Discussion:
!
!    This routine returns the cumulative probability of X for the non-central
!    Beta distribution with parameters A, B and non-centrality LAMBDA.
!
!    Note that if LAMBDA = 0, the standard Beta distribution is defined.
!
!  Modified:
!
!    10 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Russell Lenth.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Russell Lenth,
!    Algorithm AS 226:
!    Computing Noncentral Beta Probabilities,
!    Applied Statistics,
!    Volume 36, Number 2, 1987, pages 241-244.
!
!    H Frick,
!    Algorithm AS R84:
!    A Remark on Algorithm AS 226:
!    Computing Noncentral Beta Probabilities,
!    Applied Statistics,
!    Volume 39, Number 2, 1990, pages 311-312.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value defining the cumulative
!    probability lower tail.  Normally, 0 <= X <= 1, but any value
!    is allowed.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the distribution.
!    0 < A, 0 < B.
!
!    Input, real ( kind = 8 ) LAMBDA, the noncentrality parameter
!    of the distribution.  0 <= LAMBDA.  The program can produce reasonably
!    accurate results for values of LAMBDA up to about 100.
!
!    Output, integer IFAULT, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Output, real ( kind = 8 ) BETANC, the cumulative probability
!    of X.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a0
  real ( kind = 8 ) alngam
  real ( kind = 8 ) alnorm
  real ( kind = 8 ) ax
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) betain
  real ( kind = 8 ) betanc
  real ( kind = 8 ) c
  real ( kind = 8 ) errbd
  real ( kind = 8 ), parameter :: errmax = 1.0D-07
  real ( kind = 8 ) gx
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ), parameter :: itrmax = 150
  real ( kind = 8 ) lambda
  real ( kind = 8 ) q
  real ( kind = 8 ) sumq
  real ( kind = 8 ) temp
  real ( kind = 8 ), parameter :: ualpha = 5.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) xj

  ifault = 0

  if ( lambda < 0.0D+00 .or. &
       a <= 0.0D+00 .or. &
       b <= 0.0D+00 ) then
    ifault = 2
    betanc = -1.0D+00
    return
  end if

  if ( x <= 0.0D+00 ) then
    betanc = 0.0D+00
    return
  end if

  if ( 1.0D+00 <= x ) then
    betanc = 1.0D+00
    return
  end if

  c = 0.5D+00 * lambda
!
!  Initialize the series.
!
  beta = alngam ( a, ifault ) &
       + alngam ( b, ifault ) &
       - alngam ( a + b, ifault )

  temp = betain ( x, a, b, beta, ifault )

  gx = exp ( a * log ( x ) + b * log ( 1.0D+00 - x ) &
    - beta - log ( a ) )

  q = exp ( - c )

  xj = 0.0D+00
  ax = q * temp
  sumq = 1.0D+00 - q
  betanc = ax
!
!  Recur over subsequent terms until convergence is achieved.
!
  ifault = 1

  do

    xj = xj + 1.0D+00
    temp = temp - gx
    gx = x * ( a + b + xj - 1.0D+00 ) * gx / ( a + xj )
    q = q * c / xj
    sumq = sumq - q
    ax = temp * q
    betanc = betanc + ax
!
!  Check for convergence and act accordingly.
!
    errbd = abs ( ( temp - gx ) * sumq )

    if ( errbd <= errmax ) then
      ifault = 0
      exit
    end if

    if (  itrmax < int ( xj ) ) then
      exit
    end if

  end do

  return
end
function gammad ( x, p, ifault )

!*****************************************************************************80
!
!! GAMMAD computes the Incomplete Gamma Integral
!
!  Discussion:
!
!    This routine uses the auxiliary functions:
!
!      ALNGAM = logarithm of the gamma function,
!      ALNORM = algorithm AS66
!
!  Modified:
!
!    31 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by B Shea.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    B Shea,
!    Algorithm AS 239:
!    Chi-squared and Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 466-473.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, P, the parameters.
!    0 <= X, and 0 < P.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error.
!    1, X < 0 or P <= 0.
!
!    Output, real ( kind = 8 ) GAMMAD, the value of the incomplete
!    Gamma integral.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alngam
  real ( kind = 8 ) alnorm
  real ( kind = 8 ) an
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: elimit = - 88.0D+00
  real ( kind = 8 ) gammad
  integer ( kind = 4 ) ifault
  real ( kind = 8 ), parameter :: oflo = 1.0D+37
  real ( kind = 8 ) p
  real ( kind = 8 ), parameter :: plimit = 1000.0D+00
  real ( kind = 8 ) pn1
  real ( kind = 8 ) pn2
  real ( kind = 8 ) pn3
  real ( kind = 8 ) pn4
  real ( kind = 8 ) pn5
  real ( kind = 8 ) pn6
  real ( kind = 8 ) rn
  real ( kind = 8 ), parameter :: tol = 1.0D-14
  logical upper
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 1.0D+08

  gammad = 0.0D+00
!
!  Check the input.
!
  if ( x < 0.0D+00 ) then
    ifault = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMAD - Fatal error!'
    write ( *, '(a)' ) '  X < 0.'
    stop
  end if

  if ( p <= 0.0D+00 ) then
    ifault = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMAD - Fatal error!'
    write ( *, '(a)' ) '  P <= 0.'
    stop
  end if

  ifault = 0

  if ( x == 0.0D+00 ) then
    gammad = 0.0D+00
    return
  end if
!
!  If P is large, use a normal approximation.
!
  if ( plimit < p ) then

    pn1 = 3.0D+00 * sqrt ( p ) &
    * ( ( x / p )**( 1.0D+00 / 3.0D+00 ) + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )

    upper = .false.
    gammad = alnorm ( pn1, upper )
    return

  end if
!
!  If X is large set GAMMAD = 1.
!
  if ( xbig < x ) then
    gammad = 1.0D+00
    return
  end if
!
!  Use Pearson's series expansion.
!  (Note that P is not large enough to force overflow in ALOGAM).
!  No need to test IFAULT on exit since 0 < P.
!
  if ( x <= 1.0D+00 .or. x < p ) then

    arg = p * log ( x ) - x - alngam ( p + 1.0D+00, ifault )
    c = 1.0D+00
    gammad = 1.0D+00
    a = p

    do

      a = a + 1.0D+00
      c = c * x / a
      gammad = gammad + c

      if ( c <= tol ) then
        exit
      end if

    end do

    arg = arg + log ( gammad )

    if ( elimit <= arg ) then
      gammad = exp ( arg )
    else
      gammad = 0.0D+00
    end if
!
!  Use a continued fraction expansion.
!
  else

    arg = p * log ( x ) - x - alngam ( p, ifault )
    a = 1.0D+00 - p
    b = a + x + 1.0D+00
    c = 0.0D+00
    pn1 = 1.0D+00
    pn2 = x
    pn3 = x + 1.0E+00
    pn4 = x * b
    gammad = pn3 / pn4

    do

      a = a + 1.0D+00
      b = b + 2.0D+00
      c = c + 1.0D+00
      an = a * c
      pn5 = b * pn3 - an * pn1
      pn6 = b * pn4 - an * pn2

      if ( pn6 /= 0.0D+00 ) then

        rn = pn5 / pn6

        if ( abs ( gammad - rn ) <= min ( tol, tol * rn ) ) then
          exit
        end if

        gammad = rn

      end if

      pn1 = pn3
      pn2 = pn4
      pn3 = pn5
      pn4 = pn6
!
!  Re-scale terms in continued fraction if terms are large.
!
      if ( oflo <= abs ( pn5 ) ) then
        pn1 = pn1 / oflo
        pn2 = pn2 / oflo
        pn3 = pn3 / oflo
        pn4 = pn4 / oflo
      end if

    end do

    arg = arg + log ( gammad )

    if ( elimit <= arg ) then
      gammad = 1.0D+00 - exp ( arg )
    else
      gammad = 1.0D+00
    end if

  end if

  return
end
function ncbeta ( a, b, lambda, x, errmax, ifault )

!*****************************************************************************80
!
!! NCBETA computes the noncentral Beta CDF.
!
!  Discussion:
!
!    Three corrections needed to be made to the text of this routine.
!    They are noted in the comments below.
!
!    Two of these corrections were errors in transcription made when
!    producing the online copy distributed by APSTAT.
!
!    One error, an error of omission, occurred in the original printed
!    copy of the routine, and was carried over into the online copy.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by R Chattamvelli, R Shanmugam.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    R Chattamvelli, R Shanmugam,
!    Algorithm AS 310:
!    Computing the Non-central Beta Distribution Function,
!    Applied Statistics,
!    Volume 46, Number 1, 1997, pages 146-156.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the shape parameters.
!    0 <= A, 0 <= B.
!
!    Input, real ( kind = 8 ) LAMBDA, the noncentrality parameter.
!    0 <= LAMBDA.
!
!    Input, real ( kind = 8 ) X, the value at which the CDF is desired.
!
!    Input, real ( kind = 8 ) ERRMAX, the precision tolerance.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!    1, X is 0 or 1.
!    2, X < 0 or 1 < X.
!    3, A, B or LAMBDA is less than 0.
!
!    Output, double precision NCBETA, the cumulative probability of X.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), external :: alngam
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ), external :: betain
  real ( kind = 8 ), external :: betanc
  real ( kind = 8 ) c
  real ( kind = 8 ) ebd
  real ( kind = 8 ) errbd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) ftemp
  real ( kind = 8 ) fx
  real ( kind = 8 ), external :: gammad
  real ( kind = 8 ) gx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) iter1
  integer ( kind = 4 ) iter2
  integer ( kind = 4 ) iterhi
  integer ( kind = 4 ) iterlo
  integer ( kind = 4 ) j
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) m
  real ( kind = 8 ) ncbeta
  real ( kind = 8 ) mr
  real ( kind = 8 ) psum
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) s0
  real ( kind = 8 ) s1
  real ( kind = 8 ) sum
  real ( kind = 8 ) t
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  real ( kind = 8 ) temp
  real ( kind = 8 ) x
  integer ( kind = 4 ) xj

  ifault = 0
  ncbeta = x
!
!  Check parameters.
!
  if ( lambda <= 0.0D+00 ) then
    ifault = 3
    return
  end if

  if ( a <= 0.0D+00 ) then
    ifault = 3
    return
  end if

  if ( b <= 0.0D+00 ) then
    ifault = 3
    return
  end if

  if ( x <= 0.0D+00 ) then
    ncbeta = 0.0D+00
    return
  end if

  if ( 1.0D+00 <= x ) then
    ncbeta = 1.0D+00
    return
  end if

  c = 0.5D+00 * lambda
  xj = 0.0D+00
!
!  AS 226 as it stands is sufficient in this situation.
!
  if ( lambda < 54.0D+00 ) then

    ncbeta = betanc ( x, a, b, lambda, ifault )
    return

  else

    m = int ( c + 0.5D+00 )
    mr = real ( m, kind = 8 )
    iterlo = m - int ( 5.0D+00 * sqrt ( mr ) )
    iterhi = m + int ( 5.0D+00 * sqrt ( mr ) )
    t = - c + mr * log ( c ) - alngam ( mr + 1.0D+00, ifault )
    q = exp ( t )
    r = q
    psum = q

    beta = alngam ( a + mr, ifault ) &
         + alngam ( b, ifault ) &
         - alngam ( a + mr + b, ifault )

    s1 = ( a + mr ) * log ( x ) &
       + b * log ( 1.0D+00 - x ) - log ( a + mr ) - beta
    gx = exp ( s1 )
    fx = gx
    temp = betain ( x, a + mr, b, beta, ifault )
    ftemp = temp
    xj = xj + 1.0D+00
!
!  The online copy of AS 310 has "SUM = Q - TEMP" which is incorrect.
!
    sum = q * temp
    iter1 = m
!
!  The first set of iterations starts from M and goes downwards
!
    do

      if ( iter1 < iterlo ) then
        exit
      end if

      if ( q < errmax ) then
        exit
      end if
!
!  The online copy of AS 310 has "Q = Q - ITER1 / C" which is incorrect.
!
      q = q * iter1 / c
      xj = xj + 1.0D+00
      gx = ( a + iter1 ) / ( x * ( a + b + iter1 - 1.0D+00 ) ) * gx
      iter1 = iter1 - 1
      temp = temp + gx
      psum = psum + q
      sum = sum + q * temp

    end do

    t0 = alngam ( a + b, ifault ) &
       - alngam ( a + 1.0D+00, ifault ) &
       - alngam ( b, ifault )

    s0 = a * log ( x ) + b * log ( 1.0D+00 - x )
!
!  Both the online copy of AS 310 and the text printed in the reference
!  did not initialize the variable S to zero, which is incorrect.
!  JVB, 12 January 2008.
!
    s = 0.0
    do i = 1, iter1
      j = i - 1
      s = s + exp ( t0 + s0 + j * log ( x ) )
      t1 = log ( a + b + j ) - log ( a + 1.0D+00 + j ) + t0
      t0 = t1
    end do
!
!  Compute the first part of error bound.
!
    errbd = ( 1.0D+00 - gammad ( c, real ( iter1, kind = 8 ), ifault ) ) &
    * ( temp + s )

    q = r
    temp = ftemp
    gx = fx
    iter2 = m

    do

      ebd = errbd + ( 1.0D+00 - psum ) * temp

      if ( ebd < errmax .or. iterhi <= iter2 ) then
        exit
      end if

      iter2 = iter2 + 1
      xj = xj + 1.0D+00
      q = q * c / iter2
      psum = psum + q
      temp = temp - gx
      gx = x * ( a + b + iter2 - 1.0D+00 ) / ( a + iter2 ) * gx
      sum = sum + q * temp

    end do

  end if

  ncbeta = sum

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
