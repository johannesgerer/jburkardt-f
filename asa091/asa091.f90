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
subroutine chi_square_cdf_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = ChiSquareDistribution [ df ]
!      CDF [ dist, x ]
!
!  Modified:
!
!    09 August 2004
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
!    Output, integer ( kind = 4 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  integer ( kind = 4 ) a
  integer ( kind = 4 ), save, dimension ( n_max ) :: a_vec = (/ &
     1,  2,  1,  2, &
     1,  2,  3,  4, &
     1,  2,  3,  4, &
     5,  3,  3,  3, &
     3,  3, 10, 10, &
    10 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7965567455405796D-01, &
    0.4987520807317687D-02, &
    0.1124629160182849D+00, &
    0.9950166250831946D-02, &
    0.4729107431344619D+00, &
    0.1812692469220181D+00, &
    0.5975750516063926D-01, &
    0.1752309630642177D-01, &
    0.6826894921370859D+00, &
    0.3934693402873666D+00, &
    0.1987480430987992D+00, &
    0.9020401043104986D-01, &
    0.3743422675270363D-01, &
    0.4275932955291202D+00, &
    0.6083748237289110D+00, &
    0.7385358700508894D+00, &
    0.8282028557032669D+00, &
    0.8883897749052874D+00, &
    0.1721156299558408D-03, &
    0.3659846827343712D-02, &
    0.1857593622214067D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.01D+00, &
    0.01D+00, &
    0.02D+00, &
    0.02D+00, &
    0.40D+00, &
    0.40D+00, &
    0.40D+00, &
    0.40D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    3.00D+00, &
    4.00D+00, &
    5.00D+00, &
    6.00D+00, &
    1.00D+00, &
    2.00D+00, &
    3.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function gammad ( x, p, ifault )

!*****************************************************************************80
!
!! GAMMAD computes the Incomplete Gamma Integral
!
!  Auxiliary functions:
!
!    ALOGAM = logarithm of the gamma function,
!    ALNORM = algorithm AS66
!
!  Modified:
!
!    20 January 2008
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
!    Input, real ( kind = 8 ) X, P, the parameters of the incomplete
!    gamma ratio.  0 <= X, and 0 < P.
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
  real ( kind = 8 ) alnorm
  real ( kind = 8 ) alngam
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
    return
  end if

  if ( p <= 0.0D+00 ) then
    ifault = 1
    return
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

    pn1 = 3.0D+00 * sqrt ( p ) * ( ( x / p )**( 1.0D+00 / 3.0D+00 ) &
    + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )

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
!  No need to test IFAULT on exit since P > 0.
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
    pn3 = x + 1.0D+00
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
function ppchi2 ( p, v, g, ifault )

!*****************************************************************************80
!
!! PPCHI2 evaluates the percentage points of the Chi-squared PDF.
!
!  Discussion
!
!    Incorporates the suggested changes in AS R85 (vol.40(1),
!    pages 233-5, 1991) which should eliminate the need for the limited
!    range for P, though these limits have not been removed
!    from the routine.
!
!  Modified:
!
!    01 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Donald Best, DE Roberts.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Donald Best, DE Roberts,
!    Algorithm AS 91:
!    The Percentage Points of the Chi-Squared Distribution,
!    Applied Statistics,
!    Volume 24, Number 3, 1975, pages 385-390.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P,  value of the chi-squared cumulative
!    probability density function.
!    0.000002 <= P <= 0.999998.
!
!    Input, real ( kind = 8 ) V, the parameter of the chi-squared probability
!    density function.
!    0 < V.
!
!    Input, real ( kind = 8 ) G, the value of log ( Gamma ( V / 2 ) ).
!
!    Output, integer IFAULT, is nonzero if an error occurred.
!    0, no error.
!    1, P is outside the legal range.
!    2, V is not positive.
!    3, an error occurred in GAMMAD.
!    4, the result is probably as accurate as the machine will allow.
!
!    Output, real ( kind = 8 ) PPCHI2, the value of the chi-squared random
!    deviate with the property that the probability that a chi-squared random
!    deviate with parameter V is less than or equal to PPCHI2 is P.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), parameter :: aa = 0.6931471806D+00
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: c1 = 0.01D+00
  real ( kind = 8 ), parameter :: c2 = 0.222222D+00
  real ( kind = 8 ), parameter :: c3 = 0.32D+00
  real ( kind = 8 ), parameter :: c4 = 0.4D+00
  real ( kind = 8 ), parameter :: c5 = 1.24D+00
  real ( kind = 8 ), parameter :: c6 = 2.2D+00
  real ( kind = 8 ), parameter :: c7 = 4.67D+00
  real ( kind = 8 ), parameter :: c8 = 6.66D+00
  real ( kind = 8 ), parameter :: c9 = 6.73D+00
  real ( kind = 8 ), parameter :: c10 = 13.32D+00
  real ( kind = 8 ), parameter :: c11 = 60.0D+00
  real ( kind = 8 ), parameter :: c12 = 70.0D+00
  real ( kind = 8 ), parameter :: c13 = 84.0D+00
  real ( kind = 8 ), parameter :: c14 = 105.0D+00
  real ( kind = 8 ), parameter :: c15 = 120.0D+00
  real ( kind = 8 ), parameter :: c16 = 127.0D+00
  real ( kind = 8 ), parameter :: c17 = 140.0D+00
  real ( kind = 8 ), parameter :: c18 = 175.0D+00
  real ( kind = 8 ), parameter :: c19 = 210.0D+00
  real ( kind = 8 ), parameter :: c20 = 252.0D+00
  real ( kind = 8 ), parameter :: c21 = 264.0D+00
  real ( kind = 8 ), parameter :: c22 = 294.0D+00
  real ( kind = 8 ), parameter :: c23 = 346.0D+00
  real ( kind = 8 ), parameter :: c24 = 420.0D+00
  real ( kind = 8 ), parameter :: c25 = 462.0D+00
  real ( kind = 8 ), parameter :: c26 = 606.0D+00
  real ( kind = 8 ), parameter :: c27 = 672.0D+00
  real ( kind = 8 ), parameter :: c28 = 707.0D+00
  real ( kind = 8 ), parameter :: c29 = 735.0D+00
  real ( kind = 8 ), parameter :: c30 = 889.0D+00
  real ( kind = 8 ), parameter :: c31 = 932.0D+00
  real ( kind = 8 ), parameter :: c32 = 966.0D+00
  real ( kind = 8 ), parameter :: c33 = 1141.0D+00
  real ( kind = 8 ), parameter :: c34 = 1182.0D+00
  real ( kind = 8 ), parameter :: c35 = 1278.0D+00
  real ( kind = 8 ), parameter :: c36 = 1740.0D+00
  real ( kind = 8 ), parameter :: c37 = 2520.0D+00
  real ( kind = 8 ), parameter :: c38 = 5040.0D+00
  real ( kind = 8 ) ch
  real ( kind = 8 ), parameter :: e = 0.5D-06
  real ( kind = 8 ) g
  real ( kind = 8 ) gammad
  integer ( kind = 4 ) i
  integer ( kind = 4 ) if1
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ), parameter :: maxit = 20
  real ( kind = 8 ) ppchi2
  real ( kind = 8 ) p
  real ( kind = 8 ), parameter :: pmax = 0.999998D+00
  real ( kind = 8 ), parameter :: pmin = 0.000002D+00
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) ppnd
  real ( kind = 8 ) q
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) s3
  real ( kind = 8 ) s4
  real ( kind = 8 ) s5
  real ( kind = 8 ) s6
  real ( kind = 8 ) t
  real ( kind = 8 ) v
  real ( kind = 8 ) x
  real ( kind = 8 ) xx
!
!  Test arguments and initialize.
!
  ppchi2 = - 1.0D+00

  if ( p < pmin .or. pmax < p ) then
    ifault = 1
    return
  end if

  if ( v <= 0.0D+00 ) then
    ifault = 2
    return
  end if

  ifault = 0
  xx = 0.5D+00 * v
  c = xx - 1.0D+00
!
!  Starting approximation for small chi-squared
!
  if ( v < - c5 * log ( p ) ) then

    ch = ( p * xx * exp ( g + xx * aa ) ) ** ( 1.0D+00 / xx )

    if ( ch < e ) then
      ppchi2 = ch
      return
    end if
!
!  Starting approximation for V less than or equal to 0.32
!
  else if ( v <= c3 ) then

    ch = c4
    a = log ( 1.0D+00 - p )

    do

      q = ch
      p1 = 1.0D+00 + ch * ( c7 + ch )
      p2 = ch * (c9 + ch * ( c8 + ch ) )

      t = - 0.5D+00 + (c7 + 2.0D+00 * ch ) / p1 - ( c9 + ch * ( c10 + &
      3.0D+00 * ch ) ) / p2

      ch = ch - ( 1.0D+00 - exp ( a + g + 0.5D+00 * ch + c * aa ) * &
      p2 / p1) / t

      if ( abs ( q / ch - 1.0D+00 ) <= c1 ) then
        exit
      end if

    end do

  else
!
!  Call to algorithm AS 111 - note that P has been tested above.
!  AS 241 could be used as an alternative.
!
    x = ppnd ( p, ifault )
!
!  Starting approximation using Wilson and Hilferty estimate
!
    p1 = c2 / v
    ch = v * ( x * sqrt ( p1 ) + 1.0D+00 - p1)**3
!
!  Starting approximation for P tending to 1.
!
    if ( c6 * v + 6.0D+00 < ch ) then
       ch = - 2.0D+00 * ( log ( 1.0D+00 - p ) - c * log ( 0.5D+00 * ch ) + g )
    end if

  end if
!
!  Call to algorithm AS 239 and calculation of seven term
!  Taylor series
!
  do i = 1, maxit

    q = ch
    p1 = 0.5D+00 * ch
    p2 = p - gammad ( p1, xx, if1 )

    if ( if1 /= 0 ) then
      ifault = 3
      return
    end if

    t = p2 * exp ( xx * aa + g + p1 - c * log ( ch ) )
    b = t / ch
    a = 0.5D+00 * t - b * c
    s1 = ( c19 + a * ( c17 + a * ( c14 + a * ( c13 + a * ( c12 + &
    c11 * a ))))) / c24
    s2 = ( c24 + a * ( c29 + a * ( c32 + a * ( c33 + c35 * a )))) / c37
    s3 = ( c19 + a * ( c25 + a * ( c28 + c31 * a ))) / c37
    s4 = ( c20 + a * ( c27 + c34 * a) + c * ( c22 + a * ( c30 + c36 * a ))) &
    / c38
    s5 = ( c13 + c21 * a + c * ( c18 + c26 * a )) / c37
    s6 = ( c15 + c * ( c23 + c16 * c )) / c38
    ch = ch + t * ( 1.0D+00 + 0.5D+00 * t * s1 - b * c * ( s1 - b * &
    ( s2 - b * ( s3 - b * ( s4 - b * ( s5 - b * s6 ))))))

    if ( e < abs ( q / ch - 1.0D+00 ) ) then
       ppchi2 = ch
       return
    end if

  end do

 ifault = 4
 ppchi2 = ch

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
!    Original FORTRAN77 version by J Beasley, S Springer.
!    FORTRAN90 version by John Burkardt.
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
