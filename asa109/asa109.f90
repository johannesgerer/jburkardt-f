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
subroutine beta_inc_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! BETA_INC_VALUES returns some values of the incomplete Beta function.
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
function xinbta ( p, q, beta, alpha, ifault )

!*****************************************************************************80
!
!! XINBTA computes inverse of the incomplete Beta function.
!
!  Modified:
!
!    09 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by GW Cran, KJ Martin, GE Thomas.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    GW Cran, KJ Martin, GE Thomas,
!    Remark AS R19 and Algorithm AS 109:
!    A Remark on Algorithms AS 63: The Incomplete Beta Integral
!    and AS 64: Inverse of the Incomplete Beta Integeral,
!    Applied Statistics,
!    Volume 26, Number 1, 1977, pages 111-114.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, Q, the parameters of the incomplete
!    Beta function.
!
!    Input, real ( kind = 8 ) BETA, the logarithm of the value of
!    the complete Beta function.
!
!    Input, real ( kind = 8 ) ALPHA, the value of the incomplete Beta
!    function.  0 <= ALPHA <= 1.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Output, real ( kind = 8 ) XINBTA, the argument of the incomplete
!    Beta function which produces the value ALPHA.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) SAE, the most negative decimal exponent
!    which does not cause an underflow.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) acu
  real ( kind = 8 ) adj
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) betain
  real ( kind = 8 ) fpu
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) iex
  integer ( kind = 4 ) ifault
  logical indx
  real ( kind = 8 ) p
  real ( kind = 8 ) pp
  real ( kind = 8 ) prev
  real ( kind = 8 ) q
  real ( kind = 8 ) qq
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ), parameter :: sae = -37.0D+00
  real ( kind = 8 ) sq
  real ( kind = 8 ) t
  real ( kind = 8 ) tx
  real ( kind = 8 ) w
  real ( kind = 8 ) xin
  real ( kind = 8 ) xinbta
  real ( kind = 8 ) y
  real ( kind = 8 ) yprev

  fpu = 10.0D+00**sae

  ifault = 0
  xinbta = alpha
!
!  Test for admissibility of parameters.
!
  if ( p <= 0.0D+00 .or. q <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  if ( alpha < 0.0D+00 .or. 1.0D+00 < alpha ) then
    ifault = 2
    return
  end if

  if ( alpha == 0.0D+00 .or. alpha == 1.0D+00 ) then
    return
  end if
!
!  Change tail if necessary.
!
  if ( 0.5D+00 < alpha ) then
    a = 1.0D+00 - alpha
    pp = q
    qq = p
    indx = .true.
  else
    a = alpha
    pp = p
    qq = q
    indx = .false.
  end if
!
!  Calculate the initial approximation.
!
  r = sqrt ( - log ( a * a ) )

  y = r - ( 2.30753D+00 + 0.27061D+00 * r ) &
    / ( 1.0D+00 + ( 0.99229D+00 + 0.04481D+00 * r ) * r )

  if ( 1.0D+00 < pp .and. 1.0D+00 < qq ) then

    r = ( y * y - 3.0D+00 ) / 6.0D+00
    s = 1.0D+00 / ( pp + pp - 1.0D+00 )
    t = 1.0D+00 / ( qq + qq - 1.0D+00 )
    h = 2.0D+00 / ( s + t )
    w = y * sqrt ( h + r ) / h - ( t - s ) &
    * ( r + 5.0D+00 / 6.0D+00 - 2.0D+00 / ( 3.0D+00 * h ) )
    xinbta = pp / ( pp + qq * exp ( w + w ) )

  else

    r = qq + qq
    t = 1.0D+00 / ( 9.0D+00 * qq )
    t = r * ( 1.0D+00 - t + y * sqrt ( t ) )**3

    if ( t <= 0.0D+00 ) then
      xinbta = 1.0D+00 - exp ( ( log ( ( 1.0D+00 - a ) * qq ) + beta ) / qq )
    else

      t = ( 4.0D+00 * pp + r - 2.0D+00 ) / t

      if ( t <= 1.0D+00 ) then
        xinbta = exp ( ( log ( a * pp ) + beta ) / pp )
      else
        xinbta = 1.0D+00 - 2.0D+00 / ( t + 1.0D+00 )
      end if

    end if

  end if
!
!  Solve for X by a modified Newton-Raphson method,
!  using the function BETAIN.
!
  r = 1.0D+00 - pp
  t = 1.0D+00 - qq
  yprev = 0.0D+00
  sq = 1.0D+00
  prev = 1.0D+00

  if ( xinbta < 0.0001D+00 ) then
    xinbta = 0.0001D+00
  end if

  if ( 0.9999D+00 < xinbta ) then
    xinbta = 0.9999D+00
  end if

  iex = max ( - 5.0D+00 / pp**2 - 1.0D+00 / a**0.2D+00 - 13.0D+00, sae )

  acu = 10.0D+00**iex

  do

    y = betain ( xinbta, pp, qq, beta, ifault )

    if ( ifault /= 0 ) then
      ifault = 3
      return
    end if

    xin = xinbta
    y = ( y - a ) * exp ( beta + r * log ( xin ) + t * log ( 1.0D+00 - xin ) )

    if ( y * yprev <= 0.0D+00 ) then
      prev = max ( sq, fpu )
    end if

    g = 1.0D+00

    do

      do

        adj = g * y
        sq = adj * adj

        if ( sq < prev ) then

          tx = xinbta - adj

          if ( 0.0D+00 <= tx .and. tx <= 1.0D+00 ) then
            exit
          end if

        end if

        g = g / 3.0D+00

      end do

      if ( prev <= acu ) then
        if ( indx ) then
          xinbta = 1.0D+00 - xinbta
        end if
        return
      end if

      if ( y * y <= acu ) then
        if ( indx ) then
          xinbta = 1.0D+00 - xinbta
        end if
        return
      end if

      if ( tx /= 0.0D+00 .and. tx /= 1.0D+00 ) then
        exit
      end if

      g = g / 3.0D+00

    end do

    if ( tx == xinbta ) then
      exit
    end if

    xinbta = tx
    yprev = y

  end do

  if ( indx ) then
    xinbta = 1.0D+00 - xinbta
  end if

  return
end
