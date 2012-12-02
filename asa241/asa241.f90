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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
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
function r4_normal_01_cdf_inverse ( p )

!*****************************************************************************80
!
!! R4_NORMAL_01_CDF_INVERSE inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10**7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    The Percentage Points of the Normal Distribution,
!    Algorithm AS 241,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) P, the value of the cumulative probability 
!    densitity function.  0 < P < 1.  If P is outside this range, an 
!    "infinite" value will be returned. 
!
!    Output, real ( kind = 4 ) R4_NORMAL_01_CDF_INVERSE, the normal deviate 
!    value with the property that the probability of a standard normal deviate 
!    being less than or equal to this value is P.
!
  implicit none

  real ( kind = 4 ), parameter, dimension ( 4 ) :: a = (/ &
    3.3871327179E+00, 50.434271938E+00, 159.29113202E+00, 59.109374720E+00 /)
  real ( kind = 4 ), parameter, dimension ( 4 ) :: b = (/ &
    1.0E+00, 17.895169469E+00, 78.757757664E+00, 67.187563600E+00 /)
  real ( kind = 4 ), parameter, dimension ( 4 ) :: c = (/ &
    1.4234372777E+00, 2.7568153900E+00, 1.3067284816E+00, 0.17023821103E+00 /)
  real ( kind = 4 ), parameter :: const1 = 0.180625E+00
  real ( kind = 4 ), parameter :: const2 = 1.6E+00
  real ( kind = 4 ), parameter, dimension ( 3 ) :: d = (/ &
    1.0E+00, 0.73700164250E+00, 0.12021132975E+00 /)
  real ( kind = 4 ), parameter, dimension ( 4 ) :: e = (/ &
    6.6579051150E+00, 3.0812263860E+00, 0.42868294337E+00, 0.017337203997E+00 /)
  real ( kind = 4 ), parameter, dimension ( 3 ) :: f = (/ &
    1.0E+00, 0.24197894225E+00, 0.012258202635E+00 /)
  real ( kind = 4 ) p
  real ( kind = 4 ) q
  real ( kind = 4 ) r
  real ( kind = 4 ) r4_normal_01_cdf_inverse
  real ( kind = 4 ) r4poly_value
  real ( kind = 4 ), parameter :: split1 = 0.425E+00
  real ( kind = 4 ), parameter :: split2 = 5.0E+00

  if ( p <= 0.0E+00 ) then
    r4_normal_01_cdf_inverse = - huge ( p )
    return
  end if

  if ( 1.0E+00 <= p ) then
    r4_normal_01_cdf_inverse = huge ( p )
    return
  end if

  q = p - 0.5E+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    r4_normal_01_cdf_inverse = q * r4poly_value ( 4, a, r ) &
                                 / r4poly_value ( 4, b, r )

  else

    if ( q < 0.0E+00 ) then
      r = p
    else
      r = 1.0E+00 - p
    end if

    if ( r <= 0.0E+00 ) then
      r4_normal_01_cdf_inverse = -1.0E+00
      stop
    end if

    r = sqrt ( -log ( r ) )

    if ( r <= split2 ) then

      r = r - const2
      r4_normal_01_cdf_inverse = r4poly_value ( 4, c, r ) &
                               / r4poly_value ( 3, d, r )

    else

      r = r - split2
      r4_normal_01_cdf_inverse = r4poly_value ( 4, e, r ) &
                               / r4poly_value ( 3, f, r )

    end if

    if ( q < 0.0E+00 ) then
      r4_normal_01_cdf_inverse = - r4_normal_01_cdf_inverse
    end if

  end if

  return
end
function r4poly_value ( n, a, x )

!*****************************************************************************80
!
!! R4POLY_VALUE evaluates an R4POLY.
!
!  Discussion:
!
!    For sanity's sake, the value of N indicates the NUMBER of 
!    coefficients, or more precisely, the ORDER of the polynomial,
!    rather than the DEGREE of the polynomial.  The two quantities
!    differ by 1, but cause a great deal of confusion.
!
!    Given N and A, the form of the polynomial is:
!
!      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 4 ) A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real ( kind = 4 ) X, the point at which the polynomial is 
!    to be evaluated.
!
!    Output, real ( kind = 4 ) R4POLY_VALUE, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 4 ) r4poly_value
  real ( kind = 4 ) x

  r4poly_value = 0.0E+00
  do i = n, 1, -1
    r4poly_value = r4poly_value * x + a(i)
  end do

  return
end
function r8_normal_01_cdf_inverse ( p )

!*****************************************************************************80
!
!! R8_NORMAL_01_CDF_INVERSE inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10**16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    The Percentage Points of the Normal Distribution,
!    Algorithm AS 241,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the value of the cumulative probability 
!    densitity function.  0 < P < 1.  If P is outside this range,
!    an "infinite" value will be returned.
!
!    Output, real ( kind = 8 ) D_NORMAL_01_CDF_INVERSE, the normal deviate 
!    value with the property that the probability of a standard normal 
!    deviate being less than or equal to the value is P.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
    3.3871328727963666080D+00, &
    1.3314166789178437745D+02, &
    1.9715909503065514427D+03, &
    1.3731693765509461125D+04, &
    4.5921953931549871457D+04, &
    6.7265770927008700853D+04, &
    3.3430575583588128105D+04, &
    2.5090809287301226727D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: b = (/ &
    1.0D+00, &
    4.2313330701600911252D+01, &
    6.8718700749205790830D+02, &
    5.3941960214247511077D+03, &
    2.1213794301586595867D+04, &
    3.9307895800092710610D+04, &
    2.8729085735721942674D+04, &
    5.2264952788528545610D+03 /)
  real   ( kind = 8 ), parameter, dimension ( 8 ) :: c = (/ &
    1.42343711074968357734D+00, &
    4.63033784615654529590D+00, &
    5.76949722146069140550D+00, &
    3.64784832476320460504D+00, &
    1.27045825245236838258D+00, &
    2.41780725177450611770D-01, &
    2.27238449892691845833D-02, &
    7.74545014278341407640D-04 /)
  real ( kind = 8 ), parameter :: const1 = 0.180625D+00
  real ( kind = 8 ), parameter :: const2 = 1.6D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    1.0D+00, &
    2.05319162663775882187D+00, &
    1.67638483018380384940D+00, &
    6.89767334985100004550D-01, &
    1.48103976427480074590D-01, &
    1.51986665636164571966D-02, &
    5.47593808499534494600D-04, &
    1.05075007164441684324D-09 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: e = (/ &
    6.65790464350110377720D+00, &
    5.46378491116411436990D+00, &
    1.78482653991729133580D+00, &
    2.96560571828504891230D-01, &
    2.65321895265761230930D-02, &
    1.24266094738807843860D-03, &
    2.71155556874348757815D-05, &
    2.01033439929228813265D-07 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: f = (/ &
    1.0D+00, &
    5.99832206555887937690D-01, &
    1.36929880922735805310D-01, &
    1.48753612908506148525D-02, &
    7.86869131145613259100D-04, &
    1.84631831751005468180D-05, &
    1.42151175831644588870D-07, &
    2.04426310338993978564D-15 /)
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_normal_01_cdf_inverse 
  real ( kind = 8 ) r8poly_value
  real ( kind = 8 ), parameter :: split1 = 0.425D+00
  real ( kind = 8 ), parameter :: split2 = 5.0D+00

  if ( p <= 0.0D+00 ) then
    r8_normal_01_cdf_inverse = - huge ( p )
    return
  end if

  if ( 1.0D+00 <= p ) then
    r8_normal_01_cdf_inverse = huge ( p )
    return
  end if

  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    r8_normal_01_cdf_inverse = q * r8poly_value ( 8, a, r ) &
                                 / r8poly_value ( 8, b, r )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

    if ( r <= 0.0D+00 ) then
      r8_normal_01_cdf_inverse = - 1.0D+00
      stop
    end if

    r = sqrt ( -log ( r ) )

    if ( r <= split2 ) then

      r = r - const2
      r8_normal_01_cdf_inverse = r8poly_value ( 8, c, r ) &
                               / r8poly_value ( 8, d, r )

    else

      r = r - split2
      r8_normal_01_cdf_inverse = r8poly_value ( 8, e, r ) &
                               / r8poly_value ( 8, f, r )
   
    end if

    if ( q < 0.0D+00 ) then
      r8_normal_01_cdf_inverse = - r8_normal_01_cdf_inverse
    end if

  end if

  return
end
function r8poly_value ( n, a, x )

!*****************************************************************************80
!
!! R8POLY_VALUE evaluates an R8POLY
!
!  Discussion:
!
!    For sanity's sake, the value of N indicates the NUMBER of 
!    coefficients, or more precisely, the ORDER of the polynomial,
!    rather than the DEGREE of the polynomial.  The two quantities
!    differ by 1, but cause a great deal of confusion.
!
!    Given N and A, the form of the polynomial is:
!
!      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial is 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) R8POLY_VALUE, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8poly_value
  real ( kind = 8 ) x

  r8poly_value = 0.0D+00
  do i = n, 1, -1
    r8poly_value = r8poly_value * x + a(i)
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
