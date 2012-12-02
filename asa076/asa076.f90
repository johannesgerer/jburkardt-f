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
subroutine owen_values ( n_data, h, a, t )

!*****************************************************************************80
!
!! OWEN_VALUES returns some values of Owen's T function.
!
!  Discussion:
!
!    Owen's T function is useful for computation of the bivariate normal
!    distribution and the distribution of a skewed normal distribution.
!
!    Although it was originally formulated in terms of the bivariate
!    normal function, the function can be defined more directly as
!
!      T(H,A) = 1 / ( 2 * pi ) *
!        Integral ( 0 <= X <= A ) e^(-H^2*(1+X^2)/2) / (1+X^2) dX
!
!    In Mathematica, the function can be evaluated by:
!
!      fx = 1/(2*Pi) * Integrate [ E^(-h^2*(1+x^2)/2)/(1+x^2), {x,0,a} ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mike Patefield, David Tandy,
!    Fast and Accurate Calculation of Owen's T Function,
!    Journal of Statistical Software,
!    Volume 5, Number 5, 2000, pages 1-25.
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
!    Output, real ( kind = 8 ) H, a parameter.
!
!    Output, real ( kind = 8 ) A, the upper limit of the integral.
!
!    Output, real ( kind = 8 ) T, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.2500000000000000D+00, &
    0.4375000000000000D+00, &
    0.9687500000000000D+00, &
    0.0625000000000000D+00, &
    0.5000000000000000D+00, &
    0.9999975000000000D+00, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.1000000000000000D+02, &
    0.1000000000000000D+03 /)
  real ( kind = 8 ) h
  real ( kind = 8 ), save, dimension ( n_max ) :: h_vec = (/ &
    0.0625000000000000D+00, &
    6.5000000000000000D+00, &
    7.0000000000000000D+00, &
    4.7812500000000000D+00, &
    2.0000000000000000D+00, &
    1.0000000000000000D+00, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) t
  real ( kind = 8 ), save, dimension ( n_max ) :: t_vec = (/ &
    3.8911930234701366D-02, &
    2.0005773048508315D-11, &
    6.3990627193898685D-13, &
    1.0632974804687463D-07, &
    8.6250779855215071D-03, &
    6.6741808978228592D-02, &
    0.4306469112078537D-01, &
    0.6674188216570097D-01, &
    0.7846818699308410D-01, &
    0.7929950474887259D-01, &
    0.6448860284750376D-01, &
    0.1066710629614485D+00, &
    0.1415806036539784D+00, &
    0.1510840430760184D+00, &
    0.7134663382271778D-01, &
    0.1201285306350883D+00, &
    0.1666128410939293D+00, &
    0.1847501847929859D+00, &
    0.7317273327500385D-01, &
    0.1237630544953746D+00, &
    0.1737438887583106D+00, &
    0.1951190307092811D+00, &
    0.7378938035365546D-01, &
    0.1249951430754052D+00, &
    0.1761984774738108D+00, &
    0.1987772386442824D+00, &
    0.2340886964802671D+00, &
    0.2479460829231492D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    h = 0.0D+00
    a = 0.0D+00
    t = 0.0D+00
  else
    h = h_vec(n_data)
    a = a_vec(n_data)
    t = t_vec(n_data)
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
!    Original FORTRAN77 version by JC Young, Christoph Minder.
!    FORTRAN90 version by John Burkardt.
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
  real ( kind = 8 ), parameter :: tv4 = 1.0D-05
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
function tha ( h1, h2, a1, a2 )

!*****************************************************************************80
!
!! THA computes Owen's T function.
!
!  Discussion:
!
!    This function computes T(H1/H2, A1/A2) for any real numbers H1, H2,
!    A1 and A2.
!
!  Modified:
!
!    16 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by JC Young, Christoph Minder.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Boys,
!    Remark AS R80:
!    A Remark on Algorithm AS76:
!    An Integral Useful in Calculating Noncentral T and Bivariate
!    Normal Probabilities,
!    Applied Statistics,
!    Volume 38, Number 3, 1989, pages 580-582.
!
!    Youn-Min Chou,
!    Remark AS R55:
!    A Remark on Algorithm AS76:
!    An Integral Useful in Calculating Noncentral T and Bivariate
!    Normal Probabilities,
!    Applied Statistics,
!    Volume 34, Number 1, 1985, pages 100-101.
!
!    PW Goedhart, MJW Jansen,
!    Remark AS R89:
!    A Remark on Algorithm AS76:
!    An Integral Useful in Calculating Noncentral T and Bivariate
!    Normal Probabilities,
!    Applied Statistics,
!    Volume 41, Number 2, 1992, pages 496-497.
!
!    JC Young, Christoph Minder,
!    Algorithm AS 76:
!    An Algorithm Useful in Calculating Noncentral T and
!    Bivariate Normal Distributions,
!    Applied Statistics,
!    Volume 23, Number 3, 1974, pages 455-457.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) H1, H2, A1, A2, define the arguments
!    of the T function.
!
!    Output, real ( kind = 8 ) THA, the value of Owen's T function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) absa
  real ( kind = 8 ) ah
  real ( kind = 8 ) alnorm
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) ex
  real ( kind = 8 ) g
  real ( kind = 8 ) gah
  real ( kind = 8 ) gh
  real ( kind = 8 ) h
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  real ( kind = 8 ) lam
  real ( kind = 8 ) tfn
  real ( kind = 8 ) tha
  real ( kind = 8 ), parameter :: twopi = 6.2831853071795864769D+00

  if ( h2 == 0.0D+00 ) then
    tha = 0.0D+00
    return
  end if

  h = h1 / h2

  if ( a2 == 0.0D+00 ) then

    g = alnorm ( h, .false. )

    if ( h < 0.0D+00 ) then
      tha = g / 2.0D+00
    else
      tha = ( 1.0D+00 - g ) / 2.0D+00
    end if

    if ( a1 < 0.0D+00 ) then
      tha = - tha
    end if

    return
  end if

  a = a1 / a2

  if ( abs ( h ) < 0.3D+00 .and. 7.0D+00 < abs ( a ) ) then

    lam = abs ( a * h )
    ex = exp ( - lam * lam / 2.0D+00 )
    g = alnorm ( lam, .false. )
    c1 = ( ex / lam + sqrt ( twopi ) * ( g - 0.5D+00 ) ) / twopi
    c2 = ( ( lam * lam + 2.0D+00 ) * ex / lam**3 &
    + sqrt ( twopi ) * ( g - 0.5D+00 ) ) / ( 6.0D+00 * twopi )
    ah = abs ( h )
    tha = 0.25D+00 - c1 * ah + c2 * ah**3
    tha = sign ( tha, a )

  else
!
!  Correction AS R89
!
    absa = abs ( a )

    if ( absa <= 1.0D+00 ) then
      tha = tfn ( h, a )
      return
    end if

    ah = absa * h
    gh = alnorm ( h, .false. )
    gah = alnorm ( ah, .false. )
    tha = 0.5D+00 * ( gh + gah ) - gh * gah &
    - tfn ( ah, 1.0D+00 / absa )

    if ( a < 0.0D+00 ) then
      tha = - tha
    end if

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
