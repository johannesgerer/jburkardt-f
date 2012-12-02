function algdiv ( a, b )

!*****************************************************************************80
!
!! ALGDIV computes ln(gamma(b)/gamma(a+b)) when 8 <= B.
!
!  Discussion:
!
!    In this algorithm, del(x) is the function defined by
!    ln(gamma(x)) = (x - 0.5)*ln(x) - x + 0.5*ln(2*pi) + del(x).
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real algdiv
  real alnrel
  real b
  real c
  real, parameter :: c0 =  0.833333333333333E-01
  real, parameter :: c1 = -0.277777777760991E-02
  real, parameter :: c2 =  0.793650666825390E-03
  real, parameter :: c3 = -0.595202931351870E-03
  real, parameter :: c4 =  0.837308034031215E-03
  real, parameter :: c5 = -0.165322962780713E-02
  real d
  real h
  real s11
  real s3
  real s5
  real s7
  real s9
  real t
  real u
  real v
  real w
  real x
  real x2

  if ( b < a ) then
    h = b / a
    c = 1.0E+00 / ( 1.0E+00 + h )
    x = h / ( 1.0E+00 + h )
    d = a + ( b - 0.5E+00 )
  else
    h = a / b
    c = h / ( 1.0E+00 + h )
    x = 1.0E+00 / ( 1.0E+00 + h )
    d = b + ( a - 0.5E+00 )
  end if
!
!  Set sn = (1 - x**n) / ( 1 - x )
!
  x2 = x * x
  s3 = 1.0E+00 + ( x + x2)
  s5 = 1.0E+00 + ( x + x2 * s3 )
  s7 = 1.0E+00 + ( x + x2 * s5 )
  s9 = 1.0E+00 + ( x + x2 * s7 )
  s11 = 1.0E+00 + ( x + x2 * s9 )
!
!  Set w = del(b) - del(a + b)
!
  t = ( 1.0E+00 / b )**2

  w = (((( c5 * s11 &
    * t + c4 * s9 ) &
    * t + c3 * s7 ) &
    * t + c2 * s5 ) &
    * t + c1 * s3 ) &
    * t + c0

  w = w * ( c / b )
!
!  Combine the results.
!
  u = d * alnrel ( a / b )
  v = a * ( log ( b ) - 1.0E+00 )

  if ( v < u ) then
    algdiv = ( w - v ) - u
  else
    algdiv = ( w - u ) - v
  end if

  return
end
function alnrel ( a )

!*****************************************************************************80
!
!! ALNREL evaluates the function ln(1 + a).
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real alnrel
  real, parameter :: p1 = -0.129418923021993E+01
  real, parameter :: p2 =  0.405303492862024E+00
  real, parameter :: p3 = -0.178874546012214E-01
  real, parameter :: q1 = -0.162752256355323E+01
  real, parameter :: q2 =  0.747811014037616E+00
  real, parameter :: q3 = -0.845104217945565E-01
  real t
  real t2
  real w
  real x

  if ( abs ( a ) <= 0.375E+00 ) then

    t = a / ( a + 2.0E+00 )
    t2 = t * t

    w = ((( p3 * t2 + p2 ) * t2 + p1 ) * t2 + 1.0E+00 ) / &
        ((( q3 * t2 + q2 ) * t2 + q1 ) * t2 + 1.0E+00 )

    alnrel = 2.0E+00 * t * w

  else

    x = 1.0E+00 + dble ( a )
    alnrel = log ( x )

  end if

  return
end
function apser ( a, b, x, eps )

!*****************************************************************************80
!
!! APSER yields the incomplete beta ratio i(sub(1-x))(b,a) for
!     a <= min(eps,eps*b), b*x <= 1, and x <= 0.5. used when
!     a is very small. use only if above inequalities are satisfied.
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real aj
  real apser
  real b
  real bx
  real c
  real eps
  real, parameter :: g = 0.577215664901533
  real j
  real psi
  real s
  real t
  real tol
  real x

  bx = b * x
  t = x - bx

  if ( b * eps <= 2.0E-02 ) then
    c = log ( x) + psi ( b ) + g + t
  else
    c = log ( bx ) + g + t
  end if

  tol = 5.0E+00 * eps * abs ( c )
  j = 1.0E+00
  s = 0.0E+00

  do

    j = j + 1.0E+00
    t = t * ( x - bx / j )
    aj = t / j
    s = s + aj

    if ( abs ( aj ) <= tol ) then
      exit
    end if

  end do

  apser = - a * ( c + s )

  return
end
function basym ( a, b, lambda, eps )

!*****************************************************************************80
!
!! BASYM uses an asymptotic expansion for Ix(A,B) for large A and B.
!
!  Discussion:
!
!     lambda = (a + b)*y - b  and eps is the tolerance used.
!     it is assumed that lambda is nonnegative and that
!     a and b are greater than or equal to 15.
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real a0(21)
  real b
  real b0(21)
  real bcorr
  real basym
  real bsum
  real c(21)
  real d(21)
  real dsum
  real, parameter :: e0 = 1.12837916709551E+00
  real, parameter :: e1 = 0.353553390593274E+00
  real eps
  real erfc1
  real f
  real h
  real h2
  real hn
  integer i
  integer im1
  integer imj
  integer j
  real j0
  real j1
  real lambda
  integer m
  integer mm1
  integer mmj
  integer n
  integer np1
  integer, parameter :: num = 20
  real r
  real r0
  real r1
  real rlog1
  real s
  real sum2
  real t
  real t0
  real t1
  real u
  real w
  real w0
  real z
  real z0
  real z2
  real zn
  real znm1
!
!  num is the maximum value that n can take in the do loop
!  ending at statement 50. it is required that num be even.
!  the arrays a0, b0, c, d have dimension num + 1.
!
!     e0 = 2/sqrt(pi)
!     e1 = 2**(-3/2)

  basym = 0.0E+00

  if ( a < b ) then
    h = a/b
    r0 = 1.0 / ( 1.0 + h )
    r1 = ( b - a ) / b
    w0 = 1.0 / sqrt ( a * ( 1.0 + h ) )
  else
    h = b / a
    r0 = 1.0 / ( 1.0 + h )
    r1 = ( b - a ) / a
    w0 = 1.0 / sqrt ( b * ( 1.0 + h ) )
  end if

  f = a * rlog1 ( - lambda / a ) + b * rlog1 ( lambda / b )
  t = exp ( - f )

  if ( t == 0.0 ) then
    return
  end if

  z0 = sqrt ( f )
  z = 0.5 * ( z0 / e1 )
  z2 = f + f

  a0(1) = ( 2.0 / 3.0 ) * r1
  c(1) = - 0.5 * a0(1)
  d(1) = - c(1)
  j0 = ( 0.5 / e0 ) * erfc1 ( 1, z0 )
  j1 = e1
  sum2 = j0 + d(1) * w0 * j1

  s = 1.0
  h2 = h * h
  hn = 1.0
  w = w0
  znm1 = z
  zn = z2

  do n = 2, num, 2

    hn = h2*hn
    a0(n) = 2.0*r0*(1.0 + h*hn)/(n + 2.0)
    np1 = n + 1
    s = s + hn
    a0(np1) = 2.0*r1*s/(n + 3.0)

    do i = n, np1

      r = -0.5*(i + 1.0)
      b0(1) = r*a0(1)

      do m = 2, i

        bsum = 0.0
        mm1 = m - 1

        do j = 1, mm1
          mmj = m - j
          bsum = bsum + (j*r - mmj)*a0(j)*b0(mmj)
        end do

        b0(m) = r*a0(m) + bsum/m

      end do

      c(i) = b0(i)/(i + 1.0)

      dsum = 0.0
      im1 = i - 1

      do j = 1, im1
        imj = i - j
        dsum = dsum + d(imj)*c(j)
      end do

      d(i) = -(dsum + c(i))

    end do

    j0 = e1*znm1 + (n - 1.0)*j0
    j1 = e1*zn + n*j1
    znm1 = z2*znm1
    zn = z2*zn
    w = w0*w
    t0 = d(n)*w*j0
    w = w0*w
    t1 = d(np1)*w*j1
    sum2 = sum2 + (t0 + t1)

    if ((abs(t0) + abs(t1)) <= eps * sum2 ) then
      exit
    end if

  end do

  u = exp ( - bcorr ( a, b ) )
  basym = e0 * t * u * sum2

  return
end
function bcorr ( a0, b0 )

!*****************************************************************************80
!
!! BCORR evaluates  del(a0) + del(b0) - del(a0 + b0)  where
!     ln(gamma(a)) = (a - 0.5)*ln(a) - a + 0.5*ln(2*pi) + del(a).
!     it is assumed that a0 .ge. 8 and b0 .ge. 8.
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real a0
  real b
  real b0
  real bcorr
  real c
  real, parameter :: c0 =  0.833333333333333e-01
  real, parameter :: c1 = -0.277777777760991e-02
  real, parameter :: c2 =  0.793650666825390e-03
  real, parameter :: c3 = -0.595202931351870e-03
  real, parameter :: c4 =  0.837308034031215e-03
  real, parameter :: c5 = -0.165322962780713e-02
  real h
  real s11
  real s3
  real s5
  real s7
  real s9
  real t
  real w
  real x
  real x2

  a = min ( a0, b0 )
  b = max ( a0, b0 )
  h = a/b
  c = h/(1.0 + h)
  x = 1.0/(1.0 + h)
  x2 = x*x
!
!  Set sn = (1 - x**n)/(1 - x)
!
  s3 = 1.0 + (x + x2)
  s5 = 1.0 + (x + x2*s3)
  s7 = 1.0 + (x + x2*s5)
  s9 = 1.0 + (x + x2*s7)
  s11 = 1.0 + (x + x2*s9)
!
!  Set w = del(b) - del(a + b)
!
  t = (1.0/b)**2
  w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3)*t + c0
  w = w*(c/b)
!
!  Compute  del(a) + w
!
  t = (1.0/a)**2
  bcorr = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0)/a + w

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
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964.
!
!    Karl Pearson,
!    Tables of the Incomplete Beta Function,
!    Cambridge University Press, 1968.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real A, B, the parameters of the function.
!
!    Output, real X, the argument of the function.
!
!    Output, real FX, the value of the function.
!
  implicit none

  integer, parameter :: n_max = 42

  real a
  real, save, dimension ( n_max ) :: a_vec = (/ &
     0.5E+00, &
     0.5E+00, &
     0.5E+00, &
     1.0E+00, &
     1.0E+00, &
     1.0E+00, &
     1.0E+00, &
     1.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     5.5E+00, &
    10.0E+00, &
    10.0E+00, &
    10.0E+00, &
    10.0E+00, &
    20.0E+00, &
    20.0E+00, &
    20.0E+00, &
    20.0E+00, &
    20.0E+00, &
    30.0E+00, &
    30.0E+00, &
    40.0E+00, &
     0.1E+01, &
     0.1E+01, &
     0.1E+01, &
     0.1E+01, &
     0.1E+01, &
     0.1E+01, &
     0.1E+01, &
     0.1E+01, &
     0.2E+01, &
     0.3E+01, &
     0.4E+01, &
     0.5E+01 /)
  real b
  real, save, dimension ( n_max ) :: b_vec = (/ &
     0.5E+00, &
     0.5E+00, &
     0.5E+00, &
     0.5E+00, &
     0.5E+00, &
     0.5E+00, &
     0.5E+00, &
     1.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     2.0E+00, &
     5.0E+00, &
     0.5E+00, &
     5.0E+00, &
     5.0E+00, &
    10.0E+00, &
     5.0E+00, &
    10.0E+00, &
    10.0E+00, &
    20.0E+00, &
    20.0E+00, &
    10.0E+00, &
    10.0E+00, &
    20.0E+00, &
     0.5E+00, &
     0.5E+00, &
     0.5E+00, &
     0.5E+00, &
     0.2E+01, &
     0.3E+01, &
     0.4E+01, &
     0.5E+01, &
     0.2E+01, &
     0.2E+01, &
     0.2E+01, &
     0.2E+01 /)
  real fx
  real, save, dimension ( n_max ) :: fx_vec = (/ &
    0.6376856085851985E-01, &
    0.2048327646991335E+00, &
    0.1000000000000000E+01, &
    0.0000000000000000E+00, &
    0.5012562893380045E-02, &
    0.5131670194948620E-01, &
    0.2928932188134525E+00, &
    0.5000000000000000E+00, &
    0.2800000000000000E-01, &
    0.1040000000000000E+00, &
    0.2160000000000000E+00, &
    0.3520000000000000E+00, &
    0.5000000000000000E+00, &
    0.6480000000000000E+00, &
    0.7840000000000000E+00, &
    0.8960000000000000E+00, &
    0.9720000000000000E+00, &
    0.4361908850559777E+00, &
    0.1516409096347099E+00, &
    0.8978271484375000E-01, &
    0.1000000000000000E+01, &
    0.5000000000000000E+00, &
    0.4598773297575791E+00, &
    0.2146816102371739E+00, &
    0.9507364826957875E+00, &
    0.5000000000000000E+00, &
    0.8979413687105918E+00, &
    0.2241297491808366E+00, &
    0.7586405487192086E+00, &
    0.7001783247477069E+00, &
    0.5131670194948620E-01, &
    0.1055728090000841E+00, &
    0.1633399734659245E+00, &
    0.2254033307585166E+00, &
    0.3600000000000000E+00, &
    0.4880000000000000E+00, &
    0.5904000000000000E+00, &
    0.6723200000000000E+00, &
    0.2160000000000000E+00, &
    0.8370000000000000E-01, &
    0.3078000000000000E-01, &
    0.1093500000000000E-01 /)
  integer n_data
  real x
  real, save, dimension ( n_max ) :: x_vec = (/ &
    0.01E+00, &
    0.10E+00, &
    1.00E+00, &
    0.00E+00, &
    0.01E+00, &
    0.10E+00, &
    0.50E+00, &
    0.50E+00, &
    0.10E+00, &
    0.20E+00, &
    0.30E+00, &
    0.40E+00, &
    0.50E+00, &
    0.60E+00, &
    0.70E+00, &
    0.80E+00, &
    0.90E+00, &
    0.50E+00, &
    0.90E+00, &
    0.50E+00, &
    1.00E+00, &
    0.50E+00, &
    0.80E+00, &
    0.60E+00, &
    0.80E+00, &
    0.50E+00, &
    0.60E+00, &
    0.70E+00, &
    0.80E+00, &
    0.70E+00, &
    0.10E+00, &
    0.20E+00, &
    0.30E+00, &
    0.40E+00, &
    0.20E+00, &
    0.20E+00, &
    0.20E+00, &
    0.20E+00, &
    0.30E+00, &
    0.30E+00, &
    0.30E+00, &
    0.30E+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0E+00
    b = 0.0E+00
    x = 0.0E+00
    fx = 0.0E+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine beta_log_values ( n, x, y, fxy )

!*****************************************************************************80
!
!! BETA_LOG_VALUES returns some values of the Beta function for testing.
!
!  Modified:
!
!    15 February 2004
!
!  Reference:
!
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, real X, Y, the arguments of the function.
!
!    Output, real FXY, the value of the function.
!
  implicit none

  integer, parameter :: nmax = 17

  real, save, dimension ( nmax ) :: fxvec = (/ &
    1.609437912E+00,  0.9162907319E+00,  0.5108256238E+00,  0.2231435513E+00, &
    1.609437912E+00,  0.9162907319E+00,  0.0000000000E+00, -1.791759469E+00, &
   -3.401197382E+00, -4.941642423E+00,  -6.445719819E+00,  -3.737669618E+00, &
   -5.123963979E+00, -6.222576268E+00,  -7.138867000E+00,  -7.927324360E+00, &
   -9.393661429E+00 /)
  real fxy
  integer n
  real x
  real, save, dimension ( nmax ) :: xvec = (/ &
    0.2E+00, 0.4E+00, 0.6E+00, 0.8E+00, &
    1.0E+00, 1.0E+00, 1.0E+00, 2.0E+00, &
    3.0E+00, 4.0E+00, 5.0E+00, 6.0E+00, &
    6.0E+00, 6.0E+00, 6.0E+00, 6.0E+00, &
    7.0E+00 /)
  real y
  real, save, dimension ( nmax ) :: yvec = (/ &
    1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, &
    0.2E+00, 0.4E+00, 1.0E+00, 2.0E+00, &
    3.0E+00, 4.0E+00, 5.0E+00, 2.0E+00, &
    3.0E+00, 4.0E+00, 5.0E+00, 6.0E+00, &
    7.0E+00 /)

  if ( n < 0 ) then
    n = 0
  end if

  n = n + 1

  if ( nmax < n ) then
    n = 0
    x = 0.0E+00
    y = 0.0E+00
    fxy = 0.0E+00
    return
  end if

  x = xvec(n)
  y = yvec(n)
  fxy = fxvec(n)

  return
end
function betaln ( a0, b0 )

!*****************************************************************************80
!
!! BETALN evaluates the logarithm of the Beta function.
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
!  Local Parameters:
!
!    Local, real E, the value of Log ( 2 * PI ) / 2.
!
  implicit none

  real a
  real a0
  real algdiv
  real alnrel
  real b
  real b0
  real bcorr
  real betaln
  real c
  real, parameter :: e = 0.918938533204673
  real gamln
  real gsumln
  real h
  integer i
  integer n
  real u
  real v
  real w
  real z

  a = min ( a0, b0 )
  b = max ( a0, b0 )
  if (a .ge. 8.0) go to 60
  if (a .ge. 1.0) go to 20
!
!  a < 1
!
  if ( b < 8.0 ) then
    betaln = gamln ( a ) + ( gamln ( b ) - gamln ( a + b ) )
    return
  else
    betaln = gamln ( a ) + algdiv ( a, b )
    return
  end if
!
!  procedure when 1 <= a < 8
!
   20 if (a .gt. 2.0) go to 30
      if (b .gt. 2.0) go to 21
         betaln = gamln(a) + gamln(b) - gsumln(a,b)
         return
   21 w = 0.0
      if (b < 8.0) go to 40
         betaln = gamln(a) + algdiv(a,b)
         return
!
! reduction of a when b <= 1000
!
   30 if (b .gt. 1000.0) go to 50
      n = a - 1.0
      w = 1.0
      do i = 1,n
         a = a - 1.0
         h = a/b
         w = w * (h/(1.0 + h))
      end do
      w = log ( w )
      if (b < 8.0) go to 40
      betaln = w + gamln(a) + algdiv(a,b)
      return
!
!  Reduction of b when b < 8
!
   40 n = b - 1.0
      z = 1.0
      do i = 1,n
         b = b - 1.0
         z = z * (b/(a + b))
      end do
      betaln = w + log ( z ) + (gamln(a) + (gamln(b) - gsumln(a,b)))
      return
!
!  reduction of a when b .gt. 1000
!
   50 n = a - 1.0
      w = 1.0
      do i = 1,n
         a = a - 1.0
         w = w * (a/(1.0 + a/b))
      end do
      betaln = ( log ( w ) - n* log ( b ) ) + (gamln(a) + algdiv(a,b))
      return
!
!  procedure when a .ge. 8
!
   60 w = bcorr(a,b)
      h = a/b
      c = h/(1.0 + h)
      u = -(a - 0.5) * log ( c )
      v = b*alnrel(h)
      if (u <= v) go to 61
         betaln = (((-0.5* log ( b ) + e) + w) - v) - u
         return
   61 betaln = (((-0.5* log ( b ) + e) + w) - u) - v

  return
end
function bfrac ( a, b, x, y, lambda, eps )

!*****************************************************************************80
!
!! BFRAC uses a continued fraction expansion for ix(a,b) when a,b .gt. 1.
!
!     it is assumed that  lambda = (a + b)*y - b.
!
!  Modified:
!
!    28 August 2004
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real alpha
  real an
  real anp1
  real b
  real beta
  real bfrac
  real bn
  real bnp1
  real brcomp
  real c
  real c0
  real c1
  real e
  real eps
  real lambda
  real n
  real p
  real r
  real r0
  real s
  real t
  real w
  real x
  real y
  real yp1

  bfrac = brcomp ( a, b, x, y )

  if ( bfrac == 0.0 ) then
    return
  end if

  c = 1.0 + lambda
  c0 = b / a
  c1 = 1.0 + 1.0 / a
  yp1 = y + 1.0

  n = 0.0
  p = 1.0
  s = a + 1.0
  an = 0.0
  bn = 1.0
  anp1 = 1.0
  bnp1 = c / c1
  r = c1 / c
!
!  Continued fraction calculation.
!
  do

    n = n + 1.0
    t = n / a
    w = n * ( b - n ) * x
    e = a / s
    alpha = ( p * ( p + c0 ) * e * e ) * ( w * x )
    e = ( 1.0 + t ) / ( c1 + t + t )
    beta = n + w / s + e * ( c + n * yp1 )
    p = 1.0 + t
    s = s + 2.0
!
!  Update AN, BN, ANP1, and BNP1.
!
    t = alpha * an + beta * anp1
    an = anp1
    anp1 = t
    t = alpha * bn + beta * bnp1
    bn = bnp1
    bnp1 = t
    r0 = r
    r = anp1 / bnp1

    if ( abs ( r - r0 ) <= eps * r ) then
      exit
    end if
!
!  Rescale AN, BN, ANP1, and BNP1.
!
    an = an / bnp1
    bn = bn / bnp1
    anp1 = r
    bnp1 = 1.0

  end do
!
!  Termination.
!
  bfrac = bfrac * r

  return
end
subroutine bgrat ( a, b, x, y, w, eps, ierr )

!*****************************************************************************80
!
!! BGRAT uses an asymptotic expansion for Ix(a,b) when a is larger than b.
!
!  Discussion:
!
!    The result of the expansion is added to w. it is assumed
!    that a .ge. 15 and b <= 1.  eps is the tolerance used.
!    ierr is a variable that reports the status of the results.
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real algdiv
  real alnrel
  real b
  real bm1
  real bp2n
  real c(30)
  real cn
  real coef
  real d(30)
  real dj
  real eps
  real gam1
  integer i
  integer ierr
  real j
  real l
  real lnx
  integer n
  real n2
  integer nm1
  real nu
  real p
  real q
  real r
  real s
  real sum1
  real t
  real t2
  real u
  real v
  real w
  real x
  real y
  real z

  bm1 = (b - 0.5) - 0.5
  nu = a + 0.5*bm1
  if (y .gt. 0.375) go to 10
  lnx = alnrel(-y)
  go to 11
10 continue
  lnx = log ( x )
11 continue
  z = -nu*lnx
  if (b*z == 0.0) go to 100
!
!  computation of the expansion
!  set r = exp(-z)*z**b/gamma(b)
!
      r = b*(1.0 + gam1(b))*exp ( b * log ( z ) )
      r = r*exp(a*lnx)*exp(0.5*bm1*lnx)
      u = algdiv(b,a) + b* log ( nu )
      u = r*exp(-u)
      if (u == 0.0) go to 100
      call grat1(b,z,r,p,q,eps)

      v = 0.25*(1.0/nu)**2
      t2 = 0.25*lnx*lnx
      l = w/u
      j = q/r
      sum1 = j
      t = 1.0
      cn = 1.0
      n2 = 0.0
      do n = 1,30
         bp2n = b + n2
         j = (bp2n*(bp2n + 1.0)*j + (z + bp2n + 1.0)*t)*v
         n2 = n2 + 2.0
         t = t*t2
         cn = cn/(n2*(n2 + 1.0))
         c(n) = cn
         s = 0.0
         if (n == 1) go to 21
            nm1 = n - 1
            coef = b - n
            do i = 1,nm1
               s = s + coef*c(i)*d(n-i)
               coef = coef + b
            end do
   21    d(n) = bm1*cn + s/n
         dj = d(n)*j
         sum1 = sum1 + dj

         if ( sum1 <= 0.0 ) then
           go to 100
         end if

         if (abs(dj) <= eps*(sum1 + l)) then
           go to 30
         end if

      end do
!
!  Add the results to w
!
   30 ierr = 0
      w = w + u * sum1
      return
!
!  the expansion cannot be computed
!
  100 ierr = 1

  return
end
function bpser ( a, b, x, eps )

!*****************************************************************************80
!
!! BPSER uses the power series expansion for evaluating ix(a,b) when b <= 1
!     or b*x <= 0.7.  eps is the tolerance used.
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real a0
  real algdiv
  real apb
  real b
  real b0
  real betaln
  real bpser
  real c
  real eps
  real gam1
  real gamln1
  integer i
  integer m
  real n
  real sum1
  real t
  real tol
  real u
  real w
  real x
  real z

  bpser = 0.0
  if ( x == 0.0 ) then
    return
  end if
!
!  compute the factor x**a/(a*beta(a,b))
!
      a0 = min ( a, b )
      if (a0 < 1.0) go to 10
         z = a* log ( x ) - betaln(a,b)
         bpser = exp(z)/a
         go to 70
   10 b0 = max ( a, b )
      if (b0 .ge. 8.0) go to 60
      if (b0 .gt. 1.0) go to 40
!
!  procedure for a0 < 1 and b0 <= 1
!
      bpser = x**a
      if (bpser == 0.0) return

      apb = a + b
      if (apb .gt. 1.0) go to 20
         z = 1.0 + gam1(apb)
         go to 30
   20 u = dble(a) + dble(b) - 1.d0
      z = (1.0 + gam1(u))/apb

   30 c = (1.0 + gam1(a))*(1.0 + gam1(b))/z
      bpser = bpser*c*(b/apb)
      go to 70
!
!  procedure for a0 < 1 and 1 < b0 < 8
!
   40 u = gamln1(a0)
      m = b0 - 1.0
      if (m < 1) go to 50
      c = 1.0
      do i = 1,m
         b0 = b0 - 1.0
         c = c*(b0/(a0 + b0))
      end do
      u = log ( c ) + u

   50 z = a* log ( x ) - u
      b0 = b0 - 1.0
      apb = a0 + b0
      if (apb .gt. 1.0) go to 51
         t = 1.0 + gam1(apb)
         go to 52
   51 u = dble(a0) + dble(b0) - 1.d0
      t = (1.0 + gam1(u))/apb
   52 bpser = exp(z)*(a0/a)*(1.0 + gam1(b0))/t
      go to 70
!
!  procedure for a0 < 1 and b0 .ge. 8
!
   60 u = gamln1(a0) + algdiv(a0,b0)
      z = a* log ( x ) - u
      bpser = (a0/a)*exp(z)
   70 if (bpser == 0.0 .or. a <= 0.1*eps) return
!
!  Compute the series
!
  sum1 = 0.0
  n = 0.0
  c = 1.0
  tol = eps/a

  do
    n = n + 1.0
    c = c*(0.5 + (0.5 - b/n))*x
    w = c/(a + n)
    sum1 = sum1 + w
    if ( abs ( w ) <= tol ) then
      exit
    end if
  end do

  bpser = bpser*(1.0 + a*sum1)

  return
end
subroutine bratio ( a, b, x, y, w, w1, ierr )

!*****************************************************************************80
!
!! BRATIO evaluates the incomplete beta function Ix(A,B).
!
!  Discussion:
!
!    It is assumed that X <= 1
!    and Y = 1 - X.  BRATIO assigns W and W1 the values
!
!                      W  = ix(a,b)
!                      W1 = 1 - ix(a,b)
!
!     ierr is a variable that reports the status of the results.
!     if no input errors are detected then ierr is set to 0 and
!     w and w1 are computed. otherwise, if an error is detected,
!     then w and w1 are assigned the value 0 and ierr is set to
!     one of the following values ...
!
!        ierr = 1  if a or b is negative
!        ierr = 2  if a = b = 0
!        ierr = 3  if x < 0 or x .gt. 1
!        ierr = 4  if y < 0 or y .gt. 1
!        ierr = 5  if x + y /= 1
!        ierr = 6  if x = a = 0
!        ierr = 7  if y = b = 0
!
!  Author:
!
!    Alfred Morris
!    Naval Surface Warfare Center
!    Dahlgren, Virginia
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
!  Parameters:
!
!    Input, real A, B, the parameters of the function.  A and B should
!    be nonnegative.
!
  implicit none

  real a
  real a0
  real apser
  real b
  real b0
  real basym
  real bfrac
  real bpser
  real bup
  real eps
  real fpser
  integer ierr
  integer ierr1
  integer ind
  real lambda
  integer n
  real t
  real w
  real w1
  real x
  real x0
  real y
  real y0
  real z

  eps = epsilon ( eps )

  w = 0.0
  w1 = 0.0
  if (a < 0.0 .or. b < 0.0) go to 300
  if (a == 0.0 .and. b == 0.0) go to 310
  if (x < 0.0 .or. x .gt. 1.0) go to 320
      if (y < 0.0 .or. y .gt. 1.0) go to 330

  z = ((x + y) - 0.5) - 0.5

  if ( abs(z) .gt. 3.0*eps) then
    ierr = 5
    return
  end if

  ierr = 0

  if (x == 0.0) go to 200
  if (y == 0.0) go to 210
  if (a == 0.0) go to 211
  if (b == 0.0) go to 201

      eps = max ( eps, 1.e-15 )
      if ( max ( a, b ) < 1.e-3*eps) go to 230

      ind = 0
      a0 = a
      b0 = b
      x0 = x
      y0 = y
      if ( min ( a0, b0 ) .gt. 1.0) go to 30
!
!  procedure for a0 <= 1 or b0 <= 1
!
      if (x <= 0.5) go to 10
      ind = 1
      a0 = b
      b0 = a
      x0 = y
      y0 = x

   10 if (b0 < min ( eps, eps * a0 ) ) go to 80
      if (a0 < min ( eps, eps * b0 ) .and. b0*x0 <= 1.0) go to 90
      if ( max ( a0, b0 ) .gt. 1.0) go to 20
      if (a0 .ge. min ( 0.2, b0 ) ) go to 100
      if (x0**a0 <= 0.9) go to 100
      if (x0 .ge. 0.3) go to 110
      n = 20
      go to 130

   20 if (b0 <= 1.0) go to 100
      if (x0 .ge. 0.3) go to 110
      if (x0 .ge. 0.1) go to 21
      if ((x0*b0)**a0 <= 0.7) go to 100
   21 if (b0 .gt. 15.0) go to 131
      n = 20
      go to 130
!
!  procedure for a0 .gt. 1 and b0 .gt. 1
!
   30 if (a .gt. b) go to 31
         lambda = a - (a + b)*x
         go to 32
   31 lambda = (a + b)*y - b

   32 if (lambda .ge. 0.0) go to 40
      ind = 1
      a0 = b
      b0 = a
      x0 = y
      y0 = x
      lambda = abs(lambda)

   40 if (b0 < 40.0 .and. b0*x0 <= 0.7) go to 100
      if (b0 < 40.0) go to 140
      if (a0 .gt. b0) go to 50
         if (a0 <= 100.0) go to 120
         if (lambda .gt. 0.03*a0) go to 120
         go to 180
   50 if (b0 <= 100.0) go to 120
      if (lambda .gt. 0.03*b0) go to 120
      go to 180
!
!  evaluation of the appropriate algorithm
!
   80 continue
      w = fpser(a0, b0, x0, eps)
      w1 = 0.5 + (0.5 - w)
      go to 220

   90 continue
      w1 = apser(a0, b0, x0, eps)
      w = 0.5 + (0.5 - w1)
      go to 220

  100 continue
      w = bpser(a0, b0, x0, eps)
      w1 = 0.5 + (0.5 - w)
      go to 220

  110 continue
      w1 = bpser(b0, a0, y0, eps)
      w = 0.5 + (0.5 - w1)
      go to 220

  120 continue
      w = bfrac(a0, b0, x0, y0, lambda, 15.0*eps)
      w1 = 0.5 + (0.5 - w)
      go to 220

  130 continue
      w1 = bup(b0, a0, y0, x0, n, eps)
      b0 = b0 + n
  131 call bgrat(b0, a0, y0, x0, w1, 15.0*eps, ierr1)
      w = 0.5 + (0.5 - w1)
      go to 220

  140 n = b0
      b0 = b0 - n

      if (b0 == 0.0) then
         n = n - 1
         b0 = 1.0
      end if

  141 w = bup(b0, a0, y0, x0, n, eps)
      if (x0 .gt. 0.7) go to 150
      w = w + bpser(a0, b0, x0, eps)
      w1 = 0.5 + (0.5 - w)
      go to 220

  150 if ( 15.0 < a0 ) go to 151
         n = 20
         w = w + bup(a0, b0, x0, y0, n, eps)
         a0 = a0 + n
  151 call bgrat(a0, b0, x0, y0, w, 15.0*eps, ierr1)
      w1 = 0.5 + (0.5 - w)
      go to 220

  180 w = basym(a0, b0, lambda, 100.0*eps)
      w1 = 0.5 + (0.5 - w)
      go to 220
!
!  termination of the procedure
!
  200 if (a == 0.0) go to 350
  201 w = 0.0
      w1 = 1.0
      return

  210 if (b == 0.0) go to 360
  211 w = 1.0
      w1 = 0.0
      return

  220 if (ind == 0) return
      t = w
      w = w1
      w1 = t
      return
!
!  procedure for a and b < 1.e-3*eps
!
  230 w = b/(a + b)
      w1 = a/(a + b)
      return
!
!  Error return
!
  300 ierr = 1
      return
  310 ierr = 2
      return
  320 ierr = 3
      return
  330 ierr = 4
      return

  350 ierr = 6
      return
  360 ierr = 7

  return
end
function brcmp1 ( mu, a, b, x, y )

!*****************************************************************************80
!
!! BRCMP1 evaluates exp(mu) * (x**a*y**b/beta(a,b)).
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real a0
  real algdiv
  real alnrel
  real apb
  real b
  real b0
  real bcorr
  real betaln
  real brcmp1
  real c
  real, parameter :: const = 0.398942280401433
  real e
  real esum
  real gam1
  real gamln1
  real h
  integer i
  real lambda
  real lnx
  real lny
  integer mu
  integer n
  real rlog1
  real t
  real u
  real v
  real x
  real x0
  real y
  real y0
  real z

  a0 = min ( a, b )
      if (a0 .ge. 8.0) go to 100

      if (x .gt. 0.375) go to 10
         lnx = log ( x )
         lny = alnrel(-x)
         go to 20
   10 if (y .gt. 0.375) go to 11
         lnx = alnrel(-y)
         lny = log ( y )
         go to 20
   11 lnx = log ( x )
      lny = log ( y )

   20 z = a*lnx + b*lny
      if (a0 < 1.0) go to 30
      z = z - betaln(a,b)
      brcmp1 = esum(mu,z)
      return
!-----------------------------------------------------------------------
!  procedure for a < 1 or b < 1
!-----------------------------------------------------------------------
   30 b0 = max ( a, b )
      if (b0 .ge. 8.0) go to 80
      if (b0 .gt. 1.0) go to 60
!
!  algorithm for b0 <= 1
!
      brcmp1 = esum(mu,z)
      if (brcmp1 == 0.0) return

      apb = a + b
      if ( 1.0 < apb ) go to 40
         z = 1.0 + gam1(apb)
         go to 50
   40 u = dble(a) + dble(b) - 1.d0
      z = (1.0 + gam1(u))/apb

   50 c = (1.0 + gam1(a))*(1.0 + gam1(b))/z
      brcmp1 = brcmp1*(a0*c)/(1.0 + a0/b0)
      return
!
!  algorithm for 1 < b0 < 8
!
   60 u = gamln1(a0)
      n = b0 - 1.0
      if (n < 1) go to 70
      c = 1.0
      do 61 i = 1,n
         b0 = b0 - 1.0
         c = c*(b0/(a0 + b0))
   61 continue
      u = log ( c ) + u

   70 z = z - u
      b0 = b0 - 1.0
      apb = a0 + b0
      if (apb .gt. 1.0) go to 71
         t = 1.0 + gam1(apb)
         go to 72
   71 u = dble(a0) + dble(b0) - 1.d0
      t = (1.0 + gam1(u))/apb
   72 brcmp1 = a0*esum(mu,z)*(1.0 + gam1(b0))/t
      return
!
!  algorithm for b0 .ge. 8
!
   80 u = gamln1(a0) + algdiv(a0,b0)
      brcmp1 = a0*esum(mu,z - u)
      return
!-----------------------------------------------------------------------
!  procedure for a .ge. 8 and b .ge. 8
!-----------------------------------------------------------------------
  100 if (a .gt. b) go to 101
         h = a/b
         x0 = h/(1.0 + h)
         y0 = 1.0/(1.0 + h)
         lambda = a - (a + b)*x
         go to 110
  101 h = b/a
      x0 = 1.0/(1.0 + h)
      y0 = h/(1.0 + h)
      lambda = (a + b)*y - b

  110 e = -lambda/a
      if (abs(e) .gt. 0.6) go to 111
         u = rlog1(e)
         go to 120
  111 u = e - log ( x / x0 )

  120 e = lambda/b
      if (abs(e) .gt. 0.6) go to 121
         v = rlog1(e)
         go to 130
  121 v = e - log ( y / y0 )

  130 z = esum(mu,-(a*u + b*v))
      brcmp1 = const*sqrt(b*x0)*z*exp(-bcorr(a,b))

  return
end
function brcomp ( a, b, x, y )

!*****************************************************************************80
!
!! BRCOMP evaluates X**a * y**b / beta(a,b).
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real a0
  real algdiv
  real alnrel
  real apb
  real b
  real b0
  real bcorr
  real betaln
  real brcomp
  real c
  real, parameter :: const = 0.398942280401433
  real e
  real gam1
  real gamln1
  real h
  integer i
  real lambda
  real lnx
  real lny
  integer n
  real rlog1
  real t
  real u
  real v
  real x
  real x0
  real y
  real y0
  real z

  brcomp = 0.0

  if ( x == 0.0 .or. y == 0.0 ) then
    return
  end if

  a0 = min ( a, b )

      if ( 8.0 <= a0 ) go to 100

      if (x .gt. 0.375) go to 10
         lnx = log ( x )
         lny = alnrel(-x)
         go to 20
   10 if (y .gt. 0.375) go to 11
         lnx = alnrel(-y)
         lny = log ( y )
         go to 20
   11 lnx = log ( x )
      lny = log ( y )

   20 z = a*lnx + b*lny

      if ( 1.0 <= a ) then
        z = z - betaln(a,b)
        brcomp = exp(z)
        return
      end if
!-----------------------------------------------------------------------
!  procedure for a < 1 or b < 1
!-----------------------------------------------------------------------
   30 b0 = max ( a, b )
      if (b0 .ge. 8.0) go to 80
      if ( 1.0 < b0 ) go to 60
!
!  algorithm for b0 <= 1
!
      brcomp = exp(z)
      if (brcomp == 0.0) return

      apb = a + b
      if (apb .gt. 1.0) go to 40
         z = 1.0 + gam1(apb)
         go to 50
   40 u = dble(a) + dble(b) - 1.d0
      z = (1.0 + gam1(u))/apb

   50 c = (1.0 + gam1(a))*(1.0 + gam1(b))/z
      brcomp = brcomp*(a0*c)/(1.0 + a0/b0)
      return
!
!  algorithm for 1 < b0 < 8
!
   60 u = gamln1(a0)
      n = b0 - 1.0
      if (n < 1) go to 70
      c = 1.0
      do i = 1,n
         b0 = b0 - 1.0
         c = c*(b0/(a0 + b0))
      end do
      u = log ( c ) + u

   70 z = z - u
      b0 = b0 - 1.0
      apb = a0 + b0
      if (apb .gt. 1.0) go to 71
         t = 1.0 + gam1(apb)
         go to 72
   71 u = dble(a0) + dble(b0) - 1.d0
      t = (1.0 + gam1(u))/apb
   72 brcomp = a0*exp(z)*(1.0 + gam1(b0))/t
      return
!
!  algorithm for b0 .ge. 8
!
   80 u = gamln1(a0) + algdiv(a0,b0)
      brcomp = a0*exp(z - u)
      return
!-----------------------------------------------------------------------
!  procedure for a .ge. 8 and b .ge. 8
!-----------------------------------------------------------------------
  100 if (a .gt. b) go to 101
         h = a/b
         x0 = h/(1.0 + h)
         y0 = 1.0/(1.0 + h)
         lambda = a - (a + b)*x
         go to 110
  101 h = b/a
      x0 = 1.0/(1.0 + h)
      y0 = h/(1.0 + h)
      lambda = (a + b)*y - b

  110 e = -lambda/a
      if (abs(e) .gt. 0.6) go to 111
         u = rlog1(e)
         go to 120
  111 u = e - log ( x / x0 )

  120 e = lambda/b
      if (abs(e) .gt. 0.6) go to 121
         v = rlog1(e)
         go to 130
  121 v = e - log ( y / y0 )

  130 z = exp(-(a*u + b*v))
      brcomp = const*sqrt(b*x0)*z*exp(-bcorr(a,b))

  return
end
function bup ( a, b, x, y, n, eps )

!*****************************************************************************80
!
!! BUP evaluates ix(a,b) - ix(a+n,b) where n is a positive integer.
!
!     eps is the tolerance used.
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real ap1
  real apb
  real b
  real brcmp1
  real bup
  real d
  real eps
  real exparg
  integer i
  integer k
  integer kp1
  real l
  integer mu
  integer n
  integer nm1
  real r
  real t
  real w
  real x
  real y
!
!  obtain the scaling factor exp(-mu) and
!  exp(mu)*(x**a*y**b/beta(a,b))/a
!
  apb = a + b
  ap1 = a + 1.0
  mu = 0
  d = 1.0
  if (n == 1 .or. a < 1.0) go to 10
      if (apb < 1.1*ap1) go to 10
         mu = abs(exparg(1))
         k = exparg(0)
         if (k < mu) mu = k
         t = mu
         d = exp(-t)

10 continue

  bup = brcmp1(mu,a,b,x,y)/a
      if (n == 1 .or. bup == 0.0) return
      nm1 = n - 1
      w = d
!
!  let k be the index of the maximum term
!
      k = 0
      if (b <= 1.0) go to 40
      if ( 1.e-4 < y ) go to 20
         k = nm1
         go to 30
   20 r = (b - 1.0)*x/y - a
      if (r < 1.0) go to 40
      k = nm1
      t = nm1
      if (r < t) k = r
!
!  add the increasing terms of the series
!
   30 do i = 1,k
         l = i - 1
         d = ((apb + l)/(ap1 + l))*x*d
         w = w + d
      end do
      if (k == nm1) go to 50
!
!  add the remaining terms of the series
!
   40 kp1 = k + 1
      do i = kp1,nm1
         l = i - 1
         d = ((apb + l)/(ap1 + l))*x*d
         w = w + d
         if (d <= eps*w) go to 50
      end do
!
!  terminate the procedure
!
   50 bup = bup*w

  return
end
function erfc1 ( ind, x )

!*****************************************************************************80
!
!! ERFC1 evaluates the complementary error function
!
!          erfc1(ind,x) = erfc(x)            if ind = 0
!          erfc1(ind,x) = exp(x*x)*erfc(x)   otherwise
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a(5)
  real ax
  real b(3)
  real bot
  real, parameter :: c = 0.564189583547756
  real e
  real erfc1
  real exparg
  integer ind
  real p(8)
  real q(8)
  real r(5)
  real s(4)
  real t
  real top
  double precision w
  real x

      data a(1) /.771058495001320e-04/, a(2)/-.133733772997339e-02/, &
           a(3) /.323076579225834e-01/, a(4) /.479137145607681e-01/, &
           a(5) /.128379167095513e+00/
      data b(1) /.301048631703895e-02/, b(2) /.538971687740286e-01/, &
           b(3) /.375795757275549e+00/

      data p(1)/-1.36864857382717e-07/, p(2) /5.64195517478974e-01/, &
           p(3) /7.21175825088309e+00/, p(4) /4.31622272220567e+01/, &
           p(5) /1.52989285046940e+02/, p(6) /3.39320816734344e+02/, &
           p(7) /4.51918953711873e+02/, p(8) /3.00459261020162e+02/
      data q(1) /1.00000000000000e+00/, q(2) /1.27827273196294e+01/, &
           q(3) /7.70001529352295e+01/, q(4) /2.77585444743988e+02/, &
           q(5) /6.38980264465631e+02/, q(6) /9.31354094850610e+02/, &
           q(7) /7.90950925327898e+02/, q(8) /3.00459260956983e+02/
      data r(1) /2.10144126479064e+00/, r(2) /2.62370141675169e+01/, &
           r(3) /2.13688200555087e+01/, r(4) /4.65807828718470e+00/, &
           r(5) /2.82094791773523e-01/
      data s(1) /9.41537750555460e+01/, s(2) /1.87114811799590e+02/, &
           s(3) /9.90191814623914e+01/, s(4) /1.80124575948747e+01/

!
!  abs(x) <= 0.5
!
      ax = abs(x)
      if (ax .gt. 0.5) go to 10
      t = x*x
      top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0
      bot = ((b(1)*t + b(2))*t + b(3))*t + 1.0
      erfc1 = 0.5 + (0.5 - x*(top/bot))
      if (ind /= 0) erfc1 = exp(t) * erfc1
      return
!
!  0.5 < abs(x) <= 4
!
   10 if (ax .gt. 4.0) go to 20
      top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax &
                          + p(6))*ax + p(7))*ax + p(8)
      bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax &
                          + q(6))*ax + q(7))*ax + q(8)
      erfc1 = top/bot
      go to 40
!
!  abs(x) .gt. 4
!
   20 if (x <= -5.6) go to 50
      if (ind /= 0) go to 30
      if ( 100.0 < x ) go to 60
      if ( -exparg(1) < x * x ) go to 60

   30 t = (1.0/x)**2
      top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
      bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + 1.0
      erfc1 = (c - t*top/bot)/ax
!
!  final assembly
!
   40 if (ind == 0) go to 41
         if (x < 0.0) erfc1 = 2.0*exp(x*x) - erfc1
         return
   41 w = dble(x)*dble(x)
      t = w
      e = w - dble(t)
      erfc1 = ((0.5 + (0.5 - e)) * exp(-t)) * erfc1
      if (x < 0.0) erfc1 = 2.0 - erfc1
      return
!
!  limit value for large negative x
!
   50 erfc1 = 2.0
      if (ind /= 0) erfc1 = 2.0*exp(x*x)
      return
!
!  limit value for large positive x when ind = 0
!
   60 erfc1 = 0.0

  return
end
function erf ( x )

!*****************************************************************************80
!
!! ERF evaluates the real error function.
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a(5)
  real ax
  real b(3)
  real bot
  real, parameter :: c = 0.564189583547756
  real erf
  real p(8)
  real q(8)
  real r(5)
  real s(4)
  real t
  real top
  real x
  real x2

      data a(1) /.771058495001320e-04/, a(2)/-.133733772997339e-02/, &
           a(3) /.323076579225834e-01/, a(4) /.479137145607681e-01/, &
           a(5) /.128379167095513e+00/
      data b(1) /.301048631703895e-02/, b(2) /.538971687740286e-01/, &
           b(3) /.375795757275549e+00/
      data p(1)/-1.36864857382717e-07/, p(2) /5.64195517478974e-01/, &
           p(3) /7.21175825088309e+00/, p(4) /4.31622272220567e+01/, &
           p(5) /1.52989285046940e+02/, p(6) /3.39320816734344e+02/, &
           p(7) /4.51918953711873e+02/, p(8) /3.00459261020162e+02/
      data q(1) /1.00000000000000e+00/, q(2) /1.27827273196294e+01/, &
           q(3) /7.70001529352295e+01/, q(4) /2.77585444743988e+02/, &
           q(5) /6.38980264465631e+02/, q(6) /9.31354094850610e+02/, &
           q(7) /7.90950925327898e+02/, q(8) /3.00459260956983e+02/
      data r(1) /2.10144126479064e+00/, r(2) /2.62370141675169e+01/, &
           r(3) /2.13688200555087e+01/, r(4) /4.65807828718470e+00/, &
           r(5) /2.82094791773523e-01/
      data s(1) /9.41537750555460e+01/, s(2) /1.87114811799590e+02/, &
           s(3) /9.90191814623914e+01/, s(4) /1.80124575948747e+01/

  ax = abs ( x )

  if ( ax <= 0.5 ) then

    t = x*x
    top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0
    bot = ((b(1)*t + b(2))*t + b(3))*t + 1.0
    erf = x*(top/bot)
  
  else if ( ax <= 4.0 ) then

    top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax &
                        + p(6))*ax + p(7))*ax + p(8)
    bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax &
                        + q(6))*ax + q(7))*ax + q(8)
    erf = 0.5 + (0.5 - exp(-x*x)*top/bot)
    if ( x < 0.0 ) then
      erf = -erf
    end if

  else if ( ax < 5.8 ) then

    x2 = x*x
    t = 1.0/x2
    top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
    bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + 1.0
    erf = (c - top/(x2*bot)) / ax
    erf = 0.5 + (0.5 - exp(-x2)*erf)

    if ( x < 0.0 ) then
      erf = -erf
    end if

  else

   erf = sign ( 1.0, x )

  end if

  return
end
subroutine erf_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ERF_VALUES returns some values of the ERF or "error" function for testing.
!
!  Modified:
!
!    17 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and
!    N_DATA is set to the index of the test data.  On each subsequent
!    call, N_DATA is incremented and that test data is returned.  When
!    there is no more test data, N_DATA is set to 0.
!
!    Output, real X, the argument of the function.
!
!    Output, real FX, the value of the function.
!
  implicit none

  integer, parameter :: nmax = 21

  real, save, dimension ( nmax ) :: bvec = (/ &
    0.0000000000E+00, 0.1124629160E+00, 0.2227025892E+00, 0.3286267595E+00, &
    0.4283923550E+00, 0.5204998778E+00, 0.6038560908E+00, 0.6778011938E+00, &
    0.7421009647E+00, 0.7969082124E+00, 0.8427007929E+00, 0.8802050696E+00, &
    0.9103139782E+00, 0.9340079449E+00, 0.9522851198E+00, 0.9661051465E+00, &
    0.9763483833E+00, 0.9837904586E+00, 0.9890905016E+00, 0.9927904292E+00, &
    0.9953222650E+00 /)
  real fx
  integer n_data
  real x
  real, save, dimension ( nmax ) :: xvec = (/ &
    0.0E+00, 0.1E+00, 0.2E+00, 0.3E+00, &
    0.4E+00, 0.5E+00, 0.6E+00, 0.7E+00, &
    0.8E+00, 0.9E+00, 1.0E+00, 1.1E+00, &
    1.2E+00, 1.3E+00, 1.4E+00, 1.5E+00, &
    1.6E+00, 1.7E+00, 1.8E+00, 1.9E+00, &
    2.0E+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    x = 0.0E+00
    fx = 0.0E+00
  else
    x = xvec(n_data)
    fx = bvec(n_data)
  end if

  return
end
function esum ( mu, x )

!*****************************************************************************80
!
!! ESUM evaluates exp(mu + x).
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real esum
  integer mu
  real w
  real x

  if (x .gt. 0.0) go to 10

      if (mu < 0) go to 20
         w = mu + x
         if (w .gt. 0.0) go to 20
         esum = exp(w)
         return

   10 if (mu .gt. 0) go to 20
         w = mu + x
         if (w < 0.0) go to 20
         esum = exp(w)
         return

   20 w = mu

  esum = exp(w)*exp(x)

  return
end
function exparg ( l )

!*****************************************************************************80
!
!! EXPARG reports the largest safe arguments for EXP(X).
!
!    if l = 0 then  exparg(l) = the largest positive w for which
!    exp(w) can be computed.
!
!    if l is nonzero then  exparg(l) = the largest negative w for
!    which the computed value of exp(w) is nonzero.
!
!    note... only an approximate value for exparg(l) is needed.
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  integer b
  real exparg
  integer ipmpar
  integer l
  integer m
  real lnb

  b = ipmpar(4)
  if (b /= 2) go to 10
         lnb = .69314718055995
         go to 50
   10 if (b /= 8) go to 20
         lnb = 2.0794415416798
         go to 50
   20 if (b /= 16) go to 30
         lnb = 2.7725887222398
         go to 50
   30 lnb = log ( float ( b ) )

   50 if (l == 0) go to 60
         m = ipmpar(6) - 1
         exparg = 0.99999 * (m * lnb)
         return
   60 m = ipmpar(7)
      exparg = 0.99999 * (m * lnb)

  return
end
function fpser ( a, b, x, eps )

!*****************************************************************************80
!
!! FPSER evaluates Ix(A,B) for small B and X <= 0.5.
!
!          for b < min(eps,eps*a) and x <= 0.5.
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real an
  real b
  real c
  real eps
  real exparg
  real fpser
  real s
  real t
  real tol
  real x
!
!  Set FPSER = X**A.
!
  fpser = 1.0

  if ( 1.0E-03 * eps < a ) then
    fpser = 0.0
    t = a * log ( x )
    if ( t < exparg ( 1 ) ) then
      return
    end if
    fpser = exp ( t )
  end if
!
!  Note that 1/b(a,b) = b
!
  fpser = ( b / a ) * fpser
  tol = eps/a
  an = a + 1.0
  t = x
  s = t/an

  do

    an = an + 1.0
    t = x * t
    c = t / an
    s = s + c

    if ( abs ( c ) <= tol ) then
      exit
    end if

  end do

  fpser = fpser * ( 1.0E+00 + a * s )

  return
end
function gam1 ( a )

!*****************************************************************************80
!
!! GAM1 computes 1/gamma(a+1) - 1  for -0.5 <= a <= 1.5
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real bot
  real d
  real gam1
  real p(7)
  real q(5)
  real r(9)
  real s1
  real s2
  real t
  real top
  real w

      data p(1)/ .577215664901533e+00/, p(2)/-.409078193005776e+00/, &
           p(3)/-.230975380857675e+00/, p(4)/ .597275330452234e-01/, &
           p(5)/ .766968181649490e-02/, p(6)/-.514889771323592e-02/, &
           p(7)/ .589597428611429e-03/

      data q(1)/ .100000000000000e+01/, q(2)/ .427569613095214e+00/, &
           q(3)/ .158451672430138e+00/, q(4)/ .261132021441447e-01/, &
           q(5)/ .423244297896961e-02/

      data r(1)/-.422784335098468e+00/, r(2)/-.771330383816272e+00/, &
           r(3)/-.244757765222226e+00/, r(4)/ .118378989872749e+00/, &
           r(5)/ .930357293360349e-03/, r(6)/-.118290993445146e-01/, &
           r(7)/ .223047661158249e-02/, r(8)/ .266505979058923e-03/, &
           r(9)/-.132674909766242e-03/

      data s1  / .273076135303957e+00/, s2  / .559398236957378e-01/

      t = a
      d = a - 0.5
      if (d .gt. 0.0) t = d - 0.5
      if (t) 30,10,20

   10 gam1 = 0.0
      return

   20 top = (((((p(7)*t + p(6))*t + p(5))*t + p(4))*t + p(3))*t &
                        + p(2))*t + p(1)
      bot = (((q(5)*t + q(4))*t + q(3))*t + q(2))*t + 1.0
      w = top/bot
      if (d .gt. 0.0) go to 21
         gam1 = a*w
         return
   21 gam1 = (t/a)*((w - 0.5) - 0.5)
      return

   30 top = (((((((r(9)*t + r(8))*t + r(7))*t + r(6))*t + r(5))*t &
                          + r(4))*t + r(3))*t + r(2))*t + r(1)
      bot = (s2*t + s1)*t + 1.0
      w = top/bot
      if ( 0.0 < d ) go to 31
         gam1 = a*((w + 0.5) + 0.5)
         return
   31 gam1 = t*w/a

  return
end
function gamln1 ( a )

!*****************************************************************************80
!
!! GAMLN1 evaluates ln(gamma(1 + a)) for -0.2 <= A <= 1.25
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real gamln1
  real p0
  real p1
  real p2
  real p3
  real p4
  real p5
  real p6
  real q1
  real q2
  real q3
  real q4
  real q5
  real q6
  real r0
  real r1
  real r2
  real r3
  real r4
  real r5
  real s1
  real s2
  real s3
  real s4
  real s5
  real w
  real x

      data p0/ .577215664901533e+00/, p1/ .844203922187225e+00/, &
           p2/-.168860593646662e+00/, p3/-.780427615533591e+00/, &
           p4/-.402055799310489e+00/, p5/-.673562214325671e-01/, &
           p6/-.271935708322958e-02/
      data q1/ .288743195473681e+01/, q2/ .312755088914843e+01/, &
           q3/ .156875193295039e+01/, q4/ .361951990101499e+00/, &
           q5/ .325038868253937e-01/, q6/ .667465618796164e-03/

      data r0/.422784335098467e+00/,  r1/.848044614534529e+00/, &
           r2/.565221050691933e+00/,  r3/.156513060486551e+00/, &
           r4/.170502484022650e-01/,  r5/.497958207639485e-03/
      data s1/.124313399877507e+01/,  s2/.548042109832463e+00/, &
           s3/.101552187439830e+00/,  s4/.713309612391000e-02/, &
           s5/.116165475989616e-03/

  if ( a < 0.6 ) then
      w = ((((((p6*a + p5)*a + p4)*a + p3)*a + p2)*a + p1)*a + p0)/ &
          ((((((q6*a + q5)*a + q4)*a + q3)*a + q2)*a + q1)*a + 1.0)
      gamln1 = -a*w
  else
    x = (a - 0.5) - 0.5
      w = (((((r5*x + r4)*x + r3)*x + r2)*x + r1)*x + r0)/ &
          (((((s5*x + s4)*x + s3)*x + s2)*x + s1)*x + 1.0)
      gamln1 = x*w
  end if

  return
end
function gamln ( a )

!*****************************************************************************80
!
!! GAMLN evaluates ln(gamma(a)) for positive A.
!
!  Author:
!
!    Alfred Morris
!    Naval Surface Warfare Center
!    Dahlgren, Virginia
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none
!
!     d = 0.5*(ln(2*pi) - 1)
!
  real a
  real c0
  real c1
  real c2
  real c3
  real c4
  real c5
  real, parameter :: d = 0.418938533204673
  real gamln
  real gamln1
  integer i
  integer n
  real t
  real w


      data c0/.833333333333333e-01/, c1/-.277777777760991e-02/, &
           c2/.793650666825390e-03/, c3/-.595202931351870e-03/, &
           c4/.837308034031215e-03/, c5/-.165322962780713e-02/

  if ( a <= 0.8 ) then

    gamln = gamln1 ( a ) - log ( a )

  else if ( a <= 2.25 ) then

    t = (a - 0.5) - 0.5
    gamln = gamln1 ( t )

  else if ( a < 10.0 ) then

    n = a - 1.25
    t = a
    w = 1.0
    do i = 1, n
      t = t - 1.0
      w = t*w
    end do
    gamln = gamln1(t - 1.0) + log ( w )

  else

    t = (1.0/a)**2
    w = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0)/a
    gamln = (d + w) + (a - 0.5)*( log ( a ) - 1.0)

  end if

  return
end
subroutine gamma_inc_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
!
!  Discussion:
!
!    The (normalized) incomplete Gamma function P(A,X) is defined as:
!
!      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
!
!    With this definition, for all A and X,
!
!      0 <= PN(A,X) <= 1
!
!    and
!
!      PN(A,INFINITY) = 1.0
!
!    In Mathematica, the function can be evaluated by:
!
!      1 - GammaRegularized[A,X]
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
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 4 ) A, the parameter of the function.
!
!    Output, real ( kind = 4 ) X, the argument of the function.
!
!    Output, real ( kind = 4 ) FX, the value of the function.
!
  implicit none

  integer, parameter :: n_max = 20

  real ( kind = 4 ) a
  real ( kind = 4 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.10E+00, &
     0.10E+00, &
     0.10E+00, &
     0.50E+00, &
     0.50E+00, &
     0.50E+00, &
     0.10E+01, &
     0.10E+01, &
     0.10E+01, &
     0.11E+01, &
     0.11E+01, &
     0.11E+01, &
     0.20E+01, &
     0.20E+01, &
     0.20E+01, &
     0.60E+01, &
     0.60E+01, &
     0.11E+02, &
     0.26E+02, &
     0.41E+02 /)
  real ( kind = 4 ) fx
  real ( kind = 4 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.7382350532339351E+00, &
     0.9083579897300343E+00, &
     0.9886559833621947E+00, &
     0.3014646416966613E+00, &
     0.7793286380801532E+00, &
     0.9918490284064973E+00, &
     0.9516258196404043E-01, &
     0.6321205588285577E+00, &
     0.9932620530009145E+00, &
     0.7205974576054322E-01, &
     0.5891809618706485E+00, &
     0.9915368159845525E+00, &
     0.1392920235749422E+00, &
     0.7768698398515702E+00, &
     0.9990881180344455E+00, &
     0.4202103819530612E-01, &
     0.9796589705830716E+00, &
     0.9226039842296429E+00, &
     0.4470785799755852E+00, &
     0.7444549220718699E+00 /)
  integer n_data
  real ( kind = 4 ) x
  real ( kind = 4 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.30E-01, &
     0.30E+00, &
     0.15E+01, &
     0.75E-01, &
     0.75E+00, &
     0.35E+01, &
     0.10E+00, &
     0.10E+01, &
     0.50E+01, &
     0.10E+00, &
     0.10E+01, &
     0.50E+01, &
     0.15E+00, &
     0.15E+01, &
     0.70E+01, &
     0.25E+01, &
     0.12E+02, &
     0.16E+02, &
     0.25E+02, &
     0.45E+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0E+00
    x = 0.0E+00
    fx = 0.0E+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_log_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_LOG_VALUES returns some values of the Log Gamma function for testing.
!
!  Modified:
!
!    17 April 2001
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
!  Parameters:
!
!    Input/output, integer N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and
!    N_DATA is set to the index of the test data.  On each subsequent
!    call, N_DATA is incremented and that test data is returned.  When
!    there is no more test data, N_DATA is set to 0.
!
!    Output, real X, the argument of the function.
!
!    Output, real FX, the value of the function.
!
  implicit none

  integer, parameter :: nmax = 18

  real, save, dimension ( nmax ) :: bvec = (/ &
     1.524064183E+00,    0.7966780066E+00,   0.3982337117E+00,  &
     0.1520599127E+00,   0.000000000E+00,   -0.04987246543E+00, &
    -0.08537410945E+00, -0.1081747934E+00,  -0.1196128950E+00,  & 
    -0.1207822040E+00,  -0.1125917658E+00,  -0.09580771625E+00, &
    -0.07108385116E+00, -0.03898428380E+00,  0.000000000E+00,   &
    12.80182743E+00,    39.33988571E+00,    71.25704193E+00 /)
  real fx
  integer n_data
  real x
  real, save, dimension ( nmax ) :: xvec = (/ &
    0.2E+00,  0.4E+00,  0.6E+00,  0.8E+00, &
    1.0E+00,  1.1E+00,  1.2E+00,  1.3E+00, &
    1.4E+00,  1.5E+00,  1.6E+00,  1.7E+00, &
    1.8E+00,  1.9E+00,  2.0E+00, 10.0E+00, &
   20.0E+00, 30.0E+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    x = 0.0E+00
    fx = 0.0E+00
  else
    x = xvec(n_data)
    fx = bvec(n_data)
  end if

  return
end
subroutine grat1 ( a, x, r, p, q, eps )

!*****************************************************************************80
!
!! GRAT1 evaluates the incomplete Gamma ratio functions P(A,X) and Q(A,X).
!
!  Discussion:
!
!    It is assumed that a <= 1.  eps is the tolerance to be used.
!    the input argument r has the value e**(-x)*x**a/gamma(a).
!
!  Modified:
!
!    28 August 2004
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real a2n
  real a2nm1
  real am0
  real an
  real an0
  real b2n
  real b2nm1
  real c
  real cma
  real eps
  real erf
  real erfc1
  real g
  real gam1
  real h
  real j
  real l
  real p
  real q
  real r
  real rexp
  real sum2
  real t
  real tol
  real w
  real x
  real z

  if (a*x == 0.0) go to 130
  if (a == 0.5) go to 120
  if (x < 1.1) go to 10
  go to 50
!
!  Taylor series for p(a,x)/x**a
!
10 continue

      an = 3.0
      c = x
      sum2 = x/(a + 3.0)
      tol = 0.1 * eps / ( a + 1.0 )

      do
        an = an + 1.0
        c = -c*(x/an)
        t = c/(a + an)
        sum2 = sum2 + t
        if ( abs ( t ) <= tol ) then
          exit
        end if
      end do

      j = a*x*((sum2/6.0 - 0.5/(a + 2.0))*x + 1.0/(a + 1.0))

      z = a * log ( x )
      h = gam1(a)
      g = 1.0 + h
      if (x < 0.25) go to 20
         if (a < x/2.59) go to 40
         go to 30
   20 if (z .gt. -0.13394) go to 40

   30 w = exp(z)
      p = w*g*(0.5 + (0.5 - j))
      q = 0.5 + (0.5 - p)
      return

   40 l = rexp(z)
      w = 0.5 + (0.5 + l)
      q = (w*j - l)*g - h
      if (q < 0.0) go to 110
      p = 0.5 + (0.5 - q)
      return
!
!  continued fraction expansion
!
   50 a2nm1 = 1.0
      a2n = 1.0
      b2nm1 = x
      b2n = x + (1.0 - a)
      c = 1.0
   51    a2nm1 = x*a2n + c*a2nm1
         b2nm1 = x*b2n + c*b2nm1
         am0 = a2nm1/b2nm1
         c = c + 1.0
         cma = c - a
         a2n = a2nm1 + cma*a2n
         b2n = b2nm1 + cma*b2n
         an0 = a2n/b2n
         if (abs(an0 - am0) .ge. eps*an0) go to 51
      q = r*an0
      p = 0.5 + (0.5 - q)
      return
!
!  special cases
!
  100 p = 0.0
      q = 1.0
      return

  110 p = 1.0
      q = 0.0
      return

  120 if (x .ge. 0.25) go to 121
      p = erf(sqrt(x))
      q = 0.5 + (0.5 - p)
      return
  121 q = erfc1(0,sqrt(x))
      p = 0.5 + (0.5 - q)
      return

  130 if (x <= a) go to 100
      go to 110

end
function gsumln ( a, b )

!*****************************************************************************80
!
!! GSUMLN evaluates the function Log ( Gamma ( A + B ) ) in a special range.
!
!          for 1 <= a <= 2  and  1 <= b <= 2
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real a
  real alnrel
  real b
  real gamln1
  real gsumln
  real x

  x = dble ( a ) + dble ( b ) - 2.0d0

  if ( x <= 0.25 ) then
    gsumln = gamln1 ( 1.0 + x )
  else if (x <= 1.25 ) then
    gsumln = gamln1 ( x ) + alnrel ( x )
  else
    gsumln = gamln1 ( x - 1.0 ) + log ( x * ( 1.0 + x ) )
  end if

  return
end
function ipmpar ( i )

!*****************************************************************************80
!
!! IPMPAR provides the integer machine constants for the computer
!     that is used. it is assumed that the argument i is an integer
!     having one of the values 1-10. ipmpar(i) has the value ...
!
!  integers.
!
!     assume integers are represented in the n-digit, base-a form
!
!               sign ( x(n-1)*a**(n-1) + ... + x(1)*a + x(0) )
!
!               where 0 <= x(i) < a for i=0,...,n-1.
!
!     ipmpar(1) = a, the base.
!
!     ipmpar(2) = n, the number of base-a digits.
!
!     ipmpar(3) = a**n - 1, the largest magnitude.
!
!  floating-point numbers.
!
!     it is assumed that the single and double precision floating
!     point arithmetics have the same base, say b, and that the
!     nonzero numbers are represented in the form
!
!               sign (b**e) * (x(1)/b + ... + x(m)/b**m)
!
!               where x(i) = 0,1,...,b-1 for i=1,...,m,
!               x(1) .ge. 1, and emin <= e <= emax.
!
!     ipmpar(4) = b, the base.
!
!  single-precision
!
!     ipmpar(5) = m, the number of base-b digits.
!
!     ipmpar(6) = emin, the smallest exponent e.
!
!     ipmpar(7) = emax, the largest exponent e.
!
!  double-precision
!
!     ipmpar(8) = m, the number of base-b digits.
!
!     ipmpar(9) = emin, the smallest exponent e.
!
!     ipmpar(10) = emax, the largest exponent e.
!
!     to define this function for the computer being used, activate
!     the data statments for the computer by removing the c from
!     column 1. (all the other data statements should have c in
!     column 1.)
!
!     ipmpar is an adaptation of the function i1mach, written by
!     p.a. fox, a.d. hall, and n.l. schryer (bell laboratories).
!     ipmpar was formed by a.h. morris (nswc). the constants are
!     from bell laboratories, nswc, and other sources.
!
      integer imach(10)
  integer ipmpar
!
!     machine constants for amdahl machines.
!
!     data imach( 1) /   2 /
!     data imach( 2) /  31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /  16 /
!     data imach( 5) /   6 /
!     data imach( 6) / -64 /
!     data imach( 7) /  63 /
!     data imach( 8) /  14 /
!     data imach( 9) / -64 /
!     data imach(10) /  63 /
!
!     machine constants for the at&t 3b series, at&t
!     pc 7300, and at&t 6300.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     machine constants for the burroughs 1700 system.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   33 /
!     data imach( 3) / 8589934591 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -256 /
!     data imach( 7) /  255 /
!     data imach( 8) /   60 /
!     data imach( 9) / -256 /
!     data imach(10) /  255 /
!
!     machine constants for the burroughs 5700 system.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   39 /
!     data imach( 3) / 549755813887 /
!     data imach( 4) /    8 /
!     data imach( 5) /   13 /
!     data imach( 6) /  -50 /
!     data imach( 7) /   76 /
!     data imach( 8) /   26 /
!     data imach( 9) /  -50 /
!     data imach(10) /   76 /
!
!     machine constants for the burroughs 6700/7700 systems.
!
!     data imach( 1) /      2 /
!     data imach( 2) /     39 /
!     data imach( 3) / 549755813887 /
!     data imach( 4) /      8 /
!     data imach( 5) /     13 /
!     data imach( 6) /    -50 /
!     data imach( 7) /     76 /
!     data imach( 8) /     26 /
!     data imach( 9) / -32754 /
!     data imach(10) /  32780 /
!
!     machine constants for the cdc 6000/7000 series
!     60 bit arithmetic, and the cdc cyber 995 64 bit
!     arithmetic (nos operating system).
!
!     data imach( 1) /    2 /
!     data imach( 2) /   48 /
!     data imach( 3) / 281474976710655 /
!     data imach( 4) /    2 /
!     data imach( 5) /   48 /
!     data imach( 6) / -974 /
!     data imach( 7) / 1070 /
!     data imach( 8) /   95 /
!     data imach( 9) / -926 /
!     data imach(10) / 1070 /
!
!     machine constants for the cdc cyber 995 64 bit
!     arithmetic (nos/ve operating system).
!
!     data imach( 1) /     2 /
!     data imach( 2) /    63 /
!     data imach( 3) / 9223372036854775807 /
!     data imach( 4) /     2 /
!     data imach( 5) /    48 /
!     data imach( 6) / -4096 /
!     data imach( 7) /  4095 /
!     data imach( 8) /    96 /
!     data imach( 9) / -4096 /
!     data imach(10) /  4095 /
!
!     machine constants for the cray 1, xmp, 2, and 3.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    63 /
!     data imach( 3) / 9223372036854775807 /
!     data imach( 4) /     2 /
!     data imach( 5) /    47 /
!     data imach( 6) / -8189 /
!     data imach( 7) /  8190 /
!     data imach( 8) /    94 /
!     data imach( 9) / -8099 /
!     data imach(10) /  8190 /
!
!     machine constants for the data general eclipse s/200.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /   16 /
!     data imach( 5) /    6 /
!     data imach( 6) /  -64 /
!     data imach( 7) /   63 /
!     data imach( 8) /   14 /
!     data imach( 9) /  -64 /
!     data imach(10) /   63 /
!
!     machine constants for the harris 220.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   23 /
!     data imach( 3) / 8388607 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   38 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     machine constants for the honeywell 600/6000
!     and dps 8/70 series.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   63 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     machine constants for the hp 2100
!     3 word double precision option with ftn4
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   39 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     machine constants for the hp 2100
!     4 word double precision option with ftn4
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   55 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     machine constants for the hp 9000.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -126 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     machine constants for the ibm 360/370 series,
!     the icl 2900, the itel as/6, the xerox sigma
!     5/7/9 and the sel systems 85/86.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /   16 /
!     data imach( 5) /    6 /
!     data imach( 6) /  -64 /
!     data imach( 7) /   63 /
!     data imach( 8) /   14 /
!     data imach( 9) /  -64 /
!     data imach(10) /   63 /
!
!     machine constants for the ibm pc.
!
      data imach( 1) /     2 /
      data imach( 2) /    31 /
      data imach( 3) / 2147483647 /
      data imach( 4) /     2 /
      data imach( 5) /    24 /
      data imach( 6) /  -125 /
      data imach( 7) /   128 /
      data imach( 8) /    53 /
      data imach( 9) / -1021 /
      data imach(10) /  1024 /
!
!     machine constants for the macintosh ii - absoft
!     macfortran ii.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     machine constants for the microvax - vms fortran.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     machine constants for the pdp-10 (ka processor).
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   54 /
!     data imach( 9) / -101 /
!     data imach(10) /  127 /
!
!     machine constants for the pdp-10 (ki processor).
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   62 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     machine constants for the pdp-11 fortran supporting
!     32-bit integer arithmetic.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     machine constants for the sequent balance 8000.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     machine constants for the silicon graphics iris-4d
!     series (mips r3000 processor).
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     machine constants for the sun 3.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     machine constants for the univac 1100 series.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   60 /
!     data imach( 9) /-1024 /
!     data imach(10) / 1023 /
!
!     machine constants for the vax 11/780.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
      ipmpar = imach(i)

  return
end
function psi ( xx )

!*****************************************************************************80
!
!! PSI evaluates the Psi or Digamma function.
!
!     psi(xx) is assigned the value 0 when the digamma function cannot
!     be computed.
!
!     the main computation involves evaluation of rational chebyshev
!     approximations published in math. comp. 27, 123-127(1973) by
!     cody, strecok and thacher.
!
!     psi was written at argonne national laboratory for the funpack
!     package of special function subroutines. psi was modified by
!     a.h. morris (nswc).
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
!  Local parameters:
!
!    Local, double precision DX0, zero of PSI.
!
!    Local, real PIOV4 = pi/4.
!
!    Local, real XMAX1, the smallest positive floating point constant
!    with entirely integer representation.  Also used as negative of 
!    lower bound on acceptable negative arguments and as the positive 
!    argument beyond which PSI may be represented as log(x).
!
  implicit none

  real aug
  real den
  double precision, parameter :: dx0 = 1.461632144968362341262659542325721325d0
  integer i
  integer ipmpar
  integer m
  integer n
  integer nq
  real p1(7)
  real p2(4)
  real, parameter :: piov4 = 0.785398163397448e0
  real psi
  real q1(6)
  real q2(4)
  real sgn
  real upper
  real w
  real x
  real xmax1
  real xmx0
  real xsmall
  real xx
  real z
!
!  coefficients for rational approximation of
!  psi(x) / (x - x0),  0.5 <= x <= 3.0
!
      data p1(1)/.895385022981970e-02/,  p1(2)/.477762828042627e+01/, &
           p1(3)/.142441585084029e+03/,  p1(4)/.118645200713425e+04/, &
           p1(5)/.363351846806499e+04/,  p1(6)/.413810161269013e+04/, &
           p1(7)/.130560269827897e+04/
      data q1(1)/.448452573429826e+02/,  q1(2)/.520752771467162e+03/, &
           q1(3)/.221000799247830e+04/,  q1(4)/.364127349079381e+04/, &
           q1(5)/.190831076596300e+04/,  q1(6)/.691091682714533e-05/
!
!  coefficients for rational approximation of
!  psi(x) - ln(x) + 1 / (2*x),  x .gt. 3.0
!
      data p2(1)/-.212940445131011e+01/, p2(2)/-.701677227766759e+01/, &
           p2(3)/-.448616543918019e+01/, p2(4)/-.648157123766197e+00/
      data q2(1)/ .322703493791143e+02/, q2(2)/ .892920700481861e+02/, &
           q2(3)/ .546117738103215e+02/, q2(4)/ .777788548522962e+01/
!
!  machine dependent constants ...
!
!  xsmall = absolute argument below which pi*cotan(pi*x)
!  may be represented by 1/x.
!
      xmax1 = ipmpar(3)
      xmax1 =  min ( xmax1, 1.0 / epsilon ( xmax1 ) )
      xsmall = 1.e-9

      x = xx
      aug = 0.0e0
      if (x .ge. 0.5e0) go to 200
!
!  x < 0.5,  use reflection formula
!  psi(1-x) = psi(x) + pi * cotan(pi*x)
!
      if ( xsmall < abs(x) ) go to 100
      if (x == 0.0e0) go to 400
!
!  0 < abs(x) <= xsmall.  use 1/x as a substitute
!  for  pi*cotan(pi*x)
!
      aug = -1.0e0 / x
      go to 150
!
!  reduction of argument for cotan
!
  100 w = - x
      sgn = piov4
      if ( 0.0e0 < w ) go to 120
      w = - w
      sgn = -sgn
!
!  Error exit if x <= -xmax1
!
  120 if (w .ge. xmax1) go to 400
      nq = int(w)
      w = w - float(nq)
      nq = int(w*4.0e0)
      w = 4.0e0 * (w - float(nq) * .25e0)
!
!  w is now related to the fractional part of  4.0 * x.
!  adjust argument to correspond to values in first
!  quadrant and determine sign
!
      n = nq / 2
      if ((n+n) /= nq) then
        w = 1.0e0 - w
      end if
      z = piov4 * w
      m = n / 2
      if ((m+m) /= n) then 
        sgn = - sgn
      end if
!
!  Determine final value for  -pi*cotan(pi*x)
!
      n = (nq + 1) / 2
      m = n / 2
      m = m + m
      if (m /= n) go to 140
!
!  Check for singularity
!
      if (z == 0.0e0) go to 400
!
!  use cos/sin as a substitute for cotan, and
!  sin/cos as a substitute for tan
!
      aug = sgn * ((cos(z) / sin(z)) * 4.0e0)
      go to 150
  140 aug = sgn * ((sin(z) / cos(z)) * 4.0e0)
  150 x = 1.0e0 - x
  200 if (x .gt. 3.0e0) go to 300
!
!  0.5 <= x <= 3.0
!
      den = x
      upper = p1(1) * x

      do i = 1, 5
         den = (den + q1(i)) * x
         upper = (upper + p1(i+1)) * x
      end do

      den = (upper + p1(7)) / (den + q1(6))
      xmx0 = dble(x) - dx0
      psi = den * xmx0 + aug
      return
!
!  if x .ge. xmax1, psi = ln(x)
!
  300 if (x .ge. xmax1) go to 350
!
!  3.0 < x < xmax1
!
      w = 1.0e0 / (x * x)
      den = w
      upper = p2(1) * w

      do i = 1, 3
         den = (den + q2(i)) * w
         upper = (upper + p2(i+1)) * w
      end do

      aug = upper / (den + q2(4)) - 0.5e0 / x + aug
  350 psi = aug + log ( x )
      return
!
!  error return
!
  400 psi = 0.0e0

  return
end
subroutine psi_values ( n, x, fx )

!*****************************************************************************80
!
!! PSI_VALUES returns some values of the Psi or Digamma function for testing.
!
!  Discussion:
!
!    PSI(X) = d LN ( GAMMA ( X ) ) / d X = GAMMA'(X) / GAMMA(X)
!
!    PSI(1) = - Euler's constant.
!
!    PSI(X+1) = PSI(X) + 1 / X.
!
!  Modified:
!
!    17 May 2001
!
!  Reference:
!
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, real X, the argument of the function.
!
!    Output, real FX, the value of the function.
!
  implicit none

  integer, parameter :: nmax = 11

  real fx
  real, save, dimension ( nmax ) :: fxvec = (/ &
    -0.5772156649E+00, -0.4237549404E+00, -0.2890398966E+00, &
    -0.1691908889E+00, -0.0613845446E+00, -0.0364899740E+00, &
     0.1260474528E+00,  0.2085478749E+00,  0.2849914333E+00, &  
     0.3561841612E+00,  0.4227843351E+00 /)
  integer n
  real x
  real, save, dimension ( nmax ) :: xvec = (/ &
    1.0E+00,  1.1E+00,  1.2E+00,  &
    1.3E+00,  1.4E+00,  1.5E+00,  &
    1.6E+00,  1.7E+00,  1.8E+00,  &
    1.9E+00,  2.0E+00 /)

  if ( n < 0 ) then
    n = 0
  end if

  n = n + 1

  if ( nmax < n ) then
    n = 0
    x = 0.0E+00
    fx = 0.0E+00
    return
  end if

  x = xvec(n)
  fx = fxvec(n)

  return
end
function rexp ( x )

!*****************************************************************************80
!
!! REXP evaluates the function Exp(X) - 1.
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real, parameter :: p1 =  0.914041914819518e-09
  real, parameter :: p2 =  0.238082361044469e-01
  real, parameter :: q1 = -0.499999999085958e+00
  real, parameter :: q2 =  0.107141568980644e+00
  real, parameter :: q3 = -0.119041179760821e-01
  real, parameter :: q4 =  0.595130811860248e-03
  real rexp
  real w
  real x

  if ( abs ( x ) <= 0.15 ) then

    rexp = x*(((p2*x + p1)*x + 1.0) / ((((q4*x + q3)*x + q2)*x &
      + q1)*x + 1.0))

  else if ( x < 0.0 ) then

    w = exp ( x )
    rexp = ( w - 0.5 ) - 0.5

  else

    w = exp ( x )
    rexp = w * ( 0.5 + ( 0.5 - 1.0 / w ) )

  end if

  return
end
function rlog1 ( x )

!*****************************************************************************80
!
!! RLOG1 evaluates the function X - Log ( 1 + X ).
!
!  Reference:
!
!    Armido Didonato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, 1992, pages 360-373.
!
  implicit none

  real, parameter :: a = 0.566749439387324e-01
  real, parameter :: b = 0.456512608815524e-01
  real h
  real, parameter :: p0 =  0.333333333333333e+00
  real, parameter :: p1 = -0.224696413112536e+00
  real, parameter :: p2 =  0.620886815375787e-02
  real, parameter :: q1 = -0.127408923933623e+01
  real, parameter :: q2 =  0.354508718369557e+00
  real r
  real rlog1
  real t
  real w
  real w1
  real x

  if ( x < -0.39 ) then
    w = ( x + 0.5 ) + 0.5
    rlog1 = x - log ( w )
  else if ( x < -0.18 ) then
    h = ( dble ( x ) + 0.3d0 ) / 0.7
    w1 = a - h*0.3
    r = h/(h + 2.0)
    t = r * r
    w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0 )
    rlog1 = 2.0*t*(1.0/(1.0 - r) - r*w) + w1
  else if ( x <= 0.18 ) then
    h = x
    w1 = 0.0
    r = h/(h + 2.0)
    t = r*r
    w = ((p2*t + p1)*t + p0)/((q2*t + q1)*t + 1.0 )
    rlog1 = 2.0*t*(1.0/(1.0 - r) - r*w) + w1
  else if ( x <= 0.57 ) then
    h = 0.75d0*dble(x) - 0.25d0
    w1 = b + h/3.0
    r = h/(h + 2.0)
    t = r*r
    w = ((p2*t + p1)*t + p0)/((q2*t + q1)*t + 1.0)
    rlog1 = 2.0 * t * ( 1.0 / ( 1.0 - r ) - r * w ) + w1
  else
    w = ( x + 0.5 ) + 0.5
    rlog1 = x - log ( w )
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
