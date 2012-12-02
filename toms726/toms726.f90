function alga_r4 ( x )

!*****************************************************************************80
!
!! ALGA_R4 evaluates the logarithm of the gamma function.
!
!  Discussion:
!
!    This is an auxiliary routine, not optimized in any
!    sense, for evaluating the logarithm of the gamma function for positive
!    arguments X.  It is called by the routine GAMMA.  The integer m0
!    in the first executable statement is the smallest integer  m  such
!    that  1*3*5* ... *(2*m+1)/(2**m)  is greater than or equal to the
!    largest machine-representable number.
!
!    The routine is based on a rational approximation valid on [.5,1.5]
!    due to W.J. Cody and K.E. Hillstrom; see Math. Comp. 21, 1967, 198-203,
!    in particular the case  n = 7  in Table II.
!
!    For the computation of m0 it calls upon T_FUNCTION which
!    evaluates the inverse function  t  =  t(y)  of  y = t ln t.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the
!    Gamma Function,
!    Mathematics of Computation,
!    Volume 21, Number 98, April 1967, pages 198-203.
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, real X, the argument.
!
!    Output, real ALGA_R4, the logarithm of gamma(X).
!
  implicit none

  real alga_r4
  real, dimension ( 8 ) :: cnum = (/ &
       4.120843185847770E+00, &
      85.68982062831317E+00,  &
     243.175243524421E+00,    &
    -261.7218583856145E+00,   &
    -922.2613728801522E+00,   &
    -517.6383498023218E+00,   &
     -77.41064071332953E+00,  &
      -2.208843997216182E+00 /)
  real, dimension ( 8 ) :: cden = (/ &
      1.0E+00, &
     45.64677187585908E+00, &
    377.8372484823942E+00, &
    951.323597679706E+00, &
    846.0755362020782E+00, &
    262.3083470269460E+00, &
     24.43519662506312E+00, &
      0.4097792921092615E+00 /)
  integer k
  integer m
  integer m0
  integer mm1
  real p
  real sden
  real snum
  real t_function
  real x
  real xe
  real xi
!
!  The constants in the statement below are  exp(1.)  and  .5*log(8.).
!
  m0 = 2.71828E+00 * t_function ( ( log ( huge ( x ) ) - 1.03972E+00 ) &
    / 2.71828E+00 )

  xi = aint ( x )

  if ( 0.5E+00 < x - xi ) then
    xi = xi + 1.0E+00
  end if

  m = int ( xi ) - 1
!
!  Computation of log gamma on the standard interval (1/2,3/2]
!
  xe = x - real ( m, kind = 4 )
  snum = cnum(1)
  sden = cden(1)
  do k = 2, 8
    snum = xe * snum + cnum(k)
    sden = xe * sden + cden(k)
  end do
  alga_r4 = ( xe - 1.0E+00 ) * snum / sden
!
!  Computation of log gamma on (0,1/2]
!
  if ( m == -1 ) then
    alga_r4 = alga_r4 - log ( x )
    return
  else if ( m == 0 ) then
    return
  else
!
!  Computation of log gamma on (3/2,5/2]
!
    p = xe
    if ( m == 1 ) then
      alga_r4 = alga_r4 + log ( p )
      return
    else
!
!  Computation of log gamma for arguments larger than 5/2
!
      mm1 = m - 1
!
!  The else-clause in the next statement is designed to avoid possible
!  overflow in the computation of P in the if-clause, at the expense
!  of computing many logarithms.
!
      if ( m < m0 ) then
        do k = 1, mm1
          p = ( xe + real ( k, kind = 4 ) ) * p
        end do
        alga_r4 = alga_r4 + log ( p )
        return
      else
        alga_r4 = alga_r4 + log ( xe )
        do k = 1, mm1
          alga_r4 = alga_r4 + log ( xe + real ( k, kind = 4 ) )
        end do
        return
      end if
    end if
  end if

  return
end
function alga_r8 ( x )

!*****************************************************************************80
!
!! ALGA_R8 evaluates the logarithm of the gamma function.
!
!  Discussion:
!
!    A combination of recurrence and asymptotic approximation is used.
!
!    The entries in the data statement are the numerators and
!    denominators, respectively, of the quantities B[16]/(16*15),
!    B[14]/(14*13),..., B[2]/(2*1), where B[2n] are the Bernoulli
!    numbers.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the
!    Gamma Function,
!    Mathematics of Computation,
!    Volume 21, Number 98, April 1967, pages 198-203.
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) ALGA_R8, the logarithm of gamma(X).
!
  implicit none

  real ( kind = 8 ) alga_r8
  real ( kind = 8 ), dimension ( 8 ) :: dbden = (/ &
    1.224D+05, &
    1.56D+02, &
    3.6036D+05, &
    1.188D+03, &
    1.68D+03, &
    1.26D+03, &
    3.6D+02, &
    1.2D+01 /)
  real ( kind = 8 ), dimension ( 8 ) :: dbnum = (/ &
    -3.617D+03, &
     1.0D+00, &
    -6.91D+02, &
     1.0D+00, &
    -1.0D+00, &
     1.0D+00, &
    -1.0D+00, &
     1.0D+00 /)
  real ( kind = 8 ) dc
  real ( kind = 8 ) dp
  real ( kind = 8 ) dprec
  real ( kind = 8 ) ds
  real ( kind = 8 ) dt
  real ( kind = 8 ) dy
  integer i
  real ( kind = 8 ) x
  real ( kind = 4 ) y
  real ( kind = 8 ) y0
!
!  The quantity dprec in the next statement is the number of decimal
!  digits carried in double-precision floating-point arithmetic.
!
  dprec = - log10 ( epsilon ( dbnum ) )
  dc = 0.5D+00 * log ( 8.0D+00 * atan ( 1.0D+00 ) )
  dp = 1.0D+00
  dy = x
  y = real ( dy, kind = 4 )
!
!  The quantity  y0  below is the threshold value beyond which asymptotic
!  evaluation gives sufficient accuracy; see Eq. 6.1.42 in M. Abramowitz
!  and I.A. Stegun,Handbook of Mathematical Functions''. The constants
!  are .12118868...  =  ln(10)/19 and .05390522... = ln(|B[20]|/190)/19.
!
  y0 = exp ( 0.121189D+00 * dprec + 0.053905D+00 )

  do while ( y <= y0 )
    dp = dy * dp
    dy = dy + 1.0D+00
    y = real ( dy, kind = 4 )
  end do

  dt = 1.0D+00 / ( dy * dy )
!
!  The right-hand side of the next assignment statement is B[18]/(18*17).
!
  ds = 4.3867D+04 / 2.44188D+05
  do i = 1, 8
    ds = dt * ds + dbnum(i) / dbden(i)
  end do

  alga_r8 = ( dy - 0.5D+00 ) * log ( dy ) - dy + dc + ds / dy - log ( dp )

  return
end
subroutine cheb_r4 ( n, a, b, fnu, alpha, beta, s, ierr )

!*****************************************************************************80
!
!! CHEB_R4 implements the modified Chebyshev algorithm.
!
!  Discussion
!
!    This routine generates recursion coefficients ALPHA and BETA.
!
!    Given a set of polynomials  p(0),p(1),...,p(2*n-1)  satisfying
!
!      p(k+1)(x) = (x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
!      k = 0,1,...,2*n-2,
!
!      p(-1)(x) = 0,  p(0)(x)=1,
!
!    and associated modified moments
!
!      fnu(k) = integral of p(k)(x)*dlambda(x),
!      k = 0,1,...,2*n-1,
!
!    this subroutine uses the modified Chebyshev algorithm to generate the
!    recursion coefficients  alpha(k),beta(k), k = 0,1,...,n-1, for the
!    polynomials  pi(k)  orthogonal with respect to the integration
!    measure  dlambda(x), i.e.,
!
!      pi(k+1)(x) = (x-alpha(k))*pi(k)(x)-beta(k)*pi(k-1)(x),
!      k = 0,1,...,n-1,
!
!      pi(-1)(x) = 0,  pi(0)(x)=1.
!
!    On machines with limited exponent range, the occurrence of underflow
!    [overflow] in the computation of the  alpha's  and  beta's  can often
!    be avoided by multiplying all modified moments by a sufficiently large
!    [small] scaling factor and dividing the new  beta(0)  by the same
!    scaling factor.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    On Generating Orthogonal Polynomials,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 3, Number 3, 1982, pages 289-317.
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of recursion coefficients desired.
!
!    Input, real A(2*N-1), B(2*N-1), the values of A(k-1), b(k-1),
!    k = 1,2,...,2*n-1.
!
!    Input, real FNU(2*N), the values of the modified moments fnu(k-1),
!    k = 1,2,...,2*n
!
!    Output, real ALPHA(N), BETA(N), the
!    recursion coefficients  alpha(k-1),beta(k-1),
!    k = 1,2,...,n, where  beta(0)  is the total mass.
!
!    Output, real S(N), the normalization factors
!    s(k) = integral [pi(k)(x)]**2 dlambda(x), k=0,1,
!    2,...,n-1.
!
!    Output, integer IERR, an error flag.
!    0, on normal return,
!    1, if  abs(fnu(0))  is less than the machine zero,
!    2, if  n  is out of range,
!    -k  if S(K), k = 0,1,2,...,n-1, is about to underflow,
!    +k  if S(K) is about to overflow.
!
  implicit none

  integer n

  real a(2*n-1)
  real alpha(n)
  real b(2*n-1)
  real beta(n)
  real fnu(2*n)
  integer ierr
  integer k
  integer l
  integer lk
  real s(n)
  real s0(2*n)
  real s1(2*n)
  real s2(2*n)

  ierr = 0

  if ( abs ( fnu(1) ) < 10.0E+00 * tiny ( fnu(1) ) ) then
    ierr = 1
    return
  end if

  if ( n < 1 ) then
    ierr = 2
    return
  end if
!
!  Initialization.
!
  alpha(1) = a(1) + fnu(2) / fnu(1)
  beta(1) = fnu(1)

  if ( n == 1 ) then
    return
  end if

  s0(1:2*n) = 0.0E+00
  s(1) = fnu(1)
  s0(1:2*n) = 0.0E+00
  s1(1:2*n) = fnu(1:2*n)
!
!  Continuation.
!
  do k = 2, n

    lk = 2 * n - k + 1

    do l = k, 2 * n - k + 1
!
!  The quantities s2(l) for K < L are auxiliary quantities which may
!  be zero or may become so small as to underflow, without however
!  causing any harm.
!
      s2(l) = s1(l+1) - ( alpha(k-1) - a(l) ) * s1(l) &
        - beta(k-1) * s0(l) + b(l) * s1(l-1)

      if ( l == k ) then
        s(k) = s2(k)
      end if
!
!  Check impending underflow or overflow.
!
      if ( abs ( s(k) ) < 10.0E+00 * tiny ( s(k) ) ) then
        ierr = - ( k - 1 )
        return
      end if

      if ( 0.1E+00 * huge ( s(k) ) < abs ( s(k) ) ) then
        ierr = k - 1
        return
      end if

    end do
!
!  Compute the alpha and beta coefficients.
!
    alpha(k) = a(k) + ( s2(k+1) / s2(k))- ( s1(k) / s1(k-1) )
    beta(k) = s2(k) / s1(k-1)

    s0(k:lk) = s1(k:lk)
    s1(k:lk) = s2(k:lk)

  end do

  return
end
subroutine cheb_r8 ( n, da, db, dnu, dalpha, dbeta, ds, ierr )

!*****************************************************************************80
!
!! CHEB_R8 implements the modified Chebyshev algorithm.
!
!  Discussion:
!
!    The routine generates recursion coefficients ALPHA and BETA.
!
!    Given a set of polynomials  p(0),p(1),...,p(2*n-1)  satisfying
!
!      p(k+1)(x) = (x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
!      k = 0,1,...,2*n-2,
!
!      p(-1)(x) = 0,  p(0)(x)=1,
!
!    and associated modified moments
!
!      fnu(k) = integral of p(k)(x)*dlambda(x),
!      k = 0,1,...,2*n-1,
!
!    this subroutine uses the modified Chebyshev algorithm to generate the
!    recursion coefficients  alpha(k),beta(k), k = 0,1,...,n-1, for the
!    polynomials  pi(k)  orthogonal with respect to the integration
!    measure  dlambda(x), i.e.,
!
!      pi(k+1)(x) = (x-alpha(k))*pi(k)(x)-beta(k)*pi(k-1)(x),
!      k = 0,1,...,n-1,
!
!      pi(-1)(x) = 0,  pi(0)(x)=1.
!
!    On machines with limited exponent range, the occurrence of underflow
!    [overflow] in the computation of the  alpha's  and  beta's  can often
!    be avoided by multiplying all modified moments by a sufficiently large
!    [small] scaling factor and dividing the new  beta(0)  by the same
!    scaling factor.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    On Generating Orthogonal Polynomials,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 3, Number 3, 1982, pages 289-317.
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of recursion coefficients desired.
!
!    Input, real A(2*N-1), B(2*N-1), the values of A(k-1), b(k-1),
!    k = 1,2,...,2*n-1.
!
!    Input, real FNU(2*N), the values of the modified moments fnu(k-1),
!    k = 1,2,...,2*n
!
!    Output, real ALPHA(N), BETA(N), the
!    recursion coefficients  alpha(k-1),beta(k-1),
!    k = 1,2,...,n, where  beta(0)  is the total mass.
!
!    Output, real S(N), the normalization factors
!    s(k) = integral [pi(k)(x)]**2 dlambda(x), k=0,1,
!    2,...,n-1.
!
!    Output, integer IERR, an error flag.
!    0, on normal return,
!    1, if  abs(fnu(0))  is less than the machine zero,
!    2, if  n  is out of range,
!    -k  if S(K), k = 0,1,2,...,n-1, is about to underflow,
!    +k  if S(K) is about to overflow.
!
  implicit none

  integer n

  real ( kind = 8 ) da(2*n-1)
  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) db(2*n-1)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dnu(2*n)
  real ( kind = 8 ) ds(n)
  real ( kind = 8 ) ds0(2*n)
  real ( kind = 8 ) ds1(2*n)
  real ( kind = 8 ) ds2(2*n)
  integer ierr
  integer k
  integer l
  integer lk

  ierr = 0

  if ( abs ( dnu(1) ) < 10.0D+00 * tiny ( dnu(1) ) ) then
    ierr = 1
    return
  end if

  if ( n < 1 ) then
    ierr = 2
    return
  end if

  dalpha(1) = da(1) + dnu(2) / dnu(1)
  dbeta(1) = dnu(1)

  if ( n == 1 ) then
    return
  end if

  ds(1) = dnu(1)
  ds0(1:2*n) = 0.0D+00
  ds1(1:2*n) = dnu(1:2*n)

  do k = 2, n

    lk = 2 * n - k + 1

    do l = k, 2 * n - k + 1

      ds2(l) = ds1(l+1) - ( dalpha(k-1) - da(l) ) * ds1(l) - dbeta(k-1) * ds0(l) &
        + db(l) * ds1(l-1)

      if ( l == k ) then
        ds(k) = ds2(k)
      end if

    end do

    if ( abs ( ds(k) ) < 10.0D+00 * tiny ( ds(k) ) ) then
      ierr = -( k - 1 )
      return
    end if

    if ( 0.1D+00 * huge ( ds(k) ) < abs ( ds(k) ) ) then
      ierr = k - 1
      return
    end if

    dalpha(k) = da(k) + ( ds2(k+1) / ds2(k) ) - ( ds1(k) / ds1(k-1) )
    dbeta(k) = ds2(k) / ds1(k-1)

    ds0(k:lk) = ds1(k:lk)
    ds1(k:lk) = ds2(k:lk)

  end do

  return
end
subroutine chri_r4 ( n, iopt, a, b, x, y, hr, hi, alpha, beta, ierr )

!*****************************************************************************80
!
!! CHRI_R4 implements the Christoffel or generalized Christoffel theorem.
!
!  Discussion:
!
!    In all cases except  iopt = 7, it uses nonlinear recurrence
!    algorithms described in W. Gautschi,An algorithmic implementation
!    of the generalized Christoffel theorem'', Numerical Integration
!    (G. Haemmerlin, ed.), Birkhaeuser, Basel, 1982, pp. 89-106. The case
!    iopt = 7  incorporates a QR step with shift  x  in the manner of
!    J. Kautsky and G.H. Golub, On the calculation of Jacobi matrices'',
!    Linear Algebra Appl. 52/53, 1983, 439-455, using the algorithm of
!    Eq. (67.11) on p. 567 in J.H. Wilkinson,The Algebraic Eigenvalue
!    Problem'', Clarendon Press, Oxford, 1965. Given the recursion
!    coefficients  a(k),b(k), k = 0,1,...,n, for the (monic) orthogonal
!    polynomials with respect to some measure  dlambda(t), it generates
!    the recursion coefficients  alpha(k),beta(k), k = 0,1,...,n-1, for the
!    measure
!
!      (t-x)dlambda(t)               if  iopt = 1
!      [(t-x)**2+y**2]dlambda(t)     if  iopt = 2
!      (t**2+y**2)dlambda(t) with    if  iopt = 3
!      dlambda(t) and supp(dlambda) symmetric  with respect to
!      the origin
!      dlambda(t)/(t-x)              if  iopt = 4
!      dlambda(t)/[(t-x)**2+y**2]    if  iopt = 5
!      dlambda(t)/(t**2+y**2) with   if  iopt = 6
!      dlambda(t) and supp(dlambda) symmetric with respect to
!      the origin
!      [(t-x)**2]dlambda(t)          if  iopt = 7
!
!    It is assumed that  n  is larger than or equal to 2. Otherwise, the
!    routine exits immediately with the error flag  ierr  set equal to 1.
!    If  iopt  is not between 1 and 7, the routine exits with  ierr = 2.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of recurrence coefficients desired.
!
!    Input, integer IOPT, the desired weight distribution.
!
!    Input, real A(N+1), B(N+1), the recursion coefficients
!    a(k-1),b(k-1),k = 1,2,...,n+1, of the polynomials orthogonal with
!    respect to the given measure  dlambda(t)
!
!    Input, real X, Y, the linear and quadratic factors, or divisors,
!    of dlambda(t).
!
!    Input, real HR, HI, the real and imaginary part, respectively, of
!    the integral of dlambda(t)/(z-t), where z = x+iy;
!    the parameter  hr  is used only if  iopt = 4 or
!    5, the parameter  hi  only if  iopt = 5 or 6
!
!    Output, real ALPHA(N), BETA(N), the desired recursion coefficients
!    alpha(k-1), beta(k-1), k = 1,2,...,n
!
!    Output, integer IERR, an error flag.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer n

  real a(n+1)
  real alpha(n)
  real b(n+1)
  real beta(n)
  real c
  real c0
  real cm1
  real d
  real e
  real ei
  real eio
  real eioo
  real eo
  real eoo
  real eps
  real er
  real ero
  real eroo
  real gamma
  real hi
  real hr
  integer ierr
  integer iopt
  integer k
  real p2
  real q
  real s
  real so
  real t
  real u
  real x
  real y

  eps = 5.0E+00 * epsilon ( eps )
  ierr = 0

  if ( n < 2 ) then
    ierr = 1
    return
  end if
!
!  What follows implements Eq. (3.7) of W. Gautschi, op. cit.
!
  if ( iopt == 1 ) then

    e = 0.0E+00
    do k = 1, n
      q = a(k) - e - x
      beta(k) = q * e
      e = b(k+1) / q
      alpha(k) = x + q + e
    end do
!
!  Set the first beta-coefficient as discussed in Section 5.1 of the
!  companion paper.
!
    beta(1) = b(1) * ( a(1) - x )
!
!  What follows implements Eq. (4.7) of W. Gautschi, op. cit.
!
  else if ( iopt == 2 ) then

    s = x - a(1)
    t = y
    eio = 0.0E+00

    do k = 1, n
      d = s * s + t * t
      er = -b(k+1) * s / d
      ei = b(k+1) * t / d
      s = x + er - a(k+1)
      t = y + ei
      alpha(k) = x + t * er / ei - s * ei / t
      beta(k) = t * eio * ( 1.0E+00 + ( er / ei )**2 )
      eio = ei
    end do
!
!  Set the first beta-coefficient.
!
    beta(1) = b(1) * ( b(2) + ( a(1) - x )**2 + y * y )
!
!  What follows implements Eq. (4.8) of W. Gautschi, op. cit.
!
  else if ( iopt == 3 ) then

    t = y
    eio = 0.0E+00
    do k = 1, n
      ei = b(k+1) / t
      t = y + ei
      alpha(k) = 0.0E+00
      beta(k) = t * eio
      eio = ei
    end do
!
!  Set the first beta-coefficient.
!
    beta(1) = b(1) * ( b(2) + y * y )
!
!  What follows implements Eqs. (5.1),(5.2) of W. Gautschi, op. cit.
!
  else if ( iopt == 4 ) then

    alpha(1) = x - b(1) / hr
    beta(1) = -hr
    q = -b(1) / hr

    do k = 2, n
      e = a(k-1) - x - q
      beta(k) = q * e
      q = b(k) / e
      alpha(k) = q + e + x
    end do
!
!  What follows implements Eq. (5.8) of W. Gautschi, op. cit.
!
  else if ( iopt == 5 ) then

    d = hr * hr + hi * hi
    eroo = a(1) - x + b(1) * hr / d
    eioo = - b(1) * hi / d - y
    alpha(1) = x + hr * y / hi
    beta(1) = - hi / y
    alpha(2) = x - b(1) * hi * eroo / ( d * eioo ) + hr * eioo / hi
    beta(2) = y * eioo * ( 1.0E+00 + ( hr / hi )**2 )

    if ( n == 2 ) then
      return
    end if

    so = b(2) / ( eroo**2 + eioo**2 )
    ero = a(2) - x - so * eroo
    eio = so * eioo - y
    alpha(3) = x + eroo * eio / eioo + so * eioo * ero / eio
    beta(3) = - b(1) * hi * eio * ( 1.0E+00 + ( eroo / eioo )**2 ) / d

    do k = 3, n - 1
      s = b(k) / ( ero**2 + eio**2 )
      er = a(k) - x - s * ero
      ei = s * eio - y
      alpha(k+1) = x + ero * ei / eio + s * eio * er / ei
      beta(k+1) = so * eioo * ei * ( 1.0E+00 + ( ero / eio )**2 )
      eroo = ero
      eioo = eio
      ero = er
      eio = ei
      so = s
    end do
!
!  What follows implements Eq. (5.9) of W. Gautschi, op. cit.
!
  else if ( iopt == 6 ) then

    eoo = - b(1) / hi - y
    eo = b(2) / eoo - y
    alpha(1) = 0.0E+00
    beta(1) = - hi / y
    alpha(2) = 0.0E+00
    beta(2) = y * eoo

    if ( n == 2 ) then
      return
    end if

    alpha(3) = 0.0E+00
    beta(3) = -b(1) * eo / hi

    do k = 3, n - 1
      e = b(k) / eo - y
      beta(k+1) = b(k-1) * e / eoo
      alpha(k+1) = 0.0E+00
      eoo = eo
      eo = e
    end do
!
!  What follows implements a QR step with shift  x.
!
  else if ( iopt == 7 ) then

    u = 0.0E+00
    c = 1.0E+00
    c0 = 0.0E+00

    do k = 1, n

      gamma = a(k) - x - u
      cm1 = c0
      c0 = c

      if ( eps < abs ( c0 ) ) then
        p2 = ( gamma**2 ) / c0
      else
        p2 = cm1 * b(k)
      end if

      if ( 1 < k ) then
        beta(k) = s * ( p2 + b(k+1) )
      end if

      s = b(k+1) / ( p2 + b(k+1) )
      c = p2 / ( p2 + b(k+1) )
      u = s * ( gamma + a(k+1) - x )
      alpha(k) = gamma + u + x

    end do

    beta(1) = b(1) * ( b(2) + ( x - a(1) )**2 )

  else

    ierr = 2

  end if

  return
end
subroutine chri_r8 ( n, iopt, da, db, dx, dy, dhr, dhi, dalpha, dbeta, ierr )

!*****************************************************************************80
!
!! CHRI_R8 implements the Christoffel or generalized Christoffel theorem.
!
!  Discussion:
!
!    In all cases except  iopt = 7, it uses nonlinear recurrence
!    algorithms described in W. Gautschi,An algorithmic implementation
!    of the generalized Christoffel theorem'', Numerical Integration
!    (G. Haemmerlin, ed.), Birkhaeuser, Basel, 1982, pp. 89-106. The case
!    iopt = 7  incorporates a QR step with shift  x  in the manner of
!    J. Kautsky and G.H. Golub, On the calculation of Jacobi matrices'',
!    Linear Algebra Appl. 52/53, 1983, 439-455, using the algorithm of
!    Eq. (67.11) on p. 567 in J.H. Wilkinson,The Algebraic Eigenvalue
!    Problem'', Clarendon Press, Oxford, 1965. Given the recursion
!    coefficients  a(k),b(k), k = 0,1,...,n, for the (monic) orthogonal
!    polynomials with respect to some measure  dlambda(t), it generates
!    the recursion coefficients  alpha(k),beta(k), k = 0,1,...,n-1, for the
!    measure
!
!      (t-x)dlambda(t)               if  iopt = 1
!      [(t-x)**2+y**2]dlambda(t)     if  iopt = 2
!      (t**2+y**2)dlambda(t) with    if  iopt = 3
!      dlambda(t) and supp(dlambda) symmetric  with respect to
!      the origin
!      dlambda(t)/(t-x)              if  iopt = 4
!      dlambda(t)/[(t-x)**2+y**2]    if  iopt = 5
!      dlambda(t)/(t**2+y**2) with   if  iopt = 6
!      dlambda(t) and supp(dlambda) symmetric with respect to
!      the origin
!      [(t-x)**2]dlambda(t)          if  iopt = 7
!
!    It is assumed that  n  is larger than or equal to 2. Otherwise, the
!    routine exits immediately with the error flag  ierr  set equal to 1.
!    If  iopt  is not between 1 and 7, the routine exits with  ierr = 2.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of recurrence coefficients desired.
!
!    Input, integer IOPT, the desired weight distribution.
!
!    Input, real A(N+1), B(N+1), the recursion coefficients
!    a(k-1),b(k-1),k = 1,2,...,n+1, of the polynomials orthogonal with
!    respect to the given measure  dlambda(t)
!
!    Input, real X, Y, the linear and quadratic factors, or divisors,
!    of dlambda(t).
!
!    Input, real HR, HI, the real and imaginary part, respectively, of
!    the integral of dlambda(t)/(z-t), where z = x+iy;
!    the parameter  hr  is used only if  iopt = 4 or
!    5, the parameter  hi  only if  iopt = 5 or 6
!
!    Output, real ALPHA(N), BETA(N), the desired recursion coefficients
!    alpha(k-1), beta(k-1), k = 1,2,...,n
!
!    Output, integer IERR, an error flag.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer n

  real ( kind = 8 ) da(n+1)
  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) db(n+1)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dc
  real ( kind = 8 ) dc0
  real ( kind = 8 ) dcm1
  real ( kind = 8 ) dd
  real ( kind = 8 ) de
  real ( kind = 8 ) dei
  real ( kind = 8 ) deio
  real ( kind = 8 ) deioo
  real ( kind = 8 ) deo
  real ( kind = 8 ) deoo
  real ( kind = 8 ) deps
  real ( kind = 8 ) der
  real ( kind = 8 ) dero
  real ( kind = 8 ) deroo
  real ( kind = 8 ) dgam
  real ( kind = 8 ) dhi
  real ( kind = 8 ) dhr
  real ( kind = 8 ) dp2
  real ( kind = 8 ) dq
  real ( kind = 8 ) ds
  real ( kind = 8 ) dso
  real ( kind = 8 ) dt
  real ( kind = 8 ) du
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ierr
  integer iopt
  integer k

  deps = 5.0D+00 * epsilon ( deps )
  ierr = 0

  if ( n < 2 ) then
    ierr = 1
    return
  end if

  if ( iopt == 1 ) then

    de = 0.0D+00
    do k = 1, n
      dq = da(k) - de - dx
      dbeta(k) = dq * de
      de = db(k+1) / dq
      dalpha(k) = dx + dq + de
    end do

    dbeta(1) = db(1) * ( da(1) - dx )

  else if ( iopt == 2 ) then

    ds = dx - da(1)
    dt = dy
    deio = 0.0D+00

    do k = 1, n
      dd = ds * ds + dt * dt
      der = - db(k+1) * ds / dd
      dei = db(k+1) * dt / dd
      ds = dx + der - da(k+1)
      dt = dy + dei
      dalpha(k) = dx + dt * der / dei - ds * dei / dt
      dbeta(k) = dt * deio * ( 1.0D+00 + ( der / dei )**2 )
      deio = dei
    end do

    dbeta(1) = db(1) * ( db(2) + ( da(1) - dx )**2 + dy * dy )

  else if ( iopt == 3 ) then

    dt = dy
    deio = 0.0D+00

    do k = 1, n
      dei = db(k+1) / dt
      dt = dy + dei
      dalpha(k) = 0.0D+00
      dbeta(k) = dt * deio
      deio = dei
    end do

    dbeta(1) = db(1) * ( db(2) + dy * dy )

  else if ( iopt == 4 ) then

    dalpha(1) = dx - db(1) / dhr
    dbeta(1) = - dhr
    dq = - db(1) / dhr

    do k = 2, n
      de = da(k-1) - dx - dq
      dbeta(k) = dq * de
      dq = db(k) / de
      dalpha(k) = dq + de + dx
    end do

  else if ( iopt == 5 ) then

    dd = dhr * dhr + dhi * dhi
    deroo = da(1) - dx + db(1) * dhr / dd
    deioo = - db(1) * dhi / dd - dy
    dalpha(1) = dx + dhr * dy / dhi
    dbeta(1) = - dhi / dy
    dalpha(2) = dx - db(1) * dhi * deroo / ( dd * deioo ) + dhr * deioo / dhi
    dbeta(2) = dy * deioo * ( 1.0D+00 + ( dhr / dhi )**2 )

    if ( n == 2 ) then
      return
    end if

    dso = db(2) / ( deroo**2 + deioo**2 )
    dero = da(2) - dx - dso * deroo
    deio = dso * deioo - dy
    dalpha(3) = dx + deroo * deio / deioo + dso * deioo * dero / deio
    dbeta(3) = - db(1) * dhi * deio * ( 1.0D+00 + ( deroo / deioo )**2 ) / dd

    do k = 3, n - 1
      ds = db(k) / ( dero**2 + deio**2 )
      der = da(k) - dx - ds * dero
      dei = ds * deio - dy
      dalpha(k+1) = dx + dero * dei / deio + ds * deio * der / dei
      dbeta(k+1) = dso * deioo * dei * ( 1.0D+00 + ( dero / deio )**2 )
      deroo = dero
      deioo = deio
      dero = der
      deio = dei
      dso = ds
    end do

  else if ( iopt == 6 ) then

    deoo = - db(1) / dhi - dy
    deo = db(2) / deoo - dy
    dalpha(1) = 0.0D+00
    dbeta(1) = - dhi / dy
    dalpha(2) = 0.0D+00
    dbeta(2) = dy * deoo

    if ( n == 2 ) then
      return
    end if

    dalpha(3) = 0.0D+00
    dbeta(3) = - db(1) * deo / dhi

    do k = 3, n - 1
      de = db(k) / deo - dy
      dbeta(k+1) = db(k-1) * de / deoo
      dalpha(k+1) = 0.0D+00
      deoo = deo
      deo = de
    end do

  else if ( iopt == 7 ) then

    du = 0.0D+00
    dc = 1.0D+00
    dc0 = 0.0D+00

    do k = 1, n

      dgam = da(k) - dx - du
      dcm1 = dc0
      dc0 = dc

      if ( deps < abs ( dc0 ) ) then
        dp2 = ( dgam**2 ) / dc0
      else
        dp2 = dcm1*db(k)
      end if

      if ( 1 < k ) then
        dbeta(k) = ds * ( dp2 + db(k+1) )
      end if

      ds = db(k+1) / ( dp2 + db(k+1) )
      dc = dp2 / ( dp2 + db(k+1) )
      du = ds * ( dgam + da(k+1) - dx )
      dalpha(k) = dgam + du + dx

    end do

    dbeta(1) = db(1) * ( db(2) + ( dx - da(1) )**2 )

  else

    ierr = 2

  end if

  return
end
subroutine fejer_r4 ( n, x, w )

!*****************************************************************************80
!
!! FEJER_R4 generates a Fejer quadrature rule.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of quadrature nodes.
!
!    Output, real X(N), W(N), the quadrature nodes and weights.
!    The nodes are listed in increasing order.
!
  implicit none

  integer n

  real c0
  real c1
  real c2
  integer k
  integer m
  integer nh
  integer np1h
  real pi
  real t
  real total
  real w(n)
  real x(n)

  pi = 4.0E+00 * atan ( 1.0E+00 )
  nh = n / 2
  np1h = ( n + 1 ) / 2

  do k = 1, nh
    x(n+1-k) = cos ( 0.5E+00 * real ( 2 * k - 1, kind = 4 ) * pi &
      / real ( n, kind = 4 ) )
    x(k) = -x(n+1-k)
  end do

  if ( 2 * nh /= n ) then
    x(np1h) = 0.0E+00
  end if

  do k = 1, np1h

    c1 = 1.0E+00
    c0 = 2.0E+00 * x(k) * x(k) - 1.0E+00
    t = 2.0E+00 * c0
    total = c0 / 3.0E+00

    do m = 2, nh
      c2 = c1
      c1 = c0
      c0 = t * c1 - c2
      total = total + c0 / real ( 4 * m * m - 1, kind = 4 )
    end do

    w(k) = 2.0E+00 * ( 1.0E+00 - 2.0E+00 * total ) / real ( n, kind = 4 )
    w(n+1-k) = w(k)

  end do

  return
end
subroutine fejer_r8 ( n, x, w )

!*****************************************************************************80
!
!! FEJER_R8 generates a Fejer quadrature rule.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of quadrature nodes.
!
!    Output, real ( kind = 8 ) X(N), W(N), the quadrature nodes and weights.
!    The nodes are listed in increasing order.
!
  implicit none

  integer n

  real ( kind = 8 ) c0
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  integer k
  integer m
  integer nh
  integer np1h
  real ( kind = 8 ) pi
  real ( kind = 8 ) t
  real ( kind = 8 ) total
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  pi = 4.0D+00 * atan ( 1.0D+00 )
  nh = n / 2
  np1h = ( n + 1 ) / 2
  do k = 1, nh
    x(n+1-k) = cos ( 0.5D+00 * real ( 2 * k - 1, kind = 8 ) * pi &
      / real ( n, kind = 8 ) )
    x(k) = -x(n+1-k)
  end do

  if ( 2 * nh /= n ) then
    x(np1h) = 0.0D+00
  end if

  do k = 1, np1h

    c1 = 1.0D+00
    c0 = 2.0D+00 * x(k) * x(k) - 1.0D+00
    t = 2.0D+00 * c0
    total = c0 / 3.0D+00

    do m = 2, nh
      c2 = c1
      c1 = c0
      c0 = t * c1 - c2
      total = total + c0 / real ( 4 * m * m - 1, kind = 8 )
    end do

    w(k) = 2.0D+00 * ( 1.0D+00 - 2.0D+00 * total ) / real ( n, kind = 8 )
    w(n+1-k) = w(k)

  end do

  return
end
function gamma_r4 ( x, ierr )

!*****************************************************************************80
!
!! GAMMA_R4 evaluates the gamma function for real positive X.
!
!  Discussion:
!
!    The function ALGA_R4 is used. In case of overflow, the
!    routine returns the largest machine-representable number and the
!    error flag  ierr = 2.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, real X, the argument of the Gamma function.
!
!    Output, integer IERR, an error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Output, real GAMMA_R4, the value of the Gamma function at X.
!
  implicit none

  real alga_r4
  real almach
  real gamma_r4
  integer ierr
  real t
  real x

  almach = log ( huge ( almach ) )
  ierr = 0
  t = alga_r4 ( x )

  if ( almach <= t ) then
    ierr = 2
    gamma_r4 = huge ( gamma_r4 )
  else
    gamma_r4 = exp ( t )
  end if

  return
end
function gamma_r8 ( x, ierr )

!*****************************************************************************80
!
!! GAMMA_R8 evaluates the gamma function for real positive argument.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function.
!
!    Output, integer IERR, an error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Output, real ( kind = 8 ) GAMMA_R4, the value of the Gamma
!   function at X.
!
  implicit none

  real ( kind = 8 ) alga_r8
  real ( kind = 8 ) almach
  real ( kind = 8 ) gamma_r8
  integer ierr
  real ( kind = 8 ) t
  real ( kind = 8 ) x

  almach = log ( huge ( almach ) )
  ierr = 0
  t = alga_r8 ( x )

  if ( almach <= t ) then
    ierr = 2
    gamma_r8 = huge ( gamma_r8 )
  else
    gamma_r8 = exp ( t )
  end if

  return
end
subroutine gauss_r4 ( n, alpha, beta, eps, zero, weight, ierr, e )

!*****************************************************************************80
!
!! GAUSS_R4 generates an N-point Gaussian quadrature formula.
!
!  Discussion:
!
!    Given N and a measure DLAMBDA, this routine generates the N-point
!    Gaussian quadrature formula:
!
!      Integral over supp(dlambda) of f(x) dlambda(x)
!
!        =  sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
!
!    The nodes are returned as  zero(k) = x(k) and the weights as
!    weight(k) = w(k), k=1,2,...,n. The user has to supply the recursion
!    coefficients  alpha(k), beta(k), k = 0,1,2,...,n-1, for the measure
!    dlambda. The routine computes the nodes as eigenvalues, and the
!    weights in term of the first component of the respective normalized
!    eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
!    It uses a translation and adaptation of the algol procedure  imtql2,
!    Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified
!    by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for
!    Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
!    routine  imtql2.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of points in the Gaussian quadrature formula.
!
!    Input, real ALPHA(N), BETA(N), the values of alpha(k-1), beta(k-1),
!    k = 1,2,...,n
!
!    Input, real EPS, the relative accuracy desired in the nodes and weights.
!
!    Output, real ZERO(N), the Gaussian nodes in increasing order.
!
!    Output, real WEIGHT(N), the Gaussian weights.
!
!    Output, integer IERR, an error flag equal to 0 on normal return,
!    equal to I if the QR algorithm does not converge within 30 iterations on
!    evaluating the i-th eigenvalue, equal to -1 if N is not in range,
!    and equal to -2 if one of the BETA's is negative.
!
!    ?, real E(N), ?
!
  implicit none

  integer n

  real alpha(n)
  real b
  real beta(n)
  real c
  real e(n)
  real eps
  real f
  real g
  integer i
  integer ierr
  integer ii
  integer j
  integer k
  integer l
  integer m
  integer m2
  integer mml
  real p
  real r
  real s
  real weight(n)
  real zero(n)

  if ( n < 1 ) then
    ierr = -1
    return
  end if

  ierr = 0
  zero(1) = alpha(1)

  if ( beta(1) < 0.0E+00 ) then
    ierr = -2
    return
  end if

  weight(1) = beta(1)

  if ( n == 1 ) then
    return
  end if

  weight(1) = 1.0E+00
  e(n) = 0.0E+00

  do k = 2, n

    zero(k) = alpha(k)

    if ( beta(k) < 0.0E+00 ) then
      ierr = -2
      return
    end if

    e(k-1) = sqrt ( beta(k) )
    weight(k) = 0.0E+00

  end do

  do l = 1, n

    j = 0
!
!  Look for a small subdiagonal element.
!
    do

      do m2 = l, n

        m = m2

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= eps * ( abs ( zero(m) ) + abs ( zero(m+1) ) ) ) then
          exit
        end if

      end do

      p = zero(l)

      if ( m == l ) then
        exit
      end if

      if ( 30 <= j ) then
        ierr = l
        return
      end if

      j = j + 1
!
!  Form shift.
!
      g = ( zero(l+1) - p ) / ( 2.0E+00 * e(l) )
      r = sqrt ( g * g + 1.0E+00 )
      g = zero(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0E+00
      c = 1.0E+00
      p = 0.0E+00
      mml = m - l
!
!  For i = m-1 step -1 until l do ...
!
      do ii = 1, m - l

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r = sqrt ( c * c + 1.0E+00 )
          e(i+1) = f * r
          s = 1.0E+00 / r
          c = c * s
        else
          s = f / g
          r = sqrt ( s * s + 1.0E+00 )
          e(i+1) = g * r
          c = 1.0E+00 / r
          s = s * c
        end if

        g = zero(i+1) - p
        r = ( zero(i) - g ) * s + 2.0E+00 * c * b
        p = s * r
        zero(i+1) = g + p
        g = c * r - b
!
!  Form first component of vector.
!
        f = weight(i+1)
        weight(i+1) = s * weight(i) + c * f
        weight(i) = c * weight(i) - s * f

      end do

      zero(l) = zero(l) - p
      e(l) = g
      e(m) = 0.0E+00

    end do

  end do
!
!  Order the eigenvalues and eigenvectors.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = zero(i)

    do j = ii, n
      if ( zero(j) < p ) then
        k = j
        p = zero(j)
      end if
    end do

    if ( k /= i ) then
      zero(k) = zero(i)
      zero(i) = p
      p = weight(i)
      weight(i) = weight(k)
      weight(k) = p
    end if

  end do

  weight(1:n) = beta(1) * weight(1:n)**2

  return
end
subroutine gauss_r8 ( n, dalpha, dbeta, deps, dzero, dweigh, ierr, de )

!*****************************************************************************80
!
!! GAUSS_R8 generates an N-point Gaussian quadrature formula.
!
!  Discussion:
!
!    Given N and a measure DLAMBDA, this routine generates the N-point
!    Gaussian quadrature formula:
!
!      Integral over supp(dlambda) of f(x) dlambda(x)
!
!        =  sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
!
!    The nodes are returned as  zero(k) = x(k) and the weights as
!    weight(k) = w(k), k=1,2,...,n. The user has to supply the recursion
!    coefficients  alpha(k), beta(k), k = 0,1,2,...,n-1, for the measure
!    dlambda. The routine computes the nodes as eigenvalues, and the
!    weights in term of the first component of the respective normalized
!    eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
!    It uses a translation and adaptation of the algol procedure  imtql2,
!    Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified
!    by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for
!    Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
!    routine  imtql2.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of points in the Gaussian quadrature formula.
!
!    Input, real ALPHA(N), BETA(N), the values of alpha(k-1), beta(k-1),
!    k = 1,2,...,n
!
!    Input, real EPS, the relative accuracy desired in the nodes and weights.
!
!    Output, real ZERO(N), the Gaussian nodes in increasing order.
!
!    Output, real WEIGHT(N), the Gaussian weights.
!
!    Output, integer IERR, an error flag equal to 0 on normal return,
!    equal to I if the QR algorithm does not converge within 30 iterations on
!    evaluating the i-th eigenvalue, equal to -1 if N is not in range,
!    and equal to -2 if one of the BETA's is negative.
!
!    ?. real ( kind = 8 ) DE(N), ?
!
  implicit  none

  integer n

  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) db
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dc
  real ( kind = 8 ) de(n)
  real ( kind = 8 ) deps
  real ( kind = 8 ) df
  real ( kind = 8 ) dg
  real ( kind = 8 ) dp
  real ( kind = 8 ) dr
  real ( kind = 8 ) ds
  real ( kind = 8 ) dweigh(n)
  real ( kind = 8 ) dzero(n)
  integer i
  integer ierr
  integer j
  integer k
  integer l
  integer m
  integer m2
  integer mml

  if ( n < 1 ) then
    ierr = -1
    return
  end if

  ierr = 0

  dzero(1) = dalpha(1)

  if ( dbeta(1) < 0.0D+00 ) then
    ierr = -2
    return
  end if

  dweigh(1) = dbeta(1)

  if ( n == 1 ) then
    return
  end if

  dweigh(1) = 1.0D+00
  de(n) = 0.0D+00

  do k = 2, n

    dzero(k) = dalpha(k)

    if ( dbeta(k) < 0.0D+00 ) then
      ierr = -2
      return
    end if

    de(k-1) = sqrt ( dbeta(k) )
    dweigh(k) = 0.0D+00

  end do

  do l = 1, n

    j = 0

    do

      do m2 = l, n

        m = m2

        if ( m == n ) then
          exit
        end if

        if ( abs ( de(m) ) <= deps &
          * ( abs ( dzero(m) ) + abs ( dzero(m+1) ) ) ) then
          exit
        end if

      end do

      dp = dzero(l)

      if ( m == l ) then
        exit
      end if

      if ( 30 <= j ) then
        ierr = 1
        return
      end if

      j = j + 1
      dg = ( dzero(l+1) - dp ) / ( 2.0D+00 * de(l) )
      dr = sqrt ( dg * dg + 1.0D+00 )
      dg = dzero(m) - dp + de(l) / ( dg + sign ( dr, dg ) )
      ds = 1.0D+00
      dc = 1.0D+00
      dp = 0.0D+00
      mml = m - l

      do i = m - 1, 1, -1

        df = ds * de(i)
        db = dc * de(i)

        if ( abs ( df ) < abs ( dg ) ) then

          ds = df / dg
          dr = sqrt ( ds * ds + 1.0D+00 )
          de(i+1) = dg * dr
          dc = 1.0D+00 / dr
          ds = ds * dc

        else

          dc = dg / df
          dr = sqrt ( dc * dc + 1.0D+00 )
          de(i+1) = df * dr
          ds = 1.0D+00 / dr
          dc = dc * ds

        end if

        dg = dzero(i+1) - dp
        dr = ( dzero(i) - dg ) * ds + 2.0D+00 * dc * db
        dp = ds * dr
        dzero(i+1) = dg + dp
        dg = dc * dr - db
        df = dweigh(i+1)
        dweigh(i+1) = ds * dweigh(i) + dc * df
        dweigh(i) = dc * dweigh(i) - ds * df

      end do

      dzero(l) = dzero(l) - dp
      de(l) = dg
      de(m) = 0.0D+00

    end do

  end do

  do i = 1, n - 1

    k = i
    dp = dzero(i)

    do j = i + 1, n
      if ( dzero(j) < dp ) then
        k = j
        dp = dzero(j)
      end if
    end do

    if ( k /= i ) then
      dzero(k) = dzero(i)
      dzero(i) = dp
      dp = dweigh(i)
      dweigh(i) = dweigh(k)
      dweigh(k) = dp
    end if

  end do

  dweigh(1:n) = dbeta(1) * dweigh(1:n)**2

  return
end
subroutine gchri_r4 ( n, iopt, nu0, numax, eps, a, b, x, y, alpha, beta, nu, &
  ierr, ierrc, rho, rold )

!*****************************************************************************80
!
!! GCHRI_R4 implements the generalized Christoffel theorem.
!
!  Discussion:
!
!    The routine uses the method of modified moments (cf. Section 4 of
!    W. Gautschi, Minimal solutions of three-term recurrence relations and
!    orthogonal polynomials'', Math. Comp. 36, 1981, 547-554).
!
!    Given the recursion coefficients  a(k), b(k), k = 0,1,...n, for the monic
!    orthogonal polynomials with respect to some measure  dlambda(t), it
!    generates the recursion coefficients  alpha(k), beta(k), k = 0,1,2,...,n-1
!    for the measure
!
!      dlambda(t)/(t-x)        if iopt = 1
!      dlambda(t)/{(t-x)**2+y**2} if iopt = 2
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of recurrence coefficients desired.
!
!    Input, integer IOPT, the desired weight distribution.
!
!    Input, integer NU0, an estimate of the starting backward recurrence
!    index; in the absence of any better choice, take NU0 = 3 * N.
!
!    Input, integer NUMAX, controls the termination of the backward
!    recursion in case of nonconvergence; a conservative choice is
!    NUMAX = 500.
!
!    Input, real EPS, a relative error tolerance.
!
!           a,b - - arrays of dimension numax to be supplied with the
!                   recursion coefficients a(k) = alpha(k-1),b(k)=beta(k),
!                   k = 1,2,...,numax, for the measure  dlambda
!
!           x,y - - real parameters defining the linear and quadratic
!                   divisors of  dlambda
!
!   Output: alpha,beta - arrays of dimension  n  containing the desired
!                   recursion coefficients  alpha(k-1), beta(k-1), k = 1,
!                   2,...,n
!
!           nu  - - the backward recurrence index yielding convergence;
!                   in case of nonconvergence,  nu  will have the value
!                   numax
!
!           ierr  - an error flag, where
!                   ierr = 0     on normal return
!                   ierr = 1     if  iopt  is neither 1 nor 2
!                   ierr = nu0   if  numax < nu0
!                   ierr = numax if the backward recurrence algorithm does
!                              not converge
!                   ierr = -1    if  n  is not in range
!
!    Output, integer IERRC, an error flag from CHEB_R4.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer n
  integer numax

  real a(numax)
  real alpha(n)
  real b(numax)
  real beta(n)
  real eps
  real fnu(2*n)
  integer ierr
  integer ierrc
  integer iopt
  integer k
  integer nu
  integer nu0
  complex rho(2*n)
  real rold
  real s(n)
  real x
  real y
  complex z

  if ( n < 1 ) then
    ierr = -1
    return
  end if

  ierr = 0
!
!  Linear divisor.
!
  if ( iopt == 1 ) then
!
!  Generate the modified moments of dlambda.
!
    z = cmplx ( x )
    call knum_r4 ( 2 * n - 1, nu0, numax, z, eps, a, b, rho, nu, ierr )

    fnu(1:2*n) = - real ( rho(1:2*n), kind = 4 )
!
!  Compute the desired recursion coefficients by means of the modified
!  Chebyshev algorithm.
!
    call cheb_r4 ( n, a, b, fnu, alpha, beta, s, ierrc )
!
!  Quadratic divisor
!
  else if ( iopt == 2 ) then
!
!  Generate the modified moments of  dlambda.
!
    y = abs ( y )
    z = cmplx ( x, y )

    call knum_r4 ( 2 * n - 1, nu0, numax, z, eps, a, b, rho, nu, ierr )

    fnu(1:2*n) = - aimag ( rho(1:2*n) ) / y
!
!  Compute the desired recursion coefficients by means of the modified
!  Chebyshev algorithm.
!
    call cheb_r4 ( n, a, b, fnu, alpha, beta, s, ierrc )

  else

    ierr = 1

  end if

  return
end
subroutine gchri_r8 ( n, iopt, nu0, numax, deps, da, db, dx, dy, dalpha, &
  dbeta, nu, ierr, ierrc, drhor, drhoi, droldr, droldi )

!*****************************************************************************80
!
!! GCHRI_R8 implements the generalized Christoffel theorem.
!
!  Discussion:
!
!    The routine uses the method of modified moments (cf. Section 4 of
!    W. Gautschi, Minimal solutions of three-term recurrence relations and
!    orthogonal polynomials'', Math. Comp. 36, 1981, 547-554).
!
!    Given the recursion coefficients  a(k), b(k), k = 0,1,...n, for the monic
!    orthogonal polynomials with respect to some measure  dlambda(t), it
!    generates the recursion coefficients  alpha(k), beta(k), k = 0,1,2,...,n-1
!    for the measure
!
!      dlambda(t)/(t-x)        if iopt = 1
!      dlambda(t)/{(t-x)**2+y**2} if iopt = 2
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of recurrence coefficients desired.
!
!    Input, integer IOPT, selects the weight distribution.
!
!    Input, integer NU0, estimates the starting backward recurrence index.
!    If no good estimate is available, set NU0 = 3 * N.
!
!    Input, integer NUMAX, controls the termination of the backward
!    recursion in case of nonconvergence; a conservative choice is
!    NUMAX = 500.
!
!    Input, real EPS, a relative error tolerance.
!
!   Input:
!
!           a,b - - arrays of dimension numax to be supplied with the
!                   recursion coefficients a(k) = alpha(k-1),b(k)=beta(k),
!                   k = 1,2,...,numax, for the measure  dlambda
!           x,y - - real parameters defining the linear and quadratic
!                   divisors of  dlambda
!
!   Output: alpha,beta - arrays of dimension  n  containing the desired
!                   recursion coefficients  alpha(k-1), beta(k-1), k = 1,
!                   2,...,n
!           nu  - - the backward recurrence index yielding convergence;
!                   in case of nonconvergence,  nu  will have the value
!                   numax
!           ierr  - an error flag, where
!                   ierr = 0     on normal return
!                   ierr = 1     if  iopt  is neither 1 nor 2
!                   ierr = nu0   if  numax < nu0.
!                   ierr = numax if the backward recurrence algorithm does
!                              not converge
!                   ierr = -1    if  n  is not in range
!
!    Output, integer IERRC, an error flag from CHEB_R8.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer n
  integer numax

  real ( kind = 8 ) da(numax)
  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) db(numax)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) deps
  real ( kind = 8 ) dnu(2*n)
  real ( kind = 8 ) drhoi(2*n)
  real ( kind = 8 ) drhor(2*n)
  real ( kind = 8 ) droldi
  real ( kind = 8 ) droldr
  real ( kind = 8 ) ds(n)
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ierr
  integer ierrc
  integer iopt
  integer nu
  integer nu0

  if ( n < 1 ) then
    ierr = -1
    return
  end if

  ierr = 0

  if ( iopt == 1 ) then

    call knum_r8 ( 2 * n - 1, nu0, numax, dx, dy, deps, da, db, drhor, drhoi, &
      nu, ierr )

    dnu(1:2*n) = - drhor(1:2*n)

    call cheb_r8 ( n, da, db, dnu, dalpha, dbeta, ds, ierrc )

  else if ( iopt == 2 ) then

    dy = abs ( dy )

    call knum_r8 ( 2 * n - 1, nu0, numax, dx, dy, deps, da, db, drhor, drhoi, &
      nu, ierr )

    dnu(1:2*n) = - drhoi(1:2*n) / dy

    call cheb_r8 ( n, da, db, dnu, dalpha, dbeta, ds, ierrc )

  else

    ierr = 1

  end if

  return
end
subroutine kern_r4 ( n, nu0, numax, z, eps, a, b, ker, nu, ierr )

!*****************************************************************************80
!
!! KERN_R4 generates the kernels in the Gauss quadrature remainder term.
!
!  Discussion:
!
!    The kernels to be generated are
!
!      K(k)(z) = rho(k)(z)/pi(k)(z), k=0,1,2,...,n,
!
!    where rho(k) are the output quantities of the routine KNUM_R4, and
!    pi(k) the monic orthogonal polynomials.  The results are returned
!    in the array ker as ker(k) = K(k-1)(z), k=1,2,...,n+1.  All the other
!    input and output parameters have the same meaning as in the routine
!    knum.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, indicates the number of kernels to be computed.
!
!    Output, integer IERR, an error flag from KNUM_R4.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer n
  integer numax

  real a(numax)
  real b(numax)
  real eps
  integer ierr
  integer k
  complex ker(n+1)
  integer nu
  integer nu0
  complex p
  complex p0
  complex pm1
  complex z

  call knum_r4 ( n, nu0, numax, z, eps, a, b, ker, nu, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  p0 = cmplx ( 0.0E+00, 0.0E+00 )
  p = cmplx ( 1.0E+00, 0.0E+00 )

  do k = 1, n
    pm1 = p0
    p0 = p
    p = ( z - a(k) ) * p0 - b(k) * pm1
    ker(k+1) = ker(k+1) / p
  end do

  return
end
subroutine kern_r8 ( n, nu0, numax, dx, dy, deps, da, db, dkerr, dkeri, nu, &
  ierr )

!*****************************************************************************80
!
!! KERN_R8 generates the kernels in the Gauss quadrature remainder term.
!
!  Discussion:
!
!    This routine was written under the assumption that double precision complex
!    arithmetic was NOT available.
!
!    The kernels to be generated are
!
!      K(k)(z) = rho(k)(z)/pi(k)(z), k=0,1,2,...,n,
!
!    where rho(k) are the output quantities of the routine KNUM_R8, and
!    pi(k) the (monic) orthogonal polynomials.  The results are returned
!    in the array ker as ker(k) = K(k-1)(z), k=1,2,...,n+1.  All the other
!    input and output parameters have the same meaning as in the routine
!    KNUM_R8.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, indicates the number of kernels to be computed.
!
!    Output, integer IERR, an error flag from KNUM_R8.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer n
  integer numax

  real ( kind = 8 ) da(numax)
  real ( kind = 8 ) db(numax)
  real ( kind = 8 ) dden
  real ( kind = 8 ) deps
  real ( kind = 8 ) dkeri(n+1)
  real ( kind = 8 ) dkerr(n+1)
  real ( kind = 8 ) dp0i
  real ( kind = 8 ) dp0r
  real ( kind = 8 ) dpi
  real ( kind = 8 ) dpm1i
  real ( kind = 8 ) dpm1r
  real ( kind = 8 ) dpr
  real ( kind = 8 ) dt
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ierr
  integer k
  integer nu
  integer nu0

  call knum_r8 ( n, nu0, numax, dx, dy, deps, da, db, dkerr, dkeri, nu, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  dp0r = 0.0D+00
  dp0i = 0.0D+00
  dpr = 1.0D+00
  dpi = 0.0D+00

  do k = 1, n

    dpm1r = dp0r
    dpm1i = dp0i
    dp0r = dpr
    dp0i = dpi

    dpr = ( dx - da(k) ) * dp0r - dy * dp0i - db(k) * dpm1r
    dpi = ( dx - da(k) ) * dp0i + dy * dp0r - db(k) * dpm1i

    dden = dpr**2 + dpi**2
    dt = ( dkerr(k+1) * dpr + dkeri(k+1) * dpi ) / dden

    dkeri(k+1) = ( dkeri(k+1) * dpr - dkerr(k+1) * dpi ) / dden
    dkerr(k+1) = dt

  end do

  return
end
subroutine knum_r4 ( n, nu0, numax, z, eps, a, b, rho, nu, ierr )

!*****************************************************************************80
!
!! KNUM_R4 integrates certain rational polynomials.
!
!  Discussion:
!
!    The routine generates
!
!      rho(k)(z) = integral pi(k)(t)dlambda(t)/(z-t), k=0,1,2,...,n,
!
!    where  pi(k)(t)  is the (monic) k-th degree orthogonal polynomial
!    with respect to the measure  dlambda(t), and the integral is extended
!    over the support of  dlambda. It is assumed that  z  is a complex
!    number outside the smallest interval containing the support of
!    dlambda. The quantities  rho(k)(z)  are computed as the first  n+1
!    members of the minimal solution of the basic three-term recurrence
!    relation
!
!      y(k+1)(z) = (z-a(k))y(k)(z)-b(k)y(k-1)(z), k=0,1,2,...,
!
!    satisfied by the orthogonal polynomials  pi(k)(z).
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, indicates that RHO(0:N) are desired.
!
!           nu0  -  an estimate of the starting backward recurrence
!                   index; if no better estimate is known, set
!                   nu0  =  3*n/2; for Jacobi, Laguerre and Hermite
!                   weight functions, estimates of  nu0  are generated
!                   respectively by the routines  nu0jac, nu0lag  and
!                   nu0her
!
!           numax - an integer larger than  n  cutting off backward
!                   recursion in case of nonconvergence; if  nu0
!                   exceeds  numax, then the routine aborts with the
!                   error flag  ierr  set equal to  nu0
!
!           z - - - the variable in  rho(k)(z); type complex
!
!           eps - - the relative accuracy to which the  rho(k)  are
!                   desired
!
!           a,b - - arrays of dimension  numax  to be supplied with the
!                   recurrence coefficients  a(k-1), b(k-1), k = 1,2,...,
!                   numax.
!
!   Output: rho - - an array of dimension  n+1  containing the results
!                   rho(k) = rho(k-1)(z), k=1,2,...,n+1; type complex
!
!           nu  - - the starting backward recurrence index that yields
!                   convergence
!
!           ierr  - an error flag equal to zero on normal return, equal
!                   to  nu0  if numax < nu0, and equal to  numax in
!                   case of nonconvergence.
!
  implicit none

  integer n
  integer numax

  real a(numax)
  real b(numax)
  logical done
  real eps
  integer ierr
  integer j
  integer k
  integer nu
  integer nu0
  complex r
  complex rho(n+1)
  complex rold(n+1)
  complex z

  ierr = 0

  if ( numax < nu0 ) then
    ierr = nu0
    return
  end if

  if ( nu0 < n + 1 ) then
    nu0 = n + 1
  end if

  nu = nu0 - 5

  rho(1:n+1) = cmplx ( 0.0E+00, 0.0E+00 )

  do

    nu = nu + 5

    if ( numax < nu ) then
      ierr = numax
      exit
    end if

    rold(1:n+1) = rho(1:n+1)

    r = 0.0E+00

    do j = 1, nu
      r = b(nu-j+1) / ( z - a(nu-j+1) - r )
      if ( nu - j + 1 <= n + 1 ) then
        rho(nu-j+1) = r
      end if
    end do

    done = .true.

    do k = 1, n + 1
      if ( eps * abs ( rho(k) ) < abs ( rho(k) - rold(k) ) ) then
        done = .false.
        exit
      end if
    end do

    if ( done ) then
      exit
    end if

  end do

  do k = 2, n + 1
    rho(k) = rho(k) * rho(k-1)
  end do

  return
end
subroutine knum_r8 ( n, nu0, numax, dx, dy, deps, da, db, drhor, drhoi, nu, &
  ierr )

!*****************************************************************************80
!
!! KNUM_R8 is a double-precision version of the routine KNUM_R4.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, indicates that RHO(0:N) are desired.
!
  implicit none

  integer n
  integer numax

  real ( kind = 8 ) da(numax)
  real ( kind = 8 ) db(numax)
  real ( kind = 8 ) dden
  real ( kind = 8 ) deps
  logical done
  real ( kind = 8 ) drhoi(n+1)
  real ( kind = 8 ) drhor(n+1)
  real ( kind = 8 ) droldi(n+1)
  real ( kind = 8 ) droldr(n+1)
  real ( kind = 8 ) dri
  real ( kind = 8 ) drr
  real ( kind = 8 ) dt
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ierr
  integer j
  integer j1
  integer k
  integer nu
  integer nu0

  ierr = 0

  if ( numax < nu0 ) then
    ierr = nu0
    return
  end if

  nu0 = max ( nu0, n + 1 )

  nu = nu0 - 5

  drhor(1:n+1) = 0.0D+00
  drhoi(1:n+1) = 0.0D+00

  do

    nu = nu + 5

    if ( numax < nu ) then
      ierr = numax
      exit
    end if

    droldr(1:n+1) = drhor(1:n+1)
    droldi(1:n+1) = drhoi(1:n+1)

    drr = 0.0D+00
    dri = 0.0D+00

    do j = 1, nu

      j1 = nu - j + 1
      dden = ( dx - da(j1) - drr )**2 + ( dy - dri )**2
      drr = db(j1) * ( dx - da(j1) - drr ) / dden
      dri = - db(j1) * ( dy - dri ) / dden

      if ( j1 <= n + 1 ) then
        drhor(j1) = drr
        drhoi(j1) = dri
      end if

    end do

    done = .true.

    do k = 1, n + 1

      if ( ( deps**2 ) * ( drhor(k)**2 + drhoi(k)**2 ) < &
        ( drhor(k) - droldr(k) )**2 + ( drhoi(k) - droldi(k) )**2 ) then
        done = .false.
        exit
      end if

    end do

    if ( done ) then
      exit
    end if

  end do

  do k = 2, n + 1
    dt = drhor(k) * drhor(k-1) - drhoi(k) * drhoi(k-1)
    drhoi(k) = drhor(k) * drhoi(k-1) + drhoi(k) * drhor(k-1)
    drhor(k) = dt
  end do

  return
end
subroutine lancz_r4 ( n, ncap, x, w, alpha, beta, ierr, p0, p1 )

!*****************************************************************************80
!
!! LANCZ_R4 applies Stieltjes's procedure, using the Lanczos method.
!
!  Discussion:
!
!    The routine carries out the same task as the routine STI_R4,
!    but uses the more stable Lanczos method. The meaning of the input
!    and output parameters is the same as in the routine STI_R4.
!
!    This routine is adapted from the routine RKPW in W.B. Gragg and
!    W.J. Harrod,The numerically stable reconstruction of Jacobi
!    matrices from spectral data'', Numer. Math. 44, 1984, 317-335.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Output, integer IERR, an error flag.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer n
  integer ncap

  real alpha(n)
  real beta(n)
  real gam
  integer i
  integer ierr
  integer k
  real p0(ncap)
  real p1(ncap)
  real pi
  real rho
  real sig
  real t
  real tk
  real tmp
  real tsig
  real w(ncap)
  real x(ncap)
  real xlam

  if ( n <= 0 .or. ncap < n ) then
    ierr = 1
    return
  end if

  ierr = 0

  p0(1:ncap) = x(1:ncap)

  p1(1) = w(1)
  p1(2:ncap) = 0.0E+00

  do i = 1, ncap - 1

    pi = w(i+1)
    gam = 1.0E+00
    sig = 0.0E+00
    t = 0.0E+00
    xlam = x(i+1)

    do k = 1, i + 1

      rho = p1(k) + pi
      tmp = gam * rho
      tsig = sig

      if ( rho <= 0.0E+00 ) then
        gam = 1.0E+00
        sig = 0.0E+00
      else
        gam = p1(k) / rho
        sig = pi / rho
      end if

      tk = sig * ( p0(k) - xlam ) - gam * t
      p0(k) = p0(k) - ( tk - t )
      t = tk

      if ( sig <= 0.0E+00 ) then
        pi = tsig * p1(k)
      else
        pi = ( t**2 ) / sig
      end if

      tsig = sig
      p1(k) = tmp

    end do
  end do

  alpha(1:n) = p0(1:n)
  beta(1:n) = p1(1:n)

  return
end
subroutine lancz_r8 ( n, ncap, dx, dw, dalpha, dbeta, ierr, dp0, dp1 )

!*****************************************************************************80
!
!! LANCZ_R8 is a double-precision version of the routine LANCZ_R4.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Output, integer IERR, an error flag.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer n
  integer ncap

  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dgam
  real ( kind = 8 ) dp0(ncap)
  real ( kind = 8 ) dp1(ncap)
  real ( kind = 8 ) dpi
  real ( kind = 8 ) drho
  real ( kind = 8 ) dsig
  real ( kind = 8 ) dt
  real ( kind = 8 ) dtk
  real ( kind = 8 ) dtmp
  real ( kind = 8 ) dtsig
  real ( kind = 8 ) dw(ncap)
  real ( kind = 8 ) dx(ncap)
  real ( kind = 8 ) dxlam
  integer i
  integer ierr
  integer k

  if ( n <= 0 .or. ncap < n ) then
    ierr = 1
    return
  end if

  ierr = 0

  dp0(1:ncap) = dx(1:ncap)
  dp1(1) = dw(1)
  dp1(2:ncap) = 0.0D+00

  do i = 1, ncap - 1

    dpi = dw(i+1)
    dgam = 1.0D+00
    dsig = 0.0D+00
    dt = 0.0D+00
    dxlam = dx(i+1)

    do k = 1, i + 1

      drho = dp1(k) + dpi
      dtmp = dgam * drho
      dtsig = dsig

      if ( drho <= 0.0D+00 ) then
        dgam = 1.0D+00
        dsig = 0.0D+00
      else
        dgam = dp1(k) / drho
        dsig = dpi / drho
      end if

      dtk = dsig * ( dp0(k) - dxlam ) - dgam * dt
      dp0(k) = dp0(k) - ( dtk - dt )
      dt = dtk

      if ( dsig <= 0.0D+00 ) then
        dpi = dtsig * dp1(k)
      else
        dpi = ( dt**2 ) / dsig
      end if

      dtsig = dsig
      dp1(k) = dtmp

    end do
  end do

  dalpha(1:n) = dp0(1:n)
  dbeta(1:n) = dp1(1:n)

  return
end
subroutine lob_r4 ( n, alpha, beta, left, right, zero, weight, ierr )

!*****************************************************************************80
!
!! LOB_R4 generates a Gauss-Lobatto quadrature rule.
!
!  Discussion:
!
!    Given N and a measure DLAMBDA, this routine generates the
!    (N+2)-point Gauss-Lobatto quadrature formula
!
!      Integral over support(DLAMBDA) of F(X) DLAMBDA(X)
!
!      = w(0) f(x(0)) + sum from k=1 to k=n of w(k) f(x(k))
!      + w(n+1) f(x(n+1)) + R(n;f).
!
!    The nodes are returned as  zero(k) = x(k), the weights as  weight(k)
!    = w(k), k=0,1,...,n,n+1.  The user has to supply the recursion
!    coefficients  alpha(k), beta(k), k = 0,1,...,n,n+1, for the measure
!    dlambda.  The nodes and weights are computed in terms of the
!    eigenvalues and first component of the normalized eigenvectors of
!    a slightly modified Jacobi matrix of order n+2.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of interior points in the Gauss-Lobatto
!    formula.
!
!    Input, real ALPHA(N+2), BETA(N+2), the recursion coefficients
!    alpha(k-1), beta(k-1), k = 1,2,...,n+2, of the underlying measure;
!    the routine does not use alpha(n+2), beta(n+2).
!
!    Input, real LEFT, RIGHT, the prescribed left and right endpoints
!    x(0) and x(n+1) of the Gauss-Lobatto formula
!
!    Output, real ZERO(N+2), the nodes in increasing order.  zero(k)
!    = x(k), k=0,1,...,n,n+1
!
!    Output, real WEIGHT(N+2), the weights
!    weight(k) = w(k), k=0,1,...,n,n+1
!
!    Output, integer IERR, an error flag from the routine GAUSS_R4.
!
  implicit none

  integer n

  real a(n+2)
  real aleft
  real alpha(n+2)
  real b(n+2)
  real beta(n+2)
  real det
  real e(n+2)
  real epsma
  integer ierr
  integer k
  real left
  real p0l
  real p0r
  real p1l
  real p1r
  real pm1l
  real pm1r
  real right
  real weight(n+2)
  real zero(n+2)

  epsma = epsilon ( epsma )

  a(1:n+1) = alpha(1:n+1)
  b(1:n+1) = beta(1:n+1)

  p0l = 0.0E+00
  p0r = 0.0E+00
  p1l = 1.0E+00
  p1r = 1.0E+00

  do k = 1, n + 1
    pm1l = p0l
    p0l = p1l
    pm1r = p0r
    p0r = p1r
    p1l = ( left - a(k) ) * p0l - b(k) * pm1l
    p1r = ( right - a(k) ) * p0r - b(k) * pm1r
  end do

  det = p1l * p0r - p1r * p0l

  a(n+2) = ( left * p1l * p0r - right * p1r * p0l ) / det
  b(n+2) = ( right - aleft ) * p1l * p1r / det

  call gauss_r4 ( n + 2, a, b, epsma, zero, weight, ierr, e )

  return
end
subroutine lob_r8 ( n, alpha, beta, left, right, zero, weight, ierr )

!*****************************************************************************80
!
!! LOB_R8 generates a Gauss-Lobatto quadrature rule.
!
!  Discussion:
!
!    Given N and a measure DLAMBDA, this routine generates the
!    (N+2)-point Gauss-Lobatto quadrature formula
!
!      Integral over support(DLAMBDA) of F(X) DLAMBDA(X)
!
!      = w(0) f(x(0)) + sum from k=1 to k=n of w(k) f(x(k))
!      + w(n+1) f(x(n+1)) + R(n;f).
!
!    The nodes are returned as  zero(k) = x(k), the weights as  weight(k)
!    = w(k), k=0,1,...,n,n+1.  The user has to supply the recursion
!    coefficients  alpha(k), beta(k), k = 0,1,...,n,n+1, for the measure
!    dlambda.  The nodes and weights are computed in terms of the
!    eigenvalues and first component of the normalized eigenvectors of
!    a slightly modified Jacobi matrix of order  n+2.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of interior points in the Gauss-Lobatto
!    formula.
!
!    Input, real ( kind = 8 ) ALPHA(N+2), BETA(N+2), the recursion coefficients
!    alpha(k-1), beta(k-1), k = 1,2,...,n+2, of the underlying measure;
!    the routine does not use alpha(n+2), beta(n+2).
!
!    Input, real ( kind = 8 ) LEFT, RIGHT, the prescribed left and right
!    endpoints x(0) and x(n+1) of the Gauss-Lobatto formula
!
!    Output, real ( kind = 8 ) ZERO(N+2), the nodes in increasing order.  zero(k)
!    = x(k), k=0,1,...,n,n+1
!
!    Output, real ( kind = 8 ) WEIGHT(N+2), the weights
!    weight(k) = w(k), k=0,1,...,n,n+1
!
!    Output, integer IERR, an error flag from the routine GAUSS_R8.
!
  implicit none

  integer n

  real ( kind = 8 ) a(n+2)
  real ( kind = 8 ) alpha(n+2)
  real ( kind = 8 ) b(n+2)
  real ( kind = 8 ) beta(n+2)
  real ( kind = 8 ) det
  real ( kind = 8 ) e(n+2)
  real ( kind = 8 ) epsma
  integer ierr
  integer k
  real ( kind = 8 ) left
  real ( kind = 8 ) p0l
  real ( kind = 8 ) p0r
  real ( kind = 8 ) p1l
  real ( kind = 8 ) p1r
  real ( kind = 8 ) pm1l
  real ( kind = 8 ) pm1r
  real ( kind = 8 ) right
  real ( kind = 8 ) weight(n+2)
  real ( kind = 8 ) zero(n+2)

  epsma = epsilon ( epsma )

  a(1:n+1) = alpha(1:n+1)
  b(1:n+1) = beta(1:n+1)

  p0l = 0.0D+00
  p0r = 0.0D+00
  p1l = 1.0D+00
  p1r = 1.0D+00

  do k = 1, n + 1
    pm1l = p0l
    p0l = p1l
    pm1r = p0r
    p0r = p1r
    p1l = ( left - a(k) ) * p0l - b(k) * pm1l
    p1r = ( right - a(k) ) * p0r - b(k) * pm1r
  end do

  det = p1l * p0r - p1r * p0l
  a(n+2) = ( left * p1l * p0r - right * p1r * p0l ) / det
  b(n+2) = ( right - left ) * p1l * p1r / det

  call gauss_r8 ( n + 2, a, b, epsma, zero, weight, ierr, e )

  return
end
subroutine mccheb_r4 ( n, ncapm, mc, mp, xp, yp, quad_r4, eps, iq, idelta, &
  finl, finr, endl, endr, xfer, wfer, a, b, fnu, alpha, beta, ncap, kount, &
  ierr )

!*****************************************************************************80
!
!! MCCHEB_R4 is a multiple-component discretized modified Chebyshev algorithm.
!
!  Discussion:
!
!    The routine is basically a modified Chebyshev algorithm in which the
!    modified moments are discretized in the same manner as the inner
!    product in the discretization procedure MCDIS_R4.
!
!    The input and output parameters are as in MCDIS_R4.
!
!    In addition, the arrays  a,b must be filled with the recursion
!    coefficients  a(k-1),b(k-1), k = 1,2,...,2*n-1, defining the modified
!    moments.
!
!    The routine exits immediately with ierr = -1  if  n  is not in range.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Output, integer IERR, an error flag.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer mc
  integer mp
  integer n
  integer ncapm

  real a(2*n-1)
  real alpha(n)
  real b(2*n-1)
  real be(n)
  real beta(n)
  logical done
  real endl(mc)
  real endr(mc)
  real eps
  logical finl
  logical finr
  real fnu(2*n)
  integer i
  integer idelta
  integer ierr
  integer im1tn
  integer incap
  integer iq
  integer k
  integer kount
  integer l
  integer mtncap
  integer mtnpmp
  integer ncap
  integer nd
  real p
  real p1
  real pm1
  external quad_r4
  real s(n)
  real w(ncapm)
  real wfer(ncapm)
  real wm(mc*ncapm+mp)
  real x(ncapm)
  real xfer(ncapm)
  real xm(mc*ncapm+mp)
  real xp(mp)
  real yp(mp)

  nd = 2 * n

  if ( idelta <= 0 ) then
    idelta = 1
  end if

  if ( n < 1 ) then
    ierr = -1
    return
  end if
!
!  Initialization
!
  incap = 1
  kount = -1
  ierr = 0
  beta(1:n) = 0.0E+00
  ncap = ( 2 * n - 1 ) / idelta

  do

    be(1:n) = beta(1:n)

    kount = kount + 1

    if ( 1 < kount ) then
      incap = 2**( kount / 5 ) * n
    end if

    ncap = ncap + incap

    if ( ncapm < ncap ) then
      ierr = ncapm
      return
    end if
!
!  Discretization of the modified moments
!
    mtncap = mc * ncap

    do i = 1, mc

      im1tn = ( i - 1 ) * ncap

      if ( iq == 1 ) then
        call quad_r4 ( ncap, x, w, i, ierr )
      else
        call qgp_r4 ( ncap, x, w, i, ierr, mc, finl, finr, endl, endr, &
          xfer, wfer )
      end if

      if ( ierr /= 0 ) then
        ierr = i
        return
      end if

      xm(im1tn+1:im1tn+ncap) = x(1:ncap)
      wm(im1tn+1:im1tn+ncap) = w(1:ncap)

    end do

    if ( mp /= 0 ) then
      xm(mtncap+1:mtncap+mp) = xp(1:mp)
      wm(mtncap+1:mtncap+mp) = yp(1:mp)
    end if

    mtnpmp = mtncap + mp

    do k = 1, 2 * n

      fnu(k) = 0.0E+00

      do i = 1, mtnpmp

        p1 = 0.0E+00
        p = 1.0E+00

        do l = 1, k - 1
          pm1 = p1
          p1 = p
          p = ( xm(i) - a(l) ) * p1 - b(l) * pm1
        end do

        fnu(k) = fnu(k) + wm(i) * p

      end do
    end do
!
!  Computation of the desired recursion coefficients
!
    call cheb_r4 ( n, a, b, fnu, alpha, beta, s, ierr )
!
!  In the following statement, the absolute value of the beta's is
!  used to guard against failure in cases where the routine is applied
!  to variable-sign weight functions and hence the positivity of the
!  beta's is not guaranteed.
!
    done = .true.

    do k = 1, n
      if ( eps * abs ( beta(k) ) < abs ( beta(k)-be(k) ) ) then
        done = .false.
        exit
      end if
    end do

    if ( done ) then
      exit
    end if

  end do

  return
end
subroutine mccheb_r8 ( n, ncapm, mc, mp, dxp, dyp, quad_r8, deps, iq, &
  idelta, finld, finrd, dendl, dendr, dxfer, dwfer, da, db, dnu, dalpha, &
  dbeta, ncap, kount, ierrd )

!*****************************************************************************80
!
!! MCCHEB_R8 is a double-precision version of the routine MCCHEB_R4.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Output, integer IERRD, an error flag.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer mc
  integer mp
  integer n
  integer ncapm

  real ( kind = 8 ) da(2*n-1)
  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) db(2*n-1)
  real ( kind = 8 ) dbe(n)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dendl(mc)
  real ( kind = 8 ) dendr(mc)
  real ( kind = 8 ) deps
  real ( kind = 8 ) dnu(2*n)
  logical done
  real ( kind = 8 ) dp
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dpm1
  real ( kind = 8 ) ds(n)
  real ( kind = 8 ) dsum
  real ( kind = 8 ) dw(ncapm)
  real ( kind = 8 ) dwfer(ncapm)
  real ( kind = 8 ) dwm(mc*ncapm+mp)
  real ( kind = 8 ) dx(ncapm)
  real ( kind = 8 ) dxfer(ncapm)
  real ( kind = 8 ) dxm(mc*ncapm+mp)
  real ( kind = 8 ) dxp(mp)
  real ( kind = 8 ) dyp(mp)
  logical finld
  logical finrd
  integer i
  integer idelta
  integer ierr
  integer ierrd
  integer im1tn
  integer incap
  integer iq
  integer k
  integer kount
  integer l
  integer mcd
  integer mtncap
  integer mtnpmp
  integer ncap
  integer nd
  external quad_r8

  nd = 2 * n

  if ( idelta <= 0 ) then
    idelta = 1
  end if

  if ( n < 1 ) then
    ierrd = -1
    return
  end if

  incap = 1
  kount = -1
  ierrd = 0
  dbeta(1:n) = 0.0D+00
  ncap = ( 2 * n - 1 ) / idelta

  do

    dbe(1:n) = dbeta(1:n)

    kount = kount + 1

    if ( 1 < kount ) then
      incap = 2**( kount / 5 ) * n
    end if

    ncap = ncap + incap

    if ( ncapm < ncap ) then
      ierrd = ncapm
      return
    end if

    mtncap = mcd * ncap

    do i = 1, mc

      im1tn = ( i - 1 ) * ncap
      if ( iq == 1 ) then
        call quad_r8 ( ncap, dx, dw, i, ierr )
      else
        call qgp_r8 ( ncap, dx, dw, i, ierr, mc, finld, finrd, dendl, dendr, &
          dxfer, dwfer )
      end if

      if ( ierr /= 0 ) then
        ierrd = i
        return
      end if

      dxm(im1tn+1:im1tn+ncap) = dx(1:ncap)
      dwm(im1tn+1:im1tn+ncap) = dw(1:ncap)

    end do

    dxm(mtncap+1:mtncap+mp) = dxp(1:mp)
    dwm(mtncap+1:mtncap+mp) = dyp(1:mp)

    mtnpmp = mtncap + mp

    do k = 1, 2 * n

      dsum = 0.0D+00

      do i = 1, mtnpmp

        dp1 = 0.0D+00
        dp = 1.0D+00

        if ( 1 < k ) then

          do l = 1, k - 1
            dpm1 = dp1
            dp1 = dp
            dp = ( dxm(i) - da(l) ) * dp1 - db(l) * dpm1
          end do

        end if

        dsum = dsum + dwm(i) * dp

      end do

      dnu(k) = dsum

    end do

    call cheb_r8 ( n, da, db, dnu, dalpha, dbeta, ds, ierr )

    done = .true.

    do k = 1, n
      if ( deps * abs ( dbeta(k) ) < abs ( dbeta(k) - dbe(k) ) ) then
        done = .false.
        exit
      end if
    end do

    if ( done ) then
      exit
    end if

  end do

  return
end
subroutine mcdis_r4 ( n, ncapm, mc, mp, xp, yp, quad_r4, eps, iq, idelta, &
  irout, finl, finr, endl, endr, xfer, wfer, alpha, beta, ncap, kount, ierr, &
  ie )

!*****************************************************************************80
!
!! MCDIS_R4 is a multiple-component discretization procedure.
!
!  Discussion:
!
!    This routine is described in Section 4.3 of the companion paper.
!
!    It generates to a relative accuracy of EPS the recursion coefficients
!    alpha(k), beta(k), k = 0,1,...,n-1, for the polynomials orthogonal with
!    respect to a weight distribution consisting of the sum of  mc  continuous
!    components and a discrete component with  mp  points.
!
!    The continuous part of the spectrum is made up of  mc  weight functions,
!    each supported on its own interval. These intervals may or may not be
!    disjoint.
!
!    The discretization of the inner product on the i-th
!    interval is furnished either by a user-supplied subroutine QUAD_R4,
!    or by the general-purpose subroutine QGP_R4 provided in this package,
!    depending on whether  iq  is equal, or not equal, to  1, respectively.
!
!    The user-supplied routine must have the form  quad_R4(n,x,w,i,ierr)  and
!    is assumed to supply the abscissas  x(k)  and weights  w(k), k = 1,2,
!    ...,n, to be used in approximating the i-th inner product
!
!      integral of p(x)*q(x)*wf(x,i)dx
!
!    by the
!
!      sum over k from 1 to n of w(k)*p(x(k))*q(x(k)), i = 1,2,...,mc.
!
!    The desired recurrence coefficients are then approximated by the
!    recursion coefficients of the discrete orthogonal polynomials
!    belonging to the discretized inner product, which in turn are
!    computed by either the Stieltjes procedure or the Lanczos algorithm
!    according as  irout  is equal to, or not equal to  1, respectively.
!
!    Two error flags  ierr,ie  are provided which signal the occurrence
!    of an error condition in the quadrature process, or in the routine
!    STI_R4 or LANCZ_R4 (whichever is used), respectively. The point spectrum
!    is given through its abscissas  xp  and jumps  yp.
!
!    If the quadrature routine  quad  has polynomial degree of exactness
!    at least  id(n)  for each i, and if  id(n)/n  =  idelta + O(1/n)  as
!    n  goes to infinity, then the procedure is designed to converge after
!    one iteration, provided  idelta  is set with the appropriate
!    integer. Normally,  idelta = 1 (for interpolatory rules) or  idelta=2
!    (for Gaussian rules). The default value is  idelta = 1.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input:  n    - - the number of recursion coefficients desired;
!                     type integer
!
!            ncapm  - a discretization parameter indicating an upper
!                     limit of the fineness of the discretization;
!                     ncapm = 500  will usually be satisfactory; type
!                     integer
!
!            mc  - -  the number of disjoint intervals in the
!                     continuous part of the spectrum; type integer
!
!            mp  - -  the number of points in the discrete part of
!                     the spectrum; type integer. If there is no
!                     point spectrum, set  mp = 0.
!
!            xp  - -  an array of dimension  mp  containing the
!                     abscissas of the point spectrum
!
!            yp  - -  an array of dimension  mp  containing the jumps
!                     of the point spectrum
!
!            quad_r4  -  a subroutine determining the discretization of
!                     the inner product on each component interval,
!                     or a dummy routine if  iq  is not equal to  1
!                     (see below)
!
!            eps  - - the desired relative accuracy of the nonzero
!                     recursion coefficients; type real
!
!            iq   - - an integer selecting a user-supplied quadrature
!                     routine  quad  if  iq = 1  or the ORTHPOL routine
!                     qgp  otherwise
!
!            idelta - a nonzero integer, typically  1  or  2, inducing
!                     fast convergence in the case of special quadrature
!                     routines
!
!            irout  - an integer selecting the routine for generating
!                     the recursion coefficients from the discrete
!                     inner product. Specifically,  irout = 1  selects the
!                     routine STI_R4, whereas any other value selects the
!                     routine LANCZ_R4.
!
! The logical variables  finl,finr, the arrays endl,endr  of
! dimension  mc, and the arrays  xfer,wfer  of dimension  ncapm  are
! input variables to the subroutine  qgp  and are used (and hence need
! to be properly dimensioned) only if  iq  is not equal to  1.
!
!    Output:  alpha,beta - arrays of dimension n, holding as k-th
!                     element  alpha(k-1), beta(k-1), k = 1,2,...,n,
!                     respectively
!             ncap  - an integer indicating the fineness of the
!                     discretization that yields convergence within
!                     the eps-tolerance
!             kount - the number of iterations used
!             ierr  - an error flag, equal to  0  on normal return,
!                     equal to  -1  if  n  is not in the proper range,
!                     equal to  i  if there is an error condition in
!                     the discretization of the i-th interval,
!                     and equal to  ncapm  if the discretized
!                     Stieltjes procedure does not converge within the
!                     discretization resolution specified by  ncapm
!             ie - -  an error flag inherited from the routine STI_R4
!                     or LANCZ_R4 (whichever is used)
!
  implicit none

  integer mc
  integer mp
  integer n
  integer ncapm

  real alpha(n)
  real be(n)
  real beta(n)
  logical done
  real endl(mc)
  real endr(mc)
  real eps
  logical finl
  logical finr
  integer i
  integer idelta
  integer ie
  integer ierr
  integer im1tn
  integer incap
  integer iq
  integer irout
  integer k
  integer kount
  integer mtncap
  integer ncap
  real p0(mc*ncapm+mp)
  real p1(mc*ncapm+mp)
  real p2(mc*ncapm+mp)
  external quad_r4
  real w(ncapm)
  real wfer(ncapm)
  real wm(mc*ncapm+mp)
  real x(ncapm)
  real xfer(ncapm)
  real xm(mc*ncapm+mp)
  real xp(mp)
  real yp(mp)

  if ( idelta <= 0 ) then
    idelta = 1
  end if

  if ( n < 1 ) then
    ierr = - 1
    return
  end if
!
!  Initialization
!
  incap = 1
  kount = -1
  ierr = 0
  beta(1:n) = 0.0E+00
  ncap = ( 2 * n - 1 ) / idelta

  do

    be(1:n) = beta(1:n)
    kount = kount + 1

    if ( 1 < kount ) then
      incap = 2**( kount / 5 ) * n
    end if

    ncap = ncap + incap

    if ( ncapm < ncap ) then
      ierr = ncapm
      return
    end if
!
!  Discretization of the inner product
!
    mtncap = mc * ncap

    do i = 1, mc

      im1tn = ( i - 1 ) * ncap

      if ( iq == 1 ) then
        call quad_r4 ( ncap, x, w, i, ierr )
      else
        call qgp_r4 ( ncap, x, w, i, ierr, mc, finl, finr, endl, endr, &
          xfer, wfer )
      end if

      if ( ierr /= 0 ) then
        ierr = i
        return
      end if

      xm(im1tn+1:im1tn+ncap) = x(1:ncap)
      wm(im1tn+1:im1tn+ncap) = w(1:ncap)

    end do

    xm(mtncap+1:mtncap+mp) = xp(1:mp)
    wm(mtncap+1:mtncap+mp) = yp(1:mp)
!
!  Computation of the desired recursion coefficients
!
    if ( irout == 1 ) then
      call sti_r4 ( n, mtncap + mp, xm, wm, alpha, beta, ie, p0, p1, p2 )
    else
      call lancz_r4 ( n, mtncap + mp, xm, wm, alpha, beta, ie, p0, p1 )
    end if
!
!  In the following statement, the absolute value of the beta's is
!  used to guard against failure in cases where the routine is applied
!  to variable-sign weight functions and hence the positivity of the
!  beta's is not guaranteed.
!
    done = .true.

    do k = 1, n
      if ( eps * abs ( beta(k) ) < abs ( beta(k) - be(k) ) ) then
        done = .false.
        exit
      end if
    end do

    if ( done ) then
      exit
    end if

  end do

  return
end
subroutine mcdis_r8 ( n, ncapm, mc, mp, dxp, dyp, quad_r8, deps, iq, idelta, &
  irout, finld, finrd, dendl, dendr, dxfer, dwfer, dalpha, dbeta, ncap, &
  kount, ierrd, ied )

!*****************************************************************************80
!
!! MCDIS_R8 is a double-precision version of the routine MCDIS_R4.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  integer mc
  integer mp
  integer n
  integer ncapm

  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) dbe(n)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dendl(mc)
  real ( kind = 8 ) dendr(mc)
  real ( kind = 8 ) deps
  logical done
  real ( kind = 8 ) dp0(mc*ncapm+mp)
  real ( kind = 8 ) dp1(mc*ncapm+mp)
  real ( kind = 8 ) dp2(mc*ncapm+mp)
  real ( kind = 8 ) dw(ncapm)
  real ( kind = 8 ) dwfer(ncapm)
  real ( kind = 8 ) dwm(mc*ncapm+mp)
  real ( kind = 8 ) dx(ncapm)
  real ( kind = 8 ) dxfer(ncapm)
  real ( kind = 8 ) dxm(mc*ncapm+mp)
  real ( kind = 8 ) dxp(mp)
  real ( kind = 8 ) dyp(mp)
  logical finld
  logical finrd
  integer i
  integer idelta
  integer ied
  integer ierr
  integer ierrd
  integer im1tn
  integer incap
  integer iq
  integer irout
  integer k
  integer kount
  integer mtncap
  integer ncap
  external quad_r8

  if ( idelta <= 0 ) then
    idelta = 1
  end if

  if ( n < 1 ) then
    ierrd = - 1
    return
  end if

  incap = 1
  kount = -1
  ierr = 0
  dbeta(1:n) = 0.0D+00

  ncap = ( 2 * n - 1 ) / idelta

  do

    dbe(1:n) = dbeta(1:n)

    kount = kount + 1

    if ( 1 < kount ) then
      incap = 2**( kount / 5 ) * n
    end if

    ncap = ncap + incap

    if ( ncapm < ncap ) then
      ierrd = ncapm
      return
    end if

    mtncap = mc * ncap

    do i = 1, mc

      im1tn = ( i - 1 ) * ncap

      if ( iq == 1 ) then
        call quad_r8 ( ncap, dx, dw, i, ierr )
      else
        call qgp_r8 ( ncap, dx, dw, i, ierr, mc, finld, finrd, dendl, &
          dendr, dxfer, dwfer )
      end if

      if ( ierr /= 0 ) then
        ierrd = i
        return
      end if

      dxm(im1tn+1:im1tn+ncap) = dx(1:ncap)
      dwm(im1tn+1:im1tn+ncap) = dw(1:ncap)

    end do

    dxm(mtncap+1:mtncap+mp) = dxp(1:mp)
    dwm(mtncap+1:mtncap+mp) = dyp(1:mp)

    if ( irout == 1 ) then
      call sti_r8 ( n, mtncap+mp, dxm, dwm, dalpha, dbeta, ied, dp0, dp1, dp2 )
    else
      call lancz_r8 ( n, mtncap+mp, dxm, dwm, dalpha, dbeta, ied, dp0, dp1 )
    end if

    done = .true.

    do k = 1, n
      if ( deps * abs ( dbeta(k) ) < abs ( dbeta(k) - dbe(k) ) ) then
        done = .false.
        exit
      end if
    end do

    if ( done ) then
      exit
    end if

  end do

  return
end
function nu0her ( n, z, eps )

!*****************************************************************************80
!
!! NU0HER estimates a starting index for recursion with the Hermite measure.
!
!  Discussion:
!
!    This routine provides a starting backward recurrence index for the Hermite
!    measure that can be used in place of NU0 in the routines KNUM_R4 and
!    KNUM_R8.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  real eps
  integer n
  integer nu0her
  complex z

  nu0her = int ( 2.0E+00 * ( sqrt ( 0.5E+00 * real ( n + 1, kind = 4 ) ) &
    + 0.25E+00 * log ( 1.0E+00 / eps ) / abs ( aimag ( z ) ) )**2 )

  return
end
function nu0jac ( n, z, eps )

!*****************************************************************************80
!
!! NU0JAC estimates a starting index for recursion with the Jacobi measure.
!
!  Discussion:
!
!    This routine provides a starting backward recurrence index for the Jacobi
!    measure that can be used in place of NU0 in the routines KNUM_R4
!    and KNUM_R8.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  real angle
  real eps
  integer n
  integer nu0jac
  real pi
  real r
  real x
  real x2
  real y
  real y2
  complex z

  pi = 4.0E+00 * atan ( 1.0E+00 )
  x = real ( z, kind = 4 )
  y = abs ( aimag ( z ) )

  if ( x < -1.0E+00 ) then
    angle = 0.5E+00 * ( 2.0E+00 * pi + atan ( y / ( x - 1.0E+00 ) ) &
      + atan ( y / ( x + 1.0E+00 ) ) )
  else if ( x == -1.0E+00 ) then
    angle = 0.5E+00 * ( 1.5E+00 * pi - atan ( 0.5E+00 * y ) )
  else if ( x < 1.0E+00 ) then
    angle = 0.5E+00 * ( pi + atan ( y / ( x - 1.0E+00 ) ) &
      + atan ( y / ( x + 1.0E+00 ) ) )
  else if ( x == 1.0E+00 ) then
    angle = 0.5E+00 * ( 0.5E+00 * pi + atan ( 0.5E+00 * y ) )
  else if ( 1.0E+00 < x ) then
    angle = 0.5E+00 * ( atan ( y / ( x - 1.0E+00 ) ) &
      + atan ( y / ( x + 1.0E+00 ) ) )
  end if

  x2 = x * x
  y2 = y * y
  r = ( ( x2 - y2 - 1.0E+00 )**2 + 4.0E+00 * x2 * y2 )**0.25E+00
  r = sqrt ( ( x + r * cos ( angle ) )**2 + ( y + r * sin ( angle ) )**2)

  nu0jac = int ( real ( n + 1, kind = 4 ) + 0.5E+00 * log ( 1.0E+00 / eps ) &
    / log ( r ) )

  return
end
function nu0lag ( n, z, al, eps )

!*****************************************************************************80
!
!! NU0LAG estimates a starting index for recursion with the Laguerre measure.
!
!  Discussion:
!
!    This routine provides a starting backward recurrence index for the
!    Laguerre measure that can be used in place of NU0 in the routines KNUM_R4
!    and KNUM_R8.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  real al
  real eps
  integer n
  integer nu0lag
  real phi
  real pi
  real x
  real y
  complex z

  pi = 4.0E+00 * atan ( 1.0E+00 )
  x = real ( z, kind = 4 )
  y = aimag ( z )

  if ( y < 0.0E+00 ) then
    phi = 1.5E+00 * pi
  else
    phi = 0.5E+00 * pi
  end if

  if ( x /= 0.0E+00 ) then

    phi = atan ( y / x )

    if ( y <= 0.0E+00 .or. x <= 0.0E+00 ) then

      phi = phi + pi

      if ( 0.0E+00 <= x ) then
        phi = phi + pi
      end if

    end if

  end if

  nu0lag = int ( &
    ( sqrt ( real ( n + 1, kind = 4 ) + 0.5E+00 * ( al + 1.0E+00 ) ) &
    + log ( 1.0E+00 / eps ) &
    / ( 4.0E+00 * ( x * x + y * y )**0.25E+00 &
    * cos ( 0.5E+00 * ( phi - pi ) ) ) &
    )**2 - 0.5E+00 * ( al + 1.0E+00 ) )

  return
end
subroutine qgp_r4 ( n, x, w, i, ierr, mc, finl, finr, endl, endr, xfer, wfer )

!*****************************************************************************80
!
!! QGP_R4 is a general-purpose discretization routine.
!
!  Discussion:
!
!    This routine can be used as an alternative to the routine  quad  in the
!    multiple-component discretization procedure  mcdis.  It takes no account
!    of the special nature of the weight function involved and hence may result
!    in slow convergence of the discretization procedure.  This routine,
!    therefore, should be used only as a last resort, when no better, more
!    natural discretization can be found.
!
!    It is assumed that there are 1 <= MC disjoint component intervals.
!    The discretization is effected by the Fejer quadrature rule,
!    suitably transformed to the respective interval. An interval that
!    extends to minus infinity has to be indexed by  1; one that extends
!    to plus infinity has to be indexed by  mc.
!
!    The user has to supply the routine
!      function wf_r4 ( x, i ),
!    which evaluates the weight function at the point  x  on the i-th
!    component interval.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    The output variable  ierr  is given the value  0. Additional input
!    parameters and working space used by this routine are as follows:
!
!          mc      - the number of component intervals; type integer
!
!          finl    - a logical variable to be set .true. if the
!                    extreme left interval is finite, and .false.
!                    otherwise
!          finr    - a logical variable to be set .true. if the
!                    extreme right interval is finite, and .false.
!                    otherwise
!
!        endl    - an array of dimension  mc  containing the left
!                  endpoints of the component intervals; if the
!                    first of these extends to minus infinity,endl(1)
!                    can be set to an arbitrary value
!
!        endr    - an array of dimension  mc  containing the right
!                  endpoints of the component intervals; if the
!                    last of these extends to plus infinity,endr(mc)
!                    can be set to an arbitrary value
!
!          xfer,wfer-working arrays holding the Fejer nodes and
!                    weights, respectively, for the interval [-1,1].
!
  implicit none

  integer mc
  integer n

  real endl(mc)
  real endr(mc)
  logical finl
  logical finr
  integer i
  integer ierr
  integer k
  real phi
  real phi1
  real w(n)
  real wf_r4
  real wfer(*)
  real x(n)
  real xfer(*)
!
!  The arrays  xfer,wfer  are dimensioned in the routine  mcdis.
!
  ierr = 0

  if ( i == 1 ) then
    call fejer_r4 ( n, xfer, wfer )
  end if

  if ( 1 < i .and. i < mc ) go to 60

  if ( mc == 1 ) then

    if ( finl .and. finr ) go to 60
    if ( finl ) go to 20
    if ( finr ) go to 40

    do k = 1, n
      call symtr_r4 ( xfer(k), phi, phi1 )
      x(k) = phi
      w(k) = wfer(k) * wf_r4 ( phi, i ) * phi1
    end do

    return

  else

    if ( ( i == 1 .and. finl ) .or. ( i .eq. mc .and. finr ) ) then
      go to 60
    end if

    if ( i == 1 ) go to 40

  end if

   20 continue

  do k = 1, n
    call tr_r4 ( xfer(k), phi, phi1 )
    x(k) = endl(mc) + phi
    w(k) = wfer(k) * wf_r4 ( x(k), mc ) * phi1
  end do

  return

   40 continue

  do k = 1, n
    call tr_r4 ( -xfer(k), phi, phi1 )
    x(k) = endr(1) - phi
    w(k) = wfer(k) * wf_r4 ( x(k), 1 ) * phi1
  end do

  return

   60 continue

  do k = 1, n
    x(k) = 0.5E+00 * ( ( endr(i) - endl(i) ) * xfer(k) + endr(i) + endl(i) )
    w(k) = 0.5E+00 * ( endr(i) - endl(i) ) * wfer(k) * wf_r4 ( x(k), i )
  end do

  return
end
subroutine qgp_r8 ( n, dx, dw, i, ierr, mcd, finld, finrd, dendl, dendr, &
  dxfer, dwfer )

!*****************************************************************************80
!
!! QGP_R8 is a double-precision version of the routine QGP_R4.
!
!  Discussion:
!
!    The user has to supply the routine
!
!       function wf_r8 ( dx, i )
!
!    which evaluates the weight function in real ( kind = 8 ) at the
!    point  dx  on the i-th component interval.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  integer mcd
  integer n

  real ( kind = 8 ) dendl(mcd)
  real ( kind = 8 ) dendr(mcd)
  real ( kind = 8 ) dphi
  real ( kind = 8 ) dphi1
  real ( kind = 8 ) dw(n)
  real ( kind = 8 ) dwfer(*)
  real ( kind = 8 ) dx(n)
  real ( kind = 8 ) dxfer(*)
  logical finld
  logical finrd
  integer i
  integer ierr
  integer k
  real ( kind = 8 ) wf_r8
!
!  The arrays  dxfer,dwfer  are dimensioned in the routine  dmcdis.
!
  ierr = 0

  if ( i == 1 ) then
    call fejer_r8 ( n, dxfer, dwfer )
  end if

  if ( 1 < i .and. i < mcd ) then
    go to 60
  end if

  if ( mcd == 1 ) then

    if ( finld .and. finrd ) go to 60
    if ( finld ) go to 20
    if ( finrd ) go to 40

    do k = 1, n
      call symtr_r8 ( dxfer(k), dphi, dphi1 )
      dx(k) = dphi
      dw(k) = dwfer(k) * wf_r8 ( dphi, i ) * dphi1
    end do

    return

  else

    if ( (i == 1 .and. finld ) .or. ( i == mcd .and. finrd ) ) go to 60
    if ( i == 1 ) go to 40

  end if

   20 continue

  do k = 1, n
    call tr_r8 ( dxfer(k), dphi, dphi1 )
    dx(k) = dendl(mcd) + dphi
    dw(k) = dwfer(k) * wf_r8 ( dx(k), mcd ) * dphi1
  end do

  return

   40 continue

  do k = 1, n
    call tr_r8 ( -dxfer(k), dphi, dphi1 )
    dx(k) = dendr(1) - dphi
    dw(k) = dwfer(k) * wf_r8 ( dx(k), 1 ) * dphi1
  end do

  return

   60 continue

  do k = 1, n
    dx(k) = 0.5D+00 * ( ( dendr(i) - dendl(i) ) * dxfer(k) + dendr(i) + dendl(i) )
    dw(k) = 0.5D+00 * ( dendr(i) - dendl(i) ) * dwfer ( k ) * wf_r8 ( dx(k), i )
  end do

  return
end
subroutine radau_r4 ( n, alpha, beta, endl, zero, weight, ierr )

!*****************************************************************************80
!
!! RADAU_R4 generates a Gauss-Radau quadrature formula.
!
!  Discussion:
!
!    Given N and a measure  dlambda, this routine generates the
!    (N+1)-point Gauss-Radau quadrature formula
!
!      integral over supp(dlambda) of f(t)dlambda(t)
!
!      =  w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
!
!    The nodes are returned as  zero(k) = x(k), the weights as  weight(k)
!    = w(k), k=0,1,2,...,n. The user has to supply the recursion
!    coefficients  alpha(k), beta(k), k = 0,1,2,...,n, for the measure
!    dlambda. The nodes and weights are computed as eigenvalues and
!    in terms of the first component of the respective normalized
!    eigenvectors of a slightly modified Jacobi matrix of order  n+1.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of interior points in the Gauss-Radau
!    formula.
!
!    Input, real ALPHA(N+1), BETA(N+1), the recursion coefficients
!    alpha(k-1), beta(k-1), k = 1,2,...,n+1; the coefficient  alpha(n+1)  is not
!    used by the routine.
!
!    Input, real ENDL, the prescribed endpoint x(0) of the Gauss-Radau
!    formula.
!
!    Output, real ZERO(N+1), the nodes (in increasing order)  zero(k) = x(k),
!    k=0,1,2,...,n
!
!    Output, real WEIGHT(N+1), the weights weight(k) = w(k), k=0,1,2,...,n
!
!    Output, integer IERR, an error flag from GAUSS_R4.
!
  implicit none

  integer n

  real a(n+1)
  real alpha(n+1)
  real b(n+1)
  real beta(n+1)
  real e(n+1)
  real endl
  real epsma
  integer ierr
  integer k
  real p0
  real p1
  real pm1
  real weight(n+1)
  real zero(n+1)

  a(1:n) = alpha(1:n)
  b(1:n+1) = beta(1:n+1)
!
!  Determine A(N+1).
!
  p0 = 0.0E+00
  p1 = 1.0E+00

  do k = 1, n
    pm1 = p0
    p0 = p1
    p1 = ( endl - a(k) ) * p0 - b(k) * pm1
  end do

  a(n+1) = endl - b(n+1) * p0 / p1
!
!  Call GAUSS_R4.
!
  epsma = epsilon ( epsma )

  call gauss_r4 ( n+1, a, b, epsma, zero, weight, ierr, e )

  return
end
subroutine radau_r8 ( n, alpha, beta, endl, zero, weigh, ierr )

!*****************************************************************************80
!
!! RADAU_R8 generates a Gauss-Radau quadrature formula.
!
!  Discussion:
!
!    Given N and a measure  dlambda, this routine generates the
!    (N+1)-point Gauss-Radau quadrature formula
!
!      integral over supp(dlambda) of f(t)dlambda(t)
!
!      =  w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
!
!    The nodes are returned as  zero(k) = x(k), the weights as  weight(k)
!    = w(k), k=0,1,2,...,n. The user has to supply the recursion
!    coefficients  alpha(k), beta(k), k = 0,1,2,...,n, for the measure
!    dlambda. The nodes and weights are computed as eigenvalues and
!    in terms of the first component of the respective normalized
!    eigenvectors of a slightly modified Jacobi matrix of order  n+1.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of interior points in the Gauss-Radau
!    formula.
!
!    Input, real ( kind = 8 ) ALPHA(N+1), BETA(N+1), the recursion coefficients
!    alpha(k-1), beta(k-1), k = 1,2,...,n+1; the coefficient  alpha(n+1)  is not
!    used by the routine.
!
!    Input, real ( kind = 8 ) ENDL, the prescribed endpoint x(0) of the Gauss-Radau
!    formula.
!
!    Output, real ( kind = 8 ) ZERO(N+1), the nodes (in increasing order)
!    zero(k) = x(k), k=0,1,2,...,n
!
!    Output, real ( kind = 8 ) WEIGHT(N+1), the weights weight(k) = w(k),
!    k=0,1,2,...,n
!
!    Output, integer IERR, an error flag from GAUSS_R8.
!
  implicit none

  integer n

  real ( kind = 8 ) a(n+1)
  real ( kind = 8 ) alpha(n+1)
  real ( kind = 8 ) b(n+1)
  real ( kind = 8 ) beta(n+1)
  real ( kind = 8 ) e(n+1)
  real ( kind = 8 ) endl
  real ( kind = 8 ) epsma
  integer ierr
  integer k
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) pm1
  real ( kind = 8 ) weigh(n+1)
  real ( kind = 8 ) zero(n+1)

  epsma = epsilon ( epsma )

  a(1:n) = alpha(1:n)
  b(1:n+1) = beta(1:n+1)
!
!  Determine A(N+1).
!
  p0 = 0.0D+00
  p1 = 1.0D+00
  do k = 1, n
    pm1 = p0
    p0 = p1
    p1 = ( endl - a(k) ) * p0 - b(k) * pm1
  end do

  a(n+1) = endl - b(n+1) * p0 / p1
!
!  Call GAUSS_R8.
!
  call gauss_r8 ( n+1, a, b, epsma, zero, weigh, ierr, e )

  return
end
subroutine recur_r4 ( n, ipoly, al, be, a, b, ierr )

!*****************************************************************************80
!
!! RECUR_R4 generates recursion coefficients for orthogonal polynomials.
!
!  Discussion:
!
!    The routine computes the coefficients  a(k),b(k), k = 0,1,...,n-1,
!    in the recurrence relation
!
!      p(k+1)(x) = (x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
!        k = 0,1,...,n-1,
!
!      p(-1)(x) = 0,  p(0)(x)=1,
!
!    for some classical (monic) orthogonal polynomials, and sets  b(0)
!    equal to the total mass of the weight distribution. The results are
!    stored in the arrays  a,b,  which hold, respectively, the coefficients
!    a(k-1),b(k-1), k = 1,2,...,n.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer N, the number of recursion coefficients desired.
!
!    Input, integer IPOLY, chooses the polynomial as follows:
!    1, Legendre polynomial on (-1,1)
!    2, Legendre polynomial on (0,1)
!    3, Chebyshev polynomial of the first kind
!    4, Chebyshev polynomial of the second kind
!    5, Jacobi polynomial with parameters  al=-.5,be=.5
!    6, Jacobi polynomial with parameters  al,be
!    7, generalized Laguerre polynomial with parameter  al
!    8, Hermite polynomial
!
!    Input, real AL, BE, parameters for Jacobi and generalized
!    Laguerre polynomials.
!
!    Output, real A(N), B(N), the recursion coefficients
!    a(k-1),b(k-1), k = 1,2,...,n.
!
!    Output, integer IERR, an error flag, equal to 0 on normal return,
!    equal to 1 if al or be are out of range
!    when ipoly = 6 or ipoly = 7, equal to 2 if b(0)
!    overflows when ipoly = 6 or ipoly = 7, equal to 3
!    if n is out of range, and equal to 4 if ipoly
!    is not an admissible integer.  In the case ierr = 2,
!    the coefficient b(0) is set equal to the largest
!    machine-representable number.
!
  implicit none

  integer n

  real a(n)
  real al
  real al2
  real alga_r4
  real almach
  real alpbe
  real b(n)
  real be
  real be2
  real fkm1
  real, external :: gamma_r4
  integer ierr
  integer ipoly
  integer k
  real t

  if ( n < 1 ) then
    ierr = 3
    return
  end if

  almach = log ( huge ( almach ) )
  ierr = 0
  a(1:n) = 0.0E+00

  if ( ipoly == 1 ) then

    b(1) = 2.0E+00
    do k = 2, n
      fkm1 = real ( k - 1, kind = 4 )
      b(k) = 1.0E+00 / ( 4.0E+00 - 1.0E+00 / ( fkm1 * fkm1 ) )
    end do

  else if ( ipoly == 2 ) then

    a(1) = 0.5E+00
    b(1) = 1.0E+00
    do k = 2, n
      a(k) = 0.5E+00
      fkm1 = real ( k - 1, kind = 4 )
      b(k) = 0.25E+00 / ( 4.0E+00 - 1.0E+00 / ( fkm1 * fkm1 ) )
    end do

  else if ( ipoly == 3 ) then

    b(1) = 4.0E+00 * atan ( 1.0E+00 )

    if ( n == 1 ) then
      return
    end if

    b(2) = 0.5E+00
    b(3:n) = 0.25E+00

  else if ( ipoly == 4 ) then

    b(1) = 2.0E+00 * atan ( 1.0E+00 )
    b(2:n) = 0.25E+00

  else if ( ipoly == 5 ) then

    b(1) = 4.0E+00 * atan ( 1.0E+00 )
    a(1) = 0.5E+00
    b(2:n) = 0.25E+00

  else if ( ipoly == 6 ) then

    if ( al <= -1.0E+00 .or. be <= -1.0E+00 ) then
      ierr = 1
      return
    end if

    alpbe = al + be
    a(1) = ( be - al ) / ( alpbe + 2.0E+00 )
    t = ( alpbe + 1.0E+00 ) * log ( 2.0E+00 ) + alga_r4 ( al + 1.0E+00 ) &
      + alga_r4 ( be + 1.0E+00 ) - alga_r4 ( alpbe + 2.0E+00 )

    if ( almach < t ) then
      ierr = 2
      b(1) = huge ( b(1) )
    else
      b(1) = exp ( t )
    end if

    if ( n == 1 ) then
      return
    end if

    al2 = al * al
    be2 = be * be
    a(2) = ( be2 - al2 ) / ( ( alpbe + 2.0E+00 ) * ( alpbe + 4.0E+00 ) )
    b(2) = 4.0E+00 * ( al + 1.0E+00 ) * ( be + 1.0E+00 ) &
      / ( ( alpbe + 3.0E+00 ) * ( alpbe + 2.0E+00 )**2 )

    do k = 3, n

      fkm1 = real ( k - 1, kind = 4 )

      a(k) = 0.25E+00 * ( be2 - al2 ) &
        / ( fkm1 * fkm1 * ( 1.0E+00 + 0.5E+00 * alpbe / fkm1 ) &
        * ( 1.0E+00 + 0.5E+00 * ( alpbe + 2.0E+00 ) / fkm1 ) )

      b(k) = 0.25E+00 * ( 1.0E+00 + al / fkm1 ) * ( 1.0E+00 + be / fkm1 ) &
        * ( 1.0E+00 + alpbe / fkm1 ) &
        / ( ( 1.0E+00 + 0.5E+00 * ( alpbe + 1.0E+00 ) / fkm1 ) &
        * ( 1.0E+00 + 0.5E+00 * ( alpbe - 1.0E+00 ) / fkm1 ) &
        * ( 1.0E+00 + 0.5E+00 * alpbe / fkm1 )**2 )
    end do

  else if ( ipoly == 7 ) then

    if ( al <= -1.0E+00 ) then
      ierr = 1
      return
    end if

    a(1) = al + 1.0E+00
    b(1) = gamma_r4 ( al + 1.0E+00, ierr )

    if ( ierr == 2 ) then
      b(1) = huge ( b(1) )
    end if

    do k = 2, n
      fkm1 = real ( k - 1, kind = 4 )
      a(k) = 2.0E+00 * fkm1 + al + 1.0E+00
      b(k) = fkm1 * ( fkm1 + al )
    end do

  else if ( ipoly == 8 ) then

    b(1) = sqrt ( 4.0E+00 * atan ( 1.0E+00 ) )
    do k = 2, n
      b(k) = 0.5E+00 * real ( k - 1, kind = 4 )
    end do

  else

    ierr = 4

  end if

  return
end
subroutine recur_r8 ( n, ipoly, dal, dbe, da, db, ierr )

!*****************************************************************************80
!
!! RECUR_R8 is a double-precision version of the routine RECUR_R4.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  integer n

  real ( kind = 8 ) alga_r8
  real ( kind = 8 ) da(n)
  real ( kind = 8 ) dal
  real ( kind = 8 ) dal2
  real ( kind = 8 ) dalpbe
  real ( kind = 8 ) db(n)
  real ( kind = 8 ) dbe
  real ( kind = 8 ) dbe2
  real ( kind = 8 ) dkm1
  real ( kind = 8 ) dlmach
  real ( kind = 8 ) dt
  real ( kind = 8 ), external :: gamma_r8
  integer ierr
  integer ipoly
  integer k

  if ( n < 1 ) then
    ierr = 3
    return
  end if

  dlmach = log ( huge ( dlmach ) )
  ierr = 0
  da(1:n) = 0.0D+00

  if ( ipoly == 1 ) then

    db(1) = 2.0D+00
    do k = 2, n
      dkm1 = real ( k - 1, kind = 8 )
      db(k) = 1.0D+00 / ( 4.0D+00 - 1.0D+00 / ( dkm1 * dkm1 ) )
    end do

  else if ( ipoly == 2 ) then

    da(1) = 0.5D+00
    db(1) = 1.0D+00
    do k = 2, n
      da(k) = 0.5D+00
      dkm1 = real ( k - 1, kind = 8 )
      db(k) = 0.25D+00 / ( 4.0D+00 - 1.0D+00 / ( dkm1 * dkm1 ) )
    end do

  else if ( ipoly == 3 ) then

    db(1) = 4.0D+00 * atan ( 1.0D+00 )
    if ( n == 1 ) then
      return
    end if
    db(2) = 0.5D+00
    db(3:n) = 0.25D+00

  else if ( ipoly == 4 ) then

    db(1) = 2.0D+00 * atan ( 1.0D+00 )
    db(2:n) = 0.25D+00

  else if ( ipoly == 5 ) then

    db(1) = 4.0D+00 * atan ( 1.0D+00 )
    da(1) = 0.5D+00
    db(2:n) = 0.25D+00

  else if ( ipoly == 6 ) then

    if ( dal <= -1.0D+00 .or. dbe <= -1.0D+00 ) then
      ierr = 1
      return
    end if

    dalpbe = dal + dbe
    da(1) = ( dbe - dal ) / ( dalpbe + 2.0D+00 )
    dt = ( dalpbe + 1.0D+00 ) * log ( 2.0D+00 ) + alga_r8 ( dal + 1.0D+00 ) &
      + alga_r8 ( dbe + 1.0D+00 ) - alga_r8 ( dalpbe + 2.0D+00 )

    if ( dlmach < dt ) then
      ierr = 2
      db(1) = huge ( db(1) )
    else
      db(1) = exp ( dt )
    end if

    if ( n == 1 ) then
      return
    end if

    dal2 = dal * dal
    dbe2 = dbe * dbe
    da(2) = ( dbe2 - dal2 ) / ( ( dalpbe + 2.0D+00 ) * ( dalpbe + 4.0D+00 ) )
    db(2) = 4.0D+00 * ( dal + 1.0D+00 ) * ( dbe + 1.0D+00 ) &
      / ( ( dalpbe + 3.0D+00 ) * ( dalpbe + 2.0D+00 )**2 )

    do k = 3, n

      dkm1 = real ( k - 1, kind = 8 )

      da(k) = 0.25D+00 * ( dbe2 - dal2 ) &
        / ( dkm1 * dkm1 * ( 1.0D+00 + 0.5D+00 * dalpbe / dkm1 ) &
        * ( 1.0D+00 + 0.5D+00 * ( dalpbe + 2.0D+00 ) / dkm1 ) )

      db(k) = 0.25D+00 * ( 1.0D+00 + dal / dkm1 ) &
        * ( 1.0D+00 + dbe / dkm1 ) * ( 1.0D+00 + dalpbe / dkm1 ) &
        / ( ( 1.0D+00 + 0.5D+00 * ( dalpbe + 1.0D+00 ) / dkm1 ) &
        * ( 1.0D+00 + 0.5D+00 * ( dalpbe -1.0D+00 ) / dkm1 ) &
        * ( 1.0D+00 + 0.5D+00 * dalpbe / dkm1 )**2 )

    end do

  else if ( ipoly == 7 ) then

    if ( dal <= -1.0D+00 ) then
      ierr = 1
      return
    end if

    da(1) = dal + 1.0D+00
    db(1) = gamma_r8 ( dal + 1.0D+00, ierr )

    if ( ierr == 2 ) then
      db(1) = huge ( db(1) )
      return
    end if

    do k = 2, n
      dkm1 = real ( k - 1, kind = 8 )
      da(k) = 2.0D+00 * dkm1 + dal + 1.0D+00
      db(k) = dkm1 * ( dkm1 + dal )
    end do

  else if ( ipoly == 8 ) then

    db(1) = sqrt ( 4.0D+00 * atan ( 1.0D+00 ) )

    do k = 2, n
      db(k) = 0.5D+00 * real ( k - 1, kind = 8 )
    end do

  else

    ierr = 4

  end if

  return
end
subroutine sti_r4 ( n, ncap, x, w, alpha, beta, ierr, p0, p1, p2 )

!*****************************************************************************80
!
!! STI_R4 applies Stieltjes's procedure.
!
!  Discussion:
!
!    This routine generates the recursion
!    coefficients  alpha(k), beta(k) , k = 0,1,...,n-1, for the discrete
!    monic orthogonal polynomials associated with the inner product
!
!      (f,g) = sum over k from 1 to ncap of w(k)*f(x(k))*g(x(k)).
!
!    The integer  n  must be between  1  and  ncap, inclusive; otherwise,
!    there is an error exit with  ierr = 1. The results are stored in the
!    arrays  alpha, beta; the arrays  p0, p1, p2  are working arrays.
!
!    If there is a threat of underflow or overflow in the calculation
!    of the coefficients  alpha(k)  and  beta(k), the routine exits with
!    the error flag  ierr  set equal to  -k  (in the case of underflow)
!    or  +k  (in the case of overflow), where  k  is the recursion index
!    for which the problem occurs. The former [latter] can often be avoided
!    by multiplying all weights  w(k)  by a sufficiently large [small]
!    scaling factor prior to entering the routine, and, upon exit, divide
!    the coefficient  beta(0)  by the same factor.
!
!    This routine should be used with caution if  n  is relatively close
!    to ncap, since there is a distinct possibility of numerical
!    instability developing. (See W. Gautschi,Is the recurrence relation
!    for orthogonal polynomials always stable?'', BIT, 1993, to appear.)
!    In that case, the routine LANCZ_R4 should be used.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    On Generating Orthogonal Polynomials,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 3, Number 3, 1982, pages 289-317.
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  integer n
  integer ncap

  real alpha(n)
  real beta(n)
  integer ierr
  integer k
  integer m
  real p0(ncap)
  real p1(ncap)
  real p2(ncap)
  real sum0
  real sum1
  real sum2
  real t
  real w(ncap)
  real x(ncap)

  ierr = 0

  if ( n <= 0 .or. ncap < n ) then
    ierr = 1
    return
  end if
!
!  Compute the first alpha- and beta-coefficient.
!
  sum0 = sum ( w(1:ncap) )
  sum1 = dot_product ( w(1:ncap), x(1:ncap) )

  alpha(1) = sum1 / sum0
  beta(1) = sum0

  if ( n == 1 ) then
    return
  end if
!
!  Compute the remaining ALPHA's and BETA's.
!
  p1(1:ncap) = 0.0E+00
  p2(1:ncap) = 1.0E+00

  do k = 1, n - 1

    sum1 = 0.0E+00
    sum2 = 0.0E+00

    do m = 1, ncap
!
!  The following statement is designed to avoid an overflow condition
!  in the computation of p2(m) when the weights w(m) go to zero
!  faster (and underflow) than the p2(m) grow.
!
      if ( w(m) /= 0.0E+00 ) then

        p0(m) = p1(m)
        p1(m) = p2(m)
        p2(m) = ( x(m) - alpha(k) ) * p1(m) - beta(k) * p0(m)
!
!  Check for impending overflow.
!
        if ( 0.1E+00 * huge ( p2(m) ) < abs ( p2(m) ) .or. &
             0.1E+00 * huge ( sum2 ) < abs ( sum2 ) ) then
          ierr = k
          return
        end if

        t = w(m) * p2(m) * p2(m)
        sum1 = sum1 + t
        sum2 = sum2 + t * x(m)

      end if

    end do
!
!  Check for impending underflow.
!
    if ( abs ( sum1 ) < 10.0E+00 * tiny ( sum1 ) ) then
      ierr = - k
      return
    end if

    alpha(k+1) = sum2 / sum1
    beta(k+1) = sum1 / sum0
    sum0 = sum1

  end do

  return
end
subroutine sti_r8 ( n, ncap, dx, dw, dalpha, dbeta, ierr, dp0, dp1, dp2 )

!*****************************************************************************80
!
!! STI_R8 is a double-precision version of the routine STI_R4.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    On Generating Orthogonal Polynomials,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 3, Number 3, 1982, pages 289-317.
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  integer n
  integer ncap

  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dp0(ncap)
  real ( kind = 8 ) dp1(ncap)
  real ( kind = 8 ) dp2(ncap)
  real ( kind = 8 ) dsum0
  real ( kind = 8 ) dsum1
  real ( kind = 8 ) dsum2
  real ( kind = 8 ) dt
  real ( kind = 8 ) dw(ncap)
  real ( kind = 8 ) dx(ncap)
  integer ierr
  integer k
  integer m

  ierr = 0

  if ( n <= 0 .or. ncap < n ) then
    ierr = 1
    return
  end if

  dsum0 = sum ( dw(1:ncap) )
  dsum1 = dot_product ( dw(1:ncap), dx(1:ncap) )

  dalpha(1) = dsum1 / dsum0
  dbeta(1) = dsum0

  if ( n == 1 ) then
    return
  end if

  dp1(1:ncap) = 0.0D+00
  dp2(1:ncap) = 1.0D+00

  do k = 1, n - 1

    dsum1 = 0.0D+00
    dsum2 = 0.0D+00

    do m = 1, ncap

      if ( dw(m) /= 0.0D+00 ) then

        dp0(m) = dp1(m)
        dp1(m) = dp2(m)
        dp2(m) = ( dx(m) - dalpha(k) ) * dp1(m) - dbeta(k) * dp0(m)

        if ( 0.1D+00 * huge ( dp2(m) ) < abs ( dp2(m) ) .or. &
          0.1D+00 * huge ( dsum2 ) < abs ( dsum2 ) ) then
          ierr = k
          return
        end if

        dt = dw(m) * dp2(m) * dp2(m)
        dsum1 = dsum1 + dt
        dsum2 = dsum2 + dt * dx(m)

      end if

    end do

    if ( abs ( dsum1 ) < 10.0D+00 * tiny ( dsum1 ) ) then
      ierr = - k
      return
    end if

    dalpha(k+1) = dsum2 / dsum1
    dbeta(k+1) = dsum1 / dsum0
    dsum0 = dsum1

  end do

  return
end
subroutine symtr_r4 ( t, phi, phi1 )

!*****************************************************************************80
!
!! SYMTR_R4 maps T in [-1,1] to X in (-oo,oo).
!
!  Discussion:
!
!    X = T / ( 1 - T * T )
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, real T, the point in [-1,1] to be mapped.
!    -1 < T < 1.
!
!    Output, real PHI, the value of X(T).
!
!    Output, real PHI1, the derivative of the mapping.
!
  implicit none

  real phi
  real phi1
  real t
  real t2

  t2 = t * t
  phi = t / ( 1.0E+00 - t * t )
  phi1 = ( t * t + 1.0E+00 ) / ( t * t - 1.0E+00 )**2

  return
end
subroutine symtr_r8 ( t, phi, phi1 )

!*****************************************************************************80
!
!! SYMTR_R8 maps T in [-1,1] to X in (-oo,oo).
!
!  Discussion:
!
!    X = T / ( 1 - T * T )
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, real T, the point in [-1,1] to be mapped.
!    -1 < T < 1.
!
!    Output, real PHI, the value of X(T).
!
!    Output, real PHI1, the derivative of the mapping.
!
  implicit none

  real ( kind = 8 ) phi
  real ( kind = 8 ) phi1
  real ( kind = 8 ) t

  phi = t / ( 1.0D+00 - t * t )
  phi1 = ( t * t + 1.0D+00 ) / ( t * t - 1.0D+00 )**2

  return
end
function t_function ( y )

!*****************************************************************************80
!
!! T_FUNCTION solves Y = T * log ( T ) for T, given nonnegative Y.
!
!  Discussion:
!
!    This routine evaluates the function
!
!      t = t ( y )
!
!    that is the inverse of the function
!
!      y = t * ln ( t )
!
!    for nonnegative Y.  An approximation is used, and the accuracy is about
!    one percent.  For the approximation used, see the reference.
!
!  Reference:
!
!    Walter Gautschi,
!    Computational Aspects of Three-term Recurrence Relations,
!    SIAM Review,
!    Volume 9, 1967, pages 24-82.
!
!  Parameters:
!
!    Input, real Y, the value of the forward function.
!
!    Output, real T_FUNCTION, the value of the inverse function.
!
  implicit none

  real p
  real t_function
  real y
  real z

  if ( y <= 10.0E+00 ) then
    p = 0.000057941E+00 * y - 0.00176148E+00
    p = y * p + 0.0208645E+00
    p = y * p - 0.129013E+00
    p = y * p + 0.85777E+00
    t_function = y * p + 1.0125E+00
  else
    z = log ( y ) - 0.775E+00
    p = ( 0.775E+00 - log ( z ) ) / ( 1.0E+00 + z )
    p = 1.0E+00 / ( 1.0E+00 + p )
    t_function = y * p / z
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
subroutine tr_r4 ( t, phi, phi1 )

!*****************************************************************************80
!
!! TR_R4 maps T in [-1,1] to X in [0,oo).
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, real T.
!
!    Output, real PHI, the value of X(T).
!
!    Output, real PHI1, the derivative of the mapping.
!
  implicit none

  real phi
  real phi1
  real t

  phi = ( 1.0E+00 + t ) / ( 1.0E+00 - t )
  phi1 = 2.0E+00 / ( t - 1.0E+00 )**2

  return
end
subroutine tr_r8 ( t, phi, phi1 )

!*****************************************************************************80
!
!! TR_R8 maps T in [-1,1] to X in [0,oo).
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T.
!
!    Output, real ( kind = 8 ) PHI, the value of X(T).
!
!    Output, real ( kind = 8 ) PHI1, the derivative of the mapping.
!
  implicit none

  real ( kind = 8 ) phi
  real ( kind = 8 ) phi1
  real ( kind = 8 ) t

  phi = ( 1.0D+00 + t ) / ( 1.0D+00 - t )
  phi1 = 2.0D+00 / ( t - 1.0D+00 )**2

  return
end
