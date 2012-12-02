function alngam ( x, ifault )

!*****************************************************************************80
!
!! ALNGAM computes the logarithm of the Gamma function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 1999
!
!  Author:
!
!    Original FORTRAN77 version by Allan Mcleod.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Allan Mcleod,
!    Algorithm AS 245:
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Parameters:
!
!    Input, real X, the argument of the Gamma function.
!
!    Output, integer IFAULT, error flag.
!    0, no error occurred.
!    1, X is less than or equal to 0.
!    2, X is too big.
!
!    Output, real ALNGAM, the logarithm of the gamma function of X.
!
  implicit none

  real alngam
  real, parameter :: alr2pi = 0.918938533204673E+00
  integer ifault
  real, parameter, dimension ( 9 ) :: r1 = (/ &
     -2.66685511495E+00, &
     -24.4387534237E+00, &
     -21.9698958928E+00, &
      11.1667541262E+00, &
       3.13060547623E+00, &
       0.607771387771E+00, &
      11.9400905721E+00, &
      31.4690115749E+00, &
      15.2346874070E+00 /)
  real, parameter, dimension ( 9 ) :: r2 = (/ &
    -78.3359299449E+00, &
    -142.046296688E+00, &
     137.519416416E+00, &
     78.6994924154E+00, &
     4.16438922228E+00, &
     47.0668766060E+00, &
     313.399215894E+00, &
     263.505074721E+00, &
     43.3400022514E+00 /)
  real, parameter, dimension ( 9 ) :: r3 = (/ &
    -2.12159572323E+05, &
     2.30661510616E+05, &
     2.74647644705E+04, &
    -4.02621119975E+04, &
    -2.29660729780E+03, &
    -1.16328495004E+05, &
    -1.46025937511E+05, &
    -2.42357409629E+04, &
    -5.70691009324E+02 /)
  real, parameter, dimension ( 5 ) :: r4 = (/ &
     0.279195317918525E+00, &
     0.4917317610505968E+00, &
     0.0692910599291889E+00, &
     3.350343815022304E+00, &
     6.012459259764103E+00 /)
  real x
  real x1
  real x2
  real, parameter :: xlge = 5.10E+06
  real, parameter :: xlgst = 1.0E+30
  real y

  alngam = 0.0E+00
  ifault = 0
!
!  X <= 0.0
!
  if ( x <= 0.0E+00 ) then

    ifault = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ALNGAM - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
!
!  0 < X < 0.5
!
  else if ( x < 0.5E+00 ) then

    alngam = - log ( x )
    y = x + 1.0E+00

    if ( y == 1.0E+00 ) then
      return
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
!
!  0.5 <= X < 1.5
!
  else if ( x < 1.5E+00 ) then

    alngam = ( ( x - 0.5E+00 ) - 0.5E+00 ) * (((( &
           r1(5)   * x &
         + r1(4) ) * x &
         + r1(3) ) * x &
         + r1(2) ) * x &
         + r1(1) ) / (((( &
                     x &
         + r1(9) ) * x &
         + r1(8) ) * x &
         + r1(7) ) * x &
         + r1(6) )
!
!  1.5 <= X < 4.0.
!
  else if ( x < 4.0E+00 ) then

    y = ( x - 1.0E+00 ) - 1.0E+00

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
!  4.0 <= X < 12.0.
!
  else if ( x < 12.0E+00 ) then

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
!  12.0 <= X <= XLGE
!
  else if ( x <= xlge ) then

    y = log ( x )
    alngam = x * ( y - 1.0E+00 ) - 0.5E+00 * y + alr2pi

    x1 = 1.0E+00 / x
    x2 = x1**2

    alngam = alngam + x1 * ( ( &
              r4(3)   * &
         x2 + r4(2) ) * &
         x2 + r4(1) ) / ( ( &
         x2 + r4(5) ) * &
         x2 + r4(4) )
!
!  XLGE < X < XLGST
!
  else if ( x < xlgst ) then

    y = log ( x )
    alngam = x * ( y - 1.0E+00 ) - 0.5E+00 * y + alr2pi
!
!  XLGST <= X
!
  else if ( xlgst <= x ) then

    ifault = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ALNGAM - Fatal error!'
    write ( *, '(a)' ) '  X is too large.'
    stop

  end if

  return
end
function alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 1999
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
!    Input, real X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none

  real, parameter :: a1 = 5.75885480458E+00
  real, parameter :: a2 = 2.62433121679E+00
  real, parameter :: a3 = 5.92885724438E+00
  real alnorm
  real, parameter :: b1 = -29.8213557807E+00
  real, parameter :: b2 = 48.6959930692E+00
  real, parameter :: c1 = -0.000000038052E+00
  real, parameter :: c2 = 0.000398064794E+00
  real, parameter :: c3 = -0.151679116635E+00
  real, parameter :: c4 = 4.8385912808E+00
  real, parameter :: c5 = 0.742380924027E+00
  real, parameter :: c6 = 3.99019417011E+00
  real, parameter :: con = 1.28E+00
  real, parameter :: d1 = 1.00000615302E+00
  real, parameter :: d2 = 1.98615381364E+00
  real, parameter :: d3 = 5.29330324926E+00
  real, parameter :: d4 = -15.1508972451E+00
  real, parameter :: d5 = 30.789933034E+00
  real, parameter :: ltone = 7.0E+00
  real, parameter :: p = 0.398942280444E+00
  real, parameter :: q = 0.39990348504E+00
  real, parameter :: r = 0.398942280385E+00
  logical up
  logical upper
  real, parameter :: utzero = 18.66E+00
  real x
  real y
  real z

  up = upper
  z = x

  if ( z < 0.0E+00 ) then
    up = .not. up
    z = - z
  end if
!
!  TAKE ANOTHER LOOK AT THIS SET OF CONDITIONS.
!
  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0E+00
    else
      alnorm = 1.0E+00
    end if

    return

  end if

  y = 0.5E+00 * z**2

  if ( z <= con ) then

    alnorm = 0.5E+00 - z * ( p - q * y &
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
    alnorm = 1.0E+00 - alnorm
  end if

  return
end
function alogam ( x, ifault )

!*****************************************************************************80
!
!! ALOGAM computes the logarithm of the Gamma function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by Malcolm Pike, David Hill.
!    FORTRAN90 version by John Burkardt.
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
!    Input, real X, the argument of the Gamma function.
!    X should be greater than 0.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, X <= 0.
!
!    Output, real ALOGAM, the logarithm of the Gamma function of X.
!
  implicit none

  real alogam
  real f
  integer ifault
  real x
  real y
  real z

  if ( x <= 0.0E+00 ) then
    ifault = 1
    alogam = 0.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ALOGAM - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  ifault = 0
  y = x

  if ( x < 7.0E+00 ) then

    f = 1.0E+00
    z = y

    do

      if ( 7.0E+00 <= z ) then
        exit
      end if

      f = f * z
      z = z + 1.0E+00

    end do

    y = z
    f = - log ( f )

  else

    f = 0.0E+00

  end if

  z = 1.0E+00 / y / y

  alogam = f + ( y - 0.5E+00 ) * log ( y ) - y &
       + 0.918938533204673E+00 + &
       ((( &
       - 0.000595238095238E+00   * z &
       + 0.000793650793651E+00 ) * z &
       - 0.002777777777778E+00 ) * z &
       + 0.083333333333333E+00 ) / y

  return
end
function digamma ( x )

!*****************************************************************************80
!
!! DIGAMMA: the digamma or Psi function = d ( LOG ( GAMMA ( X ) ) ) / dX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    Original FORTRAN77 version by Jose Bernardo.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jose Bernardo,
!    Algorithm AS 103:
!    Psi ( Digamma ) Function,
!    Applied Statistics,
!    Volume 25, Number 3, 1976, pages 315-317.
!
!  Parameters:
!
!    Input, real X, the argument of the digamma function.
!    0 < X.
!
!    Output, real DIGAMMA, the value of the digamma function at X.
!
  implicit none

  real, parameter :: c = 8.5E+00
  real, parameter :: d1 = -0.5772156649E+00
  real digamma
  real r
  real, parameter :: s = 0.00001E+00
  real, parameter :: s3 = 0.08333333333E+00
  real, parameter :: s4 = 0.0083333333333E+00
  real, parameter :: s5 = 0.003968253968E+00
  real x
  real y
!
!  The argument must be positive.
!
  if ( x <= 0.0E+00 ) then

    digamma = 0.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGAMMA - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
!
!  Use approximation if argument <= S.
!
  else if ( x <= s ) then

    digamma = d1 - 1.0E+00 / x
!
!  Reduce the argument to DIGAMMA(X + N) where C <= (X + N).
!
  else

    digamma = 0.0E+00
    y = x

    do

      if ( c <= y ) then
        exit
      end if

      digamma = digamma - 1.0E+00 / y
      y = y + 1.0E+00

    end do
!
!  Use Stirling's (actually de Moivre's) expansion if C < argument.
!
    r = 1.0E+00 / y / y

    digamma = digamma + log ( y ) - 0.5E+00 / y &
    - r * ( s3 - r * ( s4 - r * s5 ) )

  end if

  return
end
subroutine dirichlet_check ( n, a )

!*****************************************************************************80
!
!! DIRICHLET_CHECK checks the parameters of the Dirichlet PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
  implicit none

  integer n

  real a(n)
  integer i
  logical positive

  positive = .false.

  do i = 1, n

    if ( a(i) < 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIRICHLET_CHECK - Fatal error!'
      write ( *, '(a)' ) '  A(I) < 0.'
      write ( *, '(a,i6)' ) '  For I = ', i
      write ( *, '(a,g14.6)' ) '  A(I) = ', a(i)
      stop
    else if ( 0.0E+00 < a(i) ) then
      positive = .true.
    end if

  end do

  if ( .not. positive ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRICHLET_CHECK - Fatal error!'
    write ( *, '(a)' ) '  All parameters are zero!'
    stop
  end if

  return
end
subroutine dirichlet_estimate ( k, n, x, ix, init, alpha, rlogl, v, g, &
  niter, s, eps, work, ifault )

!*****************************************************************************80
!
!! DIRICHLET_ESTIMATE estimates the parameters of a Dirichlet distribution.
!
!  Discussion:
!
!    This routine requires several auxilliary routines:
!
!      ALOGAM (CACM algorithm 291 or AS 245),
!      DIGAMA (AS 103),
!      GAMMAD (AS 239),
!      PPCHI2 (AS 91),
!      TRIGAM (AS 121).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    Original FORTRAN77 version by A Naryanan.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    A. Naryanan,
!    Algorithm AS 266:
!    Maximum Likelihood Estimation of the Parameters of the
!    Dirichlet Distribution,
!    Applied Statistics,
!    Volume 40, Number 2, 1991, pages 365-374.
!
!  Parameters:
!
!    Input, integer K, the number of parameters.
!    2 <= K.
!
!    Input, integer N, the number of observations.
!    K < N.
!
!    Input, real X(IX,K), contains the N by K array of samples
!    from the distribution.  X(I,J) is the J-th component of
!    the I-th sample.
!
!    Input, integer IX, the leading dimension of the array X.
!    N <= IX.
!
!    Input, integer INIT, specifies how the parameter estimates
!    are to be initialized:
!    1, use the method of moments;
!    2, initialize each ALPHA to the minimum of X;
!    otherwise, the input values of ALPHA already contain estimates.
!
!    Input/output, real ALPHA(K).
!    On input, if INIT is not 1 or 2, then ALPHA must contain
!    initial estimates for the parameters.
!    On output, with IFAULT = 0, ALPHA contains the computed
!    estimates for the parameters.
!
!    Output, real RLOGL, the value of the log-likelihood function
!    at the solution point.
!
!    Output, real V(K*(K+1)/2); V(J*(J-1)/2+I) contains the covariance
!    between ALPHA(I) and ALPHA(J), for I = 1 to J, J = 1 to K.
!
!    Output, real G(K), contains an estimate of the derivative of
!    the log-likelihood with respect to each component of ALPHA.
!
!    Output, integer NITER, contains the number of Newton-Raphson
!    iterations performed.
!
!    Output, real S, the value of the chi-squared statistic.
!
!    Output, real EPS, contains the probability that the chi-squared
!    statistic is less than S.
!
!    Workspace, real WORK(2*K).
!
!    Output, integer IFAULT, error indicator.
!    0, no error, the results were computed successfully;
!    1, K < 2;
!    2, N <= K;
!    3, IX < N;
!    4, if X(I,J) <= 0 for any I or J, or if
!       ABS ( Sum ( 1 <= J <= K ) X(I,J) - 1 ) >= GAMMA = 0.001;
!    5, if IFAULT is returned nonzero from the chi-square
!       routine PPCHI2;
!    6, if ALPHA(J) <= 0 for any J during any step of the iteration;
!    7, if MAXIT iterations were carried out but convergence
!       was not achieved.
!
  implicit none

  integer ix
  integer k

  real alngam
  real alpha(k)
  real, parameter :: alpha_min = 0.00001E+00
  real an
  real beta
  real chi2
  real digamma
  real eps
  real g(k)
  real, parameter :: gamma = 0.0001E+00
  real gammad
  integer i
  integer i2
  integer ifault
  integer ifault2
  integer init
  integer, parameter :: it_max = 100
  integer it_num
  integer j
  integer kk
  integer n
  integer niter
  real ppchi2
  real rk
  real rlogl
  real s
  real sum1
  real sum2
  real temp
  real trigamma
  real v((k*(k+1))/2)
  real varp1
  real work(2*k)
  real x(ix,k)
  real x_min
  real x11
  real x12

  ifault2 = 0
!
!  Check the input arguments.
!
  if ( k < 2 ) then
    ifault = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRICHLET_ESTIMATE - Fatal error!'
    write ( *, '(a)' ) '  K < 2.'
    stop
  end if

  if ( n <= k ) then
    ifault = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRICHLET_ESTIMATE - Fatal error!'
    write ( *, '(a)' ) '  N <= K.'
    stop
  end if

  if ( ix < n ) then
    ifault = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRICHLET_ESTIMATE - Fatal error!'
    write ( *, '(a)' ) '  IX < N.'
    stop
  end if

  do i = 1, n

    do j = 1, k
      if ( x(i,j) <= 0.0E+00 ) then
        niter = i
        ifault = 4
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DIRICHLET_ESTIMATE - Fatal error!'
        write ( *, '(a)' ) '  X(I,J) <= 0.'
        stop
      end if
    end do

    sum2 = sum ( x(i,1:k) )

    if ( gamma <= abs ( sum2 - 1.0E+00 ) ) then
      ifault = 4
      niter = i
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIRICHLET_ESTIMATE - Fatal error!'
      write ( *, '(a)' ) '  Row I does not sum to 1.'
      stop
    end if

  end do

  ifault = 0

  an = real ( n )
  rk = real ( k )
  niter = 0
!
!  Calculate initial estimates using the method of moments.
!
  if ( init == 1 ) then

    do j = 1, k - 1
      alpha(j) = sum ( x(1:n,j) ) / an
    end do

    alpha(k) = 1.0E+00 - sum ( alpha(1:k-1) )

    x12 = 0.0E+00
    do i = 1, n
      x12 = x12 + x(i,1) ** 2
    end do

    x12 = x12 / an
    varp1 = x12 - alpha(1) ** 2

    x11 = ( alpha(1) - x12 ) / varp1
    alpha(1:k) = x11 * alpha(1:k)
!
!  Calculate initial estimates using Ronning's suggestion.
!
  else if ( init == 2 ) then

    x_min = x(1,1)
    do j = 1, k
      do i = 1, n
        x_min = min ( x_min, x(i,j) )
      end do
    end do

    x_min = max ( x_min, alpha_min )

    alpha(1:k) = x_min

  end if
!
!  Check whether any ALPHA's are negative or zero.
!
  do j = 1, k
    if ( alpha(j) <= 0.0E+00 ) then
      ifault = 6
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIRICHLET_ESTIMATE - Fatal error!'
      write ( *, '(a,i6)' ) '  For J = ', j
      write ( *, '(a,g14.6)' ) '  ALPHA(J) = ', alpha(j)
      write ( *, '(a)' ) '  but ALPHA(J) must be positive.'
      stop
    end if
  end do
!
!  Calculate N * log ( G(J) ) for J = 1,2,...,K and store in WORK array.
!
  do j = 1, k
    work(j) = sum ( log ( x(1:n,j) ) )
  end do
!
!  Call Algorithm AS 91 to compute CHI2, the chi-squared value.
!
  chi2 = ppchi2 ( gamma, rk, ifault )

  if ( ifault /= 0 ) then
    ifault = 5
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRICHLET_ESTIMATE - Fatal error!'
    write ( *, '(a)' ) '  PPCHI2 returns error code.'
    stop
  end if
!
!  Carry out the Newton iteration.
!
  do it_num = 1, it_max

    sum2 = sum ( alpha(1:k) )

    sum1 = 0.0E+00
    do j = 1, k
      work(k+j) = trigamma ( alpha(j) )
      sum1 = sum1 + 1.0E+00 / work(k+j)
    end do

    beta = trigamma ( sum2 )
    beta = an * beta / ( 1.0E+00 - beta * sum1 )

    temp = digamma ( sum2 )

    do j = 1, k
      g(j) = an * ( temp - digamma ( alpha(j) ) ) + work(j)
    end do
!
!  Calculate the lower triangle of the Variance-Covariance matrix V.
!
    sum2 = beta / an**2
    do i = 1, k
      do j = 1, i
        kk = j + ( i * ( i - 1 ) ) / 2
        v(kk) = sum2 / ( work(k+i) * work(k+j) )
        if ( i == j ) then
          v(kk) = v(kk) + 1.0E+00 / ( an * work(k+j) )
        end if
      end do
    end do
!
!  Postmultiply the Variance-Covariance matrix V by G and store
!  in the last K elements of WORK.
!
    do i = 1, k

      sum2 = 0.0E+00
      i2 = ( i * ( i - 1 ) ) / 2
      do j = 1, i - 1
        sum2 = sum2 + v(i2+j) * g(j)
      end do
      do j = i + 1, k
        sum2 = sum2 + v(i+(j*(j-1))/2) * g(j)
      end do

      work(k+i) = sum2 + v((i*(i+1))/2) * g(i)

    end do
!
!  Update the ALPHA's.
!
    niter = it_num

    do j = 1, k
      alpha(j) = alpha(j) + work(k+j)
      alpha(j) = max ( alpha(j), alpha_min )
    end do

    do j = 1, k
      if ( alpha(j) <= 0.0E+00 ) then
        ifault = 6
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DIRICHLET_ESTIMATE - Fatal error!'
        write ( *, '(a,i6)' ) '  Newton iteration ', it_num
        write ( *, '(a)' ) '  Computed ALPHA(J) <= 0.'
        write ( *, '(a,i6)' ) '  J = ', j
        write ( *, '(a,g14.6)' ) '  ALPHA(J) = ', alpha(j)
        stop
      end if
    end do
!
!  Test for convergence.
!
    s = dot_product ( g(1:k), work(k+1:k+k) )

    if ( s < chi2 ) then
      eps = gammad ( s / 2.0E+00, rk / 2.0E+00, ifault2 )

      sum2 = sum ( alpha(1:k) )

      rlogl = 0.0E+00
      do j = 1, k
        rlogl = rlogl + ( alpha(j) - 1.0E+00 ) * work(j) - an * &
          alngam ( alpha(j), ifault2 )
      end do

      rlogl = rlogl + an * alngam ( sum2, ifault2 )

      return

    end if

  end do

  ifault = 7

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DIRICHLET_ESTIMATE - Fatal error!'
  write ( *, '(a)' ) '  No convergence.'

  stop
end
subroutine dirichlet_mean ( n, a, mean )

!*****************************************************************************80
!
!! DIRICHLET_MEAN returns the means of the Dirichlet PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
!    Output, real MEAN(N), the means of the PDF.
!
  implicit none

  integer n

  real a(n)
  real mean(n)

  call dirichlet_check ( n, a )

  mean(1:n) = a(1:n) / sum ( a(1:n) )

  return
end
subroutine dirichlet_mix_mean ( comp_max, comp_num, elem_num, a, comp_weight, &
  mean )

!*****************************************************************************80
!
!! DIRICHLET_MIX_MEAN returns the means of a Dirichlet mixture PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer COMP_MAX, the leading dimension of A, which must
!    be at least COMP_NUM.
!
!    Input, integer COMP_NUM, the number of components in the Dirichlet
!    mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    Input, integer ELEM_NUM, the number of elements of an observation.
!
!    Input, real A(COMP_MAX,ELEM_NUM), the probabilities for element ELEM_NUM
!    in component COMP_NUM.
!    Each A(I,J) should be greater than or equal to 0.0.
!
!    Input, integer COMP_WEIGHT(COMP_NUM), the mixture weights of the densities.
!    These do not need to be normalized.  The weight of a given component is
!    the relative probability that that component will be used to generate
!    the sample.
!
!    Output, real MEAN(ELEM_NUM), the means for each element.
!
  implicit none

  integer comp_max
  integer comp_num
  integer elem_num

  real a(comp_max,elem_num)
  real a_sum(comp_num)
  integer comp_i
  real comp_weight(comp_num)
  real comp_weight_sum
  integer elem_i
  real mean(elem_num)
!
!  Check.
!
  do comp_i = 1, comp_num

    do elem_i = 1, elem_num
      if ( a(comp_i,elem_i) < 0.0E+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DIRICHLET_MIX_MEAN - Fatal error!'
        write ( *, '(a)' ) '  A(COMP,ELEM) < 0.'
        write ( *, '(a,i6)' ) '  COMP = ', comp_i
        write ( *, '(a,i6)' ) '  ELEM = ', elem_i
        write ( *, '(a,g14.6)' ) '  A(COMP,ELEM) = ', a(comp_i,elem_i)
        stop
      end if
    end do

  end do

  comp_weight_sum = sum ( comp_weight(1:comp_num) )

  do comp_i = 1, comp_num
    a_sum(comp_i) = 0.0E+00
    do elem_i = 1, elem_num
      a_sum(comp_i) = a_sum(comp_i) + a(comp_i,elem_i)
    end do
  end do

  do elem_i = 1, elem_num
    mean(elem_i) = 0.0E+00
    do comp_i = 1, comp_num
      mean(elem_i) = mean(elem_i) + comp_weight(comp_i) * a(comp_i,elem_i) &
        / a_sum(comp_i)
    end do
    mean(elem_i) = mean(elem_i) / comp_weight_sum
  end do

  return
end
subroutine dirichlet_mix_sample ( comp_max, comp_num, elem_num, a, &
  comp_weight, seed, comp, x )

!*****************************************************************************80
!
!! DIRICHLET_MIX_SAMPLE samples a Dirichlet mixture PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer COMP_MAX, the leading dimension of A, which must
!    be at least COMP_NUM.
!
!    Input, integer COMP_NUM, the number of components in the Dirichlet
!    mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    Input, integer ELEM_NUM, the number of elements of an observation.
!
!    Input, real A(COMP_MAX,ELEM_NUM), the probabilities for element ELEM_NUM
!    in component COMP_NUM.
!    Each A(I,J) should be greater than or equal to 0.0.
!
!    Input, integer COMP_WEIGHT(COMP_NUM), the mixture weights of the densities.
!    These do not need to be normalized.  The weight of a given component is
!    the relative probability that that component will be used to generate
!    the sample.
!
!    Input/output, integer SEED, a seed for the random number generator.
!
!    Output, integer COMP, the index of the component of the Dirichlet
!    mixture that was chosen to generate the sample.
!
!    Output, real X(ELEM_NUM), a sample of the PDF.
!
  implicit none

  integer comp_max
  integer comp_num
  integer elem_num

  real a(comp_max,elem_num)
  integer comp
  integer comp_i
  real comp_weight(comp_num)
  real comp_weight_sum
  integer elem_i
  real r
  real r4_uniform
  integer seed
  real sum2
  real x(elem_num)
  real x_sum
!
!  Check.
!
  do comp_i = 1, comp_num

    do elem_i = 1, elem_num
      if ( a(comp_i,elem_i) < 0.0E+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DIRICHLET_MIX_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  A(COMP,ELEM) < 0.'
        write ( *, '(a,i6)' ) '  COMP = ', comp_i
        write ( *, '(a,i6)' ) '  ELEM = ', elem_i
        write ( *, '(a,g14.6)' ) '  A(COMP,ELEM) = ', a(comp_i,elem_i)
        stop
      end if
    end do

  end do
!
!  Choose a particular density MIX.
!
  comp_weight_sum = sum ( comp_weight(1:comp_num) )

  r = r4_uniform ( 0.0E+00, comp_weight_sum, seed )

  comp = 1
  sum2 = 0.0E+00

  do while ( comp < comp_num )

    sum2 = sum2 + comp_weight(comp)

    if ( r <= sum2 ) then
      exit
    end if

    comp = comp + 1

  end do
!
!  Sample density COMP.
!
  do elem_i = 1, elem_num
    call gamma_sample ( a(comp,elem_i), 1.0E+00, seed, x(elem_i) )
  end do
!
!  Normalize the result.
!
  x(1:elem_num) = x(1:elem_num) / sum ( x(1:elem_num) )

  return
end
subroutine dirichlet_sample ( n, a, seed, x )

!*****************************************************************************80
!
!! DIRICHLET_SAMPLE samples the Dirichlet PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 169.
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be
!    positive.
!
!    Input/output, integer SEED, a seed for the random number generator.
!
!    Output, real X(N), a sample of the PDF.  The entries of X should
!    sum to 1.
!
  implicit none

  integer n

  real a(n)
  integer i
  integer seed
  real x(n)

  call dirichlet_check ( n, a )

  do i = 1, n
    call gamma_sample ( a(i), 1.0E+00, seed, x(i) )
  end do
!
!  Normalize the result.
!
  x(1:n) = x(1:n) / sum ( x(1:n) )

  return
end
subroutine dirichlet_variance ( n, a, variance )

!*****************************************************************************80
!
!! DIRICHLET_VARIANCE returns the variances of the Dirichlet PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
!    Output, real VARIANCE(N), the variances of the PDF.
!
  implicit none

  integer n

  real a(n)
  real a_sum
  integer i
  real variance(n)

  call dirichlet_check ( n, a )

  a_sum = sum ( a(1:n) )

  do i = 1, n
    variance(i) = a(i) * ( a_sum - a(i) ) / ( a_sum**2 * ( a_sum + 1.0E+00 ) )
  end do

  return
end
subroutine exponential_01_sample ( seed, x )

!*****************************************************************************80
!
!! EXPONENTIAL_01_SAMPLE samples the Exponential PDF with parameters 0, 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED, a seed for the random number generator.
!
!    Output, real X, a sample of the PDF.
!
  implicit none

  real a
  real b
  real cdf
  real r4_uniform
  integer seed
  real x

  cdf = r4_uniform ( 0.0E+00, 1.0E+00, seed )

  a = 0.0E+00
  b = 1.0E+00
  call exponential_cdf_inv ( cdf, a, b, x )

  return
end
subroutine exponential_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! EXPONENTIAL_CDF_INV inverts the Exponential CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real X, the corresponding argument.
!
  implicit none

  real a
  real b
  real cdf
  real x
!
!  Check.
!
  if ( b <= 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXPONENTIAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.0'
    stop
  end if

  if ( cdf < 0.0E+00 .or. 1.0E+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXPONENTIAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a - b * log ( 1.0E+00 - cdf )

  return
end
function gamain ( x, p, ifault )

!*****************************************************************************80
!
!! GAMAIN computes the incomplete Gamma ratio.
!
!  Discussion:
!
!    A series expansion is used if X < P or X <= 1.  Otherwise, a
!    continued fraction approximation is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by G Bhattacharjee.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    G Bhattacharjee,
!    The Incomplete Gamma Integral,
!    Algorithm AS 32,
!    Applied Statistics,
!    Volume 19, Number 3, pages 285-287, 1970.
!
!  Parameters:
!
!    Input, real X, P, the parameters of the incomplete gamma ratio.
!    0 <= X, and 0 < P.
!
!    Output, integer IFAULT, error flag.
!    0, no errors.
!    1, P <= 0.
!    2, X < 0.
!    3, underflow.
!    4, error return from the Log Gamma function.
!
  implicit none

  real a
  real, parameter :: acu = 1.0E-08
  real alngam
  real an
  real arg
  real b
  real dif
  real factor
  real g
  real gamain
  real gin
  integer i
  integer ifault
  real, parameter :: oflo = 1.0E+37
  real p
  real pn(6)
  real rn
  real term
  real, parameter :: uflo = 1.0E-37
  real x
!
!  Check the input.
!
  if ( p <= 0.0E+00 ) then
    ifault = 1
    gamain = 0.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMAIN - Fatal error!'
    write ( *, '(a)' ) '  P <= 0.'
    stop
  end if

  if ( x < 0.0E+00 ) then
    ifault = 2
    gamain = 0.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMAIN - Fatal error!'
    write ( *, '(a)' ) '  X < 0.'
    stop
  end if

  if ( x == 0.0E+00 ) then
    ifault = 0
    gamain = 0.0E+00
    return
  end if

  g = alngam ( p, ifault )

  if ( ifault /= 0 ) then
    ifault = 4
    gamain = 0.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMAIN - Fatal error!'
    write ( *, '(a)' ) '  ALNGAM returns error code.'
    stop
  end if

  arg = p * log ( x ) - x - g

  if ( arg < log ( uflo ) ) then
    ifault = 3
    gamain = 0.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMAIN - Fatal error!'
    write ( *, '(a)' ) '  Underflow.'
    stop
  end if

  ifault = 0
  factor = exp ( arg )
!
!  Calculation by series expansion.
!
  if ( x <= 1.0E+00 .or. x < p ) then

    gin = 1.0E+00
    term = 1.0E+00
    rn = p

    do

      rn = rn + 1.0E+00
      term = term * x / rn
      gin = gin + term

      if ( term <= acu ) then
        exit
      end if

    end do

    gamain = gin * factor / p
    return

  end if
!
!  Calculation by continued fraction.
!
  a = 1.0E+00 - p
  b = a + x + 1.0E+00
  term = 0.0E+00

  pn(1) = 1.0E+00
  pn(2) = x
  pn(3) = x + 1.0E+00
  pn(4) = x * b

  gin = pn(3) / pn(4)

  do

    do

      a = a + 1.0E+00
      b = b + 2.0E+00
      term = term + 1.0E+00
      an = a * term
      do i = 1, 2
        pn(i+4) = b * pn(i+2) - an * pn(i)
      end do

      if ( pn(6) /= 0.0E+00 ) then

        rn = pn(5) / pn(6)
        dif = abs ( gin - rn )
!
!  Absolute error tolerance satisfied?
!
        if ( dif <= acu ) then
!
!  Relative error tolerance satisfied?
!
          if ( dif <= acu * rn ) then
            gamain = 1.0E+00 - factor * gin
            return
          end if

        end if

        gin = rn

      end if

      do i = 1, 4
        pn(i) = pn(i+2)
      end do

      if ( oflo <= abs ( pn(5) ) ) then
        exit
      end if

    end do

    do i = 1, 4
      pn(i) = pn(i) / oflo
    end do

  end do

end
function gamma_log ( x )

!*****************************************************************************80
!
!! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the
!    Gamma Function,
!    Mathematics of Computation,
!    Volume 21, Number 98, April 1967, pages 198-203.
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thacher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real X, the argument of the Gamma function.  X must be positive.
!
!    Output, real GAMMA_LOG, the logarithm of the Gamma function of X.
!    If X <= 0.0, or if overflow would occur, the program returns the
!    largest representable floating point number.
!
!
!  Explanation of machine-dependent constants
!
!  BETA   - radix for the floating-point representation.
!
!  MAXEXP - the smallest positive power of BETA that overflows.
!
!  XBIG   - largest argument for which LN(GAMMA(X)) is representable
!           in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!  FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!  Approximate values for some important machines are:
!
!                            BETA      MAXEXP         XBIG
!
!  CRAY-1        (S.P.)        2        8191       9.62E+2461
!  Cyber 180/855
!    under NOS   (S.P.)        2        1070       1.72E+319
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)        2         128       4.08E+36
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!  IBM 3033      (D.P.)       16          63       4.29D+73
!  VAX D-Format  (D.P.)        2         127       2.05D+36
!  VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                          FRTBIG
!
!  CRAY-1        (S.P.)   3.13E+615
!  Cyber 180/855
!    under NOS   (S.P.)   6.44E+79
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   1.42E+9
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)   2.25D+76
!  IBM 3033      (D.P.)   2.56D+18
!  VAX D-Format  (D.P.)   1.20D+9
!  VAX G-Format  (D.P.)   1.89D+76
!
  implicit none

  real, parameter :: d1 = - 5.772156649015328605195174E-01
  real, parameter :: d2 =   4.227843350984671393993777E-01
  real, parameter :: d4 =   1.791759469228055000094023E+00
  real, parameter :: FRTBIG = 1.42E+09
  real, parameter :: PNT68 = 0.6796875E+00
  real, parameter :: SQRTPI = 0.9189385332046727417803297E+00
  real, parameter :: XBIG = 4.08E+36

  real, parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728E-03, &
     8.4171387781295E-04, &
    -5.952379913043012E-04, &
     7.93650793500350248E-04, &
    -2.777777777777681622553E-03, &
     8.333333333333333331554247E-02, &
     5.7083835261E-03 /)
  real corr
  integer i
  real gamma_log
  real, parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888E+00, &
    2.018112620856775083915565E+02, &
    2.290838373831346393026739E+03, &
    1.131967205903380828685045E+04, &
    2.855724635671635335736389E+04, &
    3.848496228443793359990269E+04, &
    2.637748787624195437963534E+04, &
    7.225813979700288197698961E+03 /)
  real, parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064E+00, &
    5.424138599891070494101986E+02, &
    1.550693864978364947665077E+04, &
    1.847932904445632425417223E+05, &
    1.088204769468828767498470E+06, &
    3.338152967987029735917223E+06, &
    5.106661678927352456275255E+06, &
    3.074109054850539556250927E+06 /)
  real, parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062E+04, &
    2.426813369486704502836312E+06, &
    1.214755574045093227939592E+08, &
    2.663432449630976949898078E+09, &
    2.940378956634553899906876E+10, &
    1.702665737765398868392998E+11, &
    4.926125793377430887588120E+11, &
    5.606251856223951465078242E+11 /)
  real, parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036E+01, &
    1.113332393857199323513008E+03, &
    7.738757056935398733233834E+03, &
    2.763987074403340708898585E+04, &
    5.499310206226157329794414E+04, &
    6.161122180066002127833352E+04, &
    3.635127591501940507276287E+04, &
    8.785536302431013170870835E+03 /)
  real, parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942E+02, &
    7.765049321445005871323047E+03, &
    1.331903827966074194402448E+05, &
    1.136705821321969608938755E+06, &
    5.267964117437946917577538E+06, &
    1.346701454311101692290052E+07, &
    1.782736530353274213975932E+07, &
    9.533095591844353613395747E+06 /)
  real, parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843E+03, &
    6.393885654300092398984238E+05, &
    4.135599930241388052042842E+07, &
    1.120872109616147941376570E+09, &
    1.488613728678813811542398E+10, &
    1.016803586272438228077304E+11, &
    3.417476345507377132798597E+11, &
    4.463158187419713286462081E+11 /)
  real res
  real x
  real xden
  real xm1
  real xm2
  real xm4
  real xnum
  real xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0E+00 .or. XBIG < x ) then
    gamma_log = huge ( gamma_log )
    return
  end if

  if ( x <= epsilon ( x ) ) then

    res = - log ( x )

  else if ( x <= 1.5E+00 ) then

    if ( x < PNT68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0E+00
      xm1 = ( x - 0.5E+00 ) - 0.5E+00
    end if

    if ( x <= 0.5E+00 .or. PNT68 <= x ) then

      xden = 1.0E+00
      xnum = 0.0E+00

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5E+00 ) - 0.5E+00
      xden = 1.0E+00
      xnum = 0.0E+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0E+00 ) then

    xm2 = x - 2.0E+00
    xden = 1.0E+00
    xnum = 0.0E+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0E+00 ) then

    xm4 = x - 4.0E+00
    xden = - 1.0E+00
    xnum = 0.0E+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0E+00

    if ( x <= FRTBIG ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + SQRTPI - 0.5E+00 * corr
    res = res + x * ( corr - 1.0E+00 )

  end if

  gamma_log = res

  return
end
subroutine gamma_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! GAMMA_SAMPLE samples the Gamma PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 1999
!
!  Author:
!
!    Original FORTRAN77 version by Joachim Ahrens, Ulrich Dieter.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Joachim Ahrens, Ulrich Dieter,
!    Computer Methods for Sampling from Gamma, Beta, Poisson and
!    Binomial Distributions,
!    Computing,
!    Volume 12, Number 3, September 1974, pages 223-246.
!
!    Joachim Ahrens, Ulrich Dieter,
!    Generating Gamma Variates by a Modified Rejection Technique,
!    Communications of the ACM,
!    Volume 25, Number 1, January 1982, pages 47-54.
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input/output, integer SEED, a seed for the random number generator.
!
!    Output, real X, a sample of the PDF.
!
  implicit none

  real, parameter :: a1 =  0.3333333E+00
  real, parameter :: a2 = -0.2500030E+00
  real, parameter :: a3 =  0.2000062E+00
  real, parameter :: a4 = -0.1662921E+00
  real, parameter :: a5 =  0.1423657E+00
  real, parameter :: a6 = -0.1367177E+00
  real, parameter :: a7 =  0.1233795E+00
  real, parameter :: e1 = 1.0E+00
  real, parameter :: e2 = 0.4999897E+00
  real, parameter :: e3 = 0.1668290E+00
  real, parameter :: e4 = 0.0407753E+00
  real, parameter :: e5 = 0.0102930E+00
  real, parameter :: q1 =  0.04166669E+00
  real, parameter :: q2 =  0.02083148E+00
  real, parameter :: q3 =  0.00801191E+00
  real, parameter :: q4 =  0.00144121E+00
  real, parameter :: q5 = -0.00007388E+00
  real, parameter :: q6 =  0.00024511E+00
  real, parameter :: q7 =  0.00024240E+00

  real a
  real b
  real bcoef
  real c
  real d
  real e
  real p
  real q
  real q0
  real r
  real r4_uniform
  real s
  real s2
  real si
  integer seed
  real t
  real u
  real v
  real w
  real x
!
!  Allow A = 0.
!
  if ( a == 0.0E+00 ) then
    x = 0.0E+00
    return
  end if
!
!  A < 1.
!
  if ( a < 1.0E+00 ) then

    do

      p = r4_uniform ( 0.0E+00, 1.0E+00, seed )
      p = ( 1.0E+00 + 0.3678794E+00 * a ) * p

      call exponential_01_sample ( seed, e )

      if ( 1.0E+00 <= p ) then

        x = - log ( ( 1.0E+00 + 0.3678794E+00 * a - p ) / a )

        if ( ( 1.0E+00 - a ) * log ( x ) <= e ) then
          x = x / b
          return
        end if

      else

        x = exp ( log ( p ) / a )

        if ( x <= e ) then
          x = x / b
          return
        end if

      end if

    end do
!
!  1 <= A.
!
  else

    s2 = a - 0.5E+00
    s = sqrt ( a - 0.5E+00 )
    d = sqrt ( 32.0E+00 ) - 12.0E+00 * sqrt ( a - 0.5E+00 )

    call normal_01_sample ( seed, t )
    x = ( sqrt ( a - 0.5E+00 ) + 0.5E+00 * t )**2

    if ( 0.0E+00 <= t ) then
      x = x / b
      return
    end if

    u = r4_uniform ( 0.0E+00, 1.0E+00, seed )

    if ( d * u <= t**3 ) then
      x = x / b
      return
    end if

    r = 1.0E+00 / a
    q0 = ( ( ( ( ( ( &
           q7   * r &
         + q6 ) * r &
         + q5 ) * r &
         + q4 ) * r &
         + q3 ) * r &
         + q2 ) * r &
         + q1 ) * r

    if ( a <= 3.686E+00 ) then
      bcoef = 0.463E+00 + s - 0.178E+00 * s2
      si = 1.235E+00
      c = 0.195E+00 / s - 0.079E+00 + 0.016E+00 * s
    else if ( a <= 13.022E+00 ) then
      bcoef = 1.654E+00 + 0.0076E+00 * s2
      si = 1.68E+00 / s + 0.275E+00
      c = 0.062E+00 / s + 0.024E+00
    else
      bcoef = 1.77E+00
      si = 0.75E+00
      c = 0.1515E+00 / s
    end if

    if ( 0.0E+00 < sqrt ( a - 0.5E+00 ) + 0.5E+00 * t ) then

      v = 0.5E+00 * t / s

      if ( 0.25E+00 < abs ( v ) ) then
        q = q0 - s * t + 0.25E+00 * t**2 + 2.0E+00 * s2 * log ( 1.0E+00 + v )
      else
        q = q0 + 0.5E+00 * t**2 * ( ( ( ( ( ( &
               a7   * v &
             + a6 ) * v &
             + a5 ) * v &
             + a4 ) * v &
             + a3 ) * v &
             + a2 ) * v &
             + a1 ) * v
      end if

      if ( log ( 1.0E+00 - u ) <= q ) then
        x = x / b
        return
      end if

    end if

    do

      call exponential_01_sample ( seed, e )

      u = r4_uniform ( 0.0E+00, 1.0E+00, seed )
      u = 2.0E+00 * u - 1.0E+00
      t = bcoef + sign ( si * e, u )

      if ( - 0.7187449E+00 <= t ) then

        v = 0.5E+00 * t / s

        if ( 0.25E+00 < abs ( v ) ) then
          q = q0 - s * t + 0.25E+00 * t**2 + 2.0E+00 * s2 * log ( 1.0E+00 + v )
        else
          q = q0 + 0.5E+00 * t**2 * ( ( ( ( ( ( &
                 a7   * v &
               + a6 ) * v &
               + a5 ) * v &
               + a4 ) * v &
               + a3 ) * v &
               + a2 ) * v &
               + a1 ) * v
        end if

        if ( 0.0E+00 < q ) then

          if ( 0.5E+00 < q ) then
            w = exp ( q ) - 1.0E+00
          else
            w = ( ( ( ( &
                   e5   * q &
                 + e4 ) * q &
                 + e3 ) * q &
                 + e2 ) * q &
                 + e1 ) * q
          end if

          if ( c * abs ( u ) <= w * exp ( e - 0.5E+00 * t**2 ) ) then
            x = ( s + 0.5E+00 * t )**2 / b
            return
          end if

        end if

      end if

    end do

  end if

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
!      ALOGAM = logarithm of the gamma function,
!      ALNORM = algorithm AS66
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, real X, P, the parameters of the incomplete gamma ratio.
!    0 <= X, and 0 < P.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, X < 0 or P <= 0.
!
!    Output, real GAMMAD, the value of the incomplete Gamma integral.
!
  implicit none

  real, parameter :: elimit = - 88.0E+00
  real, parameter :: oflo = 1.0E+37
  real, parameter :: plimit = 1000.0E+00
  real, parameter :: tol = 1.0E-14
  real, parameter :: xbig = 1.0E+08

  real a
  real alnorm
  real alngam
  real an
  real arg
  real b
  real c
  real gammad
  integer ifault
  real p
  real pn1
  real pn2
  real pn3
  real pn4
  real pn5
  real pn6
  real rn
  logical upper
  real x

  gammad = 0.0E+00
!
!  Check the input.
!
  if ( x < 0.0E+00 ) then
    ifault = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMAD - Fatal error!'
    write ( *, '(a)' ) '  X < 0.'
    stop
  end if

  if ( p <= 0.0E+00 ) then
    ifault = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMAD - Fatal error!'
    write ( *, '(a)' ) '  P <= 0.'
    stop
  end if

  ifault = 0

  if ( x == 0.0E+00 ) then
    gammad = 0.0E+00
    return
  end if
!
!  If P is large, use a normal approximation.
!
  if ( plimit < p ) then

    pn1 = 3.0E+00 * sqrt ( p ) &
    * ( ( x / p )**( 1.0E+00 / 3.0E+00 ) + 1.0E+00 / ( 9.0E+00 * p ) - 1.0E+00 )

    upper = .false.
    gammad = alnorm ( pn1, upper )
    return

  end if
!
!  If X is large set GAMMAD = 1.
!
  if ( xbig < x ) then
    gammad = 1.0E+00
    return
  end if
!
!  Use Pearson's series expansion.
!  (Note that P is not large enough to force overflow in ALOGAM).
!  No need to test IFAULT on exit since 0 < P.
!
  if ( x <= 1.0E+00 .or. x < p ) then

    arg = p * log ( x ) - x - alngam ( p + 1.0E+00, ifault )
    c = 1.0E+00
    gammad = 1.0E+00
    a = p

    do

      a = a + 1.0E+00
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
      gammad = 0.0E+00
    end if
!
!  Use a continued fraction expansion.
!
  else

    arg = p * log ( x ) - x - alngam ( p, ifault )
    a = 1.0E+00 - p
    b = a + x + 1.0E+00
    c = 0.0E+00
    pn1 = 1.0E+00
    pn2 = x
    pn3 = x + 1.0E+00
    pn4 = x * b
    gammad = pn3 / pn4

    do

      a = a + 1.0E+00
      b = b + 2.0E+00
      c = c + 1.0E+00
      an = a * c
      pn5 = b * pn3 - an * pn1
      pn6 = b * pn4 - an * pn2

      if ( pn6 /= 0.0E+00 ) then

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
      gammad = 1.0E+00 - exp ( arg )
    else
      gammad = 1.0E+00
    end if

  end if

  return
end
function gammds ( x, p, ifault )

!*****************************************************************************80
!
!! GAMMDS computes the incomplete Gamma integral.
!
!  Discussion:
!
!    The parameters must be positive.  An infinite series is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by Chi-Leung Lau.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Chi-Leung Lau,
!    Algorithm AS 147:
!    A Simple Series for the Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 29, Number 1, 1980, pages 113-114.
!
!  Parameters:
!
!    Input, real X, P, the arguments of the incomplete Gamma integral.
!    X and P must be greater than 0.
!
!    Output, integer IFAULT, error flag.
!    0, no errors.
!    1, X <= 0 or P <= 0.
!    2, underflow during the computation.
!
!    Output, real GAMMDS, the value of the incomplete Gamma integral.
!
  implicit none

  real, parameter :: e = 1.0E-09
  real, parameter :: uflo = 1.0E-37

  real a
  real alngam
  real arg
  real c
  real f
  real gammds
  integer ifault
  integer ifault2
  real p
  real x
!
!  Check the input.
!
  if ( x <= 0.0E+00 ) then
    ifault = 1
    gammds = 0.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMDS - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  if ( p <= 0.0E+00 ) then
    ifault = 1
    gammds = 0.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMDS - Fatal error!'
    write ( *, '(a)' ) '  P <= 0.'
    stop
  end if
!
!  ALNGAM is the natural logarithm of the gamma function.
!
  ifault2 = 0
  arg = p * log ( x ) - alngam ( p + 1.0E+00, ifault2 ) - x

  if ( arg < log ( uflo ) ) then
    gammds = 0.0E+00
    ifault = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMDS - Warning!'
    write ( *, '(a)' ) '  Underflow.'
    return
  end if

  f = exp ( arg )

  if ( f == 0.0E+00 ) then
    gammds = 0.0E+00
    ifault = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMDS - Warning!'
    write ( *, '(a)' ) '  Underflow.'
    return
  end if

  ifault = 0
!
!  Series begins.
!
  c = 1.0E+00
  gammds = 1.0E+00
  a = p

  do

    a = a + 1.0E+00
    c = c * x / a
    gammds = gammds + c

    if ( c <= e * gammds ) then
      exit
    end if

  end do

  gammds = gammds * f

  return
end
function lngamma ( z, ier )

!*****************************************************************************80
!
!! LNGAMMA computes Log(Gamma(X)) using a Lanczos approximation.
!
!  Discussion:
!
!    This algorithm is not part of the Applied Statistics algorithms.
!    It is slower but gives 14 or more significant decimal digits
!    accuracy, except around X = 1 and X = 2.   The Lanczos series from
!    which this algorithm is derived is interesting in that it is a
!    convergent series approximation for the gamma function, whereas
!    the familiar series due to De Moivre (and usually wrongly called
!    Stirling's approximation) is only an asymptotic approximation, as
!    is the true and preferable approximation due to Stirling.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by Alan Miller.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Christian Lanczos,
!    A precision approximation of the gamma function,
!    SIAM Journal on Numerical Analysis, B,
!    Volume 1, 1964, pages 86-96.
!
!  Parameters:
!
!    Input, real Z, the argument of the Gamma function.
!
!    Output, integer IER, error flag.
!    0, no error occurred.
!    1, Z is less than or equal to 0.
!
!    Output, real LNGAMMA, the logarithm of the gamma function of Z.
!
  implicit none

  real, parameter :: lnsqrt2pi = 0.9189385332046727E+00

  real, parameter, dimension ( 9 ) :: a = (/ &
            0.9999999999995183E+00, &
          676.5203681218835E+00, &
       - 1259.139216722289E+00, &
          771.3234287757674E+00, &
        - 176.6150291498386E+00, &
           12.50734324009056E+00, &
          - 0.1385710331296526E+00, &
            0.9934937113930748E-05, &
            0.1659470187408462E-06 /)
  integer ier
  integer j
  real lngamma
  real tmp
  real z

  if ( z <= 0.0E+00 ) then
    ier = 1
    lngamma = 0.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LNGAMMA - Fatal error!'
    write ( *, '(a)' ) '  Z <= 0.'
    stop
  end if

  ier = 0

  lngamma = 0.0E+00
  tmp = z + 7.0E+00
  do j = 9, 2, -1
    lngamma = lngamma + a(j) / tmp
    tmp = tmp - 1.0E+00
  end do

  lngamma = lngamma + a(1)
  lngamma = log ( lngamma ) + lnsqrt2pi - ( z + 6.5E+00 ) + ( z - 0.5E+00 ) * &
    log ( z + 6.5E+00 )

  return
end
subroutine normal_01_sample ( seed, x )

!*****************************************************************************80
!
!! NORMAL_01_SAMPLE samples the standard Normal PDF.
!
!  Discussion:
!
!    The standard normal distribution has mean 0 and standard
!    deviation 1.
!
!    The Box-Muller method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED, a seed for the random number generator.
!
!    Output, real X, a sample of the PDF.
!
  implicit none

  real, parameter :: PI &
    = 3.14159265358979323846264338327950288419716939937510E+00

  integer iset
  integer seed
  real r4_uniform
  real v1
  real v2
  real x
  real xsave

  save iset
  save xsave

  data iset / 0 /
  data xsave / 0.0E+00 /

  if ( iset == 0 ) then

    v1 = r4_uniform ( 0.0E+00, 1.0E+00, seed )

    if ( v1 <= 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, '(a)' ) '  V1 <= 0.'
      write ( *, '(a,g14.6)' ) '  V1 = ', v1
      stop
    end if

    v2 = r4_uniform ( 0.0E+00, 1.0E+00, seed )

    if ( v2 <= 0.0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, '(a)' ) '  V2 <= 0.'
      write ( *, '(a,g14.6)' ) '  V2 = ', v2
      stop
    end if

    x = sqrt ( - 2.0E+00 * log ( v1 ) ) * cos ( 2.0E+00 * PI * v2 )

    xsave = sqrt ( - 2.0E+00 * log ( v1 ) ) * sin ( 2.0E+00 * PI * v2 )

    iset = 1

  else

    x = xsave
    iset = 0

  end if

  return
end
subroutine normp ( z, p, q, pdf )

!*****************************************************************************80
!
!! NORMP computes the cumulative density of the standard normal distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by Alan Miller.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thacher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real Z, divides the real line into two semi-infinite
!    intervals, over each of which the standard normal distribution
!    is to be integrated.
!
!    Output, real P, Q, the integrals of the standard normal
!    distribution over the intervals ( - Infinity, Z] and
!    [Z, + Infinity ), respectively.
!
!    Output, real PDF, the value of the standard normal distribution at Z.
!
  implicit none

  real, parameter :: cutoff = 7.071E+00
  real, parameter :: p0 = 220.2068679123761E+00
  real, parameter :: p1 = 221.2135961699311E+00
  real, parameter :: p2 = 112.0792914978709E+00
  real, parameter :: p3 = 33.91286607838300E+00
  real, parameter :: p4 = 6.373962203531650E+00
  real, parameter :: p5 = 0.7003830644436881E+00
  real, parameter :: p6 = 0.03526249659989109E+00
  real, parameter :: q0 = 440.4137358247522E+00
  real, parameter :: q1 = 793.8265125199484E+00
  real, parameter :: q2 = 637.3336333788311E+00
  real, parameter :: q3 = 296.5642487796737E+00
  real, parameter :: q4 = 86.78073220294608E+00
  real, parameter :: q5 = 16.06417757920695E+00
  real, parameter :: q6 = 1.755667163182642E+00
  real, parameter :: q7 = 0.08838834764831844E+00
  real, parameter :: root2pi = 2.506628274631001E+00

  real expntl
  real p
  real pdf
  real q
  real z
  real zabs

  zabs = abs ( z )
!
!  37 < |Z|.
!
  if ( 37.0E+00 < zabs ) then

    pdf = 0.0E+00
    p = 0.0E+00
!
!  |Z| <= 37.
!
  else

    expntl = exp ( - 0.5E+00 * zabs**2 )
    pdf = expntl / root2pi
!
!  |Z| < CUTOFF = 10 / sqrt(2).
!
    if ( zabs < cutoff ) then

      p = expntl * (((((( &
             p6   * zabs &
           + p5 ) * zabs &
           + p4 ) * zabs &
           + p3 ) * zabs &
           + p2 ) * zabs &
           + p1 ) * zabs &
           + p0 ) / ((((((( &
             q7   * zabs &
           + q6 ) * zabs &
           + q5 ) * zabs &
           + q4 ) * zabs &
           + q3 ) * zabs &
           + q2 ) * zabs &
           + q1 ) * zabs &
           + q0 )
!
!  CUTOFF <= |Z|.
!
    else

      p = pdf / ( &
           zabs + 1.0E+00 / ( &
           zabs + 2.0E+00 / ( &
           zabs + 3.0E+00 / ( &
           zabs + 4.0E+00 / ( &
           zabs + 0.65E+00 )))))

    end if

  end if

  if ( z < 0.0E+00 ) then
    q = 1.0E+00 - p
  else
    q = p
    p = 1.0E+00 - q
  end if

  return
end
subroutine nprob ( z, p, q, pdf )

!*****************************************************************************80
!
!! NPROB computes the cumulative density of the standard normal distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by AG Adams.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39:
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, 1969, pages 197-198.
!
!  Parameters:
!
!    Input, real Z, divides the real line into two semi-infinite
!    intervals, over each of which the standard normal distribution
!    is to be integrated.
!
!    Output, real P, Q, the integrals of the standard normal
!    distribution over the intervals ( - Infinity, Z] and
!    [Z, + Infinity ), respectively.
!
!    Output, real PDF, the value of the standard normal distribution at Z.
!
  implicit none

  real, parameter :: a0 = 0.5E+00
  real, parameter :: a1 = 0.398942280444E+00
  real, parameter :: a2 = 0.399903438504E+00
  real, parameter :: a3 = 5.75885480458E+00
  real, parameter :: a4 = 29.8213557808E+00
  real, parameter :: a5 = 2.62433121679E+00
  real, parameter :: a6 = 48.6959930692E+00
  real, parameter :: a7 = 5.92885724438E+00
  real, parameter :: b0 = 0.398942280385E+00
  real, parameter :: b1 = 0.000000038052E+00
  real, parameter :: b2 = 1.00000615302E+00
  real, parameter :: b3 = 0.000398064794E+00
  real, parameter :: b4 = 1.98615381364E+00
  real, parameter :: b5 = 0.151679116635E+00
  real, parameter :: b6 = 5.29330324926E+00
  real, parameter :: b7 = 4.8385912808E+00
  real, parameter :: b8 = 15.1508972451E+00
  real, parameter :: b9 = 0.742380924027E+00
  real, parameter :: b10 = 30.789933034E+00
  real, parameter :: b11 = 3.99019417011E+00

  real p
  real pdf
  real q
  real y
  real z
  real zabs

  zabs = abs ( z )
!
!  |Z| between 0 and 1.28
!
  if ( abs ( z ) <= 1.28E+00 ) then

    y = a0 * z**2
    pdf = exp ( - y ) * b0

    q = a0 - zabs * ( a1 - a2 * y &
         / ( y + a3 - a4 &
         / ( y + a5 + a6 &
         / ( y + a7 ))))
!
!  |Z| between 1.28 and 12.7
!
  else if ( abs ( z ) <= 12.7E+00 ) then

    y = a0 * z**2
    pdf = exp ( - y ) * b0

    q = pdf &
         / ( zabs - b1 + b2 &
         / ( zabs + b3 + b4 &
         / ( zabs - b5 + b6 &
         / ( zabs + b7 - b8 &
         / ( zabs + b9 + b10 &
         / ( zabs + b11 ))))))
!
!  Z far out in tail.
!
  else

    q = 0.0E+00
    pdf = 0.0E+00

  end if

  if ( z < 0.0E+00 ) then
    p = q
    q = 1.0E+00 - p
  else
    p = 1.0E+00 - q
  end if

  return
end
function ppchi2 ( p, v, ifault )

!*****************************************************************************80
!
!! PPCHI2 evaluates the percentage points of the chi-squared PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 1999
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
!    Input, real P, a value of the chi-squared cumulative probability
!    density function.
!    0.000002 <= P <= 0.999998.
!
!    Input, real V, the parameter of the chi-squared probability density
!    function.  0 < V.
!
!    Output, integer IFAULT, error flag.
!    0, no error detected.
!    1, P < PMIN or PMAX < P.
!    2, V <= 0.0.
!    3, an error occurred in the incomplete Gamma function routine.
!    4, the maximum number of iterations were taken without convergence.
!    5, an error occurred in the log Gamma routine
!
  implicit none

  real, parameter :: aa = 0.6931471806E+00
  real, parameter :: c1 = 0.01E+00
  real, parameter :: c2 = 0.222222E+00
  real, parameter :: c3 = 0.32E+00
  real, parameter :: c4 = 0.4E+00
  real, parameter :: c5 = 1.24E+00
  real, parameter :: c6 = 2.2E+00
  real, parameter :: c7 = 4.67E+00
  real, parameter :: c8 = 6.66E+00
  real, parameter :: c9 = 6.73E+00
  real, parameter :: c10 = 13.32E+00
  real, parameter :: c11 = 60.0E+00
  real, parameter :: c12 = 70.0E+00
  real, parameter :: c13 = 84.0E+00
  real, parameter :: c14 = 105.0E+00
  real, parameter :: c15 = 120.0E+00
  real, parameter :: c16 = 127.0E+00
  real, parameter :: c17 = 140.0E+00
  real, parameter :: c18 = 175.0E+00
  real, parameter :: c19 = 210.0E+00
  real, parameter :: c20 = 252.0E+00
  real, parameter :: c21 = 264.0E+00
  real, parameter :: c22 = 294.0E+00
  real, parameter :: c23 = 346.0E+00
  real, parameter :: c24 = 420.0E+00
  real, parameter :: c25 = 462.0E+00
  real, parameter :: c26 = 606.0E+00
  real, parameter :: c27 = 672.0E+00
  real, parameter :: c28 = 707.0E+00
  real, parameter :: c29 = 735.0E+00
  real, parameter :: c30 = 889.0E+00
  real, parameter :: c31 = 932.0E+00
  real, parameter :: c32 = 966.0E+00
  real, parameter :: c33 = 1141.0E+00
  real, parameter :: c34 = 1182.0E+00
  real, parameter :: c35 = 1278.0E+00
  real, parameter :: c36 = 1740.0E+00
  real, parameter :: c37 = 2520.0E+00
  real, parameter :: c38 = 5040.0E+00
  real, parameter :: e = 0.0000005E+00
  integer, parameter :: maxit = 20
  real, parameter :: pmax = 0.999998E+00
  real, parameter :: pmin = 0.000002E+00

  real a
  real alngam
  real b
  real c
  real ch
  real g
  real gammad
  integer i
  integer ifault
  integer ifault2
  real p
  real p1
  real p2
  real ppchi2
  real ppnd
  real q
  real s1
  real s2
  real s3
  real s4
  real s5
  real s6
  real t
  real v
  real x
  real xx

  ifault2 = 0
!
!  Check the input.
!
  if ( p < pmin ) then
    ifault = 1
    ppchi2 = - 1.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPCHI2 - Fatal error!'
    write ( *, '(a)' ) '  P < PMIN.'
    stop
  end if

  if ( pmax < p ) then
    ifault = 1
    ppchi2 = - 1.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPCHI2 - Fatal error!'
    write ( *, '(a)' ) '  PMAX < P.'
    stop
  end if

  if ( v <= 0.0E+00 ) then
    ifault = 2
    ppchi2 = - 1.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPCHI2 - Fatal error!'
    write ( *, '(a)' ) '  V <= 0.0.'
    stop
  end if

  ifault = 0
  xx = 0.5E+00 * v
  c = xx - 1.0E+00
!
!  Compute Log ( Gamma ( V/2 ) ).
!
  g = alngam ( v / 2.0E+00, ifault )

  if ( ifault /= 0 ) then
    ifault = 5
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPCHI2 - Fatal error!'
    write ( *, '(a)' ) '  ALNGAM returns error code.'
    stop
  end if
!
!  Starting approximation for small chi-squared.
!
  if ( v < - c5 * log ( p ) ) then

    ch = ( p * xx * exp ( g + xx * aa ) )**( 1.0E+00 / xx )

    if ( ch < e ) then
      ifault = 0
      ppchi2 = ch
      return
    end if
!
!  Starting approximation for V less than or equal to 0.32.
!
  else if ( v <= c3 ) then

    ch = c4
    a = log ( 1.0E+00 - p )

    do

      q = ch
      p1 = 1.0E+00 + ch * ( c7 + ch )
      p2 = ch * ( c9 + ch * ( c8 + ch ) )

      t = - 0.5E+00 + ( c7 + 2.0E+00 * ch ) / p1 &
        - ( c9 + ch * ( c10 + 3.0E+00 * ch ) ) / p2

      ch = ch - ( 1.0E+00 - exp ( a + g + 0.5E+00 * ch + c * aa ) * p2 / p1 ) &
        / t

      if ( abs ( q / ch - 1.0E+00 ) <= c1 ) then
        exit
      end if

    end do
!
!  Call to algorithm AS 111.
!  Note that P has been tested above.
!  AS 241 could be used as an alternative.
!
  else

    x = ppnd ( p, ifault2 )
!
!  Starting approximation using Wilson and Hilferty estimate.
!
    p1 = c2 / v
    ch = v * ( x * sqrt ( p1 ) + 1.0E+00 - p1 )**3
!
!  Starting approximation for P tending to 1.
!
    if ( c6 * v + 6.0E+00 < ch ) then
      ch = - 2.0E+00 * ( log ( 1.0E+00 - p ) - c * log ( 0.5E+00 * ch ) + g )
    end if

  end if
!
!  Call to algorithm AS 239 and calculation of seven term Taylor series.
!
  do i = 1, maxit

    q = ch
    p1 = 0.5E+00 * ch
    p2 = p - gammad ( p1, xx, ifault2 )

    if ( ifault2 /= 0 ) then
      ppchi2 = - 1.0E+00
      ifault = 3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PPCHI2 - Fatal error!'
      write ( *, '(a)' ) '  GAMMAD returns error code.'
      stop
    end if

    t = p2 * exp ( xx * aa + g + p1 - c * log ( ch ) )
    b = t / ch
    a = 0.5E+00 * t - b * c

    s1 = &
           ( c19 + a &
         * ( c17 + a &
         * ( c14 + a &
         * ( c13 + a &
         * ( c12 + a &
         *   c11 ) ) ) ) ) / c24

    s2 = &
           ( c24 + a &
         * ( c29 + a &
         * ( c32 + a &
         * ( c33 + a &
         *   c35 ) ) ) ) / c37

    s3 = &
           ( c19 + a &
         * ( c25 + a &
         * ( c28 + a * &
             c31 ) ) ) / c37

    s4 = &
           ( c20 + a &
         * ( c27 + a &
         *   c34 ) + c &
         * ( c22 + a &
         * ( c30 + a &
         *   c36 ) ) ) / c38

    s5 = ( c13 + c21 * a + c * ( c18 + c26 * a ) ) / c37

    s6 = ( c15 + c * ( c23 + c16 * c ) ) / c38

    ch = ch + t * ( 1.0E+00 + 0.5E+00 * t * s1 - b * c &
         * ( s1 - b &
         * ( s2 - b &
         * ( s3 - b &
         * ( s4 - b &
         * ( s5 - b &
         *   s6 ) ) ) ) ) )

    if ( e < abs ( q / ch - 1.0E+00 ) ) then
      ifault = 0
      ppchi2 = ch
      return
    end if

  end do

  ifault = 4
  ppchi2 = ch
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PPCHI2 - Warning!'
  write ( *, '(a)' ) '  Convergence not reached.'

  return
end
function ppnd ( p, ifault )

!*****************************************************************************80
!
!! PPND produces the normal deviate value corresponding to lower tail area = P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by JD Beasley, SG Springer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    JD Beasley, SG Springer,
!    Algorithm AS 111:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 26, Number 1, 1977, pages 118-121.
!
!  Parameters:
!
!    Input, real P, the value of the cumulative probability densitity function.
!    0 < P < 1.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, P <= 0 or 1 <= P.  PPND is returned as 0.
!
!    Output, real PPND, the normal deviate value with the property that
!    the probability of a standard normal deviate being less than or
!    equal to PPND is P.
!
  implicit none

  real, parameter :: a0 = 2.50662823884E+00
  real, parameter :: a1 = -18.61500062529E+00
  real, parameter :: a2 = 41.39119773534E+00
  real, parameter :: a3 = -25.44106049637E+00
  real, parameter :: b1 = -8.47351093090E+00
  real, parameter :: b2 = 23.08336743743E+00
  real, parameter :: b3 = -21.06224101826E+00
  real, parameter :: b4 = 3.13082909833E+00
  real, parameter :: c0 = -2.78718931138E+00
  real, parameter :: c1 = -2.29796479134E+00
  real, parameter :: c2 = 4.85014127135E+00
  real, parameter :: c3 = 2.32121276858E+00
  real, parameter :: d1 = 3.54388924762E+00
  real, parameter :: d2 = 1.63706781897E+00
  real, parameter :: split = 0.42E+00

  integer ifault
  real p
  real ppnd
  real r

  ifault = 0
!
!  0.08 < P < 0.92
!
  if ( abs ( p - 0.5E+00 ) <= split ) then

    r = ( p - 0.5E+00 )**2

    ppnd = ( p - 0.5E+00 ) * ( ( ( &
           a3   * r &
         + a2 ) * r &
         + a1 ) * r &
         + a0 ) / ( ( ( ( &
           b4   * r &
         + b3 ) * r &
         + b2 ) * r &
         + b1 ) * r &
         + 1.0E+00 )
!
!  P < 0.08 or 0.92 < P,
!  R = min ( P, 1-P )
!
  else if ( 0.0E+00 < p .and. p < 1.0E+00 ) then

    if ( 0.5E+00 < p ) then
      r = sqrt ( - log ( 1.0E+00 - p ) )
    else
      r = sqrt ( - log ( p ) )
    end if

    ppnd = ( ( ( &
           c3   * r &
         + c2 ) * r &
         + c1 ) * r &
         + c0 ) / ( ( &
           d2   * r &
         + d1 ) * r &
         + 1.0E+00 )

    if ( p < 0.5E+00 ) then
      ppnd = - ppnd
    end if
!
!  P <= 0.0 or 1.0 <= P.
!
  else

    ifault = 1
    ppnd = 0.0E+00

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPND - Warning!'
    write ( *, '(a)' ) '  P <= 0 or 1 <= P.'
    write ( *, '(a)' ) '  PPND value would be infinite.'

  end if

  return
end
function ppnd16 ( p, ifault )

!*****************************************************************************80
!
!! PPND16 produces the normal deviate value corresponding to lower tail area = P.
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
!    28 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    Algorithm AS 241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 477-484.
!
!  Parameters:
!
!    Input, double precision P, the value of the cumulative probability
!    densitity function.  0 < P < 1.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, P <= 0 or 1 <= P.
!
!    Output, double precision PPND16, the normal deviate value with the
!    property that the probability of a standard normal deviate being
!    less than or equal to PPND16 is P.
!
  implicit none

  double precision, parameter :: a0 = 3.3871328727963666080d+00
  double precision, parameter :: a1 = 1.3314166789178437745d+02
  double precision, parameter :: a2 = 1.9715909503065514427d+03
  double precision, parameter :: a3 = 1.3731693765509461125d+04
  double precision, parameter :: a4 = 4.5921953931549871457d+04
  double precision, parameter :: a5 = 6.7265770927008700853d+04
  double precision, parameter :: a6 = 3.3430575583588128105d+04
  double precision, parameter :: a7 = 2.5090809287301226727d+03
  double precision, parameter :: b1 = 4.2313330701600911252d+01
  double precision, parameter :: b2 = 6.8718700749205790830d+02
  double precision, parameter :: b3 = 5.3941960214247511077d+03
  double precision, parameter :: b4 = 2.1213794301586595867d+04
  double precision, parameter :: b5 = 3.9307895800092710610d+04
  double precision, parameter :: b6 = 2.8729085735721942674d+04
  double precision, parameter :: b7 = 5.2264952788528545610d+03
  double precision, parameter :: c0 = 1.42343711074968357734d+00
  double precision, parameter :: c1 = 4.63033784615654529590d+00
  double precision, parameter :: c2 = 5.76949722146069140550d+00
  double precision, parameter :: c3 = 3.64784832476320460504d+00
  double precision, parameter :: c4 = 1.27045825245236838258d+00
  double precision, parameter :: c5 = 2.41780725177450611770d-01
  double precision, parameter :: c6 = 2.27238449892691845833d-02
  double precision, parameter :: c7 = 7.74545014278341407640d-04
  double precision, parameter :: const1 = 0.180625d+00
  double precision, parameter :: const2 = 1.6d+00
  double precision, parameter :: d1 = 2.05319162663775882187d+00
  double precision, parameter :: d2 = 1.67638483018380384940d+00
  double precision, parameter :: d3 = 6.89767334985100004550d-01
  double precision, parameter :: d4 = 1.48103976427480074590d-01
  double precision, parameter :: d5 = 1.51986665636164571966d-02
  double precision, parameter :: d6 = 5.47593808499534494600d-04
  double precision, parameter :: d7 = 1.05075007164441684324d-09
  double precision, parameter :: e0 = 6.65790464350110377720d+00
  double precision, parameter :: e1 = 5.46378491116411436990d+00
  double precision, parameter :: e2 = 1.78482653991729133580d+00
  double precision, parameter :: e3 = 2.96560571828504891230d-01
  double precision, parameter :: e4 = 2.65321895265761230930d-02
  double precision, parameter :: e5 = 1.24266094738807843860d-03
  double precision, parameter :: e6 = 2.71155556874348757815d-05
  double precision, parameter :: e7 = 2.01033439929228813265d-07
  double precision, parameter :: f1 = 5.99832206555887937690d-01
  double precision, parameter :: f2 = 1.36929880922735805310d-01
  double precision, parameter :: f3 = 1.48753612908506148525d-02
  double precision, parameter :: f4 = 7.86869131145613259100d-04
  double precision, parameter :: f5 = 1.84631831751005468180d-05
  double precision, parameter :: f6 = 1.42151175831644588870d-07
  double precision, parameter :: f7 = 2.04426310338993978564d-15
  double precision, parameter :: split1 = 0.425d+00
  double precision, parameter :: split2 = 5.d+00

  integer ifault
  double precision p
  double precision ppnd16
  double precision q
  double precision r

  ifault = 0
  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q**2

    ppnd16 = q * ((((((( &
           a7   * r &
         + a6 ) * r &
         + a5 ) * r &
         + a4 ) * r &
         + a3 ) * r &
         + a2 ) * r &
         + a1 ) * r &
         + a0 ) / ((((((( &
           b7   * r &
         + b6 ) * r &
         + b5 ) * r &
         + b4 ) * r &
         + b3 ) * r &
         + b2 ) * r &
         + b1 ) * r &
         + 1.0D+00 )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

    if ( r <= 0.0D+00 ) then
      ifault = 1
      ppnd16 = 0.0D+00
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PPND16 - Warning!'
      write ( *, '(a)' ) '  P <= 0 or 1 <= P.'
      write ( *, '(a)' ) '  PPND16 value would be infinite.'
      return
    end if

    r = sqrt ( - log ( r ) )

    if ( r <= split2 ) then

      r = r - const2

      ppnd16 = ((((((( &
             c7   * r &
           + c6 ) * r &
           + c5 ) * r &
           + c4 ) * r &
           + c3 ) * r &
           + c2 ) * r &
           + c1 ) * r &
           + c0 ) / ((((((( &
             d7   * r &
           + d6 ) * r &
           + d5 ) * r &
           + d4 ) * r &
           + d3 ) * r &
           + d2 ) * r &
           + d1 ) * r &
           + 1.0D+00 )

    else

      r = r - split2

      ppnd16 = ((((((( &
             e7   * r &
           + e6 ) * r &
           + e5 ) * r &
           + e4 ) * r &
           + e3 ) * r &
           + e2 ) * r &
           + e1 ) * r &
           + e0 ) / ((((((( &
             f7   * r &
           + f6 ) * r &
           + f5 ) * r &
           + f4 ) * r &
           + f3 ) * r &
           + f2 ) * r &
           + f1 ) * r &
           + 1.0D+00 )

    end if

    if ( q < 0.0D+00 ) then
      ppnd16 = - ppnd16
    end if

  end if

  return
end
function ppnd7 ( p, ifault )

!*****************************************************************************80
!
!! PPND7 produces the normal deviate value corresponding to lower tail area = P.
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
!    31 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    Algorithm AS 241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 477-484.
!
!  Parameters:
!
!    Input, real P, the value of the cumulative probability densitity function.
!    0 < P < 1.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, P <= 0 or 1 <= P.
!
!    Output, real PPND7, the normal deviate value with the property that
!    the probability of a standard normal deviate being less than or
!    equal to PPND7 is P.
!
  implicit none

  real, parameter :: a0 = 3.3871327179E+00
  real, parameter :: a1 = 50.434271938E+00
  real, parameter :: a2 = 159.29113202E+00
  real, parameter :: a3 = 59.109374720E+00
  real, parameter :: b1 = 17.895169469E+00
  real, parameter :: b2 = 78.757757664E+00
  real, parameter :: b3 = 67.187563600E+00
  real, parameter :: c0 = 1.4234372777E+00
  real, parameter :: c1 = 2.7568153900E+00
  real, parameter :: c2 = 1.3067284816E+00
  real, parameter :: c3 = 0.17023821103E+00
  real, parameter :: const1 = 0.180625E+00
  real, parameter :: const2 = 1.6E+00
  real, parameter :: d1 = 0.73700164250E+00
  real, parameter :: d2 = 0.12021132975E+00
  real, parameter :: e0 = 6.6579051150E+00
  real, parameter :: e1 = 3.0812263860E+00
  real, parameter :: e2 = 0.42868294337E+00
  real, parameter :: e3 = 0.017337203997E+00
  real, parameter :: f1 = 0.24197894225E+00
  real, parameter :: f2 = 0.012258202635E+00
  real, parameter :: split1 = 0.425E+00
  real, parameter :: split2 = 5.0E+00

  integer ifault
  real p
  real ppnd7
  real q
  real r

  ifault = 0
  q = p - 0.5E+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q**2

    ppnd7 = q * ((( &
           a3   * r &
         + a2 ) * r &
         + a1 ) * r &
         + a0 ) / ((( &
           b3   * r &
         + b2 ) * r &
         + b1 ) * r &
         + 1.0E+00 )

  else

    if ( q < 0.0E+00 ) then
      r = p
    else
      r = 1.0E+00 - p
    end if

    if ( r <= 0.0E+00 ) then
      ifault = 1
      ppnd7 = 0.0E+00
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PPND7 - Fatal error!'
      write ( *, '(a)' ) '  P <= 0 or 1 <= P.'
      write ( *, '(a)' ) '  PPND7 value would be infinite.'
      return
    end if

    r = sqrt ( - log ( r ) )

    if ( r <= split2 ) then

      r = r - const2

      ppnd7 = ((( &
           c3   * r &
         + c2 ) * r &
         + c1 ) * r &
         + c0 ) / (( &
           d2   * r &
         + d1 ) * r &
         + 1.0E+00 )

    else

      r = r - split2

      ppnd7 = ((( &
           e3   * r &
         + e2 ) * r &
         + e1 ) * r &
         + e0 ) / (( &
           f2   * r &
         + f1 ) * r &
         + 1.0E+00 )

    end if

    if ( q < 0.0E+00 ) then
      ppnd7 = - ppnd7
    end if

  end if

  return
end
function psi ( xx )

!*****************************************************************************80
!
!! PSI evaluates the psi or digamma function, d/dx ln(gamma(x)).
!
!  Discussion:
!
!    The main computation involves evaluation of rational Chebyshev
!    approximations published in Math. Comp. 27, 123-127(1973) by
!    Cody, Strecok and Thacher.
!
!    PSI was written at Argonne National Laboratory for the FUNPACK
!    package of special function subroutines.  PSI was modified by
!    A. H. Morris (NSWC).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2010
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Anthony Strecok, Henry Thacher.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Anthony Strecok, Henry Thacher,
!    Chebyshev Approximations for the Psi Function,
!    Mathematics of Computation,
!    Volume 27, Number 121, January 1973, pages 123-127.
!
!  Parameters:
!
!    Input, real XX, the argument of the psi function.
!
!    Output, real PSI, the value of the psi function.  PSI is assigned
!    the value 0 when the psi function is undefined.
!
  implicit none

  real, parameter :: dx0 = 1.461632144968362341262659542325721325E+00
  real, parameter :: piov4 = 0.785398163397448E+00

  real aug
  real den
  real eps
  integer i
  integer m
  integer n
  integer nq
  real, parameter, dimension ( 7 ) :: p1 = (/&
    0.895385022981970E-02, &
    0.477762828042627E+01, &
    0.142441585084029E+03, &
    0.118645200713425E+04, &
    0.363351846806499E+04, &
    0.413810161269013E+04, &
    0.130560269827897E+04 /)
  real, parameter, dimension ( 4 ) :: p2 = (/ &
    -0.212940445131011E+01, &
    -0.701677227766759E+01, &
    -0.448616543918019E+01, &
    -0.648157123766197E+00 /)
  real psi
  real, parameter, dimension ( 6 ) :: q1 = (/&
    0.448452573429826E+02, &
    0.520752771467162E+03, &
    0.221000799247830E+04, &
    0.364127349079381E+04, &
    0.190831076596300E+04, &
    0.691091682714533E-05 /)
  real, parameter, dimension ( 4 ) :: q2 = (/ &
    0.322703493791143E+02, &
    0.892920700481861E+02, &
    0.546117738103215E+02, &
    0.777788548522962E+01 /)
  real sgn
  real temp
  real upper
  real w
  real x
  real xmax1
  real xmx0
  real xsmall
  real xx
  real z
!
!  XMAX1 is the largest positive floating point constant with entirely
!  integer representation.  It is also used as negative of lower bound
!  on acceptable negative arguments and as the positive argument beyond which
!  psi may be represented as LOG(X).
!
  xmax1 = dble ( 2147483647 )

  eps = epsilon ( eps )

  xmax1 = min ( xmax1, 1.0E+00 / eps )
!
!  XSMALL is the absolute argument below which PI*COTAN(PI*X)
!  may be represented by 1/X.
!
  xsmall = 1.0E-09

  x = xx
  aug = 0.0E+00

  if ( x == 0.0E+00 ) then
    psi = 0.0E+00
    return
  end if
!
!  X < 0.5,  Use reflection formula PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
!
  if ( x < 0.5E+00 ) then
!
!  0 < ABS(X) <= XSMALL.  Use 1/X as a substitute for PI*COTAN(PI*X)
!
    if ( abs ( x ) <= xsmall ) then
      aug = - 1.0E+00 / x
      go to 40
    end if
!
!  Reduction of argument for cotan
!
    w = - x
    sgn = piov4

    if ( w <= 0.0E+00 ) then
      w = - w
      sgn = -sgn
    end if
!
!  Make an error exit if X <= -XMAX1
!
    if ( xmax1 <= w ) then
      psi = 0.0E+00
      return
    end if

    nq = int ( w )
    w = w - real ( nq )
    nq = int ( w * 4.0E+00 )
    w = 4.0E+00 * ( w - real ( nq ) * 0.25E+00 )
!
!  W is now related to the fractional part of  4.0 * X.
!  Adjust argument to correspond to values in first
!  quadrant and determine sign
!
    n = nq / 2
    if ( ( n + n ) /= nq ) then
      w = 1.0E+00 - w
    end if

    z = piov4 * w
    m = n / 2

    if ( ( m + m ) /= n ) then
      sgn = - sgn
    end if
!
!  Determine final value for  -PI*COTAN(PI*X)
!
    n = ( nq + 1 ) / 2
    m = n / 2
    m = m + m

    if ( m == n ) then

      if ( z == 0.0E+00 ) then
        psi = 0.0E+00
       return
    end  if

      aug = 4.0E+00 * sgn * ( cos(z) / sin(z) )

    else

      aug = 4.0E+00 * sgn * ( sin(z) / cos(z) )

    end if

   40   continue

    x = 1.0E+00 - x

  end if
!
!  0.5 <= X <= 3.0
!
  if ( x <= 3.0E+00 ) then

    den = x
    upper = p1(1) * x

    do i = 1, 5
      den = ( den + q1(i) ) * x
      upper = ( upper + p1(i+1) ) * x
    end do

    den = ( upper + p1(7) ) / ( den + q1(6) )
    xmx0 = dble ( x ) - dx0
    psi = den * xmx0 + aug
!
!  3.0 < X < XMAX1
!
  else if ( x < xmax1 ) then

    w = 1.0E+00 / x**2
    den = w
    upper = p2(1) * w

    do i = 1, 3
      den = ( den + q2(i) ) * w
      upper = ( upper + p2(i+1) ) * w
    end do

    aug = upper / ( den + q2(4) ) - 0.5E+00 / x + aug
    psi = aug + log ( x )
!
!  XMAX1 <= X
!
  else

    psi = aug + log ( x )

  end if

  return
end
function r4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R4_UNIFORM returns a scaled pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 4 ) r4_uniform
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r4_uniform = a + ( b - a ) * real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
subroutine r4col_mean ( m, n, a, mean )

!*****************************************************************************80
!
!! R4COL_MEAN returns the means of columns of an R4COL.
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MEAN =
!      1.5  4.0  5.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer the leading dimension of the array, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real A(M,N), the array to be examined.
!
!    Output, real MEAN(N), the means, or averages, of the columns.
!
  implicit none

  integer m
  integer n

  real a(m,n)
  integer j
  real mean(n)

  do j = 1, n

    mean(j) = sum ( a(1:m,j) ) / real ( m )

  end do

  return
end
subroutine r4col_variance ( m, n, a, variance )

!*****************************************************************************80
!
!! R4COL_VARIANCE returns the variances of the columns of an R4COL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, real AM,N), the array whose variances are desired.
!
!    Output, real VARIANCE(N), the variances of the rows.
!
  implicit none

  integer m
  integer n

  real a(m,n)
  integer i
  integer j
  real mean
  real variance(n)

  do j = 1, n

    mean = sum ( a(1:m,j) ) / real ( m)

    variance(j) = 0.0E+00
    do i = 1, m
      variance(j) = variance(j) + ( a(i,j) - mean )**2
    end do

    if ( 1 < m ) then
      variance(j) = variance(j) / real ( m - 1 )
    else
      variance(j) = 0.0E+00
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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
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
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
function trigamma ( x )

!*****************************************************************************80
!
!! TRIGAMMA calculates trigamma(x) = d**2 log(Gamma(x)) / dx**2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    Original FORTRAN77 version by B Schneider.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    B Schneider,
!    Algorithm AS 121:
!    Trigamma Function,
!    Applied Statistics,
!    Volume 27, Number 1, 1978, page 97-99.
!
!  Parameters:
!
!    Input, real X, the argument of the trigamma function.
!    0 < X.
!
!    Output, real TRIGAMMA, the value of the trigamma function at X.
!
  implicit none

  real, parameter :: a = 0.0001E+00
  real, parameter :: b = 5.0E+00
  real, parameter :: b2 =   1.0E+00 / 6.0E+00
  real, parameter :: b4 = - 1.0E+00 / 30.0E+00
  real, parameter :: b6 =   1.0E+00 / 42.0E+00
  real, parameter :: b8 = - 1.0E+00 / 30.0E+00
  real trigamma
  real x
  real y
  real z
!
!  1): If X is not positive, fail.
!
  if ( x <= 0.0E+00 ) then

    trigamma = 0.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIGAMNA - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
!
!  2): If X is smaller than A, use a small value approximation.
!
  else if ( x <= a ) then

    trigamma = 1.0E+00 / x / x
!
!  3): Otherwise, increase the argument to B <+ ( X + I ).
!
  else

    z = x
    trigamma = 0.0E+00

    do while ( z < b )
      trigamma = trigamma + 1.0E+00 / z**2
      z = z + 1.0E+00
    end do
!
!  ...and then apply an asymptotic formula.
!
    y = 1.0E+00 / z / z

    trigamma = trigamma + 0.5E+00 * &
           y + ( 1.0E+00 &
         + y * ( b2 &
         + y * ( b4 &
         + y * ( b6 &
          + y *  b8 )))) / z

  end if

  return
end
