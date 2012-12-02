subroutine laguerre_compute ( order, xtab, weight, alpha )

!*****************************************************************************80
!
!! LAGUERRE_COMPUTE computes a Gauss-Laguerre quadrature rule.
!
!  Discussion:
!
!    In the simplest case, ALPHA is 0, and we are approximating the
!    integral from 0 to +oo of EXP(-X) * F(X).  When this is so,
!    it is easy to modify the rule to approximate the integral from
!    A to +oo as well.
!
!    If ALPHA is nonzero, then there is no simple way to extend the
!    rule to approximate the integral from A to +oo.  The simplest
!    procedures would be to approximate the integral from 0 to A.
!
!    The integration interval is [ A, +oo ) or [ 0, +oo ).
!
!    The weight function is w(x) = exp ( -x ) or exp ( -x ) * x^alpha.
!
!
!    If the integral to approximate is:
!
!        Integral ( A <= X < +oo ) exp ( - X ) * F(X) dX
!      or
!        Integral ( 0 <= X < +oo ) exp ( - X ) * X^ALPHA * F(X) dX
!
!    then the quadrature rule is:
!
!      exp ( - A ) * Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( A+XTAB(I) )
!    or
!      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!
!    If the integral to approximate is:
!
!        Integral ( A <= X < +oo ) F(X) dX
!      or
!        Integral ( 0 <= X < +oo ) X^ALPHA * F(X) dX
!
!    then the quadrature rule is:
!
!      exp ( - A ) * Sum ( 1 <= I <= ORDER ) 
!        WEIGHT(I) * EXP(A+XTAB(I)) * F ( A+XTAB(I) )
!    or
!      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * EXP(XTAB(I)) * F ( XTAB(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the quadrature rule 
!    to be computed.  ORDER must be at least 1.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    Set ALPHA = 0.0D+00 for the simplest rule.
!    ALPHA must be nonnegative.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(order)
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) cc
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) ratio
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(order)
!
!  Set the recursion coefficients.
!
  do i = 1, order
    b(i) = ( alpha + real ( 2 * i - 1, kind = 8 ) )
  end do

  do i = 1, order
    c(i) = real ( i - 1, kind = 8 ) * ( alpha + real ( i - 1, kind = 8 ) )
  end do

  cc = r8_gamma ( alpha + 1.0D+00 ) * product ( c(2:order) )

  do i = 1, order
!
!  Compute an estimate for the root.
!
    if ( i == 1 ) then

      x = ( 1.0D+00 + alpha ) * ( 3.0D+00+ 0.92 * alpha ) / &
        ( 1.0D+00 + 2.4D+00 * real ( order, kind = 8 ) + 1.8D+00 * alpha )

    else if ( i == 2 ) then

      x = x + ( 15.0D+00 + 6.25D+00 * alpha ) / &
        ( 1.0D+00 + 0.9D+00 * alpha + 2.5D+00 * real ( order, kind = 8 ) )

    else

      r1 = ( 1.0D+00 + 2.55D+00 * real ( i - 2, kind = 8 ) ) &
        / ( 1.9D+00 * real ( i - 2, kind = 8 ) )

      r2 = 1.26D+00 * real ( i - 2, kind = 8 ) * alpha / &
        ( 1.0D+00 + 3.5D+00 * real ( i - 2, kind = 8 ) )

      ratio = ( r1 + r2 ) / ( 1.0D+00 + 0.3D+00 * alpha )

      x = x + ratio * ( x - xtab(i-2) )

    end if
!
!  Use iteration to find the root.
!
    call laguerre_root ( x, order, alpha, dp2, p1, b, c )
!
!  Set the abscissa and weight.
!
    xtab(i) = x
    weight(i) = ( cc / dp2 ) / p1

  end do

  return
end
subroutine laguerre_recur ( p2, dp2, p1, x, order, alpha, b, c )

!*****************************************************************************80
!
!! LAGUERRE_RECUR finds the value and derivative of a Laguerre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 1998
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of L(ORDER)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of L'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(ORDER-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial 
!    to be computed.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor in the
!    integrand.
!
!    Input, real ( kind = 8 ) B(ORDER), C(ORDER), the recursion
!    coefficients.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(order)
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x - alpha - 1.0D+00
  dp2 = 1.0D+00

  do i = 2, order

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = ( x - b(i) ) * p1 - c(i) * p0
    dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine laguerre_root ( x, order, alpha, dp2, p1, b, c )

!*****************************************************************************80
!
!! LAGUERRE_ROOT improves an approximate root of a Laguerre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2000
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial 
!    to be computed.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!
!    Output, real ( kind = 8 ) DP2, the value of L'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(ORDER-1)(X).
!
!    Input, real ( kind = 8 ) B(ORDER), C(ORDER), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(order)
  real ( kind = 8 ) c(order)
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = epsilon ( x )

  do step = 1, step_max

    call laguerre_recur ( p2, dp2, p1, x, order, alpha, b, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine legendre_compute ( order, xtab, weight )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE computes a Gauss-Legendre quadrature rule.
!
!  Discussion:
!
!    The integration interval is [ -1, 1 ].
!
!    The weight function is w(x) = 1.0.
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    ORDER must be greater than 0.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!    The weights are positive, symmetric, and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2pn
  real ( kind = 8 ) d3pn
  real ( kind = 8 ) d4pn
  real ( kind = 8 ) dp
  real ( kind = 8 ) dpn
  real ( kind = 8 ) e1
  real ( kind = 8 ) fx
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iback
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mp1mi
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) nmove
  real ( kind = 8 ) p
  real ( kind = 8 ) :: pi = 3.1415926535897932385D+00
  real ( kind = 8 ) pk
  real ( kind = 8 ) pkm1
  real ( kind = 8 ) pkp1
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x0
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) xtemp
  real ( kind = 8 ) weight(order)

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
  end if

  e1 = real ( order * ( order + 1 ), kind = 8 )

  m = ( order + 1 ) / 2

  do i = 1, m

    mp1mi = m + 1 - i

    t = real ( 4 * i - 1, kind = 8 ) * pi &
      / real ( 4 * order + 2, kind = 8 )

    x0 = cos ( t ) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 &
      / real ( order, kind = 8 ) ) &
      / real ( 8 * order * order, kind = 8 ) )

    pkm1 = 1.0D+00
    pk = x0

    do k = 2, order
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = 8 )
      pkm1 = pk
      pk = pkp1
    end do

    d1 = real ( order, kind = 8 ) * ( pkm1 - x0 * pk )

    dpn = d1 / ( 1.0D+00 - x0 * x0 )

    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 * x0 )

    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) &
      / ( 1.0D+00 - x0 * x0 )

    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) / &
      ( 1.0D+00 - x0 * x0 )

    u = pk / dpn
    v = d2pn / dpn
!
!  Initial approximation H:
!
    h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v * v - d3pn / &
      ( 3.0D+00 * dpn ) ) ) )
!
!  Refine H using one step of Newton's method:
!
    p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 &
      * ( d3pn + 0.25D+00 * h * d4pn ) ) )

    dp = dpn + h * ( d2pn + 0.5D+00 * h * ( d3pn + h * d4pn / 3.0D+00 ) )

    h = h - p / dp

    xtemp = x0 + h

    xtab(mp1mi) = xtemp

    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )

    weight(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp * xtemp ) / ( fx * fx )

  end do

  if ( mod ( order, 2 ) == 1 ) then
    xtab(1) = 0.0D+00
  end if
!
!  Shift the data up.
!
  nmove = ( order + 1 ) / 2
  ncopy = order - nmove

  do i = 1, nmove
    iback = order + 1 - i
    xtab(iback) = xtab(iback-ncopy)
    weight(iback) = weight(iback-ncopy)
  end do
!
!  Reflect values for the negative abscissas.
!
  do i = 1, order - nmove
    xtab(i) = - xtab(order+1-i)
    weight(i) = weight(order+1-i)
  end do

  return
end
subroutine p00_alpha ( problem, alpha )

!*****************************************************************************80
!
!! P00_ALPHA returns the value of ALPHA for any problem.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!    The typical or default value is 0.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the problem.
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) problem

  if ( problem == 1 ) then
    call p01_alpha ( alpha )
  else if ( problem == 2 ) then
    call p02_alpha ( alpha )
  else if ( problem == 3 ) then
    call p03_alpha ( alpha )
  else if ( problem == 4 ) then
    call p04_alpha ( alpha )
  else if ( problem == 5 ) then
    call p05_alpha ( alpha )
  else if ( problem == 6 ) then
    call p06_alpha ( alpha )
  else if ( problem == 7 ) then
    call p07_alpha ( alpha )
  else if ( problem == 8 ) then
    call p08_alpha ( alpha )
  else if ( problem == 9 ) then
    call p09_alpha ( alpha )
  else if ( problem == 10 ) then
    call p10_alpha ( alpha )
  else if ( problem == 11 ) then
    call p11_alpha ( alpha )
  else if ( problem == 12 ) then
    call p12_alpha ( alpha )
  else if ( problem == 13 ) then
    call p13_alpha ( alpha )
  else if ( problem == 14 ) then
    call p14_alpha ( alpha )
  else if ( problem == 15 ) then
    call p15_alpha ( alpha )
  else if ( problem == 16 ) then
    call p16_alpha ( alpha )
  else if ( problem == 17 ) then
    call p17_alpha ( alpha )
  else if ( problem == 18 ) then
    call p18_alpha ( alpha )
  else if ( problem == 19 ) then
    call p19_alpha ( alpha )
  else if ( problem == 20 ) then
    call p20_alpha ( alpha )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_ALPHA - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_exact ( problem, exact )

!*****************************************************************************80
!
!! P00_EXACT returns the exact integral for any problem.
!
!  Discussion:
!
!    This routine provides a "generic" interface to the exact integral
!    routines for the various problems, and allows a problem to be called
!    by index (PROBLEM) rather than by name.
!
!    In most cases, the "exact" value of the integral is not given;
!    instead a "respectable" approximation is available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the problem.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) problem

  if ( problem == 1 ) then
    call p01_exact ( exact )
  else if ( problem == 2 ) then
    call p02_exact ( exact )
  else if ( problem == 3 ) then
    call p03_exact ( exact )
  else if ( problem == 4 ) then
    call p04_exact ( exact )
  else if ( problem == 5 ) then
    call p05_exact ( exact )
  else if ( problem == 6 ) then
    call p06_exact ( exact )
  else if ( problem == 7 ) then
    call p07_exact ( exact )
  else if ( problem == 8 ) then
    call p08_exact ( exact )
  else if ( problem == 9 ) then
    call p09_exact ( exact )
  else if ( problem == 10 ) then
    call p10_exact ( exact )
  else if ( problem == 11 ) then
    call p11_exact ( exact )
  else if ( problem == 12 ) then
    call p12_exact ( exact )
  else if ( problem == 13 ) then
    call p13_exact ( exact )
  else if ( problem == 14 ) then
    call p14_exact ( exact )
  else if ( problem == 15 ) then
    call p15_exact ( exact )
  else if ( problem == 16 ) then
    call p16_exact ( exact )
  else if ( problem == 17 ) then
    call p17_exact ( exact )
  else if ( problem == 18 ) then
    call p18_exact ( exact )
  else if ( problem == 19 ) then
    call p19_exact ( exact )
  else if ( problem == 20 ) then
    call p20_exact ( exact )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_EXACT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_exp_transform ( problem, order, result )

!*****************************************************************************80
!
!! P00_EXP_TRANSFORM applies an exponential transform and Gauss-Legendre rule.
!
!  Discussion:
!
!    To approximate:
!
!      Integral ( alpha <= x < +oo ) f(x) dx
!
!    Transform:
!
!      u = exp ( -x )
!      du = - exp ( -x ) dx
!
!      x = - log ( u )
!      dx = - du / u
!
!      x = alpha    => u = exp ( -alpha )
!      x = +oo      => u = 0
!
!    Transformed integral:
!
!      Integral ( 0 < u <= exp ( -alpha ) ) f ( -log(u) ) du / u
!
!    We apply a Gauss-Legendre rule here, but we could easily use any rule
!    that avoids evaluation at U = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the problem.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the Gauss-Legendre rule 
!    to apply.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ), allocatable, dimension ( : ) :: f_vec
  integer ( kind = 4 ) order
  integer ( kind = 4 ) problem
  real ( kind = 8 ) result
  real ( kind = 8 ), allocatable, dimension ( : ) :: u
  real ( kind = 8 ), allocatable, dimension ( : ) :: u_log
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight

  allocate ( f_vec(1:order) )
  allocate ( u(1:order) )
  allocate ( u_log(1:order) )
  allocate ( weight(1:order) )

  call p00_alpha ( problem, alpha )
!
!  Get the abscissas and weights for Gauss-Legendre quadrature.
!
  call legendre_compute ( order, u, weight )
!
!  Modify the weights from [-1,1] to [0,exp(-alpha)].
!
  weight(1:order) = exp ( -alpha ) * weight(1:order) / 2.0D+00
!
!  Linear transform of abscissas from [-1,1] to [0,exp(-alpha)].
!
  u(1:order) = ( ( 1.0D+00 + u(1:order) ) * exp ( - alpha ) &
               + ( 1.0D+00 - u(1:order) ) * 0.0D+00 )       &
               / ( 2.0D+00              )
!
!  Define U_LOG = - log ( U )
!
  u_log(1:order) = - log ( u(1:order) )
!
!  Evaluate F ( -LOG(U) ).
!
  call p00_fun ( problem, order, u_log, f_vec )
!
!  The integrand is F ( -LOG(U) ) / U
!
  f_vec(1:order) = f_vec(1:order) / u(1:order)
!
!  Sum.
!
  result = dot_product ( weight(1:order), f_vec(1:order) )

  deallocate ( f_vec )
  deallocate ( u )
  deallocate ( u_log )
  deallocate ( weight )

  return
end
subroutine p00_fun ( problem, n, x, f )

!*****************************************************************************80
!
!! P00_FUN evaluates the integrand for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the problem.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)

  if ( problem == 1 ) then
    call p01_fun ( n, x, f )
  else if ( problem == 2 ) then
    call p02_fun ( n, x, f )
  else if ( problem == 3 ) then
    call p03_fun ( n, x, f )
  else if ( problem == 4 ) then
    call p04_fun ( n, x, f )
  else if ( problem == 5 ) then
    call p05_fun ( n, x, f )
  else if ( problem == 6 ) then
    call p06_fun ( n, x, f )
  else if ( problem == 7 ) then
    call p07_fun ( n, x, f )
  else if ( problem == 8 ) then
    call p08_fun ( n, x, f )
  else if ( problem == 9 ) then
    call p09_fun ( n, x, f )
  else if ( problem == 10 ) then
    call p10_fun ( n, x, f )
  else if ( problem == 11 ) then
    call p11_fun ( n, x, f )
  else if ( problem == 12 ) then
    call p12_fun ( n, x, f )
  else if ( problem == 13 ) then
    call p13_fun ( n, x, f )
  else if ( problem == 14 ) then
    call p14_fun ( n, x, f )
  else if ( problem == 15 ) then
    call p15_fun ( n, x, f )
  else if ( problem == 16 ) then
    call p16_fun ( n, x, f )
  else if ( problem == 17 ) then
    call p17_fun ( n, x, f )
  else if ( problem == 18 ) then
    call p18_fun ( n, x, f )
  else if ( problem == 19 ) then
    call p19_fun ( n, x, f )
  else if ( problem == 20 ) then
    call p20_fun ( n, x, f )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FUN - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_gauss_laguerre ( problem, order, result )

!*****************************************************************************80
!
!! P00_GAUSS_LAGUERRE applies a Gauss-Laguerre rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the problem.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the Gauss-Laguerre rule 
!    to apply.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha2
  real ( kind = 8 ), allocatable, dimension ( : ) :: f_vec
  integer ( kind = 4 ) order
  integer ( kind = 4 ) problem
  real ( kind = 8 ) result
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab

  allocate ( f_vec(1:order) )
  allocate ( weight(1:order) )
  allocate ( xtab(1:order) )

  call p00_alpha ( problem, alpha )

  alpha2 = 0.0D+00
  call laguerre_compute ( order, xtab, weight, alpha2 )

  xtab(1:order) = xtab(1:order) + alpha

  call p00_fun ( problem, order, xtab, f_vec )
!
!  The Gauss-Laguerre rule assumes a weight of EXP(-X).
!
!  We need to multiply each F(X) by EXP(X) to implicitly 
!  adjust for this weight.
!
  f_vec(1:order) = f_vec(1:order) * exp ( xtab(1:order) )

  result = exp ( -alpha ) * dot_product ( weight(1:order), f_vec(1:order) )

  deallocate ( f_vec )
  deallocate ( weight )
  deallocate ( xtab )

  return
end
subroutine p00_problem_num ( problem_num )

!*****************************************************************************80
!
!! P00_PROBLEM_NUM returns the number of test integration problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) PROBLEM_NUM, the number of test problems.
!
  implicit none

  integer ( kind = 4 ) problem_num

  problem_num = 20

  return
end
subroutine p00_rat_transform ( problem, order, result )

!*****************************************************************************80
!
!! P00_RAT_TRANSFORM applies a rational transform and Gauss-Legendre rule.
!
!  Discussion:
!
!    To approximate:
!
!      Integral ( alpha <= x < +oo ) f(x) dx
!
!    Transform:
!
!      u = 1 / ( 1 + x )
!      du = - dx / ( 1 + x )^2
!
!      x = ( 1 - u ) / u
!      dx = - du / u^2
!
!      x = alpha    => u = 1 / ( 1 + alpha )
!      x = +oo      => u = 0
!
!    Transformed integral:
!
!      Integral ( 0 < u <= 1 / ( 1 + alpha ) ) f ( ( 1 - u ) / u ) du / u^2
!
!    We apply a Gauss-Legendre rule here, but we could easily use any rule
!    that avoids evaluation at U = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the problem.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the Gauss-Legendre rule 
!    to apply.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ), allocatable, dimension ( : ) :: f_vec
  integer ( kind = 4 ) order
  integer ( kind = 4 ) problem
  real ( kind = 8 ) result
  real ( kind = 8 ), allocatable, dimension ( : ) :: u
  real ( kind = 8 ), allocatable, dimension ( : ) :: u_rat
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight

  allocate ( f_vec(1:order) )
  allocate ( u(1:order) )
  allocate ( u_rat(1:order) )
  allocate ( weight(1:order) )

  call p00_alpha ( problem, alpha )
!
!  Get the abscissas and weights for Gauss-Legendre quadrature.
!
  call legendre_compute ( order, u, weight )
!
!  Modify the weights from [-1,1] to [0,1/(1+alpha)].
!
  weight(1:order) = weight(1:order) / 2.0D+00 / ( 1.0D+00 + alpha )
!
!  Linear transform of abscissas from [-1,1] to [0,1/(1+alpha)].
!
  u(1:order) = ( ( 1.0D+00 + u(1:order) ) / ( 1.0D+00 + alpha ) &
               + ( 1.0D+00 - u(1:order) ) * 0.0D+00 )       &
               / ( 2.0D+00              )
!
!  Define U_RAT = ( 1 - U ) / U.
!
  u_rat(1:order) = ( 1.0D+00 - u(1:order) ) / u(1:order)
!
!  Evaluate F ( ( 1 - U ) / U ).
!
  call p00_fun ( problem, order, u_rat, f_vec )
!
!  The integrand is F ( ( 1 - U ) / U ) / U^2
!
  f_vec(1:order) = f_vec(1:order) / u(1:order)**2
!
!  Sum.
!
  result = dot_product ( weight(1:order), f_vec(1:order) )

  deallocate ( f_vec )
  deallocate ( u )
  deallocate ( u_rat )
  deallocate ( weight )

  return
end
subroutine p00_title ( problem, title )

!*****************************************************************************80
!
!! P00_TITLE returns the title for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the problem.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer ( kind = 4 ) problem
  character ( len = * ) title

  if ( problem == 1 ) then
    call p01_title ( title )
  else if ( problem == 2 ) then
    call p02_title ( title )
  else if ( problem == 3 ) then
    call p03_title ( title )
  else if ( problem == 4 ) then
    call p04_title ( title )
  else if ( problem == 5 ) then
    call p05_title ( title )
  else if ( problem == 6 ) then
    call p06_title ( title )
  else if ( problem == 7 ) then
    call p07_title ( title )
  else if ( problem == 8 ) then
    call p08_title ( title )
  else if ( problem == 9 ) then
    call p09_title ( title )
  else if ( problem == 10 ) then
    call p10_title ( title )
  else if ( problem == 11 ) then
    call p11_title ( title )
  else if ( problem == 12 ) then
    call p12_title ( title )
  else if ( problem == 13 ) then
    call p13_title ( title )
  else if ( problem == 14 ) then
    call p14_title ( title )
  else if ( problem == 15 ) then
    call p15_title ( title )
  else if ( problem == 16 ) then
    call p16_title ( title )
  else if ( problem == 17 ) then
    call p17_title ( title )
  else if ( problem == 18 ) then
    call p18_title ( title )
  else if ( problem == 19 ) then
    call p19_title ( title )
  else if ( problem == 20 ) then
    call p20_title ( title )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p01_alpha ( alpha )

!*****************************************************************************80
!
!! P01_ALPHA returns ALPHA for problem 1.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 2.0D+00

  return
end
subroutine p01_exact ( exact )

!*****************************************************************************80
!
!! P01_EXACT returns the exact integral for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.19524754198276439152D+00

  return
end
subroutine p01_fun ( n, x, f )

!*****************************************************************************80
!
!! P01_FUN evaluates the integrand for problem 1.
!
!  Discussion:
!
!    D&R gives "exact" value as 0.19524753...
!    Mathematica returns        0.19524754198276439152...
!    D&R gives Laguerre(16) as  0.16623627...
!
!  Integral
!
!    exp ( -2 ) Integral ( 2 <= x < +oo ) dx / ( x * log(x)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    f(i) = exp ( -2.0D+00 ) / ( x(i) * log ( x(i) )**2 )
  end do

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns the title for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = '1 / ( x * log(x)^2 )'

  return
end
subroutine p02_alpha ( alpha )

!*****************************************************************************80
!
!! P02_ALPHA returns ALPHA for problem 2.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 2.0D+00

  return
end
subroutine p02_exact ( exact )

!*****************************************************************************80
!
!! P02_EXACT returns the exact integral for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.32510848278991335198D+00

  return
end
subroutine p02_fun ( n, x, f )

!*****************************************************************************80
!
!! P02_FUN evaluates the integrand for problem 2.
!
!  Discussion:
!
!    D&R gives "exact" value as 0.32510855...
!    Mathematica returns        0.32510848278991335198...
!    D&R gives Laguerre(16) as  0.19142399...
!
!  Integral
!
!    exp ( -2 ) Integral ( 2 <= x < +oo ) dx / ( x * log(x)^(3/2) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    f(i) = exp ( -2.0D+00 ) / ( x(i) * sqrt ( log ( x(i) )**3 ) )
  end do

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns the title for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = '1 / ( x * log(x)^(3/2) )'

  return
end
subroutine p03_alpha ( alpha )

!*****************************************************************************80
!
!! P03_ALPHA returns ALPHA for problem 3.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 2.0D+00

  return
end
subroutine p03_exact ( exact )

!*****************************************************************************80
!
!! P03_EXACT returns the exact integral for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 13.628D+00

  return
end
subroutine p03_fun ( n, x, f )

!*****************************************************************************80
!
!! P03_FUN evaluates the integrand for problem 3.
!
!  Discussion:
!
!    D&R gives "exact" value as 13.628...
!    Mathematica returns        13.440045415012575106...
!    D&R gives Laguerre(16) as   0.44996932...
!
!    This integral is "something of a numerical joke, as it is
!    scarcely distinguishable from the divergent integrand 1/x."
!
!  Integral
!
!    exp ( -2 ) Integral ( 2 <= x < +oo ) dx / ( x^1.01 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    f(i) = exp ( -2.0D+00 ) * 1.0D+00 / x(i)**1.01D+00
  end do

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns the title for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = '1 / ( x^1.01 )'

  return
end
subroutine p04_alpha ( alpha )

!*****************************************************************************80
!
!! P04_ALPHA returns ALPHA for problem 4.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 2.0D+00

  return
end
subroutine p04_exact ( exact )

!*****************************************************************************80
!
!! P04_EXACT returns the estimated integral for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = -0.0046848541335080643181D+00

  return
end
subroutine p04_fun ( n, x, f )

!*****************************************************************************80
!
!! P04_FUN evaluates the integrand for problem 4.
!
!  Discussion:
!
!    D&R gives "exact" value as -0.0046984...
!    Mathematica returns        -0.0046848541335080643181...
!    D&R gives Laguerre(16) as  -0.039258696...
!
!  Integral
!
!    exp ( -2 ) Integral ( 2 <= x < +oo ) ( sin ( x ) / x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( x(i) == 0.0D+00 ) then
      f(i) = exp ( -2.0D+00 )
    else
      f(i) = exp ( -2.0D+00 ) * sin ( x(i) ) / x(i)
    end if
  end do

  return
end
subroutine p04_title ( title )

!*****************************************************************************80
!
!! P04_TITLE returns the title for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Sine integral'

  return
end
subroutine p05_alpha ( alpha )

!*****************************************************************************80
!
!! P05_ALPHA returns ALPHA for problem 5.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 2.0D+00

  return
end
subroutine p05_exact ( exact )

!*****************************************************************************80
!
!! P05_EXACT returns the estimated integral for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.0015897286158592328774D+00

  return
end
subroutine p05_fun ( n, x, f )

!*****************************************************************************80
!
!! P05_FUN evaluates the integrand for problem 5.
!
!  Discussion:
!
!    D&R gives "exact" value as  0.00158973...
!    Mathematica returns         0.0015897286158592328774...
!    D&R gives Laguerre(16) as  -0.067859545...
!
!  Integral
!
!    exp ( -2 ) Integral ( 2 <= x < +oo ) cos ( pi * x^2 / 2 ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932385D+00
  real ( kind = 8 ) x(n)

  do i = 1, n
    f(i) = exp ( -2.0D+00 ) * cos ( 0.5D+00 * pi * x(i)**2 )
  end do

  return
end
subroutine p05_title ( title )

!*****************************************************************************80
!
!! P05_TITLE returns the title for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Fresnel integral'

  return
end
subroutine p06_alpha ( alpha )

!*****************************************************************************80
!
!! P06_ALPHA returns ALPHA for problem 6.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 2.0D+00

  return
end
subroutine p06_exact ( exact )

!*****************************************************************************80
!
!! P06_EXACT returns the exact integral for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.00056103711148387120640D+00

  return
end
subroutine p06_fun ( n, x, f )

!*****************************************************************************80
!
!! P06_FUN evaluates the integrand for problem 6.
!
!  Discussion:
!
!    D&R gives "exact" value as 0.0005610371...
!    Mathematica returns        0.00056103711148387120640...
!    D&R gives Laguerre(16) as  0.00056100775...
!
!  Integral
!
!    exp ( -2 ) Integral ( 2 <= x < +oo ) exp ( -x^2 ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: exponent_min = -80.0D+00
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( - x(i)**2 < exponent_min ) then
      f(i) = 0.0D+00
    else
      f(i) = exp ( -2.0D+00 ) * exp ( - x(i)**2 )
    end if
  end do

  return
end
subroutine p06_title ( title )

!*****************************************************************************80
!
!! P06_TITLE returns the title for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Complementary error function'

  return
end
subroutine p07_alpha ( alpha )

!*****************************************************************************80
!
!! P07_ALPHA returns ALPHA for problem 7.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 2.0D+00

  return
end
subroutine p07_exact ( exact )

!*****************************************************************************80
!
!! P07_EXACT returns the exact integral for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.16266891D+00

  return
end
subroutine p07_fun ( n, x, f )

!*****************************************************************************80
!
!! P07_FUN evaluates the integrand for problem 7.
!
!  Discussion:
!
!    D&R gives "exact" value as 0.16266891...
!    Mathematica does not return a value.
!    D&R gives Laguerre(16) as  0.097083064...
!
!  Integral
!
!    exp ( -2 ) Integral ( 2 <= x < +oo ) sin ( x - 1 ) dx
!      / sqrt ( x * ( x - 2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) == 2.0D+00 ) then
      f(i) = 0.0D+00
    else
      f(i) = exp ( -2.0D+00 ) &
        * sin ( x(i) - 1.0D+00 ) / sqrt ( x(i) * ( x(i) - 2.0D+00 ) )
    end if

  end do

  return
end
subroutine p07_title ( title )

!*****************************************************************************80
!
!! P07_TITLE returns the title for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Bessel function'

  return
end
subroutine p08_alpha ( alpha )

!*****************************************************************************80
!
!! P08_ALPHA returns ALPHA for problem 8.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 2.0D+00

  return
end
subroutine p08_exact ( exact )

!*****************************************************************************80
!
!! P08_EXACT returns the estimated integral for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.058334852497734677320D+00

  return
end
subroutine p08_fun ( n, x, f )

!*****************************************************************************80
!
!! P08_FUN evaluates the integrand for problem 8.
!
!  Discussion:
!
!    D&R gives "exact" value as 0.0583349...
!    Mathematica returns        0.058334852497734677320...
!    D&R gives Laguerre(16) as  0.05834841...
!
!  Integral
!
!    exp ( -2 ) Integral ( 2 <= x < +oo ) x / ( exp ( x ) - 1 ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: exponent_max = 80.0D+00
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( x(i) < exponent_max ) then
      f(i) = exp ( -2.0D+00 ) * x(i) / ( exp ( x(i) ) - 1.0D+00 )
    else
      f(i) = 0.0D+00
    end if
  end do

  return
end
subroutine p08_title ( title )

!*****************************************************************************80
!
!! P08_TITLE returns the title for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Debye function'

  return
end
subroutine p09_alpha ( alpha )

!*****************************************************************************80
!
!! P09_ALPHA returns ALPHA for problem 9.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.0D+00

  return
end
subroutine p09_exact ( exact )

!*****************************************************************************80
!
!! P09_EXACT returns the estimated integral for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 24.0D+00

  return
end
subroutine p09_fun ( n, x, f )

!*****************************************************************************80
!
!! P09_FUN evaluates the integrand for problem 9.
!
!  Discussion:
!
!    The integral is the definition of the Gamma function for
!    Z = 5, with exact value (Z-1)! = 24.
!
!  Integral
!
!    Integral ( 0 <= x < +oo ) x^4 exp ( -x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: exponent_min = -80.0D+00
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( -x(i) < exponent_min ) then
      f(i) = 0.0D+00
    else
      f(i) = x(i)**4 * exp ( -x(i) )
    end if
  end do

  return
end
subroutine p09_title ( title )

!*****************************************************************************80
!
!! P09_TITLE returns the title for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Gamma(Z=5) function'

  return
end
subroutine p10_alpha ( alpha )

!*****************************************************************************80
!
!! P10_ALPHA returns ALPHA for problem 10.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.0D+00

  return
end
subroutine p10_exact ( exact )

!*****************************************************************************80
!
!! P10_EXACT returns the estimated integral for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932385D+00

  exact = pi / 2.0D+00

  return
end
subroutine p10_fun ( n, x, f )

!*****************************************************************************80
!
!! P10_FUN evaluates the integrand for problem 10.
!
!  Discussion:
!
!    S&S gives exact value as pi/2 = 1.5707963267948966192...
!    S&S gives Laguerre(16) as       1.5537377347...
!    S&S gives EXP_TRANSFORM(16) as  1.4293043007...
!    S&S gives RAT_TRANSFORM(16) as  1.5707963267...
!
!  Integral
!
!    Integral ( 0 <= x < +oo ) 1/(1+x*x) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    f(i) = 1.0D+00 / ( 1.0D+00 + x(i) * x(i) )
  end do

  return
end
subroutine p10_title ( title )

!*****************************************************************************80
!
!! P10_TITLE returns the title for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = '1/(1+x*x)'

  return
end
subroutine p11_alpha ( alpha )

!*****************************************************************************80
!
!! P11_ALPHA returns ALPHA for problem 11.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.0D+00

  return
end
subroutine p11_exact ( exact )

!*****************************************************************************80
!
!! P11_EXACT returns the estimated integral for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932385D+00

  exact = pi

  return
end
subroutine p11_fun ( n, x, f )

!*****************************************************************************80
!
!! P11_FUN evaluates the integrand for problem 11.
!
!  Discussion:
!
!    S&S gives exact value as pi =  3.1415926535897932385...
!    S&S gives Laguerre(16) as      2.6652685196...
!    S&S gives EXP_TRANSFORM(16) as 2.3629036166...
!    S&S gives RAT_TRANSFORM(16) as 3.0360705907... 
!
!  Integral
!
!    Integral ( 0 <= x < +oo ) 1/((1+x)*sqrt(x)) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( x(i) == 0.0D+00 ) then
      f(i) = 0.0D+00
    else
      f(i) = 1.0D+00 / ( ( 1.0D+00 + x(i) ) * sqrt ( x(i) ) )
    end if
  end do

  return
end
subroutine p11_title ( title )

!*****************************************************************************80
!
!! P11_TITLE returns the title for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = '1 / ( (1+x) * sqrt(x) )'

  return
end
subroutine p12_alpha ( alpha )

!*****************************************************************************80
!
!! P12_ALPHA returns ALPHA for problem 12.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.0D+00

  return
end
subroutine p12_exact ( exact )

!*****************************************************************************80
!
!! P12_EXACT returns the estimated integral for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.5D+00

  return
end
subroutine p12_fun ( n, x, f )

!*****************************************************************************80
!
!! P12_FUN evaluates the integrand for problem 12.
!
!  Discussion:
!
!    S&S gives exact value as pi =  0.5
!    S&S gives Laguerre(16) as      0.5000000000...
!    S&S gives EXP_TRANSFORM(16) as 0.5019065783... 
!    S&S gives RAT_TRANSFORM(16) as 0.4988027685...
!
!  Integral
!
!    Integral ( 0 <= x < +oo ) exp ( -x ) * cos ( x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: exponent_min = -80.0D+00
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( -x(i) < exponent_min ) then
      f(i) = 0.0D+00
    else
      f(i) = exp ( -x(i) ) * cos ( x(i) )
    end if
  end do

  return
end
subroutine p12_title ( title )

!*****************************************************************************80
!
!! P12_TITLE returns the title for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'exp ( - x ) * cos ( x )'

  return
end
subroutine p13_alpha ( alpha )

!*****************************************************************************80
!
!! P13_ALPHA returns ALPHA for problem 13.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.0D+00

  return
end
subroutine p13_exact ( exact )

!*****************************************************************************80
!
!! P13_EXACT returns the estimated integral for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932385D+00

  exact = pi / 2.0D+00

  return
end
subroutine p13_fun ( n, x, f )

!*****************************************************************************80
!
!! P13_FUN evaluates the integrand for problem 13.
!
!  Discussion:
!
!    S&S gives exact value as pi/2 = 1.5707963267948966192...
!    S&S gives Laguerre(16) as       1.4399523793...
!    S&S gives EXP_TRANSFORM(16) as  1.3045186595...
!    S&S gives RAT_TRANSFORM(16) as  0.2046437026...
!
!  Integral
!
!    Integral ( 0 <= x < +oo ) sin ( x ) / x dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( x(i) == 0.0D+00 ) then
      f(i) = 1.0D+00
    else
      f(i) = sin ( x(i) ) / x(i)
    end if
  end do

  return
end
subroutine p13_title ( title )

!*****************************************************************************80
!
!! P13_TITLE returns the title for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'sin(x) / x'

  return
end
subroutine p14_alpha ( alpha )

!*****************************************************************************80
!
!! P14_ALPHA returns ALPHA for problem 14.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.0D+00

  return
end
subroutine p14_exact ( exact )

!*****************************************************************************80
!
!! P14_EXACT returns the estimated integral for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 1.0634618101722400407D+00

  return
end
subroutine p14_fun ( n, x, f )

!*****************************************************************************80
!
!! P14_FUN evaluates the integrand for problem 14.
!
!  Discussion:
!
!    S&S gives "exact" value as     1.0634618101...
!    Mathematica returns            1.0634618101722400407...
!    S&S gives Laguerre(16) as      1.0634713425...
!    S&S gives EXP_TRANSFORM(16) as 1.0634618101...
!    S&S gives RAT_TRANSFORM(16) as 1.0634574249...
!
!    The FORTRAN version of this routine, compiled with G95, was getting 
!    a floating point exception when evaluating the integrand
!    and using a Laguerre rule of order 64.  So I have had to truncate
!    the evaluation of the exponential.
!
!  Integral
!
!    Integral ( 0 <= x < +oo ) sin ( exp ( - x ) + exp ( - 4 x ) ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: exponent_min = -80.0D+00
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( - x(i) < exponent_min ) then
      f(i) = 0.0D+00
    else if ( -4.0D+00 * x(i) < exponent_min ) then
      f(i) = sin ( exp ( -x(i) ) )
    else
      f(i) = sin ( exp ( -x(i) ) + exp ( -4.0D+00 * x(i) ) )
    end if

  end do

  return
end
subroutine p14_title ( title )

!*****************************************************************************80
!
!! P14_TITLE returns the title for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'sin ( exp(-x) + exp(-4x) )'

  return
end
subroutine p15_alpha ( alpha )

!*****************************************************************************80
!
!! P15_ALPHA returns ALPHA for problem 15.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.0D+00

  return
end
subroutine p15_exact ( exact )

!*****************************************************************************80
!
!! P15_EXACT returns the estimated integral for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932385D+00

  exact = - pi * log ( 10.0D+00 ) / 20.0D+00

  return
end
subroutine p15_fun ( n, x, f )

!*****************************************************************************80
!
!! P15_FUN evaluates the integrand for problem 15.
!
!  Integral
!
!    Integral ( 0 <= x < +oo ) log(x) / (1+100*x*x) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise deDoncker-Kapenga, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983,
!    ISBN: 3540125531,
!    LC: QA299.3.Q36.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( x(i) == 0.0D+00 ) then
      f(i) = - huge ( x(i) )
    else
      f(i) = log ( x(i) ) / ( 1.0D+00 + 100.0D+00 * x(i) * x(i) )
    end if

  end do

  return
end
subroutine p15_title ( title )

!*****************************************************************************80
!
!! P15_TITLE returns the title for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'log(x) / ( 1 + 100 x^2 )'

  return
end
subroutine p16_alpha ( alpha )

!*****************************************************************************80
!
!! P16_ALPHA returns ALPHA for problem 16.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.0D+00

  return
end
subroutine p16_exact ( exact )

!*****************************************************************************80
!
!! P16_EXACT returns the estimated integral for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 1.0D+00

  return
end
subroutine p16_fun ( n, x, f )

!*****************************************************************************80
!
!! P16_FUN evaluates the integrand for problem 16.
!
!  Integral
!
!    Integral ( 0 <= x < +oo ) cos ( pi * x / 2 ) / sqrt ( x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise deDoncker-Kapenga, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983,
!    ISBN: 3540125531,
!    LC: QA299.3.Q36.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932385D+00
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( x(i) == 0.0D+00 ) then
      f(i) = huge ( x(i) )
    else
      f(i) = cos ( pi * x(i) / 2.0D+00 ) / sqrt ( x(i) )
    end if

  end do

  return
end
subroutine p16_title ( title )

!*****************************************************************************80
!
!! P16_TITLE returns the title for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'cos(pi x / 2 ) / sqrt(x)'

  return
end
subroutine p17_alpha ( alpha )

!*****************************************************************************80
!
!! P17_ALPHA returns ALPHA for problem 17.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.0D+00

  return
end
subroutine p17_exact ( exact )

!*****************************************************************************80
!
!! P17_EXACT returns the exact integral for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ), parameter :: beta = 2.0D+00
  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = sqrt ( pi ) * cos ( 0.5D+00 * atan ( 2.0D+00**beta ) ) &
    / sqrt ( sqrt ( 1.0D+00 + 0.25D+00**beta) )

  return
end
subroutine p17_fun ( n, x, fx )

!*****************************************************************************80
!
!! P17_FUN evaluates the integrand for problem 17.
!
!  Integral:
!
!    Integral ( 0 <= x < +oo ) exp ( - x / 2^beta ) * cos ( x ) / sqrt ( x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 84.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: beta = 2.0D+00
  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = exp ( - x(i) / 2.0D+00**beta ) * cos ( x(i) ) / sqrt ( x(i) )
    end if

  end do

  return
end
subroutine p17_title ( title )

!*****************************************************************************80
!
!! P17_TITLE returns the title for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'exp ( - x / 2^beta ) * cos ( x ) / sqrt ( x )'

  return
end
subroutine p18_alpha ( alpha )

!*****************************************************************************80
!
!! P18_ALPHA returns ALPHA for problem 18.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.0D+00

  return
end
subroutine p18_exact ( exact )

!*****************************************************************************80
!
!! P18_EXACT returns the exact integral for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ), parameter :: beta = 1.0D+00
  real ( kind = 8 ) exact

  exact = 2.0D+00**( 3.0D+00 * beta + 1.0D+00 )

  return
end
subroutine p18_fun ( n, x, fx )

!*****************************************************************************80
!
!! P18_FUN evaluates the integrand for problem 18.
!
!  Integral:
!
!    Integral ( 0 <= x < +oo ) x^2 * exp ( - x / 2^beta ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 84.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: beta = 1.0D+00
  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) x(n)

  fx(1:n) = x(1:n)**2 * exp ( - x(1:n) / 2**beta )

  return
end
subroutine p18_title ( title )

!*****************************************************************************80
!
!! P18_TITLE returns the title for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'x^2 * exp ( - x / 2^beta )'

  return
end
subroutine p19_alpha ( alpha )

!*****************************************************************************80
!
!! P19_ALPHA returns ALPHA for problem 19.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.0D+00

  return
end
subroutine p19_exact ( exact )

!*****************************************************************************80
!
!! P19_EXACT returns the exact integral for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ), parameter :: beta = 0.5D+00
  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  if ( beta == 1.0D+00 ) then
    exact = 1.0D+00 / 10.0D+00
  else
    exact = ( 1.0D+00 - beta ) * pi &
      / ( 10.0D+00**beta * sin ( pi * beta ) )
  end if

  return
end
subroutine p19_fun ( n, x, fx )

!*****************************************************************************80
!
!! P19_FUN evaluates the integrand for problem 19.
!
!  Integral:
!
!    Integral ( 0 <= x < +oo ) x^(beta-1) / ( 1 + 10 x )^2 dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 84.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: beta = 0.5D+00
  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( beta == 1.0D+00 ) then
      fx(i) = 1.0D+00 / ( 1.0D+00 + 10.0D+00 * x(i) )**2
    else if ( beta < 1.0D+00 .and. x(i) == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = x(i)**( beta - 1.0D+00 ) / ( 1.0D+00 + 10.0D+00 * x(i) )**2
    end if

  end do

  return
end
subroutine p19_title ( title )

!*****************************************************************************80
!
!! P19_TITLE returns the title for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'x^(beta-1) / ( 1 + 10 x )^2'

  return
end
subroutine p20_alpha ( alpha )

!*****************************************************************************80
!
!! P20_ALPHA returns ALPHA for problem 10.
!
!  Discussion:
!
!    ALPHA is the lower, finite limit of integration in the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.0D+00

  return
end
subroutine p20_exact ( exact )

!*****************************************************************************80
!
!! P20_EXACT returns the exact integral for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ), parameter :: beta = 1.0D+00
  real ( kind = 8 ) exact

  exact = &
    ( &
      log ( 1.5D+00 ) / 2.0D+00**beta &
      - 1.0D+00 / 2.0D+00**( beta + 1.0D+00 ) * &
      log ( ( 16.0D+00 + 0.25D+00**beta ) / ( 1.0D+00 + 0.25D+00**beta ) ) &
      - atan ( 2.0D+00**( beta + 2.0D+00 ) ) - atan ( 2.0D+00**beta ) &
    ) / ( 1.0D+00 + 0.25D+00**beta )

  return
end
subroutine p20_fun ( n, x, fx )

!*****************************************************************************80
!
!! P20_FUN evaluates the integrand for problem 20.
!
!  Integral:
!
!    Integral ( 0 <= x < +oo ) 
!      1 / ( 2^beta * ( ( x - 1 )^2 + (1/4)^beta ) * ( x - 2 ) ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 84.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: beta = 1.0D+00
  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( ( x(i) - 1.0D+00 )**2 + 0.25D+00**beta == 0.0D+00 .or. &
         x(i) == 2.0D+00 ) then

      fx(i) = 0.0D+00

    else

      fx(i) = 1.0D+00 / &
        ( 2.0D+00**beta &
        * ( ( x(i) - 1.0D+00 )**2 + 0.25D+00**beta ) &
        * ( x(i) - 2.0D+00 ) )

    end if

  end do

  return
end
subroutine p20_title ( title )

!*****************************************************************************80
!
!! P20_TITLE returns the title for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = '1 / ( 2^beta * ( ( x - 1 )^2 + (1/4)^beta ) * ( x - 2 ) )'

  return
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  implicit none
!
!  Coefficients for minimax approximation over (12, INF).
!
  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( kind = 8 ), parameter :: eps = 2.22D-16
  real ( kind = 8 ) fact
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) sum
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 171.624D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.0D+30
  real ( kind = 8 ), parameter :: xminin = 2.23D-308
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  parity = .false.
  fact = one
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= zero ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= zero ) then

      if ( y1 /= aint ( y1 * half ) * two ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + one

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = one / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < one ) then

      z = y
      y = y + one
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - one

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = zero
    xden = one
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + one
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + one
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - half ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= one ) then
    res = fact / res
  end if

  r8_gamma = res

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
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
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
