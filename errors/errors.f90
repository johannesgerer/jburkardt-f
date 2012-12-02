subroutine dpoly_val ( n, p, x, pval )

!*****************************************************************************80
!
!! DPOLY_VAL evaluates a real polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Input, real ( kind = 8 ) PCOF(0:N), the polynomial coefficients.
!    P(I) is the coefficient of X^I.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) PVAL, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) p(0:n)
  real ( kind = 8 ) pval
  real ( kind = 8 ) x

  pval = p(0)
  do i = 1, n
    pval = pval + p(i) * x**i
  end do

  return
end
subroutine dpoly_val_horner ( n, p, x, pval )

!*****************************************************************************80
!
!! DPOLY_VAL_HORNER evaluates a polynomial using Horner's method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Input, real ( kind = 8 ) PCOF(0:N), the polynomial coefficients.
!    P(I) is the coefficient of X**I.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) PVAL, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) p(0:n)
  real ( kind = 8 ) pval
  real ( kind = 8 ) x

  pval = p(n)
  do i = n - 1, 0, -1
    pval = pval * x + p(i)
  end do

  return
end
subroutine dpoly2_roots ( p, r )

!*****************************************************************************80
!
!! DPOLY2_ROOTS finds the roots of a quadratic polynomial.
!
!  Discussion:
!
!    The standard quadratic formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PCOF(0:2), the polynomial coefficients.
!    P(I) is the coefficient of X**I.
!
!    Output, complex R(2), the roots of the polynomial.
!
  implicit none

  real ( kind = 8 ) disc
  real ( kind = 8 ) p(0:2)
  complex ( kind = 8 ) r(2)

  if ( p(2) == 0.0D+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DPOLY2_ROOTS - Fatal error!'
    write ( *, * ) '  Quadratic coefficient is zero.'
    stop
  end if

  disc = p(1)**2 - 4.0D+00 * p(2) * p(0)

  if ( disc >= 0.0D+00 ) then

    r(1) = cmplx ( 0.5D+00 * ( - p(1) + sqrt ( disc ) ) / p(2), 0.0D+00 )
    r(2) = cmplx ( 0.5D+00 * ( - p(1) - sqrt ( disc ) ) / p(2), 0.0D+00 )

  else if ( disc < 0.0D+00 ) then

    r(1) = cmplx ( - 0.5D+00 * p(1) / p(2), 0.5D+00 * sqrt ( - disc ) / p(2) )
    r(2) = cmplx ( - 0.5D+00 * p(1) / p(2), - 0.5D+00 * sqrt ( - disc ) / p(2) )

  end if

  return
end
subroutine dpoly2_roots2 ( p, r, ierror )

!*****************************************************************************80
!
!! DPOLY2_ROOTS2 finds the roots of a quadratic polynomial.
!
!  Discussion:
!
!    An alternate form of the quadratic formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PCOF(0:2), the polynomial coefficients.
!    P(I) is the coefficient of X**I.
!
!    Output, double complex R(2), the roots of the polynomial.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error;
!    1, an error occurred.
!
  implicit none

  real ( kind = 8 ) disc
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) p(0:2)
  complex ( kind = 8 ) r(2)

  ierror = 0

  if ( p(2) == 0.0D+00 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DPOLY2_ROOTS2 - Fatal error!'
    write ( *, '(a)' ) '  Quadratic coefficient is zero.'
    stop
  end if

  disc = p(1)**2 - 4.0D+00 * p(2) * p(0)

  if ( disc >= 0.0D+00 ) then

    if ( - p(1) + sqrt ( disc ) == 0.0D+00 ) then
      ierror = 1
      r(1) = cmplx ( 0.0D+00, 0.0D+00 )
    else
      r(1) = cmplx ( 2.0D+00 * p(0) / ( - p(1) + sqrt ( disc ) ), 0.0D+00 )
    end if

    if ( - p(1) - sqrt ( disc ) == 0.0D+00 ) then
      ierror = 1
      r(2) = cmplx ( 0.0D+00, 0.0D+00 )
    else
      r(2) = cmplx ( 2.0D+00 * p(0) / ( - p(1) - sqrt ( disc ) ), 0.0D+00 )
    end if
!
!  Need to revise this part of the calculation.
!
  else if ( disc < 0.0D+00 ) then

    r(1) = cmplx ( - 0.5D+00 * p(1) / p(2), + 0.5D+00 * sqrt ( - disc ) / p(2) )
    r(2) = cmplx ( - 0.5D+00 * p(1) / p(2), - 0.5D+00 * sqrt ( - disc ) / p(2) )

  end if

  return
end
function fmin ( ax, bx, f, tol )

!*****************************************************************************80
!
!! FMIN seeks a minimizer of a scalar function of a scalar variable.
!
!  Discussion:
!
!    FMIN seeks an approximation to the point where F attains a minimum on
!    the interval (AX,BX).
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much
!    slower than that for a Fibonacci search.  If F has a continuous
!    second derivative which is positive at the minimum (which is not
!    at AX or BX), then convergence is superlinear, and usually of the
!    order of about 1.324....
!
!    The function F is never evaluated at two points closer together
!    than EPS * ABS ( FMIN ) + (TOL/3), where EPS is approximately the
!    square root of the relative machine precision.  If F is a unimodal
!    function and the computed values of F are always unimodal when
!    separated by at least EPS * ABS ( XSTAR ) + (TOL/3), then FMIN
!    approximates the abcissa of the global minimum of F on the
!    interval AX, BX with an error less than 3 * EPS * ABS ( FMIN ) + TOL.
!    If F is not unimodal, then FMIN may approximate a local, but
!    perhaps non-global, minimum to the same accuracy.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization without Derivatives,
!    Prentice Hall, 1973.
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters
!
!    Input/output, real ( kind = 4 ) AX, BX.  On input, the left and right
!    endpoints of the initial interval.  On output, the lower and upper bounds
!    for the minimizer.
!
!    Input, external F, a real function of the form
!      function f ( x )
!      real ( kind = 4 ) f
!      real ( kind = 4 ) x
!    which evaluates F(X) for any X in the interval (AX,BX).
!
!    Input, real ( kind = 4 ) TOL, the desired length of the interval of
!    uncertainty of the final result ( >= 0.0)
!
!    Output, real ( kind = 4 ) FMIN, the abcissa approximating the minimizer
!    of f.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) ax
  real ( kind = 4 ) b
  real ( kind = 4 ) bx
  real ( kind = 4 ) c
  real ( kind = 4 ) d
  real ( kind = 4 ) e
  real ( kind = 4 ) eps
  real ( kind = 4 ), external :: f
  real ( kind = 4 ) fmin
  real ( kind = 4 ) fu
  real ( kind = 4 ) fv
  real ( kind = 4 ) fw
  real ( kind = 4 ) fx
  real ( kind = 4 ) p
  real ( kind = 4 ) q
  real ( kind = 4 ) r
  real ( kind = 4 ) tol
  real ( kind = 4 ) tol1
  real ( kind = 4 ) tol2
  real ( kind = 4 ) u
  real ( kind = 4 ) v
  real ( kind = 4 ) w
  real ( kind = 4 ) x
  real ( kind = 4 ) xm

  c = 0.5E+00 * ( 3.0E+00 - sqrt ( 5.0E+00 ) )
!
!  C is the squared inverse of the golden ratio.
!
!  EPS is the square root of the relative machine precision.
!
  eps = sqrt ( epsilon ( eps ) )
!
!  Initialization.
!
  a = ax
  b = bx
  v = a + c * ( b - a )
  w = v
  x = v
  e = 0.0E+00
  fx = f(x)
  fv = fx
  fw = fx
!
!  Main loop starts here.
!
20 continue

  xm = 0.5E+00 * ( a + b )
  tol1 = eps * abs ( x ) + tol / 3.0E+00
  tol2 = 2.0E+00 * tol1
!
!  Check the stopping criterion.
!
  if ( abs ( x - xm ) <= ( tol2 - 0.5E+00 * ( b - a ) ) ) then
    fmin = x
    return
  end if
!
!  Is golden-section necessary?
!
  if ( abs ( e ) <= tol1 ) then
    go to 40
  end if
!
!  Fit a parabola.
!
  r = ( x - w ) * ( fx - fv )
  q = ( x - v ) * ( fx - fw )
  p = ( x - v ) * q - ( x - w ) * r
  q = 2.0E+00 * ( q - r )
  if ( q > 0.0E+00 ) then
    p = -p
  end if
  q = abs ( q )
  r = e
  e = d
!
!  Is a parabola acceptable?
!
   30 continue

  if ( abs ( p ) >= abs ( 0.5E+00 * q * r ) ) then
    go to 40
  end if

  if ( p <= q * ( a - x ) ) then
    go to 40
  end if

  if ( p >= q * ( b - x ) ) then
    go to 40
  end if
!
!  A parabolic interpolation step
!
  d = p / q
  u = x + d
!
!  F must not be evaluated too close to AX or BX.
!
  if ( ( u - a ) < tol2 ) then
    d = sign ( tol1, xm - x )
  end if

  if ( ( b - u ) < tol2 ) then
    d = sign ( tol1, xm - x )
  end if

  go to 50
!
!  A golden-section step.
!
   40 continue

  if ( x >= xm ) then
    e = a - x
  else
    e = b - x
  end if

  d = c * e
!
!  F must not be evaluated too close to X.
!
   50 continue

  if ( abs ( d ) >= tol1 ) then
    u = x + d
  end if

  if ( abs ( d ) < tol1 ) then
    u = x + sign ( tol1, d )
  end if

  fu = f(u)
!
!  Update  a, b, v, w, and x
!
  if ( fu <= fx ) then

    if ( u >= x ) then
      a = x
    else
      b = x
    end if

    v = w
    fv = fw
    w = x
    fw = fx
    x = u
    fx = fu
    go to 20

  end if

60 continue

  if ( u < x ) then
    a = u
  else
    b = u
  end if

  if ( fu <= fw ) then
    go to 70
  end if

  if ( w == x ) then
    go to 70
  end if

  if ( fu <= fv ) then
    go to 80
  end if

  if ( v == x ) then
    go to 80
  end if

  if ( v == w ) then
    go to 80
  end if

  go to 20

   70 continue

  v = w
  fv = fw
  w = u
  fw = fu
  go to 20
   80 continue

  v = u
  fv = fu
  go to 20

end
function isamax ( n, x, incx )

!*****************************************************************************80
!
!! ISAMAX finds the index of the vector element of maximum absolute value.
!
!  Modified:
!
!    08 April 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 4 ) X(*), the vector to be examined.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of SX.
!
!    Output, integer ( kind = 4 ) ISAMAX, the index of the element of SX of
!    maximum absolute value.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) isamax
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real ( kind = 4 ) samax
  real ( kind = 4 ) x(*)

  if ( n <= 0 ) then

    isamax = 0

  else if ( n == 1 ) then

    isamax = 1

  else if ( incx == 1 ) then

    isamax = 1
    samax = abs ( x(1) )

    do i = 2, n

      if ( abs ( x(i) ) > samax ) then
        isamax = i
        samax = abs ( x(i) )
      end if

    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    isamax = 1
    samax = abs ( x(ix) )

    ix = ix + incx

    do i = 2, n
      if ( abs ( x(ix) ) > samax ) then
        isamax = i
        samax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
function lcm_12n ( n )

!*****************************************************************************80
!
!! LCM_12N computes the least common multiple of the integers 1 through N.
!
!  Examples:
!
!    N    LCM_12N
!
!    1          1
!    2          2
!    3          3
!    4         12
!    5         60
!    6         60
!    7        420
!    8        840
!    9       2520
!   10       2520
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the value of N.
!
!    Output, integer ( kind = 4 ) LCM_12N, the least common multiple of
!    the integers 1 to N.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imult
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lcm_12n
  integer ( kind = 4 ) n

  lcm_12n = 1

  do i = 2, n

    imult = i

    do j = 1, i-1

      if ( mod ( imult, (i-j) )   ==   0 ) then
        imult = imult / ( i - j )
      end if

    end do

    lcm_12n = lcm_12n * imult

  end do

  return
end
subroutine matrix_exponential_taylor ( n, a, a_exp )

!*****************************************************************************80
!
!! MAXTRIX_EXPONENTIAL_TAYLOR uses a Taylor series for the matrix exponential.
!
!  Discussion:
!
!    This is a simple method, but NOT recommended.   It is easy to
!    find examples for which this method fails.
!
!    Successive terms in the Taylor series are added.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cleve Moler and Charles Van Loan,
!    19 Dubious Ways to Compute the Exponential of a Matrix, 25 Years Later,
!    SIAM Review,
!    Volume 45, Number 1, pages 3-49, March 2003.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real A(N,N), the matrix whose exponential is desired.
!
!    Output, real A_EXP(N,N), a Taylor estimate for the matrix exponential.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) a(n,n)
  real ( kind = 4 ) a_exp(n,n)
  real ( kind = 4 ) b_exp(n,n)
  real ( kind = 4 ) a_k(n,n)
  real ( kind = 4 ) fact_k
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 4 ), parameter :: tol = 0.0E+00

  a_exp(1:n,1:n) = 0.0E+00

  do i = 1, n
    a_exp(i,i) = 1.0E+00
  end do

  a_k(1:n,1:n) = a(1:n,1:n)

  k = 1

  do

    b_exp(1:n,1:n) = a_exp(1:n,1:n)

    a_exp(1:n,1:n) = a_exp(1:n,1:n) + a_k(1:n,1:n)

    b_exp(1:n,1:n) = abs ( b_exp(1:n,1:n) - a_exp(1:n,1:n) )

    if ( all ( b_exp(1:n,1:n) <= tol ) ) then
      exit
    end if

    k = k + 1

    a_k(1:n,1:n) = matmul ( a_k(1:n,1:n), a(1:n,1:n) ) / real ( k )

  end do

  return
end
subroutine rpoly_val ( n, p, x, pval )

!*****************************************************************************80
!
!! RPOLY_VAL evaluates a real polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Input, real ( kind = 4 ) PCOF(0:N), the polynomial coefficients.
!    P(I) is the coefficient of X**I.
!
!    Input, real ( kind = 4 ) X, the point at which the polynomial is to
!    be evaluated.
!
!    Output, real ( kind = 4 ) PVAL, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 4 ) p(0:n)
  real ( kind = 4 ) pval
  real ( kind = 4 ) x

  pval = p(0)
  do i = 1, n
    pval = pval + p(i) * x**i
  end do

  return
end
subroutine rpoly_val_horner ( n, p, x, pval )

!*****************************************************************************80
!
!! RPOLY_VAL_HORNER evaluates a real polynomial using Horner's method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Input, real ( kind = 4 ) PCOF(0:N), the polynomial coefficients.
!    P(I) is the coefficient of X**I.
!
!    Input, real ( kind = 4 ) X, the point at which the polynomial is to
!    be evaluated.
!
!    Output, rea ( kind = 4 )l PVAL, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 4 ) p(0:n)
  real ( kind = 4 ) pval
  real ( kind = 4 ) x

  pval = p(n)
  do i = n - 1, 0, -1
    pval = pval * x + p(i)
  end do

  return
end
subroutine rpoly2_roots ( p, r )

!*****************************************************************************80
!
!! RPOLY2_ROOTS finds the roots of a quadratic polynomial.
!
!  Discussion:
!
!    The standard quadratic formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) PCOF(0:2), the polynomial coefficients.
!    P(I) is the coefficient of X**I.
!
!    Output, complex ( kind = 4 ) R(2), the roots of the polynomial.
!
  implicit none

  real ( kind = 4 ) disc
  real ( kind = 4 ) p(0:2)
  complex ( kind = 4 ) r(2)

  if ( p(2) == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RPOLY2_ROOTS - Fatal error!'
    write ( *, '(a)' ) '  Quadratic coefficient is zero.'
    stop
  end if

  disc = p(1)**2 - 4.0E+00 * p(2) * p(0)

  if ( disc >= 0.0E+00 ) then

    r(1) = cmplx ( 0.5E+00 * ( - p(1) + sqrt ( disc ) ) / p(2), 0.0E+00 )
    r(2) = cmplx ( 0.5E+00 * ( - p(1) - sqrt ( disc ) ) / p(2), 0.0E+00 )

  else if ( disc < 0.0E+00 ) then

    r(1) = cmplx ( - 0.5E+00 * p(1) / p(2), 0.5E+00 * sqrt ( - disc ) / p(2) )
    r(2) = cmplx ( - 0.5E+00 * p(1) / p(2), - 0.5E+00 * sqrt ( - disc ) / p(2) )

  end if

  return
end
subroutine rpoly2_roots2 ( p, r, ierror )

!*****************************************************************************80
!
!! RPOLY2_ROOTS2 finds the roots of a quadratic polynomial.
!
!  Discussion:
!
!    An alternate form of the quadratic formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real PCOF(0:2), the polynomial coefficients.
!    P(I) is the coefficient of X**I.
!
!    Output, complex R(2), the roots of the polynomial.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error;
!    1, an error occurred.
!
  implicit none

  real ( kind = 4 ) disc
  integer ( kind = 4 ) ierror
  real ( kind = 4 ) p(0:2)
  complex ( kind = 4 ) r(2)

  ierror = 0

  if ( p(2) == 0.0E+00 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RPOLY2_ROOTS2 - Fatal error!'
    write ( *, '(a)' ) '  Quadratic coefficient is zero.'
    stop
  end if

  disc = p(1)**2 - 4.0E+00 * p(2) * p(0)

  if ( disc >= 0.0E+00 ) then

    if ( - p(1) + sqrt ( disc ) == 0.0E+00 ) then
      ierror = 1
      r(1) = cmplx ( 0.0E+00, 0.0E+00 )
    else
      r(1) = cmplx ( 2.0E+00 * p(0) / ( - p(1) + sqrt ( disc ) ), 0.0E+00 )
    end if

    if ( - p(1) - sqrt ( disc ) == 0.0E+00 ) then
      ierror = 1
      r(2) = cmplx ( 0.0E+00, 0.0E+00 )
    else
      r(2) = cmplx ( 2.0E+00 * p(0) / ( - p(1) - sqrt ( disc ) ), 0.0E+00 )
    end if
!
!  Need to revise this part of the calculation.
!
  else if ( disc < 0.0E+00 ) then

    r(1) = cmplx ( - 0.5E+00 * p(1) / p(2), + 0.5E+00 * sqrt ( - disc ) / p(2) )
    r(2) = cmplx ( - 0.5E+00 * p(1) / p(2), - 0.5E+00 * sqrt ( - disc ) / p(2) )

  end if

  return
end
function samax ( n, x, incx )

!*****************************************************************************80
!
!! SAMAX returns the maximum absolute value of the entries in a vector.
!
!  Modified:
!
!    08 April 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 4 ) X(*), the vector to be examined.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Output, real ( kind = 4 ) SAMAX, the maximum absolute value of an
!    element of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real ( kind = 4 ) samax
  real ( kind = 4 ) x(*)

  if ( n <= 0 ) then

    samax = 0.0E+00

  else if ( n == 1 ) then

    samax = abs ( x(1) )

  else if ( incx == 1 ) then

    samax = abs ( x(1) )

    do i = 2, n
      if ( abs ( x(i) ) > samax ) then
        samax = abs ( x(i) )
      end if
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    samax = abs ( x(ix) )
    ix = ix + incx

    do i = 2, n
      if ( abs ( x(ix) ) > samax ) then
        samax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
subroutine saxpy ( n, sa, x, incx, y, incy )

!*****************************************************************************80
!
!! SAXPY adds a constant times one vector to another.
!
!  Modified:
!
!    08 April 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 4 ) SA, the multiplier.
!
!    Input, real ( kind = 4 ) X(*), the vector to be scaled and added to Y.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Input/output, real ( kind = 4 ) Y(*), the vector to which a multiple
!    of X is to be added.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive
!    entries of Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) n
  real ( kind = 4 ) sa
  real ( kind = 4 ) x(*)
  real ( kind = 4 ) y(*)

  if ( n <= 0 ) then

  else if ( sa == 0.0E+00 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    y(1:n) = y(1:n) + sa * x(1:n)

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      y(iy) = y(iy) + sa * x(ix)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine scopy ( n, x, incx, y, incy )

!*****************************************************************************80
!
!! SCOPY copies one real vector into another.
!
!  Modified:
!
!    08 April 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 4 ) X(*), the vector to be copied into Y.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Output, real ( kind = 4 ) Y(*), the copy of X.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive
!    elements of Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) n
  real ( kind = 4 ) x(*)
  real ( kind = 4 ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    y(1:n) = x(1:n)

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      y(iy) = x(ix)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
function sdot ( n, x, incx, y, incy )

!*****************************************************************************80
!
!! SDOT forms the dot product of two vectors.
!
!  Modified:
!
!    02 June 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 4 ) X(*), one of the vectors to be multiplied.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Input, real ( kind = 4 ) Y(*), one of the vectors to be multiplied.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive
!    elements of Y.
!
!    Output, real ( kind = 4 ) SDOT, the dot product of X and Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) n
  real ( kind = 4 ) sdot
  real ( kind = 4 ) stemp
  real ( kind = 4 ) x(*)
  real ( kind = 4 ) y(*)

  if ( n <= 0 ) then

    sdot = 0.0E+00

  else if ( incx == 1 .and. incy == 1 ) then

    sdot = dot_product ( x(1:n), y(1:n) )

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    stemp = 0.0E+00
    do i = 1, n
      stemp = stemp + x(ix) * y(iy)
      ix = ix + incx
      iy = iy + incy
    end do

    sdot = stemp

  end if

  return
end
function sdsdot ( n, x, incx, y, incy )

!*****************************************************************************80
!
!! SDSDOT forms the dot product of two vectors using higher precision.
!
!  Modified:
!
!    02 June 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 4 ) X(*), one of the vectors to be multiplied.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Input, real ( kind = 4 ) Y(*), one of the vectors to be multiplied.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive
!    elements of Y.
!
!    Output, real ( kind = 4 ) SDSDOT, the dot product of X and Y.
!
  implicit none

  real ( kind = 8 ) dsdot
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) n
  real ( kind = 4 ) sdsdot
  real ( kind = 4 ) x(*)
  real ( kind = 4 ) y(*)

  if ( n <= 0 ) then

    dsdot = 0.0D+00

  else if ( incx == 1 .and. incy == 1 ) then

    dsdot = dot_product ( dble ( x(1:n) ), dble ( y(1:n) ) )

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    dsdot = 0.0D+00
    do i = 1, n
      dsdot = dsdot + dble ( x(ix) ) * dble ( y(iy) )
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  sdsdot = sngl ( dsdot )

  return
end
subroutine sgedi ( a, lda, n, ipvt, det, work, job )

!*****************************************************************************80
!
!! SGEDI: determinant and inverse of a matrix factored by SGECO or SGEFA.
!
!  Discussion:
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal and the inverse is requested.
!    It will not occur if the subroutines are called correctly
!    and if SGECO has set RCOND > 0.0E+00 or SGEFA has set INFO == 0.
!
!  Reference:
!
!    Dongarra, Moler, Bunch and Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(LDA,N), on input, the factored matrix.
!    as output by SGECO or SGEFA.  On output, contains the inverse
!    matrix if requested.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from SGECO or SGEFA.
!
!    Workspace, real ( kind = 4 ) WORK(N).
!
!    Output, real ( kind = 4 ) DET(2), the determinant of original matrix if
!    requested.  determinant = DET(1) * 10.0**DET(2)
!    with  1.0E+00 <= abs ( DET(1) ) < 10.0E+00
!    or DET(1) == 0.0E+00.
!
!    Input, integer ( kind = 4 ) JOB, specifies what is to be computed.
!    11, both determinant and inverse.
!    01, inverse only.
!    10, determinant only.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) l
  real ( kind = 4 ) t
  real ( kind = 4 ), parameter :: ten = 10.0E+00
  real ( kind = 4 ) work(n)
!
!  Compute the determinant.
!
  if ( job / 10 /= 0 ) then

     det(1) = 1.0E+00
     det(2) = 0.0E+00

     do i = 1, n

        if ( ipvt(i) /= i ) then
          det(1) = - det(1)
        end if

        det(1) = a(i,i) * det(1)

        if ( det(1) == 0.0E+00 ) then
          exit
        end if

        do while ( abs ( det(1) ) < 1.0E+00 )
          det(1) = ten * det(1)
          det(2) = det(2) - 1.0E+00
        end do

        do while ( abs ( det(1) ) >= ten )
          det(1) = det(1) / ten
          det(2) = det(2) + 1.0E+00
        end do

    end do

  end if
!
!  Compute inverse(U).
!
  if ( mod ( job, 10 ) /= 0 ) then

     do k = 1, n

        a(k,k) = 1.0E+00 / a(k,k)
        t = - a(k,k)
        call sscal ( k-1, t, a(1,k), 1 )

        do j = k+1, n
           t = a(k,j)
           a(k,j) = 0.0E+00
           call saxpy ( k, t, a(1,k), 1, a(1,j), 1 )
        end do

     end do
!
!  Form inverse(U) * inverse(L).
!
     do k = n-1, 1, -1

        do i = k+1, n
           work(i) = a(i,k)
           a(i,k) = 0.0E+00
        end do

        do j = k+1, n
           t = work(j)
           call saxpy ( n, t, a(1,j), 1, a(1,k), 1 )
        end do

        l = ipvt(k)
        if ( l /= k ) then
          call sswap ( n, a(1,k), 1, a(1,l), 1 )
        end if

     end do

  end if

  return
end
subroutine sgefa ( a, lda, n, ipvt, info )

!*****************************************************************************80
!
!! SGEFA factors a real matrix.
!
!  Modified:
!
!    07 March 2001
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(LDA,N).
!    On intput, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers used to
!    obtain it.  The factorization can be written A=L*U, where L is a product
!    of permutation and unit lower triangular matrices, and U is upper
!    triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity indicator.
!    0, normal value.
!    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
!    but it does indicate that SGESL or SGEDI will divide by zero if called.
!    Use RCOND in SGECO for a reliable indication of singularity.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) isamax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 4 ) t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1
!
!  Find L = pivot index.
!
    l = isamax ( n-k+1, a(k,k), 1 ) + k - 1
    ipvt(k) = l
!
!  Zero pivot implies this column already triangularized.
!
    if ( a(l,k) == 0.0E+00 ) then
      info = k
      cycle
    end if
!
!  Interchange if necessary.
!
    if ( l /= k ) then
      t = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Compute multipliers.
!
    t = -1.0E+00 / a(k,k)
    call sscal ( n-k, t, a(k+1,k), 1 )
!
!  Row elimination with column indexing.
!
    do j = k+1, n
      t = a(l,j)
      if ( l /= k ) then
        a(l,j) = a(k,j)
        a(k,j) = t
      end if
      call saxpy ( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )
    end do

  end do

  ipvt(n) = n

  if ( a(n,n) == 0.0E+00 ) then
    info = n
  end if

  return
end
subroutine sgesl ( a, lda, n, ipvt, b, job )

!*****************************************************************************80
!
!! SGESL solves a real general linear system A * X = B.
!
!  Discussion:
!
!    SGESL can solve either of the systems A * X = B or transpose ( A ) * X = B.
!
!    The system matrix must have been factored by SGECO or SGEFA.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if SGECO has set RCOND > 0.0E+00
!    or SGEFA has set INFO == 0.
!
!  Modified:
!
!    07 March 2001
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(LDA,N), the output from SGECO or SGEFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from SGECO or SGEFA.
!
!    Input/output, real ( kind = 4 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB.
!    0, solve A * X = B;
!    nonzero, solve A' * X = B.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 4 ) sdot
  real ( kind = 4 ) t
!
!  Solve A * X = B.
!
  if ( job == 0 ) then

    do k = 1, n-1

      l = ipvt(k)
      t = b(l)

      if ( l /= k ) then
        b(l) = b(k)
        b(k) = t
      end if

      call saxpy ( n-k, t, a(k+1,k), 1, b(k+1), 1 )

    end do

    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      t = -b(k)
      call saxpy ( k-1, t, a(1,k), 1, b(1), 1 )
    end do

  else
!
!  Solve A' * X = B.
!
    do k = 1, n
      t = sdot ( k-1, a(1,k), 1, b(1), 1 )
      b(k) = ( b(k) - t ) / a(k,k)
    end do

    do k = n-1, 1, -1

      b(k) = b(k) + sdot ( n-k, a(k+1,k), 1, b(k+1), 1 )
      l = ipvt(k)

      if ( l /= k ) then
        t = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

  return
end
function rvec_norm2 ( n, x, incx )

!*****************************************************************************80
!
!! RVEC_NORM2 computes the Euclidean norm of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 4 ) X(*), the vector whose norm is to be computed.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Output, real ( kind = 4 ) SNRM2, the Euclidean norm of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real ( kind = 4 ) samax
  real ( kind = 4 ) rvec_norm2
  real ( kind = 4 ) stemp
  real ( kind = 4 ) x(*)
  real ( kind = 4 ) xmax

  if ( n <= 0 ) then

    rvec_norm2 = 0.0E+00

  else

    xmax = samax ( n, x, incx )

    if ( xmax == 0.0E+00 ) then

      rvec_norm2 = 0.0E+00

    else

      if ( incx >= 0 ) then
        ix = 1
      else
        ix = ( - n + 1 ) * incx + 1
      end if

      stemp = 0.0E+00
      do i = 1, n
        stemp = stemp + ( x(ix) / xmax )**2
        ix = ix + incx
      end do

      rvec_norm2 = xmax * sqrt ( stemp )

    end if

  end if

  return
end
subroutine sqrdc ( a, lda, n, p, qraux, jpvt, work, job )

!*****************************************************************************80
!
!! SQRDC computes the QR factorization of a real rectangular matrix.
!
!  Discussion:
!
!    SQRDC uses Householder transformations.
!
!    Column pivoting based on the 2-norms of the reduced columns may be
!    performed at the user's option.
!
!  Reference:
!
!    Dongarra, Moler, Bunch and Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(LDA,P).  On input, the N by P matrix
!    whose decomposition is to be computed.  On output, A contains in
!    its upper triangle the upper triangular matrix R of the QR
!    factorization.  Below its diagonal A contains information from
!    which the orthogonal part of the decomposition can be recovered.
!    Note that if pivoting has been requested, the decomposition is not that
!    of the original matrix A but that of A with its columns permuted
!    as described by JPVT.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer ( kind = 4 ) N, the number of rows of the matrix A.
!
!    Input, integer ( kind = 4 ) P, the number of columns of the matrix A.
!
!    Output, real ( kind = 4 ) QRAUX(P), contains further information required
!    to recover the orthogonal part of the decomposition.
!
!    Input/output, integer ( kind = 4 ) JPVT(P).  On input, JPVT contains
!    integers that control the selection of the pivot columns.  The K-th
!    column A(*,K) of A is placed in one of three classes according to the
!    value of JPVT(K).
!      > 0, then A(K) is an initial column.
!      = 0, then A(K) is a free column.
!      < 0, then A(K) is a final column.
!    Before the decomposition is computed, initial columns are moved to
!    the beginning of the array A and final columns to the end.  Both
!    initial and final columns are frozen in place during the computation
!    and only free columns are moved.  At the K-th stage of the
!    reduction, if A(*,K) is occupied by a free column it is interchanged
!    with the free column of largest reduced norm.  JPVT is not referenced
!    if JOB == 0.  On output, JPVT(K) contains the index of the column of the
!    original matrix that has been interchanged into the K-th column, if
!    pivoting was requested.
!
!    Workspace, real ( kind = 4 ) WORK(P).  WORK is not referenced if JOB == 0.
!
!    Input, integer ( kind = 4 ) JOB, initiates column pivoting.
!    0, no pivoting is done.
!    nonzero, pivoting is done.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p

  real ( kind = 4 ) a(lda,p)
  integer ( kind = 4 ) jpvt(p)
  real ( kind = 4 ) qraux(p)
  real ( kind = 4 ) work(p)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) job
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lup
  integer ( kind = 4 ) maxj
  real ( kind = 4 ) maxnrm
  real ( kind = 4 ) nrmxl
  integer ( kind = 4 ) pl
  integer ( kind = 4 ) pu
  real ( kind = 4 ) sdot
  real ( kind = 4 ) rvec_norm2
  logical swapj
  real ( kind = 4 ) t
  real ( kind = 4 ) tt

  pl = 1
  pu = 0
!
!  If pivoting is requested, rearrange the columns.
!
  if ( job /= 0 ) then

    do j = 1, p

      swapj = jpvt(j) > 0

      if ( jpvt(j) < 0 ) then
        jpvt(j) = - j
      else
        jpvt(j) = j
      end if

      if ( swapj ) then

        if ( j /= pl ) then
          call sswap ( n, a(1,pl), 1, a(1,j), 1 )
        end if

        jpvt(j) = jpvt(pl)
        jpvt(pl) = j
        pl = pl + 1

      end if

    end do

    pu = p

    do j = p, 1, -1

      if ( jpvt(j) < 0 ) then

        jpvt(j) = - jpvt(j)

        if ( j /= pu ) then
          call sswap ( n, a(1,pu), 1, a(1,j), 1 )
          jp = jpvt(pu)
          jpvt(pu) = jpvt(j)
          jpvt(j) = jp
        end if

        pu = pu - 1

      end if

    end do

  end if
!
!  Compute the norms of the free columns.
!
  do j = pl, pu
    qraux(j) = rvec_norm2 ( n, a(1,j), 1 )
  end do

  work(pl:pu) = qraux(pl:pu)
!
!  Perform the Householder reduction of A.
!
  lup = min ( n, p )

  do l = 1, lup
!
!  Bring the column of largest norm into the pivot position.
!
    if ( l >= pl .and. l < pu ) then

      maxnrm = 0.0E+00
      maxj = l
      do j = l, pu
        if ( qraux(j) > maxnrm ) then
          maxnrm = qraux(j)
          maxj = j
        end if
      end do

      if ( maxj /= l ) then
        call sswap ( n, a(1,l), 1, a(1,maxj), 1 )
        qraux(maxj) = qraux(l)
        work(maxj) = work(l)
        jp = jpvt(maxj)
        jpvt(maxj) = jpvt(l)
        jpvt(l) = jp
      end if

    end if
!
!  Compute the Householder transformation for column L.
!
    qraux(l) = 0.0E+00

    if ( l /= n ) then

      nrmxl = rvec_norm2 ( n-l+1, a(l,l), 1 )

      if ( nrmxl /= 0.0E+00 ) then

        if ( a(l,l) /= 0.0E+00 ) then
          nrmxl = sign ( nrmxl, a(l,l) )
        end if

        call sscal ( n-l+1, 1.0E+00 / nrmxl, a(l,l), 1 )
        a(l,l) = 1.0E+00 + a(l,l)
!
!  Apply the transformation to the remaining columns, updating the norms.
!
        do j = l + 1, p

          t = - sdot ( n-l+1, a(l,l), 1, a(l,j), 1 ) / a(l,l)
          call saxpy ( n-l+1, t, a(l,l), 1, a(l,j), 1 )

          if ( j >= pl .and. j <= pu ) then

            if ( qraux(j) /= 0.0E+00 ) then

              tt = 1.0E+00 - ( abs ( a(l,j) ) / qraux(j) )**2
              tt = max ( tt, 0.0E+00 )
              t = tt
              tt = 1.0E+00 + 0.05E+00 * tt * ( qraux(j) / work(j) )**2

              if ( tt /= 1.0E+00 ) then
                qraux(j) = qraux(j) * sqrt ( t )
              else
                qraux(j) = rvec_norm2 ( n-l, a(l+1,j), 1 )
                work(j) = qraux(j)
              end if

            end if

          end if

        end do
!
!  Save the transformation.
!
        qraux(l) = a(l,l)
        a(l,l) = - nrmxl

      end if

    end if

  end do

  return
end
subroutine sqrsl ( a, lda, n, k, qraux, y, qy, qty, b, rsd, ab, job, info )

!*****************************************************************************80
!
!! SQRSL computes transformations, projections, and least squares solutions.
!
!  Discussion:
!
!    SQRSL requires the output of SQRDC.
!
!    For K <= min(N,P), let AK be the matrix
!
!      AK = ( A(JPVT(1)), A(JPVT(2)), ..., A(JPVT(K)) )
!
!    formed from columns JPVT(1), ..., JPVT(K) of the original
!    N by P matrix A that was input to SQRDC.  If no pivoting was
!    done, AK consists of the first K columns of A in their
!    original order.  SQRDC produces a factored orthogonal matrix Q
!    and an upper triangular matrix R such that
!
!      AK = Q * (R)
!               (0)
!
!    This information is contained in coded form in the arrays
!    A and QRAUX.
!
!    The parameters QY, QTY, B, RSD, and AB are not referenced
!    if their computation is not requested and in this case
!    can be replaced by dummy variables in the calling program.
!    To save storage, the user may in some cases use the same
!    array for different parameters in the calling sequence.  A
!    frequently occuring example is when one wishes to compute
!    any of B, RSD, or AB and does not need Y or QTY.  In this
!    case one may identify Y, QTY, and one of B, RSD, or AB, while
!    providing separate arrays for anything else that is to be
!    computed.
!
!    Thus the calling sequence
!
!      call sqrsl ( a, lda, n, k, qraux, y, dum, y, b, y, dum, 110, info )
!
!    will result in the computation of B and RSD, with RSD
!    overwriting Y.  More generally, each item in the following
!    list contains groups of permissible identifications for
!    a single calling sequence.
!
!      1. (Y,QTY,B) (RSD) (AB) (QY)
!
!      2. (Y,QTY,RSD) (B) (AB) (QY)
!
!      3. (Y,QTY,AB) (B) (RSD) (QY)
!
!      4. (Y,QY) (QTY,B) (RSD) (AB)
!
!      5. (Y,QY) (QTY,RSD) (B) (AB)
!
!      6. (Y,QY) (QTY,AB) (B) (RSD)
!
!    In any group the value returned in the array allocated to
!    the group corresponds to the last member of the group.
!
!  Reference:
!
!    Dongarra, Moler, Bunch and Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(LDA,P), contains the output of SQRDC.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the number of rows of the matrix AK.  It
!    must have the same value as N in SQRDC.
!
!    Input, integer ( kind = 4 ) K, the number of columns of the matrix AK.  K
!    must not be greater than min(N,P), where P is the same as in the
!    calling sequence to SQRDC.
!
!    Input, real ( kind = 4 ) QRAUX(P), the auxiliary output from SQRDC.
!
!    Input, real ( kind = 4 ) Y(N), a vector to be manipulated by SQRSL.
!
!    Output, real ( kind = 4 ) QY(N), contains Q * Y, if requested.
!
!    Output, real ( kind = 4 ) QTY(N), contains Q' * Y, if requested.
!
!    Output, real ( kind = 4 ) B(K), the solution of the least squares problem
!      minimize norm2 ( Y - AK * B),
!    if its computation has been requested.  Note that if pivoting was
!    requested in SQRDC, the J-th component of B will be associated with
!    column JPVT(J) of the original matrix A that was input into SQRDC.
!
!    Output, real ( kind = 4 ) RSD(N), the least squares residual Y - AK * B,
!    if its computation has been requested.  RSD is also the orthogonal
!    projection of Y onto the orthogonal complement of the column space
!    of AK.
!
!    Output, real ( kind = 4 ) AB(N), the least squares approximation Ak * B,
!    if its computation has been requested.  AB is also the orthogonal
!    projection of Y onto the column space of A.
!
!    Input, integer ( kind = 4 ) JOB, specifies what is to be computed.  JOB has
!    the decimal expansion ABCDE, with the following meaning:
!
!      if A /= 0, compute QY.
!      if B /= 0, compute QTY.
!      if C /= 0, compute QTY and B.
!      if D /= 0, compute QTY and RSD.
!      if E /= 0, compute QTY and AB.
!
!    Note that a request to compute B, RSD, or AB automatically triggers
!    the computation of QTY, for which an array must be provided in the
!    calling sequence.
!
!    Output, integer ( kind = 4 ) INFO, is zero unless the computation of B has
!    been requested and R is exactly singular.  In this case, INFO is the
!    index of the first zero diagonal element of R, and B is left unaltered.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,*)
  real ( kind = 4 ) ab(n)
  real ( kind = 4 ) b(k)
  logical cab
  logical cb
  logical cqty
  logical cqy
  logical cr
  logical cxb
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) job
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) kp1
  real ( kind = 4 ) qraux(*)
  real ( kind = 4 ) qty(n)
  real ( kind = 4 ) qy(n)
  real ( kind = 4 ) rsd(n)
  real ( kind = 4 ) sdot
  real ( kind = 4 ) t
  real ( kind = 4 ) temp
  real ( kind = 4 ) y(n)
!
!  set info flag.
!
  info = 0
!
!  Determine what is to be computed.
!
  cqy =        job / 10000         /= 0
  cqty = mod ( job,  10000 )       /= 0
  cb =   mod ( job,   1000 ) / 100 /= 0
  cr =   mod ( job,    100 ) /  10 /= 0
  cab =  mod ( job,     10 )       /= 0

  ju = min ( k, n-1 )
!
!  Special action when N = 1.
!
  if ( ju == 0 ) then

    if ( cqy ) then
      qy(1) = y(1)
    end if

    if ( cqty ) then
      qty(1) = y(1)
    end if

    if ( cab ) then
      ab(1) = y(1)
     end if

    if ( cb ) then

      if ( a(1,1) == 0.0E+00 ) then
        info = 1
      else
        b(1) = y(1) / a(1,1)
      end if

    end if

    if ( cr ) then
      rsd(1) = 0.0E+00
    end if

    return

  end if
!
!  Set up to compute QY or QTY.
!
  if ( cqy ) then
    qy(1:n) = y(1:n)
  end if

  if ( cqty ) then
    qty(1:n) = y(1:n)
  end if
!
!  Compute QY.
!
  if ( cqy ) then

    do jj = 1, ju

      j = ju - jj + 1

      if ( qraux(j) /= 0.0E+00 ) then
        temp = a(j,j)
        a(j,j) = qraux(j)
        t = - sdot ( n-j+1, a(j,j), 1, qy(j), 1 ) / a(j,j)
        call saxpy ( n-j+1, t, a(j,j), 1, qy(j), 1 )
        a(j,j) = temp
      end if

    end do

  end if
!
!  Compute Q'*Y.
!
     if ( cqty ) then

        do j = 1, ju
           if ( qraux(j) /= 0.0E+00 ) then
              temp = a(j,j)
              a(j,j) = qraux(j)
              t = - sdot ( n-j+1, a(j,j), 1, qty(j), 1 ) / a(j,j)
              call saxpy ( n-j+1, t, a(j,j), 1, qty(j), 1 )
              a(j,j) = temp
           end if
        end do

     end if
!
!  Set up to compute B, RSD, or AB.
!
     if ( cb ) then
       b(1:k) = qty(1:k)
     end if

     kp1 = k + 1

     if ( cab ) then
       ab(1:k) = qty(1:k)
     end if

     if ( cr .and. k < n ) then
       rsd(k+1:n) = qty(k+1:n)
     end if

     if ( cab .and. k+1 <= n ) then
        ab(k+1:n) = 0.0E+00
     end if

     if ( cr ) then
        rsd(1:k) = 0.0E+00
     end if
!
!  Compute B.
!
     if ( cb ) then

        do jj = 1, k

           j = k - jj + 1

           if ( a(j,j) == 0.0E+00 ) then
              info = j
              exit
           end if

           b(j) = b(j)/a(j,j)

           if ( j /= 1 ) then
              t = -b(j)
              call saxpy ( j-1, t, a(1,j), 1, b, 1 )
           end if

        end do

     end if

     if ( cr .or. cab ) then
!
!  Compute RSD or AB as required.
!
        do jj = 1, ju

           j = ju - jj + 1

           if ( qraux(j) /= 0.0E+00 ) then

              temp = a(j,j)
              a(j,j) = qraux(j)

              if ( cr ) then
                 t = - sdot ( n-j+1, a(j,j), 1, rsd(j), 1 ) / a(j,j)
                 call saxpy ( n-j+1, t, a(j,j), 1, rsd(j), 1 )
              end if

              if ( cab ) then
                 t = - sdot ( n-j+1, a(j,j), 1, ab(j), 1 ) / a(j,j)
                 call saxpy ( n-j+1, t, a(j,j), 1, ab(j), 1 )
              end if

              a(j,j) = temp

           end if

        end do

  end if

  return
end
subroutine sscal ( n, sa, x, incx )

!*****************************************************************************80
!
!! SSCAL scales a vector by a constant.
!
!  Modified:
!
!    08 April 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 4 ) SA, the multiplier.
!
!    Input/output, real ( kind = 4 ) X(*), the vector to be scaled.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ) sa
  real ( kind = 4 ) x(*)

  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    m = mod ( n, 5 )

    do i = 1, m
      x(i) = sa * x(i)
    end do

    do i = m+1, n, 5
      x(i)   = sa * x(i)
      x(i+1) = sa * x(i+1)
      x(i+2) = sa * x(i+2)
      x(i+3) = sa * x(i+3)
      x(i+4) = sa * x(i+4)
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      x(ix) = sa * x(ix)
      ix = ix + incx
    end do

  end if

  return
end
subroutine sswap ( n, x, incx, y, incy )

!*****************************************************************************80
!
!! SSWAP interchanges two vectors.
!
!  Modified:
!
!    08 April 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input/output, real ( kind = 4 ) X(*), one of the vectors to swap.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Input/output, real ( kind = 4 ) Y(*), one of the vectors to swap.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive
!    elements of Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ) stemp
  real ( kind = 4 ) x(*)
  real ( kind = 4 ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    m = mod ( n, 3 )

    do i = 1, m
      stemp = x(i)
      x(i) = y(i)
      y(i) = stemp
    end do

    do i = m+1, n, 3

      stemp = x(i)
      x(i) = y(i)
      y(i) = stemp

      stemp = x(i + 1)
      x(i + 1) = y(i + 1)
      y(i + 1) = stemp

      stemp = x(i + 2)
      x(i + 2) = y(i + 2)
      y(i + 2) = stemp

    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      stemp = x(ix)
      x(ix) = y(iy)
      y(iy) = stemp
      ix = ix + incx
      iy = iy + incy
    end do

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
