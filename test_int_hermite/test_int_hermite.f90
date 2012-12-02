subroutine hermite_compute ( order, xtab, weight )

!*****************************************************************************80
!
!! HERMITE_COMPUTE computes a Gauss-Hermite quadrature rule.
!
!  Discussion:
!
!    The abscissas are the zeros of the N-th order Hermite polynomial.
!
!    The integration interval is ( -oo, +oo ).
!
!    The weight function is w(x) = exp ( - x * x ).
!
!    The integral to approximate:
!
!      integral ( -oo < X < +oo ) exp ( - X * X ) * F(X) dX
!
!    The quadrature rule:
!
!      sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
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
!    Input, integer ( kind = 4 ) ORDER, the order of the formula.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) cc
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) s
  real ( kind = 8 ) temp
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(order)

  cc = 1.7724538509D+00 * r8_gamma ( real ( order, kind = 8 ) ) &
    / ( 2.0D+00**( order - 1 ) )

  s = ( 2.0D+00 * real ( order, kind = 8 ) + 1.0D+00 )**( 1.0D+00 / 6.0D+00 )

  do i = 1, ( order + 1 ) / 2

    if ( i == 1 ) then

      x = s**3 - 1.85575D+00 / s

    else if ( i == 2 ) then

      x = x - 1.14D+00 * ( ( real ( order, kind = 8 ) )**0.426D+00 ) / x

    else if ( i == 3 ) then

      x = 1.86D+00 * x - 0.86D+00 * xtab(1)

    else if ( i == 4 ) then

      x = 1.91D+00 * x - 0.91D+00 * xtab(2)

    else

      x = 2.0D+00 * x - xtab(i-2)

    end if

    call hermite_root ( x, order, dp2, p1 )

    xtab(i) = x
    weight(i) = ( cc / dp2 ) / p1

    xtab(order-i+1) = - x
    weight(order-i+1) = weight(i)

  end do
!
!  Reverse the order.
!
  do i = 1, order / 2
    temp            = xtab(i)
    xtab(i)         = xtab(order+1-i)
    xtab(order+1-i) = temp
  end do

  return
end
subroutine hermite_integral ( n, value )

!*****************************************************************************80
!
!! HERMITE_INTEGRAL returns the value of a Hermite polynomial integral.
!
!  Discussion:
!
!    H(n) = Integral ( -oo < x < +oo ) x^n exp(-x^2) dx
!
!    H(n) is 0 for n odd.
!
!    H(n) = (n-1)!! * sqrt(pi) / 2^(n/2) for n even.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the integral.  
!    0 <= N.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) i4_factorial2
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi_sqrt = 1.7724538509055160273D+00
  real ( kind = 8 ) value

  if ( n < 0 ) then

    value = - huge ( value )

  else if ( mod ( n, 2 ) == 1 ) then

    value = 0.0D+00

  else

    value = real ( i4_factorial2 ( n - 1 ), kind = 8 ) * pi_sqrt &
      / 2.0D+00**( n / 2 )

  end if

  return
end
subroutine hermite_recur ( p2, dp2, p1, x, order )

!*****************************************************************************80
!
!! HERMITE_RECUR finds the value and derivative of a Hermite polynomial.
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
!    Output, real ( kind = 8 ) P2, the value of H(ORDER)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of H'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of H(ORDER-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the polynomial 
!    to be computed.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) order
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x
  dp2 = 1.0D+00

  do i = 2, order

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2  = x * p1 - 0.5D+00 * ( real ( i, kind = 8 ) - 1.0D+00 ) * p0
    dp2 = x * dp1 + p1 - 0.5D+00 * ( real ( i, kind = 8 ) - 1.0D+00 ) * dp0

  end do

  return
end
subroutine hermite_root ( x, order, dp2, p1 )

!*****************************************************************************80
!
!! HERMITE_ROOT improves an approximate root of a Hermite polynomial.
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
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the Hermite polynomial.
!
!    Output, real ( kind = 8 ) DP2, the value of H'(ORDER)(X).
!
!    Output, real ( kind = 8 ) P1, the value of H(ORDER-1)(X).
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ), parameter :: eps = 1.0D-12
  integer ( kind = 4 ) order
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  do step = 1, step_max

    call hermite_recur ( p2, dp2, p1, x, order )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
function i4_factorial2 ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    N!!
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the double factorial 
!    function.  If N is less than 1, I4_FACTORIAL2 is returned as 1.
!
!    Output, integer ( kind = 4 ) I4_FACTORIAL2, the value of N!!.
!
  implicit none

  integer ( kind = 4 ) i4_factorial2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_copy

  if ( n < 1 ) then
    i4_factorial2 = 1
    return
  end if

  n_copy = n
  i4_factorial2 = 1

  do while ( 1 < n_copy )
    i4_factorial2 = i4_factorial2 * n_copy
    n_copy = n_copy - 2
  end do

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
!    31 July 2007
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
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_EXACT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_fun ( problem, option, n, x, f )

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
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the problem.
!
!    Input, integer ( kind = 4 ) OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
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
  integer ( kind = 4 ) option
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)

  if ( problem == 1 ) then
    call p01_fun ( option, n, x, f )
  else if ( problem == 2 ) then
    call p02_fun ( option, n, x, f )
  else if ( problem == 3 ) then
    call p03_fun ( option, n, x, f )
  else if ( problem == 4 ) then
    call p04_fun ( option, n, x, f )
  else if ( problem == 5 ) then
    call p05_fun ( option, n, x, f )
  else if ( problem == 6 ) then
    call p06_fun ( option, n, x, f )
  else if ( problem == 7 ) then
    call p07_fun ( option, n, x, f )
  else if ( problem == 8 ) then
    call p08_fun ( option, n, x, f )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FUN - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_gauss_hermite ( problem, order, result )

!*****************************************************************************80
!
!! P00_GAUSS_HERMITE applies a Gauss-Hermite quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2009
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

  real ( kind = 8 ), allocatable, dimension ( : ) :: f_vec
  integer ( kind = 4 ) option
  integer ( kind = 4 ) order
  integer ( kind = 4 ) problem
  real ( kind = 8 ) result
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab

  allocate ( f_vec(1:order) )
  allocate ( weight(1:order) )
  allocate ( xtab(1:order) )

  call hermite_compute ( order, xtab, weight )

  option = 1
  call p00_fun ( problem, option, order, xtab, f_vec )

  result = dot_product ( weight(1:order), f_vec(1:order) )

  deallocate ( f_vec )
  deallocate ( weight )
  deallocate ( xtab )

  return
end
subroutine p00_monte_carlo ( problem, order, result )

!*****************************************************************************80
!
!! P00_MONTE_CARLO applies a Monte Carlo procedure to Hermite integrals.
!
!  Discussion:
!
!    We wish to estimate the integral:
!
!      I(f) = integral ( -oo < x < +oo ) f(x) exp ( - x * x ) dx
!
!    We do this by a Monte Carlo sampling procedure, in which 
!    we select N points X(1:N) from a standard normal distribution,
!    and estimate
!
!      Q(f) = sum ( 1 <= I <= N ) f(x(i)) / sqrt ( pi )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2010
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

  real ( kind = 8 ), allocatable, dimension ( : ) :: f_vec
  integer ( kind = 4 ) option
  integer ( kind = 4 ) order
  real ( kind = 8 ), parameter :: pi_sqrt = 1.7724538509055160273D+00
  integer ( kind = 4 ) problem
  real ( kind = 8 ) result
  integer ( kind = 4 ) seed
  real ( kind = 8 ) weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: x_vec

  allocate ( f_vec(1:order) )
  allocate ( x_vec(1:order) )

  seed = 123456789
  call r8vec_normal_01 ( order, seed, x_vec )

  option = 2
  call p00_fun ( problem, option, order, x_vec, f_vec )

  weight = real ( order, kind = 8 ) / pi_sqrt / sqrt ( 2.0D+00 )

  result = sum ( f_vec(1:order) ) / weight

  deallocate ( f_vec )
  deallocate ( x_vec )

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
!    02 February 2010
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

  problem_num = 8

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
!    31 July 2007
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
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_turing ( problem, h, tol, n, result )

!*****************************************************************************80
!
!! P00_TURING applies the Turing quadrature rule.
!
!  Discussion:
!
!    We consider the approximation:
!
!      Integral ( -oo < x < +oo ) f(x) dx
!
!      = h * Sum ( -oo < i < +oo ) f(nh) + error term
!
!    Given H and a tolerance TOL, we start summing at I = 0, and
!    adding one more term in the positive and negative I directions,
!    until the absolute value of the next two terms being added 
!    is less than TOL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Turing,
!    A Method for the Calculation of the Zeta Function,
!    Proceedings of the London Mathematical Society,
!    Volume 48, 1943, pages 180-197.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the problem.
!
!    Input, real ( kind = 8 ) H, the spacing to use.
!
!    Input, real ( kind = 8 ) TOL, the tolerance.  
!
!    Output, integer ( kind = 4 ) N, the number of pairs of steps taken.
!    The actual number of function evaluations is 2*N+1.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  real ( kind = 8 ) f_vec(2)
  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: n_too_many = 100000
  integer ( kind = 4 ) option
  integer ( kind = 4 ) order
  integer ( kind = 4 ) problem
  real ( kind = 8 ) result
  real ( kind = 8 ) tol
  real ( kind = 8 ) xtab(2)

  option = 0
  n = 0

  result = 0.0D+00
  order = 1
  xtab(1) = 0.0D+00
  call p00_fun ( problem, option, order, xtab, f_vec )
  result = result + h * f_vec(1)

  do

    n = n + 1

    xtab(1) =   real ( n, kind = 8 ) * h
    xtab(2) = - real ( n, kind = 8 ) * h

    order = 2
    call p00_fun ( problem, option, order, xtab, f_vec )

    result = result + h * ( f_vec(1) + f_vec(2) )
!
!  Just do a simple-minded absolute error tolerance check to start with.
!
    if ( abs ( f_vec(1) ) + abs ( f_vec(2) ) <= tol ) then
      exit
    end if
!
!  Just in case things go crazy.
!
    if ( n_too_many <= n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P00_TURING - Warning!'
      write ( *, '(a,i8)' ) &
        '  Number of steps exceeded N_TOO_MANY = ', n_too_many
      exit
    end if

  end do

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
!    31 July 2007
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
  real ( kind = 8 ), parameter :: omega = 1.0D+00
  real ( kind = 8 ), parameter :: pi_sqrt = 1.7724538509055160273D+00

  exact = pi_sqrt * exp ( - omega * omega )

  return
end
subroutine p01_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P01_FUN evaluates the integrand for problem 1.
!
!  Discussion:
!
!    Squire gives exact value as sqrt(pi) * exp(-w*w).
!
!    Integral ( -oo < x < +oo ) exp(-x*x) cos(2*w*x) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Squire,
!    Comparison of Gauss-Hermite and Midpoint Quadrature with Application
!    to the Voigt Function,
!    in Numerical Integration: 
!    Recent Developments, Software and Applications,
!    edited by Patrick Keast, Graeme Fairweather,
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
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
  real ( kind = 8 ), parameter :: omega = 1.0D+00
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(n)

  f(1:n) = cos ( 2.0D+00 * omega * x(1:n) )

  if ( option == 0 ) then
    f(1:n) = f(1:n) * exp ( - x(1:n) * x(1:n) )
  else if ( option == 1 ) then

  else if ( option == 2 ) then
    f(1:n) = f(1:n) * exp ( - 0.5D+00 * x(1:n) * x(1:n) )
  end if

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
!    31 July 2007
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

  title = 'cos(2*omega*x)'

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
!    31 July 2007
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
  real ( kind = 8 ), parameter :: pi_sqrt = 1.7724538509055160273D+00

  exact = pi_sqrt

  return
end
subroutine p02_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P02_FUN evaluates the integrand for problem 2.
!
!  Discussion:
!
!    The exact value is sqrt(pi).
!
!    Integral ( -oo < x < +oo ) exp(-x*x) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
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
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(n)

  f(1:n) = 1.0D+00

  if ( option == 0 ) then
    f(1:n) = f(1:n) * exp ( - x(1:n) * x(1:n) ) 
  else if ( option == 1 ) then

  else if ( option == 2 ) then
    f(1:n) = f(1:n) * exp ( - 0.5D+00 * x(1:n) * x(1:n) )
  end if

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
!    31 July 2007
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

  title = 'exp(-x*x)'

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
!    31 July 2007
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
  real ( kind = 8 ), parameter :: p = 1.0D+00
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932385D+00
  real ( kind = 8 ), parameter :: q = 3.0D+00

  exact = pi / ( q * sin ( pi * p / q ) )

  return
end
subroutine p03_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P03_FUN evaluates the integrand for problem 3.
!
!  Discussion:
!
!    The exact value is pi / (q*sin(pi*p/q) ), assuming 0 < p < q.
!
!    Integral ( -oo < x < +oo ) exp(-px) / ( 1 + exp ( -qx) ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
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
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: p = 1.0D+00
  real ( kind = 8 ), parameter :: q = 3.0D+00
  real ( kind = 8 ) x(n)

  f(1:n) = exp ( - p * x(1:n) ) / ( 1.0D+00 + exp ( -q * x(1:n) ) )

  if ( option == 0 ) then

  else if ( option == 1 ) then
    f(1:n) = f(1:n) * exp ( + x(1:n) * x(1:n) )
  else if ( option == 2 ) then
    f(1:n) = f(1:n) * exp ( + 0.5D+00 * x(1:n) * x(1:n) )
  end if

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
!    31 July 2007
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

  title = 'exp(-px) / ( 1 + exp(-qx) )'

  return
end
subroutine p04_exact ( exact )

!*****************************************************************************80
!
!! P04_EXACT returns the exact integral for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2007
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
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932385D+00

  exact = sqrt ( pi / 2.0D+00 )

  return
end
subroutine p04_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P04_FUN evaluates the integrand for problem 4.
!
!  Discussion:
!
!    The exact value is sqrt ( pi / 2 )
!
!    Integral ( -oo < x < +oo ) sin ( x*x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
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
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(n)

  f(1:n) = sin ( x(1:n)**2 )

  if ( option == 0 ) then

  else if ( option == 1 ) then
    f(1:n) = f(1:n) * exp ( + x(1:n) * x(1:n) )
  else if ( option == 2 ) then
    f(1:n) = f(1:n) * exp ( + 0.5D+00 * x(1:n) * x(1:n) )
  end if

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
!    31 July 2007
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

  title = 'sin(x^2)'

  return
end
subroutine p05_exact ( exact )

!*****************************************************************************80
!
!! P05_EXACT returns the exact integral for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2007
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
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932385D+00

  exact = pi / 3.0D+00

  return
end
subroutine p05_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P05_FUN evaluates the integrand for problem 5.
!
!  Discussion:
!
!    The exact value is pi / 3.
!
!    Integral ( -oo < x < +oo ) dx / ( (1+x^2) sqrt(4+3x^2) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
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
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(n)

  f(1:n) = 1.0D+00 / ( ( 1.0D+00 + x(1:n)**2 ) &
    * sqrt ( 4.0D+00 + 3.0D+00 * x(1:n)**2 ) )

  if ( option == 0 ) then

  else if ( option == 1 ) then
    f(1:n) = f(1:n) * exp ( x(1:n)**2 )
  else if ( option == 2 ) then
    f(1:n) = f(1:n) * exp ( 0.5D+00 * x(1:n)**2 )
  end if

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
!    31 July 2007
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

  title = '1/( (1+x^2) sqrt(4+3x^2) )'

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
!    19 May 2009
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
  integer ( kind = 4 ) i4_factorial2
  integer ( kind = 4 ) m
  real ( kind = 8 ), parameter :: pi_sqrt = 1.7724538509055160273D+00

  call p06_param ( 'G', 'M', m )

  if ( m <= -1 ) then

    exact = - huge ( exact )

  else if ( mod ( m, 2 ) == 1 ) then

    exact = 0.0D+00

  else

    exact = real ( i4_factorial2 ( m - 1 ), kind = 8 ) * pi_sqrt &
      / 2.0D+00**( m / 2 )

  end if

  return
end
subroutine p06_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P06_FUN evaluates the integrand for problem 6.
!
!  Discussion:
!
!    The exact value is (m-1)!! * sqrt ( pi ) / sqrt ( 2^m ).
!
!    Integral ( -oo < x < +oo ) x^m exp (-x*x) dx
!
!    The parameter M is set by calling P06_PARAM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
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
  integer ( kind = 4 ) m
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(n)

  call p06_param ( 'G', 'M', m )

  f(1:n) = x(1:n)**m

  if ( option == 0 ) then
    f(1:n) = f(1:n) * exp ( - x(1:n)**2 )
  else if ( option == 1 ) then

  else if ( option == 2 ) then
    f(1:n) = f(1:n) * exp ( - 0.5D+00 * x(1:n)**2 )
  end if

  return
end
subroutine p06_param ( action, name, value )

!*****************************************************************************80
!
!! P06_PARAM gets or sets parameters for problem 6.
!
!  Discussion:
!
!    The parameter is named "M", and it represents the value of the exponent
!    in the integrand function:
!
!    Integral ( -oo < x < +oo ) x^m exp (-x*x) dx
!
!    M must be greater than -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ACTION, the action.
!    'S' to set the value,
!    'G' to get the value.
!
!    Input, character NAME, the parameter name.
!    'M', the exponent.
!
!    Input/output, integer ( kind = 4 ) VALUE, the parameter value.
!    If ACTION = 'S', then VALUE is an input quantity, and M is set to VALUE.
!    If ACTION = 'G', then VALUE is an output quantity, and VALUE is set to M.
!
  implicit none

  character action
  character name
  integer ( kind = 4 ), save :: m = 0
  integer ( kind = 4 ) value

  if ( action == 'S' .or. action == 's' ) then

    if ( value <= -1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P06_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Parameter M must be greater than -1.'
      stop
    end if

    m = value

  else if ( action == 'G' .or. action == 'g' ) then

    value = m

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P06_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized value of ACTION = "' // action // '".'
    stop

  end if

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
!    19 May 2009
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

  title = 'x^m exp(-x*x)'

  return
end
subroutine p07_exact ( exact )

!*****************************************************************************80
!
!! P07_EXACT returns the exact integral for problem 7.
!
!  Discussion:
!
!    The 20 digit values of pi^(1/2) and e^(1/4) were computed by Mathematica.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2010
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

  real ( kind = 8 ), parameter :: e_sqrt_sqrt = 1.2840254166877414841D+00
  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi_sqrt = 1.7724538509055160273D+00

  exact = 0.25D+00 * pi_sqrt / e_sqrt_sqrt

  return
end
subroutine p07_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P07_FUN evaluates the integrand for problem 7.
!
!  Discussion:
!
!    The exact value is (1/4) sqrt(pi) / sqrt(sqrt(e)).
!
!    Integral ( -oo < x < +oo ) x^2 cos(x) e^(-x^2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Prem Kythe, Michael Schaeferkotter,
!    Handbook of Computational Methods for Integration,
!    Chapman and Hall, 2004,
!    ISBN: 1-58488-428-2,
!    LC: QA299.3.K98.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
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
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(n)

  if ( option == 0 ) then
    f(1:n) = x(1:n)**2 * cos ( x(1:n) ) * exp ( - x(1:n)**2 )
  else if ( option == 1 ) then
    f(1:n) = x(1:n)**2 * cos ( x(1:n) )
  else if ( option == 2 ) then
    f(1:n) = x(1:n)**2 * cos ( x(1:n) ) * exp ( - x(1:n)**2 / 2.0D+00 )
  end if

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
!    02 February 2010
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

  title = 'x^2 cos ( x ) exp(-x*x)'

  return
end
subroutine p08_exact ( exact )

!*****************************************************************************80
!
!! P08_EXACT returns the exact integral for problem 8.
!
!  Discussion:
!
!    The 20 digit value of the answer was computed by Mathematica.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 2010
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

  exact = 3.0088235661136433510D+00

  return
end
subroutine p08_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P08_FUN evaluates the integrand for problem 8.
!
!  Discussion:
!
!    The exact value is sqrt ( 2 pi ) * HypergeometricU ( -1/2, 0, 1 ).
!
!    Integral ( -oo < x < +oo ) sqrt(1+x*x/2) * exp(-x*x/2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Prem Kythe, Michael Schaeferkotter,
!    Handbook of Computational Methods for Integration,
!    Chapman and Hall, 2004,
!    ISBN: 1-58488-428-2,
!    LC: QA299.3.K98.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
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
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(n)

  f(1:n) = sqrt ( 1.0D+00 + 0.5D+00 * x(1:n)**2 )

  if ( option == 0 ) then
    f(1:n) = f(1:n) * exp ( - 0.5D+00 * x(1:n)**2 )
  else if ( option == 1 ) then
    f(1:n) = f(1:n) * exp ( + 0.5D+00 * x(1:n)**2 )
  else if ( option == 2 ) then

  end if

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
!    30 July 2010
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

  title = 'sqrt(1+x*x/2) * exp(-x*x/2)'

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
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
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

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

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
      res = 1.0D+00 / y
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
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
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
        y = y + 1.0D+00
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
      sum = sum + ( y - 0.5D+00 ) * log ( y )
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

  if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
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
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED
!    in order to set the seed of the random number generator.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is 
!    negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer ( kind = 4 ) SAVED, is 0 or 1 depending on whether there 
!    is a single saved value left over from the previous call.
!
!    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range 
!    of entries of X that we need to compute.  This starts off as 1:N, but 
!    is adjusted if we have a saved value that can be immediately stored
!    in X(1), and so on.
!
!    Local, real ( kind = 8 ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) m
  integer ( kind = 4 ), save :: made = 0
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r(n+1)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ), save :: saved = 0
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: two = 2
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
  real ( kind = 8 ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, two ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
