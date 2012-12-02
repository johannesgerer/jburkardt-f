subroutine p00_f ( problem, x, f )

!*****************************************************************************80
!
!! P00_F evaluates the function for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, real ( kind = 8 ) X, the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  real ( kind = 8 ) f
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x

  if ( problem == 1 ) then
    call p01_f ( x, f )
  else if ( problem == 2 ) then
    call p02_f ( x, f )
  else if ( problem == 3 ) then
    call p03_f ( x, f )
  else if ( problem == 4 ) then
    call p04_f ( x, f )
  else if ( problem == 5 ) then
    call p05_f ( x, f )
  else if ( problem == 6 ) then
    call p06_f ( x, f )
  else if ( problem == 7 ) then
    call p07_f ( x, f )
  else if ( problem == 8 ) then
    call p08_f ( x, f )
  else if ( problem == 9 ) then
    call p09_f ( x, f )
  else if ( problem == 10 ) then
    call p10_f ( x, f )
  else if ( problem == 11 ) then
    call p11_f ( x, f )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_F - Fatal error!'
    write ( *, '(a,i12)' ) '  Illegal problem number PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_f1 ( problem, x, f1 )

!*****************************************************************************80
!
!! P00_F1 evaluates the first derivative for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F1, the first derivative of the 
!    objective function.
!
  implicit none

  real ( kind = 8 ) f1
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x

  if ( problem == 1 ) then
    call p01_f1 ( x, f1 )
  else if ( problem == 2 ) then
    call p02_f1 ( x, f1 )
  else if ( problem == 3 ) then
    call p03_f1 ( x, f1 )
  else if ( problem == 4 ) then
    call p04_f1 ( x, f1 )
  else if ( problem == 5 ) then
    call p05_f1 ( x, f1 )
  else if ( problem == 6 ) then
    call p06_f1 ( x, f1 )
  else if ( problem == 7 ) then
    call p07_f1 ( x, f1 )
  else if ( problem == 8 ) then
    call p08_f1 ( x, f1 )
  else if ( problem == 9 ) then
    call p09_f1 ( x, f1 )
  else if ( problem == 10 ) then
    call p10_f1 ( x, f1 )
  else if ( problem == 11 ) then
    call p11_f1 ( x, f1 )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_F1 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_f1_dif ( problem, x, f1_dif )

!*****************************************************************************80
!
!! P00_F1_DIF approximates the first derivative via finite differences.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, real ( kind = 8 ) X, the point where the gradient is to 
!    be approximated.
!
!    Output, real ( kind = 8 ) F1_DIF, the approximated gradient vector.
!
  implicit none

  real ( kind = 8 ) dx
  real ( kind = 8 ) eps
  real ( kind = 8 ) f1_dif
  real ( kind = 8 ) fminus
  real ( kind = 8 ) fplus
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x
  real ( kind = 8 ) xi

  eps = ( epsilon ( eps ) )**0.33D+00

  if ( 0.0D+00 <= x ) then
    dx = eps * ( x + 1.0D+00 )
  else
    dx = eps * ( x - 1.0D+00 )
  end if

  xi = x
  x = xi + dx
  call p00_f ( problem, x, fplus )

  x = xi - dx
  call p00_f ( problem, x, fminus )

  f1_dif = ( fplus - fminus ) / ( 2.0D+00 * dx )

  x = xi

  return
end
subroutine p00_f2 ( problem, x, f2 )

!*****************************************************************************80
!
!! P00_F2 evaluates the second derivative for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the values of the variables.
!
!    Output, real ( kind = 8 ) F2, the second derivative.
!
  implicit none

  real ( kind = 8 ) f2
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x

  if ( problem == 1 ) then
    call p01_f2 ( x, f2 )
  else if ( problem == 2 ) then
    call p02_f2 ( x, f2 )
  else if ( problem == 3 ) then
    call p03_f2 ( x, f2 )
  else if ( problem == 4 ) then
    call p04_f2 ( x, f2 )
  else if ( problem == 5 ) then
    call p05_f2 ( x, f2 )
  else if ( problem == 6 ) then
    call p06_f2 ( x, f2 )
  else if ( problem == 7 ) then
    call p07_f2 ( x, f2 )
  else if ( problem == 8 ) then
    call p08_f2 ( x, f2 )
  else if ( problem == 9 ) then
    call p09_f2 ( x, f2 )
  else if ( problem == 10 ) then
    call p10_f2 ( x, f2 )
  else if ( problem == 11 ) then
    call p11_f2 ( x, f2 )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_F2 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_f2_dif ( problem, x, f2_dif )

!*****************************************************************************80
!
!! P00_F2_DIF approximates the second derivative via finite differences.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F2_DIF, the approximate second derivative.
!
  implicit none

  real ( kind = 8 ) eps
  real ( kind = 8 ) f00
  real ( kind = 8 ) f2_dif
  real ( kind = 8 ) fmm
  real ( kind = 8 ) fpp
  integer ( kind = 4 ) problem
  real ( kind = 8 ) s
  real ( kind = 8 ) x
  real ( kind = 8 ) xi
!
!  Choose the stepsize.
!
  eps = ( epsilon ( eps ) )**0.33D+00

  s = eps * ( abs ( x ) + 1.0D+00 )

  xi = x

  call p00_f ( problem, x, f00 )

  x = xi + s
  call p00_f ( problem, x, fpp )

  x = xi - s
  call p00_f ( problem, x, fmm )

  f2_dif = ( ( fpp - f00 ) + ( fmm - f00 ) ) / s / s

  x = xi

  return
end
function p00_fmin ( a, b, problem, tol )

!*****************************************************************************80
!
!! P00_FMIN seeks a minimizer of a scalar function of a scalar variable.
!
!  Discussion:
!
!    FMIN seeks an approximation to the point where F attains a minimum on
!    the interval (A,B).
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much
!    slower than that for a Fibonacci search.  If F has a continuous
!    second derivative which is positive at the minimum (which is not
!    at A or B), then convergence is superlinear, and usually of the
!    order of about 1.324....
!
!    The function F is never evaluated at two points closer together
!    than EPS * ABS ( FMIN ) + (TOL/3), where EPS is approximately the
!    square root of the relative machine precision.  If F is a unimodal
!    function and the computed values of F are always unimodal when
!    separated by at least EPS * ABS ( XSTAR ) + (TOL/3), then FMIN
!    approximates the abcissa of the global minimum of F on the
!    interval [A, B] with an error less than 3 * EPS * ABS ( FMIN ) + TOL.
!    If F is not unimodal, then FMIN may approximate a local, but
!    perhaps non-global, minimum to the same accuracy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization without Derivatives,
!    Prentice Hall, 1973.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters
!
!    Input/output, real ( kind = 8 ) A, B.  On input, the left and right
!    endpoints of the initial interval.  On output, the lower and upper 
!    bounds for the minimizer.
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of a problem.
!
!    Input, real ( kind = 8 ) TOL, the desired length of the interval of
!    uncertainty of the final result.  TOL must not be negative.
!
!    Output, real ( kind = 8 ) P00_FMIN, the abcissa approximating the 
!    minimizer of f.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) eps
  real ( kind = 8 ) fu
  real ( kind = 8 ) fv
  real ( kind = 8 ) fw
  real ( kind = 8 ) fx
  real ( kind = 8 ) midpoint
  integer ( kind = 4 ) problem
  real ( kind = 8 ) p
  real ( kind = 8 ) p00_fmin
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) tol
  real ( kind = 8 ) tol1
  real ( kind = 8 ) tol2
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) x

  c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )
!
!  C is the squared inverse of the golden ratio.
!
!  EPS is the square root of the relative machine precision.
!
  eps = sqrt ( epsilon ( eps ) )
!
!  Initialization.
!
  v = a + c * ( b - a )
  w = v
  x = v
  e = 0.0D+00
  call p00_f ( problem, x, fx )
  fv = fx
  fw = fx
!
!  The main loop starts here.
!
  do

    midpoint = 0.5D+00 * ( a + b )
    tol1 = eps * abs ( x ) + tol / 3.0D+00
    tol2 = 2.0D+00 * tol1
!
!  Check the stopping criterion.
!
    if ( abs ( x - midpoint ) <= ( tol2 - 0.5D+00 * ( b - a ) ) ) then
      exit
    end if
!
!  Is golden-section necessary?
!
    if ( abs ( e ) <= tol1 ) then
      if ( midpoint <= x ) then
        e = a - x
      else
        e = b - x
      end if

      d = c * e
!
!  Consider fitting a parabola.
!
    else

      r = ( x - w ) * ( fx - fv )
      q = ( x - v ) * ( fx - fw )
      p = ( x - v ) * q - ( x - w ) * r
      q = 2.0D+00 * ( q - r )
      if ( 0.0D+00 < q ) then
        p = -p
      end if
      q = abs ( q )
      r = e
      e = d
!
!  Choose a golden-section step if the parabola is not advised.
!
      if ( &
        ( abs ( 0.5D+00 * q * r ) <= abs ( p ) ) .or. &
        ( p <= q * ( a - x ) ) .or. &
        ( q * ( b - x ) <= p ) ) then

        if ( midpoint <= x ) then
          e = a - x
        else
          e = b - x
        end if

        d = c * e
!
!  Choose a parabolic interpolation step.
!
      else

        d = p / q
        u = x + d

        if ( ( u - a ) < tol2 ) then
          d = sign ( tol1, midpoint - x )
        end if

        if ( ( b - u ) < tol2 ) then
          d = sign ( tol1, midpoint - x )
        end if

     end if

   end if
!
!  F must not be evaluated too close to X.
!
    if ( tol1 <= abs ( d ) ) then
      u = x + d
    end if

    if ( abs ( d ) < tol1 ) then
      u = x + sign ( tol1, d )
    end if

    call p00_f ( problem, u, fu )
!
!  Update the data.
!
    if ( fu <= fx ) then

      if ( x <= u ) then
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
      cycle

    end if

    if ( u < x ) then
      a = u
    else
      b = u
    end if

    if ( fu <= fw .or. w == x ) then
      v = w
      fv = fw
      w = u
      fw = fu
    else if ( fu <= fv .or. v == x .or. v == w ) then
      v = u
      fv = fu
    end if

  end do

  p00_fmin = x

  return
end
subroutine p00_interval ( problem, a, b )

!*****************************************************************************80
!
!! P00_INTERVAL returns a bracketing interval for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Output, real ( kind = 8 ) A, B, two points, between which a local
!    minimizer should be sought.
!
  implicit none

  integer ( kind = 4 ) problem
  real ( kind = 8 ) a
  real ( kind = 8 ) b

  if ( problem == 1 ) then
    call p01_interval ( a, b )
  else if ( problem == 2 ) then
    call p02_interval ( a, b )
  else if ( problem == 3 ) then
    call p03_interval ( a, b )
  else if ( problem == 4 ) then
    call p04_interval ( a, b )
  else if ( problem == 5 ) then
    call p05_interval ( a, b )
  else if ( problem == 6 ) then
    call p06_interval ( a, b )
  else if ( problem == 7 ) then
    call p07_interval ( a, b )
  else if ( problem == 8 ) then
    call p08_interval ( a, b )
  else if ( problem == 9 ) then
    call p09_interval ( a, b )
  else if ( problem == 10 ) then
    call p10_interval ( a, b )
  else if ( problem == 11 ) then
    call p11_interval ( a, b )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_INTERVAL - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_problem_num ( problem_num )

!*****************************************************************************80
!
!! P00_PROBLEM_NUM returns the number of problems available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!   Output, integer ( kind = 4 ) PROBLEM_NUM, the number of problems.
!
  implicit none

  integer ( kind = 4 ) problem_num

  problem_num = 11

  return
end
subroutine p00_sol ( problem, know, x )

!*****************************************************************************80
!
!! P00_SOL returns the solution for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Output, integer ( kind = 4 ) KNOW.
!    If KNOW is 0, then the solution is not known.
!    If KNOW is positive, then the solution is known, and is returned in X.
!
!    Output, real ( kind = 8 ) X, the solution, if known.
!
  implicit none

  integer ( kind = 4 ) know
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x

  if ( problem == 1 ) then
    call p01_sol ( know, x )
  else if ( problem == 2 ) then
    call p02_sol ( know, x )
  else if ( problem == 3 ) then
    call p03_sol ( know, x )
  else if ( problem == 4 ) then
    call p04_sol ( know, x )
  else if ( problem == 5 ) then
    call p05_sol ( know, x )
  else if ( problem == 6 ) then
    call p06_sol ( know, x )
  else if ( problem == 7 ) then
    call p07_sol ( know, x )
  else if ( problem == 8 ) then
    call p08_sol ( know, x )
  else if ( problem == 9 ) then
    call p09_sol ( know, x )
  else if ( problem == 10 ) then
    call p10_sol ( know, x )
  else if ( problem == 11 ) then
    call p11_sol ( know, x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SOL - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_start ( problem, x )

!*****************************************************************************80
!
!! P00_START returns a starting point for optimization for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Output, real ( kind = 8 ) X, a starting point for the optimization.
!
  implicit none

  integer ( kind = 4 ) problem
  real ( kind = 8 ) x

  if ( problem == 1 ) then
    call p01_start ( x )
  else if ( problem == 2 ) then
    call p02_start ( x )
  else if ( problem == 3 ) then
    call p03_start ( x )
  else if ( problem == 4 ) then
    call p04_start ( x )
  else if ( problem == 5 ) then
    call p05_start ( x )
  else if ( problem == 6 ) then
    call p06_start ( x )
  else if ( problem == 7 ) then
    call p07_start ( x )
  else if ( problem == 8 ) then
    call p08_start ( x )
  else if ( problem == 9 ) then
    call p09_start ( x )
  else if ( problem == 10 ) then
    call p10_start ( x )
  else if ( problem == 11 ) then
    call p11_start ( x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_START - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_title ( problem, title )

!*****************************************************************************80
!
!! P00_TITLE returns a title for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Output, character ( len = * ) TITLE, a title for the problem.
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
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p01_f ( x, f )

!*****************************************************************************80
!
!! P01_F evaluates the objective function for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) x

  f = ( x - 2.0D+00 ) * ( x - 2.0D+00 ) + 1.0D+00

  return
end
subroutine p01_f1 ( x, f1 )

!*****************************************************************************80
!
!! P01_F1 evaluates the first derivative for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F1, the first derivative of the
!    objective function.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  f1 = 2.0D+00 * ( x - 2.0D+00 )

  return
end
subroutine p01_f2 ( x, f2 )

!*****************************************************************************80
!
!! P01_F2 evaluates the second derivative for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the values of the variables.
!
!    Output, real ( kind = 8 ) F2, the second derivative.
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) x

  f2 = 2.0D+00

  return
end
subroutine p01_interval ( a, b )

!*****************************************************************************80
!
!! P01_INTERVAL returns a starting interval for optimization for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, two points defining an interval in which
!    the local minimizer should be sought.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 3.141592653589793D+00

  return
end
subroutine p01_sol ( know, x )

!*****************************************************************************80
!
!! P01_SOL returns the solution for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KNOW.
!    If KNOW is 0, then the solution is not known.
!    If KNOW is positive, then the solution is known, and is returned in X.
!
!    Output, real ( kind = 8 ) X, the solution, if known.
!
  implicit none

  integer ( kind = 4 ) know
  real ( kind = 8 ) x

  know = 1

  x = 2.0D+00

  return
end
subroutine p01_start ( x )

!*****************************************************************************80
!
!! P01_START returns a starting point for optimization for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X, a starting point for the optimization.
!
  implicit none

  real ( kind = 8 ) x

  x = 3.141592653589793D+00

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns a title for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Simple quadratic, (x-2)^2+1.'

  return
end
subroutine p02_f ( x, f )

!*****************************************************************************80
!
!! P02_F evaluates the objective function for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    LE Scales,
!    Introduction to Non-Linear Optimization,
!    Springer, 1985.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) x

  f = x * x + exp ( - x )

  return
end
subroutine p02_f1 ( x, f1 )

!*****************************************************************************80
!
!! P02_F1 evaluates the first derivative for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F1, the first derivative of the 
!    objective function.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  f1 = 2.0D+00 * x - exp ( -x )

  return
end
subroutine p02_f2 ( x, f2 )

!*****************************************************************************80
!
!! P02_F2 evaluates the second derivative for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    LE Scales,
!    Introduction to Non-Linear Optimization,
!    Springer, 1985.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the values of the variables.
!
!    Output, real ( kind = 8 ) F2, the second derivative.
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) x

  f2 = 2.0D+00 + exp ( -x )

  return
end
subroutine p02_interval ( a, b )

!*****************************************************************************80
!
!! P02_INTERVAL returns a starting interval for optimization for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, two points defining an interval in which
!    the local minimizer should be sought.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a =  0.0D+00
  b =  1.0D+00

  return
end
subroutine p02_sol ( know, x )

!*****************************************************************************80
!
!! P02_SOL returns the solution for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KNOW.
!    If KNOW is 0, then the solution is not known.
!    If KNOW is positive, then the solution is known, and is returned in X.
!
!    Output, real ( kind = 8 ) X, the solution, if known.
!
  implicit none

  integer ( kind = 4 ) know
  real ( kind = 8 ) x

  know = 1

  x = 0.351734D+00

  return
end
subroutine p02_start ( x )

!*****************************************************************************80
!
!! P02_START returns a starting point for optimization for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X, a starting point for the optimization.
!
  implicit none

  real ( kind = 8 ) x

  x = 0.8D+00

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns a title for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Quadratic plus exponential, x^2 + e^(-x).'

  return
end
subroutine p03_f ( x, f )

!*****************************************************************************80
!
!! P03_F evaluates the objective function for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    LE Scales,
!    Introduction to Non-Linear Optimization,
!    Springer, 1985.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) x

  f = ( ( x * x + 2.0D+00 ) * x + 1.0D+00 ) * x + 3.0D+00

  return
end
subroutine p03_f1 ( x, f1 )

!*****************************************************************************80
!
!! P03_F1 evaluates the first derivative for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F1, the first derivative of the 
!    objective function.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  f1 = ( 4.0D+00 * x * x + 4.0D+00 ) * x + 1.0D+00

  return
end
subroutine p03_f2 ( x, f2 )

!*****************************************************************************80
!
!! P03_F2 evaluates the second derivative for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    LE Scales,
!    Introduction to Non-Linear Optimization,
!    Springer, 1985.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the values of the variables.
!
!    Output, real ( kind = 8 ) F2, the second derivative.
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) x

  f2 = 12.0D+00 * x * x + 4.0D+00

  return
end
subroutine p03_interval ( a, b )

!*****************************************************************************80
!
!! P03_INTERVAL returns a starting interval for optimization for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, two points defining an interval in which
!    the local minimizer should be sought.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a =  -2.0D+00
  b =  +2.0D+00

  return
end
subroutine p03_sol ( know, x )

!*****************************************************************************80
!
!! P03_SOL returns the solution for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KNOW.
!    If KNOW is 0, then the solution is not known.
!    If KNOW is positive, then the solution is known, and is returned in X.
!
!    Output, real ( kind = 8 ) X, the solution, if known.
!
  implicit none

  integer ( kind = 4 ) know
  real ( kind = 8 ) x

  know = 1

  x = -0.236733D+00

  return
end
subroutine p03_start ( x )

!*****************************************************************************80
!
!! P03_START returns a starting point for optimization for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X, a starting point for the optimization.
!
  implicit none

  real ( kind = 8 ) x

  x = 1.5D+00

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns a title for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Quartic, x^4 + 2x^2 + x + 3.'

  return
end
subroutine p04_f ( x, f )

!*****************************************************************************80
!
!! P04_F evaluates the objective function for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    LE Scales,
!    Introduction to Non-Linear Optimization,
!    Springer, 1985.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) x

  f = exp ( x ) + 0.01D+00 / x

  return
end
subroutine p04_f1 ( x, f1 )

!*****************************************************************************80
!
!! P04_F1 evaluates the first derivative for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F1, the first derivative of the 
!    objective function.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  f1 = exp ( x ) - 0.01D+00 / x / x

  return
end
subroutine p04_f2 ( x, f2 )

!*****************************************************************************80
!
!! P04_F2 evaluates the second derivative for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    LE Scales,
!    Introduction to Non-Linear Optimization,
!    Springer, 1985.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the values of the variables.
!
!    Output, real ( kind = 8 ) F2, the second derivative.
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) x

  f2 = exp ( x ) + 0.02D+00 / x / x / x

  return
end
subroutine p04_interval ( a, b )

!*****************************************************************************80
!
!! P04_INTERVAL returns a starting interval for optimization for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, two points defining an interval in which
!    the local minimizer should be sought.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a =  0.0001D+00
  b =  1.0D+00

  return
end
subroutine p04_sol ( know, x )

!*****************************************************************************80
!
!! P04_SOL returns the solution for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KNOW.
!    If KNOW is 0, then the solution is not known.
!    If KNOW is positive, then the solution is known, and is returned in X.
!
!    Output, real ( kind = 8 ) X, the solution, if known.
!
  implicit none

  integer ( kind = 4 ) know
  real ( kind = 8 ) x

  know = 1

  x = 0.0953446D+00

  return
end
subroutine p04_start ( x )

!*****************************************************************************80
!
!! P04_START returns a starting point for optimization for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X, a starting point for the optimization.
!
  implicit none

  real ( kind = 8 ) x

  x = 0.95D+00

  return
end
subroutine p04_title ( title )

!*****************************************************************************80
!
!! P04_TITLE returns a title for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Steep valley, e^x + 1/(100x).'

  return
end
subroutine p05_f ( x, f )

!*****************************************************************************80
!
!! P05_F evaluates the objective function for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    LE Scales,
!    Introduction to Non-Linear Optimization,
!    Springer, 1985.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) x

  f = exp ( x ) - 2.0D+00 * x + 0.01D+00 / x - 0.000001D+00 / x / x

  return
end
subroutine p05_f1 ( x, f1 )

!*****************************************************************************80
!
!! P05_F1 evaluates the first derivative for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F1, the first derivative of the 
!    objective function.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  f1 = exp ( x ) - 2.0D+00 - 0.01D+00 / x / x + 0.000002D+00 / x / x / x

  return
end
subroutine p05_f2 ( x, f2 )

!*****************************************************************************80
!
!! P05_F2 evaluates the second derivative for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    LE Scales,
!    Introduction to Non-Linear Optimization,
!    Springer, 1985.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the values of the variables.
!
!    Output, real ( kind = 8 ) F2, the second derivative.
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) x

  f2 = exp ( x ) + 0.02D+00 / x / x / x - 0.000006D+00 / x / x / x / x

  return
end
subroutine p05_interval ( a, b )

!*****************************************************************************80
!
!! P05_INTERVAL returns a starting interval for optimization for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, two points defining an interval in which
!    the local minimizer should be sought.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a =  0.0002D+00
  b =  2.0D+00

  return
end
subroutine p05_sol ( know, x )

!*****************************************************************************80
!
!! P05_SOL returns the solution for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KNOW.
!    If KNOW is 0, then the solution is not known.
!    If KNOW is positive, then the solution is known, and is returned in X.
!
!    Output, real ( kind = 8 ) X, the solution, if known.
!
  implicit none

  integer ( kind = 4 ) know
  real ( kind = 8 ) x

  know = 1

  x = 0.703206D+00

  return
end
subroutine p05_start ( x )

!*****************************************************************************80
!
!! P05_START returns a starting point for optimization for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X, a starting point for the optimization.
!
  implicit none

  real ( kind = 8 ) x

  x = 1.5D+00

  return
end
subroutine p05_title ( title )

!*****************************************************************************80
!
!! P05_TITLE returns a title for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Steep valley, e^x - 2x + 1/(100x) - 1/(1000000x^2).'

  return
end
subroutine p06_f ( x, f )

!*****************************************************************************80
!
!! P06_F evaluates the objective function for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Prentice Hall 1973,
!    Reprinted Dover, 2002
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) x

  f = 2.0D+00 - x

  return
end
subroutine p06_f1 ( x, f1 )

!*****************************************************************************80
!
!! P06_F1 evaluates the first derivative for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F1, the first derivative of the 
!    objective function.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  f1 = -1.0D+00

  return
end
subroutine p06_f2 ( x, f2 )

!*****************************************************************************80
!
!! P06_F2 evaluates the second derivative for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    LE Scales,
!    Introduction to Non-Linear Optimization,
!    Springer, 1985.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the values of the variables.
!
!    Output, real ( kind = 8 ) F2, the second derivative.
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) x

  f2 = 0.0D+00

  return
end
subroutine p06_interval ( a, b )

!*****************************************************************************80
!
!! P06_INTERVAL returns a starting interval for optimization for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, two points defining an interval in which
!    the local minimizer should be sought.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a =  7.0D+00
  b =  9.0D+00

  return
end
subroutine p06_sol ( know, x )

!*****************************************************************************80
!
!! P06_SOL returns the solution for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KNOW.
!    If KNOW is 0, then the solution is not known.
!    If KNOW is positive, then the solution is known, and is returned in X.
!
!    Output, real ( kind = 8 ) X, the solution, if known.
!
  implicit none

  integer ( kind = 4 ) know
  real ( kind = 8 ) x

  know = 1

  x = 9.0D+00

  return
end
subroutine p06_start ( x )

!*****************************************************************************80
!
!! P06_START returns a starting point for optimization for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X, a starting point for the optimization.
!
  implicit none

  real ( kind = 8 ) x

  x = 7.2D+00

  return
end
subroutine p06_title ( title )

!*****************************************************************************80
!
!! P06_TITLE returns a title for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'line, 2 - x.'

  return
end
subroutine p07_f ( x, f )

!*****************************************************************************80
!
!! P07_F evaluates the objective function for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Prentice Hall 1973,
!    Reprinted Dover, 2002
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) x

  f = ( x + sin ( x ) ) * exp ( - x * x )

  return
end
subroutine p07_f1 ( x, f1 )

!*****************************************************************************80
!
!! P07_F1 evaluates the first derivative for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F1, the first derivative of the 
!    objective function.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  f1 = ( 1.0D+00 - 2.0D+00 * x * x + cos ( x ) &
         - 2.0D+00 * x * sin ( x ) ) * exp ( - x * x )

  return
end
subroutine p07_f2 ( x, f2 )

!*****************************************************************************80
!
!! P07_F2 evaluates the second derivative for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the values of the variables.
!
!    Output, real ( kind = 8 ) F2, the second derivative.
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) x

  f2 = ( - 4.0D+00 - 2.0D+00 * x - 4.0D+00 * x * x * x &
    - 3.0D+00 * sin ( x ) - 4.0D+00 * x * cos ( x ) & 
    + 4.0D+00 * x * x * sin ( x ) ) * exp ( - x * x )

  return
end
subroutine p07_interval ( a, b )

!*****************************************************************************80
!
!! P07_INTERVAL returns a starting interval for optimization for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, two points defining an interval in which
!    the local minimizer should be sought.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a =  -10.0D+00
  b =  +10.0D+00

  return
end
subroutine p07_sol ( know, x )

!*****************************************************************************80
!
!! P07_SOL returns the solution for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KNOW.
!    If KNOW is 0, then the solution is not known.
!    If KNOW is positive, then the solution is known, and is returned in X.
!
!    Output, real ( kind = 8 ) X, the solution, if known.
!
  implicit none

  integer ( kind = 4 ) know
  real ( kind = 8 ) x

  know = 1

  x = -0.6795786599525D+00

  return
end
subroutine p07_start ( x )

!*****************************************************************************80
!
!! P07_START returns a starting point for optimization for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X, a starting point for the optimization.
!
  implicit none

  real ( kind = 8 ) x

  x = -5.0D+00

  return
end
subroutine p07_title ( title )

!*****************************************************************************80
!
!! P07_TITLE returns a title for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The dying snake, ( x + sin(x) ) * e^(-x^2).'

  return
end
subroutine p08_f ( x, f )

!*****************************************************************************80
!
!! P08_F evaluates the objective function for problem 8.
!
!  Discussion:
!
!    This function looks positive, but has a pole at x = pi,
!    near which f -> negative infinity, and has two zeroes nearby.  
!    None of this will show up computationally.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arnold Krommer, Christoph Ueberhuber,
!    Numerical Integration on Advanced Systems,
!    Springer, 1994, pages 185-186.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( x == pi ) then
    f = - 10000.0D+00
  else
    f = 3.0D+00 * x * x + 1.0D+00 + ( log ( ( x - pi ) * ( x - pi ) ) ) / pi**4
  end if

  return
end
subroutine p08_f1 ( x, f1 )

!*****************************************************************************80
!
!! P08_F1 evaluates the first derivative for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F1, the first derivative of the 
!    objective function.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( x == pi ) then
    f1 = 0.0D+00
  else
    f1 = 6.0D+00 * x + ( 2.0D+00 / ( x - pi ) ) / pi**4
  end if

  return
end
subroutine p08_f2 ( x, f2 )

!*****************************************************************************80
!
!! P08_F2 evaluates the second derivative for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the values of the variables.
!
!    Output, real ( kind = 8 ) F2, the second derivative.
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( x == pi ) then
    f2 = 1.0D+00
  else
    f2 = 6.0D+00 + ( - 2.0D+00 / ( x - pi ) / ( x - pi ) ) / pi**4
  end if

  return
end
subroutine p08_interval ( a, b )

!*****************************************************************************80
!
!! P08_INTERVAL returns a starting interval for optimization for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, two points defining an interval in which
!    the local minimizer should be sought.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a =  2.0D+00
  b =  4.0D+00

  return
end
subroutine p08_sol ( know, x )

!*****************************************************************************80
!
!! P08_SOL returns the solution for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KNOW.
!    If KNOW is 0, then the solution is not known.
!    If KNOW is positive, then the solution is known, and is returned in X.
!
!    Output, real ( kind = 8 ) X, the solution, if known.
!
  implicit none

  integer ( kind = 4 ) know
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  know = 1

  x = pi

  return
end
subroutine p08_start ( x )

!*****************************************************************************80
!
!! P08_START returns a starting point for optimization for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X, a starting point for the optimization.
!
  implicit none

  real ( kind = 8 ) x

  x = 3.1D+00

  return
end
subroutine p08_title ( title )

!*****************************************************************************80
!
!! P08_TITLE returns a title for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The "Thin Pole", x^2+1+log((pi-x)^2)/pi^4'

  return
end
subroutine p09_f ( x, f )

!*****************************************************************************80
!
!! P09_F evaluates the objective function for problem 9.
!
!  Discussion:
!
!    This function is oscillatory, with many local minima.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) x

  f = x * x - 10.0D+00 * sin ( x * x - 3.0D+00 * x + 2.0D+00 )

  return
end
subroutine p09_f1 ( x, f1 )

!*****************************************************************************80
!
!! P09_F1 evaluates the first derivative for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F1, the first derivative of the 
!    objective function.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  f1 = 2.0D+00 * x &
    - 10.0D+00 * cos ( x * x - 3.0D+00 * x + 2.0D+00 ) &
    * ( 2.0D+00 * x - 3.0D+00 )

  return
end
subroutine p09_f2 ( x, f2 )

!*****************************************************************************80
!
!! P09_F2 evaluates the second derivative for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the values of the variables.
!
!    Output, real ( kind = 8 ) F2, the second derivative.
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) x

  f2 = 2.0D+00  &
    + 10.0D+00 * sin ( x * x - 3.0D+00 * x + 2.0D+00 ) &
    * ( 2.0D+00 * x - 3.0D+00 ) * ( 2.0D+00 * x - 3.0D+00 ) &
    - 20.0D+00 * cos ( x * x - 3.0D+00 * x + 2.0D+00 )

  return
end
subroutine p09_interval ( a, b )

!*****************************************************************************80
!
!! P09_INTERVAL returns a starting interval for optimization for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, two points defining an interval in which
!    the local minimizer should be sought.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a =  -5.0D+00
  b =  +5.0D+00

  return
end
subroutine p09_sol ( know, x )

!*****************************************************************************80
!
!! P09_SOL returns the solution for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KNOW.
!    If KNOW is 0, then the solution is not known.
!    If KNOW is positive, then the solution is known, and is returned in X.
!
!    Output, real ( kind = 8 ) X, the solution, if known.
!
  implicit none

  integer ( kind = 4 ) know
  real ( kind = 8 ) x

  know = 1

  x = 0.146621498932095D+00

  return
end
subroutine p09_start ( x )

!*****************************************************************************80
!
!! P09_START returns a starting point for optimization for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X, a starting point for the optimization.
!
  implicit none

  real ( kind = 8 ) x

  x = -2.0D+00

  return
end
subroutine p09_title ( title )

!*****************************************************************************80
!
!! P09_TITLE returns a title for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The oscillatory parabola'

  return
end
subroutine p10_f ( x, f )

!*****************************************************************************80
!
!! P10_F evaluates the objective function for problem 10.
!
!  Discussion:
!
!    This function is oscillatory.
!
!    The function has a local minimum at 1.7922 whose function value is
!    very close to the minimum value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Isabel Beichl, Dianne O'Leary, Francis Sullivan,
!    Monte Carlo Minimization and Counting: One, Two, Too Many,
!    Computing in Science and Engineering,
!    Volume 9, Number 1, January/February 2007.
!
!    Dianne O'Leary,
!    Scientific Computing with Case Studies,
!    SIAM, 2008,
!    ISBN13: 978-0-898716-66-5,
!    LC: QA401.O44.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) x

  f =           cos (           x ) &  
    + 5.0D+00 * cos ( 1.6D+00 * x ) &
    - 2.0D+00 * cos ( 2.0D+00 * x ) &
    + 5.0D+00 * cos ( 4.5D+00 * x ) &
    + 7.0D+00 * cos ( 9.0D+00 * x )

  return
end
subroutine p10_f1 ( x, f1 )

!*****************************************************************************80
!
!! P10_F1 evaluates the first derivative for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F1, the first derivative of the 
!    objective function.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  f1 = -                     sin (           x ) &  
       - 5.0D+00 * 1.6D+00 * sin ( 1.6D+00 * x ) &
       + 2.0D+00 * 2.0D+00 * sin ( 2.0D+00 * x ) &
       - 5.0D+00 * 4.5D+00 * sin ( 4.5D+00 * x ) &
       - 7.0D+00 * 9.0D+00 * sin ( 9.0D+00 * x )

  return
end
subroutine p10_f2 ( x, f2 )

!*****************************************************************************80
!
!! P10_F2 evaluates the second derivative for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the values of the variables.
!
!    Output, real ( kind = 8 ) F2, the second derivative.
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) x

  f2 = -                               cos (           x ) &  
       - 5.0D+00 * 1.6D+00 * 1.6D+00 * cos ( 1.6D+00 * x ) &
       + 2.0D+00 * 2.0D+00 * 2.0D+00 * cos ( 2.0D+00 * x ) &
       - 5.0D+00 * 4.5D+00 * 4.5D+00 * cos ( 4.5D+00 * x ) &
       - 7.0D+00 * 9.0D+00 * 9.0D+00 * cos ( 9.0D+00 * x )

  return
end
subroutine p10_interval ( a, b )

!*****************************************************************************80
!
!! P10_INTERVAL returns a starting interval for optimization for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, two points defining an interval in which
!    the local minimizer should be sought.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a =  0.0D+00
  b =  7.0D+00

  return
end
subroutine p10_sol ( know, x )

!*****************************************************************************80
!
!! P10_SOL returns the solution for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KNOW.
!    If KNOW is 0, then the solution is not known.
!    If KNOW is positive, then the solution is known, and is returned in X.
!
!    Output, real ( kind = 8 ) X, the solution, if known.
!
  implicit none

  integer ( kind = 4 ) know
  real ( kind = 8 ) x

  know = 1

  x = 5.975691087433868D+00

  return
end
subroutine p10_start ( x )

!*****************************************************************************80
!
!! P10_START returns a starting point for optimization for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X, a starting point for the optimization.
!
  implicit none

  real ( kind = 8 ) x

  x = 0.5D+00

  return
end
subroutine p10_title ( title )

!*****************************************************************************80
!
!! P10_TITLE returns a title for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The cosine combo'

  return
end
subroutine p11_f ( x, f )

!*****************************************************************************80
!
!! P11_F evaluates the objective function for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the objective function.
!
!    Output, real ( kind = 8 ) F, the value of the objective function.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) x

  f = 1.0D+00 + abs ( 3.0D+00 * x - 1.0D+00 )

  return
end
subroutine p11_f1 ( x, f1 )

!*****************************************************************************80
!
!! P11_F1 evaluates the first derivative for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) F1, the first derivative of the 
!    objective function.
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  if ( 3.0D+00 * x - 1.0D+00 < 0.0D+00 ) then
    f1 = - 3.0D+00
  else
    f1 = + 3.0D+00
  end if

  return
end
subroutine p11_f2 ( x, f2 )

!*****************************************************************************80
!
!! P11_F2 evaluates the second derivative for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the values of the variables.
!
!    Output, real ( kind = 8 ) F2, the second derivative.
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) x

  f2 = 0.0D+00

  return
end
subroutine p11_interval ( a, b )

!*****************************************************************************80
!
!! P11_INTERVAL returns a starting interval for optimization for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, two points defining an interval in which
!    the local minimizer should be sought.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a =  0.0D+00
  b =  1.0D+00

  return
end
subroutine p11_sol ( know, x )

!*****************************************************************************80
!
!! P11_SOL returns the solution for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KNOW.
!    If KNOW is 0, then the solution is not known.
!    If KNOW is positive, then the solution is known, and is returned in X.
!
!    Output, real ( kind = 8 ) X, the solution, if known.
!
  implicit none

  integer ( kind = 4 ) know
  real ( kind = 8 ) x

  know = 1

  x = 1.0D+00 / 3.0D+00

  return
end
subroutine p11_start ( x )

!*****************************************************************************80
!
!! P11_START returns a starting point for optimization for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X, a starting point for the optimization.
!
  implicit none

  real ( kind = 8 ) x

  x = 0.75D+00

  return
end
subroutine p11_title ( title )

!*****************************************************************************80
!
!! P11_TITLE returns a title for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = '1 + |3x-1|'

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
