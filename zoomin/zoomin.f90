subroutine bisect ( x, x1, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! BISECT carries out the bisection method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1, two points defining the interval
!    in which the search will take place.  F(X) and F(X1) should have
!    opposite signs.
!    On output, X is the best estimate for the root, and X1 is
!    a recently computed point such that F changes sign in the interval
!    between X and X1.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) temp
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
!
!  Initialization.
!
  ierror = 0
  k = 0
!
!  Evaluate the function at the starting points.
!
  fx = f ( x, 0 )
  fx1 = f ( x1, 0 )
!
!  Set XPOS and XNEG to the X values for which F(X) is positive
!  and negative, respectively.
!
  if ( 0.0D+00 <= fx .and. fx1 <= 0.0D+00 ) then

  else if ( fx <= 0.0D+00 .and. 0.0D+00 <= fx1 ) then
    call r8_swap ( x, x1 )
    call r8_swap ( fx, fx1 )
  else
    ierror = 1
    return
  end if
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      exit
    end if
!
!  Set the increment.
!
    dx = 0.5D+00 * ( x1 - x )
!
!  Update the iterate and function values.
!
    x2 = x + dx
    fx2 = f ( x2, 0 )

    if ( 0.0D+00 <= fx2 ) then
      x = x2
      fx = fx2
    else
      x1 = x2
      fx1 = fx2
    end if

  end do

  return
end
subroutine brent ( x, x1, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! BRENT implements the Brent bisection-based zero finder.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1.  On input, two points defining
!    the interval in which the search will take place.  F(X) and F(X1)
!    should have opposite signs.  On output, X is the best estimate for
!    the root, and X1 is a recently computed point such that F changes
!    sign in the interval between X and X1.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d
  real ( kind = 8 ) dx
  real ( kind = 8 ) e
  real ( kind = 8 ) em
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) tol
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  fx1 = f ( x1, 0 )
  x2 = x1
  fx2 = fx1
!
!  Iteration loop:
!
  do

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( 0.0D+00 < fx1 .eqv. 0.0D+00 < fx2 ) then
      x2 = x
      fx2 = fx
      d = x1 - x
      e = x1 - x
    end if

    if ( abs ( fx2 ) < abs ( fx1 ) ) then
      call r8_swap ( x1, x2 )
      call r8_swap ( fx1, fx2 )
    end if

    tol = 2.0D+00 * epsilon ( 1.0D+00 ) * abs ( x1 ) + abserr
    em = 0.5D+00 * ( x2 - x1 )

    if ( abs ( em ) <= tol .or. fx1 == 0.0D+00 ) then
      x = x1
      return
    end if

    if ( abs ( e ) < tol .or. abs ( fx ) <= abs ( fx1 ) ) then

      d = em
      e = em

    else

      s = fx1 / fx

      if ( x == x2 ) then
        p = 2.0D+00 * em * s
        q = 1.0D+00 - s
      else
        q = fx / fx2
        r = fx1 / fx2
        p = s * ( 2.0D+00 * em * q * ( q - r ) - ( x1 - x ) * ( r - 1.0D+00 ) )
        q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )
      end if

      if ( 0 < p ) then
        q = - q
      else
        p = - p
      end if

      s = e
      e = d

      if ( 2.0D+00 * p < 3.0D+00 * em * q - abs ( tol * q ) .or. &
        p < abs ( 0.5D+00 * s * q ) ) then
        d = p / q
      else
        d = em
        e = em
      end if

    end if
!
!  Set the increment.
!
    if ( tol < abs ( d ) ) then
      dx = + d
    else if ( 0.0D+00 < em ) then
      dx = + tol
    else
      dx = - tol
    end if
!
!  Remember current data for next step.
!
    x = x1
    fx = fx1
!
!  Update the iterate and function values.
!
    x1 = x1 + dx
    fx1 = f ( x1, 0 )

  end do

  return
end
subroutine c8_muller ( func, fatol, itmax, x1, x2, x3, xatol, xrtol, &
  xnew, fxnew )

!*****************************************************************************80
!
!! C8_MULLER carries out Muller's method, using C8 arithmetic.
!
!  Discussion:
!
!    "C8" arithmetic is complex double precision arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gisela Engeln-Muellges, Frank Uhlig,
!    Numerical Algorithms with C,
!    Springer, 1996,
!    ISBN: 3-540-60530-4,
!    LC: QA297.E56213.
!
!  Parameters:
!
!    Input, external FUNC, the name of the routine that evaluates the function.
!    FUNC should have the form:
!      subroutine func ( x, fx )
!      complex fx
!      complex x
!
!    Input, real ( kind = 8 ) FATOL, the absolute error tolerance for F(X).
!
!    Input, integer ( kind = 4 ) ( kind = 4 ) ITMAX, the maximum number of
!    steps allowed.
!
!    Input, complex ( kind = 8 ) X1, X2, X3, three distinct points to start the
!    iteration.
!
!    Input, real ( kind = 8 ) XATOL, XRTOL, absolute and relative
!    error tolerances for the root.
!
!    Output, complex ( kind = 8 ) XNEW, the estimated root.
!
!    Output, complex ( kind = 8 ) FXNEW, the value of the function at XNEW.
!
  implicit none

  complex ( kind = 8 ) a
  complex ( kind = 8 ) b
  complex ( kind = 8 ) c
  complex ( kind = 8 ) c8_temp
  complex ( kind = 8 ) discrm
  real ( kind = 8 ) fatol
  complex ( kind = 8 ) fminus
  complex ( kind = 8 ) fplus
  external func
  complex ( kind = 8 ) fxmid
  complex ( kind = 8 ) fxnew
  complex ( kind = 8 ) fxold
  integer ( kind = 4 ) iterate
  integer ( kind = 4 ) itmax
  real ( kind = 8 ) x_ave
  complex ( kind = 8 ) x_inc
  complex ( kind = 8 ) x1
  complex ( kind = 8 ) x2
  complex ( kind = 8 ) x3
  real ( kind = 8 ) xatol
  complex ( kind = 8 ) xlast
  complex ( kind = 8 ) xmid
  complex ( kind = 8 ) xminus
  complex ( kind = 8 ) xnew
  complex ( kind = 8 ) xold
  complex ( kind = 8 ) xplus
  real ( kind = 8 ) xrtol

  xnew = x1
  xmid = x2
  xold = x3

  call func ( xnew, fxnew )
  call func ( xmid, fxmid )
  call func ( xold, fxold )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C8_MULLER:'
  write ( *, '(a)' ) '  Muller''s method (complex root version)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Iteration     x_real              x_imag' // &
    '             ||fx||           ||disc||'
  write ( *, '(a)' ) ' '

  iterate = -2
  write ( *, '(i6,f20.10,f20.10,f20.10)' ) iterate, xold, abs ( fxold )
  iterate = -1
  write ( *, '(i6,f20.10,f20.10,f20.10)' ) iterate, xmid, abs ( fxmid )
  iterate = 0
  write ( *, '(i6,f20.10,f20.10,f20.10)' ) iterate, xnew, abs ( fxnew )

  if ( abs ( fxnew ) < fatol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C8_MULLER:'
    write ( *, '(a)' ) '  |F(X)| is below the tolerance.'
    return
  end if

  do
!
!  You may need to swap (XMID,FXMID) and (XNEW,FXNEW).
!
    if ( abs ( fxmid ) <= abs ( fxnew ) ) then

      c8_temp = xnew
      xnew = xmid
      xmid = c8_temp

      c8_temp = fxnew
      fxnew = fxmid
      fxmid = c8_temp

    end if

    xlast = xnew
    iterate = iterate + 1

    if ( itmax < iterate ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C8_MULLER:'
      write ( *, '(a)' ) '  Maximum number of steps taken.'
      exit
    end if

    a =  ( ( xmid - xnew ) * ( fxold - fxnew ) &
         - ( xold - xnew ) * ( fxmid - fxnew ) )

    b = ( ( xold - xnew )**2 * ( fxmid - fxnew ) &
        - ( xmid - xnew )**2 * ( fxold - fxnew ) )

    c = ( ( xold - xnew ) * ( xmid - xnew ) * ( xold - xmid ) * fxnew )

    xold = xmid
    xmid = xnew
!
!  Apply the quadratic formula to get roots XPLUS and XMINUS.
!
    discrm = b**2 - 4.0D+00 * a * c

    if ( a == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C8_MULLER:'
      write ( *, '(a)' ) '  The algorithm has broken down.'
      write ( *, '(a)' ) '  The quadratic coefficient A is zero.'
      exit
    end if

    xplus = xnew + ( ( - b + sqrt ( discrm ) ) / ( 2.0D+00 * a ) )

    call func ( xplus, fplus )

    xminus = xnew + ( ( - b - sqrt ( discrm ) ) / ( 2.0D+00 * a ) )

    call func ( xminus, fminus )
!
!  Choose the root with smallest function value.
!
    if ( abs ( fminus ) < abs ( fplus ) ) then
      xnew = xminus
    else
      xnew = xplus
    end if

    fxold = fxmid
    fxmid = fxnew
    call func ( xnew, fxnew )
    write ( *, '(i6,f20.10,f20.10,f20.10,f20.10)' ) &
      iterate, xnew, abs ( fxnew ), abs ( discrm )
!
!  Check for convergence.
!
    x_ave = abs ( xnew + xmid + xold ) / 3.0D+00
    x_inc = xnew - xmid

    if ( abs ( x_inc ) <= xatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C8_MULLER:'
      write ( *, '(a)' ) '  Absolute convergence of the X increment.'
      exit
    end if

    if ( abs ( x_inc ) <= xrtol * x_ave ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C8_MULLER:'
      write ( *, '(a)' ) '  Relative convergence of the X increment.'
      exit
    end if

    if ( abs ( fxnew ) <= fatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C8_MULLER:'
      write ( *, '(a)' ) '  Absolute convergence of |F(X)|.'
      exit
    end if

  end do

  return
end
subroutine cap_phi_03 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! CAP_PHI_03 implements the Traub capital Phi(0,3) method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 232.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    z = 1.0D+00 - 2.0D+00 * fx * d2fx / dfx**2
!
!  Set the increment.
!
    if ( 0.0D+00 < z ) then
      dx = - 2.0D+00 * ( fx / dfx ) / ( 1.0D+00 + sqrt ( z ) )
    else
      dx = - fx / dfx
    end if
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

  end do

  return
end
subroutine cap_phi_21 ( x, x1, x2, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! CAP_PHI_21 implements the Traub capital PHI(2,1) function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1, X2.
!    On input, three distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 and X2 are the
!    previous estimates.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) det
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = -2
  fx2 = f ( x2, 0 )
  k = -1
  fx1 = f ( x1, 0 )
  k = 0
  fx = f ( x, 0 )

  if ( x1 == x2 ) then
    ierror = 3
    return
  end if

  d1 = ( fx1 - fx2 ) / ( x1 - x2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    d2 = d1

    if ( x == x1 ) then
      ierror = 3
      return
    end if

    d1 = ( fx - fx1 ) / ( x - x1 )

    if ( x == x2 ) then
      ierror = 3
      return
    end if

    d2 = ( d1 - d2 ) / ( x - x2 )

    z = d1 + ( x - x1 ) * d2
    det = z**2 - 4.0D+00 * fx * d2
    det = max ( det, 0.0D+00 )
!
!  Set the increment.
!
    dx = - 2.0D+00 * fx / ( z + sqrt ( det ) )
!
!  Remember current data for next step.
!
    x2 = x1
    fx2 = fx1

    x1 = x
    fx1 = fx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine chebyshev ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! CHEBYSHEV implements Chebyshev's method.
!
!  Discussion:
!
!    This is also known as the Euler-Chebyshev method.
!
!    x(n+1) = x(n) - ( f(x(n)) / f'(x(n)) ) * ( 1 + 0.5 * L(x(n)) )
!
!    where
!
!      L(x) = f(x) * f''(x) / f'(x)**2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = d2fx * fx / dfx**2

    if ( 1.0D+00 + 0.5D+00 * u == 0.0D+00 ) then
      ierror = 4
      return
    end if
!
!  Set the increment.
!
    dx = - ( fx / dfx ) * ( 1.0D+00 + 0.5D+00 * u )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

  end do

  return
end
subroutine dagger_e12 ( x, x1, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! DAGGER_E12 implements the dagger E 1,2 algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 234.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1.
!    On input, two distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 is the
!    previous estimate.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfx1
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )

  dfx1 = f ( x1, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    if ( x - x1 == 0.0D+00 ) then
      ierror = 3
      return
    end if

    d = ( dfx - dfx1 ) / ( x - x1 )
!
!  Set the increment.
!
    dx = - ( fx / dfx ) - 0.5D+00 * ( fx / dfx )**2 * d / dfx
!
!  Remember current data for next step.
!
    x1 = x
    dfx1 = dfx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine e3 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! E3 implements the Traub E3 method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 232.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx
    v = d2fx / ( 2.0D+00 * dfx )
!
!  Set the increment.
!
    dx = - ( fx / dfx ) * ( 1.0D+00 + v * u )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine e4 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! E4 implements the Traub E4 method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 232.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) d3fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )
  d3fx = f ( x, 3 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx
    v = d2fx / ( 2.0D+00 * dfx )
    w = d3fx / ( 6.0D+00 * dfx )
!
!  Set the increment.
!
    dx = - u * ( 1.0D+00 + u * ( v + u * ( 2.0D+00 * v**2 - w ) ) )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )
    d3fx = f ( x, 3 )

  end do

  return
end
subroutine euler ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! EULER implements the Euler method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) bot
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    z = dfx**2 - 2.0D+00 * fx * d2fx
!
!  Set the increment.
!
    if ( 0.0D+00 < z ) then

      if ( dfx + sqrt ( z ) == 0.0D+00 ) then
        ierror = 3
        return
      end if

      dx = - 2.0D+00 * fx / ( dfx + sqrt ( z ) )

    else

      if ( dfx == 0.0D+00 ) then
        ierror = 3
        return
      end if

      dx = - fx / dfx

    end if
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

  end do

  return
end
subroutine halley1 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! HALLEY1 implements Halley's method.
!
!  Discussion:
!
!    x(n+1) = x(n) - ( f(x(n)) / f'(x(n)) ) / ( 1 - 0.5 * L(x(n)) )
!
!    where
!
!      L(x) = f(x) * f''(x) / f'(x)**2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = d2fx * fx / dfx**2

    if ( 2.0D+00 - u == 0.0D+00 ) then
      ierror = 4
      return
    end if
!
!  Set the increment.
!
    dx = - ( fx / dfx ) / ( 1.0D+00 - 0.5D+00 * u )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

  end do

  return
end
subroutine halley2 ( x, x1, x2, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! HALLEY2 implements Halley's method, with finite differences.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1, X2.
!    On input, three distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 and X2 are the
!    previous estimates.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  fx1 = f ( x1, 0 )
  fx2 = f ( x2, 0 )

  d1 = ( fx1 - fx2 ) / ( x1 - x2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    d2 = d1

    if ( x == x1 ) then
      ierror = 3
      return
    end if

    d1 = ( fx - fx1 ) / ( x - x1 )

    if ( x == x2 ) then
      ierror = 3
      return
    end if

    d2 = ( d1 - d2 ) / ( x - x2 )

    if ( d1 == 0.0D+00 ) then
      ierror = 3
      return
    end if

    z = d1 - fx1 * d2 / d1

    if ( z == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - fx / z
!
!  Remember current data for next step.
!
    x2 = x1
    fx2 = fx1
    x1 = x
    fx1 = fx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine halley_super ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! HALLEY_SUPER implements the super Halley method.
!
!  Discussion:
!
!    x(n+1) = x(n) - 0.5 * ( f(x(n)) / f'(x(n)) )
!      * ( 2 - L(x(n)) ) / ( 1 - L(x(n)) )
!
!    where
!
!      L(x) = f(x) * f''(x) / f'(x)**2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = d2fx * fx / dfx**2

    if ( 1.0D+00 - u == 0.0D+00 ) then
      ierror = 4
      return
    end if
!
!  Set the increment.
!
    dx = - 0.5D+00 * ( fx / dfx ) * ( 2.0D+00 - u ) / ( 1.0D+00 - u )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

  end do

  return
end
subroutine hansen ( x, beta, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! HANSEN implements the Hansen and Patrick method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eldon Hansen, Merrell Patrick,
!    A Family of Root Finding Methods,
!    Numerische Mathematik,
!    Volume 27, 1977, pages 257 - 269.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) BETA, ???
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) beta
  real ( kind = 8 ) bot
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )

  write ( *, '(a,i4,2g14.6)' ) 'KXFXDFX', k, x, fx
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    z = dfx**2 - ( beta + 1.0D+00 ) * fx * d2fx
    z = max ( z, 0.0D+00 )

    bot = beta * dfx + sqrt ( z )

    if ( bot == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - ( beta + 1.0D+00 ) * fx / bot
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

    write ( *, '(a,i4,3g14.6)' ) 'HANSEN - KXFXDFX', k, x, fx, dfx

  end do

  return
end
subroutine jarratt ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! JARRATT implements the Jarratt fourth order method.
!
!  Discussion:
!
!    Jarratt's method is of fourth order.
!
!    x(n+1) = x(n) - 1/2 * ( f(x(n) / f'(x(n)) )
!      + f(x(n)) / ( f'(x(n)) - 3*f'( x(n)-(2/3)*f(x(n))/f'(x(n)) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    P Jarratt,
!    Some fourth-order multipoint iterative methods for solving equations,
!    Mathematics of Computation,
!    Volume 20, Number 95, 1966, pages 434 - 437.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfz
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    z = x - 2.0D+00 * fx / ( 3.0D+00 * dfx )
    dfz = f ( z, 1 )

    if ( dfx - 3.0D+00 * dfz == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = ( - 0.5D+00 / dfx + 1.0D+00 / ( dfx - 3.0D+00 * dfz ) ) * fx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine jarratt2 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! JARRATT2 implements the inverse-free Jarratt fourth order method.
!
!  Discussion:
!
!    Jarratt's inverse-free method is of fourth order.
!
!    x(n+1) = x(n) - u(x(n)) + (3/4) * u(x(n)) * h(x(n)) * ( 1-(3/2)*h(x(n)) )
!
!    where
!
!      u(x) = f(x) / f'(x)
!      h(x) = ( f'(x-(2/3)*u(x)) - f'(x) ) / f'(x)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    P Jarratt,
!    Some fourth-order multipoint iterative methods for solving equations,
!    Mathematics of Computation,
!    Volume 20, Number 95, 1966, pages 434 - 437.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) hx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) ux
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    ux = fx / dfx

    hx = ( f ( x - 2.0D+00 / 3.0D+00 * ux, 1 ) - dfx ) / dfx
!
!  Set the increment.
!
    dx = - ux + 0.75D+00 * ux * hx * ( 1.0D+00 - 1.5D+00 * hx )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine king ( x, beta, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! KING implements a family of fourth order methods.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard King,
!    A family of fourth order methods,
!    SIAM Journal on Numerical Analysis,
!    Volume 10, 1973, pages 876 - 879.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) BETA, a parameter in the algorithm.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) beta
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fw
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) w
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    w = x - fx / dfx
    fw = f ( w, 0 )

    if ( fx + ( beta - 2.0D+00 ) * fw == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - ( fx / dfx ) - ( fw / dfx ) &
      * ( fx + beta * fw ) / ( fx + ( beta - 2.0D+00 ) * fw )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine laguerre ( x, ipoly, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! LAGUERRE implements the Laguerre rootfinding method for polynomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eldon Hansen, Merrell Patrick,
!    A Family of Root Finding Methods,
!    Numerische Mathematik,
!    Volume 27, 1977, pages 257 - 269.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, integer ( kind = 4 ) IPOLY, the polynomial degree of the function.
!    IPOLY must be at least 2.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) beta
  real ( kind = 8 ) bot
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ipoly
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Check.
!
  if ( ipoly < 2 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE - Fatal error!'
    write ( *, '(a)' ) '  IPOLY < 2.'
    return
  end if
!
!  Initialization.
!
  ierror = 0

  beta = 1.0D+00 / real ( ipoly - 1, kind = 8 )

  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    z = dfx**2 - ( beta + 1.0D+00 ) * fx * d2fx
    z = max ( z, 0.0D+00 )

    bot = beta * dfx + sqrt ( z )

    if ( bot == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - ( beta + 1.0D+00 ) * fx / bot
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

  end do

  return
end
subroutine midpoint ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! MIDPOINT implements the midpoint method.
!
!  Discussion:
!
!    The midpoint method is of order 3.
!
!    x(n+1) = x(n) - f(x(n)) / g(x(n))
!
!    where
!
!    g(x) = f'( x - f(x) / ( 2 * f'(x) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 164.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, the point that starts the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) arg
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) gx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0

  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    arg = x - fx / ( 2.0D+00 * dfx )
    gx = f ( arg, 1 )

    if ( gx == 0.0D+00 ) then
      ierror = 4
      return
    end if
!
!  Set the increment.
!
    dx = - fx / gx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine newton ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! NEWTON implements Newton's method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - fx / dfx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine newton_mod ( x, nsub, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! NEWTON_MOD implements the modified Newton method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 236.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, integer ( kind = 4 ) NSUB, the number of steps in the sub-iteration.
!    The derivative is only evaluated before the first step,
!    and after every subsequent NSUB steps.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) nsub
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( mod ( k - 1, nsub ) == 0 ) then

      dfx = f ( x, 1 )

      if ( dfx == 0.0D+00 ) then
        ierror = 3
        return
      end if

    end if
!
!  Set the increment.
!
    dx = - fx / dfx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine newton_sec ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! NEWTON_SEC implements the Newton - secant method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 236.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fxu
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    x1 = x - fx / dfx
    fxu = f ( x1, 0 )

    if ( fxu - fx == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - ( fx / dfx ) * ( 1.0D+00 - fxu / ( fxu - fx ) )
!
!  Remember current data for next step.
!
    x1 = x
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine ostrowski_sqrt ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! OSTROWSKI_SQRT implements the Ostrowski square root method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eldon Hansen, Merrell Patrick,
!    A Family of Root Finding Methods,
!    Numerische Mathematik,
!    Volume 27, 1977, pages 257 - 269.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) BETA, ???
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) bot
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    z = dfx**2 - fx * d2fx
!
!  Set the increment.
!
    if ( 0.0D+00 < z ) then
      bot = sqrt ( z )
    else
      bot = dfx
    end if

    dx = - fx / bot
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

  end do

  return
end
subroutine perp_e_12 ( x, x1, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! PERP_E_12 implements the Traub E 1,2 algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 234.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1.
!    On input, two distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 is the previous
!    estimate.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfx1
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  fx1 = f ( x1, 0 )
  dfx1 = f ( x1, 1 )

  if ( dfx1 == 0.0D+00 ) then
    ierror = 3
    return
  end if
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    if ( dfx1 == 0.0D+00 ) then
      ierror = 3
      return
    end if

    if ( fx == fx1 ) then
      ierror = 3
      return
    end if

    z = 2.0D+00 / dfx + 1.0D+00 / dfx1 - 3.0D+00 * ( x - x1 ) / ( fx - fx1 )
!
!  Set the increment.
!
    dx = - fx / dfx + fx**2 * z / ( fx - fx1 )
!
!  Remember current data for next step.
!
    x1 = x
    fx1 = fx
    dfx1 = dfx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine perp_e_21 ( x, x1, x2, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! PERP_E_21 implements the Traub perp E 21 method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 233.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1, X2.
!    On input, three distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 and X2 are the
!    previous estimates.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  fx1 = f ( x1, 0 )
  fx2 = f ( x2, 0 )

  if ( x1 == x2 ) then
    ierror = 3
    return
  end if

  d1 = ( fx1 - fx2 ) / ( x1 - x2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    d2 = d1

    if ( x == x1 ) then
      ierror = 3
      return
    end if

    if ( x == x2 ) then
      ierror = 3
      return
    end if

    d1 = ( fx - fx1 ) / ( x - x1 )
    d = ( fx - fx2 ) / ( x - x2 )
!
!  Set the increment.
!
    dx = - fx * ( 1.0D+00 / d1 + 1.0D+00 / d - 1.0D+00 / d2 )
!
!  Remember current data for next step.
!
    x2 = x1
    fx2 = fx1
    x1 = x
    fx1 = fx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine phi_12 ( x, x1, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! PHI_12 implements the Traub capital PHI(1,2) method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 233.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1.
!    On input, two distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 is the previous
!    estimate.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfx1
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) h
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  fx1 = f ( x1, 0 )
  dfx1 = f ( x1, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    c = fx - fx1

    if ( x == x1 ) then
      ierror = 3
      return
    end if

    d = ( fx - fx1 ) / ( x - x1 )

    h = 1.0D+00 / c * ( 1.0D+00 / dfx - 1.0D+00 / d ) - ( fx1 / c**2 ) * &
      ( 1.0D+00 / dfx + 1.0D+00 / dfx1 - 2.0D+00 / d )
!
!  Set the increment.
!
    dx = - ( fx / dfx ) + fx**2 * h
!
!  Remember current data for next step.
!
    x1 = x
    fx1 = fx
    dfx1 = dfx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine psi_21 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! PSI_21 implements the Traub PSI 2,1 method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 232.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) d3fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )
    d3fx = f ( x, 3 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx
    v = d2fx / ( 2.0D+00 * dfx )
    w = d3fx / ( 6.0D+00 * dfx )
!
!  Set the increment.
!
    dx = - u * ( v - u * ( v**2 - w ) ) / ( v - u * ( 2.0D+00 * v**2 - w ) )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine psi_i2 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! PSI_12 implements the Traub PSI 1,2 method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 232.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) d3fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )
    d3fx = f ( x, 3 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx
    v = d2fx / ( 2.0D+00 * dfx )
    w = d3fx / ( 6.0D+00 * dfx )
!
!  Set the increment.
!
    dx = - u / ( 1.0D+00 - u * v + ( v**2 - w ) * u**2 )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine r8_muller ( x, x1, x2, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! R8_MULLER implements Muller's method for real arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1, X2.
!    On input, three distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 and X2 are the
!    previous estimates.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) q
  real ( kind = 8 ) term
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx2 = f ( x2, 0 )
  fx1 = f ( x1, 0 )
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    q = ( x - x1 ) / ( x1 - x2 )

    a = q * fx - q * ( 1.0D+00 + q ) * fx1 + q**2 * fx2
    b = ( 2.0D+00 * q + 1.0D+00 ) * fx - ( 1.0D+00 + q )**2 * fx1 + q**2 * fx2
    c = ( 1.0D+00 + q ) * fx

    term = b**2 - 4.0D+00 * a * c
    term = max ( term, 0.0D+00 )
    term = sqrt ( term )
    if ( b < 0.0D+00 ) then
      term = - term
    end if
!
!  Set the increment.
!
    dx = - ( x - x1 ) * 2.0D+00 * c / ( b + term )
!
!  Remember current data for next step.
!
    x2 = x1
    fx2 = fx1
    x1 = x
    fx1 = fx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP switches two R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine red_cap_phi_04 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! RED_CAP_PHI_04 implements the Traub reduced capital PHI(0,4) method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 232.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) d3fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )
    d3fx = f ( x, 3 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx
    v = d2fx / ( 2.0D+00 * dfx )
    w = d3fx / ( 6.0D+00 * dfx )
    z = 1.0D+00 - 4.0D+00 * u * ( v - u * w )
!
!  Set the increment.
!
    if ( 0.0D+00 < z ) then
      dx = - 2.0D+00 * u / ( 1.0D+00 + sqrt ( z ) )
    else
      dx = - fx / dfx
    end if
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine regula ( x, x1, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! REGULA implements the Regula Falsi method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1.
!    On input, two distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 is the previous
!    estimate.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) temp
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
!
!  Initialization.
!
  ierror = 0
  k = 0

  fx = f ( x, 0 )
  fx1 = f ( x1, 0 )

  if ( fx < 0.0D+00  ) then
    call r8_swap ( x, x1 )
    call r8_swap ( fx, fx1 )
  end if
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if
!
!  Set the increment.
!
    dx = - fx1 * ( x - x1 ) / ( fx - fx1 )
!
!  Update the iterate and function values.
!
    x2 = x1 + dx
    fx2 = f ( x2, 0 )

    if ( 0.0D+00 <= fx2 ) then
      x = x2
      fx = fx2
    else
      x1 = x2
      fx1 = fx2
    end if

  end do

  return
end
subroutine rhein1 ( x, x1, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! RHEIN1 implements the Rheinboldt bisection - secant method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Werner Rheinboldt,
!    Algorithms for finding zeros of a function,
!    UMAP Journal,
!    Volume 2, Number 1, 1981, pages 43 - 72.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1.
!    On input, two distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 is the previous
!    estimate.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dx
  real ( kind = 8 ) em
  real ( kind = 8 ), external :: f
  logical forced
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
!
!  Initialization.
!
  i = 0
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  fx1 = f ( x1, 0 )

  if ( abs ( fx1 ) < abs ( fx ) ) then
    call r8_swap ( x, x1 )
    call r8_swap ( fx, fx1 )
  end if

  x2 = x1
  fx2 = fx1
  t = 0.5D+00 * abs ( x - x1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if
!
!  Force ABS ( FX ) <= ABS ( FX1 ).
!
    if ( abs ( fx1 ) < abs ( fx ) ) then

      x2 = x
      x = x1
      x1 = x2

      fx2 = fx
      fx = fx1
      fx1 = fx2

    end if

    em = 0.5D+00 * ( x1 - x )
!
!  Compute the numerator and denominator for secant step.
!
    p = ( x - x2 ) * fx
    q = fx2 - fx

    if ( p < 0.0D+00 ) then
      p = - p
      q = - q
    end if
!
!  Save the old minimum residual point.
!
    x2 = x
    fx2 = fx
!
!  Test for forced bisection.
!
    i = i + 1
    forced = .false.

    if ( 3 < i ) then
      if ( t < 8.0D+00 * abs ( em ) ) then
        forced = .true.
      else
        i = 0
        t = em
      end if
    end if
!
!  Set the increment.
!
    if ( forced ) then
      dx = em
    else if ( p <= abs ( q ) * abserr ) then
      dx = sign ( 1.0D+00, em ) * abserr
    else if ( p < q * em ) then
      dx = p / q
    else
      dx = em
    end if
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
!
!  Preserve the change of sign interval.
!
    if ( sign ( 1.0D+00, fx ) == sign ( 1.0D+00, fx1 ) ) then
      x1 = x2
      fx1 = fx2
    end if

  end do

  return
end
subroutine rhein2 ( x, x1, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! RHEIN2 implements the Rheinboldt bisection/secant/inverse quadratic method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Werner Rheinboldt,
!    Algorithms for Finding Zeros of a Function,
!    UMAP Journal,
!    Volume 2, Number 1, 1981, pages 43 - 72.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1.
!    On input, two distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 is the previous
!    estimate.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dx
  real ( kind = 8 ) em
  real ( kind = 8 ), external :: f
  logical forced
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) piq
  real ( kind = 8 ) ps
  real ( kind = 8 ) qiq
  real ( kind = 8 ) qs
  real ( kind = 8 ) stpmin
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
!
!  Initialization.
!
  i = 0
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  fx1 = f ( x1, 0 )

  if ( abs ( fx1 ) < abs ( fx ) ) then
    call r8_swap ( x, x1 )
    call r8_swap ( fx, fx1 )
  end if

  x2 = x1
  fx2 = fx1
  k = 0
  t = 0.5D+00 * abs ( x - x1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    i = i + 1

    em = 0.5D+00 * ( x1 - x )
!
!  Compute the numerator and denominator for secant step.
!
    if ( 2.0D+00 * abs ( x2 - x ) < abs ( x1 - x ) ) then
      ps = ( x - x2 ) * fx
      qs = fx2 - fx
    else
      ps = ( x - x1 ) * fx
      qs = fx1 - fx
    end if

    if ( ps < 0.0D+00 ) then
      ps = - ps
      qs = - qs
    end if
!
!  Compute the numerator and denominator for inverse quadratic.
!
    piq = 0.0D+00
    qiq = 0.0D+00

    if ( x1 /= x2 ) then
      u = fx / fx2
      v = fx2 / fx1
      w = fx / fx1
      piq = u * ( 2.0D+00 * em * v * ( v - w ) - ( x - x2 ) * ( w - 1.0D+00 ) )
      qiq = ( u - 1.0D+00 ) * ( v - 1.0D+00 ) * ( w - 1.0D+00 )

      if ( 0.0D+00 < piq ) then
        qiq =  - qiq
      end if

      piq = abs ( piq )
    end if
!
!  Save the old minimum residual point.
!
    x2 = x
    fx2 = fx

    stpmin = ( abs ( x ) + abs ( em ) + 1.0D+00 ) * abserr
!
!  Choose bisection, secant or inverse quadratic step.
!
    forced = .false.

    if ( 3 < i ) then
      if ( t < 8.0D+00 * abs ( em ) ) then
        forced = .true.
      else
        i = 0
        t = em
      end if
    end if
!
!  Set the increment.
!
    if ( forced ) then
      dx = em
    else if ( piq < 1.5D+00 * em * qiq .and. &
      abs ( qiq ) * stpmin < abs ( piq ) ) then
      dx = piq / qiq
    else if ( ps < qs * em .and. abs ( qs ) * stpmin < abs ( ps ) ) then
      dx = ps / qs
    else
      dx = em
    end if
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
!
!  Set the new X1 as either X1 or X2, depending on whether
!  F(X1) or F(X2) has the opposite sign from F(X).
!
    if ( sign ( 1.0D+00, fx ) == sign ( 1.0D+00, fx1 ) ) then
      x1 = x2
      fx1 = fx2
    end if
!
!  Force ABS ( FX ) <= ABS ( FX1 ).
!
    if ( abs ( fx1 ) < abs ( fx ) ) then
      call r8_swap ( x, x1 )
      call r8_swap ( fx, fx1 )
    end if

  end do

  return
end
subroutine script_e2 ( x, m, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! SCRIPT_E2 implements the Traub script E - 2 function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 234.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, integer ( kind = 4 ) M, the multiplicity of the root.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) m
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - real ( m, kind = 8 ) * fx / dfx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine script_e3 ( x, m, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! SCRIPT_E3 implements the Traub script E - 3 function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 235.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, integer ( kind = 4 ) M, the multiplicity of the root.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) a2
  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ) em
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) m
  real ( kind = 8 ) u
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  em = real ( m, kind = 8 )
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx
    a2 = d2fx / ( 2.0D+00 * dfx )
!
!  Set the increment.
!
    dx = - em * u * ( ( 3.0D+00 - em ) / 2.0D+00 + em * a2 * u )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine script_e4 ( x, m, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! SCRIPT_E4 implements the Traub script E - 4 function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 235.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, integer ( kind = 4 ) M, the multiplicity of the root.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) a2
  real ( kind = 8 ) a3
  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) d3fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ) em
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) m
  real ( kind = 8 ) u
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  em = real ( m, kind = 8 )
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )
    d3fx = f ( x, 3 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx
    a2 = d2fx / ( 2.0D+00 * dfx )
    a3 = d3fx / ( 6.0D+00 * dfx )
!
!  Set the increment.
!
    dx = - em * u * ( ( em**2 - 6.0D+00 * em + 11.0D+00 ) / 6.0D+00 &
       + em * ( 2.0D+00 - em ) * a2 * u &
       + em**2 * ( 2.0D+00 * a2**2 - a3 ) * u**2 )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine secant ( x, x1, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! SECANT implements the secant method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1.
!    On input, two distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 is the previous
!    estimate.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  fx1 = f ( x1, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( fx - fx1 == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - fx * ( x - x1 ) / ( fx - fx1 )
!
!  Remember current data for next step.
!
    x1 = x
    fx1 = fx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine secant_extended ( x, x1, x2, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! SECANT_EXTENDED carries out the extended secant algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1, X3.
!    On input, three distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 and X2 are two previous
!    estimates.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
!
!  Initialization.
!
  ierror = 0
  k = 0

  fx1 = f ( x1, 0 )
  fx2 = f ( x2, 0 )
  d1 = ( fx1 - fx2 ) / ( x1 - x2 )

  if ( d1 == 0.0D+00 ) then
    ierror = 3
    return
  end if

  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    d2 = d1

    if ( x == x1 ) then
      ierror = 3
      return
    end if

    d1 = ( fx - fx1 ) / ( x - x1 )

    if ( d1 == 0.0D+00 ) then
      ierror = 3
      return
    end if

    if ( fx2 == fx1 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - fx / d1 + ( fx * fx1 / ( fx - fx2 ) ) &
      * ( 1.0D+00 / d1 - 1.0D+00 / d2 )
!
!  Remember current data for next step.
!
    x2 = x1
    fx2 = fx1
    x1 = x
    fx1 = fx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine star_e12 ( x, x1, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! STAR_E12 implements the Traub *E12 method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 234.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1.
!    On input, two distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 is the previous
!    estimate.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfx1
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )

  fx1 = f ( x1, 0 )
  dfx1 = f ( x1, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx

    if ( x == x1 ) then
      ierror = 3
      return
    end if

    d = ( fx - fx1 ) / ( x - x1 )
    z = 2.0D+00 * dfx + dfx1 - 3.0D+00 * d
!
!  Set the increment.
!
    dx = - u - u**2 * z / ( dfx * ( x - x1 ) )
!
!  Remember current data for next step.
!
    x1 = x
    fx1 = fx
    dfx1 = dfx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine star_e21 ( x, x1, x2, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! STAR_E21 implements the Traub *E21 method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 234.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1, X2.
!    On input, three distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 and X2 are the
!    previous estimates.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  fx1 = f ( x1, 0 )
  fx2 = f ( x2, 0 )
  d1 = ( fx1 - fx2 ) / ( x1 - x2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    d2 = d1

    if ( x == x1 ) then
      ierror = 3
      return
    end if

    d1 = ( fx - fx1 ) / ( x - x1 )

    if ( x == x2 ) then
      ierror = 3
      return
    end if

    d = ( fx - fx2 ) / ( x - x2 )

    if ( d1 + d - d2 == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - fx / ( d1 + d - d2 )
!
!  Remember current data for next step.
!
    x2 = x1
    fx2 = fx1

    x1 = x
    fx1 = fx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine steffenson ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! STEFFENSON implements Steffenson's method.
!
!  Discussion:
!
!    Steffenson's method is of order 2.
!
!    x(n+1) = x(n) - f(x(n)) / g(x(n))
!
!    where
!
!    g(x) = ( f(x+f(x)) - f(x) ) / f(x)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 178.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, the point that starts the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) gx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0

  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    gx = ( f ( x + fx, 0 ) - fx ) / fx

    if ( gx == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - fx / gx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine stirling ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! STIRLING implements Stirling's method.
!
!  Discussion:
!
!    Stirling's method is of order 2.
!
!    x(n+1) = x(n) - f(x(n)) / f'( x(n) - f(x(n) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, the point that starts the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) gx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0

  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    gx = f ( x - fx, 1 )

    if ( gx == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - fx / gx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine t14 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! T14 implements the Traub fourteenth function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 237.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfxu
  real ( kind = 8 ) dfz
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end  if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    dfxu = f ( x - ( fx / dfx ), 1 )

    if ( dfxu == 0.0D+00 ) then
      ierror = 3
      return
    end if

    z = x - 0.25D+00 * ( fx / dfx + fx / dfxu  )
    dfz = f ( z, 1 )
!
!  Set the increment.
!
    dx = - ( fx / dfx + ( fx / dfxu ) + 4.0D+00 * ( fx / dfz ) ) / 6.0D+00
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine t15 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! T15 implements the Traub fifteenth function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 238.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfxu
  real ( kind = 8 ) dfz
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx
    dfxu = f ( x - ( fx / dfx ), 1 )

    if ( dfxu == 0.0D+00 ) then
      ierror = 3
      return
    end if

    z = x - 2.0D+00 * ( 2.0D+00 * u + fx / dfxu ) / 9.0D+00
    dfz = f ( z, 1 )

    if ( dfz == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - 0.25D+00 * ( fx / dfx + 3.0D+00 * fx / dfz )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine t16 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! T16 implements the Traub sixteenth function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 238.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfxu
  real ( kind = 8 ) dfz
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    dfxu = f ( x - ( fx / ( 3.0D+00 * dfx ) ), 1 )

    if ( dfxu == 0.0D+00 ) then
      ierror = 3
      return
    end if

    z = x - 2.0D+00 * fx / ( 3.0D+00 * dfxu )
    dfz = f ( z, 1 )

    if ( dfz == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - 0.25D+00 * ( fx / dfx + 3.0D+00 * fx / dfz )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine t_family1 ( x, c, d, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! T_FAMILY1 implements the Traub first family of iterations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, pages 236 - 237.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) C, D, ???
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfz
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    z = x - d * fx / dfx
    dfz = f ( z, 1 )

    if ( dfz == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - ( c / dfx + ( 1.0D+00 - c ) / dfz ) * fx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine t_family2 ( x, a, b, c, d, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! T_FAMILY2 implements the Traub second family of iterations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 236 - 237.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) A, B, C, D, ???
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfz
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    z = x - d * fx / dfx
    dfz = f ( z, 1 )
!
!  Set the increment.
!
    dx = - fx * ( b * dfx - c * dfz ) / ( a * dfx**2 )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine te11f ( x, x1, m, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! TE11F implements the Traub *E - 1,1(f) function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 235.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, X1.
!    On input, two distinct points that start the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR, and X1 is the previous
!    estimate.
!
!    Input, integer ( kind = 4 ) M, the multiplicity of the root.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d
  real ( kind = 8 ) dx
  real ( kind = 8 ) expon
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fnorm
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) gx
  real ( kind = 8 ) gx1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) m
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
!
!  Initialization.
!
  ierror = 0
  k = 0
  expon = 1.0D+00 / real ( m, kind = 8 )
  fx1 = f ( x1, 0 )

  fnorm = abs ( fx1 )
  gx1 = sign ( 1.0D+00, fx1 ) * fnorm**expon
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    fnorm = abs ( fx )
    gx = sign ( 1.0D+00, fx ) * fnorm**expon

    if ( x == x1 ) then
      ierror = 3
      return
    end if

    d = ( gx - gx1 ) / ( x - x1 )
!
!  Set the increment.
!
    dx = - gx / d
!
!  Remember current data for next step.
!
    gx1 = gx
    x1 = x
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine te2u ( x, em, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! TE2U implements the Traub E - 2(u) function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 235.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Output, real ( kind = 8 ) EM, an estimate of the multiplicity of the root.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ) em
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    if ( dfx**2 - fx * d2fx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    em = dfx**2 / ( dfx**2 - fx * d2fx )
    em = max ( em, 1.0D+00 )
!
!  Set the increment.
!
    dx = - em * fx / dfx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

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
subroutine tphi1u ( x, em, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! TPHI1U implements the Traub phi - 1,1(u) function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 235.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Output, real ( kind = 8 ) EM, an estimate of the multiplicity of the root.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ) em
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fu
  real ( kind = 8 ) fu1
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) u1
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  u1 = f ( x, 0 ) / f ( x, 1 )
  x = x - u1
  fx = f ( x, 0 )
  dfx = f ( x, 1 )

  fu1 = f ( u1, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx
    fu = f ( u, 0 )

    em = ( u - u1 ) / ( fu - fu1 )
    em = max ( em, 1.0D+00 )
!
!  Set the increment
!
    dx = - em * fx / dfx
!
!  Remember current data for next step.
!
    u1 = u
    fu1 = fu
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine traub1 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! TRAUB1 implements the Traub first method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 236.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfz
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    z = x - fx / dfx
    dfz = f ( z, 1 )

    if ( dfz == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - fx / dfz
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine traub4 ( x, nsub, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! TRAUB4 implements the Traub fourth method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 236.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, integer ( kind = 4 ) NSUB, the number of steps in the sub-iteration.
!    The derivatives are only evaluated before the first step,
!    and after every subsequent NSUB steps.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) nsub
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( mod ( k - 1, nsub ) == 0 ) then

      dfx = f ( x, 1 )
      d2fx = f ( x, 2 )

      if ( dfx == 0.0D+00 ) then
        ierror = 3
        return
      end if

      if ( dfx - d2fx * fx / dfx == 0.0D+00 ) then
        ierror = 3
        return
      end if

    end if
!
!  Set the increment.
!
    dx = - fx / ( dfx - d2fx * fx / dfx )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine traub8 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! TRAUB8 implements the Traub eighth function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfxu
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    dfxu = f ( x - 2.0D+00 * fx / ( 3.0D+00 * dfx ), 1 )
!
!  Set the increment.
!
    dx = - 4.0D+00 * fx / ( dfx + 3.0D+00 * dfxu )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine traub9 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! TRAUB9 implements the Traub ninth method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 237.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fxu
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) xu
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx
    xu = x - fx / dfx
    fxu = f ( xu, 0 )

    if ( 2.0D+00 * fxu - fx == 0.0D+00 ) then
      ierror = 4
      return
    end if
!
!  Set the increment.
!
    dx = - fx / dfx + u * fxu / ( 2.0D+00 * fxu - fx )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine traub_ostrowski ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! TRAUB_OSTROWSKI implements the Traub-Ostrowksi method.
!
!  Discussion:
!
!    The Traub-Ostrowski method is of order 4.
!
!    x(n+1) = x(n) - u(x(n)) * ( f ( x(n) - u(x(n)) - f(x(n)) ) /
!      ( 2 * f ( x(n) - u(x(n)) - f(x(n)) )
!
!    where
!
!    u(x) = f(x) / f'(x)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 184.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, the point that starts the method.
!    On output, X is an approximation to a root of the equation
!    which satisfies abs ( F(X) ) < ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) gx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) ux
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0

  fx = f ( x, 0 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    dfx = f ( x, 1 )

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    ux = f ( x, 0 ) / dfx

    fx1 = f ( x - ux, 0 )

    if ( 2.0D+00 * fx1 - fx == 0.0D+00 ) then
      ierror = 4
      return
    end if
!
!  Set the increment.
!
    dx = - ux * ( fx1 - fx ) / ( 2.0D+00 * fx1 - fx )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )

  end do

  return
end
subroutine tt1f ( x, a, rho, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! TT1F implements the Traub type 1 functions 10 and 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 237.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) A, ???
!
!    Input, real ( kind = 8 ) RHO, ???
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fxu
  real ( kind = 8 ) fz
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) rho
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0

  if ( rho == 0.0D+00 ) then
    ierror = 3
    return
  end if

  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    z = x + rho * fx / dfx
    fxu = f ( z, 0 )

    z = x - fxu / ( rho**2 * dfx )
    fz = f ( z, 0 )
!
!  Set the increment.
!
    dx = - ( a * fz + fxu / rho**2 ) / dfx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine tthip ( x, em, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! TTHIP implements the Traub third function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph Traub,
!    Iterative Methods for the Solution of Equations,
!    Prentice Hall, 1964, page 235.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Output, real ( kind = 8 ) EM, an estimate of the multiplicity of the root.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ) em
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    if ( fx == 0.0D+00 ) then
      em = 1.0D+00
    else if ( log ( abs ( fx / dfx ) ) == 0.0D+00 ) then
      em = 1.0D+00
    else
      em = log ( abs ( fx ) ) / log ( abs ( fx / dfx ) )
    end if

    em = max ( em, 1.0D+00 )
!
!  Set the increment.
!
    dx = - em * fx / dfx
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine vandev ( x, em, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! VANDEV implements the Van de Vel iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hugo vandeVel,
!    A method for computing a root of a single nonlinear equation,
!    including its multiplicity,
!    Computing,
!    Volume 14, 1975, pages 167 - 171.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Output, real ( kind = 8 ) EM, an estimate of the multiplicity of the root.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dfz
  real ( kind = 8 ) dx
  real ( kind = 8 ) em
  real ( kind = 8 ) em1
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fz
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) u1
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Initialization.
!
  ierror = 0
  k = 0
  em1 = 1.0D+00
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx
    z = x - em1 * fx / dfx
    fz = f ( z, 0 )
    dfz = f ( z, 1 )

    if ( dfz == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u1 = fz / dfz

    if ( u - u1 /= 0.0D+00 ) then
      em = em1 * u / ( u - u1 )
    else
      em = em1
    end if

    em = max ( em, 1.0D+00 )
!
!  Set the increment.
!
    dx = - em1 * fx / dfx - em * fz / dfz
!
!  Remember current data for next step.
!
    em1 = em
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine vandev2 ( x, em, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! VANDEV2 implements the improved Van de Vel iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard King,
!    Improving the van de Vel root-finding method,
!    Computing,
!    Volume 30, 1983, pages 373 - 378.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Output, real ( kind = 8 ) EM, an estimate of the multiplicity of the root.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  real ( kind = 8 ) em
  real ( kind = 8 ) em1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) u1
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  em1 = 1.0D+00
  fx = f ( x, 0 )
  dfx = f ( x, 1 )

  if ( dfx == 0.0D+00 ) then
    ierror = 3
    return
  end if

  u1 = fx / dfx

  x = x - fx / dfx
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = fx / dfx

    if ( u1 - u /= 0.0D+00 ) then
      em = em1 * u1 / ( u1 - u )
    else
      em = em1
    end if

    em = max ( em, 1.0D+00 )
!
!  Set the increment.
!
    dx = - em * fx / dfx
!
!  Remember current data for next step.
!
    em1 = em
    u1 = u
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )

  end do

  return
end
subroutine whittaker ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! WHITTAKER implements the convex acceleeration of Whittaker's method.
!
!  Discussion:
!
!    x(n+1) = x(n) - ( f(x(n)) / f'(x(n)) ) * ( 1 - 0.5 * L(x(n)) )
!
!    where
!
!      L(x) = f(x) * f''(x) / f'(x)**2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = d2fx * fx / dfx**2
!
!  Set the increment.
!
    dx = - ( fx / dfx ) * ( 1.0D+00 - 0.5D+00 * u )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

  end do

  return
end
subroutine whittaker2 ( x, abserr, kmax, f, ierror, k )

!*****************************************************************************80
!
!! WHITTAKER2 implements the double convex acceleeration of Whittaker's method.
!
!  Discussion:
!
!    x(n+1) = x(n) - (1/4) * ( f(x(n)) / f'(x(n)) ) *
!      ( 2 - L(x(n)) + ( 2 + L(x(n)) ) / ( 1 - L(x(n)) + 0.5 * L(x(n))**2 ) )
!
!    where
!
!      L(x) = f(x) * f''(x) / f'(x)**2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X.
!    On input, an estimate for the root of the equation.
!    On output, if IERROR = 0, X is an approximate root for which
!    abs ( F(X) ) <= ABSERR.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) d2fx
  real ( kind = 8 ) dfx
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  real ( kind = 8 ) u
  real ( kind = 8 ) x
!
!  Initialization.
!
  ierror = 0
  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    if ( dfx == 0.0D+00 ) then
      ierror = 3
      return
    end if

    u = d2fx * fx / dfx**2

    if ( 1.0D+00 - u + 0.5D+00 * u * u  == 0.0D+00 ) then
      ierror = 4
      return
    end if
!
!  Set the increment.
!
    dx = - 0.25D+00 * ( fx / dfx ) &
      * ( 2.0D+00 - u + ( 2.0D+00 + u ) / ( 1.0D+00 - u + 0.5D+00 * u * u ) )
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

  end do

  return
end
subroutine zoomin ( a, b, c, abserr, kmax, f, ipoly, mult, nder, nsub )

!*****************************************************************************80
!
!! ZOOMIN calls all the zero finders for a given problem.
!
!  Discussion:
!
!    The original code was written by Harold Deiss in BASIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, three estimates for a root of the
!    function.  The bisection routines will only be called if the function
!    evaluated at A, B and C includes both positive and negative values.
!
!    Input, real ( kind = 8 ) ABSERR, an error tolerance.
!
!    Input, real ( kind = 8 ) external F, the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!    Input, integer ( kind = 4 ) IPOLY, is 0 if the function is not known to be
!    a polynomial, or is equal to the polynomial degree otherwise.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum number of iterations allowed.
!
!    Input, integer ( kind = 4 ) MULT, the multiplicity of the root being sought.
!    If the multiplicity of the root is not known, set MULT to 1.
!
!    Input, integer ( kind = 4 ) NDER, the highest order derivative that can be
!    computed by the user function.  NDER = 0 means only the function
!    value itself is available, for instance.  NDER must be at least 0,
!    and no algorithm needs a derivative higher than 3.
!
!    Input, integer ( kind = 4 ) NSUB, the number of substeps to take.  A few
!    algorithms include a "sub-iteration".  For instance, the modified
!    Newton method evaluates the derivative function only before
!    steps 1, NSUB+1, 2*NSUB+1 and so on, and uses that derivative
!    value for NSUB iterations in a row.  This option is useful when
!    the derivative evaluation is expensive.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  logical change
  character ( len = 3 ) ctemp
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) em
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fc
  real ( kind = 8 ) ff
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) nder
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ipoly
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) mult
  character ( len = 80 ) name
  integer ( kind = 4 ) nsub
  real ( kind = 8 ) rho
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) sc
  real ( kind = 8 ) temp
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZOOMIN'
  write ( *, '(a)' ) '  A compilation of scalar zero finders,'
  write ( *, '(a)' ) '  based on the work of Joseph Traub.'
  write ( *, '(a)' ) ' '

  fa = f ( a, 0 )
  fb = f ( b, 0 )
  fc = f ( c, 0 )
!
!  Rearrange the data, if necessary, so that the pair (A,FA) has
!  the smallest function value of the three sets of data.
!
  if ( abs ( fb ) < abs ( fa ) ) then
    call r8_swap ( a, b )
    call r8_swap ( fa, fb )
  end if

  if ( abs ( fc ) < abs ( fa ) ) then
    call r8_swap ( a, c )
    call r8_swap ( fa, fc )
  end if
!
!  If necessary, switch B and C, so that FB is of opposite sign to FA,
!  if possible.
!
  sa = sign ( 1.0D+00, fa )
  sb = sign ( 1.0D+00, fb )
  sc = sign ( 1.0D+00, fc )

  if ( sa /= sb ) then
    change = .true.
  else if ( sa /= sc ) then
    call r8_swap ( b, c )
    call r8_swap ( fb, fc )
    change = .true.
  else
    change = .false.
  end if
!
!  Check other input.
!
  if ( mult < 1 ) then
    mult = 1
  end if

  if ( nder < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ZOOMIN - Fatal error!'
    write ( *, '(a,i4)' ) '  The input quantity NDER = ', nder
    write ( *, '(a)' ) '  But NDER must be at least 0.'
    return
  end if

  if ( nsub < 1 ) then
    nsub = 1
  end if
!
!  Report start up conditions.
!
  write ( *, '(a)' ) '  1 point formulas use:'
  write ( *, '(a,g14.6)' ) '    x1 = ', a
  write ( *, '(a,g14.6)' ) '    fx1=', fa
  write ( * , '(a)' ) '  2 point formulas add:'
  write ( *, '(a,g14.6)' ) '    x2 = ', b
  write ( *, '(a,g14.6)' ) '    fx2=', fb
  write ( * , '(a)' ) '  3 point formulas add:'
  write ( *, '(a,g14.6)' ) '    x3 = ', c
  write ( *, '(a,g14.6)' ) '    fx3=', fc
  write ( * , '(a)' ) ' '
  write ( * , '(a,i4)' ) '  User estimated multiplicity = ', mult
  if ( 0 <= ipoly ) then
    write ( * , '(a,i4)' ) '  Polynomial degree = ', ipoly
  else
    write ( *, '(a)' ) '  The function is not known to be polynomial.'
  end if
  write ( * , '(a,i4)' ) '  Highest derivative supplied = ', nder
  write ( * , '(a,g14.6)' ) '  Error tolerance = ', abserr
  write ( * , '(a,i4)' ) '  Maximum number of steps = ', kmax
  write ( * , '(a,i4)' ) '  Newton method substep parameter  = ', nsub
  write ( * , '(a)' ) ' '
  write ( * , '(a)' ) ' '
  write ( * , '(a)' ) '  Technique                           Root        Steps' &
     // '   Error     Multiplicity'

  write ( * , '(a)' ) ' '
  write ( * , '(a)' ) '  1.  One point iteration functions with memory:'
  write ( * , '(a)' ) ' '
!
!  Secant algorithm.
!
  name = 'Secant'
  x = a
  x1 = b
  call secant ( x, x1, abserr, kmax, f, ierror, k )

  if ( ierror /= 0 ) then
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
  else
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
  end if
!
!  Extended secant algorithm.
!
  name = 'Extended secant'
  x = a
  x1 = b
  x2 = c
  call secant_extended ( x, x1, x2, abserr, kmax, f, ierror, k )

  if ( ierror /= 0 ) then
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
  else
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
  end if
!
!  Capital Phi(2,1).
!
  if ( ipoly < 0 ) then
    write ( *, '(a)' ) '  Capital Phi(2,1) is misbehaving.'
  else
    name = 'Capital Phi(2,1)'
    x = a
    x1 = b
    x2 = c
    call cap_phi_21 ( x, x1, x2, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if
  end if
!
!  Muller's algorithm.
!
  name = 'R8_Muller'
  x = a
  x1 = b
  x2 = c
  call r8_muller ( x, x1, x2, abserr, kmax, f, ierror, k )

  if ( ierror /= 0 ) then
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
  else
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
  end if
!
!  Perp E(2,1).
!
  if ( 1 <= nder ) then

    name = 'Perp E(2,1)'
    x = a
    x1 = b
    x2 = c
    call perp_e_21 ( x, x1, x2, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Star E(2,1).
!
  name = 'Star E(2,1)'
  x = a
  x1 = b
  x2 = c
  call star_e21 ( x, x1, x2, abserr, kmax, f, ierror, k )

  if ( ierror /= 0 ) then
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
  else
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
  end if
!
!  Finite difference Halley's method.
!
  name = 'Finite difference Halley'
  x = a
  x1 = b
  x2 = c
  call halley2 ( x, x1, x2, abserr, kmax, f, ierror, k )

  if ( ierror /= 0 ) then
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
  else
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
  end if
!
!  Phi(1,2).
!
  if ( 1 <= nder ) then

    name = 'Phi(1,2)'
    x = a
    x1 = b
    call phi_12 ( x, x1, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Perp E(1,2).
!
  if ( 1 <= nder ) then

    name = 'Perp E(1,2)'
    x = a
    x1 = b
    call perp_e_12 ( x, x1, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Star E(1,2).
!
  if ( 1 <= nder ) then

    name = 'Star E(1,2)'
    x = a
    x1 = b
    call star_e12 ( x, x1, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Dagger E(1,2).
!
  if ( 1 <= nder ) then

    name = 'Dagger E(1,2)'
    x = a
    x1 = b
    call dagger_e12 ( x, x1, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  2. One point iteration functions.'
  write ( *, '(a)' ) ' '
!
!  Newton method.
!
  if ( 1 <= nder ) then

    name = 'Newton'
    x = a
    call newton ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Steffenson's method.
!
  name = 'Steffenson'
  x = a
  call steffenson ( x, abserr, kmax, f, ierror, k )

  if ( ierror /= 0 ) then
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
  else
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
  end if
!
!  Stirling's method.
!
  if ( 1 <= nder ) then

    name = 'Stirling'
    x = a
    call stirling ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Midpoint method.
!
  if ( 1 <= nder ) then

    name = 'midpoint'
    x = a
    call midpoint ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Traub-Ostrowski method.
!
  name = 'Traub-Ostrowski'
  x = a
  call traub_ostrowski ( x, abserr, kmax, f, ierror, k )

  if ( ierror /= 0 ) then
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
  else
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
  end if
!
!  Chebyshev method.
!
  if ( 2 <= nder ) then

    name = 'Chebyshev'
    x = a
    call chebyshev ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Halley Super method.
!
  if ( 2 <= nder ) then

    name = 'Halley Super'
    x = a
    call halley_super ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Whittaker convex acceleration method.
!
  if ( 2 <= nder ) then

    name = 'Whittaker'
    x = a
    call whittaker ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Whittaker double convex acceleration method.
!
  if ( 2 <= nder ) then

    name = 'Whittaker2'
    x = a
    call whittaker2 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  E3.
!
  if ( 2 <= nder ) then

    name = 'E3'
    x = a
    call e3 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  E4.
!
  if ( 3 <= nder ) then

    name = 'E4'
    x = a
    call e4 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Halley's method.
!
  if ( 2 <= nder ) then

    name = 'Halley'
    x = a
    call halley1 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  PSI(2,1).
!
  if ( 3 <= nder ) then

    name = 'Psi(2,1)'
    x = a
    call psi_21 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  PSI(1,2).
!
  if ( 3 <= nder ) then

    name = 'Psi(1,2)'
    x = a
    call psi_i2 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Capital PHI(0,3).
!
  if ( 2 <= nder ) then

    name = 'Capital Phi(0,3)'
    x = a
    call cap_phi_03 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Reduced capital PHI(0,4).
!
  if ( 3 <= nder ) then

    name = 'Reduced Capital Phi(0,4)'
    x = a
    call red_cap_phi_04 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Functions from Hansen and Patrick.
!
!  The Ostrowski square root formula is the Hansen family member with BETA = 0.
!
!     if ( ipoly < 0 ) then
!       write ( *, '(a)' ) '  Ostrowski square root is misbehaving.'
!     else
    if ( 2 <= nder ) then
      name = 'Ostrowski square root'
      x = a
      call ostrowski_sqrt ( x, abserr, kmax, f, ierror, k )

      if ( ierror /= 0 ) then
        write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
      else
        write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
      end if

    end if
!     end if
!
!  The Euler method.
!
  if ( ipoly < 0 ) then

    write ( *, '(a)' ) '  Euler is misbehaving.'

  else

    if ( 2 <= nder ) then

      name = 'Euler'
      x = a
      call euler ( x, abserr, kmax, f, ierror, k )

      if ( ierror /= 0 ) then
        write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
      else
        write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
      end if

    end if

  end if
!
!  The Laguerre method.
!
  if ( 2 <= nder .and. 2 <= ipoly ) then

    name = 'Laguerre'
    x = a
    call laguerre ( x, ipoly, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  3.  Multipoint iteration functions.'
  write ( *, '(a)' ) ' '
!
!  First function in Traub, page 236.
!
  if ( 1 <= nder ) then

    name = 'Traub first'
    x = a
    call traub1 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  A family of iteration functions including
!  forms 2, 12, and 13 of Traub, pages 236 - 237.
!
  if ( 1 <= nder ) then

    name = 'Traub second'
    e = 0.5D+00
    d = 1.0D+00
    x = a
    call t_family1 ( x, e, d, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if

  if ( 1 <= nder ) then

    name = 'Traub twelfth'
    e = 0.25D+00
    d = 2.0D+00 / 3.0D+00
    x = a
    call t_family1 ( x, e, d, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if

  if ( 1 <= nder ) then

    name = 'Traub thirteenth'
    e = 5.0D+00 / 12.0D+00
    d = 6.0D+00 / 7.0D+00
    x = a
    call t_family1 ( x, e, d, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Modified Newton method.
!
  if ( 1 <= nder .and. 1 <= nsub ) then

    write ( name, '(a,i3)' ) 'Modified Newton, NSUB=', nsub
    x = a
    call newton_mod ( x, nsub, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Fourth function in Traub.
!
  if ( 1 <= nder .and. 1 <= nsub ) then

    name = 'Traub fourth, NSUB='
    write (ctemp,'(i3)') nsub
    name(20:22) = ctemp(1:3)
    x = a
    call traub4 ( x, nsub, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Newton - secant.
!
  if ( 1 <= nder ) then

    name = 'Newton - secant'
    x = a
    call newton_sec ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  A family of iteration functions including
!  forms 6 and 7 of Traub, pages 236 - 237.
!
  if ( 1 <= nder ) then

    name = 'Traub sixth'
    e = 2.0D+00
    ff = 3.0D+00
    g = 1.0D+00
    h = 1.0D+00
    x = a
    call t_family2 ( x, e, ff, g, h, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if

  if ( 1 <= nder ) then

    name = 'Traub seventh'
    e = 4.0D+00
    ff = 7.0D+00
    g = 3.0D+00
    h = 2.0D+00 / 3.0D+00
    x = a
    call t_family2 ( x, e, ff, g, h, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Eighth function in Traub, page 237.
!
  if ( 1 <= nder ) then

    name = 'Traub eighth'
    x = a
    call traub8 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Ninth function in Traub, page 237.
!
  if ( 1 <= nder ) then

    name = 'Traub ninth'
    x = a
    call traub9 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Type 1 family of iteration functions including
!  forms 10 and 11 of Traub, pages 237.
!
  if ( 1 <= nder ) then

    name = 'Traub type 1, form 10'
    rho = ( 1.0D+00 - sqrt ( 5.0D+00 ) ) / 2.0D+00
    e = 0.0D+00
    x = a
    call tt1f ( x, e, rho, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if

  if ( 1 <= nder ) then

    name = 'Traub type 1, form 11'
    e = 1.0D+00
    x = a
    call tt1f ( x, e, rho, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Fourteenth function in Traub, page 237
!
  if ( 1 <= nder ) then

    name = 'Traub fourteenth'
    x = a
    call t14 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Fifteenth function in Traub, page 238.
!
  if ( 1 <= nder ) then

    name = 'Traub fifteenth'
    x = a
    call t15 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Sixteenth function in Traub, page 238.
!
  if ( 1 <= nder ) then

    name = 'Traub sixteenth'
    x = a
    call t16 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  A family of fourth order methods by Richard King.
!
  if ( 1 <= nder ) then

    name = 'King, BETA=0'
    beta = 0.0D+00
    x = a
    call king ( x, beta, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if

  if ( 1 <= nder ) then

    name = 'King, BETA=1'
    beta = 1.0D+00
    x = a
    call king ( x, beta, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if

  if ( 1 <= nder ) then

    name = 'King, BETA=2'
    beta = 2.0D+00
    x = a
    call king ( x, beta, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Jarratt iterative technique.
!
  if ( 1 <= nder ) then

    name = 'Jarratt'
    x = a
    call jarratt ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Jarratt inverse-free iterative technique.
!
  if ( 1 <= nder ) then

    name = 'Jarratt inverse-free'
    x = a
    call jarratt2 ( x, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  4.  Multiple root methods, multiplicity given.'
  write ( *, '(a)' ) ' '
!
!  Script E2.
!
  if ( 1 <= nder ) then

    name = 'Script E2'
    x = a
    call script_e2 ( x, mult, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Script E3.
!
  if ( 2 <= nder ) then

    name = 'Script E3'
    x = a
    call script_e3 ( x, mult, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Script E4.
!
  if ( 3 <= nder ) then

    name = 'Script E4'
    x = a
    call script_e4 ( x, mult, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if
!
!  Star E - 1,1(f) function in Traub, page 235.
!
  name = 'Star E 1,1(f)'
  x = a
  x1 = b
  call te11f ( x, x1, mult, abserr, kmax, f, ierror, k )

  if ( ierror /= 0 ) then
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
  else
    write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
  end if
!
!  Routines which determine the root and multiplicity
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  5.  Multiple root methods,'
  write ( *, '(a)' ) '      multiplicity not given.'
  write ( *, '(a)' ) ' '
!
!  E2(U).
!
  if ( 2 <= nder ) then

    name = 'E2(U)'
    x = a
    call te2u ( x, em, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, 5x, g15.6 )' ) name, x, k, em
    end if

  end if
!
!  Phi - 1,1(u) function in Traub, page 235.
!
  if ( 1 <= nder ) then

    name = 'Phi 1,1(U)'
    x = a
    call tphi1u ( x, em, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, 5x, g15.6 )' ) name, x, k, em
    end if

  end if
!
!  Third function in Traub, page 235 bottom.
!
  if ( 1 <= nder ) then

    name = 'Traub third'
    x = a
    call tthip ( x, em, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, 5x, g15.6 )' ) name, x, k, em
    end if

  end if
!
!  Van de Vel iteration.
!
  if ( 1 <= nder ) then
    name = 'Van de Vel'
    x = a
    call vandev ( x, em, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, 5x, g15.6 )' ) name, x, k, em
    end if

  end if
!
!  Improved Van de Vel iteration.
!
  if ( 1 <= nder ) then
    name = 'Improved Van de Vel'
    x = a
    call vandev2 ( x, em, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, 5x, g15.6 )' ) name, x, k, em
    end if

  end if
!
!  Bisection methods
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  6. Bisection methods'
  write ( *, '(a)' ) ' '
!
!  Bisection methods may only be used if there is a change of sign.
!
  if ( .not. change ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ZOOMIN - Note:'
    write ( *, '(a)' ) '  Bisection methods will not be called since'
    write ( *, '(a)' ) '  there is no change of sign interval.'

  else
!
!  Bisection.
!
    name = 'Bisection'
    x = a
    x1 = b
    call bisect ( x, x1, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if
!
!  Regula falsi.
!
    name = 'Regula falsi'
    x = a
    x1 = b
    call regula ( x, x1, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if
!
!  Brent.
!
    name = 'Brent'
    x = a
    x1 = b
    call brent ( x, x1, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if
!
!  Bisection + secant.
!
    name = 'Bisection + secant'
    x = a
    x1 = b
    call rhein1 ( x, x1, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if
!
!  Bisection + secant + inverse quadratic.
!
    name = 'Bisection + secant + inv quad'
    x = a
    x1 = b
    call rhein2 ( x, x1, abserr, kmax, f, ierror, k )

    if ( ierror /= 0 ) then
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k, ierror
    else
      write ( *, '(2x,a30,g15.6, i8, i5)' ) name, x, k
    end if

  end if

  return
end
