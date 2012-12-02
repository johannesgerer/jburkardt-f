subroutine bisection ( fatol, step_max, prob, xatol, xa, xb, fxa, fxb )

!*****************************************************************************80
!
!! BISECTION carries out the bisection method to seek a root of F(X) = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) FATOL, an absolute error tolerance for
!    the function value of the root.  If an approximate root X satisfies
!      ABS ( F ( X ) ) <= FATOL, then X will be accepted as the
!    root and the iteration will be terminated.
!
!    Input, integer ( kind = 4 ) STEP_MAX, the maximum number of steps
!    allowed for an iteration.
!
!    Input, integer ( kind = 4 ) PROB, the index of the function whose root is
!    to be sought.
!
!    Input, real ( kind = 8 ) XATOL, an absolute error tolerance for the root.
!
!    Input/output, real ( kind = 8 ) XA, XB, two points at which the
!    function differs in sign.  On output, these values have been adjusted
!    to a smaller interval.
!
!    Input/output, real ( kind = 8 ) FXA, FXB, the value of the function
!    at XA and XB.
!
  implicit none

  real ( kind = 8 ) fatol
  real ( kind = 8 ) fxa
  real ( kind = 8 ) fxb
  real ( kind = 8 ) fxc
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) step_max
  integer ( kind = 4 ) step_num
  real ( kind = 8 ) t
  real ( kind = 8 ) xa
  real ( kind = 8 ) xatol
  real ( kind = 8 ) xb
  real ( kind = 8 ) xc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BISECTION'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Step      XA            XB             F(XA)         F(XB)'
  write ( *, '(a)' ) ' '
!
!  Make A the root with negative F, B the root with positive F.
!
  if ( 0.0D+00 < fxa ) then
    t = xa
    xa = xb
    xb = t
    t = fxa
    fxa = fxb
    fxb = t
  end if

  step_num = 0
!
!  Loop
!
  do

    write ( *, '(2x,i4,2x,2g16.8,2g14.6)' ) step_num, xa, xb, fxa, fxb

    step_num = step_num + 1

    if ( step_max < step_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Maximum number of steps taken without convergence.'
      exit
    end if

    if ( abs ( xa - xb ) < xatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Interval small enough for convergence.'
      exit
    end if

    if ( abs ( fxa ) <= fatol .or. abs ( fxb ) <= fatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Function small enough for convergence.'
      exit
    end if
!
!  Compute the next iterate.
!
    xc = 0.5D+00 * ( xa + xb )
    call p00_fx ( prob, xc, fxc )
!
!  Replace one of the old points.
!
    if ( fxc < 0.0D+00 ) then
      xa = xc
      fxa = fxc
    else
      xb = xc
      fxb = fxc
    end if

  end do

  return
end
subroutine brent ( fatol, step_max, prob, xatol, xrtol, xa, xb, fxa, fxb )

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
!    17 January 2001
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
!  Parameters:
!
!    Input, real ( kind = 8 ) FATOL, an absolute error tolerance for the
!    function value of the root.  If an approximate root X satisfies
!      ABS ( F ( X ) ) <= FATOL, then X will be accepted as the
!    root and the iteration will be terminated.
!
!    Input, integer ( kind = 4 ) STEP_MAX, the maximum number of steps allowed
!    for an iteration.
!
!    Input, integer ( kind = 4 ) PROB, the index of the function whose root is
!    to be sought.
!
!    Input, real ( kind = 8 ) XATOL, XRTOL, absolute and relative error
!    tolerances for the root.
!
!    Input/output, real ( kind = 8 ) XA, XB, two points at which the
!    function differs in sign.  On output, these values have been adjusted
!    to a smaller interval.
!
!    Input/output, real ( kind = 8 ) FXA, FXB, the value of the function
!    at XA and XB.
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) fatol
  real ( kind = 8 ) fxa
  real ( kind = 8 ) fxb
  real ( kind = 8 ) fxc
  integer ( kind = 4 ) step_max
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) step_num
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
  real ( kind = 8 ) xc
  real ( kind = 8 ) xm
  real ( kind = 8 ) xatol
  real ( kind = 8 ) xrtol
  real ( kind = 8 ) xtol
!
!  Initialization.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BRENT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Step      XA            XB             F(XA)         F(XB)'
  write ( *, '(a)' ) ' '

  step_num = 0

  call p00_fx ( prob, xa, fxa )
  call p00_fx ( prob, xb, fxb )
!
!  Check that f(ax) and f(bx) have different signs.
!
  if ( ( fxa < 0.0D+00 .and. fxb < 0.0D+00 ) .or. &
       ( 0.0D+00 < fxa .and. 0.0D+00 < fxb ) ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BRENT - Fatal error!'
    write ( *, '(a)' ) '  F(XA) and F(XB) have same sign.'
    return

  end if

  xc = xa
  fxc = fxa
  d = xb - xa
  e = d

  do

    write ( *, '(2x,i4,2x,2g16.8,2g14.6)' ) step_num, xb, xc, fxb, fxc

    step_num = step_num + 1

    if ( step_max < step_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Maximum number of steps taken.'
      exit
    end if

    if ( abs ( fxc ) < abs ( fxb ) ) then
      xa = xb
      xb = xc
      xc = xa
      fxa = fxb
      fxb = fxc
      fxc = fxa
    end if

    xtol = 2.0D+00 * xrtol * abs ( xb ) + 0.5D+00 * xatol
!
!  XM is the halfwidth of the current change-of-sign interval.
!
    xm = 0.5D+00 * ( xc - xb )

    if ( abs ( xm ) <= xtol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Interval small enough for convergence.'
      exit
    end if

    if ( abs ( fxb ) <= fatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Function small enough for convergence.'
      exit
    end if
!
!  See if a bisection is forced.
!
    if ( abs ( e ) < xtol .or. abs ( fxa ) <= abs ( fxb ) ) then

      d = xm
      e = d

    else

      s = fxb / fxa
!
!  Linear interpolation.
!
      if ( xa == xc ) then

        p = 2.0D+00 * xm * s
        q = 1.0D+00 - s
!
!  Inverse quadratic interpolation.
!
      else

        q = fxa / fxc
        r = fxb / fxc
        p = s * ( 2.0D+00 * xm * q * ( q - r ) - ( xb - xa ) * ( r - 1.0D+00 ) )
        q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

      end if

      if ( 0.0D+00 < p ) then
        q = - q
      else
        p = - p
      end if

      s = e
      e = d

      if ( 3.0D+00 * xm * q - abs ( xtol * q ) <= 2.0D+00 * p .or. &
           abs ( 0.5D+00 * s * q ) <= p ) then
        d = xm
        e = d
      else
        d = p / q
      end if

    end if
!
!  Save in XA, FXA the previous values of XB, FXB.
!
    xa = xb
    fxa = fxb
!
!  Compute the new value of XB, and evaluate the function there.
!
    if ( xtol < abs ( d ) ) then
      xb = xb + d
    else if ( 0.0D+00 < xm ) then
      xb = xb + xtol
    else if ( xm <= 0.0D+00 ) then
      xb = xb - xtol
    end if

    call p00_fx ( prob, xb, fxb )
!
!  If the new FXB has the same sign as FXC, then replace XC by XA.
!
    if ( ( 0.0D+00 < fxb .and. 0.0D+00 < fxc ) .or. &
         ( fxb < 0.0D+00 .and. fxc < 0.0D+00 ) ) then
      xc = xa
      fxc = fxa
      d = xb - xa
      e = d
    end if

  end do

  return
end
subroutine muller ( fatol, step_max, prob, xatol, xrtol, xa, xb, xc, fxa, &
  fxb, fxc )

!*****************************************************************************80
!
!! MULLER carries out Muller's method for a real root of a nonlinear function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) FATOL, an absolute error tolerance for the
!    function value of the root.  If an approximate root X satisfies
!      ABS ( F ( X ) ) <= FATOL, then X will be accepted as the
!    root and the iteration will be terminated.
!
!    Input, integer ( kind = 4 ) STEP_MAX, the maximum number of steps allowed
!    for an iteration.
!
!    Input, integer ( kind = 4 ) PROB, the index of the function whose root is
!    to be sought.
!
!    Input, real ( kind = 8 ) XATOL, XRTOL, absolute and relative error
!    tolerances  for the root.
!
!    Input/output, real ( kind = 8 ) XA, XB, XC, three points.
!
!    Input/output, real ( kind = 8 ) FXA, FXB, FXC, the value of the
!    function at XA, XB, and XC.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) fatol
  integer ( kind = 4 ) i
  integer ( kind = 4 ) step_max
  integer ( kind = 4 ) prob
  real ( kind = 8 ) t
  real ( kind = 8 ) xa
  real ( kind = 8 ) xatol
  real ( kind = 8 ) xb
  real ( kind = 8 ) xc
  real ( kind = 8 ) xd
  real ( kind = 8 ) xrtol
  real ( kind = 8 ) fxa
  real ( kind = 8 ) fxb
  real ( kind = 8 ) fxc
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2
!
!  Initialization.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MULLER'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step      XA           XB           XC'
  write ( *, '(a)' ) '          F(XA)        F(XB)        F(XC)'
  write ( *, '(a)' ) ' '

  i = 0

  write ( *, '(2x,i4,3g14.6)' ) i,  xa,  xb,  xc
  write ( *, '(2x,4x,3g14.6)' )    fxa, fxb, fxc

  do i = 1, step_max
!
!  Determine the coefficients
!    A, B, C
!  of the polynomial
!    Y(X) = A * (X-X2)**2 + B * (X-X2) + C
!  which goes through the data:
!    (X1,Y1), (X2,Y2), (X3,Y3).
!
    a = ( ( fxa - fxc ) * ( xb - xc ) &
        - ( fxb - fxc ) * ( xa - xc ) ) / &
          ( ( xa - xc ) * ( xb - xc ) * ( xa - xb ) )

    b = ( ( fxb - fxc ) * ( xa - xc )**2 &
        - ( fxa - fxc ) * ( xb - xc )**2 ) / &
        ( ( xa - xc ) * ( xb - xc ) * ( xa - xb ) )

    c = fxc
!
!  Get the real roots of the polynomial,
!  unless A = 0, in which case the algorithm is breaking down.
!
    if ( a /= 0.0D+00 ) then

      call r8poly2_rroot ( a, b, c, z1, z2 )

    else if ( b /= 0.0D+00 ) then

      z2 = - c / b

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Polynomial fitting has failed.'
      write ( *, '(a)' ) '  Muller''s algorithm breaks down.'
      return

    end if

    xd = xc + z2
!
!  Set XA, YA, based on which of XA and XB is closer to XD.
!
    if ( abs ( xd - xb ) < abs ( xd - xa ) ) then
      t = xa
      xa = xb
      xb = t
      t = fxa
      fxa = fxb
      fxb = t
    end if
!
!  Set XB, YB, based on which of XB and XC is closer to XD.
!
    if ( abs ( xd - xc ) < abs ( xd - xb ) ) then
      t = xb
      xb = xc
      xc = t
      t = fxb
      fxb = fxc
      fxc = t
    end if
!
!  Set XC, YC.
!
    xc = xd
    call p00_fx ( prob, xc, fxc )

    write ( *, '(2x,i4,3g14.6)' ) i,  xa,  xb,  xc
    write ( *, '(2x,4x,3g14.6)' )    fxa, fxb, fxc
!
!  Estimate the relative significance of the most recent correction.
!
    if ( abs ( z2 ) <= xrtol * abs ( xc ) + xatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Stepsize small enough for convergence.'
      return
    end if

    if ( abs ( fxc ) < fatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Function small enough for convergence.'
      return
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Took maximum number of steps without convergence.'

  return
end
subroutine newton ( fatol, step_max, prob, xatol, xmin, xmax, xa, fxa )

!*****************************************************************************80
!
!! NEWTON carries out Newton's method to seek a root of F(X) = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) FATOL, an absolute error tolerance for the
!    function value of the root.  If an approximate root X satisfies
!      ABS ( F ( X ) ) <= FATOL, then X will be accepted as the
!    root and the iteration will be terminated.
!
!    Input, integer ( kind = 4 ) STEP_MAX, the maximum number of steps allowed
!    for an iteration.
!
!    Input, integer ( kind = 4 ) PROB, the index of the function whose root is
!    to be sought.
!
!    Input, real ( kind = 8 ) XATOL, an absolute error tolerance for the root.
!
!    Input, real ( kind = 8 ) XMIN, XMAX, the interval in which the root should
!    be sought.
!
!    Input/output, real ( kind = 8 ) XA.  On input, the starting point for
!    the iteration.  On output, the current approximation to the root.
!
!    Input/output, real ( kind = 8 ) FXA, the function value at XA.
!
  implicit none

  real ( kind = 8 ) fatol
  real ( kind = 8 ) fp
  real ( kind = 8 ) fxa
  integer ( kind = 4 ) step_max
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) step_num
  real ( kind = 8 ) step
  real ( kind = 8 ) xa
  real ( kind = 8 ) xatol
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  step = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NEWTON'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X             F(X)      FP(X)'
  write ( *, '(a)' ) ' '

  step_num = 0
  call p00_fx1 ( prob, xa, fp )

  write ( *, '(2x,i4,2x,g16.8,2g14.6)' ) step_num, xa, fxa, fp

  do step_num = 1, step_max

    if ( xa < xmin .or. xmax < xa ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The iterate X = ', xa
      write ( *, '(a)' ) '  has left the region [XMIN,XMAX].'
      return
    end if

    if ( abs ( fxa ) <= fatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The function norm is small enough for convergence.'
      return
    end if

    if ( 1 < step_num .and. abs ( step ) <= xatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The stepsize is small enough for convergence.'
      return
    end if

    if ( fp == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F''(X)=0, the algorithm fails.'
      return
    end if

    step = fxa / fp

    xa = xa - step

    call p00_fx ( prob, xa, fxa )
    call p00_fx1 ( prob, xa, fp )

    write ( *, '(2x,i4,2x,g16.8,2g14.6)' ) step_num, xa, fxa, fp

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Took maximum number of steps without convergence.'

  return
end
subroutine p00_fx ( prob, x, fx )

!*****************************************************************************80
!
!! P00_FX evaluates a function specified by problem number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the problem.
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) prob
  real ( kind = 8 ) x

  if ( prob == 1 ) then
    call p01_fx ( x, fx )
  else if ( prob == 2 ) then
    call p02_fx ( x, fx )
  else if ( prob == 3 ) then
    call p03_fx ( x, fx )
  else if ( prob == 4 ) then
    call p04_fx ( x, fx )
  else if ( prob == 5 ) then
    call p05_fx ( x, fx )
  else if ( prob == 6 ) then
    call p06_fx ( x, fx )
  else if ( prob == 7 ) then
    call p07_fx ( x, fx )
  else if ( prob == 8 ) then
    call p08_fx ( x, fx )
  else if ( prob == 9 ) then
    call p09_fx ( x, fx )
  else if ( prob == 10 ) then
    call p10_fx ( x, fx )
  else if ( prob == 11 ) then
    call p11_fx ( x, fx )
  else if ( prob == 12 ) then
    call p12_fx ( x, fx )
  else if ( prob == 13 ) then
    call p13_fx ( x, fx )
  else if ( prob == 14 ) then
    call p14_fx ( x, fx )
  else if ( prob == 15 ) then
    call p15_fx ( x, fx )
  else if ( prob == 16 ) then
    call p16_fx ( x, fx )
  else if ( prob == 17 ) then
    call p17_fx ( x, fx )
  else if ( prob == 18 ) then
    call p18_fx ( x, fx )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FX - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_fx1 ( prob, x, fx1 )

!*****************************************************************************80
!
!! P00_FX1: first derivative of a function specified by problem number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the problem.
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none
!
  real ( kind = 8 ) fx1
  integer ( kind = 4 ) prob
  real ( kind = 8 ) x

  if ( prob == 1 ) then
    call p01_fx1 ( x, fx1 )
  else if ( prob == 2 ) then
    call p02_fx1 ( x, fx1 )
  else if ( prob == 3 ) then
    call p03_fx1 ( x, fx1 )
  else if ( prob == 4 ) then
    call p04_fx1 ( x, fx1 )
  else if ( prob == 5 ) then
    call p05_fx1 ( x, fx1 )
  else if ( prob == 6 ) then
    call p06_fx1 ( x, fx1 )
  else if ( prob == 7 ) then
    call p07_fx1 ( x, fx1 )
  else if ( prob == 8 ) then
    call p08_fx1 ( x, fx1 )
  else if ( prob == 9 ) then
    call p09_fx1 ( x, fx1 )
  else if ( prob == 10 ) then
    call p10_fx1 ( x, fx1 )
  else if ( prob == 11 ) then
    call p11_fx1 ( x, fx1 )
  else if ( prob == 12 ) then
    call p12_fx1 ( x, fx1 )
  else if ( prob == 13 ) then
    call p13_fx1 ( x, fx1 )
  else if ( prob == 14 ) then
    call p14_fx1 ( x, fx1 )
  else if ( prob == 15 ) then
    call p15_fx1 ( x, fx1 )
  else if ( prob == 16 ) then
    call p16_fx1 ( x, fx1 )
  else if ( prob == 17 ) then
    call p17_fx1 ( x, fx1 )
  else if ( prob == 18 ) then
    call p18_fx1 ( x, fx1 )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FX1 - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_fx2 ( prob, x, fx2 )

!*****************************************************************************80
!
!! P00_FX2: second derivative of a function specified by problem number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the problem.
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  integer ( kind = 4 ) prob
  real ( kind = 8 ) x

  if ( prob == 1 ) then
    call p01_fx2 ( x, fx2 )
  else if ( prob == 2 ) then
    call p02_fx2 ( x, fx2 )
  else if ( prob == 3 ) then
    call p03_fx2 ( x, fx2 )
  else if ( prob == 4 ) then
    call p04_fx2 ( x, fx2 )
  else if ( prob == 5 ) then
    call p05_fx2 ( x, fx2 )
  else if ( prob == 6 ) then
    call p06_fx2 ( x, fx2 )
  else if ( prob == 7 ) then
    call p07_fx2 ( x, fx2 )
  else if ( prob == 8 ) then
    call p08_fx2 ( x, fx2 )
  else if ( prob == 9 ) then
    call p09_fx2 ( x, fx2 )
  else if ( prob == 10 ) then
    call p10_fx2 ( x, fx2 )
  else if ( prob == 11 ) then
    call p11_fx2 ( x, fx2 )
  else if ( prob == 12 ) then
    call p12_fx2 ( x, fx2 )
  else if ( prob == 13 ) then
    call p13_fx2 ( x, fx2 )
  else if ( prob == 14 ) then
    call p14_fx2 ( x, fx2 )
  else if ( prob == 15 ) then
    call p15_fx2 ( x, fx2 )
  else if ( prob == 16 ) then
    call p16_fx2 ( x, fx2 )
  else if ( prob == 17 ) then
    call p17_fx2 ( x, fx2 )
  else if ( prob == 18 ) then
    call p18_fx2 ( x, fx2 )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FX2 - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_prob_num ( prob_num )

!*****************************************************************************80
!
!! P00_PROB_NUM returns the number of problems available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) PROB_NUM, the number of problems available.
!
  implicit none

  integer ( kind = 4 ) prob_num

  prob_num = 18

  return
end
subroutine p00_range ( prob, range )

!*****************************************************************************80
!
!! P00_RANGE returns an interval bounding the root for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the problem.
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  integer ( kind = 4 ) prob
  real ( kind = 8 ) range(2)

  if ( prob == 1 ) then
    call p01_range ( range )
  else if ( prob == 2 ) then
    call p02_range ( range )
  else if ( prob == 3 ) then
    call p03_range ( range )
  else if ( prob == 4 ) then
    call p04_range ( range )
  else if ( prob == 5 ) then
    call p05_range ( range )
  else if ( prob == 6 ) then
    call p06_range ( range )
  else if ( prob == 7 ) then
    call p07_range ( range )
  else if ( prob == 8 ) then
    call p08_range ( range )
  else if ( prob == 9 ) then
    call p09_range ( range )
  else if ( prob == 10 ) then
    call p10_range ( range )
  else if ( prob == 11 ) then
    call p11_range ( range )
  else if ( prob == 12 ) then
    call p12_range ( range )
  else if ( prob == 13 ) then
    call p13_range ( range )
  else if ( prob == 14 ) then
    call p14_range ( range )
  else if ( prob == 15 ) then
    call p15_range ( range )
  else if ( prob == 16 ) then
    call p16_range ( range )
  else if ( prob == 17 ) then
    call p17_range ( range )
  else if ( prob == 18 ) then
    call p18_range ( range )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_RANGE - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_root ( prob, i, x )

!*****************************************************************************80
!
!! P00_ROOT returns a known root for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the problem.
!
!    Input, integer ( kind = 4 ) I, the index of the root to return.
!
!    Output, real ( kind = 8 ) X, the I-th root.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) prob
  real ( kind = 8 ) x

  if ( prob == 1 ) then
    call p01_root ( i, x )
  else if ( prob == 2 ) then
    call p02_root ( i, x )
  else if ( prob == 3 ) then
    call p03_root ( i, x )
  else if ( prob == 4 ) then
    call p04_root ( i, x )
  else if ( prob == 5 ) then
    call p05_root ( i, x )
  else if ( prob == 6 ) then
    call p06_root ( i, x )
  else if ( prob == 7 ) then
    call p07_root ( i, x )
  else if ( prob == 8 ) then
    call p08_root ( i, x )
  else if ( prob == 9 ) then
    call p09_root ( i, x )
  else if ( prob == 10 ) then
    call p10_root ( i, x )
  else if ( prob == 11 ) then
    call p11_root ( i, x )
  else if ( prob == 12 ) then
    call p12_root ( i, x )
  else if ( prob == 13 ) then
    call p13_root ( i, x )
  else if ( prob == 14 ) then
    call p14_root ( i, x )
  else if ( prob == 15 ) then
    call p15_root ( i, x )
  else if ( prob == 16 ) then
    call p16_root ( i, x )
  else if ( prob == 17 ) then
    call p17_root ( i, x )
  else if ( prob == 18 ) then
    call p18_root ( i, x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_root_num ( prob, root_num )

!*****************************************************************************80
!
!! P00_ROOT_NUM returns the number of known roots for a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the problem.
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!    This value may be zero.
!
  implicit none

  integer ( kind = 4 ) prob
  integer ( kind = 4 ) root_num

  if ( prob == 1 ) then
    call p01_root_num ( root_num )
  else if ( prob == 2 ) then
    call p02_root_num ( root_num )
  else if ( prob == 3 ) then
    call p03_root_num ( root_num )
  else if ( prob == 4 ) then
    call p04_root_num ( root_num )
  else if ( prob == 5 ) then
    call p05_root_num ( root_num )
  else if ( prob == 6 ) then
    call p06_root_num ( root_num )
  else if ( prob == 7 ) then
    call p07_root_num ( root_num )
  else if ( prob == 8 ) then
    call p08_root_num ( root_num )
  else if ( prob == 9 ) then
    call p09_root_num ( root_num )
  else if ( prob == 10 ) then
    call p10_root_num ( root_num )
  else if ( prob == 11 ) then
    call p11_root_num ( root_num )
  else if ( prob == 12 ) then
    call p12_root_num ( root_num )
  else if ( prob == 13 ) then
    call p13_root_num ( root_num )
  else if ( prob == 14 ) then
    call p14_root_num ( root_num )
  else if ( prob == 15 ) then
    call p15_root_num ( root_num )
  else if ( prob == 16 ) then
    call p16_root_num ( root_num )
  else if ( prob == 17 ) then
    call p17_root_num ( root_num )
  else if ( prob == 18 ) then
    call p18_root_num ( root_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_ROOT_NUM - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_start ( prob, i, x )

!*****************************************************************************80
!
!! P00_START returns a starting point for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the problem.
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) prob
  real ( kind = 8 ) x

  if ( prob == 1 ) then
    call p01_start ( i, x )
  else if ( prob == 2 ) then
    call p02_start ( i, x )
  else if ( prob == 3 ) then
    call p03_start ( i, x )
  else if ( prob == 4 ) then
    call p04_start ( i, x )
  else if ( prob == 5 ) then
    call p05_start ( i, x )
  else if ( prob == 6 ) then
    call p06_start ( i, x )
  else if ( prob == 7 ) then
    call p07_start ( i, x )
  else if ( prob == 8 ) then
    call p08_start ( i, x )
  else if ( prob == 9 ) then
    call p09_start ( i, x )
  else if ( prob == 10 ) then
    call p10_start ( i, x )
  else if ( prob == 11 ) then
    call p11_start ( i, x )
  else if ( prob == 12 ) then
    call p12_start ( i, x )
  else if ( prob == 13 ) then
    call p13_start ( i, x )
  else if ( prob == 14 ) then
    call p14_start ( i, x )
  else if ( prob == 15 ) then
    call p15_start ( i, x )
  else if ( prob == 16 ) then
    call p16_start ( i, x )
  else if ( prob == 17 ) then
    call p17_start ( i, x )
  else if ( prob == 18 ) then
    call p18_start ( i, x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_start_num ( prob, start_num )

!*****************************************************************************80
!
!! P00_START_NUM returns the number of starting points for a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the problem.
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) prob
  integer ( kind = 4 ) start_num

  if ( prob == 1 ) then
    call p01_start_num ( start_num )
  else if ( prob == 2 ) then
    call p02_start_num ( start_num )
  else if ( prob == 3 ) then
    call p03_start_num ( start_num )
  else if ( prob == 4 ) then
    call p04_start_num ( start_num )
  else if ( prob == 5 ) then
    call p05_start_num ( start_num )
  else if ( prob == 6 ) then
    call p06_start_num ( start_num )
  else if ( prob == 7 ) then
    call p07_start_num ( start_num )
  else if ( prob == 8 ) then
    call p08_start_num ( start_num )
  else if ( prob == 9 ) then
    call p09_start_num ( start_num )
  else if ( prob == 10 ) then
    call p10_start_num ( start_num )
  else if ( prob == 11 ) then
    call p11_start_num ( start_num )
  else if ( prob == 12 ) then
    call p12_start_num ( start_num )
  else if ( prob == 13 ) then
    call p13_start_num ( start_num )
  else if ( prob == 14 ) then
    call p14_start_num ( start_num )
  else if ( prob == 15 ) then
    call p15_start_num ( start_num )
  else if ( prob == 16 ) then
    call p16_start_num ( start_num )
  else if ( prob == 17 ) then
    call p17_start_num ( start_num )
  else if ( prob == 18 ) then
    call p18_start_num ( start_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_START_NUM - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_title ( prob, title )

!*****************************************************************************80
!
!! P00_TITLE returns the title for a given problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the problem.
!
!    Output, character ( len = * ) TITLE, the title of the given problem.
!
  implicit none

  integer ( kind = 4 ) prob
  character ( len = * ) title

  if ( prob == 1 ) then
    call p01_title ( title )
  else if ( prob == 2 ) then
    call p02_title ( title )
  else if ( prob == 3 ) then
    call p03_title ( title )
  else if ( prob == 4 ) then
    call p04_title ( title )
  else if ( prob == 5 ) then
    call p05_title ( title )
  else if ( prob == 6 ) then
    call p06_title ( title )
  else if ( prob == 7 ) then
    call p07_title ( title )
  else if ( prob == 8 ) then
    call p08_title ( title )
  else if ( prob == 9 ) then
    call p09_title ( title )
  else if ( prob == 10 ) then
    call p10_title ( title )
  else if ( prob == 11 ) then
    call p11_title ( title )
  else if ( prob == 12 ) then
    call p12_title ( title )
  else if ( prob == 13 ) then
    call p13_title ( title )
  else if ( prob == 14 ) then
    call p14_title ( title )
  else if ( prob == 15 ) then
    call p15_title ( title )
  else if ( prob == 16 ) then
    call p16_title ( title )
  else if ( prob == 17 ) then
    call p17_title ( title )
  else if ( prob == 18 ) then
    call p18_title ( title )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p01_fx ( x, fx )

!*****************************************************************************80
!
!! P01_FX evaluates sin ( x ) - x / 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  fx = sin ( x ) - 0.5D+00 * x

  return
end
subroutine p01_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P01_FX1 evaluates the derivative of the function for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  fx1 = cos ( x ) - 0.5D+00

  return
end
subroutine p01_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P01_FX2 evaluates the second derivative of the function for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = - sin ( x )

  return
end
subroutine p01_range ( range )

!*****************************************************************************80
!
!! P01_RANGE returns an interval bounding the root for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = -1000.0D+00
  range(2) =  1000.0D+00

  return
end
subroutine p01_root ( i, x )

!*****************************************************************************80
!
!! P01_ROOT returns a root for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = - 1.895494267033981D+00
  else if ( i == 2 ) then
    x = 0.0D+00
  else if ( i == 3 ) then
    x = 1.895494267033981D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P01_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p01_root_num ( root_num )

!*****************************************************************************80
!
!! P01_ROOT_NUM returns the number of known roots for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 3

  return
end
subroutine p01_start ( i, x )

!*****************************************************************************80
!
!! P01_START returns a starting point for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.5D+00 * pi
  else if ( i == 2 ) then
    x = pi
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P01_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p01_start_num ( start_num )

!*****************************************************************************80
!
!! P01_START_NUM returns the number of starting points for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 2

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns the title of problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
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

  title = 'F(X) = SIN(X) - 0.5 * X'

  return
end
subroutine p02_fx ( x, fx )

!*****************************************************************************80
!
!! P02_FX evaluates 2 * x - exp ( - x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  fx = 2.0D+00 * x - exp ( - x )

  return
end
subroutine p02_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P02_FX1 evaluates the derivative of the function for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  fx1 = 2.0D+00 + exp ( - x )

  return
end
subroutine p02_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P02_FX2 evaluates the second derivative of the function for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = - exp ( - x )

  return
end
subroutine p02_range ( range )

!*****************************************************************************80
!
!! P02_RANGE returns an interval bounding the root for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = -10.0D+00
  range(2) = 100.0D+00

  return
end
subroutine p02_root ( i, x )

!*****************************************************************************80
!
!! P02_ROOT returns a root for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.35173371124919584D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P02_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p02_root_num ( root_num )

!*****************************************************************************80
!
!! P02_ROOT_NUM returns the number of known roots for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p02_start ( i, x )

!*****************************************************************************80
!
!! P02_START returns a starting point for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.0D+00
  else if ( i == 2 ) then
    x = 1.0D+00
  else if ( i == 3 ) then
    x = -5.0D+00
  else if ( i == 4 ) then
    x = 10.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P02_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p02_start_num ( start_num )

!*****************************************************************************80
!
!! P02_START_NUM returns the number of starting points for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 4

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns the title of problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
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

  title = 'F(X) = 2 * X - EXP ( - X )'

  return
end
subroutine p03_fx ( x, fx )

!*****************************************************************************80
!
!! P03_FX evaluates x * exp ( - x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  fx = x * exp ( - x )

  return
end
subroutine p03_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P03_FX1 evaluates the derivative of the function for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  fx1 = exp ( - x ) * ( 1.0D+00 - x )

  return
end
subroutine p03_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P03_FX2 evaluates the second derivative of the function for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = exp ( - x ) * ( x - 2.0D+00 )

  return
end
subroutine p03_range ( range )

!*****************************************************************************80
!
!! P03_RANGE returns an interval bounding the root for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = -10.0D+00
  range(2) = 100.0D+00

  return
end
subroutine p03_root ( i, x )

!*****************************************************************************80
!
!! P03_ROOT returns a root for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P03_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p03_root_num ( root_num )

!*****************************************************************************80
!
!! P03_ROOT_NUM returns the number of known roots for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p03_start ( i, x )

!*****************************************************************************80
!
!! P03_START returns a starting point for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = -1.0D+00
  else if ( i == 2 ) then
    x =  0.5D+00
  else if ( i == 3 ) then
    x =  2.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P03_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p03_start_num ( start_num )

!*****************************************************************************80
!
!! P03_START_NUM returns the number of starting points for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 3

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns the title of problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
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

  title = 'F(X) = X * EXP ( - X )'

  return
end
subroutine p04_fx ( x, fx )

!*****************************************************************************80
!
!! P04_FX evaluates exp ( x ) - 1 / ( 10 * x )^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  fx = exp ( x ) - 1.0D+00 / ( 100.0D+00 * x * x )

  return
end
subroutine p04_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P04_FX1 evaluates the derivative of the function for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  fx1 = exp ( x ) + 2.0D+00 / ( 100.0D+00 * x**3 )

  return
end
subroutine p04_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P04_FX2 evaluates the second derivative of the function for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = exp ( x ) - 6.0D+00 / ( 100.0D+00 * x**4 )

  return
end
subroutine p04_range ( range )

!*****************************************************************************80
!
!! P04_RANGE returns an interval bounding the root for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) =  0.00001D+00
  range(2) = 20.0D+00

  return
end
subroutine p04_root ( i, x )

!*****************************************************************************80
!
!! P04_ROOT returns a root for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.09534461720025875D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P04_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p04_root_num ( root_num )

!*****************************************************************************80
!
!! P04_ROOT_NUM returns the number of known roots for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p04_start ( i, x )

!*****************************************************************************80
!
!! P04_START returns a starting point for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.03D+00
  else if ( i == 2 ) then
    x = 1.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P04_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p04_start_num ( start_num )

!*****************************************************************************80
!
!! P04_START_NUM returns the number of starting points for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 2

  return
end
subroutine p04_title ( title )

!*****************************************************************************80
!
!! P04_TITLE returns the title of problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
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

  title = 'F(X) = EXP ( X ) - 1 / ( 100 * X * X )'

  return
end
subroutine p05_fx ( x, fx )

!*****************************************************************************80
!
!! P05_FX evaluates ( x + 3 ) * ( x - 1 )^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  fx = ( x + 3.0D+00 ) * ( x - 1.0D+00 ) * ( x - 1.0D+00 )

  return
end
subroutine p05_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P05_FX1 evaluates the derivative of the function for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  fx1 = ( 3.0D+00 * x + 5.0D+00 ) * ( x - 1.0D+00 )

  return
end
subroutine p05_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P05_FX2 evaluates the second derivative of the function for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = 6.0D+00 * x + 2.0D+00

  return
end
subroutine p05_range ( range )

!*****************************************************************************80
!
!! P05_RANGE returns an interval bounding the root for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = -1000.0D+00
  range(2) =  1000.0D+00

  return
end
subroutine p05_root ( i, x )

!*****************************************************************************80
!
!! P05_ROOT returns a root for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = - 3.0D+00
  else if ( i == 2 ) then
    x = 1.0D+00
  else if ( i == 3 ) then
    x = 1.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P05_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p05_root_num ( root_num )

!*****************************************************************************80
!
!! P05_ROOT_NUM returns the number of known roots for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 3

  return
end
subroutine p05_start ( i, x )

!*****************************************************************************80
!
!! P05_START returns a starting point for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x =  2.0D+00
  else if ( i == 2 ) then
    x = - 5.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P05_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p05_start_num ( start_num )

!*****************************************************************************80
!
!! P05_START_NUM returns the number of starting points for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 2

  return
end
subroutine p05_title ( title )

!*****************************************************************************80
!
!! P05_TITLE returns the title of problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
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

  title = 'F(X) = ( X + 3 ) * ( X - 1 ) * ( X - 1 )'

  return
end
subroutine p06_fx ( x, fx )

!*****************************************************************************80
!
!! P06_FX evaluates exp ( x ) - 2 - 1 / ( 10 * x )^2 + 2 / ( 100 * x )^3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  fx = exp ( x ) - 2.0D+00 - 1.0D+00 / ( 10.0D+00 * x )**2 &
    + 2.0D+00 / ( 100.0D+00 * x )**3

  return
end
subroutine p06_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P06_FX1 evaluates the derivative of the function for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  fx1 = exp ( x ) + 2.0D+00 / ( 100.0D+00 * x**3 ) &
    - 6.0D+00 / ( 1000000.0D+00 * x**4 )

  return
end
subroutine p06_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P06_FX2 evaluates the second derivative of the function for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = exp ( x ) - 6.0D+00 / ( 100.0D+00 * x**4 ) &
    + 24.0D+00 / ( 1000000.0D+00 * x**5 )

  return
end
subroutine p06_range ( range )

!*****************************************************************************80
!
!! P06_RANGE returns an interval bounding the root for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) =  0.00001D+00
  range(2) = 20.0D+00

  return
end
subroutine p06_root ( i, x )

!*****************************************************************************80
!
!! P06_ROOT returns a root for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.7032048403631358D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P06_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p06_root_num ( root_num )

!*****************************************************************************80
!
!! P06_ROOT_NUM returns the number of known roots for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p06_start ( i, x )

!*****************************************************************************80
!
!! P06_START returns a starting point for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.0002D+00
  else if ( i == 2 ) then
    x = 2.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P06_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p06_start_num ( start_num )

!*****************************************************************************80
!
!! P06_START_NUM returns the number of starting points for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 2

  return
end
subroutine p06_title ( title )

!*****************************************************************************80
!
!! P06_TITLE returns the title of problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
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

  title = 'F(X) = EXP(X) - 2 - 1 / ( 10 * X )^2 - 2 / ( 100 * X )^3'

  return
end
subroutine p07_fx ( x, fx )

!*****************************************************************************80
!
!! P07_FX evaluates x^3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  fx = x * x * x

  return
end
subroutine p07_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P07_FX1 evaluates the derivative of the function for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  fx1 = 3.0D+00 * x * x

  return
end
subroutine p07_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P07_FX2 evaluates the second derivative of the function for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = 6.0D+00 * x

  return
end
subroutine p07_range ( range )

!*****************************************************************************80
!
!! P07_RANGE returns an interval bounding the root for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = -1000.0D+00
  range(2) =  1000.0D+00

  return
end
subroutine p07_root ( i, x )

!*****************************************************************************80
!
!! P07_ROOT returns a root for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P07_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p07_root_num ( root_num )

!*****************************************************************************80
!
!! P07_ROOT_NUM returns the number of known roots for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p07_start ( i, x )

!*****************************************************************************80
!
!! P07_START returns a starting point for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 1.0D+00
  else if ( i == 2 ) then
    x = -1000.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P07_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p07_start_num ( start_num )

!*****************************************************************************80
!
!! P07_START_NUM returns the number of starting points for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 2

  return
end
subroutine p07_title ( title )

!*****************************************************************************80
!
!! P07_TITLE returns the title of problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
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

  title = 'F(X) = X**3, only linear Newton convergence.'

  return
end
subroutine p08_fx ( x, fx )

!*****************************************************************************80
!
!! P08_FX evaluates cos ( x ) - x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  fx = cos ( x ) - x

  return
end
subroutine p08_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P08_FX1 evaluates the derivative of the function for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  fx1 = - sin ( x ) - 1.0D+00

  return
end
subroutine p08_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P08_FX2 evaluates the second derivative of the function for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = - cos ( x )

  return
end
subroutine p08_range ( range )

!*****************************************************************************80
!
!! P08_RANGE returns an interval bounding the root for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = - 10.0D+00
  range(2) =   10.0D+00

  return
end
subroutine p08_root ( i, x )

!*****************************************************************************80
!
!! P08_ROOT returns a root for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.7390851332151607D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P08_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p08_root_num ( root_num )

!*****************************************************************************80
!
!! P08_ROOT_NUM returns the number of known roots for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p08_start ( i, x )

!*****************************************************************************80
!
!! P08_START returns a starting point for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 1.0D+00
  else if ( i == 2 ) then
    x = 0.5D+00
  else if ( i == 3 ) then
    x = - 1.6D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P08_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p08_start_num ( start_num )

!*****************************************************************************80
!
!! P08_START_NUM returns the number of starting points for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 3

  return
end
subroutine p08_title ( title )

!*****************************************************************************80
!
!! P08_TITLE returns the title of problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
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

  title = 'F(X) = COS(X) - X'

  return
end
subroutine p09_fx ( x, fx )

!*****************************************************************************80
!
!! P09_FX evaluates the Newton Baffler.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  x2 = x - 6.25D+00

  if ( x2 < - 0.25D+00 ) then
    fx = 0.75D+00 * x2 - 0.3125D+00
  else if ( x2 < 0.25D+00 ) then
    fx = 2.0D+00 * x2
  else
    fx = 0.75D+00 * x2 + 0.3125D+00
  end if

  return
end
subroutine p09_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P09_FX1 evaluates the derivative of the function for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  x2 = x - 6.25D+00

  if ( x2 < - 0.25D+00 ) then
    fx1 = 0.75D+00
  else if ( x2 < 0.25D+00 ) then
    fx1 = 2.0D+00
  else
    fx1 = 0.75D+00
  end if

  return
end
subroutine p09_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P09_FX2 evaluates the second derivative of the function for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = 0.0D+00

  return
end
subroutine p09_range ( range )

!*****************************************************************************80
!
!! P09_RANGE returns an interval bounding the root for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = - 4.0D+00
  range(2) = 16.0D+00

  return
end
subroutine p09_root ( i, x )

!*****************************************************************************80
!
!! P09_ROOT returns a root for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 6.25D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P09_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p09_root_num ( root_num )

!*****************************************************************************80
!
!! P09_ROOT_NUM returns the number of known roots for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p09_start ( i, x )

!*****************************************************************************80
!
!! P09_START returns a starting point for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 6.25D+00 + 5.0D+00
  else if ( i == 2 ) then
    x = 6.25D+00 - 1.0D+00
  else if ( i == 3 ) then
    x = 6.25D+00 + 0.1D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P09_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p09_start_num ( start_num )

!*****************************************************************************80
!
!! P09_START_NUM returns the number of starting points for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 3

  return
end
subroutine p09_title ( title )

!*****************************************************************************80
!
!! P09_TITLE returns the title of problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
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

  title = 'The Newton Baffler'

  return
end
subroutine p10_fx ( x, fx )

!*****************************************************************************80
!
!! P10_FX evaluates the Repeller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  fx = 20.0D+00 * x / ( 100.0D+00 * x * x + 1.0D+00 )

  return
end
subroutine p10_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P10_FX1 evaluates the derivative of the function for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  fx1 = ( 1.0D+00 - 10.0D+00 * x ) * ( 1.0D+00 + 10.0D+00 * x ) &
    / ( 100.0D+00 * x * x + 1.0D+00 )**2

  return
end
subroutine p10_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P10_FX2 evaluates the second derivative of the function for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 =  - 200.0D+00 * x * ( 3.0D+00 - 100.0D+00 * x**2 ) &
    / ( 100.0D+00 * x * x + 1.0D+00 )**3

  return
end
subroutine p10_range ( range )

!*****************************************************************************80
!
!! P10_RANGE returns an interval bounding the root for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = - 10.0D+00
  range(2) = + 10.0D+00

  return
end
subroutine p10_root ( i, x )

!*****************************************************************************80
!
!! P10_ROOT returns a root for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P10_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p10_root_num ( root_num )

!*****************************************************************************80
!
!! P10_ROOT_NUM returns the number of known roots for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p10_start ( i, x )

!*****************************************************************************80
!
!! P10_START returns a starting point for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 1.0D+00
  else if ( i == 2 ) then
    x = - 0.14D+00
  else if ( i == 3 ) then
    x = 0.041D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P10_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p10_start_num ( start_num )

!*****************************************************************************80
!
!! P10_START_NUM returns the number of starting points for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 3

  return
end
subroutine p10_title ( title )

!*****************************************************************************80
!
!! P10_TITLE returns the title of problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
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

  title = 'The Repeller'

  return
end
subroutine p11_fx ( x, fx )

!*****************************************************************************80
!
!! P11_FX evaluates the Pinhead.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) epsilon
  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  epsilon = 0.00001D+00

  if ( epsilon == 0.0D+00 ) then

    fx = ( 16.0D+00 - x**4 ) / ( 16.0D+00 * x**4 )

  else

    fx = ( 16.0D+00 - x**4 ) / ( 16.0D+00 * x**4 + epsilon )

  end if

  return
end
subroutine p11_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P11_FX1 evaluates the derivative of the function for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) epsilon
  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  epsilon = 0.00001D+00

  if ( epsilon == 0.0D+00 ) then

    fx1 = - 4.0D+00 / x**5

  else

    fx1 = - 4.0D+00 * x**3 * ( epsilon + 256.0D+00 ) &
      / ( 16.0D+00 * x**4 + epsilon )**2

  end if

  return
end
subroutine p11_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P11_FX2 evaluates the second derivative of the function for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) epsilon
  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  epsilon = 0.00001D+00

  if ( epsilon == 0.0D+00 ) then

    fx2 =  20.0D+00 / x**6

  else

    fx2 = - 4.0D+00 * ( epsilon + 256.0D+00 ) &
      * ( 3.0D+00 * epsilon - 80.0D+00 * x**4 ) * x**2 &
      / ( 16.0D+00 * x**4 + epsilon )**3

  end if

  return
end
subroutine p11_range ( range )

!*****************************************************************************80
!
!! P11_RANGE returns an interval bounding the root for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) =    0.0D+00
  range(2) = + 10.0D+00

  return
end
subroutine p11_root ( i, x )

!*****************************************************************************80
!
!! P11_ROOT returns a root for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = - 2.0D+00
  else if ( i == 2 ) then
    x = 2.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P11_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p11_root_num ( root_num )

!*****************************************************************************80
!
!! P11_ROOT_NUM returns the number of known roots for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 2

  return
end
subroutine p11_start ( i, x )

!*****************************************************************************80
!
!! P11_START returns a starting point for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.25D+00
  else if ( i == 2 ) then
    x = 5.0D+00
  else if ( i == 3 ) then
    x = 1.1D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P11_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p11_start_num ( start_num )

!*****************************************************************************80
!
!! P11_START_NUM returns the number of starting points for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 3

  return
end
subroutine p11_title ( title )

!*****************************************************************************80
!
!! P11_TITLE returns the title of problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 1999
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

  title = 'The Pinhead'

  return
end
subroutine p12_fx ( x, fx )

!*****************************************************************************80
!
!! P12_FX evaluates Flat Stanley.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) factor
  real ( kind = 8 ) fx
  real ( kind = 8 ) s
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  factor = 1000.0D+00

  if ( x == 1.0D+00 ) then

    fx = 0.0D+00

  else

    y = x - 1.0D+00
    s = sign ( 1.0D+00, y )

    fx = s * exp ( log ( factor ) + log ( abs ( y ) ) - 1.0D+00 / y**2 )

  end if

  return
end
subroutine p12_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P12_FX1 evaluates the derivative of the function for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) factor
  real ( kind = 8 ) fx1
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  factor = 1000.0D+00

  if ( x == 1.0D+00 ) then
    fx1 = 0.0D+00
  else
    y = x - 1.0D+00
    fx1 = factor * exp ( - 1.0D+00 / y**2 ) * ( y**2 + 2.0D+00 ) / y**2
  end if

  return
end
subroutine p12_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P12_FX2 evaluates the second derivative of the function for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) factor
  real ( kind = 8 ) fx2
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  factor = 1000.0D+00

  if ( x == 1.0D+00 ) then
    fx2 = 0.0D+00
  else
    y = x - 1.0D+00
    fx2 = - 2.0D+00 * factor * exp ( - 1.0D+00 / y**2 ) &
      * ( y**2 - 2.0D+00 ) / y**5
  end if

  return
end
subroutine p12_range ( range )

!*****************************************************************************80
!
!! P12_RANGE returns an interval bounding the root for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = - 4.0D+00
  range(2) =   4.0D+00

  return
end
subroutine p12_root ( i, x )

!*****************************************************************************80
!
!! P12_ROOT returns a root for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 1.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P12_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p12_root_num ( root_num )

!*****************************************************************************80
!
!! P12_ROOT_NUM returns the number of known roots for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p12_start ( i, x )

!*****************************************************************************80
!
!! P12_START returns a starting point for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 2.0D+00
  else if ( i == 2 ) then
    x = 0.50D+00
  else if ( i == 3 ) then
    x = 4.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P12_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p12_start_num ( start_num )

!*****************************************************************************80
!
!! P12_START_NUM returns the number of starting points for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 3

  return
end
subroutine p12_title ( title )

!*****************************************************************************80
!
!! P12_TITLE returns the title of problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1999
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

  title = 'Flat Stanley (ALL derivatives are zero at the root.)'

  return
end
subroutine p13_fx ( x, fx )

!*****************************************************************************80
!
!! P13_FX evaluates Lazy Boy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) slope
  real ( kind = 8 ) x

  slope = 0.00000000001D+00

  fx = slope * ( x - 100.0D+00 )

  return
end
subroutine p13_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P13_FX1 evaluates the derivative of the function for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) slope
  real ( kind = 8 ) x

  slope = 0.00000000001D+00
  fx1 = slope

  return
end
subroutine p13_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P13_FX2 evaluates the second derivative of the function for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = 0.0D+00

  return
end
subroutine p13_range ( range )

!*****************************************************************************80
!
!! P13_RANGE returns an interval bounding the root for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = - 10000000000000.0D+00
  range(2) =   10000000000000.0D+00

  return
end
subroutine p13_root ( i, x )

!*****************************************************************************80
!
!! P13_ROOT returns a root for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 100.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P13_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p13_root_num ( root_num )

!*****************************************************************************80
!
!! P13_ROOT_NUM returns the number of known roots for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p13_start ( i, x )

!*****************************************************************************80
!
!! P13_START returns a starting point for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 100000000.0D+00
  else if ( i == 2 ) then
    x = 100000013.0D+00
  else if ( i == 3 ) then
    x = - 100000000000.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P13_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p13_start_num ( start_num )

!*****************************************************************************80
!
!! P13_START_NUM returns the number of starting points for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 3

  return
end
subroutine p13_title ( title )

!*****************************************************************************80
!
!! P13_TITLE returns the title of problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1999
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

  title = 'Lazy Boy (Linear function, almost flat.)'

  return
end
subroutine p14_fx ( x, fx )

!*****************************************************************************80
!
!! P14_FX evaluates the Camel.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  fx =   1.0D+00 / ( ( x - 0.3D+00 )**2 + 0.01D+00 ) &
       + 1.0D+00 / ( ( x - 0.9D+00 )**2 + 0.04D+00 ) + 2.0D+00 * x - 5.2D+00

  return
end
subroutine p14_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P14_FX1 evaluates the derivative of the function for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  fx1 = - 2.0D+00 * ( x - 0.3D+00 ) / ( ( x - 0.3D+00 )**2 + 0.01D+00 )**2 &
        - 2.0D+00 * ( x - 0.9D+00 ) / ( ( x - 0.9D+00 )**2 + 0.04D+00 )**2 &
        + 2.0D+00

  return
end
subroutine p14_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P14_FX2 evaluates the second derivative of the function for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = 0.0D+00

  return
end
subroutine p14_range ( range )

!*****************************************************************************80
!
!! P14_RANGE returns an interval bounding the root for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = - 10.0D+00
  range(2) =   10.0D+00

  return
end
subroutine p14_root ( i, x )

!*****************************************************************************80
!
!! P14_ROOT returns a root for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = - 0.1534804948126991D+00
  else if ( i == 2 ) then
    x = 1.8190323925159182D+00
  else if ( i == 3 ) then
    x = 2.1274329318603367D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P14_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p14_root_num ( root_num )

!*****************************************************************************80
!
!! P14_ROOT_NUM returns the number of known roots for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 3

  return
end
subroutine p14_start ( i, x )

!*****************************************************************************80
!
!! P14_START returns a starting point for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 3.0D+00
  else if ( i == 2 ) then
    x = - 0.5D+00
  else if ( i == 3 ) then
    x = 0.0D+00
  else if ( i == 4 ) then
    x = 2.12742D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P14_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p14_start_num ( start_num )

!*****************************************************************************80
!
!! P14_START_NUM returns the number of starting points for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 3

  return
end
subroutine p14_title ( title )

!*****************************************************************************80
!
!! P14_TITLE returns the title of problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 1999
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

  title = 'The Camel (double hump and some shallow roots.)'

  return
end
subroutine p15_fx ( x, fx )

!*****************************************************************************80
!
!! P15_FX evaluates a pathological function for Newton's method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    George Donovan, Arnold Miller, Timothy Moreland,
!    Pathological Functions for Newton's Method,
!    American Mathematical Monthly, January 1993, pages 53-58.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) r8_cube_root
  real ( kind = 8 ) x

  fx = r8_cube_root ( x ) * exp ( - x**2 )

  return
end
subroutine p15_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P15_FX1 evaluates the derivative of the function for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) r8_cube_root
  real ( kind = 8 ) x

  fx1 = ( 1.0D+00 - 6.0D+00 * x**2 ) * r8_cube_root ( x ) &
    * exp ( - x**2 ) / ( 3.0D+00 * x )

  return
end
subroutine p15_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P15_FX2 evaluates the second derivative of the function for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) r8_cube_root
  real ( kind = 8 ) x

  fx2 = ( - 2.0D+00 - 30.0D+00 * x**2 + 36.0D+00 * x**4 ) * r8_cube_root ( x ) &
    * exp ( - x**2 ) / ( 9.0D+00 * x**2 )

  return
end
subroutine p15_range ( range )

!*****************************************************************************80
!
!! P15_RANGE returns an interval bounding the root for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = - 10.0D+00
  range(2) =   10.0D+00

  return
end
subroutine p15_root ( i, x )

!*****************************************************************************80
!
!! P15_ROOT returns a root for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P15_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p15_root_num ( root_num )

!*****************************************************************************80
!
!! P15_ROOT_NUM returns the number of known roots for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p15_start ( i, x )

!*****************************************************************************80
!
!! P15_START returns a starting point for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.01D+00
  else if ( i == 2 ) then
    x = - 0.25D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P15_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p15_start_num ( start_num )

!*****************************************************************************80
!
!! P15_START_NUM returns the number of starting points for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 2

  return
end
subroutine p15_title ( title )

!*****************************************************************************80
!
!! P15_TITLE returns the title of problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2000
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

  title = 'Donovan/Miller/Moreland Pathological Function'

  return
end
subroutine p16_fx ( x, fx )

!*****************************************************************************80
!
!! P16_FX evaluates Kepler's Equation.
!
!  Discussion:
!
!    This is Kepler's equation.  The equation has the form:
!
!      X = M + E * sin ( X )
!
!    X represents the eccentric anomaly of a planet, the angle between the
!    perihelion (the point on the orbit nearest to the sun)
!    through the sun to the center of the ellipse, and the
!    line from the center of the ellipse to the planet.
!
!    There are two parameters:
!
!    E is the eccentricity of the orbit, which should be between 0 and 1.0;
!
!    M is the angle from the perihelion made by a fictitious planet traveling
!    on a circular orbit centered at the sun, and traveling at a constant
!    angular velocity equal to the average angular velocity of the true planet.
!    M is usually between 0 and 180 degrees, but can have any value.
!
!    For convenience, X and M are measured in degrees.
!
!    Sample results:
!
!    E        M      X
!    -----  ---  ----------
!    0.100    5    5.554589
!    0.200    5    6.246908
!    0.300    5    7.134960
!    0.400    5    8.313903
!    0.500    5    9.950063
!    0.600    5   12.356653
!    0.700    5   16.167990
!    0.800    5   22.656579
!    0.900    5   33.344447
!    0.990    5   45.361023
!    0.990    1   24.725822
!    0.990   33   89.722155
!    0.750   70  110.302
!    0.990    2   32.361007
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Colwell,
!    Solving Kepler's Equation Over Three Centuries,
!    Willmann-Bell, 1993
!
!    Jean Meeus,
!    Astronomical Algorithms,
!    Willman-Bell, Inc, 1991,
!    QB51.3.E43M42
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) e
  real ( kind = 8 ) fx
  real ( kind = 8 ) m
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  e = 0.8D+00
  m = 5.0D+00

  fx = ( pi * ( x - m ) / 180.0D+00 ) - e * sin ( pi * x / 180.0D+00 )

  return
end
subroutine p16_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P16_FX1 evaluates the derivative of the function for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) e
  real ( kind = 8 ) fx1
  real ( kind = 8 ) m
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  e = 0.8D+00
  m = 5.0D+00

  fx1 = ( pi / 180.0D+00 ) &
    - e * pi * cos ( pi * x / 180.0D+00  ) / 180.0D+00

  return
end
subroutine p16_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P16_FX2 evaluates the second derivative of the function for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) e
  real ( kind = 8 ) fx2
  real ( kind = 8 ) m
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  e = 0.8D+00
  m = 5.0D+00

  fx2 = e * pi * pi * sin ( pi * x / 180.0D+00  ) / 180.0D+00**2

  return
end
subroutine p16_range ( range )

!*****************************************************************************80
!
!! P16_RANGE returns an interval bounding the root for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) e
  real ( kind = 8 ) m
  real ( kind = 8 ) range(2)

  e = 0.8D+00
  m = 5.0D+00

  range(1) = m - 180.0D+00
  range(2) = m + 180.0D+00

  return
end
subroutine p16_root ( i, x )

!*****************************************************************************80
!
!! P16_ROOT returns a root for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P16_ROOT - Fatal error!'
  write ( *, '(a,i4)' ) '  Illegal root index = ', i

  stop
end
subroutine p16_root_num ( root_num )

!*****************************************************************************80
!
!! P16_ROOT_NUM returns the number of known roots for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 0

  return
end
subroutine p16_start ( i, x )

!*****************************************************************************80
!
!! P16_START returns a starting point for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  real ( kind = 8 ) e
  integer ( kind = 4 ) i
  real ( kind = 8 ) m
  real ( kind = 8 ) x

  e = 0.8D+00
  m = 5.0D+00

  if ( i == 1 ) then
    x = 0.0D+00
  else if ( i == 2 ) then
    x = m
  else if ( i == 3 ) then
    x = m + 180.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P16_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p16_start_num ( start_num )

!*****************************************************************************80
!
!! P16_START_NUM returns the number of starting points for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 3

  return
end
subroutine p16_title ( title )

!*****************************************************************************80
!
!! P16_TITLE returns the title of problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2001
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

  title = 'Kepler''s Eccentric Anomaly Equation, in degrees'

  return
end
subroutine p17_fx ( x, fx )

!*****************************************************************************80
!
!! P17_FX evaluates the function for problem 17.
!
!  Discussion:
!
!    This simple example is of historical interest, since it was used
!    by Wallis to illustrate the use of Newton's method, and has been
!    a common example ever since.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  fx = x**3 - 2.0D+00 * x - 5.0D+00

  return
end
subroutine p17_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P17_FX1 evaluates the derivative of the function for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  fx1 = 3.0D+00 * x * x - 2.0D+00

  return
end
subroutine p17_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P17_FX2 evaluates the second derivative of the function for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = 6.0D+00 * x

  return
end
subroutine p17_range ( range )

!*****************************************************************************80
!
!! P17_RANGE returns an interval bounding the root for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = 2.0D+00
  range(2) = 3.0D+00

  return
end
subroutine p17_root ( i, x )

!*****************************************************************************80
!
!! P17_ROOT returns a root for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 2.0945514815423265D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P17_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p17_root_num ( root_num )

!*****************************************************************************80
!
!! P17_ROOT_NUM returns the number of known roots for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p17_start ( i, x )

!*****************************************************************************80
!
!! P17_START returns a starting point for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 2.0D+00
  else if ( i == 2 ) then
    x = 3.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P17_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p17_start_num ( start_num )

!*****************************************************************************80
!
!! P17_START_NUM returns the number of starting points for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 2

  return
end
subroutine p17_title ( title )

!*****************************************************************************80
!
!! P17_TITLE returns the title of problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2001
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

  title = 'The Wallis example, x^3-2x-5=0'

  return
end
subroutine p18_fx ( x, fx )

!*****************************************************************************80
!
!! P18_FX evaluates the function for problem 18.
!
!  Discussion:
!
!    F(X) = 10^14 * (x-1)^7, but is written in term by term form.
!
!    This polynomial becomes difficult to evaluate accurately when 
!    written this way.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) x

  fx = 10.0D+00**14 * ( &
                   x**7 &
      -  7.0D+00 * x**6 &
      + 21.0D+00 * x**5 &
      - 35.0D+00 * x**4 &
      + 35.0D+00 * x**3 &
      - 21.0D+00 * x**2 &
      +  7.0D+00 * x    &
      -  1.0D+00 )

  return
end
subroutine p18_fx1 ( x, fx1 )

!*****************************************************************************80
!
!! P18_FX1 evaluates the derivative of the function for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX1, the first derivative of the function at X.
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) x

  fx1 = 10.0D+00**14 * ( &
          7.0D+00 * x**6 &
      -  42.0D+00 * x**5 &
      + 105.0D+00 * x**4 &
      - 140.0D+00 * x**3 &
      + 105.0D+00 * x**2 &
      -  42.0D+00 * x    &
      +   7.0D+00 )

  return
end
subroutine p18_fx2 ( x, fx2 )

!*****************************************************************************80
!
!! P18_FX2 evaluates the second derivative of the function for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the abscissa.
!
!    Output, real ( kind = 8 ) FX2, the second derivative at X.
!
  implicit none

  real ( kind = 8 ) fx2
  real ( kind = 8 ) x

  fx2 = 10.0D+00**14 * ( &
         42.0D+00 * x**5 &
      - 210.0D+00 * x**4 &
      + 420.0D+00 * x**3 &
      - 420.0D+00 * x**2 &
      + 210.0D+00 * x    &
      -  42.0D+00 )

  return
end
subroutine p18_range ( range )

!*****************************************************************************80
!
!! P18_RANGE returns an interval bounding the root for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RANGE(2), the minimum and maximum values of
!    an interval containing the root.
!
  implicit none

  real ( kind = 8 ) range(2)

  range(1) = 0.988D+00
  range(2) = 1.012D+00

  return
end
subroutine p18_root ( i, x )

!*****************************************************************************80
!
!! P18_ROOT returns a root for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the requested root.
!
!    Output, real ( kind = 8 ) X, the value of the root.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 1.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P18_ROOT - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal root index = ', i
    stop
  end if

  return
end
subroutine p18_root_num ( root_num )

!*****************************************************************************80
!
!! P18_ROOT_NUM returns the number of known roots for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ROOT_NUM, the number of known roots.
!
  implicit none

  integer ( kind = 4 ) root_num

  root_num = 1

  return
end
subroutine p18_start ( i, x )

!*****************************************************************************80
!
!! P18_START returns a starting point for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the starting point.
!
!    Output, real ( kind = 8 ) X, the starting point.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 0.990D+00
  else if ( i == 2 ) then
    x = 1.013D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P18_START - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal starting point index = ', i
    stop
  end if

  return
end
subroutine p18_start_num ( start_num )

!*****************************************************************************80
!
!! P18_START_NUM returns the number of starting points for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) START_NUM, the number of starting points.
!
  implicit none

  integer ( kind = 4 ) start_num

  start_num = 2

  return
end
subroutine p18_title ( title )

!*****************************************************************************80
!
!! P18_TITLE returns the title of problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2011
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

  title = '10^14 * (x-1)^7, written term by term.'

  return
end
function r8_cube_root ( x )

!*****************************************************************************80
!
!! R8_CUBE_ROOT returns the cube root of a real number.
!
!  Discussion:
!
!    This routine is designed to avoid the possible problems that can occur
!    when formulas like 0.0**(1/3) or (-1.0)**(1/3) are to be evaluated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose cube root is desired.
!
!    Output, real ( kind = 8 ) R8_CUBE_ROOT, the cube root of X.
!
  implicit none

  real ( kind = 8 ) r8_cube_root
  real ( kind = 8 ) x

  if ( 0.0D+00 < x ) then
    r8_cube_root = x**(1.0D+00/3.0D+00)
  else if ( x == 0.0D+00 ) then
    r8_cube_root = 0.0D+00
  else
    r8_cube_root = - ( abs ( x ) )**(1.0D+00/3.0D+00)
  end if

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real    ( kind = 8 ) r8_huge

  r8_huge = 1.0D+30

  return
end
function r8_sign ( x )

!*****************************************************************************80
!
!! R8_SIGN returns the sign of an R8.
!
!  Discussion:
!
!    value = -1 if X < 0;
!    value =  0 if X => 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose sign is desired.
!
!    Output, real ( kind = 8 ) R8_SIGN, the sign of X:
!
  implicit none

  real    ( kind = 8 ) r8_sign
  real    ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    r8_sign = -1.0D+00
  else
    r8_sign = +1.0D+00
  end if

  return
end
subroutine r8poly2_rroot ( a, b, c, r1, r2 )

!*****************************************************************************80
!
!! R8POLY2_RROOT returns the real parts of the roots of a quadratic polynomial.
!
!  Example:
!
!    A    B    C       roots              R1   R2
!   --   --   --     ------------------   --   --
!    1   -4    3     1          3          1    3
!    1    0    4         2*i      - 2*i    0    0
!    2   -6    5     3 +   i    3 -   i    3    3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the coefficients of the quadratic
!    polynomial A * X**2 + B * X + C = 0 whose roots are desired.
!    A must not be zero.
!
!    Output, real ( kind = 8 ) R1, R2, the real parts of the two roots
!    of the polynomial.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) disc
  real ( kind = 8 ) q
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_RROOT - Fatal error!'
    write ( *, '(a)' ) '  The coefficient A is zero.'
    stop
  end if

  disc = b * b - 4.0D+00 * a * c
  disc = max ( disc, 0.0D+00 )

  q = ( b + sign ( 1.0D+00, b ) * sqrt ( disc ) )
  r1 = - 0.5D+00 * q / a
  r2 = - 2.0D+00 * c / q

  return
end
subroutine regula_falsi ( fatol, step_max, prob, xatol, xa, xb, fxa, fxb )

!*****************************************************************************80
!
!! REGULA_FALSI carries out the Regula Falsi method to seek a root of F(X) = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) FATOL, an absolute error tolerance for the
!    function value of the root.  If an approximate root X satisfies
!      ABS ( F ( X ) ) <= FATOL, then X will be accepted as the
!    root and the iteration will be terminated.
!
!    Input, integer ( kind = 4 ) STEP_MAX, the maximum number of steps allowed
!    for an iteration.
!
!    Input, integer ( kind = 4 ) PROB, the index of the function whose root is
!    to be sought.
!
!    Input, real ( kind = 8 ) XATOL, an absolute error tolerance for the root.
!
!    Input/output, real ( kind = 8 ) XA, XB, two points at which the
!    function differs in sign.  On output, these values have been adjusted
!    to a smaller interval.
!
!    Input/output, real ( kind = 8 ) FXA, FXB, the value of the function
!    at XA and XB.
!
  implicit none

  real ( kind = 8 ) fatol
  real ( kind = 8 ) fxa
  real ( kind = 8 ) fxb
  real ( kind = 8 ) fxc
  integer ( kind = 4 ) step_max
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) step_num
  real ( kind = 8 ) t
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
  real ( kind = 8 ) xatol
  real ( kind = 8 ) xc
!
!  The method requires a change-of-sign interval.
!
  if ( sign ( 1.0D+00, fxa ) == sign ( 1.0D+00, fxb ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REGULA_FALSI - Fatal error!'
    write ( *, '(a)' ) '  Function does not change sign at endpoints.'
    stop
  end if
!
!  Make A the root with negative F, B the root with positive F.
!
  if ( 0.0D+00 < fxa ) then
    t = xa
    xa = xb
    xb = t
    t = fxa
    fxa = fxb
    fxb = t
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REGULA FALSI'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Step      XA            XB             F(XA)         F(XB)'
  write ( *, '(a)' ) ' '

  step_num = 0
  write ( *, '(2x,i4,2x,2g16.8,2g14.6)' ) step_num, xa, xb, fxa, fxb

  do step_num = 1, step_max

    if ( abs ( xa - xb ) < xatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Interval small enough for convergence.'
      return
    end if

    if ( abs ( fxa ) <= fatol .or. abs ( fxb ) <= fatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Function small enough for convergence.'
      return
    end if

    xc = ( fxa * xb - fxb * xa ) / ( fxa - fxb )
    call p00_fx ( prob, xc, fxc )

    if ( fxc < 0.0D+00 ) then
      xa = xc
      fxa = fxc
    else
      xb = xc
      fxb = fxc
    end if

    write ( *, '(2x,i4,2x,2g16.8,2g14.6)' ) step_num, xa, xb, fxa, fxb

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Took maximum number of steps without convergence.'

  return
end
subroutine secant ( fatol, step_max, prob, xatol, xmin, xmax, xa, xb, fxa, &
  fxb )

!*****************************************************************************80
!
!! SECANT carries out the secant method to seek a root of F(X) = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) FATOL, an absolute error tolerance for the
!    function value of the root.  If an approximate root X satisfies
!      ABS ( F ( X ) ) <= FATOL, then X will be accepted as the
!    root and the iteration will be terminated.
!
!    Input, integer ( kind = 4 ) STEP_MAX, the maximum number of steps allowed
!    for an iteration.
!
!    Input, integer ( kind = 4 ) PROB, the index of the function whose root is
!    to be sought.
!
!    Input, real ( kind = 8 ) XATOL, an absolute error tolerance for the root.
!
!    Input, real ( kind = 8 ) XMAX, XMIN, the interval in which the root should
!    be sought.
!
!    Input/output, real ( kind = 8 ) XA, XB, two points at which the
!    function differs in sign.  On output, these values have been adjusted
!    to a smaller interval.
!
!    Input/output, real ( kind = 8 ) FXA, FXB, the value of the function
!    at XA and XB.
!
  implicit none

  real ( kind = 8 ) fatol
  real ( kind = 8 ) fxa
  real ( kind = 8 ) fxb
  real ( kind = 8 ) fxc
  integer ( kind = 4 ) step_max
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) step_num
  real ( kind = 8 ) xa
  real ( kind = 8 ) xatol
  real ( kind = 8 ) xb
  real ( kind = 8 ) xc
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SECANT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X             F(X)'
  write ( *, '(a)' ) ' '

  step_num = -1
  write ( *, '(2x,i4,2x,g16.8,g14.6)' ) step_num, xa, fxa

  if ( abs ( fxa ) <= fatol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Function small enough for convergence.'
    return
  end if

  step_num = 0
  write ( *, '(2x,i4,2x,g16.8,g14.6)' ) step_num, xb, fxb

  do step_num = 1, step_max

    if ( abs ( fxb ) <= fatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Function small enough for convergence.'
      return
    end if

    if ( abs ( xa - xb ) < xatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Interval small enough for convergence.'
      return
    end if

    if ( xb < xmin .or. xmax < xb ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Iterate has left the region [XMIN,XMAX].'
      return
    end if

    if ( fxa == fxb ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(A) = F(B), algorithm fails.'
      return
    end if

    xc = ( fxa * xb - fxb * xa ) / ( fxa - fxb )

    call p00_fx ( prob, xc, fxc )

    xa = xb
    fxa = fxb
    xb = xc
    fxb = fxc
    write ( *, '(2x,i4,2x,g16.8,g14.6)' ) step_num, xb, fxb

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Took maximum number of steps.'

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
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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
