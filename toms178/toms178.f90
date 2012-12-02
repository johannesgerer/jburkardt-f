function best_nearby ( delta, point, prevbest, nvars, f, funevals )

!*****************************************************************************80
!
!! BEST_NEARBY looks for a better nearby point, one coordinate at a time.
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    Original ALGOL version by Arthur Kaupe.
!    C version by Mark Johnson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    M Bell, Malcolm Pike,
!    Remark on Algorithm 178: Direct Search,
!    Communications of the ACM,
!    Volume 9, Number 9, September 1966, page 684.
!
!    Robert Hooke, Terry Jeeves,
!    Direct Search Solution of Numerical and Statistical Problems,
!    Journal of the ACM,
!    Volume 8, Number 2, April 1961, pages 212-229.
!
!    Arthur Kaupe,
!    Algorithm 178:
!    Direct Search,
!    Communications of the ACM,
!    Volume 6, Number 6, June 1963, page 313.
!
!    FK Tomlin, LB Smith,
!    Remark on Algorithm 178: Direct Search,
!    Communications of the ACM,
!    Volume 12, Number 11, November 1969, page 637-638.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DELTA(NVARS), the size of a step in
!    each direction.
!
!    Input/output, real ( kind = 8 ) POINT(NVARS); on input, the current
!    candidate.  On output, the value of POINT may have been updated.
!
!    Input, real ( kind = 8 ) PREVBEST, the minimum value of the function seen
!    so far.
!
!    Input, integer ( kind = 4 ) NVARS, the number of variables.
!
!    Input, external real ( kind = 8 ) F, the name of the function routine,
!    which should have the form:
!      function f ( x, n )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!
!    Input/output, integer ( kind = 4 ) FUNEVALS, the number of function
!    evaluations.
!
!    Output, real ( kind = 8 ) BEST_NEARBY, the minimum value of the function
!    seen after checking the nearby neighbors.
!
  implicit none

  integer ( kind = 4 ) nvars

  real ( kind = 8 ) best_nearby
  real ( kind = 8 ) delta(nvars)
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) ftmp
  integer ( kind = 4 ) funevals
  integer ( kind = 4 ) i
  real ( kind = 8 ) minf
  real ( kind = 8 ) point(nvars)
  real ( kind = 8 ) prevbest
  real ( kind = 8 ) z(nvars)

  minf = prevbest
  z(1:nvars) = point(1:nvars)

  do i = 1, nvars

    z(i) = point(i) + delta(i)

    ftmp = f ( z, nvars )
    funevals = funevals + 1

    if ( ftmp < minf ) then

      minf = ftmp

    else

      delta(i) = - delta(i)
      z(i) = point(i) + delta(i)
      ftmp = f ( z, nvars )
      funevals = funevals + 1

      if ( ftmp < minf ) then
        minf = ftmp
      else
        z(i) = point(i)
      end if

    end if

  end do

  point(1:nvars) = z(1:nvars)
  best_nearby = minf

  return
end
function hooke ( nvars, startpt, endpt, rho, eps, itermax, f )

!*****************************************************************************80
!
!! HOOKE seeks a minimizer of a scalar function of several variables.
!
!  Discussion:
!
!    This routine find a point X where the nonlinear objective function
!    F(X) has a local minimum.  X is an N-vector and F(X) is a scalar.
!    The objective function F(X) is not required to be differentiable
!    or even continuous.  The program does not use or require derivatives
!    of the objective function.
!
!    The user supplies three things:
!    1) a subroutine that computes F(X),
!    2) an initial "starting guess" of the minimum point X,
!    3) values for the algorithm convergence parameters.
!
!    The program searches for a local minimum, beginning from the
!    starting guess, using the Direct Search algorithm of Hooke and
!    Jeeves.
!
!    This program is adapted from the Algol pseudocode found in the
!    paper by Kaupe, and includes improvements suggested by Bell and Pike,
!    and by Tomlin and Smith.
!
!    The algorithm works by taking "steps" from one estimate of
!    a minimum, to another (hopefully better) estimate.  Taking
!    big steps gets to the minimum more quickly, at the risk of
!    "stepping right over" an excellent point.  The stepsize is
!    controlled by a user supplied parameter called RHO.  At each
!    iteration, the stepsize is multiplied by RHO  (0 < RHO < 1),
!    so the stepsize is successively reduced.
!
!    Small values of rho correspond to big stepsize changes,
!    which make the algorithm run more quickly.  However, there
!    is a chance (especially with highly nonlinear functions)
!    that these big changes will accidentally overlook a
!    promising search vector, leading to nonconvergence.
!
!    Large values of RHO correspond to small stepsize changes,
!    which force the algorithm to carefully examine nearby points
!    instead of optimistically forging ahead.  This improves the
!    probability of convergence.
!
!    The stepsize is reduced until it is equal to (or smaller
!    than) EPS.  So the number of iterations performed by
!    Hooke-Jeeves is determined by RHO and EPS:
!
!      RHO^(number_of_iterations) = EPS
!
!    In general it is a good idea to set RHO to an aggressively
!    small value like 0.5 (hoping for fast convergence).  Then,
!    if the user suspects that the reported minimum is incorrect
!    (or perhaps not accurate enough), the program can be run
!    again with a larger value of RHO such as 0.85, using the
!    result of the first minimization as the starting guess to
!    begin the second minimization.
!
!    Normal use:
!    (1) Code your function F() in the C language;
!    (2) Install your starting guess;
!    (3) Run the program.
!
!    If there are doubts about the result, the computed minimizer
!    can be used as the starting point for a second minimization attempt.
!
!    To apply this method to data fitting, code your function F() to be
!    the sum of the squares of the errors (differences) between the
!    computed values and the measured values.  Then minimize F()
!    using Hooke-Jeeves.
!
!    For example, you have 20 datapoints (T(i), Y(i)) and you want to
!    find A, B and C so that:
!
!      A*t*t + B*exp(t) + C*tan(t)
!
!    fits the data as closely as possible.  Then the objective function
!    F() to be minimized is just
!
!      F(A,B,C) = sum ( 1 <= i <= 20 )
!        ( y(i) - A*t(i)*t(i) - B*exp(t(i)) - C*tan(t(i)) )^2.
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    ALGOL original by Arthur Kaupe.
!    C version by Mark Johnson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    M Bell, Malcolm Pike,
!    Remark on Algorithm 178: Direct Search,
!    Communications of the ACM,
!    Volume 9, Number 9, September 1966, page 684.
!
!    Robert Hooke, Terry Jeeves,
!    Direct Search Solution of Numerical and Statistical Problems,
!    Journal of the ACM,
!    Volume 8, Number 2, April 1961, pages 212-229.
!
!    Arthur Kaupe,
!    Algorithm 178:
!    Direct Search,
!    Communications of the ACM,
!    Volume 6, Number 6, June 1963, page 313.
!
!    FK Tomlin, LB Smith,
!    Remark on Algorithm 178: Direct Search,
!    Communications of the ACM,
!    Volume 12, Number 11, November 1969, page 637-638.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVARS, the number of spatial dimensions.
!
!    Input, real ( kind = 8 ) STARTPT(NVARS), the user-supplied
!    initial estimate for the minimizer.
!
!    Output, real ( kind = 8 ) ENDPT(NVARS), the estimate for the
!    minimizer, as calculated by the program.
!
!    Input, real ( kind = 8 ) RHO, a user-supplied convergence parameter
!    which should be set to a value between 0.0 and 1.0.  Larger values
!    of RHO give greater probability of convergence on highly nonlinear
!    functions, at a cost of more function evaluations.  Smaller
!    values of RHO reduce the number of evaluations and the program
!    running time, but increases the risk of nonconvergence.
!
!    Input, real ( kind = 8 ) EPS, the criterion for halting
!    the search for a minimum.  When the algorithm
!    begins to make less and less progress on each
!    iteration, it checks the halting criterion: if
!    the stepsize is below EPS, terminate the
!    iteration and return the current best estimate
!    of the minimum.  Larger values of EPS (such
!    as 1.0e-4) give quicker running time, but a
!    less accurate estimate of the minimum.  Smaller
!    values of EPS (such as 1.0e-7) give longer
!    running time, but a more accurate estimate of
!    the minimum.
!
!    Input, integer ( kind = 4 ) ITERMAX, a limit on the number of iterations.
!
!    Input, external real ( kind = 8 ) F, the name of the function routine,
!    which should have the form:
!      function f ( x, n )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!
!    Output, integer ( kind = 4 ) HOOKE, the number of iterations taken.
!
  implicit none

  integer ( kind = 4 ) nvars

  real ( kind = 8 ) best_nearby
  real ( kind = 8 ) delta(nvars)
  real ( kind = 8 ) endpt(nvars)
  real ( kind = 8 ) eps
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fbefore
  integer ( kind = 4 ) funevals
  integer ( kind = 4 ) hooke
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itermax
  integer ( kind = 4 ) iters
  integer ( kind = 4 ) j
  integer ( kind = 4 ) keep
  real ( kind = 8 ) newf
  real ( kind = 8 ) newx(nvars)
  real ( kind = 8 ) rho
  real ( kind = 8 ) startpt(nvars)
  real ( kind = 8 ) steplength
  real ( kind = 8 ) tmp
  logical, parameter :: verbose = .false.
  real ( kind = 8 ) xbefore(nvars)

  newx(1:nvars) = startpt(1:nvars)
  xbefore(1:nvars) = startpt(1:nvars)

  do i = 1, nvars
    if ( startpt(i) == 0.0D+00 ) then
      delta(i) = rho
    else
      delta(i) = rho * abs ( startpt(i) )
    end if
  end do

  funevals = 0
  steplength = rho
  iters = 0
  fbefore = f ( newx, nvars )
  funevals = funevals + 1
  newf = fbefore

  do while ( iters < itermax .and. eps < steplength )

    iters = iters + 1

    if ( verbose ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,g14.6)' ) &
      '  FUNEVALS, = ', funevals, '  F(X) = ', fbefore

      do j = 1, nvars
        write ( *, '(2x,i8,2x,g14.6)' ) j, xbefore(j)
      end do

    end if
!
!  Find best new point, one coordinate at a time.
!
    newx(1:nvars) = xbefore(1:nvars)

    newf = best_nearby ( delta, newx, fbefore, nvars, f, funevals )
!
!  If we made some improvements, pursue that direction.
!
    keep = 1

    do while ( newf < fbefore .and. keep == 1 )

      do i = 1, nvars
!
!  Arrange the sign of DELTA.
!
        if ( newx(i) <= xbefore(i) ) then
          delta(i) = - abs ( delta(i) )
        else
          delta(i) = abs ( delta(i) )
        end if
!
!  Now, move further in this direction.
!
        tmp = xbefore(i)
        xbefore(i) = newx(i)
        newx(i) = newx(i) + newx(i) - tmp
      end do

      fbefore = newf
      newf = best_nearby ( delta, newx, fbefore, nvars, f, funevals )
!
!  If the further (optimistic) move was bad...
!
      if ( fbefore <= newf ) then
        exit
      end if
!
!  Make sure that the differences between the new and the old points
!  are due to actual displacements; beware of roundoff errors that
!  might cause NEWF < FBEFORE.
!
      keep = 0

      do i = 1, nvars
        if ( 0.5D+00 * abs ( delta(i) ) < &
          abs ( newx(i) - xbefore(i) ) ) then
          keep = 1
          exit
        end if
      end do

    end do

    if ( eps <= steplength .and. fbefore <= newf ) then
      steplength = steplength * rho
      delta(1:nvars) = delta(1:nvars) * rho
    end if

  end do

  endpt(1:nvars) = xbefore(1:nvars)

  hooke = iters

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
