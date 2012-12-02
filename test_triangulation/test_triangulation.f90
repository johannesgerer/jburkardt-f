function angle_contains_point_2d ( p1, p2, p3, p )

!*****************************************************************************80
!
!! ANGLE_CONTAINS_POINT_2D determines if an angle contains a point, in 2D.
!
!  Discussion:
!
!    The angle is defined by the sequence of points P1, P2 and P3.
!
!    The point is "contained" by the angle if the ray P - P2
!    is between (in a counterclockwise sense) the rays P1 - P2
!    and P3 - P2.
!
!        P1
!        /
!       /   P
!      /  .
!     / .
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), P3(2), three points that define
!    the angle.  The order of these points matters!
!
!    Input, real ( kind = 8 ) P(2), the point to be checked.
!
!    Output, logical ANGLE_CONTAINS_POINT_2D, is TRUE if the point
!    is inside the angle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  logical angle_contains_point_2d
  real ( kind = 8 ) angle_rad_2d
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)

  if ( angle_rad_2d ( p1, p2, p ) <= angle_rad_2d ( p1, p2, p3 ) ) then
    angle_contains_point_2d = .true.
  else
    angle_contains_point_2d = .false.
  end if

  return
end
subroutine angle_half_2d ( p1, p2, p3, p4 )

!*****************************************************************************80
!
!! ANGLE_HALF_2D finds half an angle in 2D.
!
!  Discussion:
!
!    The original angle is defined by the sequence of points P1, P2 and P3.
!
!    The point P4 is calculated so that:
!
!      (P1,P2,P4) = (P1,P2,P3) / 2
!
!        P1
!        /
!       /   P4
!      /  .
!     / .
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), P3(2), points defining the angle.
!
!    Input, real ( kind = 8 ) P4(2), a point defining the half angle.
!    The vector P4 - P2 will have unit norm.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ) p4(dim_num)

  p4(1:2) = 0.5D+00 * ( &
      ( p1(1:2) - p2(1:2) ) / sqrt ( sum ( ( p1(1:2) - p2(1:2) )**2 ) ) &
    + ( p3(1:2) - p2(1:2) ) / sqrt ( sum ( ( p3(1:2) - p2(1:2) )**2 ) ) )

  p4(1:2) = p2(1:2) + p4(1:2) / sqrt ( sum ( p4(1:2)**2 ) )

  return
end
function angle_rad_2d ( p1, p2, p3 )

!*****************************************************************************80
!
!! ANGLE_RAD_2D returns the angle swept out between two rays in 2D.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
!
!        P1
!        /
!       /
!      /
!     /
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), P3(2), define the rays
!    P1 - P2 and P3 - P2 which define the angle.
!
!    Output, real ( kind = 8 ) ANGLE_RAD_2D, the angle swept out by the rays,
!    in radians.  0 <= ANGLE_RAD_2D < 2 * PI.  If either ray has zero
!    length, then ANGLE_RAD_2D is set to 0.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) angle_rad_2d
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)

  p(1) = ( p1(1) - p2(1) ) * ( p3(1) - p2(1) ) &
       + ( p1(2) - p2(2) ) * ( p3(2) - p2(2) )

  p(2) =   ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
         - ( p1(1) - p2(1) ) * ( p3(2) - p2(2) )

  if ( p(1) == 0.0D+00 .and. p(2) == 0.0D+00 ) then

    angle_rad_2d = 0.0D+00

  else

    angle_rad_2d = atan2 ( p(2), p(1) )

    if ( angle_rad_2d < 0.0D+00 ) then
      angle_rad_2d = angle_rad_2d + 2.0D+00 * pi
    end if

  end if

  return
end
function atan4 ( y, x )

!*****************************************************************************80
!
!! ATAN4 computes the inverse tangent of the ratio Y / X.
!
!  Discussion:
!
!    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
!    the built in functions ATAN and ATAN2 already do.
!
!    However:
!
!    * ATAN4 always returns a positive angle, between 0 and 2 PI,
!      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
!      and [-PI,+PI] respectively;
!
!    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
!     function by contrast always returns an angle in the first or fourth
!     quadrants.
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
!  Parameters:
!
!    Input, real ( kind = 8 ) Y, X, two quantities which represent the
!    tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
!
!    Output, real ( kind = 8 ) ATAN4, an angle between 0 and 2 * PI, whose
!    tangent is (Y/X), and which lies in the appropriate quadrant so that
!    the signs of its cosine and sine match those of X and Y.
!
  implicit none

  real ( kind = 8 ) abs_x
  real ( kind = 8 ) abs_y
  real ( kind = 8 ) atan4
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta_0
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Special cases:
!
  if ( x == 0.0D+00 ) then

    if ( 0.0D+00 < y ) then
      theta = pi / 2.0D+00
    else if ( y < 0.0D+00 ) then
      theta = 3.0D+00 * pi / 2.0D+00
    else if ( y == 0.0D+00 ) then
      theta = 0.0D+00
    end if

  else if ( y == 0.0D+00 ) then

    if ( 0.0D+00 < x ) then
      theta = 0.0D+00
    else if ( x < 0.0D+00 ) then
      theta = pi
    end if
!
!  We assume that ATAN2 is correct when both arguments are positive.
!
  else

    abs_y = abs ( y )
    abs_x = abs ( x )

    theta_0 = atan2 ( abs_y, abs_x )

    if ( 0.0D+00 < x .and. 0.0D+00 < y ) then
      theta = theta_0
    else if ( x < 0.0D+00 .and. 0.0D+00 < y ) then
      theta = pi - theta_0
    else if ( x < 0.0D+00 .and. y < 0.0D+00 ) then
      theta = pi + theta_0
    else if ( 0.0D+00 < x .and. y < 0.0D+00 ) then
      theta = 2.0D+00 * pi - theta_0
    end if

  end if

  atan4 = theta

  return
end
function box_contains_point_2d ( box, p )

!*****************************************************************************80
!
!! BOX_CONTAINS_POINT_2D determines if a point is inside a box in 2D.
!
!  Discussion:
!
!    A unit box in 2D is the set of points (X,Y) with the property that
!
!      0.0 <= X <= 1.0
!    and
!      0.0 <= Y <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BOX(DIM_NUM,2), the lower left and upper right
!    corners of the box.
!
!    Input, real ( kind = 8 ) P(DIM_NUM), the point to be checked.
!
!    Output, logical BOX_CONTAINS_POINT_2D, is TRUE if the point is
!    inside the box.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) box(dim_num,2)
  logical box_contains_point_2d
  real ( kind = 8 ) p(dim_num)

  box_contains_point_2d = &
    all ( box(1:dim_num,1) <= p(1:dim_num) ) .and. &
    all (                     p(1:dim_num) <= box(1:dim_num,2) )

  return
end
subroutine circle_arc_point_near_2d ( r, c, theta1, theta2, x, nearest, &
  dist )

!*****************************************************************************80
!
!! CIRCLE_ARC_POINT_NEAR_2D : nearest point on a circular arc.
!
!  Discussion:
!
!    A circular arc is defined by the portion of a circle (R,C)
!    between two angles (THETA1,THETA2).
!
!    Thus, a point (X,Y) on a circular arc satisfies
!
!      ( X - C(1) ) * ( X - C(1) ) + ( Y - C(2) ) * ( Y - C(2) ) = R * R
!
!    and
!
!      Theta1 <= Theta <= Theta2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) C(2), the center of the circle.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angles defining the arc,
!    in radians.  Normally, THETA1 < THETA2.
!
!    Input, real ( kind = 8 ) X(2), the point to be checked.
!
!    Output, real ( kind = 8 ) NEAREST(2), a point on the circular arc which is
!    nearest to the point.
!
!    Output, real ( kind = 8 ) DIST, the distance to the nearest point.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) atan4
  real ( kind = 8 ) c(dim_num)
  real ( kind = 8 ) dist
  real ( kind = 8 ) nearest(dim_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) x(dim_num)
!
!  Special case, the zero circle.
!
  if ( r == 0.0D+00 ) then
    nearest(1:2) = c(1:2)
  end if
!
!  Determine the angle made by the point.
!
  theta = atan4 ( x(2) - c(2), x(1) - c(1) )
!
!  If the angle is between THETA1 and THETA2, then you can
!  simply project the point onto the arc.
!
  if ( r8_modp ( theta  - theta1,  2.0D+00 * pi ) <= &
       r8_modp ( theta2 - theta1,  2.0D+00 * pi ) ) then

    r2 = sqrt ( sum ( ( x(1:2) - c(1:2) )**2 ) )

    nearest(1:2) = c(1:2) + ( x(1:2) - c(1:2) ) * r / r2
!
!  Otherwise, if the angle is less than the negative of the
!  average of THETA1 and THETA2, it's on the side of the arc
!  where the endpoint associated with THETA2 is closest.
!
  else if ( r8_modp ( theta - 0.5D+00 * ( theta1 + theta2 ), 2.0D+00 * pi ) &
    <= pi ) then

    nearest(1:2) = c(1:2) + r * (/ cos ( theta2 ), sin ( theta2 ) /)
!
!  Otherwise, the endpoint associated with THETA1 is closest.
  else

    nearest(1:2) = c(1:2) + r * (/ cos ( theta1 ), sin ( theta1 ) /)

  end if

  dist = sqrt ( sum ( ( x(1:2) - nearest(1:2) )**2 ) )

  return
end
function circle_imp_contains_point_2d ( r, c, x )

!*****************************************************************************80
!
!! CIRCLE_IMP_CONTAINS_POINT_2D determines if a circle contains a point in 2D.
!
!  Discussion:
!
!    An implicit circle in 2D satisfies:
!
!      ( X - C(1) ) * ( X - C(1) ) + ( Y - C(2) ) * ( Y - C(2) ) = R * R
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) C(2), the coordinates of the center of the circle.
!
!    Input, real ( kind = 8 ) X(2), the point to be checked.
!
!    Output, logical CIRCLE_IMP_CONTAINS_POINT_2D, is TRUE if the point
!    is inside or on the circle, FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) c(dim_num)
  logical circle_imp_contains_point_2d
  real ( kind = 8 ) r
  real ( kind = 8 ) x(dim_num)

  if ( ( x(1) - c(1) ) * ( x(1) - c(1) ) &
     + ( x(2) - c(2) ) * ( x(2) - c(2) ) <= r * r ) then
    circle_imp_contains_point_2d = .true.
  else
    circle_imp_contains_point_2d = .false.
  end if

  return
end
subroutine circle_imp_point_near_2d ( r, c, x, nearest, dist )

!*****************************************************************************80
!
!! CIRCLE_IMP_POINT_NEAR_2D: nearest ( implicit circle, point ) in 2D.
!
!  Discussion:
!
!    This routine finds the distance from a point to an implicitly
!    defined circle, and returns the point on the circle that is
!    nearest to the given point.
!
!    If the given point is the center of the circle, than any point
!    on the circle is "the" nearest.
!
!    An implicit circle in 2D satisfies the equation:
!
!      ( X - C(1) ) * ( X - C(1) ) + ( Y - C(2) ) * ( Y - C(2) ) = R * R
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) C(2), the coordinates of the center of the circle.
!
!    Input, real ( kind = 8 ) X(2), the point to be checked.
!
!    Output, real ( kind = 8 ) NEAREST(2), the nearest point on the circle.
!
!    Output, real ( kind = 8 ) DIST, the distance of the point to the circle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) c(dim_num)
  real ( kind = 8 ) dist
  real ( kind = 8 ) nearest(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) x(dim_num)

  if ( all ( x(1:2) == c(1:2) ) ) then
    dist = r
    nearest(1:2) = c(1:2)
    nearest(1) = nearest(1) + r
    return
  end if

  r2 = sqrt ( sum ( ( x(1:2) - c(1:2) )**2 ) )

  dist = abs (  r2 - r )

  nearest(1:2) = c(1:2) + r * ( x(1:2) - c(1:2) ) / r2

  return
end
function circle_sector_contains_point_2d ( r, c, theta1, theta2, x )

!*****************************************************************************80
!
!! CIRCLE_SECTOR_CONTAINS_POINT_2D : is a point inside a circular sector?
!
!  Discussion:
!
!    A circular sector is formed by a circular arc, and the two straight line
!    segments that join its ends to the center of the circle.
!
!    A circular sector is defined by
!
!      ( X - C(1) ) * ( X - C(1) ) + ( Y - C(2) ) * ( Y - C(2) ) <= R * R
!
!    and
!
!      Theta1 <= Theta <= Theta2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) C(2), the coordinates of the center of the circle.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angles defining the arc,
!    in radians.  Normally, THETA1 < THETA2.
!
!    Input, real ( kind = 8 ) X(2), the point to be checked.
!
!    Output, logical CIRCLE_SECTOR_CONTAINS_POINT_2D, is TRUE if the point
!    is inside or on the circular sector.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) atan4
  real ( kind = 8 ) c(dim_num)
  logical circle_sector_contains_point_2d
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  logical value
  real ( kind = 8 ) x(dim_num)

  value = .false.
!
!  Is the point inside the (full) circle?
!
  if ( ( x(1) - c(1) ) * ( x(1) - c(1) ) &
     + ( x(2) - c(2) ) * ( x(2) - c(2) ) <= r * r ) then
!
!  Is the point's angle within the arc's range?
!  Try to force the angles to lie between 0 and 2 * PI.
!
    theta = atan4 ( x(2) - c(2), x(1) - c(1) )

    if ( r8_modp ( theta  - theta1,  2.0D+00 * pi ) <= &
         r8_modp ( theta2 - theta1,  2.0D+00 * pi ) ) then

      value = .true.

    end if

  end if

  circle_sector_contains_point_2d = value

  return
end
function cos_deg ( angle )

!*****************************************************************************80
!
!! COS_DEG returns the cosine of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle, in degrees.
!
!    Output, real ( kind = 8 ) COS_DEG, the cosine of the angle.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) cos_deg
  real ( kind = 8 ), parameter :: degrees_to_radians &
    = 3.141592653589793D+00 / 180.0D+00

  cos_deg = cos ( degrees_to_radians * angle )

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC generates the next filename in a series.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Example:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  integer ( kind = 4 ) digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      digit = ichar ( c ) - 48
      digit = digit + 1

      if ( digit == 10 ) then
        digit = 0
      end if

      c = char ( digit + 48 )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

  return
end
subroutine fmin_rc ( a, b, arg, status, value )

!*****************************************************************************80
!
!! FMIN_RC seeks a minimizer of a scalar function of a scalar variable.
!
!  Discussion:
!
!    FMIN_RC seeks an approximation to the point where F attains a minimum on
!    the interval (A,B).
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much
!    slower than that for a Fibonacci search.  If F has a continuous
!    second derivative which is positive at the minimum (which is not
!    at A or B), then convergence is superlinear, and usually of the
!    order of about 1.324....
!
!    The routine is a revised version of the Brent FMIN algorithm,
!    which now uses reverse communication.
!
!    It is worth stating explicitly that this routine will NOT be
!    able to detect a minimizer that occurs at either initial endpoint
!    A or B.  If this is a concern to the user, then the user must
!    either ensure that the initial interval is larger, or to check
!    the function value at the returned minimizer against the values
!    at either endpoint.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2004
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
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A, B.  On input, the left and right
!    endpoints of the initial interval.  On output, the lower and upper
!    bounds for an interval containing the minimizer.  It is required
!    that A < B.
!
!    Output, real ( kind = 8 ) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS zero, ARG is the routine's
!    estimate for the function minimizer.
!
!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between the user
!    and the routine.  The user only sets STATUS to zero on the first call,
!    to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated minimizer.
!
!    Input, real ( kind = 8 ) VALUE, the function value at ARG, as
!    requested by the routine on the previous call.
!
!  Local parameters:
!
!    C is the squared inverse of the golden ratio.
!
!    EPS is the square root of the relative machine precision.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ), save :: c
  real ( kind = 8 ), save :: d
  real ( kind = 8 ), save :: e
  real ( kind = 8 ), save :: eps
  real ( kind = 8 ), save :: fu
  real ( kind = 8 ), save :: fv
  real ( kind = 8 ), save :: fw
  real ( kind = 8 ), save :: fx
  real ( kind = 8 ), save :: midpoint
  real ( kind = 8 ), save :: p
  real ( kind = 8 ), save :: q
  real ( kind = 8 ), save :: r
  integer ( kind = 4 ) status
  real ( kind = 8 ), save :: tol
  real ( kind = 8 ), save :: tol1
  real ( kind = 8 ), save :: tol2
  real ( kind = 8 ), save :: u
  real ( kind = 8 ), save :: v
  real ( kind = 8 ) value
  real ( kind = 8 ), save :: w
  real ( kind = 8 ), save :: x
!
!  STATUS (INPUT) = 0, startup.
!
  if ( status == 0 ) then

    if ( b <= a ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FMIN_RC - Fatal error!'
      write ( *, '(a)' ) '  A < B is required, but'
      write ( *, '(a,g14.6)' ) '  A = ', a
      write ( *, '(a,g14.6)' ) '  B = ', b
      status = -1
      stop
    end if

    c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

    eps = sqrt ( epsilon ( eps ) )
    tol = epsilon ( tol )

    v = a + c * ( b - a )
    w = v
    x = v
    e = 0.0D+00

    status = 1
    arg = x

    return
!
!  STATUS (INPUT) = 1, return with initial function value of FX.
!
  else if ( status == 1 ) then

    fx = value
    fv = fx
    fw = fx
!
!  STATUS (INPUT) = 2 or more, update the data.
!
  else if ( 2 <= status ) then

    fu = value

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

    else

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

    end if

  end if
!
!  Take the next step.
!
  midpoint = 0.5D+00 * ( a + b )
  tol1 = eps * abs ( x ) + tol / 3.0D+00
  tol2 = 2.0D+00 * tol1
!
!  If the stopping criterion is satisfied, we can exit.
!
  if ( abs ( x - midpoint ) <= ( tol2 - 0.5D+00 * ( b - a ) ) ) then
    status = 0
    return
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
!
!  Request value of F(U).
!
  arg = u
  status = status + 1

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine hex_grid_angle ( box, center, angle, h, n, r )

!*****************************************************************************80
!
!! HEX_GRID_ANGLE sets the points in an angled hex grid in a box.
!
!  Discussion:
!
!    DIM_NUM = 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BOX(DIM_NUM,2), the lower left and upper right
!    corners of the box.
!
!    Input, real ( kind = 8 ) CENTER(DIM_NUM), the center of the grid.
!    This point must be inside the unit square.
!
!    Input, real ( kind = 8 ) ANGLE, the angle, in degrees, of the grid.
!    Normally, 0 <= ANGLE <= 180, but any value is allowed.
!
!    Input, real ( kind = 8 ) H, the spacing between neighboring
!    points on a grid line.
!
!    Input, integer ( kind = 4 ) N, the number of points of the angled hex grid
!    that are within the unit square.  This value may have been computed
!    by calling HEX_GRID_ANGLE_01_SIZE.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the grid points.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle2
  real ( kind = 8 ) box(dim_num,2)
  logical box_contains_point_2d
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) cos_deg
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) layer
  integer ( kind = 4 ) layer_size
  real ( kind = 8 ) point(dim_num)
  real ( kind = 8 ) r(dim_num,n)
  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) sin_deg
  integer ( kind = 4 ) size
!
!  Ninny checks.
!
  if ( .not. box_contains_point_2d ( box, center ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE - Fatal error!'
    write ( *, '(a)' ) '  The center point of the grid is not'
    write ( *, '(a)' ) '  inside the box.'
    write ( *, '(a,2g14.6)' ) '  CENTER = ', center(1:dim_num)
    stop
  end if

  if ( h == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE - Fatal error!'
    write ( *, '(a)' ) '  The grid spacing must be nonzero.'
    write ( *, '(a,g14.6)' ) '  H = ', h
    stop
  end if

  layer = 0
  point(1:dim_num) = center(1:dim_num)

  k = 1
  if ( k <= n ) then
    r(1:dim_num,k) = center(1:dim_num)
  end if

  do

    layer = layer + 1

    layer_size = 0

    angle2 = angle
!
!  Compute the first point on the new layer.
!
    point(1:dim_num) = point(1:dim_num) &
      + h * (/ cos_deg ( angle2 ), sin_deg ( angle2 ) /)

    angle2 = r8_modp ( angle2 + 60.0D+00, 360.0D+00 )

    do i = 1, 6

      angle2 = r8_modp ( angle2 + 60.0D+00, 360.0D+00 )

      do j = 1, layer

        point(1:dim_num) = point(1:dim_num) &
          + h * (/ cos_deg ( angle2 ), sin_deg ( angle2 ) /)

        if ( box_contains_point_2d ( box, point ) ) then

          layer_size = layer_size + 1
          k = k + 1

          if ( k <= n ) then
            r(1:dim_num,k) = point(1:dim_num)
          end if

        end if

      end do

    end do

    if ( layer_size == 0 ) then
      exit
    end if

  end do

  return
end
subroutine hex_grid_angle_size ( box, center, angle, h, n )

!*****************************************************************************80
!
!! HEX_GRID_ANGLE_SIZE counts the points in an angled hex grid in a box.
!
!  Discussion:
!
!    DIM_NUM = 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BOX(DIM_NUM,2), the lower left and upper right
!    corners of the box.
!
!    Input, real ( kind = 8 ) CENTER(DIM_NUM), the center of the grid.
!    This point must be inside the box
!
!    Input, real ( kind = 8 ) ANGLE, the angle, in degrees, of the grid.
!    Normally, 0 <= ANGLE <= 180, but any value is allowed.
!
!    Input, real ( kind = 8 ) H, the spacing between neighboring
!    points on a grid line.
!
!    Output, integer ( kind = 4 ) N, the number of points of the angled hex grid
!    that are within the unit square.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle2
  real ( kind = 8 ) box(dim_num,2)
  logical box_contains_point_2d
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) cos_deg
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) layer
  integer ( kind = 4 ) layer_size
  integer ( kind = 4 ) n
  real ( kind = 8 ) point(dim_num)
  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) sin_deg
!
!  Ninny checks.
!
  if ( .not. box_contains_point_2d ( box, center ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE_SIZE - Fatal error!'
    write ( *, '(a)' ) '  The center point of the grid is not'
    write ( *, '(a)' ) '  inside the box.'
    write ( *, '(a,2g14.6)' ) '  CENTER = ', center(1:dim_num)
    stop
  end if

  if ( h == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE_SIZE - Fatal error!'
    write ( *, '(a)' ) '  The grid spacing must be nonzero.'
    write ( *, '(a,g14.6)' ) '  H = ', h
    stop
  end if

  n = 0

  layer = 0
  point(1:dim_num) = center(1:dim_num)

  n = 1

  do

    layer = layer + 1

    layer_size = 0

    angle2 = angle
!
!  Compute the first point on the new layer.
!
    point(1:dim_num) = point(1:dim_num) &
      + h * (/ cos_deg ( angle2 ), sin_deg ( angle2 ) /)

    angle2 = r8_modp ( angle2 + 60.0D+00, 360.0D+00 )

    do i = 1, 6

      angle2 = r8_modp ( angle2 + 60.0D+00, 360.0D+00 )

      do j = 1, layer

        point(1:dim_num) = point(1:dim_num) &
          + h * (/ cos_deg ( angle2 ), sin_deg ( angle2 ) /)

        if ( box_contains_point_2d ( box, point ) ) then
          layer_size = layer_size + 1
          n = n + 1
        end if

      end do

    end do

    if ( layer_size == 0 ) then
      exit
    end if

  end do

  return
end
function hexagon_contains_point_2d ( v, p )

!*****************************************************************************80
!
!! HEXAGON_CONTAINS_POINT_2D finds if a point is inside a hexagon in 2D.
!
!  Discussion:
!
!    This test is only valid if the hexagon is convex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V(2,6), the vertics, in counterclockwise order.
!
!    Input, real ( kind = 8 ) P(2), the point to be tested.
!
!    Output, logical HEXAGON_CONTAINS_POINT_2D, is TRUE if X is in the hexagon.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 6

  logical hexagon_contains_point_2d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) v(dim_num,n)
!
!  A point is inside a convex hexagon if and only if it is "inside"
!  each of the 6 halfplanes defined by lines through consecutive
!  vertices.
!
  do i = 1, n

    j = mod ( i, n ) + 1

    if (  v(1,i) * ( v(2,j) - p(2  ) ) &
        + v(1,j) * ( p(2  ) - v(2,i) ) &
        + p(1  ) * ( v(2,i) - v(2,j) ) < 0.0D+00 ) then

      hexagon_contains_point_2d = .false.
      return

    end if

  end do

  hexagon_contains_point_2d = .true.

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!  Example:
!
!        I     J     MOD  I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an integer to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the integer value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    i4_wrap = jlo
  else
    i4_wrap = jlo + i4_modp ( ival - jlo, wide )
  end if

  return
end
subroutine p00_boundary_eps ( test, h, show_nodes, eps_file_name )

!*****************************************************************************80
!
!! P00_BOUNDARY_EPS draws the boundary of a region as an EPS file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the number of the problem.
!
!    Input, real ( kind = 8 ) H, the maximum size of a segment of the boundary.
!    This controls how smoothly curved sections of the boundary will be drawn.
!
!    Input, logical SHOW_NODES, is TRUE if the boundary nodes are
!    to be shown.
!
!    Input, character ( len = * ) EPS_FILE_NAME, the name of the EPS file
!    to create.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: boundary
  character ( len = * ) eps_file_name
  integer ( kind = 4 ) eps_unit
  real ( kind = 8 ) h
  real ( kind = 8 ), allocatable, dimension ( : ) :: hi
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  real ( kind = 8 ), allocatable, dimension ( : ) :: lo
  integer ( kind = 4 ) pass
  real ( kind = 8 ) scale
  integer ( kind = 4 ) segment
  integer ( kind = 4 ) segment_length
  integer ( kind = 4 ) segment_num
  logical show_nodes
  character ( len = 40 ) string
  integer ( kind = 4 ) test
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer ( kind = 4 ) x_ps
  integer ( kind = 4 ) :: x_ps_max = 576
  integer ( kind = 4 ) :: x_ps_max_clip = 594
  integer ( kind = 4 ) :: x_ps_max_user
  integer ( kind = 4 ) :: x_ps_min = 36
  integer ( kind = 4 ) :: x_ps_min_clip = 18
  integer ( kind = 4 ) :: x_ps_min_user
  real ( kind = 8 ) x_scale
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  integer ( kind = 4 ) y_ps
  integer ( kind = 4 ) :: y_ps_max = 666
  integer ( kind = 4 ) :: y_ps_max_clip = 684
  integer ( kind = 4 ) :: y_ps_max_user
  integer ( kind = 4 ) :: y_ps_min = 126
  integer ( kind = 4 ) :: y_ps_min_clip = 108
  integer ( kind = 4 ) :: y_ps_min_user
  real ( kind = 8 ) y_scale

  call p00_boundary_segment_num ( test, segment_num )

  allocate ( lo(1:m) )
  allocate ( hi(1:m) )

  call p00_box ( test, m, lo, hi )

  x_min = lo(1) - 0.025D+00 * ( hi(1) - lo(1) )
  y_min = lo(2) - 0.025D+00 * ( hi(2) - lo(2) )
  x_max = hi(1) + 0.025D+00 * ( hi(1) - lo(1) )
  y_max = hi(2) + 0.025D+00 * ( hi(2) - lo(2) )

  x_scale = x_max - x_min
  y_scale = y_max - y_min
  scale = max ( x_scale, y_scale )
!
!  Determine the PostScript coordinates of the used box.
!
  x_ps_min_user = ( ( x_ps_max + x_ps_min ) &
    - int ( real ( x_ps_max - x_ps_min, kind = 8 ) * x_scale / scale ) ) / 2
  x_ps_max_user = ( ( x_ps_max + x_ps_min ) &
    + int ( real ( x_ps_max - x_ps_min, kind = 8 ) * x_scale / scale ) ) / 2
  y_ps_min_user = ( ( y_ps_max + y_ps_min ) &
    - int ( real ( y_ps_max - y_ps_min, kind = 8 ) * y_scale / scale ) ) / 2
  y_ps_max_user = ( ( y_ps_max + y_ps_min ) &
    + int ( real ( y_ps_max - y_ps_min, kind = 8 ) * y_scale / scale ) ) / 2

  if ( x_scale < y_scale ) then
    x_max = x_max + 0.5D+00 * ( y_scale - x_scale )
    x_min = x_min - 0.5D+00 * ( y_scale - x_scale )
  else if ( y_scale < x_scale ) then
    y_max = y_max + 0.5D+00 * ( x_scale - y_scale )
    y_min = y_min - 0.5D+00 * ( x_scale - y_scale )
  end if

  call get_unit ( eps_unit )

  open ( unit = eps_unit, file = eps_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_EPS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output EPS file:'
    write ( *, '(a)' ) '    "' // trim ( eps_file_name ) // '".'
    write ( *, '(a,i12)' ) '  IOS = ', ios
    stop
  end if

  write ( eps_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( eps_unit, '(a)' ) '%%Creator: p00_boundary_eps.f90'
  write ( eps_unit, '(a)' ) '%%Title: ' // trim ( eps_file_name )
  write ( eps_unit, '(a)' ) '%%Pages: 1'
  write ( eps_unit, '(a,i4,2x,i4,2x,i4,2x,i4)' ) '%%BoundingBox: ', &
    x_ps_min_user, y_ps_min_user, x_ps_max_user, y_ps_max_user
  write ( eps_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( eps_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( eps_unit, '(a)' ) '%%EndComments'
  write ( eps_unit, '(a)' ) '%%BeginProlog'
  write ( eps_unit, '(a)' ) '/inch {72 mul} def'
  write ( eps_unit, '(a)' ) '%%EndProlog'
  write ( eps_unit, '(a)' ) '%%Page: 1 1'
  write ( eps_unit, '(a)' ) 'save'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw a gray border around the user box.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'newpath'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) &
    '  ', x_ps_min_user, y_ps_min_user, ' moveto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) &
    '  ', x_ps_max_user, y_ps_min_user, ' lineto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) &
    '  ', x_ps_max_user, y_ps_max_user, ' lineto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) &
    '  ', x_ps_min_user, y_ps_max_user, ' lineto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) &
    '  ', x_ps_min_user, y_ps_min_user, ' lineto'
  write ( eps_unit, '(a)' ) 'stroke'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Set the RGB line color to black.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Set the font and its size.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '/Times-Roman findfont'
  write ( eps_unit, '(a)' ) '0.50 inch scalefont'
  write ( eps_unit, '(a)' ) 'setfont'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Define a clipping polygon.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'newpath'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' moveto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_min_clip, ' lineto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_max_clip, ' lineto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_max_clip, ' lineto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' lineto'
  write ( eps_unit, '(a)' ) 'clip newpath'

  do segment = 1, segment_num

    call p00_boundary_segment_length ( test, segment, h, segment_length )

    allocate ( boundary(1:m,1:segment_length) )

    call p00_boundary_segment ( test, segment, m, segment_length, boundary )

    if ( show_nodes ) then

      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '%  Set the RGB line color to green.'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '0.000  0.750  0.150 setrgbcolor'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '%  Draw the nodes.'
      write ( eps_unit, '(a)' ) '%'

      do j = 1, segment_length

        x_ps = int ( &
          ( ( x_max - boundary(1,j)         ) * real ( x_ps_min, kind = 8 )   &
          + (         boundary(1,j) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                 - x_min ) )

        y_ps = int ( &
          ( ( y_max - boundary(2,j)         ) * real ( y_ps_min, kind = 8 )   &
          + (         boundary(2,j) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                 - y_min ) )

        write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
          ' 5 0 360 arc closepath fill'

      end do
    end if

    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '%  Set the RGB line color to red.'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '%  Increase the linewidth to 3.'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '3  setlinewidth'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '%  Draw the boundary lines.'
    write ( eps_unit, '(a)' ) '%'

    do j = 1, segment_length
      x_ps = int ( &
        ( ( x_max - boundary(1,j)         ) * real ( x_ps_min, kind = 8 )   &
        + (         boundary(1,j) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                 - x_min ) )

      y_ps = int ( &
        ( ( y_max - boundary(2,j)         ) * real ( y_ps_min, kind = 8 )   &
        + (         boundary(2,j) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                 - y_min ) )

      if ( j == 1 ) then
        write ( eps_unit, '(i4,2x,i4,2x,a)' ) x_ps, y_ps, ' moveto'
      else
        write ( eps_unit, '(i4,2x,i4,2x,a)' ) x_ps, y_ps, ' lineto'
      end if

    end do

    write ( eps_unit, '(a)' ) 'stroke'

    deallocate ( boundary )

  end do

  deallocate ( hi )
  deallocate ( lo )

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'restore showpage'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% End of page'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%%Trailer'
  write ( eps_unit, '(a)' ) '%%EOF'

  close ( unit = eps_unit )

  if ( .false. ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_EPS:'
    write ( *, '(a)' ) '  An encapsulated PostScript file was created'
    write ( *, '(a)' ) '  containing an image of the boundary.'
    write ( *, '(a)' ) '  The file is named "' // trim ( eps_file_name ) // '".'
  end if

  return
end
subroutine p00_boundary_nearest ( test, dim_num, n, point, boundary )

!*****************************************************************************80
!
!! P00_BOUNDARY_NEAREST returns a nearest boundary point for any problem.
!
!  Discussion:
!
!    The given input point need not be inside the region.
!
!    In some cases, more than one boundary point may be "nearest",
!    but only one will be returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(DIM_NUM,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( dim_num, n ) :: boundary
  real ( kind = 8 ), dimension ( dim_num, n ) :: point
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 2 ) then
    call p02_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 3 ) then
    call p03_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 4 ) then
    call p04_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 5 ) then
    call p05_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 6 ) then
    call p06_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 7 ) then
    call p07_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 8 ) then
    call p08_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 9 ) then
    call p09_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 10 ) then
    call p10_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 11 ) then
    call p11_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 12 ) then
    call p12_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 13 ) then
    call p13_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 14 ) then
    call p14_boundary_nearest ( dim_num, n, point, boundary )
  else if ( test == 15 ) then
    call p15_boundary_nearest ( dim_num, n, point, boundary )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_NEAREST - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_boundary_project ( test, dim_num, n, point )

!*****************************************************************************80
!
!! P00_BOUNDARY_PROJECT projects exterior points to the boundary.
!
!  Discussion:
!
!    Interior points are not changed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(DIM_NUM,N), the coordinates
!    of the points.  Any input points that are exterior to the region
!    are replaced on output by the nearest boundary point.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( dim_num, n ) :: point
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_boundary_project ( dim_num, n, point )
  else if ( test == 2 ) then
    call p02_boundary_project ( dim_num, n, point )
  else if ( test == 3 ) then
    call p03_boundary_project ( dim_num, n, point )
  else if ( test == 4 ) then
    call p04_boundary_project ( dim_num, n, point )
  else if ( test == 5 ) then
    call p05_boundary_project ( dim_num, n, point )
  else if ( test == 6 ) then
    call p06_boundary_project ( dim_num, n, point )
  else if ( test == 7 ) then
    call p07_boundary_project ( dim_num, n, point )
  else if ( test == 8 ) then
    call p08_boundary_project ( dim_num, n, point )
  else if ( test == 9 ) then
    call p09_boundary_project ( dim_num, n, point )
  else if ( test == 10 ) then
    call p10_boundary_project ( dim_num, n, point )
  else if ( test == 11 ) then
    call p11_boundary_project ( dim_num, n, point )
  else if ( test == 12 ) then
    call p12_boundary_project ( dim_num, n, point )
  else if ( test == 13 ) then
    call p13_boundary_project ( dim_num, n, point )
  else if ( test == 14 ) then
    call p14_boundary_project ( dim_num, n, point )
  else if ( test == 15 ) then
    call p15_boundary_project ( dim_num, n, point )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_PROJECT - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_boundary_segment ( test, segment_index, m, &
  segment_length, segment )

!*****************************************************************************80
!
!! P00_BOUNDARY_SEGMENT returns a boundary segment in a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) segment(m,segment_length)
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 2 ) then
    call p02_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 3 ) then
    call p03_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 4 ) then
    call p04_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 5 ) then
    call p05_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 6 ) then
    call p06_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 7 ) then
    call p07_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 8 ) then
    call p08_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 9 ) then
    call p09_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 10 ) then
    call p10_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 11 ) then
    call p11_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 12 ) then
    call p12_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 13 ) then
    call p13_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 14 ) then
    call p14_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else if ( test == 15 ) then
    call p15_boundary_segment ( segment_index, m, segment_length, &
      segment )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_boundary_segment_length ( test, segment_index, h, &
  segment_length )

!*****************************************************************************80
!
!! P00_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the
!    boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points
!    in the segment.
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 2 ) then
    call p02_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 3 ) then
    call p03_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 4 ) then
    call p04_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 5 ) then
    call p05_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 6 ) then
    call p06_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 7 ) then
    call p07_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 8 ) then
    call p08_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 9 ) then
    call p09_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 10 ) then
    call p10_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 11 ) then
    call p11_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 12 ) then
    call p12_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 13 ) then
    call p13_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 14 ) then
    call p14_boundary_segment_length ( segment_index, h, segment_length )
  else if ( test == 15 ) then
    call p15_boundary_segment_length ( segment_index, h, segment_length )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_boundary_segment_num ( test, boundary_segment_num )

!*****************************************************************************80
!
!! P00_BOUNDARY_SEGMENT_NUM counts the boundary segments in a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_boundary_segment_num ( boundary_segment_num )
  else if ( test == 2 ) then
    call p02_boundary_segment_num ( boundary_segment_num )
  else if ( test == 3 ) then
    call p03_boundary_segment_num ( boundary_segment_num )
  else if ( test == 4 ) then
    call p04_boundary_segment_num ( boundary_segment_num )
  else if ( test == 5 ) then
    call p05_boundary_segment_num ( boundary_segment_num )
  else if ( test == 6 ) then
    call p06_boundary_segment_num ( boundary_segment_num )
  else if ( test == 7 ) then
    call p07_boundary_segment_num ( boundary_segment_num )
  else if ( test == 8 ) then
    call p08_boundary_segment_num ( boundary_segment_num )
  else if ( test == 9 ) then
    call p09_boundary_segment_num ( boundary_segment_num )
  else if ( test == 10 ) then
    call p10_boundary_segment_num ( boundary_segment_num )
  else if ( test == 11 ) then
    call p11_boundary_segment_num ( boundary_segment_num )
  else if ( test == 12 ) then
    call p12_boundary_segment_num ( boundary_segment_num )
  else if ( test == 13 ) then
    call p13_boundary_segment_num ( boundary_segment_num )
  else if ( test == 14 ) then
    call p14_boundary_segment_num ( boundary_segment_num )
  else if ( test == 15 ) then
    call p15_boundary_segment_num ( boundary_segment_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_SEGMENT_NUM - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_box ( test, m, lo, hi )

!*****************************************************************************80
!
!! P00_BOX returns a bounding box for a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), the lower and
!    upper corners of a bounding box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_box ( m, lo, hi )
  else if ( test == 2 ) then
    call p02_box ( m, lo, hi )
  else if ( test == 3 ) then
    call p03_box ( m, lo, hi )
  else if ( test == 4 ) then
    call p04_box ( m, lo, hi )
  else if ( test == 5 ) then
    call p05_box ( m, lo, hi )
  else if ( test == 6 ) then
    call p06_box ( m, lo, hi )
  else if ( test == 7 ) then
    call p07_box ( m, lo, hi )
  else if ( test == 8 ) then
    call p08_box ( m, lo, hi )
  else if ( test == 9 ) then
    call p09_box ( m, lo, hi )
  else if ( test == 10 ) then
    call p10_box ( m, lo, hi )
  else if ( test == 11 ) then
    call p11_box ( m, lo, hi )
  else if ( test == 12 ) then
    call p12_box ( m, lo, hi )
  else if ( test == 13 ) then
    call p13_box ( m, lo, hi )
  else if ( test == 14 ) then
    call p14_box ( m, lo, hi )
  else if ( test == 15 ) then
    call p15_box ( m, lo, hi )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOX - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_density ( test, m, n, point, density )

!*****************************************************************************80
!
!! P00_DENSITY returns the density for a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density
!    at each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_density ( m, n, point, density )
  else if ( test == 2 ) then
    call p02_density ( m, n, point, density )
  else if ( test == 3 ) then
    call p03_density ( m, n, point, density )
  else if ( test == 4 ) then
    call p04_density ( m, n, point, density )
  else if ( test == 5 ) then
    call p05_density ( m, n, point, density )
  else if ( test == 6 ) then
    call p06_density ( m, n, point, density )
  else if ( test == 7 ) then
    call p07_density ( m, n, point, density )
  else if ( test == 8 ) then
    call p08_density ( m, n, point, density )
  else if ( test == 9 ) then
    call p09_density ( m, n, point, density )
  else if ( test == 10 ) then
    call p10_density ( m, n, point, density )
  else if ( test == 11 ) then
    call p11_density ( m, n, point, density )
  else if ( test == 12 ) then
    call p12_density ( m, n, point, density )
  else if ( test == 13 ) then
    call p13_density ( m, n, point, density )
  else if ( test == 14 ) then
    call p14_density ( m, n, point, density )
  else if ( test == 15 ) then
    call p15_density ( m, n, point, density )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_DENSITY - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_element_size ( test, element_size )

!*****************************************************************************80
!
!! P00_ELEMENT_SIZE returns a typical element size for a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_element_size ( element_size )
  else if ( test == 2 ) then
    call p02_element_size ( element_size )
  else if ( test == 3 ) then
    call p03_element_size ( element_size )
  else if ( test == 4 ) then
    call p04_element_size ( element_size )
  else if ( test == 5 ) then
    call p05_element_size ( element_size )
  else if ( test == 6 ) then
    call p06_element_size ( element_size )
  else if ( test == 7 ) then
    call p07_element_size ( element_size )
  else if ( test == 8 ) then
    call p08_element_size ( element_size )
  else if ( test == 9 ) then
    call p09_element_size ( element_size )
  else if ( test == 10 ) then
    call p10_element_size ( element_size )
  else if ( test == 11 ) then
    call p11_element_size ( element_size )
  else if ( test == 12 ) then
    call p12_element_size ( element_size )
  else if ( test == 13 ) then
    call p13_element_size ( element_size )
  else if ( test == 14 ) then
    call p14_element_size ( element_size )
  else if ( test == 15 ) then
    call p15_element_size ( element_size )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_ELEMENT_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_fixed_num ( test, fixed_num )

!*****************************************************************************80
!
!! P00_FIXED_NUM returns the number of fixed points in a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_fixed_num ( fixed_num )
  else if ( test == 2 ) then
    call p02_fixed_num ( fixed_num )
  else if ( test == 3 ) then
    call p03_fixed_num ( fixed_num )
  else if ( test == 4 ) then
    call p04_fixed_num ( fixed_num )
  else if ( test == 5 ) then
    call p05_fixed_num ( fixed_num )
  else if ( test == 6 ) then
    call p06_fixed_num ( fixed_num )
  else if ( test == 7 ) then
    call p07_fixed_num ( fixed_num )
  else if ( test == 8 ) then
    call p08_fixed_num ( fixed_num )
  else if ( test == 9 ) then
    call p09_fixed_num ( fixed_num )
  else if ( test == 10 ) then
    call p10_fixed_num ( fixed_num )
  else if ( test == 11 ) then
    call p11_fixed_num ( fixed_num )
  else if ( test == 12 ) then
    call p12_fixed_num ( fixed_num )
  else if ( test == 13 ) then
    call p13_fixed_num ( fixed_num )
  else if ( test == 14 ) then
    call p14_fixed_num ( fixed_num )
  else if ( test == 15 ) then
    call p15_fixed_num ( fixed_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FIXED_NUM - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_fixed_points ( test, m, fixed_num, fixed )

!*****************************************************************************80
!
!! P00_FIXED_POINTS returns the fixed points in a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the
!    coordinates of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) fixed(m,fixed_num)
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_fixed_points ( m, fixed_num, fixed )
  else if ( test == 2 ) then
    call p02_fixed_points ( m, fixed_num, fixed )
  else if ( test == 3 ) then
    call p03_fixed_points ( m, fixed_num, fixed )
  else if ( test == 4 ) then
    call p04_fixed_points ( m, fixed_num, fixed )
  else if ( test == 5 ) then
    call p05_fixed_points ( m, fixed_num, fixed )
  else if ( test == 6 ) then
    call p06_fixed_points ( m, fixed_num, fixed )
  else if ( test == 7 ) then
    call p07_fixed_points ( m, fixed_num, fixed )
  else if ( test == 8 ) then
    call p08_fixed_points ( m, fixed_num, fixed )
  else if ( test == 9 ) then
    call p09_fixed_points ( m, fixed_num, fixed )
  else if ( test == 10 ) then
    call p10_fixed_points ( m, fixed_num, fixed )
  else if ( test == 11 ) then
    call p11_fixed_points ( m, fixed_num, fixed )
  else if ( test == 12 ) then
    call p12_fixed_points ( m, fixed_num, fixed )
  else if ( test == 13 ) then
    call p13_fixed_points ( m, fixed_num, fixed )
  else if ( test == 14 ) then
    call p14_fixed_points ( m, fixed_num, fixed )
  else if ( test == 15 ) then
    call p15_fixed_points ( m, fixed_num, fixed )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FIXED_POINTS - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_header ( test )

!*****************************************************************************80
!
!! P00_HEADER prints some information about a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
  implicit none

  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_header ( )
  else if ( test == 2 ) then
    call p02_header ( )
  else if ( test == 3 ) then
    call p03_header ( )
  else if ( test == 4 ) then
    call p04_header ( )
  else if ( test == 5 ) then
    call p05_header ( )
  else if ( test == 6 ) then
    call p06_header ( )
  else if ( test == 7 ) then
    call p07_header ( )
  else if ( test == 8 ) then
    call p08_header ( )
  else if ( test == 9 ) then
    call p09_header ( )
  else if ( test == 10 ) then
    call p10_header ( )
  else if ( test == 11 ) then
    call p11_header ( )
  else if ( test == 12 ) then
    call p12_header ( )
  else if ( test == 13 ) then
    call p13_header ( )
  else if ( test == 14 ) then
    call p14_header ( )
  else if ( test == 15 ) then
    call p15_header ( )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_hex_grid ( test, m, h, n, point )

!*****************************************************************************80
!
!! P00_HEX_GRID returns hex grid points in a region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, real ( kind = 8 ) H, the hexagonal spacing.
!
!    Input, integer ( kind = 4 ) N, the number of hexagonal grid points
!    that lie within the region, as computed by P00_HEX_GRID_COUNT.
!
!    Output, real POINT(M,N), the hex grid points.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) angle
  real ( kind = 8 ) box(m,2)
  real ( kind = 8 ) center(m)
  real ( kind = 8 ) h
  real ( kind = 8 ) hi(m)
  logical, allocatable, dimension ( : ) :: inside
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  real ( kind = 8 ) lo(m)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point2
  integer ( kind = 4 ) test
!
!  Get the box limits.
!
  call p00_box ( test, m, lo, hi )
!
!  How many hex points will fit in the box?
!
  box(1:2,1) = lo(1:2)
  box(1:2,2) = hi(1:2)
  center(1:2) = 0.5D+00 * ( lo(1:2) + hi(1:2) )
  angle = 0.0D+00

  call hex_grid_angle_size ( box, center, angle, h, n2 )
!
!  Generate the hex points.
!
  allocate ( inside(1:n2) )
  allocate ( point2(1:m,1:n2) )

  call hex_grid_angle ( box, center, angle, h, n2, point2 )
!
!  How many of these points are in the region?
!
  call p00_inside ( test, m, n2, point2, inside )
!
!  Copy the good hex grid points.
!
  j = 0
  do j2 = 1, n2
    if ( inside(j2) ) then
      j = j + 1
      point(1:m,j) = point2(1:m,j2)
    end if
  end do

  deallocate ( inside )
  deallocate ( point2 )

  return
end
subroutine p00_hex_grid_count ( test, m, h, n )

!*****************************************************************************80
!
!! P00_HEX_GRID_COUNT counts the number of hex grid points in a region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, real ( kind = 8 ) H, the hexagonal spacing.
!
!    Output, integer ( kind = 4 ) N, the number of hexagonal grid points
!    that lie within the region.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) angle
  real ( kind = 8 ) box(m,2)
  real ( kind = 8 ) center(m)
  real ( kind = 8 ) h
  real ( kind = 8 ) hi(m)
  logical, allocatable, dimension ( : ) :: inside
  integer ( kind = 4 ) j
  real ( kind = 8 ) lo(m)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point2
  integer ( kind = 4 ) test
!
!  Get the box limits.
!
  call p00_box ( test, m, lo, hi )
!
!  How many hex points will fit in the box?
!
  box(1:2,1) = lo(1:2)
  box(1:2,2) = hi(1:2)
  center(1:2) = 0.5D+00 * ( lo(1:2) + hi(1:2) )
  angle = 0.0D+00

  call hex_grid_angle_size ( box, center, angle, h, n2 )
!
!  Generate the hex points.
!
  allocate ( inside(1:n2) )
  allocate ( point2(1:m,1:n2) )

  call hex_grid_angle ( box, center, angle, h, n2, point2 )
!
!  How many of these points are in the region?
!
  call p00_inside ( test, m, n2, point2, inside )
!
!  Add them up.
!
  n = 0
  do j = 1, n2
    if ( inside(j) ) then
      n = n + 1
    end if
  end do

  deallocate ( inside )
  deallocate ( point2 )

  return
end
subroutine p00_hole_num ( test, hole_num )

!*****************************************************************************80
!
!! P00_HOLE_NUM counts the holes in a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of boundary segments.
!
  implicit none

  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_hole_num ( hole_num )
  else if ( test == 2 ) then
    call p02_hole_num ( hole_num )
  else if ( test == 3 ) then
    call p03_hole_num ( hole_num )
  else if ( test == 4 ) then
    call p04_hole_num ( hole_num )
  else if ( test == 5 ) then
    call p05_hole_num ( hole_num )
  else if ( test == 6 ) then
    call p06_hole_num ( hole_num )
  else if ( test == 7 ) then
    call p07_hole_num ( hole_num )
  else if ( test == 8 ) then
    call p08_hole_num ( hole_num )
  else if ( test == 9 ) then
    call p09_hole_num ( hole_num )
  else if ( test == 10 ) then
    call p10_hole_num ( hole_num )
  else if ( test == 11 ) then
    call p11_hole_num ( hole_num )
  else if ( test == 12 ) then
    call p12_hole_num ( hole_num )
  else if ( test == 13 ) then
    call p13_hole_num ( hole_num )
  else if ( test == 14 ) then
    call p14_hole_num ( hole_num )
  else if ( test == 15 ) then
    call p15_hole_num ( hole_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_HOLE_NUM - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    write ( *, '(a,i8)' ) '  TEST = ', test
    stop
  end if

  return
end
subroutine p00_hole_point ( test, hole_index, m, hole_point )

!*****************************************************************************80
!
!! P00_HOLE_POINT returns a point inside a given hole.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_hole_point ( hole_index, m, hole_point )
  else if ( test == 2 ) then
    call p02_hole_point ( hole_index, m, hole_point )
  else if ( test == 3 ) then
    call p03_hole_point ( hole_index, m, hole_point )
  else if ( test == 4 ) then
    call p04_hole_point ( hole_index, m, hole_point )
  else if ( test == 5 ) then
    call p05_hole_point ( hole_index, m, hole_point )
  else if ( test == 6 ) then
    call p06_hole_point ( hole_index, m, hole_point )
  else if ( test == 7 ) then
    call p07_hole_point ( hole_index, m, hole_point )
  else if ( test == 8 ) then
    call p08_hole_point ( hole_index, m, hole_point )
  else if ( test == 9 ) then
    call p09_hole_point ( hole_index, m, hole_point )
  else if ( test == 10 ) then
    call p10_hole_point ( hole_index, m, hole_point )
  else if ( test == 11 ) then
    call p11_hole_point ( hole_index, m, hole_point )
  else if ( test == 12 ) then
    call p12_hole_point ( hole_index, m, hole_point )
  else if ( test == 13 ) then
    call p13_hole_point ( hole_index, m, hole_point )
  else if ( test == 14 ) then
    call p14_hole_point ( hole_index, m, hole_point )
  else if ( test == 15 ) then
    call p15_hole_point ( hole_index, m, hole_point )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_HOLE_POINT - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_inside ( test, m, n, point, inside )

!*****************************************************************************80
!
!! P00_INSIDE reports if a point is inside the region in a problem.
!
!  Discussion:
!
!    For argument's sake, a point on the boundary can be considered
!    inside the region, but it is probably futile to attempt to distinguish
!    this case in general.  For more information about a point's location,
!    use P00_SDIST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical inside(n)
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_inside ( m, n, point, inside )
  else if ( test == 2 ) then
    call p02_inside ( m, n, point, inside )
  else if ( test == 3 ) then
    call p03_inside ( m, n, point, inside )
  else if ( test == 4 ) then
    call p04_inside ( m, n, point, inside )
  else if ( test == 5 ) then
    call p05_inside ( m, n, point, inside )
  else if ( test == 6 ) then
    call p06_inside ( m, n, point, inside )
  else if ( test == 7 ) then
    call p07_inside ( m, n, point, inside )
  else if ( test == 8 ) then
    call p08_inside ( m, n, point, inside )
  else if ( test == 9 ) then
    call p09_inside ( m, n, point, inside )
  else if ( test == 10 ) then
    call p10_inside ( m, n, point, inside )
  else if ( test == 11 ) then
    call p11_inside ( m, n, point, inside )
  else if ( test == 12 ) then
    call p12_inside ( m, n, point, inside )
  else if ( test == 13 ) then
    call p13_inside ( m, n, point, inside )
  else if ( test == 14 ) then
    call p14_inside ( m, n, point, inside )
  else if ( test == 15 ) then
    call p15_inside ( m, n, point, inside )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_INSIDE - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_points_eps ( test, h, m, n, points, eps_file_name )

!*****************************************************************************80
!
!! P00_POINTS_EPS draws points in a region as an EPS file.
!
!  Discussion:
!
!    The boundary of the region is also drawn.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the number of the problem.
!
!    Input, real ( kind = 8 ) H, the maximum size of a segment of the boundary.
!    This controls how smoothly curved sections of the boundary will be drawn.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) POINTS(M,N), the points to draw.
!
!    Input, character ( len = * ) EPS_FILE_NAME, the name of the EPS file
!    to create.
!
!  Local Parameters:
!
!    Local, integer CIRCLE_SIZE, controls the size of the circles
!    used to indicate points.  These are measured in PostScript points,
!    that is, 1/72 of an inch.  A value of 3, 4 or 5 is usually reasonable.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: boundary
  integer ( kind = 4 ), parameter :: circle_size = 4
  character ( len = * ) eps_file_name
  integer ( kind = 4 ) eps_unit
  real ( kind = 8 ) h
  real ( kind = 8 ), allocatable, dimension ( : ) :: hi
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  real ( kind = 8 ), allocatable, dimension ( : ) :: lo
  real ( kind = 8 ) points(m,n)
  real ( kind = 8 ) scale
  integer ( kind = 4 ) segment
  integer ( kind = 4 ) segment_length
  integer ( kind = 4 ) segment_num
  logical, parameter :: show_points = .false.
  character ( len = 40 ) string
  integer ( kind = 4 ) test
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer ( kind = 4 ) x_ps
  integer ( kind = 4 ) :: x_ps_max = 576
  integer ( kind = 4 ) :: x_ps_max_clip = 594
  integer ( kind = 4 ) :: x_ps_max_user
  integer ( kind = 4 ) :: x_ps_min = 36
  integer ( kind = 4 ) :: x_ps_min_clip = 18
  integer ( kind = 4 ) :: x_ps_min_user
  real ( kind = 8 ) x_scale
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  integer ( kind = 4 ) y_ps
  integer ( kind = 4 ) :: y_ps_max = 666
  integer ( kind = 4 ) :: y_ps_max_clip = 684
  integer ( kind = 4 ) :: y_ps_max_user
  integer ( kind = 4 ) :: y_ps_min = 126
  integer ( kind = 4 ) :: y_ps_min_clip = 108
  integer ( kind = 4 ) :: y_ps_min_user
  real ( kind = 8 ) y_scale

  call p00_boundary_segment_num ( test, segment_num )

  allocate ( lo(1:m) )
  allocate ( hi(1:m) )

  call p00_box ( test, m, lo, hi )

  x_min = lo(1) - 0.025D+00 * ( hi(1) - lo(1) )
  y_min = lo(2) - 0.025D+00 * ( hi(2) - lo(2) )
  x_max = hi(1) + 0.025D+00 * ( hi(1) - lo(1) )
  y_max = hi(2) + 0.025D+00 * ( hi(2) - lo(2) )

  x_min = min ( x_min, minval ( points(1,1:n) ) )
  x_max = max ( x_max, maxval ( points(1,1:n) ) )
  y_min = min ( y_min, minval ( points(2,1:n) ) )
  y_max = max ( y_max, maxval ( points(2,1:n) ) )

  x_scale = x_max - x_min
  y_scale = y_max - y_min
  scale = max ( x_scale, y_scale )
!
!  Determine the PostScript coordinates of the used box.
!
  x_ps_min_user = ( ( x_ps_max + x_ps_min ) &
    - int ( real ( x_ps_max - x_ps_min, kind = 8 ) * x_scale / scale ) ) / 2
  x_ps_max_user = ( ( x_ps_max + x_ps_min ) &
    + int ( real ( x_ps_max - x_ps_min, kind = 8 ) * x_scale / scale ) ) / 2
  y_ps_min_user = ( ( y_ps_max + y_ps_min ) &
    - int ( real ( y_ps_max - y_ps_min, kind = 8 ) * y_scale / scale ) ) / 2
  y_ps_max_user = ( ( y_ps_max + y_ps_min ) &
    + int ( real ( y_ps_max - y_ps_min, kind = 8 ) * y_scale / scale ) ) / 2

  if ( x_scale < y_scale ) then
    x_max = x_max + 0.5D+00 * ( y_scale - x_scale )
    x_min = x_min - 0.5D+00 * ( y_scale - x_scale )
  else if ( y_scale < x_scale ) then
    y_max = y_max + 0.5D+00 * ( x_scale - y_scale )
    y_min = y_min - 0.5D+00 * ( x_scale - y_scale )
  end if

  call get_unit ( eps_unit )

  open ( unit = eps_unit, file = eps_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_POINTS_EPS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output EPS file.'
    stop
  end if

  write ( eps_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( eps_unit, '(a)' ) '%%Creator: triangulation_plot_eps.f90'
  write ( eps_unit, '(a)' ) '%%Title: ' // trim ( eps_file_name )
  write ( eps_unit, '(a)' ) '%%Pages: 1'
  write ( eps_unit, '(a,i4,2x,i4,2x,i4,2x,i4)' ) '%%BoundingBox: ', &
    x_ps_min_user, y_ps_min_user, x_ps_max_user, y_ps_max_user
  write ( eps_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( eps_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( eps_unit, '(a)' ) '%%EndComments'
  write ( eps_unit, '(a)' ) '%%BeginProlog'
  write ( eps_unit, '(a)' ) '/inch {72 mul} def'
  write ( eps_unit, '(a)' ) '%%EndProlog'
  write ( eps_unit, '(a)' ) '%%Page: 1 1'
  write ( eps_unit, '(a)' ) 'save'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw a gray border around the user box.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'newpath'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) &
    '  ', x_ps_min_user, y_ps_min_user, ' moveto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) &
    '  ', x_ps_max_user, y_ps_min_user, ' lineto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) &
    '  ', x_ps_max_user, y_ps_max_user, ' lineto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) &
    '  ', x_ps_min_user, y_ps_max_user, ' lineto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) &
    '  ', x_ps_min_user, y_ps_min_user, ' lineto'
  write ( eps_unit, '(a)' ) 'stroke'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Set the RGB line color to black.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Set the font and its size.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '/Times-Roman findfont'
  write ( eps_unit, '(a)' ) '0.50 inch scalefont'
  write ( eps_unit, '(a)' ) 'setfont'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Define a clipping polygon.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'newpath'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' moveto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_min_clip, ' lineto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_max_clip, ' lineto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_max_clip, ' lineto'
  write ( eps_unit, '(a,i4,2x,i4,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' lineto'
  write ( eps_unit, '(a)' ) 'clip newpath'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Set the RGB line color to red.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Increase the linewidth to 3.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '3  setlinewidth'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw the boundary lines.'
  write ( eps_unit, '(a)' ) '%'

  do segment = 1, segment_num

    call p00_boundary_segment_length ( test, segment, h, segment_length )

    if ( segment_length <= 0 ) then
      cycle
    end if

    allocate ( boundary(1:m,1:segment_length) )

    call p00_boundary_segment ( test, segment, m, segment_length, boundary )

    do j = 1, segment_length

      x_ps = int ( &
        ( ( x_max - boundary(1,j)         ) * real ( x_ps_min, kind = 8 )   &
        + (         boundary(1,j) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                 - x_min ) )

      y_ps = int ( &
        ( ( y_max - boundary(2,j)         ) * real ( y_ps_min, kind = 8 )   &
        + (         boundary(2,j) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                 - y_min ) )

      if ( j == 1 ) then
        write ( eps_unit, '(i4,2x,i4,2x,a)' ) x_ps, y_ps, ' moveto'
      else
        write ( eps_unit, '(i4,2x,i4,2x,a)' ) x_ps, y_ps, ' lineto'
      end if

    end do

    write ( eps_unit, '(a)' ) 'stroke'

    deallocate ( boundary )

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Set the RGB line color to blue.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Restore the linewidth to 1.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '1  setlinewidth'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw the nodes.'
  write ( eps_unit, '(a)' ) '%'

  do j = 1, n

    x_ps = int ( &
      ( ( x_max - points(1,j)         ) * real ( x_ps_min, kind = 8 )   &
      + (         points(1,j) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
      / ( x_max               - x_min ) )

    y_ps = int ( &
      ( ( y_max - points(2,j)         ) * real ( y_ps_min, kind = 8 )   &
      + (         points(2,j) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
      / ( y_max               - y_min ) )

    write ( eps_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
      circle_size, ' 0 360 arc closepath fill'

  end do

  deallocate ( hi )
  deallocate ( lo )

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'restore showpage'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% End of page'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%%Trailer'
  write ( eps_unit, '(a)' ) '%%EOF'

  close ( unit = eps_unit )

  if ( .false. ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_POINTS_EPS:'
    write ( *, '(a)' ) '  An encapsulated PostScript file was created'
    write ( *, '(a)' ) '  containing an image of the points.'
    write ( *, '(a)' ) '  The file is named "' // trim ( eps_file_name ) // '".'
  end if

  return
end
subroutine p00_poly_write ( test, file_name )

!*****************************************************************************80
!
!! P00_POLY_WRITE collects data and writes it to a POLY file.
!
!  Discussion:
!
!    This routine collects information about the boundary for a given
!    problem, and writes that data to a POLY file that can be read by
!    TRIANGLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to create.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) boundary_segment_num
  integer ( kind = 4 ) e
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: edge_nodes
  integer ( kind = 4 ) edge_num
  character ( len = * ) file_name
  real ( kind = 8 ) h
  real ( kind = 8 ) hi(dim_num)
  integer ( kind = 4 ) hole_index
  integer ( kind = 4 ) hole_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: hole
  integer ( kind = 4 ) j
  real ( kind = 8 ) lo(dim_num)
  integer ( kind = 4 ) next
  integer ( kind = 4 ) node_num
  real ( kind = 8 ) scale
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: segment
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length
  integer ( kind = 4 ) test
!
!  Get a length scale.
!
  call p00_box ( test, dim_num, lo, hi )

  scale = maxval ( abs ( hi(1:dim_num) - lo(1:dim_num) ) )
!
!  How many boundary segments are there?
!
  call p00_boundary_segment_num ( test, boundary_segment_num )
!
!  Choose H so that we would expect to get about 100 boundary points if our
!  region were a square of any size.
!
  h = 0.04D+00 * scale

  node_num = 0
  do segment_index = 1, boundary_segment_num
    call p00_boundary_segment_length ( test, segment_index, h, segment_length )
    node_num = node_num + segment_length
  end do

  allocate ( segment(dim_num,node_num) )
!
!  Now collect all the boundary nodes into one array.
!
  next = 1

  do segment_index = 1, boundary_segment_num

    call p00_boundary_segment_length ( test, segment_index, h, segment_length )

    call p00_boundary_segment ( test, segment_index, dim_num, &
      segment_length, segment(1:dim_num,next:next+segment_length-1) )

    next = next + segment_length

  end do
!
!  How many edges are there?
!
  edge_num = 0

  do segment_index = 1, boundary_segment_num

    call p00_boundary_segment_length ( test, segment_index, h, segment_length )

    edge_num = edge_num + segment_length - 1

  end do

  allocate ( edge_nodes(2,edge_num) )
!
!  Now collect the edges.
!
  e = 0
  next = 1

  do segment_index = 1, boundary_segment_num

    call p00_boundary_segment_length ( test, segment_index, h, segment_length )

    do j = 1, segment_length - 1
      edge_nodes(1,e+j) = next + j - 1
      edge_nodes(2,e+j) = next + j
    end do

    next = next + segment_length
    e = e + segment_length - 1

  end do
!
!  Handle the holes.
!
  call p00_hole_num ( test, hole_num )

  allocate ( hole(dim_num,hole_num) )

  do hole_index = 1, hole_num
    call p00_hole_point ( test, hole_index, dim_num, &
      hole(1:dim_num,hole_index) )
  end do
!
!  Write the POLY file.
!
  call poly_write ( file_name, node_num, segment, edge_num, &
    edge_nodes, hole_num, hole )

  deallocate ( edge_nodes )
  deallocate ( hole )
  deallocate ( segment )

  return
end
subroutine p00_sample ( test, m, n, seed, point )

!*****************************************************************************80
!
!! P00_SAMPLE samples points from the region in a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_sample ( m, n, seed, point )
  else if ( test == 2 ) then
    call p02_sample ( m, n, seed, point )
  else if ( test == 3 ) then
    call p03_sample ( m, n, seed, point )
  else if ( test == 4 ) then
    call p04_sample ( m, n, seed, point )
  else if ( test == 5 ) then
    call p05_sample ( m, n, seed, point )
  else if ( test == 6 ) then
    call p06_sample ( m, n, seed, point )
  else if ( test == 7 ) then
    call p07_sample ( m, n, seed, point )
  else if ( test == 8 ) then
    call p08_sample ( m, n, seed, point )
  else if ( test == 9 ) then
    call p09_sample ( m, n, seed, point )
  else if ( test == 10 ) then
    call p10_sample ( m, n, seed, point )
  else if ( test == 11 ) then
    call p11_sample ( m, n, seed, point )
  else if ( test == 12 ) then
    call p12_sample ( m, n, seed, point )
  else if ( test == 13 ) then
    call p13_sample ( m, n, seed, point )
  else if ( test == 14 ) then
    call p14_sample ( m, n, seed, point )
  else if ( test == 15 ) then
    call p15_sample ( m, n, seed, point )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_sample_h1 ( test, m, n, h, seed, point )

!*****************************************************************************80
!
!! P00_SAMPLE_H1 samples points from the enlarged region in a problem.
!
!  Discussion:
!
!    The region is enlarged by an amount H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement amount.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) h
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_sample_h1 ( m, n, h, seed, point )
  else if ( test == 2 ) then
    call p02_sample_h1 ( m, n, h, seed, point )
  else if ( test == 3 ) then
    call p03_sample_h1 ( m, n, h, seed, point )
  else if ( test == 4 ) then
    call p04_sample_h1 ( m, n, h, seed, point )
  else if ( test == 5 ) then
   call p05_sample_h1 ( m, n, h, seed, point )
  else if ( test == 6 ) then
    call p06_sample_h1 ( m, n, h, seed, point )
  else if ( test == 7 ) then
    call p07_sample_h1 ( m, n, h, seed, point )
  else if ( test == 8 ) then
    call p08_sample_h1 ( m, n, h, seed, point )
  else if ( test == 9 ) then
    call p09_sample_h1 ( m, n, h, seed, point )
  else if ( test == 10 ) then
    call p10_sample_h1 ( m, n, h, seed, point )
  else if ( test == 11 ) then
    call p11_sample_h1 ( m, n, h, seed, point )
  else if ( test == 12 ) then
    call p12_sample_h1 ( m, n, h, seed, point )
  else if ( test == 13 ) then
    call p13_sample_h1 ( m, n, h, seed, point )
  else if ( test == 14 ) then
    call p14_sample_h1 ( m, n, h, seed, point )
  else if ( test == 15 ) then
    call p15_sample_h1 ( m, n, h, seed, point )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SAMPLE_H1 - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_sdist ( test, m, n, point, sdist )

!*****************************************************************************80
!
!! P00_SDIST returns the signed distance to the region in a problem.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance
!    of each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)
  integer ( kind = 4 ) test

  if ( test == 1 ) then
    call p01_sdist ( m, n, point, sdist )
  else if ( test == 2 ) then
    call p02_sdist ( m, n, point, sdist )
  else if ( test == 3 ) then
    call p03_sdist ( m, n, point, sdist )
  else if ( test == 4 ) then
    call p04_sdist ( m, n, point, sdist )
  else if ( test == 5 ) then
    call p05_sdist ( m, n, point, sdist )
  else if ( test == 6 ) then
    call p06_sdist ( m, n, point, sdist )
  else if ( test == 7 ) then
    call p07_sdist ( m, n, point, sdist )
  else if ( test == 8 ) then
    call p08_sdist ( m, n, point, sdist )
  else if ( test == 9 ) then
    call p09_sdist ( m, n, point, sdist )
  else if ( test == 10 ) then
    call p10_sdist ( m, n, point, sdist )
  else if ( test == 11 ) then
    call p11_sdist ( m, n, point, sdist )
  else if ( test == 12 ) then
    call p12_sdist ( m, n, point, sdist )
  else if ( test == 13 ) then
    call p13_sdist ( m, n, point, sdist )
  else if ( test == 14 ) then
    call p14_sdist ( m, n, point, sdist )
  else if ( test == 15 ) then
    call p15_sdist ( m, n, point, sdist )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SDIST - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_test_num ( test_num )

!*****************************************************************************80
!
!! P00_TEST_NUM returns the number of available tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) TEST_NUM, the number of tests.
!
  implicit none

  integer ( kind = 4 ) test_num

  test_num = 15

  return
end
subroutine p00_title ( test, title )

!*****************************************************************************80
!
!! P00_TITLE returns a title for a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the test problem
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  integer ( kind = 4 ) test
  character ( len = * ) title

  if ( test == 1 ) then
    call p01_title ( title )
  else if ( test == 2 ) then
    call p02_title ( title )
  else if ( test == 3 ) then
    call p03_title ( title )
  else if ( test == 4 ) then
    call p04_title ( title )
  else if ( test == 5 ) then
    call p05_title ( title )
  else if ( test == 6 ) then
    call p06_title ( title )
  else if ( test == 7 ) then
    call p07_title ( title )
  else if ( test == 8 ) then
    call p08_title ( title )
  else if ( test == 9 ) then
    call p09_title ( title )
  else if ( test == 10 ) then
    call p10_title ( title )
  else if ( test == 11 ) then
    call p11_title ( title )
  else if ( test == 12 ) then
    call p12_title ( title )
  else if ( test == 13 ) then
    call p13_title ( title )
  else if ( test == 14 ) then
    call p14_title ( title )
  else if ( test == 15 ) then
    call p15_title ( title )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p01_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P01_BOUNDARY_NEAREST returns a nearest boundary point in problem 01.
!
!  Discussion:
!
!    The given input point need not be inside the region.
!
!    In some cases, more than one boundary point may be "nearest",
!    but only one will be returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( m, n ) :: boundary
  real ( kind = 8 ), dimension ( 2 ) :: center =  (/ &
    0.0D+00, 0.0D+00 /)
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) :: r
  real ( kind = 8 ), parameter :: r1 = 1.0D+00

  do j = 1, n

    r = sqrt ( sum ( ( point(1:m,j) - center(1:m) )**2 ) )

    if ( r == 0.0D+00 ) then
      boundary(1,j) = center(1) + r1
      boundary(2,j) = center(2)
    else
      boundary(1:m,j) = center(1:m) &
        + ( r1 / r ) * ( point(1:m,j) - center(1:m) )
    end if

  end do

  return
end
subroutine p01_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P01_BOUNDARY_PROJECT projects exterior points to the boundary in problem 01.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ), parameter :: center =  &
    (/ 0.0D+00, 0.0D+00 /)
  logical inside(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) :: r
  real ( kind = 8 ), parameter :: r1 = 1.0D+00

  call p01_inside ( m, n, point, inside )

  do j = 1, n

    if ( inside(j) ) then
      cycle
    end if

    r = sqrt ( sum ( ( point(1:m,j) - center(1:m) )**2 ) )

    if ( r == 0.0D+00 ) then
      point(1,j) = center(1) + r1
      point(2,j) = center(2)
    else
      point(1:m,j) = center(1:m) &
        + ( r1 / r ) * ( point(1:m,j) - center(1:m) )
    end if

  end do

  return
end
subroutine p01_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P01_BOUNDARY_SEGMENT returns a boundary segment in problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) segment(m,segment_length)

  if ( segment_index == 1 ) then

    do i = 1, segment_length

      angle = 2.0D+00 * pi &
        * real (              i - 1, kind = 8 ) &
        / real ( segment_length - 1, kind = 8 )

      segment(1,i) = center(1) + r1 * cos ( angle )
      segment(2,i) = center(2) + r1 * sin ( angle )

    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P01_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p01_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P01_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the
!    boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points
!    in the segment.
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P01_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    n = nint ( 2.0D+00 * pi * r1 / h )
    n = max ( n, 5 )
    segment_length = n

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P01_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p01_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P01_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of
!    boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 1

  return
end
subroutine p01_box ( m, lo, hi )

!*****************************************************************************80
!
!! P01_BOX returns a bounding box for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00

  lo(1:m) = (/ center(1) - r1, center(2) - r1 /)
  hi(1:m) = (/ center(1) + r1, center(2) + r1 /)

  return
end
subroutine p01_density ( m, n, point, density )

!*****************************************************************************80
!
!! P01_DENSITY returns the density for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density
!    at each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p01_element_size ( element_size )

!*****************************************************************************80
!
!! P01_ELEMENT_SIZE returns a typical element size for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.2D+00

  return
end
subroutine p01_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P01_FIXED_NUM returns the number of fixed points in problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 0

  return
end
subroutine p01_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P01_FIXED_POINTS returns the fixed points in problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) fixed(m,fixed_num)

  return
end
subroutine p01_header ( )

!*****************************************************************************80
!
!! P01_HEADER prints some information about problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  integer ( kind = 4 ) boundary_segment_num
  real ( kind = 8 ), dimension ( m ) :: center = (/ 0.0D+00, 0.0D+00 /)
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num
  real ( kind = 8 ), parameter :: r1 = 1.0D+00

  call p01_boundary_segment_num ( boundary_segment_num )
  call p01_fixed_num ( fixed_num )
  call p01_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P01:'
  write ( *, '(a)' ) '  Strang and Persson example #1'
  write ( *, '(a)' ) '  The unit circle.'
  write ( *, '(a,g14.6)' ) '  Radius = ', r1
  write ( *, '(a,2g14.6)' ) '  Center = ', center(1:2)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A uniform mesh density is requested.'
  write ( *, '(a)' ) '  Element sizes tried were 0.4, 0.2, 0.1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p01_hole_num ( hole_num )

!*****************************************************************************80
!
!! P01_HOLE_NUM counts the holes in problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 0

  return
end
subroutine p01_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P01_HOLE_POINT returns a point inside a given hole in problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  return
end
subroutine p01_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P01_INSIDE reports if a point is inside the region in problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  logical inside(n)
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00

  inside(1:n) = ( ( point(1,1:n) - center(1) )**2 &
                + ( point(2,1:n) - center(2) )**2 ) <= r1**2

  return
end
subroutine p01_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P01_SAMPLE samples points from the region in problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle(n)
  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) r(n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  integer ( kind = 4 ) seed
!
!  Choose uniform random angles between 0 and 2*Pi.
!
  call r8vec_uniform_01 ( n, seed, angle )

  angle(1:n) = 2.0D+00 * pi * angle(1:n)
!
!  Choose uniform random radii, then take square root.
!
  call r8vec_uniform_01 ( n, seed, r )

  r(1:n) = sqrt ( r(1:n) )
!
!  Construct the uniformly random points in the circle of radius R1
!  centered at CENTER.
!
  point(1,1:n) = center(1) + r1 * r(1:n) * cos ( angle(1:n) )
  point(2,1:n) = center(2) + r1 * r(1:n) * sin ( angle(1:n) )

  return
end
subroutine p01_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P01_SAMPLE_H1 samples points from the enlarged region in problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement amount.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle(n)
  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) h
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) r(n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  integer ( kind = 4 ) seed
!
!  Choose uniform random angles between 0 and 2*Pi.
!
  call r8vec_uniform_01 ( n, seed, angle )

  angle(1:n) = 2.0D+00 * pi * angle(1:n)
!
!  Choose uniform random radii, then take square root.
!
  call r8vec_uniform_01 ( n, seed, r )

  r(1:n) = sqrt ( r(1:n) )
!
!  Construct the uniformly random points in the circle of radius R1 + H
!  centered at CENTER.
!
  point(1,1:n) = center(1) + ( r1 + h ) * r(1:n) * cos ( angle(1:n) )
  point(2,1:n) = center(2) + ( r1 + h ) * r(1:n) * sin ( angle(1:n) )

  return
end
subroutine p01_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P01_SDIST returns the signed distance to the region in problem 01.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ) sdist(n)

  sdist(1:n) = sqrt ( ( point(1,1:n) - center(1) )**2 &
                    + ( point(2,1:n) - center(2) )**2 ) - r1

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns a title for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = '#1: The unit circle.'

  return
end
subroutine p02_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P02_BOUNDARY_NEAREST returns a nearest boundary point in problem 02.
!
!  Discussion:
!
!    The given input point need not be inside the region.
!
!    In some cases, more than one boundary point may be "nearest",
!    but only one will be returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( m, n ) :: boundary
  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00

  do j = 1, n

    if ( point(1,j) == center(1) .and. point(2,j) == center(2) ) then

      boundary(1,j) = r2
      boundary(2,j) = center(2)

    else

      r = sqrt ( ( point(1,j) - center(1) )**2 &
               + ( point(2,j) - center(2) )**2 )

      if ( r1 - r < r - r2 ) then
        boundary(1,j) = center(1) + r1 * ( point(1,j) - center(1) ) / r
        boundary(2,j) = center(2) + r1 * ( point(2,j) - center(2) ) / r
      else
        boundary(1,j) = center(1) + r2 * ( point(1,j) - center(1) ) / r
        boundary(2,j) = center(2) + r2 * ( point(2,j) - center(2) ) / r
      end if

    end if

  end do

  return
end
subroutine p02_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P02_BOUNDARY_PROJECT projects exterior points to the boundary in problem 02.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  logical inside(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00

  call p02_inside ( m, n, point, inside )

  do j = 1, n

    if ( inside(j) ) then
      cycle
    end if

    if ( point(1,j) == center(1) .and. point(2,j) == center(2) ) then

      point(1,j) = r2
      point(2,j) = center(2)

    else

      r = sqrt ( ( point(1,j) - center(1) )**2 &
               + ( point(2,j) - center(2) )**2 )

      if ( r1 - r < r - r2 ) then
        point(1,j) = center(1) + r1 * ( point(1,j) - center(1) ) / r
        point(2,j) = center(2) + r1 * ( point(2,j) - center(2) ) / r
      else
        point(1,j) = center(1) + r2 * ( point(1,j) - center(1) ) / r
        point(2,j) = center(2) + r2 * ( point(2,j) - center(2) ) / r
      end if

    end if

  end do

  return
end
subroutine p02_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P02_BOUNDARY_SEGMENT returns a boundary segment in problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) segment(m,segment_length)

  if ( segment_index == 1 ) then

    do j = 1, segment_length
      angle = 2.0D+00 * pi &
        * real (              j - 1, kind = 8 ) &
        / real ( segment_length - 1, kind = 8 )
      segment(1,j) = center(1) + r1 * cos ( angle )
      segment(2,j) = center(2) + r1 * sin ( angle )
    end do

  else if ( segment_index == 2 ) then

    do j = 1, segment_length
      angle = 2.0D+00 * pi &
        * real ( segment_length - j, kind = 8 ) &
        / real ( segment_length - 1, kind = 8 )
      segment(1,j) = center(1) + r2 * cos ( angle )
      segment(2,j) = center(2) + r2 * sin ( angle )
    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P02_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p02_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P02_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the
!    boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points
!    in the segment.
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P02_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    n = nint ( 2.0D+00 * pi * r1 / h )
    n = max ( n, 5 )
    segment_length = n

  else if ( segment_index == 2 ) then

    n = nint ( 2.0D+00 * pi * r2 / h )
    n = max ( n, 5 )
    segment_length = n

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P02_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop
  end if

  return
end
subroutine p02_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P02_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 2

  return
end
subroutine p02_box ( m, lo, hi )

!*****************************************************************************80
!
!! P02_BOX returns a bounding box for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00

  lo(1:m) = (/ center(1) - r1, center(2) - r1 /)
  hi(1:m) = (/ center(1) + r1, center(2) + r1 /)

  return
end
subroutine p02_density ( m, n, point, density )

!*****************************************************************************80
!
!! P02_DENSITY returns the density for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p02_element_size ( element_size )

!*****************************************************************************80
!
!! P02_ELEMENT_SIZE returns a typical element size for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.1D+00

  return
end
subroutine p02_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P02_FIXED_NUM returns the number of fixed points in problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 0

  return
end
subroutine p02_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P02_FIXED_POINTS returns the fixed points in problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) fixed(m,fixed_num)

  return
end
subroutine p02_header ( )

!*****************************************************************************80
!
!! P02_HEADER prints some information about problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  integer ( kind = 4 ) boundary_segment_num
  real ( kind = 8 ), dimension ( m ) :: center = (/ 0.0D+00, 0.0D+00 /)
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00

  call p02_boundary_segment_num ( boundary_segment_num )
  call p02_fixed_num ( fixed_num )
  call p02_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P02:'
  write ( *, '(a)' ) '  Strang and Persson example #2'
  write ( *, '(a)' ) '  The unit circle, with a concentric hole.'
  write ( *, '(a,g14.6)' ) '  Radius1 = ', r1
  write ( *, '(a,g14.6)' ) '  Radius2 = ', r2
  write ( *, '(a,2g14.6)' ) '  Center = ', center(1:2)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A uniform mesh density is requested.'
  write ( *, '(a)' ) '  Element sizes tried were 0.4, 0.2, 0.1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p02_hole_num ( hole_num )

!*****************************************************************************80
!
!! P02_HOLE_NUM counts the holes in problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 1

  return
end
subroutine p02_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P02_HOLE_POINT returns a point inside a given hole in problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  hole_point(1:2) = center(1:2)

  return
end
subroutine p02_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P02_INSIDE reports if a point is inside the region in problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  logical inside(n)
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00

  inside(1:n) =                                                &
    (                                                          &
                ( point(1,1:n) - center(1) )**2                &
              + ( point(2,1:n) - center(2) )**2 <= r1**2       &
    )                                                          &
    .and.                                                      &
    (                                                          &
      r2**2 <= ( point(1,1:n) - center(1) )**2                 &
             + ( point(2,1:n) - center(2) )**2                 &
    )

  return
end
subroutine p02_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P02_SAMPLE samples points from the region in problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle(n)
  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) r(n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  integer ( kind = 4 ) seed
!
!  Choose uniform random angles between 0 and 2*Pi.
!
  call r8vec_uniform_01 ( n, seed, angle )

  angle(1:n) = 2.0D+00 * pi * angle(1:n)
!
!  Choose uniform random radii between R2**2 and R1**2, then take square root.
!
  call r8vec_uniform_01 ( n, seed, r )

  r(1:n) = r2**2 + ( r1**2 - r2**2 ) * r(1:n)
  r(1:n) = sqrt ( r(1:n) )
!
!  Construct the uniformly random points in the circle of radius 1
!  centered at 0.
!
  point(1,1:n) = center(1) + r(1:n) * cos ( angle(1:n) )
  point(2,1:n) = center(2) + r(1:n) * sin ( angle(1:n) )

  return
end
subroutine p02_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P02_SAMPLE_H1 samples points from the enlarged region in problem 02.
!
!  Discussion:
!
!    The inner circular hole can't be enlarged more than R2.
!    Normally, of course, H is intended to be small anyway.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement amount.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle(n)
  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) h
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) r(n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  integer ( kind = 4 ) seed
!
!  Choose uniform random angles between 0 and 2*Pi.
!
  call r8vec_uniform_01 ( n, seed, angle )

  angle(1:n) = 2.0D+00 * pi * angle(1:n)
!
!  Choose uniform random radii between (R1+H)**2 and (R2-H)**2,
!  then take square root.
!
  call r8vec_uniform_01 ( n, seed, r )

  r(1:n) = ( r2 - h )**2 &
    + ( ( max ( r1 + h, 0.0D+00 ) )**2 - ( r2 - h )**2 ) * r(1:n)

  r(1:n) = sqrt ( r(1:n) )
!
!  Construct the uniformly random points in the circle of radius 1
!  centered at 0.
!
  point(1,1:n) = center(1) + r(1:n) * cos ( angle(1:n) )
  point(2,1:n) = center(2) + r(1:n) * sin ( angle(1:n) )

  return
end
subroutine p02_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P02_SDIST returns the signed distance to the region in problem 02.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ) :: center = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  real ( kind = 8 ) sdist(n)

  sdist(1:n) = max                                            &
    (                                                         &
          sqrt ( ( point(1,1:n) - center(1) )**2              &
               + ( point(2,1:n) - center(2) )**2 ) - r1       &
    ,                                                         &
      - ( sqrt ( ( point(1,1:n) - center(1) )**2              &
               + ( point(2,1:n) - center(2) )**2 ) - r2    )  &
    )

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns a title for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = '#2: The circle with concentric circular hole.'

  return
end
subroutine p03_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P03_BOUNDARY_NEAREST returns a nearest boundary point in problem 03.
!
!  Discussion:
!
!    The given input point need not be inside the region.
!
!    In some cases, more than one boundary point may be "nearest",
!    but only one will be returned.
!
!    31 August 2005: Thanks to Hua Fei for pointing out that a previous
!    version of this routine gave inaccurate results for points that were
!    significantly far from the box.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( m, n ) :: boundary
  real ( kind = 8 ), dimension ( 2 ), parameter :: center =  &
    (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) dist_box
  real ( kind = 8 ) dist_circle
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) :: r
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: x1 = -1.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ) y
  real ( kind = 8 ), parameter :: y1 = -1.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00

  do j = 1, n

    x = point(1,j)
    y = point(2,j)
!
!  Special case: a point at the center of the box.
!  The closest point is ANY point on the circle.
!
    if ( point(1,j) == center(1) .and. point(2,j) == center(2) ) then
      boundary(1,j) = center(1) + r2
      boundary(2,j) = center(2)
      cycle
    end if
!
!  Is the point to the left of the box?
!
    if ( x <= x1 ) then

      boundary(1,j) = x1

      if ( y <= y1 ) then
        boundary(2,j) = y1
      else if ( y <= y2 ) then
        boundary(2,j) = y
      else if ( y2 <= y ) then
        boundary(2,j) = y2
      end if
!
!  To the right of the box?
!
    else if ( x2 <= x ) then

      boundary(1,j) = x2

      if ( y <= y1 ) then
        boundary(2,j) = y1
      else if ( y <= y2 ) then
        boundary(2,j) = y
      else if ( y2 <= y ) then
        boundary(2,j) = y2
      end if
!
!  Below the middle of the box?
!
    else if ( y <= y1 ) then

      boundary(1,j) = x
      boundary(2,j) = y1
!
!  Above the middle of the box?
!
    else if ( y2 <= y ) then

      boundary(1,j) = x
      boundary(2,j) = y2
!
!  Inside the box.
!  Figure out which side is closest by drawing the diagonal lines.
!
    else
!
!  Y is small.
!
      if (     y <= x .and.  y <= -x ) then
        boundary(1,j) = x
        boundary(2,j) = y1
!
!  X is big.
!
      else if ( y <= x .and. -y <=  x ) then
        boundary(1,j) = x2
        boundary(2,j) = y
!
!  Y is big.
!
      else if ( x <= y .and. -x <=  y ) then
        boundary(1,j) = x
        boundary(2,j) = y2
!
!  X is small.
!
      else if ( x <= y .and.  x <= -y ) then
        boundary(1,j) = x1
        boundary(2,j) = y
      end if
!
!  For points inside the box, the boundary of the circle might be closer than
!  the boundary of the box.
!
      r = sqrt ( sum ( ( point(1:m,j) - center(1:m) )**2 ) )

      dist_circle = abs ( r - r2 )
      dist_box = sqrt ( sum ( ( point(1:m,j) - boundary(1:m,j) )**2 ) )

      if ( dist_circle <  dist_box ) then
        boundary(1:m,j) = center(1:m) + r2 / r * ( point(1:m,j) - center(1:m) )
      end if

    end if

  end do

  return
end
subroutine p03_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P03_BOUNDARY_PROJECT projects exterior points to the boundary in problem 03.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!    01 September 2005: Thanks to Hua Fei for pointing out that a previous
!    version of this routine gave inaccurate results for points that were
!    significantly far from the box.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical inside(n)
  real ( kind = 8 ), dimension ( 2 ), parameter :: center =  &
    (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_min
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) :: r
  real ( kind = 8 ) temp(2)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: x1 = -1.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ) y
  real ( kind = 8 ), parameter :: y1 = -1.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00

  do j = 1, n

    x = point(1,j)
    y = point(2,j)
!
!  If the point is INSIDE the box and OUTSIDE the circle,
!  we don't do anything.
!
    if ( x1 <= x .and. x <= x2 .and. y1 <= y .and. y <= y2 .and. &
         r2**2 <= ( x - center(1) )**2 + ( y - center(2) )**2 ) then
      cycle
    end if
!
!  Special case: a point at the center of the box.
!  The closest point is ANY point on the circle.
!
    if ( x == center(1) .and. y == center(2) ) then
      point(1,j) = center(1) + r2
      point(2,j) = center(2)
      cycle
    end if
!
!  Is the point to the left of the box?
!
    if ( x <= x1 ) then

      point(1,j) = x1

      if ( y <= y1 ) then
        point(2,j) = y1
      else if ( y <= y2 ) then
        point(2,j) = y
      else if ( y2 <= y ) then
        point(2,j) = y2
      end if
!
!  To the right of the box?
!
    else if ( x2 <= x ) then

      point(1,j) = x2

      if ( y <= y1 ) then
        point(2,j) = y1
      else if ( y <= y2 ) then
        point(2,j) = y
      else if ( y2 <= y ) then
        point(2,j) = y2
      end if
!
!  Below the middle of the box?
!
    else if ( y <= y1 ) then

      point(1,j) = x
      point(2,j) = y1
!
!  Above the middle of the box?
!
    else if ( y2 <= y ) then

      point(1,j) = x
      point(2,j) = y2
!
!  Last chance: Must be inside the circle.
!
    else

      r = sqrt ( sum ( ( point(1:m,j) - center(1:m) )**2 ) )

      point(1:m,j) = center(1:m) + r2 / r * ( point(1:m,j) - center(1:m) )

    end if

  end do

  return
end
subroutine p03_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P03_BOUNDARY_SEGMENT returns a boundary segment in problem 03.
!
!  Discussion:
!
!    For boundary segment #1, the value of SEGMENT_LENGTH should be
!    at least 5.  Values of 4*N+1 will result in an "even" mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( 2 ), parameter :: center =  &
    (/ 0.0D+00, 0.0D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) s(m,4)
  real ( kind = 8 ) segment(m,segment_length)

  if ( segment_index == 1 ) then

    n1 = nint ( real ( segment_length - 1, kind = 8 ) &
              / real ( 4, kind = 8 ) )
    n2 = nint ( real ( 2 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 4, kind = 8 ) ) - n1
    n3 = nint ( real ( 3 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 4, kind = 8 ) ) - n1 - n2
    n4 = segment_length - 1 - n1 - n2 - n3

    s(1:2,1) = (/ -1.0D+00, -1.0D+00 /)
    s(1:2,2) = (/  1.0D+00, -1.0D+00 /)
    s(1:2,3) = (/  1.0D+00,  1.0D+00 /)
    s(1:2,4) = (/ -1.0D+00,  1.0D+00 /)

    j = 0

    do i = 1, n1
      j = j + 1
      segment(1:2,j) = ( real ( n1 - i + 1, kind = 8 ) * s(1:2,1)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,2) ) &
                       / real ( n1,         kind = 8 )
    end do

    do i = 1, n2
      j = j + 1
      segment(1:2,j) = ( real ( n2 - i + 1, kind = 8 ) * s(1:2,2)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,3) ) &
                       / real ( n2,         kind = 8 )
    end do

    do i = 1, n3
      j = j + 1
      segment(1:2,j) = ( real ( n3 - i + 1, kind = 8 ) * s(1:2,3)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,4) ) &
                       / real ( n3,         kind = 8 )
    end do

    do i = 1, n4
      j = j + 1
      segment(1:2,j) = ( real ( n4 - i + 1, kind = 8 ) * s(1:2,4)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,1) ) &
                       / real ( n4,         kind = 8 )
    end do

    j = j + 1
    segment(1:2,j) = s(1:2,1)

  else if ( segment_index == 2 ) then

    do i = 1, segment_length
      angle = 2.0D+00 * pi &
        * real ( segment_length - i, kind = 8 ) &
        / real ( segment_length - 1, kind = 8 )
      segment(1,i) = center(1) + r2 * cos ( angle )
      segment(2,i) = center(2) + r2 * sin ( angle )
    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P03_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p03_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P03_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the
!    boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in
!    the segment.
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P03_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    n = nint ( 8.0D+00 / h )
    n = max ( n, 5 )
    segment_length = n + mod ( 4 - mod ( n - 1, 4 ), 4 )

  else if ( segment_index == 2 ) then

    n = nint ( 2.0D+00 * pi * r2 / h )
    n = max ( n, 5 )
    segment_length = n

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P03_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop
  end if

  return
end
subroutine p03_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P03_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of
!    boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 2

  return
end
subroutine p03_box ( m, lo, hi )

!*****************************************************************************80
!
!! P03_BOX returns a bounding box for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)

  lo(1:m) = (/ -1.0D+00, -1.0D+00 /)
  hi(1:m) = (/ +1.0D+00, +1.0D+00 /)

  return
end
subroutine p03_density ( m, n, point, density )

!*****************************************************************************80
!
!! P03_DENSITY returns the density for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p03_element_size ( element_size )

!*****************************************************************************80
!
!! P03_ELEMENT_SIZE returns a typical element size for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.15D+00

  return
end
subroutine p03_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P03_FIXED_NUM returns the number of fixed points in problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 4

  return
end
subroutine p03_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P03_FIXED_POINTS returns the fixed points in problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) fixed(m,fixed_num)

  fixed(1:2,1) = (/ -1.0D+00, -1.0D+00 /)
  fixed(1:2,2) = (/  1.0D+00, -1.0D+00 /)
  fixed(1:2,3) = (/  1.0D+00,  1.0D+00 /)
  fixed(1:2,4) = (/ -1.0D+00,  1.0D+00 /)

  return
end
subroutine p03_header ( )

!*****************************************************************************80
!
!! P03_HEADER prints some information about problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num

  call p03_boundary_segment_num ( boundary_segment_num )
  call p03_fixed_num ( fixed_num )
  call p03_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P03:'
  write ( *, '(a)' ) '  Strang and Persson example #3'
  write ( *, '(a)' ) '  The unit square, with a hole.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The hole is a concentric circle of radius 0.4.'
  write ( *, '(a)' ) '  A uniform mesh density is requested.'
  write ( *, '(a)' ) '  Element sizes tried were 0.4, 0.2, 0.1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p03_hole_num ( hole_num )

!*****************************************************************************80
!
!! P03_HOLE_NUM counts the holes in problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 1

  return
end
subroutine p03_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P03_HOLE_POINT returns a point inside a given hole in problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  hole_point(1:2) = (/ 0.0D+00, 0.0D+00 /)

  return
end
subroutine p03_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P03_INSIDE reports if a point is inside the region in problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ), parameter :: center =  &
    (/ 0.0D+00, 0.0D+00 /)
  logical inside(n)
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  real ( kind = 8 ), parameter :: x1 = -1.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ), parameter :: y1 = -1.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00

  inside(1:n) =                           &
      x1 <=      point(1,1:n)       .and. &
                 point(1,1:n) <= x2 .and. &
      y1 <=      point(2,1:n)       .and. &
                 point(2,1:n) <= y2 .and. &
      r2**2 <= ( point(1,1:n) - center(1) )**2 &
             + ( point(2,1:n) - center(2) )**2

  return
end
subroutine p03_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P03_SAMPLE samples points from the region in problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ), parameter :: center =  &
    (/ 0.0D+00, 0.0D+00 /)
  integer ( kind = 4 ) have
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ), parameter :: x1 = -1.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ), parameter :: y1 = -1.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00

  have = 0
!
!  We are going to generate batches of sample points.
!
  sample_num = min ( 1000, n )

  allocate ( sample(1:m,1:sample_num) )

  do
!
!  Generate a batch of points in [0,1]x[0,1].
!
    call r8mat_uniform_01 ( m, sample_num, seed, sample )
!
!  Remap the points to the box [X1,X2] x [Y1,Y2].
!
    sample(1,1:sample_num) = x1 + sample(1,1:sample_num) * ( x2 - x1 )
    sample(2,1:sample_num) = y1 + sample(2,1:sample_num) * ( y2 - y1 )
!
!  Accept those points which are NOT in the hole.
!
    do j = 1, sample_num

      if ( r2**2 < ( sample(1,j) - center(1) )**2 &
                 + ( sample(2,j) - center(2) )**2 ) then

        have = have + 1
        point(1:m,have) = sample(1:m,j)

        if ( have == n ) then
          return
        end if

      end if

    end do

  end do

  deallocate ( sample )

  return
end
subroutine p03_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P03_SAMPLE samples points from the enlarged region in problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement amount.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ), parameter :: center =  &
    (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) h
  integer ( kind = 4 ) have
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ), parameter :: x1 = -1.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ), parameter :: y1 = -1.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00

  have = 0
!
!  We are going to generate batches of sample points.
!
  sample_num = min ( 1000, n )

  allocate ( sample(1:m,1:sample_num) )

  do
!
!  Generate a batch of points in [0,1]x[0,1].
!
    call r8mat_uniform_01 ( m, sample_num, seed, sample )
!
!  Remap the points to the box [X1,X2] x [Y1,Y2].
!
    sample(1,1:sample_num) = ( x1 - h ) &
      + sample(1,1:sample_num) * ( x2 - x1 + 2.0D+00 * h )
    sample(2,1:sample_num) = ( y1 - h ) &
      + sample(2,1:sample_num) * ( y2 - y1 + 2.0D+00 * h )
!
!  Accept those points which are NOT in the hole.
!
    if ( 0.0D+00 < r2 - h ) then

      do j = 1, sample_num

        if ( ( r2 - h )**2 < ( sample(1,j) - center(1) )**2 &
                           + ( sample(2,j) - center(2) )**2 ) then

          have = have + 1
          point(1:m,have) = sample(1:m,j)

          if ( have == n ) then
            return
          end if

        end if

      end do

    end if

  end do

  deallocate ( sample )

  return
end
subroutine p03_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P03_SDIST returns the signed distance to the region in problem 03.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!    The signed distance formula used for the rectangle is incorrect
!    for some points outside the rectangle; for the Persson and Strang
!    approach, this defect can be dealt with by fixing mesh nodes at
!    the corners in advance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ), parameter :: center =  &
    (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.4D+00
  real ( kind = 8 ) sdist(n)
  real ( kind = 8 ), parameter :: x1 = -1.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ), parameter :: xc =  0.0D+00
  real ( kind = 8 ), parameter :: y1 = -1.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00
  real ( kind = 8 ), parameter :: yc = 0.0D+00

  sdist(1:n) = max                                        &
    (                                                     &
     - min ( x2 - point(1,1:n),                           &
       min (      point(1,1:n) - x1,                      &
       min (      point(2,1:n) - y1,                      &
             y2 - point(2,1:n) ) ) )                      &
    ,                                                     &
      - ( sqrt ( ( point(1,1:n) - center(1) )**2          &
               + ( point(2,1:n) - center(2) )**2 ) - r2 ) &
    )

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns a title for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = '#3: The unit square with circular hole.'

  return
end
subroutine p04_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P04_BOUNDARY_NEAREST returns a nearest boundary point in problem 04.
!
!  Discussion:
!
!    The given input point need not be inside the region.
!
!    In some cases, more than one boundary point may be "nearest",
!    but only one will be returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( m, n ) :: boundary
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) pn(m)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  real ( kind = 8 ) t
  logical, save :: v_set = .false.
  real ( kind = 8 ), save, dimension(2,6) :: v1
  real ( kind = 8 ), save, dimension(2,6) :: v2

  if ( .not. v_set ) then

    do i = 0, 5
      angle = real ( 2 * ( i - 1 ), kind = 8 ) * pi / 6.0D+00
      v1(1,i+1) = r1 * cos ( angle )
      v1(2,i+1) = r1 * sin ( angle )
    end do

    do i = 6, 1, -1
      angle = real ( 2 * i - 1, kind = 8 ) * pi / 6.0D+00
      v2(1,7-i) = r2 * cos ( angle )
      v2(2,7-i) = r2 * sin ( angle )
    end do

    v_set = .true.

  end if

  do j = 1, n

    dist_min = huge ( dist_min )
    boundary(1:m,j) = 0.0D+00
!
!  Examine points on the outer hexagon.
!
    k = 6

    do i = 1, 6

      call segment_point_near_2d ( v1(1:2,k), v1(1:2,i), point(1:2,j), &
        pn(1:2), dist, t )

      if ( dist < dist_min ) then
        dist_min = dist
        boundary(1:2,j) = pn(1:2)
      end if

      k = i

    end do
!
!  Examine points on the inner hexagon.
!
    k = 6

    do i = 1, 6

      call segment_point_near_2d ( v2(1:2,k), v2(1:2,i), point(1:2,j), &
        pn(1:2), dist, t )

      if ( dist < dist_min ) then
        dist_min = dist
        boundary(1:2,j) = pn(1:2)
      end if

      k = i

    end do

  end do

  return
end
subroutine p04_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P04_BOUNDARY_PROJECT projects exterior points to the boundary in problem 04.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  logical inside(n)
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) pn(m)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) temp(m)
  logical, save :: v_set = .false.
  real ( kind = 8 ), save, dimension(2,6) :: v1
  real ( kind = 8 ), save, dimension(2,6) :: v2

  if ( .not. v_set ) then

    do i = 0, 5
      angle = real ( 2 * ( i - 1 ), kind = 8 ) * pi / 6.0D+00
      v1(1,i+1) = r1 * cos ( angle )
      v1(2,i+1) = r1 * sin ( angle )
    end do

    do i = 6, 1, -1
      angle = real ( 2 * i - 1, kind = 8 ) * pi / 6.0D+00
      v2(1,7-i) = r2 * cos ( angle )
      v2(2,7-i) = r2 * sin ( angle )
    end do

    v_set = .true.

  end if

  call p04_inside ( m, n, point, inside )

  do j = 1, n

    if ( inside(j) ) then
      cycle
    end if

    dist_min = huge ( dist_min )
    temp(1:m) = 0.0D+00
!
!  Examine points on the outer hexagon.
!
    k = 6

    do i = 1, 6

      call segment_point_near_2d ( v1(1:2,k), v1(1:2,i), point(1:2,j), &
        pn(1:2), dist, t )

      if ( dist < dist_min ) then
        dist_min = dist
        temp(1:m) = pn(1:2)
      end if

      k = i

    end do
!
!  Examine points on the inner hexagon.
!
    k = 6

    do i = 1, 6

      call segment_point_near_2d ( v2(1:2,k), v2(1:2,i), point(1:2,j), &
        pn(1:2), dist, t )

      if ( dist < dist_min ) then
        dist_min = dist
        temp(1:m) = pn(1:2)
      end if

      k = i

    end do

    point(1:m,j) = temp(1:m)

  end do

  return
end
subroutine p04_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P04_BOUNDARY_SEGMENT returns a boundary segment in problem 04.
!
!  Discussion:
!
!    SEGMENT_LENGTH should be no less than 7, and good values will
!    have the form 6*N+1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in
!    the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) s(m,6)
  real ( kind = 8 ) segment(m,segment_length)

  if ( segment_index == 1 ) then

    n1 = nint ( real ( segment_length - 1, kind = 8 ) &
              / real ( 6, kind = 8 ) )
    n2 = nint ( real ( 2 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1
    n3 = nint ( real ( 3 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1 - n2
    n4 = nint ( real ( 4 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1 - n2 - n3
    n5 = nint ( real ( 5 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1 - n2 - n3 - n4
    n6 = segment_length - 1 - n1 - n2 - n3 - n4 - n5

    do j = 1, 6
      angle = real ( 2 * ( j - 1 ), kind = 8 ) * pi / 6.0D+00
      s(1,j) = r1 * cos ( angle )
      s(2,j) = r1 * sin ( angle )
    end do

    j = 0

    do i = 1, n1
      j = j + 1
      segment(1:2,j) = ( real ( n1 - i + 1, kind = 8 ) * s(1:2,1)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,2) ) &
                       / real ( n1,         kind = 8 )
    end do

    do i = 1, n2
      j = j + 1
      segment(1:2,j) = ( real ( n2 - i + 1, kind = 8 ) * s(1:2,2)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,3) ) &
                       / real ( n2,         kind = 8 )
    end do

    do i = 1, n3
      j = j + 1
      segment(1:2,j) = ( real ( n3 - i + 1, kind = 8 ) * s(1:2,3)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,4) ) &
                       / real ( n3,         kind = 8 )
    end do

    do i = 1, n4
      j = j + 1
      segment(1:2,j) = ( real ( n4 - i + 1, kind = 8 ) * s(1:2,4)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,5) ) &
                       / real ( n4,         kind = 8 )
    end do

    do i = 1, n5
      j = j + 1
      segment(1:2,j) = ( real ( n5 - i + 1, kind = 8 ) * s(1:2,5)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,6) ) &
                       / real ( n5,         kind = 8 )
    end do

    do i = 1, n6
      j = j + 1
      segment(1:2,j) = ( real ( n6 - i + 1, kind = 8 ) * s(1:2,6)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,1) ) &
                       / real ( n6,         kind = 8 )
    end do

    j = j + 1
    segment(1:2,j) = s(1:2,1)

  else if ( segment_index == 2 ) then

    n1 = nint ( real ( segment_length - 1, kind = 8 ) &
              / real ( 6, kind = 8 ) )
    n2 = nint ( real ( 2 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1
    n3 = nint ( real ( 3 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1 - n2
    n4 = nint ( real ( 4 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1 - n2 - n3
    n5 = nint ( real ( 5 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1 - n2 - n3 - n4
    n6 = segment_length - 1 - n1 - n2 - n3 - n4 - n5

    do j = 1, 6
      angle = real ( 2 * ( 6 - j ) + 1, kind = 8 ) * pi / 6.0D+00
      s(1,j) = r2 * cos ( angle )
      s(2,j) = r2 * sin ( angle )
    end do

    j = 0

    do i = 1, n1
      j = j + 1
      segment(1:2,j) = ( real ( n1 - i + 1, kind = 8 ) * s(1:2,1)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,2) ) &
                       / real ( n1,         kind = 8 )
    end do

    do i = 1, n2
      j = j + 1
      segment(1:2,j) = ( real ( n2 - i + 1, kind = 8 ) * s(1:2,2)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,3) ) &
                       / real ( n2,         kind = 8 )
    end do

    do i = 1, n3
      j = j + 1
      segment(1:2,j) = ( real ( n3 - i + 1, kind = 8 ) * s(1:2,3)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,4) ) &
                       / real ( n3,         kind = 8 )
    end do

    do i = 1, n4
      j = j + 1
      segment(1:2,j) = ( real ( n4 - i + 1, kind = 8 ) * s(1:2,4)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,5) ) &
                       / real ( n4,         kind = 8 )
    end do

    do i = 1, n5
      j = j + 1
      segment(1:2,j) = ( real ( n5 - i + 1, kind = 8 ) * s(1:2,5)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,6) ) &
                       / real ( n5,         kind = 8 )
    end do

    do i = 1, n6
      j = j + 1
      segment(1:2,j) = ( real ( n6 - i + 1, kind = 8 ) * s(1:2,6)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,1) ) &
                       / real ( n6,         kind = 8 )
    end do

    j = j + 1
    segment(1:2,j) = s(1:2,1)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P04_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p04_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P04_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the
!    boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in
!    the segment.
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P04_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    n = nint ( 6.0D+00 * r1 / h )
    n = max ( n, 7 )
    segment_length = n + mod ( 6 - mod ( n - 1, 6 ), 6 )

  else if ( segment_index == 2 ) then

    n = nint ( 6.0D+00 * r2 / h )
    n = max ( n, 7 )
    segment_length = n + mod ( 6 - mod ( n - 1, 6 ), 6 )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P04_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p04_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P04_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 2

  return
end
subroutine p04_box ( m, lo, hi )

!*****************************************************************************80
!
!! P04_BOX returns a bounding box for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00

  lo(1:m) = (/ -r1, -r1 * sqrt ( 3.0D+00 ) / 2.0D+00 /)
  hi(1:m) = (/ +r1, +r1 * sqrt ( 3.0D+00 ) / 2.0D+00 /)

  return
end
subroutine p04_density ( m, n, point, density )

!*****************************************************************************80
!
!! P04_DENSITY returns the density for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p04_element_size ( element_size )

!*****************************************************************************80
!
!! P04_ELEMENT_SIZE returns a typical element size for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.1D+00

  return
end
subroutine p04_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P04_FIXED_NUM returns the number of fixed points in problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 12

  return
end
subroutine p04_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P04_FIXED_POINTS returns the fixed points in problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  real ( kind = 8 ) fixed(m,fixed_num)

  j = 0

  do i = 0, 5
    angle = real ( 2 * ( i - 1 ), kind = 8 ) * pi / 6.0D+00
    j = j + 1
    fixed(1,j) = r1 * cos ( angle )
    fixed(2,j) = r1 * sin ( angle )
  end do

  do i = 6, 1, -1
    angle = real ( 2 * i - 1, kind = 8 ) * pi / 6.0D+00
    j = j + 1
    fixed(1,j) = r2 * cos ( angle )
    fixed(2,j) = r2 * sin ( angle )
  end do

  return
end
subroutine p04_header ( )

!*****************************************************************************80
!
!! P04_HEADER prints some information about problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00

  call p04_boundary_segment_num ( boundary_segment_num )
  call p04_fixed_num ( fixed_num )
  call p04_hole_num ( hole_num )

  write ( *, '(a)' )       ' '
  write ( *, '(a)' )       'P04:'
  write ( *, '(a)' )       '  Strang and Persson example #4'
  write ( *, '(a)' )       '  The hexagon with hexagonal hole.'
  write ( *, '(a)' )       ' '
  write ( *, '(a,g14.6)' ) '  Radius of outer hexagon R1 = ', r1
  write ( *, '(a,g14.6)' ) '  Radius of outer hexagon R2 = ', r2
  write ( *, '(a)' )       ' '
  write ( *, '(a)' )       '  A uniform mesh density is requested.'
  write ( *, '(a)' )       '  Element sizes tried were ?'
  write ( *, '(a)' )       ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p04_hole_num ( hole_num )

!*****************************************************************************80
!
!! P04_HOLE_NUM counts the holes in problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 1

  return
end
subroutine p04_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P04_HOLE_POINT returns a point inside a given hole in problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  hole_point(1:2) = (/ 0.0D+00, 0.0D+00 /)

  return
end
subroutine p04_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P04_INSIDE reports if a point is inside the region in problem 04.
!
!  Discussion:
!
!    Our test asks if a point is inside the big hexagon and not inside
!    the smaller one.  For this test to work, we need to list the
!    vertices of the smaller hexagon in the same counter clockwise order
!    used for the big hexagon.  In other routines, we list the vertices
!    of the inner hexagon in clockwise order, since we think of it
!    solely as the boundary of the hexagonal annulus.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  logical hexagon_contains_point_2d
  integer ( kind = 4 ) i
  logical inside(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  logical, save :: v_set = .false.
  real ( kind = 8 ), save, dimension ( 2, 6 ) :: v1
  real ( kind = 8 ), save, dimension ( 2, 6 ) :: v2
!
!  Drat bast it!  Here, we need to list the V2 vertices in
!  the same order as V1...
!
  if ( .not. v_set ) then

    do i = 0, 5
      angle = real ( 2 * ( i - 1 ), kind = 8 ) * pi / 6.0D+00
      v1(1,i+1) = r1 * cos ( angle )
      v1(2,i+1) = r1 * sin ( angle )
    end do

    do i = 1, 6
      angle = real ( 2 * i - 1, kind = 8 ) * pi / 6.0D+00
      v2(1,i) = r2 * cos ( angle )
      v2(2,i) = r2 * sin ( angle )
    end do

    v_set = .true.

  end if

  do i = 1, n

    inside(i) = hexagon_contains_point_2d ( v1, point(1:2,i) ) &
      .and. ( .not. &
                hexagon_contains_point_2d ( v2, point(1:2,i) ) )

  end do

  return
end
subroutine p04_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P04_SAMPLE samples points from the region in problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) have
  logical, allocatable, dimension ( : ) :: inside
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  have = 0
!
!  We are going to generate batches of sample points.
!
  sample_num = min ( 1000, n )

  allocate ( inside(1:sample_num) )
  allocate ( sample(1:m,1:sample_num) )

  x1 = -r1
  x2 = +r1
  y1 = -r1
  y2 = +r1

  do
!
!  Generate a batch of points in [0,1]x[0,1].
!
    call r8mat_uniform_01 ( m, sample_num, seed, sample )
!
!  Remap the points to [x1,x2] x [y1,y2].
!
    sample(1,1:sample_num) = x1 + sample(1,1:sample_num) * ( x2 - x1 )
    sample(2,1:sample_num) = x1 + sample(2,1:sample_num) * ( y2 - y1 )
!
!  Accept each point that is inside the region.
!
    call p04_inside ( m, sample_num, sample, inside )

    do j = 1, sample_num

      if ( inside(j) ) then

        have = have + 1
        point(1:m,have) = sample(1:m,j)

        if ( have == n ) then
          return
        end if

      end if

    end do

  end do

  deallocate ( inside )
  deallocate ( sample )

  return
end
subroutine p04_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P04_SAMPLE_H1 samples points from the enlarged region in problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement value.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  real ( kind = 8 ) h
  integer ( kind = 4 ) have
  logical hexagon_contains_point_2d
  integer ( kind = 4 ) i
  logical inside
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  logical, save :: v_set = .false.
  real ( kind = 8 ), save, dimension ( 2, 6 ) :: v1
  real ( kind = 8 ), save, dimension ( 2, 6 ) :: v2
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  if ( .not. v_set ) then

    do i = 0, 5
      angle = real ( 2 * ( i - 1 ), kind = 8 ) * pi / 6.0D+00
      v1(1,i+1) = ( r1 + h ) * cos ( angle )
      v1(2,i+1) = ( r1 + h ) * sin ( angle )
    end do

    do i = 1, 6
      angle = real ( 2 * i - 1, kind = 8 ) * pi / 6.0D+00
      v2(1,i) = max ( ( r2 - h ), 0.0D+00 ) * cos ( angle )
      v2(2,i) = max ( ( r2 - h ), 0.0D+00 ) * sin ( angle )
    end do

    v_set = .true.

  end if

  have = 0
!
!  We are going to generate batches of sample points.
!
  sample_num = min ( 1000, n )

  allocate ( sample(1:m,1:sample_num) )

  x1 = - ( r1 + h )
  x2 =   ( r1 + h )
  y1 = - ( r1 + h )
  y2 =   ( r1 + h )

  do
!
!  Generate a batch of points in [0,1]x[0,1].
!
    call r8mat_uniform_01 ( m, sample_num, seed, sample )
!
!  Remap the points to [x1,x2] x [y1,y2].
!
    sample(1,1:sample_num) = x1 + sample(1,1:sample_num) * ( x2 - x1 )
    sample(2,1:sample_num) = y1 + sample(2,1:sample_num) * ( y2 - y1 )
!
!  Accept each point that is inside the enlarged region.
!
    do j = 1, sample_num

      inside = hexagon_contains_point_2d ( v1, sample(1:2,j) ) &
        .and. ( .not. &
               hexagon_contains_point_2d ( v2, sample(1:2,j) ) )

      if ( inside ) then

        have = have + 1
        point(1:m,have) = sample(1:m,j)

        if ( have == n ) then
          return
        end if

      end if

    end do

  end do

  deallocate ( sample )

  return
end
subroutine p04_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P04_SDIST returns the signed distance to the region in problem 04.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P04_SDIST - Fatal error!'
  write ( *, '(a)' ) '  The routine for this test is not written yet!'
  stop
end
subroutine p04_title ( title )

!*****************************************************************************80
!
!! P04_TITLE returns a title for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = '#4: The unit hexagon with hexagonal hole.'

  return
end
subroutine p05_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P05_BOUNDARY_NEAREST returns a nearest boundary point in problem 05.
!
!  Discussion:
!
!    The given input point need not be inside the region.
!
!    In some cases, more than one boundary point may be "nearest",
!    but only one will be returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( m, n ) :: boundary
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/  0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ -0.4D+00, 0.0D+00 /)
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_min
  integer ( kind = 4 ) j
  real ( kind = 8 ) norm
  real ( kind = 8 ), dimension ( m ) :: p1
  real ( kind = 8 ), dimension ( m ) :: p2
  real ( kind = 8 ), dimension ( m ) :: pn
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ), dimension ( m ) :: q1
  real ( kind = 8 ), dimension ( m ) :: q2
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.55D+00
  real ( kind = 8 ) t

  do j = 1, n

    dist_min = huge ( dist_min )
!
!  Distance to the semicircle #1.
!
    if ( center1(2) <= point(2,j) ) then

      norm = sqrt ( ( point(1,j) - center1(1) )**2 &
                  + ( point(2,j) - center1(2) )**2 )

      if ( 0.0D+00 == norm ) then

        pn(1) = center1(1) + r1 * sqrt ( 0.5D+00 )
        pn(2) = center1(2) + r1 * sqrt ( 0.5D+00 )

      else

        pn(1) = center1(1) + ( point(1,j) - center1(1) ) / norm
        pn(2) = center1(2) + ( point(2,j) - center1(2) ) / norm

      end if

    else if ( point(1,j) <= center1(1) ) then

      pn(1) =  center1(1) - r1
      pn(2) =  center1(2)

    else if ( center1(1) <= point(1,j) ) then

      pn(1) =  center1(1) + r1
      pn(2) =  center1(2)

    end if

    dist = sqrt ( ( point(1,j) - pn(1) )**2 + ( point(2,j) - pn(2) )**2 )

    if ( dist < dist_min ) then
      dist_min = dist
      boundary(1:2,j) = pn(1:2)
    end if
!
!  Distance to semicircle #2.
!
    if ( center2(2) <= point(2,j) ) then

      norm = sqrt ( ( point(1,j) - center2(1) )**2 &
                  + ( point(2,j) - center2(2) )**2 )

      if ( 0.0D+00 == norm ) then

        pn(1) = center2(1) + r2 * sqrt ( 0.5D+00 )
        pn(2) = center2(2) + r2 * sqrt ( 0.5D+00 )

      else

        pn(1) = center2(1) + ( point(1,j) - center2(1) ) / norm
        pn(2) = center2(2) + ( point(2,j) - center2(2) ) / norm

      end if

    else if ( point(1,j) <= center2(1) ) then

      pn(1) =  center2(1) - r2
      pn(2) =  center2(2)

    else if ( center2(1) <= point(1,j) ) then

      pn(1) =  center2(1) + r2
      pn(2) =  center2(2)

    end if

    dist = sqrt ( ( point(1,j) - pn(1) )**2 + ( point(2,j) - pn(2) )**2 )

    if ( dist < dist_min ) then
      dist_min = dist
      boundary(1:2,j) = pn(1:2)
    end if
!
!  Distance to line segment #1: (P1,P2).
!
    p1(1:2) = (/ center2(1) - r1, center2(1) /)
    p2(1:2) = (/ center2(1) - r2, center2(2) /)

    call segment_point_near_2d ( p1, p2, point(1:2,j), pn, dist, t )

    if ( dist < dist_min ) then
      dist_min = dist
      boundary(1:2,j) = pn(1:2)
    end if
!
!  Distance to line segment #2: (Q1,Q2).
!
    q1(1:2) = (/ center2(1) + r2, center2(2) /)
    q2(1:2) = (/ center1(1) + r1, center1(2) /)

    call segment_point_near_2d ( q1, q2, point(1:2,j), pn, dist, t )

    if ( dist < dist_min ) then
      dist_min = dist
      boundary(1:2,j) = pn(1:2)
    end if

  end do

  return
end
subroutine p05_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P05_BOUNDARY_PROJECT projects exterior points to the boundary in problem 05.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( m ) :: boundary
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/  0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ -0.4D+00, 0.0D+00 /)
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_min
  logical inside(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) norm
  real ( kind = 8 ), dimension ( m ) :: p1
  real ( kind = 8 ), dimension ( m ) :: p2
  real ( kind = 8 ), dimension ( m ) :: pn
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ), dimension ( m ) :: q1
  real ( kind = 8 ), dimension ( m ) :: q2
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.55D+00
  real ( kind = 8 ) t

  call p05_inside ( m, n, point, inside )

  do j = 1, n

    if ( inside(j) ) then
      cycle
    end if

    dist_min = huge ( dist_min )
!
!  Distance to the semicircle #1.
!
    if ( center1(2) <= point(2,j) ) then

      norm = sqrt ( ( point(1,j) - center1(1) )**2 &
                  + ( point(2,j) - center1(2) )**2 )

      if ( 0.0D+00 == norm ) then

        pn(1) = center1(1) + r1 * sqrt ( 0.5D+00 )
        pn(2) = center1(2) + r1 * sqrt ( 0.5D+00 )

      else

        pn(1) = center1(1) + r1 * ( point(1,j) - center1(1) ) / norm
        pn(2) = center1(2) + r1 * ( point(2,j) - center1(2) ) / norm

      end if

    else if ( point(1,j) <= center1(1) ) then

      pn(1) =  center1(1) - r1
      pn(2) =  center1(2)

    else if ( center1(1) <= point(1,j) ) then

      pn(1) =  center1(1) + r1
      pn(2) =  center1(2)

    end if

    dist = sqrt ( ( point(1,j) - pn(1) )**2 + ( point(2,j) - pn(2) )**2 )

    dist_min = dist
    boundary(1:2) = pn(1:2)
!
!  Distance to semicircle #2.
!
    if ( center2(2) <= point(2,j) ) then

      norm = sqrt ( ( point(1,j) - center2(1) )**2 &
                  + ( point(2,j) - center2(2) )**2 )

      if ( 0.0D+00 == norm ) then

        pn(1) = center2(1) + r2 * sqrt ( 0.5D+00 )
        pn(2) = center2(2) + r2 * sqrt ( 0.5D+00 )

      else

        pn(1) = center2(1) + r2 * ( point(1,j) - center2(1) ) / norm
        pn(2) = center2(2) + r2 * ( point(2,j) - center2(2) ) / norm

      end if

    else if ( point(1,j) <= center2(1) ) then

      pn(1) =  center2(1) - r2
      pn(2) =  center2(2)

    else if ( center2(1) <= point(1,j) ) then

      pn(1) =  center2(1) + r2
      pn(2) =  center2(2)

    end if

    dist = sqrt ( ( point(1,j) - pn(1) )**2 + ( point(2,j) - pn(2) )**2 )

    if ( dist < dist_min ) then
      dist_min = dist
      boundary(1:2) = pn(1:2)
    end if
!
!  Distance to line segment #1: (P1,P2).
!
    p1(1:2) = (/ center1(1) - r1, center1(2) /)
    p2(1:2) = (/ center2(1) - r2, center2(2) /)

    call segment_point_near_2d ( p1, p2, point(1:2,j), pn, dist, t )

    if ( dist < dist_min ) then
      dist_min = dist
      boundary(1:2) = pn(1:2)
    end if
!
!  Distance to line segment #2: (Q1,Q2).
!
    q1(1:2) = (/ center2(1) + r2, center2(2) /)
    q2(1:2) = (/ center1(1) + r1, center1(2) /)

    call segment_point_near_2d ( q1, q2, point(1:2,j), pn, dist, t )

    if ( dist < dist_min ) then
      dist_min = dist
      boundary(1:2) = pn(1:2)
    end if

    point(1:2,j) = boundary(1:2)

  end do

  return
end
subroutine p05_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P05_BOUNDARY_SEGMENT returns a boundary segment in problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in
!    the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/  0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ -0.4D+00, 0.0D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.55D+00
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) segment(m,segment_length)

  if ( segment_index == 1 ) then
!
!  Work out the appropriate segment lengths, and then
!  adjust N4, if necessary, to account for roundoff.
!
    n4 = nint ( real ( segment_length - 1, kind = 8 ) &
            / ( 1.0D+00 + r2 + 0.90D+00 / pi )  )

    n1 = nint ( 0.05D+00 * real ( n4, kind = 8 ) / pi )
    n2 = nint (       r2 * real ( n4, kind = 8 )            )
    n3 = nint ( 0.85D+00 * real ( n4, kind = 8 ) / pi )

    n4 = segment_length - 1 - n1 - n2 - n3

    j = 0
!
!  Piece #1, the short straight piece.
!
    do i = 1, n1
      j = j + 1
      segment(1:2,j) = &
        ( real ( n1 - i + 1, kind = 8 ) * (/ center1(1) - r1, center1(2) /)   &
        + real (      i - 1, kind = 8 ) * (/ center2(1) - r2, center2(2) /) ) &
        / real ( n1,         kind = 8 )
    end do
!
!  Piece #2, the smaller semicircle.
!
    do i = 1, n2
      angle = real ( n2 - i + 1, kind = 8 ) * pi &
            / real ( n2,         kind = 8 )
      j = j + 1
      segment(1:2,j) = (/ center2(1) + r2 * cos ( angle ), &
                          center2(2) + r2 * sin ( angle ) /)
    end do
!
!  Piece #3, the long straight piece.
!
    do i = 1, n3
      j = j + 1
      segment(1:2,j) = &
        ( real ( n3 - i + 1, kind = 8 ) * (/ center2(1) + r2, center2(2) /)   &
        + real (      i - 1, kind = 8 ) * (/ center1(1) + r1, center1(1) /) ) &
        / real ( n3,         kind = 8 )
    end do
!
!  Piece #4, the larger semicircle.
!
    do i = 1, n4
      angle = ( i - 1 ) * pi / real ( n4, kind = 8 )
      j = j + 1
      segment(1:2,j) = (/ center1(1) + r1 * cos ( angle ), &
                          center1(2) + r1 * sin ( angle ) /)
    end do

    j = j + 1
    segment(1:2,j) = (/ center1(1) - r1, center1(2) /)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P05_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p05_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P05_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the
!    boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in
!    the segment.
!
  implicit none

  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/  0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ -0.4D+00, 0.0D+00 /)
  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r1 = 1.00D+00
  real ( kind = 8 ), parameter :: r2 = 0.55D+00
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P05_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    n = nint ( ( pi * ( r1  + r2 ) &
      + ( ( center2(1) - r2 ) - ( center1(1) - r1 ) ) &
      + ( ( center1(1) + r1 ) - ( center2(1) + r2 ) ) ) / h )
    n = max ( n, 21 )
    segment_length = n

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P05_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p05_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P05_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of
!    boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 1

  return
end
subroutine p05_box ( m, lo, hi )

!*****************************************************************************80
!
!! P05_BOX returns a bounding box for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/  0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ -0.4D+00, 0.0D+00 /)
  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.55D+00

  lo(1:m) = (/ center1(1) - r1, center1(2)      /)
  hi(1:m) = (/ center1(1) + r1, center1(2) + r1 /)

  return
end
subroutine p05_density ( m, n, point, density )

!*****************************************************************************80
!
!! P05_DENSITY returns the density for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p05_element_size ( element_size )

!*****************************************************************************80
!
!! P05_ELEMENT_SIZE returns a typical element size for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.015D+00

  return
end
subroutine p05_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P05_FIXED_NUM returns the number of fixed points in problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 4

  return
end
subroutine p05_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P05_FIXED_POINTS returns the fixed points in problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/  0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ -0.4D+00, 0.0D+00 /)
  real ( kind = 8 ) fixed(m,fixed_num)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.55D+00

  fixed(1:2,1) = (/ center1(1) - r1,  center1(2) /)
  fixed(1:2,2) = (/ center2(1) - r2,  center2(2) /)
  fixed(1:2,3) = (/ center2(1) + r2,  center2(2) /)
  fixed(1:2,4) = (/ center1(1) + r1,  center1(2) /)

  return
end
subroutine p05_header ( )

!*****************************************************************************80
!
!! P05_HEADER prints some information about problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  integer ( kind = 4 ) boundary_segment_num
  real ( kind = 8 ), dimension ( m ) :: center1 = (/  0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( m ) :: center2 = (/ -0.4D+00, 0.0D+00 /)
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.55D+00

  call p05_boundary_segment_num ( boundary_segment_num )
  call p05_fixed_num ( fixed_num )
  call p05_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P05:'
  write ( *, '(a)' ) '  Strang and Persson example #5'
  write ( *, '(a)' ) '  The horn.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  Circle C1 has center = (0,0)'
  write ( *, '(a,g14.6)' ) '  Radius R1 = ', r1
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  Circle C2 has center = (-0.4,0)'
  write ( *, '(a,g14.6)' ) '  Radius R2 = ', r2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Points in the region are:'
  write ( *, '(a)' ) '    in C1 and not in C2 and have 0 <= Y.'
  write ( *, '(a)' ) '  A uniform mesh density is requested.'
  write ( *, '(a)' ) '  Element sizes tried were 0.4, 0.2, 0.1.'
  write ( *, '(a)' )       ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p05_hole_num ( hole_num )

!*****************************************************************************80
!
!! P05_HOLE_NUM counts the holes in problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 0

  return
end
subroutine p05_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P05_HOLE_POINT returns a point inside a given hole in problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  return
end
subroutine p05_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P05_INSIDE reports if a point is inside the region in problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/  0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ -0.4D+00, 0.0D+00 /)
  logical inside(n)
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 =   1.00D+00
  real ( kind = 8 ), parameter :: r2 =   0.55D+00

  inside(1:n) =                                                 &
       center1(2) <=   point(2,1:n)                             &
       .and.                                                    &
                     ( point(1,1:n) - center1(1) )**2           &
                   + ( point(2,1:n) - center1(2) )**2  <= r1**2 &
       .and.                                                    &
       r2**2 <=      ( point(1,1:n) - center2(1) )**2           &
                   + ( point(2,1:n) - center2(2) )**2


  return
end
subroutine p05_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P05_SAMPLE samples points from the region in problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/  0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ -0.4D+00, 0.0D+00 /)
  integer ( kind = 4 ) have
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 =   1.00D+00
  real ( kind = 8 ), parameter :: r2 =   0.55D+00
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  x1 = center1(1) - r1
  x2 = center1(1) + r1

  y1 = center1(2)
  y2 = center1(2) + r1

  have = 0
!
!  We are going to generate batches of sample points.
!
  sample_num = min ( 1000, n )

  allocate ( sample(1:m,1:sample_num) )

  do
!
!  Generate a batch of points in [0,1]x[0,1].
!
    call r8mat_uniform_01 ( m, sample_num, seed, sample )
!
!  Remap the points to the box [X1,X2] x [Y1,Y2].
!
    sample(1,1:sample_num) = x1 + sample(1,1:sample_num) * ( x2 - x1 )
    sample(2,1:sample_num) = y1 + sample(2,1:sample_num) * ( y2 - y1 )
!
!  Accept those points which are in the big circle and not in the
!  small circle.
!
    do j = 1, sample_num

      if (                                                          &
                   ( sample(1,j) - center1(1) )**2                  &
                 + ( sample(2,j) - center1(2) )**2 <= r1**2 .and.   &
        r2**2 <=   ( sample(1,j) - center2(1) )**2                  &
                 + ( sample(2,j) - center2(2) )**2 ) then

        have = have + 1
        point(1:m,have) = sample(1:m,j)

        if ( have == n ) then
          return
        end if

      end if

    end do

  end do

  deallocate ( sample )

  return
end
subroutine p05_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P05_SAMPLE_H1 samples points from the enlarged region in problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/  0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ -0.4D+00, 0.0D+00 /)
  real ( kind = 8 ) h
  integer ( kind = 4 ) have
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 =   1.00D+00
  real ( kind = 8 ) r1h
  real ( kind = 8 ), parameter :: r2 =   0.55D+00
  real ( kind = 8 ) r2h
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  r1h = r1 + h
  r2h = max ( r2 - h, 0.0D+00 )

  x1 = center1(1) - r1h
  x2 = center1(1) + r1h

  y1 = center1(2) - h
  y2 = center1(2) + r1h

  have = 0
!
!  We are going to generate batches of sample points.
!
  sample_num = min ( 1000, n )

  allocate ( sample(1:m,1:sample_num) )

  do
!
!  Generate a batch of points in [0,1]x[0,1].
!
    call r8mat_uniform_01 ( m, sample_num, seed, sample )
!
!  Remap the points to the box [X1,X2] x [Y1,Y2].
!
    sample(1,1:sample_num) = x1 + sample(1,1:sample_num) * ( x2 - x1 )
    sample(2,1:sample_num) = y1 + sample(2,1:sample_num) * ( y2 - y1 )

    do j = 1, sample_num
!
!  A) above the origin, and in the big circle and not in the small circle...
!
      if ( center1(2) <= sample(2,j) ) then

        if (                                                            &
                      ( sample(1,j) - center1(1) )**2                   &
                    + ( sample(2,j) - center1(2) )**2 <= r1h**2 .and.   &
          r2h**2 <=   ( sample(1,j) - center2(1) )**2                   &
                    + ( sample(2,j) - center2(2) )**2 ) then

          have = have + 1
          point(1:m,have) = sample(1:m,j)

          if ( have == n ) then
            return
          end if

        end if
!
!  B) or below the origin, but underneath one of the two segments.
!
      else if ( sample(2,j) <= center1(2) ) then

        if (                                         &
             ( center1(1) - r1h <= sample(1,j) .and. &
               sample(1,j) <= center2(1) - r2h )     &
             .or.                                    &
             ( center2(1) + r2h <= sample(1,j) .and. &
               sample(1,j) <= center1(1) + r1h )     &
          ) then

          have = have + 1
          point(1:m,have) = sample(1:m,j)

          if ( have == n ) then
            return
          end if
        end if

      end if

    end do

  end do

  deallocate ( sample )

  return
end
subroutine p05_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P05_SDIST returns the signed distance to the region in problem 05.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P05_SDIST - Fatal error!'
  write ( *, '(a)' ) '  The routine for this test is not written yet!'
  stop
end
subroutine p05_title ( title )

!*****************************************************************************80
!
!! P05_TITLE returns a title for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = '#5: The horn.'

  return
end
subroutine p06_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P06_BOUNDARY_NEAREST returns a nearest boundary point in problem 06.
!
!  Discussion:
!
!    The given input point need not be inside the region.
!
!    In some cases, more than one boundary point may be "nearest",
!    but only one will be returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) atan4
  real ( kind = 8 ), dimension ( m, n ) :: boundary
  real ( kind = 8 ) cm
  real ( kind = 8 ) cs
  real ( kind = 8 ) dstar1
  real ( kind = 8 ) dstar2
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  real ( kind = 8 ) sm
  real ( kind = 8 ) ss
  integer ( kind = 4 ) status
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) tstar1
  real ( kind = 8 ) tstar2
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  do j = 1, n

    x = point(1,j)
    y = point(2,j)

    if ( x == 0.0D+00 .and. y == 0.0D+00 ) then

      boundary(1:2,j) = (/ r2, 0.0D+00 /)
!
!  Determine the angle formed by (0,0) and the point.
!
    else

      t = atan4 ( y, x )
!
!  Find the nearest point on the superellipse x^4 + y^4 = 1^4.
!
      t1 = t - pi / 4.0D+00
      t2 = t + pi / 4.0D+00

      status = 0

      do

        call fmin_rc ( t1, t2, tstar1, status, dstar1 )

        if ( status == 0 ) then
          exit
        end if

        cm = abs (           cos ( tstar1 ) )
        cs = sign ( 1.0D+00, cos ( tstar1 ) )
        sm = abs (           sin ( tstar1 ) )
        ss = sign ( 1.0D+00, sin ( tstar1 ) )

        dstar1 = ( x - r1 * cs * sqrt ( cm ) )**2 &
               + ( y - r1 * ss * sqrt ( sm ) )**2

      end do

      boundary(1:2,j) = (/ r1 * cs * sqrt ( cm ), &
                         r1 * ss * sqrt ( sm ) /)
!
!  Find the nearest point on the superellipse x^4 + y^4 = 1/2^4.
!
      t1 = t - pi / 4.0D+00
      t2 = t + pi / 4.0D+00

      status = 0

      do

        call fmin_rc ( t1, t2, tstar2, status, dstar2 )

        if ( status == 0 ) then
          exit
        end if

        cm = abs (           cos ( tstar2 ) )
        cs = sign ( 1.0D+00, cos ( tstar2 ) )
        sm = abs (           sin ( tstar2 ) )
        ss = sign ( 1.0D+00, sin ( tstar2 ) )

        dstar2 = ( x - r2 * cs * sqrt ( cm ) )**2 &
               + ( y - r2 * ss * sqrt ( sm ) )**2

      end do

      if ( dstar2 < dstar1 ) then
        boundary(1:2,j) = (/ r2 * cs * sqrt ( cm ), &
                           r2 * ss * sqrt ( sm ) /)
      end if

    end if

  end do

  return
end
subroutine p06_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P06_BOUNDARY_NEAREST projects exterior points to the boundary in problem 06.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) atan4
  real ( kind = 8 ) cm
  real ( kind = 8 ) cs
  real ( kind = 8 ) dstar1
  real ( kind = 8 ) dstar2
  logical inside(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  real ( kind = 8 ) sm
  real ( kind = 8 ) ss
  integer ( kind = 4 ) status
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) temp(m)
  real ( kind = 8 ) tstar1
  real ( kind = 8 ) tstar2
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  call p06_inside ( m, n, point, inside )

  do j = 1, n

    if ( inside(j) ) then
      cycle
    end if

    x = point(1,j)
    y = point(2,j)

    if ( x == 0.0D+00 .and. y == 0.0D+00 ) then

      temp(1:m) = (/ r2, 0.0D+00 /)
!
!  Determine the angle formed by (0,0) and the point.
!
    else

      t = atan4 ( y, x )
!
!  Find the nearest point on the superellipse x^4 + y^4 = 1^4.
!
      t1 = t - pi / 4.0D+00
      t2 = t + pi / 4.0D+00

      status = 0

      do

        call fmin_rc ( t1, t2, tstar1, status, dstar1 )

        if ( status == 0 ) then
          exit
        end if

        cm = abs (           cos ( tstar1 ) )
        cs = sign ( 1.0D+00, cos ( tstar1 ) )
        sm = abs (           sin ( tstar1 ) )
        ss = sign ( 1.0D+00, sin ( tstar1 ) )

        dstar1 = ( x - r1 * cs * sqrt ( cm ) )**2 &
               + ( y - r1 * ss * sqrt ( sm ) )**2

      end do

      temp(1:m) = (/ r1 * cs * sqrt ( cm ), &
                        r1 * ss * sqrt ( sm ) /)
!
!  Find the nearest point on the superellipse x^4 + y^4 = 1/2^4.
!
      t1 = t - pi / 4.0D+00
      t2 = t + pi / 4.0D+00

      status = 0

      do

        call fmin_rc ( t1, t2, tstar2, status, dstar2 )

        if ( status == 0 ) then
          exit
        end if

        cm = abs (           cos ( tstar2 ) )
        cs = sign ( 1.0D+00, cos ( tstar2 ) )
        sm = abs (           sin ( tstar2 ) )
        ss = sign ( 1.0D+00, sin ( tstar2 ) )

        dstar2 = ( x - r2 * cs * sqrt ( cm ) )**2 &
               + ( y - r2 * ss * sqrt ( sm ) )**2

      end do

      if ( dstar2 < dstar1 ) then
        temp(1:m) = (/ r2 * cs * sqrt ( cm ), &
                          r2 * ss * sqrt ( sm ) /)
      end if

    end if

    point(1:m,j) = temp(1:m)

  end do

  return
end
subroutine p06_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P06_BOUNDARY_SEGMENT returns a boundary segment in problem 06.
!
!  Discussion:
!
!    The points are not uniformly spaced.  However, they correspond
!    to a uniform angular spacing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  real ( kind = 8 ) alpha
  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension(2) :: center = (/ 0.0D+00, 0.0D+00 /)
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) segment(m,segment_length)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( segment_index == 1 ) then
    r = r1
  else if ( segment_index == 2 ) then
    r = r2
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P06_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop
  end if

  do j = 1, segment_length

    angle = ( 2.0D+00 * pi * real (              j - 1, kind = 8 ) ) &
                           / real ( segment_length - 1, kind = 8 )
    x = r * cos ( angle )
    y = r * sin ( angle )

    alpha = ( r**4 / ( x**4 + y**4 ) )**0.25D+00

    segment(1,j) = center(1) + alpha * x
    segment(2,j) = center(2) + alpha * y

  end do

  return
end
subroutine p06_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P06_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 06.
!
!  Discussion:
!
!    For this region, the boundary points will not be evenly spaced,
!    and the value of SEGMENT_LENGTH returned will only approximately
!    guarantee that the maximal spacing is H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the
!    boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in
!    the segment.
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P06_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    n = nint ( 2.0D+00 * pi * r1 / h )
    n = max ( n, 5 )
    segment_length = n

  else if ( segment_index == 2 ) then

    n = nint ( 2.0D+00 * pi * r2 / h )
    n = max ( n, 5 )
    segment_length = n

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P06_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop
  end if

  return
end
subroutine p06_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P06_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of
!    boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 2

  return
end
subroutine p06_box ( m, lo, hi )

!*****************************************************************************80
!
!! P06_BOX returns a bounding box for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00

  lo(1:m) = (/ -r1, -r1 /)
  hi(1:m) = (/ +r1, +r1 /)

  return
end
subroutine p06_density ( m, n, point, density )

!*****************************************************************************80
!
!! P06_DENSITY returns the density for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = sqrt ( sqrt ( &
    (   point(1,1:n)**4 &
      + point(2,1:n)**4 ) ) )

  return
end
subroutine p06_element_size ( element_size )

!*****************************************************************************80
!
!! P06_ELEMENT_SIZE returns a typical element size for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.08D+00

  return
end
subroutine p06_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P06_FIXED_NUM returns the number of fixed points in problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 0

  return
end
subroutine p06_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P06_FIXED_POINTS returns the fixed points in problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) fixed(m,fixed_num)

  return
end
subroutine p06_header ( )

!*****************************************************************************80
!
!! P06_HEADER prints some information about problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00

  call p06_boundary_segment_num ( boundary_segment_num )
  call p06_fixed_num ( fixed_num )
  call p06_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P06:'
  write ( *, '(a)' ) '  Strang and Persson example #6'
  write ( *, '(a)' ) '  Superellipse with superelliptical hole.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Radius R1 = ', r1
  write ( *, '(a,g14.6)' ) '  Radius R2 = ', r2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A uniform mesh density is requested.'
  write ( *, '(a)' ) '  Element sizes tried were 0.4, 0.2, 0.1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p06_hole_num ( hole_num )

!*****************************************************************************80
!
!! P06_HOLE_NUM counts the holes in problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 1

  return
end
subroutine p06_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P06_HOLE_POINT returns a point inside a given hole in problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  hole_point(1:2) = (/ 0.0D+00, 0.0D+00 /)

  return
end
subroutine p06_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P06_INSIDE reports if a point is inside the region in problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical inside(n)
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00

  inside(1:n) =                                &
    (                                                  &
                point(1,1:n)**4                &
              + point(2,1:n)**4 <= r1**4       &
    )                                                  &
    .and.                                              &
    (                                                  &
      r2**4 <= point(1,1:n)**4                 &
             + point(2,1:n)**4                 &
    )

  return
end
subroutine p06_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P06_SAMPLE samples points from the region in problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) have
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  have = 0
!
!  We are going to generate batches of sample points.
!
  sample_num = min ( 1000, n )

  allocate ( sample(1:m,1:sample_num) )

  x1 = -r1
  x2 = +r1
  y1 = -r1
  y2 = +r1

  do
!
!  Generate a batch of points in [0,1]x[0,1].
!
    call r8mat_uniform_01 ( m, sample_num, seed, sample )
!
!  Remap the points to the box [x1,x2] x [y1,y2].
!
    sample(1,1:sample_num) = x1 + sample(1,1:sample_num) * ( x2 - x1 )
    sample(2,1:sample_num) = y1 + sample(2,1:sample_num) * ( y2 - y1 )
!
!  Accept those points which are in the big superellipse and not in the
!  small one.
!
    do j = 1, sample_num

      if (       sample(1,j)**4 + sample(2,j)**4 <= r1**4 .and. &
        r2**4 <= sample(1,j)**4 + sample(2,j)**4          ) then

        have = have + 1
        point(1:m,have) = sample(1:m,j)

        if ( have == n ) then
          return
        end if

      end if

    end do

  end do

  deallocate ( sample )

  return
end
subroutine p06_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P06_SAMPLE_H1 samples points from the enlarged region in problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) h
  integer ( kind = 4 ) have
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.5D+00
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ) test
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  have = 0
!
!  We are going to generate batches of sample points.
!
  sample_num = min ( 1000, n )

  allocate ( sample(1:m,1:sample_num) )

  x1 = - ( r1 + h )
  x2 = + ( r1 + h )
  y1 = - ( r1 + h )
  y2 = + ( r1 + h )

  do
!
!  Generate a batch of points in [0,1]x[0,1].
!
    call r8mat_uniform_01 ( m, sample_num, seed, sample )
!
!  Remap the points to the box [x1,x2] x [y1,y2].
!
    sample(1,1:sample_num) = x1 + sample(1,1:sample_num) * ( x2 - x1 )
    sample(2,1:sample_num) = y1 + sample(2,1:sample_num) * ( y2 - y1 )
!
!  Accept those points which are in the big superellipse and not in the
!  small one.
!
    do j = 1, sample_num

      test = sample(1,j)**4 + sample(2,j)**4

      if (                  test <= ( r1 + h )**4 .and. &
           ( r2 - h )**4 <= test                 ) then

        have = have + 1
        point(1:m,have) = sample(1:m,j)

        if ( have == n ) then
          return
        end if

      end if

    end do

  end do

  deallocate ( sample )

  return
end
subroutine p06_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P06_SDIST returns the signed distance to the region in problem 06.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P06_SDIST - Fatal error!'
  write ( *, '(a)' ) '  The routine for this test is not written yet!'
  stop
end
subroutine p06_title ( title )

!*****************************************************************************80
!
!! P06_TITLE returns a title for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = '#6: The superellipse with superelliptical hole.'

  return
end
subroutine p07_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P07_BOUNDARY_NEAREST returns a nearest boundary point in problem 07.
!
!  Discussion:
!
!    The given input point need not be inside the region.
!
!    In some cases, more than one boundary point may be "nearest",
!    but only one will be returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( m, n ) :: boundary
  real ( kind = 8 ) dstar1
  real ( kind = 8 ) dstar2
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), dimension ( m, n ) :: point
  integer ( kind = 4 ) status
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xstar
  real ( kind = 8 ) y
  real ( kind = 8 ) ystar

  do j = 1, n

    x = point(1,j)
    y = point(2,j)
!
!  Examine the upper boundary.
!
!  X can be used as the parameter.
!  The upper boundary can be broken into three parts.  Determine which one
!  we are to look at.
!
    if ( x <= -pi ) then
      x1 = -2.5D+00 * pi
      x2 = - pi
    else if ( x <= pi ) then
      x1 = - pi
      x2 =   pi
    else
      x1 =   pi
      x2 =   2.5D+00 * pi
    end if

    status = 0

    do

      call fmin_rc ( x1, x2, xstar, status, dstar1 )

      ystar = cos ( xstar )

      if ( status == 0 ) then
        exit
      end if

      dstar1 = ( x - xstar )**2 + ( y - ystar )**2

    end do

    boundary(1:2,j) = (/ xstar, ystar /)
!
!  Examine the lower boundary.
!
    x1 = -2.5D+00 * pi
    x2 =  2.5D+00 * pi

    status = 0

    do

      call fmin_rc ( x1, x2, xstar, status, dstar2 )

      ystar = 5.0D+00 * ( xstar / ( 2.5D+00 * pi ) )**4 - 5.0D+00

      if ( status == 0 ) then
        exit
      end if

      dstar2 = ( x - xstar )**2 + ( y - ystar )**2

    end do

    if ( dstar2 < dstar1 ) then
      boundary(1:2,j) = (/ xstar, ystar /)
    end if

  end do

  return
end
subroutine p07_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P07_BOUNDARY_PROJECT projects exterior points to the boundary in problem 07.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) dstar1
  real ( kind = 8 ) dstar2
  logical inside(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), dimension ( m, n ) :: point
  integer ( kind = 4 ) status
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xstar
  real ( kind = 8 ) y
  real ( kind = 8 ) ystar

  call p07_inside ( m, n, point, inside )

  do j = 1, n

    if ( inside(j) ) then
      cycle
    end if

    x = point(1,j)
    y = point(2,j)
!
!  Examine the upper boundary.
!
!  X can be used as the parameter.
!  The upper boundary can be broken into three parts.  Determine which one
!  we are to look at.
!
    if ( x <= -pi ) then
      x1 = -2.5D+00 * pi
      x2 = - pi
    else if ( x <= pi ) then
      x1 = - pi
      x2 =   pi
    else
      x1 =   pi
      x2 =   2.5D+00 * pi
    end if

    status = 0

    do

      call fmin_rc ( x1, x2, xstar, status, dstar1 )

      ystar = cos ( xstar )

      if ( status == 0 ) then
        exit
      end if

      dstar1 = ( x - xstar )**2 + ( y - ystar )**2

    end do

    point(1:2,j) = (/ xstar, ystar /)
!
!  Examine the lower boundary.
!
    x1 = -2.5D+00 * pi
    x2 =  2.5D+00 * pi

    status = 0

    do

      call fmin_rc ( x1, x2, xstar, status, dstar2 )

      ystar = 5.0D+00 * ( xstar / ( 2.5D+00 * pi ) )**4 - 5.0D+00

      if ( status == 0 ) then
        exit
      end if

      dstar2 = ( x - xstar )**2 + ( y - ystar )**2

    end do

    if ( dstar2 < dstar1 ) then
      point(1:2,j) = (/ xstar, ystar /)
    end if

  end do

  return
end
subroutine p07_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P07_BOUNDARY_SEGMENT returns a boundary segment in problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) segment(m,segment_length)
  real ( kind = 8 ) x
  real ( kind = 8 ) x_hi
  real ( kind = 8 ) x_lo
  real ( kind = 8 ) y

  if ( segment_index == 1 ) then

    n1 = ( segment_length - 3 ) / 2
    n2 = segment_length - 3 - n1

    c = 2.5D+00 * pi

    x_lo = -c
    x_hi =  c

    j = 0

    j = j + 1
    x = x_lo
    y = 0.0D+00
    segment(1:2,j) = (/ x, y /)

    do i = 1, n1
      x = ( real ( n1 - i + 1, kind = 8 ) * x_lo   &
          + real (      i,     kind = 8 ) * x_hi ) &
          / real ( n1     + 1, kind = 8 )
      y = cos ( x )
      j = j + 1
      segment(1:2,j) = (/ x, y /)
    end do

    j = j + 1
    x = x_hi
    y = 0.0D+00
    segment(1:2,j) = (/ x, y /)

    do i = 1, n2
      x = ( real ( n2 - i + 1, kind = 8 ) * x_hi   &
          + real (      i,     kind = 8 ) * x_lo ) &
          / real ( n2     + 1, kind = 8 )
      y = 5.0D+00 * ( x / c )**4 - 5.0D+00
      j = j + 1
      segment(1:2,j) = (/ x, y /)
    end do

    j = j + 1
    x = x_lo
    y = 0.0D+00
    segment(1:2,j) = (/ x, y /)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P07_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p07_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P07_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 07.
!
!  Discussion:
!
!    No attempt has been made here to accurately compute a value of N
!    which would guarantee that the boundary would be divided into pieces
!    of length no more than H.  The curve is a little too complicated
!    to make this easy to do.
!
!    Moreover, the points that will be generated will only be equally
!    spaced in their X argument, not in their arc length.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the
!    boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in
!    the segment.
!
  implicit none

  real ( kind = 8 )  h
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P07_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    n = nint ( 10.0D+00 * pi / h )
    n = max ( n, 13 )
    segment_length = n

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P07_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p07_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P07_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of
!    boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 1

  return
end
subroutine p07_box ( m, lo, hi )

!*****************************************************************************80
!
!! P07_BOX returns a bounding box for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) temp

  temp = 5.0D+00 * pi / 2.0D+00

  lo(1:m) = (/ -temp, -5.0D+00 /)
  hi(1:m) = (/ +temp, +1.0D+00 /)

  return
end
subroutine p07_density ( m, n, point, density )

!*****************************************************************************80
!
!! P07_DENSITY returns the density for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p07_element_size ( element_size )

!*****************************************************************************80
!
!! P07_ELEMENT_SIZE returns a typical element size for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.075D+00

  return
end
subroutine p07_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P07_FIXED_NUM returns the number of fixed points in problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 2

  return
end
subroutine p07_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P07_FIXED_POINTS returns the fixed points in problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) fixed(m,fixed_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) temp

  temp = 5.0D+00 * pi / 2.0D+00

  fixed(1:2,1) = (/ -temp, 0.0D+00 /)
  fixed(1:2,2) = (/  temp, 0.0D+00 /)

  return
end
subroutine p07_header ( )

!*****************************************************************************80
!
!! P07_HEADER prints some information about problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num

  call p07_boundary_segment_num ( boundary_segment_num )
  call p07_fixed_num ( fixed_num )
  call p07_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P07:'
  write ( *, '(a)' ) '  Strang and Persson example #7'
  write ( *, '(a)' ) '  Bicycle seat (implicit).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A uniform mesh density is requested.'
  write ( *, '(a)' ) '  The boundary is formed by two algebraic expressions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p07_hole_num ( hole_num )

!*****************************************************************************80
!
!! P07_HOLE_NUM counts the holes in problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 0

  return
end
subroutine p07_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P07_HOLE_POINT returns a point inside a given hole in problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  return
end
subroutine p07_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P07_INSIDE reports if a point is inside the region in problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical inside(n)
  real ( kind = 8 ) hi(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) point(m,n)

  inside(1:n) = .true.

  call p07_box ( m, lo, hi )
!
!  Check whether points are in the bounding box.
!
  do j = 1, n

    do i = 1, m
      if ( point(i,j) < lo(i) .or. hi(i) < point(i,j) ) then
        inside(j) = .false.
        cycle
      end if

    end do

  end do
!
!  Check whether points in the bounding box are in the region.
!
  do j = 1, n

    if ( .not. inside(j) ) then
      cycle
    end if

    if ( cos ( point(1,j) ) < point(2,j) ) then
      inside(j) = .false.
      cycle
    end if

    if ( point(2,j) < -5.0D+00 + 5.0D+00 * point(1,j)**4 &
         / ( 2.5D+00 * pi )**4 ) then
      inside(j) = .false.
      cycle
    end if

  end do

  return
end
subroutine p07_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P07_SAMPLE samples points from the region in problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) have
  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ) p(m)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u(m)

  call p07_box ( m, lo, hi )

  have = 0

  do

    call r8vec_uniform_01 ( m, seed, u )

    p(1:2) = ( 1.0D+00 - u(1:2) ) * lo(1:2) + u(1:2) * hi(1:2)

    if ( cos ( p(1) ) < p(2) ) then
      cycle
    end if

    if ( p(2) < -5.0D+00 + 5.0D+00 * p(1)**4 / ( 2.5D+00 * pi )**4 ) then
      cycle
    end if

    have = have + 1
    point(1:2,have) = p(1:2)

    if ( n <= have ) then
      return
    end if

  end do

  return
end
subroutine p07_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P07_SAMPLE_H1 samples points from the enlarged region in problem 07.
!
!  Discussion:
!
!    This is a crude enlargement of the region.  We simply copy the
!    upper and lower bounding curves, and shift them H units up or
!    down. No attempt is made to actually figure out the curve of points
!    which is precisely H units away from the unenlarged region.
!
!    Also, perhaps because the region is about 15 units wide and 6
!    units high, the value of H might need to be scaled (made larger)
!    in order to have the same relative effect as in problems with
!    unit sides.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) h
  real ( kind = 8 ) h2
  integer ( kind = 4 ) have
  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ) p(m)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u(m)

  call p07_box ( m, lo, hi )
!
!  Increase the size of H according to the Y range.
!
  h2 = ( hi(2) - lo(2) ) * h

  lo(2) = lo(2) - h2
  hi(2) = hi(2) + h2

  have = 0

  do

    call r8vec_uniform_01 ( m, seed, u )

    p(1:2) = ( 1.0D+00 - u(1:2) ) * lo(1:2) + u(1:2) * hi(1:2)

    if ( cos ( p(1) ) + h2 < p(2) ) then
      cycle
    end if

    if ( p(2) < -5.0D+00 - h2 + 5.0D+00 * p(1)**4 / ( 2.5D+00 * pi )**4 ) then
      cycle
    end if

    have = have + 1
    point(1:2,have) = p(1:2)

    if ( n <= have ) then
      return
    end if

  end do

  return
end
subroutine p07_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P07_SDIST returns the signed distance to the region in problem 07.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P07_SDIST - Fatal error!'
  write ( *, '(a)' ) '  The routine for this test is not written yet!'
  stop
end
subroutine p07_title ( title )

!*****************************************************************************80
!
!! P07_TITLE returns a title for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = '#7: Bicycle seat (implicit).'

  return
end
subroutine p08_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P08_BOUNDARY_NEAREST returns a nearest boundary point in problem 08.
!
!  Discussion:
!
!    The given input point need not be inside the region.
!
!    In some cases, more than one boundary point may be "nearest",
!    but only one will be returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  real ( kind = 8 ) a
  real ( kind = 8 ), dimension ( m, n ) :: boundary
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.6D+00, 0.0D+00 /)
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_min
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), dimension ( m ) :: pn
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.1D+00
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: theta1_max =  pi / 12.0D+00
  real ( kind = 8 ), parameter :: theta1_min = -pi / 12.0D+00
  real ( kind = 8 ) theta2_max
  real ( kind = 8 ) theta2_min
  real ( kind = 8 ), dimension ( 2, 4 ) :: v1 = reshape ( (/ &
    0.995435619D+00, -0.095435619D+00,   &
    0.9D+00,          0.0D+00,           &
    0.96592581D+00,   0.25881904D+00,    &
    0.0D+00,          0.0D+00 /), (/ 2, 4 /) )
  real ( kind = 8 ), dimension ( 2, 4 ) :: v2 = reshape ( (/ &
    0.9D+00,          0.0D+00,         &
    0.995435619D+00,  0.095435619D+00, &
    0.0D+00,          0.0D+00,         &
    0.96592581D+00,  -0.25881904D+00 /), (/ 2, 4 /) )

  a = ( sqrt ( 119.0D+00 ) - 9.0D+00 ) / 20.0D+00

  theta2_max = atan2 ( a, a + 0.9D+00 )
  theta2_min = -theta2_max

  do j = 1, n

    dist_min = huge ( dist_min )
    boundary(1:m,j) = 0.0D+00
!
!  Compute the distances to the line segments that form four sides
!  of the outer boundary.
!
    do k = 1, 4

      call segment_point_near_2d ( v1(1:2,k), v2(1:2,k), point(1:2,j), &
        pn(1:2), dist, t )

      if ( dist < dist_min ) then
        dist_min = dist
        boundary(1:2,j) = pn(1:2)
      end if

    end do
!
!  Compute the distances to the two arc segments that form two sides
!  of the outer boundary.
!
    call circle_arc_point_near_2d ( r1, center1, theta1_min, theta2_min, &
      point(1:2,j), pn(1:2), dist )

    if ( dist < dist_min ) then
      dist_min = dist
      boundary(1:2,j) = pn(1:2)
    end if

    call circle_arc_point_near_2d ( r1, center1, theta2_max, theta1_max, &
      point(1:2,j), pn(1:2), dist )

   if ( dist < dist_min ) then
      dist_min = dist
      boundary(1:2,j) = pn(1:2)
    end if
!
!  Compute distance to the circle that forms the inner boundary.
!
    call circle_imp_point_near_2d ( r2, center2, point(1:2,j), pn(1:2), &
      dist )

    if ( dist < dist_min ) then
      dist_min = dist
      boundary(1:2,j) = pn(1:2)
    end if

  end do

  return
end
subroutine p08_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P08_BOUNDARY_PROJECT projects exterior points to the boundary in problem 08.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  real ( kind = 8 ) a
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.6D+00, 0.0D+00 /)
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_min
  logical inside(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), dimension ( m ) :: pn
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.1D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) temp(m)
  real ( kind = 8 ), parameter :: theta1_max =  pi / 12.0D+00
  real ( kind = 8 ), parameter :: theta1_min = -pi / 12.0D+00
  real ( kind = 8 ) theta2_max
  real ( kind = 8 ) theta2_min
  real ( kind = 8 ), dimension ( 2, 4 ) :: v1 = reshape ( (/ &
    0.995435619D+00, -0.095435619D+00,   &
    0.9D+00,          0.0D+00,           &
    0.96592581D+00,   0.25881904D+00,    &
    0.0D+00,          0.0D+00 /), (/ 2, 4 /) )
  real ( kind = 8 ), dimension ( 2, 4 ) :: v2 = reshape ( (/ &
    0.9D+00,          0.0D+00,         &
    0.995435619D+00,  0.095435619D+00, &
    0.0D+00,          0.0D+00,         &
    0.96592581D+00,  -0.25881904D+00 /), (/ 2, 4 /) )

  call p08_inside ( m, n, point, inside )

  a = ( sqrt ( 119.0D+00 ) - 9.0D+00 ) / 20.0D+00

  theta2_max = atan2 ( a, a + 0.9D+00 )
  theta2_min = -theta2_max

  do j = 1, n

    if ( inside(j) ) then
      cycle
    end if

    dist_min = huge ( dist_min )
    temp(1:m) = 0.0D+00
!
!  Compute the distances to the line segments that form four sides
!  of the outer boundary.
!
    do k = 1, 4

      call segment_point_near_2d ( v1(1:2,k), v2(1:2,k), point(1:2,j), &
        pn(1:2), dist, t )

      if ( dist < dist_min ) then
        dist_min = dist
        temp(1:m) = pn(1:2)
      end if

    end do
!
!  Compute the distances to the two arc segments that form two sides
!  of the outer boundary.
!
    call circle_arc_point_near_2d ( r1, center1, theta1_min, theta2_min, &
      point(1:2,j), pn(1:2), dist )

    if ( dist < dist_min ) then
      dist_min = dist
      temp(1:m) = pn(1:2)
    end if

    call circle_arc_point_near_2d ( r1, center1, theta2_max, theta1_max, &
      point(1:2,j), pn(1:2), dist )

   if ( dist < dist_min ) then
      dist_min = dist
      temp(1:m) = pn(1:2)
    end if
!
!  Compute distance to the circle that forms the inner boundary.
!
    call circle_imp_point_near_2d ( r2, center2, point(1:2,j), pn(1:2), &
      dist )

    if ( dist < dist_min ) then
      dist_min = dist
      temp(1:m) = pn(1:2)
    end if

    point(1:m,j) = temp(1:m)

  end do

  return
end
subroutine p08_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P08_BOUNDARY_SEGMENT returns a boundary segment in problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  real ( kind = 8 ) a
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.6D+00, 0.0D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.1D+00
  real ( kind = 8 ) s(m)
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) segment(m,segment_length)
  real ( kind = 8 ) t(m)
  real ( kind = 8 ) theta
  real ( kind = 8 ), parameter :: theta1_max =  pi / 12.0D+00
  real ( kind = 8 ), parameter :: theta1_min = -pi / 12.0D+00
  real ( kind = 8 ) theta2_max
  real ( kind = 8 ) theta2_min
!
!  Segment 1: the outer boundary.
!
  if ( segment_index == 1 ) then

    a = ( sqrt ( 119.0D+00 ) - 9.0D+00 ) / 20.0D+00
    theta2_max = atan2 ( a, a + 0.9D+00 )
    theta2_min = -theta2_max
!
!  Work out the appropriate segment lengths, and then
!  adjust N6, if necessary, to account for roundoff.
!
    n1 = nint ( real ( segment_length - 1, kind = 8 ) &
            / ( 2.0D+00 + 2.0D+00 * ( theta1_max - theta2_max ) &
              + 2.0D+00 * a * sqrt ( 2.0D+00 ) ) )
    n2 = nint ( ( theta1_max - theta2_max ) * real ( n1, kind = 8 ) )
    n3 = nint (    a * sqrt ( 2.0D+00 ) * real ( n1, kind = 8 ) )

    n2 = max ( n2, 1 )
    n3 = max ( n3, 1 )

    n4 = n3
    n5 = n2
    n6 = ( segment_length - 1 - n2 - n3 - n4 - n5 ) / 2
    n1 = segment_length - 1 - n2 - n3 - n4 - n5 - n6

    j = 0

    s(1:2) = center1(1:2)

    t(1:2) = (/ center1(1) + r1 * cos ( theta1_min ), &
                center1(2) + r1 * sin ( theta1_min ) /)

    do i = 1, n1
      j = j + 1
      segment(1:2,j) = ( real ( n1 - i + 1, kind = 8 ) * s(1:2)   &
                       + real (      i - 1, kind = 8 ) * t(1:2) ) &
                       / real ( n1,         kind = 8 )
    end do

    do i = 1, n2

      theta = (   real ( n2 - i + 1, kind = 8 ) * theta1_min   &
                + real (      i - 1, kind = 8 ) * theta2_min ) &
              /   real ( n2,         kind = 8 )

      j = j + 1
      segment(1:2,j) = (/ center1(1) + r1 * cos ( theta ), &
                          center1(2) + r1 * sin ( theta ) /)
    end do

    s(1:2) = (/ center1(1) + r1 * cos ( theta2_min ), &
                center1(1) + r1 * sin ( theta2_min ) /)

    t(1:2) = (/ 0.9D+00, 0.0D+00 /)

    do i = 1, n3
      j = j + 1
      segment(1:2,j) = ( real ( n3 - i + 1, kind = 8 ) * s(1:2)   &
                       + real (      i - 1, kind = 8 ) * t(1:2) ) &
                       / real ( n3,         kind = 8 )
    end do

    s(1:2) = (/ 0.9D+00, 0.0D+00 /)

    t(1:2) = (/ center1(1) + r1 * cos ( theta2_max ), &
                center1(1) + r1 * sin ( theta2_max ) /)

    do i = 1, n4
      j = j + 1
      segment(1:2,j) = ( real ( n4 - i + 1, kind = 8 ) * s(1:2)   &
                       + real (      i - 1, kind = 8 ) * t(1:2) ) &
                       / real ( n4,         kind = 8 )
    end do

    do i = 1, n5

      theta = (   real ( n5 - i + 1, kind = 8 ) * theta2_max   &
              +   real (      i - 1, kind = 8 ) * theta1_max ) &
              /   real ( n5,         kind = 8 )

      j = j + 1
      segment(1:2,j) = (/ center1(1) + r1 * cos ( theta ), &
                          center1(1) + r1 * sin ( theta ) /)
    end do

    s(1:2) = (/ center1(1) + r1 * cos ( theta1_max ), &
                center1(1) + r1 * sin ( theta1_max ) /)

    t(1:2) = center1(1:2)

    do i = 1, n6
      j = j + 1
      segment(1:2,j) = ( real ( n6 - i + 1, kind = 8 ) * s(1:2)   &
                       + real (      i - 1, kind = 8 ) * t(1:2) ) &
                       / real ( n6,         kind = 8 )
    end do

    j = j + 1
    segment(1:2,j) = center1(1:2)
!
!  Segment 2: the circular hole.
!
  else if ( segment_index == 2 ) then

    do j = 1, segment_length
      theta = real ( segment_length - j, kind = 8 ) * 2.0D+00 * pi &
            / real ( segment_length - 1, kind = 8 )
      segment(1,j) = center2(1) + r2 * cos ( theta )
      segment(2,j) = center2(2) + r2 * sin ( theta )
    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P08_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p08_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P08_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the
!    boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in
!    the segment.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  real ( kind = 8 ) a
  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.1D+00
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length
  real ( kind = 8 ), parameter :: theta1_max =  pi / 12.0D+00
  real ( kind = 8 ), parameter :: theta1_min = -pi / 12.0D+00
  real ( kind = 8 ) theta2_max
  real ( kind = 8 ) theta2_min

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P08_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    a = ( sqrt ( 119.0D+00 ) - 9.0D+00 ) / 20.0D+00
    theta2_max = atan2 ( a, a + 0.9D+00 )
    theta2_min = -theta2_max

    n = nint ( ( 2.0D+00 + 2.0D+00 * ( theta1_max - theta2_max ) &
        + 2.0D+00 * a * sqrt ( 2.0D+00 ) ) / h )

    n = max ( n, 21 )
    segment_length = n

  else if ( segment_index == 2 ) then

    n = nint ( 2.0D+00 * pi * r2 / h )
    n = max ( n, 5 )
    segment_length = n

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P08_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop
  end if

  return
end
subroutine p08_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P08_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 2

  return
end
subroutine p08_box ( m, lo, hi )

!*****************************************************************************80
!
!! P08_BOX returns a bounding box for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.6D+00, 0.0D+00 /)
  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.1D+00
  real ( kind = 8 ) :: theta1_max =  pi / 12.0D+00
  real ( kind = 8 ) :: theta1_min = -pi / 12.0D+00

  lo(1:m) = (/  center1(1),      r1 * sin ( theta1_min ) /)
  hi(1:m) = (/  center1(1) + r1, r1 * sin ( theta1_max ) /)

  return
end
subroutine p08_density ( m, n, point, density )

!*****************************************************************************80
!
!! P08_DENSITY returns the density for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  real ( kind = 8 ) h3
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  do j = 1, n

    x = point(1,j)
    y = point(2,j)

    h1 = 0.005D+00 + 0.2D+00 *   sqrt (   x**2             + y**2 )
    h2 = 0.02D+00  + 0.2D+00 * ( sqrt ( ( x - 0.6D+00 )**2 + y**2 ) - 1.0D+00 )
    h3 = 0.005D+00 + 0.2D+00 *   sqrt ( ( x - 0.9D+00 )**2 + y**2 )

    density(j) = min ( min ( min ( h1, h2 ), h3 ), 0.03D+00 )

  end do

  return
end
subroutine p08_element_size ( element_size )

!*****************************************************************************80
!
!! P08_ELEMENT_SIZE returns a typical element size for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.005D+00

  return
end
subroutine p08_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P08_FIXED_NUM returns the number of fixed points in problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 6

  return
end
subroutine p08_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P08_FIXED_POINTS returns the fixed points in problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) a
  real ( kind = 8 ) c
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.6D+00, 0.0D+00 /)
  real ( kind = 8 ) fixed(m,fixed_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.1D+00
  real ( kind = 8 ) s
  real ( kind = 8 ) :: theta1_max =  pi / 12.0D+00
  real ( kind = 8 ) :: theta1_min = -pi / 12.0D+00

  a = ( sqrt ( 119.0D+00 ) - 9.0D+00 ) / 20.0D+00

  c = center1(1) + r1 * cos ( pi / 12.0D+00 )
  s = center1(2) + r1 * sin ( pi / 12.0D+00 )

  fixed(1:2,1) = center1(1:2)
  fixed(1:2,2) = (/ c,           -s /)
  fixed(1:2,3) = (/ 0.9D+00 + a, -a /)
  fixed(1:2,4) = (/ 0.9D+00,      0.0D+00 /)
  fixed(1:2,5) = (/ 0.9D+00 + a,  a /)
  fixed(1:2,6) = (/ c,            s /)

  return
end
subroutine p08_header ( )

!*****************************************************************************80
!
!! P08_HEADER prints some information about problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  integer ( kind = 4 ) boundary_segment_num
  real ( kind = 8 ), dimension ( m ) :: center1 = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( m ) :: center2 = (/ 0.6D+00, 0.0D+00 /)
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.1D+00

  call p08_boundary_segment_num ( boundary_segment_num )
  call p08_fixed_num ( fixed_num )
  call p08_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P08:'
  write ( *, '(a)' ) '  Strang and Persson example #8'
  write ( *, '(a)' ) '  Pie slice with notch and hole.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The pie rim is a portion of a circle C1'
  write ( *, '(a,2g14.6)' ) '  with CENTER1 = ', center1(1:2)
  write ( *, '(a,g14.6)' ) '  and radius R1 = ', r1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The interior hole is a circle C2'
  write ( *, '(a,2g14.6)' ) '  with CENTER2 = ', center2(1:2)
  write ( *, '(a,g14.6)' ) '  and radius R2 = ', r2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A uniform mesh density is requested.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p08_hole_num ( hole_num )

!*****************************************************************************80
!
!! P08_HOLE_NUM counts the holes in problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 1

  return
end
subroutine p08_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P08_HOLE_POINT returns a point inside a given hole in problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.6D+00, 0.0D+00 /)
  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  hole_point(1:2) = center2(1:2)

  return
end
subroutine p08_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P08_INSIDE reports if a point is inside the region in problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: nv = 6
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  logical angle_contains_point_2d
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.6D+00, 0.0D+00 /)
  logical inside(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  logical circle_imp_contains_point_2d
  logical circle_sector_contains_point_2d
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: r2 = 0.1D+00
  real ( kind = 8 ), dimension ( 2 ) :: p1 = (/ &
   0.995435619D+00, -0.095435619D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: p2 = (/ &
    0.9000000D+00,  0.0000000D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: p3 = (/ &
   0.995435619D+00, +0.095435619D+00 /)
  real ( kind = 8 ) :: theta1_max =  pi / 12.0D+00
  real ( kind = 8 ) :: theta1_min = -pi / 12.0D+00

  inside(1:n) = .false.

  do j = 1, n
!
!  Is the point inside the circular sector?
!
    if ( circle_sector_contains_point_2d ( r1, center1, theta1_min, &
      theta1_max, point(1:2,j) ) ) then
!
!  Is the point NOT inside the angle?
!
      if ( .not. angle_contains_point_2d ( p3, p2, p1, point(1:2,j) ) ) then
!
!  Is the point NOT inside the circle?
!
        if ( .not. &
          circle_imp_contains_point_2d ( r2, center2, point(1:2,j) ) ) then
          inside(j) = .true.
        end if
      end if
    end if

  end do

  return
end
subroutine p08_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P08_SAMPLE samples points from the region in problem 08.
!
!  Discussion:
!
!    A rejection method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: batch = 1000
  integer ( kind = 4 ) have
  real ( kind = 8 ) hi(m)
  logical, allocatable, dimension ( : ) :: inside
  integer ( kind = 4 ) j
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) reject
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed

  call p08_box ( m, lo, hi )

  have = 0
!
!  We are going to generate batches of sample points.
!
  sample_num = min ( batch, n )

  allocate ( inside (1:sample_num) )
  allocate ( sample(1:m,1:sample_num) )

  reject = 0

  do
!
!  Generate a batch of points in the bounding box.
!
    call r8mat_uniform_01 ( m, sample_num, seed, sample )
!
!  Remap the points to the box.
!
    sample(1,1:sample_num) = lo(1) + sample(1,1:sample_num) * ( hi(1) - lo(1) )
    sample(2,1:sample_num) = lo(2) + sample(2,1:sample_num) * ( hi(2) - lo(2) )

    call p08_inside ( m, sample_num, sample, inside )
!
!  Accept those points which are inside the region.
!
    do j = 1, sample_num

      if ( inside(j) ) then

        have = have + 1
        point(1:m,have) = sample(1:m,j)

        if ( have == n ) then
          return
        end if

      else

        reject = reject + 1

      end if

    end do

    if ( 10.0D+00 + 0.95D+00 * real ( have + reject, kind = 8 ) &
       < real ( reject, kind = 8 ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P08_SAMPLE - Fatal error!'
      write ( *, '(a)' ) '  Too many points rejected!'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Number generated = ', have + reject
      write ( *, '(a,i12)' ) '  Number accepted =  ', have
      write ( *, '(a,i12)' ) '  Number rejected =  ', reject
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Something appears to be wrong!'
      stop
    end if

  end do

  deallocate ( inside )
  deallocate ( sample )

  return
end
subroutine p08_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P08_SAMPLE_H1 samples points from the enlarged region in problem 08.
!
!  Discussion:
!
!    A rejection method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( m ) :: a
  real ( kind = 8 ) a_coef
  logical angle_contains_point_2d
  real ( kind = 8 ), dimension ( m ) :: b
  real ( kind = 8 ) b_coef
  integer ( kind = 4 ), parameter :: batch = 1000
  real ( kind = 8 ), dimension ( m ) :: c
  real ( kind = 8 ) c_coef
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.6D+00, 0.0D+00 /)
  logical circle_imp_contains_point_2d
  logical circle_sector_contains_point_2d
  real ( kind = 8 ), dimension ( m ) :: d
  real ( kind = 8 ), dimension ( m ) :: e
  real ( kind = 8 ), dimension ( m ) :: f
  real ( kind = 8 ) h
  integer ( kind = 4 ) have
  real ( kind = 8 ) hi(m)
  logical indent
  logical, allocatable, dimension ( : ) :: inside
  integer ( kind = 4 ) j
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ), dimension ( 2 ) :: p1 = (/ &
   0.995435619D+00, -0.095435619D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: p2 = (/ &
    0.9000000D+00,  0.0000000D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: p3 = (/ &
   0.995435619D+00, +0.095435619D+00 /)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ) r1h
  real ( kind = 8 ), parameter :: r2 = 0.1D+00
  real ( kind = 8 ) r2h
  integer ( kind = 4 ) reject
  real ( kind = 8 ) root1
  real ( kind = 8 ) root2
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension ( 2, 3 ) :: t
  real ( kind = 8 ) :: theta1_max =  pi / 12.0D+00
  real ( kind = 8 ) :: theta1_maxh
  real ( kind = 8 ) :: theta1_min = -pi / 12.0D+00
  real ( kind = 8 ) :: theta1_minh
  logical triangle_contains_point_2d
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  r1h = r1 + h
  r2h = max ( r2 - h, 0.0D+00 )
!
!  Set point A, to the left of the original angle sector vertex.
!
  a(1:2) = (/ -h / sin ( theta1_max ), 0.0D+00 /)
!
!  Set point B, on the line from A, parallel to the original angle sector
!  side, and striking the circle X^2 + Y^2 = R1H^2.
!
  a_coef = 1.0D+00 + ( sin ( theta1_max ) )**2
  b_coef = 2.0D+00 * h * sin ( theta1_max ) / cos ( theta1_max )
  c_coef = h**2 / ( cos ( theta1_max ) )**2 - r1h**2

  call r8poly2_rroot ( a_coef, b_coef, c_coef, root1, root2 )

  x = max ( root1, root2 )
  y = x * sin ( theta1_max ) + h * cos ( theta1_max )

  b(1:2) = (/ x, y /)
!
!  Set C, the reflection of B.
!
  c(1) =  b(1)
  c(2) = -b(2)
!
!  Set the triangle.
!
  t(1:2,1:3) = reshape ( (/ b(1:2), a(1:2), c(1:2) /), (/ 2, 3 /) )

  if ( h < 0.1D+00 * sqrt ( 2.0D+00 ) ) then
    indent = .true.
  else
    indent = .false.
  end if

  if ( indent ) then
!
!  Set point D, the shifted vertex of the indentation.
!
    d(1:2) = (/ 0.9D+00 + sqrt ( 2.0D+00 ) * h, 0.0D+00 /)
!
!  Set point E, the intersection of the line from D, parallel to the
!  original indentation side with the circle X^2 + Y^2 = R1H^2.
!
    a_coef = 2.0D+00
    b_coef = - 1.8D+00 - 2.0D+00 * sqrt ( 2.0D+00 ) * h
    c_coef = 0.81D+00 + 1.8D+00 * sqrt ( 2.0D+00 ) * h + 2.0D+00 * h**2 - r1h**2

    call r8poly2_rroot ( a_coef, b_coef, c_coef, root1, root2 )

    x = max ( root1, root2 )
    y = x - 0.9D+00

    e(1:2) = (/ x, y /)
!
!  Set point F, the reflection of E.
!
    f(1) =  e(1)
    f(2) = -e(2)

  end if
!
!  Set the angles that will open the circular sector to points
!  B and C.
!
  theta1_maxh = asin ( b(2) )
  theta1_minh = asin ( c(2) )
!
!  Get the ranges, and adjust them.
!
  call p08_box ( m, lo, hi )

  lo(1:2) = lo(1:2) - h
  hi(1:2) = hi(1:2) + h

  have = 0
!
!  We are going to generate batches of sample points.
!
  sample_num = min ( batch, n )

  allocate ( inside (1:sample_num) )
  allocate ( sample(1:m,1:sample_num) )

  reject = 0

  do
!
!  Generate a batch of points in the bounding box.
!
    call r8mat_uniform_01 ( m, sample_num, seed, sample )
!
!  Remap the points to the box.
!
    sample(1,1:sample_num) = lo(1) + sample(1,1:sample_num) * ( hi(1) - lo(1) )
    sample(2,1:sample_num) = lo(2) + sample(2,1:sample_num) * ( hi(2) - lo(2) )

    inside(1:sample_num) = .false.

    do j = 1, sample_num
!
!  Is the point inside the wider, thicker circular sector, or triangle?
!
      if ( &
        circle_sector_contains_point_2d ( r1h, center1, theta1_minh, &
        theta1_maxh, sample(1:2,j) ) &
      .or. &
        triangle_contains_point_2d ( t, sample(1:2,j) ) &
      ) then
!
!  Is the point NOT inside the angle?
!  I should use INDENT to tell me if the angle is an issue or not.
!
        if ( .not. angle_contains_point_2d ( e, d, f, sample(1:2,j) ) ) then
!
!  Is the point NOT inside the inner circle?
!
          if ( .not. &
            circle_imp_contains_point_2d ( r2h, center2, sample(1:2,j) ) ) then
            inside(j) = .true.
          end if
        end if
      end if

    end do
!
!  Accept those points which are inside the region.
!
    do j = 1, sample_num

      if ( inside(j) ) then

        have = have + 1
        point(1:m,have) = sample(1:m,j)

        if ( have == n ) then
          return
        end if

      else

        reject = reject + 1

      end if

    end do

    if ( 10.0D+00 + 0.95D+00 * real ( have + reject, kind = 8 ) &
       < real ( reject, kind = 8 ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P08_SAMPLE_H1 - Fatal error!'
      write ( *, '(a)' ) '  Too many points rejected!'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Number generated = ', have + reject
      write ( *, '(a,i12)' ) '  Number accepted =  ', have
      write ( *, '(a,i12)' ) '  Number rejected =  ', reject
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Something appears to be wrong!'
      stop
    end if

  end do

  deallocate ( inside )
  deallocate ( sample )

  return
end
subroutine p08_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P08_SDIST returns the signed distance to the region in problem 08.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P08_SDIST - Fatal error!'
  write ( *, '(a)' ) '  The routine for this test is not written yet!'

  stop
end
subroutine p08_title ( title )

!*****************************************************************************80
!
!! P08_TITLE returns a title for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = '#8: Pie slice with notch and hole.'

  return
end
subroutine p09_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P09_BOUNDARY_NEAREST returns a nearest boundary point in problem 04.
!
!  Discussion:
!
!    The given input point need not be inside the region.
!
!    In some cases, more than one boundary point may be "nearest",
!    but only one will be returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( m, n ) :: boundary
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.50D+00, 0.50D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.25D+00, 0.75D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center3 = (/ 0.60D+00, 0.40D+00 /)
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) pn(m)
  real ( kind = 8 ) :: r1 = 0.5D+00
  real ( kind = 8 ) :: r2 = 0.1D+00
  real ( kind = 8 ) :: r3 = 0.1D+00
  real ( kind = 8 ) t
  real ( kind = 8 ), dimension(m,4) :: v1
  real ( kind = 8 ), dimension(m,6) :: v2
  real ( kind = 8 ), dimension(m,6) :: v3

  v1(1:2,1:4) = reshape ( (/ &
    center1(1) - r1, center1(2) - r1, &
    center1(1) + r1, center1(2) - r1, &
    center1(1) + r1, center1(2) + r1, &
    center1(1) - r1, center1(2) + r1 /), (/ 2, 4 /) )

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v2(1:2,j) = (/ center2(1) + r2 * cos ( angle ), &
                   center2(2) + r2 * sin ( angle ) /)
  end do

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v3(1:2,j) = (/ center3(1) + r3 * cos ( angle ), &
                   center3(2) + r3 * sin ( angle ) /)
  end do

  do j = 1, n

    dist_min = huge ( dist_min )
    boundary(1:m,j) = 0.0D+00
!
!  Examine points on the outer square.
!
    k = 4

    do i = 1, 4

      call segment_point_near_2d ( v1(1:2,k), v1(1:2,i), point(1:2,j), &
        pn(1:2), dist, t )

      if ( dist < dist_min ) then
        dist_min = dist
        boundary(1:2,j) = pn(1:2)
      end if

      k = i

    end do
!
!  Examine points on the first hexagon.
!
    k = 6

    do i = 1, 6

      call segment_point_near_2d ( v2(1:2,k), v2(1:2,i), point(1:2,j), &
        pn(1:2), dist, t )

      if ( dist < dist_min ) then
        dist_min = dist
        boundary(1:2,j) = pn(1:2)
      end if

      k = i

    end do
!
!  Examine points on the second hexagon.
!
    k = 6

    do i = 1, 6

      call segment_point_near_2d ( v3(1:2,k), v3(1:2,i), point(1:2,j), &
        pn(1:2), dist, t )

      if ( dist < dist_min ) then
        dist_min = dist
        boundary(1:2,j) = pn(1:2)
      end if

      k = i

    end do

  end do

  return
end
subroutine p09_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P09_BOUNDARY_PROJECT projects exterior points to the boundary in problem 09.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.50D+00, 0.50D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.25D+00, 0.75D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center3 = (/ 0.60D+00, 0.40D+00 /)
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_min
  integer ( kind = 4 ) i
  logical inside(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) pn(m)
  real ( kind = 8 ) :: r1 = 0.5D+00
  real ( kind = 8 ) :: r2 = 0.1D+00
  real ( kind = 8 ) :: r3 = 0.1D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) temp(m)
  real ( kind = 8 ), dimension(m,4) :: v1
  real ( kind = 8 ), dimension(m,6) :: v2
  real ( kind = 8 ), dimension(m,6) :: v3

  v1(1:2,1:4) = reshape ( (/ &
    center1(1) - r1, center1(2) - r1, &
    center1(1) + r1, center1(2) - r1, &
    center1(1) + r1, center1(2) + r1, &
    center1(1) - r1, center1(2) + r1 /), (/ 2, 4 /) )

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v2(1:2,j) = (/ center2(1) + r2 * cos ( angle ), &
                   center2(2) + r2 * sin ( angle ) /)
  end do

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v3(1:2,j) = (/ center3(1) + r3 * cos ( angle ), &
                   center3(2) + r3 * sin ( angle ) /)
  end do

  call p09_inside ( m, n, point, inside )

  do j = 1, n

    if ( inside(j) ) then
      cycle
    end if

    dist_min = huge ( dist_min )
    temp(1:m) = 0.0D+00
!
!  Examine points on the outer square.
!
    k = 4

    do i = 1, 4

      call segment_point_near_2d ( v1(1:2,k), v1(1:2,i), point(1:2,j), &
        pn(1:2), dist, t )

      if ( dist < dist_min ) then
        dist_min = dist
        temp(1:m) = pn(1:2)
      end if

      k = i

    end do
!
!  Examine points on the first hexagon.
!
    k = 6

    do i = 1, 6

      call segment_point_near_2d ( v2(1:2,k), v2(1:2,i), point(1:2,j), &
        pn(1:2), dist, t )

      if ( dist < dist_min ) then
        dist_min = dist
        temp(1:m) = pn(1:2)
      end if

      k = i

    end do
!
!  Examine points on the second hexagon.
!
    k = 6

    do i = 1, 6

      call segment_point_near_2d ( v3(1:2,k), v3(1:2,i), point(1:2,j), &
        pn(1:2), dist, t )

      if ( dist < dist_min ) then
        dist_min = dist
        temp(1:m) = pn(1:2)
      end if

      k = i

    end do

    point(1:m,j) = temp(1:m)

  end do

  return
end
subroutine p09_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P09_BOUNDARY_SEGMENT returns a boundary segment in problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.50D+00, 0.50D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.25D+00, 0.75D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center3 = (/ 0.60D+00, 0.40D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) :: r1 = 0.5D+00
  real ( kind = 8 ) :: r2 = 0.1D+00
  real ( kind = 8 ) :: r3 = 0.1D+00
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) segment(m,segment_length)
  real ( kind = 8 ), dimension ( m, 4 ) :: v1
  real ( kind = 8 ), dimension ( m, 6 ) :: v2
  real ( kind = 8 ), dimension ( m, 6 ) :: v3

  v1(1:2,1:4) = reshape ( (/ &
    center1(1) - r1, center1(2) - r1, &
    center1(1) + r1, center1(2) - r1, &
    center1(1) + r1, center1(2) + r1, &
    center1(1) - r1, center1(2) + r1 /), (/ 2, 4 /) )

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v2(1:2,j) = (/ center2(1) + r2 * cos ( angle ), &
                   center2(2) + r2 * sin ( angle ) /)
  end do

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v3(1:2,j) = (/ center3(1) + r3 * cos ( angle ), &
                   center3(2) + r3 * sin ( angle ) /)
  end do

  if ( segment_index == 1 ) then

    n1 = nint ( real ( segment_length - 1, kind = 8 ) &
              / real ( 4, kind = 8 ) )
    n2 = nint ( real ( 2 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 4, kind = 8 ) ) - n1
    n3 = nint ( real ( 3 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 4, kind = 8 ) ) - n1 - n2
    n4 = segment_length - 1 - n1 - n2 - n3

    j = 0

    do i = 1, n1
      j = j + 1
      segment(1:2,j) = ( real ( n1 - i + 1, kind = 8 ) * v1(1:2,1)   &
                       + real (      i - 1, kind = 8 ) * v1(1:2,2) ) &
                       / real ( n1,         kind = 8 )
    end do

    do i = 1, n2
      j = j + 1
      segment(1:2,j) = ( real ( n2 - i + 1, kind = 8 ) * v1(1:2,2)   &
                       + real (      i - 1, kind = 8 ) * v1(1:2,3) ) &
                       / real ( n2,         kind = 8 )
    end do

    do i = 1, n3
      j = j + 1
      segment(1:2,j) = ( real ( n3 - i + 1, kind = 8 ) * v1(1:2,3)   &
                       + real (      i - 1, kind = 8 ) * v1(1:2,4) ) &
                       / real ( n3,         kind = 8 )
    end do

    do i = 1, n4
      j = j + 1
      segment(1:2,j) = ( real ( n4 - i + 1, kind = 8 ) * v1(1:2,4)   &
                       + real (      i - 1, kind = 8 ) * v1(1:2,1) ) &
                       / real ( n4,         kind = 8 )
    end do

    j = j + 1
    segment(1:2,j) = v1(1:2,1)

  else if ( segment_index == 2 ) then

    n1 = nint ( real ( segment_length - 1, kind = 8 ) &
              / real ( 6, kind = 8 ) )
    n2 = nint ( real ( 2 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1
    n3 = nint ( real ( 3 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1 - n2
    n4 = nint ( real ( 4 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1 - n2 - n3
    n5 = nint ( real ( 5 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1 - n2 - n3 - n4
    n6 = segment_length - 1 - n1 - n2 - n3 - n4 - n5

    j = 0

    do i = 1, n1
      j = j + 1
      segment(1:2,j) = ( real ( n1 - i + 1, kind = 8 ) * v2(1:2,1)   &
                       + real (      i - 1, kind = 8 ) * v2(1:2,2) ) &
                       / real ( n1,         kind = 8 )
    end do

    do i = 1, n2
      j = j + 1
      segment(1:2,j) = ( real ( n2 - i + 1, kind = 8 ) * v2(1:2,2)   &
                       + real (      i - 1, kind = 8 ) * v2(1:2,3) ) &
                       / real ( n2,         kind = 8 )
    end do

    do i = 1, n3
      j = j + 1
      segment(1:2,j) = ( real ( n3 - i + 1, kind = 8 ) * v2(1:2,3)   &
                       + real (      i - 1, kind = 8 ) * v2(1:2,4) ) &
                       / real ( n3,         kind = 8 )
    end do

    do i = 1, n4
      j = j + 1
      segment(1:2,j) = ( real ( n4 - i + 1, kind = 8 ) * v2(1:2,4)   &
                       + real (      i - 1, kind = 8 ) * v2(1:2,5) ) &
                       / real ( n4,         kind = 8 )
    end do

    do i = 1, n5
      j = j + 1
      segment(1:2,j) = ( real ( n5 - i + 1, kind = 8 ) * v2(1:2,5)   &
                       + real (      i - 1, kind = 8 ) * v2(1:2,6) ) &
                       / real ( n5,         kind = 8 )
    end do

    do i = 1, n6
      j = j + 1
      segment(1:2,j) = ( real ( n6 - i + 1, kind = 8 ) * v2(1:2,6)   &
                       + real (      i - 1, kind = 8 ) * v2(1:2,1) ) &
                       / real ( n6,         kind = 8 )
    end do

    j = j + 1
    segment(1:2,j) = v2(1:2,1)

  else if ( segment_index == 3 ) then

    n1 = nint ( real ( segment_length - 1, kind = 8 ) &
              / real ( 6, kind = 8 ) )
    n2 = nint ( real ( 2 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1
    n3 = nint ( real ( 3 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1 - n2
    n4 = nint ( real ( 4 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1 - n2 - n3
    n5 = nint ( real ( 5 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 6, kind = 8 ) ) - n1 - n2 - n3 - n4
    n6 = segment_length - 1 - n1 - n2 - n3 - n4 - n5

    j = 0

    do i = 1, n1
      j = j + 1
      segment(1:2,j) = ( real ( n1 - i + 1, kind = 8 ) * v3(1:2,1)   &
                       + real (      i - 1, kind = 8 ) * v3(1:2,2) ) &
                       / real ( n1,         kind = 8 )
    end do

    do i = 1, n2
      j = j + 1
      segment(1:2,j) = ( real ( n2 - i + 1, kind = 8 ) * v3(1:2,2)   &
                       + real (      i - 1, kind = 8 ) * v3(1:2,3) ) &
                       / real ( n2,         kind = 8 )
    end do

    do i = 1, n3
      j = j + 1
      segment(1:2,j) = ( real ( n3 - i + 1, kind = 8 ) * v3(1:2,3)   &
                       + real (      i - 1, kind = 8 ) * v3(1:2,4) ) &
                       / real ( n3,         kind = 8 )
    end do

    do i = 1, n4
      j = j + 1
      segment(1:2,j) = ( real ( n4 - i + 1, kind = 8 ) * v3(1:2,4)   &
                       + real (      i - 1, kind = 8 ) * v3(1:2,5) ) &
                       / real ( n4,         kind = 8 )
    end do

    do i = 1, n5
      j = j + 1
      segment(1:2,j) = ( real ( n5 - i + 1, kind = 8 ) * v3(1:2,5)   &
                       + real (      i - 1, kind = 8 ) * v3(1:2,6) ) &
                       / real ( n5,         kind = 8 )
    end do

    do i = 1, n6
      j = j + 1
      segment(1:2,j) = ( real ( n6 - i + 1, kind = 8 ) * v3(1:2,6)   &
                       + real (      i - 1, kind = 8 ) * v3(1:2,1) ) &
                       / real ( n6,         kind = 8 )
    end do

    j = j + 1
    segment(1:2,j) = v3(1:2,1)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P09_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p09_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P09_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P09_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    n = nint ( 4.0D+00 / h )
    n = max ( n, 5 )
    segment_length = n + mod ( 4 - mod ( n - 1, 4 ), 4 )

  else if ( segment_index == 2 ) then

    n = nint ( 0.6D+00 / h )
    n = max ( n, 7 )
    segment_length = n + mod ( 6 - mod ( n - 1, 6 ), 6 )

  else if ( segment_index == 3 ) then

    n = nint ( 0.6D+00 / h )
    n = max ( n, 7 )
    segment_length = n + mod ( 6 - mod ( n - 1, 6 ), 6 )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P09_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p09_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P09_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 3

  return
end
subroutine p09_box ( m, lo, hi )

!*****************************************************************************80
!
!! P09_BOX returns a bounding box for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.50D+00, 0.50D+00 /)
  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)
  real ( kind = 8 ) :: r1 = 0.5D+00

  lo(1:m) = (/ center1(1) - r1, center1(2) - r1 /)
  hi(1:m) = (/ center1(1) + r1, center1(2) + r1 /)

  return
end
subroutine p09_density ( m, n, point, density )

!*****************************************************************************80
!
!! P09_DENSITY returns the density for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p09_element_size ( element_size )

!*****************************************************************************80
!
!! P09_ELEMENT_SIZE returns a typical element size for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.005D+00

  return
end
subroutine p09_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P09_FIXED_NUM returns the number of fixed points in problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 16

  return
end
subroutine p09_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P09_FIXED_POINTS returns the fixed points in problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.50D+00, 0.50D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.25D+00, 0.75D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center3 = (/ 0.60D+00, 0.40D+00 /)
  real ( kind = 8 ), dimension(m,fixed_num) :: fixed
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) :: r1 = 0.5D+00
  real ( kind = 8 ) :: r2 = 0.1D+00
  real ( kind = 8 ) :: r3 = 0.1D+00
  real ( kind = 8 ), dimension ( m, 4 ) :: v1
  real ( kind = 8 ), dimension ( m, 6 ) :: v2
  real ( kind = 8 ), dimension ( m, 6 ) :: v3

  v1(1:2,1:4) = reshape ( (/ &
    center1(1) - r1, center1(2) - r1, &
    center1(1) + r1, center1(2) - r1, &
    center1(1) + r1, center1(2) + r1, &
    center1(1) - r1, center1(2) + r1 /), (/ 2, 4 /) )

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v2(1:2,j) = (/ center2(1) + r2 * cos ( angle ), &
                   center2(2) + r2 * sin ( angle ) /)
  end do

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v3(1:2,j) = (/ center3(1) + r3 * cos ( angle ), &
                   center3(2) + r3 * sin ( angle ) /)
  end do

  fixed(1:2,1:4) = v1(1:2,1:4)

  fixed(1:2,5:10) = v2(1:2,1:6)

  fixed(1:2,11:16) = v3(1:2,1:6)

  return
end
subroutine p09_header ( )

!*****************************************************************************80
!
!! P09_HEADER prints some information about problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  integer ( kind = 4 ) boundary_segment_num
  real ( kind = 8 ), dimension ( m ) :: center1 = (/ 0.50D+00, 0.50D+00 /)
  real ( kind = 8 ), dimension ( m ) :: center2 = (/ 0.25D+00, 0.75D+00 /)
  real ( kind = 8 ), dimension ( m ) :: center3 = (/ 0.60D+00, 0.40D+00 /)
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num
  real ( kind = 8 ) :: r1 = 0.5D+00
  real ( kind = 8 ) :: r2 = 0.1D+00
  real ( kind = 8 ) :: r3 = 0.1D+00

  call p09_boundary_segment_num ( boundary_segment_num )
  call p09_fixed_num ( fixed_num )
  call p09_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P09:'
  write ( *, '(a)' ) '  Jeff Borggaard''s example'
  write ( *, '(a)' ) '  A square with 2 hexagonal holes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  The square has "center" at ', center1(1:2)
  write ( *, '(a,g14.6)' ) '  and "radius" R1 = ', r1
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  Hexagon 1 has "center" at ', center2(1:2)
  write ( *, '(a,g14.6)' ) '  and "radius" R2 = ', r2
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  Hexagon 2 has "center" at ', center3(1:2)
  write ( *, '(a,g14.6)' ) '  and "radius" R3 = ', r3
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A uniform mesh density is requested.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p09_hole_num ( hole_num )

!*****************************************************************************80
!
!! P09_HOLE_NUM counts the holes in problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 2

  return
end
subroutine p09_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P09_HOLE_POINT returns a point inside a given hole in problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.25D+00, 0.75D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center3 = (/ 0.60D+00, 0.40D+00 /)
  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  if ( hole_index == 1 ) then
    hole_point(1:2) = center2(1:2)
  else if ( hole_index == 2 ) then
    hole_point(1:2) = center3(1:2)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P09_HOLE_POINT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of HOLE_INDEX = ', hole_index
    write ( *, '(a)' ) '  Legal values are 1 or 2.'
    stop
  end if

  return
end
subroutine p09_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P09_INSIDE reports if a point is inside the region in problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.50D+00, 0.50D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.25D+00, 0.75D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center3 = (/ 0.60D+00, 0.40D+00 /)
  logical hexagon_contains_point_2d
  logical inside(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) :: r1 = 0.5D+00
  real ( kind = 8 ) :: r2 = 0.1D+00
  real ( kind = 8 ) :: r3 = 0.1D+00
  real ( kind = 8 ), dimension ( m, 4 ) :: v1
  real ( kind = 8 ), dimension ( m, 6 ) :: v2
  real ( kind = 8 ), dimension ( m, 6 ) :: v3

  v1(1:2,1:4) = reshape ( (/ &
    center1(1) - r1, center1(2) - r1, &
    center1(1) + r1, center1(2) - r1, &
    center1(1) + r1, center1(2) + r1, &
    center1(1) - r1, center1(2) + r1 /), (/ 2, 4 /) )

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v2(1:2,j) = (/ center2(1) + r2 * cos ( angle ), &
                   center2(2) + r2 * sin ( angle ) /)
  end do

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v3(1:2,j) = (/ center3(1) + r3 * cos ( angle ), &
                   center3(2) + r3 * sin ( angle ) /)
  end do

  inside(1:n) = .true.

  do j = 1, n

    if ( point(1,j) < center1(1) - r1 .or. &
         center1(1) + r1 < point(1,j) .or. &
         point(2,j) < center1(2) - r1 .or. &
         center1(2) + r1 < point(2,j) ) then

      inside(j) = .false.

    else if ( hexagon_contains_point_2d ( v2, point(1:2,j) ) ) then

      inside(j) = .false.

    else if ( hexagon_contains_point_2d ( v3, point(1:2,j) ) ) then

      inside(j) = .false.

    end if

  end do

  return
end
subroutine p09_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P09_SAMPLE samples points from the region in problem 09.
!
!  Discussion:
!
!    A rejection method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.50D+00, 0.50D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.25D+00, 0.75D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center3 = (/ 0.60D+00, 0.40D+00 /)
  logical hexagon_contains_point_2d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) reject
  integer ( kind = 4 ) seed
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) :: r1 = 0.5D+00
  real ( kind = 8 ) :: r2 = 0.1D+00
  real ( kind = 8 ) :: r3 = 0.1D+00
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ), dimension ( 2, 4 ) :: v1
  real ( kind = 8 ), dimension ( 2, 6 ) :: v2
  real ( kind = 8 ), dimension ( 2, 6 ) :: v3
  real ( kind = 8 ) y(m)

  v1(1:2,1:4) = reshape ( (/ &
    center1(1) - r1, center1(2) - r1, &
    center1(1) + r1, center1(2) - r1, &
    center1(1) + r1, center1(2) + r1, &
    center1(1) - r1, center1(2) + r1 /), (/ 2, 4 /) )

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v2(1:2,j) = (/ center2(1) + r2 * cos ( angle ), &
                   center2(2) + r2 * sin ( angle ) /)
  end do

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v3(1:2,j) = (/ center3(1) + r3 * cos ( angle ), &
                   center3(2) + r3 * sin ( angle ) /)
  end do

  reject = 0

  do j = 1, n

    do

      do i = 1, m
        y(i) = r8_uniform_01 ( seed )
      end do

      if ( &
        ( .not. hexagon_contains_point_2d ( v2, y ) ) &
      .and. &
        ( .not. hexagon_contains_point_2d ( v3, y ) ) &
      ) then
        exit
      end if

      reject = reject + 1

      if ( 2 * n + 10 <= reject ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P09_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  (The double hexagonal hole region)'
        write ( *, '(a,i8)' ) '  Trying to generate point J = ', j
        write ( *, '(a,i8)' ) '  Number of rejections = ', reject
        write ( *, '(a,g14.6)' ) '  Rejection percentage = ', &
          real ( 100 * reject ) / real ( reject + j - 1 )
        write ( *, '(a,2g14.6)' ) '  Y = ', y(1:2)
        stop
      end if

    end do

    point(1:m,j) = y(1:m)

  end do

  return
end
subroutine p09_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P09_SAMPLE_H1 samples points from the enlarged region in problem 09.
!
!  Discussion:
!
!    A rejection method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( 2 ) :: center1 = (/ 0.50D+00, 0.50D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center2 = (/ 0.25D+00, 0.75D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: center3 = (/ 0.60D+00, 0.40D+00 /)
  real ( kind = 8 ) h
  logical hexagon_contains_point_2d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) reject
  integer ( kind = 4 ) seed
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) :: r1 = 0.5D+00
  real ( kind = 8 ) r1h
  real ( kind = 8 ) :: r2 = 0.1D+00
  real ( kind = 8 ) r2h
  real ( kind = 8 ) :: r3 = 0.1D+00
  real ( kind = 8 ) r3h
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ), dimension ( 2, 4 ) :: v1
  real ( kind = 8 ), dimension ( 2, 6 ) :: v2
  real ( kind = 8 ), dimension ( 2, 6 ) :: v3
  real ( kind = 8 ) y(m)

  r1h = r1 + h
  r2h = max ( r2 - h, 0.0D+00 )
  r3h = max ( r3 - h, 0.0D+00 )

  v1(1:2,1:4) = reshape ( (/ &
    center1(1) - r1h, center1(2) - r1h, &
    center1(1) + r1h, center1(2) - r1h, &
    center1(1) + r1h, center1(2) + r1h, &
    center1(1) - r1h, center1(2) + r1h /), (/ 2, 4 /) )

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v2(1:2,j) = (/ center2(1) + r2h * cos ( angle ), &
                   center2(2) + r2h * sin ( angle ) /)
  end do

  do j = 1, 6
    angle = real ( ( j - 1 ) * 2, kind = 8 ) * pi / 6.0D+00
    v3(1:2,j) = (/ center3(1) + r3h * cos ( angle ), &
                   center3(2) + r3h * sin ( angle ) /)
  end do

  reject = 0

  do j = 1, n

    do

      do i = 1, m
        y(i) = center1(i) &
          + ( 2.0D+00 * r8_uniform_01 ( seed ) - 1.0D+00 ) * r1h
      end do

      if ( &
        ( .not. hexagon_contains_point_2d ( v2, y ) ) &
      .and. &
        ( .not. hexagon_contains_point_2d ( v3, y ) ) &
      ) then
        exit
      end if

      reject = reject + 1

      if ( 2 * n + 10 <= reject ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P09_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  (The double hexagonal hole region)'
        write ( *, '(a,i8)' ) '  Trying to generate point J = ', j
        write ( *, '(a,i8)' ) '  Number of rejections = ', reject
        write ( *, '(a,g14.6)' ) '  Rejection percentage = ', &
          real ( 100 * reject ) / real ( reject + j - 1 )
        write ( *, '(a,2g14.6)' ) '  Y = ', y(1:2)
        stop
      end if

    end do

    point(1:m,j) = y(1:m)

  end do

  return
end
subroutine p09_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P09_SDIST returns the signed distance to the region in problem 09.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, June 2004, pages 329-345.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P09_SDIST - Fatal error!'
  write ( *, '(a)' ) '  The routine for this test is not written yet!'

  stop
end
subroutine p09_title ( title )

!*****************************************************************************80
!
!! P09_TITLE returns a title for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
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

  title = '#9: Jeff Borggaard''s Box with 2 hexagonal holes.'

  return
end
subroutine p10_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P10_BOUNDARY_NEAREST returns a nearest boundary point in problem 10.
!
!  Discussion:
!
!    The given input point need not be inside the region.
!
!    In some cases, more than one boundary point may be "nearest",
!    but only one will be returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( m, n ) :: boundary
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) temp
  real ( kind = 8 ), parameter :: x1 =  0.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ), parameter :: y1 =  0.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00

  do j = 1, n

    if ( point(1,j) <= x1 ) then

      if ( point(2,j) <= y1 ) then
        boundary(1,j) = x1
        boundary(2,j) = y1
      else if ( point(2,j) <= y2 ) then
        boundary(1,j) = x1
        boundary(2,j) = point(2,j)
      else if ( y2 < point(2,j) ) then
        boundary(1,j) = x1
        boundary(2,j) = y2
      end if

    else if ( point(1,j) <= x2 ) then

      if ( point(2,j) <= y1 ) then
        boundary(1,j) = point(1,j)
        boundary(2,j) = y1
      else if ( point(2,j) <= y2 ) then

        temp = min (            &
               point(1,j) - x1, &
          x2 - point(1,j),      &
               point(2,j) - y1, &
          y2 - point(2,j) )

        if ( point(1,j) - x1 == temp ) then
          boundary(1,j) = x1
          boundary(2,j) = point(2,j)
        else if ( point(2,j) - y1 == temp  ) then
          boundary(1,j) = point(1,j)
          boundary(2,j) = y1
        else if ( x2 - point(1,j) == temp  ) then
          boundary(1,j) = x2
          boundary(2,j) = point(2,j)
        else if ( y2 - point(1,j) == temp  ) then
          boundary(1,j) = point(1,j)
          boundary(2,j) = y2
        end if

      else if ( y2 < point(2,j) ) then
        boundary(1,j) = point(1,j)
        boundary(2,j) = y2
      end if

    else if ( x2 < point(1,j) ) then

      if ( point(2,j) <= y1 ) then
        boundary(1,j) = x2
        boundary(2,j) = y1
      else if ( point(2,j) <= y2 ) then
        boundary(1,j) = x2
        boundary(2,j) = point(2,j)
      else if ( y2 < point(2,j) ) then
        boundary(1,j) = x2
        boundary(2,j) = y2
      end if

    end if

  end do

  return
end
subroutine p10_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P10_BOUNDARY_PROJECT projects exterior points to the boundary in problem 10.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical inside(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) temp
  real ( kind = 8 ), parameter :: x1 =  0.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ), parameter :: y1 =  0.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00

  call p10_inside ( m, n, point, inside )

  do j = 1, n

    if ( inside(j) ) then
      cycle
    end if

    if ( point(1,j) <= x1 ) then

      if ( point(2,j) <= y1 ) then
        point(1,j) = x1
        point(2,j) = y1
      else if ( point(2,j) <= y2 ) then
        point(1,j) = x1
      else if ( y2 < point(2,j) ) then
        point(1,j) = x1
        point(2,j) = y2
      end if

    else if ( point(1,j) <= x2 ) then

      if ( point(2,j) <= y1 ) then
        point(2,j) = y1
      else if ( point(2,j) <= y2 ) then

        temp = min (            &
               point(1,j) - x1, &
          x2 - point(1,j),      &
               point(2,j) - y1, &
          y2 - point(2,j) )

        if ( point(1,j) - x1 == temp ) then
          point(1,j) = x1
        else if ( point(2,j) - y1 == temp  ) then
          point(2,j) = y1
        else if ( x2 - point(1,j) == temp  ) then
          point(1,j) = x2
        else if ( y2 - point(1,j) == temp  ) then
          point(2,j) = y2
        end if

      else if ( y2 < point(2,j) ) then
        point(2,j) = y2
      end if

    else if ( x2 < point(1,j) ) then

      if ( point(2,j) <= y1 ) then
        point(1,j) = x2
        point(2,j) = y1
      else if ( point(2,j) <= y2 ) then
        point(1,j) = x2
      else if ( y2 < point(2,j) ) then
        point(1,j) = x2
        point(2,j) = y2
      end if

    end if

  end do

  return
end
subroutine p10_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P10_BOUNDARY_SEGMENT returns a boundary segment in problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  real ( kind = 8 ) s(m,4)
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) segment(m,segment_length)

  if ( segment_index == 1 ) then

    n1 = nint ( real ( segment_length - 1, kind = 8 ) &
              / real ( 4, kind = 8 ) )
    n2 = nint ( real ( 2 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 4, kind = 8 ) ) - n1
    n3 = nint ( real ( 3 * ( segment_length - 1 ), kind = 8 ) &
              / real ( 4, kind = 8 ) ) - n1 - n2
    n4 = segment_length - 1 - n1 - n2 - n3

    s(1:2,1) = (/  0.0D+00,  0.0D+00 /)
    s(1:2,2) = (/  1.0D+00,  0.0D+00 /)
    s(1:2,3) = (/  1.0D+00,  1.0D+00 /)
    s(1:2,4) = (/  0.0D+00,  1.0D+00 /)

    j = 0

    do i = 1, n1
      j = j + 1
      segment(1:2,j) = ( real ( n1 - i + 1, kind = 8 ) * s(1:2,1)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,2) ) &
                       / real ( n1,         kind = 8 )
    end do

    do i = 1, n2
      j = j + 1
      segment(1:2,j) = ( real ( n2 - i + 1, kind = 8 ) * s(1:2,2)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,3) ) &
                       / real ( n2,         kind = 8 )
    end do

    do i = 1, n3
      j = j + 1
      segment(1:2,j) = ( real ( n3 - i + 1, kind = 8 ) * s(1:2,3)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,4) ) &
                       / real ( n3,         kind = 8 )
    end do

    do i = 1, n4
      j = j + 1
      segment(1:2,j) = ( real ( n4 - i + 1, kind = 8 ) * s(1:2,4)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,1) ) &
                       / real ( n4,         kind = 8 )
    end do

    j = j + 1
    segment(1:2,j) = s(1:2,1)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P10_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p10_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P10_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P10_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    n = nint ( 4.0D+00 / h )
    n = max ( n, 5 )
    segment_length = n + mod ( 4 - mod ( n - 1, 4 ), 4 )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P10_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p10_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P10_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 1

  return
end
subroutine p10_box ( m, lo, hi )

!*****************************************************************************80
!
!! P10_BOX returns a bounding box for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)

  lo(1:m) = (/  0.0D+00,  0.0D+00 /)
  hi(1:m) = (/ +1.0D+00, +1.0D+00 /)

  return
end
subroutine p10_density ( m, n, point, density )

!*****************************************************************************80
!
!! P10_DENSITY returns the density for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p10_element_size ( element_size )

!*****************************************************************************80
!
!! P10_ELEMENT_SIZE returns a typical element size for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.15D+00

  return
end
subroutine p10_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P10_FIXED_NUM returns the number of fixed points in problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 4

  return
end
subroutine p10_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P10_FIXED_POINTS returns the fixed points in problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) fixed(m,fixed_num)

  fixed(1:2,1) = (/  0.0D+00,  0.0D+00 /)
  fixed(1:2,2) = (/  1.0D+00,  0.0D+00 /)
  fixed(1:2,3) = (/  1.0D+00,  1.0D+00 /)
  fixed(1:2,4) = (/  0.0D+00,  1.0D+00 /)

  return
end
subroutine p10_header ( )

!*****************************************************************************80
!
!! P10_HEADER prints some information about problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num

  call p10_boundary_segment_num ( boundary_segment_num )
  call p10_fixed_num ( fixed_num )
  call p10_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P10:'
  write ( *, '(a)' ) '  The unit square.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p10_hole_num ( hole_num )

!*****************************************************************************80
!
!! P10_HOLE_NUM counts the holes in problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 0

  return
end
subroutine p10_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P10_HOLE_POINT returns a point inside a given hole in problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  return
end
subroutine p10_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P10_INSIDE reports if a point is inside the region in problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical inside(n)
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), parameter :: x1 =  0.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ), parameter :: y1 =  0.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00

  inside(1:n) =                         &
      x1 <=    point(1,1:n)       .and. &
               point(1,1:n) <= x2 .and. &
      y1 <=    point(2,1:n)       .and. &
               point(2,1:n) <= y2

  return
end
subroutine p10_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P10_SAMPLE samples points from the region in problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), parameter :: x1 =  0.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ), parameter :: y1 =  0.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00
!
!  Generate a batch of points in [0,1]x[0,1].
!
  call r8mat_uniform_01 ( m, n, seed, point )
!
!  Remap the points to the box [X1,X2] x [Y1,Y2].
!
  point(1,1:n) = x1 + point(1,1:n) * ( x2 - x1 )
  point(2,1:n) = y1 + point(2,1:n) * ( y2 - y1 )

  return
end
subroutine p10_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P10_SAMPLE_H1 samples points from the enlarged region in problem 10.
!
!  Discussion:
!
!    The region is enlarged by an amount H in each direction.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement amount.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) h
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), parameter :: x1 =  0.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ), parameter :: y1 =  0.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00
!
!  Generate a batch of points in [0,1]x[0,1].
!
  call r8mat_uniform_01 ( m, n, seed, point )
!
!  Remap the points to the box [X1,X2] x [Y1,Y2].
!
  point(1,1:n) = ( x1 - h ) + point(1,1:n) * ( x2 - x1 + 2.0D+00 * h )
  point(2,1:n) = ( y1 - h ) + point(2,1:n) * ( y2 - y1 + 2.0D+00 * h )

  return
end
subroutine p10_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P10_SDIST returns the signed distance to the region in problem 10.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)
  real ( kind = 8 ), parameter :: x1 =  0.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ), parameter :: y1 =  0.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00

  do j = 1, n

    if ( point(1,j) <= x1 ) then

      if ( point(2,j) <= y1 ) then
        sdist = sqrt ( ( point(1,j) - x1 )**2 + ( point(2,j) - y1 )**2 )
      else if ( point(2,j) <= y2 ) then
        sdist = x1 - point(1,j)
      else if ( y2 < point(2,j) ) then
        sdist = sqrt ( ( point(1,j) - x1 )**2 + ( point(2,j) - y2 )**2 )
      end if

    else if ( point(1,j) <= x2 ) then

      if ( point(2,j) <= y1 ) then
        sdist = y1 - point(2,j)
      else if ( point(2,j) <= y2 ) then
        sdist = - min (           &
               point(1,j) - x1, &
          x2 - point(1,j),      &
               point(2,j) - y1, &
          y2 - point(2,j) )
      else if ( y2 < point(2,j) ) then
        sdist = point(2,j) - y2
      end if

    else if ( x2 < point(1,j) ) then

      if ( point(2,j) <= y1 ) then
        sdist = sqrt ( ( point(1,j) - x2 )**2 + ( point(2,j) - y1 )**2 )
      else if ( point(2,j) <= y2 ) then
        sdist = point(1,j) - x2
      else if ( y2 < point(2,j) ) then
        sdist = sqrt ( ( point(1,j) - x2 )**2 + ( point(2,j) - y2 )**2 )
      end if

    end if

  end do

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
!    18 November 2004
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

  title = '#10: The unit square.'

  return
end
subroutine p11_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P11_BOUNDARY_NEAREST returns a nearest boundary point in problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( m, n ) :: boundary
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ), dimension(2,4) :: q1 = reshape ( (/ &
    0.00D+00, 0.00D+00, &
    1.00D+00, 0.00D+00, &
    0.75D+00, 0.25D+00, &
    0.25D+00, 0.25D+00 /), (/ 2, 4 /) )
  real ( kind = 8 ), dimension(2,4) :: q2 = reshape ( (/ &
    0.00D+00, 0.00D+00, &
    0.25D+00, 0.25D+00, &
    0.25D+00, 0.75D+00, &
    0.00D+00, 1.00D+00 /), (/ 2, 4 /) )
  real ( kind = 8 ), dimension(2,4) :: q3 = reshape ( (/ &
    0.50D+00, 0.50D+00, &
    0.50D+00, 0.25D+00, &
    0.75D+00, 0.25D+00, &
    1.00D+00, 0.50D+00 /), (/ 2, 4 /) )
  real ( kind = 8 ), dimension(2,4) :: q4 = reshape ( (/ &
    0.50D+00, 0.50D+00, &
    0.25D+00, 0.50D+00, &
    0.25D+00, 0.25D+00, &
    0.50D+00, 0.25D+00 /), (/ 2, 4 /) )
  real ( kind = 8 ), dimension(2,4) :: q5 = reshape ( (/ &
    0.50D+00, 0.50D+00, &
    0.50D+00, 1.00D+00, &
    0.25D+00, 0.75D+00, &
    0.25D+00, 0.50D+00 /), (/ 2, 4 /) )
  logical quad_contains_point_2d
  real ( kind = 8 ), dimension(2,3) :: t1 = reshape ( (/ &
    1.00D+00, 0.00D+00, &
    1.00D+00, 0.50D+00, &
    0.75D+00, 0.25D+00 /), (/ 2, 3 /) )
  real ( kind = 8 ), dimension(2,3) :: t2 = reshape ( (/ &
    0.00D+00, 1.00D+00, &
    0.25D+00, 0.25D+00, &
    0.50D+00, 1.00D+00 /), (/ 2, 3 /) )
  real ( kind = 8 ), dimension(2,3) :: t3 = reshape ( (/ &
    0.50D+00, 0.50D+00, &
    1.00D+00, 0.50D+00, &
    1.00D+00, 1.00D+00 /), (/ 2, 3 /) )
  real ( kind = 8 ), dimension(2,3) :: t4 = reshape ( (/ &
    0.50D+00, 0.50D+00, &
    1.00D+00, 1.00D+00, &
    0.50D+00, 1.00D+00 /), (/ 2, 3 /) )
  logical triangle_contains_point_2d
  real ( kind = 8 ), parameter :: x1 =  0.0D+00
  real ( kind = 8 ), parameter :: x2 =  0.5D+00
  real ( kind = 8 ), parameter :: x3 = +1.0D+00
  real ( kind = 8 ), parameter :: y1 =  0.0D+00
  real ( kind = 8 ), parameter :: y2 =  0.5D+00
  real ( kind = 8 ), parameter :: y3 = +1.0D+00

  do j = 1, n

    if ( point(1,j) <= x1 .and. point(2,j) <= y1 ) then

      boundary(1,j) = x1
      boundary(2,j) = y1

    else if ( point(1,j) <= x1 .and. point(2,j) <= y3 ) then

      boundary(1,j) = x1
      boundary(2,j) = point(2,j)

    else if ( point(1,j) <= x1 .and. y3 <= point(2,j) ) then

      boundary(1,j) = x1
      boundary(2,j) = y3

    else if ( point(1,j) <= x3 .and. point(2,j) <= y1 ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y1

    else if ( quad_contains_point_2d ( q1, point(1:2,j) ) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y1

    else if ( quad_contains_point_2d ( q2, point(1:2,j) ) ) then

      boundary(1,j) = x1
      boundary(2,j) = point(2,j)

    else if ( triangle_contains_point_2d ( t1, point(1:2,j) ) ) then

      boundary(1,j) = x3
      boundary(2,j) = point(2,j)

    else if ( quad_contains_point_2d ( q3, point(1:2,j) ) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y2

    else if ( quad_contains_point_2d ( q4, point(1:2,j) ) ) then

      boundary(1,j) = x2
      boundary(2,j) = y2

    else if ( quad_contains_point_2d ( q5, point(1:2,j) ) ) then

      boundary(1,j) = x2
      boundary(2,j) = point(2,j)

    else if ( triangle_contains_point_2d ( t2, point(1:2,j) ) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y3

    else if ( x1 <= point(1,j) .and. point(1,j) <= x2 .and. &
      y3 <= point(2,j) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y3

    else if ( triangle_contains_point_2d ( t3, point(1:2,j) ) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y2

    else if ( triangle_contains_point_2d ( t4, point(1:2,j) ) ) then

      boundary(1,j) = x2
      boundary(2,j) = point(2,j)

    else if ( x2 <= point(1,j) .and. point(1,j) <= x3 .and. &
      y3 <= point(2,j) ) then

      boundary(1,j) = x2
      boundary(2,j) = y3

    else if ( x3 <= point(1,j) .and. point(2,j) <= y1 ) then

      boundary(1,j) = x3
      boundary(2,j) = y1

    else if ( x3 <= point(1,j) .and. &
      y1 <= point(2,j) .and. point(2,j) <= y2 ) then

      boundary(1,j) = x3
      boundary(2,j) = point(2,j)

    else if ( x3 <= point(1,j) .and. &
      y2 <= point(2,j) .and. point(2,j) <= y3 ) then

      boundary(1,j) = x3
      boundary(2,j) = y2

    else if ( x3 <= point(1,j) .and. y3 <= point(2,j) .and. &
      point(2,j) <= point(1,j) ) then

      boundary(1,j) = x3
      boundary(2,j) = y2

    else if ( x3 <= point(1,j) .and. y3 <= point(2,j) .and. &
      point(1,j) <= point(2,j) ) then

      boundary(1,j) = x2
      boundary(2,j) = y3

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P11_BOUNDARY_NEAREST - Fatal error!'
      write ( *, '(a,2g14.6)' ) '  Orphan point = ', point(1:2,j)
      stop

    end if

  end do

  return
end
subroutine p11_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P11_BOUNDARY_PROJECT projects exterior points to the boundary in problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, exterior points have been projected.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ), dimension(2,4) :: q1 = reshape ( (/ &
    0.00D+00, 0.00D+00, &
    1.00D+00, 0.00D+00, &
    0.75D+00, 0.25D+00, &
    0.25D+00, 0.25D+00 /), (/ 2, 4 /) )
  real ( kind = 8 ), dimension(2,4) :: q2 = reshape ( (/ &
    0.00D+00, 0.00D+00, &
    0.25D+00, 0.25D+00, &
    0.25D+00, 0.75D+00, &
    0.00D+00, 1.00D+00 /), (/ 2, 4 /) )
  real ( kind = 8 ), dimension(2,4) :: q3 = reshape ( (/ &
    0.50D+00, 0.50D+00, &
    0.50D+00, 0.25D+00, &
    0.75D+00, 0.25D+00, &
    1.00D+00, 0.50D+00 /), (/ 2, 4 /) )
  real ( kind = 8 ), dimension(2,4) :: q4 = reshape ( (/ &
    0.50D+00, 0.50D+00, &
    0.25D+00, 0.50D+00, &
    0.25D+00, 0.25D+00, &
    0.50D+00, 0.25D+00 /), (/ 2, 4 /) )
  real ( kind = 8 ), dimension(2,4) :: q5 = reshape ( (/ &
    0.50D+00, 0.50D+00, &
    0.50D+00, 1.00D+00, &
    0.25D+00, 0.75D+00, &
    0.25D+00, 0.50D+00 /), (/ 2, 4 /) )
  real ( kind = 8 ), dimension(2,3) :: t1 = reshape ( (/ &
    1.00D+00, 0.00D+00, &
    1.00D+00, 0.50D+00, &
    0.75D+00, 0.25D+00 /), (/ 2, 3 /) )
  real ( kind = 8 ), dimension(2,3) :: t2 = reshape ( (/ &
    0.00D+00, 1.00D+00, &
    0.25D+00, 0.25D+00, &
    0.50D+00, 1.00D+00 /), (/ 2, 3 /) )
  real ( kind = 8 ), dimension(2,3) :: t3 = reshape ( (/ &
    0.50D+00, 0.50D+00, &
    1.00D+00, 0.50D+00, &
    1.00D+00, 1.00D+00 /), (/ 2, 3 /) )
  real ( kind = 8 ), dimension(2,3) :: t4 = reshape ( (/ &
    0.50D+00, 0.50D+00, &
    1.00D+00, 1.00D+00, &
    0.50D+00, 1.00D+00 /), (/ 2, 3 /) )
  logical triangle_contains_point_2d
  real ( kind = 8 ), parameter :: x1 =  0.0D+00
  real ( kind = 8 ), parameter :: x2 =  0.5D+00
  real ( kind = 8 ), parameter :: x3 = +1.0D+00
  real ( kind = 8 ), parameter :: y1 =  0.0D+00
  real ( kind = 8 ), parameter :: y2 =  0.5D+00
  real ( kind = 8 ), parameter :: y3 = +1.0D+00

  do j = 1, n

    if ( point(1,j) <= x1 .and. point(2,j) <= y1 ) then

      point(1,j) = x1
      point(2,j) = y1

    else if ( point(1,j) <= x1 .and. point(2,j) <= y3 ) then

      point(1,j) = x1

    else if ( point(1,j) <= x1 .and. y3 <= point(2,j) ) then

      point(1,j) = x1
      point(2,j) = y3

    else if ( point(1,j) <= x3 .and. point(2,j) <= y1 ) then

      point(2,j) = y1

    else if ( x1 <= point(1,j) .and. point(1,j) <= x2 .and. &
      y3 <= point(2,j) ) then

      point(2,j) = y3

    else if ( triangle_contains_point_2d ( t3, point(1:2,j) ) ) then

      point(2,j) = y2

    else if ( triangle_contains_point_2d ( t4, point(1:2,j) ) ) then

      point(1,j) = x2

    else if ( x2 <= point(1,j) .and. point(1,j) <= x3 .and. &
      y3 <= point(2,j) ) then

      point(1,j) = x2
      point(2,j) = y3

    else if ( x3 <= point(1,j) .and. point(2,j) <= y1 ) then

      point(1,j) = x3
      point(2,j) = y1

    else if ( x3 <= point(1,j) .and. &
      y1 <= point(2,j) .and. point(2,j) <= y2 ) then

      point(1,j) = x3

    else if ( x3 <= point(1,j) .and. &
      y2 <= point(2,j) .and. point(2,j) <= y3 ) then

      point(1,j) = x3
      point(2,j) = y2

    else if ( x3 <= point(1,j) .and. y3 <= point(2,j) .and. &
      point(2,j) <= point(1,j) ) then

      point(1,j) = x3
      point(2,j) = y2

    else if ( x3 <= point(1,j) .and. y3 <= point(2,j) .and. &
      point(1,j) <= point(2,j) ) then

      point(1,j) = x2
      point(2,j) = y3

    end if

  end do

  return
end
subroutine p11_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P11_BOUNDARY_SEGMENT returns a boundary segment in problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  real ( kind = 8 ) s(m,6)
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) segment(m,segment_length)

  if ( segment_index == 1 ) then

    n1 = nint ( real ( segment_length - 1, kind = 8 ) &
              / real ( 4, kind = 8 ) )

    n2 = nint ( ( 1.0D+00 + 0.5D+00 ) &
              * real ( segment_length - 1, kind = 8 ) &
              / real ( 4, kind = 8 ) ) - n1

    n3 = nint ( ( 1.0D+00 + 0.5D+00 + ( 1.0D+00 - 0.5D+00 ) ) &
                * real ( segment_length - 1, kind = 8 ) &
                / real ( 4, kind = 8 ) ) - n1 - n2

    n4 = nint ( ( 1.0D+00 + 0.5D+00 + ( 1.0D+00 - 0.5D+00 ) &
                + ( 1.0D+00 - 0.5D+00 ) ) &
                * real ( segment_length - 1, kind = 8 ) &
                / real ( 4, kind = 8 ) ) - n1 - n2 - n3

    n5 = nint ( ( 1.0D+00 + 0.5D+00 + ( 1.0D+00 - 0.5D+00 ) &
                + ( 1.0D+00 - 0.5D+00 ) + 0.5D+00 ) &
                * real ( segment_length - 1, kind = 8 ) &
                / real ( 4, kind = 8 ) ) - n1 - n2 - n3 - n4

    n6 = segment_length - 1 - n1 - n2 - n3 - n4 - n5

    s(1:2,1) = (/  0.0D+00,  0.0D+00 /)
    s(1:2,2) = (/  1.0D+00,  0.0D+00 /)
    s(1:2,3) = (/  1.0D+00,  0.5D+00 /)
    s(1:2,4) = (/  0.5D+00,  0.5D+00 /)
    s(1:2,5) = (/  0.5D+00,  1.0D+00 /)
    s(1:2,6) = (/  0.0D+00,  1.0D+00 /)

    j = 0

    do i = 1, n1
      j = j + 1
      segment(1:2,j) = ( real ( n1 - i + 1, kind = 8 ) * s(1:2,1)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,2) ) &
                       / real ( n1,         kind = 8 )
    end do

    do i = 1, n2
      j = j + 1
      segment(1:2,j) = ( real ( n2 - i + 1, kind = 8 ) * s(1:2,2)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,3) ) &
                       / real ( n2,         kind = 8 )
    end do

    do i = 1, n3
      j = j + 1
      segment(1:2,j) = ( real ( n3 - i + 1, kind = 8 ) * s(1:2,3)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,4) ) &
                       / real ( n3,         kind = 8 )
    end do

    do i = 1, n4
      j = j + 1
      segment(1:2,j) = ( real ( n4 - i + 1, kind = 8 ) * s(1:2,4)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,5) ) &
                       / real ( n4,         kind = 8 )
    end do

    do i = 1, n5
      j = j + 1
      segment(1:2,j) = ( real ( n5 - i + 1, kind = 8 ) * s(1:2,5)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,6) ) &
                       / real ( n5,         kind = 8 )
    end do

    do i = 1, n6
      j = j + 1
      segment(1:2,j) = ( real ( n6 - i + 1, kind = 8 ) * s(1:2,6)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,1) ) &
                       / real ( n6,         kind = 8 )
    end do

    j = j + 1
    segment(1:2,j) = s(1:2,1)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P11_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p11_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P11_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P11_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    n = nint ( 4.0D+00 / h )
    n = max ( n, 5 )
    segment_length = n + mod ( 4 - mod ( n - 1, 4 ), 4 )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P11_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p11_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P11_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 1

  return
end
subroutine p11_box ( m, lo, hi )

!*****************************************************************************80
!
!! P11_BOX returns a bounding box for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)

  lo(1:m) = (/  0.0D+00,  0.0D+00 /)
  hi(1:m) = (/ +1.0D+00, +1.0D+00 /)

  return
end
subroutine p11_density ( m, n, point, density )

!*****************************************************************************80
!
!! P11_DENSITY returns the density for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p11_element_size ( element_size )

!*****************************************************************************80
!
!! P11_ELEMENT_SIZE returns a typical element size for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.15D+00

  return
end
subroutine p11_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P11_FIXED_NUM returns the number of fixed points in problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 6

  return
end
subroutine p11_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P11_FIXED_POINTS returns the fixed points in problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) fixed(m,fixed_num)

  fixed(1:m,1:fixed_num) = reshape ( (/ &
    0.0D+00,  0.0D+00, &
    1.0D+00,  0.0D+00, &
    1.0D+00,  0.5D+00,    &
    0.5D+00,  0.5D+00,    &
    0.5D+00,  1.0D+00, &
    0.0D+00,  1.0D+00 /), (/ m, fixed_num /) )

  return
end
subroutine p11_header ( )

!*****************************************************************************80
!
!! P11_HEADER prints some information about problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num

  call p11_boundary_segment_num ( boundary_segment_num )
  call p11_fixed_num ( fixed_num )
  call p11_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P11:'
  write ( *, '(a)' ) '  The L-shaped region.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The lower left corner of the indentation is'
  write ( *, '(a,2g14.6)' ) '  P = ', 0.5D+00, 0.5D+00
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p11_hole_num ( hole_num )

!*****************************************************************************80
!
!! P11_HOLE_NUM counts the holes in problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 0

  return
end
subroutine p11_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P11_HOLE_POINT returns a point inside a given hole in problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  return
end
subroutine p11_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P11_INSIDE reports if a point is inside the region in problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical inside(n)
  real ( kind = 8 ) point(m,n)

  inside(1:n) =                                                      &
       ( 0.0D+00 <= point(1,1:n) .and. point(1,1:n) <= 0.5D+00 .and. &
         0.0D+00 <= point(2,1:n) .and. point(2,1:n) <= 1.0D+00    )  &
    .or.                                                             &
       ( 0.5D+00 <= point(1,1:n) .and. point(1,1:n) <= 1.0D+00 .and. &
         0.0D+00 <= point(2,1:n) .and. point(2,1:n) <= 0.5D+00 )

  return
end
subroutine p11_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P11_SAMPLE samples points from the region in problem 11.
!
!  Discussion:
!
!    With a little bit of work, we can guarantee that we don't have to
!    use a rejection method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) area
  real ( kind = 8 ) prob(n)
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  x1 = 0.0D+00
  x2 = 0.5D+00
  x3 = 1.0D+00

  y1 = 0.0D+00
  y2 = 0.5D+00
  y3 = 1.0D+00

  area = ( x3 - x1 ) * ( y3 - y1 ) - ( x3 - x2 ) * ( y3 - y2 )
!
!  Generate a batch of points in [0,1]x[0,1].
!
  call r8mat_uniform_01 ( m, n, seed, point )
!
!  Generate a batch of N probabilities.
!
  call r8vec_uniform_01 ( n, seed, prob )
!
!  Map some points into [X1,X2] x [Y1,Y3].
!
  where ( prob(1:n) < ( x2 - x1 ) * ( y3 - y1 ) / area )

    point(1,1:n) = x1 + point(1,1:n) * ( x2 - x1 )
    point(2,1:n) = y1 + point(2,1:n) * ( y3 - y1 )
!
!  Map the other points into [X2,X3] x [Y1,Y2].
!
  elsewhere

    point(1,1:n) = x2 + point(1,1:n) * ( x3 - x2 )
    point(2,1:n) = y1 + point(2,1:n) * ( y2 - y1 )

  end where

  return
end
subroutine p11_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P11_SAMPLE_H1 samples points from the enlarged region in problem 11.
!
!  Discussion:
!
!    With a little bit of work, we can guarantee that we don't have to
!    use a rejection method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) area
  real ( kind = 8 ) h
  real ( kind = 8 ) prob(n)
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  x1 = 0.0D+00 - h
  x2 = 0.5D+00 + h
  x3 = 1.0D+00 + h

  y1 = 0.0D+00 - h
  y2 = 0.5D+00 + h
  y3 = 1.0D+00 + h

  area = ( x3 - x1 ) * ( y3 - y1 ) - ( x3 - x2 ) * ( y3 - y2 )
!
!  Generate a batch of points in [0,1]x[0,1].
!
  call r8mat_uniform_01 ( m, n, seed, point )
!
!  Generate a batch of N probabilities.
!
  call r8vec_uniform_01 ( n, seed, prob )
!
!  Map some points into [X1,X2] x [Y1,Y3].
!
  where ( prob(1:n) < ( x2 - x1 ) * ( y3 - y1 ) / area )

    point(1,1:n) = x1 + point(1,1:n) * ( x2 - x1 )
    point(2,1:n) = y1 + point(2,1:n) * ( y3 - y1 )
!
!  Map the other points into [X2,X3] x [Y1,Y2].
!
  elsewhere

    point(1,1:n) = x2 + point(1,1:n) * ( x3 - x2 )
    point(2,1:n) = y1 + point(2,1:n) * ( y2 - y1 )

  end where

  return
end
subroutine p11_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P11_SDIST returns the signed distance to the region in problem 11.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P11_SDIST - Fatal error!'
  write ( *, '(a)' ) '  This routine is not written yet!'

  stop
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
!    14 January 2005
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

  title = '#11: The L-shaped region.'

  return
end
subroutine p12_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P12_BOUNDARY_NEAREST returns a nearest boundary point in problem 12.
!
!  Discussion:
!
!    Actually, we only need "reasonably accurate" values when the point
!    to be tested is outside the region, and not to far from it.  Since
!    the politics of interior propinquity are surprisingly complicated for
!    this region, we will plump for a simple, approximate scheme.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( m, n ) :: boundary
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ), dimension ( m, 4 ) :: q1
  real ( kind = 8 ), dimension ( m, 4 ) :: q2
  real ( kind = 8 ), dimension ( m, 4 ) :: q3
  real ( kind = 8 ), dimension ( m, 4 ) :: q4
  real ( kind = 8 ), dimension ( m, 4 ) :: q5
  real ( kind = 8 ), dimension ( m, 4 ) :: q6
  real ( kind = 8 ), dimension ( m, 4 ) :: q7
  real ( kind = 8 ), dimension ( m, 4 ) :: q8
  real ( kind = 8 ), dimension ( m, 4 ) :: q9
  real ( kind = 8 ), dimension ( m, 4 ) :: q10
  real ( kind = 8 ), dimension ( m, 4 ) :: q11
  real ( kind = 8 ), dimension ( m, 4 ) :: q12
  real ( kind = 8 ), dimension ( m, 4 ) :: q13
  logical quad_contains_point_2d
  real ( kind = 8 ), dimension ( m, 3 ) :: t1
  real ( kind = 8 ), dimension ( m, 3 ) :: t2
  real ( kind = 8 ), dimension ( m, 3 ) :: t3
  real ( kind = 8 ), dimension ( m, 3 ) :: t4
  real ( kind = 8 ), dimension ( m, 3 ) :: t5
  real ( kind = 8 ), dimension ( m, 3 ) :: t6
  real ( kind = 8 ), dimension ( m, 3 ) :: t7
  logical triangle_contains_point_2d
  real ( kind = 8 ), parameter :: x1 = 0.000D+00
  real ( kind = 8 ), parameter :: x2 = 0.500D+00
  real ( kind = 8 ), parameter :: x2p5 = 0.5625D+00
  real ( kind = 8 ), parameter :: x3 = 0.625D+00
  real ( kind = 8 ), parameter :: x3p5 = 0.8125D+00
  real ( kind = 8 ), parameter :: x4 = 1.000D+00
  real ( kind = 8 ), parameter :: y1 = 0.000D+00
  real ( kind = 8 ), parameter :: y2 = 0.25D+00
  real ( kind = 8 ), parameter :: y2p5 = 0.31250D+00
  real ( kind = 8 ), parameter :: y3 = 0.375D+00
  real ( kind = 8 ), parameter :: y4 = 1.000D+00

  t1(1:2,1:3) = reshape ( (/ &
    0.0000D+00, 0.0000D+00, &
    0.5000D+00, 0.0000D+00, &
    0.2500D+00, 0.2500D+00 /), (/ 2, 3 /) )

  q1(1:2,1:4) = reshape ( (/ &
    0.0000D+00, 0.0000D+00, &
    0.2500D+00, 0.2500D+00, &
    0.2500D+00, 0.7500D+00, &
    0.0000D+00, 1.0000D+00 /), (/ 2, 4 /) )

  t2(1:2,1:3) = reshape ( (/ &
    0.5000D+00, 0.0000D+00, &
    0.5000D+00, 0.2500D+00, &
    0.2500D+00, 0.2500D+00 /), (/ 2, 3 /) )

  q2(1:2,1:4) = reshape ( (/ &
    0.2500D+00, 0.2500D+00, &
    0.5000D+00, 0.2500D+00, &
    0.5000D+00, 0.3125D+00, &
    0.2500D+00, 0.3125D+00 /), (/ 2, 4 /) )

  q3(1:2,1:4) = reshape ( (/ &
    0.2500D+00, 0.3125D+00, &
    0.5000D+00, 0.2125D+00, &
    0.5000D+00, 0.3750D+00, &
    0.2500D+00, 0.375D+00 /), (/ 2, 4 /) )

  q4(1:2,1:4) = reshape ( (/ &
    0.5000D+00, 0.3750D+00, &
    0.5000D+00, 1.0000D+00, &
    0.2500D+00, 0.7500D+00, &
    0.2500D+00, 0.375D+00 /), (/ 2, 4 /) )

  t3(1:2,1:3) = reshape ( (/ &
    0.5000D+00, 1.0000D+00, &
    0.0000D+00, 1.0000D+00, &
    0.2500D+00, 0.7500D+00 /), (/ 2, 3 /) )

  q5(1:2,1:4) = reshape ( (/ &
    0.5000D+00, 0.0000D+00, &
    0.5625D+00, 0.0000D+00, &
    0.5625D+00, 0.0625D+00, &
    0.5000D+00, 0.1250D+00 /), (/ 2, 4 /) )

  q6(1:2,1:4) = reshape ( (/ &
    0.6250D+00, 0.0000D+00, &
    0.6250D+00, 0.1250D+00, &
    0.6625D+00, 0.0625D+00, &
    0.6625D+00, 0.0000D+00 /), (/ 2, 4 /) )

  t4(1:2,1:3) = reshape ( (/ &
    0.5000D+00, 0.1250D+00, &
    0.5625D+00, 0.0625D+00, &
    0.6250D+00, 0.1250D+00 /), (/ 2, 3 /) )

  t5(1:2,1:3) = reshape ( (/ &
    0.5000D+00, 0.3750D+00, &
    0.6250D+00, 0.3750D+00, &
    0.5625D+00, 0.4375D+00 /), (/ 2, 3 /) )

  q7(1:2,1:4) = reshape ( (/ &
    0.5000D+00, 0.3750D+00, &
    0.5625D+00, 0.4375D+00, &
    0.5625D+00, 1.0000D+00, &
    0.5000D+00, 1.0000D+00 /), (/ 2, 4 /) )

  q8(1:2,1:4) = reshape ( (/ &
    0.6250D+00, 0.3750D+00, &
    0.6250D+00, 1.0000D+00, &
    0.5626D+00, 1.0000D+00, &
    0.5625D+00, 0.4375D+00 /), (/ 2, 4 /) )

  t6(1:2,1:3) = reshape ( (/ &
    0.6250D+00, 0.0000D+00, &
    1.0000D+00, 0.0000D+00, &
    0.8125D+00, 0.1875D+00 /), (/ 2, 3 /) )

  q9(1:2,1:4) = reshape ( (/ &
    0.6250D+00, 0.0000D+00, &
    0.8125D+00, 0.1875D+00, &
    0.8125D+00, 0.2500D+00, &
    0.6250D+00, 0.2500D+00 /), (/ 2, 4 /) )

  q10(1:2,1:4) = reshape ( (/ &
    0.6250D+00, 0.2500D+00, &
    0.8125D+00, 0.2500D+00, &
    0.8125D+00, 0.3125D+00, &
    0.6250D+00, 0.3125D+00 /), (/ 2, 4 /) )

  q11(1:2,1:4) = reshape ( (/ &
    0.6250D+00, 0.3125D+00, &
    0.8125D+00, 0.3125D+00, &
    0.8125D+00, 0.3750D+00, &
    0.6250D+00, 0.3750D+00 /), (/ 2, 4 /) )

  q12(1:2,1:4) = reshape ( (/ &
    0.6250D+00, 0.3750D+00, &
    0.8125D+00, 0.3750D+00, &
    0.8125D+00, 0.8125D+00, &
    0.6250D+00, 1.0000D+00 /), (/ 2, 4 /) )

  q13(1:2,1:4) = reshape ( (/ &
    1.0000D+00, 0.0000D+00, &
    1.0000D+00, 1.0000D+00, &
    0.8125D+00, 0.8125D+00, &
    0.8125D+00, 0.1875D+00 /), (/ 2, 4 /) )

  t7(1:2,1:3) = reshape ( (/ &
    0.6250D+00, 1.0000D+00, &
    0.8125D+00, 0.8125D+00, &
    1.0000D+00, 1.0000D+00 /), (/ 2, 3 /) )

  do j = 1, n
!
!  Column 1
!
    if ( point(1,j) <= x1 .and. point(2,j) <= y1 ) then

      boundary(1,j) = x1
      boundary(2,j) = y1

    else if ( point(1,j) <= x1 .and. point(2,j) <= y4 ) then

      boundary(1,j) = x1
      boundary(2,j) = point(2,j)

    else if ( point(1,j) <= x1 .and. y4 <= point(2,j) ) then

      boundary(1,j) = x1
      boundary(2,j) = y4
!
!  Column 2
!
    else if ( point(1,j) <= x2 .and. point(2,j) <= y1 ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y1

    else if ( triangle_contains_point_2d ( t1, point(1:2,j) ) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y1

    else if ( quad_contains_point_2d ( q1, point(1:2,j) ) ) then

      boundary(1,j) = x1
      boundary(2,j) = point(2,j)

    else if ( quad_contains_point_2d ( t2, point(1:2,j) ) ) then

      boundary(1,j) = x2
      boundary(2,j) = point(2,j)

    else if ( quad_contains_point_2d ( q2, point(1:2,j) ) ) then

      boundary(1,j) = x2
      boundary(2,j) = y2

    else if ( quad_contains_point_2d ( q3, point(1:2,j) ) ) then

      boundary(1,j) = x2
      boundary(2,j) = y3

    else if ( quad_contains_point_2d ( q4, point(1:2,j) ) ) then

      boundary(1,j) = x2
      boundary(2,j) = point(2,j)

    else if ( triangle_contains_point_2d ( t3, point(1:2,j) ) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y4

    else if ( point(1,j) <= x2 .and. y4 <= point(2,j) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y4
!
!  Under the H crossbar.
!
    else if ( point(1,j) <= x2p5 .and. point(2,j) <= y1 ) then

      boundary(1,j) = x2
      boundary(2,j) = y1

    else if ( point(1,j) <= x3 .and. point(2,j) <= y1 ) then

      boundary(1,j) = x3
      boundary(2,j) = y1

    else if ( quad_contains_point_2d ( q5, point(1:2,j) ) ) then

      boundary(1,j) = x2
      boundary(2,j) = point(2,j)

    else if ( quad_contains_point_2d ( q6, point(1:2,j) ) ) then

      boundary(1,j) = x3
      boundary(2,j) = point(2,j)

    else if ( triangle_contains_point_2d ( t4, point(1:2,j) ) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y2
!
!  In the H crossbar.
!
    else if ( point(1,j) <= x3 .and. &
              y2 <= point(2,j) .and. point(2,j) <= y2p5 ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y2

    else if ( point(1,j) <= x3 .and. &
              y2p5 <= point(2,j) .and. point(2,j) <= y3 ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y3
!
!  Above the H crossbar.
!
    else if ( triangle_contains_point_2d ( t5, point(1:2,j) ) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y3

    else if ( quad_contains_point_2d ( q7, point(1:2,j) ) ) then

      boundary(1,j) = x2
      boundary(2,j) = point(2,j)

    else if ( quad_contains_point_2d ( q8, point(1:2,j) ) ) then

      boundary(1,j) = x3
      boundary(2,j) = point(2,j)

    else if ( point(1,j) <= x2p5 .and. y4 <= point(2,j) ) then

      boundary(1,j) = x2
      boundary(2,j) = y4

    else if ( point(1,j) <= x3 .and. y4 <= point(2,j) ) then

      boundary(1,j) = x3
      boundary(2,j) = y4
!
!  Column 4
!
    else if ( point(1,j) <= x4 .and. point(2,j) <= y1 ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y1

    else if ( triangle_contains_point_2d ( t6, point(1:2,j) ) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y1

    else if ( quad_contains_point_2d ( q9, point(1:2,j) ) ) then

      boundary(1,j) = x3
      boundary(2,j) = point(2,j)

    else if ( quad_contains_point_2d ( q10, point(1:2,j) ) ) then

      boundary(1,j) = x3
      boundary(2,j) = y2

    else if ( quad_contains_point_2d ( q11, point(1:2,j) ) ) then

      boundary(1,j) = x3
      boundary(2,j) = y3

    else if ( quad_contains_point_2d ( q12, point(1:2,j) ) ) then

      boundary(1,j) = x3
      boundary(2,j) = point(2,j)

    else if ( quad_contains_point_2d ( q13, point(1:2,j) ) ) then

      boundary(1,j) = x4
      boundary(2,j) = point(2,j)

    else if ( triangle_contains_point_2d ( t7, point(1:2,j) ) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y4

    else if ( point(1,j) <= x4 .and. y4 <= point(2,j) ) then

      boundary(1,j) = point(1,j)
      boundary(2,j) = y4
!
!  Column 5
!
    else if ( x4 <= point(1,j) .and. point(2,j) <= y1 ) then

      boundary(1,j) = x4
      boundary(2,j) = y1

    else if ( x4 <= point(1,j) .and. point(2,j) <= y4 ) then

      boundary(1,j) = x4
      boundary(2,j) = point(2,j)

    else if ( x4 <= point(1,j) .and. y4 <= point(2,j) ) then

      boundary(1,j) = x4
      boundary(2,j) = y4
    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P12_BOUNDARY_NEAREST - Fatal error!'
      write ( *, '(a,2g14.6)' ) '  Orphan point = ', point(1:2,j)
      stop

    end if

  end do

  return
end
subroutine p12_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P12_BOUNDARY_PROJECT projects exterior points to the boundary in problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, exterior points have been projected.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ), dimension ( m, 4 ) :: q5
  real ( kind = 8 ), dimension ( m, 4 ) :: q6
  real ( kind = 8 ), dimension ( m, 4 ) :: q7
  real ( kind = 8 ), dimension ( m, 4 ) :: q8
  logical quad_contains_point_2d
  real ( kind = 8 ), dimension ( m, 3 ) :: t4
  real ( kind = 8 ), dimension ( m, 3 ) :: t5
  logical triangle_contains_point_2d
  real ( kind = 8 ), parameter :: x1 = 0.000D+00
  real ( kind = 8 ), parameter :: x2 = 0.500D+00
  real ( kind = 8 ), parameter :: x2p5 = 0.5625D+00
  real ( kind = 8 ), parameter :: x3 = 0.625D+00
  real ( kind = 8 ), parameter :: x3p5 = 0.8125D+00
  real ( kind = 8 ), parameter :: x4 = 1.000D+00
  real ( kind = 8 ), parameter :: y1 = 0.000D+00
  real ( kind = 8 ), parameter :: y2 = 0.125D+00
  real ( kind = 8 ), parameter :: y2p5 = 0.250D+00
  real ( kind = 8 ), parameter :: y3 = 0.375D+00
  real ( kind = 8 ), parameter :: y4 = 1.000D+00

  q5(1:2,1:4) = reshape ( (/ &
    0.5000D+00, 0.0000D+00, &
    0.5625D+00, 0.0000D+00, &
    0.5625D+00, 0.0625D+00, &
    0.5000D+00, 0.1250D+00 /), (/ 2, 4 /) )

  q6(1:2,1:4) = reshape ( (/ &
    0.6250D+00, 0.0000D+00, &
    0.6250D+00, 0.1250D+00, &
    0.6625D+00, 0.0625D+00, &
    0.6625D+00, 0.0000D+00 /), (/ 2, 4 /) )

  t4(1:2,1:3) = reshape ( (/ &
    0.5000D+00, 0.1250D+00, &
    0.5625D+00, 0.0625D+00, &
    0.6250D+00, 0.1250D+00 /), (/ 2, 3 /) )

  t5(1:2,1:3) = reshape ( (/ &
    0.5000D+00, 0.3750D+00, &
    0.6250D+00, 0.3750D+00, &
    0.5625D+00, 0.4375D+00 /), (/ 2, 3 /) )

  q7(1:2,1:4) = reshape ( (/ &
    0.5000D+00, 0.3750D+00, &
    0.5625D+00, 0.4375D+00, &
    0.5625D+00, 1.0000D+00, &
    0.5000D+00, 1.0000D+00 /), (/ 2, 4 /) )

  q8(1:2,1:4) = reshape ( (/ &
    0.6250D+00, 0.3750D+00, &
    0.6250D+00, 1.0000D+00, &
    0.5626D+00, 1.0000D+00, &
    0.5625D+00, 0.4375D+00 /), (/ 2, 4 /) )

  do j = 1, n
!
!  Column 1
!
    if ( point(1,j) <= x1 .and. point(2,j) <= y1 ) then

      point(1,j) = x1
      point(2,j) = y1

    else if ( point(1,j) <= x1 .and. point(2,j) <= y4 ) then

      point(1,j) = x1

    else if ( point(1,j) <= x1 .and. y4 <= point(2,j) ) then

      point(1,j) = x1
      point(2,j) = y4
!
!  Column 2
!
    else if ( point(1,j) <= x2 .and. point(2,j) <= y1 ) then

      point(2,j) = y1

    else if ( point(1,j) <= x2 .and. y4 <= point(2,j) ) then

      point(2,j) = y4
!
!  Under the H crossbar.
!
    else if ( point(1,j) <= x2p5 .and. point(2,j) <= y1 ) then

      point(1,j) = x2
      point(2,j) = y1

    else if ( point(1,j) <= x3 .and. point(2,j) <= y1 ) then

      point(1,j) = x3
      point(2,j) = y1

    else if ( quad_contains_point_2d ( q5, point(1:2,j) ) ) then

      point(1,j) = x2

    else if ( quad_contains_point_2d ( q6, point(1:2,j) ) ) then

      point(1,j) = x3

    else if ( triangle_contains_point_2d ( t4, point(1:2,j) ) ) then

      point(2,j) = y2
!
!  Above the H crossbar.
!
    else if ( triangle_contains_point_2d ( t5, point(1:2,j) ) ) then

      point(2,j) = y3

    else if ( quad_contains_point_2d ( q7, point(1:2,j) ) ) then

      point(1,j) = x2

    else if ( quad_contains_point_2d ( q8, point(1:2,j) ) ) then

      point(1,j) = x3

    else if ( point(1,j) <= x2p5 .and. y4 <= point(2,j) ) then

      point(1,j) = x2
      point(2,j) = y4

    else if ( point(1,j) <= x3 .and. y4 <= point(2,j) ) then

      point(1,j) = x3
      point(2,j) = y4
!
!  Column 4
!
    else if ( point(1,j) <= x4 .and. point(2,j) <= y1 ) then

      point(2,j) = y1

    else if ( point(1,j) <= x4 .and. y4 <= point(2,j) ) then

      point(2,j) = y4
!
!  Column 5
!
    else if ( x4 <= point(1,j) .and. point(2,j) <= y1 ) then

      point(1,j) = x4
      point(2,j) = y1

    else if ( x4 <= point(1,j) .and. point(2,j) <= y4 ) then

      point(1,j) = x4

    else if ( x4 <= point(1,j) .and. y4 <= point(2,j) ) then

      point(1,j) = x4
      point(2,j) = y4

    end if

  end do

  return
end
subroutine p12_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P12_BOUNDARY_SEGMENT returns a boundary segment in problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n01
  integer ( kind = 4 ) n02
  integer ( kind = 4 ) n03
  integer ( kind = 4 ) n04
  integer ( kind = 4 ) n05
  integer ( kind = 4 ) n06
  integer ( kind = 4 ) n07
  integer ( kind = 4 ) n08
  integer ( kind = 4 ) n09
  integer ( kind = 4 ) n10
  integer ( kind = 4 ) n11
  integer ( kind = 4 ) n12
  real ( kind = 8 ) s(m,12)
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) segment(m,segment_length)

  if ( segment_index == 1 ) then

    n01 = nint ( 0.500D+00 * real ( segment_length - 1, kind = 8 ) &
              / 5.75D+00 )

    n02 = nint ( ( 0.500D+00 + 0.250D+00) &
              * real ( segment_length - 1, kind = 8 ) &
              / real ( 5.75D+00, kind = 8 ) ) - n01

    n03 = nint ( ( 0.500D+00 + 0.250D+00 + 0.125D+00 ) &
              * real ( segment_length - 1, kind = 8 ) &
              / real ( 5.75D+00, kind = 8 ) ) - n01 - n02

    n04 = nint ( ( 0.500D+00 + 0.250D+00 + 0.125D+00 + 0.250D+00 ) &
              * real ( segment_length - 1, kind = 8 ) &
              / real ( 5.75D+00, kind = 8 ) ) - n01 - n02 - n03

    n05 = nint ( ( 0.500D+00 + 0.250D+00 + 0.125D+00 + 0.250D+00 &
                 + 0.375D+00 ) &
              * real ( segment_length - 1, kind = 8 ) &
              / real ( 5.75D+00, kind = 8 ) ) - n01 - n02 - n03 - n04

    n06 = nint ( ( 0.500D+00 + 0.250D+00 + 0.125D+00 + 0.250D+00 &
                 + 0.375D+00 + 1.000D+00 ) &
              * real ( segment_length - 1, kind = 8 ) &
              / real ( 5.75D+00, kind = 8 ) ) - n01 - n02 - n03 - n04 - n05

    n07 = nint ( ( 0.500D+00 + 0.250D+00 + 0.125D+00 + 0.250D+00 &
                 + 0.375D+00 + 1.000D+00 + 0.375D+00 ) &
              * real ( segment_length - 1, kind = 8 ) &
              / real ( 5.75D+00, kind = 8 ) ) - n01 - n02 - n03 - n04 - n05 &
              - n06

    n08 = nint ( ( 0.500D+00 + 0.250D+00 + 0.125D+00 + 0.250D+00 &
                 + 0.375D+00 + 1.000D+00 + 0.375D+00 + 0.625D+00 ) &
              * real ( segment_length - 1, kind = 8 ) &
              / real ( 5.75D+00, kind = 8 ) ) - n01 - n02 - n03 - n04 - n05 &
              - n06 - n07

    n09 = nint ( ( 0.500D+00 + 0.250D+00 + 0.125D+00 + 0.250D+00 &
                 + 0.375D+00 + 1.000D+00 + 0.375D+00 + 0.625D+00 &
                 + 0.125D+00 ) &
              * real ( segment_length - 1, kind = 8 ) &
              / real ( 5.75D+00, kind = 8 ) ) - n01 - n02 - n03 - n04 - n05 &
              - n06 - n07 - n08

    n10 = nint ( ( 0.500D+00 + 0.250D+00 + 0.125D+00 + 0.250D+00 &
                 + 0.375D+00 + 1.000D+00 + 0.375D+00 + 0.625D+00 &
                 + 0.125D+00 + 0.375D+00 ) &
              * real ( segment_length - 1, kind = 8 ) &
              / real ( 5.75D+00, kind = 8 ) ) - n01 - n02 - n03 - n04 - n05 &
              - n06 - n07 - n08 - n09

    n10 = nint ( ( 0.500D+00 + 0.250D+00 + 0.125D+00 + 0.250D+00 &
                 + 0.375D+00 + 1.000D+00 + 0.375D+00 + 0.625D+00 &
                 + 0.125D+00 + 0.625D+00 ) &
              * real ( segment_length - 1, kind = 8 ) &
              / real ( 5.75D+00, kind = 8 ) ) - n01 - n02 - n03 - n04 - n05 &
              - n06 - n07 - n08 - n09

    n11 = nint ( ( 0.500D+00 + 0.250D+00 + 0.125D+00 + 0.250D+00 &
                 + 0.375D+00 + 1.000D+00 + 0.375D+00 + 0.625D+00 &
                 + 0.125D+00 + 0.625D+00 + 0.500D+00 ) &
              * real ( segment_length - 1, kind = 8 ) &
              / real ( 5.75D+00, kind = 8 ) ) - n01 - n02 - n03 - n04 - n05 &
              - n06 - n07 - n08 - n09 - n10

    n12 = segment_length - 1 - n01 - n02 - n03 - n04 - n05 &
          - n06 - n07 - n08 - n09 - n10 - n11

    s(1:2,1:12) = reshape ( (/ &
    0.000D+00,  0.000D+00, &
    0.500D+00,  0.000D+00, &
    0.500D+00,  0.250D+00, &
    0.625D+00,  0.250D+00, &
    0.625D+00,  0.000D+00, &
    1.000D+00,  0.000D+00, &
    1.000D+00,  1.000D+00, &
    0.625D+00,  1.000D+00, &
    0.625D+00,  0.375D+00, &
    0.500D+00,  0.375D+00, &
    0.500D+00,  1.000D+00, &
    0.000D+00,  1.000D+00 /), (/ m, 12 /) )

    j = 0

    do i = 1, n01
      j = j + 1
      segment(1:2,j) = ( real ( n01 - i + 1, kind = 8 ) * s(1:2,1)   &
                       + real (       i - 1, kind = 8 ) * s(1:2,2) ) &
                       / real ( n01,         kind = 8 )
    end do

    do i = 1, n02
      j = j + 1
      segment(1:2,j) = ( real ( n02 - i + 1, kind = 8 ) * s(1:2,2)   &
                       + real (       i - 1, kind = 8 ) * s(1:2,3) ) &
                       / real ( n02,         kind = 8 )
    end do

    do i = 1, n03
      j = j + 1
      segment(1:2,j) = ( real ( n03 - i + 1, kind = 8 ) * s(1:2,3)   &
                       + real (       i - 1, kind = 8 ) * s(1:2,4) ) &
                       / real ( n03,         kind = 8 )
    end do

    do i = 1, n04
      j = j + 1
      segment(1:2,j) = ( real ( n04 - i + 1, kind = 8 ) * s(1:2,4)   &
                       + real (       i - 1, kind = 8 ) * s(1:2,5) ) &
                       / real ( n04,         kind = 8 )
    end do

    do i = 1, n05
      j = j + 1
      segment(1:2,j) = ( real ( n05 - i + 1, kind = 8 ) * s(1:2,5)   &
                       + real (       i - 1, kind = 8 ) * s(1:2,6) ) &
                       / real ( n05,         kind = 8 )
    end do

    do i = 1, n06
      j = j + 1
      segment(1:2,j) = ( real ( n06 - i + 1, kind = 8 ) * s(1:2,6)   &
                       + real (       i - 1, kind = 8 ) * s(1:2,7) ) &
                       / real ( n06,         kind = 8 )
    end do

    do i = 1, n07
      j = j + 1
      segment(1:2,j) = ( real ( n07 - i + 1, kind = 8 ) * s(1:2,7)   &
                       + real (       i - 1, kind = 8 ) * s(1:2,8) ) &
                       / real ( n07,         kind = 8 )
    end do

    do i = 1, n08
      j = j + 1
      segment(1:2,j) = ( real ( n08 - i + 1, kind = 8 ) * s(1:2,8)   &
                       + real (       i - 1, kind = 8 ) * s(1:2,9) ) &
                       / real ( n08,         kind = 8 )
    end do

    do i = 1, n09
      j = j + 1
      segment(1:2,j) = ( real ( n09 - i + 1, kind = 8 ) * s(1:2,9)   &
                       + real (       i - 1, kind = 8 ) * s(1:2,10) ) &
                       / real ( n09,         kind = 8 )
    end do

    do i = 1, n10
      j = j + 1
      segment(1:2,j) = ( real ( n10 - i + 1, kind = 8 ) * s(1:2,10)   &
                       + real (       i - 1, kind = 8 ) * s(1:2,11) ) &
                       / real ( n10,         kind = 8 )
    end do

    do i = 1, n11
      j = j + 1
      segment(1:2,j) = ( real ( n11 - i + 1, kind = 8 ) * s(1:2,11)   &
                       + real (       i - 1, kind = 8 ) * s(1:2,12) ) &
                       / real ( n11,         kind = 8 )
    end do

    do i = 1, n12
      j = j + 1
      segment(1:2,j) = ( real ( n12 - i + 1, kind = 8 ) * s(1:2,12)   &
                       + real (       i - 1, kind = 8 ) * s(1:2,1) ) &
                       / real ( n12,         kind = 8 )
    end do

    j = j + 1
    segment(1:2,j) = s(1:2,1)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P12_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p12_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P12_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P12_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    n = nint ( 5.75D+00 / h )
    n = max ( n, 5 )
    segment_length = n + mod ( 4 - mod ( n - 1, 4 ), 4 )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P12_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p12_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P12_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 1

  return
end
subroutine p12_box ( m, lo, hi )

!*****************************************************************************80
!
!! P12_BOX returns a bounding box for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)

  lo(1:m) = (/  0.0D+00,  0.0D+00 /)
  hi(1:m) = (/ +1.0D+00, +1.0D+00 /)

  return
end
subroutine p12_density ( m, n, point, density )

!*****************************************************************************80
!
!! P12_DENSITY returns the density for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p12_element_size ( element_size )

!*****************************************************************************80
!
!! P12_ELEMENT_SIZE returns a typical element size for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.15D+00

  return
end
subroutine p12_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P12_FIXED_NUM returns the number of fixed points in problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 12

  return
end
subroutine p12_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P12_FIXED_POINTS returns the fixed points in problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) fixed(m,fixed_num)

  fixed(1:m,1:fixed_num) = reshape ( (/ &
    0.000D+00,  0.000D+00, &
    0.500D+00,  0.000D+00, &
    0.500D+00,  0.250D+00, &
    0.625D+00,  0.250D+00, &
    0.625D+00,  0.000D+00, &
    1.000D+00,  0.000D+00, &
    1.000D+00,  1.000D+00, &
    0.625D+00,  1.000D+00, &
    0.625D+00,  0.375D+00, &
    0.500D+00,  0.375D+00, &
    0.500D+00,  1.000D+00, &
    0.000D+00,  1.000D+00 /), (/ m, fixed_num /) )

  return
end
subroutine p12_header ( )

!*****************************************************************************80
!
!! P12_HEADER prints some information about problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num

  call p12_boundary_segment_num ( boundary_segment_num )
  call p12_fixed_num ( fixed_num )
  call p12_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P12:'
  write ( *, '(a)' ) '  John Shadid''s H-shaped region.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p12_hole_num ( hole_num )

!*****************************************************************************80
!
!! P12_HOLE_NUM counts the holes in problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 0

  return
end
subroutine p12_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P12_HOLE_POINT returns a point inside a given hole in problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  return
end
subroutine p12_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P12_INSIDE reports if a point is inside the region in problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical inside(n)
  real ( kind = 8 ) point(m,n)

  inside(1:n) =                                                          &
       ( 0.000D+00 <= point(1,1:n) .and. point(1,1:n) <= 0.500D+00 .and. &
         0.000D+00 <= point(2,1:n) .and. point(2,1:n) <= 1.000D+00    )  &
    .or.                                                                 &
       ( 0.500D+00 <= point(1,1:n) .and. point(1,1:n) <= 0.625D+00 .and. &
         0.250D+00 <= point(2,1:n) .and. point(2,1:n) <= 0.375D+00 )     &
    .or.                                                                 &
       ( 0.625D+00 <= point(1,1:n) .and. point(1,1:n) <= 1.000D+00 .and. &
         0.000D+00 <= point(2,1:n) .and. point(2,1:n) <= 1.000D+00    )

  return
end
subroutine p12_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P12_SAMPLE samples points from the region in problem 12.
!
!  Discussion:
!
!    With a little bit of work, we can guarantee that we don't have to
!    use a rejection method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area2
  real ( kind = 8 ) area3
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) prob(n)
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: x1 = 0.000D+00
  real ( kind = 8 ) :: x2 = 0.500D+00
  real ( kind = 8 ) :: x3 = 0.625D+00
  real ( kind = 8 ) :: x4 = 1.000D+00
  real ( kind = 8 ) :: y1 = 0.000D+00
  real ( kind = 8 ) :: y2 = 0.250D+00
  real ( kind = 8 ) :: y3 = 0.375D+00
  real ( kind = 8 ) :: y4 = 1.000D+00

  area1 = ( x2 - x1 ) * ( y4 - y1 )
  area2 = ( x3 - x2 ) * ( y3 - y2 )
  area3 = ( x4 - x3 ) * ( y4 - y1 )

  area = area1 + area2 + area3

  p1 =   area1                   / area
  p2 = ( area1 + area2 )         / area
!
!  Generate a batch of points in [0,1]x[0,1].
!
  call r8mat_uniform_01 ( m, n, seed, point )
!
!  Generate a batch of N probabilities.
!
  call r8vec_uniform_01 ( n, seed, prob )
!
!  Map some points into [X1,X2] x [Y1,Y4].
!
  where ( prob(1:n) <= p1 )

    point(1,1:n) = x1 + point(1,1:n) * ( x2 - x1 )
    point(2,1:n) = y1 + point(2,1:n) * ( y4 - y1 )
!
!  Map points into [X2,X3] x [Y2,Y3].
!
  elsewhere ( prob(1:n) <= p2 )

    point(1,1:n) = x2 + point(1,1:n) * ( x3 - x2 )
    point(2,1:n) = y2 + point(2,1:n) * ( y3 - y2 )
!
!  Map points into [X3,X4] x [Y1,Y4].
!
  elsewhere

    point(1,1:n) = x3 + point(1,1:n) * ( x4 - x3 )
    point(2,1:n) = y1 + point(2,1:n) * ( y4 - y1 )

  end where

  return
end
subroutine p12_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P12_SAMPLE_H1 samples points from the enlarged region in problem 12.
!
!  Discussion:
!
!    With a little bit of work, we can guarantee that we don't have to
!    use a rejection method.
!
!    We enlarge the region by a layer H.  We do not round the
!    corners of the region, which should be done if we literally only
!    want to add points within H units of the boundary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area2
  real ( kind = 8 ) area3
  real ( kind = 8 ) h
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) prob(n)
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: x1 = 0.000D+00
  real ( kind = 8 ) :: x2 = 0.500D+00
  real ( kind = 8 ) :: x3 = 0.625D+00
  real ( kind = 8 ) :: x4 = 1.000D+00
  real ( kind = 8 ) :: y1 = 0.000D+00
  real ( kind = 8 ) :: y2 = 0.250D+00
  real ( kind = 8 ) :: y3 = 0.375D+00
  real ( kind = 8 ) :: y4 = 1.000D+00
!
!  Special, simple action if H is large enough that the
!  H has become a box.
!
  if ( 1.0D+00 <= 16.0D+00 * h ) then

    x1 = 0.000D+00 - h
    x2 = 1.000D+00 + h
    y1 = 0.000D+00 - h
    y2 = 1.000D+00 + h

    call r8mat_uniform_01 ( m, n, seed, point )

    point(1,1:n) = x1 + ( x2 - x1 ) * point(1,1:n)
    point(2,1:n) = y1 + ( y2 - y1 ) * point(2,1:n)

    return
  end if

  x1 = 0.000D+00 - h
  x2 = 0.500D+00 + h
  x3 = 0.625D+00 - h
  x4 = 1.000D+00 + h
  y1 = 0.000D+00 - h
  y2 = 0.250D+00 + h
  y3 = 0.375D+00 - h
  y4 = 1.000D+00 + h

  area1 = ( x2 - x1 ) * ( y4 - y1 )
  area2 = ( x3 - x2 ) * ( y3 - y2 )
  area3 = ( x4 - x3 ) * ( y4 - y1 )

  area = area1 + area2 + area3

  p1 =   area1                   / area
  p2 = ( area1 + area2 )         / area
!
!  Generate a batch of points in [0,1]x[0,1].
!
  call r8mat_uniform_01 ( m, n, seed, point )
!
!  Generate a batch of N probabilities.
!
  call r8vec_uniform_01 ( n, seed, prob )
!
!  Map some points into [X1,X2] x [Y1,Y4].
!
  where ( prob(1:n) <= p1 )

    point(1,1:n) = x1 + point(1,1:n) * ( x2 - x1 )
    point(2,1:n) = y1 + point(2,1:n) * ( y4 - y1 )
!
!  Map points into [X2,X3] x [Y2,Y3].
!
  elsewhere ( prob(1:n) <= p2 )

    point(1,1:n) = x2 + point(1,1:n) * ( x3 - x2 )
    point(2,1:n) = y2 + point(2,1:n) * ( y3 - y2 )
!
!  Map points into [X3,X4] x [Y1,Y4].
!
  elsewhere

    point(1,1:n) = x3 + point(1,1:n) * ( x4 - x3 )
    point(2,1:n) = y1 + point(2,1:n) * ( y4 - y1 )

  end where

  return
end
subroutine p12_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P12_SDIST returns the signed distance to the region in problem 12.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P12_SDIST - Fatal error!'
  write ( *, '(a)' ) '  This routine is not written yet!'

  stop
end
subroutine p12_title ( title )

!*****************************************************************************80
!
!! P12_TITLE returns a title for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2005
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

  title = '#12: John Shadid''s H-shaped region.'

  return
end
subroutine p13_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P13_BOUNDARY_NEAREST returns a nearest boundary point in problem 13.
!
!  Discussion:
!
!    The correct computation of the distance to the boundary of the
!    region in this problem is complicated for points in the
!    exterior which must choose between the vertical shaft and
!    the semicircular annulus.  For our purposes, it is not essential
!    to get this computation exactly, so we crudely draw a 45 degree
!    dividing line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle2
  real ( kind = 8 ) angle3
  real ( kind = 8 ), dimension ( m, n ) :: boundary
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) r
  real ( kind = 8 ) y2

  angle2 = acos ( 5.0D+00 / 40.0D+00 )
  angle3 = pi - angle2

  y2 = 40.0D+00 * sin ( angle2 )

  do j = 1, n

    r = sqrt ( ( point(1,j) - 50.0D+00 )**2 + ( point(2,j) - 0.0D+00 )**2 )
    angle = atan2 ( point(2,j), point(1,j) - 50.0D+00 )
!
!  1. Left side of left foot
!
    if ( point(1,j) <= 10.0D+00 .and. point(2,j) <= 0.0D+00 ) then

      boundary(1:2,j) = (/ 10.0D+00, 0.0D+00 /)
!
!  2. Left foot
!
    else if ( point(1,j) <= 20.0D+00 .and. point(2,j) <= 0.0D+00 ) then

      boundary(1:2,j) = (/ point(1,j), 0.0D+00 /)
!
!  3. Right side of left foot
!
    else if ( point(1,j) <= 50.0D+00 .and. point(2,j) <= 0.0D+00 ) then

      boundary(1:2,j) = (/ 20.0D+00, 0.0D+00 /)
!
!  4.  Left side of right foot.
!
    else if ( point(1,j) <= 80.0D+00 .and. point(2,j) <= 0.0D+00 ) then

      boundary(1:2,j) = (/ 80.0D+00, 0.0D+00 /)
!
!  5.  Right Foot
!
    else if ( point(1,j) <= 90.0D+00 .and. point(2,j) <= 0.0D+00 ) then

      boundary(1:2,j) = (/ point(1,j), 0.0D+00 /)
!
!  6.  Right side of right foot.
!
    else if ( point(2,j) <= 0.0D+00 ) then

      boundary(1:2,j) = (/ 90.0D+00, 0.0D+00 /)
!
!  7.  Space between legs.
!
    else if ( r <= 30.0D+00 ) then

      boundary(1:2,j) = (/ 50.0D+00 + 30.0D+00 * cos ( angle ), &
                            0.0D+00 + 30.0D+00 * sin ( angle ) /)
!
!  8.  Right inside of left leg.
!
    else if ( r <= 35.0D+00 .and. point(1,j) <= 45.0D+00 ) then

      boundary(1:2,j) = (/ 50.0D+00 + 30.0D+00 * cos ( angle ), &
                            0.0D+00 + 30.0D+00 * sin ( angle ) /)
!
!  9.  Crotch.
!
    else if ( r <= 35.0D+00 .and. point(1,j) <= 55.0D+00 ) then

      boundary(1:2,j) = (/ 50.0D+00 + 30.0D+00 * cos ( angle ), &
                            0.0D+00 + 30.0D+00 * sin ( angle ) /)
!
!  10.  Left inside of right leg.
!
    else if ( r <= 35.0D+00 .and. 55.0D+00 <= point(1,j) ) then

      boundary(1:2,j) = (/ 50.0D+00 + 30.0D+00 * cos ( angle ), &
                            0.0D+00 + 30.0D+00 * sin ( angle ) /)
!
!  11.  Left inside of left leg.
!
    else if ( r <= 40.0D+00 .and. point(1,j) <= 45.0D+00 ) then

      boundary(1:2,j) = (/ 50.0D+00 + 40.0D+00 * cos ( angle ), &
                            0.0D+00 + 40.0D+00 * sin ( angle ) /)
!
!  12.  Left hip.
!
    else if ( r <= 40.0D+00 .and. point(1,j) <= 50.0D+00 &
                            .and. point(2,j) <= y2 ) then

      boundary(1:2,j) = (/ 45.0D+00, y2 /)
!
!  13.  Right hip.
!
    else if ( r <= 40.0D+00 .and. point(1,j) <= 55.0D+00 &
                            .and. point(2,j) <= y2 ) then

      boundary(1:2,j) = (/ 55.0D+00, y2 /)
!
!  14.  Right inside of right leg.
!
    else if ( r <= 40.0D+00 ) then

      boundary(1:2,j) = (/ 50.0D+00 + 40.0D+00 * cos ( angle ), &
                            0.0D+00 + 40.0D+00 * sin ( angle ) /)
!
!  15.  Left outside of left leg.
!
    else if ( angle3 < angle ) then

      boundary(1:2,j) = (/ 50.0D+00 + 30.0D+00 * cos ( angle ), &
                            0.0D+00 + 30.0D+00 * sin ( angle ) /)
!
!  16.  Right outside of right leg.
!
    else if ( angle < angle2 ) then

      boundary(1:2,j) = (/ 50.0D+00 + 40.0D+00 * cos ( angle ), &
                            0.0D+00 + 40.0D+00 * sin ( angle ) /)
!
!  17.  Left outside of trunk.
!
    else if ( point(1,j) <= 45.0D+00 .and. point(2,j) <= 90.0D+00 ) then

      boundary(1:2,j) = (/ 45.0D+00, point(2,j) /)
!
!  18.  Left inside of trunk.
!
    else if ( point(1,j) <= 50.0D+00 .and. point(2,j) <= 90.0D+00 .and. &
      point(1,j) + point(2,j) <= 135.0D+00 ) then

      boundary(1:2,j) = (/ 45.0D+00, point(2,j) /)
!
!  22.5 Throat
!
    else if ( point(1,j) <= 50.0D+00 .and. point(2,j) <= 90.0D+00 .and. &
      135.0D+00 <= point(1,j) + point(2,j) ) then

        boundary(1:2,j) = (/ point(1,j), 90.0D+00 /)

    else if ( point(1,j) <= 55.0D+00 .and. point(2,j) <= 90.0D+00 .and. &
      35.0D+00 <= point(2,j) - point(1,j) ) then

        boundary(1:2,j) = (/ point(1,j), 90.0D+00 /)
!
!  19.  Right inside of trunk.
!
    else if ( point(1,j) <= 55.0D+00 .and. point(2,j) <= 90.0D+00 .and. &
      point(2,j) - point(1,j) <= 35.0D+00 ) then

      boundary(1:2,j) = (/ 55.0D+00, point(2,j) /)
!
!  20.  Right outside of trunk.
!
    else if ( 55.0D+00 <= point(1,j) .and. point(2,j) <= 90.0D+00 ) then

      boundary(1:2,j) = (/ 55.0D+00, point(2,j) /)
!
!  21.  Left shoulder.
!
    else if ( point(1,j) <= 45.0D+00 ) then

      boundary(1:2,j) = (/ 45.0D+00, 90.0D+00 /)
!
!  22.  Neck.
!
    else if ( point(1,j) <= 55.0D+00 ) then

      boundary(1:2,j) = (/ point(1,j), 90.0D+00 /)
!
!  23.  Right Shoulder.
!
    else if ( 55.0D+00 <= point(1,j) .and. 90.0D+00 <= point(2,j) ) then

      boundary(1:2,j) = (/ 55.0D+00, 90.0D+00 /)
!
!  Missing?
!
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P13_BOUNDARY_NEAREST - Fatal error!'
      write ( *, '(a)' ) '  No nearest boundary point found for point P.'
      write ( *, '(a,2g14.6)' ) '  P = ', point(1:2,j)
      stop
    end if

  end do

  return
end
subroutine p13_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P13_BOUNDARY_PROJECT projects exterior points to the boundary in problem 13.
!
!  Discussion:
!
!    The correct computation of the distance to the boundary of the
!    region in this problem is complicated for points in the
!    exterior which must choose between the vertical shaft and
!    the semicircular annulus.  For our purposes, it is not essential
!    to get this computation exactly, so we crudely draw a 45 degree
!    dividing line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, exterior points have been projected.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle2
  real ( kind = 8 ) angle3
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) r
  real ( kind = 8 ) y2

  angle2 = acos ( 5.0D+00 / 40.0D+00 )
  angle3 = pi - angle2

  y2 = 40.0D+00 * sin ( angle2 )

  do j = 1, n

    r = sqrt ( ( point(1,j) - 50.0D+00 )**2 + ( point(2,j) - 0.0D+00 )**2 )
    angle = atan2 ( point(2,j), point(1,j) - 50.0D+00 )
!
!  1. Left side of left foot
!
    if ( point(1,j) <= 10.0D+00 .and. point(2,j) <= 0.0D+00 ) then

      point(1:2,j) = (/ 10.0D+00, 0.0D+00 /)
!
!  2. Left foot
!
    else if ( point(1,j) <= 20.0D+00 .and. point(2,j) <= 0.0D+00 ) then

      point(1:2,j) = (/ point(1,j), 0.0D+00 /)
!
!  3. Right side of left foot
!
    else if ( point(1,j) <= 50.0D+00 .and. point(2,j) <= 0.0D+00 ) then

      point(1:2,j) = (/ 20.0D+00, 0.0D+00 /)
!
!  4.  Left side of right foot.
!
    else if ( point(1,j) <= 80.0D+00 .and. point(2,j) <= 0.0D+00 ) then

      point(1:2,j) = (/ 80.0D+00, 0.0D+00 /)
!
!  5.  Right Foot
!
    else if ( point(1,j) <= 90.0D+00 .and. point(2,j) <= 0.0D+00 ) then

      point(1:2,j) = (/ point(1,j), 0.0D+00 /)
!
!  6.  Right side of right foot.
!
    else if ( point(2,j) <= 0.0D+00 ) then

      point(1:2,j) = (/ 90.0D+00, 0.0D+00 /)
!
!  7.  Space between legs.
!
    else if ( r <= 30.0D+00 ) then

      point(1:2,j) = (/ 50.0D+00 + 30.0D+00 * cos ( angle ), &
                         0.0D+00 + 30.0D+00 * sin ( angle ) /)
!
!  8.  Right inside of left leg.
!
    else if ( r <= 35.0D+00 .and. point(1,j) <= 45.0D+00 ) then
!
!  9.  Crotch.
!
    else if ( r <= 35.0D+00 .and. point(1,j) <= 55.0D+00 ) then
!
!  10.  Left inside of right leg.
!
    else if ( r <= 35.0D+00 .and. 55.0D+00 <= point(1,j) ) then
!
!  11.  Left inside of left leg.
!
    else if ( r <= 40.0D+00 .and. point(1,j) <= 45.0D+00 ) then
!
!  12.  Left hip.
!
    else if ( r <= 40.0D+00 .and. point(1,j) <= 50.0D+00 &
                            .and. point(2,j) <= y2 ) then
!
!  13.  Right hip.
!
    else if ( r <= 40.0D+00 .and. point(1,j) <= 55.0D+00 &
                            .and. point(2,j) <= y2 ) then
!
!  14.  Right inside of right leg.
!
    else if ( r <= 40.0D+00 ) then
!
!  15.  Left outside of left leg.
!
    else if ( angle3 < angle ) then

      point(1:2,j) = (/ 50.0D+00 + 30.0D+00 * cos ( angle ), &
                         0.0D+00 + 30.0D+00 * sin ( angle ) /)
!
!  16.  Right outside of right leg.
!
    else if ( angle < angle2 ) then

      point(1:2,j) = (/ 50.0D+00 + 40.0D+00 * cos ( angle ), &
                         0.0D+00 + 40.0D+00 * sin ( angle ) /)
!
!  17.  Left outside of trunk.
!
    else if ( point(1,j) <= 45.0D+00 .and. point(2,j) <= 90.0D+00 ) then

      point(1:2,j) = (/ 45.0D+00, point(2,j) /)
!
!  18.  Left inside of trunk.
!
    else if ( point(1,j) <= 50.0D+00 .and. point(2,j) <= 90.0D+00 .and. &
      point(1,j) + point(2,j) <= 135.0D+00 ) then
!
!  22.5 Throat
!
    else if ( point(1,j) <= 50.0D+00 .and. point(2,j) <= 90.0D+00 .and. &
      135.0D+00 <= point(1,j) + point(2,j) ) then

    else if ( point(1,j) <= 55.0D+00 .and. point(2,j) <= 90.0D+00 .and. &
      35.0D+00 <= point(2,j) - point(1,j) ) then
!
!  19.  Right inside of trunk.
!
    else if ( point(1,j) <= 55.0D+00 .and. point(2,j) <= 90.0D+00 .and. &
      point(2,j) - point(1,j) <= 35.0D+00 ) then
!
!  20.  Right outside of trunk.
!
    else if ( 55.0D+00 <= point(1,j) .and. point(2,j) <= 90.0D+00 ) then

      point(1:2,j) = (/ 55.0D+00, point(2,j) /)
!
!  21.  Left shoulder.
!
    else if ( point(1,j) <= 45.0D+00 ) then

      point(1:2,j) = (/ 45.0D+00, 90.0D+00 /)
!
!  22.  Neck.
!
    else if ( point(1,j) <= 55.0D+00 ) then

      point(1:2,j) = (/ point(1,j), 90.0D+00 /)
!
!  23.  Right Shoulder.
!
    else if ( 55.0D+00 <= point(1,j) .and. 90.0D+00 <= point(2,j) ) then

      point(1:2,j) = (/ 55.0D+00, 90.0D+00 /)

    end if

  end do

  return
end
subroutine p13_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P13_BOUNDARY_SEGMENT returns a boundary segment in problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle1
  real ( kind = 8 ) angle2
  real ( kind = 8 ) arc_length
  real ( kind = 8 ) arc1
  real ( kind = 8 ) arc2
  real ( kind = 8 ) arc3
  real ( kind = 8 ) arc4
  real ( kind = 8 ) arc5
  real ( kind = 8 ) arc6
  real ( kind = 8 ) arc7
  real ( kind = 8 ) arc8
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n01
  integer ( kind = 4 ) n02
  integer ( kind = 4 ) n03
  integer ( kind = 4 ) n04
  integer ( kind = 4 ) n05
  integer ( kind = 4 ) n06
  integer ( kind = 4 ) n07
  integer ( kind = 4 ) n08
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), dimension ( 2, 8 ) :: s = reshape ( (/ &
    10.000D+00,  0.000D+00, &
    20.000D+00,  0.000D+00, &
    80.000D+00,  0.000D+00, &
    90.000D+00,  0.000D+00, &
    55.000D+00, 39.686268D+00, &
    55.000D+00, 90.000D+00, &
    45.000D+00, 90.000D+00, &
    45.000D+00, 39.686268D+00 /), (/ 2, 8 /) )
  real ( kind = 8 ) segment(m,segment_length)
  integer ( kind = 4 ) segment_index

  n = segment_length - 1

  angle = acos ( 5.0D+00 / 40.0D+00 )

  arc1 = 30.0D+00 * pi
  arc2 = 10.0D+00
  arc3 = 40.0D+00 * angle
  arc4 = 90.0D+00 - 40.0D+00 * sin ( angle )
  arc5 = 10.0D+00
  arc6 = 90.0D+00 - 40.0D+00 * sin ( angle )
  arc7 = 40.0D+00 * angle
  arc8 = 10.0D+00

  arc_length = arc1 + arc2 + arc3 + arc4 + arc5 + arc6 + arc7 + arc8

  n01 = nint ( arc1 &
    * real ( n - 1, kind = 8 ) / arc_length )

  n02 = nint ( ( arc1 + arc2 ) &
    * real ( n - 1, kind = 8 ) / arc_length ) &
    - n01

  n03 = nint ( ( arc1 + arc2 + arc3 ) &
    * real ( n - 1, kind = 8 ) / arc_length ) &
    - n01 - n02

  n04 = nint ( ( arc1 + arc2 + arc3 + arc4 ) &
    * real ( n - 1, kind = 8 ) / arc_length ) &
    - n01 - n02 - n03

  n05 = nint ( ( arc1 + arc2 + arc3 + arc4 + arc5 ) &
    * real ( n - 1, kind = 8 ) / arc_length ) &
    - n01 - n02 - n03 - n04

  n06 = nint ( ( arc1 + arc2 + arc3 + arc4 + arc5 + arc6 ) &
    * real ( n - 1, kind = 8 ) / arc_length ) &
    - n01 - n02 - n03 - n04 - n05

  n07 = nint ( ( arc1 + arc2 + arc3 + arc4 + arc5 + arc6 + arc7 ) &
    * real ( n - 1, kind = 8 ) / arc_length ) &
    - n01 - n02 - n03 - n04 - n05 - n06

  n08 = n - n01 - n02 - n03 - n04 - n05 - n06 - n07

  j = 0
!
!  1: Bottom semicircle.
!
  do i = 1, n01
    j = j + 1
    angle = real ( n01 - i + 1, kind = 8 ) * pi / real ( n01, kind = 8 )
    segment(1,j) = 50.0D+00 + 30.0D+00 * cos ( angle )
    segment(2,j) =  0.0D+00 + 30.0D+00 * sin ( angle )
  end do
!
!  2: Right foot.
!
  do i = 1, n02
    j = j + 1
    segment(1:2,j) = ( real ( n02 - i + 1, kind = 8 ) * s(1:2,3)   &
                     + real (       i - 1, kind = 8 ) * s(1:2,4) ) &
                     / real ( n02,         kind = 8 )
  end do
!
!  3: Top right pen-semicircle.
!
  angle1 = 0.0D+00
  angle2 = acos ( 5.0D+00 / 40.0D+00 )

  do i = 1, n03
    j = j + 1
    angle = real ( i - 1, kind = 8 ) * angle2 / real ( n03, kind = 8 )
    segment(1,j) = 50.0D+00 + 40.0D+00 * cos ( angle )
    segment(2,j) =  0.0D+00 + 40.0D+00 * sin ( angle )
  end do
!
!  4: Right side of shaft.
!
  do i = 1, n04
    j = j + 1
    segment(1:2,j) = ( real ( n04 - i + 1, kind = 8 ) * s(1:2,5)   &
                     + real (       i - 1, kind = 8 ) * s(1:2,6) ) &
                       / real ( n04,         kind = 8 )
  end do
!
!  5: Top of shaft.
!
  do i = 1, n05
    j = j + 1
    segment(1:2,j) = ( real ( n05 - i + 1, kind = 8 ) * s(1:2,6)   &
                     + real (       i - 1, kind = 8 ) * s(1:2,7) ) &
                     / real ( n05,         kind = 8 )
  end do
!
!  6: Left side of shaft.
!
  do i = 1, n06
    j = j + 1
    segment(1:2,j) = ( real ( n06 - i + 1, kind = 8 ) * s(1:2,7)   &
                     + real (       i - 1, kind = 8 ) * s(1:2,8) ) &
                     / real ( n06,         kind = 8 )
  end do
!
!  7: Left pen-semicircle.
!
  angle1 = pi - acos ( 5.0D+00 / 40.0D+00 )
  angle2 = pi

  do i = 1, n07
    j = j + 1
    angle = ( real ( n07 - i + 1, kind = 8 ) * angle1 + &
              real (       i - 1, kind = 8 ) * angle2 ) &
            / real ( n07,         kind = 8 )
    segment(1,j) = 50.0D+00 + 40.0D+00 * cos ( angle )
    segment(2,j) =  0.0D+00 + 40.0D+00 * sin ( angle )
  end do
!
!  8: Left foot.
!
  do i = 1, n08
    j = j + 1
    segment(1:2,j) = ( real ( n08 - i + 1, kind = 8 ) * s(1:2,1)   &
                     + real (       i - 1, kind = 8 ) * s(1:2,2) ) &
                     / real ( n08,         kind = 8 )
  end do

  j = j + 1
  segment(1:2,j) = s(1:2,2)

  return
end
subroutine p13_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P13_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) arc_length
  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P13_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  angle = acos ( 5.0D+00 / 40.0D+00 )

  arc_length = &
      30.0D+00 * pi &
    + 10.0D+00 &
    + 40.0D+00 * angle &
    + 90.0D+00 - 40.0D+00 * sin ( angle ) &
    + 10.0D+00 &
    + 90.0D+00 - 40.0D+00 * sin ( angle ) &
    + 40.0D+00 * angle &
    + 10.0D+00

  n = nint ( arc_length / h )
  n = max ( n, 8 )

  segment_length = n

  return
end
subroutine p13_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P13_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 1

  return
end
subroutine p13_box ( m, lo, hi )

!*****************************************************************************80
!
!! P13_BOX returns a bounding box for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)

  lo(1:m) = (/    0.0D+00,    0.0D+00 /)
  hi(1:m) = (/ +100.0D+00, +100.0D+00 /)

  return
end
subroutine p13_density ( m, n, point, density )

!*****************************************************************************80
!
!! P13_DENSITY returns the density for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p13_element_size ( element_size )

!*****************************************************************************80
!
!! P13_ELEMENT_SIZE returns a typical element size for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.15D+00

  return
end
subroutine p13_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P13_FIXED_NUM returns the number of fixed points in problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 8

  return
end
subroutine p13_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P13_FIXED_POINTS returns the fixed points in problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) fixed(m,fixed_num)

  fixed(1:m,1:fixed_num) = reshape ( (/ &
   10.000D+00,  0.000D+00, &
   20.000D+00,  0.000D+00, &
   80.000D+00,  0.000D+00, &
   90.000D+00,  0.000D+00, &
   55.000D+00, 39.686268D+00, &
   55.000D+00, 90.000D+00, &
   45.000D+00, 90.000D+00, &
   45.000D+00, 39.686268D+00 /), (/ m, fixed_num /) )

  return
end
subroutine p13_header ( )

!*****************************************************************************80
!
!! P13_HEADER prints some information about problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num

  call p13_boundary_segment_num ( boundary_segment_num )
  call p13_fixed_num ( fixed_num )
  call p13_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P13:'
  write ( *, '(a)' ) '  The Sandia Fork region.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p13_hole_num ( hole_num )

!*****************************************************************************80
!
!! P13_HOLE_NUM counts the holes in problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 0

  return
end
subroutine p13_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P13_HOLE_POINT returns a point inside a given hole in problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  return
end
subroutine p13_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P13_INSIDE reports if a point is inside the region in problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical inside(n)
  real ( kind = 8 ) point(m,n)

  inside(1:n) =                                                       &
       (45.0D+00 <= point(1,1:n) .and. point(1,1:n) <= 55.0D+00 .and. &
        30.0D+00 <= point(2,1:n) .and. point(2,1:n) <= 90.0D+00    )  &
    .or.                                                              &
       ( 900.0D+00 <=                                                 &
         ( point(1,1:n) - 50.0D+00 )**2                               &
       + ( point(2,1:n) -  0.0D+00 )**2 .and.                         &
         ( point(1,1:n) - 50.0D+00 )**2                               &
       + ( point(2,1:n) -  0.0D+00 )**2 <= 1600.0D+00 )

  return
end
subroutine p13_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P13_SAMPLE samples points from the region in problem 13.
!
!  Discussion:
!
!    The region is contained in the box [0,100] x [0,100].
!
!    The region looks roughly like an inverted tuning fork.
!    It is the union of a rectangular strip and a partial
!    circular annulus.
!
!    If three dimensions are used, then the 2D region is simply
!    projected through the range 0 <= Z <= 20.0.
!
!    The information about this region was supplied by David
!    Crawford of Sandia National Laboratory.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) reject
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension(3) :: width = (/ &
    100.0D+00, 100.0D+00, 20.0D+00 /)
  real ( kind = 8 ) x(m)

  reject = 0

  do j = 1, n

    do

      call r8vec_uniform_01 ( m, seed, x )

      x(1:m) = width(1:m) * x(1:m)
!
!  The rectangular strip.
!
      if ( &
        45.0D+00 <= x(1)             .and. &
                    x(1) <= 55.0D+00 .and. &
        30.0D+00 <= x(2)             .and. &
                    x(2) <= 90.0D+00 ) then
        exit
!
!  The annulus.
!
      else if ( 0.0D+00 <= x(2) ) then

        c = sqrt ( ( x(1) - 50.0D+00 ) * ( x(1) - 50.0D+00 ) + x(2) * x(2) )

        if ( 30.0D+00 <= c .and. c <= 40.0D+00 ) then
          exit
        end if

      end if

      reject = reject + 1

      if ( 30 * n + 10 <= reject ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P13_SAMPLE - Fatal error!'
        write ( *, '(a,i8)' ) '  Trying to generate point J = ', j
        write ( *, '(a,i8)' ) '  Number of rejections = ', reject
        write ( *, '(a,g14.6)' ) '  Rejection percentage = ', &
          real ( 100 * reject ) / real ( reject + j - 1 )
          call r8vec_print ( m, x, '  Most recent rejected point: ' )
        stop
      end if

    end do

    point(1:m,j) = x(1:m)

  end do

  return
end
subroutine p13_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P13_SAMPLE_H1 samples points from the enlarged region in problem 13.
!
!  Discussion:
!
!    We enlarge the region by a layer H.  We do not round the
!    corners of the region, which should be done if we literally only
!    want to add points within H units of the boundary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) h
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) reject
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension(3) :: width = (/ &
    100.0D+00, 100.0D+00, 20.0D+00 /)
  real ( kind = 8 ) x(m)

  reject = 0

  do j = 1, n

    do

      call r8vec_uniform_01 ( m, seed, x )

      x(1:m) = width(1:m) * x(1:m)
!
!  The rectangular strip.
!
      if ( &
        45.0D+00 - h <= x(1)                 .and. &
                        x(1) <= 55.0D+00 + h .and. &
        30.0D+00 - h <= x(2)                 .and. &
                        x(2) <= 90.0D+00 + h ) then
        exit
!
!  The left "footprint" of the annulus.
!
    else if ( 10.0D+00 - h <= x(1)                 .and. &
                              x(1) <= 20.0D+00 + h .and. &
                        -h <= x(2)                 .and. &
                              x(2) <= 0.0D+00 ) then

      exit
!
!  The right "footprint" of the annulus.
!
    else if ( 80.0D+00 - h <= x(1)                 .and. &
                              x(1) <= 90.0D+00 + h .and. &
                        -h <= x(2)                 .and. &
                              x(2) <= 0.0D+00 ) then

      exit
!
!  The annulus.
!
      else if ( 0.0D+00 <= x(2) ) then

        c = sqrt ( ( x(1) - 50.0D+00 ) * ( x(1) - 50.0D+00 ) + x(2) * x(2) )

        if ( 30.0D+00 - h <= c .and. c <= 40.0D+00 + h ) then
          exit
        end if

      end if

      reject = reject + 1

      if ( 30 * n + 10 <= reject ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P13_SAMPLE_H1 - Fatal error!'
        write ( *, '(a,i8)' ) '  Trying to generate point J = ', j
        write ( *, '(a,i8)' ) '  Number of rejections = ', reject
        write ( *, '(a,g14.6)' ) '  Rejection percentage = ', &
          real ( 100 * reject ) / real ( reject + j - 1 )
          call r8vec_print ( m, x, '  Most recent rejected point: ' )
        stop
      end if

    end do

    point(1:m,j) = x(1:m)

  end do

  return
end
subroutine p13_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P13_SDIST returns the signed distance to the region in problem 13.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P13_SDIST - Fatal error!'
  write ( *, '(a)' ) '  This routine is not written yet!'

  stop
end
subroutine p13_title ( title )

!*****************************************************************************80
!
!! P13_TITLE returns a title for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2005
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

  title = '#13: The Sandia Fork region.'

  return
end
subroutine p14_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P14_BOUNDARY_NEAREST returns a nearest boundary point in problem 14.
!
!  Discussion:
!
!    The correct computation of the distance to the boundary of the
!    region in this problem is complicated for points in the
!    exterior which must choose between the vertical shaft and
!    the semicircular annulus.  For our purposes, it is not essential
!    to get this computation exactly, so we crudely draw a 45 degree
!    dividing line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: n1 = 90
  integer ( kind = 4 ), parameter :: n2 = 11

  real ( kind = 8 ) boundary(m,n)
  real ( kind = 8 ) dist1
  real ( kind = 8 ) dist2
  integer ( kind = 4 ) j
  real ( kind = 8 ) pn1(m)
  real ( kind = 8 ) pn2(m)
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), dimension ( 2, n1 ) :: v1 = reshape ( (/ &
    316.43027D+00, 404.47559D+00, &
    291.04946D+00, 400.70917D+00, &
    265.16504D+00, 409.77890D+00, &
    241.46794D+00, 402.40310D+00, &
    216.55145D+00, 396.52064D+00, &
    163.28492D+00, 411.37102D+00, &
    142.81752D+00, 391.16355D+00, &
    111.95404D+00, 346.70264D+00, &
    100.03538D+00, 325.72710D+00, &
    103.98723D+00, 302.51587D+00, &
    128.72978D+00, 285.72802D+00, &
    147.49111D+00, 266.23345D+00, &
    196.65261D+00, 242.24055D+00, &
    213.56835D+00, 221.67192D+00, &
    226.49969D+00, 198.09326D+00, &
    248.37126D+00, 183.50473D+00, &
    262.21952D+00, 165.39102D+00, &
    278.42330D+00, 149.91715D+00, &
    300.71846D+00, 145.82601D+00, &
    311.12698D+00, 166.71094D+00, &
    326.66315D+00, 184.58335D+00, &
    359.78574D+00, 225.48049D+00, &
    357.08892D+00, 252.88958D+00, &
    358.76685D+00, 285.34403D+00, &
    361.50834D+00, 303.71287D+00, &
    371.68926D+00, 314.92452D+00, &
    380.49890D+00, 324.58632D+00, &
    396.37634D+00, 328.88990D+00, &
    412.59116D+00, 327.25238D+00, &
    425.48394D+00, 315.28623D+00, &
    435.84305D+00, 302.44664D+00, &
    458.34025D+00, 297.55121D+00, &
    479.66439D+00, 288.99238D+00, &
    493.09812D+00, 270.20636D+00, &
    518.87309D+00, 264.56427D+00, &
    547.18014D+00, 268.18846D+00, &
    600.49708D+00, 240.62570D+00, &
    625.96183D+00, 238.40347D+00, &
    633.90530D+00, 260.70629D+00, &
    621.50451D+00, 285.88914D+00, &
    576.87224D+00, 322.14121D+00, &
    570.51915D+00, 348.85423D+00, &
    567.16400D+00, 378.24075D+00, &
    558.00668D+00, 406.86552D+00, &
    565.19008D+00, 435.75599D+00, &
    567.56437D+00, 465.33407D+00, &
    550.87626D+00, 490.96358D+00, &
    532.98174D+00, 515.84491D+00, &
    500.66817D+00, 551.89078D+00, &
    478.75120D+00, 562.17222D+00, &
    430.03371D+00, 583.94286D+00, &
    401.20454D+00, 587.69910D+00, &
    368.32214D+00, 581.10110D+00, &
    354.26303D+00, 585.86085D+00, &
    346.75200D+00, 601.10367D+00, &
    332.85137D+00, 628.74602D+00, &
    308.02188D+00, 645.84180D+00, &
    295.52344D+00, 647.18525D+00, &
    286.51519D+00, 651.60328D+00, &
    285.98846D+00, 662.07339D+00, &
    298.93455D+00, 665.66316D+00, &
    301.70226D+00, 682.79570D+00, &
    278.65857D+00, 689.63850D+00, &
    266.25737D+00, 712.11005D+00, &
    287.28701D+00, 732.77147D+00, &
    318.19548D+00, 736.85151D+00, &
    343.83067D+00, 753.60957D+00, &
    375.53164D+00, 758.35231D+00, &
    405.73444D+00, 768.98687D+00, &
    406.33873D+00, 785.59001D+00, &
    378.35436D+00, 789.44240D+00, &
    350.02151D+00, 795.02238D+00, &
    338.68030D+00, 788.87325D+00, &
    325.67930D+00, 786.10177D+00, &
    319.05995D+00, 798.04657D+00, &
    301.78158D+00, 795.34254D+00, &
    280.69272D+00, 773.86634D+00, &
    254.55844D+00, 758.02898D+00, &
    234.07759D+00, 737.42090D+00, &
    218.38337D+00, 711.41500D+00, &
    220.99086D+00, 682.17833D+00, &
    224.50640D+00, 651.96297D+00, &
    240.25971D+00, 631.36117D+00, &
    259.86174D+00, 612.60253D+00, &
    291.85381D+00, 556.70385D+00, &
    315.52139D+00, 537.56387D+00, &
    341.63663D+00, 520.12519D+00, &
    351.37130D+00, 458.75372D+00, &
    349.33183D+00, 431.31454D+00, &
    328.80465D+00, 412.43055D+00 /), (/ 2, n1 /) )
  real ( kind = 8 ), dimension ( 2, n2 ) :: v2 = reshape ( (/ &
    238.64853D+00, 266.58978D+00, &
    235.14026D+00, 287.95183D+00, &
    238.20736D+00, 303.46785D+00, &
    250.13902D+00, 303.71290D+00, &
    258.51675D+00, 297.46973D+00, &
    274.55300D+00, 291.27357D+00, &
    284.66230D+00, 280.72063D+00, &
    279.73288D+00, 267.83455D+00, &
    270.68478D+00, 255.55440D+00, &
    255.73801D+00, 249.16872D+00, &
    241.72690D+00, 256.73448D+00 /), (/ 2, n2 /) )

  do j = 1, n

    call polygon_point_near_2d ( n1, v1, point(1:m,j), pn1, dist1 )

    call polygon_point_near_2d ( n2, v2, point(1:m,j), pn2, dist2 )

    if ( dist1 < dist2 ) then
      boundary(1:m,j) = pn1(1:m)
    else
      boundary(1:m,j) = pn2(1:m)
    end if

  end do

  return
end
subroutine p14_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P14_BOUNDARY_PROJECT projects exterior points to the boundary in problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, exterior points have been projected.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical inside(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) pn(m)
  real ( kind = 8 ), dimension ( m, n ) :: point

  call p14_inside ( m, n, point, inside )

  do j = 1, n

    if ( .not. inside(j) ) then
      call p14_boundary_nearest ( m, 1, point(1:m,j), pn )
      point(1:m,j) = pn(1:m)
    end if

  end do

  return
end
subroutine p14_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P14_BOUNDARY_SEGMENT returns a boundary segment in problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: n1 = 90
  integer ( kind = 4 ), parameter :: n2 = 11
  integer ( kind = 4 ) segment_length

  real ( kind = 8 ) segment(m,segment_length)
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ), dimension ( 2, n1 ) :: v1 = reshape ( (/ &
    316.43027D+00, 404.47559D+00, &
    291.04946D+00, 400.70917D+00, &
    265.16504D+00, 409.77890D+00, &
    241.46794D+00, 402.40310D+00, &
    216.55145D+00, 396.52064D+00, &
    163.28492D+00, 411.37102D+00, &
    142.81752D+00, 391.16355D+00, &
    111.95404D+00, 346.70264D+00, &
    100.03538D+00, 325.72710D+00, &
    103.98723D+00, 302.51587D+00, &
    128.72978D+00, 285.72802D+00, &
    147.49111D+00, 266.23345D+00, &
    196.65261D+00, 242.24055D+00, &
    213.56835D+00, 221.67192D+00, &
    226.49969D+00, 198.09326D+00, &
    248.37126D+00, 183.50473D+00, &
    262.21952D+00, 165.39102D+00, &
    278.42330D+00, 149.91715D+00, &
    300.71846D+00, 145.82601D+00, &
    311.12698D+00, 166.71094D+00, &
    326.66315D+00, 184.58335D+00, &
    359.78574D+00, 225.48049D+00, &
    357.08892D+00, 252.88958D+00, &
    358.76685D+00, 285.34403D+00, &
    361.50834D+00, 303.71287D+00, &
    371.68926D+00, 314.92452D+00, &
    380.49890D+00, 324.58632D+00, &
    396.37634D+00, 328.88990D+00, &
    412.59116D+00, 327.25238D+00, &
    425.48394D+00, 315.28623D+00, &
    435.84305D+00, 302.44664D+00, &
    458.34025D+00, 297.55121D+00, &
    479.66439D+00, 288.99238D+00, &
    493.09812D+00, 270.20636D+00, &
    518.87309D+00, 264.56427D+00, &
    547.18014D+00, 268.18846D+00, &
    600.49708D+00, 240.62570D+00, &
    625.96183D+00, 238.40347D+00, &
    633.90530D+00, 260.70629D+00, &
    621.50451D+00, 285.88914D+00, &
    576.87224D+00, 322.14121D+00, &
    570.51915D+00, 348.85423D+00, &
    567.16400D+00, 378.24075D+00, &
    558.00668D+00, 406.86552D+00, &
    565.19008D+00, 435.75599D+00, &
    567.56437D+00, 465.33407D+00, &
    550.87626D+00, 490.96358D+00, &
    532.98174D+00, 515.84491D+00, &
    500.66817D+00, 551.89078D+00, &
    478.75120D+00, 562.17222D+00, &
    430.03371D+00, 583.94286D+00, &
    401.20454D+00, 587.69910D+00, &
    368.32214D+00, 581.10110D+00, &
    354.26303D+00, 585.86085D+00, &
    346.75200D+00, 601.10367D+00, &
    332.85137D+00, 628.74602D+00, &
    308.02188D+00, 645.84180D+00, &
    295.52344D+00, 647.18525D+00, &
    286.51519D+00, 651.60328D+00, &
    285.98846D+00, 662.07339D+00, &
    298.93455D+00, 665.66316D+00, &
    301.70226D+00, 682.79570D+00, &
    278.65857D+00, 689.63850D+00, &
    266.25737D+00, 712.11005D+00, &
    287.28701D+00, 732.77147D+00, &
    318.19548D+00, 736.85151D+00, &
    343.83067D+00, 753.60957D+00, &
    375.53164D+00, 758.35231D+00, &
    405.73444D+00, 768.98687D+00, &
    406.33873D+00, 785.59001D+00, &
    378.35436D+00, 789.44240D+00, &
    350.02151D+00, 795.02238D+00, &
    338.68030D+00, 788.87325D+00, &
    325.67930D+00, 786.10177D+00, &
    319.05995D+00, 798.04657D+00, &
    301.78158D+00, 795.34254D+00, &
    280.69272D+00, 773.86634D+00, &
    254.55844D+00, 758.02898D+00, &
    234.07759D+00, 737.42090D+00, &
    218.38337D+00, 711.41500D+00, &
    220.99086D+00, 682.17833D+00, &
    224.50640D+00, 651.96297D+00, &
    240.25971D+00, 631.36117D+00, &
    259.86174D+00, 612.60253D+00, &
    291.85381D+00, 556.70385D+00, &
    315.52139D+00, 537.56387D+00, &
    341.63663D+00, 520.12519D+00, &
    351.37130D+00, 458.75372D+00, &
    349.33183D+00, 431.31454D+00, &
    328.80465D+00, 412.43055D+00 /), (/ 2, n1 /) )
  real ( kind = 8 ), dimension ( 2, n2 ) :: v2 = reshape ( (/ &
    238.64853D+00, 266.58978D+00, &
    235.14026D+00, 287.95183D+00, &
    238.20736D+00, 303.46785D+00, &
    250.13902D+00, 303.71290D+00, &
    258.51675D+00, 297.46973D+00, &
    274.55300D+00, 291.27357D+00, &
    284.66230D+00, 280.72063D+00, &
    279.73288D+00, 267.83455D+00, &
    270.68478D+00, 255.55440D+00, &
    255.73801D+00, 249.16872D+00, &
    241.72690D+00, 256.73448D+00 /), (/ 2, n2 /) )

  if ( segment_index == 1 ) then
    call polyloop_points_nd ( m, n1, v1, segment_length, segment )
  else if ( segment_index == 2 ) then
    call polyloop_points_nd ( m, n2, v2, segment_length, segment )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P14_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop
  end if

  return
end
subroutine p14_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P14_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 90
  integer ( kind = 4 ), parameter :: n2 = 11

  real ( kind = 8 ) h
  real ( kind = 8 ) length
  integer ( kind = 4 ) n
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length
  real ( kind = 8 ), dimension ( 2, n1 ) :: v1 = reshape ( (/ &
    316.43027D+00, 404.47559D+00, &
    291.04946D+00, 400.70917D+00, &
    265.16504D+00, 409.77890D+00, &
    241.46794D+00, 402.40310D+00, &
    216.55145D+00, 396.52064D+00, &
    163.28492D+00, 411.37102D+00, &
    142.81752D+00, 391.16355D+00, &
    111.95404D+00, 346.70264D+00, &
    100.03538D+00, 325.72710D+00, &
    103.98723D+00, 302.51587D+00, &
    128.72978D+00, 285.72802D+00, &
    147.49111D+00, 266.23345D+00, &
    196.65261D+00, 242.24055D+00, &
    213.56835D+00, 221.67192D+00, &
    226.49969D+00, 198.09326D+00, &
    248.37126D+00, 183.50473D+00, &
    262.21952D+00, 165.39102D+00, &
    278.42330D+00, 149.91715D+00, &
    300.71846D+00, 145.82601D+00, &
    311.12698D+00, 166.71094D+00, &
    326.66315D+00, 184.58335D+00, &
    359.78574D+00, 225.48049D+00, &
    357.08892D+00, 252.88958D+00, &
    358.76685D+00, 285.34403D+00, &
    361.50834D+00, 303.71287D+00, &
    371.68926D+00, 314.92452D+00, &
    380.49890D+00, 324.58632D+00, &
    396.37634D+00, 328.88990D+00, &
    412.59116D+00, 327.25238D+00, &
    425.48394D+00, 315.28623D+00, &
    435.84305D+00, 302.44664D+00, &
    458.34025D+00, 297.55121D+00, &
    479.66439D+00, 288.99238D+00, &
    493.09812D+00, 270.20636D+00, &
    518.87309D+00, 264.56427D+00, &
    547.18014D+00, 268.18846D+00, &
    600.49708D+00, 240.62570D+00, &
    625.96183D+00, 238.40347D+00, &
    633.90530D+00, 260.70629D+00, &
    621.50451D+00, 285.88914D+00, &
    576.87224D+00, 322.14121D+00, &
    570.51915D+00, 348.85423D+00, &
    567.16400D+00, 378.24075D+00, &
    558.00668D+00, 406.86552D+00, &
    565.19008D+00, 435.75599D+00, &
    567.56437D+00, 465.33407D+00, &
    550.87626D+00, 490.96358D+00, &
    532.98174D+00, 515.84491D+00, &
    500.66817D+00, 551.89078D+00, &
    478.75120D+00, 562.17222D+00, &
    430.03371D+00, 583.94286D+00, &
    401.20454D+00, 587.69910D+00, &
    368.32214D+00, 581.10110D+00, &
    354.26303D+00, 585.86085D+00, &
    346.75200D+00, 601.10367D+00, &
    332.85137D+00, 628.74602D+00, &
    308.02188D+00, 645.84180D+00, &
    295.52344D+00, 647.18525D+00, &
    286.51519D+00, 651.60328D+00, &
    285.98846D+00, 662.07339D+00, &
    298.93455D+00, 665.66316D+00, &
    301.70226D+00, 682.79570D+00, &
    278.65857D+00, 689.63850D+00, &
    266.25737D+00, 712.11005D+00, &
    287.28701D+00, 732.77147D+00, &
    318.19548D+00, 736.85151D+00, &
    343.83067D+00, 753.60957D+00, &
    375.53164D+00, 758.35231D+00, &
    405.73444D+00, 768.98687D+00, &
    406.33873D+00, 785.59001D+00, &
    378.35436D+00, 789.44240D+00, &
    350.02151D+00, 795.02238D+00, &
    338.68030D+00, 788.87325D+00, &
    325.67930D+00, 786.10177D+00, &
    319.05995D+00, 798.04657D+00, &
    301.78158D+00, 795.34254D+00, &
    280.69272D+00, 773.86634D+00, &
    254.55844D+00, 758.02898D+00, &
    234.07759D+00, 737.42090D+00, &
    218.38337D+00, 711.41500D+00, &
    220.99086D+00, 682.17833D+00, &
    224.50640D+00, 651.96297D+00, &
    240.25971D+00, 631.36117D+00, &
    259.86174D+00, 612.60253D+00, &
    291.85381D+00, 556.70385D+00, &
    315.52139D+00, 537.56387D+00, &
    341.63663D+00, 520.12519D+00, &
    351.37130D+00, 458.75372D+00, &
    349.33183D+00, 431.31454D+00, &
    328.80465D+00, 412.43055D+00 /), (/ 2, n1 /) )
  real ( kind = 8 ), dimension ( 2, n2 ) :: v2 = reshape ( (/ &
    238.64853D+00, 266.58978D+00, &
    235.14026D+00, 287.95183D+00, &
    238.20736D+00, 303.46785D+00, &
    250.13902D+00, 303.71290D+00, &
    258.51675D+00, 297.46973D+00, &
    274.55300D+00, 291.27357D+00, &
    284.66230D+00, 280.72063D+00, &
    279.73288D+00, 267.83455D+00, &
    270.68478D+00, 255.55440D+00, &
    255.73801D+00, 249.16872D+00, &
    241.72690D+00, 256.73448D+00 /), (/ 2, n2 /) )

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P14_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then
    call polyloop_length_nd ( 2, n1, v1, length )
  else if ( segment_index == 2 ) then
    call polyloop_length_nd ( 2, n2, v2, length )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P14_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal segment index = ', segment_index
    stop
  end if

  n = nint ( length / h )

  segment_length = n

  return
end
subroutine p14_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P14_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 2

  return
end
subroutine p14_box ( m, lo, hi )

!*****************************************************************************80
!
!! P14_BOX returns a bounding box for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)

  lo(1:m) = (/  100.0D+00,  145.0D+00 /)
  hi(1:m) = (/ +634.0D+00, +799.0D+00 /)

  return
end
subroutine p14_density ( m, n, point, density )

!*****************************************************************************80
!
!! P14_DENSITY returns the density for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p14_element_size ( element_size )

!*****************************************************************************80
!
!! P14_ELEMENT_SIZE returns a typical element size for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 1.0D+00

  return
end
subroutine p14_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P14_FIXED_NUM returns the number of fixed points in problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 101

  return
end
subroutine p14_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P14_FIXED_POINTS returns the fixed points in problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) fixed(m,fixed_num)

  fixed(1:m,1:fixed_num) = reshape ( (/ &
    316.43027D+00, 404.47559D+00, &
    291.04946D+00, 400.70917D+00, &
    265.16504D+00, 409.77890D+00, &
    241.46794D+00, 402.40310D+00, &
    216.55145D+00, 396.52064D+00, &
    163.28492D+00, 411.37102D+00, &
    142.81752D+00, 391.16355D+00, &
    111.95404D+00, 346.70264D+00, &
    100.03538D+00, 325.72710D+00, &
    103.98723D+00, 302.51587D+00, &
    128.72978D+00, 285.72802D+00, &
    147.49111D+00, 266.23345D+00, &
    196.65261D+00, 242.24055D+00, &
    213.56835D+00, 221.67192D+00, &
    226.49969D+00, 198.09326D+00, &
    248.37126D+00, 183.50473D+00, &
    262.21952D+00, 165.39102D+00, &
    278.42330D+00, 149.91715D+00, &
    300.71846D+00, 145.82601D+00, &
    311.12698D+00, 166.71094D+00, &
    326.66315D+00, 184.58335D+00, &
    359.78574D+00, 225.48049D+00, &
    357.08892D+00, 252.88958D+00, &
    358.76685D+00, 285.34403D+00, &
    361.50834D+00, 303.71287D+00, &
    371.68926D+00, 314.92452D+00, &
    380.49890D+00, 324.58632D+00, &
    396.37634D+00, 328.88990D+00, &
    412.59116D+00, 327.25238D+00, &
    425.48394D+00, 315.28623D+00, &
    435.84305D+00, 302.44664D+00, &
    458.34025D+00, 297.55121D+00, &
    479.66439D+00, 288.99238D+00, &
    493.09812D+00, 270.20636D+00, &
    518.87309D+00, 264.56427D+00, &
    547.18014D+00, 268.18846D+00, &
    600.49708D+00, 240.62570D+00, &
    625.96183D+00, 238.40347D+00, &
    633.90530D+00, 260.70629D+00, &
    621.50451D+00, 285.88914D+00, &
    576.87224D+00, 322.14121D+00, &
    570.51915D+00, 348.85423D+00, &
    567.16400D+00, 378.24075D+00, &
    558.00668D+00, 406.86552D+00, &
    565.19008D+00, 435.75599D+00, &
    567.56437D+00, 465.33407D+00, &
    550.87626D+00, 490.96358D+00, &
    532.98174D+00, 515.84491D+00, &
    500.66817D+00, 551.89078D+00, &
    478.75120D+00, 562.17222D+00, &
    430.03371D+00, 583.94286D+00, &
    401.20454D+00, 587.69910D+00, &
    368.32214D+00, 581.10110D+00, &
    354.26303D+00, 585.86085D+00, &
    346.75200D+00, 601.10367D+00, &
    332.85137D+00, 628.74602D+00, &
    308.02188D+00, 645.84180D+00, &
    295.52344D+00, 647.18525D+00, &
    286.51519D+00, 651.60328D+00, &
    285.98846D+00, 662.07339D+00, &
    298.93455D+00, 665.66316D+00, &
    301.70226D+00, 682.79570D+00, &
    278.65857D+00, 689.63850D+00, &
    266.25737D+00, 712.11005D+00, &
    287.28701D+00, 732.77147D+00, &
    318.19548D+00, 736.85151D+00, &
    343.83067D+00, 753.60957D+00, &
    375.53164D+00, 758.35231D+00, &
    405.73444D+00, 768.98687D+00, &
    406.33873D+00, 785.59001D+00, &
    378.35436D+00, 789.44240D+00, &
    350.02151D+00, 795.02238D+00, &
    338.68030D+00, 788.87325D+00, &
    325.67930D+00, 786.10177D+00, &
    319.05995D+00, 798.04657D+00, &
    301.78158D+00, 795.34254D+00, &
    280.69272D+00, 773.86634D+00, &
    254.55844D+00, 758.02898D+00, &
    234.07759D+00, 737.42090D+00, &
    218.38337D+00, 711.41500D+00, &
    220.99086D+00, 682.17833D+00, &
    224.50640D+00, 651.96297D+00, &
    240.25971D+00, 631.36117D+00, &
    259.86174D+00, 612.60253D+00, &
    291.85381D+00, 556.70385D+00, &
    315.52139D+00, 537.56387D+00, &
    341.63663D+00, 520.12519D+00, &
    351.37130D+00, 458.75372D+00, &
    349.33183D+00, 431.31454D+00, &
    328.80465D+00, 412.43055D+00, &
    238.64853D+00, 266.58978D+00, &
    235.14026D+00, 287.95183D+00, &
    238.20736D+00, 303.46785D+00, &
    250.13902D+00, 303.71290D+00, &
    258.51675D+00, 297.46973D+00, &
    274.55300D+00, 291.27357D+00, &
    284.66230D+00, 280.72063D+00, &
    279.73288D+00, 267.83455D+00, &
    270.68478D+00, 255.55440D+00, &
    255.73801D+00, 249.16872D+00, &
    241.72690D+00, 256.73448D+00 /), (/ m, fixed_num /) )

  return
end
subroutine p14_header ( )

!*****************************************************************************80
!
!! P14_HEADER prints some information about problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num

  call p14_boundary_segment_num ( boundary_segment_num )
  call p14_fixed_num ( fixed_num )
  call p14_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P14:'
  write ( *, '(a)' ) '  Marcus Garvie''s Lake Alpha and Beta island.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p14_hole_num ( hole_num )

!*****************************************************************************80
!
!! P14_HOLE_NUM counts the holes in problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 1

  return
end
subroutine p14_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P14_HOLE_POINT returns a point inside a given hole in problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  hole_point(1:m) = (/ 257.068D+00, 278.225D+00 /)

  return
end
subroutine p14_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P14_INSIDE reports if a point is inside the region in problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: n1 = 90
  integer ( kind = 4 ), parameter :: n2 = 11

  logical inside(n)
  logical inside1
  logical inside2
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ), dimension ( 2, n1 ) :: v1 = reshape ( (/ &
    316.43027D+00, 404.47559D+00, &
    291.04946D+00, 400.70917D+00, &
    265.16504D+00, 409.77890D+00, &
    241.46794D+00, 402.40310D+00, &
    216.55145D+00, 396.52064D+00, &
    163.28492D+00, 411.37102D+00, &
    142.81752D+00, 391.16355D+00, &
    111.95404D+00, 346.70264D+00, &
    100.03538D+00, 325.72710D+00, &
    103.98723D+00, 302.51587D+00, &
    128.72978D+00, 285.72802D+00, &
    147.49111D+00, 266.23345D+00, &
    196.65261D+00, 242.24055D+00, &
    213.56835D+00, 221.67192D+00, &
    226.49969D+00, 198.09326D+00, &
    248.37126D+00, 183.50473D+00, &
    262.21952D+00, 165.39102D+00, &
    278.42330D+00, 149.91715D+00, &
    300.71846D+00, 145.82601D+00, &
    311.12698D+00, 166.71094D+00, &
    326.66315D+00, 184.58335D+00, &
    359.78574D+00, 225.48049D+00, &
    357.08892D+00, 252.88958D+00, &
    358.76685D+00, 285.34403D+00, &
    361.50834D+00, 303.71287D+00, &
    371.68926D+00, 314.92452D+00, &
    380.49890D+00, 324.58632D+00, &
    396.37634D+00, 328.88990D+00, &
    412.59116D+00, 327.25238D+00, &
    425.48394D+00, 315.28623D+00, &
    435.84305D+00, 302.44664D+00, &
    458.34025D+00, 297.55121D+00, &
    479.66439D+00, 288.99238D+00, &
    493.09812D+00, 270.20636D+00, &
    518.87309D+00, 264.56427D+00, &
    547.18014D+00, 268.18846D+00, &
    600.49708D+00, 240.62570D+00, &
    625.96183D+00, 238.40347D+00, &
    633.90530D+00, 260.70629D+00, &
    621.50451D+00, 285.88914D+00, &
    576.87224D+00, 322.14121D+00, &
    570.51915D+00, 348.85423D+00, &
    567.16400D+00, 378.24075D+00, &
    558.00668D+00, 406.86552D+00, &
    565.19008D+00, 435.75599D+00, &
    567.56437D+00, 465.33407D+00, &
    550.87626D+00, 490.96358D+00, &
    532.98174D+00, 515.84491D+00, &
    500.66817D+00, 551.89078D+00, &
    478.75120D+00, 562.17222D+00, &
    430.03371D+00, 583.94286D+00, &
    401.20454D+00, 587.69910D+00, &
    368.32214D+00, 581.10110D+00, &
    354.26303D+00, 585.86085D+00, &
    346.75200D+00, 601.10367D+00, &
    332.85137D+00, 628.74602D+00, &
    308.02188D+00, 645.84180D+00, &
    295.52344D+00, 647.18525D+00, &
    286.51519D+00, 651.60328D+00, &
    285.98846D+00, 662.07339D+00, &
    298.93455D+00, 665.66316D+00, &
    301.70226D+00, 682.79570D+00, &
    278.65857D+00, 689.63850D+00, &
    266.25737D+00, 712.11005D+00, &
    287.28701D+00, 732.77147D+00, &
    318.19548D+00, 736.85151D+00, &
    343.83067D+00, 753.60957D+00, &
    375.53164D+00, 758.35231D+00, &
    405.73444D+00, 768.98687D+00, &
    406.33873D+00, 785.59001D+00, &
    378.35436D+00, 789.44240D+00, &
    350.02151D+00, 795.02238D+00, &
    338.68030D+00, 788.87325D+00, &
    325.67930D+00, 786.10177D+00, &
    319.05995D+00, 798.04657D+00, &
    301.78158D+00, 795.34254D+00, &
    280.69272D+00, 773.86634D+00, &
    254.55844D+00, 758.02898D+00, &
    234.07759D+00, 737.42090D+00, &
    218.38337D+00, 711.41500D+00, &
    220.99086D+00, 682.17833D+00, &
    224.50640D+00, 651.96297D+00, &
    240.25971D+00, 631.36117D+00, &
    259.86174D+00, 612.60253D+00, &
    291.85381D+00, 556.70385D+00, &
    315.52139D+00, 537.56387D+00, &
    341.63663D+00, 520.12519D+00, &
    351.37130D+00, 458.75372D+00, &
    349.33183D+00, 431.31454D+00, &
    328.80465D+00, 412.43055D+00 /), (/ 2, n1 /) )
  real ( kind = 8 ), dimension ( 2, n2 ) :: v2 = reshape ( (/ &
    238.64853D+00, 266.58978D+00, &
    235.14026D+00, 287.95183D+00, &
    238.20736D+00, 303.46785D+00, &
    250.13902D+00, 303.71290D+00, &
    258.51675D+00, 297.46973D+00, &
    274.55300D+00, 291.27357D+00, &
    284.66230D+00, 280.72063D+00, &
    279.73288D+00, 267.83455D+00, &
    270.68478D+00, 255.55440D+00, &
    255.73801D+00, 249.16872D+00, &
    241.72690D+00, 256.73448D+00 /), (/ 2, n2 /) )

  do j = 1, n

    call polygon_contains_point_2d ( n1, v1, point(1:m,j), inside1 )

    call polygon_contains_point_2d ( n2, v2, point(1:m,j), inside2 )

    inside(j) = ( inside1 .and. ( .not. inside2 ) )

  end do

  return
end
subroutine p14_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P14_SAMPLE samples points from the region in problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2, 2 ) :: box = reshape ( (/ &
    100.0D+00,  145.0D+00, &
   +634.0D+00, +799.0D+00 /), (/ 2, 2 /) )
  logical inside
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) reject
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(m)

  reject = 0

  do j = 1, n

    do

      call r8vec_uniform_01 ( m, seed, x )

      x(1:m) = ( 1.0D+00 - x(1:m) ) * box(1:m,1) &
                         + x(1:m)   * box(1:m,2)

      call p14_inside ( m, 1, x, inside )

      if ( inside ) then
        exit
      end if

      reject = reject + 1

      if ( 30 * n + 10 <= reject ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P14_SAMPLE - Fatal error!'
        write ( *, '(a,i8)' ) '  Trying to generate point J = ', j
        write ( *, '(a,i8)' ) '  Number of rejections = ', reject
        write ( *, '(a,g14.6)' ) '  Rejection percentage = ', &
          real ( 100 * reject ) / real ( reject + j - 1 )
          call r8vec_print ( m, x, '  Most recent rejected point: ' )
        stop
      end if

    end do

    point(1:m,j) = x(1:m)

  end do

  return
end
subroutine p14_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P14_SAMPLE_H1 samples points from the enlarged region in problem 14.
!
!  Discussion:
!
!    We enlarge the region by a layer H.  We do not round the
!    corners of the region, which should be done if we literally only
!    want to add points within H units of the boundary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: n1 = 90
  integer ( kind = 4 ), parameter :: n2 = 11

  real ( kind = 8 ), dimension ( 2, 2 ) :: box = reshape ( (/ &
    100.0D+00,  145.0D+00, &
   +634.0D+00, +799.0D+00 /), (/ 2, 2 /) )
  real ( kind = 8 ) h
  logical inside
  logical inside1
  logical inside2
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) reject
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension ( 2, n1 ) :: v1 = reshape ( (/ &
    316.43027D+00, 404.47559D+00, &
    291.04946D+00, 400.70917D+00, &
    265.16504D+00, 409.77890D+00, &
    241.46794D+00, 402.40310D+00, &
    216.55145D+00, 396.52064D+00, &
    163.28492D+00, 411.37102D+00, &
    142.81752D+00, 391.16355D+00, &
    111.95404D+00, 346.70264D+00, &
    100.03538D+00, 325.72710D+00, &
    103.98723D+00, 302.51587D+00, &
    128.72978D+00, 285.72802D+00, &
    147.49111D+00, 266.23345D+00, &
    196.65261D+00, 242.24055D+00, &
    213.56835D+00, 221.67192D+00, &
    226.49969D+00, 198.09326D+00, &
    248.37126D+00, 183.50473D+00, &
    262.21952D+00, 165.39102D+00, &
    278.42330D+00, 149.91715D+00, &
    300.71846D+00, 145.82601D+00, &
    311.12698D+00, 166.71094D+00, &
    326.66315D+00, 184.58335D+00, &
    359.78574D+00, 225.48049D+00, &
    357.08892D+00, 252.88958D+00, &
    358.76685D+00, 285.34403D+00, &
    361.50834D+00, 303.71287D+00, &
    371.68926D+00, 314.92452D+00, &
    380.49890D+00, 324.58632D+00, &
    396.37634D+00, 328.88990D+00, &
    412.59116D+00, 327.25238D+00, &
    425.48394D+00, 315.28623D+00, &
    435.84305D+00, 302.44664D+00, &
    458.34025D+00, 297.55121D+00, &
    479.66439D+00, 288.99238D+00, &
    493.09812D+00, 270.20636D+00, &
    518.87309D+00, 264.56427D+00, &
    547.18014D+00, 268.18846D+00, &
    600.49708D+00, 240.62570D+00, &
    625.96183D+00, 238.40347D+00, &
    633.90530D+00, 260.70629D+00, &
    621.50451D+00, 285.88914D+00, &
    576.87224D+00, 322.14121D+00, &
    570.51915D+00, 348.85423D+00, &
    567.16400D+00, 378.24075D+00, &
    558.00668D+00, 406.86552D+00, &
    565.19008D+00, 435.75599D+00, &
    567.56437D+00, 465.33407D+00, &
    550.87626D+00, 490.96358D+00, &
    532.98174D+00, 515.84491D+00, &
    500.66817D+00, 551.89078D+00, &
    478.75120D+00, 562.17222D+00, &
    430.03371D+00, 583.94286D+00, &
    401.20454D+00, 587.69910D+00, &
    368.32214D+00, 581.10110D+00, &
    354.26303D+00, 585.86085D+00, &
    346.75200D+00, 601.10367D+00, &
    332.85137D+00, 628.74602D+00, &
    308.02188D+00, 645.84180D+00, &
    295.52344D+00, 647.18525D+00, &
    286.51519D+00, 651.60328D+00, &
    285.98846D+00, 662.07339D+00, &
    298.93455D+00, 665.66316D+00, &
    301.70226D+00, 682.79570D+00, &
    278.65857D+00, 689.63850D+00, &
    266.25737D+00, 712.11005D+00, &
    287.28701D+00, 732.77147D+00, &
    318.19548D+00, 736.85151D+00, &
    343.83067D+00, 753.60957D+00, &
    375.53164D+00, 758.35231D+00, &
    405.73444D+00, 768.98687D+00, &
    406.33873D+00, 785.59001D+00, &
    378.35436D+00, 789.44240D+00, &
    350.02151D+00, 795.02238D+00, &
    338.68030D+00, 788.87325D+00, &
    325.67930D+00, 786.10177D+00, &
    319.05995D+00, 798.04657D+00, &
    301.78158D+00, 795.34254D+00, &
    280.69272D+00, 773.86634D+00, &
    254.55844D+00, 758.02898D+00, &
    234.07759D+00, 737.42090D+00, &
    218.38337D+00, 711.41500D+00, &
    220.99086D+00, 682.17833D+00, &
    224.50640D+00, 651.96297D+00, &
    240.25971D+00, 631.36117D+00, &
    259.86174D+00, 612.60253D+00, &
    291.85381D+00, 556.70385D+00, &
    315.52139D+00, 537.56387D+00, &
    341.63663D+00, 520.12519D+00, &
    351.37130D+00, 458.75372D+00, &
    349.33183D+00, 431.31454D+00, &
    328.80465D+00, 412.43055D+00 /), (/ 2, n1 /) )
  real ( kind = 8 ), dimension ( 2, n2 ) :: v2 = reshape ( (/ &
    238.64853D+00, 266.58978D+00, &
    235.14026D+00, 287.95183D+00, &
    238.20736D+00, 303.46785D+00, &
    250.13902D+00, 303.71290D+00, &
    258.51675D+00, 297.46973D+00, &
    274.55300D+00, 291.27357D+00, &
    284.66230D+00, 280.72063D+00, &
    279.73288D+00, 267.83455D+00, &
    270.68478D+00, 255.55440D+00, &
    255.73801D+00, 249.16872D+00, &
    241.72690D+00, 256.73448D+00 /), (/ 2, n2 /) )
  real ( kind = 8 ), dimension ( 2, n1 ) :: w1
  real ( kind = 8 ), dimension ( 2, n2 ) :: w2
  real ( kind = 8 ) x(m)
!
!  It would be wise to save W1 and W2 unless the value of H changes!
!
  call polygon_expand_2d ( n1, v1, h, w1 )
  call polygon_expand_2d ( n2, v2, h, w2 )

  reject = 0

  do j = 1, n

    do

      call r8vec_uniform_01 ( m, seed, x )

      x(1:m) = ( 1.0D+00 - x(1:m) ) * box(1:m,1) &
                         + x(1:m)   * box(1:m,2)

      call polygon_contains_point_2d ( n1, w1, x(1:m), inside1 )

      call polygon_contains_point_2d ( n2, w2, x(1:m), inside2 )

      inside = ( inside1 .and. ( .not. inside2 ) )

      if ( inside ) then
        exit
      end if

      reject = reject + 1

      if ( 30 * n + 10 <= reject ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P14_SAMPLE_H1 - Fatal error!'
        write ( *, '(a,i8)' ) '  Trying to generate point J = ', j
        write ( *, '(a,i8)' ) '  Number of rejections = ', reject
        write ( *, '(a,g14.6)' ) '  Rejection percentage = ', &
          real ( 100 * reject ) / real ( reject + j - 1 )
          call r8vec_print ( m, x, '  Most recent rejected point: ' )
        stop
      end if

    end do

    point(1:m,j) = x(1:m)

  end do

  return
end
subroutine p14_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P14_SDIST returns the signed distance to the region in problem 14.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P14_SDIST - Fatal error!'
  write ( *, '(a)' ) '  This routine is not written yet!'

  stop
end
subroutine p14_title ( title )

!*****************************************************************************80
!
!! P14_TITLE returns a title for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
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

  title = '#14: Marcus Garvie''s Lake Alpha with Beta Island.'

  return
end
subroutine p15_boundary_nearest ( m, n, point, boundary )

!*****************************************************************************80
!
!! P15_BOUNDARY_NEAREST returns a nearest boundary point in problem 15.
!
!  Discussion:
!
!    The nearest boundary point assignment is incorrect for regions
!    14 and 16, and for regions 21 and 23.  The boundary between each
!    of these pairs of regions is actually a parabola, but we have not
!    bothered to model these exactly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) BOUNDARY(M,N), points on the boundary
!    that are nearest to each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( m, n ) :: boundary
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( m, n ) :: point
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  do j = 1, n

    x = point(1,j)
    y = point(2,j)
!
!  INTERIOR REGIONS.
!
!  Region 1
!
    if ( -8.0D+00 <=     x              .and. &
         -9.0D+00 <= y + x              .and. &
                     y - x <= 8.0D+00 ) then

      boundary(1,j) = -8.0D+00
      boundary(2,j) = y
!
!  Region 2
!
    else if (             y + x <= -9.0D+00 .and. &
              -1.0D+00 <=     y             .and. &
                          y - x <= -3.0D+00 .and. &
                              y <=  0.5D+00 ) then

      boundary(1,j) = x
      boundary(2,j) = -1.0D+00
!
!  Region 3
!
    else if ( -3.0D+00 <= y - x             .and. &
                              x <=  2.0D+00 .and. &
                          y     <=  0.0D+00 .and. &
               0.0D+00 <=     x             ) then

      boundary(1,j) = 2.0D+00
      boundary(2,j) = y

!
!  Region 4
!
    else if (  0.0D+00 <= y                  .and. &
                               x <=  2.0D+00 .and. &
                           y     <=  0.5D+00 .and. &
                0.0D+00 <=     x             ) then

      boundary(1,j) = 2.0D+00
      boundary(2,j) = 0.0D+00
!
!  Region 5
!
    else if (  2.0D+00 <=      x             .and. &
               0.0D+00 <=  y                 .and. &
                           y - x <= -8.0D+00 .and. &
                           y     <=  0.5D+00 ) then
      boundary(1,j) = x
      boundary(2,j) = 0.0D+00
!
!  Region 6
!
    else if ( -8.0D+00 <= y - x             .and. &
                              x <=  8.0D+00 .and. &
                          y + x <=  9.0D+00 ) then
      boundary(1,j) = 8.0D+00
      boundary(2,j) = y
!
!  Region 7
!
    else if (  9.0D+00 <= y + x             .and. &
                          y     <=  1.0D+00 .and. &
               3.0D+00 <= y - x             .and. &
               0.5D+00 <= y                 ) then

      boundary(1,j) = x
      boundary(2,j) = 1.0D+00
!
!  Region 8
!
    else if (             y - x <=  3.0D+00 .and. &
              -2.0D+00 <=     x             .and. &
               0.0D+00 <= y                 .and. &
                              x <=  0.0D+00 ) then

      boundary(1,j) = -2.0D+00
      boundary(2,j) = y
!
!  Region 9
!
    else if (              y     <=  0.0D+00 .and. &
               -2.0D+00 <=     x             .and. &
               -0.5D+00 <= y                 .and. &
                               x <=  0.0D+00 ) then

      boundary(1,j) = -2.0D+00
      boundary(2,j) =  0.0D+00
!
!  Region 10
!
    else if (                  x <= -2.0D+00 .and. &
                           y     <=  0.0D+00 .and. &
                8.0D+00 <= y - x             .and. &
                0.5D+00 <= y                 ) then

      boundary(1,j) = x
      boundary(2,j) = 0.0D+00
!
!  EXTERIOR REGIONS.
!
!  Region 11
!
    else if (  -1.0D+00 <= y                 .and. &
                               x <= -8.0D+00 .and. &
                           y     <=  0.0D+00 ) then

      boundary(1,j) = -8.0D+00
      boundary(2,j) = y
!
!  Region 12
!
    else if (              y     <= -1.0D+00 .and. &
                               x <= -8.0D+00 ) then

      boundary(1,j) = -8.0D+00
      boundary(2,j) = -1.0D+00
!
!  Region 13
!
    else if (                  x <=  2.0D+00 .and. &
                           y     <= -1.0D+00 .and. &
               -8.0D+00 <=     x             ) then

      boundary(1,j) = x
      boundary(2,j) = -1.0D+00
!
!  Region 14
!  The boundary between regions 14 and 16 is only set approximately.
!
    else if (                  x <=  3.0D+00 .and. &
                           y     <= -1.0D+00 .and. &
                2.0D+00 <=     x             ) then

      boundary(1,j) =  2.0D+00
      boundary(2,j) = -1.0D+00
!
!  Region 15
!
    else if (   2.0D+00 <=     x             .and. &
               -1.0D+00 <= y                 .and. &
                           y - x <= -2.0D+00 ) then

      boundary(1,j) =  2.0D+00
      boundary(2,j) = y
!
!  Region 16
!  The boundary between regions 14 and 16 is only set approximately.
!
    else if ( (                 x <= 3.0D+00 .and. &
                            y     <= 0.0D+00 .and. &
                -2.0D+00 <= y - x            )     &
            .or.                                   &
              (                 x <= 8.0D+00 .and. &
                            y     <= 0.0D+00 .and. &
                 3.0D+00 <=     x            ) ) then

      boundary(1,j) = x
      boundary(2,j) = 0.0D+00
!
!  Region 17
!
    else if (   8.0D+00 <=     x             .and. &
                           y     <=  0.0D+00 ) then

      boundary(1,j) = 8.0D+00
      boundary(2,j) = 0.0D+00
!
!  Region 18
!
    else if (              y     <= 1.0D+00  .and. &
                8.0D+00 <=     x             .and. &
                0.0D+00 <= y                 ) then

      boundary(1,j) = 8.0D+00
      boundary(2,j) = y
!
!  Region 19
!
    else if (   8.0D+00 <=     x             .and. &
                1.0D+00 <= y                 ) then

      boundary(1,j) = 8.0D+00
      boundary(2,j) = 1.0D+00
!
!  Region 20
!
    else if (  -2.0D+00 <=     x             .and. &
                1.0D+00 <= y                 .and. &
                               x <= 8.0D+00 ) then

      boundary(1,j) = x
      boundary(2,j) = 1.0D+00
!
!  Region 21
!  The boundary between regions 21 and 23 is only set approximately.
!
    else if (                  x <= -2.0D+00 .and. &
                1.0D+00 <= y                 .and. &
               -3.0D+00 <=     x             ) then

      boundary(1,j) = -2.0D+00
      boundary(2,j) = 1.0D+00
!
!  Region 22
!
    else if (                  x <= -2.0D+00 .and. &
                           y     <=  1.0D+00 .and. &
                2.0D+00 <= y - x             ) then

      boundary(1,j) = -2.0D+00
      boundary(2,j) = y
!
!  Region 23
!  The boundary between regions 21 and 23 is only set approximately.
!
    else if ( (             y - x <=  2.0D+00 .and. &
                 0.0D+00 <= y                 .and. &
                -3.0D+00 <=     x                   &
              ) .or. (                              &
                                x <= -3.0D+00 .and. &
                 0.0D+00 <= y                 .and. &
                -8.0D+00 <=     x ) ) then

      boundary(1,j) = x
      boundary(2,j) = 1.0D+00
!
!  Region 24
!
    else if ( x <= -8.0D+00 .and. &
              0.0D+00 <= y ) then

      boundary(1,j) = -8.0D+00
      boundary(2,j) = 0.0D+00
!
!  Somehow, a point missed all the regions.
!  This should not happen.
!
    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P15_BOUNDARY_NEAREST - Fatal error!'
      write ( *, '(a,2g14.6)' ) '  Orphan point = ', point(1:2,j)
      stop

    end if

  end do

  return
end
subroutine p15_boundary_project ( m, n, point )

!*****************************************************************************80
!
!! P15_BOUNDARY_PROJECT projects exterior points to the boundary in problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.  On output, exterior points have been projected.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical inside(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) pn(m)
  real ( kind = 8 ), dimension ( m, n ) :: point

  call p15_inside ( m, n, point, inside )

  do j = 1, n

    if ( .not. inside(j) ) then
      call p15_boundary_nearest ( m, 1, point(1:m,j), pn )
      point(1:m,j) = pn(1:m)
    end if

  end do

  return
end
subroutine p15_boundary_segment ( segment_index, m, segment_length, &
  segment )

!*****************************************************************************80
!
!! P15_BOUNDARY_SEGMENT returns a boundary segment in problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of the boundary segment.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
!    Output, real ( kind = 8 ) SEGMENT(M,SEGMENT_LENGTH), the
!    points that make up the boundary segment.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) segment_length

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  integer ( kind = 4 ) n7
  integer ( kind = 4 ) n8
  real ( kind = 8 ) s(m,8)
  integer ( kind = 4 ) segment_index
  real ( kind = 8 ) segment(m,segment_length)
  integer ( kind = 4 ) sl_extra

  if ( segment_index == 1 ) then

    sl_extra = segment_length - 9

    n1 = nint ( 10.0D+00 * real ( sl_extra, kind = 8 ) &
              / real ( 36, kind = 8 ) )

    n2 = nint (  1.0D+00 * real ( sl_extra, kind = 8 ) &
              / real ( 36, kind = 8 ) )
    n2 = max ( n2, 0 )

    n3 = nint (  6.0D+00 * real ( sl_extra, kind = 8 ) &
                / real ( 36, kind = 8 ) )
    n3 = max ( n3, 0 )

    n4 = nint (  1.0D+00 * real ( sl_extra, kind = 8 ) &
                / real ( 36, kind = 8 ) )
    n4 = max ( n4, 0 )

    n5 = nint ( 10.0D+00 * real ( sl_extra, kind = 8 ) &
                / real ( 36, kind = 8 ) )
    n5 = max ( n5, 0 )

    n6 = nint (  1.0D+00 * real ( sl_extra, kind = 8 ) &
                / real ( 36, kind = 8 ) )
    n6 = max ( n6, 0 )

    n7 = nint (  6.0D+00 * real ( sl_extra, kind = 8 ) &
                / real ( 36, kind = 8 ) )
    n7 = max ( n7, 0 )

    n8 = sl_extra - n1 - n2 - n3 - n4 - n5 - n6 - n7

    n1 = n1 + 1
    n2 = n2 + 1
    n3 = n3 + 1
    n4 = n4 + 1
    n5 = n5 + 1
    n6 = n6 + 1
    n7 = n7 + 1
    n8 = n8 + 1

    s(1:2,1) = (/  -8.0D+00,  -1.0D+00 /)
    s(1:2,2) = (/   2.0D+00,  -1.0D+00 /)
    s(1:2,3) = (/   2.0D+00,   0.0D+00 /)
    s(1:2,4) = (/   8.0D+00,   0.0D+00 /)
    s(1:2,5) = (/   8.0D+00,   1.0D+00 /)
    s(1:2,6) = (/  -2.0D+00,   1.0D+00 /)
    s(1:2,7) = (/  -2.0D+00,   0.0D+00 /)
    s(1:2,8) = (/  -8.0D+00,   0.0D+00 /)

    j = 0

    do i = 1, n1
      j = j + 1
      segment(1:2,j) = ( real ( n1 - i + 1, kind = 8 ) * s(1:2,1)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,2) ) &
                       / real ( n1,         kind = 8 )
    end do

    do i = 1, n2
      j = j + 1
      segment(1:2,j) = ( real ( n2 - i + 1, kind = 8 ) * s(1:2,2)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,3) ) &
                       / real ( n2,         kind = 8 )
    end do

    do i = 1, n3
      j = j + 1
      segment(1:2,j) = ( real ( n3 - i + 1, kind = 8 ) * s(1:2,3)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,4) ) &
                       / real ( n3,         kind = 8 )
    end do

    do i = 1, n4
      j = j + 1
      segment(1:2,j) = ( real ( n4 - i + 1, kind = 8 ) * s(1:2,4)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,5) ) &
                       / real ( n4,         kind = 8 )
    end do

    do i = 1, n5
      j = j + 1
      segment(1:2,j) = ( real ( n5 - i + 1, kind = 8 ) * s(1:2,5)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,6) ) &
                       / real ( n5,         kind = 8 )
    end do

    do i = 1, n6
      j = j + 1
      segment(1:2,j) = ( real ( n6 - i + 1, kind = 8 ) * s(1:2,6)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,7) ) &
                       / real ( n6,         kind = 8 )
    end do

    do i = 1, n7
      j = j + 1
      segment(1:2,j) = ( real ( n7 - i + 1, kind = 8 ) * s(1:2,7)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,8) ) &
                       / real ( n7,         kind = 8 )
    end do

    do i = 1, n8
      j = j + 1
      segment(1:2,j) = ( real ( n8 - i + 1, kind = 8 ) * s(1:2,8)   &
                       + real (      i - 1, kind = 8 ) * s(1:2,1) ) &
                       / real ( n8,         kind = 8 )
    end do

    j = j + 1
    segment(1:2,j) = s(1:2,1)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P15_BOUNDARY_SEGMENT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p15_boundary_segment_length ( segment_index, h, segment_length )

!*****************************************************************************80
!
!! P15_BOUNDARY_SEGMENT_LENGTH returns boundary segment lengths in problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEGMENT_INDEX, the index of one of the boundary segments.
!
!    Input, real ( kind = 8 ) H, the suggested spacing between points.
!
!    Output, integer ( kind = 4 ) SEGMENT_LENGTH, the number of points in the segment.
!
  implicit none

  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  integer ( kind = 4 ) segment_index
  integer ( kind = 4 ) segment_length

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P15_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Nonpositive H = ', h
    stop
  end if

  if ( segment_index == 1 ) then

    n = nint ( 36.0D+00 / h )
    n = max ( n, 17 )
    segment_length = n + mod ( 36 - mod ( n - 1, 36 ), 36 )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P15_BOUNDARY_SEGMENT_LENGTH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SEGMENT_INDEX = ', segment_index
    stop

  end if

  return
end
subroutine p15_boundary_segment_num ( boundary_segment_num )

!*****************************************************************************80
!
!! P15_BOUNDARY_SEGMENT_NUM counts the boundary segments in problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BOUNDARY_SEGMENT_NUM, the number of boundary segments.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num

  boundary_segment_num = 1

  return
end
subroutine p15_box ( m, lo, hi )

!*****************************************************************************80
!
!! P15_BOX returns a bounding box for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) LO(M), HI(M), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) lo(m)

  lo(1:m) = (/ -8.0D+00, -1.0D+00 /)
  hi(1:m) = (/ +8.0D+00, +1.0D+00 /)

  return
end
subroutine p15_density ( m, n, point, density )

!*****************************************************************************80
!
!! P15_DENSITY returns the density for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) DENSITY(N), the mesh density at
!    each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) density(n)
  real ( kind = 8 ) point(m,n)

  density(1:n) = 1.0D+00

  return
end
subroutine p15_element_size ( element_size )

!*****************************************************************************80
!
!! P15_ELEMENT_SIZE returns a typical element size for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_SIZE, a typical element size.
!
  implicit none

  real ( kind = 8 ) element_size

  element_size = 0.2D+00

  return
end
subroutine p15_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P15_FIXED_NUM returns the number of fixed points in problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
  implicit none

  integer ( kind = 4 ) fixed_num

  fixed_num = 8

  return
end
subroutine p15_fixed_points ( m, fixed_num, fixed )

!*****************************************************************************80
!
!! P15_FIXED_POINTS returns the fixed points in problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = 8 ) FIXED(M,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) fixed_num

  real ( kind = 8 ) fixed(m,fixed_num)

  fixed(1:m,1:fixed_num) = reshape ( (/ &
   -8.0D+00, -1.0D+00, &
    2.0D+00, -1.0D+00, &
    2.0D+00,  0.0D+00, &
    8.0D+00,  0.0D+00, &
    8.0D+00,  1.0D+00, &
   -2.0D+00,  1.0D+00, &
   -2.0D+00,  0.0D+00, &
   -8.0D+00,  0.0D+00 /), (/ m, fixed_num /) )

  return
end
subroutine p15_header ( )

!*****************************************************************************80
!
!! P15_HEADER prints some information about problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None.
!
  implicit none

  integer ( kind = 4 ) boundary_segment_num
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) hole_num

  call p15_boundary_segment_num ( boundary_segment_num )
  call p15_fixed_num ( fixed_num )
  call p15_hole_num ( hole_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P15:'
  write ( *, '(a)' ) '  Sangbum Kim''s forward step.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of boundary segments = ', boundary_segment_num
  write ( *, '(a,i8)' ) '  Number of fixed points =      ', fixed_num
  write ( *, '(a,i8)' ) '  Number of holes =             ', hole_num

  return
end
subroutine p15_hole_num ( hole_num )

!*****************************************************************************80
!
!! P15_HOLE_NUM counts the holes in problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 0

  return
end
subroutine p15_hole_point ( hole_index, m, hole_point )

!*****************************************************************************80
!
!! P15_HOLE_POINT returns a point inside a given hole in problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_INDEX, the index of the hole.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) HOLE_POINT(M), a point in the hole
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) hole_index
  real ( kind = 8 ) hole_point(m)

  return
end
subroutine p15_inside ( m, n, point, inside )

!*****************************************************************************80
!
!! P15_INSIDE reports if a point is inside the region in problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical inside(n)
  real ( kind = 8 ) point(m,n)

  inside(1:n) =                                                       &
       ( -8.0D+00 <= point(1,1:n) .and. point(1,1:n) <= 2.0D+00 .and. &
         -1.0D+00 <= point(2,1:n) .and. point(2,1:n) <= 0.0D+00    )  &
    .or.                                                              &
       ( -2.0D+00 <= point(1,1:n) .and. point(1,1:n) <= 8.0D+00 .and. &
          0.0D+00 <= point(2,1:n) .and. point(2,1:n) <= 1.0D+00 )

  return
end
subroutine p15_sample ( m, n, seed, point )

!*****************************************************************************80
!
!! P15_SAMPLE samples points from the region in problem 15.
!
!  Discussion:
!
!    With a little bit of work, we can guarantee that we don't have to
!    use a rejection method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) area
  real ( kind = 8 ) prob(n)
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4

  x1 = -8.0D+00
  x2 =  2.0D+00
  x3 = -2.0D+00
  x4 =  8.0D+00

  y1 = -1.0D+00
  y2 =  0.0D+00
  y3 =  0.0D+00
  y4 =  1.0D+00
!
!  Generate a batch of points in [0,1]x[0,1].
!
  call r8mat_uniform_01 ( m, n, seed, point )
!
!  Generate a batch of N probabilities.
!
  call r8vec_uniform_01 ( n, seed, prob )
!
!  Map half the points into [X1,X2] x [Y1,Y2].
!
  where ( prob(1:n) < 0.5D+00 )

    point(1,1:n) = x1 + point(1,1:n) * ( x2 - x1 )
    point(2,1:n) = y1 + point(2,1:n) * ( y2 - y1 )
!
!  Map the other points into [X3,X4] x [Y3,Y4].
!
  elsewhere

    point(1,1:n) = x3 + point(1,1:n) * ( x4 - x3 )
    point(2,1:n) = y3 + point(2,1:n) * ( y4 - y3 )

  end where

  return
end
subroutine p15_sample_h1 ( m, n, h, seed, point )

!*****************************************************************************80
!
!! P15_SAMPLE_H1 samples points from the enlarged region in problem 15.
!
!  Discussion:
!
!    With a little bit of work, we can guarantee that we don't have to
!    use a rejection method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) H, the enlargement.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) area
  real ( kind = 8 ) h
  real ( kind = 8 ) prob(n)
  real ( kind = 8 ) point(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4

  x1 = 0.0D+00 - h
  x2 = 0.5D+00 + h
  x3 = 1.0D+00 + h

  y1 = 0.0D+00 - h
  y2 = 0.5D+00 + h
  y3 = 1.0D+00 + h

  x1 = -8.0D+00
  x2 =  2.0D+00
  x3 = -2.0D+00
  x4 =  8.0D+00

  y1 = -1.0D+00
  y2 =  0.0D+00
  y3 =  0.0D+00
  y4 =  1.0D+00
!
!  Generate a batch of points in [0,1]x[0,1].
!
  call r8mat_uniform_01 ( m, n, seed, point )
!
!  Generate a batch of N probabilities.
!
  call r8vec_uniform_01 ( n, seed, prob )
!
!  Map half the points into [X1-H,X2+H] x [Y1-H,Y2+H].
!
  where ( prob(1:n) < 0.5D+00 )

    point(1,1:n) = x1 - h + point(1,1:n) * ( x2 + h - ( x1 - h ) )
    point(2,1:n) = y1 - h + point(2,1:n) * ( y2 + h - ( y1 - h ) )
!
!  Map the other points into [X3-H,X4+H] x [Y3-H,Y4+H].
!
  elsewhere

    point(1,1:n) = x3 - h + point(1,1:n) * ( x4 + h - ( x3 - h ) )
    point(2,1:n) = y3 - h + point(2,1:n) * ( y4 + h - ( y3 - h ) )

  end where
!
!  To be strictly correct, you need to discard HALF the points in
!  [X3-H x X2 + H ] * [ Y3 - H, Y2 + H ].
!
  return
end
subroutine p15_sdist ( m, n, point, sdist )

!*****************************************************************************80
!
!! P15_SDIST returns the signed distance to the region in problem 15.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!    For this calculation, we rely on the code from P15_BOUNDARY_NEAREST.
!    That code makes some approximations for a couple regions, but
!    for points closer than 1 unit or so, the results should be correct,
!    and we don't really need to be so fussy for what we are doing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = 8 ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sdist(n)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  do j = 1, n

    x = point(1,j)
    y = point(2,j)
!
!  INTERIOR REGIONS.
!
!  Region 1.
!
    if ( -8.0D+00 <=     x              .and. &
         -9.0D+00 <= y + x              .and. &
                     y - x <= 8.0D+00 ) then

      sdist(j) = - x - 8.0D+00
!
!  Region 2.
!
    else if (             y + x <= -9.0D+00 .and. &
              -1.0D+00 <=     y             .and. &
                          y - x <= -3.0D+00 .and. &
                              y <=  0.5D+00 ) then

      sdist(j) =  - y - 1.0D+00

!
!  Region 3.
!
    else if ( -3.0D+00 <= y - x             .and. &
                              x <=  2.0D+00 .and. &
                          y     <=  0.0D+00 .and. &
               0.0D+00 <=     x             ) then

      sdist(j) = x - 2.0D+00

!
!  Region 4.
!
    else if (  0.0D+00 <= y                  .and. &
                               x <=  2.0D+00 .and. &
                           y     <=  0.5D+00 .and. &
                0.0D+00 <=     x             ) then

      sdist(j) = - sqrt ( ( x - 2.0D+00 )**2 + y**2 )
!
!  Region 5.
!
    else if (  2.0D+00 <=      x             .and. &
               0.0D+00 <=  y                 .and. &
                           y - x <= -8.0D+00 .and. &
                           y     <=  0.5D+00 ) then
      sdist(j) = - y
!
!  Region 6.
!
    else if ( -8.0D+00 <= y - x             .and. &
                              x <=  8.0D+00 .and. &
                          y + x <=  9.0D+00 ) then
      sdist(j) = x - 8.0D+00
!
!  Region 7.
!
    else if (  9.0D+00 <= y + x             .and. &
                          y     <=  1.0D+00 .and. &
               3.0D+00 <= y - x             .and. &
               0.5D+00 <= y                 ) then

      sdist(j) = y - 1.0D+00
!
!  Region 8.
!
    else if (             y - x <=  3.0D+00 .and. &
              -2.0D+00 <=     x             .and. &
               0.0D+00 <= y                 .and. &
                              x <=  0.0D+00 ) then

      sdist(j) = - x - 2.0D+00
!
!  Region 9.
!
    else if (              y     <=  0.0D+00 .and. &
               -2.0D+00 <=     x             .and. &
               -0.5D+00 <= y                 .and. &
                               x <=  0.0D+00 ) then

      sdist(j) = - sqrt ( ( x -2.0D+00 )**2 + y**2 )
!
!  Region 10
!
    else if (                  x <= -2.0D+00 .and. &
                           y     <=  0.0D+00 .and. &
                8.0D+00 <= y - x             .and. &
                0.5D+00 <= y                 ) then

      sdist(j) = y
!
!  EXTERIOR REGIONS.
!
!  Region 11.
!
    else if (  -1.0D+00 <= y                 .and. &
                               x <= -8.0D+00 .and. &
                           y     <=  0.0D+00 ) then

      sdist(j) = - x - 8.0D+00
!
!  Region 12.
!
    else if (              y     <= -1.0D+00 .and. &
                               x <= -8.0D+00 ) then

      sdist(j) = sqrt ( ( x + 8.0D+00 )**2 + ( y + 1.0D+00 )**2 )
!
!  Region 13
!
    else if (                  x <=  2.0D+00 .and. &
                           y     <= -1.0D+00 .and. &
               -8.0D+00 <=     x             ) then

      sdist(j) = - y - 1.0D+00
!
!  Region 14.
!  The boundary between regions 14 and 16 is only set approximately.
!
    else if (                  x <=  3.0D+00 .and. &
                           y     <= -1.0D+00 .and. &
                2.0D+00 <=     x             ) then

      sdist(j) = sqrt ( ( x - 2.0D+00 )**2 + ( y + 1.0D+00 )**2 )
!
!  Region 15
!
    else if (   2.0D+00 <=     x             .and. &
               -1.0D+00 <= y                 .and. &
                           y - x <= -2.0D+00 ) then

      sdist(j) = x - 2.0D+00
!
!  Region 16
!  The boundary between regions 14 and 16 is only set approximately.
!
    else if ( (                 x <= 3.0D+00 .and. &
                            y     <= 0.0D+00 .and. &
                -2.0D+00 <= y - x            )     &
            .or.                                   &
              (                 x <= 8.0D+00 .and. &
                            y     <= 0.0D+00 .and. &
                 3.0D+00 <=     x            ) ) then

      sdist(j) = - y
!
!  Region 17
!
    else if (   8.0D+00 <=     x             .and. &
                           y     <=  0.0D+00 ) then

      sdist(j) = sqrt ( ( x - 8.0D+00 )**2 + y**2 )
!
!  Region 18.
!
    else if (              y     <= 1.0D+00  .and. &
                8.0D+00 <=     x             .and. &
                0.0D+00 <= y                 ) then

      sdist(j) = x - 8.0D+00
!
!  Region 19.
!
    else if (   8.0D+00 <=     x             .and. &
                1.0D+00 <= y                 ) then

      sdist(j) = sqrt ( ( x - 8.0D+00 )**2 + ( y - 1.0D+00 )**2 )
!
!  Region 20
!
    else if (  -2.0D+00 <=     x             .and. &
                1.0D+00 <= y                 .and. &
                               x <= 8.0D+00 ) then

      sdist(j) = y - 1.0D+00
!
!  Region 21.
!  The boundary between regions 21 and 23 is only set approximately.
!
    else if (                  x <= -2.0D+00 .and. &
                1.0D+00 <= y                 .and. &
               -3.0D+00 <=     x             ) then

      sdist(j) = sqrt ( ( x + 2.0D+00 )**2 + ( y - 1.0D+00 )**2 )
!
!  Region 22
!
    else if (                  x <= -2.0D+00 .and. &
                           y     <=  1.0D+00 .and. &
                2.0D+00 <= y - x             ) then

      sdist(j) = - x - 2.0D+00
!
!  Region 23
!  The boundary between regions 21 and 23 is only set approximately.
!
    else if ( (             y - x <=  2.0D+00 .and. &
                 0.0D+00 <= y                 .and. &
                -3.0D+00 <=     x                   &
              ) .or. (                              &
                                x <= -3.0D+00 .and. &
                 0.0D+00 <= y                 .and. &
                -8.0D+00 <=     x ) ) then

      sdist(j) = y
!
!  Region 24.
!
    else if ( x <= -8.0D+00 .and. &
              0.0D+00 <= y ) then

      sdist(j) = sqrt ( ( x + 8.0D+00 )**2 + y**2 )
!
!  Somehow, a point missed all the regions.
!  This should not happen.
!
    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P15_SDIST - Fatal error!'
      write ( *, '(a,2g14.6)' ) '  Orphan point = ', point(1:2,j)
      stop

    end if

  end do

  stop
end
subroutine p15_title ( title )

!*****************************************************************************80
!
!! P15_TITLE returns a title for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2006
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

  title = '#15: Sangbum Kim''s forward step.'

  return
end
subroutine poly_write ( file_name, node_num, segment, edge_num, &
  edge_nodes, hole_num, hole_point )

!*****************************************************************************80
!
!! POLY_WRITE writes data to a POLY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to create.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) SEGMENT(2,NODE_NUM), the nodes.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Input, integer ( kind = 4 ) EDGE_NODES(2,EDGE_NUM), the nodes that form each edge.
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of holes in the mesh.
!
!    Input, real ( kind = 8 ) HOLE_POINT(2,HOLE_NUM), a point in each hole.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num

  character ( len = 8 ) date
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_nodes(2,edge_num)
  character ( len = * ) file_name
  integer ( kind = 4 ) hole
  real ( kind = 8 ) hole_point(2,hole_num)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) node
  integer ( kind = 4 ) poly_unit
  integer ( kind = 4 ), parameter :: region_num = 0
  real ( kind = 8 ) segment(2,node_num)
  character ( len = 40 ) string

  call get_unit ( poly_unit )

  open ( unit = poly_unit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLY_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output POLY file.'
    stop
  end if

  write ( poly_unit, '(a)' ) '#  ' // trim ( file_name )
  write ( poly_unit, '(a)' ) &
    '#  Created by poly_write(test_triangulation.f90)'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Vertex  Dimension  Attribute  Marker'
  write ( poly_unit, '(a)' ) '#  Count              Count      0/1'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(2x,i8,a)' ) node_num, '  2  0  0'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Vertex  X  Y  Attributes  Marker'
  write ( poly_unit, '(a)' ) '#  Index'
  write ( poly_unit, '(a)' ) '#'
  do node = 1, node_num
    write ( poly_unit, '(2x,i8,2x,f10.6,2x,f10.6)' ) &
      node, segment(1,node), segment(2,node)
  end do
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Segment  Marker'
  write ( poly_unit, '(a)' ) '#  Count    0/1'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(2x,i8,a)' ) edge_num, '  0'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Segment_index  Node1  Node2  Marker'
  write ( poly_unit, '(a)' ) '#'
  do edge = 1, edge_num
    write ( poly_unit, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) &
      edge, edge_nodes(1,edge), edge_nodes(2,edge), 0
  end do
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Hole'
  write ( poly_unit, '(a)' ) '#  Count'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(2x,i8)' ) hole_num
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Hole_index  X  Y'
  write ( poly_unit, '(a)' ) '#'
  do hole = 1, hole_num
    write ( poly_unit, '(2x,i8,2x,f10.6,2x,f10.6)' ) &
      hole, hole_point(1,hole), hole_point(2,hole)
  end do
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Region'
  write ( poly_unit, '(a)' ) '#  Count'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(2x,i8)' ) region_num

  close ( unit = poly_unit )

  return
end
subroutine polygon_contains_point_2d ( n, v, p, inside )

!*****************************************************************************80
!
!! POLYGON_CONTAINS_POINT_2D finds if a point is inside a simple polygon in 2D.
!
!  Discussion:
!
!    A simple polygon is one whose boundary never crosses itself.
!    The polygon does not need to be convex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    M Shimrat,
!    Position of Point Relative to Polygon,
!    ACM Algorithm 112,
!    Communications of the ACM,
!    Volume 5, Number 8, page 434, August 1962.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes or vertices in the polygon.
!    N must be at least 3.
!
!    Input, real ( kind = 8 ) V(2,N), the vertices of the polygon.
!
!    Input, real ( kind = 8 ) P(2), the coordinates of the point to be tested.
!
!    Output, logical INSIDE, is TRUE if the point is inside the polygon.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) i
  logical inside
  integer ( kind = 4 ) ip1
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) v(dim_num,n)

  inside = .true.

  do i = 1, n

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    if ( ( v(2,i)   <  p(2) .and. p(2) <= v(2,ip1)   ) .or. &
         ( p(2) <= v(2,i)   .and. v(2,ip1)   < p(2) ) ) then
      if ( ( p(1) - v(1,i) ) - ( p(2) - v(2,i) ) &
         * ( v(1,ip1) - v(1,i) ) / ( v(2,ip1) - v(2,i) ) < 0 ) then
        inside = .not. inside
      end if
    end if

  end do

  inside = .not. inside

  return
end
subroutine polygon_expand_2d ( n, v, h, w )

!*****************************************************************************80
!
!! POLYGON_EXPAND_2D expands a polygon in 2D.
!
!  Discussion:
!
!    This routine simple moves each vertex of the polygon outwards
!    in such a way that the sides of the polygon advance by H.
!
!    This approach should always work if the polygon is convex, or
!    star-shaped.  But for general polygons, it is possible
!    that this procedure, for large enough H, will create a polygon
!    whose sides intersect.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of sides of the polygon.
!
!    Input, real ( kind = 8 ) V(2,N), the coordinates of the vertices.
!
!    Input, real ( kind = 8 ) H, the expansion amount.
!
!    Output, real ( kind = 8 ) W(2,N), the "expanded" coordinates.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle_rad_2d
  real ( kind = 8 ) h
  real ( kind = 8 ) h2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) ip1
  real ( kind = 8 ) p4(dim_num)
  real ( kind = 8 ) v(dim_num,n)
  real ( kind = 8 ) w(dim_num,n)
!
!  Consider each angle, formed by the nodes P(I-1), P(I), P(I+1).
!
  do i = 1, n

    im1 = i4_wrap ( i-1, 1, n )
    ip1 = i4_wrap ( i+1, 1, n )
!
!        P1
!        /
!       /   P4
!      /  .
!     / .
!    P2--------->P3
!
    call angle_half_2d ( v(1:dim_num,im1), v(1:dim_num,i), v(1:dim_num,ip1), p4 )
!
!  Compute the value of the half angle.
!
    angle = angle_rad_2d ( v(1:dim_num,im1), v(1:dim_num,i), p4(1:dim_num) )
!
!  The stepsize along the ray must be adjusted so that the sides
!  move out by H.
!
    h2 = h / sin ( angle )

    w(1:dim_num,i) = v(1:dim_num,i) - h2 * ( p4(1:dim_num) - v(1:dim_num,i) )

  end do

  return
end
subroutine polygon_point_near_2d ( n, v, p, pn, dist )

!*****************************************************************************80
!
!! POLYGON_POINT_NEAR_2D computes the nearest point on a polygon in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V(2,N), the polygon vertices.
!
!    Input, real ( kind = 8 ) P(2), the point whose nearest polygon point
!    is to be determined.
!
!    Output, real ( kind = 8 ) PN(2), the nearest point to P.
!
!    Output, real ( kind = 8 ) DIST, the distance from the point to the
!    polygon.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) dist
  real ( kind = 8 ) dist2
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) pn(dim_num)
  real ( kind = 8 ) pn2(dim_num)
  real ( kind = 8 ) tval
  real ( kind = 8 ) v(dim_num,n)
!
!  Find the distance to each of the line segments that make up the edges
!  of the polygon.
!
  dist = huge ( dist )
  pn(1:dim_num) = 0.0D+00

  do j = 1, n

    jp1 = i4_wrap ( j+1, 1, n )

    call segment_point_near_2d ( v(1:dim_num,j), v(1:dim_num,jp1), p, &
      pn2, dist2, tval )

    if ( dist2 < dist ) then
      dist = dist2
      pn(1:dim_num) = pn2(1:dim_num)
    end if

  end do

  return
end
subroutine polyloop_arclength_nd ( dim_num, nk, pk, sk )

!*****************************************************************************80
!
!! POLYLOOP_ARCLENGTH_ND computes the arclength of points on a polyloop in ND.
!
!  Discussion:
!
!    A polyloop of order NK is the geometric structure consisting of
!    the NK line segments that lie between successive elements of a list
!    of NK points, with the last point joined to the first.
!
!    Warning: I just made up the word "polyloop", so don't go repeating it!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NK, the number of points defining the polyloop.
!
!    Input, real ( kind = 8 ) PK(DIM_NUM,NK), the points defining the polyloop.
!
!    Output, real ( kind = 8 ) SK(NK+1), the arclength coordinates
!    of each point.  The first point has two arc length values,
!    namely SK(1) = 0 and SK(NK+1) = LENGTH.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nk

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) pk(dim_num,nk)
  real ( kind = 8 ) sk(nk+1)

  sk(1) = 0.0D+00

  do i = 2, nk + 1

    if ( i <= nk ) then
      j = i
    else
      j = 1
    end if

    sk(i) = sk(i-1) + sqrt ( sum ( ( pk(1:dim_num,j) - pk(1:dim_num,i-1) )**2 ) )

  end do

  return
end
subroutine polyloop_length_nd ( dim_num, nk, pk, length )

!*****************************************************************************80
!
!! POLYLOOP_LENGTH_ND computes the length of a polyloop in ND.
!
!  Discussion:
!
!    A polyloop of order NK is the geometric structure consisting of
!    the NK line segments that lie between successive elements of a list
!    of NK points, with the last point joined to the first.
!
!    Warning: I just made up the word "polyloop", so don't go repeating it!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NK, the number of points defining the polyloop.
!
!    Input, real ( kind = 8 ) PK(DIM_NUM,NK), the points defining the polyloop.
!
!    Output, real ( kind = 8 ) LENGTH, the length of the polyloop.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nk

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) length
  real ( kind = 8 ) pk(dim_num,nk)

  length = 0.0D+00

  do i = 2, nk + 1

    if ( i <= nk ) then
      j = i
    else
      j = 1
    end if

    length = length + sqrt ( sum ( ( pk(1:dim_num,j) - pk(1:dim_num,i-1) )**2 ) )

  end do

  return
end
subroutine polyloop_points_nd ( dim_num, nk, pk, nt, pt )

!*****************************************************************************80
!
!! POLYLOOP_POINTS_ND computes equally spaced points on a polyloop in ND.
!
!  Discussion:
!
!    A polyloop of order NK is the geometric structure consisting of
!    the NK line segments that lie between successive elements of a list
!    of NK points, including a segment from the last point to the first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NK, the number of points defining the polyloop.
!
!    Input, real ( kind = 8 ) PK(DIM_NUM,NK), the points defining the polyloop.
!
!    Input, integer ( kind = 4 ) NT, the number of points to be sampled.
!
!    Input, real ( kind = 8 ) PT(DIM_NUM,NT), equally spaced points
!    on the polyloop.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nk
  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  real ( kind = 8 ) pk(dim_num,nk)
  real ( kind = 8 ) pt(dim_num,nt)
  real ( kind = 8 ) sk(nk+1)
  real ( kind = 8 ) st

  call polyloop_arclength_nd ( dim_num, nk, pk, sk )

  j = 1

  do i = 1,  nt

    st = ( real ( nt - i,     kind = 8 ) * 0.0D+00 + &
           real (      i - 1, kind = 8 ) * sk(nk+1) ) &
         / real ( nt     - 1, kind = 8 )

    do

      if ( sk(j) <= st .and. st <= sk(j+1) ) then
        exit
      end if

      if ( nk <= j ) then
        exit
      end if

      j = j + 1

    end do

    jp1 = i4_wrap ( j + 1, 1, nk )

    pt(1:dim_num,i) = ( ( sk(j+1) - st         ) * pk(1:dim_num,j) &
                      + (           st - sk(j) ) * pk(1:dim_num,jp1) ) &
                      / ( sk(j+1)      - sk(j) )

  end do

  return
end
function quad_contains_point_2d ( q, p )

!*****************************************************************************80
!
!! QUAD_CONTAINS_POINT_2D finds if a point is inside a convex quadrilateral in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Q(2,4), the vertices of the quadrilateral.
!
!    Input, real ( kind = 8 ) P(2), the point to be checked.
!
!    Output, logical QUAD_CONTAINS_POINT, is TRUE if the point is
!    in the quadrilateral.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) q(dim_num,4)
  logical quad_contains_point_2d
  real ( kind = 8 ) t(dim_num,3)
  logical triangle_contains_point_2d

  t(1:2,1:3) = q(1:2,1:3)

  if ( triangle_contains_point_2d ( t, p ) ) then
    quad_contains_point_2d = .true.
    return
  end if

  t(1:2,1:3) = reshape ( (/ q(1:2,1), q(1:2,3), q(1:2,4) /), (/ 2, 3 /) )

  if ( triangle_contains_point_2d ( t, p ) ) then
    quad_contains_point_2d = .true.
    return
  end if

  quad_contains_point_2d = .false.

  return
end
function r8_modp ( x, y )

!*****************************************************************************80
!
!! R8_MODP returns the nonnegative remainder of R8 division.
!
!  Discussion:
!
!    If
!      REM = R8_MODP ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360.0) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
!
!  Example:
!
!        I         J     MOD R8_MODP  R8_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be divided.
!
!    Input, real ( kind = 8 ) Y, the number that divides X.
!
!    Output, real ( kind = 8 ) R8_MODP, the nonnegative remainder
!    when X is divided by Y.
!
  implicit none

  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( y == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MODP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R8_MODP ( X, Y ) called with Y = ', y
    stop
  end if

  r8_modp = mod ( x, y )

  if ( r8_modp < 0.0D+00 ) then
    r8_modp = r8_modp + abs ( y )
  end if

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
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
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

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
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  logical d_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i7,7x)') i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(g14.6)' ) a(i,j)

      end do

!     write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )
      write ( *, '(5x,1x,5a14)' )    ( ctemp(i), i = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_filename
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real      ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

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
!    1    0    4     2*i      - 2*i        0    0
!    2   -6    5     3 +   i    3 -   i    3    3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
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
!    Output, real ( kind = 8 ) R1, R2, the real parts of the roots
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
  r1 = -0.5D+00 * q / a
  r2 = -2.0D+00 * c / q

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
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

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

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
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) get
  integer ( kind = 4 ) put
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  put = 0
  nchar = len_trim ( s )

  do get = 1, nchar

    c = s(get:get)

    if ( c /= ' ' .and. c /= TAB ) then
      put = put + 1
      s(put:put) = c
    end if

  end do

  s(put+1:nchar) = ' '

  return
end
subroutine segment_point_near_2d ( p1, p2, p, pn, dist, t )

!*****************************************************************************80
!
!! SEGMENT_POINT_NEAR_2D finds the line segment point nearest a point in 2D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!    The nearest point will satisfy the condition
!
!      PN = (1-T) * P1 + T * P2.
!
!    T will always be between 0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), the two endpoints of the line
!    segment.  P1 should generally be different from P2, but
!    if they are equal, the program will still compute a
!    meaningful result.
!
!    Input, real ( kind = 8 ) P(2), the point whose nearest neighbor
!    on the line segment is to be determined.
!
!    Output, real ( kind = 8 ) PN(2), the point on the line segment which is
!    nearest the point P.
!
!    Output, real ( kind = 8 ) DIST, the distance from the point to the
!    nearest point on the line segment.
!
!    Output, real ( kind = 8 ) T, the relative position of the point (XN,YN)
!    to the points P1 and P2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) pn(dim_num)
  real ( kind = 8 ) t

  if ( p1(1) == p2(1) .and. p1(2) == p2(2) ) then

    t = 0.0D+00
    pn(1:dim_num) = p1(1:dim_num)

  else

    t =   ( ( p1(1) - p(1) ) * ( p1(1) - p2(1) )   &
          + ( p1(2) - p(2) ) * ( p1(2) - p2(2) ) ) &
        / ( ( p1(1) - p2(1) )**2                   &
          + ( p1(2) - p2(2) )**2 )

    if ( t < 0.0D+00 ) then
      t = 0.0D+00
      pn(1:dim_num) = p1(1:dim_num)
    else if ( 1.0D+00 < t ) then
      t = 1.0D+00
      pn(1:dim_num) = p2(1:dim_num)
    else
      pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )
    end if

  end if

  dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end
function sin_deg ( angle )

!*****************************************************************************80
!
!! SIN_DEG returns the sine of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle, in degrees.
!
!    Output, real ( kind = 8 ) SIN_DEG, the sine of the angle.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ), parameter :: degrees_to_radians &
    = 3.141592653589793D+00 / 180.0D+00
  real ( kind = 8 ) sin_deg

  sin_deg = sin ( degrees_to_radians * angle )

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
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
function triangle_contains_point_2d ( t, p )

!*****************************************************************************80
!
!! TRIANGLE_CONTAINS_POINT_2D finds if a point is inside a triangle in 2D.
!
!  Discussion:
!
!    The routine assumes that the vertices are given in counter-clockwise
!    order.
!
!    The routine determines if the point P is "to the right of" each
!    of the lines that bound the triangle.  It does this by computing
!    the cross product of vectors from a vertex to its next vertex, and to P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!    The vertices should be given in counter clockwise order.
!
!    Input, real ( kind = 8 ) P(2), the point to be checked.
!
!    Output, logical TRIANGLE_CONTAINS_POINT_2D, is TRUE if the point
!    is inside the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) t(dim_num,3)
  logical triangle_contains_point_2d

  triangle_contains_point_2d = .false.

  if ( 0.0D+00 < ( p(1)   - t(1,1) ) * ( t(2,2) - t(2,1) ) &
               - ( t(1,2) - t(1,1) ) * ( p(2)   - t(2,1) ) ) then
    return
  end if

  if ( 0.0D+00 < ( p(1)   - t(1,2) ) * ( t(2,3) - t(2,2) ) &
               - ( t(1,3) - t(1,2) ) * ( p(2)   - t(2,2) ) ) then
    return
  end if

  if ( 0.0D+00 < ( p(1)   - t(1,3) ) * ( t(2,1) - t(2,3) ) &
               - ( t(1,1) - t(1,3) ) * ( p(2)   - t(2,3) ) ) then
    return
  end if

  triangle_contains_point_2d = .true.

  return
end
