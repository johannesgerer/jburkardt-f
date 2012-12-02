program shazam

!*****************************************************************************80
!
!! SHAZAM is just a dummy program.
!
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHAZAM!'

  stop
end
block data

!*****************************************************************************80
!
!! BLOCK DATA is an untitled block data module.
!
!
  implicit none

  real x
  real y
  real z

  common / arthur / x, y, z

  data x / 1.0 /
  data y / 2.0 /
  data z / 3.0 /

  stop
end
subroutine angle_to_rgb ( angle, r, g, b )

!*****************************************************************************80
!
!! ANGLE_TO_RGB returns a color on the perimeter of the color hexagon.
!
!  Modified:
!
!    13 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, the angle in the color hexagon.  The sextants are
!    defined by the following points:
!        0 degrees, 1, 0, 0, red;
!       60 degrees, 1, 1, 0, yellow;
!      120 degrees, 0, 1, 0, green;
!      180 degrees, 0, 1, 1, cyan;
!      240 degrees, 0, 0, 1, blue;
!      300 degrees, 1, 0, 1, magenta.
!
!    Output, real R, G, B, RGB specifications for the color that lies
!    at the given angle, on the perimeter of the color hexagon.  One
!    value will be 1, and one value will be 0.
!
  implicit none

  real angle
  real angle2
  real b
  real g
  real, parameter :: degrees_to_radians = &
    3.14159265E+00 / 180.0E+00
  real r

  angle = mod ( angle, 360.0E+00 )

  if ( angle < 0.0E+00 ) then
    angle = angle + 360.0E+00
  end if

  if ( angle <= 60.0E+00 ) then

    angle2 = degrees_to_radians * 3.0E+00 * angle / 4.0E+00
    r = 1.0E+00
    g = tan ( angle2 )
    b = 0.0E+00

  else if ( angle <= 120.0E+00 ) then

    angle2 = degrees_to_radians * 3.0E+00 * angle / 4.0E+00
    r = cos ( angle2 ) / sin ( angle2 )
    g = 1.0E+00
    b = 0.0E+00

  else if ( angle <= 180.0E+00 ) then

    angle2 = degrees_to_radians * 3.0E+00 * ( angle - 120.0E+00 ) / 4.0E+00
    r = 0.0E+00
    g = 1.0E+00
    b = tan ( angle2 )

  else if ( angle <= 240.0E+00 ) then

    angle2 = degrees_to_radians * 3.0E+00 * ( angle - 120.0E+00 ) / 4.0E+00
    r = 0.0E+00
    g = cos ( angle2 ) / sin ( angle2 )
    b = 1.0E+00

  else if ( angle <= 300.0E+00 ) then

    angle2 = degrees_to_radians * 3.0E+00 * ( angle - 240.0E+00 ) / 4.0E+00
    r = tan ( angle2 )
    g = 0.0E+00
    b = 1.0E+00

  else if ( angle <= 360.0E+00 ) then

    angle2 = degrees_to_radians * 3.0E+00 * ( angle - 240.0E+00 ) / 4.0E+00
    r = 1.0E+00
    g = 0.0E+00
    b = cos ( angle2 ) / sin ( angle2 )

  end if

  return
end

!*****************************************************************************80
!
!! MAIN is a program with no PROGRAM header.
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MAIN!'

  stop
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
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real Y, X, two quantities which represent the tangent of
!    an angle.  If Y is not zero, then the tangent is (Y/X).
!
!    Output, real ATAN4, an angle between 0 and 2 * PI, whose tangent is
!    (Y/X), and which lies in the appropriate quadrant so that the signs
!    of its cosine and sine match those of X and Y.
!
  implicit none

  real abs_x
  real abs_y
  real atan4
  real, parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510E+00
  real theta
  real theta_0
  real x
  real y
!
!  Special cases:
!
  if ( x == 0.0E+00 ) then

    if ( 0.0E+00 < y ) then
      theta = pi / 2.0E+00
    else if ( y < 0.0E+00 ) then
      theta = 3.0E+00 * pi / 2.0E+00
    else if ( y == 0.0E+00 ) then
      theta = 0.0E+00
    end if

  else if ( y == 0.0E+00 ) then

    if ( 0.0E+00 < x ) then
      theta = 0.0E+00
    else if ( x < 0.0E+00 ) then
      theta = pi
    end if
!
!  We assume that ATAN2 is correct when both arguments are positive.
!
  else

    abs_y = abs ( y )
    abs_x = abs ( x )

    theta_0 = atan2 ( abs_y, abs_x )

    if ( 0.0E+00 < x .and. 0.0E+00 < y ) then
      theta = theta_0
    else if ( x < 0.0E+00 .and. 0.0E+00 < y ) then
      theta = pi - theta_0
    else if ( x < 0.0E+00 .and. y < 0.0E+00 ) then
      theta = pi + theta_0
    else if ( 0.0E+00 < x .and. y < 0.0E+00 ) then
      theta = 2.0E+00 * pi - theta_0
    end if

  end if

  atan4 = theta

  return
end
blockdata barbara

!*****************************************************************************80
!
!! BARBARA is a block data module.
!
!
  implicit none

  real u
  real v
  real w

  common / clarke / u, v, w

  data u / 1.0 /
  data v / 2.0 /
  data w / 3.0 /

  stop
end
subroutine axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, nticks )

!*****************************************************************************80
!
!! AXIS_LIMITS returns "nice" axis limits for a plot.
!
!  Discussion:
!
!    The routine is given information about the range of a variable, and
!    the number of divisions desired.  It returns suggestions for
!    labeling a plotting axis for the variable, including the
!    starting and ending points, the length of a single division,
!    and a suggested tick marking for the axis.
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
!    Input, real XMIN, XMAX, the lower and upper values that must be
!    included on the axis.  XMIN must be less than XMAX.
!
!    Input, integer NDIVS, the number of divisions desired along
!    the axis.
!
!    Output, real PXMIN, PXMAX, the recommended lower and upper axis
!    bounds.  It will be the case that PXMIN <= XMIN < XMAX <= PXMAX.
!
!    Output, real PXDIV, the recommended size of a single division.
!
!    Output, integer NTICKS, a suggested number of ticks to use,
!    if subdividing each of the NDIVS divisions of the axis.
!
  implicit none

  integer, parameter :: nsteps = 5

  real best
  real good
  integer i
  integer ihi
  integer ilo
  integer intlog
  integer iticks(5)
  integer ival
  integer j
  integer ndivs
  integer nticks
  real pxmax
  real pxmax2
  real pxmin
  real pxmin2
  real pxdiv
  real pxdiv2
  real r_log_10
  real reldif
  real steps(nsteps)
  real temp
  real xmax
  real xmin

  if ( xmin == xmax ) then
    xmin = xmin - 0.5E+00
    xmax = xmax + 0.5E+00
  else if ( xmax < xmin ) then
    temp = xmin
    xmin = xmax
    xmax = temp
  end if

  if ( ndivs <= 0 ) then
    ndivs = 5
  end if

  steps(1) =  1.0E+00
  steps(2) =  2.0E+00
  steps(3) =  4.0E+00
  steps(4) =  5.0E+00
  steps(5) = 10.0E+00

  iticks(1) = 5
  iticks(2) = 4
  iticks(3) = 4
  iticks(4) = 5
  iticks(5) = 5
!
!  Set RELDIF, the size of the X interval divided by the largest X.
!
  if ( xmax /= xmin ) then
    reldif = ( xmax - xmin ) / max ( abs ( xmax ), abs ( xmin ) )
  else
    reldif = 0.0E+00
  end if
!
!  If RELDIF tells us that XMIN and XMAX are extremely close,
!  do some simple things.
!
  if ( reldif < 0.00001E+00 ) then

    if ( xmax == 0.0E+00 ) then

      pxdiv = 1.0E+00

    else

      intlog = int ( r_log_10 ( xmax ) )

      if ( intlog < 0 ) then
        intlog = intlog - 1
      end if

      pxdiv = 10.0E+00**intlog

      if ( 1.0E+00 < pxdiv ) then
        pxdiv = 1.0E+00
      end if

    end if

    nticks = 5
    pxmin = xmax - real ( ndivs / 2 ) * pxdiv
    pxmax = xmax + real ( ndivs - ( ndivs / 2 ) ) * pxdiv
!
!  But now handle the more general case, when XMIN and XMAX
!  are relatively far apart.
!
  else

    best = - 999.0E+00
!
!  On second loop, increase INTLOG by 1.
!
    do j = 1, 2
!
!  Compute INTLOG, roughly the logarithm base 10 of the range
!  divided by the number of divisions.
!
      intlog = int ( r_log_10 ( ( xmax - xmin ) / real ( ndivs ) ) ) + ( j - 1 )

      if ( xmax - xmin  < real ( ndivs ) ) then
        intlog = intlog - 1
      end if
!
!  Now consider taking 1, 2, 4, 5 or 10 steps of size 10**INTLOG:
!
      do i = 1, nsteps
!
!  Compute the size of each step.
!
        pxdiv2 = steps(i) * 10.0E+00**intlog
!
!  Make sure NDIVS steps can reach from XMIN to XMAX, at least.
!
        if ( xmax <= xmin + ndivs * pxdiv2 ) then
!
!  Now decide where to start the axis.
!  Start the axis at PXMIN2, to the left of XMIN, and
!  representing a whole number of steps of size PXDIV2.
!
          if ( 0.0E+00 <= xmin ) then
            ival = int ( xmin / pxdiv2 )
          else
            ival = int ( xmin / pxdiv2 ) - 1
          end if

          pxmin2 = ival * pxdiv2
!
!  PXMAX2 is, of course, NDIVS steps above PXMIN2.
!
          pxmax2 = pxmin2 + ndivs * pxdiv2
!
!  Only consider going on if PXMAX2 is at least XMAX.
!
          if ( xmax <= pxmax2 ) then
!
!  Now judge this grid by the relative amount of wasted axis length.
!
            good = ( xmax - xmin ) / ( pxmax2 - pxmin2 )

            if ( best < good ) then
              best = good
              pxmax = pxmax2
              pxmin = pxmin2
              pxdiv = pxdiv2
              nticks = iticks(i)
            end if

          end if

        end if

      end do

    end do

  end if
!
!  If necessary, adjust the locations of PXMIN and PXMAX so that the
!  interval is more symmetric in containing XMIN through XMAX.
!
  do

    ilo = int ( xmin - pxmin ) / pxdiv
    ihi = int ( pxmax - xmax ) / pxdiv

    if ( ihi < ilo + 2 ) then
      exit
    end if

    pxmin = pxmin - pxdiv
    pxmax = pxmax - pxdiv

  end do

  return
end
block data charlie

!*****************************************************************************80
!
!! CHARLIE is a block data routine.
!
!
  implicit none

  real x
  real y
  real z

  common / doyle / x, y, z

  data x / 1.0 /
  data y / 2.0 /
  data z / 3.0 /

  stop
end
integer function bit_hi1_base_2 ( n )

!*****************************************************************************80
!
!! BIT_HI1_BASE_2 returns the position of the high 1 bit base 2 in an integer.
!
!  Example:
!
!       N    Binary    Hi 1
!    ----    --------  ----
!       0           0     0
!       1           1     1
!       2          10     2
!       3          11     2
!       4         100     3
!       5         101     3
!       6         110     3
!       7         111     3
!       8        1000     4
!       9        1001     4
!      10        1010     4
!      11        1011     4
!      12        1100     4
!      13        1101     4
!      14        1110     4
!      15        1111     4
!      16       10000     5
!      17       10001     5
!    1023  1111111111    10
!    1024 10000000000    11
!    1025 10000000001    11
!
!  Modified:
!
!    13 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the integer to be measured.
!    N should be nonnegative.  If N is nonpositive, BIT_HI1_BASE_2
!    will always be 0.
!
!    Output, integer BIT_HI1_BASE_2, the number of bits base 2.
!
  implicit none

  integer bit
  integer i
  integer n

  i = n
  bit = 0

  do

    if ( i <= 0 ) then
      exit
    end if

    bit = bit + 1
    i = i / 2

  end do

  bit_hi1_base_2 = bit

  return
end
real function bmi_english ( w_lb, h_ft, h_in )

!*****************************************************************************80
!
!! BMI_ENGLISH computes the body mass index given English measurements.
!
!  Modified:
!
!    20 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real W_LB, the body weight in pounds.
!
!    Input, real H_FT, H_IN, the body height in feet and inches
!
!    Output, real BMI_ENGLISH, the body mass index.
!
  implicit none

  real bmi_metric
  real feet_to_meters
  real h_ft
  real h_in
  real h_m
  real pounds_to_kilograms
  real w_kg
  real w_lb

  w_kg = pounds_to_kilograms ( w_lb )

  h_m = feet_to_meters ( h_ft + ( h_in / 12.0E+00 ) )

  bmi_english = bmi_metric ( w_kg, h_m )

  return
end
double precision function d_pi ( )

!*****************************************************************************80
!
!! D_PI returns the value of pi as a double precision quantity.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision D_PI, the value of pi.
!
  implicit none

  d_pi = 3.14159265358979323846264338327950288419716939937510D+00

  return
end
complex function c_cube_root ( x )

!*****************************************************************************80
!
!! C_CUBE_ROOT returns the principal cube root of a complex number.
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex X, the number whose cube root is desired.
!
!    Output, complex C_CUBE_ROOT, the cube root of X.
!
  implicit none

  real a
  real atan4
  real b
  real mag
  real theta
  complex x

  a = real ( x )
  b = imag ( x )
  mag = sqrt ( a * a + b * b )

  if ( mag == 0.0E+00 ) then

    c_cube_root = cmplx ( 0.0E+00, 0.0E+00 )

  else

    theta = atan4 ( b, a )
    c_cube_root = mag**( 1.0E+00 / 3.0E+00 ) &
      * cmplx ( cos ( theta / 3.0E+00 ), sin ( theta / 3.0E+00 ) )

  end if

  return
end
