subroutine angle_to_rgb ( angle, r, g, b )

!*****************************************************************************80
!
!! ANGLE_TO_RGB returns a color on the perimeter of the color hexagon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, real ( kind = 8 ) ANGLE, the angle in the color hexagon.
!    The sextants are defined by the following points:
!        0 degrees, 1, 0, 0, red;
!       60 degrees, 1, 1, 0, yellow;
!      120 degrees, 0, 1, 0, green;
!      180 degrees, 0, 1, 1, cyan;
!      240 degrees, 0, 0, 1, blue;
!      300 degrees, 1, 0, 1, magenta.
!
!    Output, real ( kind = 8 ) R, G, B, RGB specifications for the color
!    that lies at the given angle, on the perimeter of the color hexagon.  One
!    value will be 1, and one value will be 0.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle2
  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ), parameter :: degrees_to_radians = &
    3.141592653589793D+00 / 180.0D+00
  real ( kind = 8 ) r

  angle = mod ( angle, 360.0D+00 )

  if ( angle < 0.0D+00 ) then
    angle = angle + 360.0D+00
  end if

  if ( angle <= 60.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = 1.0D+00
    g = tan ( angle2 )
    b = 0.0D+00

  else if ( angle <= 120.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = cos ( angle2 ) / sin ( angle2 )
    g = 1.0D+00
    b = 0.0D+00

  else if ( angle <= 180.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = 1.0D+00
    b = tan ( angle2 )

  else if ( angle <= 240.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = cos ( angle2 ) / sin ( angle2 )
    b = 1.0D+00

  else if ( angle <= 300.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = tan ( angle2 )
    g = 0.0D+00
    b = 1.0D+00

  else if ( angle <= 360.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = 1.0D+00
    g = 0.0D+00
    b = cos ( angle2 ) / sin ( angle2 )

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
!    14 April 1999
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
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine chart_xyz_cap ( xcap, ycap, zcap, color )

!*****************************************************************************80
!
!! CHART_XYZ_CAP returns the CIE XYZ values of a 24 box color chart.
!
!  Discussion:
!
!    The chart may be drawn as a set of 4 rows of 6 squares of color:
!
!     1  2  3  4  5  6
!     7  8  9 10 11 12
!    13 14 15 16 17 18
!    19 20 21 22 23 24
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Roy Hall,
!    Illumination and Color in Computer Generated Imagery,
!    Springer Verlag, 1988.
!
!    Calvin McCamy, H Marcus and J G Davidson,
!    A Color Rendition Chart,
!    Journal of Applied Photographic Engineering,
!    Volume 11, number 3, pages 95-99.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) XCAP(24), YCAP(24), ZCAP(24), the CIE XYZ
!    color coordinates of the color squares.
!
!    Output, character ( len = * ) COLOR(24), the names of the colors.
!    The names are up to 14 characters in length.
!
  implicit none

  character ( len = * ) color(24)
  real ( kind = 8 ) xcap(24)
  real ( kind = 8 ) ycap(24)
  real ( kind = 8 ) zcap(24)

  color(1) = 'dark skin'
  xcap(1) = 0.092D+00
  ycap(1) = 0.081D+00
  zcap(1) = 0.058D+00

  color(2) = 'light skin'
  xcap(2) = 0.411D+00
  ycap(2) = 0.376D+00
  zcap(2) = 0.303D+00

  color(3) = 'blue sky'
  xcap(3) = 0.183D+00
  ycap(3) = 0.186D+00
  zcap(3) = 0.373D+00

  color(4) = 'foliage'
  xcap(4) = 0.094D+00
  ycap(4) = 0.117D+00
  zcap(4) = 0.067D+00

  color(5) = 'blue flower'
  xcap(5) = 0.269D+00
  ycap(5) = 0.244D+00
  zcap(5) = 0.503D+00

  color(6) = 'bluish green'
  xcap(6) = 0.350D+00
  ycap(6) = 0.460D+00
  zcap(6) = 0.531D+00

  color(7) = 'orange'
  xcap(7) = 0.386D+00
  ycap(7) = 0.311D+00
  zcap(7) = 0.066D+00

  color(8) = 'purplish blue'
  xcap(8) = 0.123D+00
  ycap(8) = 0.102D+00
  zcap(8) = 0.359D+00

  color(9) = 'moderate red'
  xcap(9) = 0.284D+00
  ycap(9) = 0.192D+00
  zcap(9) = 0.151D+00

  color(10) = 'purple'
  xcap(10) = 0.059D+00
  ycap(10) = 0.040D+00
  zcap(10) = 0.102D+00

  color(11) = 'yellow green'
  xcap(11) = 0.368D+00
  ycap(11) = 0.474D+00
  zcap(11) = 0.127D+00

  color(12) = 'orange yellow'
  xcap(12) = 0.497D+00
  ycap(12) = 0.460D+00
  zcap(12) = 0.094D+00

  color(13) = 'blue'
  xcap(13) = 0.050D+00
  ycap(13) = 0.035D+00
  zcap(13) = 0.183D+00

  color(14) = 'green'
  xcap(14) = 0.149D+00
  ycap(14) = 0.234D+00
  zcap(14) = 0.106D+00

  color(15) = 'red'
  xcap(15) = 0.176D+00
  ycap(15) = 0.102D+00
  zcap(15) = 0.048D+00

  color(16) = 'yellow'
  xcap(16) = 0.614D+00
  ycap(16) = 0.644D+00
  zcap(16) = 0.112D+00

  color(17) = 'magenta'
  xcap(17) = 0.300D+00
  ycap(17) = 0.192D+00
  zcap(17) = 0.332D+00

  color(18) = 'cyan'
  xcap(18) = 0.149D+00
  ycap(18) = 0.192D+00
  zcap(18) = 0.421D+00

  color(19) = 'white'
  xcap(19) = 0.981D+00
  ycap(19) = 1.000D+00
  zcap(19) = 1.184D+00

  color(20) = 'neutral 8'
  xcap(20) = 0.632D+00
  ycap(20) = 0.644D+00
  zcap(20) = 0.763D+00

  color(21) = 'neutral 6.5'
  xcap(21) = 0.374D+00
  ycap(21) = 0.381D+00
  zcap(21) = 0.451D+00

  color(22) = 'neutral 5'
  xcap(22) = 0.189D+00
  ycap(22) = 0.192D+00
  zcap(22) = 0.227D+00

  color(23) = 'neutral 3.5'
  xcap(23) = 0.067D+00
  ycap(23) = 0.068D+00
  zcap(23) = 0.080D+00

  color(24) = 'black'
  xcap(24) = 0.000D+00
  ycap(24) = 0.000D+00
  zcap(24) = 0.000D+00

  return
end
subroutine cmy_check ( c, m, y )

!*****************************************************************************80
!
!! CMY_CHECK corrects out-of-range CMY color coordinates.
!
!  Discussion:
!
!    The CMY color system describes a color based on the amounts of the
!    base colors cyan, magenta, and yellow.  Thus, a particular color
!    has three coordinates, (C,M,Y).  Each coordinate must be between
!    0 and 1.
!
!  Example:
!
!    Black   (1,1,1)
!    Blue    (1,1,0)
!    Cyan    (1,0,0)
!    Green   (1,0,1)
!    Magenta (0,1,0)
!    Red     (0,1,1)
!    White   (0,0,0)
!    Yellow  (0,0,1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) C, M, Y, the CMY color coordinates
!    to be checked.  Any value less than 0 is increased to zero.
!    Any value greater than 1 is decreased to 1.
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) m
  real ( kind = 8 ) y

  c = max ( c, 0.0D+00 )
  c = min ( c, 1.0D+00 )

  m = max ( m, 0.0D+00 )
  m = min ( m, 1.0D+00 )

  y = max ( y, 0.0D+00 )
  y = min ( y, 1.0D+00 )

  return
end
subroutine cmy_to_cmyk ( c, m, y, c2, m2, y2, k2 )

!*****************************************************************************80
!
!! CMY_TO_CMYK converts CMY to CMYK color coordinates.
!
!  Discussion:
!
!    The CMY color system describes a color based on the amounts of the
!    base colors cyan, magenta, and yellow.  Thus, a particular color
!    has three coordinates, (C,M,Y).  Each coordinate must be between
!    0 and 1.  Black is (1,1,1) and white is (0,0,0).
!
!    The CMYK color system describes a color based on the amounts of the
!    base colors cyan, magenta, yellow, and black.  The CMYK system is
!    based on the CMY system, except that equal amounts of C, M, and Y
!    are replaced by the single color K.  Thus, a particular color
!    has four coordinates, (C,M,Y,K).  Each coordinate must be between
!    0 and 1, and it must also be true that C+K, M+K and Y+K are
!    each no greater than 1.
!
!    K2 = min ( C, M, Y )
!    C2 = C - K2
!    M2 = M - K2
!    Y2 = Y - K2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, M, Y, the CMY color coordinates
!    to be converted.
!
!    Output, real ( kind = 8 ) C2, M2, Y2, K2, the corresponding CMYK
!    color coordinates.
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) c2
  real ( kind = 8 ) k2
  real ( kind = 8 ) m
  real ( kind = 8 ) m2
  real ( kind = 8 ) y
  real ( kind = 8 ) y2

  k2 = min ( c, m, y )
  c2 = c - k2
  m2 = m - k2
  y2 = y - k2

  return
end
subroutine cmy_to_rgb ( c, m, y, r, g, b )

!*****************************************************************************80
!
!! CMY_TO_RGB converts CMY to RGB color coordinates.
!
!  Discussion:
!
!    The CMY color system describes a color based on the amounts of the
!    base colors cyan, magenta, and yellow.  Thus, a particular color
!    has three coordinates, (C,M,Y).  Each coordinate must be between
!    0 and 1.  Black is (1,1,1) and white is (0,0,0).
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!    R = 1.0D+00 - C
!    G = 1.0D+00 - M
!    B = 1.0D+00 - Y
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, M, Y, the CMY color coordinates
!    to be converted.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) g
  real ( kind = 8 ) m
  real ( kind = 8 ) r
  real ( kind = 8 ) y

  r = 1.0D+00 - c
  g = 1.0D+00 - m
  b = 1.0D+00 - y

  return
end
subroutine cmyk_check ( c, m, y, k )

!*****************************************************************************80
!
!! CMYK_CHECK corrects out-of-range CMYK color coordinates.
!
!  Discussion:
!
!    The CMYK color system describes a color based on the amounts of the
!    base colors cyan, magenta, yellow, and black.  The CMYK system is
!    based on the CMY system, except that equal amounts of C, M, and Y
!    are replaced by the single color K.  Thus, a particular color
!    has four coordinates, (C,M,Y,K).  Each coordinate must be between
!    0 and 1, and it must also be true that C+K, M+K and Y+K are
!    each no greater than 1.
!
!  Examples:
!
!    Black   (0,0,0,1)
!    Blue    (1,1,0,0)
!    Cyan    (1,0,0,0)
!    Green   (1,0,1,0)
!    Magenta (0,1,0,0)
!    Red     (0,1,1,0)
!    White   (0,0,0,0)
!    Yellow  (0,0,1,0)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) C, M, Y, K, the CMYK color coordinates
!    to be checked.
!    Any value less than 0 is increased to zero.
!    Any value greater than 1 is decreased to 1.
!    Then, if any of C+K, M+K or Y+K is greater than 1, C, M or Y is reduced
!    accordingly.
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) k
  real ( kind = 8 ) m
  real ( kind = 8 ) y
!
!  1: Enforce the simple rule that C, M, Y and K must lie between 0 and 1.
!
  c = max ( c, 0.0D+00 )
  c = min ( c, 1.0D+00 )

  m = max ( m, 0.0D+00 )
  m = min ( m, 1.0D+00 )

  y = max ( y, 0.0D+00 )
  y = min ( y, 1.0D+00 )

  k = max ( k, 0.0D+00 )
  k = min ( k, 1.0D+00 )
!
!  2: Enforce C+K, M+K, Y+K each no greater than 1.
!
  c = min ( c, 1.0D+00 - k )
  m = min ( m, 1.0D+00 - k )
  y = min ( y, 1.0D+00 - k )

  return
end
subroutine cmyk_to_cmy ( c, m, y, k, c2, m2, y2 )

!*****************************************************************************80
!
!! CMYK_TO_CMY converts CMYK to CMY color coordinates.
!
!  Discussion:
!
!    The CMYK color system describes a color based on the amounts of the
!    base colors cyan, magenta, yellow, and black.  The CMYK system is
!    based on the CMY system, except that equal amounts of C, M, and Y
!    are replaced by the single color K.  Thus, a particular color
!    has four coordinates, (C,M,Y,K).  Each coordinate must be between
!    0 and 1, and it must also be true that C+K, M+K and Y+K are
!    each no greater than 1.
!
!    The CMY color system describes a color based on the amounts of the
!    base colors cyan, magenta, and yellow.  Thus, a particular color
!    has three coordinates, (C,M,Y).  Each coordinate must be between
!    0 and 1.  Black is (1,1,1) and white is (0,0,0).
!
!    C2 = C + K
!    M2 = M + K
!    Y2 = Y + K
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, M, Y, K, the CMYK color coordinates
!    to be converted.
!
!    Output, real ( kind = 8 ) C2, M2, Y2, the corresponding CMY color
!    coordinates.
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) c2
  real ( kind = 8 ) k
  real ( kind = 8 ) m
  real ( kind = 8 ) m2
  real ( kind = 8 ) y
  real ( kind = 8 ) y2

  c2 = c + k
  m2 = m + k
  y2 = y + k

  return
end
subroutine cmyk_to_rgb ( c, m, y, k, r, g, b )

!*****************************************************************************80
!
!! CMYK_TO_RGB converts CMYK to RGB color coordinates.
!
!  Discussion:
!
!    The CMYK color system describes a color based on the amounts of the
!    base colors cyan, magenta, yellow, and black.  The CMYK system is
!    based on the CMY system, except that equal amounts of C, M, and Y
!    are replaced by the single color K.  Thus, a particular color
!    has four coordinates, (C,M,Y,K).  Each coordinate must be between
!    0 and 1, and it must also be true that C+K, M+K and Y+K are
!    each no greater than 1.
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, M, Y, K, the CMYK color coordinates
!    to be converted.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) g
  real ( kind = 8 ) k
  real ( kind = 8 ) m
  real ( kind = 8 ) r
  real ( kind = 8 ) y

  r = 1.0D+00 - ( c + k )
  g = 1.0D+00 - ( m + k )
  b = 1.0D+00 - ( y + k )

  return
end
subroutine get_seed ( seed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
!
  implicit none

  integer   ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer   ( kind = 4 ) seed
  real      ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 )  today
  integer   ( kind = 4 ) values(8)
  character ( len = 5 )  zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) /  11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) /  30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) /  23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp                                    /   6.0D+00

  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( i4_huge, kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == i4_huge ) then
    seed = seed - 1
  end if

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
subroutine grayscale_luv ( n, l, u, v )

!*****************************************************************************80
!
!! GRAYSCALE_LUV returns a grayscale in the CIE LUV system.
!
!  Discussion:
!
!    Only the L (lightness) coordinate does anything interesting.
!    The U and V coordinates will all be 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of levels in the gray scale.
!
!    Output, real ( kind = 8 ) L(N), U(N), V(N), the CIE LUV color coordinates
!    of the gray scale.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) l(n)
  real ( kind = 8 ), parameter :: l_hi = 100.0D+00
  real ( kind = 8 ), parameter :: l_lo =   0.0D+00
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)

  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then

    l(1) = ( l_lo + l_hi ) / 2.0D+00

  else

    do i = 1, n

      l(i) = ( real ( n - i,     kind = 8 ) * l_lo   &
             + real (     i - 1, kind = 8 ) * l_hi ) &
             / real ( n     - 1, kind = 8 )

    end do

  end if

  u(1:n) = 0.0D+00
  v(1:n) = 0.0D+00

  return
end
subroutine grayscale_rgb ( n, r, g, b )

!*****************************************************************************80
!
!! GRAYSCALE_RGB returns a grayscale in the RGB system.
!
!  Discussion:
!
!    All we do here is make a linear gray scale, and then run it
!    through the gamma correction.  I have no evidence that this
!    is the right thing to do.  My hope is that the grays will
!    seem more equally spaced.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of levels in the gray scale.
!
!    Output, real ( kind = 8 ) R(N), G(N), B(N), the RGB color coordinates
!    of the gray scale.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) g(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) v
  real ( kind = 8 ), parameter :: v_hi = 1.0D+00
  real ( kind = 8 ), parameter :: v_lo = 0.0D+00
  real ( kind = 8 ) v_prime

  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then

    v = ( v_lo + v_hi ) / 2.0D+00
    call lin_to_nonlin ( v, v_prime )
    r(1) = v_prime

  else

    do i = 1, n

      v = ( real ( n - i,     kind = 8 ) * v_lo   &
          + real (     i - 1, kind = 8 ) * v_hi ) &
          / real ( n     - 1, kind = 8 )

      call lin_to_nonlin ( v, v_prime )

      r(i) = v_prime

    end do

  end if

  g(1:n) = r(1:n)
  b(1:n) = r(1:n)

  return
end
subroutine hls_check ( h, l, s )

!*****************************************************************************80
!
!! HLS_CHECK corrects out-of-range HLS color coordinates.
!
!  Discussion:
!
!    The HLS color system describes a color based on the qualities of
!    hue, lightness, and saturation.  A particular color has three
!    coordinates, (H,L,S).  The L and S coordinates must be between
!    0 and 1, while the H coordinate must be between 0 and 360, and
!    is interpreted as an angle.
!
!    The HLS color space is usually thought of as a double hexcone.
!    If the L coordinate is vertical, then the color space is a single
!    black point at L = 0, expands to a colorful hexagon
!    at L = 0.5, and contracts again to a white point at L = 1.
!    The colorful hexagon as the colors Red, Yellow, Green, Cyan, Blue,
!    and Magenta at its vertices.
!    The saturation coordinate varies from 0.0D+00 at the center of the
!    hexagon to 1.0D+00 at the boundary.  The corresponding color varies
!    from gray at S = 0 to the full color on the boundary at S = 1.
!
!    If the (H,S) plane is thought of as a circle, then S is the relative
!    distance from the central vertical L axis to the boundary of the hexcone.
!    Thus, even as the cone contracts to a point, S can always vary
!    from 0 to 1.  In particular, the white point can have coordinates
!    ( H, 1.0, S ) where H is any value in [ 0, 360.0D+00 ) and S is
!    any value in [ 0.0, 1.0D+00 ].
!
!
!    Some versions of the HLS model assign Blue to have HLS coordinates
!    ( 0.0, 0.5, 1.0D+00 ) instead of Red, rotating all the colors ahead
!    120 degrees.
!
!    Some versions define L as equal to the average of R, G and B, rather
!    than the average of the maximum and minimum of R, G and B.
!
!  Example:
!
!    Given RED = ( 0.0, 0.5, 1.0D+00 ) then:
!
!      Black   (   0.0, 0.0, 0.0D+00 )
!      Blue    ( 240.0, 0.5, 1.0D+00 )
!      Cyan    ( 180.0, 0.5, 1.0D+00 )
!      Green   ( 120.0, 0.5, 1.0D+00 )
!      Magenta ( 300.0, 0.5, 1.0D+00 )
!      Red     (   0.0, 0.5, 1.0D+00 )
!      White   (   0.0, 1.0, 0.0D+00 )
!      Yellow  (  60.0, 0.5, 1.0D+00 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) H, L, S, the HLS color coordinates to
!    be checked.  If H is outside the range [0, 360), it is brought back in
!    range using a MOD operation.
!    Values of L or S less than 0 are set to 0, greater than 1 are set
!    to 1.
!
  implicit none

  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) h
  real ( kind = 8 ) l
  real ( kind = 8 ) s

  h = r8_modp ( h, 360.0D+00 )

  l = max ( l, 0.0D+00 )
  l = min ( l, 1.0D+00 )

  s = max ( s, 0.0D+00 )
  s = min ( s, 1.0D+00 )

  return
end
subroutine hls_to_rgb ( h, l, s, r, g, b )

!*****************************************************************************80
!
!! HLS_TO_RGB converts HLS to RGB color coordinates.
!
!  Discussion:
!
!    The HLS color system describes a color based on the qualities of
!    hue, lightness, and saturation.  A particular color has three
!    coordinates, (H,L,S).  The L and S coordinates must be between
!    0 and 1, while the H coordinate must be between 0 and 360, and
!    is interpreted as an angle.
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) H, L, S, the HLS color coordinates to
!    be converted.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  real ( kind = 8 ) hls_value
  real ( kind = 8 ) l
  real ( kind = 8 ) m1
  real ( kind = 8 ) m2
  real ( kind = 8 ) r
  real ( kind = 8 ) s

  if ( l <= 0.5D+00 ) then
    m2 = l + l * s
  else
    m2 = l + s - l * s
  end if

  m1 = 2.0D+00 * l - m2

  if ( s == 0.0D+00 ) then
    r = l
    g = l
    b = l
  else
    r = hls_value ( m1, m2, h + 120.0D+00 )
    g = hls_value ( m1, m2, h )
    b = hls_value ( m1, m2, h - 120.0D+00 )
  end if

  return
end
function hls_value ( n1, n2, h )

!*****************************************************************************80
!
!! HLS_VALUE is a utility function used by HLS_TO_RGB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) N1, N2, H.
!
!    Output, real ( kind = 8 ) HLS_VALUE.
!
  implicit none

  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) h
  real ( kind = 8 ) hls_value
  real ( kind = 8 ) hue
  real ( kind = 8 ) n1
  real ( kind = 8 ) n2
!
!  Make sure HUE lies between 0 and 360.
!
  hue = r8_modp ( h, 360.0D+00 )

  if ( hue < 60.0D+00 ) then
    hls_value = n1 + ( n2 - n1 ) * hue / 60.0D+00
  else if ( hue < 180.0D+00 ) then
    hls_value = n2
  else if ( hue < 240.0D+00 ) then
    hls_value = n1 + ( n2 - n1 ) * ( 240.0D+00 - hue ) / 60.0D+00
  else
    hls_value = n1
  end if

  return
end
subroutine hsi_to_rgb ( h, s, i, r, g, b )

!*****************************************************************************80
!
!! HSI_TO_RGB converts HSI to RGB color coordinates.
!
!  Discussion:
!
!    The HSI color system uses coordinates of
!      Hue, an angle between 0 and 360, (0=R, 120=G, 240=B)
!      Saturation, between 0 and 1, and
!      Intensity, between 0 and 1.
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) H, S, I, the HSI color coordinates to be
!    converted.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
  implicit none

  real ( kind = 8 ) r8_modp
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  real ( kind = 8 ) h2
  real ( kind = 8 ) i
  real ( kind = 8 ) r

  real ( kind = 8 ) s
!
!  Ensure the H2 lies between 0 and 360.
!
  h2 = r8_modp ( h, 360.0D+00 )
!
!  Now turn H2 into a fraction of a full rotation, between 0 and 1.
!
  h2 = h2 / 360.0D+00

  if ( 3.0D+00 * h2 < 1.0D+00 ) then

    h2 = 2.0D+00 * pi * h2
    b = ( 1.0D+00 - s ) / 3.0D+00
    r = 1.0D+00 + s * cos ( h2 ) / cos ( pi / 3.0D+00 - h2 )
    r = r / 3.0D+00
    g = 1.0D+00 - b - r

  else if ( 3.0D+00 * h2 < 2.0D+00 ) then

    h2 = 2.0D+00 * pi * ( h2 - 1.0D+00 / 3.0D+00 )
    r = ( 1.0D+00 - s ) / 3
    g = 1.0D+00 + s * cos ( h2 ) / cos ( pi / 3.0D+00 - h2 )
    g = g / 3.0D+00
    b = 1.0D+00 - r - g

  else

    h2 = 2.0D+00 * pi * ( h2 - 2.0D+00 / 3.0D+00 )
    g = ( 1.0D+00 - s ) / 3.0D+00
    b = 1.0D+00 + s * cos ( h2 ) / cos ( pi / 3.0D+00 - h2 )
    b = b / 3.0D+00
    r = 1.0D+00 - g - b

  end if

  r = r * i * 3.0D+00
  g = g * i * 3.0D+00
  b = b * i * 3.0D+00

  return
end
subroutine hsv_check ( h, s, v )

!*****************************************************************************80
!
!! HSV_CHECK corrects out-of-range HSV color coordinates.
!
!  Discussion:
!
!    The HSV color system describes a color based on the three qualities
!    of hue, saturation, and value.  A given color will be represented
!    by three numbers, (H,S,V).  H, the value of hue, is an angle
!    between 0 and 360 degrees, with 0 representing red.  S is the
!    saturation, and is between 0 and 1.  Finally, V is the "value",
!    a measure of brightness, which goes from 0 for black, increasing
!    to a maximum of 1 for the brightest colors.  The HSV color system
!    is sometimes also called HSB, where the B stands for brightness.
!
!  Example:
!
!    Black   (   0.000   0.000   0.000 )
!    Blue    ( 240.000   1.000   1.000 )
!    Cyan    ( 180.000   1.000   1.000 )
!    Green   ( 120.000   1.000   1.000 )
!    Magenta ( 300.000   1.000   1.000 )
!    Red     (   0.000   1.000   1.000 )
!    White   (   0.000   0.000   1.000 )
!    Yellow  (  60.000   1.000   1.000 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) H, S, V, the HSV color coordinates to
!    be checked.
!
  implicit none

  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) h
  real ( kind = 8 ) s
  real ( kind = 8 ) v

  h = r8_modp ( h, 360.0D+00 )

  s = max ( s, 0.0D+00 )
  s = min ( s, 1.0D+00 )

  v = max ( v, 0.0D+00 )
  v = min ( v, 1.0D+00 )

  return
end
subroutine hsv_to_rgb ( h, s, v, r, g, b )

!*****************************************************************************80
!
!! HSV_TO_RGB converts HSV to RGB color coordinates.
!
!  Discussion:
!
!    The HSV color system describes a color based on the three qualities
!    of hue, saturation, and value.  A given color will be represented
!    by three numbers, (H,S,V).  H, the value of hue, is an angle
!    between 0 and 360 degrees, with 0 representing red.  S is the
!    saturation, and is between 0 and 1.  Finally, V is the "value",
!    a measure of brightness, which goes from 0 for black, increasing
!    to a maximum of 1 for the brightest colors.  The HSV color system
!    is sometimes also called HSB, where the B stands for brightness.
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) H, S, V, the HSV color coordinates to
!    be converted.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  real ( kind = 8 ) hue
  integer ( kind = 4 ) i
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) v

  if ( s == 0.0D+00 ) then

    r = v
    g = v
    b = v

  else
!
!  Make sure HUE lies between 0 and 360.0D+00
!
    hue = r8_modp ( h, 360.0D+00 )

    hue = hue / 60.0D+00

    i = int ( hue )
    f = hue - real ( i, kind = 8 )
    p = v * ( 1.0D+00 - s )
    q = v * ( 1.0D+00 - s * f )
    t = v * ( 1.0D+00 - s + s * f )

    if ( i == 0 ) then
      r = v
      g = t
      b = p
    else if ( i == 1 ) then
      r = q
      g = v
      b = p
    else if ( i == 2 ) then
      r = p
      g = v
      b = t
    else if ( i == 3 ) then
      r = p
      g = q
      b = v
    else if ( i == 4 ) then
      r = t
      g = p
      b = v
    else if ( i == 5 ) then
      r = v
      g = p
      b = q
    end if

  end if

  return
end
subroutine hvc_check ( h, v, c )

!*****************************************************************************80
!
!! HVC_CHECK corrects out-of-range HVC color coordinates.
!
!  Discussion:
!
!    The HVC color system, developed by Tektronix, describes a color
!    based on the three qualities of hue, value, and chroma.  A given
!    color will be represented by three numbers, (H,V,C).  H, the value
!    of hue, is an angle between 0 and 360 degrees, with 0 representing
!    red.  V, the "value", is between 0 and 100.  C, the "chroma", is
!    between 0 and 100.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) H, V, C, the HVC color coordinates
!    to be checked.
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) h
  real ( kind = 8 ) v

  h = r8_modp ( h, 360.0D+00 )

  v = max ( v, 0.0D+00 )
  v = min ( v, 100.0D+00 )

  c = max ( c, 0.0D+00 )
  c = min ( c, 100.0D+00 )

  return
end
function i4_log_2 ( i )

!*****************************************************************************80
!
!! I4_LOG_2 returns the integer part of the logarithm base 2 of |I|.
!
!  Example:
!
!     I  I4_LOG_2
!
!     0  -1
!     1,  0
!     2,  1
!     3,  1
!     4,  2
!     5,  2
!     6,  2
!     7,  2
!     8,  3
!     9,  3
!    10,  3
!   127,  6
!   128,  7
!   129,  7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 2
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_2, the integer part of the logarithm
!    base 2 of the absolute value of I.
!    For positive I4_LOG_2(I), it should be true that
!      2**I4_LOG_2(X) <= |I| < 2**(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
  implicit none

  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs

  if ( i == 0 ) then

    i4_log_2 = - huge ( i4_log_2 )

  else

    i4_log_2 = 0

    i_abs = abs ( i )

    do while ( 2 <= i_abs )
      i_abs = i_abs / 2
      i4_log_2 = i4_log_2 + 1
    end do

  end if

  return
end
subroutine i4_to_angle ( i, angle )

!*****************************************************************************80
!
!! I4_TO_ANGLE maps integers to points on a circle.
!
!  Discussion:
!
!    The angles are intended to be used to select colors on a color
!    hexagon whose 6 vertices are red, yellow, green, cyan, blue,
!    magenta.
!
!  Example:
!
!     I   X      ANGLE
!
!     0   0/3      0
!     1   1/3    120
!     2   2/3    240
!
!     3   1/6     60
!     4   3/6    180
!     5   5/6    300
!
!     6   1/12    30
!     7   3/12    90
!     8   5/12   150
!     9   7/12   210
!    10   9/12   270
!    11  11/12   330
!
!    12   1/24    15
!    13   3/24    45
!    14   5/24    75
!    etc
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the desired color.
!
!    Output, real ( kind = 8 ) ANGLE, an angle, measured in degrees,
!    between 0 and 360.
!
  implicit none

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4

  if ( 0 <= abs ( i ) .and. abs ( i ) <= 2 ) then

    angle = 120.0D+00 * real ( abs ( i ), kind = 8 )

  else

    i1 = i4_log_2 ( abs ( i ) / 3 )
    i2 = abs ( i ) + 1 - 3 * 2**i1
    i3 = 2 * ( i2 - 1 ) + 1
    i4 = 3 * 2**( i1 + 1 )

    angle = 360.0D+00 * real ( i3, kind = 8 ) / real ( i4, kind = 8 )

  end if

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
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
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine interp ( ndat, x, xdat, y, ydat )

!*****************************************************************************80
!
!! INTERP does simple linear interpolation in a table.
!
!  Discussion:
!
!    The XDAT values are assumed to be sorted in ascending order.
!
!    Input values of X less than XDAT(1) are assigned a Y value of YDAT(1).
!    Values of X above XDAT(NDAT) are handled similarly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDAT, the number of data items.
!
!    Input, real ( kind = 8 ) X, the value of the independent variable.
!
!    Input, real ( kind = 8 ) XDAT(NDAT), the tabulated X values, which should
!    be in ascending order.
!
!    Output, real ( kind = 8 ) Y, the interpolated Y value.
!
!    Input, real ( kind = 8 ) YDAT(NDAT), the tabulated Y values.
!
  implicit none

  integer ( kind = 4 ) ndat

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x
  real ( kind = 8 ) xdat(ndat)
  real ( kind = 8 ) y
  real ( kind = 8 ) ydat(ndat)

  if ( x <= xdat(1) ) then

    y = ydat(1)

  else if ( xdat(ndat) <= x ) then

    y = ydat(ndat)

  else

    do i = 1, ndat-1

      if ( xdat(i) <= x .and. x <= xdat(i+1) ) then
        j = i
        exit
      end if

    end do

    y = ( ( xdat(j+1) - x           ) * ydat(j)     &
        + (             x - xdat(j) ) * ydat(j+1) ) &
        / ( xdat(j+1)     - xdat(j) )

  end if

  return
end
subroutine lab_check ( lstar, astar, bstar )

!*****************************************************************************80
!
!! LAB_CHECK corrects out-of-range CIE LAB color coordinates.
!
!  Discussion:
!
!    The CIE LAB system describes a color based on three qualities:
!      L* is CIE lightness, similar to luminance, the amount of
!        light, but adjusted to account for human perception.
!        0 <= L* <= 100.
!      a* is the amount of red chrominance, with negative values
!        indicating green.  -500 <= a* <= 500.
!      b* is the amount of yellow chrominance, with negative values
!        indicating blue.  -200 <= b* <= 200.
!    The CIE LAB model is more suitable than the CIE LUV model
!    for situations involving reflected light.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) LSTAR, ASTAR, BSTAR, the CIE LAB
!    color coordinates to be checked.
!
  implicit none

  real ( kind = 8 ) astar
  real ( kind = 8 ) bstar
  real ( kind = 8 ) lstar

  lstar = max ( lstar, 0.0D+00 )
  lstar = min ( lstar, 100.0D+00 )

  astar = max ( astar, -500.0D+00 )
  astar = min ( astar, 500.0D+00 )

  bstar = max ( bstar, -200.0D+00 )
  bstar = min ( bstar, 200.0D+00 )

  return
end
subroutine lab_prop ( lstar, astar, bstar, chroma, hue, luminance )

!*****************************************************************************80
!
!! LAB_PROP returns certain properties of a CIE LAB color.
!
!  Discussion:
!
!    The CIE LAB system describes a color based on three qualities:
!      L* is CIE lightness, similar to luminance, the amount of
!        light, but adjusted to account for human perception.
!        0 <= L* <= 100.
!      a* is the amount of red chrominance, with negative values
!        indicating green.  -500 <= a* <= 500.
!      b* is the amount of yellow chrominance, with negative values
!        indicating blue.  -200 <= b* <= 200.
!    The CIE LAB model is more suitable than the CIE LUV model
!    for situations involving reflected light.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LSTAR, ASTAR, BSTAR, the CIE LAB coordinates
!    of the color.
!
!    Output, real ( kind = 8 ) CHROMA, HUE, LUMINANCE, the chroma, hue,
!    and relative luminance.  HUE is returned as an angle, in degrees.
!
  implicit none

  real ( kind = 8 ) astar
  real ( kind = 8 ) bstar
  real ( kind = 8 ) chroma
  real ( kind = 8 ) hue
  real ( kind = 8 ) lstar
  real ( kind = 8 ) luminance
  real ( kind = 8 ) radians_to_degrees

  chroma = sqrt ( astar**2 + bstar**2 )

  hue = atan2 ( bstar, astar )
  hue = radians_to_degrees ( hue )

  if ( lstar <= 0.0D+00 ) then
    luminance = 0.0D+00
  else if ( lstar <= 903.3D+00 * 0.008856D+00 ) then
    luminance = lstar / 903.3D+00
  else if ( lstar <= 100.0D+00 ) then
    luminance = ( ( lstar + 16.0D+00 ) / 116.0D+00 )**3
  else
    luminance = 1.0D+00
  end if

  return
end
subroutine lab_to_xyz_cap ( lstar, astar, bstar, xcap, ycap, zcap, xcapn, &
  ycapn, zcapn )

!*****************************************************************************80
!
!! LAB_TO_XYZ_CAP converts CIE LAB to CIE XYZ color coordinates.
!
!  Discussion:
!
!    The CIE LAB system describes a color based on three qualities:
!      L* is CIE lightness, similar to luminance, the amount of
!        light, but adjusted to account for human perception.
!        0 <= L* <= 100.
!      a* is the amount of red chrominance, with negative values
!        indicating green.  -500 <= a* <= 500.
!      b* is the amount of yellow chrominance, with negative values
!        indicating blue.  -200 <= b* <= 200.
!    The CIE LAB model is more suitable than the CIE LUV model
!    for situations involving reflected light.
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Giorgianni and Thomas Madden,
!    Digital Color Management, Encoding Solutions,
!    Addison Wesley, 1998.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LSTAR, ASTAR, BSTAR, the CIE LAB color
!    coordinates.
!
!    Output, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!
!    Input, real ( kind = 8 ) XCAPN, YCAPN, ZCAPN, the CIE XYZ color
!    coordinates of white.
!
  implicit none

  real ( kind = 8 ) astar
  real ( kind = 8 ) bstar
  real ( kind = 8 ) fx
  real ( kind = 8 ) fy
  real ( kind = 8 ) fz
  real ( kind = 8 ) lstar
  real ( kind = 8 ) r8_cubert
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcapn
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycapn
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcapn

  if ( lstar <= 0.0D+00 ) then
    ycap = 0.0D+00
  else if ( lstar <= 8.0D+00 ) then
    ycap = lstar * ycapn / 903.3D+00
  else if ( lstar <= 100.0D+00 ) then
    ycap = ycapn * ( ( lstar + 16.0D+00 ) / 116.0D+00 )**3
  else
    ycap = ycapn
  end if

  if ( ycap <= 0.00856D+00 * ycapn ) then
    fy = 7.787D+00 * ycap / ycapn + 16.0D+00 / 116.0D+00
  else
    fy = r8_cubert ( ycap / ycapn )
  end if

  fx = fy + ( astar / 500.0D+00 )

  if ( fx**3 <= 0.008856D+00 ) then
    xcap = xcapn * ( fx - 16.0D+00 / 116.0D+00 ) / 7.787D+00
  else
    xcap = xcapn * fx**3
  end if

  fz = fy - ( bstar / 200.0D+00 )

  if ( fz**3 <= 0.008856D+00 ) then
    zcap = zcapn * ( fz - 16.0D+00 / 116.0D+00 ) / 7.787D+00
  else
    zcap = zcapn * fz**3
  end if

  return
end
subroutine lcc_to_rgbprime ( luma, chroma1, chroma2, yr, yg, yb, rprime, &
  gprime, bprime )

!*****************************************************************************80
!
!! LCC_TO_RGBPRIME converts LCC to R'G'B' color coordinates.
!
!  Discussion:
!
!    The LCC color coordinate system records a color as three components,
!    Luma, Chroma1 and Chroma2.  The LCC color coordinates are used
!    in an intermediate calculation of the PhotoYCC color coordinates.
!    Luma is scaled luminance, Chroma1 is (B'-Luma) and Chroma2 is (R'-Luma).
!
!    The R'G'B' color system is a nonlinear video signal measurement.
!    Each coordinate must be between 0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LUMA, CHROMA1, CHROMA2, the LCC color coordinates.
!
!    Input, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
!    Output, real ( kind = 8 ) RPRIME, GPRIME, BPRIME, the R'G'B' color
!    coordinates.
!
  implicit none

  real ( kind = 8 ) chroma1
  real ( kind = 8 ) chroma2
  real ( kind = 8 ) bprime
  real ( kind = 8 ) gprime
  real ( kind = 8 ) luma
  real ( kind = 8 ) rprime
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr

  yprime = luma

  rprime = yprime + chroma2
  bprime = yprime + chroma1

  gprime = ( luma - yr * rprime - yb * bprime ) / yg

  return
end
subroutine lcc_to_ycbcr ( luma, chroma1, chroma2, yprime, cb, cr )

!*****************************************************************************80
!
!! LCC_TO_YCBCR converts LCC to Y'CbCr color coordinates.
!
!  Discussion:
!
!    The LCC color coordinate system records a color as three components,
!    Luma, Chroma1 and Chroma2.  The LCC color coordinates are used
!    in an intermediate calculation of the PhotoYCC color coordinates.
!    Luma is scaled luminance, Chroma1 is (B'-Luma) and Chroma2 is (R'-Luma).
!
!    The Y'CbCr color system is used in digital television signals.
!    The Y' component measures luma, an approximation to the luminance
!    or amount of light.  Y' is the only component displayed on black
!    and white televisions.  The Cb and Cr components contain measures
!    of the blue and red components of the color.  Y' should be between
!    0 and 1, with reference black at 16/255 and reference white at 235/255,
!    while Cb and Cr should be between -0.5 and 0.5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    C Wayne Brown and Barry Shepherd,
!    Graphics File Formats,
!    Manning Publications, 1995.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LUMA, CHROMA1, CHROMA2, the LCC color coordinates.
!
!    Output, real ( kind = 8 ) YPRIME, CB, CR, the Y'CbCr color coordinates.
!
  implicit none

  real ( kind = 8 ) cb
  real ( kind = 8 ) chroma1
  real ( kind = 8 ) chroma2
  real ( kind = 8 ) cr
  real ( kind = 8 ) luma
  real ( kind = 8 ) yprime

  yprime = ( 219.0D+00 * luma + 16.0D+00 ) / 255.0D+00
  cb = ( 224.0D+00 * 0.564D+00 * chroma1 + 128.0D+00 ) / 255.0D+00
  cr = ( 224.0D+00 * 0.713D+00 * chroma2 + 128.0D+00 ) / 255.0D+00

  return
end
subroutine lcc_to_ycc ( luma, chroma1, chroma2, yprime, c1, c2 )

!*****************************************************************************80
!
!! LCC_TO_YCC converts LCC to PhotoYCC color coordinates.
!
!  Discussion:
!
!    The LCC color coordinate system records a color as three components,
!    Luma, Chroma1 and Chroma2.  The LCC color coordinates are used
!    in an intermediate calculation of the PhotoYCC color coordinates.
!    Luma is scaled luminance, Chroma1 is (B'-Luma) and Chroma2 is (R'-Luma).
!
!    The Kodak PhotoYCC Color Interchange Space was developed for the
!    Photo CD System.  The Y' coordinate is a measure of luminance,
!    while C1 and C2 measure color difference chrominance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Giorgianni and Thomas Madden,
!    Digital Color Management, Encoding Solutions,
!    Addison Wesley, 1998.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LUMA, CHROMA1, CHROMA2, the LCC color coordinates.
!
!    Output, real ( kind = 8 ) YPRIME, C1, C2, the PhotoYCC color coordinates.
!
  implicit none

  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) chroma1
  real ( kind = 8 ) chroma2
  real ( kind = 8 ) luma
  real ( kind = 8 ) yprime

  yprime = 255.0D+00 * luma / 1.402D+00
  c1 = 111.40D+00 * chroma1 + 156.0D+00
  c2 = 135.64D+00 * chroma2 + 137.0D+00

  return
end
subroutine lin_to_nonlin ( r, rprime )

!*****************************************************************************80
!
!! LIN_TO_NONLIN converts a linear light intensity to nonlinear video signal.
!
!  Discussion:
!
!    This process is also known as gamma correction.
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
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the linear light intensity.
!
!    Output, real ( kind = 8 ) RPRIME, the nonlinear video signal.
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) rprime

  if ( r < -0.018D+00 ) then
    rprime = 0.099D+00 - 1.099D+00 * ( abs ( r ) )**0.45D+00
  else if ( abs ( r ) <= 0.018D+00 ) then
    rprime = 4.5D+00 * r
  else
    rprime = -0.099D+00 + 1.099D+00 * r**0.45D+00
  end if

  return
end
subroutine luv_check ( lstar, ustar, vstar )

!*****************************************************************************80
!
!! LUV_CHECK corrects out-of-range CIE LUV color coordinates.
!
!  Discussion:
!
!    The CIE LUV system describes a color based on three qualities:
!      L* is CIE lightness, similar to luminance, the amount of light,
!        but adjusted for human perception.  0 <= L* <= 100.
!      u* is the amount of red chrominance, with negative values
!        indicating green.
!      v* is the amount of yellow chrominance, with negative values
!        indicating blue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) LSTAR, USTAR, VSTAR, the CIE LUV
!    color coordinates to be checked.
!
  implicit none

  real ( kind = 8 ) lstar
  real ( kind = 8 ) ustar
  real ( kind = 8 ) vstar

  lstar = max ( lstar, 0.0D+00 )
  lstar = min ( lstar, 100.0D+00 )

  return
end
subroutine luv_prop ( lstar, ustar, vstar, chroma, hue, luminance, sat )

!*****************************************************************************80
!
!! LUV_PROP returns certain properties of a CIE LUV color.
!
!  Discussion:
!
!    The CIE LUV system describes a color based on three qualities:
!      L* is CIE lightness, similar to luminance, the amount of light,
!        but adjusted for human perception.  0 <= L* <= 100.
!      u* is the amount of red chrominance, with negative values
!        indicating green.
!      v* is the amount of yellow chrominance, with negative values
!        indicating blue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LSTAR, USTAR, VSTAR, the CIE LUV coordinates
!    of the color.
!
!    Output, real ( kind = 8 ) CHROMA, HUE, LUMINANCE, SAT, the chroma,
!    hue, relative luminance, and saturation.  HUE is returned as an
!    angle in degrees.
!
  implicit none

  real ( kind = 8 ) chroma
  real ( kind = 8 ) hue
  real ( kind = 8 ) lstar
  real ( kind = 8 ) luminance
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ) sat
  real ( kind = 8 ) ustar
  real ( kind = 8 ) vstar

  chroma = sqrt ( ustar**2 + vstar**2 )

  hue = atan2 ( vstar, ustar )
  hue = radians_to_degrees ( hue )

  if ( lstar == 0.0D+00 ) then
    sat = 0.0D+00
  else
    sat = chroma / lstar
  end if

  if ( lstar <= 0.0D+00 ) then
    luminance = 0.0D+00
  else if ( lstar <= 903.3D+00 * 0.008856D+00 ) then
    luminance = lstar / 903.3D+00
  else if ( lstar <= 100.0D+00 ) then
    luminance = ( ( lstar + 16.0D+00 ) / 116.0D+00 )**3
  else
    luminance = 1.0D+00
  end if

  return
end
subroutine luv_to_xyz_cap ( lstar, ustar, vstar, xcap, ycap, zcap, xcapn, &
  ycapn, zcapn )

!*****************************************************************************80
!
!! LUV_TO_XYZ_CAP converts CIE LUV to CIE XYZ color coordinates.
!
!  Discussion:
!
!    The CIE LUV system describes a color based on three qualities:
!      L* is CIE lightness, similar to luminance, the amount of light,
!        but adjusted for human perception.  0 <= L* <= 100.
!      u* is the amount of red chrominance, with negative values
!        indicating green.
!      v* is the amount of yellow chrominance, with negative values
!        indicating blue.
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LSTAR, USTAR, VSTAR, the CIE LUV color
!    coordinates.
!
!    Output, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!
!    Input, real ( kind = 8 ) XCAPN, YCAPN, ZCAPN, the CIE XYZ color
!    coordinates of white.
!
  implicit none

  real ( kind = 8 ) lstar
  real ( kind = 8 ) ustar
  real ( kind = 8 ) uprime
  real ( kind = 8 ) unprime
  real ( kind = 8 ) vstar
  real ( kind = 8 ) vprime
  real ( kind = 8 ) vnprime
  real ( kind = 8 ) wnprime
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcapn
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycapn
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcapn

  if ( lstar <= 0.0D+00 ) then

    xcap = 0.0D+00
    ycap = 0.0D+00
    zcap = 0.0D+00

  else
!
!  Compute CIE luminance from L* and YCAPN.
!
    if ( lstar <= 0.0D+00 ) then
      ycap = 0.0D+00
    else if ( lstar <= 903.3D+00 * 0.008856D+00 ) then
      ycap = lstar * ycapn / 903.3D+00
    else if ( lstar <= 100.0D+00 ) then
      ycap = ycapn * ( ( lstar + 16.0D+00 ) / 116.0D+00 )**3
    else
      ycap = ycapn
    end if
!
!  Compute (un',vn') from (XCAPN,YCAPN,ZCAPN).
!
    call xyz_cap_to_uvwprime ( xcapn, ycapn, zcapn, unprime, vnprime, wnprime )
!
!  Compute (u',v') from (un',vn') and (l*,u*,v*).
!
    uprime = unprime + ustar / ( 13.0D+00 * lstar )
    vprime = vnprime + vstar / ( 13.0D+00 * lstar )
!
!  Now compute XCAP and ZCAP from (u',v') and YCAP.
!
    call uvprimey_to_xyz_cap ( uprime, vprime, xcap, ycap, zcap )

  end if

  return
end
subroutine name_test ( itest, name )

!*****************************************************************************80
!
!! NAME_TEST supplies color names for tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ITEST, the index of the test.
!
!    Output, character ( len = * ) NAME, the name of a color to be tested.
!    The longest name is 11 characters.  If ITEST is out of the range
!    of data, then NAME is returned as ' '.
!
  implicit none

  integer ( kind = 4 ) itest
  character ( len = * ) name

  if ( itest == 1 ) then
    name = 'Red'
  else if ( itest == 2 ) then
    name = 'Green'
  else if ( itest == 3 ) then
    name = 'Blue'
  else if ( itest == 4 ) then
    name = 'Cyan'
  else if ( itest == 5 ) then
    name = 'Magenta'
  else if ( itest == 6 ) then
    name = 'Yellow'
  else if ( itest == 7 ) then
    name = 'White'
  else if ( itest == 8 ) then
    name = 'Black'
  else if ( itest == 9 ) then
    name = 'Pink'
  else if ( itest == 10 ) then
    name = 'Aquamarine'
  else if ( itest == 11 ) then
    name = 'Tan'
  else if ( itest == 12 ) then
    name = 'YellowGreen'
  else if ( itest == 13 ) then
    name = 'Maroon'
  else if ( itest == 14 ) then
    name = 'Salmon'
  else if ( itest == 15 ) then
    name = 'Mauve'
  else
    name = ' '
  end if

  return
end
subroutine name_to_primaries ( name, rx, ry, gx, gy, bx, by, wx, wy )

!*****************************************************************************80
!
!! NAME_TO_PRIMARIES returns CIE xy chromaticities of television primaries.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Martindale and Alan Paeth,
!    Television Color Encoding and Hot Broadcast Colors,
!    Graphics Gems II, edited by James Arvo,
!    Academic Press, 1991, pages 147-158.
!
!  Parameters:
!
!    Input, character ( len = * ) NAME, the name of a television system:
!    'CIE', the CIE primaries;
!    'EBU',
!    'HDTV',
!    'NTSC', National Television Systems Committee;
!    'SMPTE', Society of Motion Picture and Television Engineers.
!
!    Output, real ( kind = 8 ) RX, RY, the xy chromaticities of the R primary.
!
!    Output, real ( kind = 8 ) GX, GY, the xy chromaticities of the G primary.
!
!    Output, real ( kind = 8 ) BX, BY, the xy chromaticities of the B primary.
!
!    Output, real ( kind = 8 ) WX, WY, the xy chromaticities of the
!    reference white.
!
  implicit none

  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) gx
  real ( kind = 8 ) gy
  character ( len = * ) name
  character ( len = 30 ) name_copy
  real ( kind = 8 ) rx
  real ( kind = 8 ) ry
  logical s_eqi
  real ( kind = 8 ) wx
  real ( kind = 8 ) wy
!
!  Make a temporary copy of NAME.
!
  name_copy = adjustl ( name )
!
!  Remove blanks and underlines.
!
  call s_c_delete ( name_copy, ' ' )

  call s_c_delete ( name_copy, '_' )

  if ( s_eqi ( name_copy, 'CIE' ) ) then

    rx = 0.73467D+00
    ry = 0.26533D+00
    gx = 0.27376D+00
    gy = 0.71741D+00
    bx = 0.16658D+00
    by = 0.00886D+00
    wx = 1.0D+00 / 3.0D+00
    wy = 1.0D+00 / 3.0D+00

  else if ( s_eqi ( name_copy, 'EBU' ) ) then

    rx = 0.64D+00
    ry = 0.33D+00
    gx = 0.29D+00
    gy = 0.60D+00
    bx = 0.15D+00
    by = 0.06D+00
    wx = 0.3127D+00
    wy = 0.3291D+00

  else if ( s_eqi ( name_copy, 'HDTV' ) ) then

    rx = 0.670D+00
    ry = 0.330D+00
    gx = 0.210D+00
    gy = 0.710D+00
    bx = 0.150D+00
    by = 0.060D+00
    wx = 0.3127D+00
    wy = 0.3291D+00

  else if ( s_eqi ( name_copy, 'NTSC' ) ) then

    rx = 0.67D+00
    ry = 0.33D+00
    gx = 0.21D+00
    gy = 0.71D+00
    bx = 0.14D+00
    by = 0.08D+00
    wx = 0.3101D+00
    wy = 0.3162D+00

  else if ( s_eqi ( name_copy, 'SMPTE' ) ) then

    rx = 0.630D+00
    ry = 0.340D+00
    gx = 0.310D+00
    gy = 0.595D+00
    bx = 0.155D+00
    by = 0.070D+00
    wx = 0.3127D+00
    wy = 0.3291D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NAME_TO_PRIMARIES - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized name: "' // trim ( name ) // '"'
    stop

  end if

  return
end
subroutine name_to_rgb ( name, r, g, b )

!*****************************************************************************80
!
!! NAME_TO_RGB converts a string to RGB colors.
!
!  Discussion:
!
!    The names and information are read from the file "COLORS.TXT", a
!    modified version of the X Windows color data file "RGB.TXT".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NAME, the name of a color.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB coordinates.
!    However, these will be returned as ( -1.0, -1.0, -1.0 ) if the
!    color name was not recognized.
!
  implicit none

  real ( kind = 8 ) b
  character ( len = 20 ) :: color_file = 'colors.txt'
  real ( kind = 8 ) g
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ig
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) iunit
  character ( len = * ) name
  character ( len = 30 ) name_copy
  character ( len = 30 ) name_color
  real ( kind = 8 ) r
  logical s_eqi
!
!  Make a temporary copy of the NAME.
!  Remove blanks and underlines.
!
  name_copy = adjustl ( name )

  call s_c_delete ( name_copy, ' ' )

  call s_c_delete ( name_copy, '_' )

  call get_unit ( iunit )

  open ( unit = iunit, file = color_file, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NAME_TO_RGB - Fatal error!'
    write ( *, '(a)' ) '  Could not open the color name file:'
    write ( *, '(a)' ) trim ( color_file )
    stop
  end if

  do

    read ( iunit, '(i3,2x,i3,2x,i3,2x,a)', iostat = ios ) ir, ig, ib, name_color

    if ( ios /= 0 ) then
      r = - 1.0D+00
      g = - 1.0D+00
      b = - 1.0D+00
      exit
    end if

    name_color = adjustl ( name_color )

    call s_c_delete ( name_color, ' ' )

    call s_c_delete ( name_color, '_' )

    if ( s_eqi ( name_copy, name_color ) ) then
      r = real ( ir, kind = 8 ) / 255.0D+00
      g = real ( ig, kind = 8 ) / 255.0D+00
      b = real ( ib, kind = 8 ) / 255.0D+00
      exit
    end if

  end do

  close ( unit = iunit )

  return
end
subroutine name_to_xyz ( name, x, y, z )

!*****************************************************************************80
!
!! NAME_TO_XYZ converts a color or illuminant name to CIE xyz color coordinates.
!
!  Discussion:
!
!    Thanks to Harald Anlauf of the Technical University of Darmstadt,
!    for pointing out a programming error which meant that NAME was not
!    an input-only variable.  (30 April 2002)
!
!    In the CIE color system, the exact chromaticities of several
!    standard illuminants were defined.  These were generally chosen to
!    correspond to common lighting situations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NAME, the name of a color.  Before considering
!    the name, the routine removes all blanks and underscores, and
!    capitalizes the name.  Legal names include:
!
!      'A', the CIE illuminant A, light from a tungsten lamp, at 500 watts;
!      'B', the CIE illuminant B, direct sunlight, at 500 watts;
!      'C', the CIE illuminant C, average sunlight, at 500 watts;
!        This is used as the reference white for NTSC color encoding.
!      'D50' or 'D5000', the CIE illuminant used in graphics printing,
!        bright tunsten illumination;
!      'D55' or 'D5500', a CIE illuminant that approximates a cloudy bright day;
!      'D65' or 'D6500', the CIE illuminant that approximates daylight;
!        This is used as the reference white for SMPTE, PAL/EBU, and HDTV
!        color encoding.
!      'E', the CIE illuminant E, normalized reference source.
!
!    Output, real ( kind = 8 ) X, Y, Z, the corresponding CIE xyz color
!    coordinates, or (0,0,0) if the name is not recognized.
!
  implicit none

  character ( len = * ) name
  character ( len = 20 ) name_copy
  logical s_eqi
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!
!  Make a temporary copy of NAME.
!
  name_copy = adjustl ( name )
!
!  Remove blanks and underlines.
!
  call s_c_delete ( name_copy, ' ' )

  call s_c_delete ( name_copy, '_' )
!
!  Compare the input name to the recognized list.
!
  if ( s_eqi ( name_copy, 'A' ) ) then

    x = 0.448D+00
    y = 0.407D+00
    z = 0.145D+00

  else if ( s_eqi ( name_copy, 'B' ) ) then

    x = 0.349D+00
    y = 0.352D+00
    z = 0.299D+00

  else if ( s_eqi ( name_copy, 'C' ) ) then

    x = 0.3101D+00
    y = 0.3162D+00
    z = 0.3737D+00

  else if ( s_eqi ( name_copy, 'D50' ) .or. &
            s_eqi ( name_copy, 'D5000' ) ) then

    x =  96.42D+00 / ( 96.42D+00 + 100.00D+00 + 82.49D+00 )
    y = 100.00D+00 / ( 96.42D+00 + 100.00D+00 + 82.49D+00 )
    z =  82.49D+00 / ( 95.42D+00 + 100.00D+00 + 82.49D+00 )

  else if ( s_eqi ( name_copy, 'D55' ) .or. &
            s_eqi ( name_copy, 'D5500' ) ) then

    x = 0.3324D+00
    y = 0.3474D+00
    z = 0.3202D+00

  else if ( s_eqi ( name_copy, 'D65' ) .or. &
            s_eqi ( name_copy, 'D6500' ) ) then

    x =  95.05D+00 / ( 95.05D+00 + 100.00D+00 + 108.91D+00 )
    y = 100.00D+00 / ( 95.05D+00 + 100.00D+00 + 108.91D+00 )
    z = 108.91D+00 / ( 95.05D+00 + 100.00D+00 + 108.91D+00 )

  else if ( s_eqi ( name_copy, 'E' ) ) then

    x = 100.0D+00 / 300.0D+00
    y = 100.0D+00 / 300.0D+00
    z = 100.0D+00 / 300.0D+00

  else

    x = 0.0D+00
    y = 0.0D+00
    z = 0.0D+00

  end if

  return
end
subroutine ncs_check ( c1, c2, n, c, s )

!*****************************************************************************80
!
!! NCS_CHECK corrects out-of-range NCS color coordinates.
!
!  Discussion:
!
!    The NCS or "natural color system" describes a color based on:
!    * C1 and C2, two elementary colors from the sequence RYGB or
!      C2 = blank for a pure elementary color, or
!      C1 = N, C2 = blank for a neutral color);
!    * N, the percentage of C2;
!    * C, the colorfulness or strength, as a percentage;
!    * S, the blackness as a percentage.
!
!    The scant documentation I have seen claims that the percentages are
!    always less than 100.  I don't see why, and for now I'll let them
!    lie between 0 and 100.  The NCS designation for a color has the form
!    "CCSS C1NC2".
!
!  Example:
!
!              C1 C2  N   C   S    Designation
!              -- -- -- --- ---    -----------
!    Black   (  N  *  0   0 100 )  0099 N
!    Blue    (  B  *  0 100   0 )  9900 B
!    Cyan    (  G  B 50 100   0 )  9900 G50B
!    Green   (  G  *  0 100   0 )  9900 G
!    Magenta (  B  R 50 100   0 )  9900 B50R
!    Orange  (  R  Y 50 100   0 )  9900 R50Y
!    Red     (  R  *  0 100   0 )  9900 R
!    White   (  N  *  0   0   0 )  0000 N
!    Yellow  (  Y  *  0 100   0 )  9900 Y
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Olof Kylander and Karin Kylander,
!    GIMP: The Official Manual,
!    Coriolis Open Press, 1999.
!
!  Parameters:
!
!    Input/output, character C1, C2, integer N, C, S, the NCS color
!    coordinates to be checked.  If the color information is very
!    bad, it is replaced by the designation for black.
!
  implicit none

  integer ( kind = 4 ) c
  character c1
  character c2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
!
!  Replace C1 by its capitalized value.
!
  call ch_cap ( c1 )
!
!  We're expecting only the values R, Y, G, B, and possibly N (for neutral).
!
  if ( c1 == 'R' .or. c1 == 'Y' .or. c1 == 'G' .or. c1 == 'B' ) then

  else if ( c1 == 'N' .and. c2 == ' ' .and. c == 0 ) then

  else

    c1 = 'N'
    c2 = ' '
    n = 0
    c = 0
    s = 100
    return

  end if
!
!  Replace C2 by its capitalized value.
!  C2 may be blank, but then N should be 0.
!
  call ch_cap ( c2 )
!
!  We're expecting only the values R, Y, G, B and blank.
!
  if ( c2 == 'R' .or. c2 == 'Y' .or. c2 == 'G' .or. c2 == 'B' ) then

  else if ( c2 == ' ' .and. n == 0 ) then

  else

    c1 = 'N'
    c2 = ' '
    n = 0
    c = 0
    s = 100
    return

  end if
!
!  If necessary, swap the colors so they have a preferred order.
!
       if ( c1 == 'Y' .and. c2 == 'R' ) then
    c1 = 'R'
    c2 = 'Y'
    n = 100 - n
  else if ( c1 == 'G' .and. c2 == 'Y' ) then
    c1 = 'Y'
    c2 = 'G'
    n = 100 - n
  else if ( c1 == 'B' .and. c2 == 'G' ) then
    c1 = 'G'
    c2 = 'B'
    n = 100 - n
  else if ( c1 == 'R' .and. c2 == 'B' ) then
    c1 = 'B'
    c2 = 'R'
    n = 100 - n
  end if
!
!  Only certain pairs of colors are allowed.
!
       if ( c2 == 'R' .and. c1 == 'Y' ) then
  else if ( c2 == 'Y' .and. c1 == 'G' ) then
  else if ( c2 == 'G' .and. c1 == 'B' ) then
  else if ( c2 == 'B' .and. c1 == 'R' ) then
  else if ( c2 == ' ' .and. c1 == 'N' ) then
  else
    c1 = 'N'
    c2 = ' '
    n = 0
    c = 0
    s = 100
  end if
!
!  Only certain values of N are allowed.
!
  if ( n == 0 ) then
    c2 = ' '
  else if ( n == 100 ) then
    c1 = c2
    c2 = ' '
  else if ( n < 0 ) then
    n = 0
  else if ( 100 < n ) then
    n = 100
  else

  end if
!
!  Only certain values of C are allowed.
!  Here, we will "repair" such values.
!
  if ( c < 0 ) then
    c = 0
  else if ( 100 < c ) then
    c = 100
  else

  end if
!
!  Only certain values of S are allowed.
!  Here, we will "repair" such values.
!
  if ( s < 0 ) then
    s = 0
  else if ( 100 < s ) then
    s = 100
  else

  end if
!
!  C + S must be no more than 100.
!
  if ( 100 < c + s ) then
    c1 = 'N'
    c2 = ' '
    n = 0
    c = 0
    s = 100
    return
  end if

  return
end
subroutine ncs_to_rgb ( c1, c2, n, c, s, r, g, b )

!*****************************************************************************80
!
!! NCS_TO_RGB converts NCS to RGB color coordinates.
!
!  Discussion:
!
!    It has been difficult to find two descriptions of NCS that agree.
!    It has been impossible to find a description of the mathematical
!    details of NCS.
!    On top of that, I seem to have mixed up a few things, including
!    the ordering of the colors.  For now, I will be content that
!    this routine embodies the "flavor" of NCS, and does a reasonable
!    job of inverting RGB_TO_NCS, modulo the loss of information
!    because of real to integer percentage truncation.
!    I'll come back and fix this one day.
!
!
!    The NCS or "natural color system" describes a color based on:
!    * C1 and C2, two elementary colors from the sequence RYGB or
!      C2 = blank for a pure elementary color, or
!      C1 = N, C2 = blank for a neutral color);
!    * N, the percentage of C2;
!    * C, the colorfulness or strength, as a percentage;
!    * S, the blackness as a percentage.
!
!    The scant documentation I have seen claims that the percentages are
!    always less than 100.  I don't see why, and for now I'll let them
!    lie between 0 and 100.  The NCS designation for a color has the form
!    "CCSS C1NC2".
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Olof Kylander and Karin Kylander,
!    GIMP: The Official Manual,
!    Coriolis Open Press, 1999.
!
!  Parameters:
!
!    Input, character C1, C2, integer N, C, S, the NCS color coordinates.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  character c1
  character c2
  real ( kind = 8 ) color(3)
  real ( kind = 8 ) color1(3)
  real ( kind = 8 ) color2(3)
  real ( kind = 8 ) color3(3)
  real ( kind = 8 ) g
  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) w
!
!  Determine the colors that bracket the given color.
!
       if ( c1 == 'R' ) then
    color1(1:3) = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)
    color2(1:3) = (/ 1.0D+00, 1.0D+00, 0.0D+00 /)
  else if ( c1 == 'Y' ) then
    color1(1:3) = (/ 1.0D+00, 1.0D+00, 0.0D+00 /)
    color2(1:3) = (/ 0.0D+00, 1.0D+00, 0.0D+00 /)
  else if ( c1 == 'G' ) then
    color1(1:3) = (/ 0.0D+00, 1.0D+00, 0.0D+00 /)
    color2(1:3) = (/ 0.0D+00, 0.0D+00, 1.0D+00 /)
  else if ( c1 == 'B' ) then
    color1(1:3) = (/ 0.0D+00, 0.0D+00, 1.0D+00 /)
    color2(1:3) = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)
  end if
!
!  Apply the value of N.
!  This is trickier than it looks!
!
       if ( c1 == 'R' ) then

    if ( n < 50 ) then

      color3(1:3) = ( real ( 100 - n, kind = 8 ) * color1(1:3)   &
                    + real (       n, kind = 8 ) * color2(1:3) ) &
                    / real ( 100 - n, kind = 8 )
    else
      color3(1:3) = ( real ( 100 - n, kind = 8 ) * color1(1:3)   &
                    + real (       n, kind = 8 ) * color2(1:3) ) &
                    / real (       n, kind = 8 )
    end if

  else if ( c1 == 'Y' ) then

    if ( n < 50 ) then

      color3(1:3) = ( real ( 100 - n, kind = 8 ) * color1(1:3)   &
                    + real (       n, kind = 8 ) * color2(1:3) ) &
                    / real ( 100 - n, kind = 8 )
    else
      color3(1:3) = ( real ( 100 - n, kind = 8 ) * color1(1:3)   &
                    + real (       n, kind = 8 ) * color2(1:3) ) &
                    / real (       n, kind = 8 )
    end if

  else if ( c1 == 'G' ) then

    if ( n < 50 ) then

      color3(1:3) = ( real ( 100 - n, kind = 8 ) * color1(1:3)   &
                    + real (       n, kind = 8 ) * color2(1:3) ) &
                    / real ( 100 - n, kind = 8 )
    else
      color3(1:3) = ( real ( 100 - n, kind = 8 ) * color1(1:3)   &
                    + real (       n, kind = 8 ) * color2(1:3) ) &
                    / real (       n, kind = 8 )
    end if

  else if ( c1 == 'B' ) then

    if ( n < 50 ) then

      color3(1:3) = ( real ( 100 - n, kind = 8 ) * color1(1:3)   &
                    + real (       n, kind = 8 ) * color2(1:3) ) &
                    / real ( 100 - n, kind = 8 )
    else
      color3(1:3) = ( real ( 100 - n, kind = 8 ) * color1(1:3)   &
                    + real (       n, kind = 8 ) * color2(1:3) ) &
                    / real (       n, kind = 8 )
    end if

  end if
!
!  The color can now be computed as the triangular sum of white,
!  black, and the color from the color circle.
!
  w = 100 - c - s

  color(1:3) = ( &
      real ( c, kind = 8 ) * color3(1:3) &
    + real ( s, kind = 8 ) * (/ 0.0D+00, 0.0D+00, 0.0D+00 /) &
    + real ( w, kind = 8 ) * (/ 1.0D+00, 1.0D+00, 1.0D+00 /) ) / 100.0D+00

  r = color(1)
  g = color(2)
  b = color(3)

  return
end
subroutine nm_to_rgbcie ( w, r, g, b )

!*****************************************************************************80
!
!! NM_TO_RGBCIE converts a pure light wavelength to CIE RGB color coordinates.
!
!  Discussion:
!
!    The "RGB" used here is not the "RGB" referred to in RGB monitors.  Instead,
!    it refers to three specific monochromatic light sources:
!
!      R: red light at 700.0 nanometers;
!      G: green light at 546.1 nanometers;
!      B: blue light at 435.8 nanometers.
!
!    The relative intensities of these lights are adjusted so that
!
!      1 R + 1 G + 1 B = white.
!
!    For a light emission with given spectral power distribution SPD(W),
!    the values tabulated here can be used to produce integrated RGB
!    tristimulus values:
!
!      rval = integral ( 380 <= w <= 760 ) spd(w) * rbar(w) d w;
!
!      gval = integral ( 380 <= w <= 760 ) spd(w) * gbar(w) d w;
!
!      bval = integral ( 380 <= w <= 760 ) spd(w) * bbar(w) d w.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W, the wavelength of the monochromatic light whose
!    RGB color coordinatesare desired.  W must lie between
!    380 and 760.
!
!    Output, real ( kind = 8 ) R, G, B, the CIE RGB color coordinates for
!    the monochromatic light of the given wavelength.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndat = 39

  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( ndat ) :: bdat = (/ &
     0.00117D+00,  0.00359D+00,  0.01214D+00,  0.03707D+00,  0.11541D+00, &
     0.24769D+00,  0.31228D+00,  0.31670D+00,  0.39821D+00,  0.22991D+00, &
     0.14494D+00,  0.08257D+00,  0.04776D+00,  0.02698D+00,  0.01221D+00, &
     0.00549D+00,  0.00146D+00, -0.00058D+00,  0.00130D+00, -0.00135D+00, &
    -0.00108D+00, -0.00079D+00, -0.00049D+00, -0.00030D+00, -0.00015D+00, &
    -0.00008D+00, -0.00003D+00, -0.00001D+00,  0.00000D+00,  0.00000D+00, &
     0.00000D+00,  0.00000D+00,  0.00000D+00,  0.00000D+00,  0.00000D+00, &
     0.00000D+00,  0.00000D+00,  0.00000D+00,  0.00000D+00 /)
  real ( kind = 8 ) g
  real ( kind = 8 ), save, dimension ( ndat ) :: gdat = (/ &
    -0.00001D+00, -0.00004D+00, -0.00014D+00, -0.00041D+00, -0.00110D+00, &
    -0.00119D+00,  0.00149D+00,  0.00678D+00,  0.01485D+00,  0.02538D+00, &
     0.03914D+00,  0.05689D+00,  0.08536D+00,  0.12860D+00,  0.17468D+00, &
     0.20317D+00,  0.21466D+00,  0.21178D+00,  0.19702D+00,  0.17087D+00, &
     0.13610D+00,  0.09754D+00,  0.06246D+00,  0.03557D+00,  0.01828D+00, &
     0.00833D+00,  0.00334D+00,  0.00116D+00,  0.00037D+00,  0.00011D+00, &
     0.00003D+00,  0.00000D+00,  0.00000D+00,  0.00000D+00,  0.00000D+00, &
     0.00000D+00,  0.00000D+00,  0.00000D+00,  0.00000D+00 /)
  real ( kind = 8 ) r
  real ( kind = 8 ), save, dimension ( ndat ) :: rdat = (/ &
     0.00003D+00,  0.00010D+00,  0.00030D+00,  0.00084D+00,  0.00211D+00, &
     0.00218D+00, -0.00261D+00, -0.01213D+00, -0.02608D+00, -0.03933D+00, &
    -0.04939D+00, -0.05814D+00, -0.07173D+00, -0.08901D+00, -0.09264D+00, &
    -0.07101D+00, -0.03152D+00,  0.02279D+00,  0.09060D+00,  0.16768D+00, &
     0.24526D+00,  0.30928D+00,  0.34429D+00,  0.33971D+00,  0.29708D+00, &
     0.22677D+00,  0.15968D+00,  0.10167D+00,  0.05932D+00,  0.03149D+00, &
     0.01687D+00,  0.00819D+00,  0.00410D+00,  0.00210D+00,  0.00105D+00, &
     0.00052D+00,  0.00025D+00,  0.00012D+00,  0.00006D+00 /)
  real ( kind = 8 ) w
  real ( kind = 8 ), save, dimension ( ndat ) :: wdat = (/ &
    380.0D+00, 390.0D+00, 400.0D+00, 410.0D+00, 420.0D+00, &
    430.0D+00, 440.0D+00, 450.0D+00, 460.0D+00, 470.0D+00, &
    480.0D+00, 490.0D+00, 500.0D+00, 510.0D+00, 520.0D+00, &
    530.0D+00, 540.0D+00, 550.0D+00, 560.0D+00, 570.0D+00, &
    580.0D+00, 590.0D+00, 600.0D+00, 610.0D+00, 620.0D+00, &
    630.0D+00, 640.0D+00, 650.0D+00, 660.0D+00, 670.0D+00, &
    680.0D+00, 690.0D+00, 700.0D+00, 710.0D+00, 720.0D+00, &
    730.0D+00, 740.0D+00, 750.0D+00, 760.0D+00 /)

  if ( wdat(1) <= w .and. w <= wdat(ndat) ) then

    call interp ( ndat, w, wdat, r, rdat )
    call interp ( ndat, w, wdat, g, gdat )
    call interp ( ndat, w, wdat, b, bdat )

  else

    r = 0.0D+00
    g = 0.0D+00
    b = 0.0D+00

  end if

  return
end
subroutine nm_to_xyz ( w, x, y, z )

!*****************************************************************************80
!
!! NM_TO_XYZ converts a pure light wavelength to CIE xyz color coordinates.
!
!  Discussion:
!
!    The CIE xyz color coordinates are derived from the CIE X, Y, Z color
!    coordinates by the relationship:
!
!      x = X / ( X + Y + Z )
!      y = Y / ( X + Y + Z )
!      z = Z / ( X + Y + Z )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W, the wavelength of the pure light signal,
!    in nanometers.  W should lie between 380 nm and 780 nm.
!
!    Output, real ( kind = 8 ) X, Y, Z, the CIE xyz chromaticities.  These
!    lie between 0 and 1, and sum to 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndat = 81

  real ( kind = 8 ) w
  real ( kind = 8 ), save, dimension ( ndat ) :: wdat = (/ &
    380.0D+00, 385.0D+00, 390.0D+00, 395.0D+00, 400.0D+00, &
    405.0D+00, 410.0D+00, 415.0D+00, 420.0D+00, 425.0D+00, &
    430.0D+00, 435.0D+00, 440.0D+00, 445.0D+00, 450.0D+00, &
    455.0D+00, 460.0D+00, 465.0D+00, 470.0D+00, 475.0D+00, &
    480.0D+00, 485.0D+00, 490.0D+00, 495.0D+00, 500.0D+00, &
    505.0D+00, 510.0D+00, 515.0D+00, 520.0D+00, 525.0D+00, &
    530.0D+00, 535.0D+00, 540.0D+00, 545.0D+00, 550.0D+00, &
    555.0D+00, 560.0D+00, 565.0D+00, 570.0D+00, 575.0D+00, &
    580.0D+00, 585.0D+00, 590.0D+00, 595.0D+00, 600.0D+00, &
    605.0D+00, 610.0D+00, 615.0D+00, 620.0D+00, 625.0D+00, &
    630.0D+00, 635.0D+00, 640.0D+00, 645.0D+00, 650.0D+00, &
    655.0D+00, 660.0D+00, 665.0D+00, 670.0D+00, 675.0D+00, &
    680.0D+00, 685.0D+00, 690.0D+00, 695.0D+00, 700.0D+00, &
    705.0D+00, 710.0D+00, 715.0D+00, 720.0D+00, 725.0D+00, &
    730.0D+00, 735.0D+00, 740.0D+00, 745.0D+00, 750.0D+00, &
    755.0D+00, 760.0D+00, 765.0D+00, 770.0D+00, 775.0D+00, &
    780.0D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( ndat ) :: xdat = (/ &
    0.1741, 0.1740, 0.1738, 0.1736, 0.1733, &
    0.1730, 0.1726, 0.1721, 0.1714, 0.1703, &
    0.1689, 0.1669, 0.1644, 0.1611, 0.1566, &
    0.1510, 0.1440, 0.1355, 0.1241, 0.1096, &
    0.0913, 0.0687, 0.0454, 0.0235, 0.0082, &
    0.0039, 0.0139, 0.0389, 0.0743, 0.1142, &
    0.1547, 0.1929, 0.2296, 0.2658, 0.3016, &
    0.3373, 0.3731, 0.4087, 0.4441, 0.4788, &
    0.5125, 0.5448, 0.5752, 0.6029, 0.6270, &
    0.6482, 0.6658, 0.6801, 0.6915, 0.7006, &
    0.7079, 0.7140, 0.7190, 0.7230, 0.7260, &
    0.7283, 0.7300, 0.7311, 0.7320, 0.7327, &
    0.7334, 0.7340, 0.7344, 0.7346, 0.7347, &
    0.7347, 0.7347, 0.7347, 0.7347, 0.7347, &
    0.7347, 0.7347, 0.7347, 0.7347, 0.7347, &
    0.7347, 0.7347, 0.7347, 0.7347, 0.7347, &
    0.7347 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), save, dimension ( ndat ) :: ydat = (/ &
    0.0050, 0.0050, 0.0049, 0.0049, 0.0048, &
    0.0048, 0.0048, 0.0048, 0.0051, 0.0058, &
    0.0069, 0.0086, 0.0109, 0.0138, 0.0177, &
    0.0227, 0.0297, 0.0399, 0.0578, 0.0868, &
    0.1327, 0.2007, 0.2950, 0.4127, 0.5384, &
    0.6548, 0.7502, 0.8120, 0.8338, 0.8262, &
    0.8059, 0.7816, 0.7543, 0.7243, 0.6923, &
    0.6589, 0.6245, 0.5896, 0.5547, 0.5202, &
    0.4866, 0.4544, 0.4242, 0.3965, 0.3725, &
    0.3514, 0.3340, 0.3197, 0.3083, 0.2993, &
    0.2920, 0.2859, 0.2809, 0.2770, 0.2740, &
    0.2717, 0.2700, 0.2689, 0.2680, 0.2673, &
    0.2666, 0.2660, 0.2656, 0.2654, 0.2653, &
    0.2653, 0.2653, 0.2653, 0.2653, 0.2653, &
    0.2653, 0.2653, 0.2653, 0.2653, 0.2653, &
    0.2653, 0.2653, 0.2653, 0.2653, 0.2653, &
    0.2653 /)
  real ( kind = 8 ) z
  real ( kind = 8 ), save, dimension ( ndat ) :: zdat = (/ &
    0.8209, 0.8210, 0.8213, 0.8215, 0.8219, &
    0.8222, 0.8226, 0.8231, 0.8235, 0.8239, &
    0.8242, 0.8245, 0.8247, 0.8251, 0.8257, &
    0.8263, 0.8263, 0.8246, 0.8181, 0.8036, &
    0.7760, 0.7306, 0.6596, 0.5638, 0.4534, &
    0.3413, 0.2359, 0.1491, 0.0919, 0.0596, &
    0.0394, 0.0255, 0.0161, 0.0099, 0.0061, &
    0.0038, 0.0024, 0.0017, 0.0012, 0.0010, &
    0.0009, 0.0008, 0.0006, 0.0006, 0.0005, &
    0.0004, 0.0002, 0.0002, 0.0002, 0.0001, &
    0.0001, 0.0001, 0.0001, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000 /)

  if ( wdat(1) <= w .and. w <= wdat(ndat) ) then

    call interp ( ndat, w, wdat, x, xdat )
    call interp ( ndat, w, wdat, y, ydat )
    call interp ( ndat, w, wdat, z, zdat )

  else

    x = 0.0D+00
    y = 0.0D+00
    z = 0.0D+00

  end if

  return
end
subroutine nm_to_xyz_cap ( w, xcap, ycap, zcap )

!*****************************************************************************80
!
!! NM_TO_XYZ_CAP converts a pure light wavelength to CIE XYZ color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!    Measurements were made assuming a stimulus of one watt of pure
!    light of the indicated wavelength.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Deane Judd and Gunter Wyszecki,
!    Color in Business, Science, and Industry,
!    Wiley, 1975, pages 126-127, 130.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W, the wavelength in nanometers.  W should lie
!    between 380 nm and 780 nm.  Values outside this range will result in
!    output values of XCAP = YCAP = ZCAP = 0.
!
!    Output, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!    The values indicate the amounts of the CIE primaries X, Y and Z
!    required to match the color of the stimulus light.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndat = 90

  real ( kind = 8 ) w
  real ( kind = 8 ), save, dimension ( ndat ) :: wdat = (/ &
    380.0D+00, 385.0D+00, 390.0D+00, 395.0D+00, 400.0D+00, &
    405.0D+00, 410.0D+00, 415.0D+00, 420.0D+00, 425.0D+00, &
    430.0D+00, 435.0D+00, 440.0D+00, 445.0D+00, 450.0D+00, &
    455.0D+00, 460.0D+00, 465.0D+00, 470.0D+00, 475.0D+00, &
    480.0D+00, 485.0D+00, 490.0D+00, 495.0D+00, 500.0D+00, &
    505.0D+00, 510.0D+00, 515.0D+00, 520.0D+00, 525.0D+00, &
    530.0D+00, 535.0D+00, 540.0D+00, 545.0D+00, 550.0D+00, &
    555.0D+00, 560.0D+00, 565.0D+00, 570.0D+00, 575.0D+00, &
    580.0D+00, 585.0D+00, 590.0D+00, 595.0D+00, 600.0D+00, &
    605.0D+00, 610.0D+00, 615.0D+00, 620.0D+00, 625.0D+00, &
    630.0D+00, 635.0D+00, 640.0D+00, 645.0D+00, 650.0D+00, &
    655.0D+00, 660.0D+00, 665.0D+00, 670.0D+00, 675.0D+00, &
    680.0D+00, 685.0D+00, 690.0D+00, 695.0D+00, 700.0D+00, &
    705.0D+00, 710.0D+00, 715.0D+00, 720.0D+00, 725.0D+00, &
    730.0D+00, 735.0D+00, 740.0D+00, 745.0D+00, 750.0D+00, &
    755.0D+00, 760.0D+00, 765.0D+00, 770.0D+00, 775.0D+00, &
    780.0D+00, 785.0D+00, 790.0D+00, 795.0D+00, 800.0D+00, &
    805.0D+00, 810.0D+00, 815.0D+00, 820.0D+00, 825.0D+00 /)
  real ( kind = 8 ) xcap
  real ( kind = 8 ), save, dimension ( ndat ) :: xdat = (/ &
    0.0014, 0.0022, 0.0042, 0.0076, 0.0143, &
    0.0232, 0.0435, 0.0776, 0.1344, 0.2148, &
    0.2839, 0.3285, 0.3483, 0.3481, 0.3362, &
    0.3187, 0.2908, 0.2511, 0.1954, 0.1421, &
    0.0956, 0.0580, 0.0320, 0.0147, 0.0049, &
    0.0024, 0.0093, 0.0291, 0.0633, 0.1096, &
    0.1655, 0.2257, 0.2904, 0.3597, 0.4334, &
    0.5121, 0.5945, 0.6784, 0.7621, 0.8425, &
    0.9163, 0.9786, 1.0263, 1.0567, 1.0622, &
    1.0456, 1.0026, 0.9384, 0.8544, 0.7514, &
    0.6424, 0.5419, 0.4479, 0.3608, 0.2835, &
    0.2187, 0.1649, 0.1212, 0.0874, 0.0636, &
    0.0468, 0.0329, 0.0227, 0.0158, 0.0114, &
    0.0081, 0.0058, 0.0041, 0.0029, 0.0020, &
    0.0014, 0.0010, 0.0007, 0.0005, 0.0003, &
    0.0002, 0.0002, 0.0001, 0.0001, 0.0001, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000 /)
  real ( kind = 8 ) ycap
  real ( kind = 8 ), save, dimension ( ndat ) :: ydat = (/ &
    0.0000, 0.0001, 0.0001, 0.0002, 0.0004, &
    0.0006, 0.0012, 0.0022, 0.0040, 0.0073, &
    0.0116, 0.0168, 0.0230, 0.0298, 0.0380, &
    0.0480, 0.0600, 0.0739, 0.0910, 0.1126, &
    0.1390, 0.1693, 0.2080, 0.2586, 0.3230, &
    0.4073, 0.5030, 0.6082, 0.7100, 0.7932, &
    0.8620, 0.9149, 0.9540, 0.9803, 0.9950, &
    1.0000, 0.9950, 0.9786, 0.9520, 0.9154, &
    0.8700, 0.8163, 0.7570, 0.6949, 0.6310, &
    0.5668, 0.5030, 0.4412, 0.3810, 0.3210, &
    0.2650, 0.2170, 0.1750, 0.1382, 0.1070, &
    0.0816, 0.0610, 0.0466, 0.0320, 0.0232, &
    0.0170, 0.0119, 0.0082, 0.0057, 0.0041, &
    0.0029, 0.0021, 0.0015, 0.0010, 0.0007, &
    0.0005, 0.0004, 0.0002, 0.0002, 0.0001, &
    0.0001, 0.0001, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000 /)
  real ( kind = 8 ) zcap
  real ( kind = 8 ), save, dimension ( ndat ) :: zdat = (/ &
    0.0065, 0.0105, 0.0201, 0.0362, 0.0679, &
    0.1102, 0.2074, 0.3713, 0.6456, 1.0391, &
    1.3856, 1.6230, 1.7471, 1.7826, 1.7721, &
    1.7441, 1.6692, 1.5281, 1.2876, 1.0419, &
    0.8130, 0.6162, 0.4652, 0.3533, 0.2720, &
    0.2123, 0.1582, 0.1117, 0.0782, 0.0573, &
    0.0422, 0.0298, 0.0203, 0.0134, 0.0087, &
    0.0057, 0.0039, 0.0027, 0.0021, 0.0018, &
    0.0017, 0.0014, 0.0011, 0.0010, 0.0008, &
    0.0006, 0.0003, 0.0002, 0.0002, 0.0001, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000 /)

  if ( wdat(1) <= w .and. w <= wdat(ndat) ) then

    call interp ( ndat, w, wdat, xcap, xdat )
    call interp ( ndat, w, wdat, ycap, ydat )
    call interp ( ndat, w, wdat, zcap, zdat )

  else

    xcap = 0.0D+00
    ycap = 0.0D+00
    zcap = 0.0D+00

  end if

  return
end
subroutine nonlin_to_lin ( rprime, r )

!*****************************************************************************80
!
!! NONLIN_TO_LIN converts a nonlinear video signal to a linear light intensity.
!
!  Discussion:
!
!    This process is also known as gamma (un)correction.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RPRIME, the nonlinear video signal.
!
!    Output, real ( kind = 8 ) R, the linear light intensity.
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) rprime

  if ( rprime <= -0.081D+00 ) then
    r = - ( ( ( 0.099D+00 - rprime  ) / 1.099D+00 )**(1.0D+00/0.45D+00) )
  else if ( abs ( rprime ) <= 0.081D+00 ) then
    r = rprime / 4.5D+00
  else if ( 0.081D+00 <= rprime ) then
    r = ( ( 0.099D+00 + rprime ) / 1.099D+00 )**(1.0D+00/0.45D+00)
  end if

  return
end
subroutine primaries_to_y ( rx, ry, gx, gy, bx, by, wx, wy, yr, yg, yb )

!*****************************************************************************80
!
!! PRIMARIES_TO_Y computes the luminance function for given primaries.
!
!  Discussion:
!
!    The luminance function has the form:
!
!      Y = YR * R + YG * G + YB * B
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Martindale and Alan Paeth,
!    Television Color Encoding and Hot Broadcast Colors,
!    Graphics Gems II, edited by James Arvo,
!    Academic Press, 1991, pages 147-158.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RX, RY, the xy chromaticities of the R primary.
!
!    Input, real ( kind = 8 ) GX, GY, the xy chromaticities of the G primary.
!
!    Input, real ( kind = 8 ) BX, BY, the xy chromaticities of the B primary.
!
!    Input, real ( kind = 8 ) WX, WY, the xy chromaticities of the reference
!    white.
!
!    Output, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: nrhs = 1

  real ( kind = 8 ) a(n,n+nrhs)
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) gx
  real ( kind = 8 ) gy
  integer ( kind = 4 ) info
  real ( kind = 8 ) rx
  real ( kind = 8 ) ry
  real ( kind = 8 ) wx
  real ( kind = 8 ) wy
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yr
!
!  Set up the coefficients and right hand side of the linear system.
!
  a(1,1) = rx
  a(1,2) = gx
  a(1,3) = bx
  a(1,4) = wx / wy

  a(2,1) = ry
  a(2,2) = gy
  a(2,3) = by
  a(2,4) = wy / wy

  a(3,1) = 1.0D+00 - rx - ry
  a(3,2) = 1.0D+00 - gx - gy
  a(3,3) = 1.0D+00 - bx - by
  a(3,4) = ( 1.0D+00 - wx - wy ) / wy
!
!  Solve the linear system A * x = b.
!
  call r8mat_solve ( a, n, nrhs, info )
!
!  Extract the solution.
!
  if ( info == 0 ) then
    yr = a(1,4) * ry
    yg = a(2,4) * gy
    yb = a(3,4) * by
  else
    yr = 0.0D+00
    yg = 0.0D+00
    yb = 0.0D+00
  end if

  return
end
function r8_cubert ( x )

!*****************************************************************************80
!
!! R8_CUBERT returns the cube root of an R8.
!
!  Discussion:
!
!    R8_CUBERT is designed to avoid the possible problems that can occur
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
!    Output, real ( kind = 8 ) R8_CUBERT, the cube root of X.
!
  implicit none

  real ( kind = 8 ) r8_cubert
  real ( kind = 8 ) x

  if ( 0.0D+00 < x ) then
    r8_cubert = x**(1.0D+00/3.0D+00)
  else if ( x == 0.0D+00 ) then
    r8_cubert = 0.0D+00
  else
    r8_cubert = - ( abs ( x ) )**(1.0D+00/3.0D+00)
  end if

  return
end
function r8_modp ( x, y )

!*****************************************************************************80
!
!! R8_MODP returns the nonnegative remainder of real division.
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
!        I         J     MOD   R8_MODP  R8_MODP Factorization
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
!    29 August 1998
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
!    Output, real ( kind = 8 ) R8_MODP, the nonnegative remainder when
!    X is divided by Y.
!
  implicit none

  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( y == 0 ) then
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
subroutine r8mat_solve ( a, n, nrhs, info )

!*****************************************************************************80
!
!! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(N,N+NRHS), contains in rows and columns 1
!    to N the coefficient matrix, and in columns N+1 through
!    N+NRHS, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides.  NRHS
!    must be at least 0.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrhs

  real ( kind = 8 ) a(n,n+nrhs)
  real ( kind = 8 ) apivot
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot
  integer ( kind = 4 ) j

  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j + 1, n
      if ( abs ( apivot ) < abs ( a(i,j) ) ) then
        apivot = a(i,j)
        ipivot = i
      end if
    end do

    if ( apivot == 0.0D+00 ) then
      info = j
      return
    end if
!
!  Interchange.
!
    do i = 1, n + nrhs
      call r8_swap ( a(ipivot,i), a(j,i) )
    end do
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0D+00
    a(j,j+1:n+nrhs) = a(j,j+1:n+nrhs) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then

        factor = a(i,j)
        a(i,j) = 0.0D+00
        a(i,j+1:n+nrhs) = a(i,j+1:n+nrhs) - factor * a(j,j+1:n+nrhs)

      end if

    end do

  end do

  return
end
function radians_to_degrees ( angle )

!*****************************************************************************80
!
!! RADIANS_TO_DEGREES converts an angle from radians to degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, an angle in radians.
!
!    Output, real ( kind = 8 ) RADIANS_TO_DEGREES, the equivalent angle
!    in degrees.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  radians_to_degrees = ( angle / pi ) * 180.0D+00

  return
end
subroutine rgb709_to_xyz_cap ( r, g, b, xcap, ycap, zcap )

!*****************************************************************************80
!
!! RGB709_TO_XYZ_CAP converts RGB709 to CIE XYZ color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB 709 color coordinates to be
!    converted.
!
!    Output, real ( kind = 8 ) XCAP, YCAP, ZCAP, the corresponding CIE XYZ
!    color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) r
  real ( kind = 8 ) xcap
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap

  xcap = 0.412453D+00 * r + 0.35758D+00  * g + 0.180423D+00 * b
  ycap = 0.212671D+00 * r + 0.71516D+00  * g + 0.072169D+00 * b
  zcap = 0.019334D+00 * r + 0.119193D+00 * g + 0.950227D+00 * b

  return
end
subroutine rgb_check ( r, g, b )

!*****************************************************************************80
!
!! RGB_CHECK corrects out-of-range RGB color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!  Example:
!
!    Black   ( 0.0, 0.0, 0.0D+00 )
!    Blue    ( 0.0, 0.0, 1.0D+00 )
!    Cyan    ( 0.0, 1.0, 1.0D+00 )
!    Green   ( 0.0, 1.0, 0.0D+00 )
!    Magenta ( 1.0, 0.0, 1.0D+00 )
!    Red     ( 1.0, 0.0, 0.0D+00 )
!    White   ( 1.0, 1.0, 1.0D+00 )
!    Yellow  ( 1.0, 1.0, 0.0D+00 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) R, G, B, the RGB color coordinates to
!    be checked.  Any coordinate less than 0 is set to 0, and any
!    coordinate greater than 1 is set to 1.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) r

  r = max ( r, 0.0D+00 )
  r = min ( r, 1.0D+00 )

  g = max ( g, 0.0D+00 )
  g = min ( g, 1.0D+00 )

  b = max ( b, 0.0D+00 )
  b = min ( b, 1.0D+00 )

  return
end
subroutine rgb_uniform ( seed, r, g, b )

!*****************************************************************************80
!
!! RGB_UNIFORM returns a random RGB color.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R, G, B, the RGB values of the color.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) g
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed

  r = r8_uniform_01 ( seed )
  g = r8_uniform_01 ( seed )
  b = r8_uniform_01 ( seed )

  return
end
subroutine rgb_named_uniform ( seed, r, g, b, name )

!*****************************************************************************80
!
!! RGB_NAMED_UNIFORM returns a random named RGB color.
!
!  Discussion:
!
!    The named colors are extracted at random from those listed
!    in "colors.txt".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R, G, B, the RGB values of the color.
!
!    Output, character ( len = * ) NAME, the name of the color.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ), save, allocatable, dimension ( : ) :: b_saved
  integer ( kind = 4 ) color_index
  integer ( kind = 4 ), save :: color_num
  character ( len = 30 ), save :: color_file = 'colors.txt'
  real ( kind = 8 ) g
  real ( kind = 8 ), save, allocatable, dimension ( : ) :: g_saved
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ig
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) iunit
  logical, save :: loaded = .false.
  character ( len = * ) name
  character ( len = 30 ) name2
  character ( len = 30 ), save, allocatable, dimension ( : ) :: name_saved
  real ( kind = 8 ) r
  real ( kind = 8 ), save, allocatable, dimension ( : ) :: r_saved
  integer ( kind = 4 ) seed

  if ( .not. loaded ) then

    call get_unit ( iunit )

    open ( unit = iunit, file = color_file, status = 'old', iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RGB_NAMED_UNIFORM - Fatal error!'
      write ( *, '(a)' ) '  Could not open the color name file:'
      write ( *, '(a)' ) trim ( color_file )
      stop
    end if

    color_num = 0

    do

      read ( iunit, '(i3,2x,i3,2x,i3,2x,a)', iostat = ios ) ir, ig, ib, name

      if ( ios /= 0 ) then
        exit
      end if

      color_num = color_num + 1

    end do

    allocate ( r_saved(1:color_num) )
    allocate ( g_saved(1:color_num) )
    allocate ( b_saved(1:color_num) )
    allocate ( name_saved(1:color_num) )

    rewind ( unit = iunit )

    color_num = 0

    do

      read ( iunit, '(i3,2x,i3,2x,i3,2x,a)', iostat = ios ) ir, ig, ib, name2

      if ( ios /= 0 ) then
        exit
      end if

      color_num = color_num + 1

      r_saved(color_num) = real ( ir, kind = 8 ) / 255.0D+00
      g_saved(color_num) = real ( ig, kind = 8 ) / 255.0D+00
      b_saved(color_num) = real ( ib, kind = 8 ) / 255.0D+00
      name_saved(color_num) = name2

    end do

    close ( unit = iunit )

    loaded = .true.

  end if

  color_index = i4_uniform ( 1, color_num, seed )

  r = r_saved(color_index)
  g = g_saved(color_index)
  b = b_saved(color_index)
  name = name_saved(color_index)

  return
end
subroutine rgb_test ( itest, rtest, gtest, btest )

!*****************************************************************************80
!
!! RGB_TEST supplies RGB values for tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ITEST, the index of the test.
!
!    Output, real ( kind = 8 ) RTEST, GTEST, BTEST, sample RGB color coordinate
!    values for testing.  If ITEST is outside the range of data,
!    then RTEST = GTEST = BTEST = -1.0.
!
  implicit none

  real ( kind = 8 ) btest
  real ( kind = 8 ) gtest
  integer ( kind = 4 ) itest
  real ( kind = 8 ) rtest

  if ( itest == 1 ) then

    rtest = 0.9D+00
    gtest = 0.0D+00
    btest = 0.0D+00

  else if ( itest == 2 ) then

    rtest = 0.0D+00
    gtest = 0.8D+00
    btest = 0.0D+00

  else if ( itest == 3 ) then

    rtest = 0.0D+00
    gtest = 0.0D+00
    btest = 0.7D+00

  else if ( itest == 4 ) then

    rtest = 0.0D+00
    gtest = 0.6D+00
    btest = 0.6D+00

  else if ( itest == 5 ) then

    rtest = 0.5D+00
    gtest = 0.0D+00
    btest = 0.5D+00

  else if ( itest == 6 ) then

    rtest = 0.4D+00
    gtest = 0.4D+00
    btest = 0.0D+00
!
!  Dark gray
!
  else if ( itest == 7 ) then

    rtest = 0.3D+00
    gtest = 0.3D+00
    btest = 0.3D+00
!
!  Black
!
  else if ( itest == 8 ) then

    rtest = 0.0D+00
    gtest = 0.0D+00
    btest = 0.0D+00
!
!  White
!
  else if ( itest == 9 ) then

    rtest = 1.0D+00
    gtest = 1.0D+00
    btest = 1.0D+00

  else if ( itest == 10 ) then

    rtest = 0.1D+00
    gtest = 0.3D+00
    btest = 0.5D+00

  else if ( itest == 11 ) then

    rtest = 0.3D+00
    gtest = 0.5D+00
    btest = 0.3D+00
!
!  Signal that ITEST is out of range.
!
  else

    rtest = -1.0D+00
    gtest = -1.0D+00
    btest = -1.0D+00

  end if

  return
end
subroutine rgb_to_cmy ( r, g, b, c, m, y )

!*****************************************************************************80
!
!! RGB_TO_CMY converts RGB to CMY color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!    The CMY color system describes a color based on the amounts of the
!    base colors cyan, magenta, and yellow.  Thus, a particular color
!    has three coordinates, (C,M,Y).  Each coordinate must be between
!    0 and 1.  Black is (1,1,1) and white is (0,0,0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB color coordinates to be
!    converted.
!
!    Output, real ( kind = 8 ) C, M, Y, the corresponding CMY color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) g
  real ( kind = 8 ) m
  real ( kind = 8 ) r
  real ( kind = 8 ) y

  c = 1.0D+00 - r
  m = 1.0D+00 - g
  y = 1.0D+00 - b

  return
end
subroutine rgb_to_cmyk ( r, g, b, c, m, y, k )

!*****************************************************************************80
!
!! RGB_TO_CMYK converts RGB to CMYK color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!    The CMYK color system describes a color based on the amounts of the
!    base colors cyan, magenta, yellow, and black.  The CMYK system is
!    based on the CMY system, except that equal amounts of C, M, and Y
!    are replaced by the single color K.  Thus, a particular color
!    has four coordinates, (C,M,Y,K).  Each coordinate must be between
!    0 and 1, and it must also be true that C+K, M+K and Y+K are
!    each no greater than 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB color coordinates to be
!    converted.
!
!    Output, real ( kind = 8 ) C, M, Y, K, the corresponding CMYK color
!    coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) g
  real ( kind = 8 ) k
  real ( kind = 8 ) m
  real ( kind = 8 ) r
  real ( kind = 8 ) y
!
!  Compute the CMY equivalent colors.
!
  c = 1.0D+00 - r
  m = 1.0D+00 - g
  y = 1.0D+00 - b
!
!  Compute the black component.
!
  k = min ( c, m, y )
!
!  Subtract off the black component to complete the CMYK specification.
!
  c = c - k
  m = m - k
  y = y - k

  return
end
subroutine rgb_to_hls ( r, g, b, h, l, s )

!*****************************************************************************80
!
!! RGB_TO_HLS converts RGB to HLS color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!    The HLS color system describes a color based on the qualities of
!    hue, lightness, and saturation.  A particular color has three
!    coordinates, (H,L,S).  The L and S coordinates must be between
!    0 and 1, while the H coordinate must be between 0 and 360, and
!    is interpreted as an angle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB color coordinates
!    to be converted.
!
!    Output, real ( kind = 8 ) H, L, S, the corresponding HLS color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) bc
  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) g
  real ( kind = 8 ) gc
  real ( kind = 8 ) h
  real ( kind = 8 ) l
  real ( kind = 8 ) r
  real ( kind = 8 ) rc
  real ( kind = 8 ) rgbmax
  real ( kind = 8 ) rgbmin
  real ( kind = 8 ) s
!
!  Compute lightness.
!
  rgbmax = max ( r, g, b )
  rgbmin = min ( r, g, b )
  l = ( rgbmax + rgbmin ) / 2.0D+00
!
!  Compute saturation.
!
  if ( rgbmax == rgbmin ) then

    s = 0.0D+00

  else

    if ( l <= 0.5D+00 ) then
      s = ( rgbmax - rgbmin ) / ( rgbmax + rgbmin )
    else
      s = ( rgbmax - rgbmin ) / ( 2.0D+00 - rgbmax - rgbmin )
    end if

  end if
!
!  Compute the hue.
!
  if ( rgbmax == rgbmin ) then

    h = 0.0D+00

  else

    rc = ( rgbmax - r ) / ( rgbmax - rgbmin )
    gc = ( rgbmax - g ) / ( rgbmax - rgbmin )
    bc = ( rgbmax - b ) / ( rgbmax - rgbmin )

    if ( r == rgbmax ) then
      h = bc - gc
    else if ( g == rgbmax ) then
      h = 2.0D+00 + rc - bc
    else
      h = 4.0D+00 + gc - rc
    end if

    h = h * 60.0D+00
!
!  Make sure H lies between 0 and 360.0.
!
    h = r8_modp ( h, 360.0D+00 )

  end if

  return
end
subroutine rgb_to_hsi ( r, g, b, h, s, i )

!*****************************************************************************80
!
!! RGB_TO_HSI converts RGB to HSI color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!    The HSI color system uses coordinates of
!      Hue, an angle between 0 and 360, (0=R, 120=G, 240=B)
!      Saturation, between 0 and 1, and
!      Intensity, between 0 and 1.
!
!    Note, Grrr, that there is a typo in the formula for H in Fortner's
!    book.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!    Brand Fortner,
!    Number by Colors, A Guide to Using Color to Understand Technical Data,
!    Springer, 1997, pages 140-143.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB color coordinates to be
!    converted.
!
!    Output, real ( kind = 8 ) H, S, I, the corresponding HSI color coordinates.
!
  implicit none

  real ( kind = 8 ) atan4
  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  real ( kind = 8 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ) s

  h = atan4 ( sqrt ( 3.0D+00 ) * ( g - b ), 2.0D+00 * r - b - g )
  h = radians_to_degrees ( h )

  i = ( r + g + b ) / 3.0D+00

  if ( i == 0.0D+00 ) then
    s = 0.0D+00
  else
    s = 1.0D+00 - min ( r, g, b ) / i
  end if

  return
end
subroutine rgb_to_hsv ( r, g, b, h, s, v )

!*****************************************************************************80
!
!! RGB_TO_HSV converts RGB to HSV color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!    The HSV color system describes a color based on the three qualities
!    of hue, saturation, and value.  A given color will be represented
!    by three numbers, (H,S,V).  H, the value of hue, is an angle
!    between 0 and 360 degrees, with 0 representing red.  S is the
!    saturation, and is between 0 and 1.  Finally, V is the "value",
!    a measure of brightness, which goes from 0 for black, increasing
!    to a maximum of 1 for the brightest colors.  The HSV color system
!    is sometimes also called HSB, where the B stands for brightness.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB color coordinates
!     to be converted.
!
!    Output, real ( kind = 8 ) H, S, V, the corresponding HSV color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) bc
  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) g
  real ( kind = 8 ) gc
  real ( kind = 8 ) h
  real ( kind = 8 ) r
  real ( kind = 8 ) rc
  real ( kind = 8 ) rgbmax
  real ( kind = 8 ) rgbmin
  real ( kind = 8 ) s
  real ( kind = 8 ) v

  rgbmax = max ( r, g, b )
  rgbmin = min ( r, g, b )

  v = rgbmax
!
!  Compute the saturation.
!
  if ( rgbmax /= 0.0D+00 ) then
    s = ( rgbmax - rgbmin ) / rgbmax
  else
    s = 0.0D+00
  end if
!
!  Compute the hue.
!
  if ( s == 0.0D+00 ) then

    h = 0.0D+00

  else

    rc = ( rgbmax - r ) / ( rgbmax - rgbmin )
    gc = ( rgbmax - g ) / ( rgbmax - rgbmin )
    bc = ( rgbmax - b ) / ( rgbmax - rgbmin )

    if ( r == rgbmax ) then
      h = bc - gc
    else if ( g == rgbmax ) then
      h = 2.0D+00 + rc - bc
    else
      h = 4.0D+00 + gc - rc
    end if

    h = h * 60.0D+00
!
!  Make sure H lies between 0 and 360.0D+00
!
    h = r8_modp ( h, 360.0D+00 )

  end if

  return
end
subroutine rgb_to_hue ( r, g, b, h )

!*****************************************************************************80
!
!! RGB_TO_HUE converts (R,G,B) colors to a hue value between 0 and 1.
!
!  Discussion:
!
!    The hue computed here should be the same as the H value computed
!    for HLS and HSV, except that it ranges from 0 to 1 instead of
!    0 to 360.
!
!    A monochromatic color ( white, black, or a shade of gray) does not
!    have a hue.  This routine will return a special value of H = -1
!    for such cases.
!
!  Example:
!
!    Color    R    G    B     H
!
!    red      1.0  0.0  0.0   0.00
!    yellow   1.0  1.0  0.0   0.16
!    green    0.0  1.0  0.0   0.33
!    cyan     0.0  1.0  1.0   0.50
!    blue     0.0  0.0  1.0   0.67
!    magenta  1.0  0.0  1.0   0.83
!
!    black    0.0  0.0  0.0  -1.00
!    gray     0.5  0.5  0.5  -1.00
!    white    1.0  1.0  1.0  -1.00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the red, green and blue values of
!    the color.
!    These values should be between 0 and 1.
!
!    Output, real ( kind = 8 ) H, the corresponding hue of the color,
!    or -1.0D+00 if the color is monochromatic.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  real ( kind = 8 ) h
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) rgbmax
  real ( kind = 8 ) rgbmin
!
!  Make sure the colors are between 0 and 1.
!
  r2 = min ( max ( r, 0.0D+00 ), 1.0D+00 )
  g2 = min ( max ( g, 0.0D+00 ), 1.0D+00 )
  b2 = min ( max ( b, 0.0D+00 ), 1.0D+00 )
!
!  Compute the minimum and maximum of R, G and B.
!
  rgbmax = r2
  rgbmax = max ( rgbmax, g2 )
  rgbmax = max ( rgbmax, b2 )

  rgbmin = r2
  rgbmin = min ( rgbmin, g2 )
  rgbmin = min ( rgbmin, b2 )
!
!  If RGBMAX = RGBMIN, then the color has no hue.
!
  if ( rgbmax == rgbmin ) then

    h = - 1.0D+00
!
!  Otherwise, we need to determine the dominant color.
!
  else

    if ( r2 == rgbmax ) then
      h = ( g2 - b2 ) / ( rgbmax - rgbmin )
    else if ( g2 == rgbmax ) then
      h = 2.0D+00 + ( b2 - r2 ) / ( rgbmax - rgbmin )
    else if ( b2 == rgbmax ) then
      h = 4.0D+00 + ( r2 - g2 ) / ( rgbmax - rgbmin )
    end if

    h = h / 6.0D+00
!
!  Make sure H lies between 0 and 1.0.
!
    if ( h < 0.0D+00 ) then
      h = h + 1.0D+00
    else if ( 1.0D+00 < h ) then
      h = h - 1.0D+00
    end if

  end if

  return
end
subroutine rgb_to_name ( r, g, b, name )

!*****************************************************************************80
!
!! RGB_TO_NAME converts RGB colors to the name of the nearest color.
!
!  Discussion:
!
!    The names and information are read from the file "COLORS.TXT", a
!    modified version of the X Windows color data file "RGB.TXT".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB coordinates of the color.
!
!    Output, character ( len = * ) NAME, the name of the color that is the
!    "closest" to the given RGB coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) bcolor
  character ( len = 20 ) :: color_file = 'colors.txt'
  real ( kind = 8 ) dismin
  real ( kind = 8 ) dist
  logical first
  real ( kind = 8 ) g
  real ( kind = 8 ) gcolor
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ig
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) iunit
  character ( len = * ) name
  character ( len = 30 ) namecolor
  real ( kind = 8 ) r
  real ( kind = 8 ) rcolor

  first = .true.
  name = 'Unknown'

  call get_unit ( iunit )

  open ( unit = iunit, file = color_file, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RGB_TO_NAME - Fatal error!'
    write ( *, '(a)' ) '  Could not open the color name file:'
    write ( *, '(a)' ) trim ( color_file )
    stop
  end if

  do

    read ( iunit, '(i3,2x,i3,2x,i3,2x,a)', iostat = ios ) ir, ig, ib, namecolor

    if ( ios /= 0 ) then
      exit
    end if

    rcolor = real ( ir, kind = 8 ) / 255.0D+00
    gcolor = real ( ig, kind = 8 ) / 255.0D+00
    bcolor = real ( ib, kind = 8 ) / 255.0D+00

    dist = sqrt ( ( r - rcolor )**2 + ( g - gcolor ) **2 + ( b - bcolor ) **2 )

    if ( first ) then
      dismin = dist
      name = namecolor
      first = .false.
    else if ( dist < dismin ) then
      dismin = dist
      name = namecolor
    end if

    if ( dismin == 0.0D+00 ) then
      exit
    end if

  end do

  close ( unit = iunit )

  return
end
subroutine rgb_to_ncs ( r, g, b, c1, c2, n, c, s )

!*****************************************************************************80
!
!! RGB_TO_NCS converts RGB to NCS color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!
!    The NCS or "natural color system" describes a color based on:
!    * C1 and C2, two elementary colors from the sequence RYGB or
!      C2 = blank for a pure elementary color, or
!      C1 = N, C2 = blank for a neutral color);
!    * N, the percentage of C2;
!    * C, the colorfulness or strength, as a percentage;
!    * S, the blackness as a percentage.
!
!    The scant documentation I have seen claims that the percentages are
!    always less than 100.  I don't see why, and for now I'll let them
!    lie between 0 and 100.  The NCS designation for a color has the form
!    "CCSS C1NC2".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Olof Kylander and Karin Kylander,
!    GIMP: The Official Manual,
!    Coriolis Open Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
!    Output, character C1, C2, integer N, C, S, the NCS color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  integer ( kind = 4 ) c
  character c1
  character c2
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) rgb_min
  integer ( kind = 4 ) s
  integer ( kind = 4 ) w
  real ( kind = 8 ) y2
!
!  Copy the input data.
!
  r2 = r
  g2 = g
  b2 = b
  y2 = 0.0D+00

  call rgb_check ( r2, g2, b2 )
!
!  Extract the whiteness value.
!
  rgb_min = min ( r2, g2, b2 )
  w = int ( 100.0D+00 * rgb_min )
!
!  Subtract off the whiteness.
!
  r2 = r2 - rgb_min
  g2 = g2 - rgb_min
  b2 = b2 - rgb_min
  y2 = 0.0D+00
!
!  If R2 = G2 = B2 = zero, we have an N (neutral) color or black or white.
!
  if ( r2 == 0.0D+00 .and. g2 == 0.0D+00 .and. b2 == 0.0D+00 ) then

    if ( w == 0 ) then
      c1 = 'B'
      c2 = ' '
      n = 0
      c = 0
    else if ( w == 100 ) then
      c1 = 'W'
      c2 = ' '
      n = 0
      c = 0
    else
      c1 = 'N'
      c2 = ' '
      n = 0
      c = 0
    end if
!
!  If two colors are zero, or if G = R, then we have a pure color.
!
  else if ( r2 == 0.0D+00 .and. b2 == 0.0D+00 ) then

    c1 = 'G'
    c2 = ' '
    n = 0
    c = int ( 100.0D+00 * g2 )

  else if ( r2 == 0.0D+00 .and. g2 == 0.0D+00 ) then

    c1 = 'B'
    c2 = ' '
    n = 0
    c = int ( 100.0D+00 * b2 )

  else if ( b2 == 0.0D+00 .and. g2 == 0.0D+00 ) then

    c1 = 'R'
    c2 = ' '
    n = 0
    c = int ( 100.0D+00 * r2 )

  else if ( g2 == r2 ) then

    y2 = g2
    g2 = 0.0D+00
    r2 = 0.0D+00

    c1 = 'Y'
    c2 = ' '
    n = 0
    c = int ( 100.0D+00 * y2 )
!
!  Two colors are nonzero.
!
  else if ( r2 == 0.0D+00 ) then
    c1 = 'G'
    c2 = 'B'
    y2 = 0.0
    n = int ( 100.0D+00 * b2 / ( g2 + b2 ) )

  else if ( g2 == 0.0D+00 ) then

    c1 = 'B'
    c2 = 'R'
    y2 = 0.0D+00
    n = int ( 100.0D+00 * r2 / ( b2 + r2 ) )

  else if ( b2 == 0.0D+00 .and. g2 < r2 ) then

    c1 = 'R'
    c2 = 'Y'

    y2 = g2
    r2 = r2 - g2
    g2 = 0.0D+00

    n = int ( 100.0D+00 * y2 / ( r2 + y2 ) )

  else if ( b2 == 0.0D+00 .and. r2 < g2 ) then

    c1 = 'Y'
    c2 = 'G'

    y2 = r2
    g2 = g2 - r2
    r2 = 0.0D+00

    n = int ( 100.0D+00 * g2 / ( y2 + g2 ) )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RGB_TO_NCS - Fatal error!'
    write ( *, '(a)' ) '  Unexpected RGB combination:'
    write ( *, '(3g14.6)' ) r2, g2, b2

  end if

  c = int ( 100.0D+00 * max ( r2, g2, b2, y2 ) )
  s = 100 - c - w

  return
end
subroutine rgb_to_rgbprime ( r, g, b, rprime, gprime, bprime )

!*****************************************************************************80
!
!! RGB_TO_RGBPRIME converts RGB to R'G'B' color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!    The R'G'B' color system is a nonlinear video signal measurement.
!    Each coordinate must be between 0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB color coordinates to be
!    converted.
!
!    Output, real ( kind = 8 ) RPRIME, GPRIME, BPRIME, the corresponding
!    R'G'B' color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) bprime
  real ( kind = 8 ) g
  real ( kind = 8 ) gprime
  real ( kind = 8 ) r
  real ( kind = 8 ) rprime

  call lin_to_nonlin ( r, rprime )
  call lin_to_nonlin ( g, gprime )
  call lin_to_nonlin ( b, bprime )

  return
end
subroutine rgb_to_ycbcr ( r, g, b, yr, yg, yb, yprime, cb, cr )

!*****************************************************************************80
!
!! RGB_TO_YCBCR converts RGB to Y'CbCr color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!    The Y'CbCr color system is used in digital television signals.
!    The Y' component measures luma, an approximation to the luminance
!    or amount of light.  Y' is the only component displayed on black
!    and white televisions.  The Cb and Cr components contain measures
!    of the blue and red components of the color.  Y' should be between
!    0 and 1, while Cb and Cr should be between -0.5 and 0.5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB color coordinates to be
!    converted.
!
!    Input, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
!    Output, real ( kind = 8 ) YPRIME, CB, CR, the corresponding Y'CbCr
!    color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) cb
  real ( kind = 8 ) cr
  real ( kind = 8 ) g
  real ( kind = 8 ) r
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr

  yprime = yr * r + yg * g + yb * b
  cb = 0.5D+00 * ( b - yprime ) / ( 1.0D+00 - yb )
  cr = 0.5D+00 * ( r - yprime ) / ( 1.0D+00 - yr )

  return
end
subroutine rgb_to_yiq ( r, g, b, yr, yg, yb, yprime, i, q )

!*****************************************************************************80
!
!! RGB_TO_YIQ converts RGB to Y'IQ color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!    Y'IQ colors are used in American NTSC commercial color television
!    broadcasting.  The Y' component measures luma, an approximation to
!    luminance, or the amount of light.  Y' is the only component
!    displayed on black and white televisions.  The I and Q components
!    contain color information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB color coordinates to be
!    converted.
!
!    Input, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
!    Output, real ( kind = 8 ) YPRIME, I, Q, the corresponding Y'IQ color
!    coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) i
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr

  yprime = yr * r + yg * g + yb * b

  i = 0.7357D+00 * ( r - yprime ) - 0.2684D+00 * ( b - yprime )
  q = 0.4777D+00 * ( r - yprime ) + 0.4133D+00 * ( b - yprime )

  return
end
subroutine rgb_to_ypbpr ( r, g, b, yprime, pb, pr )

!*****************************************************************************80
!
!! RGB_TO_YPBPR converts RGB to Y'PbPr color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB color coordinates to be
!    converted.
!
!    Output, real ( kind = 8 ) YPRIME, PB, PR, the corresponding Y'PbPr
!    color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) pb
  real ( kind = 8 ) pr
  real ( kind = 8 ) r
  real ( kind = 8 ) yprime

  yprime =   0.2122 * r + 0.7013 * g + 0.0865 * b
  pb     =  -0.1162 * r - 0.3838 * g + 0.5000 * b
  pr     =   0.5000 * r - 0.4451 * g - 0.0549 * b

  return
end
subroutine rgb_to_yuv ( r, g, b, yr, yg, yb, yprime, u, v )

!*****************************************************************************80
!
!! RGB_TO_YUV converts RGB to Y'UV color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!    Y'UV colors are used in European PAL commercial color television
!    broadcasting.  The Y' component measures luma, an approximation
!    to the luminance, or amount of light.  Y' is the only component
!    displayed on black and white televisions.  The U and V components
!    contain color information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB color coordinates to be
!    converted.
!
!    Input, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
!    Output, real ( kind = 8 ) YPRIME, U, V, the corresponding Y'UV color
!    coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) r
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr

  yprime = yr * r + yg * g + yb * b

  u = ( b - yprime ) / 2.029D+00
  v = ( r - yprime ) / 1.140D+00

  return
end
subroutine rgbcie_to_xyz_cap ( r, g, b, xcap, ycap, zcap )

!*****************************************************************************80
!
!! RGBCIE_TO_XYZ_CAP converts CIE RGB to CIE XYZ color coordinates.
!
!  Discussion:
!
!    The CIE RGB color coordinates were based on color matching data
!    using three primaries, labeled R, G, and B, associated with
!    monochromatic light at the wavelengths of 700, 546.1, and 435.8
!    nanometers.  However, this set or primaries was found unsatisfactory
!    to use as the basis for color representation, since the
!    representation of some colors required negative coefficients.
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jonas Gomes and Luiz Velho,
!    Image Processing for Computer Graphics,
!    Springer, 1997.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the CIE RGB color coordinates to
!    be converted.
!
!    Output, real ( kind = 8 ) XCAP, YCAP, ZCAP, the corresponding CIE XYZ
!    color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) r
  real ( kind = 8 ) xcap
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap

  xcap = 0.489989D+00 * r + 0.310008D+00 * g + 0.20D+00    * b
  ycap = 0.176962D+00 * r + 0.81240D+00  * g + 0.01D+00    * b
  zcap = 0.0D+00      * r + 0.01D+00     * g + 0.99D+00    * b

  return
end
subroutine rgbprime_to_lcc ( rprime, gprime, bprime, yr, yg, yb, luma, &
  chroma1, chroma2 )

!*****************************************************************************80
!
!! RGBPRIME_TO_LCC converts R'G'B' to LCC color coordinates.
!
!  Discussion:
!
!    The R'G'B' color system is a nonlinear video signal measurement.
!    Each coordinate must be between 0 and 1.
!
!    The LCC color coordinate system records a color as three components,
!    Luma, Chroma1 and Chroma2.  The LCC color coordinates are used
!    in an intermediate calculation of the PhotoYCC color coordinates.
!    Luma is scaled luminance, Chroma1 is (B'-Luma) and Chroma2 is (R'-Luma).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RPRIME, GPRIME, BPRIME, the R'G'B' color
!    coordinates.
!
!    Input, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
!    Output, real ( kind = 8 ) LUMA, CHROMA1, CHROMA2, the corresponding
!    LCC color coordinates.
!
  implicit none

  real ( kind = 8 ) chroma1
  real ( kind = 8 ) chroma2
  real ( kind = 8 ) bprime
  real ( kind = 8 ) gprime
  real ( kind = 8 ) luma
  real ( kind = 8 ) rprime
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr

  yprime = yr * rprime + yg * gprime + yb * bprime

  luma = yprime
  chroma1 = bprime - luma
  chroma2 = rprime - luma

  return
end
subroutine rgbprime_to_rgb ( rprime, gprime, bprime, r, g, b )

!*****************************************************************************80
!
!! RGBPRIME_TO_RGB converts R'G'B' to RGB color coordinates.
!
!  Discussion:
!
!    The R'G'B' color system is a nonlinear video signal measurement.
!    Each coordinate must be between 0 and 1.
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RPRIME, GPRIME, BPRIME, the R'G'B' color
!    coordinates.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) bprime
  real ( kind = 8 ) g
  real ( kind = 8 ) gprime
  real ( kind = 8 ) r
  real ( kind = 8 ) rprime

  call nonlin_to_lin ( rprime, r )
  call nonlin_to_lin ( gprime, g )
  call nonlin_to_lin ( bprime, b )

  return
end
subroutine s_c_delete ( s, c )

!*****************************************************************************80
!
!! S_C_DELETE removes all occurrences of a character from a string.
!
!  Discussion:
!
!    Each time the given character is found in the string, the characters
!    to the right of the string are shifted over one position.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
!    Input, character C, the character to be removed.
!
  implicit none

  character c
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) nchar
  character ( len = * ) s

  nchar = len_trim ( s )

  iput = 1

  do iget = 1, nchar

    if ( s(iget:iget) == c ) then

    else if ( iput == iget ) then
      iput = iput + 1
    else
      s(iput:iput) = s(iget:iget)
      iput = iput + 1
    end if

  end do

  s(iput:nchar) = ' '

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine srgb_to_xyz_cap ( sr, sg, sb, xcap, ycap, zcap )

!*****************************************************************************80
!
!! SRGB_TO_XYZ_CAP converts sRGB to CIE XYZ color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!    The sRGB color space is based on the monitor characteristics expected
!    in a dimly lit office.  Each of the coordinates (SR,SG,SB) is
!    an integer between 0 and 255.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    International Electrotechnical Commission,
!    Standard IEC 61966-2-1
!    http://www.srgb.com
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SR, SG, SB, the sRGB color coordinates.
!
!    Output, real ( kind = 8 ) XCAP, YCAP, ZCAP, the corresponding CIE XYZ
!    color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b_prime
  real ( kind = 8 ) g
  real ( kind = 8 ) g_prime
  real ( kind = 8 ), parameter :: power = 2.4D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r_prime
  integer ( kind = 4 ) sb
  integer ( kind = 4 ) sg
  integer ( kind = 4 ) sr
  real ( kind = 8 ) xcap
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap

  r_prime = real ( sr, kind = 8 ) / 255.0D+00
  g_prime = real ( sg, kind = 8 ) / 255.0D+00
  b_prime = real ( sb, kind = 8 ) / 255.0D+00

  if ( r_prime <= 0.04045D+00 ) then
    r = r_prime / 12.92D+00
  else
    r = ( ( r_prime + 0.055D+00 ) / 1.055D+00 )**power
  end if

  if ( g_prime <= 0.04045D+00 ) then
    g = g_prime / 12.92D+00
  else
    g = ( ( g_prime + 0.055D+00 ) / 1.055D+00 )**power
  end if

  if ( b_prime <= 0.04045D+00 ) then
    b = b_prime / 12.92D+00
  else
    b = ( ( b_prime + 0.055D+00 ) / 1.055D+00 )**power
  end if

  xcap = 0.4124D+00 * r + 0.3576D+00 * g + 0.1805D+00 * b
  ycap = 0.2126D+00 * r + 0.7152D+00 * g + 0.0722D+00 * b
  zcap = 0.0193D+00 * r + 0.1192D+00 * g + 0.9505D+00 * b

  return
end
subroutine t_to_spd ( t, lambda, power )

!*****************************************************************************80
!
!! T_TO_SPD evaluates the black body power spectrum at a given temperature.
!
!  Discussion:
!
!    Planck's law gives the spectral power distribution function of
!    radiation from a black body, per unit volume per infinitesimal
!    increment of wavelength, as:
!
!      SPD(Lambda,T)
!        = 1/Volume * dPower / dLambda
!        = 8 * Pi * H * C /
!          ( (10**(-9))**4 * Lambda**5 * ( EXP ( P ) - 1 ) )
!
!    where
!
!      P = H * C / ( ( Lambda * 10**(-9) ) * K * T );
!      Lambda = Wavelength, in nanometers;
!      T = Temperature of the black body, in degrees Kelvin;
!      H = Planck's constant, in joule-seconds;
!      C = Speed of light, in meters/second;
!      K = Boltzmann's constant, in joules / degrees Kelvin;
!      Volume = Volume of the cavity, in cubic meters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin,
!    of the black body.
!
!    Input, real ( kind = 8 ) LAMBDA, the wavelength, in nanometers, at
!    which the spectral power distribution function (SPD) is to be evaluated.
!
!    Output, real ( kind = 8 ) POWER, the black body spectral power distribution
!    function at the given temperature and wavelength, per volume,
!    per wavelength.  The units are ( 1 / Meter**3 ) * Joules / Nanometer,
!    that is, 1/Volume * ( Energy Increment / Wavelength Increment ).
!
  implicit none

  real ( kind = 8 ), parameter :: c = 2.9979246D+08
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) expon
  real ( kind = 8 ), parameter :: h = 6.626176D-34
  real ( kind = 8 ), parameter :: k = 1.38066D-23
  real ( kind = 8 ) lambda
  real ( kind = 8 ), parameter :: nmtom = 1.0D-09
  real ( kind = 8 ) power
  real ( kind = 8 ) t

  expon = h * c / ( nmtom * lambda * k * t )

  power = 8.0D+00 * pi * h * c / &
    ( nmtom**4 * lambda**5 * ( exp ( expon ) - 1.0D+00 ) )

  return
end
subroutine t_to_xyz ( t, x, y, z )

!*****************************************************************************80
!
!! T_TO_XY returns CIE xyz color coordinates for black body radiation.
!
!  Discussion:
!
!    The CIE xyz system defines a color in terms of its normalized
!    color coordinates (x,y,z), without reference to the absolute strength
!    or luminance of the color.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Guenter Wyszecki and W Stiles,
!    Color Science, Concepts and Methods, Quantitative Data and Formulas,
!    John Wiley, 1967, page 48.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin, of the
!    black body.  Data is only available for T in the range of 1,000 to
!    30,000 degrees.  Input values of T outside this range will result in
!    output values of X and Y at the nearest endpoint of the data range.
!
!    Output, real ( kind = 8 ) X, Y, Z, the CIE xyz color coordinates of
!    the black body radiation.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndat = 53

  real ( kind = 8 ) t
  real ( kind = 8 ), save, dimension ( ndat ) :: tdat = (/ &
     1000.0D+00,  1200.0D+00,  1400.0D+00,  1500.0D+00,  1600.0D+00, &
     1700.0D+00,  1800.0D+00,  1900.0D+00,  2000.0D+00,  2100.0D+00, &
     2200.0D+00,  2300.0D+00,  2400.0D+00,  2500.0D+00,  2600.0D+00, &
     2700.0D+00,  2800.0D+00,  2900.0D+00,  3000.0D+00,  3100.0D+00, &
     3200.0D+00,  3300.0D+00,  3400.0D+00,  3500.0D+00,  3600.0D+00, &
     3700.0D+00,  3800.0D+00,  3900.0D+00,  4000.0D+00,  4100.0D+00, &
     4200.0D+00,  4300.0D+00,  4400.0D+00,  4500.0D+00,  4600.0D+00, &
     4700.0D+00,  4800.0D+00,  4900.0D+00,  5000.0D+00,  5200.0D+00, &
     5400.0D+00,  5600.0D+00,  5800.0D+00,  6000.0D+00,  6500.0D+00, &
     7000.0D+00,  7500.0D+00,  8000.0D+00,  8500.0D+00,  9000.0D+00, &
    10000.0D+00, 15000.0D+00, 30000.0D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( ndat ) :: xdat = (/ &
    0.6526, 0.6249, 0.5984, 0.5856, 0.5731, &
    0.5610, 0.5491, 0.5377, 0.5266, 0.5158, &
    0.5055, 0.4956, 0.4860, 0.4769, 0.4681, &
    0.4597, 0.4517, 0.4441, 0.4368, 0.4299, &
    0.4233, 0.4170, 0.4109, 0.4052, 0.3997, &
    0.3945, 0.3896, 0.3848, 0.3804, 0.3760, &
    0.3719, 0.3680, 0.3643, 0.3607, 0.3573, &
    0.3540, 0.3509, 0.3479, 0.3450, 0.3397, &
    0.3347, 0.3301, 0.3259, 0.3220, 0.3135, &
    0.3063, 0.3003, 0.2952, 0.2908, 0.2869, &
    0.2806, 0.2637, 0.2501 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), save, dimension ( ndat ) :: ydat = (/ &
    0.3446, 0.3676, 0.3859, 0.3932, 0.3993, &
    0.4043, 0.4083, 0.4112, 0.4133, 0.4146, &
    0.4152, 0.4152, 0.4147, 0.4137, 0.4123, &
    0.4106, 0.4086, 0.4064, 0.4041, 0.4015, &
    0.3989, 0.3962, 0.3935, 0.3907, 0.3879, &
    0.3851, 0.3822, 0.3795, 0.3767, 0.3740, &
    0.3713, 0.3687, 0.3660, 0.3635, 0.3610, &
    0.3586, 0.3562, 0.3539, 0.3516, 0.3472, &
    0.3430, 0.3391, 0.3353, 0.3318, 0.3236, &
    0.3165, 0.3103, 0.3048, 0.2999, 0.2956, &
    0.2883, 0.2673, 0.2489 /)
  real ( kind = 8 ) z

  if ( t < tdat(1) ) then
    x = xdat(1)
    y = ydat(1)
  else if ( tdat(ndat) < t ) then
    x = xdat(ndat)
    y = ydat(ndat)
  else
    call interp ( ndat, t, tdat, x, xdat )
    call interp ( ndat, t, tdat, y, ydat )
  end if

  z = 1.0D+00 - x - y

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
subroutine uvprime_to_xyz ( uprime, vprime, x, y, z )

!*****************************************************************************80
!
!! UVPRIME_TO_XYZ converts CIE u'v' to CIE xyz color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) UPRIME, VPRIME, the CIE u'v' color coordinates.
!
!    Output, real ( kind = 8 ) X, Y, Z, the CIE xyz color coordinates.
!
  implicit none

  real ( kind = 8 ) denom
  real ( kind = 8 ) uprime
  real ( kind = 8 ) vprime
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  denom = 6.0D+00 * uprime - 16.0D+00 * vprime + 12.0D+00

  if ( denom == 0.0D+00 ) then

    x = 0.0D+00
    y = 0.0D+00
    z = 0.0D+00

  else

    x = 9.0D+00 * uprime / denom
    y = 4.0D+00 * vprime / denom
    z = ( - 3.0D+00 * uprime - 20.0D+00 * vprime + 12.0D+00 ) / denom

  end if

  return
end
subroutine uvprimey_to_xyz_cap ( uprime, vprime, xcap, ycap, zcap )

!*****************************************************************************80
!
!! UVPRIMEY_TO_XYZ_CAP converts CIE u'v'Y to CIE XYZ color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) UPRIME, VPRIME, YCAP, the CIE u'v'Y color
!    coordinates.
!
!    Output, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!
  implicit none

  real ( kind = 8 ) uprime
  real ( kind = 8 ) vprime
  real ( kind = 8 ) xcap
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap

  xcap = 9.0D+00 * ycap * uprime / ( 4.0D+00 * vprime )
  zcap = - xcap / 3.0D+00 - 5.0D+00 * ycap + 3.0D+00 * ycap / vprime

  return
end
subroutine xy_to_uvwprime ( x, y, uprime, vprime, wprime )

!*****************************************************************************80
!
!! XY_TO_UVWPRIME converts CIE xy to CIE u'v'w' color coordinates.
!
!  Discussion:
!
!    The CIE xy system defines a color in terms of its normalized
!    color coordinates (x,y), without reference to the absolute strength
!    or luminance of the color.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the CIE xy color coordinates.
!
!    Output, real ( kind = 8 ) UPRIME, VPRIME, WPRIME, the CIE u'v'w'
!    color coordinates.
!
  implicit none

  real ( kind = 8 ) denom
  real ( kind = 8 ) uprime
  real ( kind = 8 ) vprime
  real ( kind = 8 ) wprime
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  denom = - 2.0D+00 * x + 12.0D+00 * y + 3.0D+00

  if ( denom == 0.0D+00 ) then

    uprime = 4.0D+00
    vprime = 9.0D+00
    wprime = 3.0D+00

  else

    uprime = 4.0D+00 * x / denom
    vprime = 9.0D+00 * y / denom
    wprime = ( - 6.0D+00 * x + 3.0D+00 * y + 3.0D+00 ) / denom

  end if

  return
end
subroutine xyy_check ( x, y, ycap )

!*****************************************************************************80
!
!! XYY_CHECK corrects out-of-range CIE xyY color coordinates.
!
!  Discussion:
!
!    The CIE xyY color system describes a color based on the normalized
!    color coordinates (x,y), and the unnormalized CIE luminance Y.
!    The luminance must be positive.  The values x and y must each be
!    nonnegative, and must sum to no more than 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y, YCAP, the CIE xyY color
!    coordinates to be checked.  YCAP must be positive.  X and Y must
!    be nonnegative, and must sum to no more than 1.
!
  implicit none

  real ( kind = 8 ) sum2
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) ycap

  ycap = max ( ycap, 0.0D+00 )

  x = max ( x, 0.0D+00 )
  y = max ( y, 0.0D+00 )

  sum2 = x + y
  if ( 1.0D+00 < sum2 ) then
    x = x / sum2
    y = y / sum2
  end if

  return
end
subroutine xyy_to_xyz_cap ( x, y, xcap, ycap, zcap )

!*****************************************************************************80
!
!! XYY_TO_XYZ_CAP converts CIE xyY to CIE XYZ color coordinates.
!
!  Discussion:
!
!    The CIE xyY system defines a color in terms of its normalized
!    color coordinates (x,y), plus the value of Y, which allows the
!    normalizing factor to be determined.
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) x, y, the CIE xy color coordinates.
!
!    Output, real ( kind = 8 ) XCAP, the CIE X color coordinate.
!
!    Input, real ( kind = 8 ) YCAP, the CIE Y color coordinate.
!
!    Output, real ( kind = 8 ) ZCAP, the CIE Z color coordinate.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) xcap
  real ( kind = 8 ) y
  real ( kind = 8 ) ycap
  real ( kind = 8 ) z
  real ( kind = 8 ) zcap

  if ( ycap == 0.0D+00 ) then

    xcap = 0.0D+00
    zcap = 0.0D+00

  else

    xcap = x * ycap / y
    xcap = max ( xcap, 0.0D+00 )

    z = 1.0D+00 - x - y
    zcap = z * ycap / y
    zcap = max ( zcap, 0.0D+00 )

  end if

  return
end
subroutine xyz_cap_to_lab ( xcap, ycap, zcap, xcapn, ycapn, zcapn, lstar, &
  astar, bstar )

!*****************************************************************************80
!
!! XYZ_CAP_TO_LAB converts CIE XYZ to CIE LAB color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!    The CIE LAB system describes a color based on three qualities:
!      L* is CIE lightness, similar to luminance, the amount of
!        light, but adjusted to account for human perception.
!        0 <= L* <= 100.
!      a* is the amount of red chrominance, with negative values
!        indicating green.  -500 <= a* <= 500.
!      b* is the amount of yellow chrominance, with negative values
!        indicating blue.  -200 <= b* <= 200.
!    The CIE LAB model is more suitable than the CIE LUV model
!    for situations involving reflected light.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!
!    Input, real ( kind = 8 ) XCAPN, YCAPN, ZCAPN, the CIE XYZ color
!    coordinates of white.
!
!    Output, real ( kind = 8 ) LSTAR, ASTAR, BSTAR, the CIE LAB color
!    coordinates.
!
  implicit none

  real ( kind = 8 ) astar
  real ( kind = 8 ) bstar
  real ( kind = 8 ) fx
  real ( kind = 8 ) fy
  real ( kind = 8 ) fz
  real ( kind = 8 ) lstar
  real ( kind = 8 ) r8_cubert
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcapn
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycapn
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcapn

  if ( xcapn == 0.0D+00 .or. ycapn == 0.0D+00 .or. zcapn == 0.0D+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYZ_CAP_TO_LAB - Fatal error!'
    write ( *, '(a)' ) '  XCAPN, YCAPN and ZCAPN cannot be zero.'
    write ( *, '(a,g14.6)' ) '  XCAPN = ', xcapn
    write ( *, '(a,g14.6)' ) '  YCAPN = ', ycapn
    write ( *, '(a,g14.6)' ) '  ZCAPN = ', zcapn
    stop

  end if
!
!  Compute the CIE lightness.
!
  if ( ycap <= 0.0D+00 ) then
    lstar = 0.0D+00
  else if ( ycap <= 0.008856D+00 * ycapn ) then
    lstar = 903.3D+00 * ycap / ycapn
  else if ( ycap <= ycapn ) then
    lstar = 116.0D+00 * r8_cubert ( ycap / ycapn ) - 16.0D+00
  else
    lstar = 100.0D+00
  end if

  if ( xcap <= 0.008856D+00 * xcapn ) then
    fx = 7.787D+00 * xcap / xcapn + 16.0D+00 / 116.0D+00
  else
    fx = r8_cubert ( xcap / xcapn )
  end if

  if ( ycap <= 0.008856D+00 * ycapn ) then
    fy = 7.787D+00 * ycap / ycapn + 16.0D+00 / 116.0D+00
  else
    fy = r8_cubert ( ycap / ycapn )
  end if

  if ( zcap <= 0.008856D+00 * zcapn ) then
    fz = 7.787D+00 * zcap / zcapn + 16.0D+00 / 116.0D+00
  else
    fz = r8_cubert ( zcap / zcapn )
  end if

  astar = 500.0D+00 * ( fx - fy )
  bstar = 200.0D+00 * ( fy - fz )

  return
end
subroutine xyz_cap_to_luv ( xcap, ycap, zcap, xcapn, ycapn, zcapn, lstar, &
  ustar, vstar )

!*****************************************************************************80
!
!! XYZ_CAP_TO_LUV converts CIE XYZ to CIE LUV color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!    The CIE LUV system describes a color based on three qualities:
!      L* is CIE lightness, similar to luminance, the amount of light,
!        but adjusted for human perception.  0 <= L* <= 100.
!      u* is the amount of red chrominance, with negative values
!        indicating green.
!      v* is the amount of yellow chrominance, with negative values
!        indicating blue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!
!    Input, real ( kind = 8 ) XCAPN, YCAPN, ZCAPN, the CIE XYZ color
!    coordinates of white.
!
!    Output, real ( kind = 8 ) LSTAR, USTAR, VSTAR, the CIE LUV
!    color coordinates.
!
  implicit none

  real ( kind = 8 ) lstar
  real ( kind = 8 ) r8_cubert
  real ( kind = 8 ) ustar
  real ( kind = 8 ) uprime
  real ( kind = 8 ) unprime
  real ( kind = 8 ) vstar
  real ( kind = 8 ) vprime
  real ( kind = 8 ) vnprime
  real ( kind = 8 ) wprime
  real ( kind = 8 ) wnprime
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcapn
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycapn
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcapn

  if ( ycap == 0.0D+00 ) then

    lstar = 0.0D+00
    ustar = 0.0D+00
    vstar = 0.0D+00

  else
!
!  Compute LSTAR, the CIE lightness from YCAP.
!
    if ( ycap <= 0.0D+00 ) then
      lstar = 0.0D+00
    else if ( ycap <= 0.008856D+00 * ycapn ) then
      lstar = 903.3D+00 * ycap / ycapn
    else if ( ycap <= ycapn ) then
      lstar = 116.0D+00 * r8_cubert ( ycap / ycapn ) - 16.0D+00
    else
      lstar = 100.0D+00
    end if
!
!  Compute (un',vn') from (XCAPN,YCAPN,ZCAPN).
!
    call xyz_cap_to_uvwprime ( xcapn, ycapn, zcapn, unprime, vnprime, wnprime )
!
!  Compute (u',v') from (XCAP,YCAP,ZCAP).
!
    call xyz_cap_to_uvwprime ( xcap, ycap, zcap, uprime, vprime, wprime )
!
!  Compute (u*,v*) from (u',v') and (un',vn').
!
    ustar = 13.0D+00 * lstar * ( uprime - unprime )
    vstar = 13.0D+00 * lstar * ( vprime - vnprime )

  end if

  return
end
subroutine xyz_cap_to_rgb709 ( xcap, ycap, zcap, r, g, b )

!*****************************************************************************80
!
!! XYZ_CAP_TO_RGB709 converts CIE XYZ to RGB709 color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!
!    Output, real ( kind = 8 ) R, G, B, the RGB 709 color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) r
  real ( kind = 8 ) xcap
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap

  r =  3.240479D+00 * xcap - 1.53715D+00  * ycap - 0.498535D+00 * zcap
  g = -0.969256D+00 * xcap + 1.875991D+00 * ycap + 0.041556D+00 * zcap
  b =  0.055648D+00 * xcap - 0.204043D+00 * ycap + 1.057311D+00 * zcap

  return
end
subroutine xyz_cap_to_rgbcie ( xcap, ycap, zcap, r, g, b )

!*****************************************************************************80
!
!! XYZ_CAP_TO_RGBCIE converts CIE XYZ to CIE RGB color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!    The CIE RGB color coordinates were based on color matching data
!    using three primaries, labeled R, G, and B, associated with
!    monochromatic light at the wavelengths of 700, 546.1, and 435.8
!    nanometers.  However, this set or primaries was found unsatisfactory
!    to use as the basis for color representation, since the
!    representation of some colors required negative coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jonas Gomes and Luiz Velho,
!    Image Processing for Computer Graphics,
!    Springer, 1997.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!
!    Input, real ( kind = 8 ) R, G, B, the CIE RGB color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) r
  real ( kind = 8 ) xcap
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap

  r =   2.3647D+00   * xcap - 0.89658D+00  * ycap - 0.468083D+00 * zcap
  g = - 0.515155D+00 * xcap + 1.426409D+00 * ycap + 0.088746D+00 * zcap
  b =   0.005203D+00 * xcap - 0.014407D+00 * ycap + 1.0092D+00   * zcap

  return
end
subroutine xyz_cap_to_srgb ( xcap, ycap, zcap, sr, sg, sb )

!*****************************************************************************80
!
!! XYZ_CAP_TO_SRGB converts CIE XYZ to sRGB color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!    The sRGB color space is based on the monitor characteristics expected
!    in a dimly lit office.  Each of the coordinates (SR,SG,SB) is an integer
!    between 0 and 255.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    International Electrotechnical Commission,
!    Standard IEC 61966-2-1
!    http://www.srgb.com
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!
!    Output, integer ( kind = 4 ) SR, SG, SB, the sRGB color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b_prime
  real ( kind = 8 ) g
  real ( kind = 8 ) g_prime
  real ( kind = 8 ), parameter :: rewop = 1.0D+00 / 2.4D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r_prime
  integer ( kind = 4 ) sb
  integer ( kind = 4 ) sg
  integer ( kind = 4 ) sr
  real ( kind = 8 ) xcap
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap

  r =  3.2406D+00 * xcap - 1.5372D+00 * ycap - 0.4986D+00 * zcap
  g = -0.9689D+00 * xcap + 1.8758D+00 * ycap + 0.0415D+00 * zcap
  b =  0.0557D+00 * xcap - 0.2040D+00 * ycap + 1.0570D+00 * zcap

  if ( r <= 0.0031308D+00 ) then
    r_prime = 12.92D+00 * r
  else
    r_prime = 1.055D+00 * r**rewop - 0.055D+00
  end if

  if ( g <= 0.0031308D+00 ) then
    g_prime = 12.92D+00 * g
  else
    g_prime = 1.055D+00 * g**rewop - 0.055D+00
  end if

  if ( b <= 0.0031308D+00 ) then
    b_prime = 12.92D+00 * b
  else
    b_prime = 1.055D+00 * b**rewop - 0.055D+00
  end if

  sr = nint ( 255.0D+00 * r_prime )
  sr = max ( sr, 0 )
  sr = min ( sr, 255 )

  sg = nint ( 255.0D+00 * g_prime )
  sg = max ( sg, 0 )
  sg = min ( sg, 255 )

  sb = nint ( 255.0D+00 * b_prime )
  sb = max ( sb, 0 )
  sb = min ( sb, 255 )

  return
end
subroutine xyz_cap_to_uvwprime ( xcap, ycap, zcap, uprime, vprime, wprime )

!*****************************************************************************80
!
!! XYZ_CAP_TO_UVWPRIME converts CIE XYZ to CIE u'v'w' color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!
!    Output, real ( kind = 8 ) UPRIME, VPRIME, WPRIME, the CIE u'v'w'
!    color coordinates.
!
  implicit none

  real ( kind = 8 ) denom
  real ( kind = 8 ) uprime
  real ( kind = 8 ) vprime
  real ( kind = 8 ) wprime
  real ( kind = 8 ) xcap
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap

  denom = xcap + 15.0D+00 * ycap + 3.0D+00 * zcap

  if ( denom == 0.0D+00 ) then

    uprime = 4.0D+00
    vprime = 9.0D+00
    wprime = 3.0D+00

  else

    uprime = 4.0D+00 * xcap / denom
    vprime = 9.0D+00 * ycap / denom
    wprime = ( - 3.0D+00 * xcap + 6.0D+00 * ycap + 3.0D+00 * zcap ) / denom

  end if

  return
end
subroutine xyz_cap_to_xyy ( xcap, ycap, zcap, x, y )

!*****************************************************************************80
!
!! XYZ_CAP_TO_XYY converts CIE XYZ to CIE xyY color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!    The CIE xyY system defines a color in terms of its normalized
!    chromaticities (x,y), plus the value of Y, which allows the
!    normalizing factor to be determined.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!
!    Output, real ( kind = 8 ) x, y, the corresponding CIE xy chromaticities.
!
  implicit none

  real ( kind = 8 ) sum2
  real ( kind = 8 ) x
  real ( kind = 8 ) xcap
  real ( kind = 8 ) y
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap

  sum2 = xcap + ycap + zcap

  if ( sum2 == 0.0D+00 ) then
    x = 0.0D+00
    y = 0.0D+00
  else
    x = xcap / sum2
    y = ycap / sum2
  end if

  return
end
subroutine xyz_cap_to_ycc ( xcap, ycap, zcap, yr, yg, yb, yprime, c1, c2 )

!*****************************************************************************80
!
!! XYZ_CAP_TO_YCC converts CIE XYZ to PhotoYCC color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!    The Kodak PhotoYCC Color Interchange Space was developed for the
!    Photo CD System.  The Y' coordinate is a measure of luminance,
!    while C1 and C2 measure color difference chrominance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Giorgianni and Thomas Madden,
!    Digital Color Management, Encoding Solutions,
!    Addison Wesley, 1998.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!
!    Input, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
!    Output, real ( kind = 8 ) YPRIME, C1, C2, the PhotoYCC color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) bprime
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) chroma1
  real ( kind = 8 ) chroma2
  real ( kind = 8 ) g
  real ( kind = 8 ) gprime
  real ( kind = 8 ) luma
  real ( kind = 8 ) r
  real ( kind = 8 ) rprime
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcap2
  real ( kind = 8 ) yb
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycap2
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcap2
!
!  Divide by 100.0D+00
!
  xcap2 = xcap / 100.0D+00
  ycap2 = ycap / 100.0D+00
  zcap2 = zcap / 100.0D+00
!
!  Convert to RGB709 values.
!
  call xyz_cap_to_rgb709 ( xcap2, ycap2, zcap2, r, g, b )
!
!  Convert to nonlinear values.
!
  call lin_to_nonlin ( r, rprime )
  call lin_to_nonlin ( g, gprime )
  call lin_to_nonlin ( b, bprime )
!
!  Compute Luma, Chroma1, Chroma2.
!
  call rgbprime_to_lcc ( rprime, gprime, bprime, yr, yg, yb, luma, chroma1, &
    chroma2 )
!
!  Now compute Y', C1, C2:
!
  call lcc_to_ycc ( luma, chroma1, chroma2, yprime, c1, c2 )

  return
end
subroutine xyzcap_check ( xcap, ycap, zcap )

!*****************************************************************************80
!
!! XYZCAP_CHECK corrects out-of-range CIE XYZ color coordinates.
!
!  Discussion:
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ
!    color coordinates to be checked.
!
  implicit none

  real ( kind = 8 ) xcap
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap

  xcap = max ( xcap, 0.0D+00 )
  ycap = max ( ycap, 0.0D+00 )
  zcap = max ( zcap, 0.0D+00 )

  return
end
subroutine ycbcr_check ( yprime, cb, cr )

!*****************************************************************************80
!
!! YCBCR_CHECK corrects out-of-range Y'CbCr color coordinates.
!
!  Discussion:
!
!    The Y'CbCr color system is used in digital television signals.
!    The Y component measures luma, an approximation to the luminance
!    or amount of light.  Y' is the only component displayed on black
!    and white televisions.  The Cb and Cr components contain measures
!    of the blue and red components of the color.  Y' should be between
!    0 and 1, while Cb and Cr should be between -0.5 and 0.5.
!
!  Example:
!
!    Black   ( 0.000   0.000   0.000 )
!    Blue    ( 0.114   0.500   0.081 )
!    Cyan    ( 0.701   0.169  -0.337 )
!    Green   ( 0.587  -0.331  -0.419 )
!    Magenta ( 0.413   0.669   0.581 )
!    Red     ( 0.299   0.169   0.500 )
!    White   ( 1.000   0.337   0.163 )
!    Yellow  ( 0.886  -0.163   0.081 )
!
!  Method:
!
!    The Y'CbCr coordinate space is "slanted", so it is not possible
!    to give simple bounds directly on the Y, Br and Cr values.
!    The routine converts the YBrCr coordinates to RGB coordinates,
!    forces each component to be between 0 and 1, and then converts
!    back to Y'CbCr.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) YPRIME, CB, CR, the Y'CbCr color
!    coordinates to be checked.
!
  implicit none

  real ( kind = 8 ) cb
  real ( kind = 8 ) cr
  real ( kind = 8 ) yprime

  yprime = max ( yprime, 0.0D+00 )
  yprime = min ( yprime, 1.0D+00 )

  cb = max ( cb, -0.5D+00 )
  cb = min ( cb, +0.5D+00 )

  cr = max ( cr, -0.5D+00 )
  cr = min ( cr, +0.5D+00 )

  return
end
subroutine ycbcr_to_lcc ( yprime, cb, cr, luma, chroma1, chroma2 )

!*****************************************************************************80
!
!! YCBCR_TO_LCC converts Y'CbCr to LCC color coordinates.
!
!  Discussion:
!
!    The Y'CbCr color system is used in digital television signals.
!    The Y' component measures luma, an approximation to the luminance
!    or amount of light.  Y' is the only component displayed on black
!    and white televisions.  The Cb and Cr components contain measures
!    of the blue and red components of the color.  Y' should be between
!    0 and 1, with reference black at 16/255 and reference white at 235/255,
!    while Cb and Cr should be between -0.5 and 0.5.
!
!    The LCC color coordinate system records a color as three components,
!    Luma, Chroma1 and Chroma2.  The LCC color coordinates are used
!    in an intermediate calculation of the PhotoYCC color coordinates.
!    Luma is scaled luminance, Chroma1 is (B'-Luma) and Chroma2 is (R'-Luma).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    C Wayne Brown and Barry Shepherd,
!    Graphics File Formats,
!    Manning Publications, 1995.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) YPRIME, CB, CR, the Y'CbCr color coordinates.
!
!    Output, real ( kind = 8 ) LUMA, CHROMA1, CHROMA2, the LCC color
!    coordinates.
!
  implicit none

  real ( kind = 8 ) cb
  real ( kind = 8 ) chroma1
  real ( kind = 8 ) chroma2
  real ( kind = 8 ) cr
  real ( kind = 8 ) luma
  real ( kind = 8 ) yprime

  luma = ( 255.0D+00 * yprime - 16.0D+00 ) / 219.0D+00
  chroma1 = ( 255.0D+00 * cb - 128.0D+00 ) / ( 224.0D+00 * 0.564D+00 )
  chroma2 = ( 255.0D+00 * cr - 128.0D+00 ) / ( 224.0D+00 * 0.713D+00 )

  return
end
subroutine ycbcr_to_rgb ( yprime, cb, cr, yr, yg, yb, r, g, b )

!*****************************************************************************80
!
!! YCBCR_TO_RGB converts Y'CbCr to RGB color coordinates.
!
!  Discussion:
!
!    The Y'CbCr color system is used in digital television signals.
!    The Y' component measures luma, an approximation to the luminance
!    or amount of light.  Y' is the only component displayed on black
!    and white televisions.  The Cb and Cr components contain measures
!    of the blue and red components of the color.  Y' should be between
!    0 and 1, while Cb and Cr should be between -0.5 and 0.5.
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) YPRIME, CB, CR, the Y'CbCr color coordinates
!    to be converted.
!
!    Input, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) cb
  real ( kind = 8 ) cr
  real ( kind = 8 ) g
  real ( kind = 8 ) r
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr

  b = yprime + 2.0D+00 * ( 1.0D+00 - yb ) * cb
  r = yprime + 2.0D+00 * ( 1.0D+00 - yr ) * cr
  g = ( yprime - yr * r - yb * b ) / yg

  return
end
subroutine ycbcr_to_ycc ( yprime, cb, cr, yprime2, c1, c2 )

!*****************************************************************************80
!
!! YCBCR_TO_YCC converts Y'CbCr to PhotoYCC color coordinates.
!
!  Discussion:
!
!    The Y'CbCr color system is used in digital television signals.
!    The Y' component measures luma, an approximation to the luminance
!    or amount of light.  Y' is the only component displayed on black
!    and white televisions.  The Cb and Cr components contain measures
!    of the blue and red components of the color.  Y' should be between
!    0 and 1, while Cb and Cr should be between -0.5 and 0.5.
!
!    The Kodak PhotoYCC Color Interchange Space was developed for the
!    Photo CD System.  The Y' coordinate is a measure of luminance,
!    while C1 and C2 measure color difference chrominance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Giorgianni and Thomas Madden,
!    Digital Color Management, Encoding Solutions,
!    Addison Wesley, 1998.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) YPRIME, CB, CR, the Y'CbCr color coordinates.
!
!    Input, real ( kind = 8 ) YPRIME2, C1, C2, the PhotoYCC color coordinates.
!
  implicit none

  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) cb
  real ( kind = 8 ) cr
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yprime2

  yprime2 = yprime / 1.402D+00
  c1 = ( cb + 73.400D+00 ) / 1.291D+00
  c2 = ( cr + 55.638D+00 ) / 1.340D+00

  return
end
subroutine ycc_test ( itest, ytest, c1test, c2test )

!*****************************************************************************80
!
!! YCC_TEST supplies PhotoYCC values for tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ITEST, the index of the test.
!
!    Output, real ( kind = 8 ) YTEST, C1TEST, C2TEST, sample PhotoYCC color
!    coordinates for testing.  If ITEST is outside the range of
!    data, then YTEST, C1TEST and C2TEST are returned as -1.
!
  implicit none

  real ( kind = 8 ) c1test
  real ( kind = 8 ) c2test
  integer ( kind = 4 ) itest
  real ( kind = 8 ) ytest

  if ( itest == 1 ) then

    ytest = 3.0D+00
    c1test = 250.0D+00
    c2test = 200.0D+00

  else if ( itest == 2 ) then

    ytest = 10.0D+00
    c1test = 200.0D+00
    c2test = 20.0D+00

  else if ( itest == 3 ) then

    ytest = 50.0D+00
    c1test = 75.0D+00
    c2test = 0.0D+00

  else if ( itest == 4 ) then

    ytest = 100.0D+00
    c1test = 30.0D+00
    c2test = 120.0D+00

  else if ( itest == 5 ) then

    ytest = 150.0D+00
    c1test = 80.0D+00
    c2test = 200.0D+00

  else

    ytest = -1.0D+00
    c1test = -1.0D+00
    c2test = -1.0D+00

  end if

  return
end
subroutine ycc_to_lcc ( yprime, c1, c2, luma, chroma1, chroma2 )

!*****************************************************************************80
!
!! YCC_TO_LCC converts PhotoYCC to LCC color coordinates.
!
!  Discussion:
!
!    The Kodak PhotoYCC Color Interchange Space was developed for the
!    Photo CD System.  The Y' coordinate is a measure of luminance,
!    while C1 and C2 measure color difference chrominance.
!
!    The LCC color coordinate system records a color as three components,
!    Luma, Chroma1 and Chroma2.  The LCC color coordinates are used
!    in an intermediate calculation of the PhotoYCC color coordinates.
!    Luma is scaled luminance, Chroma1 is (B'-Luma) and Chroma2 is (R'-Luma).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Giorgianni and Thomas Madden,
!    Digital Color Management, Encoding Solutions,
!    Addison Wesley, 1998.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) YPRIME, C1, C2, the PhotoYCC color coordinates.
!
!    Output, real ( kind = 8 ) LUMA, CHROMA1, CHROMA2, the LCC color
!    coordinates.
!
  implicit none

  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) chroma1
  real ( kind = 8 ) chroma2
  real ( kind = 8 ) luma
  real ( kind = 8 ) yprime

  luma = 1.402D+00 * yprime / 255.0D+00
  chroma1 = ( c1 - 156.0D+00 ) / 111.40D+00
  chroma2 = ( c2 - 137.0D+00 ) / 135.64D+00

  return
end
subroutine ycc_to_xyz_cap ( yprime, c1, c2, yr, yg, yb, xcap, ycap, zcap )

!*****************************************************************************80
!
!! YCC_TO_XYZ_CAP converts PhotoYCC to CIE XYZ color coordinates.
!
!  Discussion:
!
!    The Kodak PhotoYCC Color Interchange Space was developed for the
!    Photo CD System.  The Y' coordinate is a measure of luminance,
!    while C1 and C2 measure color difference chrominance.
!
!    The CIE XYZ color system describes a color in terms of its components
!    of X, Y and Z primaries.  In ordinary circumstances, all three of
!    these components must be nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Giorgianni and Thomas Madden,
!    Digital Color Management, Encoding Solutions,
!    Addison Wesley, 1998.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) YPRIME, C1, C2, the PhotoYCC color coordinates.
!
!    Input, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
!    Output, real ( kind = 8 ) XCAP, YCAP, ZCAP, the CIE XYZ color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) bprime
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) chroma1
  real ( kind = 8 ) chroma2
  real ( kind = 8 ) g
  real ( kind = 8 ) gprime
  real ( kind = 8 ) luma
  real ( kind = 8 ) r
  real ( kind = 8 ) rprime
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcap2
  real ( kind = 8 ) yb
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycap2
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcap2
!
!  Compute LUMA, CHROMA1, CHROMA2:
!
  call ycc_to_lcc ( yprime, c1, c2, luma, chroma1, chroma2 )
!
!  Compute R'G'B'.
!
  call lcc_to_rgbprime ( luma, chroma1, chroma2, yr, yg, yb, rprime, gprime, &
    bprime )
!
!  Convert to linear values.
!
  call nonlin_to_lin ( rprime, r )
  call nonlin_to_lin ( gprime, g )
  call nonlin_to_lin ( bprime, b )
!
!  Convert to CIE XYZ values.
!
  call rgb709_to_xyz_cap ( r, g, b, xcap2, ycap2, zcap2 )
!
!  Multiply by 100.0D+00
!
  xcap = xcap2 * 100.0D+00
  ycap = ycap2 * 100.0D+00
  zcap = zcap2 * 100.0D+00

  return
end
subroutine ycc_to_ycbcr ( yprime, c1, c2, yprime2, cb, cr )

!*****************************************************************************80
!
!! YCC_TO_YCBCR converts PhotoYCC to Y'CbCr color coordinates.
!
!  Discussion:
!
!    The Kodak PhotoYCC Color Interchange Space was developed for the
!    Photo CD System.  The Y' coordinate is a measure of luminance,
!    while C1 and C2 measure color difference chrominance.
!
!    The Y'CbCr color system is used in digital television signals.
!    The Y' component measures luma, an approximation to the luminance
!    or amount of light.  Y' is the only component displayed on black
!    and white televisions.  The Cb and Cr components contain measures
!    of the blue and red components of the color.  Y' should be between
!    0 and 1, while Cb and Cr should be between -0.5 and 0.5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Giorgianni and Thomas Madden,
!    Digital Color Management, Encoding Solutions,
!    Addison Wesley, 1998.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) YPRIME, C1, C2, the PhotoYCC color coordinates.
!
!    Output, real ( kind = 8 ) YPRIME2, CB, CR, the Y'CbCr color coordinates.
!
  implicit none

  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) cb
  real ( kind = 8 ) cr
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yprime2

  yprime2 = 1.402D+00 * yprime
  cb = 1.291D+00 * c1 - 73.400D+00
  cr = 1.340D+00 * c2 - 55.638D+00

  return
end
subroutine yiq_check ( yprime, i, q, yr, yg, yb )

!*****************************************************************************80
!
!! YIQ_CHECK corrects out-of-range Y'IQ color coordinates.
!
!  Discussion:
!
!    Y'IQ colors are used in American NTSC commercial color television
!    broadcasting.  The Y' component measures luma, an approximation to
!    luminance, or the amount of light.  Y' is the only component
!    displayed on black and white televisions.  The I and Q components
!    contain color information.
!
!  Example:
!
!    Black   ( 0.000,  0.000,  0.000 )
!    Blue    ( 0.114, -0.322,  0.311 )
!    Cyan    ( 0.701, -0.596, -0.212 )
!    Green   ( 0.587, -0.274, -0.523 )
!    Magenta ( 0.413,  0.274,  0.523 )
!    Red     ( 0.299,  0.596,  0.212 )
!    White   ( 1.000,  0.000,  0.000 )
!    Yellow  ( 0.886,  0.322, -0.311 )
!
!  Method:
!
!    The Y'IQ coordinate space is "slanted", so it is not possible
!    to give simple bounds directly on the Y', I and Q values.
!    The routine converts the Y'IQ coordinates to RGB coordinates,
!    forces each component to be between 0 and 1, and then converts
!    back to Y'IQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) YPRIME, I, Q, the Y'IQ color
!    coordinates to be checked.
!
!    Input, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) i
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr
!
!  Convert Y'IQ to RGB.
!
  call yiq_to_rgb ( yprime, i, q, yr, yg, yb, r, g, b )
!
!  Check RGB.
!
  call rgb_check ( r, g, b )
!
!  Convert RGB back to Y'IQ.
!
  call rgb_to_yiq ( r, g, b, yr, yg, yb, yprime, i, q )

  return
end
subroutine yiq_to_rgb ( yprime, i, q, yr, yg, yb, r, g, b )

!*****************************************************************************80
!
!! YIQ_TO_RGB converts Y'IQ to RGB color coordinates.
!
!  Discussion:
!
!    Y'IQ colors are used in American NTSC commercial color television
!    broadcasting.  The Y' component measures luma, an approximation to
!    luminance, or the amount of light.  Y' is the only component
!    displayed on black and white televisions.  The I and Q components
!    contain color information.
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) YPRIME, I, Q, the Y'IQ color coordinates
!    to be converted.
!
!    Input, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) i
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr

  r = yprime + 0.956D+00 * i + 0.621D+00 * q
  b = yprime - 1.105D+00 * i + 1.702D+00 * q

  g = ( yprime - yr * r - yb * b ) / yg

  return
end
subroutine yiq_to_yuv ( yprime, i, q, yprime2, u, v )

!*****************************************************************************80
!
!! YIQ_TO_YUV converts Y'IQ to Y'UV color coordinates.
!
!  Discussion:
!
!    Y'IQ colors are used in American NTSC commercial color television
!    broadcasting.  The Y' component measures luma, an approximation to
!    luminance, or the amount of light.  Y' is the only component
!    displayed on black and white televisions.  The I and Q components
!    contain color information.
!
!    Y'UV colors are used in European PAL commercial color television
!    broadcasting.  The Y' component measures luma, an approximation
!    to the luminance, or amount of light.  Y' is the only component
!    displayed on black and white televisions.  The U and V components
!    contain color information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) YPRIME, I, Q, the Y'IQ color coordinates
!    to be converted.
!
!    Output, real ( kind = 8 ) YPRIME2, U, V, the corresponding Y'UV
!    color coordinates.
!
  implicit none

  real ( kind = 8 ) i
  real ( kind = 8 ) q
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yprime2

  yprime2 = yprime
  u       =    - 1.1270D+00 * i + 1.8050D+00 * q
  v       =      0.9489D+00 * i + 0.6561D+00 * q

  return
end
subroutine ypbpr_to_rgb ( yprime, pb, pr, r, g, b )

!*****************************************************************************80
!
!! YPBPR_TO_RGB converts Y'PbPr to RGB color coordinates.
!
!  Discussion:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) YPRIME, PB, PR, the Y'PbPr color coordinates
!    to be converted.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) pb
  real ( kind = 8 ) pr
  real ( kind = 8 ) r
  real ( kind = 8 ) yprime

  r = yprime + 0.0000D+00 * pb + 1.5756D+00 * pr
  g = yprime - 0.2253D+00 * pb - 0.4768D+00 * pr
  b = yprime + 1.8270D+00 * pb + 0.0000D+00 * pr

  return
end
subroutine yuv_check ( yprime, u, v, yr, yg, yb )

!*****************************************************************************80
!
!! YUV_CHECK corrects out-of-range Y'UV color coordinates.
!
!  Discussion:
!
!    Y'UV colors are used in European PAL commercial color television
!    broadcasting.  The Y' component measures luma, an approximation
!    to the luminance, or amount of light.  Y' is the only component
!    displayed on black and white televisions.  The U and V components
!    contain color information.
!
!  Example:
!
!    Black   ( 0.000   0.000   0.000 )
!    Blue    ( 0.114   0.436  -0.100 )
!    Cyan    ( 0.701   0.147  -0.615 )
!    Green   ( 0.587  -0.289  -0.515 )
!    Magenta ( 0.413   0.289   0.515 )
!    Red     ( 0.299  -0.147   0.615 )
!    White   ( 1.000   0.000   0.000 )
!    Yellow  ( 0.886  -0.436   0.100 )
!
!  Method:
!
!    The Y'UV coordinate space is "slanted", so it is not possible
!    to give simple bounds directly on the Y', U and V values.
!    The routine converts the Y'UV coordinates to RGB coordinates,
!    forces each component to be between 0 and 1, and then converts
!    back to Y'UV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) YPRIME, U, V, the Y'UV color coordinates
!    to be checked.
!
!    Input, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
  implicit none
!
  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) r
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr
!
!  Convert Y'UV to RGB.
!
  call yuv_to_rgb ( yprime, u, v, yr, yg, yb, r, g, b )
!
!  Check RGB.
!
  call rgb_check ( r, g, b )
!
!  Convert RGB back to Y'UV.
!
  call rgb_to_yuv ( r, g, b, yr, yg, yb, yprime, u, v )

  return
end
subroutine yuv_to_rgb ( yprime, u, v, yr, yg, yb, r, g, b )

!*****************************************************************************80
!
!! YUV_TO_RGB converts Y'UV to RGB color coordinates.
!
!  Discussion:
!
!    Y'UV colors are used in European PAL commercial color television
!    broadcasting.  The Y' component measures luma, an approximation
!    to the luminance, or amount of light.  Y' is the only component
!    displayed on black and white televisions.  The U and V components
!    contain color information.
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) YPRIME, U, V, the Y'UV color coordinates
!    to be converted.
!
!    Input, real ( kind = 8 ) YR, YG, YB, the coefficients of the R, G and B
!    primaries in the luminance function.
!
!    Output, real ( kind = 8 ) R, G, B, the corresponding RGB color coordinates.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) r
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr

  b = yprime + 2.029D+00 * u
  r = yprime + 1.140D+00 * v

  g = ( yprime - yr * r - yb * b ) / yg

  return
end
subroutine yuv_to_yiq ( yprime, u, v, yprime2, i, q )

!*****************************************************************************80
!
!! YUV_TO_YIQ converts Y'UV to Y'IQ color coordinates.
!
!  Discussion:
!
!    Y'UV colors are used in European PAL commercial color television
!    broadcasting.  The Y' component measures luma, an approximation
!    to the luminance, or amount of light.  Y' is the only component
!    displayed on black and white televisions.  The U and V components
!    contain color information.
!
!    Y'IQ colors are used in American NTSC commercial color television
!    broadcasting.  The Y' component measures luma, an approximation to
!    luminance, or the amount of light.  Y' is the only component
!    displayed on black and white televisions.  The I and Q components
!    contain color information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) YPRIME, U, V, the Y'UV color coordinates
!    to be converted.
!
!    Output, real ( kind = 8 ) YPRIME2, I, Q, the corresponding Y'IQ
!    color coordinates.
!
  implicit none

  real ( kind = 8 ) i
  real ( kind = 8 ) q
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yprime2

  yprime2 = yprime
  i       =     - 0.2676D+00 * u + 0.7361D+00 * v
  q       =     + 0.3869D+00 * u + 0.4596D+00 * v

  return
end
