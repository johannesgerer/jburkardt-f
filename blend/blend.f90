subroutine blend_101 ( r, x0, x1, x )

!*****************************************************************************80
!
!! BLEND_101 extends scalar endpoint data to a line.
!
!  Diagram:
!
!    0-----r-----1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the coordinate where an interpolated
!    value is desired.  
!
!    Input, real ( kind = 8 ) X0, X1, the data values at the ends of the line.
!
!    Output, real ( kind = 8 ) X, the interpolated data value at (R).
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1

  x = ( 1.0D+00 - r ) * x0 + r * x1

  return
end
subroutine blend_102 ( r, s, x00, x01, x10, x11, x )

!*****************************************************************************80
!
!! BLEND_102 extends scalar point data into a square.
!
!  Diagram:
!
!    01------------11
!     |      .      |
!     |      .      |
!     |.....rs......|
!     |      .      |
!     |      .      |
!    00------------10
!
!  Formula:
!
!    Written in terms of R and S, the map has the form:
!
!      X(R,S) =
!               1     * ( + x00 )
!             + r     * ( - x00 + x10 )
!             + s     * ( - x00       + x01 )
!             + r * s * ( + x00 - x10 - x01 + x11 )
!
!    Written in terms of the coefficients, the map has the form:
!
!      X(R,S) = x00 * ( 1 - r - s + r * s )
!             + x01 * (         s - r * s )
!             + x10 * (     r     - r * s )
!             + x11 * (             r * s )
!
!             = x00 * ( 1 - r ) * ( 1 - s )
!             + x01 * ( 1 - r ) *       s
!             + x10 *       r   * ( 1 - s )
!             + x11 *       r           s
!
!    The nonlinear term ( r * s ) has an important role:
!
!      If ( x01 + x10 - x00 - x11 ) is zero, then the input data lies in
!      a plane, and the mapping is affine.  All the interpolated data 
!      will lie on the plane defined by the four corner values.  In 
!      particular, on any line through the square, data values at 
!      intermediate points will lie between the values at the endpoints.  
!
!      If ( x01 + x10 - x00 - x11 ) is not zero, then the input data does
!      not lie in a plane, and the interpolation map is nonlinear.  On
!      any line through the square, data values at intermediate points
!      may lie above or below the data values at the endpoints.  The
!      size of the coefficient of r * s will determine how severe this
!      effect is.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the coordinates where an 
!    interpolated value is desired.  
!
!    Input, real ( kind = 8 ) X00, X01, X10, X11, the data values 
!    at the corners.
!
!    Output, real ( kind = 8 ) X, the interpolated data value at (R,S).
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) x
  real ( kind = 8 ) x00
  real ( kind = 8 ) x01
  real ( kind = 8 ) x10
  real ( kind = 8 ) x11

  x =             + x00 &
      + r *     ( - x00 + x10 ) & 
      + s *     ( - x00       + x01 ) &
      + r * s * ( + x00 - x10 - x01 + x11 )

  return
end
subroutine blend_103 ( r, s, t, x000, x001, x010, x011, x100, x101, x110, &
  x111, x )

!*****************************************************************************80
!
!! BLEND_103 extends scalar point data into a cube.
!
!  Diagram:
!
!    011--------------111 
!      |               |
!      |               | 
!      |               |
!      |               |
!      |               |
!    001--------------101
!
!
!      *---------------*
!      |               |
!      |               |
!      |      rst      |
!      |               |
!      |               |
!      *---------------*
!
!
!    010--------------110
!      |               |
!      |               |
!      |               |
!      |               | 
!      |               |
!    000--------------100 
!
!
!  Formula:
!
!    Written as a polynomial in R, S and T, the interpolation map has the 
!    form:
!
!      X(R,S,T) =
!        1         * ( + x000 )
!      + r         * ( - x000 + x100 )
!      +     s     * ( - x000        + x010 )
!      +         t * ( - x000               + x001 )
!      + r * s     * ( + x000 - x100 - x010                       + x110 )
!      + r     * t * ( + x000 - x100        - x001        + x101 )
!      +     s * t * ( + x000        - x010 - x001 + x011 )
!      + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, T, the coordinates where an 
!    interpolated value is desired.
!
!    Input, real ( kind = 8 ) X000, X001, X010, X011, X100, X101, X110, 
!    X111, the data values at the corners.
!
!    Output, real ( kind = 8 ) X, the interpolated data value at (R,S,T).
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x
  real ( kind = 8 ) x000
  real ( kind = 8 ) x001
  real ( kind = 8 ) x010
  real ( kind = 8 ) x011
  real ( kind = 8 ) x100
  real ( kind = 8 ) x101
  real ( kind = 8 ) x110
  real ( kind = 8 ) x111
!
!  Interpolate the interior point.
!
  x = &
    1.0D+00     * ( + x000 ) &
    + r         * ( - x000 + x100 ) &
    +     s     * ( - x000        + x010 ) &
    +         t * ( - x000               + x001 ) &
    + r * s     * ( + x000 - x100 - x010                      + x110 ) &
    + r     * t * ( + x000 - x100        - x001        + x101 ) &
    +     s * t * ( + x000        - x010 - x001 + x011 ) &
    + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
 
  return
end
subroutine blend_112 ( r, s, x00, x01, x10, x11, xr0, xr1, x0s, x1s, x )

!*****************************************************************************80
!
!! BLEND_112 extends scalar line data into a square.
!
!  Diagram:
!
!    01-----r1-----11
!     |      .      |
!     |      .      |
!    0s.....rs.....1s
!     |      .      |
!     |      .      |
!    00-----r0-----10
!
!  Formula:
!
!    Written in terms of R and S, the interpolation map has the form:
!
!      X(R,S) =
!               1     * ( - x00       + x0s                   + xr0 )
!             + r     * (   x00       - x0s - x10       + x1s )
!             + s     * (   x00 - x01                         - xr0 + xr1 )
!             + r * s * ( - x00 + x01       + x10 - x11 )
!
!    Written in terms of the data, the map has the form:
!
!      X(R,S) =
!             - ( 1 - r ) * ( 1 - s ) * x00
!             + ( 1 - r )             * x0s
!             - ( 1 - r ) *       s   * x01
!             +             ( 1 - s ) * xr0
!             +                   s   * xr1
!             -       r   * ( 1 - s ) * x10
!             +       r               * x1s
!             -       r   *       s   * x11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the coordinates where an interpolated 
!    value is desired.  
!
!    Input, real ( kind = 8 ) X00, X01, X10, X11, the data values 
!    at the corners.
!
!    Input, real ( kind = 8 ) XR0, XR1, X0S, X1S, the data values at 
!    points along the edges corresponding to (R,0), (R,1), (0,S) and (1,S).
!
!    Output, real ( kind = 8 ) X, the interpolated data value at (R,S).
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) x
  real ( kind = 8 ) x00
  real ( kind = 8 ) x01
  real ( kind = 8 ) x10
  real ( kind = 8 ) x11
  real ( kind = 8 ) xr0
  real ( kind = 8 ) xr1
  real ( kind = 8 ) x0s
  real ( kind = 8 ) x1s

  x = - ( 1.0D+00 - r ) * ( 1.0D+00 - s ) * x00 &
      + ( 1.0D+00 - r )                   * x0s &
      - ( 1.0D+00 - r ) *             s   * x01 &
      +                   ( 1.0D+00 - s ) * xr0 &
      +                               s   * xr1 &
      -             r   * ( 1.0D+00 - s ) * x10 &
      +             r                     * x1s &
      -             r   *             s   * x11

  return
end
subroutine blend_113 ( r, s, t, x000, x001, x010, x011, x100, x101, x110, &
  x111, xr00, xr01, xr10, xr11, x0s0, x0s1, x1s0, x1s1, x00t, x01t, x10t, &
  x11t, x )

!*****************************************************************************80
!
!! BLEND_113 extends scalar line data into a cube.
!
!  Diagram:
!
!     011-----r11-----111 
!      |               |
!      |               | 
!     0s1             1s1
!      |               |
!      |               |
!     001-----r01-----101
!
!
!     01t-------------11t
!      |               |
!      |               |
!      |      rst      |
!      |               |
!      |               |
!     00t-------------10t
!
!
!     010-----r10-----110
!      |               |
!      |               |
!     0s0             1s0
!      |               | 
!      |               |
!     000-----r00-----100 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, T, the coordinates where an interpolated
!    value is desired.
!
!    Input, real ( kind = 8 ) X000, X001, X010, X011, X100, X101, X110,
!    X111, the data values at the corners.
!
!    Input, real ( kind = 8 ) XR00, XR01, XR10, XR11, X0S0, X0S1, X1S0, 
!    X1S1, X00T, X01T, X10T, X11T, the data values at points along the edges.
!
!    Output, real ( kind = 8 ) X, the interpolated data value at (R,S,T).
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x
  real ( kind = 8 ) x000
  real ( kind = 8 ) x001
  real ( kind = 8 ) x010
  real ( kind = 8 ) x011
  real ( kind = 8 ) x100
  real ( kind = 8 ) x101
  real ( kind = 8 ) x110
  real ( kind = 8 ) x111
  real ( kind = 8 ) xr00
  real ( kind = 8 ) xr01
  real ( kind = 8 ) xr0t
  real ( kind = 8 ) xr10
  real ( kind = 8 ) xr11
  real ( kind = 8 ) xr1t
  real ( kind = 8 ) xrs0
  real ( kind = 8 ) xrs1
  real ( kind = 8 ) x0s0
  real ( kind = 8 ) x0s1
  real ( kind = 8 ) x0st
  real ( kind = 8 ) x1s0
  real ( kind = 8 ) x1s1
  real ( kind = 8 ) x1st
  real ( kind = 8 ) x00t
  real ( kind = 8 ) x01t
  real ( kind = 8 ) x10t
  real ( kind = 8 ) x11t
!
!  Interpolate the points in the centers of the faces.
!
  call blend_112 ( s, t, x000, x001, x010, x011, x0s0, x0s1, x00t, x01t, x0st )
  call blend_112 ( s, t, x100, x101, x110, x111, x1s0, x1s1, x10t, x11t, x1st )
  call blend_112 ( r, t, x000, x001, x100, x101, xr00, xr01, x00t, x10t, xr0t )
  call blend_112 ( r, t, x010, x011, x110, x111, xr10, xr11, x01t, x11t, xr1t )
  call blend_112 ( r, s, x000, x010, x100, x110, xr00, xr10, x0s0, x1s0, xrs0 )
  call blend_112 ( r, s, x001, x011, x101, x111, xr01, xr11, x0s1, x1s1, xrs1 )
!
!  Interpolate the I-th coordinate component of the interior point.
!
  call blend_123 ( r, s, t, x000, x001, x010, x011, x100, x101, x110, x111, &
    xr00, xr01, xr10, xr11, x0s0, x0s1, x1s0, x1s1, x00t, x01t, x10t, x11t, &
    x0st, x1st, xr0t, xr1t, xrs0, xrs1, x )

  return
end
subroutine blend_123 ( r, s, t, x000, x001, x010, x011, x100, x101, x110, &
  x111, xr00, xr01, xr10, xr11, x0s0, x0s1, x1s0, x1s1, x00t, x01t, x10t, &
  x11t, x0st, x1st, xr0t, xr1t, xrs0, xrs1, x )
!
!*****************************************************************************80
!
!! BLEND_123 extends scalar face data into a cube.
!
!
!  Diagram:
!
!    010-----r10-----110        011-----r11-----111
!      |       .       |          |       .       |  
!      |       .       |          |       .       |
!    0s0.....rs0.....1s0        0s1.....rs1.....1s1     S
!      |       .       |          |       .       |     |
!      |       .       |          |       .       |     |
!    000-----r00-----100        001-----r01-----101     +----R
!           BOTTOM                      TOP
!
!    011-----0s1-----001        111-----1s1-----101
!      |       .       |          |       .       |  
!      |       .       |          |       .       |
!    01t.....0st.....00t        11t.....1st.....10t          T
!      |       .       |          |       .       |          |
!      |       .       |          |       .       |          |
!    010-----0s0-----000        110-----1s0-----100     S----+
!           LEFT                       RIGHT
!
!    001-----r01-----101        011-----r11-----111
!      |       .       |          |       .       |  
!      |       .       |          |       .       |
!    00t.....r0t.....100        01t.....r1t.....11t     T
!      |       .       |          |       .       |     |
!      |       .       |          |       .       |     |
!    000-----r00-----100        010-----r10-----110     +----R
!           FRONT                       BACK
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, T, the coordinates where an interpolated
!    value is desired.  
!
!    Input, real ( kind = 8 ) X000, X001, X010, X011, X100, X101, X110, 
!    X111, the data values at the corners.
!
!    Input, real ( kind = 8 ) XR00, XR01, XR10, XR11, X0S0, X0S1, X1S0, 
!    X1S1, X00T, X01T, X10T, X11T, the data values at points along the edges.
!
!    Input, real ( kind = 8 ) X0ST, X1ST, XR0T, XR1T, XRS0, XRS1, the 
!    data values at points on the faces.
!
!    Output, real ( kind = 8 ) X, the interpolated data value at (R,S,T).
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x
  real ( kind = 8 ) x000
  real ( kind = 8 ) x001
  real ( kind = 8 ) x010
  real ( kind = 8 ) x011
  real ( kind = 8 ) x100
  real ( kind = 8 ) x101
  real ( kind = 8 ) x110
  real ( kind = 8 ) x111
  real ( kind = 8 ) xr00
  real ( kind = 8 ) xr01
  real ( kind = 8 ) xr10
  real ( kind = 8 ) xr11
  real ( kind = 8 ) x0s0
  real ( kind = 8 ) x0s1
  real ( kind = 8 ) x1s0
  real ( kind = 8 ) x1s1
  real ( kind = 8 ) x00t
  real ( kind = 8 ) x01t
  real ( kind = 8 ) x10t
  real ( kind = 8 ) x11t
  real ( kind = 8 ) x0st
  real ( kind = 8 ) x1st
  real ( kind = 8 ) xr0t
  real ( kind = 8 ) xr1t
  real ( kind = 8 ) xrs0
  real ( kind = 8 ) xrs1
!
!  Interpolate the interior point.
!
  x =        ( 1.0D+00 - r ) * ( 1.0D+00 - s ) * ( 1.0D+00 - t ) * x000 &
           - ( 1.0D+00 - r ) * ( 1.0D+00 - s )                   * x00t &
           + ( 1.0D+00 - r ) * ( 1.0D+00 - s ) *             t   * x001 &
           - ( 1.0D+00 - r )                   * ( 1.0D+00 - t ) * x0s0 &
           + ( 1.0D+00 - r )                                     * x0st &
           - ( 1.0D+00 - r )                   *             t   * x0s1 &
           + ( 1.0D+00 - r ) *             s   * ( 1.0D+00 - t ) * x010 &
           - ( 1.0D+00 - r ) *             s                     * x01t &
           + ( 1.0D+00 - r ) *             s   *             t   * x011 &
           -                   ( 1.0D+00 - s ) * ( 1.0D+00 - t ) * xr00 &
           +                   ( 1.0D+00 - s )                   * xr0t &
           -                   ( 1.0D+00 - s ) *             t   * xr01 &
           +                                     ( 1.0D+00 - t ) * xrs0 &
           +                                                 t   * xrs1 &
           -                               s   * ( 1.0D+00 - t ) * xr10 &
           +                               s                     * xr1t &
           -                               s   *             t   * xr11 &
           +             r   * ( 1.0D+00 - s ) * ( 1.0D+00 - t ) * x100 &
           -             r   * ( 1.0D+00 - s )                   * x10t &
           +             r   * ( 1.0D+00 - s ) *             t   * x101 &
           -             r                     * ( 1.0D+00 - t ) * x1s0 &
           +             r                                       * x1st &
           -             r                     *             t   * x1s1 &
           +             r   *             s   * ( 1.0D+00 - t ) * x110 &
           -             r   *             s                     * x11t &
           +             r   *             s   *             t   * x111

  return
end
subroutine blend_i_0d1 ( x, m )

!*****************************************************************************80
!
!! BLEND_I_0D1 extends indexed scalar data at endpoints along a line.
!
!  Diagram:
!
!    ( X1, ..., ..., ..., ..., ..., XM )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(M).  
!    On input, X(1) and X(M) contain scalar values which are to be 
!    interpolated through the entries X(2) through X(M).  It is assumed
!    that the dependence of the data is linear in the vector index I.  
!    On output, X(2) through X(M-1) have been assigned interpolated 
!    values.
!
!    Input, integer ( kind = 4 ) M, the number of entries in X.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) x(m)

  do i = 2, m - 1

    r = real ( i - 1, kind = 8 ) &
      / real ( m - 1, kind = 8 )

    call blend_101 ( r, x(1), x(m), x(i) )

  end do

  return
end
subroutine blend_ij_0d1 ( x, m1, m2 )

!*****************************************************************************80
!
!! BLEND_IJ_0D1 extends indexed scalar data at corners into a table.
!
!  Diagram:
!
!    ( X11,  ..., ..., ..., ..., ..., X1M2  )
!    ( ...,  ..., ..., ..., ..., ..., ...   )
!    ( ...,  ..., ..., ..., ..., ..., ...   )
!    ( ...,  ..., ..., ..., ..., ..., ...   )
!    ( XM11, ..., ..., ..., ..., ..., XM1M2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(M1,M2).  
!    On input, X(1,1), X(1,M2), X(M1,1) and X(M1,M2) contain scalar 
!    values which are to be interpolated throughout the table, using 
!    the table indices I and J as independent variables. 
!    On output, all entries in X have been assigned a value.
!
!    Input, integer ( kind = 4 ) M1, M2, the number of rows and columns in X.
!
  implicit none

  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) x(m1,m2)
!
!  Interpolate values along the edges.
!
  do i = 2, m1 - 1

    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )

    call blend_101 ( r, x(1,1), x(m1,1), x(i,1) )

    call blend_101 ( r, x(1,m2), x(m1,m2), x(i,m2) )

  end do

  do j = 2, m2 - 1

    s = real (  j - 1, kind = 8 ) &
      / real ( m2 - 1, kind = 8 )

    call blend_101 ( s, x(1,1), x(1,m2), x(1,j) )

    call blend_101 ( s, x(m1,1), x(m1,m2), x(m1,j) )

  end do
!
!  Interpolate values in the interior.
!
  do i = 2, m1 - 1

    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )

    do j = 2, m2 - 1

      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )

      call blend_112 ( r, s, x(1,1), x(1,m2), x(m1,1), x(m1,m2), &
        x(i,1), x(i,m2), x(1,j), x(m1,j), x(i,j) )

    end do

  end do

  return
end
subroutine blend_ij_1d1 ( x, m1, m2 )

!*****************************************************************************80
!
!! BLEND_IJ_1D1 extends indexed scalar data along edges into a table.
!
!  Diagram:
!
!    ( X11,  X12,  X13,  X14,  X15,  X16,  X1M2  )
!    ( X21,  ...,  ...,  ...,  ...,  ...,  X2M2  )
!    ( X31,  ...,  ...,  ...,  ...,  ...,  X3M2  )
!    ( X41,  ...,  ...,  ...,  ...,  ...,  X4M2  )
!    ( XM11, XM12, XM13, XM14, XM15, XM16, XM1M2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(M1,M2).  
!    On input, data is contained in the "edge entries" X(1,J), X(I,1), 
!    X(M1,J) and X(I,M2), for I = 1 to M1, and J = 1 to M2.
!    On output, all entries in X have been assigned a value.
!
!    Input, integer ( kind = 4 ) M1, M2, the number of rows and columns in X.
!
  implicit none

  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m1,m2)
  real ( kind = 8 ) r
  real ( kind = 8 ) s
!
!  Interpolate values in the interior.
!
  do i = 2, m1 - 1

    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )

    do j = 2, m2 - 1

      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )

      call blend_112 ( r, s, x(1,1), x(1,m2), x(m1,1), x(m1,m2), &
        x(i,1), x(i,m2), x(1,j), x(m1,j), x(i,j) )

    end do

  end do

  return
end
subroutine blend_ij_w_1d1 ( x, r, s, m1, m2 )

!*****************************************************************************80
!
!! BLEND_IJ_W_1D1 extends weighted indexed scalar data along edges into a table.
!
!  Diagram:
!
!    Instead of assuming that the data in the table is equally spaced,
!    the arrays R and S are supplied, which should behave as
!    "coordinates" for the data.
!
!            S(1)  S(2)  S(3)  S(4)  S(5)  S(6)  S(M2)
!      
!    R(1)  ( X11,  X12,  X13,  X14,  X15,  X16,  X1M2  )
!    R(2)  ( X21,  ...,  ...,  ...,  ...,  ...,  X2M2  )
!    R(3)  ( X31,  ...,  ...,  ...,  ...,  ...,  X3M2  )
!    R(4)  ( X41,  ...,  ...,  ...,  ...,  ...,  X4M2  )
!    R(M1) ( XM11, XM12, XM13, XM14, XM15, XM16, XM1M2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(M1,M2).  
!    On input, data is contained in the "edge entries" X(1,J), X(I,1), 
!    X(M1,J) and X(I,M2), for I = 1 to M1, and J = 1 to M2.
!    On output, all entries in X have been assigned a value.
!
!    Input, real ( kind = 8 ) R(M1), S(M2), are "coordinates" for the rows and
!    columns of the array.  The values in R, and the values in S, should 
!    be strictly increasing or decreasing.
!
!    Input, integer ( kind = 4 ) M1, M2, the number of rows and columns in X.
!
  implicit none

  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m1,m2)
  real ( kind = 8 ) r(m1)
  real ( kind = 8 ) rr
  real ( kind = 8 ) s(m2)
  real ( kind = 8 ) ss
!
!  Interpolate values in the interior.
!
  do i = 2, m1 - 1

    rr = ( r(i) - r(1) ) / ( r(m1) - r(1) )

    do j = 2, m2 - 1

      ss = ( s(j) - s(1) ) / ( s(m2) - s(1) )

      call blend_112 ( rr, ss, x(1,1), x(1,m2), x(m1,1), x(m1,m2), &
        x(i,1), x(i,m2), x(1,j), x(m1,j), x(i,j) )

    end do

  end do

  return
end
subroutine blend_ijk_0d1 ( x, m1, m2, m3 )

!*****************************************************************************80
!
!! BLEND_IJK_0D1 extends indexed scalar corner data into a cubic table.
!
!  Diagram:
!
!    ( X111,   ...,  ...,  ...,  ...,  ...,  X1M21   )
!    ( ....,   ...,  ...,  ...,  ...,  ...,  ...     )
!    ( ....,   ...,  ...,  ...,  ...,  ...,  ...     )   First "layer"
!    ( ....,   ...,  ...,  ...,  ...,  ...,  ...     )
!    ( XM111,  ...,  ...,  ...,  ...,  ...,  XM1M21  )
!
!    ( ....,   ...,  ...,  ...,  ...,  ...,  ...     )
!    ( ....,   ...,  ...,  ...,  ...,  ...,  ...     )
!    ( ....,   ...,  ...,  ...,  ...,  ...,  ...     )   Middle "layers"
!    ( ....,   ...,  ...,  ...,  ...,  ...,  ...     )
!    ( ....,   ...,  ...,  ...,  ...,  ...,  ...     )
!
!    ( X11M3,  ...,  ...,  ...,  ...,  ...,  X1M2M3  )
!    ( ....,   ...,  ...,  ...,  ...,  ...,  ...     )
!    ( ....,   ...,  ...,  ...,  ...,  ...,  ...     )   Last "layer"
!    ( ....,   ...,  ...,  ...,  ...,  ...,  ...     )
!    ( XM11M3, ...,  ...,  ...,  ...,  ...,  XM1M2M3 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(M1,M2,M3).  
!    On input, X(1,1,1), X(1,M2,1), X(M1,1,1), X(M1,M2,1), X(1,1,M3),
!    X(1,M2,M3), X(M1,1,M3) and X(M1,M2,M3) contain scalar values 
!    which are to be interpolated throughout the table, using the table
!    indices I and J as independent variables. 
!    On output, all entries in X have been assigned a value.
!
!    Input, integer ( kind = 4 ) M1, M2, M3, the number of rows, columns, 
!    and layers in X.
!
  implicit none

  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x(m1,m2,m3)
!
!  Interpolate values along the "edges", that is, index triplets (i,j,k)
!  with exactly two of I, J, K an "extreme" value.
!
  do i = 2, m1 - 1

    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )

    call blend_101 ( r, x( 1, 1, 1), x(m1, 1, 1), x( i, 1, 1) )
    call blend_101 ( r, x( 1,m2, 1), x(m1,m2, 1), x( i,m2, 1) )
    call blend_101 ( r, x( 1, 1,m3), x(m1, 1,m3), x( i, 1,m3) )
    call blend_101 ( r, x( 1,m2,m3), x(m1,m2,m3), x( i,m2,m3) )

  end do

  do j = 2, m2 - 1

    s = real (  j - 1, kind = 8 ) &
      / real ( m2 - 1, kind = 8 )

    call blend_101 ( s, x( 1, 1, 1), x( 1,m2, 1), x( 1, j, 1) )
    call blend_101 ( s, x(m1, 1, 1), x(m1,m2, 1), x(m1, j, 1) )
    call blend_101 ( s, x( 1, 1,m3), x( 1,m2,m3), x( 1, j,m3) )
    call blend_101 ( s, x(m1, 1,m3), x(m1,m2,m3), x(m1, j,m3) )

  end do

  do k = 2, m3 - 1

    t = real (  k - 1, kind = 8 ) &
      / real ( m3 - 1, kind = 8 )

    call blend_101 ( t, x( 1, 1,1), x( 1, 1,m3), x( 1, 1,k) )
    call blend_101 ( t, x(m1, 1,1), x(m1, 1,m3), x(m1, 1,k) )
    call blend_101 ( t, x( 1,m2,1), x( 1,m2,m3), x( 1,m2,k) )
    call blend_101 ( t, x(m1,m2,1), x(m1,m2,m3), x(m1,m2,k) )
  end do
!
!  Interpolate values along the "faces", that is, index triplets (i,j,k)
!  with exactly one of I, J, K is an "extreme" value.
!
  do j = 2, m2 - 1

    s = real (  j - 1, kind = 8 ) & 
      / real ( m2 - 1, kind = 8 )

    do k = 2, m3 - 1

      t = real (  k - 1, kind = 8 ) &
        / real ( m3 - 1, kind = 8 )

      call blend_112 ( s, t, x(1,1,1), x(1,1,m3), x(1,m2,1), x(1,m2,m3), &
        x(1,j,1), x(1,j,m3), x(1,1,k), x(1,m2,k), x(1,j,k) )

      call blend_112 ( s, t, x(m1,1,1), x(m1,1,m3), x(m1,m2,1), x(m1,m2,m3), &
        x(m1,j,1), x(m1,j,m3), x(m1,1,k), x(m1,m2,k), x(m1,j,k) )

    end do
  end do 

  do i = 2, m1 - 1

    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )

    do k = 2, m3 - 1

      t = real (  k - 1, kind = 8 ) &
        / real ( m3 - 1, kind = 8 )

      call blend_112 ( r, t, x(1,1,1), x(1,1,m3), x(m1,1,1), x(m1,1,m3), &
        x(i,1,1), x(i,1,m3), x(1,1,k), x(m1,1,k), x(i,1,k) )

      call blend_112 ( r, t, x(1,m2,1), x(1,m2,m3), x(m1,m2,1), x(m1,m2,m3), &
        x(i,m2,1), x(i,m2,m3), x(1,m2,k), x(m1,m2,k), x(i,m2,k) )

    end do
  end do 

  do i = 2, m1 - 1

    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )

    do j = 2, m2 - 1

      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )

      call blend_112 ( r, s, x(1,1,1), x(1,m2,1), x(m1,1,1), x(m1,m2,1), &
        x(i,1,1), x(i,m2,1), x(1,j,1), x(m1,j,1), x(i,j,1) )

      call blend_112 ( r, s, x(1,1,m3), x(1,m2,m3), x(m1,1,m3), x(m1,m2,m3), &
        x(i,1,m3), x(i,m2,m3), x(1,j,m3), x(m1,j,m3), x(i,j,m3) )

    end do
  end do
!
!  Interpolate values in the interior.
!
  do i = 2, m1 - 1

    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )

    do j = 2, m2 - 1

      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )

      do k = 2, m3 - 1

        t = real (  k - 1, kind = 8 ) &
          / real ( m3 - 1, kind = 8 )

        call blend_123 ( r, s, t, &
          x( 1,1,1), x( 1, 1,m3), x( 1,m2,1), x( 1,m2,m3), &
          x(m1,1,1), x(m1, 1,m3), x(m1,m2,1), x(m1,m2,m3), &
          x( i,1,1), x( i, 1,m3), x( i,m2,1), x( i,m2,m3), &
          x( 1,j,1), x( 1, j,m3), x(m1, j,1), x(m1, j,m3), &
          x( 1,1,k), x( 1,m2, k), x(m1, 1,k), x(m1,m2, k), &
          x( 1,j,k), x(m1, j, k), x( i, 1,k), x( i,m2, k), &
          x( i,j,1), x( i, j,m3), x( i, j,k) )

      end do

    end do

  end do

  return
end
subroutine blend_ijk_1d1 ( x, m1, m2, m3 )

!*****************************************************************************80
!
!! BLEND_IJK_1D1 extends indexed scalar edge data into a cubic table.
!
!  Diagram:
!
!    ( X111,   X121,   X131,   X141,   X151,   X1M21   )
!    ( X211,   ...,    ...,    ...,    ...,    X2M21   )
!    ( X311,   ...,    ...,    ...,    ...,    X3M21   )   Layer 1
!    ( X411,   ...,    ...,    ...,    ...,    X4M21   )
!    ( XM111,  XM121,  XM131,  XM141,  XM151,  XM1M21  )
!
!    ( X11K,   ...,    ...,    ...,    ...,    X1M2K   )
!    ( ....,   ...,    ...,    ...,    ...,    ...     )
!    ( ....,   ...,    ...,    ...,    ...,    ...     )   Layer K
!    ( ....,   ...,    ...,    ...,    ...,    ...     )   1 < K < M3
!    ( XM11K,  ...,    ...,    ...,    ...,    XM1M2K  )
!
!    ( X11M3,  X12M3,  X13M3,  X14M3,  X15M3,  X1M2M3  )
!    ( X21M3,  ...,    ...,    ...,    ...,    X2M2M3  )
!    ( X31M3,  ...,    ...,    ...,    ...,    X3M2M3  )   Layer M3
!    ( X41M3   ...,    ...,    ...,    ...,    X4M2M3  )
!    ( XM11M3, XM12M3, XM13M3, XM14M3, XM15M3, XM1M2M3 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(M1,M2,M3).  
!    On input, there is already scalar data in the entries X(I,J,K) 
!    corresponding to "edges" of the table, that is, entries for which
!    at least two of the three indices I, J and K are equal to their
!    minimum or maximum possible values.
!    On output, all entries in X have been assigned a value, using the
!    table indices as independent variables.
!
!    Input, integer ( kind = 4 ) M1, M2, M3, the number of rows, columns, and 
!    layers in X.
!
  implicit none

  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x(m1,m2,m3)
!
!  Interpolate values along the "faces", that is, index triplets (i,j,k)
!  where exactly one of I, J, K is an "extreme" value.
!
  do j = 2, m2 - 1

    s = real (  j - 1, kind = 8 ) &
      / real ( m2 - 1, kind = 8 )

    do k = 2, m3 - 1

      t = real (  k - 1, kind = 8 ) &
        / real ( m3 - 1, kind = 8 )

      call blend_112 ( s, t, x(1,1,1), x(1,1,m3), x(1,m2,1), &
        x(1,m2,m3), x(1,j,1), x(1,j,m3), x(1,1,k), x(1,m2,k), x(1,j,k) )

      call blend_112 ( s, t, x(m1,1,1), x(m1,1,m3), x(m1,m2,1), &
        x(m1,m2,m3), x(m1,j,1), x(m1,j,m3), x(m1,1,k), x(m1,m2,k), x(m1,j,k) )

    end do
  end do 

  do i = 2, m1 - 1

    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )

    do k = 2, m3 - 1

      t = real (  k - 1, kind = 8 ) &
        / real ( m3 - 1, kind = 8 )

      call blend_112 ( r, t, x(1,1,1), x(1,1,m3), x(m1,1,1), x(m1,1,m3), &
        x(i,1,1), x(i,1,m3), x(1,1,k), x(m1,1,k), x(i,1,k) )

      call blend_112 ( r, t, x(1,m2,1), x(1,m2,m3), x(m1,m2,1), x(m1,m2,m3), &
        x(i,m2,1), x(i,m2,m3), x(1,m2,k), x(m1,m2,k), x(i,m2,k) )

    end do
  end do 

  do i = 2, m1 - 1

    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )

    do j = 2, m2 - 1

      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )

      call blend_112 ( r, s, x(1,1,1), x(1,m2,1), x(m1,1,1), x(m1,m2,1), &
        x(i,1,1), x(i,m2,1), x(1,j,1), x(m1,j,1), x(i,j,1) )

      call blend_112 ( r, s, x(1,1,m3), x(1,m2,m3), x(m1,1,m3), x(m1,m2,m3), &
        x(i,1,m3), x(i,m2,m3), x(1,j,m3), x(m1,j,m3), x(i,j,m3) )

    end do
  end do
!
!  Interpolate values in the interior.
!
  do i = 2, m1 - 1

    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )

    do j = 2, m2 - 1

      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )

      do k = 2, m3 - 1

        t = real (  k - 1, kind = 8 ) &
          / real ( m3 - 1, kind = 8 )

        call blend_123 ( r, s, t, &
          x( 1,1,1), x( 1, 1,m3), x( 1,m2,1), x( 1,m2,m3), &
          x(m1,1,1), x(m1, 1,m3), x(m1,m2,1), x(m1,m2,m3), &
          x( i,1,1), x( i, 1,m3), x( i,m2,1), x( i,m2,m3), &
          x( 1,j,1), x( 1, j,m3), x(m1, j,1), x(m1, j,m3), &
          x( 1,1,k), x( 1,m2, k), x(m1, 1,k), x(m1,m2, k), &
          x( 1,j,k), x(m1, j, k), x( i, 1,k), x( i,m2, k), &
          x( i,j,1), x( i, j,m3), x( i, j,k) )

      end do

    end do

  end do

  return
end
subroutine blend_ijk_2d1 ( x, m1, m2, m3 )

!*****************************************************************************80
!
!! BLEND_IJK_2D1 extends indexed scalar face data into a cubic table.
!
!  Diagram:
!
!    ( X111    X121    X131    X141    X151    X1M21   )
!    ( X211    X221    X231    X241    X251    X2M21   )
!    ( X311    X321    X331    X341    X351    X3M21   )   Layer 1
!    ( X411    X421    X431    X441    X451    X4M21   )
!    ( XM111   XM121   XM131   XM141   XM151   XM1M21  )
!
!    ( X11K    X12K    X13K    X14K    X15K    X1M2K   )
!    ( X21K    ...     ....    ....    ....    X2M2K   )
!    ( X31K    ...     ....    ....    ....    X3M2K   )   Layer K
!    ( X41K    ...     ....    ....    ....    X4M2K   )   1 < K < M3
!    ( XM11K   XM12K   XM13K   XM14K   XM15K   XM1M2K  )
!
!    ( X11M3   X12M3   X13M3   X14M3   X15M3   X1M2M3  )
!    ( X21M3   X22M3   X23M3   X24M3   X25M3   X2M2M3  )
!    ( X31M3   X32M3   X33M3   X34M3   X35M3   X3M2M3  )   Layer M3
!    ( X41M3   X42M3   X43M3   X44M3   X45M3   X4M2M3  )
!    ( XM11M3  XM12M3  XM13M3  XM14M3  XM15M3  XM1M2M3 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(M1,M2,M3).  
!    On input, there is already scalar data in the entries X(I,J,K) 
!    corresponding to "faces" of the table, that is, entries for which
!    at least one of the three indices I, J and K is equal to their
!    minimum or maximum possible values.
!    On output, all entries in X have been assigned a value, using the
!    table indices as independent variables.
!
!    Input, integer ( kind = 4 ) M1, M2, M3, the number of rows, columns, and 
!    layers in X.
!
  implicit none

  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x(m1,m2,m3)
!
!  Interpolate values in the interior.
!
  do i = 2, m1 - 1

    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )

    do j = 2, m2 - 1

      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )

      do k = 2, m3 - 1

        t = real (  k - 1, kind = 8 ) &
          / real ( m3 - 1, kind = 8 )

        call blend_123 ( r, s, t, &
          x( 1,1,1), x( 1, 1,m3), x( 1,m2,1), x( 1,m2,m3), &
          x(m1,1,1), x(m1, 1,m3), x(m1,m2,1), x(m1,m2,m3), &
          x( i,1,1), x( i, 1,m3), x( i,m2,1), x( i,m2,m3), &
          x( 1,j,1), x( 1, j,m3), x(m1, j,1), x(m1, j,m3), &
          x( 1,1,k), x( 1,m2, k), x(m1, 1,k), x(m1,m2, k), &
          x( 1,j,k), x(m1, j, k), x( i, 1,k), x( i,m2, k), &
          x( i,j,1), x( i, j,m3), x( i, j,k) )

      end do

    end do

  end do

  return
end
subroutine blend_r_0dn ( r, x, n, bound_r )

!*****************************************************************************80
!
!! BLEND_R_0DN extends vector data at endpoints into a line.
!
!  Diagram:
!
!    0-----r-----1
!
!  Discussion:
!
!    This is simply linear interpolation.  BLEND_R_0DN is provided
!    mainly as a "base routine" which can be compared to its 
!    generalizations, such as BLEND_RS_0DN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the (R) coordinate of the point to 
!    be evaluated.
!
!    Output, real ( kind = 8 ) X(N), the interpolated value at the point (R).
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector space.
!
!    External, BOUND_R, is a subroutine which is given (R) coordinates
!    and an component value I, and returns XI, the value of the I-th 
!    component of the N-vector at that point.  BOUND_R will only be 
!    called for "corners", that is, for values (R) where R is either 
!    0.0 or 1.0.  BOUND_R has the form:
!
!      subroutine bound_r ( r, i, xi )
!
  implicit none

  integer ( kind = 4 ) n

  external bound_r
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1

  do i = 1, n
!
!  Get the I-th coordinate component at the two corners.
!
    call bound_r ( 0.0D+00, i, x0 )
    call bound_r ( 1.0D+00, i, x1 )
!
!  Interpolate the I-th coordinate component of the interior point.
!
    call blend_101 ( r, x0, x1, x(i) )

  end do

  return
end
subroutine blend_rs_0dn ( r, s, x, n, bound_rs )

!*****************************************************************************80
!
!! BLEND_RS_0DN extends vector data at corners into a square.
!
!  Diagram:
!
!    01-----r1-----11
!     |      .      |
!     |      .      |
!    0s.....rs.....1s
!     |      .      |
!     |      .      |
!    00-----r0-----10
!
!  Discussion:
!
!    BLEND_RS_0DN should be equivalent to the use of a bilinear finite 
!    element method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the (R,S) coordinates of the point to be 
!    evaluated.
!
!    Output, real ( kind = 8 ) X(N), the interpolated value at the point (R,S).
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector space.
!
!    External, BOUND_RS, is a subroutine which is given (R,S) 
!    coordinates and an component value I, and returns XI, the value 
!    of the I-th component of the N-vector at that point.  BOUND_RS 
!    will only be called for "corners", that is, for values (R,S) where 
!    R and S are either 0.0 or 1.0.  BOUND_RS has the form:
!      subroutine bound_rs ( r, s, i, xi )
!
  implicit none

  integer ( kind = 4 ) n

  external bound_rs
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x00
  real ( kind = 8 ) x01
  real ( kind = 8 ) x10
  real ( kind = 8 ) x11
  real ( kind = 8 ) xr0
  real ( kind = 8 ) xr1
  real ( kind = 8 ) x0s
  real ( kind = 8 ) x1s

  do i = 1, n
!
!  Get the I-th coordinate component at the four corners.
!
    call bound_rs ( 0.0D+00, 0.0D+00, i, x00 )
    call bound_rs ( 0.0D+00, 1.0D+00, i, x01 )
    call bound_rs ( 1.0D+00, 0.0D+00, i, x10 )
    call bound_rs ( 1.0D+00, 1.0D+00, i, x11 )
!
!  Interpolate the I-th coordinate component at the sides.
!
    call blend_101 ( r, x00, x10, xr0 )
    call blend_101 ( r, x01, x11, xr1 )
    call blend_101 ( s, x00, x01, x0s )
    call blend_101 ( s, x10, x11, x1s )
!
!  Interpolate the I-th coordinate component of the interior point.
!
    call blend_112 ( r, s, x00, x01, x10, x11, xr0, xr1, x0s, x1s, x(i) )

  end do

  return
end
subroutine blend_rs_1dn ( r, s, x, n, bound_rs )

!*****************************************************************************80
!
!! BLEND_RS_1DN extends vector data along sides into a square.
!
!  Diagram:
!
!    01-----r1-----11
!     |      .      |
!     |      .      |
!    0s.....rs.....1s
!     |      .      |
!     |      .      |
!    00-----r0-----10
!
!  Discussion:
!
!    BLEND_RS_1DN is NOT equivalent to a bilinear finite element method,
!    since the data is sampled everywhere along the boundary lines,
!    rather than at a finite number of nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the (R,S) coordinates of the point to be 
!    evaluated.
!
!    Output, real ( kind = 8 ) X(N), the interpolated value at the point (R,S).
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector space.
!
!    External, BOUND_RS, is a subroutine which is given (R,S) 
!    coordinates and an component value I, and returns XI, the value 
!    of the I-th component of the N-vector at that point.  BOUND_RS 
!    will only be called for "sides", that is, for values (R,S) where 
!    at least one of R and S is either 0.0 or 1.0.  BOUND_RS has the 
!    form:
!      subroutine bound_rs ( r, s, i, xi )
!
  implicit none

  integer ( kind = 4 ) n

  external bound_rs
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x00
  real ( kind = 8 ) x01
  real ( kind = 8 ) x10
  real ( kind = 8 ) x11
  real ( kind = 8 ) xr0
  real ( kind = 8 ) xr1
  real ( kind = 8 ) x0s
  real ( kind = 8 ) x1s

  do i = 1, n
!
!  Get the I-th coordinate component at the four corners.
!
    call bound_rs ( 0.0D+00, 0.0D+00, i, x00 )
    call bound_rs ( 0.0D+00, 1.0D+00, i, x01 )
    call bound_rs ( 1.0D+00, 0.0D+00, i, x10 )
    call bound_rs ( 1.0D+00, 1.0D+00, i, x11 )
!
!  Get the I-th coordinate component at the sides.
!
    call bound_rs ( r, 0.0D+00, i, xr0 )
    call bound_rs ( r, 1.0D+00, i, xr1 )
    call bound_rs ( 0.0D+00, s, i, x0s )
    call bound_rs ( 1.0D+00, s, i, x1s )
!
!  Interpolate the I-th coordinate component of the interior point.
!
    call blend_112 ( r, s, x00, x01, x10, x11, xr0, xr1, x0s, x1s, x(i) )

  end do

  return
end
subroutine blend_rst_0dn ( r, s, t, x, n, bound_rst )

!*****************************************************************************80
!
!! BLEND_RST_0DN extends vector data at corners into a cube.
!
!  Diagram:
!
!    010-----r10-----110        011-----r11-----111
!      |       .       |          |       .       |  
!      |       .       |          |       .       |
!    0s0.....rs0.....1s0        0s1.....rs1.....1s1     S
!      |       .       |          |       .       |     |
!      |       .       |          |       .       |     |
!    000-----r00-----100        001-----r01-----101     +----R
!           BOTTOM                      TOP
!
!    011-----0s1-----001        111-----1s1-----101
!      |       .       |          |       .       |  
!      |       .       |          |       .       |
!    01t.....0st.....00t        11t.....1st.....10t          T
!      |       .       |          |       .       |          |
!      |       .       |          |       .       |          |
!    010-----0s0-----000        110-----1s0-----100     S----+
!           LEFT                       RIGHT
!
!    001-----r01-----101        011-----r11-----111
!      |       .       |          |       .       |  
!      |       .       |          |       .       |
!    00t.....r0t.....100        01t.....r1t.....11t     T
!      |       .       |          |       .       |     |
!      |       .       |          |       .       |     |
!    000-----r00-----100        010-----r10-----110     +----R
!           FRONT                       BACK
!
!  Discussion:
!
!    BLEND_RST_0DN is equivalent to a trilinear finite element method.
!    Data along the edges, faces, and interior of the cube is 
!    interpolated from the data at the corners.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, T, the (R,S,T) coordinates of the 
!    point to be evaluated.
!
!    Output, real ( kind = 8 ) X(N), the interpolated value at the 
!    point (R,S,T).
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector space.
!
!    External, BOUND_RST, is a subroutine which is given (R,S,T) 
!    coordinates and an component value I, and returns XI, the value 
!    of the I-th component of the N-vector at that point.  BOUND_RST 
!    will only be called for "corners", that is, for values (R,S,T) 
!    where R, S and T are either 0.0 or 1.0.  BOUND_RST has the form:
!      subroutine bound_rst ( r, s, t, i, xi )
!
  implicit none

  integer ( kind = 4 ) n

  external bound_rst
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x000
  real ( kind = 8 ) x001
  real ( kind = 8 ) x010
  real ( kind = 8 ) x011
  real ( kind = 8 ) x100
  real ( kind = 8 ) x101
  real ( kind = 8 ) x110
  real ( kind = 8 ) x111
  real ( kind = 8 ) xr00
  real ( kind = 8 ) xr01
  real ( kind = 8 ) xr10
  real ( kind = 8 ) xr11
  real ( kind = 8 ) x0s0
  real ( kind = 8 ) x0s1
  real ( kind = 8 ) x1s0
  real ( kind = 8 ) x1s1
  real ( kind = 8 ) x00t
  real ( kind = 8 ) x01t
  real ( kind = 8 ) x10t
  real ( kind = 8 ) x11t
  real ( kind = 8 ) x0st
  real ( kind = 8 ) x1st
  real ( kind = 8 ) xr0t
  real ( kind = 8 ) xr1t
  real ( kind = 8 ) xrs0
  real ( kind = 8 ) xrs1

  do i = 1, n
!
!  Get the I-th coordinate component at the corners.
!
    call bound_rst ( 0.0D+00, 0.0D+00, 0.0D+00, i, x000 )
    call bound_rst ( 0.0D+00, 0.0D+00, 1.0D+00, i, x001 )
    call bound_rst ( 0.0D+00, 1.0D+00, 0.0D+00, i, x010 )
    call bound_rst ( 0.0D+00, 1.0D+00, 1.0D+00, i, x011 )
    call bound_rst ( 1.0D+00, 0.0D+00, 0.0D+00, i, x100 )
    call bound_rst ( 1.0D+00, 0.0D+00, 1.0D+00, i, x101 )
    call bound_rst ( 1.0D+00, 1.0D+00, 0.0D+00, i, x110 )
    call bound_rst ( 1.0D+00, 1.0D+00, 1.0D+00, i, x111 )
!
!  Interpolate the I-th coordinate component at the edges.
!
    call blend_101 ( r, x000, x100, xr00 )
    call blend_101 ( r, x001, x101, xr01 )
    call blend_101 ( r, x010, x110, xr10 )
    call blend_101 ( r, x011, x111, xr11 )

    call blend_101 ( s, x000, x010, x0s0 )
    call blend_101 ( s, x001, x011, x0s1 )
    call blend_101 ( s, x100, x110, x1s0 )
    call blend_101 ( s, x101, x111, x1s1 )

    call blend_101 ( t, x000, x001, x00t )
    call blend_101 ( t, x010, x011, x01t )
    call blend_101 ( t, x100, x101, x10t )
    call blend_101 ( t, x110, x111, x11t )
!
!  Interpolate the I-th component on the faces.
!
    call blend_112 ( s, t, x000, x001, x010, x011, x0s0, x0s1, x00t, &
      x01t, x0st )

    call blend_112 ( s, t, x100, x101, x110, x111, x1s0, x1s1, x10t, &
      x11t, x1st )

    call blend_112 ( r, t, x000, x001, x100, x101, xr00, xr01, x00t, &
      x10t, xr0t )

    call blend_112 ( r, t, x010, x011, x110, x111, xr10, xr11, x01t, &
      x11t, xr1t )

    call blend_112 ( r, s, x000, x010, x100, x110, xr00, xr10, x0s0, &
      x1s0, xrs0 )

    call blend_112 ( r, s, x001, x011, x101, x111, xr01, xr11, x0s1, &
      x1s1, xrs1 )
!
!  Interpolate the I-th coordinate component of the interior point.
!
    call blend_123 ( r, s, t, x000, x001, x010, x011, x100, x101, x110, x111, &
      xr00, xr01, xr10, xr11, x0s0, x0s1, x1s0, x1s1, x00t, x01t, x10t, x11t, &
      x0st, x1st, xr0t, xr1t, xrs0, xrs1, x(i) )

  end do

  return
end
subroutine blend_rst_1dn ( r, s, t, x, n, bound_rst )

!*****************************************************************************80
!
!! BLEND_RST_1DN extends vector data on edges into a cube.
!
!  Diagram:
!
!    010-----r10-----110        011-----r11-----111
!      |       .       |          |       .       |  
!      |       .       |          |       .       |
!    0s0.....rs0.....1s0        0s1.....rs1.....1s1     S
!      |       .       |          |       .       |     |
!      |       .       |          |       .       |     |
!    000-----r00-----100        001-----r01-----101     +----R
!           BOTTOM                      TOP
!
!    011-----0s1-----001        111-----1s1-----101
!      |       .       |          |       .       |  
!      |       .       |          |       .       |
!    01t.....0st.....00t        11t.....1st.....10t          T
!      |       .       |          |       .       |          |
!      |       .       |          |       .       |          |
!    010-----0s0-----000        110-----1s0-----100     S----+
!           LEFT                       RIGHT
!
!    001-----r01-----101        011-----r11-----111
!      |       .       |          |       .       |  
!      |       .       |          |       .       |
!    00t.....r0t.....100        01t.....r1t.....11t     T
!      |       .       |          |       .       |     |
!      |       .       |          |       .       |     |
!    000-----r00-----100        010-----r10-----110     +----R
!           FRONT                       BACK
!
!  Discussion:
!
!    BLEND_RST_1D is NOT equivalent to a trilinear finite element method,
!    since the data is sampled everywhere along the corners and edges,
!    rather than at a finite number of nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, T, the (R,S,T) coordinates of the 
!    point to be evaluated.
!
!    Output, real ( kind = 8 ) X(N), the interpolated value at the 
!    point (R,S,T).
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector space.
!
!    External, BOUND_RST, is a subroutine which is given (R,S,T) 
!    coordinates and an component value I, and returns XI, the value 
!    of the I-th component of the N-vector at that point.  BOUND_RST 
!    will only be called for "edges", that is, for values (R,S,T) 
!    where at least two of R, S and T are either 0.0 or 1.0.  
!    BOUND_RST has the form:
!      subroutine bound_rst ( r, s, t, i, xi )
!
  implicit none

  integer ( kind = 4 ) n

  external bound_rst
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x000
  real ( kind = 8 ) x001
  real ( kind = 8 ) x010
  real ( kind = 8 ) x011
  real ( kind = 8 ) x100
  real ( kind = 8 ) x101
  real ( kind = 8 ) x110
  real ( kind = 8 ) x111
  real ( kind = 8 ) xr00
  real ( kind = 8 ) xr01
  real ( kind = 8 ) xr10
  real ( kind = 8 ) xr11
  real ( kind = 8 ) x0s0
  real ( kind = 8 ) x0s1
  real ( kind = 8 ) x1s0
  real ( kind = 8 ) x1s1
  real ( kind = 8 ) x00t
  real ( kind = 8 ) x01t
  real ( kind = 8 ) x10t
  real ( kind = 8 ) x11t
  real ( kind = 8 ) x0st
  real ( kind = 8 ) x1st
  real ( kind = 8 ) xr0t
  real ( kind = 8 ) xr1t
  real ( kind = 8 ) xrs0
  real ( kind = 8 ) xrs1

  do i = 1, n
!
!  Get the I-th coordinate component at the corners.
!
    call bound_rst ( 0.0D+00, 0.0D+00, 0.0D+00, i, x000 )
    call bound_rst ( 0.0D+00, 0.0D+00, 1.0D+00, i, x001 )
    call bound_rst ( 0.0D+00, 1.0D+00, 0.0D+00, i, x010 )
    call bound_rst ( 0.0D+00, 1.0D+00, 1.0D+00, i, x011 )
    call bound_rst ( 1.0D+00, 0.0D+00, 0.0D+00, i, x100 )
    call bound_rst ( 1.0D+00, 0.0D+00, 1.0D+00, i, x101 )
    call bound_rst ( 1.0D+00, 1.0D+00, 0.0D+00, i, x110 )
    call bound_rst ( 1.0D+00, 1.0D+00, 1.0D+00, i, x111 )
!
!  Get the I-th coordinate component at the edges.
!
    call bound_rst ( r, 0.0D+00, 0.0D+00, i, xr00 )
    call bound_rst ( r, 0.0D+00, 1.0D+00, i, xr01 )
    call bound_rst ( r, 1.0D+00, 0.0D+00, i, xr10 )
    call bound_rst ( r, 1.0D+00, 1.0D+00, i, xr11 )

    call bound_rst ( 0.0D+00, s, 0.0D+00, i, x0s0 )
    call bound_rst ( 0.0D+00, s, 1.0D+00, i, x0s1 )
    call bound_rst ( 1.0D+00, s, 0.0D+00, i, x1s0 )
    call bound_rst ( 1.0D+00, s, 1.0D+00, i, x1s1 )

    call bound_rst ( 0.0D+00, 0.0D+00, t, i, x00t )
    call bound_rst ( 0.0D+00, 1.0D+00, t, i, x01t )
    call bound_rst ( 1.0D+00, 0.0D+00, t, i, x10t )
    call bound_rst ( 1.0D+00, 1.0D+00, t, i, x11t )
!
!  Interpolate the I-th component on the faces.
!
    call blend_112 ( s, t, x000, x001, x010, x011, x0s0, x0s1, x00t, &
      x01t, x0st )

    call blend_112 ( s, t, x100, x101, x110, x111, x1s0, x1s1, x10t, &
      x11t, x1st )

    call blend_112 ( r, t, x000, x001, x100, x101, xr00, xr01, x00t, &
      x10t, xr0t )

    call blend_112 ( r, t, x010, x011, x110, x111, xr10, xr11, x01t, &
      x11t, xr1t )

    call blend_112 ( r, s, x000, x010, x100, x110, xr00, xr10, x0s0, &
      x1s0, xrs0 )

    call blend_112 ( r, s, x001, x011, x101, x111, xr01, xr11, x0s1, &
      x1s1, xrs1 )
!
!  Interpolate the I-th coordinate component of the interior point.
!
    call blend_123 ( r, s, t, x000, x001, x010, x011, x100, x101, x110, x111, &
      xr00, xr01, xr10, xr11, x0s0, x0s1, x1s0, x1s1, x00t, x01t, x10t, x11t, &
      x0st, x1st, xr0t, xr1t, xrs0, xrs1, x(i) )

  end do

  return
end
subroutine blend_rst_2dn ( r, s, t, x, n, bound_rst )

!*****************************************************************************80
!
!! BLEND_RST_2DN extends vector data on faces into a cube.
!
!  Diagram:
!
!    010-----r10-----110        011-----r11-----111
!      |       .       |          |       .       |  
!      |       .       |          |       .       |
!    0s0.....rs0.....1s0        0s1.....rs1.....1s1     S
!      |       .       |          |       .       |     |
!      |       .       |          |       .       |     |
!    000-----r00-----100        001-----r01-----101     +----R
!           BOTTOM                      TOP
!
!    011-----0s1-----001        111-----1s1-----101
!      |       .       |          |       .       |  
!      |       .       |          |       .       |
!    01t.....0st.....00t        11t.....1st.....10t          T
!      |       .       |          |       .       |          |
!      |       .       |          |       .       |          |
!    010-----0s0-----000        110-----1s0-----100     S----+
!           LEFT                       RIGHT
!
!    001-----r01-----101        011-----r11-----111
!      |       .       |          |       .       |  
!      |       .       |          |       .       |
!    00t.....r0t.....100        01t.....r1t.....11t     T
!      |       .       |          |       .       |     |
!      |       .       |          |       .       |     |
!    000-----r00-----100        010-----r10-----110     +----R
!           FRONT                       BACK
!
!  Discussion:
!
!    BLEND_RST_2DN is NOT equivalent to a trilinear finite element 
!    method, since the data is sampled everywhere along the corners, 
!    edges, and faces, rather than at a finite number of nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!    and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!    Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, T, the (R,S,T) coordinates of the point 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) X(N), the interpolated value at the 
!    point (R,S,T).
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector space.
!
!    External, BOUND_RST, is a subroutine which is given (R,S,T) 
!    coordinates and an component value I, and returns XI, the value 
!    of the I-th component of the N-vector at that point.  BOUND_RST 
!    will only be called for "faces", that is, for values (R,S,T) where 
!    at least one of R, S and T is either 0.0 or 1.0.  BOUND_RST has 
!    the form:
!      subroutine bound_rst ( r, s, t, i, xi )
!
  implicit none

  integer ( kind = 4 ) n

  external bound_rst
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x000
  real ( kind = 8 ) x001
  real ( kind = 8 ) x010
  real ( kind = 8 ) x011
  real ( kind = 8 ) x100
  real ( kind = 8 ) x101
  real ( kind = 8 ) x110
  real ( kind = 8 ) x111
  real ( kind = 8 ) xr00
  real ( kind = 8 ) xr01
  real ( kind = 8 ) xr10
  real ( kind = 8 ) xr11
  real ( kind = 8 ) x0s0
  real ( kind = 8 ) x0s1
  real ( kind = 8 ) x1s0
  real ( kind = 8 ) x1s1
  real ( kind = 8 ) x00t
  real ( kind = 8 ) x01t
  real ( kind = 8 ) x10t
  real ( kind = 8 ) x11t
  real ( kind = 8 ) x0st
  real ( kind = 8 ) x1st
  real ( kind = 8 ) xr0t
  real ( kind = 8 ) xr1t
  real ( kind = 8 ) xrs0
  real ( kind = 8 ) xrs1

  do i = 1, n
!
!  Get the I-th coordinate component at the corners.
!
    call bound_rst ( 0.0D+00, 0.0D+00, 0.0D+00, i, x000 )
    call bound_rst ( 0.0D+00, 0.0D+00, 1.0D+00, i, x001 )
    call bound_rst ( 0.0D+00, 1.0D+00, 0.0D+00, i, x010 )
    call bound_rst ( 0.0D+00, 1.0D+00, 1.0D+00, i, x011 )
    call bound_rst ( 1.0D+00, 0.0D+00, 0.0D+00, i, x100 )
    call bound_rst ( 1.0D+00, 0.0D+00, 1.0D+00, i, x101 )
    call bound_rst ( 1.0D+00, 1.0D+00, 0.0D+00, i, x110 )
    call bound_rst ( 1.0D+00, 1.0D+00, 1.0D+00, i, x111 )
!
!  Get the I-th coordinate component at the edges.
!
    call bound_rst ( r, 0.0D+00, 0.0D+00, i, xr00 )
    call bound_rst ( r, 0.0D+00, 1.0D+00, i, xr01 )
    call bound_rst ( r, 1.0D+00, 0.0D+00, i, xr10 )
    call bound_rst ( r, 1.0D+00, 1.0D+00, i, xr11 )

    call bound_rst ( 0.0D+00, s, 0.0D+00, i, x0s0 )
    call bound_rst ( 0.0D+00, s, 1.0D+00, i, x0s1 )
    call bound_rst ( 1.0D+00, s, 0.0D+00, i, x1s0 )
    call bound_rst ( 1.0D+00, s, 1.0D+00, i, x1s1 )

    call bound_rst ( 0.0D+00, 0.0D+00, t, i, x00t )
    call bound_rst ( 0.0D+00, 1.0D+00, t, i, x01t )
    call bound_rst ( 1.0D+00, 0.0D+00, t, i, x10t )
    call bound_rst ( 1.0D+00, 1.0D+00, t, i, x11t )
!
!  Get the I-th component on the faces.
!
    call bound_rst ( 0.0D+00, s, t, i, x0st )
    call bound_rst ( 1.0D+00, s, t, i, x1st )
    call bound_rst ( r, 0.0D+00, t, i, xr0t )
    call bound_rst ( r, 1.0D+00, t, i, xr1t )
    call bound_rst ( r, s, 0.0D+00, i, xrs0 )
    call bound_rst ( r, s, 1.0D+00, i, xrs1 )
!
!  Interpolate the I-th coordinate component of the interior point.
!
    call blend_123 ( r, s, t, x000, x001, x010, x011, x100, x101, x110, x111, &
      xr00, xr01, xr10, xr11, x0s0, x0s1, x1s0, x1s1, x00t, x01t, x10t, x11t, &
      x0st, x1st, xr0t, xr1t, xrs0, xrs1, x(i) )

  end do

  return
end
subroutine r8block_print ( l, m, n, a, title )

!*****************************************************************************80
!
!! R8BLOCK_PRINT prints a real block (a 3D matrix).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, N, the dimensions of the block.
!
!    Input, real ( kind = 8 ) A(L,M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(l,m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do k = 1, n

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  K = ', k
    write ( *, '(a)' ) ' '

    do jlo = 1, m, 5
      jhi = min ( jlo + 4, m )
      write ( *, '(a)' ) ' '
      write ( *, '(6x,5(i7,7x))' ) (j, j = jlo, jhi )
      write ( *, '(a)' ) ' '
      do i = 1, l
        write ( *, '(i6,5g14.6)' ) i, a(i,jlo:jhi,k)
      end do
    end do

  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
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
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)') j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(g14.6)' ) a(i,j)

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

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
