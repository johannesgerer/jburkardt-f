function arc_cosine ( c )

!*****************************************************************************80
!
!! ARC_COSINE computes the arc cosine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ACOS routine with an input argument that is
!    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant 
!    surprise (I did).
!
!    This routine simply truncates arguments outside the range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, the argument.
!
!    Output, real ( kind = 8 ) ARC_COSINE, an angle whose cosine is C.
!
  implicit none

  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) c
  real ( kind = 8 ) c2

  c2 = c
  c2 = max ( c2, -1.0D+00 )
  c2 = min ( c2, +1.0D+00 )

  arc_cosine = acos ( c2 )

  return
end
function arc_sine ( s )

!*****************************************************************************80
!
!! ARC_SINE computes the arc sine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ASIN routine with an input argument that is
!    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant 
!    surprise (I did).
!
!    This routine simply truncates arguments outside the range.
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
!    Input, real ( kind = 8 ) S, the argument.
!
!    Output, real ( kind = 8 ) ARC_SINE, an angle whose sine is S.
!
  implicit none

  real ( kind = 8 ) arc_sine
  real ( kind = 8 ) s
  real ( kind = 8 ) s2

  s2 = s
  s2 = max ( s2, - 1.0D+00 )
  s2 = min ( s2, + 1.0D+00 )

  arc_sine = asin ( s2 )

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
!    Output, real ( kind = 8 ) ATAN4, an angle between 0 and 2 * PI, 
!    whose tangent is (Y/X), and which lies in the appropriate quadrant so 
!    that the signs of its cosine and sine match those of X and Y.
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
      theta = PI
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
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
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

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
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
    seed = seed + i4_huge
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function r8vec_norm ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm

  r8vec_norm = sqrt ( sum ( a(1:n)**2 ) )

  return
end
subroutine r8vec_polarize ( n, a, p, a_normal, a_parallel )

!*****************************************************************************80
!
!! R8VEC_POLARIZE decomposes an R8VEC into normal and parallel components.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The (nonzero) vector P defines a direction.
!
!    The vector A can be written as the sum
!
!      A = A_normal + A_parallel
!
!    where A_parallel is a linear multiple of P, and A_normal
!    is perpendicular to P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the vector to be polarized.
!
!    Input, real ( kind = 8 ) P(N), the polarizing direction.
!
!    Output, real ( kind = 8 ) A_NORMAL(N), A_PARALLEL(N), the normal
!    and parallel components of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_dot_p
  real ( kind = 8 ) a_normal(n)
  real ( kind = 8 ) a_parallel(n)
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) p_norm

  p_norm = sqrt ( sum ( p(1:n)**2 ) )

  if ( p_norm == 0.0D+00 ) then
    a_normal(1:n) = a(1:n)
    a_parallel(1:n) = 0.0D+00
    return
  end if

  a_dot_p = dot_product ( a(1:n), p(1:n) ) / p_norm

  a_parallel(1:n) = a_dot_p * p(1:n) / p_norm

  a_normal(1:n) = a(1:n) - a_parallel(1:n)

  return
end
subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT concatenates two strings to make a third string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  if ( s1 == ' ' .and. s2 == ' ' ) then
    s3 = ' '
  else if ( s1 == ' ' ) then
    s3 = s2
  else if ( s2 == ' ' ) then
    s3 = s1
  else
    s3 = trim ( s1 ) // trim ( s2 )
  end if

  return
end
subroutine sphere01_distance_xyz ( xyz1, xyz2, dist )

!*****************************************************************************80
!
!! SPHERE01_DISTANCE_XYZ computes great circle distances on a unit sphere.
!
!  Discussion:
!
!    XYZ coordinates are used.
!
!    We assume the points XYZ1 and XYZ2 lie on the unit sphere.
!
!    This computation is a special form of the Vincenty formula.
!    It should be less sensitive to errors associated with very small 
!    or very large angular separations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    "Great-circle distance",
!    Wikipedia.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XYZ1(3), the coordinates of the first point.
!
!    Input, real ( kind = 8 ) XYZ2(3), the coordinates of the second point.
!
!    Output, real ( kind = 8 ) DIST, the great circle distance between
!    the points.
!
  implicit none

  real ( kind = 8 ) arc_sine
  real ( kind = 8 ) atan4
  real ( kind = 8 ) bot
  real ( kind = 8 ) dist
  real ( kind = 8 ) lat1
  real ( kind = 8 ) lat2
  real ( kind = 8 ) lon1
  real ( kind = 8 ) lon2
  real ( kind = 8 ) top
  real ( kind = 8 ) xyz1(3)
  real ( kind = 8 ) xyz2(3)

  lat1 = arc_sine ( xyz1(3) )
  lon1 = atan4 ( xyz1(2), xyz1(1) )

  lat2 = arc_sine ( xyz2(3) )
  lon2 = atan4 ( xyz2(2), xyz2(1) )

  top = ( cos ( lat2 ) * sin ( lon1 - lon2 ) )**2 &
      + ( cos ( lat1 ) * sin ( lat2 ) &
      -   sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) )**2

  top = sqrt ( top )

  bot = sin ( lat1 ) * sin ( lat2 ) &
      + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 )

  dist = atan2 ( top, bot )

  return
end
subroutine sphere01_sample ( n, seed, x )

!*****************************************************************************80
!
!! SPHERE01_SAMPLE picks random points on the unit sphere in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of samples.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(3,N), the sample points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arc_cosine
  integer ( kind = 4 ) j
  real ( kind = 8 ) phi
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta
  real ( kind = 8 ) vdot
  real ( kind = 8 ) x(3,n)

  do j = 1, n
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
    vdot = r8_uniform_01 ( seed )
    vdot = 2.0D+00 * vdot - 1.0D+00

    phi = arc_cosine ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
    theta = r8_uniform_01 ( seed )
    theta = 2.0D+00 * pi * theta

    x(1,j) = cos ( theta ) * sin ( phi )
    x(2,j) = sin ( theta ) * sin ( phi )
    x(3,j) = cos ( phi )

  end do

  return
end
subroutine sphere01_triangle_angles_to_area ( a, b, c, area )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_ANGLES_TO_AREA: area of a triangle on the unit sphere.
!
!  Discussion:
!
!    A unit sphere in 3D satisfies the equation:
!
!      X^2 + Y^2 + Z^2 = 1
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle on the unit sphere is:
!
!      AREA = A + B + C - PI
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the angles of the triangle.
!
!    Output, real ( kind = 8 ) AREA, the area of the sphere.
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
!
!  Apply Girard's formula.
!
  area = a + b + c - pi

  return
end
subroutine sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
  node_xyz )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_PROJECT projects from plane to spherical triangle.
!
!  Discussion:
!
!    We assume that points A, B and C lie on the unit sphere, and they
!    thus define a spherical triangle.
!
!    They also, of course, define a planar triangle.
!
!    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
!    planar triangle.
!
!    This function determines the coordinates of the point in the planar
!    triangle identified by the barycentric coordinates, and returns the
!    coordinates of the projection of that point onto the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the coordinates
!    of the points A, B, and C.
!
!    Input, integer ( kind = 4 ) F1, F2, F3, the barycentric coordinates
!    of a point in the triangle ABC.  Normally, these coordinates would
!    be real numbers, and would sum to 1.  For convenience, we allow these
!    to be integers which must be divided by F1+F2+F3.
!
!    Output, real ( kind = 8 ) NODE_XYZ(3), the coordinates of the 
!    point on the unit sphere which is the projection of the point on the plane
!    whose barycentric coordinates with respect to A, B, and C is
!    (F1,F2,F3)/(F1+F2+F3).
!
  implicit none

  real ( kind = 8 ) a_xyz(3)
  real ( kind = 8 ) b_xyz(3)
  real ( kind = 8 ) c_xyz(3)
  integer ( kind = 4 ) f1
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) f3
  real ( kind = 8 ) node_norm
  real ( kind = 8 ) node_xyz(3)
  real ( kind = 8 ) r8vec_norm

  node_xyz(1:3) = &
    ( real ( f1,           kind = 8 ) * a_xyz(1:3)   &
    + real (      f2,      kind = 8 ) * b_xyz(1:3)   &
    + real (           f3, kind = 8 ) * c_xyz(1:3) ) &
    / real ( f1 + f2 + f3, kind = 8 )

  node_norm = r8vec_norm ( 3, node_xyz(1:3) )

  node_xyz(1:3) = node_xyz(1:3) / node_norm

  return
end
subroutine sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
  node_xyz )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_PROJECT2 projects from plane to spherical triangle.
!
!  Discussion:
!
!    We assume that points A, B and C lie on the unit sphere, and they
!    thus define a spherical triangle.
!
!    They also, of course, define a planar triangle.
!
!    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
!    planar triangle.
!
!    This function determines the coordinates of the point in the planar
!    triangle identified by the barycentric coordinates, and returns the
!    coordinates of the projection of that point onto the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the coordinates
!    of the points A, B, and C.
!
!    Input, integer ( kind = 4 ) F1, F2, F3, the barycentric coordinates
!    of a point in the triangle ABC.  Normally, these coordinates would
!    be real numbers, and would sum to 1.  For convenience, we allow these
!    to be integers which must be divided by F1+F2+F3.
!
!    Output, real ( kind = 8 ) NODE_XYZ(3), the coordinates of the 
!    point on the unit sphere which is the projection of the point on the 
!    plane whose barycentric coordinates with respect to A, B, and C is
!    (F1,F2,F3)/(F1+F2+F3).
!
  implicit none

  real ( kind = 8 ) a_xyz(3)
  real ( kind = 8 ) ab(3)
  real ( kind = 8 ) ac(3)
  real ( kind = 8 ) acn(3)
  real ( kind = 8 ) acp(3)
  real ( kind = 8 ) angle
  real ( kind = 8 ) b_xyz(3)
  real ( kind = 8 ) bn(3)
  real ( kind = 8 ) bp(3)
  real ( kind = 8 ) c_xyz(3)
  real ( kind = 8 ) cn(3)
  real ( kind = 8 ) cp(3)
  integer ( kind = 4 ) f1
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) f3
  real ( kind = 8 ) node_xyz(3)
  real ( kind = 8 ) r8vec_norm
  real ( kind = 8 ) theta_ab
  real ( kind = 8 ) theta_ac
  real ( kind = 8 ) theta_bc
!
!  This check avoids 0/0 calculations later.
!
  if ( f2 == 0 .and. f3 == 0 ) then
    node_xyz(1:3) = a_xyz(1:3)
    return
  else if ( f1 == 0 .and. f3 == 0 ) then
    node_xyz(1:3) = b_xyz(1:3)
    return
  else if ( f1 == 0 .and. f2 == 0 ) then
    node_xyz(1:3) = c_xyz(1:3)
    return
  end if
!
!  Determine the angular distances (A,B) and (A,C).
!
  call sphere01_distance_xyz ( a_xyz, b_xyz, theta_ab )

  call sphere01_distance_xyz ( a_xyz, c_xyz, theta_ac )
!
!  Polarize B = BP + BN
!  Normalize BN, 
!  Same for C.
!
  call r8vec_polarize ( 3, b_xyz, a_xyz, bn, bp )
  bn(1:3) = bn(1:3) / r8vec_norm ( 3, bn )

  call r8vec_polarize ( 3, c_xyz, a_xyz, cn, cp )
  cn(1:3) = cn(1:3) / r8vec_norm ( 3, cn )
!
!  Determine AB and AC that use cos ( ( F2 + F3 ) / ( F1 + F2 + F3 ) ) of A
!  and cos ( F1 / ( F1 + F2 + F3 ) ) of B or C.
!
  angle = ( real ( f2 + f3, kind = 8 ) * theta_ab ) &
    / real ( f1 + f2 + f3, kind = 8 )
  ab(1:3) = cos ( angle ) * a_xyz(1:3) + sin ( angle ) * bn(1:3)

  angle = ( real ( f2 + f3, kind = 8 ) * theta_ac ) &
    / real ( f1 + f2 + f3, kind = 8 )
  ac(1:3) = cos ( angle ) * a_xyz(1:3) + sin ( angle ) * cn(1:3)
!
!  Determine the angular distance between AB and AC.
!
  call sphere01_distance_xyz ( ab(1:3), ac(1:3), theta_bc )
!
!  Polarize AC = ACP + ACN, normalize ACN.
!
  call r8vec_polarize ( 3, ac, ab, acn, acp )
  acn(1:3) = acn(1:3) / r8vec_norm ( 3, acn )
!
!  The interval between AB and AC is marked by F2+F3+1 vertices 0 through F2+F3.
!
  angle = ( real ( f3, kind = 8 ) * theta_bc ) / real ( f2 + f3, kind = 8 )

  node_xyz(1:3) = cos ( angle ) * ab(1:3) + sin ( angle ) * acn(1:3)

  return
end
subroutine sphere01_triangle_quad_00 ( n, v1, v2, v3, f, seed, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_00: quadrature over a triangle on the unit sphere.
!
!  Discussion:
!
!    This is a Monte Carlo approach.
!
!    The integral is approximated by averaging the values at N random points,
!    multiplied by the area.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of sample points.
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the XYZ coordinates of
!    the vertices of the triangle.
!
!    Input, real ( kind = 8 ), external :: F, evaluates the function
!    to be integrated, of the form:
!      function f ( v )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) v(3)
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) area
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) j
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  integer ( kind = 4 ) seed
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)
  real ( kind = 8 ) vc(3,n)

  call sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

  call sphere01_triangle_sample ( n, v1, v2, v3, seed, vc )

  quad = 0.0D+00
  do j = 1, n
    quad = quad + f ( vc(1:3,j) )
  end do

  result = quad * area / real ( n, kind = 8 )

  return
end
subroutine sphere01_triangle_quad_01 ( v1, v2, v3, f, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_01: quadrature over a triangle on the unit sphere.
!
!  Discussion:
!
!    The integral is approximated by the value at the centroid,
!    multiplied by the area.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the XYZ coordinates of
!    the vertices of the triangle.
!
!    Input, real ( kind = 8 ), external :: F, evaluates the function
!    to be integrated, of the form:
!      function f ( v )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) v(3)
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)
  real ( kind = 8 ) vc(3)

  call sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

  call sphere01_triangle_vertices_to_centroid ( v1, v2, v3, vc )

  quad = f ( vc )
  result = quad * area

  return
end
subroutine sphere01_triangle_quad_02 ( v1, v2, v3, f, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_02: quadrature over a triangle on the unit sphere.
!
!  Discussion:
!
!    The integral is approximated by the average of the vertex values,
!    multiplied by the area.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the XYZ coordinates of
!    the vertices of the triangle.
!
!    Input, real ( kind = 8 ), external :: F, evaluates the function
!    to be integrated, of the form:
!      function f ( v )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) v(3)
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)

  call sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

  quad = ( f ( v1 ) + f ( v2 ) + f ( v3 ) ) / 3.0D+00

  result = quad * area

  return
end
subroutine sphere01_triangle_quad_03 ( v1, v2, v3, f, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_03: quadrature over a triangle on the unit sphere.
!
!  Discussion:
!
!    The integral is approximated by the average of the midside values,
!    multiplied by the area.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the XYZ coordinates of
!    the vertices of the triangle.
!
!    Input, real ( kind = 8 ), external :: F, evaluates the function
!    to be integrated, of the form:
!      function f ( v )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) v(3)
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)
  real ( kind = 8 ) v4(3)
  real ( kind = 8 ) v5(3)
  real ( kind = 8 ) v6(3)

  call sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

  call sphere01_triangle_vertices_to_midpoints ( v1, v2, v3, v4, v5, v6 )

  quad = ( f ( v4 ) + f ( v5 ) + f ( v6 ) ) / 3.0D+00

  result = quad * area

  return
end
subroutine sphere01_triangle_quad_icos1c ( a_xyz, b_xyz, c_xyz, factor, &
  fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_ICOS1C: centroid rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over a spherical triangle on the
!    unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.  The centroids of these
!    triangles can be determined.  All of these calculations are done,
!    essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices and centroids projected to the sphere.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the vertices
!    of the spherical triangle.
!
!    Input, integer ( kind = 4 ) FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external :: FUN, evaluates the integrand, of the form:
!      subroutine fun ( n, x, v )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) v(n)
!      real ( kind = 8 ) x(3,n)
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of evaluation points.
!
!    Output, real ( kind = 8 ) RESULT, the estimated integral.
!
  implicit none

  real ( kind = 8 ) a_xyz(3)
  real ( kind = 8 ) a2_xyz(3)
  real ( kind = 8 ) area
  real ( kind = 8 ) area_total
  real ( kind = 8 ) b_xyz(3)
  real ( kind = 8 ) b2_xyz(3)
  real ( kind = 8 ) c_xyz(3)
  real ( kind = 8 ) c2_xyz(3)
  integer ( kind = 4 ) f1
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) f3
  integer ( kind = 4 ) factor
  real ( kind = 8 ), external :: fun
  integer ( kind = 4 ) node_num
  real ( kind = 8 ) node_xyz(3)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) result
  real ( kind = 8 ) v
!
!  Initialize the integral data.
!
  result = 0.0D+00
  area_total = 0.0D+00
  node_num = 0
!
!  Some subtriangles will have the same direction as the triangle.
!  Generate each in turn, by determining the barycentric coordinates
!  of the centroid (F1,F2,F3), from which we can also work out the barycentric
!  coordinates of the vertices of the subtriangle.
!
  do f3 = 1, 3 * factor - 2, 3
    do f2 = 1, 3 * factor - f3 - 1, 3

      f1 = 3 * factor - f3 - f2

      call sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
        node_xyz )

      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1 + 2, f2 - 1, f3 - 1, a2_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1 - 1, f2 + 2, f3 - 1, b2_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1 - 1, f2 - 1, f3 + 2, c2_xyz )

      call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

      v = fun ( node_xyz )    

      node_num = node_num + 1
      result = result + area * v
      area_total = area_total + area

    end do
  end do
!
!  The other subtriangles have the opposite direction from the triangle.
!  Generate each in turn, by determining the barycentric coordinates
!  of the centroid (F1,F2,F3), from which we can also work out the barycentric
!  coordinates of the vertices of the subtriangle.
!
  do f3 = 2, 3 * factor - 4, 3
    do f2 = 2, 3 * factor - f3 - 2, 3

      f1 = 3 * factor - f3 - f2

      call sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
        node_xyz )

      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1 - 2, f2 + 1, f3 + 1, a2_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1 + 1, f2 - 2, f3 + 1, b2_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1 + 1, f2 + 1, f3 - 2, c2_xyz )

      call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

      v = fun ( node_xyz )  

      node_num = node_num + 1  
      result = result + area * v
      area_total = area_total + area

    end do
  end do

  return
end
subroutine sphere01_triangle_quad_icos1m ( a_xyz, b_xyz, c_xyz, factor, &
  fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_ICOS1M: midpoint rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over a spherical triangle on the
!    unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.  The midpoints of the edges
!    of these triangles can be determined.  All of these calculations are done,
!    essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices and midpoints projected to the sphere.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the vertices
!    of the spherical triangle.
!
!    Input, integer ( kind = 4 ) FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external :: FUN, evaluates the integrand, of the form:
!      subroutine fun ( n, x, v )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) v(n)
!      real ( kind = 8 ) x(3,n)
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of evaluation points.
!
!    Output, real ( kind = 8 ) RESULT, the estimated integral.
!
  implicit none

  real ( kind = 8 ) a_xyz(3)
  real ( kind = 8 ) a2_xyz(3)
  real ( kind = 8 ) a3_xyz(3)
  real ( kind = 8 ) area
  real ( kind = 8 ) area_total
  real ( kind = 8 ) b_xyz(3)
  real ( kind = 8 ) b2_xyz(3)
  real ( kind = 8 ) b3_xyz(3)
  real ( kind = 8 ) c_xyz(3)
  real ( kind = 8 ) c2_xyz(3)
  real ( kind = 8 ) c3_xyz(3)
  integer ( kind = 4 ) f1
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) f3
  integer ( kind = 4 ) factor
  real ( kind = 8 ), external :: fun
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) result
  real ( kind = 8 ) va
  real ( kind = 8 ) vb
  real ( kind = 8 ) vc
!
!  Initialize the integral data.
!
  result = 0.0D+00
  area_total = 0.0D+00
  node_num = 0
!
!  Deal with subtriangles that have same orientation as face.
!
  do f1 = 0, factor - 1
    do f2 = 0, factor - f1 - 1
      f3 = factor - f1 - f2

      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1, a2_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1, b2_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

      call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 2 * f2 + 1, 2 * f3 - 2, a3_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, 2 * f1,     2 * f2 + 1, 2 * f3 - 1, b3_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 2 * f2,     2 * f3 - 1, c3_xyz )

      node_num = node_num + 3
      va = fun ( a3_xyz )
      vb = fun ( b3_xyz )   
      vc = fun ( c3_xyz )   
      result = result + area * ( va + vb + vc ) / 3.0D+00
      area_total = area_total + area

      end do
    end do
!
!  Deal with subtriangles that have opposite orientation as face.
!
  do f3 = 0, factor - 2
    do f2 = 1, factor - f3 - 1
      f1 = factor - f2 - f3

      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1, a2_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1, b2_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

      call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 2 * f2 - 1, 2 * f3 + 2, a3_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, 2 * f1,     2 * f2 - 1, 2 * f3 + 1, b3_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 2 * f2,     2 * f3 + 1, c3_xyz )

      node_num = node_num + 3
      va = fun ( a3_xyz )
      vb = fun ( b3_xyz )   
      vc = fun ( c3_xyz )  
      result = result + area * ( va + vb + vc ) / 3.0D+00
      area_total = area_total + area

    end do
  end do

  return
end
subroutine sphere01_triangle_quad_icos1v ( a_xyz, b_xyz, c_xyz, factor, &
  fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_ICOS1V: vertex rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over a spherical triangle on the
!    unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.   All of these calculations are 
!    done, essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices projected to the sphere.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the vertices
!    of the spherical triangle.
!
!    Input, integer ( kind = 4 ) FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external :: FUN, evaluates the integrand, of the form:
!      subroutine fun ( n, x, v )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) v(n)
!      real ( kind = 8 ) x(3,n)
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of evaluation points.
!
!    Output, real ( kind = 8 ) RESULT, the estimated integral.
!
  implicit none

  real ( kind = 8 ) a_xyz(3)
  real ( kind = 8 ) a2_xyz(3)
  real ( kind = 8 ) area
  real ( kind = 8 ) area_total
  real ( kind = 8 ) b_xyz(3)
  real ( kind = 8 ) b2_xyz(3)
  real ( kind = 8 ) c_xyz(3)
  real ( kind = 8 ) c2_xyz(3)
  integer ( kind = 4 ) f1
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) f3
  integer ( kind = 4 ) factor
  real ( kind = 8 ), external :: fun
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) result
  real ( kind = 8 ) va
  real ( kind = 8 ) vb
  real ( kind = 8 ) vc
!
!  Initialize the integral data.
!
  result = 0.0D+00
  area_total = 0.0D+00
  node_num = 0
!
!  Deal with subtriangles that have same orientation as face.
!
  do f1 = 0, factor - 1
    do f2 = 0, factor - f1 - 1
      f3 = factor - f1 - f2

      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1, a2_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1, b2_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

      call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

      node_num = node_num + 3
      va = fun ( a2_xyz )
      vb = fun ( b2_xyz )   
      vc = fun ( c2_xyz )   
      result = result + area * ( va + vb + vc ) / 3.0D+00
      area_total = area_total + area

      end do
    end do
!
!  Deal with subtriangles that have opposite orientation as face.
!
  do f3 = 0, factor - 2
    do f2 = 1, factor - f3 - 1
      f1 = factor - f2 - f3

      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1, a2_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1, b2_xyz )
      call sphere01_triangle_project ( &
        a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

      call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

      node_num = node_num + 3
      va = fun ( a2_xyz )
      vb = fun ( b2_xyz )   
      vc = fun ( c2_xyz )  
      result = result + area * ( va + vb + vc ) / 3.0D+00
      area_total = area_total + area

    end do
  end do

  return
end
subroutine sphere01_triangle_quad_icos2v ( a_xyz, b_xyz, c_xyz, factor, &
  fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_ICOS2V: vertex rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over a spherical triangle on the
!    unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.   All of these calculations are 
!    done, essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices projected to the sphere.  
!
!    This function uses a more sophisticated projection scheme than
!    SPHERE01_TRIANGLE_QUAD_ICOS1V, but this does not seem to improve
!    the results significantly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the vertices
!    of the spherical triangle.
!
!    Input, integer ( kind = 4 ) FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external :: FUN, evaluates the integrand, of the form:
!      subroutine fun ( n, x, v )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) v(n)
!      real ( kind = 8 ) x(3,n)
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of evaluation points.
!
!    Output, real ( kind = 8 ) RESULT, the estimated integral.
!
  implicit none

  real ( kind = 8 ) a_xyz(3)
  real ( kind = 8 ) a2_xyz(3)
  real ( kind = 8 ) area
  real ( kind = 8 ) area_total
  real ( kind = 8 ) b_xyz(3)
  real ( kind = 8 ) b2_xyz(3)
  real ( kind = 8 ) c_xyz(3)
  real ( kind = 8 ) c2_xyz(3)
  integer ( kind = 4 ) f1
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) f3
  integer ( kind = 4 ) factor
  real ( kind = 8 ), external :: fun
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) result
  real ( kind = 8 ) va
  real ( kind = 8 ) vb
  real ( kind = 8 ) vc
!
!  Initialize the integral data.
!
  result = 0.0D+00
  area_total = 0.0D+00
  node_num = 0
!
!  Deal with subtriangles that have same orientation as face.
!
  do f1 = 0, factor - 1
    do f2 = 0, factor - f1 - 1
      f3 = factor - f1 - f2

      call sphere01_triangle_project2 ( &
        a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1, a2_xyz )
      call sphere01_triangle_project2 ( &
        a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1, b2_xyz )
      call sphere01_triangle_project2 ( &
        a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

      call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

      node_num = node_num + 3
      va = fun ( a2_xyz )
      vb = fun ( b2_xyz )   
      vc = fun ( c2_xyz )   
      result = result + area * ( va + vb + vc ) / 3.0D+00
      area_total = area_total + area

      end do
    end do
!
!  Deal with subtriangles that have opposite orientation as face.
!
  do f3 = 0, factor - 2
    do f2 = 1, factor - f3 - 1
      f1 = factor - f2 - f3

      call sphere01_triangle_project2 ( &
        a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1, a2_xyz )
      call sphere01_triangle_project2 ( &
        a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1, b2_xyz )
      call sphere01_triangle_project2 ( &
        a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

      call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

      node_num = node_num + 3
      va = fun ( a2_xyz )
      vb = fun ( b2_xyz )   
      vc = fun ( c2_xyz )  
      result = result + area * ( va + vb + vc ) / 3.0D+00
      area_total = area_total + area

    end do
  end do

  return
end
subroutine sphere01_triangle_sample ( n, v1, v2, v3, seed, x )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_SAMPLE: sample spherical triangle on unit sphere.
!
!  Discussion:
!
!    The unit sphere has center 0 and radius 1.
!
!    A spherical triangle on the surface of the unit sphere contains those 
!    points with radius R = 1, bounded by the vertices V1, V2, V3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Arvo,
!    Stratified sampling of spherical triangles,
!    Computer Graphics Proceedings, Annual Conference Series, 
!    ACM SIGGRAPH '95, pages 437-438, 1995.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the XYZ coordinates of
!    the vertices of the spherical triangle.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(3,N), the XYZ coordinates of the 
!    sample points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) area
  real ( kind = 8 ) area_hat
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) j
  real ( kind = 8 ) q
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) r8vec_norm
  real ( kind = 8 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)
  real ( kind = 8 ) v31(3)
  real ( kind = 8 ) v4(3)
  real ( kind = 8 ) v42(3)
  real ( kind = 8 ) w
  real ( kind = 8 ) x(3,n)
  real ( kind = 8 ) xsi1
  real ( kind = 8 ) xsi2
  real ( kind = 8 ) z
!
!  Compute the sides, angles, and area of the spherical triangle;
!
  call sphere01_triangle_vertices_to_sides ( v1, v2, v3, a, b, c )

  call sphere01_triangle_sides_to_angles ( a, b, c, alpha, beta, gamma )

  call sphere01_triangle_angles_to_area ( alpha, beta, gamma, area )

  do j = 1, n
!
!  Select the new area.
!
    xsi1 = r8_uniform_01 ( seed )

    area_hat = xsi1 * area
!
!  Compute the sine and cosine of the angle phi.
!
    s = sin ( area_hat - alpha )
    t = cos ( area_hat - alpha )
!
!  Compute the pair that determines beta_hat.
!
    u = t - cos ( alpha )
    v = s + sin ( alpha ) * cos ( c )
!
!  Q is the cosine of the new edge length b_hat.
!
    q = ( ( v * t - u * s ) * cos ( alpha ) - v ) &
      / ( ( v * s + u * t ) * sin ( alpha ) )
!
!  We very occasionally get a Q value out of bounds.
!
    q = max ( q, - 1.0D+00 )
    q = min ( q, + 1.0D+00 )
!
!  V31 = normalized ( V3 - ( V3 dot V1 ) * V1 )
!
    w = dot_product ( v3, v1 )
    v31(1:3) = v3(1:3) - w * v1(1:3)
    t = r8vec_norm ( 3, v31(1:3) )
    if ( 0.0D+00 < t ) then
      v31(1:3) = v31(1:3) / t
    end if
!
!  V4 is the third vertex of the subtriangle V1, V2, V4.
!
    v4(1:3) = q * v1(1:3) + sqrt ( 1.0D+00 - q * q ) * v31(1:3)
!
!  Select cos theta, which will sample along the edge from V2 to V4.
!
    xsi2 = r8_uniform_01 ( seed )

    z = 1.0D+00 - xsi2 * ( 1.0D+00 - dot_product ( v4, v2 ) )
!
!  V42 = normalized ( V4 - ( V4 dot V2 ) * V2 )
!
    w = dot_product ( v4, v2 )

    v42(1:3) = v4(1:3) - w * v2(1:3)
    t = r8vec_norm ( 3, v42(1:3) )
    if ( 0.0D+00 < t ) then
      v42(1:3) = v42(1:3) / t
    end if
!
!  Construct the point.
!
    x(1:3,j) = z * v2(1:3) + sqrt ( 1.0D+00 - z * z ) * v42(1:3)

  end do

  return
end
subroutine sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_SIDES_TO_ANGLES: angles of triangle on unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AS, BS, CS, the (geodesic) length of the 
!    sides of the triangle.
!
!    Output, real ( kind = 8 ) A, B, C, the spherical angles of the triangle.
!    Angle A is opposite the side of length AS, and so on.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) as
  real ( kind = 8 ) asu
  real ( kind = 8 ) b
  real ( kind = 8 ) bs
  real ( kind = 8 ) bsu
  real ( kind = 8 ) c
  real ( kind = 8 ) cs
  real ( kind = 8 ) csu
  real ( kind = 8 ) ssu
  real ( kind = 8 ) tan_a2
  real ( kind = 8 ) tan_b2
  real ( kind = 8 ) tan_c2

  asu = as
  bsu = bs
  csu = cs
  ssu = ( asu + bsu + csu ) / 2.0D+00

  tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - asu )     ) )

  a = 2.0D+00 * atan ( tan_a2 )

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) )

  b = 2.0D+00 * atan ( tan_b2 )

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - csu )     ) )

  c = 2.0D+00 * atan ( tan_c2 )

  return
end
subroutine sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_AREA: area of triangle on unit sphere.
!
!  Discussion:
!
!    A unit sphere in 3D satisfies the equation:
!
!      X^2 + Y^2 + Z^2 = 1
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle on the unit sphere is:
!
!      AREA = A + B + C - PI
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real ( kind = 8 ) AREA, the area of the sphere.
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) a
  real ( kind = 8 ) as
  real ( kind = 8 ) b
  real ( kind = 8 ) bs
  real ( kind = 8 ) c
  real ( kind = 8 ) cs
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)
!
!  Compute the lengths of the sides of the spherical triangle.
!
  call sphere01_triangle_vertices_to_sides ( v1, v2, v3, as, bs, cs )
!
!  Get the spherical angles.
!
  call sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c )
!
!  Get the area.
!
  call sphere01_triangle_angles_to_area ( a, b, c, area )

  return
end
subroutine sphere01_triangle_vertices_to_centroid ( v1, v2, v3, vs )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_CENTROID: centroid of triangle on unit sphere.
!
!  Discussion:
!
!    A unit sphere in 3D satisfies the equation:
!
!      X^2 + Y^2 + Z^2 = 1
!
!    A spherical triangle is specified by three points on the sphere.
!
!    The (true) centroid of a spherical triangle is the point
!
!      VT = (XT,YT,ZT) = Integral ( X, Y, Z ) dArea / Integral 1 dArea
!
!    Note that the true centroid does NOT, in general, lie on the sphere.
!
!    The "flat" centroid VF is the centroid of the planar triangle defined by
!    the vertices of the spherical triangle.
!
!    The "spherical" centroid VS of a spherical triangle is computed by
!    the intersection of the geodesic bisectors of the triangle angles.
!    The spherical centroid lies on the sphere.
!
!    VF, VT and VS lie on a line through the center of the sphere.  We can
!    easily calculate VF by averaging the vertices, and from this determine
!    VS by normalizing.
!
!    (Of course, we still will not have actually computed VT, which lies
!    somewhere between VF and VS!)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real ( kind = 8 ) VS(3), the coordinates of the "spherical
!    centroid" of the spherical triangle.
!
  implicit none

  real ( kind = 8 ) norm
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)
  real ( kind = 8 ) vs(3)

  vs(1:3) = ( v1(1:3) + v2(1:3) + v3(1:3) ) / 3.0D+00

  norm = sqrt ( sum ( vs(1:3)**2 ) )

  vs(1:3) = vs(1:3) / norm

  return
end
subroutine sphere01_triangle_vertices_to_midpoints ( v1, v2, v3, m1, m2, m3 )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_MIDPOINTS: midsides of triangle on unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real ( kind = 8 ) M1(3), M2(3), M3(3), the coordinates of 
!    the midpoints of the sides of the spherical triangle.
!
  implicit none

  real ( kind = 8 ) m1(3)
  real ( kind = 8 ) m2(3)
  real ( kind = 8 ) m3(3)
  real ( kind = 8 ) norm
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)

  m1(1:3) = ( v1(1:3) + v2(1:3) ) / 2.0D+00
  norm = sqrt ( sum ( m1(1:3)**2 ) )
  m1(1:3) = m1(1:3) / norm

  m2(1:3) = ( v2(1:3) + v3(1:3) ) / 2.0D+00
  norm = sqrt ( sum ( m2(1:3)**2 ) )
  m2(1:3) = m2(1:3) / norm

  m3(1:3) = ( v3(1:3) + v1(1:3) ) / 2.0D+00
  norm = sqrt ( sum ( m3(1:3)**2 ) )
  m3(1:3) = m3(1:3) / norm

  return
end
subroutine sphere01_triangle_vertices_to_sides ( v1, v2, v3, as, bs, cs )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_SIDES: sides of triangle on unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the vertices of the spherical
!    triangle.
!
!    Output, real ( kind = 8 ) AS, BS, CS, the (geodesic) length of the 
!    sides of the triangle.
!
  implicit none

  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) as
  real ( kind = 8 ) bs
  real ( kind = 8 ) cs
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)

  as = arc_cosine ( dot_product ( v2(1:3), v3(1:3) ) )
  bs = arc_cosine ( dot_product ( v3(1:3), v1(1:3) ) )
  cs = arc_cosine ( dot_product ( v1(1:3), v2(1:3) ) )

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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  character ( len = 8 )  date
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  character ( len = 10 ) time
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y
  character ( len = 5 )  zone

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
