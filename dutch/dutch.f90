function angle_deg_2d ( p1, p2, p3 )

!*****************************************************************************80
!
!! ANGLE_DEG_2D returns the angle swept out between two rays in 2D.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_DEG_2D ( P1, P2, P3 ) + ANGLE_DEG_2D ( P3, P2, P1 ) = 360.0
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
!    Output, real ( kind = 8 ) ANGLE_DEG_2D, the angle swept out by the
!    rays, measured in degrees.  0 <= ANGLE_DEG_2D < 360.  If either ray
!    has zero length, then ANGLE_DEG_2D is set to 0.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) angle_deg_2d
  real ( kind = 8 ) angle_rad_2d
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( p(1) == 0.0D+00 .and. p(2) == 0.0D+00 ) then
    angle_deg_2d = 0.0D+00
    return
  end if

  angle_rad_2d = atan2 ( p(2), p(1) )

  if ( angle_rad_2d < 0.0D+00 ) then
    angle_rad_2d = angle_rad_2d + 2.0D+00 * pi
  end if

  angle_deg_2d = radians_to_degrees ( angle_rad_2d )

  return
end
function angle_rad_2d ( p1, p2, p3 )

!*****************************************************************************80
!
!! ANGLE_RAD_2D returns the angle in radians swept out between two rays in 2D.
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
!    15 January 2005
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

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( p(1) == 0.0D+00 .and. p(2) == 0.0D+00 ) then
    angle_rad_2d = 0.0D+00
    return
  end if

  angle_rad_2d = atan2 ( p(2), p(1) )

  if ( angle_rad_2d < 0.0D+00 ) then
    angle_rad_2d = angle_rad_2d + 2.0D+00 * pi
  end if

  return
end
subroutine circle_dia2imp_2d ( x1, y1, x2, y2, r, cx, cy )

!*****************************************************************************80
!
!! CIRCLE_DIA2IMP_2D converts a diameter to an implicit circle in 2D.
!
!  Discussion:
!
!    The diameter form of a circle is:
!
!      (X1,Y1) and (X2,Y2) are endpoints of a diameter.
!
!    The implicit form of a circle in 2D is:
!
!      (X-CX)**2 + (Y-CY)**2 = R**2
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
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, are the X and Y coordinates
!    of two points which form a diameter of the circle.
!
!    Output, real ( kind = 8 ) R, the computed radius of the circle.
!
!    Output, real ( kind = 8 ) CX, CY, the computed center of the circle.
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) cx
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) cy

  r = 0.5D+00 * sqrt ( ( x1 - x2 )**2 + ( y1 - y2 )**2 )

  cx = 0.5D+00 * ( x1 + x2 )
  cy = 0.5D+00 * ( y1 + y2 )

  return
end
subroutine circle_exp2imp_2d ( x1, y1, x2, y2, x3, y3, r, cx, cy )

!*****************************************************************************80
!
!! CIRCLE_EXP2IMP_2D converts a circle from explicit to implicit form in 2D.
!
!  Formula:
!
!    The explicit form of a circle in 2D is:
!
!      The circle passing through (X1,Y1), (X2,Y2), (X3,Y3).
!
!    The implicit form of a circle in 2D is:
!
!      (X-CX)**2 + (Y-CY)**2 = R**2
!
!  Discussion:
!
!    Any three points define a circle, as long as they don't lie on a straight
!    line.  (If the points do lie on a straight line, we could stretch the
!    definition of a circle to allow an infinite radius and a center at
!    some infinite point.)
!
!    Instead of the formulas used here, you can use the linear system
!    approach in the routine TRIANGLE_OUTCIRCLE_2D.
!
!    The diameter of the circle can be found by solving a 2 by 2 linear system.
!    This is because the vectors P2 - P1 and P3 - P1 are secants of the circle,
!    and each forms a right triangle with the diameter.  Hence, the dot product
!    of P2 - P1 with the diameter is equal to the square of the length
!    of P2 - P1, and similarly for P3 - P1.  These two equations determine the
!    diameter vector originating at P1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph O'Rourke,
!    Computational Geometry,
!    Cambridge University Press,
!    Second Edition, 1998, page 187.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, are the X and Y
!    coordinates of three points that lie on the circle.  These points should be
!    distinct, and not collinear.
!
!    Output, real ( kind = 8 ) R, the radius of the circle.  Normally, R will
!    be positive.  R is returned as -1 in the unlikely event that the points are
!    numerically collinear.
!
!    Output, real ( kind = 8 ) CX, CY, the center of the circle.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) r
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) cx
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) cy

  a = x2 - x1
  b = y2 - y1
  c = x3 - x1
  d = y3 - y1

  e = a * ( x1 + x2 ) + b * ( y1 + y2 )
  f = c * ( x1 + x3 ) + d * ( y1 + y3 )
!
!  Our formula is:
!
!    G = a * ( d - b ) - b * ( c - a )
!
!  but we get slightly better results using the original data.
!
  g = a * ( y3 - y2 ) - b * ( x3 - x2 )
!
!  We check for collinearity.  A more useful check would compare the
!  absolute value of G to a small quantity.
!
  if ( g == 0.0D+00 ) then
    cx = 0.0D+00
    cy = 0.0D+00
    r = -1.0D+00
    return
  end if
!
!  The center is halfway along the diameter vector from (X1,Y1).
!
  cx = 0.5D+00 * ( d * e - b * f ) / g
  cy = 0.5D+00 * ( a * f - c * e ) / g
!
!  Knowing the center, the radius is now easy to compute.
!
  r = sqrt ( ( x1 - cx )**2 + ( y1 - cy )**2 )

  return
end
subroutine circle_imp_contains_point_2d ( r, cx, cy, x, y, inside )

!*****************************************************************************80
!
!! CIRCLE_IMP_CONTAINS_POINT_2D: does an implicit circle contains a point in 2D?
!
!  Formula:
!
!    An implicit circle in 2D satisfies the equation:
!
!      ( X - CX )**2 + ( Y - CY )**2 = R**2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) CX, CY, the coordinates of the center
!    of the circle.
!
!    Input, real ( kind = 8 ) X, Y, the point to be checked.
!
!    Output, logical INSIDE, is TRUE if the point is inside or on the circle,
!    FALSE otherwise.
!
  implicit none

  logical inside
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) cx
  real ( kind = 8 ) y
  real ( kind = 8 ) cy

  if ( ( x - cx )**2 + ( y - cy )**2 <= r * r ) then
    inside = .true.
  else
    inside = .false.
  end if

  return
end
function cross0_2d ( p0, p1, p2 )

!*****************************************************************************80
!
!! CROSS0_2D finds the cross product of (P1-P0) and (P2-P0) in 2D.
!
!  Discussion:
!
!    Strictly speaking, the vectors lie in the (X,Y) plane, and
!    the cross product here is a vector in the Z direction.
!
!    The vectors are specified with respect to a basis point P0.
!    We are computing the normal to the triangle (P0,P1,P2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P0(2), P1(2), P2(2), the coordinates of
!    the three points.
!
!    Output, real ( kind = 8 ) CROSS0_2D, the Z component of the cross product
!    (P1-P0) x (P2-P0).
!
  implicit none

  real ( kind = 8 ) cross0_2d
  real ( kind = 8 ) p0(2)
  real ( kind = 8 ) p1(2)
  real ( kind = 8 ) p2(2)

  cross0_2d = ( p1(1) - p0(1) ) * ( p2(2) - p0(2) ) &
            - ( p1(2) - p0(2) ) * ( p2(1) - p0(1) )

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
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

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
!    Volume 8, Number 2, 1969, pages 136-143.
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
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
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
    seed = seed + i4_huge
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
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) wide

  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
    i4_wrap = ilo
  else
    i4_wrap = ilo + i4_modp ( ival-ilo, wide )
  end if

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = * ) title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 8 ) ctemp(incx)
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
      write ( ctemp(i2), '(i8)' ) i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a8)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine i4vec_frac ( n, a, k, iafrac )

!*****************************************************************************80
!
!! I4VEC_FRAC searches for the K-th smallest element in an N-vector.
!
!  Discussion:
!
!    Hoare's algorithm is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, integer ( kind = 4 ) A(N), array to search.  On output,
!    the elements of A have been somewhat rearranged.
!
!    Input, integer ( kind = 4 ) K, the fractile to be sought.  If K = 1, the
!    minimum entry is sought.  If K = N, the maximum is sought.
!    Other values of K search for the entry which is K-th in size.
!    K must be at least 1, and no greater than N.
!
!    Output, integer ( kind = 4 ) IAFRAC, the value of the K-th fractile of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) iafrac
  integer ( kind = 4 ) iryt
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) left

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal nonpositive value of N = ', n
    stop
  end if

  if ( k <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal nonpositive value of K = ', k
    stop
  end if

  if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal N < K, K = ', k
    stop
  end if

  left = 1
  iryt = n

  do

    if ( iryt <= left ) then
      iafrac = a(k)
      exit
    end if

    ix = a(k)
    i = left
    j = iryt

    do

      if ( j < i ) then

        if ( j < k ) then
          left = i
        end if

        if ( k < i ) then
          iryt = j
        end if

        exit

      end if
!
!  Find I so that IX <= A(I).
!
      do while ( a(i) < ix )
        i = i + 1
      end do
!
!  Find J so that A(J) <= IX.
!
      do while ( ix < a(j) )
        j = j - 1
      end do

      if ( i <= j ) then
        call i4_swap ( a(i), a(j) )
        i = i + 1
        j = j - 1
      end if

    end do

  end do

  return
end
subroutine i4vec_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_A reorders an array of integers into an ascending heap.
!
!  Discussion:
!
!    An ascending heap is an array A with the property that, for every index J,
!    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the smaller of the two values,
!  and update M if necessary.
!
        if ( a(m+1) < a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the small descendant is smaller than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( key <= a(m) ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an array of integers into an descending heap.
!
!  Discussion:
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_heap_d_extract ( n, a, val )

!*****************************************************************************80
!
!! I4VEC_HEAP_D_EXTRACT extracts the maximum value from a descending heap.
!
!  Discussion:
!
!    In other words, the routine finds the maximum value in the
!    heap, returns that value to the user, deletes that value from
!    the heap, and restores the heap to its proper form.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 150.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the heap.
!
!    Input/output, integer ( kind = 4 ) A(N), the heap.
!
!    Output, integer ( kind = 4 ) VAL, the item of maximum value, which has
!    been removed from the heap.
!
  implicit none

  integer ( kind = 4 ) a(*)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) val

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_HEAP_D_EXTRACT - Fatal error!'
    write ( *, '(a)' ) '  The heap is empty.'
    stop
  end if
!
!  Get the maximum value.
!
  val = a(1)

  if ( n == 1 ) then
    n = 0
    return
  end if
!
!  Shift the last value down.
!
  a(1) = a(n)
!
!  Restore the heap structure.
!
  n = n - 1
  call i4vec_sort_heap_d ( n, a )

  return
end
subroutine i4vec_heap_d_insert ( n, a, val )

!*****************************************************************************80
!
!! I4VEC_HEAP_D_INSERT inserts a new value into a descending heap.
!
!  Discussion:
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 150.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the heap.
!
!    Input/output, integer ( kind = 4 ) A(N), the heap.
!
!    Input, integer ( kind = 4 ) VAL, the value to be inserted.
!
  implicit none

  integer ( kind = 4 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) parent
  integer ( kind = 4 ) val

  n = n + 1
  i = n

  do while ( 1 < i )

    parent = i / 2

    if ( val <= a(parent) ) then
      exit
    end if

    a(i) = a(parent)
    i = parent

  end do

  a(i) = val

  return
end
subroutine i4vec_heap_d_max ( n, a, val_max )

!*****************************************************************************80
!
!! I4VEC_HEAP_D_MAX returns the maximum value in a descending heap of integers.
!
!  Discussion:
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 150.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the heap.
!
!    Input, integer ( kind = 4 ) A(N), the heap.
!
!    Output, integer ( kind = 4 ) VAL_MAX, the maximum value in the heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) val_max

  val_max = a(1)

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an integer vector to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_median ( n, a, median )

!*****************************************************************************80
!
!! I4VEC_MEDIAN returns the median of an unsorted integer vector.
!
!  Discussion:
!
!    Hoare's algorithm is used.  The values of the vector are
!    rearranged by this routine.
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
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to search.  On output,
!    the order of the elements of A has been somewhat changed.
!
!    Output, integer ( kind = 4 ) MEDIAN, the value of the median of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) median

  k = ( n + 1 ) / 2

  call i4vec_frac ( n, a, k, median )

  return
end
subroutine i4vec_pop ( n, x, stack1_max, stack1_num, stack1, stack2_max, &
  stack2_num, stack2 )

!*****************************************************************************80
!
!! I4VEC_POP pops an integer vector off of a stack.
!
!  Discussion:
!
!    If there are no more objects in the stack, N is returned as -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Output, integer ( kind = 4 ) X(*), the value of the vector.
!
!    Input, integer ( kind = 4 ) STACK1_MAX, the maximum size of STACK1.
!
!    Input/output, integer ( kind = 4 ) STACK1_NUM, the current size of STACK1.
!
!    Input/output, integer ( kind = 4 ) STACK1(STACK1_MAX), the vector
!    dimension stack.
!
!    Input, integer ( kind = 4 ) STACK2_MAX, the maximum size of STACK2.
!
!    Input/output, integer ( kind = 4 ) STACK2_NUM, the current size of STACK2.
!
!    Input/output, integer ( kind = 4 ) STACK2(STACK2_MAX), the vector value stack.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) stack1_max
  integer ( kind = 4 ) stack2_max

  integer ( kind = 4 ) stack1(stack1_max)
  integer ( kind = 4 ) stack1_num
  integer ( kind = 4 ) stack2(stack2_max)
  integer ( kind = 4 ) stack2_num
  integer ( kind = 4 ) x(*)

  if ( stack1_num < 1 ) then
    n = -1
    return
  end if

  n = stack1(stack1_num)
  stack1_num = stack1_num - 1

  stack2_num = stack2_num - n
  x(1:n) = stack2(stack2_num+1:stack2_num+n)

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an integer vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,i10)' ) i, a(i)
  end do

  return
end
subroutine i4vec_push ( n, x, stack1_max, stack1_num, stack1, stack2_max, &
  stack2_num, stack2 )

!*****************************************************************************80
!
!! I4VEC_PUSH pushes an integer vector onto a stack.
!
!  Discussion:
!
!    STACK1 contains a list of the dimensions of the objects stored.
!    Therefore, STACK1_MAX should be at least as big as the maximum number
!    of objects to be considered.
!
!    STACK2 contains the values of the objects.  Therefore, STACK2_MAX
!    should probably be as big as the maximum total length of the maximum
!    number of objects stored.
!
!    On first call, the user should have set STACK1_NUM and STACK2_NUM to zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.  N may be zero.
!
!    Input, integer ( kind = 4 ) X(N), the value of the vector.
!
!    Input, integer ( kind = 4 ) STACK1_MAX, the maximum size of STACK1.
!
!    Input/output, integer ( kind = 4 ) STACK1_NUM, the current size of STACK1.
!
!    Input/output, integer ( kind = 4 ) STACK1(STACK1_MAX), the vector
!    dimension stack.
!
!    Input, integer ( kind = 4 ) STACK2_MAX, the maximum size of STACK2.
!
!    Input/output, integer ( kind = 4 ) STACK2_NUM, the current size of STACK2.
!
!    Input/output, integer ( kind = 4 ) STACK2(STACK2_MAX), the vector value stack.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) stack1_max
  integer ( kind = 4 ) stack2_max

  integer ( kind = 4 ) stack1(stack1_max)
  integer ( kind = 4 ) stack1_num
  integer ( kind = 4 ) stack2(stack2_max)
  integer ( kind = 4 ) stack2_num
  integer ( kind = 4 ) x(n)

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_PUSH - Fatal error!'
    write ( *, '(a)' ) '  Input dimension N is negative.'
    stop
  end if

  if ( stack1_max < stack1_num + 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_PUSH - Fatal error!'
    write ( *, '(a)' ) '  Exceeding size of stack #1.'
    stop
  end if

  if ( stack2_max < stack2_num + n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_PUSH - Fatal error!'
    write ( *, '(a)' ) '  Exceeding size of stack #2.'
    stop
  end if

  stack1_num = stack1_num + 1
  stack1(stack1_num) = n

  stack2(stack2_num+1:stack2_num+n) = x(1:n)
  stack2_num = stack2_num + n

  return
end
subroutine i4vec_sort_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_D descending sorts an integer array using heap sort.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into ascending heap form.
!
  call i4vec_heap_a ( n, a )
!
!  2: Sort A.
!
!  The smallest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i4_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n-1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_a ( n1, a )
!
!  Take the smallest object from A(1) and move it to A(N1).
!
    call i4_swap ( a(1), a(n1) )

  end do

  return
end
subroutine i4vec_split_unsort ( n, a, split, isplit )

!*****************************************************************************80
!
!! I4VEC_SPLIT_UNSORT "splits" an unsorted I4VEC based on a splitting value.
!
!  Discussion:
!
!    If the vector is already sorted, it is simpler to do a binary search
!    on the data than to call this routine.
!
!    The vector is not assumed to be sorted before input, and is not
!    sorted during processing.  If sorting is not needed, then it is
!    more efficient to use this routine.
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
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to split.  On output,
!    all the entries of A that are less than or equal to SPLIT
!    are in A(1:ISPLIT).
!
!    Input, integer ( kind = 4 ) SPLIT, the value used to split the vector.
!    It is not necessary that any value of A actually equal SPLIT.
!
!    Output, integer ( kind = 4 ) ISPLIT, indicates the position of the last
!    entry of the split vector that is less than or equal to SPLIT.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) isplit
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) split
!
!  Partition the vector into A1, A2, A3, where
!    A1 = A(I1:J1) holds values <= SPLIT,
!    A2 = A(I2:J2) holds untested values,
!    A3 = A(I3:J3) holds values > SPLIT.
!
  i1 = 1
  j1 = 0

  i2 = 1
  j2 = n

  i3 = n+1
  j3 = n
!
!  Pick the next item from A2, and move it into A1 or A3.
!  Adjust indices appropriately.
!
  do i = 1, n

    if ( a(i2) <= split ) then
      i2 = i2 + 1
      j1 = j1 + 1
    else
      call i4_swap ( a(i2), a(i3-1) )
      i3 = i3 - 1
      j2 = j2 - 1
    end if

  end do

  isplit = j1

  return
end
subroutine ij_next ( i, j, n )

!*****************************************************************************80
!
!! IJ_NEXT returns the next matrix index.
!
!  Discussion:
!
!    For N = 3, the sequence of indices returned is:
!
!      (1,1), (1,2), (1,3), (2,1), (2,2), (2,3), (3,1), (3,2), (3,3), (0,0).
!
!    Note that once the value (N,N) is returned, the next value returned
!    will be (0,0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On input, the current pair of
!    indices.  On output, the next pair of indices.  If either index is illegal
!    on input, the output value of (I,J) will be (1,1).
!
!    Input, integer ( kind = 4 ) N, the maximum value for I and J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n

  if ( n < 1 ) then
    i = 0
    j = 0
    return
  end if

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then
    i = 1
    j = 1
    return
  end if

  if ( j < n ) then
    j = j + 1
  else if ( i < n ) then
    i = i + 1
    j = 1
  else
    i = 0
    j = 0
  end if

  return
end
subroutine ij_next_gt ( i, j, n )

!*****************************************************************************80
!
!! IJ_NEXT_GT returns the next matrix index, with the constraint that I < J.
!
!  Discussion:
!
!    For N = 5, the sequence of indices returned is:
!
!      (1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5), (4,5).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On input, the current pair of
!    indices.  On output, the next pair of indices.  If either index is illegal
!    on input, the output value of (I,J) will be (1,2).
!
!    Input, integer ( kind = 4 ) N, the maximum value for I and J.
!    A value of N less than 2 is nonsense.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n

  if ( n < 2 ) then
    i = 0
    j = 0
    return
  end if

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j .or. j <= i ) then
    i = 1
    j = 2
    return
  end if

  if ( j < n ) then
    j = j + 1
  else if ( i < n - 1 ) then
    i = i + 1
    j = i + 1
  else
    i = 0
    j = 0
  end if

  return
end
subroutine line_exp2imp_2d ( x1, y1, x2, y2, a, b, c )

!*****************************************************************************80
!
!! LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are
!    two points on the line. (X1,Y1) must be different
!    from (X2,Y2).
!
!    Output, real ( kind = 8 ) A, B, C, three coefficients which describe
!    the line that passes through (X1,Y1) and (X2,Y2).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
!
!  Take care of degenerate cases.
!
  if ( x1 == x2 .and. y1 == y2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP2IMP_2D - Fatal error!'
    write ( *, '(a)' ) '  (X1,Y1) = (X2,Y2)'
    write ( *, '(a,g14.6,2x,g14.6)' ) '  (X1,Y1) = ', x1, y1
    write ( *, '(a,g14.6,2x,g14.6)' ) '  (X2,Y2) = ', x2, y2
    stop
  end if

  a = y2 - y1
  b = x1 - x2
  c = x2 * y1 - x1 * y2

  return
end
subroutine line_exp_point_dist_2d ( x1, y1, x2, y2, x, y, dist )

!*****************************************************************************80
!
!! LINE_EXP_POINT_DIST_2D: distance ( explicit line, point ) in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
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
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are two
!    points on the line.
!
!    Input, real ( kind = 8 ) X, Y, the point whose distance from the line is
!    to be measured.
!
!    Output, real ( kind = 8 ) DIST, the distance from the point to the line.
!
  implicit none

  real ( kind = 8 ) bot
  real ( kind = 8 ) dist
  real ( kind = 8 ) dot
  real ( kind = 8 ) t
  real ( kind = 8 ) x
  real ( kind = 8 ) xn
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  real ( kind = 8 ) yn
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  bot = ( x1 - x2 )**2 + ( y1 - y2 )**2

  if ( bot == 0.0D+00 ) then

    xn = x1
    yn = y1
!
!  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
!  of the projection of (P-P1) onto (P2-P1).
!
  else

    dot = ( x - x1 ) * ( x2 - x1 ) + ( y - y1 ) * ( y2 - y1 )

    t = dot / bot

    xn = x1 + t * ( x2 - x1 )
    yn = y1 + t * ( y2 - y1 )

  end if

  dist = sqrt ( ( xn - x )**2 + ( yn - y )**2 )

  return
end
subroutine line_exp_point_dist_signed_2d ( x1, y1, x2, y2, x, y, dist_signed )

!*****************************************************************************80
!
!! LINE_EXP_POINT_DIST_SIGNED_2D: signed distance ( explicit line, point ).
!
!  Discussion:
!
!    The signed distance has two interesting properties:
!
!    *  The absolute value of the signed distance is the
!        usual (Euclidean) distance.
!
!    *  Points with signed distance 0 lie on the line,
!       points with a negative signed distance lie on one side
!         of the line,
!       points with a positive signed distance lie on the
!         other side of the line.
!
!    Assuming that C is nonnegative, then if a point is a positive
!    distance away from the line, it is on the same side of the
!    line as the point (0,0), and if it is a negative distance
!    from the line, it is on the opposite side from (0,0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, define the two points
!    (X1,Y1) and (X2,Y2) that determine the line.
!
!    Input, real ( kind = 8 ) X, Y, the point (X,Y) whose signed distance
!    is desired.
!
!    Output, real ( kind = 8 ) DIST_SIGNED, the signed distance from the
!    point to the line.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) dist_signed
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
!
!  Convert the line to implicit form.
!
  call line_exp2imp_2d ( x1, y1, x2, y2, a, b, c )
!
!  Compute the signed distance from the point to the line.
!
  dist_signed = ( a * x + b * y + c ) / sqrt ( a * a + b * b )

  return
end
subroutine line_seg_contains_point_2d ( x1, y1, x2, y2, x3, y3, u, v )

!*****************************************************************************80
!
!! LINE_SEG_CONTAINS_POINT_2D reports if a line segment contains a point in 2D.
!
!  Discussion:
!
!    In exact arithmetic, point P3 = (X3,Y3) is on the line segment between
!    P1=(X1,Y1) and P2=(X2,Y2) if and only if 0 <= U <= 1 and V = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the endpoints P1 and P2 of
!    a line segment.
!
!    Input, real ( kind = 8 ) X3, Y3, a point P3 to be tested.
!
!    Output, real ( kind = 8 ) U, the coordinate of (X3,Y3) along the axis from
!    with origin at P1 and unit at P2.
!
!    Output, real ( kind = 8 ) V, the magnitude of the off-axis portion of the
!    vector P3-P1, measured in units of (P2-P1).
!
  implicit none

  real ( kind = 8 ) u
  real ( kind = 8 ) unit
  real ( kind = 8 ) v
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  unit = sqrt ( ( x2 - x1 )**2 + ( y2 - y1 )**2 )

  if ( unit == 0.0D+00 ) then

    if ( x3 == x1 .and. y3 == y1 ) then
      u = 0.5D+00
      v = 0.0D+00
    else
      u = 0.5D+00
      v = huge ( 1.0D+00 )
    end if

  else

    u = ( ( x3 - x1 ) * ( x2 - x1 ) + ( y3 - y1 ) * ( y2 - y1 ) ) / unit**2

    v = sqrt ( ( ( u - 1.0D+00 ) * x1 - u * x2 + x3 )**2 &
             + ( ( u - 1.0D+00 ) * y1 - u * y2 + y3 )**2 ) / unit

  end if

  return
end
subroutine line_seg_vec_int_2d ( n, x1, y1, x2, y2, i, j, flag, xint, yint )

!*****************************************************************************80
!
!! LINE_SEG_VEC_INT_2D computes intersections of a set of line segments.
!
!  Discussion:
!
!    This is an implementation of the relatively inefficient algorithm
!    for computing all intersections of a set of line segments.
!
!    This is an "incremental" code, which returns as soon as it finds
!    a single intersection.  To find the next intersection, simply call
!    again.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of line segments.
!
!    Input, real ( kind = 8 ) X1(N), Y1(N), X2(N), Y2(N), the coordinates of the
!    endpoints of the line segments.
!
!    Input/output, integer ( kind = 4 ) I, J, used to keep track of the
!    computation.  On first call with a given set of line segments,
!    set I = J = 0.  On return with FLAG = 1, I and J record the indices of the
!    line segments whose intersection has just been found.  To find the
!    next intersection, simply call again, but do not alter I and J.
!
!    Output, integer ( kind = 4 ) FLAG:
!    0, no more intersections, the computation is done.
!    1, an intersection was detected between segments I and J.
!
!    Output, real ( kind = 8 ) XINT, YINT, the location of an intersection
!    of line segments I and J, if FLAG is 1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x1(n)
  real ( kind = 8 ) x2(n)
  real ( kind = 8 ) xint
  real ( kind = 8 ) y1(n)
  real ( kind = 8 ) y2(n)
  real ( kind = 8 ) yint

  do

    call ij_next_gt ( i, j, n )

    if ( i == 0 ) then
      flag = 0
      exit
    end if

    call lines_seg_int_2d ( x1(i), y1(i), x2(i), y2(i), x1(j), y1(j), &
      x2(j), y2(j), flag, xint, yint )

    if ( flag == 1 ) then
      return
    end if

  end do

  return
end
subroutine lines_exp_int_2d ( ival, x, y, x1, y1, x2, y2, x3, y3, x4, y4 )

!*****************************************************************************80
!
!! LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IVAL, reports on the intersection.
!
!     0, no intersection, the lines may be parallel or degenerate.
!     1, one intersection point, returned in X, Y.
!     2, infinitely many intersections, the lines are identical.
!
!    Output, real ( kind = 8 ) X, Y, if IVAl = 1, then X, Y contains
!    the intersection point.  Otherwise, X = 0, Y = 0.
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, define the first line.
!
!    Input, real ( kind = 8 ) X3, Y3, X4, Y4, define the second line.
!
  implicit none

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  integer ( kind = 4 ) ival
  logical point_1
  logical point_2
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4

  ival = 0
  x = 0.0D+00
  y = 0.0D+00
!
!  Check whether either line is a point.
!
  if ( x1 == x2 .and. y1 == y2 ) then
    point_1 = .true.
  else
    point_1 = .false.
  end if

  if ( x3 == x4 .and. y3 == y4 ) then
    point_2 = .true.
  else
    point_2 = .false.
  end if
!
!  Convert the lines to ABC format.
!
  if ( .not. point_1 ) then
    call line_exp2imp_2d ( x1, y1, x2, y2, a1, b1, c1 )
  end if

  if ( .not. point_2 ) then
    call line_exp2imp_2d ( x3, y3, x4, y4, a2, b2, c2 )
  end if
!
!  Search for intersection of the lines.
!
  if ( point_1 .and. point_2 ) then
    if ( x1 == x3 .and. y1 == y3 ) then
      ival = 1
      x = x1
      y = y1
    end if
  else if ( point_1 ) then
    if ( a2 * x1 + b2 * y1 == c2 ) then
      ival = 1
      x = x1
      y = y1
    end if
  else if ( point_2 ) then
    if ( a1 * x3 + b1 * y3 == c1 ) then
      ival = 1
      x = x3
      y = y3
    end if
  else
    call lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, x, y )
  end if

  return
end
subroutine lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, x, y )

!*****************************************************************************80
!
!! LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
!
!  Discussion:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1, B1, C1, define the first line.
!    At least one of A1 and B1 must be nonzero.
!
!    Input, real ( kind = 8 ) A2, B2, C2, define the second line.
!    At least one of A2 and B2 must be nonzero.
!
!    Output, integer ( kind = 4 ) IVAL, reports on the intersection.
!
!    -1, both A1 and B1 were zero.
!    -2, both A2 and B2 were zero.
!     0, no intersection, the lines are parallel.
!     1, one intersection point, returned in X, Y.
!     2, infinitely many intersections, the lines are identical.
!
!    Output, real ( kind = 8 ) X, Y, if IVAL = 1, then X, Y contains
!    the intersection point.  Otherwise, X = 0, Y = 0.
!
  implicit none

  real ( kind = 8 ) a(2,2)
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b(2,2)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) det
  integer ( kind = 4 ) ival
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  x = 0.0D+00
  y = 0.0D+00
!
!  Refuse to handle degenerate lines.
!
  if ( a1 == 0.0D+00 .and. b1 == 0.0D+00 ) then
    ival = -1
    return
  else if ( a2 == 0.0D+00 .and. b2 == 0.0D+00 ) then
    ival = -2
    return
  end if
!
!  Set up a linear system, and compute its inverse.
!
  a(1,1) = a1
  a(1,2) = b1
  a(2,1) = a2
  a(2,2) = b2

  call r8mat2_inverse ( a, b, det )
!
!  If the inverse exists, then the lines intersect.
!  Multiply the inverse times -C to get the intersection point.
!
  if ( det /= 0.0D+00 ) then

    ival = 1
    x = - b(1,1) * c1 - b(1,2) * c2
    y = - b(2,1) * c1 - b(2,2) * c2
!
!  If the inverse does not exist, then the lines are parallel
!  or coincident.  Check for parallelism by seeing if the
!  C entries are in the same ratio as the A or B entries.
!
  else

    ival = 0

    if ( a1 == 0.0D+00 ) then
      if ( b2 * c1 == c2 * b1 ) then
        ival = 2
      end if
    else
      if ( a2 * c1 == c2 * a1 ) then
        ival = 2
      end if
    end if

  end if

  return
end
subroutine lines_seg_dist_2d ( x1, y1, x2, y2, x3, y3, x4, y4, flag, x5, y5 )

!*****************************************************************************80
!
!! LINES_SEG_DIST_2D computes the distance of two line segments in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the endpoints of the first
!    segment.
!
!    Input, real ( kind = 8 ) X3, Y3, X4, Y4, the endpoints of the second
!    segment.
!
!    Output, integer ( kind = 4 ) FLAG, records the results.
!    0, the line segments do not intersect.
!    1, the line segments intersect.
!
!    Output, real ( kind = 8 ) X5, Y5.
!    If FLAG = 0, X5 = Y5 = 0.
!    If FLAG = 1, then (X5,Y5) is a point of intersection.
!
  implicit none

  integer ( kind = 4 ) flag
  integer ( kind = 4 ) ival
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) y5

  x5 = 0.0D+00
  y5 = 0.0D+00
!
!  Find the intersection of the two lines.
!
  call lines_exp_int_2d ( ival, x5, y5, x1, y1, x2, y2, x3, y3, x4, y4 )

  if ( ival == 0 ) then
    flag = 0
    return
  end if
!
!  Is the point on the first segment?
!
  call line_seg_contains_point_2d ( x1, y1, x2, y2, x5, y5, u, v )

  if ( u < 0 .or. 1.0D+00 < u .or. 0.001D+00 < v ) then
    flag = 0
    return
  end if
!
!  Is the point on the second segment?
!
  call line_seg_contains_point_2d ( x3, y3, x4, y4, x5, y5, u, v )

  if ( u < 0 .or. 1.0D+00 < u .or. 0.001D+00 < v ) then
    flag = 0
    return
  end if

  flag = 1

  return
end
subroutine lines_seg_int_1d ( x1, x2, x3, x4, flag, x5, x6 )

!*****************************************************************************80
!
!! LINES_SEG_INT_1D computes the intersection of two line segments in 1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, X2, the endpoints of the first segment.
!
!    Input, real ( kind = 8 ) X3, X4, the endpoints of the second segment.
!
!    Output, integer ( kind = 4 ) FLAG, records the results.
!    0, the line segments do not intersect.
!    1, the line segments intersect.
!
!    Output, real ( kind = 8 ) X5, X6, the endpoints of the intersection
!    segment.  If FLAG = 0, X5 = X6 = 0.
!
  implicit none

  integer ( kind = 4 ) flag
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) x6
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4

  y1 = min ( x1, x2 )
  y2 = max ( x1, x2 )
  y3 = min ( x3, x4 )
  y4 = max ( x3, x4 )

  flag = 0
  x5 = 0.0D+00
  x6 = 0.0D+00

  if ( y4 < y1 ) then
    return
  else if ( y2 < y3 ) then
    return
  end if

  flag = 1
  x5 = max ( y1, y3 )
  x6 = min ( y2, y4 )

  return
end
subroutine lines_seg_int_2d ( x1, y1, x2, y2, x3, y3, x4, y4, flag, x5, y5 )

!*****************************************************************************80
!
!! LINES_SEG_INT_2D computes the intersection of two line segments in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the endpoints of the first
!    segment.
!
!    Input, real ( kind = 8 ) X3, Y3, X4, Y4, the endpoints of the second
!    segment.
!
!    Output, integer ( kind = 4 ) FLAG, records the results.
!    0, the line segments do not intersect.
!    1, the line segments intersect.
!
!    Output, real ( kind = 8 ) X5, Y5.
!    If FLAG = 0, X5 = Y5 = 0.
!    If FLAG = 1, then (X5,Y5) is a point of intersection.
!
  implicit none

  integer ( kind = 4 ) flag
  integer ( kind = 4 ) ival
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) y5

  x5 = 0.0D+00
  y5 = 0.0D+00
!
!  Find the intersection of the two lines.
!
  call lines_exp_int_2d ( ival, x5, y5, x1, y1, x2, y2, x3, y3, x4, y4 )

  if ( ival == 0 ) then
    flag = 0
    return
  end if
!
!  Is the point on the first segment?
!
  call line_seg_contains_point_2d ( x1, y1, x2, y2, x5, y5, u, v )

  if ( u < 0 .or. 1.0D+00 < u .or. 0.001D+00 < v ) then
    flag = 0
    return
  end if
!
!  Is the point on the second segment?
!
  call line_seg_contains_point_2d ( x3, y3, x4, y4, x5, y5, u, v )

  if ( u < 0 .or. 1.0D+00 < u .or. 0.001D+00 < v ) then
    flag = 0
    return
  end if

  flag = 1

  return
end
subroutine perm_print ( n, p, title )

!*****************************************************************************80
!
!! PERM_PRINT prints a permutation.
!
!  Example:
!
!    Input:
!
!      P = 7 2 4 1 5 3 6
!
!    Printed output:
!
!      "This is the permutation:"
!
!      1 2 3 4 5 6 7
!      7 2 4 1 5 3 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation, in standard index form.
!
!    Input, character ( len = * ) TITLE, an optional title.
!    If no title is supplied, then only the permutation is printed.
!
  implicit none

  integer ( kind = 4 ), parameter :: inc = 20
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) p(n)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )

    do ilo = 1, n, inc
      ihi = min ( n, ilo + inc - 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(20i4)' ) ( i, i = ilo, ihi )
      write ( *, '(20i4)' ) p(ilo:ihi)
    end do

  else

    do ilo = 1, n, inc
      ihi = min ( n, ilo + inc - 1 )
      write ( *, '(20i4)' ) p(ilo:ihi)
    end do

  end if

  return
end
subroutine perm_random ( n, seed, p )

!*****************************************************************************80
!
!! PERM_RANDOM returns a random permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, page 78.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to permute.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) P(N), a permutation of the numbers
!    from 1 to N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  call i4vec_indicator ( n, p )

  do i = n, 2, -1
    j = i4_uniform ( 1, i, seed )
    call i4_swap ( p(i), p(j) )
  end do

  return
end
subroutine points_convex_hull_cubic_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_CONVEX_HULL_CUBIC_2D computes the convex hull of 2D points.
!
!  Discussion:
!
!    The algorithm used requires time that is cubic in the number of input
!    points.  Algorithms are available which are much more efficient!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 2-8.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of
!    the points.
!
!    Output, integer ( kind = 4 ) HULL_NUM, the number of vertices in the
!    convex hull.
!
!    Output, integer ( kind = 4 ) HULL(NODE_NUM), the HULL_NUM vertices that
!    form the convex hull, in counter clockwise order.
!
  implicit none

  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) b(node_num)
  real ( kind = 8 ) cross
  real ( kind = 8 ) cross0_2d
  integer ( kind = 4 ) e(node_num)
  integer ( kind = 4 ) hull(node_num)
  integer ( kind = 4 ) hull_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) match
  real ( kind = 8 ) node_xy(2,node_num)
  logical valid

  hull_num = 0

  do i = 1, node_num

    do j = 1,node_num
      if ( i /= j ) then
        valid = .true.
        do k = 1, node_num
          if ( k /= i .and. k /= j ) then
!
!  Compute the cross product: P(J)-P(I) x P(K)-P(I)
!  to determine if P(K) is strictly to the left of the line from P(I) to P(J).
!
!  I'll have to think if CROSS should be positive or negative...
!
            cross = cross0_2d ( node_xy(1:2,i), node_xy(1:2,j), node_xy(1:2,k) )

            if ( 0.0D+00 < cross ) then
              valid = .false.
              exit
            end if

          end if

        end do
!
!  Add the line from P(I) to P(J) to the list.
!
        if ( valid ) then
          hull_num = hull_num + 1
          b(hull_num) = i
          e(hull_num) = j
        end if

      end if

    end do

  end do
!
!  From the unordered edge list of vertex pairs,
!  construct an ordered list of vertices.
!
  hull(1) = b(1)
  match = e(1)

  do k = 2, hull_num

    hull(k) = 0

    do j = k, hull_num
      if ( b(j) == match ) then
        hull(k) = b(j)
        match = e(j)
        call i4_swap ( b(j), b(k) )
        call i4_swap ( e(j), e(k) )
        exit
      end if
    end do

    if ( hull(k) == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POINTS_2D_CONVEX_HULL_CUBIC - Fatal error!'
      write ( *, '(a,i8)' ) '  Algorithm failed for K = ', k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  B, E arrays:'
      write ( *, '(a)' ) ' '
      do i = 1, hull_num
        write ( *, '(2x,i8,2x,i8)' ) b(i), e(i)
      end do
      stop
    end if

  end do

  if ( match /= hull(1) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINTS_2D_CONVEX_HULL_CUBIC - Fatal error!'
    write ( *, '(a)' ) '  Algorithm failed to link last node to first.'
    stop
  end if
!
!  Reverse the order.
!
  hull(1:hull_num) = hull(hull_num:1:-1)

  return
end
subroutine points_convex_hull_nlogh_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_CONVEX_HULL_NLOGH_2D computes the convex hull of 2D points.
!
!  Discussion:
!
!    The work involved is N*log(H), where N is the number of points, and H is
!    the number of points that are on the hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) HULL_NUM, the number of nodes that lie
!    on the convex hull.
!
!    Output, integer ( kind = 4 ) HULL(NODE_NUM).  Entries 1 through HULL_NUM
!    contain the indices of the nodes that form the convex hull, in order.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle_max
  real ( kind = 8 ) angle_rad_2d
  real ( kind = 8 ) di
  real ( kind = 8 ) dr
  integer ( kind = 4 ) first
  integer ( kind = 4 ) hull(node_num)
  integer ( kind = 4 ) hull_num
  integer ( kind = 4 ) i
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) p_xy(2)
  integer ( kind = 4 ) q
  real ( kind = 8 ) q_xy(2)
  integer ( kind = 4 ) r
  real ( kind = 8 ) r_xy(2)

  if ( node_num < 1 ) then
    hull_num = 0
    return
  end if
!
!  If NODE_NUM = 1, the hull is the point.
!
  if ( node_num == 1 ) then
    hull_num = 1
    hull(1) = 1
    return
  end if
!
!  If NODE_NUM = 2, then the convex hull is either the two distinct points,
!  or possibly a single (repeated) point.
!
  if ( node_num == 2 ) then

    if ( node_xy(1,1) /= node_xy(1,2) .or. node_xy(2,1) /= node_xy(2,2) ) then
      hull_num = 2
      hull(1) = 1
      hull(2) = 2
    else
      hull_num = 1
      hull(1) = 1
    end if

    return

  end if
!
!  Find the leftmost point and call it "Q".
!  In case of ties, take the bottom-most.
!
  q = 1
  do i = 2, node_num
    if ( node_xy(1,i) < node_xy(1,q) .or. &
       ( node_xy(1,i) == node_xy(1,q) .and. &
         node_xy(2,i) < node_xy(2,q) ) ) then
      q = i
    end if
  end do

  q_xy(1:2) = node_xy(1:2,q)
!
!  Remember the starting point, so we know when to stop!
!
  first = q
  hull_num = 1
  hull(1) = q
!
!  For the first point, make a dummy previous point, 1 unit south,
!  and call it "P".
!
  p_xy(1) = q_xy(1)
  p_xy(2) = q_xy(2) - 1.0D+00
!
!  Now, having old point P, and current point Q, find the new point R
!  so the angle PQR is maximal.
!
!  Watch out for the possibility that the two nodes are identical.
!
  do

    r = 0
    angle_max = 0.0D+00

    do i = 1, node_num

      if ( i /= q .and. &
           ( node_xy(1,i) /= q_xy(1) .or. node_xy(2,i) /= q_xy(2) ) ) then

        angle = angle_rad_2d ( p_xy, q_xy, node_xy(1:2,i) )

        if ( r == 0 .or. angle_max < angle ) then

          r = i
          r_xy(1:2) = node_xy(1:2,r)
          angle_max = angle
!
!  In case of ties, choose the nearer point.
!
        else if ( r /= 0 .and. angle == angle_max ) then

          di = ( node_xy(1,i) - q_xy(1) )**2 + ( node_xy(2,i) - q_xy(2) )**2
          dr = ( r_xy(1)      - q_xy(1) )**2 + ( r_xy(2)      - q_xy(2) )**2

          if ( di < dr ) then
            r = i
            r_xy(1:2) = node_xy(1:2,r)
            angle_max = angle
          end if

        end if

      end if

    end do
!
!  We are done when we have returned to the first point on the convex hull.
!
    if ( r == first ) then
      exit
    end if

    hull_num = hull_num + 1

    if ( node_num < hull_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POINTS_CONVEX_HULL_NLOGH_2D - Fatal error!'
      write ( *, '(a)' ) '  The algorithm has failed.'
      stop
    end if
!
!  Add point R to convex hull.
!
    hull(hull_num) = r
!
!  Set P := Q, Q := R, and prepare to search for next point R.
!
    q = r

    p_xy(1:2) = q_xy(1:2)
    q_xy(1:2) = r_xy(1:2)

  end do
!
!  Reverse the order of the points.
!
  hull(1:hull_num) = hull(hull_num:1:-1)

  return
end
subroutine points_convex_hull_nlogn_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_CONVEX_HULL_NLOGN_2D computes the convex hull of 2D points.
!
!  Discussion:
!
!    The algorithm used requires time that is of order N * Log ( N ) in
!    N. the number of input points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 2-8.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of
!    the points.
!
!    Output, integer ( kind = 4 ) HULL_NUM, the number of vertices in the
!    convex hull.
!
!    Output, integer ( kind = 4 ) HULL(NODE_NUM), the HULL_NUM vertices that
!    form the convex hull, in counter clockwise order.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) angle_rad_2d
  integer ( kind = 4 ) hull(node_num)
  integer ( kind = 4 ) hull_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) l(node_num+100)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) p1(2)
  real ( kind = 8 ) p2(2)
  real ( kind = 8 ) p3(2)
  real ( kind = 8 ), parameter :: pi = 3.14159265D+00
!
!  Sort the (X,Y) points lexicographically.
!
  call r82vec_sort_quick_a ( node_num, node_xy )
!
!  Start with L1 = ( 1, 2 )
!
!  List L1 is stored in indices N1 through N2 of L.
!
  n1 = 1
  n2 = 2

  l(1) = 1
  p2(1:2) = node_xy(1:2,1)
  l(2) = 2
  p3(1:2) = node_xy(1:2,2)
!
!  Append I to L1.
!    Then, while L1 contains at least 3 points
!     if the last 3 points do not make a strict right turn, delete middle one.
!
  do i = 3, node_num

    n2 = n2 + 1
    l(n2) = i
    p1(1:2) = p2(1:2)
    p2(1:2) = p3(1:2)
    p3(1:2) = node_xy(1:2,i)

    do while ( 3 <= n2 + 1 - n1 )

      if ( angle_rad_2d ( p1, p2, p3 ) < pi ) then
        exit
      end if

      l(n2-1) = l(n2)
      n2 = n2 - 1
      p2(1:2) = p1(1:2)

      if ( 1 <= n2 - 2 ) then
        p1(1:2) = node_xy(1:2,l(n2-2))
      end if

    end do

  end do
!
!  Start with L2 = ( N, N-1 )
!
!  List L2 is stored in indices N2 through N3 of L.
!
  n3 = n2
  n3 = n3 + 1
  p2(1:2) = node_xy(1:2,node_num)
  l(n3) = node_num - 1
  p3(1:2) = node_xy(1:2,node_num-1)
!
!  Append I to L2.
!    Then, while L2 contains at least 3 points,
!    if the last 3 points do not make a strict right turn, delete middle one.
!
  do i = node_num-2, 1, -1

    n3 = n3 + 1
    l(n3) = i
    p1(1:2) = p2(1:2)
    p2(1:2) = p3(1:2)
    p3(1:2) = node_xy(1:2,i)

    do while ( 3 <= n3 + 1 - n2 )

      if ( angle_rad_2d ( p1, p2, p3 ) < pi ) then
        exit
      end if

      l(n3-1) = l(n3)
      n3 = n3 - 1
      p2(1:2) = p1(1:2)

      if ( n2 <= n3 - 2 ) then
        p1(1:2) = node_xy(1:2,l(n3-2))
      end if

    end do

  end do
!
!  Last entry of L is 1, so knock it off.
!
  hull_num = n3 - 1
!
!  Reverse the order.
!
  hull(1:hull_num) = l(hull_num:1:-1)

  return
end
subroutine points_minidisc1_2d ( n, px, py, qx, qy, r, cx, cy )

!*****************************************************************************80
!
!! POINTS_MINIDISC1_2D finds the smallest circle through Q containing points P.
!
!  Discussion:
!
!    It is assumed that there IS at least one disk which has the
!    points Q on its boundary, and which contains all the
!    points P.  For arbitrary points, this need not be the case.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 86-90.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points in the set P.
!
!    Input, real ( kind = 8 ) PX(N), PY(N), the X and Y coordinates of a set of
!    points in the plane.
!
!    Input, real ( kind = 8 ) QX, QY, a point in the plane.
!
!    Output, real ( kind = 8 ) R, CX, CY, the radius and center of the smallest
!    disk that encloses P and has Q on its boundary.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx
  real ( kind = 8 ) cy
  logical inside
  integer ( kind = 4 ) j
  real ( kind = 8 ) px(n)
  real ( kind = 8 ) py(n)
  real ( kind = 8 ) qx
  real ( kind = 8 ) qy
  real ( kind = 8 ) r
!
!  Determine the smallest disk with Q and P(1) on its boundary,
!  which is simply the circle whose diameter is the segment
!  from Q to P(1).
!
  call circle_dia2imp_2d ( qx, qy, px(1), py(1), r, cx, cy )
!
!  Now consider a point in the set P.  If it is not already
!  contained in the current circle, then expand the circle.
!
  do j = 2, n

    call circle_imp_contains_point_2d ( r, cx, cy, px(j), py(j), inside )

    if ( .not. inside ) then
      call points_minidisc2_2d ( j-1, px, py, px(j), py(j), qx, qy, r, cx, cy )
    end if

  end do

  return
end
subroutine points_minidisc2_2d ( n, px, py, q1x, q1y, q2x, q2y, r, cx, cy )

!*****************************************************************************80
!
!! POINTS_MINIDISC2_2D: smallest circle through Q1 and Q2 containing points P.
!
!  Discussion:
!
!    It is assumed that there IS at least one disk which has the
!    points Q1 and Q2 on its boundary, and which contains all the
!    points P.  For arbitrary points, this need not be the case.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 86-90.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points in the set P.
!
!    Input, real ( kind = 8 ) PX(N), PY(N), the X and Y coordinates of a set of
!    points in the plane.
!
!    Input, real ( kind = 8 ) Q1X, Q1Y, Q2X, Q2Y, two points in the plane.
!
!    Output, real ( kind = 8 ) R, CX, CY, the radius and center of the smallest
!    disk that encloses P and has Q1 and Q2 on its boundary.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx
  real ( kind = 8 ) cy
  logical inside
  integer ( kind = 4 ) k
  real ( kind = 8 ) px(n)
  real ( kind = 8 ) py(n)
  real ( kind = 8 ) q1x
  real ( kind = 8 ) q1y
  real ( kind = 8 ) q2x
  real ( kind = 8 ) q2y
  real ( kind = 8 ) r
!
!  Determine the smallest disk with Q1, Q2 on its boundary,
!  which is simply the circle whose diameter is the segment
!  from Q1 to Q2.
!
  call circle_dia2imp_2d ( q1x, q1y, q2x, q2y, r, cx, cy )
!
!  Now consider a point in the set P.  If it is not already
!  contained in the current circle, then expand the circle.
!
  do k = 1, n

    call circle_imp_contains_point_2d ( r, cx, cy, px(k), py(k), inside )

    if ( .not. inside ) then
      call circle_exp2imp_2d ( q1x, q1y, q2x, q2y, px(k), py(k), r, cx, cy )
    end if

  end do

  return
end
subroutine points_minidisc_2d ( n, px, py, r, cx, cy )

!*****************************************************************************80
!
!! POINTS_MINIDISC_2D finds the smallest circle containing points P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 86-90.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points in the set P.
!
!    Input, real ( kind = 8 ) PX(N), PY(N), the X and Y coordinates of a set of
!    points in the plane.
!
!    Output, real ( kind = 8 ) R, CX, CY, the radius and center of the smallest
!    disk that encloses P.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx
  real ( kind = 8 ) cy
  integer ( kind = 4 ) i
  logical inside
  real ( kind = 8 ) px(n)
  real ( kind = 8 ) py(n)
  real ( kind = 8 ) r
!
!  N = 1
!
  if ( n == 1 ) then
    r = 0.0D+00
    cx = px(1)
    cy = py(1)
    return
  end if
!
!  N = 2
!
  if ( n == 2 ) then
    r = 0.5D+00 * sqrt ( ( px(1) - px(2) )**2 + ( py(1) - py(2) )**2 )
    cx = 0.5D+00 * ( px(1) + px(2) )
    cy = 0.5D+00 * ( py(1) + py(2) )
    return
  end if
!
!  Determine the smallest disk with P(1), P(2) on its boundary,
!
  call circle_dia2imp_2d ( px(1), py(1), px(2), py(2), r, cx, cy )
!
!  Now consider a point in the set P.  If it is not already
!  contained in the current circle, then expand the circle.
!
  do i = 3, n

    call circle_imp_contains_point_2d ( r, cx, cy, px(i), py(i), inside )

    if ( .not. inside ) then
      call points_minidisc1_2d ( i-1, px, py, px(i), py(i), r, cx, cy )
    end if

  end do

  return
end
subroutine poly_triangulate_2d ( n, x, y, triang )

!*****************************************************************************80
!
!! POLY_TRIANGULATE_2D returns a triangulation of a polygon.
!
!  Discussion:
!
!    Given a polygon with N (distinct) nodes, a triangulation consists of
!    N-2 (distinct) triangles, with the property that each triangle is
!    constructed from vertices of the polygon, and that every point in the
!    polygon is contained in exactly one triangle, unless it lies on the
!    boundary of two or more triangles, and that no point exterior to the
!    polygon is contained in any triangle.  In other words, a triangulation
!    is a dissection of the area of the polygon into triangles.
!
!    This algorithm should be correct, but not efficient.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 46-49.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the polygon.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes,
!    listed in counter-clockwise order.
!
!    Output, integer ( kind = 4 ) TRIANG(3,N-2), the triangulation of
!    the polygon.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: stack2_max = 100

  integer ( kind = 4 ) best
  integer ( kind = 4 ) degree
  integer ( kind = 4 ) degree2
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_max
  logical inside
  integer ( kind = 4 ) j
  integer ( kind = 4 ) number
  integer ( kind = 4 ) poly(n)
  integer ( kind = 4 ) stack1(n-2)
  integer ( kind = 4 ) stack1_num
  integer ( kind = 4 ) stack2(stack2_max)
  integer ( kind = 4 ) stack2_num
  integer ( kind = 4 ) t
  integer ( kind = 4 ) triang(3,n-2)
  integer ( kind = 4 ) triang_num
  integer ( kind = 4 ) u
  integer ( kind = 4 ) v
  integer ( kind = 4 ) w
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4

  if ( n <= 2 ) then
    return
  end if

  degree = n
  call i4vec_indicator ( n, poly )

  call poly_reorder_nodes ( n, x, y, degree, poly )

  stack1_num = 0
  stack2_num = 0
  triang_num = 0

  do

    if ( degree == 3 ) then

      triang_num = triang_num + 1
      triang(1:3,triang_num) = poly(1:3)

      if ( stack1_num <= 0 ) then
        exit
      end if

      call i4vec_pop ( degree, poly, n-2, stack1_num, stack1, &
        stack2_max, stack2_num, stack2 )

      call poly_reorder_nodes ( n, x, y, degree, poly )

      cycle

    end if

    u = poly(degree)
    v = poly(1)
    w = poly(2)

    x1 = x(u)
    y1 = y(u)
    x2 = x(v)
    y2 = y(v)
    x3 = x(w)
    y3 = y(w)

    dist_max = 0.0D+00
    best = 0
    number = 0

    do j = 3, degree-1

      t = j
      x4 = x(poly(t))
      y4 = y(poly(t))

      call triangle_contains_point_2d ( x1, y1, x2, y2, x3, y3, x4, y4, inside )

      if ( inside ) then

        number = number + 1
        call line_exp_point_dist_2d ( x1, y1, x3, y3, x4, y4, dist )

        if ( dist_max <= dist ) then
          dist_max = dist
          best = t
        end if

      end if

    end do
!
!  If there were no polygonal nodes inside the triangle (U,V,W),
!  then triangle (U,V,W) can be used as part of the triangulation,
!  and the remaining region is simply the polygon (W,...,U).
!
    if ( number == 0 ) then

      triang_num = triang_num + 1
      triang(1,triang_num) = poly(1)
      triang(2,triang_num) = poly(2)
      triang(3,triang_num) = poly(degree)

      poly(1:n-1) = poly(2:n)
      degree = degree - 1

      call poly_reorder_nodes ( n, x, y, degree, poly )

    else
!
!  If there were polygonal nodes inside the triangle (U,V,W),
!  then POLY(T) is the node that was furthest from the line (U,W).
!
!  Split the polygon (V,W,part1,T,part2,U)
!  into (V,W,part1,T) and (V,T,part2,U)
!
      degree2 = t

      call i4vec_push ( degree2, poly(1:t), n-2, stack1_num, stack1, &
        stack2_max, stack2_num, stack2 )

      degree2 = degree + 2 - t
      poly(2:degree2) = poly(t:degree)
      degree = degree2

      call poly_reorder_nodes ( n, x, y, degree, poly )

    end if

  end do

  return
end
subroutine poly_reorder_nodes ( nxy, x, y, npoly, poly )

!*****************************************************************************80
!
!! POLY_REORDER_NODES reorders nodes of a polygon so node 1 is leftest lowest.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NXY, the number of nodes.
!
!    Input, real ( kind = 8 ) X(NXY), Y(NXY), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) NPOLY, the number of nodes of the polygon.
!
!    Input/output, POLY(NPOLY), the indices of the nodes.
!
  implicit none

  integer ( kind = 4 ) npoly
  integer ( kind = 4 ) nxy

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pmin
  integer ( kind = 4 ) poly(npoly)
  integer ( kind = 4 ) poly2(npoly)
  real ( kind = 8 ) x(nxy)
  real ( kind = 8 ) y(nxy)

  imin = 1
  pmin = poly(imin)

  do i = 2, npoly
    p = poly(i)
    if ( &
      ( x(p) < x(pmin) ) .or. &
      ( x(p) == x(pmin) .and. y(p) < y(pmin) ) ) then
      imin = i
      pmin = p
    end if
  end do

  if ( imin /= 1 ) then

    poly2(1:npoly+1-imin) = poly(imin:npoly)
    poly2(npoly+2-imin:npoly) = poly(1:imin-1)

    poly(1:npoly) = poly2(1:npoly)

  end if

  return
end
subroutine polycon_minkowski_sum_linear ( nu, ux, uy, nv, vx, vy, nw, wx, wy )

!*****************************************************************************80
!
!! POLYCON_MINKOWSKI_SUM_LINEAR: the Minkowski sum of two convex polygons.
!
!  Discussion:
!
!    For two geometric shapes U and V, the Minkowski sum W is the
!    set of all points
!      w = u + v
!    formed by adding points from the two shapes.
!
!    The Minkowski sum of two convex polygons is also a convex polygon.
!
!    The algorithm used here is only valid for convex polygons.
!
!    The vertices must be listed in counterclockwise order, with
!    the first vertex having the smallest Y coordinate.  If two
!    or more vertices have the same Y coordinate, then the one with
!    smallest X coordinate should be first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 275-281.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NU, the number of vertices of the first polygon.
!
!    Input, real ( kind = 8 ) UX(NU), UY(NU), the coordinates of the vertices
!    of the first polygon.
!
!    Input, integer ( kind = 4 ) NV, the number of vertices of the second polygon.
!
!    Input, real ( kind = 8 ) VX(NV), VY(NY), the coordinates of the vertices
!    of the second polygon.
!
!    Output, integer ( kind = 4 ) NW, the number of vertices of the sum polygon.
!    NW will be at most NV+NU.
!
!    Output, real ( kind = 8 ) WX(*), WY(*), the coordinates of the vertices
!    of the sum polygon.
!
  implicit none

  integer ( kind = 4 ) nu
  integer ( kind = 4 ) nv

  real ( kind = 8 ) angle_rad_2d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ip1w
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1w
  integer ( kind = 4 ) jw
  integer ( kind = 4 ) nw
  real ( kind = 8 ) p1(2)
  real ( kind = 8 ) p2(2)
  real ( kind = 8 ) p3(2)
  real ( kind = 8 ) q1(2)
  real ( kind = 8 ) q2(2)
  real ( kind = 8 ) q3(2)
  real ( kind = 8 ) ux(nu)
  real ( kind = 8 ) uy(nu)
  real ( kind = 8 ) u_angle
  logical u_done
  real ( kind = 8 ) vx(nv)
  real ( kind = 8 ) vy(nv)
  real ( kind = 8 ) v_angle
  logical v_done
  real ( kind = 8 ) wx(nu+nv+2)
  real ( kind = 8 ) wy(nu+nv+2)

  i = 1
  iw = 1
  u_done = .false.
  j = 1
  jw = 1
  v_done = .false.
  nw = 0

  do

    nw = nw + 1

    wx(nw) = ux(iw) + vx(jw)
    wy(nw) = uy(iw) + vy(jw)

    ip1w = i4_wrap ( i+1, 1, nu )

    p1(1) = ux(iw) + 1.0D+00
    p1(2) = uy(iw)
    p2(1) = ux(iw)
    p2(2) = uy(iw)
    p3(1) = ux(ip1w)
    p3(2) = uy(ip1w)

    u_angle = angle_rad_2d ( p1, p2, p3 )
    u_angle = u_angle + real ( i - nu, kind = 8 ) / real ( nu, kind = 8 ) &
      * 2.0D+00 * 3.141592653589793D+00

    jp1w = i4_wrap ( j+1, 1, nv )

    q1(1) = vx(jw) + 1.0D+00
    q1(2) = vy(jw)
    q2(1) = vx(jw)
    q2(2) = vy(jw)
    q3(1) = vx(jp1w)
    q3(2) = vy(jp1w)

    v_angle = angle_rad_2d ( q1, q2, q3 )
    v_angle = v_angle + real ( j - nv, kind = 8 ) / real ( nv, kind = 8 ) &
      * 2.0D+00 * 3.141592653589793D+00

    if ( u_angle < v_angle ) then
      i = i + 1
    else if ( v_angle < u_angle ) then
      j = j + 1
    else
      i = i + 1
      j = j + 1
    end if

    if ( i == nu + 1 ) then
      u_done = .true.
    end if

    if ( j == nv + 1 ) then
      v_done = .true.
    end if

    if ( u_done .and. v_done ) then
      exit
    end if

    iw = i4_wrap ( i, 1, nu )
    jw = i4_wrap ( j, 1, nv )

  end do

  return
end
subroutine polycon_minkowski_sum_n2logn2 ( nu, ux, uy, nv, vx, vy, nw, w_xy )

!*****************************************************************************80
!
!! POLYCON_MINKOWSKI_SUM_N2LOGN2 Minkowski sums two convex polygons.
!
!  Discussion:
!
!    For two geometric shapes U and V, the Minkowski sum W is the
!    set of all points
!      w = u + v
!    formed by adding points from the two shapes.
!
!    The Minkowski sum of two convex polygons is also a convex polygon.
!
!    The algorithm used here is only valid for convex polygons.
!
!    The vertices must be listed in counterclockwise order, with
!    the first vertex having the smallest Y coordinate.  If two
!    or more vertices have the same Y coordinate, then the one with
!    smallest X coordinate should be first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 275-281.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NU, the number of vertices of the first polygon.
!
!    Input, real ( kind = 8 ) UX(NU), UY(NU), the coordinates of the vertices
!    of the first polygon.
!
!    Input, integer ( kind = 4 ) NV, the number of vertices of the
!    second polygon.
!
!    Input, real ( kind = 8 ) VX(NV), VY(NY), the coordinates of the vertices
!    of the second polygon.
!
!    Output, integer ( kind = 4 ) NW, the number of vertices of the sum polygon.
!    NW will be at most NU+NV.
!
!    Output, real ( kind = 8 ) W_XY(2,*), the coordinates of the vertices
!    of the sum polygon.
!
  implicit none

  integer ( kind = 4 ) nu
  integer ( kind = 4 ) nv

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nuv
  integer ( kind = 4 ) nw
  real ( kind = 8 ) uv_xy(2,nu*nv)
  real ( kind = 8 ) ux(nu)
  real ( kind = 8 ) uy(nu)
  real ( kind = 8 ) vx(nv)
  real ( kind = 8 ) vy(nv)
  integer ( kind = 4 ) w(nu+nv)
  real ( kind = 8 ) w_xy(2,nu+nv)
!
!  Generate points from all pairs.
!
  nuv = 0
  do i = 1, nu
    do j = 1, nv
      nuv = nuv + 1
      uv_xy(1,nuv) = ux(i) + vx(j)
      uv_xy(2,nuv) = uy(i) + vy(j)
    end do
  end do
!
!  Compute the convex hull.
!  The output value of NW should be no more than NU+NV.
!
  call points_convex_hull_nlogn_2d ( nuv, uv_xy, nw, w )
!
!  Collect the points.
!
  w_xy(1:2,1:nw) = uv_xy(1:2,w(1:nw))

  return
end
function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    This routine implements the recursion
!
!      seed = ( 16807 * seed )  mod ( 2^31 - 1 )
!      r4_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
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
!    Pierre LEcuyer,
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
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 4 ) r4_uniform_01

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
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
subroutine r82vec_part_quick_a ( n, a, l, r )

!*****************************************************************************80
!
!! R82VEC_PART_QUICK_A reorders a R82VEC as part of a quick sort.
!
!  Discussion:
!
!    A R82VEC is a vector of D2's.
!    Each D2 is of type double precision, with two entries.
!    A R82VEC may be stored as 2 by N array.
!
!    The routine reorders the entries of A.  Using A(1:2,1) as a
!    key, all entries of A that are less than or equal to the key will
!    precede the key, which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      N = 8
!
!      A = ( (2,4), (8,8), (6,2), (0,2), (10,6), (10,0), (0,6), (4,8) )
!
!    Output:
!
!      L = 2, R = 4
!
!      A = ( (0,2), (0,6), (2,4), (8,8), (6,2), (10,6), (10,0), (4,8) )
!             -----------          ----------------------------------
!             LEFT          KEY    RIGHT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input/output, real ( kind = 8 ) A(2,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define the three
!    segments.
!    Let KEY = the input value of A(1:2,1).  Then
!    I <= L                 A(1:2,I) < KEY;
!         L < I < R         A(1:2,I) = KEY;
!                 R <= I    KEY < A(1:2,I).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num,n)
  logical r8vec_eq
  logical r8vec_gt
  logical r8vec_lt
  integer ( kind = 4 ) i
  real ( kind = 8 ) key(dim_num)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R82VEC_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:dim_num) = a(1:dim_num,1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( r8vec_gt ( dim_num, a(1:dim_num,l+1), key(1:dim_num) ) ) then
      r = r - 1
      call r8vec_swap ( dim_num, a(1:dim_num,r), a(1:dim_num,l+1) )
    else if ( r8vec_eq ( dim_num, a(1:dim_num,l+1), key(1:dim_num) ) ) then
      m = m + 1
      call r8vec_swap ( dim_num, a(1:dim_num,m), a(1:dim_num,l+1) )
      l = l + 1
    else if ( r8vec_lt ( dim_num, a(1:dim_num,l+1), key(1:dim_num) ) ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(1:dim_num,i) = a(1:dim_num,i+m)
  end do

  l = l - m

  do i = 1, dim_num
    a(i,l+1:l+m) = key(i)
  end do

  return
end
subroutine r82vec_sort_quick_a ( n, a )

!*****************************************************************************80
!
!! R82VEC_SORT_QUICK_A ascending sorts a R82VEC using quick sort.
!
!  Discussion:
!
!    A R82VEC is a vector of R82's.
!    Each R82 is of type double precision, with two entries.
!    A R82VEC may be stored as 2 by N array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(2,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: level_max = 25
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num,n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n_segment
  integer ( kind = 4 ) rsave(level_max)
  integer ( kind = 4 ) r_segment

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R82VEC_SORT_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  else if ( n == 1 ) then
    return
  end if

  level = 1
  rsave(level) = n + 1
  base = 1
  n_segment = n

  do
!
!  Partition the segment.
!
    call r82vec_part_quick_a ( n_segment, a(1,base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( level_max < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R82VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i8)' ) '  Exceeding recursion maximum of ', level_max
        stop
      end if

      level = level + 1
      n_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < n_segment ) then

      n_segment = n_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        n_segment = rsave(level-1) - rsave(level)
        level = level - 1

        if ( 0 < n_segment ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine r8mat_print ( lda, m, n, a, title )

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
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(LDA,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) m
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do jlo = 1, n, 5
    jhi = min ( jlo + 4, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,5(i7,7x))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,5g14.6)' ) i, a(i,jlo:jhi)
    end do
  end do

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
!    08 November 2000
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

    do i = j+1, n
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
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
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
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
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
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
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

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
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

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8mat2_inverse ( a, b, det )

!*****************************************************************************80
!
!! R8MAT2_INVERSE inverts a 2 by 2 real matrix using Cramer's rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2,2), the matrix to be inverted.
!
!    Output, real ( kind = 8 ) B(2,2), the inverse of the matrix A.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
!    If DET is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
  implicit none

  real ( kind = 8 ) a(2,2)
  real ( kind = 8 ) b(2,2)
  real ( kind = 8 ) det
!
!  Compute the determinant.
!
  det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then

    b(1:2,1:2) = 0.0D+00

    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit formula.
!
  b(1,1) = + a(2,2) / det
  b(1,2) = - a(1,2) / det
  b(2,1) = - a(2,1) / det
  b(2,2) = + a(1,1) / det

  return
end
function r8vec_eq ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_EQ is true if two R8VEC's are equal.
!
!  Discussion:
!
!    A R8VEC is an array of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), two vectors to compare.
!
!    Output, logical R8VEC_EQ, is TRUE if every pair of elements A1(I)
!    and A2(I) are equal, and FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  logical r8vec_eq

  r8vec_eq = ( all ( a1(1:n) == a2(1:n) ) )

  return
end
function r8vec_gt ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_GT == ( A1 > A2 ) for R8VEC's.
!
!  Discussion:
!
!    A R8VEC is an array of double precision real values.
!
!    The comparison is lexicographic.
!
!    A1 > A2  <=>                              A1(1) > A2(1) or
!                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
!                 ...
!                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, logical R8VEC_GT, is TRUE if and only if A1 > A2.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  logical r8vec_gt
  integer ( kind = 4 ) i

  r8vec_gt = .false.

  do i = 1, n

    if ( a2(i) < a1(i) ) then
      r8vec_gt = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      r8vec_gt = .false.
      exit
    end if

  end do

  return
end
function r8vec_lt ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_LT == ( A1 < A2 ) for R8VEC's.
!
!  Discussion:
!
!    A R8VEC is an array of double precision real values.
!
!    The comparison is lexicographic.
!
!    A1 < A2  <=>                              A1(1) < A2(1) or
!                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
!                 ...
!                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, logical R8VEC_LT, is TRUE if and only if A1 < A2.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  logical r8vec_lt
  integer ( kind = 4 ) i

  r8vec_lt = .false.

  do i = 1, n

    if ( a1(i) < a2(i) ) then
      r8vec_lt = .true.
      exit
    else if ( a2(i) < a1(i) ) then
      r8vec_lt = .false.
      exit
    end if

  end do

  return
end
subroutine r8vec_swap ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_SWAP swaps the entries of two R8VECs.
!
!  Discussion:
!
!    A R8VEC is an array of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the arrays.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N), the vectors to swap.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a3(n)

  a3(1:n) = a1(1:n)
  a1(1:n) = a2(1:n)
  a2(1:n) = a3(1:n)

  return
end
subroutine r8vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! R8VEC2_COMPARE compares elements of an R8VEC2.
!
!  Discussion:
!
!    The lexicographic ordering is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), contain the two components
!    of each item.
!
!    Input, integer ( kind = 4 ) I, J, the items to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item I > item J.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(i) > a2(j) ) then
      isgn = +1
    end if

  else if ( a1(i) > a1(j) ) then

    isgn = +1

  end if

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints a pair of real vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,2g14.6)' ) i, a1(i), a2(i)
  end do

  return
end
subroutine r8vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    Each item to be sorted is a pair (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call r8_swap ( a1(i), a1(j) )
      call r8_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call r8vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

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
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) radians_to_degrees

  radians_to_degrees = ( angle / pi ) * 180.0D+00

  return
end
subroutine rect_int_2d ( x1, y1, x2, y2, x3, y3, x4, y4, flag, x5, y5, x6, y6 )

!*****************************************************************************80
!
!! RECT_INT_2D computes the intersection of two rectangles in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the corners of the first segment.
!
!    Input, real ( kind = 8 ) X3, Y3, X4, Y4, the corners of the second segment.
!
!    Output, integer ( kind = 4 ) FLAG, records the results.
!    0, the rectangles do not intersect.
!    1, the intersection is a point.
!    2, the intersection is a line.
!    3, the intersection is a rectangle.
!
!    Output, real ( kind = 8 ) X5, Y5, X6, Y6, the corners of the intersection.
!    If FLAG = 0, X5 = Y5 = X6 = Y6 = 0.
!
  implicit none

  integer ( kind = 4 ) flag
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) x6
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) y5
  real ( kind = 8 ) y6

  x5 = 0.0D+00
  y5 = 0.0D+00
  x6 = 0.0D+00
  y6 = 0.0D+00

  call lines_seg_int_1d ( x1, x2, x3, x4, flag, x5, x6 )

  if ( flag == 0 ) then
    return
  end if

  call lines_seg_int_1d ( y1, y2, y3, y4, flag, y5, y6 )

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, real ( kind = 8 )s, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

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
subroutine triangle_contains_point_2d ( x1, y1, x2, y2, x3, y3, x, y, inside )

!*****************************************************************************80
!
!! TRIANGLE_CONTAINS_POINT_2D finds if a point is inside a triangle in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates of
!    the corners of the triangle.
!
!    Input, real ( kind = 8 ) X, Y, the point to be checked.
!
!    Output, logical INSIDE, is .TRUE. if (X,Y) is inside
!    the triangle or on its boundary, and .FALSE. otherwise.
!
  implicit none

  integer ( kind = 4 ), parameter :: N = 2
  integer ( kind = 4 ), parameter :: NRHS = 1

  real ( kind = 8 ) a(N,N+NRHS)
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  integer ( kind = 4 ) info
  logical inside
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
!
!  Set up the linear system
!
!    ( X2-X1  X3-X1 ) C1  = X-X1
!    ( Y2-Y1  Y3-Y1 ) C2    Y-Y1
!
!  which is satisfied by the barycentric coordinates of (X,Y).
!
  a(1,1) = x2 - x1
  a(1,2) = x3 - x1
  a(1,3) = x - x1

  a(2,1) = y2 - y1
  a(2,2) = y3 - y1
  a(2,3) = y - y1
!
!  Solve the linear system.
!
  call r8mat_solve ( a, N, NRHS, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_CONTAINS_POINT - Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper triangle.'
    stop
  end if

  c1 = a(1,3)
  c2 = a(2,3)
!
!  If the point is in the triangle, its barycentric coordinates
!  must both be nonnegative, and sum to no more than 1.
!
  if ( c1 < 0.0D+00 .or. c2 < 0.0D+00 ) then
    inside = .false.
  else if ( 1.0D+00 < c1 + c2 ) then
    inside = .false.
  else
    inside = .true.
  end if

  return
end
subroutine triangulate_tricolor ( node_num, triang, color )

!*****************************************************************************80
!
!! TRIANGULATE_TRICOLOR three-colors the nodes of a triangulated polygon.
!
!  Discussion:
!
!    While the data in TRIANG must represent a proper triangulation
!    of the polygon, no check for this is made.
!
!    The three-coloring is done in such a way that every triangle
!    of the triangulation has one node of each color.
!
!    This solves the art gallery problem, which asks for the minimum
!    number of guards necessary in a (polygonal) art gallery, so that
!    every point inside the gallery is seen by at least one guard.
!    The answer is that (N/3) (integer division, with truncation) guards
!    are sufficient.  Triangulate the polygon, tricolor the triangulation,
!    choose the color that shows up least often (in case N is not divisible
!    by 3) and place a guard at those points.  Since every triangle has
!    a node of that color, and the triangles make up the gallery, you're done.
!    This assumes that the gallery is a simple polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 47-48.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes in the polygon.
!
!    Input, integer ( kind = 4 ) TRIANG(3,NODE_NUM-2), the triangulation of
!    the polygon.
!
!    Output, integer ( kind = 4 ) COLOR(NODE_NUM), an assignment of the
!    "colors" 1, 2 and 3 to the triangles in such a way that every triangle
!    in the triangulation has one node of each color.
!
  implicit none

  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) color(node_num)
  integer ( kind = 4 ) color1
  integer ( kind = 4 ) color2
  integer ( kind = 4 ) color3
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  integer ( kind = 4 ) node3
  integer ( kind = 4 ) stack(6*(node_num-3))
  integer ( kind = 4 ) stack_max
  integer ( kind = 4 ) stack_num
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2
  integer ( kind = 4 ) triang(3,node_num-2)

  stack_max = 6 * ( node_num - 3 )
  t1 = 1

  node1 = triang(1,t1)
  node2 = triang(2,t1)
  node3 = triang(3,t1)

  color1 = 1
  color2 = 2
  color3 = 3

  color(node1) = color1
  color(node2) = color2
  color(node3) = color3

  stack_num = 0

  call triangulate_color_push ( t1, node1, color1, node2, color2, color3, &
    stack_max, stack_num, stack )

  call triangulate_color_push ( t1, node2, color2, node3, color3, color1, &
    stack_max, stack_num, stack )

  call triangulate_color_push ( t1, node3, color3, node1, color1, color2, &
    stack_max, stack_num, stack )

  do while ( 0 < stack_num )

    call triangulate_color_pop ( t1, node1, color1, node2, color2, color3, &
      stack_max, stack_num, stack )

    call triangulate_common_edge ( node_num-2, triang, node1, node2, t1, &
      t2, node3 )

    if ( t2 /= 0 ) then

      color(node3) = color3

      call triangulate_color_push ( t2, node1, color1, node3, color3, color2, &
        stack_max, stack_num, stack )

      call triangulate_color_push ( t2, node3, color3, node2, color2, color1, &
        stack_max, stack_num, stack )

    end if

  end do

  return
end
subroutine triangulate_color_push ( t, node1, color1, node2, color2, color3, &
  stack_max, stack_num, stack )

!*****************************************************************************80
!
!! TRIANGULATE_COLOR_PUSH pushes a side of a colored triangle onto the stack.
!
!  Discussion:
!
!    This is a utility routine used by TRIANGULATE_TRICOLOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) T, the triangle that the side belongs to.
!
!    Input, integer ( kind = 4 ) NODE1, COLOR1, the starting node of the
!    edge, and its color.
!
!    Input, integer ( kind = 4 ) NODE2, COLOR2, the end node of the edge,
!    and its color.
!
!    Input, integer ( kind = 4 ) COLOR3, the remaining color.
!
!    Input, integer ( kind = 4 ) STACK_MAX, the maximum size of the stack.
!
!    Input/output, integer ( kind = 4 ) STACK_NUM, the current size of the stack.
!
!    Input/output, integer ( kind = 4 ) STACK(STACK_MAX), the stack.
!
  implicit none

  integer ( kind = 4 ) stack_max

  integer ( kind = 4 ) color1
  integer ( kind = 4 ) color2
  integer ( kind = 4 ) color3
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  integer ( kind = 4 ) stack(stack_max)
  integer ( kind = 4 ) stack_num
  integer ( kind = 4 ) t

  stack_num = stack_num + 1
  stack(stack_num) = t
  stack_num = stack_num + 1
  stack(stack_num) = node1
  stack_num = stack_num + 1
  stack(stack_num) = color1
  stack_num = stack_num + 1
  stack(stack_num) = node2
  stack_num = stack_num + 1
  stack(stack_num) = color2
  stack_num = stack_num + 1
  stack(stack_num) = color3

  return
end
subroutine triangulate_color_pop ( t, node1, color1, node2, color2, color3, &
  stack_max, stack_num, stack )

!*****************************************************************************80
!
!! TRIANGULATE_COLOR_POP pops a side of a colored triangle from the stack.
!
!  Discussion:
!
!    This is a utility routine used by TRIANGULATE_TRICOLOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) T, the triangle that the side belongs to.
!
!    Output, integer ( kind = 4 ) NODE1, COLOR1, the starting node of the
!    edge, and its color.
!
!    Output, integer ( kind = 4 ) NODE2, COLOR2, the end node of the edge,
!    and its color.
!
!    Output, integer ( kind = 4 ) COLOR3, the remaining color.
!
!    Input, integer ( kind = 4 ) STACK_MAX, the maximum size of the stack.
!
!    Input/output, integer ( kind = 4 ) STACK_NUM, the current size of the stack.
!
!    Input/output, integer ( kind = 4 ) STACK(STACK_MAX), the stack.
!
  implicit none

  integer ( kind = 4 ) stack_max

  integer ( kind = 4 ) color1
  integer ( kind = 4 ) color2
  integer ( kind = 4 ) color3
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  integer ( kind = 4 ) stack(stack_max)
  integer ( kind = 4 ) stack_num
  integer ( kind = 4 ) t

  color3 = stack(stack_num)
  stack_num = stack_num - 1
  color2 = stack(stack_num)
  stack_num = stack_num - 1
  node2 = stack(stack_num)
  stack_num = stack_num - 1
  color1 = stack(stack_num)
  stack_num = stack_num - 1
  node1 = stack(stack_num)
  stack_num = stack_num - 1
  t = stack(stack_num)
  stack_num = stack_num - 1

  return
end
subroutine triangulate_common_edge ( triang_num, triang, node1, node2, t1, &
  t2, node3 )

!*****************************************************************************80
!
!! TRIANGULATE_COMMON_EDGE seeks the other triangle that shares an edge.
!
!  Discussion:
!
!    This is a utility routine used by TRIANGULATE_TRICOLOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TRIANG_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANG(3,TRIANG_NUM), the triangulation
!    of the polygon.
!
!    Input, integer ( kind = 4 ) NODE1, NODE2, the starting and ending nodes
!    of an edge.
!
!    Input, integer ( kind = 4 ) T1, the triangle to which the edge
!    (NODE1,NODE2) belongs.
!
!    Output, integer ( kind = 4 ) T2, is 0 if there is no triangle containing
!    the matching edge (NODE2,NODE1), or it is the index of a triangle
!    containing the edge (NODE2,NODE1).
!
!    Output, integer ( kind = 4 ) NODE3, the other node in triangle T2.
!
  implicit none

  integer ( kind = 4 ) triang_num

  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  integer ( kind = 4 ) node3
  integer ( kind = 4 ) t
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2
  integer ( kind = 4 ) triang(3,triang_num)

  t2 = 0
  node3 = 0

  do t = 1, triang_num
    do j1 = 1, 3
      j2 = i4_wrap ( j1+1, 1, 3 )
      if ( triang(j2,t) == node1 .and. triang(j1,t) == node2 ) then
        t2 = t
        j3 = i4_wrap ( j1+2, 1, 3 )
        node3 = triang(j3,t)
        return
      end if
    end do
  end do

  return
end
subroutine triangulation_boundary_count ( point_num, tri_num, bound_num )

!*****************************************************************************80
!
!! TRIANGULATION_BOUNDARY_COUNT returns the number of boundary edges.
!
!  Discussion:
!
!    We assume we are given information about a legal, maximal triangulation
!    of a set of points in the plane.  In a maximal triangulation, no more
!    edges can be added without crossing an existing edge.
!
!    Given the number of points and triangles, we are going to apply
!    Euler's formula to determine the number of edges that lie on the
!    convex hull of the set of points.
!
!    The number of faces, including the infinite face, is TRI_NUM + 1.
!
!    Let BOUND_NUM denote the number of edges on the boundary of the convex
!    hull.  Each of the TRI_NUM triangles uses three edges.  Every edge
!    occurs in two different faces, so the number of edges must be
!    ( 3 * TRI_NUM + BOUND_NUM ) / 2.
!
!    The number of points is POINT_NUM.
!
!    Euler's formula asserts that, for a simple connected figure in the
!    plane with no edge crossings, POINT_NUM points, EDGE_NUM edges and
!    FACE_NUM faces:
!
!      POINT_NUM - EDGE_NUM + FACE_NUM = 2
!
!    In our context, this becomes
!
!      POINT_NUM - ( 3 * TRI_NUM + BOUND_NUM ) / 2 + TRI_NUM + 1 = 2
!
!    or
!
!      BOUND_NUM = 2 * POINT_NUM - TRI_NUM - 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry, Section 9.1,
!    Springer, 2000.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Output, integer ( kind = 4 ) BOUND_NUM, the number of edges that lie
!    on the convex hull of the triangulation.
!
  implicit none

  integer ( kind = 4 ) bound_num
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) tri_num

  bound_num = 2 * point_num - tri_num - 2

  return
end
