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
  c2 = max ( c2, - 1.0D+00 )
  c2 = min ( c2, + 1.0D+00 )

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
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
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
!    Input, character ( len = * ) TITLE, a title.
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
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2009
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
!    Input, character ( len = * ) TITLE, a title.
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

      write ( *, '(i5,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine icos_shape ( point_num, edge_num, face_num, face_order_max, &
  point_coord, edge_point, face_order, face_point )

!*****************************************************************************80
!
!! ICOS_SHAPE describes an icosahedron.
!
!  Discussion:
!
!    The input data required for this routine can be retrieved from ICOS_SIZE.
!
!    The vertices lie on the unit sphere.
!
!    The dual of an icosahedron is a dodecahedron.
!
!    The data has been rearranged from a previous assignment.  
!    The STRIPACK program refuses to triangulate data if the first
!    three nodes are "collinear" on the sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points (12).
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges (30).
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces (20).
!
!    Input, integer ( kind = 4 ) FACE_ORDER_MAX, the maximum number of 
!    vertices per face (3).
!
!    Output, real ( kind = 8 ) POINT_COORD(3,POINT_NUM), the points.
!
!    Output, integer ( kind = 4 ) EDGE_POINT(2,EDGE_NUM), the points that 
!    make up each edge, listed in ascending order of their indexes.
!
!    Output, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Output, integer ( kind = 4 ) FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
!    FACE_POINT(I,J) is the index of the I-th point in the J-th face.  The
!    points are listed in the counter clockwise direction defined
!    by the outward normal at the face.  The nodes of each face are ordered 
!    so that the lowest index occurs first.  The faces are then sorted by
!    nodes.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ), parameter :: edge_order = 2
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) edge_point(edge_order,edge_num)
  integer ( kind = 4 ) face_order(face_num)
  integer ( kind = 4 ) face_point(face_order_max,face_num)
  real ( kind = 8 ) phi
  real ( kind = 8 ) point_coord(3,point_num)
  real ( kind = 8 ) z
!
!  Set the point coordinates.
!
  phi = 0.5D+00 * ( sqrt ( 5.0D+00 ) + 1.0D+00 )

  a = phi / sqrt ( 1.0D+00 + phi * phi )
  b = 1.0D+00 / sqrt ( 1.0D+00 + phi * phi )
  z = 0.0D+00
!
!  A*A + B*B + Z*Z = 1.
!
  point_coord(1:3,1:point_num) = reshape ( (/ &
      a,  b,  z, &
      a, -b,  z, &
      b,  z,  a, &
      b,  z, -a, &
      z,  a,  b, &
      z,  a, -b, &
      z, -a,  b, &
      z, -a, -b, &
     -b,  z,  a, &
     -b,  z, -a, &
     -a,  b,  z, &
     -a, -b,  z /), (/ 3, point_num /) )
!
!  Set the edges.
!
  edge_point(1:edge_order,1:edge_num) = reshape ( (/ &
     1,  2, &
     1,  3, &
     1,  4, &
     1,  5, &
     1,  6, &
     2,  3, &
     2,  4, &
     2,  7, &
     2,  8, &
     3,  5, &
     3,  7, &
     3,  9, &
     4,  6, &
     4,  8, &
     4, 10, &
     5,  6, &
     5,  9, &
     5, 11, &
     6, 10, &
     6, 11, &
     7,  8, &
     7,  9, &
     7, 12, &
     8, 10, &
     8, 12, &
     9, 11, &
     9, 12, &
    10, 11, &
    10, 12, &
    11, 12 /), (/ edge_order, edge_num /) )
!
!  Set the face orders.
!
  face_order(1:face_num) = (/ &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /)
!
!  Set the faces.
!
  face_point(1:face_order_max,1:face_num) = reshape ( (/ &
     1,  2,  4, &
     1,  3,  2, &
     1,  4,  6, &
     1,  5,  3, &
     1,  6,  5, &
     2,  3,  7, &
     2,  7,  8, &
     2,  8,  4, &
     3,  5,  9, &
     3,  9,  7, &
     4,  8, 10, &
     4, 10,  6, &
     5,  6, 11, &
     5, 11,  9, &
     6, 10, 11, &
     7,  9, 12, &
     7, 12,  8, &
     8, 12, 10, &
     9, 11, 12, &
    10, 12, 11 /), (/ face_order_max, face_num /) )

  return
end
subroutine icos_num ( point_num, edge_num, face_num, face_order_max )

!*****************************************************************************80
!
!! ICOS_NUM gives "sizes" for an icosahedron in 3D.
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
!    Output, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Output, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) FACE_ORDER_MAX, the maximum order of any face.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ) point_num

  point_num = 12
  edge_num = 30
  face_num = 20
  face_order_max = 3

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
!! R8MAT_TRANSPOSE_PRINT prints a R8MAT, transposed.
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
!    Input, character ( len = * ) TITLE, a title.
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
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2009
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
  character ( len = * )  title

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

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
function r8vec_diff_norm ( n, a, b )

!*****************************************************************************80
!
!! R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), B(N), the vectors
!
!    Output, real ( kind = 8 ) R8VEC_DIFF_NORM, the L2 norm of A - B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) r8vec_diff_norm

  r8vec_diff_norm = sqrt ( sum ( ( a(1:n) - b(1:n) )**2 ) )

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
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine sphere_cubed_ijk_to_xyz ( n, i, j, k, xyz )

!*****************************************************************************80
!
!! SPHERE_CUBED_IJK_TO_XYZ: cubed sphere IJK to XYZ coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of sections into which each 
!    face of the cube is to be divided.
!
!    Input, integer ( kind = 4 ) I, J, K, indices between 0 and N.  Normally,
!    at least one of the indices should have the value 0 or N.
!
!    Output, real ( kind = 8 ) XYZ(3), coordinates of the point.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) xc
  real ( kind = 8 ) xyz(3)
  real ( kind = 8 ) xyzn
  real ( kind = 8 ) yc
  real ( kind = 8 ) zc

  if ( i == 0 ) then
    xc = -1.0D+00
  else if ( i == n ) then
    xc = +1.0D+00
  else
    xc = tan ( real ( 2 * i - n, kind = 8 ) * 0.25D+00 * pi &
      / real ( n, kind = 8 ) )
  end if

  if ( j == 0 ) then
    yc = -1.0D+00
  else if ( j == n ) then
    yc = +1.0D+00
  else
    yc = tan ( real ( 2 * j - n, kind = 8 ) * 0.25D+00 * pi &
      / real ( n, kind = 8 ) )
  end if

  if ( k == 0 ) then
    zc = -1.0D+00
  else if ( k == n ) then
    zc = +1.0D+00
  else
    zc = tan ( real ( 2 * k - n, kind = 8 ) * 0.25D+00 * pi &
      / real ( n, kind = 8 ) )
  end if

  xyzn = sqrt ( xc ** 2 + yc ** 2 + zc ** 2 )

  xyz(1) = xc / xyzn
  xyz(2) = yc / xyzn
  xyz(3) = zc / xyzn

  return
end
subroutine sphere_cubed_line_num ( n, line_num )

!*****************************************************************************80
!
!! SPHERE_CUBED_LINE_NUM counts lines on a cubed sphere grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of sections into which each 
!    face of the cube is to be divided.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
  implicit none

  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) n

  line_num = 0
!
!  If N = 1, the corners form 12 lines.
!
  if ( n == 1 ) then
    line_num = 12
    return
!
!  If 1 < N, each of 8 corners connects to three neighboring edges.
!
  else
    line_num = line_num + 8 * 3
  end if
!
!  If 2 < N, then each of the 12 edges includes lines.
!
  if ( 2 < n ) then
    line_num = line_num + 12 * ( n - 2 )
  end if
!
!  Lines that belong to one of the six faces.
!
  if ( 1 < n ) then
    line_num = line_num + 6 * 2 * n * ( n - 1 )
  end if

  return
end
subroutine sphere_cubed_points ( n, ns, xyz )

!*****************************************************************************80
!
!! SPHERE_CUBED_POINTS computes the points on a cubed sphere grid.
!
!  Discussion:
!
!    For a value of N = 3, for instance, each of the 6 cube faces will
!    be divided into 3 sections, so that a single cube face will have
!    (3+1)x(3+1) points:
!
!      X---X---X---X
!      | 1 | 4 | 7 |
!      X---X---X---X
!      | 2 | 5 | 8 |
!      X---X---X---X
!      | 3 | 6 | 9 |
!      X---X---X---X
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of sections into which each 
!    face of the cube is to be divided.
!
!    Input, integer ( kind = 4 ) NS, the number of points.
!
!    Output, real ( kind = 8 ) XYZ(3,NS), distinct points on the unit sphere
!    generated by a cubed sphere grid.
!
  implicit none

  integer ( kind = 4 ) ns

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) xyz(3,ns)

  ns2 = 0
!
!  Bottom full.
!
  call sphere_cubed_points_face ( n, 0, 0, 0, n, n, 0, ns2, xyz )
!
!  To avoid repetition, draw the middles as grids of n-2 x n-1 points.
!
  call sphere_cubed_points_face ( n, 0, 0, 1, 0,   n-1, n-1, ns2, xyz )
  call sphere_cubed_points_face ( n, 0, n, 1, n-1, n,   n-1, ns2, xyz )
  call sphere_cubed_points_face ( n, n, 1, 1, n,   n,   n-1, ns2, xyz )
  call sphere_cubed_points_face ( n, 1, 0, 1, n,   0,   n-1, ns2, xyz )
!
!  Top full.
!
  call sphere_cubed_points_face ( n, 0, 0, n, n, n, n, ns2, xyz )
  
  if ( ns2 /= ns ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_CUBED_POINTS - Fatal error!'
    write ( *, '(a,i8,a)' ) '  Expected to generated NS = ', ns, ' points.'
    write ( *, '(a,i8,a)' ) '  Generated ', ns2, ' points.'
    stop
  end if

  return
end
subroutine sphere_cubed_points_face ( n, i1, j1, k1, i2, j2, k2, ns, xyz )

!*****************************************************************************80
!
!! SPHERE_CUBED_POINTS_FACE: points on one face of a cubed sphere grid.
!
!  Discussion:
!
!    This routine starts with NS = 0, and is called repeatedly to
!    add points for another face.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of sections into which each face 
!    of the cube is to be divided.
!
!    Input, integer ( kind = 4 ) I1, J1, K1, I2, J2, K2, the logical indices, 
!    between 0 and N, of two corners of the face grid.  It is guaranteed that 
!    I1 <= I2, J1 <= J2, and K1 <= K2.  
!
!    Input/output, integer ( kind = 4 ) NS, the number of points.
!
!    Input/output, real XYZ(3,NS), distinct points on the unit sphere
!    generated by a cubed sphere grid.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ns
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) xyz(3,*)
  real ( kind = 8 ) xyzn
  real ( kind = 8 ) xc
  real ( kind = 8 ) yc
  real ( kind = 8 ) zc

  do i = i1, i2

    if ( i1 < i2 ) then
      xc = tan ( real ( 2 * i - n, kind = 8 ) * 0.25D+00 * pi &
        / real ( n, kind = 8 ) )
    else if ( i1 == 0 ) then
      xc = -1.0D+00
    else if ( i1 == n ) then
      xc = +1.0D+00
    else
      xc = 0.0D+00
    end if

    do j = j1, j2

      if ( j1 < j2 ) then
        yc = tan ( real ( 2 * j - n, kind = 8 ) * 0.25D+00 * pi &
          / real ( n, kind = 8 ) )
      else if ( j1 == 0 ) then
        yc = -1.0D+00
      else if ( j1 == n ) then
        yc = +1.0D+00
      else
        yc = 0.0D+00
      end if

      do k = k1, k2

        if ( k1 < k2 ) then
          zc = tan ( real ( 2 * k - n, kind = 8 ) * 0.25D+00 * pi &
            / real ( n, kind = 8 ) );
        else if ( k1 == 0 ) then
          zc = -1.0D+00
        else if ( k1 == n ) then
          zc = +1.0D+00
        else
          zc = 0.0D+00
        end if

        xyzn = sqrt ( xc ** 2 + yc ** 2 + zc ** 2 )

        ns = ns + 1
        xyz(1,ns) = xc / xyzn
        xyz(2,ns) = yc / xyzn
        xyz(3,ns) = zc / xyzn

      end do
    end do
  end do

  return
end
subroutine sphere_cubed_point_num ( n, ns )

!*****************************************************************************80
!
!! SPHERE_CUBED_POINT_NUM counts the points on a cubed sphere grid.
!
!  Discussion:
!
!    For a value of N = 3, for instance, each of the 6 cube faces will
!    be divided into 3 sections, so that a single cube face will have
!    (3+1)x(3+1) points:
!
!      X---X---X---X
!      | 1 | 4 | 7 |
!      X---X---X---X
!      | 2 | 5 | 8 |
!      X---X---X---X
!      | 3 | 6 | 9 |
!      X---X---X---X
!
!    The number of points is simply (N+1)^3 - (N-1)^3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of sections into which 
!    each face of the cube is to be divided.
!
!    Output, integer ( kind = 4 ) NS, the number of points.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ns

  ns = ( n + 1 ) ** 3 - ( n - 1 ) ** 3

  return
end
subroutine sphere_distance_xyz ( xyz1, xyz2, dist )

!*****************************************************************************80
!
!! SPHERE_DISTANCE_XYZ computes great circle distances on a sphere.
!
!  Discussion:
!
!    XYZ coordinates are used.
!
!    We assume the points XYZ1 and XYZ2 lie on the same sphere.
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
  real ( kind = 8 ) r
  real ( kind = 8 ) r8vec_norm
  real ( kind = 8 ) top
  real ( kind = 8 ) xyz1(3)
  real ( kind = 8 ) xyz2(3)

  r = r8vec_norm ( 3, xyz1 )

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

  dist = r * atan2 ( top, bot )

  return
end
subroutine sphere_grid_q4 ( lat_num, long_num, rectangle_node )

!*****************************************************************************80
!
!! SPHERE_GRID_Q4: rectangular grid on a sphere.
!
!  Discussion:
!
!    The point numbering system is the same used in SPHERE_GRIDPOINTS,
!    and that routine may be used to compute the coordinates of the points.
!
!    A sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R * R
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LAT_NUM, the number of "rows" of rectangles to
!    be created.  LAT_NUM must be at least 2. 
!
!    Input, integer ( kind = 4 ) LONG_NUM, the number of "columns" of 
!    rectangles to be created.
!
!    Output, integer ( kind = 4 ) RECTANGLE_NODE(4,LAT_NUM*LONG_NUM), 
!    the indices of the nodes that make up each rectangle.
!
  implicit none

  integer ( kind = 4 ) lat_num
  integer ( kind = 4 ) long_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_max
  integer ( kind = 4 ) n_min
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) nw
  integer ( kind = 4 ) s
  integer ( kind = 4 ) s_max
  integer ( kind = 4 ) s_min
  integer ( kind = 4 ) se
  integer ( kind = 4 ) sw
  integer ( kind = 4 ) rectangle_node(4,lat_num*long_num)
  integer ( kind = 4 ) rectangle_num

  rectangle_num = 0
!
!  The first row.
!
  n = 1

  sw = 2
  se = sw + 1

  s_min = 2
  s_max = long_num + 1

  do j = 1, long_num

    rectangle_num = rectangle_num + 1
    rectangle_node(1:4,rectangle_num) = (/ sw, se, n, n /)

    sw = se

    if ( se == s_max ) then
      se = s_min
    else
      se = se + 1
    end if

  end do
!
!  The intermediate rows.
!
  do i = 2, lat_num - 1

    n_max = s_max
    n_min = s_min

    s_max = s_max + long_num
    s_min = s_min + long_num

    nw = n_min
    ne = nw + 1
    sw = s_min
    se = sw + 1

    do j = 1, long_num

      rectangle_num = rectangle_num + 1
      rectangle_node(1:4,rectangle_num) = (/ sw, se, ne, nw /)

      sw = se
      nw = ne

      if ( se == s_max ) then
        se = s_min
      else
        se = se + 1
      end if

      if ( ne == n_max ) then
        ne = n_min
      else
        ne = ne + 1
      end if

    end do

  end do
!
!  The last row.
!
  n_max = s_max
  n_min = s_min

  s = n_max + 1

  nw = n_min
  ne = nw + 1

  do j = 1, long_num

    rectangle_num = rectangle_num + 1
    rectangle_node(1:4,rectangle_num) = (/ ne, nw, s, s /)

    nw = ne

    if ( ne == n_max ) then
      ne = n_min
    else
      ne = ne + 1
    end if

  end do

  return
end
subroutine sphere_grid_t3 ( lat_num, long_num, triangle_node )

!*****************************************************************************80
!
!! SPHERE_GRID_T3 produces a triangle grid on a sphere.
!
!  Discussion:
!
!    The point numbering system is the same used in SPHERE_GRIDPOINTS,
!    and that routine may be used to compute the coordinates of the points.
!
!    A sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - pc(1:DIM_NUM) )^2 ) = R * R
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LAT_NUM, LONG_NUM, the number of latitude 
!    and longitude lines to draw.  The latitudes do not include the North 
!    and South poles, which will be included automatically, so LAT_NUM = 5, 
!    for instance, will result in points along 7 lines of latitude.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NODE(3,2*(LAT_NUM+1)*LONG_NUM), the
!    triangle vertices.
!
  implicit none

  integer ( kind = 4 ) lat_num
  integer ( kind = 4 ) long_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_max
  integer ( kind = 4 ) n_min
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) nw
  integer ( kind = 4 ) s
  integer ( kind = 4 ) s_max
  integer ( kind = 4 ) s_min
  integer ( kind = 4 ) se
  integer ( kind = 4 ) sw
  integer ( kind = 4 ) triangle_node(3,2*(lat_num+1)*long_num)
  integer ( kind = 4 ) triangle_num

  triangle_num = 0
!
!  The first row.
!
  n = 1

  sw = 2
  se = sw + 1

  s_min = 2
  s_max = long_num + 1

  do j = 0, long_num - 1

    triangle_num = triangle_num + 1
    triangle_node(1:3,triangle_num) = (/ sw, se, n /)

    sw = se

    if ( se == s_max ) then
      se = s_min
    else
      se = se + 1
    end if

  end do
!
!  The intermediate rows.
!
  do i = 1, lat_num

    n_max = s_max
    n_min = s_min

    s_max = s_max + long_num
    s_min = s_min + long_num

    nw = n_min
    ne = nw + 1
    sw = s_min
    se = sw + 1

    do j = 0, long_num - 1

      triangle_num = triangle_num + 1
      triangle_node(1:3,triangle_num) = (/ sw, se, nw /)

      triangle_num = triangle_num + 1
      triangle_node(1:3,triangle_num) = (/ ne, nw, se /)

      sw = se
      nw = ne

      if ( se == s_max ) then
        se = s_min
      else
        se = se + 1
      end if

      if ( ne == n_max ) then
        ne = n_min
      else
        ne = ne + 1
      end if

    end do

  end do
!
!  The last row.
!
  n_max = s_max
  n_min = s_min

  s = n_max + 1

  nw = n_min
  ne = nw + 1

  do j = 0, long_num - 1

    triangle_num = triangle_num + 1
    triangle_node(1:3,triangle_num) = (/ ne, nw, s /)

    nw = ne

    if ( ne == n_max ) then
      ne = n_min
    else
      ne = ne + 1
    end if

  end do

  return
end
subroutine sphere_icos_edge_num ( factor, edge_num )

!*****************************************************************************80
!
!! SPHERE_ICOS_EDGE_NUM sizes an icosahedral grid on a sphere.
!
!  Discussion:
!
!    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
!
!    With FACTOR = 2, each triangle of the icosahedron is subdivided into
!    2x2 subtriangles, resulting in 80 faces, 120 edges, and 
!    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
!
!    With FACTOR = 3, each triangle of the icosahedron is subdivided into
!    3x3 subtriangles, resulting in 180 faces, 270 edges and 
!    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
!
!    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
!    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
!      12 
!    + 20 * 3          * (FACTOR-1) / 2 
!    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Output, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) factor

  if ( factor < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_ICOS_EDGE_NUM - Fatal error!'
    write ( *, '(a)' ) '  Input FACTOR < 1.'
    stop
  end if

  edge_num = 30 * factor * factor

  return
end
subroutine sphere_icos_face_num ( factor, face_num )

!*****************************************************************************80
!
!! SPHERE_ICOS_FACE_NUM sizes an icosahedral grid on a sphere.
!
!  Discussion:
!
!    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
!
!    With FACTOR = 2, each triangle of the icosahedron is subdivided into
!    2x2 subtriangles, resulting in 80 faces, 120 edges, and 
!    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
!
!    With FACTOR = 3, each triangle of the icosahedron is subdivided into
!    3x3 subtriangles, resulting in 180 faces, 270 edges and 
!    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
!
!    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
!    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
!      12 
!    + 20 * 3          * (FACTOR-1) / 2 
!    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of triangles.
!
  implicit none

  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) factor

  if ( factor < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_ICOS_FACE_NUM - Fatal error!'
    write ( *, '(a)' ) '  Input FACTOR < 1.'
    stop
  end if

  face_num = 20 * factor * factor

  return
end
subroutine sphere_icos_point_num ( factor, point_num )

!*****************************************************************************80
!
!! SPHERE_ICOS_POINT_NUM sizes an icosahedral grid on a sphere.
!
!  Discussion:
!
!    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
!
!    With FACTOR = 2, each triangle of the icosahedron is subdivided into
!    2x2 subtriangles, resulting in 80 faces, 120 edges, and 
!    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
!
!    With FACTOR = 3, each triangle of the icosahedron is subdivided into
!    3x3 subtriangles, resulting in 180 faces, 270 edges and 
!    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
!
!    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
!    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
!      12 
!    + 20 * 3          * (FACTOR-1) / 2 
!    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) factor
  integer ( kind = 4 ) point_num

  if ( factor < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_ICOS_POINT_NUM - Fatal error!'
    write ( *, '(a)' ) '  Input FACTOR < 1.'
    stop
  end if

  point_num = 12                                   &
            + 10 * 3              * ( factor - 1 ) &
            + 10 * ( factor - 2 ) * ( factor - 1 )

  return
end
subroutine sphere_icos1_points ( factor, node_num, node_xyz )

!*****************************************************************************80
!
!! SPHERE_ICOS1_POINTS returns icosahedral grid points on a sphere.
!
!  Discussion:
!
!    With FACTOR = 1, the grid has 20 triangular faces and 12 nodes.
!
!    With FACTOR = 2, each triangle of the icosahedron is subdivided into
!    2x2 subtriangles, resulting in 80 faces and 
!    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
!
!    With FACTOR = 3, each triangle of the icosahedron is subdivided into
!    3x3 subtriangles, resulting in 180 faces and 
!    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
!
!    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
!    resulting in 20 * FACTOR * FACTOR faces and
!      12 
!    + 20 * 3          * (FACTOR-1) / 2 
!    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
!
!    This routine uses a simple, but only approximate, method of
!    carrying out the subdivision.  For each spherical triangle of the
!    face, we actually work in the planar triangle defined by those
!    three points.  All subdivisions are done on that planar triangle,
!    and the resulting points are then projected onto the sphere.
!    While these points are equally spaced on the planar triangle,
!    that is only approximately true on the sphere.
!
!    See SPHERE_ICOS2_POINTS for a more accurate method of subdivision
!    that works on the spherical triangle itself.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes, as reported
!    by SPHERE_GRID_ICOS_SIZE.
!
!    Output, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the node coordinates.
!
!  Local Parameters:
!
!    POINT_NUM, EDGE_NUM, FACE_NUM and FACE_ORDER_MAX are counters 
!    associated with the icosahedron, and POINT_COORD, EDGE_POINT, 
!    FACE_ORDER and FACE_POINT are data associated with the icosahedron.
!    We need to refer to this data to generate the grid.
!
!    NODE counts the number of nodes we have generated so far.  At the
!    end of the routine, it should be equal to NODE_NUM.
!
  implicit none

  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: edge_point
  integer ( kind = 4 ) f
  integer ( kind = 4 ) f1
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: face_order
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: face_point
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_norm
  real ( kind = 8 ) node_xyz(3,node_num)
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord
  integer ( kind = 4 ) point_num
  real ( kind = 8 ) r8vec_norm
!
!  Size the icosahedron.
!
  call icos_num ( point_num, edge_num, face_num, face_order_max )
!
!  Set the icosahedron.
!
  allocate ( point_coord(1:3,1:point_num) )
  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )
!
!  Generate the point coordinates.
!
!  A.  Points that are the icosahedral vertices.
!
  node = 0
  node_xyz(1:3,1:point_num) = point_coord(1:3,1:point_num)
!
!  B. Points in the icosahedral edges, at 
!  1/FACTOR, 2/FACTOR, ..., (FACTOR-1)/FACTOR.
!
  node = 12

  do edge = 1, edge_num

    a = edge_point(1,edge)
    b = edge_point(2,edge)

    do f = 1, factor - 1

      node = node + 1

      node_xyz(1:3,node) = &
        ( real ( factor - f, kind = 8 ) * point_coord(1:3,a)   &
        + real (          f, kind = 8 ) * point_coord(1:3,b) ) &
        / real ( factor,     kind = 8 )

      node_norm = r8vec_norm ( 3, node_xyz(1:3,node) )

      node_xyz(1:3,node) = node_xyz(1:3,node) / node_norm

    end do
  end do
!
!  C.  Points in the icosahedral faces.
!
  do face = 1, face_num

    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)

    do f1 = 1, factor - 1
      do f2 = 1, factor - f1 - 1

        node = node + 1

        node_xyz(1:3,node) = &
          ( real ( factor - f1 - f2, kind = 8 ) * point_coord(1:3,a)   &
          + real (          f1,      kind = 8 ) * point_coord(1:3,b)   &
          + real (               f2, kind = 8 ) * point_coord(1:3,c) ) &
          / real ( factor,           kind = 8 )

        node_norm = r8vec_norm ( 3, node_xyz(1:3,node) )

        node_xyz(1:3,node) = node_xyz(1:3,node) / node_norm

      end do
    end do

  end do
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

  return
end
subroutine sphere_icos2_points ( factor, node_num, node_xyz )

!*****************************************************************************80
!
!! SPHERE_ICOS2_POINTS returns icosahedral grid points on a sphere.
!
!  Discussion:
!
!    With FACTOR = 1, the grid has 20 triangular faces and 12 nodes.
!
!    With FACTOR = 2, each triangle of the icosahedron is subdivided into
!    2x2 subtriangles, resulting in 80 faces and 
!    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
!
!    With FACTOR = 3, each triangle of the icosahedron is subdivided into
!    3x3 subtriangles, resulting in 180 faces and 
!    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
!
!    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
!    resulting in 20 * FACTOR * FACTOR faces and
!      12 
!    + 20 * 3          * (FACTOR-1) / 2 
!    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes, as reported
!    by SPHERE_GRID_ICOS_SIZE.
!
!    Output, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the node coordinates.
!
!  Local Parameters:
!
!    POINT_NUM, EDGE_NUM, FACE_NUM and FACE_ORDER_MAX are counters 
!    associated with the icosahedron, and POINT_COORD, EDGE_POINT, 
!    FACE_ORDER and FACE_POINT are data associated with the icosahedron.
!    We need to refer to this data to generate the grid.
!
!    NODE counts the number of nodes we have generated so far.  At the
!    end of the routine, it should be equal to NODE_NUM.
!
  implicit none

  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) a
  real ( kind = 8 ) angle
  real ( kind = 8 ) ab(3)
  real ( kind = 8 ) ac(3)
  real ( kind = 8 ) acn(3)
  real ( kind = 8 ) acp(3)
  integer ( kind = 4 ) b
  real ( kind = 8 ) bn(3)
  real ( kind = 8 ) bp(3)
  integer ( kind = 4 ) c
  real ( kind = 8 ) cn(3)
  real ( kind = 8 ) cp(3)
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: edge_point
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fa
  integer ( kind = 4 ) fbc
  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: face_order
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: face_point
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xyz(3,node_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord
  integer ( kind = 4 ) point_num
  real ( kind = 8 ) r8vec_norm
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta_ab
  real ( kind = 8 ) theta_ac
  real ( kind = 8 ) theta_bc
!
!  Size the icosahedron.
!
  call icos_num ( point_num, edge_num, face_num, face_order_max )
!
!  Set the icosahedron.
!
  allocate ( point_coord(1:3,1:point_num) )
  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )
!
!  Generate the point coordinates.
!
!  A.  Points that are the icosahedral vertices.
!
  node = 0
  node_xyz(1:3,1:point_num) = point_coord(1:3,1:point_num)
!
!  B. Points in the icosahedral edges, at 
!  1/FACTOR, 2/FACTOR, ..., (FACTOR-1)/FACTOR.
!
  node = 12

  do edge = 1, edge_num

    a = edge_point(1,edge)
    b = edge_point(2,edge)
!
!  Determine the "distance" = angle between points A and B.
!
    call sphere_distance_xyz ( point_coord(1:3,a), point_coord(1:3,b), theta )
!
!  Polarize B into BP + BN and normalize BN.
!
    call r8vec_polarize ( 3, point_coord(1:3,b), point_coord(1:3,a), bn, bp )

    bn(1:3) = bn(1:3) / r8vec_norm ( 3, bn )
!
!  March from A to B, by taking equally spaced angles from 0 to THETA.
!  F = 0      => ANGLE = 0     => A
!  F = FACTOR => ANGLE = THETA => B
!
    do f = 1, factor - 1

      node = node + 1

      angle = ( real ( f, kind = 8 ) * theta ) / real ( factor, kind = 8 )
     
      node_xyz(1:3,node) = cos ( angle ) * point_coord(1:3,a) &
                         + sin ( angle ) * bn(1:3)

    end do
  end do
!
!  C.  Points in the icosahedral faces.
!
  do face = 1, face_num

    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)
!
!  Determine the "distance" = angle between points A and B, A and C.
!
    call sphere_distance_xyz ( point_coord(1:3,a), point_coord(1:3,b), &
      theta_ab )

    call sphere_distance_xyz ( point_coord(1:3,a), point_coord(1:3,c), &
      theta_ac )
!
!  Polarize B = BP + BN and normalize BN, C = CP + CN, and normalize CN.
!
    call r8vec_polarize ( 3, point_coord(1:3,b), point_coord(1:3,a), bn, bp )
    bn(1:3) = bn(1:3) / r8vec_norm ( 3, bn )

    call r8vec_polarize ( 3, point_coord(1:3,c), point_coord(1:3,a), cn, cp )
    cn(1:3) = cn(1:3) / r8vec_norm ( 3, cn )
!
!  March AB from A to B:
!    FA = 0      => ANGLE = 0        => AB = A
!    FA = FACTOR => ANGLE = THETA_AB => AB = B
!
!  March AC from A to C:
!    FA = 0      => ANGLE = 0        => AC = A
!    FA = FACTOR => ANGLE = THETA_AC => AC = C
!
    do fa = 2, factor - 1
!
!  Determine points AB and AC that use cos ( FA / FACTOR ) of A 
!  and cos ( ( FACTOR - FA ) / FACTOR ) of B or C.
!
      angle = ( real ( fa, kind = 8 ) * theta_ab ) / real ( factor, kind = 8 )
      ab(1:3) = cos ( angle ) * point_coord(1:3,a) + sin ( angle ) * bn(1:3)

      angle = ( real ( fa, kind = 8 ) * theta_ac ) / real ( factor, kind = 8 )
      ac(1:3) = cos ( angle ) * point_coord(1:3,a) + sin ( angle ) * cn(1:3)
!
!  Determine the "distance" = angle between points AB and AC.
!
      call sphere_distance_xyz ( ab(1:3), ac(1:3), theta_bc )
!
!  Polarize AC into ACP + ACN and normalize ACN.
!
      call r8vec_polarize ( 3, ac, ab, acn, acp )
      acn(1:3) = acn(1:3) / r8vec_norm ( 3, acn )
!
!  The interval between AB and AC is broken into FA intervals.
!  Go from 1 to FA - 1.
!
      do fbc = 1, fa - 1

        node = node + 1

        angle = ( real ( fbc, kind = 8 ) * theta_bc ) &
                / real ( fa, kind = 8 )
     
        node_xyz(1:3,node) = cos ( angle ) * ab(1:3) &
                           + sin ( angle ) * acn(1:3)

      end do
    end do

  end do
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

  return
end
subroutine sphere_line_project ( r, pc, n, p, maxpnt2, n2, pp, theta_min, &
  theta_max )

!*****************************************************************************80
!
!! SPHERE_LINE_PROJECT projects a line onto a sphere.
!
!  Discussion:
!
!    A sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R * R
!
!    The line to be projected is specified as a sequence of points.
!    If two successive points subtend a small angle, then the second
!    point is essentially dropped.  If two successive points subtend
!    a large angle, then intermediate points are inserted, so that
!    the projected line stays closer to the sphere.
!
!    Note that if any P coincides with the center of the sphere, then
!    its projection is mathematically undefined.  PP will
!    be returned as the center.
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
!    Input, real ( kind = 8 ) R, the radius of the sphere.  If R is
!    zero, PP will be returned as the pc, and if R is
!    negative, points will end up diametrically opposite from where
!    you would expect them for a positive R.
!
!    Input, real ( kind = 8 ) PC(3), the center of the sphere.
!
!    Input, integer ( kind = 4 ) N, the number of points on the line that is
!    to be projected.
!
!    Input, real ( kind = 8 ) P(3,N), the coordinates of
!    the points on the line that is to be projected.
!
!    Input, integer ( kind = 4 ) MAXPNT2, the maximum number of points on the
!    projected line.  Even if the routine thinks that more points are needed,
!    no more than MAXPNT2 will be generated.
!
!    Output, integer ( kind = 4 ) N2, the number of points on the projected
!    line.  N2 can be zero, if the line has an angular projection of less
!    than THETA_MIN radians.
!
!    Output, real ( kind = 8 ) PP(3,N2), the coordinates
!    of the points representing the projected line.  These points lie on the
!    sphere.  Successive points are separated by at least THETA_MIN
!    radians, and by no more than THETA_MAX radians.
!
!    Input, real ( kind = 8 ) THETA_MIN, THETA_MAX, the minimum and maximum
!    angular projections allowed between successive projected points.
!    If two successive points on the original line have projections
!    separated by more than THETA_MAX radians, then intermediate points
!    will be inserted, in an attempt to keep the line closer to the
!    sphere.  If two successive points are separated by less than
!    THETA_MIN radians, then the second point is dropped, and the
!    line from the first point to the next point is considered.
!
  implicit none

  integer ( kind = 4 ) maxpnt2
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) ang3d
  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) dot
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nfill
  integer ( kind = 4 ) n2
  real ( kind = 8 ) p(3,n)
  real ( kind = 8 ) p1(3)
  real ( kind = 8 ) p2(3)
  real ( kind = 8 ) pc(3)
  real ( kind = 8 ) pd(3)
  real ( kind = 8 ) pp(3,maxpnt2)
  real ( kind = 8 ) r
  real ( kind = 8 ) r8vec_diff_norm
  real ( kind = 8 ) r8vec_norm
  real ( kind = 8 ) theta_max
  real ( kind = 8 ) theta_min
  real ( kind = 8 ) tnorm
!
!  Check the input.
!
  if ( r == 0.0D+00 ) then
    n2 = 0
    return
  end if

  p1(1:3) = pc(1:3)
  p2(1:3) = pc(1:3)

  n2 = 0

  do i = 1, n

    if ( all ( p(1:3,i) == pc(1:3) ) ) then

    else

      p1(1:3) = p2(1:3)

      alpha = r8vec_diff_norm ( 3, p(1:3,i), pc(1:3) )

      p2(1:3) = pc(1:3) + r * ( p(1:3,i) - pc(1:3) ) / alpha
!
!  If we haven't gotten any points yet, take this point as our start.
!
      if ( n2 == 0 ) then

        n2 = n2 + 1
        pp(1:3,n2) = p2(1:3)
!
!  Compute the angular projection of P1 to P2.
!
      else if ( 1 <= n2 ) then

        dot = sum ( ( p1(1:3) - pc(1:3) ) * ( p2(1:3) - pc(1:3) ) )

        ang3d = arc_cosine (  dot / ( r * r ) )
!
!  If the angle is at least THETA_MIN, (or it's the last point),
!  then we will draw a line segment.
!
        if ( theta_min < abs ( ang3d ) .or. i == n ) then
!
!  Now we check to see if the line segment is too long.
!
          if ( theta_max < abs ( ang3d ) ) then

            nfill = int ( abs ( ang3d ) / theta_max )

            do j = 1, nfill-1

              pd(1:3) = &
                ( real ( nfill - j, kind = 8 ) * ( p1(1:3) - pc(1:3) ) &
                + real (         j, kind = 8 ) * ( p2(1:3) - pc(1:3) ) )

              tnorm = r8vec_norm ( 3, pd )

              if ( tnorm /= 0.0D+00 ) then
                pd(1:3) = pc(1:3) + r * pd(1:3) / tnorm
                n2 = n2 + 1
                pp(1:3,n2) = pd(1:3)
              end if

            end do

          end if
!
!  Now tack on the projection of point 2.
!
          n2 = n2 + 1
          pp(1:3,n2) = p2(1:3)

        end if

      end if

    end if

  end do

  return
end
subroutine sphere_ll_lines ( lat_num, long_num, line_num, line )

!*****************************************************************************80
!
!! SPHERE_LL_LINES produces lines for a latitude/longitude grid.
!
!  Discussion:
!
!    The point numbering system is the same used in SPHERE_LL_POINTS,
!    and that routine may be used to compute the coordinates of the points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LAT_NUM, LONG_NUM, the number of latitude and
!    longitude lines to draw.  The latitudes do not include the North and South
!    poles, which will be included automatically, so LAT_NUM = 5, for instance,
!    will result in points along 7 lines of latitude.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of grid lines.
!
!    Output, integer ( kind = 4 ) LINE(2,LINE_NUM), contains pairs of point 
!    indices for line segments that make up the grid.
!
  implicit none

  integer ( kind = 4 ) line_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lat_num
  integer ( kind = 4 ) l
  integer ( kind = 4 ) line(2,line_num)
  integer ( kind = 4 ) long_num
  integer ( kind = 4 ) new
  integer ( kind = 4 ) newcol
  integer ( kind = 4 ) old

  l = 0
!
!  "Vertical" lines.
!
  do j = 0, long_num - 1

    old = 1
    new = j + 2

    l = l + 1
    line(1:2,l) = (/ old, new /)

    do i = 1, lat_num - 1

      old = new
      new = old + long_num

      l = l + 1
      line(1:2,l) = (/ old, new /)

    end do

    old = new

    l = l + 1
    line(1:2,l) = (/ old, 1 + lat_num * long_num + 1 /)

  end do
!
!  "Horizontal" lines.
!
  do i = 1, lat_num

    new = 1 + ( i - 1 ) * long_num + 1

    do j = 0, long_num - 2
      old = new
      new = old + 1
      l = l + 1
      line(1:2,l) = (/ old, new /)
    end do

    old = new
    new = 1 + ( i - 1 ) * long_num + 1
    l = l + 1
    line(1:2,l) = (/ old, new /)

  end do
!
!  "Diagonal" lines.
!
  do j = 0, long_num - 1

    old = 1
    new = j + 2
    newcol = j

    do i = 1, lat_num - 1

      old = new
      new = old + long_num + 1

      newcol = newcol + 1
      if ( long_num - 1 < newcol ) then
        newcol = 0
        new = new - long_num
      end if

      l = l + 1
      line(1:2,l) = (/ old, new /)

    end do

  end do

  return
end
subroutine sphere_ll_line_num ( lat_num, long_num, line_num )

!*****************************************************************************80
!
!! SPHERE_LL_LINE_NUM counts lines for a latitude/longitude grid.
!
!  Discussion:
!
!    The number returned is the number of pairs of points to be connected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LAT_NUM, LONG_NUM, the number of latitude and
!    longitude lines to draw.  The latitudes do not include the North and South
!    poles, which will be included automatically, so LAT_NUM = 5, for instance,
!    will result in points along 7 lines of latitude.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of grid lines.
!
  implicit none

  integer ( kind = 4 ) lat_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) long_num

  line_num = long_num * ( lat_num + 1 ) &
           + lat_num * long_num &
           + long_num * ( lat_num - 1 )

  return
end
subroutine sphere_ll_points ( r, pc, lat_num, long_num, point_num, p )

!*****************************************************************************80
!
!! SPHERE_LL_POINTS produces points for a latitude/longitude grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) PC(3), the center of the sphere.
!
!    Input, integer ( kind = 4 ) LAT_NUM, LONG_NUM, the number of latitude 
!    and longitude lines to draw.  The latitudes do not include the North and 
!    South poles, which will be included automatically, so LAT_NUM = 5, for 
!    instance, will result in points along 7 lines of latitude.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Output, real ( kind = 8 ) P(3,POINT_NUM), the grid points.
!
  implicit none

  integer ( kind = 4 ) lat_num
  integer ( kind = 4 ) long_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) lat
  integer ( kind = 4 ) long
  integer ( kind = 4 ) n
  real ( kind = 8 ) p(3,point_num)
  real ( kind = 8 ) pc(3)
  real ( kind = 8 ) phi
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) theta

  n = 0
!
!  The north pole.
!
  theta = 0.0D+00
  phi = 0.0D+00
  n = n + 1
  p(1,n) = pc(1) + r * sin ( phi ) * cos ( theta )
  p(2,n) = pc(2) + r * sin ( phi ) * sin ( theta )
  p(3,n) = pc(3) + r * cos ( phi )
!
!  Do each intermediate ring of latitude.
!
  do lat = 1, lat_num

    phi = real ( lat,         kind = 8 ) * pi &
        / real ( lat_num + 1, kind = 8 )
!
!  Along that ring of latitude, compute points at various longitudes.
!
    do long = 0, long_num - 1

      theta = real ( long,     kind = 8 ) * 2.0D+00 * pi &
            / real ( long_num, kind = 8 )

      n = n + 1
      p(1,n) = pc(1) + r * sin ( phi ) * cos ( theta )
      p(2,n) = pc(2) + r * sin ( phi ) * sin ( theta )
      p(3,n) = pc(3) + r * cos ( phi )

    end do
  end do
!
!  The south pole.
!
  theta = 0.0D+00
  phi = pi
  n = n + 1
  p(1,n) = pc(1) + r * sin ( phi ) * cos ( theta )
  p(2,n) = pc(2) + r * sin ( phi ) * sin ( theta )
  p(3,n) = pc(3) + r * cos ( phi )

  return
end
subroutine sphere_ll_point_num ( lat_num, long_num, point_num )

!*****************************************************************************80
!
!! SPHERE_LL_POINT_NUM counts points for a latitude/longitude grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LAT_NUM, LONG_NUM, the number of latitude 
!    and longitude lines to draw.  The latitudes do not include the North and 
!    South poles, which will be included automatically, so LAT_NUM = 5, for 
!    instance, will result in points along 7 lines of latitude.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the number of grid points.
!
  implicit none

  integer ( kind = 4 ) lat_num
  integer ( kind = 4 ) long_num
  integer ( kind = 4 ) point_num

  point_num = 2 + lat_num * long_num

  return
end
subroutine sphere_llq_lines ( lat_num, long_num, line_num, line )

!*****************************************************************************80
!
!! SPHERE_LLQ_LINES: latitude/longitude quadrilateral grid lines.
!
!  Discussion:
!
!    The point numbering system is the same used in SPHERE_LL_POINTS,
!    and that routine may be used to compute the coordinates of the points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LAT_NUM, LONG_NUM, the number of latitude and
!    longitude lines to draw.  The latitudes do not include the North and South
!    poles, which will be included automatically, so LAT_NUM = 5, for instance,
!    will result in points along 7 lines of latitude.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of grid lines.
!
!    Output, integer ( kind = 4 ) LINE(2,LINE_NUM), contains pairs of point 
!    indices for line segments that make up the grid.
!
  implicit none

  integer ( kind = 4 ) line_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lat_num
  integer ( kind = 4 ) l
  integer ( kind = 4 ) line(2,line_num)
  integer ( kind = 4 ) long_num
  integer ( kind = 4 ) new
  integer ( kind = 4 ) newcol
  integer ( kind = 4 ) old

  l = 0
!
!  "Vertical" lines.
!
  do j = 0, long_num - 1

    old = 1
    new = j + 2

    l = l + 1
    line(1:2,l) = (/ old, new /)

    do i = 1, lat_num - 1

      old = new
      new = old + long_num

      l = l + 1
      line(1:2,l) = (/ old, new /)

    end do

    old = new

    l = l + 1
    line(1:2,l) = (/ old, 1 + lat_num * long_num + 1 /)

  end do
!
!  "Horizontal" lines.
!
  do i = 1, lat_num

    new = 1 + ( i - 1 ) * long_num + 1

    do j = 0, long_num - 2
      old = new
      new = old + 1
      l = l + 1
      line(1:2,l) = (/ old, new /)
    end do

    old = new
    new = 1 + ( i - 1 ) * long_num + 1
    l = l + 1
    line(1:2,l) = (/ old, new /)

  end do

  return
end
subroutine sphere_llq_line_num ( lat_num, long_num, line_num )

!*****************************************************************************80
!
!! SPHERE_LLQ_LINE_NUM counts lines for a latitude/longitude quadrilateral grid.
!
!  Discussion:
!
!    The number returned is the number of pairs of points to be connected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LAT_NUM, LONG_NUM, the number of latitude and
!    longitude lines to draw.  The latitudes do not include the North and South
!    poles, which will be included automatically, so LAT_NUM = 5, for instance,
!    will result in points along 7 lines of latitude.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of grid lines.
!
  implicit none

  integer ( kind = 4 ) lat_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) long_num

  line_num = long_num * ( lat_num + 1 ) &
           + lat_num * long_num

  return
end
subroutine sphere_spiralpoints ( r, pc, n, p )

!*****************************************************************************80
!
!! SPHERE_SPIRALPOINTS: spiral points on a sphere.
!
!  Discussion:
!
!    The points should be arranged on the sphere in a pleasing design.
!
!    A sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R * R
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Saff, Arno Kuijlaars,
!    Distributing Many Points on a Sphere,
!    The Mathematical Intelligencer,
!    Volume 19, Number 1, 1997, pages 5-11.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) PC(3), the center of the sphere.
!
!    Input, integer ( kind = 4 ) N, the number of points to create.
!
!    Output, real ( kind = 8 ) P(3,N), the grid points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cosphi
  integer ( kind = 4 ) i
  real ( kind = 8 ) p(3,n)
  real ( kind = 8 ) pc(3)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) sinphi
  real ( kind = 8 ) theta

  do i = 1, n

    cosphi = ( real ( n - i,     kind = 8 ) * ( -1.0D+00 ) &
             + real (     i - 1, kind = 8 ) * ( +1.0D+00 ) ) &
             / real ( n     - 1, kind = 8 )

    sinphi = sqrt ( 1.0D+00 - cosphi * cosphi )

    if ( i == 1 .or. i == n ) then
      theta = 0.0D+00
    else
      theta = theta + 3.6D+00 / ( sinphi * sqrt ( real ( n, kind = 8 ) ) )
      theta = mod ( theta, 2.0D+00 * pi )
    end if

    p(1,i) = pc(1) + r * sinphi * cos ( theta )
    p(2,i) = pc(2) + r * sinphi * sin ( theta )
    p(3,i) = pc(3) + r * cosphi

  end do

  return
end
subroutine sphere_unit_sample ( n, seed, x )

!*****************************************************************************80
!
!! SPHERE_UNIT_SAMPLE picks a random point on the unit sphere.
!
!  Discussion:
!
!    The unit sphere in 3D satisfies:
!
!      X * X + Y * Y + Z * Z = 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of samples.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
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

  end  do

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
