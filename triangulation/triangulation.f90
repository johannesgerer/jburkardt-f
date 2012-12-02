subroutine alpha_measure ( n, z, element_order, element_num, element_node, &
  alpha_min, alpha_ave, alpha_area )

!*****************************************************************************80
!
!! ALPHA_MEASURE determines the triangulated pointset quality measure ALPHA.
!
!  Discusion:
!
!    The ALPHA measure evaluates the uniformity of the shapes of the triangles
!    defined by a triangulated pointset.
!
!    We compute the minimum angle among all the triangles in the triangulated
!    dataset and divide by the maximum possible value (which, in degrees,
!    is 60).  The best possible value is 1, and the worst 0.  A good
!    triangulation should have an ALPHA score close to 1.
!
!    The code has been modified to 'allow' 6-node triangulations.
!    However, no effort is made to actually process the midside nodes.
!    Only information from the vertices is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(2,N), the points.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the triangulation.
!
!    Output, real ( kind = 8 ) ALPHA_MIN, the minimum value of ALPHA over all
!    triangles.
!
!    Output, real ( kind = 8 ) ALPHA_AVE, the value of ALPHA averaged over
!    all triangles.
!
!    Output, real ( kind = 8 ) ALPHA_AREA, the value of ALPHA averaged over
!    all triangles and weighted by area.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  real ( kind = 8 ) a_angle
  integer ( kind = 4 ) a_index
  real ( kind = 8 ) a_x
  real ( kind = 8 ) a_y
  real ( kind = 8 ) ab_len
  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha_area
  real ( kind = 8 ) alpha_ave
  real ( kind = 8 ) alpha_min
  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) area
  real ( kind = 8 ) area_total
  real ( kind = 8 ) b_angle
  integer ( kind = 4 ) b_index
  real ( kind = 8 ) b_x
  real ( kind = 8 ) b_y
  real ( kind = 8 ) bc_len
  real ( kind = 8 ) c_angle
  integer ( kind = 4 ) c_index
  real ( kind = 8 ) c_x
  real ( kind = 8 ) c_y
  real ( kind = 8 ) ca_len
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) element_node(element_order,element_num)
  real ( kind = 8 ) z(2,n)

  alpha_min = huge ( alpha )
  alpha_ave = 0.0D+00
  alpha_area = 0.0D+00
  area_total = 0.0D+00

  do triangle = 1, element_num

    a_index = element_node(1,triangle)
    b_index = element_node(2,triangle)
    c_index = element_node(3,triangle)

    a_x = z(1,a_index)
    a_y = z(2,a_index)
    b_x = z(1,b_index)
    b_y = z(2,b_index)
    c_x = z(1,c_index)
    c_y = z(2,c_index)

    area = 0.5D+00 * abs ( a_x * ( b_y - c_y ) &
                         + b_x * ( c_y - a_y ) &
                         + c_x * ( a_y - b_y ) )

    ab_len = sqrt ( ( a_x - b_x )**2 + ( a_y - b_y )**2 )
    bc_len = sqrt ( ( b_x - c_x )**2 + ( b_y - c_y )**2 )
    ca_len = sqrt ( ( c_x - a_x )**2 + ( c_y - a_y )**2 )
!
!  Take care of a ridiculous special case.
!
    if ( ab_len == 0.0D+00 .and. &
         bc_len == 0.0D+00 .and. &
         ca_len == 0.0D+00 ) then

      a_angle = 2.0D+00 * pi / 3.0D+00
      b_angle = 2.0D+00 * pi / 3.0D+00
      c_angle = 2.0D+00 * pi / 3.0D+00

    else

      if ( ca_len == 0.0D+00 .or. ab_len == 0.0D+00 ) then
        a_angle = pi
      else
        a_angle = arc_cosine ( ( ca_len**2 + ab_len**2 - bc_len**2 ) &
          / ( 2.0D+00 * ca_len * ab_len ) )
      end if

      if ( ab_len == 0.0D+00 .or. bc_len == 0.0D+00 ) then
        b_angle = pi
      else
        b_angle = arc_cosine ( ( ab_len**2 + bc_len**2 - ca_len**2 ) &
          / ( 2.0D+00 * ab_len * bc_len ) )
      end if

      if ( bc_len == 0.0D+00 .or. ca_len == 0.0D+00 ) then
        c_angle = pi
      else
        c_angle = arc_cosine ( ( bc_len**2 + ca_len**2 - ab_len**2 ) &
          / ( 2.0D+00 * bc_len * ca_len ) )
      end if

    end if

    alpha_min = min ( alpha_min, a_angle )
    alpha_min = min ( alpha_min, b_angle )
    alpha_min = min ( alpha_min, c_angle )

    alpha_ave = alpha_ave + alpha_min

    alpha_area = alpha_area + area * alpha_min

    area_total = area_total + area

  end do

  alpha_ave = alpha_ave / real ( element_num, kind = 8 )
  alpha_area = alpha_area / area_total
!
!  Normalize angles from [0,pi/3] degrees into qualities in [0,1].
!
  alpha_min = alpha_min * 3.0D+00 / pi
  alpha_ave = alpha_ave * 3.0D+00 / pi
  alpha_area = alpha_area * 3.0D+00 / pi

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
subroutine area_measure ( n, z, element_order, element_num, element_node, &
  area_min, area_max, area_ratio, area_ave, area_std )

!*****************************************************************************80
!
!! AREA_MEASURE determines the area ratio quality measure.
!
!  Discusion:
!
!    This measure computes the area of every triangle, and returns
!    the ratio of the minimum to the maximum triangle.  A value of
!    1 is "perfect", indicating that all triangles have the same area.
!    A value of 0 is the worst possible result.
!
!    The code has been modified to 'allow' 6-node triangulations.
!    However, no effort is made to actually process the midside nodes.
!    Only information from the vertices is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(2,N), the points.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the triangulation.
!
!    Output, real ( kind = 8 ) AREA_MIN, AREA_MAX, the minimum and maximum
!    areas.
!
!    Output, real ( kind = 8 ) AREA_RATIO, the ratio of the minimum to the
!    maximum area.
!
!    Output, real ( kind = 8 ) AREA_AVE, the average area.
!
!    Output, real ( kind = 8 ) AREA_STD, the standard deviation of the areas.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  real ( kind = 8 ) area
  real ( kind = 8 ) area_ave
  real ( kind = 8 ) area_max
  real ( kind = 8 ) area_min
  real ( kind = 8 ) area_ratio
  real ( kind = 8 ) area_std
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) element_node(element_order,element_num)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) z(2,n)

  area_max = 0.0D+00
  area_min = huge ( area_min )
  area_ave = 0.0

  do triangle = 1, element_num

    x1 = z(1,element_node(1,triangle))
    y1 = z(2,element_node(1,triangle))
    x2 = z(1,element_node(2,triangle))
    y2 = z(2,element_node(2,triangle))
    x3 = z(1,element_node(3,triangle))
    y3 = z(2,element_node(3,triangle))

    area = 0.5D+00 * abs ( x1 * ( y2 - y3 ) &
                         + x2 * ( y3 - y1 ) &
                         + x3 * ( y1 - y2 ) )

    area_min = min ( area_min, area )
    area_max = max ( area_max, area )

    area_ave = area_ave + area

  end do

  area_ave = area_ave / real ( element_num, kind = 8 )

  area_std = 0.0D+00
  do triangle = 1, element_num

    x1 = z(1,element_node(1,triangle))
    y1 = z(2,element_node(1,triangle))
    x2 = z(1,element_node(2,triangle))
    y2 = z(2,element_node(2,triangle))
    x3 = z(1,element_node(3,triangle))
    y3 = z(2,element_node(3,triangle))

    area = 0.5D+00 * abs ( x1 * ( y2 - y3 ) &
                         + x2 * ( y3 - y1 ) &
                         + x3 * ( y1 - y2 ) )

    area_std = area_std + ( area - area_ave )**2
  end do
  area_std = sqrt ( area_std / real ( element_num, kind = 8 ) )

  if ( 0.0D+00 < area_max ) then
    area_ratio = area_min / area_max
  else
    area_ratio = 0.0D+00
  end if

  return
end
subroutine bandwidth ( element_order, element_num, element_node, ml, mu, m )

!*****************************************************************************80
!
!! BANDWIDTH determines the bandwidth associated with a finite element mesh.
!
!  Discussion:
!
!    The quantity computed here is the "geometric" bandwidth determined
!    by the finite element mesh alone.
!
!    If a single finite element variable is associated with each node
!    of the mesh, and if the nodes and variables are numbered in the
!    same way, then the geometric bandwidth is the same as the bandwidth
!    of a typical finite element matrix.
!
!    The bandwidth M is defined in terms of the lower and upper bandwidths:
!
!      M = ML + 1 + MU
!
!    where
!
!      ML = maximum distance from any diagonal entry to a nonzero
!      entry in the same row, but earlier column,
!
!      MU = maximum distance from any diagonal entry to a nonzero
!      entry in the same row, but later column.
!
!    Because the finite element node adjacency relationship is symmetric,
!    we are guaranteed that ML = MU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Output, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths
!    of the matrix.
!
!    Output, integer ( kind = 4 ) M, the bandwidth of the matrix.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) global_i
  integer ( kind = 4 ) global_j
  integer ( kind = 4 ) local_i
  integer ( kind = 4 ) local_j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu

  ml = 0
  mu = 0

  do element = 1, element_num

    do local_i = 1, element_order
      global_i = element_node(local_i,element)

      do local_j = 1, element_order
        global_j = element_node(local_j,element)

        mu = max ( mu, global_j - global_i )
        ml = max ( ml, global_i - global_j )

      end do
    end do
  end do

  m = ml + 1 + mu

  return
end
subroutine delaunay_swap_test ( xy, swap )

!*****************************************************************************80
!
!! DELAUNAY_SWAP_TEST performs the Delaunay swap test.
!
!  Discussion:
!
!    The current triangles are formed by nodes (1,2,3) and (1,3,4).
!    if a swap is recommended, the new triangles will be (1,2,4) and (2,3,4).
!
!      4     ?     4
!     / \         /|\
!    1---3  ==>  1 | 3
!     \ /         \|/
!      2           2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Graham Carey,
!    Computational Grids:
!    Generation, Adaptation and Solution Strategies,
!    Taylor and Francis, 1997,
!    ISBN13: 978-1560326359,
!    LC: QA377.C32.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XY(2,4), the coordinates of four points.
!
!    Output, logical SWAP, is TRUE if the triangles (1,2,4) and (2,3,4)
!    are to replace triangles (1,2,3) and (1,3,4).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  logical              swap
  real ( kind = 8 ) x13
  real ( kind = 8 ) x14
  real ( kind = 8 ) x23
  real ( kind = 8 ) x24
  real ( kind = 8 ) xy(2,4)
  real ( kind = 8 ) y13
  real ( kind = 8 ) y14
  real ( kind = 8 ) y23
  real ( kind = 8 ) y24

  x13 = xy(1,1) - xy(1,3)
  x14 = xy(1,1) - xy(1,4)
  x23 = xy(1,2) - xy(1,3)
  x24 = xy(1,2) - xy(1,4)

  y13 = xy(2,1) - xy(2,3)
  y14 = xy(2,1) - xy(2,4)
  y23 = xy(2,2) - xy(2,3)
  y24 = xy(2,2) - xy(2,4)

  a = x13 * x23 + y13 * y23
  b = x24 * y14 - x14 * y24
  c = x23 * y13 - x13 * y23
  d = x24 * x14 + y14 * y24
!
!  The reference gives two initial tests before the
!  main one.  However, there seems to be an error
!  in at least one of these tests.  Since they are
!  intended to avoid error in borderline cases, but
!  instead cause real error in common cases, they are
!  omitted for now.
!
! if ( 0.0D+00 <= a .and. 0.0D+00 <= d ) then
!   swap = .true.
! else if ( a < d .and. d < 0.0D+00 ) then
!   swap = .true.
!  else if...

  if ( a * b < c * d ) then
    swap = .true.
  else
    swap = .false.
  end if

  return
end
function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! DIAEDG chooses a diagonal edge.
!
!  Discussion:
!
!    The routine determines whether 0--2 or 1--3 is the diagonal edge
!    that should be chosen, based on the circumcircle criterion, where
!    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
!    quadrilateral in counterclockwise order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!    coordinates of the vertices of a quadrilateral, given in
!    counterclockwise order.
!
!    Output, integer ( kind = 4 ) DIAEDG, chooses a diagonal:
!    +1, if diagonal edge 02 is chosen;
!    -1, if diagonal edge 13 is chosen;
!     0, if the four vertices are cocircular.
!
  implicit none

  real ( kind = 8 ) ca
  real ( kind = 8 ) cb
  integer ( kind = 4 ) diaedg
  real ( kind = 8 ) dx10
  real ( kind = 8 ) dx12
  real ( kind = 8 ) dx30
  real ( kind = 8 ) dx32
  real ( kind = 8 ) dy10
  real ( kind = 8 ) dy12
  real ( kind = 8 ) dy30
  real ( kind = 8 ) dy32
  real ( kind = 8 ) s
  real ( kind = 8 ) tol
  real ( kind = 8 ) tola
  real ( kind = 8 ) tolb
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  tol = 100.0D+00 * epsilon ( tol )

  dx10 = x1 - x0
  dy10 = y1 - y0
  dx12 = x1 - x2
  dy12 = y1 - y2
  dx30 = x3 - x0
  dy30 = y3 - y0
  dx32 = x3 - x2
  dy32 = y3 - y2

  tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
  tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

  ca = dx10 * dx30 + dy10 * dy30
  cb = dx12 * dx32 + dy12 * dy32

  if ( tola < ca .and. tolb < cb ) then

    diaedg = - 1

  else if ( ca < - tola .and. cb < - tolb ) then

    diaedg = 1

  else

    tola = max ( tola, tolb )

    s = ( dx10 * dy30 - dx30 * dy10 ) * cb &
      + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( tola < s ) then
      diaedg = - 1
    else if ( s < - tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

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
!    17 November 2004
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

  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 ) today
  integer ( kind = 4 ) values(8)
  character ( len = 5 ) zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) / 11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) / 30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) / 23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) / 59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) / 59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp / 6.0D+00
!
!  Force 0 < TEMP <= 1.
!
  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( huge ( 1 ), kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == huge ( 1 ) ) then
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
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
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
function i4_sign ( x )

!*****************************************************************************80
!
!! I4_SIGN evaluates the sign of an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
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
!    Input, integer ( kind = 4 ) X, the number whose sign is desired.
!
!    Output, integer ( kind = 4 ) I4_SIGN, the sign of X:
!
  implicit none

  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    i4_sign = -1
  else
    i4_sign = +1
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP switches two I4's.
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
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
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
!    Input, integer ( kind = 4 ) IVAL, a value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) value
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    stop
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
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

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4col_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
!
!  Discussion:
!
!    The columns of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      unique_num = unique_num + 1
      j1 = j2
    end if

  end do

  return
end
subroutine i4col_swap ( m, n, a, i, j )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input/output, integer ( kind = 4 ) A(M,N), an array of N columns
!    of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I or J is out of bounds.'
    write ( *, '(a,i8)' ) '  I =    ', i
    write ( *, '(a,i8)' ) '  J =    ', j
    write ( *, '(a,i8)' ) '  N =    ', n
    stop

  end if

  if ( i == j ) then
    return
  end if

  col(1:m) = a(1:m,i)
  a(1:m,i) = a(1:m,j)
  a(1:m,j) = col(1:m)

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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 7 ) ctemp(incx)
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
      write ( ctemp(i2), '(i7)') i
    end do

    write ( *, '(''  Row '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a7)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an I4VEC into a descending heap.
!
!  Discussion:
!
!    A descending heap is an array A with the property that, for every index J,
!    A(2*J) <= A(J) and A(2*J+1) <= A(J), (as long as the indices
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
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
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
subroutine i4vec_min ( n, a, amin )

!*****************************************************************************80
!
!! I4VEC_MIN computes the minimum element of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) AMIN, the value of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) amin

  amin = minval ( a(1:n) )

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
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
  integer ( kind = 4 ) big
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(i8,1x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(i8,1x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i8,i11)' ) i, a(i)
    end do
  end if

  return
end
subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
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
!  1: Put A into descending heap form.
!
  call i4vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
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
    call i4vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    call i4_swap ( a(1), a(n1) )

  end do

  return
end
subroutine i4vec_sorted_unique ( n, a, unique_num )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE gets the unique elements in a sorted I4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in A.
!
!    Input/output, integer ( kind = 4 ) A(N).  On input, the sorted
!    integer array.  On output, the unique elements in A.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements
!    in A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) unique_num

  unique_num = 0

  if ( n <= 0 ) then
    return
  end if

  unique_num = 1

  do itest = 2, n

    if ( a(itest) /= a(unique_num) ) then
      unique_num = unique_num + 1
      a(unique_num) = a(itest)
    end if

  end do

  return
end
subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMPARE compares pairs of integers stored in two vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), contain the two components
!    of each item.
!
!    Input, integer ( kind = 4 ) I, J, the items to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item J < item I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
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
    else if ( a2(j) < a2(i) ) then
      isgn = +1
    end if

  else if ( a1(j) < a1(i) ) then

    isgn = +1

  end if

  return
end
subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) temp

  if ( n <= 1 ) then
    return
  end if
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

      temp  = a1(i)
      a1(i) = a1(j)
      a1(j) = temp

      temp  = a2(i)
      a2(i) = a2(j)
      a2(j) = temp
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4vec2_sorted_unique ( n, a1, a2, unique_num )

!*****************************************************************************80
!
!! I4VEC2_SORTED_UNIQUE gets the unique elements in a sorted I4VEC2.
!
!  Discussion:
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of unique items.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique items.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) unique_num

  unique_num = 0

  if ( n <= 0 ) then
    return
  end if

  unique_num = 1

  do itest = 2, n

    if ( a1(itest) /= a1(unique_num) .or. a2(itest) /= a2(unique_num) ) then

      unique_num = unique_num + 1

      a1(unique_num) = a1(itest)
      a2(unique_num) = a2(itest)

    end if

  end do

  return
end
function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

!*****************************************************************************80
!
!! LRLINE determines if a point is left of, right or, or on a directed line.
!
!  Discussion:
!
!    The directed line is parallel to, and at a signed distance DV from
!    a directed base line from (XV1,YV1) to (XV2,YV2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XU, YU, the coordinates of the point whose
!    position relative to the directed line is to be determined.
!
!    Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, the coordinates of two points
!    that determine the directed base line.
!
!    Input, real ( kind = 8 ) DV, the signed distance of the directed line
!    from the directed base line through the points (XV1,YV1) and (XV2,YV2).
!    DV is positive for a line to the left of the base line.
!
!    Output, integer ( kind = 4 ) LRLINE, the result:
!    +1, the point is to the right of the directed line;
!     0, the point is on the directed line;
!    -1, the point is to the left of the directed line.
!
  implicit none

  real ( kind = 8 ) dv
  real ( kind = 8 ) dx
  real ( kind = 8 ) dxu
  real ( kind = 8 ) dy
  real ( kind = 8 ) dyu
  integer ( kind = 4 ) lrline
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolabs
  real ( kind = 8 ) xu
  real ( kind = 8 ) xv1
  real ( kind = 8 ) xv2
  real ( kind = 8 ) yu
  real ( kind = 8 ) yv1
  real ( kind = 8 ) yv2

  tol = 100.0D+00 * epsilon ( tol )

  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), &
    abs ( dyu ), abs ( dv ) )

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

  if ( tolabs < t ) then
    lrline = 1
  else if ( -tolabs <= t ) then
    lrline = 0
  else
    lrline = -1
  end if

  return
end
subroutine lvec_print ( n, a, title )

!*****************************************************************************80
!
!! LVEC_PRINT prints an LVEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, logical A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  logical a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,l1)' ) i, a(i)
  end do

  return
end
subroutine mesh_base_one ( node_num, element_order, element_num, element_node )

!*****************************************************************************80
!
!! MESH_BASE_ONE ensures that the element definition is one-based.
!
!  Discussion:
!
!    The ELEMENT_NODE array contains nodes indices that form elements.
!    The convention for node indexing might start at 0 or at 1.
!    Since a FORTRAN90 program will naturally assume a 1-based indexing, it is
!    necessary to check a given element definition and, if it is actually
!    0-based, to convert it.
!
!    This function attempts to detect 9-based node indexing and correct it.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, int NODE_NUM, the number of nodes.
!
!    Input, int ELEMENT_ORDER, the order of the elements.
!
!    Input, int ELEMENT_NUM, the number of elements.
!
!    Input/output, int ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the element
!    definitions.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_max
  integer ( kind = 4 ) node_min
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) order

  node_min = node_num + 1
  node_max = -1

  node_min = minval ( element_node(1:element_order,1:element_num) )
  node_max = maxval ( element_node(1:element_order,1:element_num) )

  if ( node_min == 0 .and. node_max == node_num - 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'MESH_BASE_ONE:'
    write ( *, '(a)' )'  The element indexing appears to be 0-based!'
    write ( *, '(a)' )'  This will be converted to 1-based.'
    element_node(1:element_order,1:element_num) = &
      element_node(1:element_order,1:element_num) + 1
  else if ( node_min == 1 .and. node_max == node_num  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'MESH_BASE_ONE:'
    write ( *, '(a)' )'  The element indexing appears to be 1-based!'
    write ( *, '(a)' )'  No conversion is necessary.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_BASE_ONE - Warning!'
    write ( *, '(a)' ) '  The element indexing is not of a recognized type.'
    write ( *, '(a,i8)' ) '  NODE_MIN = ', node_min
    write ( *, '(a,i8)' ) '  NODE_MAX = ', node_max
    write ( *, '(a,i8)' ) '  NODE_NUM = ', node_num
  end if

  return
end
subroutine node_merge ( dim_num, node_num, node_xy, tolerance, node_rep )

!*****************************************************************************80
!
!! NODE_MERGE detects nodes that should be merged.
!
!  Discussion:
!
!    Two nodes "should" be merged if they are within TOLERANCE distance
!    of each other.
!
!    With a tolerance of 0, only exactly equal nodes are counted.
!
!    With a positive tolerance, a pair of nodes inside a circle of
!    radius TOLERANCE result in a count of 1 duplicate.
!
!    However, what do we do if nodes A, B and C are arranged in a line,!
!    with A and B just within TOLERANCE of each other, and B and C just
!    within tolerance of each other?  What we do here is make a choice
!    that can be defended consistently.  A and B define an equivalence
!    class because they are closer than TOLERANCE.  C is then added to
!    this equivalence class, because it is within TOLERANCE of at least
!    on thing in that equivalence class.
!
!    Thus, if 100 nodes are separated pairwise by slightly less
!    than TOLERANCE, a total of 99 duplicates will be counted.
!
!    The program starts out by giving each node its own label.
!    If it finds that two nodes should be merged, then the index of
!    one node is used as the label for both.  This process continues
!    until all nodes have been considered.  The number of unique nodes
!    is the number of unique values in the output quantity NODE_REP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(DIM_NUM,NODE_NUM), the nodes.
!
!    Input, real ( kind = 8 ) TOLERANCE, the maximum distance between
!    two nodes regarded as duplicate.
!
!    Output, integer ( kind = 4 ) NODE_REP(NODE_NUM), the "representative" of
!    each node.  NODE_REP(NODE) is the index of a node which is within
!    TOLERANCE of node NODE, or for which a chain of nodes can be found, all
!    having the same representative, and all of which are pairwise closer
!    than TOLERANCE.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) dist
  integer ( kind = 4 ) node_rep(node_num)
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  integer ( kind = 4 ) rep
  real ( kind = 8 ) rep_dist(node_num)
  real ( kind = 8 ) tolerance

  do node1 = 1, node_num
    node_rep(node1) = node1
  end do

  do node1 = 1, node_num

    rep_dist(1:node_num) = huge ( 1.0D+00 )

    do node2 = 1, node_num

      dist = sqrt ( sum ( &
        ( node_xy(1:dim_num,node1) - node_xy(1:dim_num,node2) )**2 ) )

      rep = node_rep(node2)

      if ( dist < rep_dist(rep) ) then
        rep_dist(rep) = dist
      end if

    end do

    do node2 = 1, node_num
      rep = node_rep(node2)
      if ( rep_dist(rep) <= tolerance ) then
        node_rep(node2) = node1
      end if
    end do

  end do

  return
end
subroutine ns_adj_col_set ( node_num, element_num, variable_num, &
  element_node, element_neighbor, node_u_variable, node_v_variable, &
  node_p_variable, adj_num, adj_col )

!*****************************************************************************80
!
!! NS_ADJ_COL_SET sets the COL array in a Navier Stokes triangulation.
!
!  Discussion:
!
!    This routine also returns the value of ADJ_NUM, the number of
!    Navier-Stokes variable adjacencies, which should be identical to the
!    value that would have been computed by calling NS_ADJ_COUNT.
!
!    This routine is called to set up the ADJ_COL array, which indicates
!    the number of entries needed to store each column in the sparse
!    compressed column storage to be used for the adjacency matrix.
!
!    The triangulation is assumed to involve 6-node triangles.
!
!    Variables for the horizontal and vertical velocities are associated
!    with every node.  Variables for the pressure are associated only with
!    the vertex nodes.
!
!    We are interested in determining the number of nonzero entries in the
!    stiffness matrix of the Stokes equations, or the jacobian matrix of
!    the Navier Stokes equations.  To this end, we will say, somewhat
!    too broadly, that two variables are "adjacent" if their associated
!    nodes both occur in some common element.  This adjacency of variables
!    I and J is taken to be equivalent to the possible nonzeroness of
!    matrix entries A(I,J) and A(J,I).
!
!    A sparse compressed column format is used to store the counts for
!    the nonzeroes.  In other words, while the value ADJ_NUM reports the
!    number of adjacencies, the vector ADJ_COL is sufficient to allow us
!    to properly set up a sparse compressed matrix for the actual storage
!    of the sparse matrix, if we desire to proceed.
!
!  Local Node Numbering:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  6   5  side 2
!       |    \
!    3  |     \
!       |      \
!       1---4---2
!
!         side 1
!
!  Variable Diagram:
!
!      UVP
!       |\
!       | \
!       |  \
!      UV   UV
!       |    \
!       |     \
!       |      \
!      UVP--UV--UVP
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) VARIABLE_NUM, the number of variables.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), lists the
!    nodes that make up each triangle.  The first three nodes are the vertices,
!    in counterclockwise order.  The fourth value is the midside
!    node between nodes 1 and 2; the fifth and sixth values are
!    the other midside nodes in the logical order.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), for each
!    side of a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Input, integer ( kind = 4 ) NODE_U_VARIABLE(NODE_NUM),
!    NODE_V_VARIABLE(NODE_NUM), NODE_P_VARIABLE(NODE_NUM), the index of the
!    horizontal velocity, vertical velocity and pressure variables associated
!    with a node, or -1 if no such variable is associated with the node.
!
!    Output, integer ( kind = 4 ) ADJ_NUM, the number of
!    Navier Stokes variable adjacencies.
!
!    Output, integer ( kind = 4 ) ADJ_COL(VARIABLE_NUM+1).  Information about
!    variable J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ) variable_num

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_col(variable_num+1)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_p_variable(node_num)
  integer ( kind = 4 ) node_u_variable(node_num)
  integer ( kind = 4 ) node_v_variable(node_num)
  integer ( kind = 4 ) p1
  integer ( kind = 4 ) p2
  integer ( kind = 4 ) p3
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle2
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) u1
  integer ( kind = 4 ) u2
  integer ( kind = 4 ) u3
  integer ( kind = 4 ) u4
  integer ( kind = 4 ) u5
  integer ( kind = 4 ) u6
  integer ( kind = 4 ) v1
  integer ( kind = 4 ) v2
  integer ( kind = 4 ) v3
  integer ( kind = 4 ) v4
  integer ( kind = 4 ) v5
  integer ( kind = 4 ) v6
  integer ( kind = 4 ) variable

  adj_num = 0
!
!  Set every variable to be adjacent to itself.
!
  adj_col(1:variable_num) = 1
!
!  Set every variable to be adjacent to the other variables associated with
!  that node.
!
!  U <=> V
!  U <=> P (if there is a P variable)
!  V <=> P (if there is a P variable)
!
  do node = 1, node_num

    u1 = node_u_variable(node)
    v1 = node_v_variable(node)
    p1 = node_p_variable(node)

    adj_col(u1) = adj_col(u1) + 1
    adj_col(v1) = adj_col(v1) + 1

    if ( 0 < p1 ) then
      adj_col(u1) = adj_col(u1) + 1
      adj_col(v1) = adj_col(v1) + 1
      adj_col(p1) = adj_col(p1) + 2
    end if

  end do
!
!  Examine each triangle.
!
  do triangle = 1, element_num

    n1 = element_node(1,triangle)
    n2 = element_node(2,triangle)
    n3 = element_node(3,triangle)
    n4 = element_node(4,triangle)
    n5 = element_node(5,triangle)
    n6 = element_node(6,triangle)

    u1 = node_u_variable(n1)
    v1 = node_v_variable(n1)
    p1 = node_p_variable(n1)

    u2 = node_u_variable(n2)
    v2 = node_v_variable(n2)
    p2 = node_p_variable(n2)

    u3 = node_u_variable(n3)
    v3 = node_v_variable(n3)
    p3 = node_p_variable(n3)

    u4 = node_u_variable(n4)
    v4 = node_v_variable(n4)

    u5 = node_u_variable(n5)
    v5 = node_v_variable(n5)

    u6 = node_u_variable(n6)
    v6 = node_v_variable(n6)
!
!  For sure, we add the new adjacencies:
!
!    U5 V5 <=> U1 V1 P1
!    U6 V6 <=> U2 V2 P2
!    U4 V4 <=> U3 V3 P3
!    U5 V5 <=> U4 V4
!    U6 V6 <=> U4 V4
!    U6 V6 <=> U5 V5
!
    adj_col(u1) = adj_col(u1) + 2
    adj_col(v1) = adj_col(v1) + 2
    adj_col(p1) = adj_col(p1) + 2

    adj_col(u2) = adj_col(u2) + 2
    adj_col(v2) = adj_col(v2) + 2
    adj_col(p2) = adj_col(p2) + 2

    adj_col(u3) = adj_col(u3) + 2
    adj_col(v3) = adj_col(v3) + 2
    adj_col(p3) = adj_col(p3) + 2

    adj_col(u4) = adj_col(u4) + 7
    adj_col(v4) = adj_col(v4) + 7

    adj_col(u5) = adj_col(u5) + 7
    adj_col(v5) = adj_col(v5) + 7

    adj_col(u6) = adj_col(u6) + 7
    adj_col(v6) = adj_col(v6) + 7
!
!  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
!  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
!  Maybe add
!
!    U1 V1 P1 <=> U2 V2 P2
!    U1 V1 P1 <=> U4 V4
!    U2 V2 P2 <=> U4 V4
!
    triangle2 = element_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then

      adj_col(u1) = adj_col(u1) + 5
      adj_col(v1) = adj_col(v1) + 5
      adj_col(p1) = adj_col(p1) + 5

      adj_col(u2) = adj_col(u2) + 5
      adj_col(v2) = adj_col(v2) + 5
      adj_col(p2) = adj_col(p2) + 5

      adj_col(u4) = adj_col(u4) + 6
      adj_col(v4) = adj_col(v4) + 6

    end if
!
!  Maybe add
!
!    U2 V2 P2 <=> U3 V3 P3
!    U2 V2 P2 <=> U5 V5
!    U3 V3 P3 <=> U5 V5
!
    triangle2 = element_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then

      adj_col(u2) = adj_col(u2) + 5
      adj_col(v2) = adj_col(v2) + 5
      adj_col(p2) = adj_col(p2) + 5

      adj_col(u3) = adj_col(u3) + 5
      adj_col(v3) = adj_col(v3) + 5
      adj_col(p3) = adj_col(p3) + 5

      adj_col(u5) = adj_col(u5) + 6
      adj_col(v5) = adj_col(v5) + 6

    end if
!
!  Maybe add
!
!    U1 V1 P1 <=> U3 V3 P3
!    U1 V1 P1 <=> U6 V6
!    U3 V3 P3 <=> U6 V6
!
    triangle2 = element_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then

      adj_col(u1) = adj_col(u1) + 5
      adj_col(v1) = adj_col(v1) + 5
      adj_col(p1) = adj_col(p1) + 5

      adj_col(u3) = adj_col(u3) + 5
      adj_col(v3) = adj_col(v3) + 5
      adj_col(p3) = adj_col(p3) + 5

      adj_col(u6) = adj_col(u6) + 6
      adj_col(v6) = adj_col(v6) + 6

    end if

  end do
!
!  We used ADJ_COL to count the number of entries in each column.
!  Convert it to pointers into the ADJ array.
!
  adj_col(2:variable_num+1) = adj_col(1:variable_num)

  adj_col(1) = 1
  do variable = 2, variable_num + 1
    adj_col(variable) = adj_col(variable-1) + adj_col(variable)
  end do

  adj_num = adj_col(variable_num+1) - 1

  return
end
subroutine ns_adj_count ( node_num, element_num, variable_num, element_node, &
  element_neighbor, node_u_variable, node_v_variable, node_p_variable, &
  adj_num )

!*****************************************************************************80
!
!! NS_ADJ_COUNT counts adjacencies in a Navier Stokes triangulation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The value of ADJ_NUM computed and returned by this routine should
!    be identical to the value computed by NS_ADJ_COL_SET.
!
!    The triangulation is assumed to involve 6-node triangles.
!
!    Variables for the horizontal and vertical velocities are associated
!    with every node.  Variables for the pressure are associated only with
!    the vertex nodes.
!
!    We are interested in determining the number of nonzero entries in the
!    stiffness matrix of the Stokes equations, or the jacobian matrix of
!    the Navier Stokes equations.  To this end, we will say, somewhat
!    too broadly, that two variables are "adjacent" if their associated
!    nodes both occur in some common element.  This adjacency of variables
!    I and J is taken to be equivalent to the possible nonzeroness of
!    matrix entries A(I,J) and A(J,I).
!
!    A sparse compressed column format is used to store the counts for
!    the nonzeroes.  In other words, while the value ADJ_NUM reports the
!    number of adjacencies, the vector ADJ_COL is sufficient to allow us
!    to properly set up a sparse compressed matrix for the actual storage
!    of the sparse matrix, if we desire to proceed.
!
!  Local Node Numbering:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  6   5  side 2
!       |    \
!    3  |     \
!       |      \
!       1---4---2
!
!         side 1
!
!  Variable Diagram:
!
!      UVP
!       |\
!       | \
!       |  \
!      UV   UV
!       |    \
!       |     \
!       |      \
!      UVP--UV--UVP
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) VARIABLE_NUM, the number of variables.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), lists the
!    nodes that make up each triangle.  The first three nodes are the vertices,
!    in counterclockwise order.  The fourth value is the midside
!    node between nodes 1 and 2; the fifth and sixth values are
!    the other midside nodes in the logical order.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), for each
!    side of a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Input, integer ( kind = 4 ) NODE_U_VARIABLE(NODE_NUM),
!    NODE_V_VARIABLE(NODE_NUM), NODE_P_VARIABLE(NODE_NUM), the index of the
!    horizontal velocity, vertical velocity and pressure variables associated
!    with a node, or -1 if no such variable is associated with the node.
!
!    Output, integer ( kind = 4 ) ADJ_NUM, the value of ADJ_NUM, the number of
!    Navier Stokes variable adjacencies.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ) variable_num

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_p_variable(node_num)
  integer ( kind = 4 ) node_u_variable(node_num)
  integer ( kind = 4 ) node_v_variable(node_num)
  integer ( kind = 4 ) p1
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle2
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)

  adj_num = 0
!
!  Set every variable to be adjacent to itself.
!
  adj_num = variable_num
!
!  Set every variable to be adjacent to the other variables associated with
!  that node.
!
!  U <=> V
!  U <=> P (if there is a P variable)
!  V <=> P (if there is a P variable)
!
  do node = 1, node_num

    adj_num = adj_num + 2

    p1 = node_p_variable(node)

    if ( 0 < p1 ) then
      adj_num = adj_num + 4
    end if

  end do
!
!  Examine each triangle.
!
  do triangle = 1, element_num
!
!  For sure, we add the new adjacencies:
!
!    U5 V5 <=> U1 V1 P1
!    U6 V6 <=> U2 V2 P2
!    U4 V4 <=> U3 V3 P3
!    U5 V5 <=> U4 V4
!    U6 V6 <=> U4 V4
!    U6 V6 <=> U5 V5
!
    adj_num = adj_num + 60
!
!  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
!  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
!  Maybe add
!
!    U1 V1 P1 <=> U2 V2 P2
!    U1 V1 P1 <=> U4 V4
!    U2 V2 P2 <=> U4 V4
!
    triangle2 = element_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_num = adj_num + 42
    end if
!
!  Maybe add
!
!    U2 V2 P2 <=> U3 V3 P3
!    U2 V2 P2 <=> U5 V5
!    U3 V3 P3 <=> U5 V5
!
    triangle2 = element_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_num = adj_num + 42
    end if
!
!  Maybe add
!
!    U1 V1 P1 <=> U3 V3 P3
!    U1 V1 P1 <=> U6 V6
!    U3 V3 P3 <=> U6 V6
!
    triangle2 = element_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_num = adj_num + 42
    end if

  end do

  return
end
subroutine ns_adj_insert ( v1, v2, variable_num, adj_num, adj_col_free, &
  adj_row )

!*****************************************************************************80
!
!! NS_ADJ_INSERT inserts an adjacency into a compressed column adjacency matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) V1, V2, the indices of two items which are
!    adjacent.
!
!    Input, integer ( kind = 4 ) VARIABLE_NUM, the number of items.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of entries available
!    in ADJ_ROW.
!
!    Input/output, integer ( kind = 4 ) ADJ_COL_FREE(VARIABLE_NUM), the next
!    free location in which an entry for a given column can be stored.  On
!    output, two pointers have been updated.
!
!    Input/output, integer ( kind = 4 ) ADJ_ROW(ADJ_NUM), the row indices of
!    the Navier Stokes variable adjacency matrix.  On output, two new entries
!    have been added.
!
  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) variable_num

  integer ( kind = 4 ) adj_col_free(variable_num)
  integer ( kind = 4 ) adj_row(adj_num)
  integer ( kind = 4 ) v1
  integer ( kind = 4 ) v2

  adj_row(adj_col_free(v1)) = v2
  adj_col_free(v1) = adj_col_free(v1) + 1

  if ( v1 == v2 ) then
    return
  end if

  adj_row(adj_col_free(v2)) = v1
  adj_col_free(v2) = adj_col_free(v2) + 1

  return
end
subroutine ns_adj_row_set ( node_num, element_num, variable_num, &
  element_node, element_neighbor, node_u_variable, node_v_variable, &
  node_p_variable, adj_num, adj_col, adj_row )

!*****************************************************************************80
!
!! NS_ADJ_ROW_SET sets the Navier Stokes sparse compressed column row indices.
!
!  Discussion:
!
!    After NS_ADJ_COUNT has been called to count ADJ_NUM, the number of
!    variable adjacencies and to set up ADJ_COL, the compressed column pointer,
!    this routine can be called to assign values to ADJ_ROW, the row
!    indices of the sparse compressed column adjacency matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) VARIABLE_NUM, the number of variables.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), lists the
!    nodes that make up each triangle.  The first three nodes are the vertices,
!    in counterclockwise order.  The fourth value is the midside
!    node between nodes 1 and 2; the fifth and sixth values are
!    the other midside nodes in the logical order.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), for each
!    side of a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Input, integer ( kind = 4 ) NODE_U_VARIABLE(NODE_NUM),
!    NODE_V_VARIABLE(NODE_NUM), NODE_P_VARIABLE(NODE_NUM), the index of the
!    horizontal velocity, vertical velocity and pressure variables associated
!    with a node, or -1 if no such variable is associated with the node.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of Navier Stokes variable
!    adjacencies.
!
!    Input, integer ( kind = 4 ) ADJ_COL(VARIABLE_NUM+1).  Information about
!    variable J  is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
!
!    Output, integer ( kind = 4 ) ADJ_ROW(ADJ_NUM), the row indices of the
!    Navier Stokes variable adjacency matrix.
!
!  Local Parameters:
!
!    Local, integer ADJ_COL_FREE(VARIABLE_NUM), for each column,
!    the location in ADJ_ROW which can store the next row index.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ) variable_num

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_col(variable_num+1)
  integer ( kind = 4 ) adj_col_free(variable_num)
  integer ( kind = 4 ) adj_row(adj_num)
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_p_variable(node_num)
  integer ( kind = 4 ) node_u_variable(node_num)
  integer ( kind = 4 ) node_v_variable(node_num)
  integer ( kind = 4 ) number
  integer ( kind = 4 ) p1
  integer ( kind = 4 ) p2
  integer ( kind = 4 ) p3
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle2
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) u1
  integer ( kind = 4 ) u2
  integer ( kind = 4 ) u3
  integer ( kind = 4 ) u4
  integer ( kind = 4 ) u5
  integer ( kind = 4 ) u6
  integer ( kind = 4 ) v
  integer ( kind = 4 ) v1
  integer ( kind = 4 ) v2
  integer ( kind = 4 ) v3
  integer ( kind = 4 ) v4
  integer ( kind = 4 ) v5
  integer ( kind = 4 ) v6

  adj_row(1:adj_num) = -1

  adj_col_free(1:variable_num) = adj_col(1:variable_num)
!
!  Set every variable to be adjacent to itself.
!
  do v = 1, variable_num
    call ns_adj_insert ( v, v, variable_num, adj_num, adj_col_free, adj_row )
  end do
!
!  Set every variable to be adjacent to the other variables associated with
!  that node.
!
!  U <=> V
!  U <=> P (if there is a P variable)
!  V <=> P (if there is a P variable)
!
  do node = 1, node_num

    u1 = node_u_variable(node)
    v1 = node_v_variable(node)
    p1 = node_p_variable(node)

    call ns_adj_insert ( u1, v1, variable_num, adj_num, adj_col_free, adj_row )

    if ( 0 < p1 ) then

      call ns_adj_insert ( u1, p1, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v1, p1, variable_num, adj_num, adj_col_free, &
        adj_row )

    end if

  end do
!
!  Examine each triangle.
!
  do triangle = 1, element_num

    n1 = element_node(1,triangle)
    n2 = element_node(2,triangle)
    n3 = element_node(3,triangle)
    n4 = element_node(4,triangle)
    n5 = element_node(5,triangle)
    n6 = element_node(6,triangle)

    u1 = node_u_variable(n1)
    v1 = node_v_variable(n1)
    p1 = node_p_variable(n1)

    u2 = node_u_variable(n2)
    v2 = node_v_variable(n2)
    p2 = node_p_variable(n2)

    u3 = node_u_variable(n3)
    v3 = node_v_variable(n3)
    p3 = node_p_variable(n3)

    u4 = node_u_variable(n4)
    v4 = node_v_variable(n4)

    u5 = node_u_variable(n5)
    v5 = node_v_variable(n5)

    u6 = node_u_variable(n6)
    v6 = node_v_variable(n6)
!
!  For sure, we add the new adjacencies:
!
!    U5 V5 <=> U1 V1 P1
!    U6 V6 <=> U2 V2 P2
!    U4 V4 <=> U3 V3 P3
!    U5 V5 <=> U4 V4
!    U6 V6 <=> U4 V4
!    U6 V6 <=> U5 V5
!
    call ns_adj_insert ( u5, u1, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( u5, v1, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( u5, p1, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v5, u1, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v5, v1, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v5, p1, variable_num, adj_num, adj_col_free, adj_row )

    call ns_adj_insert ( u6, u2, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( u6, v2, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( u6, p2, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v6, u2, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v6, v2, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v6, p2, variable_num, adj_num, adj_col_free, adj_row )

    call ns_adj_insert ( u4, u3, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( u4, v3, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( u4, p3, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v4, u3, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v4, v3, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v4, p3, variable_num, adj_num, adj_col_free, adj_row )

    call ns_adj_insert ( u5, u4, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( u5, v4, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v5, u4, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v5, v4, variable_num, adj_num, adj_col_free, adj_row )

    call ns_adj_insert ( u6, u4, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( u6, v4, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v6, u4, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v6, v4, variable_num, adj_num, adj_col_free, adj_row )

    call ns_adj_insert ( u6, u5, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( u6, v5, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v6, u5, variable_num, adj_num, adj_col_free, adj_row )
    call ns_adj_insert ( v6, v5, variable_num, adj_num, adj_col_free, adj_row )
!
!  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
!  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
!  Maybe add
!
!    U1 V1 P1 <=> U2 V2 P2
!    U1 V1 P1 <=> U4 V4
!    U2 V2 P2 <=> U4 V4
!
    triangle2 = element_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then

      call ns_adj_insert ( u1, u2, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( u1, v2, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( u1, p2, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v1, u2, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v1, v2, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v1, p2, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p1, u2, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p1, v2, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p1, p2, variable_num, adj_num, adj_col_free, &
        adj_row )

      call ns_adj_insert ( u1, u4, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( u1, v4, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v1, u4, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v1, v4, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p1, u4, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p1, v4, variable_num, adj_num, adj_col_free, &
        adj_row )

      call ns_adj_insert ( u2, u4, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( u2, v4, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v2, u4, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v2, v4, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p2, u4, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p2, v4, variable_num, adj_num, adj_col_free, &
        adj_row )

    end if
!
!  Maybe add
!
!    U2 V2 P2 <=> U3 V3 P3
!    U2 V2 P2 <=> U5 V5
!    U3 V3 P3 <=> U5 V5
!
    triangle2 = element_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then

      call ns_adj_insert ( u2, u3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( u2, v3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( u2, p3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v2, u3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v2, v3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v2, p3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p2, u3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p2, v3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p2, p3, variable_num, adj_num, adj_col_free, &
        adj_row )

      call ns_adj_insert ( u2, u5, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( u2, v5, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v2, u5, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v2, v5, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p2, u5, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p2, v5, variable_num, adj_num, adj_col_free, &
        adj_row )

      call ns_adj_insert ( u3, u5, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( u3, v5, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v3, u5, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v3, v5, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p3, u5, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p3, v5, variable_num, adj_num, adj_col_free, &
        adj_row )

    end if
!
!  Maybe add
!
!    U1 V1 P1 <=> U3 V3 P3
!    U1 V1 P1 <=> U6 V6
!    U3 V3 P3 <=> U6 V6
!
    triangle2 = element_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then

      call ns_adj_insert ( u1, u3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( u1, v3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( u1, p3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v1, u3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v1, v3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v1, p3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p1, u3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p1, v3, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p1, p3, variable_num, adj_num, adj_col_free, &
        adj_row )

      call ns_adj_insert ( u1, u6, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( u1, v6, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v1, u6, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v1, v6, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p1, u6, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p1, v6, variable_num, adj_num, adj_col_free, &
        adj_row )

      call ns_adj_insert ( u3, u6, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( u3, v6, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v3, u6, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( v3, v6, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p3, u6, variable_num, adj_num, adj_col_free, &
        adj_row )
      call ns_adj_insert ( p3, v6, variable_num, adj_num, adj_col_free, &
        adj_row )

    end if

  end do
!
!  Ascending sort the entries for each variable.
!
  do v = 1, variable_num
    k1 = adj_col(v)
    k2 = adj_col(v+1)-1
    number = k2 + 1 - k1
    call i4vec_sort_heap_a ( number, adj_row(k1:k2) )
  end do

  return
end
subroutine perm_check2 ( n, p, base, ierror )

!*****************************************************************************80
!
!! PERM_CHECK2 checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from BASE to
!    to BASE+N-1 occurs among the N entries of the permutation.
!
!    Set the input quantity BASE to 0, if P is a 0-based permutation,
!    or to 1 if P is a 1-based permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) P(N), the array to check.
!
!    Input, integer ( kind = 4 ) BASE, the index base.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) find
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seek

  ierror = 0

  do seek = base, base + n - 1

    ierror = 1

    do find = 1, n
      if ( p(find) == seek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK2 - Fatal error!'
      write ( *, '(a)' ) '  The input array does not represent'
      write ( *, '(a)' ) '  a proper permutation.'
      stop
    end if

  end do

  return
end
subroutine perm_inverse ( n, p )

!*****************************************************************************80
!
!! PERM_INVERSE inverts a permutation "in place".
!
!  Discussion:
!
!    This algorithm assumes that the entries in the permutation vector are
!    strictly positive.  In particular, the value 0 must not occur.
!
!    When necessary, this function shifts the data temporarily so that
!    this requirement is satisfied.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard
!    index form.  On output, P describes the inverse permutation
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) is
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) p_min

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop
  end if
!
!  Find the least value, and shift data so it begins at 1.
!
  call i4vec_min ( n, p, p_min )
  base = 1
  p(1:n) = p(1:n) - p_min + base
!
!  Check the permutation.
!
  call perm_check2 ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if
!
!  Invert the permutation.
!
  is = 1

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    is = - i4_sign ( p(i) )
    p(i) = is * abs ( p(i) )

  end do

  do i = 1, n

    i1 = - p(i)

    if ( 0 <= i1 ) then

      i0 = i

      do

        i2 = p(i1)
        p(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do
!
!  Reverse the shift.
!
  p(1:n) = p(1:n) + p_min - base

  return
end
subroutine points_delaunay_naive_2d ( node_num, node_xy, maxtri, &
  element_num, element_node )

!*****************************************************************************80
!
!! POINTS_DELAUNAY_NAIVE_2D is a naive Delaunay triangulation scheme.
!
!  Discussion:
!
!    This routine is only suitable as a demonstration code for small
!    problems.  Its running time is of order NODE_NUM**4.  Much faster
!    algorithms are available.
!
!    Given a set of nodes in the plane, a triangulation is set of
!    triples of distinct nodes, forming triangles, so that every
!    point within the convex hull of the set of nodes is either
!    one of the nodes, or lies on an edge of one or more triangles,
!    or lies within exactly one triangle.
!
!    A Delaunay triangulation is a triangulation with additional
!    properties.
!
!    NODE_NUM must be at least 3.
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
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry,
!    Cambridge University Press,
!    Second Edition, 1998, page 187.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) MAXTRI, the maximum number of triangles.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles in
!    the triangulation.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,MAXTRI), the indices of
!    the triangle nodes.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) maxtri
  integer ( kind = 4 ) node_num

  logical              flag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) element_node(3,maxtri)
  integer ( kind = 4 ) element_num
  real ( kind = 8 ) xn
  real ( kind = 8 ) yn
  real ( kind = 8 ) z(node_num)
  real ( kind = 8 ) zn

  element_num = 0

  if ( node_num < 3 ) then
    return
  end if
!
!  Compute Z = X*X + Y*Y.
!
  z(1:node_num) = node_xy(1,1:node_num)**2 + node_xy(2,1:node_num)**2
!
!  For each triple (I,J,K):
!
  do i = 1, node_num - 2
    do j = i+1, node_num
      do k = i+1, node_num

        if ( j /= k ) then

          xn = ( node_xy(2,j) - node_xy(2,i) ) * ( z(k) - z(i) ) &
             - ( node_xy(2,k) - node_xy(2,i) ) * ( z(j) - z(i) )

          yn = ( node_xy(1,k) - node_xy(1,i) ) * ( z(j) - z(i) ) &
             - ( node_xy(1,j) - node_xy(1,i) ) * ( z(k) - z(i) )

          zn = ( node_xy(1,j) - node_xy(1,i) ) &
             * ( node_xy(2,k) - node_xy(2,i) ) &
             - ( node_xy(1,k) - node_xy(1,i) ) &
             * ( node_xy(2,j) - node_xy(2,i) )

          flag = ( zn < 0.0D+00 )

          if ( flag ) then
            do m = 1, node_num
              flag = flag .and. &
                ( ( node_xy(1,m) - node_xy(1,i) ) * xn &
                + ( node_xy(2,m) - node_xy(2,i) ) * yn &
                + ( z(m)   - z(i) )   * zn <= 0.0D+00 )
            end do
          end if

          if ( flag ) then
            if ( element_num < maxtri ) then
              element_num = element_num + 1
              element_node(1:3,element_num) = (/ i, j, k /)
            end if
          end if

        end if

      end do
    end do
  end do

  return
end
subroutine points_hull_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_HULL_2D computes the convex hull of 2D points.
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
!    Output, integer ( kind = 4 ) HULL_NUM, the number of nodes that lie on
!    the convex hull.
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
       ( node_xy(1,i) == node_xy(1,q) .and. node_xy(2,i) < node_xy(2,q) ) ) then
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
      write ( *, '(a)' ) 'POINTS_HULL_2D - Fatal error!'
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

  return
end
subroutine points_point_near_naive_nd ( dim_num, nset, pset, p, i_min, d_min )

!*****************************************************************************80
!
!! POINTS_POINT_NEAR_NAIVE_ND finds the nearest point to a given point in ND.
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(DIM_NUM,NSET), the points in the set.
!
!    Input, real ( kind = 8 ) P(DIM_NUM), the point whose nearest neighbor
!    is sought.
!
!    Output, integer ( kind = 4 ) I_MIN, the index of the nearest point in
!    PSET to P.
!
!    Output, real ( kind = 8 ) D_MIN, the distance between P(*)
!    and PSET(*,I_MIN).
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nset

  real ( kind = 8 ) d
  real ( kind = 8 ) d_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) pset(dim_num,nset)

  d_min = huge ( d_min )
  i_min = -1

  do i = 1, nset
    d = sum ( ( p(1:dim_num) - pset(1:dim_num,i) )**2 )
    if ( d < d_min ) then
      d_min = d
      i_min = i
    end if
  end do

  d_min = sqrt ( d_min )

  return
end
subroutine q_measure ( n, z, element_order, element_num, element_node, &
  q_min, q_max, q_ave, q_area )

!*****************************************************************************80
!
!! Q_MEASURE determines the triangulated pointset quality measure Q.
!
!  Discussion:
!
!    The Q measure evaluates the uniformity of the shapes of the triangles
!    defined by a triangulated pointset.
!
!    For a single triangle T, the value of Q(T) is defined as follows:
!
!      TAU_IN = radius of the inscribed circle,
!      TAU_OUT = radius of the circumscribed circle,
!
!      Q(T) = 2 * TAU_IN / TAU_OUT
!        = ( B + C - A ) * ( C + A - B ) * ( A + B - C ) / ( A * B * C )
!
!    where A, B and C are the lengths of the sides of the triangle T.
!
!    The Q measure computes the value of Q(T) for every triangle T in the
!    triangulation, and then computes the minimum of this
!    set of values:
!
!      Q_MEASURE = min ( all T in triangulation ) Q(T)
!
!    In an ideally regular mesh, all triangles would have the same
!    equilateral shape, for which Q = 1.  A good mesh would have
!    0.5 < Q.
!
!    Given the 2D coordinates of a set of N nodes, stored as Z(1:2,1:N),
!    a triangulation is a list of ELEMENT_NUM triples of node indices that form
!    triangles.  Generally, a maximal triangulation is expected, namely,
!    a triangulation whose image is a planar graph, but for which the
!    addition of any new triangle would mean the graph was no longer planar.
!    A Delaunay triangulation is a maximal triangulation which maximizes
!    the minimum angle that occurs in any triangle.
!
!    The code has been modified to 'allow' 6-node triangulations.
!    However, no effort is made to actually process the midside nodes.
!    Only information from the vertices is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Max Gunzburger and John Burkardt,
!    Uniformity Measures for Point Samples in Hypercubes.
!
!    Per-Olof Persson and Gilbert Strang,
!    A Simple Mesh Generator in MATLAB,
!    SIAM Review,
!    Volume 46, Number 2, pages 329-345, June 2004.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(2,N), the points.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the triangulation.
!
!    Output, real ( kind = 8 ) Q_MIN, Q_MAX, the minimum and maximum values
!    of Q over all triangles.
!
!    Output, real ( kind = 8 ) Q_AVE, the average value of Q.
!
!    Output, real ( kind = 8 ) Q_AREA, the average value of Q, weighted by
!    the area of each triangle.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) a_index
  real ( kind = 8 ) ab_length
  real ( kind = 8 ) area
  real ( kind = 8 ) area_total
  integer ( kind = 4 ) b_index
  real ( kind = 8 ) bc_length
  integer ( kind = 4 ) c_index
  real ( kind = 8 ) ca_length
  real ( kind = 8 ) q
  real ( kind = 8 ) q_area
  real ( kind = 8 ) q_ave
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) element_node(element_order,element_num)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) z(2,n)

  q_min =   huge ( q_min )
  q_max = - huge ( q_max )
  q_ave = 0.0D+00
  q_area = 0.0D+00
  area_total = 0.0D+00

  do triangle = 1, element_num

    a_index = element_node(1,triangle)
    b_index = element_node(2,triangle)
    c_index = element_node(3,triangle)

    ab_length = sqrt ( &
        ( z(1,a_index) - z(1,b_index) )**2 &
      + ( z(2,a_index) - z(2,b_index) )**2 )

    bc_length = sqrt ( &
        ( z(1,b_index) - z(1,c_index) )**2 &
      + ( z(2,b_index) - z(2,c_index) )**2 )

    ca_length = sqrt ( &
        ( z(1,c_index) - z(1,a_index) )**2 &
      + ( z(2,c_index) - z(2,a_index) )**2 )

    q = ( bc_length + ca_length - ab_length ) &
      * ( ca_length + ab_length - bc_length ) &
      * ( ab_length + bc_length - ca_length ) &
      / ( ab_length * bc_length * ca_length )

    x1 = z(1,element_node(1,triangle))
    y1 = z(2,element_node(1,triangle))
    x2 = z(1,element_node(2,triangle))
    y2 = z(2,element_node(2,triangle))
    x3 = z(1,element_node(3,triangle))
    y3 = z(2,element_node(3,triangle))

    area = 0.5D+00 * abs ( x1 * ( y2 - y3 ) &
                         + x2 * ( y3 - y1 ) &
                         + x3 * ( y1 - y2 ) )

    q_min = min ( q_min, q )
    q_max = max ( q_max, q )
    q_ave = q_ave + q
    q_area = q_area + q * area

    area_total = area_total + area

  end do

  q_ave = q_ave / real ( element_num, kind = 8 )

  if ( 0.0D+00 < area_total ) then
    q_area = q_area / area_total
  else
    q_area = 0.0D+00
  end if

  return
end
subroutine quad_convex_random ( seed, xy )

!*****************************************************************************80
!
!! QUAD_CONVEX_RANDOM returns a random convex quadrilateral.
!
!  Description:
!
!    The quadrilateral is constrained in that the vertices must all lie
!    with the unit square.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) XY(2,NODE_NUM), the coordinates of the
!    nodes of the quadrilateral, given in counterclockwise order.
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 4

  integer ( kind = 4 ) hull(node_num)
  integer ( kind = 4 ) hull_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) xy(2,node_num)
  real ( kind = 8 ) xy_random(2,node_num)

  do
!
!  Generate 4 random points.
!
    call r8mat_uniform_01 ( 2, node_num, seed, xy_random )
!
!  Determine the convex hull.
!
    call points_hull_2d ( node_num, xy_random, hull_num, hull )
!
!  If HULL_NUM < NODE_NUM, then our convex hull is a triangle.
!  Try again.
!
    if ( hull_num == node_num ) then
      exit
    end if

  end do
!
!  Make an ordered copy of the random points.
!
  do j = 1, node_num
    xy(1:2,j) = xy_random(1:2,hull(j))
  end do

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
subroutine r82vec_permute ( n, p, a )

!*****************************************************************************80
!
!! R82VEC_PERMUTE permutes an R82VEC in place.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of R8 values.
!
!    The same logic can be used to permute an array of objects of any
!    arithmetic type, or an array of objects of any complexity.  The only
!    temporary storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.
!
!    Input/output, real ( kind = 8 ) A(2,N), the array to be permuted.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num,n)
  real ( kind = 8 ) a_temp(dim_num)
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check2 ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp(1:dim_num) = a(1:dim_num,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(1:dim_num,iput) = a_temp(1:dim_num)
          exit
        end if

        a(1:dim_num,iput) = a(1:dim_num,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)

  return
end
subroutine r82vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R82VEC_SORT_HEAP_INDEX_A ascending index heaps an R82VEC.
!
!  Discussion:
!
!    An R82VEC is an array of R82's.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(1:2,INDX(1:N)) is sorted,
!
!    or explicitly, by the call
!
!      call r82vec_permute ( n, indx, a )
!
!    after which A(1:2,I), I = 1 to N is sorted.
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
!    Input, real ( kind = 8 ) A(2,N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(1:2,INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num,n)
  real ( kind = 8 ) aval(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval(1:dim_num) = a(1:dim_num,indxt)

    else

      indxt = indx(ir)
      aval(1:dim_num) = a(1:dim_num,indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if (   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
             ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
               a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
          j = j + 1
        end if
      end if

      if (   aval(1) <  a(1,indx(j)) .or. &
           ( aval(1) == a(1,indx(j)) .and. &
             aval(2) <  a(2,indx(j)) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

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
!    26 March 2005
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
  character ( len = * )  title

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

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

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
      write ( ctemp(i2), '(i7,7x)') i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'

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
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8's.
!
!    For now, the input quantity SEED is an integer variable.
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
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
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
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8tris2 ( node_num, node_xy, element_num, element_node, &
  element_neighbor )

!*****************************************************************************80
!
!! R8TRIS2 constructs a Delaunay triangulation of 2D vertices.
!
!  Discussion:
!
!    The routine constructs the Delaunay triangulation of a set of 2D vertices
!    using an incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (X,Y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input/output, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates
!    of the nodes.  On output, the vertices have been sorted into
!    dictionary order.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles in the
!    triangulation;  ELEMENT_NUM is equal to 2*NODE_NUM - NB - 2, where NB is
!    the number of boundary vertices.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes that
!    make up each triangle.  The elements are indices of P.  The vertices of
!    the triangles are in counter clockwise order.
!
!    Output, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbor list.  Positive elements are indices of TIL; negative
!    elements are used for links of a counter clockwise linked list of boundary
!    edges;  LINK = -(3*I + J-1) where I, J = triangle, edge index;
!    ELEMENT_NEIGHBOR(J,I) refers to the neighbor along edge from vertex J
!    to J+1 (mod 3).
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) cmax
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) indx(node_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) stack(node_num)
  integer ( kind = 4 ) t
  real ( kind = 8 ) tol
  integer ( kind = 4 ) top
  integer ( kind = 4 ) element_neighbor(3,node_num*2)
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_node(3,node_num*2)

  tol = 100.0D+00 * epsilon ( tol )

  ierr = 0
!
!  Sort the vertices by increasing (x,y).
!
  call r82vec_sort_heap_index_a ( node_num, node_xy, indx )

  call r82vec_permute ( node_num, indx, node_xy )
!
!  Make sure that the data nodes are "reasonably" distinct.
!
  m1 = 1

  do i = 2, node_num

    m = m1
    m1 = i

    k = 0

    do j = 1, dim_num

      cmax = max ( abs ( node_xy(j,m) ), abs ( node_xy(j,m1) ) )

      if ( tol * ( cmax + 1.0D+00 ) &
           < abs ( node_xy(j,m) - node_xy(j,m1) ) ) then
        k = j
        exit
      end if

    end do

    if ( k == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
      write ( *, '(a,i8)' ) '  Fails for point number I = ', i
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  M1 = ', m1
      write ( *, '(a,2g14.6)' ) '  NODE_XY(M)  = ', node_xy(1:dim_num,m)
      write ( *, '(a,2g14.6)' ) '  NODE_XY(M1) = ', node_xy(1:dim_num,m1)
      ierr = 224
      stop
    end if

  end do
!
!  Starting from nodes M1 and M2, search for a third point M that
!  makes a "healthy" triangle (M1,M2,M)
!
  m1 = 1
  m2 = 2
  j = 3

  do

    if ( node_num < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
      ierr = 225
      stop
    end if

    m = j

    lr = lrline ( node_xy(1,m), node_xy(2,m), node_xy(1,m1), &
      node_xy(2,m1), node_xy(1,m2), node_xy(2,m2), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1

  end do
!
!  Set up the triangle information for (M1,M2,M), and for any other
!  triangles you created because points were collinear with M1, M2.
!
  element_num = j - 2

  if ( lr == -1 ) then

    element_node(1,1) = m1
    element_node(2,1) = m2
    element_node(3,1) = m
    element_neighbor(3,1) = -3

    do i = 2, element_num

      m1 = m2
      m2 = i+1

      element_node(1,i) = m1
      element_node(2,i) = m2
      element_node(3,i) = m

      element_neighbor(1,i-1) = -3 * i
      element_neighbor(2,i-1) = i
      element_neighbor(3,i) = i - 1

    end do

    element_neighbor(1,element_num) = -3 * element_num - 1
    element_neighbor(2,element_num) = -5
    ledg = 2
    ltri = element_num

  else

    element_node(1,1) = m2
    element_node(2,1) = m1
    element_node(3,1) = m

    element_neighbor(1,1) = -4

    do i = 2, element_num

      m1 = m2
      m2 = i+1

      element_node(1,i) = m2
      element_node(2,i) = m1
      element_node(3,i) = m

      element_neighbor(3,i-1) = i
      element_neighbor(1,i) = -3 * i - 3
      element_neighbor(2,i) = i - 1

    end do

    element_neighbor(3,element_num) = -3 * element_num
    element_neighbor(2,1) = -3 * element_num - 2
    ledg = 2
    ltri = 1

  end if
!
!  Insert the vertices one at a time from outside the convex hull,
!  determine visible boundary edges, and apply diagonal edge swaps until
!  Delaunay triangulation of vertices (so far) is obtained.
!
  top = 0

  do i = j+1, node_num

    m = i
    m1 = element_node(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = element_node(ledg+1,ltri)
    else
      m2 = element_node(1,ltri)
    end if

    lr = lrline ( node_xy(1,m), node_xy(2,m), node_xy(1,m1), &
      node_xy(2,m1), node_xy(1,m2), node_xy(2,m2), 0.0D+00 )

    if ( 0 < lr ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -element_neighbor(ledg,ltri)
      rtri = l / 3
      redg = mod ( l, 3 ) + 1
    end if

    call vbedg ( node_xy(1,m), node_xy(2,m), node_num, node_xy, element_num, &
      element_node, element_neighbor, ltri, ledg, rtri, redg )

    n = element_num + 1
    l = -element_neighbor(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -element_neighbor(e,t)
      m2 = element_node(e,t)

      if ( e <= 2 ) then
        m1 = element_node(e+1,t)
      else
        m1 = element_node(1,t)
      end if

      element_num = element_num + 1
      element_neighbor(e,t) = element_num

      element_node(1,element_num) = m1
      element_node(2,element_num) = m2
      element_node(3,element_num) = m

      element_neighbor(1,element_num) = t
      element_neighbor(2,element_num) = element_num - 1
      element_neighbor(3,element_num) = element_num + 1

      top = top + 1

      if ( node_num < top ) then
        ierr = 8
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
        write ( *, '(a)' ) '  Stack overflow.'
        stop
      end if

      stack(top) = element_num

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    element_neighbor(ledg,ltri) = -3 * n - 1
    element_neighbor(2,n) = -3 * element_num - 2
    element_neighbor(3,element_num) = -l

    ltri = n
    ledg = 2

    call swapec ( m, top, ltri, ledg, node_num, node_xy, element_num, &
      element_node, element_neighbor, stack, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
      write ( *, '(a)' ) '  Error return from SWAPEC.'
      stop
    end if

  end do
!
!  Now account for the sorting that we did.
!
  do i = 1, 3
    do j = 1, element_num
      element_node(i,j) = indx ( element_node(i,j) )
    end do
  end do

  call perm_inverse ( node_num, indx )

  call r82vec_permute ( node_num, indx, node_xy )

  return
end
subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array that has been sorted into
!    ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer variable.
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
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
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
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

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
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I
!    and J.  (Used only when the previous call returned INDX less than 0).
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
subroutine swapec ( i, top, btri, bedg, node_num, node_xy, element_num, &
  element_node, element_neighbor, stack, ierr )

!*****************************************************************************80
!
!! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!
!  Discussion:
!
!    The routine swaps diagonal edges in a 2D triangulation, based on
!    the empty circumcircle criterion, until all triangles are Delaunay,
!    given that I is the index of the new vertex added to the triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the new vertex.
!
!    Input/output, integer ( kind = 4 ) TOP, the index of the top of the stack.
!    On output, TOP is zero.
!
!    Input/output, integer ( kind = 4 ) BTRI, BEDG; on input, if positive, are
!    the triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input/output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    triangle incidence list.  May be updated on output because of swaps.
!
!    Input/output, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbor list; negative values are used for links of the
!    counter-clockwise linked list of boundary edges;  May be updated on output
!    because of swaps.
!
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer STACK(MAXST); on input, entries 1 through TOP
!    contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer ( kind = 4 ) IERR is set to 8 for abnormal return.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bedg
  integer ( kind = 4 ) btri
  integer ( kind = 4 ) c
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) e
  integer ( kind = 4 ) ee
  integer ( kind = 4 ) em1
  integer ( kind = 4 ) ep1
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fm1
  integer ( kind = 4 ) fp1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) stack(node_num)
  integer ( kind = 4 ) swap
  integer ( kind = 4 ) t
  integer ( kind = 4 ) top
  integer ( kind = 4 ) element_node(3,element_num)
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) tt
  integer ( kind = 4 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  x = node_xy(1,i)
  y = node_xy(2,i)

  do

    if ( top <= 0 ) then
      exit
    end if

    t = stack(top)
    top = top - 1

    if ( element_node(1,t) == i ) then
      e = 2
      b = element_node(3,t)
    else if ( element_node(2,t) == i ) then
      e = 3
      b = element_node(1,t)
    else
      e = 1
      b = element_node(2,t)
    end if

    a = element_node(e,t)
    u = element_neighbor(e,t)

    if ( element_neighbor(1,u) == t ) then
      f = 1
      c = element_node(3,u)
    else if ( element_neighbor(2,u) == t ) then
      f = 2
      c = element_node(1,u)
    else
      f = 3
      c = element_node(2,u)
    end if

    swap = diaedg ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,c), &
      node_xy(2,c), node_xy(1,b), node_xy(2,b) )

    if ( swap == 1 ) then

      em1 = e - 1
      em1 = i4_wrap ( em1, 1, 3 )
      ep1 = e + 1
      ep1 = i4_wrap ( ep1, 1, 3 )
      fm1 = f - 1
      fm1 = i4_wrap ( fm1, 1, 3 )
      fp1 = f + 1
      fp1 = i4_wrap ( fp1, 1, 3 )

      element_node(ep1,t) = c
      element_node(fp1,u) = i

      r = element_neighbor(ep1,t)
      s = element_neighbor(fp1,u)

      element_neighbor(ep1,t) = u
      element_neighbor(fp1,u) = t
      element_neighbor(e,t) = s
      element_neighbor(f,u) = r

      if ( 0 < element_neighbor(fm1,u) ) then
        top = top + 1
        stack(top) = u
      end if

      if ( 0 < s ) then

        if ( element_neighbor(1,s) == u ) then
          element_neighbor(1,s) = t
        else if ( element_neighbor(2,s) == u ) then
          element_neighbor(2,s) = t
        else
          element_neighbor(3,s) = t
        end if

        top = top + 1

        if ( node_num < top ) then
          ierr = 8
          return
        end if

        stack(top) = t

      else

        if ( u == btri .and. fp1 == bedg ) then
          btri = t
          bedg = e
        end if

        l = - ( 3 * t + e - 1 )
        tt = t
        ee = em1

        do while ( 0 < element_neighbor(ee,tt) )

          tt = element_neighbor(ee,tt)

          if ( element_node(1,tt) == a ) then
            ee = 3
          else if ( element_node(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        element_neighbor(ee,tt) = l

      end if

      if ( 0 < r ) then

        if ( element_neighbor(1,r) == t ) then
          element_neighbor(1,r) = u
        else if ( element_neighbor(2,r) == t ) then
          element_neighbor(2,r) = u
        else
          element_neighbor(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( 0 < element_neighbor(ee,tt) )

          tt = element_neighbor(ee,tt)

          if ( element_node(1,tt) == b ) then
            ee = 3
          else if ( element_node(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        element_neighbor(ee,tt) = l

      end if

    end if

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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2005
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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine triangle_angles_2d ( t, angle )

!*****************************************************************************80
!
!! TRIANGLE_ANGLES_2D computes the angles of a triangle in 2D.
!
!  Discussion:
!
!    The law of cosines is used:
!
!      C^2 = A^2 + B^2 - 2 * A * B * COS ( GAMMA )
!
!    where GAMMA is the angle opposite side C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) ANGLE(3), the angles opposite
!    sides P1-P2, P2-P3 and P3-P1, in radians.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a
  real ( kind = 8 ) angle(3)
  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t(dim_num,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:dim_num,1) - t(1:dim_num,2) )**2 ) )
  b = sqrt ( sum ( ( t(1:dim_num,2) - t(1:dim_num,3) )**2 ) )
  c = sqrt ( sum ( ( t(1:dim_num,3) - t(1:dim_num,1) )**2 ) )
!
!  Take care of ridiculous special cases.
!
  if ( a == 0.0D+00 .and. b == 0.0D+00 .and. c == 0.0D+00 ) then
    angle(1:3) = 2.0D+00 * pi / 3.0D+00
    return
  end if

  if ( c == 0.0D+00 .or. a == 0.0D+00 ) then
    angle(1) = pi
  else
    angle(1) = arc_cosine ( ( c * c + a * a - b * b ) / ( 2.0D+00 * c * a ) )
  end if

  if ( a == 0.0D+00 .or. b == 0.0D+00 ) then
    angle(2) = pi
  else
    angle(2) = arc_cosine ( ( a * a + b * b - c * c ) / ( 2.0D+00 * a * b ) )
  end if

  if ( b == 0.0D+00 .or. c == 0.0D+00 ) then
    angle(3) = pi
  else
    angle(3) = arc_cosine ( ( b * b + c * c - a * a ) / ( 2.0D+00 * b * c ) )
  end if

  return
end
subroutine triangle_area_2d ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA_2D computes the area of a triangle in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) AREA, the absolute area of the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ) t(dim_num,3)

  area = 0.5D+00 * abs ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )

  return
end
subroutine triangle_circumcenter_2d ( t, center )

!*****************************************************************************80
!
!! TRIANGLE_CIRCUMCENTER_2D computes the circumcenter of a triangle in 2D.
!
!  Discussion:
!
!    The circumcenter of a triangle is the center of the circumcircle, the
!    circle that passes through the three vertices of the triangle.
!
!    The circumcircle contains the triangle, but it is not necessarily the
!    smallest triangle to do so.
!
!    If all angles of the triangle are no greater than 90 degrees, then
!    the center of the circumscribed circle will lie inside the triangle.
!    Otherwise, the center will lie outside the triangle.
!
!    The circumcenter is the intersection of the perpendicular bisectors
!    of the sides of the triangle.
!
!    In geometry, the circumcenter of a triangle is often symbolized by "O".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) CENTER(2), the circumcenter of the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) asq
  real ( kind = 8 ) bot
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) csq
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ) top(dim_num)

  asq = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2
  csq = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2

  top(1) =    ( t(2,2) - t(2,1) ) * csq - ( t(2,3) - t(2,1) ) * asq
  top(2) =  - ( t(1,2) - t(1,1) ) * csq + ( t(1,3) - t(1,1) ) * asq

  bot  =  ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) ) &
        - ( t(2,3) - t(2,1) ) * ( t(1,2) - t(1,1) )

  center(1:2) = t(1:2,1) + 0.5D+00 * top(1:2) / bot

  return
end
subroutine element_order3_physical_to_reference ( t, n, phy, ref )

!*****************************************************************************80
!
!! ELEMENT_ORDER3_PHYSICAL_TO_REFERENCE maps T3 physical to reference points.
!
!  Discussion:
!
!    Given the vertices of an order 3 physical triangle and a point
!    (X,Y) in the physical triangle, the routine computes the value
!    of the corresponding image point (XSI,ETA) in reference space.
!
!    This routine is also appropriate for an order 4 triangle, assuming
!    that the fourth node is always the centroid of the triangle.
!
!    This routine may be appropriate for an order 6
!    triangle, if the mapping between reference and physical space
!    is linear.  This implies, in particular, that the sides of the
!    image triangle are straight and that the "midside" nodes in the
!    physical triangle are halfway along the sides of
!    the physical triangle.
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the X and Y coordinates
!    of the vertices.  The vertices are assumed to be the images of
!    (0,0), (1,0) and (0,1) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) PHY(2,N), the coordinates of physical points
!    to be transformed.
!
!    Output, real ( kind = 8 ) REF(2,N), the coordinates of the corresponding
!    points in the reference space.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) phy(2,n)
  real ( kind = 8 ) ref(2,n)
  real ( kind = 8 ) t(2,3)

  ref(1,1:n) = ( ( t(2,3) - t(2,1) ) * ( phy(1,1:n) - t(1,1) )   &
               - ( t(1,3) - t(1,1) ) * ( phy(2,1:n) - t(2,1) ) ) &
             / ( ( t(2,3) - t(2,1) ) * ( t(1,2)     - t(1,1) )   &
               - ( t(1,3) - t(1,1) ) * ( t(2,2)     - t(2,1) ) )

  ref(2,1:n) = ( ( t(1,2) - t(1,1) ) * ( phy(2,1:n) - t(2,1) )   &
               - ( t(2,2) - t(2,1) ) * ( phy(1,1:n) - t(1,1) ) ) &
             / ( ( t(2,3) - t(2,1) ) * ( t(1,2)     - t(1,1) )   &
               - ( t(1,3) - t(1,1) ) * ( t(2,2)     - t(2,1) ) )

  return
end
subroutine element_order3_reference_to_physical ( t, n, ref, phy )

!*****************************************************************************80
!
!! ELEMENT_ORDER3_REFERENCE_TO_PHYSICAL maps T3 reference to physical points.
!
!  Discussion:
!
!    Given the vertices of an order 3 physical triangle and a point
!    (XSI,ETA) in the reference triangle, the routine computes the value
!    of the corresponding image point (X,Y) in physical space.
!
!    This routine is also appropriate for an order 4 triangle,
!    as long as the fourth node is the centroid of the triangle.
!
!    This routine may also be appropriate for an order 6
!    triangle, if the mapping between reference and physical space
!    is linear.  This implies, in particular, that the sides of the
!    image triangle are straight and that the "midside" nodes in the
!    physical triangle are halfway along the sides of
!    the physical triangle.
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
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
!    Input, real ( kind = 8 ) T(2,3), the coordinates of the vertices.
!    The vertices are assumed to be the images of (0,0), (1,0) and
!    (0,1) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) REF(2,N), points in the reference triangle.
!
!    Output, real ( kind = 8 ) PHY(2,N), corresponding points in the
!    physical triangle.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) phy(2,n)
  real ( kind = 8 ) ref(2,n)
  real ( kind = 8 ) t(2,3)

  do i = 1, 2
    phy(i,1:n) = t(i,1) * ( 1.0D+00 - ref(1,1:n) - ref(2,1:n) ) &
               + t(i,2) *             ref(1,1:n)                &
               + t(i,3) *                          ref(2,1:n)
  end do

  return
end
subroutine element_order6_physical_to_reference ( t, n, phy, ref )

!*****************************************************************************80
!
!! ELEMENT_ORDER6_PHYSICAL_TO_REFERENCE maps T6  physical to reference points.
!
!  Discussion:
!
!    Given the vertices of an order 6 physical triangle and a point
!    (X,Y) in the physical triangle, the routine computes the value
!    of the corresponding image point (R,S) in reference space.
!
!    The mapping from (R,S) to (X,Y) has the form:
!
!      X(R,S) = A1 * R * R + B1 * R * S + C1 * S * S
!             + D1 * R     + E1 * S     + F1
!
!      Y(R,S) = A2 * R * R + B2 * R * S + C2 * S * S
!             + D2 * R     + E2 * S     + F2
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,6), the coordinates of the vertices.
!    The vertices are assumed to be the images of (0,0), (1,0), (0,1),
!    (1/2,0), (1/2,1/2) and (0,1/2), in that order.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) PHY(2,N), the coordinates of points in the
!    physical space.
!
!    Output, real ( kind = 8 ) REF(2,N), the coordinates of the corresponding
!    points in the reference space.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b(2)
  real ( kind = 8 ) c(2)
  real ( kind = 8 ) d(2)
  real ( kind = 8 ) det
  real ( kind = 8 ) dx(2)
  real ( kind = 8 ) e(2)
  real ( kind = 8 ) f(2)
  real ( kind = 8 ) fun(2)
  real ( kind = 8 ) fun_norm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(2,2)
  integer ( kind = 4 ), parameter :: it_max = 10
  real ( kind = 8 ), parameter :: it_tol = 0.000001D+00
  real ( kind = 8 ) phy(2,n)
  real ( kind = 8 ) ref(2,n)
  real ( kind = 8 ) t(2,6)
!
!  Set iteration parameters.
!
  do i = 1, 2

    a(i) =   2.0D+00 * t(i,1) + 2.0D+00 * t(i,2)                    &
           - 4.0D+00 * t(i,4)

    b(i) =   4.0D+00 * t(i,1)                                       &
           - 4.0D+00 * t(i,4) + 4.0D+00 * t(i,5) - 4.0D+00 * t(i,6)

    c(i) =   2.0D+00 * t(i,1)                    + 2.0D+00 * t(i,3) &
                                                 - 4.0D+00 * t(i,6)

    d(i) = - 3.0D+00 * t(i,1) -           t(i,2)                    &
           + 4.0D+00 * t(i,4)

    e(i) = - 3.0D+00 * t(i,1)                    -           t(i,3) &
                                                 + 4.0D+00 * t(i,6)
    f(i) =             t(i,1)

  end do
!
!  Initialize the points by inverting the linear map.
!
  call element_order3_physical_to_reference ( t(1:2,1:3), n, phy, ref )
!
!  Carry out the Newton iteration.
!
  do j = 1, n

    do it = 1, it_max

      fun(1:2) = a(1:2) * ref(1,j) * ref(1,j) &
               + b(1:2) * ref(1,j) * ref(2,j) &
               + c(1:2) * ref(2,j) * ref(2,j) &
               + d(1:2) * ref(1,j) &
               + e(1:2) * ref(2,j) &
               + f(1:2) &
               - phy(1:2,j)

      fun_norm = sqrt ( fun(1) * fun(1) + fun(2) * fun(2) )

      if ( fun_norm <= it_tol ) then
        exit
      end if

      jac(1:2,1) = 2.0D+00 * a(1:2) * ref(1,j) &
                 +           b(1:2) * ref(2,j) + d(1:2)

      jac(1:2,2) =           b(1:2) * ref(1,j) &
                 + 2.0D+00 * c(1:2) * ref(2,j) + e(1:2)

      det = jac(1,1) * jac(2,2) - jac(1,2) * jac(2,1)

      if ( det == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) &
          'ELEMENT_ORDER6_PHYSICAL_TO_REFERENCE - Fatal error!'
        write ( *, '(a)' ) '  The jacobian of the mapping is singular.'
      end if

      dx(1) = (  jac(2,2) * fun(1) - jac(1,2) * fun(2) ) / det
      dx(2) = ( -jac(2,1) * fun(1) + jac(1,1) * fun(2) ) / det

      ref(1:2,j) = ref(1:2,j) - dx(1:2)

    end do

  end do

  return
end
subroutine element_order6_reference_to_physical ( t, n, ref, phy )

!*****************************************************************************80
!
!! ELEMENT_ORDER6_REFERENCE_TO_PHYSICAL maps T6 reference to physical points.
!
!  Discussion:
!
!    Given the vertices of an order 6 physical triangle and a point
!    (XSI,ETA) in the reference triangle, the routine computes the value
!    of the corresponding image point (X,Y) in physical space.
!
!    The mapping from (XSI,ETA) to (X,Y) has the form:
!
!      X(ETA,XSI) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
!                 + D1 * XSI    + E1 * ETA     + F1
!
!      Y(ETA,XSI) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
!                 + D2 * XSI    + E2 * ETA     + F2
!
!  Reference Element T6:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,6), the coordinates of the vertices.
!    The vertices are assumed to be the images of (0,0), (1,0),
!    (0,1),(1/2,0), (1/2,1/2) and (0,1/2) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) REF(2,N), points in the reference triangle.
!
!    Output, real ( kind = 8 ) PHY(2,N), corresponding points in the
!    physical triangle.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b(2)
  real ( kind = 8 ) c(2)
  real ( kind = 8 ) d(2)
  real ( kind = 8 ) e(2)
  real ( kind = 8 ) f(2)
  integer ( kind = 4 ) i
  real ( kind = 8 ) phy(2,n)
  real ( kind = 8 ) ref(2,n)
  real ( kind = 8 ) t(2,6)

  do i = 1, 2

    a(i) =   2.0D+00 * t(i,1) + 2.0D+00 * t(i,2)                    &
           - 4.0D+00 * t(i,4)

    b(i) =   4.0D+00 * t(i,1)                                       &
           - 4.0D+00 * t(i,4) + 4.0D+00 * t(i,5) - 4.0D+00 * t(i,6)

    c(i) =   2.0D+00 * t(i,1)                    + 2.0D+00 * t(i,3) &
                                                 - 4.0D+00 * t(i,6)

    d(i) = - 3.0D+00 * t(i,1) -           t(i,2)                    &
           + 4.0D+00 * t(i,4)

    e(i) = - 3.0D+00 * t(i,1)                    -           t(i,3) &
                                                 + 4.0D+00 * t(i,6)
    f(i) =             t(i,1)

  end do

  do i = 1, 2
    phy(i,1:n) = a(i) * ref(1,1:n) * ref(1,1:n) &
               + b(i) * ref(1,1:n) * ref(2,1:n) &
               + c(i) * ref(2,1:n) * ref(2,1:n) &
               + d(i) * ref(1,1:n) &
               + e(i) * ref(2,1:n) &
               + f(i)
  end do

  return
end
subroutine triangle_reference_sample ( n, seed, p )

!*****************************************************************************80
!
!! TRIANGLE_REFERENCE_SAMPLE returns random points in the reference triangle.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  |   \  side 2
!       |    \
!    3  |     \
!       |      \
!       1-------2
!
!         side 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points to generate.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) P(2,N), random points in the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(dim_num,n)
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  do j = 1, n

    r = r8_uniform_01 ( seed )
!
!  Interpret R as a percentage of the triangle's area.
!
!  Imagine a line L, parallel to side 1, so that the area between
!  vertex 1 and line L is R percent of the full triangle's area.
!
!  The line L will intersect sides 2 and 3 at a fraction
!  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
!
    alpha = sqrt ( r )
!
!  Now choose, uniformly at random, a point on the line L.
!
    beta = r8_uniform_01 ( seed )

    p(1,j) = ( 1.0D+00 - beta ) * alpha
    p(2,j) =             beta   * alpha

  end do

  return
end
subroutine triangle_sample ( t, n, seed, p )

!*****************************************************************************80
!
!! TRIANGLE_SAMPLE returns random points in a triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Input, integer ( kind = 4 ) N, the number of points to generate.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) P(2,N), random points in the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha(n)
  integer ( kind = 4 ) dim
  real ( kind = 8 ) p(dim_num,n)
  real ( kind = 8 ) p12(dim_num,n)
  real ( kind = 8 ) p13(dim_num,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(dim_num,3)
!
!  For comparison between F90, C++ and MATLAB codes, call R8VEC_UNIFORM_01.
!
  call r8vec_uniform_01 ( n, seed, alpha )
!
!  Interpret R as a percentage of the triangle's area.
!
!  Imagine a line L, parallel to side 1, so that the area between
!  vertex 1 and line L is R percent of the full triangle's area.
!
!  The line L will intersect sides 2 and 3 at a fraction
!  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
!
  alpha(1:n) = sqrt ( alpha(1:n) )
!
!  Determine the coordinates of the points on sides 2 and 3 intersected
!  by line L.
!
  do dim = 1, dim_num

    p12(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * t(dim,1) &
                             + alpha(1:n)   * t(dim,2)

    p13(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * t(dim,1) &
                             + alpha(1:n)   * t(dim,3)

  end do
!
!  Now choose, uniformly at random, a point on the line L.
!
  call r8vec_uniform_01 ( n, seed, alpha )

  do dim = 1, dim_num

    p(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * p12(dim,1:n) &
                           + alpha(1:n)   * p13(dim,1:n)

  end do

  return
end
function triangulation_area ( node_num, node_xy, element_order, &
  element_num, element_node )

!*****************************************************************************80
!
!! TRIANGULATION_AREA computes the area of a triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
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
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the nodes making up each triangle.
!
!    Output, real ( kind = 8 ) TRIANGULATION_AREA, the area.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) element
  real ( kind = 8 ) element_area
  integer ( kind = 4 ) element_node(element_order,element_num)
  real ( kind = 8 ) element_xy(2,3)
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) triangulation_area
  real ( kind = 8 ) value

  value = 0.0D+00

  do element = 1, element_num

    element_xy(1:2,1:3) = node_xy(1:2,element_node(1:3,element))

    call triangle_area_2d ( element_xy, element_area )

    value = value + element_area

  end do

  triangulation_area = value

  return
end
subroutine triangulation_areas ( node_num, node_xy, element_order, &
  element_num, element_node, triangle_area, triangulation_area )

!*****************************************************************************80
!
!! TRIANGULATION_AREAS computes triangle and triangulation areas.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes in the
!    triangulation.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of triangles in
!    the triangulation.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles in
!    the triangulation.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the nodes making up each triangle.
!
!    Output, real ( kind = 8 ) TRIANGLE_AREA(ELEMENT_NUM), the area of
!    the triangles.
!
!    Output, real ( kind = 8 ) TRIANGULATION_AREA, the area of
!    the triangulation.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  real ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) triangle
  real ( kind = 8 ) triangle_area(element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)
  real ( kind = 8 ) triangle_xy(2,3)
  real ( kind = 8 ) triangulation_area

  triangulation_area = 0.0D+00

  do triangle = 1, element_num

    triangle_xy(1:2,1:3) = node_xy(1:2,element_node(1:3,triangle))

    call triangle_area_2d ( triangle_xy, triangle_area(triangle) )

    triangulation_area = triangulation_area + triangle_area(triangle)

  end do

  return
end
subroutine triangulation_delaunay_discrepancy_compute ( node_num, node_xy, &
  element_order, element_num, element_node, element_neighbor, &
  angle_min, angle_min_triangle, angle_max, angle_max_triangle, value )

!*****************************************************************************80
!
!! TRIANGULATION_DELAUNAY_DISCREPANCY_COMPUTE: is a triangulation Delaunay?
!
!  Discussion:
!
!    A (maximal) triangulation is Delaunay if and only if it is locally
!    Delaunay.
!
!    A triangulation is Delaunay if the minimum angle over all triangles
!    in the triangulation is maximized.  That is, there is no other
!    triangulation of the points which has a larger minimum angle.
!
!    A triangulation is locally Delaunay if, for every pair of triangles
!    that share an edge, the minimum angle in the two triangles is larger
!    than the minimum angle in the two triangles formed by removing the
!    common edge and joining the two opposing vertices.
!
!    This function examines the question of whether a given triangulation
!    is locally Delaunay.  It does this by looking at every pair of
!    neighboring triangles and comparing the minimum angle attained
!    for the current triangle pair and the alternative triangle pair.
!
!    Let A(i,j) be the minimum angle formed by triangles T(i) and T(j),
!    which are two triangles in the triangulation which share a common edge.
!    Let B(I,J) be the minimum angle formed by triangles S(i) and S(j),
!    where S(i) and S(j) are formed by removing the common edge of T(i)
!    and T(j), and joining the opposing vertices.
!
!    Then the triangulation is Delaunay if B(i,j) <= A(i,j) for every
!    pair of neighbors T(i) and T(j).
!
!    If A(i,j) < B(i,j) for at least one pair of neighbors, the triangulation
!    is not a Delaunay triangulation.
!
!    This program returns VALUE = min ( A(i,j) - B(i,j) ) over all
!    triangle neighbors.  VALUE is scaled to be in degrees, for
!    comprehensibility.  If VALUE is negative, then at least one pair
!    of triangles violates the Delaunay condition, and so the entire
!    triangulation is not a Delaunay triangulation.  If VALUE is nonnegative,
!    then the triangulation is a Delaunay triangulation.
!
!    It is useful to return VALUE, rather than a simple True/False value,
!    because there can be cases where the Delaunay condition is only
!    "slightly" violated.  A simple example is a triangulation formed
!    by starting with a mesh of squares and dividing each square into
!    two triangles by choosing one of the diagonals of the square.
!    The Delaunay discrepancy for this mesh, if computed exactly, is 0,
!    but roundoff could easily result in discrepancies that were very
!    slightly negative.
!
!    If VALUE is positive, and not very small in magnitude, then every
!    pair of triangles in the triangulation satisfies the local Delaunay
!    condition, and so the triangulation is a Delaunay triangulation.
!
!    If VALUE is negative, and not very small in magnitude, then at least
!    one pair of triangles violates the Delaunay condition, and to a
!    significant degree.  The triangulation is not a Delaunay triangulation.
!
!    If the magnitude of VALUE is very close to zero, then the triangulation
!    is numerically ambiguous.  At least one pair of triangles violates
!    or almost violates the condition, but no triangle violates the
!    condition to a great extent.  The user must judge whether the
!    violation is significant or not.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2009
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles in
!    the triangulation.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the nodes that make up each triangle.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbor list.
!
!    Output, real ( kind = 8 ) ANGLE_MIN, the minimum angle that occurred in
!    the triangulation.
!
!    Output, integer ( kind = 4 ) ANGLE_MIN_TRIANGLE, the triangle in which
!    the minimum angle occurred.
!
!    Output, real ( kind = 8 ) ANGLE_MAX, the maximum angle that occurred in
!    the triangulation.
!
!    Output, integer ( kind = 4 ) ANGLE_MAX_TRIANGLE, the triangle in which
!    the maximum angle occurred.
!
!    Output, real ( kind = 8 ) VALUE, the minimum value of ( A(i,j) - B(i,j) ).
!    POSITIVE indicates the triangulation is Delaunay.
!    VERY NEAR ZERO is a numerically ambiguous case.
!    NEGATIVE indicates the triangulation is not Delaunay.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  real ( kind = 8 ) angle_max
  integer ( kind = 4 ) angle_max_triangle
  real ( kind = 8 ) angle_min
  integer ( kind = 4 ) angle_min_triangle
  real ( kind = 8 ) angle_min1
  real ( kind = 8 ) angle_min2
  real ( kind = 8 ) angles1(3)
  real ( kind = 8 ) angles2(3)
  integer ( kind = 8 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) neighbor
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) t(2,3)
  integer ( kind = 4 ) triangle_index
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) triangle1
  integer ( kind = 4 ) triangle2
  real ( kind = 8 ) value

  angle_max = 0.0D+00
  angle_max_triangle = - 1
  angle_min = pi
  angle_min_triangle = -1
  value = 0.0D+00
!
!  Consider triangle TRIANGLE1
!
  do triangle1 = 1, element_num
!
!  Consider the side opposite vertex NEIGHBOR.
!
    do neighbor = 1, 3

      triangle2 = element_neighbor(neighbor,triangle1)
!
!  There might be no neighbor on side NEIGHBOR.
!
      if ( triangle2 < 0 ) then
        cycle
      end if
!
!  We only need to check a pair of triangles once.
!
      if ( triangle2 < triangle1 ) then
        cycle
      end if
!
!  List the vertices of the quadrilateral in such a way
!  that the nodes of triangle 1 come first.
!
!  We rely on a property of the ELEMENT_NEIGHBOR array, namely, that
!  neighbor #1 is on the side opposite to vertex #1, and so on.
!
      i1 = i4_wrap ( neighbor + 2, 1, 3 )
      i2 = i4_wrap ( neighbor,     1, 3 )
      i3 = i4_wrap ( neighbor + 1, 1, 3 )

      n1 = element_node(i1,triangle1)
      n2 = element_node(i2,triangle1)
      n3 = element_node(i3,triangle1)
!
!  The "odd" or "opposing" node of the neighboring triangle
!  is the one which follows common node I3.
!
      n4 = -1
      do i = 1, 3
        if ( element_node(i,triangle2) == n3 ) then
          i4 = i + 1
          i4 = i4_wrap ( i4, 1, 3 )
          n4 = element_node(i4,triangle2)
          exit
        end if
      end do

      if ( n4 == -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) &
          'TRIANGULATION_DELAUNAY_DISCREPANCY_COMPUTE - Fatal error!'
        write ( *, '(a)' ) '  Could not identify the fourth node.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Triangle1 = ', triangle1
        write ( *, '(a,3i8)' ) '  Nodes =     ', element_node(1:3,triangle1)
        write ( *, '(a,3i8)' ) '  Neighbors = ', &
         element_neighbor(1:3,triangle1)
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Neighbor index = ', neighbor
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Triangle2 = ', triangle2
        write ( *, '(a,3i8)' ) '  Nodes =     ', element_node(1:3,triangle2)
        write ( *, '(a,3i8)' ) '  Neighbors = ', &
          element_neighbor(1:3,triangle2)
        stop
      end if
!
!  Compute the minimum angle for (I1,I2,I3) and (I1,I3,I4).
!
      t(1:2,1) = node_xy(1:2,n1)
      t(1:2,2) = node_xy(1:2,n2)
      t(1:2,3) = node_xy(1:2,n3)
      call triangle_angles_2d ( t, angles1 )

      t(1:2,1) = node_xy(1:2,n1)
      t(1:2,2) = node_xy(1:2,n3)
      t(1:2,3) = node_xy(1:2,n4)
      call triangle_angles_2d ( t, angles2 )

      angle_min1 = min ( minval ( angles1 ), minval ( angles2 ) )

      if ( angle_max < maxval ( angles1 ) ) then
        angle_max = maxval ( angles1 )
        angle_max_triangle = triangle1
      end if

      if ( angle_max < maxval ( angles2 ) ) then
        angle_max = maxval ( angles2 )
        angle_max_triangle = triangle2
      end if

      if ( minval ( angles1 ) < angle_min ) then
        angle_min = minval ( angles1 )
        angle_min_triangle = triangle1
      end if

      if ( minval ( angles2 ) < angle_min ) then
        angle_min = minval ( angles2 )
        angle_min_triangle = triangle2
      end if
!
!  Compute the minimum angle for (I1,I2,I4) and (I2,I3,I4).
!
      t(1:2,1) = node_xy(1:2,n1)
      t(1:2,2) = node_xy(1:2,n2)
      t(1:2,3) = node_xy(1:2,n4)
      call triangle_angles_2d ( t, angles1 )

      t(1:2,1) = node_xy(1:2,n2)
      t(1:2,2) = node_xy(1:2,n3)
      t(1:2,3) = node_xy(1:2,n4)
      call triangle_angles_2d ( t, angles2 )

      angle_min2 = min ( minval ( angles1 ), minval ( angles2 ) )
!
!  Compare this value to the current minimum.
!
      value = min ( value, angle_min1 - angle_min2 )

    end do

  end do
!
!  Scale the results to degrees.
!
  value = value * 180.0D+00 / pi
  angle_max = angle_max * 180.0D+00 / pi
  angle_min = angle_min * 180.0D+00 / pi

  return
end
subroutine triangulation_neighbor_elements ( element_order, element_num, &
  element_node, element_neighbor )

!*****************************************************************************80
!
!! TRIANGULATION_NEIGHBOR_TRIANGLES determines element neighbors.
!
!  Discussion:
!
!    A triangulation of a set of nodes can be completely described by
!    the coordinates of the nodes, and the list of nodes that make up
!    each triangular element.  However, in some cases, it is necessary to know
!    element adjacency information, that is, which element, if any,
!    is adjacent to a given element on a particular side.
!
!    This routine creates a data structure recording this information.
!
!    The primary amount of work occurs in sorting a list of 3 * ELEMENT_NUM
!    data items.
!
!    Note that ROW is a work array allocated dynamically inside this
!    routine.  It is possible, for very large values of ELEMENT_NUM,
!    that the necessary amount of memory will not be accessible, and the
!    routine will fail.  This is a limitation of the implementation of
!    dynamic arrays in FORTRAN90.  One way to get around this would be
!    to require the user to declare ROW in the calling routine
!    as an allocatable array, get the necessary memory explicitly with
!    an ALLOCATE statement, and then pass ROW into this routine.
!
!    Of course, the point of dynamic arrays was to make it easy to
!    hide these sorts of temporary work arrays from the poor user!
!
!    This routine was revised to store the edge data in a column
!    array rather than a row array.
!
!  Example:
!
!    The input information from ELEMENT_NODE:
!
!    Element    Nodes
!    --------   ---------------
!     1         3      4      1
!     2         3      1      2
!     3         3      2      8
!     4         2      1      5
!     5         8      2     13
!     6         8     13      9
!     7         3      8      9
!     8        13      2      5
!     9         9     13      7
!    10         7     13      5
!    11         6      7      5
!    12         9      7      6
!    13        10      9      6
!    14         6      5     12
!    15        11      6     12
!    16        10      6     11
!
!    The output information in ELEMENT_NEIGHBOR:
!
!    Element   Neighboring Elements
!    --------  ---------------------
!
!     1        -1     -1      2
!     2         1      4      3
!     3         2      5      7
!     4         2     -1      8
!     5         3      8      6
!     6         5      9      7
!     7         3      6     -1
!     8         5      4     10
!     9         6     10     12
!    10         9      8     11
!    11        12     10     14
!    12         9     11     13
!    13        -1     12     16
!    14        11     -1     15
!    15        16     14     -1
!    16        13     15     -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the nodes that make up each element.
!
!    Output, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the three
!    elements that are direct neighbors of a given element.
!    ELEMENT_NEIGHBOR(1,I) is the index of the element which touches side 1,
!    defined by nodes 2 and 3, and so on.  ELEMENT_NEIGHBOR(1,I) is negative
!    if there is no neighbor on that side.  In this case, that side of the
!    element lies on a boundary of the triangulation.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) col(4,3*element_num)
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) element1
  integer ( kind = 4 ) element2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) side1
  integer ( kind = 4 ) side2
!
!  Step 1.
!  From the list of nodes for element E, of the form: (I,J,K)
!  construct the three neighbor relations:
!
!    (I,J,3,E) or (J,I,3,E),
!    (J,K,1,E) or (K,J,1,E),
!    (K,I,2,E) or (I,K,2,E)
!
!  where we choose (I,J,1,E) if I < J, or else (J,I,1,E)
!
  do element = 1, element_num

    i = element_node(1,element)
    j = element_node(2,element)
    k = element_node(3,element)

    if ( i < j ) then
      col(1:4,3*(element-1)+1) = (/ i, j, 3, element /)
    else
      col(1:4,3*(element-1)+1) = (/ j, i, 3, element /)
    end if

    if ( j < k ) then
      col(1:4,3*(element-1)+2) = (/ j, k, 1, element /)
    else
      col(1:4,3*(element-1)+2) = (/ k, j, 1, element /)
    end if

    if ( k < i ) then
      col(1:4,3*(element-1)+3) = (/ k, i, 2, element /)
    else
      col(1:4,3*(element-1)+3) = (/ i, k, 2, element /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1 and 2; the routine we call here
!  sorts on rows 1 through 4 but that won't hurt us.
!
!  What we need is to find cases where two elements share an edge.
!  Say they share an edge defined by the nodes I and J.  Then there are
!  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
!  we make sure that these two columns occur consecutively.  That will
!  make it easy to notice that the triangles are neighbors.
!
  call i4col_sort_a ( 4, 3 * element_num, col )
!
!  Step 3. Neighboring elements show up as consecutive columns with
!  identical first two entries.  Whenever you spot this happening,
!  make the appropriate entries in ELEMENT_NEIGHBOR.
!
  element_neighbor(1:3,1:element_num) = -1

  icol = 1

  do

    if ( 3 * element_num <= icol ) then
      exit
    end if

    if ( col(1,icol) /= col(1,icol+1) .or. col(2,icol) /= col(2,icol+1) ) then
      icol = icol + 1
      cycle
    end if

    side1 = col(3,icol)
    element1 = col(4,icol)
    side2 = col(3,icol+1)
    element2 = col(4,icol+1)

    element_neighbor(side1,element1) = element2
    element_neighbor(side2,element2) = element1

    icol = icol + 2

  end do

  return
end
subroutine triangulation_node_order ( element_order, element_num, &
  element_node, node_num, node_order )

!*****************************************************************************80
!
!! TRIANGULATION_NODE_ORDER determines the order of nodes in a triangulation.
!
!  Discussion:
!
!    The order of a node is the number of triangles that use that node
!    as a vertex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the nodes that make up the triangles.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) NODE_ORDER(NODE_NUM), the order of each node.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) i
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_order(node_num)
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) element_node(element_order,element_num)

  node_order(1:node_num) = 0

  do triangle = 1, element_num
    do i = 1, element_order
      node = element_node(i,triangle)
      if ( node < 1 .or. node_num < node ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGULATION_NODE_ORDER - Fatal error!'
        write ( *, '(a)' ) '  Illegal entry in ELEMENT_NODE.'
        stop
      else
        node_order(node) = node_order(node) + 1
      end if
    end do
  end do

  return
end
subroutine triangulation_order3_adj_count ( node_num, element_num, &
  element_node, element_neighbor, adj_num, adj_col )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The triangulation is assumed to involve 3-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!    Also, a node is considered to be adjacent to itself.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  |   \  side 2
!       |    \
!    3  |     \
!       |      \
!       1-------2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!    A sample grid.
!
!
!    Below, we have a chart that summarizes the adjacency relationships
!    in the sample grid.  On the left, we list the node, and its neighbors,
!    with an asterisk to indicate the adjacency of the node to itself
!    (in some cases, you want to count this self adjacency and in some
!    you don't).  On the right, we list the number of adjacencies to
!    lower-indexed nodes, to the node itself, to higher-indexed nodes,
!    the total number of adjacencies for this node, and the location
!    of the first and last entries required to list this set of adjacencies
!    in a single list of all the adjacencies.
!
!    N   Adjacencies                Below  Self   Above   Total First  Last
!
!   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
!    1:  *  2  6                        0     1       2       3     1     3
!    2:  1  *  3  6  7                  1     1       3       5     4     8
!    3:  2  *  4  7  8                  1     1       3       5     9    13
!    4:  3  *  5  8  9                  1     1       3       5    14    18
!    5:  4  *  9 10                     1     1       2       4    19    22
!    6:  1  2  *  7 11                  2     1       2       5    23    27
!    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
!    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
!    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
!   10:  5  9  * 14 15                  2     1       2       5    49    53
!   11:  6  7  * 12 16                  2     1       2       5    54    58
!   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
!   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
!   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
!   15: 10 14  * 19 20                  2     1       2       5    80    84
!   16: 11 12  * 17 21                  2     1       2       5    85    89
!   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
!   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
!   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
!   20: 15 19  * 24 25                  2     1       2       5   111   115
!   21: 16 17  * 22                     2     1       1       4   116   119
!   22: 17 18 21  * 23                  3     1       1       5   120   124
!   23: 18 19 22  * 24                  3     1       1       5   125   129
!   24: 19 20 23  * 25                  3     1       1       5   130   134
!   25: 20 24  *                        2     1       0       3   135   137
!   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), lists the
!    nodes that make up each triangle, in counterclockwise order.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), for each
!    side of a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Output, integer ( kind = 4 ) ADJ_NUM, the number of adjacencies.
!
!    Output, integer ( kind = 4 ) ADJ_COL(NODE_NUM+1).  Information about
!    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_col(node_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle2
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)

  adj_num = 0
!
!  Set every node to be adjacent to itself.
!
  adj_col(1:node_num) = 1
!
!  Examine each triangle.
!
  do triangle = 1, element_num

    n1 = element_node(1,triangle)
    n2 = element_node(2,triangle)
    n3 = element_node(3,triangle)
!
!  Add edge (1,2) if this is the first occurrence,
!  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
    triangle2 = element_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_col(n1) = adj_col(n1) + 1
      adj_col(n2) = adj_col(n2) + 1
    end if
!
!  Add edge (2,3).
!
    triangle2 = element_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_col(n2) = adj_col(n2) + 1
      adj_col(n3) = adj_col(n3) + 1
    end if
!
!  Add edge (3,1).
!
    triangle2 = element_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_col(n1) = adj_col(n1) + 1
      adj_col(n3) = adj_col(n3) + 1
    end if

  end do
!
!  We used ADJ_COL to count the number of entries in each column.
!  Convert it to pointers into the ADJ array.
!
  adj_col(2:node_num+1) = adj_col(1:node_num)

  adj_col(1) = 1
  do i = 2, node_num+1
    adj_col(i) = adj_col(i-1) + adj_col(i)
  end do

  adj_num = adj_col(node_num+1) - 1

  return
end
subroutine triangulation_order3_adj_set ( node_num, element_num, &
  element_node, element_neighbor, adj_num, adj_col, adj )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_ADJ_SET sets adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to set the adjacencies, after the
!    appropriate amount of memory has been set aside.
!
!    The triangulation is assumed to involve 3-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!    Also, a node is considered to be adjacent to itself.
!
!    This routine can be used to create the compressed column storage
!    for a linear triangle finite element discretization of
!    Poisson's equation in two dimensions.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  |   \  side 2
!       |    \
!    3  |     \
!       |      \
!       1-------2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!    A sample grid
!
!
!    Below, we have a chart that summarizes the adjacency relationships
!    in the sample grid.  On the left, we list the node, and its neighbors,
!    with an asterisk to indicate the adjacency of the node to itself
!    (in some cases, you want to count this self adjacency and in some
!    you don't).  On the right, we list the number of adjacencies to
!    lower-indexed nodes, to the node itself, to higher-indexed nodes,
!    the total number of adjacencies for this node, and the location
!    of the first and last entries required to list this set of adjacencies
!    in a single list of all the adjacencies.
!
!    N   Adjacencies                Below  Self    Above  Total First  Last
!
!   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
!    1:  *  2  6                        0     1       2       3     1     3
!    2:  1  *  3  6  7                  1     1       3       5     4     8
!    3:  2  *  4  7  8                  1     1       3       5     9    13
!    4:  3  *  5  8  9                  1     1       3       5    14    18
!    5:  4  *  9 10                     1     1       2       4    19    22
!    6:  1  2  *  7 11                  2     1       2       5    23    27
!    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
!    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
!    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
!   10:  5  9  * 14 15                  2     1       2       5    49    53
!   11:  6  7  * 12 16                  2     1       2       5    54    58
!   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
!   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
!   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
!   15: 10 14  * 19 20                  2     1       2       5    80    84
!   16: 11 12  * 17 21                  2     1       2       5    85    89
!   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
!   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
!   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
!   20: 15 19  * 24 25                  2     1       2       5   111   115
!   21: 16 17  * 22                     2     1       1       4   116   119
!   22: 17 18 21  * 23                  3     1       1       5   120   124
!   23: 18 19 22  * 24                  3     1       1       5   125   129
!   24: 19 20 23  * 25                  3     1       1       5   130   134
!   25: 20 24  *                        2     1       0       3   135   137
!   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), lists the nodes
!    that make up each triangle in counterclockwise order.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), for each
!    side of a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacencies.
!
!    Input, integer ( kind = 4 ) ADJ_COL(NODE_NUM+1).  Information about
!    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
!
!    Output, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency information.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_col(node_num+1)
  integer ( kind = 4 ) adj_copy(node_num)
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) node
  integer ( kind = 4 ) number
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle2
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)

  adj(1:adj_num) = -1
  adj_copy(1:node_num) = adj_col(1:node_num)
!
!  Set every node to be adjacent to itself.
!
  do node = 1, node_num
    adj(adj_copy(node)) = node
    adj_copy(node) = adj_copy(node) + 1
  end do
!
!  Examine each triangle.
!
  do triangle = 1, element_num

    n1 = element_node(1,triangle)
    n2 = element_node(2,triangle)
    n3 = element_node(3,triangle)
!
!  Add edge (1,2) if this is the first occurrence,
!  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
    triangle2 = element_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n1)) = n2
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n2)) = n1
      adj_copy(n2) = adj_copy(n2) + 1
    end if
!
!  Add edge (2,3).
!
    triangle2 = element_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n2)) = n3
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n3)) = n2
      adj_copy(n3) = adj_copy(n3) + 1
    end if
!
!  Add edge (3,1).
!
    triangle2 = element_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n1)) = n3
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n3)) = n1
      adj_copy(n3) = adj_copy(n3) + 1
    end if

  end do
!
!  Ascending sort the entries for each node.
!
  do node = 1, node_num
    k1 = adj_col(node)
    k2 = adj_col(node+1)-1
    number = k2 + 1 - k1
    call i4vec_sort_heap_a ( number, adj(k1:k2) )
  end do

  return
end
subroutine triangulation_order3_adj_set2 ( node_num, element_num, &
  element_node, element_neighbor, adj_num, adj_col, ia, ja )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_ADJ_SET2 sets adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to set up the arrays IA and JA that
!    record which nodes are adjacent in a triangulation.
!
!    The triangulation is assumed to involve 3-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!    Also, a node is considered to be adjacent to itself.
!
!    This routine can be used to set up the sparse triplet storage
!    for a linear triangle finite element discretization of Poisson's
!    equation in two dimensions.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  |   \  side 2
!       |    \
!    3  |     \
!       |      \
!       1-------2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!    A sample grid
!
!
!    Below, we have a chart that summarizes the adjacency relationships
!    in the sample grid.  On the left, we list the node, and its neighbors,
!    with an asterisk to indicate the adjacency of the node to itself
!    (in some cases, you want to count this self adjacency and in some
!    you don't).  On the right, we list the number of adjacencies to
!    lower-indexed nodes, to the node itself, to higher-indexed nodes,
!    the total number of adjacencies for this node, and the location
!    of the first and last entries required to list this set of adjacencies
!    in a single list of all the adjacencies.
!
!    N   Adjacencies                Below  Self    Above  Total First  Last
!
!   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
!    1:  *  2  6                        0     1       2       3     1     3
!    2:  1  *  3  6  7                  1     1       3       5     4     8
!    3:  2  *  4  7  8                  1     1       3       5     9    13
!    4:  3  *  5  8  9                  1     1       3       5    14    18
!    5:  4  *  9 10                     1     1       2       4    19    22
!    6:  1  2  *  7 11                  2     1       2       5    23    27
!    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
!    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
!    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
!   10:  5  9  * 14 15                  2     1       2       5    49    53
!   11:  6  7  * 12 16                  2     1       2       5    54    58
!   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
!   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
!   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
!   15: 10 14  * 19 20                  2     1       2       5    80    84
!   16: 11 12  * 17 21                  2     1       2       5    85    89
!   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
!   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
!   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
!   20: 15 19  * 24 25                  2     1       2       5   111   115
!   21: 16 17  * 22                     2     1       1       4   116   119
!   22: 17 18 21  * 23                  3     1       1       5   120   124
!   23: 18 19 22  * 24                  3     1       1       5   125   129
!   24: 19 20 23  * 25                  3     1       1       5   130   134
!   25: 20 24  *                        2     1       0       3   135   137
!   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
!
!    For this example, the initial portion of the IA and JA arrays will be:
!
!      (1,1), (1,2), (1,6),
!      (2,1), (2,2), (2,3), (2,6), (2,7),
!      (3,2), (3,3), (3,4), (3,7), (3,8),
!      ...
!      (25,20), (25,24), (25,25)
!
!    for a total of 137 pairs of values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), lists the nodes
!    that make up each triangle in counterclockwise order.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), for each
!    side of a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacencies.
!
!    Input, integer ( kind = 4 ) ADJ_COL(NODE_NUM+1).  Information about
!    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
!
!    Output, integer ( kind = 4 ) IA(ADJ_NUM), JA(ADJ_NUM), the adjacency
!    information.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) adj_col(node_num+1)
  integer ( kind = 4 ) adj_copy(node_num)
  integer ( kind = 4 ) ia(adj_num)
  integer ( kind = 4 ) ja(adj_num)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) node
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle2
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)

  ia(1:adj_num) = -1
  ja(1:adj_num) = -1

  adj_copy(1:node_num) = adj_col(1:node_num)
!
!  Set every node to be adjacent to itself.
!
  do node = 1, node_num
    ia(adj_copy(node)) = node
    ja(adj_copy(node)) = node
    adj_copy(node) = adj_copy(node) + 1
  end do
!
!  Examine each triangle.
!
  do triangle = 1, element_num

    n1 = element_node(1,triangle)
    n2 = element_node(2,triangle)
    n3 = element_node(3,triangle)
!
!  Add edge (1,2) if this is the first occurrence,
!  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
    triangle2 = element_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then

      ia(adj_copy(n1)) = n1
      ja(adj_copy(n1)) = n2
      adj_copy(n1) = adj_copy(n1) + 1

      ia(adj_copy(n2)) = n2
      ja(adj_copy(n2)) = n1
      adj_copy(n2) = adj_copy(n2) + 1

    end if
!
!  Add edge (2,3).
!
    triangle2 = element_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then

      ia(adj_copy(n2)) = n2
      ja(adj_copy(n2)) = n3
      adj_copy(n2) = adj_copy(n2) + 1

      ia(adj_copy(n3)) = n3
      ja(adj_copy(n3)) = n2
      adj_copy(n3) = adj_copy(n3) + 1

    end if
!
!  Add edge (3,1).
!
    triangle2 = element_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then

      ia(adj_copy(n1)) = n1
      ja(adj_copy(n1)) = n3
      adj_copy(n1) = adj_copy(n1) + 1

      ia(adj_copy(n3)) = n3
      ja(adj_copy(n3)) = n1
      adj_copy(n3) = adj_copy(n3) + 1

    end if

  end do
!
!  Lexically sort the IA, JA values.
!
  call i4vec2_sort_a ( adj_num, ia, ja )

  return
end
subroutine triangulation_order3_boundary_edge_count ( element_num, &
  element_node, boundary_edge_num )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT counts the boundary edges.
!
!  Discussion:
!
!    This routine is given a triangulation, an abstract list of triples
!    of nodes.  It is assumed that the nodes in each triangle are listed
!    in a counterclockwise order, although the routine should work
!    if the nodes are consistently listed in a clockwise order as well.
!
!    It is assumed that each edge of the triangulation is either
!    * an INTERIOR edge, which is listed twice, once with positive
!      orientation and once with negative orientation, or;
!    * a BOUNDARY edge, which will occur only once.
!
!    This routine should work even if the region has holes - as long
!    as the boundary of the hole comprises more than 3 edges!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes
!    that make up the triangles.  These should be listed in counterclockwise
!    order.
!
!    Output, integer ( kind = 4 ) BOUNDARY_EDGE_NUM, the number of boundary
!    edges.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) boundary_edge_num
  integer ( kind = 4 ) e1(3*element_num)
  integer ( kind = 4 ) e2(3*element_num)
  integer ( kind = 4 ) edge(2,3*element_num)
  integer ( kind = 4 ) interior_edge_num
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) unique_num

  m = 2
  n = 3 * element_num
!
!  Set up the edge array.
!
  edge(1:2,               1:  element_num) = element_node(1:2,1:element_num)
  edge(1:2,  element_num+1:2*element_num) = element_node(2:3,1:element_num)
  edge(1  ,2*element_num+1:3*element_num) = element_node(3,  1:element_num)
  edge(2  ,2*element_num+1:3*element_num) = element_node(1,  1:element_num)
!
!  In each column, force the smaller entry to appear first.
!
  e1(1:n) = minval ( edge(1:2,1:n), dim = 1 )
  e2(1:n) = maxval ( edge(1:2,1:n), dim = 1 )

  edge(1,1:n) = e1(1:n)
  edge(2,1:n) = e2(1:n)
!
!  Ascending sort the column array.
!
  call i4col_sort_a ( m, n, edge )
!
!  Get the number of unique columns in EDGE.
!
  call i4col_sorted_unique_count ( m, n, edge, unique_num )

  interior_edge_num = 3 * element_num - unique_num

  boundary_edge_num = 3 * element_num - 2 * interior_edge_num

  return
end
subroutine triangulation_order3_boundary_edge_count_euler ( node_num, &
  element_num, hole_num, boundary_num )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER counts boundary edges.
!
!  Discussion:
!
!    We assume we are given information about a triangulation
!    of a set of nodes in the plane.
!
!    Given the number of nodes, triangles and holes, we are going to apply
!    Euler's formula to determine the number of edges that lie on the
!    boundary of the set of nodes.
!
!    The number of faces, including the infinite face and internal holes,
!    is ELEMENT_NUM + HOLE_NUM + 1.
!
!    Let BOUNDARY_NUM denote the number of edges on the boundary.
!    Each of the ELEMENT_NUM triangles uses three edges.  Every edge
!    occurs in two different faces, so the number of edges must be
!    ( 3 * ELEMENT_NUM + BOUNDARY_NUM ) / 2.
!
!    The number of nodes used in the triangulation is NODE_NUM.
!
!    Euler's formula asserts that, for a simple connected figure in the
!    plane with no edge crossings, NODE_NUM nodes, EDGE_NUM edges and
!    FACE_NUM faces:
!
!      NODE_NUM - EDGE_NUM + FACE_NUM = 2
!
!    In our context, this becomes
!
!      NODE_NUM - ( 3 * ELEMENT_NUM + BOUNDARY_NUM ) / 2
!      + ELEMENT_NUM + HOLE_NUM + 1 = 2
!
!    or
!
!      BOUNDARY_NUM = 2 * NODE_NUM + 2 * HOLE_NUM - ELEMENT_NUM - 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc Krevald, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry, Section 9.1,
!    Springer, 2000.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of internal holes.
!
!    Output, integer ( kind = 4 ) BOUNDARY_NUM, the number of edges that
!    lie on the boundary of the triangulation.
!
  implicit none

  integer ( kind = 4 ) boundary_num
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num

  boundary_num = 2 * node_num + 2 * hole_num - element_num - 2

  return
end
subroutine triangulation_order3_boundary_node ( node_num, element_num, &
  element_node, node_boundary )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_BOUNDARY_NODE indicates which nodes are on the boundary.
!
!  Discussion:
!
!    This routine is given a triangulation, an abstract list of triples
!    of nodes.  It is assumed that the nodes in each triangle are listed
!    in a counterclockwise order, although the routine should work
!    if the nodes are consistently listed in a clockwise order as well.
!
!    It is assumed that each edge of the triangulation is either
!    * an INTERIOR edge, which is listed twice, once with positive
!      orientation and once with negative orientation, or;
!    * a BOUNDARY edge, which will occur only once.
!
!    This routine should work even if the region has holes - as long
!    as the boundary of the hole comprises more than 3 edges!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes
!    that make up the triangles.  These should be listed in counterclockwise
!    order.
!
!    Output, logical NODE_BOUNDARY(NODE_NUM), is TRUE if the node
!    is on a boundary edge.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) e1(3*element_num)
  integer ( kind = 4 ) e2(3*element_num)
  integer ( kind = 4 ) edge(2,3*element_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  logical              node_boundary(node_num)
  integer ( kind = 4 ) element_node(element_order,element_num)

  m = 2
  n = 3 * element_num
!
!  Set up the edge array.
!
  edge(1:2,               1:  element_num) = element_node(1:2,1:element_num)
  edge(1:2,  element_num+1:2*element_num) = element_node(2:3,1:element_num)
  edge(1,  2*element_num+1:3*element_num) = element_node(3,  1:element_num)
  edge(2,  2*element_num+1:3*element_num) = element_node(1,  1:element_num)
!
!  In each column, force the smaller entry to appear first.
!
  e1(1:n) = minval ( edge(1:2,1:n), dim = 1 )
  e2(1:n) = maxval ( edge(1:2,1:n), dim = 1 )

  edge(1,1:n) = e1(1:n)
  edge(2,1:n) = e2(1:n)
!
!  Ascending sort the column array.
!
  call i4col_sort_a ( m, n, edge )
!
!  Records which appear twice are internal edges and can be ignored.
!
  node_boundary(1:node_num) = .false.

  j = 0

  do while ( j < 3 * element_num )

    j = j + 1

    if ( j == 3 * element_num ) then
      node_boundary(edge(1:m,j)) = .true.
    else if ( all ( edge(1:m,j) == edge(1:m,j+1) ) ) then
      j = j + 1
    else
      node_boundary(edge(1:m,j)) = .true.
    end if

  end do

  return
end
subroutine triangulation_order3_check ( node_num, element_num, &
  element_node, ierror )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_CHECK makes some simple checks on a triangulation.
!
!  Discussion:
!
!    Because this routine does not receive the physical coordinates of
!    the nodes, it cannot ensure that the triangulation is maximal,
!    that is, that no more triangles can be created.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes
!    that make up the triangles.  These should be listed in counterclockwise
!    order.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred, the triangulation is not valid.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) boundary_num
  integer ( kind = 4 ) euler
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) used(node_num)

  ierror = 0
!
!  Checks 1 and 2:
!  NODE_NUM must be at least 3.
!  ELEMENT_NUM must be at least 1.
!
  if ( node_num < 3 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
    write ( *, '(a)' ) '  The number of nodes is less than 3!'
    return
  end if

  if ( element_num < 1 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
    write ( *, '(a)' ) '  The number of triangles is less than 1!'
    return
  end if
!
!  Checks 3 and 4:
!  Verify that all node values are greater than or equal to 1
!  and less than or equal to NODE_NUM.
!
  if ( any ( element_node(1:3,1:element_num) < 1 ) ) then
    ierror = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
    write ( *, '(a)' ) '  Some nodes are less than 1!'
    return
  end if

  if ( any ( node_num < element_node(1:3,1:element_num) ) ) then
    ierror = 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
    write ( *, '(a)' ) '  Some nodes are greater than NODE_NUM!'
    return
  end if
!
!  Check 5:
!  Verify that every node is used at least once.
!
  used(1:node_num) = 0

  used(element_node(1,1:element_num)) = &
    used(element_node(1,1:element_num)) + 1
  used(element_node(2,1:element_num)) = &
    used(element_node(2,1:element_num)) + 1
  used(element_node(3,1:element_num)) = &
    used(element_node(3,1:element_num)) + 1

  if ( any ( used(1:node_num) == 0 ) ) then
    ierror = 5
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
    write ( *, '(a)' ) '  Some nodes are never used as triangle vertices!'
    return
  end if
!
!  Check 6:
!  Verify that no node is repeated in a triangle.
!
  do i = 1, element_num
    if ( element_node(1,i) == element_node(2,i) .or. &
         element_node(2,i) == element_node(3,i) .or. &
         element_node(3,i) == element_node(1,i) ) then
      ierror = 6
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
      write ( *, '(a)' ) '  A triangle contains a null edge!'
      return
    end if
  end do
!
!  Check 7:
!  Verify that no edge is repeated, and that repeated edges occur in
!  negated pairs.
!
  call triangulation_order3_edge_check ( element_num, element_node, &
    boundary_num, ierror )

  if ( ierror /= 0 ) then
    ierror = 7
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
    write ( *, '(a)' ) '  Some edges are repeated,'
    write ( *, '(a)' ) '  or given in the wrong direction!'
    return
  end if
!
!  Check 8:
!  Does the triangulation satisfy Euler's criterion?
!  If not, then the triangulation is not proper.  (For instance, there
!  might be a hole in the interior.)
!
  euler = boundary_num + element_num + 2 - 2 * node_num

  if ( euler /= 0 ) then
    ierror = 8
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
    write ( *, '(a)' ) '  The triangulation fails Euler''s criterion!'
    return
  end if

  return
end
subroutine triangulation_order3_edge_check ( element_num, element_node, &
  boundary_num, ierror )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_EDGE_CHECK checks the edges of a triangulation.
!
!  Discussion:
!
!    This routine was revised to store the edge data as columns,
!    rather than rows.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes
!    that make up each triangle.
!
!    Output, integer ( kind = 4 ) BOUNDARY_NUM, the number of edges that
!    lie on the boundary.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no errors were detected.
!    nonzero, an error occurred.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) boundary_num
  integer ( kind = 4 ) col(3,3*element_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) tri
  integer ( kind = 4 ) element_node(element_order,element_num)

  ierror = 0
!
!  Step 1.
!  From the list of nodes for triangle T, of the form: (I,J,K)
!  construct the three neighbor relations:
!
!    (I,J,+1) or (J,I,-1),
!    (J,K,+1) or (K,J,-1),
!    (K,I,+1) or (I,K,-1)
!
!  where we choose (I,J,+1) if I < J, or else (J,I,-1) and so on.
!
  do tri = 1, element_num

    i = element_node(1,tri)
    j = element_node(2,tri)
    k = element_node(3,tri)

    if ( i < j ) then
      col(1:3,3*(tri-1)+1) = (/ i, j, +1 /)
    else
      col(1:3,3*(tri-1)+1) = (/ j, i, -1 /)
    end if

    if ( j < k ) then
      col(1:3,3*(tri-1)+2) = (/ j, k, +1 /)
    else
      col(1:3,3*(tri-1)+2) = (/ k, j, -1 /)
    end if

    if ( k < i ) then
      col(1:3,3*(tri-1)+3) = (/ k, i, +1 /)
    else
      col(1:3,3*(tri-1)+3) = (/ i, k, -1 /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!
  call i4col_sort_a ( 3, 3*element_num, col )
!
!  Step 3.
!
!  Most records (/ A, B, C /) occur twice, with C being -1, then +1.
!  These records represent internal edges.
!
!  Unpaired records represent edges on the boundary.
!
  i = 0
  boundary_num = 0

  do while ( i < 3 * element_num )

    i = i + 1

    if ( i == 3 * element_num ) then

      boundary_num = boundary_num + 1

    else

      if ( col(1,i) == col(1,i+1) .and. col(2,i) == col(2,i+1) ) then

        if ( col(3,i) == col(3,i+1) ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TRIANGULATION_ORDER3_EDGE_CHECK - Warning!'
          write ( *, '(a)' ) '  An edge occurs twice.'
          return
        end if

        i = i + 1

      else

        boundary_num = boundary_num + 1

      end if

    end if

  end do

  return
end
subroutine triangulation_order3_example1 ( node_num, element_num, node_xy, &
  element_node, element_neighbor )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_EXAMPLE1 sets up a sample triangulation.
!
!  Discussion:
!
!    This triangulation is actually a Delaunay triangulation.
!
!    The appropriate input values of NODE_NUM and ELEMENT_NUM can be
!    determined by calling TRIANGULATION_ORDER3_EXAMPLE1_SIZE first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Output, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the
!    nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes
!    that make up the triangles.
!
!    Output, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbors on each side.  Negative values indicate edges that
!    lie on the exterior.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) element_neighbor(3,element_num)

  node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       2.0D+00, 2.0D+00, &
      -1.0D+00, 3.0D+00, &
      -2.0D+00, 2.0D+00, &
       8.0D+00, 2.0D+00, &
       9.0D+00, 5.0D+00, &
       7.0D+00, 4.0D+00, &
       5.0D+00, 6.0D+00, &
       6.0D+00, 7.0D+00, &
       8.0D+00, 8.0D+00, &
      11.0D+00, 7.0D+00, &
      10.0D+00, 4.0D+00, &
       6.0D+00, 4.0D+00 /), (/ dim_num, node_num /) )

  element_node(1:element_order,1:element_num ) = reshape ( (/ &
     3,   4,   1, &
     3,   1,   2, &
     3,   2,   8, &
     2,   1,   5, &
     8,   2,  13, &
     8,  13,   9, &
     3,   8,   9, &
    13,   2,   5, &
     9,  13,   7, &
     7,  13,   5, &
     6,   7,   5, &
     9,   7,   6, &
    10,   9,   6, &
     6,   5,  12, &
    11,   6,  12, &
    10,   6,  11 /), (/ element_order, element_num /) )

  element_neighbor(1:3,1:element_num) = reshape ( (/ &
       -4,  -13,    2, &
        1,    4,    3, &
        2,    5,    7, &
        2,  -43,    8, &
        3,    8,    6, &
        5,    9,    7, &
        3,    6,   -3, &
        5,    4,   10, &
        6,   10,   12, &
        9,    8,   11, &
       12,   10,   14, &
        9,   11,   13, &
      -23,   12,   16, &
       11,  -47,   15, &
       16,   14,  -50, &
       13,   15,  -39 /), (/ 3, element_num /) )

  return
end
subroutine triangulation_order3_example1_size ( node_num, element_num, &
  hole_num )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_EXAMPLE1_SIZE sets sizes for a sample triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num

  hole_num = 0
  node_num = 13
  element_num = 16

  return
end
subroutine triangulation_order3_example2 ( node_num, element_num, node_xy, &
  element_node, element_neighbor )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_EXAMPLE2 returns an example triangulation.
!
!  Discussion:
!
!    This triangulation is actually a Delaunay triangulation.
!
!    The appropriate input values of NODE_NUM and ELEMENT_NUM can be
!    determined by calling TRIANGULATION_ORDER3_EXAMPLE2_SIZE first.
!
!  Diagram:
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Output, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the
!    nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), lists the
!    nodes that make up each triangle, in counterclockwise order.
!
!    Output, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), for
!    each side of a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)

  node_xy = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    2.0D+00, 0.0D+00, &
    3.0D+00, 0.0D+00, &
    4.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, &
    2.0D+00, 1.0D+00, &
    3.0D+00, 1.0D+00, &
    4.0D+00, 1.0D+00, &
    0.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00, &
    2.0D+00, 2.0D+00, &
    3.0D+00, 2.0D+00, &
    4.0D+00, 2.0D+00, &
    0.0D+00, 3.0D+00, &
    1.0D+00, 3.0D+00, &
    2.0D+00, 3.0D+00, &
    3.0D+00, 3.0D+00, &
    4.0D+00, 3.0D+00, &
    0.0D+00, 4.0D+00, &
    1.0D+00, 4.0D+00, &
    2.0D+00, 4.0D+00, &
    3.0D+00, 4.0D+00, &
    4.0D+00, 4.0D+00  &
  /), (/ dim_num, node_num /) )

  element_node(1:element_order,1:element_num) = reshape ( (/ &
     1,  2,  6, &
     7,  6,  2, &
     2,  3,  7, &
     8,  7,  3, &
     3,  4,  8, &
     9,  8,  4, &
     4,  5,  9, &
    10,  9,  5, &
     6,  7, 11, &
    12, 11,  7, &
     7,  8, 12, &
    13, 12,  8, &
     8,  9, 13, &
    14, 13,  9, &
     9, 10, 14, &
    15, 14, 10, &
    11, 12, 16, &
    17, 16, 12, &
    12, 13, 17, &
    18, 17, 13, &
    13, 14, 18, &
    19, 18, 14, &
    14, 15, 19, &
    20, 19, 15, &
    16, 17, 21, &
    22, 21, 17, &
    17, 18, 22, &
    23, 22, 18, &
    18, 19, 23, &
    24, 23, 19, &
    19, 20, 24, &
    25, 24, 20 /), (/ element_order, element_num /) )

  element_neighbor(1:3,1:element_num) = reshape ( (/ &
    -1,  2, -1, &
     9,  1,  3, &
    -1,  4,  2, &
    11,  3,  5, &
    -1,  6,  4, &
    13,  5,  7, &
    -1,  8,  6, &
    15,  7, -1, &
     2, 10, -1, &
    17,  9, 11, &
     4, 12, 10, &
    19, 11, 13, &
     6, 14, 12, &
    21, 13, 15, &
     8, 16, 14, &
    23, 15, -1, &
    10, 18, -1, &
    25, 17, 19, &
    12, 20, 18, &
    27, 19, 21, &
    14, 22, 20, &
    29, 21, 23, &
    16, 24, 22, &
    31, 23, -1, &
    18, 26, -1, &
    -1, 25, 27, &
    20, 28, 26, &
    -1, 27, 29, &
    22, 30, 28, &
    -1, 29, 31, &
    24, 32, 30, &
    -1, 31, -1 /), (/ 3, element_num /) )

  return
end
subroutine triangulation_order3_example2_size ( node_num, element_num, &
  hole_num )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_EXAMPLE2_SIZE returns the size of an example.
!
!  Diagram:
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num

  node_num = 25
  element_num = 32
  hole_num = 0

  return
end
subroutine triangulation_order3_neighbor ( element_num, element_node, &
  t1, s1, t2, s2 )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_NEIGHBOR determines a neighbor of a given triangle.
!
!  Discussion:
!
!    A set of nodes is given.  A triangulation of the nodes has been
!    defined and recorded in TRIANGLE.  The TRIANGLE data structure records
!    triangles as sets of three nodes, N1, N2, N3, that implicitly define three
!    sides, being the line segments N1-N2, N2-N3, and N3-N1.
!
!    The nodes of the triangle are listed in counterclockwise order.
!    This means that if two triangles share a side, then the nodes
!    defining that side occur in the order (N1,N2) for one triangle,
!    and (N2,N1) for the other.
!
!    The routine is given a triangle and a side, and asked to find
!    another triangle (if any) that shares that side.  The routine
!    simply searches the ELEMENT_NODE structure for an occurrence of the
!    nodes in the opposite order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes
!    that define each triangle.
!
!    Input, integer ( kind = 4 ) T1, the index of the triangle.
!
!    Input, integer ( kind = 4 ) S1, the index of the triangle side.
!
!    Output, integer ( kind = 4 ) T2, the index of the triangle which is
!    the neighbor to T1 on side S1, or -1 if there is no such neighbor.
!
!    Output, integer ( kind = 4 ) S2, the index of the side of triangle T2
!    which is shared with triangle T1, or -1 if there is no such neighbor.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) s
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) ss
  integer ( kind = 4 ) t
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2
  integer ( kind = 4 ) element_node(element_order,element_num)

  t2 = -1
  s2 = -1

  n1 = element_node(s1,t1)
  ss = s1 + 1
  ss = i4_wrap ( ss, 1, 3 )
  n2 = element_node(ss,t1)

  do t = 1, element_num
    do s = 1, 3
      if ( element_node(s,t) == n1 ) then
        ss = s - 1
        ss = i4_wrap ( ss, 1, 3 )
        if ( element_node(ss,t) == n2 ) then
          t2 = t
          s2 = ss
          return
        end if
      end if
    end do
  end do

  return
end
subroutine triangulation_order3_neighbor_nodes ( node_num, element_num, &
  nabes_max, element_node, nabes_first, nabes_num, nabes_dim, nabes )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_NEIGHBOR_NODES determines triangulation neighbor nodes.
!
!  Example:
!
!    On input, the triangle data structure is:
!
!    Triangle  Nodes
!    --------  ----------
!     1        3,   4,   1
!     2        3,   1,   2
!     3        3,   2,   6
!     4        2,   1,   5
!     5        6,   2,   5
!
!  On output, the auxilliary neighbor arrays are:
!
!    Node  Num  First
!    ----  ---  -----
!     1     4     1
!     2     4     5
!     3     4     9
!     4     2    13
!     5     3    15
!     6     3    18
!
!  and the neighbor array is:
!
!    Position  Node
!    --------  ----
!
!     1        2
!     2        3
!     3        4
!     4        5
!    -----------
!     5        1
!     6        3
!     7        5
!     8        6
!    -----------
!     9        1
!    10        2
!    11        4
!    12        6
!    -----------
!    13        1
!    14        3
!    -----------
!    15        1
!    16        2
!    17        6
!    -----------
!    18        2
!    19        3
!    20        5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) NABES_MAX, the maximum dimension of NABES.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes
!    that make up each triangle.
!
!    Output, integer ( kind = 4 ) NABES_FIRST(NODE_NUM), the index in NABES
!    of the first neighbor in the list for each node.
!
!    Output, integer ( kind = 4 ) NABES_NUM(NODE_NUM), the number of neighbors
!    of each node.
!
!    Output, integer ( kind = 4 ) NABES_DIM, the dimension of NABES.
!
!    Output, integer ( kind = 4 ) NABES(NABES_DIM), a list of the neighbors
!    of all the nodes.  Neighbors of node 1 are listed first, and so on.
!
  implicit none

  integer ( kind = 4 ) nabes_max
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_current
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nabe
  integer ( kind = 4 ) nabes(nabes_max)
  integer ( kind = 4 ) nabes1(nabes_max)
  integer ( kind = 4 ) nabes_dim
  integer ( kind = 4 ) nabes_first(node_num)
  integer ( kind = 4 ) nabes_num(node_num)
  integer ( kind = 4 ) tri
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) unique_num
!
!  Step 1.  From the triangle list (I,J,K)
!  construct the neighbor relations: (I,J), (J,K), (K,I), (J,I), (K,J), (I,K).
!
  nabes_dim = 0
  do tri = 1, element_num
    i = element_node(1,tri)
    j = element_node(2,tri)
    k = element_node(3,tri)
    nabes1(nabes_dim+1:nabes_dim+6) = (/ i, i, j, j, k, k /)
    nabes(nabes_dim+1:nabes_dim+6) = (/ j, k, i, k, i, j /)
    nabes_dim = nabes_dim + 6
  end do
!
!  Step 2. Dictionary sort the neighbor relations.
!
  call i4vec2_sort_a ( nabes_dim, nabes1, nabes )
!
!  Step 3. Remove duplicate entries.
!
  call i4vec2_sorted_unique ( nabes_dim, nabes1, nabes, unique_num )

  nabes_dim = unique_num
!
!  Step 4. Construct the NABES_NUM and NABES_FIRST data.
!
  nabes_num(1:node_num) = 0
  nabes_first(1:node_num) = 0
  i_current = 0
  do nabe = 1, nabes_dim
    i = nabes1(nabe)
    if ( i == i_current ) then
      nabes_num(i) = nabes_num(i) + 1
    else
      i_current = i
      nabes_first(i) = nabe
      nabes_num(i) = 1
    end if
  end do

  return
end
subroutine triangulation_order3_neighbor_nodes_print ( node_num, nabes_first, &
  nabes_num, nabes_dim, nabes )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_NEIGHBOR_NODES_PRINT prints a node neighbor array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) NABES_FIRST(NODE_NUM), the index in NABES
!    of the first neighbor in the list for each node.
!
!    Input, integer ( kind = 4 ) NABES_NUM(NODE_NUM), the number of neighbors
!    of each node.
!
!    Input, integer ( kind = 4 ) NABES_DIM, the dimension of NABES.
!
!    Input, integer ( kind = 4 ) NABES(NABES_DIM), a list of the neighbors
!    of all the nodes.  Neighbors of node 1 are listed first, and so on.
!
  implicit none

  integer ( kind = 4 ) nabes_dim
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) nabes(nabes_dim)
  integer ( kind = 4 ) nabes_first(node_num)
  integer ( kind = 4 ) nabes_num(node_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node based arrays:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node  Neighbors  Index #1'
  write ( *, '(a)' ) ' '
  do i = 1, node_num
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, nabes_num(i), nabes_first(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The raw neighbor array.'
  write ( *, '(a)' ) ' '
  do i = 1, nabes_dim
    write ( *, '(2x,i8,2x,i8)' ) i, nabes(i)
  end do

  return
end
subroutine triangulation_order3_plot ( file_name, node_num, node_xy, &
  element_num, element_node, node_show, element_show )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_PLOT plots a 3-node triangulation of a set of nodes.
!
!  Discussion:
!
!    The triangulation is most usually a Delaunay triangulation,
!    but this is not necessary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), lists, for
!    each triangle, the indices of the nodes that form the vertices of the
!    triangle.
!
!    Input, integer ( kind = 4 ) NODE_SHOW,
!    0, do not show nodes;
!    1, show nodes;
!    2, show nodes and label them.
!
!    Input, integer ( kind = 4 ) ELEMENT_SHOW,
!    0, do not show triangles;
!    1, show triangles;
!    2, show triangles and label them.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  real ( kind = 8 ) ave_x
  real ( kind = 8 ) ave_y
  integer ( kind = 4 ) :: circle_size
  integer ( kind = 4 ) delta
  integer ( kind = 4 ) e
  character ( len = * )  file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_show
  real ( kind = 8 ) node_xy(2,node_num)
  character ( len = 40 ) string
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) element_show
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer ( kind = 4 ) x_ps
  integer ( kind = 4 ) :: x_ps_max = 576
  integer ( kind = 4 ) :: x_ps_max_clip = 594
  integer ( kind = 4 ) :: x_ps_min = 36
  integer ( kind = 4 ) :: x_ps_min_clip = 18
  real ( kind = 8 ) x_scale
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  integer ( kind = 4 ) y_ps
  integer ( kind = 4 ) :: y_ps_max = 666
  integer ( kind = 4 ) :: y_ps_max_clip = 684
  integer ( kind = 4 ) :: y_ps_min = 126
  integer ( kind = 4 ) :: y_ps_min_clip = 108
  real ( kind = 8 ) y_scale
!
!  We need to do some figuring here, so that we can determine
!  the range of the data, and hence the height and width
!  of the piece of paper.
!
  x_max = maxval ( node_xy(1,1:node_num) )
  x_min = minval ( node_xy(1,1:node_num) )
  x_scale = x_max - x_min

  x_max = x_max + 0.05D+00 * x_scale
  x_min = x_min - 0.05D+00 * x_scale
  x_scale = x_max - x_min

  y_max = maxval ( node_xy(2,1:node_num) )
  y_min = minval ( node_xy(2,1:node_num) )
  y_scale = y_max - y_min

  y_max = y_max + 0.05D+00 * y_scale
  y_min = y_min - 0.05D+00 * y_scale
  y_scale = y_max - y_min

  if ( x_scale < y_scale ) then

    delta = nint ( real ( x_ps_max - x_ps_min, kind = 8 ) &
      * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )

    x_ps_max = x_ps_max - delta
    x_ps_min = x_ps_min + delta

    x_ps_max_clip = x_ps_max_clip - delta
    x_ps_min_clip = x_ps_min_clip + delta

    x_scale = y_scale

  else if ( y_scale < x_scale ) then

    delta = nint ( real ( y_ps_max - y_ps_min, kind = 8 ) &
      * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )

    y_ps_max      = y_ps_max - delta
    y_ps_min      = y_ps_min + delta

    y_ps_max_clip = y_ps_max_clip - delta
    y_ps_min_clip = y_ps_min_clip + delta

    y_scale = x_scale

  end if

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_ORDER3_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Can not open output file.'
    return
  end if

  write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( file_unit, '(a)' ) '%%Creator: triangulation_order3_plot.f90'
  write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( file_unit, '(a)' ) '%%Pages: 1'
  write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%BoundingBox: ', &
    x_ps_min, y_ps_min, x_ps_max, y_ps_max
  write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( file_unit, '(a)' ) '%%EndComments'
  write ( file_unit, '(a)' ) '%%BeginProlog'
  write ( file_unit, '(a)' ) '/inch {72 mul} def'
  write ( file_unit, '(a)' ) '%%EndProlog'
  write ( file_unit, '(a)' ) '%%Page: 1 1'
  write ( file_unit, '(a)' ) 'save'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Increase line width from default 0.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '2 setlinewidth'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Draw a gray border around the page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_min, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' lineto'
  write ( file_unit, '(a)' ) 'stroke'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB color to black.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the font and its size.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '/Times-Roman findfont'
  write ( file_unit, '(a)' ) '0.50 inch scalefont'
  write ( file_unit, '(a)' ) 'setfont'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Print a title.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  210  702  moveto'
  write ( file_unit, '(a)' ) '%  (Triangulation)  show'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a)' ) 'clip newpath'
!
!  Draw the nodes.
!
  if ( node_num <= 200 ) then
    circle_size = 5
  else if ( node_num <= 500 ) then
    circle_size = 4
  else if ( node_num <= 1000 ) then
    circle_size = 3
  else if ( node_num <= 5000 ) then
    circle_size = 2
  else
    circle_size = 1
  end if

  if ( 1 <= node_show ) then
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw filled dots at the nodes.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to blue.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
    write ( file_unit, '(a)' ) '%'

    do node = 1, node_num

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
        + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                   - y_min ) )

      write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
        circle_size, '0 360 arc closepath fill'

    end do

  end if
!
!  Label the nodes.
!
  if ( 2 <= node_show ) then

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Label the nodes:'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to darker blue.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.000  0.250  0.850 setrgbcolor'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.20 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'
    write ( file_unit, '(a)' ) '%'

    do node = 1, node_num

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
        + (       + node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                   - y_min ) )

      write ( string, '(i4)' ) node
      string = adjustl ( string )

      write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, &
        ' moveto (' // trim ( string ) // ') show'

    end do

  end if
!
!  Draw the triangles.
!
  if ( 1 <= element_show ) then
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to red.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw the triangles.'
    write ( file_unit, '(a)' ) '%'

    do triangle = 1, element_num

      write ( file_unit, '(a)' ) 'newpath'

      do i = 1, 4

        e = i4_wrap ( i, 1, 3 )

        node = element_node(e,triangle)

        x_ps = int ( &
          ( ( x_max - node_xy(1,node)         ) &
          * real ( x_ps_min, kind = 8 )   &
          + (         node_xy(1,node) - x_min ) &
          * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                   - x_min ) )

        y_ps = int ( &
          ( ( y_max - node_xy(2,node)         ) &
          * real ( y_ps_min, kind = 8 )   &
          + (         node_xy(2,node) - y_min ) &
          * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                   - y_min ) )

        if ( i == 1 ) then
          write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' moveto'
        else
          write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' lineto'
        end if

      end do

      write ( file_unit, '(a)' ) 'stroke'

    end do

  end if
!
!  Label the triangles.
!
  if ( 2 <= element_show ) then

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Label the triangles:'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to darker red.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.950  0.250  0.150 setrgbcolor'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.20 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'
    write ( file_unit, '(a)' ) '%'

    do triangle = 1, element_num

      ave_x = 0.0D+00
      ave_y = 0.0D+00

      do i = 1, 3

        node = element_node(i,triangle)

        ave_x = ave_x + node_xy(1,node)
        ave_y = ave_y + node_xy(2,node)

      end do

      ave_x = ave_x / 3.0D+00
      ave_y = ave_y / 3.0D+00

      x_ps = int ( &
        ( ( x_max - ave_x         ) * real ( x_ps_min, kind = 8 )   &
        + (       + ave_x - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max         - x_min ) )

      y_ps = int ( &
        ( ( y_max - ave_y         ) * real ( y_ps_min, kind = 8 )   &
        + (         ave_y - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max         - y_min ) )

      write ( string, '(i4)' ) triangle
      string = adjustl ( string )

      write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps, ' moveto (' &
        // trim ( string ) // ') show'

    end do

  end if

  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'restore  showpage'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  End of page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%%Trailer'
  write ( file_unit, '(a)' ) '%%EOF'
  close ( unit = file_unit )

  return
end
subroutine triangulation_order3_print ( node_num, element_num, node_xy, &
  element_node, element_neighbor )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_PRINT prints information about a triangulation.
!
!  Discussion:
!
!    Triangulations created by R8TRIS2 include extra information encoded
!    in the negative values of ELEMENT_NEIGHBOR.
!
!    Because some of the nodes counted in NODE_NUM may not actually be
!    used in the triangulation, I needed to compute the true number
!    of vertices.  I added this calculation on 13 October 2001.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes
!    that make up the triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbors on each side.  If there is no triangle neighbor on
!    a particular side, the value of ELEMENT_NEIGHBOR should be negative.
!    If the triangulation data was created by R8TRIS22, then there is more
!    information encoded in the negative values.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) boundary_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) s
  logical skip
  integer ( kind = 4 ) sp1
  integer ( kind = 4 ) t
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: vertex_list
  integer ( kind = 4 ) vertex_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGULATION_ORDER3_PRINT'
  write ( *, '(a)' ) '  Information defining an order3 triangulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of nodes is ', node_num

  call r8mat_transpose_print ( dim_num, node_num, node_xy, &
    '  Node coordinates' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of triangles is ', element_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sets of three nodes are used as vertices of'
  write ( *, '(a)' ) '  the triangles.  For each triangle, the nodes'
  write ( *, '(a)' ) '  are listed in counterclockwise order.'

  call i4mat_transpose_print ( 3, element_num, element_node, &
    '  Triangle nodes:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  On each side of a given triangle, there is either'
  write ( *, '(a)' ) '  another triangle, or a piece of the convex hull.'
  write ( *, '(a)' ) '  For each triangle, we list the indices of the three'
  write ( *, '(a)' ) '  neighbors, or (if negative) the codes of the'
  write ( *, '(a)' ) '  segments of the convex hull.'

  call i4mat_transpose_print ( 3, element_num, element_neighbor, &
    '  Triangle neighbors' )
!
!  Determine the number of vertices.
!
  allocate ( vertex_list(1:3*element_num) )

  vertex_list(1:3*element_num) = reshape ( element_node(1:3,1:element_num), &
    (/ 3 * element_num /) )

  call i4vec_sort_heap_a ( 3*element_num, vertex_list )

  call i4vec_sorted_unique ( 3*element_num, vertex_list, vertex_num )

  deallocate ( vertex_list )
!
!  Determine the number of boundary points.
!
  boundary_num = 2 * vertex_num - element_num - 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of boundary points is ', boundary_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The segments that make up the convex hull can be'
  write ( *, '(a)' ) '  determined from the negative entries of the triangle'
  write ( *, '(a)' ) '  neighbor list.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     #   Tri  Side    N1    N2'
  write ( *, '(a)' ) ' '

  skip = .false.

  k = 0

  do i = 1, element_num

    do j = 1, 3

      if ( element_neighbor(j,i) < 0 ) then
        s = - element_neighbor(j,i)
        t = s / 3

        if ( t < 1 .or. element_num < t ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Sorry, this data does not use the R8TRIS2'
          write ( *, '(a)' ) '  convention for convex hull segments.'
          skip = .true.
          exit
        end if

        s = mod ( s, 3 ) + 1
        k = k + 1
        n1 = element_node(s,t)
        sp1 = s + 1
        sp1 = i4_wrap ( sp1, 1, 3 )
        n2 = element_node(sp1,t)
        write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4,2x,i4)' ) k, t, s, n1, n2
      end if

    end do

    if ( skip ) then
      exit
    end if

  end do

  return
end
subroutine triangulation_order3_quad ( node_num, node_xy, element_order, &
  element_num, element_node, quad_fun, quad_num, quad_xy, quad_w, &
  quad_value, region_area )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_QUAD approximates an integral over a triangulation.
!
!  Discussion:
!
!    The routine will accept triangulations of order higher than 3.
!    However, only the first three nodes (the vertices) of each
!    triangle will be used.  This will still produce correct results
!    for higher order triangulations, as long as the sides of the
!    triangle are straight.
!
!    We assume that the vertices of each triangle are listed first
!    in the description of higher order triangles, and we assume that
!    the vertices are listed in counterclockwise order.
!
!    The approximation of the integral is made using a quadrature rule
!    defined on the unit triangle, and supplied by the user.
!
!    The user also supplies the name of a subroutine, here called "QUAD_FUN",
!    which evaluates the integrand at a set of points.  The form of
!    this routine is:
!
!      subroutine quad_fun ( n, xy_vec, f_vec )
!      integer n
!      real ( kind = 8 ) f_vec(n)
!      real ( kind = 8 ) xy_vec(2,n)
!
!    and it returns in each entry F_VEC(1:N), the value of the integrand
!    at XY_VEC(1:2,1:N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes in the
!    triangulation.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of triangles in
!    the triangulation.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles in
!    the triangulation.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the nodes making up each triangle.
!
!    Input, external QUAD_FUN, the name of the integrand routine.
!
!    Input, integer ( kind = 4 ) QUAD_NUM, the order of the quadrature rule.
!
!    Input, real ( kind = 8 ) QUAD_XY(2,QUAD_NUM), the abscissas of the
!    quadrature rule, in the unit triangle.
!
!    Input, real ( kind = 8 ) QUAD_W(QUAD_NUM), the weights of the
!    quadrature rule.
!
!    Output, real ( kind = 8 ) QUAD_VALUE, the estimate of the integral
!    of F(X,Y) over the region covered by the triangulation.
!
!    Output, real ( kind = 8 ) REGION_AREA, the area of the region.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) quad_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) quad_f(quad_num)
  external quad_fun
  real ( kind = 8 ) quad_value
  real ( kind = 8 ) quad_w(quad_num)
  real ( kind = 8 ) quad_xy(2,quad_num)
  real ( kind = 8 ) quad2_xy(2,quad_num)
  real ( kind = 8 ) region_area
  integer ( kind = 4 ) triangle
  real ( kind = 8 ) triangle_area
  integer ( kind = 4 ) element_node(element_order,element_num)
  real ( kind = 8 ) triangle_xy(2,3)

  quad_value = 0.0D+00
  region_area = 0.0D+00

  do triangle = 1, element_num

    triangle_xy(1:2,1:3) = node_xy(1:2,element_node(1:3,triangle))

    call triangle_area_2d ( triangle_xy, triangle_area )

    call element_order3_reference_to_physical ( triangle_xy, &
      quad_num, quad_xy, quad2_xy )

    call quad_fun ( quad_num, quad2_xy, quad_f )

    quad_value = quad_value + triangle_area &
      * dot_product ( quad_w(1:quad_num), quad_f(1:quad_num) )

    region_area = region_area + triangle_area

  end do

  return
end
subroutine triangulation_order3_refine_compute ( node_num1, element_num1, &
  node_xy1, element_node1, node_num2, element_num2, edge_data, node_xy2, &
  element_node2 )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_REFINE_COMPUTE computes a refined order 3 triangulation.
!
!  Discussion:
!
!    Given a triangle defined by nodes 1, 2, 3, we need to generate
!    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
!    and T4.
!
!    The task is more complicated by the fact that we are working with
!    a mesh of triangles, so that we want to create a node only once,
!    even though it may be shared by other triangles.
!
!          3
!         / \
!        /T3 \
!      13----23
!      / \T4 / \
!     /T1 \ /T2 \
!    1----12-----2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM1, the number of triangles.
!
!    Input, real ( kind = 8 ) NODE_XY1(2,NODE_NUM1), the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE1(3,ELEMENT_NUM1), the nodes
!    that make up the triangles.  These should be listed in counterclockwise
!    order.
!
!    Input, integer ( kind = 4 ) NODE_NUM2, the number of nodes in the
!    refined mesh.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM2, the number of triangles in
!    the refined mesh.
!
!    Input, integer ( kind = 4 ) EDGE_DATA(5,3*ELEMENT_NUM1), edge information
!    computed by TRIANGULATION_ORDER3_REFINE_SIZE.
!
!    Output, real ( kind = 8 ) NODE_XY2(2,NODE_NUM2), the refined nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE2(3,ELEMENT_NUM2), the nodes
!    that make up the triangles in the refined mesh.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) element_num1
  integer ( kind = 4 ) element_num2

  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(5,3*element_num1)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_old
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_old
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xy1(2,node_num1)
  real ( kind = 8 ) node_xy2(2,node_num2)
  integer ( kind = 4 ) element_node1(3,element_num1)
  integer ( kind = 4 ) element_node2(3,element_num2)
  integer ( kind = 4 ) triangle1
  integer ( kind = 4 ) v1
  integer ( kind = 4 ) v2
!
!  Copy the old nodes.
!
  node_xy2(1:2,1:node_num1) = node_xy1(1:2,1:node_num1)

  element_node2(1:3,1:element_num2) = -1
!
!  We can assign the existing nodes to the new triangles.
!
  do triangle1 = 1, element_num1
    element_node2(1,(triangle1-1)*4+1) = element_node1(1,triangle1)
    element_node2(2,(triangle1-1)*4+2) = element_node1(2,triangle1)
    element_node2(3,(triangle1-1)*4+3) = element_node1(3,triangle1)
  end do

  node = node_num1

  n1_old = -1
  n2_old = -1

  do edge = 1, 3 * element_num1

    n1 = edge_data(1,edge)
    n2 = edge_data(2,edge)
!
!  If this edge is new, create the coordinates and index for this node.
!
    if ( n1 /= n1_old .or. n2 /= n2_old ) then
      node = node + 1

      if ( node_num2 < node ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGLE_MESH_ORDER3_REFINE - Fatal error!'
        write ( *, '(a)' ) '  Node index exceeds NODE_NUM2.'
        stop
      end if

      node_xy2(1:2,node) = &
        ( node_xy2(1:2,n1) + node_xy2(1:2,n2) ) / 2.0D+00

      n1_old = n1
      n2_old = n2

    end if
!
!  Assign the node to triangles.
!
    v1 = edge_data(3,edge)
    v2 = edge_data(4,edge)
    triangle1 = edge_data(5,edge)

    if ( v1 == 1 .and. v2 == 2 ) then

      element_node2(1,(triangle1-1)*4+2) = node
      element_node2(2,(triangle1-1)*4+1) = node
      element_node2(3,(triangle1-1)*4+4) = node

    else if ( v1 == 1 .and. v2 == 3 ) then

      element_node2(1,(triangle1-1)*4+3) = node
      element_node2(2,(triangle1-1)*4+4) = node
      element_node2(3,(triangle1-1)*4+1) = node

    else if ( v1 == 2 .and. v2 == 3 ) then

      element_node2(1,(triangle1-1)*4+4) = node
      element_node2(2,(triangle1-1)*4+3) = node
      element_node2(3,(triangle1-1)*4+2) = node

    end if

  end do

  return
end
subroutine triangulation_order3_refine_size ( node_num1, element_num1, &
  element_node1, node_num2, element_num2, edge_data )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_REFINE_SIZE sizes a refined order 3 triangulation.
!
!  Discussion:
!
!    Given a triangle defined by nodes 1, 2, 3, we need to generate
!    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
!    and T4.
!
!    The task is more complicated by the fact that we are working with
!    a mesh of triangles, so that we want to create a node only once,
!    even though it may be shared by other triangles.
!
!          3
!         / \
!        /T3 \
!      13----23
!      / \T4 / \
!     /T1 \ /T2 \
!    1----12-----2
!
!    This routine simply determines the sizes of the resulting node
!    and triangle arrays.
!
!    The primary amount of work occurs in sorting a list of 3 * ELEMENT_NUM
!    data items, one item for every edge of every triangle.  Each
!    data item records, for a given edge, the global indices
!    of the two endpoints, the local indices of the two endpoints,
!    and the index of the triangle.
!
!    Through careful sorting, it is possible to arrange this data in
!    a way that allows the proper generation of the interpolated nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes in the
!    original mesh.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM1, the number of triangles in the
!    original mesh.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE1(3,ELEMENT_NUM1), the indices
!    of the nodes that form the triangles in the input mesh.
!
!    Output, integer ( kind = 4 ) NODE_NUM2, the number of nodes in the refined
!    mesh.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM2, the number of triangles in
!    the refined mesh.
!
!    Output, integer ( kind = 4 ) EDGE_DATA(5,3*ELEMENT_NUM1), edge data that
!    will be needed by TRIANGULATION_ORDER3_REFINE_COMPUTE.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) element_num1

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(5,3*element_num1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_old
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_old
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) element_node1(3,element_num1)
  integer ( kind = 4 ) element_num2
!
!  Step 1.
!  From the list of nodes for triangle T, of the form: (I,J,K)
!  construct the edge relations:
!
!    (I,J,1,2,T)
!    (I,K,1,3,T)
!    (J,K,2,3,T)
!
!  In order to make matching easier, we reorder each pair of nodes
!  into ascending order.
!
  do triangle = 1, element_num1

    i = element_node1(1,triangle)
    j = element_node1(2,triangle)
    k = element_node1(3,triangle)

    a = min ( i, j )
    b = max ( i, j )

    edge_data(1:5,3*(triangle-1)+1) = (/ a, b, 1, 2, triangle /)

    a = min ( i, k )
    b = max ( i, k )

    edge_data(1:5,3*(triangle-1)+2) = (/ a, b, 1, 3, triangle /)

    a = min ( j, k )
    b = max ( j, k )

    edge_data(1:5,3*(triangle-1)+3) = (/ a, b, 2, 3, triangle /)

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1:2; the routine we call here
!  sorts on the full column but that won't hurt us.
!
!  What we need is to find all cases where triangles share an edge.
!  By sorting the columns of the EDGE_DATA array, we will put shared edges
!  next to each other.
!
  call i4col_sort_a ( 5, 3*element_num1, edge_data )
!
!  Step 3. All the triangles which share an edge show up as consecutive
!  columns with identical first two entries.  Figure out how many new
!  nodes there are, and allocate space for their coordinates.
!
  node_num2 = node_num1

  n1_old = -1
  n2_old = -1

  do edge = 1, 3 * element_num1
    n1 = edge_data(1,edge)
    n2 = edge_data(2,edge)
    if ( n1 /= n1_old .or. n2 /= n2_old ) then
      node_num2 = node_num2 + 1
      n1_old = n1
      n2_old = n2
    end if
  end do

  element_num2 = 4 * element_num1

  return
end
subroutine triangulation_order3_sample ( node_num, node_xy, element_num, &
  element_node, num_ran, seed, xd, td )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_SAMPLE returns random points in a triangulation.
!
!  Discussion:
!
!    It is assumed that the triangulation consists of a set of non-overlapping
!    triangles.
!
!    The point is chosen uniformly in the area covered by the triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2006
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
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes that
!    make up the triangles.
!
!    Input, integer ( kind = 4 ) NUM_RAN, the number of points to sample.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!     number generator.
!
!    Output, real ( kind = 8 ) XD(2,NUM_RAN), the sample points.
!
!    Output, integer ( kind = 4 ) TD(NUM_RAN), the triangle to which each
!    sample point belongs.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) num_ran
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  real ( kind = 8 ) area
  real ( kind = 8 ) area_cum(0:element_num)
  real ( kind = 8 ) area_total
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) left
  real ( kind = 8 ) node_xy(dim_num,node_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) right
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(dim_num,3)
  integer ( kind = 4 ) td(num_ran)
  integer ( kind = 4 ) element_node(element_order,element_num)
  real ( kind = 8 ) xd(dim_num,num_ran)
!
!  Compute the areas of the triangles.
!  Build a cumulative area vector.
!  Convert it to a relative cumulative area vector.
!
  area_cum(0) = 0.0D+00

  do i = 1, element_num

    i1 = element_node(1,i)
    i2 = element_node(2,i)
    i3 = element_node(3,i)

    t(1:2,1) = node_xy(1:2,i1)
    t(1:2,2) = node_xy(1:2,i2)
    t(1:2,3) = node_xy(1:2,i3)

    call triangle_area_2d ( t, area )

    area_cum(i) = area_cum(i-1) + area

  end do

  area_total = area_cum(element_num)

  area_cum(0:element_num) = area_cum(0:element_num) / area_total
!
!  Pick random values.  A random value R indicates the corresponding triangle
!  whose cumulative relative area contains R.
!
!  Bracket the random value in the cumulative relative areas,
!  indicating a triangle.
!
!  Pick a random point in the triangle.
!
  do i = 1, num_ran

    r = r8_uniform_01 ( seed )

    call r8vec_bracket ( element_num+1, area_cum, r, left, right )

    td(i) = right - 1

    i1 = element_node(1,td(i))
    i2 = element_node(2,td(i))
    i3 = element_node(3,td(i))

    t(1:2,1) = node_xy(1:2,i1)
    t(1:2,2) = node_xy(1:2,i2)
    t(1:2,3) = node_xy(1:2,i3)

    call triangle_sample ( t, 1, seed, xd(1:2,i) )

  end do

  return
end
subroutine triangulation_order4_plot ( file_name, node_num, node_xy, &
  element_num, element_node, node_show, element_show )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER4_PLOT plots a 4-node triangulation of a pointset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(4,ELEMENT_NUM), lists,
!    for each triangle, the indices of the points that form the vertices
!    and the centroid of the triangle.
!
!    Input, integer ( kind = 4 ) NODE_SHOW,
!    0, do not show nodes;
!    1, show nodes;
!    2, show nodes and label them.
!
!    Input, integer ( kind = 4 ) ELEMENT_SHOW,
!    0, do not show triangles;
!    1, show triangles;
!    2, show triangles and label them.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num

  real ( kind = 8 ) ave_x
  real ( kind = 8 ) ave_y
  integer ( kind = 4 ) :: circle_size
  integer ( kind = 4 ) delta
  integer ( kind = 4 ) e
  character ( len = * )  file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_show
  real ( kind = 8 ) node_xy(2,node_num)
  character ( len = 40 ) string
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) element_node(4,element_num)
  integer ( kind = 4 ) element_show
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer ( kind = 4 ) x_ps
  integer ( kind = 4 ) :: x_ps_max = 576
  integer ( kind = 4 ) :: x_ps_max_clip = 594
  integer ( kind = 4 ) :: x_ps_min = 36
  integer ( kind = 4 ) :: x_ps_min_clip = 18
  real ( kind = 8 ) x_scale
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  integer ( kind = 4 ) y_ps
  integer ( kind = 4 ) :: y_ps_max = 666
  integer ( kind = 4 ) :: y_ps_max_clip = 684
  integer ( kind = 4 ) :: y_ps_min = 126
  integer ( kind = 4 ) :: y_ps_min_clip = 108
  real ( kind = 8 ) y_scale
!
!  We need to do some figuring here, so that we can determine
!  the range of the data, and hence the height and width
!  of the piece of paper.
!
  x_max = maxval ( node_xy(1,1:node_num) )
  x_min = minval ( node_xy(1,1:node_num) )
  x_scale = x_max - x_min

  x_max = x_max + 0.05D+00 * x_scale
  x_min = x_min - 0.05D+00 * x_scale
  x_scale = x_max - x_min

  y_max = maxval ( node_xy(2,1:node_num) )
  y_min = minval ( node_xy(2,1:node_num) )
  y_scale = y_max - y_min

  y_max = y_max + 0.05D+00 * y_scale
  y_min = y_min - 0.05D+00 * y_scale
  y_scale = y_max - y_min

  if ( x_scale < y_scale ) then

    delta = nint ( real ( x_ps_max - x_ps_min, kind = 8 ) &
      * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )

    x_ps_max = x_ps_max - delta
    x_ps_min = x_ps_min + delta

    x_ps_max_clip = x_ps_max_clip - delta
    x_ps_min_clip = x_ps_min_clip + delta

    x_scale = y_scale

  else if ( y_scale < x_scale ) then

    delta = nint ( real ( y_ps_max - y_ps_min, kind = 8 ) &
      * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )

    y_ps_max      = y_ps_max - delta
    y_ps_min      = y_ps_min + delta

    y_ps_max_clip = y_ps_max_clip - delta
    y_ps_min_clip = y_ps_min_clip + delta

    y_scale = x_scale

  end if

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_ORDER4_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Can not open output file.'
    return
  end if

  write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( file_unit, '(a)' ) '%%Creator: triangulation_order4_plot.f90'
  write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( file_unit, '(a)' ) '%%Pages: 1'
  write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%BoundingBox: ', &
    x_ps_min, y_ps_min, x_ps_max, y_ps_max
  write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( file_unit, '(a)' ) '%%EndComments'
  write ( file_unit, '(a)' ) '%%BeginProlog'
  write ( file_unit, '(a)' ) '/inch {72 mul} def'
  write ( file_unit, '(a)' ) '%%EndProlog'
  write ( file_unit, '(a)' ) '%%Page: 1 1'
  write ( file_unit, '(a)' ) 'save'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Increase line width from default 0.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '2 setlinewidth'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Draw a gray border around the page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_min, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' lineto'
  write ( file_unit, '(a)' ) 'stroke'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB color to black.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the font and its size.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '/Times-Roman findfont'
  write ( file_unit, '(a)' ) '0.50 inch scalefont'
  write ( file_unit, '(a)' ) 'setfont'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Print a title.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  210  702  moveto'
  write ( file_unit, '(a)' ) '%  (Triangulation)  show'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a)' ) 'clip newpath'
!
!  Draw the nodes.
!
  if ( node_num <= 200 ) then
    circle_size = 5
  else if ( node_num <= 500 ) then
    circle_size = 4
  else if ( node_num <= 1000 ) then
    circle_size = 3
  else if ( node_num <= 5000 ) then
    circle_size = 2
  else
    circle_size = 1
  end if

  if ( 1 <= node_show ) then
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw filled dots at the nodes.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to blue.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
    write ( file_unit, '(a)' ) '%'

    do node = 1, node_num

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
        + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                   - y_min ) )

      write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
        circle_size, '0 360 arc closepath fill'

    end do

  end if
!
!  Label the nodes.
!
  if ( 2 <= node_show ) then

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Label the nodes:'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to darker blue.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.000  0.250  0.850 setrgbcolor'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.20 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'
    write ( file_unit, '(a)' ) '%'

    do node = 1, node_num

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
        + (       + node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                   - y_min ) )

      write ( string, '(i4)' ) node
      string = adjustl ( string )

      write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, &
        ' moveto (' // trim ( string ) // ') show'

    end do

  end if
!
!  Draw the triangles.
!
  if ( 1 <= element_show ) then
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to red.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw the triangles.'
    write ( file_unit, '(a)' ) '%'

    do triangle = 1, element_num

      write ( file_unit, '(a)' ) 'newpath'

      do i = 1, 4

        e = i4_wrap ( i, 1, 3 )

        node = element_node(e,triangle)

        x_ps = int ( &
          ( ( x_max - node_xy(1,node)         ) &
          * real ( x_ps_min, kind = 8 )   &
          + (         node_xy(1,node) - x_min ) &
          * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                   - x_min ) )

        y_ps = int ( &
          ( ( y_max - node_xy(2,node)         ) &
          * real ( y_ps_min, kind = 8 )   &
          + (         node_xy(2,node) - y_min ) &
          * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                   - y_min ) )

        if ( i == 1 ) then
          write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' moveto'
        else
          write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' lineto'
        end if

      end do

      write ( file_unit, '(a)' ) 'stroke'

    end do

  end if
!
!  Label the triangles.
!
  if ( 2 <= element_show ) then

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Label the triangles:'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to darker red.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.950  0.250  0.150 setrgbcolor'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.20 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'
    write ( file_unit, '(a)' ) '%'

    do triangle = 1, element_num

      ave_x = 0.0D+00
      ave_y = 0.0D+00

      do i = 1, 3

        node = element_node(i,triangle)

        ave_x = ave_x + node_xy(1,node)
        ave_y = ave_y + node_xy(2,node)

      end do

      ave_x = ave_x / 3.0D+00
      ave_y = ave_y / 3.0D+00

      x_ps = int ( &
        ( ( x_max - ave_x         ) * real ( x_ps_min, kind = 8 )   &
        + (       + ave_x - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max         - x_min ) )

      y_ps = int ( &
        ( ( y_max - ave_y         ) * real ( y_ps_min, kind = 8 )   &
        + (         ave_y - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max         - y_min ) )

      write ( string, '(i4)' ) triangle
      string = adjustl ( string )

      write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps, ' moveto (' &
        // trim ( string ) // ') show'

    end do

  end if

  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'restore  showpage'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  End of page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%%Trailer'
  write ( file_unit, '(a)' ) '%%EOF'
  close ( unit = file_unit )

  return
end
subroutine triangulation_order6_adj_count ( node_num, element_num, &
  element_node, element_neighbor, adj_num, adj_col )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_ADJ_COUNT counts adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The triangulation is assumed to involve 6-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!    Also, a node is considered to be adjacent to itself.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  6   5  side 2
!       |    \
!    3  |     \
!       |      \
!       1---4---2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\    |\    |
!    | \   | \   |
!   16 17 18 19 20
!    |   \ |   \ |
!    |    \|    \|
!   11-12-13-14-15
!    |\    |\    |
!    | \   | \   |
!    6  7  8  9 10
!    |   \ |   \ |
!    |    \|    \|
!    1--2--3--4--5
!
!    A sample grid.
!
!
!    Below, we have a chart that lists the nodes adjacent to each node, with
!    an asterisk to indicate the adjacency of the node to itself
!    (in some cases, you want to count this self adjacency and in some
!    you don't).
!
!    N   Adjacencies
!
!    1:  *  2  3  6  7 11
!    2:  1  *  3  6  7 11
!    3:  1  2  *  4  5  6  7  8  9 11 12 13
!    4:  3  *  5  8  9 13
!    5:  3  4  *  8  9 10 13 14 15
!    6:  1  2  3  *  7 11
!    7:  1  2  3  6  *  8 11 12 13
!    8:  3  4  5  7  *  9 11 12 13
!    9:  3  4  5  8  * 10 13 14 15
!   10:  5  9  * 13 14 15
!   11:  1  2  3  6  7  8  * 12 13 16 17 21
!   12:  3  7  8 11  * 13 16 17 21
!   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
!   14:  5  9 10 13  * 15 18 19 23
!   15:  5  9 10 13 14  * 18 19 20 23 24 25
!   16: 11 12 13  * 17 21
!   17: 11 12 13 16  * 18 21 22 23
!   18: 13 14 15 17  * 19 21 22 23
!   19: 13 14 15 18  * 20 23 24 25
!   20: 15 19  * 23 24 25
!   21: 11 12 13 16 17 18  * 22 23
!   22: 13 17 18 21  * 23
!   23: 13 14 15 17 18 19 20 21 22  * 24 25
!   24: 15 19 20 23  * 25
!   25: 15 19 20 23 24  *
!
!    Below, we list the number of adjacencies to lower-indexed nodes, to
!    the node itself, to higher-indexed nodes, the total number of
!    adjacencies for this node, and the location of the first and last
!    entries required to list this set of adjacencies in a single list
!    of all the adjacencies.
!
!    N   Below  Self   Above   Total First  Last
!
!   --      --    --      --      --   ---     0
!    1:      0     1       5       6     1     6
!    2:      1     1       4       6     7    12
!    3:      2     1       9      12    13    24
!    4:      1     1       4       6    25    30
!    5:      2     1       6       9    31    39
!    6:      3     1       2       6    40    45
!    7:      4     1       4       9    46    54
!    8:      4     1       4       9    55    63
!    9:      4     1       4       9    62    72
!   10:      2     1       3       6    73    78
!   11:      6     1       5      12    79    90
!   12:      4     1       4       9    91    99
!   13:      9     1       9      19   100   118
!   14:      4     1       4       9   119   127
!   15:      5     1       6      12   128   139
!   16:      3     1       2       6   140   145
!   17:      4     1       4       9   146   154
!   18:      4     1       4       9   155   163
!   19:      4     1       4       9   164   172
!   20:      2     1       3       6   173   178
!   21:      6     1       2       9   179   187
!   22:      4     1       1       6   188   193
!   23:      9     1       2      12   194   205
!   24:      4     1       1       6   206   211
!   25:      5     1       0       6   212   217
!   --      --    --      --      --   218   ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), lists the
!    nodes that make up each triangle.  The first three nodes are the vertices,
!    in counterclockwise order.  The fourth value is the midside
!    node between nodes 1 and 2; the fifth and sixth values are
!    the other midside nodes in the logical order.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), for each
!    side of a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Output, integer ( kind = 4 ) ADJ_NUM, the number of adjacencies.
!
!    Output, integer ( kind = 4 ) ADJ_COL(NODE_NUM+1).  Information about
!    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_col(node_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle2
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)

  adj_num = 0
!
!  Set every node to be adjacent to itself.
!
  adj_col(1:node_num) = 1
!
!  Examine each triangle.
!
  do triangle = 1, element_num

    n1 = element_node(1,triangle)
    n2 = element_node(2,triangle)
    n3 = element_node(3,triangle)
    n4 = element_node(4,triangle)
    n5 = element_node(5,triangle)
    n6 = element_node(6,triangle)
!
!  For sure, we add the adjacencies:
!    43 / (34)
!    51 / (15)
!    54 / (45)
!    62 / (26)
!    64 / (46)
!    65 / (56)
!
    adj_col(n3) = adj_col(n3) + 1
    adj_col(n4) = adj_col(n4) + 1
    adj_col(n1) = adj_col(n1) + 1
    adj_col(n5) = adj_col(n5) + 1
    adj_col(n4) = adj_col(n4) + 1
    adj_col(n5) = adj_col(n5) + 1
    adj_col(n2) = adj_col(n2) + 1
    adj_col(n6) = adj_col(n6) + 1
    adj_col(n4) = adj_col(n4) + 1
    adj_col(n6) = adj_col(n6) + 1
    adj_col(n5) = adj_col(n5) + 1
    adj_col(n6) = adj_col(n6) + 1
!
!  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
!  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
!  Maybe add
!    21 / 12
!    41 / 14
!    42 / 24
!
    triangle2 = element_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_col(n1) = adj_col(n1) + 1
      adj_col(n2) = adj_col(n2) + 1
      adj_col(n1) = adj_col(n1) + 1
      adj_col(n4) = adj_col(n4) + 1
      adj_col(n2) = adj_col(n2) + 1
      adj_col(n4) = adj_col(n4) + 1
    end if
!
!  Maybe add
!    32 / 23
!    52 / 25
!    53 / 35
!
    triangle2 = element_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_col(n2) = adj_col(n2) + 1
      adj_col(n3) = adj_col(n3) + 1
      adj_col(n2) = adj_col(n2) + 1
      adj_col(n5) = adj_col(n5) + 1
      adj_col(n3) = adj_col(n3) + 1
      adj_col(n5) = adj_col(n5) + 1
    end if
!
!  Maybe add
!    31 / 13
!    61 / 16
!    63 / 36
!
    triangle2 = element_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_col(n1) = adj_col(n1) + 1
      adj_col(n3) = adj_col(n3) + 1
      adj_col(n1) = adj_col(n1) + 1
      adj_col(n6) = adj_col(n6) + 1
      adj_col(n3) = adj_col(n3) + 1
      adj_col(n6) = adj_col(n6) + 1
    end if

  end do
!
!  We used ADJ_COL to count the number of entries in each column.
!  Convert it to pointers into the ADJ array.
!
  adj_col(2:node_num+1) = adj_col(1:node_num)

  adj_col(1) = 1
  do i = 2, node_num+1
    adj_col(i) = adj_col(i-1) + adj_col(i)
  end do

  adj_num = adj_col(node_num+1) - 1

  return
end
subroutine triangulation_order6_adj_set ( node_num, element_num, &
  element_node, element_neighbor, adj_num, adj_col, adj )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_ADJ_SET sets adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The triangulation is assumed to involve 6-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!    Also, a node is considered to be adjacent to itself.
!
!    This routine can be used to create the compressed column storage
!    for a quadratic triangle finite element discretization of
!    Poisson's equation in two dimensions.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  6   5  side 2
!       |    \
!    3  |     \
!       |      \
!       1---4---2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\    |\    |
!    | \   | \   |
!   16 17 18 19 20
!    |   \ |   \ |
!    |    \|    \|
!   11-12-13-14-15
!    |\    |\    |
!    | \   | \   |
!    6  7  8  9 10
!    |   \ |   \ |
!    |    \|    \|
!    1--2--3--4--5
!
!    A sample grid.
!
!
!    Below, we have a chart that lists the nodes adjacent to each node, with
!    an asterisk to indicate the adjacency of the node to itself
!    (in some cases, you want to count this self adjacency and in some
!    you don't).
!
!    N   Adjacencies
!
!    1:  *  2  3  6  7 11
!    2:  1  *  3  6  7 11
!    3:  1  2  *  4  5  6  7  8  9 11 12 13
!    4:  3  *  5  8  9 13
!    5:  3  4  *  8  9 10 13 14 15
!    6:  1  2  3  *  7 11
!    7:  1  2  3  6  *  8 11 12 13
!    8:  3  4  5  7  *  9 11 12 13
!    9:  3  4  5  8  * 10 13 14 15
!   10:  5  9  * 13 14 15
!   11:  1  2  3  6  7  8  * 12 13 16 17 21
!   12:  3  7  8 11  * 13 16 17 21
!   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
!   14:  5  9 10 13  * 15 18 19 23
!   15:  5  9 10 13 14  * 18 19 20 23 24 25
!   16: 11 12 13  * 17 21
!   17: 11 12 13 16  * 18 21 22 23
!   18: 13 14 15 17  * 19 21 22 23
!   19: 13 14 15 18  * 20 23 24 25
!   20: 15 19  * 23 24 25
!   21: 11 12 13 16 17 18  * 22 23
!   22: 13 17 18 21  * 23
!   23: 13 14 15 17 18 19 20 21 22  * 24 25
!   24: 15 19 20 23  * 25
!   25: 15 19 20 23 24  *
!
!    Below, we list the number of adjacencies to lower-indexed nodes, to
!    the node itself, to higher-indexed nodes, the total number of
!    adjacencies for this node, and the location of the first and last
!    entries required to list this set of adjacencies in a single list
!    of all the adjacencies.
!
!    N   Below  Self   Above   Total First  Last
!
!   --      --    --      --      --   ---     0
!    1:      0     1       5       6     1     6
!    2:      1     1       4       6     7    12
!    3:      2     1       9      12    13    24
!    4:      1     1       4       6    25    30
!    5:      2     1       6       9    31    39
!    6:      3     1       2       6    40    45
!    7:      4     1       4       9    46    54
!    8:      4     1       4       9    55    63
!    9:      4     1       4       9    62    72
!   10:      2     1       3       6    73    78
!   11:      6     1       5      12    79    90
!   12:      4     1       4       9    91    99
!   13:      9     1       9      19   100   118
!   14:      4     1       4       9   119   127
!   15:      5     1       6      12   128   139
!   16:      3     1       2       6   140   145
!   17:      4     1       4       9   146   154
!   18:      4     1       4       9   155   163
!   19:      4     1       4       9   164   172
!   20:      2     1       3       6   173   178
!   21:      6     1       2       9   179   187
!   22:      4     1       1       6   188   193
!   23:      9     1       2      12   194   205
!   24:      4     1       1       6   206   211
!   25:      5     1       0       6   212   217
!   --      --    --      --      --   218   ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), lists the nodes
!    that make up each triangle.  The first three nodes are the vertices,
!    in counterclockwise order.  The fourth value is the midside
!    node between nodes 1 and 2; the fifth and sixth values are
!    the other midside nodes in the logical order.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), for each
!    side of a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacencies.
!
!    Input, integer ( kind = 4 ) ADJ_COL(NODE_NUM+1).  Information about column
!    J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
!
!    Output, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency information.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_copy(node_num)
  integer ( kind = 4 ) adj_col(node_num+1)
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  integer ( kind = 4 ) node
  integer ( kind = 4 ) number
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle2
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)

  adj(1:adj_num) = -1
  adj_copy(1:node_num) = adj_col(1:node_num)
!
!  Set every node to be adjacent to itself.
!
  do node = 1, node_num
    adj(adj_copy(node)) = node
    adj_copy(node) = adj_copy(node) + 1
  end do
!
!  Examine each triangle.
!
  do triangle = 1, element_num

    n1 = element_node(1,triangle)
    n2 = element_node(2,triangle)
    n3 = element_node(3,triangle)
    n4 = element_node(4,triangle)
    n5 = element_node(5,triangle)
    n6 = element_node(6,triangle)
!
!  For sure, we add the adjacencies:
!    43 / (34)
!    51 / (15)
!    54 / (45)
!    62 / (26)
!    64 / (46)
!    65 / (56)
!
    adj(adj_copy(n3)) = n4
    adj_copy(n3) = adj_copy(n3) + 1
    adj(adj_copy(n4)) = n3
    adj_copy(n4) = adj_copy(n4) + 1

    adj(adj_copy(n1)) = n5
    adj_copy(n1) = adj_copy(n1) + 1
    adj(adj_copy(n5)) = n1
    adj_copy(n5) = adj_copy(n5) + 1

    adj(adj_copy(n4)) = n5
    adj_copy(n4) = adj_copy(n4) + 1
    adj(adj_copy(n5)) = n4
    adj_copy(n5) = adj_copy(n5) + 1

    adj(adj_copy(n2)) = n6
    adj_copy(n2) = adj_copy(n2) + 1
    adj(adj_copy(n6)) = n2
    adj_copy(n6) = adj_copy(n6) + 1

    adj(adj_copy(n4)) = n6
    adj_copy(n4) = adj_copy(n4) + 1
    adj(adj_copy(n6)) = n4
    adj_copy(n6) = adj_copy(n6) + 1

    adj(adj_copy(n5)) = n6
    adj_copy(n5) = adj_copy(n5) + 1
    adj(adj_copy(n6)) = n5
    adj_copy(n6) = adj_copy(n6) + 1
!
!  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
!  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
!  Maybe add
!    21 / 12
!    41 / 14
!    42 / 24
!
    triangle2 = element_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n1)) = n2
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n2)) = n1
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n1)) = n4
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n4)) = n1
      adj_copy(n4) = adj_copy(n4) + 1
      adj(adj_copy(n2)) = n4
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n4)) = n2
      adj_copy(n4) = adj_copy(n4) + 1
    end if
!
!  Maybe add
!    32 / 23
!    52 / 25
!    53 / 35
!
    triangle2 = element_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n2)) = n3
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n3)) = n2
      adj_copy(n3) = adj_copy(n3) + 1
      adj(adj_copy(n2)) = n5
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n5)) = n2
      adj_copy(n5) = adj_copy(n5) + 1
      adj(adj_copy(n3)) = n5
      adj_copy(n3) = adj_copy(n3) + 1
      adj(adj_copy(n5)) = n3
      adj_copy(n5) = adj_copy(n5) + 1
    end if
!
!  Maybe add
!    31 / 13
!    61 / 16
!    63 / 36
!
    triangle2 = element_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n1)) = n3
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n3)) = n1
      adj_copy(n3) = adj_copy(n3) + 1
      adj(adj_copy(n1)) = n6
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n6)) = n1
      adj_copy(n6) = adj_copy(n6) + 1
      adj(adj_copy(n3)) = n6
      adj_copy(n3) = adj_copy(n3) + 1
      adj(adj_copy(n6)) = n3
      adj_copy(n6) = adj_copy(n6) + 1
    end if

  end do
!
!  Ascending sort the entries for each node.
!
  do node = 1, node_num
    k1 = adj_col(node)
    k2 = adj_col(node+1)-1
    number = k2 + 1 - k1
    call i4vec_sort_heap_a ( number, adj(k1:k2) )
  end do

  return
end
subroutine triangulation_order6_boundary_edge_count ( element_num, &
  element_node, boundary_edge_num )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT counts the boundary edges.
!
!  Discussion:
!
!    This routine is given a triangulation, a set of 6-node triangles.
!    It is assumed that, in each list of 6 nodes, the vertices are listed
!    first, in counterclockwise order, followed by the three midside nodes,
!    in counterclockwise order, starting with the node between vertices
!    1 and 2.
!
!    It is assumed that each edge of the triangulation is either
!    * an INTERIOR edge, which is listed twice, once with positive
!      orientation and once with negative orientation, or;
!    * a BOUNDARY edge, which will occur only once.
!
!    This routine should work even if the region has holes - as long
!    as the boundary of the hole comprises more than 3 edges!
!
!    Except for the dimension of TRIANGLE, this routine is identical
!    to the routine for the order 3 case.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), the nodes that
!    make up the triangles.  These should be listed in counterclockwise order.
!
!    Output, integer ( kind = 4 ) BOUNDARY_EDGE_NUM, the number of boundary
!    edges.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6

  integer ( kind = 4 ) boundary_edge_num
  integer ( kind = 4 ) e1(3*element_num)
  integer ( kind = 4 ) e2(3*element_num)
  integer ( kind = 4 ) edge(2,3*element_num)
  integer ( kind = 4 ) interior_edge_num
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) unique_num

  m = 2
  n = 3 * element_num
!
!  Set up the edge array.
!
  edge(1:2,               1:  element_num) = element_node(1:2,1:element_num)
  edge(1:2,  element_num+1:2*element_num) = element_node(2:3,1:element_num)
  edge(1  ,2*element_num+1:3*element_num) = element_node(3,  1:element_num)
  edge(2  ,2*element_num+1:3*element_num) = element_node(1,  1:element_num)
!
!  In each column, force the smaller entry to appear first.
!
  e1(1:n) = minval ( edge(1:2,1:n), dim = 1 )
  e2(1:n) = maxval ( edge(1:2,1:n), dim = 1 )

  edge(1,1:n) = e1(1:n)
  edge(2,1:n) = e2(1:n)
!
!  Ascending sort the column array.
!
  call i4col_sort_a ( m, n, edge )
!
!  Get the number of unique columns in EDGE.
!
  call i4col_sorted_unique_count ( m, n, edge, unique_num )

  interior_edge_num = 3 * element_num - unique_num

  boundary_edge_num = 3 * element_num - 2 * interior_edge_num

  return
end
subroutine triangulation_order6_boundary_edge_count_euler ( node_num, &
  element_num, hole_num, boundary_num )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER counts boundary edges.
!
!  Discussion:
!
!    We assume we are given information about an order 6 triangulation
!    of a set of nodes in the plane.
!
!    By ignoring the midside nodes, we can determine the corresponding
!    information for an order 3 triangulation, and apply Euler's formula
!    to determine the number of edges that lie on the boundary of the
!    set of nodes.
!
!    Thus, if we have ELEMENT_NUM triangles, and NODE_NUM nodes, we
!    imagine that each triangle is replaced by 4 triangles, created
!    by adding the edges created by joining the midside nodes.
!
!    Thus, for 4 * ELEMENT_NUM triangles, we can apply Euler's formula
!    to compute the number of boundary edges.
!
!    Now, to adjust the data to our order 6 triangles, we divide the
!    number of boundary edges by 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc Krevald, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry, Section 9.1,
!    Springer, 2000.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of internal nodes.
!
!    Output, integer ( kind = 4 ) BOUNDARY_NUM, the number of edges that lie
!    on the boundary of the triangulation.
!
  implicit none

  integer ( kind = 4 ) boundary_num
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num

  boundary_num = ( 2 * node_num + 2 * hole_num - 4 * element_num - 2 ) / 2

  return
end
subroutine triangulation_order6_boundary_node ( node_num, element_num, &
  element_node, node_boundary )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_BOUNDARY_NODE indicates which nodes are on the boundary.
!
!  Discussion:
!
!    This routine is given a triangulation, an abstract list of triples
!    of nodes.  It is assumed that the nodes in each triangle are listed
!    in a counterclockwise order, although the routine should work
!    if the nodes are consistently listed in a clockwise order as well.
!
!    It is assumed that each edge of the triangulation is either
!    * an INTERIOR edge, which is listed twice, once with positive
!      orientation and once with negative orientation, or;
!    * a BOUNDARY edge, which will occur only once.
!
!    This routine should work even if the region has holes - as long
!    as the boundary of the hole comprises more than 3 edges!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), the nodes
!    that make up the triangles.
!
!    Output, logical NODE_BOUNDARY(NODE_NUM), is TRUE if the node
!    is on a boundary edge.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6

  integer ( kind = 4 ) e1(3*element_num)
  integer ( kind = 4 ) e2(3*element_num)
  integer ( kind = 4 ) edge(3,3*element_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  logical              node_boundary(node_num)
  integer ( kind = 4 ) element_node(element_order,element_num)

  m = 3
  n = 3 * element_num
!
!  Set up the edge array.  The midside node is listed last, as
!  it is not needed for the sorting process.
!
  edge(1,               1:  element_num) = element_node(1,1:element_num)
  edge(2,               1:  element_num) = element_node(4,1:element_num)
  edge(3,               1:  element_num) = element_node(2,1:element_num)

  edge(1,  element_num+1:2*element_num) = element_node(2,1:element_num)
  edge(2,  element_num+1:2*element_num) = element_node(5,1:element_num)
  edge(3,  element_num+1:2*element_num) = element_node(3,1:element_num)

  edge(1,2*element_num+1:3*element_num) = element_node(3,1:element_num)
  edge(2,2*element_num+1:3*element_num) = element_node(6,1:element_num)
  edge(3,2*element_num+1:3*element_num) = element_node(1,1:element_num)
!
!  In each column, force the smaller of the two vertices to appear first.
!
  e1(1:n) = minval ( edge(1:3:2,1:n), dim = 1 )
  e2(1:n) = maxval ( edge(1:3:2,1:n), dim = 1 )

  edge(1,1:n) = e1(1:n)
  edge(3,1:n) = e2(1:n)
!
!  Ascending sort the column array.
!
  call i4col_sort_a ( m, n, edge )
!
!  Records which appear twice are internal edges and can be ignored.
!
  node_boundary(1:node_num) = .false.

  i = 0

  do while ( i < 3 * element_num )

    i = i + 1

    if ( i == 3 * element_num ) then
      node_boundary(edge(1:m,i)) = .true.
    else if ( all ( edge(1:m,i) == edge(1:m,i+1) ) ) then
      i = i + 1
    else
      node_boundary(edge(1:m,i)) = .true.
    end if

  end do

  return
end
subroutine triangulation_order6_example1 ( node_num, element_num, node_xy, &
  element_node, element_neighbor )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_EXAMPLE1 sets up a sample triangulation.
!
!  Discussion:
!
!    This triangulation is actually a Delaunay triangulation.
!
!    The appropriate input values of NODE_NUM and ELEMENT_NUM can be
!    determined by calling TRIANGULATION_ORDER6_EXAMPLE1_SIZE first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Output, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of
!    the nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_ORDER(6,ELEMENT_NUM), the nodes
!    that make up the triangles.
!
!    Output, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbors on each side.  Negative values indicate edges that
!    lie on the exterior.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6

  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) element_neighbor(3,element_num)

  node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       1.0D+00, 0.0D+00, &
       2.0D+00, 0.0D+00, &
       3.0D+00, 0.0D+00, &
       4.0D+00, 0.0D+00, &
       5.0D+00, 0.0D+00, &
       6.0D+00, 0.0D+00, &
       0.0D+00, 1.0D+00, &
       1.0D+00, 1.0D+00, &
       2.0D+00, 1.0D+00, &
       3.0D+00, 1.0D+00, &
       4.0D+00, 1.0D+00, &
       5.0D+00, 1.0D+00, &
       6.0D+00, 1.0D+00, &
       0.0D+00, 2.0D+00, &
       1.0D+00, 2.0D+00, &
       2.0D+00, 2.0D+00, &
       3.0D+00, 2.0D+00, &
       4.0D+00, 2.0D+00, &
       5.0D+00, 2.0D+00, &
       6.0D+00, 2.0D+00, &
       0.0D+00, 3.0D+00, &
       1.0D+00, 3.0D+00, &
       2.0D+00, 3.0D+00, &
       3.0D+00, 3.0D+00, &
       5.0D+00, 3.0D+00, &
       6.0D+00, 3.0D+00, &
       0.0D+00, 4.0D+00, &
       1.0D+00, 4.0D+00, &
       2.0D+00, 4.0D+00, &
       3.0D+00, 4.0D+00, &
       4.0D+00, 4.0D+00, &
       5.0D+00, 4.0D+00, &
       6.0D+00, 4.0D+00, &
       0.0D+00, 5.0D+00, &
       1.0D+00, 5.0D+00, &
       2.0D+00, 5.0D+00, &
       3.0D+00, 5.0D+00, &
       4.0D+00, 5.0D+00, &
       5.0D+00, 5.0D+00, &
       6.0D+00, 5.0D+00, &
       0.0D+00, 6.0D+00, &
       1.0D+00, 6.0D+00, &
       2.0D+00, 6.0D+00, &
       3.0D+00, 6.0D+00, &
       4.0D+00, 6.0D+00, &
       5.0D+00, 6.0D+00, &
       6.0D+00, 6.0D+00 /), (/ dim_num, node_num /) )

  element_node(1:element_order,1:element_num ) = reshape ( (/ &
     1,  3, 15,  2,  9,  8, &
    17, 15,  3, 16,  9, 10, &
     5, 19,  3, 12, 11,  4, &
    17,  3, 19, 10, 11, 18, &
     7, 21,  5, 14, 13,  6, &
    19,  5, 21, 12, 13, 20, &
    17, 30, 15, 24, 23, 16, &
    28, 15, 30, 22, 23, 29, &
    30, 17, 32, 24, 25, 31, &
    21, 34, 19, 27, 26, 20, &
    30, 44, 28, 37, 36, 29, &
    42, 28, 44, 35, 36, 43, &
    32, 46, 30, 39, 38, 31, &
    44, 30, 46, 37, 38, 45, &
    32, 34, 46, 33, 40, 39, &
    48, 46, 34, 47, 40, 41 /), (/ element_order, element_num /) )

  element_neighbor(1:3,1:element_num) = reshape ( (/ &
    -3,   2,  -5, &
     7,   1,   4, &
     6,   4, -11, &
     2,   3, -14, &
   -15,   6, -17, &
     3,   5,  10, &
     9,   8,   2, &
   -24,   7,  11, &
     7, -28,  13, &
    27, -31,   6, &
     8,  14,  12, &
   -36,  11, -38, &
    15,  14,   9, &
    11,  13, -44, &
   -45,  16,  13, &
   -48,  15, -50 /), (/ 3, element_num /) )

  return
end
subroutine triangulation_order6_example1_size ( node_num, element_num, &
  hole_num )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_EXAMPLE1_SIZE sets sizes for a sample triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num

  node_num = 48
  element_num = 16
  hole_num = 1

  return
end
subroutine triangulation_order6_example2 ( node_num, element_num, node_xy, &
  element_node, element_neighbor )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_EXAMPLE2 returns an example triangulation.
!
!  Discussion:
!
!    This triangulation is actually a Delaunay triangulation.
!
!    The appropriate input values of NODE_NUM and ELEMENT_NUM can be
!    determined by calling TRIANGULATION_ORDER6_EXAMPLE2_SIZE first.
!
!  Diagram:
!
!   21-22-23-24-25
!    |\  6 |\  8 |
!    | \   | \   |
!   16 17 18 19 20
!    |   \ |   \ |
!    | 5  \| 7  \|
!   11-12-13-14-15
!    |\  2 |\  4 |
!    | \   | \   |
!    6  7  8  9 10
!    | 1 \ | 3 \ |
!    |    \|    \|
!    1--2--3--4--5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Output, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of
!    the nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), lists the
!    nodes that make up each triangle.  The first three nodes are the vertices,
!    in counterclockwise order.  The fourth value is the midside
!    node between nodes 1 and 2; the fifth and sixth values are
!    the other midside nodes in the logical order.
!
!    Output, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), for each
!    side of a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6

  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) element_node(element_order,element_num)

  node_xy = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    2.0D+00, 0.0D+00, &
    3.0D+00, 0.0D+00, &
    4.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, &
    2.0D+00, 1.0D+00, &
    3.0D+00, 1.0D+00, &
    4.0D+00, 1.0D+00, &
    0.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00, &
    2.0D+00, 2.0D+00, &
    3.0D+00, 2.0D+00, &
    4.0D+00, 2.0D+00, &
    0.0D+00, 3.0D+00, &
    1.0D+00, 3.0D+00, &
    2.0D+00, 3.0D+00, &
    3.0D+00, 3.0D+00, &
    4.0D+00, 3.0D+00, &
    0.0D+00, 4.0D+00, &
    1.0D+00, 4.0D+00, &
    2.0D+00, 4.0D+00, &
    3.0D+00, 4.0D+00, &
    4.0D+00, 4.0D+00  &
  /), (/ dim_num, node_num /) )

  element_node(1:element_order,1:element_num) = reshape ( (/ &
     1,  3, 11,  2,  7,  6, &
    13, 11,  3, 12,  7,  8, &
     3,  5, 13,  4,  9,  8, &
    15, 13,  5, 14,  9, 10, &
    11, 13, 21, 12, 17, 16, &
    23, 21, 13, 22, 17, 18, &
    13, 15, 23, 14, 19, 18, &
    25, 23, 15, 24, 19, 20 /), (/ element_order, element_num /) )

  element_neighbor(1:3,1:element_num) = reshape ( (/ &
    -1,  2, -1, &
     5,  1,  3, &
    -1,  4,  2, &
     7,  3, -1, &
     2,  6, -1, &
    -1,  5,  7, &
     4,  8,  6, &
    -1,  7, -1 /), (/ 3, element_num /) )

  return
end
subroutine triangulation_order6_example2_size ( node_num, element_num, &
  hole_num )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_EXAMPLE2_SIZE returns the size of an example.
!
!  Diagram:
!
!   21-22-23-24-25
!    |\  6 |\  8 |
!    | \   | \   |
!   16 17 18 19 20
!    |   \ |   \ |
!    | 5  \| 7  \|
!   11-12-13-14-15
!    |\  2 |\  4 |
!    | \   | \   |
!    6  7  8  9 10
!    | 1 \ | 3 \ |
!    |    \|    \|
!    1--2--3--4--5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num

  node_num = 25
  element_num = 8
  hole_num = 0

  return
end
subroutine triangulation_order6_neighbor ( element_num, element_node, &
  t1, s1, t2, s2 )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_NEIGHBOR determines a neighbor of a given triangle.
!
!  Discussion:
!
!    A set of nodes is given.  A triangulation of the nodes has been
!    defined and recorded in TRIANGLE.  The TRIANGLE data structure records
!    triangles as sets of six nodes, with the first three being the
!    vertices, in counterclockwise order.  The fourth node is the
!    midside node between nodes 1 and 2, and the other two are listed
!    logically.
!
!    The nodes of the triangle are listed in counterclockwise order.
!    This means that if two triangles share a side, then the nodes
!    defining that side occur in the order (N1,N2) for one triangle,
!    and (N2,N1) for the other.
!
!    The routine is given a triangle and a side, and asked to find
!    another triangle (if any) that shares that side.  The routine
!    simply searches the TRIANGLE structure for an occurrence of the
!    nodes in the opposite order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER(6,ELEMENT_NUM), the nodes
!    that define each triangle.
!
!    Input, integer ( kind = 4 ) T1, the index of the triangle.
!
!    Input, integer ( kind = 4 ) S1, the index of the triangle side.
!
!    Output, integer ( kind = 4 ) T2, the index of the triangle which is the
!    neighbor to T1 on side S1, or -1 if there is no such neighbor.
!
!    Output, integer ( kind = 4 ) S2, the index of the side of triangle T2 which
!    is shared with triangle T1, or -1 if there is no such neighbor.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6

  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) s
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) ss
  integer ( kind = 4 ) t
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2
  integer ( kind = 4 ) element_node(element_order,element_num)

  t2 = -1
  s2 = -1

  n1 = element_node(s1,t1)
  ss = s1 + 1
  ss = i4_wrap ( ss, 1, 3 )
  n2 = element_node(ss,t1)

  do t = 1, element_num
    do s = 1, 3
      if ( element_node(s,t) == n1 ) then
        ss = s - 1
        ss = i4_wrap ( ss, 1, 3 )
        if ( element_node(ss,t) == n2 ) then
          t2 = t
          s2 = ss
          return
        end if
      end if
    end do
  end do

  return
end
subroutine triangulation_order6_plot ( file_name, node_num, node_xy, &
  element_num, element_node, node_show, element_show )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_PLOT plots a 6-node triangulation of a set of nodes.
!
!  Discussion:
!
!    The triangulation is most usually a Delaunay triangulation,
!    but this is not necessary.
!
!    In a six node triangulation, it is assumed that nodes 1, 2, and 3
!    are the vertices of the triangles, and that nodes 4, 5, and 6
!    lie between 1 and 2, 2 and 3, and 3 and 1 respectively.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), lists, for
!    each triangle, the indices of the nodes that form the vertices of the
!    triangle.
!
!    Input, integer ( kind = 4 ) NODE_SHOW,
!    0, do not show nodes;
!    1, show nodes;
!    2, show nodes and label them.
!
!    Input, integer ( kind = 4 ) ELEMENT_SHOW,
!    0, do not show triangles;
!    1, show triangles;
!    2, show triangles and label them.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6

  real ( kind = 8 ) ave_x
  real ( kind = 8 ) ave_y
  integer ( kind = 4 ) :: circle_size
  integer ( kind = 4 ) delta
  character ( len = * )  file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_show
  real ( kind = 8 ) node_xy(2,node_num)
  character ( len = 40 ) string
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) element_show
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer ( kind = 4 ) x_ps
  integer ( kind = 4 ) :: x_ps_max = 576
  integer ( kind = 4 ) :: x_ps_max_clip = 594
  integer ( kind = 4 ) :: x_ps_min = 36
  integer ( kind = 4 ) :: x_ps_min_clip = 18
  real ( kind = 8 ) x_scale
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  integer ( kind = 4 ) y_ps
  integer ( kind = 4 ) :: y_ps_max = 666
  integer ( kind = 4 ) :: y_ps_max_clip = 684
  integer ( kind = 4 ) :: y_ps_min = 126
  integer ( kind = 4 ) :: y_ps_min_clip = 108
  real ( kind = 8 ) y_scale
!
!  We need to do some figuring here, so that we can determine
!  the range of the data, and hence the height and width
!  of the piece of paper.
!
  x_max = maxval ( node_xy(1,1:node_num) )
  x_min = minval ( node_xy(1,1:node_num) )
  x_scale = x_max - x_min

  x_max = x_max + 0.05D+00 * x_scale
  x_min = x_min - 0.05D+00 * x_scale
  x_scale = x_max - x_min

  y_max = maxval ( node_xy(2,1:node_num) )
  y_min = minval ( node_xy(2,1:node_num) )
  y_scale = y_max - y_min

  y_max = y_max + 0.05D+00 * y_scale
  y_min = y_min - 0.05D+00 * y_scale
  y_scale = y_max - y_min

  if ( x_scale < y_scale ) then

    delta = nint ( real ( x_ps_max - x_ps_min, kind = 8 ) &
      * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )

    x_ps_max = x_ps_max - delta
    x_ps_min = x_ps_min + delta

    x_ps_max_clip = x_ps_max_clip - delta
    x_ps_min_clip = x_ps_min_clip + delta

    x_scale = y_scale

  else if ( y_scale < x_scale ) then

    delta = nint ( real ( y_ps_max - y_ps_min, kind = 8 ) &
      * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )

    y_ps_max      = y_ps_max - delta
    y_ps_min      = y_ps_min + delta

    y_ps_max_clip = y_ps_max_clip - delta
    y_ps_min_clip = y_ps_min_clip + delta

    y_scale = x_scale

  end if
!
!  Open the file.
!
  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_ORDER6_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Can not open output file.'
    return
  end if
!
!  Write the header.
!
  write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( file_unit, '(a)' ) '%%Creator: triangulation_order6_plot.f90'
  write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( file_unit, '(a)' ) '%%Pages: 1'
  write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%BoundingBox: ', &
    x_ps_min, y_ps_min, x_ps_max, y_ps_max
  write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( file_unit, '(a)' ) '%%EndComments'
  write ( file_unit, '(a)' ) '%%BeginProlog'
  write ( file_unit, '(a)' ) '/inch {72 mul} def'
  write ( file_unit, '(a)' ) '%%EndProlog'
  write ( file_unit, '(a)' ) '%%Page: 1 1'
  write ( file_unit, '(a)' ) 'save'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Increase line width from default 0.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '2 setlinewidth'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Draw a gray border around the page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_min, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' lineto'
  write ( file_unit, '(a)' ) 'stroke'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to black.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the font and its size.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '/Times-Roman findfont'
  write ( file_unit, '(a)' ) '0.50 inch scalefont'
  write ( file_unit, '(a)' ) 'setfont'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Print a title.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  210  702  moveto'
  write ( file_unit, '(a)' ) '%  (Triangulation)  show'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a)' ) 'clip newpath'
!
!  Draw the nodes.
!
  if ( node_num <= 200 ) then
    circle_size = 5
  else if ( node_num <= 500 ) then
    circle_size = 4
  else if ( node_num <= 1000 ) then
    circle_size = 3
  else if ( node_num <= 5000 ) then
    circle_size = 2
  else
    circle_size = 1
  end if

  if ( 1 <= node_show ) then
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw filled dots at the nodes.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to blue.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
    write ( file_unit, '(a)' ) '%'

    do node = 1, node_num

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
        + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                   - y_min ) )

      write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
        circle_size, '0 360 arc closepath fill'

    end do

  end if
!
!  Label the nodes.
!
  if ( 2 <= node_show ) then

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Label the nodes:'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to darker blue.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.100  0.250  0.850 setrgbcolor'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.20 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'

    do node = 1, node_num

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
        + (       + node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                   - y_min ) )

      write ( string, '(i4)' ) node
      string = adjustl ( string )

      write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, &
        ' moveto (' // trim ( string ) // ') show'

    end do

  end if
!
!  Draw the triangles.
!
  if ( 1 <= element_show ) then
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to red.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw the triangles.'
    write ( file_unit, '(a)' ) '%'

    do triangle = 1, element_num

      write ( file_unit, '(a)' ) 'newpath'

      node = element_node(6,triangle)

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
        + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                   - y_min ) )

      write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' moveto'

      do i = 1, 3

        node = element_node(i,triangle)

        x_ps = int ( &
          ( ( x_max - node_xy(1,node)         ) &
          * real ( x_ps_min, kind = 8 )   &
          + (         node_xy(1,node) - x_min ) &
          * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                   - x_min ) )

        y_ps = int ( &
          ( ( y_max - node_xy(2,node)         ) &
          * real ( y_ps_min, kind = 8 )   &
          + (         node_xy(2,node) - y_min ) &
          * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                   - y_min ) )

        write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' lineto'

        node = element_node(i+3,triangle)

        x_ps = int ( &
          ( ( x_max - node_xy(1,node)         ) &
          * real ( x_ps_min, kind = 8 )   &
          + (         node_xy(1,node) - x_min ) &
          * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                   - x_min ) )

        y_ps = int ( &
          ( ( y_max - node_xy(2,node)         ) &
          * real ( y_ps_min, kind = 8 )   &
          + (         node_xy(2,node) - y_min ) &
          * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                   - y_min ) )

        write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' lineto'

      end do

      write ( file_unit, '(a)' ) 'stroke'

    end do

  end if
!
!  Label the triangles.
!
  if ( 2 <= element_show ) then

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Label the triangles:'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to darker red.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.950  0.250  0.150 setrgbcolor'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.20 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'

    do triangle = 1, element_num

      ave_x = 0.0D+00
      ave_y = 0.0D+00

      do i = 1, 6

        node = element_node(i,triangle)

        ave_x = ave_x + node_xy(1,node)
        ave_y = ave_y + node_xy(2,node)

      end do

      ave_x = ave_x / 6.0D+00
      ave_y = ave_y / 6.0D+00

      x_ps = int ( &
        ( ( x_max - ave_x         ) * real ( x_ps_min, kind = 8 )   &
        + (       + ave_x - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max         - x_min ) )

      y_ps = int ( &
        ( ( y_max - ave_y         ) * real ( y_ps_min, kind = 8 )   &
        + (         ave_y - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max         - y_min ) )

      write ( string, '(i4)' ) triangle
      string = adjustl ( string )

      write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps, ' moveto (' &
        // trim ( string ) // ') show'

    end do

  end if

  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'restore  showpage'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  End of page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%%Trailer'
  write ( file_unit, '(a)' ) '%%EOF'
  close ( unit = file_unit )

  return
end
subroutine triangulation_order6_print ( node_num, element_num, node_xy, &
  element_node, element_neighbor )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_PRINT prints out information defining a triangulation.
!
!  Discussion:
!
!    Triangulations created by R8TRIS2 include extra information encoded
!    in the negative values of ELEMENT_NEIGHBOR.
!
!    Because some of the nodes counted in NODE_NUM may not actually be
!    used in the triangulation, I needed to compute the true number
!    of vertices.  I added this calculation on 13 October 2001.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), the nodes
!    that make up the triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbors on each side.  If there is no triangle neighbor on
!    a particular side, the value of ELEMENT_NEIGHBOR should be negative.
!    If the triangulation data was created by R8TRIS2, then there is more
!    information encoded in the negative values.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6

  integer ( kind = 4 ) boundary_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) s
  logical              skip
  integer ( kind = 4 ) sp1
  integer ( kind = 4 ) t
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: vertex_list
  integer ( kind = 4 ) vertex_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGULATION_ORDER6_PRINT'
  write ( *, '(a)' ) '  Information defining an order6 triangulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of nodes is ', node_num

  call r8mat_transpose_print ( dim_num, node_num, node_xy, &
    '  Node coordinates' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of triangles is ', element_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sets of six nodes are used as vertices of'
  write ( *, '(a)' ) '  the triangles.  For each triangle, the vertices'
  write ( *, '(a)' ) '  are listed in counterclockwise order, followed'
  write ( *, '(a)' ) '  by the midside nodes.'

  call i4mat_transpose_print ( 6, element_num, element_node, &
    '  Triangle nodes:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  On each side of a given triangle, there is either'
  write ( *, '(a)' ) '  another triangle, or a piece of the convex hull.'
  write ( *, '(a)' ) '  For each triangle, we list the indices of the three'
  write ( *, '(a)' ) '  neighbors, or (if negative) the codes of the'
  write ( *, '(a)' ) '  segments of the convex hull.'

  call i4mat_transpose_print ( 3, element_num, element_neighbor, &
    '  Triangle neighbors' )
!
!  Determine the number of vertices.
!
  allocate ( vertex_list(1:3*element_num) )

  vertex_list(1:3*element_num) = reshape ( element_node(1:3,1:element_num), &
    (/ 3*element_num /) )

  call i4vec_sort_heap_a ( 3*element_num, vertex_list )

  call i4vec_sorted_unique ( 3*element_num, vertex_list, vertex_num )

  deallocate ( vertex_list )
!
!  Determine the number of boundary points.
!
  boundary_num = 2 * vertex_num - element_num - 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of boundary points is ', boundary_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The segments that make up the convex hull can be'
  write ( *, '(a)' ) '  determined from the negative entries of the triangle'
  write ( *, '(a)' ) '  neighbor list.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     #   Tri  Side    N1    N2    N3'
  write ( *, '(a)' ) ' '

  skip = .false.

  k = 0

  do i = 1, element_num

    do j = 1, 3

      if ( element_neighbor(j,i) < 0 ) then
        s = - element_neighbor(j,i)
        t = s / 3

        if ( t < 1 .or. element_num < t ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Sorry, this data does not use the R8TRIS2'
          write ( *, '(a)' ) '  convention for convex hull segments.'
          skip = .true.
          exit
        end if

        s = mod ( s, 3 ) + 1
        k = k + 1
        n1 = element_node(s,t)
        n2 = element_node(s+3,t)
        sp1 = s + 1
        sp1 = i4_wrap ( sp1, 1, 3 )
        n3 = element_node(sp1,t)
        write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4,2x,i4,2x,i4)' ) k, t, s, n1, n2, n3
      end if

    end do

    if ( skip ) then
      exit
    end if

  end do

  return
end
subroutine triangulation_order6_refine_compute ( node_num1, element_num1, &
  node_xy1, element_node1, node_num2, element_num2, edge_data, &
  node_xy2, element_node2 )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_REFINE_COMPUTE computes a refined order 6 triangulation.
!
!  Discussion:
!
!    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we
!    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
!    quadratic subtriangles T1, T2, T3 and T4.
!
!    The task is more complicated by the fact that we are working with
!    a mesh of triangles, so that we want to create a node only once,
!    even though it may be shared by other triangles.  (In fact, only
!    the new nodes on the edges can be shared, and then only by at most
!    one other triangle.)
!
!            3
!           / \
!          36 35
!         / T3  \
!        6--56---5
!       / \ T4  / \
!      16 46  45  25
!     / T1  \ / T2  \
!    1--14---4--24---2
!
!    This routine is given sorted information defining the edges, and uses
!    it to build the new node and triangle arrays.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM1, the number of triangles.
!
!    Input, real ( kind = 8 ) NODE_XY1(2,NODE_NUM1), the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE1(6,ELEMENT_NUM1), the nodes
!    that make up the triangles.  These should be listed in counterclockwise
!    order.
!
!    Input, integer ( kind = 4 ) NODE_NUM2, the number of nodes in the refined
!    mesh.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM2, the number of triangles in
!    the refined mesh.
!
!    Input, integer ( kind = 4 ) EDGE_DATA(5,3*ELEMENT_NUM1), edge data.
!
!    Output, real ( kind = 8 ) NODE_XY2(2,NODE_NUM2), the refined nodes.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE2(6,ELEMENT_NUM2), the nodes
!    that make up the triangles in the refined mesh.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) element_num1
  integer ( kind = 4 ) element_num2

  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(5,3*element_num1)
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) l3
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_old
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_old
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xy1(2,node_num1)
  real ( kind = 8 ) node_xy2(2,node_num2)
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2
  integer ( kind = 4 ) t3
  integer ( kind = 4 ) t4
  integer ( kind = 4 ) element_node1(6,element_num1)
  integer ( kind = 4 ) element_node2(6,element_num2)
  integer ( kind = 4 ) triangle1
  integer ( kind = 4 ) v1
  integer ( kind = 4 ) v2
  integer ( kind = 4 ) v3
  integer ( kind = 4 ) v4
  integer ( kind = 4 ) v5
  integer ( kind = 4 ) v6
!
!  Step 1:
!  Copy old nodes.
!
  node_xy2(1:2,1:node_num1) = node_xy1(1:2,1:node_num1)
!
!  Copy indices of existing nodes into new triangle array.
!
  element_node2(1:6,1:element_num2) = -1

  do triangle1 = 1, element_num1

    t1 = ( triangle1 - 1 ) * 4 + 1
    t2 = ( triangle1 - 1 ) * 4 + 2
    t3 = ( triangle1 - 1 ) * 4 + 3
    t4 = ( triangle1 - 1 ) * 4 + 4

    element_node2(1,t1) = element_node1(1,triangle1)
    element_node2(2,t1) = element_node1(4,triangle1)
    element_node2(3,t1) = element_node1(6,triangle1)

    element_node2(1,t2) = element_node1(4,triangle1)
    element_node2(2,t2) = element_node1(2,triangle1)
    element_node2(3,t2) = element_node1(5,triangle1)

    element_node2(1,t3) = element_node1(6,triangle1)
    element_node2(2,t3) = element_node1(5,triangle1)
    element_node2(3,t3) = element_node1(3,triangle1)

    element_node2(1,t4) = element_node1(5,triangle1)
    element_node2(2,t4) = element_node1(6,triangle1)
    element_node2(3,t4) = element_node1(4,triangle1)

  end do
!
!  Step 2.
!  Examine sorted edge information.  The first time an edge is encountered,
!  generate two new nodes, then assign them (usually) to the four subtriangles
!  of the two triangles that share that edge.
!
  node = node_num1

  n1_old = -1
  n2_old = -1

  do edge = 1, 3 * element_num1

    n1 = edge_data(1,edge)
    n2 = edge_data(2,edge)

    l1 = edge_data(3,edge)
    l3 = edge_data(4,edge)

    if ( l1 == 1 .and. l3 == 2 ) then
      l2 = 4
    else if ( l1 == 1 .and. l3 == 3 ) then
      l2 = 6
    else if ( l1 == 2 .and. l3 == 3 ) then
      l2 = 5
    end if

    triangle1 = edge_data(5,edge)
!
!  If this is the first time we've encountered this edge,
!  create the new nodes.
!
    if ( n1 /= n1_old .or. n2 /= n2_old ) then

      n1_old = n1
      n2_old = n2

      v1 = element_node1(l1,triangle1)
      v2 = element_node1(l2,triangle1)
      v3 = element_node1(l3,triangle1)

      node = node + 1
      v4 = node
      node_xy2(1:2,node) = 0.5D+00 * ( node_xy1(1:2,v1) + node_xy1(1:2,v2) )

      node = node + 1
      v5 = node
      node_xy2(1:2,node) = 0.5D+00 * ( node_xy1(1:2,v2) + node_xy1(1:2,v3) )

    end if

    t1 = ( triangle1 - 1 ) * 4 + 1
    t2 = ( triangle1 - 1 ) * 4 + 2
    t3 = ( triangle1 - 1 ) * 4 + 3

    if ( l1 == 1 .and. l3 == 2 ) then

      if ( element_node1(1,triangle1) == v1 ) then
        element_node2(4,t1) = v4
        element_node2(4,t2) = v5
      else
        element_node2(4,t1) = v5
        element_node2(4,t2) = v4
      end if

    else if ( l1 == 1 .and. l3 == 3 ) then

      if ( element_node1(1,triangle1) == v1 ) then
        element_node2(6,t1) = v4
        element_node2(6,t3) = v5
      else
        element_node2(6,t1) = v5
        element_node2(6,t3) = v4
      end if

    else if ( l1 == 2 .and. l3 == 3 ) then

      if ( element_node1(2,triangle1) == v1 ) then
        element_node2(5,t3) = v4
        element_node2(5,t2) = v5
      else
        element_node2(5,t3) = v5
        element_node2(5,t2) = v4
      end if

    end if

  end do
!
!  Step 3.
!  Each old triangle has a single central subtriangle, for which we now
!  need to generate three new "interior" nodes.
!
  do triangle1 = 1, element_num1

    v4 = element_node1(4,triangle1)
    v5 = element_node1(5,triangle1)
    v6 = element_node1(6,triangle1)

    t1 = ( triangle1 - 1 ) * 4 + 1
    t2 = ( triangle1 - 1 ) * 4 + 2
    t3 = ( triangle1 - 1 ) * 4 + 3
    t4 = ( triangle1 - 1 ) * 4 + 4

    node = node + 1
    node_xy2(1:2,node) = 0.5D+00 * ( node_xy1(1:2,v5) + node_xy1(1:2,v6) )
    element_node2(4,t4) = node
    element_node2(4,t3) = node

    node = node + 1
    node_xy2(1:2,node) = 0.5D+00 * ( node_xy1(1:2,v6) + node_xy1(1:2,v4) )
    element_node2(5,t4) = node
    element_node2(5,t1) = node

    node = node + 1
    node_xy2(1:2,node) = 0.5D+00 * ( node_xy1(1:2,v4) + node_xy1(1:2,v5) )
    element_node2(6,t4) = node
    element_node2(6,t2) = node

  end do

  return
end
subroutine triangulation_order6_refine_size ( node_num1, element_num1, &
  element_node1, node_num2, element_num2, edge_data )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_REFINE_SIZE sizes a refined order 6 triangulation.
!
!  Discussion:
!
!    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we
!    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
!    quadratic subtriangles T1, T2, T3 and T4.
!
!    The task is more complicated by the fact that we are working with
!    a mesh of triangles, so that we want to create a node only once,
!    even though it may be shared by other triangles.  (In fact, only
!    the new nodes on the edges can be shared, and then only by at most
!    one other triangle.)
!
!            3
!           / \
!          36 35
!         / T3  \
!        6--56---5
!       / \ T4  / \
!      16 46  45  25
!     / T1  \ / T2  \
!    1--14---4--24---2
!
!    This routine determines the sizes of the resulting node and
!    triangles, and constructs an edge array that can be used to
!    properly number the new nodes.
!
!    The primary work occurs in sorting a list related to the edges.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM1, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE1(6,ELEMENT_NUM1), the nodes
!    that make up the triangles.  These should be listed in counterclockwise
!    order.
!
!    Input, integer ( kind = 4 ) NODE_NUM2, the number of nodes in the refined
!    mesh.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM2, the number of triangles in
!    the refined mesh.
!
!    Output, integer ( kind = 4 ) EDGE_DATA(5,3*ELEMENT_NUM1), edge data
!    needed by TRIANGULATION_ORDER6_REFINE_COMPUTE.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) element_num1

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(5,3*element_num1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_old
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_old
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) element_num2
  integer ( kind = 4 ) element_node1(6,element_num1)
  integer ( kind = 4 ) triangle1
!
!  Step 1:
!  From the list of vertices for triangle T, of the form: (I,J,K),
!  construct the edge relations:
!
!    (I,J,1,2,T)
!    (I,K,1,3,T)
!    (J,K,2,3,T)
!
!  To make matching easier, we reorder each pair of nodes into
!  ascending order.
!
  do triangle1 = 1, element_num1

    i = element_node1(1,triangle1)
    j = element_node1(2,triangle1)
    k = element_node1(3,triangle1)

    a = min ( i, j )
    b = max ( i, j )

    edge_data(1:5,3*(triangle1-1)+1) = (/ a, b, 1, 2, triangle1 /)

    a = min ( i, k )
    b = max ( i, k )

    edge_data(1:5,3*(triangle1-1)+2) = (/ a, b, 1, 3, triangle1 /)

    a = min ( j, k )
    b = max ( j, k )

    edge_data(1:5,3*(triangle1-1)+3) = (/ a, b, 2, 3, triangle1 /)

  end do
!
!  Step 2: Perform an ascending dictionary sort on the relations.
!
  call i4col_sort_a ( 5, 3*element_num1, edge_data )
!
!  Step 3: Each shared edge will show up twice, consecutively,
!  in the EDGE_DATA array.  Each unique edge will generate
!  two new nodes, and each triangle will generate three new nodes.
!
  node_num2 = node_num1

  n1_old = -1
  n2_old = -1

  do edge = 1, 3 * element_num1
    n1 = edge_data(1,edge)
    n2 = edge_data(2,edge)
    if ( n1 /= n1_old .or. n2 /= n2_old ) then
      node_num2 = node_num2 + 2
      n1_old = n1
      n2_old = n2
    end if
  end do

  node_num2 = node_num2 + 3 * element_num1

  element_num2 = 4 * element_num1

  return
end
subroutine triangulation_order6_to_order3 ( element_num1, element_node1, &
  element_num2, element_node2 )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_TO_ORDER3 linearizes a quadratic triangulation.
!
!  Discussion:
!
!    A quadratic triangulation is assumed to consist of 6-node triangles,
!    as in the following:
!
!    11-12-13-14-15
!     |\    |\    |
!     | \   | \   |
!     6  7  8  9 10
!     |   \ |   \ |
!     |    \|    \|
!     1--2--3--4--5
!
!   This routine rearranges information so as to define the 3-node
!   triangulation:
!
!    11-12-13-14-15
!     |\ |\ |\ |\ |
!     | \| \| \| \|
!     6--7--8--9-10
!     |\ |\ |\ |\ |
!     | \| \| \| \|
!     1--2--3--4--5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM1, the number of triangles in
!    the quadratic triangulation.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE1(6,ELEMENT_NUM1), the quadratic
!    triangulation.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM2, the number of triangles in the
!    linear triangulation.  ELEMENT_NUM2 = 4 * ELEMENT_NUM1.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE2(3,ELEMENT_NUM2), the linear
!    triangulation.
!
  implicit none

  integer ( kind = 4 ) element_num1
  integer ( kind = 4 ) element_num2
  integer ( kind = 4 ), parameter :: element_order1 = 6
  integer ( kind = 4 ), parameter :: element_order2 = 3

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  integer ( kind = 4 ) tri1
  integer ( kind = 4 ) tri2
  integer ( kind = 4 ) element_node1(element_order1,element_num1)
  integer ( kind = 4 ) element_node2(element_order2,element_num2)

  tri2 = 0

  do tri1 = 1, element_num1

    n1 = element_node1(1,tri1)
    n2 = element_node1(2,tri1)
    n3 = element_node1(3,tri1)
    n4 = element_node1(4,tri1)
    n5 = element_node1(5,tri1)
    n6 = element_node1(6,tri1)

    tri2 = tri2 + 1
    element_node2(1:3,tri2) = (/ n1, n4, n6 /)
    tri2 = tri2 + 1
    element_node2(1:3,tri2) = (/ n2, n5, n4 /)
    tri2 = tri2 + 1
    element_node2(1:3,tri2) = (/ n3, n6, n5 /)
    tri2 = tri2 + 1
    element_node2(1:3,tri2) = (/ n4, n5, n6 /)

  end do

  return
end
subroutine triangulation_order6_vertex_count ( node_num, element_num, &
  element_node, vertex_num, midside_num )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_VERTEX_COUNT counts the vertex nodes in a triangulation.
!
!  Discussion:
!
!    In a triangulation of order 6, some nodes are midside nodes and some
!    nodes are vertex nodes.
!
!    Especially when an order 6 triangulation is used to handle the
!    Navier Stokes equations, it is useful to know the number of
!    vertex and midside nodes.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  6   5  side 2
!       |    \
!    3  |     \
!       |      \
!       1---4---2
!
!         side 1
!
!    The local node numbering.  Local nodes 1, 2 and 3 are "vertex nodes",
!    while nodes 4, 5 and 6 are "midside nodes".
!
!
!   21-22-23-24-25
!    |\    |\    |
!    | \   | \   |
!   16 17 18 19 20
!    |   \ |   \ |
!    |    \|    \|
!   11-12-13-14-15
!    |\    |\    |
!    | \   | \   |
!    6  7  8  9 10
!    |   \ |   \ |
!    |    \|    \|
!    1--2--3--4--5
!
!    A sample grid, which contains 9 vertex nodes and 16 midside nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), lists the
!    nodes that make up each triangle.  The first three nodes are the vertices,
!    in counterclockwise order.  The fourth value is the midside
!    node between nodes 1 and 2; the fifth and sixth values are
!    the other midside nodes in the logical order.
!
!    Output, integer ( kind = 4 ) VERTEX_NUM, the number of nodes which are
!    vertices.
!
!    Output, integer ( kind = 4 ) MIDSIDE_NUM, the number of nodes which are
!    midsides.  This value is inferred from NODE_NUM - VERTEX_NUM, so the value
!    of NODE_NUM needs to be correct on input!
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6

  integer ( kind = 4 ) midside_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) vertex_num
  integer ( kind = 4 ) vertices(3*element_num)

  vertices(               1:  element_num) = element_node(1,1:element_num)
  vertices(  element_num+1:2*element_num) = element_node(2,1:element_num)
  vertices(2*element_num+1:3*element_num) = element_node(3,1:element_num)

  call i4vec_sort_heap_a ( 3*element_num, vertices )

  call i4vec_sorted_unique ( 3*element_num, vertices, vertex_num )

  midside_num = node_num - vertex_num

  return
end
subroutine triangulation_search_delaunay ( node_num, node_xy, element_order, &
  element_num, element_node, element_neighbor, p, triangle_index, alpha, &
  beta, gamma, edge, step_num )

!*****************************************************************************80
!
!! TRIANGULATION_SEARCH_DELAUNAY searches a Delaunay triangulation for a point.
!
!  Discussion:
!
!    The algorithm "walks" from one triangle to its neighboring triangle,
!    and so on, until a triangle is found containing point P, or P is found
!    to be outside the convex hull.
!
!    The algorithm computes the barycentric coordinates of the point with
!    respect to the current triangle.  If all three quantities are positive,
!    the point is contained in the triangle.  If the I-th coordinate is
!    negative, then P lies on the far side of edge I, which is opposite
!    from vertex I.  This gives a hint as to where to search next.
!
!    For a Delaunay triangulation, the search is guaranteed to terminate.
!    For other triangulations, a cycle may occur.
!
!    Note the surprising fact that, even for a Delaunay triangulation of
!    a set of nodes, the nearest node to P need not be one of the
!    vertices of the triangle containing P.
!
!    The code can be called for triangulations of any order, but only
!    the first three nodes in each triangle are considered.  Thus, if
!    higher order triangles are used, and the extra nodes are intended
!    to give the triangle a polygonal shape, these will have no effect,
!    and the results obtained here might be misleading.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2012
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the nodes that make up each triangle.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbor list.
!
!    Input, real ( kind = 8 ) P(2), the coordinates of a point.
!
!    Output, integer ( kind = 4 ) TRIANGLE_INDEX, the index of the triangle
!    where the search ended.  If a cycle occurred, then TRIANGLE_INDEX = -1.
!
!    Output, real ( kind = 8 ) ALPHA, BETA, GAMMA, the barycentric
!    coordinates of the point relative to triangle TRIANGLE_INDEX.
!
!    Output, integer ( kind = 4 ) EDGE, indicates the position of the point P in
!    triangle TRIANGLE_INDEX:
!    0, the interior or boundary of the triangle;
!    -1, outside the convex hull of the triangulation, past edge 1;
!    -2, outside the convex hull of the triangulation, past edge 2;
!    -3, outside the convex hull of the triangulation, past edge 3.
!
!    Output, integer ( kind = 4 ) STEP_NUM, the number of steps.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) a
  real ( kind = 8 ) alpha
  integer ( kind = 4 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) c
  real ( kind = 8 ) det
  real ( kind = 8 ) dxp
  real ( kind = 8 ) dxa
  real ( kind = 8 ) dxb
  real ( kind = 8 ) dyp
  real ( kind = 8 ) dya
  real ( kind = 8 ) dyb
  integer ( kind = 4 ) edge
  real ( kind = 8 ) gamma
  real ( kind = 8 ) node_xy(dim_num,node_num)
  real ( kind = 8 ) p(dim_num)
  integer ( kind = 4 ) step_num
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) triangle_index
  integer ( kind = 4 ), save :: triangle_index_save = -1
  integer ( kind = 4 ) element_neighbor(3,element_num)
!
!  If possible, start with the previous successful value of TRIANGLE_INDEX.
!
  if ( triangle_index_save < 1 .or. element_num < triangle_index_save ) then
    triangle_index = ( element_num + 1 ) / 2
  else
    triangle_index = triangle_index_save
  end if

  step_num = - 1
  edge = 0

  do

    step_num = step_num + 1

    if ( element_num < step_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGULATION_SEARCH_DELAUNAY - Fatal error!'
      write ( *, '(a)' ) '  The algorithm seems to be cycling.'
      triangle_index = -1
      edge = -1
      stop
    end if
!
!  Get the nodes of triangle TRIANGLE_INDEX.
!
    a = element_node(1,triangle_index)
    b = element_node(2,triangle_index)
    c = element_node(3,triangle_index)
!
!  Using vertex C as a base, compute the distances to vertices A and B,
!  and the point P.
!
    dxa = node_xy(1,a) - node_xy(1,c)
    dya = node_xy(2,a) - node_xy(2,c)

    dxb = node_xy(1,b) - node_xy(1,c)
    dyb = node_xy(2,b) - node_xy(2,c)

    dxp = p(1)         - node_xy(1,c)
    dyp = p(2)         - node_xy(2,c)

    det = dxa * dyb - dya * dxb
!
!  Compute the barycentric coordinates of the point P with respect
!  to this triangle.
!
    alpha = ( dxp * dyb - dyp * dxb ) / det
    beta =  ( dxa * dyp - dya * dxp ) / det
    gamma = 1.0D+00 - alpha - beta
!
!  If the barycentric coordinates are all positive, then the point
!  is inside the triangle and we're done.
!
    if ( 0.0D+00 <= alpha .and. &
         0.0D+00 <= beta  .and. &
         0.0D+00 <= gamma ) then
      exit
    end if
!
!  At least one barycentric coordinate is negative.
!
!  If there is a negative barycentric coordinate for which there exists
!  an opposing triangle neighbor closer to the point, move to that triangle.
!
!  (Two coordinates could be negative, in which case we could go for the
!  most negative one, or the most negative one normalized by the actual
!  distance it represents).
!
    if ( alpha < 0.0D+00 .and. 0 < element_neighbor(2,triangle_index) ) then
      triangle_index = element_neighbor(2,triangle_index)
      cycle
    else if ( beta < 0.0D+00 .and. &
      0 < element_neighbor(3,triangle_index) ) then
      triangle_index = element_neighbor(3,triangle_index)
      cycle
    else if ( gamma < 0.0D+00 .and. &
      0 < element_neighbor(1,triangle_index) ) then
      triangle_index = element_neighbor(1,triangle_index)
      cycle
    end if
!
!  All negative barycentric coordinates correspond to vertices opposite
!  sides on the convex hull.
!
!  Note the edge and exit.
!
    if ( alpha < 0.0D+00 ) then
      edge = -2
      exit
    else if ( beta < 0.0D+00 ) then
      edge = -3
      exit
    else if ( gamma < 0.0D+00 ) then
      edge = -1
      exit
    end if

  end do

  triangle_index_save = triangle_index

  return
end
subroutine triangulation_search_naive ( node_num, node_xy, &
  element_order, element_num, element_node, p, triangle_index )

!*****************************************************************************80
!
!! TRIANGULATION_SEARCH_NAIVE naively searches a triangulation.
!
!  Discussion:
!
!    The algorithm simply checks each triangle to see if point P is
!    contained in it.  Surprisingly, this is not the fastest way to
!    do the check, at least if the triangulation is Delaunay.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2006
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
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles in
!    the triangulation.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the nodes that make up each triangle.
!
!    Input, real ( kind = 8 ) P(2), the coordinates of a point.
!
!    Output, integer ( kind = 4 ) TRIANGLE_INDEX, the index of the triangle
!    where the search ended, or -1 if no triangle was found containing
!    the point.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) a
  real ( kind = 8 ) alpha
  integer ( kind = 4 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) c
  real ( kind = 8 ) det
  real ( kind = 8 ) dxa
  real ( kind = 8 ) dxb
  real ( kind = 8 ) dxp
  real ( kind = 8 ) dya
  real ( kind = 8 ) dyb
  real ( kind = 8 ) dyp
  real ( kind = 8 ) gamma
  real ( kind = 8 ) node_xy(dim_num,node_num)
  real ( kind = 8 ) p(dim_num)
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) triangle_index

  triangle_index = -1

  do triangle = 1, element_num
!
!  Get the nodes of the triangle.
!
    a = element_node(1,triangle)
    b = element_node(2,triangle)
    c = element_node(3,triangle)
!
!  Using vertex C as a base, compute the distances to vertices A and B,
!  and the point P.
!
    dxa = node_xy(1,a) - node_xy(1,c)
    dya = node_xy(2,a) - node_xy(2,c)

    dxb = node_xy(1,b) - node_xy(1,c)
    dyb = node_xy(2,b) - node_xy(2,c)

    dxp = p(1)         - node_xy(1,c)
    dyp = p(2)         - node_xy(2,c)

    det = dxa * dyb - dya * dxb
!
!  Compute the barycentric coordinates of the point P with respect
!  to this triangle.
!
    alpha = ( dxp * dyb - dyp * dxb ) / det
    beta  = ( dxa * dyp - dya * dxp ) / det
    gamma = 1.0D+00 - alpha - beta
!
!  If the barycentric coordinates are all positive, then the point
!  is inside the triangle and we're done.
!
    if ( 0.0D+00 <= alpha .and. &
         0.0D+00 <= beta  .and. &
         0.0D+00 <= gamma ) then
      triangle_index = triangle
      return
    end if

  end do

  return
end
subroutine vbedg ( x, y, node_num, node_xy, element_num, element_node, &
  element_neighbor, ltri, ledg, rtri, redg )

!*****************************************************************************80
!
!! VBEDG determines which boundary edges are visible to a point.
!
!  Discussion:
!
!    The point (X,Y) is assumed to be outside the convex hull of the
!    region covered by the 2D triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of a point outside the
!    convex hull of the current triangulation.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the triangle
!    incidence list.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbor list; negative values are used for links of a
!    counter clockwise linked list of boundary edges;
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input/output, integer ( kind = 4 ) LTRI, LEDG.  If LTRI /= 0 then these
!    values are assumed to be already computed and are not changed, else they
!    are updated.  On output, LTRI is the index of boundary triangle to the
!    left of the leftmost boundary triangle visible from (X,Y), and LEDG is
!    the boundary edge of triangle LTRI to the left of the leftmost boundary
!    edge visible from (X,Y).  1 <= LEDG <= 3.
!
!    Input/output, integer ( kind = 4 ) RTRI.  On input, the index of the
!    boundary triangle to begin the search at.  On output, the index of the
!    rightmost boundary triangle visible from (X,Y).
!
!    Input/output, integer ( kind = 4 ) REDG, the edge of triangle RTRI that
!    is visible from (X,Y).  1 <= REDG <= 3.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  logical              ldone
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  real ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) t
  integer ( kind = 4 ) element_node(3,element_num)
  integer ( kind = 4 ) element_neighbor(3,element_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Find the rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor information.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

  do

    l = -element_neighbor(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = element_node(e,t)

    if ( e <= 2 ) then
      b = element_node(e+1,t)
    else
      b = element_node(1,t)
    end if

    lr = lrline ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,b), &
      node_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

    rtri = t
    redg = e

  end do

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

  do

    b = element_node(e,t)
    e = e - 1
    e = i4_wrap ( e, 1, 3 )

    do while ( 0 < element_neighbor(e,t) )

      t = element_neighbor(e,t)

      if ( element_node(1,t) == b ) then
        e = 3
      else if ( element_node(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = element_node(e,t)

    lr = lrline ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,b), &
      node_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end
subroutine voronoi_polygon_area ( node, neighbor_num, neighbor_index, &
  node_num, node_xy, area )

!*****************************************************************************80
!
!! VORONOI_POLYGON_AREA computes the area of a Voronoi polygon.
!
!  Discussion:
!
!    It is assumed that the Voronoi polygon is finite!  Every Voronoi
!    diagram includes some regions which are infinite, and for those,
!    this formula is not appropriate.
!
!    The routine is given the indices of the nodes that are neighbors of a
!    given "center" node.  A node is a neighbor of the center node if the
!    Voronoi polygons of the two nodes share an edge.  The triangles of the
!    Delaunay triangulation are formed from successive pairs of these neighbor
!    nodes along with the center node.
!
!    The assumption that the polygon is a Voronoi polygon is
!    used to determine the location of the boundaries of the polygon,
!    which are the perpendicular bisectors of the lines connecting
!    the center point to each of its neighbors.
!
!    The finiteness assumption is employed in part in the
!    assumption that the polygon is bounded by the finite
!    line segments from point 1 to 2, 2 to 3, ...,
!    M-1 to M, and M to 1, where M is the number of neighbors.
!
!    It is assumed that this subroutine is being called by a
!    process which has computed the Voronoi diagram of a large
!    set of nodes, so the arrays X and Y are dimensioned by
!    NODE_NUM, which may be much greater than the number of neighbor
!    nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
!    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
!    Second Edition,
!    Wiley, 2000, page 485.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the node whose Voronoi
!    polygon is to be measured.  1 <= NODE <= NODE_NUM.
!
!    Input, integer ( kind = 4 ) NEIGHBOR_NUM, the number of neighbor nodes of
!    the given node.
!
!    Input, integer ( kind = 4 ) NEIGHBOR_INDEX(NEIGHBOR_NUM), the indices
!    of the neighbor nodes (used to access X and Y).  The neighbor
!    nodes should be listed in the (counter-clockwise) order in
!    which they occur as one circles the center node.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Output, real ( kind = 8 ) AREA, the area of the Voronoi polygon.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) neighbor_num
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) a
  real ( kind = 8 ) area
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) neighbor_index(neighbor_num)
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xy(dim_num,node_num)
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) pi(dim_num)
  real ( kind = 8 ) pj(dim_num)
  real ( kind = 8 ) ui(dim_num)
  real ( kind = 8 ) uj(dim_num)

  area = 0.0D+00

  if ( node < 1 .or. node_num < node ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_POLYGON_AREA - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of input parameter NODE.'
    stop
  end if

  pc(1:dim_num) = node_xy(1:dim_num,node)

  i = neighbor_num

  pi(1:dim_num) = node_xy(1:dim_num,i)

  j = 1
  j = neighbor_index(j)

  pj(1:dim_num) = node_xy(1:dim_num,j)

  a = ( pi(1)**2 + pi(2)**2 - pc(1)**2 - pc(2)**2 )
  b = ( pj(1)**2 + pj(2)**2 - pc(1)**2 - pc(2)**2 )
  c = 2.0D+00 * ( ( pi(1) - pc(1) ) * ( pj(2) - pc(2) ) &
                - ( pj(1) - pc(1) ) * ( pi(2) - pc(2) ) )
  uj(1) = ( a * ( pj(2) - pc(2) ) - b * ( pi(2) - pc(2) )  ) / c
  uj(2) = ( a * ( pj(1) - pc(1) ) - b * ( pi(1) - pc(1) )  ) / c

  do i = 1, neighbor_num

    pi(1:dim_num) = pj(1:dim_num)

    ui(1:dim_num) = uj(1:dim_num)

    j = i + 1
    if ( neighbor_num < j ) then
      j = 1
    end if

    j = neighbor_index(j)

    pj(1:dim_num) = node_xy(1:dim_num,j)

    a = ( pi(1)**2 + pi(2)**2 - pc(1)**2 - pc(2)**2 )
    b = ( pj(1)**2 + pj(2)**2 - pc(1)**2 - pc(2)**2 )
    c = 2.0D+00 * ( ( pi(1) - pc(1) ) * ( pj(2) - pc(2) ) &
                  - ( pj(1) - pc(1) ) * ( pi(2) - pc(2) ) )
    uj(1) = ( a * ( pj(2) - pc(2) ) - b * ( pi(2) - pc(2) )  ) / c
    uj(2) = ( a * ( pj(1) - pc(1) ) - b * ( pi(1) - pc(1) )  ) / c

    area = area + uj(1) * ui(2) - ui(1) * uj(2)

  end do

  area = 0.5D+00 * area

  return
end
subroutine voronoi_polygon_centroid ( node, neighbor_num, neighbor_index, &
  node_num, node_xy, centroid )

!*****************************************************************************80
!
!! VORONOI_POLYGON_CENTROID computes the centroid of a Voronoi polygon.
!
!  Discussion:
!
!    It is assumed that the Voronoi polygon is finite!  Every Voronoi
!    diagram includes some regions which are infinite, and for those,
!    this formula is not appropriate.
!
!    The routine is given the indices of the nodes that are neighbors of a
!    given "center" node.  A node is a neighbor of the center node if the
!    Voronoi polygons of the two nodes share an edge.  The triangles of the
!    Delaunay triangulation are formed from successive pairs of these neighbor
!    nodes along with the center node.
!
!    The assumption that the polygon is a Voronoi polygon is
!    used to determine the location of the boundaries of the polygon,
!    which are the perpendicular bisectors of the lines connecting
!    the center point to each of its neighbors.
!
!    The finiteness assumption is employed in part in the
!    assumption that the polygon is bounded by the finite
!    line segments from point 1 to 2, 2 to 3, ...,
!    M-1 to M, and M to 1, where M is the number of neighbors.
!
!    It is assumed that this subroutine is being called by a
!    process which has computed the Voronoi diagram of a large
!    set of nodes, so the arrays X and Y are dimensioned by
!    NODE_NUM, which may be much greater than the number of neighbor
!    nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
!    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
!    Second Edition,
!    Wiley, 2000, page 490.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the node whose Voronoi
!    polygon is to be analyzed.  1 <= NODE <= NODE_NUM.
!
!    Input, integer ( kind = 4 ) NEIGHBOR_NUM, the number of neighbor nodes of
!    the given node.
!
!    Input, integer ( kind = 4 ) NEIGHBOR_INDEX(NEIGHBOR_NUM), the indices
!    of the neighbor nodes.  These indices are used to access the
!    X and Y arrays.  The neighbor nodes should be listed in the
!    (counter-clockwise) order in which they occur as one circles
!    the center node.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the total number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Output, real ( kind = 8 ) CENTROID(2), the coordinates of the centroid
!    of the Voronoi polygon of node NODE.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) neighbor_num
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) a
  real ( kind = 8 ) area
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) centroid(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) neighbor_index(neighbor_num)
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xy(dim_num,node_num)
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) pi(dim_num)
  real ( kind = 8 ) pj(dim_num)
  real ( kind = 8 ) ui(dim_num)
  real ( kind = 8 ) uj(dim_num)

  centroid(1:dim_num) = 0.0D+00

  if ( node < 1 .or. node_num < node ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_POLYGON_CENTROID - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of input parameter NODE.'
    stop
  end if

  pc(1:dim_num) = node_xy(1:dim_num,node)

  i = neighbor_num
  i = neighbor_index(i)

  pi(1:dim_num) = node_xy(1:dim_num,i)

  j = 1
  j = neighbor_index(j)

  pj(1:dim_num) = node_xy(1:dim_num,j)

  a = ( pi(1) * pi(1) + pi(2) * pi(2) - pc(1) * pc(1) - pc(2) * pc(2) )
  b = ( pj(1) * pj(1) + pj(2) * pj(2) - pc(1) * pc(1) - pc(2) * pc(2) )
  c = 2.0D+00 * ( ( pi(1) - pc(1) ) * ( pj(2) - pc(2) ) &
                - ( pj(1) - pc(1) ) * ( pi(2) - pc(2) ) )
  uj(1) = ( a * ( pj(2) - pc(2) ) - b * ( pi(2) - pc(2) )  ) / c
  uj(2) = ( a * ( pj(1) - pc(1) ) - b * ( pi(1) - pc(1) )  ) / c

  do i = 1, neighbor_num

    pi(1:dim_num) = pj(1:dim_num)
    ui(1:dim_num) = uj(1:dim_num)

    j = i + 1
    if ( neighbor_num < j ) then
      j = 1
    end if

    pj(1:dim_num) = node_xy(1:dim_num,j)

    a = ( pi(1) * pi(1) + pi(2) * pi(2) - pc(1) * pc(1) - pc(2) * pc(2) )
    b = ( pj(1) * pj(1) + pj(2) * pj(2) - pc(1) * pc(1) - pc(2) * pc(2) )
    c = 2.0D+00 * ( ( pi(1) - pc(1) ) * ( pj(2) - pc(2) ) &
                  - ( pj(1) - pc(1) ) * ( pi(2) - pc(2) ) )
    uj(1) = ( a * ( pj(2) - pc(2) ) - b * ( pi(2) - pc(2) )  ) / c
    uj(2) = ( a * ( pj(1) - pc(1) ) - b * ( pi(1) - pc(1) )  ) / c

    centroid(1) = centroid(1) + ( ui(2) - uj(2) ) &
      * ( ( uj(1) + ui(1) )**2 - uj(1) * ui(1) )
    centroid(2) = centroid(2) + ( ui(1) - uj(1) ) &
      * ( ( uj(2) + ui(2) )**2 - uj(2) * ui(2) )

  end do

  call voronoi_polygon_area ( node, neighbor_num, neighbor_index, &
    node_num, node_xy, area )

  centroid(1:dim_num) = centroid(1:dim_num) / ( 6.0D+00 * area )

  return
end
subroutine voronoi_polygon_vertices ( node, neighbor_num, neighbor_index, &
  node_num, node_xy, v )

!*****************************************************************************80
!
!! VORONOI_POLYGON_VERTICES computes the vertices of a Voronoi polygon.
!
!  Discussion:
!
!    This routine is only appropriate for Voronoi polygons that are finite.
!
!    The routine is given the indices of the nodes that are neighbors of a
!    given "center" node.  A node is a neighbor of the center node if the
!    Voronoi polygons of the two nodes share an edge.  The triangles of the
!    Delaunay triangulation are formed from successive pairs of these neighbor
!    nodes along with the center node.
!
!    Given only the neighbor node information, it is possible to determine
!    the location of the vertices of the polygonal Voronoi region by computing
!    the circumcenters of the Delaunay triangles.
!
!    It is assumed that this subroutine is being called by a process which has
!    computed the Voronoi diagram of a large set of nodes, so the arrays X and
!    Y are dimensioned by NODE_NUM, which may be much greater than the number
!    of neighbor nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
!    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
!    Second Edition,
!    Wiley, 2000.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the node whose Voronoi
!    polygon is to be analyzed.  1 <= NODE <= NODE_NUM.
!
!    Input, integer ( kind = 4 ) NEIGHBOR_NUM, the number of neighbor nodes of
!    the given node.
!
!    Input, integer ( kind = 4 ) NEIGHBOR_INDEX(NEIGHBOR_NUM), the indices
!    of the neighbor nodes.  These indices are used to access the
!    X and Y arrays.  The neighbor nodes should be listed in the
!    (counter-clockwise) order in which they occur as one circles
!    the center node.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Output, real ( kind = 8 ) V(2,NEIGHBOR_NUM), the coordinates of
!    the vertices of the Voronoi polygon around node NODE.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) neighbor_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) neighbor_index(neighbor_num)
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xy(dim_num,node_num)
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ) v(dim_num,neighbor_num)

  if ( node < 1 .or. node_num < node ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_POLYGON_VERTICES - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of input parameter NODE.'
    stop
  end if

  t(1:dim_num,1) = node_xy(1:dim_num,node)

  ip1 = neighbor_index(1)
  t(1:dim_num,3) = node_xy(1:dim_num,ip1)

  do i = 1, neighbor_num

    t(1:dim_num,2) = t(1:dim_num,3)

    ip1 = i + 1
    if ( neighbor_num < ip1 ) then
      ip1 = 1
    end if

    ip1 = neighbor_index(ip1)
    t(1:dim_num,3) = node_xy(1:dim_num,ip1)

    call triangle_circumcenter_2d ( t, v(1:dim_num,i) )

  end do

  return
end
