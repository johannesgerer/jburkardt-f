function alpha_measure ( n, z, triangle_order, triangle_num, triangle_node )

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
!    07 November 2005
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
!    Input, integer ( kind = 4 ) TRIANGLE_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM),
!    the triangulation.
!
!    Output, real ( kind = 8 ) ALPHA_MEASURE, the ALPHA quality measure.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) triangle_order

  real ( kind = 8 ) a_angle
  integer ( kind = 4 ) a_index
  real ( kind = 8 ) a_x
  real ( kind = 8 ) a_y
  real ( kind = 8 ) ab_len
  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha_measure
  real ( kind = 8 ) arc_cosine
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
  integer ( kind = 4 ) triangle_node(triangle_order,triangle_num)
  real ( kind = 8 ) z(2,n)

  alpha = huge ( alpha )

  do triangle = 1, triangle_num

    a_index = triangle_node(1,triangle)
    b_index = triangle_node(2,triangle)
    c_index = triangle_node(3,triangle)

    a_x = z(1,a_index)
    a_y = z(2,a_index)
    b_x = z(1,b_index)
    b_y = z(2,b_index)
    c_x = z(1,c_index)
    c_y = z(2,c_index)

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

    alpha = min ( alpha, a_angle )
    alpha = min ( alpha, b_angle )
    alpha = min ( alpha, c_angle )

  end do
!
!  Normalize angles from [0,60] degrees into qualities in [0,1].
!
  alpha_measure = alpha * 3.0D+00 / pi

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
function area_measure ( n, z, triangle_order, triangle_num, triangle_node )

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
!    07 November 2005
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
!    Input, integer ( kind = 4 ) TRIANGLE_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM),
!    the triangulation.
!
!    Output, real ( kind = 8 ) AREA_MEASURE, the AREA quality measure.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) triangle_order

  real ( kind = 8 ) area
  real ( kind = 8 ) area_max
  real ( kind = 8 ) area_measure
  real ( kind = 8 ) area_min
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle_node(triangle_order,triangle_num)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) z(2,n)

  area_max = 0.0D+00
  area_min = huge ( area_min )

  do triangle = 1, triangle_num

    x1 = z(1,triangle_node(1,triangle))
    y1 = z(2,triangle_node(1,triangle))
    x2 = z(1,triangle_node(2,triangle))
    y2 = z(2,triangle_node(2,triangle))
    x3 = z(1,triangle_node(3,triangle))
    y3 = z(2,triangle_node(3,triangle))

    area = 0.5D+00 * abs ( x1 * ( y2 - y3 ) &
                         + x2 * ( y3 - y1 ) &
                         + x3 * ( y1 - y2 ) )

    area_min = min ( area_min, area )
    area_max = max ( area_max, area )

  end do

  if ( 0.0D+00 < area_max ) then
    area_measure = area_min / area_max
  else
    area_measure = 0.0D+00
  end if

  return
end
subroutine bandwidth_mesh ( element_order, element_num, element_node, &
  ml, mu, m )

!*****************************************************************************80
!
!! BANDWIDTH_MESH: bandwidth of finite element mesh.
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
!    Output, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths of the matrix.
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
function beta_measure ( dim_num, n, z )

!*****************************************************************************80
!
!! BETA_MEASURE determines the pointset quality measure BETA.
!
!  Discussion:
!
!    The BETA measure of point distribution quality for a set Z of
!    N points in an DIM_NUM dimensional region is defined as follows:
!
!    For each point Z(I), determine the nearest distinct element of
!    the pointset by
!
!      GAMMA(I) = minimum ( 1 <= J <= N, I /= J ) distance ( Z(I), Z(J) )
!
!    Let GAMMA_AVE be the average of GAMMA(1:N).
!
!    Let GAMMA_STD be the standard deviation of the GAMMA's:
!
!      GAMMA_STD = sqrt ( 1 / ( N - 1 )
!        * sum ( 1 <= I <= N ) ( GAMMA(I) - GAMMA_AVE )**2 ) )
!
!    Then BETA is the standard deviation normalized by the average:
!
!      BETA = GAMMA_STD / GAMMA_AVE.
!
!    For an ideally regular mesh, the GAMMA(I)'s will be equal and
!    BETA will be zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 July 2005
!
!  Author:
!
!    Max Gunzburger
!    John Burkardt
!
!  Reference:
!
!    Max Gunzburger and John Burkardt,
!    Uniformity Measures for Point Samples in Hypercubes.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the points.
!
!    Output, real ( kind = 8 ) BETA_MEASURE, the BETA quality measure.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) beta_measure
  real ( kind = 8 ) gamma(n)
  real ( kind = 8 ) gamma_ave
  real ( kind = 8 ) gamma_std
  real ( kind = 8 ) z(dim_num,n)

  call pointset_spacing ( dim_num, n, z, gamma )

  gamma_ave = sum ( gamma(1:n) ) / real ( n, kind = 8 )

  if ( 1 < n ) then
    gamma_std = sqrt ( sum ( ( gamma(1:n) - gamma_ave )**2 ) &
      / real ( n - 1, kind = 8 ) )
  else
    gamma_std = 0.0D+00
  end if

  if ( 0.0D+00 < gamma_ave ) then
    beta_measure = gamma_std / gamma_ave
  else
    beta_measure = 0.0D+00
  end if

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
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer ( kind = 4 ) value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
function chi_measure ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! CHI_MEASURE determines the pointset quality measure CHI.
!
!  Discussion:
!
!    The CHI measure of point distribution quality for a set Z of
!    N points in an DIM_NUM dimensional region is defined as follows:
!
!    Assign every point X in the region to the nearest element
!    Z(I) of the point set.  For each Z(I), let H(I) be the maximum
!    distance between Z(I) and any point assigned to it by this process.
!
!    For each point Z(I), we determine the nearest distinct element of
!    the pointset by
!
!      GAMMA(I) = minimum ( 1 <= J <= N, I /= J ) distance ( Z(I), Z(J) )
!
!    Then
!
!      CHI(I) = 2 * H(I) / GAMMA(I)
!
!    and
!
!      CHI = maximum ( 1 <= I <= N ) CHI(I)
!
!    This quantity can be estimated by using sampling to pick a large
!    number of points in the region, rather than all of them.
!
!    For an ideally regular mesh, all the CHI(I)'s will be equal.
!    Any deviation from regularity increases the value of some entries
!    of CHI; thus, given two meshes, the one with a lower value of
!    CHI is to be recommended.
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
!  Reference:
!
!    Max Gunzburger and John Burkardt,
!    Uniformity Measures for Point Samples in Hypercubes.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the points.
!
!    Input, integer ( kind = 4 ) NS, the number of sample points.
!
!    Input, external SAMPLE_ROUTINE, the name of a routine which
!    is used to produce sample points in the region, of the form:
!      subroutine sample_routine ( dim_num, n, seed, x )
!      integer ( kind = 4 ) dim_num
!      integer ( kind = 4 ) n
!      integer ( kind = 4 ) seed
!      real ( kind = 8 ) x(dim_num,n)
!
!    Input, integer ( kind = 4 ) SEED_INIT, the initial value of the random
!    number seed.
!
!    Output, real ( kind = 8 ) CHI_MEASURE, the CHI quality measure.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) chi(n)
  real ( kind = 8 ) chi_measure
  integer ( kind = 4 ) closest(1)
  real ( kind = 8 ) dist
  real ( kind = 8 ) gamma(n)
  real ( kind = 8 ) h(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ns
  external sample_routine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) z(dim_num,n)

  seed = seed_init
  h(1:n) = 0.0D+00

  do k = 1, ns

    call sample_routine ( dim_num, 1, seed, x )

    call find_closest ( dim_num, n, 1, x, z, closest )

    dist = sum ( ( x(1:dim_num) - z(1:dim_num,closest(1) ) )**2 )

    h(closest(1)) = max ( h(closest(1)), dist )

  end do

  call pointset_spacing ( dim_num, n, z, gamma )
!
!  If any GAMMA is 0, then two points are identical.
!
  if ( any ( gamma(1:n) == 0.0D+00 ) ) then
    chi_measure = huge ( gamma(1) )
    return
  end if

  chi(1:n) = 2.0D+00 * sqrt ( h(1:n) ) / gamma(1:n)

  chi_measure = maxval ( chi(1:n) )

  return
end
function d_measure ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! D_MEASURE determines the pointset quality measure D.
!
!  Discussion:
!
!    The D measure of point distribution quality for a set Z of
!    N points in an DIM_NUM-dimensional region is defined as follows:
!
!    For each point Z(I) in the pointset, let V(I) be the region
!    defined by the intersection of the region with the Voronoi
!    subregion associated with Z(I).
!
!    Let D(I) be the determinant of the deviatoric tensor associated with
!    the region V(I).
!
!    Then D = maximum ( 1 <= I <= N ) D(I).
!
!    This quantity can be estimated using sampling.  A given number of
!    sample points are generated in the region, assigned to the nearest
!    element of the pointset, and used to approximate the Voronoi subregions
!    and the deviatoric tensors.
!
!    In an ideally regular mesh, each deviatoric tensor would have a
!    zero determinant, and hence D would be zero.  In general, the smaller
!    D, the better.
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
!    Max Gunzburger and John Burkardt,
!    Uniformity Measures for Point Samples in Hypercubes.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the points.
!
!    Input, integer ( kind = 4 ) NS, the number of sample points.
!
!    Input, external SAMPLE_ROUTINE, the name of a routine which
!    is used to produce sample points in the region, of the form:
!      subroutine sample_routine ( dim_num, n, seed, x )
!      integer ( kind = 4 ) dim_num
!      integer ( kind = 4 ) n
!      integer ( kind = 4 ) seed
!      real ( kind = 8 ) x(dim_num,n)
!
!    Input, integer ( kind = 4 ) SEED_INIT, the initial value of the random
!    number seed.
!
!    Output, real ( kind = 8 ) D_MEASURE, the D quality measure.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) centroid(dim_num,n)
  integer ( kind = 4 ) closest(1)
  real ( kind = 8 ) d_measure
  real ( kind = 8 ) dge_det
  real ( kind = 8 ) di
  integer ( kind = 4 ) hit(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) moment(dim_num,dim_num,n)
  integer ( kind = 4 ) ns
  external sample_routine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) tri(n)
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) z(dim_num,n)

  seed = seed_init
  centroid(1:dim_num,1:n) = 0.0D+00
  hit(1:n) = 0
  moment(1:dim_num,1:dim_num,1:n) = 0.0D+00

  do k = 1, ns

    call sample_routine ( dim_num, 1, seed, x )

    call find_closest ( dim_num, n, 1, x, z, closest )

    hit(closest(1)) = hit(closest(1)) + 1

    centroid(1:dim_num,closest(1)) = centroid(1:dim_num,closest(1)) &
      + x(1:dim_num)

    do i1 = 1, dim_num
      do i2 = 1, dim_num
        moment(i1,i2,closest(1)) = moment(i1,i2,closest(1)) + x(i1) * x(i2)
      end do
    end do

  end do

  do j = 1, n

    if ( 0 < hit(j) ) then

      centroid(1:dim_num,j) = centroid(1:dim_num,j) / real ( hit(j), kind = 8 )

      moment(1:dim_num,1:dim_num,j) = moment(1:dim_num,1:dim_num,j) &
        / real ( hit(j), kind = 8 )

      do i1 = 1, dim_num
        do i2 = 1, dim_num
          moment(i1,i2,j) = moment(i1,i2,j) - centroid(i1,j) * centroid(i2,j)
        end do
      end do

    end if

  end do

  tri(1:n) = 0.0D+00

  do j = 1, n
    do i = 1, dim_num
      tri(j) = tri(j) + moment(i,i,j)
    end do
  end do

  do j = 1, n
    do i = 1, dim_num
      moment(i,i,j) = moment(i,i,j) - tri(j) / real ( dim_num, kind = 8 )
    end do
  end do

  d_measure = 0.0D+00

  do j = 1, n

    di = dge_det ( dim_num, moment(1:dim_num,1:dim_num,j) )

    d_measure = max ( d_measure, di )

  end do

  return
end
function dge_det ( n, a )

!*****************************************************************************80
!
!! DGE_DET computes the determinant of a square matrix in DGE storage.
!
!  Discussion:
!
!    The DGE storage format is used for a general M by N matrix.  A storage
!    space is made for each logical entry.  The two dimensional logical
!    array is mapped to a vector, in which storage is by columns.
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
!    Dongarra, Bunch, Moler, Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N,N), the matrix to be analyzed.
!    On output, the matrix has been overwritten by factorization information.
!
!    Output, real ( kind = 8 ) DGE_DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) dge_det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t

  dge_det = 1.0D+00

  do k = 1, n-1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k+1, n
      if ( abs ( a(l,k) ) < abs ( a(i,k) ) ) then
        l = i
      end if
    end do

    dge_det = dge_det * a(l,k)
!
!  If the pivot value is zero, then the determinant is zero,
!  and we are done.  In fact, we do NOT want to continue, because
!  the scaling below would fail.
!
    if ( a(l,k) == 0.0D+00 ) then
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      t      = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = -a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k+1, n

      if ( l /= k ) then
        t      = a(l,j)
        a(l,j) = a(k,j)
        a(k,j) = t
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  dge_det = dge_det * a(n,n)

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
!    counter clockwise order.
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

    diaedg = -1

  else if ( ca < -tola .and. cb < -tolb ) then

    diaedg = 1

  else

    tola = max ( tola, tolb )
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( tola < s ) then
      diaedg = -1
    else if ( s < -tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

  return
end
subroutine dtris2 ( node_num, node_xy, triangle_num, triangle_node, &
  triangle_neighbor )

!*****************************************************************************80
!
!! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
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
!    Output, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles in
!    the triangulation; TRIANGLE_NUM is equal to 2*NODE_NUM - NB - 2, where
!    NB is the number of boundary vertices.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes
!    that make up each triangle.  The elements are indices of P.  The vertices
!    of the triangles are in counter clockwise order.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the
!    triangle neighbor list.  Positive elements are indices of TIL; negative
!    elements are used for links of a counter clockwise linked list of boundary
!    edges; LINK = -(3*I + J-1) where I, J = triangle, edge index;
!    TRIANGLE_NEIGHBOR(J,I) refers to the neighbor along edge from vertex J
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
  integer ( kind = 4 ) triangle_neighbor(3,node_num*2)
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) triangle_node(3,node_num*2)

  tol = 100.0D+00 * epsilon ( tol )

  ierr = 0
!
!  Sort the vertices by increasing (x,y).
!
  call r82vec_sort_heap_index_a ( node_num, node_xy, indx )

  call r82vec_permute ( node_num, node_xy, indx )
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
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
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
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
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
  triangle_num = j - 2

  if ( lr == -1 ) then

    triangle_node(1,1) = m1
    triangle_node(2,1) = m2
    triangle_node(3,1) = m
    triangle_neighbor(3,1) = -3

    do i = 2, triangle_num

      m1 = m2
      m2 = i+1

      triangle_node(1,i) = m1
      triangle_node(2,i) = m2
      triangle_node(3,i) = m

      triangle_neighbor(1,i-1) = -3 * i
      triangle_neighbor(2,i-1) = i
      triangle_neighbor(3,i) = i - 1

    end do

    triangle_neighbor(1,triangle_num) = -3 * triangle_num - 1
    triangle_neighbor(2,triangle_num) = -5
    ledg = 2
    ltri = triangle_num

  else

    triangle_node(1,1) = m2
    triangle_node(2,1) = m1
    triangle_node(3,1) = m

    triangle_neighbor(1,1) = -4

    do i = 2, triangle_num

      m1 = m2
      m2 = i+1

      triangle_node(1,i) = m2
      triangle_node(2,i) = m1
      triangle_node(3,i) = m

      triangle_neighbor(3,i-1) = i
      triangle_neighbor(1,i) = -3 * i - 3
      triangle_neighbor(2,i) = i - 1

    end do

    triangle_neighbor(3,triangle_num) = -3 * triangle_num
    triangle_neighbor(2,1) = -3 * triangle_num - 2
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
    m1 = triangle_node(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = triangle_node(ledg+1,ltri)
    else
      m2 = triangle_node(1,ltri)
    end if

    lr = lrline ( node_xy(1,m), node_xy(2,m), node_xy(1,m1), &
      node_xy(2,m1), node_xy(1,m2), node_xy(2,m2), 0.0D+00 )

    if ( 0 < lr ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -triangle_neighbor(ledg,ltri)
      rtri = l / 3
      redg = mod ( l, 3 ) + 1
    end if

    call vbedg ( node_xy(1,m), node_xy(2,m), node_num, node_xy, triangle_num, &
      triangle_node, triangle_neighbor, ltri, ledg, rtri, redg )

    n = triangle_num + 1
    l = -triangle_neighbor(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -triangle_neighbor(e,t)
      m2 = triangle_node(e,t)

      if ( e <= 2 ) then
        m1 = triangle_node(e+1,t)
      else
        m1 = triangle_node(1,t)
      end if

      triangle_num = triangle_num + 1
      triangle_neighbor(e,t) = triangle_num

      triangle_node(1,triangle_num) = m1
      triangle_node(2,triangle_num) = m2
      triangle_node(3,triangle_num) = m

      triangle_neighbor(1,triangle_num) = t
      triangle_neighbor(2,triangle_num) = triangle_num - 1
      triangle_neighbor(3,triangle_num) = triangle_num + 1

      top = top + 1

      if ( node_num < top ) then
        ierr = 8
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
        write ( *, '(a)' ) '  Stack overflow.'
        stop
      end if

      stack(top) = triangle_num

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    triangle_neighbor(ledg,ltri) = -3 * n - 1
    triangle_neighbor(2,n) = -3 * triangle_num - 2
    triangle_neighbor(3,triangle_num) = -l

    ltri = n
    ledg = 2

    call swapec ( m, top, ltri, ledg, node_num, node_xy, triangle_num, &
      triangle_node, triangle_neighbor, stack, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      write ( *, '(a)' ) '  Error return from SWAPEC.'
      stop
    end if

  end do
!
!  Now account for the sorting that we did.
!
  do i = 1, 3
    do j = 1, triangle_num
      triangle_node(i,j) = indx ( triangle_node(i,j) )
    end do
  end do

  call perm_inv ( node_num, indx )

  call r82vec_permute ( node_num, node_xy, indx )

  return
end
function e_measure ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! E_MEASURE determines the pointset quality measure E.
!
!  Discussion:
!
!    The E measure of point distribution quality for a set Z of
!    N points in an DIM_NUM dimensional region is defined as follows:
!
!    Assign every point X in the region to the nearest element
!    Z(I) of the point set.  For each point Z(I), let E_VEC(I) be the
!    integral of the distance between Z(I) and all the points assigned to
!    it:
!
!      E_VEC(I) = Integral ( all X nearest to Z(I) ) distance ( X, Z(I) )
!
!    If we let VOLUME be the volume of the region, then we define E by:
!
!      E = sum ( 1 <= I <= N ) E_VEC(I) / VOLUME
!
!    This quantity can be estimated by using sampling to pick a large
!    number of points in the region.
!
!    The E measure is minimized by a centroidal Voronoi tessellation.
!
!    Given two meshes, the one with a lower value of E is to be recommended.
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
!    Max Gunzburger and John Burkardt,
!    Uniformity Measures for Point Samples in Hypercubes.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the points.
!
!    Input, integer ( kind = 4 ) NS, the number of sample points.
!
!    Input, external SAMPLE_ROUTINE, the name of a routine which
!    is used to produce sample points in the region, of the form:
!      subroutine sample_routine ( dim_num, n, seed, x )
!      integer ( kind = 4 ) dim_num
!      integer ( kind = 4 ) n
!      integer ( kind = 4 ) seed
!      real ( kind = 8 ) x(dim_num,n)
!
!    Input, integer ( kind = 4 ) SEED_INIT, the initial value of the random 
!    number seed.
!
!    Output, real ( kind = 8 ) E_MEASURE, the E quality measure.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) closest(1)
  real ( kind = 8 ) dist
  real ( kind = 8 ) e_measure
  real ( kind = 8 ) e_vec(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ns
  external sample_routine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) z(dim_num,n)

  seed = seed_init
  e_vec(1:n) = 0.0D+00

  do k = 1, ns

    call sample_routine ( dim_num, 1, seed, x )

    call find_closest ( dim_num, n, 1, x, z, closest )

    dist = sum ( ( x(1:dim_num) - z(1:dim_num,closest(1) ) )**2 )

    e_vec(closest(1)) = e_vec(closest(1)) + dist

  end do

  e_measure = sum ( e_vec(1:n) ) / real ( ns, kind = 8 )

  return
end
subroutine file_column_count ( input_filename, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some
!    comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer ( kind = 4 ) column_num
  logical got_one
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( input_filename )
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    got_one = .true.
    exit

  end do

  if ( .not. got_one ) then

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = input_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

  return
end
subroutine file_row_count ( input_filename, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    row_num = row_num + 1

  end do

  close ( unit = input_unit )

  return
end
subroutine find_closest ( dim_num, n, sample_num, s, r, nearest )

!*****************************************************************************80
!
!! FIND_CLOSEST finds the nearest R point to each S point.
!
!  Discussion:
!
!    This routine finds the closest Voronoi cell generator by checking every
!    one.  For problems with many cells, this process can take the bulk
!    of the CPU time.  Other approaches, which group the cell generators into
!    bins, can run faster by a large factor.
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
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of cell generators.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of sample points.
!
!    Input, real ( kind = 8 ) S(DIM_NUM,SAMPLE_NUM), the points to be checked.
!
!    Input, real ( kind = 8 ) R(DIM_NUM,N), the cell generators.
!
!    Output, integer ( kind = 4 ) NEAREST(SAMPLE_NUM), the index of the nearest
!    cell generators.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) sample_num

  real ( kind = 8 ) dist_sq_min
  real ( kind = 8 ) dist_sq
  integer ( kind = 4 ) jr
  integer ( kind = 4 ) js
  integer ( kind = 4 ) nearest(sample_num)
  real ( kind = 8 ) r(dim_num,n)
  real ( kind = 8 ) s(dim_num,sample_num)

  do js = 1, sample_num

    dist_sq_min = huge ( dist_sq_min )
    nearest(js) = -1

    do jr = 1, n

      dist_sq = sum ( ( r(1:dim_num,jr) - s(1:dim_num,js) )**2 )

      if ( dist_sq < dist_sq_min ) then
        dist_sq_min = dist_sq
        nearest(js) = jr
      end if

    end do

  end do

  return
end
function gamma_measure ( dim_num, n, z )

!*****************************************************************************80
!
!! GAMMA_MEASURE determines the pointset quality measure GAMMA.
!
!  Discussion:
!
!    The GAMMA measure of point distribution quality for a set Z of
!    N points in an DIM_NUM-dimensional region is defined as follows:
!
!      GAMMA = ( GAMMA_MAX / GAMMA_MIN ),
!
!    where
!
!      GAMMA_MAX = maximum ( 1 <= I <= N ) DIST_MIN(I)
!      GAMMA_MIN = minimum ( 1 <= I <= N ) DIST_MIN(I)
!
!    and
!
!      DIST_MIN(I) = minimum ( 1 <= J <= N, I /= J ) distance ( Z(I), Z(J) )
!
!    Note that, in this code, the variable DIST_SQ_MIN is actually the square
!    of the minimum point distance, and so when we compute GAMMA, we must
!    take the square root of the given ratio.
!
!    GAMMA must be at least 1.  For an ideally regular mesh, GAMMA would
!    be equal to one.  Given two meshes, this measure recommends the one
!    with the smaller value of GAMMA.
!
!    Note, however, that GAMMA is only measuring how uniformly the
!    points are distributed, not how well.  The points could have
!    a perfect GAMMA value of 1, but be concentrated in one corner
!    of the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2004
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the points.
!
!    Output, real ( kind = 8 ) GAMMA_MEASURE, the GAMMA quality measure.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) GAMMA_SQ_MAX, the maximum, over all points,
!    of the minimum squared distance to a distinct point.
!
!    Local, real ( kind = 8 ) GAMMA_SQ_MIN, the minimum, over all points,
!    of the minimum squared distance to a distinct point.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) dist_sq
  real ( kind = 8 ) dist_sq_min
  real ( kind = 8 ) gamma_measure
  real ( kind = 8 ) gamma_sq_max
  real ( kind = 8 ) gamma_sq_min
  real ( kind = 8 ) z(1:dim_num,1:n)
!
!  Take care of ridiculous cases.
!
  if ( n <= 1 ) then
    gamma_measure = 0.0D+00
    return
  end if

  gamma_sq_max = 0.0D+00
  gamma_sq_min = huge ( gamma_sq_min )

  do j1 = 1, n

    dist_sq_min = huge ( dist_sq_min )

    do j2 = 1, n

      if ( j2 /= j1 ) then
        dist_sq = sum ( ( z(1:dim_num,j1) - z(1:dim_num,j2) )**2 )
        if ( dist_sq < dist_sq_min ) then
          dist_sq_min = dist_sq
        end if
      end if

    end do

    gamma_sq_max = max ( gamma_sq_max, dist_sq_min )
    gamma_sq_min = min ( gamma_sq_min, dist_sq_min )

  end do
!
!  Take care of the ridiculous case in which at least two points are identical.
!
  if ( gamma_sq_min == 0.0D+00 ) then

    if ( gamma_sq_max == 0.0D+00 ) then
      gamma_measure = 1.0D+00
    else
      gamma_measure = huge ( gamma_measure )
    end if
!
!  Take care of the general case when no points are identical.
!
  else
    gamma_measure = sqrt ( gamma_sq_max / gamma_sq_min )
  end if

  return
end
subroutine get_input_file_name ( input_file_name )

!*****************************************************************************80
!
!! GET_INPUT_FILE_NAME gets the input file name.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
  implicit none

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
! integer ( kind = 4 ) ierror
! integer ( kind = 4 ) ilen
  character ( len = * )  input_file_name
! integer ( kind = 4 ) ipxfargc
!
!  Get the number of command line arguments.
!
!  Old style:
!
  arg_num = iargc ( )
!
!  New style:
!
! arg_num = ipxfargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
!
!  Old style:
!
    call getarg ( iarg, input_file_name )
!
!  New style:
!
!   call pxfgetarg ( iarg, input_file_name, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'GET_INPUT_FILE_NAME - Fatal error!'
!     write ( *, '(a)' ) '  Could not read command line argument.'
!     stop
!   end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GET_INPUT_FILE_NAME:'
    write ( *, '(a)' ) '  Please enter the name of the file containing'
    write ( *, '(a)' ) '  the points to be used in the tests.'

    read ( *, '(a)' ) input_file_name

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
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
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
function h_measure ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! H_MEASURE determines the pointset quality measure H.
!
!  Discussion:
!
!    The H measure of dispersion for a set of N points in an DIM_NUM-dimensional
!    region is the maximum distance between any point in the region
!    and some point in the set.
!
!    To compute this quantity exactly, for every point X in the region,
!    find the nearest element Z of the point set and compute the distance.
!    H is then the maximum of all these distances.
!
!    To ESTIMATE this quantity, carry out the same process, but only for
!    NS sample points in the region.
!
!    Under this measure, a mesh with a smaller value of H is preferable.
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
!    Max Gunzburger and John Burkardt,
!    Uniformity Measures for Point Samples in Hypercubes.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the points.
!
!    Input, integer ( kind = 4 ) NS, the number of sample points.
!
!    Input, external SAMPLE_ROUTINE, the name of a routine which
!    is used to produce sample points in the region, of the form:
!      subroutine sample_routine ( dim_num, n, seed, x )
!      integer ( kind = 4 ) dim_num
!      integer ( kind = 4 ) n
!      integer ( kind = 4 ) seed
!      real ( kind = 8 ) x(dim_num,n)
!
!    Input, integer ( kind = 4 ) SEED_INIT, the initial value of the random 
!    number seed.
!
!    Output, real ( kind = 8 ) H_MEASURE, the H quality measure.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) k
  integer ( kind = 4 ) closest(1)
  real ( kind = 8 ) dist
  real ( kind = 8 ) h
  real ( kind = 8 ) h_measure
  integer ( kind = 4 ) ns
  external sample_routine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) z(dim_num,n)

  seed = seed_init
  h = 0.0D+00

  do k = 1, ns

    call sample_routine ( dim_num, 1, seed, x )

    call find_closest ( dim_num, n, 1, x, z, closest )

    dist = sum ( ( x(1:dim_num) - z(1:dim_num,closest(1)) )**2 )

    h = max ( h, dist )

  end do

  h_measure = sqrt ( h )

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer ( kind = 4 ) division.
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
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
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
!    Input, integer ( kind = 4 ) IVAL, the value to be wrapped.
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
!    Input, character ( len = * ) TITLE, a title.
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
function lambda_measure ( dim_num, n, z )

!*****************************************************************************80
!
!! LAMBDA_MEASURE determines the pointset quality measure LAMBDA.
!
!  Discussion:
!
!    The LAMBDA measure of point distribution quality for a set Z of
!    N points in an DIM_NUM-dimension region is defined as follows:
!
!    Let
!
!      GAMMA(I) = minimum ( 1 <= J <= N, I /= J ) distance ( Z(I), Z(J) )
!
!    and let
!
!      GAMMA_AVE = sum ( 1 <= I <= N ) GAMMA(I) / N
!
!    then
!
!      LAMBDA = sqrt ( sum ( 1 <= I <= N ) ( GAMMA(I) - GAMMA_AVE )**2 / N )
!        / GAMMA_AVE
!
!    An ideally regular mesh would have GAMMA(I) = GAMMA_AVE for all I,
!    so that LAMBDA would be 0.  Under this measure, the mesh with the
!    smaller value of LAMBDA is to be preferred.
!
!    This measure takes no account of how well the points are distributed
!    within the region, only how uniformly spaced they are.  In particular,
!    if the points were laid out equally spaced on a line, they would
!    get a perfect LAMBDA value, no matter whether they were supposed to
!    cover a larger line, or a plane or a cube or a circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2004
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the points.
!
!    Output, real ( kind = 8 ) LAMBDA_MEASURE, the LAMBDA quality measure.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) GAMMA_MAX, the maximum, over all points,
!    of the minimum distance to a distinct point.
!
!    Local, real ( kind = 8 ) GAMMA_MIN, the minimum, over all points,
!    of the minimum distance to a distinct point.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) gamma(n)
  real ( kind = 8 ) gamma_ave
  real ( kind = 8 ) lambda_measure
  real ( kind = 8 ) z(1:dim_num,1:n)
!
!  Take care of ridiculous cases.
!
  if ( n <= 1 ) then
    lambda_measure = 0.0D+00
    return
  end if
!
!  Compute the minimum spacing between distinct points of the set.
!
  call pointset_spacing ( dim_num, n, z, gamma )
!
!  Average the minimum spacing.
!
  gamma_ave = sum ( gamma(1:n) ) / real ( n, kind = 8 )
!
!  Compute a weighted variance.
!
  if ( gamma_ave <= 0.0D+00 ) then
    lambda_measure = huge ( lambda_measure )
  else
    lambda_measure = &
      sqrt ( &
             sum ( ( gamma(1:n) - gamma_ave )**2 ) / real ( n, kind = 8 ) &
           ) / gamma_ave
  end if

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
function mu_measure ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! MU_MEASURE determines the pointset quality measure MU.
!
!  Discussion:
!
!    The MU measure of dispersion, for a set of N points in an
!    DIM_NUM-dimensional region, takes the ratio of the largest and
!    smallest half-diameters of the Voronoi cells defined by a pointset.
!
!    To compute this quantity exactly, for every point X in the region,
!    find the nearest element Z of the point set and compute the distance.
!
!    Then, for each element Z(I) of the point set, define H(I) to be the
!    maximum of these distances.
!
!    MU is then the ratio of the maximum and minimum values of H.
!
!    To ESTIMATE this quantity, carry out the same process, but only for
!    NS sample points in the region.
!
!    In an ideally regular mesh, MU would be 1.  MU must be at least 1.
!    Under this measure, the mesh with the smaller value of MU is to be
!    preferred.
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
!    Max Gunzburger and John Burkardt,
!    Uniformity Measures for Point Samples in Hypercubes.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the points.
!
!    Input, integer ( kind = 4 ) NS, the number of sample points.
!
!    Input, external SAMPLE_ROUTINE, the name of a routine which
!    is used to produce sample points in the region, of the form:
!      subroutine sample_routine ( dim_num, n, seed, x )
!      integer ( kind = 4 ) dim_num
!      integer ( kind = 4 ) n
!      integer ( kind = 4 ) seed
!      real ( kind = 8 ) x(dim_num,n)
!
!    Input, integer ( kind = 4 ) SEED_INIT, the initial value of the random
!    number seed.
!
!    Output, real ( kind = 8 ) MU_MEASURE, the MU quality measure.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) k
  integer ( kind = 4 ) closest(1)
  real ( kind = 8 ) dist
  real ( kind = 8 ) h(n)
  real ( kind = 8 ) h_max
  real ( kind = 8 ) h_min
  real ( kind = 8 ) mu_measure
  integer ( kind = 4 ) ns
  external             sample_routine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) z(dim_num,n)

  seed = seed_init
  h(1:n) = 0.0D+00

  do k = 1, ns

    call sample_routine ( dim_num, 1, seed, x )

    call find_closest ( dim_num, n, 1, x, z, closest )

    dist = sum ( ( x(1:dim_num) - z(1:dim_num,closest(1)) )**2 )

    h(closest(1)) = max ( h(closest(1)), dist )

  end do

  h_max = sqrt ( maxval ( h(1:n) ) )
  h_min = sqrt ( minval ( h(1:n) ) )
!
!  If H_MIN is 0, then either
!  * some point is not the closest point to ANY sample point, (POSSIBLE), or
!  * some point is exactly equal to a sample point. (UNLIKELY)
!
  if ( h_min == 0.0D+00 ) then
    mu_measure = huge ( mu_measure )
    return
  end if

  mu_measure = h_max / h_min

  return
end
function nu_measure ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! NU_MEASURE determines the pointset quality measure NU.
!
!  Discussion:
!
!    The NU measure of dispersion for a set of N points in an
!    DIM_NUM-dimensional region is defined as follows:
!
!    For each element Z(I) of the pointset, let VOLUME(I) be the volume
!    of the corresponding Voronoi subregion, restricted to the region.
!
!    Then
!
!      NU = max ( 1 <= I <= N ) VOLUME(I) / min ( 1 <= I <= N ) VOLUME(I)
!
!    This quantity can be estimated by using a large number of sampling
!    points to estimate the Voronoi volumes.
!
!    For an ideally uniform pointset, the Voronoi volumes would be equal,
!    so that NU would be 1.  In any case, NU must be 1 or greater.  In
!    comparing two meshes, the one with the lower value of NU would be
!    preferred.
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
!    Max Gunzburger and John Burkardt,
!    Uniformity Measures for Point Samples in Hypercubes.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the points.
!
!    Input, integer ( kind = 4 ) NS, the number of sample points.
!
!    Input, external SAMPLE_ROUTINE, the name of a routine which
!    is used to produce sample points in the region, of the form:
!      subroutine sample_routine ( dim_num, n, seed, x )
!      integer ( kind = 4 ) dim_num
!      integer ( kind = 4 ) n
!      integer ( kind = 4 ) seed
!      real ( kind = 8 ) x(dim_num,n)
!
!    Input, integer ( kind = 4 ) SEED_INIT, the initial value of the random
!    number seed.
!
!    Output, real ( kind = 8 ) NU_MEASURE, the NU quality measure.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) closest(1)
  integer ( kind = 4 ) hit(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ns
  real ( kind = 8 ) nu_measure
  external             sample_routine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) volume(n)
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) z(dim_num,n)

  seed = seed_init
  hit(1:n) = 0

  do k = 1, ns

    call sample_routine ( dim_num, 1, seed, x )

    call find_closest ( dim_num, n, 1, x, z, closest )

    hit(closest(1)) = hit(closest(1)) + 1

  end do

  volume(1:n) = real ( hit(1:n), kind = 8 ) &
              / real ( ns, kind = 8 )
!
!  If the minimum volume is 0, then some point is not the closest
!  to any of the sample points.
!
  if ( minval ( volume(1:n) ) == 0.0D+00 ) then
    nu_measure = huge ( nu_measure )
    return
  end if

  nu_measure = maxval ( volume(1:n) ) / minval ( volume(1:n) )

  return
end
subroutine perm_inv ( n, p )

!*****************************************************************************80
!
!! PERM_INV inverts a permutation "in place".
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
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard
!    index form.  On output, P describes the inverse permutation
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) is
  integer ( kind = 4 ) p(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INV - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop
  end if

  is = 1

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      i2 = p(i1)
      p(i1) = - i2
      i1 = i2
    end do

    is = - sign ( 1, p(i) )
    p(i) = sign ( p(i), is )

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

  return
end
subroutine pointset_spacing ( dim_num, n, z, gamma )

!*****************************************************************************80
!
!! POINTSET_SPACING determines the minimum spacing between points in the set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 December 2002
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the point distribution.
!
!    Output, real ( kind = 8 ) GAMMA(N), the minimum distance between each
!    point and a distinct point in the set.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) gamma(n)
  real ( kind = 8 ) z(dim_num,n)

  gamma(1:n) = huge ( gamma(1) )

  do j1 = 1, n

    do j2 = 1, n

      if ( j2 /= j1 ) then
        gamma(j1) = min ( gamma(j1), &
          sum ( ( z(1:dim_num,j1) - z(1:dim_num,j2) )**2 ) )
      end if

    end do

  end do

  gamma(1:n) = sqrt ( gamma(1:n) )

  return
end
function q_measure ( n, z, triangle_order, triangle_num, triangle_node )

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
!    a triangulation is a list of TRIANGLE_NUM triples of node indices that form
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
!    07 November 2005
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
!    Input, integer ( kind = 4 ) TRIANGLE_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM),
!    the triangulation.
!
!    Output, real ( kind = 8 ) Q_MEASURE, the Q quality measure.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) triangle_order

  integer ( kind = 4 ) a_index
  real ( kind = 8 ) ab_length
  integer ( kind = 4 ) b_index
  real ( kind = 8 ) bc_length
  integer ( kind = 4 ) c_index
  real ( kind = 8 ) ca_length
  real ( kind = 8 ) q
  real ( kind = 8 ) q_min
  real ( kind = 8 ) q_measure
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle_node(triangle_order,triangle_num)
  real ( kind = 8 ) z(2,n)

  if ( triangle_num < 1 ) then
    q_measure = -1.0D+00
    return
  end if

  q_min = huge ( q_min )

  do triangle = 1, triangle_num

    a_index = triangle_node(1,triangle)
    b_index = triangle_node(2,triangle)
    c_index = triangle_node(3,triangle)

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

    q_min = min ( q_min, q )

  end do

  q_measure = q_min

  return
end
function r0_measure ( dim_num, n, z )

!*****************************************************************************80
!
!! R0_MEASURE determines the pointset quality measure R0.
!
!  Discussion:
!
!    The R0 measure of point distribution quality for a set Z of
!    N points in an DIM_NUM-dimensional region is defined as follows:
!
!      R0 = sum ( 1 <= I /= J <= N ) log ( 1 / distance ( Z(I), Z(J) ) )
!         / ( N * ( N - 1 ) )
!
!    The divisor of ( N * ( N - 1 ) ) means that R0 is
!    undefined if N < 2 or if any two points are equal.
!
!    R0 is known as the Riesz S-energy for S = 0.
!
!    Given two meshes, this measure recommends the one with the smaller
!    value of R0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Douglas Hardin, Edward Saff,
!    Discretizing Manifolds via Minimum Energy Points,
!    Notices of the AMS,
!    Volume 51, Number 10, November 2004, pages 1186-1194.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the points.
!
!    Output, real ( kind = 8 ) R0_MEASURE, the R0 quality measure.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) dist
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) r0_measure
  real ( kind = 8 ) z(1:dim_num,1:n)
!
!  Take care of ridiculous cases.
!
  if ( n <= 1 ) then
    r0_measure = huge ( z(1,1) )
    return
  end if

  r0_measure = 0.0D+00

  do j1 = 1, n

    do j2 = 1, n

      if ( j2 /= j1 ) then

        dist = sqrt ( sum ( ( z(1:dim_num,j1) - z(1:dim_num,j2) )**2 ) )

        if ( dist == 0.0D+00 ) then
          r0_measure = huge ( z(1,1) )
          return
        end if

        r0_measure = r0_measure + log ( 1.0D+00 / dist )

      end if

    end do

  end do

  if ( 1 < n ) then
    r0_measure = r0_measure / real ( n * ( n - 1 ), kind = 8 )
  else
    r0_measure = 0.0D+00
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
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer ( kind = 4 ) arithmetic never requires more than 32 bits,
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
!  Although SEED can be represented exactly as a 32 bit integer ( kind = 4 ),
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r82vec_permute ( n, a, p )

!*****************************************************************************80
!
!! R82VEC_PERMUTE permutes an R82 vector in place.
!
!  Discussion:
!
!    This routine permutes an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
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
!    11 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, real ( kind = 8 ) A(2,N), the array to be permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integer ( kind = 4 )s from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) a_temp(2)
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)
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

      a_temp(1:2) = a(1:2,istart)
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
          stop
        end if

        if ( iget == istart ) then
          a(1:2,iput) = a_temp(1:2)
          exit
        end if

        a(1:2,iput) = a(1:2,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine r82vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(1:2,INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call R82VEC_PERMUTE ( N, A, INDX )
!
!    after which A(1:2,I), I = 1 to N is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2004
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

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) aval(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then
    indx(1) = 1
    return
  end if

  call i4vec_indicator ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval(1:2) = a(1:2,indxt)

    else

      indxt = indx(ir)
      aval(1:2) = a(1:2,indxt)
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
subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_READ reads data from an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine will
!    return after reading N of them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  real ( kind = 8 ) table(m,n)
  real ( kind = 8 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine r8mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! R8MAT_HEADER_READ reads the header from an R8MAT file.
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
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
function r8mat_in_01 ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_IN_01 is TRUE if the entries of an R8MAT are in the range [0,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2004
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
!    Output, logical R8MAT_IN_01, is TRUE if every entry of A is
!    between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  logical r8mat_in_01

  if ( any ( a(1:m,1:n) < 0.0D+00 .or. 1.0D+00 < a(1:m,1:n) ) ) then
    r8mat_in_01 = .false.
  else
    r8mat_in_01 = .true.
  end if

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
!    An R8MAT is an array of real ( kind = 8 ) values.
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
subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED
!    in order to set the seed of the random number generator.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is negative,
!    then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer ( kind = 4 ) SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
!    Local, real ( kind = 8 ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) m
  integer ( kind = 4 ), save :: made = 0
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r(n+1)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ), save :: saved = 0
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
  real ( kind = 8 ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

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
subroutine radius_maximus ( dim_num, n, z, walls, radius )

!*****************************************************************************80
!
!! RADIUS_MAXIMUS finds the biggest possible nonintersecting sphere.
!
!  Discussion:
!
!    We are given a set of N points in DIM_NUM space.  We imagine that
!    at each point simultaneously, a sphere begins to expand.
!    Each sphere stops expanding as soon as it touches another sphere.
!    The radius of these spheres is to be computed.
!
!    If WALLS is true, then the spheres must not extend outside the
!    "walls" of the unit hypercube.
!
!    This routine assumes that the points are contained in the unit
!    hypercube.  It could easily be modified to handle points in an
!    arbitrary hypercube, but general regions would be much harder.
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
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the point coordinates.
!    If WALLS is TRUE, these values must be between 0 and 1.
!
!    Input, logical WALLS, is TRUE if the spheres must not extend
!    outside the unit hypercube.  If WALLS is FALSE, then this
!    restriction is not imposed.
!
!    Output, real ( kind = 8 ) RADIUS(N), the radius of the
!    maximal nonintersecting sphere around each point.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) distance_j
  integer ( kind = 4 ), parameter :: FIXED = 0
  integer ( kind = 4 ), parameter :: FREE = 1
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) next
  real ( kind = 8 ) radius(n)
  real ( kind = 8 ) radius_i
  real ( kind = 8 ) radius_min
  integer ( kind = 4 ) status(n)
  logical walls
  real ( kind = 8 ) z(dim_num,n)

  if ( walls ) then
         if ( any ( z(1:dim_num,1:n) < 0.0D+00 ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RADIUS_MAXIMUS - Fatal error!'
      write ( *, '(a)' ) '  Some coordinate is less than 0.'
      return
    else if ( any ( 1.0D+00 < z(1:dim_num,1:n) ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RADIUS_MAXIMUS - Fatal error!'
      write ( *, '(a)' ) '  Some coordinate is greater than 1.'
      return
    end if
  end if
!
!  Initially, all points are "free".
!
  radius(1:n) = 0.0D+00
  status(1:n) = FREE

  do
!
!  If all points are fixed, we're done.
!
    if ( all ( status(1:n) == FIXED ) ) then
      exit
    end if
!
!  Look at all the free points.
!  Imagine an expanding sphere at each free point, and determine
!  which such sphere will first have to stop expanding.
!
    next = -1
    radius_min = huge ( radius_min )

    do j1 = 1, n

      if ( status(j1) == FREE ) then

        if ( walls ) then
          radius_i = min ( &
            minval (           z(1:dim_num,j1) ), &
            minval ( 1.0D+00 - z(1:dim_num,j1) ) )
        else
          radius_i = huge ( radius_i )
        end if

        do j2 = 1, n

          if ( j2 /= j1 ) then

            distance_j = sqrt ( sum ( &
              ( z(1:dim_num,j1) - z(1:dim_num,j2) )**2 &
            ) )

            if ( status(j2) == FREE ) then
              radius_i = min ( radius_i, distance_j / 2.0D+00 )
            else
              radius_i = min ( radius_i, distance_j - radius(j2) )
            end if

          end if

        end do

        if ( radius_i < radius_min ) then
          next = j1
          radius_min = radius_i
        end if

      end if

    end do

    if ( next == -1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RADIUS_MAXIMUS - Fatal error!'
      write ( *, '(a)' ) '  There were points left to handle, but could'
      write ( *, '(a)' ) '  not choose the "next" one to work on.'
      stop
    end if

    radius(next) = radius_min
    status(next) = FIXED

  end do

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer ( kind = 4 ) part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer ( kind = 4 ) part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the real value that was read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 .or. nchar <= lchar+1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lchar
  real ( kind = 8 ) rvec(n)
  character ( len = * )  s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

  end do

  return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
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
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) nword
  character ( len = * ) s

  nword = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      nword = nword + 1
      blank = .false.
    end if

  end do

  return
end
subroutine sample_hypercube_uniform ( dim_num, n, seed, x )

!*****************************************************************************80
!
!! SAMPLE_HYPERCUBE_UNIFORM returns sample points in the unit hypercube.
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
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points to compute.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,N), the sample points.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(dim_num,n)

  call r8mat_uniform_01 ( dim_num, n, seed, x )

  return
end
subroutine sample_sphere_uniform ( m, n, seed, x )

!*****************************************************************************80
!
!! SAMPLE_SPHERE_UNIFORM samples points inside the unit sphere.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    We first generate a point ON the sphere, and then distribute it
!    IN the sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Russell Cheng,
!    Random Variate Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998, pages 168.
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity
!    of Queueing Networks,
!    Wiley, 1986, page 232.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the space.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) X(M,N), the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) exponent
  integer ( kind = 4 ) j
  real ( kind = 8 ) norm
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(m,n)

  exponent = 1.0D+00 / real ( m, kind = 8 )

  do j = 1, n
!
!  Fill a vector with normally distributed values.
!
    call r8vec_normal_01 ( m, seed, x(1:m,j) )
!
!  Compute the length of the vector.
!
    norm = sqrt ( sum ( x(1:m,j)**2 ) )
!
!  Normalize the vector.
!
    x(1:m,j) = x(1:m,j) / norm
!
!  Now compute a value to map the point ON the sphere INTO the sphere.
!
    r = r8_uniform_01 ( seed )

    x(1:m,j) = r**exponent * x(1:m,j)

  end do

  return
end
function sphere_measure ( dim_num, n, z )

!*****************************************************************************80
!
!! SPHERE_MEASURE determines the pointset quality measure S.
!
!  Discussion:
!
!    This routine computes a measure of even spacing in a set of N
!    points in the DIM_NUM-dimensional unit hypercube.  We will discuss
!    the program as though the space is 2-dimensional and the spheres
!    are circles, but the program may be used for general
!    DIM_NUM-dimensional data.
!
!    The points are assumed to lie in the unit square.
!
!    The program makes a circle-packing measurement on the points
!    by assuming that, at each point, a circle is centered; all
!    the circles start out with zero radius, and then expand
!    together at the same rate.  A circle stops expanding as soon
!    as it touches any other circle.
!
!    The amount of area covered by the circles is compared to the
!    area of the unit square.  This measurement has a certain amount
!    of boundary effect: some circles will naturally extend outside
!    the unit hypercube.  If this is a concern, is possible to restrict
!    the circles to remain inside the unit hypercube.  In any case,
!    this problem generally goes away as the number of points increases.
!
!    Since we are interested in the coverage of the unit hypercube,
!    it is probably best if the circles are restricted.  This way,
!    computing the area of the circles gives a measure of the even
!    coverage of the region, relative to the presumably best possible
!    covering, by the same number of circles, but of equal radius.
!
!    In the limit, the maximum relative packing density of a 2D
!    region with equal-sized circles is 0.9069.  In 3D, a density
!    of at least 0.74 can be achieved, and it is known that no
!    greater than 0.7796 is possible.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the points.
!
!    Output, real ( kind = 8 ) SPHERE_MEASURE, the amount of volume taken up
!    by the nonintersecting spheres of maximum radius around each
!    point.  Ignoring boundary effects, the "ideal" value would be
!    1 (achievable only in 1 dimension), and the maximum value
!    possible is the sphere packing density in the given spatial
!    dimension.  If boundary effects can be ignored, the value of
!    SPHERE_VOLUME reports how closely the given set of points
!    behaves like a set of close-packed spheres.
!
!  Local Parameters:
!
!    Local, logical WALLS, is TRUE if the spheres are restricted
!    to lie within the unit hypercube.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  logical r8mat_in_01
  real ( kind = 8 ), dimension ( n ) :: radius
  real ( kind = 8 ) radius_ave
  real ( kind = 8 ) radius_max
  real ( kind = 8 ) radius_min
  real ( kind = 8 ) sphere_measure
  logical, parameter :: verbose = .false.
  real ( kind = 8 ) volume
  logical, parameter :: walls = .true.
  real ( kind = 8 ), dimension ( dim_num, n ) :: z

  if ( .not. r8mat_in_01 ( dim_num, n, z ) ) then
!   write ( *, '(a)' ) ' '
!   write ( *, '(a)' ) 'SPHERE_MEASURE - Fatal error!'
!   write ( *, '(a)' ) '  Some data is not in the unit hypercube.'
    sphere_measure = huge ( 1.0D+00 )
    return
  end if

  call radius_maximus ( dim_num, n, z, walls, radius )

  sphere_measure = 0.0D+00
  do j = 1, n
    call sphere_volume_nd ( dim_num, radius(j), volume )
    sphere_measure = sphere_measure + volume
  end do

  if ( verbose ) then

    radius_ave = sum ( radius(1:n) ) / real ( n, kind = 8 )
    radius_min = minval ( radius(1:n) )
    radius_max = maxval ( radius(1:n) )

    write ( *, '(a)'      ) ' '
    write ( *, '(a,i8)'   ) '  Number of dimensions is ', dim_num
    write ( *, '(a,i8)'   ) '  Number of points is ', n
    if ( walls ) then
      write ( *, '(a)' ) &
        '  Spheres are required to stay in the unit hypercube.'
    else
      write ( *, '(a)' ) &
        '  Spheres are NOT required to stay in the unit hypercube.'
    end if
    write ( *, '(a)'      ) ' '
    write ( *, '(a,f7.4)' ) '  Average radius = ', radius_ave
    write ( *, '(a,f7.4)' ) '  Minimum radius = ', radius_min
    write ( *, '(a,f7.4)' ) '  Maximum radius = ', radius_max
    write ( *, '(a,f7.4)' ) '  Sphere volume =  ', sphere_measure
  end if

  return
end
subroutine sphere_volume_nd ( dim_num, r, volume )

!*****************************************************************************80
!
!! SPHERE_VOLUME_ND computes the volume of a sphere in ND.
!
!  Formula:
!
!    A sphere in ND satisfies the equation:
!
!      sum ( ( X(1:N) - XC(1:N) )**2 ) = R**2
!
!    where R is the radius and XC is the center.
!
!  Discussion:
!
!    N  Volume
!
!    2             PI    * R**2
!    3  (4/3)    * PI    * R**3
!    4  (1/2)    * PI**2 * R**4
!    5  (8/15)   * PI**2 * R**5
!    6  (1/6)    * PI**3 * R**6
!    7  (16/105) * PI**3 * R**7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the sphere.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) volume

  if ( mod ( dim_num, 2 ) == 0 ) then
    m = dim_num / 2
    volume = pi**m
    do i = 1, m
      volume = volume / real ( i, kind = 8 )
    end do
  else
    m = ( dim_num - 1 ) / 2
    volume = pi**m * 2.0D+00**dim_num
    do i = m+1, 2*m+1
      volume = volume / real ( i, kind = 8 )
    end do
  end if

  volume = volume * r**dim_num

  return
end
subroutine swapec ( i, top, btri, bedg, node_num, node_xy, triangle_num, &
  triangle_node, triangle_neighbor, stack, ierr )

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
!    Input/output, integer ( kind = 4 ) BTRI, BEDG; on input, if positive, are the
!    triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
!
!    Input/output, integer ( kind = 4 ) TRIANGLE_NODE(3,TRIANGLE_NUM), the
!    triangle incidence list.  May be updated on output because of swaps.
!
!    Input/output, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the
!    triangle neighbor list; negative values are used for links of the
!    counter-clockwise linked list of boundary edges;  May be updated on output
!    because of swaps.
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer ( kind = 4 ) STACK(MAXST); on input, entries 1 through TOP
!    contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer ( kind = 4 ) IERR is set to 8 for abnormal return.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) triangle_num

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
  integer ( kind = 4 ) triangle_node(3,triangle_num)
  integer ( kind = 4 ) triangle_neighbor(3,triangle_num)
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

    if ( triangle_node(1,t) == i ) then
      e = 2
      b = triangle_node(3,t)
    else if ( triangle_node(2,t) == i ) then
      e = 3
      b = triangle_node(1,t)
    else
      e = 1
      b = triangle_node(2,t)
    end if

    a = triangle_node(e,t)
    u = triangle_neighbor(e,t)

    if ( triangle_neighbor(1,u) == t ) then
      f = 1
      c = triangle_node(3,u)
    else if ( triangle_neighbor(2,u) == t ) then
      f = 2
      c = triangle_node(1,u)
    else
      f = 3
      c = triangle_node(2,u)
    end if

    swap = diaedg ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,c), &
      node_xy(2,c), node_xy(1,b), node_xy(2,b) )

    if ( swap == 1 ) then

      em1 = i4_wrap ( e - 1, 1, 3 )
      ep1 = i4_wrap ( e + 1, 1, 3 )
      fm1 = i4_wrap ( f - 1, 1, 3 )
      fp1 = i4_wrap ( f + 1, 1, 3 )

      triangle_node(ep1,t) = c
      triangle_node(fp1,u) = i

      r = triangle_neighbor(ep1,t)
      s = triangle_neighbor(fp1,u)

      triangle_neighbor(ep1,t) = u
      triangle_neighbor(fp1,u) = t
      triangle_neighbor(e,t) = s
      triangle_neighbor(f,u) = r

      if ( 0 < triangle_neighbor(fm1,u) ) then
        top = top + 1
        stack(top) = u
      end if

      if ( 0 < s ) then

        if ( triangle_neighbor(1,s) == u ) then
          triangle_neighbor(1,s) = t
        else if ( triangle_neighbor(2,s) == u ) then
          triangle_neighbor(2,s) = t
        else
          triangle_neighbor(3,s) = t
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

        do while ( 0 < triangle_neighbor(ee,tt) )

          tt = triangle_neighbor(ee,tt)

          if ( triangle_node(1,tt) == a ) then
            ee = 3
          else if ( triangle_node(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        triangle_neighbor(ee,tt) = l

      end if

      if ( 0 < r ) then

        if ( triangle_neighbor(1,r) == t ) then
          triangle_neighbor(1,r) = u
        else if ( triangle_neighbor(2,r) == t ) then
          triangle_neighbor(2,r) = u
        else
          triangle_neighbor(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( 0 < triangle_neighbor(ee,tt) )

          tt = triangle_neighbor(ee,tt)

          if ( triangle_node(1,tt) == b ) then
            ee = 3
          else if ( triangle_node(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        triangle_neighbor(ee,tt) = l

      end if

    end if

  end do

  return
end
function tau_measure ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! TAU_MEASURE determines the pointset quality measure TAU.
!
!  Discussion:
!
!    The TAU measure of point distribution quality for a set Z of
!    N points in an DIM_NUM-dimensional region is defined as follows:
!
!    For each point Z(I) in the pointset, let V(I) be the subregion
!    defined by the intersection of the region with the Voronoi
!    subregion associated with Z(I).
!
!    Let T(I) be the trace of the second moment tensor about the point
!    Z(I), associated with the region V(I).  Let T_BAR be the average
!    of the values of T(1:N).
!
!    Then TAU = maximum ( 1 <= I <= N ) abs ( T(I) - TBAR ).
!
!    This quantity can be estimated using sampling.  A given number of
!    sample points are generated in the region, assigned to the nearest
!    element of the pointset, and used to approximate the Voronoi regions
!    and the second moment tensors.
!
!    In an ideally regular mesh, the values of T would be equal, and so
!    TAU would be zero.  In general, the smaller TAU, the better.
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
!    Max Gunzburger and John Burkardt,
!    Uniformity Measures for Point Samples in Hypercubes.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) Z(DIM_NUM,N), the point distribution.
!
!    Input, integer ( kind = 4 ) NS, the number of sample points.
!
!    Input, external SAMPLE_ROUTINE, the name of a routine which
!    is used to produce sample points in the region, of the form:
!      subroutine sample_routine ( dim_num, n, seed, x )
!      integer ( kind = 4 ) dim_num
!      integer ( kind = 4 ) n
!      integer ( kind = 4 ) seed
!      real ( kind = 8 ) x(dim_num,n)
!
!    Input, integer ( kind = 4 ) SEED_INIT, the initial value of the random 
!    number seed.
!
!    Output, real ( kind = 8 ) TAU_MEASURE, a quality measure.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) centroid(dim_num,n)
  integer ( kind = 4 ) closest(1)
  integer ( kind = 4 ) hit(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) moment(dim_num,dim_num,n)
  integer ( kind = 4 ) ns
  external sample_routine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) t_bar
  real ( kind = 8 ) tau_measure
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) z(dim_num,n)

  seed = seed_init
  centroid(1:dim_num,1:n) = 0.0D+00
  hit(1:n) = 0
  moment(1:dim_num,1:dim_num,1:n) = 0.0D+00

  do k = 1, ns

    call sample_routine ( dim_num, 1, seed, x )

    call find_closest ( dim_num, n, 1, x, z, closest )

    hit(closest(1)) = hit(closest(1)) + 1

    centroid(1:dim_num,closest(1)) = centroid(1:dim_num,closest(1)) &
      + x(1:dim_num)

    do i1 = 1, dim_num
      do i2 = 1, dim_num
        moment(i1,i2,closest(1)) = moment(i1,i2,closest(1)) + x(i1) * x(i2)
      end do
    end do

  end do

  do j = 1, n

    if ( 0 < hit(j) ) then

      centroid(1:dim_num,j) = centroid(1:dim_num,j) / real ( hit(j), kind = 8 )

      moment(1:dim_num,1:dim_num,j) = moment(1:dim_num,1:dim_num,j) &
        / real ( hit(j), kind = 8 )

      do i1 = 1, dim_num
        do i2 = 1, dim_num
          moment(i1,i2,j) = moment(i1,i2,j) - centroid(i1,j) * centroid(i2,j)
        end do
      end do

    end if

  end do

  t(1:n) = 0.0D+00

  do j = 1, n
    do i = 1, dim_num
      t(j) = t(j) + moment(i,i,j)
    end do
  end do

  t_bar = sum ( t(1:n) ) / real ( n, kind = 8 )

  tau_measure = maxval ( abs ( t(1:n) - t_bar ) )

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
subroutine vbedg ( x, y, node_num, node_xy, triangle_num, triangle_node, &
  triangle_neighbor, ltri, ledg, rtri, redg )

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
!    03 November 2008
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
!    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(3,TRIANGLE_NUM), the triangle
!    incidence list.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the triangle
!    neighbor list; negative values are used for links of a counter clockwise!
!    linked list of boundary edges;
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
  integer ( kind = 4 ) triangle_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  logical ldone
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  real ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) t
  integer ( kind = 4 ) triangle_node(3,triangle_num)
  integer ( kind = 4 ) triangle_neighbor(3,triangle_num)
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

    l = -triangle_neighbor(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = triangle_node(e,t)

    if ( e <= 2 ) then
      b = triangle_node(e+1,t)
    else
      b = triangle_node(1,t)
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

    b = triangle_node(e,t)
    e = i4_wrap ( e-1, 1, 3 )

    do while ( 0 < triangle_neighbor(e,t) )

      t = triangle_neighbor(e,t)

      if ( triangle_node(1,t) == b ) then
        e = 3
      else if ( triangle_node(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = triangle_node(e,t)

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
