program main

!*****************************************************************************80
!
!! CVT_FIXED2_PRB tests the CVT library with some fixed generators.
!
!  Discussion:
!
!    The idea here is simply that we are going to generate a CVT
!    (Centroidal Voronoi Tessellation), but that some of the
!    generators are going to be fixed.  We will set this up in
!    such a way that the points can be fixed or freed at any time,
!    based on an array called CELL_FIXED.
!
!    If a cell is free, then once the cell centroid has been
!    computed, that value replaces the old cell generator.
!
!    But a fixed cell does not update its generator.
!
!    This calling program is set up to work on the "double hex" problem
!    which is a region that is a subset of the unit square, with two
!    hexagonal holes removed.
!
!
!    The region is specified by the routine REGION_SAMPLE.
!
!    The fixed nodes are listed in a file whose name is assigned to
!    INPUT_FILENAME.
!
!    The total number of nodes to generate (fixed plus free) should
!    be set as the value N.
!
!    It can be useful to adjust the value of SAMPLE_NUM, which determines
!    the accuracy of each iteration, or to change to value of IT_MAX,
!    which sets the number of iterations to take.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) CELL_CENTROID(M,N), the Voronoi cell centroids.
!
!    Local, integer CELL_FIXED(N), is 0 for "free" generators,
!    and 1 for "fixed" generators which cannot move.
!
!    Local, real ( kind = 8 ) CELL_GENERATOR(M,N), the Voronoi
!    cell generators.
!
!    Local, integer N, the number of generators.
!
!    Local, integer CENTROID_IT, a counter for the CVT centroid
!    iteration.
!
!    Local, integer CENTROID_IT_MAX, the maximum number of steps
!    of the CVT centroid iteration.
!
!    Local, real ( kind = 8 ) CUTOFF_DIST, the maximum influence distance.
!    If the nearest generator is further than this distance,
!    the point is not assigned.
!
!    Local, integer M, the spatial dimension.
!
!    Local, integer SAMPLE_NUM, the number of sample points.
!
!    Local, integer SEED, the seed for the random number generator.
!    Setting seed = 0 will cause the program to choose a seed
!    for the random number generator using the current date and time.
!
  implicit none

  integer, parameter :: m = 2
  integer, parameter :: n = 139

  real ( kind = 8 ) diff
  logical fixed(n)
  integer fixed_num
  real ( kind = 8 ) generator(m,n)
  character ( len = 80 ) :: input_filename = 'double_hex2.txt'
  integer ios
  integer, parameter :: it_max = 500
  integer it_num
  integer m2
  character ( len = 80 ) :: output_filename = 'double_hex2_cvt.txt'
  integer output_unit
  integer, parameter :: sample_num = 10000 * n
  integer :: seed = 123456789
  integer seed_init
!
!  Print introduction and options.
!
  call timestamp ( )

  seed_init = seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_FIXED2_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CVT_FIXED2 library.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate a CVT algorithm in which some'
  write ( *, '(a)' ) '  generators are held fixed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Geometry parameters:'
  write ( *, '(a)' ) '-------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension is M = ', m
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT Algorithm parameters:'
  write ( *, '(a)' ) '-------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of Voronoi cells to generate: ', n
  write ( *, '(a,i12)' ) '  Number of sampling points: ', sample_num
  write ( *, '(a,i8)' ) '  Number of centroid iterations: ', it_max
  write ( *, '(a,i12)' ) &
    '  Initial random number seed SEED = ', seed
  write ( *, '(a)' ) &
    '  Voronoi cell generators are initialized by RANDOM_NUMBER.'

  write ( *, '(a)' ) ' '
!
!  Read the fixed points from the input file.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the fixed points:'

  call dtable_header_read ( input_filename, m2, fixed_num )

  if ( m2 /= m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_READ - Fatal error!'
    write ( *, '(a)' ) '  The spatial dimension is illegal.'
    stop
  end if

  if ( n <= fixed_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_READ - Fatal error!'
    write ( *, '(a)' ) '  The number of fixed points is too large.'
    stop
  end if

  call dtable_data_read ( input_filename, m, fixed_num, generator )

  fixed(1:fixed_num) = .true.
  fixed(fixed_num+1:n) = .false.
!
!  Randomly assign the other points.
!
  call region_sample ( m, n - fixed_num, seed, &
    generator(1:m,fixed_num+1:n) )
!
!  Carry out the CVT iteration, which drives the Voronoi generators
!  and Voronoi centroids closer and closer.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Carry out the CVT iteration.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step  Change'
  write ( *, '(a)' ) ' '

  do it_num = 1, it_max

    call cvt_fixed_iteration ( m, n, sample_num, fixed, seed, generator, diff )

    write ( *, '(i8,2x,f12.6)' ) it_num, diff

  end do
!
!  Write the points to the output file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  call dtable_header_write ( output_filename, output_unit, m, n )

  write ( output_unit, '(a)' ) '#  This data was generated by CVT_FIXED2.'
  write ( output_unit, '(a)' ) '#'
  write ( output_unit, '(a,i8)' ) '#  Number of CVT iterations =    ', it_max
  write ( output_unit, '(a,i12)' ) &
    '#  Number of sampling points =   ', sample_num
  write ( output_unit, '(a,g14.6)' ) '#  L2 change in last iteration = ', diff
  write ( output_unit, '(a,i12)' ) &
    '#  Initial random number seed =  ', seed_init
  write ( output_unit, '(a,i8)' ) '#  The number of fixed nodes = ', fixed_num
  write ( output_unit, '(a)' ) '#  (The fixed points are listed first.)'
  write ( output_unit, '(a)' ) '#'

  call dtable_data_write ( output_unit, m, n, generator )

  close ( unit = output_unit )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_FIXED2_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine region_sample ( m, n, seed, x )

!******************************************************************************
!
!! REGION_SAMPLE samples the user's region.
!
!  Discussion:
!
!    This particular region comprises the unit square, less two hexagonal holes.
!
!    The purpose of this routine is to return a random point within the
!    region.  The user could also choose to apply a weighting function
!    to the point distribution, if appropriate.
!
!    The information about this region was supplied by Jeff Borggaard
!    of Virginia Tech.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points requested.
!
!    Input/output, integer SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) X(M,N), the sample points.
!
  implicit none

  integer m
  integer n

  real ( kind = 8 ) d_uniform_01
  real ( kind = 8 ), dimension ( 2, 6 ) :: hex1 = reshape ( (/ &
    0.1500D+00,  0.7500D+00, &
    0.2003D+00,  0.6634D+00, &
    0.3002D+00,  0.6640D+00, &
    0.3500D+00,  0.7500D+00, &
    0.3002D+00,  0.8360D+00, &
    0.2003D+00,  0.8366D+00 /), (/ 2, 6 /) )
  real ( kind = 8 ), dimension ( 2, 6 ) :: hex2 = reshape ( (/ &
    0.5000D+00,  0.4000D+00, &
    0.5503D+00,  0.3134D+00, &
    0.6502D+00,  0.3140D+00, &
    0.7000D+00,  0.4000D+00, &
    0.6502D+00,  0.4860D+00, &
    0.5503D+00,  0.4866D+00 /), (/ 2, 6 /) )
  logical hexagon_contains_point_2d
  integer i
  integer j
  integer reject
  integer seed
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) y(m)

  reject = 0

  do j = 1, n

    do

      do i = 1, m
        y(i) = d_uniform_01 ( seed )
      end do

      if ( &
        ( .not. hexagon_contains_point_2d ( hex1, y ) ) &
      .and. &
        ( .not. hexagon_contains_point_2d ( hex2, y ) ) &
      ) then
        exit
      end if

      reject = reject + 1

      if ( 2 * n + 10 <= reject ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DOUBLE_HEX_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  (The double hexagonal hole region)'
        write ( *, '(a,i8)' ) '  Trying to generate point J = ', j
        write ( *, '(a,i8)' ) '  Number of rejections = ', reject
        write ( *, '(a,g14.6)' ) '  Rejection percentage = ', &
          real ( 100 * reject, kind = 8 ) &
        / real ( reject + j - 1, kind = 8 )
        write ( *, '(a,2g14.6)' ) '  Y = ', y(1:2)
        stop
      end if

    end do

    x(1:m,j) = y(1:m)

  end do

  return
end
function hexagon_contains_point_2d ( h, x )

!*****************************************************************************80
!
!! HEXAGON_CONTAINS_POINT_2D finds if a point is inside a hexagon in 2D.
!
!  Discussion:
!
!    The only restriction on the hexagon is that it be convex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) H(2,6), the hexagon vertex coordinates, listed in
!    clockwise or counterclockwise order.
!
!    Input, real ( kind = 8 ) X(2), the point to be tested.
!
!    Output, logical HEXAGON_CONTAINS_POINT_2D, is TRUE if X is in the hexagon.
!
  implicit none

  integer, parameter :: n = 6

  real ( kind = 8 ) h(2,n)
  logical hexagon_contains_point_2d
  integer i
  integer j
  real ( kind = 8 ) x(2)
!
!  A point is inside a convex hexagon if and only if it is inside
!  one of the triangles formed by the first vertex and any two consecutive
!  vertices on the polygon's circumference.
!
  do i = 2, n-1

    j = i + 1

    if ( ( x(1) - h(1,1) ) * ( h(2,i) - h(2,1) ) &
       - ( x(2) - h(2,1) ) * ( h(1,i) - h(1,1) ) <= 0.0D+00 ) then

      if ( ( x(1) - h(1,i) ) * ( h(2,j) - h(2,i) ) &
         - ( x(2) - h(2,i) ) * ( h(1,j) - h(1,i) ) <= 0.0D+00 ) then

        if ( ( x(1) - h(1,j) ) * ( h(2,1) - h(2,j) ) &
           - ( x(2) - h(2,j) ) * ( h(1,1) - h(1,j) ) <= 0.0D+00 ) then

          hexagon_contains_point_2d = .true.
          return

        end if
      end if
    end if
  end do

  hexagon_contains_point_2d = .false.

  return
end
