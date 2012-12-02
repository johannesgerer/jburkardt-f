program main

!******************************************************************************
!
!! MAIN is the main program for CVT_FIXED_PRB.
!
!  Discussion:
!
!    CVT_FIXED_PRB tests the CVT library with some fixed generators.
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
!
!    A second feature has been added, a maximum influence distance.
!    This is the variable "cutoff_dist".  When sampling the region,
!    if the nearest generator is further than this distance, we
!    don't assign the point at all.
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
!  Local parameters:
!
!    Local, real ( kind = 8 ) BOX_MIN(DIM_NUM), BOX_MAX(DIM_NUM), the
!    coordinates of the two extreme corners of the bounding box.
!
!    Local, real ( kind = 8 ) CELL_CENTROID(DIM_NUM,CELL_NUM), the Voronoi
!    cell centroids.
!
!    Local, integer CELL_FIXED(CELL_NUM), is 0 for "free" generators,
!    and 1 for "fixed" generators which cannot move.
!
!    Local, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,cell_num), the Voronoi
!    cell generators.
!
!    Local, integer CELL_NUM, the number of generators.
!
!    Local, integer CELL_SAMPLE, roughly the number of sampling
!    points we'd like per cell, assuming the cells have about the
!    same area.
!
!    Local, real ( kind = 8 ) CELL_VOLUME(CELL_NUM), the (estimated) volume of
!    each cell.
!
!    Local, real ( kind = 8 ) CELL_WEIGHT(CELL_NUM), a value used to influence
!    the area of a cell.  If the weights are equal, then the cells
!    will tend to have equal area.  If one weight is 4 times
!    another weight, then roughly speaking, the corresponding cell
!    will tend to be about 4 times as large, although this is only
!    a crude control.
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
!    Local, integer DIM_NUM, the spatial dimension.
!
!    Local integer PPMB_NROW, PPMB_NCOL, the number of rows
!    and columns of pixels to use in creating a PPMB image of
!    the CVT.
!
!    Local, logical REGION_PLOT, is true if a PPMB image of the
!    CVT is to be generated after every centroid iteration.
!
!    Local, character ( len = * ) REGION_PLOT_FILE_NAME, the
!    "family name" of the files that will hold the PPMB images.
!    The numeric characters in the name will be incremented
!    to generate new names as the iteration proceeds.
!
!    Local, real ( kind = 8 ) REGION_VOLUME, the volume of the region that
!    is being divided by the CVT algorithm.
!
!    Local, integer SAMPLE_NUM, the number of sample points, which
!    is the product of cell_num, the number of cells and
!    cell_sample, the number of sample points per cell.
!
!    Local, integer SEED, the seed for the random number generator.
!    Setting seed = 0 will cause the program to choose a seed
!    for the random number generator using the current date and time.
!
  implicit none

  integer   ( kind = 4 ), parameter :: cell_num = 14
  integer   ( kind = 4 ), parameter :: dim_num = 2

  real      ( kind = 8 ), parameter, dimension ( dim_num ) :: box_min = &
    (/ 0.0D+00, 0.0D+00 /)
  real      ( kind = 8 ), parameter, dimension ( dim_num ) :: box_max = &
    (/ 10.0D+00, 10.0D+00 /)
  real      ( kind = 8 ) cell_centroid(dim_num,cell_num)
  integer   ( kind = 4 ) cell_fixed(cell_num)
  real      ( kind = 8 ) cell_generator(dim_num,cell_num)
  integer   ( kind = 4 ), parameter :: cell_sample = 5000
  real      ( kind = 8 ), dimension ( cell_num ) :: cell_volume
  real      ( kind = 8 ), dimension ( cell_num ) :: cell_weight
  integer   ( kind = 4 ) centroid_it
  integer   ( kind = 4 ), parameter :: centroid_it_max = 10
  real      ( kind = 8 ), parameter :: cutoff_dist = 4.0D+00
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) :: ppmb_ncol = 500
  integer   ( kind = 4 ) :: ppmb_nrow = 500
  logical, parameter :: region_plot = .true.
  character ( len = 80 ) :: region_plot_file_name = 'cvt_fixed_000.ppm'
  real      ( kind = 8 ), parameter :: region_volume = 100.0D+00
  integer   ( kind = 4 ) sample_num
  integer   ( kind = 4 ) :: seed = 123454321
!
!  Print introduction and options.
!
  call timestamp ( )

  sample_num = cell_num * cell_sample

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_FIXED_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CVT_FIXED library.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate a CVT algorithm in which some'
  write ( *, '(a)' ) '  generators are held fixed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Geometry parameters:'
  write ( *, '(a)' ) '-------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i3)' ) '  The spatial dimension is DIM_NUM = ', dim_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The minimum corner of the bounding box is:'
  write ( *, '(2f10.4)' ) box_min(1:dim_num)
  write ( *, '(a)' ) '  The maximum corner of the bounding box is:'
  write ( *, '(2f10.4)' ) box_max(1:dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT Algorithm parameters:'
  write ( *, '(a)' ) '-------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of Voronoi cells to generate: ', cell_num
  write ( *, '(a,i8)' ) '  Number of sampling points: ', sample_num
  write ( *, '(a,i8)' ) '  Number of centroid iterations: ', centroid_it_max
  write ( *, '(a,g14.6)' ) '  The influence cutoff distance: ', cutoff_dist
  write ( *, '(a)' ) &
    '  Voronoi cell generators are initialized by RANDOM_NUMBER.'

  write ( *, '(a)' ) ' '
!
!  Initialize the random number generator.
!
  call random_initialize ( seed )
!
!  Set the cells which are going to be fixed or free.
!
  cell_fixed(1:4) = 1
  cell_fixed(5:cell_num) = 0
!
!  Set the cell weights.
!
  cell_weight(1:cell_num) = 1.0D+00
!
!  Initialize the Voronoi cell generators.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initializing the cell generators.'

  call generator_init ( dim_num, box_min, box_max, cell_num, cell_generator )
!
!  Reset some points.
!
  cell_generator(1:2,1) = (/ 6.0D+00, 1.0D+00 /)
  cell_generator(1:2,2) = (/ 6.0D+00, 3.0D+00 /)
  cell_generator(1:2,3) = (/ 6.0D+00, 5.0D+00 /)
  cell_generator(1:2,4) = (/ 4.0D+00, 5.0D+00 /)

  cell_generator(1,5:cell_num) = 8.0D+00
!
!  Carry out the CVT iteration, which drives the Voronoi generators
!  and Voronoi centroids closer and closer.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Carry out the CVT iteration.'

  do centroid_it = 1, centroid_it_max

    if ( 1 < centroid_it ) then

      do j = 1, cell_num

        if ( cell_fixed(j) /= 1 ) then
          cell_generator(1:dim_num,j) = cell_centroid(1:dim_num,j)
        end if

      end do

    end if

    call cvt_iteration ( dim_num, box_min, box_max, cell_num, &
      cell_generator, sample_num, cell_centroid, cell_volume, &
      cutoff_dist, region_volume )

    if ( region_plot ) then

      call file_name_inc ( region_plot_file_name )

      call region_plot_ppmb ( region_plot_file_name, ppmb_nrow, ppmb_ncol, &
        cell_num, dim_num, box_min, box_max, cell_generator, &
        cell_centroid, cell_weight, cutoff_dist )

    end if

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_FIXED_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test_region ( x, dim_num, ival )

!*****************************************************************************80
!
!! TEST_REGION determines if a point is within the physical region.
!
!  Discussion:
!
!    Using a simple routine like this is only appropriate for a simple
!    region that can be easily defined by user formulas.
!
!    Computation of the "on-the-boundary" case is not considered important.
!    Only "inside" or "outside" is essential.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the point to be checked.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Output, integer ( kind = 4 ) IVAL, indicates the status of the point:
!    -1: the point is on the boundary of the region.
!     0: the point is outside the region.
!    +1: the point is inside the region.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) ival
  real    ( kind = 8 ) x(dim_num)

  if ( 0.0D+00 <= x(1) .and. x(1) <= 10.0D+00 .and. &
       0.0D+00 <= x(2) .and. x(2) <= 10.0D+00 ) then
    ival = 1
  else
    ival = 0
  end if

  return
end
