program main

!*****************************************************************************80
!
!! MAIN is the main program for CVT_REFINE_PRB.
!
!  Discussion:
!
!    CVT_REFINE_PRB tests the CVT refinement option.
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
!    28 January 2003
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
!    Local, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,CELL_NUM), the Voronoi
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
!    the point is not assigned.  When weights are used, the
!    effect of cutoff_dist is exaggerated.  Just to avoid
!    using it at all, you have to be sure to make it rather
!    big.
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
!    is the product of CELL_NUM, the number of cells and
!    CELL_SAMPLE, the number of sample points per cell.
!
!    Local, integer SEED, the seed for the random number generator.
!    Setting SEED = 0 will cause the program to choose a seed
!    for the random number generator using the current date and time.
!
  implicit none

  integer, parameter :: cell_num1 = 50
  integer, parameter :: cell_inc = 50
  integer, parameter :: dim_num = 2

  integer, parameter :: cell_num2 = cell_num1 + cell_inc

  real ( kind = 8 ), parameter, dimension ( dim_num ) :: box_min = &
    (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( dim_num ) :: box_max = &
    (/ 10.0D+00, 10.0D+00 /)
  real ( kind = 8 ) cell_centroid(dim_num,cell_num2)
  integer cell_fixed(cell_num2)
  real ( kind = 8 ) cell_generator(dim_num,cell_num2)
  integer cell_num
  integer, parameter :: cell_sample = 5000
  real ( kind = 8 ), dimension ( cell_num2 ) :: cell_size
  real ( kind = 8 ), dimension ( cell_num2 ) :: cell_volume
  real ( kind = 8 ), dimension ( cell_num2 ) :: cell_weight
  logical, parameter :: centroid_file_write = .true.
  integer centroid_it
  integer, parameter :: centroid_it_max = 100
!
!  The interpretation of CUTOFF_DIST is problematic
!  when weights are being used.  If you don't want it in
!  effect, make it ridiculously big for now.
!
  real ( kind = 8 ), parameter :: cutoff_dist = 10000.0D+00
  logical, parameter :: generator_file_write = .true.
  integer i
  integer ios
  integer iounit
  integer j
  integer :: ppmb_ncol = 500
  integer :: ppmb_nrow = 500
  logical, parameter :: region_plot = .false.
  character ( len = 80 ) :: region_plot_file_name = 'cvt_refine_000.ppmb'
  real ( kind = 8 ), parameter :: region_volume = 100.0D+00
  integer sample_num
  integer :: seed = 123456789
!
!  Print introduction and options.
!
  call timestamp ( )
!
!  Initialize the random number generator.
!
  call random_initialize ( seed )

  cell_num = cell_num1
  sample_num = cell_num * cell_sample

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_REFINE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CVT_REFINE library.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate a CVT algorithm in which a CVT'
  write ( *, '(a)' ) '  is refined.'
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
  write ( *, '(a)' ) ' '
!
!  Set the cells which are going to be fixed or free.
!
  cell_fixed(1:cell_num) = 0
!
!  Set the cell weights.
!
  cell_size(1:cell_num) = 1.0D+00

  call cvt_size_to_weight ( dim_num, cell_num, cell_size, cell_weight )

  call r8vec_print ( cell_num, cell_weight, '  Initial cell weights:' )
!
!  Initialize the Voronoi cell generators.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initializing the cell generators.'

  call generator_init ( dim_num, box_min, box_max, cell_num, &
    cell_generator, cell_fixed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Initial generators:'
  write ( *, '(a)' ) ' '
  do i = 1, cell_num
    write ( *, '(i3,2g14.6)' ) i, cell_generator(1:dim_num,i)
  end do
!
!  Carry out the CVT iteration, which drives the Voronoi generators
!  and Voronoi centroids closer and closer.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Carry out the first CVT iteration.'
  write ( *, '(a)' ) ' '

  do centroid_it = 1, centroid_it_max

    if ( 1 < centroid_it ) then

      do j = 1, cell_num

        if ( cell_fixed(j) /= 1 ) then
          cell_generator(1:dim_num,j) = cell_centroid(1:dim_num,j)
        end if

      end do

    end if

    call cvt_iteration_weight ( dim_num, box_min, box_max, cell_num, &
      cell_generator, cell_weight, sample_num, cell_centroid, cell_volume, &
      cutoff_dist, region_volume )

    if ( region_plot ) then

      call file_name_inc ( region_plot_file_name )

      call region_plot_ppmb ( region_plot_file_name, ppmb_nrow, ppmb_ncol, &
        cell_num, dim_num, box_min, box_max, cell_generator, &
        cell_centroid, cell_weight, cutoff_dist )

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Final generators for first run:'
  write ( *, '(a)' ) ' '
  do i = 1, cell_num
    write ( *, '(i3,2g14.6)' ) i, cell_generator(1:dim_num,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_REFINE_PRB'
  write ( *, '(a)' ) '  End of first CVT calculation.'

  if ( .false. ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_REFINE_PRB'
    write ( *, '(a)' ) '  Normal end of execution.'

    write ( *, '(a)' ) ' '
    call timestamp ( )
    stop

  end if
!
!  IF REQUESTED, CARRY OUT A SECOND, REFINING, ITERATION.
!
  cell_num = cell_num2
  sample_num = cell_num * cell_sample

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_REFINE_PRB'
  write ( *, '(a)' ) '  Begin second calculation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT Algorithm parameters:'
  write ( *, '(a)' ) '-------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of Voronoi cells to generate: ', cell_num2
  write ( *, '(a,i8)' ) '  Number of sampling points: ', sample_num
  write ( *, '(a,i8)' ) '  Number of centroid iterations: ', centroid_it_max
  write ( *, '(a,g14.6)' ) '  The influence cutoff distance: ', cutoff_dist
  write ( *, '(a)' ) ' '
!
!  Set the cells which are going to be fixed or free.
!
  cell_fixed(1:cell_num1) = 1
  cell_fixed(cell_num1+1:cell_num2) = 0
!
!  Set the cell weights, causing the old points to shrink a lot,
!  allowing the new points to find equilibrium.
!
  cell_size(1:cell_num1) = 1.0D+00
  cell_size(cell_num1+1:cell_num2) = 10.0D+00

  call cvt_size_to_weight ( dim_num, cell_num, cell_size, cell_weight )

  call r8vec_print ( cell_num, cell_weight, '  Initial cell weights:' )
!
!  Initialize the Voronoi cell generators.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initializing the (free) cell generators.'

  call generator_init ( dim_num, box_min, box_max, cell_num, &
    cell_generator, cell_fixed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Initial generators:'
  write ( *, '(a)' ) ' '
  do i = 1, cell_num
    write ( *, '(i3,2g14.6)' ) i, cell_generator(1:dim_num,i)
  end do
!
!  Carry out the CVT iteration, which drives the Voronoi generators
!  and Voronoi centroids closer and closer.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Carry out the second CVT iteration.'
  write ( *, '(a)' ) ' '

  do centroid_it = centroid_it_max+1, centroid_it_max+centroid_it_max
!
!  Adjust the weights at certain times to that at the end
!  the weights are equal.
!
    if ( centroid_it == 125 ) then

      cell_size(1:cell_num1) = 1.0D+00
      cell_size(cell_num1+1:cell_num2) = 5.0D+00

      call cvt_size_to_weight ( dim_num, cell_num, cell_size, cell_weight )

    else if ( centroid_it == 150 ) then

      cell_size(1:cell_num1) = 1.0D+00
      cell_size(cell_num1+1:cell_num2) = 2.0D+00

      call cvt_size_to_weight ( dim_num, cell_num, cell_size, cell_weight )

    else if ( centroid_it == 175 ) then

      cell_size(1:cell_num1) = 1.0D+00
      cell_size(cell_num1+1:cell_num2) = 1.0D+00

      call cvt_size_to_weight ( dim_num, cell_num, cell_size, cell_weight )

    end if
!
!  Update the free generators.
!
    if ( centroid_it > centroid_it_max+1 ) then

      do j = 1, cell_num

        if ( cell_fixed(j) /= 1 ) then
          cell_generator(1:dim_num,j) = cell_centroid(1:dim_num,j)
        end if

      end do

    end if
!
!  Compute the centroids.
!
    call cvt_iteration_weight ( dim_num, box_min, box_max, cell_num, &
      cell_generator, cell_weight, sample_num, cell_centroid, cell_volume, &
      cutoff_dist, region_volume )
!
!  Make a plot.
!
    if ( region_plot ) then

      call file_name_inc ( region_plot_file_name )

      call region_plot_ppmb ( region_plot_file_name, ppmb_nrow, ppmb_ncol, &
        cell_num, dim_num, box_min, box_max, cell_generator, &
        cell_centroid, cell_weight, cutoff_dist )

    end if

  end do
!
!  Write generators and centroids to files.
!
  if ( generator_file_write ) then

    call get_unit ( iounit )

    open ( unit = iounit, file = 'cvt_generators.txt', status = 'replace', &
      iostat = ios )
    do i = 1, cell_num
      write ( iounit, '(3g15.6)' ) cell_generator(1:dim_num,i)
    end do
    close ( unit = iounit )

  end if

  if ( centroid_file_write ) then

    open ( unit = iounit, file = 'cvt_centroids.txt', status = 'replace', &
      iostat = ios )
    do i = 1, cell_num
      write ( iounit, '(3g15.6)' ) cell_centroid(1:dim_num,i)
    end do
    close ( unit = iounit )

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_REFINE_PRB'
  write ( *, '(a)' ) '  End of second CVT calculation.'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_REFINE_PRB'
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
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Output, integer IVAL, indicates the status of the point:
!    -1: the point is on the boundary of the region.
!     0: the point is outside the region.
!    +1: the point is inside the region.
!
  implicit none

  integer dim_num

  integer ival
  real ( kind = 8 ) x(dim_num)

  if ( 0.0D+00 <= x(1) .and. x(1) <= 10.0D+00 .and. &
       0.0D+00 <= x(2) .and. x(2) <= 10.0D+00 ) then
    ival = 1
  else
    ival = 0
  end if

  return
end
