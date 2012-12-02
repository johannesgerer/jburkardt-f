program main

!*****************************************************************************80
!
!! MAIN is the main program for CVT_SIZE_PRB.
!
!  Discussion:
!
!    CVT_SIZE_PRB runs a test of CVT_SIZE, the sized CVT library.
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
!  Local Parameters:
!
!    Local, integer CELL_NUM, the number of generators to use.
!
!    Local, integer DIM_NUM, the spatial dimension (usually 2).
!
!    Local, integer PPMB_NCOL, the number of columns in the PPMB plots.
!
!    Local, integer PPMB_NROW, the number of rows in the PPMB plots.
!
!    Local, integer SEED, the seed for the random number generator.
!    Setting SEED = 0 will cause the program to choose a seed
!    for the random number generator using the current date and time.
!
  implicit none

  integer, parameter :: cell_num = 10
  integer, parameter :: dim_num = 2

  real ( kind = 8 ) area_disc
  real ( kind = 8 ), parameter, dimension ( dim_num ) :: box_min = &
    (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( dim_num ) :: box_max = &
    (/ 10.0D+00, 10.0D+00 /)
  real ( kind = 8 ) cell_centroid(dim_num,cell_num)
  real ( kind = 8 ) cell_generator(dim_num,cell_num)
  integer, parameter :: cell_sample = 5000
  real ( kind = 8 ) cell_volume(cell_num)
  real ( kind = 8 ) cell_volume_desired(cell_num)
  real ( kind = 8 ) cell_weight(cell_num)
  integer centroid_it
  integer, parameter :: centroid_it_max = 100
  integer i
  integer ios
  integer :: ppmb_ncol = 500
  integer :: ppmb_nrow = 500
  logical, parameter :: region_plot = .true.
  character ( len = 80 ) :: region_plot_file_name = 'cvt_000.ppmb'
  real ( kind = 8 ) region_volume
  integer sample_num
  logical, parameter :: volume_plot = .false.
  character ( len = 80 ) :: volume_file_name = 'vol_000.txt'
  integer :: seed = 123456789
  logical, parameter :: write_output = .true.
!
!  Print introduction and options.
!
  call timestamp ( )

  region_volume = 100.0D+00
!
!  Set the desired cell volumes.
!
  do i = 1, cell_num
    cell_weight(i) = real ( i, kind = 8 )
  end do
!
! do i = 1, cell_num
!   cell_weight(i) = real ( 4 * ( i - 1 ) + 1, kind = 8 )
! end do
!
! do i = 1, cell_num
!   cell_weight(i) = 2**( i - 1 )
! end do
!
! do i = 1, cell_num
!   if ( i <= cell_num / 2 ) then
!     cell_weight(i) = 1.0D+00
!   else
!     cell_weight(i) = 3.0D+00
!   end if
! end do
!
! do i = 1, cell_num
!   if ( i <= 3 * cell_num / 4 ) then
!     cell_weight(i) = 1.0D+00
!   else
!     cell_weight(i) = 4.0D+00
!   end if
! end do
!
! do i = 1, cell_num
!   if ( i <= cell_num / 2 ) then
!     cell_weight(i) = 1.0D+00
!   else
!     cell_weight(i) = 4.0D+00
!   end if
! end do

  cell_weight(1:cell_num) = cell_weight(1:cell_num) &
    / sum ( cell_weight(1:cell_num) )
!
!  We assume that the initial cell weights are in fact the desired
!  relative areas, and use this to set the desired absolute areas.
!
  cell_volume_desired(1:cell_num) = region_volume * cell_weight(1:cell_num)
!
!  Now we take the DIM_NUM-th root of CELL_WEIGHT to get a quantity
!  appropriate for treating linear distances.
!
  cell_weight(1:cell_num) = &
    cell_weight(1:cell_num)**( 1.0D+00 / real ( dim_num, kind = 8 ) )
  cell_weight(1:cell_num) = cell_weight(1:cell_num) &
    / sum ( cell_weight(1:cell_num) )

  sample_num = cell_num * cell_sample

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_SIZE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CVT_SIZE library.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A sample problem for the probabilistic'
  write ( *, '(a)' ) '  sized Centroidal Voronoi Tessellation algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Given a region in 2D, the problem is to determine'
  write ( *, '(a)' ) '  GENERATORS, a set of points which define a division'
  write ( *, '(a)' ) '  of the region into Voronoid cells, which are also'
  write ( *, '(a)' ) '  CENTROIDS of the Voronoi cells, and which have a '
  write ( *, '(a)' ) '  certain SIZE.'

  write ( *, '(a)' ) ' '
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
  write ( *, '(a,i6)' ) '  The number of Voronoi cells to generate: ', cell_num
  write ( *, '(a,i6)' ) '  Number of sampling points: ', sample_num
  write ( *, '(a,i6)' ) '  Number of centroid iterations: ', centroid_it_max
  write ( *, '(a)' )    '  (This forces generator and centroid closer.)'

  write ( *, '(a)' ) &
    '  Voronoi cell generators are initialized by RANDOM_NUMBER.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The desired CVT cell volumes and initial weights are:'
  write ( *, '(a)' ) ' '
  do i = 1, cell_num
    write ( *, '(i2,2g14.6)' ) i, cell_volume_desired(i), cell_weight(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Nearest Neighbor Search parameters:'
  write ( *, '(a)' ) '-----------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The nearest neighbor search is not speeded up.'
  write ( *, '(a)' ) '  The nearest neighbor search is done by exhaustion.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Miscellaneous parameters:'
  write ( *, '(a)' ) '------------------------'
  write ( *, '(a)' ) ' '
  if ( write_output ) then
    write ( *, '(a)' ) '  Generator and moment output files WILL be written.'
  else
    write ( *, '(a)' ) &
      '  Generator and moment output files will NOT be written.'
  end if

  write ( *, '(a)' ) ' '
!
!  Initialize the random number generator.
!
!  If SEED is zero on input, then the routine will make up a seed.
!  If SEED is nonzero, then a reproducible sequence of random numbers
!  defined by SEED will be chosen.
!
  call random_initialize ( seed )
!
!  Initialize the Voronoi cell generators.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initializing the cell generators.'

  call generator_init ( dim_num, box_min, box_max, cell_num, cell_generator )
!
!  Carry out the CVT iteration, which drives the Voronoi generators
!  and Voronoi centroids closer and closer.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Carry out the CVT iteration.'

  do centroid_it = 1, centroid_it_max

    if ( centroid_it > 1 ) then
      cell_generator(1:dim_num,1:cell_num) = &
      cell_centroid(1:dim_num,1:cell_num)
    end if

    call cvt_size_iteration ( dim_num, box_min, box_max, cell_num, &
      cell_generator, cell_weight, sample_num, cell_centroid, &
      cell_volume, region_volume )

    if ( region_plot ) then

      call file_name_inc ( region_plot_file_name )

      call region_plot_ppmb ( region_plot_file_name, ppmb_nrow, ppmb_ncol, &
        cell_num, dim_num, box_min, box_max, cell_generator, &
        cell_centroid, cell_weight )

    end if

    if ( volume_plot ) then

      call file_name_inc ( volume_file_name )

      open ( unit = 1, file = volume_file_name, status = 'replace', &
        iostat = ios )
      write ( 1, '(d15.6)' ) cell_volume(1:cell_num)
      close ( unit = 1 )

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Desired, Actual, Linear'
  write ( *, '(a)' ) '  Cell  Volume   Volume  Weight'
  write ( *, '(a)' ) ' '
  do i = 1, cell_num
    write ( *, '(i6,3g14.6)' ) &
      i, cell_volume_desired(i), cell_volume(i), cell_weight(i)

  end do
!
!  Determine the discrepancies.
!
  area_disc = 0.0D+00
  do i = 1, cell_num
    area_disc = area_disc + abs ( cell_volume_desired(i) - cell_volume(i) )
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) 'Area discrepancy = ', area_disc
!
!  Print some stuff.
!
  if ( .true. ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Cell generators and centroids:'
    write ( *, '(a)') ' '
    do i = 1, cell_num
      write ( *, '(i4,3d15.6)' ) i, cell_generator(1:dim_num,i)
      write ( *, '(4x,3d15.6)' )    cell_centroid(1:dim_num,i)
    end do
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Cell volumes:'
    write ( *, '(a)' ) ' '
    do i = 1, cell_num
      write ( *, '(i4,3d15.6)' ) i, cell_volume(i)
    end do
  end if
!
!  Write generators, weights and moments to files.
!
  if ( write_output ) then

    open ( unit = 1, file = 'cvt_generators.txt', status = 'replace', &
      iostat = ios )
    do i = 1, cell_num
      write ( 1, '(3d15.6)' ) cell_generator(1:dim_num,i)
    end do
    close ( unit = 1 )

    open ( unit = 1, file = 'cvt_volume_desired.txt', status = 'replace', &
      iostat = ios )
    write ( 1, '(d15.6)' ) cell_volume_desired(1:cell_num)
    close ( unit = 1 )

  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_SIZE_PRB'
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
