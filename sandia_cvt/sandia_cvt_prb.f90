program main

!*****************************************************************************80
!
!! MAIN is the main program for CVT_MAIN.
!
!  Discussion:
!
!    CVT_MAIN runs a test of CVT, the centroidal Voronoi tessellation library.
!
!    Typical parameters for this code are:
!
!      n = 1024
!      nbin = (/ 25, 25, 5 /)
!      ns_cvt = 5000
!      ns_mom = 5000
!
!    A "big" version of this code was run with:
!
!      n = 5000
!      nbin = (/ 25, 25, 5 /)
!      ns_cvt = 100
!      ns_mom = 100
!
!    A "huge" version of this code was run with:
!
!      n = 50000
!      nbin = (/ 50, 50, 10 /)
!      ns_cvt = 100
!      ns_mom = 100
!
!    A "monster" version of this code was run with:
!
!      n = 500000
!      nbin = (/ 100, 100, 20 /)
!      ns_cvt = 100
!      ns_mom = 100
!
!  Variables:
!
!    Note that floating point variables are declared as
!    of type "real ( kind = 8 )".  This is meant to be equivalent to
!    "double precision".  If your compiler does not recognize
!    the first kind of declaration, the second should work.
!
!    Geometry Variables:
!    ------------------
!
!    Input, integer NDIM, the spatial dimension, which should be 2 or 3.
!
!    real ( kind = 8 ) BOX_MIN(NDIM), BOX_MAX(NDIM), the minimum and maximum
!    coordinates of the bounding box.
!
!    logical USE_DIATOM, is TRUE if DIATOM is to be called to
!    determine whether a point lies in the physical region; if it is
!    FALSE than a much simplified routine is used.
!
!    real ( kind = 8 ) DR, a tolerance used by DIATOM when testing whether
!    a point is within, outside of, or on the boundary of the physical region.
!    This is typically a very small value, of the order of 0.00001 or smaller.
!
!
!    CVT Algorithm Variables:
!    -----------------------
!
!    integer N, the number of Voronoi cells to generate.
!    A typical value is 256.
!
!    integer MAXIT, the maximum number of correction iterations used in the
!    Voronoi calculation.  A typical value is 10.
!
!    integer NS_CVT, the average number of sampling points tested per
!    Voronoi cell, on each step of the correction iteration.  A typical
!    value is 5000.
!
!    Input, integer RANDOM_GENERATOR, specifies how the Voronoi cell
!    generators are to be initialized.
!    0, use the F90 RANDOM_NUMBER routine;
!    1, use the Halton sequence. (PREFERRED)
!
!    real ( kind = 8 ) CELL_GENERATOR(NDIM,N), the Voronoi cell generators
!    of the Voronoi tessellation, as approximated by the CVT algorithm.  This
!    is the output quantity of most interest.
!
!
!    Moment Calculation Variables:
!    ----------------------------
!
!    integer NS_MOM, the average number of sampling points tested per
!    Voronoi cell, for the moment calculation.  A typical value is 5000.
!
!    logical REGION_VOLUME_GIVEN,
!    0, the area or volume is already available in REGION.
!    nonzero, the area or volume needs to be calculated.
!
!    real ( kind = 8 ) REGION_VOLUME.
!    If REGION_VOLUME_GIVEN, then REGION_VOLUME must be input by the user.
!    Otherwise, REGION_VOLUME is approximated computationally.
!
!    real ( kind = 8 ) CELL_VOLUME(N), the volume of the Voronoi cells.
!
!    real ( kind = 8 ) CELL_CENTROID(NDIM,N), the centroids of the
!    Voronoi cells.
!
!    real ( kind = 8 ) CELL_MOMENT(NDIM,NDIM,N), the second moments of the
!    Voronoi cells.
!
!
!    The Nearest Neighbor Search Variables:
!    -------------------------------------
!
!    logical USE_BINS, is TRUE if the bounding box is to be divided up
!    into bins to speed up the nearest neighbor search;
!    FALSE if the nearest neighbor seach is to be done naively.
!
!    integer NBIN(3) is the number of bins to use in each direction.
!    For 2D problems, set NBIN(3) = 1.
!    For efficiency, these values should be set in such a way that the bins
!    are nearly square or cubical.
!
!    integer BIN_START(NBIN(1),NBIN(2),NBIN(3)), the index of the first
!    cell center in the bin, or -1 if none.
!
!    integer BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the index of the last
!    cell center in the bin, or -1 if none.
!
!    integer BIN_NEXT(N), the index of the next cell center in the bin
!    containing this cell center.
!
!  Miscellaneous Variables:
!  -----------------------
!
!    integer SEED, determines how to initialize the RANDOM_NUMBER routine.
!    If SEED is zero on input, then RANDOM_INITIALIZE will make up a seed
!    from the current real time clock reading.
!    If SEED is nonzero, then a reproducible sequence of random numbers
!    defined by SEED will be chosen.
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
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Reference:
!
!    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
!    User Manual and Supporting Information for Library of Codes
!    for Centroidal Voronoi Placement and Associated Zeroth,
!    First, and Second Moment Determination,
!    Sandia National Laboratories Technical Report SAND2002-0099,
!    February 2002.
!
  implicit none

  integer, parameter :: n = 1024
  integer, parameter :: ndim = 3
  integer, parameter, dimension ( ndim ) :: nbin = (/ 25, 25, 5 /)

  integer, dimension ( nbin(1),nbin(2),nbin(3) ) :: bin_last
  integer, dimension ( n ) :: bin_next
  integer, dimension ( nbin(1),nbin(2),nbin(3) ) :: bin_start
  real ( kind = 8 ), parameter, dimension ( ndim ) :: box_min = &
    (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( ndim ) :: box_max = &
    (/ 100.0D+00, 100.0D+00, 20.0D+00 /)
  real ( kind = 8 ) cell_centroid(ndim,n)
  real ( kind = 8 ) cell_generator(ndim,n)
  real ( kind = 8 ) cell_moment(ndim,ndim,n)
  real ( kind = 8 ) cell_volume(n)
  integer clock0
  integer clock1
  integer clock_max
  integer clock_rate
  real ( kind = 8 ),parameter :: dr = 0.00001D+00
  integer i
  integer ios
  integer it
  integer, parameter :: maxit = 10
  integer, parameter :: ns_cvt = 5000
  integer, parameter :: ns_mom = 5000
  logical :: quality_checks = .true.
  integer, parameter :: random_generator = 1
  real ( kind = 8 ), parameter :: region_volume = 20.0D+00 * 1700.0D+00
  logical, parameter :: region_volume_given = .true.
  integer :: seed = 0
  real tarray(2)
  real time0
  real time1
  integer, dimension ( n ) :: updates
  logical, parameter :: use_bins = .true.
  logical, parameter :: use_diatom = .false.
  logical, parameter :: write_output = .false.
!
!  Print introduction and options.
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SANDIA_CVT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SANDIA_CVT library.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A sample problem for the probabilistic'
  write ( *, '(a)' ) '  Centroidal Voronoi Tessellation algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Given a region in 2D or 3D, the problem is to determine'
  write ( *, '(a)' ) '  GENERATORS, a set of points which define a division'
  write ( *, '(a)' ) '  of the region into Voronoi cells, which are also'
  write ( *, '(a)' ) '  CENTROIDS of the Voronoi cells.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Geometry parameters:'
  write ( *, '(a)' ) '-------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension is NDIM = ', ndim
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The minimum corner of the bounding box is:'
  write ( *, '(5g14.6)' ) box_min(1:ndim)
  write ( *, '(a)' ) '  The maximum corner of the bounding box is:'
  write ( *, '(5g14.6)' ) box_max(1:ndim)
  if ( use_diatom ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DIATOM is called to determine the region.'
    write ( *, '(a,g14.6)' ) '  The DIATOM tolerance DR = ', dr
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DIATOM is not called;'
    write ( *, '(a)' ) '  a simple routine determines the region.'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT Algorithm parameters:'
  write ( *, '(a)' ) '-------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of Voronoi cells to generate: ', n
  write ( *, '(a,i8)' ) '  Number of iterations to determine CVT: ', maxit
  write ( *, '(a,i8)' ) '  Number of sampling points per Voronoi cell: ', ns_cvt

  if ( random_generator == 0 ) then
    write ( *, '(a)' ) &
      '  Voronoi cell generators are initialized by RANDOM_NUMBER.'
  else if ( random_generator == 1 ) then
    write ( *, '(a)' ) '  Voronoi cell generators are initialized by Halton.'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Moment parameters:'
  write ( *, '(a)' ) '------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of sampling points per Voronoi cell: ', ns_mom
  if ( region_volume_given ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The volume of the region is given.'
    write ( *, '(a,g14.6)' ) '  It is specified as REGION_VOLUME = ', &
      region_volume
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The volume of the region is not given.'
    write ( *, '(a)' ) '  It will be estimated.'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Nearest Neighbor Search parameters:'
  write ( *, '(a)' ) '-----------------------------------'
  write ( *, '(a)' ) ' '
  if ( use_bins ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  The nearest neighbor search is speeded up by using bins.'
    write ( *, '(a)' ) '  The bounding box is to be divided up into bins.'
    write ( *, '(a)' ) '  The number of bins is :'
    write ( *, '(5i8)' ) nbin(1:ndim)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The nearest neighbor search is not speeded up.'
    write ( *, '(a)' ) '  The nearest neighbor search is done by exhaustion.'
  end if

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
!  If DIATOM is to be used, the DIATOM setup routine must be called now.
!
  if ( use_diatom ) then
    call diatom_setup ( )
  end if
!
!  Begin timing.
!
!  The F90 CPU_TIME routine, at least on the DEC_ALPHA, can "wrap around"
!  after a relatively short time, giving negative timings.  So we also
!  call ETIME, which should be measuring the same thing, for checking.
!
!  The ETIME routine is called with an appended underscore, because
!  A) We're compiling with NOUNDERSCORE in order to interface with C;
!  B) but the UFOR library containing ETIME expects a reference to ETIME_.
!
!  The F90 SYSTEM_CLOCK routine should measure real execution time.
!
  call system_clock ( clock0, clock_rate, clock_max )
  call cpu_time ( time0 )
!
!  Initialize the Voronoi cell generators.
!
  call generator_init ( ndim, box_min, box_max, n, cell_generator, &
    use_diatom, dr, random_generator )
!
!  Carry out the CVT iteration, which drives the Voronoi generators
!  and Voronoi centroids closer and closer.
!
!  If bins are used, they must be calculated before the first call,
!  and recalculated each time the cell centers are changed.
!
  updates(1:n) = 1

  do it = 1, maxit + 1

    if ( use_bins ) then

      call bin_preprocess ( ndim, box_min, box_max, n, cell_generator, &
        nbin, bin_start, bin_last, bin_next )

    end if

    if ( it <= maxit ) then

      call cvt_iteration ( ndim, box_min, box_max, n, cell_generator, ns_cvt, &
        use_diatom, use_bins, dr, updates, nbin, bin_start, bin_last, bin_next )

    end if

  end do
!
!  Calculate cell moments.
!
  call vcm ( ndim, box_min, box_max, n, cell_generator, ns_mom, use_diatom, &
    use_bins, region_volume_given, region_volume, dr, nbin, &
    bin_start, bin_last, bin_next, cell_volume, cell_centroid, cell_moment )
!
!  End timing the "interesting" part.
!
  call cpu_time ( time1 )
  call system_clock ( clock1, clock_rate, clock_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a)' ) 'Elapsed CPU time, CPU_TIME: ', time1-time0, &
    ' seconds.'
  write ( *, '(a,g14.6,a)' ) 'Elapsed time, SYSTEM_CLOCK: ', &
    real ( clock1 - clock0 ) &
    / real ( clock_rate ), ' seconds.'
!
!  Compute quality checks.
!
  if ( quality_checks ) then
    call quality ( ndim, n, cell_moment, cell_volume, region_volume )
  end if
!
!  Write generators and moments to files.
!
  if ( write_output ) then

    open ( unit = 1, file = 'cvt_generators.txt', status = 'replace', &
      iostat = ios )
    do i = 1, n
      write ( 1, '(3d15.6)' ) cell_generator(1:ndim,i)
    end do
    close ( unit = 1 )

    open ( unit = 1, file = 'cvt_volume.txt', status = 'replace', iostat = ios )
    write ( 1, '(d15.6)' ) cell_volume(1:n)
    close ( unit = 1 )

    open ( unit = 1, file = 'cvt_centroid.txt', status = 'replace', &
      iostat = ios )
    do i = 1, n
      write ( 1,'(3d15.6)' ) cell_centroid(1:ndim,i)
    end do
    close ( unit = 1 )

    open ( unit = 1, file = 'cvt_moment.txt', status = 'replace', iostat = ios )
    do i = 1, n
      write ( 1, '(3d15.6)' ) cell_moment(1:ndim,1:ndim,i)
    end do
    close ( unit = 1 )

  end if
!
!  Terminate
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine diatom_setup ( )

!*****************************************************************************80
!
!! DIATOM_SETUP is a dummy version of the DIATOM_SETUP routine.
!
!  Discussion:
!
!    If the DIATOM library is not available, or not needed, then this
!    routine may be used to satisfy the compiler.
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
!    John Burkardt, Max Gunzburger, Janet Peterson
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DIATOM_SETUP:'
  write ( *, '(a)' ) '  Dummy version of DIATOM setup.'

  return
end
subroutine diatom_point_test2 ( x, y, z, dr, mdens, ival )

!*****************************************************************************80
!
!! DIATOM_POINT_TEST2 is a dummy version of the DIATOM_POINT_TEST2 routine.
!
!  Discussion:
!
!    If the DIATOM library is not available, then this routine may be used
!    to satisfy the compiler.
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
!    John Burkardt, Max Gunzburger, Janet Peterson
!
  implicit none

  real ( kind = 8 ) dr
  integer ival
  real ( kind = 8 ) mdens
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DIATOM_POINT_TEST2 - Fatal error!'
  write ( *, '(a)' ) '  The dummy version of DIATOM_POINT_TEST2 was called.'

  stop
end
subroutine test_region ( x, ndim, ival )

!*****************************************************************************80
!
!! TEST_REGION determines if a point is within the physical region.
!
!  Discussion:
!
!    The previous version of this routine had a flaw, because it
!    described the region as a rectangle sitting on top of half
!    an annulus.  Of course, such a rectangle only touches the
!    annulus at one point, and on either side of this contact there
!    is a small area, which should have been included in the overall
!    region, but was omitted.  This problem has been corrected.
!
!    Also, apparently, the region should extend no higher than
!    Y = 90, and no lower than Y = 0.  These constraints have
!    also been included now.
!
!    This routine should be functionally the same as the interface to
!    DIATOM, but it is much simpler, and a lot faster.  To signal that
!    the problem geometry will be determined by this routine, set
!    USE_DIATOM to FALSE in the calling program.
!
!    Using a simple routine like this is only appropriate for a simple
!    region that can be easily defined by user formulas.  This version of
!    the routine is a demonstration that implements the 2D and 3D versions
!    of the "tuning fork" test region.
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
!    12 December 2003
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, double precision X(NDIM), the point to be checked.
!
!    Input, integer NDIM, the dimension of the space.
!
!    Output, integer IVAL, indicates the status of the point:
!    -1: the point is on the boundary of the region.
!     0: the point is outside the region.
!    +1: the point is inside the region.
!
  implicit none

  integer ndim

  real ( kind = 8 ) c
  integer ival
  real ( kind = 8 ) x(ndim)
!
!  The only adjustment needed for 3D is that we require 0 <= Z <= 20.
!
  if ( ndim == 3 ) then
    if ( x(3) < 0.0D+00 .or. 20.0D+00 < x(3) ) then
      ival = 0
      return
    end if
  end if
!
!  Is the point inside the rectangular region?
!
  if ( 30.0D+00 <= x(2)             .and. &
                   x(2) <= 90.0D+00 .and. &
       45.0D+00 <= x(1)             .and. &
                   x(1) <= 55.0D+00 ) then

      ival = 1
!
!  Is the point in the annulus?
!
  else if ( 0.0D+00 <= x(2) .and. x(2) <= 40.0D+00 ) then

    c = sqrt ( ( x(1) - 50.0D+00 ) * ( x(1) - 50.0D+00 ) + x(2) * x(2) )

    if ( 30.0D+00 <= c .and. c <= 40.0D+00 ) then
      ival = 1
    else
      ival = 0
    end if

  else

    ival = 0

  end if

  return
end
