subroutine cell_volume_computation ( dim_num, box_min, box_max, cell_num, &
  cell_gen_coord, sample_num, cell_volume, region_volume )

!*****************************************************************************80
!
!! CELL_VOLUME_COMPUTATION estimates the cell volumes by sampling.
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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real BOX_MIN(DIM_NUM), BOX_MAX(DIM_NUM), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Input, integer CELL_NUM, the number of Voronoi cells.
!
!    Input/output, real CELL_GEN_COORD(DIM_NUM,CELL_NUM), the Voronoi
!    cell generators.  On output, these have been modified.
!
  implicit none

  integer cell_num
  integer dim_num
  integer sample_num

  real box_max(1:dim_num)
  real box_min(1:dim_num)
  integer cell_count(cell_num)
  real cell_gen_coord(dim_num,cell_num)
  real cell_volume(cell_num)
  logical, parameter :: debug = .false.
  integer i
  integer j
  integer nearest(sample_num)
  integer ngen
  real region_volume
  real sample_coord(dim_num,sample_num)
!
!  Generate the sampling points.
!
  call region_sampler ( dim_num, sample_num, box_min, box_max, &
    sample_coord, ngen )
!
!  Find closest generators.
!
  call find_closest ( dim_num, sample_num, sample_coord, cell_num, &
    cell_gen_coord, nearest )

  cell_count(1:cell_num) = 0

  do i = 1, sample_num
    j = nearest(i)
    if ( j < 1 .or. j > cell_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CELL_VOLUME_COMPUTATION - Fatal error!'
      write ( *, '(a,i6)' ) '  i = ', i
      write ( *, '(a,i6)' ) '  j = ', j
      write ( *, '(a,i6)' ) '  Cell_num = ', cell_num
      stop
    end if
    cell_count(j) = cell_count(j) + 1
  end do

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Volume cell counts'
    write ( *, '(a)' ) ' '
    do j = 1, cell_num
      write ( *, '(2i8)' ) j, cell_count(j)
    end do
  end if

  do j = 1, cell_num
    cell_volume(j) = region_volume * real ( cell_count(j) ) &
      / real ( sample_num )
  end do

  return
end
subroutine cvt_dense_iteration ( density_num, dim_num, box_min, box_max, &
  cell_num, cell_gen_coord, density_value, density_coord, sample_num )

!*****************************************************************************80
!
!! CVT_DENSE_ITERATION takes one step of the CVT density iteration.
!
!  Discussion:
!
!    The routine is given a set of points, called "generators", which
!    define a tessellation of the region into Voronoi cells.  Each point
!    defines a cell.  Each cell, in turn, has a centroid, but it is
!    unlikely that the centroid and the generator coincide.
!
!    Each time this CVT iteration is carried out, an attempt is made
!    to modify the generators in such a way that they are closer and
!    closer to being the centroids of the Voronoi cells they generate.
!
!    A large number of sample points are generated, and the nearest generator
!    is determined.  A count is kept of how many points were nearest to each
!    generator.  Once the sampling is completed, the location of all the
!    generators is adjusted.  This step should decrease the discrepancy
!    between the generators and the centroids.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real BOX_MIN(DIM_NUM), BOX_MAX(DIM_NUM), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Input, integer CELL_NUM, the number of Voronoi cells.
!
!    Input/output, real CELL_GEN_COORD(DIM_NUM,CELL_NUM), the Voronoi
!    cell generators.  On output, these have been modified
!
!    Input, real DENSITY_VALUE(CELL_NUM), the density values.
!
!    Input, real DENSITY_COORD(DIM_NUM,CELL_NUM), the coordinates of the
!    density definition points.
!
!    Input, integer SAMPLE_NUM, the number of sample points.
!
  implicit none

  integer cell_num
  integer density_num
  integer dim_num
  integer sample_num

  real box_max(1:dim_num)
  real box_min(1:dim_num)
  integer cell
  integer cell_count(cell_num)
  real cell_gen_coord(dim_num,cell_num)
  logical, parameter :: debug = .true.
  real density_coord(dim_num,density_num)
  real density_value(density_num)
  integer i
  integer nearest(sample_num)
  real sample_coord(dim_num,sample_num)
!
!  Generate sampling points using the current density function.
!
  call region_sampler_weighted ( cell_num, density_num, dim_num, &
    sample_num, box_min, box_max, cell_gen_coord, density_coord, &
    density_value, sample_coord )
!
!  For each sample point, find the nearest cell generator.
!
  call find_closest ( dim_num, sample_num, sample_coord, cell_num, &
    cell_gen_coord, nearest )
!
!  Add each sampling point to the averaging data for its nearest generator.
!
  cell_count(1:cell_num) = 0
  cell_gen_coord(1:dim_num,1:cell_num) = 0.0E+00

  do i = 1, sample_num

    cell = nearest(i)

    if ( cell < 1 .or. cell > cell_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DENSE_ITERATION - Fatal error!'
      write ( *, '(a,i6)' ) '  Illegal value of CELL = ', cell
      write ( *, '(a,i6)' ) '  Sample point number is I = ', i
      write ( *, '(a)' ) '  Point is '
      write ( *, '(3g14.6)' ) sample_coord(1:dim_num,i)
      stop
    end if

    cell_gen_coord(1:dim_num,cell) = cell_gen_coord(1:dim_num,cell) &
      + sample_coord(1:dim_num,i)

    cell_count(cell) = cell_count(cell) + 1

  end do

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT cell counts'
    write ( *, '(a)' ) ' '
    do cell = 1, cell_num
      write ( *, '(2i8)' ) cell, cell_count(cell)
    end do
  end if
!
!  Compute the new generators, which are the centroids of their cells.
!
!  For nonuniform densities, the centroid is not guaranteed to be
!  at the center of the cell.
!
  do cell = 1, cell_num

    if ( cell_count(cell) /= 0 ) then

      cell_gen_coord(1:dim_num,cell) = cell_gen_coord(1:dim_num,cell) &
        / real ( cell_count(cell) )

    end if

  end do

  return
end
subroutine find_closest ( dim_num, sample_num, sample_coord, cell_num, &
  cell_gen_coord, nearest )

!*****************************************************************************80
!
!! FIND_CLOSEST finds the Voronoi cell generator closest to a point.
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
!    18 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer SAMPLE_NUM, the number of sample points.
!
!    Input, real SAMPLE_COORD(DIM_NUM,SAMPLE_NUM), the points to be checked.
!
!    Input, integer CELL_NUM, the number of cell generatorrs.
!
!    Input, real CELL_GEN_COORD(DIM_NUM,CELL_NUM), the cell generators.
!
!    Output, integer NEAREST(SAMPLE_NUM), the index of the nearest
!    cell generators.
!
  implicit none

  integer cell_num
  integer dim_num
  integer sample_num

  integer cell
  real cell_gen_coord(dim_num,cell_num)
  real distance
  real dist_sq
  integer nearest(sample_num)
  integer sample
  real sample_coord(dim_num,sample_num)

  do sample = 1, sample_num

    nearest(sample) = 0
    distance = huge ( distance )

    do cell = 1, cell_num

      dist_sq = sum ( &
        ( cell_gen_coord(1:dim_num,cell) - sample_coord(1:dim_num,sample) )**2 )

      if ( dist_sq < distance ) then
        distance = dist_sq
        nearest(sample) = cell
      end if

    end do

    if ( nearest(sample) == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FIND_CLOSEST - Fatal error!'
      write ( *, '(a,i6)' ) '  Sample = ', sample
      write ( *, '(a)' ) '  SAMPLE_COORD= '
      write ( *, '(3g14.6)' ) sample_coord(1:dim_num,sample)
      stop
    end if

  end do

  return
end
subroutine random_initialize ( seed )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN 90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed ( )
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator.
!
  implicit none

  integer count
  integer count_max
  integer count_rate
  integer i
  integer seed
  integer, allocatable :: seed_vector(:)
  integer seed_size
  real t
!
!  Initialize the random number seed.
!
  call random_seed ( )
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )

  if ( seed /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i12)' ) '  Initialize RANDOM_NUMBER with user SEED = ', seed

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i12)' ) &
      '  Initialize RANDOM_NUMBER with arbitrary SEED = ', seed

  end if
!
!  Now set the seed.
!
  seed_vector(1:seed_size) = seed

  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times.
!
  do i = 1, 100
    call random_number ( harvest = t )
  end do

  return
end
subroutine region_sampler ( dim_num, sample_num, box_min, box_max, &
  sample_coord, ngen )

!*****************************************************************************80
!
!! REGION_SAMPLER returns sample points in the physical region.
!
!  Discussion:
!
!    The calculations are done in DIM_NUM dimensional space.
!
!    The physical region is enclosed in a bounding box.
!
!    A point is chosen in the bounding box by a uniform random
!    number generator.
!
!    If a user-supplied routine determines that this point is
!    within the physical region, this routine returns.  Otherwise,
!    a new random point is chosen.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer SAMPLE_NUM, the number of points to be generated.
!
!    Input, real BOX_MIN(DIM_NUM), BOX_MAX(DIM_NUM), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Output, real SAMPLE_COORD(DIM_NUM,SAMPLE_NUM), the sample points.
!
!    Output, integer NGEN, the number of points that were generated.
!    This is at least SAMPLE_NUM, but may be larger if some points
!    were rejected.
!
  implicit none

  integer dim_num
  integer sample_num

  real box_max(dim_num)
  real box_min(dim_num)
  integer ival
  integer j
  integer ngen
  real r(dim_num)
  real sample_coord(dim_num,sample_num)

  ngen = 0

  do j = 1, sample_num

    do

      ngen = ngen + 1

      if ( ngen > 10000 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'REGION_SAMPLER - Fatal error!'
        write ( *, '(a,i6)' ) '  Number of rejected points = ', ngen
        write ( *, '(a)' ) &
          '  There may be a problem with the geometry definition.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Current random value is:'
        write ( *, '(3g14.6)' ) r(1:dim_num)
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Current sample point is:'
        write ( *, '(3g14.6)' ) sample_coord(1:dim_num,j)
        stop
      end if
!
!  Generate a point at random.
!
      call random_number ( r(1:dim_num) )
!
!  Determine a point in the bounding box.
!
      sample_coord(1:dim_num,j) = &
        ( ( 1.0E+00 - r(1:dim_num) ) * box_min(1:dim_num) &
                    + r(1:dim_num)   * box_max(1:dim_num) )

      call test_region ( sample_coord(1:dim_num,j), dim_num, ival )

      if ( ival == 1 ) then
        exit
      end if

    end do

  end do

  return
end
subroutine region_sampler_weighted ( cell_num, density_num, dim_num, &
  sample_num, box_min, box_max, cell_gen_coord, density_coord, &
  density_value, sample_coord )

!*****************************************************************************80
!
!! REGION_SAMPLER_WEIGHTED returns weighted sample points
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer CELL_NUM, the number of Voronoi cells.
!
!    Input, integer DENSITY_NUM, the number of density data points.
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer SAMPLE_NUM, the number of sample points to compute.
!
!    Input, real BOX_MIN(DIM_NUM), BOX_MAX(DIM_NUM), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Input, real CELL_GEN_COORD(DIM_NUM,CELL_NUM), the coordinates of the
!    Voronoi cell generator points.
!
!    Input, real DENSITY_VALUE(DENSITY_NUM), the density values.
!
!    Input, real DENSITY_COORD(DIM_NUM,DENSITY_NUM), the coordinates of the
!    density definition points.
!
  implicit none

  integer cell_num
  integer density_num
  integer dim_num
  integer sample_num

  real box_max(dim_num)
  real box_min(dim_num)
  integer cell
  integer cell_count(cell_num)
  real cell_gen_coord(dim_num,cell_num)
  logical, parameter :: debug = .false.
  real density_coord(dim_num,density_num)
  real density_value(density_num)
  integer i
  integer nearest(1)
  integer ngen
  real r
  integer reject_num
  integer sample
  real sample_coord(dim_num,sample_num)
  real sample_dense

  reject_num = 0
!
!  Generate sample points according to the current density.
!
  do sample = 1, sample_num

    do
!
!  Pick a uniformly random point in the region.
!
      call region_sampler ( dim_num, 1, box_min, box_max, &
        sample_coord(1,sample), ngen )
!
!  Compute its density.
!
      if ( .true. ) then
        call spline_constant_val ( density_num, density_coord, density_value, &
          sample_coord(1,sample), sample_dense )
      else
        call spline_linear_val ( density_num, density_coord, density_value, &
          sample_coord(1,sample), sample_dense )
      end if
!
!  Accept the point only if a random number is less than or equal to
!  the density (which is normalized, so that the maximum value is 1).
!
      call random_number ( harvest = r )

      if ( r <= sample_dense ) then
        call find_closest ( dim_num, 1, sample_coord(1,sample), cell_num, &
          cell_gen_coord, nearest )
        cell = nearest(1)
        exit
      else
        reject_num = reject_num + 1
      end if

    end do

  end do

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REGION_SAMPLER_WEIGHTED:'
    write ( *, '(a,i8)' ) '  Reject = ', reject_num
    write ( *, '(a,i8)' ) '  Accept = ', sample_num
  end if

  return
end
subroutine rvec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! RVEC_BRACKET searches a sorted array for successive brackets of a value.
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
!    Input, integer N, length of input array.
!
!    Input, real X(N), an array sorted into ascending order.
!
!    Input, real XVAL, a value to be bracketed.
!
!    Output, integer LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      XVAL > X(N), when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none

  integer n

  integer i
  integer left
  integer right
  real x(n)
  real xval

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
subroutine rvec_even ( alo, ahi, n, a )

!*****************************************************************************80
!
!! RVEC_EVEN returns N real values, evenly spaced between ALO and AHI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the low and high values.
!
!    Input, integer N, the number of values.
!
!    Output, real A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer n

  real a(n)
  real ahi
  real alo
  integer i

  if ( n == 1 ) then

    a(1) = 0.5E+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i ) * alo + real ( i - 1 ) * ahi ) / real ( n - 1 )
    end do

  end if

  return
end
subroutine spline_constant_val ( ndata, tdata, ydata, tval, yval )

!*****************************************************************************80
!
!! SPLINE_CONSTANT_VAL evaluates a piecewise constant spline at a point.
!
!  Discussion:
!
!    NDATA-1 points TDATA define NDATA intervals, with the first
!    and last being semi-infinite.
!
!    The value of the spline anywhere in interval I is YDATA(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points defining the spline.
!
!    Input, real TDATA(NDATA-1), the breakpoints.  The values of TDATA should
!    be distinct and increasing.
!
!    Input, YDATA(NDATA), the values of the spline in the intervals
!    defined by the breakpoints.
!
!    Input, real TVAL, the point at which the spline is to be evaluated.
!
!    Output, real YVAL, the value of the spline at TVAL.
!
  implicit none

  integer ndata

  integer i
  real tdata(ndata-1)
  real tval
  real ydata(ndata)
  real yval

  do i = 1, ndata-1
    if ( tval <= tdata(i) ) then
      yval = ydata(i)
      return
    end if
  end do

  yval = ydata(ndata)

  return
end
subroutine spline_linear_val ( ndata, tdata, ydata, tval, yval )

!*****************************************************************************80
!
!! SPLINE_LINEAR_VAL evaluates a linear spline at a specific point.
!
!  Discussion:
!
!    Because of the extremely simple form of the linear spline,
!    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
!    evaluate the spline at any point.  No processing of the data
!    is required.
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
!    Input, integer NDATA, the number of data points defining the spline.
!
!    Input, real TDATA(NDATA), YDATA(NDATA), the values of the independent
!    and dependent variables at the data points.  The values of TDATA should
!    be distinct and increasing.
!
!    Input, real TVAL, the point at which the spline is to be evaluated.
!
!    Output, real YVAL, the value of the spline at TVAL.
!
  implicit none

  integer ndata

  integer left
  integer right
  real tdata(ndata)
  real tval
  real ydata(ndata)
  real ypval
  real yval
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  Now evaluate the piecewise linear function.
!
  ypval = ( ydata(right) - ydata(left) ) / ( tdata(right) - tdata(left) )

  yval = ydata(left) +  ( tval - tdata(left) ) * ypval

  return
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
!  Parameters:
!
!    Input, real X(DIM_NUM), the point to be checked.
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
  real x(dim_num)

  ival = 0

  if ( dim_num == 1 ) then

    if ( 0.0 <= x(1) .and. x(1) <= 10.0 ) then
      ival = 1
    end if

  else if ( dim_num == 2 ) then

    if ( 0.0 <= x(1) .and. x(1) <= 10.0 .and. &
         0.0 <= x(2) .and. x(2) <= 10.0 ) then
      ival = 1
    end if

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

  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone

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
