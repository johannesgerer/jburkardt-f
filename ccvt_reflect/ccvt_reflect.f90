program main

!*****************************************************************************80
!
!! MAIN is the main program for CCVT_REFLECT.
!
!  Discussion:
!
!    CCVT_REFLECT uses reflection to "stretch" a CVT onto the boundary.
!
!    This code essentially carries out a standard CVT iteration.
!    However, one goal is to adjust the mesh so that some points
!    are moved onto the boundary in a smooth way.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) AV_POINTS(2,N), used to store the running
!    total and then average of the sample points closest to each
!    generator.
!
!    Local, integer COUNT(2,N), counts the number of sample points that
!    were nearest to each generator.
!
!    Local, real ( kind = 8 ) GENERATOR(2,N), the coordinates of the
!    CVT generators.
!
!    Local, integer IT_MAX, the number of CVT iterations to carry out.
!    (the user specifies this.)
!
!    Local, integer N, the number of generators (the user specifies this).
!
!    Local, integer SEED, a seed for the random number generator.
!
  implicit none

  real ( kind = 8 ), allocatable :: av_points(:,:)
  integer batch
  real ( kind = 8 ), allocatable :: count(:)
  logical, parameter :: DEBUG = .true.
  real ( kind = 8 ) :: energy
  real ( kind = 8 ) :: energy2
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: generator
  integer init
  character ( len = 80 ) init_string
  logical initialize
  logical inside
  integer ios
  integer it_fixed
  integer it_max
  integer it_num
  integer j
  integer n
  integer dim_num
  integer nearest(1)
  integer npp
  integer outside
  real ( kind = 8 ), allocatable, dimension ( : ) :: s
  logical s_eqi
  real ( kind = 8 ), allocatable, dimension ( : ) :: s2
  integer sample
  integer sample_num
  character ( len = 80 ) sample_string
  integer seed
  integer seed_base
  integer seed_init
  integer seed_init2
  integer seed_iter

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CCVT_REFLECT:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  CVT calculation constrained to a box.'
  write ( *, '(a)' ) '  Some sample points are "reflected" to try'
  write ( *, '(a)' ) '  to drive generators to the boundary.'
!
!  Get some input from the user.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DIM_NUM is the spatial dimension.'
  write ( *, '(a)' ) '  (Try ''2'' if you have no preference.)'
  write ( *, '(a)' ) '  (Any value less than 1 terminates execution.)'

  read ( *, * ) dim_num
  write ( *, '(a,i12)' ) '  User input DIM_NUM = ', dim_num

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT'
    write ( *, '(a,i8)' ) '  The input value of DIM_NUM = ', dim_num
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N is the number of points to generate.'
  write ( *, '(a)' ) '  (Try ''100'' if you have no preference.)'
  write ( *, '(a)' ) '  (Any value less than 1 terminates execution.)'

  read ( *, * ) n
  write ( *, '(a,i12)' ) '  User input N = ', n

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NPP is the number of sample points used to'
  write ( *, '(a)' ) '  check the boundary.'
  write ( *, '(a)' ) '  (Try ''1000'' if you have no preference.)'
  write ( *, '(a)' ) '  (Any value less than 1 terminates execution.)'

  read ( *, * ) npp
  write ( *, '(a,i12)' ) '  User input NPP = ', npp

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter SEED, a seed for the random number generator:'
  write ( *, '(a)' ) '  (Try ''123456789'' if you do not have a preference.)'
  write ( *, '(a)' ) '  (Any negative value terminates execution).'

  read ( *, *, iostat = ios ) seed

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for SEED.'
    stop
  end if

  write ( *, '(a,i12)' ) '  User input SEED = ', seed
  seed_init = seed

  if ( seed < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT'
    write ( *, '(a,i12)' ) '  The input value of SEED = ', seed
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  INIT is the method of initializing the data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  file     read data from a file;'
  write ( *, '(a)' ) '  GRID     by picking points from a grid;'
  write ( *, '(a)' ) '  HALTON   from a Halton sequence;'
  write ( *, '(a)' ) '  RANDOM   using FORTRAN RANDOM function;'
  write ( *, '(a)' ) '  UNIFORM  using a simple uniform RNG;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Try ''RANDOM'' if you do not have a preference.)'
  write ( *, '(a)' ) '  (A blank value terminates execution).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter INIT:'
  read ( *, '(a)', iostat = ios ) init_string

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for INIT.'
    stop
  end if

  write ( *, '(a)' ) '  User input INIT = "' // trim ( init_string ) // '".'

  if ( s_eqi ( init_string(1:6), 'RANDOM'  ) ) then
    init = -1
  else if ( s_eqi ( init_string(1:7), 'UNIFORM' ) ) then
    init = 0
  else if ( s_eqi ( init_string(1:6), 'HALTON'  ) ) then
    init = 1
  else if ( s_eqi ( init_string(1:4), 'GRID'    ) ) then
    init = 2
  else if ( 0 < len_trim ( init_string ) ) then
    init = 3
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT'
    write ( *, '(a)' ) '  The input value of INIT '
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  if ( len_trim ( init_string ) <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT'
    write ( *, '(a)' ) '  The input value of INIT '
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) '  User input INIT = "' // trim ( init_string ) // '".'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IT_MAX is the maximum number of iterations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  An iteration carries out the following steps:'
  write ( *, '(a)' ) '  * the Voronoi region associated with each'
  write ( *, '(a)' ) '    generator is estimated by sampling;'
  write ( *, '(a)' ) '  * the centroid of each Voronoi region is estimated.'
  write ( *, '(a)' ) '  * the generator is replaced by the centroid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If "enough" sampling points are used,'
  write ( *, '(a)' ) '  and "enough" iterations are taken, this process'
  write ( *, '(a)' ) '  will converge.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Try ''50'' if you have no preference.)'
  write ( *, '(a)' ) '  (A negative value terminates execution).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter IT_MAX:  '

  read ( *, * ) it_max
  write ( *, '(a,i12)' ) '  User input IT_MAX = ', it_max

  if ( it_max < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT'
    write ( *, '(a,i12)' ) '  The input value of IT_MAX = ', it_max
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IT_FIXED is the number of consecutive iterations'
  write ( *, '(a)' ) '  to take with a fixed set of sample points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Setting IT_FIXED to 1 means a new set of sample'
  write ( *, '(a)' ) '  points is generated on every iterative step;'
  write ( *, '(a)' ) '  Setting IT_FIXED equal to IT_MAX means a single set'
  write ( *, '(a)' ) '  of sample points is used for the entire iteration.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Any value between 1 and IT_MAX is reasonable.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  (Try ', it_max, &
    ' if you do not have a preference.)'
  write ( *, '(a)' ) '  (A 0 or negative value terminates execution).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter IT_FIXED:'
  read ( *, *, iostat = ios ) it_fixed

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for IT_FIXED.'
    stop
  end if

  write ( *, '(a,i8)' ) '  User input IT_FIXED = ', it_fixed

  if ( it_fixed <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT'
    write ( *, '(a,i12)' ) '  The input value of IT_FIXED = ', it_fixed
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SAMPLE is the method of sampling the region:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GRID     by picking points from a grid;'
  write ( *, '(a)' ) '  HALTON   from a Halton sequence;'
  write ( *, '(a)' ) '  RANDOM   using FORTRAN RANDOM function;'
  write ( *, '(a)' ) '  UNIFORM  using a simple uniform RNG;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Try ''RANDOM'' if you do not have a preference.)'
  write ( *, '(a)' ) '  (A blank value terminates execution).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter SAMPLE:'
  read ( *, '(a)', iostat = ios ) sample_string

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for SAMPLE.'
    stop
  end if

  write ( *, '(a)' ) '  User input SAMPLE = "' // trim ( sample_string ) &
    // '".'

  if ( len_trim ( sample_string ) <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT'
    write ( *, '(a)' ) '  The input value of SAMPLE '
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  if ( s_eqi ( sample_string(1:6), 'RANDOM'  ) ) then
    sample = -1
  else if ( s_eqi ( sample_string(1:7), 'UNIFORM' ) ) then
    sample = 0
  else if ( s_eqi ( sample_string(1:6), 'HALTON'  ) ) then
    sample = 1
  else if ( s_eqi ( sample_string(1:4), 'GRID'    ) ) then
    sample = 2
  else
    sample = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT'
    write ( *, '(a)' ) '  The input value of SAMPLE '
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SAMPLE_NUM is the number of sample points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Voronoi regions will be explored by generating'
  write ( *, '(a)' ) '  SAMPLE_NUM points.  For each sample point, the'
  write ( *, '(a)' ) '  nearest generator is found.  Using more points'
  write ( *, '(a)' ) '  gives a better estimate of these regions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SAMPLE_NUM should be much larger than N, the'
  write ( *, '(a)' ) '  number of generators. '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Try ''10000'' if you have no preference.) '
  write ( *, '(a)' ) '  (A zero or negative value terminates execution.)'
  write ( *, '(a)' ) ' '

  read ( *, * ) sample_num
  write ( *, '(a,i12)' ) '  User input SAMPLE_NUM = ', sample_num

  if ( sample_num <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT'
    write ( *, '(a,i12)' ) '  The input value of SAMPLE_NUM = ', sample_num
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BATCH is the number of sample points to create'
  write ( *, '(a)' ) '  at one time.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BATCH should be between 1 and SAMPLE_NUM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  It is FASTER to set BATCH to SAMPLE_NUM;'
  write ( *, '(a)' ) '  setting BATCH to 1 requires the least memory.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12,a)' ) '  (Try ', min ( sample_num, 1000 ), &
    ' if you do not have a preference.)'
  write ( *, '(a)' ) '  (A zero or negative value terminates execution.)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter BATCH:'

  read ( *, *, iostat = ios ) batch

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for SAMPLE_NUM.'
    stop
  end if

  write ( *, '(a,i12)' ) '  User input BATCH = ', batch

  if ( batch <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCVT_REFLECT'
    write ( *, '(a,i12)' ) '  The input value of BATCH = ', batch
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if
!
!  Allocate space.
!
  allocate ( av_points(2,n) )
  allocate ( count(n) )
  allocate ( generator(dim_num,n) )
  allocate ( s(dim_num) )
  allocate ( s2(dim_num) )
!
!  Print the header for output.
!
  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '                            Energy    Boundary'
    write ( *, '(a)' ) '  Iteration     Seed        of CVT    Samples'
    write ( *, '(a)' ) ' '
  end if
!
!  Initialize the generators.
!
  seed_iter = seed_init

  seed = seed_init
  initialize = .true.

  call cvt_sample ( dim_num, n, n, init, initialize, seed, generator )
!
!  If the initialization and sampling steps use the same random number
!  scheme, then the sampling scheme does not have to be initialized.
!
  if ( sample == init ) then
    initialize = .false.
  else
    initialize = .true.
  end if
!
!  Compute the initial energy.
!
  seed_init2 = seed_init

  call cvt_energy ( dim_num, n, batch, sample, initialize, sample_num, seed, &
    generator, energy )

  initialize = .false.

  it_num = 0
  outside = 0

  write ( *, '(2x,i8,2x,i12,2x,g14.6,2x,i12)' ) &
    it_num, seed_init2, energy, outside
!
!  Write out the initial points.
!
  call cvt_write ( dim_num, n, batch, seed_init, seed, init_string, it_max, &
    it_fixed, it_num, energy2, sample_string, sample_num, generator, &
    'initial.txt' )

  call points_eps ( 'initial.eps', dim_num, n, generator, 'Initial generators' )
!
!   Start the iteration.
!
  do while ( it_num <= it_max )

    if ( mod ( it_num, it_fixed ) == 0 ) then
      seed_base = seed
    else
      seed = seed_base
    end if

    it_num = it_num + 1
!
!  Sample the region.
!
    seed_init2 = seed
    av_points(1:dim_num,1:n) = 0.0D+00
    count(1:n) = 0

    outside = 0

    do j = 1, sample_num

      call cvt_sample ( dim_num, sample_num, 1, sample, initialize, seed, s )

      call find_closest ( dim_num, n, 1, s, generator, nearest )

      s2(1:dim_num) = 2.0D+00 * generator(1:dim_num,nearest(1)) &
        - s(1:dim_num)

      inside = 0.0D+00 <= s2(1) .and. s2(1) <= 1.0D+00 .and. &
               0.0D+00 <= s2(2) .and. s2(2) <= 1.0D+00

      if ( inside ) then

        av_points(1:dim_num,nearest(1)) = av_points(1:dim_num,nearest(1)) &
          + s(1:dim_num)

        count(nearest(1)) = count(nearest(1)) + 1

      else

        av_points(1:dim_num,nearest(1)) = av_points(1:dim_num,nearest(1)) &
          + generator(1:dim_num,nearest(1))

        count(nearest(1)) = count(nearest(1)) + 1
        outside = outside + 1

      end if

    end do
!
!  Replace the generators by the average of the sample points.
!
    do j = 1, n

      if ( count(j) /= 0 ) then

        generator(1:dim_num,j) = av_points(1:dim_num,j) &
          / real ( count(j), kind = 8 )

      end if

    end do
!
!  Compute the energy of the updated CVT.
!
    seed = seed_iter
    seed_init2 = seed

    call cvt_energy ( dim_num, n, batch, sample, initialize, sample_num, seed, &
      generator, energy )

    write ( *, '(2x,i8,2x,i12,2x,g14.6,2x,i12)' ) &
      it_num, seed_init2, energy, outside

  end do
!
!  Write out and plot the final  (unprojected) points.
!
  call cvt_write ( dim_num, n, batch, seed_init, seed, init_string, it_max, &
    it_fixed, it_num, energy2, sample_string, sample_num, generator, &
    'final.txt' )

  call points_eps ( 'final.eps', dim_num, n, generator, 'Final generators' )
!
!  Once the iteration is done, move points onto the boundary using
!  Lili's method.
!
    call mpb ( dim_num, n, generator, npp )

    seed = seed_iter
    seed_init2 = seed

    call cvt_energy ( dim_num, n, batch, sample, initialize, sample_num, seed, &
      generator, energy )

    it_num = it_max + 1

    write ( *, '(2x,i8,2x,i12,2x,g14.6,2x,i12)' ) &
      it_num, seed_init2, energy
!
!  Write out and plot the final projected points.
!
  call cvt_write ( dim_num, n, batch, seed_init, seed, init_string, it_max, &
    it_fixed, it_num, energy2, sample_string, sample_num, generator, &
    'projected.txt' )

  call points_eps ( 'projected.eps', dim_num, n, generator, &
    'Projected generators' )
!
!   Deallocate memory.
!
  deallocate ( av_points )
  deallocate ( count )
  deallocate ( generator )
  deallocate ( s )
  deallocate ( s2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CCVT_REFLECT:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
  integer itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine cvt_energy ( dim_num, n, batch, sample, initialize, sample_num, &
  seed, r, energy )

!*****************************************************************************80
!
!! CVT_ENERGY computes the CVT energy of a dataset.
!
!  Discussion:
!
!    For a given number of generators, a CVT is a minimizer (or at least
!    a local minimizer) of the CVT energy.  During a CVT iteration,
!    it should generally be the case that the CVT energy decreases from
!    step to step, and that perturbations or adjustments of an
!    approximate CVT will almost always have higher CVT energy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of generators.
!
!    Input, integer BATCH, the maximum number of sample points to generate
!    at one time.
!
!    Input, integer SAMPLE, specifies how the sampling is done.
!    -1, 'RANDOM', using FORTRAN RANDOM function;
!     0, 'UNIFORM', using a simple uniform RNG;
!     1, 'HALTON', from a Halton sequence;
!     2, 'GRID', points from a grid;
!
!    Input, logical INITIALIZE, is TRUE if the pseudorandom process should be
!    reinitialized.
!
!    Input, integer SAMPLE_NUM, the number of sample points to use.
!
!    Input/output, integer SEED, a seed for the random number generator.
!
!    Input, real ( kind = 8 ) R(DIM_NUM,N), the coordinates of the points.
!
!    Output, real ( kind = 8 ) ENERGY, the estimated CVT energy.
!
  implicit none

  integer batch
  integer dim_num
  integer n

  real ( kind = 8 ) energy
  integer get
  integer have
  logical initialize
  integer j
  integer nearest(batch)
  real ( kind = 8 ) r(dim_num,n)
  real ( kind = 8 ) s(dim_num,batch)
  integer sample
  integer sample_num
  integer seed

  have = 0
  energy = 0.0D+00

  do while ( have < sample_num )

    get = min ( sample_num - have, batch )

    call cvt_sample ( dim_num, sample_num, get, sample, initialize, seed, s )

    have = have + get

    call find_closest ( dim_num, n, get, s, r, nearest )

    do j = 1, get

      energy = energy &
        + sum ( ( s(1:dim_num,j) - r(1:dim_num,nearest(j)) )**2 )
    end do

  end do

  energy = energy / real ( sample_num, kind = 8 )

  return
end
subroutine cvt_sample ( dim_num, n, n_now, sample, initialize, seed, r )

!******************************************************************************
!
!! CVT_SAMPLE returns sample points.
!
!  Discussion:
!
!    N sample points are to be taken from the unit box of dimension DIM_NUM.
!
!    These sample points are usually created by a pseudorandom process
!    for which the points are essentially indexed by a quantity called
!    SEED.  To get N sample points, we generate values with indices
!    SEED through SEED+N-1.
!
!    It may not be practical to generate all the sample points in a
!    single call.  For that reason, the routine allows the user to
!    request a total of N points, but to require that only N_NOW be
!    generated now (on this call).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of sample points to be generated.
!
!    Input, integer N_NOW, the number of sample points to be generated
!    on this call.  N_NOW must be at least 1.
!
!    Input, integer SAMPLE, specifies how the sampling is done.
!    -1, 'RANDOM', using FORTRAN RANDOM function;
!     0, 'UNIFORM', using a simple uniform RNG;
!     1, 'HALTON', from a Halton sequence;
!     2, 'GRID', points from a grid;
!
!    Input, logical INITIALIZE, is TRUE if the pseudorandom process should be
!    reinitialized.
!
!    Input/output, integer SEED, the random number seed.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N_NOW), the sample points.
!
  implicit none

  integer dim_num
  integer n_now

  real ( kind = 8 ) exponent
  integer, allocatable, dimension ( : ) :: halton_base
  integer, allocatable, dimension ( : ) :: halton_leap
  integer, allocatable, dimension ( : ) :: halton_seed
  integer halton_step
  integer i
  logical initialize
  integer j
  integer n
  integer ngrid
  integer prime
  real ( kind = 8 ) :: r(dim_num,n_now)
  integer rank
  integer rank_max
  integer sample
  integer seed
  integer, allocatable, dimension ( : ) :: tuple

  if ( n_now < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  N_NOW < 1.'
    stop
  end if

  if ( sample == -1 ) then

    if ( initialize ) then
      call random_initialize ( seed )
    end if

    call random_number ( harvest = r(1:dim_num,1:n_now) )

    seed = seed + n_now * dim_num

  else if ( sample == 0 ) then

    call r8mat_uniform_01 ( dim_num, n_now, seed, r )

  else if ( sample == 1 ) then

    allocate ( halton_seed(1:dim_num) )
    allocate ( halton_leap(1:dim_num) )
    allocate ( halton_base(1:dim_num) )

    halton_step = seed
    halton_seed(1:dim_num) = 0
    halton_leap(1:dim_num) = 1

    do i = 1, dim_num
      halton_base(i) = prime ( i )
    end do

    call i4_to_halton_sequence ( dim_num, n_now, halton_step, halton_seed, &
      halton_leap, halton_base, r(1:dim_num,1:n_now) )

    deallocate ( halton_seed )
    deallocate ( halton_leap )
    deallocate ( halton_base )

    seed = seed + n_now

  else if ( sample == 2 ) then

    allocate ( tuple(1:dim_num) )

    exponent = real ( 1.0, kind = 8 ) / real ( dim_num, kind = 8 )
    ngrid = int ( ( real ( n, kind = 8 ) )**exponent )
    rank_max = ngrid**dim_num

    if ( rank_max < n ) then
      ngrid = ngrid + 1
      rank_max = ngrid**dim_num
    end if

    if ( initialize ) then
      rank = -1
      call tuple_next_fast ( ngrid, dim_num, rank, tuple )
    end if

    rank = mod ( seed, rank_max )

    do j = 1, n_now
      call tuple_next_fast ( ngrid, dim_num, rank, tuple )
      rank = rank + 1
      rank = mod ( rank, rank_max )
      r(1:dim_num,j) = real ( 2 * tuple(1:dim_num) - 1, kind = 8 ) &
        / real ( 2 * ngrid, kind = 8 )
    end do

    seed = seed + n_now

    deallocate ( tuple )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_SAMPLE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  The value of SAMPLE = ', sample, ' is illegal.'
    stop

  end if

  return
end
subroutine cvt_write ( dim_num, n, batch, seed_init, seed, init_string, it_max, &
  it_fixed, it_num, energy, sample_string, sample_num, r, &
  file_out_name )

!*****************************************************************************80
!
!! CVT_WRITE writes a CVT dataset to a file.
!
!  Discussion:
!
!    The initial lines of the file are comments, which begin with a
!    '#' character.
!
!    Thereafter, each line of the file contains the DIM_NUM-dimensional
!    components of the next entry of the dataset.
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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, integer BATCH, sets the maximum number of sample points
!    generated at one time.  It is inefficient to generate the sample
!    points 1 at a time, but memory intensive to generate them all
!    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
!
!    Input, integer SEED_INIT, the initial random number seed.
!
!    Input, integer SEED, the current random number seed.
!
!    Input, character ( len = * ) INIT_STRING, specifies how the initial
!    generators are chosen:
!    filename, by reading data from a file;
!    'GRID', picking points from a grid;
!    'HALTON', from a Halton sequence;
!    'RANDOM', using FORTRAN RANDOM function;
!    'UNIFORM', using a simple uniform RNG;
!
!    Input, integer IT_MAX, the maximum number of iterations allowed.
!
!    Input, integer IT_FIXED, the number of iterations to take with a
!    fixed set of sample points.
!
!    Input, integer IT_NUM, the actual number of iterations taken.
!
!    Input, real ( kind = 8 ) ENERGY, the discrete "energy", divided
!    by the number of sample points.
!
!    Input, character ( len = * ) SAMPLE_STRING, specifies how the 
!    region is sampled:
!    'GRID', picking points from a grid;
!    'HALTON', from a Halton sequence;
!    'RANDOM', using FORTRAN RANDOM function;
!    'UNIFORM', using a simple uniform RNG;
!
!    Input, integer SAMPLE_NUM, the number of sampling points used on
!    each iteration.
!
!    Input, real ( kind = 8 ) R(DIM_NUM,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of
!    the output file.
!
  implicit none

  integer n
  integer dim_num

  integer batch
  real ( kind = 8 ) energy
  character ( len = * ) file_out_name
  integer file_out_unit
  character ( len = * ) init_string
  integer ios
  integer it_fixed
  integer it_max
  integer it_num
  integer j
  real ( kind = 8 ) r(dim_num,n)
  character ( len = * ) sample_string
  integer sample_num
  integer seed
  integer seed_init
  character ( len = 40 ) string

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, &
    status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  call timestring ( string )

  write ( file_out_unit, '(a)'       ) '#  ' // trim ( file_out_name )
  write ( file_out_unit, '(a)'       ) '#  created by routine CVT_WRITE.F90'
  write ( file_out_unit, '(a)'       ) '#  at ' // trim ( string )
  write ( file_out_unit, '(a)'       ) '#'
  write ( file_out_unit, '(a,i12)'   ) '#  Spatial dimension DIM_NUM =   ', dim_num
  write ( file_out_unit, '(a,i12)'   ) '#  Number of points N =       ', n
  write ( file_out_unit, '(a,i12)'   ) '#  Initial SEED_INIT =        ', &
    seed_init
  write ( file_out_unit, '(a,i12)'   ) '#  Current SEED =             ', seed
  write ( file_out_unit, '(a)'       ) '#  INIT =                    "' &
    // trim ( init_string ) // '".'
  write ( file_out_unit, '(a,i12)'   ) '#  Max iterations IT_MAX =    ', it_max
  write ( file_out_unit, '(a,i12)'   ) '#  IT_FIXED (fixed samples) = ', &
    it_fixed
  write ( file_out_unit, '(a,i12)'   ) '#  Iterations IT_NUM =        ', it_num
  write ( file_out_unit, '(a,g14.6)' ) '#  CVT ENERGY =               ', energy
  write ( file_out_unit, '(a)'       ) '#  SAMPLE =                  "' &
    // trim ( sample_string ) // '".'
  write ( file_out_unit, '(a,i12)'   ) '#  Samples SAMPLE_NUM    =    ', &
    sample_num
  write ( file_out_unit, '(a,i12)'   ) '#  Sampling BATCH size   =    ', batch
  write ( file_out_unit, '(a,g14.6)' ) '#  EPSILON (unit roundoff) =  ', &
    epsilon ( r(1,1) )
  write ( file_out_unit, '(a)'       ) '#'

  write ( string, '(a,i3,a)' ) '(', dim_num, 'f10.6)'
  do j = 1, n
    write ( file_out_unit, string ) r(1:dim_num,j)
  end do

  close ( unit = file_out_unit )

  return
end
subroutine find_closest ( dim_num, n, sample_num, s, r, nearest )

!******************************************************************************
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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of cell generators.
!
!    Input, integer SAMPLE_NUM, the number of sample points.
!
!    Input, real ( kind = 8 ) S(DIM_NUM,SAMPLE_NUM), the points to be checked.
!
!    Input, real ( kind = 8 ) R(DIM_NUM,N), the cell generators.
!
!    Output, integer NEAREST(SAMPLE_NUM), the index of the nearest
!    cell generators.
!
  implicit none

  integer n
  integer dim_num
  integer sample_num

  real ( kind = 8 ) dist_sq_min
  real ( kind = 8 ) dist_sq
  integer jr
  integer js
  integer nearest(sample_num)
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
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
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
function halham_leap_check ( dim_num, leap )

!*****************************************************************************80
!
!! HALHAM_LEAP_CHECK checks LEAP for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer LEAP(DIM_NUM), the leap vector.
!
!    Output, logical, HALHAM_LEAP_CHECK, true if LEAP is legal.
!
  implicit none

  integer dim_num

  logical halham_leap_check
  integer leap(dim_num)

  if ( any ( leap(1:dim_num) < 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_LEAP_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some entry of LEAP < 1!'
    write ( *, '(a)' ) ' '
    call i4vec_transpose_print ( dim_num, leap, 'LEAP:  ' )
    halham_leap_check = .false.
  else
    halham_leap_check = .true.
  end if

  return
end
function halham_n_check ( n )

!*****************************************************************************80
!
!! HALHAM_N_CHECK checks N for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the spatial dimension.
!
!    Output, logical HALHAM_N_CHECK, true if N is legal.
!
  implicit none

  logical halham_n_check
  integer n

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_N_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i12)' ) '  N = ', n
    halham_n_check = .false.
  else
    halham_n_check = .true.
  end if

  return
end
function halham_dim_num_check ( dim_num )

!*****************************************************************************80
!
!! HALHAM_DIM_NUM_CHECK checks DIM_NUM for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Output, logical HALHAM_DIM_NUM_CHECK, true if DIM_NUM is legal.
!
  implicit none

  logical halham_dim_num_check
  integer dim_num

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_DIM_NUM_CHECK - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
    halham_dim_num_check = .false.
  else
    halham_dim_num_check = .true.
  end if

  return
end
function halham_seed_check ( dim_num, seed )

!*****************************************************************************80
!
!! HALHAM_SEED_CHECK checks SEED for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer SEED(DIM_NUM), the seed vector.
!
!    Output, logical, HALHAM_SEED_CHECK, true if SEED is legal.
!
  implicit none

  integer dim_num

  logical halham_seed_check
  integer seed(dim_num)

  if ( any ( seed(1:dim_num) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_SEED_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some entry of SEED < 0!'
    write ( *, '(a)' ) ' '
    call i4vec_transpose_print ( dim_num, seed, 'SEED:  ' )
    halham_seed_check = .false.
  else
    halham_seed_check = .true.
  end if

  return
end
function halham_step_check ( step )

!*****************************************************************************80
!
!! HALHAM_STEP_CHECK checks STEP for a Halton or Hammersley sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer STEP, the index of the subsequence element.
!
!    Output, logical HALHAM_STEP_CHECK, true if STEP is legal.
!
  implicit none

  logical halham_step_check
  integer step

  if ( step < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_STEP_CHECK - Fatal error!'
    write ( *, '(a)' ) '  STEP < 0.'
    write ( *, '(a,i12)' ) '  STEP = ', step
    halham_step_check = .false.
  else
    halham_step_check = .true.
  end if

  return
end
function halton_base_check ( dim_num, base )

!*****************************************************************************80
!
!! HALTON_BASE_CHECK checks BASE for a Halton sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer BASE(DIM_NUM), the bases.
!
!    Output, logical, HALTON_BASE_CHECK, true if BASE is legal.
!
  implicit none

  integer dim_num

  integer base(dim_num)
  logical halton_base_check

  if ( any ( base(1:dim_num) <= 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALTON_BASE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some entry of BASE is <= 1!'
    write ( *, '(a)' ) ' '
    call i4vec_transpose_print ( dim_num, base, 'BASE:  ' )
    halton_base_check = .false.
  else
    halton_base_check = .true.
  end if

  return
end
subroutine i4_to_halton_sequence ( dim_num, n, step, seed, leap, base, r )

!*****************************************************************************80
!
!! I4_TO_HALTON_SEQUENCE computes N elements of a leaped Halton subsequence.
!
!  Discussion:
!
!    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
!    sequences, each generated by a particular base.
!
!    This routine selects elements of a "leaped" subsequence of the
!    Halton sequence.  The subsequence elements are indexed by a
!    quantity called STEP, which starts at 0.  The STEP-th subsequence
!    element is simply element
!
!      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
!
!    of the original Halton sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J H Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, 1960, pages 84-90.
!
!    J H Halton and G B Smith,
!    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
!    Communications of the ACM,
!    Volume 7, 1964, pages 701-702.
!
!    Ladislav Kocis and William Whiten,
!    Computational Investigations of Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 2, 1997, pages 266-294.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!    1 <= DIM_NUM is required.
!
!    Input, integer N, the number of elements of the sequence.
!
!    Input, integer STEP, the index of the subsequence element.
!    0 <= STEP is required.
!
!    Input, integer SEED(DIM_NUM), the Halton sequence index corresponding
!    to STEP = 0.
!
!    Input, integer LEAP(DIM_NUM), the succesive jumps in the Halton sequence.
!
!    Input, integer BASE(DIM_NUM), the Halton bases.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the next N elements of the
!    leaped Halton subsequence, beginning with element STEP.
!
  implicit none

  integer n
  integer dim_num

  integer base(dim_num)
  real ( kind = 8 ) base_inv
  integer digit(n)
  logical halham_leap_check
  logical halham_n_check
  logical halham_dim_num_check
  logical halham_seed_check
  logical halham_step_check
  logical halton_base_check
  integer i
  integer j
  integer leap(dim_num)
  real ( kind = 8 ) r(dim_num,n)
  integer seed(dim_num)
  integer seed2(n)
  integer step
!
!  Check the input.
!
  if ( .not. halham_dim_num_check ( dim_num ) ) then
    stop
  end if

  if ( .not. halham_n_check ( n ) ) then
    stop
  end if

  if ( .not. halham_step_check ( step ) ) then
    stop
  end if

  if ( .not. halham_seed_check ( dim_num, seed ) ) then
    stop
  end if

  if ( .not. halham_leap_check ( dim_num, leap ) ) then
    stop
  end if

  if ( .not. halton_base_check ( dim_num, base ) ) then
    stop
  end if
!
!  Calculate the data.
!
  r(1:dim_num,1:n) = 0.0D+00

  do i = 1, dim_num

    do j = 1, n
      seed2(j) = seed(i) + ( step + j - 1 ) * leap(i)
    end do

    base_inv = real ( 1.0D+00, kind = 8 ) / real ( base(i), kind = 8 )

    do while ( any ( seed2(1:n) /= 0 ) )
      digit(1:n) = mod ( seed2(1:n), base(i) )
      r(i,1:n) = r(i,1:n) + real ( digit(1:n), kind = 8 ) * base_inv
      base_inv = base_inv / real ( base(i), kind = 8 )
      seed2(1:n) = seed2(1:n) / base(i)
    end do

  end do

  return
end
subroutine i4vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
!
!  Example:
!
!    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
!    TITLE = 'My vector:  '
!
!    My vector:      1    2    3    4    5
!                    6    7    8    9   10
!                   11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  integer a(n)
  integer ihi
  integer ilo
  character ( len = 11 ) string
  character ( len = * ) title
  integer title_len

  if ( 0 < len_trim ( title ) ) then

    title_len = len ( title )

    write ( string, '(a,i3,a)' ) '(', title_len, 'x,5i12)'

    do ilo = 1, n, 5
      ihi = min ( ilo + 5 - 1, n )
      if ( ilo == 1 ) then
        write ( *, '(a, 5i12)' ) title, a(ilo:ihi)
      else
        write ( *, string      )        a(ilo:ihi)
      end if
    end do

  else

    do ilo = 1, n, 5
      ihi = min ( ilo + 5 - 1, n )
      write ( *, '(5i12)' ) a(ilo:ihi)
    end do

  end if

  return
end
subroutine mpb ( dim_num, n, generator, npp )

!*****************************************************************************80
!
!! MPB projects generators onto the boundary of the region.
!
!  Discussion:
!
!    The number NPP sets the number of subintervals into which we subdivide
!    the boundary.  It does NOT specify how many points will be pulled onto
!    the boundary.  The reason for this is that, after the first boundary
!    subinterval has had a generator pulled into it, on every subsequent
!    subinterval the nearest generator is likely to be the one in the
!    previous subinterval!  Unless an interior generator is closer than
!    some small distance, this process will simply drag some unfortunate
!    generator onto the boundary, and then around from interval to interval
!    for a considerable time.
!
!    The algorithm could be changed, if desired, so that points snapped
!    to the boundary are guaranteed not to move, at least not twice in 
!    one application of this routine!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2004
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of generators.
!
!    Input/output, real ( kind = 8 ) GENERATOR(DIM_NUM,N), the coordinates of
!    the generators.  On output, some generators will have been moved.
!
!    Input, integer NPP, the number of subintervals into which the
!    perimeter is divided.
!
  implicit none

  integer n
  integer dim_num

  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) generator(dim_num,n)
  real ( kind = 8 ) hh
  integer i
  integer nearest(1)
  integer npp
  real ( kind = 8 ) u
  real ( kind = 8 ) sample(dim_num)

  dx = 1.0D+00
  dy = 1.0D+00
!
!  HH is the length of an individual segment of the perimeter of the region.
!
!  U is set in such a way that on step I, it measures the distance from
!  the lower left corner of the box to the midpoint of the I-th subinterval 
!  on the perimeter of the box.
!  
  hh = 2.0D+00 * ( dx + dy ) / real ( npp, kind = 8 )

  u = -0.5D+00 * hh

  do i = 1, npp

    u = u + hh
!
!  The portion of the bottom perimeter from (0,0) to (1,0).
!
    if ( u < dx ) then
      sample(1:2) = (/ u, 0.0D+00 /)
      call find_closest ( dim_num, n, 1, sample, generator, nearest )
      generator(2,nearest(1)) = 0.0D+00
!
!  The portion of the right perimeter from (1,0) to (1,1).
!
    else if ( dx < u .and. u < dx + dy ) then
      sample(1:2) = (/ 1.0D+00, ( u - dx ) /)
      call find_closest ( dim_num, n, 1, sample, generator, nearest )
      generator(1,nearest(1)) = 1.0D+00
!
!  The portion of the top perimeter from (1,1) to (0,1).
!
    else if ( dx + dy < u .and. u < dx + dy + dx ) then
      sample(1:2) = (/ 1.0D+00 - ( u - dx - dy ), 1.0D+00 /)
      call find_closest ( dim_num, n, 1, sample, generator, nearest )
      generator(2,nearest(1)) = 1.0D+00
!
!  The portion of the left perimeter from (0,1) to (0,0).
!
    else if ( dx + dy + dx < u ) then
      sample(1:2) = (/ 0.0D+00, 1.0D+00 - ( u - dx - dy - dx ) /)
      call find_closest ( dim_num, n, 1, sample, generator, nearest )
      generator(1,nearest(1)) = 0.0D+00
    end if  
 
  end do

  return           
end
subroutine points_eps ( file_name, dim_num, node_num, node_xy, title )

!*****************************************************************************80
!
!! POINTS_EPS creates an EPS file image of a set of points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to create.
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates
!    of the nodes.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
  implicit none

  integer node_num
  integer dim_num

  integer, parameter :: circle_size = 3
  logical, parameter :: DEBUG = .false.
  real ( kind = 8 ) dif
  integer eps_unit
  integer eps_x
  integer eps_y
  character ( len = * ) file_name
  integer ios
  integer node
  real ( kind = 8 ) node_xy(dim_num,node_num)
  real ( kind = 8 ) node_x_max
  real ( kind = 8 ) node_x_min
  real ( kind = 8 ) node_y_max
  real ( kind = 8 ) node_y_min
  real ( kind = 8 ) scale
  character ( len = 40 ) string
  character ( len = * ) title
!
!  Determine the range of the points.
!
  node_x_min =  huge ( node_x_min )
  node_x_max = -huge ( node_x_max )
  node_y_min =  huge ( node_y_min )
  node_y_max = -huge ( node_y_max )

  do node = 1, node_num
    node_x_min = min ( node_x_min, node_xy(1,node) )
    node_x_max = max ( node_x_max, node_xy(1,node) )
    node_y_min = min ( node_y_min, node_xy(2,node) )
    node_y_max = max ( node_y_max, node_xy(2,node) )
  end do

  if ( node_y_max - node_y_min < node_x_max - node_x_min ) then
    scale = node_x_max - node_x_min
    dif = ( node_x_max - node_x_min ) - ( node_y_max - node_y_min )
    node_y_max = node_y_max + 0.5D+00 * dif
    node_y_min = node_y_min - 0.5D+00 * dif
  else
    scale = node_y_max - node_y_min
    dif = ( node_y_max - node_y_min ) - ( node_x_max - node_x_min )
    node_x_max = node_x_max + 0.5D+00 * dif
    node_x_min = node_x_min - 0.5D+00 * dif
  end if

  call get_unit ( eps_unit )

  open ( unit = eps_unit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINTS_EPS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output EPS file.'
    stop
  end if

  call timestring ( string )

  write ( eps_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( eps_unit, '(a)' ) &
    '%%Creator: points_eps(ccvt_box.f90)'
  write ( eps_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( eps_unit, '(a)' ) '%%CreationDate: ' // trim ( string )
  write ( eps_unit, '(a)' ) '%%Pages: 1'
  write ( eps_unit, '(a)' ) '%%BoundingBox:    36    36   576   756'
  write ( eps_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( eps_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( eps_unit, '(a)' ) '%%EndComments'
  write ( eps_unit, '(a)' ) '%%BeginProlog'
  write ( eps_unit, '(a)' ) '/inch {72 mul} def'
  write ( eps_unit, '(a)' ) '%%EndProlog'
  write ( eps_unit, '(a)' ) '%%Page:      1     1'
  write ( eps_unit, '(a)' ) 'save'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Set RGB line color.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.9000 0.9000 0.9000 setrgbcolor'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Draw a gray border around the page.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'newpath'
  write ( eps_unit, '(a)' ) '    36   126 moveto'
  write ( eps_unit, '(a)' ) '   576   126 lineto'
  write ( eps_unit, '(a)' ) '   576   666 lineto'
  write ( eps_unit, '(a)' ) '    36   666 lineto'
  write ( eps_unit, '(a)' ) '    36   126 lineto'
  write ( eps_unit, '(a)' ) 'stroke'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Set RGB line color.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Label the plot:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'
  write ( eps_unit, '(a)' ) '/Times-Roman findfont 0.50 inch scalefont setfont'
  write ( eps_unit, '(a)' ) '    36   666 moveto'
  write ( eps_unit, '(a)' ) '(' // trim ( title ) // ') show'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Define a clipping polygon'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '    36   126 moveto'
  write ( eps_unit, '(a)' ) '   576   126 lineto'
  write ( eps_unit, '(a)' ) '   576   666 lineto'
  write ( eps_unit, '(a)' ) '    36   666 lineto'
  write ( eps_unit, '(a)' ) '    36   126 lineto'
  write ( eps_unit, '(a)' ) 'clip newpath'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw the boundary in red:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.9000 0.0000 0.0000 setrgbcolor'
  write ( eps_unit, '(a)' ) 'newpath'
  write ( eps_unit, '(a)' ) '    61   151 moveto'
  write ( eps_unit, '(a)' ) '   551   151 lineto'
  write ( eps_unit, '(a)' ) '   551   641 lineto'
  write ( eps_unit, '(a)' ) '    61   641 lineto'
  write ( eps_unit, '(a)' ) '    61   151 lineto'
  write ( eps_unit, '(a)' ) 'stroke'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw filled dots at each node:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.9000 setrgbcolor'

  do node = 1, node_num

    eps_x = int ( &
      ( ( node_x_max - node_xy(1,node)              ) *  61.0D+00   &
      + (            + node_xy(1,node) - node_x_min ) * 551.0D+00 ) &
      / scale )

    eps_y = int ( &
      ( ( node_y_max - node_xy(2,node)              ) * 151.0D+00   &
      + (              node_xy(2,node) - node_y_min ) * 641.0D+00 ) &
      / scale )

    write ( eps_unit, '(a,i4,2x,i4,2x,i4,a)' ) &
      'newpath  ', eps_x, eps_y, circle_size, ' 0 360 arc closepath fill'

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'restore showpage'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% End of page'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%%Trailer'
  write ( eps_unit, '(a)' ) '%%EOF'

  close ( unit = eps_unit )

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINTS_EPS:'
    write ( *, '(a)' ) '  An encapsulated PostScript file was created'
    write ( *, '(a)' ) '  containing an image of the points.'
    write ( *, '(a)' ) '  The file is named "' // trim ( file_name ) // '".'
  end if

  return
end
function prime ( n )

!*****************************************************************************80
!
!! PRIME returns any of the first PRIME_MAX prime numbers.
!
!  Discussion:
!
!    PRIME_MAX is 1600, and the largest prime stored is 13499.
!
!    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
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
!  Reference:
!
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964, pages 870-873.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 95-98.
!
!  Parameters:
!
!    Input, integer N, the index of the desired prime number.
!    In general, is should be true that 0 <= N <= PRIME_MAX.
!    N = -1 returns PRIME_MAX, the index of the largest prime available.
!    N = 0 is legal, returning PRIME = 1.
!
!    Output, integer PRIME, the N-th prime.  If N is out of range,
!    PRIME is returned as -1.
!
  implicit none

  integer, parameter :: prime_max = 1600

  integer, save :: icall = 0
  integer n
  integer, save, dimension ( prime_max ) :: npvec
  integer prime

  if ( icall == 0 ) then

    icall = 1

    npvec(1:100) = (/ &
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /)

    npvec(101:200) = (/ &
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /)

    npvec(201:300) = (/ &
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /)

    npvec(301:400) = (/ &
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /)

    npvec(401:500) = (/ &
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /)

    npvec(501:600) = (/ &
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /)

    npvec(601:700) = (/ &
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /)

    npvec(701:800) = (/ &
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /)

    npvec(801:900) = (/ &
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /)

    npvec(901:1000) = (/ &
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)

    npvec(1001:1100) = (/ &
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /)

    npvec(1101:1200) = (/ &
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /)

    npvec(1201:1300) = (/ &
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /)

    npvec(1301:1400) = (/ &
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949, &
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /)

    npvec(1401:1500) = (/ &
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, &
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, &
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, &
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, &
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, &
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, &
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, &
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, &
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, &
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /)

   npvec(1501:1600) = (/ &
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, &
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, &
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, &
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, &
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, &
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, &
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, &
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, &
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, &
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 /)

  end if

  if ( n == -1 ) then
    prime = prime_max
  else if ( n == 0 ) then
    prime = 1
  else if ( n <= prime_max ) then
    prime = npvec(n)
  else
    prime = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRIME - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i8)' ) '  N should be between 1 and PRIME_MAX =', prime_max
    stop
  end if

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
subroutine random_initialize ( seed )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
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
!    27 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer SEED, the value used to initialize the random 
!    number generator.
!
  implicit none

  logical, parameter :: DEBUG = .false.
  integer seed
  integer, allocatable :: seed_vector(:)
  integer seed_size
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

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, user SEED = ', seed
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

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer i
  integer len1
  integer len2
  integer lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

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
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
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
  character ( len = * ) string
  character ( len = 10 ) time
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

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tuple_next_fast ( m, n, rank, x )

!*****************************************************************************80
!
!! TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
!
!  Discussion:
!
!    The elements are N vectors.  Each entry is constrained to lie
!    between 1 and M.  The elements are produced one at a time.
!    The first element is
!      (1,1,...,1)
!    and the last element is
!      (M,M,...,M)
!    Intermediate elements are produced in lexicographic order.
!
!    This code was written as a possibly faster version of TUPLE_NEXT.
!
!  Example:
!
!    N = 2,
!    M = 3
!
!    INPUT        OUTPUT
!    -------      -------
!    Rank          X
!    ----          ----
!   -1            -1 -1
!
!    0             1  1
!    1             1  2
!    2             1  3
!    3             2  1
!    4             2  2
!    5             2  3
!    6             3  1
!    7             3  2
!    8             3  3
!    9             1  1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the maximum entry in any component.
!    M must be greater than 0.
!
!    Input, integer N, the number of components.
!    N must be greater than 0.
!
!    Input, integer RANK, indicates the rank of the tuple.
!    Typically, 0 <= RANK < N**M.  Values of RANK greater than
!    N**M are legal and meaningful; they are equivalent to the
!    corresponding value mod (N**M).  If RANK < 0, this indicates
!    that this is the first call for the given values of (M,N).
!    Initialization is done, and X is set to a dummy value.
!
!    Output, integer X(N), the next tuple, or a dummy value if
!    initialization has just been done.
!
  implicit none

  integer n

  integer, save, allocatable, dimension ( : ) :: base
  integer i
  integer m
  integer rank
  integer x(n)

  if ( rank < 0 ) then

    if ( m <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
      write ( *, '(a)' ) '  The value M <= 0 is not allowed.'
      write ( *, '(a,i8)' ) '  M = ', m
      stop
    end if

    if ( n <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
      write ( *, '(a)' ) '  The value N <= 0 is not allowed.'
      write ( *, '(a,i8)' ) '  N = ', n
      stop
    end if

    if ( allocated ( base ) ) then
      deallocate ( base )
    end if
    allocate ( base(1:n) )

    base(n) = 1
    do i = n-1, 1, -1
      base(i) = base(i+1) * m
    end do

    x(1:n) = -1

  else

    x(1:n) = mod ( rank / base(1:n), m ) + 1

  end if

  return
end
