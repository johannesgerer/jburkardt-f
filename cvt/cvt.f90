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
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
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
subroutine cvt ( dim_num, n, batch, init, sample, sample_num, it_max, &
  it_fixed, seed, r, it_num, it_diff, energy )

!*****************************************************************************80
!
!! CVT computes a Centroidal Voronoi Tessellation.
!
!  Discussion:
!
!    This routine carries out the CVT iteration.
!
!    It initializes the CVT generators, unless the user indicates that
!    this has already been done.
!
!    It sets a flag INITIALIZE that indicates whether the random number
!    generator needs to be initialized.
!
!    It decides whether to use a new seed on each iteration, or to
!    reuse a seed value.
!
!    It then repeatedly calls a routine that takes another step 
!    of the iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Qiang Du, Vance Faber, Max Gunzburger,
!    Centroidal Voronoi Tessellations: Applications and Algorithms,
!    SIAM Review, Volume 41, 1999, pages 637-676.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of Voronoi cells.
!
!    Input, integer ( kind = 4 ) BATCH, sets the maximum number of sample points
!    generated at one time.  It is inefficient to generate the sample
!    points 1 at a time, but memory intensive to generate them all
!    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
!    BATCH must be at least 1.
!
!    Input, integer ( kind = 4 ) INIT, specifies how the points are to be 
!    initialized.
!    -1, 'RANDOM', using FORTRAN RANDOM function;
!     0, 'UNIFORM', using a simple uniform RNG;
!     1, 'HALTON', from a Halton sequence;
!     2, 'GRID', points from a grid;
!     3, 'USER', call "user" routine;
!     4, points are already initialized on input.
!
!    Input, integer ( kind = 4 ) SAMPLE, specifies how the sampling is done.
!    -1, 'RANDOM', using FORTRAN RANDOM function;
!     0, 'UNIFORM', using a simple uniform RNG;
!     1, 'HALTON', from a Halton sequence;
!     2, 'GRID', points from a grid;
!     3, 'USER', call "user" routine.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of sample points.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Input, integer ( kind = 4 ) IT_FIXED, the maximum number of iterations to 
!    take with a fixed set of sample points.
!
!    Input/output, integer ( kind = 4 ) SEED, the current random number seed.
!
!    Input/output, real ( kind = 8 ) R(DIM_NUM,N), the approximate CVT points.
!    If INIT = 4 on input, then it is assumed that these values have been
!    initialized.  On output, the CVT iteration has been applied to improve
!    the value of the points.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.  
!    Generally, this will be equal to IT_MAX, unless the iteration tolerance was
!    satisfied early.
!
!    Output, real ( kind = 8 ) IT_DIFF, the L2 norm of the difference
!    between the iterates.
!
!    Output, real ( kind = 8 ) ENERGY, the discrete "energy", divided
!    by the number of sample points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) batch
  logical, parameter :: DEBUG = .true.
  real    ( kind = 8 ) energy
  integer ( kind = 4 ) init
  logical              initialize
  real    ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real    ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_base
  integer ( kind = 4 ) seed_init

  if ( batch < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT - Fatal error!'
    write ( *, '(a)' ) '  The input value BATCH < 1.'
    stop
  end if

  if ( seed <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT - Fatal error!'
    write ( *, '(a)' ) '  The input value SEED <= 0.'
    stop
  end if

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Step    SEED        L2-Change         Energy'
    write ( *, '(a)' ) ' '
  end if

  it_num = 0
  it_diff = 0.0D+00
  energy = 0.0D+00
  seed_init = seed
!
!  Initialize the data unless the user has already done that.
!
  if ( init /= 4 ) then

    initialize = .true.

    call cvt_sample ( dim_num, n, n, init, initialize, seed, r )

  end if
!
!  NOTE THIS IS NEW:...
!
  call cvt_energy ( dim_num, n, batch, sample, initialize, sample_num, &
    seed, r, energy )

  if ( DEBUG ) then
    write ( *, '(2x,i4,2x,i12,2x,14x,2x,g14.6)' ) &
      it_num, seed_init, energy
  end if
!
!  If the initialization and sampling steps use the same random number
!  scheme, then the sampling scheme does not have to be initialized.
!
  if ( init == sample ) then
    initialize = .false.
  else
    initialize = .true.
  end if
!
!  Carry out the iteration.
!
  do while ( it_num < it_max )
!
!  If it's time to update the seed, save its current value
!  as the starting value for all iterations in this cycle.
!  If it's not time to update the seed, restore it to its initial
!  value for this cycle.
!
    if ( mod ( it_num, it_fixed ) == 0 ) then
      seed_base = seed
    else
      seed = seed_base
    end if

    it_num = it_num + 1

    seed_init = seed

    call cvt_iterate ( dim_num, n, batch, sample, initialize, sample_num, &
      seed, r, it_diff, energy )

    initialize = .false.

    if ( DEBUG ) then
      write ( *, '(2x,i4,2x,i12,2x,g14.6,2x,g14.6)' ) &
        it_num, seed_init, it_diff, energy
    end if

  end do

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
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of generators.
!
!    Input, integer ( kind = 4 ) BATCH, the maximum number of sample points 
!    to generate at one time.
!
!    Input, integer ( kind = 4 ) SAMPLE, specifies how the sampling is done.
!    -1, 'RANDOM', using FORTRAN RANDOM function;
!     0, 'UNIFORM', using a simple uniform RNG;
!     1, 'HALTON', from a Halton sequence;
!     2, 'GRID', points from a grid;
!     3, 'USER', call "user" routine.
!
!    Input, logical INITIALIZE, is TRUE if the pseudorandom process should be
!    reinitialized.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of sample points to use.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Input, real ( kind = 8 ) R(DIM_NUM,N), the coordinates of the points.
!
!    Output, real ( kind = 8 ) ENERGY, the estimated CVT energy.
!
  implicit none

  integer ( kind = 4 ) batch
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real    ( kind = 8 ) energy
  integer ( kind = 4 ) get
  integer ( kind = 4 ) have
  logical              initialize
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nearest(batch)
  real    ( kind = 8 ) r(dim_num,n)
  real    ( kind = 8 ) s(dim_num,batch)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed

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
subroutine cvt_iterate ( dim_num, n, batch, sample, initialize, sample_num, &
  seed, r, it_diff, energy )

!*****************************************************************************80
!
!! CVT_ITERATE takes one step of the CVT iteration.
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
!    The centroidal Voronoi tessellation minimizes the "energy",
!    defined to be the integral, over the region, of the square of
!    the distance between each point in the region and its nearest generator.
!    The sampling technique supplies a discrete estimate of this
!    energy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Qiang Du, Vance Faber, Max Gunzburger,
!    Centroidal Voronoi Tessellations: Applications and Algorithms,
!    SIAM Review, Volume 41, 1999, pages 637-676.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points to generate.
!
!    Input, integer ( kind = 4 ) BATCH, sets the maximum number of sample points
!    generated at one time.  It is inefficient to generate the sample
!    points 1 at a time, but memory intensive to generate them all
!    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
!    BATCH must be at least 1.
!
!    Input, integer ( kind = 4 ) SAMPLE, specifies how the sampling is done.
!    -1, 'RANDOM', using FORTRAN RANDOM function;
!     0, 'UNIFORM', using a simple uniform RNG;
!     1, 'HALTON', from a Halton sequence;
!     2, 'GRID', points from a grid;
!     3, 'USER', call "user" routine.
!
!    Input, logical INITIALIZE, is TRUE if the random number generator
!    should be initialized, because this is the first call to it.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of sample points.
!
!    Input/output, integer ( kind = 4 ) SEED, the random number seed.
!
!    Input/output, real ( kind = 8 ) R(DIM_NUM,N), the Voronoi
!    cell generators.  On output, these have been modified
!
!    Output, real ( kind = 8 ) IT_DIFF, the L2 norm of the difference
!    between the iterates.
!
!    Output, real ( kind = 8 ) ENERGY, the discrete "energy", divided
!    by the number of sample points.
!
  implicit none

  integer ( kind = 4 ) batch
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n
  integer ( kind = 4 ) sample_num

  integer ( kind = 4 ) count(n)
  real    ( kind = 8 ) energy
  integer ( kind = 4 ) get
  integer ( kind = 4 ) have
  logical              initialize
  real    ( kind = 8 ) it_diff
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nearest(batch)
  real    ( kind = 8 ) r(dim_num,n)
  real    ( kind = 8 ) r2(dim_num,n)
  real    ( kind = 8 ) s(dim_num,batch)
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) seed
!
!  Take each generator as the first sample point for its region.
!  This can slightly slow the convergence, but it simplifies the
!  algorithm by guaranteeing that no region is completely missed
!  by the sampling.
!
  energy = 0.0D+00
  r2(1:dim_num,1:n) = r(1:dim_num,1:n)
  count(1:n) = 1
!
!  Generate the sampling points S in batches.
!
  have = 0

  do while ( have < sample_num )

    get = min ( sample_num - have, batch )

    call cvt_sample ( dim_num, sample_num, get, sample, initialize, seed, s )

    initialize = .false.
    have = have + get
!
!  Find the index N of the nearest cell generator to each sample point S.
!
    call find_closest ( dim_num, n, get, s, r, nearest )
!
!  Add S to the centroid associated with generator N.
!
    do j = 1, get
      r2(1:dim_num,nearest(j)) = r2(1:dim_num,nearest(j)) + s(1:dim_num,j)
      energy = energy + sum ( ( r(1:dim_num,nearest(j)) - s(1:dim_num,j) )**2 )
      count(nearest(j)) = count(nearest(j)) + 1
    end do

  end do
!
!  Estimate the centroids.
!
  do j = 1, n
    r2(1:dim_num,j) = r2(1:dim_num,j) / real ( count(j), kind = 8 )
  end do
!
!  Determine the sum of the distances between the old generators 
!  and the estimated centroids.
!
  it_diff = 0.0D+00
  do j = 1, n
    it_diff = it_diff + sqrt ( sum ( ( r2(1:dim_num,j) - r(1:dim_num,j) )**2 ) )
  end do
!
!  Replace the generators by the centroids.
!
  r(1:dim_num,1:n) = r2(1:dim_num,1:n)
!
!  Normalize the discrete energy estimate.
!
  energy = energy / real ( sample_num, kind = 8 ) 

  return
end
subroutine cvt_sample ( dim_num, n, n_now, sample, initialize, seed, r )

!*****************************************************************************80
!
!! CVT_SAMPLE returns sample points.
!
!  Discussion:
!
!    N sample points are to be taken from the region, of spatial
!    dimension DIM_NUM.  In most cases, this region is a unit box.
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
!    23 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of sample points to be generated.
!
!    Input, integer ( kind = 4 ) N_NOW, the number of sample points to be 
!    generated on this call.  N_NOW must be at least 1.
!
!    Input, integer ( kind = 4 ) SAMPLE, specifies how the sampling is done.
!    -1, 'RANDOM', using FORTRAN RANDOM function;
!     0, 'UNIFORM', using a simple uniform RNG;
!     1, 'HALTON', from a Halton sequence;
!     2, 'GRID', points from a grid;
!     3, 'USER', from the "user" routine.
!
!    Input, logical INITIALIZE, is TRUE if the pseudorandom process should be
!    reinitialized.
!
!    Input/output, integer ( kind = 4 ) SEED, the random number seed.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N_NOW), the sample points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n_now

  real    ( kind = 8 ) exponent
  integer ( kind = 4 ), allocatable, dimension ( : ) :: halton_base
  integer ( kind = 4 ), allocatable, dimension ( : ) :: halton_leap
  integer ( kind = 4 ), allocatable, dimension ( : ) :: halton_seed
  integer ( kind = 4 ) halton_step
  integer ( kind = 4 ) i
  logical              initialize
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ngrid
  integer ( kind = 4 ) prime
  real    ( kind = 8 ) r(dim_num,n_now)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_max
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable, dimension ( : ) :: tuple

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

    exponent = 1.0D+00 / real ( dim_num, kind = 8 )
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

  else if ( sample == 3 ) then

    call user ( dim_num, n_now, seed, r )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_SAMPLE - Fatal error!'
    write ( *, '(a,i6,a)' ) '  The value of SAMPLE = ', sample, ' is illegal.'
    stop

  end if

  return
end
subroutine data_read ( file_in_name, dim_num, n, r, success )

!*****************************************************************************80
!
!! DATA_READ reads generator coordinate data from a file.
!
!  Discussion:
!
!    The file is assumed to contain one record per line.
!
!    Records beginning with the '#' character are comments, and are ignored.
!    Blank lines are also ignored.
!
!    Each line that is not ignored is assumed to contain exactly (or at least)
!    M real numbers, representing the coordinates of a point.
!
!    There are assumed to be exactly (or at least) N such records.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) N, the number of points.  The program
!    will stop reading data once N values have been read.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the point coordinates.
!
!    Output, logical SUCCESS, is TRUE if the data was read properly.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  real ( kind = 8 ) r(dim_num,n)
  logical success
  real ( kind = 8 ) x(dim_num)

  success = .true.

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    success = .false.
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( file_in_name )
    return
  end if

  i = 0

  do while ( i < n )

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      success = .false.
      return
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, dim_num, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    i = i + 1

    r(1:dim_num,i) = x(1:dim_num)

  end do

  close ( unit = file_in_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DATA_READ:'
  write ( *, '(a,i6)' ) '  Read coordinate data from file.'

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

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n
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
!    02 August 2004
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

  temp = temp + real ( values(2) - 1, kind = 8 ) /  11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) /  30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) /  23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp                                    /   6.0D+00

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
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the leap vector.
!
!    Output, logical, HALHAM_LEAP_CHECK, true if LEAP is legal.
!
  implicit none

  integer ( kind = 4 ) dim_num

  logical halham_leap_check
  integer ( kind = 4 ) leap(dim_num)

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
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, logical HALHAM_N_CHECK, true if N is legal.
!
  implicit none

  logical halham_n_check
  integer ( kind = 4 ) n

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
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, logical HALHAM_DIM_NUM_CHECK, true if DIM_NUM is legal.
!
  implicit none

  logical halham_dim_num_check
  integer ( kind = 4 ) dim_num

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
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the seed vector.
!
!    Output, logical, HALHAM_SEED_CHECK, true if SEED is legal.
!
  implicit none

  integer ( kind = 4 ) dim_num

  logical halham_seed_check
  integer ( kind = 4 ) seed(dim_num)

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
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!
!    Output, logical HALHAM_STEP_CHECK, true if STEP is legal.
!
  implicit none

  logical halham_step_check
  integer ( kind = 4 ) step

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
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases.
!
!    Output, logical, HALTON_BASE_CHECK, true if BASE is legal.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) base(dim_num)
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
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!    1 <= DIM_NUM is required.
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence.
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!    0 <= STEP is required.
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the Halton sequence index 
!    corresponding to STEP = 0.
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the succesive jumps in the 
!    Halton sequence.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the Halton bases.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the next N elements of the
!    leaped Halton subsequence, beginning with element STEP.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) base(dim_num)
  real ( kind = 8 ) base_inv
  integer ( kind = 4 ) digit(n)
  logical halham_leap_check
  logical halham_n_check
  logical halham_dim_num_check
  logical halham_seed_check
  logical halham_step_check
  logical halton_base_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leap(dim_num)
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) seed2(n)
  integer ( kind = 4 ) step
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
!! I4VEC_TRANSPOSE_PRINT prints an integer vector "transposed".
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
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = 11 ) string
  character ( len = * ) title
  integer ( kind = 4 ) title_len

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
!    Input, integer ( kind = 4 ) N, the index of the desired prime number.
!    In general, is should be true that 0 <= N <= PRIME_MAX.
!    N = -1 returns PRIME_MAX, the index of the largest prime available.
!    N = 0 is legal, returning PRIME = 1.
!
!    Output, integer ( kind = 4 ) PRIME, the N-th prime.  If N is out of range,
!    PRIME is returned as -1.
!
  implicit none

  integer ( kind = 4 ), parameter :: prime_max = 1600

  integer ( kind = 4 ), save :: icall = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( prime_max ) :: npvec
  integer ( kind = 4 ) prime

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
    write ( *, '(a,i6)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i6)' ) '  N should be between 1 and PRIME_MAX =', prime_max
    stop
  end if

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints a real matrix, transposed.
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
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of a real matrix, transposed.
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
!! R8MAT_UNIFORM_01 fills an array with pseudorandom numbers.
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
!  Reference:
!
!    Paul Bratley, Bennett Fox, L E Schrage,
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
!  Parameters:
!
!    Input, integer  ( kind = 4 )M, N, the number of rows and columns in 
!    the array.
!
!    Input/output, integer  ( kind = 4 )SEED, the "seed" value, which should 
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer  ( kind = 4 )m
  integer  ( kind = 4 )n

  integer  ( kind = 4 )i
  integer  ( kind = 4 )j
  integer  ( kind = 4 )k
  integer  ( kind = 4 )seed
  real     ( kind = 8 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875E-10

    end do
  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
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
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_filename
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real      ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

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
!    Input, integer ( kind = 4 ) SEED, the value used to initialize the random 
!    number generator.
!
  implicit none

  logical, parameter :: DEBUG = .false.
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
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
subroutine s_cap ( s )

!*****************************************************************************80
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer nchar
  character ( len = * ) s

  nchar = len_trim ( s )

  do i = 1, nchar

    c = s(i:i)
    call ch_cap ( c )
    s(i:i) = c

  end do

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
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
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
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
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character of S used.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

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
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
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
!    0, no errors occurred.
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
  r = 0.0
  lchar = - 1
  isgn = 1
  rtop = 0.0
  rbot = 1.0
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
        rtop = real ( 10.0, kind = 8 ) * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = real ( 10.0, kind = 8 ) * rtop + real ( ndig, kind = 8 )
        rbot = real ( 10.0, kind = 8 ) * rbot
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
    rexp = 1.0
  else

    if ( jbot == 1 ) then
      rexp = real ( 10.0, kind = 8 )**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = real ( 10.0, kind = 8 )**rexp
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
  character ( len = * ) s

  i = 0

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
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
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
!    Input, integer ( kind = 4 ) M, the maximum entry in any component.
!    M must be greater than 0.
!
!    Input, integer ( kind = 4 ) N, the number of components.
!    N must be greater than 0.
!
!    Input, integer ( kind = 4 ) RANK, indicates the rank of the tuple.
!    Typically, 0 <= RANK < N**M.  Values of RANK greater than
!    N**M are legal and meaningful; they are equivalent to the
!    corresponding value mod (N**M).  If RANK < 0, this indicates 
!    that this is the first call for the given values of (M,N).  
!    Initialization is done, and X is set to a dummy value.
!
!    Output, integer ( kind = 4 ) X(N), the next tuple, or a dummy value if
!    initialization has just been done.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), save, allocatable, dimension ( : ) :: base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  if ( rank < 0 ) then

    if ( m <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
      write ( *, '(a)' ) '  The value M <= 0 is not allowed.'
      write ( *, '(a,i6)' ) '  M = ', m
      stop
    end if

    if ( n <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
      write ( *, '(a)' ) '  The value N <= 0 is not allowed.'
      write ( *, '(a,i6)' ) '  N = ', n
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
subroutine user ( dim_num, n, seed, r )

!*****************************************************************************80
!
!! USER samples points in a user-specified region with given density.
!
!  Discussion:
!
!    This routine can be used to 
!
!    * specify an interesting initial configuration for the data,
!      by specifing that USER be used for initialization (INIT = 3);
!
!    * specify the shape of the computational region, by specifying
!      that sample points are to be generated by this routine, 
!      (SAMPLE = 3) and then returning sample points uniformly at random.
!
!    * specify the distribution or density function, by specifying
!      that sample points are to be generated by this routine, 
!      (SAMPLE = 3 ) and then returning sample points according to a 
!      given probability density function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of sample points desired.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value.  The meaning and
!    use of this variable is up to the user.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), an array of sample points from
!    the region.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) angle(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(dim_num,n)
  real ( kind = 8 ) radius(n)
!
!  We sample points in the unit circle uniformly.
!
  call random_number ( harvest = angle(1:n) )
  angle(1:n) = 2.0D+00 * pi * angle(1:n)

  call random_number ( harvest = radius(1:n) )
  radius(1:n) = sqrt ( radius(1:n) )

  r(1,1:n) = radius(1:n) * cos ( angle(1:n) )
  r(2,1:n) = radius(1:n) * sin ( angle(1:n) )

  return
end
