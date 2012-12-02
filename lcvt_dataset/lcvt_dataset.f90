program main

!*****************************************************************************80
!
!! MAIN is the main program for LCVT_DATASET.
!
!  Discussion:
!
!    LCVT_DATASET generates a Latinized CVT dataset and writes it to a file.
!
!    This program is meant to be used interactively.  It's also
!    possible to prepare a simple input file beforehand and use it
!    in batch mode.
!
!    The program requests input values from the user:
!
!    * DIM_NUM, the spatial dimension;
!    * N, the number of points to generate;
!    * SEED_INIT, a seed to use for random number generation;
!    * INIT, initialize the points:
!      ** file, by reading data from file;
!      ** GRID, picking points from a grid;
!      ** HALTON, from a Halton sequence;
!      ** RANDOM, using FORTRAN RANDOM function;
!      ** UNIFORM, using a simple uniform RNG;
!      ** USER, call the "user" routine;
!    * CVT_IT_NUM, the maximum number of iterations;
!    * SAMPLE, how to conduct the sampling:
!      ** GRID, picking points from a grid;
!      ** HALTON, from a Halton sequence;
!      ** RANDOM, using FORTRAN RANDOM function;
!      ** UNIFORM, using a simple uniform RNG;
!      ** USER, call the "user" routine.
!    * SAMPLE_NUM, the number of sampling points;
!    * BATCH, the number of sampling points to create at one time.
!    * LAT_IT_NUM, the maximum number of iterations;
!    * OUTPUT, a file in which to store the data.
!
!    To indicate that no further computations are desired, it is 
!    enough to input a nonsensical value, such as -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Qiang Du, Vance Faber, Max Gunzburger,
!    Centroidal Voronoi Tessellations: Applications and Algorithms,
!    SIAM Review, 
!    Volume 41, 1999, pages 637-676.
!
  implicit none

  integer ( kind = 4 ) batch
  real ( kind = 8 ) cvt_energy
  real ( kind = 8 ) cvt_it_diff
  integer ( kind = 4 ) cvt_it
  integer ( kind = 4 ) cvt_it_num
  logical, parameter :: DEBUG = .true.
  integer ( kind = 4 ) dim_num
  character ( len = 80 ) :: input_file_name
  integer ( kind = 4 ) init
  character ( len = 80 ) :: init_string
  integer ( kind = 4 ) ios
  real ( kind = 8 ) lat_energy
  integer ( kind = 4 ) lat_it
  integer ( kind = 4 ) lat_it_num
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_total
  character ( len = 80 ) :: output_file_name
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r
  logical reset
  logical s_eqi
  integer ( kind = 4 ) sample
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) :: sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LCVT_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Generate a Latinized CVT dataset.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program is meant to be used interactively.'
  write ( *, '(a)' ) '  It is also possible to prepare a simple input '
  write ( *, '(a)' ) '  file beforehand and use it in batch mode.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program requests input values from the user:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  * DIM_NUM, the spatial dimension,'
  write ( *, '(a)' ) '  * N, the number of points to generate,'
  write ( *, '(a)' ) '  * SEED_INIT, a seed to use for random number generation,'
  write ( *, '(a)' ) '  * INIT, initialize the points:'
  write ( *, '(a)' ) '    ** file, read data from a file;'
  write ( *, '(a)' ) '    ** GRID, by picking points from a grid;'
  write ( *, '(a)' ) '    ** HALTON, from a Halton sequence;'
  write ( *, '(a)' ) '    ** RANDOM, using FORTRAN RANDOM function;'
  write ( *, '(a)' ) '    ** UNIFORM, using a simple uniform RNG;'
  write ( *, '(a)' ) '    ** USER, call the "user" routine;'
  write ( *, '(a)' ) '  * CVT_IT_NUM, the maximum number of CVT iterations.'
  write ( *, '(a)' ) '  * SAMPLE, how to conduct the sampling.'
  write ( *, '(a)' ) '    ** GRID, by picking points from a grid;'
  write ( *, '(a)' ) '    ** HALTON, from a Halton sequence;'
  write ( *, '(a)' ) '    ** RANDOM, using FORTRAN RANDOM function;'
  write ( *, '(a)' ) '    ** UNIFORM, using a simple uniform RNG;'
  write ( *, '(a)' ) '    ** USER, call the "user" routine;'
  write ( *, '(a)' ) '  * SAMPLE_NUM, the number of sample points.'
  write ( *, '(a)' ) '  * BATCH, number of sample points to create at one time.'
  write ( *, '(a)' ) '  * LAT_IT_NUM, the number of Latinizing iterations.'
  write ( *, '(a)' ) '  * OUTPUT, a file in which to store the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  To indicate that no further computations are '
  write ( *, '(a)' ) '  desired, it is enough to input a nonsensical value, '
  write ( *, '(a)' ) '  such as -1.'

  write ( *, '(a)' ) '  *'
  write ( *, '(a)' ) ' *'
  write ( *, '(a)' ) '*  Ready to generate a new dataset:'
  write ( *, '(a)' ) ' *'
  write ( *, '(a)' ) '  *'
  write ( *, '(a)' ) '  Enter DIM_NUM, the spatial dimension:'
  write ( *, '(a)' ) '  (Try "2" if you do not have a preference.)'
  write ( *, '(a)' ) '  (0 or any negative value terminates execution).'

  read ( *, *, iostat = ios ) dim_num

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for DIM_NUM.'
    stop
  end if

  write ( *, '(a,i12)' ) '  User input DIM_NUM = ', dim_num

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET'
    write ( *, '(a,i12)' ) '  The input value of DIM_NUM = ', dim_num
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter N, the number of points to generate:'
  write ( *, '(a)' ) '  (Try "25" if you do not have a preference.)'
  write ( *, '(a)' ) '  (0 or any negative value terminates execution).'

  read ( *, *, iostat = ios ) n

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for N.'
    stop
  end if

  write ( *, '(a,i12)' ) '  User input N = ', n

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET'
    write ( *, '(a,i12)' ) '  The input value of N = ', n
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter SEED_INIT, a seed for the random number generator:'
  write ( *, '(a)' ) '  (Try "123456789" if you do not have a preference.)'
  write ( *, '(a)' ) '  (Any negative value terminates execution).'

  read ( *, *, iostat = ios ) seed_init

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for SEED_INIT.'
    stop
  end if

  write ( *, '(a,i12)' ) '  User input SEED_INIT = ', seed_init

  if ( seed_init < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET'
    write ( *, '(a,i12)' ) '  The input value of SEED_INIT = ', seed_init
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  INIT is the method of initializing the data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  file     read data from a file;'
  write ( *, '(a)' ) '  GRID     by picking points from a grid;'
  write ( *, '(a)' ) '  HALTON   from a Halton sequence;'
  write ( *, '(a)' ) '  RANDOM   using FORTRAN RANDOM function;'
  write ( *, '(a)' ) '  UNIFORM  using a simple uniform RNG;'
  write ( *, '(a)' ) '  USER     call the "user" routine;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Try "RANDOM" if you do not have a preference.)'
  write ( *, '(a)' ) '  (A blank value terminates execution).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter INIT:'
  read ( *, '(a)', iostat = ios ) init_string

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for INIT.'
    stop
  end if

  if ( s_eqi ( init_string(1:6), 'RANDOM'  ) ) then
    init = -1
    write ( *, '(a)' ) '  User input INIT = "RANDOM".'
  else if ( s_eqi ( init_string(1:7), 'UNIFORM' ) ) then
    init = 0
    write ( *, '(a)' ) '  User input INIT = "UNIFORM".'
  else if ( s_eqi ( init_string(1:6), 'HALTON'  ) ) then
    init = 1
    write ( *, '(a)' ) '  User input INIT = "HALTON".'
  else if ( s_eqi ( init_string(1:4), 'GRID'    ) ) then
    init = 2
    write ( *, '(a)' ) '  User input INIT = "GRID".'
  else if ( s_eqi ( init_string(1:4), 'USER'    ) ) then
    init = 3
    write ( *, '(a)' ) '  User input INIT = "USER".'
  else if ( 0 < len_trim ( init_string ) ) then
    init = 4
    write ( *, '(a)' ) '  User input INIT = FILE_NAME = "RANDOM".'
    input_file_name = init_string
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET'
    write ( *, '(a)' ) '  The input value of INIT '
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CVT_IT_NUM is the number of CVT iterations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A CVT iteration carries out the following steps:'
  write ( *, '(a)' ) '  * the Voronoi region associated with each'
  write ( *, '(a)' ) '    generator is estimated by sampling;'
  write ( *, '(a)' ) '  * the centroid of each Voronoi region is estimated.'
  write ( *, '(a)' ) '  * the generator is replaced by the centroid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If "enough" sampling points are used,'
  write ( *, '(a)' ) '  and "enough" iterations are taken, this process'
  write ( *, '(a)' ) '  will converge.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Try "50" if you do not have a preference.)'
  write ( *, '(a)' ) '  (A negative value terminates execution).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter CVT_IT_NUM:'
  read ( *, *, iostat = ios ) cvt_it_num

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for CVT_IT_NUM.'
    stop
  end if

  write ( *, '(a,i6)' ) '  User input CVT_IT_NUM = ', cvt_it_num

  if ( cvt_it_num < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET'
    write ( *, '(a,i12)' ) '  The input value of CVT_IT_NUM = ', cvt_it_num
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
  write ( *, '(a)' ) '  USER     call the "user" routine;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Try "RANDOM" if you do not have a preference.)'
  write ( *, '(a)' ) '  (A blank value terminates execution).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter SAMPLE:'
  read ( *, '(a)', iostat = ios ) sample_string

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for SAMPLE.'
    stop
  end if

  if ( s_eqi ( sample_string(1:6), 'RANDOM'  ) ) then
    write ( *, '(a)' ) '  User input INIT = "RANDOM".'
    sample = -1
  else if ( s_eqi ( sample_string(1:7), 'UNIFORM' ) ) then
    write ( *, '(a)' ) '  User input INIT = "UNIFORM".'
    sample = 0
  else if ( s_eqi ( sample_string(1:6), 'HALTON'  ) ) then
    write ( *, '(a)' ) '  User input INIT = "HALTON".'
    sample = 1
  else if ( s_eqi ( sample_string(1:4), 'GRID'    ) ) then
    write ( *, '(a)' ) '  User input INIT = "GRID".'
    sample = 2
  else if ( s_eqi ( sample_string(1:4), 'USER' ) ) then
    write ( *, '(a)' ) '  User input INIT = "USER".'
    sample = 3
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET'
    write ( *, '(a)' ) '  The input value of SAMPLE '
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SAMPLE_NUM is the number of sample points for CVT.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Voronoi regions will be explored by generating'
  write ( *, '(a)' ) '  SAMPLE_NUM points.  For each sample point, the'
  write ( *, '(a)' ) '  nearest generator is found.  Using more points'
  write ( *, '(a)' ) '  gives a better estimate of these regions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SAMPLE_NUM should be much larger than N, the'
  write ( *, '(a)' ) '  number of generators.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Try "10000" if you do not have a preference.)'
  write ( *, '(a)' ) '  (A zero or negative value terminates execution.)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter SAMPLE_NUM:'

  read ( *, *, iostat = ios ) sample_num

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for SAMPLE_NUM.'
    stop
  end if

  write ( *, '(a,i12)' ) '  User input SAMPLE_NUM = ', sample_num

  if ( sample_num <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET'
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
    write ( *, '(a)' ) 'LCVT_DATASET - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for SAMPLE_NUM.'
    stop
  end if

  write ( *, '(a,i12)' ) '  User input BATCH = ', batch

  if ( batch <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET'
    write ( *, '(a,i12)' ) '  The input value of BATCH = ', batch
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LAT_IT_NUM is the number of Latinizing iterations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Each step of the latinizing iteration begins'
  write ( *, '(a)' ) '  by carrying out CVT_IT_NUM steps of CVT iteration,'
  write ( *, '(a)' ) '  after which the data is "latinized".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Often, one latinizing step is enough.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In some cases, it may be worth while to carry'
  write ( *, '(a)' ) '  out several latinizing steps; that is, the'
  write ( *, '(a)' ) '  Latinized data is smoothed by another series'
  write ( *, '(a)' ) '  of CVT steps, then latinized, and so on.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Try "1" if you do not have a preference.)'
  write ( *, '(a)' ) '  (A negative value terminates execution).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter LAT_IT_NUM:'
  read ( *, *, iostat = ios ) lat_it_num

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for LAT_IT_NUM.'
    stop
  end if

  write ( *, '(a,i6)' ) '  User input LAT_IT_NUM = ', lat_it_num

  if ( cvt_it_num < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET'
    write ( *, '(a,i12)' ) '  The input value of LAT_IT_NUM = ', lat_it_num
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  OUTPUT is the name of a file into which'
  write ( *, '(a)' ) '  the computed data may be stored.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Try "lcvt.txt" if you do not have a preference.)'
  write ( *, '(a)' ) '  (A blank value terminates execution).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter OUTPUT:'
  read ( *, '(a)', iostat = ios ) output_file_name

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET - Warning!'
    write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
    write ( *, '(a)' ) '  while expecting input for OUTPUT.'
    stop
  end if

  write ( *, '(a)' ) '  User input OUTPUT = "' // trim ( output_file_name ) &
    // '".'

  if ( len_trim ( output_file_name ) <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_DATASET'
    write ( *, '(a)' ) '  The input value of OUTPUT '
    write ( *, '(a)' ) '  is interpreted as a request for termination.'
    write ( *, '(a)' ) '  Normal end of execution.'
    stop
  end if
!
!  Initialize the data.
!
  allocate ( r(1:dim_num,1:n) )

  if ( init == 4 ) then

    call r8table_data_read ( input_file_name, dim_num, n, r )

  else

    n_total = n
    reset = .true.

    call region_sampler ( dim_num, n, n_total, r, init, reset, seed )

  end if

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Latin IT      CVT Energy    Latin Energy'
    write ( *, '(a)' ) ' '
  end if

  do lat_it = 1, lat_it_num

    if ( .false. ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    CVT IT  Change'
      write ( *, '(a)' ) ' '
    end if

    do cvt_it = 1, cvt_it_num

      call cvt_iteration ( dim_num, n, r, sample_num, &
        sample, seed, cvt_it_diff )

      if ( .false. ) then
        write ( *, '(2x,i8,2x,g14.6)' ) cvt_it, cvt_it_diff
      end if

    end do

    if ( .false. ) then
      write ( *, '(a)' ) ' '
    end if

    call cluster_energy ( dim_num, n, r, sample_num, &
      sample, seed, cvt_energy )

    call r8mat_latinize ( dim_num, n, r )

    call cluster_energy ( dim_num, n, r, sample_num, &
      sample, seed, lat_energy )

    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) &
      lat_it, cvt_energy, lat_energy

  end do
!
!  Write the data to a file.
!
  call lcvt_write ( dim_num, n, seed_init, init, &
    input_file_name, sample, sample_num, &
    cvt_it_num, cvt_energy, lat_it_num, lat_energy, &
    r, output_file_name )

  deallocate ( r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data was written to the file "' &
     // trim ( output_file_name ) // '".'


  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LCVT_DATASET:'
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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
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
subroutine cluster_energy ( dim_num, n, cell_generator, sample_num, &
  sample_function_cvt, seed, energy )

!******************************************************************************
!
!! CLUSTER_ENERGY returns the energy of a dataset.
!
!  Discussion:
!
!    The energy is the integral of the square of the distance from each point
!    in the region to its nearest generator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of Voronoi cells.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,N), the Voronoi
!    cell generators.  On output, these have been modified
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of sample points.
!
!    Input, integer ( kind = 4 ) SAMPLE_FUNCTION_CVT, specifies how the region is sampled:
!    -1, the sampling function is RANDOM_NUMBER (Fortran90 intrinsic),
!    0, the sampling function is UNIFORM,
!    1, the sampling function is HALTON,
!    2, the sampling function is GRID.
!
!    Input/output, integer ( kind = 4 ) SEED, the random number seed.
!
!    Output, real ( kind = 8 ) ENERGY, the (estimated) energy of the dataset.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) cell_generator(dim_num,n)
  real ( kind = 8 ) energy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nearest
  logical reset
  integer ( kind = 4 ) sample_function_cvt
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(dim_num)

  energy = 0.0D+00
  reset = .true.

  do j = 1, sample_num
!
!  Generate a sampling point X.
!
    call region_sampler ( dim_num, 1, sample_num, x, &
     sample_function_cvt, reset, seed )

    reset = .false.
!
!  Find the nearest cell generator.
!
    call find_closest ( dim_num, n, x, cell_generator, nearest )
!
!  Add the contribution to the energy.
!
    energy = energy &
      + sum ( cell_generator(1:dim_num,nearest) - x(1:dim_num) )**2

  end do

  energy = energy / real ( sample_num, kind = 8 )

  return
end
subroutine cvt_iteration ( dim_num, n, r, sample_num, &
  sample_function_cvt, seed, change_l2 )

!******************************************************************************
!
!! CVT_ITERATION takes one step of the CVT iteration.
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
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of Voronoi cells.
!
!    Input/output, real ( kind = 8 ) R(DIM_NUM,N), the Voronoi
!    cell generators.  On output, these have been modified
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of sample points.
!
!    Input, integer ( kind = 4 ) SAMPLE_FUNCTION_CVT, specifies how the region is sampled:
!    -1, the sampling function is RANDOM_NUMBER (Fortran90 intrinsic),
!    0, the sampling function is UNIFORM,
!    1, the sampling function is HALTON,
!    2, the sampling function is GRID,
!
!    Input/output, integer ( kind = 4 ) SEED, the random number seed.
!
!    Output, real ( kind = 8 ) CHANGE_L2, the L2 norm of the difference between
!    the input and output data.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) change_l2
  integer ( kind = 4 ) count(n)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nearest
  real ( kind = 8 ) r(dim_num,n)
  real ( kind = 8 ) r2(dim_num,n)
  logical reset
  integer ( kind = 4 ) sample_function_cvt
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(dim_num)

  r2(1:dim_num,1:n) = 0.0D+00
  count(1:n) = 0
  reset = .true.

  do j = 1, sample_num
!
!  Generate a sampling point X.
!
    call region_sampler ( dim_num, 1, sample_num, x, sample_function_cvt, &
      reset, seed )

    reset = .false.
!
!  Find the nearest cell generator.
!
    call find_closest ( dim_num, n, x, r, nearest )
!
!  Add X to the averaging data for R(*,NEAREST).
!
    r2(1:dim_num,nearest) = r2(1:dim_num,nearest) + x(1:dim_num)

    count(nearest) = count(nearest) + 1

  end do
!
!  Compute the new generators.
!
  do j = 1, n
    if ( count(j) /= 0 ) then
      r2(1:dim_num,j) = r2(1:dim_num,j) / real ( count(j), kind = 8 )
    end if
  end do
!
!  Determine the change.
!
  change_l2 = sqrt ( sum ( ( r2(1:dim_num,1:n) - r(1:dim_num,1:n) )**2 ) )
!
!  Update.
!
  r(1:dim_num,1:n) = r2(1:dim_num,1:n)

  return
end
subroutine file_column_count ( file_in_name, ncolumn )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of NCOLUMN words, separated
!    by spaces.  There may also be some blank lines, and some comment lines,
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
!    Input, character ( len = * ) FILE_IN_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NCOLUMN, the number of columns assumed to be in the file.
!
  implicit none

  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  logical got_one
  integer ( kind = 4 ) ios
  character ( len = 256 ) line
  integer ( kind = 4 ) ncolumn
!
!  Open the file.
!
  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    ncolumn = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( file_in_name )
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( file_in_unit, '(a)', iostat = ios ) line

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

    rewind ( file_in_unit )

    do

      read ( file_in_unit, '(a)', iostat = ios ) line

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

  close ( unit = file_in_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    ncolumn = 0
    return
  end if

  call s_word_count ( line, ncolumn )

  return
end
subroutine file_row_count ( file_in_name, row_num )

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
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 100 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( file_in_name )
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( file_in_unit, '(a)', iostat = ios ) line

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

  close ( unit = file_in_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_ROW_COUNT:'
  write ( *, '(a,i8)' ) '  Number of records:         ', record_num
  write ( *, '(a,i8)' ) '  Number of data records:    ', row_num
  write ( *, '(a,i8)' ) '  Number of comment records: ', comment_num

  return
end
subroutine find_closest ( dim_num, n, x, generator, nearest )

!******************************************************************************
!
!! FIND_CLOSEST finds the Voronoi cell generator closest to a point X.
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
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of cell generatorrs.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the point to be checked.
!
!    Input, real ( kind = 8 ) GENERATOR(DIM_NUM,N), the cell generators.
!
!    Output, integer ( kind = 4 ) NEAREST, the index of the nearest cell generators.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) dist_sq
  real ( kind = 8 ) dist_sq_min
  real ( kind = 8 ) generator(dim_num,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nearest
  real ( kind = 8 ) x(dim_num)

  nearest = -1
  dist_sq_min = huge ( dist_sq_min )

  do i = 1, n

    dist_sq = sum ( ( generator(1:dim_num,i) - x(1:dim_num) )**2 )

    if ( dist_sq < dist_sq_min ) then
      dist_sq_min = dist_sq
      nearest = i
    end if

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

  temp = temp + real ( values(2) - 1, kind = 8 ) / real (  11.0, kind = 8 )
  temp = temp + real ( values(3) - 1, kind = 8 ) / real (  30.0, kind = 8 )
  temp = temp + real ( values(5),     kind = 8 ) / real (  23.0, kind = 8 )
  temp = temp + real ( values(6),     kind = 8 ) / real (  59.0, kind = 8 )
  temp = temp + real ( values(7),     kind = 8 ) / real (  59.0, kind = 8 )
  temp = temp + real ( values(8),     kind = 8 ) / real ( 999.0, kind = 8 )
  temp = temp                                    / real (   6.0, kind = 8 )

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
subroutine i4_to_halton ( seed, base, m, r )

!*****************************************************************************80
!
!! I4_TO_HALTON computes an element of a vector Halton sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J H Halton,
!    Numerische Mathematik,
!    Volume 2, pages 84-90.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED, the index of the desired element.
!    Only the absolute value of SEED is considered.  SEED = 0 is allowed,
!    and returns R = 0.
!
!    Input, integer ( kind = 4 ) BASE(M), the Halton bases, which should be
!    distinct prime numbers.  This routine only checks that each base
!    is greater than 1.
!
!    Input, integer ( kind = 4 ) M, the dimension of the sequence.
!
!    Output, real ( kind = 8 ) R(M), the SEED-th element of the Halton
!    sequence for the given bases.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) base(m)
  real ( kind = 8 ) base_inv(m)
  integer ( kind = 4 ) digit(m)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r(m)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed2(m)

  seed2(1:m) = abs ( seed )

  r(1:m) = 0.0E+00

  if ( any ( base(1:m) <= 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
    write ( *, '(a)' ) '  An input base BASE is <= 1!'
    do i = 1, m
      write ( *, '(2i8)' ) i, base(i)
    end do
    stop
  end if

  base_inv(1:m) = 1.0E+00 / real ( base(1:m), kind = 8 )

  do while ( any ( seed2(1:m) /= 0 ) )
    digit(1:m) = mod ( seed2(1:m), base(1:m) )
    r(1:m) = r(1:m) + real ( digit(1:m), kind = 8 ) * base_inv(1:m)
    base_inv(1:m) = base_inv(1:m) / real ( base(1:m), kind = 8 )
    seed2(1:m) = seed2(1:m) / base(1:m)
  end do

  return
end
subroutine lcvt_write ( dim_num, n, seed_init, sample_function_init, &
  file_in_name, sample_function_cvt, sample_num, cvt_it_num, &
  cvt_energy, latin_it_num, latin_energy, cell_generator, file_out_name )

!*****************************************************************************80
!
!! LCVT_WRITE writes a Latinized CVT dataset to a file.
!
!  Discussion:
!
!    The initial lines of the file are comments, which begin with a
!    '#' character.
!
!    Thereafter, each line of the file contains the M-dimensional
!    components of a Latinized CVT generator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) SEED_INIT, the initial random number seed.
!
!    Input, integer ( kind = 4 ) SAMPLE_FUNCTION_INIT, specifies how the initial
!    generators are chosen:
!    -1, the initialization function is RANDOM_NUMBER (Fortran90 intrinsic),
!    0, the initialization function is UNIFORM,
!    1, the initialization function is HALTON,
!    2, the initialization function is GRID,
!    3, the initial values are read in from a file.
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the file
!    from which initialization values were read for the generators,
!    if SAMPLE_FUNCTION_INIT = 3.
!
!    Input, integer ( kind = 4 ) SAMPLE_FUNCTION_CVT, specifies how the region is sampled:
!    -1, the sampling function is RANDOM_NUMBER (Fortran90 intrinsic),
!    0, the sampling function is UNIFORM,
!    1, the sampling function is HALTON,
!    2, the sampling function is GRID.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of sampling points used on
!    each CVT iteration.
!
!    Input, integer ( kind = 4 ) CVT_IT_NUM, the number of CVT iterations.
!
!    Input, real ( kind = 8 ) CVT_ENERGY, the energy of the final CVT dataset.
!
!    Input, integer ( kind = 4 ) LATIN_IT_NUM, the number of Latin iterations.
!
!    Input, real ( kind = 8 ) LATIN_ENERGY, the energy of the Latinized 
!    CVT dataset.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of
!    the output file.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) cell_generator(dim_num,n)
  real ( kind = 8 ) cvt_energy
  integer ( kind = 4 ) cvt_it_num
  character ( len = * ) file_in_name
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  real ( kind = 8 ) latin_energy
  integer ( kind = 4 ) latin_it_num
  integer ( kind = 4 ) sample_function_cvt
  integer ( kind = 4 ) sample_function_init
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed_init
  character ( len = 40 ) string

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, &
    status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCVT_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  call timestring ( string )

  write ( file_out_unit, '(a)' ) '#  ' // trim ( file_out_name )
  write ( file_out_unit, '(a)' ) &
    '#  created by routine LCVT_WRITE.F90'
  write ( file_out_unit, '(a)' ) '#  at ' // trim ( string )
  write ( file_out_unit, '(a)' ) '#'
  write ( file_out_unit, '(a,i8)' ) '#  Spatial dimension DIM_NUM = ', dim_num
  write ( file_out_unit, '(a,i8)' ) '#  Number of points N = ', n
  write ( file_out_unit, '(a,g14.6)' ) '#  EPSILON (unit roundoff) = ', &
    epsilon ( cell_generator(1,1) )

  if ( sample_function_init == 0 .or. &
       sample_function_init == 1 .or. &
       sample_function_cvt == 0 .or. &
       sample_function_cvt == 1 ) then
    write ( file_out_unit, '(a)' ) '#'
    write ( file_out_unit, '(a,i12)' ) '#  Initial SEED = ', seed_init
  end if

  write ( file_out_unit, '(a)' ) '#'
  if ( sample_function_init == -1 ) then
    write ( file_out_unit, '(a)' ) &
      '#  Initialization by RANDOM_NUMBER (Fortran90 intrinsic).'
  else if ( sample_function_init == 0 ) then
    write ( file_out_unit, '(a)' ) '#  Initialization by UNIFORM.'
  else if ( sample_function_init == 1 ) then
    write ( file_out_unit, '(a)' ) '#  Initialization by HALTON.'
  else if ( sample_function_init == 2 ) then
    write ( file_out_unit, '(a)' ) '#  Initialization by GRID.'
  else if ( sample_function_init == 3 ) then
    write ( file_out_unit, '(a)' ) '#  Initialization from file: "' // &
      trim ( file_in_name ) // '".'
  end if

  if ( sample_function_cvt == -1 ) then
    write ( file_out_unit, '(a)' ) &
      '#  Sampling by RANDOM_NUMBER (Fortran90 intrinsic).'
  else if ( sample_function_cvt == 0 ) then
    write ( file_out_unit, '(a)' ) '#  Sampling by UNIFORM.'
  else if ( sample_function_cvt == 1 ) then
    write ( file_out_unit, '(a)' ) '#  Sampling by HALTON.'
  else if ( sample_function_cvt == 2 ) then
    write ( file_out_unit, '(a)' ) '#  Sampling by GRID.'
  end if

  write ( file_out_unit, '(a,i12)' ) &
    '#  Number of sample points = ', sample_num
  write ( file_out_unit, '(a,i8)' ) &
    '#  Number of CVT iterations = ', cvt_it_num
  write ( file_out_unit, '(a,g14.6)' ) &
    '#  Energy of CVT dataset = ', cvt_energy
  write ( file_out_unit, '(a,i8)' ) &
    '#  Number of Latin iterations = ', latin_it_num
  write ( file_out_unit, '(a,g14.6)' ) &
    '#  Energy of Latinized CVT dataset = ', latin_energy
  write ( file_out_unit, '(a)' ) '#'

  write ( string, '(a,i3,a)' ) '(', dim_num, '(2x,f10.6))'
  do j = 1, n
    write ( file_out_unit, string ) cell_generator(1:dim_num,j)
  end do

  close ( unit = file_out_unit )

  return
end
subroutine param_print ( dim_num, n, cvt_it_num, latin_it_num, seed, &
  seed_init, sample_function_cvt, sample_function_init, sample_num )

!*****************************************************************************80
!
!! PARAM_PRINT prints the program parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of Voronoi cells.
!
!    Input, integer ( kind = 4 ) CVT_IT_NUM, the number of CVT iterations.
!
!    Input, integer ( kind = 4 ) LATIN_IT_NUM, the number of Latin iterations.
!
!    Input, integer ( kind = 4 ) SEED, the current random number seed.
!
!    Input, integer ( kind = 4 ) SEED_INIT, the initial random number seed.
!
!    Input, integer ( kind = 4 ) SAMPLE_FUNCTION_CVT, specifies how the region is sampled:
!    -1, the sampling function is RANDOM_NUMBER (Fortran90 intrinsic),
!    0, the sampling function is UNIFORM,
!    1, the sampling function is HALTON,
!    2, the sampling function is GRID.
!
!    Input, integer ( kind = 4 ) SAMPLE_FUNCTION_INIT, specifies how the initial
!    generators are chosen:
!    0, the initialization function is UNIFORM,
!    1, the initialization function is HALTON,
!    2, the initialization function is GRID,
!    3, the initial values are read in from a file.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of sample points for the
!    CVT iteration.
!
  implicit none

  integer ( kind = 4 ) cvt_it_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) latin_it_num
  integer ( kind = 4 ) n
  integer ( kind = 4 ) sample_function_cvt
  integer ( kind = 4 ) sample_function_init
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Geometry parameters:'
  write ( *, '(a)' ) '-------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The spatial dimension is DIM_NUM = ', dim_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data lies in the unit hypercube [0,1]^DIM_NUM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Latin hypercube parameters:'
  write ( *, '(a)' ) '-------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Number of Latin iterations: ', latin_it_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT Algorithm parameters:'
  write ( *, '(a)' ) '-------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of Voronoi cells to generate: ', n

  write ( *, '(a,i12)' ) '  Number of CVT iterations: ', cvt_it_num
  write ( *, '(a,i12)' ) '  Number of sampling points for CVT iteration: ', &
    sample_num

  if ( sample_function_init == -1 ) then
    write ( *, '(a)' ) &
      '  The generators are initialized by RANDOM_NUMBER (Fortran90 intrinsic).'
  else if ( sample_function_init == 0 ) then
    write ( *, '(a)' ) &
      '  The generators are initialized by UNIFORM.'
  else if ( sample_function_init == 1 ) then
    write ( *, '(a)' ) &
      '  The generators are initialized by HALTON.'
  else if ( sample_function_init == 2 ) then
    write ( *, '(a)' ) &
      '  The generators are initialized by GRID.'
  else if ( sample_function_init == 3 ) then
    write ( *, '(a)' ) &
      '  The generators are initialized from a file.'
  end if

  if ( sample_function_cvt == -1 ) then
    write ( *, '(a)' ) &
      '  The region is sampled by RANDOM_NUMBER (Fortran90 intrinsic).'
  else if ( sample_function_cvt == 0 ) then
    write ( *, '(a)' ) '  The region is sampled by UNIFORM.'
  else if ( sample_function_cvt == 1 ) then
    write ( *, '(a)' ) '  The region is sampled by HALTON.'
  else if ( sample_function_cvt == 2 ) then
    write ( *, '(a)' ) '  The region is sampled by GRID.'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Miscellaneous parameters:'
  write ( *, '(a)' ) '------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial random number seed was ', seed_init
  write ( *, '(a,i12)' ) '  The current random number seed is  ', seed


  write ( *, '(a)' ) ' '

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
    write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i8)' ) '  N should be between 1 and PRIME_MAX =', prime_max
    stop
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
!    Input/output, integer ( kind = 4 ) ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) ( kind = 4 ) seed

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
subroutine r8mat_latinize ( dim_num, n, table )

!******************************************************************************
!
!! R8MAT_LATINIZE "Latinizes" an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of values of type "real ( kind = 8 )".
!
!    It is assumed, though not necessary, that the input dataset
!    has points that lie in the unit hypercube.
!
!    In any case, the output dataset will have this property.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of cells.
!
!    Input/output, real ( kind = 8 ) TABLE(DIM_NUM,N).  On input, the dataset to
!    be "Latinized".  On output, the Latinized dataset.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) table(dim_num,n)

  do i = 1, dim_num
    call r8vec_sort_heap_index_a ( n, table(i,1:n), indx )
    do j = 1, n
      table(i,indx(j)) = real ( 2 * j - 1, kind = 8 ) / real ( 2 * n, kind = 8 )
    end do
  end do

  return
end
subroutine r8table_data_read ( input_file_name, m, n, table )

!*****************************************************************************80
!
!! R8TABLE_DATA_READ reads data from an R8TABLE file.
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
!    26 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
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
  character ( len = * ) input_file_name
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  real ( kind = 8 ) table(m,n)
  real ( kind = 8 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8TABLE_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_file_name ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8TABLE_DATA_READ - Fatal error!'
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
subroutine r8vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of values of type "real ( kind = 8 )".
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call DVEC_PERMUTE ( N, A, INDX )
!
!    after which A(I), I = 1 to N is sorted.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
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

  do i = 1, n
    indx(i) = i
  end do

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
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
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
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
subroutine region_sampler ( dim_num, n, n_total, r, sample_function, reset, &
  seed )

!******************************************************************************
!
!! REGION_SAMPLER returns a sample point in the physical region.
!
!  Discussion:
!
!    This routine original interfaced with a lower routine called
!    TEST_REGION, which tested whether the points generated in the
!    bounding box were actually inside a possibly smaller physical
!    region of interest.  It's been a long time since that option
!    was actually used, so it's been dropped.
!
!    A point is chosen in the bounding box, either by a uniform random
!    number generator, or from a vector Halton sequence.
!
!    The original coding for this routine only supported a Halton
!    sequence of dimension 3 or less.  This restriction has been removed.
!
!    Note that RESET was made an input-only quantity, in part to match
!    the behavior of the routine in MATLAB, where it's cumbersome to
!    support an input/output variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points to generate now.
!
!    Input, integer ( kind = 4 ) N_TOTAL, the total number of points to generate.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the sample points.
!
!    Input, integer ( kind = 4 ) SAMPLE_FUNCTION, specifies how the region is sampled:
!    -1, the sampling function is RANDOM_NUMBER (Fortran90 intrinsic),
!    0, the sampling function is UNIFORM,
!    1, the sampling function is HALTON,
!    2, the sampling function is GRID,
!    3, sample points are generated elsewhere, and this routine is skipped.
!
!    Input, logical RESET, if TRUE, then internal data should be reset
!    for a new problem.  Set RESET to TRUE on the first call to this
!    routine.  Also set RESET to TRUE if you are using the Halton sampler,
!    and you want to reset the Halton seed to 1.
!
!    Input/output, integer ( kind = 4 ) SEED, the random number seed.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) exponent
  integer ( kind = 4 ), save, allocatable, dimension ( : ) :: halton_base
  integer ( kind = 4 ), save :: halton_seed = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: dim_num_save = 0
  integer ( kind = 4 ) n_total
  integer ( kind = 4 ), save :: ngrid
  integer ( kind = 4 ) prime
  real ( kind = 8 ) r(dim_num,n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ), save :: rank
  logical reset
  integer ( kind = 4 ) sample_function
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) tuple(dim_num)

  if ( sample_function == -1 ) then

    call random_number ( harvest = r(1:dim_num,1:n) )

  else if ( sample_function == 0 ) then

    call random_number ( harvest = r(1:dim_num,1:n) )

    do j = 1, n
      do i = 1, dim_num
!       r(i,j) = r8_uniform_01 ( seed )
      end do
    end do

  else if ( sample_function == 1 ) then

    if ( reset ) then
      halton_seed = 1
    end if

    if ( dim_num_save < dim_num ) then

      if ( allocated ( halton_base ) ) then
        deallocate ( halton_base )
      end if
      allocate ( halton_base(1:dim_num) )

      dim_num_save = dim_num

      do i = 1, dim_num
        halton_base(i) = prime ( i )
      end do

    end if

    do j = 1, n
      call i4_to_halton ( halton_seed, halton_base, dim_num, r(1:dim_num,j) )
      halton_seed = halton_seed + 1
    end do

  else if ( sample_function == 2 ) then

    if ( reset ) then

      rank = 0
      exponent = 1.0D+00 / real ( dim_num, kind = 8 )
      ngrid = int ( ( real ( n_total, kind = 8 ) )**exponent )

      if ( ngrid**dim_num < n_total ) then
        ngrid = ngrid + 1
      end if

    end if

    do j = 1, n
      call tuple_next_fast ( ngrid, dim_num, rank, tuple )
      rank = rank + 1
      r(1:dim_num,j) = real ( 2 * tuple(1:dim_num) - 1, kind = 8 ) &
                     / real ( 2 * ngrid, kind = 8 ) 
    end do

  else if ( sample_function == 3 ) then

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REGION_SAMPLER - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal SAMPLE_FUNCTION = ', sample_function
    stop

  end if

  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 July 1998
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
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  iput = 0
  nchar = len_trim ( s )

  do iget = 1, nchar

    c = s(iget:iget)

    if ( c /= ' ' .and. c /= TAB ) then
      iput = iput + 1
      s(iput:iput) = c
    end if

  end do

  s(iput+1:nchar) = ' '

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
  integer ( kind = 4 ) nchar
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
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 value from a string.
!
!  Discussion:
!
!    An "R8" value is simply a real number to be stored as a
!    variable of type "real ( kind = 8 )".
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
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
!    07 September 2004
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
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( s_length < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

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
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
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
!  Scientific notation exponent marker.
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
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

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
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
  if ( iterm /= 1 .and. length+1 == s_length ) then
    length = s_length
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, r8vec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Discussion:
!
!    An R8VEC is a vector of real values, of type "real ( kind = 8 )".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2005
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
!    Output, real ( kind = 8 ) R8VEC(N), the values read from the string.
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
  real ( kind = 8 ) r8vec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), r8vec(i), ierror, lchar )

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
subroutine test_region ( x, dim_num, value )

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
!    05 March 2003
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
!    Output, integer ( kind = 4 ) VALUE, indicates the status of the point:
!    -1: the point is on the boundary of the region.
!     0: the point is outside the region.
!    +1: the point is inside the region.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) value
  real ( kind = 8 ) x(dim_num)

  value = 1

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
  integer ( kind = 4 )   ( kind = 4 ) d
  integer ( kind = 4 )   ( kind = 4 ) h
  integer ( kind = 4 )   ( kind = 4 ) m
  integer ( kind = 4 )   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 )   ( kind = 4 ) n
  integer ( kind = 4 )   ( kind = 4 ) s
  integer ( kind = 4 )   ( kind = 4 ) values(8)
  integer ( kind = 4 )   ( kind = 4 ) y

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
  character ( len = * ) string
  character ( len = 10 ) time
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
!    Rank  X      X
!    ----  ---    ---
!    0     * *    1 1
!    1     1 1    1 2
!    2     1 2    1 3
!    3     1 3    2 1
!    4     2 1    2 2
!    5     2 2    2 3
!    6     2 3    3 1
!    7     3 1    3 2
!    8     3 2    3 3
!    9     3 3    1 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 April 2003
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the maximum entry.
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, integer ( kind = 4 ) RANK, indicates the rank of the tuples.
!    On the very first call, it is necessary that the user set RANK = 0.  
!
!    Output, integer ( kind = 4 ) X(N), the tuple of rank RANK.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), save, allocatable, dimension ( : ) :: base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  if ( rank == 0 ) then

    if ( allocated ( base ) ) then
      deallocate ( base )
    end if
    allocate ( base(1:n) )
    base(n) = 1
    do i = n-1, 1, -1
      base(i) = base(i+1) * m
    end do

  end if

  x(1:n) = mod ( rank / base(1:n), m ) + 1

  return
end
