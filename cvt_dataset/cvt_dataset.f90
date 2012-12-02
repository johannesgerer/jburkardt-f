program main

!*****************************************************************************80
!
!! MAIN is the main program for CVT_DATASET.
!
!  Discussion:
!
!    CVT_DATASET generates a CVT dataset and writes it to a file.
!
!    This program is meant to be used interactively.  It's also
!    possible to prepare a simple input file beforehand and use it
!    in batch mode.
!
!    The program requests input values from the user:
!
!    * DIM_NUM, the spatial dimension;
!    * N, the number of points to generate;
!    * SEED, a seed to use for random number generation;
!    * INIT, initialize the points:
!      ** file, by reading data from file;
!      ** GRID, picking points from a grid;
!      ** HALTON, from a Halton sequence;
!      ** RANDOM, using FORTRAN RANDOM function;
!      ** UNIFORM, using a simple uniform RNG;
!      ** USER, call the "user" routine;
!    * IT_MAX, the maximum number of iterations;
!    * IT_FIXED, the number of iterative steps to take
!      using a fixed set of sampling points.
!    * SAMPLE, how to conduct the sampling:
!      ** GRID, picking points from a grid;
!      ** HALTON, from a Halton sequence;
!      ** RANDOM, using FORTRAN RANDOM function;
!      ** UNIFORM, using a simple uniform RNG;
!      ** USER, call the "user" routine.
!    * SAMPLE_NUM, the number of sampling points;
!    * BATCH, the number of sampling points to create at one time.
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
!    04 July 2005
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

  integer batch
  logical, parameter :: comment = .false.
  logical, parameter :: DEBUG = .true.
  real ( kind = 8 ) energy
  character ( len = 80 ) :: file_out_name
  integer init
  character ( len = 80 ) :: init_string
  integer ios
  real ( kind = 8 ) it_diff
  integer it_fixed
  integer it_max
  integer it_num
  integer n
  integer dim_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r
  logical s_eqi
  integer sample
  integer sample_num
  character ( len = 80 ) :: sample_string
  integer seed
  integer seed_init
  logical success

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Generate a CVT dataset.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program is meant to be used interactively.'
  write ( *, '(a)' ) '  It is also possible to prepare a simple input '
  write ( *, '(a)' ) '  file beforehand and use it in batch mode.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program requests input values from the user:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  * DIM_NUM, the spatial dimension,'
  write ( *, '(a)' ) '  * N, the number of points to generate,'
  write ( *, '(a)' ) '  * SEED, a seed to use for random number generation,'
  write ( *, '(a)' ) '  * INIT, initialize the points:'
  write ( *, '(a)' ) '    ** file, read data from a file;'
  write ( *, '(a)' ) '    ** GRID, by picking points from a grid;'
  write ( *, '(a)' ) '    ** HALTON, from a Halton sequence;'
  write ( *, '(a)' ) '    ** RANDOM, using FORTRAN RANDOM function;'
  write ( *, '(a)' ) '    ** UNIFORM, using a simple uniform RNG;'
  write ( *, '(a)' ) '    ** USER, call the "user" routine;'
  write ( *, '(a)' ) '  * IT_MAX, the maximum number of iterations.'
  write ( *, '(a)' ) '  * IT_FIXED, the number of iterative steps to take'
  write ( *, '(a)' ) '    using a fixed set of sampling points.'
  write ( *, '(a)' ) '  * SAMPLE, how to conduct the sampling.'
  write ( *, '(a)' ) '    ** GRID, by picking points from a grid;'
  write ( *, '(a)' ) '    ** HALTON, from a Halton sequence;'
  write ( *, '(a)' ) '    ** RANDOM, using FORTRAN RANDOM function;'
  write ( *, '(a)' ) '    ** UNIFORM, using a simple uniform RNG;'
  write ( *, '(a)' ) '    ** USER, call the "user" routine;'
  write ( *, '(a)' ) '  * SAMPLE_NUM, the number of sample points.'
  write ( *, '(a)' ) '  * BATCH, number of sample points to create at one time.'
  write ( *, '(a)' ) '  * OUTPUT, a file in which to store the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  To indicate that no further computations are '
  write ( *, '(a)' ) '  desired, it is enough to input a nonsensical value, '
  write ( *, '(a)' ) '  such as -1.'

  do

    write ( *, '(a)' ) '  *'
    write ( *, '(a)' ) ' *'
    write ( *, '(a)' ) '*  Ready to generate a new dataset:'
    write ( *, '(a)' ) ' *'
    write ( *, '(a)' ) '  *'
    write ( *, '(a)' ) '  Enter DIM_NUM, the spatial dimension:'
    write ( *, '(a)' ) '  (Try ''2'' if you do not have a preference.)'
    write ( *, '(a)' ) '  (0 or any negative value terminates execution).'

    read ( *, *, iostat = ios ) dim_num

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for DIM_NUM.'
      exit
    end if

    write ( *, '(a,i12)' ) '  User input DIM_NUM = ', dim_num

    if ( dim_num < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET'
      write ( *, '(a,i12)' ) '  The input value of DIM_NUM = ', dim_num
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter N, the number of points to generate:'
    write ( *, '(a)' ) '  (Try ''25'' if you do not have a preference.)'
    write ( *, '(a)' ) '  (0 or any negative value terminates execution).'

    read ( *, *, iostat = ios ) n

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for N.'
      exit
    end if

    write ( *, '(a,i12)' ) '  User input N = ', n

    if ( n < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET'
      write ( *, '(a,i12)' ) '  The input value of N = ', n
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter SEED, a seed for the random number generator:'
    write ( *, '(a)' ) '  (Try ''123456789'' if you do not have a preference.)'
    write ( *, '(a)' ) '  (Any negative value terminates execution).'

    read ( *, *, iostat = ios ) seed

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for SEED.'
      exit
    end if

    write ( *, '(a,i12)' ) '  User input SEED = ', seed

    if ( seed < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET'
      write ( *, '(a,i12)' ) '  The input value of SEED = ', seed
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

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
    write ( *, '(a)' ) '  (Try ''RANDOM'' if you do not have a preference.)'
    write ( *, '(a)' ) '  (A blank value terminates execution).'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter INIT:'
    read ( *, '(a)', iostat = ios ) init_string

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for INIT.'
      exit
    end if

    write ( *, '(a)' ) '  User input INIT = "' // trim ( init_string ) // '".'

    if ( s_eqi ( init_string, 'RANDOM'  ) ) then
      init = -1
    else if ( s_eqi ( init_string, 'UNIFORM' ) ) then
      init = 0
    else if ( s_eqi ( init_string, 'HALTON'  ) ) then
      init = 1
    else if ( s_eqi ( init_string, 'GRID'    ) ) then
      init = 2
    else if ( s_eqi ( init_string, 'USER'    ) ) then
      init = 3
    else if ( 0 < len_trim ( init_string ) ) then
      init = 4
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET'
      write ( *, '(a)' ) '  The input value of INIT '
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

    if ( len_trim ( init_string ) <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET'
      write ( *, '(a)' ) '  The input value of INIT '
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

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
    write ( *, '(a)' ) '  (Try ''50'' if you do not have a preference.)'
    write ( *, '(a)' ) '  (A negative value terminates execution).'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter IT_MAX:'
    read ( *, *, iostat = ios ) it_max

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for IT_MAX.'
      exit
    end if

    write ( *, '(a,i6)' ) '  User input IT_MAX = ', it_max

    if ( it_max < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET'
      write ( *, '(a,i12)' ) '  The input value of IT_MAX = ', it_max
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
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
    write ( *, '(a,i6,a)' ) '  (Try ', it_max, &
      ' if you do not have a preference.)'
    write ( *, '(a)' ) '  (A 0 or negative value terminates execution).'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter IT_FIXED:'
    read ( *, *, iostat = ios ) it_fixed

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for IT_FIXED.'
      exit
    end if

    write ( *, '(a,i6)' ) '  User input IT_FIXED = ', it_fixed

    if ( it_fixed <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET'
      write ( *, '(a,i12)' ) '  The input value of IT_FIXED = ', it_fixed
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
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
    write ( *, '(a)' ) '  (Try ''RANDOM'' if you do not have a preference.)'
    write ( *, '(a)' ) '  (A blank value terminates execution).'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter SAMPLE:'
    read ( *, '(a)', iostat = ios ) sample_string

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for SAMPLE.'
      exit
    end if

    write ( *, '(a)' ) '  User input SAMPLE = "' // trim ( sample_string ) &
      // '".'

    if ( len_trim ( sample_string ) <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET'
      write ( *, '(a)' ) '  The input value of SAMPLE '
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

    if ( s_eqi ( sample_string, 'RANDOM'  ) ) then
      sample = -1
    else if ( s_eqi ( sample_string, 'UNIFORM' ) ) then
      sample = 0
    else if ( s_eqi ( sample_string, 'HALTON'  ) ) then
      sample = 1
    else if ( s_eqi ( sample_string, 'GRID'    ) ) then
      sample = 2
    else if ( s_eqi ( sample_string, 'USER' ) ) then
      sample = 3
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET'
      write ( *, '(a)' ) '  The input value of SAMPLE '
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
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
    write ( *, '(a)' ) '  number of generators.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (Try ''10000'' if you do not have a preference.)'
    write ( *, '(a)' ) '  (A zero or negative value terminates execution.)'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter SAMPLE_NUM:'

    read ( *, *, iostat = ios ) sample_num

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for SAMPLE_NUM.'
      exit
    end if

    write ( *, '(a,i12)' ) '  User input SAMPLE_NUM = ', sample_num

    if ( sample_num <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET'
      write ( *, '(a,i12)' ) '  The input value of SAMPLE_NUM = ', sample_num
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
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
      write ( *, '(a)' ) 'CVT_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for SAMPLE_NUM.'
      exit
    end if

    write ( *, '(a,i12)' ) '  User input BATCH = ', batch

    if ( batch <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET'
      write ( *, '(a,i12)' ) '  The input value of BATCH = ', batch
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  OUTPUT is the name of a file into which'
    write ( *, '(a)' ) '  the computed data may be stored.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (Try ''cvt.txt'' if you do not have a preference.)'
    write ( *, '(a)' ) '  (A blank value terminates execution).'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter OUTPUT:'
    read ( *, '(a)', iostat = ios ) file_out_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for OUTPUT.'
      exit
    end if

    write ( *, '(a)' ) '  User input OUTPUT = "' // trim ( file_out_name ) &
      // '".'

    if ( len_trim ( file_out_name ) <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_DATASET'
      write ( *, '(a)' ) '  The input value of OUTPUT '
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if
!
!  Initialize the data.
!
    allocate ( r(1:dim_num,1:n) )

    if ( init == 4 ) then

      call data_read ( init_string, dim_num, n, r, success )

      if ( .not. success ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CVT_DATASET - Fatal error!'
        write ( *, '(a)' ) '  The data could not be read from the file.'
        stop
      end if

    end if

    seed_init = seed

    call cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed, &
      seed, r, it_num, it_diff, energy )
!
!  Write the data to a file.
!
    call cvt_write ( dim_num, n, batch, seed_init, seed, init_string, it_max, &
      it_fixed, it_num, it_diff, energy, sample_string, sample_num, r, &
      file_out_name, comment )

    deallocate ( r )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The data was written to the file "' &
       // trim ( file_out_name ) // '".'

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Final value of SEED = ', seed

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
