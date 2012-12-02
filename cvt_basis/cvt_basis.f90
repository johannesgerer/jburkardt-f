program main

!*****************************************************************************80
!
!! MAIN is the main routine for the CVT_BASIS program.
!
!  Discussion:
!
!    What we really want to do is take a thousand points and put
!    them into 8 or 9 clusters.  Each cluster will be represented
!    by a special generator point (which need not be one of the original
!    data points).  The cost of a particular clustering is the
!    sum of the squares of the distances of each data point to its
!    generator point.
!
!    The goal of this particular application is then to use the
!    generator points as a good basis for representing the dynamics
!    of the system that generated the original data.
!
!    Or, in other words, we're going to solve a PDE, save a bunch
!    of solutions, pick out a few representative ones, and assume
!    that the most interesting behavior of solutions to the PDE
!    is incorporated in the representatives.
!
!
!    The data to be examined is assumed to be stored in a file.
!
!    The file is assumed to contain a number of records, with each
!    record stored on its own line.
!
!    Each record, in turn, contains a fixed number of data values
!    that describe a particular gene expression experiment.
!
!    Each record will be regarded as a point in N dimensional space.
!
!    The program will try to cluster the data, that is, to organize
!    the data by defining a number of cluster centers, which are
!    also points in N dimensional space, and assigning each record
!    to the cluster associated with a particular center.
!
!    The method of assigning data aims to minimize the cluster energy,
!    which is taken to be the sum of the squares of the distances of
!    each data point from its cluster center.
!
!    In some contexts, it makes sense to use the usual Euclidean sort
!    of distance.  In others, it may make more sense to replace each
!    data record by a normalized version, and to assign distance
!    by computing angles between the unit vectors.
!
!  Diary:
!
!    20 July 2005: starting from CVT_BASIS_FLOW, make a "flow free"
!    version, that is, one which drops all the flow features.
!    This will stand in relationship to CVT_BASIS_FLOW as 
!    SVD_BASIS is related to POD_BASIS_FLOW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: v0_file_max = 10

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) a_dot_p
  integer ( kind = 4 ), allocatable, dimension ( : ) :: cluster
  real ( kind = 8 ), allocatable, dimension (:,:) :: cluster_center
  integer ( kind = 4 ), parameter :: cluster_max = 10
  integer ( kind = 4 ) cluster_hi
  integer ( kind = 4 ) cluster_it_max
  integer ( kind = 4 ) cluster_lo
  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) cluster_range
  real ( kind = 8 ), allocatable, dimension ( : ) :: cluster_size
  real ( kind = 8 ), allocatable, dimension ( : ) :: cluster_size_inv
  logical comment
  character comment_char
  integer ( kind = 4 ) comp_num
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) r8vec_norm2
  character ( len = 255 ) element_file_name
  integer ( kind = 4 ) element_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: energy
  integer ( kind = 4 ) energy_it_max
  logical file_exist
  character ( len = 255 ) gen_file_name
  integer ( kind = 4 ) gen_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: node
  integer ( kind = 4 ) node_num
  real ( kind = 8 ) norm
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ), parameter :: null_cluster_policy = 1
  real ( kind = 8 ), allocatable, dimension ( :, :) :: point
  integer ( kind = 4 ) point_num
  logical s_eqi
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) swap_num
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) time_cvt
  real ( kind = 8 ), allocatable, dimension ( : ) :: v
  character ( len = 255 ) v_file_name
  integer ( kind = 4 ) v_file_num
  character ( len = 255 ) v0_file_name(v0_file_max)
  integer ( kind = 4 ) v0_file_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: v1
  logical, parameter :: verbose = .false.

  call timestamp ( )

  time_cvt = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_BASIS'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Arrange a set abstract data vectors into clusters.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This version of the code handles up to'
  write ( *, '(i6,a)' ) v0_file_max, ' families of data.'

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Null cluster policy:'
    write ( *, '(a)' ) '  0, do nothing, accept null clusters;'
    write ( *, '(a)' ) '  1, reset center to a random data point;'
    write ( *, '(a)' ) '  2, reset center to random point in hull;'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  NULL_CLUSTER_POLICY = ', null_cluster_policy
  end if

  seed = 123456789

  call random_initialize ( seed )

  comp_num = 0

  v_file_num = 0
  v0_file_num = 0
  i = 0
!
!  Get the V0 file name.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  It is time to read sets of solution files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We assume that a series of solution files exists,'
  write ( *, '(a)' ) '  with "consecutive" names, like'
  write ( *, '(a)' ) '    fred001.txt, fred002,txt, ...'
  write ( *, '(a)' ) '  Just specify the FIRST name in the series, and'
  write ( *, '(a)' ) '  the program will read them all.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program can read another series of data if'
  write ( *, '(a)' ) '  you specify yet another first name, or you can'
  write ( *, '(a)' ) '  type "none" when there are no more file series.'
  write ( *, '(a)' ) ' '

  do while ( v0_file_num < v0_file_max )

    i = i + 1

    if ( i == 1 ) then
      call s_input ( 'What is the first solution file (in the first series)?', &
        v0_file_name(i), ierror )
    else
      call s_input ( &
        'What is the first solution file (in the NEXT series) or "NONE"?', &
        v0_file_name(i), ierror )
    end if

    if ( s_eqi ( v0_file_name(i), 'NONE' ) ) then
      exit
    end if

    v0_file_num = v0_file_num + 1
!
!  Presumably, all the solution files have the same name as the first
!  solution file, but with a numerical increment.  To begin with, simply count
!  the number of files.
!
    v_file_name = v0_file_name(i)

    do

      if ( .not. file_exist ( v_file_name ) ) then
        exit
      end if

      v_file_num = v_file_num + 1

      call file_name_inc ( v_file_name )

    end do

  end do

  if ( v_file_num == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_BASIS - Fatal error!'
    write ( *, '(a)' ) '  There do not seem to be any solution files.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) &
    '  The number of initial solution files is ', v0_file_num
  write ( *, '(a,i6)' ) &
    '  The total number of solution files is ', v_file_num
!
!  If there was no steady state file, then you need to read the
!  first data file in order to determine the number of components.
!
  if ( comp_num == 0 ) then

    call r8mat_header_read ( v0_file_name(1), comp_num, node_num )

    allocate ( v(comp_num*node_num) )

  end if
!
!  Now we have enough information to set up a data structure.
!
!  Determine the spatial dimension (columns) and number of points (rows).
!
  dim_num = comp_num * node_num

  point_num = v_file_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data is stored in an M by N matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The "spatial" dimension M is   ', dim_num
  write ( *, '(a,i6)' ) '  The number of data points N is ', point_num
!
!  Allocate space for the POINT array.
!
  allocate ( point(1:dim_num,1:point_num) )
!
!  Now read the data from the individual files, process it if necessary,
!  and gather it into a single array called POINT.
!
  write ( *, '(a)' ) ' '
  j = 0

  do i = 1, v0_file_num

    v_file_name = v0_file_name(i)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Processing files starting with "' &
      // trim ( v0_file_name(i) ) // '".'

    do

      if ( .not. file_exist ( v_file_name ) ) then
        exit
      end if

      j = j + 1

      if ( .false. ) then
        write ( *, '(2x,i6,2x,a)' ) j, trim ( v_file_name )
      end if

      call r8mat_data_read ( v_file_name, comp_num, node_num, v )

      point(1:dim_num,j) = v(1:dim_num)

      call file_name_inc ( v_file_name )

    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  All the data has been read.'
!
!  Get the range of cluster sizes to check.
!
  call i4_range_input ( 'Enter lower and upper number of clusters', &
    cluster_lo, cluster_hi, ierror )

  if ( point_num < cluster_hi ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_BASIS - Warning!'
    write ( *, '(a)' ) '  CLUSTER_HI exceeds number of points.'
    write ( *, '(a,i12)' ) '  CLUSTER_HI = ', cluster_hi
    write ( *, '(a,i12)' ) '  POINT_NUM =  ', point_num
    cluster_hi = point_num
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program will try to determine the minimum energy'
  write ( *, '(a)' ) '  of a clustering, for cluster sizes in the range:'
  write ( *, '(2x,2i6)' ) cluster_lo, cluster_hi

  allocate ( cluster(1:point_num) )
  allocate ( cluster_center(1:dim_num,1:cluster_hi) )
  allocate ( energy(1:cluster_hi) )
!
!  Get the number of different random starting cluster configurations.
!
  call i4_input ( &
    'Enter the number of different random cluster configurations to check', &
    cluster_it_max, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For each number of clusters, the number of'
  write ( *, '(a)' ) '  distinct initial random configurations to be checked'
  write ( *, '(a,i6)' ) '  will be  ', cluster_it_max
!
!  Get the number of energy iterations for a particular random start:
!
  call i4_input ( 'Enter the number of energy iterations', &
    energy_it_max, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For each initial random configuration, the number of'
  write ( *, '(a)' ) '  times the program will recompute the cluster centers,'
  write ( *, '(a,i6)' ) '  cluster components, and energy is ', energy_it_max

  call cpu_time ( t1 )

  call analysis_raw ( dim_num, point_num, point, cluster_lo, &
    cluster_hi, cluster_it_max, energy_it_max, energy, cluster_center, &
    null_cluster_policy, seed )

  call cpu_time ( t2 )
  time_cvt = time_cvt + t2 - t1
!
!  Multiple cluster values, presumably we want an energy plot.
!
  if ( cluster_lo < cluster_hi ) then

    allocate ( cluster_size(1:cluster_hi) )
    allocate ( cluster_size_inv(1:cluster_hi) )

    do i = 1, cluster_hi
      cluster_size(i) = real ( i, kind = 8 )
    end do

    cluster_size_inv(1:cluster_hi) = 1.0D+00 / cluster_size(1:cluster_hi)
      cluster_range = cluster_hi + 1 - cluster_lo

    call data_to_gnuplot ( 'raw.txt', cluster_range, &
      cluster_size(cluster_lo:cluster_hi), energy(cluster_lo:cluster_hi) )

    call data_to_gnuplot ( 'raw2.txt', cluster_range, &
      cluster_size_inv(cluster_lo:cluster_hi), energy(cluster_lo:cluster_hi) )

    deallocate ( cluster_size )
    deallocate ( cluster_size_inv )
!
!  If there was only one cluster value, presumably we want to
!  * print out the population of each cluster;
!  * write out the generators.
!
  else
!
!  Why am I recomputing the cluster assignments here?
!
    cluster_num = cluster_hi

    cluster(1:point_num) = 1

    call nearest_cluster_raw ( dim_num, point_num, cluster_num, &
      cluster_center, point, cluster, swap_num )

    call cluster_census ( dim_num, point_num, cluster_num, cluster_center, &
      point, cluster )

    if ( .false. ) then
      call cluster_list ( point_num, cluster )
    end if
!
!  Write vectors to files.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_BASIS:'
    write ( *, '(a)' ) '  Ready to write the cluster generators to files.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Do you want comments in the header of the file?'
    write ( *, '(a)' ) '  (These begin with the "#" character.) (Y/N)'

    call s_input ( '  Enter Y or N:', comment_char, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_BASIS_FLOW - Fatal error!'
      write ( *, '(a)' ) '  Input error reading the comment option.'
      write ( *, '(a)' ) '  We will assume comments are acceptable.'
      comment_char = 'Y'
    end if

    if ( comment_char == 'Y' .or. comment_char == 'y' ) then
      comment = .true.
      write ( *, * ) 'The output files will include a commented header.'
    else
      comment = .false.
      write ( *, * ) 'The output files will NOT include a commented header.'
    end if

    gen_file_name = 'gen_000.txt'
    call get_unit ( gen_unit )

    do j = 1, cluster_hi

      call file_name_inc ( gen_file_name )

      if ( j == 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Writing first file ' // trim ( gen_file_name )
      end if

      if ( j == cluster_hi ) then
        write ( *, '(a)' ) '  Writing last file  ' // trim ( gen_file_name )
      end if

      v(1:comp_num*node_num) = cluster_center(1:comp_num*node_num,j)
      
      call r8mat_write ( gen_file_name, comp_num, node_num, v )

    end do

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CPU_TIME (seconds):'
  write ( *, '(a,g14.6)' ) '  CVT:         ', time_cvt
!
!  Free memory.
!
  deallocate ( cluster )
  deallocate ( cluster_center )
  deallocate ( energy )
  deallocate ( point )
  deallocate ( v )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_BASIS'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine analysis_raw ( dim_num, point_num, point, cluster_lo, cluster_hi, &
  cluster_it_max, energy_it_max, energy, cluster_center, null_cluster_policy, &
  seed )

!*****************************************************************************80
!
!! ANALYSIS_RAW computes the energy for a range of number of clusters.
!
!  Discussion:
!
!    This version of the analysis routine is for "raw" data.  That is,
!    all the solution vectors are treated as points in Euclidean space
!    with the usual distance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NORMAL:
!    0, analyze the raw data.
!    1, analyze the normalized data.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, integer ( kind = 4 ) CLUSTER_LO, CLUSTER_HI, the low and high
!    cluster numbers that define a range of clusters to check.
!
!    Input, integer ( kind = 4 ) CLUSTER_IT_MAX, the number of different random
!    startup configurations to try.
!
!    Input, integer ( kind = 4 ) ENERGY_IT_MAX, the maximum number of energy iterations.
!
!    Output, real ( kind = 8 ) ENERGY(CLUSTER_HI), contains in entries
!    CLUSTER_LO through CLUSTER_HI the estimated minimum energy over 
!    all clusterings of this size.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_HI), contains the
!    generators for the minimum cluster energy calculation, but only
!    for the last cluster size considered, namely, CLUSTER_HI.
!
!    Input, integer ( kind = 4 ) NULL_CLUSTER_POLICY, specifies what to do if a
!    null cluster is encountered.
!    0, do nothing.
!    1, reset center of null cluster to a random data point;
!    2, reset center of null cluster to a random point in the data hull;
!
!    Input/output, integer SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) cluster_hi
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension (dim_num,cluster_hi) :: cluster_center
  real ( kind = 8 ), dimension (dim_num,cluster_hi) :: cluster_center_test
  integer ( kind = 4 ) cluster_it
  integer ( kind = 4 ) cluster_it_max
  integer ( kind = 4 ) cluster_lo
  integer ( kind = 4 ) cluster_num
  logical, parameter :: debug = .true.
  real ( kind = 8 ), dimension (cluster_hi) :: energy
  integer ( kind = 4 ) energy_it_max
  real ( kind = 8 ) energy_max(cluster_hi)
  real ( kind = 8 ) energy_min(cluster_hi)
  real ( kind = 8 ) energy_test
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) null_cluster_policy
  real ( kind = 8 ), dimension (dim_num,point_num) :: point
  real ( kind = 8 ), dimension (dim_num) :: r_min
  real ( kind = 8 ), dimension (dim_num) :: r_max
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
!
!  Compute the minimum and maximum component values.
!
  do i = 1, dim_num
    r_min(i) = minval ( point(i,1:point_num) )
    r_max(i) = maxval ( point(i,1:point_num) )
  end do

  if ( dim_num <= 20 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Minimum and maximum data values:'
    write ( *, '(a)' ) ' '
    do i = 1, dim_num
      write ( *, '(i6,2f10.4)' ) i, r_min(i), r_max(i)
    end do
  end if
!
!  Consider a range of clusters.
!
  do cluster_num = cluster_lo, cluster_hi

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'Number of clusters allowed: ', cluster_num
!
!  For each cluster size, try several random starting configurations.
!
    energy_min(cluster_num) = huge ( energy_min )
    energy_max(cluster_num) = 0.0D+00

    do cluster_it = 1, cluster_it_max

      write ( *, '(a)' ) ' '
      write ( *, '(i6)' ) cluster_it

      call cluster_initialize_raw ( dim_num, point_num, cluster_num, point, &
        cluster, cluster_center_test, energy_test )

      it = 0
      write ( *, '(a,g14.6,i6)' ) 'Initial_RAW  ', energy_test, it

      call hmeans_raw ( dim_num, point_num, cluster_num, energy_it_max, &
        it, point, cluster, cluster_center_test, energy_test, &
        null_cluster_policy, seed )

      write ( *, '(a,g14.6,i6)' ) 'HMEANS_RAW   ', energy_test, it

      call kmeans_raw ( dim_num, point_num, cluster_num, energy_it_max, &
        it, point, cluster, cluster_center_test, energy_test, &
        null_cluster_policy, seed )

      write ( *, '(a,g14.6,i6)' ) 'KMEANS_RAW   ', energy_test, it

      if ( energy_test < energy_min(cluster_num) ) then
        energy_min(cluster_num) = energy_test
        cluster_center(1:dim_num,1:cluster_hi) = &
          cluster_center_test(1:dim_num,1:cluster_hi)
      end if

      energy_max(cluster_num) = max ( energy_max(cluster_num), energy_test )

    end do

    energy(cluster_num) = energy_min(cluster_num)

  end do
!
!  Report energy ranges for the various cluster sizes.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ANALYSIS_RAW:'
  write ( *, '(a)' ) '  Computed energy range for given cluster size:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (The minimum and maximum should be close if'
  write ( *, '(a)' ) '  we''re taking enough iterations.)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Number'
  write ( *, '(a)' ) '  of       Minimum      Maximum'
  write ( *, '(a)' ) '  Clusters Energy       Energy'
  write ( *, '(a)' ) ' '

  do cluster_num = cluster_lo, cluster_hi
    write ( *, '(i7,2f14.4)' ) &
      cluster_num, energy_min(cluster_num), energy_max(cluster_num)
  end do
!
!  Report best energy for the various cluster sizes.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Energy table:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Number'
  write ( *, '(a)' ) 'of                   Energy'
  write ( *, '(a)' ) 'Clusters  Energy     /point      Sqrt(E/Pt)'
  write ( *, '(a)' ) ' '

  do cluster_num = cluster_lo, cluster_hi
    t = energy(cluster_num) / dble ( point_num )
    write ( *, '(i7,3f14.4)' ) cluster_num, energy(cluster_num), t, sqrt ( t )
  end do

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
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

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

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
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!
  implicit none

  character c
  logical ch_is_digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
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
subroutine cluster_census ( dim_num, point_num, cluster_num, cluster_center, &
  point, cluster )

!*****************************************************************************80
!
!! CLUSTER_CENSUS computes and prints the population of each cluster.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the cluster
!    generators.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster to which each
!    point is assigned.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) cluster(point_num)
  real ( kind = 8 ), dimension (dim_num,cluster_num) :: cluster_center
  real ( kind = 8 ), dimension ( cluster_num ) :: cluster_energy
  integer ( kind = 4 ), dimension ( cluster_num ) :: cluster_max
  integer ( kind = 4 ), dimension ( cluster_num ) :: cluster_min
  integer ( kind = 4 ), dimension ( cluster_num ) :: cluster_population
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) percent1
  integer ( kind = 4 ) percent2
  real ( kind = 8 ), dimension (dim_num,point_num) :: point
  integer ( kind = 4 ) swap_num
  real ( kind = 8 ) total_energy

  call nearest_cluster_raw ( dim_num, point_num, cluster_num, &
    cluster_center, point, cluster, swap_num )

  cluster_energy(1:cluster_num) = 0.0D+00
  cluster_population(1:cluster_num) = 0

  do i = 1, point_num

    j = cluster(i)

    cluster_population(j) = cluster_population(j) + 1

    cluster_energy(j) = cluster_energy(j) &
      + sum ( ( cluster_center(1:dim_num,j) - point(1:dim_num,i) )**2 )

  end do

  total_energy = sum ( cluster_energy(1:cluster_num) )

  cluster_min(1:cluster_num) = point_num + 1
  cluster_max(1:cluster_num) = 0

  do i = 1, point_num
    j = cluster(i)
    cluster_min(j) = min ( cluster_min(j), i )
    cluster_max(j) = max ( cluster_max(j), i )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CLUSTER_CENSUS'
  write ( *, '(a)' ) '  Individual cluster population and energy'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Index    Population   Percentage   Energy  Percentage  Min  Max'
  write ( *, '(a)' ) ' '
  do i = 1, cluster_num
    percent1 = int ( dble ( 100 * cluster_population(i) ) / dble ( point_num ) )
    percent2 = int ( 100.0D+00 * cluster_energy(i) ) / total_energy
    write ( *, '(1x,i6,8x,i6,10x,i3,g14.6,4x,i3,2i5)' ) &
      i, cluster_population(i), &
      percent1, cluster_energy(i), percent2, cluster_min(i), cluster_max(i)
  end do

  percent1 = 100
  percent2 = 100
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               ------          ---  ------------    ---'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a5,8x,i6,10x,i3,g14.6,4x,i3,2i5)' ) 'Total', &
  sum ( cluster_population(1:cluster_num) ), &
    percent1, sum ( cluster_energy(1:cluster_num) ), percent2, 1, point_num

  return
end
subroutine cluster_initialize_raw ( dim_num, point_num, cluster_num, point, &
  cluster, cluster_center, energy )

!*****************************************************************************80
!
!! CLUSTER_INITIALIZE_RAW initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, each cluster center is a random convex combination
!    of the data points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates of 
!    the points.
!
!    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the clusters to which points
!    are assigned.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), 
!    the coordinates of the cluster centers.
!
!    Output, real ( kind = 8 ) ENERGY, the energy of the clustering.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) cluster(point_num)
  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  real ( kind = 8 ) column_sum
  real ( kind = 8 ) energy
  real ( kind = 8 ) factor(point_num,cluster_num)
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(dim_num,point_num)
  integer ( kind = 4 ) swap_num
!
!  Get a PxC block of random factors.
!
  call random_number ( harvest = factor(1:point_num,1:cluster_num) )
!
!  Make each column of factors have unit sum.
!
  do j = 1, cluster_num
    column_sum = sum ( factor(1:point_num,j) )
    factor(1:point_num,j) = factor(1:point_num,j) / column_sum
  end do
!
!  Set centers = points * factors.
!
  cluster_center(1:dim_num,1:cluster_num) = &
    matmul ( point(1:dim_num,1:point_num), factor(1:point_num,1:cluster_num) )
!
!  Assign points to the nearest centers.
!
  call nearest_cluster_raw ( dim_num, point_num, cluster_num, &
    cluster_center, point, cluster, swap_num )
!
!  Determine the energy of the clustering.
!
  call energy_raw ( dim_num, point_num, cluster_num, point, &
    cluster_center, cluster, energy )

  return
end
subroutine cluster_list ( point_num, cluster )

!*****************************************************************************80
!
!! CLUSTER_LIST prints out the assignments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster to which each
!    point is assigned.
!
  implicit none

  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) cluster(point_num)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I  Cluster(I)'
  write ( *, '(a)' ) ' '
  do i = 1, point_num
    write ( *, '(i4,8x,i4)' ) i, cluster(i)
  end do

  return
end
subroutine data_to_gnuplot ( file_name, n, x, y )

!*****************************************************************************80
!
!! DATA_TO_GNUPLOT writes data to a file suitable for processing by GNUPLOT.
!
!  Discussion:
!
!    Once the data file is written, the GNUPLOT program can be used
!    to display the data, using commands like:
!
!      set term post default
!      set grid
!      set yrnage [0:*]
!      set title "Number of Clusters vs Total Energy, Normalized Data"
!      plot 'file_name'
!      quit
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file
!    to which the data is written.
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the data.
!
  integer ( kind = 4 ) n

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace' )

  do i = 1, n
    write ( iunit, '(g16.8,2x,g16.8)' ) x(i), y(i)
  end do

  close ( unit = iunit )

  return
end
subroutine digit_inc ( c )

!*****************************************************************************80
!
!! DIGIT_INC increments a decimal digit.
!
!  Example:
!
!    Input  Output
!    -----  ------
!    '0'    '1'
!    '1'    '2'
!    ...
!    '8'    '9'
!    '9'    '0'
!    'A'    'A'
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
!    Input/output, character C, a digit to be incremented.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  call ch_to_digit ( c, digit )

  if ( digit == -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, c )

  return
end
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
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
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine energy_raw ( dim_num, point_num, cluster_num, point, &
  cluster_center, cluster, energy )

!*****************************************************************************80
!
!! ENERGY_RAW computes the total energy of a given clustering.
!
!  Discussion:
!
!    This routine is used with the raw data.  No normalization is
!    done to the data. 
!
!    The total energy function is the sum of the cluster energies.
!
!    The energy of a cluster is the sum of the squares of the distances of
!    each point in the cluster to the center point of the cluster.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the 
!    center points.
!
!    Input, integer ( kind = 4 ) CLUSTER(POINT_NUM), indicates the cluster 
!    to which each data point belongs.
!
!    Output, real ( kind = 8 ) ENERGY, the total energy of the clustering.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  real ( kind = 8 ), dimension ( cluster_num ) :: cluster_energy
  real ( kind = 8 ) energy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point

  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_energy(j) = cluster_energy(j) + sum ( &
      ( cluster_center(1:dim_num,j) - point(1:dim_num,i) )**2 )
  end do

  energy = sum ( cluster_energy(1:cluster_num) )

  return
end
subroutine file_column_count ( input_filename, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some 
!    comment lines,
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
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in 
!    the file.
!
  implicit none

  integer ( kind = 4 ) column_num
  logical got_one
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( input_filename )
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = ios ) line

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

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = ios ) line

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

  close ( unit = input_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

  return
end
function file_exist ( file_name )

!*****************************************************************************80
!
!! FILE_EXIST reports whether a file exists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, logical FILE_EXIST, is TRUE if the file exists.
!
  implicit none

  character ( len = * ) file_name
  logical file_exist

  inquire ( file = file_name, exist = file_exist )

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC generates the next filename in a series.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Example:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  logical ch_is_digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( ch_is_digit ( c ) ) then

      call digit_inc ( c )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

  return
end
subroutine file_row_count ( input_filename, row_num )

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
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = ios ) line

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

  close ( unit = input_unit )

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
subroutine hmeans_raw ( dim_num, point_num, cluster_num, it_max, it, &
  point, cluster, cluster_center, energy, null_cluster_policy, seed )

!*****************************************************************************80
!
!! HMEANS_RAW seeks the minimal energy of a cluster of a given size.
!
!  Discussion:
!
!    This routine works with the raw data, and does not do any
!    normalization.
!
!    The data for the H-Means problem is a set of N points X in
!    M-dimensions, and a desired number of clusters K.
!
!    The goal is to determine K points Z, called cluster centers, so that
!    if we associate each point X with its nearest Z value, we minimize
!    the standard deviation or cluster energy.  Writing CLUSTER(I) to
!    indicate the index of the nearest cluster center to point X(I), the
!    energy can be written as:
!
!      Energy = Sum ( 1 <= I <= N ) distance ( X(I), Z(CLUSTER(I)) )**2
!
!    where
!
!      distance ( X - Z ) = Sqrt ( Sum ( 1 <= J <= M ) ( X(J) - Z(J) )**2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez and Angel Martinez,
!    Computational Statistics Handbook with MATLAB,
!    pages 373-376,
!    Chapman and Hall / CRC, 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations allowed.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster to which each
!    data point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), 
!    the centers associated with the minimal energy clustering.
!
!    Output, real ( kind = 8 ) ENERGY, the total energy associated with 
!    the minimal energy clustering.
!
!    Input, integer ( kind = 4 ) NULL_CLUSTER_POLICY, specifies what to do if a
!    null cluster is encountered.
!    0, do nothing.
!    1, reset center of null cluster to a random data point;
!    2, reset center of null cluster to a random point in the data hull;
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  integer ( kind = 4 ), dimension ( cluster_num ) :: cluster_population
  real ( kind = 8 ) column_sum
  logical, parameter :: debug = .false.
  real ( kind = 8 ) energy
  real ( kind = 8 ) factor(point_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) null_cluster_num
  integer ( kind = 4 ) null_cluster_policy
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) swap_num

  do it = 1, it_max
!
!  #1: Assign each point to the cluster of its nearest center.
!
    call nearest_cluster_raw ( dim_num, point_num, cluster_num, &
      cluster_center, point, cluster, swap_num )
!
!  #2: Determine the energy of the new clustering with the current centers.
!
    call energy_raw ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, energy )

    if ( .false. ) then
      write ( *, '(i6,3x,14x,3x,g14.6)' ) it, energy
    end if
!
!  #3: Determine the populations of the new clusters.
!  If null clusters are not OK, we need to handle any empty clusters.
!  If any cluster is empty, set its center to a random point within
!  the convex hull of the data, reassign points, and repeat the entire
!  process.
!
    do

      cluster_population(1:cluster_num) = 0

      do i = 1, point_num
        j = cluster(i)
        cluster_population(j) = cluster_population(j) + 1
      end do

      null_cluster_num = 0

      do j = 1, cluster_num
        if ( cluster_population(j) == 0 ) then
          null_cluster_num = null_cluster_num + 1
        end if
      end do

      if ( null_cluster_num == 0 ) then
        exit
      end if

      if ( null_cluster_policy == 0 ) then
        exit
      end if

      if ( debug ) then
        write ( *, '(a,i6)' ) &
          'HMEANS_RAW, number of empty clusters = ', null_cluster_num
      end if

      if ( null_cluster_policy == 1 ) then

        do j = 1, cluster_num
          if ( cluster_population(j) == 0 ) then
            k = i4_uniform ( 1, point_num, seed )
            cluster_center(1:dim_num,j) = point(1:dim_num,k)
          end if
        end do

      else if ( null_cluster_policy == 2 ) then

        do j = 1, cluster_num
          if ( cluster_population(j) == 0 ) then
            call random_number ( harvest = factor(1:point_num) )
            column_sum = sum ( factor(1:point_num) )
            factor(1:point_num) = factor(1:point_num) / column_sum
            cluster_center(1:dim_num,j) = &
              matmul ( point(1:dim_num,1:point_num), factor(1:point_num) )
          end if
        end do

      end if
!
!  Resort the points.
!
      call nearest_cluster_raw ( dim_num, point_num, cluster_num, &
        cluster_center, point, cluster, swap_num )

    end do
!
!  #4: Recompute the cluster centers as the centroids of the points
!  in the cluster.
!
    cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

    do i = 1, point_num
      j = cluster(i)
      cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
        + point(1:dim_num,i)
    end do

    do i = 1, cluster_num
      cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
        / dble ( max ( cluster_population(i), 1 ) )
    end do
!
!  #5: Determine the energy of the current clustering with the new centers.
!
    call energy_raw ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, energy )

    if ( .false. ) then
      write ( *, '(i6,3x,g14.6)' ) it, energy
    end if

    if ( swap_num == 0 .and. 1 < it ) then
      exit
    end if

  end do

  return
end
subroutine i4_input ( string, value, ierror )

!*****************************************************************************80
!
!! I4_INPUT prints a prompt string and reads an integer from the user.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#') or is
!    blank, the routine ignores that line, and tries to read the next one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, integer ( kind = 4 ) VALUE, the value input by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  character ( len = 255 ) line
  character ( len = * ) string
  integer ( kind = 4 ) value

  ierror = 0
  value = huge ( value )
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) line

    if ( ierror /= 0 ) then
      return
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if
!
!  Extract integer information from the string.
!
    call s_to_i4 ( line, value, ierror, last )

    if ( ierror /= 0 ) then
      value = huge ( value )
      return
    end if

    exit

  end do

  return
end
subroutine i4_range_input ( string, value1, value2, ierror )

!*****************************************************************************80
!
!! I4_RANGE_INPUT reads a pair of integers from the user, representing a range.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#') or is blank,
!    the routine ignores that line, and tries to read the next one.
!
!    The pair of integers may be separated by spaces or a comma or both.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, integer ( kind = 4 ) VALUE1, VALUE2, the values entered by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  character, parameter :: comma = ','
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  integer ( kind = 4 ) last2
  character ( len = 255 ) line
  character, parameter :: space = ' '
  character ( len = * ) string
  integer ( kind = 4 ) value1
  integer ( kind = 4 ) value2

  ierror = 0
  value1 = huge ( value1 )
  value2 = huge ( value2 )
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) line

    if ( ierror /= 0 ) then
      return
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if
!
!  Remove commas.
!
    call s_rep_ch ( line, comma, space )
!
!  Extract integer information from the string.
!
    call s_to_i4 ( line, value1, ierror, last )

    if ( ierror /= 0 ) then
      value1 = huge ( value1 )
      return
    end if

    call s_to_i4 ( line(last+1:), value2, ierror, last2 )

    if ( ierror /= 0 ) then
      value2 = huge ( value2 )
      return
    end if

    exit

  end do

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2006
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
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an integer vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2000
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
  integer ( kind = 4 ) big
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(i6,1x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(i6,1x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i6,i11)' ) i, a(i)
    end do
  end if

  return
end
subroutine kmeans_raw ( dim_num, point_num, cluster_num, it_max, it, &
  point, cluster, cluster_center, energy, null_cluster_policy, seed )

!*****************************************************************************80
!
!! KMEANS_RAW tries to improve a partition of points.
!
!  Discussion:
!
!    This routine works with the raw data, and does not do any
!    normalization.
!
!    The data for the K-Means problem is a set of N points X in
!    M-dimensions, and a desired number of clusters K.
!
!    The goal is to determine K points Z, called cluster centers, so that
!    if we associate each point X with its nearest Z value, we minimize
!    the standard deviation or cluster energy.  Writing CLUSTER(I) to
!    indicate the index of the nearest cluster center to point X(I), the
!    energy can be written as:
!
!      Energy = Sum ( 1 <= I <= N ) distance ( X(I), Z(CLUSTER(I)) )^2
!
!    where
!
!      distance ( X - Z ) = Sqrt ( Sum ( 1 <= J <= M ) ( X(J) - Z(J) )^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez and Angel Martinez,
!    Computational Statistics Handbook with MATLAB,
!    pages 373-376,
!    Chapman and Hall / CRC, 2002
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations allowed.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster to which
!    each point belongs.  On output, these may have been altered.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), 
!    the centers associated with the clustering.  On output, these may 
!    have been altered.
!
!    Input, integer ( kind = 4 ) NULL_CLUSTER_POLICY, specifies what to do if a
!    null cluster is encountered.
!    0, do nothing.
!    1, reset center of null cluster to a random data point;
!    2, reset center of null cluster to a random point in the data hull;
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) ci
  integer ( kind = 4 ) cj
  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  integer ( kind = 4 ), dimension ( cluster_num ) :: cluster_population
  real ( kind = 8 ) column_sum
  logical, parameter :: debug = .false.
  real ( kind = 8 ), dimension ( cluster_num ) :: distsq
  real ( kind = 8 ) energy
  real ( kind = 8 ), dimension ( point_num ) :: factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) list(1)
  integer ( kind = 4 ) null_cluster_num
  integer ( kind = 4 ) null_cluster_policy
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) swap
  integer ( kind = 4 ) swap_num
  integer ( kind = 4 ) swap_total
!
!  For each observation, calculate the distance from each cluster
!  center, and assign to the nearest.
!
  do 

    call nearest_cluster_raw ( dim_num, point_num, cluster_num, &
      cluster_center, point, cluster, swap_num )
!
!  Determine the cluster populations.
!
    cluster_population(1:cluster_num) = 0
    do i = 1, point_num
      j = cluster(i)
      cluster_population(j) = cluster_population(j) + 1
    end do

    if ( debug ) then
      call i4vec_print ( cluster_num, cluster_population, &
        '  KMEANS_RAW: Input Pops:' )
    end if
!
!  If a cluster is empty, give it a new cluster center, and restart
!  the process.
!  WARNING: This can take a long time!
!
    null_cluster_num = 0

    do j = 1, cluster_num
      if ( cluster_population(j) == 0 ) then
        null_cluster_num = null_cluster_num + 1
      end if
    end do

    if ( null_cluster_num == 0 ) then
      exit
    end if

    if ( null_cluster_policy == 0 ) then
      exit
    end if

    if ( debug ) then
      write ( *, '(a,i6)' ) &
        'KMEANS_RAW - Number of null clusters = ', null_cluster_num
    end if

    if ( null_cluster_policy == 1 ) then

      do j = 1, cluster_num
        if ( cluster_population(j) == 0 ) then
          k = i4_uniform ( 1, point_num, seed )
          cluster_center(1:dim_num,j) = point(1:dim_num,k)
        end if
      end do

    else if ( null_cluster_policy == 2 ) then

      do j = 1, cluster_num
        if ( cluster_population(j) == 0 ) then
          call random_number ( harvest = factor(1:point_num) )
          column_sum = sum ( factor(1:point_num) )
          factor(1:point_num) = factor(1:point_num) / column_sum
          cluster_center(1:dim_num,j) = &
            matmul ( point(1:dim_num,1:point_num), factor(1:point_num) )
        end if
      end do

    end if

  end do
!
!  Calculate the new cluster centers.
!
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00 

  do i = 1, point_num
    j = cluster(i)
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
      + point(1:dim_num,i)
  end do

  do i = 1, cluster_num
    cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
      / dble ( max ( cluster_population(i), 1 ) )
  end do
!
!  Carry out the iteration.
!
  it = 0
  swap_total = 0

  do
!
!  Determine the energy.
!
    call energy_raw ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, energy )

    if ( .false. ) then
      write ( *, '(i6,3x,14x,3x,g14.6)' ) it, energy
    end if

    if ( it_max <= it ) then
      exit
    end if

    it = it + 1

    swap = 0

    do i = 1, point_num

      ci = cluster(i)

      if ( cluster_population(ci) <= 1 ) then
        cycle
      end if

      do cj = 1, cluster_num

        if ( cj == ci ) then

          distsq(cj) = sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,cj) )**2 ) &
            * dble ( cluster_population(cj) ) &
            / dble ( cluster_population(cj) - 1 )

        else if ( cluster_population(cj) == 0 ) then

          cluster_center(1:dim_num,cj) = point(1:dim_num,i)
          distsq(cj) = 0.0D+00

        else

          distsq(cj) = sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,cj) )**2 ) &
            * dble ( cluster_population(cj) ) &
            / dble ( cluster_population(cj) + 1 )

        end if

      end do
!
!  Find the index of the minimum value of DISTSQ.
!
      list = minloc ( distsq(1:cluster_num) )
!
!  If that's not the cluster to which point I now belongs, move it there.
!
      if ( list(1) == ci ) then
        cycle
      end if

      cj = list(1)

      cluster_center(1:dim_num,ci) = &
        ( dble ( cluster_population(ci) ) * cluster_center(1:dim_num,ci) &
        - point(1:dim_num,i) ) / dble ( cluster_population(ci) - 1 )

      cluster_center(1:dim_num,cj) = &
        ( dble ( cluster_population(cj) ) * cluster_center(1:dim_num,cj) &
        + point(1:dim_num,i) ) / dble ( cluster_population(cj) + 1 )

      cluster_population(ci) = cluster_population(ci) - 1

      cluster_population(cj) = cluster_population(cj) + 1

      cluster(i) = cj

      swap = swap + 1
      swap_total = swap_total + 1

    end do

    if ( swap == 0 ) then
      exit
    end if

  end do
!
!  Determine the cluster populations.
!
  if ( debug ) then

    cluster_population(1:cluster_num) = 0
    do i = 1, point_num
      j = cluster(i)
      cluster_population(j) = cluster_population(j) + 1
    end do

    call i4vec_print ( cluster_num, cluster_population, &
      '  KMEANS_RAW: Output Pops:' )

  end if

  return
end
subroutine nearest_cluster_raw ( dim_num, point_num, cluster_num, &
  cluster_center, point, cluster, swap_num )

!*****************************************************************************80
!
!! NEAREST_CLUSTER_RAW finds the cluster nearest to a data point.
!
!  Discussion:
!
!    This routine uses the "raw" data.  Data is not normalized.
!    Distance is the usual Euclidean distance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of center points.
!
!    Input, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!    coordinates of the center points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points 
!    to be checked.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM).  On input, the cluster to
!    which each point was assigned.  On output, the cluster to which
!    each point has been reassigned.
!
!    Output, integer ( kind = 4 ) SWAP_NUM, the number of times a point was moved
!    from its input cluster to a different cluster.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_new
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nearest
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  integer ( kind = 4 ) swap_num

  swap_num = 0

  do i = 1, point_num

    dist = huge ( dist )
    nearest = 0

    do j = 1, cluster_num

      dist_new = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

      if ( dist_new < dist ) then
        dist = dist_new
        nearest = j
      end if

    end do

    if ( nearest /= cluster(i) ) then
      swap_num = swap_num + 1
    end if

    cluster(i) = nearest

  end do

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
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2*^1 - 1 )
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
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

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
subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_READ reads data from an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
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
!    18 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) TABLE(M,N), the data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  real ( kind = 8 ) table(m,n)
  real ( kind = 8 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
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
subroutine r8mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! R8MAT_HEADER_READ reads the header from an R8MAT file.
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
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

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
!    Input, real ( kind = 8 ) TABLE(M,N), the data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
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
function r8vec_norm2 ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM2 returns the 2-norm of a vector.
!
!  Discussion:
!
!    The vector 2-norm is defined as:
!
!      R8VEC_NORM2 = Sqrt ( Sum ( 1 <= I <= N ) A(I)**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose 2-norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM2, the 2-norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm2

  r8vec_norm2 = sqrt ( sum ( a(1:n)**2 ) )

  return
end
subroutine r8vec_range_input ( string, dim_num, value1, value2, ierror )

!*****************************************************************************80
!
!! R8VEC_RANGE_INPUT reads two DP vectors from the user, representing a range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of dimensions.
!
!    Output, real ( kind = 8 ) VALUE1(DIM_NUM), VALUE2(DIM_NUM), the values
!    entered by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) ierror
  character ( len = * ) string
  real ( kind = 8 ) value1(dim_num)
  real ( kind = 8 ) value2(dim_num)

  ierror = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )
  read ( *, *, iostat = ierror ) value1(1:dim_num), value2(1:dim_num)

  return
end
subroutine r8vec_unit_euclidean ( n, a )

!*****************************************************************************80
!
!! R8VEC_UNIT_EUCLIDEAN normalizes a N-vector in the Euclidean norm.
!
!  Discussion;
!
!    The euclidean norm is also sometimes called the l2 or
!    least squares norm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input/output, real ( kind = 8 ) A(N), the vector to be normalized.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) norm

  norm = sqrt ( sum ( a(1:n)**2 ) )

  if ( norm /= 0.0D+00 ) then
    a(1:n) = a(1:n) / norm
  end if

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
!    Input/output, integer ( kind = 4 ) SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator.
!
  implicit none

  integer ( kind = 4 ) count
  integer ( kind = 4 ) count_max
  integer ( kind = 4 ) count_rate
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
  real ( kind = 8 ) t
  logical, parameter :: verbose = .false.
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

    if ( verbose ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i12)' ) &
        '  Initialize RANDOM_NUMBER with user SEED = ', seed
    end if

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    if ( verbose ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i12)' ) &
        '  Initialize RANDOM_NUMBER with arbitrary SEED = ', seed
    end if

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
  integer ( kind = 4 ) get
  integer ( kind = 4 ) put
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  put = 0
  nchar = len_trim ( s )

  do get = 1, nchar

    c = s(get:get)

    if ( c /= ' ' .and. c /= TAB ) then
      put = put + 1
      s(put:put) = c
    end if

  end do

  s(put+1:nchar) = ' '

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
subroutine s_input ( string, value, ierror )

!*****************************************************************************80
!
!! S_INPUT prints a prompt string and reads a string from the user.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#'), or is blank,
!    the routine ignores that line, and tries to read the next one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, character ( len = * ) VALUE, the value input by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = * ) string
  character ( len = * ) value

  ierror = 0
  value = ' '
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) value

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'S_INPUT: Fatal error!'
      write ( *, '(a)' ) '  Input error!'
      stop
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( value(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( value ) == 0 ) then
      cycle
    end if

    exit

  end do

  return
end
subroutine s_rep_ch ( s, c1, c2 )

!*****************************************************************************80
!
!! S_REP_CH replaces all occurrences of one character by another.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string.
!
!    Input, character C1, C2, the character to be replaced, and the
!    replacement character.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  character ( len = * ) s

  do i = 1, len ( s )
    if ( s(i:i) == c1 ) then
      s(i:i) = c2
    end if
  end do

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
!    part of the real ( kind = 8 ) number.
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
!    ' 2D-1'           0.2
!    '23.45'           23.45
!    '-4.2D+2'         -420.0
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
!    characters can be read to form a legal number.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the value that was read from the string.
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

  character c
  logical ch_eqi
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
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
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
        rtop = 10.0D+00 * rtop + dble ( ndig )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + dble ( ndig )
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
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
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
!    07 September 2004
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
!    Output, integer IERROR, error flag.
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
  ierror = 0
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
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
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
!    Output, integer LAST, the last character of S used to make IVAL.
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
subroutine s_to_i4vec ( s, n, ivec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an I4VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2003
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
!    Output, integer ( kind = 4 ) IVEC(N), the values read from the string.
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
  integer ( kind = 4 ) ivec(n)
  integer ( kind = 4 ) length
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

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
