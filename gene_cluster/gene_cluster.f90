program main

!*****************************************************************************80
!
!! MAIN is the main program for GENE_CLUSTER.
!
!  Discussion:
!
!    GENE_CLUSTER groups gene expression data into clusters.
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
  implicit none

  integer ( kind = 4 ), parameter :: cluster_max = 10
  integer ( kind = 4 ) cluster_hi
  integer ( kind = 4 ) cluster_it_max
  integer ( kind = 4 ) cluster_lo
  integer ( kind = 4 ) cluster_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: cluster_size
  real ( kind = 8 ), allocatable, dimension ( : ) :: cluster_size_inv
  character ( len = 255 ) data_file
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: energy
  integer ( kind = 4 ) energy_it_max
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) inunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) normal
  real ( kind = 8 ), allocatable, dimension ( :, :) :: point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) seed

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GENE_CLUSTER'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Arrange gene expression data into clusters.'
  write ( *, '(a)' ) ' '

  seed = 1234567
!
!  Get the name of the input data file.
!
  call s_input ( 'Enter the name of the data file', data_file, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENE_CLUSTER'
    write ( *, '(a)' ) '  Error getting name of the data file.'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data will be read from the file ' // trim ( data_file )
!
!  Determine the spatial dimension (columns) and number of points (rows).
!
  call file_column_count ( data_file, dim_num )
  call file_line_count ( data_file, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension is     ', dim_num
  write ( *, '(a,i8)' ) '  Number of data points is ', point_num
!
!  Allocate space for the POINT array.
!
  allocate ( point(1:dim_num,1:point_num) )
!
!  Read the data.
!
  call get_unit ( inunit )

  open ( unit = inunit, file = data_file )
  do j = 1, point_num
    read ( inunit, * ) point(1:dim_num,j)
  end do
  close ( unit = inunit )
!
!  Get the range of cluster values to check.
!
  call i4_range_input ( 'Enter lower and upper number of clusters', &
    cluster_lo, cluster_hi, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENE_CLUSTER'
    write ( *, '(a)' ) '  Error getting cluster sizes.'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    stop
  end if

  if ( point_num < cluster_hi ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  CLUSTER_HI exceeds number of points.'
    cluster_hi = point_num
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program will try to determine the minimum energy'
  write ( *, '(a)' ) '  of a clustering, for cluster sizes in the range:'
  write ( *, '(2x,2i12)' ) cluster_lo, cluster_hi

  allocate ( cluster_size(1:cluster_hi) )
  allocate ( cluster_size_inv(1:cluster_hi) )
  allocate ( energy(1:cluster_hi) )

  do i = 1, cluster_hi
    cluster_size(i) = real ( i, kind = 8 )
  end do
  cluster_size_inv(1:cluster_hi) = 1.0D+00 / cluster_size(1:cluster_hi)
!
!  Get the number of different random starts.
!
  call i4_input ( 'Enter number of different random configurations to check', &
    cluster_it_max, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENE_CLUSTER'
    write ( *, '(a)' ) '  Error getting number of random configurations.'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For each number of clusters, the number of'
  write ( *, '(a)' ) '  distinct initial random configurations to be checked'
  write ( *, '(a,i8)' ) '  will be  ', cluster_it_max
!
!  Get the number of energy iterations for a particular random start:
!
  call i4_input ( 'Enter the number of energy iterations', &
    energy_it_max, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENE_CLUSTER'
    write ( *, '(a)' ) '  Error getting number of energy iterations.'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    stop
  end if
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For each initial random configuration, the number of'
  write ( *, '(a)' ) '  times the program will recompute the cluster centers,'
  write ( *, '(a,i8)' ) '  cluster components, and energy is ', energy_it_max
!
!  Analyze the un-normalized data.
!
  normal = 0

  call analysis ( normal, dim_num, point_num, point, cluster_lo, &
    cluster_hi, cluster_it_max, energy_it_max, energy, seed )

  cluster_num = cluster_hi + 1 - cluster_lo

  call data_to_gnuplot ( 'unnormal.txt', cluster_num, &
    cluster_size(cluster_lo:cluster_hi), energy(cluster_lo:cluster_hi) )

  call data_to_gnuplot ( 'unnormal2.txt', cluster_num, &
    cluster_size_inv(cluster_lo:cluster_hi), energy(cluster_lo:cluster_hi) )
!
!  Analyze the normalized data.
!
  normal = 1

  call analysis ( normal, dim_num, point_num, point, cluster_lo, & 
    cluster_hi, cluster_it_max, energy_it_max, energy, seed )

  call data_to_gnuplot ( 'normal.txt', cluster_num, &
    cluster_size(cluster_lo:cluster_hi), energy(cluster_lo:cluster_hi) )

  call data_to_gnuplot ( 'normal2.txt', cluster_num, &
    cluster_size_inv(cluster_lo:cluster_hi), energy(cluster_lo:cluster_hi) )
!
!  Free memory.
!
  deallocate ( cluster_size )
  deallocate ( cluster_size_inv )
  deallocate ( energy )
  deallocate ( point )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GENE_CLUSTER'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine analysis ( normal, dim_num, point_num, point, &
  cluster_lo, cluster_hi, cluster_it_max, energy_it_max, energy, seed )

!*****************************************************************************80
!
!! ANALYSIS computes the energy for a range of number of clusters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NORMAL:
!    0, analyze the unnormalized data.
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
!    Input, integer ( kind = 4 ) ENERGY_IT, the number of energy iterations.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) cluster_hi
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ), dimension (dim_num,cluster_hi) :: cluster_center
  integer ( kind = 4 ) cluster_it
  integer ( kind = 4 ) cluster_it_max
  integer ( kind = 4 ) cluster_lo
  integer ( kind = 4 ) cluster_num
  real ( kind = 8 ), dimension (cluster_hi) :: energy
  real ( kind = 8 ) energy_it
  real ( kind = 8 ) energy_max
  real ( kind = 8 ) energy_min
  integer ( kind = 4 ) energy_it_max
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) norm
  integer ( kind = 4 ) normal
  real ( kind = 8 ), dimension (dim_num,point_num) :: point
  real ( kind = 8 ), dimension (dim_num) :: r_min
  real ( kind = 8 ), dimension (dim_num) :: r_max
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
!
!  If requested, normalize the data.
!
  if ( normal == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Data is to be normalized.'
    write ( *, '(a)' ) ' '
    do j = 1, point_num
      norm = sqrt ( sum ( point(1:dim_num,j)**2 ) )
      if ( norm /= 0.0D+00 ) then
        point(1:dim_num,j) = point(1:dim_num,j) / norm
      end if
    end do
  end if

  if ( .false. ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Point data:'
    write ( *, '(a)' ) ' '
    call point_print ( dim_num, point_num, point )
  end if
!
!  Compute the minimum and maximum component values.
!
  do i = 1, dim_num
    r_min(i) = minval ( point(i,1:point_num) )
    r_max(i) = maxval ( point(i,1:point_num) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Minimum and maximum data values:'
  write ( *, '(a)' ) ' '
  do i = 1, dim_num
    write ( *, '(i3,2f10.4)' ) i, r_min(i), r_max(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cluster  Minimum    Maximum'
  write ( *, '(a)' ) '  Size     Energy     Energy'
  write ( *, '(a)' ) ' '
!
!  Consider a range of clusters.
!
  do cluster_num = cluster_lo, cluster_hi

    energy_min = huge ( energy_min )
    energy_max = 0.0D+00

    do cluster_it = 1, cluster_it_max

      call cluster_iteration ( dim_num, point_num, cluster_num, point, &
        r_min, r_max, seed, cluster_center, energy_it, energy_it_max  )

      energy_min = min ( energy_min, energy_it )
      energy_max = max ( energy_max, energy_it )

    end do

    energy(cluster_num) = energy_min

    write ( *, '(i9,2f12.4)' ) cluster_num, energy_min, energy_max

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Energy table:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Cluster              Energy'
  write ( *, '(a)' ) 'Size      Energy     /point    Sqrt(E/Pt)'
  write ( *, '(a)' ) ' '

  do cluster_num = cluster_lo, cluster_hi
    t = energy(cluster_num) / real ( point_num, kind = 8 )
    write ( *, '(i9,3f12.4)' ) cluster_num, energy(cluster_num), t, sqrt ( t )
  end do

  return
end
subroutine cluster_iteration ( dim_num, point_num, cluster_num, point, &
  r_min, r_max, seed, center, energy, energy_it_max )

!*****************************************************************************80
!
!! CLUSTER_ITERATION seeks the minimal energy of a cluster of a given size.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2005
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
!    Input, real ( kind = 8 ) R_MIN(DIM_NUM), R_MAX(DIM_NUM), the
!    minimum and maximum corners of the region.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) CENTER(DIM_NUM,CLUSTER_NUM), the centers
!    associated with the minimal energy clustering.
!
!    Output, real ( kind = 8 ) ENERGY, the total energy associated with the
!    minimal energy clustering.
!
!    Input, integer ( kind = 4 ) ENERGY_IT, the number of energy iterations.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: center
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: centroid
  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ) energy
  integer ( kind = 4 ) energy_it
  integer ( kind = 4 ) energy_it_max
  real ( kind = 8 ) energy_new
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) missed
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  real ( kind = 8 ), dimension ( dim_num ) :: r
  real ( kind = 8 ), dimension ( dim_num ) :: r_max
  real ( kind = 8 ), dimension ( dim_num ) :: r_min
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), dimension ( cluster_num ) :: subs
!
!  Initialize the centers randomly.
!
  do j = 1, cluster_num
    call random_number ( harvest = r(1:dim_num) )
    center(1:dim_num,j) = ( 1.0D+00 - r(1:dim_num) ) * r_min(1:dim_num) &
                                    + r(1:dim_num)   * r_max(1:dim_num)
  end do
!
!  Initialize the clusters randomly.
!
  subs(1:cluster_num) = 0
  do i = 1, point_num
    cluster(i) = i4_uniform ( 1, cluster_num, seed )
    j = cluster(i)
    subs(j) = subs(j) + 1
  end do

  call energy_computation ( dim_num, point_num, cluster_num, point, &
    center, cluster, energy )

  if ( .false. ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) &
      '  Initial total energy =                        ', energy
  end if

  if ( .false. ) then
    do j = 1, cluster_num
      write ( *, '(2i5,14f10.4)' ) j, subs(j), center(1:dim_num,j)
    end do
  end if

  do energy_it = 1, energy_it_max
!
!  #1:
!  Assign each point to the cluster of its nearest center.
!
    do i = 1, point_num
      call nearest_point ( dim_num, cluster_num, center, point(1,i), j )
      cluster(i) = j
    end do

    call energy_computation ( dim_num, point_num, cluster_num, point, &
      center, cluster, energy )

    if ( .false. ) then
      write ( *, '(a,g14.6)' ) &
        '  Total energy with old centers, new clusters = ', energy
    end if
!
!  #2:
!  Determine the centroids of the clusters.
!
    centroid(1:dim_num,1:cluster_num) = 0.0D+00
    subs(1:cluster_num) = 0

    do i = 1, point_num
      j = cluster(i)
      subs(j) = subs(j) + 1
      centroid(1:dim_num,j) = centroid(1:dim_num,j) + point(1:dim_num,i)
    end do

    missed = 0

    do j = 1, cluster_num

      if ( subs(j) /= 0 ) then
        centroid(1:dim_num,j) = centroid(1:dim_num,j) &
          / real ( subs(j), kind = 8 )
      else
        missed = missed + 1

        centroid(1:dim_num,j) = point(1:dim_num,missed)

!       call random_number ( harvest = r(1:dim_num) )
!       centroid(1:dim_num,j) = ( 1.0D+00 - r(1:dim_num) ) * r_min(1:dim_num) &
!                                         + r(1:dim_num)   * r_max(1:dim_num)

      end if

    end do

    if ( missed /= 0 ) then
      if ( .false. ) then
        write ( *, '(a,i8)' ) '  Number of empty clusters was ', missed
      end if
    end if

    call energy_computation ( dim_num, point_num, cluster_num, point, &
      centroid, cluster, energy_new )

    if ( .false. ) then
      write ( *, '(a,g14.6)' ) &
        '  Total energy with new centers, new clusters = ', energy_new
    end if

    if ( .false. ) then
      do j = 1, cluster_num
        write ( *, '(2i8,14f10.4)' ) j, subs(j), center(1:dim_num,j)
        write ( *, '(16x,14f10.4)' )             centroid(1:dim_num,j)
      end do
    end if
!
!  Update the centers.
!
    center(1:dim_num,1:cluster_num) = centroid(1:dim_num,1:cluster_num)

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

  character ( len = * )  file_name
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
subroutine energy_computation ( dim_num, point_num, cluster_num, point, &
  center, cluster, energy )

!*****************************************************************************80
!
!! ENERGY_COMPUTATION computes the total energy of a given clustering.
!
!  Discussion:
!
!    The total energy of a clustering is the sum of the energies of
!    each cluster.
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
!    09 August 2001
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
!    Input, real ( kind = 8 ) CENTER(DIM_NUM,CLUSTER_NUM), the center points.
!
!    Input, integer ( kind = 4 ) CLUSTER(POINT_NUM), indicates the cluster to which
!    each data point belongs.
!
!    Output, real ( kind = 8 ) ENERGY, the total energy of the clustering.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: center
  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ) energy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point

  energy = 0.0D+00
  do i = 1, point_num
    j = cluster(i)
    energy = energy + sum ( ( center(1:dim_num,j) - point(1:dim_num,i) )**2 )
  end do

  return
end
subroutine file_column_count ( file_name, ncolumn )

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
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NCOLUMN, the number of columns assumed 
!    to be in the file.
!
  implicit none

  character ( len = * )   file_name
  logical                 got_one
  integer ( kind = 4 )  ios
  integer ( kind = 4 )  iunit
  character ( len = 255 ) line
  integer ( kind = 4 )  ncolumn
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    ncolumn = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( file_name )
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( iunit, '(a)', iostat = ios ) line

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

    rewind ( iunit )

    do

      read ( iunit, '(a)', iostat = ios ) line

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

  close ( unit = iunit )

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
subroutine file_line_count ( file_name, nline )

!*****************************************************************************80
!
!! FILE_LINE_COUNT counts the number of lines in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Blank lines and comment lines, which begin with '#', are not counted.
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
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NLINE, the number of lines found in the file.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 256 ) line
  integer ( kind = 4 ) nline

  nline = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    nline = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:' // trim ( file_name )
    return
  end if
!
!  Count the lines.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    nline = nline + 1

  end do

  close ( unit = iunit )

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
subroutine i4_input ( string, value, ierror )

!*****************************************************************************80
!
!! I4_INPUT prints a prompt string and reads an integer from the user.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#') or is blank, 
!    the routine ignores that line, and tries to read the next one.
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
  character ( len = 80 ) line
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
  character ( len = 80 ) line
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
!    An I4 is an integer value.
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
subroutine nearest_point ( dim_num, cluster_num, center, point, nearest )

!*****************************************************************************80
!
!! NEAREST_POINT finds the center point nearest a data point.
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
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of center points.
!
!    Input, real ( kind = 8 ) CENTER(DIM_NUM,CLUSTER_NUM), the 
!    center points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM), the data point to be checked.
!
!    Output, integer ( kind = 4 ) NEAREST, the index of the center point closest to
!    the data point.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: center
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_new
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nearest
  real ( kind = 8 ), dimension ( dim_num ) :: point

  dist = huge ( dist )
  nearest = 0

  do j = 1, cluster_num

    dist_new = sum ( ( point(1:dim_num) - center(1:dim_num,j) )**2 )

    if ( dist_new < dist ) then
      dist = dist_new
      nearest = j
    end if

  end do

  return
end
subroutine point_generate ( point_dist, dim_num, r_min, r_max, point_num, &
  point )

!*****************************************************************************80
!
!! POINT_GENERATE generates data points for the problem.
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
!    Input, integer ( kind = 4 ) POINT_DIST, the point distribution to use.
!    1, use a uniform random distribution.
!    2, use a uniform grid of points.  (This hasn't been set up properly
!       except for 1D!).
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, real ( kind = 8 ) R_MIN(DIM_NUM), R_MAX(DIM_NUM), the
!    minimum and maximum corners of the region.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points to generate.
!
!    Output, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) i
  real ( kind = 8 ) point(dim_num,point_num)
  integer ( kind = 4 ) point_dist
  real ( kind = 8 ) r
  real ( kind = 8 ) r_max(dim_num)
  real ( kind = 8 ) r_min(dim_num)

  if ( point_dist == 1 ) then

    call random_number ( harvest = point(1:dim_num,1:point_num) )

    do i = 1, dim_num
      point(i,1:point_num) = r_min(i) + &
        point(i,1:point_num) * ( r_max(i) - r_min(i) )
    end do

  else if ( point_dist == 2 ) then

    do i = 1, point_num

      if ( 1 < point_num ) then
        r = real (         i - 1, kind = 8 ) &
          / real ( point_num - 1, kind = 8 )
      else
        r = 0.5D+00
      end if

      point(1:dim_num,i) = ( 1.0D+00 - r ) * r_min(1:dim_num) &
                                   + r   * r_max(1:dim_num)
    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINT_GENERATE - Fatal error!'
    write ( *, '(a)' ) '  Meaningless input value of point distribution,'
    write ( *, '(a,i8)' ) '  POINT_DIST = ', point_dist
    stop

  end if

  return
end
subroutine point_print ( dim_num, point_num, point )

!*****************************************************************************80
!
!! POINT_PRINT prints out the values of the data points.
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
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) i
  real ( kind = 8 ) point(dim_num,point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data points:'
  write ( *, '(a)' ) ' '

  do i = 1, point_num
    write ( *, '(2x,i8)' ) i
    write ( *, '(2x,8f10.4)' ) point(1:dim_num,i)
  end do

  return
end
subroutine s_input ( string, value, ierror )

!*****************************************************************************80
!
!! S_INPUT prints a prompt string and reads a string from the user.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#') or is blank, 
!    the routine ignores that line, and tries to read the next one.
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
!    Output, character ( len = * ) VALUE, the value input by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = 80 ) line
  character ( len = * ) string
  character ( len = * ) value

  ierror = 0
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) line

    if ( ierror /= 0 ) then
      value = 'S_INPUT: Input error!'
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

    value = line
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
!    Output, integer ( kind = 4 ) LAST, the last character of S used to make IVAL.
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
