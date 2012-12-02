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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding value.  If C was
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
subroutine cluster_energy_compute ( dim_num, point_num, cluster_num, point, &
  cluster, cluster_center, cluster_energy )

!*****************************************************************************80
!
!! CLUSTER_ENERGY_COMPUTE computes the energy of the clusters.
!
!  Discussion:
!
!    The cluster energy is defined as the sum of the distance
!    squared from each point to its cluster center.  It is the goal
!    of the H-means and K-means algorithms to find, for a fixed number
!    of clusters, a clustering that minimizes this energy
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 March 2002
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
!    Input, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster to which each
!    data point belongs.
!
!    Input, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the 
!    centers associated with the minimal energy clustering.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy
!    associated with each cluster.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  real ( kind = 8 ), dimension ( cluster_num ) :: cluster_energy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  real ( kind = 8 ) point_energy

  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    j = cluster(i)

    point_energy = sum ( &
      ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

    cluster_energy(j) = cluster_energy(j) + point_energy

  end do

  return
end
subroutine cluster_initialize_1 ( dim_num, point_num, cluster_num, point, &
  cluster_center )

!*****************************************************************************80
!
!! CLUSTER_INITIALIZE_1 initializes the clusters to data points.
!
!  Discussion:
!
!    The cluster centers are simply chosen to be the first data points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2002
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
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  integer ( kind = 4 ) i
  real ( kind = 8 ) point(dim_num,point_num)

  do i = 1, cluster_num
    cluster_center(1:dim_num,i) = point(1:dim_num,i)
  end do

  return
end
subroutine cluster_initialize_2 ( dim_num, point_num, cluster_num, point, &
  seed, cluster_center )

!*****************************************************************************80
!
!! CLUSTER_INITIALIZE_2 initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, the hyperbox containing the data is computed.
!
!    Then the cluster centers are chosen uniformly at random within
!    this hyperbox.
!
!    Of course, if the data is not smoothly distributed throughout
!    the box, many cluster centers will be isolated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2006
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
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), 
!    the coordinates of the cluster centers.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  integer ( kind = 4 ) i
  real ( kind = 8 ) point(dim_num,point_num)
  real ( kind = 8 ) r(dim_num)
  real ( kind = 8 ) r_max(dim_num)
  real ( kind = 8 ) r_min(dim_num)
  integer ( kind = 4 ) seed

  r_min = minval ( point, 2 )
  r_max = maxval ( point, 2 )

  do i = 1, cluster_num

    call r8vec_uniform_01 ( dim_num, seed, r )

    cluster_center(1:dim_num,i) = &
      ( 1.0D+00 - r(1:dim_num) ) * r_min(1:dim_num) &
                + r(1:dim_num)   * r_max(1:dim_num)
  end do
    
  return
end
subroutine cluster_initialize_3 ( dim_num, point_num, cluster_num, point, &
  seed, cluster_center )

!*****************************************************************************80
!
!! CLUSTER_INITIALIZE_3 initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, each point is randomly assigned to a cluster, and
!    the cluster centers are then computed as the centroids of the points 
!    in the cluster.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2006
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
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), 
!    the coordinates of the cluster centers.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  integer ( kind = 4 ) cluster_population(cluster_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(dim_num,point_num)
  integer ( kind = 4 ) seed
!
!  Assign one point to each cluster center.
!
  do i = 1, cluster_num
    cluster_center(1:dim_num,i) = point(1:dim_num,i)
  end do

  cluster_population(1:cluster_num) = 1
!
!  The rest of the points get assigned randomly.
!
  do i = cluster_num+1, point_num
    j = i4_uniform ( 1, cluster_num, seed )
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) + &
      point(1:dim_num,i)
    cluster_population(j) = cluster_population(j) + 1
  end do
!
!  Now average the points to get the centroid.
!
  do i = 1, cluster_num
    cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) / &
      real ( cluster_population(i), kind = 8 )
  end do
    
  return
end
subroutine cluster_initialize_4 ( dim_num, point_num, cluster_num, point, &
  seed, cluster_center )

!*****************************************************************************80
!
!! CLUSTER_INITIALIZE_4 initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, each data point is divided randomly among the
!    the cluster centers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2002
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
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), 
!    the coordinates of the cluster centers.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  real ( kind = 8 ) cluster_factor(cluster_num)
  real ( kind = 8 ) cluster_weight(cluster_num)
  real ( kind = 8 ) divisor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(dim_num,point_num)
  integer ( kind = 4 ) seed

  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00
  cluster_weight(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    call r8vec_uniform_01 ( cluster_num, seed, cluster_factor )

    divisor = sum ( cluster_factor(1:cluster_num) )
    cluster_factor(1:cluster_num) = cluster_factor(1:cluster_num) / divisor

    do j = 1, cluster_num
      cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
        + cluster_factor(j) * point(1:dim_num,i)
    end do

    cluster_weight(1:cluster_num) = cluster_weight(1:cluster_num) &
      + cluster_factor(1:cluster_num)

  end do
!
!  Now normalize,  so that each cluster center is now a convex 
!  combination of the points.
!
  do i = 1, cluster_num
    cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
      / cluster_weight(i)
  end do

  return
end
subroutine cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
  seed, cluster_center )

!*****************************************************************************80
!
!! CLUSTER_INITIALIZE_5 initializes the cluster centers to random values.
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
!    06 June 2006
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
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  real ( kind = 8 ) column_sum
  real ( kind = 8 ) factor(point_num,cluster_num)
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(dim_num,point_num)
  integer ( kind = 4 ) seed
!
!  Get a PxC block of random factors.
!
  call r8mat_uniform_01 ( point_num, cluster_num, seed, factor )
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

  return
end
subroutine cluster_print_summary ( point_num, cluster_num, &
  cluster_population, cluster_energy, cluster_variance )

!*****************************************************************************80
!
!! CLUSTER_PRINT_SUMMARY prints a summary of data about a clustering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number of
!    points assigned to each cluster.
!
!    Input, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy of 
!    the clusters.
!
!    Input, real ( kind = 8 ) CLUSTER_VARIANCE(CLUSTER_NUM), the variance of 
!    the clusters.
!
  implicit none

  integer ( kind = 4 ) cluster_num

  real ( kind = 8 ) ce
  integer ( kind = 4 ) cep
  real ( kind = 8 ) ce_total
  real ( kind = 8 ) cluster_energy(cluster_num)
  integer ( kind = 4 ) cluster_population(cluster_num)
  real ( kind = 8 ) cluster_variance(cluster_num)
  integer ( kind = 4 ) cp
  integer ( kind = 4 ) cpp
  real ( kind = 8 ) cv
  integer ( kind = 4 ) i
  integer ( kind = 4 ) point_num

  ce_total = sum ( cluster_energy(1:cluster_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Clustering statistics:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '    Number of clusters is ', cluster_num
  write ( *, '(a,i8)' ) '    Number of points is   ', point_num
  write ( *, '(a,g14.6)' ) '    Total energy is       ', ce_total
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Cluster   Population        Energy          Variance'
  write ( *, '(a)' ) &
    '    -------  -----------  -----------------  --------------'
  write ( *, '(a)' ) '                  #    %     value        %'
  write ( *, '(a)' ) ' '

  do i = 1, cluster_num
    cp = cluster_population(i)
    cpp = int ( real ( 100 * cp, kind = 8 ) / real ( point_num, kind = 8 ) )
    ce = cluster_energy(i)
    cep = int ( ( ce * 100.0D+00 ) / ce_total )
    cv = cluster_variance(i)
    write ( *, '(6x,i3,2x,i8,2x,i3,g14.6,2x,i3,2x,g14.6)' ) &
      i, cp, cpp, ce, cep, cv
  end do

  cp = sum ( cluster_population(1:cluster_num) )
  cpp = 100
  ce = sum ( cluster_energy(1:cluster_num) )
  cep = 100
  cv = sum ( cluster_population(1:cluster_num) &
    * cluster_variance(1:cluster_num) ) / cp

  write ( *, '(a)' ) ' '
  write ( *, '(a9,2x,i8,2x,i3,g14.6,2x,i3,2x,g14.6)' ) &
    '    Total', cp, cpp, ce, cep, cv

  return
end
subroutine cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
  cluster, cluster_center, cluster_variance )

!*****************************************************************************80
!
!! CLUSTER_VARIANCE_COMPUTE computes the variance of the clusters.
!
!  Discussion:
!
!    The cluster variance (from the cluster center) is the average of the 
!    sum of the squares of the distances of each point in the cluster to the 
!    cluster center.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2009
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
!    Input, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster to which each
!    data point belongs.
!
!    Input, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the 
!    centers associated with the minimal energy clustering.
!
!    Output, real ( kind = 8 ) CLUSTER_VARIANCE(CLUSTER_NUM), the variance
!    associated with each cluster.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) cluster(point_num)
  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  integer ( kind = 4 ) cluster_population(cluster_num)
  real ( kind = 8 ) cluster_variance(cluster_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(dim_num,point_num)
  real ( kind = 8 ) point_variance

  cluster_population(1:cluster_num) = 0
  cluster_variance(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    j = cluster(i)

    point_variance = sum ( &
      ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

    cluster_variance(j) = cluster_variance(j) + point_variance
    cluster_population(j) = cluster_population(j) + 1

  end do

  cluster_variance(1:cluster_num) = cluster_variance(1:cluster_num) &
    / cluster_population(1:cluster_num)

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
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer   ( kind = 4 )   column_num
  logical                  got_one
  character ( len = * )    input_filename
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  character ( len = 255 )  line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = input_status )

  if ( input_status /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' &
      // trim ( input_filename ) // '" on unit ', input_unit
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
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

      read ( input_unit, '(a)', iostat = input_status ) line

      if ( input_status /= 0 ) then
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

  integer   ( kind = 4 )   bad_num
  integer   ( kind = 4 )   comment_num
  integer   ( kind = 4 )   ierror
  character ( len = * )    input_filename
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  character ( len = 255 )  line
  integer   ( kind = 4 )   record_num
  integer   ( kind = 4 )   row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
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
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
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
subroutine hmeans_01 ( dim_num, point_num, cluster_num, it_max, it_num, point, &
  cluster, cluster_center, cluster_population, cluster_energy )

!*****************************************************************************80
!
!! HMEANS_01 applies the H-Means algorithm.
!
!  Discussion:
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
!      Energy = Sum ( 1 <= I <= N ) || X(I) - Z(CLUSTER(I)) ||**2
!
!    where
!
!      || X - Z ||**2 = Sum ( 1 <= J <= M ) ( X(J) - Z(J) )**2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez, Angel Martinez,
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
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM).  On input, the user 
!    may specify an initial cluster for each point, or leave all entrie of
!    CLUSTER set to 0.  On output, CLUSTER contains the index of the
!    cluster to which each data point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the 
!    centers associated with the minimal energy clustering.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM),
!    the populuation of each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy
!    associated with each cluster.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) c
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: centroid
  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  real ( kind = 8 ), dimension ( cluster_num ) :: cluster_energy
  integer ( kind = 4 ), dimension ( cluster_num ) :: cluster_population
  logical, parameter :: debug = .true.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) missed
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  real ( kind = 8 ) point_energy
  real ( kind = 8 ) point_energy_min
  integer ( kind = 4 ) swap
!
!  Data checks.
!
  if ( cluster_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_01 - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_01 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if

  if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_01 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop
  end if

  if ( it_max < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_01 - Fatal error!'
    write ( *, '(a)' ) '  IT_MAX < 0.'
    stop
  end if
!
!  On input, legal entries in CLUSTER are preserved, but
!  otherwise, each point is assigned to its nearest cluster.
!
  do i = 1, point_num
    if ( cluster(i) <= 0 .or. cluster_num < cluster(i) ) then

      point_energy_min = huge ( point_energy_min )

      do j = 1, cluster_num

        point_energy = sum ( &
          ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

        if ( point_energy < point_energy_min ) then
          point_energy_min = point_energy
          cluster(i) = j
        end if

      end do

    end if
  end do
  it_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1
!
!  #1:
!  Assign each point to the cluster of its nearest center.
!
    swap = 0

    do i = 1, point_num

      point_energy_min = huge ( point_energy_min )
      c = cluster(i)

      do j = 1, cluster_num

        point_energy = sum ( &
          ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

        if ( point_energy < point_energy_min ) then
          point_energy_min = point_energy
          cluster(i) = j
        end if

      end do

      if ( c /= cluster(i) ) then
        swap = swap + 1
      end if
      
    end do
!
!  Terminate if no points were swapped.
!
    if ( 1 < it_num ) then
      if ( swap == 0 ) then
        exit
      end if
    end if
!
!  #2:
!  Determine the total energy of the new clustering with current centroids.
!
    cluster_energy(1:cluster_num) = 0.0D+00

    do i = 1, point_num

      c = cluster(i)

      cluster_energy(c) = cluster_energy(c) + sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,c) )**2 )

    end do
!
!  #3:
!  Determine the centroids of the clusters.
!
    centroid(1:dim_num,1:cluster_num) = 0.0D+00
    cluster_population(1:cluster_num) = 0

    do i = 1, point_num
      j = cluster(i)
      cluster_population(j) = cluster_population(j) + 1
      centroid(1:dim_num,j) = centroid(1:dim_num,j) + point(1:dim_num,i)
    end do
!
!  Now divide by the population to get the centroid.
!  But if a center has no population, pick a point at random.
!
    missed = 0

    do j = 1, cluster_num

      if ( cluster_population(j) /= 0 ) then
        centroid(1:dim_num,j) = centroid(1:dim_num,j) &
          / real ( cluster_population(j), kind = 8 )
      else
        missed = missed + 1

        centroid(1:dim_num,j) = point(1:dim_num,missed)

      end if

    end do

    cluster_center(1:dim_num,1:cluster_num) = centroid(1:dim_num,1:cluster_num)
!
!  #4:
!  Determine the total energy of the current clustering with new centroids.
!
    cluster_energy(1:cluster_num) = 0.0D+00

    do i = 1, point_num

      c = cluster(i)

      cluster_energy(c) = cluster_energy(c) + sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,c) )**2 )

    end do
 
  end do

  return
end
subroutine hmeans_02 ( dim_num, point_num, cluster_num, it_max, it_num, point, &
  cluster, cluster_center, cluster_population, cluster_energy, seed )

!*****************************************************************************80
!
!! HMEANS_02 applies the H-Means algorithm.
!
!  Discussion:
!
!    This is a simple routine to group a set of points into K clusters,
!    each with a center point, in such a way that the total cluster 
!    energy is minimized.  The total cluster energy is the sum of the
!    squares of the distances of each point to the center of its cluster.
!
!    The algorithm begins with an initial estimate for the cluster centers:
!
!    1. The points are assigned to the nearest cluster centers.
!
!    2. The iteration stops if the total energy has not changed 
!        significantly, or we have reached the maximum number of iterations.
!
!    3. Each cluster center is replaced by the centroid of the points
!       in the cluster.
!
!    4. Return to step 1.
!
!    The algorithm may fail to find the best solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2005
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
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM).  On input, the user 
!    may specify an initial cluster for each point, or leave all entrie of
!    CLUSTER set to 0.  On output, CLUSTER contains the index of the
!    cluster to which each data point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number of
!    points assigned to each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy of 
!    the clusters.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) c
  integer ( kind = 4 ) cluster(point_num)
  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  real ( kind = 8 ) cluster_energy(cluster_num)
  integer ( kind = 4 ) cluster_population(cluster_num)
  logical, parameter :: debug = .false.
  real ( kind = 8 ) energy(cluster_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) list(1)
  real ( kind = 8 ) point(dim_num,point_num)
  real ( kind = 8 ) point_energy
  real ( kind = 8 ) point_energy_min
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) swap
!
!  Data checks.
!
  if ( cluster_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_02 - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_02 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if

  if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_02 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop
  end if

  if ( it_max < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_02 - Fatal error!'
    write ( *, '(a)' ) '  IT_MAX < 0.'
    stop
  end if
!
!  On input, legal entries in CLUSTER are preserved, but
!  otherwise, each point is assigned to its nearest cluster.
!
  do i = 1, point_num
    if ( cluster(i) <= 0 .or. cluster_num < cluster(i) ) then

      point_energy_min = huge ( point_energy_min )

      do j = 1, cluster_num

        point_energy = sum ( &
          ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

        if ( point_energy < point_energy_min ) then
          point_energy_min = point_energy
          cluster(i) = j
        end if

      end do

    end if
  end do

  it_num = 0

  do
!
!  Given centers, assign points to nearest center.
!
    cluster_population(1:cluster_num) = 0
    cluster_energy(1:cluster_num) = 0.0D+00

    swap = 0

    do i = 1, point_num

      do j = 1, cluster_num
        energy(j) = sum ( &
          ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
      end do

      list = minloc ( energy(1:cluster_num) )

      c = list(1)

      if ( c /= cluster(i) ) then
        swap = swap + 1
      end if

      cluster(i) = c
      cluster_energy(c) = cluster_energy(c) + energy(c) 
      cluster_population(c) = cluster_population(c) + 1

    end do

    if ( debug ) then
      write ( *, '(i3,g14.6)' ) it_num, sum ( cluster_energy(1:cluster_num) )
    end if

    if ( 0 < it_num ) then
      if ( swap == 0 ) then
        exit
      end if
    end if

    if ( it_max <= it_num ) then
      exit
    end if

    it_num = it_num + 1
!
!  Given points in cluster, replace center by centroid.
!
    cluster_center(1:dim_num,1:cluster_num) = 0.0D+00
  
    do i = 1, point_num
      c = cluster(i)
      cluster_center(1:dim_num,c) = cluster_center(1:dim_num,c) &
        + point(1:dim_num,i)
    end do

    do i = 1, cluster_num
      if ( cluster_population(i) /= 0 ) then
        cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) / &
          real ( cluster_population(i), kind = 8 )
      else
        j = i4_uniform ( 1, point_num, seed )
        cluster_center(1:dim_num,i) = point(1:dim_num,j)
      end if
    end do

  end do
!
!  Compute the energy based on the final value of the cluster centers.
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    c = cluster(i)

    cluster_energy(c) = cluster_energy(c) + sum ( &
      ( point(1:dim_num,i) - cluster_center(1:dim_num,c) )**2 )

  end do

  return
end
subroutine hmeans_w_01 ( dim_num, point_num, cluster_num, it_max, it_num, &
  point, weight, cluster, cluster_center, cluster_population, cluster_energy )

!*****************************************************************************80
!
!! HMEANS_W_01 applies the weighted H-Means algorithm. 
!
!  Discussion:
!
!    The input data for the weight H-Means problem includes:
!    * a set of N data points X in M dimensions, 
!    * a set of N nonnegative weights W,
!    * a desired number of clusters K.
!    * an initial set of cluster centers Z,
!    * an (optional) initial set of cluster assignments.
!
!    The goal is to determine K points Z, called cluster centers, and
!    to assign each point X(I) to some cluster Z(J), so that we minimize
!    the weighted standard deviation of the distance of each data point
!    to the center of its cluster.  Writing J = CLUSTER(I) to
!    indicate the index of the nearest cluster center Z(J) to the 
!    point X(I), the quantity we are trying to minimize is the sum
!    of the weighted cluster energies E(J), where:
!
!      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||**2
!
!    Here, we assume that we are using the Euclidean norm, so that
!    
!      || X(I) - Z(J) ||**2 = Sum ( 1 <= K <= M )
!        ( X(I)(K) - Z(J)(K) )**2
!
!    In this notation, X(I)(K) is the K-th spatial component of the
!    I-th data point.
!
!    Note that this routine should give the same results as HMEANS_01
!    in any case in which all the entries of the WEIGHT vector are equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez, Angel Martinez,
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
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, real ( kind = 8 ) WEIGHT(POINT_NUM), the weights
!    assigned to the data points.  These must be nonnegative, and
!    at least one must be strictly positive.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM).  On input, the user 
!    may specify an initial cluster for each point, or leave all entrie of
!    CLUSTER set to 0.  On output, CLUSTER contains the index of the
!    cluster to which each data point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!    centers associated with the minimal energy clustering.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM),
!    the populuation of each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy
!    associated with each cluster.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) c
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: centroid
  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  real ( kind = 8 ), dimension ( cluster_num ) :: cluster_energy
  integer ( kind = 4 ), dimension ( cluster_num ) :: cluster_population
  real ( kind = 8 ), dimension ( cluster_num ) :: cluster_weight
  logical, parameter :: debug = .true.
  real ( kind = 8 ) energy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) missed
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  real ( kind = 8 ) point_energy
  real ( kind = 8 ) point_energy_min
  integer ( kind = 4 ) swap
  real ( kind = 8 ) weight(point_num)
!
!  Data checks.
!
  if ( cluster_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if

  if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop
  end if

  if ( it_max < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  IT_MAX < 0.'
    stop
  end if

  if ( any ( weight(1:point_num) < 0.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  Some weight entry is negative.'
    stop
  end if

  if ( all ( weight(1:point_num) <= 0.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  No weight entry is positive.'
    stop
  end if
!
!  On input, legal entries in CLUSTER are preserved, but
!  otherwise, each point is assigned to its nearest cluster.
!
  do i = 1, point_num
    if ( cluster(i) <= 0 .or. cluster_num < cluster(i) ) then

      point_energy_min = huge ( point_energy_min )

      do j = 1, cluster_num

        point_energy = sum ( &
          ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

        if ( point_energy < point_energy_min ) then
          point_energy_min = point_energy
          cluster(i) = j
        end if

      end do

    end if
  end do

  it_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1
!
!  #1:
!  Reassign points to clusters:
!  Assign each point to the cluster whose center is nearest;
!  Count the number of points whose cluster assignment is changed.
!
    swap = 0

    do i = 1, point_num

      point_energy_min = huge ( point_energy_min )
      c = cluster(i)

      do j = 1, cluster_num

        point_energy = sum ( &
          ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

        if ( point_energy < point_energy_min ) then
          point_energy_min = point_energy
          cluster(i) = j
        end if

      end do

      if ( c /= cluster(i) ) then
        swap = swap + 1
      end if

    end do
!
!  If no point changed its cluster assignment, the algorithm can make no 
!  more improvements, so terminate.
!
    if ( 1 < it_num ) then
      if ( swap == 0 ) then
        exit
      end if
    end if
!
!  Determine the current energy.
!
    energy = 0.0
    do i = 1, point_num
      energy = energy + weight(i) * sum ( &
          ( point(1:dim_num,i) - cluster_center(1:dim_num,cluster(i)) )**2 )
    end do
    write ( *, * ) it_num, energy
!
!  #2:
!  Determine the centroids of the clusters, and set the 
!  cluster center to the cluster centroid.
!
    centroid(1:dim_num,1:cluster_num) = 0.0D+00
    cluster_population(1:cluster_num) = 0
    cluster_weight(1:cluster_num) = 0.0D+00

    do i = 1, point_num
      c = cluster(i)
      cluster_population(c) = cluster_population(c) + 1
      cluster_weight(c) = cluster_weight(c) + weight(i)
      centroid(1:dim_num,c) = centroid(1:dim_num,c) &
        + weight(i) * point(1:dim_num,i)
    end do

    missed = 0

    do c = 1, cluster_num

      if ( cluster_weight(c) /= 0.0D+00 ) then
        centroid(1:dim_num,c) = centroid(1:dim_num,c) / cluster_weight(c)
      else
        missed = missed + 1
        centroid(1:dim_num,c) = point(1:dim_num,missed)
      end if

    end do

    cluster_center(1:dim_num,1:cluster_num) = centroid(1:dim_num,1:cluster_num)

  end do
!
!  Compute the energy based on the final value of the cluster centers.
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    c = cluster(i)

    cluster_energy(c) = cluster_energy(c) + weight(i) * sum ( &
      ( point(1:dim_num,i) - cluster_center(1:dim_num,c) )**2 )

  end do

  return
end
subroutine hmeans_w_02 ( dim_num, point_num, cluster_num, it_max, it_num, &
  point, weight, cluster, cluster_center, cluster_population, &
  cluster_energy, seed )

!*****************************************************************************80
!
!! HMEANS_W_02 applies the weighted H-Means algorithm.
!
!  Discussion:
!
!    The input data for the weight H-Means problem includes:
!    * a set of N data points X in M dimensions, 
!    * a set of N nonnegative weights W,
!    * a desired number of clusters K.
!    * an initial set of cluster centers Z,
!    * an (optional) initial set of cluster assignments.
!
!    The goal is to determine K points Z, called cluster centers, and
!    to assign each point X(I) to some cluster Z(J), so that we minimize
!    the weighted standard deviation of the distance of each data point
!    to the center of its cluster.  Writing J = CLUSTER(I) to
!    indicate the index of the nearest cluster center Z(J) to the 
!    point X(I), the quantity we are trying to minimize is the sum
!    of the weighted cluster energies E(J), where:
!
!      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||**2
!
!    Here, we assume that we are using the Euclidean norm, so that
!    
!      || X(I) - Z(J) ||**2 = Sum ( 1 <= K <= M )
!        ( X(I)(K) - Z(J)(K) )**2
!
!    In this notation, X(I)(K) is the K-th spatial component of the
!    I-th data point.
!
!    Note that this routine should give the same results as HMEANS_02
!    in any case in which all the entries of the WEIGHT vector are equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2006
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
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM).  On input, the user 
!    may specify an initial cluster for each point, or leave all entrie of
!    CLUSTER set to 0.  On output, CLUSTER contains the index of the
!    cluster to which each data point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number of
!    points assigned to each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy of 
!    the clusters.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) c
  integer ( kind = 4 ) cluster(point_num)
  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  real ( kind = 8 ) cluster_energy(cluster_num)
  integer ( kind = 4 ) cluster_population(cluster_num)
  real ( kind = 8 ) cluster_weight(cluster_num)
  logical, parameter :: debug = .false.
  real ( kind = 8 ) energy(cluster_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) list(1)
  real ( kind = 8 ) point(dim_num,point_num)
  real ( kind = 8 ) point_energy
  real ( kind = 8 ) point_energy_min
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) swap
  real ( kind = 8 ) weight(point_num)
!
!  Data checks.
!
  if ( cluster_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_W_02 - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_W_02 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if

  if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_W_02 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop
  end if

  if ( it_max < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_W_02 - Fatal error!'
    write ( *, '(a)' ) '  IT_MAX < 0.'
    stop
  end if

  if ( any ( weight(1:point_num) < 0.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_W_02 - Fatal error!'
    write ( *, '(a)' ) '  Some weight entry is negative.'
    stop
  end if

  if ( all ( weight(1:point_num) <= 0.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HMEANS_W_02 - Fatal error!'
    write ( *, '(a)' ) '  No weight entry is positive.'
    stop
  end if
!
!  On input, legal entries in CLUSTER are preserved, but
!  otherwise, each point is assigned to its nearest cluster.
!
  do i = 1, point_num
    if ( cluster(i) <= 0 .or. cluster_num < cluster(i) ) then

      point_energy_min = huge ( point_energy_min )

      do j = 1, cluster_num

        point_energy = sum ( &
          ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

        if ( point_energy < point_energy_min ) then
          point_energy_min = point_energy
          cluster(i) = j
        end if

      end do

    end if
  end do

  it_num = 0

  do
!
!  Given centers, assign points to nearest center.
!
    cluster_population(1:cluster_num) = 0
    cluster_weight(1:cluster_num) = 0.0D+00
    cluster_energy(1:cluster_num) = 0.0D+00

    swap = 0

    do i = 1, point_num

      do j = 1, cluster_num
        energy(j) = sum ( &
          ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
      end do

      list = minloc ( energy(1:cluster_num) )

      c = list(1)

      if ( c /= cluster(i) ) then
        swap = swap + 1
      end if

      cluster(i) = c
      cluster_energy(c) = cluster_energy(c) + weight(i) * energy(c) 
      cluster_population(c) = cluster_population(c) + 1
      cluster_weight(c) = cluster_weight(c) + weight(i)

    end do

    if ( debug ) then
      write ( *, '(i3,g14.6)' ) it_num, sum ( cluster_energy(1:cluster_num) )
    end if

    if ( 0 < it_num ) then
      if ( swap == 0 ) then
        exit
      end if
    end if

    if ( it_max <= it_num ) then
      exit
    end if

    it_num = it_num + 1
!
!  Given points in cluster, replace center by weighted centroid.
!
    cluster_center(1:dim_num,1:cluster_num) = 0.0D+00
  
    do i = 1, point_num
      c = cluster(i)
      cluster_center(1:dim_num,c) = cluster_center(1:dim_num,c) &
        + weight(i) * point(1:dim_num,i)
    end do

    do i = 1, cluster_num
      if ( cluster_weight(i) /= 0 ) then
        cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) / &
          cluster_weight(i)
      else
        j = i4_uniform ( 1, point_num, seed )
        cluster_center(1:dim_num,i) = point(1:dim_num,j)
      end if
    end do

  end do
!
!  Compute the energy based on the final value of the cluster centers.
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    c = cluster(i)

    cluster_energy(c) = cluster_energy(c) + weight(i) * sum ( &
      ( point(1:dim_num,i) - cluster_center(1:dim_num,c) )**2 )

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
subroutine i4mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_WRITE writes an I4MAT file.
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
!    Input, integer ( kind = 4 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_filename
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  integer   ( kind = 4 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
  write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine kmeans_01 ( dim_num, point_num, cluster_num, it_max, it_num, point, &
  cluster, cluster_center, cluster_population, cluster_energy )

!*****************************************************************************80
!
!! KMEANS_01 applies the K-Means algorithm.
!
!  Discussion:
!
!    Given a matrix of POINT_NUM observations on DIM_NUM variables, the
!    observations are to be allocated to CLUSTER_NUM clusters in such 
!    a way that the within-cluster sum of squares is minimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    Original FORTRAN77 version by David Sparks.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Sparks,
!    Algorithm AS 58: 
!    Euclidean Cluster Analysis,
!    Applied Statistics,
!    Volume 22, Number 1, 1973, pages 126-130.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the points.
!
!    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), indicates which cluster
!    each point belongs to.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number 
!    of points in each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the 
!    cluster energies.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) cluster(point_num)
  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  real ( kind = 8 ) cluster_energy(cluster_num)
  integer ( kind = 4 ) cluster_population(cluster_num)
  real ( kind = 8 ) dc
  real ( kind = 8 ) de
  real ( kind = 8 ) f(point_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) list(1)
  real ( kind = 8 ) point(dim_num,point_num)
  integer ( kind = 4 ) swap

  it_num = 0
!
!  Idiot checks.
!
  if ( cluster_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if

  if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop
  end if
!
!  For each observation, calculate the distance from each cluster
!  center, and assign to the nearest.
!
  do i = 1, point_num

    do j = 1, cluster_num
      cluster_energy(j) = sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
    end do

    list = minloc ( cluster_energy(1:cluster_num) )
    cluster(i) = list(1)

  end do
!
!  Determine the cluster population counts.
!
  cluster_population(1:cluster_num) = 0

  do i = 1, point_num
    j = cluster(i)
    cluster_population(j) = cluster_population(j) + 1
  end do
!
!  Calculate the mean and sum of squares for each cluster.
!
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
      + point(1:dim_num,i)
  end do

  do i = 1, cluster_num
    if ( 0 < cluster_population(i) ) then
      cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
        / real ( cluster_population(i), kind = 8 )
    end if
  end do
!
!  Set the point energies.
!
  f(1:point_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    f(i) = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
  end do
!
!  Set the cluster energies.
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_energy(j) = cluster_energy(j) + f(i)
  end do
!
!  Adjust the point energies by a weight factor.
!
  do i = 1, point_num
    j = cluster(i)
    if ( 1 < cluster_population(j) ) then
      f(i) = f(i) * real ( cluster_population(j), kind = 8 ) &
        / real ( cluster_population(j) - 1, kind = 8 )
    end if
  end do
!
!  Examine each observation in turn to see if it should be
!  reassigned to a different cluster.
!
  it_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1

    swap = 0

    do i = 1, point_num

      il = cluster(i)
      ir = il

      if ( cluster_population(il) <= 1 ) then
        cycle
      end if

      dc = f(i)

      do j = 1, cluster_num

        if ( j /= il ) then

          de = sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 ) &
            * real ( cluster_population(j), kind = 8 ) &
            / real ( cluster_population(j) + 1, kind = 8 )

          if ( de < dc ) then
            dc = de
            ir = j
          end if

        end if

      end do
!
!  If the lowest value was obtained by staying in the current cluster,
!  then cycle.
!
      if ( ir == il ) then
        cycle
      end if
!
!  Reassign the point from cluster IL to cluster IR.
!
      cluster_center(1:dim_num,il) = &
        ( cluster_center(1:dim_num,il) &
        * real ( cluster_population(il), kind = 8 ) &
        - point(1:dim_num,i) ) / real ( cluster_population(il) - 1, kind = 8 )

      cluster_center(1:dim_num,ir) = &
        ( cluster_center(1:dim_num,ir) &
        * real ( cluster_population(ir), kind = 8 ) &
        + point(1:dim_num,i) ) / real ( cluster_population(ir) + 1, kind = 8 )

      cluster_energy(il) = cluster_energy(il) - f(i)
      cluster_energy(ir) = cluster_energy(ir) + dc
      cluster_population(ir) = cluster_population(ir) + 1
      cluster_population(il) = cluster_population(il) - 1

      cluster(i) = ir
!
!  Adjust the value of F for points in clusters IL and IR.
!
      do j = 1, point_num

        k = cluster(j)

        if ( k == il .or. k == ir ) then

          f(j) = sum ( &
           ( point(1:dim_num,j) - cluster_center(1:dim_num,k) )**2 )

          if ( 1 < cluster_population(k) ) then
            f(j) = f(j) * real ( cluster_population(k), kind = 8 ) &
              / ( real ( cluster_population(k) - 1, kind = 8 ) )
          end if

        end if

      end do

      swap = swap + 1

    end do
!
!  Exit if no reassignments were made during this iteration.
!
    if ( swap == 0 ) then
      exit
    end if

  end do
!
!  Compute the cluster energies.
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    j = cluster(i)

    cluster_energy(j) = cluster_energy(j) + sum ( &
      ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

  end do

  return
end
subroutine kmeans_02 ( dim_num, point_num, cluster_num, it_max, it_num, point, &
  cluster, cluster_center, cluster_population, cluster_energy )

!*****************************************************************************80
!
!! KMEANS_02 applies the K-Means algorithm.
!
!  Discussion:
!
!    The routine attempts to divide POINT_NUM points in 
!    DIM_NUM-dimensional space into CLUSTER_NUM clusters so that the within 
!    cluster sum of squares is minimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2002
!
!  Author:
!
!    Original FORTRAN77 by John Hartigan, Manchek Wong.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hartigan, Manchek Wong,
!    Algorithm AS 136:
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster each 
!    point belongs to.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number 
!    of points in each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the 
!    within-cluster sum of squares.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) an1(cluster_num)
  real ( kind = 8 ) an2(cluster_num)
  integer ( kind = 4 ) cluster(point_num)
  integer ( kind = 4 ) cluster2(point_num)
  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  real ( kind = 8 ) cluster_energy(cluster_num)
  integer ( kind = 4 ) cluster_population(cluster_num)
  real ( kind = 8 ) d(point_num)
  real ( kind = 8 ) db
  real ( kind = 8 ) dt(2)
  integer ( kind = 4 ) i 
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) il
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) itran(cluster_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) live(cluster_num)
  integer ( kind = 4 ) ncp(cluster_num)
  real ( kind = 8 ) point(dim_num,point_num)
  real ( kind = 8 ) temp

  it_num = 0

  if ( cluster_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_02 - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop
  end if

  if ( point_num <= cluster_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_02 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM <= CLUSTER_NUM.'
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_02 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if

  if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_02 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop
  end if
!
!  For each point I, find its two closest centers, CLUSTER(I) and CLUSTER2(I).
!  Assign it to CLUSTER(I).
!
  do i = 1, point_num

    cluster(i) = 1
    cluster2(i) = 2

    do il = 1, 2
      dt(il) = sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,il) )**2 )
    end do

    if ( dt(2) < dt(1) ) then
      cluster(i) = 2
      cluster2(i) = 1
      temp = dt(1)
      dt(1) = dt(2)
      dt(2) = temp
    end if

    do l = 3, cluster_num

      db = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,l) )**2 )

      if ( db < dt(1) ) then
        dt(2) = dt(1)
        cluster2(i) = cluster(i)
        dt(1) = db
        cluster(i) = l
      else if ( db < dt(2) ) then
        dt(2) = db
        cluster2(i) = l
      end if

    end do

  end do
!
!  Update cluster centers to be the average of points contained
!  within them.
!
  cluster_population(1:cluster_num) = 0
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

  do i = 1, point_num
    l = cluster(i)
    cluster_population(l) = cluster_population(l) + 1
    cluster_center(1:dim_num,l) = cluster_center(1:dim_num,l) &
      + point(1:dim_num,i)
  end do
!
!  Check to see if there is any empty cluster.
!
  do l = 1, cluster_num

    if ( cluster_population(l) == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'KMEANS_02 - Fatal error!'
      write ( *, '(a)' ) '  There is an empty cluster.'
      stop
    end if

    cluster_center(1:dim_num,l) = cluster_center(1:dim_num,l) &
      / real ( cluster_population(l), kind = 8 )
!
!  Initialize AN1, AN2, ITRAN and NCP
!  AN1(L) = CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) - 1)
!  AN2(L) = CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) + 1)
!  ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
!           = 0 otherwise
!  In the optimal-transfer stage, NCP(L) stores the step at which
!  cluster L is last updated.
!  In the quick-transfer stage, NCP(L) stores the step at which
!  cluster L is last updated plus POINT_NUM.
!
    an2(l) = real ( cluster_population(l), kind = 8 ) &
      / real ( cluster_population(l) + 1, kind = 8 )

    if ( 1 < cluster_population(l) ) then
      an1(l) = real ( cluster_population(l), kind = 8 ) &
        / real ( cluster_population(l) - 1, kind = 8 )
    else
      an1(l) = huge ( an1(l) )
    end if

    itran(l) = 1
    ncp(l) = -1

  end do

  indx = 0
  ifault = 2
  it_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1
!
!  In this stage, there is only one pass through the data.   Each
!  point is re-allocated, if necessary, to the cluster that will
!  induce the maximum reduction in within-cluster sum of squares.
!
    call kmeans_02_optra ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, cluster2, cluster_population, an1, an2, &
      ncp, d, itran, live, indx )
!
!  Stop if no transfer took place in the last POINT_NUM optimal transfer steps.
!
    if ( indx == point_num ) then
      ifault = 0
      exit
    end if
!
!  Each point is tested in turn to see if it should be re-allocated
!  to the cluster to which it is most likely to be transferred,
!  CLUSTER2(I), from its present cluster, CLUSTER(I).   Loop through the
!  data until no further change is to take place.
!
    call kmeans_02_qtran ( dim_num, point_num, cluster_num, point, &
      cluster_center, cluster, cluster2, cluster_population, an1, an2, &
      ncp, d, itran, indx )
!
!  If there are only two clusters, there is no need to re-enter the
!  optimal transfer stage.
!
    if ( cluster_num == 2 ) then
      ifault = 0
      exit
    end if
!
!  NCP has to be set to 0 before entering OPTRA.
!
    ncp(1:cluster_num) = 0

  end do

  if ( ifault == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_02 - Warning!'
    write ( *, '(a)' ) '  Maximum number of iterations reached'
    write ( *, '(a)' ) '  without convergence.'
  end if
!
!  Compute the within-cluster sum of squares for each cluster.
!
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

  do i = 1, point_num
    cluster_center(1:dim_num,cluster(i)) = &
      cluster_center(1:dim_num,cluster(i)) + point(1:dim_num,i)
  end do

  do j = 1, dim_num
    cluster_center(j,1:cluster_num) = cluster_center(j,1:cluster_num) &
      / real ( cluster_population(1:cluster_num), kind = 8 )
  end do
!
!  Compute the cluster energies.
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    j = cluster(i)

    cluster_energy(j) = cluster_energy(j) + sum ( &
      ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

  end do

  return
end
subroutine kmeans_02_optra ( dim_num, point_num, cluster_num, point, &
  cluster_center, cluster, cluster2, cluster_population, an1, an2, ncp, &
  d, itran, live, indx )

!*****************************************************************************80
!
!! KMEANS_02_OPTRA carries out the optimal transfer stage.
!
!  Discussion:
!
!    Each point is re-allocated, if necessary, to the cluster that
!    will induce a maximum reduction in the within-cluster sum of
!    squares.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 March 2002
!
!  Author:
!
!    Original FORTRAN77 by John Hartigan, Manchek Wong.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hartigan, Manchek Wong,
!    Algorithm AS 136:
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
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
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster 
!    each point belongs to.
!
!    Input/output, integer ( kind = 4 ) CLUSTER2(POINT_NUM), the cluster 
!    to which each point is most likely to be transferred to.
!
!    Input/output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), 
!    the number of points in each cluster.
!
!    Input/output, real ( kind = 8 ) AN1(CLUSTER_NUM), 
!    CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) - 1)
!
!    Input/output, real ( kind = 8 ) AN2(CLUSTER_NUM), 
!    CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) + 1)
!
!    Input/output, integer NCP(CLUSTER_NUM), ?
!
!    Input/output, real ( kind = 8 ) D(POINT_NUM), ?
!
!    Input/output, integer ( kind = 4 ) ITRAN(CLUSTER_NUM), 
!    1 if cluster L is updated in the quick-transfer stage,
!    0 otherwise.  Reset to zero on output.
!
!    Input/output, integer ( kind = 4 ) LIVE(CLUSTER_NUM), ?
!
!    Input/output, integer ( kind = 4 ) INDX, ?
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) al1
  real ( kind = 8 ) al2
  real ( kind = 8 ) alt
  real ( kind = 8 ) alw
  real ( kind = 8 ) an1(cluster_num)
  real ( kind = 8 ) an2(cluster_num)
  integer ( kind = 4 ) cluster(point_num)
  integer ( kind = 4 ) cluster2(point_num)
  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  integer ( kind = 4 ) cluster_population(cluster_num)
  real ( kind = 8 ) d(point_num)
  real ( kind = 8 ) dc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) itran(cluster_num)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) live(cluster_num)
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) ncp(cluster_num)
  real ( kind = 8 ) point(dim_num,point_num)
  real ( kind = 8 ) r2
  real ( kind = 8 ) rr
!
!  If cluster L is updated in the last quick-transfer stage, it
!  belongs to the live set throughout this stage.   Otherwise, at
!  each step, it is not in the live set if it has not been updated
!  in the last POINT_NUM optimal transfer steps.
!
  do l = 1, cluster_num
    if ( itran(l) == 1 ) then
      live(l) = point_num + 1
    end if
  end do

  do i = 1, point_num

    indx = indx + 1
    l1 = cluster(i)
    l2 = cluster2(i)
    ll = l2
!
!  If point I is the only member of cluster L1, no transfer.
!
    if ( 1 < cluster_population(l1) ) then
!
!  If L1 has been updated in this stage, re-compute D(I).
!
      if ( ncp(l1) /= 0 ) then
        d(i) = an1(l1) * sum ( &
          ( point(1:dim_num,i) - cluster_center(1:dim_num,l1) )**2 )
      end if
!
!  Find the cluster with minimum R2.
!
      r2 = an2(l2) * sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,l2) )**2 )

      do l = 1, cluster_num
!
!  If LIVE(L1) <= I, then L1 is not in the live set.   If this is
!  true, we only need to consider clusters that are in the live set
!  for possible transfer of point I.   
!
!  Otherwise, we need to consider all possible clusters.
!
        if ( ( i < live(l1) .or. i < live(l) ) .and. &
             l /= l1 .and. &
             l /= ll ) then

          rr = r2 / an2(l)

          dc = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,l) )**2 )

          if ( dc < rr ) then
            r2 = dc * an2(l)
            l2 = l
          end if

        end if

      end do
!
!  If no transfer is necessary, L2 is the new CLUSTER2(I).
! 
      if ( d(i) <= r2 ) then

        cluster2(i) = l2

      else
!
!  Update cluster centers, LIVE, NCP, AN1 and AN2 for clusters L1 and
!  L2, and update CLUSTER(I) and CLUSTER2(I).
!
        indx = 0
        live(l1) = point_num + i
        live(l2) = point_num + i
        ncp(l1) = i
        ncp(l2) = i
        al1 = cluster_population(l1)
        alw = al1 - 1.0D+00
        al2 = cluster_population(l2)
        alt = al2 + 1.0D+00

        cluster_center(1:dim_num,l1) = ( cluster_center(1:dim_num,l1) * al1 &
          - point(1:dim_num,i) ) / alw

        cluster_center(1:dim_num,l2) = ( cluster_center(1:dim_num,l2) * al2 &
          + point(1:dim_num,i) ) / alt

        cluster_population(l1) = cluster_population(l1) - 1
        cluster_population(l2) = cluster_population(l2) + 1
        an2(l1) = alw / al1

        if ( 1.0D+00 < alw ) then
          an1(l1) = alw / ( alw - 1.0D+00 )
        else
          an1(l1) = huge ( an1(l1) )
        end if

        an1(l2) = alt / al2
        an2(l2) = alt / ( alt + 1.0D+00 )
        cluster(i) = l2
        cluster2(i) = l1

      end if

    end if

    if ( indx == point_num ) then
      return
    end if

  end do
!
!  ITRAN(L) = 0 before entering QTRAN.
!
  itran(1:cluster_num) = 0
!
!  LIVE(L) has to be decreased by POINT_NUM before re-entering OPTRA.
!
  live(1:cluster_num) = live(1:cluster_num) - point_num

  return
end
subroutine kmeans_02_qtran ( dim_num, point_num, cluster_num, point, &
  cluster_center, cluster, cluster2, cluster_population, an1, an2, ncp, &
  d, itran, indx )

!*****************************************************************************80
!
!! KMEANS_02_QTRAN carries out the quick transfer stage.
!
!  Discussion:
!
!    For each point I, CLUSTER(I) and CLUSTER2(I) are switched, if necessary, 
!    to reduce within-cluster sum of squares.  The cluster centers are
!    updated after each step.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2002
!
!  Author:
!
!    Original FORTRAN77 by John Hartigan, Manchek Wong.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hartigan, Manchek Wong,
!    Algorithm AS 136:
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster 
!    each point belongs to.
!
!    Input/output, integer ( kind = 4 ) CLUSTER2(POINT_NUM), the cluster to 
!    which each point is most likely to be transferred to.
!
!    Input/output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), 
!    the number of points in each cluster.
!
!    Input/output, real ( kind = 8 ) AN1(CLUSTER_NUM), 
!    CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) - 1).
!
!    Input/output, real ( kind = 8 ) AN2(CLUSTER_NUM), 
!    CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) + 1).
!
!    Input/output, integer NCP(CLUSTER_NUM), ?
!
!    Input/output, real ( kind = 8 ) D(POINT_NUM), ?
!
!    Input/output, integer ITRAN(CLUSTER_NUM), 
!    1 if cluster L is updated in the quick-transfer stage,
!    0 otherwise.
!
!    Input/output, integer ( kind = 4 ) INDX, is set to 0 if any 
!    updating occurs.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) al1
  real ( kind = 8 ) al2
  real ( kind = 8 ) alt
  real ( kind = 8 ) alw
  real ( kind = 8 ) an1(cluster_num)
  real ( kind = 8 ) an2(cluster_num)
  integer ( kind = 4 ) cluster(point_num)
  integer ( kind = 4 ) cluster2(point_num)
  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  integer ( kind = 4 ) cluster_population(cluster_num)
  integer ( kind = 4 ) count
  real ( kind = 8 ) d(point_num)
  real ( kind = 8 ) dd
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) itran(cluster_num)
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) ncp(cluster_num)
  real ( kind = 8 ) point(dim_num,point_num)
  real ( kind = 8 ) r2
  integer ( kind = 4 ) step
!
!  In the optimal transfer stage, NCP(L) indicates the step at which
!  cluster L is last updated.   In the quick transfer stage, NCP(L)
!  is equal to the step at which cluster L is last updated plus POINT_NUM.
!
  count = 0
  step = 0

  do

    do i = 1, point_num

      count = count + 1
      step = step + 1
      l1 = cluster(i)
      l2 = cluster2(i)
!
!  If point I is the only member of cluster L1, no transfer.
!
      if ( 1 < cluster_population(l1) ) then
!
!  If NCP(L1) < STEP, no need to re-compute distance from point I to
!  cluster L1.   Note that if cluster L1 is last updated exactly POINT_NUM
!  steps ago, we still need to compute the distance from point I to
!  cluster L1.
!
        if ( step <= ncp(l1) ) then
          d(i) = an1(l1) * sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,l1) )**2 )
        end if
!
!  If STEP >= both NCP(L1) and NCP(L2) there will be no transfer of
!  point I at this step.
!
        if ( step < ncp(l1) .or. step < ncp(l2) ) then

          r2 = d(i) / an2(l2)

          dd = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,l2) )**2 )
!
!  Update cluster centers, NCP, CLUSTER_POPULATION, ITRAN, AN1 and AN2 
!  for clusters L1 and L2.   Also update CLUSTER(I) and CLUSTER2(I).   
!
!  Note that if any updating occurs in this stage, INDX is set back to 0.
!
          if ( dd < r2 ) then

            count = 0
            indx = 0
            itran(l1) = 1
            itran(l2) = 1
            ncp(l1) = step + point_num
            ncp(l2) = step + point_num
            al1 = cluster_population(l1)
            alw = al1 - 1.0D+00
            al2 = cluster_population(l2)
            alt = al2 + 1.0D+00

            cluster_center(1:dim_num,l1) = &
              ( cluster_center(1:dim_num,l1) * al1 &
              - point(1:dim_num,i) ) / alw

            cluster_center(1:dim_num,l2) = &
              ( cluster_center(1:dim_num,l2) * al2 &
              + point(1:dim_num,i) ) / alt

            cluster_population(l1) = cluster_population(l1) - 1
            cluster_population(l2) = cluster_population(l2) + 1
            an2(l1) = alw / al1

            if ( 1.0D+00 < alw ) then
              an1(l1) = alw / ( alw - 1.0D+00 )
            else
              an1(l1) = huge ( an1(l1) )
            end if

            an1(l2) = alt / al2
            an2(l2) = alt / ( alt + 1.0D+00 )
            cluster(i) = l2
            cluster2(i) = l1

          end if

        end if

      end if
!
!  If no re-allocation took place in the last POINT_NUM steps, return.
!
      if ( count == point_num ) then
        return
      end if

    end do

  end do

  return
end
subroutine kmeans_03 ( dim_num, point_num, cluster_num, it_max, it_num, point, &
  cluster, cluster_center, cluster_population, cluster_energy )

!*****************************************************************************80
!
!! KMEANS_03 applies the K-Means algorithm.
!
!  Discussion:
!
!    It is possible for a straightforward K-Means algorithm to
!    halt at a non-optimal partition of the points.  This routine
!    tries to improve the input partition if possible.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez, Angel Martinez,
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
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster to which
!    each point belongs.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the 
!    centers associated with the clustering.  On output, these may 
!    have been altered.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number
!    of points in each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy of 
!    the clusters.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) ci
  integer ( kind = 4 ) cj
  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  real ( kind = 8 ), dimension ( cluster_num ) :: cluster_energy
  integer ( kind = 4 ), dimension ( cluster_num ) :: cluster_population
  logical, parameter :: debug = .true.
  real ( kind = 8 ), dimension ( cluster_num ) :: distsq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) list(1)
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  integer ( kind = 4 ) swap
!
!  Check the input.
!
  if ( cluster_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_03 - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_03 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if

  if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_03 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop
  end if

  if ( it_max < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_03 - Fatal error!'
    write ( *, '(a)' ) '  IT_MAX < 0.'
    stop
  end if
!
!  Assign each point to the nearest cluster center.
!
  do i = 1, point_num

    do j = 1, cluster_num
      cluster_energy(j) = sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
    end do

    list = minloc ( cluster_energy(1:cluster_num) )
    cluster(i) = list(1)

  end do
!
!  Determine the cluster populations.
!
  cluster_population(1:cluster_num) = 0
  do i = 1, point_num
    ci = cluster(i)
    cluster_population(ci) = cluster_population(ci) + 1
  end do
!
!  Average the points in each cluster to get a new cluster center.
!
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
      + point(1:dim_num,i)
  end do

  do i = 1, cluster_num
    cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
      / real ( cluster_population(i), kind = 8 )
  end do
!
!  Carry out the iteration.
!
  it_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1

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
            * real ( cluster_population(cj), kind = 8 ) &
            / real ( cluster_population(cj) - 1, kind = 8 )

        else if ( cluster_population(cj) == 0 ) then

          cluster_center(1:dim_num,cj) = point(1:dim_num,i)
          distsq(cj) = 0.0D+00

        else

          distsq(cj) = sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,cj) )**2 ) &
            * real ( cluster_population(cj), kind = 8 ) &
            / real ( cluster_population(cj) + 1, kind = 8 )

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
        ( real ( cluster_population(ci), kind = 8 ) &
        * cluster_center(1:dim_num,ci) &
        - point(1:dim_num,i) ) / real ( cluster_population(ci) - 1, kind = 8 )

      cluster_center(1:dim_num,cj) = &
        ( real ( cluster_population(cj), kind = 8 ) &
        * cluster_center(1:dim_num,cj) &
        + point(1:dim_num,i) ) / real ( cluster_population(cj) + 1, kind = 8 )

      cluster_population(ci) = cluster_population(ci) - 1
      cluster_population(cj) = cluster_population(cj) + 1

      cluster(i) = cj

      swap = swap + 1

    end do
!
!  Exit if no reassignments were made during this iteration.
!
    if ( swap == 0 ) then
      exit
    end if

  end do
!
!  Compute the cluster energies.
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    j = cluster(i)

    cluster_energy(j) = cluster_energy(j) + sum ( &
      ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

  end do

  return
end
subroutine kmeans_w_01 ( dim_num, point_num, cluster_num, it_max, it_num, &
  point, weight, cluster, cluster_center, cluster_population, cluster_energy )

!*****************************************************************************80
!
!! KMEANS_W_01 applies the weighted K-Means algorithm.
!
!  Discussion:
!
!    The input data for the weight K-Means problem includes:
!    * a set of N data points X in M dimensions, 
!    * a set of N nonnegative weights W,
!    * a desired number of clusters K.
!    * an initial set of cluster centers Z,
!    * an (optional) initial set of cluster assignments.
!
!    The goal is to determine K points Z, called cluster centers, and
!    to assign each point X(I) to some cluster Z(J), so that we minimize
!    the weighted standard deviation of the distance of each data point
!    to the center of its cluster.  Writing J = CLUSTER(I) to
!    indicate the index of the nearest cluster center Z(J) to the 
!    point X(I), the quantity we are trying to minimize is the sum
!    of the weighted cluster energies E(J), where:
!
!      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||^2
!
!    Here, we assume that we are using the Euclidean norm, so that
!    
!      || X(I) - Z(J) ||^2 = Sum ( 1 <= K <= M )
!         ( X(I)(K) - Z(J)(K) )^2
!
!    In this notation, X(I)(K) is the K-th spatial component of the
!    I-th data point.
!
!    Note that this routine should give the same results as KMEANS_01
!    in any case in which all the entries of the WEIGHT vector are equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Sparks,
!    Algorithm AS 58: Euclidean Cluster Analysis,
!    Applied Statistics,
!    Volume 22, Number 1, 1973, pages 126-130.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the points.
!
!    Input, real ( kind = 8 ) WEIGHT(POINT_NUM), the weights.
!
!    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), indicates which cluster
!    each point belongs to.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number 
!    of points in each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the 
!    cluster energies.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) c
  integer ( kind = 4 ) cluster(point_num)
  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  real ( kind = 8 ) cluster_energy(cluster_num)
  integer ( kind = 4 ) cluster_population(cluster_num)
  real ( kind = 8 ) cluster_weight(cluster_num)
  real ( kind = 8 ) dc
  real ( kind = 8 ) de
  real ( kind = 8 ) f(point_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) list(1)
  real ( kind = 8 ) point(dim_num,point_num)
  integer ( kind = 4 ) swap
  real ( kind = 8 ) weight(point_num)

  it_num = 0
!
!  Idiot checks.
!
  if ( cluster_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if

  if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop
  end if

  if ( it_max < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  IT_MAX < 0.'
    stop
  end if

  if ( any ( weight(1:point_num) < 0.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  Some weight entry is negative.'
    stop
  end if

  if ( all ( weight(1:point_num) <= 0.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  No weight entry is positive.'
    stop
  end if
!
!  Assign each point to the nearest cluster.
!
  do i = 1, point_num

    do j = 1, cluster_num
      cluster_energy(j) = sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
    end do

    list = minloc ( cluster_energy(1:cluster_num) )
    cluster(i) = list(1)

  end do
!
!  Determine the cluster populations and weights.
!
  cluster_population(1:cluster_num) = 0
  cluster_weight(1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_population(j) = cluster_population(j) + 1
    cluster_weight(j) = cluster_weight(j) + weight(i)
  end do
!
!  Calculate the mean and sum of squares for each cluster.
!
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
      + weight(i) * point(1:dim_num,i)
  end do

  do i = 1, cluster_num
    if ( 0.0D+00 < cluster_weight(i) ) then
      cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
        / cluster_weight(i)
    end if
  end do
!
!  Set the point energies.
!
  f(1:point_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    f(i) = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
  end do
!
!  Set the cluster energies.
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_energy(j) = cluster_energy(j) + weight(i) * f(i)
  end do
!
!  Adjust the point energies by a weight factor.
!
  do i = 1, point_num
    j = cluster(i)
    if ( weight(i) < cluster_weight(j) ) then
      f(i) = f(i) * cluster_weight(j) / ( cluster_weight(j) - weight(i) )
    end if
  end do
!
!  Examine each observation in turn to see if it should be
!  reassigned to a different cluster.
!
  it_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1

    swap = 0

    do i = 1, point_num

      il = cluster(i)
      ir = il

      if ( cluster_population(il) <= 1 ) then
        cycle
      end if

      dc = f(i)

      do j = 1, cluster_num

        if ( j /= il ) then

          de = sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 ) &
            * cluster_weight(j) / ( cluster_weight(j) + weight(i) )

          if ( de < dc ) then
            dc = de
            ir = j
          end if

        end if

      end do
!
!  If the lowest value was obtained by staying in the current cluster,
!  then cycle.
!
      if ( ir == il ) then
        cycle
      end if
!
!  Reassign the point from cluster IL to cluster IR.
!
      cluster_center(1:dim_num,il) = &
        ( cluster_weight(il) * cluster_center(1:dim_num,il) &
        - weight(i) * point(1:dim_num,i) ) &
        / ( cluster_weight(il) - weight(i) )

      cluster_center(1:dim_num,ir) = &
        ( cluster_weight(ir) * cluster_center(1:dim_num,ir) &
        + weight(i) * point(1:dim_num,i) ) &
        / ( cluster_weight(ir) + weight(i) )

      cluster_weight(il) = cluster_weight(il) - weight(i)
      cluster_weight(ir) = cluster_weight(ir) + weight(i)

      cluster_energy(il) = cluster_energy(il) - f(i)
      cluster_energy(ir) = cluster_energy(ir) + dc

      cluster_population(ir) = cluster_population(ir) + 1
      cluster_population(il) = cluster_population(il) - 1

      cluster(i) = ir
!
!  Adjust the value of F for all points in clusters IL and IR.
!
      do j = 1, point_num

        k = cluster(j)

        if ( k == il .or. k == ir ) then

          f(j) = sum ( &
           ( point(1:dim_num,j) - cluster_center(1:dim_num,k) )**2 )

          if ( weight(j) < cluster_weight(k) ) then
            f(j) = f(j) * cluster_weight(k) &
              / ( cluster_weight(k) - weight(j) )
          end if

        end if

      end do

      swap = swap + 1

    end do
!
!  Exit if no reassignments were made during this iteration.
!
    if ( swap == 0 ) then
      exit
    end if

  end do
!
!  Compute the energy based on the final value of the cluster centers.
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    c = cluster(i)

    cluster_energy(c) = cluster_energy(c) + weight(i) * sum ( &
      ( point(1:dim_num,i) - cluster_center(1:dim_num,c) )**2 )

  end do

  return
end
subroutine kmeans_w_03 ( dim_num, point_num, cluster_num, it_max, it_num, &
  point, weight, cluster, cluster_center, cluster_population, cluster_energy )

!*****************************************************************************80
!
!! KMEANS_W_03 applies the weighted K-Means algorithm.
!
!  Discussion:
!
!    The input data for the weight K-Means problem includes:
!    * a set of N data points X in M dimensions, 
!    * a set of N nonnegative weights W,
!    * a desired number of clusters K.
!    * an initial set of cluster centers Z,
!    * an (optional) initial set of cluster assignments.
!
!    The goal is to determine K points Z, called cluster centers, and
!    to assign each point X(I) to some cluster Z(J), so that we minimize
!    the weighted standard deviation of the distance of each data point
!    to the center of its cluster.  Writing J = CLUSTER(I) to
!    indicate the index of the nearest cluster center Z(J) to the 
!    point X(I), the quantity we are trying to minimize is the sum
!    of the weighted cluster energies E(J), where:
!
!      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||**2
!
!    Here, we assume that we are using the Euclidean norm, so that
!    
!      || X(I) - Z(J) ||**2 = Sum ( 1 <= K <= M )
!        ( X(I)(K) - Z(J)(K) )**2
!
!    In this notation, X(I)(K) is the K-th spatial component of the
!    I-th data point.
!
!    Note that this routine should give the same results as KMEANS_01
!    in any case in which all the entries of the WEIGHT vector are equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez, Angel Martinez,
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
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, real ( kind = 8 ) WEIGHT(POINT_NUM), the weights.
!
!    Input/output, integer ( kind = 4 ) CLUSTER(POINT_NUM), the cluster 
!    to which each point belongs.  On output, these may have been altered.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!    centers associated with the clustering.  On output, these may
!    have been altered.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number
!    of points in each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the energy of
!    the clusters.
!
  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) ci
  integer ( kind = 4 ) cj
  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ), dimension ( dim_num, cluster_num ) :: cluster_center
  real ( kind = 8 ), dimension ( cluster_num ) :: cluster_energy
  integer ( kind = 4 ), dimension ( cluster_num ) :: cluster_population
  real ( kind = 8 ), dimension ( cluster_num ) :: cluster_weight
  logical, parameter :: debug = .true.
  real ( kind = 8 ), dimension ( cluster_num ) :: distsq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) list(1)
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  integer ( kind = 4 ) swap
  real ( kind = 8 ) weight(point_num)
!
!  Check the input.
!
  if ( cluster_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_03 - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_03 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if

  if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_03 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop
  end if

  if ( it_max < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_03 - Fatal error!'
    write ( *, '(a)' ) '  IT_MAX < 0.'
    stop
  end if

  if ( any ( weight(1:point_num) < 0.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_03 - Fatal error!'
    write ( *, '(a)' ) '  Some weight entry is negative.'
    stop
  end if

  if ( all ( weight(1:point_num) <= 0.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_03 - Fatal error!'
    write ( *, '(a)' ) '  No weight entry is positive.'
    stop
  end if
!
!  Assign each observation to the nearest cluster center.
!
  do i = 1, point_num

    do j = 1, cluster_num
      cluster_energy(j) = sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
    end do

    list = minloc ( cluster_energy(1:cluster_num) )
    cluster(i) = list(1)

  end do
!
!  Determine the cluster populations and weights.
!
  cluster_population(1:cluster_num) = 0
  cluster_weight(1:cluster_num) = 0.0D+00

  do i = 1, point_num
    ci = cluster(i)
    cluster_population(ci) = cluster_population(ci) + 1
    cluster_weight(ci) = cluster_weight(ci) + weight(i)
  end do
!
!  Average the points in each cluster to get a new cluster center.
!
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
      + weight(i) * point(1:dim_num,i)
  end do

  do i = 1, cluster_num
    if ( 0.0D+00 < cluster_weight(i) ) then
      cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
        / cluster_weight(i)
    end if
  end do
!
!  Carry out the iteration.
!
  it_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1

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
            * cluster_weight(cj)  &
            / ( cluster_weight(cj) - weight(i) )

        else if ( cluster_population(cj) == 0 ) then

          cluster_center(1:dim_num,cj) = point(1:dim_num,i)
          distsq(cj) = 0.0D+00

        else

          distsq(cj) = sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,cj) )**2 ) &
            * cluster_weight(cj) &
            / ( cluster_weight(cj) + weight(i) )

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
        ( cluster_weight(ci) * cluster_center(1:dim_num,ci) &
        - weight(i) * point(1:dim_num,i) ) &
        / ( cluster_weight(ci) - weight(i) )

      cluster_center(1:dim_num,cj) = &
        ( cluster_weight(cj) * cluster_center(1:dim_num,cj) &
        + weight(i) * point(1:dim_num,i) ) &
        / ( cluster_weight(cj) + weight(i) )

      cluster_population(ci) = cluster_population(ci) - 1
      cluster_population(cj) = cluster_population(cj) + 1

      cluster_weight(ci) = cluster_weight(ci) - weight(i)
      cluster_weight(cj) = cluster_weight(cj) + weight(i)

      cluster(i) = cj

      swap = swap + 1

    end do
!
!  Exit if no reassignments were made during this iteration.
!
    if ( swap == 0 ) then
      exit
    end if

  end do
!
!  Compute the energy based on the final value of the cluster centers.
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    ci = cluster(i)

    cluster_energy(ci) = cluster_energy(ci) + weight(i) * sum ( &
      ( point(1:dim_num,i) - cluster_center(1:dim_num,ci) )**2 )

  end do

  return
end
function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r4_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
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
!    11 August 2004
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
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, 
!    which should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r4_uniform_01

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

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
!    Output, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 )   m
  integer   ( kind = 4 )   n

  integer   ( kind = 4 )   ierror
  character ( len = * )    input_filename
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  integer   ( kind = 4 )   j
  character ( len = 255 )  line
  real      ( kind = 8 )   table(m,n)
  real      ( kind = 8 )   x(m)

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

  character ( len = * )  input_filename
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

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
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
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
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
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
!  For greater precision in the output file, try:
!
!                                            '(', m, 'g', 24, '.', 16, ')'
!
  write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
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
!    Input, integer ( kind = 4 ) N, the number of entries in 
!    the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

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
!    19 December 2001
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
!    use the input value of SEED to initialize the random number generator,
!    and SEED is not changed on output.
!
  implicit none

  integer ( kind = 4 ) count
  integer ( kind = 4 ) count_max
  integer ( kind = 4 ) count_rate
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
  real ( kind = 8 ) t
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

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, user SEED = ', seed
    end if

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, arbitrary SEED = ', &
        seed
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
!
!    0, no errors occurred.
!
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
