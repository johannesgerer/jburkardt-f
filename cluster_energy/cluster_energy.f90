program main

!*****************************************************************************80
!
!! MAIN is the main program for CLUSTER_ENERGY.
!
!  Discussion:
!
!    CLUSTER_ENERGY computes the minimal cluster energy for a set of data.
!
!    The current code only generates the data internally, and does
!    not read an input file of data, which would be more useful.
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
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) c_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: center
  integer ( kind = 4 ) c_hi
  integer ( kind = 4 ) c_lo
  real ( kind = 8 ), allocatable, dimension ( : ) :: e
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) it_hi
  integer ( kind = 4 ) it_lo
  real ( kind = 8 ), allocatable, dimension ( :, :) :: point
  integer ( kind = 4 ) point_dist
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: r_max
  real ( kind = 8 ), allocatable, dimension ( : ) :: r_min
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CLUSTER_ENERGY'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Compute the cluster energy of a set of points.'
  write ( *, '(a)' ) ' '

  seed = 123456789

  call random_initialize ( seed )

  call i4_input ( 'Enter spatial dimension DIM_NUM', dim_num, ierror )

  allocate ( r_min(1:dim_num) )
  allocate ( r_max(1:dim_num) )

  call input_get ( point_dist, dim_num, r_min, r_max, point_num, c_lo, &
    c_hi, it_lo, it_hi )

  call input_check ( point_dist, dim_num, r_min, r_max, point_num, c_lo, &
    c_hi, it_lo, it_hi  )

  call input_print ( point_dist, dim_num, r_min, r_max, point_num, c_lo, &
    c_hi, it_lo, it_hi  )

  allocate ( point(1:dim_num,1:point_num) )

  call point_generate ( point_dist, dim_num, r_min, r_max, point_num, point )

  if ( .false. ) then
    call point_print ( dim_num, point_num, point )
  end if

  allocate ( e(c_lo:c_hi) )
  allocate ( center(1:dim_num,1:c_hi) )

  do c_num = c_lo, c_hi

    call cluster_iteration ( dim_num, point_num, c_num, point, r_min, &
      r_max, seed, center, e(c_num), it_lo, it_hi  )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Energy table:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Number    Total        Energy         Square'
  write ( *, '(a)' ) '    of      Energy        per            Root'
  write ( *, '(a)' ) '  Clusters               point'
  write ( *, '(a)' ) ' '

  do c_num = c_lo, c_hi
    t = e(c_num) / real ( point_num, kind = 8 )
    write ( *, '(2x,i4,2x,g14.6,g14.6,g14.6)' ) c_num, e(c_num), t, sqrt ( t )
  end do
!
!  Free memory.
!
  deallocate ( center )
  deallocate ( e )
  deallocate ( point )
  deallocate ( r_max )
  deallocate ( r_min )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CLUSTER_ENERGY'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine cluster_iteration ( dim_num, point_num, c_num, point, r_min, &
  r_max, seed, center, e, it_lo, it_hi  )

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
!    20 February 2005
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
!    Input, integer ( kind = 4 ) C_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) R_MIN(DIM_NUM), R_MAX(DIM_NUM), the coordinates
!    of the minimum and maximum corners of the region.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) CENTER(DIM_NUM,C_NUM), the centers associated
!    with the minimal energy clustering.
!
!    Output, real ( kind = 8 ) E, the energy associated with the minimal 
!    energy clustering.
!
!    Input, integer ( kind = 4 ) IT_LO, IT_HI, the least and most number of energy
!    iterations.
!
  implicit none

  integer ( kind = 4 ) c_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ), dimension ( dim_num, c_num ) :: center
  real ( kind = 8 ), dimension ( dim_num, c_num ) :: centroid
  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ) e
  real ( kind = 8 ) e_new
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it_hi
  integer ( kind = 4 ) it_lo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) missed
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point
  real ( kind = 8 ), dimension ( dim_num ) :: r
  real ( kind = 8 ), dimension ( dim_num ) :: r_max
  real ( kind = 8 ), dimension ( dim_num ) :: r_min
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), dimension ( c_num ) :: subs

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of clusters = ', c_num
!
!  Initialize the centers randomly.
!
  do j = 1, c_num
    call random_number ( harvest = r(1:dim_num) )
    center(1:dim_num,j) = ( 1.0D+00 - r(1:dim_num) ) * r_min(1:dim_num) &
                                    + r(1:dim_num)   * r_max(1:dim_num)
  end do
!
!  Initialize the clusters randomly.
!
  subs(1:c_num) = 0
  do i = 1, point_num
    cluster(i) = i4_uniform ( 1, c_num, seed )
    j = cluster(i)
    subs(j) = subs(j) + 1
  end do

  call energy ( dim_num, point_num, c_num, point, center, cluster, e )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Initial energy =                        ', e

  do j = 1, c_num
    write ( *, '(2i5,6f10.4)' ) j, subs(j), center(1:dim_num,j)
  end do

  do it = 1, it_hi
!
!  #1:
!  Assign each point to the cluster of its nearest center.
!
    do i = 1, point_num
      call nearest_point ( dim_num, c_num, center, point(1,i), j )
      cluster(i) = j
    end do

    call energy ( dim_num, point_num, c_num, point, center, cluster, e )

    write ( *, '(a,g14.6)' ) '  Energy with old centers, new clusters = ', e
!
!  #2:
!  Determine the centroids of the clusters.
!
    centroid(1:dim_num,1:c_num) = 0.0D+00
    subs(1:c_num) = 0

    do i = 1, point_num
      j = cluster(i)
      subs(j) = subs(j) + 1
      centroid(1:dim_num,j) = centroid(1:dim_num,j) + point(1:dim_num,i)
    end do

    missed = 0

    do j = 1, c_num

      if ( subs(j) /= 0 ) then
        centroid(1:dim_num,j) = centroid(1:dim_num,j) &
          / real ( subs(j), kind = 8 )
      else
        missed = missed + 1
        call random_number ( harvest = r(1:dim_num) )
        centroid(1:dim_num,j) = ( 1.0D+00 - r(1:dim_num) ) * r_min(1:dim_num) &
                                          + r(1:dim_num)   * r_max(1:dim_num)

      end if

    end do

    if ( missed /= 0 ) then
      write ( *, '(a,i8)' ) 'Number of empty clusters was ', missed
    end if

    call energy ( dim_num, point_num, c_num, point, centroid, cluster, e_new )

    write ( *, '(a,g14.6)' ) '  Energy with new centers, new clusters = ', e_new

    do j = 1, c_num
      write ( *, '(2i8,6f10.4)' ) j, subs(j), center(1:dim_num,j)
      write ( *, '(16x,6f10.4)' )              centroid(1:dim_num,j)
    end do
!
!  Update the centers.
!
    center(1:dim_num,1:c_num) = centroid(1:dim_num,1:c_num)

  end do

  return
end
subroutine energy ( dim_num, point_num, c_num, point, center, cluster, e )

!*****************************************************************************80
!
!! ENERGY computes the energy of a given clustering.
!
!  Discussion:
!
!    The energy of the clustering is the sum of the energy of each cluster.
!
!    The sum of a cluster is the sum of the squares of the distances of 
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
!    Input, integer ( kind = 4 ) C_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, real ( kind = 8 ) CENTER(DIM_NUM,C_NUM), the center points.
!
!    Input, integer ( kind = 4 ) CLUSTER(POINT_NUM), indicates the cluster to which
!    each data point belongs.
!
!    Output, real ( kind = 8 ) E, the energy of the clustering.
!
  implicit none

  integer ( kind = 4 ) c_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ), dimension ( dim_num, c_num ) :: center
  integer ( kind = 4 ), dimension ( point_num ) :: cluster
  real ( kind = 8 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: point

  e = 0.0D+00
  do i = 1, point_num
    j = cluster(i)
    e = e + sum ( ( center(1:dim_num,j) - point(1:dim_num,i) )**2 )
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
!    Output, integer ( kind = 4 ) VALUE, the value input by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = * ) string
  integer ( kind = 4 ) value

  ierror = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  read ( *, *, iostat = ierror ) value

  return
end
subroutine i4_range_input ( string, value1, value2, ierror )

!*****************************************************************************80
!
!! I4_RANGE_INPUT reads a pair of integers from the user, representing a range.
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
!    Output, integer ( kind = 4 ) VALUE1, VALUE2, the values entered by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = * ) string
  integer ( kind = 4 ) value1
  integer ( kind = 4 ) value2

  ierror = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )
  read ( *, *, iostat = ierror ) value1, value2

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
subroutine input_check ( point_dist, dim_num, r_min, r_max, point_num, &
  c_lo, c_hi, it_lo, it_hi )

!*****************************************************************************80
!
!! INPUT_CHECK performs some simple checks on the problem data.
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
!    Input, real ( kind = 8 ) R_MIN(DIM_NUM), R_MAX(DIM_NUM), the coordinates
!    of the minimum and maximum corners of the region.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) C_LO, C_HI, the minimum and maximum number of
!    clusters to consider.
!
!    Input, integer ( kind = 4 ) IT_LO, IT_HI, the least and most number of energy
!    iterations.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) c_hi
  integer ( kind = 4 ) c_lo
  integer ( kind = 4 ) it_hi
  integer ( kind = 4 ) it_lo
  integer ( kind = 4 ) point_dist
  integer ( kind = 4 ) point_num
  real ( kind = 8 ) r_max(dim_num)
  real ( kind = 8 ) r_min(dim_num)

  if ( point_dist < 1 .or. 2 < point_dist ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INPUT_CHECK - Fatal error!'
    write ( *, '(a)' ) &
      '  Point distribution option less than 1 or greater than 2.'
    write ( *, '(a,i8)' ) '  POINT_DIST = ', point_dist
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INPUT_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Number of dimensions DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( any ( r_max(1:dim_num) <= r_min(1:dim_num) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INPUT_CHECK - Fatal error!'
    write ( *, '(a)' ) '  For at least one index I, R_MAX(I) <= R_MIN(I).'
    stop
  end if

  if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INPUT_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Number of points POINT_NUM < 1.'
    write ( *, '(a,i8)' ) '  POINT_NUM = ', point_num
    stop
  end if

  if ( c_hi < c_lo ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INPUT_CHECK - Fatal error!'
    write ( *, '(a)' ) '  CHI_NUM < CLO_NUM.'
    write ( *, '(a,i8)' ) '  CLO_NUM = ', c_lo
    write ( *, '(a,i8)' ) '  CHI_NUM = ', c_hi
    stop
  end if

  if ( point_num < c_hi ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INPUT_CHECK - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < CHI_NUM.'
    write ( *, '(a,i8)' ) '  CHI_NUM = ', c_hi
    write ( *, '(a,i8)' ) '  POINT_NUM = ', point_num
    stop
  end if

  if ( it_lo < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INPUT_CHECK - Fatal error!'
    write ( *, '(a)' ) '  IT_LO < 1.'
    write ( *, '(a,i8)' ) '  IT_LO = ', it_lo
    stop
  end if

  if ( it_hi < it_lo ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INPUT_CHECK - Fatal error!'
    write ( *, '(a)' ) '  IT_HI < IT_LO.'
    write ( *, '(a,i8)' ) '  IT_LO = ', it_lo
    write ( *, '(a,i8)' ) '  IT_HI = ', it_hi
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INPUT_CHECK - The input data passed all checks.'

  return
end
subroutine input_get ( point_dist, dim_num, r_min, r_max, point_num, &
  c_lo, c_hi, it_lo, it_hi )

!*****************************************************************************80
!
!! INPUT_GET gets the value of certain input quantities from the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) POINT_DIST, the point distribution to use.
!    1, use a uniform random distribution.
!    2, use a uniform grid of points.  (This hasn't been set up properly
!       except for 1D!).
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Output, real ( kind = 8 ) R_MIN(DIM_NUM), R_MAX(DIM_NUM), the coordinates
!    of the minimum and maximum corners of the region.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Output, integer ( kind = 4 ) C_LO, C_HI, the minimum and maximum number of
!    clusters to consider.
!
!    Output, integer ( kind = 4 ) IT_LO, IT_HI, the least and most number of energy
!    iterations.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) c_hi
  integer ( kind = 4 ) c_lo
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) it_hi
  integer ( kind = 4 ) it_lo
  integer ( kind = 4 ) point_dist
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), dimension ( dim_num ) :: r_max
  real ( kind = 8 ), dimension ( dim_num ) :: r_min

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The region is assumed to be an interval [A,B], or'
  write ( *, '(a)' ) '  a rectangle [A(1),B(1)] x [A(2),B(2)], or'
  write ( *, '(a)' ) '  a generalized box, described by lower and upper bounds.'

  call r8vec_range_input ( &
    'Enter all the lower bounds for the region, then the upper bounds:', &
    dim_num, r_min, r_max, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program will place data points in the region.'

  call i4_input ( 'Enter the number of data points POINT_NUM', &
    point_num, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data points can be randomly or regularly placed.'

  call i4_input ( 'Enter 1 for uniform random data, 2 for uniform grid', &
    point_dist, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For a given number of clusters, the program seeks the'
  write ( *, '(a)' ) '  minimal energy.  Usually, we do this calculation for a'
  write ( *, '(a)' ) '  range of clusters, and tabulate their energy.'

  call i4_range_input ( 'Enter the lower and upper number of clusters', c_lo, &
    c_hi, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  To seek the minimal energy, the program carries out an iteration.'
  write ( *, '(a)' ) &
    '  Usually, we specify lower and upper limits for the number of iterations.'

  call i4_range_input ( &
    'Enter the lower and upper number of energy iterations', &
    it_lo, it_hi, ierror )

  return
end
subroutine input_print ( point_dist, dim_num, r_min, r_max, point_num, &
  c_lo, c_hi, it_lo, it_hi )

!*****************************************************************************80
!
!! INPUT_PRINT prints the value of the input quantities.
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
!    Input, real ( kind = 8 ) R_MIN(DIM_NUM), R_MAX(DIM_NUM), the coordinates
!    of the minimum and maximum corners of the region.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) C_LO, C_HI, the minimum and maximum number of
!    clusters to consider.
!
!    Input, integer ( kind = 4 ) IT_LO, IT_HI, the least and most number of energy
!    iterations.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) c_hi
  integer ( kind = 4 ) c_lo
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_hi
  integer ( kind = 4 ) it_lo
  integer ( kind = 4 ) point_dist
  integer ( kind = 4 ) point_num
  real ( kind = 8 ) r_max(dim_num)
  real ( kind = 8 ) r_min(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INPUT_PRINT:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of spatial dimensions = ', dim_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lower and upper region bounds:'
  write ( *, '(a)' ) ' '
  do i = 1, dim_num
    write ( *, '(4x,2f10.6)' ) r_min(i), r_max(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of data points =              ', point_num
  write ( *, '(a,i8)' ) '  Lowest number of clusters =  	', c_lo
  write ( *, '(a,i8)' ) '  Highest number of clusters =         ', c_hi
  write ( *, '(a,i8)' ) '  Least number of energy iterations =  ', it_lo
  write ( *, '(a,i8)' ) '  Most number of energy iterations =   ', it_hi
  write ( *, '(a)' ) ' '
  if ( point_dist == 1 ) then
    write ( *, '(a)' ) '  Data points will be uniform random values.'
  else if ( point_dist == 2 ) then
    write ( *, '(a)' ) '  Data points will be a uniform grid.'
  end if

  return
end
subroutine nearest_point ( dim_num, c_num, center, point, nearest )

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
!    Input, integer ( kind = 4 ) C_NUM, the number of center points.
!
!    Input, real ( kind = 8 ) CENTER(DIM_NUM,C_NUM), the coordinates of the center points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM), the data point to be checked.
!
!    Output, integer ( kind = 4 ) NEAREST, the index of the center point closest to
!    the data point.
!
  implicit none

  integer ( kind = 4 ) c_num
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_new
  real ( kind = 8 ), dimension ( dim_num, c_num ) :: center
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nearest
  real ( kind = 8 ), dimension ( dim_num ) :: point

  dist = huge ( dist )
  nearest = 0

  do j = 1, c_num

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
!    Input, real ( kind = 8 ) R_MIN(DIM_NUM), R_MAX(DIM_NUM), the coordinates of the
!    minimum and maximum corners of the region.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points to generate.
!
!    Output, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates of
!    the points.
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
        r = real ( i - 1, kind = 8 ) / real ( point_num - 1, kind = 8 )
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
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates of
!    the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) i
  real ( kind = 8 ) point(dim_num,point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Data points:'
  write ( *, '(a)' ) ' '

  do i = 1, point_num
    write ( *, '(i8,6f10.4)' ) i, point(1:dim_num,i)
  end do

  return
end
subroutine r8vec_range_input ( string, dim_num, value1, value2, ierror )

!*****************************************************************************80
!
!! R8VEC_RANGE_INPUT reads two real vectors from the user, representing a range.
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
