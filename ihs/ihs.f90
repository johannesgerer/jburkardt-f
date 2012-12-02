subroutine covariance ( dim_num, n, x, average, std, covc )

!*****************************************************************************80
!
!! COVARIANCE does a covariance calculation for IHS solutions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points to be generated.
!
!    Input, integer ( kind = 4 ) X(DIM_NUM,N), the points.
!
!    Output, real ( kind = 8 ) AVERAGE, the average minimum distance.
!
!    Output, real ( kind = 8 ) STD, the standard deviation of the
!    minimum distances.
!
!    Output, real ( kind = 8 ) COVC, the covariance of the minimum distances.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) average
  real ( kind = 8 ) covc
  real ( kind = 8 ) dist
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) mindist(n)
  real ( kind = 8 ), parameter :: r8_huge = 1.0D+30
  real ( kind = 8 ) std
  real ( kind = 8 ) vec(dim_num)
  integer ( kind = 4 ) x(dim_num,n)
!
!  Set up the distance matrix.
!
  do i = 1, n
    mindist(i) = r8_huge
    do j = 1, n
      if ( i /= j ) then
        vec(1:dim_num) = real ( x(1:dim_num,i) - x(1:dim_num,j), kind = 8 )
        dist = sqrt ( dot_product ( vec(1:dim_num), vec(1:dim_num) ) )
        mindist(i) = min ( mindist(i), dist )
      end if
    end do
  end do
!
!  Find the average minimum distance.
!
  average = sum ( mindist(1:n) ) / real ( n, kind = 8 )
!
!  Compute the standard deviation of the distances.
!
  call r8vec_std ( n, mindist, std )
!
!  Compute the covariance.
!
  covc = std / average

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

  integer   ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer   ( kind = 4 ) seed
  real ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8  ) today
  integer   ( kind = 4 ) values(8)
  character ( len = 5  ) zone

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

  seed = int ( real ( i4_huge, kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == i4_huge ) then
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
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
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
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
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
    seed = seed + i4_huge
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
subroutine i4vec_uniform ( n, a, b, seed, x )

!*****************************************************************************80
!
!! I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer ( kind = 4 ) values.
!
!    The pseudorandom numbers should be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) X(N), a vector of numbers between A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
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

    x(i) = value

  end do

  return
end
subroutine ihs ( dim_num, n, duplication, seed, x )

!*****************************************************************************80
!
!! IHS implements the improved distributed hypercube sampling algorithm.
!
!  Discussion:
!
!    N Points in an DIM_NUM dimensional Latin hypercube are to be selected.
!    Each of the coordinate dimensions is discretized to the values
!    1 through N.  The points are to be chosen in such a way that
!    no two points have any coordinate value in common.  This is
!    a standard Latin hypercube requirement, and there are many
!    solutions.
!
!    This algorithm differs in that it tries to pick a solution
!    which has the property that the points are "spread out"
!    as evenly as possible.  It does this by determining an optimal
!    even spacing, and using the DUPLICATION factor to allow it
!    to choose the best of the various options available to it.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Brian Beachkofski, Ramana Grandhi,
!    Improved Distributed Hypercube Sampling,
!    American Institute of Aeronautics and Astronautics Paper 2002-1274.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points to be generated.
!
!    Input, integer ( kind = 4 ) DUPLICATION, the duplication factor.  This must
!    be at least 1.  A value of 5 is reasonable.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) X(DIM_NUM,N), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) duplication
  integer ( kind = 4 ) n

  integer ( kind = 4 ) avail(dim_num,n)
  integer ( kind = 4 ) best
  integer ( kind = 4 ) count
  real ( kind = 8 ) dist
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) list(duplication*n)
  real ( kind = 8 ) min_all
  real ( kind = 8 ) min_can
  real ( kind = 8 ) opt
  integer ( kind = 4 ) point(dim_num,duplication*n)
  integer ( kind = 4 ) point_index
  real ( kind = 8 ), parameter :: r8_huge = 1.0D+30
  integer ( kind = 4 ) seed
  real ( kind = 8 ) vec(dim_num)
  integer ( kind = 4 ) x(dim_num,n)

  opt = real ( n, kind = 8 ) / &
    ( real ( n, kind = 8 ) )**( 1.0D+00 / real ( dim_num, kind = 8 ) )
!
!  Pick the first point.
!
  call i4vec_uniform ( dim_num, 1, n, seed, x(1:dim_num,n) )
!
!  Initialize AVAIL,
!  and set an entry in a random row of each column of AVAIL to N.
!
  do j = 1, n
    avail(1:dim_num,j) = j
  end do

  do i = 1, dim_num
    avail(i,x(i,n)) = n
  end do
!
!  Main loop:
!  Assign a value to X(1:DIM_NUM,COUNT) for COUNT = N-1 down to 2.
!
  do count = n-1, 2, -1
!
!  Generate valid points.
!
    do i = 1, dim_num

      do k = 1, duplication
        list(count*(k-1)+1:k*count) = avail(i,1:count)
      end do

      do k = count*duplication, 1, -1
        point_index = i4_uniform ( 1, k, seed )
        point(i,k) = list(point_index)
        list(point_index) = list(k)
      end do

    end do
!
!  For each candidate, determine the distance to all the
!  points that have already been selected, and save the minimum value.
!
    min_all = r8_huge
    best = 0

    do k = 1, duplication*count

      min_can = r8_huge

      do j = count+1, n
        vec(1:dim_num) = real ( point(1:dim_num,k) - x(1:dim_num,j), kind = 8 )
        dist = sqrt ( dot_product ( vec(1:dim_num), vec(1:dim_num) ) )
        min_can = min ( min_can, dist )
      end do

      if ( abs ( min_can - opt ) < min_all ) then
        min_all = abs ( min_can - opt )
        best = k
      end if

    end do

    x(1:dim_num,count) = point(1:dim_num,best)
!
!  Having chosen X(*,COUNT), update AVAIL.
!
    do i = 1, dim_num

      do j = 1, n
        if ( avail(i,j) == x(i,count) ) then
          avail(i,j) = avail(i,count)
        end if
      end do

    end do

  end do
!
!  For the last point, there's only one choice.
!
  x(1:dim_num,1) = avail(1:dim_num,1)

  return
end
subroutine ihs_write ( dim_num, n, d, seed_init, seed, x, file_out_name )

!*****************************************************************************80
!
!! IHS_WRITE writes an IHS dataset to a file.
!
!  Discussion:
!
!    The initial lines of the file are comments, which begin with a
!    '#' character.
!
!    Thereafter, each line of the file contains the M-dimensional
!    components of the next entry of the sequence.
!
!    Note that the actual values of the data are integers between 1
!    and N.  For our convenience, these are rescaled by the
!    mapping
!
!      I -> ( 2 * I - 1 )/ ( 2 * N ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 June 2004
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
!    Input, integer ( kind = 4 ) D, the duplication factor.
!
!    Input, integer ( kind = 4 ) SEED_INIT, the initial random number seed.
!
!    Input, integer ( kind = 4 ) SEED, the current random number seed.
!
!    Input, integer ( kind = 4 ) X(DIM_NUM,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of
!    the output file.
!
  implicit none

  integer   ( kind = 4 ) dim_num
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) d
  character ( len = * ) file_out_name
  integer   ( kind = 4 ) file_out_unit
  integer   ( kind = 4 ) ios
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) seed
  integer   ( kind = 4 ) seed_init
  character ( len = 80 ) string
  integer   ( kind = 4 ) x(dim_num,n)

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IHS_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  write ( file_out_unit, '(a)'       ) '#  ' // trim ( file_out_name )
  write ( file_out_unit, '(a)'       ) '#  created by IHS_DATASET.F90.'
  write ( file_out_unit, '(a)'       ) '#'
  write ( file_out_unit, '(a)'       ) '#'
  write ( file_out_unit, '(a,i8)'    ) &
    '#  Spatial dimension DIM_NUM = ', dim_num
  write ( file_out_unit, '(a,i8)'    ) '#  Number of points N =  ', n
  write ( file_out_unit, '(a,g14.6)' ) '#  EPSILON (unit roundoff) = ', &
    epsilon ( 1.0D+00 )
  write ( file_out_unit, '(a,i8)'    ) '#  Duplication factor D = ', d
  write ( file_out_unit, '(a,i12)'   ) '#  Initial SEED_INIT = ', seed_init
  write ( file_out_unit, '(a,i12)'   ) '#  Current SEED =      ', seed
  write ( file_out_unit, '(a)'       ) '#'

  write ( string, '(a,i3,a)' ) '(', dim_num, '(2x,f10.6))'
  do j = 1, n
    write ( file_out_unit, string ) &
      real ( 2 * x(1:dim_num,j) - 1 ) / real ( 2 * n )
  end do

  close ( unit = file_out_unit )

  return
end
subroutine r8vec_std ( n, a, std )

!*****************************************************************************80
!
!! R8VEC_STD returns the standard deviation of a real vector.
!
!  Discussion:
!
!    The standard deviation of a vector X of length N is defined as
!
!      mean ( X(1:n) ) = sum ( X(1:n) ) / n
!
!      std ( X(1:n) ) = sqrt ( sum ( ( X(1:n) - mean )**2 ) / ( n - 1 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!    N should be at least 2.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, real ( kind = 8 ) STD, the standard deviation of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) mean
  real ( kind = 8 ) std

  if ( n < 2 ) then

    std = 0.0D+00

  else

    mean = sum ( a(1:n) ) / real ( n, kind = 8 )

    std = sum ( ( a(1:n) - mean )**2 )

    std = sqrt ( std / real ( n - 1, kind = 8 ) )

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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  character ( len = 8  ) date
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  character ( len = 10 )  time
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y
  character ( len = 5  ) zone

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
