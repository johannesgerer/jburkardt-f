program main

!*****************************************************************************80
!
!! MAIN is the main program for BALL_VOLUME_MONTE_CARLO.
!
!  Discussion:
!
!    DIM_NUM = 6 is a reasonable test.
!
!    N_LOG2_MAX = 25 puts a strain on the system, since we generate that
!    many temporary points at once.  To solve bigger problems, it would
!    be better to compute the new points in batches whose maximum size
!    is limited.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) arg_num
  character ( len = 80 ) arg_string
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) estimate
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: fx
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_more
  integer ( kind = 4 ) n_log2
  integer ( kind = 4 ), parameter :: n_log2_max = 25
  real ( kind = 8 ) quad
  real ( kind = 8 ) quad_more
  integer ( kind = 4 ) seed
  real ( kind = 8 ) volume
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BALL_VOLUME_MONTE_CARLO:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Estimate the volume of the unit sphere using'
  write ( *, '(a)' ) '  a Monte Carlo procedure.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the spatial dimension.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, arg_string )
    call s_to_i4 ( arg_string, dim_num, ierror, last )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BALL_VOLUME_MONTE_CARLO:'
    write ( *, '(a)' ) '  Enter the spatial dimension of the sphere'

    read ( *, '(a)' ) dim_num

  end if
!
!  Get the random number seed if it was supplied on the command line.
!
  if ( 2 <= arg_num ) then

    iarg = 2
    call getarg ( iarg, arg_string )
    call s_to_i4 ( arg_string, seed, ierror, last )

  else

    seed = 123456789
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BALL_VOLUME_MONTE_CARLO:'
    write ( *, '(a)' ) '  Using default seed for random number generator.'

  end if
!
!  Report user input.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' )  '  The spatial dimension is  ', dim_num
  write ( *, '(a,i12)' ) '  The random number seed is ', seed
!
!  Begin computation.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Log(N)         N      Estimate         Error'
  write ( *, '(a)' ) ' '

  quad = 0.0D+00
  volume = 2.0D+00**dim_num

  do n_log2 = 0, n_log2_max

    if ( n_log2 == 0 ) then
      quad = 0.0D+00
      n_more = 1
      n = 0
    else if ( n_log2 == 1 ) then
      n_more = 1
    else
      n_more = 2 * n_more
    end if

    allocate ( x(1:dim_num,1:n_more) )
    allocate ( fx(1:n_more) )

    call r8mat_uniform_01 ( dim_num, n_more, seed, x )
!
!  Rescale X from [0,1] to [-1,1].
!
    x(1:dim_num,1:n_more) = 2.0D+00 * x(1:dim_num,1:n_more) - 1.0D+00

    call sphere_indicator ( dim_num, n_more, x, fx )

    quad_more = sum ( fx(1:n_more) )

    deallocate ( fx )
    deallocate ( x )
!
!  Incorporate the new data into the totals.
!
    n = n + n_more
    quad = quad + quad_more

    estimate = volume * quad / real ( n, kind = 8 )
    call sphere_unit_volume_nd ( dim_num, exact )
    error = abs ( exact - estimate )
    write ( *, '(2x,i8,2x,i8,2x,g16.8,2x,g10.2)' ) n_log2, n, estimate, error

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(8x,a2,8x,a2,2x,g16.8,2x,g10.2)' ) 'oo', 'oo', exact, 0.0D+00
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BALL_VOLUME_MONTE_CARLO:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8's.
!
!    For now, the input quantity SEED is an integer variable.
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
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
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
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

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
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S 
!    used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) length
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
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end
subroutine sphere_indicator ( dim_num, point_num, x, fx )

!*****************************************************************************80
!
!! SPHERE_INDICATOR evaluates the unit sphere indicator function.
!
!  Discussion:
!
!    F(X) = 1 if X is on or inside the unit sphere, and 0 elsewhere.
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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer POINT_NUM, the number of points to evaluate.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the points.
!
!    Output, real ( kind = 8 ) FX(POINT_NUM), the unit sphere indicator 
!    function value.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) fx(point_num)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(dim_num,point_num)
!
!  Use the F95 WHERE if possible.
!
  if ( .false. ) then

    do j = 1, point_num
      if ( sum ( x(1:dim_num,j)**2 ) <= 1.0D+00 ) then
        fx(j) = 1.0D+00
      else
        fx(j) = 0.0D+00
      end if
    end do

  else

    where ( sum ( x**2, dim = 1 ) <= 1.0D+00 )
      fx = 1.0D+00
    elsewhere
      fx = 0.0D+00
    end where

  end if

  return
end
subroutine sphere_unit_volume_nd ( dim_num, volume )

!*****************************************************************************80
!
!! SPHERE_UNIT_VOLUME_ND computes the volume of a unit sphere in M-dimensions.
!
!  Discussion:
!
!    DIM_NUM  Volume
!
!    2             PI
!    3  (4/3)    * PI
!    4  (1/2)    * PI^2
!    5  (8/15)   * PI^2
!    6  (1/6)    * PI^3
!    7  (16/105) * PI^3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the sphere.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) volume

  volume = 1.0D+00

  if ( mod ( dim_num, 2 ) == 0 ) then
    m = dim_num / 2
    do i = 1, m
      volume = volume * pi / real ( i, kind = 8 )
    end do
  else
    m = ( dim_num - 1 ) / 2
    do i = 1, m
      volume = volume * pi * 2.0D+00
    end do
    do i = m + 1, 2 * m + 1
      volume = volume * 2.0D+00 / real ( i, kind = 8 )
    end do
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
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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
