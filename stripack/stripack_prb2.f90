program main

!*****************************************************************************80
!
!! MAIN is the main program for STRIPACK_PRB2.
!
!  Discussion:
!
!    STRIPACK_PRB2 demonstrates how to use STRIPACK data.
!
!    STRIPACK can compute the Voronoi diagram for data on a sphere.
!
!    This routine has STRIPACK compute the Voronoi diagram, then
!    takes just a few of the "interesting" arrays, and uses them
!    to visit every Voronoi polygon.  Just to prove that is what
!    we are doing, we compute the area of each subtriangle of the
!    polygons, and sum them.  This should be equal to the total area
!    of the sphere, 4 * PI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Local Parameters:
!
!    Local parameter, integer ( kind = 4 ) N, the number of generators of
!    Voronoi cells.
!
!    Local parameter, real ( kind = 8 ) X(N), Y(N), Z(N), randomly chosen 
!    coordinates for the generators.  These points are normalized so that they 
!    lie on the unit sphere.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) area(n)
  real    ( kind = 8 ) centroid(3,n)
  integer ( kind = 4 ) first(n+1)
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ), allocatable :: list(:)
  integer ( kind = 4 ) listc(6*(n-2))
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) order(n)
  integer ( kind = 4 ) order_sum
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) xc(2*(n-2))
  real    ( kind = 8 ) yc(2*(n-2))
  real    ( kind = 8 ) zc(2*(n-2))
  real    ( kind = 8 ) xyz(3,n)
  real    ( kind = 8 ), allocatable :: xyzv(:,:)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPACK_PRB2'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A sample application of the data produced by'
  write ( *, '(a)' ) '  STRIPACK.  Here, we have STRIPACK compute the'
  write ( *, '(a)' ) '  Voronoi diagram of a set of points on the unit'
  write ( *, '(a)' ) '  sphere, and then we do a simple check by computing'
  write ( *, '(a)' ) '  the polygonal area and summing.'
!
!  Choose a set of random values for X, Y, Z.
!
  seed = 123456789
  call uniform_on_sphere01_map ( 3, n, seed, xyz )

  call r8mat_transpose_print ( 3, n, xyz, '  Data points:' )
!
!  Now we compute the Voronoi information on the sphere.
!
  call voronoi_get ( n, xyz(1,:), xyz(2,:), xyz(3,:), nt, xc, yc, zc, &
    lend, listc, lptr )

  allocate ( xyzv(1:3,1:nt) )

  xyzv(1,1:nt) = xc(1:nt)
  xyzv(2,1:nt) = yc(1:nt)
  xyzv(3,1:nt) = zc(1:nt)

  call r8mat_transpose_print ( 3, nt, xyzv, '  Voronoi vertices' )
!
!  Get the order of each Voronoi polygon.
!
  call voronoi_order ( n, lend, lptr, order )

  call i4vec_print ( n, order, '  Voronoi polygon orders:' )
!
!  Get the Voronoi polygons as a list.
!
  order_sum = sum ( order(1:n) )
  allocate ( list(1:order_sum) )

  call voronoi_polygons ( n, order_sum, lend, listc, lptr, first, list )

  call i4list_print ( n, first, order_sum, list, '  Voronoi polygons:' )
!
!  Compute the area of each Voronoi polygon.
!
  call voronoi_area ( n, xyz(1,:), xyz(2,:), xyz(3,:), xc, yc, zc, lend, &
    listc, lptr, area )

  call r8vec_print ( n, area, '  Voronoi polygon areas:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) &
    '  Sphere area from Voronoi polygons =  ', sum ( area(1:n) )
  write ( *, '(a,f12.6)' ) &
    '  Exact area from spherical geometry = ', 4.0D+00 * pi
!
!  Compute the centroid of each Voronoi polygon.
!
  call voronoi_centroids ( n, xyz(1,:), xyz(2,:), xyz(3,:), xc, yc, zc, lend, &
    listc, lptr, centroid )

  call r8mat_transpose_print ( 3, n, centroid, '  Voronoi polygon centroids:' )
!
!  "Visit" each Voronoi polygon, and each triangle of that polygon.
!
  call voronoi_traverse ( n, xyz(1,:), xyz(2,:), xyz(3,:), xc, yc, zc, lend, &
    listc, lptr )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPACK_PRB2:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  deallocate ( list )
  deallocate ( xyzv )

  stop
end
subroutine i4list_print ( n, first, list_num, list, title )

!*****************************************************************************80
!
!! I4LIST_PRINT prints an I4LIST.
!
!  Discussion:
!
!    An I4LIST is a list of integers grouped into N segments.
!    An index vector locates the first entry of each segment.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of segments.
!
!    Input, integer ( kind = 4 ) FIRST(N+1), indexes the first entry
!    of each segment.
!
!    Input, integer ( kind = 4 ) LIST_NUM, the number of entries.
!
!    Input, integer ( kind = 4 ) LIST(LIST_NUM), the data.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) list_num
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) first(n+1)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  integer   ( kind = 4 ) list(list_num)
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n

    do jlo = first(i), first(i+1) - 1, 5
      jhi = min ( jlo + 4, first(i+1) - 1 )
      if ( jlo == first(i) ) then
        write ( *, '(i5,a,5(2x,i8))' ) i, ':', list(jlo:jhi)
      else
        write ( *, '(6x,  5(2x,i8))' )         list(jlo:jhi)
      end if
    end do

  end do

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(m,n)
  character ( len = * )  title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 10
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(m,n)
  character ( len = 8 )  ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

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
      write ( ctemp(i2), '(i8)' ) i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i5,a,i10)' ) i, ':', a(i)
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
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

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
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
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
!    10 September 2009
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

  integer   ( kind = 4 ), parameter :: incx = 5
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

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
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED
!    in order to set the seed of the random number generator.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is 
!    negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer ( kind = 4 ) SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range 
!    of entries of X that we need to compute.  This starts off as 1:N, but 
!    is adjusted if we have a saved value that can be immediately stored
!    in X(1), and so on.
!
!    Local, real ( kind = 8 ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) m
  integer ( kind = 4 ), save :: made = 0
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) r(n+1)
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ), save :: saved = 0
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: two = 2
  real    ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
  real    ( kind = 8 ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, two ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i5,a,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
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
  real    ( kind = 8 ) r(n)

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
subroutine stri_vertices_to_centroid ( r, v1, v2, v3, vs )

!*****************************************************************************80
!
!! STRI_VERTICES_TO_CENTROID gets a spherical triangle centroid in 3D.
!
!  Discussion:
!
!    A sphere centered at 0 in 3D satisfies the equation:
!
!      X*X + Y*Y + Z*Z = R*R
!
!    A spherical triangle is specified by three points on the sphere.
!
!    The (true) centroid of a spherical triangle is the point
!
!      VT = (XT,YT,ZT) = Integral ( X, Y, Z ) dArea / Integral 1 dArea
!
!    Note that the true centroid does NOT, in general, lie on the sphere.  
!
!    The "flat" centroid VF is the centroid of the planar triangle defined by
!    the vertices of the spherical triangle.
!
!    The "spherical" centroid VS of a spherical triangle is computed by
!    the intersection of the geodesic bisectors of the triangle angles.
!    The spherical centroid lies on the sphere.
!
!    VF, VT and VS lie on a line through the center of the sphere.  We can
!    easily calculate VF by averaging the vertices, and from this determine
!    VS by normalizing.
!
!    Of course, we still will not have actually computed VT, which lies
!    somewhere between VF and VS!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real ( kind = 8 ) VS(3), the coordinates of the "spherical
!    centroid" of the spherical triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real    ( kind = 8 ) norm
  real    ( kind = 8 ) r
  real    ( kind = 8 ) v1(dim_num)
  real    ( kind = 8 ) v2(dim_num)
  real    ( kind = 8 ) v3(dim_num)
  real    ( kind = 8 ) vs(dim_num)

  vs(1:dim_num) = ( v1(1:dim_num) + v2(1:dim_num) + v3(1:dim_num) ) / 3.0D+00

  norm = sqrt ( sum ( vs(1:dim_num)**2 ) )

  vs(1:dim_num) = r * vs(1:dim_num) / norm

  return
end
subroutine uniform_on_sphere01_map ( dim_num, n, seed, x )

!*****************************************************************************80
!
!! UNIFORM_ON_SPHERE01_MAP maps uniform points onto the unit sphere.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    This procedure is valid for any spatial dimension DIM_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Russell Cheng,
!    Random Variate Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998, pages 168.
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity 
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,N), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real    ( kind = 8 ) norm
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) x(dim_num,n)

  do j = 1, n
!
!  Fill a vector with normally distributed values.
!
    call r8vec_normal_01 ( dim_num, seed, x(1:dim_num,j) )
!
!  Compute the length of the vector.
!
    norm = sqrt ( sum ( x(1:dim_num,j)**2 ) )
!
!  Normalize the vector.
!
    x(1:dim_num,j) = x(1:dim_num,j) / norm

  end do

  return
end
subroutine voronoi_area ( n, x, y, z, xc, yc, zc, lend, listc, lptr, area )

!*****************************************************************************80
!
!! VORONOI_AREA computes the area of each polygon in a Voronoi diagram.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes, and Voronoi polygons.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XC(6*(N-2)), YC(6*(N-2)), ZC(6*(N-2)), the 
!    coordinates of the vertices.
!
!    Input, integer ( kind = 4 ) LEND(N), points to the "first" vertex in the
!    Voronoi polygon around a particular node.
!
!    Input, integer ( kind = 4 ) LISTC(6*(N-2)), the Voronoi vertex indices.
!
!    Input, integer ( kind = 4 ) LPTR(6*(N-2)), given a vertex, returns the next
!    vertex in the Voronoi polygon.  (The vertex numbering is done
!    in such a way that the physical vertex has three distince indices,
!    depending on which polygon we are considering.  Thus, it is always
!    possible to answer the question "which is the next vertex from this
!    one?" because the vertex index also tells us what polygon we are in.)
!
!    Output, real ( kind = 8 ) AREA(N), the area of each polygon.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) area(n)
  real    ( kind = 8 ) areas
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) listc(6*(n-2))
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_last
  integer ( kind = 4 ) node_new
  integer ( kind = 4 ) node_stop
  real    ( kind = 8 ) v1(3)
  real    ( kind = 8 ) v2(3)
  real    ( kind = 8 ) v3(3)
  integer ( kind = 4 ) vertex_last
  integer ( kind = 4 ) vertex_new
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xc(6*(n-2))
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) yc(6*(n-2))
  real    ( kind = 8 ) z(n)
  real    ( kind = 8 ) zc(6*(n-2))

  do node = 1, n

    area(node) = 0.0D+00

    node_stop = lend(node)

    node_new = node_stop
    vertex_new = listc(node_new)

    do

      node_last = node_new
      node_new = lptr(node_last)

      vertex_last = vertex_new
      vertex_new = listc(node_new)

      v1(1:3) = (/ x(node),         y(node),         z(node)         /)
      v2(1:3) = (/ xc(vertex_last), yc(vertex_last), zc(vertex_last) /)
      v3(1:3) = (/ xc(vertex_new),  yc(vertex_new),  zc(vertex_new)  /)

      area(node) = area(node) + areas ( v1, v2, v3 )

      if ( node_new == node_stop ) then
        exit
      end if

    end do

  end do

  return
end
subroutine voronoi_centroids ( n, x, y, z, xc, yc, zc, lend, listc, lptr, &
  centroid )

!*****************************************************************************80
!
!! VORONOI_CENTROIDS computes the centroids of each polygon in a Voronoi diagram.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes, and Voronoi polygons.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XC(6*(N-2)), YC(6*(N-2)), ZC(6*(N-2)), the 
!    coordinates of the vertices.
!
!    Input, integer ( kind = 4 ) LEND(N), points to the "first" vertex in the
!    Voronoi polygon around a particular node.
!
!    Input, integer ( kind = 4 ) LISTC(6*(N-2)), the Voronoi vertex indices.
!
!    Input, integer ( kind = 4 ) LPTR(6*(N-2)), given a vertex, returns the next
!    vertex in the Voronoi polygon.  (The vertex numbering is done
!    in such a way that the physical vertex has three distince indices,
!    depending on which polygon we are considering.  Thus, it is always
!    possible to answer the question "which is the next vertex from this
!    one?" because the vertex index also tells us what polygon we are in.)
!
!    Output, real ( kind = 8 ) CENTROID(3,N), the centroid of each Voronoi polygon.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) centroid(3,n)
  real    ( kind = 8 ) areas
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) listc(6*(n-2))
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_last
  integer ( kind = 4 ) node_new
  integer ( kind = 4 ) node_stop
  real    ( kind = 8 ) norm
  real    ( kind = 8 ) r
  real    ( kind = 8 ) stri_area
  real    ( kind = 8 ) stri_centroid(3)
  real    ( kind = 8 ) v1(3)
  real    ( kind = 8 ) v2(3)
  real    ( kind = 8 ) v3(3)
  integer ( kind = 4 ) vertex_last
  integer ( kind = 4 ) vertex_new
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xc(6*(n-2))
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) yc(6*(n-2))
  real    ( kind = 8 ) z(n)
  real    ( kind = 8 ) zc(6*(n-2))

  centroid(1:3,1:n) = 0.0D+00
  r = 1.0D+00

  do node = 1, n

    node_stop = lend(node)

    node_new = node_stop
    vertex_new = listc(node_new)

    do

      node_last = node_new
      node_new = lptr(node_last)

      vertex_last = vertex_new
      vertex_new = listc(node_new)

      v1(1:3) = (/ x(node),         y(node),         z(node)         /)
      v2(1:3) = (/ xc(vertex_last), yc(vertex_last), zc(vertex_last) /)
      v3(1:3) = (/ xc(vertex_new),  yc(vertex_new),  zc(vertex_new)  /)

      stri_area = areas ( v1, v2, v3 )

      call stri_vertices_to_centroid ( r, v1, v2, v3, stri_centroid )

      centroid(1:3,node) = centroid(1:3,node) + stri_area * stri_centroid(1:3)

      if ( node_new == node_stop ) then
        exit
      end if

    end do

    norm = sqrt ( sum ( centroid(1:3,node)**2 ) )
    centroid(1:3,node) = centroid(1:3,node) / norm

  end do

  return
end
subroutine voronoi_get ( n, x, y, z, nt, xc, yc, zc, lend, listc, lptr )

!*****************************************************************************80
!
!! VORONOI_GET calls STRIPACK routines to get Voronoi information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the coordinates of points 
!    on the sphere.
!
!    Output, integer ( kind = 4 ) NT, the number of Delaunay triangles
!    and Voronoi vertices.
!
!    Output, real ( kind = 8 ) XC(6*(N-2)), YC(6*(N-2)), ZC(6*(N-2)), the 
!    coordinates of the vertices.
!
!    Output, integer ( kind = 4 ) LEND(N), points to the "first" vertex in the
!    Voronoi polygon around a particular node.
!
!    Output, integer ( kind = 4 ) LISTC(6*(N-2)), the Voronoi vertex indices.
!
!    Output, integer ( kind = 4 ) LPTR(6*(N-2)), given a vertex, returns the 
!    next vertex in the Voronoi polygon.  (The vertex numbering is done
!    in such a way that the physical vertex has three distince indices,
!    depending on which polygon we are considering.  Thus, it is always
!    possible to answer the question "which is the next vertex from this
!    one?" because the vertex index also tells us what polygon we are in.)
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: nrow = 9

  real    ( kind = 8 ) ds(n)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iwk(2*n)
  integer ( kind = 4 ) lbtri(6,n)
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) listc(6*(n-2)) 
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) ltri(nrow,2*(n-2))
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nt
  real    ( kind = 8 ) rc(2*(n-2))
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xc(2*(n-2))
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) yc(2*(n-2))
  real    ( kind = 8 ) z(n)
  real    ( kind = 8 ) zc(2*(n-2))
!
!  Create the triangulation.
!
  call trmesh ( n, x, y, z, list, lptr, lend, lnew, iwk, iwk(n+1), ds, ierror )

  if ( ierror == -2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_GET - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  The first 3 nodes are collinear.'
    stop
  end if

  if ( 0 < ierror ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_GET - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  Duplicate nodes encountered.'
    stop
  end if
!
!  Create a triangle list.
!
  call trlist ( n, list, lptr, lend, nrow, nt, ltri, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_GET - Fatal error!'
    write ( *, '(a)' ) '  Error in TRLIST.'
    stop
  end if

  call i4mat_transpose_print ( nrow, nt, ltri, '  Vertices/Triangles/Arcs:' )
!
!  Construct the Voronoi diagram.
!
!  Note that the triangulation data structure is altered if NB > 0.
!
  call crlist ( n, n, x, y, z, list, lend, lptr, lnew, &
    lbtri, listc, nb, xc, yc, zc, rc, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_GET - Fatal error!'
    write ( *, '(a)' ) '  Error in CRLIST.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  return
end
subroutine voronoi_order ( n, lend, lptr, order )

!*****************************************************************************80
!
!! VORONOI_ORDER computes the order of each polygon in a Voronoi diagram.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes, and Voronoi polygons.
!
!    Input, integer ( kind = 4 ) LEND(N), points to the "first" vertex in the
!    Voronoi polygon around a particular node.
!
!    Input, integer ( kind = 4 ) LPTR(6*(N-2)), given a vertex, returns the next
!    vertex in the Voronoi polygon.  (The vertex numbering is done
!    in such a way that the physical vertex has three distince indices,
!    depending on which polygon we are considering.  Thus, it is always
!    possible to answer the question "which is the next vertex from this
!    one?" because the vertex index also tells us what polygon we are in.)
!
!    Output, integer ( kind = 4 ) ORDER(N), the order of each polygon.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_last
  integer ( kind = 4 ) node_new
  integer ( kind = 4 ) node_stop
  integer ( kind = 4 ) order(n)

  do node = 1, n

    order(node) = 0;

    node_stop = lend(node)
    node_new = node_stop

    do

      node_last = node_new
      node_new = lptr(node_last)

      order(node) = order(node) + 1

      if ( node_new == node_stop ) then
        exit
      end if

    end do

  end do

  return
end
subroutine voronoi_polygons ( n, list_num, lend, listc, lptr, first, list )

!*****************************************************************************80
!
!! VORONOI_POLYGONS creates a list of Voronoi polygons.
!
!  Discussion:
!
!    STRIPACK defines a data structure recording the location of
!    the vertices of the Voronoi diagram, and their connectivity.
!    The purpose of this routine is to construct a simplified data structure
!    that lists the indices of the Voronoi vertices that form each 
!    Voronoi polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes, and Voronoi polygons.
!
!    Input, integer ( kind = 4 ) LIST_NUM, the number of entries to be stored
!    in LIST.
!
!    Input, integer ( kind = 4 ) LEND(N), points to the "first" vertex in the
!    Voronoi polygon around a particular node.
!
!    Input, integer ( kind = 4 ) LISTC(6*(N-2)), the Voronoi vertex indices.
!
!    Input, integer ( kind = 4 ) LPTR(6*(N-2)), given a vertex, returns the 
!    next vertex in the Voronoi polygon.  (The vertex numbering is done
!    in such a way that the physical vertex has three distince indices,
!    depending on which polygon we are considering.  Thus, it is always
!    possible to answer the question "which is the next vertex from this
!    one?" because the vertex index also tells us what polygon we are in.)
!
!    Output, integer FIRST(N+1), for each polygon, points to the location
!    in LIST of the index.
!
!    Output, integer LIST(LIST_NUM), the list of vertices that form each 
!    polygon.
!
  implicit none

  integer ( kind = 4 ) list_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) first(n+1)
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(list_num)
  integer ( kind = 4 ) listc(6*(n-2))
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_new
  integer ( kind = 4 ) node_stop
  integer ( kind = 4 ) used

  used = 0

  do node = 1, n

    first(node) = used + 1

    node_stop = lend(node)
    node_new = node_stop

    used = used + 1
    list(used) = listc(node_new)

    do

      node_new = lptr(node_new)

      if ( node_new == node_stop ) then
        exit
      end if

      used = used + 1
      list(used) = listc(node_new)

    end do

  end do

  first(n+1) = used + 1

  return
end
subroutine voronoi_traverse ( n, x, y, z, xc, yc, zc, lend, listc, lptr )

!*****************************************************************************80
!
!! VORONOI_TRAVERSE traverses the polygons in a Voronoi diagram.
!
!  Discussion:
!
!    STRIPACK defines a data structure recording the location of
!    the vertices of the Voronoi diagram, and their connectivity.
!    The purpose of this routine is to "visit" each polygon, and,
!    in fact, each subtriangle of each polygon.  Such a procedure
!    would be done when estimating an integral by quadrature, for instance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes, and Voronoi polygons.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XC(6*(N-2)), YC(6*(N-2)), ZC(6*(N-2)), the 
!    coordinates of the vertices.
!
!    Input, integer ( kind = 4 ) LEND(N), points to the "first" vertex in the
!    Voronoi polygon around a particular node.
!
!    Input, integer ( kind = 4 ) LISTC(6*(N-2)), the Voronoi vertex indices.
!
!    Input, integer ( kind = 4 ) LPTR(6*(N-2)), given a vertex, returns the 
!    next vertex in the Voronoi polygon.  (The vertex numbering is done
!    in such a way that the physical vertex has three distince indices,
!    depending on which polygon we are considering.  Thus, it is always
!    possible to answer the question "which is the next vertex from this
!    one?" because the vertex index also tells us what polygon we are in.)
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) area_polygon
  real    ( kind = 8 ) area_triangle
  real    ( kind = 8 ) areas
  integer ( kind = 4 ) index_polygon
  integer ( kind = 4 ) index_triangle
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) listc(6*(n-2))
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_last
  integer ( kind = 4 ) node_new
  integer ( kind = 4 ) node_stop
  integer ( kind = 4 ) order
  real    ( kind = 8 ) r
  real    ( kind = 8 ) v1(3)
  real    ( kind = 8 ) v2(3)
  real    ( kind = 8 ) v3(3)
  integer ( kind = 4 ) vertex_last
  integer ( kind = 4 ) vertex_new
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xc(6*(n-2))
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) yc(6*(n-2))
  real    ( kind = 8 ) z(n)
  real    ( kind = 8 ) zc(6*(n-2))

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VORONOI_TRAVERSE'
  write ( *, '(a)' ) '  Visit each Voronoi polygon.'
  write ( *, '(a)' ) '  Compute the (spherical) area of each subtriangle'
  write ( *, '(a)' ) '  Add up to get the area of each polygon.'
!
!  To access every polygon, start by accessing a particular node.
!
!  The Voronoi polygon around a node NODE has a pointer LEND(NODE) to the 
!  first (or last) vertex of the Voronoi polygon around NODE.
!
!  To access all the vertices of the polygon in order, start at the
!  special vertex, and then repeatedly use the LPTR array to get the
!  next vertex on the polygon.  Stop when you return to LEND(NODE).
!
!  To subdivide the polygon into triangles, use NODE, VERTEX_LAST,
!  and VERTEX.
!
!  To get the coordinates of these points:
!
!    NODE ==>        X(NODE),         Y(NODE),         Z(NODE).
!
!    VERTEX_LAST ==> XC(VERTEX_LAST), YC(VERTEX_LAST), ZC(VERTEX_LAST)
!    VERTEX      ==> XC(VERTEX     ), YC(VERTEX     ), ZC(VERTEX     ) 
!
  index_polygon = 0

  do node = 1, n

    area_polygon = 0.0D+00
    index_triangle = 0
    order = 0
    
    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Polygon ', node

    node_stop = lend(node)

    node_new = node_stop

    vertex_new = listc(node_new)
!
!  Each iteration of this DO walks along one side of the polygon,
!  considering the subtriangle NODE --> VERTEX_LAST --> VERTEX.
!
    do

      index_triangle = index_triangle + 1
      order = order + 1

      node_last = node_new
      node_new = lptr(node_last)

      vertex_last = vertex_new
      vertex_new = listc(node_new)
!
!  Here is a good place to process information about the polygon side 
!
!   VERTEX_LAST --> VERTEX 
!
!  or about the subtriangle
!
!   NODE --> VERTEX_LAST --> VERTEX.
!
      r = 1.0D+00
      v1(1:3) = (/ x(node),         y(node),         z(node)         /)
      v2(1:3) = (/ xc(vertex_last), yc(vertex_last), zc(vertex_last) /)
      v3(1:3) = (/ xc(vertex_new),  yc(vertex_new),  zc(vertex_new)      /)

      area_triangle = areas ( v1, v2, v3 )

      area_polygon = area_polygon + area_triangle

      write ( *, '(a,2x,i8,2x,a,2x,g14.6)' ) &
        '  Subtriangle ', index_triangle, '  area = ', area_triangle
!
!  Now if we have reached the vertex where we started, we are done with
!  this polygon.
!
      if ( node_new == node_stop ) then
        exit
      end if

    end do

    write ( *, '(a,2x,19x,g14.6)' ) '  Polygon area =', area_polygon

  end do

  return
end
