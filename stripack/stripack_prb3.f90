program main

!*****************************************************************************80
!
!! MAIN is the main program for STRIPACK_PRB3.
!
!  Discussion:
!
!    STRIPACK_PRB3 looks at the Delaunay triangulation created by TRMESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 April 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: nrow = 6

  real    ( kind = 8 ), allocatable :: area ( : )
  real    ( kind = 8 ) area_sum
  real    ( kind = 8 ) dist(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) ltri(nrow,2*(n-2))
  integer ( kind = 4 ) near(n)
  integer ( kind = 4 ) next(n)
  integer ( kind = 4 ) nt
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ), parameter :: r = 1.0D+00
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) xyz(3,n)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPACK_PRB3'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Have TRMESH create the Delaunay triangulation of'
  write ( *, '(a)' ) '  points on a sphere.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Have TRLIST return the triangle vertices,'
  write ( *, '(a)' ) '  and the triangle neighbors.'
!
!  Choose a set of random values for X, Y, Z.
!
  seed = 123456789
  call uniform_on_sphere01_map ( 3, n, seed, xyz )

  call r8mat_transpose_print ( 3, n, xyz, '  Data points:' )
!
!  Compute the Delaunay triangulation.
!
  call trmesh ( n, xyz(1,:), xyz(2,:), xyz(3,:), list, lptr, lend, lnew, near, &
    next, dist, ier )

  if ( ier /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_PRB3 - Fatal error!'
    write ( *, '(a,i8)' ) '  TRMESH returns IER = ', ier
    stop
  end if
!
!  Create the triangle list.
!
  call trlist ( n, list, lptr, lend, nrow, nt, ltri, ier )

  if ( ier /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_PRB3 - Fatal error!'
    write ( *, '(a,i8)' ) '  TRLIST returns IER = ', ier
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of triangles = ', nt
!
!  Euler's Formula.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Check Euler''s formula:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Faces =    ', nt
  write ( *, '(a,i8)' ) '  Vertices = ', n
  write ( *, '(a,i8)' ) '  Edges =    ', ( 3 * nt ) / 2
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  F+V-E-2 =  ', nt + n - ( ( 3 * nt ) / 2 ) - 2
!
!  Print the triangulation.
!
  call i4mat_transpose_print ( 3, nt, ltri(1:3,:), '  Delaunay triangles:' )
!
!  Print the triangle neighbors.
!
  call i4mat_transpose_print ( 3, nt, ltri(4:6,:), '  Triangle neighbors:' )
!
!  Compute areas.
!
  allocate ( area(1:nt) )
  do i = 1, nt
    i1 = ltri(1,i)
    i2 = ltri(2,i)
    i3 = ltri(3,i)
    call stri_vertices_to_area ( r, xyz(1:3,i1), xyz(1:3,i2), xyz(1:3,i3), area(i) )
  end do

  call r8vec_print ( nt, area, '  Spherical area of triangles:' )

  area_sum = sum ( area (1:nt) )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Area sum = ', area_sum
  write ( *, '(a,g14.6)' ) '  4*PI =     ', 4.0D+00 * pi

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPACK_PRB3:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  deallocate ( area )

  stop
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
!7
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
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
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

      write ( *, '(2x,i3,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

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

      write ( *, '(2x,i3,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

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
    write ( *, '(2x,i3,a,g16.8)' ) i, ':', a(i)
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
subroutine stri_angles_to_area_3d ( r, a, b, c, area )

!*****************************************************************************80
!
!! STRI_ANGLES_TO_AREA_3D computes the area of a spherical triangle.
!
!  Discussion:
!
!    A sphere centered at 0 in 3D satisfies the equation:
!
!      X*X + Y*Y + Z*Z = R*R
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle is:
!
!      AREA = ( A + B + C - PI ) * R*R
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) A, B, C, the angles of the triangle.
!
!    Output, real ( kind = 8 ) AREA, the area of the spherical triangle.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) area
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) r
!
!  Apply Girard's formula.
!
  area = r * r * ( a + b + c - pi )

  return
end
subroutine stri_contains_point ( v1, v2, v3, p, contains )

!*****************************************************************************80
!
!! STRI_CONTAINS_POINT determines if a spherical triangle contains a point.
!
!  Discussion:
!
!    A sphere centered at 0 in 3D satisfies the equation:
!
!      X*X + Y*Y + Z*Z = R*R
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.  The inside of the triangle is defined by the fact
!    that the three points are listed in counterclockwise order. 
!    Here "counterclockwise" is with reference to an observer standing
!    outside the sphere.
!
!    If P is a point on the sphere, we say that the spherical triangle
!    "contains" P if P is in the interior of the spherical triangle.
!    We do not actually require that P be a point on the sphere.  Instead,
!    we consider the ray defined from the origin through P, which intersects
!    the sphere.  It is essentially this point of intersection we are
!    considering.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Input, real ( kind = 8 ) P(3), a point on the sphere, or the point on
!    the sphere determined by the ray from the origin through P.  P must
!    not be zero.
!
!    Output, real ( kind = 8 ) CONTAINS, is positive if the spherical triangle
!    contains P, zero if P is exactly on the boundary of the triangle, and
!    negative if P is outside the triangle.  
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real    ( kind = 8 ) contains
  real    ( kind = 8 ) p(dim_num)
  real    ( kind = 8 ) p_direction(dim_num)
  real    ( kind = 8 ) p_norm
  real    ( kind = 8 ) normal_direction(dim_num)
  real    ( kind = 8 ) normal_norm
  real    ( kind = 8 ) v1(dim_num)
  real    ( kind = 8 ) v2(dim_num)
  real    ( kind = 8 ) v3(dim_num)
!
!  Determine the normal vector to (V1,V2,V3), which is (V2-V1)x(V3-V13)..
!
  normal_direction(1) = ( v2(2) - v1(2) ) * ( v3(3) - v1(3) ) &
                      - ( v2(3) - v1(3) ) * ( v3(2) - v1(2) )

  normal_direction(2) = ( v2(3) - v1(3) ) * ( v3(1) - v1(1) ) &
                      - ( v2(1) - v1(1) ) * ( v3(3) - v1(3) )

  normal_direction(3) = ( v2(1) - v1(1) ) * ( v3(2) - v1(2) ) &
                      - ( v2(2) - v1(2) ) * ( v3(1) - v1(1) )

  normal_norm = sqrt ( sum ( normal_direction(1:dim_num)**2 ) )

  if ( normal_norm == 0.0D+00 ) then
    contains = - huge ( contains )
    return
  end if

  normal_direction(1:dim_num) = normal_direction(1:dim_num) / normal_norm
!
!  Determine the length of P.
!
  p_norm = sqrt ( sum ( p(1:dim_num)**2 ) )

  if ( p_norm == 0.0D+00 ) then
    contains = - huge ( contains )
    return
  end if

  p_direction(1:dim_num) = p_direction(1:dim_num) / p_norm
!
!  CONTAINS is the dot product of the normal vector to (V1,V2,V3)
!  against the unit direction vector defined by P.
!
  contains = dot_product ( normal_direction, p_direction )

  return
end
subroutine stri_sides_to_angles_3d ( r, as, bs, cs, a, b, c )

!*****************************************************************************80
!
!! STRI_SIDES_TO_ANGLES_3D computes spherical triangle angles in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) AS, BS, CS, the (geodesic) length of the 
!    sides of the triangle.
!
!    Output, real ( kind = 8 ) A, B, C, the spherical angles of the triangle.
!    Angle A is opposite the side of length AS, and so on.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) as
  real    ( kind = 8 ) asu
  real    ( kind = 8 ) b
  real    ( kind = 8 ) bs
  real    ( kind = 8 ) bsu
  real    ( kind = 8 ) c
  real    ( kind = 8 ) cs
  real    ( kind = 8 ) csu
  real    ( kind = 8 ) r
  real    ( kind = 8 ) ssu
  real    ( kind = 8 ) tan_a2
  real    ( kind = 8 ) tan_b2
  real    ( kind = 8 ) tan_c2

  asu = as / r
  bsu = bs / r
  csu = cs / r
  ssu = ( asu + bsu + csu ) / 2.0D+00

  tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) / & 
                  ( sin ( ssu ) * sin ( ssu - asu )     ) )

  a = 2.0D+00 * atan ( tan_a2 )

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) / & 
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) )

  b = 2.0D+00 * atan ( tan_b2 )

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) / & 
                  ( sin ( ssu ) * sin ( ssu - csu )     ) )

  c = 2.0D+00 * atan ( tan_c2 )

  return
end
subroutine stri_vertices_to_area ( r, v1, v2, v3, area )

!*****************************************************************************80
!
!! STRI_VERTICES_TO_AREA computes the area of a spherical triangle.
!
!  Discussion:
!
!    A sphere centered at 0 in 3D satisfies the equation:
!
!      X*X + Y*Y + Z*Z = R*R
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle is:
!
!      AREA = ( A + B + C - PI ) * R*R
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
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
!    Output, real ( kind = 8 ) AREA, the area of the spherical triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real    ( kind = 8 ) a
  real    ( kind = 8 ) area
  real    ( kind = 8 ) as
  real    ( kind = 8 ) b
  real    ( kind = 8 ) bs
  real    ( kind = 8 ) c
  real    ( kind = 8 ) cs
  real    ( kind = 8 ) r
  real    ( kind = 8 ) v1(dim_num)
  real    ( kind = 8 ) v2(dim_num)
  real    ( kind = 8 ) v3(dim_num)
!
!  Compute the lengths of the sides of the spherical triangle.
!
  call stri_vertices_to_sides ( r, v1, v2, v3, as, bs, cs )
!
!  Get the spherical angles.
!
  call stri_sides_to_angles_3d ( r, as, bs, cs, a, b, c )
!
!  Get the area
!
  call stri_angles_to_area_3d ( r, a, b, c, area )

  return
end
subroutine stri_vertices_to_centroid_3d ( r, v1, v2, v3, vs )

!*****************************************************************************80
!
!! STRI_VERTICES_TO_CENTROID_3D gets a spherical triangle centroid in 3D.
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
subroutine stri_vertices_to_sides ( r, v1, v2, v3, as, bs, cs )

!*****************************************************************************80
!
!! STRI_VERTICES_TO_SIDES_3D computes spherical triangle sides in 3D.
!
!  Discussion:
!
!    We can use the ACOS system call here, but the ARC_COSINE routine
!    will automatically take care of cases where the input argument is
!    (usually slightly) out of bounds.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the vertices of the spherical
!    triangle.
!
!    Output, real ( kind = 8 ) AS, BS, CS, the (geodesic) length of the sides
!    of the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real    ( kind = 8 ) arc_cosine
  real    ( kind = 8 ) as
  real    ( kind = 8 ) bs
  real    ( kind = 8 ) cs
  real    ( kind = 8 ) r
  real    ( kind = 8 ) v1(dim_num)
  real    ( kind = 8 ) v2(dim_num)
  real    ( kind = 8 ) v3(dim_num)

  as = r * arc_cosine ( dot_product ( v2(1:dim_num), v3(1:dim_num) ) / r**2 )
  bs = r * arc_cosine ( dot_product ( v3(1:dim_num), v1(1:dim_num) ) / r**2 )
  cs = r * arc_cosine ( dot_product ( v1(1:dim_num), v2(1:dim_num) ) / r**2 )

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
