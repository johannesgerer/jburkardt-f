program main

!*****************************************************************************80
!
!! MAIN is the main program for LAPACK_TEST.
!
!  Discussion:
!
!    LAPACK_TEST tests some real symmetric eigenproblem routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAPACK_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test some of the LAPACK routines for'
  write ( *, '(a)' ) '  real symmetric eigenproblems.'

  n = 1

  do i = 1, 4
    n = n * 4
    call test01 ( n )
  end do

  n = 1

  do i = 1, 4
    n = n * 4
    call test02 ( n )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAPACK_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)') ' '
  call timestamp ( )

  stop
end
subroutine test01 ( n )

!*****************************************************************************80
!
!! TEST01 tests DSYEV.
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
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) aq(n,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  character jobz
  real ( kind = 8 ) lambda(n)
  real ( kind = 8 ) lambda2(n)
  real ( kind = 8 ), parameter :: lambda_dev = 1.0D+00
  real ( kind = 8 ) lambda_max
  real ( kind = 8 ), parameter :: lambda_mean = 1.0D+00
  real ( kind = 8 ) lambda_min
  integer ( kind = 4 ) lwork
  real ( kind = 8 ) q(n,n)
  real ( kind = 8 ) q2(n,n)
  real ( kind = 8 ) r(n,n)
  integer ( kind = 4 ) :: seed = 12345
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) time_setup
  real ( kind = 8 ) time_solve
  character uplo
  real ( kind = 8 ) work(3*n-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For a double precision real matrix (D)'
  write ( *, '(a)' ) '  in symmetric storage mode (SY)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DSYEV computes the eigenvalues and eigenvectors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Matrix order = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) &
    '  Initialize random number generator using SEED = ', seed

  call random_initialize ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8SYMM_TEST will give us a symmetric matrix'
  write ( *, '(a)' ) '  with known eigenstructure.'

  call cpu_time ( t1 )

  call r8symm_test ( n, lambda_mean, lambda_dev, seed, a, q, lambda )

  call cpu_time ( t2 )

  time_setup = t2 - t1

  if ( n <= 5 ) then

    call r8mat_print ( n, n, a, '  The matrix A:' )

    call r8mat_print ( n, n, q, '  The eigenvector matrix Q:' )

  end if

  lambda_min = minval ( lambda(1:n) )
  lambda_max = maxval ( lambda(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  LAMBDA_MIN = ', lambda_min
  write ( *, '(a,g14.6)' ) '  LAMBDA_MAX = ', lambda_max

  if ( n <= 10 ) then

    call r8vec_print ( n, lambda, '  The eigenvalue vector LAMBDA:' )

  end if
!
!  Verify the claim that A*Q = Q * LAMBDA.
!
  if ( n <= 5 ) then

    aq(1:n,1:n) = matmul ( a(1:n,1:n), q(1:n,1:n) )

    do j = 1, n
      lambda2(j) = sqrt ( sum ( aq(1:n,j)**2 ) )
    end do

    call r8vec_print ( n, lambda2, '  The column norms of A*Q:' )

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now call DSYEV'
  write ( *, '(a)' ) '  and see if it can recover Q and LAMBDA.'
!
!  Copy the matrix.
!
  q2(1:n,1:n) = a(1:n,1:n)

  jobz = 'V'
  uplo = 'U'
  lwork = 3 * n - 1

  call cpu_time ( t1 )

  call dsyev ( jobz, uplo, n, q2, n, lambda2, work, lwork, info )
 
  call cpu_time ( t2 )

  time_solve = t2 - t1 

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Warning!'
    write ( *, '(a,i6)' ) '  DSYEV returned an error flag INFO = ', info
    return
  end if

  lambda_min = minval ( lambda2(1:n) )
  lambda_max = maxval ( lambda2(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  LAMBDA_MIN = ', lambda_min
  write ( *, '(a,g14.6)' ) '  LAMBDA_MAX = ', lambda_max

  if ( n <= 10 ) then
    call r8vec_print ( n, lambda2, '  The computed eigenvalues Lambda:' )
  end if

  if ( jobz == 'V' ) then

    if ( n <= 5 ) then
      call r8mat_print ( n, n, q2, '  The eigenvector matrix:' )

      r(1:n,1:n) = matmul ( a(1:n,1:n), q2(1:n,1:n) )

      do j = 1, n
        r(1:n,j) = r(1:n,j) - lambda2(j) * q2(1:n,j)
      end do

      call r8mat_print ( n, n, r, '  The residual (A-Lambda*I)*X:' )

    end if

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Setup time = ', time_setup
  write ( *, '(a,g14.6)' ) '  Solve time = ', time_solve

  return
end
subroutine test02 ( n )

!*****************************************************************************80
!
!! TEST02 tests DSYEVD.
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
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) aq(n,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iwork(3+5*n)
  integer ( kind = 4 ) j
  character jobz
  real ( kind = 8 ) lambda(n)
  real ( kind = 8 ) lambda2(n)
  real ( kind = 8 ), parameter :: lambda_dev = 1.0D+00
  real ( kind = 8 ) lambda_max
  real ( kind = 8 ), parameter :: lambda_mean = 1.0D+00
  real ( kind = 8 ) lambda_min
  integer ( kind = 4 ) liwork
  integer ( kind = 4 ) lwork
  real ( kind = 8 ) q(n,n)
  real ( kind = 8 ) q2(n,n)
  real ( kind = 8 ) r(n,n)
  integer ( kind = 4 ) :: seed = 12345
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) time_setup
  real ( kind = 8 ) time_solve
  character uplo
  real ( kind = 8 ) work(1+6*n+2*n*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For a double precision real matrix (D)'
  write ( *, '(a)' ) '  in symmetric storage mode (SY)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DSYEVD computes the eigenvalues and eigenvectors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Matrix order = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) &
    '  Initialize random number generator using SEED = ', seed

  call random_initialize ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8SYMM_TEST will give us a symmetric matrix'
  write ( *, '(a)' ) '  with known eigenstructure.'

  call cpu_time ( t1 )

  call r8symm_test ( n, lambda_mean, lambda_dev, seed, a, q, lambda )

  call cpu_time ( t2 )

  time_setup = t2 - t1

  if ( n <= 5 ) then

    call r8mat_print ( n, n, a, '  The matrix A:' )

    call r8mat_print ( n, n, q, '  The eigenvector matrix Q:' )

  end if

  lambda_min = minval ( lambda(1:n) )
  lambda_max = maxval ( lambda(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  LAMBDA_MIN = ', lambda_min
  write ( *, '(a,g14.6)' ) '  LAMBDA_MAX = ', lambda_max

  if ( n <= 10 ) then

    call r8vec_print ( n, lambda, '  The eigenvalue vector LAMBDA:' )

  end if
!
!  Verify the claim that A*Q = Q * LAMBDA.
!
  if ( n <= 5 ) then

    aq(1:n,1:n) = matmul ( a(1:n,1:n), q(1:n,1:n) )

    do j = 1, n
      lambda2(j) = sqrt ( sum ( aq(1:n,j)**2 ) )
    end do

    call r8vec_print ( n, lambda2, '  The column norms of A*Q:' )

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now call DSYEVD'
  write ( *, '(a)' ) '  and see if it can recover Q and LAMBDA.'
!
!  Copy the matrix.
!
  q2(1:n,1:n) = a(1:n,1:n)

  jobz = 'V'
  uplo = 'U'
  lwork = 1 + 6 * n + 2 * n * n
  liwork = 3 + 5 * n

  call cpu_time ( t1 )

  call dsyevd ( jobz, uplo, n, q2, n, lambda2, work, lwork, &
    iwork, liwork, info )
 
  call cpu_time ( t2 )

  time_solve = t2 - t1 

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Warning!'
    write ( *, '(a,i6)' ) '  DSYEVD returned an error flag INFO = ', info
    return
  end if

  lambda_min = minval ( lambda2(1:n) )
  lambda_max = maxval ( lambda2(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  LAMBDA_MIN = ', lambda_min
  write ( *, '(a,g14.6)' ) '  LAMBDA_MAX = ', lambda_max

  if ( n <= 10 ) then
    call r8vec_print ( n, lambda2, '  The computed eigenvalues Lambda:' )
  end if

  if ( jobz == 'V' ) then

    if ( n <= 5 ) then
      call r8mat_print ( n, n, q2, '  The eigenvector matrix:' )

      r(1:n,1:n) = matmul ( a(1:n,1:n), q2(1:n,1:n) )

      do j = 1, n
        r(1:n,j) = r(1:n,j) - lambda2(j) * q2(1:n,j)
      end do

      call r8mat_print ( n, n, r, '  The residual (A-Lambda*I)*X:' )

    end if

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Setup time = ', time_setup
  write ( *, '(a,g14.6)' ) '  Solve time = ', time_solve

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
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
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
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
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
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
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input/output, integer SEED, the "seed" value, which should NOT be 0.
!    On output, SEED has been updated.
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

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + huge ( seed )
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end
subroutine random_initialize ( seed_input )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!
!    The G95 compiler is currently generating internal compiler errors
!    when it tries to compile this routine.  Just one more exasperation
!    on the mountain of complications because of the ragged interface with
!    the nonstandard random number generator standard!
!
!    If you don't initialize the FORTRAN90 random number generator
!    routine RANDOM_NUMBER, which is used by calls like
!
!      call random_number ( harvest = x )
!
!    then its behavior is not specified.  That may be OK for you.  But
!    you may want to be able to force the same sequence to be generated
!    each time, or to force a different sequence.  
!
!    To control the sequence of random numbers, you need to set the seed.
!    In FORTRAN90, this is done by calling the RANDOM+SEED routine.
!    You can call it with no arguments, in fact.  But if you call
!    it with no arguments:
!
!      call random_seed ( )
!
!    then its behavior (or more particularly, the behavior of RANDOM_NUMBER)
!    is still not specified.  You might hope that the system will go out
!    and pick a nice random seed for you, but there's no guarantee. 
!
!
!    For example, on the DEC ALPHA, if you compile a program that calls
!    RANDOM_NUMBER, then every time you run it, you get the same sequence
!    of "random" values.  If you compile a program that calls RANDOM_SEED
!    with no arguments, and then calls RANDOM_NUMBER, you still get the
!    same sequence each time you run the program.
!
!    In order to actually try to scramble up the random number generator 
!    a bit, this routine goes through the tedious process of getting the 
!    size of the random number seed, making up values based on the current 
!    time, and setting the random number seed.
!
!    Unfortunately, the RANDOM_SEED routine has a very elastic definition.
!    It does not use a single scalar integer SEED.  Instead, it communicates
!    with the user through an integer vector whose size is not specified.
!    You actually have to "ask" the routine to tell you the size of this
!    vector.  Then you can fill up the vector with values to be used to
!    influence the seeding of the random number routine.  The details of
!    how the seed affects the sequence are also unspecified, but the only
!    thing we are reasonably confident about is that specifying the same
!    seed should result in the same sequence, and specifying different
!    seeds should result in different sequences!
!
!    I assume this is far more than you wanted to know.  (It's certainly
!    more than I wanted to know!)
!
!    The argument SEED is an input quantity only, so it is legal
!    to type
!
!      call random_initialize ( 0 )
!
!    or
!
!      call random_initialize ( 18867 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED_INPUT, the user's "suggestion" for a seed.
!    However, if the input value is 0, the routine will come up with
!    its own "suggestion", based on the system clock.
!
  implicit none

  integer ( kind = 4 ) count
  integer ( kind = 4 ) count_max
  integer ( kind = 4 ) count_rate
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_input
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
  real ( kind = 4 ) t
  integer ( kind = 4 ), parameter :: warm_up = 100

  seed = seed_input
!
!  Initialize the random seed routine.
!
  call random_seed ( )
!
!  Determine the size of the random number seed vector.
!
  call random_seed ( size = seed_size )
!
!  Allocate a vector of the right size to be used as a random seed.
!
  allocate ( seed_vector(seed_size) )
!
!  If the user supplied a SEED value, use that.
!
!  Otherwise, use the system clock value to make up a value that is
!  likely to change based on when this routine is called.
!
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
!  Set the seed vector.  We don't know the significance of the
!  individual entries of the internal seed vector, so we'll just set 
!  all entries to SEED.
!
  seed_vector(1:seed_size) = seed
!
!  Now call RANDOM_SEED, and tell it to use this seed vector.
!
  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times just to "warm it up".
!
  do i = 1, warm_up
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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
