program main

!*****************************************************************************80
!
!! MAIN is the main program for EISPACK_PRB2.
!
!  Discussion:
!
!    EISPACK_PRB2 does some symmetric eigenproblem tests on EISPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2008
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
  write ( *, '(a)' ) 'EISPACK_PRB2'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the EISPACK library.'
  write ( *, '(a)' ) '  Do some symmetric eigenproblem tests.'

  n = 1

  do i = 1, 4
    n = n * 4
    call test01 ( n )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EISPACK_PRB2'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)') ' '
  call timestamp ( )

  stop
end
subroutine test01 ( n )

!*****************************************************************************80
!
!! TEST01 tests RS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) aq(n,n)
  real    ( kind = 8 ) a2(n,n)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  real    ( kind = 8 ) lambda(n)
  real    ( kind = 8 ) lambda2(n)
  real    ( kind = 8 ), parameter :: lambda_dev = 1.0D+00
  real    ( kind = 8 ) lambda_max
  real    ( kind = 8 ), parameter :: lambda_mean = 1.0D+00
  real    ( kind = 8 ) lambda_min
  integer ( kind = 4 ) matz
  real    ( kind = 8 ) q(n,n)
  real    ( kind = 8 ) q2(n,n)
  real    ( kind = 8 ) r(n,n)
  integer ( kind = 4 ) :: seed = 12345
  real    ( kind = 8 ) t1
  real    ( kind = 8 ) t2
  real    ( kind = 8 ) time_setup
  real    ( kind = 8 ) time_solve

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  RS computes the eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) &
    '  Initialize random number generator using SEED = ', seed

  call random_initialize ( seed )

  write ( *, '(a,i12)' ) '  RANDOM_INITIALIZE returns suggested seed = ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SYMM_TEST will give us a symmetric matrix'
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
  write ( *, '(a)' ) '  Now call EISPACK routine RS'
  write ( *, '(a)' ) '  and see if it can recover Q and LAMBDA.'
!
!  Copy the matrix.
!
  a2(1:n,1:n) = a(1:n,1:n)

  matz = 1

  call cpu_time ( t1 )

  call rs ( n, a2, lambda2, matz, q2, ierror )

  call cpu_time ( t2 )

  time_solve = t2 - t1

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Warning!'
    write ( *, '(a,i8)' ) '  RS returned an error flag IERROR = ', ierror
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

  if ( matz /= 0 ) then

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
