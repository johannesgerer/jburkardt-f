program main

!*****************************************************************************80
!
!! MAIN is the main program for LINPACK_D_PRB.
!
!  Discussion:
!
!    LINPACK_D_PRB calls the double precision real LINPACK test routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINPACK_D_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LINPACK_D library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )

  call test20 ( )
  call test21 ( )
  call test22 ( )
  call test23 ( )
  call test24 ( )
  call test25 ( )
  call test26 ( )
  call test27 ( )
  call test28 ( )
  call test29 ( )

  call test30 ( )
  call test31 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINPACK_D_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DCHDC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real ( kind = 8 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For real ( kind = 8 ), general storage,'
  write ( *, '(a)' ) '  DCHDC computes the Cholesky decomposition.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of equations is N = ', n
!
!  Set the values of the matrix A.
!
  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = 2.0D+00
  end do

  do i = 1, n-1
    a(i,i+1) = -1.0D+00
  end do

  do i = 2, n
    a(i-1,i) = -1.0D+00
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) a(i,1:n)
  end do
!
!  Decompose the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Decompose the matrix.'

  job = 0
  ipvt(1:n) = 0

  call dchdc ( a, lda, n, work, ipvt, job, info )

  if ( info /= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  DCHDC returned INFO = ', info
    write ( *, '(a)' ) '  This means the matrix is not positive definite.'
    return
  end if
!
!  Zero out the lower diagonal.
!
  do i = 2, n
    do j = 1, i-1
      a(i,j) = 0.0D+00
    end do
  end do
!
!  Print the factorization.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Cholesky factor U:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) a(i,1:n)
  end do
!
!  Compute the Cholesky product.
!
  a(1:n,1:n) = matmul ( transpose ( a(1:n,1:n) ), a(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product U'' * U: '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) a(i,1:n)
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests DCHEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n
  integer ( kind = 4 ), parameter :: nz = 1

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) s(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) work(n)
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For double precision real general storage,'
  write ( *, '(a)' ) '  DCHEX can shift columns in a Cholesky factorization.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of equations is N = ', n
!
!  Set the values of the matrix A.
!
  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = 2.0D+00
  end do

  do i = 1, n-1
    a(i,i+1) = -1.0D+00
  end do

  do i = 2, n
    a(i-1,i) = -1.0D+00
  end do

  do i = 1, n
    z(i) = real ( i, kind = 8 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) a(i,1:n)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The vector Z:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,g14.6)' ) z(i)
  end do
!
!  Decompose the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Decompose the matrix.'

  job = 0
  ipvt(1:n) = 0

  call dchdc ( a, lda, n, work, ipvt, job, info )

  if ( info /= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  DCHDC returned INFO = ', info
    write ( *, '(a)' ) '  This means the matrix is not positive definite.'
    return
  end if
!
!  Zero out the lower diagonal.
!
  do i = 2, n
    do j = 1, i-1
      a(i,j) = 0.0D+00
    end do
  end do
!
!  Print the factorization.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Cholesky factor U:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) a(i,1:n)
  end do
!
!  Right circular shift columns L through K.
!
  k = 1
  l = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  Right circular shift columns K  = ', k, &
    ' through L = ', l

  job = 1
  call dchex ( a, lda, n, k, l, z, n, nz, c, s, job )
!
!  Left circular shift columns K+1 through L.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  Left circular shift columns K+1 = ', k+1, &
    ' through L = ', l

  job = 2
  call dchex ( a, lda, n, k+1, l, z, n, nz, c, s, job )
!
!  Print the factorization.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The shifted Cholesky factor U:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) a(i,1:n)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The shifted vector Z:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,g14.6)' ) z(i)
  end do
!
!  Compute the Cholesky product.
!
  a(1:n,1:n) = matmul ( transpose ( a(1:n,1:n) ), a(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The shifted product U'' * U: '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) a(i,1:n)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests DCHUD and DTRSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: p = 20
  integer ( kind = 4 ), parameter :: ldr = p
  integer ( kind = 4 ), parameter :: nz = 1

  real ( kind = 8 ) b(p)
  real ( kind = 8 ) beta(p)
  real ( kind = 8 ) c(p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real ( kind = 8 ) r(ldr,p)
  real ( kind = 8 ) rho(nz)
  real ( kind = 8 ) row(p)
  real ( kind = 8 ) s(p)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(p)
  real ( kind = 8 ) y(nz)
  real ( kind = 8 ) z(p,nz)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For double precision real general storage,'
  write ( *, '(a)' ) '  DCHUD updates a Cholesky decomposition.'
  write ( *, '(a)' ) '  DTRSL can solve a triangular linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we use DCHUD to solve a'
  write ( *, '(a)' ) '  least squares problem R * b = z.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of equations is P = ', p
!
!  Initialize.
!
  r(1:p,1:p) = 0.0D+00
  z(1:p,1:nz) = 0.0D+00
  do i = 1, p
    x(i) = real ( i, kind = 8 )
  end do
!
!  Use DCHUD to form R, Z and RHO by adding X and Y a row at a time.
!  X is a row of the least squares matrix and Y the right hand side.
!
  seed = 123456789

  do i = 1, p
    call r8mat_uniform_01 ( 1, p, seed, row )
    y(1) = dot_product ( row(1:p), x(1:p) )
    rho(1) = 0.0D+00
    call dchud ( r, ldr, p, row, z, p, nz, y, rho, c, s )
  end do
!
!  Generate the least squares solution, b = inverse ( R ) * Z.
!
  do j = 1, nz

    b(1:p) = z(1:p,j)
    job = 01

    call dtrsl ( r, ldr, p, b, job, info )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Solution vector # ', j
    write ( *, '(a)' ) '  (Should be (1,2,3...,n))'
    write ( *, '(a)' ) ' '

    do i = 1, p
      if ( i <= 5 .or. p-5 < i ) then
        write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
      end if
      if ( i == 5 ) then
        write ( *, '(a)' ) '  ......  ..............'
      end if
    end do

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests DGBCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: lda = 2 * ml + mu + 1

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) rcond
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For a banded matrix in general format,'
  write ( *, '(a)' ) '  DGBCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  m = ml + mu + 1
  write ( *, '(a,i8)' ) '  The bandwidth of the matrix is ', m

  do j = 1, n
    a(m-1,j) = -1.0D+00
    a(m,j) = 2.0D+00
    a(m+1,j) = -1.0D+00
  end do
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call dgbco ( a, lda, n, ml, mu, ipivot, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition = ', rcond

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests DGBFA and DGBSL.
!
!  Discussion:
!
!    The problem solved here is a larger version of this one:
!
!    Matrix A is ( 2 -1  0  0  0)    right hand side B is  (1)
!                (-1  2 -1  0  0)                          (0)
!                ( 0 -1  2 -1  0)                          (0)
!                ( 0  0 -1  2 -1)                          (0)
!                ( 0  0  0 -1  2)                          (1)
!
!
!    Solution is   (1)
!                  (1)
!                  (1)
!                  (1)
!                  (1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Local Parameters:
!
!    N is the number of equations.
!
!    ML is the number of subdiagonals,
!    MU the number of superdiagonals.
!
!    LDA is the leading dimension of the array used to store the
!    matrix, which must be at least 2*ML+MU+1.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: lda = 2 * ml + mu + 1

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For a banded matrix in general format,'
  write ( *, '(a)' ) '  DGBFA factors the matrix,'
  write ( *, '(a)' ) '  DGBSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the right hand side B.
!
  b(1) = 1.0D+00
  b(2:n-1) = 0.0D+00
  b(n) = 1.0D+00
!
!  Set the matrix A.
!
  m = ml + mu + 1
  write ( *, '(a,i8)' ) '  The bandwidth of the matrix is ', m

  do j = 1, n
    a(m-1,j) = -1.0D+00
    a(m,j) = 2.0D+00
    a(m+1,j) = -1.0D+00
  end do
!
!  Call DGBFA to factor the matrix A.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call dgbfa ( a, lda, n, ml, mu, ipivot, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  DGBFA returns INFO = ', info
    return
  end if
!
!  Call DGBSL to solve the linear system.  The solution
!  is returned in B, that is, it overwrites the right hand side.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call dgbsl ( a, lda, n, ml, mu, ipivot, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 solution entries:'
  write ( *, '(a)' ) '  (All should be 1):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    if ( i <= 5 .or. n-5 < i ) then
      write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
    end if
    if ( i == 5 ) then
      write ( *, '(a)' ) '  ......  ..............'
    end if
  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests DGBFA and DGBDI.
!
!  Discussion:
!
!    Matrix A is ( 2 -1  0  0  0)
!                (-1  2 -1  0  0)
!                ( 0 -1  2 -1  0)
!                ( 0  0 -1  2 -1)
!                ( 0  0  0 -1  2)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Local Parameters:
!
!    N is the number of equations.
!
!    ML is the number of subdiagonals,
!    MU the number of superdiagonals.
!
!    LDA is the leading dimension of the array used to store the
!    matrix, which must be at least 2*ML+MU+1.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 128
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: lda = 2 * ml + mu + 1

  real ( kind = 8 ) a(lda,n_max)
  real ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n_max)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  For a banded matrix in general format,'
  write ( *, '(a)' ) '  DGBFA factors the matrix,'
  write ( *, '(a)' ) '  DGBDI computes the determinant as'
  write ( *, '(a)' ) '    det = MANTISSA * 10**EXPONENT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Find the determinant of the -1,2,-1 matrix'
  write ( *, '(a)' ) '  for N = 2, 4, 8, 16, 32, 64, 128.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (For this matrix, det ( A ) = N + 1.)'
!
!  Set the matrix A.
!
  m = ml + mu + 1
  write ( *, '(a,i8)' ) '  The bandwidth of the matrix is ', m
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N    Mantissa       Exponent'
  write ( *, '(a)' ) ' '

  n = 1

  do n_log = 1, 7

    n = 2 * n

    a(1:lda,1:n) = 0.0D+00

    do j = 1, n
      a(m-1,j) = -1.0D+00
      a(m,j) = 2.0D+00
      a(m+1,j) = -1.0D+00
    end do

    call dgbfa ( a, lda, n, ml, mu, ipivot, info )

    if ( info /= 0 ) then
      write ( *, '(a,i8)' ) '  Error!  DGBFA returns INFO = ', info
      return
    end if

    call dgbdi ( a, lda, n, ml, mu, ipivot, det )

    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, det(1), det(2)

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests DGBFA and DGBSL.
!
!  Discussion:
!
!    DGBFA and DGBSL are for general banded matrices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: ml = 25
  integer ( kind = 4 ), parameter :: mu = 25
  integer ( kind = 4 ), parameter :: lda = 2 * ml + mu + 1

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) m
  real ( kind = 8 ) temp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  For a banded matrix in general format,'
  write ( *, '(a)' ) '  DGBFA factors the matrix,'
  write ( *, '(a)' ) '  DGBSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to matrix A and right hand side B.
!
!  We want to try a problem with a significant bandwidth.
!
  m = ml + mu + 1
  write ( *, '(a,i8)' ) '  The bandwidth of the matrix is ', m

  do j = 1, n

    ilo = max ( 1, j - mu )
    ihi = min ( n, j + ml )

    temp = 0.0D+00
    do i = ilo, ihi
      a(i-j+m,j) = -1.0D+00
      temp = temp - 1.0D+00
    end do

    temp = temp + 1.0D+00
    a(m,j) = 4.0D+00 - temp
    b(j) = 4.0D+00

  end do
!
!  Factor the matrix A.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call dgbfa ( a, lda, n, ml, mu, ipivot, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  DGBFA returns INFO = ', info
    return
  end if
!
!  Call DGBSL to solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call dgbsl ( a, lda, n, ml, mu, ipivot, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 solution entries:'
  write ( *, '(a)' ) '  (All should be 1):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    if ( i <= 5 .or. n-5 < i ) then
      write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
    end if
    if ( i == 5 ) then
      write ( *, '(a)' ) '  ......  ..............'
    end if
  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 calls DGECO and DGESL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Local Parameters:
!
!    LDA defines the maximum matrix size we will use.
!
  implicit none

  integer ( kind = 4 ), parameter :: lda = 10

  real ( kind = 8 ) a(lda,lda)
  real ( kind = 8 ) b(lda)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(lda)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) n
  real ( kind = 8 ) rcond
  real ( kind = 8 ) z(lda)

  n = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  DGECO factors a general matrix and computes'
  write ( *, '(a)' ) '  its reciprocal condition number;'
  write ( *, '(a)' ) '  DGESL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = 1.0D+00
  a(1,2) = 2.0D+00
  a(1,3) = 3.0D+00

  a(2,1) = 4.0D+00
  a(2,2) = 5.0D+00
  a(2,3) = 6.0D+00

  a(3,1) = 7.0D+00
  a(3,2) = 8.0D+00
  a(3,3) = 0.0D+00
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call dgeco ( a, lda, n, ipvt, rcond, z )

  write ( *, '(a,g14.6)' ) '  The reciprocal matrix condition number = ', rcond

  if ( rcond + 1.0D+00 == 1.0D+00 ) then
    write ( *, '(a)' ) '  Error!  The matrix is nearly singular!'
    return
  end if
!
!  Set a right hand side.
!
  b(1:3) = (/ 14.0D+00, 32.0D+00, 23.0D+00 /)
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call dgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution returned by DGESL'
  write ( *, '(a)' ) '  (Should be (1,2,3))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) b(i)
  end do
!
!  A second right hand side can be solved without refactoring a.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call DGESL for a new right hand '
  write ( *, '(a)' ) '  side for the same, factored matrix.'
!
!  Set the right hand side.
!
  b(1:3) = (/ 1.0D+00, 4.0D+00, 7.0D+00 /)
!
!  Solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve a linear system.'

  job = 0
  call dgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution returned by DGESL'
  write ( *, '(a)' ) '  (should be (1,0,0))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) b(i)
  end do
!
!  The transposed problem A'*X = B can be solved by DGESL
!  as well, without any refactoring.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call DGESL for transposed problem.'
!
!  Set the right hand side.
!
  b(1:3) = (/ 6.0D+00, 6.0D+00, -3.0D+00 /)
!
!  Solve the transposed system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call DGESL to solve a transposed linear system.'

  job = 1
  call dgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution returned by DGESL'
  write ( *, '(a)' ) '  (should be (-1,0,1))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) b(i)
  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests DGEFA and DGEDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  real ( kind = 8 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  DGEFA factors a general matrix;'
  write ( *, '(a)' ) '  DGEDI computes the inverse and determinant'
  write ( *, '(a)' ) '  of a factored matrix.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = 1.0D+00
  a(1,2) = 2.0D+00
  a(1,3) = 3.0D+00

  a(2,1) = 4.0D+00
  a(2,2) = 5.0D+00
  a(2,3) = 6.0D+00

  a(3,1) = 7.0D+00
  a(3,2) = 8.0D+00
  a(3,3) = 0.0D+00
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix'

  call dgefa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) '  Error!  The matrix is nearly singular!'
    return
  end if
!
!  Get the inverse and determinant.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get the inverse and determinant'

  job = 11
  call dgedi ( a, lda, n, ipvt, det, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  The determinant = ', det(1), ' * 10 ** ', det(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The inverse matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') a(i,1:n)
  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests DGEFA and DGESL.
!
!  Discussion:
!
!    Solve A*x = b where A is a given matrix, and B a right hand side.
!
!    We will also assume that A is stored in the simplest
!    possible way.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  DGEFA factors a general matrix;'
  write ( *, '(a)' ) '  DGESL solves a factored linear system;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of equations is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = 1.0D+00
  a(1,2) = 2.0D+00
  a(1,3) = 3.0D+00

  a(2,1) = 4.0D+00
  a(2,2) = 5.0D+00
  a(2,3) = 6.0D+00

  a(3,1) = 7.0D+00
  a(3,2) = 8.0D+00
  a(3,3) = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) a(i,1:n)
  end do
!
!  Set the values of the right hand side vector B.
!
  b(1:3) = (/ 14.0D+00, 32.0D+00, 23.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The right hand side B is '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) b(i)
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix'

  call dgefa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  DGEFA returned an error flag INFO = ', info
    return
  end if
!
!  If no error occurred, now use DGESL to solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call dgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DGESL returns the solution:'
  write ( *, '(a)' ) '  (Should be (1,2,3))'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,g14.6)' ) b(i)
  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests DGEFA and DGESL.
!
!  Discussion:
!
!    In this example, we solve a relatively large linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  DGEFA factors a general matrix;'
  write ( *, '(a)' ) '  DGESL solves a factored linear system;'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to the matrix A and the right hand side B.
!
!  The problem is just an enlarged version of the
!  problem for N = 5, which is:
!
!  Matrix A is ( n -1 -1 -1 -1)    Right hand side B is  (1)
!              (-1  n -1 -1 -1)                          (1)
!              (-1 -1  n -1 -1)                          (1)
!              (-1 -1 -1  n -1)                          (1)
!              (-1 -1 -1 -1  n)                          (1)
!
!  Solution is   (1)
!                (1)
!                (1)
!                (1)
!                (1)
!
  b(1:n) = 1.0D+00

  a(1:n,1:n) = -1.0D+00
  do i = 1, n
    a(i,i) = real ( n, kind = 8 )
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call dgefa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  DGEFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call dgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last five solution entries:'
  write ( *, '(a)' ) '  (All of them should be 1.)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    if ( i <= 5 .or. n-5 < i ) then
      write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
    end if
    if ( i == 5 ) then
      write ( *, '(a)' ) '  ......  ..............'
    end if
  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests DGTSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  For a general tridiagonal matrix,'
  write ( *, '(a)' ) '  DGTSL factors and solves a linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
  write ( *, '(a)' ) ' '
!
!  Set up the linear system, by storing the values of the
!  subdiagonal, diagonal, and superdiagonal in C, D, and E,
!  and the right hand side in B.
!
  c(1) = 0.0D+00
  c(2:n) = -1.0D+00

  d(1:n) = 2.0D+00

  e(1:n-1) = -1.0D+00
  e(n) = 0.0D+00

  b(1:n-1) = 0.0D+00
  b(n) = real ( n + 1, kind = 8 )
!
!  Factor and solve the system in one step.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix and solve the system.'

  call dgtsl ( n, c, d, e, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  DGTSL returns nonzero INFO = ', info
    return
  end if
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 solution entries:'
  write ( *, '(a)' ) '  (Should be (1,2,3,4,5,...,n-1,n))'
  write ( *, '(a)' ) ' '

  do i = 1, n
    if ( i <= 5 .or. n-5 < i ) then
      write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
    end if
    if ( i == 5 ) then
      write ( *, '(a)' ) '  ......  ..............'
    end if
  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests DPBCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ), parameter :: lda = m+1

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real ( kind = 8 ) rcond
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  For a positive definite symmetric band matrix,'
  write ( *, '(a)' ) '  DPBCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the value of the superdiagonal and diagonal.
!
  a(1,1)   =  0.0D+00
  a(1,2:n) = -1.0D+00
  a(2,1:n) =  2.0D+00
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call dpbco ( a, lda, n, m, rcond, z, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Reciprocal condition  = ', rcond

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests DPBDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 128
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ), parameter :: lda = m+1

  real ( kind = 8 ) a(lda,n_max)
  real ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  For a positive definite symmetric band matrix,'
  write ( *, '(a)' ) '  DPBDI computes the determinant as'
  write ( *, '(a)' ) '    det = MANTISSA * 10**EXPONENT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Find the determinant of the -1,2,-1 matrix'
  write ( *, '(a)' ) '  for N = 2, 4, 8, 16, 32, 64, 128.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (For this matrix, det ( A ) = N + 1.)'
  write ( *, '(a)' ) ' '

  write ( *, '(a,i8)' ) '  The bandwidth of the matrix is ', 2 * m + 1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N    Mantissa       Exponent'
  write ( *, '(a)' ) ' '

  n = 1

  do n_log = 1, 7

    n = 2 * n

    a(1:lda,1:n) = 0.0D+00

    a(1,1)   =  0.0D+00
    a(1,2:n) = -1.0D+00
    a(2,1:n) =  2.0D+00

    call dpbfa ( a, lda, n, m, info )

    if ( info /= 0 ) then
      write ( *, '(a,i8)' ) '  Error!  DPBFA returns INFO = ', info
      return
    end if

    call dpbdi ( a, lda, n, m, det )

    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, det(1), det(2)

  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests DPBFA and DPBSL.
!
!  Discussion:
!
!    DPBFA and DPBSL are for a positive definite symmetric band matrix.
!
!    The problem is just an enlarged version of the
!    problem for N = 5, which is:
!
!    Matrix A is ( 2 -1  0  0  0)    right hand side B is  (1)
!                (-1  2 -1  0  0)                          (0)
!                ( 0 -1  2 -1  0)                          (0)
!                ( 0  0 -1  2 -1)                          (0)
!                ( 0  0  0 -1  2)                          (1)
!
!    solution is   (1)
!                  (1)
!                  (1)
!                  (1)
!                  (1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ), parameter :: lda = m+1

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  For a positive definite symmetric band matrix,'
  write ( *, '(a)' ) '  DPBFA computes the LU factors.'
  write ( *, '(a)' ) '  DPBSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to matrix A and right hand side B.
!
!  Set the right hand side.
!
  b(1) = 1.0D+00
  b(2:n-1) = 0.0D+00
  b(n) = 1.0D+00
!
!  Set the value of the superdiagonal and diagonal.
!
  a(1,1)   =  0.0D+00
  a(1,2:n) = -1.0D+00
  a(2,1:n) =  2.0D+00
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call dpbfa ( a, lda, n, m, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  DPBFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call dpbsl ( a, lda, n, m, b )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 solution entries:'
  write ( *, '(a)' ) '  (All should be 1):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    if ( i <= 5 .or. n-5 < i ) then
      write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
    end if
    if ( i == 5 ) then
      write ( *, '(a)' ) '  ......  ..............'
    end if
  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests DPOCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real ( kind = 8 ) rcond
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  For a positive definite symmetric matrix,'
  write ( *, '(a)' ) '  DPOCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = 2.0D+00
    if ( 1 < i ) then
      a(i,i-1) = -1.0D+00
    end if
    if ( i < n ) then
      a(i,i+1) = -1.0D+00
    end if
  end do
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call dpoco ( a, lda, n, rcond, z, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Reciprocal condition  = ', rcond

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests DPOFA and DPODI.
!
!  Discussion:
!
!    DPOFA factors a positive definite symmetric matrix,
!    and DPODI can compute the determinant or the inverse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  For a positive definite symmetric matrix,'
  write ( *, '(a)' ) '  DPOFA computes the LU factors,'
  write ( *, '(a)' ) '  DPODI computes the inverse or determinant.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = 2.0D+00
    if ( 1 < i ) then
      a(i,i-1) = -1.0D+00
    end if
    if ( i < n ) then
      a(i,i+1) = -1.0D+00
    end if
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call dpofa ( a, lda, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, DPOFA returns INFO = ', info
    return
  end if
!
!  Get the determinant and inverse.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get the determinant and inverse.'

  job = 11
  call dpodi ( a, lda, n, det, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant  = ', det(1), ' * 10 ** ', det(2)
!
!  DPODI produces only the 'upper half triangle' of the inverse,
!  which is actually symmetric.  Thus, the lower half could be
!  produced by copying from the upper half.  However, the first row
!  of A, as returned, is exactly the first row of the inverse.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First row of inverse:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) a(1,1:n)

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests DPOFA and DPOSL.
!
!  Discussion:
!
!    DPOFA factors a positive definite symmetric matrix,
!    and DPOSL can solve a factored linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  For a positive definite symmetric matrix,'
  write ( *, '(a)' ) '  DPOFA computes the LU factors.'
  write ( *, '(a)' ) '  DPOSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = 2.0D+00
    if ( 1 < i ) then
      a(i,i-1) = -1.0D+00
    end if
    if ( i < n ) then
      a(i,i+1) = -1.0D+00
    end if
  end do
!
!  Set the right hand side.
!
  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call dpofa ( a, lda, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, DPOFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call dposl ( a, lda, n, b )
!
!  Print the result.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 solution entries:'
  write ( *, '(a)' ) '  (Should be 1,2,3,4,5,...,n-1,n):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    if ( i <= 5 .or. n-5 < i ) then
      write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
    end if
    if ( i == 5 ) then
      write ( *, '(a)' ) '  ......  ..............'
    end if
  end do

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests DPPCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a((n*(n+1))/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) rcond
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  For a positive definite symmetric packed matrix,'
  write ( *, '(a)' ) '  DPPCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      if ( i == j - 1 ) then
        a(k) = -1.0D+00
      else if ( i == j ) then
        a(k) = 2.0D+00
      else
        a(k) = 0.0D+00
      end if
    end do
  end do
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition number.'

  call dppco ( a, n, rcond, z, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Reciprocal condition number = ', rcond

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests DPPFA and DPPDI.
!
!  Discussion:
!
!    DPPFA factors a packed positive definite symmetric matrix,
!    and DPPDI can compute the determinant or the inverse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  For a positive definite symmetric packed matrix,'
  write ( *, '(a)' ) '  DPPFA factors the matrix.'
  write ( *, '(a)' ) '  DPPDI computes the inverse or determinant.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      if ( i == j - 1 ) then
        a(k) = -1.0D+00
      else if ( i == j ) then
        a(k) = 2.0D+00
      else
        a(k) = 0.0D+00
      end if
    end do
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call dppfa ( a, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, DPPFA returns INFO = ', info
    return
  end if
!
!  Invert the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get the determinant and inverse.'

  job = 11
  call dppdi ( a, n, det, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant  = ', det(1), ' * 10 ** ', det(2)
!
!  DPPDI produces only the 'upper half triangle' of the inverse,
!  which is actually symmetric.  Thus, the lower half could be
!  produced by copying from the upper half.  However, the first row
!  of A, as returned, is exactly the first row of the inverse.
!
  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      b(i,j) = a(k)
      b(j,i) = a(k)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Inverse:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,5g14.6)' ) b(i,1:n)
  end do

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests DPPFA and DPPSL.
!
!  Discussion:
!
!    DPOFA factors a positive definite symmetric matrix,
!    and DPOSL can solve a factored linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  For a positive definite symmetric packed matrix,'
  write ( *, '(a)' ) '  DPPFA factors the matrix.'
  write ( *, '(a)' ) '  DPPSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  b(1:n) = 0.0D+00

  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      if ( i == j - 1 ) then
        a(k) = -1.0D+00
        b(i) = b(i) + a(k) * x(j)
        b(j) = b(j) + a(k) * x(i)
      else if ( i == j ) then
        a(k) = 2.0D+00
        b(i) = b(i) + a(k) * x(i)
      else
        a(k) = 0.0D+00
      end if
    end do
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call dppfa ( a, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, DPPFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call dppsl ( a, n, b )
!
!  Print the result.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 solution entries:'
  write ( *, '(a)' ) '  (Should be 1,2,3,4,5,...,n-1,n):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    if ( i <= 5 .or. n-5 < i ) then
      write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
    end if
    if ( i == 5 ) then
      write ( *, '(a)' ) '  ......  ..............'
    end if
  end do

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests DPTSL.
!
!  Discussion:
!
!    DPTSL factors and solves a positive definite symmetric tridiagonal system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  For a positive definite symmetric tridiagonal matrix,'
  write ( *, '(a)' ) '  DPTSL factors and solves a linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  b(1:n) = 0.0D+00
  d(1:n) = 2.0D+00
  e(1:n-1) = -1.0D+00
  e(n) = 0.0D+00

  do i = 1, n

    if ( 1 < i ) then
      b(i) = b(i) + e(i-1) * x(i-1)
    end if
    b(i) = b(i) + d(i) * x(i)
    if ( i < n ) then
      b(i) = b(i) + e(i) * x(i+1)
    end if

  end do
!
!  Factor and solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix and solve the system.'

  call dptsl ( n, d, e, b )
!
!  Print the result.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 solution entries:'
  write ( *, '(a)' ) '  (Should be 1,2,3,4,5,...,n-1,n):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    if ( i <= 5 .or. n-5 < i ) then
      write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
    end if
    if ( i == 5 ) then
      write ( *, '(a)' ) '  ......  ..............'
    end if
  end do

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests DQRDC and DQRSL.
!
!  Discussion:
!
!    DQRDC and DQRSL compute the QR factorization, and use it
!    to solve linear systems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: p = 3
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,p)
  real ( kind = 8 ) b(lda,p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(p)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real ( kind = 8 ) q(n,n)
  real ( kind = 8 ) qraux(p)
  real ( kind = 8 ) qty(n)
  real ( kind = 8 ) qy(n)
  real ( kind = 8 ) r(n,p)
  real ( kind = 8 ) rsd(n)
  real ( kind = 8 ) work(p)
  real ( kind = 8 ) xb(n)
  real ( kind = 8 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  For a general matrix,'
  write ( *, '(a)' ) '  DQRDC computes the QR decomposition of a '
  write ( *, '(a)' ) '  matrix, but does not return Q and R explicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Show how Q and R can be recovered using DQRSL.'
!
!  Set the matrix A.
!
  a(1,1) = 1.0D+00
  a(2,1) = 1.0D+00
  a(3,1) = 0.0D+00

  a(1,2) = 1.0D+00
  a(2,2) = 0.0D+00
  a(3,2) = 1.0D+00

  a(1,3) = 0.0D+00
  a(2,3) = 1.0D+00
  a(3,3) = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The original matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) a(i,1:p)
  end do
!
!  Decompose the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Decompose the matrix.'

  job = 0
  ipvt(1:p) = 0

  call dqrdc ( a, lda, n, p, qraux, ipvt, work, job )
!
!  Print out what DQRDC has stored in A...
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The packed matrix A which describes Q and R:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) a(i,1:p)
  end do
!
!  ...and in QRAUX.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The QRAUX vector, containing some additional'
  write ( *, '(a)' ) '  information defining Q:'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,5g14.6)' )  qraux(1:n)
!
!  Print out the resulting R factor.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The R factor:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    do j = 1, p

      if ( j < i ) then
        r(i,j) = 0.0D+00
      else
        r(i,j) = a(i,j)
      end if

    end do

    write ( *, '(2x,5g14.6)' ) r(i,1:p)

  end do
!
!  Call DQRSL to extract the information about the Q matrix.
!  We do this, essentially, by asking DQRSL to tell us the
!  value of Q*Y, where Y is a column of the identity matrix.
!
  job = 10000

  do i = 1, n
!
!  Set the vector Y.
!
    y(1:n) = 0.0D+00

    y(i) = 1.0D+00
!
!  Ask DQRSL to tell us what Q*Y is.
!
    call dqrsl ( a, lda, n, p, qraux, y, qy, qty, b, rsd, xb, job, info )

    if ( info /= 0 ) then
      write ( *, '(a,i8)' ) '  Error!  DQRSL returns INFO = ', info
      return
    end if
!
!  Copy QY into the appropriate column of Q.
!
    q(1:n,i) = qy(1:n)

  end do
!
!  Now print out the Q matrix we have extracted.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Q factor:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) q(i,1:n)
  end do
!
!  Compute Q*R to verify that it equals A.
!
  b(1:n,1:p) = matmul ( q(1:n,1:n), r(1:n,1:p) )
!
!  Print the result.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product Q * R:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) b(i,1:p)
  end do

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests DSICO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  real ( kind = 8 ) rcond
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  For a symmetric indefinite matrix,'
  write ( *, '(a)' ) '  DSICO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to the matrix A.
!
  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = 2.0D+00
    if ( i < n ) then
      a(i,i+1) = -1.0D+00
    end if
  end do
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call dsico ( a, lda, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition = ', rcond

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests DSIFA and DSISL.
!
!  Discussion:
!
!    DSIFA and DSISL are for symmetric indefinite matrices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  For a symmetric indefinite matrix,'
  write ( *, '(a)' ) '  DSIFA factors the matrix,'
  write ( *, '(a)' ) '  DSISL solves a factored linear system,'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to the matrix A and the right hand side B.
!
  b(1:n-1) = 0.0D+00
  b(n)= real ( n + 1, kind = 8 )

  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = 2.0D+00
    if ( i < n ) then
      a(i,i+1) = -1.0D+00
    end if
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call dsifa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  DSIFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call dsisl ( a, lda, n, ipvt, b )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 solution entries:'
  write ( *, '(a)' ) '  (Should be (1,2,3,4,5,...,n-1,n))'
  write ( *, '(a)' ) ' '

  do i = 1, n
    if ( i <= 5 .or. n-5 < i ) then
      write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
    end if
    if ( i == 5 ) then
      write ( *, '(a)' ) '  ......  ..............'
    end if
  end do

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26 tests DSPCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100

  real ( kind = 8 ) a((n*(n+1))/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) rcond
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  For a symmetric indefinite packed matrix,'
  write ( *, '(a)' ) '  DSPCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to the matrix A.
!
  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      if ( i == j ) then
        a(k) = 2.0D+00
      else if ( j == i+1 ) then
        a(k) = -1.0D+00
      else
        a(k) = 0.0D+00
      end if
    end do
  end do
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call dspco ( a, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition = ', rcond

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests DSPFA and DSPSL.
!
!  Discussion:
!
!    DSPFA and DSPSL are for packed symmetric indefinite matrices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  For a symmetric indefinite packed matrix,'
  write ( *, '(a)' ) '  DSPFA factors the matrix,'
  write ( *, '(a)' ) '  DSPSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to the matrix A and the right hand side B.
!
  b(1:n-1) = 0.0D+00
  b(n)= real ( n + 1, kind = 8 )

  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      if ( i == j ) then
        a(k) = 2.0D+00
      else if ( j == i+1 ) then
        a(k) = -1.0D+00
      else
        a(k) = 0.0D+00
      end if
    end do
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call dspfa ( a, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  DSPFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call dspsl ( a, n, ipvt, b )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 solution entries:'
  write ( *, '(a)' ) '  (Should be (1,2,3,4,5,...,n-1,n))'
  write ( *, '(a)' ) ' '

  do i = 1, n
    if ( i <= 5 .or. n-5 < i ) then
      write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
    end if
    if ( i == 5 ) then
      write ( *, '(a)' ) '  ......  ..............'
    end if
  end do

  return
end
subroutine test28 ( )

!*****************************************************************************80
!
!! TEST28 tests DSVDC.
!
!  Discussion:
!
!    DSVDC computes the singular value decomposition:
!
!      A = U * S * V'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 6
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) e(max(m+1,n))
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) job
  real ( kind = 8 ) s(max(m+1,n))
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sigma(m,n)
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) work(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  For an MxN matrix A in general storage,'
  write ( *, '(a)' ) '  DSVDC computes the singular value decomposition:'
  write ( *, '(a)' ) '    A = U * S * V'''
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set A.
!
  seed = 123456789

  call r8mat_uniform_01 ( m, n, seed, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,7f10.4)' ) a(i,1:n)
  end do
!
!  Decompose the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Decompose the matrix.'

  job = 11
  lda = m
  ldu = m
  ldv = n

  call dsvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Warning:'
    write ( *, '(a,i8)' ) '  DSVDC returned nonzero INFO = ', info
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Singular values:'
  write ( *, '(a)' ) ' '

  do i = 1, min ( m, n )
    write ( *, '(2x,i4,2x,g14.6)' ) i, s(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Left Singular Vector Matrix U:'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,7f10.4)' ) u(i,1:m)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Right Singular Vector Matrix V:'
  write ( *, '(a)' ) ' '

  do i = 1,  n
    write ( *, '(2x,7f10.4)' ) v(i,1:n)
  end do

  sigma(1:m,1:n) = 0.0D+00
  do i = 1, min ( m, n )
    sigma(i,i) = s(i)
  end do

  a(1:m,1:n) = matmul ( u(1:m,1:m), &
    matmul ( sigma(1:m,1:n), transpose ( v(1:n,1:n) ) ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product U * S * V'' (should equal A):'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,7f10.4)' ) a(i,1:n)
  end do

  return
end
subroutine test29 ( )

!*****************************************************************************80
!
!! TEST29 tests DTRCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real ( kind = 8 ) rcond
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  For a triangular matrix,'
  write ( *, '(a)' ) '  DTRCO computes the LU factors and'
  write ( *, '(a)' ) '  computes its reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Lower triangular matrix A.
!
  call r8mat_uniform_01 ( n, n, seed, a )

  do i = 1, n
    do j = i+1, n
      a(i,j) = 0.0D+00
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lower triangular matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') a(i,1:n)
  end do

  job = 0
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition:'

  call dtrco ( a, lda, n, rcond, z, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The reciprocal condition number = ', rcond
!
!  Upper triangular matrix A.
!
  call r8mat_uniform_01 ( n, n, seed, a )

  do i = 1, n
    do j = 1, i - 1
      a(i,j) = 0.0D+00
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Upper triangular matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') a(i,1:n)
  end do
!
!  Estimate the condition.
!
  job = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition:'

  call dtrco ( a, lda, n, rcond, z, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The reciprocal condition number = ', rcond

  return
end
subroutine test30 ( )

!*****************************************************************************80
!
!! TEST30 tests DTRDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  For a triangular matrix,'
  write ( *, '(a)' ) '  DTRDI computes the determinant or inverse.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Lower triangular matrix A.
!
  call r8mat_uniform_01 ( n, n, seed, a )

  do i = 1, n
    do j = i+1, n
      a(i,j) = 0.0D+00
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lower triangular matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') a(i,1:n)
  end do

  job = 110

  call dtrdi ( a, lda, n, det, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  The determinant = ', det(1), ' * 10 ** ', det(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The inverse matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') a(i,1:n)
  end do
!
!  Upper triangular matrix A.
!
  call r8mat_uniform_01 ( n, n, seed, a )

  do i = 1, n
    do j = 1, i - 1
      a(i,j) = 0.0D+00
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Upper triangular matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') a(i,1:n)
  end do

  job = 111

  call dtrdi ( a, lda, n, det, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  The determinant = ', det(1), ' * 10 ** ', det(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The inverse matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') a(i,1:n)
  end do

  return
end
subroutine test31 ( )

!*****************************************************************************80
!
!! TEST31 tests DTRSL.
!
!  Discussion:
!
!    DTRSL solves triangular linear systems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  For a triangular matrix,'
  write ( *, '(a)' ) '  DTRSL solves a linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Lower triangular matrix A.
!
  call r8mat_uniform_01 ( n, n, seed, a )

  do i = 1, n
    do j = i+1, n
      a(i,j) = 0.0D+00
    end do
  end do

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For a lower triangular matrix A,'
  write ( *, '(a)' ) '  solve A * x = b'

  job = 00

  call dtrsl ( a, lda, n, b, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The solution (should be 1,2,3,4,5):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
  end do

  b(1:n) = matmul ( transpose ( a(1:n,1:n) ), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For a lower triangular matrix A,'
  write ( *, '(a)' ) '  solve A'' * x = b'

  job = 10

  call dtrsl ( a, lda, n, b, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The solution (should be 1,2,3,4,5):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
  end do
!
!  Upper triangular matrix A.
!
  call r8mat_uniform_01 ( n, n, seed, a )

  do i = 1, n
    do j = 1, i - 1
      a(i,j) = 0.0D+00
    end do
  end do

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For an upper triangular matrix A,'
  write ( *, '(a)' ) '  solve A * x = b'

  job = 01

  call dtrsl ( a, lda, n, b, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The solution (should be 1,2,3,4,5):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
  end do

  b(1:n) = matmul ( transpose ( a(1:n,1:n) ), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For an upper triangular matrix A,'
  write ( *, '(a)' ) '  solve A'' * x = b'

  job = 11

  call dtrsl ( a, lda, n, b, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The solution (should be 1,2,3,4,5):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
  end do

  return
end
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
!    use I4_HUGE() and HUGE interchangeably.
!
!    However, when using the G95, the values returned by HUGE were
!    not equal to 2147483647, apparently, and were causing severe
!    and obscure errors in my random number generator, which needs to
!    add I4_HUGE to the seed whenever the seed is negative.  So I
!    am backing away from invoking HUGE, whereas I4_HUGE is under
!    my control.
!
!    Explanation: because under G95 the default integer type is 64 bits!
!    So HUGE ( 1 ) = a very very huge integer indeed, whereas
!    I4_HUGE ( ) = the same old 32 bit big value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer I4_HUGE, a "huge" I4.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

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
!    For now, the input quantity SEED is an integer variable.
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
!    Input, integer M, N, the number of rows and columns
!    in the array.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
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
        seed = seed + i4_huge ( )
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
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
!    26 February 2005
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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
