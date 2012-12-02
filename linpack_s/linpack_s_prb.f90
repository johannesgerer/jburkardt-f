program main

!*****************************************************************************80
!
!! MAIN is the main program for LINPACK_S_PRB.
!
!  Discussion:
!
!    LINPACK_S_PRB calls the single precision real LINPACK test routines.
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

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINPACK_S_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LINPACK_S library.'

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
  write ( *, '(a)' ) 'LINPACK_S_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests SCHDC.
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

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real ( kind = 4 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  SCHDC computes the Cholesky decomposition.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of equations is N = ', n
!
!  Set the values of the matrix A.
!
  a(1:n,1:n) = 0.0E+00

  do i = 1, n
    a(i,i) = 2.0E+00
  end do

  do i = 1, n-1
    a(i,i+1) = -1.0E+00
  end do

  do i = 2, n
    a(i-1,i) = -1.0E+00
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

  call schdc ( a, lda, n, work, ipvt, job, info )

  if ( info /= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  SCHDC returned INFO = ', info
    write ( *, '(a)' ) '  This means the matrix is not positive definite.'
    return
  end if
!
!  Zero out the lower diagonal.
!
  do i = 2, n
    do j = 1, i-1
      a(i,j) = 0.0E+00
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
!! TEST02 tests SCHEX.
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

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n
  integer ( kind = 4 ), parameter :: nz = 1

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 4 ) s(n)
  integer ( kind = 4 ) seed
  real ( kind = 4 ) work(n)
  real ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  SCHEX can shift columns in a Cholesky factorization.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of equations is N = ', n
!
!  Set the values of the matrix A.
!
  a(1:n,1:n) = 0.0E+00

  do i = 1, n
    a(i,i) = 2.0E+00
  end do

  do i = 1, n-1
    a(i,i+1) = -1.0E+00
  end do

  do i = 2, n
    a(i-1,i) = -1.0E+00
  end do

  do i = 1, n
    z(i) = real ( i )
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

  call schdc ( a, lda, n, work, ipvt, job, info )

  if ( info /= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  SCHDC returned INFO = ', info
    write ( *, '(a)' ) '  This means the matrix is not positive definite.'
    return
  end if
!
!  Zero out the lower diagonal.
!
  do i = 2, n
    do j = 1, i-1
      a(i,j) = 0.0E+00
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
  call schex ( a, lda, n, k, l, z, n, nz, c, s, job )
!
!  Left circular shift columns K+1 through L.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  Left circular shift columns K+1 = ', k+1, &
    ' through L = ', l

  job = 2
  call schex ( a, lda, n, k+1, l, z, n, nz, c, s, job )
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
!! TEST03 tests SCHUD and STRSL.
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

  integer ( kind = 4 ), parameter :: p = 20
  integer ( kind = 4 ), parameter :: ldr = p
  integer ( kind = 4 ), parameter :: nz = 1

  real ( kind = 4 ) b(p)
  real ( kind = 4 ) beta(p)
  real ( kind = 4 ) c(p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real ( kind = 4 ) r(ldr,p)
  real ( kind = 4 ) rho(nz)
  real ( kind = 4 ) row(p)
  real ( kind = 4 ) s(p)
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x(p)
  real ( kind = 4 ) y(nz)
  real ( kind = 4 ) z(p,nz)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  SCHUD updates a Cholesky decomposition.'
  write ( *, '(a)' ) '  STRSL solves a triangular linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we use SCHUD to solve a'
  write ( *, '(a)' ) '  least squares problem R * b = z.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of equations is P = ', p
!
!  Initialize.
!
  r(1:p,1:p) = 0.0E+00
  z(1:p,1:nz) = 0.0E+00
  do i = 1, p
    x(i) = real ( i )
  end do
!
!  Use SCHUD to form R, Z and RHO by adding X and Y a row at a time.
!  X is a row of the least squares matrix and Y the right hand side.
!
  seed = 123456789

  do i = 1, p
    call r4mat_uniform_01 ( 1, p, seed, row )
    y(1) = dot_product ( row(1:p), x(1:p) )
    rho(1) = 0.0E+00
    call schud ( r, ldr, p, row, z, p, nz, y, rho, c, s )
  end do
!
!  Generate the least squares solution, b = inverse ( R ) * Z.
!
  do j = 1, nz

    b(1:p) = z(1:p,j)
    job = 01

    call strsl ( r, ldr, p, b, job, info )

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
!! TEST04 tests SGBCO.
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

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: lda = 2 * ml + mu + 1

  real ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 4 ) rcond
  real ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and general band storage (GB),'
  write ( *, '(a)' ) '  SGBCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  m = ml + mu + 1
  write ( *, '(a,i8)' ) '  The bandwidth of the matrix is ', m

  do j = 1, n
    a(m-1,j) = -1.0E+00
    a(m,j) = 2.0E+00
    a(m+1,j) = -1.0E+00
  end do
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call sgbco ( a, lda, n, ml, mu, ipivot, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition = ', rcond

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests SGBFA and SGBSL.
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
!    03 May 2007
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

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and general band storage (GB),'
  write ( *, '(a)' ) '  SGBFA factors the matrix,'
  write ( *, '(a)' ) '  SGBSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the right hand side B.
!
  b(1) = 1.0E+00
  b(2:n-1) = 0.0E+00
  b(n) = 1.0E+00
!
!  Set the matrix A.
!
  m = ml + mu + 1
  write ( *, '(a,i8)' ) '  The bandwidth of the matrix is ', m

  do j = 1, n
    a(m-1,j) = -1.0E+00
    a(m,j) = 2.0E+00
    a(m+1,j) = -1.0E+00
  end do
!
!  Call SGBFA to factor the matrix A.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call sgbfa ( a, lda, n, ml, mu, ipivot, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  SGBFA returns INFO = ', info
    return
  end if
!
!  Call SGBSL to solve the linear system.  The solution
!  is returned in B, that is, it overwrites the right hand side.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call sgbsl ( a, lda, n, ml, mu, ipivot, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 entries of the solution:'
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
!! TEST06 tests SGBFA and SGBDI.
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

  real ( kind = 4 ) a(lda,n_max)
  real ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n_max)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and general band storage (GB),'
  write ( *, '(a)' ) '  SGBFA factors the matrix,'
  write ( *, '(a)' ) '  SGBDI computes the determinant as'
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

    a(1:lda,1:n) = 0.0E+00

    do j = 1, n
      a(m-1,j) = -1.0E+00
      a(m,  j) =  2.0E+00
      a(m+1,j) = -1.0E+00
    end do

    call sgbfa ( a, lda, n, ml, mu, ipivot, info )

    if ( info /= 0 ) then
      write ( *, '(a,i8)' ) '  Error!  SGBFA returns INFO = ', info
      return
    end if

    call sgbdi ( a, lda, n, ml, mu, ipivot, det )

    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, det(1), det(2)

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests SGBFA and SGBSL.
!
!  Discussion:
!
!    SGBFA and SGBSL are for general banded matrices.
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

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: ml = 25
  integer ( kind = 4 ), parameter :: mu = 25
  integer ( kind = 4 ), parameter :: lda = 2 * ml + mu + 1

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) m
  real ( kind = 4 ) temp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and general band storage (GB),'
  write ( *, '(a)' ) '  SGBFA factors the matrix,'
  write ( *, '(a)' ) '  SGBSL solves a factored linear system.'
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

    temp = 0.0E+00
    do i = ilo, ihi
      a(i-j+m,j) = -1.0E+00
      temp = temp - 1.0E+00
    end do

    temp = temp + 1.0E+00
    a(m,j) = 4.0E+00 - temp
    b(j) = 4.0E+00

  end do
!
!  Factor the matrix A.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call sgbfa ( a, lda, n, ml, mu, ipivot, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  SGBFA returns INFO = ', info
    return
  end if
!
!  Call SGBSL to solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call sgbsl ( a, lda, n, ml, mu, ipivot, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 entries of the solution:'
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
!! TEST08 calls SGECO and SGESL.
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
!  Local Parameters:
!
!    LDA defines the maximum matrix size we will use.
!
  implicit none

  integer ( kind = 4 ), parameter :: lda = 10

  real ( kind = 4 ) a(lda,lda)
  real ( kind = 4 ) b(lda)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(lda)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) n
  real ( kind = 4 ) rcond
  real ( kind = 4 ) z(lda)

  n = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and general matrix storage (GE),'
  write ( *, '(a)' ) '  SGECO factors the matrix and computes'
  write ( *, '(a)' ) '  its reciprocal condition number;'
  write ( *, '(a)' ) '  SGESL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = 1.0E+00
  a(1,2) = 2.0E+00
  a(1,3) = 3.0E+00

  a(2,1) = 4.0E+00
  a(2,2) = 5.0E+00
  a(2,3) = 6.0E+00

  a(3,1) = 7.0E+00
  a(3,2) = 8.0E+00
  a(3,3) = 0.0E+00
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call sgeco ( a, lda, n, ipvt, rcond, z )

  write ( *, '(a,g14.6)' ) '  The reciprocal matrix condition number = ', rcond

  if ( rcond + 1.0E+00 == 1.0E+00 ) then
    write ( *, '(a)' ) '  Error!  The matrix is nearly singular!'
    return
  end if
!
!  Set a right hand side.
!
  b(1:3) = (/ 14.0E+00, 32.0E+00, 23.0E+00 /)
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call sgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution returned by SGESL'
  write ( *, '(a)' ) '  (Should be (1,2,3))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) b(i)
  end do
!
!  A second right hand side can be solved without refactoring a.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call SGESL for a new right hand '
  write ( *, '(a)' ) '  side for the same, factored matrix.'
!
!  Set the right hand side.
!
  b(1:3) = (/ 1.0E+00, 4.0E+00, 7.0E+00 /)
!
!  Solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve a linear system.'

  job = 0
  call sgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution returned by SGESL'
  write ( *, '(a)' ) '  (should be (1,0,0))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) b(i)
  end do
!
!  The transposed problem A'*X = B can be solved by SGESL
!  as well, without any refactoring.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call SGESL for transposed problem.'
!
!  Set the right hand side.
!
  b(1:3) = (/ 6.0E+00, 6.0E+00, -3.0E+00 /)
!
!  Solve the transposed system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call SGESL to solve a transposed linear system.'

  job = 1
  call sgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution returned by SGESL'
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
!! TEST09 tests SGEFA and SGEDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(n,n)
  real ( kind = 4 ) a_save(n,n)
  real ( kind = 4 ) c(n,n)
  real ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  real ( kind = 4 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and general matrix storage (GE),'
  write ( *, '(a)' ) '  SGEFA factors the matrix;'
  write ( *, '(a)' ) '  SGEDI computes the inverse and determinant'
  write ( *, '(a)' ) '  of a factored matrix.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = 1.0E+00
  a(1,2) = 2.0E+00
  a(1,3) = 3.0E+00

  a(2,1) = 4.0E+00
  a(2,2) = 5.0E+00
  a(2,3) = 6.0E+00

  a(3,1) = 7.0E+00
  a(3,2) = 8.0E+00
  a(3,3) = 0.0E+00

  a_save(1:n,1:n) = a(1:n,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') a(i,1:n)
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix'

  call sgefa ( a, lda, n, ipvt, info )

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
  call sgedi ( a, lda, n, ipvt, det, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  The determinant = ', det(1), ' * 10 ** ', det(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The inverse matrix inverse(A):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') a(i,1:n)
  end do

  c(1:n,1:n) = matmul ( a(1:n,1:n), a_save(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product inverse(A) * A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') c(i,1:n)
  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests SGEFA and SGESL.
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
!    03 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and general matrix storage (GE),'
  write ( *, '(a)' ) '  SGEFA factors the matrix;'
  write ( *, '(a)' ) '  SGESL solves a factored linear system;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of equations is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = 1.0E+00
  a(1,2) = 2.0E+00
  a(1,3) = 3.0E+00

  a(2,1) = 4.0E+00
  a(2,2) = 5.0E+00
  a(2,3) = 6.0E+00

  a(3,1) = 7.0E+00
  a(3,2) = 8.0E+00
  a(3,3) = 0.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)' ) a(i,1:n)
  end do
!
!  Set the values of the right hand side vector B.
!
  b(1:3) = (/ 14.0E+00, 32.0E+00, 23.0E+00 /)

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

  call sgefa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  SGEFA returned an error flag INFO = ', info
    return
  end if
!
!  If no error occurred, now use SGESL to solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call sgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SGESL returns the solution:'
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
!! TEST11 tests SGEFA and SGESL.
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
!    03 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and general matrix storage (GE),'
  write ( *, '(a)' ) '  SGEFA factors a general matrix;'
  write ( *, '(a)' ) '  SGESL solves a factored linear system;'
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
  b(1:n) = 1.0E+00

  a(1:n,1:n) = -1.0E+00
  do i = 1, n
    a(i,i) = real ( n )
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call sgefa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  SGEFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  job = 0
  call sgesl ( a, lda, n, ipvt, b, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last five entries of the solution:'
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
!! TEST12 tests SGTSL.
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

  integer ( kind = 4 ), parameter :: n = 100

  real ( kind = 4 ) b(n)
  real ( kind = 4 ) c(n)
  real ( kind = 4 ) d(n)
  real ( kind = 4 ) e(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and general tridiagonal storage (GT),'
  write ( *, '(a)' ) '  SGTSL factors and solves a linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
  write ( *, '(a)' ) ' '
!
!  Set up the linear system, by storing the values of the
!  subdiagonal, diagonal, and superdiagonal in C, D, and E,
!  and the right hand side in B.
!
  c(1) = 0.0E+00
  c(2:n) = -1.0E+00

  d(1:n) = 2.0E+00

  e(1:n-1) = -1.0E+00
  e(n) = 0.0E+00

  b(1:n-1) = 0.0E+00
  b(n) = real ( n + 1 )
!
!  Factor and solve the system in one step.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix and solve the system.'

  call sgtsl ( n, c, d, e, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  SGTSL returns nonzero INFO = ', info
    return
  end if
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 entries of solution:'
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
!! TEST13 tests SPBCO.
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

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ), parameter :: lda = m+1

  real ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real ( kind = 4 ) rcond
  real ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and positive definite symmetric band storage (PB),'
  write ( *, '(a)' ) '  SPBCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the value of the superdiagonal and diagonal.
!
  a(1,1)   =  0.0E+00
  a(1,2:n) = -1.0E+00
  a(2,1:n) =  2.0E+00
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call spbco ( a, lda, n, m, rcond, z, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Reciprocal condition  = ', rcond

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests SPBDI.
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

  integer ( kind = 4 ), parameter :: n_max = 128
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ), parameter :: lda = m+1

  real ( kind = 4 ) a(lda,n_max)
  real ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and positive definite symmetric band storage (PB),'
  write ( *, '(a)' ) '  SPBDI computes the determinant as'
  write ( *, '(a)' ) '    det = MANTISSA * 10**EXPONENT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Find the determinant of the -1,2,-1 matrix'
  write ( *, '(a)' ) '  for N = 2, 4, 8, 16, 32, 64, 128.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (For this matrix, det ( A ) = N + 1.)'
  write ( *, '(a)' ) ' '
!
!  Set the number of nonzero diagonals.
!
  write ( *, '(a,i8)' ) '  The bandwidth of the matrix is ', 2 * m + 1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N    Mantissa       Exponent'
  write ( *, '(a)' ) ' '

  n = 1

  do n_log = 1, 7

    n = 2 * n

    a(1:lda,1:n) = 0.0E+00

    a(1,1)   =  0.0E+00
    a(1,2:n) = -1.0E+00
    a(2,1:n) =  2.0E+00

    call spbfa ( a, lda, n, m, info )

    if ( info /= 0 ) then
      write ( *, '(a,i8)' ) '  Error!  SPBFA returns INFO = ', info
      return
    end if

    call spbdi ( a, lda, n, m, det )

    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, det(1), det(2)

  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests SPBFA and SPBSL.
!
!  Discussion:
!
!    SPBFA and SPBSL are for a positive definite symmetric band matrix.
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

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ), parameter :: lda = m+1

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and positive definite symmetric band storage (PB),'
  write ( *, '(a)' ) '  SPBFA computes the LU factors.'
  write ( *, '(a)' ) '  SPBSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to matrix A and right hand side B.
!
!  The problem is just an enlarged version of the
!  problem for N = 5, which is:
!
!  Matrix A is ( 2 -1  0  0  0)    right hand side B is  (1)
!              (-1  2 -1  0  0)                          (0)
!              ( 0 -1  2 -1  0)                          (0)
!              ( 0  0 -1  2 -1)                          (0)
!              ( 0  0  0 -1  2)                          (1)
!
!
!  solution is   (1)
!                (1)
!                (1)
!                (1)
!                (1)
!
!  Set the right hand side.
!
  b(1) = 1.0E+00
  b(2:n-1) = 0.0E+00
  b(n) = 1.0E+00
!
!  Set the value of the superdiagonal and diagonal.
!
  a(1,1)   =  0.0E+00
  a(1,2:n) = -1.0E+00
  a(2,1:n) =  2.0E+00
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call spbfa ( a, lda, n, m, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  SPBFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call spbsl ( a, lda, n, m, b )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 entries of the solution:'
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
!! TEST16 tests SPOCO.
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

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real ( kind = 4 ) rcond
  real ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and positive definite symmetric storage (PO),'
  write ( *, '(a)' ) '  SPOCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  a(1:n,1:n) = 0.0E+00

  do i = 1, n
    a(i,i) = 2.0E+00
    if ( 1 < i ) then
      a(i,i-1) = -1.0E+00
    end if
    if ( i < n ) then
      a(i,i+1) = -1.0E+00
    end if
  end do
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call spoco ( a, lda, n, rcond, z, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Reciprocal condition  = ', rcond

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests SPOFA and SPODI.
!
!  Discussion:
!
!    SPOFA factors a positive definite symmetric matrix,
!    and SPODI can compute the determinant or the inverse.
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

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and positive definite symmetric storage (PO),'
  write ( *, '(a)' ) '  SPOFA computes the LU factors,'
  write ( *, '(a)' ) '  SPODI computes the inverse or determinant.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  a(1:n,1:n) = 0.0E+00

  do i = 1, n
    a(i,i) = 2.0E+00
    if ( 1 < i ) then
      a(i,i-1) = -1.0E+00
    end if
    if ( i < n ) then
      a(i,i+1) = -1.0E+00
    end if
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call spofa ( a, lda, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, SPOFA returns INFO = ', info
    return
  end if
!
!  Get the determinant and inverse.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get the determinant and inverse.'

  job = 11
  call spodi ( a, lda, n, det, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant  = ', det(1), ' * 10 ** ', det(2)
!
!  SPODI produces only the 'upper half triangle' of the inverse,
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
!! TEST18 tests SPOFA and SPOSL.
!
!  Discussion:
!
!    SPOFA factors a positive definite symmetric matrix,
!    and SPOSL can solve a factored linear system.
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

  integer ( kind = 4 ), parameter :: n = 20
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and positive definite symmetric storage (PO),'
  write ( *, '(a)' ) '  SPOFA computes the LU factors.'
  write ( *, '(a)' ) '  SPOSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  a(1:n,1:n) = 0.0E+00

  do i = 1, n
    a(i,i) = 2.0E+00
    if ( 1 < i ) then
      a(i,i-1) = -1.0E+00
    end if
    if ( i < n ) then
      a(i,i+1) = -1.0E+00
    end if
  end do
!
!  Set the right hand side.
!
  do i = 1, n
    x(i) = real ( i )
  end do

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call spofa ( a, lda, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, SPOFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call sposl ( a, lda, n, b )
!
!  Print the result.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 entries of the solution:'
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
!! TEST19 tests SPPCO.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 4 ) a((n*(n+1))/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) rcond
  real ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and positive definite symmetric packed storage (PP),'
  write ( *, '(a)' ) '  SPPCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      if ( i == j - 1 ) then
        a(k) = -1.0E+00
      else if ( i == j ) then
        a(k) = 2.0E+00
      else
        a(k) = 0.0E+00
      end if
    end do
  end do
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition number.'

  call sppco ( a, n, rcond, z, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Reciprocal condition number = ', rcond

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests SPPFA and SPPDI.
!
!  Discussion:
!
!    SPPFA factors a packed positive definite symmetric matrix,
!    and SPPDI can compute the determinant or the inverse.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 4 ) a((n*(n+1))/2)
  real ( kind = 4 ) b(n,n)
  real ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and positive definite symmetric packed storage (PP),'
  write ( *, '(a)' ) '  SPPFA factors the matrix.'
  write ( *, '(a)' ) '  SPPDI computes the inverse or determinant.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      if ( i == j - 1 ) then
        a(k) = -1.0E+00
      else if ( i == j ) then
        a(k) = 2.0E+00
      else
        a(k) = 0.0E+00
      end if
    end do
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call sppfa ( a, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, SPPFA returns INFO = ', info
    return
  end if
!
!  Invert the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get the determinant and inverse.'

  job = 11
  call sppdi ( a, n, det, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant  = ', det(1), ' * 10 ** ', det(2)
!
!  SPPDI produces only the 'upper half triangle' of the inverse,
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
!! TEST21 tests SPPFA and SPPSL.
!
!  Discussion:
!
!    SPOFA factors a positive definite symmetric matrix,
!    and SPOSL can solve a factored linear system.
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

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 4 ) a((n*(n+1))/2)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and positive definite symmetric packed storage (PP),'
  write ( *, '(a)' ) '  SPPFA factors the matrix.'
  write ( *, '(a)' ) '  SPPSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  do i = 1, n
    x(i) = real ( i )
  end do

  b(1:n) = 0.0E+00

  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      if ( i == j - 1 ) then
        a(k) = -1.0E+00
        b(i) = b(i) + a(k) * x(j)
        b(j) = b(j) + a(k) * x(i)
      else if ( i == j ) then
        a(k) = 2.0E+00
        b(i) = b(i) + a(k) * x(i)
      else
        a(k) = 0.0E+00
      end if
    end do
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call sppfa ( a, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, SPPFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call sppsl ( a, n, b )
!
!  Print the result.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 entries of the solution:'
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
!! TEST22 tests SPTSL.
!
!  Discussion:
!
!    SPTSL factors and solves a positive definite symmetric tridiagonal system.
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

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 4 ) b(n)
  real ( kind = 4 ) d(n)
  real ( kind = 4 ) e(n)
  integer ( kind = 4 ) i
  real ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and positive definite symmetric tridiagonal (PT),'
  write ( *, '(a)' ) '  SPTSL factors and solves a linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the matrix A.
!
  do i = 1, n
    x(i) = real ( i )
  end do

  b(1:n) = 0.0E+00
  d(1:n) = 2.0E+00
  e(1:n-1) = -1.0E+00
  e(n) = 0.0E+00

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

  call sptsl ( n, d, e, b )
!
!  Print the result.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 entries of the solution:'
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
!! TEST23 tests SQRDC and SQRSL.
!
!  Discussion:
!
!    SQRDC and SQRSL compute the QR factorization, and use it
!    to solve linear systems.
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

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: p = 3
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,p)
  real ( kind = 4 ) b(lda,p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(p)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real ( kind = 4 ) q(n,n)
  real ( kind = 4 ) qraux(p)
  real ( kind = 4 ) qty(n)
  real ( kind = 4 ) qy(n)
  real ( kind = 4 ) r(n,p)
  real ( kind = 4 ) rsd(n)
  real ( kind = 4 ) work(p)
  real ( kind = 4 ) xb(n)
  real ( kind = 4 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  SQRDC computes the QR decomposition of a '
  write ( *, '(a)' ) '  matrix, but does not return Q and R explicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Show how Q and R can be recovered using SQRSL.'
!
!  Set the matrix A.
!
  a(1,1) = 1.0E+00
  a(2,1) = 1.0E+00
  a(3,1) = 0.0E+00

  a(1,2) = 1.0E+00
  a(2,2) = 0.0E+00
  a(3,2) = 1.0E+00

  a(1,3) = 0.0E+00
  a(2,3) = 1.0E+00
  a(3,3) = 1.0E+00

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

  call sqrdc ( a, lda, n, p, qraux, ipvt, work, job )
!
!  Print out what SQRDC has stored in A...
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
        r(i,j) = 0.0E+00
      else
        r(i,j) = a(i,j)
      end if

    end do

    write ( *, '(2x,5g14.6)' ) r(i,1:p)

  end do
!
!  Call SQRSL to extract the information about the Q matrix.
!  We do this, essentially, by asking SQRSL to tell us the
!  value of Q*Y, where Y is a column of the identity matrix.
!
  job = 10000

  do i = 1, n
!
!  Set the vector Y.
!
    y(1:n) = 0.0E+00

    y(i) = 1.0E+00
!
!  Ask SQRSL to tell us what Q*Y is.
!
    call sqrsl ( a, lda, n, p, qraux, y, qy, qty, b, rsd, xb, job, info )

    if ( info /= 0 ) then
      write ( *, '(a,i8)' ) '  Error!  SQRSL returns INFO = ', info
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
!! TEST24 tests SSICO.
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

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  real ( kind = 4 ) rcond
  real ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and symmetric indefinite storage (SI),'
  write ( *, '(a)' ) '  SSICO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to the matrix A.
!
  a(1:n,1:n) = 0.0E+00

  do i = 1, n
    a(i,i) = 2.0E+00
    if ( i < n ) then
      a(i,i+1) = -1.0E+00
    end if
  end do
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call ssico ( a, lda, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition = ', rcond

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests SSIFA and SSISL.
!
!  Discussion:
!
!    SSIFA and SSISL are for symmetric indefinite matrices.
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

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and symmetric indefinite storage (SI),'
  write ( *, '(a)' ) '  SSIFA factors the matrix,'
  write ( *, '(a)' ) '  SSISL solves a factored linear system,'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to the matrix A and the right hand side B.
!
  b(1:n-1) = 0.0E+00
  b(n)= real ( n + 1 )

  a(1:n,1:n) = 0.0E+00

  do i = 1, n
    a(i,i) = 2.0E+00
    if ( i < n ) then
      a(i,i+1) = -1.0E+00
    end if
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call ssifa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  SSIFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call ssisl ( a, lda, n, ipvt, b )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 entries of the solution:'
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
!! TEST26 tests SSPCO.
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

  integer ( kind = 4 ), parameter :: n = 100

  real ( kind = 4 ) a((n*(n+1))/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) rcond
  real ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and symmetric indefinite packed storage (SP),'
  write ( *, '(a)' ) '  SSPCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to the matrix A.
!
  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      if ( i == j ) then
        a(k) = 2.0E+00
      else if ( j == i+1 ) then
        a(k) = -1.0E+00
      else
        a(k) = 0.0E+00
      end if
    end do
  end do
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call sspco ( a, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition = ', rcond

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests SSPFA and SSPSL.
!
!  Discussion:
!
!    SSPFA and SSPSL are for packed symmetric indefinite matrices.
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

  integer ( kind = 4 ), parameter :: n = 100

  real ( kind = 4 ) a((n*(n+1))/2)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and symmetric indefinite packed storage (SP),'
  write ( *, '(a)' ) '  SSPFA factors the matrix,'
  write ( *, '(a)' ) '  SSPSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Assign values to the matrix A and the right hand side B.
!
  b(1:n-1) = 0.0E+00
  b(n)= real ( n + 1 )

  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      if ( i == j ) then
        a(k) = 2.0E+00
      else if ( j == i+1 ) then
        a(k) = -1.0E+00
      else
        a(k) = 0.0E+00
      end if
    end do
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call sspfa ( a, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  SSPFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call sspsl ( a, n, ipvt, b )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first and last 5 entries of the solution:'
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
!! TEST28 tests SSVDC.
!
!  Discussion:
!
!    SSVDC computes the singular value decomposition:
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

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) e(max(m+1,n))
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) job
  real ( kind = 4 ) s(max(m+1,n))
  integer ( kind = 4 ) seed
  real ( kind = 4 ) sigma(m,n)
  real ( kind = 4 ) u(m,m)
  real ( kind = 4 ) v(n,n)
  real ( kind = 4 ) work(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and an MxN matrix A in general storage,'
  write ( *, '(a)' ) '  SSVDC computes the singular value decomposition:'
  write ( *, '(a)' ) '    A = U * S * V'''
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set A.
!
  seed = 123456789

  call r4mat_uniform_01 ( m, n, seed, a )

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

  call ssvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Warning:'
    write ( *, '(a,i8)' ) '  SSVDC returned nonzero INFO = ', info
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

  sigma(1:m,1:n) = 0.0E+00
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
!! TEST29 tests STRCO.
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

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real ( kind = 4 ) rcond
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and triangular storage (TR),'
  write ( *, '(a)' ) '  STRCO computes the LU factors and'
  write ( *, '(a)' ) '  computes its reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Lower triangular matrix A.
!
  call r4mat_uniform_01 ( n, n, seed, a )

  do i = 1, n
    do j = i+1, n
      a(i,j) = 0.0E+00
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

  call strco ( a, lda, n, rcond, z, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The reciprocal condition number = ', rcond
!
!  Upper triangular matrix A.
!
  call r4mat_uniform_01 ( n, n, seed, a )

  do i = 1, n
    do j = 1, i - 1
      a(i,j) = 0.0E+00
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

  call strco ( a, lda, n, rcond, z, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The reciprocal condition number = ', rcond

  return
end
subroutine test30 ( )

!*****************************************************************************80
!
!! TEST30 tests STRDI.
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

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and triangular storage (TR),'
  write ( *, '(a)' ) '  STRDI computes the determinant or inverse.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Lower triangular matrix A.
!
  call r4mat_uniform_01 ( n, n, seed, a )

  do i = 1, n
    do j = i+1, n
      a(i,j) = 0.0E+00
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lower triangular matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') a(i,1:n)
  end do

  job = 110

  call strdi ( a, lda, n, det, job, info )

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
  call r4mat_uniform_01 ( n, n, seed, a )

  do i = 1, n
    do j = 1, i - 1
      a(i,j) = 0.0E+00
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Upper triangular matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5g14.6)') a(i,1:n)
  end do

  job = 111

  call strdi ( a, lda, n, det, job, info )

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
!! TEST31 tests STRSL.
!
!  Discussion:
!
!    STRSL solves triangular linear systems.
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

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  For single precision real arithmetic (S)'
  write ( *, '(a)' ) '  and triangular storage (TR),'
  write ( *, '(a)' ) '  STRSL solves a linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Lower triangular matrix A.
!
  call r4mat_uniform_01 ( n, n, seed, a )

  do i = 1, n
    do j = i+1, n
      a(i,j) = 0.0E+00
    end do
  end do

  do i = 1, n
    x(i) = real ( i )
  end do

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For a lower triangular matrix A,'
  write ( *, '(a)' ) '  solve A * x = b'

  job = 00

  call strsl ( a, lda, n, b, job, info )

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

  call strsl ( a, lda, n, b, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The solution (should be 1,2,3,4,5):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
  end do
!
!  Upper triangular matrix A.
!
  call r4mat_uniform_01 ( n, n, seed, a )

  do i = 1, n
    do j = 1, i - 1
      a(i,j) = 0.0E+00
    end do
  end do

  do i = 1, n
    x(i) = real ( i )
  end do

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For an upper triangular matrix A,'
  write ( *, '(a)' ) '  solve A * x = b'

  job = 01

  call strsl ( a, lda, n, b, job, info )

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

  call strsl ( a, lda, n, b, job, info )

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
subroutine r4mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R4MAT_UNIFORM_01 returns a unit pseudorandom R4MAT.
!
!  Discussion:
!
!    An R4MAT is an array of real ( kind = 4 ) values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2005
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
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input/output, integer SEED, the "seed" value, which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r(m,n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4MAT_UNIFORM_01 - Fatal error!'
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

      r(i,j) = real ( seed, kind = 4 ) * 4.656612875E-10

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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
