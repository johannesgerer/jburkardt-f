program main

!*****************************************************************************80
!
!! MAIN is the main program for LINPACK_Z_PRB.
!
!  Discussion:
!
!    LINPACK_Z_PRB calls the double precision complex LINPACK test routines.
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
  write ( *, '(a)' ) 'LINPACK_Z_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LINPACK_Z library.'

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
  call test32 ( )
  call test33 ( )
  call test34 ( )
  call test345 ( )
  call test35 ( )
  call test36 ( )
  call test37 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINPACK_Z_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests ZCHDC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  complex ( kind = 8 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) &
    '  For a double complex Hermitian positive definite matrix,'
  write ( *, '(a)' ) '  ZCHDC computes the Cholesky decomposition.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of equations is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = cmplx ( 2.5281D+00,  0.0000D+00, kind = 8 )
  a(2,1) = cmplx ( 2.1341D+00,  0.2147D+00, kind = 8 )
  a(3,1) = cmplx ( 2.4187D+00, -0.2932D+00, kind = 8 )

  a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00, kind = 8 )
  a(2,2) = cmplx ( 3.0371D+00,  0.0000D+00, kind = 8 )
  a(3,2) = cmplx ( 2.0905D+00, -1.1505D+00, kind = 8 )

  a(1,3) = cmplx ( 2.4187D+00,  0.2932D+00, kind = 8 )
  a(2,3) = cmplx ( 2.0905D+00,  1.1505D+00, kind = 8 )
  a(3,3) = cmplx ( 2.7638D+00,  0.0000D+00, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,6f10.4)' ) a(i,1:n)
  end do
!
!  Decompose the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Decompose the matrix.'

  job = 0
  ipvt(1:n) = 0

  call zchdc ( a, lda, n, work, ipvt, job, info )

  if ( info /= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZCHDC returned INFO = ', info
    write ( *, '(a)' ) '  The matrix is not Hermitian positive definite.'
    return
  end if
!
!  Zero out the lower diagonal.
!
  do i = 2, n
    do j = 1, i-1
      a(i,j) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    end do
  end do
!
!  Print the factorization.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Cholesky factor U:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,6f10.4)' ) a(i,1:n)
  end do
!
!  Compute the Cholesky product.
!
  a(1:n,1:n) = matmul ( conjg ( transpose ( a(1:n,1:n) ) ), a(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product U^H * U: '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,6f10.4)' ) a(i,1:n)
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests ZCHEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n
  integer ( kind = 4 ), parameter :: ldz = n
  integer ( kind = 4 ), parameter :: nz = 1

  complex ( kind = 8 ) a(lda,n)
  real    ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  complex ( kind = 8 ) s(n)
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) work(n)
  complex ( kind = 8 ) z(ldz,nz)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) &
    '  For a double complex Hermitian positive definite matrix,'
  write ( *, '(a)' ) '  ZCHEX can shift rows and columns in a Cholesky factorization.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of equations is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = cmplx ( 2.5281D+00,  0.0000D+00, kind = 8 )
  a(2,1) = cmplx ( 2.1341D+00,  0.2147D+00, kind = 8 )
  a(3,1) = cmplx ( 2.4187D+00, -0.2932D+00, kind = 8 )

  a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00, kind = 8 )
  a(2,2) = cmplx ( 3.0371D+00,  0.0000D+00, kind = 8 )
  a(3,2) = cmplx ( 2.0905D+00, -1.1505D+00, kind = 8 )

  a(1,3) = cmplx ( 2.4187D+00,  0.2932D+00, kind = 8 )
  a(2,3) = cmplx ( 2.0905D+00,  1.1505D+00, kind = 8 )
  a(3,3) = cmplx ( 2.7638D+00,  0.0000D+00, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,6f10.4)' ) a(i,1:n)
  end do

  do i = 1, n
    z(i,1) = cmplx ( i, 0.0D+00, kind = 8 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The vector Z:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,2g14.6)' ) z(i,1)
  end do
!
!  Decompose the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Decompose the matrix.'

  job = 0
  ipvt(1:n) = 0

  call zchdc ( a, lda, n, work, ipvt, job, info )

  if ( info /= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZCHDC returned INFO = ', info
    write ( *, '(a)' ) '  This means the matrix is not positive definite.'
    return
  end if
!
!  Zero out the lower diagonal.
!
  do i = 2, n
    do j = 1, i-1
      a(i,j) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    end do
  end do
!
!  Print the factorization.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Cholesky factor U:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,6f10.4)' ) a(i,1:n)
  end do
!
!  Right circular shift columns L through K.
!
  k = 1
  l = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  Right circular shift rows and columns K  = ', k, &
    ' through L = ', l
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Logical matrix is now:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  33 31 32'
  write ( *, '(a)' ) '  13 11 12'
  write ( *, '(a)' ) '  23 21 22'

  job = 1
  call zchex ( a, lda, n, k, l, z, ldz, nz, c, s, job )
!
!  Left circular shift columns K+1 through L.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  Left circular shift rows and columns K+1 = ', k+1, &
    ' through L = ', l
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Logical matrix is now:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  33 32 31'
  write ( *, '(a)' ) '  23 22 21'
  write ( *, '(a)' ) '  13 12 11'

  job = 2
  call zchex ( a, lda, n, k+1, l, z, ldz, nz, c, s, job )
!
!  Print the factorization.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The shifted Cholesky factor UU:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,6f10.4)' ) a(i,1:n)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The shifted vector ZZ:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,2g14.6)' ) z(i,1)
  end do
!
!  Compute the Cholesky product.
!
  a(1:n,1:n) = matmul ( conjg ( transpose ( a(1:n,1:n) ) ), a(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The shifted product AA = UU'' * UU: '
  write ( *, '(a)' ) '  The rows and columns of the original matrix A reappear,'
  write ( *, '(a)' ) '  but in reverse order.'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,6f10.4)' ) a(i,1:n)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests ZCHUD and ZTRSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: p = 20
  integer ( kind = 4 ), parameter :: ldr = p
  integer ( kind = 4 ), parameter :: ldz = p
  integer ( kind = 4 ), parameter :: nz = 1

  complex ( kind = 8 ) b(p)
  complex ( kind = 8 ) beta(p)
  real    ( kind = 8 ) c(p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  complex ( kind = 8 ) r(ldr,p)
  real    ( kind = 8 ) rho(nz)
  complex ( kind = 8 ) row(p)
  complex ( kind = 8 ) s(p)
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) x(p)
  complex ( kind = 8 ) y(nz)
  complex ( kind = 8 ) z(ldz,nz)
  complex ( kind = 8 ) zdotu

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For a double complex Hermitian matrix'
  write ( *, '(a)' ) '  ZCHUD updates a Cholesky decomposition.'
  write ( *, '(a)' ) '  ZTRSL solves a triangular linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we use ZCHUD to solve a'
  write ( *, '(a)' ) '  least squares problem R * b = z.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of equations is P = ', p
!
!  Initialize.
!
  r(1:p,1:p) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
  z(1:p,1:nz) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  do i = 1, p
    x(i) = cmplx ( i, mod ( i, 2 ), kind = 8 )
  end do
!
!  Use ZCHUD to form R, Z and RHO by adding X and Y a row at a time.
!  X is a row of the least squares matrix and Y the right hand side.
!
  seed = 123456789

  do i = 1, p
    call c8vec_uniform_01 ( p, seed, row )
    y(1) = zdotu ( p, row, 1, x, 1 )
    rho(1) = 0.0D+00
    call zchud ( r, ldr, p, row, z, ldz, nz, y, rho, c, s )
  end do
!
!  Generate the least squares solution, b = inverse ( R ) * Z.
!
  do j = 1, nz

    b(1:p) = z(1:p,j)
    job = 01

    call ztrsl ( r, ldr, p, b, job, info )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Solution vector # ', j
    write ( *, '(a)' ) '  (Should be (1,1) (2,0), (3,1) (4,0) ...)'
    write ( *, '(a)' ) ' '

    do i = 1, p
      if ( i <= 5 .or. p-5 < i ) then
        write ( *, '(2x,i8,2x,2g14.6)' ) i, b(i)
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
!! TEST04 tests ZGBCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ), parameter :: lda = 2*ml+mu+1

  complex ( kind = 8 ) a(lda,n)
  complex ( kind = 8 ) a_save(n,n)
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real    ( kind = 8 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For a double complex general band storage matrix:'
  write ( *, '(a)' ) '  ZGBCO factors the matrix and estimates the'
  write ( *, '(a)' ) '  reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
  write ( *, '(a,i8)' ) '  The lower band is ML =  ', ml
  write ( *, '(a,i8)' ) '  The upper band is MU =  ', mu
!
!  Set the values of the matrix A.
!
  a_save(1:n,1:n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  m = ml + mu + 1

  seed = 123456789

  do j = 1, n
    i1 = max ( 1, j - mu )
    i2 = min ( n, j + ml )
    do i = i1, i2
      k = i - j + m
      a(k,j) = c8_uniform_01 ( seed )
      a_save(i,j) = a(k,j)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a_save(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call zgbco ( a, lda, n, ml, mu, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition RCOND = ', rcond

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests ZGBFA and ZGBSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ), parameter :: lda = 2*ml+mu+1

  complex ( kind = 8 ) a(lda,n)
  complex ( kind = 8 ) a_save(n,n)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For a double complex general band storage matrix:'
  write ( *, '(a)' ) '  ZGBFA factors the matrix;'
  write ( *, '(a)' ) '  ZGBSL solves a factored linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
  write ( *, '(a,i8)' ) '  The lower band is ML =  ', ml
  write ( *, '(a,i8)' ) '  The upper band is MU =  ', mu
!
!  Set the values of the matrix A.
!
  a_save(1:n,1:n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  m = ml + mu + 1

  seed = 123456789

  do j = 1, n
    i1 = max ( 1, j - mu )
    i2 = min ( n, j + ml )
    do i = i1, i2
      k = i - j + m
      a(k,j) = c8_uniform_01 ( seed )
      a_save(i,j) = a(k,j)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a_save(i,1:n)
  end do
!
!  Set the values of the right hand side vector B.
!
  call c8vec_uniform_01 ( n, seed, x )

  b(1:n) = matmul ( a_save(1:n,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The right hand side B is '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2f8.4)' ) b(i)
  end do
!
!  Factor the matrix A.
!
  call zgbfa ( a, lda, n, ml, mu, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  ZGBFA returned INFO = ', info
    return
  end if
!
!  Solve the system.
!
  job = 0
  call zgbsl ( a, lda, n, ml, mu, ipvt, b, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed                     Exact'
  write ( *, '(a)' ) '  Solution                     Solution'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(4g14.6)' ) b(i), x(i)
  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests ZGBFA and ZGBDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ), parameter :: lda = 2*ml+mu+1

  complex ( kind = 8 ) a(lda,n)
  complex ( kind = 8 ) a_save(n,n)
  complex ( kind = 8 ) c8_uniform_01
  complex ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  For a double complex general band storage matrix:'
  write ( *, '(a)' ) '  ZGBFA factors the matrix.'
  write ( *, '(a)' ) '  ZGBDI computes the determinant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
  write ( *, '(a,i8)' ) '  The lower band is ML =  ', ml
  write ( *, '(a,i8)' ) '  The upper band is MU =  ', mu
!
!  Set the values of the matrix A.
!
  a_save(1:n,1:n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  m = ml + mu + 1

  seed = 123456789

  do j = 1, n
    i1 = max ( 1, j - mu )
    i2 = min ( n, j + ml )
    do i = i1, i2
      k = i - j + m
      a(k,j) = c8_uniform_01 ( seed )
      a_save(i,j) = a(k,j)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a_save(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call zgbfa ( a, lda, n, ml, mu, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  ZGBFA returned INFO = ', info
    return
  end if
!
!  Get the determinant.
!
  call zgbdi ( a, lda, n, ml, mu, ipvt, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', real ( det(2),  kind = 8 )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests ZGECO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  real    ( kind = 8 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) z(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  For a double complex general storage matrix:'
  write ( *, '(a)' ) '  ZGECO factors the matrix and estimates the'
  write ( *, '(a)' ) '  reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  call c8mat_uniform_01 ( n, n, seed, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call zgeco ( a, lda, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition RCOND = ', rcond

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests ZGEFA and ZGESL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) ai
  real    ( kind = 8 ) ar
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) x(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  For a double complex general storage matrix:'
  write ( *, '(a)' ) '  ZGEFA factors the matrix.'
  write ( *, '(a)' ) '  ZGESL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  call c8mat_uniform_01 ( n, n, seed, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a(i,1:n)
  end do
!
!  Set the values of the right hand side vector B.
!
  call c8vec_uniform_01 ( n, seed, x )

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The right hand side B is '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2f8.4)' ) b(i)
  end do
!
!  Factor the matrix A.
!
  call zgefa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZGEFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the system.
!
  job = 0
  call zgesl ( a, lda, n, ipvt, b, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed                     Exact'
  write ( *, '(a)' ) '  Solution                     Solution'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(4g14.6)' ) b(i), x(i)
  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests ZGEFA and ZGEDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a(n,n)
  complex ( kind = 8 ) a_save(n,n)
  real    ( kind = 8 ) ai
  real    ( kind = 8 ) ar
  complex ( kind = 8 ) c(n,n)
  complex ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) work(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  For a double complex general storage matrix:'
  write ( *, '(a)' ) '  ZGEFA factors the matrix.'
  write ( *, '(a)' ) '  ZGEDI computes the determinant or inverse.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  call c8mat_uniform_01 ( n, n, seed, a )

  a_save(1:n,1:n) = a(1:n,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call zgefa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZGEFA returned an error flag INFO = ', info
    return
  end if
!
!  Get the determinant.
!
  job = 10
  call zgedi ( a, lda, n, ipvt, det, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', real ( det(2), kind = 8 )
!
!  Get the inverse.
!
  job = 01
  call zgedi ( a, lda, n, ipvt, det, work, job )

  c(1:n,1:n) = matmul ( a(1:n,1:n), a_save(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product inv(A) * A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) c(i,1:n)
  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests ZGTSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) c(n)
  complex ( kind = 8 ) d(n)
  complex ( kind = 8 ) e(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  For a double complex tridiagonal matrix:'
  write ( *, '(a)' ) '  ZGTSL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  seed = 123456789

  c(1) = cmplx ( 0.0E+00, 0.0E+00, kind = 8 )
  call c8vec_uniform_01 ( n-1, seed, c(2) )

  call c8vec_uniform_01 ( n-1, seed, e(1) )
  e(n) = cmplx ( 0.0E+00, 0.0E+00, kind = 8 )

  d(1:n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  d(1:n-1) = d(1:n-1) - 2.0D+00 * e(1:n-1)
  d(2:n)   = d(2:n) - 2.0D+00 * c(2:n)
!
!  Set the desired solution
!
  do i = 1, n
    x(i) = cmplx ( i, 10 * i, kind = 8 )
  end do
!
!  Compute the corresponding right hand side.
!
  b(1) = d(1) * x(1) + e(1) * x(2)
  do i = 2, n-1
    b(i) = c(i) * x(i-1) + d(i) * x(i) + e(i) * x(i+1)
  end do
  b(n) = c(n) * x(n-1) + d(n) * x(n)
!
!  Solve the tridiagonal system.
!
  call zgtsl ( n, c, d, e, b, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed                     Exact'
  write ( *, '(a)' ) '  Solution                     Solution'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(4g14.6)' ) b(i), x(i)
  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests ZHICO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) ai
  real    ( kind = 8 ) ar
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  real    ( kind = 8 ) r8_uniform_01
  real    ( kind = 8 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) z(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  For a double complex Hermitian matrix:'
  write ( *, '(a)' ) '  ZHICO factors the matrix and estimates'
  write ( *, '(a)' ) '  the reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  do i = 1, n
    a(i,i) = cmplx ( r8_uniform_01 ( seed ), 0.0D+00 )
    do j = i+1, n
      a(i,j) = c8_uniform_01 ( seed )
      a(j,i) = conjg ( a(i,j) )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call zhico ( a, lda, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition RCOND = ', rcond

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests ZHIFA and ZHISL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a(n,n)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) x(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  For a double complex Hermitian matrix:'
  write ( *, '(a)' ) '  ZHIFA factors the matrix.'
  write ( *, '(a)' ) '  ZHISL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  do i = 1, n
    a(i,i) = cmplx ( r8_uniform_01 ( seed ), 0.0D+00 )
    do j = i+1, n
      a(i,j) = c8_uniform_01 ( seed )
      a(j,i) = conjg ( a(i,j) )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a(i,1:n)
  end do
!
!  Set the values of the right hand side vector B.
!
  call c8vec_uniform_01 ( n, seed, x )

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The right hand side B is '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2f8.4)' ) b(i)
  end do
!
!  Factor the matrix A.
!
  call zhifa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZHIFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the system.
!
  call zhisl ( a, lda, n, ipvt, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed                     Exact'
  write ( *, '(a)' ) '  Solution                     Solution'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(4g14.6)' ) b(i), x(i)
  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests ZHIFA and ZHIDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a(n,n)
  complex ( kind = 8 ) a_save(n,n)
  complex ( kind = 8 ) c(n,n)
  complex ( kind = 8 ) c8_uniform_01
  real    ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inert(3)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) work(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  For a double complex hermitian matrix:'
  write ( *, '(a)' ) '  ZHIFA factors the matrix.'
  write ( *, '(a)' ) '  ZHIDI computes the determinant, inverse,'
  write ( *, '(a)' ) '  or inertia.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  do i = 1, n
    a(i,i) = cmplx ( r8_uniform_01 ( seed ), 0.0D+00 )
    do j = i+1, n
      a(i,j) = c8_uniform_01 ( seed )
      a(j,i) = conjg ( a(i,j) )
    end do
  end do

  a_save(1:n,1:n) = a(1:n,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call zhifa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZHIFA returned an error flag INFO = ', info
    return
  end if
!
!  Get the determinant.
!
  job = 010
  call zhidi ( a, lda, n, ipvt, det, inert, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', det(2)
!
!  Get the inertia.
!
  job = 100
  call zhidi ( a, lda, n, ipvt, det, inert, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The inertia:'
  write ( *, '(a)' ) ' '

  do i = 1, 3
    write ( *, '(2x,i8)' ) inert(i)
  end do
!
!  Get the inverse.
!
  job = 001
  call zhidi ( a, lda, n, ipvt, det, inert, work, job )
!
!  Only the upper triangle is set, so the user must set up the
!  lower triangle:
!
  do i = 1, n
    do j = 1, i-1
      a(i,j) = conjg ( a(j,i) )
    end do
  end do

  c(1:n,1:n) = matmul ( a(1:n,1:n), a_save(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product inv(A) * A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) c(i,1:n)
  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests ZHPCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a((n*(n+1))/2)
  complex ( kind = 8 ) a_save(n,n)
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r8_uniform_01
  real    ( kind = 8 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  For a double complex Hermitian matrix'
  write ( *, '(a)' ) '  using packed storage,'
  write ( *, '(a)' ) '  ZHPCO factors the matrix and estimates'
  write ( *, '(a)' ) '  the reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  k = 0
  seed = 123456789

  do j = 1, n

    do i = 1, j-1
      k = k + 1
      a(k) = c8_uniform_01 ( seed )
      a_save(i,j) = a(k)
      a_save(j,i) = conjg ( a(k) )
    end do

    k = k + 1
    a(k) = cmplx ( r8_uniform_01 ( seed ), 0.0D+00, kind = 8 )
    a_save(j,j) = a(k)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a_save(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call zhpco ( a, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition RCOND = ', rcond

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests ZHPFA and ZHPSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a((n*(n+1))/2)
  complex ( kind = 8 ) a_save(n,n)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  For a double complex Hermitian matrix,'
  write ( *, '(a)' ) '  using packed storage,'
  write ( *, '(a)' ) '  ZHPFA factors the matrix.'
  write ( *, '(a)' ) '  ZHPSL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  k = 0
  seed = 123456789

  do j = 1, n

    do i = 1, j-1
      k = k + 1
      a(k) = c8_uniform_01 ( seed )
      a_save(i,j) = a(k)
      a_save(j,i) = conjg ( a(k) )
    end do

    k = k + 1
    a(k) = cmplx ( r8_uniform_01 ( seed ), 0.0D+00, kind = 8 )
    a_save(j,j) = a(k)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a_save(i,1:n)
  end do
!
!  Set the values of the right hand side vector B.
!
  call c8vec_uniform_01 ( n, seed, x )

  b(1:n) = matmul ( a_save(1:n,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The right hand side B is '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2f8.4)' ) b(i)
  end do
!
!  Factor the matrix A.
!
  call zhpfa ( a, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZHPFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the system.
!
  call zhpsl ( a, n, ipvt, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed                     Exact'
  write ( *, '(a)' ) '  Solution                     Solution'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(4g14.6)' ) b(i), x(i)
  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests ZHPFA and ZHPDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a((n*(n+1))/2)
  complex ( kind = 8 ) a_save(n,n)
  complex ( kind = 8 ) b(n,n)
  complex ( kind = 8 ) c(n,n)
  complex ( kind = 8 ) c8_uniform_01
  real    ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inert(3)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  For a double complex hermitian matrix,'
  write ( *, '(a)' ) '  using packed storage,'
  write ( *, '(a)' ) '  ZHPFA factors the matrix.'
  write ( *, '(a)' ) '  ZHPDI computes the determinant, inverse,'
  write ( *, '(a)' ) '  or inertia.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  k = 0
  seed = 123456789

  do j = 1, n

    do i = 1, j-1
      k = k + 1
      a(k) = c8_uniform_01 ( seed )
      a_save(i,j) = a(k)
      a_save(j,i) = conjg ( a(k) )
    end do

    k = k + 1
    a(k) = cmplx ( r8_uniform_01 ( seed ), 0.0D+00, kind = 8 )
    a_save(j,j) = a(k)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a_save(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call zhpfa ( a, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZHPFA returned an error flag INFO = ', info
    return
  end if
!
!  Get the determinant.
!
  job = 010
  call zhpdi ( a, n, ipvt, det, inert, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', det(2)
!
!  Get the inertia.
!
  job = 100
  call zhpdi ( a, n, ipvt, det, inert, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The inertia:'
  write ( *, '(a)' ) ' '

  do i = 1, 3
    write ( *, '(2x,i8)' ) inert(i)
  end do
!
!  Get the inverse.
!
  job = 001
  call zhpdi ( a, n, ipvt, det, inert, work, job )
!
!  Only the upper triangle is set, so the user must set up the
!  lower triangle:
!
  k = 0
  do j = 1, n
    do i = 1, j-1
      k = k + 1
      b(i,j) = a(k)
      b(j,i) = conjg ( a(k) )
    end do
    k = k + 1
    b(j,j) = a(k)
  end do

  c(1:n,1:n) = matmul ( b(1:n,1:n), a_save(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product inv(A) * A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) c(i,1:n)
  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests ZPBCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 1
  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ), parameter :: lda = m + 1

  complex ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real    ( kind = 8 ) rcond
  complex ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  For a double complex '
  write ( *, '(a)' ) '  positive definite hermitian band matrix,'
  write ( *, '(a)' ) '  ZPBCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the value of the superdiagonal and diagonal.
!
  a(1,1) = cmplx ( 0.0000D+00,  0.0000D+00, kind = 8 )
  a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00, kind = 8 )
  a(1,3) = cmplx ( 2.0905D+00,  1.1505D+00, kind = 8 )

  a(2,1) = cmplx ( 4.5281D+00,  0.0000D+00, kind = 8 )
  a(2,2) = cmplx ( 5.0371D+00,  0.0000D+00, kind = 8 )
  a(2,3) = cmplx ( 4.7638D+00,  0.0000D+00, kind = 8 )
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call zpbco ( a, lda, n, m, rcond, z, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZPBCO returned INFO = ', info
    write ( *, '(a)' ) '  The factorization was not completed.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Reciprocal condition  = ', rcond

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests ZPBDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ), parameter :: lda = m+1

  complex ( kind = 8 ) a(lda,n)
  real    ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  For a double complex '
  write ( *, '(a)' ) '  positive definite hermitian band matrix,'
  write ( *, '(a)' ) '  ZPBDI computes the determinant as'
  write ( *, '(a)' ) '    det = MANTISSA * 10**EXPONENT'
  write ( *, '(a)' ) ' '
!
!  Set the value of the superdiagonal and diagonal.
!
  a(1,1) = cmplx ( 0.0000D+00,  0.0000D+00, kind = 8 )
  a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00, kind = 8 )
  a(1,3) = cmplx ( 2.0905D+00,  1.1505D+00, kind = 8 )

  a(2,1) = cmplx ( 4.5281D+00,  0.0000D+00, kind = 8 )
  a(2,2) = cmplx ( 5.0371D+00,  0.0000D+00, kind = 8 )
  a(2,3) = cmplx ( 4.7638D+00,  0.0000D+00, kind = 8 )

  call zpbfa ( a, lda, n, m, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  ZPBFA returns INFO = ', info
    return
  end if

  call zpbdi ( a, lda, n, m, det )

  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', det(2)

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests ZPBFA and ZPBSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ), parameter :: lda = m+1

  complex ( kind = 8 ) a(lda,n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  For a double complex'
  write ( *, '(a)' ) '  positive definite hermitian band matrix,'
  write ( *, '(a)' ) '  ZPBFA computes the LU factors.'
  write ( *, '(a)' ) '  ZPBSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the value of the superdiagonal and diagonal.
!
  a(1,1) = cmplx ( 0.0000D+00,  0.0000D+00, kind = 8 )
  a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00, kind = 8 )
  a(1,3) = cmplx ( 2.0905D+00,  1.1505D+00, kind = 8 )

  a(2,1) = cmplx ( 4.5281D+00,  0.0000D+00, kind = 8 )
  a(2,2) = cmplx ( 5.0371D+00,  0.0000D+00, kind = 8 )
  a(2,3) = cmplx ( 4.7638D+00,  0.0000D+00, kind = 8 )
!
!  Set the right hand side.
!
  b(1) = cmplx (  8.7963D+00, -0.4294D+00, kind = 8 )
  b(2) = cmplx ( 18.4798D+00,  3.6662D+00, kind = 8 )
  b(3) = cmplx ( 18.4724D+00, -2.3010D+00, kind = 8 )
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call zpbfa ( a, lda, n, m, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  ZPBFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call zpbsl ( a, lda, n, m, b )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The solution:'
  write ( *, '(a)' ) '  (Should be roughly (1,2,3)):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,2g14.6)' ) i, b(i)
  end do

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests ZPOCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real    ( kind = 8 ) rcond
  complex ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) &
    '  For a double complex Hermitian positive definite matrix,'
  write ( *, '(a)' ) '  ZPOCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = cmplx ( 2.5281D+00,  0.0000D+00, kind = 8 )
  a(2,1) = cmplx ( 2.1341D+00,  0.2147D+00, kind = 8 )
  a(3,1) = cmplx ( 2.4187D+00, -0.2932D+00, kind = 8 )

  a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00, kind = 8 )
  a(2,2) = cmplx ( 3.0371D+00,  0.0000D+00, kind = 8 )
  a(3,2) = cmplx ( 2.0905D+00, -1.1505D+00, kind = 8 )

  a(1,3) = cmplx ( 2.4187D+00,  0.2932D+00, kind = 8 )
  a(2,3) = cmplx ( 2.0905D+00,  1.1505D+00, kind = 8 )
  a(3,3) = cmplx ( 2.7638D+00,  0.0000D+00, kind = 8 )
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call zpoco ( a, lda, n, rcond, z, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Reciprocal condition  = ', rcond

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests ZPOFA and ZPODI.
!
!  Discussion:
!
!    ZPOFA factors a positive definite symmetric matrix,
!    and ZPODI can compute the determinant or the inverse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 8 ) a(lda,n)
  real    ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) &
    '  For a double complex Hermitian positive definite matrix,'
  write ( *, '(a)' ) '  ZPOFA computes the LU factors,'
  write ( *, '(a)' ) '  ZPODI computes the inverse or determinant.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = cmplx ( 2.5281D+00,  0.0000D+00, kind = 8 )
  a(2,1) = cmplx ( 2.1341D+00,  0.2147D+00, kind = 8 )
  a(3,1) = cmplx ( 2.4187D+00, -0.2932D+00, kind = 8 )

  a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00, kind = 8 )
  a(2,2) = cmplx ( 3.0371D+00,  0.0000D+00, kind = 8 )
  a(3,2) = cmplx ( 2.0905D+00, -1.1505D+00, kind = 8 )

  a(1,3) = cmplx ( 2.4187D+00,  0.2932D+00, kind = 8 )
  a(2,3) = cmplx ( 2.0905D+00,  1.1505D+00, kind = 8 )
  a(3,3) = cmplx ( 2.7638D+00,  0.0000D+00, kind = 8 )
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call zpofa ( a, lda, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, ZPOFA returns INFO = ', info
    return
  end if
!
!  Get the determinant and inverse.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get the determinant and inverse.'

  job = 11
  call zpodi ( a, lda, n, det, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant  = ', det(1), ' * 10 ** ', det(2)
!
!  ZPODI produces only the 'upper half triangle' of the inverse,
!  which is actually symmetric.  Thus, the lower half could be
!  produced by copying from the upper half.  However, the first row
!  of A, as returned, is exactly the first row of the inverse.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First row of inverse:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,6f10.4)' ) a(1,1:n)

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests ZPOFA and ZPOSL.
!
!  Discussion:
!
!    ZPOFA factors a Hermitian positive definite matrix,
!    and ZPOSL can solve a factored linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 8 ) a(lda,n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) &
    '  For a double complex Hermitian positive definite matrix,'
  write ( *, '(a)' ) '  ZPOFA computes the LU factors.'
  write ( *, '(a)' ) '  ZPOSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = cmplx ( 2.5281D+00,  0.0000D+00, kind = 8 )
  a(2,1) = cmplx ( 2.1341D+00,  0.2147D+00, kind = 8 )
  a(3,1) = cmplx ( 2.4187D+00, -0.2932D+00, kind = 8 )

  a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00, kind = 8 )
  a(2,2) = cmplx ( 3.0371D+00,  0.0000D+00, kind = 8 )
  a(3,2) = cmplx ( 2.0905D+00, -1.1505D+00, kind = 8 )

  a(1,3) = cmplx ( 2.4187D+00,  0.2932D+00, kind = 8 )
  a(2,3) = cmplx ( 2.0905D+00,  1.1505D+00, kind = 8 )
  a(3,3) = cmplx ( 2.7638D+00,  0.0000D+00, kind = 8 )
!
!  Set the right hand side.
!
  do i = 1, n
    x(i) = cmplx ( 2 * i - 1, 2 * i, kind = 8  )
  end do

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call zpofa ( a, lda, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, ZPOFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call zposl ( a, lda, n, b )
!
!  Print the result.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The solution:'
  write ( *, '(a)' ) '  (Should be (1+2i),(3+4i),(5+6i):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,2g14.6)' ) i, b(i)
  end do

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests ZPPCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a((n*(n+1))/2)
  integer ( kind = 4 ) info
  real    ( kind = 8 ) rcond
  complex ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) &
    '  For a double complex Hermitian positive definite packed matrix,'
  write ( *, '(a)' ) '  ZPPCO estimates the reciprocal condition number.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the values of the matrix A.
!
  a(1) = cmplx ( 2.5281D+00,  0.0000D+00, kind = 8 )

  a(2) = cmplx ( 2.1341D+00, -0.2147D+00, kind = 8 )
  a(3) = cmplx ( 3.0371D+00,  0.0000D+00, kind = 8 )

  a(4) = cmplx ( 2.4187D+00,  0.2932D+00, kind = 8 )
  a(5) = cmplx ( 2.0905D+00,  1.1505D+00, kind = 8 )
  a(6) = cmplx ( 2.7638D+00,  0.0000D+00, kind = 8 )
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition number.'

  call zppco ( a, n, rcond, z, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Reciprocal condition number = ', rcond

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests ZPPFA and ZPPDI.
!
!  Discussion:
!
!    ZPPFA factors a Hermitian positive definite packed matrix,
!    and ZPPDI can compute the determinant or the inverse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a((n*(n+1))/2)
  complex ( kind = 8 ) b(n,n)
  real    ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) &
    '  For a double complex Hermitian positive definite packed matrix,'
  write ( *, '(a)' ) '  ZPPFA factors the matrix.'
  write ( *, '(a)' ) '  ZPPDI computes the inverse or determinant.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the values of the matrix A.
!
  a(1) = cmplx ( 2.5281D+00,  0.0000D+00, kind = 8 )

  a(2) = cmplx ( 2.1341D+00, -0.2147D+00, kind = 8 )
  a(3) = cmplx ( 3.0371D+00,  0.0000D+00, kind = 8 )

  a(4) = cmplx ( 2.4187D+00,  0.2932D+00, kind = 8 )
  a(5) = cmplx ( 2.0905D+00,  1.1505D+00, kind = 8 )
  a(6) = cmplx ( 2.7638D+00,  0.0000D+00, kind = 8 )
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call zppfa ( a, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, ZPPFA returns INFO = ', info
    return
  end if
!
!  Invert the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get the determinant and inverse.'

  job = 11
  call zppdi ( a, n, det, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant  = ', det(1), ' * 10 ** ', det(2)
!
!  ZPPDI produces only the 'upper half triangle' of the inverse,
!  which is actually symmetric.  Thus, the lower half could be
!  produced by copying from the upper half.
!
  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      b(i,j) = a(k)
      b(j,i) = conjg ( a(k) )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Inverse:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,6f10.4)' ) b(i,1:n)
  end do

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests ZPPFA and ZPPSL.
!
!  Discussion:
!
!    ZPOFA factors a Hermitian positive definite packed matrix,
!    and ZPOSL can solve a factored linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a((n*(n+1))/2)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) &
    '  For a double complex Hermitian positive definite packed matrix,'
  write ( *, '(a)' ) '  ZPPFA factors the matrix.'
  write ( *, '(a)' ) '  ZPPSL solves a factored linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the values of the matrix A.
!
  a(1) = cmplx ( 2.5281D+00,  0.0000D+00, kind = 8 )

  a(2) = cmplx ( 2.1341D+00, -0.2147D+00, kind = 8 )
  a(3) = cmplx ( 3.0371D+00,  0.0000D+00, kind = 8 )

  a(4) = cmplx ( 2.4187D+00,  0.2932D+00, kind = 8 )
  a(5) = cmplx ( 2.0905D+00,  1.1505D+00, kind = 8 )
  a(6) = cmplx ( 2.7638D+00,  0.0000D+00, kind = 8 )

  b(1) = cmplx ( 20.12350D+00, 28.92670D+00, kind = 8 )
  b(2) = cmplx ( 14.36550D+00, 34.92680D+00, kind = 8 )
  b(3) = cmplx ( 27.69760D+00, 26.03750D+00, kind = 8 )
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call zppfa ( a, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, ZPPFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call zppsl ( a, n, b )
!
!  Print the result.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The solution:'
  write ( *, '(a)' ) '  (Should be (1+2i),(3+4i),(5+6i):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,2g14.6)' ) i, b(i)
  end do

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26 tests ZPTSL.
!
!  Discussion:
!
!    ZPTSL factors and solves a Hermitian positive definite
!    tridiagonal system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) d(n)
  complex ( kind = 8 ) e(n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) &
    '  For a double complex Hermitian positive definite tridiagonal matrix,'
  write ( *, '(a)' ) '  ZPTSL factors and solves a linear system.'
  write ( *, '(a,i8)' ) '  The matrix size is N = ', n
!
!  Set the value of the superdiagonal and diagonal.
!
  e(1) = cmplx ( 2.1341D+00, -0.2147D+00, kind = 8 )
  e(2) = cmplx ( 2.0905D+00,  1.1505D+00, kind = 8 )
  e(3) = cmplx ( 0.0000D+00,  0.0000D+00, kind = 8 )

  d(1) = cmplx ( 4.5281D+00,  0.0000D+00, kind = 8 )
  d(2) = cmplx ( 5.0371D+00,  0.0000D+00, kind = 8 )
  d(3) = cmplx ( 4.7638D+00,  0.0000D+00, kind = 8 )
!
!  Set the right hand side.
!
  b(1) = cmplx (  8.7963D+00, -0.4294D+00, kind = 8 )
  b(2) = cmplx ( 18.4798D+00,  3.6662D+00, kind = 8 )
  b(3) = cmplx ( 18.4724D+00, -2.3010D+00, kind = 8 )
!
!  Factor and solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix and solve the system.'

  call zptsl ( n, d, e, b )
!
!  Print the result.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The solution:'
  write ( *, '(a)' ) '  (Should be roughly (1,2,3)):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,2g14.6)' ) i, b(i)
  end do

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests ZQRDC and ZQRSL.
!
!  Discussion:
!
!    ZQRDC and ZQRSL compute the QR factorization, and use it
!    to solve linear systems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: p = 3
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 8 ) a(lda,p)
  complex ( kind = 8 ) b(lda,p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(p)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  complex ( kind = 8 ) q(n,n)
  complex ( kind = 8 ) qraux(p)
  complex ( kind = 8 ) qty(n)
  complex ( kind = 8 ) qy(n)
  complex ( kind = 8 ) r(n,p)
  complex ( kind = 8 ) rsd(n)
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) work(p)
  complex ( kind = 8 ) xb(n)
  complex ( kind = 8 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  For a double complex general matrix,'
  write ( *, '(a)' ) '  ZQRDC computes the QR decomposition of a '
  write ( *, '(a)' ) '  matrix, but does not return Q and R explicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Show how Q and R can be recovered using ZQRSL.'
!
!  Set the values of the matrix A.
!
  seed = 123456789

  call c8mat_uniform_01 ( n, p, seed, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(6f8.4)' ) a(i,1:p)
  end do
!
!  Decompose the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Decompose the matrix.'

  job = 0
  ipvt(1:p) = 0

  call zqrdc ( a, lda, n, p, qraux, ipvt, work, job )
!
!  Print out what ZQRDC has stored in A...
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The packed matrix A which describes Q and R:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(6f8.4)' ) a(i,1:p)
  end do
!
!  ...and in QRAUX.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The QRAUX vector, containing some additional'
  write ( *, '(a)' ) '  information defining Q:'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,6f8.4)' )  qraux(1:n)
!
!  Print out the resulting R factor.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The R factor:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    do j = 1, p

      if ( j < i ) then
        r(i,j) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
      else
        r(i,j) = a(i,j)
      end if

    end do

    write ( *, '(2x,6f8.4)' ) r(i,1:p)

  end do
!
!  Call ZQRSL to extract the information about the Q matrix.
!  We do this, essentially, by asking ZQRSL to tell us the
!  value of Q*Y, where Y is a column of the identity matrix.
!
  job = 10000

  do j = 1, n
!
!  Set the vector Y.
!
    y(1:n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

    y(j) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
!
!  Ask ZQRSL to tell us what Q*Y is.
!
    call zqrsl ( a, lda, n, p, qraux, y, qy, qty, b, rsd, xb, job, info )

    if ( info /= 0 ) then
      write ( *, '(a,i8)' ) '  Error!  ZQRSL returns INFO = ', info
      return
    end if
!
!  Copy QY into the appropriate column of Q.
!
    q(1:n,j) = qy(1:n)

  end do
!
!  Now print out the Q matrix we have extracted.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Q factor:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,6f8.4)' ) q(i,1:n)
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
    write ( *, '(2x,6f8.4)' ) b(i,1:p)
  end do

  return
end
subroutine test28 ( )

!*****************************************************************************80
!
!! TEST28 tests ZSICO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a(n,n)
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  real    ( kind = 8 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) z(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  For a double complex symmetric matrix:'
  write ( *, '(a)' ) '  ZSICO factors the matrix and estimates'
  write ( *, '(a)' ) '  the reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  do i = 1, n
    a(i,i) = c8_uniform_01 ( seed )
    do j = i+1, n
      a(i,j) = c8_uniform_01 ( seed )
      a(j,i) = a(i,j)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call zsico ( a, lda, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition RCOND = ', rcond

  return
end
subroutine test29 ( )

!*****************************************************************************80
!
!! TEST29 tests ZSIFA and ZSISL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a(n,n)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) x(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  For a double complex symmetric matrix:'
  write ( *, '(a)' ) '  ZSIFA factors the matrix.'
  write ( *, '(a)' ) '  ZSISL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  do i = 1, n
    a(i,i) = c8_uniform_01 ( seed )
    do j = i+1, n
      a(i,j) = c8_uniform_01 ( seed )
      a(j,i) = a(i,j)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a(i,1:n)
  end do
!
!  Set the values of the right hand side vector B.
!
  call c8vec_uniform_01 ( n, seed, x )

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The right hand side B is '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2f8.4)' ) b(i)
  end do
!
!  Factor the matrix A.
!
  call zsifa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZSIFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the system.
!
  call zsisl ( a, lda, n, ipvt, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed                     Exact'
  write ( *, '(a)' ) '  Solution                     Solution'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(4g14.6)' ) b(i), x(i)
  end do

  return
end
subroutine test30 ( )

!*****************************************************************************80
!
!! TEST30 tests ZSIFA and ZSIDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a(n,n)
  complex ( kind = 8 ) a_save(n,n)
  complex ( kind = 8 ) c(n,n)
  complex ( kind = 8 ) c8_uniform_01
  complex ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) work(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  For a double complex symmetric matrix:'
  write ( *, '(a)' ) '  ZSIFA factors the matrix.'
  write ( *, '(a)' ) '  ZSIDI computes the determinant or inverse.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  do i = 1, n
    a(i,i) = c8_uniform_01 ( seed )
    do j = i+1, n
      a(i,j) = c8_uniform_01 ( seed )
      a(j,i) = a(i,j)
    end do
  end do

  a_save(1:n,1:n) = a(1:n,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call zsifa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZSIFA returned an error flag INFO = ', info
    return
  end if
!
!  Get the determinant.
!
  job = 10
  call zsidi ( a, lda, n, ipvt, det, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', real ( det(2), kind = 8 )
!
!  Get the inverse.
!
  job = 01
  call zsidi ( a, lda, n, ipvt, det, work, job )
!
!  Only the upper triangle is set, so the user must set up the
!  lower triangle:
!
  do i = 1, n
    do j = 1, i-1
      a(i,j) = a(j,i)
    end do
  end do

  c(1:n,1:n) = matmul ( a(1:n,1:n), a_save(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product inv(A) * A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) c(i,1:n)
  end do

  return
end
subroutine test31 ( )

!*****************************************************************************80
!
!! TEST31 tests ZSPCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a((n*(n+1))/2)
  complex ( kind = 8 ) a_save(n,n)
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  For a double complex symmetric matrix'
  write ( *, '(a)' ) '  in packed storage,'
  write ( *, '(a)' ) '  ZSPCO factors the matrix and estimates'
  write ( *, '(a)' ) '  the reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the packed matrix A.
!
  k = 0
  seed = 123456789

  do j = 1, n

    do i = 1, j-1
      k = k + 1
      a(k) = c8_uniform_01 ( seed )
    end do

    k = k + 1
    a(k) = c8_uniform_01 ( seed )

  end do
!
!  Copy the packed matrix into a "normal" matrix.
!
  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      a_save(i,j) = a(k)
    end do
  end do

  do j = 1, n
    do i = j+1, n
      a_save(i,j) = a_save(j,i)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a_save(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call zspco ( a, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition RCOND = ', rcond

  return
end
subroutine test32 ( )

!*****************************************************************************80
!
!! TEST32 tests ZSPFA and ZSPSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a((n*(n+1))/2)
  complex ( kind = 8 ) a_save(n,n)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32'
  write ( *, '(a)' ) '  For a double complex symmetric matrix'
  write ( *, '(a)' ) '  in packed storage,'
  write ( *, '(a)' ) '  ZSPFA factors the matrix.'
  write ( *, '(a)' ) '  ZSPSL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the packed matrix A.
!
  k = 0
  seed = 123456789

  do j = 1, n

    do i = 1, j-1
      k = k + 1
      a(k) = c8_uniform_01 ( seed )
    end do

    k = k + 1
    a(k) = c8_uniform_01 ( seed )

  end do
!
!  Copy the packed matrix into a "normal" matrix.
!
  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      a_save(i,j) = a(k)
    end do
  end do

  do j = 1, n
    do i = j+1, n
      a_save(i,j) = a_save(j,i)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a_save(i,1:n)
  end do
!
!  Set the values of the right hand side vector B.
!
  call c8vec_uniform_01 ( n, seed, x )

  b(1:n) = matmul ( a_save(1:n,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The right hand side B is '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2f8.4)' ) b(i)
  end do
!
!  Factor the matrix A.
!
  call zspfa ( a, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZSPFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the system.
!
  call zspsl ( a, n, ipvt, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed                     Exact'
  write ( *, '(a)' ) '  Solution                     Solution'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(4g14.6)' ) b(i), x(i)
  end do

  return
end
subroutine test33 ( )

!*****************************************************************************80
!
!! TEST33 tests ZSPFA and ZSPDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a((n*(n+1))/2)
  complex ( kind = 8 ) a_save(n,n)
  complex ( kind = 8 ) b_save(n,n)
  complex ( kind = 8 ) c(n,n)
  complex ( kind = 8 ) c8_uniform_01
  complex ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33'
  write ( *, '(a)' ) '  For a double complex symmetric matrix'
  write ( *, '(a)' ) '  in packed storage,'
  write ( *, '(a)' ) '  ZSPFA factors the matrix.'
  write ( *, '(a)' ) '  ZSPDI computes the determinant or inverse.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the packed matrix A.
!
  k = 0
  seed = 123456789

  do j = 1, n

    do i = 1, j-1
      k = k + 1
      a(k) = c8_uniform_01 ( seed )
    end do

    k = k + 1
    a(k) = c8_uniform_01 ( seed )

  end do
!
!  Copy the packed matrix into a "normal" matrix.
!
  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      a_save(i,j) = a(k)
    end do
  end do

  do j = 1, n
    do i = j+1, n
      a_save(i,j) = a_save(j,i)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a_save(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call zspfa ( a, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ZSPFA returned an error flag INFO = ', info
    return
  end if
!
!  Get the determinant.
!
  job = 10
  call zspdi ( a, n, ipvt, det, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', real ( det(2), kind = 8 )
!
!  Get the inverse.
!
  job = 01
  call zspdi ( a, n, ipvt, det, work, job )
!
!  Copy the packed matrix into a "normal" matrix.
!
  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      b_save(i,j) = a(k)
    end do
  end do

  do j = 1, n
    do i = j+1, n
      b_save(i,j) = b_save(j,i)
    end do
  end do

  c(1:n,1:n) = matmul ( b_save(1:n,1:n), a_save(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product inv(A) * A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) c(i,1:n)
  end do

  return
end
subroutine test34 ( )

!*****************************************************************************80
!
!! TEST34 tests ZSVDC.
!
!  Discussion:
!
!    ZSVDC computes the singular value decomposition:
!
!      A = U * S * conjg-transpose ( V )
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

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 8 ) a(m,n)
  complex ( kind = 8 ) e(max(m+1,n))
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) job
  complex ( kind = 8 ) s(max(m+1,n))
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) sigma(m,n)
  complex ( kind = 8 ) u(m,m)
  complex ( kind = 8 ) v(n,n)
  complex ( kind = 8 ) work(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST34'
  write ( *, '(a)' ) '  For an MxN matrix A in double complex general storage,'
  write ( *, '(a)' ) '  ZSVDC computes the singular value decomposition:'
  write ( *, '(a)' ) '    A = U * S * V^H'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set A.
!
  seed = 123456789

  call c8mat_uniform_01 ( m, n, seed, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,6f10.4)' ) a(i,1:n)
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

  call zsvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Warning:'
    write ( *, '(a,i8)' ) '  ZSVDC returned nonzero INFO = ', info
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Singular values:'
  write ( *, '(a)' ) ' '

  do i = 1, min ( m, n )
    write ( *, '(2x,i4,2x,2g14.6)' ) i, s(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Left Singular Vector Matrix U:'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,8f10.4)' ) u(i,1:m)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Right Singular Vector Matrix V:'
  write ( *, '(a)' ) ' '

  do i = 1,  n
    write ( *, '(2x,6f10.4)' ) v(i,1:n)
  end do

  sigma(1:m,1:n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
  do i = 1, min ( m, n )
    sigma(i,i) = s(i)
  end do

  a(1:m,1:n) = matmul ( u(1:m,1:m), &
    matmul ( sigma(1:m,1:n), transpose ( conjg ( v(1:n,1:n) ) ) ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product U * S * V^H (should equal A):'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,6f10.4)' ) a(i,1:n)
  end do

  return
end
subroutine test345 ( )

!*****************************************************************************80
!
!! TEST345 tests ZSVDC.
!
!  Discussion:
!
!    ZSVDC computes the singular value decomposition:
!
!      A = U * S * conjg-transpose ( V )
!
!    on the matrix
!
!       1  1 1  1
!      -i -1 1  i
!      -1 -1 1 -1
!       i  1 1 -i
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 4

  complex ( kind = 8 ) a(m,n)
  complex ( kind = 8 ) e(max(m+1,n))
  complex ( kind = 8 ) :: eye = cmplx ( 0.0D+00, 1.0D+00 )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  complex ( kind = 8 ) :: one = cmplx ( 1.0D+00, 0.0D+00 )
  complex ( kind = 8 ) s(max(m+1,n))
  complex ( kind = 8 ) sigma(m,n)
  complex ( kind = 8 ) u(m,m)
  complex ( kind = 8 ) v(n,n)
  complex ( kind = 8 ) work(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST345'
  write ( *, '(a)' ) '  For an MxN matrix A in double complex general storage,'
  write ( *, '(a)' ) '  ZSVDC computes the singular value decomposition:'
  write ( *, '(a)' ) '    A = U * S * V^H'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set A.
!
  a = reshape ( (/ &
    one, - eye, - one,   eye, &
    one, - one, - one,   one, &
    one,   one,   one,   one, &
    one,   eye, - one, - eye /), (/ 4, 4 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,8f10.4)' ) a(i,1:n)
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

  call zsvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Warning:'
    write ( *, '(a,i8)' ) '  ZSVDC returned nonzero INFO = ', info
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Singular values:'
  write ( *, '(a)' ) ' '

  do i = 1, min ( m, n )
    write ( *, '(2x,i4,2x,2g14.6)' ) i, s(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Left Singular Vector Matrix U:'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,8f10.4)' ) u(i,1:m)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Right Singular Vector Matrix V:'
  write ( *, '(a)' ) ' '

  do i = 1,  n
    write ( *, '(2x,8f10.4)' ) v(i,1:n)
  end do

  sigma(1:m,1:n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
  do i = 1, min ( m, n )
    sigma(i,i) = s(i)
  end do

  a(1:m,1:n) = matmul ( u(1:m,1:m), &
    matmul ( sigma(1:m,1:n), transpose ( conjg ( v(1:n,1:n) ) ) ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product U * S * V^H (should equal A):'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,8f10.4)' ) a(i,1:n)
  end do

  return
end
subroutine test35 ( )

!*****************************************************************************80
!
!! TEST35 tests ZTRCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 8 ) a(n,n)
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real    ( kind = 8 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST35'
  write ( *, '(a)' ) '  For a double complex triangular matrix,'
  write ( *, '(a)' ) '  ZTRCO estimates the condition.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  seed = 123456789

  do i = 1, n
    do j = 1, i
      a(i,j) = c8_uniform_01 ( seed )
    end do
    do j = i+1, n
      a(i,j) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    end do
  end do
!
!  Get the condition of the lower triangular matrix.
!
  job = 0
  call ztrco ( a, lda, n, rcond, z, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated reciprocal condition RCOND = ', rcond

  return
end
subroutine test36 ( )

!*****************************************************************************80
!
!! TEST36 tests ZTRDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 8 ) a(n,n)
  complex ( kind = 8 ) a_save(n,n)
  complex ( kind = 8 ) c(n,n)
  complex ( kind = 8 ) c8_uniform_01
  complex ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST36'
  write ( *, '(a)' ) '  For a double complex triangular matrix,'
  write ( *, '(a)' ) '  ZTRDI computes the determinant or inverse.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  seed = 123456789

  do i = 1, n
    do j = 1, i
      a(i,j) = c8_uniform_01 ( seed )
    end do
    do j = i+1, n
      a(i,j) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    end do
  end do

  a_save(1:n,1:n) = a(1:n,1:n)
!
!  Get the determinant of the lower triangular matrix.
!
  job = 100
  call ztrdi ( a, lda, n, det, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', real ( det(2), kind = 8 )
!
!  Get the inverse of the lower triangular matrix.
!
  job = 010
  call ztrdi ( a, lda, n, det, job, info )

  c(1:n,1:n) = matmul ( a(1:n,1:n), a_save(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product inv(A) * A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) c(i,1:n)
  end do

  return
end
subroutine test37 ( )

!*****************************************************************************80
!
!! TEST37 tests ZTRSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 8 ) a(n,n)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) seed
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST37'
  write ( *, '(a)' ) '  For a double complex triangular matrix,'
  write ( *, '(a)' ) '  ZTRSL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  seed = 123456789

  do i = 1, n
    do j = 1, i
      a(i,j) = c8_uniform_01 ( seed )
    end do
    do j = i+1, n
      a(i,j) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    end do
  end do
!
!  Set the desired solution
!
  do i = 1, n
    x(i) = cmplx ( i, 10 * i, kind = 8 )
  end do
!
!  Compute the corresponding right hand side.
!
  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )
!
!  Solve the lower triangular system.
!
  job = 0
  call ztrsl ( a, lda, n, b, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed                     Exact'
  write ( *, '(a)' ) '  Solution                     Solution'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(4g14.6)' ) b(i), x(i)
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
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
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
!    Input/output, integer SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r8_uniform_01
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
    seed = seed + i4_huge ( )
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
function c8_uniform_01 ( seed )

!*****************************************************************************80
!
!! C8_UNIFORM_01 returns a unit pseudorandom C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    The angle should be uniformly distributed between 0 and 2 * PI,
!    the square root of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 8 ) C8_UNIFORM_01, a pseudorandom complex value.
!
  implicit none

  real ( kind = 8 ) r
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta
  complex ( kind = 8 ) c8_uniform_01

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + huge ( seed )
  end if

  r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + huge ( seed )
  end if

  theta = 2.0D+00 * pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

  c8_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

  return
end
subroutine c8mat_uniform_01 ( m, n, seed, c )

!*****************************************************************************80
!
!! C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of complex ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer variable.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns
!    in the matrix.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 8 ) C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C8MAT_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do j = 1, n
    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + huge ( seed )
      end if

      r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + huge ( seed )
      end if

      theta = 2.0D+00 * pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

    end do

  end do

  return
end
subroutine c8vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of complex ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer variable.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values to compute.
!
!    Input/output, integer SEED, the "seed" value, which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, complex ( kind = 8 ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) r
  integer ( kind = 4 ) k
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) theta

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + huge ( seed )
    end if

    r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + huge ( seed )
    end if

    theta = 2.0D+00 * pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

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
