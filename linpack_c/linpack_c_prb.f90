program main

!*****************************************************************************80
!
!! MAIN is the main program for LINPACK_C_PRB.
!
!  Discussion:
!
!    LINPACK_C_PRB tests the single precision complex LINPACK routines.
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
  write ( *, '(a)' ) 'LINPACK_C_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LINPACK_C library.'

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
  write ( *, '(a)' ) 'LINPACK_C_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CCHDC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  complex ( kind = 4 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian positive definite matrix,'
  write ( *, '(a)' ) '  CCHDC computes the Cholesky decomposition.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = cmplx ( 2.5281E+00,  0.0000E+00 )
  a(2,1) = cmplx ( 2.1341E+00,  0.2147E+00 )
  a(3,1) = cmplx ( 2.4187E+00, -0.2932E+00 )

  a(1,2) = cmplx ( 2.1341E+00, -0.2147E+00 )
  a(2,2) = cmplx ( 3.0371E+00,  0.0000E+00 )
  a(3,2) = cmplx ( 2.0905E+00, -1.1505E+00 )

  a(1,3) = cmplx ( 2.4187E+00,  0.2932E+00 )
  a(2,3) = cmplx ( 2.0905E+00,  1.1505E+00 )
  a(3,3) = cmplx ( 2.7638E+00,  0.0000E+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix:'
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

  call cchdc ( a, lda, n, work, ipvt, job, info )

  if ( info /= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CCHDC returned INFO = ', info
    write ( *, '(a)' ) '  The matrix is not Hermitian positive definite.'
    return
  end if
!
!  Zero out the lower diagonal.
!
  do i = 2, n
    do j = 1, i-1
      a(i,j) = cmplx ( 0.0E+00, 0.0E+00 )
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
!! TEST02 tests CCHEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 May 2006
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

  complex ( kind = 4 ) a(lda,n)
  real    ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  complex ( kind = 4 ) s(n)
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) work(n)
  complex ( kind = 4 ) z(ldz,nz)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian positive definite matrix,'
  write ( *, '(a)' ) '  CCHEX can shift columns in a Cholesky factorization.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = cmplx ( 2.5281E+00,  0.0000E+00 )
  a(2,1) = cmplx ( 2.1341E+00,  0.2147E+00 )
  a(3,1) = cmplx ( 2.4187E+00, -0.2932E+00 )

  a(1,2) = cmplx ( 2.1341E+00, -0.2147E+00 )
  a(2,2) = cmplx ( 3.0371E+00,  0.0000E+00 )
  a(3,2) = cmplx ( 2.0905E+00, -1.1505E+00 )

  a(1,3) = cmplx ( 2.4187E+00,  0.2932E+00 )
  a(2,3) = cmplx ( 2.0905E+00,  1.1505E+00 )
  a(3,3) = cmplx ( 2.7638E+00,  0.0000E+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,6f10.4)' ) a(i,1:n)
  end do

  do i = 1, n
    z(i,1) = cmplx ( i, 0.0E+00 )
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

  call cchdc ( a, lda, n, work, ipvt, job, info )

  if ( info /= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CCHDC returned INFO = ', info
    write ( *, '(a)' ) '  This means the matrix is not positive definite.'
    return
  end if
!
!  Zero out the lower diagonal.
!
  do i = 2, n
    do j = 1, i-1
      a(i,j) = cmplx ( 0.0E+00, 0.0E+00 )
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
  write ( *, '(a,i8,a,i8)' ) '  Right circular shift columns K  = ', k, &
    ' through L = ', l

  job = 1
  call cchex ( a, lda, n, k, l, z, ldz, nz, c, s, job )
!
!  Left circular shift columns K+1 through L.
!
  k = 2
  l = 3
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  Left circular shift columns K = ', k, &
    ' through L = ', l

  job = 2
  call cchex ( a, lda, n, k, l, z, ldz, nz, c, s, job )
!
!  Print the factorization.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The shifted Cholesky factor U:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,6f10.4)' ) a(i,1:n)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The shifted vector Z:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,2g14.6)' ) z(i,1)
  end do
!
!  Compute the Cholesky product.
!
  a(1:n,1:n) = matmul ( conjg ( transpose ( a(1:n,1:n) ) ), a(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The shifted product U'' * U: '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,6f10.4)' ) a(i,1:n)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests CCHUD and CTRSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2007
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

  complex ( kind = 4 ) b(p)
  complex ( kind = 4 ) beta(p)
  real    ( kind = 4 ) c(p)
  complex ( kind = 4 ) cdotu
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  complex ( kind = 4 ) r(ldr,p)
  real    ( kind = 4 ) rho(nz)
  complex ( kind = 4 ) row(p)
  complex ( kind = 4 ) s(p)
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) x(p)
  complex ( kind = 4 ) y(nz)
  complex ( kind = 4 ) z(ldz,nz)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian matrix'
  write ( *, '(a)' ) '  CCHUD updates a Cholesky decomposition.'
  write ( *, '(a)' ) '  CTRSL solves a triangular linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we use CCHUD to solve a'
  write ( *, '(a)' ) '  least squares problem R * b = z.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is P = ', p
!
!  Initialize.
!
  r(1:p,1:p) = cmplx ( 0.0E+00, 0.0E+00 )
  z(1:p,1:nz) = cmplx ( 0.0E+00, 0.0E+00 )

  do i = 1, p
    x(i) = cmplx ( i, mod ( i, 2 ) )
  end do
!
!  Use CCHUD to form R, Z and RHO by adding X and Y a row at a time.
!  X is a row of the least squares matrix and Y the right hand side.
!
  seed = 123456789

  do i = 1, p
    call c4vec_uniform_01 ( p, seed, row )
    y(1) = cdotu ( p, row, 1, x, 1 )
    rho(1) = 0.0E+00
    call cchud ( r, ldr, p, row, z, ldz, nz, y, rho, c, s )
  end do
!
!  Generate the least squares solution, b = inverse ( R ) * Z.
!
  do j = 1, nz

    b(1:p) = z(1:p,j)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  RHS #', j
    write ( *, '(a)' ) ' '

    do i = 1, p
      if ( i <= 5 .or. p-5 < i ) then
        write ( *, '(2x,i8,2x,2g14.6)' ) i, b(i)
      end if
      if ( i == 5 ) then
        write ( *, '(a)' ) '  ......  ..............'
      end if
    end do

    job = 01

    call ctrsl ( r, ldr, p, b, job, info )

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
!! TEST04 tests CGBCO.
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

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real    ( kind = 4 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  general band storage matrix (GB):'
  write ( *, '(a)' ) '  CGBCO factors the matrix and estimates the'
  write ( *, '(a)' ) '  reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
  write ( *, '(a,i8)' ) '  The lower band is ML =  ', ml
  write ( *, '(a,i8)' ) '  The upper band is MU =  ', mu
!
!  Set the values of the matrix A.
!
  a_save(1:n,1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  m = ml + mu + 1

  seed = 123456789

  do j = 1, n
    i1 = max ( 1, j - mu )
    i2 = min ( n, j + ml )
    do i = i1, i2
      k = i - j + m
      a(k,j) = c4_uniform_01 ( seed )
      a_save(i,j) = a(k,j)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a_save(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call cgbco ( a, lda, n, ml, mu, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  Estimated reciprocal condition number = ', rcond

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests CGBFA and CGBSL.
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

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) c4_uniform_01
  complex ( kind = 4 ) b(n)
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
  complex ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  general band storage matrix (GB):'
  write ( *, '(a)' ) '  CGBFA factors the matrix;'
  write ( *, '(a)' ) '  CGBSL solves a factored linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
  write ( *, '(a,i8)' ) '  The lower band is ML =  ', ml
  write ( *, '(a,i8)' ) '  The upper band is MU =  ', mu
!
!  Set the values of the matrix A.
!
  a_save(1:n,1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  m = ml + mu + 1

  seed = 123456789

  do j = 1, n
    i1 = max ( 1, j - mu )
    i2 = min ( n, j + ml )
    do i = i1, i2
      k = i - j + m
      a(k,j) = c4_uniform_01 ( seed )
      a_save(i,j) = a(k,j)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a_save(i,1:n)
  end do
!
!  Set the values of the right hand side vector B.
!
  call c4vec_uniform_01 ( n, seed, x )

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
  call cgbfa ( a, lda, n, ml, mu, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CGBFA returned INFO = ', info
    return
  end if
!
!  Solve the system.
!
  job = 0
  call cgbsl ( a, lda, n, ml, mu, ipvt, b, job )

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
!! TEST06 tests CGBFA and CGBDI.
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

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) c4_uniform_01
  complex ( kind = 4 ) det(2)
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
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  general band storage matrix (GB):'
  write ( *, '(a)' ) '  CGBFA factors the matrix.'
  write ( *, '(a)' ) '  CGBDI computes the determinant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
  write ( *, '(a,i8)' ) '  The lower band is ML =  ', ml
  write ( *, '(a,i8)' ) '  The upper band is MU =  ', mu
!
!  Set the values of the matrix A.
!
  a_save(1:n,1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  m = ml + mu + 1

  seed = 123456789

  do j = 1, n
    i1 = max ( 1, j - mu )
    i2 = min ( n, j + ml )
    do i = i1, i2
      k = i - j + m
      a(k,j) = c4_uniform_01 ( seed )
      a_save(i,j) = a(k,j)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a_save(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call cgbfa ( a, lda, n, ml, mu, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CGBFA returned INFO = ', info
    return
  end if
!
!  Get the determinant.
!
  call cgbdi ( a, lda, n, ml, mu, ipvt, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', real ( det(2) )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests CGECO.
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

  complex ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  real    ( kind = 4 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) z(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  general storage matrix (GE):'
  write ( *, '(a)' ) '  CGECO factors the matrix and estimates the'
  write ( *, '(a)' ) '  reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  call c4mat_uniform_01 ( n, n, seed, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call cgeco ( a, lda, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  Estimated reciprocal condition number = ', rcond

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests CGEFA and CGESL.
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

  complex ( kind = 4 ) a(n,n)
  complex ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) x(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  general storage matrix (GE):'
  write ( *, '(a)' ) '  CGEFA factors the matrix.'
  write ( *, '(a)' ) '  CGESL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  call c4mat_uniform_01 ( n, n, seed, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a(i,1:n)
  end do
!
!  Set the values of the right hand side vector B.
!
  call c4vec_uniform_01 ( n, seed, x )

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The right hand side:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2f8.4)' ) b(i)
  end do
!
!  Factor the matrix A.
!
  call cgefa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CGEFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the system.
!
  job = 0
  call cgesl ( a, lda, n, ipvt, b, job )

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
!! TEST09 tests CGEFA and CGEDI.
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

  complex ( kind = 4 ) a(n,n)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) c(n,n)
  complex ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) work(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  general storage matrix (GE):'
  write ( *, '(a)' ) '  CGEFA factors the matrix.'
  write ( *, '(a)' ) '  CGEDI computes the determinant or inverse.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  call c4mat_uniform_01 ( n, n, seed, a )

  a_save(1:n,1:n) = a(1:n,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) a(i,1:n)
  end do
!
!  Factor the matrix A.
!
  call cgefa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CGEFA returned an error flag INFO = ', info
    return
  end if
!
!  Get the determinant.
!
  job = 10
  call cgedi ( a, lda, n, ipvt, det, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', real ( det(2) )
!
!  Get the inverse.
!
  job = 01
  call cgedi ( a, lda, n, ipvt, det, work, job )

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
!! TEST10 tests CGTSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2006
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) c(n)
  complex ( kind = 4 ) d(n)
  complex ( kind = 4 ) e(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  general tridiagonal matrix (GT):'
  write ( *, '(a)' ) '  CGTSL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the matrix.
!
  seed = 123456789

  c(1) = cmplx ( 0.0E+00, 0.0E+00 )
  call c4vec_uniform_01 ( n-1, seed, c(2) )

  call c4vec_uniform_01 ( n-1, seed, e(1) )
  e(n) = cmplx ( 0.0E+00, 0.0E+00 )

  d(1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  d(1:n-1) = d(1:n-1) - 2.0E+00 * e(1:n-1)
  d(2:n)   = d(2:n) - 2.0E+00 * c(2:n)
!
!  Set the desired solution
!
  do i = 1, n
    x(i) = cmplx ( i, 10 * i )
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
  call cgtsl ( n, c, d, e, b, info )

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
!! TEST11 tests CHICO.
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

  complex ( kind = 4 ) a(n,n)
  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  real    ( kind = 4 ) r4_uniform_01
  real    ( kind = 4 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) z(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian matrix (HI):'
  write ( *, '(a)' ) '  CHICO factors the matrix and estimates'
  write ( *, '(a)' ) '  the reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  do i = 1, n
    a(i,i) = cmplx ( r4_uniform_01 ( seed ), 0.0E+00 )
    do j = i+1, n
      a(i,j) = c4_uniform_01 ( seed )
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
  call chico ( a, lda, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  Estimated reciprocal condition number = ', rcond

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests CHIFA and CHISL.
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

  complex ( kind = 4 ) a(n,n)
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  real    ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) x(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian matrix (HI):'
  write ( *, '(a)' ) '  CHIFA factors the matrix.'
  write ( *, '(a)' ) '  CHISL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  do i = 1, n
    a(i,i) = cmplx ( r4_uniform_01 ( seed ), 0.0E+00 )
    do j = i+1, n
      a(i,j) = c4_uniform_01 ( seed )
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
  call c4vec_uniform_01 ( n, seed, x )

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
  call chifa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CHIFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the system.
!
  call chisl ( a, lda, n, ipvt, b )

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
!! TEST13 tests CHIFA and CHIDI.
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

  complex ( kind = 4 ) a(n,n)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) c(n,n)
  complex ( kind = 4 ) c4_uniform_01
  real    ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inert(3)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  real    ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) work(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian matrix (HI):'
  write ( *, '(a)' ) '  CHIFA factors the matrix.'
  write ( *, '(a)' ) '  CHIDI computes the determinant, inverse,'
  write ( *, '(a)' ) '  or inertia.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  do i = 1, n
    a(i,i) = cmplx ( r4_uniform_01 ( seed ), 0.0E+00 )
    do j = i+1, n
      a(i,j) = c4_uniform_01 ( seed )
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
  call chifa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CHIFA returned an error flag INFO = ', info
    return
  end if
!
!  Get the determinant.
!
  job = 010
  call chidi ( a, lda, n, ipvt, det, inert, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', det(2)
!
!  Get the inertia.
!
  job = 100
  call chidi ( a, lda, n, ipvt, det, inert, work, job )

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
  call chidi ( a, lda, n, ipvt, det, inert, work, job )
!
!  Only the upper triangle is set, so the user must set up the
!  lower triangle:
!
  do i = 1, n
    a(i,1:i-1) = conjg ( a(1:i-1,i) )
  end do

  c(1:n,1:n) = matmul ( a(1:n,1:n), a_save(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product inverse(A) * A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) c(i,1:n)
  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests CHPCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 4 ) a((n*(n+1))/2)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r4_uniform_01
  real    ( kind = 4 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian matrix using packed storage (HP),'
  write ( *, '(a)' ) '  CHPCO factors the matrix and estimates'
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
      a(k) = c4_uniform_01 ( seed )
      a_save(i,j) = a(k)
      a_save(j,i) = conjg ( a(k) )
    end do

    k = k + 1
    a(k) = cmplx ( r4_uniform_01 ( seed ), 0.0E+00 )
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
  call chpco ( a, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  Estimated reciprocal condition number = ', rcond

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests CHPFA and CHPSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 4 ) a((n*(n+1))/2)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian matrix using packed storage (HP),'
  write ( *, '(a)' ) '  CHPFA factors the matrix.'
  write ( *, '(a)' ) '  CHPSL solves a linear system.'
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
      a(k) = c4_uniform_01 ( seed )
      a_save(i,j) = a(k)
      a_save(j,i) = conjg ( a(k) )
    end do

    k = k + 1
    a(k) = cmplx ( r4_uniform_01 ( seed ), 0.0E+00 )
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
  call c4vec_uniform_01 ( n, seed, x )

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
  call chpfa ( a, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CHPFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the system.
!
  call chpsl ( a, n, ipvt, b )

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
!! TEST16 tests CHPFA and CHPDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 4 ) a((n*(n+1))/2)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) b(n,n)
  complex ( kind = 4 ) c(n,n)
  complex ( kind = 4 ) c4_uniform_01
  real    ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inert(3)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian matrix using packed storage (HP),'
  write ( *, '(a)' ) '  CHPFA factors the matrix.'
  write ( *, '(a)' ) '  CHPDI computes the determinant, inverse,'
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
      a(k) = c4_uniform_01 ( seed )
      a_save(i,j) = a(k)
      a_save(j,i) = conjg ( a(k) )
    end do

    k = k + 1
    a(k) = cmplx ( r4_uniform_01 ( seed ), 0.0E+00 )
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
  call chpfa ( a, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CHPFA returned an error flag INFO = ', info
    return
  end if
!
!  Get the determinant.
!
  job = 010
  call chpdi ( a, n, ipvt, det, inert, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', det(2)
!
!  Get the inertia.
!
  job = 100
  call chpdi ( a, n, ipvt, det, inert, work, job )

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
  call chpdi ( a, n, ipvt, det, inert, work, job )
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
!! TEST17 tests CPBCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 1
  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ), parameter :: lda = m + 1

  complex ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real    ( kind = 4 ) rcond
  complex ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  positive definite hermitian band matrix (PB),'
  write ( *, '(a)' ) '  CPBCO estimates the reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the value of the superdiagonal and diagonal.
!
  a(1,1) = cmplx ( 0.0000E+00,  0.0000E+00 )
  a(1,2) = cmplx ( 2.1341E+00, -0.2147E+00 )
  a(1,3) = cmplx ( 2.0905E+00,  1.1505E+00 )

  a(2,1) = cmplx ( 4.5281E+00,  0.0000E+00 )
  a(2,2) = cmplx ( 5.0371E+00,  0.0000E+00 )
  a(2,3) = cmplx ( 4.7638E+00,  0.0000E+00 )
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call cpbco ( a, lda, n, m, rcond, z, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CPBCO returned INFO = ', info
    write ( *, '(a)' ) '  The factorization was not completed.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  Estimated reciprocal condition number = ', rcond

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests CPBDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ), parameter :: lda = m+1

  complex ( kind = 4 ) a(lda,n)
  real    ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  positive definite hermitian band matrix (PB),'
  write ( *, '(a)' ) '  CPBDI computes the determinant as'
  write ( *, '(a)' ) '    det = MANTISSA * 10**EXPONENT'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the value of the superdiagonal and diagonal.
!
  a(1,1) = cmplx ( 0.0000E+00,  0.0000E+00 )
  a(1,2) = cmplx ( 2.1341E+00, -0.2147E+00 )
  a(1,3) = cmplx ( 2.0905E+00,  1.1505E+00 )

  a(2,1) = cmplx ( 4.5281E+00,  0.0000E+00 )
  a(2,2) = cmplx ( 5.0371E+00,  0.0000E+00 )
  a(2,3) = cmplx ( 4.7638E+00,  0.0000E+00 )

  call cpbfa ( a, lda, n, m, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  CPBFA returns INFO = ', info
    return
  end if

  call cpbdi ( a, lda, n, m, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', det(2)

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests CPBFA and CPBSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ), parameter :: lda = m+1

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  positive definite hermitian band matrix (PB),'
  write ( *, '(a)' ) '  CPBFA computes the LU factors.'
  write ( *, '(a)' ) '  CPBSL solves a factored linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the value of the superdiagonal and diagonal.
!
  a(1,1) = cmplx ( 0.0000E+00,  0.0000E+00 )
  a(1,2) = cmplx ( 2.1341E+00, -0.2147E+00 )
  a(1,3) = cmplx ( 2.0905E+00,  1.1505E+00 )

  a(2,1) = cmplx ( 4.5281E+00,  0.0000E+00 )
  a(2,2) = cmplx ( 5.0371E+00,  0.0000E+00 )
  a(2,3) = cmplx ( 4.7638E+00,  0.0000E+00 )
!
!  Set the right hand side.
!
  b(1) = cmplx (  8.7963E+00, -0.4294E+00 )
  b(2) = cmplx ( 18.4798E+00,  3.6662E+00 )
  b(3) = cmplx ( 18.4724E+00, -2.3010E+00 )
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call cpbfa ( a, lda, n, m, info )

  if ( info /= 0 ) then
    write ( *, '(a,i8)' ) '  Error!  CPBFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call cpbsl ( a, lda, n, m, b )
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
!! TEST20 tests CPOCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real    ( kind = 4 ) rcond
  complex ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian positive definite matrix (PO),'
  write ( *, '(a)' ) '  CPOCO estimates the reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = cmplx ( 2.5281E+00,  0.0000E+00 )
  a(2,1) = cmplx ( 2.1341E+00,  0.2147E+00 )
  a(3,1) = cmplx ( 2.4187E+00, -0.2932E+00 )

  a(1,2) = cmplx ( 2.1341E+00, -0.2147E+00 )
  a(2,2) = cmplx ( 3.0371E+00,  0.0000E+00 )
  a(3,2) = cmplx ( 2.0905E+00, -1.1505E+00 )

  a(1,3) = cmplx ( 2.4187E+00,  0.2932E+00 )
  a(2,3) = cmplx ( 2.0905E+00,  1.1505E+00 )
  a(3,3) = cmplx ( 2.7638E+00,  0.0000E+00 )
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition.'

  call cpoco ( a, lda, n, rcond, z, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  Estimated reciprocal condition number = ', rcond

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests CPOFA and CPODI.
!
!  Discussion:
!
!    CPOFA factors a positive definite symmetric matrix,
!    and CPODI can compute the determinant or the inverse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) a_save(lda,n)
  complex ( kind = 4 ) b(lda,n)
  real    ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian positive definite matrix (PO),'
  write ( *, '(a)' ) '  CPOFA computes the LU factors,'
  write ( *, '(a)' ) '  CPODI computes the inverse or determinant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = cmplx ( 2.5281E+00,  0.0000E+00 )
  a(2,1) = cmplx ( 2.1341E+00,  0.2147E+00 )
  a(3,1) = cmplx ( 2.4187E+00, -0.2932E+00 )

  a(1,2) = cmplx ( 2.1341E+00, -0.2147E+00 )
  a(2,2) = cmplx ( 3.0371E+00,  0.0000E+00 )
  a(3,2) = cmplx ( 2.0905E+00, -1.1505E+00 )

  a(1,3) = cmplx ( 2.4187E+00,  0.2932E+00 )
  a(2,3) = cmplx ( 2.0905E+00,  1.1505E+00 )
  a(3,3) = cmplx ( 2.7638E+00,  0.0000E+00 )

  do i = 1, n
    do j = 1, n
      a_save(i,j) = a(i,j)
    end do
  end do
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call cpofa ( a, lda, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, CPOFA returns INFO = ', info
    return
  end if
!
!  Get the determinant and inverse.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get the determinant and inverse.'

  job = 11
  call cpodi ( a, lda, n, det, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant  = ', det(1), ' * 10 ** ', det(2)
!
!  CPODI produces only the 'upper half triangle' of the inverse,
!  which is actually symmetric.  Thus, the lower half could be
!  produced by copying from the upper half.
!
  do i = 1, n
    a(i,1:i-1) = conjg ( a(1:i-1,i) )
  end do

  b(1:n,1:n) = matmul ( a(1:n,1:n), a_save(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product inverse(A) * A is '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(10f8.4)' ) b(i,1:n)
  end do

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests CPOFA and CPOSL.
!
!  Discussion:
!
!    CPOFA factors a Hermitian positive definite matrix,
!    and CPOSL can solve a factored linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  complex ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian positive definite matrix (PO),'
  write ( *, '(a)' ) '  CPOFA computes the LU factors.'
  write ( *, '(a)' ) '  CPOSL solves a factored linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  a(1,1) = cmplx ( 2.5281E+00,  0.0000E+00 )
  a(2,1) = cmplx ( 2.1341E+00,  0.2147E+00 )
  a(3,1) = cmplx ( 2.4187E+00, -0.2932E+00 )

  a(1,2) = cmplx ( 2.1341E+00, -0.2147E+00 )
  a(2,2) = cmplx ( 3.0371E+00,  0.0000E+00 )
  a(3,2) = cmplx ( 2.0905E+00, -1.1505E+00 )

  a(1,3) = cmplx ( 2.4187E+00,  0.2932E+00 )
  a(2,3) = cmplx ( 2.0905E+00,  1.1505E+00 )
  a(3,3) = cmplx ( 2.7638E+00,  0.0000E+00 )
!
!  Set the right hand side.
!
  do i = 1, n
    x(i) = cmplx ( 2 * i - 1, 2 * i  )
  end do

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call cpofa ( a, lda, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, CPOFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call cposl ( a, lda, n, b )
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
!! TEST23 tests CPPCO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 4 ) a((n*(n+1))/2)
  integer ( kind = 4 ) info
  real    ( kind = 4 ) rcond
  complex ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian positive definite packed matrix (PP),'
  write ( *, '(a)' ) '  CPPCO estimates the reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  a(1) = cmplx ( 2.5281E+00,  0.0000E+00 )

  a(2) = cmplx ( 2.1341E+00, -0.2147E+00 )
  a(3) = cmplx ( 3.0371E+00,  0.0000E+00 )

  a(4) = cmplx ( 2.4187E+00,  0.2932E+00 )
  a(5) = cmplx ( 2.0905E+00,  1.1505E+00 )
  a(6) = cmplx ( 2.7638E+00,  0.0000E+00 )
!
!  Estimate the condition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the condition number.'

  call cppco ( a, n, rcond, z, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  Estimated reciprocal condition number = ', rcond

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests CPPFA and CPPDI.
!
!  Discussion:
!
!    CPPFA factors a Hermitian positive definite packed matrix,
!    and CPPDI can compute the determinant or the inverse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 4 ) a((n*(n+1))/2)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) b(n,n)
  complex ( kind = 4 ) c(n,n)
  real    ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian positive definite packed matrix (PP),'
  write ( *, '(a)' ) '  CPPFA factors the matrix.'
  write ( *, '(a)' ) '  CPPDI computes the inverse or determinant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  a(1) = cmplx ( 2.5281E+00,  0.0000E+00 )

  a(2) = cmplx ( 2.1341E+00, -0.2147E+00 )
  a(3) = cmplx ( 3.0371E+00,  0.0000E+00 )

  a(4) = cmplx ( 2.4187E+00,  0.2932E+00 )
  a(5) = cmplx ( 2.0905E+00,  1.1505E+00 )
  a(6) = cmplx ( 2.7638E+00,  0.0000E+00 )

  a_save(1,1) = a(1)

  a_save(1,2) = a(2)
  a_save(2,1) = conjg ( a(2) )
  a_save(2,2) = a(3)

  a_save(1,3) = a(4)
  a_save(3,1) = conjg ( a(4) )
  a_save(2,3) = a(5)
  a_save(3,2) = conjg ( a(5) )
  a_save(3,3) = a(6)
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call cppfa ( a, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, CPPFA returns INFO = ', info
    return
  end if
!
!  Invert the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get the determinant and inverse.'

  job = 11
  call cppdi ( a, n, det, job )
!
!  Print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  Determinant  = ', det(1), ' * 10 ** ', det(2)
!
!  CPPDI produces only the 'upper half triangle' of the inverse,
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
  write ( *, '(a)' ) '  Matrix Inverse(A):'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,6f10.4)' ) b(i,1:n)
  end do

  c(1:n,1:n) = matmul ( b(1:n,1:n), a_save(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix Inverse(A) * A:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,6f10.4)' ) c(i,1:n)
  end do

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests CPPFA and CPPSL.
!
!  Discussion:
!
!    CPOFA factors a Hermitian positive definite packed matrix,
!    and CPOSL can solve a factored linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 4 ) a((n*(n+1))/2)
  complex ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  complex ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian positive definite packed matrix (PP),'
  write ( *, '(a)' ) '  CPPFA factors the matrix.'
  write ( *, '(a)' ) '  CPPSL solves a factored linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  a(1) = cmplx ( 2.5281E+00,  0.0000E+00 )

  a(2) = cmplx ( 2.1341E+00, -0.2147E+00 )
  a(3) = cmplx ( 3.0371E+00,  0.0000E+00 )

  a(4) = cmplx ( 2.4187E+00,  0.2932E+00 )
  a(5) = cmplx ( 2.0905E+00,  1.1505E+00 )
  a(6) = cmplx ( 2.7638E+00,  0.0000E+00 )

  b(1) = cmplx ( 20.12350E+00, 28.92670E+00 )
  b(2) = cmplx ( 14.36550E+00, 34.92680E+00 )
  b(3) = cmplx ( 27.69760E+00, 26.03750E+00 )
!
!  Factor the matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix.'

  call cppfa ( a, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Error, CPPFA returns INFO = ', info
    return
  end if
!
!  Solve the linear system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the linear system.'

  call cppsl ( a, n, b )
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
!! TEST26 tests CPTSL.
!
!  Discussion:
!
!    CPTSL factors and solves a Hermitian positive definite
!    tridiagonal system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) d(n)
  complex ( kind = 4 ) e(n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  Hermitian positive definite tridiagonal matrix (PT),'
  write ( *, '(a)' ) '  CPTSL factors and solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the value of the superdiagonal and diagonal.
!
  e(1) = cmplx ( 2.1341E+00, -0.2147E+00 )
  e(2) = cmplx ( 2.0905E+00,  1.1505E+00 )
  e(3) = cmplx ( 0.0000E+00,  0.0000E+00 )

  d(1) = cmplx ( 4.5281E+00,  0.0000E+00 )
  d(2) = cmplx ( 5.0371E+00,  0.0000E+00 )
  d(3) = cmplx ( 4.7638E+00,  0.0000E+00 )
!
!  Set the right hand side.
!
  b(1) = cmplx (  8.7963E+00, -0.4294E+00 )
  b(2) = cmplx ( 18.4798E+00,  3.6662E+00 )
  b(3) = cmplx ( 18.4724E+00, -2.3010E+00 )
!
!  Factor and solve the system.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Factor the matrix and solve the system.'

  call cptsl ( n, d, e, b )
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
!! TEST27 tests CQRDC and CQRSL.
!
!  Discussion:
!
!    CQRDC and CQRSL compute the QR factorization, and use it
!    to solve linear systems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: p = 3
  integer ( kind = 4 ), parameter :: lda = n

  complex ( kind = 4 ) a(lda,p)
  complex ( kind = 4 ) b(lda,p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(p)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  complex ( kind = 4 ) q(n,n)
  complex ( kind = 4 ) qraux(p)
  complex ( kind = 4 ) qty(n)
  complex ( kind = 4 ) qy(n)
  complex ( kind = 4 ) r(n,p)
  complex ( kind = 4 ) rsd(n)
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) work(p)
  complex ( kind = 4 ) xb(n)
  complex ( kind = 4 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  general matrix,'
  write ( *, '(a)' ) '  CQRDC computes the QR decomposition of a '
  write ( *, '(a)' ) '  matrix, but does not return Q and R explicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix row order is N    = ', n
  write ( *, '(a,i8)' ) '  The matrix column order is P = ', p
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Show how Q and R can be recovered using CQRSL.'
!
!  Set the values of the matrix A.
!
  seed = 123456789

  call c4mat_uniform_01 ( n, p, seed, a )

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

  call cqrdc ( a, lda, n, p, qraux, ipvt, work, job )
!
!  Print out what CQRDC has stored in A...
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

  do i = 1, n
    write ( *, '(2x,2f8.4)' )  qraux(i)
  end do
!
!  Print out the resulting R factor.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The R factor:'
  write ( *, '(a)' ) ' '

  do i = 1, n

    r(i,1:i-1) = cmplx ( 0.0E+00, 0.0E+00 )
    r(i,i:p) = a(i,i:p)

    write ( *, '(2x,6f8.4)' ) r(i,1:p)

  end do
!
!  Call CQRSL to extract the information about the Q matrix.
!  We do this, essentially, by asking CQRSL to tell us the
!  value of Q*Y, where Y is a column of the identity matrix.
!
  job = 10000

  do j = 1, n
!
!  Set the vector Y.
!
    y(1:n) = cmplx ( 0.0E+00, 0.0E+00 )

    y(j) = cmplx ( 1.0E+00, 0.0E+00 )
!
!  Ask CQRSL to tell us what Q*Y is.
!
    call cqrsl ( a, lda, n, p, qraux, y, qy, qty, b, rsd, xb, job, info )

    if ( info /= 0 ) then
      write ( *, '(a,i8)' ) '  Error!  CQRSL returns INFO = ', info
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
!! TEST28 tests CSICO.
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

  complex ( kind = 4 ) a(n,n)
  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  real    ( kind = 4 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) z(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  symmetric matrix (SI):'
  write ( *, '(a)' ) '  CSICO factors the matrix and estimates'
  write ( *, '(a)' ) '  the reciprocal condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  do i = 1, n
    a(i,i) = c4_uniform_01 ( seed )
    do j = i+1, n
      a(i,j) = c4_uniform_01 ( seed )
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
  call csico ( a, lda, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  Estimated reciprocal condition number = ', rcond

  return
end
subroutine test29 ( )

!*****************************************************************************80
!
!! TEST29 tests CSIFA and CSISL.
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

  complex ( kind = 4 ) a(n,n)
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) x(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  symmetric matrix (SI):'
  write ( *, '(a)' ) '  CSIFA factors the matrix.'
  write ( *, '(a)' ) '  CSISL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  do i = 1, n
    a(i,i) = c4_uniform_01 ( seed )
    do j = i+1, n
      a(i,j) = c4_uniform_01 ( seed )
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
  call c4vec_uniform_01 ( n, seed, x )

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
  call csifa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CSIFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the system.
!
  call csisl ( a, lda, n, ipvt, b )

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
!! TEST30 tests CSIFA and CSIDI.
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

  complex ( kind = 4 ) a(n,n)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) c(n,n)
  complex ( kind = 4 ) c4_uniform_01
  complex ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) work(n)

  lda = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  symmetric matrix (SI):'
  write ( *, '(a)' ) '  CSIFA factors the matrix.'
  write ( *, '(a)' ) '  CSIDI computes the determinant or inverse.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the values of the matrix A.
!
  seed = 123456789

  do i = 1, n
    a(i,i) = c4_uniform_01 ( seed )
    do j = i+1, n
      a(i,j) = c4_uniform_01 ( seed )
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
  call csifa ( a, lda, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CSIFA returned an error flag INFO = ', info
    return
  end if
!
!  Get the determinant.
!
  job = 10
  call csidi ( a, lda, n, ipvt, det, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', real ( det(2) )
!
!  Get the inverse.
!
  job = 01
  call csidi ( a, lda, n, ipvt, det, work, job )
!
!  Only the upper triangle is set, so the user must set up the
!  lower triangle:
!
  do i = 1, n
    a(i,1:i-1) = a(1:i-1,i)
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
!! TEST31 tests CSPCO.
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

  complex ( kind = 4 ) a((n*(n+1))/2)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 4 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  symmetric matrix in packed storage (SP),'
  write ( *, '(a)' ) '  CSPCO factors the matrix and estimates'
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
      a(k) = c4_uniform_01 ( seed )
    end do

    k = k + 1
    a(k) = c4_uniform_01 ( seed )

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
    a_save(j+1:n,j) = a_save(j,j+1:n)
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
  call cspco ( a, n, ipvt, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  Estimated reciprocal condition number = ', rcond

  return
end
subroutine test32 ( )

!*****************************************************************************80
!
!! TEST32 tests CSPFA and CSPSL.
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

  complex ( kind = 4 ) a((n*(n+1))/2)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  symmetric matrix in packed storage (SP),'
  write ( *, '(a)' ) '  CSPFA factors the matrix.'
  write ( *, '(a)' ) '  CSPSL solves a linear system.'
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
      a(k) = c4_uniform_01 ( seed )
    end do

    k = k + 1
    a(k) = c4_uniform_01 ( seed )

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
    a_save(j+1:n,j) = a_save(j,j+1:n)
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
  call c4vec_uniform_01 ( n, seed, x )

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
  call cspfa ( a, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CSPFA returned an error flag INFO = ', info
    return
  end if
!
!  Solve the system.
!
  call cspsl ( a, n, ipvt, b )

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
!! TEST33 tests CSPFA and CSPDI.
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

  complex ( kind = 4 ) a((n*(n+1))/2)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) b_save(n,n)
  complex ( kind = 4 ) c(n,n)
  complex ( kind = 4 ) c4_uniform_01
  complex ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  symmetric matrix in packed storage (SP),'
  write ( *, '(a)' ) '  CSPFA factors the matrix.'
  write ( *, '(a)' ) '  CSPDI computes the determinant or inverse.'
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
      a(k) = c4_uniform_01 ( seed )
    end do

    k = k + 1
    a(k) = c4_uniform_01 ( seed )

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
    a_save(j+1:n,j) = a_save(j,j+1:n)
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
  call cspfa ( a, n, ipvt, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  CSPFA returned an error flag INFO = ', info
    return
  end if
!
!  Get the determinant.
!
  job = 10
  call cspdi ( a, n, ipvt, det, work, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', real ( det(2) )
!
!  Get the inverse.
!
  job = 01
  call cspdi ( a, n, ipvt, det, work, job )
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
    b_save(j+1:n,j) = b_save(j,j+1:n)
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

!***********************************************************************
!
!! TEST34 tests CSVDC.
!
!  Discussion:
!
!    CSVDC computes the singular value decomposition:
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

  complex ( kind = 4 ) a(m,n)
  complex ( kind = 4 ) e(max(m+1,n))
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) job
  complex ( kind = 4 ) s(max(m+1,n))
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) sigma(m,n)
  complex ( kind = 4 ) u(m,m)
  complex ( kind = 4 ) v(n,n)
  complex ( kind = 4 ) work(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST34'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  general storage matrix,'
  write ( *, '(a)' ) '  CSVDC computes the singular value decomposition:'
  write ( *, '(a)' ) '    A = U * S * V^H'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix row order is M =    ', m
  write ( *, '(a,i8)' ) '  The matrix column order is N = ', n
!
!  Set A.
!
  seed = 123456789

  call c4mat_uniform_01 ( m, n, seed, a )

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

  call csvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Warning:'
    write ( *, '(a,i8)' ) '  CSVDC returned nonzero INFO = ', info
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

  sigma(1:m,1:n) = cmplx ( 0.0E+00, 0.0E+00 )
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
!! TEST345 tests CSVDC.
!
!  Discussion:
!
!    CSVDC computes the singular value decomposition:
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

  complex ( kind = 4 ) a(m,n)
  complex ( kind = 4 ) e(max(m+1,n))
  complex ( kind = 4 ) :: eye = cmplx ( 0.0E+00, 1.0E+00 )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  complex ( kind = 4 ) :: one = cmplx ( 1.0E+00, 0.0E+00 )
  complex ( kind = 4 ) s(max(m+1,n))
  complex ( kind = 4 ) sigma(m,n)
  complex ( kind = 4 ) u(m,m)
  complex ( kind = 4 ) v(n,n)
  complex ( kind = 4 ) work(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST345'
  write ( *, '(a)' ) '  For an MxN matrix A in complex general storage,'
  write ( *, '(a)' ) '  CSVDC computes the singular value decomposition:'
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

  call csvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Warning:'
    write ( *, '(a,i8)' ) '  CSVDC returned nonzero INFO = ', info
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

  sigma(1:m,1:n) = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )
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
!! TEST35 tests CTRCO.
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

  complex ( kind = 4 ) a(n,n)
  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real    ( kind = 4 ) rcond
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST35'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  triangular matrix (TR),'
  write ( *, '(a)' ) '  CTRCO estimates the condition.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the matrix.
!
  seed = 123456789

  do i = 1, n
    do j = 1, i
      a(i,j) = c4_uniform_01 ( seed )
    end do
    a(i,i+1:n) = cmplx ( 0.0E+00, 0.0E+00 )
  end do
!
!  Get the condition of the lower triangular matrix.
!
  job = 0
  call ctrco ( a, lda, n, rcond, z, job )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) '  Estimated reciprocal condition number = ', rcond

  return
end
subroutine test36 ( )

!*****************************************************************************80
!
!! TEST36 tests CTRDI.
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

  complex ( kind = 4 ) a(n,n)
  complex ( kind = 4 ) a_save(n,n)
  complex ( kind = 4 ) c(n,n)
  complex ( kind = 4 ) c4_uniform_01
  complex ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST36'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  triangular matrix (TR),'
  write ( *, '(a)' ) '  CTRDI computes the determinant or inverse.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the matrix.
!
  seed = 123456789

  do i = 1, n
    do j = 1, i
      a(i,j) = c4_uniform_01 ( seed )
    end do
    a(i,i+1:n) = cmplx ( 0.0E+00, 0.0E+00 )
  end do

  a_save(1:n,1:n) = a(1:n,1:n)
!
!  Get the determinant of the lower triangular matrix.
!
  job = 100
  call ctrdi ( a, lda, n, det, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6,a,g14.6)' ) &
    '  Determinant = ', det(1), ' * 10** ', real ( det(2) )
!
!  Get the inverse of the lower triangular matrix.
!
  job = 010
  call ctrdi ( a, lda, n, det, job, info )

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
!! TEST37 tests CTRSL.
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

  complex ( kind = 4 ) a(n,n)
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) seed
  complex ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST37'
  write ( *, '(a)' ) '  For a single precision complex (C)'
  write ( *, '(a)' ) '  triangular matrix (TR),'
  write ( *, '(a)' ) '  CTRSL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The matrix order is N = ', n
!
!  Set the matrix.
!
  seed = 123456789

  do i = 1, n
    do j = 1, i
      a(i,j) = c4_uniform_01 ( seed )
    end do
    a(i,i+1:n) = cmplx ( 0.0E+00, 0.0E+00 )
  end do
!
!  Set the desired solution
!
  do i = 1, n
    x(i) = cmplx ( i, 10 * i )
  end do
!
!  Compute the corresponding right hand side.
!
  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )
!
!  Solve the lower triangular system.
!
  job = 0
  call ctrsl ( a, lda, n, b, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed                     Exact'
  write ( *, '(a)' ) '  Solution                     Solution'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(4g14.6)' ) b(i), x(i)
  end do

  return
end
function c4_uniform_01 ( seed )

!*****************************************************************************80
!
!! C4_UNIFORM_01 returns a unit pseudorandom C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
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
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 4 ) C4_UNIFORM_01, a pseudorandom complex value.
!
  implicit none

  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  real    ( kind = 4 ), parameter :: pi = 3.1415926E+00
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  real    ( kind = 4 ) theta

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C4_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  r = sqrt ( real ( seed, kind = 4 ) * 4.656612875E-10 )

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  theta = 2.0E+00 * pi * ( real ( seed, kind = 4 ) * 4.656612875E-10 )

  c4_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ), kind = 4 )

  return
end
subroutine c4mat_uniform_01 ( m, n, seed, c )

!*****************************************************************************80
!
!! C4MAT_UNIFORM_01 returns a unit pseudorandom C4MAT.
!
!  Discussion:
!
!    A C4MAT is a matrix of complex ( kind = 4 ) values.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 4 ) C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) c(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) j
  real    ( kind = 4 ) r
  integer ( kind = 4 ) k
  real    ( kind = 4 ), parameter :: pi = 3.1415926E+00
  integer ( kind = 4 ) seed
  real    ( kind = 4 ) theta

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C4MAT_UNIFORM_01 - Fatal error!'
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

      r = sqrt ( real ( seed, kind = 4 ) * 4.656612875E-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge ( )
      end if

      theta = 2.0D+00 * pi * ( real ( seed, kind = 4 ) * 4.656612875E-10 )

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 4 )

    end do

  end do

  return
end
subroutine c4vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C4VEC_UNIFORM_01 returns a unit pseudorandom C4VEC.
!
!  Discussion:
!
!    A C4VEC is a vector of complex ( kind = 4 ) values.
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
!    Input, integer ( kind = 4 ) N, the number of values to compute.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 4 ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  real    ( kind = 4 ), parameter :: pi = 3.1415926E+00
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  real    ( kind = 4 ) theta

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C4VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge ( )
    end if

    r = sqrt ( real ( seed, kind = 4 ) * 4.656612875E-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge ( )
    end if

    theta = 2.0E+00 * pi * ( real ( seed, kind = 4 ) * 4.656612875E-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 4 )

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
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end
function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r4_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer ( kind = 4 ) arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
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
!    11 August 2004
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
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 4 ) r4_uniform_01

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

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
