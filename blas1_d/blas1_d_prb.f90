program main

!*****************************************************************************80
!
!! MAIN is the main program for BLAS1_D_PRB.
!
!  Discussion:
!
!    BLAS1_D_PRB tests the BLAS1 double precision real routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLAS1_D_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BLAS1_D library.'

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
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLAS1_D_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DASUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: lda = 6
  integer ( kind = 4 ), parameter :: ma = 5
  integer ( kind = 4 ), parameter :: na = 4
  integer ( kind = 4 ), parameter :: nx = 10

  real ( kind = 8 ) a(lda,na)
  real ( kind = 8 ) dasum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(nx)

  do i = 1, nx
    x(i) = ( -1.0D+00 )**i * real ( 2 * i, kind = 8 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  DASUM adds the absolute values of elements'
  write ( *, '(a)' ) '  of a double precision real vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, nx
    write ( *, '(2x,i6,g14.6)' ) i, x(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  DASUM ( NX, X, 1 ) =   ', dasum ( nx, x, 1 )
  write ( *, '(a,g14.6)' ) '  DASUM ( NX/2, X, 2 ) = ', dasum ( nx/2, x, 2 )
  write ( *, '(a,g14.6)' ) '  DASUM ( 2, X, NX/2 ) = ', dasum ( 2, x, nx/2 )

  a(1:lda,1:na) = 0.0D+00

  do i = 1, ma
    do j = 1, na
      a(i,j) = ( -1.0D+00 )**( i + j ) * real ( 10 * i + j, kind = 8 )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate with a matrix A:'
  write ( *, '(a)' ) ' '
  do i = 1, ma
    write ( *, '(2x,5g14.6)' ) a(i,1:na)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  DASUM(MA,A(1,2),1) =   ', dasum ( ma, a(1,2), 1 )
  write ( *, '(a,g14.6)' ) &
    '  DASUM(NA,A(2,1),LDA) = ', dasum ( na, a(2,1), lda )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests DAXPY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ) da
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    y(i) = real ( 100 * i, kind = 8 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  DAXPY adds a multiple of a double precision real'
  write ( *, '(a)' ) '  vector X to vector Y.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, x(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Y = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  da = 1.0D+00
  call daxpy ( n, da, x, 1, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  DAXPY ( N, ', da, ', X, 1, Y, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  do i = 1, n
    y(i) = real ( 100 * i, kind = 8 )
  end do

  da = -2.0D+00
  call daxpy ( n, da, x, 1, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  DAXPY ( N, ', da, ', X, 1, Y, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  do i = 1, n
    y(i) = real ( 100 * i, kind = 8 )
  end do

  da = +3.0D+00
  call daxpy ( 3, da, x, 2, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  DAXPY ( 3, ', da, ', X, 2, Y, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  do i = 1, n
    y(i) = real ( 100 * i, kind = 8 )
  end do

  da = -4.0D+00
  call daxpy ( 3, da, x, 1, y, 2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  DAXPY ( 3, ', da, ', X, 1, Y, 2 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests DCOPY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(5,5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(10)
  real ( kind = 8 ) y(10)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  DCOPY copies one double precision real vector'
  write ( *, '(a)' ) '  into another.'

  do i = 1, 10
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, 10
    y(i) = real ( 10 * i, kind = 8 )
  end do

  do i = 1, 5
    do j = 1, 5
      a(i,j) = real ( 10 * i + j, kind = 8 )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i6,g14.6)' ) i, x(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Y = '
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A = '
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,5f8.2)' ) a(i,1:5)
  end do

  call dcopy ( 5, x, 1, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  DCOPY ( 5, X, 1, Y, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  do i = 1, 10
    y(i) = real ( 10 * i, kind = 8 )
  end do

  call dcopy ( 3, x, 2, y, 3 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  DCOPY ( 3, X, 2, Y, 3 )'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  call dcopy ( 5, x, 1, a, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  DCOPY ( 5, X, 1, A, 1 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A = '
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,5f8.2)' ) a(i,1:5)
  end do

  do i = 1, 5
    do j = 1, 5
      a(i,j) = real ( 10 * i + j, kind = 8 )
    end do
  end do

  call dcopy ( 5, x, 2, a, 5 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  DCOPY ( 5, X, 2, A, 5 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A = '
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,5f8.2)' ) a(i,1:5)
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests DDOT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = 10
  integer ( kind = 4 ), parameter :: ldb = 7
  integer ( kind = 4 ), parameter :: ldc = 6

  real ( kind = 8 ) a(lda,lda)
  real ( kind = 8 ) b(ldb,ldb)
  real ( kind = 8 ) c(ldc,ldc)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) ddot
  real ( kind = 8 ) sum1
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  DDOT computes the dot product of two'
  write ( *, '(a)' ) '  double precision real vectors.'

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    y(i) = - real ( i, kind = 8 )
  end do

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( i + j, kind = 8 )
    end do
  end do

  do i = 1, n
    do j = 1, n
      b(i,j) = real ( i - j, kind = 8 )
    end do
  end do
!
!  To compute a simple dot product of two vectors, use a
!  call like this:
!
  sum1 = ddot ( n, x, 1, y, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Dot product of X and Y is ', sum1
!
!  To multiply a ROW of a matrix A times a vector X, we need to
!  specify the increment between successive entries of the row of A:
!
  sum1 = ddot ( n, a(2,1), lda, x, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Product of row 2 of A and X is ', sum1
!
!  Product of a column of A and a vector is simpler:
!
  sum1 = ddot ( n, a(1,2), 1, x, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Product of column 2 of A and X is ', sum1
!
!  Here's how matrix multiplication, c = a*b, could be done
!  with DDOT:
!
  do i = 1, n
    do j = 1, n
      c(i,j) = ddot ( n, a(i,1), lda, b(1,j), 1 )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix product computed with DDOT:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,5g14.6)' ) c(i,1:n)
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests DMACH.
!
!  Discussion:
!
!    The DMACH routine is not part of the official BLAS release.
!    It was used for the testing routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) dmach
  integer ( kind = 4 ) job

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  DMACH computes several machine-dependent'
  write ( *, '(a)' ) '  double precision real arithmetic parameters.'

  write ( *, '(a)' ) ' '
  write ( *, * ) '  DMACH(1)  = machine epsilon = ', dmach ( 1 )
  write ( *, * ) '  DMACH(2)  = a tiny value    = ', dmach ( 2 )
  write ( *, * ) '  DMACH(3)  = a huge value    = ', dmach ( 3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FORTRAN90 parameters:'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  EPSILON() = machine epsilon = ', epsilon ( 1.0D+00 )
  write ( *, * ) '  TINY()    = a tiny value    = ', tiny ( 1.0D+00 )
  write ( *, * ) '  HUGE()    = a huge value    = ', huge ( 1.0D+00 )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests DNRM2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lda = n + 5
!
!  These parameters illustrate the fact that matrices are typically
!  dimensioned with more space than the user requires.
!
  real ( kind = 8 ) a(lda,lda)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) j
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) sum1
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  DNRM2 computes the Euclidean norm of a '
  write ( *, '(a)' ) '  double precision real vector.'
!
!  Compute the euclidean norm of a vector:
!
  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The vector X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,f8.4)' ) i, x(i)
  end do
  incx = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The 2-norm of X is ', dnrm2 ( n, x, incx )
!
!  Compute the euclidean norm of a row or column of a matrix:
!
  do i = 1, n
    do j = 1, n
      a(i,j) = real ( i + j, kind = 8 )
    end do
  end do

  incx = lda
  sum1 = dnrm2 ( n, a(2,1), incx )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The 2-norm of row 2 of A is ', sum1

  incx = 1
  sum1 = dnrm2 ( n, a(1,2), incx )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The 2-norm of column 2 of A is ', sum1

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests DROT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  real ( kind = 8 ) s
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    y(i) = real ( i * i - 12, kind = 8 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) &
    '  DROT carries out a double precision real Givens rotation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X and Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
  end do

  c = 0.5D+00
  s = sqrt ( 1.0D+00 - c * c )
  call drot ( n, x, 1, y, 1, c, s )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a,f8.4,a)' ) '  DROT ( N, X, 1, Y, 1, ', c, ',', s, ' )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
  end do

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    y(i) = real ( i * i - 12, kind = 8 )
  end do

  c = x(1) / sqrt ( x(1) * x(1) + y(1) * y(1) )
  s = y(1) / sqrt ( x(1) * x(1) + y(1) * y(1) )
  call drot ( n, x, 1, y, 1, c, s )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a,f8.4,a)' ) '  DROT ( N, X, 1, Y, 1, ', c, ',', s, ' )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests DROTG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) s
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  DROTG generates a double precision real Givens rotation'
  write ( *, '(a)' ) '    (  C  S ) * ( A ) = ( R )'
  write ( *, '(a)' ) '    ( -S  C )   ( B )   ( 0 )'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    a = r8_uniform_01 ( seed )
    b = r8_uniform_01 ( seed )

    sa = a
    sb = b

    call drotg ( sa, sb, c, s )

    r = sa
    z = sb

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6,a,g14.6)' ) '  A =  ', a,  '  B =  ', b
    write ( *, '(a,g14.6,a,g14.6)' ) '  C =  ', c,  '  S =  ', s
    write ( *, '(a,g14.6,a,g14.6)' ) '  R =  ', r,  '  Z =  ', z
    write ( *, '(a,g14.6)' ) '   C*A+S*B = ',  c * a + s * b
    write ( *, '(a,g14.6)' ) '  -S*A+C*B = ', -s * a + c * b

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests DSCAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ) da
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  DSCAL multiplies a double precision real scalar '
  write ( *, '(a)' ) '  times a double precision real vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, x(i)
  end do

  da = 5.0D+00
  call dscal ( n, da, x, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  DSCAL ( N, ', da, ', X, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, x(i)
  end do

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  da = -2.0D+00
  call dscal ( 3, da, x, 2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  DSCAL ( 3, ', da, ', X, 2 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, x(i)
  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests DSWAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    y(i) = real ( 100 * i, kind = 8 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  DSWAP swaps two double precision real vectors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X and Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
  end do

  call dswap ( n, x, 1, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DSWAP ( N, X, 1, Y, 1 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X and Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
  end do

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    y(i) = real ( 100 * i, kind = 8 )
  end do

  call dswap ( 3, x, 2, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  DSWAP ( 3, X, 2, Y, 1 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X and Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests IDAMAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) idamax
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  IDAMAX returns the index of the entry of'
  write ( *, '(a)' ) '  maximum magnitude in a double precision real vector.'

  do i = 1, n
    x(i) = real ( mod ( 7 * i, 11 ), kind = 8 ) - real ( n / 2, kind = 8 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The vector X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,f8.4)' ) i, x(i)
  end do

  incx = 1

  i1 = idamax ( n, x, incx )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The index of maximum magnitude = ', i1

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests IDAMAX, DAXPY and DSCAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  Use IDAMAX, DAXPY and DSCAL'
  write ( *, '(a)' ) '  in a Gauss elimination routine.'
!
!  Set the matrix.
!
  do i = 1, n
    do j = 1, n

      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( i == j + 1 ) then
        a(i,j) = - 1.0D+00
      else if ( i == j - 1 ) then
        a(i,j) = - 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do
!
!  Set the right hand side.
!
  b(1:n-1) = 0.0D+00
  b(n) = real ( n, kind = 8 ) + 1.0D+00

  info = 0

  do k = 1, n - 1

    l = idamax ( n-k+1, a(k,k), 1 ) + k - 1
    ipvt(k) = l

    if ( a(l,k) == 0.0D+00 ) then

      info = k

    else

      if ( l /= k ) then
        t = a(l,k)
        a(l,k) = a(k,k)
        a(k,k) = t
      end if

      t = -1.0D+00 / a(k,k)
      call dscal ( n-k, t, a(k+1,k), 1 )

      do j = k+1, n
        t = a(l,j)
        if ( l /= k ) then
          a(l,j) = a(k,j)
          a(k,j) = t
        end if
        call daxpy ( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )
      end do

    end if

  end do

  ipvt(n) = n
  if ( a(n,n) == 0.0D+00 ) then
    info = n
  end if

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The matrix is singular.'
    return
  end if

  do k = 1, n-1
    l = ipvt(k)
    t = b(l)
    if ( l /= k ) then
      b(l) = b(k)
      b(k) = t
    end if
    call daxpy ( n-k, t, a(k+1,k), 1, b(k+1), 1 )
  end do

  do k = n, 1, -1
    b(k) = b(k) / a(k,k)
    t = - b(k)
    call daxpy ( k-1, t, a(1,k), 1, b(1), 1 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First five entries of solution:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) b(1:5)

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
!    Philip Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r8_uniform_01

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
