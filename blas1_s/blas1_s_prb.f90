program main

!*****************************************************************************80
!
!! MAIN is the main program for BLAS1_PRB.
!
!  Discussion:
!
!    BLAS1_S_PRB tests the BLAS1 single precision real routines.
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
  write ( *, '(a)' ) 'BLAS1_S_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BLAS1_S library.'

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
  write ( *, '(a)' ) 'BLAS1_S_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests ISAMAX.
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
  integer ( kind = 4 ) isamax
  real ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  ISAMAX returns the index of the entry of '
  write ( *, '(a)' ) '  maximum magnitude in a single precision real vector.'

  do i = 1, n
    x(i) = real ( mod ( 7 * i, 11 ), kind = 4 ) - real ( n / 2, kind = 4 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The vector X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,f8.4)' ) i, x(i)
  end do

  incx = 1

  i1 = isamax ( n, x, incx )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The index of maximum magnitude = ', i1

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests ISAMAX, SAXPY and SSCAL.
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

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isamax
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 4 ) t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Use ISAMAX, SAXPY and SSCAL'
  write ( *, '(a)' ) '  in a Gauss elimination routine.'
!
!  Set the matrix.
!
  do i = 1, n
    do j = 1, n

      if ( i == j ) then
        a(i,j) = 2.0E+00
      else if ( i == j + 1 ) then
        a(i,j) = - 1.0E+00
      else if ( i == j - 1 ) then
        a(i,j) = - 1.0E+00
      else
        a(i,j) = 0.0E+00
      end if

    end do
  end do
!
!  Set the right hand side.
!
  b(1:n-1) = 0.0E+00
  b(n) = real ( n, kind = 4 ) + 1.0E+00

  info = 0

  do k = 1, n - 1

    l = isamax ( n-k+1, a(k,k), 1 ) + k - 1
    ipvt(k) = l

    if ( a(l,k) == 0.0E+00 ) then

      info = k

    else

      if ( l /= k ) then
        t = a(l,k)
        a(l,k) = a(k,k)
        a(k,k) = t
      end if

      t = -1.0E+00 / a(k,k)
      call sscal ( n-k, t, a(k+1,k), 1 )

      do j = k+1, n
        t = a(l,j)
        if ( l /= k ) then
          a(l,j) = a(k,j)
          a(k,j) = t
        end if
        call saxpy ( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )
      end do

    end if

  end do

  ipvt(n) = n
  if ( a(n,n) == 0.0E+00 ) then
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
    call saxpy ( n-k, t, a(k+1,k), 1, b(k+1), 1 )
  end do

  do k = n, 1, -1
    b(k) = b(k) / a(k,k)
    t = - b(k)
    call saxpy ( k-1, t, a(1,k), 1, b(1), 1 )
  end do

  write ( *, '(a,g14.6)' ) ' '
  write ( *, '(a,g14.6)' ) '  First five entries of solution:'
  write ( *, '(a,g14.6)' ) ' '
  write ( *, '(2x,5g14.6)' ) b(1:5)

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests SASUM.
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

  real ( kind = 4 ) a(lda,na)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 4 ) sasum
  real ( kind = 4 ) x(nx)

  do i = 1, nx
    x(i) = ( -1.0E+00 )**i * real ( 2 * i, kind = 4 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  SASUM adds the absolute values of elements '
  write ( *, '(a)' ) '  of a single precision real vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, nx
    write ( *, '(2x,i6,g14.6)' ) i, x(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  SASUM ( NX, X, 1 ) =   ', sasum ( nx, x, 1 )
  write ( *, '(a,g14.6)' ) '  SASUM ( NX/2, X, 2 ) = ', sasum ( nx/2, x, 2 )
  write ( *, '(a,g14.6)' ) '  SASUM ( 2, X, NX/2 ) = ', sasum ( 2, x, nx/2 )

  a(1:lda,1:na) = 0.0E+00

  do i = 1, ma
    do j = 1, na
      a(i,j) = ( -1.0E+00 )**( i + j ) * real ( 10 * i + j, kind = 4 )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate with a matrix A:'
  write ( *, '(a)' ) ' '
  do i = 1, ma
    write ( *, '(2x,5g14.6)' ) a(i,1:na)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  SASUM(MA,A(1,2),1) =   ', &
    sasum ( ma, a(1,2), 1 )
  write ( *, '(a,g14.6)' ) '  SASUM(NA,A(2,1),LDA) = ', &
    sasum ( na, a(2,1), lda )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests SAXPY.
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

  real ( kind = 4 ) da
  integer ( kind = 4 ) i
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(n)

  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do

  do i = 1, n
    y(i) = real ( 100 * i, kind = 4 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  SAXPY adds a multiple of '
  write ( *, '(a)' ) '  one single precision real vector to another.'
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

  da = 1.0E+00
  call saxpy ( n, da, x, 1, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  SAXPY ( N, ', da, ', X, 1, Y, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  do i = 1, n
    y(i) = real ( 100 * i, kind = 4 )
  end do

  da = -2.0E+00
  call saxpy ( n, da, x, 1, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  SAXPY ( N, ', da, ', X, 1, Y, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  do i = 1, n
    y(i) = real ( 100 * i, kind = 4 )
  end do

  da = +3.0E+00
  call saxpy ( 3, da, x, 2, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  SAXPY ( 3, ', da, ', X, 2, Y, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  do i = 1, n
    y(i) = real ( 100 * i, kind = 4 )
  end do

  da = -4.0E+00
  call saxpy ( 3, da, x, 1, y, 2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  SAXPY ( 3, ', da, ', X, 1, Y, 2 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests SCOPY.
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

  real ( kind = 4 ) a(5,5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 4 ) x(10)
  real ( kind = 4 ) y(10)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  SCOPY copies a single precision real vector.'

  do i = 1, 10
    x(i) = real ( i, kind = 4 )
  end do

  do i = 1, 10
    y(i) = real ( 10 * i, kind = 4 )
  end do

  do i = 1, 5
    do j = 1, 5
      a(i,j) = real ( 10 * i + j, kind = 4 )
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

  call scopy ( 5, x, 1, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  SCOPY ( 5, X, 1, Y, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  do i = 1, 10
    y(i) = real ( 10 * i, kind = 4 )
  end do

  call scopy ( 3, x, 2, y, 3 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  SCOPY ( 3, X, 2, Y, 3 )'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i6,g14.6)' ) i, y(i)
  end do

  call scopy ( 5, x, 1, a, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  SCOPY ( 5, X, 1, A, 1 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A = '
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,5f8.2)' ) a(i,1:5)
  end do

  do i = 1, 5
    do j = 1, 5
      a(i,j) = real ( 10 * i + j, kind = 4 )
    end do
  end do

  call scopy ( 5, x, 2, a, 5 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  SCOPY ( 5, X, 2, A, 5 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A = '
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,5f8.2)' ) a(i,1:5)
  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests SDOT.
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

  real ( kind = 4 ) a(lda,lda)
  real ( kind = 4 ) b(ldb,ldb)
  real ( kind = 4 ) c(ldc,ldc)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 4 ) sdot
  real ( kind = 4 ) sum1
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  SDOT computes the dot product of '
  write ( *, '(a)' ) '  single precision real vectors.'

  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do

  do i = 1, n
    y(i) = - real ( i, kind = 4 )
  end do

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( i + j, kind = 4 )
    end do
  end do

  do i = 1, n
    do j = 1, n
      b(i,j) = real ( i - j, kind = 4 )
    end do
  end do
!
!  To compute a simple dot product of two vectors, use a
!  call like this:
!
  sum1 = sdot ( n, x, 1, y, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Dot product of X and Y is ', sum1
!
!  To multiply a ROW of a matrix A times a vector X, we need to
!  specify the increment between successive entries of the row of A:
!
  sum1 = sdot ( n, a(2,1), lda, x, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Product of row 2 of A and X is ', sum1
!
!  Product of a column of A and a vector is simpler:
!
  sum1 = sdot ( n, a(1,2), 1, x, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Product of column 2 of A and X is ', sum1
!
!  Here's how matrix multiplication, c = a*b, could be done
!  with SDOT:
!
  do i = 1, n
    do j = 1, n
      c(i,j) = sdot ( n, a(i,1), lda, b(1,j), 1 )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix product computed with SDOT:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,5g14.6)' ) c(i,1:n)
  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests SMACH.
!
!  Discussion:
!
!    The SMACH routine is not part of the official BLAS release.
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
  implicit none

  integer ( kind = 4 ) job
  real ( kind = 4 ) smach

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  SMACH computes several machine-dependent'
  write ( *, '(a)' ) '  single precision real arithmetic parameters.'

  write ( *, '(a)' ) ' '
  write ( *, * ) '  SMACH(1)  = machine epsilon = ', smach ( 1 )
  write ( *, * ) '  SMACH(2)  = a tiny value    = ', smach ( 2 )
  write ( *, * ) '  SMACH(3)  = a huge value    = ', smach ( 3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FORTRAN90 parameters:'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  EPSILON() = machine epsilon = ', epsilon ( 1.0E+00 )
  write ( *, * ) '  TINY()    = a tiny value    = ', tiny ( 1.0E+00 )
  write ( *, * ) '  HUGE()    = a huge value    = ', huge ( 1.0E+00 )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests SNRM2.
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
  real ( kind = 4 ) a(lda,lda)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) j
  real ( kind = 4 ) snrm2
  real ( kind = 4 ) sum1
  real ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  SNRM2 computes the Euclidean norm of '
  write ( *, '(a)' ) '  a single precision real vector.'
!
!  Compute the euclidean norm of a vector:
!
  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The vector X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,f8.4)' ) i, x(i)
  end do
  incx = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The 2-norm of X is ', snrm2 ( n, x, incx )
!
!  Compute the euclidean norm of a row or column of a matrix:
!
  do i = 1, n
    do j = 1, n
      a(i,j) = real ( i + j, kind = 4 )
    end do
  end do

  incx = lda
  sum1 = snrm2 ( n, a(2,1), incx )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The 2-norm of row 2 of A is ', sum1

  incx = 1
  sum1 = snrm2 ( n, a(1,2), incx )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The 2-norm of column 2 of A is ', sum1

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests SROT.
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

  real ( kind = 4 ) c
  integer ( kind = 4 ) i
  real ( kind = 4 ) s
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(n)

  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do

  do i = 1, n
    y(i) = real ( i * i - 12, kind = 4 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) &
    '  SROT carries out a single precision real Givens rotation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X and Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
  end do

  c = 0.5E+00
  s = sqrt ( 1.0E+00 - c * c )
  call srot ( n, x, 1, y, 1, c, s )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a,f8.4,a)' ) '  SROT ( N, X, 1, Y, 1, ', c, ',', s, ' )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
  end do

  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do

  do i = 1, n
    y(i) = real ( i * i - 12, kind = 4 )
  end do

  c = x(1) / sqrt ( x(1) * x(1) + y(1) * y(1) )
  s = y(1) / sqrt ( x(1) * x(1) + y(1) * y(1) )
  call srot ( n, x, 1, y, 1, c, s )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a,f8.4,a)' ) '  SROT ( N, X, 1, Y, 1, ', c, ',', s, ' )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests SROTG.
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

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  real ( kind = 4 ) r
  real ( kind = 4 ) r4_uniform_01
  real ( kind = 4 ) s
  real ( kind = 4 ) sa
  real ( kind = 4 ) sb
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5
  real ( kind = 4 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  SROTG generates a single precision real Givens rotation'
  write ( *, '(a)' ) '    (  C  S ) * ( A ) = ( R )'
  write ( *, '(a)' ) '    ( -S  C )   ( B )   ( 0 )'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    a = r4_uniform_01 ( seed )
    b = r4_uniform_01 ( seed )

    sa = a
    sb = b

    call srotg ( sa, sb, c, s )

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
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests SSCAL.
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

  real ( kind = 4 ) da
  integer ( kind = 4 ) i
  real ( kind = 4 ) x(n)

  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  SSCAL multiplies a single precision real scalar times'
  write ( *, '(a)' ) '  a single precision real vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X = '
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, x(i)
  end do

  da = 5.0E+00
  call sscal ( n, da, x, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  SSCAL ( N, ', da, ', X, 1 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, x(i)
  end do

  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do

  da = -2.0E+00
  call sscal ( 3, da, x, 2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  SSCAL ( 3, ', da, ', X, 2 )'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6)' ) i, x(i)
  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests SSWAP.
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

  integer ( kind = 4 ) i
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(n)

  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do

  do i = 1, n
    y(i) = real ( 100 * i, kind = 4 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  SSWAP swaps two single precision real vectors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X and Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
  end do

  call sswap ( n, x, 1, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SSWAP ( N, X, 1, Y, 1 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X and Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
  end do

  do i = 1, n
    x(i) = real ( i, kind = 4 )
  end do

  do i = 1, n
    y(i) = real ( 100 * i, kind = 4 )
  end do

  call sswap ( 3, x, 2, y, 1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  SSWAP ( 3, X, 2, Y, 1 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X and Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
  end do

  return
end
function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r4_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
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
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r4_uniform_01

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r4_uniform_01 = real ( real ( seed, kind = 8 ) * 4.656612875D-10, kind = 4 )

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
