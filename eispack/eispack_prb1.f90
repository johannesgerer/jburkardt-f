program main

!*****************************************************************************80
!
!! MAIN is the main program for EISPACK_PRB1.
!
!  Discussion:
!
!    EISPACK_PRB1 calls the EISPACK sample programs.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EISPACK_PRB1'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the EISPACK library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test065 ( )
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
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EISPACK_PRB1'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CG.
!
!  Discussion:
!
!    CG is for the eigenvalues of a complex general matrix.
!
!    eigenvalues and eigenvectors of a complex general matrix
!    note that the eigenvalues of such a matrix are in general complex.
!    however, we will use the same example we used before, namely
!    a hermitian matrix, so the eigenvalues will in fact be real.
!
!    (3     1     0     0+2i)
!    (1     3     0-2i  0   )
!    (0     0+2i  1     1   )
!    (0-2i  0     1     1   )
!
!    The eigenvalues are 2+2*sqrt(2), 2-2*sqrt(2), 4 and 0
!
!    The eigenvector matrix is
!
!    (  1+sqrt(2),  1,                -1,          1)
!    (  1+sqrt(2),  1,                 1,         -1)
!    (     i,       -(1+sqrt(2))*i,    i,          i)
!    (    -i,        (1+sqrt(2))*i,    i,          i)
!
!    Note that the actual eigenvector matrix from EISPACK could
!    be scaled by a real value, or by i, and the columns may
!    appear in any order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) ai(n,n)
  real ( kind = 8 ) ar(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) matz
  real ( kind = 8 ) wi(n)
  real ( kind = 8 ) wr(n)
  real ( kind = 8 ) xi(n,n)
  real ( kind = 8 ) xr(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  CG computes the eigenvalues and eigenvectors of '
  write ( *, '(a)' ) '  a complex general matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  Set the values of the matrix.
!
  ar(1,1) = 3.0D+00
  ar(1,2) = 1.0D+00
  ar(1,3) = 0.0D+00
  ar(1,4) = 0.0D+00

  ar(2,1) = 1.0D+00
  ar(2,2) = 3.0D+00
  ar(2,3) = 0.0D+00
  ar(2,4) = 0.0D+00

  ar(3,1) = 0.0D+00
  ar(3,2) = 0.0D+00
  ar(3,3) = 1.0D+00
  ar(3,4) = 1.0D+00

  ar(4,1) = 0.0D+00
  ar(4,2) = 0.0D+00
  ar(4,3) = 1.0D+00
  ar(4,4) = 1.0D+00

  ai(1,1) = 0.0D+00
  ai(1,2) = 0.0D+00
  ai(1,3) = 0.0D+00
  ai(1,4) = 2.0D+00

  ai(2,1) = 0.0D+00
  ai(2,2) = 0.0D+00
  ai(2,3) = -2.0D+00
  ai(2,4) = 0.0D+00

  ai(3,1) = 0.0D+00
  ai(3,2) = 2.0D+00
  ai(3,3) = 0.0D+00
  ai(3,4) = 0.0D+00

  ai(4,1) = -2.0D+00
  ai(4,2) = -0.0D+00
  ai(4,3) = -0.0D+00
  ai(4,4) = 0.0D+00
!
!  matz = 0 for eigenvalues only,
!  matz = 1 for eigenvalues and eigenvectors.
!
  matz = 1
  call cg ( n, ar, ai, wr, wi, matz, xr, xi, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec2_print ( n, wr, wi, '  Real and imaginary parts of eigenvalues:' )

  if ( matz /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'The eigenvectors are:'
    do i = 1, n
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Eigenvector ', i
      write ( *, '(a)' ) ' '
      do j = 1, n
        write ( *, '(2g14.6)' ) xr(i,j), xi(i,j)
      end do
    end do
  end if

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests CH.
!
!  Discussion:
!
!    CH is for the eigenvalues of a complex hermitian matrix.
!
!    Eigenvalues and eigenvectors of a complex hermitian matrix
!
!    Note that the eigenvalues (though not the eigenvectors) of
!    a hermitian matrix are real.
!
!    (3     1     0     0+2i)
!    (1     3     0-2i  0   )
!    (0     0+2i  1     1   )
!    (0-2i  0     1     1   )
!
!    The eigenvalues are 2+2*sqrt(2), 2-2*sqrt(2), 4 and 0
!
!    The eigenvector matrix is
!
!    (  1+sqrt(2),  1,                -1,          1)
!    (  1+sqrt(2),  1,                 1,         -1)
!    (     i,       -(1+sqrt(2))*i,    i,          i)
!    (    -i,        (1+sqrt(2))*i,    i,          i)
!
!    Note that the actual eigenvector matrix from EISPACK could
!    be scaled by a real value, or by i, and the columns may
!    appear in any order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) ar(n,n)
  real ( kind = 8 ) ai(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) matz
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) xr(n,n)
  real ( kind = 8 ) xi(n,n)
!
!  Set the values of the matrix.
!
  ar(1,1) = 3.0D+00
  ar(1,2) = 1.0D+00
  ar(1,3) = 0.0D+00
  ar(1,4) = 0.0D+00

  ar(2,1) = 1.0D+00
  ar(2,2) = 3.0D+00
  ar(2,3) = 0.0D+00
  ar(2,4) = 0.0D+00

  ar(3,1) = 0.0D+00
  ar(3,2) = 0.0D+00
  ar(3,3) = 1.0D+00
  ar(3,4) = 1.0D+00

  ar(4,1) = 0.0D+00
  ar(4,2) = 0.0D+00
  ar(4,3) = 1.0D+00
  ar(4,4) = 1.0D+00

  ai(1,1) = 0.0D+00
  ai(1,2) = 0.0D+00
  ai(1,3) = 0.0D+00
  ai(1,4) = 2.0D+00

  ai(2,1) = 0.0D+00
  ai(2,2) = 0.0D+00
  ai(2,3) = -2.0D+00
  ai(2,4) = 0.0D+00

  ai(3,1) = 0.0D+00
  ai(3,2) = 2.0D+00
  ai(3,3) = 0.0D+00
  ai(3,4) = 0.0D+00

  ai(4,1) = -2.0D+00
  ai(4,2) = -0.0D+00
  ai(4,3) = -0.0D+00
  ai(4,4) = 0.0D+00

  matz = 1

  call ch ( n, ar, ai, w, matz, xr, xi, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  CH computes the eigenvalues and eigenvectors of'
  write ( *, '(a)' ) '  a complex hermitian matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Error flag = ', ierr
  write ( *, '(a)' ) ' '

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  write ( *, '(a)' ) ' '

  if ( matz /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Eigenvectors are:'
    do i = 1, n
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Eigenvector ', i
      write ( *, '(a)' ) ' '
      do j = 1, n
        write ( *, '(2g14.6)' ) xr(i,j), xi(i,j)
      end do
    end do
  end if

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests MINFIT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: nb = 1
  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) acopy(m,n)
  real ( kind = 8 ) b(m,nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(m)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  a(1,1) =   1.00D+00
  a(2,1) =   2.05D+00
  a(3,1) =   3.06D+00
  a(4,1) = - 1.02D+00
  a(5,1) =   4.08D+00

  a(1,2) =   1.00D+00
  a(2,2) = - 1.00D+00
  a(3,2) =   1.00D+00
  a(4,2) =   2.00D+00
  a(5,2) = - 1.00D+00

  acopy(1:m,1:n) = a(1:m,1:n)

  b(1,1) = 1.98D+00
  b(2,1) = 0.95D+00
  b(3,1) = 3.98D+00
  b(4,1) = 0.92D+00
  b(5,1) = 2.90D+00

  r(1:m) = - b(1:m,1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  MINFIT solves an overdetermined linear system'
  write ( *, '(a)' ) '  using least squares methods.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows = ', m
  write ( *, '(a,i8)' ) '  Matrix columns = ', n

  call r8mat_print ( m, n, a, '  The matrix A:' )

  call r8mat_print ( m, nb, b, '  The right hand side B:' )

  call minfit ( m, m, n, a, w, nb, b, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03 - Warning!'
    write ( *, '(a,i8)' ) '  MINFIT error code IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The singular values:' )
!
!  B now contains U' * B.
!  We need to divide by the singular values, and multiply by V.
!
  b(1:n,1) = b(1:n,1) / w(1:n)

  do i = 1, n
    x(i) = 0.0D+00
    do j = 1, n
      x(i) = x(i) + a(i,j) * b(j,1)
    end do
  end do

  call r8vec_print ( n, x, '  The least squares solution X:' )

  do i = 1, m
    do j = 1, n
      r(i) = r(i) + acopy(i,j) * x(j)
    end do
  end do

  call r8vec_print ( m, r, '  The residual A * X - B:' )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests RG.
!
!  Discussion:
!
!    RG is for the eigenvalues of a general real matrix.
!
!    The matrix A is nonsymmetric.  The eigenvalues may therefore be
!    complex numbers.
!
!    ( 33  16  72)
!    (-24 -10 -57)
!    ( -8  -4 -17)
!
!    The eigenvalues of A are (1,2,3)
!
!    The eigenvectors of A are
!
!    (-15 -16 -4)
!    ( 12  13  3)
!    (  4   4  1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) acopy(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) matz
  real ( kind = 8 ) sum3
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) wi(n)
  real ( kind = 8 ) wr(n)
  real ( kind = 8 ) x(n,n)
!
!  Set the values of the matrix.
!
  a(1,1) = 33.0D+00
  a(1,2) = 16.0D+00
  a(1,3) = 72.0D+00

  a(2,1) = -24.0D+00
  a(2,2) = -10.0D+00
  a(2,3) = -57.0D+00

  a(3,1) = -8.0D+00
  a(3,2) = -4.0D+00
  a(3,3) = -17.0D+00
!
!  RG overwrites A with garbage, so save a copy now!
!
  acopy(1:n,1:n) = a(1:n,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  RG computes the eigenvalues and eigenvectors of'
  write ( *, '(a)' ) '  a real general matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, a, '  The matrix A:' )

  matz = 1

  call rg ( n, a, wr, wi, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag was IERR = ', ierr
    return
  end if

  call r8vec2_print ( n, wr, wi, '  Real and imaginary parts of eigenvalues:' )

  if ( matz /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The eigenvectors may be complex:'
    do j = 1, n
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Eigenvector ', j
      write ( *, '(a)' ) ' '
      do i = 1, n
        if ( wi(j) == 0.0D+00 ) then
          write ( *, '(g14.6)' ) x(i,j)
        else if ( 0.0D+00 < wi(j) ) then
          write ( *, '(2g14.6)' ) x(i,j), x(i,j+1)
        else if ( wi(j) < 0.0D+00 ) then
          write ( *, '(2g14.6)' ) x(i,j-1), -x(i,j)
        end if
      end do
    end do
!
!  Check.
!  First, restore the original values of A.
!
    a(1:n,1:n) = acopy(1:n,1:n)

    do k = 1, n

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Residuals (A*x-Lambda*x) for eigenvalue ', k
      write ( *, '(a)' ) ' '

      if ( wi(k) == 0.0D+00 ) then

        do i = 1, n
          sum3 = 0.0D+00
          do j = 1, n
            sum3 = sum3 + a(i,j) * x(j,k)
          end do
          sum3 = sum3 - wr(k) * x(i,k)
          write ( *, '(g14.6)' ) sum3
        end do

      else if ( 0.0D+00 < wi(k) ) then

        do i = 1, n
          sum1 = 0.0D+00
          sum2 = 0.0D+00
          do j = 1, n
            sum1 = sum1 + a(i,j) * x(j,k)
            sum2 = sum2 + a(i,j) * x(j,k+1)
          end do
          sum1 = sum1 - wr(k) * x(i,k) + wi(k) * x(i,k+1)
          sum2 = sum2 - wi(k) * x(i,k) - wr(k) * x(i,k+1)
          write ( *, '(2g14.6)' ) sum1, sum2
        end do

      else if ( wi(k) < 0.0D+00 ) then

        do i = 1, n
          sum1 = 0.0D+00
          sum2 = 0.0D+00
          do j = 1, n
            sum1 = sum1 + a(i,j) * x(j,k-1)
            sum2 = sum2 - a(i,j) * x(j,k)
          end do
          sum1 = sum1 - wr(k) * x(i,k-1) - wi(k) * x(i,k)
          sum2 = sum2 - wi(k) * x(i,k-1) + wr(k) * x(i,k)
          write ( *, '(2g14.6)' ) sum1, sum2
        end do

      end if

    end do

  end if

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests RGG.
!
!  Discussion:
!
!    RGG is for a real generalized general eigenvalue problem.
!
!    A generalized eigenvalue problem.  Given matrices A and B, find
!    N numbers LAMBDA, and for each LAMBDA a vector X, so that
!
!      A*x = lambda*B*x
!
!    The matrix A is
!
!    ( -7 7  6  6)
!    (-10 8 10  8)
!    ( -8 3 10 11)
!    ( -4 0  4 12)
!
!    The matrix B is
!
!    (2 1 0 0)
!    (1 2 1 0)
!    (0 1 2 1)
!    (0 0 1 2)
!
!    The correct eigenvalues LAMBDA are
!
!    (1,2,3,4)
!
!    The correct eigenvectors X are
!
!    (4 3 2 1)
!    (3 3 2 1)
!    (2 2 2 1)
!    (1 1 1 1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) acopy(n,n)
  real ( kind = 8 ) alfi(n)
  real ( kind = 8 ) alfr(n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) bcopy(n,n)
  real ( kind = 8 ) beta(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) matz
  real ( kind = 8 ) sum3
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) x(n,n)
!
!  Set the values in the A matrix.
!
  a(1,1) = -7.0D+00
  a(1,2) = 7.0D+00
  a(1,3) = 6.0D+00
  a(1,4) = 6.0D+00

  a(2,1) = -10.0D+00
  a(2,2) = 8.0D+00
  a(2,3) = 10.0D+00
  a(2,4) = 8.0D+00

  a(3,1) = -8.0D+00
  a(3,2) = 3.0D+00
  a(3,3) = 10.0D+00
  a(3,4) = 11.0D+00

  a(4,1) = -4.0D+00
  a(4,2) = 0.0D+00
  a(4,3) = 4.0D+00
  a(4,4) = 12.0D+00
!
!  Save a copy of A.
!
  acopy(1:n,1:n) = a(1:n,1:n)
!
!  Set the values in the B matrix.
!
  b(1,1) = 2.0D+00
  b(1,2) = 1.0D+00
  b(1,3) = 0.0D+00
  b(1,4) = 0.0D+00

  b(2,1) = 1.0D+00
  b(2,2) = 2.0D+00
  b(2,3) = 1.0D+00
  b(2,4) = 0.0D+00

  b(3,1) = 0.0D+00
  b(3,2) = 1.0D+00
  b(3,3) = 2.0D+00
  b(3,4) = 1.0D+00

  b(4,1) = 0.0D+00
  b(4,2) = 0.0D+00
  b(4,3) = 1.0D+00
  b(4,4) = 2.0D+00
!
!  Save a copy of B.
!
  bcopy(1:n,1:n) = b(1:n,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  RGG for real generalized problem.'
  write ( *, '(a)' ) '  Find scalars LAMBDA and vectors X so that'
  write ( *, '(a)' ) '    A*X = LAMBDA * B * X'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call r8mat_print ( n, n, b, '  The matrix B:' )

  matz = 1

  call rgg ( n, a, b, alfr, alfi, beta, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  alfr(1:n) = alfr(1:n) / beta(1:n)
  alfi(1:n) = alfi(1:n) / beta(1:n)

  call r8vec2_print ( n, alfr, alfi, '  Real and imaginary parts of eigenvalues:' )

  if ( matz /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The eigenvectors are:'
    do i = 1, n
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Eigenvector ', i
      write ( *, '(a)' ) ' '
      do j = 1, n
        write ( *, '(g14.6)' ) x(i,j)
      end do
    end do
  end if
!
!  Check.
!  First, restore the original values of A and B.
!
  if ( matz /= 0 ) then

    a(1:n,1:n) = acopy(1:n,1:n)
    b(1:n,1:n) = bcopy(1:n,1:n)

    do k = 1, n
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) &
        '  Residuals (A*x-(Alfr+Alfi*I)*B*x) for eigenvalue ', k
      write ( *, '(a)' ) ' '

      if ( alfi(k) == 0.0D+00 ) then

        do i = 1, n

          sum3 = dot_product ( a(i,1:n), x(1:n,k) )

          do j = 1, n
            sum3 = sum3 - alfr(k) * b(i,j) * x(j,k)
          end do

          write ( *, '(g14.6)' ) sum3
        end do

      else if ( 0.0D+00 < alfi(k) ) then

        do i = 1, n

          sum1 = 0.0D+00
          sum2 = 0.0D+00
          do j = 1, n
            sum1 = sum1 + a(i,j) * x(j,k)
            sum2 = sum2 + a(i,j) * x(j,k+1)
          end do

          do j = 1, n
            sum1 = sum1 - alfr(k) * b(i,j) * x(j,k) &
                        + alfi(k) * b(i,j) * x(j,k+1)

            sum2 = sum2 - alfi(k) * b(i,j) * x(j,k) &
                        - alfr(k) * b(i,j) * x(j,k+1)

          end do

          write ( *, '(2g14.6)' ) sum1, sum2
        end do

      else if ( alfi(k) < 0.0D+00 ) then

        do i = 1, n

          sum1 = 0.0D+00
          sum2 = 0.0D+00
          do j = 1, n
            sum1 = sum1 + a(i,j) * x(j,k-1)
            sum2 = sum2 - a(i,j) * x(j,k)
          end do

          do j = 1, n

            sum1 = sum1 - alfr(k) * b(i,j) * x(j,k-1) &
                        - alfi(k) * b(i,j) * x(j,k)

            sum2 = sum2 - alfi(k) * b(i,j) * x(j,k-1) &
                        + alfr(k) * b(i,j) * x(j,k)
          end do

          write ( *, '(2g14.6)' ) sum1, sum2
        end do

      end if

    end do

  end if

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests RS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a2(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) matz
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n,n)
!
!  Set the values in the matrix.
!
  a(1,1) = 5.0D+00
  a(1,2) = 4.0D+00
  a(1,3) = 1.0D+00
  a(1,4) = 1.0D+00

  a(2,1) = 4.0D+00
  a(2,2) = 5.0D+00
  a(2,3) = 1.0D+00
  a(2,4) = 1.0D+00

  a(3,1) = 1.0D+00
  a(3,2) = 1.0D+00
  a(3,3) = 4.0D+00
  a(3,4) = 2.0D+00

  a(4,1) = 1.0D+00
  a(4,2) = 1.0D+00
  a(4,3) = 2.0D+00
  a(4,4) = 4.0D+00
!
!  Save a copy of the matrix.
!
  a2(1:n,1:n) = a(1:n,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  RS computes the eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, a, '  The matrix A:' )

  matz = 1

  call rs ( n, a, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix:' )

    r(1:n,1:n) = matmul ( a2(1:n,1:n), x(1:n,1:n) )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual (A-Lambda*I)*X:' )

  end if

  return
end
subroutine test065 ( )

!*****************************************************************************80
!
!! TEST065 tests RS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a2(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) matz
  real ( kind = 8 ) r(n,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST065'
  write ( *, '(a)' ) '  RS computes the eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  seed = 123456789

  call r8mat_uniform_01 ( n, n, seed, a )

  do i = 1, n - 1
    do j = i + 1, n
      t = ( a(i,j) + a(j,i) ) / 2.0D+00
      a(i,j) = t
      a(j,i) = t
    end do
  end do
!
!  Save a copy of the matrix.
!
  a2(1:n,1:n) = a(1:n,1:n)

  call r8mat_print ( n, n, a, '  The matrix A:' )

  matz = 1

  call rs ( n, a, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST065 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix:' )

    r(1:n,1:n) = matmul ( a2(1:n,1:n), x(1:n,1:n) )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual (A-Lambda*I)*X:' )

  end if

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests RSB.
!
!  Discussion:
!
!    RSB solves the eigenvalue problem for a symmetric banded matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: mb = 2

  real ( kind = 8 ) a(n,mb)
  real ( kind = 8 ) a2(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) matz
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n,n)
!
!  A contains the band matrix in banded storage.
!
  a(1:n,1:mb) = 0.0D+00
  a(1:n,mb) = 2.0D+00
  a(2:n,1) = -1.0D+00
!
!  A2 contains the band matrix in full storage.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a2(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a2(i,j) = - 1.0D+00
      else
        a2(i,j) = 0.0D+00
      end if
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  RSB computes the eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric band matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, a2, '  The matrix A:' )

  matz = 1

  call rsb ( n, mb, a, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST07 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual (A-Lambda*I)*X:' )

  end if

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests RSG.
!
!  Discussion:
!
!    RGG is for a real generalized eigenvalue problem of the form
!
!      A*x = lambda*B*x
!
!    with A symmetric and B positive definite symmetric.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a2(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) b2(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) matz
  real ( kind = 8 ) sum3
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n,n)

  do i = 1, n
    do j = 1, n
      a(i,j) = abs ( i - j )
    end do
  end do

  a2(1:n,1:n) = a(1:n,1:n)

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        b(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        b(i,j) = - 1.0D+00
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  b2(1:n,1:n) = b(1:n,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  RSG for real symmetric generalized problem.'
  write ( *, '(a)' ) '  Find scalars LAMBDA and vectors X so that'
  write ( *, '(a)' ) '    A*X = LAMBDA * B * X'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call r8mat_print ( n, n, b, '  The matrix B:' )

  matz = 1

  call rsg ( n, a, b, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    a(1:n,1:n) = a2(1:n,1:n)
    b(1:n,1:n) = b2(1:n,1:n)

    do k = 1, n

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) &
        'Residuals (A*x-(w*I)*B*x) for eigenvalue ', k
      write ( *, '(a)' ) ' '

        do i = 1, n

          sum3 = 0.0D+00
          do j = 1, n
            sum3 = sum3 + a(i,j) * x(j,k)
          end do

          do j = 1, n
            sum3 = sum3 - w(k) * b(i,j) * x(j,k)
          end do

          write ( *, '(g14.6)' ) sum3
        end do

    end do

  end if

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests RSGAB.
!
!  Discussion:
!
!    RGGAB is for a real generalized eigenvalue problem of the form
!
!      A*B*x = lambda*x
!
!    with A symmetric and B positive definite symmetric.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a2(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) b2(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) matz
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n,n)

  do i = 1, n
    do j = 1, n
      a(i,j) = abs ( i - j )
    end do
  end do

  a2(1:n,1:n) = a(1:n,1:n)

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        b(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        b(i,j) = - 1.0D+00
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  b2(1:n,1:n) = b(1:n,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09:'
  write ( *, '(a)' ) '  RSGAB for real symmetric generalized problem.'
  write ( *, '(a)' ) '  Find scalars LAMBDA and vectors X so that'
  write ( *, '(a)' ) '    A*B*X = LAMBDA * X'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call r8mat_print ( n, n, b, '  The matrix B:' )

  matz = 1

  call rsgab ( n, a, b, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST09 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  EThe eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( b2, x )

    r = matmul ( a2, r )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual matrix (A*B-Lambda*I)*X:' )

  end if

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests RSGBA.
!
!  Discussion:
!
!    RGGBA is for a real generalized eigenvalue problem of the form
!
!      B*A*x = lambda*x
!
!    with A symmetric and B positive definite symmetric.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a2(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) b2(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) matz
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n,n)

  do i = 1, n
    do j = 1, n
      a(i,j) = abs ( i - j )
    end do
  end do

  a2(1:n,1:n) = a(1:n,1:n)

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        b(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        b(i,j) = - 1.0D+00
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  b2(1:n,1:n) = b(1:n,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10:'
  write ( *, '(a)' ) '  RSGBA for real symmetric generalized problem.'
  write ( *, '(a)' ) '  Find scalars LAMBDA and vectors X so that'
  write ( *, '(a)' ) '    B*A*X = LAMBDA * X'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call r8mat_print ( n, n, b, '  The matrix B:' )

  matz = 1

  call rsgba ( n, a, b, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST10 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    r = matmul ( b2, r )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual matrix (B*A-Lambda*I)*X:' )

  end if

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests RSM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: m = n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a2(n,n)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) matz
  real ( kind = 8 ) r(n,m)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n,m)

  a(1,1) = 5.0D+00
  a(1,2) = 4.0D+00
  a(1,3) = 1.0D+00
  a(1,4) = 1.0D+00

  a(2,1) = 4.0D+00
  a(2,2) = 5.0D+00
  a(2,3) = 1.0D+00
  a(2,4) = 1.0D+00

  a(3,1) = 1.0D+00
  a(3,2) = 1.0D+00
  a(3,3) = 4.0D+00
  a(3,4) = 2.0D+00

  a(4,1) = 1.0D+00
  a(4,2) = 1.0D+00
  a(4,3) = 2.0D+00
  a(4,4) = 4.0D+00

  a2(1:n,1:n) = a(1:n,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  RSM computes some eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n
  write ( *, '(a,i8)' ) '  Number of eigenvectors desired = ', m

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call rsm ( n, a, w, m, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( 0 < m ) then

    call r8mat_print ( n, m, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, m
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, m, r, '  The residual (A-Lambda*I)*X:' )

  end if

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests RSP.
!
!  Discussion:
!
!    RSP is for the eigenvalues of a real symmetric packed matrix.
!
!    A is symmetric.  Because of this, we know that the eigenvalues
!    of A must be real (rather than complex) numbers.
!
!
!    The entries of A are
!
!    (5 4 1 1)
!    (4 5 1 1)
!    (1 1 4 2)
!    (1 1 2 4)
!
!    The eigenvalues of A are (10, 5, 2, 1)
!
!    One set of eigenvectors of A is:
!
!    ( 2 -1  0 -1)
!    ( 2 -1  0  1)
!    ( 1  2 -1  0)
!    ( 1  2  1  0)
!
!    However, this set is not orthonormal, and EISPACK will compute
!    a different set of values.
!
!    Note that the I-th eigenvector corresponding to the I-th eigenvalue
!    consists of the I-th column of the above matrix of eigenvectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: nv = ( n * ( n + 1 ) ) / 2

  real ( kind = 8 ) a(nv)
  real ( kind = 8 ) a2(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) matz
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n,n)
!
!  Set the values in the matrix.
!
  a(1) = 5.0D+00

  a(2) = 4.0D+00
  a(3) = 5.0D+00

  a(4) = 1.0D+00
  a(5) = 1.0D+00
  a(6) = 4.0D+00

  a(7) = 1.0D+00
  a(8) = 1.0D+00
  a(9) = 2.0D+00
  a(10) = 4.0D+00

  a2(1,1) = 5.0D+00
  a2(1,2) = 4.0D+00
  a2(1,3) = 1.0D+00
  a2(1,4) = 1.0D+00

  a2(2,1) = 4.0D+00
  a2(2,2) = 5.0D+00
  a2(2,3) = 1.0D+00
  a2(2,4) = 1.0D+00

  a2(3,1) = 1.0D+00
  a2(3,2) = 1.0D+00
  a2(3,3) = 4.0D+00
  a2(3,4) = 2.0D+00

  a2(4,1) = 1.0D+00
  a2(4,2) = 1.0D+00
  a2(4,3) = 2.0D+00
  a2(4,4) = 4.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  RSP computes the eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric packed matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, a2, '  The matrix A:' )

  matz = 1

  call rsp ( n, nv, a, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST12 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag was IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual matrix (A-Lambda*I)*X:' )

  end if

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests RSPP.
!
!  Discussion:
!
!    RSPP is for some eigenvalues of a real symmetric packed matrix.
!
!    A is symmetric.  Because of this, we know that the eigenvalues
!    of A must be real (rather than complex) numbers.
!
!
!    The entries of A are
!
!    (5 4 1 1)
!    (4 5 1 1)
!    (1 1 4 2)
!    (1 1 2 4)
!
!    The eigenvalues of A are (10, 5, 2, 1)
!
!    One set of eigenvectors of A is:
!
!    ( 2 -1  0 -1)
!    ( 2 -1  0  1)
!    ( 1  2 -1  0)
!    ( 1  2  1  0)
!
!    However, this set is not orthonormal, and EISPACK will compute
!    a different set of values.
!
!    Note that the I-th eigenvector corresponding to the I-th eigenvalue
!    consists of the I-th column of the above matrix of eigenvectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: m = n
  integer ( kind = 4 ), parameter :: nv = ( n * ( n + 1 ) ) / 2

  real ( kind = 8 ) a(nv)
  real ( kind = 8 ) a2(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) matz
  real ( kind = 8 ) r(n,m)
  logical              type
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(n,m)
!
!  Set the values in the matrix.
!
  a(1) = 5.0D+00

  a(2) = 4.0D+00
  a(3) = 5.0D+00

  a(4) = 1.0D+00
  a(5) = 1.0D+00
  a(6) = 4.0D+00

  a(7) = 1.0D+00
  a(8) = 1.0D+00
  a(9) = 2.0D+00
  a(10) = 4.0D+00

  a2(1,1) = 5.0D+00
  a2(1,2) = 4.0D+00
  a2(1,3) = 1.0D+00
  a2(1,4) = 1.0D+00

  a2(2,1) = 4.0D+00
  a2(2,2) = 5.0D+00
  a2(2,3) = 1.0D+00
  a2(2,4) = 1.0D+00

  a2(3,1) = 1.0D+00
  a2(3,2) = 1.0D+00
  a2(3,3) = 4.0D+00
  a2(3,4) = 2.0D+00

  a2(4,1) = 1.0D+00
  a2(4,2) = 1.0D+00
  a2(4,3) = 2.0D+00
  a2(4,4) = 4.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  RSPP finds some eigenvalues and eigenvectors of'
  write ( *, '(a)' ) '  a real symmetric packed matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, a2, '  The matrix A:' )
!
!  Set MATZ = 0 for no eigenvectors, 1 for eigenvectors.
!
  matz = 1
!
!  TYPE = TRUE to find smallest eigenvalues, FALSE for largest.
!
  type = .true.

  call rspp ( n, nv, a, w, matz, x, ierr, m, type )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST13 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag was IERR = ', ierr
    return
  end if

  call r8vec_print ( m, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call r8mat_print ( n, m, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, m
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, m, r, '  The residual matrix (A-Lambda*I)*X:' )

  end if

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests RST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) matz
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n,n)
!
!  Here is where the matrix is defined.
!
  w(1:n) = 2.0D+00

  e(1) = 0.0D+00
  e(2:n) = -1.0D+00
!
!  We only set up and store the matrix A this way in order to make it easy
!  to compute the residual.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a(i,j) = - 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  RST computes the eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric tridiagonal matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, a, '  The matrix A:' )

  matz = 1

  call rst ( n, w, e, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST14 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( a, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual matrix (A-Lambda*I)*X:' )

  end if

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests RT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,3)
  real ( kind = 8 ) a2(n,n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) matz
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n,n)
!
!  Here is where the matrix is defined.
!
  a(2:n,  1) = - 1.0D+00
  a(1:n,  2) =   2.0D+00
  a(1:n-1,3) = - 1.0D+00
!
!  We only set up and store the matrix A this way in order to make it easy
!  to compute the residual.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a2(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a2(i,j) = - 1.0D+00
      else
        a2(i,j) = 0.0D+00
      end if
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  RT computes the eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real sign-symmetric tridiagonal matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, a2, '  The matrix A:' )

  matz = 1

  call rt ( n, a, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual matrix (A-Lambda*I)*X:' )

  end if

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests SVD.
!
!  Discussion:
!
!    In our special example, the matrix is square and symmetric.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: m = n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical              matu
  logical              matv
  real ( kind = 8 ) r(m,n)
  real ( kind = 8 ) u(m,n)
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) w(n)
!
!  Set the values of the matrix.
!
  a(1,1) = 0.9900D+00
  a(1,2) = 0.0020D+00
  a(1,3) = 0.0060D+00
  a(1,4) = 0.0020D+00

  a(2,1) = 0.0020D+00
  a(2,2) = 0.9900D+00
  a(2,3) = 0.0020D+00
  a(2,4) = 0.0060D+00

  a(3,1) = 0.0060D+00
  a(3,2) = 0.0020D+00
  a(3,3) = 0.9900D+00
  a(3,4) = 0.0020D+00

  a(4,1) = 0.0020D+00
  a(4,2) = 0.0060D+00
  a(4,3) = 0.0020D+00
  a(4,4) = 0.9900D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  SVD computes the singular value decomposition'
  write ( *, '(a)' ) '  of a real general matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( m, n, a, '  The matrix A:' )

  matu = .true.
  matv = .true.

  call svd ( m, n, a, w, matu, u, matv, v, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST16 - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The singular values S' )

  call r8mat_print ( m, n, u, '  The U matrix:' )

  call r8mat_print ( n, n, v, '  The V matrix:' )

  do j = 1, n
    v(1:n,j) = w(j) * v(1:n,j)
  end do

  r(1:m,1:n) = matmul ( u(1:m,1:n), transpose ( v(1:n,1:n) ) )

  call r8mat_print ( m, n, r, '  The product U * S * Transpose(V):' )

  return
end
