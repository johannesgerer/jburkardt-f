program band_qr_prb

!*****************************************************************************80
!
!! BAND_QR_PRB tests BAND_QR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
! Modified:
!
!    02 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BAND_QR:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BAND_QR library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BAND_QR:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DGEQRF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nrhs = 1

  integer ( kind = 4 ), parameter :: lda = m
  integer ( kind = 4 ), parameter :: ldb = max ( m, n )
  integer ( kind = 4 ), parameter :: lwork = n
  integer ( kind = 4 ), parameter :: tau_size = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(ldb,nrhs)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  real ( kind = 8 ) tau(tau_size)
  real ( kind = 8 ) work(lwork)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Call the standard LAPACK routine'
  write ( *, '(a)' ) '  DGEQRF to get the QR factorization of'
  write ( *, '(a)' ) '  a matrix A stored in GE format, but which'
  write ( *, '(a)' ) '  is actually banded.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Then solve the linear system A*X=B using an'
  write ( *, '(a)' ) '  explicit version of DGEQRS_TWO.'
  write ( *, '(a)' ) '  DGEQRS_TWO is a version of '
  write ( *, '(a)' ) '  DGEQRS with no calls to subroutines.'

  do j = 1, n
    do i = 1, m
      if ( j == i - 1 ) then
        a(i,j) = - 1.0D+00
      else if ( j == i ) then
        a(i,j) = 2.0D+00
      else if ( j == i + 1 ) then
        a(i,j) = -1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  call r8ge_print ( m, n, a, '  Input matrix:' )
!
!  QR factor the matrix.
!
  call dgeqrf ( m, n, a, lda, tau, work, lwork, info )

  if ( info == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGEQRF called successfully.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGEQRF returned error flag.'
    write ( *, '(a,i8)' ) '  INFO = ', info
    return
  end if

  call r8ge_print ( m, n, a, '  Factored matrix:' )

  call r8vec_print ( tau_size, tau, '  Tau:' );
!
!  Set up and solve a linear system using DQEQRS_TWO.
!
  b(1:n-1,1) = 0.0D+00
  b(n,1) = real ( n + 1, kind = 8 )

  call dgeqrs_two ( m, n, nrhs, a, lda, tau, b, ldb, work, lwork, info )

  if ( info == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGEQRS_TWO called successfully.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGEQRS_TWO returned error flag.'
    write ( *, '(a,i8)' ) '  INFO = ', info
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution from DGEQRS_TWO:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(2x,g14.6)' ) b(i,1)
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests DGEBQR2.
!
!  Discussion:
!
!    The same operation is being carried out as in Test #1, except
!    that DGEBQR2 knows that the matrix A is banded.  The matrix A
!    is still stored as a full storage "GE" matrix, though.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nrhs = 1

  integer ( kind = 4 ), parameter :: lda = m
  integer ( kind = 4 ), parameter :: ldb = max ( m, n )
  integer ( kind = 4 ), parameter :: lwork = min ( ml + mu, n )
  integer ( kind = 4 ), parameter :: tau_size = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(ldb,nrhs)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  real ( kind = 8 ) tau(tau_size)
  real ( kind = 8 ) work(lwork)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Call the LAPACK-like routine'
  write ( *, '(a)' ) '  DGEBQR2 to get the QR factorization of'
  write ( *, '(a)' ) '  a matrix A stored in GE format, but which'
  write ( *, '(a)' ) '  is actually banded.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Then solve the linear system A*X=B using an'
  write ( *, '(a)' ) '  explicit version of DGEQRS_TWO.'
  write ( *, '(a)' ) '  DGEQRS_TWO is a version of '
  write ( *, '(a)' ) '  DGEQRS with no calls to subroutines.'

  do j = 1, n
    do i = 1, m
      if ( j == i - 1 ) then
        a(i,j) = - 1.0D+00
      else if ( j == i ) then
        a(i,j) = 2.0D+00
      else if ( j == i + 1 ) then
        a(i,j) = -1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  call r8ge_print ( m, n, a, '  Input matrix:' )
!
!  QR factor the matrix.
!
  call dgebqr2 ( m, n, ml, mu, a, lda, tau, work, info )

  if ( info == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGEBQR2 called successfully.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGEBQR2 returned error flag.'
    write ( *, '(a,i8)' ) '  INFO = ', info
    return
  end if

  call r8ge_print ( m, n, a, '  Factored matrix:' )

  call r8vec_print ( tau_size, tau, '  Tau:' );
!
!  Set up and solve a linear system using DQEQRS_TWO.
!
  b(1:n-1,1) = 0.0D+00
  b(n,1) = real ( n + 1, kind = 8 )

  call dgeqrs_two ( m, n, nrhs, a, lda, tau, b, ldb, work, lwork, info )

  if ( info == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGEQRS_TWO called successfully.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGEQRS_TWO returned error flag.'
    write ( *, '(a,i8)' ) '  INFO = ', info
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution from DGEQRS_TWO:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(2x,g14.6)' ) b(i,1)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests DGBBQR2.
!
!  Discussion:
!
!    The same calculation is carried out as in Test #2.  However,
!    DGBBQR2 knows that the matrix A is banded, and that it is stored
!    in the "GB" format.  The QR factorization data overwrites the
!    values of A, and fits in the GB format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nrhs = 1

  integer ( kind = 4 ), parameter :: lda = 2 * ml + mu + 1
  integer ( kind = 4 ), parameter :: ldb = max ( m, n )
  integer ( kind = 4 ), parameter :: lwork = min ( ml + mu, n )
  integer ( kind = 4 ), parameter :: tau_size = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(ldb,nrhs)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  real ( kind = 8 ) tau(tau_size)
  real ( kind = 8 ) work(lwork)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Call the LAPACK-like routine'
  write ( *, '(a)' ) '  DGBBQR2 to get the QR factorization of'
  write ( *, '(a)' ) '  a matrix A stored in GB format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call DGBBQRS to solve the linear system'
  write ( *, '(a)' ) '  A*X=B.'
!
!  Zero out the matrix.
!
  a(1:lda,1:n) = 0.0D+00
!
!  Superdiagonal,
!  Diagonal,
!  Subdiagonal.
!
  do j = 2, n
    a(ml + mu + 1 - 1,j) = -1.0D+00
  end do

  do j = 1, n
    a(ml + mu + 1,j) = 2.0D+00
  end do

  do j = 1, n - 1
    a(ml + mu + 1 + 1,j) = -1.0D+00
  end do

  call r8gb_print ( m, n, ml, mu, a, '  Input matrix:' )
!
!  QR factor the matrix.
!
  call dgbbqr2 ( m, n, ml, mu, a, lda, tau, work, info )

  if ( info == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGBBQR2 called successfully.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGBBQR2 returned error flag.'
    write ( *, '(a,i8)' ) '  INFO = ', info
    return
  end if

  call r8gb_print ( m, n, ml, mu, a, '  Factored matrix:' )

  call r8vec_print ( tau_size, tau, '  Tau:' );
!
!  Set up and solve a linear system using DQBBQRS.
!
  b(1:n-1,1) = 0.0D+00
  b(n,1) = real ( n + 1, kind = 8 )

  call dgbbqrs ( m, n, ml, mu, nrhs, a, lda, tau, b, ldb, work, lwork, info )

  if ( info == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGBBQRS called successfully.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGBBQRS returned error flag.'
    write ( *, '(a,i8)' ) '  INFO = ', info
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution from DGBBQRS:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(2x,g14.6)' ) b(i,1)
  end do

  return
end

