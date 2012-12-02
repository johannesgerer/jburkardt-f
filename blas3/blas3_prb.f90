program main

!*****************************************************************************80
!
!! MAIN is the main program for BLAS3_PRB.
!
!  Discussion:
!
!    BLAS3_PRB tests the BLAS3 routines.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLAS3_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BLAS3 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLAS3_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests a Cholesky solver.
!
!  Discussion:
!
!    Solve a positive definite symmetric linear system A*X = B.
!
  implicit none

  integer, parameter :: n = 20

  real ( kind = 8 ) a(n,n)
  integer i
  integer info
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Demonstrate Cholesky factor/solve routines'
  write ( *, '(a)' ) '  built out of level 1, 2 and 3 BLAS.'
  write ( *, '(a)' ) ' '
!
!  Set the entries of the matrix.
!
  do i = 1, n

    a(i,1:n) = 0.0D+00

    a(i,i) = 2.0D+00

    if ( 1 < i ) then
      a(i,i-1) = -1.0D+00
    end if

    if ( i < n ) then
      a(i,i+1) = -1.0D+00
    end if

  end do
!
!  Set the right hand side vector.
!
  x(1:n-1) = 0.0D+00
  x(n) = real ( n + 1, kind = 8 )
!
!  Factor the matrix.
!
  call dlltb ( n, a, n, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The matrix is singular!'
    return
  end if
!
!  Solve the system.
!
  call dllts ( n, a, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 10 entries of solution:'
  write ( *, '(a)' ) ' '
  write ( *, '(5g14.6)' ) x(1:10)

  return
end
subroutine dlltb ( n, a, lda, info )

!*****************************************************************************80
!
!! DLLTB Cholesky factors a positive definite symmetric matrix.
!
!  Discussion:
!
!    The factorization has the form A = L * L'.
!
!    The parameter NB determines the 'blocking factor' used
!    by the level 3 BLAS routines.
!
!    DLLT is an equivalent routine, but it solves the problem
!    at a lower level.
!
!    DLLTB solves the problem in chunks,
!    which is hoped to allow for greater optimization.
!
  implicit none

  integer lda
  integer n
  integer, parameter :: nb = 64

  real ( kind = 8 ) a(lda,*)
  integer info
  integer j
  integer jb

  info = 0

  do j = 1, n, nb

    jb = min ( nb, n-j+1 )
!
!  Update diagonal block
!
    call dsyrk ( 'lower', 'no transpose', jb, j-1, -1.0D+00, a(j,1), lda, &
      1.0D+00, a(j,j), lda )
!
!  Factorize diagonal block and test for non-positive-definiteness.
!
    call dllt ( jb, a(j,j), lda, info )

    if ( info /= 0 )then
      info = info + j - 1
      return
    end if

    if ( j + jb <= n ) then
!
!  Update subdiagonal block.
!
      call dgemm ( 'no transpose', 'transpose', n-j-jb+1, jb, j-1, &
        -1.0D+00, a(j+jb,1), lda, a(j,1), lda, 1.0D+00, a(j+jb,j), lda )
!
!  Compute the subdiagonal block of L.
!
      call dtrsm ( 'right', 'lower', 'transpose', 'non-unit', n-j-jb+1, &
        jb, 1.0D+00, a(j,j), lda, a(j+jb,j), lda )

    end if

  end do

  return
end
subroutine dllt ( n, a, lda, info )

!*****************************************************************************80
!
!! DLLT computes an L*L' factorization of an SPD matrix A.
!
!  Discussion:
!
!    DLLT uses level 2 and level 1 BLAS routines.
!
  implicit none

  integer lda
  integer n

  real ( kind = 8 ) a(lda,n)
  integer info
  integer j
  real ( kind = 8 ) ddot

  info = 0

  do j = 1, n
!
!  Update A(J,J).
!
    a(j,j) = a(j,j) - ddot ( j-1, a(j,1), lda, a(j,1), lda )
!
!  Compute L(J,J) and test for non-positive-definiteness.
!
    if ( a(j,j) <= 0.0D+00 ) then
      info = j
      return
    end if

    a(j,j) = sqrt ( a(j,j) )
!
!  Update elements J+1 to N of J-th column.
!
    if ( j  <  n ) then
      call dgemv ( 'no transpose', n-j, j-1, -1.0D+00, a(j+1,1), lda, &
        a(j,1), lda, 1.0D+00, a(j+1,j), 1 )
!
!  Compute elements J+1 to N of J-th column of L.
!
      call dscal ( n-j, 1.0D+00/a(j,j), a(j+1,j), 1 )

    end if

  end do

  return
end
subroutine dllts ( n, a, lda, b )

!*****************************************************************************80
!
!! DLLTS solves A*X=B, for a positive definite symmetric matrix.
!
!  Discussion:
!
!    The matrix A should have been factored by DLLT or DLLTB.
!
  implicit none

  integer lda
  integer n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer k
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t

  do k = 1, n
    t = ddot ( k-1, a(k,1), lda, b(1), 1 )
    b(k) = ( b(k) - t ) / a(k,k)
  end do

  do k = n, 1, -1
    b(k) = b(k) / a(k,k)
    t = - b(k)
    call daxpy ( k-1, t, a(k,1), lda, b(1), 1 )
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
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
