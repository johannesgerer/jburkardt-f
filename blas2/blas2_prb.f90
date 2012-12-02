program main

!*****************************************************************************80
!
!! MAIN is the main program for BLAS2_PRB.
!
!  Discussion:
!
!    BLAS2_PRB tests the BLAS2 routines.
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
  write ( *, '(a)' ) 'BLAS2_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BLAS2 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLAS2_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DGEMV.
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 5

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer i
  integer incx
  integer incy
  integer j
  integer lda
  character trans
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For a general matrix A,'
  write ( *, '(a)' ) '  DGEMV computes y := alpha * A * x + beta * y'

  trans = 'N'
  alpha = 2.0D+00
  lda = m
  incx = 1
  beta = 3.0D+00
  incy = 1

  do i = 1, m
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( i == j - 1 .or. i == j + 1 ) then
        a(i,j) = -1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, m
    y(i) = real ( 10 * i, kind = 8 )
  end do

  call dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Result vector Y = '
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,g14.6)' ) y(i)
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests DGBMV.
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: kl = 1
  integer, parameter :: ku = 1
  integer, parameter :: lda = kl + 1 + ku

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer i
  integer incx
  integer incy
  integer j
  integer jhi
  integer jlo
  character trans
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For a general band matrix A,'
  write ( *, '(a)' ) '  DGBMV computes y := alpha * A * x + beta * y'

  trans = 'N'
  alpha = 2.0D+00
  incx = 1
  beta = 3.0D+00
  incy = 1

  do i = 1, m

    jlo = max ( 1, i - kl )
    jhi = min ( n, i + ku )

    do j = jlo, jhi

      if ( i == j ) then
        a(ku+1+i-j,j) = 2.0D+00
      else if ( i == j - 1 .or. i == j + 1 ) then
        a(ku+1+i-j,j) = -1.0D+00
      else
        a(ku+1+i-j,j) = 0.0D+00
      end if

    end do
  end do

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, m
    y(i) = real ( 10 * i, kind = 8 )
  end do

  call dgbmv ( trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Result vector Y = '
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,g14.6)' ) y(i)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests DSYMV.
!
  implicit none

  integer, parameter :: n = 5
  integer, parameter :: lda = n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer i
  integer incx
  integer incy
  integer j
  character uplo
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For a general symmetric matrix A,'
  write ( *, '(a)' ) '  DSYMV computes y := alpha * A * x + beta * y'

  uplo = 'U'
  alpha = 2.0D+00
  incx = 1
  beta = 3.0D+00
  incy = 1

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( i == j - 1 ) then
        a(i,j) = -1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    y(i) = real ( 10 * i, kind = 8 )
  end do

  call dsymv ( uplo, n, alpha, a, lda, x, incx, beta, y, incy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Result vector Y = '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,g14.6)' ) y(i)
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests DSBMV.
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: k = 1
  integer, parameter :: lda = k + 1

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer i
  integer incx
  integer incy
  integer j
  integer jhi
  integer jlo
  character uplo
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For a symmetric band matrix A,'
  write ( *, '(a)' ) '  DSBMV computes y := alpha * A * x + beta * y'

  uplo = 'U'
  alpha = 2.0D+00
  incx = 1
  beta = 3.0D+00
  incy = 1

  do i = 1, m

    jhi = min ( n, i + k )

    do j = i, jhi

      if ( i == j ) then
        a(k+1+i-j,j) = 2.0D+00
      else if ( i == j - 1 ) then
        a(k+1+i-j,j) = -1.0D+00
      else
        a(k+1+i-j,j) = 0.0D+00
      end if

    end do
  end do

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  do i = 1, m
    y(i) = real ( 10 * i, kind = 8 )
  end do

  call dsbmv ( uplo, n, k, alpha, a, lda, x, incx, beta, y, incy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Result vector Y = '
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,g14.6)' ) y(i)
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
