program main

!*****************************************************************************80
!
!! MAIN is the main program for WHERE_TEST.
!
!  Discussion:
!
!    WHERE_TEST tests the FORTRAN90 WHERE statement.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 10
  integer, parameter :: n = 10

  real, dimension(m,n) :: a
  real ahi
  real alo

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WHERE_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate the use of the WHERE statement.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,2x,i8)' ) '  The matrix is of dimension : ', m, n
  write ( *, '(a)' ) '  In the printout, no more than 5 rows and columns '
  write ( *, '(a)' ) '  will be shown.'
 
  alo = -3.0E+00
  ahi =  3.0E+00

  call r4mat_random ( m, m, n, a, alo, ahi )

  a = nint ( a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initial matrix A:'

  call r4mat_print_some ( m, n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Where A nonzero, A = 2 / A:'

  where ( a /= 0.0E+00 ) 
    a = 2.0E+00 / a
  end where

  call r4mat_print_some ( m, n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Where A positive, A = Sqrt ( A ),'
  write ( *, '(a)' ) '  Elsewhere         A = - A**2.'

  where ( 0.0E+00 < a )
    a = sqrt ( a )
  elsewhere
    a = -(a**2)
  end where

  call r4mat_print_some ( m, n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Where |A| <= 1, A = arcsin ( A ):'

  where ( abs ( a ) <= 1.0E+00 ) 
    a = asin ( a )
  end where

  call r4mat_print_some ( m, n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WHERE_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r4mat_print_some ( m, n, a )

!*****************************************************************************80
!
!! R4MAT_PRINT prints some of a real matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the dimensions of the matrix.
!
!    Input, real A(M,N), the matrix to be printed.
!
  implicit none

  integer, intent ( in ) :: m
  integer, intent ( in ) :: n

  integer :: i
  integer :: ihi
  integer :: j
  integer :: jhi
  real, dimension ( m, n ), intent ( in ) :: a

  ihi = min ( m, 5 )
  jhi = min ( n, 5 )

  write ( *, '(a)' ) ' '
  write ( *, '(5x,5(i7,7x))' ) ( j, j = 1, jhi )
  write ( *, '(a)' ) ' '
  do i = 1, ihi
    write ( *, '(i3,2x,5g14.6)' ) i, a(i,1:jhi)
  end do

  return
end
subroutine r4mat_random ( lda, m, n, a, alo, ahi )

!*****************************************************************************80
!
!! R4MAT_RANDOM returns a matrix of uniform random values between AHI and ALO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer M, N, the number of rows and columns of A.
!
!    Output, real A(LDA,N), the random matrix.
!
!    Input, real ALO, AHI, the minimum and maximum values that
!    the matrix entries can have.
!
  implicit none

  integer lda
  integer n

  real a(lda,n)
  real ahi
  real alo
  integer m

  call random_number ( harvest = a(1:m,1:n) )

  a(1:m,1:n) = alo + ( ahi - alo ) * a(1:m,1:n)

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
