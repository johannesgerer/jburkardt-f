program main

!*****************************************************************************80
!
!! MAIN is the main program for MATRIX_FUNCTION_TEST.
!
!  Discussion:
!
!    MATRIX_FUNCTION_TEST tests the ability to define a matrix function.
!
!    This INTERFACE is the key to how a function can return a vector or
!    array value.  The INTERFACE has to be included in the source code of
!    any routine that uses the function.
!
!    If there are many interfaces, or many places where an interface is used,
!    it may be more convenient to put all the interfaces into a module, and
!    then invoke it where needed by a simple USE statement.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  interface
    function symmetric_part ( n, a )
      integer n
      real, dimension ( n, n ) :: a
      real, dimension ( n, n ) :: symmetric_part
    end function symmetric_part
  end interface

  integer, parameter :: n = 5

  real, dimension(n,n) :: a
  real ahi
  real alo
  real, dimension(n,n) :: b
  integer i

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATRIX_FUNCTION_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate the use of the matrix function facility.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The matrix is of order : ', n

  alo = - 3.0E+00
  ahi =   3.0E+00

  call random_number ( harvest = a(1:n,1:n) )

  a(1:n,1:n) = alo + ( ahi - alo ) * a(1:n,1:n)

  a(1:n,1:n) = nint ( a(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initial matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,5f10.4)' ) a(i,1:n)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Symmetric part B:'
  write ( *, '(a)' ) ' '

  b(1:n,1:n) = symmetric_part ( n, a )

  do i = 1, n
    write ( *, '(2x,5f10.4)' ) b(i,1:n)
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATRIX_FUNCTION_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function symmetric_part ( n, a )

!*****************************************************************************80
!
!! SYMMETRIC_PART returns the symmetric part of a matrix.
!
!  Discussion:
!
!    The symmetric part of a (square) matrix A is defined as the matrix B
!    made by averaging the original matrix and its transpose.  That is,
!
!      B(I,J) = 0.5 * ( A(I,J) + A(J,I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N,N), a matrix whose symmetric part is desired.
!
!    Output, real SYMMETRIC_PART(N,N), the symmetric part of the matrix.
!
  implicit none

  integer n

  real a(n,n)
  real symmetric_part(n,n)

  symmetric_part(1:n,1:n) = 0.5E+00 * (  a(1:n,1:n) + transpose ( a(1:n,1:n) ) )

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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
