program main

!*****************************************************************************80
!
!! MAIN is the main program for SGE_MOD_PRB.
!
!  Discussion:
!
!    SGE_MOD_PRB tests the SGE_MODULE example.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2007
!
!  Author:
!
!    John Burkardt
!
  use sge_module

  implicit none

  integer, parameter :: n = 5

  real, dimension ( n, n ) :: a
  real, dimension ( n ) :: b
  real det
  integer i
  real sge_det
  real, dimension ( n ) :: x

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGE_MOD_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A calling program for a demonstration of the'
  write ( *, '(a)' ) '  use of modules.'
!
!  Assign values to A.
!
  call random_number ( harvest = a(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,5f12.4)' ) a(i,1:n)
  end do
!
!  Set up the desired solution vector.
!
  do i = 1, n
    x(i) = real ( i )
  end do
!
!  Compute the right hand side B = A*X.
!
  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

  call r4vec_print ( n, b, '  The right hand side vector:' )
!
!  Now "register" the matrix A.
!
  call sge_create ( n, a )
!
!  Factor it.
!
  call sge_fa ( )
!
!  Compute its determinant.
!
  det = sge_det ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The matrix determinant is ', det
!
!  Solve the linear system.
!
  call sge_sl ( b, x )

  call r4vec_print ( n, x, '  The solution vector:' )
!
!  Discard it.
!
  call sge_delete ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGE_MOD_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
