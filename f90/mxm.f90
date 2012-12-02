module global_data

!*****************************************************************************80
!
!! GLOBAL_DATA is a module which stores the data shared by the routines.
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
  implicit none

  integer, parameter :: n1 = 100
  integer, parameter :: n2 = 100
  integer, parameter :: n3 = 100

  real, dimension (n1,n2) :: a
  real, dimension (n2,n3) :: b
  real, dimension (n1,n3) :: c
  integer, dimension(6) :: count
  integer :: count_max
  integer :: count_rate
  real, parameter :: pi = 3.1415926535E+00

end module global_data
program main

!*****************************************************************************80
!
!! MAIN is the main program for MXM.
!
!  Discussion:
!
!    MXM carries out a timing of matrix multiplication.
!
!    Two matrices A and B are defined, and then multiplied to
!    get the product matrix C.  Various operations are timed
!    and reported.  It is possible that some operations might
!    be carried out in parallel.  The F90 intrinsic functions
!    TRANSPOSE and MATMUL are used to perform the calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2000
!
!  Author:
!
!    John Burkardt
!
  use global_data

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MXM'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A program to multiply two matrices.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A * B = C'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Dimensions:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Matrix    Rows    Columns    Entries'
  write ( *, '(a)' ) ' '
  write ( *, '(9x,a1,i8,3x,i8,i11)' ) 'A', n1, n2, n1*n2
  write ( *, '(9x,a1,i8,3x,i8,i11)' ) 'B', n2, n3, n2*n3
  write ( *, '(9x,a1,i8,3x,i8,i11)' ) 'C', n1, n3, n1*n3
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Floating point operations:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '    Adds:       ', n1*n2*n3
  write ( *, '(a,i12)' ) '    Multiplies: ', n1*n2*n3
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Begin execution.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Set A by calling A_SET.'

  call system_clock ( count(1), count_rate, count_max )

  call a_set ( )

  call system_clock ( count(2), count_rate, count_max )

  write ( *, '(a)' ) 'Set B = Transpose(A) by calling B_SET.'

  call system_clock ( count(3), count_rate, count_max )

  call b_set ( )

  call system_clock ( count(4), count_rate, count_max )

  write ( *, '(a)' ) 'Compute C by calling A_TIMES_B.'

  call system_clock ( count(5), count_rate, count_max )

  call a_times_b ( )

  call system_clock ( count(6), count_rate, count_max )

  write ( *, '(a)' ) 'Call REPORT.'

  call report ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MXM:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine a_set ( )

!*****************************************************************************80
!
!! A_SET sets the matrix A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2000
!
!  Author:
!
!    John Burkardt
!
  use global_data

  real angle
  integer i
  integer j
  integer n

  n = max ( n1, n2 )

  do i = 1, n1
    do j = 1, n2

      angle = 2.0E+00 * real ( i * j ) * pi / real ( 2 * n + 1 )
      a(i,j) = 2.0E+00 * sin ( angle ) / sqrt ( real ( 2 * n + 1 ) )

    end do
  end do

  return
end
subroutine b_set ( )

!*****************************************************************************80
!
!! B_SET sets the matrix B, by transposing A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2000
!
!  Author:
!
!    John Burkardt
!
  use global_data

  b = transpose ( a )

  return
end
subroutine a_times_b ( )

!*****************************************************************************80
!
!! A_TIMES_B carries out the multiplication, using MATMUL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2000
!
!  Author:
!
!    John Burkardt
!
  use global_data

  c = matmul ( a, b )

  return
end
subroutine report ( )

!*****************************************************************************80
!
!! REPORT prints a report of the results.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2000
!
!  Author:
!
!    John Burkardt
!
  use global_data

  integer :: i
  integer :: ihi
  integer :: j
  integer :: jhi
  integer :: temp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Task                     Ticks   Seconds'
  write ( *, '(a)' ) ' '
  temp = 1
  write ( *, '(a,i12,g14.6)' ) '  1 Tick            ', temp, &
    real ( temp ) / count_rate
  temp = count(2) - count(1)
  write ( *, '(a,i12,g14.6)' ) '  Set A             ', temp, &
    real ( temp ) / count_rate
  temp = count(4) - count(3)
  write ( *, '(a,i12,g14.6)' ) '  Set B = Trans(A): ', temp, &
    real ( temp ) / count_rate
  temp = count(6) - count(5)
  write ( *, '(a,i12,g14.6)' ) '  Set C = A * B:    ', temp, &
    real ( temp ) / count_rate

  ihi = min ( n1, 5 )
  jhi = min ( n3, 5 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Initial block of the matrix:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,5(i7,7x))' ) ( j, j = 1, jhi )
  write ( *, '(a)' ) ' '
  do i = 1, ihi
    write ( *, '(i4,5g14.6)' ) i, c(i,1:jhi)
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
