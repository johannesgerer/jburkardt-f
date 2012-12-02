program main

!*****************************************************************************80
!
!! MAIN is the main program for MXM_SERIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ) angle
  real ( kind = 8 ), allocatable :: b(:,:)
  real ( kind = 8 ), allocatable :: c(:,:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) s

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MXM_SERIAL:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Compute matrix product C = A * B.'

  n = 2000
  write ( *, '(a,i8)' ) '  The matrix order N                 = ', n
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
!
!  Loop 1: Evaluate A.
!
  s = 1.0D+00 / sqrt ( real ( n, kind = 8 ) )

  do i = 1, n
    do j = 1, n
      angle = 2.0D+00 * pi * ( i - 1 ) * ( j - 1 ) / real ( n, kind = 8 )
      a(i,j) = s * ( sin ( angle ) + cos ( angle ) ) 
    end do
  end do
!
!  Loop 2: Copy A into B.
!
  do i = 1, n
    do j = 1, n
      b(i,j) = a(i,j)
    end do
  end do
!
!  Loop 3: Compute C = A * B.
!
  do i = 1, n
    do j = 1, n
      c(i,j) = 0.0D+00
      do k = 1, n
        c(i,j) = c(i,j) + a(i,k) * b(k,j)
      end do
    end do
  end do

  write ( *, '(a,g14.6)' ) '  C(100,100)  = ', c(100,100)
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MXM_SERIAL:'
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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2001
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
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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
