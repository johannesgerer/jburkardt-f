program main

!*****************************************************************************80
!
!! MAIN is the main program for ALLOCATABLE.
!
!  Discussion:
!
!    ALLOCATABLE attempts to allocate an allocatable array in a
!    subroutine, and use it in the main program.
!
!    Allocatable arrays were introduced with FORTRAN90, but this particular
!    usage was not allowed until AFTER the FORTRAN95 standard was approved.
!    ISO/IEC TR:15581.
!
!    I have tried compiling this with various versions of Gnu Fortran,
!    but have encountered a runtime error every time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: a(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ALLOCATABLE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Try to send an allocatable (but unallocated)'
  write ( *, '(a)' ) '  array to a subroutine, where it is allocated'
  write ( *, '(a)' ) '  assigned and returned.'

  call test01 ( n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Array size is N = ', n
  do i = 1, n
    write ( *, '(a,i8,2x,i8)' ) i, a(i)
  end do

  deallocate ( a )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ALLOCATABLE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( n, a )

!*****************************************************************************80
!
!! TEST01 allocates and assigns the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: a(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real    ( kind = 4 ) x

  call random_number ( harvest = x )

  n = int ( 10.0E+00 * x ) + 1

  allocate ( a(1:n) )

  do i = 1, n
    a(i) = 101 * i
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
