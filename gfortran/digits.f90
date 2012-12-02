program main

!*****************************************************************************80
!
!! MAIN is the main program for DIGITS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DIGITS'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  How many digits are appropriate for data?'

  call test01
  call test02

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DIGITS:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 looks at data stored as REAL ( KIND = 4 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: digit_max = 20

  real    ( kind = 4 ) c
  real    ( kind = 4 ) d
  integer ( kind = 4 ) digit
  real    ( kind = 4 ), dimension ( digit_max ) :: pi = (/ &
    3.0E+00, &
    3.1E+00, &
    3.14E+00, &
    3.141E+00, &
    3.1415E+00, &
    3.14159E+00, &
    3.141592E+00, &
    3.1415926E+00, &
    3.14159265E+00, &
    3.141592653E+00, &
    3.1415926535E+00, &
    3.14159265358E+00, &
    3.141592653589E+00, &
    3.1415926535897E+00, &
    3.14159265358979E+00, &
    3.141592653589793E+00, &
    3.1415926535897932E+00, &
    3.14159265358979323E+00, &
    3.141592653589793238E+00, &
    3.1415926535897932384E+00 /)
  real    ( kind = 4 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  REAL ( KIND = 4 ) arithmetic'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Store PI using 1 through 20 digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Digits     sin(PI(I))   1+cos(PI(I))     (PI(I)-PI(20))'
  write ( *, '(a)' ) ' '

  do digit = 1, digit_max
    s = sin ( pi(digit) )
    c = 1.0E+00 + cos ( pi(digit) )
    d = pi(digit) - pi(digit_max)
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) digit, s, c, d
  end do

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 looks at data stored as REAL ( KIND = 8 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: digit_max = 20

  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  integer ( kind = 4 ) digit
  real    ( kind = 8 ), dimension ( digit_max ) :: pi = (/ &
    3.0D+00, &
    3.1D+00, &
    3.14D+00, &
    3.141D+00, &
    3.1415D+00, &
    3.14159D+00, &
    3.141592D+00, &
    3.1415926D+00, &
    3.14159265D+00, &
    3.141592653D+00, &
    3.1415926535D+00, &
    3.14159265358D+00, &
    3.141592653589D+00, &
    3.1415926535897D+00, &
    3.14159265358979D+00, &
    3.141592653589793D+00, &
    3.1415926535897932D+00, &
    3.14159265358979323D+00, &
    3.141592653589793238D+00, &
    3.1415926535897932384D+00 /)
  real    ( kind = 8 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  REAL ( KIND = 8 ) arithmetic'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Store PI using 1 through 20 digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Digits     sin(PI(I))   1+cos(PI(I))     (PI(I)-PI(20))'
  write ( *, '(a)' ) ' '

  do digit = 1, digit_max
    s = sin ( pi(digit) )
    c = 1.0D+00 + cos ( pi(digit) )
    d = pi(digit) - pi(digit_max)
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) digit, s, c, d
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
