program main

!*****************************************************************************80
!
!! MAIN is the main program for EXPONENT_FORMAT_OVERFLOW.
!
!  Discussion:
!
!    EXPONENT_FORMAT_OVERFLOW examines the format used to print real values
!    in exponential format, in cases where the exponent has large magnitude.
!
!    It has been observed that, at least for some compilers, the case in
!    which the exponent has three digits is handled in a very bad way
!    that is misleading and liable to result in errors, particularly if
!    one program writes out the data and another program is to read it
!    back in.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXPONENT_FORMAT_OVERFLOW:'
  write ( *, '(a)' ) '  FORTRAN90 version.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXPONENT_FORMAT_OVERFLOW:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 prints some large values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) x6
  real ( kind = 8 ) x7
  real ( kind = 8 ) x8
  real ( kind = 8 ) x9
  real ( kind = 8 ) x10
  real ( kind = 8 ) x11
  real ( kind = 8 ) x12

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Real numbers can have exponents greater than 99.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  It is unsettling to discover that there is a'
  write ( *, '(a)' ) '  common FLAW in certain output formats, which means'
  write ( *, '(a)' ) '  that the printing of real numbers with exponents'
  write ( *, '(a)' ) '  of magnitude more than 99 is handled poorly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The "E" marker, used to indicate scientific notation,'
  write ( *, '(a)' ) '  is simply suppressed, as though to make room for'
  write ( *, '(a)' ) '  the leading digit of the exponent.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  But this inconsistency can be deadly.  In particular,'
  write ( *, '(a)' ) '  if you write out such data, it CANNOT be read back in'
  write ( *, '(a)' ) '  properly!  (THAT example is easy to write, but I will'
  write ( *, '(a)' ) '  do that later.)'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Define some big numbers:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X1 = huge ( X1 )'
  write ( *, '(a)' ) '  X2 = 0.99D+00 * X1'
  write ( *, '(a)' ) '  X3 = X1 / 100.0D+00'
  write ( *, '(a)' ) '  X4 = 1.0D+101.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Define some small numbers:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X5 = tiny ( X5 )'
  write ( *, '(a)' ) '  X6 = 1.01D+00 * X5'
  write ( *, '(a)' ) '  X7 = X5 * 100.0D+00'
  write ( *, '(a)' ) '  X8 = 1.0D-101.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Define some comparison numbers:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X9  = 1.0D+98'
  write ( *, '(a)' ) '  X10 = 1.0D-98'
  write ( *, '(a)' ) '  X11 = 123456789.0'
  write ( *, '(a)' ) '  X12 = 0.123456789'

  x1 = huge ( x1 )
  x2 = 0.99D+00 * x1
  x3 = x1 / 10.0D+00
  x4 = 1.0D+101

  x5 = tiny ( x5 )
  x6 = 1.01D+00 * x5
  x7 = 100.0D+00 * x5
  x8 = 1.0D-101

  x9 = 1.0D+98
  x10 = 1.0D-98
  x11 = 123456789.0D+00
  x12 = 0.123456789D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Print with a WRITE(*,*) format:'
  write ( *, '(a)' ) '  This seems to work OK:'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  X1 = ', x1
  write ( *, * ) '  X2 = ', x2
  write ( *, * ) '  X3 = ', x3
  write ( *, * ) '  X4 = ', x4
  write ( *, * ) '  X5 = ', x5
  write ( *, * ) '  X6 = ', x6
  write ( *, * ) '  X7 = ', x7
  write ( *, * ) '  X8 = ', x8
  write ( *, * ) '  X9 = ', x9
  write ( *, * ) '  X10 = ', x10
  write ( *, * ) '  X11 = ', x11
  write ( *, * ) '  X12 = ', x12

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Print with a WRITE(*,''(G24.10)'') format:'
  write ( *, '(a)' ) '  Notice the disaster:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.10)' ) '  X1 = ', x1
  write ( *, '(a,g24.10)' ) '  X2 = ', x2
  write ( *, '(a,g24.10)' ) '  X3 = ', x3
  write ( *, '(a,g24.10)' ) '  X4 = ', x4
  write ( *, '(a,g24.10)' ) '  X5 = ', x5
  write ( *, '(a,g24.10)' ) '  X6 = ', x6
  write ( *, '(a,g24.10)' ) '  X7 = ', x7
  write ( *, '(a,g24.10)' ) '  X8 = ', x8
  write ( *, '(a,g24.10)' ) '  X9 = ', x9
  write ( *, '(a,g24.10)' ) '  X10 = ', x10
  write ( *, '(a,g24.10)' ) '  X11 = ', x11
  write ( *, '(a,g24.10)' ) '  X12 = ', x12

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Print with a WRITE(*,''(D24.10)'') format:'
  write ( *, '(a)' ) '  Same disaster:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,d24.10)' ) '  X1 = ', x1
  write ( *, '(a,d24.10)' ) '  X2 = ', x2
  write ( *, '(a,d24.10)' ) '  X3 = ', x3
  write ( *, '(a,d24.10)' ) '  X4 = ', x4
  write ( *, '(a,d24.10)' ) '  X5 = ', x5
  write ( *, '(a,d24.10)' ) '  X6 = ', x6
  write ( *, '(a,d24.10)' ) '  X7 = ', x7
  write ( *, '(a,d24.10)' ) '  X8 = ', x8
  write ( *, '(a,d24.10)' ) '  X9 = ', x9
  write ( *, '(a,d24.10)' ) '  X10 = ', x10
  write ( *, '(a,d24.10)' ) '  X11 = ', x11
  write ( *, '(a,d24.10)' ) '  X12 = ', x12

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Print with a WRITE(*,''(E24.10)'') format:'
  write ( *, '(a)' ) '  Same disaster:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,e24.10)' ) '  X1 = ', x1
  write ( *, '(a,e24.10)' ) '  X2 = ', x2
  write ( *, '(a,e24.10)' ) '  X3 = ', x3
  write ( *, '(a,e24.10)' ) '  X4 = ', x4
  write ( *, '(a,e24.10)' ) '  X5 = ', x5
  write ( *, '(a,e24.10)' ) '  X6 = ', x6
  write ( *, '(a,e24.10)' ) '  X7 = ', x7
  write ( *, '(a,e24.10)' ) '  X8 = ', x8
  write ( *, '(a,e24.10)' ) '  X9 = ', x9
  write ( *, '(a,e24.10)' ) '  X10 = ', x10
  write ( *, '(a,e24.10)' ) '  X11 = ', x11
  write ( *, '(a,e24.10)' ) '  X12 = ', x12

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
