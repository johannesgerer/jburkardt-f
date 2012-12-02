program main

!*****************************************************************************80
!
!! MAIN is the main program for REAL_KIND_PRECISION24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REAL_KIND_PRECISION24'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Inquire about the properties of real data that'
  write ( *, '(a)' ) '  is requested with a precision of 24 decimal digits'
  write ( *, '(a)' ) '  and a max exponent 30.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REAL_KIND_PRECISION24'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 examines the default real numbers.
!
!  Discussion:
!
!    We use SELECTED_REAL_KIND to request a real type with
!    a precision of 24 decimal digits and a maximum exponent of 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: precision24 = selected_real_kind ( 24, 30 )
  real ( kind = precision24 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use SELECTED_REAL_KIND to request real numbers'
  write ( *, '(a)' ) '  with a precision of 24 decimal digits.'

  if ( precision24 == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Warning!'
    write ( *, '(a)' ) '  There is no real kind with the requested'
    write ( *, '(a)' ) '  precision.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HUGE returns the largest value of the given type.'
  write ( *, '(a)' ) '  TINY returns the smallest value of the given type.'
  write ( *, '(a)' ) '  EPSILON returns the precision of a real type.'
  write ( *, '(a)' ) ' '
  write ( *, * ) '    HUGE(R) =    ', huge ( r )
  write ( *, * ) '    TINY(R) =    ', tiny ( r )
  write ( *, '(a,g14.6)' ) '    EPSILON(R) = ', epsilon ( r )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DIGITS counts the significant binary digits.'
  write ( *, '(a)' ) '  RANGE provides the decimal exponent range.'
  write ( *, '(a)' ) '  PRECISION provides the decimal precision.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '    DIGITS(R) =    ', digits ( r )
  write ( *, '(a,i12)' ) '    RANGE(R) =     ', range ( r )
  write ( *, '(a,i12)' ) '    PRECISION(R) = ', precision ( r )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RADIX provides the base of the model.'
  write ( *, '(a)' ) '  MAXEXPONENT returns the maximum exponent of a variable.'
  write ( *, '(a)' ) '  MINEXPONENT returns the minimum exponent of a variable.'
  write ( *, '(a)' ) '  (These are exponents of the RADIX)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '    RADIX(R) =       ', radix ( r )
  write ( *, '(a,i12)' ) '    MAXEXPONENT(R) = ', maxexponent ( r )
  write ( *, '(a,i12)' ) '    MINEXPONENT(R) = ', minexponent ( r )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KIND returns the "kind" of a variable.'
  write ( *, '(a,i12)' ) '    KIND(R) = ', kind ( r )

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
