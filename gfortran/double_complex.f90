program main

!*****************************************************************************80
!
!! MAIN is the main program for DOUBLE_COMPLEX.
!
!  Discussion:
!
!    DOUBLE_COMPLEX runs the double precision complex demonstration examples.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DOUBLE_COMPLEX'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate double precision complex arithmetic.'

  call test01
  call test02
  call test03

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DOUBLE_COMPLEX'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 shows what you get by default.
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

  complex a
  complex b
  complex c

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use complex values of the default type.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RANGE provides the decimal exponent range.'
  write ( *, '(a,i12)' ) '    RANGE(I1) = ', range ( a )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KIND returns the "kind" of a given integer.'
  write ( *, '(a)' ) '  Do you expect a value of 4, or of 8?'
  write ( *, '(a,i12)' ) '    KIND(I1) = ', kind ( a )

  a = cmplx ( 1.0E+00, 1.0E+00 )
  b = sqrt ( a )
  c = a - b * b

  write ( *, '(a)' ) ' '
  write ( *, * ) '  A = ', a
  write ( *, * ) '  B = sqrt ( A ) = ', b
  write ( *, * ) '  C = A - B * B = ', c

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 shows what you get with KIND = 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 8 ) a
  complex ( kind = 8 ) b
  complex ( kind = 8 ) c

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Use complex values of KIND = 8.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RANGE provides the decimal exponent range.'
  write ( *, '(a,i12)' ) '    RANGE(A) = ', range ( a )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KIND returns the "kind" of a given integer.'
  write ( *, '(a,i12)' ) '    KIND(A) = ', kind ( a )

  a = cmplx ( 1.0D+00, 1.0D+00, kind = 8 )
  b = sqrt ( a )
  c = a - b * b

  write ( *, '(a)' ) ' '
  write ( *, * ) '  A = ', a
  write ( *, * ) '  B = sqrt ( A ) = ', b
  write ( *, * ) '  C = A - B * B = ', c

  return
end
subroutine test03

!*****************************************************************************80
!
!! TEST03 shows what happens if you ask for "DOUBLE COMPLEX".
!
!  Discussion:
!
!    "Double complex" is a nonstandard data type that was implemented
!    on some machines for FORTRAN77.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  double complex a
  double complex b
  double complex c

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Use "DOUBLE COMPLEX" values.'

  a = cmplx ( 1.0D+00, 1.0D+00 )
  b = sqrt ( a )
  c = a - b * b

  write ( *, '(a)' ) ' '
  write ( *, * ) '  A = ', a
  write ( *, * ) '  B = sqrt ( A ) = ', b
  write ( *, * ) '  C = A - B * B = ', c

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
