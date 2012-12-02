program main

!*****************************************************************************80
!
!! MAIN is the main program for DIVISION.
!
!  Discussion:
!
!    DIVISION demonstrates that division should be done to full precision.
!
!    Especially when carrying out division, but also in other arithmetic
!    operations, it is important that constants be specified as double
!    precision values, by using the "D" identifier in the exponent field.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ) a(4,4)
  real    ( kind = 8 ) b(4,4)
  integer ( kind = 4 ) i

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DIVISION:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate 16 ways to compute 1/3, expecting'
  write ( *, '(a)' ) '  a double precision real result.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The results, stored as a table, are:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  1       / 3   1       / 3.0   1       / 3.0E+00  1       / 3.0D+00'
  write ( *, '(a)' ) &
    '  1.0     / 3   1.0     / 3.0   1.0     / 3.0E+00  1.0     / 3.0D+00'
  write ( *, '(a)' ) &
    '  1.0E+00 / 3   1.0E+00 / 3.0   1.0E+00 / 3.0E+00  1.0E+00 / 3.0D+00'
  write ( *, '(a)' ) &
    '  1.0D+00 / 3   1.0D+00 / 3.0   1.0E+00 / 3.0E+00  1.0D+00 / 3.0D+00'

  a(1,1) = 1 / 3
  a(1,2) = 1 / 3.0
  a(1,3) = 1 / 3.0E+00
  a(1,4) = 1 / 3.0D+00

  a(2,1) = 1.0 / 3
  a(2,2) = 1.0 / 3.0
  a(2,3) = 1.0 / 3.0E+00
  a(2,4) = 1.0 / 3.0D+00

  a(3,1) = 1.0E+00 / 3
  a(3,2) = 1.0E+00 / 3.0
  a(3,3) = 1.0E+00 / 3.0E+00
  a(3,4) = 1.0E+00 / 3.0D+00

  a(4,1) = 1.0D+00 / 3
  a(4,2) = 1.0D+00 / 3.0
  a(4,3) = 1.0D+00 / 3.0E+00
  a(4,4) = 1.0D+00 / 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix of results:'
  write ( *, '(a)' ) ' '

  do i = 1, 4
    write ( *, '(4g20.14)' ) a(i,1:4)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Defects might be more obvious if we multiply'
  write ( *, '(a)' ) '  by 3:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  B = 3 * A:'
  write ( *, '(a)' ) ' '

  b(1:4,1:4) = 3 * a(1:4,1:4)

  do i = 1, 4
    write ( *, '(4g20.14)' ) b(i,1:4)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  B = 3.0 * A:'
  write ( *, '(a)' ) ' '

  b(1:4,1:4) = 3.0 * a(1:4,1:4)

  do i = 1, 4
    write ( *, '(4g20.14)' ) b(i,1:4)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  B = 3.0E+00 * A:'
  write ( *, '(a)' ) ' '

  b(1:4,1:4) = 3.0E+00 * a(1:4,1:4)

  do i = 1, 4
    write ( *, '(4g20.14)' ) b(i,1:4)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  B = 3.0D+00 * A:'
  write ( *, '(a)' ) ' '

  b(1:4,1:4) = 3.0D+00 * a(1:4,1:4)

  do i = 1, 4
    write ( *, '(4g20.14)' ) b(i,1:4)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DIVISION:'
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
