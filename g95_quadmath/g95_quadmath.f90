program main

!*****************************************************************************80
!
!! MAIN is the main program for G95_QUADMATH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'G95_QUADMATH'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the G95 quadmath facility.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'G95_QUADMATH'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses REAL ( KIND = 4 ) arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer divs
  real ( kind = 4 ) x
  real ( kind = 4 ) x_old
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Using REAL ( KIND = 4 ) arithmetic:'
  write ( *, '(a)' ) '  Compute smallest 1/2^DIV that can be added to 1.'

  x = 1.0
  z = 1.0
  divs = 0

  do 
    x_old = x
    x = x / 2.0
    y = 1.0 + x
    if ( y <= z ) then
      exit
    end if
    divs = divs + 1
  end do

  write ( *, '(a,i4)' ) '  Number of divisions DIV = ', divs
  write ( *, '(a,g14.6)' ) '  1/2^DIV =         ', x_old
  write ( *, '(a,g14.6)' ) '  Machine epsilon = ', epsilon ( x_old )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 uses REAL ( KIND = 8 ) arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer divs
  real ( kind = 8 ) x
  real ( kind = 8 ) x_old
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Using REAL ( KIND = 8 ) arithmetic:'
  write ( *, '(a)' ) '  Compute smallest 1/2^DIV that can be added to 1.'

  x = 1.0
  z = 1.0
  divs = 0

  do
    x_old = x
    x = x / 2.0
    y = 1.0 + x
    if ( y <= z ) then
      exit
    end if
    divs = divs + 1
  end do

  write ( *, '(a,i4)' ) '  Number of divisions DIV = ', divs
  write ( *, '(a,g14.6)' ) '  1/2^DIV =         ', x_old
  write ( *, '(a,g14.6)' ) '  Machine epsilon = ', epsilon ( x_old )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 uses REAL ( KIND = 10 ) arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer divs
  real ( kind = 10 ) x
  real ( kind = 10 ) x_old
  real ( kind = 10 ) y
  real ( kind = 10 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Using REAL ( KIND = 10 ) arithmetic:'
  write ( *, '(a)' ) '  Compute smallest 1/2^DIV that can be added to 1.'

  x = 1.0
  z = 1.0
  divs = 0

  do 
    x_old = x
    x = x / 2.0
    y = 1.0 + x
    if ( y <= z ) then
      exit
    end if
    divs = divs + 1
  end do

  write ( *, '(a,i4)' ) '  Number of divisions DIV = ', divs
  write ( *, '(a,g14.6)' ) '  1/2^DIV =         ', x_old
  write ( *, '(a,g14.6)' ) '  Machine epsilon = ', epsilon ( x_old )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 uses REAL ( KIND = 16 ) arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer divs
  real ( kind = 16 ) x
  real ( kind = 16 ) x_old
  real ( kind = 16 ) y
  real ( kind = 16 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Using REAL ( KIND = 16 ) arithmetic:'
  write ( *, '(a)' ) '  Compute smallest 1/2^DIV that can be added to 1.'

  x = 1.0
  z = 1.0
  divs = 0

  do 
    x_old = x
    x = x / 2.0
    y = 1.0 + x
    if ( y <= z ) then
      exit
    end if
    divs = divs + 1
  end do

  write ( *, '(a,i4)' ) '  Number of divisions DIV = ', divs
  write ( *, '(a,g14.6)' ) '  1/2^DIV =         ', x_old
  write ( *, '(a,g14.6)' ) '  Machine epsilon = ', epsilon ( x_old )

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
