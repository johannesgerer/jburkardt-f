program main

!*****************************************************************************80
!
!! MAIN is the main program for SUM_MILLION.
!
!  Discussion:
!
!    This code estimates the power of a computer by summing the integers
!    from 1 to 1,000,000.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 1000000

  real ( kind = 8 ) ctime
  real ( kind = 8 ) error
  real ( kind = 8 ), parameter :: exact = 500000500000.0D+00
  integer ( kind = 4 ) i
  real ( kind = 8 ) mflops
  real ( kind = 8 ) total
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SUM_MILLION'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Sum the integers from 1 to 1,000,000.'
  write ( *, '(a)' ) '  Correct answer is 500000500000.'

  call set_up ( n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N      CPU time        MFLOPS         ERROR'
  write ( *, '(a)' ) '                (seconds)'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call sum_up ( n, x, total, ctime )

    mflops = real ( n, kind = 8 ) / 1000000.0D+00 / ctime

    error = total - exact
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) n, ctime, mflops, error

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SUM_MILLION:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine set_up ( n, x )

!*****************************************************************************80
!
!! SET_UP sets up the data for the SUM_MILLION program.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values to define.
!
!    Output, real ( kind = 8 ) X(N), a vector containing the values 1 through N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  return
end
subroutine sum_up ( n, x, total, ctime )

!*****************************************************************************80
!
!! SUM_UP carries out the sum for the SUM_MILLION program.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values to define.
!
!    Input, real ( kind = 8 ) X(N), the data to be summed.
!
!    Output, real ( kind = 8 ) TOTAL, the sum of the data.
!
!    Output, real ( kind = 8 ) CTIME, the cpu time required to sum the data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ctime
  real ( kind = 8 ) ctime1
  real ( kind = 8 ) ctime2
  real ( kind = 8 ) total
  real ( kind = 8 ) x(n)

  call cpu_time ( ctime1 )

  total = sum ( x(1:n) )

  call cpu_time ( ctime2 )

  ctime = ctime2 - ctime1

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
  intege ( kind = 4 )r d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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
