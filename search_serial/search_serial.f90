program main

!*****************************************************************************80
!
!! MAIN is the main program for SEARCH_SERIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fj
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  real ( kind = 8 ) u
  real ( kind = 8 ) wtime
  real ( kind = 8 ) wtime1
  real ( kind = 8 ) wtime2

  a = 1
  b = i4_huge
  c = 45

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SEARCH_SERIAL:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Search the integers from A to B'
  write ( *, '(a)' ) '  for a value J such that F(J) = C.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  A           = ', a
  write ( *, '(a,i12)' ) '  B           = ', b
  write ( *, '(a,i12)' ) '  C           = ', c

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   J     C=F(J)      U=C/(maximum integer)'
  write ( *, '(a)' ) ' '
  do j = 0, 20
    c = f ( j )
    u = real ( c, kind = 8 ) / real ( i4_huge, kind = 8 )
    write ( *, '(2x,i2,2x,i12,2x,g14.6)' ) j, c, u
  end do

  call cpu_time ( wtime1 )

  call search ( a, b, c, j )

  call cpu_time ( wtime2 )
  wtime = wtime2 - wtime

  if ( j == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  No solution was found.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  Found     J = ', j
    write ( *, '(a,i12)' ) '  Verify F(J) = ', f ( j )
  end if

  write ( *, '(a,g14.6)' ) '  Elapsed CPU time is ', wtime
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SEARCH_SERIAL:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine search ( a, b, c, j )

!*****************************************************************************80
!
!! SEARCH searches integers in [A,B] for a J so that F(J) = C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the search range.
!
!    Input, integer ( kind = 4 ) C, the desired function value.
!
!    Output, integer ( kind = 4 ) J, the computed solution, or -1
!    if no solution was found.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  j = -1

  do i = a, b

    fi = f ( i )

    if ( fi == c ) then
      j = i
!     return
    end if

  end do

  return
end
function f ( i )

!*****************************************************************************80
!
!! F is the function we are analyzing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the argument.
!
!    Input, integer ( kind = 4 ) F, the value.
!
  implicit none

  integer ( kind = 4 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) value

  value = i

  do j = 1, 5

    k = value / 127773

    value = 16807 * ( value - k * 127773 ) - k * 2836

    if ( value < 0 ) then
      value = value + i4_huge
    end if

  end do

  f = value

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
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
