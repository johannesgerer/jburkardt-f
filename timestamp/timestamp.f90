subroutine hms_current_hms ( h, m, s, mm )

!*****************************************************************************80
!
!! HMS_CURRENT_HMS returns the current HMS time as integers.
!
!  Example:
!
!    If the current time is 9:45:54.872 AM, then
!
!    H = 9 
!    M = 45
!    S = 54
!    MM = 872
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
!    Output, integer ( kind = 4 ) H, M, S, MM, the current hour, minute, 
!    second, and thousandths of a second.
!
  implicit none

  integer ( kind = 4 ) h
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)

  call date_and_time ( values = values )

  h = values(5)
  m = values(6)
  s = values(7)
  mm = values(8)

  return
end
subroutine hms_current_print ( string )

!*****************************************************************************80
!
!! HMS_CURRENT_PRINT prints the current HMS time, and a user specified string.
!
!  Example:
!
!     Wallclock:  9:45:54.872 AM  Started determinant calculation.
!     Wallclock:  9:47:32.738 AM  Finished determinant calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string to be printed.
!
  implicit none

  character ( len = 15 ) string2
  character ( len = *  ) string

  call hms_current_string ( string2 )

  write ( *, '(a,2x,a,2x,a)' ) 'Wallclock:', string2, trim ( string )

  return
end
subroutine hms_current_string ( string )

!*****************************************************************************80
!
!! HMS_CURRENT_STRING writes the current HMS data into a string.
!
!  Example:
!
!    STRING = ' 9:45:54.872 AM'
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
!    Output, character ( len = * ) STRING, contains the HMS information.
!    A character length of 15 should always be sufficient.
!
  implicit none

  character ( len = 2  ) ampm
  integer ( kind = 4 ) h
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = *  ) string
  integer ( kind = 4 ) values(8)

  call date_and_time ( values = values )

  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Nn'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Md'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( string, '(i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine hms_delta_print ( string )

!*****************************************************************************80
!
!! HMS_DELTA_PRINT prints the change in HMS time, and a user specified string.
!
!  Example:
!
!    Delta Wallclock:  0:01:37.966 AM  Determinant calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string to be printed.
!
  implicit none

  integer ( kind = 4 ), save :: h = -1
  integer ( kind = 4 ) h_del
  integer ( kind = 4 ) h_old
  integer ( kind = 4 ), save :: m = 0
  integer ( kind = 4 ) m_del
  integer ( kind = 4 ) m_old
  integer ( kind = 4 ), save :: mm = 0
  integer ( kind = 4 ) mm_del
  integer ( kind = 4 ) mm_old
  integer ( kind = 4 ), save :: s = 0
  integer ( kind = 4 ) s_del
  integer ( kind = 4 ) s_old
  character ( len = *  ) string
!
!  Back up the previous time.
!
  if ( h == -1 ) then

    call hms_current_hms ( h, m, s, mm )

    h_old = h
    m_old = m
    s_old = s
    mm_old = mm

  else

    h_old = h
    m_old = m
    s_old = s
    mm_old = mm

    call hms_current_hms ( h, m, s, mm )

  end if

  h_del = h - h_old
  m_del = m - m_old
  s_del = s - s_old
  mm_del = mm - mm_old

  if ( mm_del < 0 ) then
    s_del = s_del - 1
    mm_del = mm_del + 1000
  end if

  if ( s_del < 0 ) then
    m_del = m_del - 1
    s_del = s_del + 60
  end if

  if ( m_del < 0 ) then
    m_del = m_del + 60
    h_del = h_del - 1
  end if

  if ( h_del < 0 ) then
    h_del = h_del + 24
  end if

  write ( *, '(a,i2,a1,i2.2,a1,i2.2,a1,i3.3,2x,a)' ) &
    'Delta Wallclock: ', h_del, ':', m_del, ':', s_del, '.', mm_del, &
    trim ( string )

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
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = '31 May 2001   9:45:54.872 AM'
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
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  implicit none

  character ( len = 8  ) ampm
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
  character ( len = *  ) string
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

  write ( string, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
