function i4_div_rounded ( a, b )

!*****************************************************************************80
!
!! I4_DIV_ROUNDED computes the rounded result of I4 division.
!
!  Discussion:
!
!    This routine computes C = A / B, where A, B and C are integers
!    and C is the closest integer value to the exact real result.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number to be divided.
!
!    Input, integer ( kind = 4 ) B, the divisor, or the number of parts.
!
!    Output, integer ( kind = 4 ) I4_DIV_ROUNDED, the rounded result
!    of the division.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) a_abs
  integer ( kind = 4 ) b
  integer ( kind = 4 ) b_abs
  integer ( kind = 4 ) i4_div_rounded
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) value

  if ( a == 0 .and. b == 0 ) then

    value = i4_huge
 
  else if ( a == 0 ) then

    value = 0

  else if ( b == 0 ) then

    if ( a < 0 ) then
      value = - i4_huge
    else
      value = + i4_huge
    end if

  else

    a_abs = abs ( a )
    b_abs = abs ( b )

    value = a_abs / b_abs
!
!  Round the value.
!
    if ( ( 2 * value + 1 ) * b_abs < 2 * a_abs ) then
      value = value + 1
    end if
!
!  Set the sign.
!
    if ( ( a < 0 .and. 0 < b ) .or. ( 0 < a .and. b < 0 ) ) then
      value = - value
    end if

  end if

  i4_div_rounded = value

  return
end
subroutine task_division ( task_number, proc_first, proc_last )

!*****************************************************************************80
!
!! TASK_DIVISION divides tasks among processors.
!
!  Discussion:
!
!    This routine assigns each of T tasks to P processors, assuming that 
!    the assignment is to be beforehand.
!
!    In that case, we just want to make sure that we assign each task
!    to a processor, that we assign about the same number of tasks
!    to each processor, and that we assign each processor a contiguous
!    range of tasks, say tasks I_LO to I_HI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TASK_NUMBER, the number of tasks.
!
!    Input, integer PROC_FIRST, PROC_LAST, the first and last processors.
!
  implicit none

  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  integer ( kind = 4 ) i4_div_rounded
  integer ( kind = 4 ) p
  integer ( kind = 4 ) proc
  integer ( kind = 4 ) proc_first
  integer ( kind = 4 ) proc_last
  integer ( kind = 4 ) proc_number
  integer ( kind = 4 ) proc_remain
  integer ( kind = 4 ) task_number
  integer ( kind = 4 ) task_proc
  integer ( kind = 4 ) task_remain

  p = proc_last + 1 - proc_first

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TASK_DIVISION'
  write ( *, '(a)' ) '  Divide T tasks among P processors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of tasks T = ', task_number
  write ( *, '(a,i8)' ) '  Number of processors P = ', p
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  P_FIRST = ', proc_first
  write ( *, '(a,i8)' ) '  P_LAST =  ', proc_last
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             Number of   First      Last'
  write ( *, '(a)' ) ' Processor     Tasks     Task       Task'
  write ( *, '(a)' ) ' '

  i_hi = 0

  task_remain = task_number
  proc_remain = p

  do proc = proc_first, proc_last

    task_proc = i4_div_rounded ( task_remain, proc_remain )

    proc_remain = proc_remain - 1
    task_remain = task_remain - task_proc

    i_lo = i_hi + 1
    i_hi = i_hi + task_proc

    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) proc, task_proc, i_lo, i_hi

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
  integer ( kind = 4 ) d
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
