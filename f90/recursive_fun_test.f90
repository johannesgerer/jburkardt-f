program main

!*****************************************************************************80
!
!! MAIN is the main program for RECURSIVE_FUN_TEST.
!
!  Discussion;
!
!    RECURSIVE_FUN_TEST demonstrates the use of recursive functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer f
  integer f_hofstadter
  integer i

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RECURSIVE_FUN_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate recursive function definitions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F_HOFSTADTER evaluates Hofstadter''s recursive'
  write ( *, '(a)' ) '  F function, and does so using an F90 recursive'
  write ( *, '(a)' ) '  function.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N   F(N)'
  write ( *, '(a)' ) ' '

  do i = 0, 30
    f = f_hofstadter ( i )
    write ( *, '(2x,i8,2x,i8)' ) i, f
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RECURSIVE_FUN_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
recursive function f_hofstadter ( n ) result ( value )

!*****************************************************************************80
!
!! F_HOFSTADTER computes the Hofstadter F sequence.
!
!  Discussion:
!
!    F(N) = 0                if N = 0
!         = N - F ( N - 1 ), otherwise.
!
!    F(N) is defined for all nonnegative integers, and turns out
!    to be equal to int ( ( N + 1 ) / 2 ).
!
!    In nonrecursive FORTRAN functions, the name of the function
!    is used as the value of the function.  For a recursive FORTRAN90
!    function, the value is given a separate name, and that name
!    is specified using the extra
!
!      result ( "NAME" )
!
!    field in the function declaration.  In the body of the function,
!    the assignment is made to the variable defined in the "result"
!    field, and NOT to a variable that shares its name with that
!    of the function.
!
!  Table:
!
!     N  F(N)
!    --  ----
!
!     0     0
!     1     1
!     2     1
!     3     2
!     4     2
!     5     3
!     6     3
!     7     4
!     8     4
!     9     5
!    10     5
!    11     6
!    12     6
!    13     7
!    14     7
!    15     8
!    16     8
!    17     9
!    18     9
!    19    10
!    20    10
!    21    11
!    22    11
!    23    12
!    24    12
!    25    13
!    26    13
!    27    14
!    28    14
!    29    15
!    30    15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Douglas Hofstadter,
!    Goedel, Escher, Bach,
!    Basic Books, 1979.
!
!  Parameters:
!
!    Input, integer N, the argument of the function.
!
!    Output, integer VALUE, the value of the function.
!
  implicit none

  integer n
  integer value

  if ( n <= 0 ) then
    value = 0
  else
    value = n - f_hofstadter ( n - 1 )
  end if

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
