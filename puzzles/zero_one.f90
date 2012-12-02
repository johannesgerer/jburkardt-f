program main

!*****************************************************************************80
!
!! ZERO_ONE tests some ways of defining a function returning 1 if I == J.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZERO_ONE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tests of functions designed to return 1 if I == J.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The challenge is to design a one line arithmetic'
  write ( *, '(a)' ) '  function which accepts two integer arguments I and J,'
  write ( *, '(a)' ) '  and returns 1 if they are equal and 0 otherwise.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  While solutions are possible involving logical '
  write ( *, '(a)' ) '  operations, conversions, and other sophisticated means,'
  write ( *, '(a)' ) '  a solution is preferred that is as simple and close to '
  write ( *, '(a)' ) '  "number-theoretical" as possible.'

  call test00
  call test01
! call test02
  call test03
  call test04
  call test05
  call test06
  call test07

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZERO_ONE'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test00

!*****************************************************************************80
!
!! TEST00 simply assigns the desired values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: min_val = -10
  integer, parameter :: max_val = +10

  integer a(min_val:max_val,min_val:max_val)
  integer i
  integer j

  do i = min_val, max_val
    do j = min_val, max_val
      if ( i == j ) then
        a(i,j) = 1
      else
        a(i,j) = 0
      end if
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST00'
  write ( *, '(a)' ) '  For comparison only.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  if ( i == j ) then'
  write ( *, '(a)' ) '    a(i,j) = 1'
  write ( *, '(a)' ) '  else'
  write ( *, '(a)' ) '    a(i,j) = 0'
  write ( *, '(a)' ) '  end if'
  write ( *, '(a)' ) ' '

  do i = min_val, max_val
    write ( *, '(4x,41i1)' ) a(i,min_val:max_val)
  end do

  return
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 is the first test.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: min_val = -10
  integer, parameter :: max_val = +10

  integer a(min_val:max_val,min_val:max_val)
  integer i
  integer j

  do i = min_val, max_val
    do j = min_val, max_val
      a(i,j) = 1 - int ( abs ( sign ( 1, i - j ) - sign ( 1, j - i ) ) ) / 2
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  a(i,j) = 1 - int ( abs ( sign ( 1, i - j ) - sign ( 1, j - i ) ) ) / 2'
  write ( *, '(a)' ) ' '

  do i = min_val, max_val
    write ( *, '(4x,41i1)' ) a(i,min_val:max_val)
  end do

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 is the second test.  It has been cancelled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: min_val = -10
  integer, parameter :: max_val = +10

  integer a(min_val:max_val,min_val:max_val)
  integer i
  integer j

  do i = min_val, max_val
    do j = min_val, max_val
!     a(i,j) = int ( i == j )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  a(i,j) = int ( i == j )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TEST CANCELLED:' 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FORTRAN90 will not allow the INT operator to be used'
  write ( *, '(a)' ) '  except on a numeric argument.'
  write ( *, '(a)' ) ' '

  do i = min_val, max_val
    write ( *, '(4x,41i1)' ) a(i,min_val:max_val)
  end do

  return
end
subroutine test03

!*****************************************************************************80
!
!! TEST03 is the third test.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: min_val = 1
  integer, parameter :: max_val = +10

  integer a(min_val:max_val,min_val:max_val)
  integer i
  integer j

  do i = min_val, max_val
    do j = min_val, max_val
      a(i,j) = (i/j) * (j/i)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  a(i,j) = (i/j) * (j/i)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test restricted to I, J greater than 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Thanks to David Levine, St Bonaventure University.'
  write ( *, '(a)' ) ' '

  do i = min_val, max_val
    write ( *, '(4x,41i1)' ) a(i,min_val:max_val)
  end do

  return
end
subroutine test04

!*****************************************************************************80
!
!! TEST04 is the fourth test.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: min_val = -10
  integer, parameter :: max_val = +10

  integer a(min_val:max_val,min_val:max_val)
  integer i
  integer j

  do i = min_val, max_val
    do j = min_val, max_val
      a(i,j) = 1 / ( 1 + ( i - j ) * ( i - j ) )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  a(i,j) = 1 / ( 1 + ( i - j ) * ( i - j ) )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Thanks to David Faden, Iowa State University.'
  write ( *, '(a)' ) ' '

  do i = min_val, max_val
    write ( *, '(4x,41i1)' ) a(i,min_val:max_val)
  end do

  return
end
subroutine test05

!*****************************************************************************80
!
!! TEST05 is the fifth test.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: min_val = -10
  integer, parameter :: max_val = +10

  integer a(min_val:max_val,min_val:max_val)
  integer i
  integer j

  do i = min_val, max_val
    do j = min_val, max_val
      a(i,j) = 1 - min ( abs ( i - j ), 1 )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  a(i,j) = 1 - min ( abs ( i - j ), 1 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Thanks to David Faden, Iowa State University.'
  write ( *, '(a)' ) ' '

  do i = min_val, max_val
    write ( *, '(4x,41i1)' ) a(i,min_val:max_val)
  end do

  return
end
subroutine test06

!*****************************************************************************80
!
!! TEST06 is the sixth test.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: min_val = -10
  integer, parameter :: max_val = +10

  integer a(min_val:max_val,min_val:max_val)
  integer i
  integer j

  do i = min_val, max_val
    do j = min_val, max_val
      a(i,j) = 1 - ( abs ( i - j ) + 1 - abs ( abs ( i - j ) - 1 ) ) / 2
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  a(i,j) = 1 - ( abs ( i - j ) + '
  write ( *, '(a)' ) '           1 - abs ( abs ( i - j ) - 1 ) ) / 2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Thanks to David Faden, Iowa State University.'
  write ( *, '(a)' ) ' '

  do i = min_val, max_val
    write ( *, '(4x,41i1)' ) a(i,min_val:max_val)
  end do

  return
end
subroutine test07

!*****************************************************************************80
!
!! TEST07 is the seventh test.
!
!  Discussion:
!
!    Note that conversion and rounding are used in this example.
!
!    This can be seen by writing out explicitly what is happening:
!
!      a(i,j) = int ( real ( 1 ) - real ( i - j ) / ( real ( i - j ) + 0.1 ) )
!
!    Note, also, that it is important here that the final conversion
!    to an integer always rounds towards 0, so you get the same
!    result for a(i,j) and a(j,i).  For instance,
!
!      a(2,1) = int ( 1.0 - 1.0 / 1.1 ) 
!             = int ( 0.09090909... )
!             = 0
!
!      a(1,2) = int ( 1.0 - ( - 1.0 ) / ( -0.9 ) )
!             = int ( 1.0 - 1.0 / 0.9 ) 
!             = int ( -0.09090909... ) 
!             = 0
!
!    Moreover, the discrepency, that is, the amount that has to be
!    rounded, is at most 1/9, the value achieved when I and J differ
!    by 1.  As the difference between I and J increases, the discrepancy
!    goes to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2007
!
!  Author:
! 
!    John Burkardt
!
  implicit none

  integer, parameter :: min_val = -10
  integer, parameter :: max_val = +10

  integer a(min_val:max_val,min_val:max_val)
  integer i
  integer j

  do i = min_val, max_val
    do j = min_val, max_val
      a(i,j) = 1 - ( i - j ) / ( i - j + 0.1 )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  a(i,j) = 1 - ( i - j ) / ( i - j + 0.1 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Thanks to Li Wang, Washington State University'
  write ( *, '(a)' ) ' '

  do i = min_val, max_val
    write ( *, '(4x,41i1)' ) a(i,min_val:max_val)
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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2001
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
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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
