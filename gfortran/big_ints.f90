program main

!*****************************************************************************80
!
!! MAIN is the main program for BIG_INTS.
!
!  Discussion:
!
!    BIG_INTS demonstrates how to set up integers with a bigger range.
!
!    Although you can try specifying
!
!      integer ( kind = 8 ) :: i
!
!    to get a 64 bit integer, the FORTRAN90 standard apparently
!    doesn't specify the actual values associated with KIND, leaving
!    that up to individual compiler writers.  The natural choice
!    seems to be to let KIND represent the number of bytes in
!    the representation, but this is not required.
!
!    What is "guaranteed" to work is something like
!
!      integer, parameter :: kind_val = selected_int_kind ( 20 )
!
!    which asks that KIND_VAL be set to the appropriate KIND value
!    to be able to represent integers with up to 20 decimal digits.
!    (Of course, if there is no integer kind that can handle 20
!    digits, you're out of luck.)
!
!    Then, you use KIND_VAL in your declarations:
!
!      integer ( kind_val ) :: i
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BIG_INTS'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate the use of "big" integers.'

  call test01
  call test02
  call test03
  call test04

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BIG_INTS'
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
!    27 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use integers of the default type.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HUGE() returns the largest value of the given type.'
  write ( *, '(a,i24)' ) '    HUGE(I1) = ', huge ( i1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RANGE provides the decimal exponent range.'
  write ( *, '(a,i24)' ) '    RANGE(I1) = ', range ( i1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KIND returns the "kind" of a given integer.'
  write ( *, '(a,i24)' ) '    KIND(I1) = ', kind ( i1 )

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 shows what you get with KIND = 4.
!
!  Discussion:
!
!    Using "4" or "8" as the kind of an integer is just an "inspired guess".
!    It is not part of the standard.  We're just guessing it makes
!    sense to identify a big integer by the number of bytes used to
!    store it.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Use integers of KIND = 4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HUGE() returns the largest value of the given type.'
  write ( *, '(a,i24)' ) '    HUGE(I1) = ', huge ( i1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RANGE provides the decimal exponent range.'
  write ( *, '(a,i24)' ) '    RANGE(I1) = ', range ( i1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KIND returns the "kind" of a given integer.'
  write ( *, '(a,i24)' ) '    KIND(I1) = ', kind ( i1 )

  return
end
subroutine test03

!*****************************************************************************80
!
!! TEST03 shows what you get with KIND = 8.
!
!  Discussion:
!
!    Using "8" as the kind of a big integer is just an "inspired guess".
!    It is not part of the standard.  We're just guessing it makes
!    sense to identify a big integer by the number of bytes used to
!    store it.  
!
!    On the other hand, the correct way of determining the argument
!    for KIND, (which will turn out to be 8, almost surely) follows
!    in the next subroutine, and is not particularly memorable,
!    although correct.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 8 ) i1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Use integers of KIND = 8.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HUGE() returns the largest value of the given type.'
  write ( *, '(a,i24)' ) '    HUGE(I1) = ', huge ( i1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RANGE provides the decimal exponent range.'
  write ( *, '(a,i24)' ) '    RANGE(I1) = ', range ( i1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KIND returns the "kind" of a given integer.'
  write ( *, '(a,i24)' ) '    KIND(I1) = ', kind ( i1 )

  return
end
subroutine test04

!*****************************************************************************80
!
!! TEST04 shows what you get using SELECTED_INT_KIND.
!
!  Discussion:
!
!    This example shows the recommended procedure for getting
!    a big integer, which is to decide how many integer digits 
!    you'd like to use, and then politely asking for the integer
!    type, if any, that would have that many digits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: digits = 18

  integer, parameter :: kind_val = selected_int_kind ( digits )

  integer ( kind_val ) :: i1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a,i3,a)' ) '  Set KIND_VAL = SELECTED_INT_KIND ( ', digits, ' )'
  write ( *, '(a)' ) '  and declare I1 of type " integer ( KIND_VAL )"'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i24)' ) '  KIND_VAL = ', kind_val
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HUGE() returns the largest value of the given type.'
  write ( *, '(a,i24)' ) '    HUGE(I1) = ', huge ( i1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RANGE provides the decimal exponent range.'
  write ( *, '(a,i24)' ) '    RANGE(I1) = ', range ( i1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KIND returns the "kind" of a given integer.'
  write ( *, '(a,i24)' ) '    KIND(I1) = ', kind ( i1 )

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
