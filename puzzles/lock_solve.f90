program main

!*****************************************************************************80
!
!! MAIN is the main program for LOCK_SOLVE.
!
!  Discussion:
!
!    A random 4 digit TARGET is given.
!
!    The user starts at another random 4 digit starting point.
!
!    The user can indicate which digit to be advanced, and by how much.
!    That digit is advanced, modulo 10.
!
!    When the target is achieved, the user has won.
!
!    The problem is that the digits are "sticky".  So advancing digit 1
!    will cause digit 2 to move the same amount.  If digit 2 is advanced,
!    it will cause digits 1 and 3 to move the same.  Similarly, digit 4
!    will cause digit 3 to move as well.
!
!    Although the problem is still solvable, and solvable in at most 4 steps,
!    the complications of stickiness make it harder to see what to do.
!
!    This program applies a formula to determine, given the starting and
!    desired combinations, the amount that each digit must be advanced.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ) digit(4)
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) solve(4)
  integer ( kind = 4 ) target(4)
  integer ( kind = 4 ) x(4)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOCK_SOLVE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Solve the lock combination problem.'

  write ( *, '(a)') ' '
  write ( *, * ) 'Enter the target combination'
  read ( *, * ) target(1:4)
  write ( *, * ) 'Enter the current combination:'
  read ( *, * ) digit(1:4)

  x(1:4) = target(1:4) - digit(1:4)

  solve(1) = i4_modp (   x(1)        - x(3) + x(4), 10 )
  solve(2) = i4_modp (                 x(3) - x(4), 10 )
  solve(3) = i4_modp ( - x(1) + x(2),               10 )
  solve(4) = i4_modp (   x(1) - x(2)        + x(4), 10 )

  write ( *, * ) solve(1:4)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOCK_SOLVE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
