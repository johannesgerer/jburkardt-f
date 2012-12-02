program main

!*****************************************************************************80
!
!! MAIN is the main program for BISECTION_INTEGER_PRB.
!
!  Discussion:
!
!    BISECTION_INTEGER_PRB tests the BISECTION_INTEGER library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BISECTION_INTEGER_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BISECTION_INTEGER library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BISECTION_INTEGER_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BISECTION_INTEGER;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ), external :: f01
  integer ( kind = 4 ) fc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  BISECTION_INTEGER attempts to locate an integer root C'
  write ( *, '(a)' ) '  of an equation F(C) = 0.'
  write ( *, '(a)' ) '  The user supplies a change of sign interval [A,B].'
  write ( *, '(a)' ) '  The function considered here has two real roots'
  write ( *, '(a)' ) '  as well as an integer root, so the algorithm can'
  write ( *, '(a)' ) '  fail depending on how the change of sign interval is chosen.'

  a = 4
  b = 100

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The initial change of sign interval is:'
  write ( *, '(a,i8,a,i8)' ) '  F(', a, ') = ', f01 ( a )
  write ( *, '(a,i8,a,i8)' ) '  F(', b, ') = ', f01 ( b )

  call bisection_integer ( f01, a, b, c, fc )

  if ( fc == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  An exact root was found at C = ', c
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  An exact root was NOT found.'
    write ( *, '(a,i8)' ) '  The change of sign interval is now:'
    write ( *, '(a,i8,a,i8)' ) '  F(', a, ') = ', f01 ( a )
    write ( *, '(a,i8,a,i8)' ) '  F(', b, ') = ', f01 ( b )
  end if

  a = -10
  b = 15

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The initial change of sign interval is:'
  write ( *, '(a,i8,a,i8)' ) '  F(', a, ') = ', f01 ( a )
  write ( *, '(a,i8,a,i8)' ) '  F(', b, ') = ', f01 ( b )

  call bisection_integer ( f01, a, b, c, fc )

  if ( fc == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  An exact root was found at C = ', c
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  An exact root was NOT found.'
    write ( *, '(a,i8)' ) '  The change of sign interval is now:'
    write ( *, '(a,i8,a,i8)' ) '  F(', a, ') = ', f01 ( a )
    write ( *, '(a,i8,a,i8)' ) '  F(', b, ') = ', f01 ( b )
  end if


  return
end
function f01 ( n )

!*****************************************************************************80
!
!! F01 is a test function.
!
!  Discussion:
!
!    The polynomial has roots 1/2, 7/2, and 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument.
!
!    Output, integer F01, the function value.
!
  implicit none

  integer ( kind = 4 ) f01
  integer ( kind = 4 ) n

  f01 = ( 2 * n - 7 ) * ( 2 * n - 1 ) * ( n - 10 )

  return
end
