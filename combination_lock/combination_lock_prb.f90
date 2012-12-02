program main

!*****************************************************************************80
!
!! MAIN is the main program for COMBINATION_LOCK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMBINATION_LOCK'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the COMBINATION_LOCK libary.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMBINATION_LOCK'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BICYCLE_LOCK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) step
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  A bicycle combination lock consists of 3 dials,'
  write ( *, '(a)' ) '  each having 10 symbols, 0 through 9.'
  write ( *, '(a)' ) '  We seek to determine the combination C.'
!
!  Report on the problem data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The number of dials is M = ', m
  write ( *, '(a,i4)' ) '  The number of symbols is N = ', n
  write ( *, '(a,i8)' ) &
    '  The number of possible combinations is M^N = ', n ** m

  call get_seed ( seed )
  c = i4_uniform ( 0, 999, seed )

  write ( *, '(a,i3)' ) '  The "secret" combination is ', c

  call bicycle_lock ( c, step )

  if ( step == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The combination was not found!'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The combination was found on step ', step
  end if

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests COMBINATION_LOCK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of dials.
!
!    Input, integer N, the number of symbols on each dial.
!    We assume the symbols are the integers 0 to N-1.
!
!    Input, integer C(M), the combination.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4

  integer ( kind = 4 ), dimension ( m ) :: c = (/ 1, 2, 3, 4 /)
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ) step
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  A combination lock consists of M dials,'
  write ( *, '(a)' ) '  each having N symbols.'
  write ( *, '(a)' ) '  We seek to determine the combination C.'
!
!  Report on the problem data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The number of dials is M = ', m
  write ( *, '(a,i4)' ) '  The number of symbols is N = ', n
  write ( *, '(a,i8)' ) &
    '  The number of possible combinations is M^N = ', n ** m

  call i4vec_print ( m, c, '  The "secret" combination:' );

  call combination_lock ( m, n, c, step )

  if ( step == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The combination was not found!'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The combination was found on step ', step
  end if

  return
end
