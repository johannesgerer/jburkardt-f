program main

!*****************************************************************************80
!
!! MAIN is the main program for UNICYCLE_PRB.
!
!  Discussion:
!
!    UNICYCLE_PRB tests the UNICYCLE library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNICYCLE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the UNICYCLE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNICYCLE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests PERM_IS_UNICYCLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: test_num = 10

  integer ( kind = 4 ) p(n)
  logical perm_is_unicycle
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) u(n)
  logical value

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  PERM_IS_UNICYCLE determines whether a permutation'
  write ( *, '(a)' ) '  is a unicyle'

  do test = 1, test_num

    call perm_random ( n, seed, p )

    value = perm_is_unicycle ( n, p )

    if ( value ) then

      call perm_print ( n, p, '  This permutation is a unicycle' )
      call unicycle_index_to_sequence ( n, p, u )
      call unicycle_print ( n, u, '  The permutation in sequence form' )

    else

      call perm_print ( n, p, '  This permutation is NOT a unicycle' )

    end if

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests UNICYCLE_ENUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10

  integer ( kind = 4 ) n
  integer ( kind = 4 ) num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  UNICYCLE_ENUM enumerates the unicycles of N objects.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N    Number'
  write ( *, '(a)' ) ' '

  do n = 0, n_max

    call unicycle_enum ( n, num )
    write ( *, '(2x,i3,2x,i8)' ) n, num

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests UNICYCLE_INVERSE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ), dimension ( n ) :: u = (/ 1, 7, 6, 2, 4, 3, 5 /)
  integer ( kind = 4 ) u_inverse(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  UNICYCLE_INVERSE inverts a unicycle;'

  call unicycle_print ( n, u, '  The original unicycle:' )
 
  call unicycle_inverse ( n, u, u_inverse )
 
  call unicycle_print ( n, u_inverse, '  The inverse unicycle:' )
 
  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests UNICYCLE_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer rank
  integer ( kind = 4 ) u(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  UNICYCLE_NEXT generates unicycles in lex order.'
  write ( *, '(a)' ) ' '
  rank = -1
 
  do

    call unicycle_next ( n, u, rank )

    if ( rank == - 1 ) then
      exit
    end if

    write ( *, '(2x,i3,a1,2x,10i2)' ) rank, ':', u(1:n)

  end do
 
  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests UNICYCLE_RANDOM;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) u(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  UNICYCLE_RANDOM produces a random unicyle;'
  write ( *, '(a,i8)' ) '  For this test, N = ', n
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 5
    call unicycle_random ( n, seed, u )
    call unicycle_print ( n, u, ' ' )
  end do
 
  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests UNICYCLE_RANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ), save, dimension ( n ) :: u = (/ 1, 5, 2, 3, 4 /)
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  UNICYCLE_RANK ranks a unicycle.'

  call unicycle_print ( n, u, '  The unicycle:' )
 
  call unicycle_rank ( n, u, rank )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The rank is:', rank
 
  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests UNICYCLE_UNRANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) u(n)
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  UNICYCLE_UNRANK, given a rank, computes the'
  write ( *, '(a)' ) '  corresponding unicycle.'
  write ( *, '(a)' ) ' '
  rank = 6
  write ( *, '(a,i8)' ) '  The requested rank is ', rank
 
  call unicycle_unrank ( n, rank, u )
 
  call unicycle_print ( n, u, '  The unicycle:' )
 
  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests UNICYCLE_INDEX, UNICYCLE_INDEX_TO_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ) u(n)
  integer ( kind = 4 ) u_index(n)
  integer ( kind = 4 ) u2(n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) :: test_num = 5

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  UNICYCLE_INDEX converts a unicycle to index form.'
  write ( *, '(a)' ) '  UNICYCLE_INDEX_TO_SEQUENCE converts an index to unicycle form.'

  do test = 1, test_num 

    call unicycle_random ( n, seed, u )

    call unicycle_print ( n, u, '  The unicycle:' )

    call unicycle_index ( n, u, u_index )
    
    call unicycle_index_print ( n, u_index, '  The index form:' )

    call unicycle_index_to_sequence ( n, u_index, u2 )

    call unicycle_print ( n, u2, '  The unicycle recovered:' )

  end do

  return
end
