program main

!*****************************************************************************80
!
!! MAIN is the main program for GRID_PRB.
!
!  Discussion:
!
!    GRID_PRB tests the grid routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) center
  integer ( kind = 4 ) seed

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRID_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the GRID library.'

  center = 1
  seed = 123456789

  call test01 ( center, seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat TEST01 with a different seed from the first run.'

  seed = 987654321
  call test01 ( center, seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat TEST01 with the same seed as the first run.'

  seed = 123456789
  call test01 ( center, seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat TEST01 with different centering values.'

  do center = 1, 5

    seed = 123456789
    call test01 ( center, seed )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRID_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( center, seed )

!*****************************************************************************80
!
!! TEST01 tests GRID_GENERATE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: point_num = 10

  integer ( kind = 4 ) center
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(dim_num,point_num)

  write ( *, '(a)'     ) ' '
  write ( *, '(a)'     ) 'TEST01'
  write ( *, '(a)'     ) '  GRID_GENERATE randomly chooses a given number of'
  write ( *, '(a)'     ) '  points on a uniform grid.'
  write ( *, '(a)'     ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension =  ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points =   ', point_num
  write ( *, '(a,i12)' ) '  Random number SEED = ', seed
  write ( *, '(a,i8)'  ) '  Centering option =   ', center

  call grid_generate ( dim_num, point_num, center, seed, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The grid points:'
  write ( *, '(a)' ) ' '

  do j = 1, point_num
    write ( *, '(2f10.4)' ) x(1:dim_num,j)
  end do

  return
end
