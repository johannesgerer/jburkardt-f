program main

!*****************************************************************************80
!
!! TETRAHEDRON_GRID_TEST tests TETRAHEDRON_GRID.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TETRAHEDRON_GRID_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TETRAHEDRON_GRID library.'

  call tetrahedron_grid_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TETRAHEDRON_GRID_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine tetrahedron_grid_test01 ( )

!*****************************************************************************80
!
!! TETRAHEDRON_GRID_TEST01 tests TETRAHEDRON_GRID.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) filename
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ng
  real ( kind = 8 ) :: t(3,4) = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00 /), (/ 3, 4 /) )
  real ( kind = 8 ), allocatable :: tg(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  TETRAHEDRON_GRID can define a tetrahedral grid'
  write ( *, '(a)' ) '  with N+1 points on a side, based on any tetrahedron.'

  n = 10
  write ( *, '(a,i8)' ) '  N = ', n

  call tetrahedron_grid_count ( n, ng )

  call r8mat_print ( 3, 4, t, '  Tetrahedron vertices:' )

  allocate ( tg(1:3,1:ng) )

  call tetrahedron_grid ( n, t, ng, tg )

  call r83vec_print_part ( ng, tg, 20, '  Part of the grid point array:' )

  filename = 'tetrahedron_grid_test01.xyz'

  call r8mat_write ( filename, 3, ng, tg )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to the file "' // trim ( filename ) // '".'

  deallocate ( tg )

  return
end
