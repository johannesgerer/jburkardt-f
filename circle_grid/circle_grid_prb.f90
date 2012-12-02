program main

!*****************************************************************************80
!
!! CIRCLE_GRID_TEST tests CIRCLE_GRID.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIRCLE_GRID_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CIRCLE_GRID library.'

  call circle_grid_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIRCLE_GRID_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine circle_grid_test01 ( )

!*****************************************************************************80
!
!! CIRCLE_GRID_TEST01 tests CIRCLE_GRID.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) c(2)
  real ( kind = 8 ), allocatable :: cg(:,:)
  character ( len = 80 ) filename
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ng
  real ( kind = 8 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  CIRCLE_GRID can define a grid of points'
  write ( *, '(a)' ) '  with N+1 points on a horizontal or vertical radius,'
  write ( *, '(a)' ) '  based on any circle.'

  n = 20
  r = 2.0D+00
  c(1) = 1.0D+00
  c(2) = 5.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  We use N = ', n
  write ( *, '(a,g14.6)' ) '  Radius R = ', r
  write ( *, '(a,g14.6,a,g14.6,a)' ) '  Center C = (', c(1), ',', c(2), ')'

  call circle_grid_count ( n, r, c, ng )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of grid points will be ', ng

  allocate ( cg(1:2,1:ng) )

  call circle_grid ( n, r, c, ng, cg )

  call r82vec_print_part ( ng, cg, 20, '  Part of the grid point array:' )

  filename = 'circle_grid_test01.xy'

  call r8mat_write ( filename, 2, ng, cg )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to the file "' // trim ( filename ) // '".'

  deallocate ( cg )

  return
end
