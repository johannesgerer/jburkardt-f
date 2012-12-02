program main

!*****************************************************************************80
!
!! ELLIPSE_GRID_TEST tests ELLIPSE_GRID.
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

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELLIPSE_GRID_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ELLIPSE_GRID library.'

  call ellipse_grid_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELLIPSE_GRID_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine ellipse_grid_test01 ( )

!*****************************************************************************80
!
!! ELLIPSE_GRID_TEST01 tests ELLIPSE_GRID.
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

  real ( kind = 8 ) c(2)
  character ( len = 80 ) filename
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ng
  real ( kind = 8 ) r(2)
  real ( kind = 8 ), allocatable :: xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  ELLIPSE_GRID can define a grid of points'
  write ( *, '(a)' ) '  with N+1 points on the minor half axis,'
  write ( *, '(a)' ) '  based on any ellipse.'

  n = 8
  r(1) = 2.0D+00
  r(2) = 1.0D+00
  c(1) = 1.0D+00
  c(2) = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  We use N = ', n
  write ( *, '(a,g14.6,a,g14.6,a)' ) '  Radius R = (', r(1), ',', r(2), ')'
  write ( *, '(a,g14.6,a,g14.6,a)' ) '  Center C = (', c(1), ',', c(2), ')'

  call ellipse_grid_count ( n, r, c, ng )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of grid points will be ', ng

  allocate ( xy(2,ng) )

  call ellipse_grid ( n, r, c, ng, xy )

  call r82vec_print_part ( ng, xy, 20, '  Part of the grid point array:' )

  filename = 'ellipse_grid_test01.xy'

  call r8mat_write ( filename, 2, ng, xy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to the file "' // trim ( filename ) // '".'

  deallocate ( xy )

  return
end
