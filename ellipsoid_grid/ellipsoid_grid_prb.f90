program main

!*****************************************************************************80
!
!! ELLIPSOID_GRID_TEST tests ELLIPSOID_GRID.
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
  write ( *, '(a)' ) 'ELLIPSOID_GRID_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ELLIPSOID_GRID library.'

  call ellipsoid_grid_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELLIPSOID_GRID_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine ellipsoid_grid_test01 ( )

!*****************************************************************************80
!
!! ELLIPSOID_GRID_TEST01 tests ELLIPSOID_GRID.
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

  real ( kind = 8 ) c(3)
  character ( len = 80 ) filename
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ng
  real ( kind = 8 ) r(3)
  real ( kind = 8 ), allocatable :: xyz(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  ELLIPSOID_GRID can define a grid of points'
  write ( *, '(a)' ) '  with N+1 points on the minor half axis,'
  write ( *, '(a)' ) '  based on any ellipsoid.'

  n = 4
  r(1) = 2.0D+00
  r(2) = 1.0D+00
  r(3) = 1.5D+00
  c(1) = 1.0D+00
  c(2) = 2.0D+00
  c(3) = 1.5D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use N = ', n
  write ( *, '(a,g14.6,a,g14.6,a,g14.6,a)' ) &
    '  Radius R = (', r(1), ',', r(2), ',', r(3), ')'
  write ( *, '(a,g14.6,a,g14.6,a,g14.6,a)' ) &
    '  Center C = (', c(1), ',', c(2), ',', c(3), ')'

  call ellipsoid_grid_count ( n, r, c, ng )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of grid points will be ', ng

  allocate ( xyz(1:3,1:ng) )

  call ellipsoid_grid ( n, r, c, ng, xyz )

  call r83vec_print_part ( ng, xyz, 20, '  Part of the grid point array:' )

  filename = 'ellipsoid_grid_test01.xyz'

  call r8mat_write ( filename, 3, ng, xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to the file "' // trim ( filename ) // '".'

  deallocate ( xyz )

  return
end
