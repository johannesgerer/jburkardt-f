program main

!*****************************************************************************80
!
!! BROWNIAN_MOTION_SIMULATION_PRB tests the BROWNIAN_MOTION_SIMULATION library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ), allocatable :: dsq(:,:)
  character ( len = 80 ) header
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable :: x(:,:)

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BROWNIAN_MOTION_SIMULATION_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BROWNIAN_MOTION_SIMULATION library.'
!
!  Compute the path of a particle undergoing Brownian motion.
!
  do m = 1, 2

    n = 1001
    d = 10.0D+00
    t = 1.0D+00
    seed = 123456789
    allocate ( x(1:m,1:n) )
    call brownian_motion_simulation ( m, n, d, t, seed, x )
    if ( m == 1 ) then
      header = 'motion_1d'
    else if ( m == 2 ) then
      header = 'motion_2d'
    end if
    call brownian_motion_display ( m, n, x, header )
    deallocate ( x )

  end do
!
!  Estimate the average displacement of the particle from the origin
!  as a function of time.
!
  do m = 1, 3

    k = 40
    n = 1001
    d = 10.0D+00
    t = 1.0D+00
    seed = 123456789

    allocate ( dsq(1:k,1:n) )
    call brownian_displacement_simulation ( k, n, m, d, t, seed, dsq )
    if ( m == 1 ) then
      header = 'displacement_1d'
    else if ( m == 2 ) then
      header = 'displacement_2d'
    else if ( m == 3 ) then
      header = 'displacement_3d'
    end if
    call brownian_displacement_display ( k, n, d, t, dsq, header )
    deallocate ( dsq )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BROWNIAN_MOTION_SIMULATION_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( );

  return
end
