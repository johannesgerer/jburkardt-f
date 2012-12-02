program main

!*****************************************************************************80
!
!! POLYGON_MOMENTS_PRB tests POLYGON_MOMENTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLYGON_MOMENTS_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test POLYGON_MOMENTS library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLYGON_MOMENTS_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 carries out a test on a rectangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ), dimension ( 6 ) :: alpha_exact = (/ &
    1.0D+00, &
    5.0D+00, 4.0D+00, &
    30.66666666666667D+00, 22.0D+00, 18.66666666666666D+00 /)
  real ( kind = 8 ) alpha_pq
  integer ( kind = 4 ) k
  real ( kind = 8 ), dimension ( 6 ) :: mu_exact = (/ &
    1.0D+00, &
    0.0D+00, 0.0D+00, &
    5.666666666666667D+00, 2.0D+00, 2.666666666666667D+00 /)
  real ( kind = 8 ) mu_pq
  real ( kind = 8 ), dimension ( 6 ) :: nu_exact = (/ &
    40.0D+00, &
    200.0D+00, 160.0D+00, &
    1226.66666666666667D+00, 880.0D+00, 746.66666666666666D+00 /)
  real ( kind = 8 ) nu_pq
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) s
  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    2.0D+00, 10.0D+00, 8.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( n ) :: y = (/ &
    0.0D+00,  4.0D+00, 8.0D+00, 4.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Check normalized moments of a rectangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   P   Q             Nu(P,Q)'
  write ( *, '(a)' ) '            Computed         Exact'
  write ( *, '(a)' ) ' '
  k = 0
  do s = 0, 2
    do p = s, 0, -1
      q = s - p
      k = k + 1
      call moment ( n, x, y, p, q, nu_pq )
      write ( *, '(2x,i2,2x,i2,2x,g14.6,2x,g14.6)' ) p, q, nu_pq, nu_exact(k)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   P   Q           Alpha(P,Q)'
  write ( *, '(a)' ) '            Computed         Exact'
  write ( *, '(a)' ) ' '
  k = 0
  do s = 0, 2
    do p = s, 0, -1
      q = s - p
      k = k + 1
      call moment_normalized ( n, x,y, p, q, alpha_pq )
      write ( *, '(2x,i2,2x,i2,2x,g14.6,2x,g14.6)' ) p, q, alpha_pq, alpha_exact(k)
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   P   Q             Mu(P,Q)'
  write ( *, '(a)' ) '            Computed         Exact'
  write ( *, '(a)' ) ' '
  k = 0
  do s = 0, 2
    do p = s, 0, -1
      q = s - p
      k = k + 1
      call moment_central ( n, x, y , p, q, mu_pq )
      write ( *, '(2x,i2,2x,i2,2x,g14.6,2x,g14.6)' ) p, q, mu_pq, mu_exact(k)
    end do
  end do

  return
end
