program main

!*****************************************************************************80
!
!! SHEPARD_INTERP_2D_PRB tests SHEPARD_INTERP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) p_test_num
  parameter ( p_test_num = 4 )

  integer ( kind = 4 ) g
  integer ( kind = 4 ) j
  real ( kind = 8 ) p
  real ( kind = 8 ), dimension ( p_test_num ) :: p_test = (/ &
    1.0D+00, 2.0D+00, 4.0D+00, 8.0D+00 /)
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SHEPARD_INTERP_2D_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SHEPARD_INTERP_2D library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  This test also needs the TEST_INTERP_2D library.'

  call f00_num ( prob_num )
  g = 1

  do prob = 1, prob_num

    do j = 1, p_test_num
      p = p_test(j)
      call test01 ( prob, g, p )
    end do
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SHEPARD_INTERP_2D_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( prob, g, p )

!*****************************************************************************80
!
!! TEST01 tests SHEPARD_INTERP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem number.
!
!    Input, integer ( kind = 4 ) G, the grid number.
!
!    Input, real ( kind = 8 ) P, the power used in the distance weighting.
!
  implicit none

  logical, parameter :: debug = .false.
  integer ( kind = 4 ) g
  real ( kind = 8 ) int_error
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  real ( kind = 8 ) p
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yi(:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a,i6)' ) '  Interpolate data from TEST_INTERP_2D problem #', prob
  write ( *, '(a,i6)' ) '  using grid #', g
  write ( *, '(a,g14.6)' ) '  using Shepard interpolation with P = ', p

  call g00_size ( g, nd )
  write ( *, '(a,i6)' ) '  Number of data points = ', nd

  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  call g00_xy ( g, nd, xd, yd )
  
  allocate ( zd(1:nd) )
  call f00_f0 ( prob, nd, xd, yd, zd )

  if ( debug ) then
    call r8vec3_print ( nd, xd, yd, zd, '  X, Y, Z data:' )
  end if
!
!  #1:  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)
  yi(1:ni) = yd(1:ni)

  allocate ( zi(1:ni) )
  call shepard_interp_2d ( nd, xd, yd, zd, p, ni, xi, yi, zi )

  int_error = r8vec_norm_affine ( nd, zi, zd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', &
    int_error

  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )
  deallocate ( zd )
  deallocate ( zi )

  return
end
