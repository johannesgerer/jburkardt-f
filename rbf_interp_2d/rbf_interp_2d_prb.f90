program main

!*****************************************************************************80
!
!! RBF_INTERP_2D_TEST tests RBF_INTERP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) g
  external phi1
  external phi2
  external phi3
  external phi4
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RBF_INTERP_2D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the RBF_INTERP_2D library.'
  write ( *, '(a)' ) '  The R8LIB library is required.'
  write ( *, '(a)' ) '  This test also needs the TEST_INTERP_2D library.'

  call f00_num ( prob_num )
  g = 1

  do prob = 1, prob_num
    call test01 ( prob, g, phi1, 'phi1' )
    call test01 ( prob, g, phi2, 'phi2' )
    call test01 ( prob, g, phi3, 'phi3' )
    call test01 ( prob, g, phi4, 'phi4' )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RBF_INTERP_2D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( prob, g, phi, phi_name )

!*****************************************************************************80
!
!! RBF_INTERP_2D_TEST01 tests RBF_INTERP_2D.
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
!    Input, integer ( kind = 4 ) PROB, the index of the problem.
!
!    Input, integer ( kind = 4 ) G, the index of the grid.
!
!    Input, external PHI, the radial basis function.
!
!    Input, character ( len = * ) PHI_NAME, the name of the radial basis function.
!
  implicit none

  logical, parameter :: debug = .false.
  real ( kind = 8 ) e
  integer ( kind = 4 ) g
  integer ( kind = 4 ) i
  real ( kind = 8 ) int_error
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  external phi
  integer ( kind = 4 ) prob
  character ( len = * ) phi_name
  real ( kind = 8 ) r0
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ) volume
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ), allocatable :: xyd(:,:)
  real ( kind = 8 ), allocatable :: xyi(:,:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RBF_INTERP_2D_TEST01:'
  write ( *, '(a,i6)' ) '  Interpolate data from TEST_INTERP_2D problem #', prob
  write ( *, '(a,i6)' ) '  using grid #', g
  write ( *, '(a)' ) '  using radial basis function "' // trim ( phi_name ) &
    // '".'

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

  m = 2
  allocate ( xyd(2,nd) )

  do i = 1, nd
    xyd(1,i) = xd(i)
    xyd(2,i) = yd(i)
  end do

  xmax = maxval ( xd(1:nd) )
  xmin = minval ( xd(1:nd) )
  ymax = maxval ( yd(1:nd) )
  ymin = minval ( yd(1:nd) )
  volume = ( xmax - xmin ) * ( ymax - ymin )

  e = 1.0D+00 / real ( m, kind = 8 )
  r0 = ( volume / nd ) ** e

  write ( *, '(a,g14.6)' ) '  Setting R0 = ', r0

  allocate ( w(1:nd) )

  call rbf_weight ( m, nd, xyd, r0, phi, zd, w )
!
!  #1:  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xyi(1:2,1:ni) )
  xyi(1:2,1:ni) = xyd(1:2,1:ni)
  allocate ( zi(1:ni) )

  call rbf_interp ( m, nd, xyd, r0, phi, w, ni, xyi, zi )

  int_error = r8vec_norm_affine ( ni, zi, zd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( w )
  deallocate ( xd )
  deallocate ( xyd )
  deallocate ( xyi )
  deallocate ( yd )
  deallocate ( zd )
  deallocate ( zi )

  return
end
