program main

!*****************************************************************************80
!
!! RBF_INTERP_1D_TEST tests RBF_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) nd
  external phi1
  external phi2
  external phi3
  external phi4
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ) r0
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ), allocatable :: xy(:,:)

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RBF_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the RBF_INTERP_1D library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  The test needs the TEST_INTERP library.'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num
!
!  Determine an appropriate value of R0, the spacing parameter.
!
    call p00_data_num ( prob, nd )
    allocate ( xy(1:2,1:nd) )
    call p00_data ( prob, 2, nd, xy )
    allocate ( xd(1:nd) )
    xd(1:nd) = xy(1,1:nd)
    xmax = maxval ( xd )
    xmin = minval ( xd )
    r0 = ( xmax - xmin ) / real ( nd - 1, kind = 8 )
    deallocate ( xd )
    deallocate ( xy )

    call test01 ( prob, phi1, 'phi1', r0 )
    call test01 ( prob, phi2, 'phi2', r0 )
    call test01 ( prob, phi3, 'phi3', r0 )
    call test01 ( prob, phi4, 'phi4', r0 )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RBF_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( prob, phi, phi_name, r0 )

!*****************************************************************************80
!
!! TEST01 tests RBF_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the index of the problem.
!
!    Input, external PHI, the name of the radial basis function.
!
!    Input, character ( len = * ) PHI_NAME, the name of the radial basis function.
!
!    Input, real ( kind = 8 ) R0, the scale factor.  Typically, this might be
!    a small multiple of the average distance between points.
!
  implicit none

  real ( kind = 8 ) int_error
  real ( kind = 8 ) ld
  real ( kind = 8 ) li
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  external phi
  character ( len = * ) phi_name
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r0
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ), allocatable :: xy(:,:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yi(:)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a,i4)' ) '  Interpolate data from TEST_INTERP problem #', prob
  write ( *, '(a)' ) '  using radial basis function "' // trim ( phi_name ) // '".'
  write ( *, '(a,g14.6)' ) '  Scale factor R0 = ', r0

  call p00_data_num ( prob, nd )
  write ( *, '(a,i4)' ) '  Number of data points = ', nd

  allocate ( xy(1:2,1:nd) )
  call p00_data ( prob, 2, nd, xy )
  
  call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )

  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  xd(1:nd) = xy(1,1:nd)
  yd(1:nd) = xy(2,1:nd)

  m = 1
  allocate ( w(1:nd) )
  call rbf_weight ( m, nd, xd, r0, phi, yd, w )
!
!  #1:  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:nd) = xd(1:nd)
  call rbf_interp ( m, nd, xd, r0, phi, w, ni, xi, yi )

  int_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( xi )
  deallocate ( yi )
!
!  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
!  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
!  (YMAX-YMIN).
!
  xmin = minval ( xd(1:nd) )
  xmax = maxval ( xd(1:nd) )
  ymin = minval ( yd(1:nd) )
  ymax = maxval ( yd(1:nd) )

  ni = 501
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  call r8vec_linspace ( ni, xmin, xmax, xi )
  call rbf_interp ( m, nd, xd, r0, phi, w, ni, xi, yi )

  ld = sum ( sqrt ( ( ( xd(2:nd) - xd(1:nd-1) ) / ( xmax - xmin ) ) ** 2 &
                  + ( ( yd(2:nd) - yd(1:nd-1) ) / ( ymax - ymin ) ) ** 2 ) )

  li = sum ( sqrt ( ( ( xi(2:ni) - xi(1:ni-1) ) / ( xmax - xmin ) ) ** 2 &
                  + ( ( yi(2:ni) - yi(1:ni-1) ) / ( ymax - ymin ) ) ** 2 ) )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Normalized length of piecewise linear interpolant = ', ld
  write ( *, '(a,g14.6)' ) '  Normalized length of polynomial interpolant       = ', li

  deallocate ( w )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( xy )
  deallocate ( yd )
  deallocate ( yi )

  return
end
