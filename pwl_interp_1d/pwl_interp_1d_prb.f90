program main

!*****************************************************************************80
!
!! PWL_INTERP_1D_TEST tests PWL_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PWL_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PWL_INTERP_1D library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  The test needs the TEST_INTERP library.'

  call p00_prob_num ( prob_num )
  do prob = 1, prob_num
    call pwl_interp_1d_test01 ( prob )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PWL_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine pwl_interp_1d_test01 ( prob )

!*****************************************************************************80
!
!! PWL_INTERP_1D_TEST01 tests PWL_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
  implicit none

  real ( kind = 8 ) int_error
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ), allocatable :: xy(:,:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PWL_INTERP_1D_TEST01:'
  write ( *, '(a)' ) '  PWL_INTERP_1D evaluates the piecewise linear interpolant.'
  write ( *, '(a,i2)' ) '  Interpolate data from TEST_INTERP problem #', prob

  call p00_data_num ( prob, nd )
  write ( *, '(a,i4)' ) '  Number of data points = ', nd

  allocate ( xy(1:2,1:nd) )

  call p00_data ( prob, 2, nd, xy )
  
  call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )

  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  
  xd(1:nd) = xy(1,1:nd)
  yd(1:nd) = xy(2,1:nd)
!
!  #1:  Does interpolant match function at interpolation points?
!
  ni = nd

  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )

  xi(1:ni) = xd(1:nd)
  call pwl_interp_1d ( nd, xd, yd, ni, xi, yi )

  int_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( xd )
  deallocate ( xi )
  deallocate ( xy )
  deallocate ( yd )
  deallocate ( yi )

  return
end
