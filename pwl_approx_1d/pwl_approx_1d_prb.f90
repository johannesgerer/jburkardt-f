program main

!*****************************************************************************80
!
!! PWL_APPROX_1D_TEST tests PWL_APPROX_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nc_test_num = 4
  integer ( kind = 4 ), parameter :: nd_test_num = 2

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nc
  integer ( kind = 4 ), dimension(nc_test_num) :: nc_test = (/ 2, 4, 8, 16 /)
  integer ( kind = 4 ) nd
  integer ( kind = 4 ), dimension(nd_test_num) :: nd_test = (/ 16, 64 /)
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PWL_APPROX_1D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PWL_APPROX_1D library.'
  write ( *, '(a)' ) '  The QR_SOLVE library is needed.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  The test also needs the TEST_INTERP_1D library.'
  call p00_prob_num ( prob_num )

  do prob = 1, prob_num
    do j = 1, nc_test_num
      nc = nc_test(j)
      do k = 1, nd_test_num
        nd = nd_test(k)
        call test01 ( prob, nc, nd )
      end do
    end do
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PWL_APPROX_1D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( prob, nc, nd )

!*****************************************************************************80
!
!! TEST01 tests PWL_APPROX_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) NC, the number of control points.
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
  implicit none

  real ( kind = 8 ) app_error
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ), allocatable :: xc(:)
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ), allocatable :: yc(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yi(:)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a,i4)' ) '  Approximate data from TEST_INTERP_1D problem #', prob
  write ( *, '(a,i4)' ) '  Number of control points = ', nc
  write ( *, '(a,i4)' ) '  Number of data points = ', nd

  xmin = 0.0D+00
  xmax = 1.0D+00

  allocate ( xd(1:nd) )
  call r8vec_linspace ( nd, xmin, xmax, xd )

  allocate ( yd(1:nd) )
  call p00_f ( prob, nd, xd, yd )

  if ( nd < 10 ) then
    call r8vec2_print ( nd, xd, yd, '  Data array:' )
  end if
!
!  Determine control values.
!
  allocate ( xc(1:nc) )
  call r8vec_linspace ( nc, xmin, xmax, xc )

  allocate ( yc(1:nc) )
  call pwl_approx_1d ( nd, xd, yd, nc, xc, yc )
!
!  #1:  Does approximant come close to function at data points?
!
  ni = nd
  allocate ( xi(1:ni) )
  xi(1:ni) = xd(1:ni)

  allocate ( yi(1:ni) )
  call pwl_interp_1d ( nc, xc, yc, ni, xi, yi )

  app_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 approximation error averaged per data node = ', app_error

  deallocate ( xc )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yc )
  deallocate ( yd )
  deallocate ( yi )

  return
end
