program main

!*****************************************************************************80
!
!! NEAREST_INTERP_1D_TEST tests NEAREST_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NEAREST_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the NEAREST_INTERP_1D library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  The test needs the TEST_INTERP library.'

  call p00_prob_num ( prob_num )

  ni = 11
  do prob = 1, prob_num
    call nearest_interp_1d_test01 ( prob, ni )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NEAREST_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine nearest_interp_1d_test01 ( prob, ni )

!*****************************************************************************80
!
!! NEAREST_INTERP_1D_TEST01 tests NEAREST_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the index of the problem.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
  implicit none

  real ( kind = 8 ), allocatable :: d(:,:)
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) prob
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ) xd_max
  real ( kind = 8 ) xd_min
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yi(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NEAREST_INTERP_1D_TEST01'
  write ( *, '(a,i2)' ) &
    '  Sample the nearest neighbor interpolant for problem # ', prob

  call p00_data_num ( prob, nd )

  allocate ( d(1:2,1:nd) )
  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )

  call p00_data ( prob, 2, nd, d )

  xd(1:nd) = d(1,1:nd)
  yd(1:nd) = d(2,1:nd)

  xd_min = minval ( xd(1:nd) )
  xd_max = maxval ( xd(1:nd) )

  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )

  call r8vec_linspace ( ni, xd_min, xd_max, xi )
  call nearest_interp_1d ( nd, xd, yd, ni, xi, yi )

  write ( title, '(a,i2)' ) 'X, Y for problem ', prob

  call r8vec2_print ( ni, xi, yi, title )

  deallocate ( d )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )

  return
end

