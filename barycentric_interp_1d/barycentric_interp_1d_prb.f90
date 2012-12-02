program main

!*****************************************************************************80
!
!! BARYCENTRIC_INTERP_1D_TEST tests BARYCENTRIC_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nd_test_num = 6

  integer ( kind = 4 ) i
  integer ( kind = 4 ) nd
  integer ( kind = 4 ), dimension ( nd_test_num ) :: nd_test = (/ &
    4, 8, 16, 32, 64, 1000 /)
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BARYCENTRIC_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BARYCENTRIC_INTERP_1D library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  The tests need the TEST_INTERP_1D library.'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num
    do i = 1, nd_test_num
      nd = nd_test(i)
      call test01 ( prob, nd )
    end do
  end do

  do prob = 1, prob_num
    do i = 1, nd_test_num
      nd = nd_test(i)
      call test02 ( prob, nd )
    end do
  end do

  do prob = 1, prob_num
    do i = 1, nd_test_num
      nd = nd_test(i)
      call test03 ( prob, nd )
    end do
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BARYCENTRIC_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( prob, nd )

!*****************************************************************************80
!
!! BARYCENTRIC_INTERP_1D_TEST01 tests LAGCHEBY1_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) ND, the number of data points to use.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) int_error
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BARYCENTRIC_INTERP_1D_TEST01:'
  write ( *, '(a,i6)' ) '  Interpolate data from TEST_INTERP_1D problem #', prob
  write ( *, '(a)' ) '  Use Chebyshev Type 1 spacing for data points.'
  write ( *, '(a,i6)' ) '  Number of data points = ', nd
!
!  Define the data.
!
  a =  0.0D+00
  b = +1.0D+00
  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  call r8vec_cheby1space ( nd, a, b, xd )
  call p00_f ( prob, nd, xd, yd )

  if ( nd < 10 ) then
    call r8vec2_print ( nd, xd, yd, '  Data array:' )
  end if
!
!  #1:  Does the interpolant match the function at the interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)
  call lagcheby1_interp_1d ( nd, xd, yd, ni, xi, yi )

  int_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )

  return
end
subroutine test02 ( prob, nd )

!*****************************************************************************80
!
!! BARYCENTRIC_INTERP_1D_TEST02 tests LAGCHEBY2_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) ND, the number of data points to use.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) int_error
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BARYCENTRIC_INTERP_1D_TEST02:'
  write ( *, '(a,i6)' ) '  Interpolate data from TEST_INTERP_1D problem #', prob
  write ( *, '(a)' ) '  Use Chebyshev Type 2 spacing for data points.'
  write ( *, '(a,i6)' ) '  Number of data points = ', nd
!
!  Define the data.
!
  a =  0.0D+00
  b = +1.0D+00
  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  call r8vec_cheby2space ( nd, a, b, xd )
  call p00_f ( prob, nd, xd, yd )

  if ( nd < 10 ) then
    call r8vec2_print ( nd, xd, yd, '  Data array:' )
  end if
!
!  #1:  Does the interpolant match the function at the interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)
  call lagcheby2_interp_1d ( nd, xd, yd, ni, xi, yi )

  int_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )

  return
end
subroutine test03 ( prob, nd )

!*****************************************************************************80
!
!! BARYCENTRIC_INTERP_1D_TEST03 tests LAGEVEN_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) ND, the number of data points to use.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) int_error
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BARYCENTRIC_INTERP_1D_TEST03:'
  write ( *, '(a,i6)' ) '  Interpolate data from TEST_INTERP_1D problem #', prob
  write ( *, '(a)' ) '  Use even spacing for data points.'
  write ( *, '(a,i6)' ) '  Number of data points = ', nd
!
!  Define the data.
!
  a =  0.0D+00
  b = +1.0D+00
  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  call r8vec_midspace ( nd, a, b, xd )
  call p00_f ( prob, nd, xd, yd )

  if ( nd < 10 ) then
    call r8vec2_print ( nd, xd, yd, '  Data array:' )
  end if
!
!  #1:  Does the interpolant match the function at the interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)
  call lageven_interp_1d ( nd, xd, yd, ni, xi, yi )

  int_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )

  return
end
