program main

!*****************************************************************************80
!
!! LAGRANGE_INTERP_1D_TEST tests LAGRANGE_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nd_test_num = 6

  integer ( kind = 4 ) j
  integer ( kind = 4 ) nd
  integer ( kind = 4 ), dimension ( nd_test_num ) :: nd_test = (/ &
    4, 8, 16, 32, 64, 256 /)
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'LAGRANGE_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LAGRANGE_INTERP_1D library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  The tests needs the TEST_INTERP_1D library.'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num
    do j = 1, nd_test_num
      nd = nd_test(j)
      call test02 ( prob, nd )
    end do
  end do

  do prob = 1, prob_num
    do j = 1, nd_test_num
      nd = nd_test(j)
      call test03 ( prob, nd )
    end do
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'LAGRANGE_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test02 ( prob, nd )

!*****************************************************************************80
!
!! TEST02 tests LAGRANGE_VALUE_1D with evenly spaced data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2012
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

  integer ( kind = 4 ) nd

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) int_error
  real ( kind = 8 ) ld
  real ( kind = 8 ) li
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ), allocatable :: yi(:)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a,i4)' ) '  Interpolate data from TEST_INTERP_1D problem #', prob
  write ( *, '(a)' ) '  Use even spacing for data points.'
  write ( *, '(a,i4)' ) '  Number of data points = ', nd

  a = 0.0D+00
  b = 1.0D+00
  
  call r8vec_linspace ( nd, a, b, xd )

  call p00_f ( prob, nd, xd, yd )

  if ( nd < 10 ) then
    call r8vec2_print ( nd, xd, yd, '  Data array:' )
  end if
!
!  #1:  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)

  call lagrange_value_1d ( nd, xd, yd, ni, xi, yi )

  int_error = r8vec_norm_affine ( nd, yi, yd ) / real ( ni, kind = 8 )

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
  ymin = minval ( yd(1:nd) )
  ymax = maxval ( yd(1:nd) )

  ni = 501
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  call r8vec_linspace ( ni, a, b, xi )
  call lagrange_value_1d ( nd, xd, yd, ni, xi, yi )

  ld = 0.0D+00
  do i = 1, nd - 1
    ld = ld + sqrt ( ( ( xd(i+1) - xd(i) ) / ( b - a ) )**2 &
                   + ( ( yd(i+1) - yd(i) ) / ( ymax - ymin ) )**2 ) 
  end do

  li = 0.0D+00
  do i = 1, ni - 1
    li = li + sqrt ( ( ( xi(i+1) - xi(i) ) / ( b - a ) )**2 &
                   + ( ( yi(i+1) - yi(i) ) / ( ymax - ymin ) )**2 )
  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  Normalized length of piecewise linear interpolant = ', ld
  write ( *, '(a,g14.6)' ) &
    '  Normalized length of polynomial interpolant       = ', li

  deallocate ( xi )
  deallocate ( yi )

  return
end
subroutine test03 ( prob, nd )

!*****************************************************************************80
!
!! TEST03 tests LAGRANGE_VALUE_1D with Chebyshev spaced data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2012
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

  integer ( kind = 4 ) nd

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) int_error
  real ( kind = 8 ) ld
  real ( kind = 8 ) li
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ), allocatable :: yi(:)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a,i4)' ) '  Interpolate data from TEST_INTERP_1D problem #', prob
  write ( *, '(a)' ) '  Use Chebyshev spacing for data points.'
  write ( *, '(a,i4)' ) '  Number of data points = ', nd

  a = 0.0D+00
  b = 1.0D+00
  call r8vec_chebyspace ( nd, a, b, xd )

  call p00_f ( prob, nd, xd, yd )

  if ( nd < 10 ) then
    call r8vec2_print ( nd, xd, yd, '  Data array:' )
  end if
!
!  #1:  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)

  call lagrange_value_1d ( nd, xd, yd, ni, xi, yi )

  int_error = r8vec_norm_affine ( nd, yi, yd ) / real ( ni, kind = 8 )

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
  ymin = minval ( yd(1:nd) )
  ymax = maxval ( yd(1:nd) )

  ni = 501
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  call r8vec_linspace ( ni, a, b, xi )
  call lagrange_value_1d ( nd, xd, yd, ni, xi, yi )

  ld = 0.0D+00
  do i = 1, nd - 1
    ld = ld + sqrt ( ( ( xd(i+1) - xd(i) ) / ( b - a ) )**2 &
                   + ( ( yd(i+1) - yd(i) ) / ( ymax - ymin ) )**2 ) 
  end do

  li = 0.0D+00
  do i = 1, ni - 1
    li = li + sqrt ( ( ( xi(i+1) - xi(i) ) / ( b - a ) )**2 &
                   + ( ( yi(i+1) - yi(i) ) / ( ymax - ymin ) )**2 )
  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  Normalized length of piecewise linear interpolant = ', ld
  write ( *, '(a,g14.6)' ) &
    '  Normalized length of polynomial interpolant       = ', li

  deallocate ( xi )
  deallocate ( yi )

  return
end
