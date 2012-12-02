program main

!*****************************************************************************80
!
!! LAGRANGE_APPROX_1D_TEST tests LAGRANGE_APPROX_1D.
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

  integer ( kind = 4 ), parameter :: m_test_num = 7
  integer ( kind = 4 ), parameter :: nd_test_num = 3

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ), dimension ( m_test_num) :: m_test = (/ &
    0, 1, 2, 3, 4, 8, 16 /)
  integer ( kind = 4 ) nd
  integer ( kind = 4 ), dimension ( nd_test_num) :: nd_test = (/ &
    16, 64, 1000 /)
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'LAGRANGE_APPROX_1D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LAGRANGE_APPROX_1D library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  The QR_SOLVE library is needed.'
  write ( *, '(a)' ) '  These tests need the TEST_INTERP_1D library.'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num
    do j = 1, m_test_num
      m = m_test(j)
      do k = 1, nd_test_num
        nd = nd_test(k)
        call test02 ( prob, m, nd )
      end do
    end do
  end do

  do prob = 1, prob_num
    do j = 1, m_test_num
      m = m_test(j)
      do k = 1, nd_test_num
        nd = nd_test(k)
        call test03 ( prob, m, nd )
      end do
    end do
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'LAGRANGE_APPROX_1D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test02 ( prob, m, nd )

!*****************************************************************************80
!
!! TEST02 tests LAGRANGE_APPROX_1D with evenly spaced data
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) M, the polynomial approximant degree.
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) int_error
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a,i4)' ) &
    '  Approximate evenly spaced data from TEST_INTERP_1D problem #', prob
  write ( *, '(a,i4)' ) '  Use polynomial approximant of degree ', m
  write ( *, '(a,i4)' ) '  Number of data points = ', nd

  a = 0.0D+00
  b = 1.0D+00
  allocate ( xd(1:nd) )
  call r8vec_linspace ( nd, a, b, xd )

  allocate ( yd(1:nd) )
  call p00_f ( prob, nd, xd, yd )

  if ( nd < 10 ) then
    call r8vec2_print ( nd, xd, yd, '  Data array:' )
  end if
!
!  #1:  Does approximant come close to function at data points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)
  call lagrange_approx_1d ( m, nd, xd, yd, ni, xi, yi )

  int_error = r8vec_norm_affine ( nd, yi, yd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 approximation error averaged per data node = ', int_error

  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )

  return
end
subroutine test03 ( prob, m, nd )

!*****************************************************************************80
!
!!  TEST03 tests LAGRANGE_APPROX_1D with Chebyshev spaced data.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) M, the polynomial approximant degree.
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) int_error
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a,i4)' ) '  Approximate Chebyshev-spaced data from TEST_INTERP_1D problem #', prob
  write ( *, '(a,i4)' ) '  Use polynomial approximant of degree ', m
  write ( *, '(a,i4)' ) '  Number of data points = ', nd

  a = 0.0D+00
  b = 1.0D+00
  allocate ( xd(1:nd) )
  call r8vec_chebyspace ( nd, a, b, xd )

  allocate ( yd(1:nd) )
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
  call lagrange_approx_1d ( m, nd, xd, yd, ni, xi, yi )

  int_error = r8vec_norm_affine ( nd, yi, yd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 approximation error averaged per data node = ', int_error

  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )

  return
end
