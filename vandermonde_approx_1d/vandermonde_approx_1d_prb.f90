program main

!*****************************************************************************80
!
!! VANDERMONDE_APPROX_1D_TEST tests VANDERMONDE_APPROX_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m_test_num = 8

  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ), dimension(m_test_num) :: m_test = (/ &
    0, 1, 2, 3, 4, 5, 9, 12 /)
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'VANDERMONDE_APPROX_1D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the VANDERMONDE_APPROX_1D library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  The QR_SOLVE library is needed.'
  write ( *, '(a)' ) '  The test needs the CONDITION library.'
  write ( *, '(a)' ) '  The test needs the TEST_INTERP libary.'

  call p00_prob_num ( prob_num )
  do prob = 1, prob_num
    do j = 1, m_test_num
      m = m_test(j)
      call test01 ( prob, m )
    end do
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'VANDERMONDE_APPROX_1D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( prob, m )

!*****************************************************************************80
!
!! TEST01 tests VANDERMONDE_APPROX_1D_MATRIX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem number.
!
!    Input, integer ( kind = 4 ) M, the polynomial degree.
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: c(:)
  logical, parameter :: debug = .false.
  real ( kind = 8 ) ld
  real ( kind = 8 ) li
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
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
  write ( *, '(a,i4)' ) '  Approximate data from TEST_INTERP problem #', prob

  call p00_data_num ( prob, nd )
  write ( *, '(a,i4)' ) '  Number of data points = ', nd

  allocate ( xy(1:2,1:nd) )
  call p00_data ( prob, 2, nd, xy )
  
  if ( debug ) then
    call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )
  end if

  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  xd = xy(1,1:nd)
  yd = xy(2,1:nd)
!
!  Compute the Vandermonde matrix.
!
  write ( *, '(a,i4)' ) '  Using polynomial approximant of degree ', m

  allocate ( a(1:nd,0:m) )
  call vandermonde_approx_1d_matrix ( nd, m, xd, a )
!
!  Solve linear system.
!
  allocate ( c(0:m) )
  call qr_solve ( nd, m + 1, a, yd, c )
!
!  #1:  Does approximant match function at data points?
!
  ni = nd
  allocate ( xi(1:ni) )
  xi(1:ni) = xd(1:ni)

  allocate ( yi(1:ni) )
  call r8poly_value ( m, c, ni, xi, yi )

  app_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 data approximation error = ', app_error

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
  call r8vec_linspace ( ni, xmin, xmax, xi )

  allocate ( yi(1:ni) )
  call r8poly_value ( m, c, ni, xi, yi )

  ld = sum ( sqrt ( ( ( xd(2:nd) - xd(1:nd-1) ) / ( xmax - xmin ) ) ** 2 &
                  + ( ( yd(2:nd) - yd(1:nd-1) ) / ( ymax - ymin ) ) ** 2 ) )

  li = sum ( sqrt ( ( ( xi(2:ni) - xi(1:ni-1) ) / ( xmax - xmin ) ) ** 2 &
                  + ( ( yi(2:ni) - yi(1:ni-1) ) / ( ymax - ymin ) ) ** 2 ) )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Normalized length of piecewise linear interpolant = ', ld
  write ( *, '(a,g14.6)' ) '  Normalized length of polynomial approximant       = ', li

  deallocate ( a )
  deallocate ( c )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( xy )
  deallocate ( yd )
  deallocate ( yi )

  return
end
