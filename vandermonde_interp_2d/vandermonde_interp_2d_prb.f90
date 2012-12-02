program main

!*****************************************************************************80
!
!! VANDERMONDE_INTERP_2D_PRB tests VANDERMONDE_INTERP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m_test_num = 5

  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ), dimension ( m_test_num) :: m_test = (/ &
    1, 2, 3, 4, 8 /)
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'VANDERMONDE_INTERP_2D_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the VANDERMONDE_INTERP_2D library.'
  write ( *, '(a)' ) '  The QR_SOLVE library is needed.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  This test needs the TEST_INTERP_2D library.'

  call f00_num ( prob_num )
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
  write ( *, '(a)' ) 'VANDERMONDE_INTERP_2D_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( prob, m )

!*****************************************************************************80
!
!! TEST01 tests VANDERMONDE_INTERP_2D_MATRIX.
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
!    Input, integer ( kind = 4 ) M, the degree of interpolation.
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: c(:)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) tmp1
  integer ( kind = 4 ) triangle_num
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yi(:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a,i6)' ) '  Interpolate data from TEST_INTERP_2D problem #', prob
  write ( *, '(a,i6)' ) '  Create an interpolant of total degree ', m
  tmp1 = triangle_num ( m + 1 )
  write ( *, '(a,i6)' ) '  Number of data values needed is', tmp1

  nd = tmp1

  seed = 123456789

  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )

  call r8vec_uniform_01 ( nd, seed, xd )
  call r8vec_uniform_01 ( nd, seed, yd )

  allocate ( zd(1:nd) )
  call f00_f0 ( prob, nd, xd, yd, zd )

  if ( debug ) then
    call r8vec3_print ( nd, xd, yd, zd, '  X, Y, Z data:' )
  end if
!
!  Compute the Vandermonde matrix.
!
  allocate ( a(1:nd,1:nd) )
  call vandermonde_interp_2d_matrix ( nd, m, xd, yd, a )
!
!  Solve linear system.
!
  allocate ( c(1:nd) )
  call qr_solve ( nd, nd, a, zd, c )
!
!  #1:  Does interpolant match function at data points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)
  yi(1:ni) = yd(1:ni)
  allocate ( zi(1:ni) )
  call r8poly_value_2d ( m, c, ni, xi, yi, zi )

  app_error = r8vec_norm_affine ( zi - zd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 data interpolation error = ', app_error

  deallocate ( a )
  deallocate ( c )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )
  deallocate ( zd )
  deallocate ( zi )

  return
end
