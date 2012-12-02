program main

!*****************************************************************************80
!
!! VANDERMONDE_APPROX_2D_TEST tests VANDERMONDE_APPROX_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m_test_num = 5

  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ), dimension(m_test_num) :: m_test = (/ 0, 1, 2, 4, 8 /)
  integer ( kind = 4 ) grid
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'VANDERMONDE_APPROX_2D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the VANDERMONDE_APPROX_2D library.'
  write ( *, '(a)' ) '  This test also needs the TEST_INTERP_2D library.'

  call f00_num ( prob_num )

  do prob = 1, prob_num
    grid = 1
    do j = 1, m_test_num
      m = m_test(j)
      call test01 ( prob, grid, m )
    end do
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'VANDERMONDE_APPROX_2D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( prob, grd, m )

!*****************************************************************************80
!
!! VANDERMONDE_APPROX_2D_TEST01 tests VANDERMONDE_APPROX_2D_MATRIX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem number.
!
!    Input, integer ( kind = 4 ) GRD, the grid number.
!    (Can't use GRID as the name because that's also a plotting function.)
!
!    Input, integer ( kind = 4 ) M, the total polynomial degree.
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) grd
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) tm
  integer ( kind = 4 ) triangle_num
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yi(:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a,i4)' ) '  Approximate data from TEST_INTERP_2D problem #', prob
  write ( *, '(a,i4)' ) '  Use grid from TEST_INTERP_2D with index #', grd
  write ( *, '(a,i4)' ) '  Using polynomial approximant of total degree ', m

  call g00_size ( grd, nd )
  write ( *, '(a,i6)' ) '  Number of data points = ', nd

  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  call g00_xy ( grd, nd, xd, yd )

  allocate ( zd(1:nd) )
  call f00_f0 ( prob, nd, xd, yd, zd )

  if ( nd < 10 ) then
    call r8vec3_print ( nd, xd, yd, zd, '  X, Y, Z data:' )
  end if
!
!  Compute the Vandermonde matrix.
!
  tm = triangle_num ( m + 1 );
  allocate ( a(1:nd,1:tm) )
  call vandermonde_approx_2d_matrix ( nd, m, tm, xd, yd, a )
!
!  Solve linear system.
!
  allocate ( c(1:tm) )
  call qr_solve ( nd, tm, a, zd, c )
!
!  #1:  Does approximant match function at data points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)
  yi(1:ni) = yd(1:ni)

  allocate ( zi(1:ni) )
  call r8poly_value_2d ( m, c, ni, xi, yi, zi )

  app_error = r8vec_norm_affine ( ni, zi, zd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 data approximation error = ', app_error

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
