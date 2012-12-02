program main

!*****************************************************************************80
!
!! LAGRANGE_INTERP_2D_TEST tests LAGRANGE_INTERP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m_test_num = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ), dimension ( m_test_num ) :: m_test = (/ &
    1, 2, 3, 4, 8 /)
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'LAGRANGE_INTERP_2D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LAGRANGE_INTERP_2D library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  The test needs the TEST_INTERP_2D library.'

  call f00_num ( prob_num )
!
!  Numerical tests.
!
  do prob = 1, prob_num
    do i = 1, m_test_num
      m = m_test(i)
      call test01 ( prob, m )
    end do
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'LAGRANGE_INTERP_2D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( prob, m )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_2D_TEST01 tests LAGRANGE_INTERP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem number.
!
!    Input, integer ( kind = 4 ) M, the polynomial degree in each dimension.
!
  implicit none

  real ( kind = 8 ) app_error
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  real ( kind = 8 ) int_error
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) my
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) prob
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xd_1d(:)
  real ( kind = 8 ), allocatable :: xi(:)
  real ( kind = 8 ), allocatable :: xi_1d(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: yd_1d(:)
  real ( kind = 8 ), allocatable :: yi(:)
  real ( kind = 8 ), allocatable :: yi_1d(:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: zdm(:)
  real ( kind = 8 ), allocatable :: zi(:)

  mx = m
  my = m

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'LAGRANGE_INTERP_2D_TEST01:'
  write ( *, '(a,i2)' ) '  Interpolate data from TEST_INTERP_2D problem #', prob
  write ( *, '(a,i2,a,i2)' ) '  Using polynomial interpolant of product degree ', mx, ' x ', my

  nd = ( mx + 1 ) * ( my + 1 )
  write ( *, '(a,i6)' ) '  Number of data points = ', nd

  allocate ( xd_1d(1:mx+1) )
  allocate ( yd_1d(1:my+1) )

  call r8vec_chebyspace ( mx + 1, 0.0D+00, 1.0D+00, xd_1d )
  call r8vec_chebyspace ( my + 1, 0.0D+00, 1.0D+00, yd_1d )

  allocate ( xd((mx+1)*(my+1)) )
  allocate ( yd((mx+1)*(my+1)) )
  allocate ( zd((mx+1)*(my+1)) )

  ij = 0
  do j = 1, my + 1
    do i = 1, mx + 1
      ij = ij + 1
      xd(ij) = xd_1d(i)
      yd(ij) = yd_1d(j)
    end do
  end do

  call f00_f0 ( prob, nd, xd, yd, zd )

  if ( nd <= 20 ) then
    call r8vec3_print ( nd, xd, yd, zd, '  X, Y, Z data:' )
  end if
!
!  #1:  Does interpolant match function at data points?
!
  ni = nd

  allocate ( xi(ni) )
  allocate ( yi(ni) )
  allocate ( zi(ni) )

  xi(1:ni) = xd(1:ni)
  yi(1:ni) = yd(1:ni)

  call lagrange_interp_2d ( mx, my, xd_1d, yd_1d, zd, ni, xi, yi, zi )

  if ( ni <= 20 ) then
    call r8vec3_print ( ni, xi, yi, zi, '  X, Y, Z interpolation:' )
  end if

  int_error = r8vec_norm_affine ( ni, zi, zd ) / real ( ni, kind = 8 )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  RMS data interpolation error = ', int_error

  deallocate ( xi )
  deallocate ( yi )
  deallocate ( zi )
!
!  #2:  Does interpolant approximate data at midpoints?
!
  if ( 1 < nd ) then

    allocate ( xi_1d(1:mx) )
    allocate ( yi_1d(1:my) )

    xi_1d(1:mx) = 0.5D+00 * ( xd_1d(1:mx) + xd_1d(2:mx+1) )
    yi_1d(1:my) = 0.5D+00 * ( yd_1d(1:my) + yd_1d(2:my+1) )

    ni = mx * my

    allocate ( xi(ni) )
    allocate ( yi(ni) )
    allocate ( zi(ni) )
    
    ij = 0
    do j = 1, my
      do i = 1, mx
        ij = ij + 1
        xi(ij) = xi_1d(i)
        yi(ij) = yi_1d(j)
      end do
    end do

    allocate ( zdm(ni) )

    call f00_f0 ( prob, ni, xi, yi, zdm )

    call lagrange_interp_2d ( mx, my, xd_1d, yd_1d, zd, ni, xi, yi, zi )

    app_error = r8vec_norm_affine ( ni, zi, zdm ) / real ( ni, kind = 8 )

    write ( *, '(a)' ) ''
    write ( *, '(a,g14.6)' ) '  RMS data approximation error = ', app_error

    deallocate ( xi )
    deallocate ( xi_1d )
    deallocate ( yi )
    deallocate ( yi_1d )
    deallocate ( zdm )
    deallocate ( zi )

  end if

  deallocate ( xd )
  deallocate ( xd_1d )
  deallocate ( yd )
  deallocate ( yd_1d )
  deallocate ( zd )

  return
end
