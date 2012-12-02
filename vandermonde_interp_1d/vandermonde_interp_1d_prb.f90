program main

!*****************************************************************************80
!
!! VANDERMONDE_INTERP_1D_TEST tests VANDERMONDE_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2012
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
  write ( *, '(a)' ) 'VANDERMONDE_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the VANDERMONDE_INTERP_1D library.'
  write ( *, '(a)' ) '  The QR_SOLVE library is needed.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  This test needs the CONDITION library.'
  write ( *, '(a)' ) '  This test needs the TEST_INTERP library.'

  call p00_prob_num ( prob_num )
  do prob = 1, prob_num
    call test01 ( prob )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'VANDERMONDE_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( prob )

!*****************************************************************************80
!
!! TEST01 tests VANDERMONDE_INTERP_1D_MATRIX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: c(:)
  real ( kind = 8 ) condition
  logical, parameter :: debug = .false.
  real ( kind = 8 ) int_error
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
  write ( *, '(a,i2)' ) '  Interpolate data from TEST_INTERP problem #', prob

  call p00_data_num ( prob, nd )
  write ( *, '(a,i2)' ) '  Number of data points = ', nd

  allocate ( xy(1:2,1:nd) )

  call p00_data ( prob, 2, nd, xy )
  
  if ( debug ) then
    call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )
  end if

  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )

  xd(1:nd) = xy(1,1:nd)
  yd(1:nd) = xy(2,1:nd)
!
!  Choose the degree of the polynomial to be ND - 1.
!
  m = nd - 1
!
!  Compute Vandermonde matrix and get condition number.
!
  allocate ( a(1:nd,1:nd) )

  call vandermonde_interp_1d_matrix ( nd, xd, a )

  call condition_hager ( nd, a, condition )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Condition of Vandermonde matrix is ', condition
!
!  Solve linear system.
!
  allocate ( c(1:nd) )

  call qr_solve ( nd, nd, a, yd, c )
!
!  #1:  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)
  call r8poly_value ( m, c, ni, xi, yi )

  int_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = 8 )

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
  xmin = minval ( xd(1:nd) )
  xmax = maxval ( xd(1:nd) )
  ymin = minval ( yd(1:nd) )
  ymax = maxval ( yd(1:nd) )

  ni = 501
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  call r8vec_linspace ( ni, xmin, xmax, xi )
  call r8poly_value ( m, c, ni, xi, yi )

  ld = sum ( sqrt ( ( ( xd(2:nd) - xd(1:nd-1) ) / ( xmax - xmin ) )**2 &
                  + ( ( yd(2:nd) - yd(1:nd-1) ) / ( ymax - ymin ) )**2 ) )

  li = sum ( sqrt ( ( ( xi(2:ni) - xi(1:ni-1) ) / ( xmax - xmin ) )**2 &
                  + ( ( yi(2:ni) - yi(1:ni-1) ) / ( ymax - ymin ) )**2 ) )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Normalized length of piecewise linear interpolant = ', ld
  write ( *, '(a,g14.6)' ) '  Normalized length of polynomial interpolant       = ', li

  deallocate ( a )
  deallocate ( c )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( xy )
  deallocate ( yd )
  deallocate ( yi )

  return
end
