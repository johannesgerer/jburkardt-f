program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_MAT_PRB.
!
!  Discussion:
!
!    TEST_MAT_PRB calls the TEST_MAT test routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MAT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_MAT library.'

  call test_cond ( )
  call test_determinant ( )
  call test_eigen ( )
  call test_inverse ( )
  call test_null ( )
  call test_plu ( )
  call test_solution ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MAT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test_cond ( )

!*****************************************************************************80
!
!! TEST_COND tests the condition number computations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) cond
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  character ( len = 20 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_COND'
  write ( *, '(a)' ) '  Compute the condition number of an example of each'
  write ( *, '(a)' ) '  test matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix title             N      COND'
  write ( *, '(a)' ) ' '
!
!  AEGERTER matrix.
!
  title = 'AEGERTER'
  n = 5
  call aegerter_condition ( n, cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
!
!  BAB matrix.
!
  title = 'BAB'
  seed = 123456789
  alpha = anint ( 50.0D+00 * r8_uniform_01 ( seed ) - 25.0D+00 ) / 5.0D+00
  beta = anint ( 50.0D+00 * r8_uniform_01 ( seed ) - 25.0D+00 ) / 5.0D+00
  call bab_condition ( n, alpha, beta, cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
!
!  BODEWIG matrix.
!
  title = 'BODEWIG'
  n = 4
  call bodewig_condition ( cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
!
!  COMBIN matrix.
!
  title = 'COMBIN'
  n = 3
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta ) / 5.0D+00
  call combin_condition ( alpha, beta, n, cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
!
!  CONEX3 matrix.
!
  title = 'CONEX3'
  n = 5
  call conex3_condition ( n, cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
!
!  RUTIS5 matrix.
!
  title = 'RUTIS5'
  n = 4
  call rutis5_condition ( cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
!
!  SUMMATION matrix.
!
  title = 'SUMMATION'
  n = 5
  call summation_condition ( n, cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
!
!  TRI_UPPER matrix.
!
  title = 'TRI_UPPER'
  n = 5
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call tri_upper_condition ( alpha, n, cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
!
!  WILK03 matrix.
!
  title = 'WILK03'
  n = 3
  call wilk03_condition ( cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
!
!  WILSON matrix.
!
  title = 'WILSON'
  n = 4
  call wilson_condition ( cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond

  return
end
subroutine test_determinant ( )

!*****************************************************************************80
!
!! TEST_DETERMINANT tests the determinant computations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) col_num
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) d3
  real ( kind = 8 ) d4
  real ( kind = 8 ) d5
  real ( kind = 8 ) da
  real ( kind = 8 ) determ1
  real ( kind = 8 ) determ2
  real ( kind = 8 ) di
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm_frobenius
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p
  real ( kind = 8 ) perturb
  integer ( kind = 4 ), allocatable, dimension ( : ) :: pivot
  real ( kind = 8 ) prob
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) row_num
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_save
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: u
  real ( kind = 8 ), allocatable, dimension ( : ) :: v1
  real ( kind = 8 ), allocatable, dimension ( : ) :: v2
  real ( kind = 8 ), allocatable, dimension ( : ) :: v3
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  integer ( kind = 4 ) x_n
  real ( kind = 8 ), allocatable, dimension ( : ) :: y
  integer ( kind = 4 ) y_n
  real ( kind = 8 ) y_sum
  real ( kind = 8 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DETERMINANT'
  write ( *, '(a)' ) '  Compute the determinants of an example of each'
  write ( *, '(a)' ) '  test matrix; compare with the determinant routine,'
  write ( *, '(a)' ) '  if available.  Print the matrix Frobenius norm'
  write ( *, '(a)' ) '  for an estimate of magnitude.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix title             N      ' // &
                     'Determ          Determ         ||A||'
  write ( *, '(a)' ) ' '
!
!  AEGERTER matrix.
!
  title = 'AEGERTER'
  n = 5
  allocate ( a(1:n,1:n) )
  call aegerter ( n, a )
  call aegerter_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ANTICIRCULANT matrix.
!
  title = 'ANTICIRCULANT'
  n = 3
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = ( anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) ) / 5.0D+00
  call anticirculant ( n, n, x, a )
  call anticirculant_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  ANTICIRCULANT matrix.
!
  title = 'ANTICIRCULANT'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = ( anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) ) / 5.0D+00
  call anticirculant ( n, n, x, a )
  call anticirculant_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  ANTICIRCULANT matrix.
!
  title = 'ANTICIRCULANT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = ( anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) ) / 5.0D+00
  call anticirculant ( n, n, x, a )
  call anticirculant_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  ANTIHADAMARD matrix.
!
  title = 'ANTIHADAMARD'
  n = 5
  allocate ( a(1:n,1:n) )
  call antihadamard ( n, a )
  call antihadamard_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ANTISYMM_RANDOM matrix.
!
  title = 'ANTISYMM_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  call antisymm_random ( n, seed, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  ANTISYMM_RANDOM matrix.
!
  title = 'ANTISYMM_RANDOM'
  n = 6
  allocate ( a(1:n,1:n) )
  seed = 123456789
  call antisymm_random ( n, seed, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  BAB matrix.
!
  title = 'BAB'
  n = 5
  seed = 123456789
  alpha = anint ( 50.0D+00 * r8_uniform_01 ( seed ) - 25.0D+00 ) / 5.0D+00
  beta = anint ( 50.0D+00 * r8_uniform_01 ( seed ) - 25.0D+00 ) / 5.0D+00
  allocate ( a(1:n,1:n) )
  call bab ( n, alpha, beta, a )
  call bab_determinant ( n, alpha, beta, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  BIMARKOV_RANDOM matrix.
!
  title = 'BIMARKOV_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  call bimarkov_random ( n, seed, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  BIS matrix.
!
  title = 'BIS'
  n = 5
  seed = 123456789
  alpha = anint ( 50.0D+00 * r8_uniform_01 ( seed ) ) / 5.0D+00
  beta = anint ( 50.0D+00 * r8_uniform_01 ( seed ) ) / 5.0D+00
  allocate ( a(1:n,1:n) )
  call bis ( alpha, beta, n, n, a )
  call bis_determinant ( alpha, beta, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  BODEWIG matrix.
!
  title = 'BODEWIG'
  n = 4
  allocate ( a(1:n,1:n) )
  call bodewig ( a )
  call bodewig_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  BOOTHROYD matrix.
!
  title = 'BOOTHROYD'
  n = 5
  allocate ( a(1:n,1:n) )
  call boothroyd ( n, a )
  call boothroyd_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  BORDERBAND matrix.
!
  title = 'BORDERBAND'
  n = 5
  allocate ( a(1:n,1:n) )
  call borderband ( n, a )
  call borderband_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CARRY matrix.
!
  title = 'CARRY'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  k = i4_uniform ( 2, 20, seed )
  call carry ( k, n, a )
  call carry_determinant ( k, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CAUCHY matrix.
!
  title = 'CAUCHY'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  call r8vec_uniform_01 ( n, seed, y )
  call cauchy ( n, x, y, a )
  call cauchy_determinant ( n, x, y, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  CHEBY_DIFF1 matrix.
!
  title = 'CHEBY_DIFF1'
  n = 5
  allocate ( a(1:n,1:n) )
  call cheby_diff1 ( n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  CHEBY_DIFF1 matrix.
!
  title = 'CHEBY_DIFF1'
  n = 6
  allocate ( a(1:n,1:n) )
  call cheby_diff1 ( n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  CHEBY_T matrix.
!
  title = 'CHEBY_T'
  n = 5
  allocate ( a(1:n,1:n) )
  call cheby_t ( n, a )
  call cheby_t_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CHEBY_U matrix.
!
  title = 'CHEBY_U'
  n = 5
  allocate ( a(1:n,1:n) )
  call cheby_u ( n, a )
  call cheby_u_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CHEBY_VAN1 matrix.
!
  title = 'CHEBY_VAN1'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 10.0D+00 * x(1:n) - 5.0D+00 )
  call cheby_van1 ( n, x, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  CHEBY_VAN2 matrix.
!
  do n = 2, 10
    title = 'CHEBY_VAN2'
    allocate ( a(1:n,1:n) )
    call cheby_van2 ( n, a )
    call cheby_van2_determinant ( n, determ1 )
    call r8mat_determinant ( n, a, determ2 )
    norm_frobenius = r8mat_norm_fro ( n, n, a )
    write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      title, n, determ1, determ2, norm_frobenius
    deallocate ( a )
  end do
!
!  CHEBY_VAN3 matrix.
!
  title = 'CHEBY_VAN3'
  n = 5
  allocate ( a(1:n,1:n) )
  call cheby_van3 ( n, a )
  call cheby_van3_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CHOW matrix.
!
  title = 'CHOW'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed );
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta ) / 5.0D+00
  call chow ( alpha, beta, n, n, a )
  call chow_determinant ( alpha, beta, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CIRCULANT matrix.
!
  title = 'CIRCULANT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789;
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call circulant ( n, n, x, a )
  call circulant_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  CIRCULANT2 matrix.
!
  title = 'CIRCULANT2'
  n = 3
  allocate ( a(1:n,1:n) )
  call circulant2 ( n, a )
  call circulant2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CIRCULANT2 matrix.
!
  title = 'CIRCULANT2'
  n = 4
  allocate ( a(1:n,1:n) )
  call circulant2 ( n, a )
  call circulant2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CIRCULANT2 matrix.
!
  title = 'CIRCULANT2'
  n = 5
  allocate ( a(1:n,1:n) )
  call circulant2 ( n, a )
  call circulant2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CLEMENT1 matrix.
!
  title = 'CLEMENT1'
  n = 5
  allocate ( a(1:n,1:n) )
  call clement1 ( n, a )
  call clement1_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CLEMENT1 matrix.
!
  title = 'CLEMENT1'
  n = 6
  allocate ( a(1:n,1:n) )
  call clement1 ( n, a )
  call clement1_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CLEMENT2 matrix.
!
  title = 'CLEMENT2'
  n = 5
  allocate ( a(1:n,1:n) )
  call clement2 ( n, a )
  call clement2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CLEMENT2 matrix.
!
  title = 'CLEMENT2'
  n = 6
  allocate ( a(1:n,1:n) )
  call clement2 ( n, a )
  call clement2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CLEMENT3.
!
  title = 'CLEMENT3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n-1) )
  allocate ( y(1:n-1) )
  seed = 123456789
  call r8vec_uniform_01 ( n - 1, seed, x )
  x(1:n-1) = anint ( 10.0D+00 * x(1:n-1) - 5.0D+00 )
  call r8vec_uniform_01 ( n - 1, seed, y )
  y(1:n-1) = anint ( 10.0D+00 * y(1:n-1) - 5.0D+00 )
  call clement3 ( n, x, y, a )
  call clement3_determinant ( n, x, y, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  CLEMENT3.
!
  title = 'CLEMENT3'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n-1) )
  allocate ( y(1:n-1) )
  seed = 123456789
  call r8vec_uniform_01 ( n - 1, seed, x )
  x(1:n-1) = anint ( 10.0D+00 * x(1:n-1) - 5.0D+00 )
  call r8vec_uniform_01 ( n - 1, seed, y )
  y(1:n-1) = anint ( 10.0D+00 * y(1:n-1) - 5.0D+00 )
  call clement3 ( n, x, y, a )
  call clement3_determinant ( n, x, y, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  COMBIN matrix.
!
  title = 'COMBIN'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta ) / 5.0D+00
  call combin ( alpha, beta, n, a )
  call combin_determinant ( alpha, beta, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  COMPANION.
!
  title = 'COMPANION'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 10.0D+00 * x(1:n) - 5.0D+00 )
  call companion ( n, x, a )
  call companion_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  COMPLEX_I
!
  title = 'COMPLEX_I'
  n = 2
  allocate ( a(1:n,1:n) )
  call complex_i ( a )
  call complex_i_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CONEX1 matrix.
!
  title = 'CONEX1'
  n = 4
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call conex1 ( alpha, a )
  call conex1_determinant ( alpha, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CONEX2 matrix.
!
  title = 'CONEX2'
  n = 3
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call conex2 ( alpha, a )
  call conex2_determinant ( alpha, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CONEX3 matrix.
!
  title = 'CONEX3'
  n = 5
  allocate ( a(1:n,1:n) )
  call conex3 ( n, a )
  call conex3_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CONFERENCE matrix.
!
  title = 'CONFERENCE'
  n = 6
  allocate ( a(1:n,1:n) )
  call conference ( n, a )
  call conference_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  CREATION matrix.
!
  title = 'CREATION'
  n = 5
  allocate ( a(1:n,1:n) )
  call creation ( n, n, a )
  call creation_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DAUB2 matrix.
!
  title = 'DAUB2'
  n = 4
  allocate ( a(1:n,1:n) )
  call daub2 ( n, a )
  call daub2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DAUB4 matrix.
!
  title = 'DAUB4'
  n = 8
  allocate ( a(1:n,1:n) )
  call daub4 ( n, a )
  call daub4_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DAUB6 matrix.
!
  title = 'DAUB6'
  n = 12
  allocate ( a(1:n,1:n) )
  call daub6 ( n, a )
  call daub6_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DAUB8 matrix.
!
  title = 'DAUB8'
  n = 16
  allocate ( a(1:n,1:n) )
  call daub8 ( n, a )
  call daub8_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DAUB10 matrix.
!
  title = 'DAUB10'
  n = 20
  allocate ( a(1:n,1:n) )
  call daub10 ( n, a )
  call daub10_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DAUB12 matrix.
!
  title = 'DAUB12'
  n = 24
  allocate ( a(1:n,1:n) )
  call daub12 ( n, a )
  call daub12_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DIAGONAL.
!
  title = 'DIAGONAL'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 10.0D+00 * x(1:n) - 5.0D+00 )
  call diagonal ( n, n, x, a )
  call diagonal_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  DIF1 matrix.
!
  title = 'DIF1'
  n = 5
  allocate ( a(1:n,1:n) )
  call dif1 ( n, n, a )
  call dif1_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DIF1CYCLIC matrix.
!
  title = 'DIF1CYCLIC'
  n = 5
  allocate ( a(1:n,1:n) )
  call dif1cyclic ( n, a )
  call dif1cyclic_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DIF2 matrix.
!
  title = 'DIF2'
  n = 5
  allocate ( a(1:n,1:n) )
  call dif2 ( n, n, a )
  call dif2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DIF2CYCLIC matrix.
!
  title = 'DIF2CYCLIC'
  n = 5
  allocate ( a(1:n,1:n) )
  call dif2cyclic ( n, a )
  call dif2cyclic_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  DORR matrix.
!
  title = 'DORR'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call dorr ( alpha, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  DOWNSHIFT matrix.
!
  title = 'DOWNSHIFT'
  n = 5
  allocate ( a(1:n,1:n) )
  call downshift ( n, a )
  call downshift_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  EBERLEIN matrix.
!
  title = 'EBERLEIN'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call eberlein ( alpha, n, a )
  call eberlein_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  EULERIAN matrix.
!
  title = 'EULERIAN'
  n = 5
  allocate ( a(1:n,1:n) )
  call eulerian ( n, n, a )
  call eulerian_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  EXCHANGE matrix.
!
  title = 'EXCHANGE'
  n = 5
  allocate ( a(1:n,1:n) )
  call exchange ( n, n, a )
  call exchange_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FIBONACCI1 matrix.
!
  title = 'FIBONACCI1'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta ) / 5.0D+00
  call fibonacci1 ( n, alpha, beta, a )
  call fibonacci1_determinant ( n, alpha, beta, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FIBONACCI2 matrix.
!
  title = 'FIBONACCI2'
  n = 5
  allocate ( a(1:n,1:n) )
  call fibonacci2 ( n, a )
  call fibonacci2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FIBONACCI3 matrix.
!
  title = 'FIBONACCI3'
  n = 5
  allocate ( a(1:n,1:n) )
  call fibonacci3 ( n, a )
  call fibonacci3_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FIEDLER.
!
  title = 'FIEDLER'
  n = 7
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call fiedler ( n, n, x, a )
  call fiedler_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  FORSYTHE matrix.
!
  title = 'FORSYTHE'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta ) / 5.0D+00
  call forsythe ( alpha, beta, n, a )
  call forsythe_determinant ( alpha, beta, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FOURIER_COSINE matrix.
!
  title = 'FOURIER_COSINE'
  n = 5
  allocate ( a(1:n,1:n) )
  call fourier_cosine ( n, a )
  call fourier_cosine_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FOURIER_SINE matrix.
!
  title = 'FOURIER_SINE'
  n = 5
  allocate ( a(1:n,1:n) )
  call fourier_sine ( n, a )
  call fourier_sine_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  FRANK matrix.
!
  title = 'FRANK'
  n = 5
  allocate ( a(1:n,1:n) )
  call frank ( n, a )
  call frank_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  GEAR matrix.
!
  do n = 4, 8
    title = 'GEAR'
    allocate ( a(1:n,1:n) )
    seed = 123456789
    ii = i4_uniform ( -n, n, seed )
    jj = i4_uniform ( -n, n, seed )
    call gear ( ii, jj, n, a )
    call gear_determinant ( ii, jj, n, determ1 )
    call r8mat_determinant ( n, a, determ2 )
    norm_frobenius = r8mat_norm_fro ( n, n, a )
    write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      title, n, determ1, determ2, norm_frobenius
    deallocate ( a )
  end do
!
!  GFPP matrix.
!
  title = 'GFPP'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call gfpp ( n, alpha, a )
  call gfpp_determinant ( n, alpha, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  GIVENS matrix.
!
  title = 'GIVENS'
  n = 5
  allocate ( a(1:n,1:n) )
  call givens ( n, n, a )
  call givens_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  GK316 matrix.
!
  title = 'GK316'
  n = 5
  allocate ( a(1:n,1:n) )
  call gk316 ( n, a )
  call gk316_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  GK323 matrix.
!
  title = 'GK323'
  n = 5
  allocate ( a(1:n,1:n) )
  call gk323 ( n, n, a )
  call gk323_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  GK324 matrix.
!
  title = 'GK324'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call gk324 ( n, n, x, a )
  call gk324_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  GRCAR matrix.
!
  title = 'GRCAR'
  n = 5
  seed = 123456789
  allocate ( a(1:n,1:n) )
  k = i4_uniform ( 1, n - 1, seed )
  call grcar ( n, n, k, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  HADAMARD matrix.
!
  title = 'HADAMARD'
  n = 5
  allocate ( a(1:n,1:n) )
  call hadamard ( n, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  HANKEL matrix.
!
  title = 'HANKEL'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:2*n-1) )
  seed = 123456789
  call r8vec_uniform_01 ( 2 * n - 1, seed, x )
  x(1:2*n-1) = anint ( 50.0D+00 * x(1:2*n-1) - 25.0D+00 ) / 5.0D+00
  call hankel ( n, x, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  HANOWA matrix.
!
  title = 'HANOWA'
  n = 6
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call hanowa ( alpha, n, a )
  call hanowa_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HARMAN matrix.
!
  title = 'HARMAN'
  n = 8
  allocate ( a(1:n,1:n) )
  call harman ( a )
  call harman_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HARTLEY matrix.
!
  title = 'HARTLEY'
  do n = 5, 8
    allocate ( a(1:n,1:n) )
    call hartley ( n, a )
    call hartley_determinant ( n, determ1 )
    call r8mat_determinant ( n, a, determ2 )
    norm_frobenius = r8mat_norm_fro ( n, n, a )
    write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      title, n, determ1, determ2, norm_frobenius
    deallocate ( a )
  end do
!
!  HELMERT matrix.
!
  title = 'HELMERT'
  n = 5
  allocate ( a(1:n,1:n) )
  call helmert ( n, a )
  call helmert_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HELMERT2 matrix.
!
  title = 'HELMERT2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call helmert2 ( n, x, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  HERMITE matrix.
!
  title = 'HERMITE'
  n = 5
  allocate ( a(1:n,1:n) )
  call hermite ( n, a )
  call hermite_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HERNDON matrix.
!
  title = 'HERNDON'
  n = 5
  allocate ( a(1:n,1:n) )
  call herndon ( n, a )
  call herndon_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HILBERT matrix.
!
  title = 'HILBERT'
  n = 5
  allocate ( a(1:n,1:n) )
  call hilbert ( n, n, a )
  call hilbert_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  HOUSEHOLDER matrix.
!
  title = 'HOUSEHOLDER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call householder ( n, x, a )
  call householder_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  IDEM_RANDOM matrix.
!
  title = 'IDEM_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  rank = i4_uniform ( 0, n, seed )
  call idem_random ( n, rank, seed, a )
  call idem_random_determinant ( n, rank, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  IDENTITY matrix.
!
  title = 'IDENTITY'
  n = 5
  allocate ( a(1:n,1:n) )
  call identity ( n, n, a )
  call identity_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  IJFACT1 matrix.
!
  title = 'IJFACT1'
  n = 5
  allocate ( a(1:n,1:n) )
  call ijfact1 ( n, a )
  call ijfact1_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  IJFACT2 matrix.
!
  title = 'IJFACT2'
  n = 5
  allocate ( a(1:n,1:n) )
  call ijfact2 ( n, a )
  call ijfact2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ILL3 matrix.
!
  title = 'ILL3'
  n = 3
  allocate ( a(1:n,1:n) )
  call ill3 ( a )
  call ill3_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  INTEGRATION matrix.
!
  title = 'INTEGRATION'
  n = 6
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call integration ( alpha, n, a )
  call integration_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  INVOL matrix.
!
  title = 'INVOL'
  n = 5
  allocate ( a(1:n,1:n) )
  call invol ( n, a )
  call invol_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  INVOL_RANDOM matrix.
!
  title = 'INVOL_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  rank = i4_uniform ( 0, n, seed )
  call invol_random ( n, rank, seed, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  JACOBI matrix.
!
  title = 'JACOBI'
  n = 5
  allocate ( a(1:n,1:n) )
  call jacobi ( n, n, a )
  call jacobi_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  JACOBI matrix.
!
  title = 'JACOBI'
  n = 6
  allocate ( a(1:n,1:n) )
  call jacobi ( n, n, a )
  call jacobi_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  JORDAN matrix.
!
  title = 'JORDAN'
  n = 6
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call jordan ( alpha, n, n, a )
  call jordan_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  KAHAN matrix.
!
  title = 'KAHAN'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call kahan ( alpha, n, n, a )
  call kahan_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  KERSHAW matrix.
!
  title = 'KERSHAW'
  n = 4
  allocate ( a(1:n,1:n) )
  call kershaw ( a )
  call kershaw_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  KERSHAWTRI matrix.
!
  title = 'KERSHAWTRI'
  n = 5
  x_n = ( n + 1 ) / 2
  allocate ( a(1:n,1:n) )
  allocate ( x(1:x_n) )
  seed = 123456789
  call r8vec_uniform_01 ( x_n, seed, x )
  x(1:x_n) = anint ( 50.0D+00 * x(1:x_n) - 25.0D+00 ) / 5.0D+00
  call kershawtri ( n, x, a )
  call kershawtri_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  KMS matrix.
!
  title = 'KMS'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call kms ( alpha, n, n, a )
  call kms_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LAGUERRE matrix.
!
  title = 'LAGUERRE'
  n = 5
  allocate ( a(1:n,1:n) )
  call laguerre ( n, a )
  call laguerre_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LEHMER matrix.
!
  title = 'LEHMER'
  n = 5
  allocate ( a(1:n,1:n) )
  call lehmer ( n, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  LESLIE matrix.
!
  title = 'LESLIE'
  n = 4
  allocate ( a(1:n,1:n) )
  b =  0.025D+00
  di = 0.010D+00
  da = 0.100D+00
  call leslie ( b, di, da, a )
  call leslie_determinant ( b, di, da, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LESP matrix.
!
  title = 'LESP'
  n = 5
  allocate ( a(1:n,1:n) )
  call lesp ( n, n, a )
  call lesp_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LIETZKE matrix.
!
  title = 'LIETZKE'
  n = 5
  allocate ( a(1:n,1:n) )
  call lietzke ( n, a )
  call lietzke_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LIGHTS_OUT matrix.
!
  title = 'LIGHTS_OUT'
  row_num = 5
  col_num = 5
  n = row_num * col_num;
  allocate ( a(1:row_num*col_num,1:row_num*col_num) )
  call lights_out ( row_num, col_num, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  LINE_ADJ matrix.
!
  title = 'LINE_ADJ'
  n = 5
  allocate ( a(1:n,1:n) )
  call line_adj ( n, a )
  call line_adj_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LINE_LOOP_ADJ matrix.
!
  title = 'LINE_LOOP_ADJ'
  n = 5
  allocate ( a(1:n,1:n) )
  call line_loop_adj ( n, a )
  call line_loop_adj_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  LOEWNER matrix.
!
  title = 'LOEWNER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( w(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )
  allocate ( z(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, w )
  call r8vec_uniform_01 ( n, seed, x )
  call r8vec_uniform_01 ( n, seed, y )
  call r8vec_uniform_01 ( n, seed, z )
  call loewner ( w, x, y, z, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( w )
  deallocate ( x )
  deallocate ( y )
  deallocate ( z )
!
!  LOTKIN matrix.
!
  title = 'LOTKIN'
  n = 5
  allocate ( a(1:n,1:n) )
  call lotkin ( n, n, a )
  call lotkin_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  MARKOV_RANDOM matrix.
!
  title = 'MARKOV_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  call markov_random ( n, seed, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  MAXIJ matrix.
!
  title = 'MAXIJ'
  n = 5
  allocate ( a(1:n,1:n) )
  call maxij ( n, n, a )
  call maxij_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  MILNES matrix.
!
  title = 'MILNES'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call milnes ( n, n, x, a )
  call milnes_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  MINIJ matrix.
!
  title = 'MINIJ'
  n = 5
  allocate ( a(1:n,1:n) )
  call minij ( n, n, a )
  call minij_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  MOLER1 matrix.
!
  title = 'MOLER1'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call moler1 ( alpha, n, n, a )
  call moler1_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  MOLER2 matrix.
!
  title = 'MOLER2'
  n = 5
  allocate ( a(1:n,1:n) )
  call moler2 ( a )
  call moler2_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  MOLER3 matrix.
!
  title = 'MOLER3'
  n = 5
  allocate ( a(1:n,1:n) )
  call moler3 ( n, n, a )
  call moler3_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  NEUMANN matrix.
!
  title = 'NEUMANN'
  row_num = 5
  col_num = 5
  allocate ( a(1:row_num*col_num,1:row_num*col_num) )
  call neumann ( row_num, col_num, n, a )
  call neumann_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ONE matrix.
!
  title = 'ONE'
  n = 5
  allocate ( a(1:n,1:n) )
  call one ( n, n, a )
  call one_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ORTEGA matrix.
!
  title = 'ORTEGA'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( v1(1:n) )
  allocate ( v2(1:n) )
  allocate ( v3(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, v1 )
  v1(1:n) = anint ( 50.0D+00 * v1(1:n) - 25.0D+00  ) / 5.0D+00
  call r8vec_uniform_01 ( n, seed, v2 )
  v2(1:n) = anint ( 50.0D+00 * v2(1:n) - 25.0D+00  ) / 5.0D+00
  call r8vec_uniform_01 ( n, seed, v3 )
  v3(1:n) = anint ( 50.0D+00 * v3(1:n) - 25.0D+00  ) / 5.0D+00
  call ortega ( n, v1, v2, v3, a )
  call ortega_determinant ( n, v1, v2, v3, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( v1 )
  deallocate ( v2 )
  deallocate ( v3 )
!
!  ORTH_RANDOM matrix.
!
  title = 'ORTH_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  call orth_random ( n, seed, a )
  call orth_random_determinant ( n, seed, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ORTH_SYMM matrix.
!
  title = 'ORTH_SYMM'
  n = 5
  allocate ( a(1:n,1:n) )
  call orth_symm ( n, a )
  call orth_symm_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  OTO matrix.
!
  title = 'OTO'
  n = 5
  allocate ( a(1:n,1:n) )
  call oto ( n, n, a )
  call oto_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PARTER matrix.
!
  title = 'PARTER'
  n = 5
  allocate ( a(1:n,1:n) )
  call parter ( n, n, a )
  call parter_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PASCAL1 matrix.
!
  title = 'PASCAL1'
  n = 5
  allocate ( a(1:n,1:n) )
  call pascal1 ( n, a )
  call pascal1_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PASCAL2 matrix.
!
  title = 'PASCAL2'
  n = 5
  allocate ( a(1:n,1:n) )
  call pascal2 ( n, a )
  call pascal2_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PASCAL3 matrix.
!
  title = 'PASCAL3'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call pascal3 ( n, alpha, a )
  call pascal3_determinant ( n, alpha, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PDS_RANDOM matrix.
!
  title = 'PDS_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  seed_save = seed
  call pds_random ( n, seed, a )
  seed = seed_save
  call pds_random_determinant ( n, seed, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PEI matrix.
!
  title = 'PEI'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call pei ( alpha, n, a )
  call pei_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PERMUTATION_RANDOM matrix.
!
  title = 'PERMUTATION_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  seed_save = seed
  call permutation_random ( n, seed, a )
  seed = seed_save
  call permutation_random_determinant ( n, seed, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PLU matrix.
!
  title = 'PLU'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( l(1:n,1:n) )
  allocate ( p(1:n,1:n) )
  allocate ( pivot(n) )
  allocate ( u(1:n,1:n) )
  do i = 1, n
    pivot(i) = i
  end do
  call plu ( n, pivot, p, l, u, a )
  call plu_determinant ( n, p, l, u, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( pivot )
  deallocate ( u )
!
!  POISSON matrix.
!
  title = 'POISSON'
  row_num = 5
  col_num = 5
  n = row_num * col_num
  allocate ( a(1:n,1:n) )
  call poisson ( row_num, col_num, n, a )
  call poisson_determinant ( row_num, col_num, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  PROLATE matrix.
!
  title = 'PROLATE'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call prolate ( alpha, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  RECTANGLE_ADJ matrix.
!
  title = 'RECTANGLE_ADJ'
  row_num = 5
  col_num = 5
  n = row_num * col_num
  allocate ( a(1:n,1:n) )
  call rectangle_adj ( row_num, col_num, n, a )
  call rectangle_adj_determinant ( row_num, col_num, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  REDHEFFER matrix.
!
  title = 'REDHEFFER'
  n = 5
  allocate ( a(1:n,1:n) )
  call redheffer ( n, a )
  call redheffer_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  REF_RANDOM matrix.
!
  title = 'REF_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  prob = 0.65D+00
  seed_save = 123456789
  seed = seed_save
  call ref_random ( n, n, prob, seed, a )
  seed = seed_save
  call ref_random_determinant ( n, prob, seed, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  REF_RANDOM matrix.
!
  title = 'REF_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  prob = 0.85D+00
  seed_save = 123456789
  seed = seed_save
  call ref_random ( n, n, prob, seed, a )
  seed = seed_save
  call ref_random_determinant ( n, prob, seed, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  RIEMANN matrix.
!
  title = 'RIEMANN'
  n = 5
  allocate ( a(1:n,1:n) )
  call riemann ( n, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  RING_ADJ matrix.
!
  do n = 1, 8
    title = 'RING_ADJ'
    allocate ( a(1:n,1:n) )
    call ring_adj ( n, a )
    call ring_adj_determinant ( n, determ1 )
    call r8mat_determinant ( n, a, determ2 )
    norm_frobenius = r8mat_norm_fro ( n, n, a )
    write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      title, n, determ1, determ2, norm_frobenius
    deallocate ( a )
  end do
!
!  RIS matrix.
!
  title = 'RIS'
  n = 5
  allocate ( a(1:n,1:n) )
  call ris ( n, a )
  call ris_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  RODMAN matrix.
!
  title = 'RODMAN'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call rodman ( alpha, n, n, a )
  call rodman_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ROSSER1 matrix.
!
!  Note that while the correct determinant of this matrix is 0,
!  that value is very difficult to calculate correctly.  MATLAB
!  returns det ( A ) = -10611, for instance.
!
  title = 'ROSSER1'
  n = 8
  allocate ( a(1:n,1:n) )
  call rosser1 ( a )
  call rosser1_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ROUTH matrix.
!
  title = 'ROUTH'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n ) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) ) / 5.0D+00
  call routh ( n, x, a )
  call routh_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  RUTIS1 matrix.
!
  title = 'RUTIS1'
  n = 4
  allocate ( a(1:n,1:n) )
  call rutis1 ( a )
  call rutis1_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  RUTIS2 matrix.
!
  title = 'RUTIS2'
  n = 4
  allocate ( a(1:n,1:n) )
  call rutis2 ( a )
  call rutis2_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  RUTIS3 matrix.
!
  title = 'RUTIS3'
  n = 4
  allocate ( a(1:n,1:n) )
  call rutis3 ( a )
  call rutis3_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  RUTIS4 matrix.
!
  title = 'RUTIS4'
  n = 5
  allocate ( a(1:n,1:n) )
  call rutis4 ( n, a )
  call rutis4_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  RUTIS5 matrix.
!
  title = 'RUTIS5'
  n = 4
  allocate ( a(1:n,1:n) )
  call rutis5 ( a )
  call rutis5_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  SCHUR_BLOCK matrix.
!
  title = 'SCHUR_BLOCK'
  n = 5
  x_n = ( n + 1 ) / 2
  y_n = n / 2
  allocate ( a(1:n,1:n) )
  allocate ( x(1:x_n) )
  allocate ( y(1:y_n) )
  seed = 123456789
  call r8vec_uniform_01 ( x_n, seed, x )
  x(1:x_n) = anint ( 50.0D+00 * x(1:x_n) - 25.0D+00 ) / 5.0D+00
  call r8vec_uniform_01 ( y_n, seed, y )
  y(1:y_n) = anint ( 50.0D+00 * y(1:y_n) - 25.0D+00 ) / 5.0D+00
  call schur_block ( n, x, y, a )
  call schur_block_determinant ( n, x, y, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  SKEW_CIRCULANT matrix.
!
  title = 'SKEW_CIRCULANT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call skew_circulant ( n, n, x, a )
  call skew_circulant_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  SPLINE matrix.
!
  title = 'SPLINE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call spline ( n, x, a )
  call spline_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  STIRLING matrix.
!
  title = 'STIRLING'
  n = 5
  allocate ( a(1:n,1:n) )
  call stirling ( n, n, a )
  call stirling_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  STRIPE matrix.
!
  title = 'STRIPE'
  n = 5
  allocate ( a(1:n,1:n) )
  call stripe ( n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  SUMMATION matrix.
!
  title = 'SUMMATION'
  n = 5
  allocate ( a(1:n,1:n) )
  call summation ( n, n, a )
  call summation_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  SWEET1 matrix.
!
  title = 'SWEET1'
  n = 6
  allocate ( a(1:n,1:n) )
  perturb = 0.0D+00
  call sweet1 ( perturb, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  SWEET2 matrix.
!
  title = 'SWEET2'
  n = 6
  allocate ( a(1:n,1:n) )
  perturb = 0.0D+00
  call sweet2 ( perturb, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  SWEET3 matrix.
!
  title = 'SWEET3'
  n = 6
  allocate ( a(1:n,1:n) )
  perturb = 0.0D+00
  call sweet3 ( perturb, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  SWEET4 matrix.
!
  title = 'SWEET4'
  n = 13
  allocate ( a(1:n,1:n) )
  perturb = 0.0D+00
  call sweet4 ( perturb, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  SYLVESTER matrix.
!
  title = 'SYLVESTER'
  n = 5
  x_n = 3 + 1
  y_n = 2 + 1
  allocate ( a(1:n,1:n) )
  allocate ( x(1:x_n) )
  allocate ( y(1:y_n) )
  seed = 123456789
  call r8vec_uniform_01 ( x_n, seed, x )
  x(1:x_n) = anint ( 50.0D+00 * x(1:x_n) - 25.0D+00 ) / 5.0D+00
  call r8vec_uniform_01 ( y_n, seed, y )
  y(1:y_n) = anint ( 50.0D+00 * y(1:y_n) - 25.0D+00 ) / 5.0D+00
  call sylvester ( n, x_n - 1, x, y_n - 1, y, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  SYMM_RANDOM matrix.
!
  title = 'SYMM_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call symm_random ( n, x, seed, a )
  call symm_random_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  TOEPLITZ matrix.
!
  title = 'TOEPLITZ'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:2*n-1) )
  seed = 123456789
  call r8vec_uniform_01 ( 2 * n - 1, seed, x )
  x(1:2*n-1) = anint ( 50.0D+00 * x(1:x_n) - 25.0D+00 ) / 5.0D+00
  call toeplitz ( n, x, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  TOEPLITZ_5DIAG matrix.
!
  title = 'TOEPLITZ_5DIAG'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  d1 = r8_uniform_01 ( seed )
  d1 = anint ( 50.0D+00 * d1 - 25.0D+00 ) / 5.0D+00
  d2 = r8_uniform_01 ( seed )
  d2 = anint ( 50.0D+00 * d2 - 25.0D+00 ) / 5.0D+00
  d3 = r8_uniform_01 ( seed )
  d3 = anint ( 50.0D+00 * d3 - 25.0D+00 ) / 5.0D+00
  d4 = r8_uniform_01 ( seed )
  d4 = anint ( 50.0D+00 * d4 - 25.0D+00 ) / 5.0D+00
  d5 = r8_uniform_01 ( seed )
  d5 = anint ( 50.0D+00 * d5 - 25.0D+00 ) / 5.0D+00
  call toeplitz_5diag ( n, d1, d2, d3, d4, d5, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  TOEPLITZ_5S matrix.
!
  title = 'TOEPLITZ_5S'
  row_num = 5
  col_num = 5
  n = row_num * col_num
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta - 25.0D+00 ) / 5.0D+00
  gamma = r8_uniform_01 ( seed )
  gamma = anint ( 50.0D+00 * gamma - 25.0D+00 ) / 5.0D+00
  call toeplitz_5s ( row_num, col_num, alpha, beta, gamma, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  TOEPLITZ_PDS matrix.
!
  title = 'TOEPLITZ_PDS'
  m = 3
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:m) )
  allocate ( y(1:m) )
  seed = 123456789
  call r8vec_uniform_01 ( m, seed, x )
  call r8vec_uniform_01 ( m, seed, y )
  y_sum = sum ( y(1:m) )
  y(1:m) = y(1:m) / y_sum
  call toeplitz_pds ( m, n, x, y, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  TOURNAMENT_RANDOM matrix.
!
  title = 'TOURNAMENT_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  seed_save = 1123456789
  seed = seed_save
  call tournament_random ( n, seed, a )
  seed = seed_save
  call tournament_random_determinant ( n, seed, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  TRANSITION_RANDOM matrix.
!
  title = 'TRANSITION_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  call transition_random ( n, seed, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  TRENCH matrix.
!
  title = 'TRENCH'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call trench ( alpha, n, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  TRI_UPPER matrix.
!
  title = 'TRI_UPPER'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call tri_upper ( alpha, n, a )
  call tri_upper_determinant ( alpha, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  TRIS matrix.
!
  title = 'TRIS'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta - 25.0D+00 ) / 5.0D+00
  gamma = r8_uniform_01 ( seed )
  gamma = anint ( 50.0D+00 * gamma - 25.0D+00 ) / 5.0D+00
  call tris ( n, n, alpha, beta, gamma, a )
  call tris_determinant ( n, alpha, beta, gamma, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  TRIV matrix.
!
  title = 'TRIV'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n-1) )
  allocate ( y(1:n) )
  allocate ( z(1:n-1) )
  seed = 123456789
  call r8vec_uniform_01 ( n - 1, seed, x )
  x(1:n-1) = anint ( 50.0D+00 * x(1:n-1) - 25.0D+00 ) / 5.0D+00
  call r8vec_uniform_01 ( n, seed, y )
  y(1:n) = anint ( 50.0D+00 * y(1:n) - 25.0D+00 ) / 5.0D+00
  call r8vec_uniform_01 ( n - 1, seed, z )
  z(1:n-1) = anint ( 50.0D+00 * z(1:n-1) - 25.0D+00 ) / 5.0D+00
  call triv ( n, x, y, z, a )
  call triv_determinant ( n, x, y, z, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
  deallocate ( z )
!
!  TRIW matrix.
!
  title = 'TRIW'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  k = i4_uniform ( 0, n - 1, seed )
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call triw ( alpha, k, n, a )
  call triw_determinant ( alpha, k, n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  UPSHIFT matrix.
!
  title = 'UPSHIFT'
  n = 5
  allocate ( a(1:n,1:n) )
  call upshift ( n, a )
  call upshift_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  VAND1 matrix.
!
  title = 'VAND1'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call vand1 ( n, x, a )
  call vand1_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  VAND2 matrix.
!
  title = 'VAND2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call vand2 ( n, x, a )
  call vand2_determinant ( n, x, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
  deallocate ( x )
!
!  WATHEN matrix.
!
  title = 'WATHEN'
  row_num = 5
  col_num = 5
  call wathen_order ( row_num, col_num, n )
  allocate ( a(1:n,1:n) )
  call wathen ( row_num, col_num, n, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  WILK03 matrix.
!
  title = 'WILK03'
  n = 3
  allocate ( a(1:n,1:n) )
  call wilk03 ( a )
  call wilk03_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  WILK04 matrix.
!
  title = 'WILK04'
  n = 4
  allocate ( a(1:n,1:n) )
  call wilk04 ( a )
  call wilk04_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  WILK05 matrix.
!
  title = 'WILK05'
  n = 5
  allocate ( a(1:n,1:n) )
  call wilk05 ( a )
  call wilk05_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  WILK12 matrix.
!
  title = 'WILK12'
  n = 12
  allocate ( a(1:n,1:n) )
  call wilk12 ( a )
  call wilk12_determinant ( determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  WILK20 matrix.
!
  title = 'WILK20'
  n = 20
  allocate ( a(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call wilk20 ( alpha, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )
!
!  WILK21 matrix.
!
  title = 'WILK21'
  n = 21
  allocate ( a(1:n,1:n) )
  call wilk21 ( n, a )
  call wilk21_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  WILSON matrix.
!
  title = 'WILSON'
  n = 4
  allocate ( a(1:n,1:n) )
  call wilson ( a )
  call wilson_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ZERO matrix.
!
  title = 'ZERO'
  n = 5
  allocate ( a(1:n,1:n) )
  call zero ( n, n, a )
  call zero_determinant ( n, determ1 )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, determ1, determ2, norm_frobenius
  deallocate ( a )
!
!  ZIELKE matrix.
!
  title = 'ZIELKE'
  n = 5
  allocate ( a(1:n,1:n) )
  seed = 123456789
  d1 = r8_uniform_01 ( seed )
  d1 = anint ( 50.0D+00 * d1 - 25.0D+00 ) / 5.0D+00
  d2 = r8_uniform_01 ( seed )
  d2 = anint ( 50.0D+00 * d2 - 25.0D+00 ) / 5.0D+00
  d3 = r8_uniform_01 ( seed )
  d3 = anint ( 50.0D+00 * d3 - 25.0D+00 ) / 5.0D+00
  call zielke ( n, d1, d2, d3, a )
  call r8mat_determinant ( n, a, determ2 )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) &
    title, n,          determ2, norm_frobenius
  deallocate ( a )

  return
end
subroutine test_eigen ( )

!*****************************************************************************80
!
!! TEST_EIGEN tests the eigenvalue computations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) error_frobenius
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 8 ), allocatable, dimension ( : ) :: lambda
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm_frobenius
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_save
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: v1
  real ( kind = 8 ), allocatable, dimension ( : ) :: v2
  real ( kind = 8 ), allocatable, dimension ( : ) :: v3
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_EIGEN'
  write ( *, '(a)' ) '  Compute the Frobenius norm of the eigenvalue error:'
  write ( *, '(a)' ) '    A * X - X * LAMBDA'
  write ( *, '(a)' ) '  given a set of K eigenvectors X and eigenvalues LAMBDA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix title             N     K      ||A||          ||(A-Lambda*I)*X||'
  write ( *, '(a)' ) ' '
!
!  BODEWIG matrix.
!
  title = 'BODEWIG'
  n = 4
  k = 4
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )

  call bodewig ( a )
  call bodewig_eigenvalues ( lambda )
  call bodewig_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  CARRY matrix.
!
  title = 'CARRY'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  seed = 123456789
  i1 = i4_uniform ( 2, 20, seed )
  call carry ( i1, n, a )
  call carry_eigenvalues ( i1, n, lambda )
  call carry_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  CHOW matrix.
!
  title = 'CHOW'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta ) / 5.0D+00
  call chow ( alpha, beta, n, n, a )
  call chow_eigenvalues ( alpha, beta, n, lambda )
  call chow_right ( alpha, beta, n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  COMBIN matrix.
!
  title = 'COMBIN'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta ) / 5.0D+00
  call combin ( alpha, beta, n, a )
  call combin_eigenvalues ( alpha, beta, n, lambda )
  call combin_right ( alpha, beta, n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  DIF2 matrix.
!
  title = 'DIF2'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call dif2 ( n, n, a )
  call dif2_eigenvalues ( n, lambda )
  call dif2_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  EXCHANGE matrix.
!
  title = 'EXCHANGE'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call exchange ( n, n, a )
  call exchange_eigenvalues ( n, lambda )
  call exchange_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  IDEM_RANDOM matrix.
!
  title = 'IDEM_RANDOM'
  n = 5
  k = 5
  rank = 3
  seed_save = 987654321
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  seed = seed_save
  call idem_random ( n, rank, seed, a )
  call idem_random_eigenvalues ( n, rank, lambda )
  seed = seed_save
  call idem_random_right ( n, rank, seed, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  IDENTITY matrix.
!
  title = 'IDENTITY'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call identity ( n, n, a )
  call identity_eigenvalues ( n, lambda )
  call identity_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  ILL3 matrix.
!
  title = 'ILL3'
  n = 3
  k = 3
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call ill3 ( a )
  call ill3_eigenvalues ( lambda )
  call ill3_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  KMS matrix.
!
  title = 'KMS'
  n = 5
  k = 5
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call kms ( alpha, n, n, a )
  call kms_eigenvalues ( alpha, n, lambda )
  call kms_right ( alpha, n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  ONE matrix.
!
  title = 'ONE'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call one ( n, n, a )
  call one_eigenvalues ( n, lambda )
  call one_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  ORTEGA matrix.
!
  title = 'ORTEGA'
  n = 5
  k = n
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:n) )
  allocate ( v1(1:n) )
  allocate ( v2(1:n) )
  allocate ( v3(1:n) )
  allocate ( x(1:n,1:k) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, v1 )
  v1(1:n) = anint ( 50.0D+00 * v1(1:n) - 25.0D+00  ) / 5.0D+00
  call r8vec_uniform_01 ( n, seed, v2 )
  v2(1:n) = anint ( 50.0D+00 * v2(1:n) - 25.0D+00  ) / 5.0D+00
  call r8vec_uniform_01 ( n, seed, v3 )
  v3(1:n) = anint ( 50.0D+00 * v3(1:n) - 25.0D+00  ) / 5.0D+00
  call ortega ( n, v1, v2, v3, a )
  call ortega_eigenvalues ( n, v1, v2, v3, lambda )
  call ortega_right ( n, v1, v2, v3, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( v1 )
  deallocate ( v2 )
  deallocate ( v3 )
  deallocate ( x )
!
!  OTO matrix.
!
  title = 'OTO'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call oto ( n, n, a )
  call oto_eigenvalues ( n, lambda )
  call oto_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  PDS_RANDOM matrix.
!
  title = 'PDS_RANDOM'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  seed_save = 123456789
  seed = seed_save
  call pds_random ( n, seed, a )
  seed = seed_save
  call pds_random_eigenvalues ( n, seed, lambda )
  seed = seed_save
  call pds_random_right ( n, seed, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  PEI matrix.
!
  title = 'PEI'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call pei ( alpha, n, a )
  call pei_eigenvalues ( alpha, n, lambda )
  call pei_right ( alpha, n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  RODMAN matrix.
!
  title = 'RODMAN'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call rodman ( alpha, n, n, a )
  call rodman_eigenvalues ( alpha, n, lambda )
  call rodman_right ( alpha, n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  ROSSER1 matrix.
!
  title = 'ROSSER1'
  n = 8
  k = 8
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call rosser1 ( a )
  call rosser1_eigenvalues ( lambda )
  call rosser1_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  RUTIS1 matrix.
!
  title = 'RUTIS1'
  n = 4
  k = 4
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call rutis1 ( a )
  call rutis1_eigenvalues ( lambda )
  call rutis1_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  RUTIS2 matrix.
!
  title = 'RUTIS2'
  n = 4
  k = 4
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call rutis2 ( a )
  call rutis2_eigenvalues ( lambda )
  call rutis2_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  RUTIS3 matrix.
!  COMPLEX eigenvalues cannot be handled yet!
!
  if ( .false. ) then
    title = 'RUTIS3'
    n = 4
    k = 4
    allocate ( a(1:n,1:n) )
    allocate ( lambda(1:k) )
    allocate ( x(1:n,1:k) )
    call rutis3 ( a )
    call rutis3_eigenvalues ( lambda )
    call rutis3_right ( x )
!   call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
    norm_frobenius = r8mat_norm_fro ( n, n, a )
    write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
      title, n, k, norm_frobenius, error_frobenius
    deallocate ( a )
    deallocate ( lambda )
    deallocate ( x )

  end if
!
!  RUTIS5 matrix.
!
  title = 'RUTIS5'
  n = 4
  k = 4
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call rutis5 ( a )
  call rutis5_eigenvalues ( lambda )
  call rutis5_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  WILK12 matrix.
!
  title = 'WILK12'
  n = 12
  k = 12
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call wilk12 ( a )
  call wilk12_eigenvalues ( lambda )
  call wilk12_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  WILSON matrix.
!
  title = 'WILSON'
  n = 4
  k = 4
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call wilson ( a )
  call wilson_eigenvalues ( lambda )
  call wilson_right ( x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )
!
!  ZERO matrix.
!
  title = 'ZERO'
  n = 5
  k = 5
  allocate ( a(1:n,1:n) )
  allocate ( lambda(1:k) )
  allocate ( x(1:n,1:k) )
  call zero ( n, n, a )
  call zero_eigenvalues ( n, lambda )
  call zero_right ( n, x )
  call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( lambda )
  deallocate ( x )

  return
end
subroutine test_inverse ( )

!*****************************************************************************80
!
!! TEST_INVERSE tests the inverse computations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 June 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: b
  real ( kind = 8 ) beta
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: c
  real ( kind = 8 ) error_ab
  real ( kind = 8 ) error_ac
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: l
  integer ( kind = 4 ) n
  real ( kind = 8 ) norma_frobenius
  real ( kind = 8 ) normc_frobenius
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p
  integer ( kind = 4 ), allocatable, dimension ( : ) :: pivot
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_save
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: u
  real ( kind = 8 ), allocatable, dimension ( : ) :: v1
  real ( kind = 8 ), allocatable, dimension ( : ) :: v2
  real ( kind = 8 ), allocatable, dimension ( : ) :: v3
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  integer ( kind = 4 ) x_n
  real ( kind = 8 ), allocatable, dimension ( : ) :: y
  integer ( kind = 4 ) y_n
  real ( kind = 8 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INVERSE'
  write ( *, '(a)' ) '  A = a test matrix of order N;'
  write ( *, '(a)' ) '  B = inverse as computed by a routine.'
  write ( *, '(a)' ) '  C = inverse as computed by R8MAT_INVERSE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ||A||    = Frobenius norm of A.'
  write ( *, '(a)' ) '  ||C||    = Frobenius norm of C.'
  write ( *, '(a)' ) '  ||I-AC|| = Frobenius norm of I-A*C.'
  write ( *, '(a)' ) '  ||I-AB|| = Frobenius norm of I-A*B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix title             N        ' // &
                     '||A||          ||C||      ||I-AC||        ||I-AB||'
  write ( *, '(a)' ) ' '
!
!  AEGERTER matrix.
!
  title = 'AEGERTER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call aegerter ( n, a )
  call aegerter_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BAB matrix.
!
  title = 'BAB'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta - 25.0D+00 ) / 5.0D+00
  call bab ( n, alpha, beta, a );
  call bab_inverse ( n, alpha, beta, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BERNSTEIN matrix.
!
  title = 'BERNSTEIN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call bernstein ( n, a )
  call bernstein_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BIS matrix.
!
  title = 'BIS'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta ) / 5.0D+00
  call bis ( alpha, beta, n, n, a );
  call bis_inverse ( alpha, beta, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BODEWIG matrix.
!
  title = 'BODEWIG'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call bodewig ( a )
  call bodewig_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BOOTHROYD matrix.
!
  title = 'BOOTHROYD'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call boothroyd ( n, a )
  call boothroyd_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  BORDERBAND matrix.
!
  title = 'BORDERBAND'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call borderband ( n, a )
  call borderband_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CARRY matrix.
!
  title = 'CARRY'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  k = i4_uniform ( 2, 20, seed )
  call carry ( k, n, a )
  call carry_inverse ( k, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CAUCHY matrix.
!
  title = 'CAUCHY'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  call r8vec_uniform_01 ( n, seed, y )
  call cauchy ( n, x, y, a )
  call cauchy_inverse ( n, x, y, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
  deallocate ( y )
!
!  CHEBY_T matrix.
!
  title = 'CHEBY_T'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call cheby_t ( n, a )
  call cheby_t_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CHEBY_U matrix.
!
  title = 'CHEBY_U'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call cheby_u ( n, a )
  call cheby_u_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CHEBY_VAN2 matrix.
!
  title = 'CHEBY_VAN2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call cheby_van2 ( n, a )
  call cheby_van2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CHEBY_VAN3 matrix.
!
  title = 'CHEBY_VAN3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call cheby_van3 ( n, a )
  call cheby_van3_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CHOW matrix.
!
  title = 'CHOW'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta ) / 5.0D+00
  call chow ( alpha, beta, n, n, a )
  call chow_inverse ( alpha, beta, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CIRCULANT matrix.
!
  title = 'CIRCULANT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call circulant ( n, n, x, a )
  call circulant_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  CIRCULANT2 matrix.
!
  title = 'CIRCULANT2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call circulant2 ( n, a )
  call circulant2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CLEMENT1 matrix.
!
  title = 'CLEMENT1'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call clement1 ( n, a )
  call clement1_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CLEMENT2 matrix.
!
  title = 'CLEMENT2'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call clement2 ( n, a )
  call clement2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CLEMENT3.
!
  title = 'CLEMENT3'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n-1) )
  allocate ( y(1:n-1) )
  seed = 123456789
  call r8vec_uniform_01 ( n - 1, seed, x )
  x(1:n-1) = anint ( 50.0D+00 * x(1:n-1) - 25.0D+00 ) / 5.0D+00
  call r8vec_uniform_01 ( n - 1, seed, y )
  y(1:n-1) = anint ( 50.0D+00 * y(1:n-1) - 25.0D+00 ) / 5.0D+00
  call clement3 ( n, x, y, a )
  call clement3_inverse ( n, x, y, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
  deallocate ( y )
!
!  COMBIN matrix.
!
  title = 'COMBIN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta ) / 5.0D+00
  call combin ( alpha, beta, n, a )
  call combin_inverse ( alpha, beta, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  COMPANION.
!
  title = 'COMPANION'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 10.0D+00 * x(1:n) - 5.0D+00 )
  call companion ( n, x, a )
  call companion_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  COMPLEX_I
!
  title = 'COMPLEX_I'
  n = 2
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call complex_i ( a )
  call complex_i_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CONEX1 matrix.
!
  title = 'CONEX1'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call conex1 ( alpha, a )
  call conex1_inverse ( alpha, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CONEX2 matrix.
!
  title = 'CONEX2'
  n = 3
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call conex2 ( alpha, a )
  call conex2_inverse ( alpha, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CONEX3 matrix.
!
  title = 'CONEX3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call conex3 ( n, a )
  call conex3_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  CONFERENCE matrix.
!
  title = 'CONFERENCE'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call conference ( n, a )
  call conference_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DAUB2.
!
  title = 'DAUB2'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call daub2 ( n, a )
  call daub2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DAUB4.
!
  title = 'DAUB4'
  n = 8
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call daub4 ( n, a )
  call daub4_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DAUB6.
!
  title = 'DAUB6'
  n = 12
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call daub6 ( n, a )
  call daub6_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DAUB8.
!
  title = 'DAUB8'
  n = 16
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call daub8 ( n, a )
  call daub8_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DAUB10.
!
  title = 'DAUB10'
  n = 20
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call daub10 ( n, a )
  call daub10_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DAUB12.
!
  title = 'DAUB12'
  n = 24
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call daub12 ( n, a )
  call daub12_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DIAGONAL.
!
  title = 'DIAGONAL'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call diagonal ( n, n, x, a )
  call diagonal_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  DIF2 matrix.
!
  title = 'DIF2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call dif2 ( n, n, a )
  call dif2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DOWNSHIFT matrix.
!
  title = 'DOWNSHIFT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call downshift ( n, a )
  call downshift_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  DRMAC
!

!
!  EULERIAN matrix.
!
  title = 'EULERIAN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call eulerian ( n, n, a )
  call eulerian_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  EXCHANGE matrix.
!
  title = 'EXCHANGE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call exchange ( n, n, a )
  call exchange_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  FIBONACCI2 matrix.
!
  title = 'FIBONACCI2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call fibonacci2 ( n, a )
  call fibonacci2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  FIBONACCI3 matrix.
!
  title = 'FIBONACCI3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call fibonacci3 ( n, a )
  call fibonacci3_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  FIEDLER.
!  The FIEDLER_INVERSE routine assumes the X vector is sorted.
!
  title = 'FIEDLER'
  n = 7
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call r8vec_sort_bubble_a ( n, x )
  call fiedler ( n, n, x, a )
  call fiedler_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  FORSYTHE matrix.
!
  title = 'FORSYTHE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta ) / 5.0D+00
  call forsythe ( alpha, beta, n, a )
  call forsythe_inverse ( alpha, beta, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  FOURIER_COSINE matrix.
!
  title = 'FOURIER_COSINE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call fourier_cosine ( n, a )
  call fourier_cosine_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  FOURIER_SINE matrix.
!
  title = 'FOURIER_SINE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call fourier_sine ( n, a )
  call fourier_sine_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  FRANK matrix.
!
  title = 'FRANK'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call frank ( n, a )
  call frank_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  GFPP matrix.
!
  title = 'GFPP'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  call gfpp ( n, alpha, a )
  call gfpp_inverse ( n, alpha, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  GIVENS matrix.
!
  title = 'GIVENS'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call givens ( n, n, a )
  call givens_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  GK316 matrix.
!
  title = 'GK316'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call gk316 ( n, a )
  call gk316_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  GK323 matrix.
!
  title = 'GK323'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call gk323 ( n, n, a )
  call gk323_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  GK324 matrix.
!
  title = 'GK324'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call gk324 ( n, n, x, a )
  call gk324_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  HANOWA matrix.
!
  title = 'HANOWA'
  n = 8
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call hanowa ( alpha, n, a )
  call hanowa_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HARMAN matrix.
!
  title = 'HARMAN'
  n = 8
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call harman ( a )
  call harman_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HARTLEY matrix.
!
  title = 'HARTLEY'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call hartley ( n, a )
  call hartley_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HELMERT matrix.
!
  title = 'HELMERT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call helmert ( n, a )
  call helmert_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HELMERT2 matrix.
!
  title = 'HELMERT2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call helmert2 ( n, x, a )
  call helmert2_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  HERMITE matrix.
!
  title = 'HERMITE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call hermite ( n, a )
  call hermite_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HERNDON matrix.
!
  title = 'HERNDON'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call herndon ( n, a )
  call herndon_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HILBERT matrix.
!
  title = 'HILBERT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call hilbert ( n, n, a )
  call hilbert_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  HOUSEHOLDER matrix.
!
  title = 'HOUSEHOLDER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call householder ( n, x, a )
  call householder_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  IDENTITY matrix.
!
  title = 'IDENTITY'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call identity ( n, n, a )
  call identity_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  ILL3 matrix.
!
  title = 'ILL3'
  n = 3
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call ill3 ( a )
  call ill3_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  INTEGRATION matrix.
!
  title = 'INTEGRATION'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call integration ( alpha, n, a )
  call integration_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  INVOL matrix.
!
  title = 'INVOL'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call invol ( n, a )
  call invol_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  JORDAN matrix.
!
  title = 'JORDAN'
  n = 6
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call jordan ( alpha, n, n, a )
  call jordan_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  KAHAN matrix.
!
  title = 'KAHAN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call kahan ( alpha, n, n, a )
  call kahan_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  KERSHAW matrix.
!
  title = 'KERSHAW'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call kershaw ( a )
  call kershaw_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  KERSHAWTRI matrix.
!
  title = 'KERSHAWTRI'
  n = 5
  x_n = ( n + 1 ) / 2
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:x_n) )
  seed = 123456789
  call r8vec_uniform_01 ( x_n, seed, x )
  x(1:x_n) = anint ( 50.0D+00 * x(1:x_n) - 25.0D+00 ) / 5.0D+00
  call kershawtri ( n, x, a )
  call kershawtri_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  KMS matrix.
!
  title = 'KMS'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call kms ( alpha, n, n, a )
  call kms_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  LAGUERRE matrix.
!
  title = 'LAGUERRE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call laguerre ( n, a )
  call laguerre_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  LEGENDRE matrix.
!
  title = 'LEGENDRE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call legendre ( n, a )
  call legendre_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  LEHMER matrix.
!
  title = 'LEHMER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call lehmer ( n, n, a )
  call lehmer_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  LIETZKE matrix.
!
  title = 'LIETZKE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call lietzke ( n, a )
  call lietzke_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  LOTKIN matrix.
!
  title = 'LOTKIN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call lotkin ( n, n, a )
  call lotkin_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  MAXIJ matrix.
!
  title = 'MAXIJ'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call maxij ( n, n, a )
  call maxij_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  MILNES matrix.
!
  title = 'MILNES'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call milnes ( n, n, x, a )
  call milnes_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  MINIJ matrix.
!
  title = 'MINIJ'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call minij ( n, n, a )
  call minij_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  MOLER1 matrix.
!
  title = 'MOLER1'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call moler1 ( alpha, n, n, a )
  call moler1_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  MOLER3 matrix.
!
  title = 'MOLER3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call moler3 ( n, n, a )
  call moler3_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  ORTEGA matrix.
!
  title = 'ORTEGA'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( v1(1:n) )
  allocate ( v2(1:n) )
  allocate ( v3(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, v1 )
  v1(1:n) = anint ( 50.0D+00 * v1(1:n) - 25.0D+00  ) / 5.0D+00
  call r8vec_uniform_01 ( n, seed, v2 )
  v2(1:n) = anint ( 50.0D+00 * v2(1:n) - 25.0D+00  ) / 5.0D+00
  call r8vec_uniform_01 ( n, seed, v3 )
  v3(1:n) = anint ( 50.0D+00 * v3(1:n) - 25.0D+00  ) / 5.0D+00
  call ortega ( n, v1, v2, v3, a )
  call ortega_inverse ( n, v1, v2, v3, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( v1 )
  deallocate ( v2 )
  deallocate ( v3 )
!
!  ORTH_SYMM matrix.
!
  title = 'ORTH_SYMM'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call orth_symm ( n, a )
  call orth_symm_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  OTO matrix.
!
  title = 'OTO'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call oto ( n, n, a )
  call oto_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PARTER matrix.
!
  title = 'PARTER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call parter ( n, n, a )
  call parter_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PASCAL1 matrix.
!
  title = 'PASCAL1'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call pascal1 ( n, a )
  call pascal1_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PASCAL2 matrix.
!
  title = 'PASCAL2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call pascal2 ( n, a )
  call pascal2_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PASCAL3 matrix.
!
  title = 'PASCAL3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call pascal3 ( n, alpha, a )
  call pascal3_inverse ( n, alpha, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PDS_RANDOM matrix.
!
  title = 'PDS_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed_save = 123456789
  seed = seed_save
  call pds_random ( n, seed, a )
  seed = seed_save
  call pds_random_inverse ( n, seed, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PEI matrix.
!
  title = 'PEI'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call pei ( alpha, n, a )
  call pei_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PERMUTATION_RANDOM matrix.
!
  title = 'PERMUTATION_RANDOM'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  seed_save = seed
  call permutation_random ( n, seed, a )
  seed = seed_save
  call permutation_random_inverse ( n, seed, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  PLU matrix.
!
  title = 'PLU'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( l(1:n,1:n) )
  allocate ( p(1:n,1:n) )
  allocate ( pivot(1:n) )
  allocate ( u(1:n,1:n) )
  do i = 1, n
    pivot(i) = i
  end do
  call plu ( n, pivot, p, l, u, a )
  call plu_inverse ( n, p, l, u, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( l )
  deallocate ( p )
  deallocate ( pivot )
  deallocate ( u )
!
!  RIS matrix.
!
  title = 'RIS'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call ris ( n, a )
  call ris_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  RODMAN matrix.
!
  title = 'RODMAN'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call rodman ( alpha, n, n, a )
  call rodman_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  RUTIS1 matrix.
!
  title = 'RUTIS1'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call rutis1 ( a )
  call rutis1_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  RUTIS2 matrix.
!
  title = 'RUTIS2'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call rutis2 ( a )
  call rutis2_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  RUTIS3 matrix.
!
  title = 'RUTIS3'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call rutis3 ( a )
  call rutis3_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  RUTIS4 matrix.
!
  title = 'RUTIS4'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call rutis4 ( n, a )
  call rutis4_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  RUTIS5 matrix.
!
  title = 'RUTIS5'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call rutis5 ( a )
  call rutis5_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  SCHUR_BLOCK matrix.
!
  title = 'SCHUR_BLOCK'
  n = 5
  x_n = ( n + 1 ) / 2
  y_n = n / 2
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:x_n) )
  allocate ( y(1:y_n) )
  seed = 123456789
  call r8vec_uniform_01 ( x_n, seed, x )
  x(1:x_n) = anint ( 50.0D+00 * x(1:x_n) - 25.0D+00 ) / 5.0D+00
  call r8vec_uniform_01 ( y_n, seed, y )
  y(1:y_n) = anint ( 50.0D+00 * y(1:y_n) - 25.0D+00 ) / 5.0D+00
  call schur_block ( n, x, y, a )
  call schur_block_inverse ( n, x, y, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
  deallocate ( y )
!
!  SPLINE matrix.
!
  title = 'SPLINE'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n-1) )
  seed = 123456789
  call r8vec_uniform_01 ( n - 1, seed, x )
  x(1:n - 1) = anint ( 50.0D+00 * x(1:n - 1) - 25.0D+00 ) / 5.0D+00
  call spline ( n, x, a )
  call spline_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  STIRLING matrix.
!
  title = 'STIRLING'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call stirling ( n, n, a )
  call stirling_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  SUMMATION matrix.
!
  title = 'SUMMATION'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call summation ( n, n, a )
  call summation_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  TRI_UPPER matrix.
!
  title = 'TRI_UPPER'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call tri_upper ( alpha, n, a )
  call tri_upper_inverse ( alpha, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  TRIS matrix.
!
  title = 'TRIS'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  beta = r8_uniform_01 ( seed )
  beta = anint ( 50.0D+00 * beta - 25.0D+00 ) / 5.0D+00
  gamma = r8_uniform_01 ( seed )
  gamma = anint ( 50.0D+00 * gamma - 25.0D+00 ) / 5.0D+00
  call tris ( n, n, alpha, beta, gamma, a )
  call tris_inverse ( n, alpha, beta, gamma, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  TRIV matrix.
!
  title = 'TRIV'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n-1) )
  allocate ( y(1:n) )
  allocate ( z(1:n-1) )
  seed = 123456789
  call r8vec_uniform_01 ( n - 1, seed, x )
  x(1:n-1) = anint ( 50.0D+00 * x(1:n-1) - 25.0D+00 ) / 5.0D+00
  call r8vec_uniform_01 ( n, seed, y )
  y(1:n) = anint ( 50.0D+00 * y(1:n) - 25.0D+00 ) / 5.0D+00
  call r8vec_uniform_01 ( n - 1, seed, z )
  z(1:n-1) = anint ( 50.0D+00 * z(1:n-1) - 25.0D+00 ) / 5.0D+00
  call triv ( n, x, y, z, a )
  call triv_inverse ( n, x, y, z, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
  deallocate ( y )
  deallocate ( z )
!
!  TRIW matrix.
!
  title = 'TRIW'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  seed = 123456789
  k = i4_uniform ( 0, n - 1, seed )
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call triw ( alpha, k, n, a )
  call triw_inverse ( alpha, k, n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  UPSHIFT matrix.
!
  title = 'UPSHIFT'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call upshift ( n, a )
  call upshift_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  VAND1 matrix.
!
  title = 'VAND1'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call vand1 ( n, x, a )
  call vand1_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  VAND2 matrix.
!
  title = 'VAND2'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )
  x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
  call vand2 ( n, x, a )
  call vand2_inverse ( n, x, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )
!
!  WILK03 matrix.
!
  title = 'WILK03'
  n = 3
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call wilk03 ( a )
  call wilk03_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  WILK04 matrix.
!
  title = 'WILK04'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call wilk04 ( a )
  call wilk04_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  WILK05 matrix.
!
  title = 'WILK05'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call wilk05 ( a )
  call wilk05_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  WILK21 matrix.
!
  title = 'WILK21'
  n = 21
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call wilk21 ( n, a )
  call wilk21_inverse ( n, b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  WILSON matrix.
!
  title = 'WILSON'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
  call wilson ( a )
  call wilson_inverse ( b )
  call r8mat_inverse ( n, a, c )
  call r8mat_is_inverse ( n, a, b, error_ab )
  call r8mat_is_inverse ( n, a, c, error_ac )
  norma_frobenius = r8mat_norm_fro ( n, n, a )
  normc_frobenius = r8mat_norm_fro ( n, n, c )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )

  return
end
subroutine test_null ( )

!*****************************************************************************80
!
!! TEST_NULL tests the null vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: at
  real ( kind = 8 ) alpha
  integer ( kind = 4 ) col_num
  real ( kind = 8 ) error_l2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mt
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nt
  real ( kind = 8 ) norm_a_frobenius
  real ( kind = 8 ) norm_x_l2
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) r8mat_norm_fro
  real ( kind = 8 ) r8vec_norm_l2
  integer ( kind = 4 ) row_num
  integer ( kind = 4 ) seed
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_NULL'
  write ( *, '(a)' ) '  A = a test matrix of order M by N'
  write ( *, '(a)' ) '  x = an N vector, candidate for a null vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ||A|| = Frobenius norm of A.'
  write ( *, '(a)' ) '  ||x|| = L2 norm of x.'
  write ( *, '(a)' ) '  ||A*x||/||x|| = L2 norm of A*x over L2 norm of x.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix title	           M     N      ' // &
               '||A||            ||x||        ||A*x||/||x||'
  write ( *, '(a)' ) ' '
!
!  ARCHIMEDES matrix.
!
  title = 'ARCHIMEDES'
  m = 7
  n = 8
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call archimedes ( a )
  call archimedes_null ( x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  CHEBY_DIFF1 matrix.
!
  title = 'CHEBY_DIFF1'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call cheby_diff1 ( n, a )
  call cheby_diff1_null ( n, x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  CREATION matrix.
!
  title = 'CREATION'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call creation ( m, n, a )
  call creation_null ( n, x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  DIF1 matrix.
!  Only has null vectors for N odd.
!
  title = 'DIF1'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call dif1 ( m, n, a )
  call dif1_null ( n, x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  DIF1CYCLIC matrix.
!
  title = 'DIF1CYCLIC'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call dif1cyclic ( n, a )
  call dif1cyclic_null ( n, x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  DIF2CYCLIC matrix.
!
  title = 'DIF2CYCLIC'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call dif2cyclic ( n, a )
  call dif2cyclic_null ( n, x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  EBERLEIN matrix.
!  We have a LEFT null vector.
!
  title = 'EBERLEIN (left)'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( at(1:n,1:m) )
  allocate ( x(1:m) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call eberlein ( alpha, n, a )
  mt = n
  nt = m
  at(1:mt,1:nt) = transpose ( a(1:m,1:n) )
  call eberlein_null_left ( n, x )
  call r8mat_is_null_vector ( mt, nt, at, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( mt, nt, at )
  norm_x_l2 = r8vec_norm_l2 ( nt, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( at )
  deallocate ( x )
!
!  FIBONACCI1 matrix.
!
  title = 'FIBONACCI1'
  m = 5
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  seed = 123456789
  f1 = r8_uniform_01 ( seed )
  f1 = anint ( 50.0D+00 * f1 - 25.0D+00 ) / 5.0D+00
  f2 = r8_uniform_01 ( seed )
  f2 = anint ( 50.0D+00 * f2 - 25.0D+00 ) / 5.0D+00
  call fibonacci1 ( n, f1, f2, a )
  call fibonacci1_null ( n, f1, f2, x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  LAUCHLI matrix.
!  We have a LEFT null vector of a RECTANGULAR matrix.
!
  title = 'LAUCHLI (left)'
  m = 6
  n = m - 1
  allocate ( a(1:m,1:n) )
  allocate ( at(1:n,1:m) )
  allocate ( x(1:m) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
  call lauchli ( alpha, m, n, a )
  mt = n
  nt = m
  at(1:mt,1:nt) = transpose ( a(1:m,1:n) )
  call lauchli_null_left ( alpha, m, n, x )
  call r8mat_is_null_vector ( mt, nt, at, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( mt, nt, at )
  norm_x_l2 = r8vec_norm_l2 ( nt, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( at )
  deallocate ( x )
!
!  LINE_ADJ matrix.
!
  title = 'LINE_ADJ'
  m = 7
  n = m
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call line_adj ( n, a )
  call line_adj_null ( n, x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  MOLER2 matrix.
!
  title = 'MOLER2'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call moler2 ( a )
  call moler2_null ( x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  NEUMANN matrix.
!
  title = 'NEUMANN'
  row_num = 5
  col_num = 5
  m = row_num * col_num
  n = row_num * col_num
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call neumann ( row_num, col_num, n, a )
  call neumann_null ( row_num, col_num, x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  ONE matrix.
!
  title = 'ONE'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call one ( m, n, a )
  call one_null ( n, x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  RING_ADJ matrix.
!  N must be a multiple of 4 for there to be a null vector.
!
  title = 'RING_ADJ'
  m = 12
  n = 12
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call ring_adj ( n, a )
  call ring_adj_null ( n, x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  ROSSER1 matrix.
!
  title = 'ROSSER1'
  m = 8
  n = 8
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call rosser1 ( a )
  call rosser1_null ( x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )
!
!  ZERO matrix.
!
  title = 'ZERO'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  call zero ( m, n, a )
  call zero_null ( n, x )
  call r8mat_is_null_vector ( m, n, a, x, error_l2 )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  norm_x_l2 = r8vec_norm_l2 ( n, x )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
  deallocate ( a )
  deallocate ( x )

  return
end
subroutine test_plu ( )

!*****************************************************************************80
!
!! TEST_PLU tests the PLU factors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) error_frobenius
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm_a_frobenius
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) seed
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: u

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_PLU'
  write ( *, '(a)' ) '  A = a test matrix of order M by N'
  write ( *, '(a)' ) '  P, L, U are the PLU factors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ||A|| = Frobenius norm of A.'
  write ( *, '(a)' ) '  ||A-PLU|| = Frobenius norm of A-P*L*U.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix title	           M     N      ' // &
               '||A||            ||A-PLU||'
  write ( *, '(a)' ) ' '
!
!  BODEWIG matrix.
!
  title = 'BODEWIG'
  m = 4
  n = 4
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call bodewig ( a )
  call bodewig_plu ( p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  BORDERBAND matrix.
!
  title = 'BORDERBAND'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call borderband ( n, a )
  call borderband_plu ( n, p, l, u )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  DIF2 matrix.
!
  title = 'DIF2'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call dif2 ( m, n, a )
  call dif2_plu ( n, p, l, u )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  GFPP matrix.
!
  title = 'GFPP'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call gfpp ( n, alpha, a )
  call gfpp_plu ( n, alpha, p, l, u )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  GIVENS matrix.
!
  title = 'GIVENS'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call givens ( n, n, a )
  call givens_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  KMS matrix.
!
  title = 'KMS'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call kms ( alpha, m, n, a )
  call kms_plu ( alpha, n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  MAXIJ matrix.
!
  title = 'MAXIJ'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call maxij ( n, n, a )
  call maxij_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  MINIJ matrix.
!
  title = 'MINIJ'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call minij ( n, n, a )
  call minij_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  MOLER1 matrix.
!
  title = 'MOLER1'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  seed = 123456789
  alpha = r8_uniform_01 ( seed )
  alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
  call moler1 ( alpha, n, n, a )
  call moler1_plu ( alpha, n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  MOLER3 matrix.
!
  title = 'MOLER3'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call moler3 ( m, n, a )
  call moler3_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  OTO matrix.
!
  title = 'OTO'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call oto ( m, n, a )
  call oto_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  PASCAL2 matrix.
!
  title = 'PASCAL2'
  m = 5
  n = 5
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call pascal2 ( n, a )
  call pascal2_plu ( n, p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )
!
!  WILSON matrix.
!
  title = 'WILSON'
  m = 4
  n = 4
  allocate ( a(1:m,1:n) )
  allocate ( p(1:m,1:m) )
  allocate ( l(1:m,1:m) )
  allocate ( u(1:m,1:n) )
  call wilson ( a )
  call wilson_plu ( p, l, u  )
  call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
  norm_a_frobenius = r8mat_norm_fro ( m, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    title, m, n, norm_a_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( l )
  deallocate ( p )
  deallocate ( u )

  return
end
subroutine test_solution ( )

!*****************************************************************************80
!
!! TEST_SOLUTION tests the linear solution computations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: b
  real ( kind = 8 ) beta
  real ( kind = 8 ) error_frobenius
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncol
  real ( kind = 8 ) norm_frobenius
  integer ( kind = 4 ) nrow
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) r8mat_norm_fro
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_save
  character ( len = 20 ) title
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SOLUTION'
  write ( *, '(a)' ) '  Compute the Frobenius norm of the solution error:'
  write ( *, '(a)' ) '    A * X - B'
  write ( *, '(a)' ) '  given MxN matrix A, NxK solution X, MxK right hand side B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix title             M     N     K      ||A||         ||A*X-B||'
  write ( *, '(a)' ) ' '
!
!  BODEWIG matrix.
!
  title = 'BODEWIG'
  m = 4
  n = 4
  k = 1
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call bodewig ( a )
  call bodewig_rhs ( b )
  call bodewig_solution ( x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  DIF2 matrix.
!
  title = 'DIF2'
  m = 10
  n = 10
  k = 2
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call dif2 ( m, n, a )
  call dif2_rhs ( m, k, b )
  call dif2_solution ( n, k, x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  FRANK matrix.
!
  title = 'FRANK'
  m = 10
  n = 10
  k = 2
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call frank ( n, a )
  call frank_rhs ( m, k, b )
  call frank_solution ( n, k, x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  POISSON matrix.
!
  title = 'POISSON'
  nrow = 4
  ncol = 5
  m = nrow * ncol
  n = nrow * ncol
  k = 1
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call poisson ( nrow, ncol, n, a )
  call poisson_rhs ( nrow, ncol, n, b )
  call poisson_solution ( nrow, ncol, n, x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  WILK03 matrix.
!
  title = 'WILK03'
  m = 3
  n = 3
  k = 1
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call wilk03 ( a )
  call wilk03_rhs ( b )
  call wilk03_solution ( x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  WILK04 matrix.
!
  title = 'WILK04'
  m = 4
  n = 4
  k = 1
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call wilk04 ( a )
  call wilk04_rhs ( b )
  call wilk04_solution ( x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )
!
!  WILSON matrix.
!
  title = 'WILSON'
  m = 4
  n = 4
  k = 1
  allocate ( a(1:m,1:n) )
  allocate ( b(1:m,1:k) )
  allocate ( x(1:n,1:k) )
  call wilson ( a )
  call wilson_rhs ( b )
  call wilson_solution ( x )
  call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
  norm_frobenius = r8mat_norm_fro ( n, n, a )
  write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) &
    title, m, n, k, norm_frobenius, error_frobenius
  deallocate ( a )
  deallocate ( b )
  deallocate ( x )

  return
end
