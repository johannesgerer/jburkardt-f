program main

!*****************************************************************************80
!
!! CONDITION_PRB tests the CONDITION library.
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
  write ( *, '(a)' ) 'CONDITION_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CONDITION library.'
  write ( *, '(a)' ) '  The R8LIB library must also be available.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CONDITION_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CONDITION_LINPACK.
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

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: a_inverse(:,:)
  real ( kind = 8 ) a_inverse_norm_l1
  real ( kind = 8 ) a_norm_l1
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) cond
  real ( kind = 8 ) cond_l1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) n
  character ( len = 80 ) name
  integer ( kind = 4 ), allocatable :: pivot(:)
  real ( kind = 8 ) r8mat_norm_l1
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: z(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  CONDITION_LINPACK estimates the L1 condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix               Order   Condition         Linpack'
  write ( *, '(a)' ) ' '
!
!  Combinatorial matrix.
!
  name = 'Combinatorial'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  allocate ( pivot(1:n) )
  allocate ( z(1:n) )
  alpha = 2.0D+00
  beta = 3.0D+00
  call combin ( alpha, beta, n, a )
  call combin_inverse ( alpha, beta, n, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  call condition_linpack ( n, a, pivot, cond, z )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
  deallocate ( a )
  deallocate ( a_inverse )
  deallocate ( pivot )
  deallocate ( z )
!
!  CONEX1
!
  name = 'CONEX1'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  allocate ( pivot(1:n) )
  allocate ( z(1:n) )
  alpha = 100.0D+00
  call conex1 ( alpha, a )
  call conex1_inverse ( alpha, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  call condition_linpack ( n, a, pivot, cond, z )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
  deallocate ( a )
  deallocate ( a_inverse )
  deallocate ( pivot )
  deallocate ( z )
!
!  CONEX2
!
  name = 'CONEX2'
  n = 3
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  allocate ( pivot(1:n) )
  allocate ( z(1:n) )
  alpha = 100.0D+00
  call conex2 ( alpha, a )
  call conex2_inverse ( alpha, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  call condition_linpack ( n, a, pivot, cond, z )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
  deallocate ( a )
  deallocate ( a_inverse )
  deallocate ( pivot )
  deallocate ( z )
!
!  CONEX3
!
  name = 'CONEX3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  allocate ( pivot(1:n) )
  allocate ( z(1:n) )
  call conex3 ( n, a )
  call conex3_inverse ( n, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  call condition_linpack ( n, a, pivot, cond, z )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
  deallocate ( a )
  deallocate ( a_inverse )
  deallocate ( pivot )
  deallocate ( z )
!
!  CONEX4
!
  name = 'CONEX4'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  allocate ( pivot(1:n) )
  allocate ( z(1:n) )
  call conex4 ( a )
  call conex4_inverse ( a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  call condition_linpack ( n, a, pivot, cond, z )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
  deallocate ( a )
  deallocate ( a_inverse )
  deallocate ( pivot )
  deallocate ( z )
!
!  KAHAN
!
  name = 'KAHAN'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  allocate ( pivot(1:n) )
  allocate ( z(1:n) )
  alpha = 0.25D+00
  call kahan ( alpha, n, n, a )
  call kahan_inverse ( alpha, n, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  call condition_linpack ( n, a, pivot, cond, z )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
  deallocate ( a )
  deallocate ( a_inverse )
  deallocate ( pivot )
  deallocate ( z )
!
!  Random
!
  seed = 123456789

  do i = 1, 5
    name = 'RANDOM'
    n = 4
    allocate ( a(1:n,1:n) )
    allocate ( a_inverse(1:n,1:n) )
    allocate ( pivot(1:n) )
    allocate ( z(1:n) )
    call r8mat_uniform_01 ( n, n, seed, a )
    a_inverse(1:n,1:n) = a(1:n,1:n)
    call r8ge_fa ( n, a_inverse, pivot, info )
    call r8ge_inverse ( n, a_inverse, pivot )
    a_norm_l1         = r8mat_norm_l1 ( n, n, a )
    a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
    cond_l1           = a_norm_l1 * a_inverse_norm_l1
    call condition_linpack ( n, a, pivot, cond, z )
    write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
    deallocate ( a )
    deallocate ( a_inverse )
    deallocate ( pivot )
    deallocate ( z )
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests CONDITION_SAMPLE1.
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

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: a_inverse(:,:)
  real ( kind = 8 ) a_inverse_norm_l1
  real ( kind = 8 ) a_norm_l1
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) cond
  real ( kind = 8 ) cond_l1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) :: m_test(3) = (/ 10, 1000, 100000 /)
  integer ( kind = 4 ) n
  character ( len = 80 ) name
  integer ( kind = 4 ), allocatable :: pivot(:)
  real ( kind = 8 ) r8mat_norm_l1
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  CONDITION_SAMPLE1 estimates the L1 condition number using sampling.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix                 Samples Order   Condition        Estimate'
!
!  Combinatorial matrix.
!
  name = 'Combinatorial'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  alpha = 2.0D+00
  beta = 3.0D+00
  call combin ( alpha, beta, n, a )
  call combin_inverse ( alpha, beta, n, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  write ( *, '(a)' ) ' '
  do i = 1, 3
    m = m_test(i)
    call condition_sample1 ( n, a, m, cond )
    write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) name, m, n, cond_l1, cond
  end do
  deallocate ( a )
  deallocate ( a_inverse )
!
!  CONEX1
!
  name = 'CONEX1'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  alpha = 100.0D+00
  call conex1 ( alpha, a )
  call conex1_inverse ( alpha, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  write ( *, '(a)' ) ' '
  do i = 1, 3
    m = m_test(i)
    call condition_sample1 ( n, a, m, cond )
    write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) name, m, n, cond_l1, cond
  end do
  deallocate ( a )
  deallocate ( a_inverse )
!
!  CONEX2
!
  name = 'CONEX2'
  n = 3
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  alpha = 100.0D+00
  call conex2 ( alpha, a )
  call conex2_inverse ( alpha, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  write ( *, '(a)' ) ' '
  do i = 1, 3
    m = m_test(i)
    call condition_sample1 ( n, a, m, cond )
    write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) name, m, n, cond_l1, cond
  end do
  deallocate ( a )
  deallocate ( a_inverse )
!
!  CONEX3
!
  name = 'CONEX3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  call conex3 ( n, a )
  call conex3_inverse ( n, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  write ( *, '(a)' ) ' '
  do i = 1, 3
    m = m_test(i)
    call condition_sample1 ( n, a, m, cond )
    write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) name, m, n, cond_l1, cond
  end do
  deallocate ( a )
  deallocate ( a_inverse )
!
!  CONEX4
!
  name = 'CONEX4'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  call conex4 ( a )
  call conex4_inverse ( a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  write ( *, '(a)' ) ' '
  do i = 1, 3
    m = m_test(i)
    call condition_sample1 ( n, a, m, cond )
    write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) name, m, n, cond_l1, cond
  end do
  deallocate ( a )
  deallocate ( a_inverse )
!
!  KAHAN
!
  name = 'KAHAN'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  alpha = 0.25D+00
  call kahan ( alpha, n, n, a )
  call kahan_inverse ( alpha, n, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  write ( *, '(a)' ) ' '
  do i = 1, 3
    m = m_test(i)
    call condition_sample1 ( n, a, m, cond )
    write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) name, m, n, cond_l1, cond
  end do
  deallocate ( a )
  deallocate ( a_inverse )
!
!  Random
!
  seed = 123456789

  do j = 1, 5
    name = 'RANDOM'
    n = 4
    allocate ( a(1:n,1:n) )
    allocate ( a_inverse(1:n,1:n) )
    allocate ( pivot(1:n) )
    call r8mat_uniform_01 ( n, n, seed, a )
    a_inverse(1:n,1:n) = a(1:n,1:n)
    call r8ge_fa ( n, a_inverse, pivot, info )
    call r8ge_inverse ( n, a_inverse, pivot )
    a_norm_l1         = r8mat_norm_l1 ( n, n, a )
    a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
    cond_l1           = a_norm_l1 * a_inverse_norm_l1
    write ( *, '(a)' ) ' '
    do i = 1, 3
      m = m_test(i)
      call condition_sample1 ( n, a, m, cond )
      write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) name, m, n, cond_l1, cond
    end do
    deallocate ( a )
    deallocate ( a_inverse )
    deallocate ( pivot )
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests CONDITION_HAGER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: a_inverse(:,:)
  real ( kind = 8 ) a_inverse_norm_l1
  real ( kind = 8 ) a_norm_l1
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) cond
  real ( kind = 8 ) cond_l1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) n
  character ( len = 80 ) name
  integer ( kind = 4 ), allocatable :: pivot(:)
  real ( kind = 8 ) r8mat_norm_l1
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  CONDITION_HAGER estimates the L1 condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix               Order   Condition         Hager'
  write ( *, '(a)' ) ' '
!
!  Combinatorial matrix.
!
  name = 'Combinatorial'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  alpha = 2.0D+00
  beta = 3.0D+00
  call combin ( alpha, beta, n, a )
  call combin_inverse ( alpha, beta, n, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  call condition_hager ( n, a, cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
  deallocate ( a )
  deallocate ( a_inverse )
!
!  CONEX1
!
  name = 'CONEX1'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  alpha = 100.0D+00
  call conex1 ( alpha, a )
  call conex1_inverse ( alpha, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  call condition_hager ( n, a, cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
  deallocate ( a )
  deallocate ( a_inverse )
!
!  CONEX2
!
  name = 'CONEX2'
  n = 3
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  alpha = 100.0D+00
  call conex2 ( alpha, a )
  call conex2_inverse ( alpha, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  call condition_hager ( n, a, cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
  deallocate ( a )
  deallocate ( a_inverse )
!
!  CONEX3
!
  name = 'CONEX3'
  n = 5
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  call conex3 ( n, a )
  call conex3_inverse ( n, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  call condition_hager ( n, a, cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
  deallocate ( a )
  deallocate ( a_inverse )
!
!  CONEX4
!
  name = 'CONEX4'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  call conex4 ( a )
  call conex4_inverse ( a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  call condition_hager ( n, a, cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
  deallocate ( a )
  deallocate ( a_inverse )
!
!  KAHAN
!
  name = 'KAHAN'
  n = 4
  allocate ( a(1:n,1:n) )
  allocate ( a_inverse(1:n,1:n) )
  alpha = 0.25D+00
  call kahan ( alpha, n, n, a )
  call kahan_inverse ( alpha, n, a_inverse )
  a_norm_l1         = r8mat_norm_l1 ( n, n, a )
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
  cond_l1           = a_norm_l1 * a_inverse_norm_l1
  call condition_hager ( n, a, cond )
  write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
  deallocate ( a )
  deallocate ( a_inverse )
!
!  Random
!
  seed = 123456789

  do i = 1, 5
    name = 'RANDOM'
    n = 4
    allocate ( a(1:n,1:n) )
    allocate ( a_inverse(1:n,1:n) )
    allocate ( pivot(1:n) )
    call r8mat_uniform_01 ( n, n, seed, a )
    a_inverse(1:n,1:n) = a(1:n,1:n)
    call r8ge_fa ( n, a_inverse, pivot, info )
    call r8ge_inverse ( n, a_inverse, pivot )
    a_norm_l1         = r8mat_norm_l1 ( n, n, a )
    a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
    cond_l1           = a_norm_l1 * a_inverse_norm_l1
    call condition_hager ( n, a, cond )
    write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) name, n, cond_l1, cond
    deallocate ( a )
    deallocate ( a_inverse )
    deallocate ( pivot )
  end do

  return
end
