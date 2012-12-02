program main

!*****************************************************************************80
!
!! MAIN is the main program for QR_SOLVE_PRB.
!
!  Discussion:
!
!    QR_SOLVE_PRB tests the QR_SOLVE library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QR_SOLVE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the QR_SOLVE library.'
  write ( *, '(a)' ) '  QR_SOLVE needs the R8LIB library.'
  write ( *, '(a)' ) '  This test also needs the TEST_LS library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QR_SOLVE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests NORMAL_SOLVE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ) b_norm
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ), allocatable :: r1(:)
  real ( kind = 8 ) r1_norm
  real ( kind = 8 ), allocatable :: r2(:)
  real ( kind = 8 ) r2_norm
  real ( kind = 8 ) r8vec_norm
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ) x_diff_norm
  real ( kind = 8 ), allocatable :: x1(:)
  real ( kind = 8 ) x1_norm
  real ( kind = 8 ), allocatable :: x2(:)
  real ( kind = 8 ) x2_norm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  NORMAL_SOLVE is a function with a simple interface which'
  write ( *, '(a)' ) '  solves a linear system A*x = b in the least squares sense.'
  write ( *, '(a)' ) '  Compare a tabulated solution X1 to the NORMAL_SOLVE result X2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NORMAL_SOLVE cannot be applied when N < M,'
  write ( *, '(a)' ) '  or if the matrix does not have full column rank.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of problems = ', prob_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index     M     N     ||B||         ||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||'
  write ( *, '(a)' ) ' '

  do prob = 1, prob_num
!
!  Get problem size.
!
    call p00_m ( prob, m )
    call p00_n ( prob, n )
!
!  Allocate space.
!
    allocate ( a(1:m,1:n) )
    allocate ( b(1:m) )
    allocate ( r1(1:m) )
    allocate ( r2(1:m) )
    allocate ( x1(1:n) )
    allocate ( x2(1:n) )
!
!  Retrieve problem data.
!
    call p00_a ( prob, m, n, a )
    call p00_b ( prob, m, b )
    call p00_x ( prob, n, x1 )

    b_norm = r8vec_norm ( m, b )
    x1_norm = r8vec_norm ( n, x1 )
    r1(1:m) = matmul ( a(1:m,1:n), x1(1:n) ) - b(1:m)
    r1_norm = r8vec_norm ( m, r1 )
!
!  Use NORMAL_SOLVE on the problem.
!
    call normal_solve ( m, n, a, b, x2, flag )

    if ( flag /= 0 ) then

      write ( *, '(2x,i5,2x,i4,2x,i4,2x,g12.4,2x,a,2x,g12.4,2x,a,2x,g12.4,2x,a)' ) &
        prob, m, n, b_norm, '------------', x1_norm, '------------', r1_norm, '------------'

    else

      x2_norm = r8vec_norm ( n, x2 )
      r2(1:m) = matmul ( a(1:m,1:n), x2(1:n) ) - b(1:m)
      r2_norm = r8vec_norm ( m, r2 )
!
!  Compare tabulated and computed solutions.
!
      x_diff_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report results for this problem.
!
      write ( *, '(2x,i5,2x,i4,2x,i4,2x,g12.4,2x,g12.4,2x,g12.4,2x,g12.4,2x,g12.4,2x,g12.4)' ) &
        prob, m, n, b_norm, x_diff_norm, x1_norm, x2_norm, r1_norm, r2_norm

    end if
!
!  Deallocate memory.
!
    deallocate ( a )
    deallocate ( b )
    deallocate ( r1 )
    deallocate ( r2 )
    deallocate ( x1 )
    deallocate ( x2 )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests QR_SOLVE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ) b_norm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ), allocatable :: r1(:)
  real ( kind = 8 ) r1_norm
  real ( kind = 8 ), allocatable :: r2(:)
  real ( kind = 8 ) r2_norm
  real ( kind = 8 ) r8vec_norm
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ) x_diff_norm
  real ( kind = 8 ), allocatable :: x1(:)
  real ( kind = 8 ) x1_norm
  real ( kind = 8 ), allocatable :: x2(:)
  real ( kind = 8 ) x2_norm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  QR_SOLVE is a function with a simple interface which'
  write ( *, '(a)' ) '  solves a linear system A*x = b in the least squares sense.'
  write ( *, '(a)' ) '  Compare a tabulated solution X1 to the QR_SOLVE result X2.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of problems = ', prob_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index     M     N     ||B||         ||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||'
  write ( *, '(a)' ) ' '

  do prob = 1, prob_num
!
!  Get problem size.
!
    call p00_m ( prob, m )
    call p00_n ( prob, n )
!
!  Allocate space.
!
    allocate ( a(1:m,1:n) )
    allocate ( b(1:m) )
    allocate ( r1(1:m) )
    allocate ( r2(1:m) )
    allocate ( x1(1:n) )
    allocate ( x2(1:n) )
!
!  Retrieve problem data.
!
    call p00_a ( prob, m, n, a )
    call p00_b ( prob, m, b )
    call p00_x ( prob, n, x1 )

    b_norm = r8vec_norm ( m, b )
    x1_norm = r8vec_norm ( n, x1 )
    r1(1:m) = matmul ( a(1:m,1:n), x1(1:n) ) - b(1:m)
    r1_norm = r8vec_norm ( m, r1 )
!
!  Use QR_SOLVE on the problem.
!
    call qr_solve ( m, n, a, b, x2 )

    x2_norm = r8vec_norm ( n, x2 )
    r2(1:m) = matmul ( a(1:m,1:n), x2(1:n) ) - b(1:m)
    r2_norm = r8vec_norm ( m, r2 )
!
!  Compare tabulated and computed solutions.
!
    x_diff_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report results for this problem.
!
    write ( *, '(2x,i5,2x,i4,2x,i4,2x,g12.4,2x,g12.4,2x,g12.4,2x,g12.4,2x,g12.4,2x,g12.4)' ) &
      prob, m, n, b_norm, x_diff_norm, x1_norm, x2_norm, r1_norm, r2_norm
!
!  Deallocate memory.
!
    deallocate ( a )
    deallocate ( b )
    deallocate ( r1 )
    deallocate ( r2 )
    deallocate ( x1 )
    deallocate ( x2 )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests SVD_SOLVE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ) b_norm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  real ( kind = 8 ), allocatable :: r1(:)
  real ( kind = 8 ) r1_norm
  real ( kind = 8 ), allocatable :: r2(:)
  real ( kind = 8 ) r2_norm
  real ( kind = 8 ) r8vec_norm
  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ) x_diff_norm
  real ( kind = 8 ), allocatable :: x1(:)
  real ( kind = 8 ) x1_norm
  real ( kind = 8 ), allocatable :: x2(:)
  real ( kind = 8 ) x2_norm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  SVD_SOLVE is a function with a simple interface which'
  write ( *, '(a)' ) '  solves a linear system A*x = b in the least squares sense.'
  write ( *, '(a)' ) '  Compare a tabulated solution X1 to the SVD_SOLVE result X2.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of problems = ', prob_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index     M     N     ||B||         ||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||'
  write ( *, '(a)' ) ' '

  do prob = 1, prob_num
!
!  Get problem size.
!
    call p00_m ( prob, m )
    call p00_n ( prob, n )
!
!  Allocate space.
!
    allocate ( a(1:m,1:n) )
    allocate ( b(1:m) )
    allocate ( r1(1:m) )
    allocate ( r2(1:m) )
    allocate ( x1(1:n) )
    allocate ( x2(1:n) )
!
!  Retrieve problem data.
!
    call p00_a ( prob, m, n, a )
    call p00_b ( prob, m, b )
    call p00_x ( prob, n, x1 )

    b_norm = r8vec_norm ( m, b )
    x1_norm = r8vec_norm ( n, x1 )
    r1(1:m) = matmul ( a(1:m,1:n), x1(1:n) ) - b(1:m)
    r1_norm = r8vec_norm ( m, r1 )
!
!  Use SVD_SOLVE on the problem.
!
    call svd_solve ( m, n, a, b, x2 )

    x2_norm = r8vec_norm ( n, x2 )
    r2(1:m) = matmul ( a(1:m,1:n), x2(1:n) ) - b(1:m)
    r2_norm = r8vec_norm ( m, r2 )
!
!  Compare tabulated and computed solutions.
!
    x_diff_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report results for this problem.
!
    write ( *, '(2x,i5,2x,i4,2x,i4,2x,g12.4,2x,g12.4,2x,g12.4,2x,g12.4,2x,g12.4,2x,g12.4)' ) &
      prob, m, n, b_norm, x_diff_norm, x1_norm, x2_norm, r1_norm, r2_norm
!
!  Deallocate memory.
!
    deallocate ( a )
    deallocate ( b )
    deallocate ( r1 )
    deallocate ( r2 )
    deallocate ( x1 )
    deallocate ( x2 )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests DQRLS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) itask
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jpvt(n)
  integer ( kind = 4 ) kr
  real ( kind = 8 ) qraux(n)
  real ( kind = 8 ) tol
  real ( kind = 8 ) work(n)
  real ( kind = 8 ) x(n)
!
!  Set up least-squares problem
!  quadratic model, equally-spaced points
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  DQRLS solves a linear system A*x = b in the least squares sense.'

  do i = 1, m
    a(i,1) = 1.0D+00
    do j = 2, n
      a(i,j) = a(i,j-1) * real ( i, kind = 8 )
    end do
  end do

  b(1:5) = (/ 1.0D+00, 2.3D+00, 4.6D+00, 3.1D+00, 1.2D+00 /)

  tol = 1.0D-06

  call r8mat_print ( m, n, a, '  Coefficient matrix A:' )

  call r8vec_print ( m, b, '  Right hand side b:' )
!
!  Solve least-squares problem
!
  itask = 1
  call dqrls ( a, m, m, n, tol, kr, b, x, b, jpvt, qraux, work, itask, ind )
!
!  Print results
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Error code =', ind
  write ( *, '(a,i4)' ) '  Estimated matrix rank =', kr

  call r8vec_print ( n, x, '  Least squares solution x:' )

  call r8vec_print ( m, b, '  Residuals A*x-b' )

  return
end
