program main

!*****************************************************************************80
!
!! REGRESSION_PRB runs the REGRESSION tests.
!
!  Modified:
!
!    13 January 2011
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REGRESSION_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the REGRESSION library.'

  call test01 ( 'x01.txt' )
  call test01 ( 'x03.txt' )
  call test01 ( 'x60.txt' )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( 'x06.txt' )
  call test08 ( 'x06.txt' )
  call test09 ( 'x03.txt' )
!
!  Some L2 solvers.
!
  call test10 ( 'x03.txt' )
  call test11 ( 'x03.txt' )
  call test12 ( 'x03.txt' )
  call test13 ( 'x03.txt' )
  call test14 ( 'x03.txt' )
  call test15 ( 'x03.txt' )
!
!  Some L1 solvers.
!
  call test16 ( 'x10.txt' )
  call test17 ( 'x10.txt' )
  call test18 ( 'x10.txt' )
!
!  Some LI solvers.
!
  call test19 ( 'x03.txt' )
  call test20 ( 'x03.txt' )
  call test21 ( 'x03.txt' )
!
!  Orthogonal regression.
!
  call test30 ( 'x02.txt' )
  call test305 ( 'x02.txt' )
  call test306 ( 'x02.txt' )
  call test307 ( 'x02.txt' )
!
!  Minimize the norm of A*X-B, subject to linear equalities and inequalities.
!
  call test33 ( 'x54.txt' )
  call test33 ( 'x61.txt' )
  call test34 ( 'x54.txt' )
  call test35 ( 'x54.txt' )
!
!  Other stuff.
!
  call test22 ( 'x02.txt' )
  call test225 ( 'x54.txt' )

  call test23 ( 'x10.txt' )
  call test24 ( 'x43.txt' )
  call test25 ( 'x04.txt' )
  call test26 ( 'x54_03.txt' )
  call test27 ( 'x03.txt' )
  call test28 ( 'x03.txt' )
  call test29 ( 'x03.txt' )
  call test31 ( 'x03.txt' )
  call test32 ( 'x03.txt' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REGRESSION_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( file_name )

!*****************************************************************************80
!
!! TEST01 tests EXAMPLE_READ, EXAMPLE_PRINT.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  character ( len = * ) file_name
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  EXAMPLE_SIZE gets the size of an example file.'
  write ( *, '(a)' ) '  EXAMPLE_READ reads an example file.'
  write ( *, '(a)' ) '  EXAMPLE_PRINT prints an example file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open and print file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  allocate ( a(m,n) )
  allocate ( b(m) )
  allocate ( a_title(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )

  call example_print ( m, n, a, b, a_title, b_title )

  deallocate ( a )
  deallocate ( b )
  deallocate ( a_title )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests GEN, EXAMPLE_PRINT.
!
  implicit none

  integer ( kind = 4 ), parameter :: ntest = 3

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iseed
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter, dimension ( ntest ) :: m_test = (/ 5,8, 4 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter, dimension ( ntest ) :: n_test = (/ 3, 4, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  GEN generates a random example.'

  iseed = 1483

  do i = 1, ntest

    m = m_test(i)
    n = n_test(i)

    allocate ( a(1:m,1:n) )
    allocate ( a_title(1:n) )
    allocate ( b(1:m) )

    call gen ( m, n, a, b, iseed )

    do j = 1, n
      a_title(j) = 'A'
      write ( a_title(j)(2:2), '(i1)' ) j
    end do

    b_title = 'B'

    call example_print ( m, n, a, b, a_title, b_title )

    deallocate ( a )
    deallocate ( a_title )
    deallocate ( b )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests SCR.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nu = 3

  real ( kind = 4 ) a(m,n)
  real ( kind = 4 ) aa(m,nu)
  logical first
  integer ( kind = 4 ) kbit(n)
  integer ( kind = 4 ) na
  integer ( kind = 4 ), parameter :: nl = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  SCR returns all the M by NA submatrices of an'
  write ( *, '(a)' ) '  M by N matrix A, where NA is between NL and NU.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For our problem:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  M =  ', m
  write ( *, '(a,i6)' ) '  N =  ', n
  write ( *, '(a,i6)' ) '  NL = ', nl
  write ( *, '(a,i6)' ) '  NU = ', nu

  call r4mat_indicator ( m, n, a )

  call r4mat_print ( m, n, a, '  The M by N matrix A:' )

  kbit(1:n) = 0
  first = .true.

  do

    call scr ( a, m, n, nl, nu, first, na, aa, kbit )

    if ( first ) then
      exit
    end if

    call r4mat_print ( m, na, aa, ' ' )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests QRBD.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 4 ) bd(n,n)
  real ( kind = 4 ) e(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  real ( kind = 4 ) q(n)
  real ( kind = 4 ) s(n,n)
  real ( kind = 4 ) u(n,n)
  real ( kind = 4 ) v(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  QRBD computes the singular values S of a bidiagonal'
  write ( *, '(a)' ) '    matrix BD, and can also compute the decomposition'
  write ( *, '(a)' ) '    factors U and V, so that'
  write ( *, '(a)' ) '      S = U * BD * V.'

  call random_number ( harvest = q(1:n) )
  call random_number ( harvest = e(2:n) )

  bd(1:n,1:n) = 0.0E+00
  do i = 1, n
    bd(i,i) = q(i)
  end do
  do i = 2, n
    bd(i-1,i) = e(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The bidiagonal matrix BD:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(5f12.4)' ) bd(i,1:n)
  end do

  u(1:n,1:n) = 0.0E+00
  v(1:n,1:n) = 0.0E+00
  do i = 1, n
    u(i,i) = 1.0E+00
    v(i,i) = 1.0E+00
  end do

  call qrbd ( q, e, n, v, n, n, u, n, n,iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04 - Warning!'
    write ( *, '(a,i6)' ) '  QRBD returned IFLAG = ', iflag
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The singular values of BD:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, q(i)
  end do

  s(1:n,1:n) = 0.0E+00
  do i = 1, n
    s(i,i) = q(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The factor U:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(5f12.4)' ) u(i,1:n)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The factor V:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(5f12.4)' ) v(i,1:n)
  end do

  bd(1:n,1:n) = matmul ( matmul ( transpose ( u(1:n,1:n) ), s(1:n,1:n) ), &
    transpose ( v(1:n,1:n) ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product U'' * S * V'' = BD:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(5f12.4)' ) bd(i,1:n)
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests SVD.
!
!  Discussion:
!
!    In our special example, the matrix is square and symmetric.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: m = n
  integer ( kind = 4 ), parameter :: nm = n

  real ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical matu
  logical matv
  real ( kind = 4 ) r(m,n)
  real ( kind = 4 ) s(n)
  real ( kind = 4 ) u(m,n)
  real ( kind = 4 ) v(n,n)
!
!  Set the values of the matrix.
!
  a(1,1) = 0.9900E+00
  a(1,2) = 0.0020E+00
  a(1,3) = 0.0060E+00
  a(1,4) = 0.0020E+00

  a(2,1) = 0.0020E+00
  a(2,2) = 0.9900E+00
  a(2,3) = 0.0020E+00
  a(2,4) = 0.0060E+00

  a(3,1) = 0.0060E+00
  a(3,2) = 0.0020E+00
  a(3,3) = 0.9900E+00
  a(3,4) = 0.0020E+00

  a(4,1) = 0.0020E+00
  a(4,2) = 0.0060E+00
  a(4,3) = 0.0020E+00
  a(4,4) = 0.9900E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  SVD computes the singular value decomposition'
  write ( *, '(a)' ) '  of a real general matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Matrix order = ', n

  call r4mat_print ( m, n, a, '  The matrix A:' )

  matu = .true.
  matv = .true.

  call svd ( m, n, a, s, matu, u, matv, v, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Warning!'
    write ( *, '(a,i6)' ) '  SVD returned error flag IFLAG = ', iflag
  end if

  call r4vec_print ( n, s, '  The singular values S' )

  call r4mat_print ( m, n, u, '  The U matrix:' )

  call r4mat_print ( n, n, v, '  The V matrix:' )

  do j = 1, n
    v(1:n,j) = s(j) * v(1:n,j)
  end do

  r(1:m,1:n) = matmul ( u(1:m,1:n), transpose ( v(1:n,1:n) ) )

  call r4mat_print ( m, n, r, '  The product U * S * Transpose(V):' )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests C01M.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 8

  integer ( kind = 4 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  C01M generates all subsets of size M from a'
  write ( *, '(a)' ) '  set of size N, one at a time, trying to use'
  write ( *, '(a)' ) '  simple exchanges.'
  write ( *, '(a)' ) ' '

  m = 2
  d(1:n) = 0
  i = 0

  do

    call c01m ( n, m, d )

    i = i + 1
    write ( *, '(i3,5x,10i2)' ) i, d(1:n)
!
!  We're done when the first M things are 1.
!
    if ( sum ( d(1:m) ) == m ) then
      exit
    end if

  end do

  return
end
subroutine test07 ( file_name )

!*****************************************************************************80
!
!! TEST07 tests CWLR_L2.
!
!  Discussion:
!
!    Setting up values of ML and MU, and obtaining a partition that
!    satisfies the appropriate constraints can be a tedious task.
!
  implicit none

  integer ( kind = 4 ), parameter :: s = 4

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ) d
  real ( kind = 4 ), allocatable, dimension ( : ) :: e
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: x
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  CWLR_L2 uses clustering techniques.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  ml = n + 1
  mu = m

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations, M =   ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables, N =   ', n
  write ( *, '(a,i6)' ) '  Number of clusters,          S =   ', s
  write ( *, '(a,i6)' ) '  Minimum cluster population,  ML =  ', ml
  write ( *, '(a,i6)' ) '  Maximum cluster population,  MU =  ', mu

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( e(s) )
  allocate ( r(m) )
  allocate ( x(s,n) )
  allocate ( z(m) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Compute a random initial partition with occupancy limits.
!
  call random_partition2 ( m, s, ml, mu, z )

  call cwlr_l2 ( a, m, n, b, s, eps, z, ml, mu, x, e, d, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST07 - Warning!'
    write ( *, '(a,i6)' ) '  CWLR_L2 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call r4mat_print ( s, n, x, '  Solution column vectors X:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Objective function = ', d

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( e )
  deallocate ( r )
  deallocate ( x )
  deallocate ( z )

  return
end
subroutine test08 ( file_name )

!*****************************************************************************80
!
!! TEST08 tests CWLR_LI.
!
  implicit none

  integer ( kind = 4 ), parameter :: s = 4

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ) d
  real ( kind = 4 ), allocatable, dimension ( : ) :: e
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: x
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  CWLR_LI uses clustering techniques.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  ml = n + 1
  mu = m

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations, M =   ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables, N =   ', n
  write ( *, '(a,i6)' ) '  Number of clusters,          S =   ', s
  write ( *, '(a,i6)' ) '  Minimum cluster population,  ML =  ', ml
  write ( *, '(a,i6)' ) '  Maximum cluster population,  MU =  ', mu

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( e(s) )
  allocate ( r(m) )
  allocate ( x(s,n) )
  allocate ( z(m) )

  write ( *, * ) 'DEBUG: Call example_read'
  call example_read ( file_name, m, n, a, b, a_title, b_title )
  write ( *, * ) 'DEBUG: Returned from example_read'
!
!  Compute a random initial partition with occupancy limits.
!
  write ( *, * ) 'DEBUG: Call random_partition2'
  call random_partition2 ( m, s, ml, mu, z )
  write ( *, * ) 'DEBUG: Call cwlr_li'
  write ( *, * ) ' '
  write ( *, * ) '  ERRORS ARE OCCURRING IN CALL To CWKL_LI.'
  write ( *, * ) '  WE WILL SKIP THIS CALL FOR NOW'

  if ( .false. ) then

  call cwlr_li ( a, m, n, b, s, eps, z, ml, mu, x, e, d, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Warning!'
    write ( *, '(a,i6)' ) '  CWLR_LI returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call r4mat_print ( s, n, x, '  Solution column vectors X:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Objective function = ', d

  end if

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( e )
  deallocate ( r )
  deallocate ( x )
  deallocate ( z )

  return
end
subroutine test09 ( file_name )

!*****************************************************************************80
!
!! TEST09 tests REGR_LP.
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 8

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ) cond
  real ( kind = 4 ), parameter :: eps1 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps2 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps3 = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), parameter :: itmax = 30
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ) p
  real ( kind = 4 ), dimension ( test_num ) :: p_test = (/ &
    1.1E+00, 1.2E+00, 1.4E+00, 1.7E+00, 2.0E+00, 2.5E+00, 4.0E+00, 6.5E+00 /)
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) rank
  real ( kind = 4 ) r4vec_norm_lp
  integer ( kind = 4 ) test
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  REGR_LP minimizes the LP norm of the residual'
  write ( *, '(a)' ) '  A*X-B, where 1 < P < Infinity.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )

  do test = 1, test_num

    p = p_test(test)
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Using P = ', p
!
!  Solve the linear system.
!
    call regr_lp ( a, m, n, b, p, eps1, eps2, eps3, itmax, x, iflag )

    if ( iflag /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09 - Warning!'
      write ( *, '(a,i6)' ) '  REGR_LP returned IFLAG = ', iflag
      cycle
    end if
!
!  Print the solution and residual.
!
!  The routine destroys the matrix A, so we need to recover the
!  original data in order to compute the residual.
!
    call example_read ( file_name, m, n, a, b, a_title, b_title )

    call r4vec_print ( n, x, '  Solution:' )

    r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  LP norm of residual = ', &
      r4vec_norm_lp ( m, r, p )

  end do

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test10 ( file_name )

!*****************************************************************************80
!
!! TEST10 tests NORMAL_L2.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  real ( kind = 4 ), parameter :: lambda = 0.0E+00
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  NORMAL_L2 solves the normal equations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call normal_l2 ( a, m, n, b, lambda, x )
!
!  Print the solution and residual.
!
  call r4vec_print ( n, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of residual = ', r4vec_norm_l2 ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test11 ( file_name )

!*****************************************************************************80
!
!! TEST11 tests MGS_L2.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  logical bnew
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  MGS_L2 uses the modified Gram Schmidt method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  bnew = .false.

  call mgs_l2 ( a, m, n, b, eps, bnew, x, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Warning!'
    write ( *, '(a,i6)' ) '  MGS_L2 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
!  The MGS_L2 algorithm destroys the information in A and B, so
!  we have to recover it from the example file.
!
  call example_read ( file_name, m, n, a, b, a_title, b_title )

  call r4vec_print ( n, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of residual = ', r4vec_norm_l2 ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test12 ( file_name )

!*****************************************************************************80
!
!! TEST12 tests ICMGS_L2.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a2
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  logical bnew
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  ICMGS_L2 uses the modified Gram Schmidt method.'
  write ( *, '(a)' ) '  Note that this code should get the same answer'
  write ( *, '(a)' ) '  as MGS_L2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n
!
!  Note that X is of size N+1 for this problem!
!
  allocate ( a(m,n) )
  allocate ( a2(m,n-1) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Our input data matrix A includes an initial column of 1's.
!  In order to test ICMGS_L2, we need to make a copy of A without
!  this initial column.
!
  n2 = n - 1
  a2(1:m,1:n2) = a(1:m,2:n)
!
!  Solve the linear system.
!
  call icmgs_l2 ( a2, m, n2, b, eps, x, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST12 - Warning!'
    write ( *, '(a,i6)' ) '  ICMGS_L2 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!  We still have our old matrix A around, because ICMGS_L2 destroyed A2,
!  not A.  So we don't need to read A back in again to compute
!  the residual.
!
  call r4vec_print ( n, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of residual = ', r4vec_norm_l2 ( m, r )

  deallocate ( a )
  deallocate ( a2 )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test13 ( file_name )

!*****************************************************************************80
!
!! TEST13 tests GIVR_L2.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  GIVR_L2 uses the Givens Rotation method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call givr_l2 ( a, m, n, b, eps, x, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST13 - Warning!'
    write ( *, '(a,i6)' ) '  GIVR_L2 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
!  Because GIVR_L2 modifies the matrix A, we need to recover the
!  original data in order to compute the residual.
!
  call example_read ( file_name, m, n, a, b, a_title, b_title )

  call r4vec_print ( n, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of residual = ', r4vec_norm_l2 ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test14 ( file_name )

!*****************************************************************************80
!
!! TEST14 tests HFTI_L2.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  HFTI_L2 uses Householder transformations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call hfti_l2 ( a, m, n, b, eps, x )
!
!  Print the solution and residual.
!
!  HFTI_L2 destroys the matrix A, so we need to recover the
!  original data in order to compute the residual.
!
  call example_read ( file_name, m, n, a, b, a_title, b_title )

  call r4vec_print ( n, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of residual = ', r4vec_norm_l2 ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test15 ( file_name )

!*****************************************************************************80
!
!! TEST15 tests SVDR_L2.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ) cond
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) rank
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  SVDR_L2 uses the singular value decomposition to'
  write ( *, '(a)' ) '  solve a least squares problem.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call svdr_l2 ( a, m, n, b, eps, x, cond, rank, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Warning!'
    write ( *, '(a,i6)' ) '  SVDR_L2 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
!  SVDR_L2 destroys the matrix A, so we need to recover the
!  original data in order to compute the residual.
!
  call example_read ( file_name, m, n, a, b, a_title, b_title )

  call r4vec_print ( n, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of residual = ', r4vec_norm_l2 ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test16 ( file_name )

!*****************************************************************************80
!
!! TEST16 tests A478_L1.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) rank
  real ( kind = 4 ) r4vec_norm_l1
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  A478_L1 minimizes the L1 norm of the residual'
  write ( *, '(a)' ) '  A*X-B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call a478_l1 ( a, m, n, b, eps, rank, x, r, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST16 - Warning!'
    write ( *, '(a,i6)' ) '  A478_L1 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call r4vec_print ( n, x, '  Solution:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Numerical estimate of rank = ', rank
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L1 norm of computed residual =   ', &
    r4vec_norm_l1 ( m, r )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a,g14.6)' ) '  L1 norm of recomputed residual = ', &
    r4vec_norm_l1 ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test17 ( file_name )

!*****************************************************************************80
!
!! TEST17 tests AFK_L1.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l1
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  AFK_L1 minimizes the L1 norm of the residual'
  write ( *, '(a)' ) '  A*X-B.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call afk_l1 ( a, m, n, b, eps, x, r, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Warning!'
    write ( *, '(a,i6)' ) '  AFK_L1 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call example_read ( file_name, m, n, a, b, a_title, b_title )

  call r4vec_print ( n, x, '  Solution:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L1 norm of computed residual =   ', &
    r4vec_norm_l1 ( m, r )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a,g14.6)' ) '  L1 norm of recomputed residual = ', &
    r4vec_norm_l1 ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test18 ( file_name )

!*****************************************************************************80
!
!! TEST18 tests BLOD_L1.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) r4vec_norm_l1
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  BLOD_L1 minimizes the L1 norm of the residual'
  write ( *, '(a)' ) '  A*X-B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call blod_l1 ( a, m, n, b, x, res_norm )
!
!  Print the solution and residual.
!
  call example_read ( file_name, m, n, a, b, a_title, b_title )

  call r4vec_print ( n, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Computed L1 norm of residual =   ', res_norm
  write ( *, '(a,g14.6)' ) '  Recomputed L1 norm of residual = ', &
    r4vec_norm_l1 ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test19 ( file_name )

!*****************************************************************************80
!
!! TEST19 tests A328_LI.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) r4vec_norm_li
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  A328_LI minimizes the L-infinity norm of the residual'
  write ( *, '(a)' ) '  A*X-B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call a328_li ( a, m, n, b, x, res_norm, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST19 - Warning!'
    write ( *, '(a,i6)' ) '  A328_LI returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call example_read ( file_name, m, n, a, b, a_title, b_title )

  call r4vec_print ( n, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Computed L-infinity norm of residual =   ', &
    res_norm
  write ( *, '(a,g14.6)' ) '  Recomputed L-infinity norm of residual = ', &
    r4vec_norm_li ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test20 ( file_name )

!*****************************************************************************80
!
!! TEST20 tests A495_LI.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps1 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps2 = 0.0E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) rank
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) r4vec_norm_li
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  A495_LI minimizes the L-infinity norm of the residual'
  write ( *, '(a)' ) '  A*X-B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call a495_li ( a, m, n, b, eps1, eps2, x, rank, res_norm, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST20 - Warning!'
    write ( *, '(a,i6)' ) '  A495_LI returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call example_read ( file_name, m, n, a, b, a_title, b_title )

  call r4vec_print ( n, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Computed L-infinity norm of residual =   ', &
    res_norm
  write ( *, '(a,g14.6)' ) '  Recomputed L-infinity norm of residual = ', &
    r4vec_norm_li ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test21 ( file_name )

!*****************************************************************************80
!
!! TEST21 tests ABD_LI.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) rank
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) r4vec_norm_li
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  ABD_LI minimizes the L-infinity norm of the residual'
  write ( *, '(a)' ) '  A*X-B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call abd_li ( a, m, n, b, eps, rank, x, r, res_norm, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST21 - Warning!'
    write ( *, '(a,i6)' ) '  ABD_LI returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call example_read ( file_name, m, n, a, b, a_title, b_title )

  call r4vec_print ( n, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Computed L-infinity norm of residual =   ', &
    res_norm
  write ( *, '(a,g14.6)' ) '  Recomputed L-infinity norm of residual = ', &
    r4vec_norm_li ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test22 ( file_name )

!*****************************************************************************80
!
!! TEST22 tests NN_L2.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  NN_L2 minimizes the L2 norm of the residual'
  write ( *, '(a)' ) '  A*X-B for nonnegative X.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call nn_l2 ( a, m, n, b, x, res_norm, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST22 - Warning!'
    write ( *, '(a,i6)' ) '  NN_L2 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call example_read ( file_name, m, n, a, b, a_title, b_title )

  call r4vec_print ( n, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Computed L2 norm of residual =   ', res_norm
  write ( *, '(a,g14.6)' ) '  Recomputed L2 norm of residual = ', &
    r4vec_norm_l2 ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test225 ( file_name )

!*****************************************************************************80
!
!! TEST225 tests NN_L1.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: c
  real ( kind = 4 ), allocatable, dimension ( : ) :: d
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: e
  real ( kind = 4 ), parameter :: eps = 0.00001E+00
  character ( len = * ) file_name
  character ( len = 80 ) file_name2
  real ( kind = 4 ), allocatable, dimension ( : ) :: h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), parameter :: itmax = 15
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r1
  real ( kind = 4 ), allocatable, dimension ( : ) :: r2
  real ( kind = 4 ), allocatable, dimension ( : ) :: r3
  real ( kind = 4 ), allocatable, dimension ( : ) :: res
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) r4vec_norm_l1
  integer ( kind = 4 ) s
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST225'
  write ( *, '(a)' ) '  NN_L1 minimizes the L1 norm of the residual'
  write ( *, '(a)' ) '  A*X-B for nonnegative X.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Find an N vector X which minimizes the residual:'
  write ( *, '(a)' ) '    || A * X - B ||'
  write ( *, '(a)' ) '  and satisifies the linear equalities:'
  write ( *, '(a)' ) '    C * X = D'
  write ( *, '(a)' ) '  and the linear inequalities:'
  write ( *, '(a)' ) '    E * X >= H'
  write ( *, '(a)' ) '  and the nonnegativity constraint:'
  write ( *, '(a)' ) '    X >= 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_multi_size ( file_name, m, n, s, file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n
  write ( *, '(a,i6)' ) '  Number of subsystems =        ', s

  allocate ( a_title(n) )
!
!  Read A and B.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name2 )

  call example_size ( file_name2, m, n )

  allocate ( a(m,n) )
  allocate ( b(m) )
  allocate ( r1(m) )

  call example_read ( file_name2, m, n, a, b, a_title, b_title )
!
!  Read C and D.
!
  call file_name_inc ( file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name2 )

  call example_size ( file_name2, l, n )

  allocate ( c(l,n) )
  allocate ( d(l) )
  allocate ( r2(l) )

  call example_read ( file_name2, l, n, c, d, a_title, b_title )
!
!  Read E and H.
!
  call file_name_inc ( file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name2 )

  call example_size ( file_name2, k, n )

  allocate ( e(k,n) )
  allocate ( h(k) )
  allocate ( r3(k) )

  call example_read ( file_name2, k, n, e, h, a_title, b_title )
!
!  Solve the linear system.
!
  allocate ( res(m+l+k) )
  allocate ( x(n) )

  call nn_l1 ( m, l, k, n, a, b, c, d, e, h, eps, itmax, x, res, &
    res_norm, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST225 - Warning!'
    write ( *, '(a,i6)' ) '  NN_L1 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call r4vec_print ( n, x, '  Solution:' )

  r1(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)
  r2(1:l) = matmul ( c(1:l,1:n), x(1:n) ) - d(1:l)
  r3(1:k) = matmul ( e(1:k,1:n), x(1:n) ) - h(1:k)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L1 norms of residuals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  NN_L1 claims A*X-B residual norm = ', res_norm
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A*X-B = ', r4vec_norm_l1 ( m, r1 )
  write ( *, '(a,g14.6)' ) '  C*X-D = ', r4vec_norm_l1 ( l, r2 )
  write ( *, '(a,g14.6)' ) '  E*X-H = ', r4vec_norm_l1 ( k, r3 )
  call r4vec_print ( k, r3, '  E*X-H:' )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( c )
  deallocate ( d )
  deallocate ( e )
  deallocate ( h )
  deallocate ( r1 )
  deallocate ( r2 )
  deallocate ( r3 )
  deallocate ( res )
  deallocate ( x )

  return
end
subroutine test23 ( file_name )

!*****************************************************************************80
!
!! TEST23 tests SCRF_L1.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ), allocatable, dimension ( : ) :: kbit
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: na = 2
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) rank
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) r4vec_norm_l1
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  SCRF_L1 minimizes the L1 norm of the residual'
  write ( *, '(a)' ) '  A*X-B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( kbit(n) )
  allocate ( r(m) )
  allocate ( x(na) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call scrf_l1 ( a, m, n, b, na, eps, x, kbit, res_norm, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST23 - Warning!'
    write ( *, '(a,i6)' ) '  SCRF_L1 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call i4vec_print ( n, kbit, '  Index vector' )
  call r4vec_print ( na, x, '  Solution:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Computed residual =   ', res_norm

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( kbit )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test24 ( file_name )

!*****************************************************************************80
!
!! TEST24 tests AVLLSQ.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps1 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps2 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps3 = 0.000001E+00
  character ( len = * ) file_name
  character ( len = 80 ) file_name2
  real ( kind = 4 ) fx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), parameter :: itmax = 15
  integer ( kind = 4 ) m
  integer ( kind = 4 ), allocatable, dimension ( : ) :: m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  integer ( kind = 4 ) s
  real ( kind = 4 ), allocatable, dimension ( : ) :: w
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  AVLLSQ carries out average linear regression.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_multi_size ( file_name, m, n, s, file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n
  write ( *, '(a,i6)' ) '  Number of clusters =          ', s
  write ( *, '(a)' ) '  Base filename = ' // trim ( file_name2 )

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( m2(s) )
  allocate ( r(m) )
  allocate ( w(s) )
  allocate ( x(n) )

  i1 = 0
  i2 = 0

  do i = 1, s

    write ( *, '(a)' ) '  Open subsystem data file ' // trim ( file_name2 )

    call example_size ( file_name2, m2(i), n2 )

    if ( n2 /= n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Fatal error!'
      write ( *, '(a)' ) '  Subsystem and system have different values of N.'
      write ( *, '(a,i6)' ) '  System N    = ', n
      write ( *, '(a,i6)' ) '  Subsystem N = ', n2
      return
    end if

    i1 = i2 + 1
    i2 = i1 + m2(i) - 1

    write ( *, '(a,i6)' ) '  Number of rows/observations = ', m2(i)

    call example_read ( file_name2, m2(i), n2, a(i1:i2,1:n), b(i1:i2), &
      a_title, b_title )

    call file_name_inc ( file_name2 )

  end do
!
!  Verify that we got everything in correctly.
!
  call example_print ( m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  w(1:s) = 1.0E+00 / real ( s )

  call avllsq ( a, m, n, b, m2, s, w, itmax, eps1, eps2, eps3, x, &
    fx, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST24 - Warning!'
    write ( *, '(a,i6)' ) '  AVLLSQ returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call r4vec_print ( n, x, '  Solution:' )

  x = (/ 0.539414E+00, 2.39305E+00 /)

  call r4vec_print ( n, x, '  Spaeth''s solution, page 190:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Value of objective function = ', fx

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( m2 )
  deallocate ( r )
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test25 ( file_name )

!*****************************************************************************80
!
!! TEST25 tests ROBUST.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps1 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps2 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps3 = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), parameter :: itmax = 15
  integer ( kind = 4 ) m
  integer ( kind = 4 ) method
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) sr
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  ROBUST minimizes an objective function based on'
  write ( *, '(a)' ) '  the residuals of A*X-B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  do method = 1, 8

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'Method = ', method

    call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
    call robust ( a, m, n, b, method, eps1, eps2, eps3, itmax, &
      x, sr, iflag )

    if ( iflag /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST25 - Warning!'
      write ( *, '(a,i6)' ) '  ROBUST returned IFLAG = ', iflag
      return
    end if
!
!  Print the solution and residual.
!
    call r4vec_print ( n, x, '  Solution:' )

    r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Value of SR =                    ', sr
    write ( *, '(a,g14.6)' ) '  Recomputed L2 norm of residual = ', &
      r4vec_norm_l2 ( m, r )

  end do

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test26 ( file_name )

!*****************************************************************************80
!
!! TEST26 tests LDP_L2.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  real ( kind = 4 ), parameter :: lambda = 0.0E+00
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ), allocatable, dimension ( : ) :: x
  real ( kind = 4 ) xnorm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  LDP_L2 solves a least distance programming problem.'
  write ( *, '(a)' ) '  Find the vector X of minimum L2 norm which'
  write ( *, '(a)' ) '  satisfies the linear inequalities E * X >= H.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )

  b(1:m) = - b(1:m)
!
!  Solve the linear system.
!
  call ldp_l2 ( a, m, n, b, x, xnorm, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST26 - Warning!'
    write ( *, '(a,i6)' ) '  LDP_L2 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call r4vec_print ( n, x, '  Solution X:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  call r4vec_print ( m, r, '  Residual E * X - H:' )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test27 ( file_name )

!*****************************************************************************80
!
!! TEST27 tests RR_L2.
!
  implicit none

  integer ( kind = 4 ), parameter :: lambda_num = 3

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  real ( kind = 4 ) lambda
  real ( kind = 4 ), dimension ( lambda_num ) :: lambda_test = &
    (/ 0.0E+00, 1.0E+00, 100.0E+00 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  RR_L2 minimizes the L2 norm of:'
  write ( *, '(a)' ) '    [ A          ] * X - [ B ]'
  write ( *, '(a)' ) '    [ LAMBDA * I ]       [ 0 ] '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m+n) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  do i = 1, lambda_num

    lambda = lambda_test(i)
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Solution for LAMBDA = ', lambda

    call rr_l2 ( a, m, n, b, eps, x, lambda, iflag )

    if ( iflag /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST27 - Warning!'
      write ( *, '(a,i6)' ) '  RR_L2 returned IFLAG = ', iflag
      return
    end if
!
!  Print the solution and residual.
!
    call r4vec_print ( n, x, '  Solution:' )

    r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)
    r(m+1:m+n) = lambda * x(1:n)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  L2 norm of AX-B =     ', &
      r4vec_norm_l2 ( m, r(1:m) )
    write ( *, '(a,g14.6)' ) '  L2 norm of LAMBDA*X = ', &
      r4vec_norm_l2 ( n, r(m+1:m+n) )
    write ( *, '(a,g14.6)' ) '  L2 norm of residual = ', &
      r4vec_norm_l2 ( m+n, r )

  end do

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test28 ( file_name )

!*****************************************************************************80
!
!! TEST28 tests RR_L1.
!
  implicit none

  integer ( kind = 4 ), parameter :: lambda_num = 3

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  real ( kind = 4 ) lambda
  real ( kind = 4 ), dimension ( lambda_num ) :: lambda_test = &
    (/ 0.0E+00, 1.0E+00, 100.0E+00 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l1
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  RR_L1 minimizes the L1 norm of:'
  write ( *, '(a)' ) '    [ A          ] * X - [ B ]'
  write ( *, '(a)' ) '    [ LAMBDA * I ]       [ 0 ] '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m+n) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  do i = 1, lambda_num

    lambda = lambda_test(i)
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Solution for LAMBDA = ', lambda

    call rr_l1 ( a, m, n, b, eps, x, lambda, iflag )

    if ( iflag /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST28 - Warning!'
      write ( *, '(a,i6)' ) '  RR_L1 returned IFLAG = ', iflag
      return
    end if
!
!  Print the solution and residual.
!
    call r4vec_print ( n, x, '  Solution:' )

    r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)
    r(m+1:m+n) = lambda * x(1:n)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  L1 norm of AX-B =     ', &
      r4vec_norm_l1 ( m, r(1:m) )
    write ( *, '(a,g14.6)' ) '  L1 norm of LAMBDA*X = ', &
      r4vec_norm_l1 ( n, r(m+1:m+n) )
    write ( *, '(a,g14.6)' ) '  L1 norm of residual = ', &
      r4vec_norm_l1 ( m+n, r )

  end do

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test29 ( file_name )

!*****************************************************************************80
!
!! TEST29 tests RR_LI.
!
  implicit none

  integer ( kind = 4 ), parameter :: lambda_num = 3

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps1 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps2 = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  real ( kind = 4 ) lambda
  real ( kind = 4 ), dimension ( lambda_num ) :: lambda_test = &
    (/ 0.0E+00, 1.0E+00, 100.0E+00 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_li
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  RR_LI minimizes the L-infinity norm of:'
  write ( *, '(a)' ) '    [ A          ] * X - [ B ]'
  write ( *, '(a)' ) '    [ LAMBDA * I ]       [ 0 ] '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m+n) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  do i = 1, lambda_num

    lambda = lambda_test(i)
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Solution for LAMBDA = ', lambda

    call rr_li ( a, m, n, b, eps1, eps2, x, lambda, iflag )

    if ( iflag /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST29 - Warning!'
      write ( *, '(a,i6)' ) '  RR_LI returned IFLAG = ', iflag
      return
    end if
!
!  Print the solution and residual.
!
    call r4vec_print ( n, x, '  Solution:' )

    r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)
    r(m+1:m+n) = lambda * x(1:n)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  L-infinity norm of AX-B =     ', &
      r4vec_norm_li ( m, r(1:m) )
    write ( *, '(a,g14.6)' ) '  L-infinity norm of LAMBDA*X = ', &
      r4vec_norm_li ( n, r(m+1:m+n) )
    write ( *, '(a,g14.6)' ) '  L-infinity norm of residual = ', &
      r4vec_norm_li ( m+n, r )

  end do

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test30 ( file_name )

!*****************************************************************************80
!
!! TEST30 tests ORTH_L2.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps1 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps2 = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), parameter :: itmax = 15
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  logical start
  real ( kind = 4 ), allocatable, dimension ( : ) :: x
  real ( kind = 4 ), allocatable, dimension ( : ) :: x_start

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  ORTH_L2 minimizes the L2 norm of:'
  write ( *, '(a)' ) '    A*X - X(N+1)*B'
  write ( *, '(a)' ) '  with L2 norm of X equal to 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n+1) )
  allocate ( x_start(n+1) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  start = .false.

  call orth_l2 ( a, m, n, b, eps1, eps2, itmax, x, start, x_start, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST30 - Warning!'
    write ( *, '(a,i6)' ) '  ORTH_L2 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call r4vec_print ( n+1, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - x(n+1) * b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of AX-X(N+1)*B =     ', &
    r4vec_norm_l2 ( m, r )
  write ( *, '(a,g14.6)' ) '  L2 norm of X =               ', &
    r4vec_norm_l2 ( n+1, x )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )
  deallocate ( x_start )

  return
end
subroutine test305 ( file_name )

!*****************************************************************************80
!
!! TEST305 tests ORTH_L1.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps1 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps2 = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), parameter :: itmax = 15
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l1
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST305'
  write ( *, '(a)' ) '  ORTH_L1 minimizes the L1 norm of:'
  write ( *, '(a)' ) '    A*X - X(N+1)*B'
  write ( *, '(a)' ) '  with L2 norm of X equal to 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n+1) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call orth_l1 ( a, m, n, b, eps1, eps2, itmax, x, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST305 - Warning!'
    write ( *, '(a,i6)' ) '  ORTH_L1 returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call r4vec_print ( n+1, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - x(n+1) * b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L1 norm of AX-X(N+1)*B =     ', &
    r4vec_norm_l1 ( m, r )
  write ( *, '(a,g14.6)' ) '  L2 norm of X =               ', &
    r4vec_norm_l2 ( n+1, x )

  x(1:n+1) = (/ 0.730791E+00, -0.138925E+00, 0.668315E+00 /)

  call r4vec_print ( n+1, x, '  Spaeth solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - x(n+1) * b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L1 norm of AX-X(N+1)*B =     ', &
    r4vec_norm_l1 ( m, r )
  write ( *, '(a,g14.6)' ) '  L2 norm of X =               ', &
    r4vec_norm_l2 ( n+1, x )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test306 ( file_name )

!*****************************************************************************80
!
!! TEST306 tests ORTH_LI.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps1 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps2 = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), parameter :: itmax = 15
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_li
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: x
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST306'
  write ( *, '(a)' ) '  ORTH_LI minimizes the LI norm of:'
  write ( *, '(a)' ) '    A*X - X(N+1)*B'
  write ( *, '(a)' ) '  with L2 norm of X equal to 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n+1) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call orth_li ( a, m, n, b, eps1, eps2, itmax, x, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST306 - Warning!'
    write ( *, '(a,i6)' ) '  ORTH_LI returned IFLAG = ', iflag
    return
  end if
!
!  Print the solution and residual.
!
  call r4vec_print ( n+1, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - x(n+1) * b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  LI norm of AX-X(N+1)*B =     ', &
    r4vec_norm_li ( m, r )
  write ( *, '(a,g14.6)' ) '  L2 norm of X =               ', &
    r4vec_norm_l2 ( n+1, x )

  x(1:n+1) = (/ 0.699750E+00, -0.0845856E+00, 0.709362E+00 /)

  call r4vec_print ( n+1, x, '  Spaeth solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - x(n+1) * b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  LI norm of AX-X(N+1)*B =     ', &
    r4vec_norm_li ( m, r )
  write ( *, '(a,g14.6)' ) '  L2 norm of X =               ', &
    r4vec_norm_l2 ( n+1, x )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test307 ( file_name )

!*****************************************************************************80
!
!! TEST307 tests ORTH_LP.
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 8

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps1 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps2 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps3 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps4 = 0.000001E+00
  real ( kind = 4 ), parameter :: eps5 = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), parameter :: itmax1 = 30
  integer ( kind = 4 ), parameter :: itmax2 = 50
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ) p
  real ( kind = 4 ), dimension ( test_num ) :: p_test = (/ &
    1.1E+00, 1.2E+00, 1.4E+00, 1.7E+00, 2.0E+00, 2.5E+00, 4.0E+00, 6.5E+00 /)
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ) r4vec_norm_lp
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST307'
  write ( *, '(a)' ) '  ORTH_LP minimizes the LP norm of:'
  write ( *, '(a)' ) '    A*X - X(N+1)*B'
  write ( *, '(a)' ) '  with L2 norm of X equal to 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( x(n+1) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  do i = 1, test_num

    p = p_test(i)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  P = ', p

    call orth_lp ( a, m, n, b, p, eps1, eps2, eps3, eps4, eps5, &
      itmax1, itmax2, x, iflag )

    if ( iflag /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST307 - Warning!'
      write ( *, '(a,i6)' ) '  ORTH_LP returned IFLAG = ', iflag
    end if
!
!  Print the solution and residual.
!
    call r4vec_print ( n+1, x, '  Solution:' )

    r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - x(n+1) * b(1:m)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  LP norm of AX-X(N+1)*B =     ', &
      r4vec_norm_lp ( m, r, p )
    write ( *, '(a,g14.6)' ) '  L2 norm of X =               ', &
      r4vec_norm_l2 ( n+1, x )

  end do

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )

  return
end
subroutine test31 ( file_name )

!*****************************************************************************80
!
!! TEST31 tests SVDRS.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: s
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  SVDRS uses the singular value decomposition to'
  write ( *, '(a)' ) '  solve a least squares problem.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( s(n) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  call svdrs ( a, m, m, n, b, m, 1, s )

  call r4vec_print ( n, s, '  Singular values S:' )

  do i = 1, m
    if ( i <= n ) then
      if ( s(i) /= 0.0E+00 ) then
        b(i) = b(i) / s(i)
      else
        b(i) = 0.0
      end if
    else
      b(i) = 0.0E+00
    end if
  end do

  x(1:n) = matmul ( a(1:n,1:n), b(1:n) )

  call r4vec_print ( n, x, '  Solution = V * inv(S) * U'' * B:' )
!
!  Print the solution and residual.
!
!  SVDRS destroys the matrix A, so we need to recover the
!  original data in order to compute the residual.
!
  call example_read ( file_name, m, n, a, b, a_title, b_title )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of residual = ', r4vec_norm_l2 ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( s )
  deallocate ( x )

  return
end
subroutine test32 ( file_name )

!*****************************************************************************80
!
!! TEST32 tests SVD.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), parameter :: eps = 0.000001E+00
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) m
  logical matu
  logical matv
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r
  real ( kind = 4 ) r4vec_norm_l2
  real ( kind = 4 ), allocatable, dimension ( : ) :: s
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: u
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: v
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32'
  write ( *, '(a)' ) '  SVD applies the singular value decomposition to'
  write ( *, '(a)' ) '  a matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n

  allocate ( a(m,n) )
  allocate ( a_title(n) )
  allocate ( b(m) )
  allocate ( r(m) )
  allocate ( s(n) )
  allocate ( u(m,n) )
  allocate ( v(n,n) )
  allocate ( x(n) )

  call example_read ( file_name, m, n, a, b, a_title, b_title )
!
!  Solve the linear system.
!
  matu = .true.
  matv = .true.
  call svd ( m, n, a, s, matu, u, matv, v, iflag )

  call r4vec_print ( n, s, '  Singular values S:' )

  x(1:n) = matmul ( transpose ( u(1:m,1:n) ), b(1:m) )

  do i = 1, n
    if ( s(i) /= 0.0E+00 ) then
      x(i) = x(i) / s(i)
    end if
  end do

  x(1:n) = matmul ( v(1:n,1:n), x(1:n) )
!
!  Print the solution and residual.
!
!  SVD destroys the matrix A, so we need to recover the
!  original data in order to compute the residual.
!
  call example_read ( file_name, m, n, a, b, a_title, b_title )

  call r4vec_print ( n, x, '  Solution:' )

  r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of residual = ', r4vec_norm_l2 ( m, r )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( r )
  deallocate ( s )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )

  return
end
subroutine test33 ( file_name )

!*****************************************************************************80
!
!! TEST33 tests CON_L1.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: c
  integer ( kind = 4 ) code
  real ( kind = 4 ), allocatable, dimension ( : ) :: d
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: e
  real ( kind = 4 ), parameter :: eps = 0.0001E+00
  real ( kind = 4 ), allocatable, dimension ( : ) :: f
  character ( len = * ) file_name
  character ( len = 80 ) file_name2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), parameter :: itmax = 20
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: res
  real ( kind = 4 ), allocatable, dimension ( : ) :: r1
  real ( kind = 4 ), allocatable, dimension ( : ) :: r2
  real ( kind = 4 ), allocatable, dimension ( : ) :: r3
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) r4vec_norm_l1
  integer ( kind = 4 ) s
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33'
  write ( *, '(a)' ) '  CON_L1 solves a constrained minimization problem'
  write ( *, '(a)' ) '  in the L1 norm:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Find an N vector X which minimizes the residual:'
  write ( *, '(a)' ) '    || A * X - B ||'
  write ( *, '(a)' ) '  and satisifies the linear equalities:'
  write ( *, '(a)' ) '    C * X = D'
  write ( *, '(a)' ) '  and the linear inequalities:'
  write ( *, '(a)' ) '    E * X >= F.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_multi_size ( file_name, m, n, s, file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n
  write ( *, '(a,i6)' ) '  Number of subsystems =        ', s

  allocate ( a_title(n) )
!
!  Read A and B.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name2 )

  call example_size ( file_name2, m, n )

  allocate ( a(m,n) )
  allocate ( b(m) )
  allocate ( r1(m) )

  call example_read ( file_name2, m, n, a, b, a_title, b_title )
!
!  Read C and D.
!
  call file_name_inc ( file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name2 )

  call example_size ( file_name2, l, n )

  allocate ( c(l,n) )
  allocate ( d(l) )
  allocate ( r2(l) )

  call example_read ( file_name2, l, n, c, d, a_title, b_title )
!
!  Read E and F.
!
  call file_name_inc ( file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name2 )

  call example_size ( file_name2, k, n )

  allocate ( e(k,n) )
  allocate ( f(k) )
  allocate ( r3(k) )

  call example_read ( file_name2, k, n, e, f, a_title, b_title )
!
!  Some more preparations.
!
  allocate ( res(m+l+k) )
  allocate ( x(n) )
  code = 0
  x(1:n) = 0.0E+00
!
!  Solve the linear system.
!
  call con_l1 ( m, l, k, n, a, b, c, d, e, f, code, eps, itmax, x, &
    res, res_norm, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST33 - Warning!'
    write ( *, '(a,i6)' ) '  CON_L1 returned error flag IFLAG = ', iflag
  end if
!
!  Print the solution and residual.
!
  call r4vec_print ( n, x, '  Solution:' )

  r1(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)
  r2(1:l) = matmul ( c(1:l,1:n), x(1:n) ) - d(1:l)
  r3(1:k) = matmul ( e(1:k,1:n), x(1:n) ) - f(1:k)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L1 norms of residuals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CON_L1 claims A*X-B residual norm = ', res_norm
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A*X-B = ', r4vec_norm_l1 ( m, r1 )
  write ( *, '(a,g14.6)' ) '  C*X-D = ', r4vec_norm_l1 ( l, r2 )
  write ( *, '(a,g14.6)' ) '  E*X-F = ', r4vec_norm_l1 ( k, r3 )
  call r4vec_print ( k, r3, '  E*X-F:' )

  if ( file_name == 'x54.txt' ) then

    x(1:n) = &
      (/ 0.0816327E+00, 0.0E+00, 0.0E+00, 0.0782313E+00, -0.0646259E+00 /)

    call r4vec_print ( n, x, '  Spaeth''s solution:' )

    r1(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)
    r2(1:l) = matmul ( c(1:l,1:n), x(1:n) ) - d(1:l)
    r3(1:k) = matmul ( e(1:k,1:n), x(1:n) ) - f(1:k)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  L1 norms of residuals:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  A*X-B = ', r4vec_norm_l1 ( m, r1 )
    write ( *, '(a,g14.6)' ) '  C*X-D = ', r4vec_norm_l1 ( l, r2 )
    write ( *, '(a,g14.6)' ) '  E*X-F = ', r4vec_norm_l1 ( k, r3 )
    call r4vec_print ( k, r3, '  E*X-F:' )

  else if ( file_name == 'x61.txt' ) then

    x(1:n) = &
      (/ 0.0, 1.737705, 0.0, -0.2377049, -0.1885246 /)

    call r4vec_print ( n, x, '  Barrodale and Roberts''s solution:' )

    r1(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)
    r2(1:l) = matmul ( c(1:l,1:n), x(1:n) ) - d(1:l)
    r3(1:k) = matmul ( e(1:k,1:n), x(1:n) ) - f(1:k)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  L1 norms of residuals:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  A*X-B = ', r4vec_norm_l1 ( m, r1 )
    write ( *, '(a,g14.6)' ) '  C*X-D = ', r4vec_norm_l1 ( l, r2 )
    write ( *, '(a,g14.6)' ) '  E*X-F = ', r4vec_norm_l1 ( k, r3 )
    call r4vec_print ( k, r3, '  E*X-F:' )

  end if

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( c )
  deallocate ( d )
  deallocate ( e )
  deallocate ( f )
  deallocate ( r1 )
  deallocate ( r2 )
  deallocate ( r3 )
  deallocate ( res )
  deallocate ( x )

  return
end
subroutine test34 ( file_name )

!*****************************************************************************80
!
!! TEST34 tests CON_L2.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: c
  real ( kind = 4 ), allocatable, dimension ( : ) :: d
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: e
  real ( kind = 4 ), parameter :: eps = 0.0001E+00
  real ( kind = 4 ), allocatable, dimension ( : ) :: f
  character ( len = * ) file_name
  character ( len = 80 ) file_name2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mlk
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rank1
  integer ( kind = 4 ) rank2
  real ( kind = 4 ), allocatable, dimension ( : ) :: res
  real ( kind = 4 ), allocatable, dimension ( : ) :: r1
  real ( kind = 4 ), allocatable, dimension ( : ) :: r2
  real ( kind = 4 ), allocatable, dimension ( : ) :: r3
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) r4vec_norm_l2
  integer ( kind = 4 ) s
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST34'
  write ( *, '(a)' ) '  CON_L2 solves a constrained minimization problem'
  write ( *, '(a)' ) '  in the L2 norm:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Find an N vector X which minimizes the residual:'
  write ( *, '(a)' ) '    || A * X - B ||'
  write ( *, '(a)' ) '  and satisifies the linear equalities:'
  write ( *, '(a)' ) '    C * X = D'
  write ( *, '(a)' ) '  and the linear inequalities:'
  write ( *, '(a)' ) '    E * X >= F.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_multi_size ( file_name, mlk, n, s, file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', mlk
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n
  write ( *, '(a,i6)' ) '  Number of subsystems =        ', s

  allocate ( a_title(n) )
!
!  Read A and B.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name2 )

  call example_size ( file_name2, m, n )
  write ( *, '(a,i6)' ) '  Number of minimizine equations M = ', m

  allocate ( a(m,n) )
  allocate ( b(m) )
  allocate ( r1(m) )

  call example_read ( file_name2, m, n, a, b, a_title, b_title )
!
!  Read C and D.
!
  call file_name_inc ( file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name2 )

  call example_size ( file_name2, l, n )

  write ( *, '(a,i6)' ) '  Number of equality constraints L = ', l

  allocate ( c(l,n) )
  allocate ( d(l) )
  allocate ( r2(l) )

  call example_read ( file_name2, l, n, c, d, a_title, b_title )
!
!  Read E and F.
!
  call file_name_inc ( file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name2 )

  call example_size ( file_name2, k, n )

  write ( *, '(a,i6)' ) '  Number of inequality constraints K = ', k

  allocate ( e(k,n) )
  allocate ( f(k) )
  allocate ( r3(k) )

  call example_read ( file_name2, k, n, e, f, a_title, b_title )
!
!  Some more preparations.
!
  allocate ( x(n) )

  x(1:n) = 0.0E+00
!
!  Solve the linear system.
!
  call con_l2 ( m, l, k, n, a, b, c, d, e, f, eps, x, rank1, rank2, &
    res_norm, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST34 - Warning!'
    write ( *, '(a,i6)' ) '  CON_L2 returned error flag IFLAG = ', iflag
  end if
!
!  Print the solution and residual.
!
  call r4vec_print ( n, x, '  Solution:' )

  r1(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)
  r2(1:l) = matmul ( c(1:l,1:n), x(1:n) ) - d(1:l)
  r3(1:k) = matmul ( e(1:k,1:n), x(1:n) ) - f(1:k)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norms of residuals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CON_L2 claims A*X-B residual norm = ', res_norm
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A*X-B = ', r4vec_norm_l2 ( m, r1 )
  write ( *, '(a,g14.6)' ) '  C*X-D = ', r4vec_norm_l2 ( l, r2 )
  write ( *, '(a,g14.6)' ) '  E*X-F = ', r4vec_norm_l2 ( k, r3 )
  call r4vec_print ( k, r3, '  E*X-F:' )

  x(1:n) = (/ -0.0411582E+00, 0.0411582E+00, 0.252699E-08, 0.101372E+00, &
    -0.0414852E+00 /)

  call r4vec_print ( n, x, '  Spaeth''s solution:' )

  r1(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)
  r2(1:l) = matmul ( c(1:l,1:n), x(1:n) ) - d(1:l)
  r3(1:k) = matmul ( e(1:k,1:n), x(1:n) ) - f(1:k)

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norms of residuals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A*X-B = ', r4vec_norm_l2 ( m, r1 )
  write ( *, '(a,g14.6)' ) '  C*X-D = ', r4vec_norm_l2 ( l, r2 )
  write ( *, '(a,g14.6)' ) '  E*X-F = ', r4vec_norm_l2 ( k, r3 )
  call r4vec_print ( k, r3, '  E*X-F:' )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( c )
  deallocate ( d )
  deallocate ( e )
  deallocate ( f )
  deallocate ( r1 )
  deallocate ( r2 )
  deallocate ( r3 )
  deallocate ( x )

  return
end
subroutine test35 ( file_name )

!*****************************************************************************80
!
!! TEST35 tests CON_LI.
!
  implicit none

  real ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  character ( len = 20 ), allocatable, dimension ( : ) :: a_title
  real ( kind = 4 ), allocatable, dimension ( : ) :: b
  character ( len = 20 ) b_title
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: c
  real ( kind = 4 ), allocatable, dimension ( : ) :: d
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: e
  real ( kind = 4 ), parameter :: eps = 0.0001E+00
  real ( kind = 4 ), allocatable, dimension ( : ) :: f
  character ( len = * ) file_name
  character ( len = 80 ) file_name2
  real ( kind = 4 ), allocatable, dimension ( : ) :: g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 4 ), allocatable, dimension ( : ) :: r1
  real ( kind = 4 ), allocatable, dimension ( : ) :: r2
  real ( kind = 4 ), allocatable, dimension ( : ) :: r3
  real ( kind = 4 ), allocatable, dimension ( : ) :: r4
  real ( kind = 4 ) res_norm
  real ( kind = 4 ) r4vec_norm_li
  integer ( kind = 4 ) s
  real ( kind = 4 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST35'
  write ( *, '(a)' ) '  CON_LI solves a constrained minimization problem'
  write ( *, '(a)' ) '  in the LI norm:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Find an N vector X which minimizes the residual:'
  write ( *, '(a)' ) '    || A * X - B ||'
  write ( *, '(a)' ) '  and satisifies the linear equalities:'
  write ( *, '(a)' ) '    C * X = D'
  write ( *, '(a)' ) '  and the linear inequalities:'
  write ( *, '(a)' ) '    F <= E * X <= G.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name )

  call example_multi_size ( file_name, m, n, s, file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of rows/observations = ', m
  write ( *, '(a,i6)' ) '  Number of columns/variables = ', n
  write ( *, '(a,i6)' ) '  Number of subsystems =        ', s

  allocate ( a_title(n) )
!
!  Read A and B.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name2 )

  call example_size ( file_name2, m, n )

  allocate ( a(m,n) )
  allocate ( b(m) )
  allocate ( r1(m) )

  call example_read ( file_name2, m, n, a, b, a_title, b_title )
!
!  Read C and D.
!
  call file_name_inc ( file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name2 )

  call example_size ( file_name2, l, n )

  allocate ( c(l,n) )
  allocate ( d(l) )
  allocate ( r2(l) )

  call example_read ( file_name2, l, n, c, d, a_title, b_title )
!
!  Read E and F.
!  We don't care about G, so set it large.
!
  call file_name_inc ( file_name2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Open data file ' // trim ( file_name2 )

  call example_size ( file_name2, k, n )

  allocate ( e(k,n) )
  allocate ( f(k) )
  allocate ( g(k) )
  allocate ( r3(k) )
  allocate ( r4(k) )

  call example_read ( file_name2, k, n, e, f, a_title, b_title )
  g(1:k) = 100.0E+00
!
!  Some more preparations.
!
  allocate ( x(n) )

  x(1:n) = 0.0E+00
!
!  Solve the linear system.
!
  call con_li ( m, l, k, n, a, b, c, d, e, f, g, eps, x, res_norm, &
    iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST35 - Warning!'
    write ( *, '(a,i6)' ) '  CON_LI returned error flag IFLAG = ', iflag
  end if
!
!  Print the solution and residual.
!
  call r4vec_print ( n, x, '  Solution:' )

  r1(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)
  r2(1:l) = matmul ( c(1:l,1:n), x(1:n) ) - d(1:l)
  r3(1:k) = matmul ( e(1:k,1:n), x(1:n) ) - f(1:k)
  r4(1:k) = g(1:k) - matmul ( e(1:k,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  LI norms of residuals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CON_LI claims A*X-B residual norm = ', res_norm
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A*X-B = ', r4vec_norm_li ( m, r1 )
  write ( *, '(a,g14.6)' ) '  C*X-D = ', r4vec_norm_li ( l, r2 )
  call r4vec_print ( k, r3, '  E*X-F:' )
  call r4vec_print ( k, r4, '  G-E*X:' )

  x(1:n) = (/ 0.0408165E+00, 0.0E+00, 0.0E+00, 0.0891156E+00, -0.0537416E+00 /)

  call r4vec_print ( n, x, '  Spaeth''s solution:' )

  r1(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)
  r2(1:l) = matmul ( c(1:l,1:n), x(1:n) ) - d(1:l)
  r3(1:k) = matmul ( e(1:k,1:n), x(1:n) ) - f(1:k)
  r3(1:k) = g(1:k) - matmul ( e(1:k,1:n), x(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  LI norms of residuals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A*X-B = ', r4vec_norm_li ( m, r1 )
  write ( *, '(a,g14.6)' ) '  C*X-D = ', r4vec_norm_li ( l, r2 )
  call r4vec_print ( k, r3, '  E*X-F:' )
  call r4vec_print ( k, r4, '  G-E*X:' )

  deallocate ( a )
  deallocate ( a_title )
  deallocate ( b )
  deallocate ( c )
  deallocate ( d )
  deallocate ( e )
  deallocate ( f )
  deallocate ( g )
  deallocate ( r1 )
  deallocate ( r2 )
  deallocate ( r3 )
  deallocate ( r4 )
  deallocate ( x )

  return
end
