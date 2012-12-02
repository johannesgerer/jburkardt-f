program main

!*****************************************************************************80
!
!! MAIN is the main program for TOEPLITZ_PRB.
!
!  Discussion:
!
!    TOEPLITZ_PRB tests the TOEPLITZ routines.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOEPLITZ_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOEPLITZ library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOEPLITZ_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CCI_SL.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  complex a(n)
  complex x(n)
  complex x2(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  CCI_SL solves a complex circulant system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call cci_random ( n, a )

  call cci_print ( n, a, '  The circulant matrix:' )
!
!  Set the desired solution.
!
  call c4vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call cci_mxv ( n, a, x, x2 )
!
!  Solve the linear system.
!
  call cci_sl ( n, a, x2 )

  call c4vec_print ( n, x2, '  Solution:' )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests CTO_SL.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  complex a(2*n-1)
  complex b(n)
  integer ( kind = 4 ) job
  complex x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  CTO_SL solves a complex Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call cto_random ( n, a )

  call cto_print ( n, a, '  The Toeplitz matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call c4vec_indicator ( n, x )

    write ( *, '(a)' ) ' '
    if ( job == 0 ) then
      write ( *, '(a)' ) '  Desired solution:'
    else
      write ( *, '(a)' ) '  Desired solution to transposed system:'
    end if

    write ( *, '(a)' ) ' '
    call c4vec_print_some ( n, x, 10 )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call cto_mxv ( n, a, x, b )
    else
      call cto_vxm ( n, a, x, b )
    end if

    write ( *, '(a)' ) ' '
    if ( job == 0 ) then
      write ( *, '(a)' ) '  Right Hand Side:'
    else
      write ( *, '(a)' ) '  Right Hand Side of transposed system:'
    end if

    write ( *, '(a)' ) ' '
    call c4vec_print_some ( n, b, 10 )
!
!  Solve the linear system.
!
    call cto_sl ( n, a, b, x, job )

    write ( *, '(a)' ) ' '
    if ( job == 0 ) then
      write ( *, '(a)' ) '  Solution:'
    else
      write ( *, '(a)' ) '  Solution to transposed system:'
    end if

    write ( *, '(a)' ) ' '
    call c4vec_print_some ( n, x, 10 )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests STO_SL.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 4 ) a(2*n-1)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) job
  real ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  STO_SL solves a real Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call sto_random ( n, a )

  call sto_print ( n, a, '  The Toeplitz matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call r4vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call sto_mxv ( n, a, x, b )
    else
      call sto_vxm ( n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call sto_sl ( n, a, b, x, job )

    if ( job == 0 ) then
      call r4vec_print_some ( n, x, 10, '  Solution:' )
    else
      call r4vec_print_some ( n, x, 10, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests SBTO_MXV, SBTO_VXM.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 3
  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 4 ), dimension ( m, m, l ) ::  a1 = reshape ( (/ &
    1.0E+00, 5.0E+00, 2.0E+00, 5.0E+00, &
    3.0E+00, 6.0E+00, 4.0E+00, 6.0E+00, &
    5.0E+00, 7.0E+00, 6.0E+00, 7.0E+00 /), (/ m, m, l /) )

  real ( kind = 4 ), dimension ( m, m, l-1 ) :: a2 = reshape ( (/ &
    7.0E+00, 8.0E+00, 8.0E+00, 8.0E+00, &
    9.0E+00, 9.0E+00, 0.0E+00, 9.0E+00 /), (/ m, m, l-1 /) )
  real ( kind = 4 ) b(m,l)
  real ( kind = 4 ) x(m,l)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For a real block Toeplitz matrix,'
  write ( *, '(a)' ) '  SBTO_MXV computes A * x.'
  write ( *, '(a)' ) '  SBTO_VXM computes x * A.'

  call sbto_print ( m, l, a1, a2, '  The block Toeplitz matrix:' )

  call r4vec_indicator ( m*l, x )

  call r4vec_print ( m*l, x, '  The vector x:' )

  call sbto_mxv ( m, l, a1, a2, x, b )

  call r4vec_print ( m*l, b, '  The product A*x:' )

  call sbto_vxm ( m, l, a1, a2, x, b )

  call r4vec_print ( m*l, b, '  The product x*A:' )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests SBTO_SL.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: l = 3
  integer ( kind = 4 ), parameter :: n = m * l

  real ( kind = 4 ), dimension ( m, m, l ) ::  a1 = reshape ( (/ &
    9.0E+00, 2.0E+00, 1.0E+00, 8.0E+00, &
    3.0E+00, 6.0E+00, 4.0E+00, 6.0E+00, &
    5.0E+00, 7.0E+00, 6.0E+00, 7.0E+00 /), (/ m, m, l /) )

  real ( kind = 4 ), dimension ( m, m, l-1 ) :: a2 = reshape ( (/ &
    7.0E+00, 8.0E+00, 8.0E+00, 8.0E+00, &
    9.0E+00, 9.0E+00, 0.0E+00, 9.0E+00 /), (/ m, m, l-1 /) )
  real ( kind = 4 ) b(n)
  real ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  SBTO_SL solves a block Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call sbto_print ( m, l, a1, a2, '  The block Toeplitz matrix:' )
!
!  Set the desired solution.
!
  call r4vec_indicator ( n, x )
!
!  Compute the right hand side.
!
  call sbto_mxv ( m, l, a1, a2, x, b )

  call r4vec_print ( n, b, '  Right hand side:' )

  call sbto_sl ( m, l, a1, a2, b, x )

  call r4vec_print ( n, x, '  Computed solution:' )

  return
end
