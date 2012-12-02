program main

!*****************************************************************************80
!
!! MAIN is the main program for MGMRES_PRB.
!
!  Discussion:
!
!    MGMRES_PRB runs the quick checks for the MGMRES code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 August 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MGMRES_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the MGMRES library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MGMRES_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests MGMRES_ST on the simple -1,2-1 matrix.
!
!  Discussion:
!
!    This is a very weak test, since the matrix has such a simple
!    structure, is diagonally dominant (though not strictly), 
!    and is symmetric.
!
!    To make the matrix bigger, simply increase the value of N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ), parameter :: nz_num = 3 * n - 2

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nz_num)
  integer ( kind = 4 ) itr_max
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mr
  real ( kind = 8 ) rhs(n)
  integer ( kind = 4 ) test
  real ( kind = 8 ) tol_abs
  real ( kind = 8 ) tol_rel
  real ( kind = 8 ) x_error
  real ( kind = 8 ) x_estimate(n)
  real ( kind = 8 ) x_exact(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test MGMRES_ST on the simple -1,2-1 matrix.'
!
!  Set the matrix.
!
  k = 0

  do i = 1, n

    if ( 1 < i ) then
      k = k + 1
      ia(k) = i
      ja(k) = i - 1
      a(k) = -1.0D+00
    end if

    k = k + 1
    ia(k) = i
    ja(k) = i
    a(k) = 2.0D+00

    if ( i < n ) then
      k = k + 1
      ia(k) = i
      ja(k) = i + 1
      a(k) = -1.0D+00
    end if

  end do
!
!  Set the right hand side:
!
  rhs(1:n-1) = 0.0D+00
  rhs(n) = real ( n + 1, kind = 8 )
!
!  Set the exact solution.
!
  do i = 1, n
    x_exact(i) = real ( i, kind = 8 )
  end do

  do test = 1, 3
!
!  Set the initial solution estimate.
!
    x_estimate(1:n) = 0.0D+00

    x_error = sqrt ( sum ( ( x_exact(1:n) - x_estimate(1:n) )**2 ) )

    if ( test == 1 ) then
      itr_max = 1
      mr = 20
    else if ( test == 2 ) then
      itr_max = 2
      mr = 10
    else if ( test == 3 ) then
      itr_max = 5
      mr = 4
    end if

    tol_abs = 1.0D-08
    tol_rel = 1.0D-08

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Test ', test
    write ( *, '(a,i8)' ) '  Matrix order N = ', n
    write ( *, '(a,i8)' ) '  Inner iteration limit = ', mr
    write ( *, '(a,i8)' ) '  Outer iteration limit = ', itr_max
    write ( *, '(a,g14.6)' ) '  Initial X_ERROR = ', x_error

    call mgmres_st ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr, &
      tol_abs, tol_rel )

    x_error = sqrt ( sum ( ( x_exact(1:n) - x_estimate(1:n) )**2 ) )

    write ( *, '(a,g14.6)' ) '  Final X_ERROR = ', x_error

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests MGMRES_ST on a 9 by 9 matrix.
!
!  Discussion:
!
!    A = 
!      2  0  0 -1  0  0  0  0  0
!      0  2 -1  0  0  0  0  0  0
!      0 -1  2  0  0  0  0  0  0
!     -1  0  0  2 -1  0  0  0  0
!      0  0  0 -1  2 -1  0  0  0
!      0  0  0  0 -1  2 -1  0  0
!      0  0  0  0  0 -1  2 -1  0
!      0  0  0  0  0  0 -1  2 -1
!      0  0  0  0  0  0  0 -1  2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9
  integer ( kind = 4 ), parameter :: nz_num = 23

  real ( kind = 8 ), dimension(nz_num) :: a = (/ &
     2.0D+00, -1.0D+00, &
     2.0D+00, -1.0D+00, &
    -1.0D+00,  2.0D+00, &
    -1.0D+00,  2.0D+00, -1.0D+00, &
    -1.0D+00,  2.0D+00, -1.0D+00, &
    -1.0D+00,  2.0D+00, -1.0D+00, &
    -1.0D+00,  2.0D+00, -1.0D+00, &
    -1.0D+00,  2.0D+00, -1.0D+00, &
    -1.0D+00,  2.0D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension(nz_num) :: ia = (/ &
    1, 1, &
    2, 2, &
    3, 3, &
    4, 4, 4, &
    5, 5, 5, &
    6, 6, 6, &
    7, 7, 7, &
    8, 8, 8, &
    9, 9 /)
  integer ( kind = 4 ) itr_max
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension(nz_num) :: ja = (/ &
    1, 4, &
    2, 3, &
    2, 3, &
    1, 4, 5, &
    4, 5, 6, &
    5, 6, 7, &
    6, 7, 8, &
    7, 8, 9, &
    8, 9 /)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mr
  real ( kind = 8 ), dimension(n) :: rhs = (/ &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00 /)
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  real ( kind = 8 ) tol_abs
  real ( kind = 8 ) tol_rel
  real ( kind = 8 ) x_error
  real ( kind = 8 ) x_estimate(n)
  real ( kind = 8 ), dimension(n) :: x_exact = (/ &
    3.5D+00, &
    1.0D+00, &
    1.0D+00, &
    6.0D+00, &
    7.5D+00, &
    8.0D+00, &
    7.5D+00, &
    6.0D+00, &
    3.5D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Test MGMRES_ST on a matrix that is not quite '
  write ( *, '(a,i8)' ) '  the -1,2,-1 matrix, of order N = ', n

  do test = 1, 2

    if ( test == 1 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First try, use zero initial vector:'
      write ( *, '(a)' ) ' '

      x_estimate(1:n) = 0.0D+00

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Second try, use random initial vector:'
      write ( *, '(a)' ) ' '

      call r8vec_uniform_01 ( n, seed, x_estimate )

    end if

    x_error = sqrt ( sum ( ( x_exact(1:n) - x_estimate(1:n) )**2 ) )

    write ( *, '(a,g14.6)' ) '  Before solving, X_ERROR = ', x_error

    itr_max = 20
    mr = n - 1
    tol_abs = 1.0D-08
    tol_rel = 1.0D-08

    call mgmres_st ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr, &
      tol_abs, tol_rel )

    x_error = sqrt ( sum ( ( x_exact(1:n) - x_estimate(1:n) )**2 ) )

    write ( *, '(a,g14.6)' ) '  After solving, X_ERROR = ', x_error

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Final solution estimate:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i8,2x,g14.6)' ) i, x_estimate(i)
    end do

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests PMGMRES_ILU_CR on the simple -1,2-1 matrix.
!
!  Discussion:
!
!    This is a very weak test, since the matrix has such a simple
!    structure, is diagonally dominant (though not strictly), 
!    and is symmetric.
!
!    To make the matrix bigger, simply increase the value of N.
!
!    Note that PGMRES_ILU_CR expects the matrix to be stored using the
!    sparse compressed row format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20
  integer ( kind = 4 ), parameter :: nz_num = ( 3 * n - 2 )

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) itr_max
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mr
  real ( kind = 8 ) rhs(n)
  integer ( kind = 4 ) test
  real ( kind = 8 ) tol_abs
  real ( kind = 8 ) tol_rel
  real ( kind = 8 ) x_error
  real ( kind = 8 ) x_estimate(n)
  real ( kind = 8 ) x_exact(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Test PMGMRES_ILU_CR on the simple -1,2-1 matrix.'
!
!  Set the matrix.
!  Note that we use 1-based index values in IA and JA.
!
  k = 1
  ia(1) = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a,i4)' ) '  ia(', 1, ') = ', ia(1)

  do i = 1, n

    ia(i+1) = ia(i)

    if ( 1 < i ) then
      ia(i+1) = ia(i+1) + 1
      ja(k) = i - 1
      a(k) = -1.0D+00
      k = k + 1
    end if

    ia(i+1) = ia(i+1) + 1
    ja(k) = i
    a(k) = 2.0D+00
    k = k + 1

    if ( i < n ) then
      ia(i+1) = ia(i+1) + 1
      ja(k) = i + 1
      a(k) = -1.0D+00
      k = k + 1
    end if
    write ( *, '(a,i4,a,i4)' ) '  ia(', i + 1, ') = ', ia(i+1)
  end do
!
!  Set the right hand side:
!
  rhs(1:n-1) = 0.0D+00
  rhs(n) = real ( n + 1, kind = 8 )
!
!  Set the exact solution.
!
  do i = 1, n
    x_exact(i) = real ( i, kind = 8 )
  end do

  do test = 1, 3
!
!  Set the initial solution estimate.
!
    x_estimate(1:n) = 0.0D+00
    x_error = 0.0D+00
    do i = 1, n
      x_error = x_error + ( x_exact(i) - x_estimate(i) ) ** 2
    end do
    x_error = sqrt ( x_error )

    if ( test == 1 ) then
      itr_max = 1
      mr = 20
    else if ( test == 2 ) then
      itr_max = 2
      mr = 10
    else if ( test == 3 ) then
      itr_max = 5
      mr = 4
    end if

    tol_abs = 1.0D-08
    tol_rel = 1.0D-08

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Test ', test
    write ( *, '(a,i4)' ) '  Matrix order N = ', n
    write ( *, '(a,i4)' ) '  Inner iteration limit = ', mr
    write ( *, '(a,i4)' ) '  Outer iteration limit = ', itr_max
    write ( *, '(a,g14.6)' ) '  Initial X_ERROR = ', x_error

    call pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, &
      mr, tol_abs, tol_rel )

    x_error = 0.0D+00
    do i = 1, n
      x_error = x_error + ( x_exact(i) - x_estimate(i) ) ** 2
    end do
    x_error = sqrt ( x_error )

    write ( *, '(a,g14.6)' ) '  Final X_ERROR = ', x_error

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests PMGMRES_ILU_CR on a simple 5 by 5 matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 August 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 9

  real ( kind = 8 ), dimension ( nz_num ) :: a = (/ &
     1.0, 2.0, 1.0, &
     2.0, &
     3.0, 3.0, &
     4.0, &
     1.0, 5.0 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( n + 1 ) :: ia = (/ 1, 4, 5, 7, 8, 10 /)
  integer ( kind = 4 ) itr_max
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension ( nz_num ) :: ja = (/ &
    1, 4, 5, &
    2, &
    1, 3, &
    4, &
    2, 5 /)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mr
  real ( kind = 8 ), dimension ( n ) :: rhs = (/ &
    14.0, 4.0, 12.0, 16.0, 27.0 /)
  integer ( kind = 4 ) test
  real ( kind = 8 ) tol_abs
  real ( kind = 8 ) tol_rel
  real ( kind = 8 ) x_error
  real ( kind = 8 ) x_estimate(n)
  real ( kind = 8 ), dimension ( n ) :: x_exact = (/ 1.0, 2.0, 3.0, 4.0, 5.0 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Test PMGMRES_ILU_CR on a simple 5 x 5 matrix.'

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(a,i2,a,i2)' ) '  ia(', i, ') = ', ia(i)
  end do
 
  do test = 1, 3
!
!  Set the initial solution estimate.
!
    x_estimate(1:n) = 0.0D+00
    x_error = 0.0D+00
    do i = 1, n
      x_error = x_error + ( x_exact(i) - x_estimate(i) ) ** 2
    end do
    x_error = sqrt ( x_error )

    if ( test == 1 ) then
      itr_max = 1
      mr = 20
    else if ( test == 2 ) then
      itr_max = 2
      mr = 10
    else if ( test == 3 ) then
      itr_max = 5
      mr = 4
    end if

    tol_abs = 1.0D-08
    tol_rel = 1.0D-08

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Test ', test
    write ( *, '(a,i4)' ) '  Matrix order N = ', n
    write ( *, '(a,i4)' ) '  Inner iteration limit = ', mr
    write ( *, '(a,i4)' ) '  Outer iteration limit = ', itr_max
    write ( *, '(a,g14.6)' ) '  Initial X_ERROR = ', x_error

    call pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, &
      mr, tol_abs, tol_rel )

    x_error = 0.0D+00
    do i = 1, n
      x_error = x_error + ( x_exact(i) - x_estimate(i) ) ** 2
    end do
    x_error = sqrt ( x_error )

    write ( *, '(a,g14.6)' ) '  Final X_ERROR = ', x_error

  end do

  return
end
