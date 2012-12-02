program main

!*******************************************************************************
!
!! MAIN is the main program for SGE_MOD_PRB.
!
!  Discussion:
!
!    SGE_MOD_PRB tests the SGE_MODULE example.
!
  use sge_module

  implicit none

  integer, parameter :: n = 5

  real, dimension ( n, n ) :: a
  real, dimension ( n ) :: b
  real det
  integer i
  real sge_det
  real, dimension ( n ) :: x

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGE_MOD_PRB'
  write ( *, '(a)' ) '  A calling program for a demonstration of the'
  write ( *, '(a)' ) '  use of modules.'
!
!  Assign values to A.
!
  call random_number ( harvest = a(1:n,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix A:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,5f12.4)' ) a(i,1:n)
  end do
!
!  Set up the desired solution vector.
!
  do i = 1, n
    x(i) = real ( i )
  end do
!
!  Compute the right hand side B = A*X.
!
  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

  call rvec_print ( n, b, '  The right hand side vector:' )
!
!  Now "register" the matrix A.
!
  call sge_create ( n, a )
!
!  Factor it.
!
  call sge_fa ( )
!
!  Compute its determinant.
!
  det = sge_det ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The matrix determinant is ', det
!
!  Solve the linear system.
!
  call sge_sl ( b, x )

  call rvec_print ( n, x, '  The solution vector:' )
!
!  Discard it.
!
  call sge_delete ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGE_MOD_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
