program main

!*****************************************************************************80
!
!! MAIN is the main program for SOR_GE_PRB.
!
!  Discussion:
!
!    SOR_GE_PRB tests the TEMPLATES SOR GE solver.
!
  implicit none

  integer, parameter :: n = 5

  real a(n,n)
  real b(n)
  integer i
  integer info
  integer iter
  integer iter_max
  integer iter_total
  integer j
  real omega
  real r(n)
  real resid
  real restol
  real s
  real x(n)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SOR_GE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEMPLATES SOR GE solver.'
!
!  Set the matrix A.
!
  call dif ( n, n, a )
!
!  Set the solution X.
!
  s = 1.0
  do i = 1, n
    x(i) = s * real ( i )
    s = - s
  end do
!
!  Set the right hand side B.
!
  b(1:n) = 0.0

  call matvec_ge ( n, a, 1.0, x, 0.0, b, b )

  call matvec_ge ( n, a, 1.0, x, -1.0, b, r )

  call r4vec_print_some ( n, b, 10, '  Right hand side B:' )

  call r4vec_print_some ( n, x, 10, '  Exact solution X:' )

  x(1:n) = 0.0
!
!  Call SOR_GE
!
  iter = 5
  iter_max = 30
  restol = 0.000001
  omega = 1.2

  iter_total = 0

  do

    if ( iter_max < iter_total ) then
      write ( *, * ) ' '
      write ( *, * ) 'The iteration was halted.'
      write ( *, * ) 'The iteration limit was reached.'
      exit
    end if

    call sor_ge ( n, a, b, x, omega, iter, restol, resid, info )

    iter_total = iter_total + iter

    write ( *, * ) ' '
    write ( *, * ) 'The total number of iterations was ', iter_total
    write ( *, * ) 'Convergence measure is ', resid

    call matvec_ge ( n, a, 1.0, x, -1.0, b, r )

    call r4vec_print_some ( n, x, 10, '  Computed solution X:' )

    call r4vec_print_some ( n, r, 10, '  Residual R:' )

    if ( info == 0 ) then
      exit
      write ( *, * ) ' '
      write ( *, * ) '  The convergence tolerance was reached.'
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'SOR_GE_PRB'
  write ( *, * ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
