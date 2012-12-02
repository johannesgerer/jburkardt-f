program main

!*****************************************************************************80
!
!! MAIN is the main program for JACOBI_GE_PRB.
!
!  Discussion:
!
!    JACOBI_GE_PRB tests the TEMPLATES Jacobi GE solver.
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
  real r(n)
  real resid
  real restol
  real s
  real x(n)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'JACOBI_GE_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEMPLATES JACOBI_GE solver.'
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

  call r4vec_print_some ( n, x, 10, '  Exact solution X:' )
!
!  Set the right hand side B, by calling the multiply routine!
!
  b(1:n) = 0.0

  call matvec_ge ( n, a, 1.0, x, 0.0, b, b )

  call r4vec_print_some ( n, b, 10, '  Right hand side B:' )
!
!  Wipe out the solution.
!
  x(1:n) = 0.0
  call matvec_ge ( n, a, 1.0, x, -1.0, b, r )

  call r4vec_print_some ( n, r, 10, '  Residual R = A*X-B:' )
!
!  Call JACOBI to solve the system..
!
  iter = 5
  iter_max = 50
  iter_total = 0

  restol = 0.000001

  do 

    if ( iter_max < iter_total ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'JACOBI_GE_PRB:'
      write ( *, '(a)' ) '  The iteration limit was reached.'
      exit
    end if

    call jacobi_ge ( n, a, b, x, iter, restol, resid, info )

    iter_total = iter_total + iter

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_GE_PRB:'
    write ( *, '(a,i6)' ) '  Number of iterations:  ', iter_total
    write ( *, '(a,g14.6)' ) '  Convergence measure is ', resid

    call matvec_ge ( n, a, 1.0, x, -1.0, b, r )

    call r4vec_print_some ( n, x, 10, '  Computed solution X:' )

    call r4vec_print_some ( n, r, 10, '  Residual R:' )

    if ( info == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'JACOBI_GE_PRB:'
      write ( *, '(a)' ) '  The convergence tolerance was reached.'
      exit
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'JACOBI_GE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
