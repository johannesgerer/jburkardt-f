program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_NONLIN_PRB.
!
!  Discussion:
!
!    TEST_NONLIN_PRB demonstrates the use of the test functions in TEST_NONLIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_NONLIN_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_NONLIN library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_NONLIN_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 prints the title of each problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Print the title of each problem.'
  write ( *, '(a)' ) ' '
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num
    call p00_title ( problem, title )
    write ( *, '(2x,i2,2x,a)' ) problem, '"' // trim ( title ) // '"'
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 uses a simple Newton method on all the problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) abstol
  logical converged
  real ( kind = 8 ) damp
  real ( kind = 8 ) damptol
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  real ( kind = 8 ) fnrm
  real ( kind = 8 ) fnrm_old
  real ( kind = 8 ) fnrm_zero
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: fprime
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) iknow
  integer ( kind = 4 ) info
  integer ( kind = 4 ) it
  integer ( kind = 4 ) maxit
  integer ( kind = 4 ) n
  integer ( kind = 4 ), allocatable, dimension ( : ) :: pivot
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  real ( kind = 8 ) reltol
  real ( kind = 8 ) r8vec_norm2
  real ( kind = 8 ) snrm
  real ( kind = 8 ), allocatable, dimension ( : ) :: step
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) ::  x
  real ( kind = 8 ), allocatable, dimension ( : ) ::  xnew
  real ( kind = 8 ) xnrm
  real ( kind = 8 ) xold

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Seek roots of the nonlinear functions'
  write ( *, '(a)' ) '  in TEST_NONLIN using Newton''s method.'
  write ( *, '(a)' ) ' '
!
!  Initialization.
!
  abstol = 1.0D-06
  damptol = 1.0D-06
  reltol = 1.0D-06
  maxit = 40
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )
!
!  Solve each problem.
!
  do problem = 1, problem_num
!
!  Get the default problem size.
!  If none, use N = 10.
!
    call p00_n ( problem, n )

    if ( n < 0 ) then
      n = 10
    end if
!
!  Allocate memory.
!
    allocate ( f(n) )
    allocate ( fprime(n,n) )
    allocate ( pivot(n) )
    allocate ( step(n) )
    allocate ( x(n) )
    allocate ( xnew(n) )
!
!  Get the problem title.
!
    call p00_title ( problem, title )

    do ijac = 0, 1

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Problem index ', problem
      write ( *, '(3x,a)' ) trim ( title )
      write ( *, '(a,i6)' ) '  Problem size = ', n
      if ( ijac == 0 ) then
        write ( *, '(a)' ) '  Use analytic jacobian'
      else
        write ( *, '(a)' ) '  Use finite difference jacobian'
      end if
!
!  Get the problem starting point.
!
      call p00_start ( problem, n, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Initial X:'
      write ( *, '(a)' ) ' '
      write ( *, '(5g14.6)' ) x(1:n)
      xnrm = r8vec_norm2 ( n, x )
!
!  Evaluate function at starting point.
!
      call p00_fx ( problem, f, n, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Function value at initial X:'
      write ( *, '(a)' ) ' '
      write ( *, '(5g14.6)' ) f(1:n)
      fnrm = r8vec_norm2 ( n, f )
      fnrm_old = fnrm
      fnrm_zero = fnrm

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Iteration   ||X||         ||F(X)||   Damping   ||dX||'
      write ( *, '(a)' ) ' '
      it = 0
      write ( *, '(7x,i3,g14.6,g14.6)' ) it, xnrm, fnrm

      converged = .false.

      do it = 1, maxit
!
!  Get the Jacobian matrix.
!
        if ( ijac == 0 ) then
          call p00_jac ( problem, n, fprime, x )
        else
          call p00_dif ( problem, n, fprime, x )
        end if
!
!  Factor the jacobian matrix.
!
        call r8ge_fa ( n, fprime, pivot, info )

        if ( info /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TEST02 - Warning:'
          write ( *, '(a)' ) '  The iteration must be halted.'
          write ( *, '(a)' ) '  The jacobian matrix is singular!'
          exit
        end if
!
!  Solve Fprime * DeltaX = F.
!
        step(1:n) = f(1:n)

        call r8ge_sl ( n, fprime, pivot, step, 0 )

        snrm = r8vec_norm2 ( n, step )

        if ( 10.0D+00 * ( xnrm + 10.0D+00 ) < snrm .and. 1 < it ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TEST02'
          write ( *, '(a,g14.6)' ) '  Excessive correction norm = ', snrm
          write ( *, '(a,g14.6)' ) '  Solution norm was           ', xnrm
          write ( *, '(a)' ) '  The iteration is terminated.'
          exit
        end if
!
!  Update XNEW := X - DeltaX.
!
        if ( snrm < 2.0D+00 * ( xnrm + 1.0D+00 ) ) then
          damp = 1.0D+00
        else
          damp = 2.0D+00 * ( xnrm + 1.0D+00 ) / snrm
        end if

        damp = max ( damp, damptol )

10      continue

        xnew(1:n) = x(1:n) - damp * step(1:n)

        xold = xnrm
        xnrm = r8vec_norm2 ( n, xnew )
        call p00_fx ( problem, f, n, xnew )
        fnrm = r8vec_norm2 ( n, f )
!
!  Should we reject this step?
!
        if ( 1.2D+00 * ( fnrm_old + abstol ) < fnrm ) then

          if ( damptol < damp ) then
            damp = 0.50D+00 * damp
            go to 10
          else
            write ( *, '(7x,3x,14x,14x,14x,g14.6)' ) damp * snrm
            if ( damp == 1.0D+00 ) then
              write ( *, '(7x,i3,g14.6,g14.6)' ) it, xnrm, fnrm
            else
              write ( *, '(7x,i3,g14.6,g14.6,g14.6)' ) it, xnrm, fnrm, damp
            end if
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) '  Norm of function increased!'
            write ( *, '(a)' ) '  Iteration halted.'
            exit
          end if

        end if

        if ( 1000.0D+00 * ( xold + 1.0D+00 ) < xnrm ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Norm of X blew up!'
          write ( *, '(a)' ) '  Iteration halted.'

          if ( damptol < damp ) then
            damp = 0.25D+00 * damp
            go to 10
          else
            write ( *, '(7x,3x,14x,14x,14x,g14.6)' ) damp * snrm
            if ( damp == 1.0D+00 ) then
              write ( *, '(7x,i3,g14.6,g14.6)' ) it, xnrm, fnrm
            else
              write ( *, '(7x,i3,g14.6,g14.6,g14.6)' ) it, xnrm, fnrm, damp
            end if
            exit
          end if

        end if

        write ( *, '(7x,3x,14x,14x,14x,g14.6)' ) damp * snrm
        if ( damp == 1.0D+00 ) then
          write ( *, '(7x,i3,g14.6,g14.6)' ) it, xnrm, fnrm
        else
          write ( *, '(7x,i3,g14.6,g14.6,g14.6)' ) it, xnrm, fnrm, damp
        end if
!
!  Accept the point as the next iterate.
!
        x(1:n) = xnew(1:n)

        fnrm_old = fnrm
!
!  Should we stop at this step?
!
        if ( snrm < reltol * ( xnrm + 1.0D+00 ) .and. fnrm < abstol ) then
          converged = .true.
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Convergence criteria satisfied.'
          exit
        end if

      end do

      if ( .not. converged ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The Newton iteration did not converge.'
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Final X:'
      write ( *, '(a)' ) ' '
      write ( *, '(5g14.6)' ) x(1:n)
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Final F(X):'
      write ( *, '(a)' ) ' '
      write ( *, '(5g14.6)' ) f(1:n)
!
!  Check against the exact solution, if it is known.
!
      call p00_sol ( problem, iknow, n, x )

      if ( 0 < iknow ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The exact solution X:'
        write ( *, '(a)' ) ' '
        write ( *, '(5g14.6)' ) x(1:n)
        write ( *, '(a)' ) ' '
        call p00_fx ( problem, f, n, x )
        write ( *, '(a)' ) '  F(X):'
        write ( *, '(a)' ) ' '
        write ( *, '(5g14.6)' ) f(1:n)
      else if ( iknow == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The exact solution is not known.'
      else if ( iknow == -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  An exact solution is not given.'
      end if

    end do
!
!  Deallocate memory.
!
    deallocate ( f )
    deallocate ( fprime )
    deallocate ( pivot )
    deallocate ( step )
    deallocate ( x )
    deallocate ( xnew )

  end do

  return
end
