subroutine p00_dif ( problem, n, fjac, x )

!*****************************************************************************80
!
!! P00_DIF approximates the jacobian via finite differences.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Output, real ( kind = 8 ) FJAC(N,N), the approximante N by N
!    jacobian matrix.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be approximated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dxj
  real ( kind = 8 ), parameter :: eps = 0.0001D+00
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) fplus(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xsave

  call p00_fx ( problem, f, n, x )

  do j = 1, n

    if ( 0.0D+00 <= x(j) ) then
      dxj = eps * ( x(j) + 1.0D+00 )
    else
      dxj = eps * ( x(j) - 1.0D+00 )
    end if

    xsave = x(j)
    x(j) = xsave + dxj

    call p00_fx ( problem, fplus, n, x )

    fjac(1:n,j) = ( fplus(1:n) - f(1:n) ) / dxj

    x(j) = xsave

  end do

  return
end
subroutine p00_fx ( problem, f, n, x )

!*****************************************************************************80
!
!! P00_FX evaluates the function for any problem.
!
!  Discussion:
!
!    Most of the problems were originally part of ACM algorithm 566,
!    and were used to test the MINPACK package.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Dennis, David Gay, Phuong Vu,
!    A new nonlinear equations test problem,
!    Technical Report 83-16,
!    Mathematical Sciences Department,
!    Rice University, 1983.
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    Noel deVilliers, David Glasser,
!    A continuation method for nonlinear regression,
!    SIAM Journal on Numerical Analysis,
!    Volume 18, Number 6, December 1981, pages 1139-1154.
!
!    Chris Fraley,
!    Solution of nonlinear least-squares problems,
!    Technical Report STAN-CS-1165,
!    Computer Science Department,
!    Stanford University, 1987.
!
!    Chris Fraley,
!    Software performance on nonlinear least-squares problems,
!    Technical Report SOL 88-17,
!    Systems Optimization Laboratory,
!    Department of Operations Research,
!    Stanford University, 1988.
!
!    JJ McKeown,
!    Specialized versus general-purpose algorithms for functions
!    that are sums of squared terms,
!    Mathematical Programming,
!    Volume 9, 1975, pages 57-68.
!
!    JJ McKeown,
!    On algorithms for sums of squares problems,
!    in Towards Global Optimisation,
!    edited by Laurence Dixon, Gabor Szego,
!    North-Holland, 1975, pages 229-257,
!    ISBN: 0444109552,
!    LC: QA402.5.T7.
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    Algorithm 566:
!    Testing unconstrained optimization software,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 1, March 1981, pages 17-41.
!
!    Douglas Salane,
!    A continuation approach for solving large residual nonlinear
!    least squares problems,
!    SIAM Journal of Scientific and Statistical Computing,
!    Volume 8, Number 4, July 1987, pages 655-671.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the number of the problem.
!
!    Output, real ( kind = 8 ) F(N), the value of the function at X.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which F is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)

  if ( problem == 1 ) then
    call p01_fx ( n, x, f )
  else if ( problem == 2 ) then
    call p02_fx ( n, x, f )
  else if ( problem  ==  3 ) then
    call p03_fx ( n, x, f )
  else if ( problem == 4 ) then
    call p04_fx ( n, x, f )
  else if ( problem == 5 ) then
    call p05_fx ( n, x, f )
  else if ( problem == 6 ) then
    call p06_fx ( n, x, f )
  else if ( problem == 7 ) then
    call p07_fx ( n, x, f )
  else if ( problem == 8 ) then
    call p08_fx ( n, x, f )
  else if ( problem == 9 ) then
    call p09_fx ( n, x, f )
  else if ( problem == 10 ) then
    call p10_fx ( n, x, f )
  else if ( problem == 11 ) then
    call p11_fx ( n, x, f )
  else if ( problem == 12 ) then
    call p12_fx ( n, x, f )
  else if ( problem == 13 ) then
    call p13_fx ( n, x, f )
  else if ( problem == 14 ) then
    call p14_fx ( n, x, f )
  else if ( problem == 15 ) then
    call p15_fx ( n, x, f )
  else if ( problem == 16 ) then
    call p16_fx ( n, x, f )
  else if ( problem == 17 ) then
    call p17_fx ( n, x, f )
  else if ( problem == 18 ) then
    call p18_fx ( n, x, f )
  else if ( problem == 19 ) then
    call p19_fx ( n, x, f )
  else if ( problem == 20 ) then
    call p20_fx ( n, x, f )
  else if ( problem == 21 ) then
    call p21_fx ( n, x, f )
  else if ( problem == 22 ) then
    call p22_fx ( n, x, f )
  else if ( problem == 23 ) then
    call p23_fx ( n, x, f )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FX - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_jac ( problem, n, fjac, x )

!*****************************************************************************80
!
!! P00_JAC evaluates the jacobian for any problem.
!
!  Discussion:
!
!    P00_JAC evaluates the matrix FJAC(I,J)  =  D F(I) / D X(J)
!    given the problem chosen, the number of equations, and the value
!    of the point X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)

  if ( problem == 1 ) then
    call p01_jac ( fjac, n, x )
  else if ( problem == 2 ) then
    call p02_jac ( fjac, n, x )
  else if ( problem  ==  3 ) then
    call p03_jac ( fjac, n, x )
  else if ( problem == 4 ) then
    call p04_jac ( fjac, n, x )
  else if ( problem == 5 ) then
    call p05_jac ( fjac, n, x )
  else if ( problem == 6 ) then
    call p06_jac ( fjac, n, x )
  else if ( problem == 7 ) then
    call p07_jac ( fjac, n, x )
  else if ( problem == 8 ) then
    call p08_jac ( fjac, n, x )
  else if ( problem == 9 ) then
    call p09_jac ( fjac, n, x )
  else if ( problem == 10 ) then
    call p10_jac ( fjac, n, x )
  else if ( problem == 11 ) then
    call p11_jac ( fjac, n, x )
  else if ( problem == 12 ) then
    call p12_jac ( fjac, n, x )
  else if ( problem == 13 ) then
    call p13_jac ( fjac, n, x )
  else if ( problem == 14 ) then
    call p14_jac ( fjac, n, x )
  else if ( problem == 15 ) then
    call p15_jac ( fjac, n, x )
  else if ( problem == 16 ) then
    call p16_jac ( fjac, n, x )
  else if ( problem == 17 ) then
    call p17_jac ( fjac, n, x )
  else if ( problem == 18 ) then
    call p18_jac ( fjac, n, x )
  else if ( problem == 19 ) then
    call p19_jac ( fjac, n, x )
  else if ( problem == 20 ) then
    call p20_jac ( fjac, n, x )
  else if ( problem == 21 ) then
    call p21_jac ( fjac, n, x )
  else if ( problem == 22 ) then
    call p22_jac ( fjac, n, x )
  else if ( problem == 23 ) then
    call p23_jac ( fjac, n, x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_JAC - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_n ( problem, n )

!*****************************************************************************80
!
!! P00_N returns the number of equations for a problem.
!
!  Discussion:
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Output, integer ( kind = 4 ) N, the number of equations.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) problem

  if ( problem == 1 ) then
    call p01_n ( n )
  else if ( problem == 2 ) then
    call p02_n ( n )
  else if ( problem == 3 ) then
    call p03_n ( n )
  else if ( problem == 4 ) then
    call p04_n ( n )
  else if ( problem == 5 ) then
    call p05_n ( n )
  else if ( problem == 6 ) then
    call p06_n ( n )
  else if ( problem == 7 ) then
    call p07_n ( n )
  else if ( problem == 8 ) then
    call p08_n ( n )
  else if ( problem == 9 ) then
    call p09_n ( n )
  else if ( problem == 10 ) then
    call p10_n ( n )
  else if ( problem == 11 ) then
    call p11_n ( n )
  else if ( problem == 12 ) then
    call p12_n ( n )
  else if ( problem == 13 ) then
    call p13_n ( n )
  else if ( problem == 14 ) then
    call p14_n ( n )
  else if ( problem == 15 ) then
    call p15_n ( n )
  else if ( problem == 16 ) then
    call p16_n ( n )
  else if ( problem == 17 ) then
    call p17_n ( n )
  else if ( problem == 18 ) then
    call p18_n ( n )
  else if ( problem == 19 ) then
    call p19_n ( n )
  else if ( problem == 20 ) then
    call p20_n ( n )
  else if ( problem == 21 ) then
    call p21_n ( n )
  else if ( problem == 22 ) then
    call p22_n ( n )
  else if ( problem == 23 ) then
    call p23_n ( n )
  else
    n = 0
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_N - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_problem_num ( problem_num )

!*****************************************************************************80
!
!! P00_PROBLEM_NUM returns the number of problems available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) PROBLEM_NUM, the number of problems available.
!
  implicit none

  integer ( kind = 4 ) problem_num

  problem_num = 23

  return
end
subroutine p00_sol ( problem, iknow, n, x )

!*****************************************************************************80
!
!! P00_SOL returns the solution of any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!   if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)
!
!  Call the appropriate function-specific routine.
!
  if ( problem == 1 ) then
    call p01_sol ( iknow, n, x )
  else if ( problem == 2 ) then
    call p02_sol ( iknow, n, x )
  else if ( problem == 3 ) then
    call p03_sol ( iknow, n, x )
  else if ( problem == 4 ) then
    call p04_sol ( iknow, n, x )
  else if ( problem == 5 ) then
    call p05_sol ( iknow, n, x )
  else if ( problem == 6 ) then
    call p06_sol ( iknow, n, x )
  else if ( problem == 7 ) then
    call p07_sol ( iknow, n, x )
  else if ( problem == 8 ) then
    call p08_sol ( iknow, n, x )
  else if ( problem == 9 ) then
    call p09_sol ( iknow, n, x )
  else if ( problem == 10 ) then
    call p10_sol ( iknow, n, x )
  else if ( problem == 11 ) then
    call p11_sol ( iknow, n, x )
  else if ( problem == 12 ) then
    call p12_sol ( iknow, n, x )
  else if ( problem == 13 ) then
    call p13_sol ( iknow, n, x )
  else if ( problem == 14 ) then
    call p14_sol ( iknow, n, x )
  else if ( problem == 15 ) then
    call p15_sol ( iknow, n, x )
  else if ( problem == 16 ) then
    call p16_sol ( iknow, n, x )
  else if ( problem == 17 ) then
    call p17_sol ( iknow, n, x )
  else if ( problem == 18 ) then
    call p18_sol ( iknow, n, x )
  else if ( problem == 19 ) then
    call p19_sol ( iknow, n, x )
  else if ( problem == 20 ) then
    call p20_sol ( iknow, n, x )
  else if ( problem == 21 ) then
    call p21_sol ( iknow, n, x )
  else if ( problem == 22 ) then
    call p22_sol ( iknow, n, x )
  else if ( problem == 23 ) then
    call p23_sol ( iknow, n, x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SOL - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_start ( problem, n, x )

!*****************************************************************************80
!
!! P00_START specifies a standard approximate solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, is the problem number.
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(n)

  if ( problem == 1 ) then
    call p01_start ( n, x )
  else if ( problem == 2 ) then
    call p02_start ( n, x )
  else if ( problem == 3 ) then
    call p03_start ( n, x )
  else if ( problem == 4 ) then
    call p04_start ( n, x )
  else if ( problem == 5 ) then
    call p05_start ( n, x )
  else if ( problem == 6 ) then
    call p06_start ( n, x )
  else if ( problem == 7 ) then
    call p07_start ( n, x )
  else if ( problem == 8 ) then
    call p08_start ( n, x )
  else if ( problem == 9 ) then
    call p09_start ( n, x )
   else if ( problem == 10 ) then
    call p10_start ( n, x )
  else if ( problem == 11 ) then
    call p11_start ( n, x )
  else if ( problem == 12 ) then
    call p12_start ( n, x )
  else if ( problem == 13 ) then
    call p13_start ( n, x )
  else if ( problem == 14 ) then
    call p14_start ( n, x )
  else if ( problem == 15 ) then
    call p15_start ( n, x )
  else if ( problem == 16 ) then
    call p16_start ( n, x )
  else if ( problem == 17 ) then
    call p17_start ( n, x )
  else if ( problem == 18 ) then
    call p18_start ( n, x )
  else if ( problem == 19 ) then
    call p19_start ( n, x )
  else if ( problem == 20 ) then
    call p20_start ( n, x )
  else if ( problem == 21 ) then
    call p21_start ( n, x )
  else if ( problem == 22 ) then
    call p22_start ( n, x )
  else if ( problem == 23 ) then
    call p23_start ( n, x )
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_START - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal problem number = ', problem
    stop

  end if

  return
end
subroutine p00_title ( problem, title )

!*****************************************************************************80
!
!! P00_TITLE returns the title of the problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer ( kind = 4 ) problem
  character ( len = * ) title

  if ( problem == 1 ) then
    call p01_title ( title )
  else if ( problem == 2 ) then
    call p02_title ( title )
  else if ( problem == 3 ) then
    call p03_title ( title )
  else if ( problem == 4 ) then
    call p04_title ( title )
  else if ( problem == 5 ) then
    call p05_title ( title )
  else if ( problem == 6 ) then
    call p06_title ( title )
  else if ( problem == 7 ) then
    call p07_title ( title )
  else if ( problem == 8 ) then
    call p08_title ( title )
  else if ( problem == 9 ) then
    call p09_title ( title )
  else if ( problem == 10 ) then
    call p10_title ( title )
  else if ( problem == 11 ) then
    call p11_title ( title )
  else if ( problem == 12 ) then
    call p12_title ( title )
  else if ( problem == 13 ) then
    call p13_title ( title )
  else if ( problem == 14 ) then
    call p14_title ( title )
  else if ( problem == 15 ) then
    call p15_title ( title )
  else if ( problem == 16 ) then
    call p16_title ( title )
  else if ( problem == 17 ) then
    call p17_title ( title )
  else if ( problem == 18 ) then
    call p18_title ( title )
  else if ( problem == 19 ) then
    call p19_title ( title )
  else if ( problem == 20 ) then
    call p20_title ( title )
  else if ( problem == 21 ) then
    call p21_title ( title )
  else if ( problem == 22 ) then
    call p22_title ( title )
  else if ( problem == 23 ) then
    call p23_title ( title )
  else
    title = ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p01_fx ( n, x, f )

!*****************************************************************************80
!
!! P01_FX evaluates the function for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  f(1) = 1.0D+00 - x(1)
  f(2:n) = 10.0D+00 * ( x(2:n) - x(1:n-1) * x(1:n-1) )

  return
end
subroutine p01_n ( n )

!*****************************************************************************80
!
!! P01_N returns the number of equations for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = -2

  return
end
subroutine p01_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P01_JAC sets the jacobian for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point where the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fjac(1:n,1:n) = 0.0D+00

  fjac(1,1) = - 1.0D+00

  do i = 2, n
    fjac(i,i-1) = - 20.0D+00 * x(i-1)
    fjac(i,i) = 10.0D+00
  end do

  return
end
subroutine p01_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P01_SOL returns the solution of problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1:n) = 1.0D+00

  return
end
subroutine p01_start ( n, x )

!*****************************************************************************80
!
!! P01_START specifies a standard approximate solution for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1) = - 1.2D+00
  x(2:n) = 1.0D+00

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns the title of problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Generalized Rosenbrock function.'

  return
end
subroutine p02_fx ( n, x, f )

!*****************************************************************************80
!
!! P02_FX evaluates the function for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)

  f(1) = x(1) + 10.0D+00 * x(2)
  f(2) = sqrt ( 5.0D+00 ) * ( x(3) - x(4) )
  f(3) = ( x(2) - 2.0D+00 * x(3) )**2
  f(4) = sqrt ( 10.0D+00 ) * ( x(1) - x(4) ) * ( x(1) - x(4) )

  return
end
subroutine p02_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P02_JAC sets the jacobian for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) x(n)

  fjac(1,1) = 1.0D+00
  fjac(1,2) = 10.0D+00
  fjac(1,3) = 0.0D+00
  fjac(1,4) = 0.0D+00

  fjac(2,1) = 0.0D+00
  fjac(2,2) = 0.0D+00
  fjac(2,3) = sqrt ( 5.0D+00 )
  fjac(2,4) = - sqrt ( 5.0D+00 )

  fjac(3,1) = 0.0D+00
  fjac(3,2) = 2.0D+00 * ( x(2) - 2.0D+00 * x(3) )
  fjac(3,3) = - 4.0D+00 * ( x(2) - 2.0D+00 * x(3) )
  fjac(3,4) = 0.0D+00

  fjac(4,1) = 2.0D+00 * sqrt ( 10.0D+00 ) * ( x(1) - x(4) )
  fjac(4,2) = 0.0D+00
  fjac(4,3) = 0.0D+00
  fjac(4,4) = - 2.0D+00 * sqrt ( 10.0D+00 ) * ( x(1) - x(4) )

  return
end
subroutine p02_n ( n )

!*****************************************************************************80
!
!! P02_N returns the number of equations for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = 4

  return
end
subroutine p02_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P02_SOL returns the solution of problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1:n) = 0.0D+00

  return
end
subroutine p02_start ( n, x )

!*****************************************************************************80
!
!! P02_START specifies a standard approximate solution for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:4) =   (/ 3.0D+00, -1.0D+00, 0.0D+00, 1.0D+00 /)

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns the title of problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Powell singular function.'

  return
end
subroutine p03_fx ( n, x, f )

!*****************************************************************************80
!
!! P03_FX evaluates the function for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)

  f(1) = 10000.0D+00 * x(1) * x(2) - 1.0D+00
  f(2) = exp ( - x(1) ) + exp ( - x(2) ) - 1.0001D+00

  return
end
subroutine p03_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P03_JAC sets the jacobian for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) x(n)

  fjac(1,1) = 10000.0D+00 * x(2)
  fjac(1,2) = 10000.0D+00 * x(1)

  fjac(2,1) = - exp ( - x(1) )
  fjac(2,2) = - exp ( - x(2) )

  return
end
subroutine p03_n ( n )

!*****************************************************************************80
!
!! P03_N returns the number of equations for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p03_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P03_SOL returns the solution of problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1:2) = (/ 1.098159D-05, 9.106146D+00 /)

  return
end
subroutine p03_start ( n, x )

!*****************************************************************************80
!
!! P03_START specifies a standard approximate solution for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ 0.0D+00, 1.0D+00 /)

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns the title of problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Powell badly scaled function.'

  return
end
subroutine p04_fx ( n, x, f )

!*****************************************************************************80
!
!! P04_FX evaluates the function for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  real ( kind = 8 ) x(n)

  temp1 = x(2) - x(1) * x(1)
  temp2 = x(4) - x(3) * x(3)

  f(1) = - 200.0D+00 * x(1) * temp1 - ( 1.0D+00 - x(1) )

  f(2) = 200.0D+00 * temp1 + 20.2D+00 * ( x(2) - 1.0D+00 ) &
    + 19.8D+00 * ( x(4) - 1.0D+00 )

  f(3) = - 180.0D+00 * x(3) * temp2 - ( 1.0D+00 - x(3) )

  f(4) = 180.0D+00 * temp2 + 20.2D+00 * ( x(4) - 1.0D+00 ) &
    + 19.8D+00 * ( x(2) - 1.0D+00 )

  return
end
subroutine p04_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P04_JAC sets the jacobian for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) x(n)

  fjac(1,1) = - 200.0D+00 * ( x(2) - 3.0D+00 * x(1)**2 ) + 1.0D+00
  fjac(1,2) = - 200.0D+00 * x(1)
  fjac(1,3) = 0.0D+00
  fjac(1,4) = 0.0D+00

  fjac(2,1) = - 400.0D+00 * x(1)
  fjac(2,2) =   220.2D+00
  fjac(2,3) =     0.0D+00
  fjac(2,4) =    19.8D+00

  fjac(3,1) = 0.0D+00
  fjac(3,2) = 0.0D+00
  fjac(3,3) = - 180.0D+00 * ( x(4) - 3.0D+00 * x(3)**2 ) + 1.0D+00
  fjac(3,4) = - 180.0D+00 * x(3)

  fjac(4,1) =     0.0D+00
  fjac(4,2) =    19.8D+00
  fjac(4,3) = - 360.0D+00 * x(3)
  fjac(4,4) =   300.2D+00

  return
end
subroutine p04_n ( n )

!*****************************************************************************80
!
!! P04_N returns the number of equations for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = 4

  return
end
subroutine p04_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P04_SOL returns the solution of problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1:4) = 1.0D+00

  return
end
subroutine p04_start ( n, x )

!*****************************************************************************80
!
!! P04_START specifies a standard approximate solution for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:4) = (/ -3.0D+00, -1.0D+00, -3.0D+00, -1.0D+00 /)

  return
end
subroutine p04_title ( title )

!*****************************************************************************80
!
!! P04_TITLE returns the title of problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Wood function.'

  return
end
subroutine p05_fx ( n, x, f )

!*****************************************************************************80
!
!! P05_FX evaluates the function for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) temp
  real ( kind = 8 ) x(n)

  if ( 0.0D+00 < x(1) ) then
    temp = atan ( x(2) / x(1) ) / ( 2.0D+00 * pi )
  else if (x(1) < 0.0D+00 ) then
    temp = atan ( x(2) / x(1) ) / ( 2.0D+00 * pi ) + 0.5D+00
  else
    temp = sign ( 0.25D+00, x(2) )
  end if

  f(1) = 10.0D+00 * ( x(3) - 10.0D+00 * temp )

  f(2) = 10.0D+00 * ( sqrt ( x(1)**2 + x(2)**2 ) - 1.0D+00 )

  f(3) = x(3)

  return
end
subroutine p05_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P05_JAC sets the jacobian for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  fjac(1,1) =  100.0D+00 * x(2) / ( 2.0D+00 * pi * ( x(1)**2 + x(2)**2 ) )
  fjac(1,2) = -100.0D+00 * x(1) / ( 2.0D+00 * pi * ( x(1)**2 + x(2)**2 ) )
  fjac(1,3) =   10.0D+00

  fjac(2,1) = 10.0D+00 * x(1) / sqrt ( x(1)**2 + x(2)**2 )
  fjac(2,2) = 10.0D+00 * x(2) / sqrt ( x(1)**2 + x(2)**2 )
  fjac(2,3) = 0.0D+00

  fjac(3,1) = 0.0D+00
  fjac(3,2) = 0.0D+00
  fjac(3,3) = 1.0D+00

  return
end
subroutine p05_n ( n )

!*****************************************************************************80
!
!! P05_N returns the number of equations for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = 3

  return
end
subroutine p05_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P05_SOL returns the solution of problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1:3) = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)

  return
end
subroutine p05_start ( n, x )

!*****************************************************************************80
!
!! P05_START specifies a standard approximate solution for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:3) = (/ -1.0D+00, 0.0D+00, 0.0D+00 /)

  return
end
subroutine p05_title ( title )

!*****************************************************************************80
!
!! P05_TITLE returns the title of problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Helical valley function.'

  return
end
subroutine p06_fx ( n, x, f )

!*****************************************************************************80
!
!! P06_FX evaluates the function for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp
  real ( kind = 8 ) ti
  real ( kind = 8 ) x(n)

  f(1:n) = 0.0D+00

  do i = 1, 29

    ti = real ( i, kind = 8 ) / 29.0D+00

    sum1 = 0.0D+00
    temp = 1.0D+00
    do j = 2, n
      sum1 = sum1 + real ( j - 1, kind = 8 ) * temp * x(j)
      temp = ti * temp
    end do

    sum2 = 0.0D+00
    temp = 1.0D+00
    do j = 1, n
      sum2 = sum2 + temp * x(j)
      temp = ti * temp
    end do

    temp = 1.0D+00 / ti

    do k = 1, n
      f(k) = f(k) + temp * ( sum1 - sum2 * sum2 - 1.0D+00 ) &
        * ( real ( k - 1, kind = 8 ) - 2.0D+00 * ti * sum2 )
      temp = ti * temp
    end do

  end do

  f(1) = f(1) + 3.0D+00 * x(1) - 2.0D+00 * x(1) * x(2) + 2.0D+00 * x(1)**3
  f(2) = f(2) + x(2) - x(1)**2 - 1.0D+00

  return
end
subroutine p06_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P06_JAC sets the jacobian for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp
  real ( kind = 8 ) temp1
  real ( kind = 8 ) ti
  real ( kind = 8 ) tj
  real ( kind = 8 ) tk
  real ( kind = 8 ) x(n)

  fjac(1:n,1:n) = 0.0D+00

  do i = 1, 29

    ti = real ( i, kind = 8 ) / 29.0D+00

    sum1 = 0.0D+00
    temp = 1.0D+00
    do j = 2, n
      sum1 = sum1 + real ( j - 1, kind = 8 ) * temp * x(j)
      temp = ti * temp
    end do

    sum2 = 0.0D+00
    temp = 1.0D+00
    do j = 1, n
      sum2 = sum2 + temp * x(j)
      temp = ti * temp
    end do

    temp1 = 2.0D+00 * ( sum1 - sum2 * sum2 - 1.0D+00 )
    tk = 1.0D+00

    do k = 1, n
      tj = tk
      do j = k, n
        fjac(k,j) = fjac(k,j) &
          + tj * ( ( real ( k - 1, kind = 8 ) / ti - 2.0D+00 * sum2 ) &
          * ( real ( j - 1, kind = 8 ) / ti - 2.0D+00 * sum2 ) - temp1 )
        tj = ti * tj
      end do
      tk = ti**2 * tk
    end do

  end do

  fjac(1,1) = fjac(1,1) + 3.0D+00 - 2.0D+00 * x(2) + 6.0D+00 * x(1)**2
  fjac(1,2) = fjac(1,2) - 2.0D+00 * x(1)

  fjac(2,1) = fjac(2,1) - 2.0D+00 * x(1)
  fjac(2,2) = fjac(2,2) + 1.0D+00

  do k = 1, n
    fjac(k:n,k) = fjac(k,k:n)
  end do

  return
end
subroutine p06_n ( n )

!*****************************************************************************80
!
!! P06_N returns the number of equations for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 2

  return
end
subroutine p06_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P06_SOL returns the solution of problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 0

  x(1:n) = 0.0D+00

  return
end
subroutine p06_start ( n, x )

!*****************************************************************************80
!
!! P06_START specifies a standard approximate solution for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = 0.0D+00

  return
end
subroutine p06_title ( title )

!*****************************************************************************80
!
!! P06_TITLE returns the title of problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Watson function.'

  return
end
subroutine p07_fx ( n, x, f )

!*****************************************************************************80
!
!! P07_FX evaluates the function for problem 7.
!
!  Discussion:
!
!    The Chebyquad function is related to the computation of the
!    abscissas for Chebyshev quadrature.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) x(n)

  f(1:n) = 0.0D+00

  do j = 1, n

    t1 = 1.0D+00
    t2 = x(j)

    do i = 1, n

      f(i) = f(i) + t2

      t3 = 2.0D+00 * x(j) * t2 - t1
      t1 = t2
      t2 = t3

    end do

  end do

  f(1:n) = f(1:n) / real ( n, kind = 8 )

  do i = 1, n

    if ( mod ( i, 2 ) == 0 ) then
      f(i) = f(i) + 1.0D+00 / real ( i * i - 1, kind = 8 )
    end if

  end do

  return
end
subroutine p07_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P07_JAC sets the jacobian for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ) t5
  real ( kind = 8 ) t6
  real ( kind = 8 ) x(n)

  do j = 1, n

    t1 = 1.0D+00
    t2 = x(j)

    t4 = 0.0D+00
    t5 = 1.0D+00

    do i = 1, n

      fjac(i,j) = t5

      t6 = 2.0D+00 * t2 + 2.0D+00 * t5 * x(j) - t4
      t4 = t5
      t5 = t6

      t3 = 2.0D+00 * x(j) * t2 - t1
      t1 = t2
      t2 = t3

    end do

  end do

  fjac(1:n,1:n) = fjac(1:n,1:n) / real ( n, kind = 8 )

  return
end
subroutine p07_n ( n )

!*****************************************************************************80
!
!! P07_N returns the number of equations for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n
!
!  Only the values N = 1, 2, 3, 4, 5, 6, 7 and 9 may be used.
!
  n = 9

  return
end
subroutine p07_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P07_SOL returns the solution of problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then
    iknow = 1
    x(1) = 0.5000000000000000D+00
  else if ( n == 2 ) then
    iknow = 1
    x(1) = 0.2113248654051871D+00
    x(2) = 0.7886751345948129D+00
  else if ( n == 3 ) then
    iknow = 1
    x(1) = 0.1464466094067263D+00
    x(2) = 0.5000000000000000D+00
    x(3) = 0.8535533905932737D+00
  else if ( n == 4 ) then
    iknow = 1
    x(1) = 0.1026727638500000D+00
    x(2) = 0.4062037629500000D+00
    x(3) = 0.5937962370500000D+00
    x(4) = 0.8973272361500000D+00
  else if ( n == 5 ) then
    iknow = 1
    x(1) = 8.3751256499999982D-02
    x(2) = 0.3127292952000000D+00
    x(3) = 0.5000000000000000D+00
    x(4) = 0.6872707048000000D+00
    x(5) = 0.9162487435000000D+00
  else if ( n == 6 ) then
    iknow = 1
    x(1) = 6.6876590949999981D-02
    x(2) = 0.2887406731000000D+00
    x(3) = 0.3666822992500000D+00
    x(4) = 0.6333177007499999D+00
    x(5) = 0.7112593269000000D+00
    x(6) = 0.9331234090500000D+00
  else if ( n == 7 ) then
    iknow = 1
    x(1) = 5.8069149599999981D-02
    x(2) = 0.2351716123500000D+00
    x(3) = 0.3380440947500000D+00
    x(4) = 0.5000000000000000D+00
    x(5) = 0.6619559052500000D+00
    x(6) = 0.7648283876499999D+00
    x(7) = 0.9419308504000000D+00
  else if ( n == 9 ) then
    iknow = 1
    x(1) = 4.4205346149999991D-02
    x(2) = 0.1994906723000000D+00
    x(3) = 0.2356191084500000D+00
    x(4) = 0.4160469079000000D+00
    x(5) = 0.5000000000000000D+00
    x(6) = 0.5839530921000000D+00
    x(7) = 0.7643808915500000D+00
    x(8) = 0.8005093276999999D+00
    x(9) = 0.9557946538500000D+00
  else
    iknow = - 1
    x(1:n) = 0.5D+00
  end if

  x(1:n) = 2.0D+00 * x(1:n) - 1.0D+00

  return
end
subroutine p07_start ( n, x )

!*****************************************************************************80
!
!! P07_START specifies a standard approximate solution for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = real ( 2 * i - 1 - n, kind = 8 ) / real ( n + 1, kind = 8 )
  end do

  return
end
subroutine p07_title ( title )

!*****************************************************************************80
!
!! P07_TITLE returns the title of problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Chebyquad function.'

  return
end
subroutine p08_fx ( n, x, f )

!*****************************************************************************80
!
!! P08_FX evaluates the function for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)

  f(1:n-1) = x(1:n-1) + sum ( x(1:n) ) - real ( n + 1, kind = 8 )

  f(n) = product ( x(1:n) ) - 1.0D+00

  return
end
subroutine p08_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P08_JAC sets the jacobian for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) prod
  real ( kind = 8 ) x(n)

  do j = 1, n
    fjac(1:n-1,j) = 1.0D+00
  end do

  do i = 1, n - 1
    fjac(i,i) = 2.0D+00
  end do
!
!  Last row:
!
  do j = 1, n
    prod = 1.0D+00
    do k = 1, n
      if ( k /= j ) then
        prod = x(k) * prod
      end if
    end do
    fjac(n,j) = prod
  end do

  return
end
subroutine p08_n ( n )

!*****************************************************************************80
!
!! P08_N returns the number of equations for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p08_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P08_SOL returns the solution of problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1:n) = 1.0D+00

  return
end
subroutine p08_start ( n, x )

!*****************************************************************************80
!
!! P08_START specifies a standard approximate solution for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = 0.5D+00

  return
end
subroutine p08_title ( title )

!*****************************************************************************80
!
!! P08_TITLE returns the title of problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Brown almost linear function.'

  return
end
subroutine p09_fx ( n, x, f )

!*****************************************************************************80
!
!! P09_FX evaluates the function for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) h
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  h = 1.0D+00 / real ( n + 1, kind = 8 )

  do k = 1, n

    f(k) = 2.0D+00 * x(k) + 0.5D+00 * h * h &
      * ( x(k) + real ( k, kind = 8 ) * h + 1.0D+00 )**3

    if ( 1 < k ) then
      f(k) = f(k) - x(k-1)
    end if

    if ( k < n ) then
      f(k) = f(k) - x(k+1)
    end if

  end do

  return
end
subroutine p09_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P09_JAC sets the jacobian for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fjac(1:n,1:n) = 0.0D+00

  do i = 1, n

    fjac(i,i) = 2.0D+00 + 1.5D+00 * ( x(i) + 1.0D+00 + real ( i, kind = 8 ) &
      / real ( n + 1, kind = 8 )  )**2 / real ( n + 1, kind = 8 )**2

    if ( 1 < i ) then
      fjac(i,i-1) = - 1.0D+00
    end if

    if ( i < n ) then
      fjac(i,i+1) = - 1.0D+00
   end if

  end do

  return
end
subroutine p09_n ( n )

!*****************************************************************************80
!
!! P09_N returns the number of equations for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p09_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P09_SOL returns the solution of problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 0

  x(1:n) = 0.0D+00

  return
end
subroutine p09_start ( n, x )

!*****************************************************************************80
!
!! P09_START specifies a standard approximate solution for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = real ( i * ( i - n - 1 ), kind = 8 ) / real ( n + 1, kind = 8 )**2
  end do

  return
end
subroutine p09_title ( title )

!*****************************************************************************80
!
!! P09_TITLE returns the title of problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Discrete boundary value function.'

  return
end
subroutine p10_fx ( n, x, f )

!*****************************************************************************80
!
!! P10_FX evaluates the function for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) tj
  real ( kind = 8 ) tk
  real ( kind = 8 ) x(n)

  h = 1.0D+00 / real ( n + 1, kind = 8 )

  do k = 1, n

    tk = real ( k, kind = 8 ) / real ( n + 1, kind = 8 )

    sum1 = 0.0D+00
    do j = 1, k
      tj = real ( j, kind = 8 ) * h
      sum1 = sum1 + tj * ( x(j) + tj + 1.0D+00 )**3
    end do

    sum2 = 0.0D+00
    do j = k+1, n
      tj = real ( j, kind = 8 ) * h
      sum2 = sum2 + ( 1.0D+00 - tj ) * ( x(j) + tj + 1.0D+00 )**3
    end do

    f(k) = x(k) + h * ( ( 1.0D+00 - tk ) * sum1 + tk * sum2 ) / 2.0D+00

  end do

  return
end
subroutine p10_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P10_JAC sets the jacobian for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  real ( kind = 8 ) ti
  real ( kind = 8 ) tj
  real ( kind = 8 ) x(n)

  do i = 1, n

    ti = real ( i, kind = 8 ) / real ( n + 1, kind = 8 )

    do j = 1, n
      tj = real ( j, kind = 8 ) / real ( n + 1, kind = 8 )
      temp1 = ( x(j) + tj + 1.0D+00 )**2
      temp2 = min ( ti, tj ) - ti * tj
      fjac(i,j) = 1.5D+00 * temp2 * temp1 / real ( n + 1, kind = 8 )
    end do

    fjac(i,i) = fjac(i,i) + 1.0D+00

  end do

  return
end
subroutine p10_n ( n )

!*****************************************************************************80
!
!! P10_N returns the number of equations for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p10_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P10_SOL returns the solution of problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 0

  x(1:n) = 0.0D+00

  return
end
subroutine p10_start ( n, x )

!*****************************************************************************80
!
!! P10_START specifies a standard approximate solution for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = real ( i * ( i - n - 1 ), kind = 8 ) / real ( n + 1, kind = 8 )**2
  end do

  return
end
subroutine p10_title ( title )

!*****************************************************************************80
!
!! P10_TITLE returns the title of problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Discrete integral equation function.'

  return
end
subroutine p11_fx ( n, x, f )

!*****************************************************************************80
!
!! P11_FX evaluates the function for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c_sum
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  c_sum = sum ( cos ( x(1:n) ) )

  do k = 1, n
    f(k) =  real ( n, kind = 8 ) &
      - c_sum + real ( k, kind = 8 ) * ( 1.0D+00 - cos ( x(k) ) ) &
      - sin ( x(k) )
  end do

  return
end
subroutine p11_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P11_JAC sets the jacobian for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  do i = 1, n

    do j = 1, n

      if ( i /= j ) then
        fjac(i,j) = sin ( x(j) )
      else
        fjac(i,j) = real ( j + 1, kind = 8 ) * sin ( x(j) ) - cos ( x(j) )
      end if

    end do

  end do

  return
end
subroutine p11_n ( n )

!*****************************************************************************80
!
!! P11_N returns the number of equations for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p11_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P11_SOL returns the solution of problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 0

  x(1:n) = 0.0D+00

  return
end
subroutine p11_start ( n, x )

!*****************************************************************************80
!
!! P11_START specifies a standard approximate solution for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = 1.0D+00 / real ( n, kind = 8 )

  return
end
subroutine p11_title ( title )

!*****************************************************************************80
!
!! P11_TITLE returns the title of problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Trigonometric function.'

  return
end
subroutine p12_fx ( n, x, f )

!*****************************************************************************80
!
!! P12_FX evaluates the function for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) sum1
  real ( kind = 8 ) x(n)

  sum1 = 0.0D+00
  do j = 1, n
    sum1 = sum1 + real ( j, kind = 8 ) * ( x(j) - 1.0D+00 )
  end do

  do k = 1, n
    f(k) = x(k) - 1.0D+00 + real ( k, kind = 8 ) * sum1 &
      * ( 1.0D+00 + 2.0D+00 * sum1 * sum1 )
  end do

  return
end
subroutine p12_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P12_JAC sets the jacobian for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) sum1
  real ( kind = 8 ) x(n)

  sum1 = 0.0D+00
  do i = 1, n
    sum1 = sum1 + real ( i, kind = 8 ) * ( x(i) - 1.0D+00 )
  end do

  do i = 1, n
    do j = 1, n

      fjac(i,j) = real ( i * j, kind = 8 ) * ( 1.0D+00 + 6.0D+00 * sum1**2 )

      if ( i == j ) then
        fjac(i,j) = fjac(i,j) + 1.0D+00
      end if

    end do
  end do

  return
end
subroutine p12_n ( n )

!*****************************************************************************80
!
!! P12_N returns the number of equations for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p12_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P12_SOL returns the solution of problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1:n) = 1.0D+00

  return
end
subroutine p12_start ( n, x )

!*****************************************************************************80
!
!! P12_START specifies a standard approximate solution for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = 1.0D+00 - real ( i, kind = 8 ) / real ( n, kind = 8 )
  end do

  return
end
subroutine p12_title ( title )

!*****************************************************************************80
!
!! P12_TITLE returns the title of problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Variably dimensioned function.'

  return
end
subroutine p13_fx ( n, x, f )

!*****************************************************************************80
!
!! P13_FX evaluates the function for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  do k = 1, n

    f(k) = ( 3.0D+00 - 2.0D+00 * x(k) ) * x(k) + 1.0D+00

    if ( 1 < k ) then
      f(k) = f(k) - x(k-1)
    end if

    if ( k < n ) then
      f(k) = f(k) - 2.0D+00 * x(k+1)
    end if

  end do

  return
end
subroutine p13_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P13_JAC sets the jacobian for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  fjac(1:n,1:n) = 0.0D+00

  do k = 1, n

    fjac(k,k) = 3.0D+00 - 4.0D+00 * x(k)

    if ( 1 < k ) then
      fjac(k,k-1) = - 1.0D+00
    end if

    if ( k < n ) then
      fjac(k,k+1) = - 2.0D+00
    end if

  end do

  return
end
subroutine p13_n ( n )

!*****************************************************************************80
!
!! P13_N returns the number of equations for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p13_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P13_SOL returns the solution of problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 0

  x(1:n) = 0.0D+00

  return
end
subroutine p13_start ( n, x )

!*****************************************************************************80
!
!! P13_START specifies a standard approximate solution for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = - 1.0D+00

  return
end
subroutine p13_title ( title )

!*****************************************************************************80
!
!! P13_TITLE returns the title of problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Broyden tridiagonal function.'

  return
end
subroutine p14_fx ( n, x, f )

!*****************************************************************************80
!
!! P14_FX evaluates the function for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  real ( kind = 8 ) temp
  real ( kind = 8 ) x(n)

  ml = 5
  mu = 1

  do k = 1, n

    k1 = max ( 1, k - ml )
    k2 = min ( n, k + mu )

    temp = 0.0D+00
    do j = k1, k2
      if ( j /= k ) then
        temp = temp + x(j) * ( 1.0D+00 + x(j) )
      end if
    end do

    f(k) = x(k) * ( 2.0D+00 + 5.0D+00 * x(k) * x(k) ) + 1.0D+00 - temp

  end do

  return
end
subroutine p14_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P14_JAC sets the jacobian for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  real ( kind = 8 ) x(n)

  fjac(1:n,1:n) = 0.0D+00

  ml = 5
  mu = 1

  do k = 1, n

    k1 = max ( 1, k - ml )
    k2 = min ( n, k + mu )

    do j = k1, k2
      if ( j /= k ) then
        fjac(k,j) = - ( 1.0D+00 + 2.0D+00 * x(j) )
      else
        fjac(k,j) = 2.0D+00 + 15.0D+00 * x(k) * x(k)
      end if
    end do

  end do

  return
end
subroutine p14_n ( n )

!*****************************************************************************80
!
!! P14_N returns the number of equations for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = - 1

  return
end
subroutine p14_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P14_SOL returns the solution of problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 0

  x(1:n) = 0.0D+00

  return
end
subroutine p14_start ( n, x )

!*****************************************************************************80
!
!! P14_START specifies a standard approximate solution for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = - 1.0D+00

  return
end
subroutine p14_title ( title )

!*****************************************************************************80
!
!! P14_TITLE returns the title of problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Broyden banded function.'

  return
end
subroutine p15_fx ( n, x, f )

!*****************************************************************************80
!
!! P15_FX evaluates the function for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)

  f(1) = ( x(1) * x(1) + x(2) * x(3) ) - 0.0001D+00
  f(2) = ( x(1) * x(2) + x(2) * x(4) ) - 1.0D+00
  f(3) = ( x(3) * x(1) + x(4) * x(3) ) - 0.0D+00
  f(4) = ( x(3) * x(2) + x(4) * x(4) ) - 0.0001D+00

  return
end
subroutine p15_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P15_JAC sets the jacobian for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) x(n)

  fjac(1,1) = 2.0D+00 * x(1)
  fjac(1,2) = x(3)
  fjac(1,3) = x(2)
  fjac(1,4) = 0.0D+00

  fjac(2,1) = x(2)
  fjac(2,2) = x(1) + x(4)
  fjac(2,3) = 0.0D+00
  fjac(2,4) = x(2)

  fjac(3,1) = x(3)
  fjac(3,2) = 0.0D+00
  fjac(3,3) = x(1) + x(4)
  fjac(3,4) = x(3)

  fjac(4,1) = 0.0D+00
  fjac(4,2) = x(3)
  fjac(4,3) = x(2)
  fjac(4,4) = 2.0D+00 * x(4)

  return
end
subroutine p15_n ( n )

!*****************************************************************************80
!
!! P15_N returns the number of equations for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = 4

  return
end
subroutine p15_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P15_SOL returns the solution of problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1) = 0.01D+00
  x(2) = 50.0D+00
  x(3) = 0.0D+00
  x(4) = 0.01D+00

  return
end
subroutine p15_start ( n, x )

!*****************************************************************************80
!
!! P15_START specifies a standard approximate solution for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1) = 1.0D+00
  x(2) = 0.0D+00
  x(3) = 0.0D+00
  x(4) = 1.0D+00

  return
end
subroutine p15_title ( title )

!*****************************************************************************80
!
!! P15_TITLE returns the title of problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Hammarling 2 by 2 matrix square root problem.'

  return
end
subroutine p16_fx ( n, x, f )

!*****************************************************************************80
!
!! P16_FX evaluates the function for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)

  f(1) = ( x(1) * x(1) + x(2) * x(4) + x(3) * x(7) ) - 0.0001D+00
  f(2) = ( x(1) * x(2) + x(2) * x(5) + x(3) * x(8) ) - 1.0D+00
  f(3) = ( x(1) * x(3) + x(2) * x(6) + x(3) * x(9) )

  f(4) = ( x(4) * x(1) + x(5) * x(4) + x(6) * x(7) )
  f(5) = ( x(4) * x(2) + x(5) * x(5) + x(6) * x(8) ) - 0.0001D+00
  f(6) = ( x(4) * x(3) + x(5) * x(6) + x(6) * x(9) )

  f(7) = ( x(7) * x(1) + x(8) * x(4) + x(9) * x(7) )
  f(8) = ( x(7) * x(2) + x(8) * x(5) + x(9) * x(8) )
  f(9) = ( x(7) * x(3) + x(8) * x(6) + x(9) * x(9) ) - 0.0001D+00

  return
end
subroutine p16_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P16_JAC sets the jacobian for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) x(n)

  fjac(1:n,1:n) = 0.0D+00

  fjac(1,1) = 2.0D+00 * x(1)
  fjac(1,2) = x(4)
  fjac(1,3) = x(7)
  fjac(1,4) = x(2)
  fjac(1,7) = x(3)

  fjac(2,1) = x(2)
  fjac(2,2) = x(1) + x(5)
  fjac(2,3) = x(8)
  fjac(2,5) = x(2)
  fjac(2,8) = x(3)

  fjac(3,1) = x(3)
  fjac(3,2) = x(6)
  fjac(3,3) = x(1) + x(9)
  fjac(3,6) = x(2)
  fjac(3,9) = x(3)

  fjac(4,1) = x(4)
  fjac(4,4) = x(1) + x(5)
  fjac(4,5) = x(4)
  fjac(4,6) = x(7)
  fjac(4,7) = x(6)

  fjac(5,2) = x(4)
  fjac(5,4) = x(2)
  fjac(5,5) = 2.0D+00 * x(5)
  fjac(5,6) = x(8)
  fjac(5,8) = x(6)

  fjac(6,3) = x(4)
  fjac(6,4) = x(3)
  fjac(6,5) = x(6)
  fjac(6,6) = x(5) + x(9)
  fjac(6,9) = x(6)

  fjac(7,1) = x(7)
  fjac(7,4) = x(8)
  fjac(7,7) = x(1) + x(9)
  fjac(7,8) = x(4)
  fjac(7,9) = x(7)

  fjac(8,2) = x(7)
  fjac(8,5) = x(8)
  fjac(8,7) = x(2)
  fjac(8,8) = x(5) + x(9)
  fjac(8,9) = x(8)

  fjac(9,3) = x(7)
  fjac(9,6) = x(8)
  fjac(9,7) = x(3)
  fjac(9,8) = x(6)
  fjac(9,9) = 2.0D+00 * x(9)

  return
end
subroutine p16_n ( n )

!*****************************************************************************80
!
!! P16_N returns the number of equations for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = 9

  return
end
subroutine p16_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P16_SOL returns the solution of problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1) = 0.01D+00
  x(2) = 50.0D+00
  x(3) = 0.0D+00

  x(4) = 0.0D+00
  x(5) = 0.01D+00
  x(6) = 0.0D+00

  x(7) = 0.0D+00
  x(8) = 0.0D+00
  x(9) = 0.01D+00

  return
end
subroutine p16_start ( n, x )

!*****************************************************************************80
!
!! P16_START specifies a standard approximate solution for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1) = 1.0D+00
  x(2) = 0.0D+00
  x(3) = 0.0D+00
  x(4) = 0.0D+00
  x(5) = 1.0D+00
  x(6) = 0.0D+00
  x(7) = 0.0D+00
  x(8) = 0.0D+00
  x(9) = 1.0D+00

  return
end
subroutine p16_title ( title )

!*****************************************************************************80
!
!! P16_TITLE returns the title of problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Hammarling 3 by 3 matrix square root problem.'

  return
end
subroutine p17_fx ( n, x, f )

!*****************************************************************************80
!
!! P17_FX evaluates the function for problem 17.
!
!  Discussion:
!
!    The equations are:
!
!      F1(X) = X(1) + X(2) - 3
!      F2(X) = X(1)^2 + X(2)^2 - 9
!
!    with roots (3,0) and (0,3).
!
!    Using a starting point of (1,5), here are the iterates for Broyden's
!    method and Newton's method:
!
!    Broyden's Method
!
!    0   1.0              5.0
!    1  -0.625            3.625
!    2  -0.0757575757575  3.0757575757575
!    3  -0.0127942681679  3.0127942681679
!    4  -0.0003138243387  3.0003138243387
!    5  -0.0000013325618  3.0000013325618
!    6  -0.0000000001394  3.0000000001394
!    7   0.0              3.0
!
!    Newton's Method
!
!    0   1.0              5.0
!    1  -0.625            3.625
!    2  -0.0919117647059  3.0919117647059
!    3  -0.0026533419372  3.0026533419372
!    4  -0.0000023425973  3.0000023425973
!    5  -0.0000000000018  3.0000000000018
!    6   0.0              3.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)

  f(1) = x(1) + x(2) - 3.0D+00
  f(2) = x(1)**2 + x(2)**2 - 9.0D+00

  return
end
subroutine p17_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P17_JAC sets the jacobian for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) x(n)

  fjac(1,1) = 1.0D+00
  fjac(1,2) = 1.0D+00

  fjac(2,1) = 2.0D+00 * x(1)
  fjac(2,2) = 2.0D+00 * x(2)

  return
end
subroutine p17_n ( n )

!*****************************************************************************80
!
!! P17_N returns the number of equations for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p17_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P17_SOL returns the solution of problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1) = 0.0D+00
  x(2) = 3.0D+00

  return
end
subroutine p17_start ( n, x )

!*****************************************************************************80
!
!! P17_START specifies a standard approximate solution for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1) = 1.0D+00
  x(2) = 5.0D+00

  return
end
subroutine p17_title ( title )

!*****************************************************************************80
!
!! P17_TITLE returns the title of problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Dennis and Schnabel 2 by 2 example.'

  return
end
subroutine p18_fx ( n, x, f )

!*****************************************************************************80
!
!! P18_FX evaluates the function for problem 18.
!
!  Discussion:
!
!    This problem has as roots any point (x,y) with x or y equal to
!    zero.  The jacobian is of rank one at all roots, except at the
!    origin, where it completely vanishes.
!
!    Newton iterates can converge to the origin, or to points on
!    the X or Y axes, or can even "converge" to a point at infinity,
!    generating a sequence of ever larger points with ever smaller
!    function values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)

  if ( x(1) /= 0.0D+00 ) then
    f(1) = x(2)**2 * ( 1.0D+00 - exp ( - x(1) * x(1) ) ) / x(1)
  else
    f(1) = 0.0D+00
  end if

  if ( x(2) /= 0.0D+00 ) then
    f(2) = x(1) * ( 1.0D+00 - exp ( - x(2)**2 ) ) / x(2)
  else
    f(2) = 0.0D+00
  end if

  return
end
subroutine p18_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P18_JAC sets the jacobian for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) x(n)

  if ( x(1) /= 0.0D+00 ) then
    fjac(1,1) = x(2)**2 * ( 2.0D+00 * exp ( - x(1)**2 ) - &
      ( 1.0D+00 - exp ( - x(1)**2 ) ) / x(1)**2)

    fjac(1,2) = 2.0D+00 * x(2) * ( 1.0D+00 - exp ( - x(1)**2 ) ) / x(1)

  else

    fjac(1,1) = x(2)**2

    fjac(1,2) = 0.0D+00

  end if

  if ( x(2) /= 0.0D+00 ) then

    fjac(2,1) = ( 1.0D+00 - exp ( - x(2)**2 ) ) / x(2)

    fjac(2,2) = x(1) * ( 2.0D+00 * exp ( - x(2)**2 ) - &
      ( 1.0D+00 - exp ( - x(2)**2 ) ) / x(2)**2 )

  else

    fjac(2,1) = 0.0D+00

    fjac(2,2) = x(1)

  end if

  return
end
subroutine p18_n ( n )

!*****************************************************************************80
!
!! P18_N returns the number of equations for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p18_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P18_SOL returns the solution of problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1) = 0.0D+00
  x(2) = 0.0D+00

  return
end
subroutine p18_start ( n, x )

!*****************************************************************************80
!
!! P18_START specifies a standard approximate solution for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1) = 2.0D+00
  x(2) = 2.0D+00

  return
end
subroutine p18_title ( title )

!*****************************************************************************80
!
!! P18_TITLE returns the title of problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Sample problem 18.'

  return
end
subroutine p19_fx ( n, x, f )

!*****************************************************************************80
!
!! P19_FX evaluates the function for problem 19.
!
!  Discussion:
!
!    This problem has a single root at the origin.  Convergence of the
!    Newton iterates should be monotonic, but only linear in rate,
!    since the jacobian is singular at the origin.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)

  f(1) = x(1) * ( x(1)**2 + x(2)**2 )
  f(2) = x(2) * ( x(1)**2 + x(2)**2 )

  return
end
subroutine p19_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P19_JAC sets the jacobian for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) x(n)

  fjac(1,1) = 3.0D+00 * x(1)**2 + x(2)**2
  fjac(1,2) = 2.0D+00 * x(1) * x(2)

  fjac(2,1) = 2.0D+00 * x(1) * x(2)
  fjac(2,2) = x(1)**2 + 3.0D+00 * x(2)**2

  return
end
subroutine p19_n ( n )

!*****************************************************************************80
!
!! P19_N returns the number of equations for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p19_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P19_SOL returns the solution of problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1) = 0.0D+00
  x(2) = 0.0D+00

  return
end
subroutine p19_start ( n, x )

!*****************************************************************************80
!
!! P19_START specifies a standard approximate solution for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ 3.0D+00, 3.0D+00 /)

  return
end
subroutine p19_title ( title )

!*****************************************************************************80
!
!! P19_TITLE returns the title of problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Sample problem 19.'

  return
end
subroutine p20_fx ( n, x, f )

!*****************************************************************************80
!
!! P20_FX evaluates the function for problem 20.
!
!  Discussion:
!
!    This problem has a single root at the origin, and a multiple root
!    at x = 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)

  f(1) = x(1) * ( x(1) - 5.0D+00 ) * ( x(1) - 5.0D+00 )

  return
end
subroutine p20_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P20_JAC sets the jacobian for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) x(n)

  fjac(1,1) = ( 3.0D+00 * x(1) - 5.0D+00 ) * ( x(1) - 5.0D+00 )

  return
end
subroutine p20_n ( n )

!*****************************************************************************80
!
!! P20_N returns the number of equations for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = 1

  return
end
subroutine p20_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P20_SOL returns the solution of problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), save :: icall = -1
  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  icall = icall + 1
  icall = mod ( icall, iknow )

  iknow = 2

  if ( icall == 0 ) then
    x(1) = 0.0D+00
  else if ( icall == 1 ) then
    x(1) = 5.0D+00
  end if

  return
end
subroutine p20_start ( n, x )

!*****************************************************************************80
!
!! P20_START specifies a standard approximate solution for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1) = 1.0D+00

  return
end
subroutine p20_title ( title )

!*****************************************************************************80
!
!! P20_TITLE returns the title of problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Scalar problem f(x) = x * ( x - 5 ) * ( x - 5 ).'

  return
end
subroutine p21_fx ( n, x, f )

!*****************************************************************************80
!
!! P21_FX evaluates the function for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) x(n)

  f(1) = x(1) - x(2)**3 + 5.0D+00 * x(2)**2 -  2.0D+00 * x(2) - 13.0D+00
  f(2) = x(1) + x(2)**3 +           x(2)**2 - 14.0D+00 * x(2) - 29.0D+00

  return
end
subroutine p21_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P21_JAC sets the jacobian for problem 21
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) x(n)

  fjac(1,1) = 1.0D+00
  fjac(1,2) = - 3.0D+00 * x(2)**2 + 10.0D+00 * x(2) -  2.0D+00

  fjac(2,1) = 1.0D+00
  fjac(2,2) =   3.0D+00 * x(2)**2 +  2.0D+00 * x(2) - 14.0D+00

  return
end
subroutine p21_n ( n )

!*****************************************************************************80
!
!! P21_N returns the number of equations for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p21_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P21_SOL returns the solution of problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1) = 5.0D+00
  x(2) = 4.0D+00

  return
end
subroutine p21_start ( n, x )

!*****************************************************************************80
!
!! P21_START specifies a standard approximate solution for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1) =   0.5D+00
  x(2) = - 2.0D+00

  return
end
subroutine p21_title ( title )

!*****************************************************************************80
!
!! P21_TITLE returns the title of problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Freudenstein-Roth function.'

  return
end
subroutine p22_fx ( n, x, f )

!*****************************************************************************80
!
!! P22_FX evaluates the function for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  f(1) = x(1) * x(1) - x(2) + 1.0D+00
  f(2) = x(1) - cos ( 0.5D+00 * pi * x(2) )

  return
end
subroutine p22_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P22_JAC sets the jacobian for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(n,n)
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  fjac(1,1) = 2.0D+00 * x(1)
  fjac(1,2) = - 1.0D+00

  fjac(2,1) = 1.0D+00
  fjac(2,2) = 0.5D+00 * pi * sin ( 0.5D+00 * pi * x(2) )

  return
end
subroutine p22_n ( n )

!*****************************************************************************80
!
!! P22_N returns the number of equations for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = 2

  return
end
subroutine p22_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P22_SOL returns the solution of problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = 1

  x(1) = 0.0D+00
  x(2) = 1.0D+00

  return
end
subroutine p22_start ( n, x )

!*****************************************************************************80
!
!! P22_START specifies a standard approximate solution for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:2) = (/ 1.0D+00, 0.0D+00 /)

  return
end
subroutine p22_title ( title )

!*****************************************************************************80
!
!! P22_TITLE returns the title of problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Boggs function.'

  return
end
subroutine p23_fx ( n, x, f )

!*****************************************************************************80
!
!! P23_FX evaluates the function for problem 23.
!
!    The Chandrasekhar H function is
!
!      F(H)(mu) = H(mu)
!        - 1 / ( 1 - c/2
!        * Integral ( 0 <= nu <= 1 ) [ mu * H(mu) / ( mu + nu ) ] dnu )
!
!    Applying the composite midpoint rule, and setting, for i = 1 to N,
!
!      mu(i) = ( i - 0.5 ) / N
!
!    we have the discrete function
!
!      F(h)(i) = h(i)
!        - 1 / ( 1 - c/(2N)
!        * sum ( 1 <= j <= N ) ( mu(i) * h(j) / ( mu(i) + mu(j) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Subramanyan Chandrasekhar,
!    Radiative Transfer,
!    Dover, 1960,
!    ISBN13: 978-0486605906,
!    LC: QB461.C46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the function.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) C, the constant, which must be between
!    0 and 1.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: c = 0.9D+00
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mu(n)
  real ( kind = 8 ) term
  real ( kind = 8 ) x(n)

  f(1:n) = x(1:n)

  do i = 1, n
    mu(i) = real ( 2 * i - 1, kind = 8 ) / real ( 2 * n, kind = 8 )
  end do

  do i = 1, n

    term = 1.0D+00 - c * sum ( mu(i) * x(1:n) / ( mu(i) + mu(1:n) ) ) &
      / real ( 2 * n,  kind = 8 )

    f(i) = f(i) - 1.0D+00 / term

  end do

  return
end
subroutine p23_n ( n )

!*****************************************************************************80
!
!! P23_N returns the number of equations for problem 23.
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
!  Parameters:
!
!    Output, integer ( kind = 4 ) N.
!    If N is positive, then N is the only value for the number
!    of equations for this problem.
!    If N is negative, then the absolute value of N is the
!    MINIMUM possible value for the number of equations for
!    this problem, and all larger values may also be used.
!
  implicit none

  integer ( kind = 4 ) n

  n = -1

  return
end
subroutine p23_jac ( fjac, n, x )

!*****************************************************************************80
!
!! P23_JAC sets the jacobian for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJAC(N,N), the N by N jacobian matrix.
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input, real ( kind = 8 ) X(N), the point where the jacobian
!    is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: c = 0.9D+00
  real ( kind = 8 ) fjac(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) mu(n)
  real ( kind = 8 ) term
  real ( kind = 8 ) termp
  real ( kind = 8 ) x(n)

  fjac(1:n,1:n) = 0.0D+00

  do i = 1, n
    mu(i) = real ( 2 * i - 1, kind = 8 ) / real ( 2 * n, kind = 8 )
  end do

  do i = 1, n

    term = 1.0D+00 - c * sum ( mu(i) * x(1:n) / ( mu(i) + mu(1:n) ) ) &
      / real ( 2 * n,  kind = 8 )

    do j = 1, n

      termp = c * mu(i) / ( mu(i) + mu(j) ) / real ( 2 * n, kind = 8 )

      fjac(i,j) = ( termp / term ) / term

    end do

  end do

  return
end
subroutine p23_sol ( iknow, n, x )

!*****************************************************************************80
!
!! P23_SOL returns the solution of problem 23.
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
!  Parameters:
!
!    Output, integer IKNOW.
!    -1, there is no solution to the system.
!    0, the solution of the system is not known.
!    positive, IKNOW solutions are known.
!
!    Input, integer ( kind = 4 ) N, the order of the system.  N is primarily
!    needed for problems where N may vary.
!
!    Output, real ( kind = 8 ) X(N), the value of a solution of the system,
!    if known.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iknow
  real ( kind = 8 ) x(n)

  iknow = -1
  x(1:n) = 0.0D+00

  return
end
subroutine p23_start ( n, x )

!*****************************************************************************80
!
!! P23_START specifies a standard approximate solution for problem 23.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the number of equations for the
!    problem.
!
!    Output, real ( kind = 8 ) X(N), the starting point.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) x(n)

  x(1:n) = 1.0D+00

  return
end
subroutine p23_title ( title )

!*****************************************************************************80
!
!! P23_TITLE returns the title of problem 23.
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
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Chandrasekhar function.'

  return
end
function r8vec_norm2 ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM2 returns the 2-norm of a vector.
!
!  Definition:
!
!    The vector 2-norm is defined as:
!
!      value = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose 2-norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM2, the 2-norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm2

  r8vec_norm2 = sqrt ( sum ( a(1:n)**2 ) )

  return
end
subroutine r8ge_fa ( n, a, pivot, info )

!*****************************************************************************80
!
!! R8GE_FA performs a LINPACK style PLU factorization of an R8GE matrix.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!    R8GE_FA is a simplified version of the LINPACK routine SGEFA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N,N), the matrix to be factored.
!    On output, A contains an upper triangular matrix and the multipliers
!    which were used to obtain it.  The factorization can be written
!    A = L * U, where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Output, integer ( kind = 4 ) PIVOT(N), a vector of pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t

  info = 0

  do k = 1, n-1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k+1, n
      if ( abs ( a(l,k) ) < abs ( a(i,k) ) ) then
        l = i
      end if
    end do

    pivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      t      = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = -a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k + 1, n

      if ( l /= k ) then
        t      = a(l,j)
        a(l,j) = a(k,j)
        a(k,j) = t
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  pivot(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop
  end if

  return
end
subroutine r8ge_sl ( n, a_lu, pivot, b, job )

!*****************************************************************************80
!
!! R8GE_SL solves a system factored by R8GE_FA.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!    R8GE_SL is a simplified version of the LINPACK routine SGESL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A_LU(N,N), the LU factors from R8GE_FA.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector from R8GE_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB, specifies the operation.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve PL * Y = B.
!
    do k = 1, n - 1

      l = pivot(k)

      if ( l /= k ) then
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      end if

      b(k+1:n) = b(k+1:n) + a_lu(k+1:n,k) * b(k)

    end do
!
!  Solve U * X = Y.
!
    do k = n, 1, -1
      b(k) = b(k) / a_lu(k,k)
      b(1:k-1) = b(1:k-1) - a_lu(1:k-1,k) * b(k)
    end do
!
!  Solve A' * X = B.
!
  else
!
!  Solve U' * Y = B.
!
    do k = 1, n
      b(k) = ( b(k) - sum ( b(1:k-1) * a_lu(1:k-1,k) ) ) / a_lu(k,k)
    end do
!
!  Solve ( PL )' * X = Y.
!
    do k = n - 1, 1, -1

      b(k) = b(k) + sum ( b(k+1:n) * a_lu(k+1:n,k) )

      l = pivot(k)

      if ( l /= k ) then
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
