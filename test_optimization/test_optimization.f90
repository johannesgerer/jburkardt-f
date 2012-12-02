subroutine p00_ab ( problem, m, a, b )

!*****************************************************************************80
!
!! P00_AB evaluates the limits of the optimization region for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) problem

  if ( problem == 1 ) then
    call p01_ab ( m, a, b )
  else if ( problem == 2 ) then
    call p02_ab ( m, a, b )
  else if ( problem == 3 ) then
    call p03_ab ( m, a, b )
  else if ( problem == 4 ) then
    call p04_ab ( m, a, b )
  else if ( problem == 5 ) then
    call p05_ab ( m, a, b )
  else if ( problem == 6 ) then
    call p06_ab ( m, a, b )
  else if ( problem == 7 ) then
    call p07_ab ( m, a, b )
  else if ( problem == 8 ) then
    call p08_ab ( m, a, b )
  else if ( problem == 9 ) then
    call p09_ab ( m, a, b )
  else if ( problem == 10 ) then
    call p10_ab ( m, a, b )
  else if ( problem == 11 ) then
    call p11_ab ( m, a, b )
  else if ( problem == 12 ) then
    call p12_ab ( m, a, b )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_AB - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_compass_search ( problem, m, x0, delta_tol, delta_init, &
  k_max, x, fx, k )

!*****************************************************************************80
!
!! P00_COMPASS_SEARCH carries out a direct search minimization algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Tamara Kolda, Robert Michael Lewis, Virginia Torczon,
!    Optimization by Direct Search: New Perspectives on Some Classical 
!    and Modern Methods,
!    SIAM Review,
!    Volume 45, Number 3, 2003, pages 385-482. 
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X0(M), a starting estimate for the minimizer.
!
!    Input, real ( kind = 8 ) DELTA_TOL, the smallest step size that is allowed.
!
!    Input, real ( kind = 8 ) DELTA_INIT, the starting stepsize.  
!
!    Input, integer ( kind = 4 ) K_MAX, the maximum number of steps allowed.
!
!    Output, real ( kind = 8 ) X(M), the estimated minimizer.
!
!    Output, real ( kind = 8 ) FX, the function value at X.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: n = 1

  logical decrease
  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_init
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ) fx
  real ( kind = 8 ) fxd
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) problem
  real ( kind = 8 ) s
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)
  real ( kind = 8 ) xd(m)

  k = 0
  x(1:m) = x0(1:m)
  call p00_f ( problem, m, n, x, fx )

  if ( delta_tol <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_COMPASS_SEARCH - Fatal error!'
    write ( *, '(a)' ) '  DELTA_TOL <= 0.0.'
    write ( *, '(a,g14.6)' ) '  DELTA_TOL = ', delta_tol
    stop
  end if

  if ( delta_init <= delta_tol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_COMPASS_SEARCH - Fatal error!'
    write ( *, '(a)' ) '  DELTA_INIT < DELTA_TOL.'
    write ( *, '(a,g14.6)' ) '  DELTA_INIT = ', delta_init
    write ( *, '(a,g14.6)' ) '  DELTA_TOL = ', delta_tol
    stop
  end if

  delta = delta_init

  do while ( k < k_max )

    k = k + 1
!
!  For each coordinate direction I, seek a lower function value
!  by increasing or decreasing X(I) by DELTA.
!
    decrease = .false.
    s = + 1.0D+00
    i = 1

    do ii = 1, 2 * m

      xd = x
      xd(i) = xd(i) + s * delta
      call p00_f ( problem, m, n, xd, fxd )
!
!  As soon as a decrease is noticed, accept the new point.
!
      if ( fxd < fx ) then
        x = xd
        fx = fxd
        decrease = .true.
        exit
      end if

      s = - s
      if ( s == + 1.0D+00 ) then
        i = i + 1
      end if

    end do
!
!  If no decrease occurred, reduce DELTA.
!
    if ( .not. decrease ) then
      delta = delta / 2.0D+00
      if ( delta < delta_tol ) then
        exit
      end if
    end if

  end do

  return
end
subroutine p00_f ( problem, m, n, x, f )

!*****************************************************************************80
!
!! P00_F evaluates the objective function for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the argument of the objective function.
!
!    Output, real ( kind = 8 ) F(N), the objective function evaluated at
!    each argument.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(m,n)

  if ( problem == 1 ) then
    call p01_f ( m, n, x, f )
  else if ( problem == 2 ) then
    call p02_f ( m, n, x, f )
  else if ( problem == 3 ) then
    call p03_f ( m, n, x, f )
  else if ( problem == 4 ) then
    call p04_f ( m, n, x, f )
  else if ( problem == 5 ) then
    call p05_f ( m, n, x, f )
  else if ( problem == 6 ) then
    call p06_f ( m, n, x, f )
  else if ( problem == 7 ) then
    call p07_f ( m, n, x, f )
  else if ( problem == 8 ) then
    call p08_f ( m, n, x, f )
  else if ( problem == 9 ) then
    call p09_f ( m, n, x, f )
  else if ( problem == 10 ) then
    call p10_f ( m, n, x, f )
  else if ( problem == 11 ) then
    call p11_f ( m, n, x, f )
  else if ( problem == 12 ) then
    call p12_f ( m, n, x, f )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_F - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number PROBLEM = ', problem
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
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!   Output, integer ( kind = 4 ) PROBLEM_NUM, the number of problems available.
!
  implicit none

  integer ( kind = 4 ) problem_num

  problem_num = 12

  return
end
subroutine p00_sol ( problem, m, know, x )

!*****************************************************************************80
!
!! P00_SOL returns the solution for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(m)

  if ( problem == 1 ) then
    call p01_sol ( m, know, x )
  else if ( problem == 2 ) then
    call p02_sol ( m, know, x )
  else if ( problem == 3 ) then
    call p03_sol ( m, know, x )
  else if ( problem == 4 ) then
    call p04_sol ( m, know, x )
  else if ( problem == 5 ) then
    call p05_sol ( m, know, x )
  else if ( problem == 6 ) then
    call p06_sol ( m, know, x )
  else if ( problem == 7 ) then
    call p07_sol ( m, know, x )
  else if ( problem == 8 ) then
    call p08_sol ( m, know, x )
  else if ( problem == 9 ) then
    call p09_sol ( m, know, x )
  else if ( problem == 10 ) then
    call p10_sol ( m, know, x )
  else if ( problem == 11 ) then
    call p11_sol ( m, know, x )
  else if ( problem == 12 ) then
    call p12_sol ( m, know, x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SOL - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p00_title ( problem, title )

!*****************************************************************************80
!
!! P00_TITLE returns a title for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the number of the problem.
!
!    Output, character ( len = * ) TITLE, a title for the problem.
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
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of PROBLEM = ', problem
    stop
  end if

  return
end
subroutine p01_ab ( m, a, b )

!*****************************************************************************80
!
!! P01_AB evaluates the limits of the optimization region for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = - 5.0D+00
  b(1:m) = + 5.0D+00

  return
end
subroutine p01_f ( m, n, x, f )

!*****************************************************************************80
!
!! P01_F evaluates the objective function for problem 01.
!
!  Discussion:
!
!    The function is continuous, convex, and unimodal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hugues Bersini, Marco Dorigo, Stefan Langerman, Gregory Seront, 
!    Luca Gambardella,
!    Results of the first international contest on evolutionary optimisation,
!    In Proceedings of 1996 IEEE International Conference on Evolutionary 
!    Computation,
!    IEEE Press, pages 611-615, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function evaluated at the arguments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,n)

  do j = 1, n
    f(j) = sum ( ( x(1:m,j) - 1.0D+00 ) ** 2 )
  end do

  return
end
subroutine p01_sol ( m, know, x )

!*****************************************************************************80
!
!! P01_SOL returns the solution for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = 1.0D+00
  else
    know = 0
  end if

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns a title for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The sphere model.'

  return
end
subroutine p02_ab ( m, a, b )

!*****************************************************************************80
!
!! P02_AB evaluates the limits of the optimization region for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = - 5.12D+00
  b(1:m) = + 5.12D+00

  return
end
subroutine p02_f ( m, n, x, f )

!*****************************************************************************80
!
!! P02_F evaluates the objective function for problem 02.
!
!  Discussion:
!
!    This function is also known as the weighted sphere model.
!
!    The function is continuous, convex, and unimodal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function evaluated at the arguments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,m)
  real ( kind = 8 ) y(m)

  call r8vec_indicator ( m, y )

  do j = 1, n
    f(j) = sum ( y(1:m) * x(1:m,j) ** 2 )
  end do

  return
end
subroutine p02_sol ( m, know, x )

!*****************************************************************************80
!
!! P02_SOL returns the solution for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = 0.0D+00
  else
    know = 0
  end if

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns a title for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The axis-parallel hyper-ellipsoid function.'

  return
end
subroutine p03_ab ( m, a, b )

!*****************************************************************************80
!
!! P03_AB evaluates the limits of the optimization region for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = - 65.536D+00
  b(1:m) = + 65.536D+00

  return
end
subroutine p03_f ( m, n, x, f )

!*****************************************************************************80
!
!! P03_F evaluates the objective function for problem 03.
!
!  Discussion:
!
!    This function is also known as the weighted sphere model.
!
!    The function is continuous, convex, and unimodal.
!
!     There is a typographical error in Molga and Smutnicki, so that the
!     formula for this function is given incorrectly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function evaluated at the arguments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) x_sum

  do j = 1, n

    f(j) = 0.0D+00
    x_sum = 0.0D+00

    do i = 1, m
      x_sum = x_sum + x(i,j)
      f(j) = f(j) + x_sum**2
    end do

  end do

  return
end
subroutine p03_sol ( m, know, x )

!*****************************************************************************80
!
!! P03_SOL returns the solution for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = 0.0D+00
  else
    know = 0
  end if

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns a title for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The rotated hyper-ellipsoid function.'

  return
end
subroutine p04_ab ( m, a, b )

!*****************************************************************************80
!
!! P04_AB evaluates the limits of the optimization region for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = - 2.048D+00
  b(1:m) = + 2.048D+00

  return
end
subroutine p04_f ( m, n, x, f )

!*****************************************************************************80
!
!! P04_F evaluates the objective function for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Howard Rosenbrock,
!    An Automatic Method for Finding the Greatest or Least Value of a Function,
!    Computer Journal,
!    Volume 3, 1960, pages 175-184.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function evaluated at the arguments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,n)

  do j = 1, n
    f(j) = sum ( ( 1.0D+00 - x(1:m,j) )**2 ) &
         + sum ( ( x(2:m,j) - x(1:m-1,j) )**2 )
  end do

  return
end
subroutine p04_sol ( m, know, x )

!*****************************************************************************80
!
!! P04_SOL returns the solution for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = 1.0D+00
  else
    know = 0
  end if

  return
end
subroutine p04_title ( title )

!*****************************************************************************80
!
!! P04_TITLE returns a title for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Rosenbrock''s valley.'

  return
end
subroutine p05_ab ( m, a, b )

!*****************************************************************************80
!
!! P05_AB evaluates the limits of the optimization region for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = - 5.12D+00
  b(1:m) = + 5.12D+00

  return
end
subroutine p05_f ( m, n, x, f )

!*****************************************************************************80
!
!! P05_F evaluates the objective function for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function evaluated at the arguments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(m,n)

  do j = 1, n

    f(j) = real ( 10 * m, kind = 8 )

    do i = 1, m
      f(j) = f(j) + x(i,j) ** 2 - 10.0D+00 * cos ( 2.0D+00 * pi * x(i,j) )
    end do

  end do

  return
end
subroutine p05_sol ( m, know, x )

!*****************************************************************************80
!
!! P05_SOL returns the solution for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = 0.0D+00
  else
    know = 0
  end if

  return
end
subroutine p05_title ( title )

!*****************************************************************************80
!
!! P05_TITLE returns a title for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Rastrigin''s function.'

  return
end
subroutine p06_ab ( m, a, b )

!*****************************************************************************80
!
!! P06_AB evaluates the limits of the optimization region for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = - 500.0D+00
  b(1:m) = + 500.0D+00

  return
end
subroutine p06_f ( m, n, x, f )

!*****************************************************************************80
!
!! P06_F evaluates the objective function for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hans-Paul Schwefel,
!    Numerical optimization of computer models,
!    Wiley, 1981,
!    ISBN13: 978-0471099888,
!    LC: QA402.5.S3813.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function evaluated at the arguments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,n)

  do j = 1, n
    f(j) = - sum ( x(1:m,j) * sin ( sqrt ( abs ( x(1:m,j) ) ) ) )
  end do

  return
end
subroutine p06_sol ( m, know, x )

!*****************************************************************************80
!
!! P06_SOL returns the solution for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = 420.9687D+00
  else
    know = 0
  end if

  return
end
subroutine p06_title ( title )

!*****************************************************************************80
!
!! P06_TITLE returns a title for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Schwefel''s function.'

  return
end
subroutine p07_ab ( m, a, b )

!*****************************************************************************80
!
!! P07_AB evaluates the limits of the optimization region for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = - 600.0D+00
  b(1:m) = + 600.0D+00

  return
end
subroutine p07_f ( m, n, x, f )

!*****************************************************************************80
!
!! P07_F evaluates the objective function for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function evaluated at the arguments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) y(m)

  call r8vec_indicator ( m, y )
  y(1:m) = sqrt ( y(1:m) )

  do j = 1, n
    f(j) = sum ( x(1:m,j) ** 2 ) / 4000.0D+00 &
      - product ( cos ( x(1:m,j) / y(1:m) ) ) + 1.0D+00
  end do

  return
end
subroutine p07_sol ( m, know, x )

!*****************************************************************************80
!
!! P07_SOL returns the solution for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = 0.0D+00
  else
    know = 0
  end if

  return
end
subroutine p07_title ( title )

!*****************************************************************************80
!
!! P07_TITLE returns a title for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Griewank''s function.'

  return
end
subroutine p08_ab ( m, a, b )

!*****************************************************************************80
!
!! P08_AB evaluates the limits of the optimization region for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = - 1.0D+00
  b(1:m) = + 1.0D+00

  return
end
subroutine p08_f ( m, n, x, f )

!*****************************************************************************80
!
!! P08_F evaluates the objective function for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function evaluated at the arguments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) y(m)

  call r8vec_indicator ( m, y )
  y(1:m) = y(1:m) + 1.0D+00

  do j = 1, n
    f(j) = sum ( abs ( x(1:m,j) ) ** y(1:m) )
  end do

  return
end
subroutine p08_sol ( m, know, x )

!*****************************************************************************80
!
!! P08_SOL returns the solution for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = 0.0D+00
  else
    know = 0
  end if

  return
end
subroutine p08_title ( title )

!*****************************************************************************80
!
!! P08_TITLE returns a title for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The power sum function.'

  return
end
subroutine p09_ab ( m, a, b )

!*****************************************************************************80
!
!! P09_AB evaluates the limits of the optimization region for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = - 32.768D+00
  b(1:m) = + 32.768D+00

  return
end
subroutine p09_f ( m, n, x, f )

!*****************************************************************************80
!
!! P09_F evaluates the objective function for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function evaluated at the arguments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: a = 20.0D+00
  real ( kind = 8 ), parameter :: b = 0.2D+00
  real ( kind = 8 ), parameter :: c = 0.2D+00
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(m,n)

  do j = 1, n
    f(j) = - a * exp ( - b * sqrt ( sum ( x(1:m,j)**2 ) &
      / real ( m, kind = 8 ) ) ) &
      - exp ( sum ( cos ( c * pi * x(1:m,j) ) ) / real ( m, kind = 8 ) ) &
      + a + exp ( 1.0D+00 )
  end do

  return
end
subroutine p09_sol ( m, know, x )

!*****************************************************************************80
!
!! P09_SOL returns the solution for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = 0.0D+00
  else
    know = 0
  end if

  return
end
subroutine p09_title ( title )

!*****************************************************************************80
!
!! P09_TITLE returns a title for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Ackley''s function.'

  return
end
subroutine p10_ab ( m, a, b )

!*****************************************************************************80
!
!! P10_AB evaluates the limits of the optimization region for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  a(1:m) = 0.0D+00
  b(1:m) = pi

  return
end
subroutine p10_f ( m, n, x, f )

!*****************************************************************************80
!
!! P10_F evaluates the objective function for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function evaluated at the arguments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: p = 10
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) y(m)

  call r8vec_indicator ( m, y )

  do j = 1, n
    f(j) = - sum ( &
      sin ( x(1:m,j) ) * ( sin ( x(1:m,j)**2 * y(1:m) / pi ) ) ** ( 2 * p ) &
    )
  end do

  return
end
subroutine p10_sol ( m, know, x )

!*****************************************************************************80
!
!! P10_SOL returns the solution for problem 10.
!
!  Discussion:
!
!    The minimum value is - 0.966 * M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  know = 0
  x(1:m) = 0.0D+00

  return
end
subroutine p10_title ( title )

!*****************************************************************************80
!
!! P10_TITLE returns a title for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Michalewicz''s function.'

  return
end
subroutine p11_ab ( m, a, b )

!*****************************************************************************80
!
!! P11_AB evaluates the limits of the optimization region for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = - 5.12D+00
  b(1:m) = + 5.12D+00

  return
end
subroutine p11_f ( m, n, x, f )

!*****************************************************************************80
!
!! P11_F evaluates the objective function for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function evaluated at the arguments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) rsq
  real ( kind = 8 ) x(m,n)

  do j = 1, n

    rsq = sum ( x(1:m,j)**2 )

    f(j) = - ( 1.0D+00 + cos ( 12.0D+00 * sqrt ( rsq ) ) ) &
      / ( 0.5D+00 * rsq + 2.0D+00 )

  end do

  return
end
subroutine p11_sol ( m, know, x )

!*****************************************************************************80
!
!! P11_SOL returns the solution for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    x(1:m) = 0.0D+00
    know = 1
  else
    know = 0
  end if

  return
end
subroutine p11_title ( title )

!*****************************************************************************80
!
!! P11_TITLE returns a title for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Drop wave function.'

  return
end
subroutine p12_ab ( m, a, b )

!*****************************************************************************80
!
!! P12_AB evaluates the limits of the optimization region for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), the lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  return
end
subroutine p12_f ( m, n, x, f )

!*****************************************************************************80
!
!! P12_F evaluates the objective function for problem 12.
!
!  Discussion:
!
!    In dimension I, the function is a piecewise linear function with
!    local minima at 0 and 1.0, and a global minimum at ALPHA(I) = I/(M+1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marcin Molga, Czeslaw Smutnicki,
!    Test functions for optimization needs.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function evaluated at the arguments.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha(m)
  real ( kind = 8 ), parameter :: beta = 2.0D+00
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,n)
!
!  I'm just choosing ALPHA in [0,1] arbitrarily.
!
  do i = 1, m
    alpha(i) = real ( i, kind = 8 ) / real ( m + 1, kind = 8 )
  end do

  do j = 1, n

    f(j) = 0.0D+00

    do i = 1, m

      if ( x(i,j) <= 0.0D+00 ) then
        g = x(i,j)
      else if ( x(i,j) <= 0.8D+00 * alpha(i) ) then
        g = 0.8D+00 - x(i,j) / alpha(i)
      else if ( x(i,j) <= alpha(i) ) then
        g = 5.0D+00 * x(i,j) / alpha(i) - 4.0D+00
      else if ( x(i,j) <= ( 1.0D+00 + 4.0D+00 * alpha(i) ) / 5.0D+00 ) then
        g = 1.0D+00 + 5.0D+00 * ( x(i,j) - alpha(i) ) / ( alpha(i) - 1.0D+00 )
      else if ( x(i,j) <= 1.0D+00 ) then
        g = 0.8D+00 + ( x(i,j) - 1.0D+00 ) / ( 1.0D+00 - alpha(i) )
      else
        g = x(i,j) - 1.0D+00
      end if

      f(j) = f(j) + g

    end do

    f(j) = f(j) / real ( m, kind = 8 )
    f(j) = - ( f(j) ** beta )

  end do

  return
end
subroutine p12_sol ( m, know, x )

!*****************************************************************************80
!
!! P12_SOL returns the solution for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution, if known.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) alpha(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  do i = 1, m
    alpha(i) = real ( i, kind = 8 ) / real ( m + 1, kind = 8 )
  end do

  if ( know == 0 ) then
    know = 1
    x(1:m) = alpha(1:m)
  else
    know = 0
  end if

  return
end
subroutine p12_title ( title )

!*****************************************************************************80
!
!! P12_TITLE returns a title for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'The deceptive function.'

  return
end
subroutine r8col_uniform ( m, n, a, b, seed, r )

!*****************************************************************************80
!
!! R8COL_UNIFORM fills an R8COL with scaled pseudorandom numbers.
!
!  Discussion:
!
!    An R8COL is an array of R8 values, regarded as a set of column vectors.
!
!    The user specifies a minimum and maximum value for each row.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = a(i) &
        + ( b(i) - a(i) ) * real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8vec_indicator ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR sets an R8VEC to the indicator vector.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, real ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = real ( i, kind = 8 )
  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
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
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
