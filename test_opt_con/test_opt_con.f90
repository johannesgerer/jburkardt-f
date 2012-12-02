subroutine p00_ab ( problem, m, a, b )

!*****************************************************************************80
!
!! P00_AB returns bounds for a problem specified by index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), lower and upper bounds.
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
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_AB - Fatal error!'
    write ( *, '(a)' ) '  Problem index out of bounds.'
    stop
  end if

  return
end
subroutine p00_f ( problem, m, n, x, f )

!*****************************************************************************80
!
!! P00_F returns the objective function value for a problem specified by index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function values.
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
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_F - Fatal error!'
    write ( *, '(a)' ) '  Problem index out of bounds.'
    stop
  end if

  return
end
subroutine p00_m ( problem, m )

!*****************************************************************************80
!
!! P00_M returns the spatial dimension for a problem specified by index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Output, integer ( kind = 4 ) M, the spatial dimension.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) problem

  if ( problem == 1 ) then
    call p01_m ( m )
  else if ( problem == 2 ) then
    call p02_m ( m )
  else if ( problem == 3 ) then
    call p03_m ( m )
  else if ( problem == 4 ) then
    call p04_m ( m )
  else if ( problem == 5 ) then
    call p05_m ( m )
  else if ( problem == 6 ) then
    call p06_m ( m )
  else if ( problem == 7 ) then
    call p07_m ( m )
  else if ( problem == 8 ) then
    call p08_m ( m )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_M - Fatal error!'
    write ( *, '(a)' ) '  Problem index out of bounds.'
    stop
  end if

  return
end
subroutine p00_problem_num ( problem_num )

!*****************************************************************************80
!
!! P00_PROBLEM_NUM returns the number of problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) PROBLEM_NUM, the number of defined problems.
!
  implicit none

  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num

  problem_num = 8

  return
end
subroutine p00_sol ( problem, m, know, x )

!*****************************************************************************80
!
!! P00_SOL returns known solutions for a problem specified by index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) M, the order of the problem.
!
!    Input/output, integer ( kind = 4 ) KNOW.
!    On input, KNOW is 0, or the index of the previously returned solution.
!    On output, KNOW is 0 if there are no more solutions, or it is the
!    index of the next solution.
!
!    Output, real ( kind = 8 ) X(M), the solution.
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
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SOL - Fatal error!'
    write ( *, '(a)' ) '  Problem index out of bounds.'
    stop
  end if

  return
end
subroutine p00_title ( problem, title )

!*****************************************************************************80
!
!! P00_TITLE returns a title for a problem specified by index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
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
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a)' ) '  Problem number out of bounds.'
    stop
  end if

  return
end
subroutine p01_ab ( m, a, b )

!*****************************************************************************80
!
!! P01_AB returns bounds for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  return
end
subroutine p01_f ( m, n, x, f )

!*****************************************************************************80
!
!! P01_F returns the objective function value for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) p
  real ( kind = 8 ) s
  real ( kind = 8 ) x(m,n)

  do j = 1, n
    p = product ( x(1:m,j) )
    s = sum ( x(1:m,j) )
    f(j) = - exp ( p ) * sin ( s )
  end do

  return
end
subroutine p01_m ( m )

!*****************************************************************************80
!
!! P01_M returns the spatial dimension for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the spatial dimension.
!
  implicit none

  integer ( kind = 4 ) m

  m = 4

  return
end
subroutine p01_sol ( m, know, x )

!*****************************************************************************80
!
!! P01_SOL returns known solutions for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
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
!    Output, real ( kind = 8 ) X(M), the solution.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = (/ &
      0.409887209247642D+00, &
      0.409887209247642D+00, &
      0.409887209247642D+00, &
      0.409887209247642D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns a title for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x) = - exp(prod(x)) * sin(sum(x)).'

  return
end
subroutine p02_ab ( m, a, b )

!*****************************************************************************80
!
!! P02_AB returns bounds for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  return
end
subroutine p02_f ( m, n, x, f )

!*****************************************************************************80
!
!! P02_F returns the objective function value for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) p
  real ( kind = 8 ) s
  real ( kind = 8 ) x(m,n)

  do j = 1, n
    p = x(1,j) * x(2,j)**2 * x(3,j)**3 * x(4,j)**4
    s = sum ( x(1:m,j) )
    f(j) = - exp ( p ) * sin ( s )
  end do

  return
end
subroutine p02_m ( m )

!*****************************************************************************80
!
!! P02_M returns the spatial dimension for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the spatial dimension.
!
  implicit none

  integer ( kind = 4 ) m

  m = 4

  return
end
subroutine p02_sol ( m, know, x )

!*****************************************************************************80
!
!! P02_SOL returns known solutions for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
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
!    Output, real ( kind = 8 ) X(M), the solution.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = (/ &
      0.390500591228663D+00, &
      0.392051909813608D+00, &
      0.393601661544812D+00, &
      0.395149843840982D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns a title for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x) = - exp(x(1)*x(2)^2*x(3)^3*x(4)^4) * sin(sum(x)).'

  return
end
subroutine p03_ab ( m, a, b )

!*****************************************************************************80
!
!! P03_AB returns bounds for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  return
end
subroutine p03_f ( m, n, x, f )

!*****************************************************************************80
!
!! P03_F returns the objective function value for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) s
  real ( kind = 8 ) x(m,n)

  do j = 1, n
    s = - x(1,j) - 2.0D+00 * x(2,j) - 3.0D+00 * x(3,j) - 4.0D+00 * x(4,j)
    f(j) = - 10000.0D+00 * product ( x(1:m,j) ) * exp ( s )
  end do

  return
end
subroutine p03_m ( m )

!*****************************************************************************80
!
!! P03_M returns the spatial dimension for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the spatial dimension.
!
  implicit none

  integer ( kind = 4 ) m

  m = 4

  return
end
subroutine p03_sol ( m, know, x )

!*****************************************************************************80
!
!! P03_SOL returns known solutions for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
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
!    Output, real ( kind = 8 ) X(M), the solution.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = (/ &
      0.999980569087140D+00, &
      0.500000721280566D+00, &
      0.333341891834645D+00, &
      0.249997266604697D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns a title for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x) = -1000 * product(x) * exp(-x(1)-2x(2)-3x(3)-4x(4)).'

  return
end
subroutine p04_ab ( m, a, b )

!*****************************************************************************80
!
!! P04_AB returns bounds for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  return
end
subroutine p04_f ( m, n, x, f )

!*****************************************************************************80
!
!! P04_F returns the objective function value for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,n)

  do j = 1, n
    f(j) = - 100.0D+00 * product ( x(1:m,j) ) * exp ( - x(4,j) ) &
      / ( 1.0D+00 + x(1,j) * x(2,j) * x(3,j) )**2
  end do

  return
end
subroutine p04_m ( m )

!*****************************************************************************80
!
!! P04_M returns the spatial dimension for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the spatial dimension.
!
  implicit none

  integer ( kind = 4 ) m

  m = 4

  return
end
subroutine p04_sol ( m, know, x )

!*****************************************************************************80
!
!! P04_SOL returns known solutions for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
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
!    Output, real ( kind = 8 ) X(M), the solution.
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
!! P04_TITLE returns a title for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x) = -100 * product(x) * exp(-x(4)) / (1+x(1)+x(2)+x(3)).'

  return
end
subroutine p05_ab ( m, a, b )

!*****************************************************************************80
!
!! P05_AB returns bounds for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  return
end
subroutine p05_f ( m, n, x, f )

!*****************************************************************************80
!
!! P05_F returns the objective function value for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,n)

  do j = 1, n
    f(j) = ( x(1,j) -  3.0D+00 / 11.0D+00 )**2 &
         + ( x(2,j) -  6.0D+00 / 13.0D+00 )**2 &
         + ( x(3,j) - 12.0D+00 / 23.0D+00 )**4 &
         + ( x(4,j) -  8.0D+00 / 37.0D+00 )**6
  end do

  return
end
subroutine p05_m ( m )

!*****************************************************************************80
!
!! P05_M returns the spatial dimension for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the spatial dimension.
!
  implicit none

  integer ( kind = 4 ) m

  m = 4

  return
end
subroutine p05_sol ( m, know, x )

!*****************************************************************************80
!
!! P05_SOL returns known solutions for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
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
!    Output, real ( kind = 8 ) X(M), the solution.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = (/ &
       3.0D+00 / 11.0D+00, &
       6.0D+00 / 13.0D+00, &
      12.0D+00 / 23.0D+00, &
       8.0D+00 / 37.0D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p05_title ( title )

!*****************************************************************************80
!
!! P05_TITLE returns a title for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x) = (x(1)-3/11)^2+(x(2)-6/13)^2+(x(3)-12/23)^4+(x(4)-8/37)^6'

  return
end
subroutine p06_ab ( m, a, b )

!*****************************************************************************80
!
!! P06_AB returns bounds for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

  return
end
subroutine p06_f ( m, n, x, f )

!*****************************************************************************80
!
!! P06_F returns the objective function value for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,n)

  do j = 1, n
    arg = &
        1.0D+00 / x(1,j) &
      + 1.0D+00 / x(2,j) &
      + 1.0D+00 / x(3,j) &
      + 1.0D+00 / x(4,j)
    f(j) = - sin ( arg )
  end do

  return
end
subroutine p06_m ( m )

!*****************************************************************************80
!
!! P06_M returns the spatial dimension for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the spatial dimension.
!
  implicit none

  integer ( kind = 4 ) m

  m = 4

  return
end
subroutine p06_sol ( m, know, x )

!*****************************************************************************80
!
!! P06_SOL returns known solutions for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
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
!    Output, real ( kind = 8 ) X(M), the solution.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  if ( know == 0 ) then
    know = 1
    x(1:m) = (/ &
      0.509282516910744D+00, &
      0.509282516910744D+00, &
      0.509282516910746D+00, &
      0.509282516910744D+00 /)
  else
    know = 0
  end if

  return
end
subroutine p06_title ( title )

!*****************************************************************************80
!
!! P06_TITLE returns a title for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Harald Niederreiter, Kevin McCurley,
!    Optimization of functions by quasi-random search methods,
!    Computing,
!    Volume 22, Number 2, 1979, pages 119-123.
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x) = - sin(1/x(1)+1/x(2)+1/x(3)+1/x(4))'

  return
end
subroutine p07_ab ( m, a, b )

!*****************************************************************************80
!
!! P07_AB returns bounds for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = 0.0D+00
  b(1:m) = 10.0D+00

  return
end
subroutine p07_f ( m, n, x, f )

!*****************************************************************************80
!
!! P07_F returns the objective function value for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Langerman10 reference?
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2, 5 ) :: a = reshape ( (/ &
    3.0D+00, 5.0D+00, &
    5.0D+00, 2.0D+00, &
    2.0D+00, 1.0D+00, &
    1.0D+00, 4.0D+00, &
    7.0D+00, 9.0D+00 /), &
    (/ 2, 5 /) )
  real ( kind = 8 ) arg
  real ( kind = 8 ) :: c(5) = (/ &
    1.0D+00, 2.0D+00, 5.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(m,n)

  f(1:n) = 0.0D+00

  do j = 1, n
    do k = 1, 5
      arg = dot_product ( x(1:m,j) - a(1:m,k), x(1:m,j) - a(1:m,k) )
      f(j) = f(j) - c(k) * exp ( - arg / pi ) * cos ( pi * arg )
    end do
  end do

  return
end
subroutine p07_m ( m )

!*****************************************************************************80
!
!! P07_M returns the spatial dimension for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the spatial dimension.
!
  implicit none

  integer ( kind = 4 ) m

  m = 2

  return
end
subroutine p07_sol ( m, know, x )

!*****************************************************************************80
!
!! P07_SOL returns known solutions for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2012
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
!    Output, real ( kind = 8 ) X(M), the solution.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  know = 0

  return
end
subroutine p07_title ( title )

!*****************************************************************************80
!
!! P07_TITLE returns a title for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2012
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

  title = 'f(x) = Langerman2 function'

  return
end
subroutine p08_ab ( m, a, b )

!*****************************************************************************80
!
!! P08_AB returns bounds for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, real ( kind = 8 ) A(M), B(M), lower and upper bounds.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)

  a(1:m) = 0.0D+00
  b(1:m) = 10.0D+00

  return
end
subroutine p08_f ( m, n, x, f )

!*****************************************************************************80
!
!! P08_F returns the objective function value for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Langerman10 reference?
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(M,N), the arguments.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 10, 30 ) :: a = reshape ( (/ &
    9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020, &
    9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374, &
    8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982, &
    2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426, &
    8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567, &
    7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208, &
    1.256, 3.605, 8.623, 6.905, 4.584, 8.133, 6.071, 6.888, 4.187, 5.448, &
    8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762, &
    0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637, &
    7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247, &
    0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016, &
    2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789, &
    8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109, &
    2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564, &
    4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670, &
    8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826, &
    8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591, &
    4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740, &
    2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675, &
    6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258, &
    0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070, &
    5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234, &
    3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027, &
    8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064, &
    1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224, &
    0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644, &
    0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229, &
    4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506, &
    9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732, &
    4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699, 6.500 /), &
    (/ 10, 30 /) )
  real ( kind = 8 ) arg
  real ( kind = 8 ) :: c(30) = (/ &
    0.806, 0.517, 1.500, 0.908, 0.965, &
    0.669, 0.524, 0.902, 0.531, 0.876, &
    0.462, 0.491, 0.463, 0.714, 0.352, &
    0.869, 0.813, 0.811, 0.828, 0.964, &
    0.789, 0.360, 0.369, 0.992, 0.332, &
    0.817, 0.632, 0.883, 0.608, 0.326 /)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(m,n)

  f(1:n) = 0.0D+00

  do j = 1, n
    do k = 1, 30
      arg = dot_product ( x(1:m,j) - a(1:m,k), x(1:m,j) - a(1:m,k) )
      f(j) = f(j) - c(k) * exp ( - arg / pi ) * cos ( pi * arg )
    end do
  end do

  return
end
subroutine p08_m ( m )

!*****************************************************************************80
!
!! P08_M returns the spatial dimension for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, the spatial dimension.
!
  implicit none

  integer ( kind = 4 ) m

  m = 10

  return
end
subroutine p08_sol ( m, know, x )

!*****************************************************************************80
!
!! P08_SOL returns known solutions for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2012
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
!    Output, real ( kind = 8 ) X(M), the solution.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) know
  real ( kind = 8 ) x(m)

  know = 0

  return
end
subroutine p08_title ( title )

!*****************************************************************************80
!
!! P08_TITLE returns a title for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2012
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

  title = 'f(x) = Langerman10 function'

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
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
