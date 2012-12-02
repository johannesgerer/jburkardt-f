subroutine p00_f ( prob, n, x, value )

!*****************************************************************************80
!
!! P00_F evaluates the function for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the desired test problem.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) prob
  real ( kind = 8 ) value(n)
  real ( kind = 8 ) x(n)

  if ( prob == 1 ) then
    call p01_f ( n, x, value )
  else if ( prob == 2 ) then
    call p02_f ( n, x, value )
  else if ( prob == 3 ) then
    call p03_f ( n, x, value )
  else if ( prob == 4 ) then
    call p04_f ( n, x, value )
  else if ( prob == 5 ) then
    call p05_f ( n, x, value )
  else if ( prob == 6 ) then
    call p06_f ( n, x, value )
  else if ( prob == 7 ) then
    call p07_f ( n, x, value )
  else if ( prob == 8 ) then
    call p08_f ( n, x, value )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_F - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_prob_num ( prob_num )

!*****************************************************************************80
!
!! P00_PROB_NUM returns the number of problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) PROB_NUM, the number of problems.
!
  implicit none

  integer ( kind = 4 ) prob_num

  prob_num = 8

  return
end
subroutine p00_title ( prob, title )

!*****************************************************************************80
!
!! P00_TITLE returns the title of any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the desired test problem.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer ( kind = 4 ) prob
  character ( len = * ) title

  if ( prob == 1 ) then
    call p01_title ( title )
  else if ( prob == 2 ) then
    call p02_title ( title )
  else if ( prob == 3 ) then
    call p03_title ( title )
  else if ( prob == 4 ) then
    call p04_title ( title )
  else if ( prob == 5 ) then
    call p05_title ( title )
  else if ( prob == 6 ) then
    call p06_title ( title )
  else if ( prob == 7 ) then
    call p07_title ( title )
  else if ( prob == 8 ) then
    call p08_title ( title )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p01_f ( n, x, value )

!*****************************************************************************80
!
!! P01_F evaluates the function for problem p01.
!
!  Discussion:
!
!    Step function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) value(n)
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( x(i) <= 1.0D+00 / 3.0D+00 ) then
      value(i) = -1.0D+00
    else if ( x(i) <= 4.0D+00 / 5.0D+00 ) then
      value(i) = 2.0D+00
    else
      value(i) = 1.0D+00
    end if
  end do

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns the title of problem p01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
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

  title = 'f(x) = steps -1/2/1 at [0,1/3], [1/3,4/5], [4/5,1].'

  return
end
subroutine p02_f ( n, x, value )

!*****************************************************************************80
!
!! P02_F evaluates the function for problem p01.
!
!  Discussion:
!
!    Nondifferentiable function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) value(n)
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( x(i) <= 1.0D+00 / 3.0D+00 ) then
      value(i) = 1.0D+00 - 3.0D+00 * x(i)
    else
      value(i) = 6.0D+00 * x(i) - 2.0D+00
    end if
  end do

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns the title of problem p02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
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

  title = 'f(x) = (1-3x), x < 1/3 (6x-2) if 1/3 < x'

  return
end
subroutine p03_f ( n, x, value )

!*****************************************************************************80
!
!! P03_F evaluates the function for problem p03.
!
!  Discussion:
!
!    Step function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) value(n)
  real ( kind = 8 ) x(n)

  value(1:n) = x(1:n) * ( 10.0D+00 * x(1:n) - 1.0D+00 ) &
    * ( 5.0D+00 * x(1:n) - 2.0D+00 ) * ( 5.0D+00 * x(1:n) - 2.0D+00 ) &
    * ( 4.0D+00 * x(1:n) - 3.4D+00 ) * ( x(1:n) - 1.0D+00 )

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns the title of problem p03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
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

  title = 'f(x) = x (10*x-1) (5x-2) (5x-2) (4x-3.4) (x-1)'

  return
end
subroutine p04_f ( n, x, value )

!*****************************************************************************80
!
!! P04_F evaluates the function for problem p04.
!
!  Discussion:
!
!    Step function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) value(n)
  real ( kind = 8 ) x(n)

  value(1:n) = atan ( 40.0D+00 * x(1:n) - 15.0D+00 )

  return
end
subroutine p04_title ( title )

!*****************************************************************************80
!
!! P04_TITLE returns the title of problem p04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
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

  title = 'f(x) = atan ( 40 * x - 15 )'

  return
end
subroutine p05_f ( n, x, value )

!*****************************************************************************80
!
!! P05_F evaluates the function for problem p05.
!
!  Discussion:
!
!    Step function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) value(n)
  real ( kind = 8 ) x(n)

  value(1:n) =           cos (  7.0D+00 * x(1:n) ) &
             + 5.0D+00 * cos ( 11.2D+00 * x(1:n) ) &
             - 2.0D+00 * cos ( 14.0D+00 * x(1:n) ) &
             + 5.0D+00 * cos ( 31.5D+00 * x(1:n) ) &
             + 7.0D+00 * cos ( 63.0D+00 * x(1:n) )

  return
end
subroutine p05_title ( title )

!*****************************************************************************80
!
!! P05_TITLE returns the title of problem p05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
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

  title = 'f(x) = cos(7*x)+5*cos(11.2*x)-2*cos(14*x)+5*cos(31.5*x)+7*cos(63*x).'

  return
end
subroutine p06_f ( n, x, value )

!*****************************************************************************80
!
!! P06_F evaluates the function for problem p06.
!
!  Discussion:
!
!    f(x) = exp ( - (4 * x - 1)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz,
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D Reidel, 1987, pages 337-340,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) value(n)
  real ( kind = 8 ) x(n)

  value(1:n) = exp ( - ( 4.0D+00 * x(1:n) - 1.0D+00 ) ** 2 )

  return
end
subroutine p06_title ( title )

!*****************************************************************************80
!
!! P06_TITLE returns the title of problem p06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
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

  title = 'f(x) = exp ( - ( 4*x-1 )^2 )'

  return
end
subroutine p07_f ( n, x, value )

!*****************************************************************************80
!
!! P07_F evaluates the function for problem p07.
!
!  Discussion:
!
!    f(x) = exp ( 4 * x ) if x <= 1/2
!           0                  otherwise
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz,
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D Reidel, 1987, pages 337-340,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) value(n)
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( x(i) < 0.5D+00 ) then
      value(i) = exp ( 4.0D+00 * x(i) )
    else
      value(i) = 0.0D+00
    end if
  end do

  return
end
subroutine p07_title ( title )

!*****************************************************************************80
!
!! P07_TITLE returns the title of problem p07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
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

  title = 'f(x) = exp ( 2 x ) if x < 0.5, 0 otherwise'

  return
end
subroutine p08_f ( n, x, value )

!*****************************************************************************80
!
!! P08_F evaluates the function for problem p08.
!
!  Discussion:
!
!    This is a famous example, due to Runge.  If equally spaced
!    abscissas are used, the sequence of interpolating polynomials Pn(X)
!    diverges, in the sense that the max norm of the difference
!    between Pn(X) and F(X) becomes arbitrarily large as N increases.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) value(n)
  real ( kind = 8 ) x(n)

  value(1:n) = 1.0D+00 / ( ( 10.0D+00 * ( x(1:n) - 0.5D+00 ) ) ** 2 + 1.0D+00 )

  return
end
subroutine p08_title ( title )

!*****************************************************************************80
!
!! P08_TITLE returns the title of problem p08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2012
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

  title = 'f(x) = 1 / ( 1 + ( 10 * (x-1/2) )^2 )'

  return
end
