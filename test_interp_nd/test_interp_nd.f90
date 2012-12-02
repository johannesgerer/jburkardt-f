function csevl ( x, cs, n )

!*****************************************************************************80
!
!! CSEVL evaluates an N term Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, 1973, pages 254-256.
!
!    Leslie Fox, Ian Parker,
!    Chebyshev Polynomials in Numerical Analysis,
!    Oxford Press, 1968,
!    LC: QA297.F65.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument at which the series is 
!    to be evaluated.
!    X must satisfy -1.0 <= X <= 1.0.
!
!    Input, real ( kind = 8 ) CS(N), the array of N terms of a Chebyshev series.
!    In evaluating CS, only half the first coefficient is summed.
!
!    Input, integer ( kind = 4 ) N, the number of terms in array CS.
!    N must be at least 1, and no more than 1000.
!
!    Output, real ( kind = 8 ) CSEVL, the value of the Chebyshev series.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) cs(n)
  real ( kind = 8 ) csevl
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSEVL - Fatal error!'
    write ( *, '(a)' ) '  Number of terms N is less than 1.'
    stop
  end if

  if ( 1000 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSEVL - Fatal error!'
    write ( *, '(a)' ) '  The number of terms is more than 1000.'
    stop
  end if

  if ( x < -1.0D+00 .or. 1.0D+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSEVL - Fatal error!'
    write ( *, '(a)' ) '  The input argument X is outside the interval [-1,1].'
    stop
  end if

  b1 = 0.0D+00
  b0 = 0.0D+00

  do i = n, 1, -1
    b2 = b1
    b1 = b0
    b0 = 2.0D+00 * x * b1 - b2 + cs(i)
  end do

  csevl = 0.5D+00 * ( b0 - b2 )

  return
end
function inits ( os, nos, eta )

!*****************************************************************************80
!
!! INITS estimates the order of an orthogonal series for a given accuracy.
!
!  Discussion:
!
!    Because this is a function, it is not possible to print out
!    warning error messages.  Therefore, if an error condition is
!    detected, a bogus value of INITS is returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, 1973, pages 254-256.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) OS(NOS), the coefficients in the series.
!
!    Input, integer ( kind = 4 ) NOS, the number of coefficients.  NOS must be
!    at least 1, or an error condition arises.
!
!    Input, real ( kind = 8 ) ETA, the requested accuracy of the series.
!    Ordinarily, ETA will be chosen to be one-tenth machine precision.
!
!    Output, integer ( kind = 4 ) INITS, the order of the series guaranteeing 
!    the given accuracy.  However, on error, INITS will be returned
!    as a negative number.
!
  implicit none

  integer ( kind = 4 ) nos

  real ( kind = 8 ) err
  real ( kind = 8 ) eta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inits
  real ( kind = 8 ) os(nos)
!
!  Fatal error.  Number of coefficients less than 1.
!  But an error message cannot be printed out because this is a function.
!
  if ( nos < 1 ) then
    inits = - huge ( i )
    return
  end if

  err = 0.0D+00

  i = 0

  do ii = 1, nos

    i = nos + 1 - ii
    err = err + abs ( os(i) )

    if ( eta < err ) then
      i = nos + 1 - ii
      exit
    end if

  end do
!
!  Eta may be too small.  Accuracy cannot be guaranteed.
!  But an error message cannot be printed out because this is a function.
!
  if ( i == 0 ) then
    i = - nos
  end if

  inits = i

  return
end
subroutine p00_c ( prob, m, seed, c )

!*****************************************************************************80
!
!! P00_CW computes a random C parameter vector for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem number.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) C(M), the parameter vector.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ), dimension ( 6 ) :: b = (/ &
    1.5D+00, 0.0D+00, 1.85D+00, 7.03D+00, 20.4D+00, 4.3D+00 /)
  real ( kind = 8 ) c(m)
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) seed

  b(2) = real ( m, kind = 8 )

  call r8vec_uniform_01 ( m, seed, c )
  c(1:m) = b(prob) * c(1:m) / sum ( c )

  return
end
subroutine p00_d ( prob, m, id, c, w, n, x, d )

!*****************************************************************************80
!
!! P00_D returns a derivative component of any function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the index of the function.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) ID, the spatial coordinate to differentiate.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) D(N), the ID-th derivative component.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) id
  integer ( kind = 4 ) prob
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)  

  if ( id < 0 .or. m < id ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_D - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal spatial coordinate ID = ', id
    stop
  end if

  if ( prob == 1 ) then
    call p01_d ( m, id, c, w, n, x, d )
  else if ( prob == 2 ) then
    call p02_d ( m, id, c, w, n, x, d )
  else if ( prob == 3 ) then
    call p03_d ( m, id, c, w, n, x, d )
  else if ( prob == 4 ) then
    call p04_d ( m, id, c, w, n, x, d )
  else if ( prob == 5 ) then
    call p05_d ( m, id, c, w, n, x, d )
  else if ( prob == 6 ) then
    call p06_d ( m, id, c, w, n, x, d )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_D - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal function index PROB = ', prob
    stop
  end if

  return
end
subroutine p00_f ( prob, m, c, w, n, x, f )

!*****************************************************************************80
!
!! P00_F returns the value of any function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the index of the function.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) prob
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)  

  if ( prob == 1 ) then
    call p01_f ( m, c, w, n, x, f )
  else if ( prob == 2 ) then
    call p02_f ( m, c, w, n, x, f )
  else if ( prob == 3 ) then
    call p03_f ( m, c, w, n, x, f )
  else if ( prob == 4 ) then
    call p04_f ( m, c, w, n, x, f )
  else if ( prob == 5 ) then
    call p05_f ( m, c, w, n, x, f )
  else if ( prob == 6 ) then
    call p06_f ( m, c, w, n, x, f )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_F - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal function index PROB = ', prob
    stop
  end if

  return
end
subroutine p00_prob_num ( prob_num )

!*****************************************************************************80
!
!! P00_PROB_NUM returns the number of test functions available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!   Output, integer ( kind = 4 ) PROB_NUM, the number of test functions.
!
  implicit none

  integer ( kind = 4 ) prob_num

  prob_num = 6

  return
end
subroutine p00_q ( prob, m, c, w, q )

!*****************************************************************************80
!
!! P00_Q returns the integral of any function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the index of the function.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Output, real ( kind = 8 ) Q, the integral.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) c(m)
  integer ( kind = 4 ) prob
  real ( kind = 8 ) q
  real ( kind = 8 ) w(m)

  if ( prob == 1 ) then
    call p01_q ( m, c, w, q )
  else if ( prob == 2 ) then
    call p02_q ( m, c, w, q )
  else if ( prob == 3 ) then
    call p03_q ( m, c, w, q )
  else if ( prob == 4 ) then
    call p04_q ( m, c, w, q )
  else if ( prob == 5 ) then
    call p05_q ( m, c, w, q )
  else if ( prob == 6 ) then
    call p06_q ( m, c, w, q )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_Q - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal function index PROB = ', prob
    stop
  end if

  return
end
subroutine p00_title ( prob, title )

!*****************************************************************************80
!
!! P00_TITLE returns the title for any function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the index of the function.
!
!    Output, character ( len = * ) TITLE, the function title.
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
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal function index PROB = ', prob
    stop
  end if

  return
end
subroutine p00_w ( prob, m, seed, w )

!*****************************************************************************80
!
!! P00_W computes a random W parameter vector for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem number.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) W(M), the parameter vector.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) prob
  integer ( kind = 4 ) seed
  real ( kind = 8 ) w(m)

  call r8vec_uniform_01 ( m, seed, w )

  return
end
subroutine p01_d ( m, id, c, w, n, x, d )

!*****************************************************************************80
!
!! P01_D evaluates any derivative component for problem p01.
!
!  Discussion:
!
!    f(x) = cos ( 2 * pi * w(1) + sum ( c(1:m) * x(1:m) ) )
!
!    Default values are:
!
!    c(1:m) = 1/m
!    w(1) = 0.3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
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
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) ID, the spatial coordinate to differentiate.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) D(N), the ID-th derivative component.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  d(1:n) = 2.0D+00 * pi * w(1)
  do i = 1, m
    d(1:n) = d(1:n) + c(i) * x(i,1:n)
  end do
  d(1:n) = - c(id) * sin ( d(1:n) )

  return
end
subroutine p01_f ( m, c, w, n, x, f )

!*****************************************************************************80
!
!! P01_F evaluates the function for problem p01.
!
!  Discussion:
!
!    f(x) = cos ( 2 * pi * w(1) + sum ( c(1:m) * x(1:m) ) )
!
!    Default values are:
!
!    c(1:m) = 1/m
!    w(1) = 0.3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
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
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  f(1:n) = 2.0D+00 * pi * w(1)
  do i = 1, m
    f(1:n) = f(1:n) + c(i) * x(i,1:n)
  end do
  f(1:n) = cos ( f(1:n) )

  return
end
subroutine p01_q ( m, c, w, q )

!*****************************************************************************80
!
!! P01_Q evaluates the integral for problem p01.
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
!  Reference:
!
!    Alan Genz,
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Output, real ( kind = 8 ) Q, the integral.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) c(m)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) q
  real ( kind = 8 ) w(m)

  q = 2.0D+00 ** m * &
    cos ( 2.0D+00 * pi * w(1) + 0.5D+00 * sum ( c(1:m) ) ) * &
    product ( sin ( 0.5D+00 * c(1:m) ) / c(1:m) )

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns the name of problem p01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 2012
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

  title = 'Oscillatory'

  return
end
subroutine p02_d ( m, id, c, w, n, x, d )

!*****************************************************************************80
!
!! P02_D evaluates an derivative component for problem p02.
!
!  Discussion:
!
!    f(x) = 1 / product ( c(1:m)^(-2) + ( x(1:m) - w(1:m) )^2 )
!
!    Default values are:
!
!    c(1:m) = 1
!    w(1:m) = 0.5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
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
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) ID, the spatial coordinate to differentiate.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) D(N), the ID-th derivative component.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  d(1:n) = 1.0D+00
  do i = 1, m
    d(1:n) = d(1:n) * ( c(i)**(-2) + ( x(i,1:n) - w(i) )**2 )
  end do
  d(1:n) = - 2.0D+00 / d(1:n) * ( x(id,1:n) - w(id) ) / &
    ( c(id) ** ( -2 ) + ( x(id,1:n) - w(id) ) ** 2 )

  return
end
subroutine p02_f ( m, c, w, n, x, f )

!*****************************************************************************80
!
!! P02_F evaluates the function for problem p02.
!
!  Discussion:
!
!    f(x) = 1 / product ( c(1:m)^(-2) + ( x(1:m) - w(1:m) )^2 )
!
!    Default values are:
!
!    c(1:m) = 1
!    w(1:m) = 0.5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
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
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  f(1:n) = 1.0D+00
  do i = 1, m
    f(1:n) = f(1:n) * ( c(i)**(-2) + ( x(i,1:n) - w(i) )**2 )
  end do
  f(1:n) = 1.0D+00 / f(1:n)

  return
end
subroutine p02_q ( m, c, w, q )

!*****************************************************************************80
!
!! P02_Q evaluates the integral for problem p02.
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
!  Reference:
!
!    Alan Genz,
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Output, real ( kind = 8 ) Q, the integral.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) q
  real ( kind = 8 ) w(m)

  q = product ( &
               (   atan ( ( 1.0D+00 - w(1:m) ) * c(1:m) ) &
                 + atan (             w(1:m)   * c(1:m) ) &
               ) * c(1:m) &
             )
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
!    04 August 2012
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

  title = 'Product Peak'

  return
end
subroutine p03_d ( m, id, c, w, n, x, d )

!*****************************************************************************80
!
!! P03_D evaluates any derivative component for problem p03.
!
!  Discussion:
!
!    f(x) = 1 / ( 1 + sum ( c(1:m) * x(1:m) ) ) ^ ( m + 1 )
!
!    Default values are:
!
!    c(1:m) = 1/m
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
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
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) ID, the spatial coordinate to differentiate.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) D(N), the ID-th derivative component.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  d(1:n) = 1.0D+00
  do i = 1, m
    d(1:n) = d(1:n) + c(i) * x(i,1:n)
  end do
  d(1:n) = - c(id) * real ( m + 1, kind = 8 ) / d(1:n) ** ( m + 2 )

  return
end
subroutine p03_f ( m, c, w, n, x, f )

!*****************************************************************************80
!
!! P03_F evaluates the function for problem p03.
!
!  Discussion:
!
!    f(x) = 1 / ( 1 + sum ( c(1:m) * x(1:m) ) ) ^ ( m + 1 )
!
!    Default values are:
!
!    c(1:m) = 1/m
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
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
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  f(1:n) = 1.0D+00
  do i = 1, m
    f(1:n) = f(1:n) + c(i) * x(i,1:n)
  end do
  f(1:n) = 1.0D+00 / f(1:n) ** ( m + 1 )

  return
end
subroutine p03_q ( m, c, w, q )

!*****************************************************************************80
!
!! P03_Q evaluates the integral for problem p03.
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
!  Reference:
!
!    Alan Genz,
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Output, real ( kind = 8 ) Q, the integral.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) c(m)
  integer ( kind = 4 ) ivec(m)
  real ( kind = 8 ) q
  integer ( kind = 4 ) rank
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) r8vec_i4vec_dot_product
  integer ( kind = 4 ) s
  real ( kind = 8 ) w(m)
!
!  Here, we need to generate all possible DIM_NUM tuples with
!  values of 0 or 1.
! 
  q = 0.0D+00
  rank = 0

  do

    call tuple_next ( 0, 1, m, rank, ivec )

    if ( rank == 0 ) then
      exit
    end if

    s = sum ( ivec(1:m) )

    q = q + r8_mop ( s ) / ( 1.0D+00 + r8vec_i4vec_dot_product ( m, c, ivec ) )

  end do

  q = q / ( r8_factorial ( m ) * product ( c(1:m) ) )

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
!    04 August 2012
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

  title = 'Corner Peak'

  return
end
subroutine p04_d ( m, id, c, w, n, x, d )

!*****************************************************************************80
!
!! P04_D evaluates any derivative component for problem p04.
!
!  Discussion:
!
!    f(x) = exp ( - sum ( c(1:m)^2 * ( x(1:m) - w(1:m) )^2 )
!
!    Default values are:
!
!    c(1:m) = 1 / m
!    w(1:m) = 0.5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
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
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) ID, the spatial coordinate to differentiate.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) D(N), the ID-th derivative component.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  d(1:n) = 0.0D+00
  do i = 1, m
    d(1:n) = d(1:n) + ( c(i) * ( x(i,1:n) - w(i) ) )**2
  end do
  d(1:n) = exp ( - d(1:n) ) * c(id)**2 * ( - 2.0D+00 ) * ( x(id,1:n) - w(id) )

  return
end
subroutine p04_f ( m, c, w, n, x, f )

!*****************************************************************************80
!
!! P04_F evaluates the function for problem p04.
!
!  Discussion:
!
!    f(x) = exp ( - sum ( c(1:m)^2 * ( x(1:m) - w(1:m) )^2 )
!
!    Default values are:
!
!    c(1:m) = 1 / m
!    w(1:m) = 0.5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
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
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  f(1:n) = 0.0D+00
  do i = 1, m
    f(1:n) = f(1:n) + ( c(i) * ( x(i,1:n) - w(i) ) )**2
  end do
  f(1:n) = exp ( - f(1:n) )

  return
end
subroutine p04_q ( m, c, w, q )

!*****************************************************************************80
!
!! P04_Q evaluates the integral for problem p04.
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
!  Reference:
!
!    Alan Genz,
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Output, real ( kind = 8 ) Q, the integral.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) c(m)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) q
  real ( kind = 8 ) r8_error
  real ( kind = 8 ) w(m)

  q = 1.0D+00

  do i = 1, m

    q = q * sqrt ( pi ) &
      * ( r8_error ( c(i) * ( 1.0D+00 - w(i) ) ) &
        + r8_error ( c(i) *             w(i) ) ) &
      / ( 2.0D+00 * c(i) )

  end do

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
!    04 August 2012
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

  title = 'Gaussian'

  return
end
subroutine p05_d ( m, id, c, w, n, x, d )

!*****************************************************************************80
!
!! P05_D evaluates any derivative component for problem p05.
!
!  Discussion:
!
!    f(x) = exp ( - sum ( c(1:m) * abs ( x(1:m) - w(1:m) ) ) )
!
!    Default values are:
!
!    c(1:m) = 2.0
!    w(1:m) = 0.5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
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
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) ID, the spatial coordinate to differentiate.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) D(N), the ID-th derivative component.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) j
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  d(1:n) = 0.0D+00
  do i = 1, m
    d(1:n) = d(1:n) + c(i) * abs ( x(i,1:n) - w(i) )
  end do
  d(1:n) = exp ( - d(1:n) )

  do j = 1, n
    if ( x(id,j) - w(id) <= 0.0D+00 ) then
      d(j) = d(j) * c(id)
    else
      d(j) = - d(j) * c(id)
    end if
  end do

  return
end
subroutine p05_f ( m, c, w, n, x, f )

!*****************************************************************************80
!
!! P05_F evaluates the function for problem p05.
!
!  Discussion:
!
!    f(x) = exp ( - sum ( c(1:m) * abs ( x(1:m) - w(1:m) ) ) )
!
!    Default values are:
!
!    c(1:m) = 2.0
!    w(1:m) = 0.5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
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
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  f(1:n) = 0.0D+00
  do i = 1, m
    f(1:n) = f(1:n) + c(i) * abs ( x(i,1:n) - w(i) )
  end do
  f(1:n) = exp ( - f(1:n) )

  return
end
subroutine p05_q ( m, c, w, q )

!*****************************************************************************80
!
!! P05_Q evaluates the integral for problem p05.
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
!  Reference:
!
!    Alan Genz,
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Output, real ( kind = 8 ) Q, the integral.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) c(m)
  integer ( kind = 4 ) i
  real ( kind = 8 ) q
  real ( kind = 8 ) w(m)

  q = 1.0D+00

  do i = 1, m
!
!  W < 0 < 1
!
!  | X - W | = X - W from 0 to 1.
!
    if ( w(i) < 0.0D+00 ) then

      q = q * &
      ( exp ( - c(i) * (         - w(i) ) ) &
      - exp ( - c(i) * ( 1.0D+00 - w(i) ) ) ) / c(i)
!
!  0 < W < 1
!
!  | X - W | = W - X from 0 to Z, 
!            = X - W from      Z to 1.
!
    else if ( w(i) < 1.0D+00 ) then

      q = q * ( 2.0D+00 &
          - exp ( - c(i) * (           w(i) ) ) &
          - exp ( - c(i) * ( 1.0D+00 - w(i) ) ) ) / c(i)
!
!  0 < 1 < W
!
!  | X - W | = W - X from 0 to 1.
!
    else

      q = q * &
      ( exp ( - c(i) * ( w(i) - 1.0D+00 ) ) &
      - exp ( - c(i) * ( w(i)           ) ) ) / c(i)

    end if

  end do

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
!    04 August 2012
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

  title = 'Continuous'

  return
end
subroutine p06_d ( m, id, c, w, n, x, d )

!*****************************************************************************80
!
!! P06_D evaluates any derivative component for problem p06.
!
!  Discussion:
!
!    f(x) = exp ( c(1:m) * x(1:m) ) if x(1) <= w(1) and x(2) <= w(2).
!           0                          otherwise
!
!    Default values are:
!
!    c(1:m) = 0.5^(1/m)
!    w(1:2) = 0.5^(1/m)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
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
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) ID, the spatial coordinate to differentiate.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) D(N), the ID-th derivative component.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) j
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  if ( m == 1 ) then

    d(1:n) = c(1) * exp ( c(1) * x(1,1:n) )

    do j = 1, n
      if ( w(1) < x(1,j) ) then
        d(j) = 0.0D+00
      end if
    end do

  else

    d(1:n) = 0.0D+00
    do i = 1, m
      d(1:n) = d(1:n) + c(i) * x(i,1:n)
    end do
    d(1:n) = c(id) * exp ( d(1:n) )

    do j = 1, n
      if ( w(1) < x(1,j) .or. w(2) < x(2,j) ) then
        d(j) = 0.0D+00
      end if
    end do

  end if

  return
end
subroutine p06_f ( m, c, w, n, x, f )

!*****************************************************************************80
!
!! P06_F evaluates the function for problem p06.
!
!  Discussion:
!
!    f(x) = exp ( c(1:m) * x(1:m) ) if x(1) <= w(1) and x(2) <= w(2).
!           0                          otherwise
!
!    Default values are:
!
!    c(1:m) = 0.5^(1/m)
!    w(1:2) = 0.5^(1/m)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2012
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
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) w(m)
  real ( kind = 8 ) x(m,n)

  if ( m == 1 ) then

    f(1:n) = exp ( c(1) * x(1,1:n) )

    do j = 1, n
      if ( w(1) < x(1,j) ) then
        f(j) = 0.0D+00
      end if
    end do

  else

    f(1:n) = 0.0D+00
    do i = 1, m
      f(1:n) = f(1:n) + c(i) * x(i,1:n)
    end do
    f(1:n) = exp ( f(1:n) )

    do j = 1, n
      if ( w(1) < x(1,j) .or. w(2) < x(2,j) ) then
        f(j) = 0.0D+00
      end if
    end do

  end if

  return
end
subroutine p06_q ( m, c, w, q )

!*****************************************************************************80
!
!! P06_Q evaluates the integral for problem p06.
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
!  Reference:
!
!    Alan Genz,
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the argument.
!
!    Input, real ( kind = 8 ) C(M), W(M), the problem parameters.
!
!    Output, real ( kind = 8 ) Q, the integral.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) c(m)
  integer ( kind = 4 ) i
  real ( kind = 8 ) q
  real ( kind = 8 ) w(m)
!
!  To simplify the calculation, force W(3:M) to be at least 1.0.
!
  w(3:m) = 1.0D+00

  q = 1.0D+00

  do i = 1, m

    if ( w(i) <= 0.0D+00 ) then

      q = q * 0.0D+00

    else if ( w(i) <= 1.0D+00 ) then

      if ( c(i) == 0.0D+00 ) then
        q = q * w(i)
      else
        q = q * ( exp ( c(i) * w(i) ) - 1.0D+00 ) / c(i)
      end if

    else

      if ( c(i) /= 0.0D+00 ) then
        q = q * ( exp ( c(i) * w(i) ) - 1.0D+00 ) / c(i)
      end if

    end if

  end do

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
!    04 August 2012
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

  title = 'Discontinuous'

  return
end
function r8_error ( x )

!*****************************************************************************80
!
!! R8_ERROR computes the error function.
!
!  Discussion:
!
!    This function was renamed "R8_ERROR" from "ERF", to avoid a conflict
!    with the name of a corresponding routine often, but not always,
!    supplied as part of the math support library.
!
!    The definition of the error function is:
!
!      erf(x) = ( 2 / sqrt ( pi ) ) * Integral ( 0 <= t <= x ) exp ( -t^2 ) dt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the error function.
!
!    Output, real ( kind = 8 ) R8_ERROR, the value of the error function at X.
!
  implicit none

  real ( kind = 8 ) csevl
  real ( kind = 8 ) r8_error
  real ( kind = 8 ) r8_errorc
  real ( kind = 8 ), parameter, dimension ( 13 ) :: erfcs = (/ &
    -0.049046121234691808D+00, -0.14226120510371364D+00, &
     0.010035582187599796D+00, -0.000576876469976748D+00, &
     0.000027419931252196D+00, -0.000001104317550734D+00, &
     0.000000038488755420D+00, -0.000000001180858253D+00, &
     0.000000000032334215D+00, -0.000000000000799101D+00, &
     0.000000000000017990D+00, -0.000000000000000371D+00, &
     0.000000000000000007D+00 /)
  integer ( kind = 4 ) inits
  integer ( kind = 4 ), save :: nterf = 0
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ), parameter :: sqrtpi = 1.7724538509055160D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xbig = 0.0D+00
  real ( kind = 8 ) y
!
!  Initialize the Chebyshev series.
!
  if ( nterf == 0 ) then
    nterf = inits ( erfcs, 13, 0.1D+00 * epsilon ( erfcs ) )
    xbig = sqrt ( - log ( sqrtpi * epsilon ( xbig ) ) )
    sqeps = sqrt ( 2.0D+00 * epsilon ( sqeps ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then
    r8_error = 2.0D+00 * x / sqrtpi
  else if ( y <= 1.0D+00 ) then
    r8_error = x &
      * ( 1.0D+00 + csevl ( 2.0D+00 * x**2 - 1.0D+00, erfcs, nterf ) )
  else if ( y <= xbig ) then
    r8_error = sign ( 1.0D+00 - r8_errorc ( y ), x )
  else
    r8_error = sign ( 1.0D+00, x )
  end if

  return
end
function r8_errorc ( x )

!*****************************************************************************80
!
!! R8_ERRORC computes the complementary error function.
!
!  Discussion:
!
!    This function was renamed "R8_ERRORC" from "ERFC", to avoid a conflict
!    with the name of a corresponding routine often, but not always,
!    supplied as part of the math support library.
!
!    The definition of the complementary error function is:
!
!      ERFC(X) = 1 - ERF(X)
!
!    where ERF(X) is the error function.
!
!  Modified:
!
!    26 August 2001
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the error function.
!
!    Output, real ( kind = 8 ) R8_ERRORC, the value of the complementary 
!    error function at X.
!
  implicit none

  real ( kind = 8 ) csevl
  real ( kind = 8 ), parameter, dimension ( 13 ) :: erfcs = (/ &
   -0.049046121234691808D+00, &
   -0.14226120510371364D+00,  &
    0.010035582187599796D+00, &
   -0.000576876469976748D+00, &
    0.000027419931252196D+00, &
   -0.000001104317550734D+00, &
    0.000000038488755420D+00, &
   -0.000000001180858253D+00, &
    0.000000000032334215D+00, &
   -0.000000000000799101D+00, &
    0.000000000000017990D+00, &
   -0.000000000000000371D+00, &
    0.000000000000000007D+00 /)
  real ( kind = 8 ), parameter, dimension ( 24 ) :: erfccs = (/ &
    0.0715179310202925D+00, &
   -0.026532434337606719D+00, &
    0.001711153977920853D+00, &
   -0.000163751663458512D+00, &
    0.000019871293500549D+00, &
   -0.000002843712412769D+00, &
    0.000000460616130901D+00, &
   -0.000000082277530261D+00, &
    0.000000015921418724D+00, &
   -0.000000003295071356D+00, &
    0.000000000722343973D+00, &
   -0.000000000166485584D+00, &
    0.000000000040103931D+00, &
   -0.000000000010048164D+00, &
    0.000000000002608272D+00, &
   -0.000000000000699105D+00, &
    0.000000000000192946D+00, &
   -0.000000000000054704D+00, &
    0.000000000000015901D+00, &
   -0.000000000000004729D+00, &
    0.000000000000001432D+00, &
   -0.000000000000000439D+00, &
    0.000000000000000138D+00, &
   -0.000000000000000048D+00 /)
  real ( kind = 8 ), parameter, dimension ( 23 ) :: erc2cs = (/ &
   -0.069601346602309501D+00, &
   -0.041101339362620893D+00, &
    0.003914495866689626D+00, &
   -0.000490639565054897D+00, &
    0.000071574790013770D+00, &
   -0.000011530716341312D+00, &
    0.000001994670590201D+00, &
   -0.000000364266647159D+00, &
    0.000000069443726100D+00, &
   -0.000000013712209021D+00, &
    0.000000002788389661D+00, &
   -0.000000000581416472D+00, &
    0.000000000123892049D+00, &
   -0.000000000026906391D+00, &
    0.000000000005942614D+00, &
   -0.000000000001332386D+00, &
    0.000000000000302804D+00, &
   -0.000000000000069666D+00, &
    0.000000000000016208D+00, &
   -0.000000000000003809D+00, &
    0.000000000000000904D+00, &
   -0.000000000000000216D+00, &
    0.000000000000000052D+00 /)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) inits
  integer ( kind = 4 ), save :: nterc2 = 0
  integer ( kind = 4 ), save :: nterf = 0
  integer ( kind = 4 ), save :: nterfc = 0
  real ( kind = 8 ) r8_errorc
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ), parameter :: sqrtpi = 1.7724538509055160D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xmax = 0.0D+00
  real ( kind = 8 ), save :: xsml = 0.0D+00
  real ( kind = 8 ) y

  if ( nterf == 0 ) then

    eta = 0.1D+00 * epsilon ( eta )
    nterf = inits ( erfcs, 13, eta )
    nterfc = inits ( erfccs, 24, eta )
    nterc2 = inits ( erc2cs, 23, eta )

    xsml = -sqrt ( - log ( sqrtpi * epsilon ( xsml ) ) )
    xmax = sqrt ( - log ( sqrtpi * tiny ( xmax ) ) )
    xmax = xmax - 0.5D+00 * log ( xmax ) / xmax - 0.01D+00
    sqeps = sqrt ( 2.0D+00 * epsilon ( sqeps ) )

  end if

  if ( x <= xsml ) then
    r8_errorc = 2.0D+00
    return
  end if
!
!  X so big that ERFC will underflow.
!
  if ( xmax < x ) then
    r8_errorc = 0.0D+00
    return
  end if

  y = abs ( x )
!
!  erfc(x) = 1.0D+00 - erf(x) for -1 <= x <= 1.
!
  if ( y <= 1.0D+00 ) then

    if ( y < sqeps ) then
      r8_errorc = 1.0D+00 - 2.0D+00 * x / sqrtpi
    else if ( sqeps <= y ) then
      r8_errorc = 1.0D+00 - x * ( 1.0D+00 + &
        csevl ( 2.0D+00 * x * x - 1.0D+00, erfcs, nterf ) )
    end if

    return

  end if
!
!  For 1 < |x| <= xmax, erfc(x) = 1.0D+00 - erf(x)
!
  y = y * y

  if ( y <= 4.0D+00 ) then
    r8_errorc = exp ( -y ) / abs ( x ) * ( 0.5D+00 &
      + csevl ( ( 8.0D+00 / y - 5.0D+00 ) / 3.0D+00, erc2cs, nterc2 ) )
  else
    r8_errorc = exp ( -y ) / abs ( x ) * ( 0.5D+00 &
      + csevl ( 8.0D+00 / y - 1.0D+00, erfccs, nterfc ) )
  end if

  if ( x < 0.0D+00 ) then
    r8_errorc = 2.0D+00 - r8_errorc
  end if

  return
end
subroutine tuple_next ( m1, m2, n, rank, x )

!*****************************************************************************80
!
!! TUPLE_NEXT computes the next element of a tuple space.
!
!  Discussion:
!
!    The elements are N vectors.  Each entry is constrained to lie
!    between M1 and M2.  The elements are produced one at a time.
!    The first element is
!      (M1,M1,...,M1),
!    the second element is
!      (M1,M1,...,M1+1),
!    and the last element is
!      (M2,M2,...,M2)
!    Intermediate elements are produced in lexicographic order.
!
!    For example, if our scalar input was
!
!      N = 2,
!      M1 = 1,
!      M2 = 3
!
!    then the sequence of input and output vectors would be:
!
!      INPUT        OUTPUT
!      -------      -------
!      Rank  X      Rank   X
!      ----  ---    -----  ---
!      0     * *    1      1 1
!      1     1 1    2      1 2
!      2     1 2    3      1 3
!      3     1 3    4      2 1
!      4     2 1    5      2 2
!      5     2 2    6      2 3
!      6     2 3    7      3 1
!      7     3 1    8      3 2
!      8     3 2    9      3 3
!      9     3 3    0      0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M1, M2, the minimum and maximum entries.
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input/output, integer ( kind = 4 ) RANK, counts the elements.
!    On first call, set RANK to 0.  Thereafter, the output value of RANK
!    will indicate the order of the element returned.  When there are no
!    more elements, RANK will be returned as 0.
!
!    Input/output, integer ( kind = 4 ) X(N), on input the previous tuple.
!    On output, the next tuple.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  if ( m2 < m1 ) then
    rank = 0
    return
  end if

  if ( rank <= 0 ) then

    x(1:n) = m1
    rank = 1

  else

    rank = rank + 1
    i = n

    do

      if ( x(i) < m2 ) then
        x(i) = x(i) + 1
        exit
      end if

      x(i) = m1

      if ( i == 1 ) then
        rank = 0
        x(1:n) = m1
        exit
      end if

      i = i - 1

    end do

  end if

  return
end
