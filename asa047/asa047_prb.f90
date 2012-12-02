program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA047_PRB.
!
!  Discussion:
!
!    ASA047_PRB calls the ASA047 routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA047_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA047 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA047_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of NELMIN on ROSENBROCK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) kcount
  integer ( kind = 4 ) konvge
  integer ( kind = 4 ) numres
  real ( kind = 8 ) reqmin
  real ( kind = 8 ), external :: rosenbrock
  real ( kind = 8 ) start(n)
  real ( kind = 8 ) step(n)
  real ( kind = 8 ) xmin(n)
  real ( kind = 8 ) ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Apply NELMIN to the ROSENBROCK function.'

  start(1:n) = (/ -1.2D+00, 1.0D+00 /)

  reqmin = 1.0D-08

  step(1:n) = (/ 1.0D+00, 1.0D+00 /)

  konvge = 10
  kcount = 500

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting point X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) start(i)
  end do

  ynewlo = rosenbrock ( start )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

  call nelmin ( rosenbrock, n, start, xmin, ynewlo, reqmin, step, &
    konvge, kcount, icount, numres, ifault )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate of minimizing value X*:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) xmin(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations = ', icount
  write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

  return
end
function rosenbrock ( x )

!*****************************************************************************80
!
!! ROSENBROCK evaluates the Rosenbrock parabolic value function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(2), the argument.
!
!    Output, real ( kind = 8 ) ROSENBROCK, the value of the function.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  real ( kind = 8 ) rosenbrock
  real ( kind = 8 ) x(3)

  fx1 = x(2) - x(1) * x(1)
  fx2 = 1.0D+00 - x(1)

  fx = 100.0D+00 * fx1 * fx1 &
     +             fx2 * fx2

  rosenbrock = fx

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates the use of NELMIN on POWELL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) kcount
  integer ( kind = 4 ) konvge
  integer ( kind = 4 ) numres
  real ( kind = 8 ), external :: powell
  real ( kind = 8 ) reqmin
  real ( kind = 8 ) start(n)
  real ( kind = 8 ) step(n)
  real ( kind = 8 ) xmin(n)
  real ( kind = 8 ) ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Apply NELMIN to POWELL quartic function.'

  start(1:n) = (/  3.0D+00, - 1.0D+00,   0.0D+00,   1.0D+00 /)

  reqmin = 1.0D-08

  step(1:n) = (/ 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /)

  konvge = 10
  kcount = 500

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting point X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) start(i)
  end do

  ynewlo = powell ( start )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

  call nelmin ( powell, n, start, xmin, ynewlo, reqmin, step, &
    konvge, kcount, icount, numres, ifault )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate of minimizing value X*:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) xmin(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations = ', icount
  write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

  return
end
function powell ( x )

!*****************************************************************************80
!
!! POWELL evaluates the Powell quartic function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(4), the argument.
!
!    Output, real ( kind = 8 ) POWELL, the value of the function.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  real ( kind = 8 ) fx3
  real ( kind = 8 ) fx4
  real ( kind = 8 ) powell
  real ( kind = 8 ) x(4)

  fx1 = x(1) + 10.0D+00 * x(2)
  fx2 = x(3) - x(4)
  fx3 = x(2) - 2.0D+00 * x(3)
  fx4 = x(1) - x(4)

  fx =            fx1 * fx1 &
     +  5.0D+00 * fx2 * fx2 &
     +            fx3 * fx3 * fx3 * fx3 &
     + 10.0D+00 * fx4 * fx4 * fx4 * fx4

  powell = fx

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 demonstrates the use of NELMIN on HELICAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), external :: helical
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) kcount
  integer ( kind = 4 ) konvge
  integer ( kind = 4 ) numres
  real ( kind = 8 ) reqmin
  real ( kind = 8 ) start(n)
  real ( kind = 8 ) step(n)
  real ( kind = 8 ) xmin(n)
  real ( kind = 8 ) ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Apply NELMIN to the HELICAL function.'

  start(1:n) = (/ - 1.0D+00,   0.0D+00,   0.0D+00 /)

  reqmin = 1.0D-08

  step(1:n) = (/ 1.0D+00, 1.0D+00, 1.0D+00 /)

  konvge = 10
  kcount = 500

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting point X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) start(i)
  end do

  ynewlo = helical ( start )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

  call nelmin ( helical, n, start, xmin, ynewlo, reqmin, step, &
    konvge, kcount, icount, numres, ifault )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate of minimizing value X*:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) xmin(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations = ', icount
  write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

  return
end
function helical ( x )

!*****************************************************************************80
!
!! HELICAL evaluates the Fletcher-Powell helical valley function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(3), the argument.
!
!    Output, real ( kind = 8 ) HELICAL, the value of the function.
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  real ( kind = 8 ) fx3
  real ( kind = 8 ) helical
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(3)

  if ( 0.0D+00 < x(1) ) then
    theta = atan2 ( x(2), x(1) ) / 2.0D+00 / pi
  else if ( x(1) < 0.0D+00 ) then
    theta = 0.5D+00 + atan2 ( x(2), x(1) ) / 2.0D+00 / pi
  else if ( x(1) == 0.0D+00 ) then
    theta = 0.25D+00
  end if

  fx1 = x(3) - 10.0D+00 * theta
  fx2 = sqrt ( x(1) * x(1) + x(2) * x(2) )
  fx3 = x(3)

  fx = 100.0D+00 * fx1 * fx1 &
     +             fx2 * fx2 &
     +             fx3 * fx3

  helical = fx

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 demonstrates the use of NELMIN on QUARTIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) kcount
  integer ( kind = 4 ) konvge
  integer ( kind = 4 ) numres
  real ( kind = 8 ), external :: quartic
  real ( kind = 8 ) reqmin
  real ( kind = 8 ) start(n)
  real ( kind = 8 ) step(n)
  real ( kind = 8 ) xmin(n)
  real ( kind = 8 ) ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Apply NELMIN to the QUARTIC function.'

  start(1:n) = 1.0D+00

  reqmin = 1.0D-08

  step(1:n) = 1.0D+00

  konvge = 10
  kcount = 2000

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting point X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) start(i)
  end do

  ynewlo = quartic ( start )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

  call nelmin ( quartic, n, start, xmin, ynewlo, reqmin, step, &
    konvge, kcount, icount, numres, ifault )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate of minimizing value X*:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) xmin(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations = ', icount
  write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

  return
end
function quartic ( x )

!*****************************************************************************80
!
!! QUARTIC evaluates a function defined by a sum of fourth powers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(10), the argument.
!
!    Output, real ( kind = 8 ) QUARTIC, the value of the function.
!
  implicit none

  real ( kind = 8 ) quartic
  real ( kind = 8 ) x(10)

  quartic = sum ( x(1:10)**4 )

  return
end
