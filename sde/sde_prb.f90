program main

!*****************************************************************************80
!
!! MAIN is the main program for SDE_PRB.
!
!  Discussion:
!
!    SDE_PRB demonstrates the use of SDE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SDE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SDE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )
  call test11 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SDE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BPATH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 500

  integer ( kind = 4 ) seed
  real ( kind = 8 ) w(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  BPATH generates a sample Brownian motion path'

  seed = 123456789

  call bpath ( seed, n, w )

  call bpath_gnuplot ( n, w )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests BPATH_AVERAGE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 1000
  integer ( kind = 4 ), parameter :: n = 500

  real ( kind = 8 ) error
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u(m,0:n)
  real ( kind = 8 ) umean(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  BPATH_AVERAGE generates many Brownian paths'
  write ( *, '(a)' ) '  and averages them.'

  seed = 123456789

  call bpath_average ( seed, m, n, u, umean, error )

  call bpath_average_gnuplot ( m, n, u, umean )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests CHAIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 200

  real ( kind = 8 ) diff
  integer ( kind = 4 ) seed
  real ( kind = 8 ) vem(0:n)
  real ( kind = 8 ) xem(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  CHAIN solves a stochastic differential equation for'
  write ( *, '(a)' ) '  a function of a stochastic variable X.'
  write ( *, '(a)' ) '  We can solve for X(t), and then evaluate V(X(t)).'
  write ( *, '(a)' ) '  Or, we can apply the stochastic chain rule to derive an'
  write ( *, '(a)' ) '  an SDE for V, and solve that.'

  seed = 123456789

  call chain ( seed, n, xem, vem, diff )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum | Sqrt(X) - V | = ', diff

  call chain_gnuplot ( n, xem, vem )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests EM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 256

  real ( kind = 8 ) diff
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(0:n)
  real ( kind = 8 ) t2(0:n/4)
  real ( kind = 8 ) xtrue(0:n)
  real ( kind = 8 ) xem(0:n/4)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  EM solves a stochastic differential equation'
  write ( *, '(a)' ) '  using the Euler-Maruyama method.'

  seed = 123456789

  call em ( seed, n, t, xtrue, t2, xem, diff )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  | Exact X(T) - EM X(T) | = ', diff

  call em_gnuplot ( n, t, xtrue, t2, xem )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests EMSTRONG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 100
  integer ( kind = 4 ), parameter :: n = 512
  integer ( kind = 4 ), parameter :: p_max = 6

  real ( kind = 8 ) dtvals(p_max)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) xerr(p_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  EMSTRONG investigates the strong convergence'
  write ( *, '(a)' ) '  of the Euler-Maruyama method.'

  seed = 123456789

  call emstrong ( seed, m, n, p_max, dtvals, xerr )

  call emstrong_gnuplot ( p_max, dtvals, xerr )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests EMWEAK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 50000
  integer ( kind = 4 ), parameter :: p_max = 5

  real ( kind = 8 ) dtvals(p_max)
  integer ( kind = 4 ) method
  integer ( kind = 4 ) seed
  real ( kind = 8 ) xerr(p_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  EMWEAK investigates the weak convergence'
  write ( *, '(a)' ) '  of the Euler-Maruyama method.'

  seed = 123456789
  method = 0

  call emweak ( seed, method, m, p_max, dtvals, xerr )

  call emweak_gnuplot ( p_max, dtvals, xerr, method )

  seed = 123456789
  method = 1

  call emweak ( seed, method, m, p_max, dtvals, xerr )

  call emweak_gnuplot ( p_max, dtvals, xerr, method )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests MILSTRONG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: p_max = 4

  real ( kind = 8 ) dtvals(p_max)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) xerr(p_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  MILSTRONG investigates the strong convergence'
  write ( *, '(a)' ) '  of the Milstein method.'

  seed = 123456789

  call milstrong ( seed, p_max, dtvals, xerr )

  call milstrong_gnuplot ( p_max, dtvals, xerr )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests STAB_ASYMPTOTIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 1000
  integer ( kind = 4 ), parameter :: p_max = 3
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  STAB_ASYMPTOTIC investigates the asymptotic'
  write ( *, '(a)' ) '  stability of the Euler-Maruyama method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For technical reasons, the plotting is done'
  write ( *, '(a)' ) '  in the same routine as the computations.'

  seed = 123456789

  call stab_asymptotic ( seed, n, p_max )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests STAB_MEANSQUARE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09:'
  write ( *, '(a)' ) '  STAB_MEANSQUARE investigates the mean square'
  write ( *, '(a)' ) '  stability of the Euler-Maruyama method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For technical reasons, the plotting is done'
  write ( *, '(a)' ) '  in the same routine as the computations.'

  seed = 123456789

  call stab_meansquare ( seed )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests STOCHASTIC_INTEGRAL_ITO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10:'
  write ( *, '(a)' ) '  Estimate the Ito integral of W(t) dW over [0,1].'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '                                                 Abs          Rel'
  write ( *, '(a)' ) &
    '         N        Exact        Estimate          Error        Error'
  write ( *, '(a)' ) ' '

  n = 100
  seed = 123456789

  do i = 1, 7

    call stochastic_integral_ito ( n, seed, estimate, exact, error )

    write ( *, '(2x,i8,2x,g16.8,2x,g16.8,2x,g10.2,2x,g10.2)' ) &
      n, exact, estimate, error, error / exact

    n = n * 4

  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests STOCHASTIC_INTEGRAL_STRAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11:'
  write ( *, '(a)' ) &
    '  Estimate the Stratonovich integral of W(t) dW over [0,1].'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '                                                 Abs          Rel'
  write ( *, '(a)' ) &
    '         N        Exact        Estimate          Error        Error'
  write ( *, '(a)' ) ' '

  n = 100
  seed = 123456789

  do i = 1, 7

    call stochastic_integral_strat ( n, seed, estimate, exact, error )

    write ( *, '(2x,i8,2x,g16.8,2x,g16.8,2x,g10.2,2x,g10.2)' ) &
      n, exact, estimate, error, error / exact

    n = n * 4

  end do

  return
end


