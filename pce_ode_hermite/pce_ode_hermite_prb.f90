program main

!*****************************************************************************80
!
!! MAIN is the main program for PCE_ODE_HERMITE_TEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( );
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PCE_ODE_HERMITE_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test PCE_ODE_HERMITE.'

  call pce_ode_hermite_test01 ( )
  call pce_ode_hermite_test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PCE_ODE_HERMITE_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine pce_ode_hermite_test01 ( )

!*****************************************************************************80
!
!! PCE_ODE_HERMITE_TEST01 runs a test problem with PCE_ODE_HERMITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: np = 4
  integer ( kind = 4 ), parameter :: nt = 200

  real ( kind = 8 ) alpha_mu
  real ( kind = 8 ) alpha_sigma
  integer ( kind = 4 ) i
  real ( kind = 8 ) t(0:nt)
  real ( kind = 8 ) tf
  real ( kind = 8 ) ti
  real ( kind = 8 ) u(0:nt,0:np)
  real ( kind = 8 ) uex(0:nt)
  real ( kind = 8 ) ui

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PCE_ODE_HERMITE_TEST01:'
  write ( *, '(a)' ) '  Call PCE_ODE_HERMITE to compute a polynomial chaos expansion'
  write ( *, '(a)' ) '  for the ODE:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    u'' = - alpha * u,'
  write ( *, '(a)' ) '    u(0) = 1.'

  ti = 0.0D+00
  tf = 2.0D+00
  ui = 1.0D+00
  alpha_mu = 0.0D+00
  alpha_sigma = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Initial time         TI = ', ti
  write ( *, '(a,g14.6)' ) '  Final time           TF = ', tf
  write ( *, '(a,i6)' ) '  Number of time steps NT = ', nt
  write ( *, '(a,g14.6)' ) '  Initial condition    UI = ', ui
  write ( *, '(a,i6)' ) '  Expansion degree     NP = ', np
  write ( *, '(a,g14.6)' ) '  E(ALPHA)       ALPHA_MU = ', alpha_mu
  write ( *, '(a,g14.6)' ) '  STD(ALPHA)  ALPHA_SIGMA = ', alpha_sigma

  call pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, alpha_sigma, t, u )
!
!  Evaluate the exact expected value function.
!
  uex(0:nt) = ui * exp ( t(0:nt)**2 / 2.0D+00 )
!
!  Compare the first computed component against the exact expected value.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' i  T(i)  E(U(T(i)))    U(T(i),0)'
  write ( *, '(a)' ) ' '
  do i = 0, nt, 10
    write ( *, '(2x,i4,2x,f6.3,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, t(i), uex(i), u(i,0), abs ( uex(i) - u(i,0) )
  end do

  return
end
subroutine pce_ode_hermite_test02 ( )

!*****************************************************************************80
!
!! PCE_ODE_HERMITE_TEST02 looks at convergence behavior for a fixed time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nt = 2000

  real ( kind = 8 ) alpha_mu
  real ( kind = 8 ) alpha_sigma
  real ( kind = 8 ) ep(0:5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) np
  real ( kind = 8 ) t(0:nt)
  real ( kind = 8 ) tf
  real ( kind = 8 ) ti
  real ( kind = 8 ), allocatable :: u(:,:)
  real ( kind = 8 ) uex(0:nt)
  real ( kind = 8 ) uexf
  real ( kind = 8 ) ui

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PCE_ODE_HERMITE_TEST02:'
  write ( *, '(a)' ) '  Examine convergence behavior as the PCE degree increases:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    u'' = - alpha * u,'
  write ( *, '(a)' ) '    u(0) = 1.'

  ti = 0.0D+00
  tf = 2.0D+00
  ui = 1.0D+00
  alpha_mu = 0.0D+00
  alpha_sigma = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Initial time         TI = ', ti
  write ( *, '(a,g14.6)' ) '  Final time           TF = ', tf
  write ( *, '(a,i6)' ) '  Number of time steps NT = ', nt
  write ( *, '(a,g14.6)' ) '  Initial condition    UI = ', ui
  write ( *, '(a,g14.6)' ) '  E(ALPHA)       ALPHA_MU = ', alpha_mu
  write ( *, '(a,g14.6)' ) '  STD(ALPHA)  ALPHA_SIGMA = ', alpha_sigma

  uexf = ui * exp ( tf**2 / 2.0D+00 )

  do np = 0, 5

    allocate ( u(0:nt,0:np) )

    call pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, alpha_sigma, t, u )

    ep(np) = abs ( uexf - u(nt,0) )

    deallocate ( u )

  end do
!
!  Print error in expected value as a function of the PCE degree.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    NP     Error(NP)     Log(Error(NP))'
  write ( *, '(a)' ) ' '
  do np = 0, 5
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) np, ep(np), log ( ep(np) )
  end do

  return
end
