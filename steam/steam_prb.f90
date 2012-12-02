program main

!*****************************************************************************80
!
!! MAIN is the main program for STEAM_PRB.
!
!  Discussion:
!
!    STEAM_PRB tests the STEAM routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STEAM_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the STEAM library.'

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
  call test12 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STEAM_PRB'
  write ( *, '(a)' ) '  Normal end of STEAM tests.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests the computation of CP.
!
!  Discussion:
!
!    See the table of page 229 of the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cjth
  real ( kind = 8 ) cjtt
  real ( kind = 8 ) cp
  real ( kind = 8 ) cp2
  real ( kind = 8 ) cv
  real ( kind = 8 ) dpdr
  real ( kind = 8 ) dpdt
  real ( kind = 8 ) error
  real ( kind = 8 ) g
  real ( kind = 8 ) gascon
  real ( kind = 8 ) h
  integer ( kind = 4 ) n
  real ( kind = 8 ) pbar
  real ( kind = 8 ) pmpa
  real ( kind = 8 ) ps
  real ( kind = 8 ) rho
  real ( kind = 8 ) rhol
  real ( kind = 8 ) rhov
  real ( kind = 8 ) rho_start
  real ( kind = 8 ) s
  real ( kind = 8 ) tc
  real ( kind = 8 ) tk
  real ( kind = 8 ), parameter :: t_crit = 647.1260000001D+00
  real ( kind = 8 ) u

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  THERM computes CP, the specific heat capacity'
  write ( *, '(a)' ) '  at constant pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     P            T          Rho         CP          CP'
  write ( *, '(a)' ) &
    '     (bar)        (C)        g/cm3       tabulated   computed'
  write ( *, '(a)' ) ' '

  error = 0.0D+00
  n = 0

  do

    call cp_values ( n, tc, pbar, cp )

    if ( n == 0 ) then
      exit
    end if

    pmpa = pbar / 10.0D+00
    tk = tc + 273.15D+00
!
!  Determine the density RHO.
!
    ps = 20000.0D+00
    rhol = 0.0D+00
    rhov = 0.0D+00

    if ( tk < t_crit ) then
      call psat ( tk, ps, rhol, rhov )
    end if

    if ( pmpa > ps ) then
      rho_start = rhol
    else
      rho_start = pmpa / ( gascon() * tk )
    end if

    call dense ( pmpa, tk, rho_start, rho, dpdr )
!
!  Get CP (and all the other thermodynamic quantities).
!
    call therm ( tk, rho, a, cjth, cjtt, cp2, cv, dpdr, dpdt, g, h, pmpa, s, u )

    error = max ( error, abs ( cp - cp2 ) )

    write ( *, '(3f12.6,2f14.8)' ) pbar, tc, rho, cp, cp2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum absolute error was ', error

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests DENSE at roughly atmospheric pressure.
!
!  Discussion:
!
!    See the table of page 34 of the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
  implicit none

  integer ( kind = 4 ), parameter :: ntest = 39

  real ( kind = 8 ) dpdr
  integer ( kind = 4 ) i
  real ( kind = 8 ) pbar
  real ( kind = 8 ) pmpa
  real ( kind = 8 ) rho
  real ( kind = 8 ) rhol
  real ( kind = 8 ) rhov
  real ( kind = 8 ) rho_start
  real ( kind = 8 ) t2
  real ( kind = 8 ) tc(ntest)
  real ( kind = 8 ) tk

  pbar = 1.0D+00
  pmpa = pbar / 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  DENSE computes fluid density.'
  write ( *, '(a,g14.6,a)' ) '  Pressure is fixed at ', pmpa, ' MegaPascals.'

  call tsat ( pmpa, t2, rhol, rhov )
  write ( *, '(a,g14.6,a)' ) '  The saturation temperature is ', t2, &
    ' degrees Kelvin'

  tc(1) = 0.0D+00
  tc(2) = 5.0D+00
  tc(3) = 10.0D+00
  tc(4) = 15.0D+00
  tc(5) = 20.0D+00
  tc(6) = 25.0D+00
  tc(7) = 30.0D+00
  tc(8) = 35.0D+00
  tc(9) = 40.0D+00
  tc(10) = 45.0D+00
  tc(11) = 50.0D+00
  tc(12) = 55.0D+00
  tc(13) = 60.0D+00
  tc(14) = 65.0D+00
  tc(15) = 70.0D+00
  tc(16) = 75.0D+00
  tc(17) = 80.0D+00
  tc(18) = 85.0D+00
  tc(19) = 90.0D+00
  tc(20) = 95.0D+00
  tc(21) = 100.0D+00
  tc(22) = 110.0D+00
  tc(23) = 120.0D+00
  tc(24) = 130.0D+00
  tc(25) = 140.0D+00
  tc(26) = 150.0D+00
  tc(27) = 160.0D+00
  tc(28) = 170.0D+00
  tc(29) = 180.0D+00
  tc(30) = 190.0D+00
  tc(31) = 200.0D+00
  tc(32) = 300.0D+00
  tc(33) = 400.0D+00
  tc(34) = 500.0D+00
  tc(35) = 600.0D+00
  tc(36) = 700.0D+00
  tc(37) = 800.0D+00
  tc(38) = 900.0D+00
  tc(39) = 1000.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T (C)         T (K)   RHO (KG/M3)          dPdR'
  write ( *, '(a)' ) ' '

  do i = 1, ntest

    tk = tc(i) + 273.15D+00

    if ( tk < t2 ) then
      rho_start = 1.9D+00
    else
      rho_start = 0.01D+00
    end if

    call dense ( pmpa, tk, rho_start, rho, dpdr )

    rho = rho * 1000.0D+00
    dpdr = dpdr / 1000.0D+00
    write ( *, '(2g14.6,2g16.8)' ) tc(i), tk, rho, dpdr

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests DIELEC.
!
!  Discussion:
!
!    See the table of page 266 of the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
  implicit none

  real ( kind = 8 ) dpdr
  real ( kind = 8 ) eps
  real ( kind = 8 ) eps2
  real ( kind = 8 ) error
  real ( kind = 8 ) gascon
  integer ( kind = 4 ) n
  real ( kind = 8 ) pbar
  real ( kind = 8 ) pmpa
  real ( kind = 8 ) ps
  real ( kind = 8 ) rho
  real ( kind = 8 ) rhol
  real ( kind = 8 ) rhov
  real ( kind = 8 ) rho_start
  real ( kind = 8 ), parameter :: t_crit = 647.1260000001D+00
  real ( kind = 8 ) tc
  real ( kind = 8 ) tk

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  DIELEC computes the dielectric constant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     P            T          Rho         EPS         EPS'
  write ( *, '(a)' ) &
    '     (bar)        (C)        g/cm3       tabulated   computed'
  write ( *, '(a)' ) ' '

  error = 0.0D+00
  n = 0

  do

    call dielectric_values ( n, tc, pbar, eps )

    if ( n == 0 ) then
      exit
    end if

    pmpa = pbar / 10.0D+00
    tk = tc + 273.15D+00
!
!  Determine the density RHO.
!
    ps = 20000.0D+00
    rhol = 0.0D+00
    rhov = 0.0D+00

    if ( tk < t_crit ) then
      call psat ( tk, ps, rhol, rhov )
    end if

    if ( pmpa > ps ) then
      rho_start = rhol
    else
      rho_start = pmpa / ( gascon() * tk )
    end if

    call dense ( pmpa, tk, rho_start, rho, dpdr )

    call dielectric ( tk, rho, eps2 )

    error = max ( error, abs ( eps - eps2 ) )

    write ( *, '(3f12.6,2f14.8)' ) pbar, tc, rho, eps, eps2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum absolute error = ', error

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests PRANDTL.
!
!  Discussion:
!
!    See the table of page 267 of the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
  implicit none

  real ( kind = 8 ) error
  integer ( kind = 4 ) n
  real ( kind = 8 ) pbar
  real ( kind = 8 ) pmpa
  real ( kind = 8 ) pr
  real ( kind = 8 ) pr2
  real ( kind = 8 ) tc
  real ( kind = 8 ) tk

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  PRANDTL computes the Prandtl number'
  write ( *, '(a)' ) '  given the pressure and temperature.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     P            T          PR          PR'
  write ( *, '(a)' ) &
    '     (bar)        (C)        tabulated   computed'
  write ( *, '(a)' ) ' '

  error = 0.0D+00
  n = 0

  do

    call prandtl_values ( n, tc, pbar, pr )

    if ( n == 0 ) then
      exit
    end if

    pmpa = pbar / 10.0D+00
    tk = tc + 273.15D+00
    call prandtl ( tk, pmpa, pr2 )

    error = max ( error, abs ( pr - pr2 ) )

    write ( *, '(2f12.6,2f14.8)' ) pbar, tc, pr, pr2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum absolute error = ', error

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests PSAT.
!
!  Discussion:
!
!    See the table of pages 9-15 of the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
  implicit none

  real ( kind = 8 ) error
  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ) p2
  real ( kind = 8 ) rhol
  real ( kind = 8 ) rhov
  real ( kind = 8 ) tc
  real ( kind = 8 ) tk

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  PSAT computes the saturation pressure'
  write ( *, '(a)' ) '  given the temperature.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     T          P           P'
  write ( *, '(a)' ) '     (C)        (bar)       (bar)'
  write ( *, '(a)' ) '                tabulated   computed'
  write ( *, '(a)' ) ' '

  error = 0.0D+00
  n = 0

  do

    call psat_values ( n, tc, p )

    if ( n == 0 ) then
      exit
    end if

    tk = tc + 273.15D+00

    call psat ( tk, p2, rhol, rhov )
    p2 = p2 * 10.0D+00

    error = max ( error, abs ( p - p2 ) )

    write ( *, '(f12.6,2f14.8)' ) tc, p, p2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum absolute error = ', error

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests SECVIR.
!
!  Discussion:
!
!    See the table on pages 24-25 of the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
  implicit none

  integer ( kind = 4 ), parameter :: ntest = 45

  real ( kind = 8 ) error
  integer ( kind = 4 ) n
  real ( kind = 8 ) tc
  real ( kind = 8 ) tk
  real ( kind = 8 ) vir
  real ( kind = 8 ) vir2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  SECVIR computes the second virial coefficient'
  write ( *, '(a)' ) '  as a function of temperature.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     T          VIR         VIR'
  write ( *, '(a)' ) '     (C)        tabulated   computed'
  write ( *, '(a)' ) ' '

  error = 0.0D+00
  n = 0

  do

    call secvir_values ( n, tc, vir )

    if ( n == 0 ) then
      exit
    end if

    tk = tc + 273.15D+00

    call secvir ( tk, vir2 )

    error = max ( error, abs ( vir - vir2 ) )

    write ( *, '(f12.6,2f14.8)' ) tc, vir, vir2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum absolute error = ', error

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests SOUND.
!
!  Discussion:
!
!    See the table on page 238-246 of the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
  implicit none

  integer ( kind = 4 ), parameter :: nptest = 3
  integer ( kind = 4 ), parameter :: nttest = 20

  real ( kind = 8 ) c
  real ( kind = 8 ) c2
  real ( kind = 8 ) error
  integer ( kind = 4 ) n
  real ( kind = 8 ) pbar
  real ( kind = 8 ) pmpa
  real ( kind = 8 ) tc
  real ( kind = 8 ) tk

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  SOUND computes the speed of sound.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     P            T          C           C'
  write ( *, '(a)' ) &
    '     (bar)        (C)        tabulated   computed'
  write ( *, '(a)' ) ' '

  error = 0.0D+00
  n = 0

  do

    call sound_values ( n, tc, pbar, c )

    if ( n == 0 ) then
      exit
    end if

    pmpa = pbar / 10.0D+00
    tk = tc + 273.15D+00

    call sound ( tk, pmpa, c2 )

    error = max ( error, abs ( c - c2 ) )

    write ( *, '(2f12.6,2f14.8)' ) pbar, tc, c, c2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' )  '  Maximum absolute error = ', error

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests SURTEN.
!
!  Discussion:
!
!    See the table of page 267 of the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
  implicit none

  real ( kind = 8 ) error
  integer ( kind = 4 ) n
  real ( kind = 8 ) sigma
  real ( kind = 8 ) sigma2
  real ( kind = 8 ) tc
  real ( kind = 8 ) tk

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  SURTEN computes the surface tension'
  write ( *, '(a)' ) '  given the temperature.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     T          SIGMA       SIGMA'
  write ( *, '(a)' ) '     (C)        tabulated   computed'
  write ( *, '(a)' ) ' '

  error = 0.0D+00
  n = 0

  do

    call surten_values ( n, tc, sigma )

    if ( n == 0 ) then
      exit
    end if

    tk = tc + 273.15D+00

    call surten ( tk, sigma2 )

    error = max ( error, abs ( sigma - sigma2 ) )

    write ( *, '(f12.6,2f14.8)' ) tc, sigma, sigma2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum absolute error = ', error

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests THERCON.
!
!  Discussion:
!
!    See the table of page 264 of the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
  implicit none

  real ( kind = 8 ) dpdr
  real ( kind = 8 ) error
  real ( kind = 8 ) gascon
  integer ( kind = 4 ) n
  real ( kind = 8 ) lambda
  real ( kind = 8 ) lambda2
  real ( kind = 8 ) pbar
  real ( kind = 8 ) pmpa
  real ( kind = 8 ) ps
  real ( kind = 8 ) rho
  real ( kind = 8 ) rho_start
  real ( kind = 8 ) rhol
  real ( kind = 8 ) rhov
  real ( kind = 8 ), parameter :: t_crit = 647.1260000001D+00
  real ( kind = 8 ) tc
  real ( kind = 8 ) tk

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  THERCON computes thermal conductivity.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     P            T          LAMBDA      LAMBDA'
  write ( *, '(a)' ) &
    '     (bar)        (C)        tabulated   computed'
  write ( *, '(a)' ) ' '

  error = 0.0D+00
  n = 0

  do

    call thercon_values ( n, tc, pbar, lambda )

    if ( n == 0 ) then
      exit
    end if

    pmpa = pbar / 10.0D+00
    tk = tc + 273.15D+00

    rho_start = 1.9D+00
    call dense ( pmpa, tk, rho_start, rho, dpdr )
!
!  Determine the density RHO.
!
    ps = 20000.0D+00
    rhol = 0.0D+00
    rhov = 0.0D+00

    if ( tk < t_crit ) then
      call psat ( tk, ps, rhol, rhov )
    end if

    if ( pmpa > ps ) then
      rho_start = rhol
    else
      rho_start = pmpa / ( gascon() * tk )
    end if

    call dense ( pmpa, tk, rho_start, rho, dpdr )

    call thercon ( tk, rho, lambda2 )

    error = max ( error, abs ( lambda - lambda2 ) )

    write ( *, '(2f12.6,2f14.8)' ) pbar, tc, lambda, lambda2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum absolute error = ', error

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests TSAT.
!
!  Discussion:
!
!    See the table of pages 16-22 of the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
  implicit none

  real ( kind = 8 ) error
  integer ( kind = 4 ) n
  real ( kind = 8 ) pbar
  real ( kind = 8 ) pmpa
  real ( kind = 8 ) rhol
  real ( kind = 8 ) rhov
  real ( kind = 8 ) tc
  real ( kind = 8 ) tc2
  real ( kind = 8 ) tk

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  TSAT computes the saturation temperature'
  write ( *, '(a)' ) '  given the pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P          T           T'
  write ( *, '(a)' ) '     (bar)      (C)         (C)'
  write ( *, '(a)' ) '                tabulated   computed'
  write ( *, '(a)' ) ' '

  error = 0.0D+00
  n = 0

  do

    call tsat_values ( n, pbar, tc )

    if ( n == 0 ) then
      exit
    end if

    tk = tc + 273.15D+00

    pmpa = pbar / 10.0D+00

    call tsat ( pmpa, tk, rhol, rhov )
    tc2 = tk - 273.15D+00

    error = max ( error, abs ( tc - tc2 ) )

    write ( *, '(f12.6,2f14.8)' ) pbar, tc, tc2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum absolute error = ', error

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests VISCOSITY.
!
!  Discussion:
!
!    See the table on page 263 of the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
  implicit none

  real ( kind = 8 ) dpdr
  real ( kind = 8 ) error
  real ( kind = 8 ) eta
  real ( kind = 8 ) eta2
  integer ( kind = 4 ) n
  real ( kind = 8 ) pbar
  real ( kind = 8 ) pmpa
  real ( kind = 8 ) rho
  real ( kind = 8 ) rho_start
  real ( kind = 8 ) tc
  real ( kind = 8 ) tk

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  VISCOSITY computes the viscosity.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     P            T          RHO         ETA         ETA'
  write ( *, '(a)' ) &
    '     (bar)        (C)        g/cm3       tabulated   computed'
  write ( *, '(a)' ) ' '

  error = 0.0D+00
  n = 0

  do

    call viscosity_values ( n, tc, pbar, eta )

    if ( n == 0 ) then
      exit
    end if

    pmpa = pbar / 10.0D+00
    tk = tc + 273.15D+00
!
!  The call to VISCOSITY is producing a few bad results,
!  for P = 1, T = 100, 200.
!
!  Changing the value of RHO_START seems to affect the behavior
!  of the algorithm...and for the worse.  Try a value of 0.2, for instance.
!
!  I am most discouraged that, again, a single variable is being used
!  to do two things.  I'll get back to this some time later.
!
    rho_start = 1.9D+00
    call dense ( pmpa, tk, rho_start, rho, dpdr )
    call viscosity ( tk, rho, eta2 )

    error = max ( error, abs ( eta - eta2 ) )

    write ( *, '(3f12.6,2f14.8)' ) pbar, tc, rho, eta, eta2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum absolute error = ', error

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests VOLUME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!    for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
  implicit none

  integer ( kind = 4 ), parameter :: nptest = 26
  integer ( kind = 4 ), parameter :: nttest = 22

  real ( kind = 8 ) dpdr
  real ( kind = 8 ) dvdr
  real ( kind = 8 ) dvdt
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) pbar(nptest)
  real ( kind = 8 ) pmpa
  real ( kind = 8 ) rho
  real ( kind = 8 ) rho_start
  real ( kind = 8 ) tc(nttest)
  real ( kind = 8 ) tk
  real ( kind = 8 ) v

  tc(1) = 0.0D+00
  tc(2) = 25.0D+00
  tc(3) = 50.0D+00
  tc(4) = 75.0D+00
  tc(5) = 100.0D+00
  tc(6) = 150.0D+00
  tc(7) = 200.0D+00
  tc(8) = 250.0D+00
  tc(9) = 300.0D+00
  tc(10) = 350.0D+00
  tc(11) = 375.0D+00
  tc(12) = 400.0D+00
  tc(13) = 425.0D+00
  tc(14) = 450.0D+00
  tc(15) = 475.0D+00
  tc(16) = 500.0D+00
  tc(17) = 550.0D+00
  tc(18) = 600.0D+00
  tc(19) = 650.0D+00
  tc(20) = 700.0D+00
  tc(21) = 750.0D+00
  tc(22) = 800.0D+00

  pbar(1) = 1.0D+00
  pbar(2) = 5.0D+00
  pbar(3) = 10.0D+00
  pbar(4) = 25.0D+00
  pbar(5) = 50.0D+00
  pbar(6) = 75.0D+00
  pbar(7) = 100.0D+00
  pbar(8) = 125.0D+00
  pbar(9) = 150.0D+00
  pbar(10) = 175.0D+00
  pbar(11) = 200.0D+00
  pbar(12) = 225.0D+00
  pbar(13) = 250.0D+00
  pbar(14) = 275.0D+00
  pbar(15) = 300.0D+00
  pbar(16) = 350.0D+00
  pbar(17) = 400.0D+00
  pbar(18) = 450.0D+00
  pbar(19) = 500.0D+00
  pbar(20) = 550.0D+00
  pbar(21) = 600.0D+00
  pbar(22) = 650.0D+00
  pbar(23) = 700.0D+00
  pbar(24) = 800.0D+00
  pbar(25) = 900.0D+00
  pbar(26) = 1000.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  VOLUME computes the specific volume.'

  do i = 1, 2
    tk = tc(i) + 273.15D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a,g16.8)' ) '  T(C) = ', tc(i)
    write ( *, '(a,g16.8)' ) '  T(K) = ', tk
    write ( *, '(a)' ) '      RHO             V          dVdT       dVdRHO  '
    write ( *, '(a)' ) ' '
    do j = 1, nptest
      pmpa = pbar(j) / 10.0D+00
      rho_start = 1.9D+00
      call dense ( pmpa, tk, rho_start, rho, dpdr )
      call volume ( tk, rho, v, dvdt, dvdr )
      write ( *, '(g14.6,3g16.8)' ) rho, v, dvdt, dvdr
    end do
  end do

  return
end
