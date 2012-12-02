program main

!*****************************************************************************80
!
!! MAIN is the main program for STEAM_INTERACT.
!
!  Discussion:
!
!    STEAM_INTERACT interactively queries the NBS/NRC Steam Table routines.
!
!    The units used internally are:
!
!    Density      G/cm**3   grams per cubic centimeter
!    Enthalpy     J/G       joules per gram
!    Entropy      J/(GK)    joules per gram degree Kelvin
!    Length       M         meters
!    Pressure     MPa       megapascals
!    Temperature  K         degrees Kelvin
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) c
  real ( kind = 8 ) cjth
  real ( kind = 8 ) cjtt
  real ( kind = 8 ) cp
  real ( kind = 8 ) cv
  logical ok
  real ( kind = 8 ) dpdr
  real ( kind = 8 ) dpdt
  real ( kind = 8 ) dvdr
  real ( kind = 8 ) dvdt
  real ( kind = 8 ) eps
  real ( kind = 8 ) eta
  real ( kind = 8 ) g
  real ( kind = 8 ) gascon
  real ( kind = 8 ) gl
  real ( kind = 8 ) gv
  real ( kind = 8 ) h
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) itmax
  character ( len = 72 ) label
  real ( kind = 8 ) lambda
  real ( kind = 8 ) p
  real ( kind = 8 ), parameter :: p_crit = 22.055D+00
  real ( kind = 8 ), parameter :: p_max = 22.00D+00
  real ( kind = 8 ), parameter :: p_min = 0.0D00
  real ( kind = 8 ) pl
  real ( kind = 8 ) pr
  real ( kind = 8 ) ps
  real ( kind = 8 ) pv
  real ( kind = 8 ) rho
  real ( kind = 8 ) rhol
  real ( kind = 8 ), parameter :: rho_max = 1.9D+00
  real ( kind = 8 ), parameter :: rho_min = 0.0D+00
  real ( kind = 8 ) rho_start
  real ( kind = 8 ) rhov
  real ( kind = 8 ) s
  real ( kind = 8 ) sigma
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: t_crit = 647.1260000001D+00
  real ( kind = 8 ), parameter :: t_max = 647.0D+00
  real ( kind = 8 ), parameter :: t_min = 0.0D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) vir

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STEAM_INTERACT:'
  write ( *, '(a)' ) '  Interactive queries of the NBS/NRC Steam Table Program.'
!
!  Read the next choice.
!
  do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter option'
    write ( *, '(a)' ) '  1, table for given density and temperature'
    write ( *, '(a)' ) '  2, table for given pressure and temperature'
    write ( *, '(a)' ) '  3, saturation table for given temperature'
    write ( *, '(a)' ) '  4, saturation table for given pressure'
    write ( *, '(a)' ) '  7, reset error tolerance and iteration max.'
    write ( *, '(a)' ) '  0, to quit.'

    read ( *, *, iostat = ios ) iopt

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Input error!'
      exit
    end if
!
!***********************************************************************
!  Option 0, Quit.
!***********************************************************************
!
    if ( iopt == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Shutdown request'
      exit
!
!***********************************************************************
!  Option 1, properties for given density and temperature.
!***********************************************************************
!
    else if ( iopt == 1 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Compute table for given density and temperature.'
      write ( *, '(a)' ) ' '

      do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter the density in grams/cubic centimeter,'
        write ( *, '(a,g14.6)' ) '  greater than ', rho_min
        write ( *, '(a,g14.6)' ) '  and less than or equal to ', rho_max

        read ( *, *, iostat = ios ) rho

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'Input error!'
          exit
        end if

        call d_gtle_range ( rho, rho_min, rho_max, ok )

        if ( ok ) then
          exit
        end if

      end do

      if ( ios /= 0 ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 'Specified density is ', rho

      do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter the temperature in degrees Kelvin,'
        write ( *, '(a,g14.6)' ) '  greater than ', t_min
        write ( *, '(a,g14.6)' ) '  and less than or equal to ', t_max

        read ( *, *, iostat = ios ) t

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'Input error!'
          exit
        end if

        call d_gtle_range ( t, t_min, t_max, ok )

        if ( ok ) then
          exit
        end if

      end do

      if ( ios /= 0 ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 'Specified temperature is ', t
!
!  Compute thermodynamic quantities as function of T and D.
!
      call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

      call dielectric ( t, rho, eps )
      call prandtl ( t, p, pr )
      call secvir ( t, vir )
      call sound ( t, p, c )
      call surten ( t, sigma )
      call thercon ( t, rho, lambda )
      call viscosity ( t, rho, eta )
      call volume ( t, rho, v, dvdt, dvdr )

      label = 'Temperature and density input.'
      call output ( a, c, cjth, cjtt, cp, cv, rho, dpdr, dpdt, dvdr, &
        dvdt, eps, eta, g, h, label, lambda, p, pr, s, sigma, t, u, v, vir )
!
!  Check for metastability.
!
!  Compute the saturation pressure, and liquid and gas saturation densities.
!  Compute residual and base functions, and liquid and gas pressures.
!
      if ( t < t_crit ) then

        write ( *, * ) 'Call PSAT'
        rhol = 0.0D+00
        rhov = 0.0D+00
        call psat ( t, ps, rhol, rhov )

        rho = rhol

        call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

        call dielectric ( t, rho, eps )
        call prandtl ( t, p, pr )
        call secvir ( t, vir )
        call sound ( t, p, c )
        call surten ( t, sigma )
        call thercon ( t, rho, lambda )
        call viscosity ( t, rho, eta )
        call volume ( t, rho, v, dvdt, dvdr )

        label = 'Liquid phase:'

        call output ( a, c, cjth, cjtt, cp, cv, rho, dpdr, dpdt, dvdr, &
          dvdt, eps, eta, g, h, label, lambda, p, pr, s, sigma, t, u, v, vir )

        gl = g
        pl = p
        rho = rhov

        call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

        call dielectric ( t, rho, eps )
        call prandtl ( t, p, pr )
        call secvir ( t, vir )
        call sound ( t, p, c )
        call surten ( t, sigma )
        call thercon ( t, rho, lambda )
        call viscosity ( t, rho, eta )
        call volume ( t, rho, v, dvdt, dvdr )

        label = 'Vapor phase:'

        call output ( a, c, cjth, cjtt, cp, cv, rho, dpdr, dpdt, dvdr, &
          dvdt, eps, eta, g, h, label, lambda, p, pr, s, sigma, t, u, v, vir )

        pv = p
        gv = g

        if ( rhov < rho .and. rho < rhol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'Caution!'
          write ( *, '(a)' ) &
            '  The density and temperature define a metastable state.'
          write ( *, '(a,2g14.6)' ) 'Pressure (MPa)  ', pl, pv
          write ( *, '(a,2g14.6)' ) 'Gibbs (KJ/Kg)   ', gl, gv
          write ( *, '(a,2g14.6)' ) 'Density (g/cm^3)', rhol, rhov
        end if

      end if
!
!***********************************************************************
!  Option 2, properties for given pressure and temperature.
!***********************************************************************
!
    else if ( iopt == 2 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Compute table for given pressure and temperature.'
      write ( *, '(a)' ) ' '

      do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter the pressure in MPa'
        write ( *, '(a,g14.6)' ) 'greater than ', p_min
        write ( *, '(a,g14.6)' ) 'and less than or equal to ', p_max

        read ( *, *, iostat = ios ) p

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'Input error!'
          exit
        end if

        call d_gtle_range ( p, p_min, p_max, ok )

        if ( ok ) then
          exit
        end if

      end do

      if ( ios /= 0 ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Specified pressure is ', p

      do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter the temperature in degrees Kelvin'
        write ( *, '(a,g14.6)' ) 'greater than ', t_min
        write ( *, '(a,g14.6)' ) 'and less than or equal to ', t_max

        read ( *, *, iostat = ios ) t

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'Input error!'
          exit
        end if

        call d_gtle_range ( t, t_min, t_max, ok )

        if ( ok ) then
          exit
        end if

      end do

      if ( ios /= 0 ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Specified temperature is ', t
!
!  Get the density as function of P and T.
!
      ps = 20000.0D+00
      rhol = 0.0D+00
      rhov = 0.0D+00

      if ( t < t_crit ) then
        call psat ( t, ps, rhol, rhov )
      end if

      if ( p > ps ) then
        rho_start = rhol
      else
        rho_start = p / ( gascon() * t )
      end if

      call dense ( p, t, rho_start, rho, dpdr )
!
!  Compute properties as function of T and D.
!
      call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

      call dielectric ( t, rho, eps )
      call prandtl ( t, p, pr )
      call secvir ( t, vir )
      call sound ( t, p, c )
      call surten ( t, sigma )
      call thercon ( t, rho, lambda )
      call viscosity ( t, rho, eta )
      call volume ( t, rho, v, dvdt, dvdr )

      label = 'Pressure and temperature were input.'

      call output ( a, c, cjth, cjtt, cp, cv, rho, dpdr, dpdt, dvdr, &
        dvdt, eps, eta, g, h, label, lambda, p, pr, s, sigma, t, u, v, vir )
!
!***********************************************************************
!  Option 3, Saturation properties as a function of T.
!***********************************************************************
!
    else if ( iopt == 3 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Compute saturation table for given temperature.'
      write ( *, '(a)' ) ' '

      do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter the temperature in degrees Kelvin'
        write ( *, '(a,g14.6)' ) 'greater than ', t_min
        write ( *, '(a,g14.6)' ) 'and less than or equal to ', t_crit

        read ( *, *, iostat = ios ) t

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Input error while reading T!'
          exit
        end if

        call d_gtle_range ( t, t_min, t_crit, ok )

        if ( ok ) then
          exit
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The specified temperature is ', t

      if ( ios /= 0 ) then
        exit
      end if

      rhol = 0.0D+00
      rhov = 0.0D+00
      call psat ( t, p, rhol, rhov )
!
!  Calculate saturated liquid values as a function of T and DLL.
!
      rho = rhol

      call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

      call dielectric ( t, rho, eps )
      call prandtl ( t, p, pr )
      call secvir ( t, vir )
      call sound ( t, p, c )
      call surten ( t, sigma )
      call thercon ( t, rho, lambda )
      call viscosity ( t, rho, eta )
      call volume ( t, rho, v, dvdt, dvdr )

      label = &
        'Saturated liquid properties as function of temperature.'

      call output ( a, c, cjth, cjtt, cp, cv, rho, dpdr, dpdt, dvdr, &
        dvdt, eps, eta, g, h, label, lambda, p, pr, s, sigma, t, u, v, vir )
!
!  Calculate saturated vapor values, as a function of T and DVV.
!
      rho = rhov

      call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

      call dielectric ( t, rho, eps )
      call prandtl ( t, p, pr )
      call secvir ( t, vir )
      call sound ( t, p, c )
      call surten ( t, sigma )
      call thercon ( t, rho, lambda )
      call viscosity ( t, rho, eta )
      call volume ( t, rho, v, dvdt, dvdr )

      label = &
        'Saturated vapor properties as a function of temperature.'

      call output ( a, c, cjth, cjtt, cp, cv, rho, dpdr, dpdt, dvdr, &
        dvdt, eps, eta, g, h, label, lambda, p, pr, s, sigma, t, u, v, vir )
!
!***********************************************************************
!  Option 4, Saturation properties as a function of P.
!***********************************************************************
!
    else if ( iopt == 4 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Compute saturation table for given pressure.'
      write ( *, '(a)' ) ' '

      do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter the pressure in MPa'
        write ( *, '(a,g14.6)' ) 'greater than ', p_min
        write ( *, '(a,g14.6)' ) 'and less than or equal to ', p_crit

        read ( *, *, iostat = ios ) p

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'Input error!'
          exit
        end if

        call d_gtle_range ( p, p_min, p_crit, ok )

        if ( ok ) then
          exit
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Specified pressure is P = ', p

      if ( ios /= 0 ) then
        exit
      end if

      rhol = 0.0D+00
      rhov = 0.0D+00
      call tsat ( t, p, rhol, rhov )
!
!  Calculate saturated liquid values as a function of T and DLL.
!
      rho = rhol

      call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

      call dielectric ( t, rho, eps )
      call prandtl ( t, p, pr )
      call secvir ( t, vir )
      call sound ( t, p, c )
      call surten ( t, sigma )
      call thercon ( t, rho, lambda )
      call viscosity ( t, rho, eta )
      call volume ( t, rho, v, dvdt, dvdr )

      label = 'Saturated liquid properties as a function of pressure.'

      call output ( a, c, cjth, cjtt, cp, cv, rho, dpdr, dpdt, dvdr, &
        dvdt, eps, eta, g, h, label, lambda, p, pr, s, sigma, t, u, v, vir )
!
!  Calculate saturated vapor values as a function of T and DVV.
!
      rho = rhov

      call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

      call dielectric ( t, rho, eps )
      call prandtl ( t, p, pr )
      call secvir ( t, vir )
      call sound ( t, p, c )
      call surten ( t, sigma )
      call thercon ( t, rho, lambda )
      call viscosity ( t, rho, eta )
      call volume ( t, rho, v, dvdt, dvdr )

      label = 'Saturated vapor properties as a function of pressure.'

      call output ( a, c, cjth, cjtt, cp, cv, rho, dpdr, dpdt, dvdr, &
        dvdt, eps, eta, g, h, label, lambda, p, pr, s, sigma, t, u, v, vir )
!
!***********************************************************************
!  Option 7, Specify convergence tolerances.
!***********************************************************************
!
    else if ( iopt == 7 ) then

      itmax = 10

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Current iteration maximum is ', itmax
      write ( *, '(a)' ) '  Enter a new value.'

      read ( *, *, iostat = ios ) itmax

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Input error!'
        exit
      end if

!     call machin ( itmax, 1 )
!
!***********************************************************************
!  Unrecognized option.
!***********************************************************************
!
    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Unrecognized option!'
      write ( *, '(a,i6)' ) '  IOPT = ', iopt

    end if

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STEAM_INTERACT:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine output ( a, c, cjth, cjtt, cp, cv, rho, dpdr, dpdt, dvdr, &
  dvdt, eps, eta, g, h, label, lambda, p, pr, s, sigma, t, u, v, vir )

!***********************************************************************
!
!! OUTPUT prints the values of the thermodynamic quantities.
!
!  Modified:
!
!    26 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the Helmholtz function.
!
!    Input, real ( kind = 8 ) C, the speed of sound.
!
!    Input, real ( kind = 8 ) CJTH, the Joule-Thomson coefficient.
!
!    Input, real ( kind = 8 ) CJTT, the isothermal Joule-Thomson coefficient.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) c
  real ( kind = 8 ) cjth
  real ( kind = 8 ) cjtt
  real ( kind = 8 ) cp
  real ( kind = 8 ) cv
  real ( kind = 8 ) dpdr
  real ( kind = 8 ) dpdt
  real ( kind = 8 ) dvdr
  real ( kind = 8 ) dvdt
  real ( kind = 8 ) eps
  real ( kind = 8 ) eta
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  character ( len = 72 ) label
  real ( kind = 8 ) lambda
  real ( kind = 8 ) p
  real ( kind = 8 ) pr
  real ( kind = 8 ) rho
  real ( kind = 8 ) s
  real ( kind = 8 ) sigma
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) vir

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( label )
  write ( *, '(a)' )' '
  write ( *, '(a,g16.8,a)' ) &
    'A =      ', a,       ' KJ/kg        (Helmholtz function)'
  write ( *, '(a,g16.8,a)' ) &
    'C =      ', c,       ' M/s          (speed of sound)'
  write ( *, '(a,g16.8,a)' ) &
    'CJTH =   ', cjth,    ' K/MPa        (Joule-Thomson coefficient)'
  write ( *, '(a,g16.8,a)' ) &
    'CJTT =   ', cjtt,    ' CM3/g        (Isothermal Joule-Thomson coefficient)'
  write ( *, '(a,g16.8,a)' ) &
    'CP =     ', cp,      ' KJ/(kg K)    (Specific heat at constant pressure)'
  write ( *, '(a,g16.8,a)' ) &
    'CV =     ', cv,      ' KJ/(kg K)    (Specific heat at constant volume)'
  write ( *, '(a,g16.8,a)' ) &
    'RHO =    ', rho,     ' g/CM3        (density)'
  write ( *, '(a,g16.8,a)' ) &
    'DPDR =   ', dpdr,    ' MPa CM3/g    (D Pressure/D Density)'
  write ( *, '(a,g16.8,a)' ) &
    'DPDT =   ', dpdt,    ' MPa/K        (D Pressure/D Temperature)'
  write ( *, '(a,g16.8,a)' ) &
    'DVDR =   ', dvdr,    ' CM6/g2       (d Specific Volume/d Density)'
  write ( *, '(a,g16.8,a)' ) &
    'DVDT =   ', dvdt,    ' CM3/(g K)    (d Specific Volume/ d T)'
  write ( *, '(a,g16.8,a)' ) &
    'EPS =    ', eps,     ' 1            (static dielectric constant)'
  write ( *, '(a,g16.8,a)' ) &
    'ETA =    ', eta,     ' MPa s        (viscosity)'
  write ( *, '(a,g16.8,a)' ) &
    'G =      ', g,       ' KJ/kg        (Gibbs specific energy)'
  write ( *, '(a,g16.8,a)' ) &
    'H =      ', h,       ' KJ/kg        (enthalpy)'
  write ( *, '(a,g16.8,a)' ) &
    'LAMBDA = ', lambda,  ' mW/(m K)     (thermal conductivity)'
  write ( *, '(a,g16.8,a)' ) &
    'P =      ', p,       ' MPa          (pressure)'
  write ( *, '(a,g16.8,a)' ) &
    'PR =     ', pr,      ' 1            (Prandtl number)'
  write ( *, '(a,g16.8,a)' ) &
    'S =      ', s,       ' KJ/(kg K)    (entropy)'
  write ( *, '(a,g16.8,a)' ) &
    'SIGMA =  ', sigma,   ' Pa M         (surface tension)'
  write ( *, '(a,g16.8,a)' ) &
    'T =      ', t,       ' K            (temperature)'
  write ( *, '(a,g16.8,a)' ) &
    'U =      ', u,       ' KJ/kg        (internal energy)'
  write ( *, '(a,g16.8,a)' ) &
    'V =      ', v,       ' CM3/g        (specific volume)'
  write ( *, '(a,g16.8,a)' ) &
    'VIR =    ', vir,     ' CM3/g        (second virial coefficient)'

  return
end
subroutine d_gtle_range ( x, a, b, ok )

!***********************************************************************
!
!! D_GTLE_RANGE checks that A < X <= B.
!
!  Modified:
!
!    31 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value to be checked.
!
!    Input, real ( kind = 8 ) A, B, the left and right endpoints.
!
!    Output, logical OK, is TRUE if A < X <= B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical ok
  real ( kind = 8 ) x

  if ( x <= a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The input quantity is too small.'
    write ( *, '(a,g14.6)' ) '  The minimum legal value is ', a
    write ( *, '(a,g14.6)' ) '  The input value is         ', x
    ok = .false.
  else if ( b < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The input quantity is too large.'
    write ( *, '(a,g14.6)' ) '  The maximum legal value is ', a
    write ( *, '(a,g14.6)' ) '  The input value is         ', x
    ok = .false.
  else
    ok = .true.
  end if

  return
end
