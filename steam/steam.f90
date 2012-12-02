subroutine base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )

!*****************************************************************************80
!
!! BASE calculates quantities associated with the base Helmholtz function.
!
!  Discussion:
!
!    The equation for the base Helmholtz function AB(T,RHO) is:
!
!      AB(T,RHO) = R * T * (
!        - ln ( 1 - y )
!        - ( beta - 1 ) / ( 1 - y )
!        + ( alpha + beta + 1 ) / ( 2 * ( 1 - y )**2 )
!        + 4 * y * ( ( Bbar / b ) - gamma )
!        - 0.5 * ( alpha - beta + 3 )
!        + ln ( RHO * R * T / P0 ) )
!                                                      (Equation 2)
!   where
!
!     y = b * rho / 4,
!     alpha = 11,
!     beta = 133/3,
!     gamma = 7/2,
!     P0 = 0.101325 MegaPascals = 1 atm
!
!   and
!
!     b(T) = b1 * ln(T/T0) + sum(j=0,1,3,5) b(j)*(T0/T)**j  (Equation 3)
!
!     Bbar(T) = sum(j=0,1,2,4) B(j)*(T0/T)**j               (Equation 4).
!
!   where
!
!     T0=647.073 K and the coefficients b(j) and B(j) are
!
!     j    b(j)                         B(j)
!    --    -----------                  ----------
!     0    0.7478629                    1.1278334
!     1   -0.3540782                   -0.5944001
!     2    0                           -5.010996
!     3    0.007159876                  0
!     4    0                            0.63684256
!     5   -0.003528426                  0
!
!  For the derived quantities, the following relations are used:
!
!    Pressure:                  PB      = RHO**2 * dAB/dRHO
!    Density derivative:        DPDRB   = 2*PB/RHO + RHO**2 * d2AB/dRHO2
!    Temperature derivative:    DPDTB   = RHO**2 * d2AB/(dRHO dT)
!    Specific entropy:          SB      = ( UB - AB ) / T
!    Specific internal energy:  UB      = AB + T * SB
!    Specific enthalpy:         HB      = UB + PB / RHO
!    Specific heat capacity
!      at constant volume:      CVB     = - T * d2AB/dT2
!    Specific Gibbs function:   GB      = AB + PB / RHO
!
!  Modified:
!
!    03 February 2002
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Input, real ( kind = 8 ) RHO, the density, in G/CM3.
!
!    Output, real ( kind = 8 ) AB, the base value of the Helmholtz function,
!    in KJ/kg.
!
!    Output, real ( kind = 8 ) CVB, the base value of the isochoric (constant
!    volume) heat capacity, in KJ/(kg degrees Kelvin).
!
!    Output, real ( kind = 8 ) DPDRB, the base value of the partial
!    derivative dP(T,RHO)/dRHO, with T held fixed, in (MegaPascals CM3)/G.
!
!    Output, real ( kind = 8 ) DPDTB, the base value of the partial
!    derivative dP(T,RHO)/dT, with RHO held fixed, in
!    MegaPascals/degrees Kelvin.
!
!    Output, real ( kind = 8 ) GB, the base value of the Gibbs free energy,
!    in KJ/kg.
!
!    Output, real ( kind = 8 ) HB, the base value of enthalpy, in KJ/kg.
!
!    Output, real ( kind = 8 ) PB, the base pressure, in MegaPascals.
!
!    Output, real ( kind = 8 ) SB, the base value of entropy,
!    in KJ/(kg degrees Kelvin).
!
!    Output, real ( kind = 8 ) UB, the base value of internal energy,
!    in KJ/kg.
!
  implicit none

  real ( kind = 8 ) ab
  real ( kind = 8 ), parameter :: alpha = 11.0D+00
  real ( kind = 8 ) b1
  real ( kind = 8 ) b1t
  real ( kind = 8 ) b1tt
  real ( kind = 8 ) b2
  real ( kind = 8 ) b2t
  real ( kind = 8 ) b2tt
  real ( kind = 8 ), parameter :: beta = 44.333333333333D+00
  real ( kind = 8 ) cvb
  real ( kind = 8 ) dpdrb
  real ( kind = 8 ) dpdtb
  real ( kind = 8 ) dz
  real ( kind = 8 ) dz0
  real ( kind = 8 ), parameter :: gamma = 3.5D+00
  real ( kind = 8 ) gascon
  real ( kind = 8 ) gb
  real ( kind = 8 ) hb
  real ( kind = 8 ), parameter :: p_zero = 0.101325D+00
  real ( kind = 8 ) pb
  real ( kind = 8 ) rho
  real ( kind = 8 ) sb
  real ( kind = 8 ) t
  real ( kind = 8 ) ub
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  real ( kind = 8 ) z0
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASE - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASE - Fatal error!'
    write ( *, '(a)' ) '  The input density RHO must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was RHO = ', rho
    stop
  end if
!
!  Compute auxilliary quantities for Equation 2.
!
  call bb ( t, b1, b2, b1t, b2t, b1tt, b2tt )

  y = 0.25D+00 * b1 * rho

  x = 1.0D+00 - y
!
!  Evaluate Equation 2.
!
  ab =   - log ( 1.0D+00 - y ) &
         - ( beta - 1.0D+00 ) / ( 1.0D+00 - y ) &
         + ( alpha + beta + 1.0D+00 ) / ( 2.0D+00 * ( 1.0D+00 - y )**2 ) &
         + 4.0D+00 * y * ( ( b2 / b1 ) - gamma ) &
         - 0.5D+00 * ( alpha - beta + 3.0D+00 ) &
         + log ( rho * gascon() * t / p_zero )
!
!  Determine quantities defined in terms of AB.
!
  pb = ( 1.0D+00 + alpha * y + beta * y**2 ) / ( 1.0D+00 - y )**3 &
    + 4.0D+00 * y * ( b2 / b1 - gamma )

  z0 = ( 1.0D+00 + alpha * y + beta * y**2 ) / ( 1.0D+00 - y )**3

  z = z0 + 4.0D+00 * y * ( b2 / b1 - gamma )

  dz0 = ( alpha + 2.0D+00 * beta * y ) / ( 1.0D+00 - y )**3 &
    + 3.0D+00 * ( 1.0D+00 + alpha * y + beta * y**2 ) / ( 1.0D+00 - y )**4

  dz = dz0 + 4.0D+00 * ( b2 / b1 - gamma )

  gb = ab + pb

  ub = - t * b1t * ( pb - 1.0D+00 - rho * b2 ) / b1 - rho * t * b2t

  hb = pb + ub
!
!  An incorrect version of this equation began:
!
!    cvb = 2.0D+00 * ub + ( pb - 1.0D+00 ) &
!
!  and caused me no end of trouble.  My fault, JVB, 03 February 2002
!
  cvb = 2.0D+00 * ub + ( z0 - 1.0D+00 ) &
    * ( ( t * b1t / b1 )**2 - t**2 * b1tt / b1 ) &
    - rho * t**2 * ( b2tt - gamma * b1tt ) - ( t * b1t / b1 )**2 * y * dz0

  dpdtb = pb / t + rho * ( 0.25D+00 * ( dz0 + 4.0D+00 * ( b2 / b1 - gamma ) ) &
    * b1t + b2t - b2 / b1 * b1t )

  sb = ub - ab

  dpdrb = pb + y * ( dz0 + 4.0D+00 * ( b2 / b1 - gamma ) )
!
!  Assign dimensions.
!
  ab =    gascon() * t       * ab
  cvb =   gascon()           * cvb
  dpdrb = gascon() * t       * dpdrb
  dpdtb = gascon() * t * rho * dpdtb
  gb =    gascon() * t       * gb
  hb =    gascon() * t       * hb
  pb =    gascon() * t * rho * pb
  sb =    gascon()           * sb
  ub =    gascon() * t       * ub

  return
end
subroutine bb ( t, b1, b2, b1t, b2t, b1tt, b2tt )

!*****************************************************************************80
!
!! BB calculates the B's of equations 3 and 4.
!
!  Discussion:
!
!    Here
!
!      b(T) = b1 * ln(T/T0) + sum(j=0,1,3,5) b(j)*(T0/T)**j  (Equation 3)
!
!      Bbar(T) = sum(j=0,1,2,4) B(j)*(T0/T)**j               (Equation 4).
!
!    where
!
!      T0 = 647.073 K
!
!    and the coefficients b(j) and B(j) are
!
!      j    b(j)                         B(j)
!     --    -----------                  ----------
!      0    0.7478629                    1.1278334
!      1   -0.3540782                   -0.5944001
!      2    0                           -5.010996
!      3    0.007159876                  0
!      4    0                            0.63684256
!      5   -0.003528426                  0
!
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Output, real ( kind = 8 ) B1, the coefficient b from equation 3,
!    in CM3/G.
!
!    Output, real ( kind = 8 ) B2, the coefficient Bbar from equation 4,
!    in CM3/G.
!
!    Output, real ( kind = 8 ) B1T, the derivative dB1/dT,
!    in (CM3)/(G Degrees Kelvin).
!
!    Output, real ( kind = 8 ) B2T, the derivative dB2/dT,
!    in (CM3)/(G Degrees Kelvin).
!
!    Output, real ( kind = 8 ) B1TT, the second derivative of B1 with
!    respect to T, in (CM3)/(G (Degrees Kelvin)**2 ).
!
!    Output, real ( kind = 8 ) B2TT, the second derivative of B2 with
!    respect to T, in (CM3)/(G (Degrees Kelvin)**2 ).
!
  implicit none

  real ( kind = 8 ) b1
  real ( kind = 8 ) b1t
  real ( kind = 8 ) b1tt
  real ( kind = 8 ) b2
  real ( kind = 8 ) b2t
  real ( kind = 8 ) b2tt
  real ( kind = 8 ), parameter, dimension ( 10 ) :: bp = (/ &
    0.7478629D+00,   -0.3540782D+00,    0.0D+00,           0.0D+00, &
    0.007159876D+00,  0.0D+00,         -0.003528426D+00,   0.0D+00, &
    0.0D+00,          0.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 10 ) :: bq = (/ &
    1.1278334D+00,    0.0D+00,         -0.5944001D+00,   -5.010996D+00, &
    0.0D+00,          0.63684256D+00,   0.0D+00,          0.0D+00, &
    0.0D+00,          0.0D+00 /)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: t_ref = 647.073D+00
  real ( kind = 8 ) v(10)
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BB - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Set V(I) = ( T_REF / T )**(I-1).
!
  v(1) = 1.0D+00
  do i = 2, 10
    v(i) = v(i-1) * t_ref / t
  end do
!
!  Set B1, B1T, B1TT.
!
  b1 = bp(1) + bp(2) * log ( 1.0D+00 / v(2) )
  b1t = bp(2) * v(2) / t_ref
  b1tt = 0.0D+00
  do i = 3, 10
    b1 = b1 + bp(i) * v(i-1)
    b1t = b1t - dble ( i - 2 ) * bp(i) * v(i-1) / t
    b1tt = b1tt + bp(i) * dble ( i - 2 )**2 * v(i-1) / t**2
  end do

  b1tt = b1tt -  ( b1t / t )
!
!  Set B2, B2T, B2TT.
!
  b2 = bq(1)
  b2t = 0.0D+00
  b2tt = 0.0D+00
  do i = 3, 10
    b2 = b2 + bq(i) * v(i-1)
    b2t = b2t - dble ( i - 2 ) * bq(i) * v(i-1) / t
    b2tt = b2tt + bq(i) * dble ( i - 2 )**2 * v(i-1) / t**2
  end do

  b2tt = b2tt - ( b2t / t )

  return
end
subroutine corr ( t, p, p_consistent, rhol, rhov, delg )

!*****************************************************************************80
!
!! CORR evaluates an adjustment to the Gibbs function.
!
!  Discussion:
!
!    CORR is given T and P at or near the vapor pressure and evaluates
!    the corresponding liquid and vapor densities, and the residual
!    function DELG = (GL-GV)/(R*T) where GL and GV are the Gibbs functions
!    for the liquid and vapor phases, respectively.
!
!    These quantities are used to calculate a correction to the vapor
!    pressure or the vapor temperature.
!
!    The states corresponding to the coexisting phases of liquid
!    and vapor for the temperature range from the triple point
!    to within 0.5 C of the critical point 0.01 <= t <= tk-0.5 C
!    have been determined in exact accord with the Gibbs condition
!    of phase equilibrium: DELG = G(g)-G(l) = 0, P, t constant,
!    where G(g) and G(l) are the values of the Gibbs function
!    for saturated gas and liquid respectively.
!
!    For the region (tk-t)<=0.5 C, an exact solution for the
!    Helmholtz function yields values of density for the saturated
!    liquid that are shifted to lower values.  Also, the isotherms
!    in the pressure-density plane and the Gibbs function-density
!    plane are nearly flat, so that it is difficult to obtain
!    solutions.  As an alternative to exact solution, the power
!    law equation is used to define states:
!
!      rho(gas) = 0.322 - 0.657 * (1 - T/647.126)**0.325   (g/cm3).
!      rho(liq) = 0.322 + 0.657 * (1 - T/647.126)**0.325   (g/cm3).
!
!    In a poor instance of programming, the input pressure was
!    originally overwritten on output by a value consistent with
!    the computed densities.  This causes no end of misunderstandings,
!    since other routines expect the value of pressure to be input
!    only.  The code is now revised so that there is an input P
!    and an output P.  In a huff, JVB 05 February 2002.
!
!  Modified:
!
!    05 February 2002
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the vapor temperature, in degrees Kelvin.
!
!    Input, real ( kind = 8 ) P, the vapor pressure, in MegaPascals.
!
!    Output, real ( kind = 8 ) P_CONSISTENT, the vapor pressure, in MegaPascals,
!    consistent with RHOL and RHOV.  This is equal to the input value of
!    P unless 646.3 <= T.
!
!    Input/output, real ( kind = 8 ) RHOL, the liquid density, in G/CM3.
!    On input, if RHOL is positive, it is used as an initial
!    estimate for the iteration.
!
!    Input/output, real ( kind = 8 ) RHOV, the vapor density, in G/CM3.
!    On input, if RHOV is positive, it is used as an initial
!    estimate for the iteration.
!
!    Output, real ( kind = 8 ) DELG, the residual function (GL-GV)/(R*T),
!    where GL is the liquid Gibbs function, GV the vapor Gibbs function,
!    dimensionless.  If 646.3 < T, DELG is 0.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ab
  real ( kind = 8 ) ar
  real ( kind = 8 ) cd
  real ( kind = 8 ) cjth
  real ( kind = 8 ) cjtt
  real ( kind = 8 ) cp
  real ( kind = 8 ) cv
  real ( kind = 8 ) cvb
  real ( kind = 8 ) cvr
  logical, parameter :: debug = .false.
  real ( kind = 8 ) delg
  real ( kind = 8 ) dpdr
  real ( kind = 8 ) dpdrb
  real ( kind = 8 ) dpdrr
  real ( kind = 8 ) dpdt
  real ( kind = 8 ) dpdtb
  real ( kind = 8 ) dpdtr
  real ( kind = 8 ) g
  real ( kind = 8 ) gascon
  real ( kind = 8 ) gb
  real ( kind = 8 ) gl
  real ( kind = 8 ) gr
  real ( kind = 8 ) gv
  real ( kind = 8 ) h
  real ( kind = 8 ) hb
  real ( kind = 8 ) hr
  real ( kind = 8 ) p
  real ( kind = 8 ) p_consistent
  real ( kind = 8 ), parameter :: p_crit = 22.055D+00
  real ( kind = 8 ) pb
  real ( kind = 8 ) pr
  real ( kind = 8 ) rho
  real ( kind = 8 ) rhol
  real ( kind = 8 ), parameter :: rho_min = 1.0D-08
  real ( kind = 8 ) rhov
  real ( kind = 8 ) rho_start
  real ( kind = 8 ) s
  real ( kind = 8 ) sb
  real ( kind = 8 ) sr
  real ( kind = 8 ) t
  real ( kind = 8 ) tau
  real ( kind = 8 ), parameter :: t_crit = 647.1260000001D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) ub
  real ( kind = 8 ) ur

  p_consistent = p
!
!  Initialize output quantities.
!
  delg = 0.0D+00
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CORR - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  The input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative pressures.
!
  if ( p <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CORR - Fatal error!'
    write ( *, '(a)' ) '  The input pressure P must be positive.'
    write ( *, '(a,g14.6)' ) '  The input value was P = ', p
    stop
  end if

  if ( t <= 646.3D+00 ) then

    if ( rhol <= 0.0D+00 ) then
      rho_start = 1.11D+00 - 0.0004D+00 * t
    else
      rho_start = rhol
    end if

    call dense ( p_consistent, t, rho_start, rho, dpdr )

    call therm ( t, rho, a, cjth, cjtt, cd, cv, dpdr, dpdt, g, h, &
      p_consistent, s, u )

    rhol = rho
    gl = g

    if ( rhov <= 0.0D+00 ) then
      rho_start = p_consistent / ( gascon() * t )
    else
      rho_start = rhov
    end if

    call dense ( p_consistent, t, rho_start, rho, dpdr )

    rho = max ( rho, rho_min )

    call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, &
      p_consistent, s, u )

    rhov = rho
    gv = g
    delg = ( gl - gv ) / ( gascon() * t )

    p_consistent = p

    if ( debug ) then
      write ( *, '(a,g14.6)' ) '  CORR - RHOL = ', rhol
      write ( *, '(a,g14.6)' ) '  RHOV = ', rhov
    end if

  else if ( t <= t_crit ) then

    if ( debug ) then
      write ( *, '(a)' ) 'CORR - Twilight zone'
    end if

    delg = 0.0D+00
    tau = 0.657128D+00 * ( 1.0D+00 - t / t_crit )**0.325D+00
    rhol = 0.322D+00 + tau
    rhov = 0.322D+00 - tau
    rho = rhov
    call base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )
    call resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )
    p_consistent = pb + pr

  else

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CORR - Weirdo zone'
    end if

    rhol = 0.322D+00
    rhov = 0.322D+00
    p_consistent = p_crit
    delg = 0.0D+00

  end if

  return
end
subroutine cp_values ( n, tc, p, cp )

!*****************************************************************************80
!
!! CP_VALUES returns some values of the specific heat at constant pressure.
!
!  Modified:
!
!    04 February 2002
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
!    TJ270.H3, pages 229-237.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) CP, the specific heat at constant pressure,
!    in KJ/(kg K).
!
  implicit none

  integer ( kind = 4 ), parameter :: n_data = 24

  real ( kind = 8 ) cp
  real ( kind = 8 ), save, dimension ( n_data ) :: cpvec = (/ &
    4.228D+00, 2.042D+00, 1.975D+00, 2.013D+00, 2.040D+00, &
    2.070D+00, 2.135D+00, 2.203D+00, 2.378D+00, 2.541D+00, &
    2.792D+00, 2.931D+00, 4.226D+00, 4.223D+00, 4.202D+00, &
    4.177D+00, 4.130D+00, 4.089D+00, 4.053D+00, 4.021D+00, &
    3.909D+00, 3.844D+00, 3.786D+00, 2.89D+00 /)
  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_data ) :: pvec = (/ &
       1.0D+00,    1.0D+00,    1.0D+00,    1.0D+00,   1.0D+00, &
       1.0D+00,    1.0D+00,    1.0D+00,    1.0D+00,   1.0D+00, &
       1.0D+00,    1.0D+00,    5.0D+00,   10.0D+00,  50.0D+00, &
     100.0D+00,  200.0D+00,  300.0D+00,  400.0D+00, 500.0D+00, &
    1000.0D+00, 1500.0D+00, 2000.0D+00, 5000.0D+00  /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_data ) :: tcvec = (/ &
       0.0D+00,  100.0D+00,  200.0D+00,  300.0D+00,  350.0D+00, &
     400.0D+00,  500.0D+00,  600.0D+00,  850.0D+00, 1100.0D+00, &
    1600.0D+00, 2000.0D+00,    0.0D+00,    0.0D+00,    0.0D+00, &
       0.0D+00,    0.0D+00,    0.0D+00,    0.0D+00,    0.0D+00, &
       0.0D+00,    0.0D+00,    0.0D+00,    0.0D+00 /)

  if ( n < 0 ) then
    n = 0
  end if

  if ( n_data <= n ) then
    n = 0
    tc = 0.0D+00
    p = 0.0D+00
    cp = 0.0D+00
  else
    n = n + 1
    tc = tcvec(n)
    p = pvec(n)
    cp = cpvec(n)
  end if

  return
end
subroutine dense ( p, t, rho_start, rho, dpdr )

!*****************************************************************************80
!
!! DENSE computes the density for a given pressure and temperature.
!
!  Discussion:
!
!    The use of the variable RHO_START for two opposing purposes is
!    poor practice and will be corrected one of these days.  Meanwhile,
!    the algorithm's behavior, particularly in the two-phase region,
!    is very suspect.
!
!  Modified:
!
!    19 November 1998
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the pressure, in MegaPascals.
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Input, real ( kind = 8 ) RHO_START, an initial guess for the density,
!    in G/CM3.  The value of RHO_START also signals whether a vapor or liquid
!    calculation is to be done.  If DPDR is computed negative, then for
!    0.2967 <= RHO_START, liquid is assumed, otherwise gas.
!
!    Output, real ( kind = 8 ) RHO, the density for the given
!    pressure and temperature, in G/CM3.
!
!    Output, real ( kind = 8 ) DPDR, the partial derivative
!    dP(T,RHO)/dRHO, with T held fixed, in (MegaPascals CM3)/G.
!
  implicit none

  real ( kind = 8 ) ab
  real ( kind = 8 ) ar
  real ( kind = 8 ) cvb
  real ( kind = 8 ) cvr
  real ( kind = 8 ) dp
  real ( kind = 8 ) dpdr
  real ( kind = 8 ) dpdrb
  real ( kind = 8 ) dpdrr
  real ( kind = 8 ) dpdtb
  real ( kind = 8 ) dpdtr
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) errtol
  real ( kind = 8 ) gb
  real ( kind = 8 ) gr
  real ( kind = 8 ) hb
  real ( kind = 8 ) hr
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 50
  real ( kind = 8 ) p
  real ( kind = 8 ) pb
  real ( kind = 8 ) pp
  real ( kind = 8 ) pr
  real ( kind = 8 ) rho
  real ( kind = 8 ), parameter :: rho_max = 1.9D+00
  real ( kind = 8 ), parameter :: rho_min = 1.0D-08
  real ( kind = 8 ) rho_start
  real ( kind = 8 ) sb
  real ( kind = 8 ) sr
  real ( kind = 8 ) t
  real ( kind = 8 ) ub
  real ( kind = 8 ) ur
  real ( kind = 8 ) x

  errtol = sqrt ( epsilon ( errtol ) )
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DENSE - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative pressures.
!
  if ( p <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DENSE - Fatal error!'
    write ( *, '(a)' ) '  The input pressure P must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was P = ', p
    stop
  end if

  rho = rho_start
  rho = max ( rho, rho_min )
  rho = min ( rho, rho_max )

  do it = 1, it_max

    call resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )

    call base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )

    pp = pb + pr
    dpdr = dpdrb + dpdrr
!
!  Check for negative DP/DRho, which characterizes the two-phase region.
!
    if ( dpdr <= 0.0D+00 ) then

      if ( 0.2967D+00 <= rho_start ) then
        rho = rho * 1.02D+00
      else
        rho = rho * 0.98D+00
      end if

      if ( it <= 10 ) then
        cycle
      end if

    end if

    dpdx = 1.1D+00 * dpdr
    dpdx = max ( dpdx, 0.01D+00 )

    dp = abs ( 1.0D+00 - pp / p )

    if ( dp <= errtol .or. &
       ( 0.3D+00 < rho .and. dp <= errtol ) .or. &
       ( 0.7D+00 < rho .and. dp <= 10.0D+00 * errtol ) ) then
      return
    end if

    x = ( p - pp ) / dpdx
    if ( 0.1D+00 < abs ( x ) ) then
      x = x * 0.1D+00 / abs ( x )
    end if

    rho = rho + x

    rho = max ( rho, rho_min )
    rho = min ( rho, rho_max )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DENSE - Warning!'
  write ( *, '(a)' ) '  The iteration did not converge.'
  write ( *, '(a,i6)' ) '  Number of iterations was ', it_max
  write ( *, '(a,g14.6)' ) '  Last iterate was ', rho

  return
end
subroutine dielectric ( t, rho, eps )

!*****************************************************************************80
!
!! DIELECTRIC returns the static dielectric constant.
!
!  Discussion:
!
!    According to the IAPS, the equation used is valid in the range
!
!      273.15 degrees Kelvin <= T <= 823.15 degrees K
!      0 MegaPascals         <= P <= 500 MegaPascals.
!
!  Modified:
!
!    02 February 2002
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
!    TJ270.H3, page 266.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Input, real ( kind = 8 ) RHO, the density, in G/CM3.
!
!    Output, real ( kind = 8 ) EPS, the dielectric constant, dimensionless.
!
  implicit none

  integer ( kind = 4 ), parameter :: npol_c = 4

  real ( kind = 8 ), parameter, dimension ( 10 ) :: a = (/ &
    7.62571D+00,  244.003D+00,  -140.569D+00,   27.7841D+00, -96.2805D+00, &
   41.7909D+00,   -10.2099D+00,  -45.2059D+00,  84.6395D+00, -35.8644D+00 /)
  real ( kind = 8 ) c(0:npol_c)
  real ( kind = 8 ) eps
  real ( kind = 8 ) rho
  real ( kind = 8 ) t
  real ( kind = 8 ) t_copy
  real ( kind = 8 ), parameter :: t_max = 823.15D+00
  real ( kind = 8 ), parameter :: t_min = 273.15D+00
  real ( kind = 8 ), parameter :: t_ref = 298.15D+00
  real ( kind = 8 ) t_star
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIELECTRIC - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIELECTRIC - Fatal error!'
    write ( *, '(a)' ) '  The input density RHO must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was RHO = ', rho
    stop
  end if

  t_copy = t
  t_copy = min ( t_copy, t_max )
  t_copy = max ( t_copy, t_min )

  t_star = t_copy / t_ref

  c(0) = 1.0D+00
  c(1) = a(1) / t_star
  c(2) = ( a(2) / t_star ) + a(3) + a(4) * t_star
  c(3) = ( a(5) / t_star ) + a(6) * t_star + a(7) * t_star**2
  c(4) = ( a(8) / t_star**2 ) + ( a(9) / t_star ) + a(10)

  call dpoly_val_horner ( npol_c, c, rho, eps )

  return
end
subroutine dielectric_values ( n, tc, p, eps )

!*****************************************************************************80
!
!! DIELECTRIC_VALUES returns some values of the static dielectric constant.
!
!  Modified:
!
!    03 February 2002
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
!    TJ270.H3, page 266.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) EPS, the dielectric constant, dimensionless.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_data = 15

  real ( kind = 8 ) eps
  real ( kind = 8 ), save, dimension ( n_data ) :: epsvec = (/ &
    88.29D+00, 90.07D+00, 92.02D+00, 95.14D+00, 100.77D+00, &
    78.85D+00, 70.27D+00, 62.60D+00, 55.78D+00,  44.31D+00, &
    35.11D+00, 20.40D+00,  1.17D+00,  1.11D+00,   1.08D+00 /)
  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_data ) :: pvec = (/ &
    100.0D+00, 500.0D+00, 1000.0D+00, 2000.0D+00, 5000.0D+00, &
    100.0D+00, 100.0D+00, 100.0D+00, 100.0D+00, 100.0D+00, &
    100.0D+00, 100.0D+00, 100.0D+00, 100.0D+00, 100.0D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_data ) :: tcvec = (/ &
      0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
     25.0D+00,  50.0D+00,  75.0D+00, 100.0D+00, 150.0D+00, &
    200.0D+00, 300.0D+00, 400.0D+00, 500.0D+00, 600.0D+00 /)

  if ( n < 0 ) then
    n = 0
  end if

  if ( n_data <= n ) then
    n = 0
    tc = 0.0D+00
    p = 0.0D+00
    eps = 0.0D+00
  else
    n = n + 1
    tc = tcvec(n)
    p = pvec(n)
    eps = epsvec(n)
  end if

  return
end
subroutine dpoly_val_horner ( n, c, x, cx )

!*****************************************************************************80
!
!! DPOLY_VAL_HORNER evaluates a polynomial using Horner's method.
!
!  Modified:
!
!    08 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of C.
!
!    Input, real ( kind = 8 ) C(0:N), the polynomial coefficients.
!    C(I) is the coefficient of X**I.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) CX, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n)
  real ( kind = 8 ) cx
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  cx = c(n)
  do i = n - 1, 0, -1
    cx = cx * x + c(i)
  end do

  return
end
function gascon ( )

!*****************************************************************************80
!
!! GASCON returns the value of the specific gas constant.
!
!  Note:
!
!    The specific gas constant R is related to the universal gas
!    constant R-bar = 8.31441 J/(mol degrees Kelvin) by the molar mass
!    M = 18.0152 g/mol:
!
!      R = R-bar / M.
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) GASCON, the value of the specific gas
!    constant, in J/(g degrees Kelvin).
!
  implicit none

  real ( kind = 8 ) gascon

  gascon = 0.461522D+00

  return
end
subroutine ideal ( t, ai, cpi, cvi, gi, hi, si, ui )

!*****************************************************************************80
!
!! IDEAL computes ideal gas thermodynamic properties of water.
!
!  Discussion:
!
!    Values for thermodynamic properties of water in the ideal
!    gas state were reported by Woolley.  The formula for the ideal gas
!    term of the Helmholtz function approximates a term by term summation of
!    contributions from each of the rotation and vibration states.
!    The formula, equation #6 in the reference, is:
!
!    A(ideal)(T) = -R * T * ( 1 + ( C(1)/Tr + C(2) ) * ln(Tr)
!      + Sum ( 3 <= I <= 18) C(I) * Tr**(I-6)
!
!    where Tr=T/100 K.  The C(i) are tabulated coefficients.  Equation
!    6 can be used for temperatures below 3000 K, and is accurate to
!    within the tolerance of the gas constant for 50<=T<=2000 K.
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
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Output, real ( kind = 8 ) AI, the Helmholtz function, in KJ/kg.
!
!    Output, real ( kind = 8 ) CPI, the heat capacity at constant pressure,
!    in KJ/(kg degrees Kelvin).
!
!    Output, real ( kind = 8 ) CVI, the heat capacity at constant volume,
!    in KJ/(kg degrees Kelvin).
!
!    Output, real ( kind = 8 ) GI, the Gibbs free energy, in KJ/kg.
!
!    Output, real ( kind = 8 ) HI, the enthalpy, in KJ/kg.
!
!    Output, real ( kind = 8 ) SI, the entropy, in KJ/(kg degrees Kelvin).
!
!    Output, real ( kind = 8 ) UI, the internal energy, in KJ/kg.
!
  implicit none

  real ( kind = 8 ) ai
  real ( kind = 8 ), parameter, dimension ( 18 ) :: c = (/ &
   19.730271018D+00,      20.9662681977D+00,     -0.483429455355D+00, &
    6.05743189245D+00,    22.56023885D+00,       -9.87532442D+00, &
   -4.3135538513D+00,      0.458155781D+00,      -0.047754901883D+00, &
    0.0041238460633D+00,  -0.00027929052852D+00,  0.14481695261D-04, &
   -0.56473658748D-06,     0.16200446D-07,       -0.3303822796D-09, &
    0.451916067368D-11,   -0.370734122708D-13,    0.137546068238D-15 /)
  real ( kind = 8 ) cpi
  real ( kind = 8 ) cvi
  real ( kind = 8 ) gascon
  real ( kind = 8 ) gi
  real ( kind = 8 ) hi
  integer ( kind = 4 ) i
  real ( kind = 8 ) si
  real ( kind = 8 ) t
  real ( kind = 8 ) tt
  real ( kind = 8 ) ui
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IDEAL - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if

  tt = t / 100.0D+00

  gi = - ( c(1) / tt + c(2) ) * log ( tt )
  do i = 3, 18
    gi = gi - c(i) * tt**(i-6)
  end do

  hi = c(2) + c(1) * ( 1.0D+00 - log ( tt ) ) / tt
  do i = 3, 18
    hi = hi + dble ( i - 6 ) * c(i) * tt**(i-6)
  end do

  cpi = c(2) - c(1) / tt
  do i = 3, 18
    cpi = cpi + dble ( ( i - 6 ) * ( i - 5 ) ) * c(i) * tt**(i-6)
  end do

  ai = gi - 1.0D+00
  ui = hi - 1.0D+00
  cvi = cpi - 1.0D+00
  si = hi - gi
!
!  Assign dimensions.
!
  ai =  gascon() * t * ai
  cpi = gascon()     * cpi
  cvi = gascon()     * cvi
  gi =  gascon() * t * gi
  hi =  gascon() * t * hi
  si =  gascon()     * si
  ui =  gascon() * t * ui

  return
end
subroutine prandtl ( t, p, pr )

!*****************************************************************************80
!
!! PRANDTL computes the Prandtl number.
!
!  Discussion:
!
!    This routine was NOT working properly for large pressures,
!    because the routine CORR was changing the input value of P.
!
!  Modified:
!
!    17 November 1998
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Input, real ( kind = 8 ) P, the pressure, in MegaPascals.
!
!    Output, real ( kind = 8 ) PR, the Prandtl number, dimensionless.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cjth
  real ( kind = 8 ) cjtt
  real ( kind = 8 ) cp
  real ( kind = 8 ) cv
  logical, parameter :: debug = .false.
  real ( kind = 8 ) dpdr
  real ( kind = 8 ) dpdt
  real ( kind = 8 ) eta
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  real ( kind = 8 ) lambda
  real ( kind = 8 ) p
  real ( kind = 8 ) pr
  real ( kind = 8 ) rho
  real ( kind = 8 ) rhol
  real ( kind = 8 ) rho_start
  real ( kind = 8 ) rhov
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) t2
  real ( kind = 8 ) u
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRANDTL - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  The input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative pressures.
!
  if ( p <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRANDTL - Fatal error!'
    write ( *, '(a)' ) '  The input pressure P must be positive.'
    write ( *, '(a,g14.6)' ) '  The input value was P = ', p
    stop
  end if
!
!  Compute the density.
!
  if ( debug ) then
    write ( *, * ) 'PRANDTL - Call TSAT, with P = ', p
  end if

  rhol = 0.0D+00
  rhov = 0.0D+00

  call tsat ( p, t2, rhol, rhov )

  if ( debug ) then
    write ( *, * ) 'PRANDTL - T2 = ', t2
  end if

  if ( t < t2 ) then
    rho_start = 1.9D+00
  else
    rho_start = 0.01D+00
  end if

  call dense ( p, t, rho_start, rho, dpdr )

  if ( debug ) then
    write ( *, * ) 'PRANDTL - RHO = ', rho
  end if
!
!  Now from T and RHO, compute CP, ETA and LAMBDA.
!
  call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

  call viscosity ( t, rho, eta )

  call thercon ( t, rho, lambda )

  if ( debug ) then
    write ( *, '(7f10.4)' ) t, p, rho, eta, cp, lambda, pr
  end if

  pr = eta * cp / lambda

  return
end
subroutine prandtl_values ( n, tc, p, pr )

!*****************************************************************************80
!
!! PRANDTL_VALUES returns some values of the Prandtl number for testing.
!
!  Modified:
!
!    04 February 2002
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
!    TJ270.H3, page 265.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) PR, the Prandtl number, dimensionless.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_data = 35

  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ) pr
  real ( kind = 8 ), save, dimension ( n_data ) :: prvec = (/ &
    13.50D+00, 13.48D+00, 13.46D+00, 13.39D+00, 13.27D+00, &
    13.15D+00, 13.04D+00, 12.93D+00, 12.83D+00, 12.73D+00, &
    12.63D+00, 12.53D+00, 12.43D+00, 12.34D+00, 12.25D+00, &
    12.08D+00, 11.92D+00, 11.77D+00, 11.62D+00, 11.48D+00, &
    11.36D+00, 11.23D+00, 11.12D+00, 10.91D+00, 10.72D+00, &
    10.55D+00,  6.137D+00, 3.555D+00, 2.378D+00, 1.000D+00, &
     0.974D+00, 0.960D+00, 0.924D+00, 0.899D+00, 0.882D+00 /)
  real ( kind = 8 ), save, dimension ( n_data ) :: pvec = (/ &
       1.0D+00,    5.0D+00,   10.0D+00,   25.0D+00,   50.0D+00, &
      75.0D+00,  100.0D+00,  125.0D+00,  150.0D+00,  175.0D+00, &
     200.0D+00,  225.0D+00,  250.0D+00,  275.0D+00,  300.0D+00, &
     350.0D+00,  400.0D+00,  450.0D+00,  500.0D+00,  550.0D+00, &
     600.0D+00,  650.0D+00,  700.0D+00,  800.0D+00,  900.0D+00, &
    1000.0D+00,    1.0D+00,    1.0D+00,    1.0D+00,    1.0D+00, &
       1.0D+00,    1.0D+00,    1.0D+00,    1.0D+00,    1.0D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_data ) :: tcvec = (/ &
    0.0D+00,     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
    0.0D+00,     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
    0.0D+00,     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
    0.0D+00,     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
    0.0D+00,     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
    0.0D+00,    25.0D+00,  50.0D+00,  75.0D+00, 100.0D+00, &
    150.0D+00, 200.0D+00, 400.0D+00, 600.0D+00, 800.0D+00 /)

  if ( n < 0 ) then
    n = 0
  end if

  if ( n_data <= n ) then
    n = 0
    tc = 0.0D+00
    p = 0.0D+00
    pr = 0.0D+00
  else
    n = n + 1
    tc = tcvec(n)
    p = pvec(n)
    pr = prvec(n)
  end if

  return
end
subroutine psat ( t, p, rhol, rhov )

!*****************************************************************************80
!
!! PSAT calculates the vapor pressure, and the liquid and vapor densities.
!
!  Discussion:
!
!    These quantities correspond to the input temperature T, corrected
!    so that the Gibbs functions for liquid and vapor phase are
!    equal to within a tolerance.
!
!  Modified:
!
!    04 February 2002
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the vapor temperature, in degrees Kelvin.
!
!    Output, real ( kind = 8 ) P, the vapor pressure, in MegaPascals.
!
!    Output, real ( kind = 8 ) RHOL, the liquid density, in G/CM3.
!
!    Output, real ( kind = 8 ) RHOV, the vapor density, in G/CM3.
!
  implicit none

  real ( kind = 8 ) bot
  real ( kind = 8 ) delg
  real ( kind = 8 ) dp
  real ( kind = 8 ) errtol
  real ( kind = 8 ) gascon
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  real ( kind = 8 ) p
  real ( kind = 8 ) p_consistent
  real ( kind = 8 ) p_old
  real ( kind = 8 ) rhol
  real ( kind = 8 ) rhov
  real ( kind = 8 ) t
!
!  Ensure that output quantities are initialized,, obliterating any
!  input values.
!
  p = 0.0D+00
  rhol = 0.0D+00
  rhov = 0.0D+00
!
!  Set the error tolerance.
!
  errtol = 100.0D+00 * sqrt ( epsilon ( errtol ) )
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PSAT - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Get an estimate for the saturation pressure.
!
  call psat_est ( t, p )

  dp = 0.0D+00

  do it = 1, it_max

    call corr ( t, p, p_consistent, rhol, rhov, delg )

    bot = ( rhol - rhov ) / ( rhol * rhov )

    if ( abs ( bot ) < errtol ) then
      write ( *, * ) 'PSAT - Warning, what is this?'
      bot = sign ( errtol, bot )
    end if

    dp = delg * gascon() * t / bot

    p_old = p
    p = p + dp

    if ( abs ( dp ) <= errtol * ( abs ( p ) + 1.0D+00 ) ) then
      return
    end if

    if ( p <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PSAT - Warning!'
      write ( *, '(a)' ) '  The iterates have become nonpositive.'
      write ( *, '(a,i6)' ) '  Iteration number = ', it
      write ( *, '(a,g14.6)' ) '  Last iterate was ', p
      write ( *, '(a,g14.6)' ) '  Previous iterate was ', p_old
      write ( *, '(a,g14.6)' ) '  Last correction was ', dp
      write ( *, '(a)' ) '  Trying to recover...'
      p = 0.5D+00 * p_old
    end if

    if ( p <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PSAT - Fatal error!'
      write ( *, '(a)' ) '  The iterates have become nonpositive.'
      write ( *, '(a,i6)' ) '  Iteration number = ', it
      write ( *, '(a,g14.6)' ) '  Last iterate was ', p
      write ( *, '(a,g14.6)' ) '  Previous iterate was ', p_old
      write ( *, '(a,g14.6)' ) '  Last correction was ', dp
      stop
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PSAT - Warning!'
  write ( *, '(a)' ) '  The iteration did not converge.'
  write ( *, '(a,i6)' ) '  The number of iterations was ', it_max
  write ( *, '(a,g14.6)' ) '  Convergence tolerance was ', errtol
  write ( *, '(a,g14.6)' ) '  Last iterate was ', p
  write ( *, '(a,g14.6)' ) '  Last correction was ', dp

  return
end
subroutine psat_est ( t, p )

!*****************************************************************************80
!
!! PSAT_EST makes a rough estimate of the vapor pressure.
!
!  Discussion:
!
!    The calculation agrees with tabulated data to within
!    0.02% for temperature to within a degree or so of the critical
!    temperature.  The approximate vapor pressure can be refined
!    by imposing the condition that the Gibbs functions of the vapor
!    and liquid phases be equal.
!
!  Modified:
!
!    21 November 1998
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Output, real ( kind = 8 ) P, the vapor pressure, in MegaPascals.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
    -7.8889166D+00,   2.5514255D+00,   -6.716169D+00, 33.239495D+00, &
    -105.38479D+00,   174.35319D+00,  -148.39348D+00, 48.631602D+00 /)
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: t_ref = 647.25D+00
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) z
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PSAT_EST - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if

  if ( t <= 314.0D+00 ) then

    p = 0.1D+00 * exp ( 6.3573118D+00 - 8858.843D+00 / t &
      + 607.56335D+00 * t**( -0.6D+00 ) )

  else

    v = t / t_ref
    w = abs ( 1.0D+00 - v )
    b = 0.0D+00
    do i = 1, 8
      z = i
      b = b + a(i) * w**( ( z + 1.0D+00 ) / 2.0D+00 )
    end do

    q = b / v
    p = 22.093D+00 * exp ( q )

  end if

  return
end
subroutine psat_values ( n, tc, p )

!*****************************************************************************80
!
!! PSAT_VALUES returns some values of the saturation pressure.
!
!  Modified:
!
!    04 February 2002
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
!    TJ270.H3, pages 9-15.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the saturation pressure, in bar.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_data = 12

  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_data ) :: pvec = (/ &
    0.0061173D+00,  0.0065716D+00, 0.0087260D+00,  0.12344D+00,  1.0132D+00, &
    2.3201D+00,     4.7572D+00,   15.537D+00,     39.737D+00,   85.838D+00, &
    165.21D+00,   220.55D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_data ) :: tcvec = (/ &
        0.01D+00,   1.0D+00,    5.0D+00,   50.0D+00,   100.0D+00, &
      125.0D+00,  150.0D+00,  200.0D+00,  250.0D+00,  300.0D+00, &
      350.0D+00,  373.976D+00 /)

  if ( n < 0 ) then
    n = 0
  end if

  if ( n_data <= n ) then
    n = 0
    tc = 0.0D+00
    p = 0.0D+00
  else
    n = n + 1
    tc = tcvec(n)
    p = pvec(n)
  end if

  return
end
subroutine resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )

!*****************************************************************************80
!
!! RESID calculates residual contributions to thermodynamic quantities.
!
!  Discussion:
!
!    The residual function consists of 40 terms.  The first 36 are
!    used in a global least squares fit to experimental data.
!
!    Three terms were added that contribute only in the immediate
!    neighborhood of the critical point
!      (tk-5) <= T <= (tk+5) C
!      0.20   <= rho <= 0.44 g/cm3,
!
!    A single term was added for the region of high pressure and
!    low temperature: T < 75 C, 300 MPa < P.
!
!    Except in these limited regions, the residual function is
!    given by the first 36 terms.  The equation is
!
!      A(residual)(rho,T)=
!        sum(i=1 to 36) (g(i)/k(i)) * (T0/T)**(l(i)) (1-exp(-rho))**(k(i))
!      + sum(i=37 to 40) g(i)*delta(i)**(k(i))
!        * exp(-alpha(i)*delta(i)**(k(i)) - beta(i)*tau(i)**2)
!                                                     (Equation 5)
!
!    where
!
!      g(i) are coefficients determined by fits to data,
!      delta(i) are reduced densities (delta(i)=((rho-rho(i))/rho(i))
!      tau(i) are reduced temperatures (tau(i)=((T-tau(i))/tau(i))
!      rho(i) are specified densities.
!      tau(i) are specified temperatures.
!      The k(i) and l(i) are specified integers.
!
!  Modified:
!
!    22 November 1998
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Input, real ( kind = 8 ) RHO, the density, in G/CM3.
!
!    Output, real ( kind = 8 ) AR, the residual contribution to the
!    Helmholtz function, in KJ/kg.
!
!    Output, real ( kind = 8 ) CVR, the residual contribution to the
!    isochoric (constant volume) heat capacity, in KJ/(kg degrees Kelvin).
!
!    Output, real ( kind = 8 ) DPDRR, the residual contribution to
!    the partial derivative dP(T,RHO)/dRHO, with T held fixed, in
!    (MegaPascals CM3)/G.
!
!    Output, real ( kind = 8 ) DPDTR, the residual contribution to
!    the partial derivative dP(T,RHO)/dT, with RHO held fixed,
!    in MegaPascals/degrees Kelvin.
!
!    Output, real ( kind = 8 ) GR, the residual contribution to the Gibbs
!    function, in KJ/kg.
!
!    Output, real ( kind = 8 ) HR, the residual contribution to the
!    enthalpy, in KJ/kg.
!
!    Output, real ( kind = 8 ) PR, the residual contribution to the pressure,
!    in MegaPascals.
!
!    Output, real ( kind = 8 ) SR, the residual contribution to the entropy,
!    in KJ/(kg degrees Kelvin).
!
!    Output, real ( kind = 8 ) UR, the residual contribution to the
!    internal energy, in KJ/kg.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 4 ) :: aad = (/ &
    34.0D+00, 40.0D+00, 30.0D+00, 1050.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 4 ) :: aat = (/ &
    20000.0D+00, 20000.0D+00, 40000.0D+00, 25.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 4 ) :: adz = (/ &
    0.319D+00, 0.319D+00, 0.319D+00, 1.55D+00 /)
  real ( kind = 8 ) ar
  real ( kind = 8 ) att
  real ( kind = 8 ), parameter, dimension ( 4 ) ::  atz = (/ &
    640.0D+00, 640.0D+00, 641.6D+00, 270.0D+00 /)
  real ( kind = 8 ) cvr
  real ( kind = 8 ) dadt
  real ( kind = 8 ) ddz
  real ( kind = 8 ) del
  real ( kind = 8 ) dex
  real ( kind = 8 ) dfdt
  real ( kind = 8 ) dpdrr
  real ( kind = 8 ) dpdtr
  real ( kind = 8 ) e
  real ( kind = 8 ) errtol
  real ( kind = 8 ) ex0
  real ( kind = 8 ) ex1
  real ( kind = 8 ) ex2
  real ( kind = 8 ) fct
  real ( kind = 8 ), parameter, dimension ( 40 ) :: g = (/ &
    -530.62968529023D+00,  0.22744901424408D+04, 0.78779333020687D+03, &
    -69.830527374994D+00,  0.17863832875422D+05,-0.39514731563338D+05, &
    0.33803884280753D+05, -0.13855050202703D+05,-0.25637436613260D+06, &
    0.48212575981415D+06, -0.34183016969660D+06, 0.12223156417448D+06, &
    0.11797433655832D+07, -0.21734810110373D+07, 0.10829952168620D+07, &
   -0.25441998064049D+06, -0.31377774947767D+07, 0.52911910757704D+07, &
   -0.13802577177877D+07, -0.25109914369001D+06, 0.46561826115608D+07, &
   -0.72752773275387D+07,  0.41774246148294D+06, 0.14016358244614D+07, &
   -0.31555231392127D+07,  0.47929666384584D+07, 0.40912664781209D+06, &
   -0.13626369388386D+07,  0.69625220862664D+06,-0.10834900096447D+07, &
   -0.22722827401688D+06,  0.38365486000660D+06, 0.68833257944332D+04, &
    0.21757245522644D+05, -0.26627944829770D+04,-0.70730418082074D+05, &
   -0.225D+00, -1.68D+00, 0.055D+00, -93.0D+00 /)
  real ( kind = 8 ) gascon
  real ( kind = 8 ) gr
  real ( kind = 8 ) hr
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter, dimension ( 40 ) :: ii = (/ &
    0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6, &
    8,8,8,8,2,2,0,4,2,2,2,4 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter, dimension ( 40 ) :: jj = (/ &
    2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,&
    2,3,5,7,1,4,4,4,0,2,0,0 /)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nc
  real ( kind = 8 ) pr
  real ( kind = 8 ) q10
  real ( kind = 8 ) q20
  real ( kind = 8 ) q2a
  real ( kind = 8 ) q5t
  real ( kind = 8 ) qm
  real ( kind = 8 ) qp
  real ( kind = 8 ) qr(11)
  real ( kind = 8 ) qt(10)
  real ( kind = 8 ) rho
  real ( kind = 8 ) sr
  real ( kind = 8 ), parameter :: s_ref = 7.6180720166752D+00
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: t_ref = 647.073D+00
  real ( kind = 8 ) tau
  real ( kind = 8 ) tx
  real ( kind = 8 ), parameter :: u_ref = - 4328.4549774261D+00
  real ( kind = 8 ) ur
  real ( kind = 8 ) v

  errtol = sqrt ( epsilon ( errtol ) )
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RESID - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RESID - Fatal error!'
    write ( *, '(a)' ) '  The input density RHO must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was RHO = ', rho
    stop
  end if

  nc = 36
  dpdrr = 0.0D+00
  pr = 0.0D+00
  ar = 0.0D+00
  dadt = 0.0D+00
  cvr = 0.0D+00
  dpdtr = 0.0D+00

  ex0 = - rho
  ex0 = max ( ex0, -225.0D+00 )
  ex0 = min ( ex0,  225.0D+00 )
  e = exp ( ex0 )

  q10 = rho * rho * e
  q20 = 1.0D+00 - e

  qr(1) = 0.0D+00
  qr(2) = q10
  do i = 2, 10
    qr(i+1) = qr(i) * q20
  end do

  v = t_ref / t
  qt(1) = t / t_ref
  do i = 2, 10
    qt(i) = qt(i-1) * v
  end do

  do i = 1, nc

    k = ii(i) + 1
    l = jj(i)
    qp = g(i) * qr(k+1) * qt(l+1)
    pr = pr + qp

    dpdrr = dpdrr + ( 2.0D+00 / rho - ( 1.0D+00 - e * dble ( k - 1 ) / &
      ( 1.0D+00 - e ) ) ) * qp

    ar = ar + g(i) * qr(k+2) * qt(l+1) / ( rho**2 * e * dble ( k ) &
      * gascon ( ) * t )

    dfdt = ( 1.0D+00 - e )**k * dble ( 1 - l ) * qt(l+2) / t_ref / dble ( k )

    dadt = dadt + g(i) * dfdt

    dpdtr = dpdtr + g(i) * dfdt * rho**2 * e * dble ( k ) / ( 1.0D+00 - e )

    cvr = cvr + g(i) * dble ( l ) * dfdt / gascon()

  end do

  qp = 0.0D+00
  q2a = 0.0D+00

  do j = 37, 40

    k = ii(j)
    ddz = adz(j-36)
    del = rho / ddz - 1.0D+00

    if ( abs ( del ) < errtol ) then
      del = errtol
    end if

    ex1 = - aad(j-36) * del**k
    ex1 = max ( ex1, -225.0D+00 )
    ex1 = min ( ex1,  225.0D+00 )
    dex = exp ( ex1 ) * del**jj(j)

    att = aat(j-36)
    tx = atz(j-36)
    tau = ( t / tx ) - 1.0D+00

    ex2 = - att * tau * tau

    ex2 = max ( ex2, -225.0D+00 )
    ex2 = min ( ex2,  225.0D+00 )
    q10 = dex * exp ( ex2 )

    qm = dble ( jj(j) ) / del - dble ( k ) * aad(j-36) * del**(k-1)
    fct = qm * rho**2 * q10 / ddz

    q5t = fct * ( 2.0D+00 / rho + qm / ddz ) - ( rho / ddz )**2 * q10 * &
      ( dble ( jj(j) ) / del**2 + dble ( k * ( k - 1 ) ) * aad(j-36) * &
      del**(k-2) )

    dpdrr = dpdrr + q5t * g(j)
    qp = qp + g(j) * fct
    dadt = dadt - 2.0D+00 * g(j) * att * tau * q10 / tx
    dpdtr = dpdtr - 2.0D+00 * g(j) * att * tau * fct / tx

    q2a = q2a + t * g(j) * att * ( 4.0D+00 * ex2 + 2.0D+00 ) * q10 / tx**2

    ar = ar + q10 * g(j) / ( gascon() * t )

  end do

  cvr = cvr + q2a / gascon()
  pr = pr + qp
  sr = - dadt / gascon()
  ur = ar + sr
!
!  Assign dimensions.
!
  ar =  gascon() * t *  ar
  cvr = gascon() *     cvr
  sr =  gascon() *      sr
  ur =  gascon() * t *  ur
!
!  Adjust energies.
!
  ar = ar + gascon ( ) * t * s_ref - gascon ( ) * u_ref
  sr = sr - gascon ( ) * s_ref
  ur = ur - gascon ( ) * u_ref

  gr = ar + pr / rho
  hr = ur + pr / rho

  return
end
subroutine secvir ( t, vir )

!*****************************************************************************80
!
!! SECVIR calculates the second virial coefficient at a given temperature.
!
!  Discussion:
!
!    The second virial coefficient VIR appears in the first correction term
!    to the ideal gas equation of state:
!
!      P = R * T / volume + VIR / volume**2 + ...
!
!  Modified:
!
!    28 November 1998
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Output, real ( kind = 8 ) VIR, the second virial coefficient, in CM3/G.
!
  implicit none

  real ( kind = 8 ) b1
  real ( kind = 8 ) b1t
  real ( kind = 8 ) b1tt
  real ( kind = 8 ) b2
  real ( kind = 8 ) b2t
  real ( kind = 8 ) b2tt
  real ( kind = 8 ), parameter, dimension ( 5 ) :: g = (/ &
    -0.53062968529023D+03,  0.22744901424408D+04, -0.26627944829770D+04, &
     0.78779333020687D+03, -0.69830527374994D+02 /)
  real ( kind = 8 ) gascon
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: t_ref = 647.073D+00
  real ( kind = 8 ) v
  real ( kind = 8 ) vir
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SECVIR - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if

  call bb ( t, b1, b2, b1t, b2t, b1tt, b2tt )

  v = t_ref / t

  vir = b2 + ( &
      v    * ( g(1) &
    + v    * ( g(2) &
    + v    * ( g(3) &
    + v    * ( g(4) &
    + v**2 *   g(5) ))))) &
    / ( gascon ( ) * t )

  return
end
subroutine secvir_values ( n, tc, vir )

!*****************************************************************************80
!
!! SECVIR_VALUES returns some values of the second virial coefficient.
!
!  Modified:
!
!    03 February 2002
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
!    TJ270.H3, pages 24-25.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) VIR, the second virial coefficient, in
!    m^3/kg.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_data = 19

  integer ( kind = 4 ) n
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_data ) :: tcvec = (/ &
        0.0D+00,    5.0D+00,   10.0D+00,   20.0D+00,   30.0D+00, &
       40.0D+00,   60.0D+00,   90.0D+00,  120.0D+00,  150.0D+00, &
      180.0D+00,  210.0D+00,  240.0D+00,  300.0D+00,  400.0D+00, &
      500.0D+00,  700.0D+00, 1000.0D+00, 2000.0D+00 /)
  real ( kind = 8 ) vir
  real ( kind = 8 ), save, dimension ( n_data ) :: virvec = (/ &
    -98.96D+00, -90.08D+00, -82.29D+00, -69.36D+00, -59.19D+00, &
    -51.07D+00, -39.13D+00, -27.81D+00, -20.83D+00, -16.21D+00, &
    -12.98D+00, -10.63D+00,  -8.85D+00,  -6.39D+00,  -4.03D+00, &
     -2.71D+00,  -1.32D+00,  -0.39D+00,   0.53D+00 /)

  if ( n < 0 ) then
    n = 0
  end if

  if ( n_data <= n ) then
    n = 0
    tc = 0.0D+00
    vir = 0.0D+00
  else
    n = n + 1
    tc = tcvec(n)
    vir = virvec(n)
  end if

  return
end
subroutine sound ( t, p, c )

!*****************************************************************************80
!
!! SOUND computes the speed of sound given temperature and pressure.
!
!  Modified:
!
!    22 November 1998
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Input, real ( kind = 8 ) P, the pressure, in MegaPascals.
!
!    Output, real ( kind = 8 ) C, the speed of sound, in m/s.
!
  implicit none

  real ( kind = 8 ) ab
  real ( kind = 8 ) ai
  real ( kind = 8 ) ar
  real ( kind = 8 ) c
  real ( kind = 8 ) cp
  real ( kind = 8 ) cpi
  real ( kind = 8 ) cv
  real ( kind = 8 ) cvb
  real ( kind = 8 ) cvi
  real ( kind = 8 ) cvr
  real ( kind = 8 ) dpdr
  real ( kind = 8 ) dpdrb
  real ( kind = 8 ) dpdrr
  real ( kind = 8 ) dpdt
  real ( kind = 8 ) dpdtb
  real ( kind = 8 ) dpdtr
  real ( kind = 8 ) gb
  real ( kind = 8 ) gi
  real ( kind = 8 ) gr
  real ( kind = 8 ) hb
  real ( kind = 8 ) hi
  real ( kind = 8 ) hr
  real ( kind = 8 ) p
  real ( kind = 8 ) pb
  real ( kind = 8 ) pr
  real ( kind = 8 ) rho
  real ( kind = 8 ) rhol
  real ( kind = 8 ) rho_start
  real ( kind = 8 ) rhov
  real ( kind = 8 ) sb
  real ( kind = 8 ) si
  real ( kind = 8 ) sr
  real ( kind = 8 ) t
  real ( kind = 8 ) t2
  real ( kind = 8 ) ub
  real ( kind = 8 ) ui
  real ( kind = 8 ) ur
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SOUND - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative pressures.
!
  if ( p <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SOUND - Fatal error!'
    write ( *, '(a)' ) '  The input pressure P must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was P = ', p
    stop
  end if
!
!  For the given pressure, compute the saturation temperature.
!
  rhol = 0.0D+00
  rhov = 0.0D+00

  call tsat ( p, t2, rhol, rhov )
!
!  Depending on whether the temperature is above or below the
!  saturation temperature, we expect to compute the density of
!  a liquid or vapor.
!
  if ( t < t2 ) then
    rho_start = 1.9D+00
  else
    rho_start = 0.01D+00
  end if

  call dense ( p, t, rho_start, rho, dpdr )
!
!  From T and RHO, compute the thermodynamic properties.
!
  call ideal ( t, ai, cpi, cvi, gi, hi, si, ui )

  call resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )

  call base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )

  cv =     cvb +   cvr + cvi
  dpdr = dpdrb + dpdrr
  dpdt = dpdtb + dpdtr

  cp = cv + t * dpdt**2 / ( dpdr * rho**2 )

  c = sqrt ( 1000.0D+00 * cp * dpdr / cv )

  return
end
subroutine sound_values ( n, tc, p, c )

!*****************************************************************************80
!
!! SOUND_VALUES returns some values of the speed of sound.
!
!  Modified:
!
!    04 February 2002
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
!    TJ270.H3, page 238-246.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) C, the speed of sound, in m/s.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_data = 20

  real ( kind = 8 ) c
  real ( kind = 8 ), save, dimension ( n_data ) :: cvec = (/ &
    1401.0D+00,  472.8D+00, 533.7D+00, 585.7D+00, 609.5D+00, &
     632.2D+00,  674.6D+00, 713.9D+00, 802.0D+00, 880.1D+00, &
    1017.8D+00, 1115.9D+00, 1401.7D+00,1402.6D+00, 1409.6D+00, &
    1418.1D+00, 1443.1D+00, 1484.6D+00, 1577.1D+00, 1913.4D+00 /)
  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_data ) :: pvec = (/ &
      1.0D+00,   1.0D+00,   1.0D+00,    1.0D+00,    1.0D+00, &
      1.0D+00,   1.0D+00,   1.0D+00,    1.0D+00,    1.0D+00, &
      1.0D+00,   1.0D+00,   5.0D+00,  10.0D+00,    50.0D+00, &
    100.0D+00, 250.0D+00, 500.0D+00, 1000.0D+00, 2500.0D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_data ) :: tcvec = (/ &
       0.0D+00,  100.0D+00,  200.0D+00,  300.0D+00,  350.0D+00, &
     400.0D+00,  500.0D+00,  600.0D+00,  850.0D+00, 1100.0D+00, &
    1600.0D+00, 2000.0D+00,    0.0D+00,    0.0D+00,    0.0D+00, &
       0.0D+00,    0.0D+00,    0.0D+00,    0.0D+00,    0.0D+00 /)

  if ( n < 0 ) then
    n = 0
  end if

  if ( n_data <= n ) then
    n = 0
    tc = 0.0D+00
    p = 0.0D+00
    c = 0.0D+00
  else
    n = n + 1
    tc = tcvec(n)
    p = pvec(n)
    c = cvec(n)
  end if

  return
end
subroutine surten ( t, sigma )

!*****************************************************************************80
!
!! SURTEN returns the surface tension as a function of temperature.
!
!  Discussion:
!
!    SURTEN uses an equation that yields values of the surface tension to
!    within the accuracy of measurements from the triple point to the
!    critical point.
!
!      Sigma = B * ( (TSTAR-T)/TSTAR)**Mu * (1+b*(TSTAR-T)/TSTAR)
!
!    where:
!
!      TSTAR = 647.15 Degrees Kelvin,
!      B = 0.2358 Pascals * Meters
!      b = -0.625,
!      Mu = 1.256.
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Output, real ( kind = 8 ) SIGMA, the surface tension,
!    in Pascal * m = Newton / m.
!
  implicit none

  real ( kind = 8 ), parameter :: b_cap = 0.2358D+00
  real ( kind = 8 ), parameter :: b_small = -0.625D+00
  real ( kind = 8 ), parameter :: mu = 1.256D+00
  real ( kind = 8 ) sigma
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: t_star = 647.15D+00
  real ( kind = 8 ) term
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SURTEN - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if

  term = ( t_star - t ) / t_star
  sigma = b_cap * term**mu * ( 1.0D+00 + b_small * term )
!
!  Need this conversion to match the table, but justification is there none.
!
  sigma = 1000.0D+00 * sigma

  return
end
subroutine surten_values ( n, tc, sigma )

!*****************************************************************************80
!
!! SURTEN_VALUES returns some values of the surface tension.
!
!  Modified:
!
!    04 February 2002
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
!    TJ270.H3, pages 267.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) SIGMA, the surface tension,
!    in Pascal * m = Newton / m.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_data = 14

  integer ( kind = 4 ) n
  real ( kind = 8 ) sigma
  real ( kind = 8 ), save, dimension ( n_data ) :: sigmavec = (/ &
    74.22D+00, 72.74D+00, 71.20D+00, 69.60D+00, 67.95D+00, &
    58.92D+00, 48.75D+00, 37.68D+00, 26.05D+00, 14.37D+00, &
     8.78D+00,  3.67D+00,  0.40D+00,  0.0D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_data ) :: tcvec = (/ &
     10.0D+00,  20.0D+00,  30.0D+00,  40.0D+00,   50.0D+00, &
    100.0D+00, 150.0D+00, 200.0D+00, 250.0D+00,  300.0D+00, &
    325.0D+00, 350.0D+00, 370.0D+00, 373.976D+00 /)

  if ( n < 0 ) then
    n = 0
  end if

  if ( n_data <= n ) then
    n = 0
    tc = 0.0D+00
    sigma = 0.0D+00
  else
    n = n + 1
    tc = tcvec(n)
    sigma = sigmavec(n)
  end if

  return
end
subroutine tdpsdt ( t, dp )

!*****************************************************************************80
!
!! TDPSDT computes the quantity T * dP(Sat)/dT.
!
!  Discussion:
!
!    Here T is the temperature and P(Sat) is the vapor pressure.
!    It is used by TSAT_EST and TSAT.
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Output, real ( kind = 8 ) DP, the value T*(dP(Sat)/dT),
!    in MegaPascals.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
      -7.8889166D+00,   2.5514255D+00,   -6.716169D+00, 33.239495D+00, &
    -105.38479D+00,   174.35319D+00,   -148.39348D+00,  48.631602D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) dp
  integer ( kind = 4 ) i
  real ( kind = 8 ) q
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: t_ref = 647.25D+00
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TDPSDT - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if

  v = t / t_ref
  w = 1.0D+00 - v
  b = 0.0D+00
  c = 0.0D+00
  do i = 1, 8
    z = dble ( i + 1 ) / 2.0D+00
    y = a(i) * w**z
    c = c + ( y / w ) * ( 0.5D+00 - 0.5D+00 * dble ( i ) - 1.0D+00 / v )
    b = b + y
  end do

  q = b / v
  dp = 22.093D+00 * exp ( q ) * c

  return
end
subroutine thercon ( t, rho, lambda )

!*****************************************************************************80
!
!! THERCON calculates thermal conductivity for given temperature and density.
!
!  Modified:
!
!    20 November 1998
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Input, real ( kind = 8 ) RHO, the density, in G/CM3.
!
!    Output, real ( kind = 8 ) LAMBDA, the thermal conductivity,
!    in mW/(m degrees Kelvin).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), parameter, dimension ( 0:3 ) :: acof = (/ &
    2.02223D+00, 14.11166D+00, 5.25597D+00, -2.01870D+00 /)
  real ( kind = 8 ), parameter :: a_con = 18.66D+00
  real ( kind = 8 ) b(0:4,0:5)
  real ( kind = 8 ), parameter :: b_con = 1.00D+00
  real ( kind = 8 ), parameter :: c_con = 3.7711D-08
  real ( kind = 8 ) chi
  real ( kind = 8 ) cjth
  real ( kind = 8 ) cjtt
  real ( kind = 8 ) cp
  real ( kind = 8 ) cv
  real ( kind = 8 ) dpdr
  real ( kind = 8 ) dpdr2
  real ( kind = 8 ) dpdt
  real ( kind = 8 ) eta
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lambda
  real ( kind = 8 ) lambda0
  real ( kind = 8 ) lambda_del
  real ( kind = 8 ), parameter :: omega = 0.4678D+00
  real ( kind = 8 ) p
  real ( kind = 8 ), parameter :: p_ref = 22.115D+00
  real ( kind = 8 ) power
  real ( kind = 8 ) rho
  real ( kind = 8 ), parameter :: rho_ref = 317.763D+00
  real ( kind = 8 ) rho2
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: t_ref = 647.27D+00
  real ( kind = 8 ) total
  real ( kind = 8 ) u

  data b / &
    1.3293046D+00,    1.7018363D+00,   5.2246158D+00,  &
    8.7127675D+00, -1.8525999D+00, &
   -0.40452437D+00,  -2.2156845D+00, -10.124111D+00,   &
   -9.5000611D+00,  0.93404690D+00, &
    0.24409490D+00,   1.6511057D+00,   4.9874687D+00, &
    4.3786606D+00,  0.0D+00, &
    0.018660751D+00, -0.76736002D+00, -0.27297694D+00, &
    -0.91783782D+00, 0.0D+00, &
    -0.12961068D+00,  0.37283344D+00, -0.43083393D+00,  &
     0.0D+00,        0.0D+00, &
    0.044809953D+00, -0.11203160D+00,  0.13333849D+00,  &
    0.0D+00,        0.0D+00 /
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'THERCON - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'THERCON - Fatal error!'
    write ( *, '(a)' ) '  The input density RHO must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was RHO = ', rho
    stop
  end if
!
!  Compute DPDR, DPDT, ETA.
!
  call therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

  call viscosity ( t, rho, eta )
!
!  Convert RHO from G/CM3 to kg/M3,
!  Convert DPDR from ? to ?.
!
  rho2 = 1000.0D+00 * rho
  dpdr2 = dpdr / 1000.0D+00
!
!  Compute LAMBDA0.
!
  total = 0.0D+00
  do i = 0, 3
    total = total + acof(i) * ( t_ref / t )**i
  end do

  lambda0 = sqrt ( t / t_ref ) / total
!
!  Compute CHI.
!
  chi = rho2 * p_ref / ( rho_ref**2 * dpdr2 )
!
!  Compute delta_Lambda
!
  power = - a_con * ( ( t_ref - t ) / t )**2 - b_con * ( ( rho2 - rho_ref ) &
    / rho_ref )**4

  lambda_del = ( c_con / eta ) * ( ( t * rho_ref ) / ( t_ref * rho ) )**2 &
    * ( t_ref / p_ref )**2 * dpdt**2 * chi**omega * sqrt ( rho2 / rho_ref ) &
    * exp ( power )
!
!  Compute LAMBDA.
!
  total = 0.0D+00
  do i = 0, 4
    do j = 0, 5
      total = total + b(i,j) * ( ( t_ref - t ) / t )**i * &
        ( ( rho2 - rho_ref ) / rho_ref )**j
    end do
  end do

  lambda = lambda0 * exp ( ( rho2 / rho_ref ) * total ) + lambda_del
!
!  Temporary fix.
!
  lambda = 1000.0D+00 * lambda

  return
end
subroutine thercon_values ( n, tc, p, lambda )

!*****************************************************************************80
!
!! THERCON_VALUES returns some values of the thermal conductivity.
!
!  Modified:
!
!    04 February 2002
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
!    TJ270.H3, page 264.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) LAMBDA, the thermal conductivity, in
!    mW/(m degrees Kelvin).
!
  implicit none

  integer ( kind = 4 ), parameter :: n_data = 35

  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_data ) :: lambdavec = (/ &
    561.0D+00, 561.3D+00, 561.5D+00, 562.4D+00,  563.7D+00, &
    565.1D+00, 566.5D+00, 567.9D+00, 569.3D+00,  570.6D+00, &
    572.0D+00, 573.4D+00, 574.8D+00, 576.1D+00,  577.5D+00, &
    580.2D+00, 582.9D+00, 585.5D+00, 588.1D+00,  590.7D+00, &
    593.3D+00, 595.8D+00, 598.3D+00, 603.1D+00,  607.8D+00, &
    612.2D+00, 607.2D+00, 643.6D+00, 666.8D+00,   25.08D+00, &
    28.85D+00,  33.28D+00, 54.76D+00, 79.89D+00, 107.3D+00 /)
  real ( kind = 8 ), save, dimension ( n_data ) :: pvec = (/ &
       1.0D+00,    5.0D+00,   10.0D+00,   25.0D+00,   50.0D+00, &
      75.0D+00,  100.0D+00,  125.0D+00,  150.0D+00,  175.0D+00, &
     200.0D+00,  225.0D+00,  250.0D+00,  275.0D+00,  300.0D+00, &
     350.0D+00,  400.0D+00,  450.0D+00,  500.0D+00,  550.0D+00, &
     600.0D+00,  650.0D+00,  700.0D+00,  800.0D+00,  900.0D+00, &
    1000.0D+00,    1.0D+00,    1.0D+00,    1.0D+00,    1.0D+00, &
       1.0D+00,    1.0D+00,    1.0D+00,    1.0D+00,    1.0D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_data ) :: tcvec = (/ &
    0.0D+00,     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
    0.0D+00,     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
    0.0D+00,     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
    0.0D+00,     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
    0.0D+00,     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
    0.0D+00,    25.0D+00,  50.0D+00,  75.0D+00, 100.0D+00, &
    150.0D+00, 200.0D+00, 400.0D+00, 600.0D+00, 800.0D+00 /)

  if ( n < 0 ) then
    n = 0
  end if

  if ( n_data <= n ) then
    n = 0
    tc = 0.0D+00
    p = 0.0D+00
    lambda = 0.0D+00
  else
    n = n + 1
    tc = tcvec(n)
    p = pvec(n)
    lambda = lambdavec(n)
  end if

  return
end
subroutine therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )

!*****************************************************************************80
!
!! THERM calculates thermodynamic functions given temperature and density.
!
!  Discussion:
!
!    Thermodynamic values were calculated from an analytic equation
!    that approximates the Helmholtz function (specific Helmholtz
!    energy) for ordinary water and steam, of the form A=A(rho,T)
!    where A is the Helmholtz function, rho the density, and T
!    the absolute (thermodynamic) temperature.  Any thermodynamic
!    value for any state, liquid, vapor or metastable, may be
!    calculated by differentiation of this equation in accord with
!    the first and second laws of thermodynamics.
!
!    The International Association for the Properties of Steam
!    has provisionally accepted this formulation for the range
!    273.15 <= T <= 1273.15 degrees Kelvin, where, for 423.15 <= T,
!    the maximum pressure is Pmax = 1500 MPa = 15000 bar, and for
!    273.15 <= T < 423.15, the maximum pressure is
!    Pmax = 100 * (5 + (T-273.15)/15) MPa.
!
!    Close to the critical point, a small region is excluded:
!    Abs(T-Tk) < 1, abs((rho-rhok)/rhok) < 0.3.
!
!    The equation has a wider useful range, namely, fluid states
!    of pure, undissociated water and steam defined by
!    260 <= T <= 2500 K and 0 <= P <= 3000 MPa.
!
!    Thermodynamic property values for specific volume, density,
!    specific internal energy, specific enthalpy, and specific
!    entropy of water and steam were tabulated over the range
!    0 <= t <= 2000 C, 0 <= P <= 3000 MPa.  The reference
!    state is the liquid at the triple point, for which the
!    internal energy and entropy have been assigned the value zero.
!
!    Thermodynamic quantities are determined from the Helmholtz function
!    A(rho,T), which is computed as the sum of three terms:
!
!      A(rho,T) = A(base)(rho,T) + A(residual)(rho,T) + A(ideal)(T)
!                                                       (Equation 1)
!
!    Because A(rho,T) is everywhere single valued and analytic,
!    we can derive closed form relations for all other properties.
!    In the following, unless otherwise indicated, the independent
!    variables are temperature T and density RHO, and differentiation
!    with respect to one variable is to imply that the other is fixed.
!
!    Pressure:                  P       = RHO**2 * dA/dRHO
!    Density derivative:        dP/dRHO = 2*P/RHO + RHO**2 * d2A/dRHO2
!    Temperature derivative:    dP/dT   = RHO**2 * d2A/(dRHO dT)
!    Specific entropy:          S       = - dA/dT
!    Specific internal energy:  U       = A + T*S
!    Specific enthalpy:         H       = U + P/RHO
!    Specific heat capacity
!      at constant volume:      Cv      = - T * d2A/dT2
!    Specific Gibbs function:   G       = A + P/RHO
!    Specific heat capacity
!      at constant pressure:    Cp      = Cv + (T*(dP/dT)**2)/(RHO**2*dP/dRHO)
!    Speed of sound:            Omega   = Sqrt ((Cp/Cv) * dP/dRHO)
!    Second virial coefficient: B       = 1/(2*R*T) * (d2P/dRHO2) (at RHO=0)
!    Isothermal Joule-Thomson
!      coefficient:             DeltaT  = (dH/dP) (fixed T) =
!                                         (1/RHO)-(T*dP/dT)/(RHO**2*dP/dRHO)
!    Joule-Thomson coefficient: Mu      = (dT/dP) (fixed H) = DeltaT/Cp
!    Isentropic temperature-
!      pressure coefficient:    BetaS   = (dT/dP) (fixed S) =
!                                         (DeltaT - 1/RHO)/Cp
!
!  Modified:
!
!    19 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Input, real ( kind = 8 ) RHO, the fluid density, in G/CM3.
!
!    Output, real ( kind = 8 ) A, the Helmholtz function, in KJ/kg.
!
!    Output, real ( kind = 8 ) CJTH, the Joule-Thomson coefficient,
!    in K/MegaPascals.
!
!    Output, real ( kind = 8 ) CJTT, the isothermal Joule-Thomson coefficient,
!    in CM3/G.
!
!    Output, real ( kind = 8 ) CP, the isobaric (constant pressure) heat
!    capacity, in KJ/(kg degrees Kelvin).
!
!    Output, real ( kind = 8 ) CV, the isochoric (constant volume) heat capacity,
!    in KJ/(kg degrees Kelvin).
!
!    Output, real ( kind = 8 ) DPDR, the partial derivative
!    dP(T,RHO)/dRHO, with T held fixed, in MegaPascals*CM3/G.
!
!    Output, real ( kind = 8 ) DPDT, the partial derivative
!    dP(T,RHO)/dT, with RHO held fixed, in MegaPascals/degrees Kelvin.
!
!    Output, real ( kind = 8 ) G, the Gibbs free energy, in KJ/kg.
!
!    Output, real ( kind = 8 ) H, the enthalpy, in KJ/kg.
!
!    Output, real ( kind = 8 ) P, the pressure, in MegaPascals.
!
!    Output, real ( kind = 8 ) S, the entropy, in KJ/(kg degrees Kelvin).
!
!    Output, real ( kind = 8 ) U, the internal energy, in KJ/kg.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ab
  real ( kind = 8 ) ai
  real ( kind = 8 ) ar
  real ( kind = 8 ) cjth
  real ( kind = 8 ) cjtt
  real ( kind = 8 ) cp
  real ( kind = 8 ) cpi
  real ( kind = 8 ) cv
  real ( kind = 8 ) cvb
  real ( kind = 8 ) cvi
  real ( kind = 8 ) cvr
  logical, parameter :: debug = .false.
  real ( kind = 8 ) dpdr
  real ( kind = 8 ) dpdrb
  real ( kind = 8 ) dpdrr
  real ( kind = 8 ) dpdt
  real ( kind = 8 ) dpdtb
  real ( kind = 8 ) dpdtr
  real ( kind = 8 ) g
  real ( kind = 8 ) gb
  real ( kind = 8 ) gi
  real ( kind = 8 ) gr
  real ( kind = 8 ) h
  real ( kind = 8 ) hb
  real ( kind = 8 ) hi
  real ( kind = 8 ) hr
  real ( kind = 8 ) p
  real ( kind = 8 ) pb
  real ( kind = 8 ) pr
  real ( kind = 8 ) rho
  real ( kind = 8 ) s
  real ( kind = 8 ) sb
  real ( kind = 8 ) si
  real ( kind = 8 ) sr
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) ub
  real ( kind = 8 ) ui
  real ( kind = 8 ) ur
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'THERM - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'THERM - Fatal error!'
    write ( *, '(a)' ) '  The input density RHO must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was RHO = ', rho
    stop
  end if

  call ideal ( t, ai, cpi, cvi, gi, hi, si, ui )

  call resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )

  call base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )

  a =       ab +    ar +  ai
  cv =     cvb +   cvr + cvi

  if ( debug ) then
    write ( *, * ) ' '
    write ( *, * ) 'THERM:'
    write ( *, * ) '  CVB = ', cvb
    write ( *, * ) '  CVR = ', cvr
    write ( *, * ) '  CVI = ', cvi
    write ( *, * ) '  CV  = ', cv
  end if

  dpdr = dpdrb + dpdrr
  dpdt = dpdtb + dpdtr
  p =       pb +    pr
  s =       sb +    sr +  si
  u =       ub +    ur +  ui

  if ( debug ) then
    write ( *, * ) ' '
    write ( *, * ) 'THERM:'
    write ( *, * ) '  UB = ', ub
    write ( *, * ) '  UR = ', ur
    write ( *, * ) '  UI = ', ui
  end if

  g = a + p / rho
  h = u + p / rho
  cp = cv + t * dpdt**2 / ( dpdr * rho**2 )
  cjtt = 1.0D+00 / rho - t * dpdt / ( dpdr * rho**2 )
  cjth = - cjtt / cp

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
!
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
subroutine tsat ( p, t, rhol, rhov )

!*****************************************************************************80
!
!! TSAT calculates the saturation temperature for a given pressure.
!
!  Discussion:
!
!    The corresponding liquid and vapor densities are also computed.
!    The saturation temperature is also known as the "vapor temperature".
!
!  Modified:
!
!    04 February 2002
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the vapor pressure, in MegaPascals.
!
!    Output, real ( kind = 8 ) T, the vapor temperature, in degrees Kelvin.
!
!    Output, real ( kind = 8 ) RHOL, the liquid density, in G/CM3.
!
!    Output, real ( kind = 8 ) RHOV, the vapor density, in G/CM3.
!
  implicit none

  logical, parameter :: debug = .false.
  real ( kind = 8 ) delg
  real ( kind = 8 ) dp
  real ( kind = 8 ) dp2
  real ( kind = 8 ) errtol
  real ( kind = 8 ) gascon
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 50
  real ( kind = 8 ) p
  real ( kind = 8 ) p_consistent
  real ( kind = 8 ) rhol
  real ( kind = 8 ) rhov
  real ( kind = 8 ) t
!
!  Initialize output quantities, obliterating any input value.
!
  t = 0.0D+00
  rhol = 0.0D+00
  rhov = 0.0D+00
!
!  Set the error tolerance.
!
  errtol = sqrt ( epsilon ( errtol ) )
!
!  Refuse to handle zero or negative pressure.
!
  if ( p <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TSAT - Fatal error!'
    write ( *, '(a)' ) '  The input pressure must be positive!'
    write ( *, '(a,g14.6)' ) '  Your value was P = ', p
    stop
  end if
!
!  Estimate the saturation temperature.
!
  call tsat_est ( p, t )

  if ( debug ) then
    write ( *, * ) ' '
    write ( *, * ) 'TSAT:'
    write ( *, '(2g14.6)' ) p, t, rhol, rhov
  end if

  do it = 1, it_max

    call corr ( t, p, p_consistent, rhol, rhov, delg )

    dp = delg * gascon ( ) * t  * rhol * rhov / ( rhol - rhov )

    call tdpsdt ( t, dp2 )

    t = t * ( 1.0D+00 - dp / dp2 )

    if ( debug ) then
      write ( *, '(2g14.6)' ) p, t, rhol, rhov
    end if

    if ( abs ( delg ) < errtol ) then
      return
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TSAT - Warning!'
  write ( *, '(a)' ) '  The iteration did not converge.'
  write ( *, '(a,i6)' ) '  Number of iterations was ', it_max
  write ( *, '(a,g14.6)' ) '  Last iterate was ', t
  write ( *, '(a,g14.6)' ) '  Last DELG was ', delg

  return
end
subroutine tsat_est ( p, t )

!*****************************************************************************80
!
!! TSAT_EST makes a rough estimate of the saturation temperature.
!
!  Discussion:
!
!    The saturation temperature is also called the vapor temperature.
!
!  Modified:
!
!    02 February 2002
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the pressure, in MegaPascals.  The tabulated
!    range for P is
!      0.00061173 MegaPascals <= P <= P_CRIT = 22.055 MegaPascals.
!    The input value of P must be positive.
!
!    Output, real ( kind = 8 ) T, the saturation temperature,
!    in degrees Kelvin.  This value will always be in the range
!    [ 273.15, 647.126 ].
!
  implicit none

  integer ( kind = 4 ), parameter :: npol = 4

  real ( kind = 8 ), parameter, dimension ( 0:npol ) :: c = (/ &
    372.83D+00, 27.7589D+00, 2.3819D+00, 0.24834D+00, 0.0193855D+00 /)
  real ( kind = 8 ) dp
  real ( kind = 8 ) dt
  real ( kind = 8 ) errtol
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 8
  real ( kind = 8 ) p
  real ( kind = 8 ) pl
  real ( kind = 8 ), parameter :: p_crit = 22.055D+00
  real ( kind = 8 ) pp
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: t_crit = 647.126D+00
  real ( kind = 8 ), parameter :: t_min = 273.15D+00
  real ( kind = 8 ) t_old

  errtol = sqrt ( epsilon ( errtol ) )
!
!  Refuse to handle zero or negative pressure.
!
  if ( p <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TSAT_EST - Fatal error!'
    write ( *, '(a)' ) '  The input pressure must be positive!'
    write ( *, '(a,g14.6)' ) '  Your value was P = ', p
    stop
  end if

  if ( p_crit < p ) then
    t = t_crit
    return
  end if
!
!  The initial estimate for T uses a polyonmial in the logarithm of P.
!
  pl = 2.302585D+00 + log ( p )

  call dpoly_val_horner ( npol, c, pl, t )

  t = min ( t, t_crit )
  t = max ( t, t_min )

  dt = 0.0D+00

  do it = 1, it_max

    call psat_est ( t, pp )

    call tdpsdt ( t, dp )

    if ( abs ( p - pp ) < errtol * p ) then
      return
    end if

    dt = t * ( p - pp ) / dp

    t_old = t
    t = t * ( 1.0D+00 + ( p - pp ) / dp )
    t = min ( t, t_crit )
    t = max ( t, t_min )

    if ( abs ( dt ) < errtol * ( abs ( t ) + 1.0D+00 ) ) then
      return
    else if ( abs ( t - t_old ) < errtol ) then
      return
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TSAT_EST - Warning!'
  write ( *, '(a)' ) '  The iteration did not converge.'
  write ( *, '(a,i6)' ) '  Number of iterations was ', it_max
  write ( *, '(a,g14.6)' ) '  Convergence tolerance was ', errtol
  write ( *, '(a,g14.6)' ) '  Last iterate was ', t
  write ( *, '(a,g14.6)' ) '  Last correction was ', dt

  return
end
subroutine tsat_values ( n, p, tc )

!*****************************************************************************80
!
!! TSAT_VALUES returns some values of the saturation temperature.
!
!  Modified:
!
!    05 February 2002
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
!    TJ270.H3, pages 16-22.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) TC, the saturation temperature, in
!    degrees Celsius.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_data = 20

  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_data ) :: pvec = (/ &
    0.0061173D+00,  0.012D+00,   0.025D+00,  0.055D+00,  0.080D+00, &
    0.11D+00,       0.16D+00,    0.25D+00,   0.50D+00,   0.75D+00, &
    1.0D+00,        1.5D+00,   2.0D+00,      5.0D+00,   10.0D+00, &
    20.0D+00,      50.0D+00, 100.0D+00,    200.0D+00,  220.55D+00 /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_data ) :: tcvec = (/ &
      0.010D+00,  9.655D+00,   21.080D+00,  34.589D+00,  41.518D+00, &
     47.695D+00, 55.327D+00,   64.980D+00,  81.339D+00,  91.783D+00, &
     99.632D+00, 111.378D+00, 120.443D+00, 151.866D+00, 179.916D+00, &
    212.417D+00, 263.977D+00, 311.031D+00, 365.800D+00, 373.976D+00 /)

  if ( n < 0 ) then
    n = 0
  end if

  if ( n_data <= n ) then
    n = 0
    p = 0.0D+00
    tc = 0.0D+00

  else
    n = n + 1
    p = pvec(n)
    tc = tcvec(n)
  end if

  return
end
subroutine viscosity ( t, rho, eta )

!*****************************************************************************80
!
!! VISCOSITY calculates the viscosity for given temperature and density.
!
!  Discussion:
!
!    On 02 February 2002, I discovered that the Haar/Gallagher/Kell
!    reference apparently reversed the sign on the A3 coefficient.
!    That made the results better, but still off.
!
!    Apparently Haar/Gallagher/Kell had a transcription error in
!    the value of B(4,1), which they list as -0.273093, but which
!    should be -0.253093.
!
!    These two corrections courtesy of Meyer/McClintock/Silvestri/Spencer.
!
!    Now the results look proper!  And just 12 years late...
!
!  Modified:
!
!    02 February 2002
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
!    International Association for the Properties of Steam,
!    Release on Dynamic Viscosity of Water Substance,
!    National Bureau of Standards, Washington DC, 1975, revised 1983.
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers,
!    Fifth Edition, 1983,
!    TJ270.A75.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Input, real ( kind = 8 ) RHO, the density, in G/CM3.
!
!    Output, real ( kind = 8 ) ETA, the viscosity, in MegaPascal seconds.
!
  implicit none

  integer ( kind = 4 ), parameter :: npol_t = 3

  real ( kind = 8 ), parameter, dimension ( 0:npol_t ) :: a = (/ &
    0.0181583D+00, 0.0177624D+00, 0.0105287D+00, -0.0036744D+00 /)
  real ( kind = 8 ) arg
  real ( kind = 8 ), dimension(0:5,0:4) :: b = reshape ( (/ &
    0.501938D+00,   0.162888D+00,  -0.130356D+00,  &
    0.907919D+00,  -0.551119D+00,   0.146543D+00,  &
    0.235622D+00,   0.789393D+00,   0.673665D+00,  &
    1.207552D+00,   0.0670665D+00, -0.0843370D+00, &
   -0.274637D+00,  -0.743539D+00,  -0.959456D+00,  &
   -0.687343D+00,  -0.497089D+00,   0.195286D+00,  &
    0.145831D+00,   0.263129D+00,   0.347247D+00,  &
    0.213486D+00,   0.100754D+00,  -0.032932D+00,  &
   -0.0270448D+00, -0.0253093D+00, -0.0267758D+00, &
   -0.0822904D+00,  0.0602253D+00, -0.0202595D+00 /), &
    (/ 6,5 /) )
  logical, parameter :: debug = .false.
  real ( kind = 8 ) eta
  real ( kind = 8 ) eta0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) rho
  real ( kind = 8 ), parameter :: rho_max = 1.050D+00
  real ( kind = 8 ), parameter :: rho_ref = 0.317763D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) total
  real ( kind = 8 ), parameter :: t_max = 800.00D+00
  real ( kind = 8 ), parameter :: t_ref = 647.27D+00
  logical, save :: warning1 = .false.
  logical, save :: warning2 = .false.
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VISCOSITY - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  The input value was T = ', t
    stop
  end if

  if ( t_max < t .and. .not. warning1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VISCOSITY - Warning (once only)!'
    write ( *, '(a,g14.6)' ) &
      '  The input temperature T should be no more than ', t_max
    write ( *, '(a,g14.6)' ) '  The input value was T = ', t
    write ( *, '(a)' ) ' '
    warning1 = .true.
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VISCOSITY - Fatal error!'
    write ( *, '(a)' ) '  The input density RHO must be positive.'
    write ( *, '(a,g14.6)' ) '  The input value was RHO = ', rho
    stop
  end if

  if ( rho_max < rho .and. .not. warning2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VISCOSITY - Warning (once only)!'
    write ( *, '(a,g14.6)' ) &
      '  The input density RHO should be no more than ', rho_max
    write ( *, '(a,g14.6)' ) '  The input value was RHO = ', rho
    write ( *, '(a)' ) ' '
    warning2 = .true.
  end if
!
!  Compute ETA0.
!
  arg = t_ref / t

  call dpoly_val_horner ( npol_t, a, arg, total )

  eta0 = sqrt ( t / t_ref ) / total
!
!  Compute ETA.
!
  total = 0.0D+00
  do i = 0, 5
    do j = 0, 4
      total = total + b(i,j) * ( ( t_ref - t ) / t )**i * &
        ( ( rho - rho_ref ) / rho_ref )**j
    end do
  end do

  eta = eta0 * exp ( ( rho / rho_ref ) * total )

  return
end
subroutine viscosity_values ( n, tc, p, eta )

!*****************************************************************************80
!
!! VISCOSITY_VALUES returns some values of the viscosity function for testing.
!
!  Modified:
!
!    04 February 2002
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
!    TJ270.H3, page 263.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, real ( kind = 8 ) TC, the temperature, in degrees Celsius.
!
!    Output, real ( kind = 8 ) P, the pressure, in bar.
!
!    Output, real ( kind = 8 ) ETA, the viscosity, in MegaPascal seconds.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_data = 34

  real ( kind = 8 ) eta
  real ( kind = 8 ), save, dimension ( n_data ) :: etavec = (/ &
    1792.0D+00, 1791.0D+00, 1790.0D+00, 1786.0D+00, 1780.0D+00,  &
    1775.0D+00, 1769.0D+00, 1764.0D+00, 1759.0D+00, 1754.0D+00,  &
    1749.0D+00, 1744.0D+00, 1739.0D+00, 1735.0D+00, 1731.0D+00,  &
    1722.0D+00, 1714.0D+00, 1707.0D+00, 1700.0D+00, 1694.0D+00,  &
    1687.0D+00, 1682.0D+00, 1676.0D+00, 1667.0D+00, 1659.0D+00,  &
    1653.0D+00,  890.8D+00,  547.1D+00,  378.4D+00,   12.28D+00, &
      16.18D+00,  24.45D+00,  32.61D+00,  40.38D+00 /)
  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_data ) :: pvec = (/ &
       1.0D+00,    5.0D+00,   10.0D+00,   25.0D+00,   50.0D+00, &
      75.0D+00,  100.0D+00,  125.0D+00,  150.0D+00,  175.0D+00, &
     200.0D+00,  225.0D+00,  250.0D+00,  275.0D+00,  300.0D+00, &
     350.0D+00,  400.0D+00,  450.0D+00,  500.0D+00,  550.0D+00, &
     600.0D+00,  650.0D+00,  700.0D+00,  800.0D+00,  900.0D+00, &
    1000.0D+00,    1.0D+00,    1.0D+00,    1.0D+00,    1.0D+00, &
       1.0D+00,    1.0D+00,    1.0D+00,    1.0D+00  /)
  real ( kind = 8 ) tc
  real ( kind = 8 ), save, dimension ( n_data ) :: tcvec = (/ &
     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
     0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00,   0.0D+00, &
     0.0D+00,  25.0D+00,  50.0D+00,  75.0D+00, 100.0D+00, &
   200.0D+00, 400.0D+00, 600.0D+00, 800.0D+00 /)

  if ( n < 0 ) then
    n = 0
  end if

  if ( n_data <= n ) then
    n = 0
    tc = 0.0D+00
    p = 0.0D+00
    eta = 0.0D+00
  else
    n = n + 1
    tc = tcvec(n)
    p = pvec(n)
    eta = etavec(n)
  end if

  return
end
subroutine volume ( t, rho, v, dvdt, dvdr )

!*****************************************************************************80
!
!! VOLUME computes specific volume derivatives given temperature and density.
!
!  Discussion:
!
!    Because A(rho,T) is everywhere single valued and analytic,
!    we can derive closed form relations for all other properties.
!
!    The independent variables are temperature T and density RHO,
!    and differentiation with respect to one variable is to imply that
!    the other is held at a fixed value.
!
!  Modified:
!
!    28 November 1998
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
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the temperature, in degrees Kelvin.
!
!    Input, real ( kind = 8 ) RHO, the fluid density, in G/CM3.
!
!    Output, real ( kind = 8 ) V, the specific volume, in CM3/G.
!
!    Output, real ( kind = 8 ) DVDT, the partial derivative dV(T,RHO)/dT,
!    where V is the specific volume, in CM3 / (G * degrees Kelvin).
!
!    Output, real ( kind = 8 ) DVDR, the partial derivative dV(T,RHO)/dRHO,
!    where V is the specific volume, in CM3**2 / ( G**2 ).
!
  implicit none

  real ( kind = 8 ) ab
  real ( kind = 8 ) ar
  real ( kind = 8 ) cvb
  real ( kind = 8 ) cvr
  real ( kind = 8 ) dpdr
  real ( kind = 8 ) dpdrb
  real ( kind = 8 ) dpdrr
  real ( kind = 8 ) dpdt
  real ( kind = 8 ) dpdtb
  real ( kind = 8 ) dpdtr
  real ( kind = 8 ) dvdr
  real ( kind = 8 ) dvdt
  real ( kind = 8 ) gb
  real ( kind = 8 ) gr
  real ( kind = 8 ) hb
  real ( kind = 8 ) hr
  real ( kind = 8 ) pb
  real ( kind = 8 ) pr
  real ( kind = 8 ) rho
  real ( kind = 8 ) sb
  real ( kind = 8 ) sr
  real ( kind = 8 ) t
  real ( kind = 8 ) ub
  real ( kind = 8 ) ur
  real ( kind = 8 ) v
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VOLUME - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  The input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VOLUME - Fatal error!'
    write ( *, '(a)' ) '  The input density RHO must be positive.'
    write ( *, '(a,g14.6)' ) '  The input value was RHO = ', rho
    stop
  end if

  call resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )

  call base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )

  dpdr = dpdrb + dpdrr
  dpdt = dpdtb + dpdtr

  dvdt = dpdt / ( dpdr * rho**2 )
!
!  Because V = 1/Rho, dV/dRho = -1/Rho**2
!
  dvdr = - 1.0D+00 / rho**2

  v = 1.0D+00 / rho

  return
end
