subroutine rk1_ti_step ( x, t, h, q, fi, gi, seed, xstar )

!*****************************************************************************80
!
!! RK1_TI_STEP takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is first-order, and suitable for time-invariant
!    systems in which F and G do not depend explicitly on time.
!
!    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value at the current time.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) H, the time step.
!
!    Input, real ( kind = 8 ) Q, the spectral density of the input white noise.
!
!    Input, external real ( kind = 8 ) FI, the name of the deterministic
!    right hand side function.
!
!    Input, external real ( kind = 8 ) GI, the name of the stochastic
!    right hand side function.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) XSTAR, the value at time T+H.
!
  implicit none

  real ( kind = 8 ) a21
  real ( kind = 8 ), external :: fi
  real ( kind = 8 ), external :: gi
  real ( kind = 8 ) h
  real ( kind = 8 ) k1
  real ( kind = 8 ) q
  real ( kind = 8 ) q1
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) w1
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) xstar

  a21 =   1.0D+00

  q1 = 1.0D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
  k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1 

  xstar = x1 + a21 * k1

  return
end
subroutine rk2_ti_step ( x, t, h, q, fi, gi, seed, xstar )

!*****************************************************************************80
!
!! RK2_TI_STEP takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is second-order, and suitable for time-invariant
!    systems.
!
!    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value at the current time.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) H, the time step.
!
!    Input, real ( kind = 8 ) Q, the spectral density of the input white noise.
!
!    Input, external real ( kind = 8 ) FI, the name of the deterministic
!    right hand side function.
!
!    Input, external real ( kind = 8 ) GI, the name of the stochastic
!    right hand side function.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) XSTAR, the value at time T+H.
!
  implicit none

  real ( kind = 8 ) a21
  real ( kind = 8 ) a31
  real ( kind = 8 ) a32
  real ( kind = 8 ), external :: fi
  real ( kind = 8 ), external :: gi
  real ( kind = 8 ) h
  real ( kind = 8 ) k1
  real ( kind = 8 ) k2
  real ( kind = 8 ) q
  real ( kind = 8 ) q1
  real ( kind = 8 ) q2
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xstar

  a21 =   1.0D+00
  a31 =   0.5D+00
  a32 =   0.5D+00

  q1 = 2.0D+00
  q2 = 2.0D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
  k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1 

  t2 = t1 + a21 * h
  x2 = x1 + a21 * k1
  w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h )
  k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2

  xstar = x1 + a31 * k1 + a32 * k2

  return
end
subroutine rk3_ti_step ( x, t, h, q, fi, gi, seed, xstar )

!*****************************************************************************80
!
!! RK3_TI_STEP takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is third-order, and suitable for time-invariant
!    systems in which F and G do not depend explicitly on time.
!
!    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value at the current time.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) H, the time step.
!
!    Input, real ( kind = 8 ) Q, the spectral density of the input white noise.
!
!    Input, external real ( kind = 8 ) FI, the name of the deterministic
!    right hand side function.
!
!    Input, external real ( kind = 8 ) GI, the name of the stochastic
!    right hand side function.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) XSTAR, the value at time T+H.
!
  implicit none

  real ( kind = 8 ) a21
  real ( kind = 8 ) a31
  real ( kind = 8 ) a32
  real ( kind = 8 ) a41
  real ( kind = 8 ) a42
  real ( kind = 8 ) a43
  real ( kind = 8 ), external :: fi
  real ( kind = 8 ), external :: gi
  real ( kind = 8 ) h
  real ( kind = 8 ) k1
  real ( kind = 8 ) k2
  real ( kind = 8 ) k3
  real ( kind = 8 ) q
  real ( kind = 8 ) q1
  real ( kind = 8 ) q2
  real ( kind = 8 ) q3
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xstar

  a21 =   1.52880952525675D+00
  a31 =   0.0D+00
  a32 =   0.51578733443615D+00
  a41 =   0.53289582961739D+00
  a42 =   0.25574324768195D+00
  a43 =   0.21136092270067D+00

  q1 = 1.87653936176981D+00
  q2 = 3.91017166264989D+00
  q3 = 4.73124353935667D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
  k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1 

  t2 = t1 + a21 * h
  x2 = x1 + a21 * k1
  w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h )
  k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2

  t3 = t1 + a31 * h  + a32 * h
  x3 = x1 + a31 * k1 + a32 * k2
  w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h )
  k3 = h * fi ( x3 ) + h * gi ( x3 ) * w3

  xstar = x1 + a41 * k1 + a42 * k2 + a43 * k3

  return
end
subroutine rk4_ti_step ( x, t, h, q, fi, gi, seed, xstar )

!*****************************************************************************80
!
!! RK4_TI_STEP takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is fourth-order, and suitable for time-invariant
!    systems in which F and G do not depend explicitly on time.
!
!    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value at the current time.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) H, the time step.
!
!    Input, real ( kind = 8 ) Q, the spectral density of the input white noise.
!
!    Input, external real ( kind = 8 ) FI, the name of the deterministic
!    right hand side function.
!
!    Input, external real ( kind = 8 ) GI, the name of the stochastic
!    right hand side function.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) XSTAR, the value at time T+H.
!
  implicit none

  real ( kind = 8 ) a21
  real ( kind = 8 ) a31
  real ( kind = 8 ) a32
  real ( kind = 8 ) a41
  real ( kind = 8 ) a42
  real ( kind = 8 ) a43
  real ( kind = 8 ) a51
  real ( kind = 8 ) a52
  real ( kind = 8 ) a53
  real ( kind = 8 ) a54
  real ( kind = 8 ), external :: fi
  real ( kind = 8 ), external :: gi
  real ( kind = 8 ) h
  real ( kind = 8 ) k1
  real ( kind = 8 ) k2
  real ( kind = 8 ) k3
  real ( kind = 8 ) k4
  real ( kind = 8 ) q
  real ( kind = 8 ) q1
  real ( kind = 8 ) q2
  real ( kind = 8 ) q3
  real ( kind = 8 ) q4
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) w4
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) xstar

  a21 =   2.71644396264860D+00
  a31 = - 6.95653259006152D+00
  a32 =   0.78313689457981D+00
  a41 =   0.0D+00
  a42 =   0.48257353309214D+00
  a43 =   0.26171080165848D+00
  a51 =   0.47012396888046D+00
  a52 =   0.36597075368373D+00
  a53 =   0.08906615686702D+00
  a54 =   0.07483912056879D+00

  q1 =   2.12709852335625D+00
  q2 =   2.73245878238737D+00
  q3 =  11.22760917474960D+00
  q4 =  13.36199560336697D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
  k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1 

  t2 = t1 + a21 * h
  x2 = x1 + a21 * k1
  w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h )
  k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2

  t3 = t1 + a31 * h  + a32 * h
  x3 = x1 + a31 * k1 + a32 * k2
  w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h )
  k3 = h * fi ( x3 ) + h * gi ( x3 ) * w3

  t4 = t1 + a41 * h  + a42 * h + a43 * h
  x4 = x1 + a41 * k1 + a42 * k2
  w4 = r8_normal_01 ( seed ) * sqrt ( q4 * q / h )
  k4 = h * fi ( x4 ) + h * gi ( x4 ) * w4

  xstar = x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4

  return
end
subroutine rk1_tv_step ( x, t, h, q, fv, gv, seed, xstar )

!*****************************************************************************80
!
!! RK1_TV_STEP takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is first-order, and suitable for time-varying
!    systems.
!
!    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value at the current time.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) H, the time step.
!
!    Input, real ( kind = 8 ) Q, the spectral density of the input white noise.
!
!    Input, external real ( kind = 8 ) FV, the name of the deterministic
!    right hand side function.
!
!    Input, external real ( kind = 8 ) GV, the name of the stochastic
!    right hand side function.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) XSTAR, the value at time T+H.
!
  implicit none

  real ( kind = 8 ) a21
  real ( kind = 8 ), external :: fv
  real ( kind = 8 ), external :: gv
  real ( kind = 8 ) h
  real ( kind = 8 ) k1
  real ( kind = 8 ) q
  real ( kind = 8 ) q1
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) w1
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) xstar

  a21 =   1.0D+00

  q1 = 1.0D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
  k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1 

  xstar = x1 + a21 * k1

  return
end
subroutine rk2_tv_step ( x, t, h, q, fv, gv, seed, xstar )

!*****************************************************************************80
!
!! RK2_TV_STEP takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is second-order, and suitable for time-varying
!    systems.
!
!    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value at the current time.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) H, the time step.
!
!    Input, real ( kind = 8 ) Q, the spectral density of the input white noise.
!
!    Input, external real ( kind = 8 ) FV, the name of the deterministic
!    right hand side function.
!
!    Input, external real ( kind = 8 ) GV, the name of the stochastic
!    right hand side function.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) XSTAR, the value at time T+H.
!
  implicit none

  real ( kind = 8 ) a21
  real ( kind = 8 ) a31
  real ( kind = 8 ) a32
  real ( kind = 8 ), external :: fv
  real ( kind = 8 ), external :: gv
  real ( kind = 8 ) h
  real ( kind = 8 ) k1
  real ( kind = 8 ) k2
  real ( kind = 8 ) q
  real ( kind = 8 ) q1
  real ( kind = 8 ) q2
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xstar

  a21 =   1.0D+00
  a31 =   0.5D+00
  a32 =   0.5D+00

  q1 = 2.0D+00
  q2 = 2.0D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
  k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1 

  t2 = t1 + a21 * h
  x2 = x1 + a21 * k1
  w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h )
  k2 = h * fv ( t2, x2 ) + h * gv ( t2, x2 ) * w2

  xstar = x1 + a31 * k1 + a32 * k2

  return
end
subroutine rk4_tv_step ( x, t, h, q, fv, gv, seed, xstar )

!*****************************************************************************80
!
!! RK4_TV_STEP takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is fourth-order, and suitable for time-varying
!    systems.
!
!    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value at the current time.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) H, the time step.
!
!    Input, real ( kind = 8 ) Q, the spectral density of the input white noise.
!
!    Input, external real ( kind = 8 ) FV, the name of the deterministic
!    right hand side function.
!
!    Input, external real ( kind = 8 ) GV, the name of the stochastic
!    right hand side function.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) XSTAR, the value at time T+H.
!
  implicit none

  real ( kind = 8 ) a21
  real ( kind = 8 ) a31
  real ( kind = 8 ) a32
  real ( kind = 8 ) a41
  real ( kind = 8 ) a42
  real ( kind = 8 ) a43
  real ( kind = 8 ) a51
  real ( kind = 8 ) a52
  real ( kind = 8 ) a53
  real ( kind = 8 ) a54
  real ( kind = 8 ), external :: fv
  real ( kind = 8 ), external :: gv
  real ( kind = 8 ) h
  real ( kind = 8 ) k1
  real ( kind = 8 ) k2
  real ( kind = 8 ) k3
  real ( kind = 8 ) k4
  real ( kind = 8 ) q
  real ( kind = 8 ) q1
  real ( kind = 8 ) q2
  real ( kind = 8 ) q3
  real ( kind = 8 ) q4
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) w4
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) xstar

  a21 =   0.66667754298442D+00
  a31 =   0.63493935027993D+00
  a32 =   0.00342761715422D+00
  a41 = - 2.32428921184321D+00
  a42 =   2.69723745129487D+00
  a43 =   0.29093673271592D+00
  a51 =   0.25001351164789D+00
  a52 =   0.67428574806272D+00
  a53 = - 0.00831795169360D+00
  a54 =   0.08401868181222D+00

  q1 = 3.99956364361748D+00
  q2 = 1.64524970733585D+00
  q3 = 1.59330355118722D+00
  q4 = 0.26330006501868D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
  k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1 

  t2 = t1 + a21 * h
  x2 = x1 + a21 * k1
  w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h )
  k2 = h * fv ( t2, x2 ) + h * gv ( t2, x2 ) * w2

  t3 = t1 + a31 * h  + a32 * h
  x3 = x1 + a31 * k1 + a32 * k2
  w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h )
  k3 = h * fv ( t3, x3 ) + h * gv ( t3, x3 ) * w3

  t4 = t1 + a41 * h  + a42 * h  + a43 * h
  x4 = x1 + a41 * k1 + a42 * k2 + a43 * k3
  w4 = r8_normal_01 ( seed ) * sqrt ( q4 * q / h )
  k4 = h * fv ( t4, x4 ) + h * gv ( t4, x4 ) * w4

  xstar = x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4

  return
end
function r8_normal_01 ( seed )

!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    Because this routine uses the Box Muller method, it requires pairs
!    of uniform random values to generate a pair of normal random values.
!    This means that on every other call, the code can use the second
!    value that it calculated.
!
!    However, if the user has changed the SEED value between calls,
!    the routine automatically resets itself and discards the saved data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a normally distributed
!    random value.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_01
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: seed1 = 0
  integer ( kind = 4 ), save :: seed2 = 0
  integer ( kind = 4 ), save :: seed3 = 0
  integer ( kind = 4 ), parameter :: two = 2
  integer ( kind = 4 ), save :: used = 0
  real ( kind = 8 ) v1
  real ( kind = 8 ), save :: v2 = 0.0D+00
!
!  If USED is odd, but the input SEED does not match
!  the output SEED on the previous call, then the user has changed
!  the seed.  Wipe out internal memory.
!
  if ( mod ( used, two ) == 1 ) then

    if ( seed /= seed2 ) then
      used = 0
      seed1 = 0
      seed2 = 0
      seed3 = 0
      v2 = 0.0D+00
    end if

  end if
!
!  If USED is even, generate two uniforms, create two normals,
!  return the first normal and its corresponding seed.
!
  if ( mod ( used, two ) == 0 ) then

    seed1 = seed

    r1 = r8_uniform_01 ( seed )

    if ( r1 == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    seed2 = seed

    r2 = r8_uniform_01 ( seed )

    seed3 = seed

    v1 = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    v2 = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )

    r8_normal_01 = v1
    seed = seed2
!
!  If USED is odd (and the input SEED matched the output value from
!  the previous call), return the second normal and its corresponding seed.
!
  else

    r8_normal_01 = v2
    seed = seed3

  end if

  used = used + 1

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

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
