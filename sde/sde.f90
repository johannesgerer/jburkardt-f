subroutine bpath ( seed, n, w )

!*****************************************************************************80
!
!! BPATH performs a Brownian path simulation.
!
!  Discussion:
!
!    This routine computes one simulation of discretized Brownian 
!    motion over the time interval [0,1] using N time steps.
!    The user specifies a random number seed.  Different values of
!    the seed will result in different realizations of the path.
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
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Input, integer ( kind = 4 ) N, the number of steps.
!
!    Output, real ( kind = 8 ) W(0:N), the Brownian path.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dt
  real ( kind = 8 ) dw(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tmax
  real ( kind = 8 ) w(0:n)

  tmax = 1.0D+00
  dt = tmax / real ( n, kind = 8 )
!
!  Define the increments dW.
!
  call r8vec_normal_01 ( n, seed, dw )

  dw(1:n) = sqrt ( dt ) * dw(1:n)
!
!  W is the sum of the previous increments.
!
  w(0) = 0.0D+00
  do j = 1, n
    w(j) = w(j-1) + dw(j)
  end do

  return
end
subroutine bpath_gnuplot ( n, w )

!*****************************************************************************80
!
!! BPATH_GNUPLOT writes a GNUPLOT input file to plot BPATH data.
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
!    John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of steps.
!
!    Input, real ( kind = 8 ) W(0:N), the Brownian path.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  integer ( kind = 4 ) i
  real ( kind = 8 ) t
  real ( kind = 8 ) w(0:n)
!
!  Create the data file.
!
  call get_unit ( data_unit )

  data_filename = 'bpath_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 0, n
    t = real ( i, kind = 8 ) / real ( n, kind = 8 )
    write ( data_unit, '(2x,g14.6,2x,g14.6)' ) t, w(i)
  end do
  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BPATH data stored in "' &
    // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = 'bpath_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# bpath_commands.txt'
  write ( command_unit, '(a)' ) '# created by sde::bpath_gnuplot.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < bpath_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "bpath.png"'
  write ( command_unit, '(a)' ) 'set xlabel "t"'
  write ( command_unit, '(a)' ) 'set ylabel "W(t)"'
  write ( command_unit, '(a)' ) 'set title "Brownian motion by BPATH"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "bpath_data.txt" using 1:2'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  BPATH plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine bpath_average ( seed, m, n, u, umean, error )

!*****************************************************************************80
!
!! BPATH_AVERAGE: displays the average of 1000 Brownian paths.
!
!  Discussion:
!
!    This routine computes M simulations of discretized Brownian 
!    motion W(t) over the time interval [0,1] using N time steps.
!    The user specifies a random number seed.  Different values of
!    the seed will result in a different set of realizations of the path.
!
!    Actually, we are interested in a function u(W(t)):
!
!      u(W(t)) = exp ( t + W(t)/2 )
!
!    The routine plots 5 of the simulations, as well as the average
!    of all the simulations.  
!
!    The plot of the average should be quite smooth.  Its expected
!    value is exp ( 9 * t / 8 ), and we compute the 'error', that is,
!    the difference between the averaged value and this expected
!    value.  This 'error' should decrease as the number of simulation
!    is increased.
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
!    Original Matlab version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Input, integer ( kind = 4 ) M, the number of simulations to compute 
!    and average.  A typical value is 1000.
!
!    Input, integer ( kind = 4 ) N, the number of steps.  A typical value
!    is 500.
!
!    Output, real ( kind = 8 ) U(M,0:N), the M paths.
!
!    Output, real ( kind = 8 ) UMEAN(0:N), the averaged path.
!
!    Output, real ( kind = 8 ) ERROR, the maximum difference between the
!    averaged path and the exact expected value.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) dt
  real ( kind = 8 ) dw(n)
  real ( kind = 8 ) error
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(0:n)
  real ( kind = 8 ) tmax
  real ( kind = 8 ) u(1:m,0:n)
  real ( kind = 8 ) umean(0:n)
  real ( kind = 8 ) w(0:n)

  tmax = 1.0D+00
  dt = tmax / real ( n, kind = 8 )
  do j = 0, n
    t(j) = real ( j, kind = 8 ) * tmax / real ( n, kind = 8 )
  end do

  do i = 1, m
!
!  Define the increments dW.
!
    call r8vec_normal_01 ( n, seed, dw )

    dw(1:n) = sqrt ( dt ) * dw(1:n)
!
!  W is the sum of the previous increments.
!
    w(0) = 0.0D+00
    do j = 1, n
      w(j) = w(j-1) + dw(j)
    end do

    u(i,0:n) = exp ( t(0:n) + 0.5D+00 * w(0:n) )

  end do
!
!  Average the M estimates of the path.
!
  umean(0:n) = sum ( u(1:m,0:n), 1 ) / real ( m, kind = 8 )

  error = maxval ( abs ( umean(0:n) - exp ( 9.0D+00 * t(0:n) / 8.0D+00 ) ) )

  return
end
subroutine bpath_average_gnuplot ( m, n, u, umean )

!*****************************************************************************80
!
!! BPATH_AVERAGE_GNUPLOT writes a GNUPLOT input file to plot BPATH_AVERAGE data.
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
!    John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of simulations.
!
!    Input, integer ( kind = 4 ) N, the number of steps. 
!
!    Input, real ( kind = 8 ) U(M,0:N), the M paths.
!
!    Input, real ( kind = 8 ) UMEAN(0:N), the averaged path.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  integer ( kind = 4 ) i
  real ( kind = 8 ) t
  real ( kind = 8 ) u(m,0:n)
  real ( kind = 8 ) umean(0:n)
!
!  Create the data file.
!
  call get_unit ( data_unit )

  data_filename = 'bpath_average_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 0, n
    t = real ( i, kind = 8 ) / real ( n, kind = 8 )
    write ( data_unit, '(7(2x,g14.6))' ) t, u(1:5,i), umean(i)
  end do
  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BPATH_AVERAGE data stored in "' &
    // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = 'bpath_average_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# bpath_average_commands.txt'
  write ( command_unit, '(a)' ) '# created by sde::bpath_average_gnuplot.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < bpath_average_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "bpath_average.png"'
  write ( command_unit, '(a)' ) 'set xlabel "t"'
  write ( command_unit, '(a)' ) 'set ylabel "W(t)"'
  write ( command_unit, '(a)' ) 'set title "Averaged Brownian paths"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) &
    'plot "bpath_average_data.txt" using 1:2 title "sample 1", \'
  write ( command_unit, '(a)' ) &
    '     "bpath_average_data.txt" using 1:3 title "sample 2", \'
  write ( command_unit, '(a)' ) &
    '     "bpath_average_data.txt" using 1:4 title "sample 3", \'
  write ( command_unit, '(a)' ) &
    '     "bpath_average_data.txt" using 1:5 title "sample 4", \'
  write ( command_unit, '(a)' ) &
    '     "bpath_average_data.txt" using 1:6 title "sample 5", \'
  write ( command_unit, '(a)' ) &
    '     "bpath_average_data.txt" using 1:7 title "average" lw 3'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  BPATH_AVERAGE plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine chain ( seed, n, xem, vem, diff )

!*****************************************************************************80
!
!! CHAIN tests the stochastic Chain Rule.
!
!  Discussion:
!
!    This function solves a stochastic differential equation for 
!
!      V = sqrt(X) 
!
!    where X satisfies the stochastic differential equation:
! 
!      dX = ( alpha - X ) * dt + beta * sqrt(X) dW,
!      X(0) = Xzero,
!
!    with 
!
!      alpha = 2,
!      beta = 1,
!      Xzero = 1.
!
!    From the stochastic Chain Rule, the SDE for V is therefore:
!
!      dV = ( ( 4 * alpha - beta^2 ) / ( 8 * V ) - 1/2 V ) dt + 1/2 beta dW
!      V(0) = sqrt ( Xzero ).
!
!    Xem is the Euler-Maruyama solution for X. 
!
!    Vem is the Euler-Maruyama solution of the SDE for V from
!    the stochastic Chain Rule.
!
!    Hence, we compare sqrt(Xem) and Vem.
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
!    Original Matlab version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Input, integer ( kind = 4 ) N, the number of time steps.
!
!    Output, real ( kind = 8 ) XEM(0:N), the computed value of X.
!
!    Output, real ( kind = 8 ) VEM(0:N), the computed value of V.
!
!    Output, real ( kind = 8 ) DIFF, the maximum value of |sqrt(XEM)-V|.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) diff
  real ( kind = 8 ) dt
  real ( kind = 8 ) dt2
  real ( kind = 8 ) dw(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tmax
  real ( kind = 8 ) vem(0:n)
  real ( kind = 8 ) xem(0:n)
!
!  Set problem parameters.
!
  alpha = 2.0D+00
  beta = 1.0D+00
!
!  Stepping parameters.
!  dt2 is the size of the Euler-Maruyama steps.
!
  tmax = 1.0D+00
  dt = tmax / real ( n, kind = 8 )
  dt2 = dt
!
!  Define the increments dW.
!
  call r8vec_normal_01 ( n, seed, dw )

  dw(1:n) = sqrt ( dt ) * dw(1:n)
!
!  Solve for X(t).
!
  xem(0) = 1.0D+00
  do j = 1, n
    xem(j) = xem(j-1) + ( alpha - xem(j-1) ) * dt2 &
                      + beta * sqrt ( xem(j-1) ) * dw(j)
  end do
!
!  Solve for V(t).
!
  vem(0) = sqrt ( xem(0) )
  do j = 1, n
    vem(j) = vem(j-1) &
      + ( ( 4.0D+00 * alpha - beta ** 2 ) / ( 8.0D+00 * vem(j-1) ) &
      - 0.5D+00 * vem(j-1) ) * dt2 &
      + 0.5D+00 * beta * dw(j)
  end do
!
!  Compare sqrt(X) and V.
!
  diff = 0.0D+00
  do i = 0, n
    diff = max ( diff, abs ( sqrt ( xem(i) ) - vem(i) ) )
  end do

  return
end
subroutine chain_gnuplot ( n, x, v )

!*****************************************************************************80
!
!! CHAIN_GNUPLOT writes a GNUPLOT input file to plot CHAIN data.
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
!    John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of steps.
!
!    Input, real ( kind = 8 ) X(0:N), the value of X.
!
!    Input, real ( kind = 8 ) V(0:N), the value of V.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  integer ( kind = 4 ) i
  real ( kind = 8 ) t
  real ( kind = 8 ) v(0:n)
  real ( kind = 8 ) x(0:n)
!
!  Create the data file.
!
  call get_unit ( data_unit )

  data_filename = 'chain_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 0, n
    t = real ( i, kind = 8 ) / real ( n, kind = 8 )
    write ( data_unit, '(3(2x,g14.6))' ) t, sqrt ( x(i) ), v(i)
  end do
  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CHAIN data stored in "' &
    // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = 'chain_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# chain_commands.txt'
  write ( command_unit, '(a)' ) '# created by sde::chain_gnuplot.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < chain_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "chain.png"'
  write ( command_unit, '(a)' ) 'set xlabel "t"'
  write ( command_unit, '(a)' ) 'set ylabel "Sqrt(X(t)) vs V(X(t))"'
  write ( command_unit, '(a)' ) &
    'set title "V(X(t)) from X(t) and from Chain Rule"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) &
    'plot "chain_data.txt" using 1:2 title "Sqrt(X(t))", \'
  write ( command_unit, '(a)' ) &
    '     "chain_data.txt" using 1:3 title "V(X(t))"'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  CHAIN plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine em ( seed, n, t, xtrue, t2, xem, error )

!*****************************************************************************80
!
!! EM applies the Euler-Maruyama method to a linear SDE.
!
!  Discussion:
!
!    The SDE is 
!
!      dX = lambda * X dt + mu * X dW,   
!      X(0) = Xzero,
!
!    where 
!
!      lambda = 2,
!      mu = 1,
!      Xzero = 1.
!
!    The discretized Brownian path over [0,1] uses
!    a stepsize dt = 2^(-8).
!
!    The Euler-Maruyama method uses a larger timestep Dt = R*dt,
!    where R is an integer.  For an SDE of the form
!
!      dX = f(X(t)) dt + g(X(t)) dW(t)
!
!    it has the form
!
!      X(j) = X(j-1) + f(X(j-1)) * Dt + g(X(j-1)) * ( W(j*Dt) - W((j-1)*Dt) )
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
!    Original Matlab version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Input, integer ( kind = 4 ) N, the number of time steps.  A typical
!    value is 2^8.  N should be a multiple of 4.
!
!    Output, real ( kind = 8 ) T(0:N), the time values for the exact solution.
!
!    Output, real ( kind = 8 ) XTRUE(0:N), the exact solution.
!
!    Output, real ( kind = 8 ) T2(0:N/4), the time values for the 
!    Euler-Maruyama solution.
!
!    Output, real ( kind = 8 ) XEM(0:N/4), the Euler-Maruyama solution.
!
!    Output, real ( kind = 8 ) ERROR, the value of | XEM(T) - XTRUE(T) |.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dt
  real ( kind = 8 ) dt2
  real ( kind = 8 ) dw(n)
  real ( kind = 8 ) dw2
  real ( kind = 8 ) error
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  real ( kind = 8 ) lambda
  real ( kind = 8 ) mu
  integer ( kind = 4 ) r
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(0:n)
  real ( kind = 8 ) t2(0:n/4)
  real ( kind = 8 ) tmax
  real ( kind = 8 ) w(0:n)
  real ( kind = 8 ) xem(0:n/4)
  real ( kind = 8 ) xtrue(0:n)
  real ( kind = 8 ) xzero
!
!  Set problem parameters.
!
  lambda = 2.0D+00
  mu = 1.0D+00
  xzero = 1.0D+00
!
!  Set stepping parameters.
!
  tmax = 1.0D+00
  dt = tmax / real ( n, kind = 8 )
!
!  Define the increments dW.
!
  call r8vec_normal_01 ( n, seed, dw )

  dw(1:n) = sqrt ( dt ) * dw(1:n)
!
!  Sum the Brownian increments.
!
  w(0) = 0.0D+00
  do j = 1, n
    w(j) = w(j-1) + dw(j)
  end do

  do j = 0, n
    t(j) = real ( j, kind = 8 ) * tmax / real ( n, kind = 8 )
  end do
!
!  Compute the discretized Brownian path.
!
  xtrue(0:n) = xzero * exp ( ( lambda - 0.5D+00 * mu**2 ) &
    * ( t(0:n) + mu * w(0:n) ) )
!
!  Set:
!  R, the multiplier for the EM step, 
!  Dt, the EM stepsize,
!  L, the number of EM steps (we need N to be a multiple of R!)
!
  r = 4
  dt2 = real ( r, kind = 8 ) * dt
  l = n / r

  do j = 0, l
    t2(j) = real ( j, kind = 8 ) * tmax / real ( l, kind = 8 )
  end do
!
!  Preallocate Xem for efficiency.
!
  xem(0) = xzero
  do j = 1, l
    dw2 = sum ( dw ( r*(j-1)+1 : r*j ) )
    xem(j) = xem(j-1) + dt2 * lambda * xem(j-1) + mu * xem(j-1) * dw2
  end do

  error = abs ( xem(l) - xtrue(n) )

  return
end
subroutine em_gnuplot ( n, t, xtrue, t2, xem )

!*****************************************************************************80
!
!! EM_GNUPLOT writes a GNUPLOT input file to plot EM data.
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
!    John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of steps.
!
!    Input, real ( kind = 8 ) T(0:N), the time values for the exact solution.
!
!    Input, real ( kind = 8 ) XTRUE(0:N), the exact solution.
!
!    Input, real ( kind = 8 ) T2(0:N/4), the time values for the 
!    Euler-Maruyama solution.
!
!    Input, real ( kind = 8 ) XEM(0:N/4), the Euler-Maruyama solution.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  integer ( kind = 4 ) i
  real ( kind = 8 ) t(0:n)
  real ( kind = 8 ) t2(0:n/4)
  real ( kind = 8 ) xem(0:n/4)
  real ( kind = 8 ) xtrue(0:n)
!
!  Create data file #1.
!
  call get_unit ( data_unit )

  data_filename = 'em1_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 0, n
    write ( data_unit, '(3(2x,g14.6))' ) t(i), xtrue(i)
  end do
  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EM data #1 stored in "' &
    // trim ( data_filename ) // '".'
!
!  Create data file #2.
!
  call get_unit ( data_unit )

  data_filename = 'em2_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 0, n / 4
    write ( data_unit, '(3(2x,g14.6))' ) t2(i), xem(i)
  end do
  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EM data #2 stored in "' &
    // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = 'em_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# em_commands.txt'
  write ( command_unit, '(a)' ) '# created by sde::em_gnuplot.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < em_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "em.png"'
  write ( command_unit, '(a)' ) 'set xlabel "t"'
  write ( command_unit, '(a)' ) 'set ylabel "X(t)"'
  write ( command_unit, '(a)' ) &
    'set title "Exact X(t) and Euler-Maruyama Estimate"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) &
    'plot "em1_data.txt" using 1:2 title "Exact X(t))", \'
  write ( command_unit, '(a)' ) &
    '     "em2_data.txt" using 1:2 title "EM X(t)"'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  EM plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine emstrong ( seed, m, n, p_max, dtvals, xerr )

!*****************************************************************************80
!
!! EMSTRONG tests the strong convergence of the EM method.
!
!  Discussion:
!
!    The SDE is 
!
!      dX = lambda * X dt + mu * X dW,   
!      X(0) = Xzero,
!
!    where 
!
!      lambda = 2,
!      mu = 1,
!      Xzero = 1.
!
!    The discretized Brownian path over [0,1] has dt = 2^(-9).
!
!    The Euler-Maruyama method uses 5 different timesteps: 
!      16*dt, 8*dt, 4*dt, 2*dt, dt.
!
!    We are interested in examining strong convergence at T=1,
!    that is
!
!      E | X_L - X(T) |.
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
!    Original Matlab version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Input, integer ( kind = 4 ) M, the number of simulations to perform.
!    A typical value is M = 1000.
!
!    Input, integer ( kind = 4 ) N, the number of time steps to take.
!    A typical value is N = 512.
!
!    Input, integer ( kind = 4 ) P_MAX, the number of time step sizes to use.
!    A typical value is 5.
!
!    Output, real ( kind = 8 ) DTVALS(P_MAX), the time steps used.
!
!    Output, real ( kind = 8 ) XERR(P_MAX), the averaged absolute error in the
!    solution estimate at the final time.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p_max

  real ( kind = 8 ) a(p_max,2)
  real ( kind = 8 ) dt
  real ( kind = 8 ) dt2
  real ( kind = 8 ) dtvals(p_max)
  real ( kind = 8 ) dw(n)
  real ( kind = 8 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  real ( kind = 8 ) lambda
  real ( kind = 8 ) mu
  integer ( kind = 4 ) p
  real ( kind = 8 ) q
  integer ( kind = 4 ) r
  real ( kind = 8 ) resid
  real ( kind = 8 ) rhs(p_max)
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sol(2)
  real ( kind = 8 ) tmax
  real ( kind = 8 ) w(0:n)
  real ( kind = 8 ) winc
  real ( kind = 8 ) xerr(p_max)
  real ( kind = 8 ) xtemp
  real ( kind = 8 ) xtrue
  real ( kind = 8 ) xzero
!
!  Set problem parameters.
!
  lambda = 2.0D+00
  mu = 1.0D+00
  xzero = 1.0D+00
!
!  Set stepping parameters.
!
  tmax = 1.0D+00
  dt = tmax / real ( n, kind = 8 )

  do p = 1, p_max
    dtvals(p) = dt * 2.0D+00 ** ( p - 1 )
  end do
!
!  Sample over discrete Brownian paths.
!
  xerr(1:p_max) = 0.0D+00

  do s = 1, m
!
!  Define the increments dW.
!
    call r8vec_normal_01 ( n, seed, dw )

    dw(1:n) = sqrt ( dt ) * dw(1:n)
!
!  Sum the increments to get the Brownian path.
!
    w(0) = 0.0D+00
    do j = 1, n
      w(j) = w(j-1) + dw(j)
    end do
!
!  Determine the true solution.
!
    xtrue = xzero * exp ( ( lambda - 0.5 * mu ** 2 ) + mu * w(n) )
!
!  Use the Euler-Maruyama method with 5 different time steps dt2 = r * dt
!  to estimate the solution value at time TMAX.
!
    do p = 1, p_max                       
      dt2 = dtvals(p)
      r = 2 ** ( p - 1 )
      l = n / r
      xtemp = xzero
      do j = 1, l
        winc = sum ( dw(r*(j-1)+1:r*j) )
        xtemp = xtemp + dt2 * lambda * xtemp + mu * xtemp * winc
      end do
      xerr(p) = xerr(p) + abs ( xtemp - xtrue )
    end do

  end do

  xerr(1:p_max) = xerr(1:p_max) / real ( m, kind = 8 )
!
!  Least squares fit of error = c * dt^q.
!
  a(1:p_max,1) = 1.0D+00
  a(1:p_max,2) = log ( dtvals(1:p_max) )
  rhs(1:p_max) = log ( xerr(1:p_max) )

  call qr_solve ( p_max, 2, a, rhs, sol )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EMSTRONG:'
  write ( *, '(a)' ) '  Least squares solution to Error = c * dt ^ q'
  write ( *, '(a)' ) '  (Expecting Q to be about 1/2.)'
  write ( *, '(a,g14.6)' ) '  Computed Q = ', sol(2)

  resid = 0.0D+00
  do i = 1, p_max
    e = a(i,1) * sol(1) + a(i,2) * sol(2) - rhs(i)
    resid = resid + e * e
  end do
  resid = sqrt ( resid )
  write ( *, '(a,g14.6)' ) '  Residual is ', resid

  return
end
subroutine emstrong_gnuplot ( p_max, dtvals, xerr )

!*****************************************************************************80
!
!! EMSTRONG_GNUPLOT writes a GNUPLOT input file to plot EMSTRONG data.
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
!    John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_MAX, the number of time step sizes to use.
!
!    Input, real ( kind = 8 ) DTVALS(P_MAX), the time steps used.
!
!    Input, real ( kind = 8 ) XERR(P_MAX), the averaged absolute error in the
!    solution estimate at the final time.
!
  implicit none

  integer ( kind = 4 ) p_max

  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  integer ( kind = 4 ) i
  real ( kind = 8 ) dtvals(p_max)
  real ( kind = 8 ) xerr(p_max)
!
!  Create data file.
!
  call get_unit ( data_unit )

  data_filename = 'emstrong_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 1, p_max
    write ( data_unit, '(3(2x,g14.6))' ) dtvals(i), xerr(i), sqrt ( dtvals(i) )
  end do
  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EMSTRONG data stored in "' &
    // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = 'emstrong_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# emstrong_commands.txt'
  write ( command_unit, '(a)' ) '# created by sde::emstrong_gnuplot.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < emstrong_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "emstrong.png"'
  write ( command_unit, '(a)' ) 'set xlabel "Log(dt)"'
  write ( command_unit, '(a)' ) 'set ylabel "Log(Averaged Error at final T)"'
  write ( command_unit, '(a)' ) 'set logscale xy 10'
  write ( command_unit, '(a)' ) &
    'set title "Euler-Maruyama Error as function of DT"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data linespoints'
  write ( command_unit, '(a)' ) &
    'plot "emstrong_data.txt" using 1:2 title "Error", \'
  write ( command_unit, '(a)' ) &
    '     "emstrong_data.txt" using 1:3 title "Slope = 1/2"'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  EMSTRONG plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine emweak ( seed, method, m, p_max, dtvals, xerr )

!*****************************************************************************80
!
!! EMWEAK tests the weak convergence of the Euler-Maruyama method.
!
!  Discussion:
!
!    The SDE is 
!
!      dX = lambda * X dt + mu * X dW,   
!      X(0) = Xzero,
!
!    where 
!
!      lambda = 2,
!      mu = 1,
!      Xzero = 1.
!
!    The discretized Brownian path over [0,1] has dt = 2^(-9).
!
!    The Euler-Maruyama method will use 5 different timesteps:
!
!      2^(p-10),  p = 1,2,3,4,5.
!
!    We examine weak convergence at T=1:
!
!      | E (X_L) - E (X(T)) |.
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
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Input, integer ( kind = 4 ) METHOD.
!    0, use the standard Euler-Maruyama method;
!    1, use the weak Euler-Maruyama method.
!
!    Input, integer ( kind = 4 ) M, the number of simulations to perform.
!    A typical value is M = 1000.
!
!    Input, integer ( kind = 4 ) P_MAX, the number of time step sizes to use.
!    A typical value is 5.
!
!    Output, real ( kind = 8 ) DTVALS(P_MAX), the time steps used.
!
!    Output, real ( kind = 8 ) XERR(P_MAX), the averaged absolute error in the
!    solution estimate at the final time.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) p_max

  real ( kind = 8 ) a(p_max,2)
  real ( kind = 8 ) dt
  real ( kind = 8 ) dt2
  real ( kind = 8 ) dtvals(p_max)
  real ( kind = 8 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) method
  real ( kind = 8 ) mu
  integer ( kind = 4 ) p
  real ( kind = 8 ) q
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_sign
  real ( kind = 8 ) resid
  real ( kind = 8 ) rhs(p_max)
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sol(2)
  real ( kind = 8 ) tmax
  real ( kind = 8 ) winc(m)
  real ( kind = 8 ) xem(p_max)
  real ( kind = 8 ) xerr(p_max)
  real ( kind = 8 ) xtemp(m)
  real ( kind = 8 ) xtrue
  real ( kind = 8 ) xzero
!
!  Problem parameters;
!
  lambda = 2.0D+00
  mu = 0.1D+00
  xzero = 1.0D+00
!
!  Stepping parameters.
!
  tmax = 1.0D+00
  do p = 1, p_max
    dtvals(p) = 2.0D+00 ** ( p - 10 )
  end do
!
!  Take various Euler timesteps.
!  For stepsize dt, we will need to take L Euler steps to reach time TMAX.
!
  do p = 1, p_max

    l = 2 ** ( 10 - p )
    dt = dtvals(p)

    xtemp(1:m) = xzero          

    do j = 1, l
    
      if ( method == 0 ) then
        call r8vec_normal_01 ( m, seed, winc )
        winc(1:m) = sqrt ( dt ) * winc (1:m)
      else
        call r8vec_normal_01 ( m, seed, winc )
        do i = 1, m
          winc(i) = sqrt ( dt ) * r8_sign ( winc(i) )
        end do
      end if

      xtemp(1:m) = xtemp(1:m) + dt * lambda * xtemp(1:m) &
        + mu * xtemp(1:m) * winc

    end do
!
!  Average the M results for this stepsize.
!
    call r8vec_mean ( m, xtemp, xem(p) )

  end do
!
!  Compute the error in the estimates for each stepsize.
!
  xerr(1:p_max) = abs ( xem(1:p_max) - exp ( lambda ) )
!
!  Least squares fit of error = c * dt^q.
!
  a(1:p_max,1) = 1.0D+00
  a(1:p_max,2) = log ( dtvals(1:p_max) )
  rhs(1:p_max) = log ( xerr(1:p_max) )

  call qr_solve ( p_max, 2, a, rhs, sol )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EMWEAK:'
  if ( method == 0 ) then
    write ( *, '(a)' ) '  Using standard Euler-Maruyama method.'
  else
    write ( *, '(a)' ) '  Using weak Euler-Maruyama method.'
  end if
  write ( *, '(a)' ) '  Least squares solution to Error = c * dt ^ q'
  write ( *, '(a)' ) '  (Expecting Q to be about 1.)'
  write ( *, '(a,g14.6)' ) '  Computed Q = ', sol(2)

  resid = 0.0D+00
  do i = 1, p_max
    e = a(i,1) * sol(1) + a(i,2) * sol(2) - rhs(i)
    resid = resid + e * e
  end do
  resid = sqrt ( resid )
  write ( *, '(a,g14.6)' ) '  Residual is ', resid

  return
end
subroutine emweak_gnuplot ( p_max, dtvals, xerr, method )

!*****************************************************************************80
!
!! EMWEAK_GNUPLOT writes a GNUPLOT input file to plot EMWEAK data.
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
!    John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_MAX, the number of time step sizes to use.
!
!    Input, real ( kind = 8 ) DTVALS(P_MAX), the time steps used.
!
!    Input, real ( kind = 8 ) XERR(P_MAX), the averaged absolute error in the
!    solution estimate at the final time.
!
!    Input, integer ( kind = 4 ) METHOD.
!    0, use the standard Euler-Maruyama method;
!    1, use the weak Euler-Maruyama method.
!
  implicit none

  integer ( kind = 4 ) p_max

  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  integer ( kind = 4 ) i
  real ( kind = 8 ) dtvals(p_max)
  integer ( kind = 4 ) method
  real ( kind = 8 ) xerr(p_max)
!
!  Create data file.
!
  call get_unit ( data_unit )

  if ( method == 0 ) then
    data_filename = 'emweak0_data.txt'
  else 
    data_filename = 'emweak1_data.txt'
  end if

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 1, p_max
    write ( data_unit, '(3(2x,g14.6))' ) dtvals(i), xerr(i)
  end do
  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EMWEAK data stored in "' &
    // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  if ( method == 0 ) then
    command_filename = 'emweak0_commands.txt'
  else
    command_filename = 'emweak1_commands.txt'
  end if

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  if ( method == 0 ) then
    write ( command_unit, '(a)' ) '# emweak0_commands.txt'
  else
    write ( command_unit, '(a)' ) '# emweak1_commands.txt'
  end if

  write ( command_unit, '(a)' ) '# created by sde::emweak_gnuplot.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  if ( method == 0 ) then
    write ( command_unit, '(a)' ) '#  gnuplot < emweak0_commands.txt'
  else
    write ( command_unit, '(a)' ) '#  gnuplot < emweak1_commands.txt'
  end if
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  if ( method == 0 ) then
    write ( command_unit, '(a)' ) 'set output "emweak0.png"'
  else
    write ( command_unit, '(a)' ) 'set output "emweak1.png"'
  end if
  write ( command_unit, '(a)' ) 'set xlabel "Log(dt)"'
  write ( command_unit, '(a)' ) 'set ylabel "Log(Averaged Error at final T)"'
  write ( command_unit, '(a)' ) 'set logscale xy 10'
  if ( method == 0 ) then
    write ( command_unit, '(a)' ) &
      'set title "Standard Euler-Maruyama Error as function of DT"'
  else
    write ( command_unit, '(a)' ) &
      'set title "Weak Euler-Maruyama Error as function of DT"'
  end if
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data linespoints'
  if ( method == 0 ) then
    write ( command_unit, '(a)' ) &
      'plot "emweak0_data.txt" using 1:2 title "Error", \'
    write ( command_unit, '(a)' ) &
      '     "emweak0_data.txt" using 1:1 title "Slope = 1"'
  else
    write ( command_unit, '(a)' ) &
      'plot "emweak1_data.txt" using 1:2 title "Error", \'
    write ( command_unit, '(a)' ) &
      '     "emweak1_data.txt" using 1:1 title "Slope = 1"'
  end if

  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  EMWEAK plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine filename_inc ( filename )

!*****************************************************************************80
!
!! FILENAME_INC increments a partially numeric filename.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected.
!
!    If the name is empty, then the routine stops.
!
!    If the name contains no digits, the empty string is returned.
!
!  Example:
!
!      Input            Output
!      -----            ------
!      'a7to11.txt'     'a7to12.txt'
!      'a7to99.txt'     'a8to00.txt'
!      'a9to99.txt'     'a0to00.txt'
!      'cat.txt'        ' '
!      ' '              STOP!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILENAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  integer ( kind = 4 ) change
  integer ( kind = 4 ) digit
  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( filename )

  if ( lens <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILENAME_INC - Fatal error!'
    write ( *, '(a)' ) '  The input string is empty.'
    stop
  end if

  change = 0

  do i = lens, 1, -1

    c = filename(i:i)

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      change = change + 1

      digit = ichar ( c ) - 48
      digit = digit + 1

      if ( digit == 10 ) then
        digit = 0
      end if

      c = char ( digit + 48 )

      filename(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do
!
!  No digits were found.  Return blank.
!
  if ( change == 0 ) then
    filename = ' '
    return
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine milstrong ( seed, p_max, dtvals, xerr )

!*****************************************************************************80
!
!! MILSTRONG tests the strong convergence of the Milstein method.
!
!  Discussion:
!
!    This function solves the stochastic differential equation
!
!      dX = sigma * X * ( k - X ) dt + beta * X dW,  
!      X(0) = Xzero,
!
!    where 
!
!       sigma = 2, 
!       k = 1, 
!       beta = 1,
!       Xzero = 0.5.
!
!    The discretized Brownian path over [0,1] has dt = 2^(-11).
!
!    The Milstein method uses timesteps 128*dt, 64*dt, 32*dt, 16*dt 
!    (also dt for reference).
!
!    We examine strong convergence at T=1:  
!
!      E | X_L - X(T) |.
!
!    The code is vectorized: all paths computed simultaneously.
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
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Input, integer ( kind = 4 ) P_MAX, the number of time step sizes to use.
!    A typical value is 4.
!
!    Output, real ( kind = 8 ) DTVALS(P_MAX), the time steps used.
!
!    Output, real ( kind = 8 ) XERR(P_MAX), the averaged absolute error in the
!    solution estimate at the final time.
!
  implicit none

  integer ( kind = 4 ) p_max

  real ( kind = 8 ) a(p_max,2)
  real ( kind = 8 ) beta
  real ( kind = 8 ) dt
  real ( kind = 8 ) dtp
  real ( kind = 8 ) dtvals(1:p_max)
  real ( kind = 8 ), allocatable :: dw(:,:)
  real ( kind = 8 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  integer ( kind = 4 ) r
  real ( kind = 8 ) resid
  real ( kind = 8 ) rhs(p_max)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sigma
  real ( kind = 8 ) sol(2)
  real ( kind = 8 ) tmax
  real ( kind = 8 ), allocatable :: winc(:)
  real ( kind = 8 ) xerr(p_max)
  real ( kind = 8 ), allocatable :: xref(:)
  real ( kind = 8 ), allocatable :: xtemp(:)
  real ( kind = 8 ) xzero
!
!  Set problem parameters.
!
  sigma = 2.0D+00
  k = 1.0D+00
  beta = 0.25D+00
  xzero = 0.5D+00
!
!  Set stepping parameters.
!
  tmax = 1.0D+00
  n = 2 ** 11
  dt = tmax / real ( n, kind = 8 )
!
!  Number of paths sampled.
!
  m = 500
!
!  Define the increments dW.
!
  allocate ( dw(1:m,1:n) )
  call r8mat_normal_01 ( m, n, seed, dw )
  dw(1:m,1:n) = sqrt ( dt ) * dw(1:m,1:n)
!
!  Estimate the reference solution at time T M times.
!
  allocate ( xref(1:m) )

  xref(1:m) = xzero

  do j = 1, n
    xref(1:m) = xref(1:m) &
      + dt * sigma * xref(1:m) * ( k - xref(1:m) ) &
      + beta * xref(1:m) * dw(1:m,j) &
      + 0.5D+00 * beta ** 2 * xref(1:m) * ( dw(1:m,j) ** 2 - dt )
  end do
!
!  Now compute M Milstein approximations at each of 4 timesteps,
!  and record the average errors.
!
  do p = 1, p_max
    dtvals(p) = dt * 8.0D+00 * 2.0D+00 ** p
  end do

  xerr(1:p_max) = 0.0D+00

  allocate ( winc(1:m) )
  allocate ( xtemp(1:m) )

  do p = 1, p_max

    r = 8 * 2 ** p
    dtp = dtvals(p)
    l = n / r
    xtemp(1:m) = xzero

    do j = 1, l
      winc(1:m) = sum ( dw(1:m,r*(j-1)+1:r*j), 2 )
      xtemp(1:m) = xtemp(1:m) &
        + dtp * sigma * xtemp(1:m) * ( k - xtemp(1:m) ) &
        + beta * xtemp(1:m) * winc(1:m) &
        + 0.5D+00 * beta ** 2 * xtemp(1:m) * ( winc(1:m) ** 2 - dtp )
    end do

    xerr(p) = sum ( abs ( xtemp(1:m) - xref(1:m) ) ) / real ( m, kind = 8 )

  end do
!
!  Least squares fit of error = C * dt^q
!
  a(1:p_max,1) = 1.0D+00
  a(1:p_max,2) = log ( dtvals(1:p_max) )
  rhs(1:p_max) = log ( xerr(1:p_max) )

  call qr_solve ( p_max, 2, a, rhs, sol )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MILSTEIN:'
  write ( *, '(a)' ) '  Least squares solution to Error = c * dt ^ q'
  write ( *, '(a)' ) '  Expecting Q to be about 1.'
  write ( *, '(a,g14.6)' ) '  Computed Q = ', sol(2)

  resid = 0.0D+00
  do i = 1, p_max
    e = a(i,1) * sol(1) + a(i,2) * sol(2) - rhs(i)
    resid = resid + e * e
  end do
  resid = sqrt ( resid )
  write ( *, '(a,g14.6)' ) '  Residual is ', resid

  deallocate ( dw )
  deallocate ( winc )
  deallocate ( xref )
  deallocate ( xtemp )

  return
end
subroutine milstrong_gnuplot ( p_max, dtvals, xerr )

!*****************************************************************************80
!
!! MILSTRONG_GNUPLOT writes a GNUPLOT input file to plot MILSTRONG data.
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
!    John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_MAX, the number of time step sizes to use.
!
!    Input, real ( kind = 8 ) DTVALS(P_MAX), the time steps used.
!
!    Input, real ( kind = 8 ) XERR(P_MAX), the averaged absolute error in the
!    solution estimate at the final time.
!
  implicit none

  integer ( kind = 4 ) p_max

  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  integer ( kind = 4 ) i
  real ( kind = 8 ) dtvals(p_max)
  real ( kind = 8 ) xerr(p_max)
!
!  Create data file.
!
  call get_unit ( data_unit )

  data_filename = 'milstrong_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 1, p_max
    write ( data_unit, '(2(2x,g14.6))' ) dtvals(i), xerr(i)
  end do
  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  MILSTRONG data stored in "' &
    // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = 'milstrong_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# milstrong_commands.txt'
  write ( command_unit, '(a)' ) '# created by sde::milstrong_gnuplot.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < milstrong_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "milstrong.png"'
  write ( command_unit, '(a)' ) 'set xlabel "Log(dt)"'
  write ( command_unit, '(a)' ) 'set ylabel "Log(Averaged Error at final T)"'
  write ( command_unit, '(a)' ) 'set logscale xy 10'
  write ( command_unit, '(a)' ) &
    'set title "Milstein Error as function of DT"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data linespoints'
  write ( command_unit, '(a)' ) &
    'plot "milstrong_data.txt" using 1:2 title "Error", \'
  write ( command_unit, '(a)' ) &
    '     "milstrong_data.txt" using 1:1 title "Slope = 1"'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  MILSTRONG plot commands stored in "' &
    // trim ( command_filename ) // '".'

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
!    This means that on every other call, essentially, the input value of
!    SEED is ignored, since the code saves the second normal random value.
!
!    If you didn't know this, you might be confused since, usually, the
!    output of a random number generator can be completely controlled by
!    the input value of the SEED.  If I were more careful, I could rewrite
!    this routine so that it would distinguish between cases where the input
!    value of SEED is the output value from the previous call (all is well)
!    and those cases where it is not (the user has decided to do something
!    new.  Restart the uniform random number sequence.)  But I'll leave
!    that for later.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a sample of the standard
!    normal PDF.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_01
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: seed2 = 0
  integer ( kind = 4 ), save :: used = 0
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: y = 0.0D+00
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal and its corresponding seed.
!
  if ( mod ( used, 2 ) == 0 ) then

    r1 = r8_uniform_01 ( seed )

    if ( r1 == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    seed2 = seed
    r2 = r8_uniform_01 ( seed2 )

    x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    y = sqrt ( - 2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
!
!  On odd calls, return the second normal and its corresponding seed.
!
  else

    seed = seed2
    x = y

  end if

  used = used + 1

  r8_normal_01 = x

  return
end
function r8_sign ( x )

!*****************************************************************************80
!
!! R8_SIGN returns the sign of an R8.
!
!  Discussion:
!
!    value = -1 if X < 0;
!    value = +1 if X => 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose sign is desired.
!
!    Output, real ( kind = 8 ) R8_SIGN, the sign of X:
!
  implicit none

  real ( kind = 8 ) r8_sign
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    r8_sign = -1.0D+00
  else
    r8_sign = +1.0D+00
  end if

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
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
!    05 July 2006
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
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8mat_normal_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_NORMAL_01 returns a unit pseudonormal R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2010
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudonormal values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  call r8vec_normal_01 ( m * n, seed, r )

  return
end
subroutine r8vec_mean ( n, a, mean )

!*****************************************************************************80
!
!! R8VEC_MEAN returns the mean of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector whose mean is desired.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) mean

  mean = sum ( a(1:n) ) / real ( n, kind = 8 )

  return
end
subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED
!    in order to set the seed of the random number generator.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is
!    negative,then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
!    Local, real ( kind = 8 ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) m
  integer ( kind = 4 ), save :: made = 0
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r(n+1)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ), save :: saved = 0
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
  real ( kind = 8 ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
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
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine stab_asymptotic ( seed, n, p_max )

!*****************************************************************************80
!
!! STAB_ASYMPTOTIC examines asymptotic stability.
!
!  Discussion:
!
!    The function tests the asymptotic stability
!    of the Euler-Maruyama method applied to a stochastic differential
!    equation (SDE).
!
!    The SDE is
!
!      dX = lambda*X dt + mu*X dW,
!      X(0) = Xzero,
!
!    where 
!
!      lambda is a constant,
!      mu is a constant,
!      Xzero = 1.
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
!    Original Matlab version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Input, integer ( kind = 4 ) N, the number of time steps for the
!    first solution.
!
!    Input, integer ( kind = 4 ) P_MAX, the number of time step sizes.
!
  implicit none

  integer ( kind = 4 ) p_max

  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  real ( kind = 8 ) dt
  real ( kind = 8 ) dtvals(p_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lambda
  real ( kind = 8 ) mu
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nval
  integer ( kind = 4 ) p
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) test
  real ( kind = 8 ) tmax
  real ( kind = 8 ) u(1000)
  real ( kind = 8 ) winc
  real ( kind = 8 ), allocatable :: xemabs(:)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xtemp
  real ( kind = 8 ) xzero

  data_filename = 'stab_asymptotic0_data.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STAB_ASYMPTOTIC:'
  write ( *, '(a)' ) '  Investigate asymptotic stability of Euler-Maruyama'
  write ( *, '(a)' ) '  solution with stepsize DT and MU.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SDE is asymptotically stable if'
  write ( *, '(a)' ) '    Real ( lambda - 1/2 mu^2 ) < 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EM with DT is asymptotically stable if'
  write ( *, '(a)' ) '    E log ( | 1 + lambda dt - sqrt(dt) mu n(0,1) | ) < 0.'
  write ( *, '(a)' ) '  where n(0,1) is a normal random value.'
!
!  Problem parameters.
!
  lambda = 0.5D+00
  mu = sqrt ( 6.0D+00 )
  xzero = 1.0D+00
!
!  Test the SDE.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Lambda = ', lambda
  write ( *, '(a,g14.6)' ) '  Mu =     ', mu
  test = lambda - 0.5D+00 * mu ** 2
  write ( *, '(a,g14.6)' ) '  SDE asymptotic stability test = ', test 
!
!  Step parameters.
!
  tmax = 500.0D+00
!
!  For each stepsize, compute the Euler-Maruyama solution.
!
  do p = 1, p_max

    nval = n * 2 ** ( p - 1 )
    dt = tmax / real ( nval, kind = 8 )
    dtvals(p) = dt
!
!  Test the EM for this DT.
! 
    write ( *, '(a)' ) ' '           
    write ( *, '(a,g14.6)' ) '  dt = ', dt
    call r8vec_normal_01 ( 1000, seed, u )
    u(1:1000) = &
      log ( abs ( 1.0D+00 + lambda * dt - sqrt ( dt ) * mu * u(1:1000) ) )
    test = sum ( u(1:1000) ) / 1000.0D+00
    write ( *, '(a,g14.6)' ) '  EM asymptotic test = ', test

    xtemp = xzero
    allocate ( xemabs(0:nval) )
    xemabs(0) = xtemp

    do j = 1, nval
      winc = sqrt ( dt ) * r8_normal_01 ( seed )
      xtemp = xtemp + dt * lambda * xtemp + mu * xtemp * winc
      xemabs(j) = abs ( xtemp )
    end do
!
!  Write this data to a file.
!
    call get_unit ( data_unit )

    call filename_inc ( data_filename )

    open ( unit = data_unit, file = data_filename, status = 'replace' )
!
!  We have to impose a tiny lower bound on the values because we
!  will end up plotting their logs.
!
    xmin = exp ( -200.0D+00 )
    do i = 0, nval
      t = tmax * real ( i, kind = 8 ) / real ( nval, kind = 8 )
      write ( data_unit, '(2x,g14.6,2x,g14.6)' ) t, max ( xemabs(i), xmin )
    end do
    close ( unit = data_unit )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6,a)' ) '  Data for DT = ', dt, ' stored in "' &
      // trim ( data_filename ) // '".'

    deallocate ( xemabs )

  end do
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = 'stab_asymptotic_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# stab_asymptotic_commands.txt'
  write ( command_unit, '(a)' ) '# created by sde::stab_asymptotic.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < stab_asymptotic_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "stab_asymptotic.png"'
  write ( command_unit, '(a)' ) 'set xlabel "t"'
  write ( command_unit, '(a)' ) 'set ylabel "|X(t)|"'
  write ( command_unit, '(a)' ) 'set title "Absolute value of EM Solution"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set logscale y 10'
  write ( command_unit, '(a)' ) 'set style data lines'

  data_filename = 'stab_asymptotic0_data.txt'

  call filename_inc ( data_filename )
  write ( command_unit, '(a)' ) &
    'plot "' // trim ( data_filename ) // '" using 1:2, \'

  do p = 2, p_max - 1
    call filename_inc ( data_filename )
    write ( command_unit, '(a)' ) &
      '     "' // trim ( data_filename ) // '" using 1:2, \'
  end do
  call filename_inc ( data_filename )
  write ( command_unit, '(a)' ) &
    '     "' // trim ( data_filename ) // '" using 1:2'

  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  STAB_ASYMPTOTIC plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine stab_meansquare ( seed )

!*****************************************************************************80
!
!! STAB_MEANSQUARE examines mean-square stability.
!
!  Discussion:
!
!    The function tests the mean-square stability
!    of the Euler-Maruyama method applied to a stochastic differential
!    equation (SDE).
!
!    The SDE is
!
!      dX = lambda*X dt + mu*X dW,
!      X(0) = Xzero,
!
!    where 
!
!      lambda is a constant,
!      mu is a constant,
!      Xzero = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2006
!
!  Author:
!
!    Original Matlab version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED, a seed for the random number generator.
!    In the reference, this value is set to 100.
!
  implicit none

  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  real ( kind = 8 ) dt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) m
  real ( kind = 8 ) mu
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) test
  real ( kind = 8 ) tmax
  real ( kind = 8 ), allocatable :: winc(:)
  real ( kind = 8 ), allocatable :: xms(:)
  real ( kind = 8 ), allocatable :: xtemp(:)
  real ( kind = 8 ) xzero

  data_filename = 'stab_meansquare0_data.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STAB_MEANSQUARE:'
  write ( *, '(a)' ) '  Investigate mean square stability of Euler-Maruyama'
  write ( *, '(a)' ) '  solution with stepsize DT and MU.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SDE is mean square stable if'
  write ( *, '(a)' ) '    Real ( lambda + 1/2 |mu|^2 ) < 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EM with DT is mean square stable if'
  write ( *, '(a)' ) '    |1+dt^2| + dt * |mu|^2 - 1.0 < 0.'
!
!  Set problem parameters.
!
  tmax = 20.0D+00
  m = 50000
  xzero = 1.0D+00 
!
!  Problem parameters.
!
  lambda = -3.0D+00
  mu = sqrt ( 3.0D+00 )
!
!  Test the SDE.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Lambda = ', lambda
  write ( *, '(a,g14.6)' ) '  Mu =     ', mu
  test = lambda + 0.5D+00 * mu ** 2
  write ( *, '(a,g14.6)' ) '  SDE mean square stability test = ', test
!
!  XMS is the mean square estimate of M paths.
!
  allocate ( winc(1:m) )

  do k = 1, 3

    dt = 2.0D+00 ** ( 1 - k )                      
    n = 20 * 2 ** ( k - 1 )
!
!  Test the EM for this DT.
! 
    write ( *, '(a)' ) ' '           
    write ( *, '(a,g14.6)' ) '  dt = ', dt
    test = ( 1.0D+00 + dt * lambda ) ** 2 + dt * mu ** 2 - 1.0D+00
    write ( *, '(a,g14.6)' ) '  EM mean square stability test = ', test

    allocate ( xms(0:n) )
    allocate ( xtemp(1:m) )
    xtemp(1:m) = xzero
    xms(0) = xzero

    do j = 1, n
      call r8vec_normal_01 ( m, seed, winc )
      winc(1:m) = sqrt ( dt ) * winc(1:m) 
      xtemp(1:m) = xtemp(1:m) &
        + dt * lambda * xtemp(1:m) &
        + mu * xtemp(1:m) * winc(1:m)
      xms(j) = sum ( xtemp(1:m) ** 2 ) / real ( m, kind = 8 )
    end do
!
!  Write this data to a file.
!
    call get_unit ( data_unit )

    call filename_inc ( data_filename )

    open ( unit = data_unit, file = data_filename, status = 'replace' )

    do j = 0, n
      t = tmax * real ( j, kind = 8 ) / real ( n, kind = 8 )
      write ( data_unit, '(2x,g14.6,2x,g14.6)' ) t, xms(j)
    end do
    close ( unit = data_unit )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6,a)' ) '  Data for DT = ', dt, ' stored in "' &
      // trim ( data_filename ) // '".'

    deallocate ( xtemp )
    deallocate ( xms )

  end do

  deallocate ( winc )
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = 'stab_meansquare_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# stab_meansquare_commands.txt'
  write ( command_unit, '(a)' ) '# created by sde::stab_meansquare.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < stab_meansquare_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "stab_meansquare.png"'
  write ( command_unit, '(a)' ) 'set xlabel "t"'
  write ( command_unit, '(a)' ) 'set ylabel "E|X^2(t)|"'
  write ( command_unit, '(a)' ) 'set title "Mean Square of EM Solution"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set logscale y 10'
  write ( command_unit, '(a)' ) 'set style data lines'

  data_filename = 'stab_meansquare0_data.txt'

  call filename_inc ( data_filename )
  write ( command_unit, '(a)' ) &
    'plot "' // trim ( data_filename ) // '" using 1:2, \'

  do k = 2, 2
    call filename_inc ( data_filename )
    write ( command_unit, '(a)' ) &
      '     "' // trim ( data_filename ) // '" using 1:2, \'
  end do
  call filename_inc ( data_filename )
  write ( command_unit, '(a)' ) &
    '     "' // trim ( data_filename ) // '" using 1:2'

  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  STAB_MEANSQUARE plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine stochastic_integral_ito ( n, seed, estimate, exact, error )

!*****************************************************************************80
!
!! STOCHASTIC_INTEGRAL_ITO approximates the Ito integral of W(t) dW.
!
!  Discussion:
!
!    This function estimates the Ito integral of W(t) dW over 
!    the interval [0,1].
!
!    The estimates is made by taking N steps.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2012
!
!  Author:
!
!    Original Matlab version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of steps to take.
!
!    Input, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) ESTIMATE, the estimate of the integral.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
!    Output, real ( kind = 8 ) ERROR, the error in the integral estimate.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dt
  real ( kind = 8 ) dw(n)
  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tmax
  real ( kind = 8 ) w(0:n)
!
!  Set step parameters.
!
  tmax = 1.0D+00
  dt = tmax / real ( n, kind = 8 )
!
!  Define the increments dW.
!
  call r8vec_normal_01 ( n, seed, dw )

  dw(1:n) = sqrt ( dt ) * dw(1:n)
!
!  Sum the increments to get the Brownian path.
!
  w(0) = 0.0D+00
  do j = 1, n
    w(j) = w(j-1) + dw(j)
  end do
!
!  Approximate the Ito integral.
!
  estimate = dot_product ( w(0:n-1), dw(1:n) )
!
!  Compare with the exact solution.
!
  exact = 0.5D+00 * ( w(n) ** 2 - tmax )
  error = abs ( estimate - exact )

  return
end
subroutine stochastic_integral_strat ( n, seed, estimate, exact, error )

!*****************************************************************************80
!
!! STOCHASTIC_INTEGRAL_STRAT approximates the Stratonovich integral of W(t) dW.
!
!  Discussion:
!
!    This function estimates the Stratonovich integral of W(t) dW over 
!    the interval [0,1].
!
!    The estimates is made by taking N steps.
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
!    Original Matlab version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of steps to take.
!
!    Input, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) ESTIMATE, the estimate of the integral.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
!    Output, real ( kind = 8 ) ERROR, the error in the integral estimate.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dt
  real ( kind = 8 ) dw(n)
  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exact
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tmax
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) w(0:n)
!
!  Set step parameters.
!
  tmax = 1.0D+00
  dt = tmax / real ( n, kind = 8 )
!
!  Define the increments dW.
!
  call r8vec_normal_01 ( n, seed, dw )

  dw(1:n) = sqrt ( dt ) * dw(1:n)
!
!  Sum the increments to get the Brownian path.
!
  w(0) = 0.0D+00
  do j = 1, n
    w(j) = w(j-1) + dw(j)
  end do
!
!  Approximate the Stratonovich integral.
!
  call r8vec_normal_01 ( n, seed, u )
  v(1:n) = 0.5D+00 * ( w(0:n-1) + w(1:n) ) + 0.5D+00 * sqrt ( dt ) * u(1:n)

  estimate = dot_product ( v, dw )
!
!  Compare with the exact solution.
!
  exact = 0.5D+00 * w(n) ** 2
  error = abs ( estimate - exact )

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
