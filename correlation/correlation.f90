subroutine correlation_besselj ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_BESSELJ evaluates the Bessel J correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_besj0
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhohat

  do i = 1, n
    rhohat = abs ( rho(i) ) / rho0
    c(i) = r8_besj0 ( rhohat )
  end do

  return
end
subroutine correlation_besselk ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_BESSELK evaluates the Bessel K correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_besk1
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhohat

  do i = 1, n

    if ( rho(i) == 0.0D+00 ) then
      c(i) = 1.0D+00
    else
      rhohat = abs ( rho(i) ) / rho0
      c(i) = rhohat * r8_besk1 ( rhohat )
    end if

  end do

  return
end
subroutine correlation_brownian ( m, n, s, t, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_BROWNIAN computes the Brownian correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of arguments.
!
!    Input, real ( kind = 8 ) S(M), T(N), two samples.
!    0 <= S(*), T(*).
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(M,N), the correlations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) rho0
  real ( kind = 8 ) s(m)
  real ( kind = 8 ) t(n)

  do j = 1, n
    do i = 1, m
      if ( 0.0D+00 < max ( s(i), t(j) ) ) then
        c(i,j) = sqrt ( min ( s(i), t(j) ) / max ( s(i), t(j) ) )
      else
        c(i,j) = 1.0D+00
      end if
    end do
  end do

  return
end
subroutine correlation_brownian_display ( )

!*****************************************************************************80
!
!! CORRELATION_BROWNIAN_DISPLAY displays 4 slices of the Brownian Correlation.
!
!  Discussion:
!
!    The correlation function is C(S,T) = sqrt ( min ( s, t ) / max ( s, t ) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 101
  integer ( kind = 4 ), parameter :: n2 = 4

  real ( kind = 8 ) c(n,n2)
  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) t(n2)

  call r8vec_linspace ( n, 0.0D+00, 5.0D+00, s )

  t(1:n2) = (/ 0.25D+00, 1.50D+00, 2.50D+00, 3.75D+00 /)

  do i = 1, n
    do j = 1, n2
      c(i,j) = sqrt ( min ( s(i), t(j) ) / max ( s(i), t(j) ) )
    end do
  end do

  call get_unit ( data_unit )
  data_filename = 'brownian_plots_data.txt'
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 1, n
    write ( data_unit, '(5(2x,g14.6))' ) s(i), c(i,1:n2)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) '  Created data file "' // trim ( data_filename ) // '".'

  call get_unit ( command_unit )
  command_filename = 'brownian_plots_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set key off'
  write ( command_unit, '(a)' ) &
    'set output "brownian_plots.png"'
  write ( command_unit, '(a)' ) &
    'set title "Brownian correlation C(S,T), S = 0.25, 1.5, 2.5, 3.75"'
  write ( command_unit, '(a)' ) 'set xlabel "S"'
  write ( command_unit, '(a)' ) 'set ylabel "C(s,t)"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 lw 3 linecolor rgb "blue", \'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:3 lw 3 linecolor rgb "blue", \'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:4 lw 3 linecolor rgb "blue", \'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:5 lw 3 linecolor rgb "blue"'
  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
end
subroutine correlation_circular ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_CIRCULAR evaluates the circular correlation function.
!
!  Discussion:
!
!    This correlation is based on the area of overlap of two circles
!    of radius RHO0 and separation RHO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhohat

  do i = 1, n

    rhohat = min ( abs ( rho(i) ) / rho0, 1.0D+00 )

    c(i) = ( 1.0D+00 - ( 2.0D+00 / pi ) &
      * ( rhohat * sqrt ( 1.0D+00 - rhohat ** 2 ) + asin ( rhohat ) ) )

  end do

  return
end
subroutine correlation_constant ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_CONSTANT evaluates the constant correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0

  c(1:n) = 1.0D+00

  return
end
subroutine correlation_cubic ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_CUBIC evaluates the cubic correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhohat

  do i = 1, n

    rhohat = min ( abs ( rho(i) ) / rho0, 1.0D+00 )

    c(i) = 1.0D+00 &
         - 7.0D+00  * rhohat ** 2 &
         + 8.75D+00 * rhohat ** 3 &
         - 3.5D+00  * rhohat ** 5 &
         + 0.75D+00 * rhohat ** 7

  end do

  return
end
subroutine correlation_damped_cosine ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_DAMPED_COSINE evaluates the damped cosine correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0

  c(1:n) = exp ( - abs ( rho(1:n) ) / rho0 ) * cos ( abs ( rho(1:n) ) / rho0 )

  return
end
subroutine correlation_damped_sine ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_DAMPED_SINE evaluates the damped sine correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhohat

  do i = 1, n

    if ( rho(i) == 0.0D+00 ) then
      c(i) = 1.0D+00
    else
      rhohat = abs ( rho(i) ) / rho0
      c(i) = sin ( rhohat ) / rhohat
    end if

  end do

  return
end
subroutine correlation_exponential ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_EXPONENTIAL evaluates the exponential correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0

  c(1:n) = exp ( - abs ( rho(1:n) ) / rho0 )

  return
end
subroutine correlation_gaussian ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_GAUSSIAN evaluates the Gaussian correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0

  c(1:n) = exp ( - ( ( rho(1:n) / rho0 ) ** 2 ) )

  return
end
subroutine correlation_hole ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_HOLE evaluates the hole correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0

  c(1:n) = ( 1.0D+00 - abs ( rho(1:n) ) / rho0 ) &
    * exp ( - abs ( rho(1:n) ) / rho0 )

  return
end
subroutine correlation_linear ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_LINEAR evaluates the linear correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0

  do i = 1, n
    if ( rho0 < abs ( rho(i) ) ) then
      c(i) = 0.0D+00
    else
      c(i) = ( rho0 - abs ( rho(i) ) ) / rho0
    end if
  end do

  return
end
subroutine correlation_matern ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_MATERN evaluates the Matern correlation function.
!
!  Discussion:
!
!    In order to call this routine under a dummy name, I had to drop NU from
!    the parameter list.
!
!    The Matern correlation is
!
!      rho1 = 2 * sqrt ( nu ) * rho / rho0
!
!      c(rho) = ( rho1 )^nu * BesselK ( nu, rho1 ) 
!               / gamma ( nu ) / 2 ^ ( nu - 1 )
!
!    The Matern covariance has the form:
!
!      K(rho) = sigma^2 * c(rho)
!
!    A Gaussian process with Matern covariance has sample paths that are
!    differentiable (nu - 1) times.
!
!    When nu = 0.5, the Matern covariance is the exponential covariance.
!
!    As nu goes to +oo, the correlation converges to exp ( - (rho/rho0)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!    0.0 <= RHO.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!    0.0 < RHO0.
!
!    Input, real ( kind = 8 ) NU, the smoothness parameter.
!    NU has a default value of 2.5;
!    0 < NU.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) nu
  real ( kind = 8 ) r8_besk
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rho1

  nu = 2.5D+00

  do i = 1, n

    rho1 = 2.0D+00 * sqrt ( nu ) * abs ( rho(i) ) / rho0

    if ( rho1 == 0.0D+00 ) then
      c(i) = 1.0D+00
    else
      c(i) = rho1 ** nu * r8_besk ( nu, rho1 ) / r8_gamma ( nu ) &
        / 2.0 ** ( nu - 1.0D+00 )
    end if

  end do

  return
end
subroutine correlation_pentaspherical ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_PENTASPHERICAL evaluates the pentaspherical correlation function.
!
!  Discussion:
!
!    This correlation is based on the volume of overlap of two spheres
!    of radius RHO0 and separation RHO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhohat

  do i = 1, n

    rhohat = min ( abs ( rho(i) ) / rho0, 1.0 )

    c(i) = 1.0D+00 - 1.875D+00 * rhohat + 1.25D+00 * rhohat ** 3 &
      - 0.375D+00 * rhohat ** 5

  end do

  return
end
subroutine correlation_plot ( n, rho, c, header, title )

!*****************************************************************************80
!
!! CORRELATION_PLOT makes a plot of a correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) C(N), the correlations.
!
!    Input, character ( len = * ) HEADER, an identifier for the files.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  character ( len = * ) header
  integer ( kind = 4 ) i
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0
  character ( len = * ) title

  call get_unit ( data_unit )
  data_filename = trim ( header ) // '_data.txt'
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 1, n
    write ( data_unit, '(2x,g14.6,2x,g14.6)' ) rho(i), c(i)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) '  Created data file "' // trim ( data_filename ) // '".'

  call get_unit ( command_unit )
  command_filename = trim ( header ) // '_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) &
    'set output "' // trim ( header ) // '_plot.png"'
  write ( command_unit, '(a)' ) 'set xlabel "Distance Rho"'
  write ( command_unit, '(a)' ) 'set ylabel "Correlation C(Rho)"'
  write ( command_unit, '(a)' ) 'set title "' // trim ( title ) // '"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 lw 3 linecolor rgb "blue"'
  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
end
subroutine correlation_plots ( n, n2, rho, rho0, c, header, title )

!*****************************************************************************80
!
!! CORRELATION_PLOTS plots correlations for a range of correlation lengths.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values of RHO.
!
!    Input, integer ( kind = 4 ) N2, the number of values of RHO0.
!
!    Input, real ( kind = 8 ) RHO(N), the independent value.
!
!    Input, real ( kind = 8 ) RHO0(N2), the correlation lengths.
!
!    Input, real ( kind = 8 ) C(N,N2), the correlations.
!
!    Input, character ( len = * ) HEADER, an identifier for the files.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2

  real ( kind = 8 ) c(n,n2)
  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  character ( len = 40 ) format_string
  character ( len = * ) header
  integer ( kind = 4 ) i
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0(n2)
  character ( len = * ) title

  write ( format_string, '(a1,i8,a1,i8,a1,i8,a1)' ) &
    '(', n2 + 1, 'g', 14, '.', 6, ')'

  call get_unit ( data_unit )
  data_filename = trim ( header ) // '_plots_data.txt'
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 1, n
    write ( data_unit, format_string ) rho(i), c(i,1:n2)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) '  Created data file "' // trim ( data_filename ) // '".'

  call get_unit ( command_unit )
  command_filename = trim ( header ) // '_plots_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) &
    'set output "' // trim ( header ) // '_plots.png"'
  write ( command_unit, '(a)' ) 'set xlabel "Rho"'
  write ( command_unit, '(a)' ) 'set ylabel "Correlation(Rho)"'
  write ( command_unit, '(a)' ) 'set title "' // trim ( title ) // '"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'set key off'
  if ( n2 == 1 ) then
    write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 lw 3'
  else
    write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 lw 3, \'
    do i = 2, n2 - 1
      write ( command_unit, '(a,i3,a)' ) '     "' // trim ( data_filename ) // &
      '" using 1:', i + 1, ' lw 3, \'
    end do
    write ( command_unit, '(a,i3,a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:', n2 + 1, ' lw 3'
  end if
  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
end
subroutine correlation_power ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_POWER evaluates the power correlation function.
!
!  Discussion:
!
!    In order to be able to call this routine under a dummy name, I had
!    to drop E from the argument list.
!
!    The power correlation is
!
!      C(rho) = ( 1 - |rho| )^e  if 0 <= |rho| <= 1
!             = 0                otherwise
!
!      The constraint on the exponent is 2 <= e.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!    0.0 <= RHO.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!    0.0 < RHO0.
!
!    Input, real ( kind = 8 ) E, the exponent.
!    E has a default value of 2.0;
!    2.0 <= E.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) e
  integer ( kind = 4 ) i
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhohat

  e = 2.0D+00

  do i = 1, n
    rhohat = abs ( rho(i) ) / rho0
    if ( rhohat <= 1.0D+00 ) then
      c(i) = ( 1.0D+00 - rhohat ) ** e
    else
      c(i) = 0.0D+00
    end if
  end do

  return
end
subroutine correlation_rational_quadratic ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_RATIONAL_QUADRATIC: rational quadratic correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0

  c(1:n) = 1.0D+00 / ( 1.0D+00 + ( rho(1:n) / rho0 ) ** 2 )

  return
end
subroutine correlation_spherical ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_SPHERICAL evaluates the spherical correlation function.
!
!  Discussion:
!
!    This correlation is based on the volume of overlap of two spheres
!    of radius RHO0 and separation RHO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhohat

  do i = 1, n
    rhohat = min ( abs ( rho(i) ) / rho0, 1.0D+00 )
    c(i) = 1.0D+00 - 1.5D+00 * rhohat + 0.5D+00 * rhohat ** 3
  end do

  return
end
subroutine correlation_to_covariance ( n, c, sigma, k )

!*****************************************************************************80
!
!! CORRELATION_TO_COVARIANCE: covariance matrix from a correlation matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) C(N,N), the correlation matrix.
!
!    Input, real ( kind = 8 ) SIGMA(N), the standard deviations.
!
!    Output, real ( kind = 8 ) K(N,N), the covariance matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n,n)
  real ( kind = 8 ) c_max
  real ( kind = 8 ) c_min
  real ( kind = 8 ) e
  real ( kind = 8 ) error_frobenius
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) k(n,n)
  real ( kind = 8 ) sigma(n)
  real ( kind = 8 ) tol

  tol = sqrt ( epsilon ( tol ) )
!
!  C must be symmetric.
!
  call r8mat_is_symmetric ( n, n, c, error_frobenius )

  if ( tol < error_frobenius ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CORRELATION_TO_COVARIANCE - Fatal error!'
    write ( *, '(a,g14.6)' ) &
      '  Input matrix C fails symmetry test with error ', error_frobenius
    stop
  end if
!
!  The diagonal must be 1.
!
  do i = 1, n
    e = abs ( c(i,i) - 1.0D+00 )
    if ( tol < e ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'CORRELATION_TO_COVARIANCE - Fatal error!'
      write ( *, '(a)' ) '  Input matrix C has non-unit diagonal entries.'
      write ( *, '(a,i4,a,g14.6)' ) '  Error on row ', i, ' is ', e
      stop
    end if
  end do
!
!  Off-diagonals must be between -1 and 1.
!
  c_min = minval ( c(1:n,1:n) )

  if ( c_min < - 1.0D+00 - tol ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CORRELATION_TO_COVARIANCE - Fatal error!'
    write ( *, '(a)' ) '  Input matrix C has entries less than -1.0'
    stop
  end if

  c_max = maxval ( c(1:n,1:n) )

  if ( 1.0D+00 + tol < c_max ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CORRELATION_TO_COVARIANCE - Fatal error!'
    write ( *, '(a)' ) '  Input matrix C has entries greater than +1.0'
    stop
  end if
!
!  Form K.
!
  do j = 1, n
    do i = 1, n
      k(i,j) = sigma(i) * c(i,j) * sigma(j)
    end do
  end do

  return
end
subroutine correlation_white_noise ( n, rho, rho0, c )

!*****************************************************************************80
!
!! CORRELATION_WHITE_NOISE evaluates the white noise correlation function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Petter Abrahamsen,
!    A Review of Gaussian Random Fields and Correlation Functions,
!    Norwegian Computing Center, 1997.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) RHO(N), the arguments.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Output, real ( kind = 8 ) C(N), the correlations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0

  do i = 1, n
    if ( rho(i) == 0.0D+00 ) then
      c(i) = 1.0D+00
    else
      c(i) = 0.0D+00
    end if
  end do

  return
end
subroutine covariance_to_correlation ( n, k, c, sigma )

!*****************************************************************************80
!
!! COVARIANCE_TO_CORRELATION: correlation matrix from a covariance matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) K(N,N), the covariance matrix.
!
!    Output, real ( kind = 8 ) C(N,N), the correlation matrix.
!
!    Output, real ( kind = 8 ) SIGMA(N), the standard deviations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n,n)
  real ( kind = 8 ) e
  real ( kind = 8 ) error_frobenius
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) k(n,n)
  real ( kind = 8 ) sigma(n)
  real ( kind = 8 ) sigma_min
  real ( kind = 8 ) tol

  tol = sqrt ( epsilon ( tol ) )
!
!  K must be symmetric.
!
  call r8mat_is_symmetric ( n, n, k, error_frobenius )

  if ( tol < error_frobenius ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COVARIANCE_TO_CORRELATION - Fatal error!'
    write ( *, '(a,g14.6)' ) &
      '  Input matrix K fails symmetry test with error ', error_frobenius
    stop
  end if
!
!  It must be the case that K(I,J)^2 <= K(I,I) * K(J,J).
!
  e = 0.0D+00
  do i = 1, n
    do j = i + 1, n
      e = max ( e, k(i,j) ** 2 - k(i,i) * k(j,j) )
    end do
  end do

  if ( tol < e ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COVARIANCE_TO_CORRELATION - Fatal error!'
    write ( *, '(a)' ) '  Input matrix K fails K(I,J)^2 <= K(I,I)*K(J,J)'
    stop
  end if
!
!  Get the diagonal.
!
  do i = 1, n
    sigma(i) = k(i,i)
  end do
!
!  Ensure the diagonal is positive.
!
  sigma_min = minval ( sigma(1:n) )

  if ( sigma_min <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COVARIANCE_TO_CORRELATION - Fatal error!'
    write ( *, '(a,g14.6)' ) &
      '  Input matrix K has nonpositive diagonal entry = ', sigma_min
    stop
  end if
!
!  Convert from variance to standard deviation.
!
  sigma(1:n) = sqrt ( sigma(1:n) )
!
!  Form C.
!
  do j = 1, n
    do i = 1, n
      c(i,j) = k(i,j) / sigma(i) / sigma(j)
    end do
  end do

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
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, a value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) value
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
subroutine minij ( m, n, a )

!*****************************************************************************80
!
!! MINIJ returns the MINIJ matrix.
!
!  Formula:
!
!    A(I,J) = min ( I, J )
!
!  Example:
!
!    N = 5
!
!    1 1 1 1 1
!    1 2 2 2 2
!    1 2 3 3 3
!    1 2 3 4 4
!    1 2 3 4 5
!
!  Properties:
!
!    A is integral, therefore det ( A ) is integral, and 
!    det ( A ) * inverse ( A ) is integral.
!
!    A is positive definite.
!
!    A is symmetric: A' = A.
!
!    Because A is symmetric, it is normal.
!
!    Because A is normal, it is diagonalizable.
!
!    The inverse of A is tridiagonal.
!
!    The eigenvalues of A are
!
!      LAMBDA(I) = 0.5 / ( 1 - cos ( ( 2 * I - 1 ) * pi / ( 2 * N + 1 ) ) ),
!
!    (N+1)*ONES(N) - A also has a tridiagonal inverse.
!
!    Gregory and Karney consider the matrix defined by
!
!      B(I,J) = N + 1 - MAX(I,J)
!
!    which is equal to the MINIJ matrix, but with the rows and
!    columns reversed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Gregory, David Karney,
!    Example 3.12, Example 4.14,
!    A Collection of Matrices for Testing Computational Algorithms,
!    Wiley, 1969, page 41, page 74, 
!    LC: QA263.G68.
!
!    Daniel Rutherford,
!    Some continuant determinants arising in physics and chemistry II,
!    Proceedings of the Royal Society Edinburgh,
!    Volume 63, A, 1952, pages 232-241.
!
!    John Todd,
!    Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
!    Academic Press, 1977, page 158.
!
!    Joan Westlake,
!    A Handbook of Numerical Matrix Inversion and Solution of 
!    Linear Equations,
!    John Wiley, 1968,
!    ISBN13: 978-0471936756,
!    LC: QA263.W47.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n
    do i = 1, m
      a(i,j) = real ( min ( i, j ), kind = 8 )
    end do
  end do

  return
end
subroutine paths_plot ( n, n2, rho, x, header, title )

!*****************************************************************************80
!
!! PATHS_PLOT plots a sequence of paths or simulations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points in each path.
!
!    Input, integer ( kind = 4 ) N2, the number of paths.
!
!    Input, real ( kind = 8 ) RHO(N), the independent value.
!
!    Input, real ( kind = 8 ) X(N,N2), the path values.
!
!    Input, character ( len = * ) HEADER, an identifier for the files.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2

  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  character ( len = 40 ) format_string
  character ( len = * ) header
  integer ( kind = 4 ) i
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) rho0
  character ( len = * ) title
  real ( kind = 8 ) x(n,n2)

  write ( format_string, '(a1,i8,a1,i8,a1,i8,a1)' ) &
    '(', n2 + 1, 'g', 14, '.', 6, ')'

  call get_unit ( data_unit )
  data_filename = trim ( header ) // '_path_data.txt'
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 1, n
    write ( data_unit, format_string ) rho(i), x(i,1:n2)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) '  Created data file "' // trim ( data_filename ) // '".'

  call get_unit ( command_unit )
  command_filename = trim ( header ) // '_path_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) &
    'set output "' // trim ( header ) // '_paths.png"'
  write ( command_unit, '(a)' ) 'set xlabel "Rho"'
  write ( command_unit, '(a)' ) 'set ylabel "X(Rho)"'
  write ( command_unit, '(a)' ) 'set title "' // trim ( title ) // '"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'set key off'
  if ( n2 == 1 ) then
    write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 lw 3'
  else
    write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2, \'
    do i = 2, n2 - 1
      write ( command_unit, '(a,i3,a)' ) '     "' // trim ( data_filename ) // &
      '" using 1:', i + 1, ', \'
    end do
    write ( command_unit, '(a,i3)' ) '     "' // trim ( data_filename ) // &
    '" using 1:', n2 + 1
  end if
  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
end
function pythag ( a, b )

!*****************************************************************************80
!
!! PYTHAG computes sqrt ( A * A + B * B ) carefully.
!
!  Discussion:
!
!    The formula
!
!      PYTHAG = sqrt ( A * A + B * B )
!
!    is reasonably accurate, but can fail if, for example, A^2 is larger
!    than the machine overflow.  The formula can lose most of its accuracy
!    if the sum of the squares is very large or very small.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Modified:
!
!    04 February 2003
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the two legs of a right triangle.
!
!    Output, real ( kind = 8 ) PYTHAG, the length of the hypotenuse.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) p
  real ( kind = 8 ) pythag
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) u

  p = max ( abs ( a ), abs ( b ) )

  if ( p /= 0.0D+00 ) then

    r = ( min ( abs ( a ), abs ( b ) ) / p ) ** 2

    do

      t = 4.0D+00 + r

      if ( t == 4.0D+00 ) then
        exit
      end if

      s = r / t
      u = 1.0D+00 + 2.0D+00 * s
      p = u * p
      r = ( s / u ) ** 2 * r

    end do

  end if

  pythag = p

  return
end
subroutine r8_b0mp ( x, ampl, theta )

!*****************************************************************************80
!
!! R8_B0MP evaluates the modulus and phase for the Bessel J0 and Y0 functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) AMPL, THETA, the modulus and phase.
!
  implicit none

  real ( kind = 8 ) ampl
  real ( kind = 8 ) bm0cs(37)
  real ( kind = 8 ) bm02cs(40)
  real ( kind = 8 ) bt02cs(39)
  real ( kind = 8 ) bth0cs(44)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) nbm0
  integer ( kind = 4 ) nbm02
  integer ( kind = 4 ) nbt02
  integer ( kind = 4 ) nbth0
  real ( kind = 8 ) pi4
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) theta
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) z

  save bm0cs
  save bm02cs
  save bt02cs
  save bth0cs
  save nbm0
  save nbm02
  save nbt02
  save nbth0
  save xmax

  data bm0cs(  1) / +0.9211656246827742712573767730182D-01/
  data bm0cs(  2) / -0.1050590997271905102480716371755D-02/
  data bm0cs(  3) / +0.1470159840768759754056392850952D-04/
  data bm0cs(  4) / -0.5058557606038554223347929327702D-06/
  data bm0cs(  5) / +0.2787254538632444176630356137881D-07/
  data bm0cs(  6) / -0.2062363611780914802618841018973D-08/
  data bm0cs(  7) / +0.1870214313138879675138172596261D-09/
  data bm0cs(  8) / -0.1969330971135636200241730777825D-10/
  data bm0cs(  9) / +0.2325973793999275444012508818052D-11/
  data bm0cs( 10) / -0.3009520344938250272851224734482D-12/
  data bm0cs( 11) / +0.4194521333850669181471206768646D-13/
  data bm0cs( 12) / -0.6219449312188445825973267429564D-14/
  data bm0cs( 13) / +0.9718260411336068469601765885269D-15/
  data bm0cs( 14) / -0.1588478585701075207366635966937D-15/
  data bm0cs( 15) / +0.2700072193671308890086217324458D-16/
  data bm0cs( 16) / -0.4750092365234008992477504786773D-17/
  data bm0cs( 17) / +0.8615128162604370873191703746560D-18/
  data bm0cs( 18) / -0.1605608686956144815745602703359D-18/
  data bm0cs( 19) / +0.3066513987314482975188539801599D-19/
  data bm0cs( 20) / -0.5987764223193956430696505617066D-20/
  data bm0cs( 21) / +0.1192971253748248306489069841066D-20/
  data bm0cs( 22) / -0.2420969142044805489484682581333D-21/
  data bm0cs( 23) / +0.4996751760510616453371002879999D-22/
  data bm0cs( 24) / -0.1047493639351158510095040511999D-22/
  data bm0cs( 25) / +0.2227786843797468101048183466666D-23/
  data bm0cs( 26) / -0.4801813239398162862370542933333D-24/
  data bm0cs( 27) / +0.1047962723470959956476996266666D-24/
  data bm0cs( 28) / -0.2313858165678615325101260800000D-25/
  data bm0cs( 29) / +0.5164823088462674211635199999999D-26/
  data bm0cs( 30) / -0.1164691191850065389525401599999D-26/
  data bm0cs( 31) / +0.2651788486043319282958336000000D-27/
  data bm0cs( 32) / -0.6092559503825728497691306666666D-28/
  data bm0cs( 33) / +0.1411804686144259308038826666666D-28/
  data bm0cs( 34) / -0.3298094961231737245750613333333D-29/
  data bm0cs( 35) / +0.7763931143074065031714133333333D-30/
  data bm0cs( 36) / -0.1841031343661458478421333333333D-30/
  data bm0cs( 37) / +0.4395880138594310737100799999999D-31/

  data bth0cs(  1) / -0.24901780862128936717709793789967D+00/
  data bth0cs(  2) / +0.48550299609623749241048615535485D-03/
  data bth0cs(  3) / -0.54511837345017204950656273563505D-05/
  data bth0cs(  4) / +0.13558673059405964054377445929903D-06/
  data bth0cs(  5) / -0.55691398902227626227583218414920D-08/
  data bth0cs(  6) / +0.32609031824994335304004205719468D-09/
  data bth0cs(  7) / -0.24918807862461341125237903877993D-10/
  data bth0cs(  8) / +0.23449377420882520554352413564891D-11/
  data bth0cs(  9) / -0.26096534444310387762177574766136D-12/
  data bth0cs( 10) / +0.33353140420097395105869955014923D-13/
  data bth0cs( 11) / -0.47890000440572684646750770557409D-14/
  data bth0cs( 12) / +0.75956178436192215972642568545248D-15/
  data bth0cs( 13) / -0.13131556016891440382773397487633D-15/
  data bth0cs( 14) / +0.24483618345240857495426820738355D-16/
  data bth0cs( 15) / -0.48805729810618777683256761918331D-17/
  data bth0cs( 16) / +0.10327285029786316149223756361204D-17/
  data bth0cs( 17) / -0.23057633815057217157004744527025D-18/
  data bth0cs( 18) / +0.54044443001892693993017108483765D-19/
  data bth0cs( 19) / -0.13240695194366572724155032882385D-19/
  data bth0cs( 20) / +0.33780795621371970203424792124722D-20/
  data bth0cs( 21) / -0.89457629157111779003026926292299D-21/
  data bth0cs( 22) / +0.24519906889219317090899908651405D-21/
  data bth0cs( 23) / -0.69388422876866318680139933157657D-22/
  data bth0cs( 24) / +0.20228278714890138392946303337791D-22/
  data bth0cs( 25) / -0.60628500002335483105794195371764D-23/
  data bth0cs( 26) / +0.18649748964037635381823788396270D-23/
  data bth0cs( 27) / -0.58783732384849894560245036530867D-24/
  data bth0cs( 28) / +0.18958591447999563485531179503513D-24/
  data bth0cs( 29) / -0.62481979372258858959291620728565D-25/
  data bth0cs( 30) / +0.21017901684551024686638633529074D-25/
  data bth0cs( 31) / -0.72084300935209253690813933992446D-26/
  data bth0cs( 32) / +0.25181363892474240867156405976746D-26/
  data bth0cs( 33) / -0.89518042258785778806143945953643D-27/
  data bth0cs( 34) / +0.32357237479762298533256235868587D-27/
  data bth0cs( 35) / -0.11883010519855353657047144113796D-27/
  data bth0cs( 36) / +0.44306286907358104820579231941731D-28/
  data bth0cs( 37) / -0.16761009648834829495792010135681D-28/
  data bth0cs( 38) / +0.64292946921207466972532393966088D-29/
  data bth0cs( 39) / -0.24992261166978652421207213682763D-29/
  data bth0cs( 40) / +0.98399794299521955672828260355318D-30/
  data bth0cs( 41) / -0.39220375242408016397989131626158D-30/
  data bth0cs( 42) / +0.15818107030056522138590618845692D-30/
  data bth0cs( 43) / -0.64525506144890715944344098365426D-31/
  data bth0cs( 44) / +0.26611111369199356137177018346367D-31/

  data bm02cs(  1) / +0.9500415145228381369330861335560D-01/
  data bm02cs(  2) / -0.3801864682365670991748081566851D-03/
  data bm02cs(  3) / +0.2258339301031481192951829927224D-05/
  data bm02cs(  4) / -0.3895725802372228764730621412605D-07/
  data bm02cs(  5) / +0.1246886416512081697930990529725D-08/
  data bm02cs(  6) / -0.6065949022102503779803835058387D-10/
  data bm02cs(  7) / +0.4008461651421746991015275971045D-11/
  data bm02cs(  8) / -0.3350998183398094218467298794574D-12/
  data bm02cs(  9) / +0.3377119716517417367063264341996D-13/
  data bm02cs( 10) / -0.3964585901635012700569356295823D-14/
  data bm02cs( 11) / +0.5286111503883857217387939744735D-15/
  data bm02cs( 12) / -0.7852519083450852313654640243493D-16/
  data bm02cs( 13) / +0.1280300573386682201011634073449D-16/
  data bm02cs( 14) / -0.2263996296391429776287099244884D-17/
  data bm02cs( 15) / +0.4300496929656790388646410290477D-18/
  data bm02cs( 16) / -0.8705749805132587079747535451455D-19/
  data bm02cs( 17) / +0.1865862713962095141181442772050D-19/
  data bm02cs( 18) / -0.4210482486093065457345086972301D-20/
  data bm02cs( 19) / +0.9956676964228400991581627417842D-21/
  data bm02cs( 20) / -0.2457357442805313359605921478547D-21/
  data bm02cs( 21) / +0.6307692160762031568087353707059D-22/
  data bm02cs( 22) / -0.1678773691440740142693331172388D-22/
  data bm02cs( 23) / +0.4620259064673904433770878136087D-23/
  data bm02cs( 24) / -0.1311782266860308732237693402496D-23/
  data bm02cs( 25) / +0.3834087564116302827747922440276D-24/
  data bm02cs( 26) / -0.1151459324077741271072613293576D-24/
  data bm02cs( 27) / +0.3547210007523338523076971345213D-25/
  data bm02cs( 28) / -0.1119218385815004646264355942176D-25/
  data bm02cs( 29) / +0.3611879427629837831698404994257D-26/
  data bm02cs( 30) / -0.1190687765913333150092641762463D-26/
  data bm02cs( 31) / +0.4005094059403968131802476449536D-27/
  data bm02cs( 32) / -0.1373169422452212390595193916017D-27/
  data bm02cs( 33) / +0.4794199088742531585996491526437D-28/
  data bm02cs( 34) / -0.1702965627624109584006994476452D-28/
  data bm02cs( 35) / +0.6149512428936330071503575161324D-29/
  data bm02cs( 36) / -0.2255766896581828349944300237242D-29/
  data bm02cs( 37) / +0.8399707509294299486061658353200D-30/
  data bm02cs( 38) / -0.3172997595562602355567423936152D-30/
  data bm02cs( 39) / +0.1215205298881298554583333026514D-30/
  data bm02cs( 40) / -0.4715852749754438693013210568045D-31/

  data bt02cs(  1) / -0.24548295213424597462050467249324D+00/
  data bt02cs(  2) / +0.12544121039084615780785331778299D-02/
  data bt02cs(  3) / -0.31253950414871522854973446709571D-04/
  data bt02cs(  4) / +0.14709778249940831164453426969314D-05/
  data bt02cs(  5) / -0.99543488937950033643468850351158D-07/
  data bt02cs(  6) / +0.85493166733203041247578711397751D-08/
  data bt02cs(  7) / -0.86989759526554334557985512179192D-09/
  data bt02cs(  8) / +0.10052099533559791084540101082153D-09/
  data bt02cs(  9) / -0.12828230601708892903483623685544D-10/
  data bt02cs( 10) / +0.17731700781805131705655750451023D-11/
  data bt02cs( 11) / -0.26174574569485577488636284180925D-12/
  data bt02cs( 12) / +0.40828351389972059621966481221103D-13/
  data bt02cs( 13) / -0.66751668239742720054606749554261D-14/
  data bt02cs( 14) / +0.11365761393071629448392469549951D-14/
  data bt02cs( 15) / -0.20051189620647160250559266412117D-15/
  data bt02cs( 16) / +0.36497978794766269635720591464106D-16/
  data bt02cs( 17) / -0.68309637564582303169355843788800D-17/
  data bt02cs( 18) / +0.13107583145670756620057104267946D-17/
  data bt02cs( 19) / -0.25723363101850607778757130649599D-18/
  data bt02cs( 20) / +0.51521657441863959925267780949333D-19/
  data bt02cs( 21) / -0.10513017563758802637940741461333D-19/
  data bt02cs( 22) / +0.21820381991194813847301084501333D-20/
  data bt02cs( 23) / -0.46004701210362160577225905493333D-21/
  data bt02cs( 24) / +0.98407006925466818520953651199999D-22/
  data bt02cs( 25) / -0.21334038035728375844735986346666D-22/
  data bt02cs( 26) / +0.46831036423973365296066286933333D-23/
  data bt02cs( 27) / -0.10400213691985747236513382399999D-23/
  data bt02cs( 28) / +0.23349105677301510051777740800000D-24/
  data bt02cs( 29) / -0.52956825323318615788049749333333D-25/
  data bt02cs( 30) / +0.12126341952959756829196287999999D-25/
  data bt02cs( 31) / -0.28018897082289428760275626666666D-26/
  data bt02cs( 32) / +0.65292678987012873342593706666666D-27/
  data bt02cs( 33) / -0.15337980061873346427835733333333D-27/
  data bt02cs( 34) / +0.36305884306364536682359466666666D-28/
  data bt02cs( 35) / -0.86560755713629122479172266666666D-29/
  data bt02cs( 36) / +0.20779909972536284571238399999999D-29/
  data bt02cs( 37) / -0.50211170221417221674325333333333D-30/
  data bt02cs( 38) / +0.12208360279441714184191999999999D-30/
  data bt02cs( 39) / -0.29860056267039913454250666666666D-31/

  data nbm0 / 0 /
  data nbm02 / 0 /
  data nbt02 / 0 /
  data nbth0 / 0 /
  data pi4 / 0.785398163397448309615660845819876D+00 /
  data xmax / 0.0D+00 /

  if ( nbm0 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    nbm0 = r8_inits ( bm0cs, 37, eta )
    nbt02 = r8_inits ( bt02cs, 39, eta )
    nbm02 = r8_inits ( bm02cs, 40, eta )
    nbth0 = r8_inits ( bth0cs, 44, eta )
    xmax = 1.0D+00 / r8_mach ( 4 )
  end if

  if ( x < 4.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_B0MP - Fatal error!'
    write ( *, '(a)' ) '  X < 4.'
    stop
  else if ( x <= 8.0D+00 ) then
    z = ( 128.0D+00 / x / x - 5.0D+00 ) / 3.0D+00
    ampl = ( 0.75D+00 + r8_csevl ( z, bm0cs, nbm0 ) ) / sqrt ( x )
    theta = x - pi4 + r8_csevl ( z, bt02cs, nbt02 ) / x
  else
    z = 128.0D+00 / x / x - 1.0D+00
    ampl = ( 0.75D+00 + r8_csevl ( z, bm02cs, nbm02) ) / sqrt ( x )
    theta = x - pi4 + r8_csevl ( z, bth0cs, nbth0 ) / x
  end if

  return
end
function r8_besi1 ( x )

!*****************************************************************************80
!
!! R8_BESI1 evaluates the Bessel function I of order 1 of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESI1, the Bessel function I of order 1 of X.
!
  implicit none

  real ( kind = 8 ) bi1cs(17)
  integer ( kind = 4 ) nti1
  real ( kind = 8 ) r8_besi1
  real ( kind = 8 ) r8_besi1e
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save bi1cs
  save nti1
  save xmax
  save xmin
  save xsml

  data bi1cs(  1) / -0.19717132610998597316138503218149D-02 /
  data bi1cs(  2) / +0.40734887667546480608155393652014D+00 /
  data bi1cs(  3) / +0.34838994299959455866245037783787D-01 /
  data bi1cs(  4) / +0.15453945563001236038598401058489D-02 /
  data bi1cs(  5) / +0.41888521098377784129458832004120D-04 /
  data bi1cs(  6) / +0.76490267648362114741959703966069D-06 /
  data bi1cs(  7) / +0.10042493924741178689179808037238D-07 /
  data bi1cs(  8) / +0.99322077919238106481371298054863D-10 /
  data bi1cs(  9) / +0.76638017918447637275200171681349D-12 /
  data bi1cs( 10) / +0.47414189238167394980388091948160D-14 /
  data bi1cs( 11) / +0.24041144040745181799863172032000D-16 /
  data bi1cs( 12) / +0.10171505007093713649121100799999D-18 /
  data bi1cs( 13) / +0.36450935657866949458491733333333D-21 /
  data bi1cs( 14) / +0.11205749502562039344810666666666D-23 /
  data bi1cs( 15) / +0.29875441934468088832000000000000D-26 /
  data bi1cs( 16) / +0.69732310939194709333333333333333D-29 /
  data bi1cs( 17) / +0.14367948220620800000000000000000D-31 /

  data nti1 / 0 /
  data xmax / 0.0D+00 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( nti1 == 0 ) then
    nti1 = r8_inits ( bi1cs, 17, 0.1D+00 * r8_mach ( 3 ) )
    xmin = 2.0D+00 * r8_mach ( 1 )
    xsml = sqrt ( 8.0D+00 * r8_mach ( 3 ) )
    xmax = log ( r8_mach ( 2 ) )
  end if

  y = abs ( x )

  if ( y <= xmin ) then
    r8_besi1 = 0.0D+00
  else if ( y <= xsml ) then
    r8_besi1 = 0.5D+00 * x
  else if ( y <= 3.0D+00 ) then
    r8_besi1 = x * ( 0.875D+00 &
      + r8_csevl ( y * y / 4.5D+00 - 1.0D+00, bi1cs, nti1 ) )
  else if ( y <= xmax ) then
    r8_besi1 = exp ( y ) * r8_besi1e ( x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESI1 - Fatal error!'
    write ( *, '(a)' ) '  Result overflows.'
    stop
  end if

  return
end
function r8_besi1e ( x )

!*****************************************************************************80
!
!! R8_BESI1E evaluates the exponentially scaled Bessel function I1(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESI1E, the exponentially scaled Bessel 
!    function I1(X).
!
  implicit none

  real ( kind = 8 ) ai12cs(69)
  real ( kind = 8 ) ai1cs(46)
  real ( kind = 8 ) bi1cs(17)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) ntai1
  integer ( kind = 4 ) ntai12
  integer ( kind = 4 ) nti1
  real ( kind = 8 ) r8_besi1e
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save ai12cs
  save ai1cs
  save bi1cs
  save ntai1
  save ntai12
  save nti1
  save xmin
  save xsml

  data bi1cs(  1) / -0.19717132610998597316138503218149D-02 /
  data bi1cs(  2) / +0.40734887667546480608155393652014D+00 /
  data bi1cs(  3) / +0.34838994299959455866245037783787D-01 /
  data bi1cs(  4) / +0.15453945563001236038598401058489D-02 /
  data bi1cs(  5) / +0.41888521098377784129458832004120D-04 /
  data bi1cs(  6) / +0.76490267648362114741959703966069D-06 /
  data bi1cs(  7) / +0.10042493924741178689179808037238D-07 /
  data bi1cs(  8) / +0.99322077919238106481371298054863D-10 /
  data bi1cs(  9) / +0.76638017918447637275200171681349D-12 /
  data bi1cs( 10) / +0.47414189238167394980388091948160D-14 /
  data bi1cs( 11) / +0.24041144040745181799863172032000D-16 /
  data bi1cs( 12) / +0.10171505007093713649121100799999D-18 /
  data bi1cs( 13) / +0.36450935657866949458491733333333D-21 /
  data bi1cs( 14) / +0.11205749502562039344810666666666D-23 /
  data bi1cs( 15) / +0.29875441934468088832000000000000D-26 /
  data bi1cs( 16) / +0.69732310939194709333333333333333D-29 /
  data bi1cs( 17) / +0.14367948220620800000000000000000D-31 /

  data ai1cs(  1) / -0.2846744181881478674100372468307D-01 /
  data ai1cs(  2) / -0.1922953231443220651044448774979D-01 /
  data ai1cs(  3) / -0.6115185857943788982256249917785D-03 /
  data ai1cs(  4) / -0.2069971253350227708882823777979D-04 /
  data ai1cs(  5) / +0.8585619145810725565536944673138D-05 /
  data ai1cs(  6) / +0.1049498246711590862517453997860D-05 /
  data ai1cs(  7) / -0.2918338918447902202093432326697D-06 /
  data ai1cs(  8) / -0.1559378146631739000160680969077D-07 /
  data ai1cs(  9) / +0.1318012367144944705525302873909D-07 /
  data ai1cs( 10) / -0.1448423418183078317639134467815D-08 /
  data ai1cs( 11) / -0.2908512243993142094825040993010D-09 /
  data ai1cs( 12) / +0.1266388917875382387311159690403D-09 /
  data ai1cs( 13) / -0.1664947772919220670624178398580D-10 /
  data ai1cs( 14) / -0.1666653644609432976095937154999D-11 /
  data ai1cs( 15) / +0.1242602414290768265232168472017D-11 /
  data ai1cs( 16) / -0.2731549379672432397251461428633D-12 /
  data ai1cs( 17) / +0.2023947881645803780700262688981D-13 /
  data ai1cs( 18) / +0.7307950018116883636198698126123D-14 /
  data ai1cs( 19) / -0.3332905634404674943813778617133D-14 /
  data ai1cs( 20) / +0.7175346558512953743542254665670D-15 /
  data ai1cs( 21) / -0.6982530324796256355850629223656D-16 /
  data ai1cs( 22) / -0.1299944201562760760060446080587D-16 /
  data ai1cs( 23) / +0.8120942864242798892054678342860D-17 /
  data ai1cs( 24) / -0.2194016207410736898156266643783D-17 /
  data ai1cs( 25) / +0.3630516170029654848279860932334D-18 /
  data ai1cs( 26) / -0.1695139772439104166306866790399D-19 /
  data ai1cs( 27) / -0.1288184829897907807116882538222D-19 /
  data ai1cs( 28) / +0.5694428604967052780109991073109D-20 /
  data ai1cs( 29) / -0.1459597009090480056545509900287D-20 /
  data ai1cs( 30) / +0.2514546010675717314084691334485D-21 /
  data ai1cs( 31) / -0.1844758883139124818160400029013D-22 /
  data ai1cs( 32) / -0.6339760596227948641928609791999D-23 /
  data ai1cs( 33) / +0.3461441102031011111108146626560D-23 /
  data ai1cs( 34) / -0.1017062335371393547596541023573D-23 /
  data ai1cs( 35) / +0.2149877147090431445962500778666D-24 /
  data ai1cs( 36) / -0.3045252425238676401746206173866D-25 /
  data ai1cs( 37) / +0.5238082144721285982177634986666D-27 /
  data ai1cs( 38) / +0.1443583107089382446416789503999D-26 /
  data ai1cs( 39) / -0.6121302074890042733200670719999D-27 /
  data ai1cs( 40) / +0.1700011117467818418349189802666D-27 /
  data ai1cs( 41) / -0.3596589107984244158535215786666D-28 /
  data ai1cs( 42) / +0.5448178578948418576650513066666D-29 /
  data ai1cs( 43) / -0.2731831789689084989162564266666D-30 /
  data ai1cs( 44) / -0.1858905021708600715771903999999D-30 /
  data ai1cs( 45) / +0.9212682974513933441127765333333D-31 /
  data ai1cs( 46) / -0.2813835155653561106370833066666D-31 /

  data ai12cs(  1) / +0.2857623501828012047449845948469D-01  /
  data ai12cs(  2) / -0.9761097491361468407765164457302D-02  /
  data ai12cs(  3) / -0.1105889387626237162912569212775D-03  /
  data ai12cs(  4) / -0.3882564808877690393456544776274D-05  /
  data ai12cs(  5) / -0.2512236237870208925294520022121D-06  /
  data ai12cs(  6) / -0.2631468846889519506837052365232D-07  /
  data ai12cs(  7) / -0.3835380385964237022045006787968D-08  /
  data ai12cs(  8) / -0.5589743462196583806868112522229D-09  /
  data ai12cs(  9) / -0.1897495812350541234498925033238D-10 /
  data ai12cs( 10) / +0.3252603583015488238555080679949D-10 /
  data ai12cs( 11) / +0.1412580743661378133163366332846D-10 /
  data ai12cs( 12) / +0.2035628544147089507224526136840D-11 /
  data ai12cs( 13) / -0.7198551776245908512092589890446D-12 /
  data ai12cs( 14) / -0.4083551111092197318228499639691D-12 /
  data ai12cs( 15) / -0.2101541842772664313019845727462D-13 /
  data ai12cs( 16) / +0.4272440016711951354297788336997D-13 /
  data ai12cs( 17) / +0.1042027698412880276417414499948D-13 /
  data ai12cs( 18) / -0.3814403072437007804767072535396D-14 /
  data ai12cs( 19) / -0.1880354775510782448512734533963D-14 /
  data ai12cs( 20) / +0.3308202310920928282731903352405D-15 /
  data ai12cs( 21) / +0.2962628997645950139068546542052D-15 /
  data ai12cs( 22) / -0.3209525921993423958778373532887D-16 /
  data ai12cs( 23) / -0.4650305368489358325571282818979D-16 /
  data ai12cs( 24) / +0.4414348323071707949946113759641D-17 /
  data ai12cs( 25) / +0.7517296310842104805425458080295D-17 /
  data ai12cs( 26) / -0.9314178867326883375684847845157D-18 /
  data ai12cs( 27) / -0.1242193275194890956116784488697D-17 /
  data ai12cs( 28) / +0.2414276719454848469005153902176D-18 /
  data ai12cs( 29) / +0.2026944384053285178971922860692D-18 /
  data ai12cs( 30) / -0.6394267188269097787043919886811D-19 /
  data ai12cs( 31) / -0.3049812452373095896084884503571D-19 /
  data ai12cs( 32) / +0.1612841851651480225134622307691D-19 /
  data ai12cs( 33) / +0.3560913964309925054510270904620D-20 /
  data ai12cs( 34) / -0.3752017947936439079666828003246D-20 /
  data ai12cs( 35) / -0.5787037427074799345951982310741D-22 /
  data ai12cs( 36) / +0.7759997511648161961982369632092D-21 /
  data ai12cs( 37) / -0.1452790897202233394064459874085D-21 /
  data ai12cs( 38) / -0.1318225286739036702121922753374D-21 /
  data ai12cs( 39) / +0.6116654862903070701879991331717D-22 /
  data ai12cs( 40) / +0.1376279762427126427730243383634D-22 /
  data ai12cs( 41) / -0.1690837689959347884919839382306D-22 /
  data ai12cs( 42) / +0.1430596088595433153987201085385D-23 /
  data ai12cs( 43) / +0.3409557828090594020405367729902D-23 /
  data ai12cs( 44) / -0.1309457666270760227845738726424D-23 /
  data ai12cs( 45) / -0.3940706411240257436093521417557D-24 /
  data ai12cs( 46) / +0.4277137426980876580806166797352D-24 /
  data ai12cs( 47) / -0.4424634830982606881900283123029D-25 /
  data ai12cs( 48) / -0.8734113196230714972115309788747D-25 /
  data ai12cs( 49) / +0.4045401335683533392143404142428D-25 /
  data ai12cs( 50) / +0.7067100658094689465651607717806D-26 /
  data ai12cs( 51) / -0.1249463344565105223002864518605D-25 /
  data ai12cs( 52) / +0.2867392244403437032979483391426D-26 /
  data ai12cs( 53) / +0.2044292892504292670281779574210D-26 /
  data ai12cs( 54) / -0.1518636633820462568371346802911D-26 /
  data ai12cs( 55) / +0.8110181098187575886132279107037D-28 /
  data ai12cs( 56) / +0.3580379354773586091127173703270D-27 /
  data ai12cs( 57) / -0.1692929018927902509593057175448D-27 /
  data ai12cs( 58) / -0.2222902499702427639067758527774D-28 /
  data ai12cs( 59) / +0.5424535127145969655048600401128D-28 /
  data ai12cs( 60) / -0.1787068401578018688764912993304D-28 /
  data ai12cs( 61) / -0.6565479068722814938823929437880D-29 /
  data ai12cs( 62) / +0.7807013165061145280922067706839D-29 /
  data ai12cs( 63) / -0.1816595260668979717379333152221D-29 /
  data ai12cs( 64) / -0.1287704952660084820376875598959D-29 /
  data ai12cs( 65) / +0.1114548172988164547413709273694D-29 /
  data ai12cs( 66) / -0.1808343145039336939159368876687D-30 /
  data ai12cs( 67) / -0.2231677718203771952232448228939D-30 /
  data ai12cs( 68) / +0.1619029596080341510617909803614D-30 /
  data ai12cs( 69) / -0.1834079908804941413901308439210D-31 /

  data ntai1 / 0 /
  data ntai12 / 0 /
  data nti1 / 0 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( nti1 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    nti1 = r8_inits ( bi1cs, 17, eta )
    ntai1 = r8_inits ( ai1cs, 46, eta )
    ntai12 = r8_inits ( ai12cs, 69, eta )
    xmin = 2.0D+00 * r8_mach ( 1 )
    xsml = sqrt ( 8.0D+00 * r8_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= xmin ) then
    r8_besi1e = 0.0D+00
  else if ( y <= xsml ) then
    r8_besi1e = 0.5D+00 * x * exp ( - y )
  else if ( y <= 3.0D+00 ) then
    r8_besi1e = x * ( 0.875D+00 &
      + r8_csevl ( y * y / 4.5D+00 - 1.0D+00, bi1cs, nti1 ) ) &
      * exp ( - y ) 
  else if ( y <= 8.0D+00 ) then
    r8_besi1e = ( 0.375D+00 &
      + r8_csevl ( ( 48.0D+00 / y - 11.0D+00) / 5.0D+00, &
      ai1cs, ntai1 ) ) / sqrt ( y )
    if ( x < 0.0D+00 ) then
      r8_besi1e = - r8_besi1e
    end if
  else
    r8_besi1e = ( 0.375D+00 &
      + r8_csevl ( 16.0D+00 / y - 1.0D+00, ai12cs, ntai12 ) ) &
      / sqrt ( y )
    if ( x < 0.0D+00 ) then
      r8_besi1e = - r8_besi1e
    end if
  end if

  return
end
function r8_besj0 ( x )

!*****************************************************************************80
!
!! R8_BESJ0 evaluates the Bessel function J of order 0 of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESJ0, the Bessel function J of order 0 of X.
!
  implicit none

  real ( kind = 8 ) ampl
  real ( kind = 8 ) arg
  real ( kind = 8 ) bj0cs(19)
  integer ( kind = 4 ) ntj0
  real ( kind = 8 ) r8_besj0
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) theta
  real ( kind = 8 ) x
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save bj0cs
  save ntj0
  save xsml

  data bj0cs(  1) / +0.10025416196893913701073127264074D+00 /
  data bj0cs(  2) / -0.66522300776440513177678757831124D+00 /
  data bj0cs(  3) / +0.24898370349828131370460468726680D+00 /
  data bj0cs(  4) / -0.33252723170035769653884341503854D-01 /
  data bj0cs(  5) / +0.23114179304694015462904924117729D-02 /
  data bj0cs(  6) / -0.99112774199508092339048519336549D-04 /
  data bj0cs(  7) / +0.28916708643998808884733903747078D-05 /
  data bj0cs(  8) / -0.61210858663032635057818407481516D-07 /
  data bj0cs(  9) / +0.98386507938567841324768748636415D-09 /
  data bj0cs( 10) / -0.12423551597301765145515897006836D-10 /
  data bj0cs( 11) / +0.12654336302559045797915827210363D-12 /
  data bj0cs( 12) / -0.10619456495287244546914817512959D-14 /
  data bj0cs( 13) / +0.74706210758024567437098915584000D-17 /
  data bj0cs( 14) / -0.44697032274412780547627007999999D-19 /
  data bj0cs( 15) / +0.23024281584337436200523093333333D-21 /
  data bj0cs( 16) / -0.10319144794166698148522666666666D-23 /
  data bj0cs( 17) / +0.40608178274873322700800000000000D-26 /
  data bj0cs( 18) / -0.14143836005240913919999999999999D-28 /
  data bj0cs( 19) / +0.43910905496698880000000000000000D-31 /

  data ntj0 / 0 /
  data xsml / 0.0D+00 /

  if ( ntj0 == 0 ) then
    ntj0 = r8_inits ( bj0cs, 19, 0.1D+00 * r8_mach ( 3 ) )
    xsml = sqrt ( 4.0D+00 * r8_mach ( 3 ) )
  end if

  y = abs ( x )

  if ( y <= xsml ) then
    r8_besj0 = 1.0D+00
  else if ( y <= 4.0D+00 ) then
    arg = 0.125D+00 * y * y - 1.0D+00
    r8_besj0 = r8_csevl ( arg, bj0cs, ntj0 )
  else
    call r8_b0mp ( y, ampl, theta )
    r8_besj0 = ampl * cos ( theta )
  end if

  return
end
function r8_besk ( nu, x )

!*****************************************************************************80
!
!! R8_BESK evaluates the Bessel function K of order NU of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 November 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) NU, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESK, the Bessel function K of order NU at X.
!
  implicit none

  real ( kind = 8 ), allocatable :: bke(:)
  integer ( kind = 4 ) nin
  real ( kind = 8 ) nu
  real ( kind = 8 ) r8_besk
  real ( kind = 8 ) x
  real ( kind = 8 ) xnu

  xnu = nu - int ( nu )
  nin = int ( nu ) + 1
  allocate ( bke(1:nin) )

  call r8_besks ( xnu, x, nin, bke )

  r8_besk = bke(nin)

  deallocate ( bke )

  return
end
function r8_besk1 ( x )

!*****************************************************************************80
!
!! R8_BESK1 evaluates the Bessel function K of order 1 of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESK1, the Bessel function K of order 1 of X.
!
  implicit none

  real ( kind = 8 ) bk1cs(16)
  integer ( kind = 4 ) ntk1
  real ( kind = 8 ) r8_besi1
  real ( kind = 8 ) r8_besk1
  real ( kind = 8 ) r8_besk1e
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save bk1cs
  save ntk1
  save xmax
  save xmin
  save xsml

  data bk1cs(  1) / +0.25300227338947770532531120868533D-01 /
  data bk1cs(  2) / -0.35315596077654487566723831691801D+00 /
  data bk1cs(  3) / -0.12261118082265714823479067930042D+00 /
  data bk1cs(  4) / -0.69757238596398643501812920296083D-02 /
  data bk1cs(  5) / -0.17302889575130520630176507368979D-03 /
  data bk1cs(  6) / -0.24334061415659682349600735030164D-05 /
  data bk1cs(  7) / -0.22133876307347258558315252545126D-07 /
  data bk1cs(  8) / -0.14114883926335277610958330212608D-09 /
  data bk1cs(  9) / -0.66669016941993290060853751264373D-12 /
  data bk1cs( 10) / -0.24274498505193659339263196864853D-14 /
  data bk1cs( 11) / -0.70238634793862875971783797120000D-17 /
  data bk1cs( 12) / -0.16543275155100994675491029333333D-19 /
  data bk1cs( 13) / -0.32338347459944491991893333333333D-22 /
  data bk1cs( 14) / -0.53312750529265274999466666666666D-25 /
  data bk1cs( 15) / -0.75130407162157226666666666666666D-28 /
  data bk1cs( 16) / -0.91550857176541866666666666666666D-31 /

  data ntk1 / 0 /
  data xmax / 0.0D+00 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( ntk1 == 0 ) then
    ntk1 = r8_inits ( bk1cs, 16, 0.1D+00 * r8_mach ( 3 ) )
    xmin = exp ( max ( log ( r8_mach ( 1 ) ), &
      - log ( r8_mach ( 2 ) ) ) + 0.01D+00 )
    xsml = sqrt ( 4.0D+00 * r8_mach ( 3 ) )
    xmax = - log ( r8_mach ( 1 ) )
    xmax = xmax - 0.5D+00 * xmax * log ( xmax ) &
      / ( xmax + 0.5D+00 ) - 0.01D+00
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESK1 = Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0D+00
    r8_besk1 = log ( 0.5D+00 * x ) * r8_besi1 ( x ) + ( 0.75D+00 &
      + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk1cs, ntk1 ) ) / x
  else if ( x <= 2.0D+00 ) then
    y = x * x
    r8_besk1 = log ( 0.5D+00 * x ) * r8_besi1 ( x ) + ( 0.75D+00 &
      + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk1cs, ntk1 ) ) / x
  else if ( x <= xmax ) then
    r8_besk1 = exp ( - x ) * r8_besk1e ( x )
  else
    r8_besk1 = 0.0D+00
  end if

  return
end
function r8_besk1e ( x )

!*****************************************************************************80
!
!! R8_BESK1E evaluates the exponentially scaled Bessel function K1(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_BESK1E, the exponentially scaled Bessel 
!    function K1(X).
!
  implicit none

  real ( kind = 8 ) ak12cs(33)
  real ( kind = 8 ) ak1cs(38)
  real ( kind = 8 ) bk1cs(16)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) ntak1
  integer ( kind = 4 ) ntak12
  integer ( kind = 4 ) ntk1
  real ( kind = 8 ) r8_besi1
  real ( kind = 8 ) r8_besk1e
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save ak12cs
  save ak1cs
  save bk1cs
  save ntak1
  save ntak12
  save ntk1
  save xmin
  save xsml

  data bk1cs(  1) / +0.25300227338947770532531120868533D-01 /
  data bk1cs(  2) / -0.35315596077654487566723831691801D+00 /
  data bk1cs(  3) / -0.12261118082265714823479067930042D+00 /
  data bk1cs(  4) / -0.69757238596398643501812920296083D-02 /
  data bk1cs(  5) / -0.17302889575130520630176507368979D-03 /
  data bk1cs(  6) / -0.24334061415659682349600735030164D-05 /
  data bk1cs(  7) / -0.22133876307347258558315252545126D-07 /
  data bk1cs(  8) / -0.14114883926335277610958330212608D-09 /
  data bk1cs(  9) / -0.66669016941993290060853751264373D-12 /
  data bk1cs( 10) / -0.24274498505193659339263196864853D-14 /
  data bk1cs( 11) / -0.70238634793862875971783797120000D-17 /
  data bk1cs( 12) / -0.16543275155100994675491029333333D-19 /
  data bk1cs( 13) / -0.32338347459944491991893333333333D-22 /
  data bk1cs( 14) / -0.53312750529265274999466666666666D-25 /
  data bk1cs( 15) / -0.75130407162157226666666666666666D-28 /
  data bk1cs( 16) / -0.91550857176541866666666666666666D-31 /

  data ak1cs(  1) / +0.27443134069738829695257666227266D+00 /
  data ak1cs(  2) / +0.75719899531993678170892378149290D-01 /
  data ak1cs(  3) / -0.14410515564754061229853116175625D-02 /
  data ak1cs(  4) / +0.66501169551257479394251385477036D-04 /
  data ak1cs(  5) / -0.43699847095201407660580845089167D-05 /
  data ak1cs(  6) / +0.35402774997630526799417139008534D-06 /
  data ak1cs(  7) / -0.33111637792932920208982688245704D-07 /
  data ak1cs(  8) / +0.34459775819010534532311499770992D-08 /
  data ak1cs(  9) / -0.38989323474754271048981937492758D-09 /
  data ak1cs( 10) / +0.47208197504658356400947449339005D-10 /
  data ak1cs( 11) / -0.60478356628753562345373591562890D-11 /
  data ak1cs( 12) / +0.81284948748658747888193837985663D-12 /
  data ak1cs( 13) / -0.11386945747147891428923915951042D-12 /
  data ak1cs( 14) / +0.16540358408462282325972948205090D-13 /
  data ak1cs( 15) / -0.24809025677068848221516010440533D-14 /
  data ak1cs( 16) / +0.38292378907024096948429227299157D-15 /
  data ak1cs( 17) / -0.60647341040012418187768210377386D-16 /
  data ak1cs( 18) / +0.98324256232648616038194004650666D-17 /
  data ak1cs( 19) / -0.16284168738284380035666620115626D-17 /
  data ak1cs( 20) / +0.27501536496752623718284120337066D-18 /
  data ak1cs( 21) / -0.47289666463953250924281069568000D-19 /
  data ak1cs( 22) / +0.82681500028109932722392050346666D-20 /
  data ak1cs( 23) / -0.14681405136624956337193964885333D-20 /
  data ak1cs( 24) / +0.26447639269208245978085894826666D-21 /
  data ak1cs( 25) / -0.48290157564856387897969868800000D-22 /
  data ak1cs( 26) / +0.89293020743610130180656332799999D-23 /
  data ak1cs( 27) / -0.16708397168972517176997751466666D-23 /
  data ak1cs( 28) / +0.31616456034040694931368618666666D-24 /
  data ak1cs( 29) / -0.60462055312274989106506410666666D-25 /
  data ak1cs( 30) / +0.11678798942042732700718421333333D-25 /
  data ak1cs( 31) / -0.22773741582653996232867840000000D-26 /
  data ak1cs( 32) / +0.44811097300773675795305813333333D-27 /
  data ak1cs( 33) / -0.88932884769020194062336000000000D-28 /
  data ak1cs( 34) / +0.17794680018850275131392000000000D-28 /
  data ak1cs( 35) / -0.35884555967329095821994666666666D-29 /
  data ak1cs( 36) / +0.72906290492694257991679999999999D-30 /
  data ak1cs( 37) / -0.14918449845546227073024000000000D-30 /
  data ak1cs( 38) / +0.30736573872934276300799999999999D-31 /

  data ak12cs(  1) / +0.6379308343739001036600488534102D-01 /
  data ak12cs(  2) / +0.2832887813049720935835030284708D-01 /
  data ak12cs(  3) / -0.2475370673905250345414545566732D-03 /
  data ak12cs(  4) / +0.5771972451607248820470976625763D-05 /
  data ak12cs(  5) / -0.2068939219536548302745533196552D-06 /
  data ak12cs(  6) / +0.9739983441381804180309213097887D-08 /
  data ak12cs(  7) / -0.5585336140380624984688895511129D-09 /
  data ak12cs(  8) / +0.3732996634046185240221212854731D-10 /
  data ak12cs(  9) / -0.2825051961023225445135065754928D-11 /
  data ak12cs( 10) / +0.2372019002484144173643496955486D-12 /
  data ak12cs( 11) / -0.2176677387991753979268301667938D-13 /
  data ak12cs( 12) / +0.2157914161616032453939562689706D-14 /
  data ak12cs( 13) / -0.2290196930718269275991551338154D-15 /
  data ak12cs( 14) / +0.2582885729823274961919939565226D-16 /
  data ak12cs( 15) / -0.3076752641268463187621098173440D-17 /
  data ak12cs( 16) / +0.3851487721280491597094896844799D-18 /
  data ak12cs( 17) / -0.5044794897641528977117282508800D-19 /
  data ak12cs( 18) / +0.6888673850418544237018292223999D-20 /
  data ak12cs( 19) / -0.9775041541950118303002132480000D-21 /
  data ak12cs( 20) / +0.1437416218523836461001659733333D-21 /
  data ak12cs( 21) / -0.2185059497344347373499733333333D-22 /
  data ak12cs( 22) / +0.3426245621809220631645388800000D-23 /
  data ak12cs( 23) / -0.5531064394246408232501248000000D-24 /
  data ak12cs( 24) / +0.9176601505685995403782826666666D-25 /
  data ak12cs( 25) / -0.1562287203618024911448746666666D-25 /
  data ak12cs( 26) / +0.2725419375484333132349439999999D-26 /
  data ak12cs( 27) / -0.4865674910074827992378026666666D-27 /
  data ak12cs( 28) / +0.8879388552723502587357866666666D-28 /
  data ak12cs( 29) / -0.1654585918039257548936533333333D-28 /
  data ak12cs( 30) / +0.3145111321357848674303999999999D-29 /
  data ak12cs( 31) / -0.6092998312193127612416000000000D-30 /
  data ak12cs( 32) / +0.1202021939369815834623999999999D-30 /
  data ak12cs( 33) / -0.2412930801459408841386666666666D-31 /

  data ntak1 / 0 /
  data ntak12 / 0 /
  data ntk1 / 0 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( ntk1 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    ntk1 = r8_inits ( bk1cs, 16, eta )
    ntak1 = r8_inits ( ak1cs, 38, eta )
    ntak12 = r8_inits ( ak12cs, 33, eta )
    xmin = exp ( max ( log ( r8_mach ( 1 ) ), &
      - log ( r8_mach ( 2 ) ) ) + 0.01D+00 )
    xsml = sqrt ( 4.0D+00 * r8_mach ( 3 ) )
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESK1E = Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= xsml ) then
    y = 0.0D+00
    r8_besk1e = exp ( x ) * ( log ( 0.5D+00 * x ) * r8_besi1 ( x ) &
      + ( 0.75D+00 &
      + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk1cs, ntk1 ) ) / x )
  else if ( x <= 2.0D+00 ) then
    y = x * x
    r8_besk1e = exp ( x ) * ( log ( 0.5D+00 * x ) * r8_besi1 ( x ) &
      + ( 0.75D+00 &
      + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk1cs, ntk1 ) ) / x )
  else if ( x <= 8.0D+00 ) then
    r8_besk1e = ( 1.25D+00 &
      + r8_csevl ( ( 16.0D+00 / x - 5.0D+00 ) / 3.0D+00, ak1cs, &
      ntak1 ) ) / sqrt ( x )
  else
    r8_besk1e = ( 1.25D+00 + &
      r8_csevl ( 16.0D+00 / x - 1.0D+00, ak12cs, ntak12 ) ) &
      / sqrt ( x )
  end if

  return
end
subroutine r8_beskes ( xnu, x, nin, bke )

!*****************************************************************************80
!
!! R8_BESKES: a sequence of exponentially scaled K Bessel functions at X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XNU, ?
!    |XNU| < 1.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, integer ( kind = 4 ) NIN, indicates the number of terms to compute.
!
!    Output, real ( kind = 8 ) BKE(abs(NIN)), the exponentially scaled 
!    K Bessel functions.
!
  implicit none

  real ( kind = 8 ) alnbig
  real ( kind = 8 ) bke(*)
  real ( kind = 8 ) bknu1
  real ( kind = 8 ) direct
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iswtch
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nin
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) v
  real ( kind = 8 ) vend
  real ( kind = 8 ) vincr
  real ( kind = 8 ) x
  real ( kind = 8 ) xnu

  save alnbig

  data alnbig / 0.0D+00 /

  if ( alnbig == 0.0D+00 ) then
    alnbig = log ( r8_mach ( 2 ) )
  end if

  v = abs ( xnu )
  n = abs ( nin )

  if ( 1.0D+00 <= v ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESKES - Fatal error!'
    write ( *, '(a)' ) '  |XNU| must be less than 1.'
    stop
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESKES - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  if ( n == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_BESKES - Fatal error!'
    write ( *, '(a)' ) '  N = 0.'
    stop
  end if

  call r8_knus ( v, x, bke(1), bknu1, iswtch )

  if ( n == 1 ) then
    return
  end if

  if ( nin < 0 ) then
    vincr = - 1.0D+00
  else
    vincr = + 1.0D+00
  end if

  if ( xnu < 0.0D+00 ) then
    direct = - vincr
  else
    direct = vincr
  end if

  bke(2) = bknu1

  if ( direct < 0.0D+00 ) then
    call r8_knus ( abs ( xnu + vincr ), x, bke(2), bknu1, iswtch )
  end if

  if ( n == 2 ) then
    return
  end if

  vend = abs ( xnu + real ( nin, kind = 8 ) ) - 1.0D+00

  v = xnu
  do i = 3, n
    v = v + vincr
    bke(i) = 2.0D+00 * v * bke(i-1) / x + bke(i-2)
  end do

  return
end
subroutine r8_besks ( xnu, x, nin, bk )

!*****************************************************************************80
!
!! R8_BESKS evaluates a sequence of K Bessel functions at X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XNU, ?
!    |XNU| < 1.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, integer ( kind = 4 ) NIN, indicates the number of terms to compute.
!
!    Output, real ( kind = 8 ) BK(abs(NIN)), the K Bessel functions.
!
  implicit none

  integer ( kind = 4 ) nin

  real ( kind = 8 ) bk(abs(nin))
  real ( kind = 8 ) expxi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xnu

  save xmax

  data xmax / 0.0D+00 /

  if ( xmax == 0.0D+00 ) then
    xmax = - log ( r8_mach ( 1 ) )
    xmax = xmax + 0.5D+00 * log ( 3.14D+00 * 0.5D+00 / xmax )
  end if

  call r8_beskes ( xnu, x, nin, bk )

  expxi = exp ( - x )
  n = abs ( nin )

  do i = 1, n
    bk(i) = expxi * bk(i)
  end do

  return
end
function r8_csevl ( x, a, n )

!*****************************************************************************80
!
!! R8_CSEVL evaluates a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, Number 4, April 1973, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Input, real ( kind = 8 ) CS(N), the Chebyshev coefficients.
!
!    Input, integer ( kind = 4 ) N, the number of Chebyshev coefficients.
!
!    Output, real ( kind = 8 ) R8_CSEVL, the Chebyshev series evaluated at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) twox
  real ( kind = 8 ) x

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
    write ( *, '(a)' ) '  Number of terms <= 0.'
    stop
  end if

  if ( 1000 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
    write ( *, '(a)' ) '  Number of terms > 1000.'
    stop
  end if

  if ( x < -1.1D+00 .or. 1.1D+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
    write ( *, '(a)' ) '  X outside (-1,+1)'
    write ( *, '(a,g14.6)' ) '  X = ', x
    stop
  end if

  twox = 2.0D+00 * x
  b1 = 0.0D+00
  b0 = 0.0D+00

  do i = n, 1, -1
    b2 = b1
    b1 = b0
    b0 = twox * b1 - b2 + a(i)
  end do

  r8_csevl = 0.5D+00 * ( b0 - b2 )

  return
end
subroutine r8_gaml ( xmin, xmax )

!*****************************************************************************80
!
!! R8_GAML evaluates bounds for an R8 argument of the gamma function.
!
!  Discussion:
!
!    This function calculates the minimum and maximum legal bounds 
!    for X in the evaluation of GAMMA ( X ).
!
!    XMIN and XMAX are not the only bounds, but they are the only 
!    non-trivial ones to calculate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) XMIN, XMAX, the bounds.
!
  implicit none

  real ( kind = 8 ) alnbig
  real ( kind = 8 ) alnsml
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) xln
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xold

  alnsml = log ( r8_mach ( 1 ) )
  xmin = - alnsml

  do i = 1, 10

    xold = xmin
    xln = log ( xmin )
    xmin = xmin - xmin * ( ( xmin + 0.5D+00 ) * xln - xmin &
      - 0.2258D+00 + alnsml ) / ( xmin * xln + 0.5D+00 )

    if ( abs ( xmin - xold ) < 0.005D+00 ) then

      xmin = - xmin + 0.01D+00

      alnbig = log ( r8_mach ( 2 ) )
      xmax = alnbig

      do j = 1, 10

        xold = xmax
        xln = log ( xmax )
        xmax = xmax - xmax * ( ( xmax - 0.5D+00 ) * xln - xmax &
          + 0.9189D+00 - alnbig ) / ( xmax * xln - 0.5D+00 )

        if ( abs ( xmax - xold ) < 0.005D+00 ) then
          xmax = xmax - 0.01D+00
          xmin = max ( xmin, - xmax + 1.0D+00 )
          return
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_GAML - Fatal error!'
      write ( *, '(a)' ) '  Unable to find XMAX.'
      stop

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_GAML - Fatal error!'
  write ( *, '(a)' ) '  Unable to find XMIN.'

  stop
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates the gamma function of an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the gamma function of X.
!
  implicit none

  real ( kind = 8 ) dxrel
  real ( kind = 8 ) gcs(42)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ngcs
  real ( kind = 8 ) pi
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) r8_gamma
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_lgmc
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) sinpiy
  real ( kind = 8 ) sq2pil
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xsml
  real ( kind = 8 ) y

  save dxrel
  save gcs
  save ngcs
  save pi
  save sq2pil
  save xmax
  save xmin
  save xsml

  data gcs(  1) / +0.8571195590989331421920062399942D-02 /
  data gcs(  2) / +0.4415381324841006757191315771652D-02 /
  data gcs(  3) / +0.5685043681599363378632664588789D-01 /
  data gcs(  4) / -0.4219835396418560501012500186624D-02 /
  data gcs(  5) / +0.1326808181212460220584006796352D-02 /
  data gcs(  6) / -0.1893024529798880432523947023886D-03 /
  data gcs(  7) / +0.3606925327441245256578082217225D-04 /
  data gcs(  8) / -0.6056761904460864218485548290365D-05 /
  data gcs(  9) / +0.1055829546302283344731823509093D-05 /
  data gcs( 10) / -0.1811967365542384048291855891166D-06 /
  data gcs( 11) / +0.3117724964715322277790254593169D-07 /
  data gcs( 12) / -0.5354219639019687140874081024347D-08 /
  data gcs( 13) / +0.9193275519859588946887786825940D-09 /
  data gcs( 14) / -0.1577941280288339761767423273953D-09 /
  data gcs( 15) / +0.2707980622934954543266540433089D-10 /
  data gcs( 16) / -0.4646818653825730144081661058933D-11 /
  data gcs( 17) / +0.7973350192007419656460767175359D-12 /
  data gcs( 18) / -0.1368078209830916025799499172309D-12 /
  data gcs( 19) / +0.2347319486563800657233471771688D-13 /
  data gcs( 20) / -0.4027432614949066932766570534699D-14 /
  data gcs( 21) / +0.6910051747372100912138336975257D-15 /
  data gcs( 22) / -0.1185584500221992907052387126192D-15 /
  data gcs( 23) / +0.2034148542496373955201026051932D-16 /
  data gcs( 24) / -0.3490054341717405849274012949108D-17 /
  data gcs( 25) / +0.5987993856485305567135051066026D-18 /
  data gcs( 26) / -0.1027378057872228074490069778431D-18 /
  data gcs( 27) / +0.1762702816060529824942759660748D-19 /
  data gcs( 28) / -0.3024320653735306260958772112042D-20 /
  data gcs( 29) / +0.5188914660218397839717833550506D-21 /
  data gcs( 30) / -0.8902770842456576692449251601066D-22 /
  data gcs( 31) / +0.1527474068493342602274596891306D-22 /
  data gcs( 32) / -0.2620731256187362900257328332799D-23 /
  data gcs( 33) / +0.4496464047830538670331046570666D-24 /
  data gcs( 34) / -0.7714712731336877911703901525333D-25 /
  data gcs( 35) / +0.1323635453126044036486572714666D-25 /
  data gcs( 36) / -0.2270999412942928816702313813333D-26 /
  data gcs( 37) / +0.3896418998003991449320816639999D-27 /
  data gcs( 38) / -0.6685198115125953327792127999999D-28 /
  data gcs( 39) / +0.1146998663140024384347613866666D-28 /
  data gcs( 40) / -0.1967938586345134677295103999999D-29 /
  data gcs( 41) / +0.3376448816585338090334890666666D-30 /
  data gcs( 42) / -0.5793070335782135784625493333333D-31 /

  data dxrel / 0.0D+00 /
  data ngcs / 0 /
  data pi / 3.14159265358979323846264338327950D+00 /
  data sq2pil / 0.91893853320467274178032973640562D+00 /
  data xmax / 0.0D+00 /
  data xmin / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( ngcs == 0 ) then
    ngcs = r8_inits ( gcs, 42, 0.1D+00 * r8_mach ( 3 ) )
    call r8_gaml ( xmin, xmax )
    xsml = exp ( max ( log ( r8_mach ( 1 ) ), &
      - log ( r8_mach ( 2 ) ) ) + 0.01D+00 )
    dxrel = sqrt ( r8_mach ( 4 ) )
  end if

  y = abs ( x )

  if ( y <= 10.0D+00 ) then

    n = int ( x )
    if ( x < 0.0D+00 ) then
      n = n - 1
    end if
    y = x - real ( n, kind = 8 )
    n = n - 1
    r8_gamma = 0.9375D+00 + r8_csevl ( 2.0D+00 * y - 1.0D+00, gcs, ngcs )

    if ( n == 0 ) then

      return

    else if ( n < 0 ) then

      n = - n

      if ( x == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
        write ( *, '(a)' ) '  X is 0.'
        stop
      end if

      if ( x < 0.0D+00 .and. x + real ( n - 2, kind = 8 ) == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
        write ( *, '(a)' ) '  X is a negative integer.'
        stop
      end if

      if ( x < - 0.5D+00 .and. &
        abs ( ( x - aint ( x - 0.5D+00 ) ) / x ) < dxrel ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_GAMMA - Warning!'
        write ( *, '(a)' ) '  X too near a negative integer,'
        write ( *, '(a)' ) '  answer is half precision.'
      end if

      if ( y < xsml ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
        write ( *, '(a)' ) '  X is so close to zero that Gamma overflows.'
        stop
      end if

      do i = 1, n
        r8_gamma = r8_gamma / ( x + real ( i - 1, kind = 8 ) )
      end do

    else if ( n == 0 ) then

    else

      do i = 1, n
        r8_gamma = ( y + real ( i, kind = 8 ) ) * r8_gamma
      end do

    end if

  else

    if ( xmax < x ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
      write ( *, '(a)' ) '  X so big that Gamma overflows.'
      stop
    end if
!
!  Underflow.
!
    if ( x < xmin ) then
      r8_gamma = 0.0D+00
      return
    end if

    r8_gamma = exp ( ( y - 0.5D+00 ) * log ( y ) - y + sq2pil + r8_lgmc ( y ) )

    if ( 0.0D+00 < x ) then
      return
    end if

    if ( abs ( ( x - aint ( x - 0.5D+00 ) ) / x ) < dxrel ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_GAMMA - Warning!'
      write ( *, '(a)' ) '  X too near a negative integer,'
      write ( *, '(a)' ) '  answer is half precision.'
    end if

    sinpiy = sin ( pi * y )

    if ( sinpiy == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
      write ( *, '(a)' ) '  X is a negative integer.'
      stop
    end if

    r8_gamma = - pi / ( y * sinpiy * r8_gamma )

  end if

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real ( kind = 8 ) r8_huge

  r8_huge = 1.0D+30

  return
end
function r8_inits ( dos, nos, eta )

!*****************************************************************************80
!
!! R8_INITS initializes a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, Number 4, April 1973, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DOS(NOS), the Chebyshev coefficients.
!
!    Input, integer ( kind = 4 ) NOS, the number of coefficients.
!
!    Input, real ( kind = 8 ) ETA, the desired accuracy.
!
!    Output, integer ( kind = 4 ) R8_INITS, the number of terms of the 
!    series needed to ensure the requested accuracy.
!
  implicit none

  integer ( kind = 4 ) nos

  real ( kind = 8 ) dos(nos)
  real ( kind = 8 ) err 
  real ( kind = 8 ) eta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) r8_inits

  if ( nos < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_INITS - Fatal error!'
    write ( *, '(a)' ) '  Number of coefficients < 1.'
    stop
  end if

  err = 0.0D+00

  do i = nos, 1, -1
    err = err + abs ( dos(i) )
    if ( eta < err ) then
      r8_inits = i
      return
    end if
  end do

  r8_inits = nos
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_INITS - Warning!'
  write ( *, '(a)' ) '  ETA may be too small.'

  return
end
subroutine r8_knus ( xnu, x, bknu, bknu1, iswtch )

!*****************************************************************************80
!
!! R8_KNUS computes a sequence of K Bessel functions.
!
!  Discussion:
!
!    This routine computes Bessel functions 
!      exp(x) * k-sub-xnu (x)  
!    and
!      exp(x) * k-sub-xnu+1 (x) 
!    for 0.0 <= xnu < 1.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XNU, the order parameter.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) BKNU, BKNU1, the two K Bessel functions.
!
!    Output, integer ( kind = 4 ) ISWTCH, ?
!
  implicit none

  real ( kind = 8 ) a(32)
  real ( kind = 8 ) a0
  real ( kind = 8 ) aln2
  real ( kind = 8 ) alnbig
  real ( kind = 8 ) alneps
  real ( kind = 8 ) alnsml
  real ( kind = 8 ) alnz
  real ( kind = 8 ) alpha(32)
  real ( kind = 8 ) an
  real ( kind = 8 ) b0
  real ( kind = 8 ) beta(32)
  real ( kind = 8 ) bknu
  real ( kind = 8 ) bknu0
  real ( kind = 8 ) bknu1
  real ( kind = 8 ) bknud
  real ( kind = 8 ) bn
  real ( kind = 8 ) c0
  real ( kind = 8 ) c0kcs(29)
  real ( kind = 8 ) eta
  real ( kind = 8 ) euler
  real ( kind = 8 ) expx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inu
  integer ( kind = 4 ) iswtch
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ntc0k
  integer ( kind = 4 ) nterms
  integer ( kind = 4 ) ntznu1
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) p3
  real ( kind = 8 ) qq
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) r8_gamma
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) result
  real ( kind = 8 ) sqpi2
  real ( kind = 8 ) sqrtx
  real ( kind = 8 ) v
  real ( kind = 8 ) vlnz
  real ( kind = 8 ) x
  real ( kind = 8 ) x2n
  real ( kind = 8 ) x2tov
  real ( kind = 8 ) xi
  real ( kind = 8 ) xmu
  real ( kind = 8 ) xnu
  real ( kind = 8 ) xnusml
  real ( kind = 8 ) xsml
  real ( kind = 8 ) z
  real ( kind = 8 ) znu1cs(20)
  real ( kind = 8 ) ztov

  save aln2
  save alnbig
  save alneps
  save alnsml
  save c0kcs
  save euler
  save ntc0k
  save ntznu1
  save sqpi2
  save xnusml
  save xsml
  save znu1cs

  data c0kcs(  1) / +0.60183057242626108387577445180329D-01     /
  data c0kcs(  2) / -0.15364871433017286092959755943124D+00     /
  data c0kcs(  3) / -0.11751176008210492040068229226213D-01     /
  data c0kcs(  4) / -0.85248788891979509827048401550987D-03     /
  data c0kcs(  5) / -0.61329838767496791874098176922111D-04     /
  data c0kcs(  6) / -0.44052281245510444562679889548505D-05     /
  data c0kcs(  7) / -0.31631246728384488192915445892199D-06     /
  data c0kcs(  8) / -0.22710719382899588330673771793396D-07     /
  data c0kcs(  9) / -0.16305644608077609552274620515360D-08     /
  data c0kcs( 10) / -0.11706939299414776568756044043130D-09     /
  data c0kcs( 11) / -0.84052063786464437174546593413792D-11    /
  data c0kcs( 12) / -0.60346670118979991487096050737198D-12    /
  data c0kcs( 13) / -0.43326960335681371952045997366903D-13    /
  data c0kcs( 14) / -0.31107358030203546214634697772237D-14    /
  data c0kcs( 15) / -0.22334078226736982254486133409840D-15    /
  data c0kcs( 16) / -0.16035146716864226300635791528610D-16    /
  data c0kcs( 17) / -0.11512717363666556196035697705305D-17    /
  data c0kcs( 18) / -0.82657591746836959105169479089258D-19    /
  data c0kcs( 19) / -0.59345480806383948172333436695984D-20    /
  data c0kcs( 20) / -0.42608138196467143926499613023976D-21    /
  data c0kcs( 21) / -0.30591266864812876299263698370542D-22    /
  data c0kcs( 22) / -0.21963541426734575224975501815516D-23    /
  data c0kcs( 23) / -0.15769113261495836071105750684760D-24    /
  data c0kcs( 24) / -0.11321713935950320948757731048056D-25    /
  data c0kcs( 25) / -0.81286248834598404082792349714433D-27    /
  data c0kcs( 26) / -0.58360900893453226552829349315949D-28    /
  data c0kcs( 27) / -0.41901241623610922519452337780905D-29    /
  data c0kcs( 28) / -0.30083737960206435069530504212862D-30    /
  data c0kcs( 29) / -0.21599152067808647728342168089832D-31    /

  data znu1cs(  1) / +0.203306756994191729674444001216911D+00    /
  data znu1cs(  2) / +0.140077933413219771062943670790563D+00    /
  data znu1cs(  3) / +0.791679696100161352840972241972320D-02    /
  data znu1cs(  4) / +0.339801182532104045352930092205750D-03    /
  data znu1cs(  5) / +0.117419756889893366664507228352690D-04    /
  data znu1cs(  6) / +0.339357570612261680333825865475121D-06    /
  data znu1cs(  7) / +0.842594176976219910194629891264803D-08    /
  data znu1cs(  8) / +0.183336677024850089184748150900090D-09    /
  data znu1cs(  9) / +0.354969844704416310863007064469557D-11   /
  data znu1cs( 10) / +0.619032496469887332205244342078407D-13   /
  data znu1cs( 11) / +0.981964535680439424960346115456527D-15   /
  data znu1cs( 12) / +0.142851314396490474211473563005985D-16   /
  data znu1cs( 13) / +0.191894921887825298966162467488436D-18   /
  data znu1cs( 14) / +0.239430979739498914162313140597128D-20   /
  data znu1cs( 15) / +0.278890246815347354835870465474995D-22   /
  data znu1cs( 16) / +0.304606650633033442582845214092865D-24   /
  data znu1cs( 17) / +0.313173237042191815771564260932089D-26   /
  data znu1cs( 18) / +0.304133098987854951645174908005034D-28   /
  data znu1cs( 19) / +0.279840384636833084343185097659733D-30   /
  data znu1cs( 20) / +0.244637186274497596485238794922666D-32   /

  data aln2 / 0.69314718055994530941723212145818D+00 /
  data alnbig / 0.0D+00 /
  data alneps / 0.0D+00 /
  data alnsml / 0.0D+00 /
  data euler / 0.57721566490153286060651209008240D+00 /
  data ntc0k / 0 /
  data ntznu1 / 0 /
  data sqpi2 / +1.2533141373155002512078826424055D+00 /
  data xnusml / 0.0D+00 /
  data xsml / 0.0D+00 /

  if ( ntc0k == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    ntc0k = r8_inits ( c0kcs, 29, eta )
    ntznu1 = r8_inits ( znu1cs, 20, eta )
    xnusml = sqrt ( r8_mach ( 3 ) / 8.0D+00 )
    xsml = 0.1D+00 * r8_mach ( 3 )
    alnsml = log ( r8_mach ( 1 ) )
    alnbig = log ( r8_mach ( 2 ) )
    alneps = log ( 0.1D+00 * r8_mach ( 3 ) )
  end if

  if ( xnu < 0.0D+00 .or. 1.0D+00 <= xnu ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_KNUS - Fatal error!'
    write ( *, '(a)' ) '  XNU < 0 or. 1 <= XNU.'
    stop
  end if

  if ( x <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_KNUS - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  end if

  iswtch = 0
!
!  X is small.  Compute k-sub-xnu (x) and the derivative of k-sub-xnu (x)
!  then find k-sub-xnu+1 (x).  xnu is reduced to the interval (-0.5,+0.5)
!  then to (0., .5), because k of negative order (-nu) = k of positive
!  order (+nu).
!
  if ( x <= 2.0D+00 ) then

    if ( xnu <= 0.5D+00 ) then
      v = xnu
    else
      v = 1.0D+00 - xnu
    end if
!
!  Carefully find (x/2)^xnu and z^xnu where z = x*x/4.
!
    alnz = 2.0D+00 * ( log ( x ) - aln2 )

    if ( x <= xnu ) then

      if ( alnbig < - 0.5D+00 * xnu * alnz - aln2 - log ( xnu ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_KNUS - Fatal error!'
        write ( *, '(a)' ) '  Small X causing overflow.'
        stop
      end if

    end if

    vlnz = v * alnz
    x2tov = exp ( 0.5D+00 * vlnz )

    if ( vlnz <= alnsml ) then
      ztov = 0.0D+00
    else
      ztov = x2tov * x2tov
    end if

    a0 = 0.5D+00 * r8_gamma ( 1.0D+00 + v )
    b0 = 0.5D+00 * r8_gamma ( 1.0D+00 - v )
    c0 = - euler
    if ( 0.5D+00 <= ztov .and. xnusml < v ) then
      c0 = - 0.75D+00 + &
        r8_csevl ( ( 8.0D+00 * v ) * v - 1.0D+00, c0kcs, ntc0k )
    end if

    if ( ztov <= 0.5D+00 ) then
      alpha(1) = ( a0 - ztov * b0 ) / v
    else
      alpha(1) = c0 - alnz * ( 0.75D+00 + &
        r8_csevl ( vlnz / 0.35D+00 + 1.0D+00, znu1cs, ntznu1 ) ) * b0
    end if

    beta(1) = - 0.5D+00 * ( a0 + ztov * b0 )

    if ( x <= xsml ) then
      z = 0.0D+00
    else
      z = 0.25D+00 * x * x
    end if

    nterms = max ( 2, int ( 11.0D+00 &
      + ( 8.0D+00 * alnz - 25.19D+00 - alneps ) &
      / ( 4.28D+00 - alnz ) ) )

    do i = 2, nterms
      xi = real ( i - 1, kind = 8 )
      a0 = a0 / ( xi * ( xi - v ) )
      b0 = b0 / ( xi * ( xi + v ) )
      alpha(i) = ( alpha(i-1) + 2.0D+00 * xi * a0 ) / ( xi * ( xi + v ) )
      beta(i) = ( xi - 0.5D+00 * v ) * alpha(i) - ztov * b0
    end do

    bknu = alpha(nterms)
    bknud = beta(nterms)
    do ii = 2, nterms
      i = nterms + 1 - ii
      bknu = alpha(i) + bknu * z
      bknud = beta(i) + bknud * z
    end do

    expx = exp ( x )
    bknu = expx * bknu / x2tov

    if ( alnbig < - 0.5D+00 * ( xnu + 1.0D+00 ) * alnz - 2.0D+00 * aln2 ) then
      iswtch = 1
      return
    end if

    bknud = expx * bknud * 2.0D+00 / ( x2tov * x )

    if ( xnu <= 0.5D+00 ) then
      bknu1 = v * bknu / x - bknud
      return
    end if

    bknu0 = bknu
    bknu = - v * bknu / x - bknud
    bknu1 = 2.0D+00 * xnu * bknu / x + bknu0
!
!  x is large.  find k-sub-xnu (x) and k-sub-xnu+1 (x) with y. l. luke-s
!  rational expansion.
!
  else

    sqrtx = sqrt ( x )

    if ( 1.0D+00 / xsml < x ) then
      bknu = sqpi2 / sqrtx
      bknu1 = bknu
      return
    end if

    an = - 0.60D+00 - 1.02D+00 / x
    bn = - 0.27D+00 - 0.53D+00 / x
    nterms = min ( 32, max ( 3, int ( an + bn * alneps ) ) )

    do inu = 1, 2

      if ( inu == 1 ) then
        if ( xnu <= xnusml ) then
          xmu = 0.0D+00
        else
          xmu = ( 4.0D+00 * xnu ) * xnu
        end if
      else
        xmu = 4.0D+00 * ( abs ( xnu ) + 1.0D+00 ) ** 2
      end if

      a(1) = 1.0D+00 - xmu
      a(2) = 9.0D+00 - xmu
      a(3) = 25.0D+00 - xmu

      if ( a(2) == 0.0D+00 ) then

        result = sqpi2 * ( 16.0D+00 * x + xmu + 7.0D+00 ) &
          / ( 16.0D+00 * x * sqrtx )

      else

        alpha(1) = 1.0D+00
        alpha(2) = ( 16.0D+00 * x + a(2) ) / a(2)
        alpha(3) = ( ( 768.0D+00 * x + 48.0D+00 * a(3) ) * x &
          + a(2) * a(3) ) / ( a(2) * a(3) )

        beta(1) = 1.0D+00
        beta(2) = ( 16.0D+00 * x + ( xmu + 7.0D+00 ) ) / a(2)
        beta(3) = ( ( 768.0D+00 * x &
          + 48.0D+00 * ( xmu + 23.0D+00 ) ) * x + &
          ( ( xmu + 62.0D+00 ) * xmu + 129.0D+00 ) ) &
          / ( a(2) * a(3) )

        do i = 4, nterms

          n = i - 1
          x2n = real ( 2 * n - 1, kind = 8 )

          a(i) = ( x2n + 2.0D+00 ) ** 2 - xmu
          qq = 16.0D+00 * x2n / a(i)
          p1 = - x2n * ( real ( 12 * n * n - 20 * n, kind = 8 ) - a(1) ) &
            / ( ( x2n - 2.0D+00 ) * a(i) ) - qq * x
          p2 = ( real ( 12 * n * n - 28 * n + 8, kind = 8 ) - a(1) ) / a(i) &
            - qq * x
          p3 = - x2n * a(i-3) / ( ( x2n - 2.0D+00 ) * a(i))

          alpha(i) = - p1 * alpha(i-1) &
                     - p2 * alpha(i-2) &
                     - p3 * alpha(i-3)

          beta(i) =  - p1 * beta(i-1) &
                     - p2 * beta(i-2) &
                     - p3 * beta(i-3)

        end do

        result = sqpi2 * beta(nterms) / ( sqrtx * alpha(nterms) )

      end if

      if ( inu == 1 ) then
        bknu = result
      else
        bknu1 = result
      end if

    end do

  end if

  return
end
function r8_lgmc ( x )

!*****************************************************************************80
!
!! R8_LGMC evaluates the log gamma correction factor for an R8 argument.
!
!  Discussion:
!
!    For 10 <= X, compute the log gamma correction factor so that
!
!      log ( gamma ( x ) ) = log ( sqrt ( 2 * pi ) ) 
!                          + ( x - 0.5 ) * log ( x ) - x 
!                          + r8_lgmc ( x )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_LGMC, the correction factor.
!
  implicit none

  real ( kind = 8 ) algmcs(15)
  integer ( kind = 4 ) nalgm
  real ( kind = 8 ) r8_csevl
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_lgmc
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xbig
  real ( kind = 8 ) xmax

  save algmcs
  save nalgm
  save xbig
  save xmax

  data algmcs(  1) / +0.1666389480451863247205729650822D+00 /
  data algmcs(  2) / -0.1384948176067563840732986059135D-04 /
  data algmcs(  3) / +0.9810825646924729426157171547487D-08 /
  data algmcs(  4) / -0.1809129475572494194263306266719D-10 /
  data algmcs(  5) / +0.6221098041892605227126015543416D-13 /
  data algmcs(  6) / -0.3399615005417721944303330599666D-15 /
  data algmcs(  7) / +0.2683181998482698748957538846666D-17 /
  data algmcs(  8) / -0.2868042435334643284144622399999D-19 /
  data algmcs(  9) / +0.3962837061046434803679306666666D-21 /
  data algmcs( 10) / -0.6831888753985766870111999999999D-23 /
  data algmcs( 11) / +0.1429227355942498147573333333333D-24 /
  data algmcs( 12) / -0.3547598158101070547199999999999D-26 /
  data algmcs( 13) / +0.1025680058010470912000000000000D-27 /
  data algmcs( 14) / -0.3401102254316748799999999999999D-29 /
  data algmcs( 15) / +0.1276642195630062933333333333333D-30 /

  data nalgm / 0 /
  data xbig / 0.0D+00 /
  data xmax / 0.0D+00 /

  if ( nalgm == 0 ) then
    nalgm = r8_inits ( algmcs, 15, r8_mach ( 3 ) )
    xbig = 1.0D+00 / sqrt ( r8_mach ( 3 ) )
    xmax = exp ( min ( log ( r8_mach ( 2 ) / 12.0D+00 ), &
      - log ( 12.0D+00 * r8_mach ( 1 ) ) ) )
  end if

  if ( x < 10.0D+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_LGMC - Fatal error!'
    write ( *, '(a)' ) '  X must be at least 10.'
    stop

  else if ( x < xbig ) then

    r8_lgmc = r8_csevl ( 2.0D+00 * ( 10.0D+00 / x ) &
      * ( 10.0D+00 / x ) - 1.0D+00, algmcs, nalgm ) / x

  else if ( x < xmax ) then

    r8_lgmc = 1.0D+00 / ( 12.0D+00 * x )

  else

    r8_lgmc = 0.0D+00

  end if

  return
end
function r8_mach ( i )

!*****************************************************************************80
!
!! R8_MACH returns real ( kind = 8 ) real machine-dependent constants.
!
!  Discussion:
!
!    R8_MACH can be used to obtain machine-dependent parameters
!    for the local machine environment.  It is a function
!    with one input argument, and can be called as follows:
!
!      D = R8_MACH ( I )
!
!    where I=1,...,5.  The output value of D above is
!    determined by the input value of I:.
!
!    R8_MACH ( 1) = B^(EMIN-1), the smallest positive magnitude.
!    R8_MACH ( 2) = B^EMAX*(1 - B^(-T)), the largest magnitude.
!    R8_MACH ( 3) = B^(-T), the smallest relative spacing.
!    R8_MACH ( 4) = B^(1-T), the largest relative spacing.
!    R8_MACH ( 5) = LOG10(B)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528:
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the desired constant.
!
!    Output, real ( kind = 8 ) R8_MACH, the value of the constant.
!
  implicit none

  real ( kind = 8 ) r8_mach
  integer ( kind = 4 ) i

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    r8_mach = 0.0D+00
    stop
  else if ( i == 1 ) then
    r8_mach = 4.450147717014403D-308
  else if ( i == 2 ) then
    r8_mach = 8.988465674311579D+307
  else if ( i == 3 ) then
    r8_mach = 1.110223024625157D-016
  else if ( i == 4 ) then
    r8_mach = 2.220446049250313D-016
  else if ( i == 5 ) then
    r8_mach = 0.301029995663981D+000
  else if ( 5 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    r8_mach = 0.0D+00
    stop
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
subroutine r8mat_cholesky_factor ( n, a, c, flag )

!*****************************************************************************80
!
!! R8MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric matrix.
!
!  Discussion:
!
!    The matrix must be symmetric and positive semidefinite.
!
!    For a positive semidefinite symmetric matrix A, the Cholesky factorization
!    is a lower triangular matrix L such that:
!
!      A = L * L'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 8 ) C(N,N), the N by N lower triangular
!    Cholesky factor.
!
!    Output, integer ( kind = 4 ) FLAG:
!    0, no error occurred.
!    1, the matrix is not positive definite.
!    2, the matrix is not nonnegative definite.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) sum2
  real ( kind = 8 ) tol

  flag = 0
  tol = sqrt ( epsilon ( tol ) )

  c(1:n,1:n) = a(1:n,1:n)

  do j = 1, n

    c(1:j-1,j) = 0.0D+00

    do i = j, n

      sum2 = c(j,i) - dot_product ( c(j,1:j-1), c(i,1:j-1) )

      if ( i == j ) then

        if ( 0.0D+00 < sum2 ) then
          c(i,j) = sqrt ( sum2 )
        else if ( sum2 < - tol ) then
          flag = 2
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8MAT_CHOLESKY_FACTOR - Fatal error!'
          write ( *, '(a)' ) '  Matrix is not nonnegative definite.'
          write ( *, '(a,i4)' ) '  Diagonal I = ', i
          write ( *, '(a,g14.6)' ) '  SUM2 = ', sum2
          return
        else
          flag = 1
          c(i,j) = 0.0D+00
        end if

      else

        if ( c(j,j) /= 0.0D+00 ) then
          c(i,j) = sum2 / c(j,j)
        else
          c(i,j) = 0.0D+00
        end if
      end if

    end do

  end do

  return
end
subroutine r8mat_is_symmetric ( m, n, a, error_frobenius )

!*****************************************************************************80
!
!! R8MAT_IS_SYMMETRIC checks an R8MAT for symmetry.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = 8 ) values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) ERROR_FROBENIUS, measures the 
!    Frobenius norm of ( A - A' ), which would be zero if the matrix
!    were exactly symmetric.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) error_frobenius
  real ( kind = 8 ) r8_huge

  if ( m /= n ) then
    error_frobenius = r8_huge ( )
    return
  end if

  error_frobenius = sqrt ( &
                      sum ( &
                        ( &
                          abs ( a(1:m,1:n) - transpose ( a(1:m,1:n) ) ) &
                         )**2 &
                       ) &
                     )

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
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = 8 ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = 8 ) * a   &
             + real (     i - 1, kind = 8 ) * b ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

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
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
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
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

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
subroutine sample_paths_cholesky ( n, n2, rhomax, rho0, correlation, seed, x )

!*****************************************************************************80
!
!! SAMPLE_PATHS_CHOLESKY: sample paths for stationary correlation functions.
!
!  Discussion:
!
!    This method uses the Cholesky factorization of the correlation matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points on each path.
!
!    Input, integer ( kind = 4 ) N2, the number of paths.
!
!    Input, real ( kind = 8 ) RHOMAX, the maximum value of RHO.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Input, external CORRELATION, the name of the subroutine which evaluates
!    the correlation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X(N,N2), the sample paths.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2

  real ( kind = 8 ) cor(n,n)
  real ( kind = 8 ) cor_vec(n)
  external correlation
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) l(n,n)
  real ( kind = 8 ) r(n,n2)
  real ( kind = 8 ) rho_vec(n)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhomax
  real ( kind = 8 ) rhomin
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n,n2)
!
!  Choose N equally spaced sample points from 0 to RHOMAX.
!
  rhomin = 0.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho_vec )
!
!  Evaluate the correlation function.
!
  call correlation ( n, rho_vec, rho0, cor_vec )
!
!  Construct the correlation matrix;
!
!  From the vector 
!    [ C(0), C(1), C(2), ... C(N-1) ]
!  construct the vector
!    [ C(N-1), ..., C(2), C(1), C(0), C(1), C(2), ...  C(N-1) ]
!  Every row of the correlation matrix can be constructed by a subvector
!  of this vector.
!
  do j = 1, n
    do i = 1, n
      k = i4_wrap ( abs ( j - i ) + 1, 1, n )
      cor(i,j) = cor_vec(k)
    end do
  end do
!
!  Get the Cholesky factorization of COR:
!
!    COR = L * L'.
!
  call r8mat_cholesky_factor ( n, cor, l, flag )
!
!  The matrix might not be nonnegative definite.
!
  if ( flag == 2 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'SAMPLE_PATHS_CHOLESKY - Fatal error!'
    write ( *, '(a)' ) '  The correlation matrix is not'
    write ( *, '(a)' ) '  symmetric nonnegative definite.'
    stop
  end if
!
!  Compute a matrix of N by N2 normally distributed values.
!
  call r8mat_normal_01 ( n, n2, seed, r )
!
!  Compute the sample path.
!
  x(1:n,1:n2) = matmul ( l(1:n,1:n), r(1:n,1:n2) )

  return
end
subroutine sample_paths_eigen ( n, n2, rhomax, rho0, correlation, seed, x )

!*****************************************************************************80
!
!! SAMPLE_PATHS_EIGEN: sample paths for stationary correlation functions.
!
!  Discussion:
!
!    This method uses the eigen-decomposition of the correlation matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points on each path.
!
!    Input, integer ( kind = 4 ) N2, the number of paths.
!
!    Input, real ( kind = 8 ) RHOMAX, the maximum value of RHO.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Input, external CORRELATION, the name of the subroutine which evaluates
!    the correlation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X(N,N2), the sample paths.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2

  real ( kind = 8 ) c(n,n)
  real ( kind = 8 ) cor(n,n)
  real ( kind = 8 ) cor_vec(n)
  external correlation
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) dmin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r(n,n2)
  real ( kind = 8 ) rho_vec(n)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhomax
  real ( kind = 8 ) rhomin
  integer ( kind = 4 ) seed
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n,n2)
!
!  Choose N equally spaced sample points from 0 to RHOMAX.
!
  rhomin = 0.0D+00
  call r8vec_linspace ( n, rhomin, rhomax, rho_vec )
!
!  Evaluate the correlation function.
!
  call correlation ( n, rho_vec, rho0, cor_vec )
!
!  Construct the correlation matrix;
!
!  From the vector 
!    [ C(0), C(1), C(2), ... C(N-1) ]
!  construct the vector
!    [ C(N-1), ..., C(2), C(1), C(0), C(1), C(2), ...  C(N-1) ]
!  Every row of the correlation matrix can be constructed by a subvector
!  of this vector.
!
  do j = 1, n
    do i = 1, n
      k = i4_wrap ( abs ( i - j ), 0, n - 1 ) + 1
      cor(i,j) = cor_vec(k)
    end do
  end do
!
!  Get the eigendecomposition of COR:
!
!    COR = V * D * V'.
!
!  Because COR is symmetric, V is orthogonal.
!
  call tred2 ( n, cor, d, w, v )

  call tql2 ( n, d, w, v, ierr )
!
!  We assume COR is non-negative definite, and hence that there
!  are no negative eigenvalues.
!
  dmin = minval ( d )

  if ( dmin < - sqrt ( epsilon ( dmin ) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SAMPLE_PATHS_EIGEN - Warning!'
    write ( *, '(a,g14.6)' ) '  Negative eigenvalues observed as low as ', dmin
  end if

  do i = 1, n
    d(i) = max ( d(i), 0.0D+00 )
  end do
!
!  Compute the eigenvalues of the factor C.
!
  d(1:n) = sqrt ( d(1:n) )
!
!  Compute C, such that C' * C = COR.
!
  do j = 1, n
    do i = 1, n
      c(i,j) = 0.0D+00
      do k = 1, n
        c(i,j) = c(i,j) + d(k) * v(i,k) * v(j,k)
      end do
    end do
  end do
!
!  Compute N by N2 independent random normal values.
!
  call r8mat_normal_01 ( n, n2, seed, r )
!
!  Multiply to get the variables X which have correlation COR.
!
  x(1:n,1:n2) = matmul ( c(1:n,1:n), r(1:n,1:n2) )

  return
end
subroutine sample_paths2_cholesky ( n, n2, rhomin, rhomax, rho0, correlation2, &
  seed, x )

!*****************************************************************************80
!
!! SAMPLE_PATHS2_CHOLESKY: sample paths for stationary correlation functions.
!
!  Discussion:
!
!    This method uses the Cholesky factorization of the correlation matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points on each path.
!
!    Input, integer ( kind = 4 ) N2, the number of paths.
!
!    Input, real ( kind = 8 ) RHOMIN, RHOMAX, the range of RHO.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Input, external CORRELATION2, the name of the subroutine which evaluates
!    the correlation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X(N,N2), the sample paths.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2

  real ( kind = 8 ) cor(n,n)
  external correlation2
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) l(n,n)
  real ( kind = 8 ) r(n,n2)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhomax
  real ( kind = 8 ) rhomin
  real ( kind = 8 ) s(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n,n2)
!
!  Choose N equally spaced sample points from RHOMIN to RHOMAX.
!
  call r8vec_linspace ( n, rhomin, rhomax, s )
!
!  Evaluate the correlation matrix.
!
  call correlation2 ( n, n, s, s, rho0, cor )
!
!  Get the Cholesky factorization of COR:
!
!    COR = L * L'.
!
  call r8mat_cholesky_factor ( n, cor, l, flag )
!
!  The matrix might not be nonnegative definite.
!
  if ( flag == 2 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'SAMPLE_PATHS2_CHOLESKY - Fatal error!'
    write ( *, '(a)' ) '  The correlation matrix is not'
    write ( *, '(a)' ) '  symmetric nonnegative definite.'
    stop
  end if
!
!  Compute a matrix of N by N2 normally distributed values.
!
  call r8mat_normal_01 ( n, n2, seed, r )
!
!  Compute the sample path.
!
  x(1:n,1:n2) = matmul ( l(1:n,1:n), r(1:n,1:n2) )

  return
end
subroutine sample_paths2_eigen ( n, n2, rhomin, rhomax, rho0, correlation2, &
  seed, x )

!*****************************************************************************80
!
!! SAMPLE_PATHS2_EIGEN: sample paths for stationary correlation functions.
!
!  Discussion:
!
!    This method uses the eigen-decomposition of the correlation matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points on each path.
!
!    Input, integer ( kind = 4 ) N2, the number of paths.
!
!    Input, real ( kind = 8 ) RHOMIN, RHOMAX, the range of RHO.
!
!    Input, real ( kind = 8 ) RHO0, the correlation length.
!
!    Input, external CORRELATION2, the name of the subroutine which evaluates
!    the correlation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X(N,N2), the sample paths.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2

  real ( kind = 8 ) c(n,n)
  real ( kind = 8 ) cor(n,n)
  external correlation2
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) dmin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r(n,n2)
  real ( kind = 8 ) rho0
  real ( kind = 8 ) rhomax
  real ( kind = 8 ) rhomin
  real ( kind = 8 ) s(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n,n2)
!
!  Choose N equally spaced sample points from RHOMIN to RHOMAX.
!
  call r8vec_linspace ( n, rhomin, rhomax, s )
!
!  Evaluate the correlation function.
!
  call correlation2 ( n, n, s, s, rho0, cor )
!
!  Get the eigendecomposition of COR:
!
!    COR = V * D * V'.
!
!  Because COR is symmetric, V is orthogonal.
!
  call tred2 ( n, cor, d, w, v )

  call tql2 ( n, d, w, v, ierr )
!
!  We assume COR is non-negative definite, and hence that there
!  are no negative eigenvalues.
!
  dmin = minval ( d )

  if ( dmin < - sqrt ( epsilon ( dmin ) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SAMPLE_PATHS2_EIGEN - Warning!'
    write ( *, '(a,g14.6)' ) '  Negative eigenvalues observed as low as ', dmin
  end if

  do i = 1, n
    d(i) = max ( d(i), 0.0D+00 )
  end do
!
!  Compute the eigenvalues of the factor C.
!
  d(1:n) = sqrt ( d(1:n) )
!
!  Compute C, such that C' * C = COR.
!
  do j = 1, n
    do i = 1, n
      c(i,j) = 0.0D+00
      do k = 1, n
        c(i,j) = c(i,j) + d(k) * v(i,k) * v(j,k)
      end do
    end do
  end do
!
!  Compute N by N2 independent random normal values.
!
  call r8mat_normal_01 ( n, n2, seed, r )
!
!  Multiply to get the variables X which have correlation COR.
!
  x(1:n,1:n2) = matmul ( c(1:n,1:n), r(1:n,1:n2) )

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
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
subroutine tql2 ( n, d, e, z, ierr )

!*****************************************************************************80
!
!! TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This subroutine finds the eigenvalues and eigenvectors of a symmetric
!    tridiagonal matrix by the QL method.  The eigenvectors of a full
!    symmetric matrix can also be found if TRED2 has been used to reduce this
!    full matrix to tridiagonal form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Bowdler, Martin, Reinsch, Wilkinson,
!    TQL2,
!    Numerische Mathematik,
!    Volume 11, pages 293-306, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal elements of
!    the matrix.  On output, the eigenvalues in ascending order.  If an error
!    exit is made, the eigenvalues are correct but unordered for indices
!    1,2,...,IERR-1.
!
!    Input/output, real ( kind = 8 ) E(N).  On input, E(2:N) contains the
!    subdiagonal elements of the input matrix, and E(1) is arbitrary.
!    On output, E has been destroyed.
!
!    Input, real ( kind = 8 ) Z(N,N).  On input, the transformation matrix
!    produced in the reduction by TRED2, if performed.  If the eigenvectors of
!    the tridiagonal matrix are desired, Z must contain the identity matrix.
!    On output, Z contains the orthonormal eigenvectors of the symmetric
!    tridiagonal (or full) matrix.  If an error exit is made, Z contains
!    the eigenvectors associated with the stored eigenvalues.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, normal return,
!    J, if the J-th eigenvalue has not been determined after
!    30 iterations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) dl1
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) el1
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) pythag
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) s2
  real ( kind = 8 ) t
  real ( kind = 8 ) tst1
  real ( kind = 8 ) tst2
  real ( kind = 8 ) z(n,n)

  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do

  f = 0.0D+00
  tst1 = 0.0D+00
  e(n) = 0.0D+00

  do l = 1, n

    j = 0
    h = abs ( d(l) ) + abs ( e(l) )
    tst1 = max ( tst1, h )
!
!  Look for a small sub-diagonal element.
!
    do m = l, n
      tst2 = tst1 + abs ( e(m) )
      if ( tst2 == tst1 ) then
        exit
      end if
    end do

    if ( m /= l ) then

      do

        if ( 30 <= j ) then
          ierr = l
          return
        end if

        j = j + 1
!
!  Form shift.
!
        l1 = l + 1
        l2 = l1 + 1
        g = d(l)
        p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
        r = pythag ( p, 1.0D+00 )
        d(l) = e(l) / ( p + sign ( r, p ) )
        d(l1) = e(l) * ( p + sign ( r, p ) )
        dl1 = d(l1)
        h = g - d(l)
        d(l2:n) = d(l2:n) - h
        f = f + h
!
!  QL transformation.
!
        p = d(m)
        c = 1.0D+00
        c2 = c
        el1 = e(l1)
        s = 0.0D+00
        mml = m - l

        do ii = 1, mml

          c3 = c2
          c2 = c
          s2 = s
          i = m - ii
          g = c * e(i)
          h = c * p
          r = pythag ( p, e(i) )
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * ( c * g + s * d(i) )
!
!  Form vector.
!
          do k = 1, n
            h = z(k,i+1)
            z(k,i+1) = s * z(k,i) + c * h
            z(k,i) = c * z(k,i) - s * h
          end do

        end do

        p = - s * s2 * c3 * el1 * e(l) / dl1
        e(l) = s * p
        d(l) = c * p
        tst2 = tst1 + abs ( e(l) )

        if ( tst2 <= tst1 ) then
          exit
        end if

      end do

    end if

    d(l) = d(l) + f

  end do
!
!  Order eigenvalues and eigenvectors.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n

      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if

    end do

    if ( k /= i ) then

      d(k) = d(i)
      d(i) = p

      do j = 1, n
        t      = z(j,i)
        z(j,i) = z(j,k)
        z(j,k) = t
      end do

    end if

  end do

  return
end
subroutine tred2 ( n, a, d, e, z )

!*****************************************************************************80
!
!! TRED2 transforms a real symmetric matrix to symmetric tridiagonal form.
!
!  Discussion:
!
!    This subroutine reduces a real symmetric matrix to a
!    symmetric tridiagonal matrix using and accumulating
!    orthogonal similarity transformations.
!
!    A and Z may coincide, in which case a single storage area is used
!    for the input of A and the output of Z.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Martin, Reinsch, Wilkinson,
!    TRED2,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the real symmetric input matrix.  Only the
!    lower triangle of the matrix need be supplied.
!
!    Output, real ( kind = 8 ) D(N), the diagonal elements of the tridiagonal
!    matrix.
!
!    Output, real ( kind = 8 ) E(N), contains the subdiagonal elements of the
!    tridiagonal matrix in E(2:N).  E(1) is set to zero.
!
!    Output, real ( kind = 8 ) Z(N,N), the orthogonal transformation matrix
!    produced in the reduction.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  real ( kind = 8 ) hh
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) scale
  real ( kind = 8 ) z(n,n)

  do i = 1, n
    z(i:n,i) = a(i:n,i)
  end do

  d(1:n) = a(n,1:n)

  do ii = 2, n

    i = n + 2 - ii
    l = i - 1
    h = 0.0D+00
!
!  Scale row.
!
    scale = sum ( abs ( d(1:l) ) )

    if ( scale == 0.0D+00 ) then

      e(i) = d(l)

      do j = 1, l
        d(j) = z(l,j)
        z(i,j) = 0.0D+00
        z(j,i) = 0.0D+00
      end do

      d(i) = 0.0D+00

      cycle

    end if

    d(1:l) = d(1:l) / scale

    h = h + dot_product ( d(1:l), d(1:l) )

    f = d(l)
    g = - sign ( sqrt ( h ), f )
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g
!
!  Form A*U.
!
    e(1:l) = 0.0D+00

    do j = 1, l

      f = d(j)
      z(j,i) = f
      g = e(j) + z(j,j) * f

      do k = j + 1, l
        g = g + z(k,j) * d(k)
        e(k) = e(k) + z(k,j) * f
      end do

      e(j) = g

    end do
!
!  Form P.
!
    e(1:l) = e(1:l) / h

    f = dot_product ( e(1:l), d(1:l) )

    hh = 0.5D+00 * f / h
!
!  Form Q.
!
    e(1:l) = e(1:l) - hh * d(1:l)
!
!  Form reduced A.
!
    do j = 1, l

      f = d(j)
      g = e(j)

      z(j:l,j) = z(j:l,j) - f * e(j:l) - g * d(j:l)

      d(j) = z(l,j)
      z(i,j) = 0.0D+00

    end do

    d(i) = h

  end do
!
!  Accumulation of transformation matrices.
!
  do i = 2, n

    l = i - 1
    z(n,l) = z(l,l)
    z(l,l) = 1.0D+00
    h = d(i)

    if ( h /= 0.0D+00 ) then

      d(1:l) = z(1:l,i) / h

      do j = 1, l

        g = dot_product ( z(1:l,i), z(1:l,j) )

        do k = 1, l
          z(k,j) = z(k,j) - g * d(k)
        end do

      end do

    end if

    z(1:l,i) = 0.0D+00

  end do

  d(1:n) = z(n,1:n)

  z(n,1:n-1) = 0.0D+00
  z(n,n) = 1.0D+00

  e(1) = 0.0D+00

  return
end
