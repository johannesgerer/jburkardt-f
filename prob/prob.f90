subroutine angle_cdf ( x, n, cdf )

!*****************************************************************************80
!
!! ANGLE_CDF evaluates the Angle CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation and Sensitivity of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N must be at least 2.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) n
  real ( kind = 8 ) n_real
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) sin_power_int
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ANGLE_CDF - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop
  end if

  if ( x <= 0.0D+00 ) then
    cdf = 0.0D+00
  else if ( pi <= x ) then
    cdf = 1.0D+00
  else if ( n == 2 ) then
    cdf = x / pi
  else
    n_real = real ( n, kind = 8 )
    cdf = sin_power_int ( zero, x, n - 2 ) * r8_gamma ( n_real / 2.0D+00 ) &
      / ( sqrt ( pi ) * r8_gamma ( ( n_real - 1.0D+00 ) / 2.0D+00 ) )
  end if

  return
end
subroutine angle_mean ( n, mean )

!*****************************************************************************80
!
!! ANGLE_MEAN returns the mean of the Angle PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N must be at least 2.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) mean
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  mean = pi / 2.0D+00

  return
end
subroutine angle_pdf ( x, n, pdf )

!*****************************************************************************80
!
!! ANGLE_PDF evaluates the Angle PDF.
!
!  Discussion:
!
!    X is an angle between 0 and PI, corresponding to the angle
!    made in an N-dimensional space, between a fixed line passing
!    through the origin, and an arbitrary line that also passes
!    through the origin, which is specified by a choosing any point
!    on the N-dimensional sphere with uniform probability.
!
!    The formula is
!
!      PDF(X) = ( sin ( X ) )**(N-2) * Gamma ( N / 2 )
!               / ( sqrt ( PI ) * Gamma ( ( N - 1 ) / 2 ) )
!
!      PDF(X) = 1 / PI if N = 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation and Sensitivity of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N must be at least 2.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) x

  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ANGLE_PDF - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop
  end if

  if ( x < 0.0D+00 .or. pi < x ) then
    pdf = 0.0D+00
  else if ( n == 2 ) then
    pdf = 1.0D+00 / pi
  else
    pdf = ( sin ( x ) )**( n - 2 ) * r8_gamma ( real ( n, kind = 8 ) / 2.0D+00 ) &
      / ( sqrt ( pi ) * r8_gamma ( real ( n - 1, kind = 8 ) / 2.0D+00 ) )
  end if

  return
end
subroutine anglit_cdf ( x, cdf )

!*****************************************************************************80
!
!! ANGLIT_CDF evaluates the Anglit CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( x <  - 0.25D+00 * pi ) then
    cdf = 0.0D+00
  else if ( x < 0.25D+00 * pi ) then
    cdf = 0.5D+00 - 0.5D+00 * cos ( 2.0D+00 * x + pi / 2.0D+00 )
  else
    cdf = 1.0D+00
  end if

  return
end
subroutine anglit_cdf_inv ( cdf, x )

!*****************************************************************************80
!
!! ANGLIT_CDF_INV inverts the Anglit CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ANGLIT_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = 0.5D+00 * ( acos ( 1.0D+00 - 2.0D+00 * cdf ) - pi / 2.0D+00 )

  return
end
subroutine anglit_mean ( mean )

!*****************************************************************************80
!
!! ANGLIT_MEAN returns the mean of the Anglit PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) mean

  mean = 0.0D+00

  return
end
subroutine anglit_pdf ( x, pdf )

!*****************************************************************************80
!
!! ANGLIT_PDF evaluates the Anglit PDF.
!
!  Discussion:
!
!    PDF(X) = sin ( 2 * X + PI / 2 ) for -PI/4 <= X <= PI/4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( x <= - 0.25D+00 * pi .or. 0.25D+00 * pi <= x ) then
    pdf = 0.0D+00
  else
    pdf = sin ( 2.0D+00 * x + 0.25D+00 * pi )
  end if

  return
end
subroutine anglit_sample ( seed, x )

!*****************************************************************************80
!
!! ANGLIT_SAMPLE samples the Anglit PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 1999
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
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call anglit_cdf_inv ( cdf, x )

  return
end
subroutine anglit_variance ( variance )

!*****************************************************************************80
!
!! ANGLIT_VARIANCE returns the variance of the Anglit PDF.
!
!  Discussion:
!
!    Variance =
!      Integral ( -PI/4 <= X <= PI/4 ) X^2 * sin ( 2 * X + PI / 2 )
!
!    Antiderivative =
!      0.5D+00 * X * sin ( 2 * X + PI / 2 )
!      + ( 0.25 - 0.5D+00 * X^2 ) * cos ( 2 * X + PI / 2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) variance

  variance = 0.0625D+00 * pi * pi - 0.5D+00

  return
end
subroutine arcsin_cdf ( x, a, cdf )

!*****************************************************************************80
!
!! ARCSIN_CDF evaluates the Arcsin CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, the parameter of the CDF.
!    A must be positive.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( x <= -a ) then
    cdf = 0.0D+00
  else if ( x < a ) then
    cdf = 0.5D+00 + asin ( x / a ) / pi
  else
    cdf = 1.0D+00
  end if

  return
end
subroutine arcsin_cdf_inv ( cdf, a, x )

!*****************************************************************************80
!
!! ARCSIN_CDF_INV inverts the Arcsin CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, the parameter of the CDF.
!    A must be positive.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ARCSIN_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a * sin ( pi * ( cdf - 0.5D+00 ) )

  return
end
function arcsin_check ( a )

!*****************************************************************************80
!
!! ARCSIN_CHECK checks the parameter of the Arcsin CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0 < A.
!
!    Output, logical ARCSIN_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  logical arcsin_check

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ARCSIN_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    arcsin_check = .false.
    return
  end if

  arcsin_check = .true.

  return
end
subroutine arcsin_mean ( a, mean )

!*****************************************************************************80
!
!! ARCSIN_MEAN returns the mean of the Arcsin PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the CDF.
!    A must be positive.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) mean

  mean = 0.0D+00

  return
end
subroutine arcsin_pdf ( x, a, pdf )

!*****************************************************************************80
!
!! ARCSIN_PDF evaluates the Arcsin PDF.
!
!  Discussion:
!
!    The LOGISTIC EQUATION has the form:
!
!      X(N+1) = 4.0D+00 * LAMBDA * ( 1.0D+00 - X(N) ).
!
!    where 0 < LAMBDA <= 1.  This nonlinear difference equation maps
!    the unit interval into itself, and is a simple example of a system
!    exhibiting chaotic behavior.  Ulam and von Neumann studied the
!    logistic equation with LAMBDA = 1, and showed that iterates of the
!    function generated a sequence of pseudorandom numbers with
!    the Arcsin probability density function.
!
!    The derived sequence
!
!      Y(N) = ( 2 / PI ) * Arcsin ( SQRT ( X(N) ) )
!
!    is a pseudorandom sequence with the uniform probability density
!    function on [0,1].  For certain starting values, such as X(0) = 0, 0.75,
!    or 1.0D+00, the sequence degenerates into a constant sequence, and for
!    values very near these, the sequence takes a while before becoming
!    chaotic.
!
!    The formula is:
!
!      PDF(X) = 1 / ( pi * sqrt ( A^2 - X^2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, Stephen Kokoska,
!    CRC Standard Probability and Statistics Tables and Formulae,
!    Chapman and Hall/CRC, 2000, pages 114-115.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    -A < X < A.
!
!    Input, real ( kind = 8 ) A, the parameter of the CDF.
!    A must be positive.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ARCSIN_PDF - Fatal error!'
    write ( *, '(a)' ) '  Parameter A must be positive.'
    stop
  end if

  if ( x <= -a .or. a <= x ) then
    pdf = 0.0D+00
  else
    pdf = 1.0D+00 / ( pi * sqrt ( a * a - x * x ) )
  end if

  return
end
subroutine arcsin_sample ( a, seed, x )

!*****************************************************************************80
!
!! ARCSIN_SAMPLE samples the Arcsin PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the CDF.
!    A must be positive.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call arcsin_cdf_inv ( cdf, a, x )

  return
end
subroutine arcsin_variance ( a, variance )

!*****************************************************************************80
!
!! ARCSIN_VARIANCE returns the variance of the Arcsin PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the CDF.
!    A must be positive.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) variance

  variance = a * a / 2.0D+00

  return
end
subroutine benford_pdf ( x, pdf )

!*****************************************************************************80
!
!! BENFORD_PDF returns the Benford PDF.
!
!  Discussion:
!
!    Benford's law is an empirical formula explaining the observed
!    distribution of initial digits in lists culled from newspapers,
!    tax forms, stock market prices, and so on.  It predicts the observed
!    high frequency of the initial digit 1, for instance.
!
!    Note that the probabilities of digits 1 through 9 are guaranteed
!    to add up to 1, since
!      LOG10 ( 2/1 ) + LOG10 ( 3/2) + LOG10 ( 4/3 ) + ... + LOG10 ( 10/9 )
!      = LOG10 ( 2/1 * 3/2 * 4/3 * ... * 10/9 ) = LOG10 ( 10 ) = 1.
!
!    The formula is:
!
!      PDF(X) = LOG10 ( ( X + 1 ) / X ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Frank Benford,
!    The Law of Anomalous Numbers,
!    Proceedings of the American Philosophical Society,
!    Volume 78, pages 551-572, 1938.
!
!    Ted Hill,
!    The First Digit Phenomenon,
!    American Scientist,
!    Volume 86, July/August 1998, pages 358 - 363.
!
!    Ralph Raimi,
!    The Peculiar Distribution of First Digits,
!    Scientific American,
!    December 1969, pages 109-119.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the string of significant digits to be checked.
!    If X is 1, then we are asking for the Benford probability that
!    a value will have first digit 1.  If X is 123, we are asking for
!    the probability that the first three digits will be 123, and so on.
!
!    Output, real ( kind = 8 ) PDF, the Benford probability that an item taken
!    from a real world distribution will have the initial digits X.
!
  implicit none

  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  if ( x <= 0 ) then
    pdf = 0.0D+00
  else
    pdf = log10 ( real ( x + 1, kind = 8 ) / real ( x, kind = 8 ) )
  end if

  return
end
subroutine bessel_ix_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_IX_VALUES returns some values of the Ix Bessel function.
!
!  Discussion:
!
!    This set of data considers the less common case in which the
!    index of the Bessel function In is actually not an integer.
!    We may suggest this case by occasionally replacing the symbol
!    "In" by "Ix".
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselI[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.3592084175833614D+00,  &
    0.9376748882454876D+00,  &
    2.046236863089055D+00,   &
    3.053093538196718D+00,   &
    4.614822903407601D+00,   &
    26.47754749755907D+00,   &
    2778.784603874571D+00,   &
    4.327974627242893D+07,   &
    0.2935253263474798D+00,  &
    1.099473188633110D+00,   &
    21.18444226479414D+00,   &
    2500.906154942118D+00,   &
    2.866653715931464D+20,   &
    0.05709890920304825D+00, &
    0.3970270801393905D+00,  &
    13.76688213868258D+00,   &
    2028.512757391936D+00,   &
    2.753157630035402D+20,   &
    0.4139416015642352D+00,  &
    1.340196758982897D+00,   &
    22.85715510364670D+00,   &
    2593.006763432002D+00,   &
    2.886630075077766D+20,   &
    0.03590910483251082D+00, &
    0.2931108636266483D+00,  &
    11.99397010023068D+00,   &
    1894.575731562383D+00,   &
    2.716911375760483D+20 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) nu
  real ( kind = 8 ), save, dimension ( n_max ) :: nu_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     1.0D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     5.0D+00, &
    10.0D+00, &
    20.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine birthday_cdf ( n, cdf )

!*****************************************************************************80
!
!! BIRTHDAY_CDF returns the Birthday Concurrence CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of people whose birthdays have
!    been disclosed.
!
!    Output, real ( kind = 8 ) CDF, the probability that at least
!    two of the N people have matching birthays..
!
  implicit none

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  if ( n < 1 ) then
    cdf = 0.0D+00
    return
  else if ( 365 < n ) then
    cdf = 1.0D+00
    return
  end if
!
!  Compute the probability that N people have distinct birthdays.
!
  cdf = 1.0D+00
  do i = 1, n
    cdf = cdf * real ( 365 + 1 - i, kind = 8 ) / 365.0D+00
  end do
!
!  Compute the probability that it is NOT the case that N people
!  have distinct birthdays.  This is the cumulative probability
!  that person 2 matches person 1, or person 3 matches 1 or 2,
!  etc.
!
  cdf = 1.0D+00 - cdf

  return
end
subroutine birthday_cdf_inv ( cdf, n )

!*****************************************************************************80
!
!! BIRTHDAY_CDF_INV inverts the Birthday Concurrence CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the probability that at least
!    two of the N people have matching birthays.
!
!    Output, integer ( kind = 4 ) N, the corresponding number of people whose
!    birthdays need to be disclosed.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf_not
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  if ( cdf <= 0.0D+00 ) then
    n = 0
    return
  else if ( 1.0D+00 <= cdf ) then
    n = 365
    return
  end if
!
!  Compute the probability that N people have distinct birthdays.
!
  cdf_not = 1.0D+00

  do i = 1, 365
    cdf_not = cdf_not * real ( 365 + 1 - i, kind = 8 ) / 365.0D+00
    if ( cdf <= 1.0D+00 - cdf_not ) then
      n = i
      return
    end if
  end do

  n = 365

  return
end
subroutine birthday_pdf ( n, pdf )

!*****************************************************************************80
!
!! BIRTHDAY_PDF returns the Birthday Concurrence PDF.
!
!  Discussion:
!
!    The probability is the probability that the N-th person is the
!    first one to match a birthday with someone earlier.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of people whose birthdays have
!    been disclosed.
!
!    Output, real ( kind = 8 ) PDF, the probability that the N-th person
!    is the first to match a birthday with someone earlier.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) pdf

  if ( n < 1 .or. 365 < n ) then
    pdf = 0.0D+00
    return
  end if

  pdf = 1.0D+00
!
!  Compute the probability that N-1 people have distinct birthdays.
!
  do i = 1, n-1
    pdf = pdf * real ( 365 + 1 - i, kind = 8 ) / 365.0D+00
  end do
!
!  Compute the probability that person N has one of those N-1 birthdays.
!
  pdf = pdf * real ( n - 1, kind = 8 ) / 365.0D+00

  return
end
subroutine bernoulli_cdf ( x, a, cdf )

!*****************************************************************************80
!
!! BERNOULLI_CDF evaluates the Bernoulli CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the number of successes on a single trial.
!    X = 0 or 1.
!
!    Input, real ( kind = 8 ) A, the probability of success on one trial.
!    0.0D+00 <= A <= 1.0.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    cdf = 0.0D+00
  else if ( x == 0 ) then
    cdf = 1.0D+00 - a
  else
    cdf = 1.0D+00
  end if

  return
end
subroutine bernoulli_cdf_inv ( cdf, a, x )

!*****************************************************************************80
!
!! BERNOULLI_CDF_INV inverts the Bernoulli CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 <= A <= 1.0.
!
!    Output, integer ( kind = 4 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BERNOULLI_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf <= 1.0D+00 - a ) then
    x = 0
  else
    x = 1
  end if

  return
end
function bernoulli_check ( a )

!*****************************************************************************80
!
!! BERNOULLI_CHECK checks the parameter of the Bernoulli CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0 <= A <= 1.0.
!
!    Output, logical BERNOULLI_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  logical bernoulli_check

  if ( a < 0.0D+00 .or. 1.0D+00 < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BERNOULLI_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A < 0 or 1 < A.'
    bernoulli_check = .false.
    return
  end if

  bernoulli_check = .true.

  return
end
subroutine bernoulli_mean ( a, mean )

!*****************************************************************************80
!
!! BERNOULLI_MEAN returns the mean of the Bernoulli PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the probability of success.
!    0.0D+00 <= A <= 1.0.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine bernoulli_pdf ( x, a, pdf )

!*****************************************************************************80
!
!! BERNOULLI_PDF evaluates the Bernoulli PDF.
!
!  Discussion:
!
!    PDF(A;X) = A^X * ( 1.0D+00 - A )^( X - 1 )
!
!    X = 0 or 1.
!
!    The Bernoulli PDF describes the simple case in which a single trial
!    is carried out, with two possible outcomes, called "success" and
!    "failure"; the probability of success is A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the number of successes on a single trial.
!    X = 0 or 1.
!
!    Input, real ( kind = 8 ) A, the probability of success on one trial.
!    0.0 <= A <= 1.0.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    pdf = 0.0D+00
  else if ( x == 0 ) then
    pdf = 1.0D+00 - a
  else if ( x == 1 ) then
    pdf = a
  else
    pdf = 0.0D+00
  end if

  return
end
subroutine bernoulli_sample ( a, seed, x )

!*****************************************************************************80
!
!! BERNOULLI_SAMPLE samples the Bernoulli PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the probability of success on one trial.
!    0.0D+00 <= A <= 1.0.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x

  cdf = r8_uniform_01 ( seed )

  call bernoulli_cdf_inv ( cdf, a, x )

  return
end
subroutine bernoulli_variance ( a, variance )

!*****************************************************************************80
!
!! BERNOULLI_VARIANCE returns the variance of the Bernoulli PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the probability of success on one trial.
!    0.0D+00 <= A <= 1.0.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) variance

  variance = a * ( 1.0D+00 - a )

  return
end
function bessel_i0 ( arg )

!*****************************************************************************80
!
!! BESSEL_I0 evaluates the modified Bessel function I0(X).
!
!  Discussion:
!
!    The main computation evaluates slightly modified forms of
!    minimax approximations generated by Blair and Edwards, Chalk
!    River (Atomic Energy of Canada Limited) Report AECL-4928,
!    October, 1974.  This transportable program is patterned after
!    the machine dependent FUNPACK packet NATSI0, but cannot match
!    that version for efficiency or accuracy.  This version uses
!    rational functions that theoretically approximate I-SUB-0(X)
!    to at least 18 significant decimal digits.
!
!  Machine dependent constants:
!
!    beta   = Radix for the floating-point system
!    maxexp = Smallest power of beta that overflows
!    XMAX =   Largest argument acceptable to BESI0;  Solution to
!             equation:
!               W(X) * (1+1/(8*X)+9/(128*X^2) = beta^maxexp
!             where  W(X) = EXP(X)/sqrt(2*PI*X)
!
!    Approximate values for some important machines are:
!
!                             beta       maxexp       XMAX
!
!    CRAY-1        (S.P.)       2         8191       5682.810
!    Cyber 180/855
!      under NOS   (S.P.)       2         1070        745.893
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)       2          128         91.900
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)       2         1024        713.986
!    IBM 3033      (D.P.)      16           63        178.182
!    VAX           (S.P.)       2          127         91.203
!    VAX D-Format  (D.P.)       2          127         91.203
!    VAX G-Format  (D.P.)       2         1023        713.293
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.
!
!    Output, real ( kind = 8 ) BESSEL_I0, the value of the modified
!    Bessel function of the first kind.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) bessel_i0
  real ( kind = 8 ), parameter :: exp40 = 2.353852668370199854D+17
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter, dimension ( 15 ) :: p = (/ &
    -5.2487866627945699800D-18, &
    -1.5982226675653184646D-14, &
    -2.6843448573468483278D-11, &
    -3.0517226450451067446D-08, &
    -2.5172644670688975051D-05, &
    -1.5453977791786851041D-02, &
    -7.0935347449210549190D+00, &
    -2.4125195876041896775D+03, &
    -5.9545626019847898221D+05, &
    -1.0313066708737980747D+08, &
    -1.1912746104985237192D+10, &
    -8.4925101247114157499D+11, &
    -3.2940087627407749166D+13, &
    -5.5050369673018427753D+14, &
    -2.2335582639474375249D+15 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: pp = (/ &
    -3.9843750000000000000D-01, &
     2.9205384596336793945D+00, &
    -2.4708469169133954315D+00, &
     4.7914889422856814203D-01, &
    -3.7384991926068969150D-03, &
    -2.6801520353328635310D-03, &
     9.9168777670983678974D-05, &
    -2.1877128189032726730D-06 /)
  real ( kind = 8 ), parameter, dimension ( 5 ) :: q = (/ &
    -3.7277560179962773046D+03, &
     6.5158506418655165707D+06, &
    -6.5626560740833869295D+09, &
     3.7604188704092954661D+12, &
    -9.7087946179594019126D+14 /)
  real ( kind = 8 ), parameter, dimension ( 7 ) :: qq = (/ &
    -3.1446690275135491500D+01, &
     8.5539563258012929600D+01, &
    -6.0228002066743340583D+01, &
     1.3982595353892851542D+01, &
    -1.1151759188741312645D+00, &
     3.2547697594819615062D-02, &
    -5.5194330231005480228D-04 /)
  real ( kind = 8 ), parameter :: rec15 = 6.6666666666666666666D-02
  real ( kind = 8 ) sump
  real ( kind = 8 ) sumq
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xmax = 91.9D+00
  real ( kind = 8 ) xx

  x = abs ( arg )

  if ( x < epsilon ( arg ) ) then
    value = 1.0D+00
  else if ( x < 15.0D+00 ) then
!
!  EPSILON ( ARG ) <= ABS(ARG) < 15.0D+00
!
    xx = x * x
    sump = p(1)
    do i = 2, 15
      sump = sump * xx + p(i)
    end do

    xx = xx - 225.0D+00
    sumq = (((( &
           xx + q(1) ) &
         * xx + q(2) ) &
         * xx + q(3) ) &
         * xx + q(4) ) &
         * xx + q(5)

    value = sump / sumq

  else if ( 15.0D+00 <= x ) then

    if ( xmax < x ) then
      value = huge ( value )
    else
!
!  15.0D+00 <= ABS(ARG)
!
      xx = 1.0D+00 / x - rec15

      sump = ((((((  &
                  pp(1) &
           * xx + pp(2) ) &
           * xx + pp(3) ) &
           * xx + pp(4) ) &
           * xx + pp(5) ) &
           * xx + pp(6) ) &
           * xx + pp(7) ) &
           * xx + pp(8)

      sumq = (((((( &
             xx + qq(1) ) &
           * xx + qq(2) ) &
           * xx + qq(3) ) &
           * xx + qq(4) ) &
           * xx + qq(5) ) &
           * xx + qq(6) ) &
           * xx + qq(7)

      value = sump / sumq
!
!  Calculation reformulated to avoid premature overflow.
!
      if ( x <= xmax - 15.0D+00 ) then
        a = exp ( x )
        b = 1.0D+00
      else
        a = exp ( x - 40.0D+00 )
        b = exp40
      end if

      value = ( ( value * a - pp(1) * a ) / sqrt ( x ) ) * b

    end if

  end if

  bessel_i0 = value

  return
end
subroutine bessel_i0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_I0_VALUES returns some values of the I0 Bessel function.
!
!  Discussion:
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    The modified Bessel function I0(Z) corresponds to N = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselI[0,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.1000000000000000D+01, &
    0.1010025027795146D+01, &
    0.1040401782229341D+01, &
    0.1092045364317340D+01, &
    0.1166514922869803D+01, &
    0.1266065877752008D+01, &
    0.1393725584134064D+01, &
    0.1553395099731217D+01, &
    0.1749980639738909D+01, &
    0.1989559356618051D+01, &
    0.2279585302336067D+01, &
    0.3289839144050123D+01, &
    0.4880792585865024D+01, &
    0.7378203432225480D+01, &
    0.1130192195213633D+02, &
    0.1748117185560928D+02, &
    0.2723987182360445D+02, &
    0.6723440697647798D+02, &
    0.4275641157218048D+03, &
    0.2815716628466254D+04 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.20D+00, &
    0.40D+00, &
    0.60D+00, &
    0.80D+00, &
    0.10D+01, &
    0.12D+01, &
    0.14D+01, &
    0.16D+01, &
    0.18D+01, &
    0.20D+01, &
    0.25D+01, &
    0.30D+01, &
    0.35D+01, &
    0.40D+01, &
    0.45D+01, &
    0.50D+01, &
    0.60D+01, &
    0.80D+01, &
    0.10D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function bessel_i1 ( arg )

!*****************************************************************************80
!
!! BESSEL_I1 evaluates the Bessel I function of order I.
!
!  Discussion:
!
!    The main computation evaluates slightly modified forms of
!    minimax approximations generated by Blair and Edwards.
!    This transportable program is patterned after the machine-dependent
!    FUNPACK packet NATSI1, but cannot match that version for efficiency
!    or accuracy.  This version uses rational functions that theoretically
!    approximate I-SUB-1(X) to at least 18 significant decimal digits.
!    The accuracy achieved depends on the arithmetic system, the compiler,
!    the intrinsic functions, and proper selection of the machine-dependent
!    constants.
!
!  Machine-dependent constants:
!
!    beta   = Radix for the floating-point system.
!    maxexp = Smallest power of beta that overflows.
!    XMAX =   Largest argument acceptable to BESI1;  Solution to
!             equation:
!               EXP(X) * (1-3/(8*X)) / SQRT(2*PI*X) = beta**maxexp
!
!
!    Approximate values for some important machines are:
!
!                            beta       maxexp    XMAX
!
!    CRAY-1        (S.P.)       2         8191    5682.810
!    Cyber 180/855
!      under NOS   (S.P.)       2         1070     745.894
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)       2          128      91.906
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)       2         1024     713.987
!    IBM 3033      (D.P.)      16           63     178.185
!    VAX           (S.P.)       2          127      91.209
!    VAX D-Format  (D.P.)       2          127      91.209
!    VAX G-Format  (D.P.)       2         1023     713.293
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2004
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Blair, Edwards,
!    Chalk River Report AECL-4928,
!    Atomic Energy of Canada, Limited,
!    October, 1974.
!
!  Parameters:
!
!    Input, real (kind = 8 ) ARG, the argument.
!
!    Output, real ( kind = 8 ) BESSEL_I1, the value of the Bessel
!    I1 function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) bessel_i1
  real ( kind = 8 ), parameter :: exp40 = 2.353852668370199854D+17
  real ( kind = 8 ), parameter :: forty = 40.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: one5 = 15.0D+00
  real ( kind = 8 ), dimension(15) :: p = (/ &
    -1.9705291802535139930D-19, &
    -6.5245515583151902910D-16, &
    -1.1928788903603238754D-12, &
    -1.4831904935994647675D-09, &
    -1.3466829827635152875D-06, &
    -9.1746443287817501309D-04, &
    -4.7207090827310162436D-01, &
    -1.8225946631657315931D+02, &
    -5.1894091982308017540D+04, &
    -1.0588550724769347106D+07, &
    -1.4828267606612366099D+09, &
    -1.3357437682275493024D+11, &
    -6.9876779648010090070D+12, &
    -1.7732037840791591320D+14, &
    -1.4577180278143463643D+15 /)
  real ( kind = 8 ) :: pbar = 3.98437500D-01
  real ( kind = 8 ), dimension(8) :: pp = (/ &
    -6.0437159056137600000D-02, &
     4.5748122901933459000D-01, &
    -4.2843766903304806403D-01, &
     9.7356000150886612134D-02, &
    -3.2457723974465568321D-03, &
    -3.6395264712121795296D-04, &
     1.6258661867440836395D-05, &
    -3.6347578404608223492D-07 /)
  real ( kind = 8 ), dimension(5) :: q = (/ &
    -4.0076864679904189921D+03, &
     7.4810580356655069138D+06, &
    -8.0059518998619764991D+09, &
     4.8544714258273622913D+12, &
    -1.3218168307321442305D+15 /)
  real ( kind = 8 ), dimension(6) :: qq = (/ &
    -3.8806586721556593450D+00, &
     3.2593714889036996297D+00, &
    -8.5017476463217924408D-01, &
     7.4212010813186530069D-02, &
    -2.2835624489492512649D-03, &
     3.7510433111922824643D-05 /)
  real ( kind = 8 ), parameter :: rec15 = 6.6666666666666666666D-02
  real ( kind = 8 ) sump
  real ( kind = 8 ) sumq
  real ( kind = 8 ), parameter :: two25 = 225.0D+00
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xmax = 713.987D+00
  real ( kind = 8 ) xx
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  x = abs ( arg )
!
!  ABS(ARG) < EPSILON ( ARG )
!
  if ( x < epsilon ( x ) ) then

    value = half * x
!
!  EPSILON ( ARG ) <= ABS(ARG) < 15.0
!
  else if ( x < one5 ) then

    xx = x * x
    sump = p(1)
    do j = 2, 15
      sump = sump * xx + p(j)
    end do

    xx = xx - two25

    sumq = (((( &
          xx + q(1) &
      ) * xx + q(2) &
      ) * xx + q(3) &
      ) * xx + q(4) &
      ) * xx + q(5)

    value = ( sump / sumq ) * x

  else if ( xmax < x ) then

    value = huge ( x )
!
!  15.0 <= ABS(ARG)
!
  else

    xx = one / x - rec15

    sump = ((((((    &
               pp(1) &
        * xx + pp(2) &
      ) * xx + pp(3) &
      ) * xx + pp(4) &
      ) * xx + pp(5) &
      ) * xx + pp(6) &
      ) * xx + pp(7) &
      ) * xx + pp(8)

    sumq = (((((     &
          xx + qq(1) &
      ) * xx + qq(2) &
      ) * xx + qq(3) &
      ) * xx + qq(4) &
      ) * xx + qq(5) &
      ) * xx + qq(6)

    value = sump / sumq

    if ( xmax - one5 < x ) then
      a = exp ( x - forty )
      b = exp40
    else
      a = exp ( x )
      b = one
    end if

    value = ( ( value * a + pbar * a ) / sqrt ( x ) ) * b

  end if

  if ( arg < zero ) then
    value = -value
  end if

  bessel_i1 = value

  return
end
subroutine bessel_i1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_I1_VALUES returns some values of the I1 Bessel function.
!
!  Discussion:
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselI[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.1005008340281251D+00, &
    0.2040267557335706D+00, &
    0.3137040256049221D+00, &
    0.4328648026206398D+00, &
    0.5651591039924850D+00, &
    0.7146779415526431D+00, &
    0.8860919814143274D+00, &
    0.1084810635129880D+01, &
    0.1317167230391899D+01, &
    0.1590636854637329D+01, &
    0.2516716245288698D+01, &
    0.3953370217402609D+01, &
    0.6205834922258365D+01, &
    0.9759465153704450D+01, &
    0.1538922275373592D+02, &
    0.2433564214245053D+02, &
    0.6134193677764024D+02, &
    0.3998731367825601D+03, &
    0.2670988303701255D+04 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.20D+00, &
    0.40D+00, &
    0.60D+00, &
    0.80D+00, &
    0.10D+01, &
    0.12D+01, &
    0.14D+01, &
    0.16D+01, &
    0.18D+01, &
    0.20D+01, &
    0.25D+01, &
    0.30D+01, &
    0.35D+01, &
    0.40D+01, &
    0.45D+01, &
    0.50D+01, &
    0.60D+01, &
    0.80D+01, &
    0.10D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function beta ( a, b )

!*****************************************************************************80
!
!! BETA returns the value of the Beta function.
!
!  Discussion:
!
!    The Beta function is defined as
!
!      BETA(A,B) = ( GAMMA ( A ) * GAMMA ( B ) ) / GAMMA ( A + B )
!                = Integral ( 0 <= T <= 1 ) T**(A-1) (1-T)**(B-1) dT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the function.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) BETA, the value of the function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma_log

  if ( a <= 0.0D+00 .or. b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA - Fatal error!'
    write ( *, '(a)' ) '  Both A and B must be greater than 0.'
    stop
  end if

  beta = exp ( gamma_log ( a ) + gamma_log ( b ) - gamma_log ( a + b ) )

  return
end
subroutine beta_binomial_cdf ( x, a, b, c, cdf )

!*****************************************************************************80
!
!! BETA_BINOMIAL_CDF evaluates the Beta Binomial CDF.
!
!  Discussion:
!
!    A simple summing approach is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Input, integer ( kind = 4 ) C, a parameter of the PDF.
!    0 <= C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y

  if ( x < 0 ) then

    cdf = 0.0D+00

  else if ( x < c ) then

    cdf = 0.0D+00
    do y = 0, x
      pdf = beta ( a + real ( y, kind = 8 ), &
        b + real ( c - y, kind = 8 ) ) / ( real ( c + 1, kind = 8 ) &
        * beta ( real ( y + 1 , kind = 8), &
        real ( c - y + 1, kind = 8 ) ) * beta ( a, b ) )
      cdf = cdf + pdf
    end do

  else if ( c <= x ) then

    cdf = 1.0D+00

  end if

  return
end
subroutine beta_binomial_cdf_inv ( cdf, a, b, c, x )

!*****************************************************************************80
!
!! BETA_BINOMIAL_CDF_INV inverts the Beta Binomial CDF.
!
!  Discussion:
!
!    A simple discrete approach is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, B, parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Input, integer ( kind = 4 ) C, a parameter of the PDF.
!    0 <= C.
!
!    Output, integer ( kind = 4 ) X, the smallest X whose cumulative density 
!    function is greater than or equal to CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cum
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_BINOMIAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  cum = 0.0D+00

  do y = 0, c

    pdf = beta ( a + real ( y, kind = 8 ), &
      b + real ( c - y, kind = 8 ) ) / ( real ( c + 1, kind = 8 ) &
      * beta ( real ( y + 1, kind = 8 ), &
      real ( c - y + 1, kind = 8 ) ) * beta ( a, b ) )

    cum = cum + pdf

    if ( cdf <= cum ) then
      x = y
      return
    end if

  end do

  x = c

  return
end
function beta_binomial_check ( a, b, c )

!*****************************************************************************80
!
!! BETA_BINOMIAL_CHECK checks the parameters of the Beta Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Input, integer ( kind = 4 ) C, a parameter of the PDF.
!    0 <= C.
!
!    Output, logical BETA_BINOMIAL_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical beta_binomial_check
  integer ( kind = 4 ) c

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_BINOMIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    beta_binomial_check = .false.
    return
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_BINOMIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    beta_binomial_check = .false.
    return
  end if

  if ( c < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_BINOMIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C < 0.'
    beta_binomial_check = .false.
    return
  end if

  beta_binomial_check = .true.

  return
end
subroutine beta_binomial_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! BETA_BINOMIAL_MEAN returns the mean of the Beta Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Input, integer ( kind = 4 ) C, a parameter of the PDF.
!    0 <= N.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) mean

  mean = real ( c, kind = 8 ) * a / ( a + b )

  return
end
subroutine beta_binomial_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! BETA_BINOMIAL_PDF evaluates the Beta Binomial PDF.
!
!  Discussion:
!
!    The PDF is defined as:
!
!      PDF(A,B,C;X) = Beta(A+X,B+C-X)
!        / ( (C+1) * Beta(X+1,C-X+1) * Beta(A,B) )  for 0 <= X <= C.
!
!    This PDF can be reformulated as:
!
!      The beta binomial probability density function for X successes
!      out of N trials is
!
!      PDF2(X)( N, MU, THETA ) =
!        C(N,X) * Product ( 0 <= R <= X - 1 ) ( MU + R * THETA )
!               * Product ( 0 <= R <= N - X - 1 ) ( 1 - MU + R * THETA )
!               / Product ( 0 <= R <= N - 1 )  ( 1 + R * THETA )
!
!      where
!
!        C(N,X) is the combinatorial coefficient;
!        MU is the expectation of the underlying Beta distribution;
!        THETA is a shape parameter.
!
!      A THETA value of 0 ( or A+B --> Infinity ) results in the binomial
!      distribution:
!
!        PDF2(X) ( N, MU, 0 ) = C(N,X) * MU**X * ( 1 - MU )**(N-X)
!
!    Given A, B, C for PDF, then the equivalent PDF2 has:
!
!      N     = C
!      MU    = A / ( A + B )
!      THETA = 1 / ( A + B )
!
!    Given N, MU, THETA for PDF2, the equivalent PDF has:
!
!      A = MU / THETA
!      B = ( 1 - MU ) / THETA
!      C = N
!
!    BETA_BINOMIAL_PDF(1,1,C;X) = UNIFORM_DISCRETE_PDF(0,C-1;X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Input, integer ( kind = 4 ) C, a parameter of the PDF.
!    0 <= C.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) c
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  if ( x < 0 ) then

    pdf = 0.0D+00

  else if ( x <= c ) then

    pdf = beta ( a + real ( x, kind = 8 ), b + real ( c - x, kind = 8 ) ) &
      / ( real ( c + 1, kind = 8 ) &
      * beta ( real ( x + 1, kind = 8 ), &
      real ( c - x + 1, kind = 8 ) ) * beta ( a, b ) )

  else if ( c < x ) then

    pdf = 0.0D+00

  end if

  return
end
subroutine beta_binomial_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! BETA_BINOMIAL_SAMPLE samples the Beta Binomial CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Input, integer ( kind = 4 ) C, a parameter of the PDF.
!    0 <= C.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x

  cdf = r8_uniform_01 ( seed )

  call beta_binomial_cdf_inv ( cdf, a, b, c, x )

  return
end
subroutine beta_binomial_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! BETA_BINOMIAL_VARIANCE returns the variance of the Beta Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Input, integer ( kind = 4 ) C, a parameter of the PDF.
!    0 <= C.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) variance

  variance = ( real ( c, kind = 8 ) * a * b ) &
    * ( a + b + real ( c, kind = 8 ) ) &
    / ( ( a + b )**2 * ( a + b + 1.0D+00 ) )

  return
end
subroutine beta_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! BETA_CDF evaluates the Beta CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta_inc
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then
    cdf = 0.0D+00
  else if ( x <= 1.0D+00 ) then
    cdf = beta_inc ( a, b, x )
  else
    cdf = 1.0D+00
  end if

  return
end
subroutine beta_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! BETA_CDF_INV inverts the Beta CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2001
!
!  Author:
!
!    Original FORTRAN77 version by Roger Abernathy, Robert Smith.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Roger Abernathy, Robert Smith,
!    Algorithm 724,
!    Program to Calculate F Percentiles,
!    ACM Transactions on Mathematical Software,
!    Volume 19, Number 4, December 1993, pages 481-483.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the argument of the CDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxk = 20

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bcoeff
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf_x
  real ( kind = 8 ) d(2:maxk,0:maxk-2)
  real ( kind = 8 ), parameter :: error = 0.0001D+00
  real ( kind = 8 ), parameter :: errapp = 0.01D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) loopct
  real ( kind = 8 ) pdf_x
  real ( kind = 8 ) q
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) sum2
  real ( kind = 8 ) t
  real ( kind = 8 ) tail
  real ( kind = 8 ) x
  real ( kind = 8 ) xold

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if
!
!  Estimate the solution.
!
  x = a / ( a + b )

  xold = 0.0D+00
  loopct = 2

  do while ( errapp <= abs ( ( x - xold ) / x ) .and. loopct /= 0 )

    xold = x
    loopct = loopct - 1
!
!  CDF_X = PROB { BETA(A,B) <= X }.
!  Q = ( CDF - CDF_X ) / PDF_X.
!
    call beta_cdf ( x, a, b, cdf_x )

    call beta_pdf ( x, a, b, pdf_x )

    q = ( cdf - cdf_x ) / pdf_x
!
!  D(N,K) = C(N,K) * Q**(N+K-1) / (N-1)!
!
    t = 1.0D+00 - x
    s1 = q * ( b - 1.0D+00 ) / t
    s2 = q * ( 1.0D+00 - a ) / x
    d(2,0) = s1 + s2
    tail = d(2,0) * q / 2.0D+00
    x = x + q + tail

    k = 3

    do while ( error < abs ( tail / x ) .and. k <= maxk )
!
!  Find D(2,K-2).
!
      s1 = q * ( real ( k, kind = 8 ) - 2.0D+00 ) * s1 / t
      s2 = q * ( 2.0D+00 - real ( k, kind = 8 ) ) * s2 / x
      d(2,k-2) = s1 + s2
!
!  Find D(3,K-3), D(4,K-4), D(5,K-5), ... , D(K-1,1).
!
      do i = 3, k-1
        sum2 = d(2,0) * d(i-1,k-i)
        bcoeff = 1.0D+00
        do j = 1, k-i
          bcoeff = ( bcoeff * real ( k - i - j + 1, kind = 8 ) ) &
            / real ( j, kind = 8 )
          sum2 = sum2 + bcoeff * d(2,j) * d(i-1,k-i-j)
        end do
        d(i,k-i) = sum2 + d(i-1,k-i+1) / real ( i - 1, kind = 8 )
      end do
!
!  Compute D(K,0) and use it to expand the series.
!
      d(k,0) = d(2,0) * d(k-1,0) + d(k-1,1) / real ( k - 1, kind = 8 )
      tail = d(k,0) * q / real ( k, kind = 8 )
      x = x + tail
!
!  Check for divergence.
!
      if ( x <= 0.0D+00 .or. 1.0D+00 <= x )  then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BETA_CDF_INV - Fatal error!'
        write ( *, '(a)' ) '  The series has diverged.'
        write ( *, '(a,g14.6)' ) '  X = ', x
        x = - 1.0D+00
        return
      end if

      k = k + 1

    end do

  end do

  return
end
subroutine beta_cdf_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! BETA_CDF_VALUES returns some values of the Beta CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = BetaDistribution [ a, b ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Karl Pearson,
!    Tables of the Incomplete Beta Function,
!    Cambridge University Press, 1968.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the parameters of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.10D+01, &
     0.10D+01, &
     0.10D+01, &
     0.10D+01, &
     0.10D+01, &
     0.10D+01, &
     0.10D+01, &
     0.10D+01, &
     0.20D+01, &
     0.30D+01, &
     0.40D+01, &
     0.50D+01 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     0.50D+00, &
     0.50D+00, &
     0.50D+00, &
     0.50D+00, &
     0.20D+01, &
     0.30D+01, &
     0.40D+01, &
     0.50D+01, &
     0.20D+01, &
     0.20D+01, &
     0.20D+01, &
     0.20D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5131670194948620D-01, &
    0.1055728090000841D+00, &
    0.1633399734659245D+00, &
    0.2254033307585166D+00, &
    0.3600000000000000D+00, &
    0.4880000000000000D+00, &
    0.5904000000000000D+00, &
    0.6723200000000000D+00, &
    0.2160000000000000D+00, &
    0.8370000000000000D-01, &
    0.3078000000000000D-01, &
    0.1093500000000000D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function beta_check ( a, b )

!*****************************************************************************80
!
!! BETA_CHECK checks the parameters of the Beta PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, logical BETA_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical beta_check

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    beta_check = .false.
    return
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    beta_check = .false.
    return
  end if

  beta_check = .true.

  return
end
function beta_inc ( a, b, x )

!*****************************************************************************80
!
!! BETA_INC returns the value of the incomplete Beta function.
!
!  Discussion:
!
!    This calculation requires an iteration.  In some cases, the iteration
!    may not converge rapidly, or may become inaccurate.
!
!    The formula is:
!
!      BETA_INC(A,B,X)
!
!        =   Integral ( 0 <= T <= X ) T**(A-1) (1-T)**(B-1) dT
!          / Integral ( 0 <= T <= 1 ) T**(A-1) (1-T)**(B-1) dT
!
!        =   Integral ( 0 <= T <= X ) T**(A-1) (1-T)**(B-1) dT
!          / BETA(A,B)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    KL Majumder, GP Bhattacharjee,
!    Algorithm AS63,
!    Applied Statistics,
!    1973, volume 22, number 3.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the function.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!    Normally, 0.0D+00 <= X <= 1.0.
!
!    Output, real ( kind = 8 ) BETA_INC, the value of the function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) beta_inc
  real ( kind = 8 ) cx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 1000
  logical indx
  integer ( kind = 4 ) ns
  real ( kind = 8 ) pp
  real ( kind = 8 ) psq
  real ( kind = 8 ) qq
  real ( kind = 8 ) rx
  real ( kind = 8 ) temp
  real ( kind = 8 ) term
  real ( kind = 8 ), parameter :: tol = 1.0D-07
  real ( kind = 8 ) x
  real ( kind = 8 ) xx

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_INC - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    stop
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_INC - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    stop
  end if

  if ( x <= 0.0D+00 ) then
    beta_inc = 0.0D+00
    return
  else if ( 1.0D+00 <= x ) then
    beta_inc = 1.0D+00
    return
  end if
!
!  Change tail if necessary and determine S.
!
  psq = a + b

  if ( a < ( a + b ) * x ) then
    xx = 1.0D+00 - x
    cx = x
    pp = b
    qq = a
    indx = .true.
  else
    xx = x
    cx = 1.0D+00 - x
    pp = a
    qq = b
    indx = .false.
  end if

  term = 1.0D+00
  i = 1
  beta_inc = 1.0D+00

  ns = int ( qq + cx * ( a + b ) )
!
!  Use Soper's reduction formulas.
!
  rx = xx / cx

  temp = qq - real ( i, kind = 8 )
  if ( ns == 0 ) then
    rx = xx
  end if

  it = 0

  do

    it = it + 1

    if ( it_max < it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BETA_INC - Fatal error!'
      write ( *, '(a)' ) '  Maximum number of iterations exceeded!'
      write ( *, '(a,i8)' ) '  IT_MAX = ', it_max
      stop
    end if

    term = term * temp * rx / ( pp + real ( i, kind = 8 ) )
    beta_inc = beta_inc + term
    temp = abs ( term )

    if ( temp <= tol .and. temp <= tol * beta_inc ) then
      exit
    end if

    i = i + 1
    ns = ns - 1

    if ( 0 <= ns ) then
      temp = qq - real ( i, kind = 8 )
      if ( ns == 0 ) then
        rx = xx
      end if
    else
      temp = psq
      psq = psq + 1.0D+00
    end if

  end do
!
!  Finish calculation.
!
  beta_inc = beta_inc * exp ( pp * log ( xx ) &
    + ( qq - 1.0D+00 ) * log ( cx ) ) / ( beta ( a, b ) * pp )

  if ( indx ) then
    beta_inc = 1.0D+00 - beta_inc
  end if

  return
end
subroutine beta_inc_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! BETA_INC_VALUES returns some values of the incomplete Beta function.
!
!  Discussion:
!
!    The incomplete Beta function may be written
!
!      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
!                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
!
!    Thus,
!
!      BETA_INC(A,B,0.0) = 0.0D+00
!      BETA_INC(A,B,1.0) = 1.0
!
!    In Mathematica, the function can be evaluated by:
!
!      BETA[X,A,B] / BETA[A,B]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Karl Pearson,
!    Tables of the Incomplete Beta Function,
!    Cambridge University Press, 1968.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the parameters of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.5D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    30.0D+00, &
    30.0D+00, &
    40.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.0D+00, &
     0.5D+00, &
     5.0D+00, &
     5.0D+00, &
    10.0D+00, &
     5.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6376856085851985D-01, &
    0.2048327646991335D+00, &
    0.1000000000000000D+01, &
    0.0000000000000000D+00, &
    0.5012562893380045D-02, &
    0.5131670194948620D-01, &
    0.2928932188134525D+00, &
    0.5000000000000000D+00, &
    0.2800000000000000D-01, &
    0.1040000000000000D+00, &
    0.2160000000000000D+00, &
    0.3520000000000000D+00, &
    0.5000000000000000D+00, &
    0.6480000000000000D+00, &
    0.7840000000000000D+00, &
    0.8960000000000000D+00, &
    0.9720000000000000D+00, &
    0.4361908850559777D+00, &
    0.1516409096347099D+00, &
    0.8978271484375000D-01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.4598773297575791D+00, &
    0.2146816102371739D+00, &
    0.9507364826957875D+00, &
    0.5000000000000000D+00, &
    0.8979413687105918D+00, &
    0.2241297491808366D+00, &
    0.7586405487192086D+00, &
    0.7001783247477069D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.01D+00, &
    0.10D+00, &
    1.00D+00, &
    0.00D+00, &
    0.01D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    0.50D+00, &
    0.90D+00, &
    0.50D+00, &
    1.00D+00, &
    0.50D+00, &
    0.80D+00, &
    0.60D+00, &
    0.80D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.70D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine beta_mean ( a, b, mean )

!*****************************************************************************80
!
!! BETA_MEAN returns the mean of the Beta PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a / ( a + b )

  return
end
subroutine beta_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! BETA_PDF evaluates the Beta PDF.
!
!  Discussion:
!
!    The formula for the PDF is:
!
!      PDF(A,B;X) = X**(A-1) * (1-X)**(B-1) / BETA(A,B).
!
!    A = B = 1 yields the Uniform distribution on [0,1].
!    A = B = 1/2 yields the Arcsin distribution.
!        B = 1 yields the power function distribution.
!    A = B -> Infinity tends to the Normal distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
    pdf = 0.0D+00
  else
    pdf = x**( a - 1.0D+00 ) * ( 1.0D+00 - x )**( b - 1.0D+00 ) / beta ( a, b )
  end if

  return
end
subroutine beta_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! BETA_SAMPLE samples the Beta PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Kennedy, James Gentle,
!    Algorithm BN,
!    Statistical Computing,
!    Dekker, 1980.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) mu
  integer ( kind = 4 ) seed
  real ( kind = 8 ) stdev
  real ( kind = 8 ) test
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  mu = ( a - 1.0D+00 ) / ( a + b - 2.0D+00 )
  stdev = 0.5D+00 / sqrt ( a + b - 2.0D+00 )

  do

    call normal_01_sample ( seed, y )

    x = mu + stdev * y

    if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
      cycle
    end if

    u = r8_uniform_01 ( seed )

    test =     ( a - 1.0D+00 )     * log (         x   / ( a - 1.0D+00 ) ) &
             + ( b - 1.0D+00 )     * log ( ( 1.0D+00 - x ) / ( b - 1.0D+00 ) ) &
             + ( a + b - 2.0D+00 ) * log ( a + b - 2.0D+00 ) + 0.5D+00 * y * y

    if ( log ( u ) <= test ) then
      exit
    end if

  end do

  return
end
subroutine beta_variance ( a, b, variance )

!*****************************************************************************80
!
!! BETA_VARIANCE returns the variance of the Beta PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = ( a * b ) / ( ( a + b )**2 * ( 1.0D+00 + a + b ) )

  return
end
subroutine binomial_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! BINOMIAL_CDF evaluates the Binomial CDF.
!
!  Discussion:
!
!    CDF(X)(A,B) is the probability of at most X successes in A trials,
!    given that the probability of success on a single trial is B.
!
!    A sequence of trials with fixed probability of success on
!    any trial is known as a sequence of Bernoulli trials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the desired number of successes.
!    0 <= X <= A.
!
!    Input, integer ( kind = 4 ) A, the number of trials.
!    1 <= A.
!
!    Input, real ( kind = 8 ) B, the probability of success on one trial.
!    0.0D+00 <= B <= 1.0.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) cnk
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) j
  real ( kind = 8 ) pr
  integer ( kind = 4 ) x

  if ( x < 0 ) then

    cdf = 0.0D+00

  else if ( a <= x ) then

    cdf = 1.0D+00

  else if ( b == 0.0D+00 ) then

    cdf = 1.0D+00

  else if ( b == 1.0D+00 ) then

    cdf = 0.0D+00

  else

    cdf = 0.0D+00

    do j = 0, x

      call binomial_coef ( a, j, cnk )

      pr = real ( cnk, kind = 8 ) * b**j * ( 1.0D+00 - b )**( a - j )

      cdf = cdf + pr

    end do

  end if

  return
end
subroutine binomial_cdf_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
!
!  Discussion:
!
!    CDF(X)(A,B) is the probability of at most X successes in A trials,
!    given that the probability of success on a single trial is B.
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`DiscreteDistributions`]
!      dist = BinomialDistribution [ n, p ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition, CRC Press, 1996, pages 651-652.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) A, a parameter of the function.
!
!    Output, real ( kind = 8 ) B, a parameter of the function.
!
!    Output, integer ( kind = 4 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  integer ( kind = 4 ) a
  integer ( kind = 4 ), save, dimension ( n_max ) :: a_vec = (/ &
     2,  2,  2,  2, &
     2,  4,  4,  4, &
     4, 10, 10, 10, &
    10, 10, 10, 10, &
    10 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
    0.05D+00, &
    0.05D+00, &
    0.05D+00, &
    0.50D+00, &
    0.50D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.05D+00, &
    0.10D+00, &
    0.15D+00, &
    0.20D+00, &
    0.25D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.9025000000000000D+00, &
    0.9975000000000000D+00, &
    0.1000000000000000D+01, &
    0.2500000000000000D+00, &
    0.7500000000000000D+00, &
    0.3164062500000000D+00, &
    0.7382812500000000D+00, &
    0.9492187500000000D+00, &
    0.9960937500000000D+00, &
    0.9999363101685547D+00, &
    0.9983650626000000D+00, &
    0.9901259090013672D+00, &
    0.9672065024000000D+00, &
    0.9218730926513672D+00, &
    0.8497316674000000D+00, &
    0.6331032576000000D+00, &
    0.3769531250000000D+00 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) x
  integer ( kind = 4 ), save, dimension ( n_max ) :: x_vec = (/ &
     0, 1, 2, 0, &
     1, 0, 1, 2, &
     3, 4, 4, 4, &
     4, 4, 4, 4, &
     4 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0
    b = 0.0D+00
    x = 0
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine binomial_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! BINOMIAL_CDF_INV inverts the Binomial CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, integer ( kind = 4 ) A, the number of trials.
!    1 <= A.
!
!    Input, real ( kind = 8 ) B, the probability of success on one trial.
!    0.0D+00 <= B <= 1.0.
!
!    Output, integer ( kind = 4 ) X, the corresponding argument.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BINOMIAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  cdf2 = 0.0D+00

  do x2 = 0, a

    call binomial_pdf ( x2, a, b, pdf )

    cdf2 = cdf2 + pdf

    if ( cdf <= cdf2 ) then
      x = x2
      return
    end if

  end do

  return
end
function binomial_check ( a, b )

!*****************************************************************************80
!
!! BINOMIAL_CHECK checks the parameter of the Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of trials.
!    1 <= A.
!
!    Input, real ( kind = 8 ) B, the probability of success on one trial.
!    0.0D+00 <= B <= 1.0.
!
!    Output, logical BINOMIAL_CHECK, is true if the parameters are legal.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  logical binomial_check

  if ( a < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BINOMIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A < 1.'
    binomial_check = .false.
    return
  end if

  if ( b < 0.0D+00 .or. 1.0D+00 < b ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BINOMIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B < 0 or 1 < B.'
    binomial_check = .false.
    return
  end if

  binomial_check = .true.

  return
end
subroutine binomial_coef ( n, k, cnk )

!*****************************************************************************80
!
!! BINOMIAL_COEF computes the Binomial coefficient C(N,K).
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula is:
!
!      CNK = C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Combinatorial of M Things Taken N at a Time,
!    ACM Algorithm 160,
!    Communications of the ACM,
!    April, 1963.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, integer ( kind = 4 ) CNK, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) cnk
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    cnk = 0

  else if ( mn == 0 ) then

    cnk = 1

  else

    mx = max ( k, n - k )
    cnk = mx + 1

    do i = 2, mn
      cnk = ( cnk * ( mx + i ) ) / i
    end do

  end if

  return
end
subroutine binomial_coef_log ( n, k, cnk_log )

!*****************************************************************************80
!
!! BINOMIAL_COEF_LOG computes the logarithm of the Binomial coefficient.
!
!  Discussion:
!
!    The formula is:
!
!      CNK_LOG = LOG ( C(N,K) ) = LOG ( N! / ( K! * (N-K)! ) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, real ( kind = 8 ) CNK_LOG, the logarithm of C(N,K).
!
  implicit none

  real ( kind = 8 ) cnk_log
  real ( kind = 8 ) factorial_log
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  cnk_log = factorial_log ( n ) - factorial_log ( k ) - factorial_log ( n - k )

  return
end
subroutine binomial_mean ( a, b, mean )

!*****************************************************************************80
!
!! BINOMIAL_MEAN returns the mean of the Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of trials.
!    1 <= A.
!
!    Input, real ( kind = 8 ) B, the probability of success on one trial.
!    0.0D+00 <= B <= 1.0.
!
!    Output, real ( kind = 8 ) MEAN, the expected value of the number of
!    successes in A trials.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = real ( a, kind = 8 ) * b

  return
end
subroutine binomial_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! BINOMIAL_PDF evaluates the Binomial PDF.
!
!  Discussion:
!
!    PDF(A,B;X) is the probability of exactly X successes in A trials,
!    given that the probability of success on a single trial is B.
!
!    The formula is:
!
!      PDF(A,B;X) = C(N,X) * B**X * ( 1.0D+00 - B )**( A - X )
!
!    Binomial_PDF(1,B;X) = Bernoulli_PDF(B;X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the desired number of successes.
!    0 <= X <= A.
!
!    Input, integer ( kind = 4 ) A, the number of trials.
!    1 <= A.
!
!    Input, real ( kind = 8 ) B, the probability of success on one trial.
!    0.0D+00 <= B <= 1.0.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) cnk
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  if ( a < 1 ) then

    pdf = 0.0D+00

  else if ( x < 0 .or. a < x ) then

    pdf = 0.0D+00

  else if ( b == 0.0D+00 ) then

    if ( x == 0 ) then
      pdf = 1.0D+00
    else
      pdf = 0.0D+00
    end if

  else if ( b == 1.0D+00 ) then

    if ( x == a ) then
      pdf = 1.0D+00
    else
      pdf = 0.0D+00
    end if

  else

    call binomial_coef ( a, x, cnk )

    pdf = real ( cnk, kind = 8 ) * b**x * ( 1.0D+00 - b )**( a - x )

  end if

  return
end
subroutine binomial_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! BINOMIAL_SAMPLE samples the Binomial PDF.
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
!  Reference:
!
!    William Kennedy, James Gentle,
!    Algorithm BU,
!    Statistical Computing,
!    Dekker, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of trials.
!    1 <= A.
!
!    Input, real ( kind = 8 ) B, the probability of success on one trial.
!    0.0D+00 <= B <= 1.0.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u
  integer ( kind = 4 ) x

  x = 0

  do i = 1, a

    u = r8_uniform_01 ( seed )

    if ( u <= b ) then
      x = x + 1
    end if

  end do

  return
end
subroutine binomial_variance ( a, b, variance )

!*****************************************************************************80
!
!! BINOMIAL_VARIANCE returns the variance of the Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of trials.
!    1 <= A.
!
!    Input, real ( kind = 8 ) B, the probability of success on one trial.
!    0.0D+00 <= B <= 1.0.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = real ( a, kind = 8 ) * b * ( 1.0D+00 - b )

  return
end
subroutine bradford_cdf ( x, a, b, c, cdf )

!*****************************************************************************80
!
!! BRADFORD_CDF evaluates the Bradford CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x <= a ) then
    cdf = 0.0D+00
  else if ( x <= b ) then
    cdf = log ( 1.0D+00 + c * ( x - a ) / ( b - a ) ) / log ( c + 1.0D+00 )
  else if ( b < x ) then
    cdf = 1.0D+00
  end if

  return
end
subroutine bradford_cdf_inv ( cdf, a, b, c, x )

!*****************************************************************************80
!
!! BRADFORD_CDF_INV inverts the Bradford CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BRADFORD_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf <= 0.0D+00 ) then
    x = a
  else if ( cdf < 1.0D+00 ) then
    x = a + ( b - a ) * ( ( c + 1.0D+00 )**cdf - 1.0D+00 ) / c
  else if ( 1.0D+00 <= cdf ) then
    x = b
  end if

  return
end
function bradford_check ( a, b, c )

!*****************************************************************************80
!
!! BRADFORD_CHECK checks the parameters of the Bradford PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A < B,
!    0.0D+00 < C.
!
!    Output, logical BRADFORD_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical bradford_check
  real ( kind = 8 ) c

  if ( b <= a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BRADFORD_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= A.'
    bradford_check = .false.
    return
  end if

  if ( c <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BRADFORD_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C <= 0.'
    bradford_check = .false.
    return
  end if

  bradford_check = .true.

  return
end
subroutine bradford_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! BRADFORD_MEAN returns the mean of the Bradford PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) mean

  mean = ( c * ( b - a ) + log ( c + 1.0D+00 ) * ( a * ( c + 1.0D+00 ) - b ) ) &
    / ( c * log ( c + 1.0D+00 ) )

  return
end
subroutine bradford_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! BRADFORD_PDF evaluates the Bradford PDF.
!
!  Discussion:
!
!    The formula is:
!
!      PDF(A,B,C;X) =
!        C / ( ( C * ( X - A ) + B - A ) * log ( C + 1 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A <= X
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x <= a ) then
    pdf = 0.0D+00
  else if ( x <= b ) then
    pdf = c / ( ( c * ( x - a ) + b - a ) * log ( c + 1.0D+00 ) )
  else if ( b < x ) then
    pdf = 0.0D+00
  end if

  return
end
subroutine bradford_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! BRADFORD_SAMPLE samples the Bradford PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A < B,
!    0.0D+00 < C.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  x = a + ( b - a ) * ( ( c + 1.0D+00 )**cdf - 1.0D+00 ) / c

  return
end
subroutine bradford_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! BRADFORD_VARIANCE returns the variance of the Bradford PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) variance

  variance = ( b - a )**2 * &
    ( c * ( log ( c + 1.0D+00 ) - 2.0D+00 ) + 2.0D+00 * log ( c + 1.0D+00 ) ) &
    / ( 2.0D+00 * c * ( log ( c + 1.0D+00 ) )**2 )

  return
end
subroutine buffon_laplace_pdf ( a, b, l, pdf )

!*****************************************************************************80
!
!! BUFFON_LAPLACE_PDF evaluates the Buffon-Laplace PDF.
!
!  Discussion:
!
!    In the Buffon-Laplace needle experiment, we suppose that the plane has been
!    tiled into a grid of rectangles of width A and height B, and that a
!    needle of length L is dropped "at random" onto this grid.
!
!    We may assume that one end, the "eye" of the needle falls at the point
!    (X1,Y1), taken uniformly at random in the cell [0,A]x[0,B].
!
!    ANGLE, the angle that the needle makes is taken to be uniformly random.
!    The point of the needle, (X2,Y2), therefore lies at
!
!      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
!
!    The needle will have crossed at least one grid line if any of the
!    following are true:
!
!      X2 <= 0, A <= X2, Y2 <= 0, B <= Y2.
!
!    If L is larger than sqrt ( A*A + B*B ), then the needle will
!    cross every time, and the computation is uninteresting.  However, if
!    L is smaller than this limit, then the probability of a crossing on
!    a single trial is
!
!      P(L,A,B) = ( 2 * L * ( A + B ) - L * L ) / ( PI * A * B )
!
!    and therefore, a record of the number of hits for a given number of
!    trials can be used as a very roundabout way of estimating PI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Sudarshan Raghunathan,
!    Making a Supercomputer Do What You Want: High Level Tools for
!    Parallel Programming,
!    Computing in Science and Engineering,
!    Volume 8, Number 5, September/October 2006, pages 70-80.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the horizontal and vertical dimensions
!    of each cell of the grid.  0 <= A, 0 <= B.
!
!    Input, real ( kind = 8 ) L, the length of the needle.
!    0 <= L <= min ( A, B ).
!
!    Output, real ( kind = 8 ) PDF, the Buffon-Laplace PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) l
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793238462643D+00

  if ( a < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BUFFON_LAPLACE_PDF - Fatal error!'
    write ( *, '(a)' ) '  Input A < 0.'
    stop
  else if ( a == 0.0D+00 ) then
    pdf = 1.0D+00
    return
  end if

  if ( b < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BUFFON_LAPLACE_PDF - Fatal error!'
    write ( *, '(a)' ) '  Input B < 0.'
    stop
  else if ( b == 0.0D+00 ) then
    pdf = 1.0D+00
    return
  end if

  if ( l < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BUFFON_LAPLACE_PDF - Fatal error!'
    write ( *, '(a)' ) '  Input L < 0.'
    stop
  else if ( l == 0.0D+00 ) then
    pdf = 0.0D+00
    return
  else if ( min ( a, b ) < l ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BUFFON_LAPLACE_PDF - Fatal error!'
    write ( *, '(a)' ) '  min ( A, B ) < L.'
    stop
  end if

  pdf = l * ( 2.0D+00 * ( a + b ) - l ) / ( pi * a * b )

  return
end
function buffon_laplace_simulate ( a, b, l, trial_num )

!*****************************************************************************80
!
!! BUFFON_LAPLACE_SIMULATE simulates a Buffon-Laplace needle experiment.
!
!  Discussion:
!
!    In the Buffon-Laplace needle experiment, we suppose that the plane has 
!    been tiled into a grid of rectangles of width A and height B, and that a
!    needle of length L is dropped "at random" onto this grid.
!
!    We may assume that one end, the "eye" of the needle falls at the point
!    (X1,Y1), taken uniformly at random in the cell [0,A]x[0,B].
!
!    ANGLE, the angle that the needle makes is taken to be uniformly random.
!    The point of the needle, (X2,Y2), therefore lies at
!
!      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
!
!    The needle will have crossed at least one grid line if any of the
!    following are true:
!
!      X2 <= 0, A <= X2, Y2 <= 0, B <= Y2.
!
!    This routine simulates the tossing of the needle, and returns the number
!    of times that the needle crossed at least one grid line.
!
!    If L is larger than sqrt ( A*A + B*B ), then the needle will
!    cross every time, and the computation is uninteresting.  However, if
!    L is smaller than this limit, then the probability of a crossing on
!    a single trial is
!
!      P(L,A,B) = ( 2 * L * ( A + B ) - L * L ) / ( PI * A * B )
!
!    and therefore, a record of the number of hits for a given number of
!    trials can be used as a very roundabout way of estimating PI.
!    (Particularly roundabout, since we actually will use a good value of
!    PI in order to pick the random angles!)
!
!    Note that this routine will try to generate 5 * TRIAL_NUM random
!    double precision values at one time, using automatic arrays.
!    When I tried this with TRIAL_NUM = 1,000,000, the program failed,
!    because of internal system limits on such arrays.
!
!    Such a problem could be avoided by using a DO loop running through
!    each trial individually, but this tend to run much more slowly than
!    necessary.
!
!    Since this routine invokes the FORTRAN90 random number generator,
!    the user should initialize the random number generator, particularly
!    if it is desired to control whether the sequence is to be varied
!    or repeated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Sudarshan Raghunathan,
!    Making a Supercomputer Do What You Want: High Level Tools for
!    Parallel Programming,
!    Computing in Science and Engineering,
!    Volume 8, Number 5, September/October 2006, pages 70-80.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the horizontal and vertical dimensions
!    of each cell of the grid.  0 <= A, 0 <= B.
!
!    Input, real ( kind = 8 ) L, the length of the needle.
!    0 <= L <= min ( A, B ).
!
!    Input, integer ( kind = 4 ) TRIAL_NUM, the number of times the needle is
!    to be dropped onto the grid.
!
!    Output, integer ( kind = 4 ) BUFFON_LAPLACE_SIMULATE, the number of times
!    the needle crossed at least one line of the grid of cells.
!
!  Local Parameters:
!
!    Local, integer BATCH_SIZE, specifies the number of trials to be done
!    in a single batch.  Setting BATCH_SIZE to 1 will be very slow.
!    Replacing it by TRIAL_NUM would be fine except that your system
!    may have a limit on the size of automatic arrays.  We have set a default
!    value of 10,000 here which should be large enough to be efficient
!    but small enough not to annoy the system.
!
  implicit none

  integer ( kind = 4 ), parameter :: batch_size = 10000
  integer ( kind = 4 ) trial_num

  real ( kind = 8 ) a
  real ( kind = 8 ) angle(batch_size)
  real ( kind = 8 ) b
  integer ( kind = 4 ) batch
  integer ( kind = 4 ) buffon_laplace_simulate
  integer ( kind = 4 ) hits
  real ( kind = 8 ) l
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793238462643D+00
  real ( kind = 8 ) x1(batch_size)
  real ( kind = 8 ) x2(batch_size)
  real ( kind = 8 ) y1(batch_size)
  real ( kind = 8 ) y2(batch_size)

  hits = 0

  do batch = 1, trial_num, batch_size

    n = min ( batch_size, trial_num + 1 - batch )
!
!  Randomly choose the location of the eye of the needle in [0,0]x[A,B],
!  and the angle the needle makes.
!
    call random_number ( harvest = x1(1:n) )
    call random_number ( harvest = y1(1:n) )
    call random_number ( harvest = angle(1:n) )

    x1(1:n) = a * x1(1:n)
    y1(1:n) = b * y1(1:n)
    angle(1:n) = 2.0D+00 * pi * angle(1:n)
!
!  Compute the location of the point of the needle.
!
    x2(1:n) = x1(1:n) + l * cos ( angle(1:n) )
    y2(1:n) = y1(1:n) + l * sin ( angle(1:n) )
!
!  Count the end locations that lie outside the cell.
!
    hits = hits + count (      x2(1:n) <= 0.0 .or. &
                          a <= x2(1:n)        .or. &
                               y2(1:n) <= 0.0 .or. &
                          b <= y2(1:n) )

  end do

  buffon_laplace_simulate = hits

  return
end
subroutine buffon_pdf ( a, l, pdf )

!*****************************************************************************80
!
!! BUFFON_PDF evaluates the Buffon PDF.
!
!  Discussion:
!
!    In the Buffon needle experiment, we suppose that the plane has been
!    ruled by vertical lines with a spacing of A units, and that a
!    needle of length L is dropped "at random" onto this grid.
!
!    Because of the various symmetries, we may assume that this eye of
!    this needle lands in the first infinite strip, and we may further
!    assume that its Y coordinate is 0.  Thus, we have
!    the eye as (X1,Y1) with 0 <= X1 <= A and Y1 = 0.
!
!    ANGLE, the angle that the needle makes is taken to be uniformly random.
!    The point of the needle, (X2,Y2), therefore lies at
!
!      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
!
!    The needle will have crossed at least one grid line if any of the
!    following are true:
!
!      X2 <= 0, A <= X2.
!
!    The probability of a crossing on a single trial is
!
!      P(A,L) = ( 2 * L ) / ( PI * A )
!
!    and therefore, a record of the number of hits for a given number of
!    trials can be used as a very roundabout way of estimating PI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the horizontal spacing between the
!    vertical grid lines.  0 <= A.
!
!    Input, real ( kind = 8 ) L, the length of the needle.
!
!    Output, real ( kind = 8 ) PDF, the Buffon PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) l
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793238462643D+00

  if ( a < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BUFFON_PDF - Fatal error!'
    write ( *, '(a)' ) '  Input A < 0.'
    stop
  else if ( a == 0.0D+00 ) then
    pdf = 1.0D+00
    return
  end if

  if ( l < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BUFFON_PDF - Fatal error!'
    write ( *, '(a)' ) '  Input L < 0.'
    stop
  else if ( l == 0.0D+00 ) then
    pdf = 0.0D+00
    return
  end if

  pdf = ( 2.0D+00 * l ) / ( pi * a )

  return
end
function buffon_simulate ( a, l, trial_num )

!*****************************************************************************80
!
!! BUFFON_SIMULATE simulates a Buffon needle experiment.
!
!  Discussion:
!
!    In the Buffon needle experiment, we suppose that the plane has been
!    ruled by vertical lines with a spacing of A units, and that a
!    needle of length L is dropped "at random" onto this grid.
!
!    Because of the various symmetries, we may assume that this eye of
!    this needle lands in the first infinite strip, and we may further
!    assume that its Y coordinate is 0.  Thus, we have
!    the eye as (X1,Y1) with 0 <= X1 <= A and Y1 = 0.
!
!    ANGLE, the angle that the needle makes is taken to be uniformly random.
!    The point of the needle, (X2,Y2), therefore lies at
!
!      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
!
!    The needle will have crossed at least one grid line if any of the
!    following are true:
!
!      X2 <= 0, A <= X2.
!
!    The probability of a crossing on a single trial is
!
!      P(A,L) = ( 2 * L ) / ( PI * A )
!
!    and therefore, a record of the number of hits for a given number of
!    trials can be used as a very roundabout way of estimating PI.
!
!    Note that this routine will try to generate 4 * TRIAL_NUM random
!    double precision values at one time, using automatic arrays.
!    When I tried this with TRIAL_NUM = 1,000,000, the program failed,
!    because of internal system limits on such arrays.
!
!    Such a problem could be avoided by using a DO loop running through
!    each trial individually, but this tend to run much more slowly than
!    necessary.
!
!    Since this routine invokes the FORTRAN90 random number generator,
!    the user should initialize the random number generator, particularly
!    if it is desired to control whether the sequence is to be varied
!    or repeated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the horizontal spacing between the
!    vertical grid lines.  0 <= A.
!
!    Input, real ( kind = 8 ) L, the length of the needle.
!
!    Input, integer ( kind = 4 ) TRIAL_NUM, the number of times the needle is
!    to be dropped onto the grid.
!
!    Output, integer ( kind = 4 ) BUFFON_SIMULATE, the number of times the needle
!    crossed at least one line of the grid of cells.
!
!  Local Parameters:
!
!    Local, integer BATCH_SIZE, specifies the number of trials to be done
!    in a single batch.  Setting BATCH_SIZE to 1 will be very slow.
!    Replacing it by TRIAL_NUM would be fine except that your system
!    may have a limit on the size of automatic arrays.  We have set a default
!    value of 10,000 here which should be large enough to be efficient
!    but small enough not to annoy the system.
!
  implicit none

  integer ( kind = 4 ), parameter :: batch_size = 10000
  integer ( kind = 4 ) trial_num

  real ( kind = 8 ) a
  real ( kind = 8 ) angle(batch_size)
  integer ( kind = 4 ) batch
  integer ( kind = 4 ) buffon_simulate
  integer ( kind = 4 ) hits
  real ( kind = 8 ) l
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793238462643D+00
  real ( kind = 8 ) x1(batch_size)
  real ( kind = 8 ) x2(batch_size)

  hits = 0

  do batch = 1, trial_num, batch_size

    n = min ( batch_size, trial_num + 1 - batch )
!
!  Randomly choose the location (X1,Y1) of the eye of the needle
!  in [0,0]x[A,0], and the angle the needle makes.
!
    call random_number ( harvest = x1(1:n) )
    call random_number ( harvest = angle(1:n) )

    x1(1:n) = a * x1(1:n)
    angle(1:n) = 2.0D+00 * pi * angle(1:n)
!
!  Compute the location of the point of the needle.
!  We only need to know the value of X2, not Y2!
!
    x2(1:n) = x1(1:n) + l * cos ( angle(1:n) )
!
!  Count the end locations that lie outside the cell.
!
    hits = hits + count (      x2(1:n) <= 0.0 .or. &
                          a <= x2(1:n) )

  end do

  buffon_simulate = hits

  return
end
subroutine burr_cdf ( x, a, b, c, d, cdf )

!*****************************************************************************80
!
!! BURR_CDF evaluates the Burr CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) d
  real ( kind = 8 ) x

  if ( x <= a ) then

    cdf = 0.0D+00

  else

    cdf = 1.0D+00 / ( 1.0D+00 + ( b / ( x - a ) )**c  )**d

  end if

  return
end
subroutine burr_cdf_inv ( cdf, a, b, c, d, x )

!*****************************************************************************80
!
!! BURR_CDF_INV inverts the Burr CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) d
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BURR_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a + b / ( ( 1.0D+00 / cdf )**(1.0D+00 / d ) - 1.0D+00 )**( 1.0D+00 / c )

  return
end
function burr_check ( a, b, c, d )

!*****************************************************************************80
!
!! BURR_CHECK checks the parameters of the Burr CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, logical BURR_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical burr_check
  real ( kind = 8 ) c
  real ( kind = 8 ) d

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BURR_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    burr_check = .false.
    return
  end if

  if ( c <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BURR_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C <= 0.'
    burr_check = .false.
    return
  end if

  burr_check = .true.

  return
end
subroutine burr_mean ( a, b, c, d, mean )

!*****************************************************************************80
!
!! BURR_MEAN returns the mean of the Burr PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) mean
  real ( kind = 8 ) r8_gamma

  mean = a + b * r8_gamma ( 1.0D+00 - 1.0D+00 / c ) &
    * r8_gamma ( d + 1.0D+00 / c ) / r8_gamma ( d )

  return
end
subroutine burr_pdf ( x, a, b, c, d, pdf )

!*****************************************************************************80
!
!! BURR_PDF evaluates the Burr PDF.
!
!  Discussion:
!
!    The formula is:
!
!      PDF(A,B,C,D;X) = ( C * D / B ) * ( ( X - A ) / B )**( - C - 1 )
!        * ( 1 + ( ( X - A ) / B )**( - C ) )**( - D - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Johnson,
!    Multivariate Statistical Simulation,
!    Wiley, 1987.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A <= X
!
!    Input, real ( kind = 8 ) A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x <= a ) then
    pdf = 0.0D+00
  else

    y = ( x - a ) / b

    pdf = ( c * d / b ) * y**( - c - 1.0D+00 ) &
      * ( 1.0D+00 + y**( - c ) )**( - d - 1.0D+00 )

  end if

  return
end
subroutine burr_sample ( a, b, c, d, seed, x )

!*****************************************************************************80
!
!! BURR_SAMPLE samples the Burr PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) d
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call burr_cdf_inv ( cdf, a, b, c, d, x )

  return
end
subroutine burr_variance ( a, b, c, d, variance )

!*****************************************************************************80
!
!! BURR_VARIANCE returns the variance of the Burr PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) k
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) variance

  if ( c <= 2.0D+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BURR_VARIANCE - Warning!'
    write ( *, '(a)' ) '  Variance undefined for C <= 2.'
    variance = huge ( variance )

  else

    k = r8_gamma ( d ) * r8_gamma ( 1.0D+00 - 2.0D+00 / c ) &
      * r8_gamma ( d + 2.0D+00 / c ) &
      - ( r8_gamma ( 1.0D+00 - 1.0D+00 / c ) * r8_gamma ( d + 1.0D+00 / c ) )**2

    variance = k * b * b / ( r8_gamma ( d ) )**2

  end if

  return
end
subroutine c4_normal_01_sample ( seed, x )

!*****************************************************************************80
!
!! C4_NORMAL_01_SAMPLE samples the complex Normal 01 PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 2005
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
!    Output, complex ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 4 ) r4_uniform_01
  real ( kind = 4 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 4 ) v1
  real ( kind = 4 ) v2
  complex ( kind = 4 ) x
  real ( kind = 4 ) x_c
  real ( kind = 4 ) x_r

  v1 = r4_uniform_01 ( seed )
  v2 = r4_uniform_01 ( seed )

  x_r = sqrt ( - 2.0E+00 * log ( v1 ) ) * cos ( 2.0E+00 * pi * v2 )
  x_c = sqrt ( - 2.0E+00 * log ( v1 ) ) * sin ( 2.0E+00 * pi * v2 )

  x = cmplx ( x_r, x_c, kind = 4 )

  return
end
subroutine cardioid_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! CARDIOID_CDF evaluates the Cardioid CDF.
!
!  Discussion:
!
!    The angle X is assumed to lie between A - PI and A + PI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!    A - PI <= X <= A + PI.
!
!    Input, real ( kind = 8 ) A, B, the parameters.
!    -0.5 <= B <= 0.5.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( x <= a - pi ) then
    cdf = 0.0D+00
  else if ( x < a + pi ) then
    cdf = ( pi + x - a + 2.0D+00 * b * sin ( x - a ) ) / ( 2.0D+00 * pi )
  else
    cdf = 1.0D+00
  end if

  return
end
subroutine cardioid_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! CARDIOID_CDF_INV inverts the Cardioid CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0 <= CDF <= 1.
!
!    Input, real ( kind = 8 ) A, B, the parameters.
!    -0.5 <= B <= 0.5.
!
!    Output, real ( kind = 8 ) X, the argument with the given CDF.
!    A - PI <= X <= A + PI.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) fp
  real ( kind = 8 ) fx
  integer ( kind = 4 ) it
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: tol = 0.000001D+00
  real ( kind = 8 ) x

  if ( cdf <= 0.0D+00 ) then

    x = a - pi

  else if ( cdf < 1.0D+00 ) then

    x = a

    it = 0

    do

      fx = cdf - ( pi + x - a + 2.0D+00 * b * sin ( x - a ) ) / ( 2.0D+00 * pi )

      if ( abs ( fx ) < tol ) then
        exit
      end if

      if ( 10 < it ) then
        stop
      end if

      fp = - ( 1.0D+00 + 2.0D+00 * b * cos ( x - a ) ) / ( 2.0D+00 * pi )

      x = x - fx / fp
      x = max ( x, a - pi )
      x = min ( x, a + pi )

      it = it + 1

    end do

  else

    x = a + pi

  end if

  return
end
function cardioid_check ( a, b )

!*****************************************************************************80
!
!! CARDIOID_CHECK checks the parameters of the Cardioid CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    -0.5 <= B <= 0.5.
!
!    Output, logical CARDIOID_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical cardioid_check

  if ( b < -0.5D+00 .or. 0.5D+00 < b ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CARDIOID_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B < -0.5 or 0.5 < B.'
    cardioid_check = .false.
    return
  end if

  cardioid_check = .true.

  return
end
subroutine cardioid_mean ( a, b, mean )

!*****************************************************************************80
!
!! CARDIOID_MEAN returns the mean of the Cardioid PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    -0.5 <= B <= 0.5.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine cardioid_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! CARDIOID_PDF evaluates the Cardioid PDF.
!
!  Discussion:
!
!    The cardioid PDF can be thought of as being applied to points on
!    a circle.  Compare this distribution with the "Cosine PDF".
!
!    PDF(A,B;X) = ( 1 / ( 2 * PI ) ) * ( 1 + 2 * B * COS ( X - A ) )
!    for  A - PI <= X <= A + PI, -1/2 <= B <= 1/2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Nicholas Fisher,
!    Statistical Analysis of Circular Data,
!    Cambridge, 1993.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A - PI <= X <= A + PI.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    -0.5 <= B <= 0.5.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  pdf = ( 1.0D+00 + 2.0D+00 * b * cos ( x - a ) ) / ( 2.0D+00 * pi )

  return
end
subroutine cardioid_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! CARDIOID_SAMPLE samples the Cardioid PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    -0.5 <= B <= 0.5.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!    A - PI <= X <= A + PI.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call cardioid_cdf_inv ( cdf, a, b, x )

  return
end
subroutine cardioid_variance ( a, b, variance )

!*****************************************************************************80
!
!! CARDIOID_VARIANCE returns the variance of the Cardioid PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    -0.5 <= B <= 0.5.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = 0.0D+00

  return
end
subroutine cauchy_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! CAUCHY_CDF evaluates the Cauchy CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  cdf = 0.5D+00 + atan2 ( x - a, b ) / PI

  return
end
subroutine cauchy_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! CAUCHY_CDF_INV inverts the Cauchy CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CAUCHY_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a + b * tan ( pi * ( cdf - 0.5D+00 ) )

  return
end
subroutine cauchy_cdf_values ( n_data, mu, sigma, x, fx )

!*****************************************************************************80
!
!! CAUCHY_CDF_VALUES returns some values of the Cauchy CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = CauchyDistribution [ mu, sigma ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) MU, the mean of the distribution.
!
!    Output, real ( kind = 8 ) SIGMA, the variance of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.8524163823495667D+00, &
    0.9220208696226307D+00, &
    0.9474315432887466D+00, &
    0.6475836176504333D+00, &
    0.6024163823495667D+00, &
    0.5779791303773693D+00, &
    0.5628329581890012D+00, &
    0.6475836176504333D+00, &
    0.5000000000000000D+00, &
    0.3524163823495667D+00, &
    0.2500000000000000D+00 /)
  real ( kind = 8 ) mu
  real ( kind = 8 ), save, dimension ( n_max ) :: mu_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ), save, dimension ( n_max ) :: sigma_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    mu = 0.0D+00
    sigma = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    mu = mu_vec(n_data)
    sigma = sigma_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function cauchy_check ( a, b )

!*****************************************************************************80
!
!! CAUCHY_CHECK checks the parameters of the Cauchy CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, logical CAUCHY_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical cauchy_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CAUCHY_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    cauchy_check = .false.
    return
  end if

  cauchy_check = .true.

  return
end
subroutine cauchy_mean ( a, b, mean )

!*****************************************************************************80
!
!! CAUCHY_MEAN returns the mean of the Cauchy PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine cauchy_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! CAUCHY_PDF evaluates the Cauchy PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = 1 / ( PI * B * ( 1 + ( ( X - A ) / B )^2 ) )
!
!    The Cauchy PDF is also known as the Breit-Wigner PDF.  It
!    has some unusual properties.  In particular, the integrals for the
!    expected value and higher order moments are "singular", in the
!    sense that the limiting values do not exist.  A result can be
!    obtained if the upper and lower limits of integration are set
!    equal to +T and -T, and the limit as T=>INFINITY is taken, but
!    this is a very weak and unreliable sort of limit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  y = ( x - a ) / b

  pdf = 1.0D+00 / ( pi * b * ( 1.0D+00 + y * y ) )

  return
end
subroutine cauchy_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! CAUCHY_SAMPLE samples the Cauchy PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call cauchy_cdf_inv ( cdf, a, b, x )

  return
end
subroutine cauchy_variance ( a, b, variance )

!*****************************************************************************80
!
!! CAUCHY_VARIANCE returns the variance of the Cauchy PDF.
!
!  Discussion:
!
!    The variance of the Cauchy PDF is not well defined.  This routine
!    is made available for completeness only, and simply returns
!    a "very large" number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = huge ( variance )

  return
end
subroutine chi_cdf ( x, a, b, c, cdf )

!*****************************************************************************80
!
!! CHI_CDF evaluates the Chi CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) gamma_inc
  real ( kind = 8 ) p2
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) y

  if ( x <= a ) then

    cdf = 0.0D+00

  else

    y = ( x - a ) / b
    x2 = 0.5D+00 * y * y
    p2 = 0.5D+00 * c

    cdf = gamma_inc ( p2, x2 )

  end if

  return
end
subroutine chi_cdf_inv ( cdf, a, b, c, x )

!*****************************************************************************80
!
!! CHI_CDF_INV inverts the Chi CDF.
!
!  Discussion:
!
!    A simple bisection method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf1
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) cdf3
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  real ( kind = 8 ), parameter :: tol = 0.0001D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHI_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = a
    return
  else if ( 1.0D+00 == cdf ) then
    x = huge ( x )
    return
  end if

  x1 = a
  cdf1 = 0.0D+00

  x2 = a + 1.0D+00

  do

    call chi_cdf ( x2, a, b, c, cdf2 )

    if ( cdf < cdf2 ) then
      exit
    end if

    x2 = a + 2.0D+00 * ( x2 - a )

  end do
!
!  Now use bisection.
!
  it = 0

  do

    it = it + 1

    x3 = 0.5D+00 * ( x1 + x2 )
    call chi_cdf ( x3, a, b, c, cdf3 )

    if ( abs ( cdf3 - cdf ) < tol ) then
      x = x3
      return
    end if

    if ( it_max < it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CHI_CDF_INV - Fatal error!'
      write ( *, '(a)' ) '  Iteration limit exceeded.'
      return
    end if

    if ( sign ( 1.0D+00, cdf3 - cdf ) == sign ( 1.0D+00, cdf1 - cdf ) ) then
      x1 = x3
      cdf1 = cdf3
    else
      x2 = x3
      cdf2 = cdf3
    end if

  end do

  return
end
function chi_check ( a, b, c )

!*****************************************************************************80
!
!! CHI_CHECK checks the parameters of the Chi CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, logical CHI_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical chi_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHI_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.0.'
    chi_check = .false.
    return
  end if

  if ( c <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHI_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C <= 0.0.'
    chi_check = .false.
    return
  end if

  chi_check = .true.

  return
end
subroutine chi_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! CHI_MEAN returns the mean of the Chi PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) MEAN, the mean value.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) mean
  real ( kind = 8 ) r8_gamma

  mean = a + sqrt ( 2.0D+00 ) * b * r8_gamma ( 0.5D+00 * ( c + 1.0D+00 ) ) &
    / r8_gamma ( 0.5D+00 * c )

  return
end
subroutine chi_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! CHI_PDF evaluates the Chi PDF.
!
!  Discussion:
!
!    PDF(A,B,C;X) = EXP ( - 0.5D+00 * ( ( X - A ) / B )^2 )
!      * ( ( X - A ) / B )^( C - 1 ) /
!      ( 2^( 0.5D+00 * C - 1 ) * B * GAMMA ( 0.5D+00 * C ) )
!
!    CHI(A,B,1) is the Half Normal PDF;
!    CHI(0,B,2) is the Rayleigh PDF;
!    CHI(0,B,3) is the Maxwell PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A <= X
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) pdf
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x <= a ) then

    pdf = 0.0D+00

  else

    y = ( x - a ) / b

    pdf = exp ( - 0.5D+00 * y * y ) * y**( c - 1.0D+00 ) / &
      ( 2.0D+00**( 0.5D+00 * c - 1.0D+00 ) * b * r8_gamma ( 0.5D+00 * c ) )

  end if

  return
end
subroutine chi_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! CHI_SAMPLE samples the Chi PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  call chi_square_sample ( c, seed, x )

  x = a + b * sqrt ( x )

  return
end
subroutine chi_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! CHI_VARIANCE returns the variance of the Chi PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) variance

  variance = b * b * ( c - 2.0D+00 * &
    ( r8_gamma ( 0.5D+00 * ( c + 1.0D+00 ) ) / r8_gamma ( 0.5D+00 * c ) )**2 )

  return
end
subroutine chi_square_cdf ( x, a, cdf )

!*****************************************************************************80
!
!! CHI_SQUARE_CDF evaluates the Chi squared CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the random deviate.
!
!    Input, real ( kind = 8 ) A, the parameter of the distribution, usually
!    the number of degrees of freedom.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) b2
  real ( kind = 8 ) c2
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  x2 = 0.5D+00 * x

  a2 = 0.0D+00
  b2 = 1.0D+00
  c2 = 0.5D+00 * a

  call gamma_cdf ( x2, a2, b2, c2, cdf )

  return
end
subroutine chi_square_cdf_inv ( cdf, a, x )

!*****************************************************************************80
!
!! CHI_SQUARE_CDF_INV inverts the Chi squared PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2004
!
!  Author:
!
!    Original FORTAN77 version by Donald Best, Roberts.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Donald Best, Roberts,
!    The Percentage Points of the Chi-Squared Distribution,
!    Algorithm AS 91,
!    Applied Statistics,
!    Volume 24, Number ?, pages 385-390, 1975.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, a value of the chi-squared cumulative
!    probability density function.
!    0.000002 <= CDF <= 0.999998.
!
!    Input, real ( kind = 8 ) A, the parameter of the chi-squared
!    probability density function.  0 < A.
!
!    Output, real ( kind = 8 ) X, the value of the chi-squared random deviate
!    with the property that the probability that a chi-squared random
!    deviate with parameter A is less than or equal to X is CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ), parameter :: aa = 0.6931471806D+00
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: c1 = 0.01D+00
  real ( kind = 8 ), parameter :: c2 = 0.222222D+00
  real ( kind = 8 ), parameter :: c3 = 0.32D+00
  real ( kind = 8 ), parameter :: c4 = 0.4D+00
  real ( kind = 8 ), parameter :: c5 = 1.24D+00
  real ( kind = 8 ), parameter :: c6 = 2.2D+00
  real ( kind = 8 ), parameter :: c7 = 4.67D+00
  real ( kind = 8 ), parameter :: c8 = 6.66D+00
  real ( kind = 8 ), parameter :: c9 = 6.73D+00
  real ( kind = 8 ), parameter :: c10 = 13.32D+00
  real ( kind = 8 ), parameter :: c11 = 60.0D+00
  real ( kind = 8 ), parameter :: c12 = 70.0D+00
  real ( kind = 8 ), parameter :: c13 = 84.0D+00
  real ( kind = 8 ), parameter :: c14 = 105.0D+00
  real ( kind = 8 ), parameter :: c15 = 120.0D+00
  real ( kind = 8 ), parameter :: c16 = 127.0D+00
  real ( kind = 8 ), parameter :: c17 = 140.0D+00
  real ( kind = 8 ), parameter :: c18 = 175.0D+00
  real ( kind = 8 ), parameter :: c19 = 210.0D+00
  real ( kind = 8 ), parameter :: c20 = 252.0D+00
  real ( kind = 8 ), parameter :: c21 = 264.0D+00
  real ( kind = 8 ), parameter :: c22 = 294.0D+00
  real ( kind = 8 ), parameter :: c23 = 346.0D+00
  real ( kind = 8 ), parameter :: c24 = 420.0D+00
  real ( kind = 8 ), parameter :: c25 = 462.0D+00
  real ( kind = 8 ), parameter :: c26 = 606.0D+00
  real ( kind = 8 ), parameter :: c27 = 672.0D+00
  real ( kind = 8 ), parameter :: c28 = 707.0D+00
  real ( kind = 8 ), parameter :: c29 = 735.0D+00
  real ( kind = 8 ), parameter :: c30 = 889.0D+00
  real ( kind = 8 ), parameter :: c31 = 932.0D+00
  real ( kind = 8 ), parameter :: c32 = 966.0D+00
  real ( kind = 8 ), parameter :: c33 = 1141.0D+00
  real ( kind = 8 ), parameter :: c34 = 1182.0D+00
  real ( kind = 8 ), parameter :: c35 = 1278.0D+00
  real ( kind = 8 ), parameter :: c36 = 1740.0D+00
  real ( kind = 8 ), parameter :: c37 = 2520.0D+00
  real ( kind = 8 ), parameter :: c38 = 5040.0D+00
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: cdf_max = 0.999998D+00
  real ( kind = 8 ), parameter :: cdf_min = 0.000002D+00
  real ( kind = 8 ) ch
  real ( kind = 8 ), parameter :: e = 0.0000005D+00
  real ( kind = 8 ) g
  real ( kind = 8 ) gamma_inc
  real ( kind = 8 ) gamma_log
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: it_max = 20
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) q
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) s3
  real ( kind = 8 ) s4
  real ( kind = 8 ) s5
  real ( kind = 8 ) s6
  real ( kind = 8 ) t
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) xx

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHI_SQUARE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    write ( *, '(a,g14.6)' ) '  CDF = ', cdf
    stop
  end if

  if ( cdf < cdf_min ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHI_SQUARE_CDF_INV - Warning!'
    write ( *, '(a)' ) '  CDF < CDF_MIN.'
    write ( *, '(a,g14.6)' ) '  CDF = ', cdf
    write ( *, '(a,g14.6)' ) '  CDF_MIN = ', cdf_min
  end if

  if ( cdf_max < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHI_SQUARE_CDF_INV - Warning!'
    write ( *, '(a)' ) '  CDF_MAX < CDF.'
    write ( *, '(a,g14.6)' ) '  CDF = ', cdf
    write ( *, '(a,g14.6)' ) '  CDF_MAX = ', cdf_max
  end if

  xx = 0.5D+00 * a
  c = xx - 1.0D+00
!
!  Compute Log ( Gamma ( A/2 ) ).
!
  g = gamma_log ( a / 2.0D+00 )
!
!  Starting approximation for small chi-squared.
!
  if ( a < - c5 * log ( cdf ) ) then

    ch = ( cdf * xx * exp ( g + xx * aa ) )**( 1.0D+00 / xx )

    if ( ch < e ) then
      x = ch
      return
    end if
!
!  Starting approximation for A less than or equal to 0.32.
!
  else if ( a <= c3 ) then

    ch = c4
    a2 = log ( 1.0D+00 - cdf )

    do

      q = ch
      p1 = 1.0D+00 + ch * ( c7 + ch )
      p2 = ch * ( c9 + ch * ( c8 + ch ) )

      t = - 0.5D+00 + ( c7 + 2.0D+00 * ch ) / p1 &
        - ( c9 + ch * ( c10 + 3.0D+00 * ch ) ) / p2

      ch = ch - ( 1.0D+00 - exp ( a2 + g + 0.5D+00 * ch + c * aa ) &
        * p2 / p1 ) / t

      if ( abs ( q / ch - 1.0D+00 ) <= c1 ) then
        exit
      end if

    end do
!
!  Call to algorithm AS 111.
!  Note that P has been tested above.
!  AS 241 could be used as an alternative.
!
  else

    call normal_01_cdf_inv ( cdf, x2 )
!
!  Starting approximation using Wilson and Hilferty estimate.
!
    p1 = c2 / a
    ch = a * ( x2 * sqrt ( p1 ) + 1.0D+00 - p1 )**3
!
!  Starting approximation for P tending to 1.
!
    if ( c6 * a + 6.0D+00 < ch ) then
      ch = - 2.0D+00 * ( log ( 1.0D+00 - cdf ) - c * log ( 0.5D+00 * ch ) + g )
    end if

  end if
!
!  Call to algorithm AS 239 and calculation of seven term Taylor series.
!
  do i = 1, it_max

    q = ch
    p1 = 0.5D+00 * ch
    p2 = cdf - gamma_inc ( xx, p1 )
    t = p2 * exp ( xx * aa + g + p1 - c * log ( ch ) )
    b = t / ch
    a2 = 0.5D+00 * t - b * c

    s1 = ( c19 + a2 &
       * ( c17 + a2 &
       * ( c14 + a2 &
       * ( c13 + a2 &
       * ( c12 + a2 &
       *   c11 ) ) ) ) ) / c24

    s2 = ( c24 + a2 &
       * ( c29 + a2 &
       * ( c32 + a2 &
       * ( c33 + a2 &
       *   c35 ) ) ) ) / c37

    s3 = ( c19 + a2 &
       * ( c25 + a2 &
       * ( c28 + a2 &
       *   c31 ) ) ) / c37

    s4 = ( c20 + a2 &
       * ( c27 + a2 &
       *   c34 ) + c &
       * ( c22 + a2 &
       * ( c30 + a2 &
       *   c36 ) ) ) / c38

    s5 = ( c13 + c21 * a2 + c * ( c18 + c26 * a2 ) ) / c37

    s6 = ( c15 + c * ( c23 + c16 * c ) ) / c38

    ch = ch + t * ( 1.0D+00 + 0.5D+00 * t * s1 - b * c &
      * ( s1 - b &
      * ( s2 - b &
      * ( s3 - b &
      * ( s4 - b &
      * ( s5 - b &
      *   s6 ) ) ) ) ) )

    if ( e < abs ( q / ch - 1.0D+00 ) ) then
      x = ch
      return
    end if

  end do

  x = ch
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHI_SQUARE_CDF_INV - Warning!'
  write ( *, '(a)' ) '  Convergence not reached.'

  return
end
subroutine chi_square_cdf_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = ChiSquareDistribution [ df ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  integer ( kind = 4 ) a
  integer ( kind = 4 ), save, dimension ( n_max ) :: a_vec = (/ &
     1,  2,  1,  2, &
     1,  2,  3,  4, &
     1,  2,  3,  4, &
     5,  3,  3,  3, &
     3,  3, 10, 10, &
    10 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7965567455405796D-01, &
    0.4987520807317687D-02, &
    0.1124629160182849D+00, &
    0.9950166250831946D-02, &
    0.4729107431344619D+00, &
    0.1812692469220181D+00, &
    0.5975750516063926D-01, &
    0.1752309630642177D-01, &
    0.6826894921370859D+00, &
    0.3934693402873666D+00, &
    0.1987480430987992D+00, &
    0.9020401043104986D-01, &
    0.3743422675270363D-01, &
    0.4275932955291202D+00, &
    0.6083748237289110D+00, &
    0.7385358700508894D+00, &
    0.8282028557032669D+00, &
    0.8883897749052874D+00, &
    0.1721156299558408D-03, &
    0.3659846827343712D-02, &
    0.1857593622214067D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.01D+00, &
    0.01D+00, &
    0.02D+00, &
    0.02D+00, &
    0.40D+00, &
    0.40D+00, &
    0.40D+00, &
    0.40D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    3.00D+00, &
    4.00D+00, &
    5.00D+00, &
    6.00D+00, &
    1.00D+00, &
    2.00D+00, &
    3.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function chi_square_check ( a )

!*****************************************************************************80
!
!! CHI_SQUARE_CHECK checks the parameter of the central Chi squared PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the distribution.
!    1 <= A.
!
!    Output, logical CHI_SQUARE_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  logical chi_square_check

  if ( a < 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHI_SQUARE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A < 1.0.'
    chi_square_check = .false.
    return
  end if

  chi_square_check = .true.

  return
end
subroutine chi_square_mean ( a, mean )

!*****************************************************************************80
!
!! CHI_SQUARE_MEAN returns the mean of the central Chi squared PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the distribution.
!    1 <= A.
!
!    Output, real ( kind = 8 ) MEAN, the mean value.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine chi_square_pdf ( x, a, pdf )

!*****************************************************************************80
!
!! CHI_SQUARE_PDF evaluates the central Chi squared PDF.
!
!  Discussion:
!
!    PDF(A;X) =
!      EXP ( - X / 2 ) * X**((A-2)/2) / ( 2**(A/2) * GAMMA ( A/2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    1 <= A.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    pdf = 0.0D+00
  else
    b = a / 2.0D+00
    pdf = exp ( - 0.5D+00 * x ) * x**( b - 1.0D+00 ) &
      / ( 2.0D+00**b * r8_gamma ( b ) )
  end if

  return
end
subroutine chi_square_sample ( a, seed, x )

!*****************************************************************************80
!
!! CHI_SQUARE_SAMPLE samples the central Chi squared PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    1 <= A.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) b2
  real ( kind = 8 ) c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: it_max = 100
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  n = int ( a )

  if ( real ( n, kind = 8 ) == a .and. n <= it_max ) then

    x = 0.0D+00
    do i = 1, n
      call normal_01_sample ( seed, x2 )
      x = x + x2 * x2
    end do

  else

    a2 = 0.0D+00
    b2 = 1.0D+00
    c2 = a / 2.0D+00

    call gamma_sample ( a2, b2, c2, seed, x )

    x = 2.0D+00 * x

  end if

  return
end
subroutine chi_square_variance ( a, variance )

!*****************************************************************************80
!
!! CHI_SQUARE_VARIANCE returns the variance of the central Chi squared PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the distribution.
!    1 <= A.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) variance

  variance = 2.0D+00 * a

  return
end
function chi_square_noncentral_check ( a, b )

!*****************************************************************************80
!
!! CHI_SQUARE_NONCENTRAL_CHECK check parameters of noncentral Chi Squared PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the parameter of the PDF.
!    1.0D+00 <= A.
!
!    Input, real ( kind = 8 ) B, the noncentrality parameter of the PDF.
!    0.0D+00 <= B.
!
!    Output, logical CHI_SQUARE_NONCENTRAL_CHECK, is true if the parameters
!    are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical chi_square_noncentral_check

  if ( a < 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHI_SQUARE_NONCENTRAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A < 1.'
    chi_square_noncentral_check = .false.
    return
  end if

  if ( b < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHI_SQUARE_NONCENTRAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B < 0.'
    chi_square_noncentral_check = .false.
    return
  end if

  chi_square_noncentral_check = .true.

  return
end
subroutine chi_square_noncentral_cdf_values ( n_data, df, lambda, x, cdf )

!*****************************************************************************80
!
!! CHI_SQUARE_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NoncentralChiSquareDistribution [ df, lambda ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) DF, the number of degrees of freedom.
!
!    Output, real ( kind = 8 ) LAMBDA, the noncentrality parameter.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) CDF, the noncentral chi CDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) cdf
  real ( kind = 8 ), save, dimension ( n_max ) :: cdf_vec = (/ &
    0.8399444269398261D+00, &
    0.6959060300435139D+00, &
    0.5350879697078847D+00, &
    0.7647841496310313D+00, &
    0.6206436532195436D+00, &
    0.4691667375373180D+00, &
    0.3070884345937569D+00, &
    0.2203818092990903D+00, &
    0.1500251895581519D+00, &
    0.3071163194335791D-02, &
    0.1763982670131894D-02, &
    0.9816792594625022D-03, &
    0.1651753140866208D-01, &
    0.2023419573950451D-03, &
    0.4984476352854074D-06, &
    0.1513252400654827D-01, &
    0.2090414910614367D-02, &
    0.2465021206048452D-03, &
    0.2636835050342939D-01, &
    0.1857983220079215D-01, &
    0.1305736595486640D-01, &
    0.5838039534819351D-01, &
    0.4249784402463712D-01, &
    0.3082137716021596D-01, &
    0.1057878223400849D+00, &
    0.7940842984598509D-01, &
    0.5932010895599639D-01, &
    0.2110395656918684D+00 /)
  integer ( kind = 4 ) df
  integer ( kind = 4 ), save, dimension ( n_max ) :: df_vec = (/ &
      1,   2,   3, &
      1,   2,   3, &
      1,   2,   3, &
      1,   2,   3, &
     60,  80, 100, &
      1,   2,   3, &
     10,  10,  10, &
     10,  10,  10, &
     10,  10,  10, &
      8 /)
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_max ) :: lambda_vec = (/ &
      0.5D+00, &
      0.5D+00, &
      0.5D+00, &
      1.0D+00, &
      1.0D+00, &
      1.0D+00, &
      5.0D+00, &
      5.0D+00, &
      5.0D+00, &
     20.0D+00, &
     20.0D+00, &
     20.0D+00, &
     30.0D+00, &
     30.0D+00, &
     30.0D+00, &
      5.0D+00, &
      5.0D+00, &
      5.0D+00, &
      2.0D+00, &
      3.0D+00, &
      4.0D+00, &
      2.0D+00, &
      3.0D+00, &
      4.0D+00, &
      2.0D+00, &
      3.0D+00, &
      4.0D+00, &
      0.5D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
     3.000D+00, &
    60.000D+00, &
    60.000D+00, &
    60.000D+00, &
     0.050D+00, &
     0.050D+00, &
     0.050D+00, &
     4.000D+00, &
     4.000D+00, &
     4.000D+00, &
     5.000D+00, &
     5.000D+00, &
     5.000D+00, &
     6.000D+00, &
     6.000D+00, &
     6.000D+00, &
     5.000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    lambda = 0.0D+00
    df = 0
    cdf = 0.0D+00
  else
    x = x_vec(n_data)
    lambda = lambda_vec(n_data)
    df = df_vec(n_data)
    cdf = cdf_vec(n_data)
  end if

  return
end
subroutine chi_square_noncentral_mean ( a, b, mean )

!*****************************************************************************80
!
!! CHI_SQUARE_NONCENTRAL_MEAN returns the mean of the noncentral Chi squared PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the parameter of the PDF.
!    1.0D+00 <= A.
!
!    Input, real ( kind = 8 ) B, the noncentrality parameter of the PDF.
!    0.0D+00 <= B.
!
!    Output, real ( kind = 8 ) MEAN, the mean value.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a + b

  return
end
subroutine chi_square_noncentral_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! CHI_SQUARE_NONCENTRAL_SAMPLE samples the noncentral Chi squared PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the parameter of the PDF.
!    1.0D+00 <= A.
!
!    Input, real ( kind = 8 ) B, the noncentrality parameter of the PDF.
!    0.0D+00 <= B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  a1 = a - 1.0D+00

  call chi_square_sample ( a1, seed, x1 )

  a2 = sqrt ( b )
  b2 = 1.0D+00
  call normal_sample ( a2, b2, seed, x2 )

  x = x1 + x2 * x2

  return
end
subroutine chi_square_noncentral_variance ( a, b, variance )

!*****************************************************************************80
!
!! CHI_SQUARE_NONCENTRAL_VARIANCE: variance of the noncentral Chi squared PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    1 <= A.
!
!    Input, real ( kind = 8 ) B, the noncentrality parameter of the PDF.
!    0.0D+00 <= B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance value.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = 2.0D+00 * ( a + 2.0D+00 * b )

  return
end
subroutine circle_sample ( a, b, c, seed, x1, x2 )

!*****************************************************************************80
!
!! CIRCLE_SAMPLE samples points from a circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the circle.
!    The circle is centered at (A,B) and has radius C.
!    0.0D+00 < C.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X1, X2, a sampled point of the circle.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) angle
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) radius_frac
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  radius_frac = r8_uniform_01 ( seed )
  radius_frac = sqrt ( radius_frac )

  angle = 2.0D+00 * pi * r8_uniform_01 ( seed )

  x1 = a + c * radius_frac * cos ( angle )
  x2 = b + c * radius_frac * sin ( angle )

  return
end
subroutine circular_normal_01_mean ( mean )

!*****************************************************************************80
!
!! CIRCULAR_NORMAL_01_MEAN returns the mean of the Circular Normal 01 PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MEAN(2), the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) mean(2)

  mean(1:2) = 0.0D+00

  return
end
subroutine circular_normal_01_pdf ( x, pdf )

!*****************************************************************************80
!
!! CIRCULAR_NORMAL_01_PDF evaluates the Circular Normal 01 PDF.
!
!  Discussion:
!
!    PDF(X) = EXP ( - 0.5D+00 * ( X(1)^2 + X(2)^2 ) ) / ( 2 * PI )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(2), the argument of the PDF.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(2)

  pdf = exp ( - 0.5D+00 * ( x(1)**2 + x(2)**2 ) ) / ( 2.0D+00 * pi )

  return
end
subroutine circular_normal_01_sample ( seed, x )

!*****************************************************************************80
!
!! CIRCULAR_NORMAL_01_SAMPLE samples the Circular Normal 01 PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 1999
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
!    Output, real ( kind = 8 ) X(2), a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) x(2)

  v1 = r8_uniform_01 ( seed )
  v2 = r8_uniform_01 ( seed )

  x(1) = sqrt ( - 2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * pi * v2 )
  x(2) = sqrt ( - 2.0D+00 * log ( v1 ) ) * sin ( 2.0D+00 * pi * v2 )

  return
end
subroutine circular_normal_01_variance ( variance )

!*****************************************************************************80
!
!! CIRCULAR_NORMAL_01_VARIANCE returns the variance of the Circular Normal 01 PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VARIANCE(2), the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) variance(2)

  variance(1) = 1.0D+00
  variance(2) = 1.0D+00

  return
end
subroutine circular_normal_mean ( a, b, mean )

!*****************************************************************************80
!
!! CIRCULAR_NORMAL_MEAN returns the mean of the Circular Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2), a parameter of the PDF, the mean value.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF, the standard deviation.
!
!    Output, real ( kind = 8 ) MEAN(2), the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b
  real ( kind = 8 ) mean(2)

  mean(1:2) = a(1:2)

  return
end
subroutine circular_normal_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! CIRCULAR_NORMAL_PDF evaluates the Circular Normal PDF.
!
!  Discussion:
!
!    PDF(X) = EXP ( - 0.5D+00 * ( ( (X(1)-A(1))^2 + (X(2)-A(2))^2 ) / B^2 ) 
!      / ( 2 * PI * B^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(2), the argument of the PDF.
!
!    Input, real ( kind = 8 ) A(2), a parameter of the PDF, the mean value.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF, the standard deviation.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b
  real ( kind = 8 ) d
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(2)

  d = ( ( x(1) - a(1) )**2 + ( x(2) - a(2) )**2 ) / b**2

  pdf = exp ( - 0.5D+00 * d ) / ( 2.0D+00 * b**2 * pi )

  return
end
subroutine circular_normal_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! CIRCULAR_NORMAL_SAMPLE samples the Circular Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2), a parameter of the PDF, the mean value.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF, the standard deviation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(2), a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) x(2)

  v1 = r8_uniform_01 ( seed )
  v2 = r8_uniform_01 ( seed )

  r = sqrt ( - 2.0D+00 * log ( v1 ) )

  x(1) = a(1) + b * r * cos ( 2.0D+00 * pi * v2 )
  x(2) = a(2) + b * r * sin ( 2.0D+00 * pi * v2 )

  return
end
subroutine circular_normal_variance ( a, b, variance )

!*****************************************************************************80
!
!! CIRCULAR_NORMAL_VARIANCE returns the variance of the Circular Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2), a parameter of the PDF, the mean value.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF, the standard deviation.
!
!    Output, real ( kind = 8 ) VARIANCE(2), the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b
  real ( kind = 8 ) variance(2)

  variance(1) = b**2
  variance(2) = b**2

  return
end
function combinatorial ( n, k )

!*****************************************************************************80
!
!! COMBINATORIAL computes the binomial coefficient C(N,K).
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    April, 1963.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, integer ( kind = 4 ) COMBINATORIAL, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) combinatorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    mx = max ( k, n - k )
    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  combinatorial = value

  return
end
subroutine cosine_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! COSINE_CDF evaluates the Cosine CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameter of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x <= a - pi * b ) then

    cdf = 0.0D+00

  else if ( x <= a + pi * b ) then

    y = ( x - a ) / b

    cdf = 0.5D+00 + ( y + sin ( y ) ) / ( 2.0D+00 * pi )

  else if ( a + pi * b < x ) then

    cdf = 1.0D+00

  end if

  return
end
subroutine cosine_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! COSINE_CDF_INV inverts the Cosine CDF.
!
!  Discussion:
!
!    A simple bisection method is used on the interval
!    [ A - PI * B, A + PI * B ].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf1
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) cdf3
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: tol = 0.0001D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COSINE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = a - pi * b
    return
  else if ( 1.0D+00 == cdf ) then
    x = a + pi * b
    return
  end if

  x1 = a - pi * b
  cdf1 = 0.0D+00

  x2 = a + pi * b
  cdf2 = 1.0D+00
!
!  Now use bisection.
!
  it = 0

  do it = 1, it_max

    x3 = 0.5D+00 * ( x1 + x2 )
    call cosine_cdf ( x3, a, b, cdf3 )

    if ( abs ( cdf3 - cdf ) < tol ) then
      x = x3
      return
    end if

    if ( sign ( 1.0D+00, cdf3 - cdf ) == sign ( 1.0D+00, cdf1 - cdf ) ) then
      x1 = x3
      cdf1 = cdf3
    else
      x2 = x3
      cdf2 = cdf3
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COSINE_CDF_INV - Fatal error!'
  write ( *, '(a)' ) '  Iteration limit exceeded.'

  stop
end
function cosine_check ( a, b )

!*****************************************************************************80
!
!! COSINE_CHECK checks the parameters of the Cosine CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameter of the PDF.
!    0.0D+00 < B.
!
!    Output, logical COSINE_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical cosine_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COSINE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.0'
    cosine_check = .false.
    return
  end if

  cosine_check = .true.

  return
end
subroutine cosine_mean ( a, b, mean )

!*****************************************************************************80
!
!! COSINE_MEAN returns the mean of the Cosine PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine cosine_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! COSINE_PDF evaluates the Cosine PDF.
!
!  Discussion:
!
!    The cosine PDF can be thought of as being applied to points on
!    a circle.
!
!    PDF(A,B;X) = ( 1 / ( 2 * PI * B ) ) * COS ( ( X - A ) / B )
!    for A - PI * B <= X <= A + PI * B
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x < a - pi * b ) then
    pdf = 0.0D+00

  else if ( x <= a + pi * b ) then

    y = ( x - a ) / b

    pdf = 1.0D+00 / ( 2.0D+00 * pi * b ) * cos ( y )

  else if ( a + pi * b < x ) then

    pdf = 0.0D+00

  end if

  return
end
subroutine cosine_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! COSINE_SAMPLE samples the Cosine PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call cosine_cdf_inv ( cdf, a, b, x )

  return
end
subroutine cosine_variance ( a, b, variance )

!*****************************************************************************80
!
!! COSINE_VARIANCE returns the variance of the Cosine PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) variance

  variance = ( pi * pi / 3.0D+00 - 2.0D+00 ) * b * b

  return
end
subroutine coupon_complete_pdf ( type_num, box_num, pdf )

!*****************************************************************************80
!
!! COUPON_COMPLETE_PDF evaluates the Complete Coupon Collection PDF.
!
!  Discussion:
!
!    PDF(TYPE_NUM;BOX_NUM) is the probability that, given an inexhaustible
!    supply of boxes, inside each of which there is one of TYPE_NUM distinct
!    coupons, which are uniformly distributed among the boxes, that it will
!    require opening exactly BOX_NUM boxes to achieve at least one of each
!    kind of coupon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Herbert Wilf,
!    Some New Aspects of the Coupon Collector's Problem,
!    SIAM Review,
!    Volume 48, Number 3, September 2006, pages 549-565.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BOX_NUM, the number of boxes that had to be
!    opened in order to just get at least one of each coupon.
!    0 <= BOX_NUM.  If BOX_NUM < TYPE_NUM, then PDF is surely 0.
!
!    Input, integer ( kind = 4 ) TYPE_NUM, the number of distinct coupons.
!    1 <= TYPE_NUM.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) box_num
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) stirling2_value
  integer ( kind = 4 ) type_num
!
!  Nonsense cases.
!
  if ( box_num < 0 ) then

    pdf = 0.0D+00

  else if ( type_num < 1 ) then

    pdf = 0.0D+00
!
!  Degenerate but meaningful case.
!
  else if ( type_num == 1 ) then

    if ( box_num == 1 ) then
      pdf = 1.0D+00
    else
      pdf = 0.0D+00
    end if
!
!  Easy cases.
!
  else if ( box_num < type_num ) then

    pdf = 0.0D+00
!
!  General case.
!
  else

    factor = 1.0D+00
    do i = 1, type_num
      factor = factor * real ( i ) / real ( type_num )
    end do
    do i = type_num+1, box_num
      factor = factor / real ( type_num )
    end do

    pdf = factor * real ( &
      stirling2_value ( box_num - 1, type_num - 1 ), kind = 8 )

  end if

  return
end
subroutine coupon_mean ( j, type_num, mean )

!*****************************************************************************80
!
!! COUPON_MEAN returns the mean of the Coupon PDF.
!
!  Discussion:
!
!    In this version of the coupon collector's problem, we assume
!    that each box contains 1 coupon, that there are TYPE_NUM distinct types
!    of coupon, uniformly distributed among an inexhaustible supply
!    of boxes, and that the collector's goal is to get J distinct
!    types of coupons by opening one box after another.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J, the number of distinct coupons to be collected.
!    J must be between 1 and TYPE_NUM.
!
!    Input, integer ( kind = 4 ) TYPE_NUM, the number of distinct coupons.
!
!    Output, real ( kind = 8 ) MEAN, the mean number of boxes that
!    must be opened in order to just get J distinct kinds.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) mean
  integer ( kind = 4 ) type_num

  if ( type_num < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COUPON_MEAN - Fatal error!'
    write ( *, '(a)' ) '  Number of distinct coupons desired must be no more'
    write ( *, '(a)' ) '  than the total number of distinct coupons.'
    stop
  end if

  mean = 0.0D+00
  do i = 1, j
    mean = mean + 1.0D+00 / real ( type_num - i + 1, kind = 8 )
  end do
  mean = mean * real ( type_num, kind = 8 )

  return
end
subroutine coupon_simulate ( type_num, seed, coupon, box_num )

!*****************************************************************************80
!
!! COUPON_SIMULATE simulates the coupon collector's problem.
!
!  Discussion:
!
!    The coupon collector needs to collect one of each of TYPE_NUM
!    coupons.  The collector may draw one coupon (or, in some settings,
!    open one box) on each trial, and takes as many trials as necessary
!    to complete the task.  On each trial, the probability of picking
!    any particular type of coupon is always 1 / TYPE_NUM.
!
!    Interesting questions include;
!
!    * what is the expected number of drawings necessary to complete
!      the collection?
!
!    * How does the expected number of drawings necessary to complete
!      the collection vary as TYPE_NUM increases?
!
!    * What is the distribution of the numbers of each type of coupon
!      in a typical collection when it is just completed?
!
!    As TYPE_NUM increases, the number of coupons necessary to be
!    collected in order to get a complete set in any simulation
!    strongly tends to the value TYPE_NUM * LOG ( TYPE_NUM ).
!
!    If TYPE_NUM is 1, the simulation ends with a single drawing.
!
!    If TYPE_NUM is 2, then we may call the coupon taken on the first drawing
!    a "Head", say, and the process then is similar to the question of the
!    length, plus one, of a run of Heads or Tails in coin flipping.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TYPE_NUM, the number of types of coupons.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) COUPON(TYPE_NUM), the number of coupons
!    of each type that were collected during the simulation.
!
!    Output, integer ( kind = 4 ) BOX_NUM, the total number of boxes opened.
!
  implicit none

  integer ( kind = 4 ) type_num

  integer ( kind = 4 ), parameter :: box_max = 2000
  integer ( kind = 4 ) box_num
  integer ( kind = 4 ) coupon(type_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) straight

  coupon(1:type_num) = 0

  straight = 0
  box_num = 0
!
!  Draw another coupon.
!
  do while ( box_num < box_max )

    i = i4_uniform ( 1, type_num, seed )
!
!  Increment the number of I coupons.
!
    coupon(i) = coupon(i) + 1
    box_num = box_num + 1
!
!  If I is the next one we needed, increase STRAIGHT by 1.
!
    if ( i == straight + 1 ) then

      do

        straight = straight + 1
!
!  If STRAIGHT = TYPE_NUM, we have all of them.
!
        if ( type_num <= straight ) then
          return
        end if
!
!  If the next coupon has not been collected, our straight is over.
!
        if ( coupon(straight+1) <= 0 ) then
          exit
        end if

      end do

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COUPON_SIMULATE - Fatal error!'
  write ( *, '(a)' ) '  Maximum number of coupons drawn without success.'

  stop
end
subroutine coupon_variance ( j, type_num, variance )

!*****************************************************************************80
!
!! COUPON_VARIANCE returns the variance of the Coupon PDF.
!
!  Discussion:
!
!    In this version of the coupon collector's problem, we assume
!    that each box contains 1 coupon, that there are TYPE_NUM distinct types
!    of coupon, uniformly distributed among an inexhaustible supply
!    of boxes, and that the collector's goal is to get J distinct
!    types of coupons by opening one box after another.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J, the number of distinct coupons to be
!    collected.  1 <= J <= TYPE_NUM.
!
!    Input, integer ( kind = 4 ) TYPE_NUM, the number of types of coupons.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the number of
!    boxes that must be opened in order to just get J distinct kinds
!    of coupons.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) type_num
  real ( kind = 8 ) variance

  if ( type_num < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COUPON_VARIANCE - Fatal error!'
    write ( *, '(a)' ) '  Number of distinct coupons desired must be no more'
    write ( *, '(a)' ) '  than the total number of distinct coupons.'
    stop
  end if

  variance = 0.0D+00
  do i = 1, j
    variance = variance + real ( i - 1, kind = 8 ) &
      / real ( type_num - i + 1, kind = 8 )**2
  end do
  variance = variance * real ( type_num, kind = 8 )

  return
end
function csc ( theta )

!*****************************************************************************80
!
!! CSC returns the cosecant of X.
!
!  Discussion:
!
!    CSC ( THETA ) = 1.0 / SIN ( THETA )
!
!    CSC is not a built-in function in FORTRAN, and occasionally it
!    is handier, or more concise, to be able to refer to it directly
!    rather than through its definition in terms of the sine function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) THETA, the angle, in radians, whose
!    cosecant is desired.  It must be the case that SIN ( THETA ) is not zero.
!
!    Output, real ( kind = 8 ) CSC, the cosecant of THETA.
!
  implicit none

  real ( kind = 8 ) csc
  real ( kind = 8 ) theta

  csc = sin ( theta )

  if ( csc == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSC - Fatal error!'
    write ( *, '(a,g14.6)' ) '  CSC undefined for THETA = ', theta
    stop
  end if

  csc = 1.0D+00 / csc

  return
end
subroutine deranged_cdf ( x, a, cdf )

!*****************************************************************************80
!
!! DERANGED_CDF evaluates the Deranged CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the maximum number of items in
!    their correct places.
!    0 <= X <= A.
!
!    Input, integer ( kind = 4 ) A, the number of items.
!    1 <= A.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) cnk
  integer ( kind = 4 ) deranged_enum
  integer ( kind = 4 ) dnmk
  real ( kind = 8 ) i4_factorial
  integer ( kind = 4 ) sum2
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  if ( x < 0 .or. a < x ) then
    cdf = 0.0D+00
  else
    sum2 = 0
    do x2 = 0, x
      call binomial_coef ( a, x2, cnk )
      dnmk = deranged_enum ( a-x2 )
      sum2 = sum2 + cnk * dnmk
    end do
    cdf = real ( sum2, kind = 8 ) / i4_factorial ( a )
  end if

  return
end
subroutine deranged_cdf_inv ( cdf, a, x )

!*****************************************************************************80
!
!! DERANGED_CDF_INV inverts the Deranged CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, integer ( kind = 4 ) A, the number of items.
!    1 <= A.
!
!    Output, integer ( kind = 4 ) X, the corresponding argument.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DERANGED_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  cdf2 = 0.0D+00

  do x2 = 0, a

    call deranged_pdf ( x2, a, pdf )

    cdf2 = cdf2 + pdf

    if ( cdf <= cdf2 ) then
      x = x2
      return
    end if

  end do

  x = a

  return
end
function deranged_check ( a )

!*****************************************************************************80
!
!! DERANGED_CHECK checks the parameter of the Deranged PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the total number of items.
!    1 <= A.
!
!    Output, logical DERANGED_CHECK, is true if the parameters are legal.
!
  implicit none

  integer ( kind = 4 ) a
  logical deranged_check

  if ( a < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DERANGED_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A < 1.'
    deranged_check = .false.
    return
  end if

  deranged_check = .true.

  return
end
function deranged_enum ( n )

!*****************************************************************************80
!
!! DERANGED_ENUM returns the number of derangements of N objects.
!
!  Discussion:
!
!    A derangement of N objects is a permutation with no fixed
!    points.  If we symbolize the permutation operation by "P",
!    then for a derangment, P(I) is never equal to I.
!
!      D(0) = 1
!      D(1) = 0
!      D(2) = 1
!      D(N) = (N-1) * ( D(N-1) + D(N-2) )
!
!    or
!
!      D(0) = 1
!      D(1) = 0
!      D(N) = N * D(N-1) + (-1)^N
!
!    D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
!
!    Based on the inclusion/exclusion law.
!
!    D(N) is the number of ways of placing N non-attacking rooks on
!    an N by N chessboard with one diagonal deleted.
!
!    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
!
!    The number of permutations with exactly K items in the right
!    place is COMB(N,K) * D(N-K).
!
!  First values:
!
!     N         D(N)
!     0           1
!     1           0
!     2           1
!     3           2
!     4           9
!     5          44
!     6         265
!     7        1854
!     8       14833
!     9      133496
!    10     1334961
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Output, integer ( kind = 4 ) DERANGED_ENUM, the number of derangements
!    of N objects.
!
  implicit none

  integer ( kind = 4 ) deranged_enum
  integer ( kind = 4 ) dn
  integer ( kind = 4 ) dnm1
  integer ( kind = 4 ) dnm2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  if ( n < 0 ) then

    dn = 0

  else if ( n == 0 ) then

    dn = 1

  else if ( n == 1 ) then

    dn = 0

  else if ( n == 2 ) then

    dn = 1

  else

    dnm1 = 0
    dn = 1

    do i = 3, n
      dnm2 = dnm1
      dnm1 = dn
      dn = ( i - 1 ) * ( dnm1 + dnm2 )
    end do

  end if

  deranged_enum = dn

  return
end
subroutine deranged_mean ( a, mean )

!*****************************************************************************80
!
!! DERANGED_MEAN returns the mean of the Deranged CDF.
!
!  Discussion:
!
!    The mean is computed by straightforward summation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of items.
!    1 <= A.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) mean
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  mean = 0.0D+00
  do x = 1, a
    call deranged_pdf ( x, a, pdf )
    mean = mean + pdf * x
  end do

  return
end
subroutine deranged_pdf ( x, a, pdf )

!*****************************************************************************80
!
!! DERANGED_PDF evaluates the Deranged PDF.
!
!  Discussion:
!
!    PDF(A;X) is the probability that exactly X items will occur in
!    their proper place after a random permutation of A items.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the number of items in their
!    correct places.  0 <= X <= A.
!
!    Input, integer ( kind = 4 ) A, the total number of items.
!    1 <= A.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) cnk
  integer ( kind = 4 ) deranged_enum
  integer ( kind = 4 ) dnmk
  real ( kind = 8 ) i4_factorial
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  if ( x < 0 .or. a < x ) then
    pdf = 0.0D+00
  else
    call binomial_coef ( a, x, cnk )
    dnmk = deranged_enum ( a-x )
    pdf = real ( cnk * dnmk, kind = 8 ) / i4_factorial ( a )
  end if

  return
end
subroutine deranged_sample ( a, seed, x )

!*****************************************************************************80
!
!! DERANGED_SAMPLE samples the Deranged PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of items.
!    1 <= A.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x

  cdf = r8_uniform_01 ( seed )

  call deranged_cdf_inv ( cdf, a, x )

  return
end
subroutine deranged_variance ( a, variance )

!*****************************************************************************80
!
!! DERANGED_VARIANCE returns the variance of the Deranged CDF.
!
!  Discussion:
!
!    The variance is computed by straightforward summation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of items.
!    1 <= A.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) mean
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  real ( kind = 8 ) variance

  call deranged_mean ( a, mean )

  variance = 0.0D+00
  do x = 1, a
    call deranged_pdf ( x, a, pdf )
    variance = variance + pdf * ( x - mean )**2
  end do

  return
end
function digamma ( x )

!*****************************************************************************80
!
!! DIGAMMA calculates the digamma or Psi function.
!
!  Discussion:
!
!    DiGamma ( X ) = d ( log ( Gamma ( X ) ) ) / dX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    Original FORTRAN77 version by Jose Bernardo.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jose Bernardo,
!    Algorithm AS 103:
!    Psi ( Digamma ) Function,
!    Applied Statistics,
!    Volume 25, Number 3, pages 315-317, 1976.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the digamma function.
!    0 < X.
!
!    Output, real ( kind = 8 ) DIGAMMA, the value of the digamma function at X.
!
  implicit none

  real ( kind = 8 ), parameter :: c = 8.5D+00
  real ( kind = 8 ), parameter :: d1 = -0.5772156649D+00
  real ( kind = 8 ) digamma
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: s = 0.00001D+00
  real ( kind = 8 ), parameter :: s3 = 0.08333333333D+00
  real ( kind = 8 ), parameter :: s4 = 0.0083333333333D+00
  real ( kind = 8 ), parameter :: s5 = 0.003968253968D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  The argument must be positive.
!
  if ( x <= 0.0D+00 ) then

    digamma = 0.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGAMMA - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
!
!  Use approximation if argument <= S.
!
  else if ( x <= s ) then

    digamma = d1 - 1.0D+00 / x
!
!  Reduce the argument to DIGAMMA(X + N) where C <= (X + N).
!
  else

    digamma = 0.0D+00
    y = x

    do while ( y < c )
      digamma = digamma - 1.0D+00 / y
      y = y + 1.0D+00
    end do
!
!  Use Stirling's (actually de Moivre's) expansion if C < argument.
!
    r = 1.0D+00 / ( y * y )
    digamma = digamma + log ( y ) - 0.5D+00 / y &
      - r * ( s3 - r * ( s4 - r * s5 ) )

  end if

  return
end
subroutine dipole_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! DIPOLE_CDF evaluates the Dipole CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI
!    is interesting, and -1.0D+00 <= B <= 1.0.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  cdf = 0.5D+00 + ( 1.0D+00 / pi ) * atan ( x ) &
    + b * b * ( x * cos ( 2.0D+00 * a ) &
    - sin ( 2.0D+00 * a ) ) / ( pi * ( 1.0D+00 + x * x ) )

  return
end
subroutine dipole_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! DIPOLE_CDF_INV inverts the Dipole CDF.
!
!  Discussion:
!
!    A simple bisection method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    -1.0D+00 <= B <= 1.0.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf1
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) cdf3
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  real ( kind = 8 ), parameter :: tol = 0.0001D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIPOLE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = - huge ( x )
    return
  else if ( 1.0D+00 == cdf ) then
    x = huge ( x )
    return
  end if
!
!  Seek X1 < X < X2.
!
  x1 = - 1.0D+00

  do

    call dipole_cdf ( x1, a, b, cdf1 )

    if ( cdf1 <= cdf ) then
      exit
    end if

    x1 = 2.0D+00 * x1

  end do

  x2 = 1.0D+00

  do

    call dipole_cdf ( x2, a, b, cdf2 )

    if ( cdf <= cdf2 ) then
      exit
    end if

    x2 = 2.0D+00 * x2

  end do
!
!  Now use bisection.
!
  it = 0

  do

    it = it + 1

    x3 = 0.5D+00 * ( x1 + x2 )
    call dipole_cdf ( x3, a, b, cdf3 )

    if ( abs ( cdf3 - cdf ) < tol ) then
      x = x3
      exit
    end if

    if ( it_max < it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIPOLE_CDF_INV - Fatal error!'
      write ( *, '(a)' ) '  Iteration limit exceeded.'
      stop
    end if

    if ( sign ( 1.0D+00, cdf3 - cdf ) == sign ( 1.0D+00, cdf1 - cdf ) ) then
      x1 = x3
      cdf1 = cdf3
    else
      x2 = x3
      cdf2 = cdf3
    end if

  end do

  return
end
function dipole_check ( a, b )

!*****************************************************************************80
!
!! DIPOLE_CHECK checks the parameters of the Dipole CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI
!    is interesting, and -1.0D+00 <= B <= 1.0.
!
!    Output, logical DIPOLE_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical dipole_check

  if ( b < -1.0D+00 .or. 1.0D+00 < b ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIPOLE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  -1.0D+00 <= B <= 1.0D+00 is required.'
    write ( *, '(a,g14.6)' ) '  The input B = ', b
    dipole_check = .false.
    return
  end if

  dipole_check = .true.

  return
end
subroutine dipole_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! DIPOLE_PDF evaluates the Dipole PDF.
!
!  Discussion:
!
!    PDF(A,B;X) =
!        1 / ( PI * ( 1 + X^2 ) )
!      + B^2 * ( ( 1 - X^2 ) * cos ( 2 * A ) + 2 * X * sin ( 2 * A ) )
!      / ( PI * ( 1 + X^2 )^2 )
!
!    Densities of this kind commonly occur in the analysis of resonant
!    scattering of elementary particles.
!
!    DIPOLE_PDF(A,0;X) = CAUCHY_PDF(A;X)
!
!    A = 0, B = 1 yields the single channel dipole distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Knop,
!    Algorithm 441,
!    Random Deviates from the Dipole Distribution,
!    ACM Transactions on Mathematical Software,
!    Volume 16, Number 1, 1973, page 51.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI
!      is interesting,
!    and -1.0D+00 <= B <= 1.0.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  pdf = 1.0D+00 / ( pi * ( 1.0D+00 + x * x ) ) &
    + b * b * ( ( 1.0D+00 - x * x ) * cos ( 2.0D+00 * a ) &
    + 2.0D+00 * x * sin ( 2.0D+00 * x ) ) / ( pi * ( 1.0D+00 + x * x )**2 )

  return
end
subroutine dipole_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! DIPOLE_SAMPLE samples the Dipole PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Knop,
!    Algorithm 441,
!    Random Deviates from the Dipole Distribution,
!    ACM Transactions on Mathematical Software,
!    Volume 16, Number 1, 1973, page 51.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI
!      is interesting,
!    and -1.0D+00 <= B <= 1.0D+00.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) c2
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
!
!  Find (X1,X2) at random in a circle.
!
  a2 = b * sin ( a )
  b2 = b * cos ( a )
  c2 = 1.0D+00

  call circle_sample ( a2, b2, c2, seed, x1, x2 )
!
!  The dipole variate is the ratio X1 / X2.
!
  x = x1 / x2

  return
end
function dirichlet_check ( n, a )

!*****************************************************************************80
!
!! DIRICHLET_CHECK checks the parameters of the Dirichlet PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, real ( kind = 8 ) A(N), the probabilities for each component.
!    Each A(I) should be positive.
!
!    Output, logical DIRICHLET_CHECK, is true if the parameters are legal.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical dirichlet_check
  integer ( kind = 4 ) i
  logical positive

  positive = .false.

  do i = 1, n

    if ( a(i) <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIRICHLET_CHECK - Fatal error!'
      write ( *, '(a)' ) '  A(I) <= 0.'
      write ( *, '(a,i8)' ) '  For I = ', i
      write ( *, '(a,g14.6)' ) '  A(I) = ', a(i)
      dirichlet_check = .false.
      return
    else if ( 0.0D+00 < a(i) ) then
      positive = .true.
    end if

  end do

  if ( .not. positive ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRICHLET_CHECK - Fatal error!'
    write ( *, '(a)' ) '  All parameters are zero!'
    dirichlet_check = .false.
    return
  end if

  dirichlet_check = .true.

  return
end
subroutine dirichlet_mean ( n, a, mean )

!*****************************************************************************80
!
!! DIRICHLET_MEAN returns the means of the Dirichlet PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, real ( kind = 8 ) A(N), the probabilities for each component.
!    Each A(I) should be positive.
!
!    Output, real ( kind = 8 ) MEAN(N), the means of the PDF.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) mean(n)

  mean(1:n) = a(1:n)

  call r8vec_unit_sum ( n, mean )

  return
end
function dirichlet_mix_check ( comp_num, elem_num, a, comp_weight )

!*****************************************************************************80
!
!! DIRICHLET_MIX_CHECK checks the parameters of a Dirichlet mixture PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in the
!    Dirichlet
!    mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    Input, integer ( kind = 4 ) ELEM_NUM, the number of elements of an
!    observation.
!
!    Input, real ( kind = 8 ) A(ELEM_NUM,COMP_NUM), the probabilities
!    for element ELEM_NUM in component COMP_NUM.
!    Each A(I,J) should be positive.
!
!    Input, real ( kind = 8 ) COMP_WEIGHT(COMP_NUM), the mixture weights of
!    the densities.  These do not need to be normalized.  The weight of a
!    given component is the relative probability that that component will
!    be used to generate the sample.
!
!    Output, logical DIRICHLET_MIX_CHECK, is true if the parameters are legal.
!
  implicit none

  integer ( kind = 4 ) comp_num
  integer ( kind = 4 ) elem_num

  real ( kind = 8 ) a(elem_num,comp_num)
  integer ( kind = 4 ) comp_i
  real ( kind = 8 ) comp_weight(comp_num)
  logical dirichlet_mix_check
  integer ( kind = 4 ) elem_i
  logical positive

  do comp_i = 1, comp_num

    do elem_i = 1, elem_num
      if ( a(elem_i,comp_i) <= 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DIRICHLET_MIX_CHECK - Fatal error!'
        write ( *, '(a)' ) '  A(ELEM,COMP) <= 0.'
        write ( *, '(a,i8)' ) '  COMP = ', comp_i
        write ( *, '(a,i8)' ) '  ELEM = ', elem_i
        write ( *, '(a,g14.6)' ) '  A(COMP,ELEM) = ', a(elem_i,comp_i)
        dirichlet_mix_check = .false.
        return
      end if
    end do

  end do

  positive = .false.

  do comp_i = 1, comp_num

    if ( comp_weight(comp_i) < 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIRICHLET_MIX_CHECK - Fatal error!'
      write ( *, '(a)' ) '  COMP_WEIGHT(COMP) < 0.'
      write ( *, '(a,i8)' ) '  COMP = ', comp_i
      write ( *, '(a,g14.6)' ) '  COMP_WEIGHT(COMP) = ', comp_weight(comp_i)
      dirichlet_mix_check = .false.
      return
    else if ( 0.0D+00 < comp_weight(comp_i) ) then
      positive = .true.
    end if

  end do

  if ( .not. positive ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRICHLET_MIX_CHECK - Fatal error!'
    write ( *, '(a)' ) '  All component weights are zero.'
    dirichlet_mix_check = .false.
    return
  end if

  dirichlet_mix_check = .true.

  return
end
subroutine dirichlet_mix_mean ( comp_num, elem_num, a, comp_weight, &
  mean )

!*****************************************************************************80
!
!! DIRICHLET_MIX_MEAN returns the means of a Dirichlet mixture PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in the
!    Dirichlet mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    Input, integer ( kind = 4 ) ELEM_NUM, the number of elements of an observation.
!
!    Input, real ( kind = 8 ) A(ELEM_NUM,COMP_NUM), the probabilities for
!    element ELEM_NUM in component COMP_NUM.
!    Each A(I,J) should be positive.
!
!    Input, real ( kind = 8 ) COMP_WEIGHT(COMP_NUM), the mixture weights of
!    the densities.  These do not need to be normalized.  The weight of a
!    given component is the relative probability that that component will
!    be used to generate the sample.
!
!    Output, real ( kind = 8 ) MEAN(ELEM_NUM), the means for each element.
!
  implicit none

  integer ( kind = 4 ) comp_num
  integer ( kind = 4 ) elem_num

  real ( kind = 8 ) a(elem_num,comp_num)
  integer ( kind = 4 ) comp_i
  real ( kind = 8 ) comp_mean(elem_num)
  real ( kind = 8 ) comp_weight(comp_num)
  real ( kind = 8 ) comp_weight_sum
  real ( kind = 8 ) mean(elem_num)

  comp_weight_sum = sum ( comp_weight )

  mean(1:elem_num) = 0.0D+00

  do comp_i = 1, comp_num
    call dirichlet_mean ( elem_num, a(1,comp_i), comp_mean )
    mean(1:elem_num) = mean(1:elem_num) &
      + comp_weight(comp_i) * comp_mean(1:elem_num)
  end do

  mean(1:elem_num) = mean(1:elem_num) / comp_weight_sum

  return
end
subroutine dirichlet_mix_pdf ( x, comp_num, elem_num, a, &
  comp_weight, pdf )

!*****************************************************************************80
!
!! DIRICHLET_MIX_PDF evaluates a Dirichlet mixture PDF.
!
!  Discussion:
!
!    The PDF is a weighted sum of Dirichlet PDF's.  Each PDF is a
!    "component", with an associated weight.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(ELEM_NUM), the argument of the PDF.
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in the Dirichlet
!    mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    Input, integer ( kind = 4 ) ELEM_NUM, the number of elements of an observation.
!
!    Input, real ( kind = 8 ) A(ELEM_NUM,COMP_NUM), the probabilities for
!    element ELEM_NUM in component COMP_NUM.
!    Each A(I,J) should be positive.
!
!    Input, real ( kind = 8 ) COMP_WEIGHT(COMP_NUM), the mixture weights of
!    the densities.  These do not need to be normalized.  The weight of a
!    given component is the relative probability that that component will
!    be used to generate the sample.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) comp_num
  integer ( kind = 4 ) elem_num

  real ( kind = 8 ) a(elem_num,comp_num)
  integer ( kind = 4 ) comp_i
  real ( kind = 8 ) comp_pdf
  real ( kind = 8 ) comp_weight(comp_num)
  real ( kind = 8 ) comp_weight_sum
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x(elem_num)

  comp_weight_sum = sum ( comp_weight )

  pdf = 0.0D+00
  do comp_i = 1, comp_num

    call dirichlet_pdf ( x, elem_num, a(1,comp_i), comp_pdf )

    pdf = pdf + comp_weight(comp_i) * comp_pdf / comp_weight_sum

  end do

  return
end
subroutine dirichlet_mix_sample ( comp_num, elem_num, a, &
  comp_weight, seed, comp, x )

!*****************************************************************************80
!
!! DIRICHLET_MIX_SAMPLE samples a Dirichlet mixture PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) COMP_NUM, the number of components in the 
!    Dirichlet mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    Input, integer ( kind = 4 ) ELEM_NUM, the number of elements of 
!    an observation.
!
!    Input, real ( kind = 8 ) A(ELEM_NUM,COMP_NUM), the probabilities for
!    element ELEM_NUM in component COMP_NUM.
!    Each A(I,J) should be positive.
!
!    Input, real ( kind = 8 ) COMP_WEIGHT(COMP_NUM), the mixture weights of
!    the densities.  These do not need to be normalized.  The weight of a
!    given component is the relative probability that that component will
!    be used to generate the sample.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) COMP, the index of the component of the 
!    Dirichlet mixture that was chosen to generate the sample.
!
!    Output, real ( kind = 8 ) X(ELEM_NUM), a sample of the PDF.
!
  implicit none

  integer ( kind = 4 ) comp_num
  integer ( kind = 4 ) elem_num

  real ( kind = 8 ) a(elem_num,comp_num)
  integer ( kind = 4 ) comp
  real ( kind = 8 ) comp_weight(comp_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(elem_num)
!
!  Choose a particular density component COMP.
!
  call discrete_sample ( comp_num, comp_weight, seed, comp )
!
!  Sample the density number COMP.
!
  call dirichlet_sample ( elem_num, a(1,comp), seed, x )

  return
end
subroutine dirichlet_moment2 ( n, a, m2 )

!*****************************************************************************80
!
!! DIRICHLET_MOMENT2 returns the second moments of the Dirichlet PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, real ( kind = 8 ) A(N), the probabilities for each component.
!    Each A(I) should be positive.
!
!    Output, real ( kind = 8 ) M2(N,N), the second moments of the PDF.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_sum
  real ( kind = 8 ) m2(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  a_sum = sum ( a(1:n) )

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        m2(i,j) = a(i) * ( a(i) + 1.0D+00 ) / ( a_sum * ( a_sum + 1.0D+00 ) )
      else
        m2(i,j) = a(i) * a(j) / ( a_sum * ( a_sum + 1.0D+00 ) )
      end if
    end do
  end do

  return
end
subroutine dirichlet_multinomial_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! DIRICHLET_MULTINOMIAL_PDF evaluates a Dirichlet Multinomial PDF.
!
!  Discussion:
!
!    PDF(A,B,C;X) = Comb(A,B,X) * ( Gamma(C_Sum) / Gamma(C_Sum+A) )
!      Product ( 1 <= I <= B ) Gamma(C(I)+X(I)) / Gamma(C(I))
!
!    where:
!
!      Comb(A,B,X) is the multinomial coefficient C( A; X(1), X(2), ..., X(B) ),
!      C_Sum = Sum ( 1 <= I <= B ) C(I)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kenneth Lange,
!    Mathematical and Statistical Methods for Genetic Analysis,
!    Springer, 1997, page 45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X(B); X(I) counts the number of occurrences of
!    outcome I, out of the total of A trials.
!
!    Input, integer ( kind = 4 ) A, the total number of trials.
!
!    Input, integer ( kind = 4 ) B, the number of different possible outcomes on
!    one trial.
!
!    Input, real ( kind = 8 ) C(B); C(I) is the Dirichlet parameter associated
!    with outcome I.
!
!    Output, real ( kind = 8 ) PDF, the value of the Dirichlet multinomial PDF.
!
  implicit none

  integer ( kind = 4 ) b

  integer ( kind = 4 ) a
  real ( kind = 8 ) c(b)
  real ( kind = 8 ) c_sum
  real ( kind = 8 ) gamma_log
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  real ( kind = 8 ) pdf_log
  integer ( kind = 4 ) x(b)

  c_sum = sum ( c(1:b) )

  pdf_log = - gamma_log ( c_sum + real ( a, kind = 8 ) ) + gamma_log ( c_sum ) &
            + gamma_log ( real ( a + 1, kind = 8 ) )

  do i = 1, b
    pdf_log = pdf_log + gamma_log ( c(i) + real ( x(i), kind = 8 ) ) &
      - gamma_log ( c(i) ) - gamma_log ( real ( x(i) + 1, kind = 8 ) )
  end do

  pdf = exp ( pdf_log )

  return
end
subroutine dirichlet_pdf ( x, n, a, pdf )

!*****************************************************************************80
!
!! DIRICHLET_PDF evaluates the Dirichlet PDF.
!
!  Discussion:
!
!    PDF(N,A;X) = Product ( 1 <= I <= N ) X(I)**( A(I) - 1 )
!      * Gamma ( A_SUM ) / A_PROD
!
!    where
!
!      0 <= A(I) for all I;
!      0 <= X(I) for all I;
!      Sum ( 1 <= I <= N ) X(I) = 1;
!      A_SUM = Sum ( 1 <= I <= N ) A(I).
!      A_PROD = Product ( 1 <= I <= N ) Gamma ( A(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the PDF.  Each X(I) should
!    be greater than 0.0D+00, and the X(I)'s must add up to 1.0.
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, real ( kind = 8 ) A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be
!    positive.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_prod
  real ( kind = 8 ) a_sum
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ), parameter :: tol = 0.0001D+00
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_sum

  do i = 1, n
    if ( x(i) <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIRICHLET_PDF - Fatal error!'
      write ( *, '(a)' ) '  X(I) <= 0.'
    end if
  end do

  x_sum = sum ( x(1:n) )

  if ( tol < abs ( x_sum - 1.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRICHLET_PDF - Fatal error!'
    write ( *, '(a)' ) '  SUM X(I) =/= 1.'
  end if

  a_sum = sum ( a(1:n) )

  a_prod = 1.0D+00
  do i = 1, n
    a_prod = a_prod * r8_gamma ( a(i) )
  end do

  pdf = r8_gamma ( a_sum ) / a_prod
  do i = 1, n
    pdf = pdf * x(i)**( a(i) - 1.0D+00 )
  end do

  return
end
subroutine dirichlet_sample ( n, a, seed, x )

!*****************************************************************************80
!
!! DIRICHLET_SAMPLE samples the Dirichlet PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 169.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, real ( kind = 8 ) A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be
!    positive.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the PDF.  The entries
!    of X should sum to 1.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a2
  real ( kind = 8 ) b2
  real ( kind = 8 ) c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  a2 = 0.0D+00
  b2 = 1.0D+00

  do i = 1, n
    c2 = a(i)
    call gamma_sample ( a2, b2, c2, seed, x(i) )
  end do
!
!  Rescale the vector to have unit sum.
!
  call r8vec_unit_sum ( n, x )

  return
end
subroutine dirichlet_variance ( n, a, variance )

!*****************************************************************************80
!
!! DIRICHLET_VARIANCE returns the variances of the Dirichlet PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, real ( kind = 8 ) A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
!    Output, real ( kind = 8 ) VARIANCE(N), the variances of the PDF.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_sum
  integer ( kind = 4 ) i
  real ( kind = 8 ) variance(n)

  a_sum = sum ( a(1:n) )

  do i = 1, n
    variance(i) = a(i) * ( a_sum - a(i) ) / ( a_sum**2 * ( a_sum + 1.0D+00 ) )
  end do

  return
end
subroutine discrete_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! DISCRETE_CDF evaluates the Discrete CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the item whose probability is desired.
!
!    Input, integer ( kind = 4 ) A, the number of probabilities assigned.
!
!    Input, real ( kind = 8 ) B(A), the relative probabilities of outcomes
!    1 through A.  Each entry must be nonnegative.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) x

  if ( x < 1 ) then
    cdf = 0.0D+00
  else if ( x < a ) then
    cdf = sum ( b(1:x) ) / sum ( b(1:a) )
  else if ( a <= x ) then
    cdf = 1.0D+00
  end if

  return
end
subroutine discrete_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! DISCRETE_CDF_INV inverts the Discrete CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, integer ( kind = 4 ) A, the number of probabilities assigned.
!
!    Input, real ( kind = 8 ) B(A), the relative probabilities of outcomes
!    1 through A.  Each entry must be nonnegative.
!
!    Output, integer ( kind = 4 ) X, the corresponding argument for which
!    CDF(X-1) < CDF <= CDF(X)
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) b_sum
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cum
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DISCRETE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  b_sum = sum ( b(1:a) )

  cum = 0.0D+00

  do j = 1, a

    cum = cum + b(j) / b_sum

    if ( cdf <= cum ) then
      x = j
      return
    end if

  end do

  x = a

  return
end
function discrete_check ( a, b )

!*****************************************************************************80
!
!! DISCRETE_CHECK checks the parameters of the Discrete CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of probabilities assigned.
!
!    Input, real ( kind = 8 ) B(A), the relative probabilities of
!    outcomes 1 through A.  Each entry must be nonnegative.
!
!    Output, logical DISCRETE_CHECK, is true if the parameters are legal.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) b_sum
  logical discrete_check
  integer ( kind = 4 ) j

  do j = 1, a
    if ( b(j) < 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISCRETE_CHECK - Fatal error!'
      write ( *, '(a)' ) '  Negative probabilities not allowed.'
      discrete_check = .false.
      return
    end if
  end do

  b_sum = sum ( b(1:a) )

  if ( b_sum == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DISCRETE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Total probablity is zero.'
    discrete_check = .false.
    return
  end if

  discrete_check = .true.

  return
end
subroutine discrete_mean ( a, b, mean )

!*****************************************************************************80
!
!! DISCRETE_MEAN evaluates the mean of the Discrete PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of probabilities assigned.
!
!    Input, real ( kind = 8 ) B(A), the relative probabilities of
!    outcomes 1 through A.  Each entry must be nonnegative.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) b_sum
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean

  b_sum = sum ( b(1:a) )

  mean = 0.0D+00
  do j = 1, a
    mean = mean + real ( j, kind = 8 ) * b(j)
  end do

  mean = mean / b_sum

  return
end
subroutine discrete_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! DISCRETE_PDF evaluates the Discrete PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = B(X) if 1 <= X <= A
!                = 0    otherwise
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the item whose probability is desired.
!
!    Input, integer ( kind = 4 ) A, the number of probabilities assigned.
!
!    Input, real ( kind = 8 ) B(A), the relative probabilities of
!    outcomes 1 through A.  Each entry must be nonnegative.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) b_sum
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  b_sum = sum ( b(1:a) )

  if ( 1 <= x .and. x <= a ) then
    pdf = b(x) / b_sum
  else
    pdf = 0.0D+00
  end if

  return
end
subroutine discrete_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! DISCRETE_SAMPLE samples the Discrete PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of probabilities assigned.
!
!    Input, real ( kind = 8 ) B(A), the relative probabilities of
!    outcomes 1 through A.  Each entry must be nonnegative.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) b_sum
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x

  b_sum = sum ( b(1:a) )

  cdf = r8_uniform_01 ( seed )

  call discrete_cdf_inv ( cdf, a, b, x )

  return
end
subroutine discrete_variance ( a, b, variance )

!*****************************************************************************80
!
!! DISCRETE_VARIANCE evaluates the variance of the Discrete PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of probabilities assigned.
!
!    Input, real ( kind = 8 ) B(A), the relative probabilities of
!    outcomes 1 through A.  Each entry must be nonnegative.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) b_sum
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance

  b_sum = sum ( b(1:a) )

  mean = 0.0D+00
  do j = 1, a
    mean = mean + real ( j, kind = 8 ) * b(j)
  end do

  mean = mean / b_sum

  variance = 0.0D+00
  do j = 1, a
    variance = variance + b(j) * ( j - mean )**2
  end do

  variance = variance / b_sum

  return
end
function e_constant ( )

!*****************************************************************************80
!
!! E_CONSTANT returns the value of E.
!
!  Discussion:
!
!   "E" was named in honor of Euler.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) E_CONSTANT, the base of the natural
!    logarithm system.
!
  implicit none

  real ( kind = 8 ) e_constant

  e_constant = 2.71828182845904523536028747135266249775724709369995D+00

  return
end
subroutine empirical_discrete_cdf ( x, a, b, c, cdf )

!*****************************************************************************80
!
!! EMPIRICAL_DISCRETE_CDF evaluates the Empirical Discrete CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, integer ( kind = 4 ) A, the number of values.
!    0 < A.
!
!    Input, real ( kind = 8 ) B(A), the weights of each value.
!    0 <= B(1:A) and at least one value is nonzero.
!
!    Input, real ( kind = 8 ) C(A), the values.
!    The values must be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) bsum
  real ( kind = 8 ) c(a)
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  cdf = 0.0D+00

  bsum = sum ( b(1:a) )

  do i = 1, a

    if ( x < c(i) ) then
      return
    end if

    cdf = cdf + b(i) / bsum

  end do

  return
end
subroutine empirical_discrete_cdf_inv ( cdf, a, b, c, x )

!*****************************************************************************80
!
!! EMPIRICAL_DISCRETE_CDF_INV inverts the Empirical Discrete CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, integer ( kind = 4 ) A, the number of values.
!    0 < A.
!
!    Input, real ( kind = 8 ) B(A), the weights of each value.
!    0 <= B(1:A) and at least one value is nonzero.
!
!    Input, real ( kind = 8 ) C(A), the values.
!    The values must be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) X, the smallest argument whose CDF is greater
!    than or equal to CDF.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) bsum
  real ( kind = 8 ) c(a)
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf2
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EMPIRICAL_DISCRETE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  bsum = sum ( b(1:a) )

  x = c(1)
  cdf2 = b(1) / bsum

  do i = 2, a

    if ( cdf <= cdf2 ) then
      return
    end if

    x = c(i)
    cdf2 = cdf2 + b(i) / bsum

  end do

  return
end
function empirical_discrete_check ( a, b, c )

!*****************************************************************************80
!
!! EMPIRICAL_DISCRETE_CHECK checks the parameters of the Empirical Discrete CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of values.
!    0 < A.
!
!    Input, real ( kind = 8 ) B(A), the weights of each value.
!    0 <= B(1:A) and at least one value is nonzero.
!
!    Input, real ( kind = 8 ) C(A), the values.
!    The values must be distinct and in ascending order.
!
!    Output, logical EMPIRICAL_DISCRETE_CHECK, is true if the parameters
!    are legal.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) c(a)
  logical empirical_discrete_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( a <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EMPIRICAL_DISCRETE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A must be positive.'
    write ( *, '(a,i12)' ) '  Input A = ', a
    write ( *, '(a)' ) '  A is the number of weights.'
    empirical_discrete_check = .false.
    return
  end if

  if ( any ( b(1:a) < 0.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EMPIRICAL_DISCRETE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some B(*) < 0.'
    write ( *, '(a)' ) '  But all B values must be nonnegative.'
    empirical_discrete_check = .false.
    return
  end if

  if ( all ( b(1:a) == 0.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EMPIRICAL_DISCRETE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  All B(*) = 0.'
    write ( *, '(a)' ) '  But at least one B values must be nonzero.'
    empirical_discrete_check = .false.
    return
  end if

  do i = 1, a
    do j = i+1, a
      if ( c(i) == c(j) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EMPIRICAL_DISCRETE_CHECK - Fatal error!'
        write ( *, '(a)' ) '  All values C must be unique.'
        write ( *, '(a)' ) '  But at least two values are identical.'
        empirical_discrete_check = .false.
        return
      end if
    end do
  end do

  do i = 1, a-1
    if ( c(i+1) < c(i) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EMPIRICAL_DISCRETE_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The values in C must be in ascending order.'
      empirical_discrete_check = .false.
      return
    end if
  end do

  empirical_discrete_check = .true.

  return
end
subroutine empirical_discrete_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! EMPIRICAL_DISCRETE_MEAN returns the mean of the Empirical Discrete PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of values.
!    0 < A.
!
!    Input, real ( kind = 8 ) B(A), the weights of each value.
!    0 <= B(1:A) and at least one value is nonzero.
!
!    Input, real ( kind = 8 ) C(A), the values.
!    The values must be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) c(a)
  real ( kind = 8 ) mean

  mean = dot_product ( b(1:a), c(1:a) ) / sum ( b(1:a) )

  return
end
subroutine empirical_discrete_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! EMPIRICAL_DISCRETE_PDF evaluates the Empirical Discrete PDF.
!
!  Discussion:
!
!    A set of A values C(1:A) are assigned nonnegative weights B(1:A),
!    with at least one B nonzero.  The probability of C(I) is the
!    value of B(I) divided by the sum of the weights.
!
!    The C's must be distinct, and given in ascending order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, integer ( kind = 4 ) A, the number of values.
!    0 < A.
!
!    Input, real ( kind = 8 ) B(A), the weights of each value.
!    0 <= B(1:A) and at least one value is nonzero.
!
!    Input, real ( kind = 8 ) C(A), the values.
!    The values must be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) c(a)
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  do i = 1, a
    if ( x == c(i) ) then
      pdf = b(i) / sum ( b(1:a) )
      return
    end if
  end do

  pdf = 0.0D+00

  return
end
subroutine empirical_discrete_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! EMPIRICAL_DISCRETE_SAMPLE samples the Empirical Discrete PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of values.
!    0 < A.
!
!    Input, real ( kind = 8 ) B(A), the weights of each value.
!    0 <= B(1:A) and at least one value is nonzero.
!
!    Input, real ( kind = 8 ) C(A), the values.
!    The values must be distinct and in ascending order.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) c(a)
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call empirical_discrete_cdf_inv ( cdf, a, b, c, x )

  return
end
subroutine empirical_discrete_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! EMPIRICAL_DISCRETE_VARIANCE returns the variance of the Empirical Discrete PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of values.
!    0 < A.
!
!    Input, real ( kind = 8 ) B(A), the weights of each value.
!    0 <= B(1:A) and at least one value is nonzero.
!
!    Input, real ( kind = 8 ) C(A), the values.
!    The values must be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  integer ( kind = 4 ) a

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) bsum
  real ( kind = 8 ) c(a)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance

  bsum = sum ( b(1:a) )

  call empirical_discrete_mean ( a, b, c, mean )

  variance = 0.0D+00

  do i = 1, a
    variance = variance + ( b(i) / bsum ) * ( c(i) - mean )**2
  end do

  return
end
subroutine english_sentence_length_cdf ( x, cdf )

!*****************************************************************************80
!
!! ENGLISH_SENTENCE_LENGTH_CDF evaluates the English Sentence Length CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry Kucera, Winthrop Francis,
!    Computational Analysis of Present-Day American English,
!    Brown University Press, 1967.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the sentence length whose CDF is desired.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: sentence_length_max = 79

  real ( kind = 8 ) cdf
  real ( kind = 8 ), dimension ( sentence_length_max ) :: pdf_vec = (/ &
    0.00806D+00, &
    0.01370D+00, &
    0.01862D+00, &
    0.02547D+00, &
    0.03043D+00, &
    0.03189D+00, &
    0.03516D+00, &
    0.03545D+00, &
    0.03286D+00, &
    0.03533D+00, &
    0.03562D+00, &
    0.03788D+00, &
    0.03669D+00, &
    0.03751D+00, &
    0.03518D+00, &
    0.03541D+00, &
    0.03434D+00, &
    0.03305D+00, &
    0.03329D+00, &
    0.03103D+00, &
    0.02867D+00, &
    0.02724D+00, &
    0.02647D+00, &
    0.02526D+00, &
    0.02086D+00, &
    0.02178D+00, &
    0.02128D+00, &
    0.01801D+00, &
    0.01690D+00, &
    0.01556D+00, &
    0.01512D+00, &
    0.01326D+00, &
    0.01277D+00, &
    0.01062D+00, &
    0.01051D+00, &
    0.00901D+00, &
    0.00838D+00, &
    0.00764D+00, &
    0.00683D+00, &
    0.00589D+00, &
    0.00624D+00, &
    0.00488D+00, &
    0.00477D+00, &
    0.00406D+00, &
    0.00390D+00, &
    0.00350D+00, &
    0.00318D+00, &
    0.00241D+00, &
    0.00224D+00, &
    0.00220D+00, &
    0.00262D+00, &
    0.00207D+00, &
    0.00174D+00, &
    0.00174D+00, &
    0.00128D+00, &
    0.00121D+00, &
    0.00103D+00, &
    0.00117D+00, &
    0.00124D+00, &
    0.00082D+00, &
    0.00088D+00, &
    0.00061D+00, &
    0.00061D+00, &
    0.00075D+00, &
    0.00063D+00, &
    0.00056D+00, &
    0.00052D+00, &
    0.00057D+00, &
    0.00031D+00, &
    0.00029D+00, &
    0.00021D+00, &
    0.00017D+00, &
    0.00021D+00, &
    0.00034D+00, &
    0.00031D+00, &
    0.00011D+00, &
    0.00011D+00, &
    0.00008D+00, &
    0.00006D+00 /)
  real, parameter :: pdf_sum = 0.99768D+00
  integer ( kind = 4 ) x

  if ( x < 1 ) then
    cdf = 0.0D+00
  else if ( x < sentence_length_max ) then
    cdf = sum ( pdf_vec(1:x) ) / pdf_sum
  else if ( sentence_length_max <= x ) then
    cdf = 1.0D+00
  end if

  return
end
subroutine english_sentence_length_cdf_inv ( cdf, x )

!*****************************************************************************80
!
!! ENGLISH_SENTENCE_LENGTH_CDF_INV inverts the English Sentence Length CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry Kucera, Winthrop Francis,
!    Computational Analysis of Present-Day American English,
!    Brown University Press, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Output, integer ( kind = 4 ) X, the corresponding sentence length for which
!    CDF(X-1) < CDF <= CDF(X)
!
  implicit none

  integer ( kind = 4 ), parameter :: sentence_length_max = 79

  real ( kind = 8 ) cdf
  real ( kind = 8 ) cum
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( sentence_length_max ) :: pdf_vec = (/ &
    0.00806D+00, &
    0.01370D+00, &
    0.01862D+00, &
    0.02547D+00, &
    0.03043D+00, &
    0.03189D+00, &
    0.03516D+00, &
    0.03545D+00, &
    0.03286D+00, &
    0.03533D+00, &
    0.03562D+00, &
    0.03788D+00, &
    0.03669D+00, &
    0.03751D+00, &
    0.03518D+00, &
    0.03541D+00, &
    0.03434D+00, &
    0.03305D+00, &
    0.03329D+00, &
    0.03103D+00, &
    0.02867D+00, &
    0.02724D+00, &
    0.02647D+00, &
    0.02526D+00, &
    0.02086D+00, &
    0.02178D+00, &
    0.02128D+00, &
    0.01801D+00, &
    0.01690D+00, &
    0.01556D+00, &
    0.01512D+00, &
    0.01326D+00, &
    0.01277D+00, &
    0.01062D+00, &
    0.01051D+00, &
    0.00901D+00, &
    0.00838D+00, &
    0.00764D+00, &
    0.00683D+00, &
    0.00589D+00, &
    0.00624D+00, &
    0.00488D+00, &
    0.00477D+00, &
    0.00406D+00, &
    0.00390D+00, &
    0.00350D+00, &
    0.00318D+00, &
    0.00241D+00, &
    0.00224D+00, &
    0.00220D+00, &
    0.00262D+00, &
    0.00207D+00, &
    0.00174D+00, &
    0.00174D+00, &
    0.00128D+00, &
    0.00121D+00, &
    0.00103D+00, &
    0.00117D+00, &
    0.00124D+00, &
    0.00082D+00, &
    0.00088D+00, &
    0.00061D+00, &
    0.00061D+00, &
    0.00075D+00, &
    0.00063D+00, &
    0.00056D+00, &
    0.00052D+00, &
    0.00057D+00, &
    0.00031D+00, &
    0.00029D+00, &
    0.00021D+00, &
    0.00017D+00, &
    0.00021D+00, &
    0.00034D+00, &
    0.00031D+00, &
    0.00011D+00, &
    0.00011D+00, &
    0.00008D+00, &
    0.00006D+00 /)
  real, parameter :: pdf_sum = 0.99768D+00
  integer ( kind = 4 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ENGLISH_SENTENCE_LENGTH_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  cum = 0.0D+00

  do j = 1, sentence_length_max

    cum = cum + pdf_vec(j)

    if ( cdf <= cum / pdf_sum ) then
      x = j
      return
    end if

  end do

  x = sentence_length_max

  return
end
subroutine english_sentence_length_mean ( mean )

!*****************************************************************************80
!
!! ENGLISH_SENTENCE_LENGTH_MEAN: mean of the English Sentence Length PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry Kucera, Winthrop Francis,
!    Computational Analysis of Present-Day American English,
!    Brown University Press, 1967.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: sentence_length_max = 79

  integer ( kind = 4 ) j
  real ( kind = 8 ) mean
  real ( kind = 8 ), dimension ( sentence_length_max ) :: pdf_vec = (/ &
    0.00806D+00, &
    0.01370D+00, &
    0.01862D+00, &
    0.02547D+00, &
    0.03043D+00, &
    0.03189D+00, &
    0.03516D+00, &
    0.03545D+00, &
    0.03286D+00, &
    0.03533D+00, &
    0.03562D+00, &
    0.03788D+00, &
    0.03669D+00, &
    0.03751D+00, &
    0.03518D+00, &
    0.03541D+00, &
    0.03434D+00, &
    0.03305D+00, &
    0.03329D+00, &
    0.03103D+00, &
    0.02867D+00, &
    0.02724D+00, &
    0.02647D+00, &
    0.02526D+00, &
    0.02086D+00, &
    0.02178D+00, &
    0.02128D+00, &
    0.01801D+00, &
    0.01690D+00, &
    0.01556D+00, &
    0.01512D+00, &
    0.01326D+00, &
    0.01277D+00, &
    0.01062D+00, &
    0.01051D+00, &
    0.00901D+00, &
    0.00838D+00, &
    0.00764D+00, &
    0.00683D+00, &
    0.00589D+00, &
    0.00624D+00, &
    0.00488D+00, &
    0.00477D+00, &
    0.00406D+00, &
    0.00390D+00, &
    0.00350D+00, &
    0.00318D+00, &
    0.00241D+00, &
    0.00224D+00, &
    0.00220D+00, &
    0.00262D+00, &
    0.00207D+00, &
    0.00174D+00, &
    0.00174D+00, &
    0.00128D+00, &
    0.00121D+00, &
    0.00103D+00, &
    0.00117D+00, &
    0.00124D+00, &
    0.00082D+00, &
    0.00088D+00, &
    0.00061D+00, &
    0.00061D+00, &
    0.00075D+00, &
    0.00063D+00, &
    0.00056D+00, &
    0.00052D+00, &
    0.00057D+00, &
    0.00031D+00, &
    0.00029D+00, &
    0.00021D+00, &
    0.00017D+00, &
    0.00021D+00, &
    0.00034D+00, &
    0.00031D+00, &
    0.00011D+00, &
    0.00011D+00, &
    0.00008D+00, &
    0.00006D+00 /)
  real, parameter :: pdf_sum = 0.99768D+00

  mean = 0.0D+00
  do j = 1, sentence_length_max
    mean = mean + real ( j, kind = 8 ) * pdf_vec(j)
  end do

  mean = mean / pdf_sum

  return
end
subroutine english_sentence_length_pdf ( x, pdf )

!*****************************************************************************80
!
!! ENGLISH_SENTENCE_LENGTH_PDF evaluates the English Sentence Length PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = B(X) if 1 <= X <= A
!                = 0    otherwise
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry Kucera, Winthrop Francis,
!    Computational Analysis of Present-Day American English,
!    Brown University Press, 1967.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the sentence length whose probability 
!    is desired.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: sentence_length_max = 79

  real ( kind = 8 ) pdf
  real ( kind = 8 ), dimension ( sentence_length_max ) :: pdf_vec = (/ &
    0.00806D+00, &
    0.01370D+00, &
    0.01862D+00, &
    0.02547D+00, &
    0.03043D+00, &
    0.03189D+00, &
    0.03516D+00, &
    0.03545D+00, &
    0.03286D+00, &
    0.03533D+00, &
    0.03562D+00, &
    0.03788D+00, &
    0.03669D+00, &
    0.03751D+00, &
    0.03518D+00, &
    0.03541D+00, &
    0.03434D+00, &
    0.03305D+00, &
    0.03329D+00, &
    0.03103D+00, &
    0.02867D+00, &
    0.02724D+00, &
    0.02647D+00, &
    0.02526D+00, &
    0.02086D+00, &
    0.02178D+00, &
    0.02128D+00, &
    0.01801D+00, &
    0.01690D+00, &
    0.01556D+00, &
    0.01512D+00, &
    0.01326D+00, &
    0.01277D+00, &
    0.01062D+00, &
    0.01051D+00, &
    0.00901D+00, &
    0.00838D+00, &
    0.00764D+00, &
    0.00683D+00, &
    0.00589D+00, &
    0.00624D+00, &
    0.00488D+00, &
    0.00477D+00, &
    0.00406D+00, &
    0.00390D+00, &
    0.00350D+00, &
    0.00318D+00, &
    0.00241D+00, &
    0.00224D+00, &
    0.00220D+00, &
    0.00262D+00, &
    0.00207D+00, &
    0.00174D+00, &
    0.00174D+00, &
    0.00128D+00, &
    0.00121D+00, &
    0.00103D+00, &
    0.00117D+00, &
    0.00124D+00, &
    0.00082D+00, &
    0.00088D+00, &
    0.00061D+00, &
    0.00061D+00, &
    0.00075D+00, &
    0.00063D+00, &
    0.00056D+00, &
    0.00052D+00, &
    0.00057D+00, &
    0.00031D+00, &
    0.00029D+00, &
    0.00021D+00, &
    0.00017D+00, &
    0.00021D+00, &
    0.00034D+00, &
    0.00031D+00, &
    0.00011D+00, &
    0.00011D+00, &
    0.00008D+00, &
    0.00006D+00 /)
  real, parameter :: pdf_sum = 0.99768D+00
  integer ( kind = 4 ) x

  if ( 1 <= x .and. x <= sentence_length_max ) then
    pdf = pdf_vec(x) / pdf_sum
  else
    pdf = 0.0D+00
  end if

  return
end
subroutine english_sentence_length_sample ( seed, x )

!*****************************************************************************80
!
!! ENGLISH_SENTENCE_LENGTH_SAMPLE samples the English Sentence Length PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry Kucera, Winthrop Francis,
!    Computational Analysis of Present-Day American English,
!    Brown University Press, 1967.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x

  cdf = r8_uniform_01 ( seed )

  call english_sentence_length_cdf_inv ( cdf, x )

  return
end
subroutine english_sentence_length_variance ( variance )

!*****************************************************************************80
!
!! ENGLISH_SENTENCE_LENGTH_VARIANCE: variance of the English Sentence Length PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry Kucera, Winthrop Francis,
!    Computational Analysis of Present-Day American English,
!    Brown University Press, 1967.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: sentence_length_max = 79

  integer ( kind = 4 ) j
  real ( kind = 8 ) mean
  real ( kind = 8 ), dimension ( sentence_length_max ) :: pdf_vec = (/ &
    0.00806D+00, &
    0.01370D+00, &
    0.01862D+00, &
    0.02547D+00, &
    0.03043D+00, &
    0.03189D+00, &
    0.03516D+00, &
    0.03545D+00, &
    0.03286D+00, &
    0.03533D+00, &
    0.03562D+00, &
    0.03788D+00, &
    0.03669D+00, &
    0.03751D+00, &
    0.03518D+00, &
    0.03541D+00, &
    0.03434D+00, &
    0.03305D+00, &
    0.03329D+00, &
    0.03103D+00, &
    0.02867D+00, &
    0.02724D+00, &
    0.02647D+00, &
    0.02526D+00, &
    0.02086D+00, &
    0.02178D+00, &
    0.02128D+00, &
    0.01801D+00, &
    0.01690D+00, &
    0.01556D+00, &
    0.01512D+00, &
    0.01326D+00, &
    0.01277D+00, &
    0.01062D+00, &
    0.01051D+00, &
    0.00901D+00, &
    0.00838D+00, &
    0.00764D+00, &
    0.00683D+00, &
    0.00589D+00, &
    0.00624D+00, &
    0.00488D+00, &
    0.00477D+00, &
    0.00406D+00, &
    0.00390D+00, &
    0.00350D+00, &
    0.00318D+00, &
    0.00241D+00, &
    0.00224D+00, &
    0.00220D+00, &
    0.00262D+00, &
    0.00207D+00, &
    0.00174D+00, &
    0.00174D+00, &
    0.00128D+00, &
    0.00121D+00, &
    0.00103D+00, &
    0.00117D+00, &
    0.00124D+00, &
    0.00082D+00, &
    0.00088D+00, &
    0.00061D+00, &
    0.00061D+00, &
    0.00075D+00, &
    0.00063D+00, &
    0.00056D+00, &
    0.00052D+00, &
    0.00057D+00, &
    0.00031D+00, &
    0.00029D+00, &
    0.00021D+00, &
    0.00017D+00, &
    0.00021D+00, &
    0.00034D+00, &
    0.00031D+00, &
    0.00011D+00, &
    0.00011D+00, &
    0.00008D+00, &
    0.00006D+00 /)
  real, parameter :: pdf_sum = 0.99768D+00
  real ( kind = 8 ) variance

  mean = 0.0D+00
  do j = 1, sentence_length_max
    mean = mean + real ( j, kind = 8 ) * pdf_vec(j)
  end do

  mean = mean / pdf_sum

  variance = 0.0D+00
  do j = 1, sentence_length_max
    variance = variance + pdf_vec(j) * ( j - mean )**2
  end do

  variance = variance / pdf_sum

  return
end
subroutine english_word_length_cdf ( x, cdf )

!*****************************************************************************80
!
!! ENGLISH_WORD_LENGTH_CDF evaluates the English Word Length CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry Kucera, Winthrop Francis,
!    Computational Analysis of Present-Day American English,
!    Brown University Press, 1967.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the word length whose CDF is desired.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: word_length_max = 27

  real ( kind = 8 ) cdf
  real ( kind = 8 ), dimension ( word_length_max ) :: pdf_vec = (/ &
    0.03160D+00, &
    0.16975D+00, &
    0.21192D+00, &
    0.15678D+00, &
    0.10852D+00, &
    0.08524D+00, &
    0.07724D+00, &
    0.05623D+00, &
    0.04032D+00, &
    0.02766D+00, &
    0.01582D+00, &
    0.00917D+00, &
    0.00483D+00, &
    0.00262D+00, &
    0.00099D+00, &
    0.00050D+00, &
    0.00027D+00, &
    0.00022D+00, &
    0.00011D+00, &
    0.00006D+00, &
    0.00005D+00, &
    0.00002D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00 /)
  real, parameter :: pdf_sum = 0.99997D+00
  integer ( kind = 4 ) x

  if ( x < 1 ) then
    cdf = 0.0D+00
  else if ( x < word_length_max ) then
    cdf = sum ( pdf_vec(1:x) ) / pdf_sum
  else if ( word_length_max <= x ) then
    cdf = 1.0D+00
  end if

  return
end
subroutine english_word_length_cdf_inv ( cdf, x )

!*****************************************************************************80
!
!! ENGLISH_WORD_LENGTH_CDF_INV inverts the English Word Length CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry Kucera, Winthrop Francis,
!    Computational Analysis of Present-Day American English,
!    Brown University Press, 1967.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Output, integer ( kind = 4 ) X, the corresponding word length for which
!    CDF(X-1) < CDF <= CDF(X)
!
  implicit none

  integer ( kind = 4 ), parameter :: word_length_max = 27

  real ( kind = 8 ) cdf
  real ( kind = 8 ) cum
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( word_length_max ) :: pdf_vec = (/ &
    0.03160D+00, &
    0.16975D+00, &
    0.21192D+00, &
    0.15678D+00, &
    0.10852D+00, &
    0.08524D+00, &
    0.07724D+00, &
    0.05623D+00, &
    0.04032D+00, &
    0.02766D+00, &
    0.01582D+00, &
    0.00917D+00, &
    0.00483D+00, &
    0.00262D+00, &
    0.00099D+00, &
    0.00050D+00, &
    0.00027D+00, &
    0.00022D+00, &
    0.00011D+00, &
    0.00006D+00, &
    0.00005D+00, &
    0.00002D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00 /)
  real, parameter :: pdf_sum = 0.99997D+00
  integer ( kind = 4 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ENGLISH_WORD_LENGTH_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  cum = 0.0D+00

  do j = 1, word_length_max

    cum = cum + pdf_vec(j)

    if ( cdf <= cum / pdf_sum ) then
      x = j
      return
    end if

  end do

  x = word_length_max

  return
end
subroutine english_word_length_mean ( mean )

!*****************************************************************************80
!
!! ENGLISH_WORD_LENGTH_MEAN evaluates the mean of the English Word Length PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry Kucera, Winthrop Francis,
!    Computational Analysis of Present-Day American English,
!    Brown University Press, 1967.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: word_length_max = 27

  integer ( kind = 4 ) j
  real ( kind = 8 ) mean
  real ( kind = 8 ), dimension ( word_length_max ) :: pdf_vec = (/ &
    0.03160D+00, &
    0.16975D+00, &
    0.21192D+00, &
    0.15678D+00, &
    0.10852D+00, &
    0.08524D+00, &
    0.07724D+00, &
    0.05623D+00, &
    0.04032D+00, &
    0.02766D+00, &
    0.01582D+00, &
    0.00917D+00, &
    0.00483D+00, &
    0.00262D+00, &
    0.00099D+00, &
    0.00050D+00, &
    0.00027D+00, &
    0.00022D+00, &
    0.00011D+00, &
    0.00006D+00, &
    0.00005D+00, &
    0.00002D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00 /)
  real, parameter :: pdf_sum = 0.99997D+00

  mean = 0.0D+00
  do j = 1, word_length_max
    mean = mean + real ( j, kind = 8 ) * pdf_vec(j)
  end do

  mean = mean / pdf_sum

  return
end
subroutine english_word_length_pdf ( x, pdf )

!*****************************************************************************80
!
!! ENGLISH_WORD_LENGTH_PDF evaluates the English Word Length PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = B(X) if 1 <= X <= A
!                = 0    otherwise
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry Kucera, Winthrop Francis,
!    Computational Analysis of Present-Day American English,
!    Brown University Press, 1967.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the word length whose probability is desired.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: word_length_max = 27

  real ( kind = 8 ) pdf
  real ( kind = 8 ), dimension ( word_length_max ) :: pdf_vec = (/ &
    0.03160D+00, &
    0.16975D+00, &
    0.21192D+00, &
    0.15678D+00, &
    0.10852D+00, &
    0.08524D+00, &
    0.07724D+00, &
    0.05623D+00, &
    0.04032D+00, &
    0.02766D+00, &
    0.01582D+00, &
    0.00917D+00, &
    0.00483D+00, &
    0.00262D+00, &
    0.00099D+00, &
    0.00050D+00, &
    0.00027D+00, &
    0.00022D+00, &
    0.00011D+00, &
    0.00006D+00, &
    0.00005D+00, &
    0.00002D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00 /)
  real, parameter :: pdf_sum = 0.99997D+00
  integer ( kind = 4 ) x

  if ( 1 <= x .and. x <= word_length_max ) then
    pdf = pdf_vec(x) / pdf_sum
  else
    pdf = 0.0D+00
  end if

  return
end
subroutine english_word_length_sample ( seed, x )

!*****************************************************************************80
!
!! ENGLISH_WORD_LENGTH_SAMPLE samples the English Word Length PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry Kucera, Winthrop Francis,
!    Computational Analysis of Present-Day American English,
!    Brown University Press, 1967.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x

  cdf = r8_uniform_01 ( seed )

  call english_word_length_cdf_inv ( cdf, x )

  return
end
subroutine english_word_length_variance ( variance )

!*****************************************************************************80
!
!! ENGLISH_WORD_LENGTH_VARIANCE: variance of the English Word Length PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry Kucera, Winthrop Francis,
!    Computational Analysis of Present-Day American English,
!    Brown University Press, 1967.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: word_length_max = 27

  integer ( kind = 4 ) j
  real ( kind = 8 ) mean
  real ( kind = 8 ), dimension ( word_length_max ) :: pdf_vec = (/ &
    0.03160D+00, &
    0.16975D+00, &
    0.21192D+00, &
    0.15678D+00, &
    0.10852D+00, &
    0.08524D+00, &
    0.07724D+00, &
    0.05623D+00, &
    0.04032D+00, &
    0.02766D+00, &
    0.01582D+00, &
    0.00917D+00, &
    0.00483D+00, &
    0.00262D+00, &
    0.00099D+00, &
    0.00050D+00, &
    0.00027D+00, &
    0.00022D+00, &
    0.00011D+00, &
    0.00006D+00, &
    0.00005D+00, &
    0.00002D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00, &
    0.00001D+00 /)
  real, parameter :: pdf_sum = 0.99997D+00
  real ( kind = 8 ) variance

  mean = 0.0D+00
  do j = 1, word_length_max
    mean = mean + real ( j, kind = 8 ) * pdf_vec(j)
  end do

  mean = mean / pdf_sum

  variance = 0.0D+00
  do j = 1, word_length_max
    variance = variance + pdf_vec(j) * ( j - mean )**2
  end do

  variance = variance / pdf_sum

  return
end
function error_f ( x )

!*****************************************************************************80
!
!! ERROR_F evaluates the error function ERF.
!
!  Discussion:
!
!    Since some compilers already supply a routine named ERF which evaluates
!    the error function, this routine has been given a distinct, if
!    somewhat unnatural, name.
!
!    The function is defined by:
!
!      ERF(X) = ( 2 / sqrt ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( - T^2 ) dT.
!
!    Properties of the function include:
!
!      Limit ( X -> -Infinity ) ERF(X) =          -1.0;
!                               ERF(0) =           0.0;
!                               ERF(0.476936...) = 0.5;
!      Limit ( X -> +Infinity ) ERF(X) =          +1.0.
!
!      0.5D+00 * ( ERF(X/sqrt(2)) + 1 ) = Normal_01_CDF(X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    Rational Chebyshev Approximations for the Error Function,
!    Mathematics of Computation,
!    1969, pages 631-638.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the error function.
!
!    Output, real ( kind = 8 ) ERROR_F, the value of the error function.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 5 ) :: a = (/ &
    3.16112374387056560D+00, &
    1.13864154151050156D+02, &
    3.77485237685302021D+02, &
    3.20937758913846947D+03, &
    1.85777706184603153D-01 /)
  real ( kind = 8 ), parameter, dimension ( 4 ) :: b = (/ &
    2.36012909523441209D+01, &
    2.44024637934444173D+02, &
    1.28261652607737228D+03, &
    2.84423683343917062D+03 /)
  real ( kind = 8 ), parameter, dimension ( 9 ) :: c = (/ &
    5.64188496988670089D-01, &
    8.88314979438837594D+00, &
    6.61191906371416295D+01, &
    2.98635138197400131D+02, &
    8.81952221241769090D+02, &
    1.71204761263407058D+03, &
    2.05107837782607147D+03, &
    1.23033935479799725D+03, &
    2.15311535474403846D-08 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    1.57449261107098347D+01, &
    1.17693950891312499D+02, &
    5.37181101862009858D+02, &
    1.62138957456669019D+03, &
    3.29079923573345963D+03, &
    4.36261909014324716D+03, &
    3.43936767414372164D+03, &
    1.23033935480374942D+03 /)
  real ( kind = 8 ) del
  real ( kind = 8 ) error_f
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter, dimension ( 6 ) :: p = (/ &
    3.05326634961232344D-01, &
    3.60344899949804439D-01, &
    1.25781726111229246D-01, &
    1.60837851487422766D-02, &
    6.58749161529837803D-04, &
    1.63153871373020978D-02 /)
  real ( kind = 8 ), parameter, dimension ( 5 ) :: q = (/ &
    2.56852019228982242D+00, &
    1.87295284992346047D+00, &
    5.27905102951428412D-01, &
    6.05183413124413191D-02, &
    2.33520497626869185D-03 /)
  real ( kind = 8 ), parameter :: sqrpi = 0.56418958354775628695D+00
  real ( kind = 8 ), parameter :: thresh = 0.46875D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xabs
  real ( kind = 8 ), parameter :: xbig = 26.543D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ) xnum
  real ( kind = 8 ) xsq

  xabs = abs ( ( x ) )
!
!  Evaluate ERF(X) for |X| <= 0.46875.
!
  if ( xabs <= thresh ) then

    if ( epsilon ( xabs ) < xabs ) then
      xsq = xabs * xabs
    else
      xsq = 0.0D+00
    end if

    xnum = a(5) * xsq
    xden = xsq
    do i = 1, 3
      xnum = ( xnum + a(i) ) * xsq
      xden = ( xden + b(i) ) * xsq
    end do

    error_f = x * ( xnum + a(4) ) / ( xden + b(4) )
!
!  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
!
  else if ( xabs <= 4.0D+00 ) then

    xnum = c(9) * xabs
    xden = xabs
    do i = 1, 7
      xnum = ( xnum + c(i) ) * xabs
      xden = ( xden + d(i) ) * xabs
    end do

    error_f = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = real ( int ( xabs * 16.0D+00 ), kind = 8 ) / 16.0D+00
    del = ( xabs - xsq ) * ( xabs + xsq )
    error_f = exp ( - xsq * xsq ) * exp ( - del ) * error_f

    error_f = ( 0.5D+00 - error_f ) + 0.5D+00

    if ( x < 0.0D+00 ) then
      error_f = - error_f
    end if
!
!  Evaluate ERFC(X) for 4.0D+00 < |X|.
!
  else

    if ( xbig <= xabs ) then

      if ( 0.0D+00 < x ) then
        error_f = 1.0D+00
      else
        error_f = - 1.0D+00
      end if

    else

      xsq = 1.0D+00 / ( xabs * xabs )

      xnum = p(6) * xsq
      xden = xsq
      do i = 1, 4
        xnum = ( xnum + p(i) ) * xsq
        xden = ( xden + q(i) ) * xsq
      end do

      error_f = xsq * ( xnum + p(5) ) / ( xden + q(5) )
      error_f = ( sqrpi - error_f ) / xabs
      xsq = real ( int ( xabs * 16.0D+00 ), kind = 8 ) / 16.0D+00
      del = ( xabs - xsq ) * ( xabs + xsq )
      error_f = exp ( - xsq * xsq ) * exp ( - del ) * error_f

      error_f = ( 0.5D+00 - error_f ) + 0.5D+00

      if ( x < 0.0D+00 ) then
        error_f = - error_f
      end if

    end if

  end if

  return
end
function error_f_inverse ( y )

!*****************************************************************************80
!
!! ERROR_F_INVERSE inverts the error function ERF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Y, the value of the error function.
!
!    Output, real ( kind = 8 ) ERROR_F_INVERSE, the value X such that
!    ERROR_F(X) = Y.
!
  implicit none

  real ( kind = 8 ) error_f_inverse
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = ( y + 1.0D+00 ) / 2.0D+00

  call normal_01_cdf_inv ( z, x )

  error_f_inverse = x / sqrt ( 2.0D+00 )

  return
end
subroutine erlang_cdf ( x, a, b, c, cdf )

!*****************************************************************************80
!
!! ERLANG_CDF evaluates the Erlang CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, integer C, the parameters of the PDF.
!    0.0D+00 < B.
!    0 < C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) gamma_inc
  real ( kind = 8 ) p2
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  if ( x < a ) then

    cdf = 0.0D+00

  else

    x2 = ( x - a ) / b
    p2 = real ( c, kind = 8 )

    cdf = gamma_inc ( p2, x2 )

  end if

  return
end
subroutine erlang_cdf_inv ( cdf, a, b, c, x )

!*****************************************************************************80
!
!! ERLANG_CDF_INV inverts the Erlang CDF.
!
!  Discussion:
!
!    A simple bisection method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, B, integer C, the parameters of the PDF.
!    0.0D+00 < B.
!    0 < C.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf1
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) cdf3
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  real ( kind = 8 ), parameter :: tol = 0.0001D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ERLANG_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = a
    return
  else if ( 1.0D+00 == cdf ) then
    x = huge ( x )
    return
  end if

  x1 = a
  cdf1 = 0.0D+00

  x2 = a + 1.0D+00

  do

    call erlang_cdf ( x2, a, b, c, cdf2 )

    if ( cdf < cdf2 ) then
      exit
    end if

    x2 = a + 2.0D+00 * ( x2 - a )

  end do
!
!  Now use bisection.
!
  it = 0

  do

    it = it + 1

    x3 = 0.5D+00 * ( x1 + x2 )
    call erlang_cdf ( x3, a, b, c, cdf3 )

    if ( abs ( cdf3 - cdf ) < tol ) then
      x = x3
      exit
    end if

    if ( it_max < it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ERLANG_CDF_INV - Fatal error!'
      write ( *, '(a)' ) '  Iteration limit exceeded.'
      return
    end if

    if ( sign ( 1.0D+00, cdf3 - cdf ) == sign ( 1.0D+00, cdf1 - cdf ) ) then
      x1 = x3
      cdf1 = cdf3
    else
      x2 = x3
      cdf2 = cdf3
    end if

  end do

  return
end
function erlang_check ( a, b, c )

!*****************************************************************************80
!
!! ERLANG_CHECK checks the parameters of the Erlang PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, integer C, the parameters of the PDF.
!    0.0D+00 < B.
!    0 < C.
!
!    Output, logical ERLANG_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  logical erlang_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ERLANG_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.0'
    erlang_check = .false.
    return
  end if

  if ( c <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ERLANG_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C <= 0.'
    erlang_check = .false.
    return
  end if

  erlang_check = .true.

  return
end
subroutine erlang_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! ERLANG_MEAN returns the mean of the Erlang PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, integer C, the parameters of the PDF.
!    0.0D+00 < B.
!    0 < C.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) mean

  mean =  a + b * real ( c, kind = 8 )

  return
end
subroutine erlang_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! ERLANG_PDF evaluates the Erlang PDF.
!
!  Discussion:
!
!    PDF(A,B,C;X) = ( ( X - A ) / B )**( C - 1 )
!      / ( B * Gamma ( C ) * EXP ( ( X - A ) / B ) )
!
!    for 0 < B, 0 < C integer, A <= X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, integer C, the parameters of the PDF.
!    0.0D+00 < B.
!    0 < C.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) i4_factorial
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x <= a ) then

    pdf = 0.0D+00

  else

    y = ( x - a ) / b

    pdf = y**( c - 1 ) / ( b * i4_factorial ( c - 1 ) * exp ( y ) )

  end if

  return
end
subroutine erlang_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! ERLANG_SAMPLE samples the Erlang PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, integer C, the parameters of the PDF.
!    0.0D+00 < B.
!    0 < C.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  a2 = 0.0D+00
  b2 = b
  x = a
  do i = 1, c
    call exponential_sample ( a2, b2, seed, x2 )
    x = x + x2
  end do

  return
end
subroutine erlang_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! ERLANG_VARIANCE returns the variance of the Erlang PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, integer C, the parameters of the PDF.
!    0.0D+00 < B.
!    0 < C.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) variance

  variance =  b * b * real ( c )

  return
end
function euler_constant ( )

!*****************************************************************************80
!
!! EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
!
!  Discussion:
!
!    The Euler-Mascheroni constant is often denoted by a lower-case
!    Gamma.  Gamma is defined as
!
!      Gamma = limit ( M -> Infinity )
!        ( Sum ( 1 <= N <= M ) 1 / N ) - Log ( M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EULER_CONSTANT, the value of the
!    Euler-Mascheroni constant.
!
  implicit none

  real ( kind = 8 ) euler_constant

  euler_constant = 0.577215664901532860606512090082402431042D+00

  return
end
subroutine exponential_01_cdf ( x, cdf )

!*****************************************************************************80
!
!! EXPONENTIAL_01_CDF evaluates the Exponential 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then
    cdf = 0.0D+00
  else
    cdf = 1.0D+00 - exp ( - x )
  end if

  return
end
subroutine exponential_01_cdf_inv ( cdf, x )

!*****************************************************************************80
!
!! EXPONENTIAL_01_CDF_INV inverts the Exponential 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXPONENTIAL_01_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = - log ( 1.0D+00 - cdf )

  return
end
subroutine exponential_01_mean ( mean )

!*****************************************************************************80
!
!! EXPONENTIAL_01_MEAN returns the mean of the Exponential 01 PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) mean

  mean = 1.0D+00

  return
end
subroutine exponential_01_pdf ( x, pdf )

!*****************************************************************************80
!
!! EXPONENTIAL_01_PDF evaluates the Exponential 01 PDF.
!
!  Discussion:
!
!    PDF(X) = EXP ( - X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    pdf = 0.0D+00
  else
    pdf = exp ( - x )
  end if

  return
end
subroutine exponential_01_sample ( seed, x )

!*****************************************************************************80
!
!! EXPONENTIAL_01_SAMPLE samples the Exponential PDF with parameter 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2003
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
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  x = - log ( 1.0D+00 - cdf )

  return
end
subroutine exponential_01_variance ( variance )

!*****************************************************************************80
!
!! EXPONENTIAL_01_VARIANCE returns the variance of the Exponential 01 PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) variance

  variance = 1.0D+00

  return
end
subroutine exponential_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! EXPONENTIAL_CDF evaluates the Exponential CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameter of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x <= a ) then
    cdf = 0.0D+00
  else
    cdf = 1.0D+00 - exp ( ( a - x ) / b )
  end if

  return
end
subroutine exponential_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! EXPONENTIAL_CDF_INV inverts the Exponential CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXPONENTIAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a - b * log ( 1.0D+00 - cdf )

  return
end
subroutine exponential_cdf_values ( n_data, lambda, x, fx )

!*****************************************************************************80
!
!! EXPONENTIAL_CDF_VALUES returns some values of the Exponential CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = ExponentialDistribution [ lambda ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) LAMBDA, the parameter of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 9

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.3934693402873666D+00, &
    0.6321205588285577D+00, &
    0.7768698398515702D+00, &
    0.8646647167633873D+00, &
    0.8646647167633873D+00, &
    0.9816843611112658D+00, &
    0.9975212478233336D+00, &
    0.9996645373720975D+00, &
    0.9999546000702375D+00 /)
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_max ) :: lambda_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    lambda = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    lambda = lambda_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function exponential_check ( a, b )

!*****************************************************************************80
!
!! EXPONENTIAL_CHECK checks the parameters of the Exponential CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameter of the PDF.
!    0.0D+00 < B.
!
!    Output, logical EXPONENTIAL_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical exponential_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXPONENTIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.0'
    exponential_check = .false.
    return
  end if

  exponential_check = .true.

  return
end
subroutine exponential_mean ( a, b, mean )

!*****************************************************************************80
!
!! EXPONENTIAL_MEAN returns the mean of the Exponential PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a + b

  return
end
subroutine exponential_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! EXPONENTIAL_PDF evaluates the Exponential PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = ( 1 / B ) * EXP ( ( A - X ) / B )
!
!    The time interval between two Poisson events is a random
!    variable with the Exponential PDF.  The parameter B is the
!    average interval between events.
!
!    In another context, the Exponential PDF is related to
!    the Boltzmann distribution, which describes the relative
!    probability of finding a system, which is in thermal equilibrium
!    at absolute temperature T, in a given state having energy E.
!    The relative probability is
!
!      Boltzmann_Relative_Probability(E,T) = exp ( - E / ( k * T ) ),
!
!    where k is the Boltzmann constant,
!
!      k = 1.38 * 10**(-23) joules / degree Kelvin
!
!    and normalization requires a determination of the possible
!    energy states of the system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A <= X
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x < a ) then
    pdf = 0.0D+00
  else
    pdf = ( 1.0D+00 / b ) * exp ( ( a - x ) / b )
  end if

  return
end
subroutine exponential_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! EXPONENTIAL_SAMPLE samples the Exponential PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call exponential_cdf_inv ( cdf, a, b, x )

  return
end
subroutine exponential_variance ( a, b, variance )

!*****************************************************************************80
!
!! EXPONENTIAL_VARIANCE returns the variance of the Exponential PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = b * b

  return
end
subroutine extreme_values_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! EXTREME_VALUES_CDF evaluates the Extreme Values CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  y = ( x - a ) / b

  cdf = exp ( - exp ( - y ) )

  return
end
subroutine extreme_values_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! EXTREME_VALUES_CDF_INV inverts the Extreme Values CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXTREME_VALUES_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a - b * log ( - log ( cdf ) )

  return
end
subroutine extreme_values_cdf_values ( n_data, alpha, beta, x, fx )

!*****************************************************************************80
!
!! EXTREME_VALUES_CDF_VALUES returns some values of the Extreme Values CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = ExtremeValuesDistribution [ alpha, beta ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) ALPHA, the first parameter of the distribution.
!
!    Output, real ( kind = 8 ) BETA, the second parameter of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) alpha
  real ( kind = 8 ), save, dimension ( n_max ) :: alpha_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  real ( kind = 8 ) beta
  real ( kind = 8 ), save, dimension ( n_max ) :: beta_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.3678794411714423D+00, &
    0.8734230184931166D+00, &
    0.9818510730616665D+00, &
    0.9975243173927525D+00, &
    0.5452392118926051D+00, &
    0.4884435800065159D+00, &
    0.4589560693076638D+00, &
    0.4409910259429826D+00, &
    0.5452392118926051D+00, &
    0.3678794411714423D+00, &
    0.1922956455479649D+00, &
    0.6598803584531254D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    alpha = 0.0D+00
    beta = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    alpha = alpha_vec(n_data)
    beta = beta_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function extreme_values_check ( a, b )

!*****************************************************************************80
!
!! EXTREME_VALUES_CHECK checks the parameters of the Extreme Values CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, logical EXTREME_VALUES_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical extreme_values_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXTREME_VALUES_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    extreme_values_check = .false.
    return
  end if

  extreme_values_check = .true.

  return
end
subroutine extreme_values_mean ( a, b, mean )

!*****************************************************************************80
!
!! EXTREME_VALUES_MEAN returns the mean of the Extreme Values PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) euler_constant
  real ( kind = 8 ) mean

  mean = a + b * euler_constant ( )

  return
end
subroutine extreme_values_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! EXTREME_VALUES_PDF evaluates the Extreme Values PDF.
!
!  Discussion:
!
!    PDF(A,B;X) =
!      ( 1 / B ) * exp ( ( A - X ) / B ) * exp ( - exp ( ( A - X ) / B  ) ).
!
!    The Extreme Values PDF is also known as the Fisher-Tippet PDF
!    and the Log-Weibull PDF.
!
!    The special case A = 0 and B = 1 is the Gumbel PDF.
!
!    The Extreme Values PDF is the limiting distribution for the
!    smallest or largest value in a large sample drawn from
!    any of a great variety of distributions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein, editor,
!    CRC Concise Encylopedia of Mathematics,
!    CRC Press, 1998.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  pdf = ( 1.0D+00 / b ) * exp ( ( a - x ) / b - exp ( ( a - x ) / b ) )

  return
end
subroutine extreme_values_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! EXTREME_VALUES_SAMPLE samples the Extreme Values PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call extreme_values_cdf_inv ( cdf, a, b, x )

  return
end
subroutine extreme_values_variance ( a, b, variance )

!*****************************************************************************80
!
!! EXTREME_VALUES_VARIANCE returns the variance of the Extreme Values PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) variance

  variance = pi * pi * b * b / 6.0D+00

  return
end
subroutine f_cdf ( x, m, n, cdf )

!*****************************************************************************80
!
!! F_CDF evaluates the F central CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the PDF.
!    1 <= M,
!    1 <= N.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) beta_inc
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then

    cdf = 0.0D+00

  else

    arg1 = 0.5D+00 * real ( n, kind = 8 )
    arg2 = 0.5D+00 * real ( m, kind = 8 )
    arg3 = real ( n, kind = 8 ) &
      / ( real ( n, kind = 8 ) + real ( m, kind = 8 ) * x )

    cdf = 1.0D+00 - beta_inc ( arg1, arg2, arg3 )

  end if

  return
end
subroutine f_cdf_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! F_CDF_VALUES returns some values of the F CDF test function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = FRatioDistribution [ dfn, dfd ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) A, integer B, the parameters of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) a
  integer ( kind = 4 ), save, dimension ( n_max ) :: a_vec = (/ &
    1, &
    1, &
    5, &
    1, &
    2, &
    4, &
    1, &
    6, &
    8, &
    1, &
    3, &
    6, &
    1, &
    1, &
    1, &
    1, &
    2, &
    3, &
    4, &
    5 /)
  integer ( kind = 4 ) b
  integer ( kind = 4 ), save, dimension ( n_max ) :: b_vec = (/ &
     1, &
     5, &
     1, &
     5, &
    10, &
    20, &
     5, &
     6, &
    16, &
     5, &
    10, &
    12, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5, &
     5 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.4999714850534485D+00, &
    0.4996034370170990D+00, &
    0.7496993658293228D+00, &
    0.7504656462757382D+00, &
    0.7514156325324275D+00, &
    0.8999867031372156D+00, &
    0.8997127554259699D+00, &
    0.9002845660853669D+00, &
    0.9500248817817622D+00, &
    0.9500574946122442D+00, &
    0.9501926400000000D+00, &
    0.9750133887312993D+00, &
    0.9900022327445249D+00, &
    0.9949977837872073D+00, &
    0.9989999621122122D+00, &
    0.5687988496283079D+00, &
    0.5351452100063650D+00, &
    0.5143428032407864D+00, &
    0.5000000000000000D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     1.00D+00, &
     0.528D+00, &
     1.89D+00, &
     1.69D+00, &
     1.60D+00, &
     1.47D+00, &
     4.06D+00, &
     3.05D+00, &
     2.09D+00, &
     6.61D+00, &
     3.71D+00, &
     3.00D+00, &
    10.01D+00, &
    16.26D+00, &
    22.78D+00, &
    47.18D+00, &
     1.00D+00, &
     1.00D+00, &
     1.00D+00, &
     1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0
    b = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function f_check ( m, n )

!*****************************************************************************80
!
!! F_CHECK checks the parameters of the F PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the PDF.
!    1 <= M,
!    1 <= N.
!
!    Output, logical F_CHECK, is TRUE if the parameters are legal.
!
  implicit none

  logical f_check
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( m < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F_CHECK - Fatal error!'
    write ( *, '(a)' ) '  M < 1.'
    f_check = .false.
    return
  end if

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    f_check = .false.
    return
  end if

  f_check = .true.

  return
end
subroutine f_mean ( m, n, mean )

!*****************************************************************************80
!
!! F_MEAN returns the mean of the F central PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the PDF.
!    1 <= M,
!    1 <= N.
!    Note, however, that the mean is not defined unless 3 <= N.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  integer ( kind = 4 ) m
  real ( kind = 8 ) mean
  integer ( kind = 4 ) n

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F_MEAN - Fatal error!'
    write ( *, '(a)' ) '  The mean is not defined for N < 3.'
    stop
  end if

  mean = real ( n, kind = 8 ) / real ( n - 2, kind = 8 )

  return
end
subroutine f_pdf ( x, m, n, pdf )

!*****************************************************************************80
!
!! F_PDF evaluates the F central PDF.
!
!  Discussion:
!
!    PDF(M,N;X) = M**(M/2) * X**((M-2)/2)
!      / ( Beta(M/2,N/2) * N**(M/2) * ( 1 + (M/N) * X )**((M+N)/2)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the PDF.
!    1 <= M,
!    1 <= N.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) bot1
  real ( kind = 8 ) bot2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) pdf
  real ( kind = 8 ) top
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then

    pdf = 0.0D+00

  else

    a = real ( m, kind = 8 )
    b = real ( n, kind = 8 )

    top = sqrt ( a**m * b**n * x**( m - 2 ) )
    bot1 = beta ( a / 2.0D+00, b / 2.0D+00 )
    bot2 =  sqrt ( ( b + a * x )**( m + n ) )

    pdf = top / ( bot1 * bot2 )

  end if

  return
end
subroutine f_sample ( m, n, seed, x )

!*****************************************************************************80
!
!! F_SAMPLE samples the F central PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the PDF.
!    1 <= M,
!    1 <= N.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  real ( kind = 8 ) xm
  real ( kind = 8 ) xn

  a = real ( m, kind = 8 )
  call chi_square_sample ( a, seed, xm )

  a = real ( n, kind = 8 )
  call chi_square_sample ( a, seed, xn )

  x = real ( n, kind = 8 ) * xm / ( real ( m, kind = 8 ) * xn )

  return
end
subroutine f_variance ( m, n, variance )

!*****************************************************************************80
!
!! F_VARIANCE returns the variance of the F central PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the PDF.
!    1 <= M,
!    1 <= N.
!    Note, however, that the variance is not defined unless 5 <= N.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) variance

  if ( n < 5 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F_VARIANCE - Fatal error!'
    write ( *, '(a)' ) '  The variance is not defined for N < 5.'
    stop
  end if

  variance = real ( 2 * n * n * ( m + n - 2 ), kind = 8 ) / &
    real ( m * ( n - 2 )**2 * ( n - 4 ), kind = 8 )

  return
end
subroutine f_noncentral_cdf_values ( n_data, n1, n2, lambda, x, fx )

!*****************************************************************************80
!
!! F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NoncentralFRatioDistribution [ n1, n2, lambda ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N1, integer N2, the numerator and denominator
!    degrees of freedom.
!
!    Output, real ( kind = 8 ) LAMBDA, the noncentrality parameter.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 22

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.6367825323508774D+00, &
    0.5840916116305482D+00, &
    0.3234431872392788D+00, &
    0.4501187879813550D+00, &
    0.6078881441188312D+00, &
    0.7059275551414605D+00, &
    0.7721782003263727D+00, &
    0.8191049017635072D+00, &
    0.3170348430749965D+00, &
    0.4327218008454471D+00, &
    0.4502696915707327D+00, &
    0.4261881186594096D+00, &
    0.6753687206341544D+00, &
    0.4229108778879005D+00, &
    0.6927667261228938D+00, &
    0.3632174676491226D+00, &
    0.4210054012695865D+00, &
    0.4266672258818927D+00, &
    0.4464016600524644D+00, &
    0.8445888579504827D+00, &
    0.4339300273343604D+00 /)
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_max ) :: lambda_vec = (/ &
    0.00D+00, &
    0.00D+00, &
    0.25D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    1.00D+00, &
    1.00D+00, &
    0.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) n1
  integer ( kind = 4 ), save, dimension ( n_max ) :: n1_vec = (/ &
     1,  1,  1,  1, &
     1,  1,  1,  1, &
     1,  1,  2,  2, &
     3,  3,  4,  4, &
     5,  5,  6,  6, &
     8, 16 /)
  integer ( kind = 4 ) n2
  integer ( kind = 4 ), save, dimension ( n_max ) :: n2_vec = (/ &
     1,  5,  5,  5, &
     5,  5,  5,  5, &
     5,  5,  5, 10, &
     5,  5,  5,  5, &
     1,  5,  6, 12, &
    16,  8 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    0.50D+00, &
    1.00D+00, &
    2.00D+00, &
    3.00D+00, &
    4.00D+00, &
    5.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    2.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n1 = 0
    n2 = 0
    lambda = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    n1 = n1_vec(n_data)
    n2 = n2_vec(n_data)
    lambda = lambda_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function f_noncentral_check ( a, m, n )

!*****************************************************************************80
!
!! F_NONCENTRAL_CHECK checks the parameters of the F noncentral PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the PDF.
!    1 <= M,
!    1 <= N.
!
!    Output, logical F_NONCENTRAL_CHECK, is TRUE if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  logical f_noncentral_check
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F_NONCENTRAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    f_noncentral_check = .false.
    return
  end if

  if ( m < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F_NONCENTRAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  M < 1.'
    f_noncentral_check = .false.
    return
  end if

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F_NONCENTRAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    f_noncentral_check = .false.
    return
  end if

  f_noncentral_check = .true.

  return
end
subroutine f_noncentral_mean ( a, m, n, mean )

!*****************************************************************************80
!
!! F_NONCENTRAL_MEAN returns the mean of the F noncentral PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!
!    Input, integer ( kind = 4 ) M, N, parameters of the PDF.
!    1 <= M,
!    1 <= N.
!    Note, however, that the mean is not defined unless 3 <= N.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  integer ( kind = 4 ) m
  real ( kind = 8 ) mean
  integer ( kind = 4 ) n

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F_NONCENTRAL_MEAN - Fatal error!'
    write ( *, '(a)' ) '  The mean is not defined for N < 3.'
    stop
  end if

  mean = ( real ( m, kind = 8 ) + a ) * real ( n, kind = 8 ) &
    / ( real ( m, kind = 8 ) * real ( n - 2, kind = 8 ) )

  return
end
subroutine f_noncentral_variance ( a, m, n, variance )

!*****************************************************************************80
!
!! F_NONCENTRAL_VARIANCE returns the variance of the F noncentral PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!
!    Input, integer ( kind = 4 ) M, N, parameters of the PDF.
!    1 <= M,
!    1 <= N.
!    Note, however, that the variance is not defined unless 5 <= N.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  integer ( kind = 4 ) m
  real ( kind = 8 ) mr
  integer ( kind = 4 ) n
  real ( kind = 8 ) nr
  real ( kind = 8 ) variance

  if ( n < 5 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F_NONCENTRAL_VARIANCE - Fatal error!'
    write ( *, '(a)' ) '  The variance is not defined for N < 5.'
    stop
  end if

  mr = real ( m, kind = 8 )
  nr = real ( n, kind = 8 )

  variance = ( ( mr + a )**2 + 2.0D+00 * ( mr + a ) * nr**2 ) &
    / ( ( nr - 2.0D+00 ) * ( nr - 4.0D+00 ) * mr**2 ) - &
    ( mr + a )**2 * nr**2 / ( ( nr - 2.0D+00 )**2 * mr**2 )

  return
end
function factorial_log ( n )

!*****************************************************************************80
!
!! FACTORIAL_LOG returns the logarithm of N factorial.
!
!  Discussion:
!
!    N! = Product ( 1 <= I <= N ) I
!
!    N! = Gamma(N+1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the function.
!    0 <= N.
!
!    Output, real ( kind = 8 ) FACTORIAL_LOG, the logarithm of N!.
!
  implicit none

  real ( kind = 8 ) factorial_log
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACTORIAL_LOG - Fatal error!'
    write ( *, '(a)' ) '  N < 0.'
    stop
  end if

  factorial_log = 0.0D+00

  do i = 2, n
    factorial_log = factorial_log + log ( real ( i, kind = 8 ) )
  end do

  return
end
function factorial_stirling ( n )

!*****************************************************************************80
!
!! FACTORIAL_STIRLING computes Stirling's approximation to N!.
!
!  Discussion:
!
!    This routine returns the raw approximation for all nonnegative
!    values of N.  If N is less than 0, the value is returned as 0,
!    and if N is 0, the value of 1 is returned.  In all other cases,
!    Stirling's formula is used.
!
!    The factorial function N! is defined by
!
!      N! = Product ( 1 <= I <= N ) I
!
!    Stirling's approximation to N! is
!
!      Stirling ( N ) = sqrt ( 2 * PI * N ) * ( N / E )**N * E**(1/(12*N) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, real ( kind = 8 ) FACTORIAL_STIRLING, an approximation to N!.
!
  implicit none

  real ( kind = 8 ), parameter :: e_natural = &
    2.71828182845904523536028747135266249775724709369995
  real ( kind = 8 ) factorial_stirling
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) value

  if ( n < 0 ) then

    value = 0.0D+00

  else if ( n == 0 ) then

    value = 1.0D+00

  else

    value = sqrt ( 2.0D+00 * pi * real ( n, kind = 8 ) ) &
      * ( real ( n, kind = 8 ) / e_natural )**n &
      * exp ( 1.0D+00 / real ( 12 * n, kind = 8 ) )

  end if

  factorial_stirling = value

  return
end
subroutine fermi_dirac_sample ( u, v, seed, z )

!*****************************************************************************80
!
!! FERMI_DIRAC_SAMPLE returns a sample from a (continuous) Fermi-Dirac distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    Original BASIC version by Frederick Ruckdeschel.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Frederick Ruckdeschel,
!    BASIC Scientific Subroutines,
!    Volume I,
!    McGraw Hill, 1980,
!    ISBN: 0-07-054201-5,
!    LC: QA76.95.R82.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U, V, the parameters of the distribution.
!    The value of U represents the halfway point for the distribution.
!    Half the probability is to the left, and half to the right, of
!    the value U.  The value of V controls the shape of the distribution.
!    The ratio U/V determines the relative shape of the distribution.
!    Values of U/V in excess of 100 will risk overflow.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) Z, a sample from the Fermi-Dirac distribution.
!    Output values will be nonnegative, and roughly half of them should
!    be less than or equal to U.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ), parameter :: iter_max = 1000
  integer ( kind = 4 ) iter_num
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) z

  x = r8_uniform_01 ( seed )
  y = 1.0D+00
  a = exp ( 4.0D+00 * u / v )
  b = ( x - 1.0D+00 ) * log ( 1.0D+00 + a )

  iter_num = 0

  do

    y1 = b + log ( a + exp ( y ) )

    if ( abs ( y - y1 ) < 0.001D+00 ) then
      exit
    end if

    y = y1

    iter_num = iter_num + 1

    if ( iter_max < iter_num ) then
      exit
    end if

  end do

  z = v * y1 / 4.0D+00

  return
end
subroutine fisher_pdf ( x, kappa, mu, pdf )

!*****************************************************************************80
!
!! FISHER_PDF evaluates the Fisher PDF.
!
!  Discussion:
!
!    The formula for the PDF is:
!
!      PDF(KAPPA,MU;X) = C(KAPPA) * exp ( KAPPA * MU' * X )
!
!    where:
!
!      0 <= KAPPA is the concentration parameter,
!      MU is a point on the unit sphere, the mean direction,
!      X is any point on the unit sphere,
!      and C(KAPPA) is a normalization factor:
!
!      C(KAPPA) = sqrt ( KAPPA ) / ( ( 2 * pi )**(3/2) * I(0.5,KAPPA) )
!
!    where
!
!      I(nu,X) is the Bessel function of order NU and argument X.
!
!    For a fixed value of MU, the value of KAPPA determines the
!    tendency of sample points to tend to be near MU.  In particular,
!    KAPPA = 0 corresponds to a uniform distribution of points on the
!    sphere, but as KAPPA increases, the sample points will tend to
!    cluster more closely to MU.
!
!    The Fisher distribution for points on the unit sphere is
!    analogous to the normal distribution of points on a line,
!    and, more precisely, to the von Mises distribution on a circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kanti Mardia, Peter Jupp,
!    Directional Statistics,
!    Wiley, 2000,
!    LC: QA276.M335
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(3), the argument of the PDF.
!    X should have unit Euclidean norm, but this routine will
!    automatically work with a normalized version of X.
!
!    Input, real ( kind = 8 ) KAPPA, the concentration parameter.
!    0 <= KAPPA is required.
!
!    Input, real ( kind = 8 ) MU(3), the mean direction.
!    MU should have unit Euclidean norm, but this routine will
!    automatically work with a normalized version of MU.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: nb = 1

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) b(nb)
  real ( kind = 8 ) cf
  integer ( kind = 4 ) ize
  real ( kind = 8 ) kappa
  real ( kind = 8 ) mu(3)
  real ( kind = 8 ) mu_norm
  integer ( kind = 4 ) ncalc
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(3)
  real ( kind = 8 ) x_norm

  if ( kappa < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FISHER_PDF - Fatal error!'
    write ( *, '(a)' ) '  KAPPA must be nonnegative.'
    write ( *, '(a,g14.6)' ) '  Input KAPPA = ', kappa
    stop
  end if

  if ( kappa == 0.0D+00 ) then
    pdf = 1.0D+00 / ( 4.0D+00 * pi )
    return
  end if
!
!  Compute the normalization factor CF.
!
  alpha = 0.5D+00
  ize = 1

  call ribesl ( kappa, alpha, nb, ize, b, ncalc )

  cf = sqrt ( kappa ) / ( sqrt ( ( 2.0D+00 * pi )**3 ) * b(1) )
!
!  Normalize MU.
!
  mu_norm = sqrt ( sum ( mu(1:3)**2 ) )

  if ( mu_norm == 0.0D+00 ) then
    pdf = cf
    return
  end if
!
!  Normalize X.
!
  x_norm = sqrt ( sum ( x(1:3)**2 ) )

  if ( x_norm == 0.0D+00 ) then
    pdf = cf
    return
  end if
!
!  Evaluate the PDF.
!
  arg = kappa * dot_product ( x(1:3), mu(1:3) ) / ( x_norm * mu_norm )

  pdf = cf * exp ( arg )

  return
end
subroutine fisher_sample ( kappa, mu, n, seed, xyz )

!*****************************************************************************80
!
!! FISHER_SAMPLE samples the Fisher distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Nicholas Fisher, Toby Lewis, Brian Embleton,
!    Statistical Analysis of Spherical Data,
!    Cambridge, 2003,
!    ISBN13: 978-0521456999,
!    LC: QA276.F489.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) KAPPA, the concentration parameter.
!
!    Input, real ( kind = 8 ) MU(3), the mean direction.
!    MU should have unit Euclidean norm, but this routine will
!    automatically work with a normalized version of MU.
!
!    Input, integer ( kind = 4 ) N, the number of samples to choose.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) XYZ(3,N), a sample of the Fisher distribution.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) ALPHA, BETA, the colatitude (theta) and
!    longitude (phi) of the mean direction.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,3)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) kappa
  real ( kind = 8 ) lambda
  real ( kind = 8 ) mu(3)
  real ( kind = 8 ) mu_norm
  real ( kind = 8 ) phi(1:n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta(1:n)
  real ( kind = 8 ) xyz(3,n)

  mu_norm = sqrt ( sum ( mu(1:3)**2 ) )

  if ( mu_norm == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FISHER_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  MU(1:3) = 0.'
    stop
  end if

  alpha = - acos ( mu(3) / mu_norm )
  beta = atan2 ( mu(2), mu(1) )

  lambda = exp ( - 2.0D+00 * kappa )

  call r8vec_uniform_01 ( n, seed, theta )

  if ( kappa == 0.0D+00 ) then

    theta(1:n) = 2.0D+00 * asin ( sqrt ( 1.0D+00 - theta(1:n) ) )

  else

    theta(1:n) = 2.0D+00 * asin ( sqrt ( &
      - log ( theta(1:n) * ( 1.0D+00 - lambda ) + lambda ) &
      / ( 2.0D+00 * kappa ) ) )

  end if

  call r8vec_uniform_01 ( n, seed, phi )

  phi(1:n) = 2.0D+00 * pi * phi(1:n)
!
!  Compute the unrotated points.
!
  xyz(1,1:n) = sin ( theta(1:n) ) * cos ( phi(1:n) )
  xyz(2,1:n) = sin ( theta(1:n) ) * sin ( phi(1:n) )
  xyz(3,1:n) = cos ( theta(1:n) )
!
!  Compute the rotation matrix.
!
  a(1,1) =   cos ( alpha ) * cos ( beta )
  a(2,1) =                 - sin ( beta )
  a(3,1) =   sin ( alpha ) * cos ( beta )

  a(1,2) =   cos ( alpha ) * sin ( beta )
  a(2,2) =                 + cos ( beta )
  a(3,2) =   sin ( alpha ) * sin ( beta )

  a(1,3) = - sin ( alpha )
  a(2,3) =   0.0D+00
  a(3,3) =   cos ( alpha )
!
!  Rotate the points.
!
  xyz(1:3,1:n) = matmul ( a(1:3,1:3), xyz(1:3,1:n) )

  return
end
subroutine fisk_cdf ( x, a, b, c, cdf )

!*****************************************************************************80
!
!! FISK_CDF evaluates the Fisk CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x <= a ) then
    cdf = 0.0D+00
  else
    cdf = 1.0D+00 / ( 1.0D+00 + ( b / ( x - a ) )**c )
  end if

  return
end
subroutine fisk_cdf_inv ( cdf, a, b, c, x )

!*****************************************************************************80
!
!! FISK_CDF_INV inverts the Fisk CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FISK_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf <= 0.0D+00 ) then
    x = a
  else if ( cdf < 1.0D+00 ) then
    x = a + b * ( cdf / ( 1.0D+00 - cdf ) )**( 1.0D+00 / c )
  else if ( 1.0D+00 <= cdf ) then
    x = huge ( x )
  end if

  return
end
function fisk_check ( a, b, c )

!*****************************************************************************80
!
!! FISK_CHECK checks the parameters of the Fisk PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, logical FISK_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical fisk_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FISK_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    fisk_check = .false.
    return
  end if

  if ( c <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FISK_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C <= 0.'
    fisk_check = .false.
    return
  end if

  fisk_check = .true.

  return
end
subroutine fisk_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! FISK_MEAN returns the mean of the Fisk PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) csc
  real ( kind = 8 ) mean
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  if ( c <= 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FISK_MEAN - Fatal error!'
    write ( *, '(a)' ) '  No mean defined for C <= 1.0'
    stop
  end if

  mean = a + pi * ( b / c ) * csc ( pi / c )

  return
end
subroutine fisk_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! FISK_PDF evaluates the Fisk PDF.
!
!  Discussion:
!
!    The Fisk PDF is also known as the Log Logistic PDF.
!
!    The formula for the PDF is:
!
!    PDF(A,B,C;X) =
!      ( C / B ) * ( ( X - A ) / B )^( C - 1 ) /
!      ( 1 + ( ( X - A ) / B )^C )^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A <= X
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x <= a ) then

    pdf = 0.0D+00

  else

    y = ( x - a ) / b

    pdf = ( c / b ) * y**( c - 1.0D+00 ) / ( 1.0D+00 + y**c )**2

  end if

  return
end
subroutine fisk_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! FISK_SAMPLE samples the Fisk PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call fisk_cdf_inv ( cdf, a, b, c, x )

  return
end
subroutine fisk_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! FISK_VARIANCE returns the variance of the Fisk PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) csc
  real ( kind = 8 ) g
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) variance

  if ( c <= 2.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FISK_VARIANCE - Fatal error!'
    write ( *, '(a)' ) '  No variance defined for C <= 2.0'
    stop
  end if

  g = pi / c

  variance = b * b &
    * ( 2.0D+00 * g * csc ( 2.0D+00 * g ) - ( g * csc ( g ) )**2 )

  return
end
subroutine folded_normal_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! FOLDED_NORMAL_CDF evaluates the Folded Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!    0.0D+00 <= X.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf1
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  if ( x < 0.0D+00 ) then
    cdf = 0.0D+00
  else
    x1 = ( x - a ) / b
    call normal_01_cdf ( x1, cdf1 )
    x2 = ( - x - a ) / b
    call normal_01_cdf ( x2, cdf2 )
    cdf = cdf1 - cdf2
  end if

  return
end
subroutine folded_normal_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! FOLDED_NORMAL_CDF_INV inverts the Folded Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the argument of the CDF.
!    0.0D+00 <= X.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf1
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) cdf3
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  real ( kind = 8 ), parameter :: tol = 0.0001D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FOLDED_NORMAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = 0.0D+00
    return
  else if ( 1.0D+00 == cdf ) then
    x = huge ( x )
    return
  end if
!
!  Find X1, for which the value of CDF will be too small.
!
  if ( 0.0D+00 <= a ) then
    call normal_cdf_inv ( cdf, a, b, x1 )
  else
    call normal_cdf_inv ( cdf, -a, b, x1 )
  end if
  x1 = max ( x1, 0.0D+00 )
  call folded_normal_cdf ( x1, a, b, cdf1 )
!
!  Find X2, for which the value of CDF will be too big.
!
  cdf2 = ( 1.0D+00 - cdf ) / 2.0D+00

  call normal_cdf_inv ( cdf2, a, b, xa )
  call normal_cdf_inv ( cdf2, -a, b, xb )
  x2 = max ( abs ( xa ), abs ( xb ) )
  call folded_normal_cdf ( x2, a, b, cdf2 )
!
!  Now use bisection.
!
  it = 0

  do

    it = it + 1

    x3 = 0.5D+00 * ( x1 + x2 )
    call folded_normal_cdf ( x3, a, b, cdf3 )

    if ( abs ( cdf3 - cdf ) < tol ) then
      x = x3
      exit
    end if

    if ( it_max < it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FOLDED_NORMAL_CDF_INV - Fatal error!'
      write ( *, '(a)' ) '  Iteration limit exceeded.'
      stop
    end if

    if ( sign ( 1.0D+00, cdf3 - cdf ) == sign ( 1.0D+00, cdf1 - cdf ) ) then
      x1 = x3
      cdf1 = cdf3
    else
      x2 = x3
      cdf2 = cdf3
    end if

  end do

  return
end
function folded_normal_check ( a, b )

!*****************************************************************************80
!
!! FOLDED_NORMAL_CHECK checks the parameters of the Folded Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A,
!    0.0D+00 < B.
!
!    Output, logical FOLDED_NORMAL_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical folded_normal_check

  if ( a < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FOLDED_NORMAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A < 0.'
    folded_normal_check = .false.
    return
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FOLDED_NORMAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    folded_normal_check = .false.
    return
  end if

  folded_normal_check = .true.

  return
end
subroutine folded_normal_mean ( a, b, mean )

!*****************************************************************************80
!
!! FOLDED_NORMAL_MEAN returns the mean of the Folded Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) mean
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  a2 = a / b

  call normal_01_cdf ( a2, cdf )

  mean = b * sqrt ( 2.0D+00 / PI ) * exp ( - 0.5D+00 * a2 * a2 ) &
    - a * ( 1.0D+00 - 2.0D+00 * cdf )

  return
end
subroutine folded_normal_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! FOLDED_NORMAL_PDF evaluates the Folded Normal PDF.
!
!  Discussion:
!
!    The formula for the PDF is:
!
!    PDF(A;X) = sqrt ( 2 / PI ) * ( 1 / B ) * cosh ( A * X / B^2 )
!      * exp ( - 0.5D+00 * ( X^2 + A^2 ) / B^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    pdf = 0.0D+00
  else
    pdf = sqrt ( 2.0D+00 / PI ) * ( 1.0D+00 / b ) * cosh ( a * x / b**2 ) &
      * exp ( - 0.5D+00 * ( x * x + a * a ) / b**2 )
  end if

  return
end
subroutine folded_normal_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! FOLDED_NORMAL_SAMPLE samples the Folded Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A,
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call folded_normal_cdf_inv ( cdf, a, b, x )

  return
end
subroutine folded_normal_variance ( a, b, variance )

!*****************************************************************************80
!
!! FOLDED_NORMAL_VARIANCE returns the variance of the Folded Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance

  call folded_normal_mean ( a, b, mean )

  variance = a * a + b * b - mean * mean

  return
end
subroutine frechet_cdf ( x, alpha, cdf )

!*****************************************************************************80
!
!! FRECHET_CDF evaluates the Frechet CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    It is required that 0.0 < ALPHA.
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( alpha <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FRECHET_CDF - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= 0.0.'
    stop
  end if

  if ( x <= 0.0D+00 ) then
    cdf = 0.0D+00
  else
    cdf = exp ( - 1.0D+00 / x**alpha )
  end if

  return
end
subroutine frechet_cdf_inv ( cdf, alpha, x )

!*****************************************************************************80
!
!! FRECHET_CDF_INV inverts the Frechet CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    It is required that 0.0 < ALPHA.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FRECHET_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( alpha <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FRECHET_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= 0.0.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = 0.0D+00
  else
    x =  ( - 1.0D+00 / log ( cdf ) ) ** ( 1.0D+00 / alpha )
  end if

  return
end
subroutine frechet_mean ( alpha, mean )

!*****************************************************************************80
!
!! FRECHET_MEAN returns the mean of the Frechet PDF.
!
!  Discussion:
!
!    The distribution does not have a mean value unless 1 < ALPHA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    It is required that 1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) mean
  real ( kind = 8 ) r8_gamma

  if ( alpha <= 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FRECHET_MEAN - Fatal error!'
    write ( *, '(a)' ) '  Mean does not exist if ALPHA <= 1.'
    stop
  end if

  mean = r8_gamma ( ( alpha - 1.0D+00 ) / alpha )

  return
end
subroutine frechet_pdf ( x, alpha, pdf )

!*****************************************************************************80
!
!! FRECHET_PDF evaluates the Frechet PDF.
!
!  Discussion:
!
!    PDF(X) = ALPHA * exp ( -1 / X^ALPHA ) / X^(ALPHA+1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    It is required that 0.0 < ALPHA.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( alpha <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FRECHET_PDF - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= 0.0.'
    stop
  end if

  pdf = alpha * exp ( - 1.0D+00 / x**alpha ) / x**( alpha + 1.0D+00 )

  return
end
subroutine frechet_sample ( alpha, seed, x )

!*****************************************************************************80
!
!! FRECHET_SAMPLE samples the Frechet PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    It is required that 0.0 < ALPHA.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  if ( alpha <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FRECHET_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= 0.0.'
    stop
  end if

  cdf = r8_uniform_01 ( seed )

  call frechet_cdf_inv ( cdf, alpha, x )

  return
end
subroutine frechet_variance ( alpha, variance )

!*****************************************************************************80
!
!! FRECHET_VARIANCE returns the variance of the Frechet PDF.
!
!  Discussion:
!
!    The PDF does not have a variance unless 2 < ALPHA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    It is required that 2.0 < ALPHA.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) mean
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) variance

  if ( alpha <= 2.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FRECHET_VARIANCE - Fatal error!'
    write ( *, '(a)' ) '  Variance does not exist if ALPHA <= 2.'
    stop
  end if

  mean = r8_gamma ( ( alpha - 1.0D+00 ) / alpha )

  variance = r8_gamma ( ( alpha - 2.0D+00 ) / alpha ) - mean * mean

  return
end
subroutine gamma_cdf ( x, a, b, c, cdf )

!*****************************************************************************80
!
!! GAMMA_CDF evaluates the Gamma CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A <= X
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) gamma_inc
  real ( kind = 8 ) p2
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  x2 = ( x - a ) / b
  p2 = c

  cdf = gamma_inc ( p2, x2 )

  return
end
subroutine gamma_cdf_values ( n_data, mu, sigma, x, fx )

!*****************************************************************************80
!
!! GAMMA_CDF_VALUES returns some values of the Gamma CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = GammaDistribution [ mu, sigma ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) MU, the mean of the distribution.
!
!    Output, real ( kind = 8 ) SIGMA, the variance of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.8646647167633873D+00, &
    0.9816843611112658D+00, &
    0.9975212478233336D+00, &
    0.9996645373720975D+00, &
    0.6321205588285577D+00, &
    0.4865828809674080D+00, &
    0.3934693402873666D+00, &
    0.3296799539643607D+00, &
    0.4421745996289254D+00, &
    0.1911531694619419D+00, &
    0.6564245437845009D-01, &
    0.1857593622214067D-01 /)
  real ( kind = 8 ) mu
  real ( kind = 8 ), save, dimension ( n_max ) :: mu_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ), save, dimension ( n_max ) :: sigma_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    mu = 0.0D+00
    sigma = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    mu = mu_vec(n_data)
    sigma = sigma_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function gamma_check ( a, b, c )

!*****************************************************************************80
!
!! GAMMA_CHECK checks the parameters of the Gamma PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, logical GAMMA_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical gamma_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMA_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    write ( *, '(a,g14.6)' ) '  B = ', b
    gamma_check = .false.
    return
  end if

  if ( c <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMA_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C <= 0.'
    write ( *, '(a,g14.6)' ) '  C = ', c
    gamma_check = .false.
    return
  end if

  gamma_check = .true.

  return
end
function gamma_inc ( p, x )

!*****************************************************************************80
!
!! GAMMA_INC computes the incomplete Gamma function.
!
!  Discussion:
!
!    GAMMA_INC(P,       0) = 0,
!    GAMMA_INC(P,Infinity) = 1.
!
!    GAMMA_INC(P,X) = Integral ( 0 <= T <= X ) T**(P-1) EXP(-T) DT / GAMMA(P).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    Original FORTRAN77 version by B L Shea.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    BL Shea,
!    Chi-squared and Incomplete Gamma Integral,
!    Algorithm AS239,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 466-473.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the exponent parameter.
!    0.0D+00 < P.
!
!    Input, real ( kind = 8 ) X, the integral limit parameter.
!    If X is less than or equal to 0, GAMMA_INC is returned as 0.
!
!    Output, real ( kind = 8 ) GAMMA_INC, the value of the function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: exp_arg_min = -88.0D+00
  real ( kind = 8 ) gamma_inc
  real ( kind = 8 ) gamma_log
  real ( kind = 8 ), parameter :: overflow = 1.0D+37
  real ( kind = 8 ) p
  real ( kind = 8 ), parameter :: plimit = 1000.0D+00
  real ( kind = 8 ) pn1
  real ( kind = 8 ) pn2
  real ( kind = 8 ) pn3
  real ( kind = 8 ) pn4
  real ( kind = 8 ) pn5
  real ( kind = 8 ) pn6
  real ( kind = 8 ) rn
  real ( kind = 8 ), parameter :: tol = 1.0D-07
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 1.0D+08

  gamma_inc = 0.0D+00

  if ( p <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMA_INC - Fatal error!'
    write ( *, '(a)' ) '  Parameter P <= 0.'
    stop
  end if

  if ( x <= 0.0D+00 ) then
    gamma_inc = 0.0D+00
    return
  end if
!
!  Use a normal approximation if PLIMIT < P.
!
  if ( plimit < p ) then
    pn1 = 3.0D+00 * sqrt ( p ) * ( ( x / p ) ** ( 1.0D+00 / 3.0D+00 ) &
      + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )
    call normal_01_cdf ( pn1, cdf )
    gamma_inc = cdf
    return
  end if
!
!  Is X extremely large compared to P?
!
  if ( xbig < x ) then
    gamma_inc = 1.0D+00
    return
  end if
!
!  Use Pearson's series expansion.
!  (P is not large enough to force overflow in the log of Gamma.
!
  if ( x <= 1.0D+00 .or. x < p ) then

    arg = p * log ( x ) - x - gamma_log ( p + 1.0D+00 )
    c = 1.0D+00
    gamma_inc = 1.0D+00
    a = p

    do

      a = a + 1.0D+00
      c = c * x / a
      gamma_inc = gamma_inc + c

      if ( c <= tol ) then
        exit
      end if

    end do

    arg = arg + log ( gamma_inc )

    if ( exp_arg_min <= arg ) then
      gamma_inc = exp ( arg )
    else
      gamma_inc = 0.0D+00
    end if

  else
!
!  Use a continued fraction expansion.
!
    arg = p * log ( x ) - x - gamma_log ( p )
    a = 1.0D+00 - p
    b = a + x + 1.0D+00
    c = 0.0D+00
    pn1 = 1.0D+00
    pn2 = x
    pn3 = x + 1.0D+00
    pn4 = x * b
    gamma_inc = pn3 / pn4

    do

      a = a + 1.0D+00
      b = b + 2.0D+00
      c = c + 1.0D+00
      pn5 = b * pn3 - a * c * pn1
      pn6 = b * pn4 - a * c * pn2

      if ( 0.0D+00 < abs ( pn6 ) ) then

        rn = pn5 / pn6

        if ( abs ( gamma_inc - rn ) <= min ( tol, tol * rn ) ) then

          arg = arg + log ( gamma_inc )

          if ( exp_arg_min <= arg ) then
            gamma_inc = 1.0D+00 - exp ( arg )
          else
            gamma_inc = 1.0D+00
          end if

          return

        end if

        gamma_inc = rn

      end if

      pn1 = pn3
      pn2 = pn4
      pn3 = pn5
      pn4 = pn6
!
!  Rescale terms in continued fraction if terms are large.
!
      if ( overflow <= abs ( pn5 ) ) then
        pn1 = pn1 / overflow
        pn2 = pn2 / overflow
        pn3 = pn3 / overflow
        pn4 = pn4 / overflow
      end if

    end do

  end if

  return
end
subroutine gamma_inc_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
!
!  Discussion:
!
!    The (normalized) incomplete Gamma function P(A,X) is defined as:
!
!      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
!
!    With this definition, for all A and X,
!
!      0 <= PN(A,X) <= 1
!
!    and
!
!      PN(A,INFINITY) = 1.0
!
!    In Mathematica, the function can be evaluated by:
!
!      1 - GammaRegularized[A,X]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+01, &
    0.10D+01, &
    0.10D+01, &
    0.11D+01, &
    0.11D+01, &
    0.11D+01, &
    0.20D+01, &
    0.20D+01, &
    0.20D+01, &
    0.60D+01, &
    0.60D+01, &
    0.11D+02, &
    0.26D+02, &
    0.41D+02 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7382350532339351D+00, &
    0.9083579897300343D+00, &
    0.9886559833621947D+00, &
    0.3014646416966613D+00, &
    0.7793286380801532D+00, &
    0.9918490284064973D+00, &
    0.9516258196404043D-01, &
    0.6321205588285577D+00, &
    0.9932620530009145D+00, &
    0.7205974576054322D-01, &
    0.5891809618706485D+00, &
    0.9915368159845525D+00, &
    0.1018582711118352D-01, &
    0.4421745996289254D+00, &
    0.9927049442755639D+00, &
    0.4202103819530612D-01, &
    0.9796589705830716D+00, &
    0.9226039842296429D+00, &
    0.4470785799755852D+00, &
    0.7444549220718699D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.30D-01, &
    0.30D+00, &
    0.15D+01, &
    0.75D-01, &
    0.75D+00, &
    0.35D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.15D+00, &
    0.15D+01, &
    0.70D+01, &
    0.25D+01, &
    0.12D+02, &
    0.16D+02, &
    0.25D+02, &
    0.45D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function gamma_log ( x )

!*****************************************************************************80
!
!! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ).
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in references 1 and 2.
!    The program uses rational functions that theoretically approximate
!    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
!    approximation for 12 < X is from Hart et al, while approximations
!    for X < 12.0D+00 are similar to those in Cody and Hillstrom,
!    but are unpublished.
!
!    The accuracy achieved depends on the arithmetic system, the compiler,
!    intrinsic functions, and proper selection of the machine dependent
!    constants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Mathematics of Computation,
!    Volume 21, 1967, pages 198-203.
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thacher, Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function.
!    X must be positive.
!
!    Output, real ( kind = 8 ) GAMMA_LOG, the logarithm of the Gamma
!    function of X.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) BETA, the radix for the floating-point
!    representation.
!
!    Local, integer MAXEXP, the smallest positive power of BETA that overflows.
!
!    Local, real ( kind = 8 ) XBIG, the largest argument for which
!    LN(GAMMA(X)) is representable in the machine, the solution to the equation
!      LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!    Local, real ( kind = 8 ) FRTBIG, a rough estimate of the fourth root
!    of XBIG.
!
!  Approximate values for some important machines are:
!
!                            BETA      MAXEXP         XBIG     FRTBIG
!
!  CRAY-1        (S.P.)        2        8191       9.62D+2461  3.13D+615
!  Cyber 180/855 (S.P.)        2        1070       1.72D+319   6.44D+79
!  IEEE (IBM/XT) (S.P.)        2         128       4.08D+36    1.42D+9
!  IEEE (IBM/XT) (D.P.)        2        1024       2.55D+305   2.25D+76
!  IBM 3033      (D.P.)       16          63       4.29D+73    2.56D+18
!  VAX D-Format  (D.P.)        2         127       2.05D+36    1.20D+9
!  VAX G-Format  (D.P.)        2        1023       1.28D+305   1.89D+76
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = 8 ) corr
  real ( kind = 8 ), parameter :: d1 = -5.772156649015328605195174D-01
  real ( kind = 8 ), parameter :: d2 =  4.227843350984671393993777D-01
  real ( kind = 8 ), parameter :: d4 =  1.791759469228055000094023D+00
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: frtbig = 1.42D+09
  real ( kind = 8 ) gamma_log
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real ( kind = 8 ), parameter :: pnt68 = 0.6796875D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 4.08D+36
  real ( kind = 8 ) xden
  real ( kind = 8 ) xm1
  real ( kind = 8 ) xm2
  real ( kind = 8 ) xm4
  real ( kind = 8 ) xnum
  real ( kind = 8 ) xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0D+00 .or. xbig < x ) then
    gamma_log = huge ( gamma_log )
    return
  end if

  if ( x <= epsilon ( x ) ) then

    res = -log ( x )

  else if ( x <= 1.5D+00 ) then

    if ( x < pnt68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0D+00
      xm1 = ( x - 0.5D+00 ) - 0.5D+00
    end if

    if ( x <= 0.5D+00 .or. pnt68 <= x ) then

      xden = 1.0D+00
      xnum = 0.0D+00

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5D+00 ) - 0.5D+00
      xden = 1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0D+00 ) then

    xm2 = x - 2.0D+00
    xden = 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0D+00 ) then

    xm4 = x - 4.0D+00
    xden = - 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0D+00

    if ( x <= frtbig ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + sqrtpi - 0.5D+00 * corr
    res = res + x * ( corr - 1.0D+00 )

  end if

  gamma_log = res

  return
end
function gamma_log_int ( n )

!*****************************************************************************80
!
!! GAMMA_LOG_INT computes the logarithm of Gamma of an integer N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the logarithm of the
!    Gamma function.  0 < N.
!
!    Output, real ( kind = 8 ) GAMMA_LOG_INT, the logarithm of
!    the Gamma function of N.
!
  implicit none

  real ( kind = 8 ) gamma_log
  real ( kind = 8 ) gamma_log_int
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMA_LOG_INT - Fatal error!'
    write ( *, '(a,i12)' ) '  Illegal input value of N = ', n
    write ( *, '(a)' ) '  But N must be strictly positive.'
    stop
  end if

  gamma_log_int = gamma_log ( real ( n, kind = 8 ) )

  return
end
subroutine gamma_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! GAMMA_MEAN returns the mean of the Gamma PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) mean

  mean = a + b * c

  return
end
subroutine gamma_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! GAMMA_PDF evaluates the Gamma PDF.
!
!  Discussion:
!
!    PDF(A,B,C;X) = exp ( - ( X - A ) / B ) * ( ( X - A ) / B )^(C-1)
!      / ( B * GAMMA ( C ) )
!
!    GAMMA_PDF(A,B,C;X), where C is an integer, is the Erlang PDF.
!    GAMMA_PDF(A,B,1;X) is the Exponential PDF.
!    GAMMA_PDF(0,2,C/2;X) is the Chi Squared PDF with C degrees of freedom.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A <= X
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A controls the location of the peak;  A is often chosen to be 0.0.
!    B is the "scale" parameter; 0.0 < B, and is often 1.0.
!    C is the "shape" parameter; 0.0 < C, and is often 1.0.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) pdf
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x <= a ) then

    pdf = 0.0D+00

  else

    y = ( x - a ) / b

    pdf = y**( c - 1.0D+00 ) / ( b * r8_gamma ( c ) * exp ( y ) )

  end if

  return
end
subroutine gamma_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! GAMMA_SAMPLE samples the Gamma PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 2000
!
!  Author:
!
!    Original FORTRAN77 version by Joachim Ahrens, Ulrich Dieter.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Joachim Ahrens, Ulrich Dieter,
!    Generating Gamma Variates by a Modified Rejection Technique,
!    Communications of the ACM,
!    Volume 25, Number 1, January 1982, pages 47 - 54.
!
!    Joachim Ahrens, Ulrich Dieter,
!    Computer Methods for Sampling from Gamma, Beta, Poisson and
!    Binomial Distributions.
!    Computing,
!    Volume 12, 1974, pages 223 - 246.
!
!    Joachim Ahrens, KD Kohrt, Ulrich Dieter,
!    Algorithm 599,
!    ACM Transactions on Mathematical Software,
!    Volume 9, Number 2, June 1983, pages 255-257.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), parameter :: a1 =   0.3333333D+00
  real ( kind = 8 ), parameter :: a2 = - 0.2500030D+00
  real ( kind = 8 ), parameter :: a3 =   0.2000062D+00
  real ( kind = 8 ), parameter :: a4 = - 0.1662921D+00
  real ( kind = 8 ), parameter :: a5 =   0.1423657D+00
  real ( kind = 8 ), parameter :: a6 = - 0.1367177D+00
  real ( kind = 8 ), parameter :: a7 =   0.1233795D+00
  real ( kind = 8 ) b
  real ( kind = 8 ) bcoef
  real ( kind = 8 ) c
  real ( kind = 8 ) co
  real ( kind = 8 ) d
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) e
  real ( kind = 8 ), parameter :: e1 = 1.0D+00
  real ( kind = 8 ), parameter :: e2 = 0.4999897D+00
  real ( kind = 8 ), parameter :: e3 = 0.1668290D+00
  real ( kind = 8 ), parameter :: e4 = 0.0407753D+00
  real ( kind = 8 ), parameter :: e5 = 0.0102930D+00
  real ( kind = 8 ), parameter :: euler = 2.71828182845904D+00
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) q0
  real ( kind = 8 ), parameter :: q1 =   0.04166669D+00
  real ( kind = 8 ), parameter :: q2 =   0.02083148D+00
  real ( kind = 8 ), parameter :: q3 =   0.00801191D+00
  real ( kind = 8 ), parameter :: q4 =   0.00144121D+00
  real ( kind = 8 ), parameter :: q5 = - 0.00007388D+00
  real ( kind = 8 ), parameter :: q6 =   0.00024511D+00
  real ( kind = 8 ), parameter :: q7 =   0.00024240D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ) si
  real ( kind = 8 ) s2
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) x
!
!  Allow C = 0.
!
  if ( c == 0.0D+00 ) then
    x = a
    return
  end if
!
!  C < 1.
!
  if ( c < 1.0D+00 ) then

    do

      u = r8_uniform_01 ( seed )
      t = 1.0D+00 + c / euler
      p = u * t

      call exponential_01_sample ( seed, s )

      if ( p < 1.0D+00 ) then
        x = exp ( log ( p ) / c )
        if ( x <= s ) then
          exit
        end if
      else
        x = - log ( ( t - p ) / c )
        if ( ( 1.0D+00 - c ) * log ( x ) <= s ) then
          exit
        end if
      end if

    end do

    x = a + b * x
    return
!
!  1 <= C.
!
  else

    s2 = c - 0.5D+00
    s = sqrt ( c - 0.5D+00 )
    d = sqrt ( 32.0D+00 ) - 12.0D+00 * sqrt ( c - 0.5D+00 )

    call normal_01_sample ( seed, t )
    x = ( sqrt ( c - 0.5D+00 ) + 0.5D+00 * t )**2

    if ( 0.0D+00 <= t ) then
      x = a + b * x
      return
    end if

    u = r8_uniform_01 ( seed )

    if ( d * u <= t**3 ) then
      x = a + b * x
      return
    end if

    r = 1.0D+00 / c

    q0 = ( ( ( ( ( ( &
           q7   * r &
         + q6 ) * r &
         + q5 ) * r &
         + q4 ) * r &
         + q3 ) * r &
         + q2 ) * r &
         + q1 ) * r

    if ( c <= 3.686D+00 ) then
      bcoef = 0.463D+00 + s - 0.178D+00 * s2
      si = 1.235D+00
      co = 0.195D+00 / s - 0.079D+00 + 0.016D+00 * s
    else if ( c <= 13.022D+00 ) then
      bcoef = 1.654D+00 + 0.0076D+00 * s2
      si = 1.68D+00 / s + 0.275D+00
      co = 0.062D+00 / s + 0.024D+00
    else
      bcoef = 1.77D+00
      si = 0.75D+00
      co = 0.1515D+00 / s
    end if

    if ( 0.0D+00 < sqrt ( c - 0.5D+00 ) + 0.5D+00 * t ) then

      v = 0.5D+00 * t / s

      if ( 0.25D+00 < abs ( v ) ) then
        q = q0 - s * t + 0.25D+00 * t * t + 2.0D+00 * s2 * log ( 1.0D+00 + v )
      else
        q = q0 + 0.5D+00 * t**2 * ( ( ( ( ( ( &
               a7   * v &
             + a6 ) * v &
             + a5 ) * v &
             + a4 ) * v &
             + a3 ) * v &
             + a2 ) * v &
             + a1 ) * v
      end if

      if ( log ( 1.0D+00 - u ) <= q ) then
        x = a + b * x
        return
      end if

    end if

    do

      call exponential_01_sample ( seed, e )

      u = r8_uniform_01 ( seed )

      u = 2.0D+00 * u - 1.0D+00
      t = bcoef + sign ( si * e, u )

      if ( -0.7187449D+00 <= t ) then

        v = 0.5D+00 * t / s

        if ( 0.25D+00 < abs ( v ) ) then
          q = q0 - s * t + 0.25D+00 * t**2 + 2.0D+00 * s2 * log ( 1.0D+00 + v )
        else
          q = q0 + 0.5D+00 * t**2 * ( ( ( ( ( ( &
               a7   * v &
             + a6 ) * v &
             + a5 ) * v &
             + a4 ) * v &
             + a3 ) * v &
             + a2 ) * v &
             + a1 ) * v
        end if

        if ( 0.0D+00 < q ) then

          if ( 0.5D+00 < q ) then
            w = exp ( q ) - 1.0D+00
          else
            w = ( ( ( ( &
                    e5   * q &
                  + e4 ) * q &
                  + e3 ) * q &
                  + e2 ) * q &
                  + e1 ) * q
          end if

          if ( co * abs ( u ) <= w * exp ( e - 0.5D+00 * t**2 ) ) then
            x = a + b * ( s + 0.5D+00 * t )**2
            return
          end if

        end if

      end if

    end do

  end if

  return
end
subroutine gamma_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! GAMMA_VARIANCE returns the variance of the Gamma PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) variance

  variance = b * b * c

  return
end
subroutine genlogistic_cdf ( x, a, b, c, cdf )

!*****************************************************************************80
!
!! GENLOGISTIC_CDF evaluates the Generalized Logistic CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  y = ( x - a ) / b

  cdf = 1.0D+00 / ( 1.0D+00 + exp ( - y ) )**c

  return
end
subroutine genlogistic_cdf_inv ( cdf, a, b, c, x )

!*****************************************************************************80
!
!! GENLOGISTIC_CDF_INV inverts the Generalized Logistic CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENLOGISTIC_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = - huge ( x )
  else if ( cdf < 1.0D+00 ) then
    x = a - b * log ( cdf**( - 1.0D+00 / c ) - 1.0D+00 )
  else if ( 1.0D+00 == cdf ) then
    x = huge ( x )
  end if

  return
end
function genlogistic_check ( a, b, c )

!*****************************************************************************80
!
!! GENLOGISTIC_CHECK checks the parameters of the Generalized Logistic CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, logical GENLOGISTIC_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical genlogistic_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENLOGISTIC_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    genlogistic_check = .false.
    return
  end if

  if ( c <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENLOGISTIC_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C <= 0.'
    genlogistic_check = .false.
    return
  end if

  genlogistic_check = .true.

  return
end
subroutine genlogistic_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! GENLOGISTIC_MEAN returns the mean of the Generalized Logistic PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) digamma
  real ( kind = 8 ) euler_constant
  real ( kind = 8 ) mean

  mean = a + b * ( euler_constant ( ) + digamma ( c ) )

  return
end
subroutine genlogistic_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! GENLOGISTIC_PDF evaluates the Generalized Logistic PDF.
!
!  Discussion:
!
!    PDF(A,B,C;X) = ( C / B ) * exp ( ( A - X ) / B ) /
!      ( ( 1 + exp ( ( A - X ) / B ) )**(C+1) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  y = ( x - a ) / b

  pdf = ( c / b ) * exp ( - y ) / ( 1.0D+00 + exp ( - y ) )**( c + 1.0D+00 )

  return
end
subroutine genlogistic_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! GENLOGISTIC_SAMPLE samples the Generalized Logistic PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call genlogistic_cdf_inv ( cdf, a, b, c, x )

  return
end
subroutine genlogistic_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! GENLOGISTIC_VARIANCE returns the variance of the Generalized Logistic PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) trigamma
  real ( kind = 8 ) variance

  variance = b * b * ( pi * pi / 6.0D+00 + trigamma ( c ) )

  return
end
subroutine geometric_cdf ( x, a, cdf )

!*****************************************************************************80
!
!! GEOMETRIC_CDF evaluates the Geometric CDF.
!
!  Discussion:
!
!    CDF(X,P) is the probability that there will be at least one
!    successful trial in the first X Bernoulli trials, given that
!    the probability of success in a single trial is P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the maximum number of trials.
!
!    Input, real ( kind = 8 ) A, the probability of success on one trial.
!    0.0D+00 <= A <= 1.0.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) x

  if ( x <= 0 ) then
    cdf = 0.0D+00
  else if ( a == 0.0D+00 ) then
    cdf = 0.0D+00
  else if ( a == 1.0D+00 ) then
    cdf = 1.0D+00
  else
    cdf = 1.0D+00 - ( 1.0D+00 - a )**x
  end if

  return
end
subroutine geometric_cdf_inv ( cdf, a, x )

!*****************************************************************************80
!
!! GEOMETRIC_CDF_INV inverts the Geometric CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0D+00
!
!    Input, real ( kind = 8 ) A, the probability of success on one trial.
!    0.0D+00 <= A <= 1.0.
!
!    Output, integer ( kind = 4 ) X, the corresponding value of X.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GEOMETRIC_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( a == 1.0D+00 ) then
    x = 1
  else if ( a == 0.0D+00 ) then
    x = huge ( x )
  else
    x = 1 + int ( log ( 1.0D+00 - cdf ) / log ( 1.0D+00 - a ) )
  end if

  return
end
subroutine geometric_cdf_values ( n_data, x, p, cdf )

!*****************************************************************************80
!
!! GEOMETRIC_CDF_VALUES returns values of the geometric CDF.
!
!  Discussion:
!
!    The geometric or Pascal probability density function gives the
!    probability that the first success will happen on the X-th Bernoulli
!    trial, given that the probability of a success on a single trial is P.
!
!    The value of CDF ( X, P ) is the probability that the first success
!    will happen on or before the X-th trial.
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`DiscreteDistributions`]
!      dist = GeometricDistribution [ p ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!    Daniel Zwillinger, Stephen Kokoska,
!    CRC Standard Probability and Statistics Tables and Formulae,
!    Chapman and Hall / CRC Press, 2000.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) X, the number of trials.
!
!    Output, real ( kind = 8 ) P, the probability of success
!    on one trial.
!
!    Output, real ( kind = 8 ) CDF, the cumulative density function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 14

  real ( kind = 8 ) cdf
  real ( kind = 8 ), save, dimension ( n_max ) :: cdf_vec = (/ &
    0.1900000000000000D+00, &
    0.2710000000000000D+00, &
    0.3439000000000000D+00, &
    0.6861894039100000D+00, &
    0.3600000000000000D+00, &
    0.4880000000000000D+00, &
    0.5904000000000000D+00, &
    0.9141006540800000D+00, &
    0.7599000000000000D+00, &
    0.8704000000000000D+00, &
    0.9375000000000000D+00, &
    0.9843750000000000D+00, &
    0.9995117187500000D+00, &
    0.9999000000000000D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_max ) :: p_vec = (/ &
    0.1D+00, &
    0.1D+00, &
    0.1D+00, &
    0.1D+00, &
    0.2D+00, &
    0.2D+00, &
    0.2D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    0.9D+00 /)
  integer ( kind = 4 ) x
  integer ( kind = 4 ), save, dimension ( n_max ) :: x_vec = (/ &
    1,  2,  3, 10, 1, &
    2,  3, 10,  3, 3, &
    3,  5, 10,  3 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0
    p = 0.0D+00
    cdf = 0.0D+00
  else
    x = x_vec(n_data)
    p = p_vec(n_data)
    cdf = cdf_vec(n_data)
  end if

  return
end
function geometric_check ( a )

!*****************************************************************************80
!
!! GEOMETRIC_CHECK checks the parameter of the Geometric CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the probability of success on one trial.
!    0.0D+00 <= A <= 1.0.
!
!    Output, logical GEOMETRIC_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  logical geometric_check

  if ( a < 0.0D+00 .or. 1.0D+00 < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GEOMETRIC_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A < 0 or 1 < A.'
    geometric_check = .false.
    return
  end if

  geometric_check = .true.

  return
end
subroutine geometric_mean ( a, mean )

!*****************************************************************************80
!
!! GEOMETRIC_MEAN returns the mean of the Geometric PDF.
!
!  Discussion:
!
!    MEAN is the expected value of the number of trials required
!    to obtain a single success.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the probability of success on one trial.
!    0.0D+00 <= A <= 1.0.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) mean

  mean = 1.0D+00 / a

  return
end
subroutine geometric_pdf ( x, a, pdf )

!*****************************************************************************80
!
!! GEOMETRIC_PDF evaluates the Geometric PDF.
!
!  Discussion:
!
!    PDF(A;X) = A * ( 1 - A )**(X-1)
!
!    PDF(A;X) is the probability that exactly X Bernoulli trials, each
!    with probability of success A, will be required to achieve
!    a single success.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the number of trials.
!    0 < X
!
!    Input, real ( kind = 8 ) A, the probability of success on one trial.
!    0.0D+00 <= A <= 1.0.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
!
!  Special cases.
!
  if ( x < 1 ) then

    pdf = 0.0D+00

  else if ( a == 0.0D+00 ) then

    pdf = 0.0D+00

  else if ( a == 1.0D+00 ) then

    if ( x == 1 ) then
      pdf = 1.0D+00
    else
      pdf = 0.0D+00
    end if

  else

    pdf = a * ( 1.0D+00 - a )**( x - 1 )

  end if

  return
end
subroutine geometric_sample ( a, seed, x )

!*****************************************************************************80
!
!! GEOMETRIC_SAMPLE samples the Geometric PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the probability of success on one trial.
!    0.0D+00 <= A <= 1.0.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x

  cdf = r8_uniform_01 ( seed )

  call geometric_cdf_inv ( cdf, a, x )

  return
end
subroutine geometric_variance ( a, variance )

!*****************************************************************************80
!
!! GEOMETRIC_VARIANCE returns the variance of the Geometric PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the probability of success on one trial.
!    0.0D+00 <= A <= 1.0.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) variance

  variance = ( 1.0D+00 - a ) / ( a * a )

  return
end
subroutine gompertz_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! GOMPERTZ_CDF evaluates the Gompertz CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Norman Johnson, Samuel Kotz, Balakrishnan,
!    Continuous Univariate Distributions, second edition,
!    Wiley, 1994,
!    QA273.6.J6
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    1 < A, 0 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then
    cdf = 0.0D+00
  else
    cdf = 1.0D+00 - exp ( - b * ( a**x - 1.0D+00 ) / log ( a ) )
  end if

  return
end
subroutine gompertz_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! GOMPERTZ_CDF_INV inverts the Gompertz CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Norman Johnson, Samuel Kotz, Balakrishnan,
!    Continuous Univariate Distributions, second edition,
!    Wiley, 1994,
!    QA273.6.J6
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    1 < A, 0 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GOMPERTZ_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf < 1.0D+00 ) then
    x = log ( 1.0D+00 - log ( 1.0D+00 - cdf ) * log ( a ) / b  ) / log ( a )
  else
    x = huge ( x )
  end if

  return
end
function gompertz_check ( a, b )

!*****************************************************************************80
!
!! GOMPERTZ_CHECK checks the parameters of the Gompertz PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Norman Johnson, Samuel Kotz, Balakrishnan,
!    Continuous Univariate Distributions, second edition,
!    Wiley, 1994,
!    QA273.6.J6
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    1 < A, 0 < B.
!
!    Output, logical GOMPERTZ_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical gompertz_check

  if ( a <= 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GOMPERTZ_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 1.0!'
    gompertz_check = .false.
    return
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GOMPERTZ_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.0!'
    gompertz_check = .false.
    return
  end if

  gompertz_check = .true.

  return
end
subroutine gompertz_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! GOMPERTZ_PDF evaluates the Gompertz PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = B * A**X / exp ( B * ( A**X - 1 ) / log ( A ) )
!
!    for
!
!      0.0 <= X
!      1.0 <  A
!      0.0 <  B
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Norman Johnson, Samuel Kotz, Balakrishnan,
!    Continuous Univariate Distributions, second edition,
!    Wiley, 1994,
!    QA273.6.J6
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    1 < A, 0 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then

    pdf = 0.0D+00

  else if ( 1.0D+00 < a ) then

    pdf = exp ( log ( b ) + x * log ( a ) &
      - ( b / log ( a ) ) * ( a**x - 1.0D+00 ) )

  end if

  return
end
subroutine gompertz_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! GOMPERTZ_SAMPLE samples the Gompertz PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    1 < A, 0 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call gompertz_cdf_inv ( cdf, a, b, x )

  return
end
subroutine gumbel_cdf ( x, cdf )

!*****************************************************************************80
!
!! GUMBEL_CDF evaluates the Gumbel CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  cdf = exp ( - exp ( - x ) )

  return
end
subroutine gumbel_cdf_inv ( cdf, x )

!*****************************************************************************80
!
!! GUMBEL_CDF_INV inverts the Gumbel CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GUMBEL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x =  - log ( - log ( cdf ) )

  return
end
subroutine gumbel_mean ( mean )

!*****************************************************************************80
!
!! GUMBEL_MEAN returns the mean of the Gumbel PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) euler_constant
  real ( kind = 8 ) mean

  mean = euler_constant ( )

  return
end
subroutine gumbel_pdf ( x, pdf )

!*****************************************************************************80
!
!! GUMBEL_PDF evaluates the Gumbel PDF.
!
!  Discussion:
!
!    PDF(X) = exp ( -X ) * exp ( - exp ( -X  ) ).
!
!    GUMBEL_PDF(X) = EXTREME_PDF(0,1;X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein, editor,
!    CRC Concise Encylopedia of Mathematics,
!    CRC Press, 1998.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  pdf = exp ( - x - exp ( - x ) )

  return
end
subroutine gumbel_sample ( seed, x )

!*****************************************************************************80
!
!! GUMBEL_SAMPLE samples the Gumbel PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
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
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call gumbel_cdf_inv ( cdf, x )

  return
end
subroutine gumbel_variance ( variance )

!*****************************************************************************80
!
!! GUMBEL_VARIANCE returns the variance of the Gumbel PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) variance

  variance = pi * pi / 6.0D+00

  return
end
subroutine half_normal_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! HALF_NORMAL_CDF evaluates the Half Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) x

  if ( x <= a ) then
    cdf = 0.0D+00
  else
    call normal_cdf ( x, a, b, cdf2 )
    cdf = 2.0D+00 * cdf2 - 1.0D+00
  end if

  return
end
subroutine half_normal_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! HALF_NORMAL_CDF_INV inverts the Half Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALF_NORMAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  cdf2 = 0.5D+00 * ( cdf + 1.0D+00 )

  call normal_cdf_inv ( cdf2, a, b, x )

  return
end
function half_normal_check ( a, b )

!*****************************************************************************80
!
!! HALF_NORMAL_CHECK checks the parameters of the Half Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, logical HALF_NORMAL_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical half_normal_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALF_NORMAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    half_normal_check = .false.
    return
  end if

  half_normal_check = .true.

  return
end
subroutine half_normal_mean ( a, b, mean )

!*****************************************************************************80
!
!! HALF_NORMAL_MEAN returns the mean of the Half Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  mean = a + b * sqrt ( 2.0D+00 / pi )

  return
end
subroutine half_normal_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! HALF_NORMAL_PDF evaluates the Half Normal PDF.
!
!  Discussion:
!
!    PDF(A,B;X) =
!      sqrt ( 2 / PI ) * ( 1 / B ) * exp ( - 0.5D+00 * ( ( X - A ) / B )^2 )
!
!    for A <= X
!
!    The Half Normal PDF is a special case of both the Chi PDF and the
!    Folded Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A <= X
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x <= a ) then

    pdf = 0.0D+00

  else

    y = ( x - a ) / b

    pdf = sqrt ( 2.0D+00 / pi ) * ( 1.0D+00 / b ) * exp ( - 0.5D+00 * y * y )

  end if

  return
end
subroutine half_normal_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! HALF_NORMAL_SAMPLE samples the Half Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call half_normal_cdf_inv ( cdf, a, b, x )

  return
end
subroutine half_normal_variance ( a, b, variance )

!*****************************************************************************80
!
!! HALF_NORMAL_VARIANCE returns the variance of the Half Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) variance

  variance = b * b * ( 1.0D+00 - 2.0D+00 / pi )

  return
end
subroutine hypergeometric_cdf ( x, n, m, l, cdf )

!*****************************************************************************80
!
!! HYPERGEOMETRIC_CDF evaluates the Hypergeometric CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the CDF.
!
!    Input, integer ( kind = 4 ) N, the number of balls selected.
!    0 <= N <= L.
!
!    Input, integer ( kind = 4 ) M, the number of white balls.
!    0 <= M <= L.
!
!    Input, integer ( kind = 4 ) L, the number of balls to select from.
!    0 <= L.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) c1_log
  real ( kind = 8 ) c2_log
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  call binomial_coef_log ( l - m, n, c1_log )
  call binomial_coef_log ( l, n, c2_log )

  pdf = exp ( c1_log - c2_log )
  cdf = pdf

  do x2 = 0, x - 1

    pdf = pdf * real ( ( m - x2 ) * ( n - x2 ), kind = 8 ) &
      / real ( ( x2 + 1 ) * ( l - m - n + x2 + 1 ), kind = 8 )

    cdf = cdf + pdf

  end do

  return
end
subroutine hypergeometric_cdf_values ( n_data, sam, suc, pop, n, fx )

!*****************************************************************************80
!
!! HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric CDF.
!
!  Discussion:
!
!    CDF(X)(A,B) is the probability of at most X successes in A trials,
!    given that the probability of success on a single trial is B.
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`DiscreteDistributions`]
!      dist = HypergeometricDistribution [ sam, suc, pop ]
!      CDF [ dist, n ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition, CRC Press, 1996, pages 651-652.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) SAM, integer SUC, integer POP, the sample
!    size, success size, and population parameters of the function.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6001858177500578D-01, &
    0.2615284665839845D+00, &
    0.6695237889132748D+00, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.5332595856827856D+00, &
    0.1819495964117640D+00, &
    0.4448047017527730D-01, &
    0.9999991751316731D+00, &
    0.9926860896560750D+00, &
    0.8410799901444538D+00, &
    0.3459800113391901D+00, &
    0.0000000000000000D+00, &
    0.2088888139634505D-02, &
    0.3876752992448843D+00, &
    0.9135215248834896D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     7,  8,  9, 10, &
     6,  6,  6,  6, &
     6,  6,  6,  6, &
     0,  0,  0,  0 /)
  integer ( kind = 4 ) pop
  integer ( kind = 4 ), save, dimension ( n_max ) :: pop_vec = (/ &
    100, 100, 100, 100, &
    100, 100, 100, 100, &
    100, 100, 100, 100, &
    90,  200, 1000, 10000 /)
  integer ( kind = 4 ) sam
  integer ( kind = 4 ), save, dimension ( n_max ) :: sam_vec = (/ &
    10, 10, 10, 10, &
     6,  7,  8,  9, &
    10, 10, 10, 10, &
    10, 10, 10, 10 /)
  integer ( kind = 4 ) suc
  integer ( kind = 4 ), save, dimension ( n_max ) :: suc_vec = (/ &
    90, 90, 90, 90, &
    90, 90, 90, 90, &
    10, 30, 50, 70, &
    90, 90, 90, 90 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    sam = 0
    suc = 0
    pop = 0
    n = 0
    fx = 0.0D+00
  else
    sam = sam_vec(n_data)
    suc = suc_vec(n_data)
    pop = pop_vec(n_data)
    n = n_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function hypergeometric_check ( n, m, l )

!*****************************************************************************80
!
!! HYPERGEOMETRIC_CHECK checks the parameters of the Hypergeometric CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of balls selected.
!    0 <= N <= L.
!
!    Input, integer ( kind = 4 ) M, the number of white balls in the population.
!    0 <= M <= L.
!
!    Input, integer ( kind = 4 ) L, the number of balls to select from.
!    0 <= L.
!
!    Output, logical HYPERGEOMETRIC_CHECK, is true if the parameters are legal.
!
  implicit none

  logical hypergeometric_check
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( n < 0 .or. l < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HYPERGEOMETRIC_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Input N is out of range.'
    hypergeometric_check = .false.
    return
  end if

  if ( m < 0 .or. l < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HYPERGEOMETRIC_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Input M is out of range.'
    hypergeometric_check = .false.
    return
  end if

  if ( l < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HYPERGEOMETRIC_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Input L is out of range.'
    hypergeometric_check = .false.
    return
  end if

  hypergeometric_check = .true.

  return
end
subroutine hypergeometric_mean ( n, m, l, mean )

!*****************************************************************************80
!
!! HYPERGEOMETRIC_MEAN returns the mean of the Hypergeometric PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of balls selected.
!    0 <= N <= L.
!
!    Input, integer ( kind = 4 ) M, the number of white balls in the population.
!    0 <= M <= L.
!
!    Input, integer ( kind = 4 ) L, the number of balls to select from.
!    0 <= L.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real ( kind = 8 ) mean
  integer ( kind = 4 ) n

  mean = real ( n * m, kind = 8 ) / real ( l, kind = 8 )

  return
end
subroutine hypergeometric_pdf ( x, n, m, l, pdf )

!*****************************************************************************80
!
!! HYPERGEOMETRIC_PDF evaluates the Hypergeometric PDF.
!
!  Discussion:
!
!    PDF(N,M,L;X) = C(M,X) * C(L-M,N-X) / C(L,N).
!
!    PDF(N,M,L;X) is the probability of drawing X white balls in a
!    single random sample of size N from a population containing
!    M white balls and a total of L balls.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the desired number of white balls.
!    0 <= X <= N, usually, although any value of X can be given.
!
!    Input, integer ( kind = 4 ) N, the number of balls selected.
!    0 <= N <= L.
!
!    Input, integer ( kind = 4 ) M, the number of white balls.
!    0 <= M <= L.
!
!    Input, integer ( kind = 4 ) L, the number of balls to select from.
!    0 <= L.
!
!    Output, real ( kind = 8 ) PDF, the probability of exactly K white balls.
!
  implicit none

  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) pdf
  real ( kind = 8 ) pdf_log
  integer ( kind = 4 ) x
!
!  Special cases.
!
  if ( x < 0 ) then

    pdf = 1.0D+00

  else if ( n < x ) then

    pdf = 0.0D+00

  else if ( m < x ) then

    pdf = 0.0D+00

  else if ( l < x ) then

    pdf = 0.0D+00

  else if ( n == 0 ) then

    if ( x == 0 ) then
      pdf = 1.0D+00
    else
      pdf = 0.0D+00
    end if

  else

    call binomial_coef_log ( m, x, c1 )
    call binomial_coef_log ( l-m, n-x, c2 )
    call binomial_coef_log ( l, n, c3 )

    pdf_log = c1 + c2 - c3

    pdf = exp ( pdf_log )

  end if

  return
end
subroutine hypergeometric_sample ( n, m, l, seed, x )

!*****************************************************************************80
!
!! HYPERGEOMETRIC_SAMPLE samples the Hypergeometric PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 165.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of balls selected.
!    0 <= N <= L.
!
!    Input, integer ( kind = 4 ) M, the number of white balls.
!    0 <= M <= L.
!
!    Input, integer ( kind = 4 ) L, the number of balls to select from.
!    0 <= L.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c1_log
  real ( kind = 8 ) c2_log
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u
  integer ( kind = 4 ) x

  call binomial_coef_log ( l - m, n, c1_log )
  call binomial_coef_log ( l, n, c2_log )

  a = exp ( c1_log - c2_log )
  b = a

  u = r8_uniform_01 ( seed )

  x = 0

  do while ( a < u )

    b = b * real ( ( m - x ) * ( n - x ), kind = 8 ) &
      / real ( ( x + 1 ) * ( l - m - n + x + 1 ), kind = 8 )

    a = a + b

    x = x + 1

  end do

  return
end
subroutine hypergeometric_variance ( n, m, l, variance )

!*****************************************************************************80
!
!! HYPERGEOMETRIC_VARIANCE returns the variance of the Hypergeometric PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of balls selected.
!    0 <= N <= L.
!
!    Input, integer ( kind = 4 ) M, the number of white balls.
!    0 <= M <= L.
!
!    Input, integer ( kind = 4 ) L, the number of balls to select from.
!    0 <= L.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) variance

  variance = real ( n * m * ( l - m ) * ( l - n ), kind = 8 ) &
    / real ( l * l * ( l - 1 ), kind = 8 )

  return
end
function i4_factorial ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL returns N!.
!
!  Discussion:
!
!    N! = Product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the function.
!    0 <= N.
!
!    Output, real ( kind = 8 ) I4_FACTORIAL, the factorial of N.
!
  implicit none

  real ( kind = 8 ) i4_factorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_FACTORIAL - Fatal error!'
    write ( *, '(a)' ) '  N < 0.'
    stop
  end if

  i4_factorial = 1.0D+00

  do i = 2, n
    i4_factorial = i4_factorial * real ( i, kind = 8 )
  end do

  return
end
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
!    use I4_HUGE() and HUGE interchangeably.
!
!    However, when using the G95, the values returned by HUGE were
!    not equal to 2147483647, apparently, and were causing severe
!    and obscure errors in my random number generator, which needs to
!    add I4_HUGE to the seed whenever the seed is negative.  So I
!    am backing away from invoking HUGE, whereas I4_HUGE is under
!    my control.
!
!    Explanation: because under G95 the default integer type is 64 bits!
!    So HUGE ( 1 ) = a very very huge integer indeed, whereas
!    I4_HUGE ( ) = the same old 32 bit big value.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
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
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine i4row_max ( m, n, x, ixmax, xmax )

!*****************************************************************************80
!
!! I4ROW_MAX returns the maximums of rows of an I4ROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) X(M,N), the array to be examined.
!
!    Output, integer ( kind = 4 ) IXMAX(M); IXMAX(I) is the column of X in
!    which the maximum for row I occurs.
!
!    Output, integer ( kind = 4 ) XMAX(M), the maximums of the rows of X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ixmax(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x(m,n)
  integer ( kind = 4 ) xmax(m)

  do i = 1, m

    ixmax(i) = 1
    xmax(i) = x(i,1)
    do j = 2, n
      if ( xmax(i) < x(i,j) ) then
        ixmax(i) = j
        xmax(i) = x(i,j)
      end if
    end do

  end do

  return
end
subroutine i4row_mean ( m, n, a, mean )

!*****************************************************************************80
!
!! I4ROW_MEAN returns the means of the rows of an I4ROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) A(M,N), the array.
!
!    Output, real ( kind = 8 ) MEAN(M), the mean of each row.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mean(m)

  do i = 1, m
    mean(i) = sum ( a(i,1:n) ) / real ( n, kind = 8 )
  end do

  return
end
subroutine i4row_min ( m, n, x, ixmin, xmin )

!*****************************************************************************80
!
!! I4ROW_MIN returns the minimums of rows of an I4ROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the array.
!
!    Input, integer ( kind = 4 ) X(M,N), the array to be examined.
!
!    Output, integer ( kind = 4 ) IXMIN(M); IXMIN(I) is the column of X in which
!    the minimum for row I occurs.
!
!    Output, integer ( kind = 4 ) XMIN(M), the minimums of the rows of X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ixmin(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x(m,n)
  integer ( kind = 4 ) xmin(m)

  do i = 1, m

    ixmin(i) = 1
    xmin(i) = x(i,1)
    do j = 2, n
      if ( x(i,j) < xmin(i) ) then
        ixmin(i) = j
        xmin(i) = x(i,j)
      end if
    end do

  end do

  return
end
subroutine i4row_variance ( m, n, a, variance )

!*****************************************************************************80
!
!! I4ROW_VARIANCE returns the variances of the rows of an I4ROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) A(M,N), the array.
!
!    Output, real ( kind = 8 ) VARIANCE(M), the variance of each row.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance(m)

  do i = 1, m

    mean = real ( sum ( a(i,1:n) ), kind = 8 ) / real ( n, kind = 8 )

    variance(i) = 0.0D+00
    do j = 1, n
      variance(i) = variance(i) + ( real ( a(i,j), kind = 8 ) - mean )**2
    end do

    if ( 1 < n ) then
      variance(i) = variance(i) / real ( n - 1, kind = 8 )
    else
      variance(i) = 0.0D+00
    end if

  end do

  return
end
subroutine i4vec_max ( n, iarray, index, imax )

!*****************************************************************************80
!
!! I4VEC_MAX computes the maximum element of an I4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) IARRAY(N), the array.
!
!    Output, integer ( kind = 4 ) INDEX, the index of the largest entry.
!
!    Output, integer ( kind = 4 ) IMAX, the value of the largest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) index
  integer ( kind = 4 ) imax

  if ( n <= 0 ) then

    imax = 0
    index = 0

  else

    imax = iarray(1)
    index = 1
    do i = 2, n

      if ( imax < iarray(i) ) then
        imax = iarray(i)
        index = i
      end if

    end do

  end if

  return
end
subroutine i4vec_mean ( n, x, mean )

!*****************************************************************************80
!
!! I4VEC_MEAN returns the mean of an I4VEC.
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
!    Input, integer ( kind = 4 ) X(N), the vector whose mean is desired.
!
!    Output, real ( kind = 8 ) MEAN, the mean, or average, of
!    the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) mean
  integer ( kind = 4 ) x(n)

  mean = real ( sum ( x(1:n) ), kind = 8 ) / real ( n, kind = 8 )

  return
end
subroutine i4vec_min ( n, iarray, index, imin )

!*****************************************************************************80
!
!! I4VEC_MIN computes the minimum element of an I4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) IARRAY(N), the array.
!
!    Output, integer ( kind = 4 ) INDEX, the index of the smallest entry.
!
!    Output, integer ( kind = 4 ) IMIN, the value of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) index

  if ( n <= 0 ) then

    imin = 0
    index = 0

  else

    imin = iarray(1)
    index = 1
    do i = 2, n
      if ( iarray(i) < imin ) then
        imin = iarray(i)
        index = i
      end if
    end do

  end if

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) big
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(i8,1x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(i8,1x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i8,i11)' ) i, a(i)
    end do
  end if

  return
end
subroutine i4vec_run_count ( n, a, run_count )

!*****************************************************************************80
!
!! I4VEC_RUN_COUNT counts runs of equal values in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    A run is a sequence of equal values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be examined.
!
!    Output, integer ( kind = 4 ) RUN_COUNT, the number of runs.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) run_count
  integer ( kind = 4 ) test

  run_count = 0

  if ( n < 1 ) then
    return
  end if

  test = 0

  do i = 1, n

    if ( i == 1 .or. a(i) /= test ) then
      run_count = run_count + 1
      test = a(i)
    end if

  end do

  return
end
subroutine i4vec_variance ( n, x, variance )

!*****************************************************************************80
!
!! I4VEC_VARIANCE returns the variance of an I4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) X(N), the vector whose variance is desired.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) mean
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(n)

  call i4vec_mean ( n, x, mean )

  variance = sum ( ( real ( x(1:n), kind = 8 ) - mean )**2 )

  if ( 1 < n ) then
    variance = variance / real ( n - 1, kind = 8 )
  else
    variance = 0.0D+00
  end if

  return
end
subroutine inverse_gaussian_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! INVERSE_GAUSSIAN_CDF evaluates the Inverse Gaussian CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!    0.0D+00 < X.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf1
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  if ( x <= 0.0D+00 ) then

    cdf = 0.0D+00

  else

    x1 = sqrt ( b / x ) * ( x - a ) / a
    call normal_01_cdf ( x1, cdf1 )

    x2 = - sqrt ( b / x ) * ( x + a ) / a
    call normal_01_cdf ( x2, cdf2 )

    cdf = cdf1 + exp ( 2.0D+00 * b / a ) * cdf2

  end if

  return
end
function inverse_gaussian_check ( a, b )

!*****************************************************************************80
!
!! INVERSE_GAUSSIAN_CHECK checks the parameters of the Inverse Gaussian CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, logical INVERSE_GAUSSIAN_CHECK, is true if the parameters
!    are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical inverse_gaussian_check

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INVERSE_GAUSSIAN_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    inverse_gaussian_check = .false.
    return
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INVERSE_GAUSSIAN_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    inverse_gaussian_check = .false.
    return
  end if

  inverse_gaussian_check = .true.

  return
end
subroutine inverse_gaussian_mean ( a, b, mean )

!*****************************************************************************80
!
!! INVERSE_GAUSSIAN_MEAN returns the mean of the Inverse Gaussian PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine inverse_gaussian_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! INVERSE_GAUSSIAN_PDF evaluates the Inverse Gaussian PDF.
!
!  Discussion:
!
!    The Inverse Gaussian PDF is also known as the Wald PDF
!    and the Inverse Normal PDF.
!
!    PDF(A,B;X)
!      = sqrt ( B / ( 2 * PI * X^3 ) )
!        * exp ( - B * ( X - A )^2 / ( 2.0D+00 * A^2 * X ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 < X
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then
    pdf = 0.0D+00
  else
    pdf = sqrt ( b / ( 2.0D+00 * pi * x**3 ) ) * &
      exp ( - b * ( x - a )**2 / ( 2.0D+00 * a * a * x ) )
  end if

  return
end
subroutine inverse_gaussian_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! INVERSE_GAUSSIAN_SAMPLE samples the Inverse Gaussian PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) phi
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  phi = b / a
  call normal_01_sample ( seed, z )
  y = z * z

  t = 1.0D+00 + 0.5D+00 * ( y - sqrt ( 4.0D+00 * phi * y + y * y ) ) / phi
  u = r8_uniform_01 ( seed )

  if ( u * ( 1.0D+00 + t ) <= 1.0D+00 ) then
    x = a * t
  else
    x = a / t
  end if

  return
end
subroutine inverse_gaussian_variance ( a, b, variance )

!*****************************************************************************80
!
!! INVERSE_GAUSSIAN_VARIANCE returns the variance of the Inverse Gaussian PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = a**3 / b

  return
end
subroutine laplace_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! LAPLACE_CDF evaluates the Laplace CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  y = ( x - a ) / b

  if ( x <= a ) then
    cdf = 0.5D+00 * exp ( y )
  else
    cdf = 1.0D+00 - 0.5D+00 * exp ( - y )
  end if

  return
end
subroutine laplace_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! LAPLACE_CDF_INV inverts the Laplace CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAPLACE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf <= 0.5D+00 ) then
    x = a + b * log ( 2.0D+00 * cdf )
  else
    x = a - b * log ( 2.0D+00 * ( 1.0D+00 - cdf ) )
  end if

  return
end
subroutine laplace_cdf_values ( n_data, mu, beta, x, fx )

!*****************************************************************************80
!
!! LAPLACE_CDF_VALUES returns some values of the Laplace CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = LaplaceDistribution [ mu, beta ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) MU, the mean of the distribution.
!
!    Output, real ( kind = 8 ) BETA, the shape parameter.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) beta
  real ( kind = 8 ), save, dimension ( n_max ) :: beta_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.8160602794142788D+00, &
    0.9323323583816937D+00, &
    0.9751064658160680D+00, &
    0.6967346701436833D+00, &
    0.6417343447131054D+00, &
    0.6105996084642976D+00, &
    0.5906346234610091D+00, &
    0.5000000000000000D+00, &
    0.3032653298563167D+00, &
    0.1839397205857212D+00, &
    0.1115650800742149D+00 /)
  real ( kind = 8 ) mu
  real ( kind = 8 ), save, dimension ( n_max ) :: mu_vec = (/ &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.0000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    mu = 0.0D+00
    beta = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    mu = mu_vec(n_data)
    beta = beta_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function laplace_check ( a, b )

!*****************************************************************************80
!
!! LAPLACE_CHECK checks the parameters of the Laplace PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, logical LAPLACE_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical laplace_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAPLACE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    laplace_check = .false.
    return
  end if

  laplace_check = .true.

  return
end
subroutine laplace_mean ( a, b, mean )

!*****************************************************************************80
!
!! LAPLACE_MEAN returns the mean of the Laplace PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine laplace_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! LAPLACE_PDF evaluates the Laplace PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = exp ( - abs ( X - A ) / B ) / ( 2 * B )
!
!    The Laplace PDF is also known as the Double Exponential PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  pdf = exp ( - abs ( x - a ) / b ) / ( 2.0D+00 * b )

  return
end
subroutine laplace_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! LAPLACE_SAMPLE samples the Laplace PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call laplace_cdf_inv ( cdf, a, b, x )

  return
end
subroutine laplace_variance ( a, b, variance )

!*****************************************************************************80
!
!! LAPLACE_VARIANCE returns the variance of the Laplace PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = 2.0D+00 * b * b

  return
end
function lerch ( a, b, c )

!*****************************************************************************80
!
!! LERCH estimates the Lerch transcendent function.
!
!  Discussion:
!
!    The Lerch transcendent function is defined as:
!
!      LERCH ( A, B, C ) = Sum ( 0 <= K < Infinity ) A**K / ( C + K )**B
!
!    excluding any term with ( C + K ) = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein, editor,
!    CRC Concise Encylopedia of Mathematics,
!    CRC Press, 1998.
!
!  Thanks:
!
!    Oscar van Vlijmen
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the function.
!
!    Output, real ( kind = 8 ) LERCH, an approximation to the Lerch
!    transcendent function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a_k
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) k
  real ( kind = 8 ) lerch
  real ( kind = 8 ) sum2
  real ( kind = 8 ) sum2_old

  sum2 = 0.0D+00
  k = 0
  a_k = 1.0D+00

  do

    sum2_old = sum2

    if ( c + real ( k, kind = 8 ) == 0.0D+00 ) then
      k = k + 1
      a_k = a_k * a
      cycle
    end if

    sum2 = sum2 + a_k / ( c + real ( k, kind = 8 ) )**b

    if ( sum2 <= sum2_old ) then
      exit
    end if

    k = k + 1
    a_k = a_k * a

  end do

  lerch = sum2

  return
end
subroutine levy_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! LEVY_CDF evaluates the Levy CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!    Normally, A <= X.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) error_f
  real ( kind = 8 ) x

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEVY_CDF - Fatal error!'
    write ( *, '(a)' ) '  Input parameter B <= 0.0'
    stop
  end if

  if ( x <= a ) then
    cdf = 0.0D+00
  else
    cdf = 1.0D+00 - error_f ( sqrt ( b / ( 2.0D+00 * ( x - a ) ) ) )
  end if

  return
end
subroutine levy_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! LEVY_CDF_INV inverts the Levy CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf1
  real ( kind = 8 ) x
  real ( kind = 8 ) x1

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEVY_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEVY_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  Input parameter B <= 0.0'
    stop
  end if

  cdf1 = 1.0D+00 - 0.5D+00 * cdf
  call normal_01_cdf_inv ( cdf1, x1 )
  x = a + b / ( x1 * x1 )

  return
end
subroutine levy_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! LEVY_PDF evaluates the Levy PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = sqrt ( B / ( 2 * PI ) )
!               * exp ( - B / ( 2 * ( X - A ) )
!               / ( X - A )^(3/2)
!
!    for A <= X.
!
!    Note that the Levy PDF does not have a finite mean or variance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    Normally, A <= X.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEVY_PDF - Fatal error!'
    write ( *, '(a)' ) '  Input parameter B <= 0.0'
    stop
  end if

  if ( x <= a ) then
    pdf = 0.0D+00
  else
    pdf = sqrt ( b / ( 2.0D+00 * pi ) ) &
        * exp ( - b / ( 2.0D+00 * ( x - a ) ) ) &
        / sqrt ( ( x - a )**3 )
  end if

  return
end
subroutine levy_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! LEVY_SAMPLE samples the Levy PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call levy_cdf_inv ( cdf, a, b, x )

  return
end
subroutine logistic_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! LOGISTIC_CDF evaluates the Logistic CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  cdf = 1.0D+00 / ( 1.0D+00 + exp ( ( a - x ) / b ) )

  return
end
subroutine logistic_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! LOGISTIC_CDF_INV inverts the Logistic CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOGISTIC_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a - b * log ( ( 1.0D+00 - cdf ) / cdf )

  return
end
subroutine logistic_cdf_values ( n_data, mu, beta, x, fx )

!*****************************************************************************80
!
!! LOGISTIC_CDF_VALUES returns some values of the Logistic CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = LogisticDistribution [ mu, beta ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) MU, the mean of the distribution.
!
!    Output, real ( kind = 8 ) BETA, the shape parameter of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) beta
  real ( kind = 8 ), save, dimension ( n_max ) :: beta_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.8807970779778824D+00, &
    0.9820137900379084D+00, &
    0.9975273768433652D+00, &
    0.6224593312018546D+00, &
    0.5825702064623147D+00, &
    0.5621765008857981D+00, &
    0.5498339973124779D+00, &
    0.6224593312018546D+00, &
    0.5000000000000000D+00, &
    0.3775406687981454D+00, &
    0.2689414213699951D+00 /)
  real ( kind = 8 ) mu
  real ( kind = 8 ), save, dimension ( n_max ) :: mu_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    mu = 0.0D+00
    beta = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    mu = mu_vec(n_data)
    beta = beta_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function logistic_check ( a, b )

!*****************************************************************************80
!
!! LOGISTIC_CHECK checks the parameters of the Logistic CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, logical LOGISTIC_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical logistic_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOGISTIC_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    logistic_check = .false.
    return
  end if

  logistic_check = .true.

  return
end
subroutine logistic_mean ( a, b, mean )

!*****************************************************************************80
!
!! LOGISTIC_MEAN returns the mean of the Logistic PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine logistic_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! LOGISTIC_PDF evaluates the Logistic PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = exp ( ( A - X ) / B ) /
!      ( B * ( 1 + exp ( ( A - X ) / B ) )^2 )
!
!    The Logistic PDF is also known as the Sech-Squared PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) temp
  real ( kind = 8 ) x

  temp = exp ( ( a - x ) / b )

  pdf = temp / ( b * ( 1.0D+00 + temp )**2 )

  return
end
subroutine logistic_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! LOGISTIC_SAMPLE samples the Logistic PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call logistic_cdf_inv ( cdf, a, b, x )

  return
end
subroutine logistic_variance ( a, b, variance )

!*****************************************************************************80
!
!! LOGISTIC_VARIANCE returns the variance of the Logistic PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) variance

  variance = pi * pi * b * b / 3.0D+00

  return
end
subroutine log_normal_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! LOG_NORMAL_CDF evaluates the Lognormal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 < X.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) logx
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then

    cdf = 0.0D+00

  else

    logx = log ( x )

    call normal_cdf ( logx, a, b, cdf )

  end if

  return
end
subroutine log_normal_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! LOG_NORMAL_CDF_INV inverts the Lognormal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) logx
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOG_NORMAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  call normal_cdf_inv ( cdf, a, b, logx )

  x = exp ( logx )

  return
end
subroutine log_normal_cdf_values ( n_data, mu, sigma, x, fx )

!*****************************************************************************80
!
!! LOG_NORMAL_CDF_VALUES returns some values of the Log Normal CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = LogNormalDistribution [ mu, sigma ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) MU, the mean of the distribution.
!
!    Output, real ( kind = 8 ) SIGMA, the shape parameter of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.2275013194817921D-01, &
    0.2697049307349095D+00, &
    0.5781741008028732D+00, &
    0.7801170895122241D+00, &
    0.4390310097476894D+00, &
    0.4592655190218048D+00, &
    0.4694258497695908D+00, &
    0.4755320473858733D+00, &
    0.3261051056816658D+00, &
    0.1708799040927608D+00, &
    0.7343256357952060D-01, &
    0.2554673736161761D-01 /)
  real ( kind = 8 ) mu
  real ( kind = 8 ), save, dimension ( n_max ) :: mu_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ), save, dimension ( n_max ) :: sigma_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    mu = 0.0D+00
    sigma = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    mu = mu_vec(n_data)
    sigma = sigma_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function log_normal_check ( a, b )

!*****************************************************************************80
!
!! LOG_NORMAL_CHECK checks the parameters of the Lognormal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, logical LOG_NORMAL_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical log_normal_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOG_NORMAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    log_normal_check = .false.
    return
  end if

  log_normal_check = .true.

  return
end
subroutine log_normal_mean ( a, b, mean )

!*****************************************************************************80
!
!! LOG_NORMAL_MEAN returns the mean of the Lognormal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = exp ( a + 0.5D+00 * b * b )

  return
end
subroutine log_normal_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! LOG_NORMAL_PDF evaluates the Lognormal PDF.
!
!  Discussion:
!
!    PDF(A,B;X)
!      = exp ( - 0.5 * ( ( log ( X ) - A ) / B )^2 )
!        / ( B * X * sqrt ( 2 * PI ) )
!
!    The Lognormal PDF is also known as the Cobb-Douglas PDF,
!    and as the Antilog_normal PDF.
!
!    The Lognormal PDF describes a variable X whose logarithm
!    is normally distributed.
!
!    The special case A = 0, B = 1 is known as Gilbrat's PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 < X
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then
    pdf = 0.0D+00
  else
    pdf = exp ( - 0.5D+00 * ( ( log ( x ) - a ) / b )**2 ) &
      / ( b * x * sqrt ( 2.0D+00 * pi ) )
  end if

  return
end
subroutine log_normal_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! LOG_NORMAL_SAMPLE samples the Lognormal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call log_normal_cdf_inv ( cdf, a, b, x )

  return
end
subroutine log_normal_variance ( a, b, variance )

!*****************************************************************************80
!
!! LOG_NORMAL_VARIANCE returns the variance of the Lognormal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = exp ( 2.0D+00 * a + b * b ) * ( exp ( b * b ) - 1.0D+00 )

  return
end
subroutine log_series_cdf ( x, a, cdf )

!*****************************************************************************80
!
!! LOG_SERIES_CDF evaluates the Logarithmic Series CDF.
!
!  Discussion:
!
!    Simple summation is used, with a recursion to generate successive
!    values of the PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Thanks:
!
!    Oscar van Vlijmen
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the PDF.
!    0 < X
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A < 1.0.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  cdf = 0.0D+00

  do x2 = 1, x

    if ( x2 == 1 ) then
      pdf = - a / log ( 1.0D+00 - a )
    else
      pdf = real ( x2 - 1, kind = 8 ) * a * pdf / real ( x2, kind = 8 )
    end if

    cdf = cdf + pdf

  end do

  return
end
subroutine log_series_cdf_inv ( cdf, a, x )

!*****************************************************************************80
!
!! LOG_SERIES_CDF_INV inverts the Logarithmic Series CDF.
!
!  Discussion:
!
!    Simple summation is used.  The only protection against an
!    infinite loop caused by roundoff is that X cannot be larger
!    than 1000.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A < 1.0.
!
!    Output, real ( kind = 8 ) X, the argument of the CDF for which
!    CDF(X-1) <= CDF <= CDF(X).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  integer ( kind = 4 ), parameter :: xmax = 1000

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOG_SERIES_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  cdf2 = 0.0D+00
  x = 1

  do while ( cdf2 < cdf .and. x < xmax )

    if ( x == 1 ) then
      pdf = - a / log ( 1.0D+00 - a )
    else
      pdf = real ( x - 1, kind = 8 ) * a * pdf / real ( x, kind = 8 )
    end if

    cdf2 = cdf2 + pdf

    x = x + 1

  end do

  return
end
subroutine log_series_cdf_values ( n_data, t, n, fx )

!*****************************************************************************80
!
!! LOG_SERIES_CDF_VALUES returns some values of the log series CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`DiscreteDistributions`]
!      dist = LogSeriesDistribution [ t ]
!      CDF [ dist, n ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) T, the parameter of the function.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 29

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.9491221581029903D+00, &
    0.9433541128559735D+00, &
    0.9361094611773272D+00, &
    0.9267370278044118D+00, &
    0.9141358246245129D+00, &
    0.8962840235449100D+00, &
    0.8690148741955517D+00, &
    0.8221011541254772D+00, &
    0.7213475204444817D+00, &
    0.6068261510845583D+00, &
    0.5410106403333613D+00, &
    0.4970679476476894D+00, &
    0.4650921887927060D+00, &
    0.4404842934597863D+00, &
    0.4207860535926143D+00, &
    0.4045507673897055D+00, &
    0.3908650337129266D+00, &
    0.2149757685421097D+00, &
    0.0000000000000000D+00, &
    0.2149757685421097D+00, &
    0.3213887739704539D+00, &
    0.3916213575531612D+00, &
    0.4437690508633213D+00, &
    0.4850700239649681D+00, &
    0.5191433267738267D+00, &
    0.5480569580144867D+00, &
    0.5731033910767085D+00, &
    0.5951442521714636D+00, &
    0.6147826594068904D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     1, 1, 1, 1, 1, &
     1, 1, 1, 1, 1, &
     1, 1, 1, 1, 1, &
     1, 1, 1, 0, 1, &
     2, 3, 4, 5, 6, &
     7, 8, 9, 10 /)
  real ( kind = 8 ) t
  real ( kind = 8 ), save, dimension ( n_max ) :: t_vec = (/ &
    0.1000000000000000D+00, &
    0.1111111111111111D+00, &
    0.1250000000000000D+00, &
    0.1428571428571429D+00, &
    0.1666666666666667D+00, &
    0.2000000000000000D+00, &
    0.2500000000000000D+00, &
    0.3333333333333333D+00, &
    0.5000000000000000D+00, &
    0.6666666666666667D+00, &
    0.7500000000000000D+00, &
    0.8000000000000000D+00, &
    0.8333333333333333D+00, &
    0.8571485714857149D+00, &
    0.8750000000000000D+00, &
    0.8888888888888889D+00, &
    0.9000000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00, &
    0.9900000000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    t = 0.0D+00
    n = 0
    fx = 0.0D+00
  else
    t = t_vec(n_data)
    n = n_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function log_series_check ( a )

!*****************************************************************************80
!
!! LOG_SERIES_CHECK checks the parameter of the Logarithmic Series PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A < 1.0.
!
!    Output, logical LOG_SERIES_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  logical log_series_check

  if ( a <= 0.0D+00 .or. 1.0D+00 <= a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOG_SERIES_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.0D+00 or 1.0D+00 <= A'
    log_series_check = .false.
    return
  end if

  log_series_check = .true.

  return
end
subroutine log_series_mean ( a, mean )

!*****************************************************************************80
!
!! LOG_SERIES_MEAN returns the mean of the Logarithmic Series PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A < 1.0.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) mean

  mean = - a / ( ( 1.0D+00 - a ) * log ( 1.0D+00 - a ) )

  return
end
subroutine log_series_pdf ( x, a, pdf )

!*****************************************************************************80
!
!! LOG_SERIES_PDF evaluates the Logarithmic Series PDF.
!
!  Discussion:
!
!    PDF(A;X) = - A**X / ( X * log ( 1 - A ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the PDF.
!    0 < X
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A < 1.0.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  if ( x <= 0 ) then
    pdf = 0.0D+00
  else
    pdf = - a**x / ( real ( x, kind = 8 ) * log ( 1.0D+00 - a ) )
  end if

  return
end
subroutine log_series_sample ( a, seed, x )

!*****************************************************************************80
!
!! LOG_SERIES_SAMPLE samples the Logarithmic Series PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Luc Devroye,
!    Non-Uniform Random Variate Generation,
!    Springer-Verlag, 1986, page 547.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A < 1.0.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  integer ( kind = 4 ) x

  u = r8_uniform_01 ( seed )
  v = r8_uniform_01 ( seed )

  x = int ( 1.0D+00 + log ( v ) / ( log ( 1.0D+00 - ( 1.0D+00 - a )**u ) ) )

  return
end
subroutine log_series_variance ( a, variance )

!*****************************************************************************80
!
!! LOG_SERIES_VARIANCE returns the variance of the Logarithmic Series PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A < 1.0.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) variance

  alpha = - 1.0D+00 / log ( 1.0D+00 - a )

  variance = a * alpha * ( 1.0D+00 - alpha * a ) / ( 1.0D+00 - a )**2

  return
end
subroutine log_uniform_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! LOG_UNIFORM_CDF evaluates the Log Uniform CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x <= a ) then
    cdf = 0.0D+00
  else if ( x < b ) then
    cdf = ( log ( x ) - log ( a ) ) / ( log ( b ) - log ( a ) )
  else
    cdf = 1.0D+00
  end if

  return
end
subroutine log_uniform_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! LOG_UNIFORM_CDF_INV inverts the Log Uniform CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOG_UNIFORM_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a * exp ( ( log ( b ) - log ( a ) ) * cdf )

  return
end
function log_uniform_check ( a, b )

!*****************************************************************************80
!
!! LOG_UNIFORM_CHECK checks the parameters of the Log Uniform CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    1.0 < A < B.
!
!    Output, logical LOG_UNIFORM_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical log_uniform_check

  if ( a <= 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOG_UNIFORM_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 1.'
    log_uniform_check = .false.
    return
  end if

  if ( b <= a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOG_UNIFORM_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= A.'
    log_uniform_check = .false.
    return
  end if

  log_uniform_check = .true.

  return
end
subroutine log_uniform_mean ( a, b, mean )

!*****************************************************************************80
!
!! LOG_UNIFORM_MEAN returns the mean of the Log Uniform PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    1.0 < A < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = ( b - a ) / ( log ( b ) - log ( a ) )

  return
end
subroutine log_uniform_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! LOG_UNIFORM_PDF evaluates the Log Uniform PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = 1 / ( X * ( log ( B ) - log ( A ) ) ) for A <= X <= B
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    1.0 < A < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x < a ) then
    pdf = 0.0D+00
  else if ( x <= b ) then
    pdf = 1.0D+00 / ( x * ( log ( b ) - log ( a ) ) )
  else
    pdf = 0.0D+00
  end if

  return
end
subroutine log_uniform_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! LOG_UNIFORM_SAMPLE samples the Log Uniform PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    1.0 < A < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call log_uniform_cdf_inv ( cdf, a, b, x )

  return
end
subroutine lorentz_cdf ( x, cdf )

!*****************************************************************************80
!
!! LORENTZ_CDF evaluates the Lorentz CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  cdf = 0.5D+00 + atan ( x ) / pi

  return
end
subroutine lorentz_cdf_inv ( cdf, x )

!*****************************************************************************80
!
!! LORENTZ_CDF_INV inverts the Lorentz CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LORENTZ_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = tan ( pi * ( cdf - 0.5D+00 ) )

  return
end
subroutine lorentz_mean ( mean )

!*****************************************************************************80
!
!! LORENTZ_MEAN returns the mean of the Lorentz PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) mean

  mean = 0.0D+00

  return
end
subroutine lorentz_pdf ( x, pdf )

!*****************************************************************************80
!
!! LORENTZ_PDF evaluates the Lorentz PDF.
!
!  Discussion:
!
!    PDF(X) = 1 / ( PI * ( 1 + X^2 ) )
!
!    The chief interest of the Lorentz PDF is that it is easily
!    inverted, and can be used to dominate other PDF's in an
!    acceptance/rejection method.
!
!    LORENTZ_PDF(X) = CAUCHY_PDF(0,1;X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  pdf = 1.0D+00 / ( pi * ( 1.0D+00 + x * x ) )

  return
end
subroutine lorentz_sample ( seed, x )

!*****************************************************************************80
!
!! LORENTZ_SAMPLE samples the Lorentz PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
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
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call lorentz_cdf_inv ( cdf, x )

  return
end
subroutine lorentz_variance ( variance )

!*****************************************************************************80
!
!! LORENTZ_VARIANCE returns the variance of the Lorentz PDF.
!
!  Discussion:
!
!    The variance of the Lorentz PDF is not well defined.  This routine
!    is made available for completeness only, and simply returns
!    a "very large" number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VARIANCE, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) variance

  variance = huge ( variance )

  return
end
subroutine maxwell_cdf ( x, a, cdf )

!*****************************************************************************80
!
!! MAXWELL_CDF evaluates the Maxwell CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0 < A.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) gamma_inc
  real ( kind = 8 ) p2
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  if ( x <= 0.0D+00 ) then

    cdf = 0.0D+00

  else

    x2 = x / a
    p2 = 1.5D+00

    cdf = gamma_inc ( p2, x2 )

  end if

  return
end
subroutine maxwell_cdf_inv ( cdf, a, x )

!*****************************************************************************80
!
!! MAXWELL_CDF_INV inverts the Maxwell CDF.
!
!  Discussion:
!
!    A simple bisection method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0 < A.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf1
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) cdf3
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  real ( kind = 8 ), parameter :: tol = 0.0001D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAXWELL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = 0.0D+00
    return
  else if ( 1.0D+00 == cdf ) then
    x = huge ( x )
    return
  end if

  x1 = 0.0D+00
  cdf1 = 0.0D+00

  x2 = 1.0D+00

  do

    call maxwell_cdf ( x2, a, cdf2 )

    if ( cdf < cdf2 ) then
      exit
    end if

    x2 = 2.0D+00 * x2

    if ( 1000000.0D+00 < x2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MAXWELL_CDF_INV - Fatal error!'
      write ( *, '(a)' ) '  Initial bracketing effort fails.'
      stop
    end if

  end do
!
!  Now use bisection.
!
  it = 0

  do

    it = it + 1

    x3 = 0.5D+00 * ( x1 + x2 )
    call maxwell_cdf ( x3, a, cdf3 )

    if ( abs ( cdf3 - cdf ) < tol ) then
      x = x3
      exit
    end if

    if ( it_max < it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MAXWELL_CDF_INV - Fatal error!'
      write ( *, '(a)' ) '  Iteration limit exceeded.'
      stop
    end if

    if ( sign ( 1.0D+00, cdf3 - cdf ) == sign ( 1.0D+00, cdf1 - cdf ) ) then
      x1 = x3
      cdf1 = cdf3
    else
      x2 = x3
      cdf2 = cdf3
    end if

  end do

  return
end
function maxwell_check ( a )

!*****************************************************************************80
!
!! MAXWELL_CHECK checks the parameters of the Maxwell CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0 < A.
!
!    Output, logical MAXWELL_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  logical maxwell_check

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAXWELL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.0.'
    maxwell_check = .false.
    return
  end if

  maxwell_check = .true.

  return
end
subroutine maxwell_mean ( a, mean )

!*****************************************************************************80
!
!! MAXWELL_MEAN returns the mean of the Maxwell PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0 < A.
!
!    Output, real ( kind = 8 ) MEAN, the mean value.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) mean
  real ( kind = 8 ) r8_gamma

  mean = sqrt ( 2.0D+00 ) * a * r8_gamma ( 2.0D+00 ) / r8_gamma ( 1.5D+00 )

  return
end
subroutine maxwell_pdf ( x, a, pdf )

!*****************************************************************************80
!
!! MAXWELL_PDF evaluates the Maxwell PDF.
!
!  Discussion:
!
!    PDF(A;X) = exp ( - 0.5D+00 * ( X / A )^2 ) * ( X / A )^2 /
!      ( sqrt ( 2 ) * A * GAMMA ( 1.5D+00 ) )
!
!    MAXWELL_PDF(A;X) = CHI_PDF(0,A,3;X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0 < X
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0 < A.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) pdf
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x <= 0.0D+00 ) then

    pdf = 0.0D+00

  else

    y = x / a

    pdf = exp ( - 0.5D+00 * y * y ) * y * y &
      / ( sqrt ( 2.0D+00 ) * a * r8_gamma ( 1.5D+00 ) )

  end if

  return
end
subroutine maxwell_sample ( a, seed, x )

!*****************************************************************************80
!
!! MAXWELL_SAMPLE samples the Maxwell PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0 < A.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  a2 = 3.0D+00
  call chi_square_sample ( a2, seed, x )

  x = a * sqrt ( x )

  return
end
subroutine maxwell_variance ( a, variance )

!*****************************************************************************80
!
!! MAXWELL_VARIANCE returns the variance of the Maxwell PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0 < A.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) variance

  variance = a * a * ( 3.0D+00 - 2.0D+00 &
    * ( r8_gamma ( 2.0D+00 ) / r8_gamma ( 1.5D+00 ) )**2 )

  return
end
function multicoef_check ( nfactor, factor )

!*****************************************************************************80
!
!! MULTICOEF_CHECK checks the parameters of the multinomial coefficient.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NFACTOR, the number of factors.
!    1 <= NFACTOR.
!
!    Input, integer ( kind = 4 ) FACTOR(NFACTOR), contains the factors.
!    0.0D+00 <= FACTOR(I).
!
!    Output, logical MULTICOEF_CHECK, is true if the parameters are legal.
!
  implicit none

  integer ( kind = 4 ) nfactor

  integer ( kind = 4 ) factor(nfactor)
  integer ( kind = 4 ) i
  logical multicoef_check

  if ( nfactor < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MULTICOEF_CHECK - Fatal error!'
    write ( *, '(a)' ) '  NFACTOR < 1.'
    multicoef_check = .false.
    return
  end if

  do i = 1, nfactor

    if ( factor(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MULTICOEF_CHECK - Fatal error'
      write ( *, '(a,i8)' ) '  Factor ', I
      write ( *, '(a,i8)' ) '  = ', factor(i)
      write ( *, '(a)' ) '  But this value must be nonnegative.'
      multicoef_check = .false.
      return
    end if

  end do

  multicoef_check = .true.

  return
end
subroutine multinomial_coef1 ( nfactor, factor, ncomb )

!*****************************************************************************80
!
!! MULTINOMIAL_COEF1 computes a Multinomial coefficient.
!
!  Discussion:
!
!    The multinomial coefficient is a generalization of the binomial
!    coefficient.  It may be interpreted as the number of combinations of
!    N objects, where FACTOR(1) objects are indistinguishable of type 1,
!    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
!    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
!
!    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
!
!    The log of the gamma function is used, to avoid overflow.
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
!    Input, integer ( kind = 4 ) NFACTOR, the number of factors.
!    1 <= NFACTOR.
!
!    Input, integer ( kind = 4 ) FACTOR(NFACTOR), contains the factors.
!    0.0D+00 <= FACTOR(I).
!
!    Output, integer ( kind = 4 ) NCOMB, the value of the multinomial coefficient.
!
  implicit none

  integer ( kind = 4 ) nfactor

  logical              check
  real ( kind = 8 ) facn
  integer ( kind = 4 ) factor(nfactor)
  real ( kind = 8 ) factorial_log
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  logical              multicoef_check
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncomb

  check = multicoef_check ( nfactor, factor )

  if ( .not. check ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MULTINOMIAL_COEF1 - Fatal error!'
    write ( *, '(a)' ) '  MULTICOEF_CHECK failed.'
    ncomb = - i4_huge ( )
    return
  end if
!
!  The factors sum to N.
!
  n = sum ( factor(1:nfactor) )

  facn = factorial_log ( n )

  do i = 1, nfactor

    facn = facn - factorial_log ( factor(i) )

  end do

  ncomb = nint ( exp ( facn ) )

  return
end
subroutine multinomial_coef2 ( nfactor, factor, ncomb )

!*****************************************************************************80
!
!! MULTINOMIAL_COEF2 computes a Multinomial coefficient.
!
!  Discussion:
!
!    The multinomial coefficient is a generalization of the binomial
!    coefficient.  It may be interpreted as the number of combinations of
!    N objects, where FACTOR(1) objects are indistinguishable of type 1,
!    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
!    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
!
!    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
!
!    A direct method is used, which should be exact.  However, there
!    is a possibility of intermediate overflow of the result.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NFACTOR, the number of factors.
!    1 <= NFACTOR.
!
!    Input, integer ( kind = 4 ) FACTOR(NFACTOR), contains the factors.
!    0.0D+00 <= FACTOR(I).
!
!    Output, integer ( kind = 4 ) NCOMB, the value of the multinomial coefficient.
!
  implicit none

  integer ( kind = 4 ) nfactor

  logical              check
  integer ( kind = 4 ) factor(nfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical              multicoef_check
  integer ( kind = 4 ) ncomb

  check = multicoef_check ( nfactor, factor )

  if ( .not. check ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MULTINOMIAL_COEF2 - Fatal error!'
    write ( *, '(a)' ) '  MULTICOEF_CHECK failed.'
    ncomb = - i4_huge ( )
    return
  end if

  ncomb = 1
  k = 0

  do i = 1, nfactor

    do j = 1, factor(i)
      k = k + 1
      ncomb = ( ncomb * k ) / j
    end do

  end do

  return
end
function multinomial_check ( a, b, c )

!*****************************************************************************80
!
!! MULTINOMIAL_CHECK checks the parameters of the Multinomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of trials.
!
!    Input, integer ( kind = 4 ) B, the number of outcomes possible on one trial.
!    1 <= B.
!
!    Input, real ( kind = 8 ) C(B).  C(I) is the probability of outcome I on
!    any trial.
!    0.0D+00 <= C(I) <= 1.0D+00,
!    Sum ( 1 <= I <= B ) C(I) = 1.0.
!
!    Output, logical MULTINOMIAL_CHECK, is true if the parameters are legal.
!
  implicit none

  integer ( kind = 4 ) b

  integer ( kind = 4 ) a
  real ( kind = 8 ) c(b)
  real ( kind = 8 ) c_sum
  integer ( kind = 4 ) i
  logical multinomial_check

  if ( b < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MULTINOMIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B < 1.'
    multinomial_check = .false.
    return
  end if

  do i = 1, b

    if ( c(i) < 0.0D+00 .or. 1.0D+00 < c(i) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MULTINOMIAL_CHECK - Fatal error!'
      write ( *, '(a)' ) '  Input C(I) is out of range.'
      multinomial_check = .false.
      return
    end if

  end do

  c_sum = sum ( c )

  if ( 0.0001D+00 < abs ( 1.0D+00 - c_sum ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MULTINOMIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  The probabilities do not sum to 1.'
    multinomial_check = .false.
    return
  end if

  multinomial_check = .true.

  return
end
subroutine multinomial_covariance ( a, b, c, covariance )

!*****************************************************************************80
!
!! MULTINOMIAL_COVARIANCE returns the covariances of the Multinomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of trials.
!
!    Input, integer ( kind = 4 ) B, the number of outcomes possible on one trial.
!    1 <= B.
!
!    Input, real ( kind = 8 ) C(B).  C(I) is the probability of outcome I on
!    any trial.
!    0.0D+00 <= C(I) <= 1.0D+00,
!    SUM ( 1 <= I <= B) C(I) = 1.0.
!
!    Output, real ( kind = 8 ) COVARIANCE(B,B), the covariance matrix.
!
  implicit none

  integer ( kind = 4 ) b

  integer ( kind = 4 ) a
  real ( kind = 8 ) c(b)
  real ( kind = 8 ) covariance(b,b)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, b
    do j = 1, b

      if ( i == j ) then
        covariance(i,j) = real ( a, kind = 8 ) * c(i) * ( 1.0D+00 - c(i) )
      else
        covariance(i,j) = - real ( a, kind = 8 ) * c(i) * c(j)
      end if

    end do
  end do

  return
end
subroutine multinomial_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! MULTINOMIAL_MEAN returns the means of the Multinomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of trials.
!
!    Input, integer ( kind = 4 ) B, the number of outcomes possible on one trial.
!    1 <= B.
!
!    Input, real ( kind = 8 ) C(B).  C(I) is the probability of outcome I on
!    any trial.
!    0.0D+00 <= C(I) <= 1.0D+00,
!    SUM ( 1 <= I <= B) C(I) = 1.0.
!
!    Output, real ( kind = 8 ) MEAN(B), MEAN(I) is the expected value of the
!    number of outcome I in N trials.
!
  implicit none

  integer ( kind = 4 ) b

  integer ( kind = 4 ) a
  real ( kind = 8 ) c(b)
  real ( kind = 8 ) mean(b)

  mean(1:b) = real ( a, kind = 8 ) * c(1:b)

  return
end
subroutine multinomial_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! MULTINOMIAL_PDF computes a Multinomial PDF.
!
!  Discussion:
!
!    PDF(A,B,C;X) = Comb(A,B,X) * Product ( 1 <= I <= B ) C(I)**X(I)
!
!    where Comb(A,B,X) is the multinomial coefficient
!      C( A; X(1), X(2), ..., X(B) )
!
!    PDF(A,B,C;X) is the probability that in A trials there
!    will be exactly X(I) occurrences of event I, whose probability
!    on one trial is C(I), for I from 1 to B.
!
!    As soon as A or B gets large, the number of possible X's explodes,
!    and the probability of any particular X can become extremely small.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X(B); X(I) counts the number of occurrences of
!    outcome I, out of the total of A trials.
!
!    Input, integer ( kind = 4 ) A, the total number of trials.
!
!    Input, integer ( kind = 4 ) B, the number of different possible outcomes on
!    one trial.
!
!    Input, real ( kind = 8 ) C(B); C(I) is the probability of outcome I on
!    any one trial.
!
!    Output, real ( kind = 8 ) PDF, the value of the multinomial PDF.
!
  implicit none

  integer ( kind = 4 ) b

  integer ( kind = 4 ) a
  real ( kind = 8 ) c(b)
  real ( kind = 8 ) gamma_log
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  real ( kind = 8 ) pdf_log
  integer ( kind = 4 ) x(b)
!
!  To try to avoid overflow, do the calculation in terms of logarithms.
!  Note that Gamma(A+1) = A factorial.
!
  pdf_log = gamma_log ( real ( a + 1, kind = 8 ) )

  do i = 1, b
    pdf_log = pdf_log + x(i) * log ( c(i) ) &
      - gamma_log ( real ( x(i) + 1, kind = 8 ) )
  end do

  pdf = exp ( pdf_log )

  return
end
subroutine multinomial_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! MULTINOMIAL_SAMPLE samples the Multinomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Luc Devroye,
!    Non-Uniform Random Variate Generation,
!    Springer-Verlag, New York, 1986, page 559.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the total number of trials.
!    0 <= A.
!
!    Input, integer ( kind = 4 ) B, the number of outcomes possible on
!    one trial.  1 <= B.
!
!    Input, real ( kind = 8 ) C(B).  C(I) is the probability of outcome I on
!    any trial.
!    0.0D+00 <= C(I) <= 1.0D+00,
!    sum ( 1 <= I <= B) C(I) = 1.0.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) X(B); X(I) is the number of
!    occurrences of event I during the N trials.
!
  implicit none

  integer ( kind = 4 ) b

  integer ( kind = 4 ) a
  real ( kind = 8 ) c(b)
  integer ( kind = 4 ) ifactor
  integer ( kind = 4 ) ntot
  real ( kind = 8 ) prob
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sum2
  integer ( kind = 4 ) x(b)

  ntot = a

  sum2 = 1.0D+00

  x(1:b) = 0

  do ifactor = 1, b - 1

    prob = c(ifactor) / sum2
!
!  Generate a binomial random deviate for NTOT trials with
!  single trial success probability PROB.
!
    call binomial_sample ( ntot, prob, seed, x(ifactor) )

    ntot = ntot - x(ifactor)
    if ( ntot <= 0 ) then
      return
    end if

    sum2 = sum2 - c(ifactor)

  end do
!
!  The last factor gets what's left.
!
  x(b) = ntot

  return
end
subroutine multinomial_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! MULTINOMIAL_VARIANCE returns the variances of the Multinomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number of trials.
!
!    Input, integer ( kind = 4 ) B, the number of outcomes possible on one trial.
!    1 <= B.
!
!    Input, real ( kind = 8 ) C(B).  C(I) is the probability of outcome I on
!    any trial.
!    0.0D+00 <= C(I) <= 1.0D+00,
!    sum ( 1 <= I <= B ) C(I) = 1.0.
!
!    Output, real ( kind = 8 ) VARIANCE(B), VARIANCE(I) is the variance of the
!    total number of events of type I.
!
  implicit none

  integer ( kind = 4 ) b

  integer ( kind = 4 ) a
  real ( kind = 8 ) c(b)
  integer ( kind = 4 ) i
  real ( kind = 8 ) variance(b)

  do i = 1, b
    variance(i) = real ( a, kind = 8 ) * c(i) * ( 1.0D+00 - c(i) )
  end do

  return
end
subroutine multivariate_normal_sample ( n, mean, covar_factor, seed, x )

!*****************************************************************************80
!
!! MULTIVARIATE_NORMAL_SAMPLE samples the Multivariate Normal PDF.
!
!  Discussion:
!
!    PDF ( Mean(1:N), S(1:N,1:N); X(1:N) ) = 1 / ( 2 * pi * det ( S ) )^(N/2)
!      * exp ( - ( X - Mean )' * inverse ( S ) * ( X - Mean ) / 2 )
!
!    Here,
!
!      X is the argument vector of length N,
!      Mean is the mean vector of length N,
!      S is an N by N positive definite symmetric covariance matrix.
!
!    The properties of S guarantee that it has a lower triangular
!    matrix L, the Cholesky factor, such that S = L * L'.  It is the
!    matrix L, rather than S, that is required by this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 167.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) MEAN(N), the mean vector.
!
!    Input, real ( kind = 8 ) COVAR_FACTOR(N,N), the lower triangular Cholesky
!    factor L of the covariance matrix S.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample point of the distribution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) covar_factor(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) z

  do i = 1, n

    call normal_01_sample ( seed, z )

    x(i) = mean(i)

    do j = 1, i
      x(i) = x(i) + covar_factor(i,j) * z
    end do

  end do

  return
end
subroutine nakagami_cdf ( x, a, b, c, cdf )

!*****************************************************************************80
!
!! NAKAGAMI_CDF evaluates the Nakagami CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) gamma_inc
  real ( kind = 8 ) p2
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) y

  if ( x <= 0.0D+00 ) then

    cdf = 0.0D+00

  else if ( 0.0D+00 < x ) then

    y = ( x - a ) / b
    x2 = c * y * y
    p2 = c

    cdf = gamma_inc ( p2, x2 )

  end if

  return
end
function nakagami_check ( a, b, c )

!*****************************************************************************80
!
!! NAKAGAMI_CHECK checks the parameters of the Nakagami PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, logical NAKAGAMI_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical nakagami_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NAKAGAMI_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    nakagami_check = .false.
    return
  end if

  if ( c <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NAKAGAMI_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C <= 0.'
    nakagami_check = .false.
    return
  end if

  nakagami_check = .true.

  return
end
subroutine nakagami_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! NAKAGAMI_MEAN returns the mean of the Nakagami PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B
!    0.0D+00 < C
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) mean
  real ( kind = 8 ) r8_gamma

  mean = a + b * r8_gamma ( c + 0.5D+00 ) / ( sqrt ( c ) * r8_gamma ( c ) )

  return
end
subroutine nakagami_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! NAKAGAMI_PDF evaluates the Nakagami PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) pdf
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x <= 0.0D+00 ) then

    pdf = 0.0D+00

  else if ( 0.0D+00 < x ) then

    y = ( x - a ) / b

    pdf = 2.0D+00 * c**c / ( b * r8_gamma ( c ) ) &
      * y**( 2.0D+00 * c - 1.0D+00 ) &
      * exp ( - c * y * y )

  end if

  return
end
subroutine nakagami_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! NAKAGAMI_VARIANCE returns the variance of the Nakagami PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B
!    0.0D+00 < C
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) variance

  t1 = r8_gamma ( c + 0.5D+00 )
  t2 = r8_gamma ( c )

  variance = b * b * ( 1.0D+00 - t1 * t1 / ( c * t2 * t2 ) )

  return
end
subroutine negative_binomial_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! NEGATIVE_BINOMIAL_CDF evaluates the Negative Binomial CDF.
!
!  Discussion:
!
!    A simple summing approach is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the CDF.
!
!    Input, integer ( kind = 4 ) A, a parameter of the PDF.
!    0 <= A.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    0 < B <= 1.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) cnk
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y

  cdf = 0.0D+00

  do y = a, x

    call binomial_coef ( y - 1, a - 1, cnk )

    pdf = real ( cnk, kind = 8 ) * b**a * ( 1.0D+00 - b )**( y - a )

    cdf = cdf + pdf

  end do

  return
end
subroutine negative_binomial_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! NEGATIVE_BINOMIAL_CDF_INV inverts the Negative Binomial CDF.
!
!  Discussion:
!
!    A simple discrete approach is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, integer ( kind = 4 ) A, a parameter of the PDF.
!    0 <= A.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    0 < B <= 1.
!
!    Output, integer ( kind = 4 ) X, the smallest X whose cumulative density
!    function is greater than or equal to CDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cum
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  integer ( kind = 4 ), parameter :: x_max = 1000

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NEGATIVE_BINOMIAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if


  cum = 0.0D+00

  x = a

  do

    call negative_binomial_pdf ( x, a, b, pdf )

    cum = cum + pdf

    if ( cdf <= cum .or. x_max <= x ) then
      exit
    end if

    x = x + 1

  end do

  return
end
subroutine negative_binomial_cdf_values ( n_data, f, s, p, cdf )

!*****************************************************************************80
!
!! NEGATIVE_BINOMIAL_CDF_VALUES returns values of the negative binomial CDF.
!
!  Discussion:
!
!    Assume that a coin has a probability P of coming up heads on
!    any one trial.  Suppose that we plan to flip the coin until we
!    achieve a total of S heads.  If we let F represent the number of
!    tails that occur in this process, then the value of F satisfies
!    a negative binomial PDF:
!
!      PDF(F,S,P) = Choose ( F from F+S-1 ) * P**S * (1-P)**F
!
!    The negative binomial CDF is the probability that there are F or
!    fewer failures upon the attainment of the S-th success.  Thus,
!
!      CDF(F,S,P) = sum ( 0 <= G <= F ) PDF(G,S,P)
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`DiscreteDistributions`]
!      dist = NegativeBinomialDistribution [ s, p ]
!      CDF [ dist, f ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    FC Powell,
!    Statistical Tables for Sociology, Biology and Physical Sciences,
!    Cambridge University Press, 1982.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) F, the maximum number of failures.
!
!    Output, integer ( kind = 4 ) S, the number of successes.
!
!    Output, real ( kind = 8 ) P, the probability of a success on one trial.
!
!    Output, real ( kind = 8 ) CDF, the probability of at most F failures
!    before the S-th success.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 27

  real ( kind = 8 ) cdf
  real ( kind = 8 ), save, dimension ( n_max ) :: cdf_vec = (/ &
    0.6367187500000000D+00, &
    0.3632812500000000D+00, &
    0.1445312500000000D+00, &
    0.5000000000000000D+00, &
    0.2265625000000000D+00, &
    0.6250000000000000D-01, &
    0.3437500000000000D+00, &
    0.1093750000000000D+00, &
    0.1562500000000000D-01, &
    0.1792000000000000D+00, &
    0.4096000000000000D-01, &
    0.4096000000000000D-02, &
    0.7047000000000000D-01, &
    0.1093500000000000D-01, &
    0.7290000000000000D-03, &
    0.9861587127990000D+00, &
    0.9149749500510000D+00, &
    0.7471846521450000D+00, &
    0.8499053647030009D+00, &
    0.5497160941090026D+00, &
    0.2662040052146710D+00, &
    0.6513215599000000D+00, &
    0.2639010709000000D+00, &
    0.7019082640000000D-01, &
    0.1000000000000000D+01, &
    0.1990000000000000D-01, &
    0.1000000000000000D-03 /)
  integer ( kind = 4 ) f
  integer ( kind = 4 ), save, dimension ( n_max ) :: f_vec = (/ &
     4,  3,  2, &
     3,  2,  1, &
     2,  1,  0, &
     2,  1,  0, &
     2,  1,  0, &
    11, 10,  9, &
    17, 16, 15, &
     9,  8,  7, &
     2,  1,  0 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ), save, dimension ( n_max ) :: p_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.40D+00, &
    0.40D+00, &
    0.40D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.10D-01, &
    0.10D-01, &
    0.10D-01 /)
  integer ( kind = 4 ) s
  integer ( kind = 4 ), save, dimension ( n_max ) :: s_vec = (/ &
    4, 5, 6, &
    4, 5, 6, &
    4, 5, 6, &
    4, 5, 6, &
    4, 5, 6, &
    1, 2, 3, &
    1, 2, 3, &
    1, 2, 3, &
    0, 1, 2 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    f = 0
    s = 0
    p = 0.0D+00
    cdf = 0.0D+00
  else
    f = f_vec(n_data)
    s = s_vec(n_data)
    p = p_vec(n_data)
    cdf = cdf_vec(n_data)
  end if

  return
end
function negative_binomial_check ( a, b )

!*****************************************************************************80
!
!! NEGATIVE_BINOMIAL_CHECK checks the parameters of the Negative Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, a parameter of the PDF.
!    0 <= A.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    0 < B <= 1.
!
!    Output, logical NEGATIVE_BINOMIAL_CHECK, is true if the
!    parameters are legal.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  logical negative_binomial_check

  if ( a < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NEGATIVE_BINOMIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A < 0.'
    negative_binomial_check = .false.
    return
  end if

  if ( b <= 0.0D+00 .or. 1.0D+00 < b ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NEGATIVE_BINOMIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0 or 1 < B.'
    negative_binomial_check = .false.
    return
  end if

  negative_binomial_check = .true.

  return
end
subroutine negative_binomial_mean ( a, b, mean )

!*****************************************************************************80
!
!! NEGATIVE_BINOMIAL_MEAN returns the mean of the Negative Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, a parameter of the PDF.
!    0 <= A.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    0 < B <= 1.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = real ( a, kind = 8 ) / b

  return
end
subroutine negative_binomial_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! NEGATIVE_BINOMIAL_PDF evaluates the Negative Binomial PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = C(X-1,A-1) * B**A * ( 1 - B )**(X-A)
!
!    PDF(A,B;X) is the probability that the A-th success will
!    occur on the X-th trial, given that the probability
!    of a success on a single trial is B.
!
!    The Negative Binomial PDF is also known as the Pascal PDF or
!    the "Polya" PDF.
!
!    NEGATIVE_BINOMIAL_PDF(1,B;X) = GEOMETRIC_PDF(B;X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the number of trials.
!    A <= X.
!
!    Input, integer ( kind = 4 ) A, the number of successes required.
!    0 <= A <= X, normally.
!
!    Input, real ( kind = 8 ) B, the probability of a success on a single trial.
!    0.0 < B <= 1.0.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) cnk
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  if ( x < a ) then

    pdf = 0.0D+00

  else

    call binomial_coef ( x - 1, a - 1, cnk )

    pdf = real ( cnk, kind = 8 ) * b**a * ( 1.0D+00 - b )**( x - a )

  end if

  return
end
subroutine negative_binomial_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! NEGATIVE_BINOMIAL_SAMPLE samples the Negative Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, a parameter of the PDF.
!    0 <= A.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    0 < B <= 1.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) num_success
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x

  if ( b == 1.0D+00 ) then
    x = a
    return
  else if ( b == 0.0D+00 ) then
    x = i4_huge ( )
    return
  end if

  x = 0
  num_success = 0

  do while ( num_success < a )

    x = x + 1
    r = r8_uniform_01 ( seed )

    if ( r <= b ) then
      num_success = num_success + 1
    end if

  end do

  return
end
subroutine negative_binomial_variance ( a, b, variance )

!*****************************************************************************80
!
!! NEGATIVE_BINOMIAL_VARIANCE returns the variance of the Negative Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, a parameter of the PDF.
!    0 <= A.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    0 < B <= 1.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = real ( a, kind = 8 ) * ( 1.0D+00 - b ) / ( b * b )

  return
end
subroutine normal_01_cdf ( x, cdf )

!*****************************************************************************80
!
!! NORMAL_01_CDF evaluates the Normal 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39,
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, pages 197-198, 1969.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
  real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
  real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
  real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
  real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
  real ( kind = 8 ), parameter :: b1 = 3.8052D-08
  real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
  real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
  real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
  real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
  real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
  real ( kind = 8 ) cdf
  real ( kind = 8 ) q
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  |X| <= 1.28.
!
  if ( abs ( x ) <= 1.28D+00 ) then

    y = 0.5D+00 * x * x

    q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
!
!  1.28 < |X| <= 12.7
!
  else if ( abs ( x ) <= 12.7D+00 ) then

    y = 0.5D+00 * x * x

    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
      + b2 / ( abs ( x ) + b3 &
      + b4 / ( abs ( x ) - b5 &
      + b6 / ( abs ( x ) + b7 &
      - b8 / ( abs ( x ) + b9 &
      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!
!  12.7 < |X|
!
  else

    q = 0.0D+00

  end if
!
!  Take account of negative X.
!
  if ( x < 0.0D+00 ) then
    cdf = q
  else
    cdf = 1.0D+00 - q
  end if

  return
end
subroutine normal_01_cdf_inv ( p, x )

!*****************************************************************************80
!
!! NORMAL_01_CDF_INV inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10**16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 2007
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    Algorithm AS241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the value of the cumulative probability
!    densitity function.  0 < P < 1.  If P is outside this range, an
!    "infinite" value will be returned.
!
!    Output, real ( kind = 8 ) X, the normal deviate value
!    with the property that the probability of a standard normal deviate being
!    less than or equal to the value is P.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
    3.3871328727963666080D+00, &
    1.3314166789178437745D+02, &
    1.9715909503065514427D+03, &
    1.3731693765509461125D+04, &
    4.5921953931549871457D+04, &
    6.7265770927008700853D+04, &
    3.3430575583588128105D+04, &
    2.5090809287301226727D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: b = (/ &
    1.0D+00, &
    4.2313330701600911252D+01, &
    6.8718700749205790830D+02, &
    5.3941960214247511077D+03, &
    2.1213794301586595867D+04, &
    3.9307895800092710610D+04, &
    2.8729085735721942674D+04, &
    5.2264952788528545610D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: c = (/ &
    1.42343711074968357734D+00, &
    4.63033784615654529590D+00, &
    5.76949722146069140550D+00, &
    3.64784832476320460504D+00, &
    1.27045825245236838258D+00, &
    2.41780725177450611770D-01, &
    2.27238449892691845833D-02, &
    7.74545014278341407640D-04 /)
  real ( kind = 8 ), parameter :: const1 = 0.180625D+00
  real ( kind = 8 ), parameter :: const2 = 1.6D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    1.0D+00, &
    2.05319162663775882187D+00, &
    1.67638483018380384940D+00, &
    6.89767334985100004550D-01, &
    1.48103976427480074590D-01, &
    1.51986665636164571966D-02, &
    5.47593808499534494600D-04, &
    1.05075007164441684324D-09 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: e = (/ &
    6.65790464350110377720D+00, &
    5.46378491116411436990D+00, &
    1.78482653991729133580D+00, &
    2.96560571828504891230D-01, &
    2.65321895265761230930D-02, &
    1.24266094738807843860D-03, &
    2.71155556874348757815D-05, &
    2.01033439929228813265D-07 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: f = (/ &
    1.0D+00, &
    5.99832206555887937690D-01, &
    1.36929880922735805310D-01, &
    1.48753612908506148525D-02, &
    7.86869131145613259100D-04, &
    1.84631831751005468180D-05, &
    1.42151175831644588870D-07, &
    2.04426310338993978564D-15 /)
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) r8poly_value
  real ( kind = 8 ), parameter :: split1 = 0.425D+00
  real ( kind = 8 ), parameter :: split2 = 5.0D+00
  real ( kind = 8 ) x

  if ( p <= 0.0D+00 ) then
    x = - huge ( x )
    return
  end if

  if ( 1.0D+00 <= p ) then
    x = huge ( x )
    return
  end if

  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    x = q * r8poly_value ( 8, a, r ) / r8poly_value ( 8, b, r )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

    if ( r <= 0.0D+00 ) then

      x = huge ( x )

    else

      r = sqrt ( - log ( r ) )

      if ( r <= split2 ) then

        r = r - const2
        x = r8poly_value ( 8, c, r ) / r8poly_value ( 8, d, r )

      else

        r = r - split2
        x = r8poly_value ( 8, e, r ) / r8poly_value ( 8, f, r )

      end if

    end if

    if ( q < 0.0D+00 ) then
      x = -x
    end if

  end if

  return
end
subroutine normal_01_cdf_values ( n_data, x, fx )

!*****************************************************************************80
!
!! NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NormalDistribution [ 0, 1 ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.5398278372770290D+00, &
    0.5792597094391030D+00, &
    0.6179114221889526D+00, &
    0.6554217416103242D+00, &
    0.6914624612740131D+00, &
    0.7257468822499270D+00, &
    0.7580363477769270D+00, &
    0.7881446014166033D+00, &
    0.8159398746532405D+00, &
    0.8413447460685429D+00, &
    0.9331927987311419D+00, &
    0.9772498680518208D+00, &
    0.9937903346742239D+00, &
    0.9986501019683699D+00, &
    0.9997673709209645D+00, &
    0.9999683287581669D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000D+00, &
    0.1000000000000000D+00, &
    0.2000000000000000D+00, &
    0.3000000000000000D+00, &
    0.4000000000000000D+00, &
    0.5000000000000000D+00, &
    0.6000000000000000D+00, &
    0.7000000000000000D+00, &
    0.8000000000000000D+00, &
    0.9000000000000000D+00, &
    0.1000000000000000D+01, &
    0.1500000000000000D+01, &
    0.2000000000000000D+01, &
    0.2500000000000000D+01, &
    0.3000000000000000D+01, &
    0.3500000000000000D+01, &
    0.4000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine normal_01_mean ( mean )

!*****************************************************************************80
!
!! NORMAL_01_MEAN returns the mean of the Normal 01 PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) mean

  mean = 0.0D+00

  return
end
subroutine normal_01_pdf ( x, pdf )

!*****************************************************************************80
!
!! NORMAL_01_PDF evaluates the Normal 01 PDF.
!
!  Discussion:
!
!    The Normal 01 PDF is also called the "Standard Normal" PDF, or
!    the Normal PDF with 0 mean and variance 1.
!
!    PDF(X) = exp ( - 0.5 * X^2 ) / sqrt ( 2 * PI )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  pdf = exp ( -0.5D+00 * x * x ) / sqrt ( 2.0D+00 * pi )

  return
end
subroutine normal_01_sample ( seed, x )

!*****************************************************************************80
!
!! NORMAL_01_SAMPLE samples the standard normal probability distribution.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    The Box-Muller method is used, which is efficient, but
!    generates two values at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2002
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
!    Output, real ( kind = 8 ) X, a sample of the standard normal PDF.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: used = -1
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: y = 0.0D+00

  if ( used == -1 ) then
    used = 0
  end if
!
!  If we've used an even number of values so far, generate two more,
!  return one and save one.
!
  if ( mod ( used, 2 ) == 0 ) then

    do

      r1 = r8_uniform_01 ( seed )

      if ( r1 /= 0.0D+00 ) then
        exit
      end if

    end do

    r2 = r8_uniform_01 ( seed )

    x = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    y = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
!
!  Otherwise, return the second, saved, value.
!
  else

    x = y

  end if

  used = used + 1

  return
end
subroutine normal_01_variance ( variance )

!*****************************************************************************80
!
!! NORMAL_01_VARIANCE returns the variance of the Normal 01 PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) variance

  variance = 1.0D+00

  return
end
subroutine normal_01_vector ( n, seed, x )

!*****************************************************************************80
!
!! NORMAL_01_VECTOR samples the standard normal probability distribution.
!
!  Discussion:
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
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is
!    negative, then the code will flush its internal memory; in particular,
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
!    Local, integer MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real R(N+1), is used to store some uniform random values.
!    Its dimension is N+1, but really it is only needed to be the
!    smallest even number greater than or equal to N.
!
!    Local, integer SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
!    Local, real Y, the value saved from the previous call, if
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

    call r8vec_uniform_01 ( 2 * m, seed, r )

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

    call r8vec_uniform_01 ( 2 * m, seed, r )

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
subroutine normal_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! NORMAL_CDF evaluates the Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  y = ( x - a ) / b

  call normal_01_cdf ( y, cdf )

  return
end
subroutine normal_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! NORMAL_CDF_INV inverts the Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NORMAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  call normal_01_cdf_inv ( cdf, x2 )

  x = a + b * x2

  return
end
subroutine normal_cdf_values ( n_data, mu, sigma, x, fx )

!*****************************************************************************80
!
!! NORMAL_CDF_VALUES returns some values of the Normal CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NormalDistribution [ mu, sigma ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) MU, the mean of the distribution.
!
!    Output, real ( kind = 8 ) SIGMA, the variance of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.9772498680518208D+00, &
    0.9999683287581669D+00, &
    0.9999999990134124D+00, &
    0.6914624612740131D+00, &
    0.6305586598182364D+00, &
    0.5987063256829237D+00, &
    0.5792597094391030D+00, &
    0.6914624612740131D+00, &
    0.5000000000000000D+00, &
    0.3085375387259869D+00, &
    0.1586552539314571D+00 /)
  real ( kind = 8 ) mu
  real ( kind = 8 ), save, dimension ( n_max ) :: mu_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ), save, dimension ( n_max ) :: sigma_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    mu = 0.0D+00
    sigma = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    mu = mu_vec(n_data)
    sigma = sigma_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function normal_check ( a, b )

!*****************************************************************************80
!
!! NORMAL_CHECK checks the parameters of the Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, logical NORMAL_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical normal_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NORMAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    normal_check = .false.
    return
  end if

  normal_check = .true.

  return
end
subroutine normal_mean ( a, b, mean )

!*****************************************************************************80
!
!! NORMAL_MEAN returns the mean of the Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine normal_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! NORMAL_PDF evaluates the Normal PDF.
!
!  Discussion:
!
!    PDF(A,B;X)
!      = exp ( - 0.5D+00 * ( ( X - A ) / B )^2 ) / ( B * sqrt ( 2 * PI ) )
!
!    The normal PDF is also known as the Gaussian PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  y = ( x - a ) / b

  pdf = exp ( - 0.5D+00 * y * y )  / ( b * sqrt ( 2.0D+00 * pi ) )

  return
end
subroutine normal_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! NORMAL_SAMPLE samples the Normal PDF.
!
!  Discussion:
!
!    The Box-Muller method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  call normal_01_sample ( seed, x )

  x = a + b * x

  return
end
subroutine normal_variance ( a, b, variance )

!*****************************************************************************80
!
!! NORMAL_VARIANCE returns the variance of the Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = b * b

  return
end
subroutine normal_vector ( n, mean, dev, seed, x )

!*****************************************************************************80
!
!! NORMAL_VECTOR samples the normal probability distribution.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) has
!    a user-specified mean and standard deviation.
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
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is
!    negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input, real ( kind = 8 ) MEAN, the desired mean value.
!
!    Input, real ( kind = 8 ) DEV, the desired standard deviation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) dev
  real ( kind = 8 ) mean

  call normal_01_vector ( n, seed, x )

  x(1:n) = mean + dev * x(1:n)

  return
end
subroutine owen_values ( n_data, h, a, t )

!*****************************************************************************80
!
!! OWEN_VALUES returns some values of Owen's T function.
!
!  Discussion:
!
!    Owen's T function is useful for computation of the bivariate normal
!    distribution and the distribution of a skewed normal distribution.
!
!    Although it was originally formulated in terms of the bivariate
!    normal function, the function can be defined more directly as
!
!      T(H,A) = 1 / ( 2 * pi ) *
!        Integral ( 0 <= X <= A ) e^(-H^2*(1+X^2)/2) / (1+X^2) dX
!
!    In Mathematica, the function can be evaluated by:
!
!      fx = 1/(2*Pi) * Integrate [ E^(-h^2*(1+x^2)/2)/(1+x^2), {x,0,a} ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) H, a parameter.
!
!    Output, real ( kind = 8 ) A, the upper limit of the integral.
!
!    Output, real ( kind = 8 ) T, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 22

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.1000000000000000D+02, &
    0.1000000000000000D+03 /)
  real ( kind = 8 ) h
  real ( kind = 8 ), save, dimension ( n_max ) :: h_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) t
  real ( kind = 8 ), save, dimension ( n_max ) :: t_vec = (/ &
    0.4306469112078537D-01, &
    0.6674188216570097D-01, &
    0.7846818699308410D-01, &
    0.7929950474887259D-01, &
    0.6448860284750376D-01, &
    0.1066710629614485D+00, &
    0.1415806036539784D+00, &
    0.1510840430760184D+00, &
    0.7134663382271778D-01, &
    0.1201285306350883D+00, &
    0.1666128410939293D+00, &
    0.1847501847929859D+00, &
    0.7317273327500385D-01, &
    0.1237630544953746D+00, &
    0.1737438887583106D+00, &
    0.1951190307092811D+00, &
    0.7378938035365546D-01, &
    0.1249951430754052D+00, &
    0.1761984774738108D+00, &
    0.1987772386442824D+00, &
    0.2340886964802671D+00, &
    0.2479460829231492D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    h = 0.0D+00
    a = 0.0D+00
    t = 0.0D+00
  else
    h = h_vec(n_data)
    a = a_vec(n_data)
    t = t_vec(n_data)
  end if

  return
end
subroutine pareto_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! PARETO_CDF evaluates the Pareto CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x < a ) then
    cdf = 0.0D+00
  else
    cdf = 1.0D+00 - ( a / x )**b
  end if

  return
end
subroutine pareto_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! PARETO_CDF_INV inverts the Pareto CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARETO_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a / ( 1.0D+00 - cdf )**( 1.0D+00 / b )

  return
end
function pareto_check ( a, b )

!*****************************************************************************80
!
!! PARETO_CHECK checks the parameters of the Pareto CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, logical PARETO_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical pareto_check

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARETO_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    pareto_check = .false.
    return
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARETO_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    pareto_check = .false.
    return
  end if

  pareto_check = .true.

  return
end
subroutine pareto_mean ( a, b, mean )

!*****************************************************************************80
!
!! PARETO_MEAN returns the mean of the Pareto PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  if ( b <= 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARETO_MEAN - Fatal error!'
    write ( *, '(a)' ) '  For B <= 1, the mean does not exist.'
    mean = 0.0D+00
    return
  end if

  mean = b * a / ( b - 1.0D+00 )

  return
end
subroutine pareto_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! PARETO_PDF evaluates the Pareto PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = B * A**B / X**(B+1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A <= X
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x < a ) then
    pdf = 0.0D+00
  else
    pdf = b * a**b / x**( b + 1.0D+00 )
  end if

  return
end
subroutine pareto_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! PARETO_SAMPLE samples the Pareto PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call pareto_cdf_inv ( cdf, a, b, x )

  return
end
subroutine pareto_variance ( a, b, variance )

!*****************************************************************************80
!
!! PARETO_VARIANCE returns the variance of the Pareto PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  if ( b <= 2.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARETO_VARIANCE - Warning!'
    write ( *, '(a)' ) '  For B <= 2, the variance does not exist.'
    variance = 0.0D+00
    return
  end if

  variance = a * a * b / ( ( b - 1.0D+00 )**2 * ( b - 2.0D+00 ) )

  return
end
function pearson_05_check ( a, b, c )

!*****************************************************************************80
!
!! PEARSON_05_CHECK checks the parameters of the Pearson 5 PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < A, 0.0D+00 < B.
!
!    Output, logical PEARSON_05_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical pearson_05_check

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PEARSON_05_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    pearson_05_check = .false.
    return
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PEARSON_05_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    pearson_05_check = .false.
    return
  end if

  pearson_05_check = .true.

  return
end
subroutine pearson_05_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! PEARSON_05_MEAN evaluates the mean of the Pearson 5 PDF.
!
!  Discussion:
!
!    The mean is undefined for B <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < A, 0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) mean

  if ( b <= 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PEARSON_05_MEAN - Warning!'
    write ( *, '(a)' ) '  MEAN undefined for B <= 1.'
    mean = c
    return
  end if

  mean = c + a / ( b - 1.0D+00 )

  return
end
subroutine pearson_05_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! PEARSON_05_PDF evaluates the Pearson 5 PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = A**B * ( X - C )**(-B-1)
!      * exp ( - A / ( X - C ) ) / Gamma ( B )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    C < X
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < A, 0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) pdf
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) x

  if ( x <= c ) then
    pdf = 0.0D+00
  else
    pdf = a**b * ( x - c )**( - b - 1.0D+00 ) &
      * exp ( - a / ( x - c ) ) / r8_gamma ( b )
  end if

  return
end
subroutine pearson_05_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! PEARSON_05_SAMPLE samples the Pearson 5 PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < A, 0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) c
  real ( kind = 8 ) c2
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  a2 = 0.0D+00
  b2 = b
  c2 = 1.0D+00 / a

  call gamma_sample ( a2, b2, c2, seed, x2 )

  x = c + 1.0D+00 / x2

  return
end
function planck_check ( a, b )

!*****************************************************************************80
!
!! PLANCK_CHECK checks the parameters of the Planck PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, logical PLANCK_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical planck_check

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANCK_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    planck_check = .false.
    return
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANCK_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    planck_check = .false.
    return
  end if

  planck_check = .true.

  return
end
subroutine planck_mean ( a, b, mean )

!*****************************************************************************80
!
!! PLANCK_MEAN returns the mean of the Planck PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean
  real ( kind = 8 ) zeta

  mean = ( b + 1.0D+00 ) * zeta ( b + 2.0D+00 ) / zeta ( b + 1.0D+00 )

  return
end
subroutine planck_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! PLANCK_PDF evaluates the Planck PDF.
!
!  Discussion:
!
!    The Planck PDF has the form
!
!      PDF(A,B;X) = A**(B+1) * X**B / ( exp ( A * X ) - 1 ) / K
!
!    where K is the normalization constant, and has the value
!
!      K = Gamma ( B + 1 ) * Zeta ( B + 1 ).
!
!    The original Planck distribution governed the frequencies in
!    blackbody radiation at a given temperature T, and has the form
!
!      PDF(A;X) = K * X**3 / ( exp ( A * X ) - 1 )
!
!    with
!
!      K = 15 / PI**4.
!
!    Thus, in terms of the Planck PDF, the original Planck distribution
!    has A = 1, B = 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Norman Johnson, Samuel Kotz, Balakrishnan,
!    Continuous Univariate Distributions, second edition,
!    Wiley, 1994,
!    QA273.6.J6
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0 <= X
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) k
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) x
  real ( kind = 8 ) zeta

  if ( x < 0.0D+00 ) then
    pdf = 0.0D+00
  else
    k = r8_gamma ( b + 1.0D+00 ) * zeta ( b + 1.0D+00 )
    pdf = a**( b + 1.0D+00 ) * x**b / ( exp ( a * x ) - 1.0D+00 ) / k
  end if

  return
end
subroutine planck_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! PLANCK_SAMPLE samples the Planck PDF.
!
!  Discussion:
!
!    The Planck sampling seems to be giving incorrect results.
!    I suspect this has to do with a possible problem in the
!    ZIPF_SAMPLE routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Luc Devroye,
!    Non-Uniform Random Variate Generation,
!    Springer Verlag, 1986, pages 552.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) c2
  real ( kind = 8 ) g
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  integer ( kind = 4 ) z

  a2 = 0.0D+00
  b2 = 1.0D+00
  c2 = b + 1.0D+00

  call gamma_sample ( a2, b2, c2, seed, g )

  call zipf_sample ( c2, seed, z )

  x = g / ( a * real ( z, kind = 8 ) )

  return
end
subroutine planck_variance ( a, b, variance )

!*****************************************************************************80
!
!! PLANCK_VARIANCE returns the variance of the Planck PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance
  real ( kind = 8 ) zeta

  call planck_mean ( a, b, mean )

  variance = ( b + 1.0D+00 ) * ( b + 2.0D+00 ) &
    * zeta ( b + 3.0D+00 ) / zeta ( b + 1.0D+00 ) - mean * mean

  return
end
subroutine point_distance_1d_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! POINT_DISTANCE_1D_PDF evaluates the point distance PDF in 1D.
!
!  Discussion:
!
!    It is assumed that a set of points has been generated in 1D
!    according to a Poisson process.  The number of points in a region
!    of size LENGTH is a Poisson variate with mean value B * LENGTH.
!
!    For a point chosen at random, we may now find the nearest
!    Poisson point, the second nearest and so on.  We are interested
!    in the PDF that governs the expected behavior of the distances
!    of rank A = 1, 2, 3, ... with Poisson density B.
!
!    Note that this PDF is a form of the Gamma PDF.???
!
!    PDF(A,B;X) = B**A * X**( A - 1 ) * exp ( - B * X ) / ( A - 1 )!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X.
!
!    Input, integer ( kind = 4 ) A, indicates the degree of nearness of the
!    point.  A = 1 means the nearest point, A = 2 the second nearest, and so on.
!    0 < A.
!
!    Input, real ( kind = 8 ) B, the point density.  0.0 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) i4_factorial
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( a < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINT_DISTANCE_1D_PDF - Fatal error!'
    write ( *, '(a)' ) '  Input parameter A < 1.'
    stop
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINT_DISTANCE_1D_PDF - Fatal error!'
    write ( *, '(a)' ) '  Input parameter B <= 0.0.'
    stop
  end if

  if ( x < 0.0D+00 ) then
    pdf = 0.0D+00
  else
    pdf = b**a * x**( a - 1 ) * exp ( - b * x ) / i4_factorial ( a - 1 )
  end if

  return
end
subroutine point_distance_2d_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! POINT_DISTANCE_2D_PDF evaluates the point distance PDF in 2D.
!
!  Discussion:
!
!    It is assumed that a set of points has been generated in 2D
!    according to a Poisson process.  The number of points in a region
!    of size AREA is a Poisson variate with mean value B * AREA.
!
!    For a point chosen at random, we may now find the nearest
!    Poisson point, the second nearest and so on.  We are interested
!    in the PDF that governs the expected behavior of the distances
!    of rank A = 1, 2, 3, ... with Poisson density B.
!
!    PDF(A,B;X) = 2 * ( B * PI )**A * X**( 2 * A - 1 )
!      * EXP ( - B * PI * X * X ) / ( A - 1 )!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 579.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X.
!
!    Input, integer ( kind = 4 ) A, indicates the degree of nearness of the
!    point.  A = 1 means the nearest point, A = 2 the second nearest, and so on.
!    0 < A.
!
!    Input, real ( kind = 8 ) B, the point density.  0.0 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) i4_factorial
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( a < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINT_DISTANCE_2D_PDF - Fatal error!'
    write ( *, '(a)' ) '  Input parameter A < 1.'
    stop
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINT_DISTANCE_2D_PDF - Fatal error!'
    write ( *, '(a)' ) '  Input parameter B <= 0.0.'
    stop
  end if

  if ( x < 0.0D+00 ) then
    pdf = 0.0D+00
  else
    pdf = 2.0D+00 * ( b * pi )**a * x**( 2 * a - 1 ) &
      * exp ( - b * pi * x * x ) / i4_factorial ( a - 1 )
  end if

  return
end
subroutine point_distance_3d_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! POINT_DISTANCE_3D_PDF evaluates the point distance PDF in the 3D.
!
!  Discussion:
!
!    It is assumed that a set of points has been generated in 3D
!    according to a Poisson process.  The number of points in a region
!    of size VOLUME is a Poisson variate with mean value B * VOLUME.
!
!    For a point chosen at random, we may now find the nearest
!    Poisson point, the second nearest and so on.  We are interested
!    in the PDF that governs the expected behavior of the distances
!    of rank A = 1, 2, 3, ... with Poisson density B.
!
!    PDF(A,B;X) = 3 * ( (4/3) * B * PI )**A * X**( 3 * A - 1 )
!      * EXP ( - (4/3) * B * PI * X * X * X ) / ( A - 1 )!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 580.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X.
!
!    Input, integer ( kind = 4 ) A, indicates the degree of nearness of the
!    point.  A = 1 means the nearest point, A = 2 the second nearest, and so on.
!    0 < A.
!
!    Input, real ( kind = 8 ) B, the Poisson point density.  0.0 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) i4_factorial
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( a < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINT_DISTANCE_3D_PDF - Fatal error!'
    write ( *, '(a)' ) '  Input parameter A < 1.'
    stop
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINT_DISTANCE_3D_PDF - Fatal error!'
    write ( *, '(a)' ) '  Input parameter B <= 0.0.'
    stop
  end if

  if ( x < 0.0D+00 ) then
    pdf = 0.0D+00
  else
    pdf = 3.0D+00 * ( ( 4.0D+00 / 3.0D+00 ) * b * pi )**a &
      * x**( 3 * a - 1 ) * exp ( - ( 4.0D+00 / 3.0D+00 ) * b * pi * x**3 ) &
      / i4_factorial ( a - 1 )
  end if

  return
end
subroutine poisson_cdf ( x, a, cdf )

!*****************************************************************************80
!
!! POISSON_CDF evaluates the Poisson CDF.
!
!  Discussion:
!
!    CDF(X,A) is the probability that the number of events observed
!    in a unit time period will be no greater than X, given that the
!    expected number of events in a unit time period is A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the CDF.
!    0 <= X.
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) last
  real ( kind = 8 ) new
  real ( kind = 8 ) sum2
  integer ( kind = 4 ) x

  if ( x < 0 ) then

    cdf = 0.0D+00

  else

    new = exp ( - a )
    sum2 = new

    do i = 1, x
      last = new
      new = last * a / real ( i, kind = 8 )
      sum2 = sum2 + new
    end do

    cdf = sum2

  end if

  return
end
subroutine poisson_cdf_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! POISSON_CDF_VALUES returns some values of the Poisson CDF.
!
!  Discussion:
!
!    CDF(X)(A) is the probability of at most X successes in unit time,
!    given that the expected mean number of successes is A.
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`DiscreteDistributions`]
!      dist = PoissonDistribution [ a ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition, CRC Press, 1996, pages 653-658.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, the parameter of the function.
!
!    Output, integer ( kind = 4 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.02D+00, &
    0.10D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    1.00D+00, &
    2.00D+00, &
    2.00D+00, &
    2.00D+00, &
    2.00D+00, &
    5.00D+00, &
    5.00D+00, &
    5.00D+00, &
    5.00D+00, &
    5.00D+00, &
    5.00D+00, &
    5.00D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.9801986733067553D+00, &
    0.9048374180359596D+00, &
    0.9953211598395555D+00, &
    0.6065306597126334D+00, &
    0.9097959895689501D+00, &
    0.9856123220330293D+00, &
    0.3678794411714423D+00, &
    0.7357588823428846D+00, &
    0.9196986029286058D+00, &
    0.9810118431238462D+00, &
    0.1353352832366127D+00, &
    0.4060058497098381D+00, &
    0.6766764161830635D+00, &
    0.8571234604985470D+00, &
    0.6737946999085467D-02, &
    0.4042768199451280D-01, &
    0.1246520194830811D+00, &
    0.2650259152973617D+00, &
    0.4404932850652124D+00, &
    0.6159606548330631D+00, &
    0.7621834629729387D+00 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) x
  integer ( kind = 4 ), save, dimension ( n_max ) :: x_vec = (/ &
     0, 0, 1, 0, &
     1, 2, 0, 1, &
     2, 3, 0, 1, &
     2, 3, 0, 1, &
     2, 3, 4, 5, &
     6 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine poisson_cdf_inv ( cdf, a, x )

!*****************************************************************************80
!
!! POISSON_CDF_INV inverts the Poisson CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, a value of the CDF.
!    0 <= CDF < 1.
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, integer ( kind = 4 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) last
  real ( kind = 8 ) new
  real ( kind = 8 ) sum2
  real ( kind = 8 ) sumold
  integer ( kind = 4 ) x
  integer ( kind = 4 ), parameter :: xmax = 100

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POISSON_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if
!
!  Now simply start at X = 0, and find the first value for which
!  CDF(X-1) <= CDF <= CDF(X).
!
  sum2 = 0.0D+00

  do i = 0, xmax

    sumold = sum2

    if ( i == 0 ) then
      new = exp ( - a )
      sum2 = new
    else
      last = new
      new = last * a / real ( i, kind = 8 )
      sum2 = sum2 + new
    end if

    if ( sumold <= cdf .and. cdf <= sum2 ) then
      x = i
      return
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POISSON_CDF_INV - Warning!'
  write ( *, '(a,i8)' ) '  Exceeded XMAX = ', xmax

  x = xmax

  return
end
function poisson_check ( a )

!*****************************************************************************80
!
!! POISSON_CHECK checks the parameter of the Poisson PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, logical POISSON_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  logical poisson_check

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POISSON_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    poisson_check = .false.
    return
  end if

  poisson_check = .true.

  return
end
subroutine poisson_mean ( a, mean )

!*****************************************************************************80
!
!! POISSON_MEAN returns the mean of the Poisson PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine poisson_kernel ( r, n, c, x, y, p )

!*****************************************************************************80
!
!! POISSON_KERNEL evaluates the Poisson kernel.
!
!  Discussion:
!
!    P(X,Y) = ( R^2 - |X-C|^2 ) / ( R * A * |X-Y|^N )
!
!    where the N-dimensional ball has radius R and center C,
!    and A is the area of the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the ball.
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) C(N), the center of the ball.
!
!    Input, real ( kind = 8 ) X(N), a point inside the ball.
!
!    Input, real ( kind = 8 ) Y(N), a point on the surface of the ball.
!
!    Output, real ( kind = 8 ) P, the Poisson kernel function P(X,Y).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) area
  real ( kind = 8 ) b
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) p
  real ( kind = 8 ) r
  real ( kind = 8 ) r8vec_diff_norm
  real ( kind = 8 ) sphere_unit_area_nd
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xc_diff_norm
  real ( kind = 8 ) xy_diff_norm
  real ( kind = 8 ) y(n)

  xc_diff_norm = r8vec_diff_norm ( n, x, c )
  xy_diff_norm = r8vec_diff_norm ( n, x, y )
  area = sphere_unit_area_nd ( n )

  t = ( r + xc_diff_norm ) * ( r - xc_diff_norm )
  b = r * area * ( xy_diff_norm ) ** n
  p = t / b

  return
end
subroutine poisson_pdf ( x, a, pdf )

!*****************************************************************************80
!
!! POISSON_PDF evaluates the Poisson PDF.
!
!  Discussion:
!
!    PDF(A;X) = EXP ( - A ) * A**X / X!
!
!    PDF(A;X) is the probability that the number of events observed
!    in a unit time period will be X, given the expected number
!    of events in a unit time.
!
!    The parameter A is the expected number of events per unit time.
!
!    The Poisson PDF is a discrete version of the Exponential PDF.
!
!    The time interval between two Poisson events is a random
!    variable with the Exponential PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the PDF.
!    0 <= X
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) i4_factorial
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    pdf = 0.0D+00
  else
    pdf = exp ( - a ) * a**x / i4_factorial ( x )
  end if

  return
end
subroutine poisson_sample ( a, seed, x )

!*****************************************************************************80
!
!! POISSON_SAMPLE samples the Poisson PDF.
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
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x

  cdf = r8_uniform_01 ( seed )

  call poisson_cdf_inv ( cdf, a, x )

  return
end
subroutine poisson_variance ( a, variance )

!*****************************************************************************80
!
!! POISSON_VARIANCE returns the variance of the Poisson PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) variance

  variance = a

  return
end
subroutine power_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! POWER_CDF evaluates the Power CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A, 0.0D+00 < B,
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then
    cdf = 0.0D+00
  else if ( x <= b ) then
    cdf = ( x / b )**a
  else
    cdf = 1.0D+00
  end if

  return
end
subroutine power_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! POWER_CDF_INV inverts the Power CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A, 0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POWER_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = 0.0D+00
  else if ( cdf < 1.0D+00 ) then
    x = b * exp ( log ( cdf ) / a )
  else
    x = b
  end if

  return
end
function power_check ( a, b )

!*****************************************************************************80
!
!! POWER_CHECK checks the parameter of the Power PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A, 0.0D+00 < B.
!
!    Output, logical POWER_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical power_check

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POWER_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    power_check = .false.
    return
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POWER_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    power_check = .false.
    return
  end if

  power_check = .true.

  return
end
subroutine power_mean ( a, b, mean )

!*****************************************************************************80
!
!! POWER_MEAN returns the mean of the Power PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A, 0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a * b / ( a + 1.0D+00 )

  return
end
subroutine power_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! POWER_PDF evaluates the Power PDF.
!
!  Discussion:
!
!    PDF(A;X) = (A/B) * (X/B)**(A-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, Stephen Kokoska,
!    CRC Standard Probability and Statistics Tables and Formulae,
!    Chapman and Hall/CRC, 2000, pages 152-153.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X <= B.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A, 0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x < 0.0D+00 .or. b < x ) then
    pdf = 0.0D+00
  else
    pdf = ( a / b ) * ( x / b )**( a - 1.0D+00 )
  end if

  return
end
subroutine power_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! POWER_SAMPLE samples the Power PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A, 0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call power_cdf_inv ( cdf, a, b, x )

  return
end
subroutine power_variance ( a, b, variance )

!*****************************************************************************80
!
!! POWER_VARIANCE returns the variance of the Power PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A, 0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = b * b * a / ( ( a + 1.0D+00 )**2 * ( a + 2.0D+00 ) )

  return
end
subroutine psi_values ( n_data, x, fx )

!*****************************************************************************80
!
!! PSI_VALUES returns some values of the Psi or Digamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      PolyGamma[x]
!
!    or
!
!      PolyGamma[0,x]
!
!    PSI(X) = d ln ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
!
!    PSI(1) = -Euler's constant.
!
!    PSI(X+1) = PSI(X) + 1 / X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.5772156649015329D+00, &
    -0.4237549404110768D+00, &
    -0.2890398965921883D+00, &
    -0.1691908888667997D+00, &
    -0.6138454458511615D-01, &
     0.3648997397857652D-01, &
     0.1260474527734763D+00, &
     0.2085478748734940D+00, &
     0.2849914332938615D+00, &
     0.3561841611640597D+00, &
     0.4227843350984671D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine quasigeometric_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! QUASIGEOMETRIC_CDF evaluates the Quasigeometric CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the maximum number of trials.
!
!    Input, real ( kind = 8 ) A, the probability of 0 successes.
!    0.0 <= A <= 1.0.
!
!    Input, real ( kind = 8 ) B, the depreciation constant.
!    0.0 <= B < 1.0.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    cdf = 0.0D+00
  else if ( x == 0 ) then
    cdf = a
  else if ( b == 0.0D+00 ) then
    cdf = 1.0D+00
  else
    cdf = a + ( 1.0D+00 - a ) * ( 1.0D+00 - b**x )
  end if

  return
end
subroutine quasigeometric_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! QUASIGEOMETRIC_CDF_INV inverts the Quasigeometric CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0D+00
!
!    Input, real ( kind = 8 ) A, the probability of 0 successes.
!    0.0 <= A <= 1.0.
!
!    Input, real ( kind = 8 ) B, the depreciation constant.
!    0.0 <= B < 1.0.
!
!    Output, integer ( kind = 4 ) X, the corresponding value of X.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUASIGEOMETRIC_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf < a ) then
    x = 0
  else if ( b == 0.0D+00 ) then
    x = 1
  else
    x = 1 + int ( ( log ( 1.0D+00 - cdf ) - log ( 1.0D+00 - a ) ) / log ( b ) )
  end if

  return
end
function quasigeometric_check ( a, b )

!*****************************************************************************80
!
!! QUASIGEOMETRIC_CHECK checks the parameters of the Quasigeometric CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the probability of 0 successes.
!    0.0 <= A <= 1.0.
!
!    Input, real ( kind = 8 ) B, the depreciation constant.
!    0.0 <= B < 1.0.
!
!    Output, logical QUASIGEOMETRIC_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical              quasigeometric_check

  if ( a < 0.0D+00 .or. 1.0D+00 < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUASIGEOMETRIC_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A < 0 or 1 < A.'
    quasigeometric_check = .false.
    return
  end if

  if ( b < 0.0D+00 .or. 1.0D+00 <= b ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUASIGEOMETRIC_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B < 0 or 1 <= B.'
    quasigeometric_check = .false.
    return
  end if

  quasigeometric_check = .true.

  return
end
subroutine quasigeometric_mean ( a, b, mean )

!*****************************************************************************80
!
!! QUASIGEOMETRIC_MEAN returns the mean of the Quasigeometric PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the probability of 0 successes.
!    0.0 <= A <= 1.0.
!
!    Input, real ( kind = 8 ) B, the depreciation constant.
!    0.0 <= B < 1.0.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = ( 1.0D+00 - a  ) / ( 1.0D+00 - b )

  return
end
subroutine quasigeometric_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! QUASIGEOMETRIC_PDF evaluates the Quasigeometric PDF.
!
!  Discussion:
!
!    PDF(A,B;X) =    A                     if 0  = X;
!               = (1-A) * (1-B) * B^(X-1)  if 1 <= X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Darren Glass, Philip Lowry,
!    Quasiquasigeometric Distributions and Extra Inning Baseball Games,
!    Mathematics Magazine,
!    Volume 81, Number 2, April 2008, pages 127-137.
!
!    Paul Nahin,
!    Digital Dice: Computational Solutions to Practical Probability Problems,
!    Princeton University Press, 2008,
!    ISBN13: 978-0-691-12698-2,
!    LC: QA273.25.N34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the independent variable.
!    0 <= X
!
!    Input, real ( kind = 8 ) A, the probability of 0 successes.
!    0.0 <= A <= 1.0.
!
!    Input, real ( kind = 8 ) B, the depreciation constant.
!    0.0 <= B < 1.0.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  if ( x < 0 ) then

    pdf = 0.0D+00

  else if ( x == 0 ) then

    pdf = a

  else if ( b == 0.0D+00 ) then

    if ( x == 1 ) then
      pdf = 1.0D+00
    else
      pdf = 0.0D+00
    end if

  else

    pdf = ( 1.0D+00 - a ) * ( 1.0D+00 - b ) * b**( x - 1 )

  end if

  return
end
subroutine quasigeometric_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! QUASIGEOMETRIC_SAMPLE samples the Quasigeometric PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the probability of 0 successes.
!    0.0 <= A <= 1.0.
!
!    Input, real ( kind = 8 ) B, the depreciation constant.
!    0.0 <= B < 1.0.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x

  cdf = r8_uniform_01 ( seed )

  call quasigeometric_cdf_inv ( cdf, a, b, x )

  return
end
subroutine quasigeometric_variance ( a, b, variance )

!*****************************************************************************80
!
!! QUASIGEOMETRIC_VARIANCE returns the variance of the Quasigeometric PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the probability of 0 successes.
!    0.0 <= A <= 1.0.
!
!    Input, real ( kind = 8 ) B, the depreciation constant.
!    0.0 <= B < 1.0.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = ( 1.0D+00 - a ) * ( a + b ) / ( 1.0D+00 - b ) / ( 1.0D+00 - b )

  return
end
function r4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R4_UNIFORM returns a scaled real ( kind = 4 ) pseudorandom number.
!
!  Discussion:
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
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
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  real ( kind = 4 ) r4_uniform
  integer ( kind = 4 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  r4_uniform = a + ( b - a ) * real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit real ( kind = 4 ) pseudorandom number.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r4_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
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
!    11 August 2004
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
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r4_uniform_01

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
function r8_ceiling ( r )

!*****************************************************************************80
!
!! R8_CEILING rounds an R8 "up" to the nearest integer.
!
!  Example:
!
!    R     Value
!
!   -1.1  -1
!   -1.0  -1
!   -0.9   0
!    0.0   0
!    5.0   5
!    5.1   6
!    5.9   6
!    6.0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the real value to be rounded up.
!
!    Output, integer ( kind = 4 ) R8_CEILING, the rounded value.
!
  implicit none

  integer ( kind = 4 ) r8_ceiling
  real ( kind = 8 ) r
  integer ( kind = 4 ) value

  value = int ( r )
  if ( real ( value, kind = 8 ) < r ) then
    value = value + 1
  end if

  r8_ceiling = value

  return
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  implicit none
!
!  Coefficients for minimax approximation over (12, INF).
!
  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( kind = 8 ), parameter :: eps = 2.22D-16
  real ( kind = 8 ) fact
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) sum
  real ( kind = 8 ), parameter :: twelve = 12.0D+00
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 171.624D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.0D+30
  real ( kind = 8 ), parameter :: xminin = 2.23D-308
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  parity = .false.
  fact = one
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= zero ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= zero ) then

      if ( y1 /= aint ( y1 * half ) * two ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + one

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = one / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < twelve ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < one ) then

      z = y
      y = y + one
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - one

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = zero
    xden = one
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + one
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + one
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - half ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= one ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end
function r8_is_int ( r )

!*****************************************************************************80
!
!! R8_IS_INT determines if an R8 represents an integer value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the number to be checked.
!
!    Output, logical R8_IS_INT, is TRUE if R is an integer value.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  real ( kind = 8 ) r
  logical r8_is_int

  if ( real ( i4_huge ( ), kind = 8 ) < r ) then
    r8_is_int = .false.
  else if ( r < - real ( i4_huge ( ) , kind = 8 ) ) then
    r8_is_int = .false.
  else if ( r == real ( int ( r ), kind = 8 ) ) then
    r8_is_int = .true.
  else
    r8_is_int = .false.
  end if

  return
end
function r8_pi ( )

!*****************************************************************************80
!
!! R8_PI returns the value of pi to 16 decimal places.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision R8_PI, the value of pi.
!
  implicit none

  double precision r8_pi

  r8_pi = 3.141592653589793D+00

  return
end
function r8_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R8_UNIFORM returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
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
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 )seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  r8_uniform = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

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
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
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

  integer ( kind = 4 ) i4_huge
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
    seed = seed + i4_huge ( )
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
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
!    Input, character ( len = * ) TITLE, a title to be printed.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  logical                r8_is_int
  character ( len = * )  title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)') j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( r8_is_int ( a(i,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
function r8poly_value ( n, a, x )

!*****************************************************************************80
!
!! R8POLY_VALUE evaluates an R8POLY
!
!  Discussion:
!
!    For sanity's sake, the value of N indicates the NUMBER of
!    coefficients, or more precisely, the ORDER of the polynomial,
!    rather than the DEGREE of the polynomial.  The two quantities
!    differ by 1, but cause a great deal of confusion.
!
!    Given N and A, the form of the polynomial is:
!
!      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) R8POLY_VALUE, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8poly_value
  real ( kind = 8 ) x

  r8poly_value = a(n)
  do i = n - 1, 1, -1
    r8poly_value = r8poly_value * x + a(i)
  end do

  return
end
subroutine r8row_max ( m, n, x, ixmax, xmax )

!*****************************************************************************80
!
!! R8ROW_MAX returns the maximums of rows of an R8ROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the array.
!
!    Input, real ( kind = 8 ) X(M,N), the array to be examined.
!
!    Output, integer ( kind = 4 ) IXMAX(M); IXMAX(I) is the column of X in which
!    the maximum for row I occurs.
!
!    Output, real ( kind = 8 ) XMAX(M), the maximums of the rows of X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ixmax(m)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) xmax(m)

  do i = 1, m

    ixmax(i) = 1
    xmax(i) = x(i,1)
    do j = 2, n
      if ( xmax(i) < x(i,j) ) then
        ixmax(i) = j
        xmax(i) = x(i,j)
      end if
    end do

  end do

  return
end
subroutine r8row_mean ( m, n, x, mean )

!*****************************************************************************80
!
!! R8ROW_MEAN returns the means of rows of an R8ROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the array.
!
!    Input, real ( kind = 8 ) X(M,N), the array whose row means are desired.
!
!    Output, real ( kind = 8 ) MEAN(M), the means, or averages,
!    of the rows of X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) mean(m)
  real ( kind = 8 ) x(m,n)

  do i = 1, m
    mean(i) = sum ( x(i,1:n) ) / real ( n )
  end do

  return
end
subroutine r8row_min ( m, n, x, ixmin, xmin )

!*****************************************************************************80
!
!! R8ROW_MIN returns the minimums of rows of an R8ROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the array.
!
!    Input, real ( kind = 8 ) X(M,N), the array to be examined.
!
!    Output, integer ( kind = 4 ) IXMIN(M); IXMIN(I) is the column of X in which
!    the minimum for row I occurs.
!
!    Output, real ( kind = 8 ) XMIN(M), the minimums of the rows of X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ixmin(m)
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) xmin(m)

  do i = 1, m

    ixmin(i) = 1
    xmin(i) = x(i,1)
    do j = 2, n
      if ( x(i,j) < xmin(i) ) then
        ixmin(i) = j
        xmin(i) = x(i,j)
      end if
    end do

  end do

  return
end
subroutine r8row_variance ( m, n, x, variance )

!*****************************************************************************80
!
!! R8ROW_VARIANCE returns the variances of the rows of an R8ROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the array.
!
!    Input, real ( kind = 8 ) X(M,N), the array whose row means are desired.
!
!    Output, real ( kind = 8 ) VARIANCE(M), the variances of the rows of X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance(m)
  real ( kind = 8 ) x(m,n)

  do i = 1, m

    mean = sum ( x(i,1:n) ) / real ( n, kind = 8 )

    variance(i) = sum ( ( x(i,1:n) - mean )**2 )

    if ( 1 < n ) then
      variance(i) = variance(i) / real ( n - 1, kind = 8 )
    else
      variance(i) = 0.0D+00
    end if

  end do

  return
end
subroutine r8vec_circular_variance ( n, x, circular_variance )

!*****************************************************************************80
!
!! R8VEC_CIRCULAR_VARIANCE returns the circular variance of an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose variance is desired.
!
!    Output, real ( kind = 8 ) CIRCULAR VARIANCE, the circular variance
!    of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) circular_variance
  real ( kind = 8 ) mean
  real ( kind = 8 ) x(n)

  call r8vec_mean ( n, x, mean )

  circular_variance = &
      ( sum ( cos ( x(1:n) - mean ) ) )**2 &
    + ( sum ( sin ( x(1:n) - mean ) ) )**2

  circular_variance = sqrt ( circular_variance ) / real ( n, kind = 8 )

  circular_variance = 1.0D+00 - circular_variance

  return
end
function r8vec_diff_norm ( n, a, b )

!*****************************************************************************80
!
!! R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), B(N), the vectors
!
!    Output, real ( kind = 8 ) R8VEC_DIFF_NORM, the L2 norm of A - B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) r8vec_diff_norm

  r8vec_diff_norm = sqrt ( sum ( ( a(1:n) - b(1:n) )**2 ) )

  return
end
subroutine r8vec_max ( n, x, index, xmax )

!*****************************************************************************80
!
!! R8VEC_MAX returns the maximum value in an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) X(N), the array.
!
!    Output, integer ( kind = 4 ) INDEX, the index of the largest entry.
!
!    Output, real ( kind = 8 ) XMAX, the value of the largest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax

  if ( n <= 0 ) then

    index = 0
    xmax = 0.0D+00

  else

    index = 1
    xmax = x(1)

    do i = 2, n
      if ( xmax < x(i) ) then
        xmax = x(i)
        index = i
      end if
    end do

  end if

  return
end
subroutine r8vec_mean ( n, x, mean )

!*****************************************************************************80
!
!! R8VEC_MEAN returns the mean of an R8VEC.
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
!    Input, real ( kind = 8 ) X(N), the vector whose mean is desired.
!
!    Output, real ( kind = 8 ) MEAN, the mean, or average,
!    of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) mean
  real ( kind = 8 ) x(n)

  mean = sum ( x(1:n) ) / real ( n, kind = 8 )

  return
end
subroutine r8vec_min ( n, x, index, xmin )

!*****************************************************************************80
!
!! R8VEC_MIN returns the minimum value of an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) X(N), the array.
!
!    Output, integer ( kind = 4 ) INDEX, the index of the smallest entry.
!
!    Output, real ( kind = 8 ) XMIN, the value of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmin

  if ( n <= 0 ) then

    index = 0
    xmin = 0.0D+00

  else

    xmin = x(1)
    index = 1
    do i = 2, n
      if ( x(i) < xmin ) then
        xmin = x(i)
        index = i
      end if
    end do

  end if

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 December 1999
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
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) i, a(i)
  end do

  return
end
subroutine r8vec_uniform ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
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
!    Input, integer ( kind = 4 ) M, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge ( )
    end if

    r(i) = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

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
!    An R8VEC is a vector of real ( kind = 8 ) values.
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
  integer ( kind = 4 ) i4_huge
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
      seed = seed + i4_huge ( )
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec_unit_sum ( n, a )

!*****************************************************************************80
!
!! R8VEC_UNIT_SUM normalizes an R8VEC to have unit sum.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, real A(N), the vector to be normalized.  On output,
!    the entries of A should have unit sum.  However, if the input vector
!    has zero sum, the routine halts.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_sum

  a_sum = sum ( a(1:n) )

  if ( a_sum == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIT_SUM - Fatal error!'
    write ( *, '(a)' ) '  The vector entries sum to 0.'
    stop
  end if

  a(1:n) = a(1:n) / a_sum

  return
end
subroutine r8vec_variance ( n, x, variance )

!*****************************************************************************80
!
!! R8VEC_VARIANCE returns the variance of an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose variance is desired.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) mean
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(n)

  call r8vec_mean ( n, x, mean )

  variance = sum ( ( x(1:n) - mean )**2 )

  if ( 1 < n ) then
    variance = variance / real ( n - 1, kind = 8 )
  else
    variance = 0.0D+00
  end if

  return
end
subroutine random_initialize ( seed_input )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the FORTRAN90 random number generator
!    routine RANDOM_NUMBER, which is used by calls like
!
!      call random_number ( harvest = x )
!
!    then its behavior is not specified.  That may be OK for you.  But
!    you may want to be able to force the same sequence to be generated
!    each time, or to force a different sequence.
!
!    To control the sequence of random numbers, you need to set the seed.
!    In FORTRAN90, this is done by calling the RANDOM+SEED routine.
!    You can call it with no arguments, in fact.  But if you call
!    it with no arguments:
!
!      call random_seed ( )
!
!    then its behavior (or more particularly, the behavior of RANDOM_NUMBER)
!    is still not specified.  You might hope that the system will go out
!    and pick a nice random seed for you, but there's no guarantee.
!
!
!    For example, on the DEC ALPHA, if you compile a program that calls
!    RANDOM_NUMBER, then every time you run it, you get the same sequence
!    of "random" values.  If you compile a program that calls RANDOM_SEED
!    with no arguments, and then calls RANDOM_NUMBER, you still get the
!    same sequence each time you run the program.
!
!    In order to actually try to scramble up the random number generator
!    a bit, this routine goes through the tedious process of getting the
!    size of the random number seed, making up values based on the current
!    time, and setting the random number seed.
!
!    Unfortunately, the RANDOM_SEED routine has a very elastic definition.
!    It does not use a single scalar integer SEED.  Instead, it communicates
!    with the user through an integer vector whose size is not specified.
!    You actually have to "ask" the routine to tell you the size of this
!    vector.  Then you can fill up the vector with values to be used to
!    influence the seeding of the random number routine.  The details of
!    how the seed affects the sequence are also unspecified, but the only
!    thing we are reasonably confident about is that specifying the same
!    seed should result in the same sequence, and specifying different
!    seeds should result in different sequences!
!
!    I assume this is far more than you wanted to know.  (It's certainly
!    more than I wanted to know!)
!
!    The argument SEED is an input quantity only, so it is legal
!    to type
!
!      call random_initialize ( 0 )
!
!    or
!
!      call random_initialize ( 18867 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED_INPUT, the user's "suggestion" for a seed.
!    However, if the input value is 0, the routine will come up with
!    its own "suggestion", based on the system clock.
!
  implicit none

  integer ( kind = 4 ) count
  integer ( kind = 4 ) count_max
  integer ( kind = 4 ) count_rate
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_input
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
  real ( kind = 8 ) t
  integer ( kind = 4 ), parameter :: warm_up = 100

  seed = seed_input
!
!  Initialize the random seed routine.
!
  call random_seed ( )
!
!  Determine the size of the random number seed vector.
!
  call random_seed ( size = seed_size )
!
!  Allocate a vector of the right size to be used as a random seed.
!
  allocate ( seed_vector(seed_size) )
!
!  If the user supplied a SEED value, use that.
!
!  Otherwise, use the system clock value to make up a value that is
!  likely to change based on when this routine is called.
!
  if ( seed /= 0 ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, user SEED = ', seed
    end if

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, arbitrary SEED = ', &
        seed
    end if

  end if
!
!  Set the seed vector.  We don't know the significance of the
!  individual entries of the internal seed vector, so we'll just set
!  all entries to SEED.
!
  seed_vector(1:seed_size) = seed
!
!  Now call RANDOM_SEED, and tell it to use this seed vector.
!
  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times just to "warm it up".
!
  do i = 1, warm_up
    call random_number ( harvest = t )
  end do

  return
end
subroutine rayleigh_cdf ( x, a, cdf )

!*****************************************************************************80
!
!! RAYLEIGH_CDF evaluates the Rayleigh CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!    0.0D+00 <= X.
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    cdf = 0.0D+00
  else
    cdf = 1.0D+00 - exp ( - x**2 / ( 2.0D+00 * a**2 ) )
  end if

  return
end
subroutine rayleigh_cdf_inv ( cdf, a, x )

!*****************************************************************************80
!
!! RAYLEIGH_CDF_INV inverts the Rayleigh CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAYLEIGH_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = sqrt ( - 2.0D+00 * a * a * log ( 1.0D+00 - cdf ) )

  return
end
subroutine rayleigh_cdf_values ( n_data, sigma, x, fx )

!*****************************************************************************80
!
!! RAYLEIGH_CDF_VALUES returns some values of the Rayleigh CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = RayleighDistribution [ sigma ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) SIGMA, the shape parameter of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 9

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.8646647167633873D+00, &
    0.9996645373720975D+00, &
    0.9999999847700203D+00, &
    0.999999999999987D+00, &
    0.8646647167633873D+00, &
    0.3934693402873666D+00, &
    0.1992625970831920D+00, &
    0.1175030974154046D+00, &
    0.7688365361336422D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ), save, dimension ( n_max ) :: sigma_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    sigma = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    sigma = sigma_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function rayleigh_check ( a )

!*****************************************************************************80
!
!! RAYLEIGH_CHECK checks the parameter of the Rayleigh PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, logical RAYLEIGH_CHECK, is true if the parameter is legal.
!
  implicit none

  real ( kind = 8 ) a
  logical rayleigh_check

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAYLEIGH_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    rayleigh_check = .false.
    return
  end if

  rayleigh_check = .true.

  return
end
subroutine rayleigh_mean ( a, mean )

!*****************************************************************************80
!
!! RAYLEIGH_MEAN returns the mean of the Rayleigh PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) mean
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  mean = a * sqrt ( 0.5D+00 * pi )

  return
end
subroutine rayleigh_pdf ( x, a, pdf )

!*****************************************************************************80
!
!! RAYLEIGH_PDF evaluates the Rayleigh PDF.
!
!  Discussion:
!
!    PDF(A;X) = ( X / A^2 ) * EXP ( - X^2 / ( 2 * A^2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0 < A.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    pdf = 0.0D+00
  else
    pdf = ( x / a**2 ) * exp ( - x**2 / ( 2.0D+00 * a**2 ) )
  end if

  return
end
subroutine rayleigh_sample ( a, seed, x )

!*****************************************************************************80
!
!! RAYLEIGH_SAMPLE samples the Rayleigh PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call rayleigh_cdf_inv ( cdf, a, x )

  return
end
subroutine rayleigh_variance ( a, variance )

!*****************************************************************************80
!
!! RAYLEIGH_VARIANCE returns the variance of the Rayleigh PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameters of the PDF.
!    0.0D+00 < A.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) variance

  variance = 2.0D+00 * a**2 * ( 1.0D+00 - 0.25D+00 * pi )

  return
end
subroutine reciprocal_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! RECIPROCAL_CDF evaluates the Reciprocal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A <= B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then

    cdf = 0.0D+00

  else if ( 0.0D+00 < x ) then

    cdf = log ( a / x ) / log ( a / b )

  end if

  return
end
subroutine reciprocal_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! RECIPROCAL_CDF_INV inverts the Reciprocal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A <= B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RECIPROCAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = 0.0D+00
  else if ( 0.0D+00 < cdf ) then
    x = b**cdf / a**( cdf - 1.0D+00 )
  end if

  return
end
function reciprocal_check ( a, b )

!*****************************************************************************80
!
!! RECIPROCAL_CHECK checks the parameters of the Reciprocal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A <= B.
!
!    Output, logical RECIPROCAL_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical reciprocal_check

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RECIPROCAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.0'
    reciprocal_check = .false.
    return
  end if

  if ( b < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RECIPROCAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B < A'
    reciprocal_check = .false.
    return
  end if

  reciprocal_check = .true.

  return
end
subroutine reciprocal_mean ( a, b, mean )

!*****************************************************************************80
!
!! RECIPROCAL_MEAN returns the mean of the Reciprocal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A <= B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = ( a - b ) / log ( a / b )

  return
end
subroutine reciprocal_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! RECIPROCAL_PDF evaluates the Reciprocal PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = 1.0D+00 / ( X * LOG ( B / A ) )
!    for 0.0D+00 <= X
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A <= B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then
    pdf = 0.0D+00
  else if ( 0.0D+00 < x ) then
    pdf = 1.0D+00 / ( x * log ( b / a ) )
  end if

  return
end
subroutine reciprocal_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! RECIPROCAL_SAMPLE samples the Reciprocal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A <= B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  x = b**cdf / a**( cdf - 1.0D+00 )

  return
end
subroutine reciprocal_variance ( a, b, variance )

!*****************************************************************************80
!
!! RECIPROCAL_VARIANCE returns the variance of the Reciprocal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A <= B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d
  real ( kind = 8 ) variance

  d = log ( a / b )

  variance = ( a - b )* ( a * ( d - 2.0D+00 ) &
    + b * ( d + 2.0D+00 ) ) / ( 2.0D+00 * d**2 )

  return
end
subroutine ribesl ( x, alpha, nb, ize, b, ncalc )

!*****************************************************************************80
!
!! RIBESL calculates I Bessel function with non-integer orders.
!
!  Discussion:
!
!    This routine calculates Bessel functions I SUB(N+ALPHA) (X)
!    for non-negative argument X, and non-negative order N+ALPHA,
!    with or without exponential scaling.
!
!    This program is based on a program written by David
!    Sookne that computes values of the Bessel functions J or
!    I of real argument and integer order.  Modifications include
!    the restriction of the computation to the I Bessel function
!    of non-negative real argument, the extension of the computation
!    to arbitrary positive order, the inclusion of optional
!    exponential scaling, and the elimination of most underflow.
!
!    In case of an error, NCALC will not equal NB, and not all I's are
!    calculated to the desired accuracy.
!
!    If NCALC < 0:  An argument is out of range. For example,
!    NB <= 0, IZE is not 1 or 2, or IZE = 1 and EXPARG <= ABS(X)
!    In this case, the B-vector is not calculated, and NCALC is
!    set to MIN(NB,0)-1 so that NCALC /= NB.
!
!    If 0 < NCALC < NB, then not all requested function values could
!    be calculated accurately.  This usually occurs because NB is
!    much larger than ABS(X).  In this case, B(N) is calculated
!    to the desired accuracy for N <= NCALC, but precision
!    is lost for NCALC < N <= NB.  If B(N) does not vanish
!    for NCALC < N (because it is too small to be represented),
!    and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
!    significant figures of B(N) can be trusted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Frank Olver, David Sookne,
!    A Note on Backward Recurrence Algorithms,
!    Mathematics of Computation,
!    Volume 26, 1972, pages 941-947.
!
!    David Sookne,
!    Bessel Functions of Real Argument and Integer Order,
!    NBS Journal of Research B,
!    Volume 77B, 1973, pages 125-132.
!
!    William Cody,
!    Algorithm 597:
!    Sequence of Modified Bessel Functions of the First Kind,
!    ACM Transactions of Mathematical Software,
!    Volume 9, Number 2, June 1983, pages 242-245.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument for which the functions
!    are to be calculated.
!
!    Input, real ( kind = 8 ) ALPHA,the fractional part of the order
!    for which the functions are to be calculated.
!    0 <= ALPHA < 1.0.
!
!    Input, integer ( kind = 4 ) NB, the number of functions to be calculated.
!    The first function calculated is of order ALPHA, and the
!    last is of order (NB - 1 + ALPHA).  1 <= NB.
!
!    Input, integer ( kind = 4 ) IZE, scaling option.
!    1, unscaled I's are to calculated,
!    2, exponentially scaled I's are to be calculated.
!
!    Output, real ( kind = 8 ) B(NB), the values of the functions
!    I(ALPHA,X) through I(NB-1+ALPHA,X), with scaling if requested.
!
!    Output, integer ( kind = 4 ) NCALC, error indicator.
!    If NCALC = NB, then all the requested values were calculated
!    to the desired accuracy.
!
!  Local Parameeters:
!
!    BETA, the radix for the floating-point system.
!
!    MINEXP, smallest representable power of BETA.
!
!    MAXEXP, smallest power of BETA that overflows
!
!    IT, number of bits in the mantissa of a working precision variable.
!
!    NSIG, decimal significance desired.  Should be set to
!    INT(LOG10(2)*IT+1).  Setting NSIG lower will result
!    in decreased accuracy while setting NSIG higher will
!    increase CPU time without increasing accuracy.  The
!    truncation error is limited to a relative error of
!    T=.5*10^(-NSIG).
!
!    ENTEN, 10.0^K, where K is the largest integer such that
!    ENTEN is machine-representable in working precision
!
!    ENSIG, 10.0^NSIG
!
!    RTNSIG, 10.0^(-K) for the smallest integer K such that
!    NSIG/4 <= K.
!
!    ENMTEN, smallest ABS(X) such that X/4 does not underflow
!
!    XLARGE, upper limit on the magnitude of X when IZE=2.  Bear
!    in mind that if ABS(X)=N, then at least N iterations
!    of the backward recursion will be executed.  The value
!    of 10.0^4 is used on every machine.
!
!    EXPARG, largest working precision argument that the library
!    EXP routine can handle and upper limit on the
!    magnitude of X when IZE=1; approximately log(BETA^MAXEXP).
!
!    Approximate values for some important machines are:
!
!                        beta       minexp      maxexp       it
!
!  CRAY-1        (S.P.)    2        -8193        8191        48
!  Cyber 180/855
!    under NOS   (S.P.)    2         -975        1070        48
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)    2         -126         128        24
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)    2        -1022        1024        53
!  IBM 3033      (D.P.)   16          -65          63        14
!  VAX           (S.P.)    2         -128         127        24
!  VAX D-Format  (D.P.)    2         -128         127        56
!  VAX G-Format  (D.P.)    2        -1024        1023        53
!
!
!                        NSIG       ENTEN       ENSIG      RTNSIG
!
! CRAY-1        (S.P.)    15       1.0E+2465   1.0E+15     1.0E-4
! Cyber 180/855
!   under NOS   (S.P.)    15       1.0E+322    1.0E+15     1.0E-4
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)    16       1.0D+308    1.0D+16     1.0D-4
! IBM 3033      (D.P.)     5       1.0D+75     1.0D+5      1.0D-2
! VAX           (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
! VAX D-Format  (D.P.)    17       1.0D+38     1.0D+17     1.0D-5
! VAX G-Format  (D.P.)    16       1.0D+307    1.0D+16     1.0D-4
!
!
!                         ENMTEN      XLARGE   EXPARG
!
! CRAY-1        (S.P.)   1.84E-2466   1.0E+4    5677
! Cyber 180/855
!   under NOS   (S.P.)   1.25E-293    1.0E+4     741
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)   4.70E-38     1.0E+4      88
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)   8.90D-308    1.0D+4     709
! IBM 3033      (D.P.)   2.16D-78     1.0D+4     174
! VAX           (S.P.)   1.17E-38     1.0E+4      88
! VAX D-Format  (D.P.)   1.17D-38     1.0D+4      88
! VAX G-Format  (D.P.)   2.22D-308    1.0D+4     709
!
  implicit none

  integer ( kind = 4 ) nb

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(nb)
  real ( kind = 8 ), parameter :: const = 1.585D+00
  real ( kind = 8 ) em
  real ( kind = 8 ) empal
  real ( kind = 8 ) emp2al
  real ( kind = 8 ) en
  real ( kind = 8 ), parameter :: enmten = 8.9D-308
  real ( kind = 8 ), parameter :: ensig = 1.0D+16
  real ( kind = 8 ), parameter :: enten = 1.0D+308
  real ( kind = 8 ), parameter :: exparg = 709.0D+00
  logical flag
  real ( kind = 8 ), parameter :: half = 0.5D+00
  real ( kind = 8 ) halfx
  integer ( kind = 4 ) ize
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) magx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbmx
  integer ( kind = 4 ) ncalc
  integer ( kind = 4 ) nend
  integer ( kind = 4 ), parameter :: nsig = 16
  integer ( kind = 4 ) nstart
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) p
  real ( kind = 8 ) plast
  real ( kind = 8 ) pold
  real ( kind = 8 ) psave
  real ( kind = 8 ) psavel
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ), parameter :: rtnsig = 1.0D-04
  real ( kind = 8 ) tempa
  real ( kind = 8 ) tempb
  real ( kind = 8 ) tempc
  real ( kind = 8 ) test
  real ( kind = 8 ) total
  real ( kind = 8 ) tover
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xlarge = 1.0D+04
  real ( kind = 8 ), parameter :: zero = 0.0D+00
!
!  Check for X, NB, OR IZE out of range.
!
  if ( nb <= 0 ) then
    ncalc = min ( nb, 0 ) - 1
    return
  end if

  if ( x < 0.0D+00 ) then
    ncalc = min ( nb, 0 ) - 1
    return
  end if

  if ( alpha < 0.0D+00 ) then
    ncalc = min ( nb, 0 ) - 1
    return
  end if

  if ( 1.0D+00 <= alpha ) then
    ncalc = min ( nb, 0 ) - 1
    return
  end if

  if ( ize == 1 .and. exparg < x ) then
    ncalc = min ( nb, 0 ) - 1
    return
  end if

  if ( ize == 2 .and. xlarge < x ) then
    ncalc = min ( nb, 0 ) - 1
    return
  end if
!
!  Use 2-term ascending series for small X.
!
  ncalc = nb
  magx = int ( x )
!
!  Initialize the forward sweep, the P-sequence of Olver.
!
  if ( rtnsig <= x ) then

    nbmx = nb - magx
    n = magx + 1
    en = real ( n + n, kind = 8 ) + ( alpha + alpha )
    plast = one
    p = en / x
!
!  Calculate general significance test.
!
    test = ensig + ensig

    if ( 5 * nsig < 2 * magx ) then
      test = sqrt ( test * p )
    else
      test = test / const**magx
    end if
!
!  Calculate P-sequence until N = NB-1.  Check for possible overflow.
!
    flag = .false.

    if ( 3 <= nbmx ) then

      tover = enten / ensig
      nstart = magx + 2
      nend = nb - 1

      do k = nstart, nend

        n = k
        en = en + two
        pold = plast
        plast = p
        p = en * plast / x + pold
!
!  To avoid overflow, divide P-sequence by TOVER.  Calculate
!  P-sequence until 1 < ABS(P).
!
        if ( tover < p ) then

          tover = enten
          p = p / tover
          plast = plast / tover
          psave = p
          psavel = plast
          nstart = n + 1

          do

            n = n + 1
            en = en + two
            pold = plast
            plast = p
            p = en * plast / x + pold

            if ( 1.0D+00 < p ) then
              exit
            end if

          end do

          tempb = en / x
!
!  Calculate backward test, and find NCALC, the highest N
!  such that the test is passed.
!
          test = pold * plast / ensig
          test = test * ( half - half / ( tempb * tempb ) )
          p = plast * tover
          n = n - 1
          en = en - two
          nend = min ( nb, n )

          ncalc = nend + 1

          do l = nstart, nend

            pold = psavel
            psavel = psave
            psave = en * psavel / x + pold

            if ( test < psave * psavel ) then
              ncalc = l
              exit
            end if

          end do

          ncalc = ncalc - 1
          flag = .true.
          exit

        end if

      end do

      if ( .not. flag ) then

        n = nend
        en = real ( n + n, kind = 8 ) + ( alpha + alpha )
!
!  Calculate special significance test for 2 < NBMX.
!
        test = max ( test, sqrt ( plast * ensig ) * sqrt ( p + p ) )

      end if

    end if
!
!  Calculate P-sequence until significance test passed.
!
    if ( .not. flag ) then

      do

        n = n + 1
        en = en + two
        pold = plast
        plast = p
        p = en * plast / x + pold

        if ( test <= p ) then
          exit
        end if

      end do

    end if
!
!  Initialize the backward recursion and the normalization sum.
!
    n = n + 1
    en = en + two
    tempb = zero
    tempa = one / p
    em = real ( n, kind = 8 ) - one
    empal = em + alpha
    emp2al = ( em - one ) + ( alpha + alpha )
    total = tempa * empal * emp2al / em
    nend = n - nb
!
!  N < NB, so store B(N) and set higher orders to zero.
!
    if ( nend < 0 ) then

      b(n) = tempa
      nend = -nend

      do l = 1, nend
        b(n+l) = zero
      end do

      nend = n - 2
!
!  Calculate via difference equation and store B(N), until N = 2.
!
      if ( 0 < nend ) then

        do l = 1, nend
          n = n - 1
          en = en - two
          b(n) = ( en * b(n+1) ) / x + b(n+2)
          em = em - one
          emp2al = emp2al - one
          if ( n == 2 ) then
            emp2al = one
          end if
          empal = empal - one
          total = ( total + b(n) * empal ) * emp2al / em
        end do

      end if
!
!  Calculate B(1).
!
      b(1) = two * empal * b(2) / x + b(3)

      total = ( total + total ) + b(1)
!
!  Recur backward via difference equation, calculating (but
!  not storing) B(N), until N = NB.
!
    else

      if ( 0 < nend ) then

        do l = 1, nend

          n = n - 1
          en = en - two
          tempc = tempb
          tempb = tempa
          tempa = ( en * tempb ) / x + tempc
          em = em - one
          emp2al = emp2al - one

          if ( n == 1 ) then
            exit
          end if

          if ( n == 2 ) then
            emp2al = one
          end if

          empal = empal - one
          total = ( total + tempa * empal ) * emp2al / em

        end do

      end if
!
!  Store B(NB).
!
      b(n) = tempa

      if ( nb <= 1 ) then

        total = ( total + total ) + tempa
!
!  Calculate and Store B(NB-1).
!
      else

        n = n - 1
        en = en - two
        b(n) = ( en * tempa ) / x + tempb

        if ( 1 < n  ) then

          em = em - one
          emp2al = emp2al - one

          if ( n == 2 ) then
            emp2al = one
          end if

          empal = empal - one
          total = ( total + b(n) * empal ) * emp2al / em

          nend = n - 2
!
!  Calculate via difference equation and store B(N), until N = 2.
!
          if ( 0 < nend ) then

            do l = 1, nend
              n = n - 1
              en = en - two
              b(n) = ( en * b(n+1) ) / x + b(n+2)
              em = em - one
              emp2al = emp2al - one
              if ( n == 2 ) then
                emp2al = one
              end if
              empal = empal - one
              total = ( total + b(n) * empal ) * emp2al / em
            end do

          end if
!
!  Calculate B(1).
!
          b(1) = two * empal * b(2) / x + b(3)

        end if

        total = ( total + total ) + b(1)

      end if

    end if
!
!  Normalize.  Divide all B(N) by TOTAL.
!


    if ( alpha /= zero ) then
       total = total * r8_gamma ( one + alpha ) * ( x * half )**( -alpha )
    end if

    if ( ize == 1 ) then
      total = total * exp ( -x )
    end if

    tempa = enmten

    if ( 1.0D+00 < total ) then
      tempa = tempa * total
    end if

    do n = 1, nb
      if ( b(n) < tempa ) then
        b(n) = zero
      end if
      b(n) = b(n) / total
    end do

    return
!
!  Two-term ascending series for small X.
!
  else

    tempa = one
    empal = one + alpha
    halfx = zero

    if ( enmten < x ) then
      halfx = half * x
    end if

    if ( alpha /= zero ) then
      tempa = halfx**alpha / r8_gamma ( empal )
    end if

    if ( ize == 2 ) then
      tempa = tempa * exp ( - x )
    end if

    tempb = zero

    if ( one < x + one ) then
      tempb = halfx * halfx
    end if

    b(1) = tempa + tempa * tempb / empal

    if ( x /= zero .and. b(1) == zero ) then
      ncalc = 0
    end if

    if ( 1 < nb ) then

      if ( x == zero ) then

        b(2:nb) = zero
!
!  Calculate higher-order functions.
!
      else

        tempc = halfx
        tover = ( enmten + enmten ) / x

        if ( tempb /= zero ) then
          tover = enmten / tempb
        end if

        do n = 2, nb

          tempa = tempa / empal
          empal = empal + one
          tempa = tempa * tempc

          if ( tempa <= tover * empal ) then
            tempa = zero
          end if

          b(n) = tempa + tempa * tempb / empal

          if ( b(n) == zero .and. n < ncalc ) then
            ncalc = n - 1
          end if

        end do

      end if

    end if

  end if

  return
end
subroutine runs_mean ( m, n, mean )

!*****************************************************************************80
!
!! RUNS_MEAN returns the mean of the Runs PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the PDF.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  integer ( kind = 4 ) m
  real ( kind = 8 ) mean
  integer ( kind = 4 ) n

  mean = real ( m + 2 * m * n + n, kind = 8 ) &
       / real ( m             + n, kind = 8 )

  return
end
subroutine runs_pdf ( m, n, r, pdf )

!*****************************************************************************80
!
!! RUNS_PDF evaluates the Runs PDF.
!
!  Discussion:
!
!    Suppose we have M symbols of one type and N of another, and we consider
!    the various possible permutations of these symbols.
!
!    Let "R" be the number of runs in a given permutation.  By a "run", we
!    mean a maximal sequence of identical symbols.  Thus, for instance,
!    the permutation
!
!      ABBBAAAAAAAA
!
!    has three runs.
!
!    The probability that a permutation of M+N symbols, with M of one kind
!    and N of another, will have exactly R runs is:
!
!      PDF(M,N)(R) = 2 * C(M-1,R/2-1) * C(N-1,R/2-1)
!                    / C(M+N,N) for R even;
!
!                  = ( C(M-1,(R-1)/2) * C(N-1,(R-3)/2 )
!                    + C(M-1,(R-3)/2) * C(N-1,(R-1)/2 )
!                    ) / C(M+N,N) for R odd.
!
!    The minimum number of runs is:
!
!      1 if M or N is 0,
!      2 otherwise.
!
!    The maximum number of runs is:
!
!      M + N,                if M = N
!      2 * min ( M, N ) + 1  otherwise
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kalimutha Krishnamoorthy,
!    Handbook of Statistical Distributions with Applications,
!    Chapman and Hall, 2006,
!    ISBN: 1-58488-635-8,
!    LC: QA273.6.K75.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the PDF.
!
!    Input, integer ( kind = 4 ) R, the number of runs.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) combinatorial
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) r

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RUN_PDF - Fatal error!'
    write ( *, '(a)' ) '  M must be at least 0.'
    write ( *, '(a,i8)' ) '  The input value of M = ', m
    stop
  end if

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RUN_PDF - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 0.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop
  end if

  if ( n + m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RUN_PDF - Fatal error!'
    write ( *, '(a)' ) '  M+N must be at least 1.'
    write ( *, '(a,i8)' ) '  The input value of M+N = ', m + n
    stop
  end if
!
!  If all the symbols are of one type, there is always 1 run.
!
  if ( m == 0 .or. n == 0 ) then
    if ( r == 1 ) then
      pdf = 1.0D+00
    else
      pdf = 0.0D+00
    end if
    return
  end if
!
!  Take care of extreme values of R.
!
  if ( r < 2 .or. m + n < r ) then
    pdf = 0.0D+00
    return
  end if
!
!  The normal cases.
!
  if ( mod ( r, 2 ) == 0 ) then

    pdf = real ( 2 * combinatorial ( m - 1, ( r / 2 ) - 1 ) &
                   * combinatorial ( n - 1, ( r / 2 ) - 1 ), kind = 8 ) &
        / real (     combinatorial ( m + n, n ), kind = 8 )

  else

    pdf = real (   combinatorial ( m - 1, ( r - 1 ) / 2 ) &
                 * combinatorial ( n - 1, ( r - 3 ) / 2 ) &
                 + combinatorial ( m - 1, ( r - 3 ) / 2 ) &
                 * combinatorial ( n - 1, ( r - 1 ) / 2 ), kind = 8 ) &
        / real (   combinatorial ( m + n, n ), kind = 8 )

  end if

  return
end
subroutine runs_sample ( m, n, seed, r )

!*****************************************************************************80
!
!! RUNS_SAMPLE samples the Runs PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) R, the number of runs.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m+n)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) seed

  call runs_simulate ( m, n, seed, a )

  call i4vec_run_count ( m+n, a, r )

  return
end
subroutine runs_simulate ( m, n, seed, a )

!*****************************************************************************80
!
!! RUNS_SIMULATE simulates a case governed by the Runs PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) A(M+N), a sequence of M 0's and N 1's chosen
!    uniformly at random.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m+n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  a(1:m) = 0
  a(m+1:m+n) = 1

  do i = 1, m+n-1

    j = i4_uniform ( i, m+n, seed )

    k    = a(i)
    a(i) = a(j)
    a(j) = k

  end do

  return
end
subroutine runs_variance ( m, n, variance )

!*****************************************************************************80
!
!! RUNS_VARIANCE returns the variance of the Runs PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the PDF.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) variance

  variance = real ( 2 * m * n * ( 2 * m * n - m - n ), kind = 8 ) &
           / real ( ( m + n ) * ( m + n ) * ( m + n - 1 ), kind = 8 )

  return
end
function sech ( x )

!*****************************************************************************80
!
!! SECH returns the hyperbolic secant.
!
!  Discussion:
!
!    SECH ( X ) = 1.0D+00 / COSH ( X ) = 2.0D+00 / ( EXP ( X ) + EXP ( - X ) )
!
!    SECH is not a built-in function in FORTRAN, and occasionally it
!    is handier, or more concise, to be able to refer to it directly
!    rather than through its definition in terms of the sine function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) SECH, the hyperbolic secant of X.
!
  implicit none

  real ( kind = 8 ) sech
  real ( kind = 8 ) x

  sech = 1.0D+00 / cosh ( x )

  return
end
subroutine sech_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! SECH_CDF evaluates the Hyperbolic Secant CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameter of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  y = ( x - a ) / b

  cdf = 2.0D+00 * atan ( exp ( y ) ) / pi

  return
end
subroutine sech_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! SECH_CDF_INV inverts the Hyperbolic Secant CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SECH_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = - huge ( x )
  else if ( cdf < 1.0D+00 ) then
    x = a + b * log ( tan ( 0.5D+00 * pi * cdf ) )
  else if ( 1.0D+00 == cdf ) then
    x = huge ( x )
  end if

  return
end
function sech_check ( a, b )

!*****************************************************************************80
!
!! SECH_CHECK checks the parameters of the Hyperbolic Secant CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameter of the PDF.
!    0.0D+00 < B.
!
!    Output, logical SECH_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical sech_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SECH_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.0'
    sech_check = .false.
    return
  end if

  sech_check = .true.

  return
end
subroutine sech_mean ( a, b, mean )

!*****************************************************************************80
!
!! SECH_MEAN returns the mean of the Hyperbolic Secant PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine sech_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! SECH_PDF evaluates the Hypebolic Secant PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = sech ( ( X - A ) / B ) / ( PI * B )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) sech
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  y = ( x - a ) / b

  pdf = sech ( y ) / ( pi * b )

  return
end
subroutine sech_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! SECH_SAMPLE samples the Hyperbolic Secant PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  x = a + b * log ( tan ( 0.5D+00 * pi * cdf ) )

  return
end
subroutine sech_variance ( a, b, variance )

!*****************************************************************************80
!
!! SECH_VARIANCE returns the variance of the Hyperbolic Secant PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) variance

  variance = 0.25D+00 * ( pi * b )**2

  return
end
subroutine semicircular_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! SEMICIRCULAR_CDF evaluates the Semicircular CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameter of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x <= a - b ) then

    cdf = 0.0D+00

  else if ( x <= a + b ) then

    y = ( x - a ) / b

    cdf = 0.5D+00 + ( y * sqrt ( 1.0D+00 - y**2 ) + asin ( y ) ) / pi

  else if ( a + b < x ) then

    cdf = 1.0D+00

  end if

  return
end
subroutine semicircular_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! SEMICIRCULAR_CDF_INV inverts the Semicircular CDF.
!
!  Discussion:
!
!    A simple bisection method is used on the interval [ A - B, A + B ].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf1
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) cdf3
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  real ( kind = 8 ), parameter :: tol = 0.0001D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SEMICIRCULAR_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = a - b
    return
  else if ( 1.0D+00 == cdf ) then
    x = a + b
    return
  end if

  x1 = a - b
  cdf1 = 0.0D+00

  x2 = a + b
  cdf2 = 1.0D+00
!
!  Now use bisection.
!
  it = 0

  do

    it = it + 1

    x3 = 0.5D+00 * ( x1 + x2 )
    call semicircular_cdf ( x3, a, b, cdf3 )

    if ( abs ( cdf3 - cdf ) < tol ) then
      x = x3
      exit
    end if

    if ( it_max < it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SEMICIRCULAR_CDF_INV - Fatal error!'
      write ( *, '(a)' ) '  Iteration limit exceeded.'
      stop
    end if

    if ( sign ( 1.0D+00, cdf3 - cdf ) == sign ( 1.0D+00, cdf1 - cdf ) ) then
      x1 = x3
      cdf1 = cdf3
    else
      x2 = x3
      cdf2 = cdf3
    end if

  end do

  return
end
function semicircular_check ( a, b )

!*****************************************************************************80
!
!! SEMICIRCULAR_CHECK checks the parameters of the Semicircular CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameter of the PDF.
!    0.0D+00 < B.
!
!    Output, logical SEMICIRCULAR_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical semicircular_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SEMICIRCULAR_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.0'
    semicircular_check = .false.
    return
  end if

  semicircular_check = .true.

  return
end
subroutine semicircular_mean ( a, b, mean )

!*****************************************************************************80
!
!! SEMICIRCULAR_MEAN returns the mean of the Semicircular PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine semicircular_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! SEMICIRCULAR_PDF evaluates the Semicircular PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = ( 2 / ( B * PI ) ) * SQRT ( 1 - ( ( X - A ) / B )^2 )
!    for A - B <= X <= A + B
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x < a - b ) then

    pdf = 0.0D+00

  else if ( x <= a + b ) then

    y = ( x - a ) / b

    pdf = 2.0D+00 / ( b * pi ) * sqrt ( 1.0D+00 - y**2 )

  else if ( a + b < x ) then

    pdf = 0.0D+00

  end if

  return
end
subroutine semicircular_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! SEMICIRCULAR_SAMPLE samples the Semicircular PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) angle
  real ( kind = 8 ) b
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) radius
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  radius = r8_uniform_01 ( seed )
  radius = b * sqrt ( radius )
  angle = pi * r8_uniform_01 ( seed )
  x = a + radius * cos ( angle )

  return
end
subroutine semicircular_variance ( a, b, variance )

!*****************************************************************************80
!
!! SEMICIRCULAR_VARIANCE returns the variance of the Semicircular PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = b * b / 4.0D+00

  return
end
function sin_power_int ( a, b, n )

!*****************************************************************************80
!
!! SIN_POWER_INT evaluates the sine power integral.
!
!  Discussion:
!
!    The function is defined by
!
!      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin ( t ))^n dt
!
!    The algorithm uses the following fact:
!
!      Integral sin^n ( t ) = (1/n) * (
!        sin^(n-1)(t) * cos(t) + ( n-1 ) * Integral sin^(n-2) ( t ) dt )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, integer ( kind = 4 ) N, the power of the sine function.
!
!    Output, real ( kind = 8 ) SIN_POWER_INT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) ca
  real ( kind = 8 ) cb
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mlo
  integer ( kind = 4 ) n
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) sin_power_int
  real ( kind = 8 ) value

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIN_POWER_INT - Fatal error!'
    write ( *, '(a)' ) '  Power N < 0.'
    value = 0.0D+00
    stop
  end if

  sa = sin ( a )
  sb = sin ( b )
  ca = cos ( a )
  cb = cos ( b )

  if ( mod ( n, 2 ) == 0 ) then

    value = b - a
    mlo = 2
  else
    value = ca - cb
    mlo = 3
  end if

  do m = mlo, n, 2
    value = ( real ( m - 1, kind = 8 ) * value &
              + sa**(m-1) * ca - sb**(m-1) * cb ) &
      / real ( m, kind = 8 )
  end do

  sin_power_int = value

  return
end
function sphere_unit_area_nd ( dim_num )

!*****************************************************************************80
!
!! SPHERE_UNIT_AREA_ND computes the surface area of a unit sphere in ND.
!
!  Discussion:
!
!    The unit sphere in ND satisfies:
!
!      sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
!
!    Results for the first few values of N are:
!
!    DIM_NUM   Area
!
!     2    2        * PI
!     3    4        * PI
!     4  ( 2 /   1) * PI^2
!     5  ( 8 /   3) * PI^2
!     6  ( 1 /   1) * PI^3
!     7  (16 /  15) * PI^3
!     8  ( 1 /   3) * PI^4
!     9  (32 / 105) * PI^4
!    10  ( 1 /  12) * PI^5
!
!    For the unit sphere, Area(DIM_NUM) = DIM_NUM * Volume(DIM_NUM)
!
!    Sphere_Unit_Area ( DIM_NUM ) = 2 * PI**(DIM_NUM/2) / Gamma ( DIM_NUM / 2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Output, real ( kind = 8 ) SPHERE_UNIT_AREA_ND, the area of the sphere.
!
  implicit none

  real ( kind = 8 ) area
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) sphere_unit_area_nd

  if ( mod ( dim_num, 2 ) == 0 ) then
    m = dim_num / 2
    area = 2.0D+00 * ( pi )**m
    do i = 1, m - 1
      area = area / real ( i, kind = 8 )
    end do
  else
    m = ( dim_num - 1 ) / 2
    area = ( pi )**m * 2.0D+00**dim_num
    do i = m + 1, 2 * m
      area = area / real ( i,  kind = 8 )
    end do
  end if

  sphere_unit_area_nd = area

  return
end
function stirling2_value ( n, m )

!*****************************************************************************80
!
!! STIRLING2_VALUE computes a Stirling number of the second kind.
!
!  Discussion:
!
!    S2(N,M) represents the number of distinct partitions of N elements
!    into M nonempty sets.  For a fixed N, the sum of the Stirling
!    numbers S2(N,M) is represented by B(N), called "Bell's number",
!    and represents the number of distinct partitions of N elements.
!
!    For example, with 4 objects, there are:
!
!    1 partition into 1 set:
!
!      (A,B,C,D)
!
!    7 partitions into 2 sets:
!
!      (A,B,C) (D)
!      (A,B,D) (C)
!      (A,C,D) (B)
!      (A) (B,C,D)
!      (A,B) (C,D)
!      (A,C) (B,D)
!      (A,D) (B,C)
!
!    6 partitions into 3 sets:
!
!      (A,B) (C) (D)
!      (A) (B,C) (D)
!      (A) (B) (C,D)
!      (A,C) (B) (D)
!      (A,D) (B) (C)
!      (A) (B,D) (C)
!
!    1 partition into 4 sets:
!
!      (A) (B) (C) (D)
!
!    So S2(4,1) = 1, S2(4,2) = 7, S2(4,3) = 6, S2(4,4) = 1, and B(4) = 15.
!
!
!  First terms:
!
!    N/M: 1    2    3    4    5    6    7    8
!
!    1    1    0    0    0    0    0    0    0
!    2    1    1    0    0    0    0    0    0
!    3    1    3    1    0    0    0    0    0
!    4    1    7    6    1    0    0    0    0
!    5    1   15   25   10    1    0    0    0
!    6    1   31   90   65   15    1    0    0
!    7    1   63  301  350  140   21    1    0
!    8    1  127  966 1701 1050  266   28    1
!
!  Recursion:
!
!    S2(N,1) = 1 for all N.
!    S2(I,I) = 1 for all I.
!    S2(I,J) = 0 if I < J.
!
!    S2(N,M) = M * S2(N-1,M) + S2(N-1,M-1)
!
!  Properties:
!
!    sum ( 1 <= K <= M ) S2(I,K) * S1(K,J) = Delta(I,J)
!
!    X**N = sum ( 0 <= K <= N ) S2(N,K) X_K
!    where X_K is the falling factorial function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows of the table.
!
!    Input, integer ( kind = 4 ) M, the number of columns of the table.
!
!    Output, integer ( kind = 4 ) STIRLING2_VALUE, the value of S2(N,M).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s2(n,m)
  integer ( kind = 4 ) stirling2_value

  if ( n <= 0 ) then
    stirling2_value = 0
    return
  end if

  if ( m <= 0 ) then
    stirling2_value = 0
    return
  end if

  s2(1,1) = 1
  s2(1,2:m) = 0

  do i = 2, n

    s2(i,1) = 1

    do j = 2, m
      s2(i,j) = j * s2(i-1,j) + s2(i-1,j-1)
    end do

  end do

  stirling2_value = s2(n,m)

  return
end
subroutine student_cdf ( x, a, b, c, cdf )

!*****************************************************************************80
!
!! STUDENT_CDF evaluates the central Student T CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, shape parameters of the PDF,
!    used to transform the argument X to a shifted and scaled
!    value Y = ( X - A ) / B.  It is required that B be nonzero.
!    For the standard distribution, A = 0 and B = 1.
!
!    Input, real ( kind = 8 ) C, is usually called the number of
!    degrees of freedom of the distribution.  C is typically an
!    integer, but that is not essential.  It is required that
!    C be strictly positive.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) beta_inc
  real ( kind = 8 ) c
  real ( kind = 8 ) c2
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  y = ( x - a ) / b

  a2 = 0.5D+00 * c
  b2 = 0.5D+00
  c2 = c / ( c + y * y )

  if ( y <= 0.0D+00 ) then
    cdf = 0.5D+00 * beta_inc ( a2, b2, c2 )
  else
    cdf = 1.0D+00 - 0.5D+00 * beta_inc ( a2, b2, c2 )
  end if

  return
end
subroutine student_cdf_values ( n_data, c, x, fx )

!*****************************************************************************80
!
!! STUDENT_CDF_VALUES returns some values of the Student CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = StudentTDistribution [ c ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) C, is usually called the number of
!    degrees of freedom of the distribution.  C is typically an
!    integer, but that is not essential.  It is required that
!    C be strictly positive.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) c
  real ( kind = 8 ), save, dimension ( n_max ) :: c_vec = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, &
    5.0D+00, 2.0D+00, 5.0D+00, 2.0D+00, &
    5.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, &
    5.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6000231200328521D+00, &
    0.6001080279134390D+00, &
    0.6001150934648930D+00, &
    0.6000995134721354D+00, &
    0.5999341989834830D+00, &
    0.7498859393137811D+00, &
    0.7500879487671045D+00, &
    0.9500004222186464D+00, &
    0.9499969138365968D+00, &
    0.9900012348724744D+00, &
    0.9900017619355059D+00, &
    0.9900004567580596D+00, &
    0.9900007637471291D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.325D+00, &
    0.289D+00, &
    0.277D+00, &
    0.271D+00, &
    0.267D+00, &
    0.816D+00, &
    0.727D+00, &
    2.920D+00, &
    2.015D+00, &
    6.965D+00, &
    4.541D+00, &
    3.747D+00, &
    3.365D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    c = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    c = c_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function student_check ( a, b, c )

!*****************************************************************************80
!
!! STUDENT_CHECK checks the parameter of the central Student T CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, shape parameters of the PDF,
!    used to transform the argument X to a shifted and scaled
!    value Y = ( X - A ) / B.  It is required that B be nonzero.
!    For the standard distribution, A = 0 and B = 1.
!
!    Input, real ( kind = 8 ) C, is usually called the number of
!    degrees of freedom of the distribution.  C is typically an
!    integer, but that is not essential.  It is required that
!    C be strictly positive.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical student_check

  if ( b == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STUDENT_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B must be nonzero.'
    student_check = .false.
    return
  end if

  if ( c <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STUDENT_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C must be greater than 0.'
    student_check = .false.
    return
  end if

  student_check = .true.

  return
end
subroutine student_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! STUDENT_MEAN returns the mean of the central Student T PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, shape parameters of the PDF,
!    used to transform the argument X to a shifted and scaled
!    value Y = ( X - A ) / B.  It is required that B be nonzero.
!    For the standard distribution, A = 0 and B = 1.
!
!    Input, real ( kind = 8 ) C, is usually called the number of
!    degrees of freedom of the distribution.  C is typically an
!    integer, but that is not essential.  It is required that
!    C be strictly positive.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine student_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! STUDENT_PDF evaluates the central Student T PDF.
!
!  Discussion:
!
!    PDF(A,B,C;X) = Gamma ( (C+1)/2 ) /
!      ( Gamma ( C / 2 ) * Sqrt ( PI * C )
!      * ( 1 + ((X-A)/B)^2/C )^(C + 1/2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, shape parameters of the PDF,
!    used to transform the argument X to a shifted and scaled
!    value Y = ( X - A ) / B.  It is required that B be nonzero.
!    For the standard distribution, A = 0 and B = 1.
!
!    Input, real ( kind = 8 ) C, is usually called the number of
!    degrees of freedom of the distribution.  C is typically an
!    integer, but that is not essential.  It is required that
!    C be strictly positive.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  y = ( x - a ) / b

  pdf = r8_gamma ( 0.5D+00 * ( c + 1.0D+00 ) ) / ( sqrt ( pi * c ) &
    * r8_gamma ( 0.5D+00 * c ) &
    * sqrt ( ( 1.0D+00 + y * y / c )**( 2 * c + 1.0D+00 ) ) )

  return
end
subroutine student_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! STUDENT_SAMPLE samples the central Student T PDF.
!
!  Discussion:
!
!    For the sampling algorithm, it is necessary that 2 < C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, shape parameters of the PDF,
!    used to transform the argument X to a shifted and scaled
!    value Y = ( X - A ) / B.  It is required that B be nonzero.
!    For the standard distribution, A = 0 and B = 1.
!
!    Input, real ( kind = 8 ) C, is usually called the number of
!    degrees of freedom of the distribution.  C is typically an
!    integer, but that is not essential.  It is required that
!    C be strictly positive.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) c
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  if ( c < 3.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STUDENT_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Sampling fails for C <= 2.'
    return
  end if

  a2 = 0.0D+00
  b2 = c / ( c - 2 )

  call normal_sample ( a2, b2, seed, x2 )

  call chi_square_sample ( c, seed, x3 )
  x3 = x3 * c / ( c - 2.0D+00 )

  x = a + b * x2 * sqrt ( c ) / x3

  return
end
subroutine student_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! STUDENT_VARIANCE returns the variance of the central Student T PDF.
!
!  Discussion:
!
!    The variance is not defined unless 2 < C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, shape parameters of the PDF,
!    used to transform the argument X to a shifted and scaled
!    value Y = ( X - A ) / B.  It is required that B be nonzero.
!    For the standard distribution, A = 0 and B = 1.
!
!    Input, real ( kind = 8 ) C, is usually called the number of
!    degrees of freedom of the distribution.  C is typically an
!    integer, but that is not essential.  It is required that
!    C be strictly positive.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) variance

  if ( c <= 2.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STUDENT_VARIANCE - Fatal error!'
    write ( *, '(a)' ) '  Variance not defined for C <= 2.'
    stop
  end if

  variance = b * b * c / ( c - 2.0D+00 )

  return
end
subroutine student_noncentral_cdf ( x, idf, d, cdf )

!*****************************************************************************80
!
!! STUDENT_NONCENTRAL_CDF evaluates the noncentral Student T CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by B E Cooper.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    BE Cooper,
!    Algorithm AS 5:
!    The Integral of the Non-Central T-Distribution,
!    Applied Statistics,
!    Volume 17, 1968, page 193.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, integer ( kind = 4 ) IDF, the number of degrees of freedom.
!
!    Input, real ( kind = 8 ) D, the noncentrality parameter.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  integer ( kind = 4 ), parameter :: a_max = 100
  real ( kind = 8 ) ak
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) d
  real ( kind = 8 ) drb
  real ( kind = 8 ), parameter :: emin = 12.5D+00
  real ( kind = 8 ) f
  real ( kind = 8 ) fk
  real ( kind = 8 ) fmkm1
  real ( kind = 8 ) fmkm2
  real ( kind = 8 ) gamma_log
  integer ( kind = 4 ) idf
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp
  real ( kind = 8 ) tfn
  real ( kind = 8 ) x

  f = real ( idf, kind = 8 )

  if ( idf == 1 ) then

    a = x / sqrt ( f )
    b = f / ( f + x**2 )
    drb = d * sqrt ( b )

    call normal_01_cdf ( drb, cdf2 )
    cdf = 1.0D+00 - cdf2 + 2.0D+00 * tfn ( drb, a )

  else if ( idf <= a_max ) then

    a = x / sqrt ( f )
    b = f / ( f + x * x )
    drb = d * sqrt ( b )
    sum2 = 0.0D+00

    fmkm2 = 0.0D+00
    if ( abs ( drb ) < emin ) then
      call normal_01_cdf ( a * drb, cdf2 )
      fmkm2 = a * sqrt ( b ) * exp ( - 0.5D+00 * drb**2 ) * cdf2 &
        / sqrt ( 2.0D+00 * pi )
    end if

    fmkm1 = b * d * a * fmkm2
    if ( abs ( d ) < emin ) then
      fmkm1 = fmkm1 + 0.5D+00 * b * a * exp ( - 0.5D+00 * d**2 ) / pi
    end if

    if ( mod ( idf, 2 ) == 0 ) then
      sum2 = fmkm2
    else
      sum2 = fmkm1
    end if

    ak = 1.0D+00

    do k = 2, idf - 2, 2

      fk = real ( k, kind = 8 )

      fmkm2 = b * ( d * a * ak * fmkm1 + fmkm2 ) * ( fk - 1.0D+00 ) / fk

      ak = 1.0D+00 / ( ak * ( fk - 1.0D+00 ) )
      fmkm1 = b * ( d * a * ak * fmkm2 + fmkm1 ) * fk / ( fk + 1.0D+00 )

      if ( mod ( idf, 2 ) == 0 ) then
        sum2 = sum2 + fmkm2
      else
        sum2 = sum2 + fmkm1
      end if

      ak = 1.0D+00 / ( ak * fk )

    end do

    if ( mod ( idf, 2 ) == 0 ) then
      call normal_01_cdf ( d, cdf2 )
      cdf = 1.0D+00 - cdf2 + sum2 * sqrt ( 2.0D+00 * pi )
    else
      call normal_01_cdf ( drb, cdf2 )
      cdf = 1.0D+00 - cdf2 + 2.0D+00 * ( sum2 + tfn ( drb, a ) )
    end if
!
!  Normal approximation.
!
  else

    a = sqrt ( 0.5D+00 * f ) * exp ( gamma_log ( 0.5D+00 * ( f - 1.0D+00 ) ) &
      - gamma_log ( 0.5D+00 * f ) ) * d

    temp = ( x - a ) / sqrt ( f * ( 1.0D+00 + d**2 ) / ( f - 2.0D+00 ) - a**2 )

    call normal_01_cdf ( temp, cdf2 )
    cdf = cdf2

  end if

  return
end
subroutine student_noncentral_cdf_values ( n_data, df, lambda, x, fx )

!*****************************************************************************80
!
!! STUDENT_NONCENTRAL_CDF_VALUES returns values of the noncentral Student CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NoncentralStudentTDistribution [ df, lambda ]
!      CDF [ dist, x ]
!
!    Mathematica seems to have some difficulty computing this function
!    to the desired number of digits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) DF, real ( kind = 8 ) LAMBDA, the parameters of the
!    function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30

  integer ( kind = 4 ) df
  integer ( kind = 4 ), save, dimension ( n_max ) :: df_vec = (/ &
     1,  2,  3, &
     1,  2,  3, &
     1,  2,  3, &
     1,  2,  3, &
     1,  2,  3, &
    15, 20, 25, &
     1,  2,  3, &
    10, 10, 10, &
    10, 10, 10, &
    10, 10, 10 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.8975836176504333D+00, &
    0.9522670169D+00, &
    0.9711655571887813D+00, &
    0.8231218864D+00, &
    0.9049021510D+00, &
    0.9363471834D+00, &
    0.7301025986D+00, &
    0.8335594263D+00, &
    0.8774010255D+00, &
    0.5248571617D+00, &
    0.6293856597D+00, &
    0.6800271741D+00, &
    0.20590131975D+00, &
    0.2112148916D+00, &
    0.2074730718D+00, &
    0.9981130072D+00, &
    0.9994873850D+00, &
    0.9998391562D+00, &
    0.168610566972D+00, &
    0.16967950985D+00, &
    0.1701041003D+00, &
    0.9247683363D+00, &
    0.7483139269D+00, &
    0.4659802096D+00, &
    0.9761872541D+00, &
    0.8979689357D+00, &
    0.7181904627D+00, &
    0.9923658945D+00, &
    0.9610341649D+00, &
    0.8688007350D+00 /)
  real ( kind = 8 ) lambda
  real ( kind = 8 ), save, dimension ( n_max ) :: lambda_vec = (/ &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00, &
    4.0D+00, &
    4.0D+00, &
    4.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
    15.00D+00, &
    15.00D+00, &
    15.00D+00, &
     0.05D+00, &
     0.05D+00, &
     0.05D+00, &
     4.00D+00, &
     4.00D+00, &
     4.00D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00, &
     6.00D+00, &
     6.00D+00, &
     6.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    df = 0
    lambda = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    df = df_vec(n_data)
    lambda = lambda_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function tfn ( h, a )

!*****************************************************************************80
!
!! TFN calculates the T function of Owen.
!
!  Discussion:
!
!    Owen's T function is useful for computation of the bivariate normal
!    distribution and the distribution of a skewed normal distribution.
!
!    Although it was originally formulated in terms of the bivariate
!    normal function, the function can be defined more directly as
!
!      T(H,A) = 1 / ( 2 * pi ) *
!        Integral ( 0 <= X <= A ) e^( -H^2 * (1+X^2) / 2 ) / (1+X^2) dX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by J C Young, C E Minder.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Donald Owen,
!    Tables for computing the bivariate normal distribution,
!    Annals of Mathematical Statistics,
!    Volume 27, pages 1075-1090, 1956.
!
!    JC Young, CE Minder,
!    Algorithm AS 76,
!    An Algorithm Useful in Calculating Non-Central T and
!    Bivariate Normal Distributions,
!    Applied Statistics,
!    Volume 23, Number 3, 1974, pages 455-457.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) H, A, the arguments of the T function.
!
!    Output, real ( kind = 8 ) TFN, the value of the T function.
!
  implicit none

  integer ( kind = 4 ), parameter :: ngauss = 10

  real ( kind = 8 ) a
  real ( kind = 8 ) as
  real ( kind = 8 ) h
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  real ( kind = 8 ) hs
  integer ( kind = 4 ) i
  real ( kind = 8 ) rt
  real ( kind = 8 ) tfn
  real ( kind = 8 ), parameter :: two_pi_inverse = 0.1591549430918953D+00
  real ( kind = 8 ), parameter :: tv1 = 1.0D-35
  real ( kind = 8 ), parameter :: tv2 = 15.0D+00
  real ( kind = 8 ), parameter :: tv3 = 15.0D+00
  real ( kind = 8 ), parameter :: tv4 = 1.0D-05
  real ( kind = 8 ), parameter, dimension ( ngauss ) :: weight = (/ &
    0.666713443086881375935688098933D-01, &
    0.149451349150580593145776339658D+00, &
    0.219086362515982043995534934228D+00, &
    0.269266719309996355091226921569D+00, &
    0.295524224714752870173892994651D+00, &
    0.295524224714752870173892994651D+00, &
    0.269266719309996355091226921569D+00, &
    0.219086362515982043995534934228D+00, &
    0.149451349150580593145776339658D+00, &
    0.666713443086881375935688098933D-01 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter, dimension ( ngauss ) :: xtab = (/ &
   -0.973906528517171720077964012084D+00, &
   -0.865063366688984510732096688423D+00, &
   -0.679409568299024406234327365115D+00, &
   -0.433395394129247190799265943166D+00, &
   -0.148874338981631210884826001130D+00, &
    0.148874338981631210884826001130D+00, &
    0.433395394129247190799265943166D+00, &
    0.679409568299024406234327365115D+00, &
    0.865063366688984510732096688423D+00, &
    0.973906528517171720077964012084D+00 /)
!
!  Test for H near zero.
!
  if ( abs ( h ) < tv1 ) then
    tfn = atan ( a ) * two_pi_inverse
!
!  Test for large values of abs(H).
!
  else if ( tv2 < abs ( h ) ) then
    tfn = 0.0D+00
!
!  Test for A near zero.
!
  else if ( abs ( a ) < tv1 ) then
    tfn = 0.0D+00
!
!  Test whether abs(A) is so large that it must be truncated.
!  If so, the truncated value of A is H2.
!
  else

    hs = - 0.5D+00 * h * h
    h2 = a
    as = a * a
!
!  Computation of truncation point by Newton iteration.
!
    if ( tv3 <= log ( 1.0D+00 + as ) - hs * as ) then

      h1 = 0.5D+00 * a
      as = 0.25D+00 * as

      do

        rt = as + 1.0D+00
        h2 = h1 + ( hs * as + tv3 - log ( rt ) ) &
          / ( 2.0D+00 * h1 * ( 1.0D+00 / rt - hs ) )
        as = h2 * h2

        if ( abs ( h2 - h1 ) < tv4 ) then
          exit
        end if

        h1 = h2

      end do

    end if
!
!  Gaussian quadrature on the interval [0,H2].
!
    rt = 0.0D+00
    do i = 1, ngauss
      x = 0.5D+00 * h2 * ( xtab(i) + 1.0D+00 )
      rt = rt + weight(i) * exp ( hs * ( 1.0D+00 + x * x ) ) &
        / ( 1.0D+00 + x * x )
    end do

    tfn = rt * ( 0.5D+00 * h2 ) * two_pi_inverse

  end if

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2003
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

  character ( len = 40 ) string

  call timestring ( string )

  write ( *, '(a)' ) trim ( string )

  return
end
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
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
  character ( len = * ) string
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

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

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine triangle_cdf ( x, a, b, c, cdf )

!*****************************************************************************80
!
!! TRIANGLE_CDF evaluates the Triangle CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A <= B <= C and A < C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x <= a ) then

    cdf = 0.0D+00

  else if ( x <= b ) then

    if ( a == b ) then
      cdf = 0.0D+00
    else
      cdf = ( x - a ) * ( x - a ) / ( b - a ) / ( c - a )
    end if

  else if ( x <= c ) then

    cdf = ( b - a ) / ( c - a ) &
        + ( 2.0D+00 * c - b - x ) * ( x - b ) / ( c - b ) / ( c - a )

  else

    cdf = 1.0D+00

  end if

  return
end
subroutine triangle_cdf_inv ( cdf, a, b, c, x )

!*****************************************************************************80
!
!! TRIANGLE_CDF_INV inverts the Triangle CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A <= B <= C and A < C.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf_mid
  real ( kind = 8 ) d
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  d = 2.0D+00 / ( c - a )
  cdf_mid = 0.5D+00 * d * ( b - a )

  if ( cdf <= cdf_mid ) then
    x = a + sqrt ( cdf * ( b - a ) * ( c - a ) )
  else
    x = c - sqrt ( ( c - b ) * ( ( c - b ) - ( cdf - cdf_mid ) * ( c - a ) ) )
  end if

  return
end
function triangle_check ( a, b, c )

!*****************************************************************************80
!
!! TRIANGLE_CHECK checks the parameters of the Triangle CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A <= B <= C and A < C.
!
!    Output, logical TRIANGLE_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical triangle_check

  if ( b < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B < A.'
    triangle_check = .false.
    return
  end if

  if ( c < b ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C < B.'
    triangle_check = .false.
    return
  end if

  if ( a == c ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A == C.'
    triangle_check = .false.
    return
  end if

  triangle_check = .true.

  return
end
subroutine triangle_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! TRIANGLE_MEAN returns the mean of the Triangle PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A <= B <= C and A < C.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the discrete uniform PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) mean

  mean = a + ( c + b - 2.0D+00 * a ) / 3.0D+00

  return
end
subroutine triangle_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! TRIANGLE_PDF evaluates the Triangle PDF.
!
!  Discussion:
!
!    Given points A <= B <= C, the probability is 0 to the left of A,
!    rises linearly to a maximum of 2/(C-A) at B, drops linearly to zero
!    at C, and is zero for all values greater than C.
!
!    PDF(A,B,C;X)
!      = 2 * ( X - A ) / ( B - A ) / ( C - A ) for A <= X <= B
!      = 2 * ( C - X ) / ( C - B ) / ( C - A ) for B <= X <= C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A <= B <= C and A < C.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x <= a ) then

    pdf = 0.0D+00

  else if ( x <= b ) then

    if ( a == b ) then
      pdf = 0.0D+00
    else
      pdf = 2.0D+00 * ( x - a ) / ( b - a ) / ( c - a )
    end if

  else if ( x <= c ) then

    if ( b == c ) then
      pdf = 0.0D+00
    else
      pdf = 2.0D+00 * ( c - x ) / ( c - b ) / ( c - a )
    end if

  else
    pdf = 0.0D+00
  end if

  return
end
subroutine triangle_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! TRIANGLE_SAMPLE samples the Triangle PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A <= B <= C and A < C.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call triangle_cdf_inv ( cdf, a, b, c, x )

  return
end
subroutine triangle_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! TRIANGLE_VARIANCE returns the variance of the Triangle PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    A <= B <= C and A < C.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) variance

  variance = ( ( c - a ) * ( c - a ) &
             - ( c - a ) * ( b - a ) &
             + ( b - a ) * ( b - a ) ) / 18.0D+00

  return
end
subroutine triangular_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! TRIANGULAR_CDF evaluates the Triangular CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x <= a ) then
    cdf = 0.0D+00
  else if ( x <= 0.5D+00 * ( a + b ) ) then
    cdf = 2.0D+00 * ( x**2 - 2.0D+00 * a * x + a**2 ) / ( b - a )**2
  else if ( x <= b ) then
    cdf = 0.5D+00 + ( - 2.0D+00 * x**2 + 4.0D+00 * b * x + 0.5D+00 * a**2 &
      - a * b - 1.5D+00 * b**2 ) / ( b - a )**2
  else
    cdf = 1.0D+00
  end if

  return
end
subroutine triangular_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! TRIANGULAR_CDF_INV inverts the Triangular CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULAR_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf <= 0.5D+00 ) then
    x = a + 0.5D+00 * ( b - a ) * sqrt ( 2.0D+00 * cdf )
  else
    x = b - 0.5D+00 * ( b - a ) * sqrt ( 2.0D+00 * ( 1.0D+00 - cdf ) )
  end if

  return
end
function triangular_check ( a, b )

!*****************************************************************************80
!
!! TRIANGULAR_CHECK checks the parameters of the Triangular CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Output, logical TRIANGULAR_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical triangular_check

  if ( b <= a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULAR_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= A.'
    triangular_check = .false.
    return
  end if

  triangular_check = .true.

  return
end
subroutine triangular_mean ( a, b, mean )

!*****************************************************************************80
!
!! TRIANGULAR_MEAN returns the mean of the Triangular PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = 0.5D+00 * ( a + b )

  return
end
subroutine triangular_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! TRIANGULAR_PDF evaluates the Triangular PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = 4 * ( X - A ) / ( B - A )^2 for A <= X <= (A+B)/2
!               = 4 * ( B - X ) / ( B - A )^2 for (A+B)/2 <= X <= B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x <= a ) then
    pdf = 0.0D+00
  else if ( x <= 0.5D+00 * ( a + b ) ) then
    pdf = 4.0D+00 * ( x - a ) / ( b - a )**2
  else if ( x <= b ) then
    pdf = 4.0D+00 * ( b - x ) / ( b - a )**2
  else
    pdf = 0.0D+00
  end if

  return
end
subroutine triangular_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! TRIANGULAR_SAMPLE samples the Triangular PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call triangular_cdf_inv ( cdf, a, b, x )

  return
end
subroutine triangular_variance ( a, b, variance )

!*****************************************************************************80
!
!! TRIANGULAR_VARIANCE returns the variance of the Triangular PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = ( b - a )**2 / 24.0D+00

  return
end
function trigamma ( x )

!*****************************************************************************80
!
!! TRIGAMMA calculates the TriGamma function.
!
!  Discussion:
!
!    TriGamma(x) = d^2 log ( Gamma ( x ) ) / dx^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    FORTRAN77 original version by B Schneider
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    BE Schneider,
!    Algorithm AS 121:
!    Trigamma Function,
!    Applied Statistics,
!    Volume 27, Number 1, page 97-99, 1978.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the trigamma function.
!    0 < X.
!
!    Output, real ( kind = 8 ) TRIGAMMA, the value of the
!    trigamma function at X.
!
  implicit none

  real ( kind = 8 ), parameter :: a = 0.0001D+00
  real ( kind = 8 ), parameter :: b = 5.0D+00
  real ( kind = 8 ), parameter :: b2 =   1.0D+00 / 6.0D+00
  real ( kind = 8 ), parameter :: b4 = - 1.0D+00 / 30.0D+00
  real ( kind = 8 ), parameter :: b6 =   1.0D+00 / 42.0D+00
  real ( kind = 8 ), parameter :: b8 = - 1.0D+00 / 30.0D+00
  real ( kind = 8 ) trigamma
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!
!  1): If X is not positive, fail.
!
  if ( x <= 0.0D+00 ) then

    trigamma = 0.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIGAMMA - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
!
!  2): If X is smaller than A, use a small value approximation.
!
  else if ( x <= a ) then

    trigamma = 1.0D+00 / x**2
!
!  3): Otherwise, increase the argument to B <= ( X + I ).
!
  else

    z = x
    trigamma = 0.0D+00

    do while ( z < b )
      trigamma = trigamma + 1.0D+00 / z**2
      z = z + 1.0D+00
    end do
!
!  ...and then apply an asymptotic formula.
!
    y = 1.0D+00 / z**2

    trigamma = trigamma + 0.5D+00 * &
            y + ( 1.0D+00 &
          + y * ( b2 &
          + y * ( b4 &
          + y * ( b6 &
          + y *   b8 )))) / z

  end if

  return
end
subroutine uniform_01_cdf ( x, cdf )

!*****************************************************************************80
!
!! UNIFORM_01_CDF evaluates the Uniform 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    cdf = 0.0D+00
  else if ( 1.0D+00 < x ) then
    cdf = 1.0D+00
  else
    cdf = x
  end if

  return
end
subroutine uniform_01_cdf_inv ( cdf, x )

!*****************************************************************************80
!
!! UNIFORM_01_CDF_INV inverts the Uniform 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNIFORM_01_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = cdf

  return
end
subroutine uniform_01_mean ( mean )

!*****************************************************************************80
!
!! UNIFORM_01_MEAN returns the mean of the Uniform 01 PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) mean

  mean = 0.5D+00

  return
end
subroutine uniform_01_order_sample ( n, seed, x )

!*****************************************************************************80
!
!! UNIFORM_01_ORDER_SAMPLE samples the Uniform 01 Order PDF.
!
!  Discussion:
!
!    In effect, this routine simply generates N samples of the
!    Uniform 01 PDF; but it generates them in order.  (Actually,
!    it generates them in descending order, but stores them in
!    the array in ascending order).  This saves the work of
!    sorting the results.  Moreover, if the order statistics
!    for another PDF are desired, and the inverse CDF is available,
!    then the desired values may be generated, presorted, by
!    calling this routine and using the results as input to the
!    inverse CDF routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 168.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the sample.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), N samples of the Uniform 01 PDF, in
!    ascending order.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x(n)

  v = 1.0D+00
  do i = n, 1, -1
    u = r8_uniform_01 ( seed )
    v = v * u**( 1.0D+00 / real ( i, kind = 8 ) )
    x(i) = v
  end do

  return
end
subroutine uniform_01_pdf ( x, pdf )

!*****************************************************************************80
!
!! UNIFORM_01_PDF evaluates the Uniform 01 PDF.
!
!  Discussion:
!
!    PDF(X) = 1 for 0 <= X <= 1
!           = 0 otherwise
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X <= 1.0.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
    pdf = 0.0D+00
  else
    pdf = 1.0D+00
  end if

  return
end
function uniform_01_sample ( seed )

!*****************************************************************************80
!
!! UNIFORM_01_SAMPLE is a portable random number generator.
!
!  Discussion:
!
!    SEED = SEED * (7**5) mod ( 2**31 - 1 )
!    UNIFORM_01_SAMPLE = SEED * / ( 2**31 - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the integer "seed" used to generate
!    the output random number, and updated in preparation for the
!    next one.  SEED should not be zero.
!
!    Output, real ( kind = 8 ) UNIFORM_01_SAMPLE, a random value between 0 and 1.
!
!  Local parameters:
!
!    IA = 7**5
!    IB = 2**15
!    IB16 = 2**16
!    IP = 2**31-1
!
  implicit none

  integer ( kind = 4 ), parameter :: ia = 16807
  integer ( kind = 4 ), parameter :: ib15 = 32768
  integer ( kind = 4 ), parameter :: ib16 = 65536
  integer ( kind = 4 ), parameter :: ip = 2147483647
  integer ( kind = 4 ) iprhi
  integer ( kind = 4 ) ixhi
  integer ( kind = 4 ) k
  integer ( kind = 4 ) leftlo
  integer ( kind = 4 ) loxa
  integer ( kind = 4 ) seed
  real ( kind = 8 ) uniform_01_sample
!
!  Don't let SEED be 0 or IP
!
  if ( seed == 0 .or. seed == ip ) then
    seed = ip / 2
  end if
!
!  Get the 15 high order bits of SEED.
!
  ixhi = seed / ib16
!
!  Get the 16 low bits of SEED and form the low product.
!
  loxa = ( seed - ixhi * ib16 ) * ia
!
!  Get the 15 high order bits of the low product.
!
  leftlo = loxa / ib16
!
!  Form the 31 highest bits of the full product.
!
  iprhi = ixhi * ia + leftlo
!
!  Get overflow past the 31st bit of full product.
!
  k = iprhi / ib15
!
!  Assemble all the parts and presubtract IP.  The parentheses are essential.
!
  seed = ( ( ( loxa - leftlo * ib16 ) - ip ) + ( iprhi - k * ib15 ) * ib16 ) &
    + k
!
!  Add IP back in if necessary.
!
  if ( seed < 0 ) then
    seed = seed + ip
  end if
!
!  Multiply by 1 / (2**31-1).
!
  uniform_01_sample = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine uniform_01_variance ( variance )

!*****************************************************************************80
!
!! UNIFORM_01_VARIANCE returns the variance of the Uniform 01 PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) variance

  variance = 1.0D+00 / 12.0D+00

  return
end
subroutine uniform_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! UNIFORM_CDF evaluates the Uniform CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( x < a ) then
    cdf = 0.0D+00
  else if ( b < x ) then
    cdf = 1.0D+00
  else
    cdf = ( x - a ) / ( b - a )
  end if

  return
end
subroutine uniform_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! UNIFORM_CDF_INV inverts the Uniform CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNIFORM_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a + ( b - a ) * cdf

  return
end
function uniform_check ( a, b )

!*****************************************************************************80
!
!! UNIFORM_CHECK checks the parameters of the Uniform CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Output, logical UNIFORM_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical uniform_check

  if ( b <= a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNIFORM_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= A.'
    uniform_check = .false.
    return
  end if

  uniform_check = .true.

  return
end
subroutine uniform_mean ( a, b, mean )

!*****************************************************************************80
!
!! UNIFORM_MEAN returns the mean of the Uniform PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = 0.5D+00 * ( a + b )

  return
end
subroutine uniform_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! UNIFORM_PDF evaluates the Uniform PDF.
!
!  Discussion:
!
!    The Uniform PDF is also known as the "Rectangular" or "de Moivre" PDF.
!
!    PDF(A,B;X) = 1 / ( B - A ) for A <= X <= B
!               = 0 otherwise
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  if ( x < a .or. b < x ) then
    pdf = 0.0D+00
  else
    pdf = 1.0D+00 / ( b - a )
  end if

  return
end
subroutine uniform_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! UNIFORM_SAMPLE samples the Uniform PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call uniform_cdf_inv ( cdf, a, b, x )

  return
end
subroutine uniform_variance ( a, b, variance )

!*****************************************************************************80
!
!! UNIFORM_VARIANCE returns the variance of the Uniform PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) variance

  variance = ( b - a )**2 / 12.0D+00

  return
end
subroutine uniform_discrete_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! UNIFORM_DISCRETE_CDF evaluates the Uniform Discrete CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the CDF.
!
!    Input, integer ( kind = 4 ) A, B, the parameters of the PDF.
!    A <= B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) x

  if ( x < a ) then
    cdf = 0.0D+00
  else if ( b < x ) then
    cdf = 1.0D+00
  else
    cdf = real ( x + 1 - a, kind = 8 ) / real ( b + 1 - a, kind = 8 )
  end if

  return
end
subroutine uniform_discrete_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! UNIFORM_DISCRETE_CDF_INV inverts the Uniform Discrete CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, integer ( kind = 4 ) A, B, the parameters of the PDF.
!    A <= B.
!
!    Output, integer ( kind = 4 ) X, the smallest argument whose CDF is greater
!    than or equal to CDF.
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) a2
  integer ( kind = 4 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) x
  real ( kind = 8 ) x2

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNIFORM_DISCRETE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  a2 = real ( a, kind = 8 ) - 0.5D+00
  b2 = real ( b, kind = 8 ) + 0.5D+00
  x2 = a + cdf * ( b2 - a2 )

  x = nint ( x2 )

  x = max ( x, a )
  x = min ( x, b )

  return
end
function uniform_discrete_check ( a, b )

!*****************************************************************************80
!
!! UNIFORM_DISCRETE_CHECK checks the parameters of the Uniform discrete CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the parameters of the PDF.
!    A <= B.
!
!    Output, logical UNIFORM_DISCRETE_CHECK, is true if the parameters
!    are legal.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  logical uniform_discrete_check

  if ( b < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNIFORM_DISCRETE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B < A.'
    uniform_discrete_check = .false.
    return
  end if

  uniform_discrete_check = .true.

  return
end
subroutine uniform_discrete_mean ( a, b, mean )

!*****************************************************************************80
!
!! UNIFORM_DISCRETE_MEAN returns the mean of the Uniform discrete PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the parameters of the PDF.
!    A <= B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) mean

  mean = 0.5D+00 * real ( a + b, kind = 8 )

  return
end
subroutine uniform_discrete_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! UNIFORM_DISCRETE_PDF evaluates the Uniform discrete PDF.
!
!  Discussion:
!
!    The Uniform Discrete PDF is also known as the "Rectangular"
!    Discrete PDF.
!
!    PDF(A,B;X) = 1 / ( B + 1 - A ) for A <= X <= B.
!
!    The parameters define the interval of integers
!    for which the PDF is nonzero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the PDF.
!
!    Input, integer ( kind = 4 ) A, B, the parameters of the PDF.
!    A <= B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  if ( x < a .or. b < x ) then
    pdf = 0.0D+00
  else
    pdf = 1.0D+00 / real ( b + 1 - a, kind = 8 )
  end if

  return
end
subroutine uniform_discrete_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! UNIFORM_DISCRETE_SAMPLE samples the Uniform discrete PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the parameters of the PDF.
!    A <= B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x

  cdf = r8_uniform_01 ( seed )

  call uniform_discrete_cdf_inv ( cdf, a, b, x )

  return
end
subroutine uniform_discrete_variance ( a, b, variance )

!*****************************************************************************80
!
!! UNIFORM_DISCRETE_VARIANCE returns the variance of the Uniform discrete PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the parameters of the PDF.
!    A <= B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) variance

  variance = real ( ( b + 1 - a )**2 - 1, kind = 8 ) / 12.0D+00

  return
end
subroutine uniform_nsphere_sample ( n, seed, x )

!*****************************************************************************80
!
!! UNIFORM_NSPHERE_SAMPLE samples the Uniform Unit Sphere PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 168.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the sphere.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a point on the unit N sphere, chosen
!    with a uniform probability.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  do i = 1, n
    call normal_01_sample ( seed, x(i) )
  end do

  x(1:n) = x(1:n) / sqrt ( sum ( x(1:n)**2 ) )

  return
end
subroutine von_mises_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! VON_MISES_CDF evaluates the von Mises CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    Original FORTRAN77 version by Geoffrey Hill.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Geoffrey Hill,
!    Algorithm 518,
!    Incomplete Bessel Function I0: The von Mises Distribution,
!    ACM Transactions on Mathematical Software,
!    Volume 3, Number 3, September 1977, pages 279-284.
!
!    Kanti Mardia, Peter Jupp,
!    Directional Statistics,
!    Wiley, 2000,
!    QA276.M335
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!    A - PI <= X <= A + PI.
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0D+00 <= B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), parameter :: a1 = 12.0D+00
  real ( kind = 8 ), parameter :: a2 = 0.8D+00
  real ( kind = 8 ), parameter :: a3 = 8.0D+00
  real ( kind = 8 ), parameter :: a4 = 1.0D+00
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: c1 = 56.0D+00
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: ck = 10.5D+00
  real ( kind = 8 ) cn
  real ( kind = 8 ) error_f
  real ( kind = 8 ) erfx
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) sn
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!
!  We expect -PI <= X - A <= PI.
!
  if ( x - a <= -pi ) then
    cdf = 0.0D+00
    return
  end if

  if ( pi <= x - a ) then
    cdf = 1.0D+00
    return
  end if
!
!  Convert the angle (X - A) modulo 2 PI to the range ( 0, 2 * PI ).
!
  z = b

  u = mod ( x - a + pi, 2.0D+00 * pi )

  if ( u < 0.0D+00 ) then
    u = u + 2.0D+00 * pi
  end if

  y = u - pi
!
!  For small B, sum IP terms by backwards recursion.
!
  if ( z <= ck ) then

    v = 0.0D+00

    if ( 0.0D+00 < z ) then

      ip = int ( z * a2 - a3 / ( z + a4 ) + a1 )
      p = real ( ip, kind = 8 )
      s = sin ( y )
      c = cos ( y )
      y = p * y
      sn = sin ( y )
      cn = cos ( y )
      r = 0.0D+00
      z = 2.0D+00 / z

      do n = 2, ip
        p = p - 1.0D+00
        y = sn
        sn = sn * c - cn * s
        cn = cn * c + y * s
        r = 1.0D+00 / ( p * z + r )
        v = ( sn / p + v ) * r
      end do

    end if

    cdf = ( u * 0.5D+00 + v ) / pi
!
!  For large B, compute the normal approximation and left tail.
!
  else

    c = 24.0D+00 * z
    v = c - c1
    r = sqrt ( ( 54.0D+00 / ( 347.0D+00 / v + 26.0D+00 - c ) - 6.0D+00 + c ) &
      / 12.0D+00 )
    z = sin ( 0.5D+00 * y ) * r
    s = 2.0D+00 * z**2
    v = v - s + 3.0D+00
    y = ( c - s - s - 16.0D+00 ) / 3.0D+00
    y = ( ( s + 1.75D+00 ) * s + 83.5D+00 ) / v - y
    arg = z * ( 1.0D+00 - s / y**2 )
    erfx = error_f ( arg )
    cdf = 0.5D+00 * erfx + 0.5D+00

  end if

  cdf = max ( cdf, 0.0D+00 )
  cdf = min ( cdf, 1.0D+00 )

  return
end
subroutine von_mises_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! VON_MISES_CDF_INV inverts the von Mises CDF.
!
!  Discussion:
!
!    A simple bisection method is used on the interval [ A - PI, A + PI ].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0D+00 <= B.
!
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!    A - PI <= X <= A + PI.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf1
  real ( kind = 8 ) cdf2
  real ( kind = 8 ) cdf3
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: tol = 0.000001D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VON_MISES_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf == 0.0D+00 ) then
    x = a - pi
    return
  else if ( 1.0D+00 == cdf ) then
    x = a + pi
    return
  end if

  x1 = a - pi
  cdf1 = 0.0D+00

  x2 = a + pi
  cdf2 = 1.0D+00
!
!  Now use bisection.
!
  it = 0

  do

    it = it + 1

    x3 = 0.5D+00 * ( x1 + x2 )
    call von_mises_cdf ( x3, a, b, cdf3 )

    if ( abs ( cdf3 - cdf ) < tol ) then
      x = x3
      exit
    end if

    if ( it_max < it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'VON_MISES_CDF_INV - Fatal error!'
      write ( *, '(a)' ) '  Iteration limit exceeded.'
      stop
    end if

    if ( sign ( 1.0D+00, cdf3 - cdf ) == sign ( 1.0D+00, cdf1 - cdf ) ) then
      x1 = x3
      cdf1 = cdf3
    else
      x2 = x3
      cdf2 = cdf3
    end if

  end do

  return
end
subroutine von_mises_cdf_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! VON_MISES_CDF_VALUES returns some values of the von Mises CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kanti Mardia, Peter Jupp,
!    Directional Statistics,
!    Wiley, 2000,
!    QA276.M335
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Output, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0D+00 <= B.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 23

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.0D+00, &
     0.0D+00, &
     0.0D+00, &
     0.0D+00, &
     0.0D+00, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
    -0.2D+01, &
    -0.1D+01, &
     0.0D+01, &
     0.1D+01, &
     0.2D+01, &
     0.3D+01, &
     0.0D+00, &
     0.0D+00, &
     0.0D+00, &
     0.0D+00, &
     0.0D+00, &
     0.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.1D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01, &
     0.2D+01, &
     0.3D+01, &
     0.3D+01, &
     0.3D+01, &
     0.3D+01, &
     0.3D+01, &
     0.3D+01, &
     0.0D+00, &
     0.1D+01, &
     0.2D+01, &
     0.3D+01, &
     0.4D+01, &
     0.5D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.2535089956281180D-01, &
    0.1097539041177346D+00, &
    0.5000000000000000D+00, &
    0.8043381312498558D+00, &
    0.9417460124555197D+00, &
    0.5000000000000000D+00, &
    0.6018204118446155D+00, &
    0.6959356933122230D+00, &
    0.7765935901304593D+00, &
    0.8410725934916615D+00, &
    0.8895777369550366D+00, &
    0.9960322705517925D+00, &
    0.9404336090170247D+00, &
    0.5000000000000000D+00, &
    0.5956639098297530D-01, &
    0.3967729448207649D-02, &
    0.2321953958111930D-03, &
    0.6250000000000000D+00, &
    0.7438406999109122D+00, &
    0.8369224904294019D+00, &
    0.8941711407897124D+00, &
    0.9291058600568743D+00, &
    0.9514289900655436D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -0.2617993977991494D+01, &
    -0.1570796326794897D+01, &
     0.0000000000000000D+00, &
     0.1047197551196598D+01, &
     0.2094395102393195D+01, &
     0.1000000000000000D+01, &
     0.1200000000000000D+01, &
     0.1400000000000000D+01, &
     0.1600000000000000D+01, &
     0.1800000000000000D+01, &
     0.2000000000000000D+01, &
     0.0000000000000000D+00, &
     0.0000000000000000D+00, &
     0.0000000000000000D+00, &
     0.0000000000000000D+00, &
     0.0000000000000000D+00, &
     0.0000000000000000D+00, &
     0.7853981633974483D+00, &
     0.7853981633974483D+00, &
     0.7853981633974483D+00, &
     0.7853981633974483D+00, &
     0.7853981633974483D+00, &
     0.7853981633974483D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function von_mises_check ( a, b )

!*****************************************************************************80
!
!! VON_MISES_CHECK checks the parameters of the von Mises PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0D+00 <= B.
!
!    Output, logical VON_MISES_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  logical von_mises_check

  if ( a < -pi .or. pi < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VON_MISES_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A < -PI or PI < A.'
    von_mises_check = .false.
    return
  end if

  if ( b < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VON_MISES_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B < 0.0'
    von_mises_check = .false.
    return
  end if

  von_mises_check = .true.

  return
end
subroutine von_mises_circular_variance ( a, b, circular_variance )

!*****************************************************************************80
!
!! VON_MISES_CIRCULAR_VARIANCE returns the circular variance of the von Mises PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0D+00 <= B.
!
!    Output, real ( kind = 8 ) CIRCULAR_VARIANCE, the circular variance
!    of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bessel_i0
  real ( kind = 8 ) bessel_i1
  real ( kind = 8 ) circular_variance

  circular_variance = 1.0D+00 - bessel_i1 ( b ) / bessel_i0 ( b )

  return
end
subroutine von_mises_mean ( a, b, mean )

!*****************************************************************************80
!
!! VON_MISES_MEAN returns the mean of the von Mises PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0D+00 <= B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) mean

  mean = a

  return
end
subroutine von_mises_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! VON_MISES_PDF evaluates the von Mises PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = EXP ( B * COS ( X - A ) ) / ( 2 * PI * I0(B) )
!
!    where:
!
!      I0(*) is the modified Bessel function of the first
!      kind of order 0.
!
!    The von Mises distribution for points on the unit circle is
!    analogous to the normal distribution of points on a line.
!    The variable X is interpreted as a deviation from the angle A,
!    with B controlling the amount of dispersion.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 160.
!
!    Donald Best, Nicholas Fisher,
!    Efficient Simulation of the von Mises Distribution,
!    Applied Statistics,
!    Volume 28, Number 2, pages 152-157.
!
!    Merran Evans, Nicholas Hastings, Brian Peacock,
!    Statistical Distributions,
!    Wiley, 2000,
!    LC: QA273.6.E92, pages 189-191.
!
!    Kanti Mardia, Peter Jupp,
!    Directional Statistics,
!    Wiley, 2000,
!    LC: QA276.M335
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A - PI <= X <= A + PI.
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0D+00 <= B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bessel_i0
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( x < a - pi ) then
    pdf = 0.0D+00
  else if ( x <= a + pi ) then
    pdf = exp ( b * cos ( x - a ) ) / ( 2.0D+00 * pi * bessel_i0 ( b ) )
  else
    pdf = 0.0D+00
  end if

  return
end
subroutine von_mises_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! VON_MISES_SAMPLE samples the von Mises PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Best, Nicholas Fisher,
!    Efficient Simulation of the von Mises Distribution,
!    Applied Statistics,
!    Volume 28, Number 2, pages 152-157.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0D+00 <= B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) f
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) rho
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tau
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) u3
  real ( kind = 8 ) x
  real ( kind = 8 ) z

  tau = 1.0D+00 + sqrt ( 1.0D+00 + 4.0D+00 * b * b )
  rho = ( tau - sqrt ( 2.0D+00 * tau ) ) / ( 2.0D+00 * b )
  r = ( 1.0D+00 + rho**2 ) / ( 2.0D+00 * rho )

  do

    u1 = r8_uniform_01 ( seed )
    z = cos ( pi * u1 )
    f = ( 1.0D+00 + r * z ) / ( r + z )
    c = b * ( r - f )

    u2 = r8_uniform_01 ( seed )

    if ( u2 < c * ( 2.0D+00 - c ) ) then
      exit
    end if

    if ( c <= log ( c / u2 ) + 1.0D+00 ) then
      exit
    end if

  end do

  u3 = r8_uniform_01 ( seed )

  x = a + sign ( 1.0D+00, u3 - 0.5D+00 ) * acos ( f )

  return
end
subroutine weibull_cdf ( x, a, b, c, cdf )

!*****************************************************************************80
!
!! WEIBULL_CDF evaluates the Weibull CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!    A <= X.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x < a ) then
    cdf = 0.0D+00
  else
    y = ( x - a ) / b
    cdf = 1.0D+00 - 1.0D+00 / exp ( y**c )
  end if

  return
end
subroutine weibull_cdf_inv ( cdf, a, b, c, x )

!*****************************************************************************80
!
!! WEIBULL_CDF_INV inverts the Weibull CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 < CDF < 1.0.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEIBULL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a + b * ( - log ( 1.0D+00 - cdf ) )**( 1.0D+00 / c )

  return
end
subroutine weibull_cdf_values ( n_data, alpha, beta, x, fx )

!*****************************************************************************80
!
!! WEIBULL_CDF_VALUES returns some values of the Weibull CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = WeibullDistribution [ alpha, beta ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    LC: QA47.A34,
!    ISBN: 0-486-61272-4.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) ALPHA, the first parameter of the distribution.
!
!    Output, real ( kind = 8 ) BETA, the second parameter of the distribution.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) alpha
  real ( kind = 8 ), save, dimension ( n_max ) :: alpha_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01 /)
  real ( kind = 8 ) beta
  real ( kind = 8 ), save, dimension ( n_max ) :: beta_vec = (/ &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.5000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.8646647167633873D+00, &
    0.9816843611112658D+00, &
    0.9975212478233336D+00, &
    0.9996645373720975D+00, &
    0.6321205588285577D+00, &
    0.4865828809674080D+00, &
    0.3934693402873666D+00, &
    0.3296799539643607D+00, &
    0.8946007754381357D+00, &
    0.9657818816883340D+00, &
    0.9936702845725143D+00, &
    0.9994964109502630D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.4000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01, &
    0.3000000000000000D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    alpha = 0.0D+00
    beta = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    alpha = alpha_vec(n_data)
    beta = beta_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function weibull_check ( a, b, c )

!*****************************************************************************80
!
!! WEIBULL_CHECK checks the parameters of the Weibull CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, logical WEIBULL_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical weibull_check

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEIBULL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    weibull_check = .false.
    return
  end if

  if ( c <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEIBULL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C <= 0.'
    weibull_check = .false.
    return
  end if

  weibull_check = .true.

  return
end
subroutine weibull_mean ( a, b, c, mean )

!*****************************************************************************80
!
!! WEIBULL_MEAN returns the mean of the Weibull PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) mean
  real ( kind = 8 ) r8_gamma

  mean = b * r8_gamma ( ( c + 1.0D+00 ) / c ) + a

  return
end
subroutine weibull_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! WEIBULL_PDF evaluates the Weibull PDF.
!
!  Discussion:
!
!    PDF(A,B,C;X) = ( C / B ) * ( ( X - A ) / B )**( C - 1 )
!     * EXP ( - ( ( X - A ) / B )**C ).
!
!    The Weibull PDF is also known as the Frechet PDF.
!
!    WEIBULL_PDF(A,B,1;X) is the Exponential PDF.
!
!    WEIBULL_PDF(0,1,2;X) is the Rayleigh PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A <= X
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x < a ) then

    pdf = 0.0D+00

  else

    y = ( x - a ) / b

    pdf = ( c / b ) * y**( c - 1.0D+00 )  / exp ( y**c )

  end if

  return
end
subroutine weibull_sample ( a, b, c, seed, x )

!*****************************************************************************80
!
!! WEIBULL_SAMPLE samples the Weibull PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  cdf = r8_uniform_01 ( seed )

  call weibull_cdf_inv ( cdf, a, b, c, x )

  return
end
subroutine weibull_variance ( a, b, c, variance )

!*****************************************************************************80
!
!! WEIBULL_VARIANCE returns the variance of the Weibull PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) variance

  g1 = r8_gamma ( ( c + 2.0D+00 ) / c )
  g2 = r8_gamma ( ( c + 1.0D+00 ) / c )

  variance = b * b * ( g1 - g2 * g2 )

  return
end
subroutine weibull_discrete_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! WEIBULL_DISCRETE_CDF evaluates the Discrete Weibull CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the CDF.
!    0 <= X.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A <= 1.0D+00,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    cdf = 0.0D+00
  else
    cdf = 1.0D+00 - ( 1.0D+00 - a )**((x+1)**b)
  end if

  return
end
subroutine weibull_discrete_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! WEIBULL_DISCRETE_CDF_INV inverts the Discrete Weibull CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A <= 1.0D+00,
!    0.0D+00 < B.
!
!    Output, integer ( kind = 4 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) r8_ceiling
  integer ( kind = 4 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEIBULL_DISCRETE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = r8_ceiling ( &
    ( log ( 1.0D+00 - cdf ) / log ( 1.0D+00 - a ) )**( 1.0D+00 / b ) - 1.0D+00 )

  return
end
function weibull_discrete_check ( a, b )

!*****************************************************************************80
!
!! WEIBULL_DISCRETE_CHECK checks the parameters of the discrete Weibull CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A <= 1.0D+00,
!    0.0D+00 < B.
!
!    Output, logical WEIBULL_DISCRETE_CHECK, is true if the parameters
!    are legal.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical weibull_discrete_check

  if ( a < 0.0D+00 .or. 1.0D+00 < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEIBULL_DISCRETE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A < 0 or 1 < A.'
    weibull_discrete_check = .false.
    return
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEIBULL_DISCRETE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    weibull_discrete_check = .false.
    return
  end if

  weibull_discrete_check = .true.

  return
end
subroutine weibull_discrete_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! WEIBULL_DISCRETE_PDF evaluates the discrete Weibull PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = ( 1 - A )**X**B - ( 1 - A )**(X+1)**B.
!
!    WEIBULL_DISCRETE_PDF(A,1;X) = GEOMETRIC_PDF(A;X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the PDF.
!    0 <= X
!
!    Input, real ( kind = 8 ) A, B, the parameters that define the PDF.
!    0 <= A <= 1,
!    0 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    pdf = 0.0D+00
  else
    pdf = ( 1.0D+00 - a )**(x**b) - ( 1.0D+00 - a )**((x+1)**b)
  end if

  return
end
subroutine weibull_discrete_sample ( a, b, seed, x )

!*****************************************************************************80
!
!! WEIBULL_DISCRETE_SAMPLE samples the discrete Weibull PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A <= 1.0D+00,
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x

  cdf = r8_uniform_01 ( seed )

  call weibull_discrete_cdf_inv ( cdf, a, b, x )

  return
end
function zeta ( p )

!*****************************************************************************80
!
!! ZETA estimates the Riemann Zeta function.
!
!  Discussion:
!
!    For 1 < P, the Riemann Zeta function is defined as:
!
!      ZETA ( P ) = Sum ( 1 <= N < Infinity ) 1 / N**P
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the power to which the integers are raised.
!    P must be greater than 1.  For integral P up to 20, a
!    precomputed value of ZETA is returned; otherwise the infinite
!    sum is approximated.
!
!    Output, real ( kind = 8 ) ZETA, an approximation to the Riemann
!    Zeta function.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) zsum
  real ( kind = 8 ) zsum_old
  real ( kind = 8 ) zeta

  if ( p <= 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ZETA - Fatal error!'
    write ( *, '(a)' ) '  Exponent P <= 1.0.'
    stop
  else if ( p == 2.0D+00 ) then
    zeta = pi**2 / 6.0D+00
  else if ( p == 3.0D+00 ) then
    zeta = 1.2020569032D+00
  else if ( p == 4.0D+00 ) then
    zeta = pi**4 / 90.0D+00
  else if ( p == 5.0D+00 ) then
    zeta = 1.0369277551D+00
  else if ( p == 6.0D+00 ) then
    zeta = pi**6 / 945.0D+00
  else if ( p == 7.0D+00 ) then
    zeta = 1.0083492774D+00
  else if ( p == 8.0D+00 ) then
    zeta = pi**8 / 9450.0D+00
  else if ( p == 9.0D+00 ) then
    zeta = 1.0020083928D+00
  else if ( p == 10.0D+00 ) then
    zeta = pi**10 / 93555.0D+00
  else if ( p == 11.0D+00 ) then
    zeta = 1.0004941886D+00
  else if ( p == 12.0D+00 ) then
    zeta = 1.0002460866D+00
  else if ( p == 13.0D+00 ) then
    zeta = 1.0001227133D+00
  else if ( p == 14.0D+00 ) then
    zeta = 1.0000612482D+00
  else if ( p == 15.0D+00 ) then
    zeta = 1.0000305882D+00
  else if ( p == 16.0D+00 ) then
    zeta = 1.0000152823D+00
  else if ( p == 17.0D+00 ) then
    zeta = 1.0000076372D+00
  else if ( p == 18.0D+00 ) then
    zeta = 1.0000038173D+00
  else if ( p == 19.0D+00 ) then
    zeta = 1.0000019082D+00
  else if ( p == 20.0D+00 ) then
    zeta = 1.0000009540D+00
  else

    zsum = 0.0D+00
    n = 0

    do

      n = n + 1
      zsum_old = zsum
      zsum = zsum + 1.0D+00 / ( real ( n, kind = 8 ) )**p
      if ( zsum <= zsum_old ) then
        exit
      end if

    end do

    zeta = zsum

  end if

  return
end
subroutine zipf_cdf ( x, a, cdf )

!*****************************************************************************80
!
!! ZIPF_CDF evaluates the Zipf CDF.
!
!  Discussion:
!
!    Simple summation is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the PDF.
!    1 <= N
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    1.0D+00 < A.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), save :: asave = 0.0D+00
  real ( kind = 8 ), save :: c = 0.0D+00
  real ( kind = 8 ) cdf
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
  real ( kind = 8 ) zeta

  if ( x < 1 ) then

    cdf = 0.0D+00

  else

    if ( a /= asave ) then

      c = zeta ( a )
      asave = a

    end if

    cdf = 0.0D+00
    do y = 1, x
      pdf = ( 1.0D+00 / real ( y, kind = 8 )**a ) / c
      cdf = cdf + pdf
    end do

  end if

  return
end
function zipf_check ( a )

!*****************************************************************************80
!
!! ZIPF_CHECK checks the parameter of the Zipf PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    1.0D+00 < A.
!
!    Output, logical ZIPF_CHECK, is true if the parameters are legal.
!
  implicit none

  real ( kind = 8 ) a
  logical zipf_check

  if ( a <= 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ZIPF_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 1.'
    zipf_check = .false.
    return
  end if

  zipf_check = .true.

  return
end
subroutine zipf_mean ( a, mean )

!*****************************************************************************80
!
!! ZIPF_MEAN returns the mean of the Zipf PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    1.0D+00 < A.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!    The mean is only defined for 2 < A.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) mean
  real ( kind = 8 ) zeta

  if ( a <= 2.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ZIPF_MEAN - Fatal error!'
    write ( *, '(a)' ) '  No mean defined for A <= 2.'
    stop
  end if

  mean = zeta ( a - 1.0D+00 ) / zeta ( a )

  return
end
subroutine zipf_pdf ( x, a, pdf )

!*****************************************************************************80
!
!! ZIPF_PDF evaluates the Zipf PDF.
!
!  Discussion:
!
!    PDF(A;X) = ( 1 / X**A ) / C
!
!    where the normalizing constant is chosen so that
!
!    C = Sum ( 1 <= I < Infinity ) 1 / I**A.
!
!    From observation, the frequency of different words in long
!    sequences of text seems to follow the Zipf PDF, with
!    parameter A slightly greater than 1.  The Zipf PDF is sometimes
!    known as the "discrete Pareto" PDF.
!
!    Lotka's law is a version of the Zipf PDF in which A is 2 or approximately
!    2.  Lotka's law describes the frequency of publications by authors in a
!    given field, and estimates that the number of authors with X papers is
!    about 1/X^A of the number of authors with 1 paper.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alfred Lotka,
!    The frequency distribution of scientific productivity,
!    Journal of the Washington Academy of Sciences,
!    Volume 16, Number 12, 1926, pages 317-324.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the PDF.
!    1 <= N
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    1.0D+00 < A.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), save :: asave = 0.0D+00
  real ( kind = 8 ), save :: c = 0.0D+00
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  real ( kind = 8 ) zeta

  if ( x < 1 ) then

    pdf = 0.0D+00

  else

    if ( a /= asave ) then

      c = zeta ( a )
      asave = a

    end if

    pdf = ( 1.0D+00 / real ( x, kind = 8 )**a ) / c

  end if

  return
end
subroutine zipf_sample ( a, seed, x )

!*****************************************************************************80
!
!! ZIPF_SAMPLE samples the Zipf PDF.
!
!  Discussion:
!
!    I am concerned that there may be a discrepancy in the results
!    of this routine, which do not seem to have the predicted variances.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Luc Devroye,
!    Non-Uniform Random Variate Generation,
!    Springer Verlag, 1986, pages 550-551.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    1.0D+00 < A.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i4_huge
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) test
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  integer ( kind = 4 ) x

  test = real ( i4_huge ( ) )

  b = 2.0D+00**( a - 1.0D+00 )

  do

    u = r8_uniform_01 ( seed )
    v = r8_uniform_01 ( seed )
    w = aint ( 1.0D+00 / u**( 1.0D+00 / ( a - 1.0D+00 ) ) )
!
!  Very small values of U can cause W to be very large,
!  bigger than the largest integer...
!
    if ( test < w ) then
      cycle
    end if

    t = ( ( w + 1.0D+00 ) / w )**( a - 1.0D+00 )

    if ( v * w * ( t - 1.0D+00 ) * b <= t * ( b - 1.0D+00 ) ) then
      exit
    end if

  end do

  x = int ( w )

  return
end
subroutine zipf_variance ( a, variance )

!*****************************************************************************80
!
!! ZIPF_VARIANCE returns the variance of the Zipf PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    1.0D+00 < A.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!    The variance is only defined for 3 < A.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance
  real ( kind = 8 ) zeta

  if ( a <= 3.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ZIPF_VARIANCE - Fatal error!'
    write ( *, '(a)' ) '  No variance defined for A <= 3.0.'
    stop
  end if

  call zipf_mean ( a, mean )

  variance = zeta ( a - 2.0D+00 ) / zeta ( a ) - mean * mean

  return
end
