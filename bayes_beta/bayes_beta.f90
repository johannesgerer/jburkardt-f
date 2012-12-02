program main

!*****************************************************************************80
!
!! MAIN is the main program for BAYES_BETA.
!
!  Discussion:
!
!    BAYES_BETA does a simple demonstration of Bayesian statistics.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   06 December 1999
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ntoss_obs
  integer ( kind = 4 ) obs_num
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) theta1_true

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BAYES_BETA'
  write ( *, '(a)' ) '  Simple Bayesian Statistics demonstrations.'

  theta1_true = 0.75D+00
  ntoss_obs = 5
  obs_num = 4
  seed = 123456789

  call test01 ( ntoss_obs, obs_num, theta1_true, seed )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BAYES_BETA:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( ntoss_obs, obs_num, theta1_true, seed )

!*****************************************************************************80
!
!! TEST01 does a simple demonstration of Bayesian statistics.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   18 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTOSS_OBS, the number of coin tosses
!    per observation.
!
!    Input, integer ( kind = 4 ) OBS_NUM, the number of observations.
!
!    Input, real ( kind = 8 ) THETA1_TRUE, the value of THETA1, which the
!    program will try to estimate.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
  implicit none

  real    ( kind = 8 ) alpha1_init
  real    ( kind = 8 ) alpha1_post
  real    ( kind = 8 ) alpha1_prior
  real    ( kind = 8 ) alpha2_init
  real    ( kind = 8 ) alpha2_post
  real    ( kind = 8 ) alpha2_prior
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nhead_obs
  integer ( kind = 4 ) nhead_post
  integer ( kind = 4 ) nhead_prior
  integer ( kind = 4 ) ntail_obs
  integer ( kind = 4 ) ntail_post
  integer ( kind = 4 ) ntail_prior
  integer ( kind = 4 ) ntoss_obs
  integer ( kind = 4 ) ntoss_post
  integer ( kind = 4 ) ntoss_prior
  integer ( kind = 4 ) obs_i
  integer ( kind = 4 ) obs_num
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) step_num
  real    ( kind = 8 ) theta1_init
  real    ( kind = 8 ) theta1_post
  real    ( kind = 8 ) theta1_prior
  real    ( kind = 8 ) theta1_true
  real    ( kind = 8 ) theta2_init
  real    ( kind = 8 ) theta2_post
  real    ( kind = 8 ) theta2_prior
  real    ( kind = 8 ) theta2_true
  real    ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BAYES_BETA'
  write ( *, '(a)' ) '  Simple Bayesian Statistics demonstrations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Suppose we''re watching a "system" and trying'
  write ( *, '(a)' ) '  to analyze its behavior.  Each time we observe'
  write ( *, '(a)' ) '  the system, it flips a coin a certain number of'
  write ( *, '(a)' ) '  times, and reports the number of heads and'
  write ( *, '(a)' ) '  tails.  We want to estimate THETA1 and THETA2,'
  write ( *, '(a)' ) '  the probabilities of heads and of tails.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We treat the values of THETA1 and THETA2 as'
  write ( *, '(a)' ) '  random variables themselves, controlled by a'
  write ( *, '(a)' ) '  Beta probability density function, which has'
  write ( *, '(a)' ) '  parameters ALPHA1 and ALPHA2.  We make an '
  write ( *, '(a)' ) '  arbitrary or informed guess for initial values'
  write ( *, '(a)' ) '  of ALPHA1 and ALPHA2.  We observe the system,'
  write ( *, '(a)' ) '  and adjust ALPHA1 and ALPHA2 using Bayes''s'
  write ( *, '(a)' ) '  Law.  We continue until we are satisfied'
  write ( *, '(a)' ) '  that our estimates seem to have converged.'
!
!  Get the initial estimates of the parameters ALPHA1 and ALPHA2.
!
  alpha1_init = 0.5D+00
  alpha2_init = 0.5D+00

  call beta_max ( alpha1_init, alpha2_init, theta1_init, theta2_init )

  alpha1_post = alpha1_init
  alpha2_post = alpha2_init
  nhead_post = 0
  ntail_post = 0
  ntoss_post = 0
  theta1_post = theta1_init
  theta2_post = theta2_init
  theta1_true = max ( theta1_true, 0.0D+00 )
  theta1_true = min ( theta1_true, 1.0D+00 )
  theta2_true = 1.0D+00 - theta1_true
!
!  Report on run parameters:
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Run parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  THETA1/THETA2:    ', theta1_true, theta2_true
  write ( *, '(a,i6)' ) '  Number of observations to make =     ', obs_num
  write ( *, '(a,i6)' ) '  Number of tosses per observation =   ', ntoss_obs
!
!  Report on initial data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initial data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  ALPHA1/ALPHA2     ', alpha1_init, alpha2_init
  write ( *, '(a,2g14.6)' ) '  THETA1/THETA2:    ', theta1_init, theta2_init

  step_num = 20
  call beta_plot ( step_num, alpha1_init, alpha2_init )
!
!  Make observations of the system.
!
  do obs_i = 1, obs_num

    alpha1_prior = alpha1_post
    alpha2_prior = alpha2_post
    nhead_prior = nhead_post
    ntail_prior = ntail_post
    ntoss_prior = ntoss_post
    theta1_prior = theta1_post
    theta2_prior = theta2_post

    nhead_obs = 0
    ntail_obs = 0
!
!  Flip the coin NTOSS_OBS times.
!
    do i = 1, ntoss_obs

      x = r8_uniform_01 ( seed )

      if ( x <= theta1_true ) then
        nhead_obs = nhead_obs + 1
      else
        ntail_obs = ntail_obs + 1
    end if

  end do
!
!  Use the observations to adjust our estimates of the system.
!
    alpha1_post = alpha1_prior + nhead_obs
    alpha2_post = alpha2_prior + ntail_obs

    nhead_post = nhead_prior + nhead_obs
    ntail_post = ntail_prior + ntail_obs
    ntoss_post = ntoss_prior + ntoss_obs

    call beta_max ( alpha1_post, alpha2_post, theta1_post, theta2_post )
!
!  Print out the data.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'BAYES_BETA - After observation ', obs_i
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Prior observation data:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Tosses:           ', ntoss_prior
    write ( *, '(a,2i6)' ) '  Heads/Tails:      ', nhead_prior, ntail_prior
    write ( *, '(a,2g14.6)' ) '  ALPHA1/ALPHA2:    ', alpha1_prior, &
      alpha2_prior
    write ( *, '(a,2g14.6)' ) '  THETA1/THETA2:    ', theta1_prior, &
      theta2_prior
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Observation data:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Tosses:           ', ntoss_obs
    write ( *, '(a,2i6)' ) '  Heads/Tails:      ', nhead_obs, ntail_obs
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Post observation data:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Tosses:           ', ntoss_post
    write ( *, '(a,2i6)' ) '  Heads/Tails:      ', nhead_post, ntail_post
    write ( *, '(a,2g14.6)' ) '  ALPHA1/ALPHA2:    ', alpha1_post, alpha2_post
    write ( *, '(a,2g14.6)' ) '  THETA1/THETA2:    ', theta1_post, theta2_post

  end do

  step_num = 20
  call beta_plot ( step_num, alpha1_post, alpha2_post )

  return
end
function beta ( x, y )

!*****************************************************************************80
!
!! BETA returns the value of the Beta function.
!
!  Discussion:
!
!    BETA(X,Y) = ( GAMMA(X) * GAMMA(Y) ) / GAMMA(X+Y)
!    BETA(X,Y) = BETA(Y,X).
!    BETA(X,Y) = Integral ( 0 <= T <= 1 ) T**(X-1) (1-T)**(Y-1) dT.
!
!    Both X and Y must be greater than 0.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the two parameters that define the
!    Beta function.
!
!    Output, real ( kind = 8 ) BETA, the value of the Beta function.
!
  implicit none

  real    ( kind = 8 ) beta
  real    ( kind = 8 ) gamma_log
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y

  if ( x <= 0.0D+00 .or. y <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA - Fatal error!'
    write ( *, '(a)' ) '  Both X and Y must be greater than 0.'
    stop
  end if

  beta = exp ( gamma_log ( x ) + gamma_log ( y ) - gamma_log ( x + y ) )

  return
end
subroutine beta_max ( alpha1, alpha2, t1, t2 )

!*****************************************************************************80
!
!! BETA_MAX returns the most likely values of T1 and T2 for a Beta PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA1, ALPHA2, the parameter values for the Beta
!    probability density function.
!
!    Output, real ( kind = 8 ) T1, T2, the most likely values of T1 and T2,
!    which are controlled by the Beta probability density function.
!
  implicit none

  real    ( kind = 8 ) alpha1
  real    ( kind = 8 ) alpha2
  real    ( kind = 8 ) t1
  real    ( kind = 8 ) t2

  t1 = alpha1 / ( alpha1 + alpha2 )
  t2 = alpha2 / ( alpha1 + alpha2 )

  return
end
function beta_pdf ( t, x, y )

!*****************************************************************************80
!
!! BETA_PDF returns the value of the Beta probability density function.
!
!  Discussion:
!
!    BETA_PDF(T,X,Y) = T**(X-1) * (1-T)**(Y-1) / BETA(X,Y).
!
!    BETA_PDF(T,X,Y) = BETA_PDF(T,Y,X).
!    BETA_PDF(0,X,Y) = 0.
!    BETA_PDF(1,X,Y) = 0.
!    Integral ( 0 <= T <= 1 ) BETA_PDF(T,X,Y) dT = 1.
!    Expected value of T is X / ( X + Y ).
!
!    T must be in [0,1], X and Y must be greater than 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   16 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the value at which the Beta probability density
!    function is to be evaluated.
!
!    Input, real ( kind = 8 ) X, Y, the two parameters that define
!    the Beta function.  X and Y must be greater than 0.
!
!    Output, real ( kind = 8 ) BETA_PDF, the value of the Beta function.
!
  implicit none

  real    ( kind = 8 ) beta_log
  real    ( kind = 8 ) beta_pdf
  real    ( kind = 8 ) gamma_log
  real    ( kind = 8 ) t
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y

  if ( x <= 0.0D+00 .or. y <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_PDF - Fatal error!'
    write ( *, '(a)' ) '  Both X and Y must be greater than 0.'
    stop
  end if

  if ( t < 0.0D+00 .or. 1.0D+00 < t ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_PDF - Fatal error!'
    write ( *, '(a)' ) '  0 <= T <= 1 is required.'
    write ( *, '(a,g14.6)' ) '  Your value of T is ', t
    stop
  end if
!
!  If you use the straightforward formula for the Beta PDF, you
!  will get intermediate overflow with values of X and Y as small
!  as 200, although the values of BETA_PDF are reasonable.
!
  if ( t == 0.0D+00 .or. t == 1.0D+00 ) then
    beta_pdf = 0.0D+00
  else
    beta_log = ( x - 1.0D+00 ) * log ( t ) &
      + ( y - 1.0D+00 ) * log ( 1.0D+00 - t ) &
      + gamma_log ( x + y ) - gamma_log ( x ) - gamma_log ( y )

    beta_pdf = exp ( beta_log )

  end if

  return
end
subroutine beta_plot ( step_num, alpha1, alpha2 )

!*****************************************************************************80
!
!! BETA_PLOT "plots" the Beta distribution for given parameters.
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
!    Input, integer ( kind = 4 ) STEP_NUM, the number of steps to take
!    from 0 to 1 in the plotting.
!
!    Input, real ( kind = 8 ) ALPHA1, ALPHA2, the parameters of
!    the distribution.
!
  implicit none

  real    ( kind = 8 ) alpha1
  real    ( kind = 8 ) alpha2
  real    ( kind = 8 ) beta_pdf
  real    ( kind = 8 ) pdf
  integer ( kind = 4 ) step_i
  integer ( kind = 4 ) step_num
  real    ( kind = 8 ) t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Plot of Beta distribution,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  ALPHA1 = ', alpha1
  write ( *, '(a,g14.6)' ) '  ALPHA2 = ', alpha2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            BETA(T,ALPHA1,ALPHA2)'
  write ( *, '(a)' ) ' '

  do step_i = 0, step_num
    t = real ( step_i, kind = 8 ) / real ( step_num, kind = 8 )
    pdf = beta_pdf ( t, alpha1, alpha2 )
    write ( *, '(2g14.6)' ) t, pdf
  end do

  return
end
function gamma_log ( x )

!*****************************************************************************80
!
!! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!  Discussion:
!
!    The program uses rational functions that theoretically approximate
!    log ( GAMMA(X) ) to at least 18 significant decimal digits.
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
!    Hart, Ward Cheney, Charles Lawson, Maehly, Charles Mesztenyi,
!    John Rice, Thacher, Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function.
!    X must be positive.
!
!    Output, real ( kind = 8 ) GAMMA_LOG, the logarithm of the Gamma
!    function of X.  If X <= 0.0, or if overflow would occur, the
!    program returns the value HUGE().
!
!  Machine-dependent constants:
!
!    BETA   - radix for the floating-point representation.
!
!    MAXEXP - the smallest positive power of BETA that overflows.
!
!    XBIG   - largest argument for which LN(GAMMA(X)) is representable
!             in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!    XINF   - largest machine representable floating-point number;
!             approximately BETA**MAXEXP.
!
!    FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!    Approximate values for some important machines are:
!
!                              BETA      MAXEXP         XBIG
!
!    CRAY-1        (S.P.)        2        8191       9.62D+2461
!    Cyber 180/855
!      under NOS   (S.P.)        2        1070       1.72D+319
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)        2         128       4.08D+36
!    IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!    IBM 3033      (D.P.)       16          63       4.29D+73
!    VAX D-Format  (D.P.)        2         127       2.05D+36
!    VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                            FRTBIG
!
!    CRAY-1        (S.P.)   3.13D+615
!    Cyber 180/855
!      under NOS   (S.P.)   6.44D+79
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)   1.42D+9
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)   2.25D+76
!    IBM 3033      (D.P.)   2.56D+18
!    VAX D-Format  (D.P.)   1.20D+9
!    VAX G-Format  (D.P.)   1.89D+76
!
  implicit none

  real    ( kind = 8 ), parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real    ( kind = 8 ) corr
  real    ( kind = 8 ), parameter :: d1 = - 5.772156649015328605195174D-01
  real    ( kind = 8 ), parameter :: d2 =   4.227843350984671393993777D-01
  real    ( kind = 8 ), parameter :: d4 =   1.791759469228055000094023D+00
  real    ( kind = 8 ) eps
  real    ( kind = 8 ), parameter :: frtbig = 1.42D+09
  integer ( kind = 4 ) i
  real    ( kind = 8 ) gamma_log
  real    ( kind = 8 ), parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real    ( kind = 8 ), parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real    ( kind = 8 ), parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real    ( kind = 8 ), parameter :: pnt68 = 0.6796875D+00
  real    ( kind = 8 ), parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real    ( kind = 8 ), parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real    ( kind = 8 ), parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real    ( kind = 8 ) res
  real    ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real    ( kind = 8 ) x
  real    ( kind = 8 ), parameter :: xbig = 4.08D+36
  real    ( kind = 8 ) xden
  real    ( kind = 8 ) xm1
  real    ( kind = 8 ) xm2
  real    ( kind = 8 ) xm4
  real    ( kind = 8 ) xnum
  real    ( kind = 8 ) xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0D+00 .or. xbig < x ) then
    gamma_log = huge ( gamma_log )
    return
  end if

  eps = epsilon ( eps )

  if ( x <= eps ) then

    res = - log ( x )

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
    xden = -1.0D+00
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
!      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
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
!    Volume 8, Number 2, 1969, pages 136-143.
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
  real    ( kind = 8 ) r8_uniform_01
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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
