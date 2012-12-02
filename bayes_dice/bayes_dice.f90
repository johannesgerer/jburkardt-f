program main

!*****************************************************************************80
!
!! MAIN is the main program for BAYES_DICE.
!
!  Discussion:
!
!    BAYES_DICE does a simple demonstration of Bayesian statistics.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   15 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: elem_num = 6

  integer ( kind = 4 ) obs_num
  real ( kind = 8 ) theta_true(elem_num)
  integer ( kind = 4 ) rolls_per_obs

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BAYES_DICE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Simple Bayesian Statistics demonstrations.'

  theta_true(1:6) = (/ 1.0D+00, 2.0D+00, 4.0D+00, 6.0D+00, 1.0D+00, 2.0D+00 /) &
    / 16.0D+00

  rolls_per_obs = 500
  obs_num = 4

  call test01 ( elem_num, obs_num, rolls_per_obs, theta_true )

  theta_true(1:6) = (/ 1.0D+00, 2.0D+00, 4.0D+00, 6.0D+00, 1.0D+00, 2.0D+00 /) &
    / 16.0D+00

  rolls_per_obs = 5000
  obs_num = 4

  call test01 ( elem_num, obs_num, rolls_per_obs, theta_true )

  theta_true(1:6) = (/ 1.0D+00, 2.0D+00, 4.0D+00, 6.0D+00, 1.0D+00, 2.0D+00 /) &
    / 16.0D+00

  rolls_per_obs = 500
  obs_num = 40

  call test01 ( elem_num, obs_num, rolls_per_obs, theta_true )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BAYES_DICE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( elem_num, obs_num, rolls_per_obs, theta_true )

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
!   30 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEM_NUM, the number of sides on the dice.
!
!    Input, integer OBS_NUM, the number of observations.
!
!    Input, integer ROLLS_PER_OBS, the number of dice rolls in each observation.
!
!    Input, real ( kind = 8 ) THETA_TRUE(ELEM_NUM), the true weights
!    of the dice outcomes.
!
  implicit none

  integer ( kind = 4 ) elem_num

  real ( kind = 8 ) alpha_init(elem_num)
  real ( kind = 8 ) alpha_post(elem_num)
  real ( kind = 8 ) alpha_prior(elem_num)
  integer ( kind = 4 ) count_obs(elem_num)
  integer ( kind = 4 ) count_total(elem_num)
  integer ( kind = 4 ) elem_i
  integer ( kind = 4 ) obs_i
  integer ( kind = 4 ) obs_num
  integer ( kind = 4 ) roll_i
  integer ( kind = 4 ) rolls_per_obs
  integer ( kind = 4 ) rolls_total
  real ( kind = 8 ) theta_init(elem_num)
  real ( kind = 8 ) theta_post(elem_num)
  real ( kind = 8 ) theta_true(elem_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BAYES_DICE'
  write ( *, '(a)' ) '  Simple Bayesian Statistics demonstrations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Suppose we''re watching a "system" and trying'
  write ( *, '(a)' ) '  to analyze its behavior.  Each time we observe'
  write ( *, '(a)' ) '  the system, it rolls a die a certain number of'
  write ( *, '(a)' ) '  times, and reports the results.  We want to'
  write ( *, '(a)' ) '  estimate THETA(1) through and THETA(6),'
  write ( *, '(a)' ) '  the probabilities of each result.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We treat the values of THETA as'
  write ( *, '(a)' ) '  random variables themselves, controlled by a'
  write ( *, '(a)' ) '  Dirichlet probability density function with'
  write ( *, '(a)' ) '  parameters ALPHA(1) through ALPHA(6).  We make'
  write ( *, '(a)' ) '  an arbitrary guess for initial ALPHA values.'
  write ( *, '(a)' ) '  We observe the system and adjust the ALPHA''s'
  write ( *, '(a)' ) '  using Bayes''s Law.'
!
!  Get the initial estimates of the parameters.
!
  rolls_total = 0
  alpha_init(1:elem_num) = 1.0D+00

  call dirichlet_mean ( elem_num, alpha_init, theta_init )

  alpha_post(1:elem_num) = alpha_init(1:elem_num)
  theta_post(1:elem_num) = theta_init(1:elem_num)
  count_total(1:elem_num) = 0
!
!  Report run parameters:
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Run parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, THETA_TRUE(I):'
  write ( *, '(a)' ) ' '
  do elem_i = 1, elem_num
    write ( *, '(2x,i6,g14.6)' ) elem_i, theta_true(elem_i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of observations =   ', obs_num
  write ( *, '(a,i6)' ) '  Rolls per observation =    ', rolls_per_obs
!
!  Report initial parameters:
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Initial parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, ALPHA_INIT(I), THETA_INIT(I):'
  write ( *, '(a)' ) ' '
  do elem_i = 1, elem_num
    write ( *, '(2x,i6,2g14.6)' ) elem_i, alpha_init(elem_i), theta_init(elem_i)
  end do
!
!  Make an observation of the system.
!
  do obs_i = 1, obs_num

    alpha_prior(1:elem_num) = alpha_post(1:elem_num)
    count_obs(1:elem_num) = 0
!
!  Roll the dice.
!
    do roll_i = 1, rolls_per_obs

      rolls_total = rolls_total + 1

      call discrete_sample ( elem_num, theta_true, elem_i )

      count_obs(elem_i) = count_obs(elem_i) + 1

    end do
!
!  Use the observations to adjust our estimates of the system.
!
    alpha_post(1:elem_num) = alpha_prior(1:elem_num) + count_obs(1:elem_num)
    count_total(1:elem_num) = count_total(1:elem_num) + count_obs(1:elem_num)

    call dirichlet_mean ( elem_num, alpha_post, theta_post )
!
!  Print out the data.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'BAYES_DICE - Observation: ', obs_i
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Observed data:'
    write ( *, '(a)' ) ' '
    do elem_i = 1, elem_num
      write ( *, '(2x,i6,2x,i6)' ) elem_i, count_obs(elem_i)
    end do
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Post observation data:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Total rolls:      ', rolls_total
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Estimated parameters:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  I, COUNT(I), ALPHA_POST(I), THETA_POST(I):'
    write ( *, '(a)' ) ' '
    do elem_i = 1, elem_num
      write ( *, '(2x,i6,i6,2g14.6)' ) elem_i, count_total(elem_i), &
        alpha_post(elem_i), theta_post(elem_i)
    end do

  end do

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
!    Each A(I) should be nonnegative, and at least one should be positive.
!
!    Output, real ( kind = 8 ) MEAN(N), the means of the PDF.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_sum
  integer ( kind = 4 ) i
  real ( kind = 8 ) mean(n)
!
!  Check.
!
  do i = 1, n
    if ( a(i) < 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIRICHLET_MEAN - Fatal error!'
      write ( *, '(a)' ) '  A(I) < 0.'
      stop
    end if
  end do

  a_sum = sum ( a(1:n) )

  mean(1:n) = a(1:n) / a_sum

  return
end
subroutine discrete_cdf_inv ( cdf, n, p, i )

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
!    0.0 <= CDF <= 1.0.
!
!    Input, integer ( kind = 4 ) N, the number of probabilities assigned.
!
!    Input, real ( kind = 8 ) P(N), the relative probabilities of
!    outcomes 1 through N.  Each entry must be nonnegative.
!
!    Output, integer ( kind = 4 ) I, the corresponding argument for which
!    CDF(I-1) < CDF <= CDF(I)
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cdf
  real ( kind = 8 ) cum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) p_sum
!
!  Check.
!
  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DISCRETE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  do j = 1, n
    if ( p(j) < 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISCRETE_CDF_INV - Fatal error!'
      write ( *, '(a)' ) '  Negative probabilities not allowed.'
      stop
    end if
  end do

  p_sum = sum ( p(1:n) )

  if ( p_sum == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DISCRETE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  Total probablity is zero.'
    stop
  end if

  cum = 0.0D+00

  do j = 1, n

    cum = cum + p(j) / p_sum

    if ( cdf <= cum ) then
      i = j
      return
    end if

  end do

  i = n

  return
end
subroutine discrete_sample ( n, p, i )

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
!    Input, integer ( kind = 4 ) N, the number of probabilities assigned.
!
!    Input, real ( kind = 8 ) P(N), the relative probabilities of
!    outcomes 1 through N.  Each entry must be nonnegative.
!
!    Output, integer ( kind = 4 ) I, a sample of the PDF.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cdf
  real ( kind = 8 ) cdf_max
  real ( kind = 8 ) cdf_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) p_sum
  real ( kind = 8 ) uniform_01_sample
!
!  Check.
!
  do j = 1, n
    if ( p(j) < 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISCRETE_SAMPLE - Fatal error!'
      write ( *, '(a)' ) '  Negative probabilities not allowed.'
      stop
    end if
  end do

  p_sum = sum ( p(1:n) )

  if ( p_sum == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DISCRETE_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Total probablity is zero.'
    stop
  end if

  cdf_min = 0.0D+00
  cdf_max = 1.0D+00

  call r8_random ( cdf_min, cdf_max, cdf )

  call discrete_cdf_inv ( cdf, n, p, i )

  return
end
subroutine r8_random ( rlo, rhi, r )

!*****************************************************************************80
!
!! R8_RANDOM returns a random R8 in a given range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RLO, RHI, the minimum and maximum values.
!
!    Output, real ( kind = 8 ) R, the randomly chosen value.
!
  implicit none

  logical, save :: seed = .false.
  real ( kind = 8 ) r
  real ( kind = 8 ) rhi
  real ( kind = 8 ) rlo
  real ( kind = 8 ) t

  if ( .not. seed ) then
    call random_seed ( )
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R.
!
  r = ( 1.0D+00 - t ) * rlo + t * rhi

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
