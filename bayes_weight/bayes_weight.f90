program main

!*****************************************************************************80
!
!! MAIN is the main program for BAYES_WEIGHT.
!
!  Discussion:
!
!    BAYES_WEIGHT does a simple demonstration of Bayesian statistics.
!
!    Choose one of two dice each time and roll it.  The probabilities
!    of choosing die 1 or 2 are the unknown weights W1 and W2.  The
!    PDF's of each die are known.  Estimate the values of W1 and W2
!    by observing and analyzing the results of a sequence of rolls.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   17 December 1999
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: comp_num = 2
  integer ( kind = 4 ), parameter :: elem_num = 6

  real ( kind = 8 ) alpha_init(comp_num)
  real ( kind = 8 ) alpha_post(comp_num)
  real ( kind = 8 ) alpha_prior(comp_num)
  real ( kind = 8 ) coef(comp_num)
  integer ( kind = 4 ) comp
  integer ( kind = 4 ) comp_i
  real ( kind = 8 ) comp_weight(comp_num)
  real ( kind = 8 ) comp_weight_est(comp_num)
  integer ( kind = 4 ) count_comp(comp_num)
  integer ( kind = 4 ) count_obs(elem_num)
  integer ( kind = 4 ) count_total(elem_num)
  integer ( kind = 4 ) elem_i
  integer ( kind = 4 ) obs_i
  integer ( kind = 4 ) obs_num
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) roll_i
  integer ( kind = 4 ) rolls_per_obs
  integer ( kind = 4 ) rolls_total
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta_true(elem_num,comp_num)
!
!  Blather
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BAYES_WEIGHT'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Simple Bayesian Statistics demonstrations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Suppose an event is a roll of a die.'
  write ( *, '(a)' ) '  We have two dice to choose from, each with'
  write ( *, '(a)' ) '  a different characteristic probability'
  write ( *, '(a)' ) '  density function (PDF).  We know the form'
  write ( *, '(a)' ) '  of both PDF''s.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For each event, the probability we will '
  write ( *, '(a)' ) '  choose die #1 or die #2 is W1 or W2.'
  write ( *, '(a)' ) '  But W1 and W2 are unknown to us, and we must'
  write ( *, '(a)' ) '  deduce them from the results of the rolls.'
!
!  Initialize.
!
  comp_weight(1) = 0.75D+00
  comp_weight(2) = 0.25D+00

  theta_true(1,1) = 1.0D+00 / 16.0D+00
  theta_true(2,1) = 2.0D+00 / 16.0D+00
  theta_true(3,1) = 4.0D+00 / 16.0D+00
  theta_true(4,1) = 6.0D+00 / 16.0D+00
  theta_true(5,1) = 1.0D+00 / 16.0D+00
  theta_true(6,1) = 2.0D+00 / 16.0D+00

  theta_true(1,2) = 9.0D+00 / 16.0D+00
  theta_true(2,2) = 1.0D+00 / 16.0D+00
  theta_true(3,2) = 1.0D+00 / 16.0D+00
  theta_true(4,2) = 0.0D+00 / 16.0D+00
  theta_true(5,2) = 3.0D+00 / 16.0D+00
  theta_true(6,2) = 2.0D+00 / 16.0D+00

  seed = 123456789

  rolls_per_obs = 10
  obs_num = 1000

  rolls_total = 0
  alpha_init(1:comp_num) = 1.0D+00

  call dirichlet_mean ( comp_num, alpha_init, comp_weight_est )

  alpha_post(1:comp_num) = alpha_init(1:comp_num)
  count_comp(1:comp_num) = 0
  count_total(1:elem_num) = 0
!
!  Report run parameters:
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Run parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of components is ', comp_num

  call r8vec_print ( comp_num, comp_weight, '  Exact weights:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, THETA_TRUE(I,1)  THETA_TRUE(I,2):'
  write ( *, '(a)' ) ' '
  do elem_i = 1, elem_num
    write ( *, '(2x,i6,2g14.6)' ) elem_i, theta_true(elem_i,1), &
      theta_true(elem_i,2)
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
  write ( *, '(a)' ) '  I, ALPHA(I), WEIGHT_EST(I):'
  write ( *, '(a)' ) ' '
  do comp_i = 1, comp_num
    write ( *, '(2x,i6,2g14.6)' ) comp_i, alpha_init(comp_i), &
      comp_weight_est(comp_i)
  end do
!
!  Observe the system.
!
  do obs_i = 1, obs_num

    alpha_prior(1:comp_num) = alpha_post(1:comp_num)

    count_obs(1:elem_num) = 0
!
!  Roll the dice a bunch of times.
!
    do roll_i = 1, rolls_per_obs

      rolls_total = rolls_total + 1

      call discrete_sample ( comp_num, comp_weight, seed, comp )

      count_comp(comp) = count_comp(comp) + 1

      call discrete_sample ( elem_num, theta_true(1,comp), elem_i )

      count_obs(elem_i) = count_obs(elem_i) + 1
!
!  Here, I analyze each roll individually.  Do I have to do this?
!  Can I do it at the end of the loop, using only the summed
!  results?
!
      do comp_i = 1, comp_num
        call discrete_pdf ( elem_i, elem_num, theta_true(1,comp_i), pdf )
        coef(comp_i) = pdf * alpha_prior(comp_i)
      end do

      coef(1:comp_num) = coef(1:comp_num) / sum ( coef(1:comp_num) )

      alpha_post(1:comp_num) = alpha_prior(1:comp_num) + coef(1:comp_num)

    end do

    count_total(1:elem_num) = count_total(1:elem_num) + count_obs(1:elem_num)

    call dirichlet_mean ( comp_num, alpha_post, comp_weight_est )
!
!  Print out the data.
!
    if ( obs_i <= 10 .or. mod ( obs_i, obs_num / 10 ) == 0 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'BAYES_BETA - Observation: ', obs_i
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Current and summary events'
      write ( *, '(a)' ) ' '
      do elem_i = 1, elem_num
        write ( *, '(2x,i6,2x,2i8)' ) elem_i, count_obs(elem_i), &
          count_total(elem_i)
      end do
      write ( *, '(2x,a6,2x,2i8)' ) 'Total', rolls_per_obs, rolls_total
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  (Unknown) component calls and probability:'
      write ( *, '(a)' ) ' '
      do comp_i = 1, comp_num
        write ( *, '(2x,i6,2x,i8,g14.6)' ) comp_i, count_comp(comp_i), &
          real ( count_comp(comp_i), kind = 8 ) / real ( rolls_total, kind = 8 )
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimated parameters:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I, ALPHA(I), WEIGHT_EST(I):'
      write ( *, '(a)') ' '
      do comp_i = 1, comp_num
        write ( *, '(2x,i6,2g14.6)' ) comp_i, alpha_post(comp_i), &
          comp_weight_est(comp_i)
      end do

    end if

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BAYES_WEIGHT:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
  real ( kind = 8 ) mean(n)

  mean(1:n) = a(1:n) / sum ( a(1:n) )

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
!    0.0 <= CDF <= 1.0.
!
!    Input, integer ( kind = 4 ) A, the number of probabilities assigned.
!
!    Input, real ( kind = 8 ) B(A), the relative probabilities of outcomes 1 through A.
!    Each entry must be nonnegative.
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
subroutine discrete_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! DISCRETE_PDF evaluates the Discrete PDF.
!
!  Formula:
!
!    PDF(X)(A,B) = B(X) if 1 <= X <= A
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
!    Input, real ( kind = 8 ) B(A), the relative probabilities of outcomes 1 through A.
!    Each entry must be nonnegative.
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
!    Input, real ( kind = 8 ) B(A), the relative probabilities of outcomes 1 through A.
!    Each entry must be nonnegative.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
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
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC, with an optional title.
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
    write ( *, '(a)' ) title
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
  end do

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
