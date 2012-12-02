program main

!*****************************************************************************80
!
!! MAIN is the main program for BDMLIB_PRB.
!
!  Discussion:
!
!    BDMLIB_PRB tests the BDMLIB Bayes Dirichlet Mixture routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: acid_num = 20
  integer ( kind = 4 ), parameter :: comp_max = 9
  integer ( kind = 4 ), parameter :: iunit = 1

  character acid_sym(acid_num)
  real ( kind = 8 ) beta(acid_num,comp_max)
  real ( kind = 8 ) beta_sum(comp_max)
  integer ( kind = 4 ) comp_label(comp_max)
  integer ( kind = 4 ) comp_num
  real ( kind = 8 ) comp_weight(comp_max)
  integer ( kind = 4 ) ierror
  character ( len = 30 ) mixture_file_name

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BDMLIB_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BDMLIB library.'
!
!  Read information about the mixture.
!
  mixture_file_name = 'mixture.txt'

  open ( unit = iunit, file = mixture_file_name, form = 'formatted' )

  call mixture_read ( acid_num, acid_sym, beta, beta_sum, &
    comp_label, comp_max, comp_num, comp_weight, ierror, iunit )

  close ( unit = iunit )
!
!  Print the amino acid parameters.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Numeric key for amino acid abbreviations.'
  write ( *, '(a)' ) ' '

  call amino_print ( acid_num, acid_sym )
!
!  Print the component parameters.
!
  call comp_param_print ( acid_num, acid_sym, comp_max, comp_num, &
    beta, beta_sum, comp_weight )
!
!  Check the parameters.
!
  call dirichlet_mix_check ( comp_num, acid_num, acid_num, beta, comp_weight )
!
!  Perform the simple test of generating 10 isoleucines in a row.
!
  call test01 ( acid_num, beta, comp_num, comp_weight )
!
!  Now test a random sampling.
!
  call test02 ( acid_num, beta, comp_num, comp_weight )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BDMLIB_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( acid_num, beta, comp_num, comp_weight )

!*****************************************************************************80
!
!! TEST01 generates 10 isoleucine events in a row.
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
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_num

  integer ( kind = 4 ) acid_i
  real ( kind = 8 ) alpha(comp_num)
  real ( kind = 8 ) alpha_0(comp_num)
  real ( kind = 8 ) beta(acid_num,comp_num)
  integer ( kind = 4 ) comp_i
  real ( kind = 8 ) comp_weight(comp_num)
  real ( kind = 8 ) comp_weight_est(comp_num)
  integer ( kind = 4 ) event_i
  integer ( kind = 4 ) event_num
  real ( kind = 8 ) p(comp_num)
  real ( kind = 8 ) p_hat(comp_num)
  integer ( kind = 4 ) site_num
  integer ( kind = 4 ) x_sample(acid_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Generate a (nonrandom) sequence of'
  write ( *, '(a)' ) '  10 isoleucine results in a row.'
  write ( *, '(a)' ) ' '
!
!  Initialize information about the Bayesian process.
!
  alpha_0(1:comp_num) = 1.0D+00

  alpha(1:comp_num) = alpha_0(1:comp_num)

  p(1:comp_num) = 1.0D+00 / real ( comp_num, kind = 8 )
  p_hat(1:comp_num) = 1.0D+00 / real ( comp_num, kind = 8 )

  event_num = 10

  call r8vec_print ( comp_num, comp_weight, 'Exact component weights:' )

  call r8vec_print ( comp_num, alpha, 'Initial ALPHA:' )
!
!  Based on the current ALPHA's, compute the mean/expected value/estimate 
!  for the weights.
!
  call dirichlet_mean ( comp_num, alpha, comp_weight_est )

  call r8vec_print ( comp_num, comp_weight_est, &
    'Initial estimated component weights:' )

  site_num = 1

  do event_i = 1, event_num
!
!  Observe a single isoleucine.
!
    x_sample(1:acid_num) = 0
    x_sample(8) = site_num
!
!  Update ALPHA, the estimated weight parameters, based on X.
!
    call event_process ( acid_num, alpha, beta, comp_num, p, p_hat, site_num, &
      x_sample )

    call r8vec_print ( comp_num, alpha, 'Current ALPHA:' )
!
!  Based on the current ALPHA's, compute the mean/expected value/estimate 
!  for the weights.
!
    call dirichlet_mean ( comp_num, alpha, comp_weight_est )

    call r8vec_print ( comp_num, comp_weight_est, &
      'Estimated component weights:' )

  end do

  return
end
subroutine test02 ( acid_num, beta, comp_num, comp_weight )

!*****************************************************************************80
!
!! TEST02 generates random events.
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
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_num

  integer ( kind = 4 ) acid_i
  real ( kind = 8 ) alpha(comp_num)
  real ( kind = 8 ) alpha_0(comp_num)
  real ( kind = 8 ) beta(acid_num,comp_num)
  integer ( kind = 4 ) comp_i
  integer ( kind = 4 ) comp_sample
  real ( kind = 8 ) comp_weight(comp_num)
  real ( kind = 8 ) comp_weight_est(comp_num)
  integer ( kind = 4 ) event_i
  integer ( kind = 4 ) event_num
  real ( kind = 8 ) p(comp_num)
  real ( kind = 8 ) p_hat(comp_num)
  real ( kind = 8 ) p_sample(acid_num)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) site_num
  integer ( kind = 4 ) x_sample(acid_num)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Generate many random events.'
  write ( *, '(a)' ) '  We should be able to approximate the'
  write ( *, '(a)' ) '  exact component weights.'
  write ( *, '(a)' ) ' '
!
!  Initialize information about the Bayesian process.
!
  alpha_0(1:comp_num) = 1.0D+00

  alpha(1:comp_num) = alpha_0(1:comp_num)

  p(1:comp_num) = 1.0D+00 / real ( comp_num, kind = 8 )
  p_hat(1:comp_num) = 1.0D+00 / real ( comp_num, kind = 8 )

  event_num = 1000

  call r8vec_print ( comp_num, comp_weight, 'Exact component weights:' )

  call r8vec_print ( comp_num, alpha, 'Initial ALPHA:' )
!
!  Based on the current ALPHA's, compute the mean/expected value/estimate 
!  for the weights.
!
  call dirichlet_mean ( comp_num, alpha, comp_weight_est )

  call r8vec_print ( comp_num, comp_weight_est, &
    'Initial estimated component weights:' )

  site_num = 10

  do event_i = 1, event_num
!
!  Randomly choose COMP_SAMPLE, the component PDF to sample.
!
    call discrete_sample ( comp_num, comp_weight, seed, comp_sample )
!
!  Now generate the probabilities P_SAMPLE for a multinomial PDF by sampling
!  the COMP_SAMPLE-th Dirichlet distribution.
!
    call dirichlet_sample ( acid_num, beta(1,comp_sample), seed, p_sample )
!
!  Now generate the event X_SAMPLE by sampling the multinomial PDF with
!  the given P_SAMPLE parameters.
!
    call multinomial_sample ( site_num, acid_num, p_sample, seed, x_sample )
!
!  Update ALPHA, the estimated weight parameters, based on X.
!
    call event_process ( acid_num, alpha, beta, comp_num, p, p_hat, site_num, &
      x_sample )
!
!  Based on the current ALPHA's, compute the mean/expected value/estimate 
!  for the weights.
!
    call dirichlet_mean ( comp_num, alpha, comp_weight_est )

    if ( event_i <= 10 .or. mod ( event_i, event_num/10 ) == 0 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'Event ', event_i

      call r8vec_print ( comp_num, alpha, 'Current ALPHA:' )

      call r8vec_print ( comp_num, comp_weight_est, &
        'Estimated component weights:' )

    end if

  end do

  return
end
