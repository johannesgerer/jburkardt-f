program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA189_PRB.
!
!  Discussion:
!
!    ASA189_PRB controls the test of ASA189.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a_est
  real ( kind = 8 ) b
  real ( kind = 8 ) b_est
  integer ( kind = 4 ) c
  real ( kind = 8 ) mu
  real ( kind = 8 ) mu_est
  integer ( kind = 4 ) sample_log
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta_est

  call timestamp ( )

  a = 2.0D+00
  b = 3.0D+00
  c = 4

  call beta_binomial_check ( a, b, c )

  mu = a / ( a + b )
  theta = 1.0D+00 / ( a + b )

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA189_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA189 library.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstration of "ASA189",'
  write ( *, '(a)' ) '  Applied Statistics Algorithm 189'
  write ( *, '(a)' ) '  Estimate the parameters A and B of a '
  write ( *, '(a)' ) '  beta binomial distribution.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          A             B                MU         THETA'
  write ( *, '(a)' ) ' '
  write ( *, '(6x,4g14.6)' ) a, b, mu, theta
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimated values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Samples   A_est         B_est        MU_est     THETA_est'
  write ( *, '(a)' ) ' '

  do sample_log = 2, 5

    sample_num = 10**sample_log

    call analyze ( sample_num, a, b, c, seed, mu_est, theta_est )
!
!  Convert the ASA189 "THETA" and "MU" parameters to "A" and "B".
!
    a_est = mu_est / theta_est
    b_est = ( 1.0D+00 - mu_est ) / theta_est

    write ( *, '(i6,4g14.6)' ) sample_num, a_est, b_est, mu_est, theta_est

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA189_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine analyze ( sample_num, a, b, c, seed, mu_est, theta_est )

!*****************************************************************************80
!
!! ANALYZE generates data and analyzes it with ASA189.
!
!  Discussion:
!
!    The routine to generate the samples uses parameter A, B and C,
!    while ASA189 prefers a related form MU, THETA.  The calling routine
!    has to figure out how these are related.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of samples to generate.
!
!    Input, real ( kind = 8 ) A, real ( kind = 8 ) B, integer ( kind = 4 ) C,
!    the parameters to use in the beta binomial distribution.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real  ( kind = 8 ) MU_EST, THETA_EST, the estimates of MU and THETA
!    produced by ASA189.
!
  implicit none

  integer ( kind = 4 ), parameter :: mrl = 4
  integer ( kind = 4 ) sample_num

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) ccrit
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) in(sample_num)
  integer ( kind = 4 ) iter
  real ( kind = 8 ) mu_est
  real ( kind = 8 ) mu_se
  integer ( kind = 4 ) rl(mrl,3)
  real ( kind = 8 ) rnl
  integer ( kind = 4 ) sample_i
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta_est
  real ( kind = 8 ) theta_se
  integer ( kind = 4 ) x(sample_num)
!
!  Generate the sample data using the exact parameters A, B.
!
  do sample_i = 1, sample_num
    call beta_binomial_sample ( a, b, c, seed, x(sample_i) )
  end do
!
!  Analyze the sample data, trying to estimate the parameters
!  in the form "MU", "THETA".  Note that ASA189 expects to be told
!  the value of C (although C could vary from one sample to the next).
!
  in(1:sample_num) = c

  iter = 10
  ccrit = 0.001D+00

  call bbml ( sample_num, x, in, rl, mrl, iter, ccrit, mu_est, &
    theta_est, mu_se, theta_se, rnl, ifault )

  return
end
