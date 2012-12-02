program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA266_PRB.
!
!  Discussion:
!
!    ASA266_PRB calls the ASA266 tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA266_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA266 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA266_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests ALNORM, NORMP, NPROB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 16

  real alnorm
  real ccdf1
  real ccdf2
  real ccdf3
  real cdf1
  real cdf2
  real cdf3
  integer i
  real pdf2
  real pdf3
  logical upper
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  ALNORM,'
  write ( *, '(a)' ) '  NORMP, and'
  write ( *, '(a)' ) '  NPROB are routines that compute the cumulative'
  write ( *, '(a)' ) '  density function for the normal distribution.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X  CDF1  1-CDF1'
  write ( *, '(a)' ) '     CDF2  1-CDF2  PDF2'
  write ( *, '(a)' ) '     CDF3  1-CDF3  PDF3'
  write ( *, '(a)' ) ' '

  do i = 1, ntest

    x = 3.0E+00 * real ( i - 1 ) / real ( ntest - 1 )

    upper = .false.
    cdf1 = alnorm ( x, upper )

    upper = .true.
    ccdf1 = alnorm ( x, upper )

    call normp ( x, cdf2, ccdf2, pdf2 )

    call nprob ( x, cdf3, ccdf3, pdf3 )

    write ( *, '(a)' ) ' '
    write ( *,     '(3g14.6)' ) x, cdf1, ccdf1
    write ( *, '(14x,3g14.6)' )    cdf2, ccdf2, pdf2
    write ( *, '(14x,3g14.6)' )    cdf3, ccdf3, pdf3

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests PPND, PPND7, PPND16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 9

  real cdf
  double precision cdf2
  integer i
  integer ifault
  real ppnd
  double precision ppnd16
  real ppnd7
  real x1
  real x2
  double precision x3

  ifault = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  PPND, '
  write ( *, '(a)' ) '  PPND7 and '
  write ( *, '(a)' ) '  PPND16 compute the percentage '
  write ( *, '(a)' ) '  points of the normal distribution.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CDF, PPND(CDF), PPND7(CDF), PPND16(CDF)'
  write ( *, '(a)' ) ' '

  do i = 1, ntest

    cdf = real ( i ) / real ( ntest + 1 )
    cdf2 = dble ( cdf )
    x1 = ppnd ( cdf, ifault )
    x2 = ppnd7 ( cdf, ifault )
    x3 = ppnd16 ( cdf2, ifault )

    write ( *, '(4g14.6)' ) cdf, x1, x2, x3

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests DIGAMMA, PSI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 10

  real digamma
  integer i
  real psi
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  digamma(X) = d ( Log ( Gamma ( X ) ) ) / dX.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DIGAMMA and'
  write ( *, '(a)' ) '  PSI compute the digamma function:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X  DIGAMMA   PSI'
  write ( *, '(a)' ) ' '

  do i = 1, ntest

    x = real ( i ) / real ( ntest )

    write ( *, '(3g14.6)' ) x, digamma ( x ), psi ( x )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests TRIGAMMA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 10

  integer i
  real t
  real x
  real trigamma

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  TRIGAMMA computes the trigamma function:'
  write ( *, '(a)' ) '    trigamma(X) = d**2 ( Log ( Gamma ( X ) ) ) / dX**2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X  TRIGAMMA'
  write ( *, '(a)' ) ' '

  do i = 1, ntest

    x = real ( i ) / real ( ntest )

    t = trigamma ( x )
    write ( *, '(2g14.6)' ) x, t

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests ALNGAM, ALOGAM, GAMMA_LOG, LNGAMMA;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 10

  real alngam
  real alogam
  real gamma_log
  integer i
  integer ifault
  real lngamma
  real log1
  real log2
  real log3
  real log4
  real log5
  real x

  ifault = 0
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  ALNGAM,'
  write ( *, '(a)' ) '  ALOGAM,'
  write ( *, '(a)' ) '  GAMMA_LOG, and'
  write ( *, '(a)' ) '  LNGAMMA compute the logarithm of the gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X  ALNGAM  ALOGAM  GAMMA_LOG LNGAMMA'
  write ( *, '(a)' ) ' '

  do i = 1, ntest

    x = real ( i ) / real ( ntest )

    log1 = alngam ( x, ifault )
    log2 = alogam ( x, ifault )
    log3 = gamma_log ( x )
    log4 = lngamma ( x, ifault )

    write ( *, '(5g14.6)' ) x, log1, log2, log3, log4

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests GAMAIN, GAMMDS, GAMMAD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 10

  real g1
  real g2
  real g3
  real gamain
  real gammad
  real gammds
  integer i
  integer ifault
  integer j
  real p
  real x

  ifault = 0
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  GAMAIN, '
  write ( *, '(a)' ) '  GAMMDS and '
  write ( *, '(a)' ) '  GAMMAD compute the incomplete Gamma integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X  P  GAMMDS  GAMMAD  GAMAIN'
  write ( *, '(a)' ) ' '

  do i = 1, ntest

    x = real ( i ) / real ( ntest )

    write ( *, '(a)' ) ' '

    do j = 1, ntest

      p = real ( j ) / real ( ntest )
      g1 = gammds ( x, p, ifault )
      if ( ifault /= 0 ) then
        g1 = -99.0E+00
      end if

      g2 = gammad ( x, p, ifault )
      if ( ifault /= 0 ) then
        g2 = -99.0E+00
      end if

      g3 = gamain ( x, p, ifault )
      if ( ifault /= 0 ) then
        g3 = - 99.0E+00
      end if

      write ( *, '(5g14.6)' ) x, p, g1, g2, g3

    end do

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests PPCHI2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nitest = 9
  integer, parameter :: njtest = 9

  real cdf
  integer i
  integer ifault
  integer j
  real ppchi2
  real v
  real x1

  ifault = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  PPCHI2 computes the percentage points'
  write ( *, '(a)' ) '  of the chi squared distribution.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CDF, PPCHI2(CDF)'
  write ( *, '(a)' ) ' '

  do j = 1, njtest

    v = real ( j )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  For Chi**2 parameter value ', v
    write ( *, '(a)' ) ' '

    do i = 1, nitest

      cdf = real ( i ) / real ( nitest + 1 )
      x1 = ppchi2 ( cdf, v, ifault )

      write ( *, '(4g14.6)' ) cdf, x1

    end do

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests DIRICHLET_ESTIMATE, DIRICHLET_MEAN, DIRICHLET_VARIANCE.
!
!  Discussion:
!
!    Canned data is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: elem_num = 3
  integer, parameter :: sample_num = 23

  real alpha(elem_num)
  real alpha_sum
  real aminus
  real aplus
  integer elem_i
  real eps
  real g(elem_num)
  integer i
  integer ifault
  integer init
  real mean(elem_num)
  integer niter
  real rlogl
  real s
  integer sample_i
  real v((elem_num*(elem_num+1))/2)
  real vari
  real variance(elem_num)
  real work(2*elem_num)
  real x(sample_num,elem_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  For samples of a Dirichlet PDF,'
  write ( *, '(a)' ) '  DIRICHLET_ESTIMATE estimates the parameters.'
  write ( *, '(a)' ) '  DIRICHLET_MEAN finds the means;'
  write ( *, '(a)' ) '  DIRICHLET_VARIANCE finds the variances;'
!
!  Set the data
!
  x(1,1) = 0.178E+00
  x(1,2) = 0.346E+00
  x(1,3) = 0.476E+00

  x(2,1) = 0.162E+00
  x(2,2) = 0.307E+00
  x(2,3) = 0.531E+00

  x(3,1) = 0.083E+00
  x(3,2) = 0.448E+00
  x(3,3) = 0.469E+00

  x(4,1) = 0.087E+00
  x(4,2) = 0.474E+00
  x(4,3) = 0.439E+00

  x(5,1) = 0.078E+00
  x(5,2) = 0.503E+00
  x(5,3) = 0.419E+00

  x(6,1) = 0.040E+00
  x(6,2) = 0.456E+00
  x(6,3) = 0.504E+00

  x(7,1) = 0.049E+00
  x(7,2) = 0.363E+00
  x(7,3) = 0.588E+00

  x(8,1) = 0.100E+00
  x(8,2) = 0.317E+00
  x(8,3) = 0.583E+00

  x(9,1) = 0.075E+00
  x(9,2) = 0.394E+00
  x(9,3) = 0.531E+00

  x(10,1) = 0.084E+00
  x(10,2) = 0.445E+00
  x(10,3) = 0.471E+00

  x(11,1) = 0.060E+00
  x(11,2) = 0.435E+00
  x(11,3) = 0.505E+00

  x(12,1) = 0.089E+00
  x(12,2) = 0.418E+00
  x(12,3) = 0.493E+00

  x(13,1) = 0.050E+00
  x(13,2) = 0.485E+00
  x(13,3) = 0.465E+00

  x(14,1) = 0.073E+00
  x(14,2) = 0.378E+00
  x(14,3) = 0.549E+00

  x(15,1) = 0.064E+00
  x(15,2) = 0.562E+00
  x(15,3) = 0.374E+00

  x(16,1) = 0.085E+00
  x(16,2) = 0.465E+00
  x(16,3) = 0.450E+00

  x(17,1) = 0.094E+00
  x(17,2) = 0.388E+00
  x(17,3) = 0.518E+00

  x(18,1) = 0.014E+00
  x(18,2) = 0.449E+00
  x(18,3) = 0.537E+00

  x(19,1) = 0.060E+00
  x(19,2) = 0.544E+00
  x(19,3) = 0.396E+00

  x(20,1) = 0.031E+00
  x(20,2) = 0.569E+00
  x(20,3) = 0.400E+00

  x(21,1) = 0.025E+00
  x(21,2) = 0.491E+00
  x(21,3) = 0.484E+00

  x(22,1) = 0.045E+00
  x(22,2) = 0.613E+00
  x(22,3) = 0.342E+00

  x(23,1) = 0.0195E+00
  x(23,2) = 0.526E+00
  x(23,3) = 0.4545E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sampled data:'
  write ( *, '(a)' ) ' '

  do sample_i = 1, sample_num
    write ( *, '(i6,3g14.6)' ) sample_i, x(sample_i,1:elem_num)
  end do
!
!  Compute the observed averages.
!
  call r4col_mean ( sample_num, elem_num, x, mean )

  call r4col_variance ( sample_num, elem_num, x, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Observed means, variances are:'
  write ( *, '(a)' ) ' '
  do elem_i = 1, elem_num
    write ( *, '(i6,2g14.6)' ) elem_i, mean(elem_i), variance(elem_i)
  end do

  init = 1

  call dirichlet_estimate ( elem_num, sample_num, x, sample_num, &
    init, alpha, rlogl, v, g, niter, s, eps, work, ifault )

  if ( ifault /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WARNING!'
    write ( *, '(a)' ) '  DIRICHLET_ESTIMATE error code:'
    write ( *, '(a,i6)' ) '  IFAULT = ', ifault
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index, Estimate, Lower Limit, Upper Limit:'
  write ( *, '(a)' ) ' '

  do elem_i = 1, elem_num
    vari = v((elem_i*(elem_i-1))/2+elem_i)
    aminus = alpha(elem_i) - 1.96E+00 * sqrt ( vari )
    aplus = alpha(elem_i) + 1.96E+00 * sqrt ( vari )
    write ( *, '(i6,3g14.6)' ) elem_i, alpha(elem_i), aminus, aplus
  end do

  call dirichlet_mean ( elem_num, alpha, mean )

  call dirichlet_variance ( elem_num, alpha, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Expected means, variances are:'
  write ( *, '(a)' ) ' '
  do elem_i = 1, elem_num
    write ( *, '(i6,2g14.6)' ) elem_i, mean(elem_i), variance(elem_i)
  end do

  alpha_sum = 0.0E+00
  do elem_i = 1, elem_num
    alpha_sum = alpha_sum + alpha(elem_i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Alpha sum is ', alpha_sum
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NORMALIZED VALUES:'
  write ( *, '(a)' ) '  Index, Estimate, Lower Limit, Upper Limit:'
  write ( *, '(a)' ) ' '

  do elem_i = 1, elem_num
    vari = v((elem_i*(elem_i-1))/2+elem_i)
    aminus = ( alpha(elem_i) - 1.96E+00 * sqrt ( vari ) ) / alpha_sum
    aplus = ( alpha(elem_i) + 1.96E+00 * sqrt ( vari ) ) / alpha_sum
    write ( *, '(i6,3g14.6)' ) elem_i, alpha(elem_i)/alpha_sum, aminus, aplus
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) 'Log likelikhood function = ', rlogl

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests DIRICHLET_ESTIMATE, DIRICHLET_MEAN, DIRICHLET_VARIANCE, DIRICHLET_SAMPLE.
!
!  Discussion:
!
!    Data is generated by sampling a distribution with known parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: elem_num = 3
  integer, parameter :: sample_num = 1000

  real alpha(elem_num)
  real alpha_sum
  real aminus
  real aplus
  integer elem_i
  real eps
  real g(elem_num)
  integer ifault
  integer init
  real mean(elem_num)
  integer niter
  real rlogl
  real s
  integer sample_i
  integer seed
  real v((elem_num*(elem_num+1))/2)
  real vari
  real variance(elem_num)
  real work(2*elem_num)
  real x_sample(sample_num,elem_num)
  real x(elem_num)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  For a Dirichlet distribution,'
  write ( *, '(a)' ) '  DIRICHLET_SAMPLE samples;'
  write ( *, '(a)' ) '  DIRICHLET_MEAN finds the means;'
  write ( *, '(a)' ) '  DIRICHLET_VARIANCE finds the variances;'
  write ( *, '(a)' ) '  DIRICHLET_ESTIMATE estimates the parameters.'
!
!  Report.
!
  alpha(1:3) = (/ 3.22E+00, 20.38E+00, 21.68E+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Distribution parameters are:'
  write ( *, '(a)' ) ' '
  do elem_i = 1, elem_num
    write ( *, '(i6,g14.6)' ) elem_i, alpha(elem_i)
  end do

  call dirichlet_mean ( elem_num, alpha, mean )

  call dirichlet_variance ( elem_num, alpha, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Distribution means, variances are:'
  write ( *, '(a)' ) ' '
  do elem_i = 1, elem_num
    write ( *, '(i6,2g14.6)' ) elem_i, mean(elem_i), variance(elem_i)
  end do
!
!  Sample the distribution.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of samples is ', sample_num

  do sample_i = 1, sample_num

    call dirichlet_sample ( elem_num, alpha, seed, x )

    do elem_i = 1, elem_num
      x_sample(sample_i,elem_i) = x(elem_i)
    end do

  end do
!
!  Print some results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First few samples:'
  write ( *, '(a)' ) ' '

  do sample_i = 1, min ( sample_num, 10 )
    write ( *, '(i6,3g14.6)' ) sample_i, x_sample(sample_i,1:elem_num)
  end do
!
!  Compute means, variances.
!
  call r4col_mean ( sample_num, elem_num, x_sample, mean )

  call r4col_variance ( sample_num, elem_num, x_sample, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Observed means, variances are:'
  write ( *, '(a)' ) ' '
  do elem_i = 1, elem_num
    write ( *, '(i6,2g14.6)' ) elem_i, mean(elem_i), variance(elem_i)
  end do
!
!  Destroy the values of ALPHA.
!
  alpha(1:elem_num) = 0.0E+00
!
!  Try to recover the values of ALPHA.
!
  init = 1

  call dirichlet_estimate ( elem_num, sample_num, x_sample, sample_num, &
    init, alpha, rlogl, v, g, niter, s, eps, work, ifault )

  if ( ifault /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WARNING!'
    write ( *, '(a)' ) '  DIRICHLET_ESTIMATE error code:'
    write ( *, '(a,i6)' ) 'IFAULT = ', ifault
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index, Estimate, Lower Limit, Upper Limit:'
  write ( *, '(a)' ) ' '

  do elem_i = 1, elem_num

    vari = v((elem_i*(elem_i-1))/2+elem_i)
    aminus = alpha(elem_i) - 1.96E+00 * sqrt ( vari )
    aplus = alpha(elem_i) + 1.96E+00 * sqrt ( vari )
    write ( *, '(i6,3g14.6)' ) elem_i, alpha(elem_i), aminus, aplus

  end do

  alpha_sum = 0.0E+00
  do elem_i = 1, elem_num
    alpha_sum = alpha_sum + alpha(elem_i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Alpha sum is ', alpha_sum
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NORMALIZED VALUES:'
  write ( *, '(a)' ) '  Index, Estimate, Lower Limit, Upper Limit:'
  write ( *, '(a)' ) ' '

  do elem_i = 1, elem_num
    vari = v((elem_i*(elem_i-1))/2+elem_i)
    aminus = ( alpha(elem_i) - 1.96 * sqrt ( vari ) ) / alpha_sum
    aplus = ( alpha(elem_i) + 1.96 * sqrt ( vari ) ) / alpha_sum
    write ( *, '(i6,3g14.6)' ) elem_i, alpha(elem_i)/alpha_sum, aminus, aplus
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) 'Log likelikhood function = ', rlogl

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests DIRICHLET_MIX_SAMPLE, DIRICHLET_MIX_MEAN, DIRICHLET_MIX_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: comp_num = 3
  integer, parameter :: comp_max = comp_num
  integer, parameter :: elem_num = 3
  integer, parameter :: sample_num = 200

  real a(elem_num)
  real alpha(comp_max,elem_num)
  integer comp_sample(sample_num)
  integer comp_i
  real comp_weight(comp_num)
  integer elem_i
  real mean(elem_num)
  integer sample_i
  integer seed
  real variance(elem_num)
  real x(elem_num)
  real x_sample(sample_num,elem_num)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  For a Dirichlet mixture distribution,'
  write ( *, '(a)' ) '  DIRICHLET_MIX_SAMPLE samples;'
  write ( *, '(a)' ) '  DIRICHLET_MIX_MEAN computes means;'
  write ( *, '(a)' ) '  DIRICHLET_MIX_VARIANCE computes variances.'
!
!  Report.
!
  alpha(1,1) = 0.05E+00
  alpha(1,2) = 0.20E+00
  alpha(1,3) = 0.75E+00

  alpha(2,1) = 0.85E+00
  alpha(2,2) = 0.10E+00
  alpha(2,3) = 0.05E+00

  alpha(3,1) = 0.00E+00
  alpha(3,2) = 0.50E+00
  alpha(3,3) = 0.50E+00

  comp_weight(1:3) = (/ 3.0E+00, 2.0E+00, 1.0E+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Component Weight'
  write ( *, '(a)' ) ' '
  do comp_i = 1, comp_num
    write ( *, '(i6,1x,f10.6,4x,3f10.6)' ) comp_i, comp_weight(comp_i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Component  Parameters Means Variances'
  write ( *, '(a)' ) ' '
  do comp_i = 1, comp_num
    write ( *, * ) ' '
    write ( *, '(i6)' ) comp_i
    do elem_i = 1, elem_num
      a(elem_i) = alpha(comp_i,elem_i)
    end do
    call dirichlet_mean ( elem_num, a, mean )
    call dirichlet_variance ( elem_num, a, variance )
    do elem_i = 1, elem_num
      write ( *, '(i6,1x,3f10.6)' ) elem_i, alpha(comp_i,elem_i), &
        mean(elem_i), variance(elem_i)
    end do
  end do

  call dirichlet_mix_mean ( comp_max, comp_num, elem_num, alpha, &
    comp_weight, mean )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Element  Mean'
  write ( *, '(a)' ) ' '
  do elem_i = 1, elem_num
    write ( *, '(i6,1x,2f10.6)' ) elem_i, mean(elem_i)
  end do
!
!  Sample the distribution.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of samples is ', sample_num

  do sample_i = 1, sample_num

    call dirichlet_mix_sample ( comp_max, comp_num, elem_num, alpha, &
      comp_weight, seed, comp_i, x )

    x_sample(sample_i,1:elem_num) = x(1:elem_num)

    comp_sample(sample_i) = comp_i

  end do
!
!  Print some results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'First few samples:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Sample  Component  X'
  write ( *, '(a)' ) ' '

  do sample_i = 1, min ( sample_num, 10 )

    write ( *, '(2x,i2,2x,i2,2x,3f10.6)' ) sample_i, comp_sample(sample_i), &
      x_sample(sample_i,1:elem_num)

  end do
!
!  Compute the observed averages.
!
  call r4col_mean ( sample_num, elem_num, x_sample, mean )

  call r4col_variance ( sample_num, elem_num, x_sample, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Element  Observed mean, variance'
  write ( *, '(a)' ) ' '
  do elem_i = 1, elem_num
    write ( *, '(i6,1x,2f10.6)' ) elem_i, mean(elem_i), variance(elem_i)
  end do

  return
end







