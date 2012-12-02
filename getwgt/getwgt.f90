subroutine getwgt ( aacnts, pseudocount )

!*****************************************************************************80
!
!! GETWGT updates the Dirichlet mixture weights based on a set of counts.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer AACNTS(ACID_NUM), the number of counts of each amino acid.
!    The implicit order used is: 
!
!      1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
!      A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
!
!    Output, real PSEUDOCOUNT(ACID_NUM), the estimated pseudocount vector.
!
!  Internals:
!
!    Internal, integer ACID_NUM, the number of amino acids, which is 20.
!
!    Internal, integer COMP_NUM, the number of mixture components,
!    which is assumed to be 10.
!
!    Internal, real ALPHA(COMP_NUM), the estimated Dirichlet parameters
!    for the mixture.
!
!    Internal, real BETA(ACID_NUM,COMP_NUM); BETA(I,J) is the parameter
!    for the J-th acid in the I-th Dirichlet mixture component.
!
!    Internal, COMP_WEIGHT_ESTIMATE(COMP_NUM); the estimated value of the
!    weights for the components of the mixture.
!
!    Internal, integer NCALL, the number of times this routine has been called.
!
!    Internal, real P(COMP_NUM), P(I) is the Bayesian posterior probability 
!    of component I, given the observation of the most recent event, which is 
!    proportional to the probability of the event under the component I PDF,
!    times the prior probability of component I.
!
!    Internal, real P_HAT(COMP_NUM), the prior probabilities of the
!    components.
!
  implicit none

  integer ( kind = 4 ), parameter :: acid_num = 20
  integer ( kind = 4 ), parameter :: comp_max = 10

  integer ( kind = 4 ) aacnts(acid_num)
  integer ( kind = 4 ) acid_i
  character acid_sym(acid_num)
  real ( kind = 4 ), save, dimension ( comp_max ) :: alpha
  real ( kind = 4 ), save, dimension ( acid_num, comp_max ) :: beta
  real ( kind = 4 ) beta_sum(comp_max)
  integer ( kind = 4 ) comp_i
  integer ( kind = 4 ) comp_label(comp_max)
  integer ( kind = 4 ), save :: comp_num = 0
  real ( kind = 4 ) comp_weight(comp_max)
  real ( kind = 4 ), save, dimension ( comp_max ) :: comp_weight_estimate
  character ( len = 80 ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ), save :: ncall = 0
  real ( kind = 4 ), save, dimension ( comp_max ) :: p
  real ( kind = 4 ), save, dimension ( comp_max ) :: p_hat
  real ( kind = 4 ) pseudocount(acid_num)
  integer ( kind = 4 ) site_num

  ncall = ncall + 1
!
!  On first call, we need to retrieve some information.
!
  if ( ncall == 1 ) then

    file_name = 'weights.txt'
    iunit = 1

    open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
      err = 10 )

    call mixture_read ( acid_num, acid_sym, beta, beta_sum, comp_label, &
      comp_max, comp_num, comp_weight, ierror, iunit )

    close ( unit = iunit )

    if ( comp_num > comp_max ) then
      write ( *, * ) ' '
      write ( *, * ) 'GETWGT - Fatal error!'
      write ( *, * ) '  Number of components exceeds internal max.'
      stop
    end if

    call comp_param_print ( acid_num, acid_sym, comp_max, comp_num, beta, &
      beta_sum, comp_weight )

  end if
!
!  Initialize the ALPHA's.
!
!  Options:
!    Use the weights stored in the file as ALPHA's.
!    Use equal ALPHA's of 1 / comp_num.
!    Multiply these values to increase the importance of the initial values.
!
  alpha(1:comp_num) = comp_weight(1:comp_num)

  call r4vec_print ( comp_num, alpha, 'Old Alphas' )

  call r4vec_print ( comp_num, comp_weight, 'Old Weights' )
!
!  Initialize the prior probabilities.
!
  p(1:comp_num) = comp_weight(1:comp_num)
  p_hat(1:comp_num) = comp_weight(1:comp_num)

  site_num = sum ( aacnts(1:acid_num) )
!
!  Process the new information.
!
  call event_process ( acid_num, alpha, beta, comp_max, comp_num, p, p_hat, &
    site_num, aacnts )
!
!  From the new ALPHA's we update the estimated weights.
!
  call r4vec_print ( comp_num, alpha, 'New Alphas' )

  call dirichlet_mean ( comp_num, alpha, comp_weight_estimate )

  call r4vec_print ( comp_num, comp_weight_estimate, 'New Weights' )
!
!  From the updated weight estimates, compute the corresponding pseudo count.
!
  do acid_i = 1, acid_num
    pseudocount(acid_i) = 0.0
    do comp_i = 1, comp_num
      pseudocount(acid_i) = pseudocount(acid_i) &
        + comp_weight_estimate(comp_i) * beta(acid_i,comp_i)
    end do
  end do

  return
!
!  Error: could not open the data file.
!
10    continue

  write ( *, * ) ' '
  write ( *, * ) 'GETWGT - Fatal error!'
  write ( *, * ) '  Could not open the data file:'
  write ( *, '(a)' ) file_name

  ncall = 0

  stop
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c2
  character cc1
  character cc2

  cc1 = c1
  cc2 = c2

  call ch_cap ( cc1 )
  call ch_cap ( cc2 )

  if ( cc1 == cc2 ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_next ( line, cval, done )

!*****************************************************************************80
!
!! CH_NEXT "reads" space-separated characters from a string, one at a time.
!
!  Example:
!
!    Input:
!
!      LINE = ' A  B, C    DE  F'
!
!    Output:
!
!      'A', 'B', 'C', 'D', 'E', 'F', and then blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing
!    characters, possibly separated by spaces or commas.
!
!    Output, character CVAL.  If DONE is FALSE, then CVAL contains the
!    "next" character read from LINE.  If DONE is TRUE, then
!    CVAL is blank.
!
!    Input/output, logical DONE.
!    On input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!    On output, the routine sets DONE to FALSE if another character
!    was read, or TRUE if no more characters could be read.
!
  implicit none

  character cval
  logical done
  integer ( kind = 4 ) i
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1

  if ( done ) then
    next = 1
    done = .false.
  end if

  do i = next, len_trim ( line )

    if ( line(i:i) /= ' ' .and. line(i:i) /= ',' ) then
      cval = line(i:i)
      next = i + 1
      return
    end if

  end do

  done = .true.
  next = 1
  cval = ' '

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine comp_param_print ( acid_num, acid_sym, comp_max, comp_num, beta, &
  beta_sum, comp_weight )

!*****************************************************************************80
!
!! COMP_PARAM_PRINT prints the parameters for the mixture components.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ACID_NUM, the number of amino acids.
!
!    Input, character ACID_SYM(ACID_NUM), the one letter amino acid codes.
!
!    Input, integer COMP_MAX, the maximum number of Dirichlet mixture
!    components.
!
!    Input, integer COMP_NUM, the number of components in the Dirichlet mixture.
!
!    Input, real BETA(ACID_NUM,COMP_MAX); BETA(I,J) is the parameter
!    for the J-th acid in the I-th Dirichlet mixture component.
!
!    Input, real BETA_SUM(COMP_MAX), the sum of the values of
!    BETA(ACID_I,COMP_I) for a given component COMP_I.
!
!    Input, real COMP_WEIGHT(COMP_NUM), the mixture weight of each component.
!    These values should be nonnegative, and sum to 1.  They represent the
!    relative proportion of each component in the mixture.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_max

  integer ( kind = 4 ) acid_i
  character acid_sym(acid_num)
  integer ( kind = 4 ) comp_i
  real ( kind = 4 ) beta(acid_num,comp_max)
  real ( kind = 4 ) beta_sum(comp_max)
  integer ( kind = 4 ) comp_num
  real ( kind = 4 ) comp_weight(comp_max)

  write ( *, * ) ' '
  write ( *, * ) '  Number of components = ', comp_num
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, '(''Compon:'',20i8)' ) ( comp_i, comp_i = 1, comp_num )
  write ( *, '(''Weight:'',20f8.4)' ) comp_weight(1:comp_num)
  write ( *, * ) ' '

  do acid_i = 1, acid_num
    write ( *, '(i2,2x,a1,2x,20f8.4)' ) acid_i, acid_sym(acid_i), &
      beta(acid_i,1:comp_num)
  end do

  write ( *, * ) ' '
  write ( *, '(a3,4x,20f8.4)' ) 'Sum', beta_sum(1:comp_num)

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
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
!    Output, real MEAN(N), the means of the PDF.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) a(n)
  real ( kind = 4 ) mean(n)

  mean(1:n) = a(1:n)

  call r4vec_unit_sum ( n, mean )

  return
end
subroutine dirichlet_multinomial_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! DIRICHLET_MULTINOMIAL_PDF evaluates a Dirichlet Multinomial PDF.
!
!  Discussion:
!
!    PDF(X)(A,B,C) = Comb(A,B,X) * ( Gamma(C_Sum) / Gamma(C_Sum+A) )
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
!    Input, integer X(B); X(I) counts the number of occurrences of
!    outcome I, out of the total of A trials.
!
!    Input, integer A, the total number of trials.
!
!    Input, integer B, the number of different possible outcomes on
!    one trial.
!
!    Input, integer C(B); C(I) is the Dirichlet parameter associated
!    with outcome I.
!
!    Output, real PDF, the value of the Dirichlet multinomial PDF.
!
  implicit none

  integer ( kind = 4 ) b

  integer ( kind = 4 ) a
  real ( kind = 4 ) c(b)
  real ( kind = 4 ) c_sum
  real ( kind = 4 ) gamma_log
  integer ( kind = 4 ) i
  real ( kind = 4 ) pdf
  real ( kind = 4 ) pdf_log
  integer ( kind = 4 ) x(b)

  c_sum = sum ( c(1:b) )

  pdf_log = - gamma_log ( c_sum + real ( a ) ) + gamma_log ( c_sum ) &
    + gamma_log ( real ( a + 1 ) )

  do i = 1, b

    pdf_log = pdf_log + gamma_log ( c(i) + real ( x(i) ) ) &
      - gamma_log ( c(i) ) - gamma_log ( real ( x(i) + 1 ) )
  end do

  pdf = exp ( pdf_log )

  return
end
subroutine event_process ( acid_num, alpha, beta, comp_max, comp_num, p, &
  p_hat, site_num, x_sample )

!*****************************************************************************80
!
!! EVENT_PROCESS updates the mixture weight distribution parameters.
!
!  Discussion:
!
!    This routine updates the values of ALPHA.  It does this by
!    considering the results of the most recent event.  If we knew
!    which component PDF had generated the event, then we would 
!    simply add 1 to the ALPHA for that component.  Instead, we
!    use Bayesian analysis to estimate the proportion of the event
!    that is to be added to each ALPHA.
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
!  Reference:
!
!    BS Everitt, DJ Hand,
!    Finite Mixture Distributions,
!    Chapman and Hall, 1981.
!
!    AFM Smith, UE Makov,
!    A Quasi-Bayes Sequential Procedure for Mixtures,
!    Journal of the Royal Statistical Society,
!    Volume 40, Number 1, B, 1978, pages 106-112.
!
!  Parameters:
!
!    Input, integer ACID_NUM, the number of amino acids.
!
!    Input/output, real ALPHA(COMP_NUM), the Dirichlet parameters for the 
!    weights.
!
!    Input, real BETA(ACID_NUM,COMP_MAX); BETA(I,J) is the multinomial Dirichlet
!    parameter for the J-th acid in the I-th Dirichlet mixture component.
!
!    Input, integer COMP_NUM, the number of components in the Dirichlet mixture.
!
!    Input/output, real P(COMP_NUM); P(I) is the Bayesian posterior probability 
!    of component I, given the observation of the most recent event, which is 
!    proportional to the probability of the event under the component I PDF,
!    times the prior probability of component I.
!
!    Input/output, real P_HAT(COMP_NUM), the prior probabilities of the
!    components.
!
!    Input, integer SITE_NUM, the number of sites observed for this event.
!    This value might change from call to call, although in the demonstration
!    I'm keeping it fixed.
!
!    Input, integer X_SAMPLE(ACID_NUM), the "current event", namely, the count
!    vector for the number of occurrences of each acid out of the total
!    of SITE_NUM sites analyzed.  This is the evidence used to update the
!    "theory" for the value of ALPHA.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_max

  real ( kind = 4 ) alpha(comp_max)
  real ( kind = 4 ) alpha_sum
  real ( kind = 4 ) beta(acid_num,comp_max)
  integer ( kind = 4 ) comp_i
  integer ( kind = 4 ) comp_num
  real ( kind = 4 ) comp_pdf
  real ( kind = 4 ) comp_sum
  real ( kind = 4 ) p(comp_max)
  real ( kind = 4 ) p_hat(comp_max)
  integer ( kind = 4 ) site_num
  real ( kind = 4 ) sum
  integer ( kind = 4 ) x_sample(acid_num)
!
!  Sum the parameters.
!
  alpha_sum = sum ( alpha(1:comp_num) )
!
!  Update P_HAT.
!
  p_hat(1:comp_num) = ( ( alpha_sum - 1.0 ) * p_hat(1:comp_num) &
    + p(1:comp_num) ) / alpha_sum
!
!  Generate the new P's.
!  P(COMP_I) = the Bayesian posterior probability of component I,
!  given the observation of event EVENT_I, which is proportional
!  to the probability of event EVENT_I in the component I PDF,
!  times the prior probability of component I.
!
  do comp_i = 1, comp_num
!
!  Compute the probability of this event, for a given component.
!
    call dirichlet_multinomial_pdf ( x_sample, site_num, acid_num, &
      beta(1,comp_i), comp_pdf )
!
!  Multiply by the probability of that component to get the relative
!  probability of the event.
!
    p(comp_i) = comp_pdf * p_hat(comp_i)

  end do

  write ( *, * ) ' '
  write ( *, * ) 'BAYES'
  write ( *, * ) ' '
  comp_sum = 0.0E+00
  do comp_i = 1, comp_num
    comp_sum = comp_sum + p(comp_i)
  end do

  do comp_i = 1, comp_num
    call dirichlet_multinomial_pdf ( x_sample, site_num, acid_num, &
      beta(1,comp_i), comp_pdf )
    write ( *, '(i4,4g14.6)' ) comp_i, comp_pdf, p_hat(comp_i), p(comp_i), &
      p(comp_i) / comp_sum
  end do
!
!  Normalize the P's to get the absolute Bayesian probability.
!
  call r4vec_unit_sum ( comp_num, p )
!
!  Update the alpha's by adding adding appropriate portions of
!  the most recent event to each component's parameter.
!
  alpha(1:comp_num) = alpha(1:comp_num) + p(1:comp_num)

  return
end
function gamma_log ( x )

!*****************************************************************************80
!
!! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 June 1999
!
!  Authors:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    W. J. Cody and K. E. Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Math. Comp.
!    Volume 21, 1967, pages 198-203.
!
!    K. E. Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    Hart, Et. Al.,
!    Computer Approximations,
!    Wiley and sons, New York, 1968.
!
!  Parameters:
!
!    Input, real X, the argument of the Gamma function.  X must be positive.
!
!    Output, real GAMMA_LOG, the logarithm of the Gamma function of X.
!    If X <= 0.0, or if overflow would occur, the program returns the
!    value XINF, the largest representable floating point number.
!
!
!  Explanation of machine-dependent constants
!
!  BETA   - radix for the floating-point representation.
!
!  MAXEXP - the smallest positive power of BETA that overflows.
!
!  XBIG   - largest argument for which LN(GAMMA(X)) is representable
!           in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!  XINF   - largest machine representable floating-point number;
!           approximately BETA**MAXEXP.
!
!  EPS    - The smallest positive floating-point number such that
!           1.0+EPS .GT. 1.0
!
!  FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!  Approximate values for some important machines are:
!
!                            BETA      MAXEXP         XBIG
!
!  CRAY-1        (S.P.)        2        8191       9.62E+2461
!  Cyber 180/855
!    under NOS   (S.P.)        2        1070       1.72E+319
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)        2         128       4.08E+36
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!  IBM 3033      (D.P.)       16          63       4.29D+73
!  VAX D-Format  (D.P.)        2         127       2.05D+36
!  VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                            XINF        EPS        FRTBIG
!
!  CRAY-1        (S.P.)   5.45E+2465   7.11E-15    3.13E+615
!  Cyber 180/855
!    under NOS   (S.P.)   1.26E+322    3.55E-15    6.44E+79
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   3.40E+38     1.19E-7     1.42E+9
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)   1.79D+308    2.22D-16    2.25D+76
!  IBM 3033      (D.P.)   7.23D+75     2.22D-16    2.56D+18
!  VAX D-Format  (D.P.)   1.70D+38     1.39D-17    1.20D+9
!  VAX G-Format  (D.P.)   8.98D+307    1.11D-16    1.89D+76
!
  implicit none

  real ( kind = 4 ) EPS
  real ( kind = 4 ) FRTBIG
  real ( kind = 4 ) PNT68
  real ( kind = 4 ) SQRTPI
  real ( kind = 4 ) XBIG
  real ( kind = 4 ) XINF

  parameter ( EPS = 1.19e-7 )
  parameter ( FRTBIG = 1.42e9 )
  parameter ( PNT68 = 0.6796875 )
  parameter ( SQRTPI = 0.9189385332046727417803297 )
  parameter ( XBIG = 4.08e36 )
  parameter ( XINF = 3.401e38 )

  real ( kind = 4 ), parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728e-03, &
     8.4171387781295e-04, &
    -5.952379913043012e-04, &
     7.93650793500350248e-04, &
    -2.777777777777681622553e-03, &
     8.333333333333333331554247e-02, &
     5.7083835261e-03 /)
  real ( kind = 4 ) corr
  integer ( kind = 4 ) i
  real ( kind = 4 ), parameter :: d1 = -5.772156649015328605195174e-1
  real ( kind = 4 ), parameter :: d2 = 4.227843350984671393993777e-1
  real ( kind = 4 ), parameter :: d4 = 1.791759469228055000094023e0
  real ( kind = 4 ) gamma_log
  real ( kind = 4 ), parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888e0, &
    2.018112620856775083915565e2, &
    2.290838373831346393026739e3, &
    1.131967205903380828685045e4, &
    2.855724635671635335736389e4, &
    3.848496228443793359990269e4, &
    2.637748787624195437963534e4, &
    7.225813979700288197698961e3 /)
  real ( kind = 4 ), parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064e0, &
    5.424138599891070494101986e2, &
    1.550693864978364947665077e4, &
    1.847932904445632425417223e5, &
    1.088204769468828767498470e6, &
    3.338152967987029735917223e6, &
    5.106661678927352456275255e6, &
    3.074109054850539556250927e6 /)
  real ( kind = 4 ), parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062e4, &
    2.426813369486704502836312e6, &
    1.214755574045093227939592e8, &
    2.663432449630976949898078e9, &
    2.940378956634553899906876e10, &
    1.702665737765398868392998e11, &
    4.926125793377430887588120e11, &
    5.606251856223951465078242e11 /)
  real ( kind = 4 ), parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036e1, &
    1.113332393857199323513008e3, &
    7.738757056935398733233834e3, &
    2.763987074403340708898585e4, &
    5.499310206226157329794414e4, &
    6.161122180066002127833352e4, &
    3.635127591501940507276287e4, &
    8.785536302431013170870835e3 /)
  real ( kind = 4 ), parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942e2, &
    7.765049321445005871323047e3, &
    1.331903827966074194402448e5, &
    1.136705821321969608938755e6, &
    5.267964117437946917577538e6, &
    1.346701454311101692290052e7, &
    1.782736530353274213975932e7, &
    9.533095591844353613395747e6 /)
  real ( kind = 4 ), parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843e3, &
    6.393885654300092398984238e5, &
    4.135599930241388052042842e7, &
    1.120872109616147941376570e9, &
    1.488613728678813811542398e10, &
    1.016803586272438228077304e11, &
    3.417476345507377132798597e11, &
    4.463158187419713286462081e11 /)
  real ( kind = 4 ) res
  real ( kind = 4 ) x
  real ( kind = 4 ) xden
  real ( kind = 4 ) xm1
  real ( kind = 4 ) xm2
  real ( kind = 4 ) xm4
  real ( kind = 4 ) xnum
  real ( kind = 4 ) xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0 .or. x > XBIG ) then
    gamma_log = XINF
    return
  end if

  if ( x <= EPS ) then

    res = - log ( x )

  else if ( x <= 1.5 ) then

    if ( x < PNT68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0
      xm1 = ( x - 0.5 ) - 0.5
    end if

    if ( x <= 0.5 .or. x >= PNT68 ) then

      xden = 1.0
      xnum = 0.0

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5 ) - 0.5
      xden = 1.0
      xnum = 0.0
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0 ) then

    xm2 = x - 2.0
    xden = 1.0
    xnum = 0.0
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0 ) then

    xm4 = x - 4.0
    xden = - 1.0
    xnum = 0.0
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0

    if ( x <= FRTBIG ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + SQRTPI - 0.5 * corr
    res = res + x * ( corr - 1.0 )

  end if

  gamma_log = res

  return
end
subroutine i4_next ( line, ival, done )

!*****************************************************************************80
!
!! I4_NEXT "reads" integers from a string, one at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing
!    integers.  These may be separated by spaces or commas.
!
!    Output, integer IVAL.  If DONE is FALSE, then IVAL contains the
!    "next" integer read from LINE.  If DONE is TRUE, then
!    IVAL is zero.
!
!    Input/output, logical DONE.
!    On input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!    On output, the routine sets DONE to FALSE if another integer
!    was read, or TRUE if no more integers could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) lchar
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1

  ival = 0

  if ( done ) then
    next = 1
    done = .false.
  end if

  if ( next > len_trim ( line ) ) then
    done = .true.
    return
  end if

  call s_to_i4 ( line(next:), ival, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    done = .true.
    next = 1
  else
    done = .false.
    next = next + lchar
  end if

  return
end
subroutine mixture_read ( acid_num, acid_sym, beta, beta_sum, comp_label, &
  comp_max, comp_num, comp_weight, ierror, iunit )

!*****************************************************************************80
!
!! MIXTURE_READ reads the Dirichlet mixture parameters from a file.
!
!  Discussion:
!
!    The data in the file is delimited by keywords.
!
!    The first lines (not necessarily in order!) may include
!
!      ClassName = string
!      NumDistr = N           the number of components in the mixture.
!      Alphabet = string
!      Order = A C D E ...    the order of the amino acids.
!      AlphaChar = 20
!      NumDistr = 9           the number of distributions
!      EndClassName = string
!
!    For each component, there are four lines:
!
!      Number= N              the component number, starting with 0
!      Mixture= N             the mixture weight, out of a total of 1.0
!      Alpha=  |A| A1 A2 ...  the parameter sum, and individual parameters
!      Comment=               a comment, which describes the frequencies.
!
!    In the comment, the symbol "><" indicates the mean background frequency;
!    residues to the left of that symbol occur more frequently
!    than background, residues to the right less frequently.  Commas separate 
!    residues differing in frequency by a factor of 2.
!
!    For example, the comment
!      S A T , C G P >< N V M , Q H R I K F L D W , E Y
!    indicates that for this component, the frequency of
!    proline is just above the mean, and serine, alanine and
!    threonine are twice as frequent in this component than they
!    are on average.  By contrast, tyrosine and glutamic acid are
!    between 4 and 8 times less likely in this component than on
!    average.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ACID_NUM, the number of amino acids.
!
!    Output, character ACID_SYM(ACID_NUM), the one letter amino acid codes.
!
!    Output, real BETA(ACID_NUM,COMP_MAX); BETA(I,J) is the parameter
!    for the J-th acid in the I-th Dirichlet mixture component.
!
!    Output, real BETA_SUM(COMP_MAX), the sum of the values of
!    BETA(ACID_I,COMP_I) for a given component COMP_I.
!
!    Output, integer COMP_LABEL(COMP_NUM), the label of each component.
!    Normally, component I has label I.
!
!    Input, integer COMP_MAX, the maximum number of Dirichlet mixture
!    components.
!
!    Output, integer COMP_NUM, the number of components in the Dirichlet
!    mixture.
!
!    Output, real COMP_WEIGHT(COMP_NUM), the mixture weight of each component.
!    These values should be nonnegative, and sum to 1.  They represent the
!    relative proportion of each component in the mixture.
!
!    Output, integer IERROR, error indicator.
!    0: no error occurred; nonzero: an error occurred.
!
!    Input, integer IUNIT, the FORTRAN unit from which the data is to
!    be read.
!
  implicit none

  integer ( kind = 4 ) acid_num
  integer ( kind = 4 ) comp_max

  integer ( kind = 4 ) acid_i
  character acid_sym(acid_num)
  real ( kind = 4 ) beta(acid_num,comp_max)
  real ( kind = 4 ) beta_sum(comp_max)
  integer ( kind = 4 ) comp_i
  integer ( kind = 4 ) comp_label(comp_max)
  integer ( kind = 4 ) comp_num
  real ( kind = 4 ) comp_weight(comp_max)
  logical done
  integer ( kind = 4 ) iequal
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ngoofy
  integer ( kind = 4 ) nrec
  logical s_begin
  character ( len = 500 ) string

  ierror = 0
  comp_i = 0
  comp_num = 0
  nrec = 0
  ngoofy = 0

10    continue

  read ( iunit, '(a)', end = 20 ) string
  nrec = nrec + 1
!
!  Ignore blank lines.
!
  if ( string == ' ' ) then
!
!  Ignore the CLASSNAME field.
!
  else if ( s_begin ( string, 'CLASSNAME' ) ) then
!
!  Ignore the ENDCLASSNAME field.
!
  else if ( s_begin ( string, 'ENDCLASSNAME' ) ) then
!
!  Ignore the NAME field.
!
  else if ( s_begin ( string, 'NAME' ) ) then
!
!  Ignore the ALPHABET field.
!
  else if ( s_begin ( string, 'ALPHABET' ) ) then
!
!  Read the ORDER field, since it tells us how to interpret the ALPHA's.
!
  else if ( s_begin ( string, 'ORDER' ) ) then

    iequal = index ( string, '=' )
    done = .true.
    do acid_i = 1, acid_num
      call ch_next ( string(iequal+1:), acid_sym(acid_i), done )
    end do
!
!  Ignore the ALPHACHAR field.
!
  else if ( s_begin ( string, 'ALPHACHAR' ) ) then
!
!  Read the NUMDISTR field.
!
  else if ( s_begin ( string, 'NUMDISTR' ) ) then

    iequal = index ( string, '=' )
    done = .true.
    call i4_next ( string(iequal+1:), comp_num, done )

    if ( comp_num < 1 ) then
      ierror = 1
      return
    else if ( comp_num > comp_max ) then
      ierror = 2
      return
    end if
!
!  Read the NUMBER field.
!
  else if ( s_begin ( string, 'NUMBER' ) ) then

    comp_i = comp_i + 1

    if ( comp_i > comp_num ) then
      write ( *, * ) ' '
      write ( *, * ) 'MIXTURE_READ - Fatal error!'
      write ( *, * ) '  Number of components = ', comp_i
      write ( *, * ) '  exceeding reported value of ', comp_num
      stop
    end if

    iequal = index ( string, '=' )
    done = .true.
    call i4_next ( string(iequal+1:), comp_label(comp_i), done )
!
!  Read the MIXTURE field.
!
  else if ( s_begin ( string, 'MIXTURE' ) ) then

    iequal = index ( string, '=' )
    done = .true.
    call r4_next ( string(iequal+1:), comp_weight(comp_i), done )
!
!  Read the ALPHA field.
!
  else if ( s_begin ( string, 'ALPHA' ) ) then

    iequal = index ( string, '=' )
    done = .true.
    call r4_next ( string(iequal+1:), beta_sum(comp_i), done )
    do acid_i = 1, acid_num
      call r4_next ( string(iequal+1:), beta(acid_i,comp_i), done )
    end do
!
!  Ignore the COMMENT field.
!
  else if ( s_begin ( string, 'COMMENT' ) ) then
!
!  Unexpected field:
!
  else

    write ( *, * ) ' '
    write ( *, * ) 'MIXTURE_READ - Warning!'
    write ( *, * ) '  Goofy record: '
    write ( *, '(a)' ) string(1:20)

    ngoofy = ngoofy + 1

  end if

  go to 10

20    continue

  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) 'MIXTURE_READ - Note:'
  write ( *, * ) '  Number of records read was ', nrec
  write ( *, * ) '  Number of goofy records was ', ngoofy

  return
end
subroutine r4_next ( line, rval, done )

!*****************************************************************************80
!
!! R4_NEXT "reads" real numbers from a string, one at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing real
!    numbers.  These may be separated by spaces or commas.
!
!    Output, real RVAL.  If DONE is FALSE, then RVAL contains the
!    "next" real value read from LINE.  If DONE is TRUE, then
!    RVAL is zero.
!
!    Input/output, logical DONE.
!    On input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!    On output, the routine sets DONE to FALSE if another real
!    value was read, or TRUE if no more reals could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) lchar
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1
  real ( kind = 4 ) rval

  rval = 0.0

  if ( done ) then
    next = 1
    done = .false.
  end if

  if ( next > len_trim ( line ) ) then
    done = .true.
    return
  end if

  call s_to_r4 ( line(next:), rval, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    done = .true.
    next = 1
  else
    done = .false.
    next = next + lchar
  end if

  return
end
subroutine r4vec_print ( n, a, title )

!*****************************************************************************80
!
!! R4VEC_PRINT prints an R4VEC, with an optional title.
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
!    Input, integer N, the number of components of the vector.
!
!    Input, real A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
  end do

  return
end
subroutine r4vec_unit_sum ( n, a )

!*****************************************************************************80
!
!! R4VEC_UNIT_SUM normalizes an R4VEC to have unit sum.
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
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, real A(N), the vector to be normalized.  On output,
!    the entries of A should have unit sum.  However, if the input vector
!    has zero sum, the routine halts.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) a(n)

  a(1:n) = a(1:n) / sum ( a(1:n) )

  return
end
function s_begin ( s1, s2 )

!*****************************************************************************80
!
!! S_BEGIN is TRUE if one string matches the beginning of the other.
!
!  Example:
!
!     S1              S2      S_BEGIN
!
!    'Bob'          'BOB'     TRUE
!    '  B  o b '    ' bo b'   TRUE
!    'Bob'          'Bobby'   TRUE
!    'Bobo'         'Bobb'    FALSE
!    ' '            'Bob'     FALSE    (Do not allow a blank to match
!                                       anything but another blank string.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to be compared.
!
!    Output, logical S_BEGIN, .TRUE. if the strings match up to
!    the end of the shorter string, ignoring case,
!    FALSE otherwise.
!
  implicit none

  logical ch_eqi
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  logical s_begin
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len_trim ( s1 )
  len2 = len_trim ( s2 )
!
!  If either string is blank, then both must be blank to match.
!  Otherwise, a blank string matches anything, which is not 
!  what most people want.
!
  if ( len1 == 0 .or. len2 == 0 ) then

    if ( len1 == 0 .and. len2 == 0 ) then
      s_begin = .true.
    else
      s_begin = .false.
    end if

    return

  end if

  i1 = 0
  i2 = 0
!
!  Find the next nonblank in S1.
!
10    continue

  i1 = i1 + 1

  if ( i1 > len1 ) then
    s_begin = .true.
    return
  end if

  if ( s1(i1:i1) == ' ' ) then
    go to 10
  end if
!
!  Find the next nonblank in S2.
!
20    continue

  i2 = i2 + 1

  if ( i2 > len2 ) then
    s_begin = .true.
    return
  end if

  if ( s2(i2:i2) == ' ' ) then
    go to 20
  end if
!
!  If the characters match, get the next pair.
!
  if ( ch_eqi ( s1(i1:i1), s2(i2:i2) ) ) then
    go to 10
  end if

  s_begin = .false.

  return
end
subroutine s_to_i4 ( string, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, a string to be examined.
!
!    Output, integer IVAL, the integer value read from the string.
!    If STRING is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LAST, the last character in STRING that was
!    part of the representation of IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lens
  character ( len = * ) string

  ierror = 0
  istate = 0

  isgn = 1
  ival = 0

  lens = len_trim ( string )

  i = 0

10    continue

  i = i + 1

  c = string(i:i)

  if ( istate == 0 ) then

    if ( c == ' ' ) then

    else if ( c == '-' ) then
      istate = 1
      isgn = -1
    else if ( c == '+' ) then
      istate = 1
      isgn = + 1
    else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
      istate = 2
      ival = ichar ( c ) - ichar ( '0' )
    else
      ierror = 1
      return
    end if

  else if ( istate == 1 ) then

    if ( c == ' ' ) then

    else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
      istate = 2
      ival = ichar ( c ) - ichar ( '0' )
    else
      ierror = 1
      return
    end if

  else if ( istate == 2 ) then

    if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
      ival = 10 * ival + ichar ( c ) - ichar ( '0' )
    else
      istate = 3
    end if

  end if
!
!  Continue or exit?
!
  if ( istate == 3 ) then
    ival = isgn * ival
    last = i - 1
    return
  else if ( i >= lens ) then
    if ( istate == 2 ) then
      ival = isgn * ival
      last = lens
    else
      ierror = 1
      last = 0
    end if
    return
  end if

  go to 10
end
subroutine s_to_r4 ( string, rval, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R4 reads an R4 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    STRING            RVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real RVAL, the real value that was read from the string.
!
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    STRING to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character chrtmp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 4 ) rbot
  real ( kind = 4 ) rexp
  real ( kind = 4 ) rtop
  real ( kind = 4 ) rval
  character ( len = * ) string
  character TAB

  nchar = len_trim ( string )
  TAB = char(9)
  ierror = 0
  rval = 0.0
  lchar = - 1
  isgn = 1
  rtop = 0.0
  rbot = 1.0
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

10    continue

  lchar = lchar + 1
  chrtmp = string(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
  if ( chrtmp == ' ' .or. chrtmp == TAB ) then
!
!  20 November 1993
!
!  I would like to allow input like "+ 2", where there is a space
!  between the plus and the number.  So I am going to comment out
!  this line, because I think that's all that's keeping me from
!  doing this.
!
!       if ( ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    if ( ihave == 2 ) then

    else if ( ihave == 6 .or. ihave == 7 ) then
      iterm = 1
    else if ( ihave > 1 ) then
      ihave = 11
    end if
!
!  Comma.
!
  else if ( chrtmp == ',' .or. chrtmp == ';' ) then

    if ( ihave /= 1 ) then
      iterm = 1
      ihave = 12
      lchar = lchar + 1
    end if
!
!  Minus sign.
!
  else if ( chrtmp == '-' ) then

    if ( ihave == 1 ) then
      ihave = 2
      isgn = - 1
    else if ( ihave == 6 ) then
      ihave = 7
      jsgn = - 1
    else
      iterm = 1
    end if
!
!  Plus sign.
!
  else if ( chrtmp == '+' ) then

    if ( ihave == 1 ) then
      ihave = 2
    else if ( ihave == 6 ) then
      ihave = 7
    else
      iterm = 1
    end if
!
!  Decimal point.
!
  else if ( chrtmp == '.' ) then

    if ( ihave < 4 ) then
      ihave = 4
    else if ( ihave >= 6 .and. ihave <= 8 ) then
      ihave = 9
    else
      iterm = 1
    end if
!
!  Exponent marker.
!
  else if ( ch_eqi ( chrtmp, 'E' ) .or. ch_eqi ( chrtmp, 'D' ) ) then

    if ( ihave < 6 ) then
      ihave = 6
    else
      iterm = 1
    end if
!
!  Digit.
!
  else if ( ihave < 11 .and. lge ( chrtmp, '0' ) .and. lle ( chrtmp, '9' ) ) then

    if ( ihave <= 2 ) then
      ihave = 3
    else if ( ihave == 4 ) then
      ihave = 5
    else if ( ihave == 6 .or. ihave == 7 ) then
      ihave = 8
    else if ( ihave == 9 ) then
      ihave = 10
    end if

    call ch_to_digit ( chrtmp, ndig )

    if ( ihave == 3 ) then
      rtop = 10.0 * rtop + real ( ndig )
    else if ( ihave == 5 ) then
      rtop = 10.0 * rtop + real ( ndig )
      rbot = 10.0 * rbot
    else if ( ihave == 8 ) then
      jtop = 10 * jtop + ndig
    else if ( ihave == 10 ) then
      jtop = 10 * jtop + ndig
      jbot = 10 * jbot
    end if
!
!  Anything else is regarded as a terminator.
!
  else
    iterm = 1
  end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
  if ( iterm /= 1 .and. lchar+1 < nchar ) then
    go to 10
  end if
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0
  else

    if ( jbot == 1 ) then
      rexp = 10.0**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0**rexp
    end if

  end if

  rval = isgn * rexp * rtop / rbot

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
