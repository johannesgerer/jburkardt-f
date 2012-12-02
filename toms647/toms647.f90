subroutine infaur ( flag, dimen, atmost )

!*****************************************************************************80
!
!! INFAUR initializes the Faure quasirandom number generator.
!
!  Discussion:
!
!    INFAUR checks whether the user-supplied dimension DIMEN of the
!    quasirandom vectors is acceptable (between 2 and 40).
!
!    Then it calculates an upper summation limit HISUM based on DIMEN 
!    and the user-supplied number ATMOST of quasirandom vectors required.
!
!    Then the routine computes values to be stored in the common block
!    /FAURE/ for use by GOFAUR.
!
!    Thanks to Michael Baudin for pointing out an error in a previous
!    version of this function.
!
!  Modified:
!
!    07 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Bennett Fox.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Henri Faure,
!    Discrepance de suites associees a un systeme de numeration
!    (en dimension s),
!    Acta Arithmetica,
!    Volume XLI, 1982, pages 337-351, especially page 342.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom 
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!  Parameters:
!
!    Output, logical FLAG(2), error flags.
!    FLAG(1) is FALSE if the input value of DIMEN is unacceptable.
!    FLAG(2) is FALSE if the input value of ATMOST is unacceptable.
!
!    Input, integer DIMEN, the spatial dimension.  DIMEN should
!    satisfy: 2 <= DIMEN <= 40.
!
!    Input, integer ATMOST, the maximum number of quasirandom
!    vectors to be computed.
!
!  Global Parameters:
!
!    Global, integer DIG_MAX, a bound for ( the number of digits in the 
!    base QS representation of ATMOST + TESTN ) - 1.
!
!    Global, integer S, the spatial dimension (a copy of DIMEN).
!
!    Global, integer QS, the smallest prime greater than or equal to S.
!
!    Global, integer COEF(0:19,0:19), a table of binomial coefficients.
!
!    Global, integer NEXTN, the index of the next quasirandom vector,
!    initialized here to (QS**4 - 1).
!
!    Global, integer TESTN, initialized to QS**4.
!
!    Global, integer HISUM, initialized to 3.
!
!    Global, real RQS, the value ( 1.0 / QS ).
!
  implicit none

  integer, parameter :: dig_max = 19
  integer, parameter :: dim_max = 40

  integer atmost
  integer coef(0:dig_max,0:dig_max)
  integer dimen
  logical flag(2)
  integer hisum
  integer i
  integer j
  integer nextn
  integer, dimension(dim_max) :: primes = (/ &
     1,  2,  3,  5,  5,  7,  7, 11, 11, 11, &
    11, 13, 13, 17, 17, 17, 17, 19, 19, 23, &
    23, 23, 23, 29, 29, 29, 29, 29, 29, 31, &
    31, 37, 37, 37, 37, 37, 37, 41, 41, 41 /)
  integer qs
  real rqs
  integer s
  integer testn

  common /faure/ s, qs, coef, nextn, testn, hisum, rqs

  save /faure/
!
!  Check the spatial dimension.
!
  s = dimen

  flag(1) = .true.
  flag(2) = .true.

  if ( s < 2 .or. dim_max < s ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INFAUR - Fatal error!'
    write ( *, '(a)' ) '  The spatial dimension S should satisfy:'
    write ( *, '(a,i6)' ) '    2 <= S <= ', dim_max
    write ( *, '(a,i6)' ) '  But this input value is S = ', s
    flag(1) = .false.
    return
  end if

  qs = primes(s)
  testn = qs**4
!
!  Compute log ( ATMOST + TESTN ) in base QS using a ratio of natural logs 
!  to get an upper bound on ( the number of digits in the base QS 
!  representation of ATMOST + TESTN ) - 1.
!
  hisum = nint ( log ( real ( atmost + testn ) ) / log ( real ( qs ) ) )

  if ( dig_max < hisum ) then
    flag(2) = .false.
    return
  end if
!
!  Now find binomial coefficients mod QS in a lower-triangular matrix COEF
!  using the recursions:
!
!    binom(i,j) = binom(i-1,j) + binom(i-1,j-1)
!
!  and
!
!    a = b + c implies mod(a,d) = mod ( mod(b,d) + mod(c,d), d )
!
  coef(0:hisum,0:hisum) = 0

  coef(0:hisum,0) = 1

  do j = 1, hisum
    coef(j,j) = 1
  end do

  do j = 1, hisum
    do i = j + 1, hisum
      coef(i,j) = mod ( coef(i-1,j) + coef(i-1,j-1), qs )
    end do
  end do
!
!  Calculating these coefficients mod QS avoids possible overflow
!  problems with raw binomial coefficients.
!
!  Complete the initialization as described in section 2.
!  NEXTN has 4 digits in base QS, so HISUM is set to 3.
!
  nextn = testn - 1
  hisum = 3
  rqs = 1.0E+00 / real ( qs )

  return
end
subroutine gofaur ( quasi )

!*****************************************************************************80
!
!! GOFAUR generates a new quasirandom Faure vector with each call.
!
!  Discussion:
!
!    This routine implements a method of H. Faure for computing
!    quasirandom numbers.
!
!    The routine INFAUR must be called once for a particular
!    sequence before using GOFAUR.
!
!    All inputs come from INFAUR via labelled common /FAURE/; for 
!    their definitions, see INFAUR.
!
!  Modified:
!
!    17 November 2008
!
!  Author:
!
!    Original FORTRAN77 version by Bennett Fox.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    H Faure,
!    Discrepance de suites associees a un systeme de numeration
!    (en dimension s),
!    Acta Arithmetica,
!    Volume XLI, 1982, pages 337-351, especially page 342.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom 
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!  Parameters:
!
!    Output, real QUASI(DIMEN), the next quasirandom vector.
!
  implicit none

  integer, parameter :: dig_max = 19

  integer coef(0:dig_max,0:dig_max)
  integer hisum
  integer i
  integer j
  integer k
  integer ktemp
  integer ltemp
  integer mtemp
  integer nextn
  integer qs
  real quasi(*)
  real r
  real rqs
  integer s
  integer testn
  integer ytemp(0:dig_max)
  integer ztemp

  common /faure/ s, qs, coef, nextn, testn, hisum, rqs

  save /faure/
!
!  Find QUASI(1) using the method of Faure (section 3.3).
!
!  NEXTN has a representation in base QS of the form: 
!
!    Sum ( 0 <= J <= HISUM ) YTEMP(J) * QS**J
!
!  We now compute the YTEMP(J)'s.
!
  ktemp = testn
  ltemp = nextn
  do i = hisum, 0, -1
    ktemp = ktemp / qs
    mtemp = mod ( ltemp, ktemp )
    ytemp(i) = ( ltemp - mtemp ) / ktemp
    ltemp = mtemp
  end do
!
!  QUASI(K) has the form
!
!    Sum ( 0 <= J <= HISUM ) YTEMP(J) / QS**(J+1)
!
!  Ready to compute QUASI(1) using nested multiplication.
!
  r = real ( ytemp(hisum) )
  do i = hisum - 1, 0, -1
    r = real ( ytemp(i) ) + rqs * r
  end do

  quasi(1) = r * rqs
!
!  Find components QUASI(2:S) using the Faure method 
!  (sections 3.2 and 3.3).
!
  do k = 2, s

    quasi(k) = 0.0E+00
    r = rqs

    do j = 0, hisum

      ztemp = dot_product ( ytemp(j:hisum), coef(j:hisum,j) )
!
!  No apparent alternative one-dimensional coefficient array
!  except via subscript address computations and equivalencing.
!
!  New YTEMP(J) is:
!
!    Sum ( J <= I <= HISUM ) ( old ytemp(i) * binom(i,j) ) mod QS.
!
      ytemp(j) = mod ( ztemp, qs )
      quasi(k) = quasi(k) + real ( ytemp(j) ) * r
      r = r * rqs

    end do

  end do
!
!  Update NEXTN and, if needed, TESTN and HISUM.
!  HISUM is guaranteed to be no greater than DIG_MAX.
!
  nextn = nextn + 1

  if ( nextn == testn ) then
    testn = testn * qs
    hisum = hisum + 1
  end if

  return
end
subroutine inhalt ( flag, dimen, atmost, quasi )

!*****************************************************************************80
!
!! INHALT initializes the Halton quasirandom number generator.
!
!  Discussion:
!
!    INHALT first checks whether the user-supplied dimension DIMEN of
!    the quasirandom vectors is acceptable (between 2 and 40).
!
!    INHALT then calculates a tolerance parameter E to make the program work
!    correctly in finite precision arithmetic and a parameter DELTA
!    to check that E works.  If the test is not passed, then ATMOST
!    is too big relative to the machine precision.
!
!    Otherwise, INHALT computes and returns the first vector QUASI.
!    For the following values of QUASI, it is necessary to call GOHALT.
!
!  Modified:
!
!    18 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Bennett Fox.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom 
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    J H Halton and G B Smith,
!    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
!    Communications of the ACM,
!    Volume 7, 1964, pages 701-702.
!
!  Parameters:
!
!    Output, logical FLAG(2), error flags.
!    FLAG(1) is FALSE if the input value of DIMEN is unacceptable.
!    FLAG(2) is FALSE if the input value of ATMOST is unacceptable.
!
!    Input, integer DIMEN, the spatial dimension.  DIMEN should
!    satisfy: 2 <= DIMEN <= 40.
!
!    Input, integer ATMOST, the maximum number of quasirandom
!    vectors to be computed.
!
!    Output, real ( kind = 8 ) QUASI(DIMEN), the first element of
!    the Halton sequence.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) PRIME(40), the first 40 primes.
!
!  Global Parameters:
!
!    Stored in common block /HALTON/:
!
!    Global, real ( kind = 8 ) E, a tolerance.
!
!    Global, real ( kind = 8 ) PRIME_INV(40), the reciprocals of the
!    first 40 primes.
!
!    Global, integer S, the spatial dimension.
!
  implicit none

  integer dimen
  integer, parameter :: dim_max = 40

  integer atmost
  real ( kind = 8 ) delta
  real ( kind = 8 ) e
  logical flag(2)
  integer prime(dim_max)
  real ( kind = 8 ) prime_inv(dim_max)
  real ( kind = 8 ) quasi(dimen)
  integer s
  real ( kind = 8 ) small

  common /halton/ e, prime_inv, s

  save /halton/
!
!  Check DIMEN.
!
  flag(1) = .true.
  flag(2) = .true.

  s = dimen

  if ( s < 2 .or. dim_max < s ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INHALT - Fatal error!'
    write ( *, '(a)' ) '  The spatial dimension S should satisfy:'
    write ( *, '(a,i6)' ) '    2 <= S <= ', dim_max
    write ( *, '(a,i6)' ) '  But this input value is S = ', s
    flag(1) = .false.
    return
  end if
!
!  Set the primes.
!
  prime(1:dim_max) = (/ &
      2,   3,   5,   7,  11,  13,  17,  19,  23,  29, &
     31,  37,  41,  43,  47,  53,  59,  61,  67,  71, &
     73,  79,  83,  89,  97, 101, 103, 107, 109, 113, &
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173 /)
!
!  Compute the tolerance and make the check.
!
  small = epsilon ( small )

  e = 0.9D+00 * ( 1.0D+00 / ( dble ( atmost * prime(s) ) ) - 10.0D+00 * small )

  delta = 100.0D+00 * small * dble ( atmost + 1 ) * log10 ( dble ( atmost ) )

  if ( 0.09D+00 * ( e - 10.0D+00 * small ) < delta ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INHALT - Fatal error!'
    write ( *, '(a)' ) '  The value of ATMOST is too great.'
    flag(2) = .false.
    return
  end if
!
!  Set the inverse primes.
!
  prime_inv(1:dim_max) = 1.0D+00 / dble ( prime(1:dim_max) )
!
!  Compute the first vector.
!
  quasi(1:s) = prime_inv(1:s)

  return
end
subroutine gohalt ( quasi )

!*****************************************************************************80
!
!! GOHALT generates a new quasirandom Halton vector with each call.
!
!  Discussion:
!
!    The routine adapts key ideas from Halton and Smith.
!
!    The routine INHALT must be called once before using
!    this routine.
!
!  Modified:
!
!    17 November 2008
!
!  Author:
!
!    Original FORTRAN77 version by Bennett Fox.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom 
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    J H Halton and G B Smith,
!    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
!    Communications of the ACM,
!    Volume 7, 1964, pages 701-702.
!
!  Parameters:
!
!    Input/output, real QUASI(DIMEN), on input, the previous
!    quasirandom vector; on output, the next quasirandom vector.
!    On the first call, the input value should be the output
!    value given by INHALT.
!
!  Global Parameters:
!
!    In labelled common block /HALTON/:
!
!    Global, real ( kind = 8 ) E, a tolerance.
!
!    Global, real ( kind = 8 ) PRIME_INV(40), the reciprocal of
!    the first 40 prime numbers.
!
!    Global, integer S, the spatial dimension.
!
  implicit none

  integer, parameter :: dim_max = 40

  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer i
  real ( kind = 8 ) prime_inv(dim_max)
  real ( kind = 8 ) quasi(*)
  integer s
  real ( kind = 8 ) t

  common /halton/ e, prime_inv, s

  save /halton/
!
!  Generate QUASI one component at a time, using radix 1/PRIME(K) for 
!  component K.
!
  do i = 1, s

    t = prime_inv(i)
    f = 1.0D+00 - quasi(i)
    g = 1.0D+00
    h = prime_inv(i)

    do

      if ( e <= f - h ) then
        exit
      end if
!
!  This checks whether Q + H > 1 - E.
!
      g = h
      h = h * t
!
!  If this is the (K-1)-st time this statement is reached, check whether
!  QUASI(I) + R**(-K) > 1-E.
!
    end do
!
!  For the appropriate I (depending on how many times the loop above
!  is executed), add H**(I+1) + H**(I) - 1
!  to the old QUASI(I) to get the next QUASI(I).
!
    quasi(i) = g + h - f

  end do

  return
end
subroutine insobl ( flag, dimen, atmost, taus )

!*****************************************************************************80
!
!! INSOBL initializes the Sobol quasirandom number generator.
!
!  Discussion:
!
!    INSOBL first checks whether the user-supplied dimension DIMEN 
!    of the quasi-random vectors is strictly between 0 and 41.
!
!    Next it checks the size of ATMOST, an upper bound on the number
!    of calls the user intends to make on GOSOBL.  
!
!    The leading elements of each row of V were initialized by BDSOBL.
!    Each row corresponds to a primitive polynomial.  If the polynomial has
!    degree M, elements after the first M are now calculated.
!    The numbers in V are actually binary fractions.
!
!    INSOBL implicitly computes the first Sobol vector, but since it
!    is all zero, it does not return it to the calling program.  
!    Subsequent vectors come from GOSOBL.
!
!    The array POLY gives successive primitive polynomials coded in 
!    binary, e.g.
!
!      45 = 100101
!
!    has bits 5, 2, and 0 set (counting from the right) and therefore 
!    represents the polynomial
!
!      x**5 + x**2 + x**0.
!
!    These polynomials are listed in the order used by Sobol. 
!
!    A more complete table is given in Sobol and Levitan.
!
!    The initialization of the array V is from the
!    latter paper.  For a polynomial of degree M, M initial
!    values are needed.  These are the values given in the initialization
!    statements.  Subsequent values of V are calculated.
!
!  Modified:
!
!    17 November 2008
!
!  Author:
!
!    Original FORTRAN77 version by Bennett Fox.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom 
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Paul Bratley and Bennett Fox,
!    Algorithm 659:
!    Implementing Sobol's Quasirandom Sequence Generator,
!    ACM Transactions on Mathematical Software,
!    Volume 14, Number 1, pages 88-100, 1988.
!
!    I Sobol,
!    USSR Computational Mathematics and Mathematical Physics,
!    Volume 16, pages 236-242, 1977.
!
!    I Sobol and Levitan, 
!    The Production of Points Uniformly Distributed in a Multidimensional 
!    Cube (in Russian),
!    Preprint IPM Akad. Nauk SSSR, 
!    Number 40, Moscow 1976.
!
!  Parameters:
!
!    Output, logical FLAG(2), error flags.
!    FLAG(1) is FALSE if the input value of DIMEN is unacceptable.
!    FLAG(2) is FALSE if the input value of ATMOST is unacceptable.
!
!    Input, integer DIMEN, the spatial dimension.  DIMEN should
!    satisfy: 2 <= DIMEN <= 40.
!
!    Input, integer ATMOST, the maximum number of quasirandom
!    vectors to be computed.  ATMOST should be positive, and less
!    than 2**30.
!
!    Output, integer TAUS, is for determining "favorable" values.  As
!    discussed in Bratley and Fox, these have the form N = 2**K 
!    where (TAUS+DIMEN-1) <= K for integration and TAUS < K for global 
!    optimization.  Useful values of TAUS are available for
!    2 <= DIMEN <= 13.
!
!  Global Parameters:
!
!    In common block /SOBOL/:
!
!    Global, integer V(40,30), table of direction numbers.
!
!    Global, integer S, the spatial dimension.
!
!    Global, integer MAXCOL, the last column of V to be used.
!
!    Global, integer COUNT, the sequence number of this call.
!
!    Global, integer LASTQ(40), the numerators for the last vector generated.
!
!    Global, real RECIPD, (1/denominator) for these numerators.
!
  implicit none

  integer dimen
  integer, parameter :: dim_max = 40

  integer atmost
  integer count
  logical flag(2)
  integer i
  logical includ(8)
  integer j
  integer j2
  integer k
  integer l
  integer lastq(40)
  integer m
  integer maxcol
  integer newv
  integer poly(40)
  real recipd
  integer s
  integer, save, dimension ( 13 ) :: tau = (/ &
    0, 0, 1, 3, 5, 8, 11, 15, 19, 23, 27, 31, 35 /)
  integer taus
  integer, dimension ( 40, 30 ) :: v

  common /sobol/ v, s, maxcol, count, poly, lastq, recipd

  save /sobol/
!
!  Initialize (part of) V.
!
  v(1:40,1) = (/ &
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)

  v(3:40,2) = (/ &
          1, 3, 1, 3, 1, 3, 3, 1, &
    3, 1, 3, 1, 3, 1, 1, 3, 1, 3, &
    1, 3, 1, 3, 3, 1, 3, 1, 3, 1, &
    3, 1, 1, 3, 1, 3, 1, 3, 1, 3 /)

  v(4:40,3) = (/ &
             7, 5, 1, 3, 3, 7, 5, & 
    5, 7, 7, 1, 3, 3, 7, 5, 1, 1, &
    5, 3, 3, 1, 7, 5, 1, 3, 3, 7, &
    5, 1, 1, 5, 7, 7, 5, 1, 3, 3 /)

  v(6:40,4) = (/ &
                   1, 7, 9,13,11, &
    1, 3, 7, 9, 5,13,13,11, 3,15, &
    5, 3,15, 7, 9,13, 9, 1,11, 7, &
    5,15, 1,15,11, 5, 3, 1, 7, 9 /)

  v(8:40,5) = (/ &
                         9, 3,27, &
   15,29,21,23,19,11,25, 7,13,17, &
    1,25,29, 3,31,11, 5,23,27,19, &
   21, 5, 1,17,13, 7,15, 9,31, 9 /)

  v(14:40,6) = (/ &
            37,33, 7, 5,11,39,63, &
   27,17,15,23,29, 3,21,13,31,25, &
    9,49,33,19,29,11,19,27,15,25 /)

  v(20:40,7) = (/ &
                                       13, &
   33,115, 41, 79, 17, 29,119, 75, 73,105, &
    7, 59, 65, 21,  3,113, 61, 89, 45,107 /)

  v(38:40,8) = (/ &
                                7, 23, 39 /)
!
!  Set POLY.
!
  poly(1:40)= (/ &
      1,   3,   7,  11,  13,  19,  25,  37,  59,  47, &
     61,  55,  41,  67,  97,  91, 109, 103, 115, 131, &
    193, 137, 145, 143, 241, 157, 185, 167, 229, 171, &
    213, 191, 253, 203, 211, 239, 247, 285, 369, 299 /)
!
!  Check parameters.
!
  flag(1) = .true.
  flag(2) = .true.

  s = dimen

  if ( s < 2 .or. dim_max < s ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INSOBL - Fatal error!'
    write ( *, '(a)' ) '  The spatial dimension S should satisfy:'
    write ( *, '(a,i6)' ) '    2 <= S <= ', dim_max
    write ( *, '(a,i6)' ) '  But this input value is S = ', s
    flag(1) = .false.
    return
  end if

  if ( atmost < 1 .or. 2**30 <= atmost ) then
    flag(2) = .false.
    return
  end if

  if ( s <= 13 ) then
    taus = tau(s)
  else
    taus = -1
  end if
!
!  Find the number of bits in ATMOST.
!
  i = atmost
  maxcol = 0

  do

    maxcol = maxcol + 1
    i = i / 2
    if ( i <= 0 ) then
      exit
    end if

  end do
!
!  Initialize row 1 of V.
!
  v(1,1:maxcol) = 1
!
!  Initialize the remaining rows of V.
!
  do i = 2, s
!
!  The bit pattern of the integer POLY(I) gives the form
!  of polynomial I.
!
!  Find the degree of polynomial I from binary encoding.
!
    j = poly(i)
    m = 0

    do

      j = j / 2

      if ( j <= 0 ) then
        exit
      end if

      m = m + 1

    end do
!
!  We expand this bit pattern to separate components
!  of the logical array INCLUD.
!
    j = poly(i)
    do k = m, 1, -1
      j2 = j / 2
      includ(k) = ( j /= 2 * j2 )
      j = j2
    end do
!
!  Calculate the remaining elements of row I as explained
!  in Bratley and Fox, section 2.
!
    do j = m + 1, maxcol

      newv = v(i,j-m)
      l = 1

      do k = 1, m

        l = 2 * l

        if ( includ(k) ) then
          newv = ieor ( newv, l * v(i,j-k) )
        end if

      end do

      v(i,j) = newv

    end do

  end do
!
!  Multiply columns of V by appropriate power of 2.
!
  l = 1
  do j = maxcol - 1, 1, -1
    l = 2 * l
    v(1:s,j) = v(1:s,j) * l
  end do
!
!  RECIPD is 1/(common denominator of the elements in V).
!
  recipd = 1.0E+00 / real ( 2 * l )
!
!  Set up first vector and values for GOSOBL.
!
  count = 0
  lastq(1:s) = 0

  return
end
subroutine gosobl ( quasi )

!*****************************************************************************80
!
!! GOSOBL generates a new quasirandom Sobol vector with each call.
!
!  Discussion:
!
!    The routine adapts the ideas of Antonov and Saleev.
!
!    The routine INSOBL must be called once, for a particular
!    set of parameters, before calling GOSOBL.
!
!    GOSOBL checks that the user does not make more calls than
!    the value ATMOST that was specified in the call to INSOBL.
!
!  Modified:
!
!    17 November 2008
!
!  Author:
!
!    Original FORTRAN77 version by Bennett Fox.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Antonov and Saleev,
!    USSR Computational Mathematics and Mathematical Physics,
!    Volume 19, 1980, pages 252 - 256.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom 
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!  Parameters:
!
!    Output, real QUASI(DIMEN), the next quasirandom vector.
!
!  Global Parameters:
!
!    In labelled common /SOBOL/:
!
!    Global, integer V(40,30), table of direction numbers.
!
!    Global, integer S, the spatial dimension.
!
!    Global, integer MAXCOL, the last column of V to be used.
!
!    Global, integer COUNT, the sequence number of this call.
!
!    Global, integer LASTQ(40), the numerators for the last vector generated.
!
!    Global, real RECIPD, (1/denominator) for these numerators.
!
  implicit none

  integer count
  integer i
  integer i2
  integer l
  integer lastq(40)
  integer maxcol
  integer poly(40)
  real quasi(*)
  real recipd
  integer s
  integer v(40,30)

  common /sobol/ v, s, maxcol, count, poly, lastq, recipd

  save /sobol/
!
!  Find the position of the right-hand zero in COUNT.
!
  l = 0
  i = count

  do

    l = l + 1
    i2 = i / 2

    if ( i == 2 * i2 ) then
      exit
    end if

    i = i2

  end do
!
!  Check that the user is not calling too many times!
!
  if ( maxcol < l ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GOSOBL - Fatal error!'
    write ( *, '(a)' ) '  Too many calls!'
    write ( *, '(a,i12)' ) '  MAXCOL = ', maxcol
    write ( *, '(a,i12)' ) '  L =      ', l
    stop
  end if
!
!  Calculate the new components of QUASI.
!  RECIPD is used to normalize the values.
!
  do i = 1, s

     lastq(i) = ieor ( lastq(i), v(i,l) )

     quasi(i) = real ( lastq(i) ) * recipd

  end do

  count = count + 1

  return
end
function exor ( iin, jin )

!*****************************************************************************80
!
!! EXOR calculates the exclusive OR of two integers.
!
!  Discussion:
!
!    The FORTRAN90 intrinsic IEOR ( IIN, JIN ) should be used
!    instead of this routine!
!
!  Modified:
!
!    17 November 2008
!
!  Author:
!
!    Original FORTRAN77 version by Bennett Fox.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer IIN, JIN, two values whose exclusive OR is needed.
!
!    Output, integer EXOR, the exclusive OR of IIN and JIN.
!
  implicit none

  integer exor
  integer i
  integer i2
  integer iin
  integer j
  integer j2
  integer jin
  integer k
  integer l

  i = iin
  j = jin
  k = 0
  l = 1

  do

    if ( i == 0 .and. j == 0 ) then
      exit
    end if
!
!  Check the current right-hand bits of i and j.
!  If they differ, set the appropriate bit of k.
!
    i2 = i / 2
    j2 = j / 2

    if ( ( i == 2 * i2 ) .neqv. ( j == 2 * j2 ) ) then
      k = k + l
    end if

    i = i2
    j = j2
    l = 2 * l

  end do

  exor = k

  return
end
function unif ( seed )

!*****************************************************************************80
!
!! UNIF is a portable pseudorandom number generator.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      unif = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!  Modified:
!
!    17 November 2008
!
!  Author:
!
!    Original FORTRAN77 version by Bennett Fox.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, L E Schrage,
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
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which must
!    be greater than 0 and less than 2**31 - 1.  On output,
!    SEED has been updated.
!
!    Output, real UNIF, a new pseudorandom variate, strictly between
!    0 and 1.
!
  implicit none

  integer k1
  integer seed
  real unif

  k1 = seed / 127773

  seed = 16807 * ( seed - k1 * 127773 ) - k1 * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  unif = real ( seed ) * 4.656612875e-10

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
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
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
