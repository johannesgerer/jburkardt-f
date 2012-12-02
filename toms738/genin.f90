program main

!*****************************************************************************80
!
!! MAIN is the main program for GENIN.
!
!  Discussion:
!
!    GENIN is the driver for the interactive Niederreiter sequence program.
!
!    This program tests the accuracy of numerical integration
!    using the low-discrepancy binary sequences of
!    Harald Niederreiter (1988) as implemented in INLO, GOLO, and
!    related programs.  Various possible test integrals are
!    provided by the function TESTF.
!
!    For a prime-power base, arithmetic tables must be read from
!    the file "gftabs.txt".  For any base, a set of irreducible
!    polynomials is read from "gfplys.txt".
!
!    Both the general-base and base-2 programs assume that
!    your computer's word length is 31 bits, excluding sign.
!    If this is not the case, modify the parameter NBITS
!    throughout the PARAMETER statements in this set of
!    programs accordingly.
!
!
!    There are theoretical reasons to believe that BASE ** E,
!    where E is defined for example in Bratley, Fox, and
!    Niederreiter (1991), would be a good choice.  
!
!    However, we don't want to come anywhere near the largest possible
!    machine-representable integer; hence, the compromise
!    exponents above.  
!
!    Note: Subject to this conditon, it can't hurt to take an exponent 
!    greater than E, because warm-up skipping of initial values is done 
!    implicitly in O(1) time.  The maximum value of E for a fixed dimension
!    S grows like log(S).  We allow some "fat" for the implicit
!    constant in our choice of POWER.
!
!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederrieter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, pages 494-495, 1994.
!
!    Harald Niederreiter,
!    Low-discrepancy and low-dispersion sequences,
!    Journal of Number Theory,
!    Volume 30, 1988, pages 51-70.
!
!  Local Parameters:
!
!    Local, integer MAXBAS, the maximum asymptotically-optimal base 
!    used up to MAXDIM.
!
!    Local, integer MAXDIM, the maximum spatial dimension that will be used.
!
!    Local, integer OPTBAS(2:MAXDIM), gives the asymptotically-optimal
!    base for the respective dimension.
!
!    Local, integer OUNIT, the unit to save the output.
!
!    Local, integer POWER(2:MAXBAS), values used in a possible
!    warm-up calculation.
!
  implicit none

  integer, parameter :: maxbas = 13
  integer, parameter :: maxdim = 20

  integer base
  real ( kind = 8 ) correct
  integer dimen
  real ( kind = 8 ) error
  real ( kind = 8 ) estimate
  real ( kind = 8 ) exactf
  character ( len = 80 ) :: file_out_name = 'genin_prb_output.txt'
  integer file_out_unit
  integer i
  integer i4_characteristic
  real ( kind = 8 ) integral_est
  integer num
  integer, dimension(2:12) :: optbas = (/ &
    2, 3, 3, 5, 7, 7, 9, 9, 11, 11, 13 /)
  integer pbase
  integer, dimension(2:maxbas) :: power = (/ &
    12, 8, 8, 6, 6, 6, 4, 4, 4, 4, 4, 4 /)
  real ( kind = 8 ) quasi(maxdim)
  integer seqlen
  integer skip
  integer step
  real ( kind = 8 ) testf

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GENIN'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program demonstrates the Niederreiter'
  write ( *, '(a)' ) '  quasirandom sequence by applying it to the problem'
  write ( *, '(a)' ) '  of approximating a multidimensional integral.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) '  Output will be stored in "' // trim ( file_out_name ) // '".'

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace' )

  do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Choose a test integral (1 to 4) or 0 to quit :'
    read ( *, * ) num

    if ( num <= 0 ) then
      close ( unit = file_out_unit )
      exit
    end if

    if ( 4 < num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GENIN - Error!'
      write ( *, '(a)' ) '  No such test integral.'
      cycle
    end if
!
!  Each test integral is parameterized by its dimension.
!
    do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Enter the spatial dimension, '
      write ( *, '(a,i8)' ) '  between 1 and ', maxdim

      read ( *, * ) dimen

      if ( 1 <= dimen .and. dimen <= maxdim ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GENIN - Error!'
      write ( *, '(a,i8)' ) '  The dimension must be between 1 and ', maxdim
 
    end do

    do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Choose a prime or prime-power base.'
      if ( dimen <= 12 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The asymptotically-optimal base '
        write ( *, '(a,i8)' ) '  for this dimension is ', optbas(dimen)
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (This base might not be empirically optimal.)'
        write ( *, '(a)' ) ' '
      end if
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Enter base: '

      read ( *, * ) base

      if ( i4_characteristic ( base ) /= 0 ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GENIN - Error!'
      write ( *, '(a)' ) '  Out of range or bad value:  try again'

    end do
!
!  The sequence length is the number of quasi-random points used to 
!  estimate the integral, excluding warm-up.  The number of initial 
!  quasi-random points deleted during warm-up is given by SKIP,
!  chosen below.
!
!  Some users may wish to rewrite the driver to test a [heuristic] 
!  "convergence" criterion, stopping the generation of points when it 
!  is passed or when a specified number of points have been generated,
!  whichever occurs first.
!
    do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Choose sequence length :'
      write ( *, '(a)' ) '  A power of the base is recommended, such as:'
      write ( *, '(a)' ) ' '

      write ( *, '(i12)' ) base ** power(base)
      write ( *, '(i12)' ) base ** ((power(base) + 1))
      write ( *, '(i12)' ) base ** ((power(base) + 2))
      write ( *, '(i12)' ) base ** ((power(base) + 3))

      write ( *, '(a)' ) ' '
      write ( *, '(a)') '  Enter your choice for the sequence length:'

      read ( *, * ) seqlen

      if ( 0 < seqlen ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GENIN - Error!'
      write ( *, '(a)' ) '  The sequence length must be strictly positive.'

    end do

    do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Choose SKIP, the number of initially computed '
      write ( *, '(a)' ) '  sequence values to skip.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A heuristic formula is'
      write ( *, '(a)' ) '    SKIP = BASE ** power(BASE)'
      write ( *, '(a)' ) '  when BASE <= MAXBAS, otherwise SKIP = 10000.'

      if ( base <= maxbas ) then
        skip = base ** power(base)
      else
        skip = 10000
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Numerically, this heuristic value is ', skip
      write ( *, '(a)' ) '  Enter SKIP (not necessarily the value above):'            
      read ( *, * ) skip

      if ( 0 <= skip ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GENIN - Error!'
      write ( *, '(a)' ) '  The value of SKIP must be >= 0.'
 
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENIN:'
    write ( *, '(a)' ) '  Call INLO for initialization.'

    call inlo ( dimen, base, skip )
!
!  Write title and the exact value of the integral.
!
    correct = exactf ( num, dimen )

    write ( file_out_unit, '(a)' ) ' '
    write ( file_out_unit, '(a,i8)' ) 'Test integral:     ', num
    write ( file_out_unit, '(a)' ) ' '
    write ( file_out_unit, '(a,i8)' ) 'Spatial dimension: ', dimen
    write ( file_out_unit, '(a,i8)' ) 'Base:              ', base
    write ( file_out_unit, '(a,i8)' ) 'Sequence length:   ', seqlen
    write ( file_out_unit, '(a,i8)' ) 'Skipped values:    ', skip
    write ( file_out_unit, '(a)' ) ' '
    write ( file_out_unit, '(a,g14.6)' ) 'Correct value is ', correct
    write ( file_out_unit, '(a)' ) ' '
    write ( file_out_unit, '(a)' ) &
      '      Iteration     Estimate of integral       Error'
    write ( file_out_unit, '(a)' ) ' '
!
!  Now estimate the integral
!
    pbase = base
    integral_est = 0.0D+00
    step = 500

    do i = 1, seqlen

      call golo ( quasi )

      integral_est = integral_est + testf ( num, dimen, quasi )

      if ( mod ( i, step ) == 0 ) then
        estimate = integral_est / real ( i, kind = 8 )
        error = abs ( estimate - correct )
        write ( file_out_unit, '(i12,2g24.7)' ) i, estimate, error
      end if
!
!  This finds the next power of the base.  There is reason to believe
!  that convergence properties of the sequence of estimates is
!  better along the subsequence corrsponding to powers of the base.
!
      if ( mod ( i, pbase ) == 0 )  then
        estimate = integral_est / real ( i, kind = 8 )
        error = abs ( estimate - correct )
        write ( file_out_unit, '(i12,2g24.7)' ) i, estimate, error
        pbase = pbase * base
      end if

      if ( i == 5000 ) then
        step = 1000
      else if ( i == 10000 ) then
        step = 5000
      end if

    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENIN:'
    write ( *, '(a)' ) '  End of iteration for this problem.'

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GENIN:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine calcc ( )

!*****************************************************************************80
!
!! CALCC calculates the value of the constants C(I,J,R).
!
!  Discussion:
!
!    This routine calculates the values of the constants C(I,J,R).
!    As far as possible, we use Niederreiter's notation.
!    We calculate the values of C for each I in turn.
!    When all the values of C have been calculated, we return
!    this array to the calling program.
!
!    Irreducible polynomials are read from file "gfplys.txt"
!    This file should have been created earlier by running the
!    GFPLYS program.
!
!    Polynomials stored as arrays have the coefficient of degree n 
!    in POLY(N), and the degree of the polynomial in POLY(-1).  
!    The parameter DEG is just to remind us of this last fact.  
!    A polynomial which is identically 0 is given degree -1.
!
!    Thanks to Michael Baudin for pointing out that MAXE should
!    be increased from 5 to 7, since one of the irreducible polynomials
!    that must be stored in PX has degree 7, 07 June 2010.
!
!  Modified:
!
!    07 June 2010
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederrieter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, pages 494-495, 1994.
!
!    Harald Niederreiter,
!    Low-discrepancy and low-dispersion sequences,
!    Journal of Number Theory,
!    Volume 30, 1988, pages 51-70.
!
!  Local Parameters:
!
!    Local, integer MAXDEG, the highest degree of polynomial
!    to be handled.  
!
!    Local, integer MAXDIM, the maximum dimension that will be used.
!
!    Local, integer MAXE; we need MAXDIM irreducible polynomials over GF(Q).
!    MAXE is the highest degree among these.
!
!    Local, integer MAXFIG, the maximum number of base Q digits we can handle.
!
!    Local, integer MAXINT, the largest fixed point integer we can represent.
!
!    Local, integer MAXQ, the order of the largest field to be handled.
!
!    Local, integer MAXV, the maximum index used in V.
!
!    Local, integer NBITS, the number of bits in a fixed-point integer, not
!    counting the sign.
!
!    Local, integer NPOLS, the number of precalculated irreducible polynomials.
!
!  Global parameters in /COMM/:
!
!    Global, integer C(MAXDIM,MAXFIG,0:MAXFIG-1), the values of 
!    Niederreiter's C(I,J,R).
!
!    Global, integer COUNT(0:MAXFIG-1), the index of the current item 
!    in the sequence, expressed as an array of base-Q digits.  COUNT(R)
!    is the same as Niederreiter's AR(N) (page 54) except that N is implicit.
!
!    Global, integer D(MAXDIM,MAXFIG).
!
!    Global, integer NEXTQ(MAXDIM), the numerators of the next item in 
!    the series.  These are like Niederreiter's XI(N) (page 54) except that
!    N is implicit, and the NEXTQ are integers.  To obtain the values of 
!    XI(N), multiply by RECIP.
!
!    Global, integer QPOW(MAXFIG), to speed things up a bit. 
!    QPOW(I) = Q ** (NFIGS-I).
!
!    Global, integer DIMEN, the dimension of the sequence to be generated.
!
!    Global, integer NFIGS, the number of base Q digits we are using.
!
!    Global, real ( kind = 8 ) RECIP = 1.0 / Q**NFIGS.
!
  implicit none

  integer, parameter :: deg = -1
  integer, parameter :: maxdeg = 50
  integer, parameter :: maxdim = 20
! integer, parameter :: maxe = 5
  integer, parameter :: maxe = 7
  integer, parameter :: maxfig = 20
  integer, parameter :: maxq = 50
  integer, parameter :: nbits = 31
  integer, parameter :: npols = 25

  integer, parameter :: maxv = maxfig + maxe

  integer add(0:maxq-1,0:maxq-1)
  integer b(-1:maxdeg)
  integer c(maxdim,maxfig,0:maxfig-1)
  integer count(0:maxfig-1)
  integer d(maxdim,maxfig)
  integer dimen
  integer e
  character ( len = 80 ) :: file_in_name = 'gfplys.txt'
  integer file_in_unit
  integer i
  integer ios
  integer j
  integer mul(0:maxq-1,0:maxq-1)
  integer nextq(maxdim)
  integer nfigs
  integer p
  integer px(-1:maxe)
  integer q
  integer qpow(maxfig)
  integer r
  real ( kind = 8 ) recip
  integer sub(0:maxq-1,0:maxq-1)
  integer u
  integer v(0:maxv)

  common /comm/ c, count, d, nextq, qpow, dimen, nfigs, recip
  common /field/ p, q, add, mul, sub

  save /field/
  save /comm/
!
!  Read the irreducible polynomials.
!
  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old' )

  do

    read ( file_in_unit, '(20i3)', iostat = ios ) i

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CALCC - Fatal error!' 
      write ( *, '(a,i8)' ) '  Could not find tables for Q = ', q
      stop
    end if

    if ( i == q ) then
      exit
    end if

    do j = 1, npols
      read ( file_in_unit, '(20i3)' )
    end do

  end do

  do i = 1, dimen
!
!  For each dimension, we need to calculate powers of an
!  appropriate irreducible polynomial.  See Niederreiter
!  page 65, just below equation (19).
!
!  Read the appropriate irreducible polynomial into PX,
!  and its degree into E.  Set polynomial B = PX ** 0 = 1.
!  M is the degree of B.  Subsequently B will hold higher
!  powers of PX.
!
!  The polynomial PX is stored in 'gfplys.txt' in the format
!
!    n  a0  a1  a2  ... an
!
!  where n is the degree of the polynomial and the ai are
!  its coefficients.
!
    read ( file_in_unit, '(20i3)' ) e, px(0:e)

    px(deg) = e
    b(deg) = 0
    b(0) = 1
!
!  Niederreiter (page 56, after equation (7), defines two variables 
!  Q and U.  We do not need Q explicitly, but we do need U.
!
    u = 0

    do j = 1, nfigs
!
!  If U = 0, we need to set B to the next power of PX
!  and recalculate V.  This is done by subroutine CALCV.
!
      if ( u == 0 ) then
        call calcv ( px, b, v, maxv )
      end if
!
!  Now C is obtained from V.  Neiderreiter obtains A from V 
!  (page 65, near the bottom), and then gets C from A (page 56,
!  equation (7)).  However this can be done in one step.
!
      do r = 0, nfigs - 1
        c(i,j,r) = v(r+u)
      end do
!
!  Increment U.  If U = E, then U = 0 and in Niederreiter's
!  paper Q = Q + 1.  Here, however, Q is not used explicitly.
!
      u = u + 1
      if ( u == e ) then
        u = 0
      end if

    end do

  end do

  close ( unit = file_in_unit )

  return
end
subroutine calcv ( px, b, v, maxv )

!*****************************************************************************80
!
!! CALCV calculates the constants V(J,R).
!
!  Discussion:
!
!    This program calculates the values of the constants V(J,R) as
!    described in Bratley, Fox and Niederreiter, section 3.3.  It 
!    is called from either CALCC or CALCC2.  The values transmitted 
!    through common /FIELD/ determine which field we are working in.
!
!    Polynomials stored as arrays have the coefficient of degree n 
!    in POLY(N), and the degree of the polynomial in POLY(-1).  The 
!    parameter DEG is just to remind us of this last fact.  A polynomial 
!    which is identically 0 is given degree -1.
!
!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederrieter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, pages 494-495, 1994.
!
!  Parameters:
!
!    Input, integer PX(-1:MAXDEG), the appropriate irreducible polynomial 
!    for the dimension currently being considered.  The degree of PX will 
!    be called E.
!
!    Input/output, integer B(-1:MAXDEG).  On input, B is the polynomial 
!    defined in section 2.3 of BFN.  The degree of B implicitly defines 
!    the parameter J of section 3.3, by degree(B) = E*(J-1).  On output, 
!    B has been multiplied by PX, so its degree is now E*J.
!
!    Input, integer V(0:MAXV), contains the values required.
!
!    Input, integer MAXV, the dimension of the array V.
!
!  Local Parameters:
!
!    Local, integer ARBIT(), indicates where the user can place
!    an arbitrary element of the field of order Q.  For the code,
!    this means 0 <= ARBIT < Q.  Within these limits, the user can 
!    do what he likes.  ARBIT could be declared as a function 
!    that returned a different arbitrary value each time it is referenced.
!
!    Local, integer BIGM, is the M used in section 3.3.  It differs from 
!    the [little] m used in section 2.3, denoted here by M.
!
!    Local, integer MAXDEG, the highest degree of polynomial to be
!    handled. 
! 
!    Local, integer MAXQ, the order of the largest field to be handled.
!
!    Local, integer NONZER shows where the user must put an arbitrary 
!    non-zero element of the same field.  For the code, this means 
!    0 < NONZER < Q.  Within these limits, the user can do what he likes.  
!
  implicit none

  integer, parameter :: deg = -1
  integer, parameter :: maxdeg = 50
  integer, parameter :: maxq = 50

  integer add(0:maxq-1,0:maxq-1)
  integer, parameter :: arbit = 1
  integer b(-1:maxdeg)
  integer bigm
  integer e
  integer h(-1:maxdeg)
  integer i
  integer j
  integer kj
  integer m
  integer maxv
  integer mul(0:maxq-1,0:maxq-1)
  integer, parameter :: nonzer = 1
  integer p
  integer px(-1:maxdeg)
  integer q
  integer r
  integer sub(0:maxq-1,0:maxq-1)
  integer term
  integer v(0:maxv)

  common /field/ p, q, add, mul, sub

  save /field/

  e = px(deg)
!
!  The polynomial H is PX**(J-1), which is the value of B on arrival.
!
!  In section 3.3, the values of Hi are defined with a minus sign:
!  don't forget this if you use them later!
!
  do i = -1, b(deg)
    h(i) = b(i)
  end do

  bigm = h(deg)
!
!  Now multiply B by PX so B becomes PX**J.
!
!  In section 2.3, the values of Bi are defined with a minus sign:
!  don't forget this if you use them later!
!
  call plymul ( px, b, b )
  m = b(deg)
!
!  We don't use J explicitly anywhere, but here it is just in case.
!
  j = m / e
!
!  Now choose a value of Kj as defined in section 3.3.
!  We must have 0 <= Kj < E*J = M.
!  The limit condition on Kj does not seem very relevant
!  in this program.
!
  kj = bigm
!
!  Now choose values of V in accordance with the conditions in
!  section 3.3
!
  v(0:kj-1) = 0
  v(kj) = 1

  if ( kj < bigm ) then

    term = sub ( 0, h(kj) )

    do r = kj + 1, bigm - 1

      v(r) = arbit
!
!  Check the condition of section 3.3,
!  remembering that the H's have the opposite sign.
!
      term = sub ( term, mul ( h(r), v(r) ) )

    end do
!
!  Now V(BIGM) is anything but TERM.
!
    v(bigm) = add ( nonzer, term )
    v(bigm+1:m-1) = arbit

  else

    v(kj+1:m-1) = arbit

  end if
!
!  Calculate the remaining V's using the recursion of section 2.3,
!  remembering that the B's have the opposite sign.
!
  do r = 0, maxv - m
    term = 0
    do i = 0, m - 1
      term = sub ( term, mul (b(i), v(r+i)) )
    end do
    v(r+m) = term
  end do

  return
end
function exactf ( n, dimen )

!*****************************************************************************80
!
!! EXACTF returns the exact integral of a scalar test integrand.
!
!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederrieter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Implementation and Tests of Low-Discrepancy Sequences,
!    ACM Transactions on Modeling and Computer Simulation,
!    Volume 2, Number 3, pages 195-213, 1992.
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, pages 494-495, 1994.
!
!  Parameters:
!
!    Input, integer N, the index of the test integral.
!
!    Input, integer DIMEN, the spatial dimension.
!
!    Output, real EXACTF, the exact value of the integral.
!
  implicit none

  integer dimen
  real ( kind = 8 ) exactf
  integer i
  integer n
  real ( kind = 8 ) x
!
!  Integral over [0,+1]**DIMEN Product ( 1 <= i <= DIMEN ) 
!  abs ( 4 * x(i) - 2 ) dx
! 
  if ( n == 1 ) then

    exactf = 1.0D+00
!
!  Integral over [0,+1]**DIMEN Product ( 1 <= i <= DIMEN  )
!  i * cos ( i * x(i) ) dx
!
  else if ( n == 2 ) then

    exactf = 1.0D+00
    do i = 1, dimen
      exactf = exactf * sin ( real ( i, kind = 8 ) )
    end do
!
!  Integral over [0,+1]**DIMEN Product ( 1 <= i <= DIMEN )
!  T(N(i))( 2 * x(i) - 1 ) dx
!
!  where T(N)(X) is the Chebyshev polynomial of degree N,
!  and N(i) = mod ( i, 4 ) + 1.
!
  else if ( n == 3 ) then

    exactf = 0.0D+00
!
!  Integral over [0,+1]**DIMEN Sum ( 1 <= i <= DIMEN )
!  (-1)**i * Product ( 1 <= j <= i ) x(j) dx
!
  else if ( n == 4 ) then

    x = 1.0D+00 / real ( 2**dimen, kind = 8 )
    if ( mod ( dimen, 2 ) == 0 ) then
      exactf = ( x - 1.0D+00 ) / 3.0D+00
    else
      exactf = ( x + 1.0D+00 ) / 3.0D+00
    end if

  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
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
!    Output, integer IUNIT.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

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
subroutine golo ( quasi )

!*****************************************************************************80
!
!! GOLO generates a new quasi-random vector on each call.
!
!  Discussion:
!
!    Before the first call to this routine, a call must be made
!    to subroutine INLO to carry out some initializations.
!
!    Polynomials stored as arrays have the coefficient of degree n 
!    in POLY(N), and the degree of the polynomial in POLY(-1).  
!    The parameter DEG is just to remind us of this last fact.  
!    A polynomial which is identically 0 is given degree -1.
!
!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederrieter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, pages 494-495, 1994.
!
!    Harald Niederreiter,
!    Low-discrepancy and low-dispersion sequences,
!    Journal of Number Theory,
!    Volume 30, 1988, pages 51-70.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) QUASI(*), the next vector in the sequence.
!
!  Local Parameters:
!
!    Local, integer MAXDEG, the highest degree of polynomial to be handled. 
! 
!    Local, integer MAXDIM, the maximum dimension that will be used.
!
!    Local, integer MAXFIG, the maximum number of base-Q digits we can handle.
!
!    Local, integer MAXINT, the largest fixed point integer we can represent.
!
!    Local, integer MAXQ, the order of the largest field to be handled.
!
!    Local, integer NBITS, the number of bits in a fixed-point integer, not
!    counting the sign.
!
!  Global parameters in /COMM/:
!
!    Global, integer C(MAXDIM,MAXFIG,0:MAXFIG-1), the values of 
!    Niederreiter's C(I,J,R).
!
!    Global, integer COUNT(0:MAXFIG-1), the index of the current item 
!    in the sequence, expressed as an array of base-Q digits.  COUNT(R)
!    is the same as Niederreiter's AR(N) (page 54) except that N is implicit.
!
!    Global, integer D(MAXDIM,MAXFIG).
!
!    Global, integer NEXTQ(MAXDIM), the numerators of the next item in 
!    the series.  These are like Niederreiter's XI(N) (page 54) except that
!    N is implicit, and the NEXTQ are integers.  To obtain the values of 
!    XI(N), multiply by RECIP.
!
!    Global, integer QPOW(MAXFIG), to speed things up a bit. 
!    QPOW(I) = Q ** (NFIGS-I).
!
!    Global, integer DIMEN, the dimension of the sequence to be generated.
!
!    Global, integer NFIGS, the number of base Q digits we are using.
!
!    Global, real ( kind = 8 ) RECIP = 1.0 / Q**NFIGS.
!
  implicit none

  integer, parameter :: deg = -1
  integer, parameter :: maxdeg = 50
  integer, parameter :: maxdim = 20
  integer, parameter :: maxfig = 20
  integer, parameter :: maxq = 50
  integer, parameter :: nbits = 31

  integer add(0:maxq-1,0:maxq-1)
  integer c(maxdim,maxfig,0:maxfig-1)
  integer count(0:maxfig-1)
  integer d(maxdim, maxfig)
  integer diff
  integer dimen
  integer i
  integer j
  integer mul(0:maxq-1,0:maxq-1)
  integer nextq(maxdim)
  integer nfigs
  integer nq
  integer oldcnt
  integer p
  integer q
  integer qpow(maxfig)
  real ( kind = 8 ) quasi(*)
  integer r
  real ( kind = 8 ) recip
  integer sub(0:maxq-1,0:maxq-1)

  common /comm/ c, count, d, nextq, qpow, dimen, nfigs, recip
  common /field/ p, q, add, mul, sub

  save /comm/
  save /field/
!
!  Multiply the numerators in NEXTQ by RECIP to get the next
!  quasi-random vector.
!
  quasi(1:dimen) = real ( nextq(1:dimen), kind = 8 ) * recip
!
!  Update COUNT, treated as a base-Q integer.  Instead of
!  recalculating the values of D from scratch, we update
!  them for each digit of COUNT which changes.  In terms of
!  Niederreiter page 54, NEXTQ(I) corresponds to XI(N), with
!  N being implicit, and D(I,J) corresponds to XI(N,J), again
!  with N implicit.  Finally COUNT(R) corresponds to AR(N).
!
  r = 0

  do

    if ( nfigs <= r ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GOLO - Fatal error!'
      write ( *, '(a)' ) '  Too many calls!'
      stop
    end if

    oldcnt = count(r)

    if ( count(r) < q - 1 ) then
      count(r) = count(r) + 1
    else
      count(r) = 0
    end if

    diff = sub ( count(r), oldcnt )
!
!  Digit R has just changed.  DIFF says how much it changed
!  by.  We use this to update the values of array D.
!
    do i = 1, dimen
      do j = 1, nfigs
        d(i,j) = add ( d(i,j), mul ( c(i,j,r), diff ) )
      end do
    end do
!
!  If COUNT(R) is now zero, we must propagate the carry.
!
    if ( count(r) /= 0 ) then
      exit
    end if

    r = r + 1

  end do
!
!  Now use the updated values of D to calculate NEXTQ.
!  Array QPOW helps to speed things up a little:
!  QPOW(J) is Q ** (NFIGS-J).
!
  do i = 1, dimen
    nq = 0
    do j = 1, nfigs
      nq = nq + d(i,j) * qpow(j)
    end do
    nextq(i) = nq
  end do

  return
end
function i4_characteristic ( q )

!*****************************************************************************80
!
!! I4_CHARACTERISTIC gives the characteristic for an integer.
!
!  Discussion:
!
!    For any positive integer Q, the characteristic is:
!
!    Q, if Q is a prime;
!    P, if Q = P**N for some prime P and some integer N;
!    0, otherwise, that is, if Q is negative, 0, 1, or the product
!       of more than one distinct prime.
!
!    A faster code would only consider prime values of I,
!    but that entails storing a table of primes and limiting the
!    size of Q.  Simplicity and flexibility for now!
!
!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederrieter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738:
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, pages 494-495, 1994.
!
!  Parameters:
!
!    Input, integer Q, the value to be tested.
!
!    Output, integer I4_CHARACTERISTIC, the characteristic of Q.
!
  implicit none

  integer i
  integer i4_characteristic
  integer i_max
  integer q
  integer q_copy

  if ( q <= 1 ) then
    i4_characteristic = 0
    return
  end if
!
!  If Q is not prime, then there is at least one prime factor
!  of Q no greater than SQRT(Q)+1.
!
  i_max = int ( sqrt ( real ( q, kind = 8 ) ) ) + 1
  q_copy = q

  do i = 2, i_max

    if ( mod ( q_copy, i ) == 0 ) then

      do while ( mod ( q_copy, i ) == 0 )
        q_copy = q_copy / i
      end do

      if ( q_copy == 1 ) then
        i4_characteristic = i
      else
        i4_characteristic = 0
      end if

      return

    end if

  end do
!
!  If no factor was found, then Q is prime.
!
  i4_characteristic = q

  return
end
subroutine inlo ( dim, base, skip )

!*****************************************************************************80
!
!! INLO calculates the values of C(I,J,R).
!
!  Discussion:
!
!    This subroutine calculates the values of Niederreiter's
!    C(I,J,R) and performs other initialization necessary
!    before calling GOLO.
!
!    Polynomials stored as arrays have the coefficient of degree n 
!    in POLY(N), and the degree of the polynomial in POLY(-1).  
!    The parameter DEG is just to remind us of this last fact.  
!    A polynomial which is identically 0 is given degree -1.
!
!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederrieter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, pages 494-495, 1994.
!
!    Harald Niederreiter,
!    Low-discrepancy and low-dispersion sequences,
!    Journal of Number Theory,
!    Volume 30, 1988, pages 51-70.
!
!  Parameters:
!
!    Input, integer DIM, the dimension of the sequence to be generated.
!    The value of DIM is copied into DIMEN in the common block.
!
!    Input, integer BASE, the prime or prime-power base to be used.
!
!    Input, integer SKIP, the number of values to throw away at the 
!    beginning of the sequence.
!
!  Local Parameters:
!
!    Local, integer MAXQ, the order of the largest field to be handled.
!
!    Local, integer MAXDEG, the highest degree of polynomial
!    to be handled. 
!
!    Local, integer MAXDIM, the maximum dimension that will be used.
!
!    Local, integer MAXFIG, the maximum number of base-Q digits we can handle.
!
!    Local, integer MAXINT, the largest fixed point integer we can represent.
!
!    Local, integer NBITS, the number of bits in a fixed-point integer, not
!    counting the sign.
!
!  Global parameters in /COMM/:
!
!    Global, integer C(MAXDIM,MAXFIG,0:MAXFIG-1), the values of 
!    Niederreiter's C(I,J,R).
!
!    Global, integer COUNT(0:MAXFIG-1), the index of the current item 
!    in the sequence, expressed as an array of base-Q digits.  COUNT(R)
!    is the same as Niederreiter's AR(N) (page 54) except that N is implicit.
!
!    Global, integer D(MAXDIM,MAXFIG).
!
!    Global, integer NEXTQ(MAXDIM), the numerators of the next item in 
!    the series.  These are like Niederreiter's XI(N) (page 54) except that
!    N is implicit, and the NEXTQ are integers.  To obtain the values of 
!    XI(N), multiply by RECIP.
!
!    Global, integer QPOW(MAXFIG), to speed things up a bit. 
!    QPOW(I) = Q ** (NFIGS-I).
!
!    Global, integer DIMEN, the dimension of the sequence to be generated.
!
!    Global, integer NFIGS, the number of base Q digits we are using.
!
!    Global, real ( kind = 8 ) RECIP = 1.0 / Q**NFIGS.
!
  implicit none

  integer, parameter :: deg = -1
  integer, parameter :: maxdeg = 50
  integer, parameter :: maxdim = 20
  integer, parameter :: maxfig = 20
  integer, parameter :: maxq = 50
  integer, parameter :: nbits = 31

  integer add(0:maxq-1,0:maxq-1)
  integer base
  integer c(maxdim, maxfig,0:maxfig-1)
  integer count(0:maxfig-1)
  integer d(maxdim, maxfig)
  integer dim
  integer dimen
  integer i
  integer i4_characteristic
  integer j
  integer mul(0:maxq-1,0:maxq-1)
  integer nextq(maxdim)
  integer nfigs
  integer nq
  integer p
  integer q
  integer qpow(maxfig)
  integer r
  real ( kind = 8 ) recip
  integer skip
  integer sub(0:maxq-1,0:maxq-1)
  real ( kind = 8 ) temp

  common /comm/ c, count, d, nextq, qpow, dimen, nfigs, recip
  common /field/ p, q, add, mul, sub

  save /comm/
  save /field/

  dimen = dim

  if ( dimen <= 0 .or. maxdim < dimen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INLO - Fatal error!'
    write ( *, '(a)' ) '  Bad spatial dimension.'
    stop
  end if

  if ( i4_characteristic ( base ) == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INLO - Fatal error!'
    write ( *, '(a)' ) '  Base not prime power or out of range.'
    stop
  end if

  call setfld ( base )
!
!  Calculate how many figures to use in base Q = BASE
!
  temp = log ( 2.0D+00**nbits - 1.0D+00 ) / log ( real ( q, kind = 8 ) )

  nfigs = min ( maxfig, int ( temp ) )
!
!  Calculate the C array.
!
  call calcc
!
!  Set RECIP.
!
  recip = 1.0D+00 / real ( q**nfigs, kind = 8 )
!
!  Set QPOW(I) = Q ** (NFIGS-I).
!
  qpow(nfigs) = 1
  do i = nfigs-1, 1, -1
    qpow(i) = q * qpow(i+1)
  end do
!
!  Initialize COUNT.
!
  i = skip

  do r = 0, nfigs-1
    count(r) = mod ( i, q )
    i = i / q
  end do

  if ( i /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INLO - Fatal error!'
    write ( *, '(a)' ) '  SKIP is too long'
    stop
  end if
!
!  Initialize D.
!
  d(1:dimen,1:nfigs) = 0

  do r = 0, nfigs-1
    if ( count(r) /= 0 ) then
      do i = 1, dimen
        do j = 1, nfigs
          d(i,j) = add ( d(i,j), mul ( c(i,j,r), count(r) ) )
        end do
      end do
    end if
  end do
!
!  Initialize NEXTQ.
!
  do i = 1, dimen
    nq = 0
    do j = 1, nfigs
      nq = nq + d(i,j) * qpow(j)
    end do
    nextq(i) = nq
  end do

  return
end
subroutine plymul ( pa, pb, pc )

!*****************************************************************************80
!
!! PLYMUL computes the product of two polynomials in the field Q.
!
!  Discussion:
!
!    The polynomials have coefficients in the field of order Q.
!
!    Polynomials stored as arrays have the coefficient of degree N
!    in POLY(N), and the degree of the polynomial in POLY(-1).  
!
!    The parameter DEG is just to remind us of this last fact.  
!    A polynomial which is identically 0 is given degree -1.
!
!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederrieter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, pages 494-495, 1994.
!
!  Parameters:
!
!    Input, integer PA(-1:MAXDEG), PB(-1:MAXDEG), two polynomials
!    to be multiplied.
!
!    Output, integer PC(-1:MAXDEG), the product of PA and PB, 
!    in the field Q.
!
!  Local Parameters:
!
!    Local, integer MAXDEG, the highest degree of polynomial
!    to be handled. 
!
!    Local, integer MAXQ, the order of the largest field to be handled.
!
  implicit none

  integer, parameter :: deg = -1
  integer, parameter :: maxdeg = 50
  integer, parameter :: maxq = 50

  integer add(0:maxq-1,0:maxq-1)
  integer dega
  integer degb
  integer degc
  integer i
  integer j
  integer mul(0:maxq-1,0:maxq-1)
  integer p
  integer pa(-1:maxdeg)
  integer pb(-1:maxdeg)
  integer pc(-1:maxdeg)
  integer pt(-1:maxdeg)
  integer q
  integer sub(0:maxq-1,0:maxq-1)
  integer term

  common /field/ p, q, add, mul, sub

  save /field/ 

  dega = pa(deg)
  degb = pb(deg)

  if ( dega == -1 .or. degb == -1 ) then
    degc = -1
  else
    degc = dega + degb
  end if

  if ( maxdeg < degc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLYMUL - Fatal error!'
    write ( *, '(a)' ) '  The degree of the product exceeds MAXDEG.'
    stop
  end if

  do i = 0, degc
    term = 0
    do j = max ( 0, i - dega ), min ( degb, i )
      term = add ( term, mul ( pa(i-j), pb(j) ) )
    end do
    pt(i) = term
  end do

  pc(deg) = degc
  do i = 0, degc
    pc(i) = pt(i)
  end do

  do i = degc + 1, maxdeg
    pc(i) = 0
  end do

  return
end
subroutine setfld ( qin )

!*****************************************************************************80
!
!! SETFLD sets up arithmetic tables for the field of order QIN.
!
!  Discussion:
!
!    If necessary, it reads precalculated tables from the file
!    "gfarit.txt".  These precalculated tables
!    are supposed to have been put there by GFARIT.
!
!    Polynomials stored as arrays have the coefficient of degree n 
!    in POLY(N), and the degree of the polynomial in POLY(-1).  
!    The parameter DEG is just to remind us of this last fact.  
!    A polynomial which is identically 0 is given degree -1.
!
!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederrieter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, pages 494-495, 1994.
!
!  Parameters:
!
!    Input, integer QIN, the order of the field.
!
!  Local Parameters:
!
!    Local, integer MAXDEG, the highest degree of polynomial to be handled.  
!
!    Local, integer MAXQ, the order of the largest field to be handled.
!
  implicit none

  integer, parameter :: deg = -1
  integer, parameter :: maxdeg = 50
  integer, parameter :: maxq = 50

  integer add(0:maxq-1,0:maxq-1)
  character ( len = 80 ) :: file_in_name = 'gfarit.txt'
  integer file_in_unit
  integer i
  integer i4_characteristic
  integer ios
  integer j
  integer mul(0:maxq-1,0:maxq-1)
  integer n
  integer p
  integer q
  integer qin
  integer sub(0:maxq-1,0:maxq-1)

  common /field/ p, q, add, mul, sub

  save /field/

  if ( qin <= 1 .or. maxq < qin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETFLD - Fatal error!'
    write ( *, '(a)' ) '  Bad value of Q!'
    stop
  end if

  q = qin
  p = i4_characteristic ( q )

  if ( p == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETFLD - Fatal error!'
    write ( *, '(a,i8)' ) '  There is no field of order Q = ', q
    stop
  end if
!
!  Set up to handle a field of prime order: calculate ADD and MUL.
!
  if ( p == q ) then
    do i = 0, q-1
      do j = 0, q-1
        add(i,j) = mod ( i + j, p)
        mul(i,j) = mod ( i * j, p)
      end do
    end do
!
!  Set up to handle a field of prime-power order: tables for
!  ADD and MUL are in the file 'gfarit.txt'.
!
  else

    call get_unit ( file_in_unit )

    open ( unit = file_in_unit, file = file_in_name, status = 'old' )

    do

      read ( file_in_unit, '(20i3)', iostat = ios ) n

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETFLD - Fatal error!'
        write ( *, '(a,i8)' ) '  Could not find tables for Q =', q
        stop
      end if

      do i = 0, n-1
        read ( file_in_unit, '(20i3)' ) add(i,0:n-1)
      end do

      do i = 0, n-1
        read ( file_in_unit, '(20i3)' ) mul(i,0:n-1)
      end do

      if ( n == q ) then
        exit
      end if

    end do

    close ( unit = file_in_unit )

  end if
!
!  Use the addition table to set the subtraction table.
!
  do i = 0, q-1
    do j = 0, q-1
      sub(add(i,j), i) = j
    end do
  end do

  return
end
function testf ( n, dimen, quasi )

!*****************************************************************************80
!
!! TESTF evaluates a scalar test integrand of a multidimensional argument.
!!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederrieter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Implementation and Tests of Low-Discrepancy Sequences,
!    ACM Transactions on Modeling and Computer Simulation,
!    Volume 2, Number 3, pages 195-213, 1992.
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, pages 494-495, 1994.
!
!  Parameters:
!
!    Input, integer N, the index of the test function.
!
!    Input, integer DIMEN, the spatial dimension.
!
!    Input, real ( kind = 8 ) QUASI(DIMEN), the argument.
!
!    Output, real ( kind = 8 ) TESTF, the integrand at QUASI(*).
!
  implicit none

  integer dimen

  integer i
  integer ii
  integer n
  real ( kind = 8 ) quasi(dimen)
  real ( kind = 8 ) testf
  real ( kind = 8 ) x
!
!  Integral over [0,+1]**DIMEN Product ( 1 <= i <= DIMEN ) 
!  abs ( 4 * x(i) - 2 ) dx
! 
  if ( n == 1 ) then

    testf = product ( abs ( 4.0D+00 * quasi(1:dimen) - 2.0D+00 ) )
!
!  Integral over [0,+1]**DIMEN Product ( 1 <= i <= DIMEN  )
!  i * cos ( i * x(i) ) dx
!
  else if ( n == 2 ) then

    testf = 1.0D+00
    do i = 1, dimen
      testf = testf * real ( i, kind = 8 ) &
        * cos ( real ( i, kind = 8 ) * quasi(i) )
    end do
!
!  Integral over [0,+1]**DIMEN Product ( 1 <= i <= DIMEN )
!  T(N(i))( 2 * x(i) - 1 ) dx
!
!  where T(N)(X) is the Chebyshev polynomial of degree N,
!  and N(i) = mod ( i, 4 ) + 1.
!
  else if ( n == 3 ) then

    testf = 1.0D+00

    do i = 1, dimen

      x = 2.0D+00 * quasi(i) - 1.0D+00

      ii = mod ( i, 4 )

      if ( ii == 1 ) then
        testf = testf * x
      else if ( ii == 2 ) then
        testf = testf * ( 2.0D+00 * x**2 - 1.0D+00 )
      else if ( ii == 3 ) then
        testf = testf * ( 4.0D+00 * x**2 - 3.0D+00 ) * x
      else if ( ii == 0 ) then
        testf = testf * ( 8.0D+00 * x**4 - 8.0D+00 * x**2 + 1.0D+00 )
      end if

    end do
!
!  Integral over [0,+1]**DIMEN Sum ( 1 <= i <= DIMEN )
!  (-1)**i * Product ( 1 <= j <= i ) x(j) dx
!
  else if ( n == 4 ) then

    testf = 0.0D+00
    x = 1.0D+00
    do i = 1, dimen
      x = - x * quasi(i)
      testf = testf + x
    end do

  else

    testf = 0.0D+00

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
