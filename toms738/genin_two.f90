program main

!*****************************************************************************80
!
!! MAIN is the main program for GENIN_TWO.
!
!  Discussion:
!
!    GENIN_TWO is the driver for the base-2 program.
!
!    This program tests the accuracy of numerical integration
!    using the low-discrepancy binary sequences of
!    H. Niederreiter (1988) as implemented in INLO2, GOLO2, and
!    related programs.  Various possible test integrals are
!    provided by the function TESTF.  GENIN2 generates only
!    sequences with base 2.
!
!    These program assumes that your computer's word length
!    is 31 bits, excluding sign.  If this is not the case,
!    modify the parameter NBITS throughout accordingly.
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
!    Local, integer MAXDIM, the maximum dimension that will be used.  
!
!    Local, integer POWER, is used in a possible warm-up formula.
!
  implicit none

  integer, parameter :: maxdim = 20
  integer, parameter :: power = 12

  real correct
  real error
  real estimate
  integer dimen
  real exactf
  character ( len = 80 ) :: file_out_name = 'genin_two_prb_output.txt'
  integer file_out_unit
  integer i
  integer num
  integer seqlen
  integer skip
  integer step
  integer pbase
  real quasi(maxdim)
  real testf
  double precision total

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GENIN_TWO'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Generate a Niederreiter sequence for a given'
  write ( *, '(a)' ) '  spatial dimension.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Special version of GENIN for the case BASE = 2.'

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace' )

  do

    do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Choose a test integral (1 to 4) or 0 to quit:'

      read ( *, * ) num

      if ( num <= 4 ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GENIN_TWO - Error!'
      write ( *, '(a)' ) '  There is no such test integral.'

    end do

    if ( num <= 0 ) then
      exit
    end if
!
!  Each test integral is parameterized by its dimension.
!
    do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Enter the spatial dimension'
      write ( *, '(a,i6)' ) '  between 1 and ', maxdim

      read ( *, * ) dimen

      if ( 0 < dimen .and. dimen <= maxdim ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GENIN_TWO - Error!'
      write ( *, '(a,i6)' ) '  The dimension may not exceed ', maxdim

    end do
!
!  The sequence length is the number of quasi-random points used to 
!  estimate the integral, excluding warm-up.
!
!  Some users may wish to rewrite the driver to test a [heuristic] 
!  "convergence" criterion, stopping the generation of points
!  when it is passed or when a specified number of points have been 
!  generated, whichever occurs first.
!
    do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Choose the sequence length:'
!
!  Except when comparing results with other bases, we suggest taking 
!  SEQLEN to be a power of 2.
!
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  If you do not have a preference, we'
      write ( *, '(a)' ) '  suggest using a large power of two, such as:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  2**10 = ', 2**10
      write ( *, '(a,i12)' ) '  2**15 = ', 2**15
      write ( *, '(a,i12)' ) '  2**20 = ', 2**20
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Enter the sequence length:'

      read ( *, * ) seqlen

      if ( 0 < seqlen ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GENIN_TWO - Error!'
      write ( *, '(a)' ) '  The sequence length must be strictly positive.'

    end do

    do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Choose the number of values to skip:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  There is reason to believe that BASE**E, '
      write ( *, '(a)' ) '  where E is defined for example in '
      write ( *, '(a)' ) '  Bratley, Fox, and Niederreiter [1991], '
      write ( *, '(a)' ) '  would be a good choice.  Our formula has '
      write ( *, '(a)' ) '  has the form:'
      write ( *, '(a)' ) '    SKIP = 2 ** POWER,'
      write ( *, '(a)' ) '  where POWER is chosen so that SKIP comes nowhere '
      write ( *, '(a)' ) '  near the largest possible machine-representable'
      write ( *, '(a)' ) '  integer.  It does not hurt to choose '
      write ( *, '(a)' ) '  POWER larger than E, because skipping is '
      write ( *, '(a)' ) '  done implicitly in O(1) time. '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The maximum value of E for a fixed dimension '
      write ( *, '(a)' ) '  S grows like log S.  We allow some "fat" for '
      write ( *, '(a)' ) '  the implicit constant. '
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Numerically, 2**POWER = ', 2 ** power
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Enter SKIP (not necessarily the value above)'

      read ( *, * ) skip

      if ( 0 <= skip ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GENIN_TWO - Error!'
      write ( *, '(a)' ) '  The number must be nonnegative.'

    end do
!
!  Calculate the values of C(I,J,R).
!
    call inlo2 ( dimen, skip )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENIN_TWO:'
    write ( *, '(a)' ) '  Initialization completed.'
!
!  Write titles and the exact value of the integral
!
    correct = exactf ( num, dimen )

    write ( file_out_unit, '(a)' ) ' '
    write ( file_out_unit, '(a,i2)' ) '  Test integral ', num
    write ( file_out_unit, '(a)' ) ' '
    write ( file_out_unit, '(a,i6)' ) '  Spatial dimension = ', dimen
    write ( file_out_unit, '(a,i6)' ) '  Base =              ', 2
    write ( file_out_unit, '(a,i6)' ) '  Sequence length =   ', seqlen
    write ( file_out_unit, '(a,i6)' ) '  Skip size =         ', skip
    write ( file_out_unit, '(a)' ) ' '
    write ( file_out_unit, '(a,g16.7)' ) '  Correct value is ', correct
    write ( file_out_unit, '(a)' ) ' '
    write ( file_out_unit, '(a)' ) &
      '      Iteration     Estimate of integral      Error'
    write ( file_out_unit, '(a)' ) ' '
!
!  Now estimate the integral
!
    total = 0.0

    pbase = 2**6

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Odd-looking iteration numbers are powers of 2.'

    step = 500

    do i = 1, seqlen

      call golo2 ( quasi )

      total = total + testf ( num, dimen, quasi )

      if ( mod ( i, step ) == 0 ) then
        estimate = total / real ( i )
        error = abs ( estimate - correct )
        write ( file_out_unit, '(i12,3g24.7)' ) i, estimate, error
      end if

      if ( mod ( i, pbase ) == 0 ) then
        estimate = total / real ( i )
        error = abs ( estimate - correct )
        write ( file_out_unit, '(i12,3g24.7)' ) i, estimate, error
        pbase = 2 * pbase
      end if
!
!  There is reason to believe that the subsequence
!  of estimates along powers of the base [here 2]
!  converges faster than the original sequence or
!  the subsequence corresponding to STEP.
!
      if ( i == 5000 ) then
        step = 1000
      end if

      if ( i == 10000 ) then
        step = 5000
      end if

    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENIN_TWO:'
    write ( *, '(a)' ) '  Iteration completed.'

  end do

  close ( unit = file_out_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GENIN_TWO:'
  write ( *, '(a)' ) '  Output file name is ' // trim ( file_out_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GENIN_TWO:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine calcc2

!*****************************************************************************80
!
!! CALCC2 computes values of the constants C(I,J,R).
!
!  Discussion:
!
!    This program calculates the values of the constants C(I,J,R).
!    As far as possible, we use Niederreiter's notation.
!    For each value of I, we first calculate all the corresponding
!    values of C.  These are held in the array CI.  All these
!    values are either 0 or 1.  Next we pack the values into the
!    array CJ, in such a way that CJ(I,R) holds the values of C
!    for the indicated values of I and R and for every value of
!    J from 1 to NBITS.  The most significant bit of CJ(I,R)
!    (not counting the sign bit) is C(I,1,R) and the least
!    significant bit is C(I,NBITS,R).
!
!  Reference:
!
!    Rudolf Lidl and Harald Niederreiter, 
!    Finite Fields,
!    Cambridge University Press, 1984, page 553.
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
!    Local, integer MAXE, we need MAXDIM irreducible polynomials over Z2.
!    MAXE is the highest degree among these.
!
!    Local, integer MAXQ, the order of the largest field to be handled.
!
!    Local, integer MAXV, the maximum possible index used in V.
!
!    Local, integer NBITS, the number of bits (not counting the sign) in a
!    fixed-point integer.
!
!  Global Parameters:
!
!    In common block /COMM2/:
!
!    Global, integer CJ(MAXDIM,0:NBITS-1), the packed values of 
!    Niederreiter's C(I,J,R)
!
!    Global, integer DIMEN, the dimension of the sequence to be generated
!
!    Global, integer COUNT, the index of the current item in the sequence,
!    expressed as an array of bits.  COUNT(R) is the same as Niederreiter's
!    AR(N) (page 54) except that N is implicit.
!
!    Global, integer NEXTQ(MAXDIM), the numerators of the next item 
!    in the series.  These are like Niederreiter's XI(N) (page 54) except 
!    that N is implicit, and the NEXTQ are integers.  To obtain
!    the values of XI(N), multiply by RECIP.
!
  implicit none

  integer, parameter :: deg = -1
  integer, parameter :: maxdeg = 50
  integer, parameter :: maxdim = 20
  integer, parameter :: maxe = 5
  integer, parameter :: maxq = 50
  integer, parameter :: nbits = 31

  integer, parameter :: maxv = nbits + maxe

  integer add(0:maxq-1,0:maxq-1)
  integer b(-1:maxdeg)
  integer ci(nbits,0:nbits-1)
  integer cj(maxdim,0:nbits-1)
  integer count
  integer dimen
  integer e
  integer i
  integer, save, dimension ( maxdim, -1:maxe ) :: irred
  integer j
  integer mul(0:maxq-1,0:maxq-1)
  integer nextq(maxdim)
  integer p
  integer px(-1:maxdeg)
  integer q
  integer r
  integer sub(0:maxq-1,0:maxq-1)
  integer term
  integer u
  integer v(0:maxv)

  common /comm2/ cj, dimen, count, nextq
  common /field/ p, q, add, mul, sub

  save /comm2/
  save /field/
!
!  Here we supply the coefficients and the
!  degrees of the first 12 irreducible polynomials over Z2.
!
!  They are taken from Lidl and Niederreiter.
!
!  The order of these polynomials is the same as the order in
!  file 'gfplys.dat' used by the general program.
!
!  In this block PX(I, -1) is the degree of the Ith polynomial,
!  and PX(I, N) is the coefficient of x**n in the Ith polynomial.
!
  irred( 1,-1:1) = (/ 1,0,1 /)
  irred( 2,-1:1) = (/ 1,1,1 /)
  irred( 3,-1:2) = (/ 2,1,1,1 /)
  irred( 4,-1:3) = (/ 3,1,1,0,1 /)
  irred( 5,-1:3) = (/ 3,1,0,1,1 /)
  irred( 6,-1:4) = (/ 4,1,1,0,0,1 /)
  irred( 7,-1:4) = (/ 4,1,0,0,1,1 /)
  irred( 8,-1:4) = (/ 4,1,1,1,1,1 /)
  irred( 9,-1:5) = (/ 5,1,0,1,0,0,1 /)
  irred(10,-1:5) = (/ 5,1,0,0,1,0,1 /)
  irred(11,-1:5) = (/ 5,1,1,1,1,0,1 /)
  irred(12,-1:5) = (/ 5,1,1,1,0,1,1 /)
!
!  Prepare to work in Z2.
!
  call setfld ( 2 )

  do i = 1, dimen
!
!  For each dimension, we need to calculate powers of an
!  appropriate irreducible polynomial:  see Niederreiter
!  page 65, just below equation (19).
!
!  Copy the appropriate irreducible polynomial into PX,
!  and its degree into E.  Set polynomial B = PX ** 0 = 1.
!  M is the degree of B.  Subsequently B will hold higher
!  powers of PX.
!
    e = irred(i,deg)

    do j = -1, e
      px(j) = irred(i,j)
    end do

    b(deg) = 0
    b(0) = 1
!
!  Niederreiter (page 56, after equation (7), defines two
!  variables Q and U.  We do not need Q explicitly, but we do need U.
!
    u = 0

    do j = 1, nbits
!
!  If U = 0, we need to set B to the next power of PX
!  and recalculate V.  This is done by subroutine CALCV.
!
      if ( u == 0 ) then
        call calcv ( px, b, v, maxv )
      end if
!
!  Now C is obtained from V.  Niederreiter
!  obtains A from V (page 65, near the bottom), and then gets
!  C from A (page 56, equation (7)).  However this can be done
!  in one step.  Here CI(J,R) corresponds to Niederreiter's C(I,J,R).
!
      do r = 0, nbits-1
        ci(j,r) = v(r+u)
      end do
!
!  Increment U.  
!
!  If U = E, then U = 0 and in Niederreiter's
!  paper Q = Q + 1.  Here, however, Q is not used explicitly.
!
      u = u + 1
      if ( u == e ) then
        u = 0
      end if

    end do
!
!  The array CI now holds the values of C(I,J,R) for this value
!  of I.  We pack them into array CJ so that CJ(I,R) holds all
!  the values of C(I,J,R) for J from 1 to NBITS.
!
    do r = 0, nbits-1
      term = 0
      do j = 1, nbits
        term = 2 * term + ci(j,r)
      end do
      cj(i,r) = term
    end do

  end do

  return
end
subroutine calcv ( px, b, v, maxv )

!*****************************************************************************80
!
!! CALCV calculates the value of the constants V(J,R).
!
!  Discussion:
!
!    This program calculates the values of the constants V(J,R) as
!    described in BFN section 3.3.  It is called from either CALCC or
!    CALCC2.  The values transmitted through common /FIELD/ determine
!    which field we are working in.
!
!    Polynomials stored as arrays have the coefficient of degree n 
!    in POLY(N), and the degree of the polynomial in POLY(-1).  The 
!    parameter DEG is just to remind us of this last fact.  A 
!    polynomial which is identically 0 is given degree -1.
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
!    for the dimension currently being considered.  The degree of PX 
!    will be called E.
!
!    Input/output, integer B(-1:MAXDEG).  On input, B is the polynomial 
!    defined in section 2.3 of BFN.  The degree of B implicitly defines 
!    the parameter J of section 3.3, by degree(B) = E*(J-1).  On output,
!    B has been multiplied by PX, so its degree is now E * J.
!
!    Output, integer V(0:MAXV), the computed V array.
!
!    Input, integer MAXV gives the dimension of the array V.
!
!  Local Parameters:
!
!    Local, integer ARBIT, indicates where the user can place
!    an arbitrary element of the field of order Q.  This means 
!    0 <= ARBIT < Q.  
!
!    Local, integer BIGM, is the M used in section 3.3.
!    It differs from the [little] m used in section 2.3,
!    denoted here by M.
!
!    Local, integer MAXDEG, the highest degree of polynomial
!    to be handled. 
!
!    Local, integer MAXQ, the order of the largest field to be handled.
!
!    Local, integer NONZER, shows where the must put an arbitrary 
!    non-zero element of the field.  For the code, this means 
!    0 < NONZER < Q.
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
!  Multiply B by PX so B becomes PX**J.
!  In section 2.3, the values of Bi are defined with a minus sign:
!  don't forget this if you use them later!
!
  call plymul (px, b, b)
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
!  Choose values of V in accordance with the conditions in section 3.3
!
  do r = 0, kj-1
    v(r) = 0
  end do
  v(kj) = 1

  if ( kj < bigm ) then

    term = sub ( 0, h(kj) )

    do r = kj+1, bigm-1

      v(r) = arbit
!
!  Check the condition of section 3.3,
!  remembering that the H's have the opposite sign.
!
      term = sub ( term, mul ( h(r), v(r) ) )

    end do
!
!  Now V(BIGM) is anything but TERM
!
    v(bigm) = add (nonzer, term)

    do r = bigm+1, m-1
      v(r) = arbit
    end do

  else

    do r = kj+1, m-1
      v(r) = arbit
    end do

  end if
!
!  Calculate the remaining V's using the recursion of section 2.3,
!  remembering that the B's have the opposite sign.
!
  do r = 0, maxv-m
    term = 0
    do i = 0, m-1
      term = sub ( term, mul ( b(i), v(r+i) ) )
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
  real exactf
  integer i
  integer n
  real x
!
!  Integral over [0,+1]**DIMEN Product ( 1 <= i <= DIMEN )
!  abs ( 4 * x(i) - 2 ) dx
!
  if ( n == 1 ) then

    exactf = 1.0E+00
!
!  Integral over [0,+1]**DIMEN Product ( 1 <= i <= DIMEN  )
!  i * cos ( i * x(i) ) dx
!
  else if ( n == 2 ) then

    exactf = 1.0E+00
    do i = 1, dimen
      exactf = exactf * sin ( real ( i ) )
    end do
!
!  Integral over [0,+1]**DIMEN Product ( 1 <= i <= DIMEN )
!  T(N(i))( 2 * x(i) - 1 ) dx
!
!  where T(N)(X) is the Chebyshev polynomial of degree N,
!  and N(i) = mod ( i, 4 ) + 1.
!
  else if ( n == 3 ) then

    exactf = 0.0E+00
!
!  Integral over [0,+1]**DIMEN Sum ( 1 <= i <= DIMEN )
!  (-1)**i * Product ( 1 <= j <= i ) x(j) dx
!
  else if ( n == 4 ) then

    x = 1.0E+00 / real ( 2**dimen )
    if ( mod ( dimen, 2 ) == 0 ) then
      exactf = ( x - 1.0E+00 ) / 3.0E+00
    else
      exactf = ( x + 1.0E+00 ) / 3.0E+00
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
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
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
subroutine golo2 ( quasi )

!*****************************************************************************80
!
!! GOLO2 generates a new quasi-random vector on each call.
!
!  Reference:
!
!    Harald Niederreiter,
!    Low-discrepancy and low-dispersion sequences,
!    Journal of Number Theory,
!    Volume 30, 1988, pages 51-70.
!
!  Parameters:
!
!    Output, real QUASI(DIMEN), the next quasirandom vector.
!
!  Local Parameters:
!
!    Local, integer MAXDIM, the maximum dimension that will be used.
!
!    Local, integer NBITS, the number of bits (not counting the sign) in a
!    fixed-point integer.
!
!    Local, real RECIP is the multiplier which changes the
!    integers in NEXTQ into the required real values in QUASI.
!
!  Global Parameters:
!
!    In the common block /COMM2/:
!
!     CJ    - Contains the packed values of Niederreiter's C(I,J,R)
!
!     DIMEN   - The dimension of the sequence to be generated
!
!     COUNT - The index of the current item in the sequence,
!             expressed as an array of bits.  COUNT(R) is the same
!             as Niederreiter's AR(N) (page 54) except that N is
!             implicit.
!
!     NEXTQ - The numerators of the next item in the series.  These
!             are like Niederreiter's XI(N) (page 54) except that
!             N is implicit, and the NEXTQ are integers.  To obtain
!             the values of XI(N), multiply by RECIP (see GOLO2).
!
  implicit none

  integer, parameter :: maxdim = 20
  integer, parameter :: nbits = 31

  integer cj(maxdim,0:nbits-1)
  integer count
  integer dimen
  integer i
  integer nextq(maxdim)
  real quasi(*)
  integer r
  real, parameter :: recip = 2.0**(-nbits)

  common /comm2/ cj, dimen, count, nextq

  save /comm2/
!
!  Multiply the numerators in NEXTQ by RECIP to get the next
!  quasi-random vector.
!
  quasi(1:dimen) = real ( nextq(1:dimen) ) * recip
!
!  Find the position of the right-hand zero in COUNT.  This
!  is the bit that changes in the Gray-code representation as
!  we go from COUNT to COUNT+1.
!
  r = 0
  i = count

  do while ( mod ( i, 2 ) /= 0 )
    r = r + 1
    i = i / 2
  end do
!
!  Check that we have not passed 2**NBITS calls on GOLO2.
!
  if ( nbits <= r ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GOLO2 - Fatal error!'
    write ( *, '(a)' ) '  Too many calls!'
    stop
  end if
!
!  Compute the new numerators in vector NEXTQ.
!
  do i = 1, dimen
    nextq(i) = ieor ( nextq(i), cj(i,r) )
  end do

  count = count + 1

  return
end
function i_characteristic ( q )

!*****************************************************************************80
!
!! I_CHARACTERISTIC gives the characteristic for an integer.
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
!    Output, integer I_CHARACTERISTIC, the characteristic of Q.
!
  implicit none

  integer i
  integer i_characteristic
  integer i_max
  integer q
  integer q_copy

  if ( q <= 1 ) then
    i_characteristic = 0
    return
  end if
!
!  If Q is not prime, then there is at least one prime factor
!  of Q no greater than SQRT(Q)+1.
!
  i_max = int ( sqrt ( real ( q ) ) ) + 1
  q_copy = q

  do i = 2, i_max

    if ( mod ( q_copy, i ) == 0 ) then

      do while ( mod ( q_copy, i ) == 0 )
        q_copy = q_copy / i
      end do

      if ( q_copy == 1 ) then
        i_characteristic = i
      else
        i_characteristic = 0
      end if

      return

    end if

  end do
!
!  If no factor was found, then Q is prime.
!
  i_characteristic = q

  return
end
subroutine inlo2 ( dim, skip )

!*****************************************************************************80
!
!! INLO2 calculates the values of C(I,J,R).
!
!  Discussion:
!
!    This subroutine calculates the values of Niederreiter's
!    C(I,J,R) and performs other initialization necessary
!    before calling the user may call GOLO2.
!
!  Reference:
!
!    Harald Niederreiter,
!    Low-discrepancy and low-dispersion sequences,
!    Journal of Number Theory,
!    Volume 30, 1988, pages 51-70.
!
!  Parameters:
!
!    Input, integer DIM, the dimension of the sequence to be generated.
!
!    Input, integer SKIP, the number of values to throw away at the 
!    beginning of the sequence.
!
!  Local Parameters:
!
!    Local, integer MAXDIM, the maximum dimension that will be used.
!
!    Local, integer NBITS, the number of bits (not counting the sign) in a
!    fixed-point integer.
!
!  Global Parameters:
!
!    In the common block /COMM2/:
!
!    Global, integer CJ(MAXDIM,0:NBITS-1), the packed values of 
!    Niederreiter's C(I,J,R).
!
!    Global, integer DIMEN, the spatial dimension of the sequence.
!
!    Global, integer COUNT, the index of the current item in the sequence,
!    expressed as an array of bits.  COUNT(R) is the same as Niederreiter's
!    AR(N) (page 54) except that N is implicit.
!
!    Global, integer NEXTQ(MAXDIM), the numerators of the next item in the
!    series.  These are like Niederreiter's XI(N) (page 54) except that
!    N is implicit, and the NEXTQ are integers.  To obtain
!    the values of XI(N), multiply by RECIP.
!
  implicit none

  integer, parameter :: maxdim = 20
  integer, parameter :: nbits = 31

  integer cj(maxdim,0:nbits-1)
  integer count
  integer dim
  integer dimen
  integer gray
  integer i
  integer nextq(maxdim)
  integer r
  integer skip

  common /comm2/ cj, dimen, count, nextq

  save /comm2/

  dimen = dim

  if ( dimen <= 0 .or. maxdim < dimen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INLO2 - Fatal error!'
    write ( *, '(a)' ) '  Bad spatial dimension.'
    stop
  end if
!
!  Calculate the C array.
!
  call calcc2
!
!  Translate SKIP into Gray code.
!
  gray = ieor ( skip, skip / 2 )
!
!  Set up NEXTQ appropriately for this value of the Gray code.
!
  nextq(1:dimen) = 0

  r = 0

  do while ( gray /= 0 )

    if ( mod ( gray, 2 ) /= 0 ) then
      do i = 1, dimen
        nextq(i) = ieor ( nextq(i), cj(i,r) )
      end do
    end if

    gray = gray / 2
    r = r + 1

  end do

  count = skip

  return
end
subroutine plymul ( pa, pb, pc )

!*****************************************************************************80
!
!! PLYMUL multiplies two polynomials in the field of order Q.
!
!  Discussion:
!
!    Polynomials stored as arrays have the
!    coefficient of degree n in POLY(N), and the degree of the
!    polynomial in POLY(-1).  The parameter DEG is just to remind
!    us of this last fact.  A polynomial which is identically 0
!    is given degree -1.
!
!  Parameters:
!
!    Input, integer PA(-1:MAXDEG), PB(-1:MAXDEG), two polynomials
!    to be multiplied.
!
!    Output, integer PC(-1:MAXDEG), the product polynomial.
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
!
  common /field/ p, q, add, mul, sub
!
  save /field/
!
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
    write ( *, '(a)' ) '  Degree of the product exceeds MAXDEG.'
    stop
  end if

  do i = 0, degc
    term = 0
    do j = max ( 0, i-dega ), min ( degb, i )
      term = add ( term, mul ( pa(i-j), pb(j) ) )
    end do
    pt(i) = term
  end do

  pc(deg) = degc

  do i = 0, degc
    pc(i) = pt(i)
  end do

  do i = degc+1, maxdeg
    pc(i) = 0
  end do

  return
end
subroutine setfld ( qin )

!*****************************************************************************80
!
!! SETFLD sets up arithmetic tables for the finite field.
!
!  Discussion:
!
!    SETFLD sets up addition, multiplication, and subtraction tables 
!    for the finite field of order QIN.
!
!    If necessary, it reads precalculated tables from the file
!    'gfarit.dat' using unit 1.  These precalculated tables
!    are supposed to have been put there by GFARIT.
!
!  Parameters:
!
!    Input, integer QIN, the order of the field.
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
  integer file_in_unit
  integer i
  integer i_characteristic
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
    write ( *, '(a)' ) '  Bad value of Q.'
    stop
  end if

  q = qin

  p = i_characteristic ( q )

  if ( p == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETFLD - Fatal error!'
    write ( *, '(a,i6)' ) '  There is no field of order Q = ', q
    stop
  end if
!
!  Handle a field of prime order: calculate ADD and MUL.
!
  if ( p == q ) then

    do i = 0, q-1
      do j = 0, q-1
        add(i,j) = mod ( i + j, p )
        mul(i,j) = mod ( i * j, p )
      end do
    end do
!
!  Handle a field of prime-power order: tables for
!  ADD and MUL are in the file 'gfarit.dat'.
!
  else

    call get_unit ( file_in_unit )

    open ( unit = file_in_unit, file = 'gfarit.dat', status = 'old')

    do

      read ( file_in_unit, '(20i3)', iostat = ios ) n

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETFLD - Fatal error!'
        write ( *, '(a,i6)' ) '  Could not find tables for Q =', q
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
!    Input, real QUASI(DIMEN), the argument.
!
!    Output, real TESTF, the integrand at QUASI(*).
!
  implicit none

  integer dimen

  integer i
  integer ii
  integer n
  real quasi(dimen)
  real testf
  real x
!
!  Integral over [0,+1]**DIMEN Product ( 1 <= i <= DIMEN )
!  abs ( 4 * x(i) - 2 ) dx
!
  if ( n == 1 ) then

    testf = product ( abs ( 4.0E+00 * quasi(1:dimen) - 2.0E+00 ) )
!
!  Integral over [0,+1]**DIMEN Product ( 1 <= i <= DIMEN  )
!  i * cos ( i * x(i) ) dx
!
  else if ( n == 2 ) then

    testf = 1.0E+00
    do i = 1, dimen
      testf = testf * real ( i ) * cos ( real ( i ) * quasi(i) )
    end do
!
!  Integral over [0,+1]**DIMEN Product ( 1 <= i <= DIMEN )
!  T(N(i))( 2 * x(i) - 1 ) dx
!
!  where T(N)(X) is the Chebyshev polynomial of degree N,
!  and N(i) = mod ( i, 4 ) + 1.
!
  else if ( n == 3 ) then

    testf = 1.0E+00

    do i = 1, dimen

      x = 2.0E+00 * quasi(i) - 1.0E+00

      ii = mod ( i, 4 )

      if ( ii == 1 ) then
        testf = testf * x
      else if ( ii == 2 ) then
        testf = testf * ( 2.0E+00 * x**2 - 1.0E+00 )
      else if ( ii == 3 ) then
        testf = testf * ( 4.0E+00 * x**2 - 3.0E+00 ) * x
      else if ( ii == 0 ) then
        testf = testf * ( 8.0E+00 * x**4 - 8.0E+00 * x**2 + 1.0E+00 )
      end if

    end do
!
!  Integral over [0,+1]**DIMEN Sum ( 1 <= i <= DIMEN )
!  (-1)**i * Product ( 1 <= j <= i ) x(j) dx
!
  else if ( n == 4 ) then

    testf = 0.0E+00
    x = 1.0E+00
    do i = 1, dimen
      x = - x * quasi(i)
      testf = testf + x
    end do

  else

    testf = 0.0E+00

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
