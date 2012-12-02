program main

!*****************************************************************************80
!
!! MAIN is the main program for NIEDERREITER2_DATASET.
!
!  Discussion:
!
!    NIEDERREITER2_DATASET generates a Niederreiter2 dataset and writes it out.
!
!    These program assumes that your computer's word length
!    is 31 bits, excluding sign.  If this is not the case,
!    modify the parameter NBITS throughout accordingly.
!
!  Usage:
!
!    niederreiter2_dataset m n skip
!
!    where
!
!    * M, the spatial dimension,
!    * N, the number of points to generate,
!    * SKIP, number of initial values to skip.
!
!    creates an M by N dataset and writes it to the
!    file "niederreiter2_M_N.txt".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, pages 494-495, 1994.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) DIM_MAX, the maximum dimension that 
!    will be used.  
!
!    Local, integer ( kind = 4 ) POWER, is used in a possible warm-up formula.
!
  implicit none

  integer   ( kind = 4 ), parameter :: dim_max = 20

  integer   ( kind = 4 ), parameter :: base = 2
  integer   ( kind = 4 )  m
  character ( len = 255 ) output_filename
  integer   ( kind = 4 )  ios
  integer   ( kind = 4 )  n
  integer   ( kind = 4 ), parameter :: power = 12
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r
  integer   ( kind = 4 )  skip

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NIEDERREITER2_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Create a Niederreiter sequence for a given'
  write ( *, '(a)' ) '  spatial dimension, using BASE = 2.'
!
!  Get M.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the spatial dimension'
  write ( *, '(a,i6)' ) '  between 1 and ', dim_max

  read ( *, *, iostat = ios ) m

  if ( dim_max < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NIEDERREITER2_DATASET - Fatal error!'
    write ( *, '(a,i6)' ) '  The dimension may not exceed ', dim_max
    stop
  end if
!
!  Get N.
!
!  The sequence length is the number of quasi-random points used to 
!  estimate the integral, excluding warm-up.
!
!  Some users may wish to rewrite the driver to test a [heuristic] 
!  "convergence" criterion, stopping the generation of points
!  when it is passed or when a specified number of points have been 
!  generated, whichever occurs first.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Choose the sequence length:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If you do not have a preference, we'
  write ( *, '(a)' ) '  suggest using a large power of two, such as:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  2**10 = ', 2**10
  write ( *, '(a,i12)' ) '  2**15 = ', 2**15
  write ( *, '(a,i12)' ) '  2**20 = ', 2**20
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the sequence length:'

  read ( *, *, iostat = ios ) n
!
!  Read SKIP.
!
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
  write ( *, '(a)' ) '  integer ( kind = 4 ).  It does not hurt to choose '
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

  read ( *, *, iostat = ios ) skip
!
!  Compute the data.
!
  allocate ( r(1:m,1:n) )

  call niederreiter2_generate ( m, n, skip, r )
!
!  Write the data to a file.
!
  write ( output_filename, '(a,i2.2,a,i5.5,a)' ) &
    'niederreiter2_', m, '_', n, '.txt'

  call r8mat_write ( output_filename, m, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NIEDERREITER2_DATASET:'
  write ( *, '(a)' ) '  Output file name is ' // trim ( output_filename )
!
!  Free memory.
!
  deallocate ( r )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NIEDERREITER2_DATASET:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine calcc2 ( dim_num, cj )

!*****************************************************************************80
!
!! CALCC2 computes the constants C(I,J,R).
!
!  Discussion:
!
!    As far as possible, Niederreiter's notation is used.
!
!    For each value of I, we first calculate all the corresponding
!    values of C.  These are held in the array CI.  All these
!    values are either 0 or 1.  
!
!    Next we pack the values into the
!    array CJ, in such a way that CJ(I,R) holds the values of C
!    for the indicated values of I and R and for every value of
!    J from 1 to NBITS.  The most significant bit of CJ(I,R)
!    (not counting the sign bit) is C(I,1,R) and the least
!    significant bit is C(I,NBITS,R).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
!    Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Rudolf Lidl, Harald Niederreiter, 
!    Finite Fields,
!    Second Edition,
!    Cambridge University Press, 1997,
!    ISBN: 0521392314,
!    LC: QA247.3.L53
!
!    Harald Niederreiter,
!    Low-discrepancy and low-dispersion sequences,
!    Journal of Number Theory,
!    Volume 30, 1988, pages 51-70.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the sequence to 
!    be generated.
!
!    Output, integer ( kind = 4 ) CJ(DIM_MAX,0:NBITS-1), the packed values of 
!    Niederreiter's C(I,J,R)
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) MAXDEG, the highest degree of polynomial
!    to be handled. 
!
!    Local, integer ( kind = 4 ) DIM_MAX, the maximum dimension that will 
!    be used.
!
!    Local, integer ( kind = 4 ) MAXE; we need DIM_MAX irreducible polynomials 
!    over Z2.  MAXE is the highest degree among these.
!
!    Local, integer ( kind = 4 ) MAXV, the maximum possible index used in V.
!
!    Local, integer ( kind = 4 ) NBITS, the number of bits (not counting 
!    the sign) in a fixed-point integer.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxdeg = 50
  integer ( kind = 4 ), parameter :: dim_max = 20
  integer ( kind = 4 ), parameter :: maxe = 6
  integer ( kind = 4 ), parameter :: nbits = 31

  integer ( kind = 4 ), parameter :: maxv = nbits + maxe

  integer ( kind = 4 ) add(0:1,0:1)
  integer ( kind = 4 ) b(-1:maxdeg)
  integer ( kind = 4 ) ci(nbits,0:nbits-1)
  integer ( kind = 4 ) cj(dim_max,0:nbits-1)
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( dim_max, -1:maxe ) :: irred
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mul(0:1,0:1)
  integer ( kind = 4 ) px(-1:maxdeg)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) sub(0:1,0:1)
  integer ( kind = 4 ) term
  integer ( kind = 4 ) u
  integer ( kind = 4 ) v(0:maxv)
!
!  Here we supply the coefficients and the
!  degrees of the first DIM_MAX irreducible polynomials over Z2.
!
!  They are taken from Lidl and Niederreiter.
!
!  The order of these polynomials is the same as the order in
!  file 'gfplys.dat' used by the general program.
!
!  In this block PX(I, -1) is the degree of the Ith polynomial,
!  and PX(I, N) is the coefficient of x**n in the Ith polynomial.
!
  irred(1:dim_max,-1:maxe) = 0

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
  irred(13,-1:5) = (/ 5,1,1,0,1,1,1 /)
  irred(14,-1:5) = (/ 5,1,0,1,1,1,1 /)
  irred(15,-1:6) = (/ 6,1,1,0,0,0,0,1 /)
  irred(16,-1:6) = (/ 6,1,0,0,1,0,0,1 /)
  irred(17,-1:6) = (/ 6,1,1,1,0,1,0,1 /)
  irred(18,-1:6) = (/ 6,1,1,0,1,1,0,1 /)
  irred(19,-1:6) = (/ 6,1,0,0,0,0,1,1 /)
  irred(20,-1:6) = (/ 6,1,1,1,0,0,1,1 /)
!
!  Prepare to work in Z2.
!
  call setfld2 ( add, mul, sub )

  do i = 1, dim_num
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
    e = irred(i,-1)

    do j = -1, e
      px(j) = irred(i,j)
    end do

    b(-1) = 0
    b(0) = 1
!
!  Niederreiter (page 56, after equation (7), defines two
!  variables Q and U.  We do not need Q explicitly, but we do need U.
!
    u = 0

    do j = 1, nbits
!
!  If U = 0, we need to set B to the next power of PX
!  and recalculate V.  This is done by subroutine CALCV2.
!
      if ( u == 0 ) then
        call calcv2 ( maxv, px, add, mul, sub, b, v )
      end if
!
!  Now C is obtained from V.  Niederreiter obtains A from V (page 65, 
!  near the bottom), and then gets C from A (page 56, equation (7)).  
!  However this can be done in one step.  Here CI(J,R) corresponds to
!  Niederreiter's C(I,J,R).
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
subroutine calcv2 ( maxv, px, add, mul, sub, b, v )

!*****************************************************************************80
!
!! CALCV2 calculates the constants V(J,R).
!
!  Discussion:
!
!    This program calculates the values of the constants V(J,R) as
!    described in the reference (BFN) section 3.3.  It is called from 
!    either CALCC or CALCC2.  
!
!    Polynomials stored as arrays have the coefficient of degree N 
!    in POLY(N), and the degree of the polynomial in POLY(-1).  
!
!    A polynomial which is identically 0 is given degree -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
!    Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, 1994, pages 494-495.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXV, the dimension of the array V.
!
!    Input, integer ( kind = 4 ) PX(-1:MAXDEG), the appropriate irreducible 
!    polynomial for the dimension currently being considered.  The degree of PX 
!    will be called E.
!
!    Input, integer ( kind = 4 ) ADD(0:1,0:1), MUL(0:1,0:1), SUB(0:1,0:1), the 
!    addition, multiplication, and subtraction tables, mod 2.
!
!    Input/output, integer ( kind = 4 ) B(-1:MAXDEG).  On input, B is the 
!    polynomial defined in section 2.3 of BFN.  The degree of B implicitly 
!    defines the parameter J of section 3.3, by degree(B) = E*(J-1).  On output,
!    B has been multiplied by PX, so its degree is now E * J.
!
!    Output, integer ( kind = 4 ) V(0:MAXV), the computed V array.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) ARBIT, indicates where the user can place
!    an arbitrary element of the field of order 2.  This means 
!    0 <= ARBIT < 2.  
!
!    Local, integer ( kind = 4 ) BIGM, is the M used in section 3.3.
!    It differs from the [little] m used in section 2.3,
!    denoted here by M.
!
!    Local, integer ( kind = 4 ) MAXDEG, the highest degree of polynomial
!    to be handled. 
!
!    Local, integer ( kind = 4 ) NONZER, shows where the user must put an 
!    arbitrary non-zero element of the field.  For the code, this means 
!    0 < NONZER < 2.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxdeg = 50

  integer ( kind = 4 ) add(0:1,0:1)
  integer ( kind = 4 ), parameter :: arbit = 1
  integer ( kind = 4 ) b(-1:maxdeg)
  integer ( kind = 4 ) bigm
  integer ( kind = 4 ) e
  integer ( kind = 4 ) h(-1:maxdeg)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kj
  integer ( kind = 4 ) m
  integer ( kind = 4 ) maxv
  integer ( kind = 4 ) mul(0:1,0:1)
  integer ( kind = 4 ), parameter :: nonzer = 1
  integer ( kind = 4 ) p
  integer ( kind = 4 ) px(-1:maxdeg)
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  integer ( kind = 4 ) sub(0:1,0:1)
  integer ( kind = 4 ) term
  integer ( kind = 4 ) v(0:maxv)

  p = 2
  q = 2
  e = px(-1)
!
!  The polynomial H is PX**(J-1), which is the value of B on arrival.
!
!  In section 3.3, the values of Hi are defined with a minus sign:
!  don't forget this if you use them later!
!
  do i = -1, b(-1)
    h(i) = b(i)
  end do

  bigm = h(-1)
!
!  Multiply B by PX so B becomes PX**J.
!  In section 2.3, the values of Bi are defined with a minus sign:
!  don't forget this if you use them later!
!
  call plymul2 ( add, mul, px, b, b )
  m = b(-1)
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
!  Choose values of V in accordance with the conditions in section 3.3.
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
!  Now V(BIGM) is anything but TERM.
!
    v(bigm) = add ( nonzer, term )

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
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical              lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

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
subroutine niederreiter2 ( dim_num, seed, quasi )

!*****************************************************************************80
!
!! NIEDERREITER2 returns an element of the Niederreiter sequence base 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
!    Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
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
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the sequence to 
!    be generated.
!
!    Input/output, integer ( kind = 4 ) SEED, the index of the element entry to
!    compute.  On output, SEED is typically reset by this routine
!    to SEED+1.
!
!    Output, real ( kind = 8 ) QUASI(DIM_NUM), the next quasirandom vector.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) DIM_MAX, the maximum dimension that will 
!    be used.
!
!    Local, integer ( kind = 4 ) NBITS, the number of bits (not counting the 
!    sign) in a fixed-point integer.
!
!    Local, real ( kind = 8 ) RECIP, the multiplier which changes the
!    integers in NEXTQ into the required real values in QUASI.
!
!    Local, integer ( kind = 4 ) CJ(DIM_MAX,0:NBITS-1), the packed values of 
!    Niederreiter's C(I,J,R).
!
!    Local, integer ( kind = 4 ) DIM_SAVE, the spatial dimension of the sequence
!    as specified on an initialization call.
!
!    Local, integer ( kind = 4 ) NEXTQ(DIM_MAX), the numerators of the next item 
!    in the series.  These are like Niederreiter's XI(N) (page 54) except that
!    N is implicit, and the NEXTQ are integers.  To obtain
!    the values of XI(N), multiply by RECIP.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), parameter :: dim_max = 20
  integer ( kind = 4 ), parameter :: nbits = 31

  integer ( kind = 4 ), save, dimension (dim_max,0:nbits-1) :: cj
  integer ( kind = 4 ), save :: dim_save = 0
  integer ( kind = 4 ) gray
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save, dimension ( dim_max ) :: nextq
  real ( kind = 8 ) quasi(dim_num)
  integer ( kind = 4 ) r
  real ( kind = 8 ), parameter :: recip = 2.0D+00**(-nbits)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: seed_save = 0
!
!  Initialization.
!
  if ( dim_save < 1 .or. dim_num /= dim_save .or. seed <= 0 ) then

    if ( dim_num <= 0 .or. dim_max < dim_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NIEDERREITER2 - Fatal error!'
      write ( *, '(a)' ) '  Bad spatial dimension.'
      stop
    end if

    dim_save = dim_num

    if ( seed < 0 ) then
      seed = 0
    end if

    seed_save = seed
!
!  Calculate the C array.
!
    call calcc2 ( dim_save, cj )

  end if
!
!  Set up NEXTQ appropriately, depending on the Gray code of SEED.
!
!  You can do this every time, starting NEXTQ back at 0,
!  or you can do it once, and then carry the value of NEXTQ
!  around from the previous computation.
!
  if ( seed /= seed_save + 1 ) then

    gray = ieor ( seed, seed / 2 )

    nextq(1:dim_save) = 0

    r = 0

    do while ( gray /= 0 )

      if ( mod ( gray, 2 ) /= 0 ) then
        do i = 1, dim_save
          nextq(i) = ieor ( nextq(i), cj(i,r) )
        end do
      end if

      gray = gray / 2
      r = r + 1

    end do

  end if
!
!  Multiply the numerators in NEXTQ by RECIP to get the next
!  quasi-random vector.
!
  quasi(1:dim_save) = real ( nextq(1:dim_save), kind = 8 ) * recip
!
!  Find the position of the right-hand zero in SEED.  This
!  is the bit that changes in the Gray-code representation as
!  we go from SEED to SEED+1.
!
  r = 0
  i = seed

  do while ( mod ( i, 2 ) /= 0 )
    r = r + 1
    i = i / 2
  end do
!
!  Check that we have not passed 2**NBITS calls.
!
  if ( nbits <= r ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NIEDERREITER2 - Fatal error!'
    write ( *, '(a)' ) '  Too many calls!'
    stop
  end if
!
!  Compute the new numerators in vector NEXTQ.
!
  do i = 1, dim_save
    nextq(i) = ieor ( nextq(i), cj(i,r) )
  end do

  seed_save = seed
  seed = seed + 1

  return
end
subroutine niederreiter2_generate ( dim_num, n, seed, r )

!*****************************************************************************80
!
!! NIEDERREITER2_GENERATE generates a set of Niederreiter values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points desired.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed

  do j = 1, n
    call niederreiter2 ( dim_num, seed, r(1:dim_num,j) )
  end do

  return
end
subroutine plymul2 ( add, mul, pa, pb, pc )

!*****************************************************************************80
!
!! PLYMUL2 multiplies two polynomials in the field of order 2.
!
!  Discussion:
!
!    Polynomials stored as arrays have the coefficient of degree N in 
!    POLY(N), and the degree of the polynomial in POLY(-1).  
!
!    A polynomial which is identically 0 is given degree -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
!    Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, 1994, pages 494-495.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADD(0:1,0:1), MUL(0:1,0:1),
!    the addition and multiplication tables, mod 2.
!
!    Input, integer ( kind = 4 ) PA(-1:MAXDEG), PB(-1:MAXDEG), two polynomials
!    to be multiplied.
!
!    Output, integer ( kind = 4 ) PC(-1:MAXDEG), the product polynomial.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) MAXDEG, the highest degree of polynomial
!    to be handled. 
!
  implicit none

  integer ( kind = 4 ), parameter :: maxdeg = 50

  integer ( kind = 4 ) add(0:1,0:1)
  integer ( kind = 4 ) dega
  integer ( kind = 4 ) degb
  integer ( kind = 4 ) degc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mul(0:1,0:1)
  integer ( kind = 4 ) pa(-1:maxdeg)
  integer ( kind = 4 ) pb(-1:maxdeg)
  integer ( kind = 4 ) pc(-1:maxdeg)
  integer ( kind = 4 ) pt(-1:maxdeg)
  integer ( kind = 4 ) term

  dega = pa(-1)
  degb = pb(-1)

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

  pc(-1) = degc

  do i = 0, degc
    pc(i) = pt(i)
  end do

  do i = degc+1, maxdeg
    pc(i) = 0
  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_filename
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine setfld2 ( add, mul, sub )

!*****************************************************************************80
!
!! SETFLD2 sets up arithmetic tables for the finite field of order 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
!    Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, 1994, pages 494-495.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADD(0:1,0:1), MUL(0:1,0:1), SUB(0:1,0:1), the 
!    addition, multiplication, and subtraction tables, mod 2.
!
  implicit none

  integer ( kind = 4 ) add(0:1,0:1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mul(0:1,0:1)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) sub(0:1,0:1)

  q = 2

  p = 2

  do i = 0, q-1
    do j = 0, q-1
      add(i,j) = mod ( i + j, p )
    end do
  end do

  do i = 0, q-1
    do j = 0, q-1
      mul(i,j) = mod ( i * j, p )
    end do
  end do

  do i = 0, q-1
    do j = 0, q-1
      sub(add(i,j), i) = j
    end do
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
