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
!    Polynomials stored as arrays have the coefficient of degree N 
!    in POLY(N), and the degree of the polynomial in POLY(-1).  
!    The parameter DEG is just to remind us of this last fact.  
!    A polynomial which is identically 0 is given degree -1.
!
!    Thanks to Michael Baudin for pointing out that MAXE should
!    be increased from 5 to 7, since one of the irreducible polynomials
!    that must be stored in PX has degree 7, 07 June 2010.
!
!    In fact, the size of MAXE depends on the highest degree polynomial
!    computed by GFPLYS, which in turn depends in part on the 
!    value NPOL in that program.  To allow DIM_MAX = 50, we increased
!    NPOL to 50 in GFPLYS and here, and hence MAXE to 8, JVB, 
!    and DIM_MAX from 20 to 50 in CALCC, GOLO and INLO.  07 June 2010.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 June 2010
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
!    Harald Niederreiter,
!    Low-discrepancy and low-dispersion sequences,
!    Journal of Number Theory,
!    Volume 30, 1988, pages 51-70.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) DEG_MAX, the highest degree of polynomial
!    to be handled.  
!
!    Local, integer ( kind = 4 ) DIM_MAX, the maximum dimension that will 
!    be used.
!
!    Local, integer ( kind = 4 ) MAXE; we need DIM_MAX irreducible polynomials 
!    over GF(Q).  MAXE is the highest degree among these.
!
!    Local, integer ( kind = 4 ) FIG_MAX, the maximum number of base-Q digits 
!    we can handle.  BASE**FIG_MAX - 1 must be representable in FORTRAN.
!    For base 2, this implies that FIG_MAX could be as high as 31.
!    In the original version of the program, FIG_MAX was set to 20.
!
!    Local, integer ( kind = 4 ) Q_MAX, the order of the largest field 
!    to be handled.
!
!    Local, integer ( kind = 4 ) V_MAX, the maximum index used in V.
!
!    Local, integer ( kind = 4 ) NPOLS, the number of precalculated 
!    irreducible polynomials.
!
!  Global parameters in /COMM/:
!
!    Global, integer ( kind = 4 ) C(DIM_MAX,FIG_MAX,0:FIG_MAX-1), the values of 
!    Niederreiter's C(I,J,R).
!
!    Global, integer ( kind = 4 ) COUNT(0:FIG_MAX-1), the index of the current 
!    item in the sequence, expressed as an array of base-Q digits.  COUNT(R)
!    is the same as Niederreiter's AR(N) (page 54) except that N is implicit.
!
!    Global, integer ( kind = 4 ) D(DIM_MAX,FIG_MAX).
!
!    Global, integer ( kind = 4 ) NEXTQ(DIM_MAX), the numerators of the next
!    item in the series.  These are like Niederreiter's XI(N) (page 54) except 
!    that N is implicit, and the NEXTQ are integers.  To obtain the values of 
!    XI(N), multiply by RECIP.
!
!    Global, integer ( kind = 4 ) QPOW(FIG_MAX), to speed things up a bit. 
!    QPOW(I) = Q ** (NFIGS-I).
!
!    Global, integer ( kind = 4 ) DIMEN, the dimension of the sequence to be 
!    generated.
!
!    Global, integer ( kind = 4 ) NFIGS, the number of base Q digits we 
!    are using.
!
!    Global, real ( kind = 8 ) RECIP = 1.0 / Q^NFIGS.
!
  implicit none

  integer ( kind = 4 ), parameter :: deg = -1
  integer ( kind = 4 ), parameter :: deg_max = 50
  integer ( kind = 4 ), parameter :: dim_max = 50
! integer ( kind = 4 ), parameter :: maxe = 5
  integer ( kind = 4 ), parameter :: maxe = 8
! integer ( kind = 4 ), parameter :: fig_max = 20
  integer ( kind = 4 ), parameter :: fig_max = 31
  integer ( kind = 4 ), parameter :: q_max = 50
! integer ( kind = 4 ), parameter :: npols = 25
  integer ( kind = 4 ), parameter :: npols = 50

  integer ( kind = 4 ), parameter :: v_max = fig_max + maxe

  integer ( kind = 4 ) add(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) b(-1:deg_max)
  integer ( kind = 4 ) c(dim_max,fig_max,0:fig_max-1)
  integer ( kind = 4 ) count(0:fig_max-1)
  integer ( kind = 4 ) d(dim_max,fig_max)
  integer ( kind = 4 ) dimen
  integer ( kind = 4 ) e
  character ( len = 80 ) :: file_in_name = 'gfplys.txt'
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mul(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) nextq(dim_max)
  integer ( kind = 4 ) nfigs
  integer ( kind = 4 ) p
  integer ( kind = 4 ) px(-1:maxe)
  integer ( kind = 4 ) q
  integer ( kind = 4 ) qpow(fig_max)
  integer ( kind = 4 ) r
  real ( kind = 8 ) recip
  integer ( kind = 4 ) sub(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) u
  integer ( kind = 4 ) v(0:v_max)

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
!  and its degree into E.  Set polynomial B = PX^0 = 1.
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
        call calcv ( px, b, v, v_max )
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
subroutine calcv ( px, b, v, v_max )

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 2007
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
!    Input, integer ( kind = 4 ) PX(-1:DEG_MAX), the appropriate irreducible 
!    polynomial for the dimension currently being considered.  The degree of 
!    PX will be called E.
!
!    Input/output, integer ( kind = 4 ) B(-1:DEG_MAX).  On input, B is the 
!    polynomial defined in section 2.3 of the BFN paper.  The degree of B 
!    implicitly defines the parameter J of section 3.3, by degree(B) = E*(J-1).  
!    On output, B has been multiplied by PX, so its degree is now E*J.
!
!    Input, integer ( kind = 4 ) V(0:V_MAX), contains the values required.
!
!    Input, integer ( kind = 4 ) V_MAX, the dimension of the array V.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) ARBIT, indicates where the user can place
!    an arbitrary element of the field of order Q.  For the code,
!    this means 0 <= ARBIT < Q.  Within these limits, the user can 
!    do what he likes.  ARBIT could be declared as a function 
!    that returned a different arbitrary value each time it is referenced.
!
!    Local, integer ( kind = 4 ) BIGM, is the M used in section 3.3.  It differs
!    from the [little] m used in section 2.3, denoted here by M.
!
!    Local, integer ( kind = 4 ) DEG_MAX, the highest degree of polynomial to be
!    handled. 
! 
!    Local, integer ( kind = 4 ) Q_MAX, the order of the largest field to 
!    be handled.
!
!    Local, integer ( kind = 4 ) NONZER shows where the user must put an 
!    arbitrary non-zero element of the same field.  For the code, this means 
!    0 < NONZER < Q.  Within these limits, the user's choice is free.
!
  implicit none

  integer ( kind = 4 ), parameter :: deg = -1
  integer ( kind = 4 ), parameter :: deg_max = 50
  integer ( kind = 4 ), parameter :: q_max = 50

  integer ( kind = 4 ) add(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ), parameter :: arbit = 1
  integer ( kind = 4 ) b(-1:deg_max)
  integer ( kind = 4 ) bigm
  integer ( kind = 4 ) e
  integer ( kind = 4 ) h(-1:deg_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kj
  integer ( kind = 4 ) m
  integer ( kind = 4 ) v_max
  integer ( kind = 4 ) mul(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ), parameter :: nonzer = 1
  integer ( kind = 4 ) p
  integer ( kind = 4 ) px(-1:deg_max)
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  integer ( kind = 4 ) sub(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) term
  integer ( kind = 4 ) v(0:v_max)

  common /field/ p, q, add, mul, sub

  save /field/

  e = px(deg)
!
!  The polynomial H is PX^(J-1), which is the value of B on arrival.
!
!  In section 3.3, the values of Hi are defined with a minus sign:
!  don't forget this if you use them later!
!
  do i = -1, b(deg)
    h(i) = b(i)
  end do

  bigm = h(deg)
!
!  Now multiply B by PX so B becomes PX^J.
!
!  In section 2.3, the values of Bi are defined with a minus sign:
!  don't forget this if you use them later!
!
  call plymul ( px, h, b )
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
  do r = 0, v_max - m
    term = 0
    do i = 0, m - 1
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
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Output, integer ( kind = 4 ) IUNIT.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2009
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
!    Local, integer ( kind = 4 ) DEG_MAX, the highest degree of polynomial 
!    to be handled. 
! 
!    Local, integer ( kind = 4 ) DIM_MAX, the maximum dimension that will 
!    be used.
!
!    Local, integer ( kind = 4 ) FIG_MAX, the maximum number of base-Q digits 
!    we can handle.  BASE**FIG_MAX - 1 must be representable in FORTRAN.
!    For base 2, this implies that FIG_MAX will be 31.
!    In the original version of the program, FIG_MAX was set to 20.
!
!    Local, integer ( kind = 4 ) Q_MAX, the order of the largest field to 
!    be handled.
!
!    Local, integer ( kind = 4 ) NBITS, the number of bits in a fixed-point 
!    integer, not counting the sign.
!
!  Global parameters in /COMM/:
!
!    Global, integer ( kind = 4 ) C(DIM_MAX,FIG_MAX,0:FIG_MAX-1), the values of 
!    Niederreiter's C(I,J,R).
!
!    Global, integer ( kind = 4 ) COUNT(0:FIG_MAX-1), the index of the current 
!    item in the sequence, expressed as an array of base-Q digits.  COUNT(R)
!    is the same as Niederreiter's AR(N) (page 54) except that N is implicit.
!
!    Global, integer ( kind = 4 ) D(DIM_MAX,FIG_MAX).
!
!    Global, integer ( kind = 4 ) NEXTQ(DIM_MAX), the numerators of the next 
!    item in the series.  These are like Niederreiter's XI(N) (page 54) except 
!    that N is implicit, and the NEXTQ are integers.  To obtain the values of 
!    XI(N), multiply by RECIP.
!
!    Global, integer ( kind = 4 ) QPOW(FIG_MAX), to speed things up a bit. 
!    QPOW(I) = Q^(NFIGS-I).
!
!    Global, integer ( kind = 4 ) DIMEN, the dimension of the sequence to 
!    be generated.
!
!    Global, integer ( kind = 4 ) NFIGS, the number of base Q digits we 
!    are using.
!
!    Global, real ( kind = 8 ) RECIP = 1.0 / Q^NFIGS.
!
  implicit none

  integer ( kind = 4 ), parameter :: deg = -1
  integer ( kind = 4 ), parameter :: deg_max = 50
! integer ( kind = 4 ), parameter :: dim_max = 20
  integer ( kind = 4 ), parameter :: dim_max = 50
! integer ( kind = 4 ), parameter :: fig_max = 20
  integer ( kind = 4 ), parameter :: fig_max = 31
  integer ( kind = 4 ), parameter :: q_max = 50

  integer ( kind = 4 ) add(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) c(dim_max,fig_max,0:fig_max-1)
  integer ( kind = 4 ) count(0:fig_max-1)
  integer ( kind = 4 ) d(dim_max, fig_max)
  integer ( kind = 4 ) diff
  integer ( kind = 4 ) dimen
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mul(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) nextq(dim_max)
  integer ( kind = 4 ) nfigs
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) oldcnt
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) qpow(fig_max)
  real ( kind = 8 ) quasi(*)
  integer ( kind = 4 ) r
  real ( kind = 8 ) recip
  integer ( kind = 4 ) sub(0:q_max-1,0:q_max-1)

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
      write ( *, '(a)' ) '  NFIGS <= R.'
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
!  QPOW(J) is Q^(NFIGS-J).
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
!    P, if Q = P^N for some prime P and some integer N;
!    0, otherwise, that is, if Q is negative, 0, 1, or the product
!       of more than one distinct prime.
!
!    A faster code would only consider prime values of I,
!    but that entails storing a table of primes and limiting the
!    size of Q.  Simplicity and flexibility for now!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 2007
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
!    Input, integer ( kind = 4 ) Q, the value to be tested.
!
!    Output, integer ( kind = 4 ) I4_CHARACTERISTIC, the characteristic of Q.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_characteristic
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) q
  integer ( kind = 4 ) q_copy

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2009
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
!    Harald Niederreiter,
!    Low-discrepancy and low-dispersion sequences,
!    Journal of Number Theory,
!    Volume 30, 1988, pages 51-70.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM, the dimension of the sequence to be 
!    generated.  The value of DIM is copied into DIMEN in the common block.
!
!    Input, integer ( kind = 4 ) BASE, the prime or prime-power base to be used.
!
!    Input, integer ( kind = 4 ) SKIP, the number of values to throw away at 
!    the beginning of the sequence.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) Q_MAX, the order of the largest field to be 
!    handled.
!
!    Local, integer ( kind = 4 ) DEG_MAX, the highest degree of polynomial
!    to be handled. 
!
!    Local, integer ( kind = 4 ) DIM_MAX, the maximum dimension that will 
!    be used.
!
!    Local, integer ( kind = 4 ) FIG_MAX, the maximum number of base-Q digits 
!    we can handle.  BASE**FIG_MAX - 1 must be representable in FORTRAN.
!    For base 2, this implies that FIG_MAX will be 31.
!    In the original version of the program, FIG_MAX was set to 20.
!
!    Local, integer ( kind = 4 ) NBITS, the number of bits in a fixed-point 
!    integer, not counting the sign.
!
!  Global parameters in /COMM/:
!
!    Global, integer ( kind = 4 ) C(DIM_MAX,FIG_MAX,0:FIG_MAX-1), the values of 
!    Niederreiter's C(I,J,R).
!
!    Global, integer ( kind = 4 ) COUNT(0:FIG_MAX-1), the index of the current 
!    item in the sequence, expressed as an array of base-Q digits.  COUNT(R)
!    is the same as Niederreiter's AR(N) (page 54) except that N is implicit.
!
!    Global, integer ( kind = 4 ) D(DIM_MAX,FIG_MAX).
!
!    Global, integer ( kind = 4 ) NEXTQ(DIM_MAX), the numerators of the next 
!    item in the series.  These are like Niederreiter's XI(N) (page 54) except 
!    that N is implicit, and the NEXTQ are integers.  To obtain the values of 
!    XI(N), multiply by RECIP.
!
!    Global, integer ( kind = 4 ) QPOW(FIG_MAX), to speed things up a bit. 
!    QPOW(I) = Q^(NFIGS-I).
!
!    Global, integer ( kind = 4 ) DIMEN, the dimension of the sequence to 
!    be generated.
!
!    Global, integer ( kind = 4 ) NFIGS, the number of base Q digits we 
!    are using.
!
!    Global, real ( kind = 8 ) RECIP = 1.0 / Q^NFIGS.
!
  implicit none

  integer ( kind = 4 ), parameter :: deg = -1
  integer ( kind = 4 ), parameter :: deg_max = 50
! integer ( kind = 4 ), parameter :: dim_max = 20
  integer ( kind = 4 ), parameter :: dim_max = 50
! integer ( kind = 4 ), parameter :: fig_max = 20
  integer ( kind = 4 ), parameter :: fig_max = 31
  integer ( kind = 4 ), parameter :: q_max = 50
  integer ( kind = 4 ), parameter :: nbits = 31

  integer ( kind = 4 ) add(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) c(dim_max, fig_max,0:fig_max-1)
  integer ( kind = 4 ) count(0:fig_max-1)
  integer ( kind = 4 ) d(dim_max, fig_max)
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dimen
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_characteristic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mul(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) nextq(dim_max)
  integer ( kind = 4 ) nfigs
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) qpow(fig_max)
  integer ( kind = 4 ) r
  real ( kind = 8 ) recip
  integer ( kind = 4 ) skip
  integer ( kind = 4 ) sub(0:q_max-1,0:q_max-1)
  real ( kind = 8 ) temp

  common /comm/ c, count, d, nextq, qpow, dimen, nfigs, recip
  common /field/ p, q, add, mul, sub

  save /comm/
  save /field/

  dimen = dim

  if ( dimen <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INLO - Fatal error!'
    write ( *, '(a)' ) '  DIM <= 0.'
    stop
  end if

  if ( dim_max < dimen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INLO - Fatal error!'
    write ( *, '(a)' ) '  DIM_MAX < DIM.'
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

  nfigs = min ( fig_max, int ( temp ) )
!
!  Calculate the C array.
!
  call calcc ( )
!
!  Set RECIP.
!
  recip = 1.0D+00 / real ( q**nfigs, kind = 8 )
!
!  Set QPOW(I) = Q^(NFIGS-I).
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
subroutine niederreiter ( dim_num, base, seed, r )

!*****************************************************************************80
!
!! NIEDERREITER returns an element of a Niederreiter sequence for base BASE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) BASE, the base to use for the Niederreiter 
!    sequence.  The base should be a prime, or a power of a prime.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) R(DIM_NUM), the element of the sequence.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) base
  integer ( kind = 4 ), save :: dim_num_save = -1
  real ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) skip

  if ( dim_num_save < 1 .or. dim_num /= dim_num_save .or. seed <= 0 ) then
!
!  It seems as though INLO wants the "SKIP" argument to be nonzero.
!
    skip = 1

    call inlo ( dim_num, base, skip )

    dim_num_save = dim_num

  end if

  call golo ( r )

  seed = seed + 1

  return
end
subroutine niederreiter_generate ( dim_num, n, base, seed, r )

!*****************************************************************************80
!
!! NIEDERREITER_GENERATE generates a set of Niederreiter values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2007
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
!    Input, integer ( kind = 4 ) BASE, the base to use for the Niederreiter 
!    sequence.  The base should be a prime, or a power of a prime.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) R(DIM_NUM,N), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(dim_num,n)
  integer ( kind = 4 ) seed

  do j = 1, n
    call niederreiter ( dim_num, base, seed, r(1:dim_num,j) )
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 2007
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
!    Input, integer ( kind = 4 ) PA(-1:DEG_MAX), PB(-1:DEG_MAX), two polynomials
!    to be multiplied.
!
!    Output, integer ( kind = 4 ) PC(-1:DEG_MAX), the product of PA and PB, 
!    in the field Q.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) DEG_MAX, the highest degree of polynomial
!    to be handled. 
!
!    Local, integer ( kind = 4 ) Q_MAX, the order of the largest field to 
!    be handled.
!
  implicit none

  integer ( kind = 4 ), parameter :: deg = -1
  integer ( kind = 4 ), parameter :: deg_max = 50
  integer ( kind = 4 ), parameter :: q_max = 50

  integer ( kind = 4 ) add(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) dega
  integer ( kind = 4 ) degb
  integer ( kind = 4 ) degc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mul(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pa(-1:deg_max)
  integer ( kind = 4 ) pb(-1:deg_max)
  integer ( kind = 4 ) pc(-1:deg_max)
  integer ( kind = 4 ) pt(-1:deg_max)
  integer ( kind = 4 ) q
  integer ( kind = 4 ) sub(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) term

  common /field/ p, q, add, mul, sub

  save /field/ 

  dega = pa(deg)
  degb = pb(deg)

  if ( dega == -1 .or. degb == -1 ) then
    degc = -1
  else
    degc = dega + degb
  end if

  if ( deg_max < degc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLYMUL - Fatal error!'
    write ( *, '(a)' ) '  The degree of the product exceeds DEG_MAX.'
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

  do i = degc + 1, deg_max
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

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 2007
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
!    Input, integer ( kind = 4 ) QIN, the order of the field.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) DEG_MAX, the highest degree of polynomial to 
!    be handled.  
!
!    Local, integer ( kind = 4 ) Q_MAX, the order of the largest field to 
!    be handled.
!
  implicit none

  integer ( kind = 4 ), parameter :: deg = -1
  integer ( kind = 4 ), parameter :: deg_max = 50
  integer ( kind = 4 ), parameter :: q_max = 50

  integer ( kind = 4 ) add(0:q_max-1,0:q_max-1)
  character ( len = 80 ) :: file_in_name = 'gfarit.txt'
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_characteristic
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mul(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) qin
  integer ( kind = 4 ) sub(0:q_max-1,0:q_max-1)

  common /field/ p, q, add, mul, sub

  save /field/

  if ( qin <= 1 .or. q_max < qin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETFLD - Fatal error!'
    write ( *, '(a)' ) '  Bad value of Q!'
    write ( *, '(a)' ) '  Q <= 1 or Q_MAX < Q.'
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
    do i = 0, q - 1
      do j = 0, q - 1
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

      do i = 0, n - 1
        read ( file_in_unit, '(20i3)' ) add(i,0:n-1)
      end do

      do i = 0, n - 1
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
  do i = 0, q - 1
    do j = 0, q - 1
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

  character ( len = 8 )  ampm
  integer ( kind = 4 ) d
  character ( len = 8 )  date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 )  zone

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
