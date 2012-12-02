program main

!*****************************************************************************80
!
!! MAIN is the main program for GFPLYS.
!
!  Discussion:
!
!    GFPLYS writes out data about irreducible polynomials.
!
!    The program calculates irreducible polynomials for various
!    finite fields, and writes them out to the file "gfplys.txt".
!
!    Finite field arithmetic is carried out with the help of
!    precalculated addition and multiplication tables found on
!    the file "gfarit.txt".  This file should have been computed
!    and written by the program GFARIT.
!
!    The format of the irreducible polynomials on the output file is
!
!      Q
!      d1   a(1)  a(2) ... a(d1)
!      d2   b(1)  b(2) ... b(d2)
!      ...
!
!    where 
!
!      Q is the order of the field, 
!      d1 is the degree of the first irreducible polynomial, 
!      a(1), a(2), ..., a(d1) are its coefficients.
!
!    Polynomials stored as arrays have the coefficient of degree N in 
!    POLY(N), and the degree of the polynomial in POLY(-1).  The parameter
!    DEG is just to remind us of this last fact.  A polynomial which is
!    identically 0 is given degree -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2007
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
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) DEG_MAX, the highest degree of polynomial 
!    to be handled.  
!
!    Local, integer ( kind = 4 ) Q_MAX, the order of the largest field to 
!    be handled.
!
!  Global Parameters:
!
!    The common block /FIELD/ contains:
!
!    Common, integer ( kind = 4 ) Q, the order of the field.
!
!    Common, integer ( kind = 4 ) P, the characteristic of the field.
!
!    Common, integer ( kind = 4 ) ADD(0:Q_MAX-1,0:Q_MAX-1), the field addition table. 
!
!    Common, integer ( kind = 4 ) MUL(0:Q_MAX-1,0:Q_MAX-1), the field multiplication table. 
!
!    Common, integer ( kind = 4 ) SUB(0:Q_MAX-1,0:Q_MAX-1), the field 
!    subtraction table.
!
  implicit none

  integer ( kind = 4 ), parameter :: deg = -1
  integer ( kind = 4 ), parameter :: deg_max = 50
  integer ( kind = 4 ), parameter :: q_max = 50

  integer ( kind = 4 ) add(0:q_max-1,0:q_max-1)
  character ( len = 80 ) :: file_out_name = 'gfplys.txt'
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) mul(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) q_init
  integer ( kind = 4 ) sub(0:q_max-1,0:q_max-1)

  common /field/ p, q, add, mul, sub

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GFPLYS:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A program to compute a set of irreducible'
  write ( *, '(a)' ) '  polynomials over fields of certain orders Q.'

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace' )

  do q_init = 2, q_max
    call irred ( file_out_unit, q_init )
  end do

  close ( unit = file_out_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GFPLYS:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function find ( n, tab, i, tab_max )

!*****************************************************************************80
!
!! FIND seeks the value N in the range TAB(I) to TAB(TAB_MAX).
!
!  Discussion:
!
!    The vector TAB does not have to be sorted or have any other
!    special properties.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2007
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
!    Input, integer ( kind = 4 ) N, the value being sought.
!
!    Input, integer ( kind = 4 ) TAB(*), the table to be searched.
!
!    Input, integer ( kind = 4 ) I, TAB_MAX, the first and last entries of
!    TAB to be examined.
!
!    Output, integer ( kind = 4 ) FIND, is the index ( between I and TAB_MAX) of the 
!    entry in TAB that is equal to N, or else -1 if no such value
!    was found.
!
  implicit none

  integer ( kind = 4 ) tab_max

  integer ( kind = 4 ) find
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) tab(tab_max)

  find = -1

  if ( tab(tab_max) < n ) then
    return
  end if

  do j = i, tab_max
    if ( tab(j) == n ) then
      find = j
      return
    end if
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2007
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
!  A faster code would only consider prime values of I,
!  but that entails storing a table of primes and limiting the
!  size of Q.  Simplicity and flexibility for now!
!
  i_max = int ( sqrt ( real ( q ) ) ) + 1
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
subroutine irred ( ounit, q_init )

!*****************************************************************************80
!
!! IRRED computes and writes out a set of irreducible polynomials.
!
!  Discussion:
!
!    We find the irreducible polynomials using a sieve.  
!
!    Polynomials stored as arrays have the coefficient of degree n in 
!    POLY(N), and the degree of the polynomial in POLY(-1).  The parameter
!    DEG is just to remind us of this last fact.  A polynomial which is
!    identically 0 is given degree -1.
!
!    Note that the value of NPOL controls the number of polynomials
!    computed, and hence the maximum spatial dimension for the
!    subsequence Niederreiter sequences, JVB, 07 June 2010.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUNIT, the unit number of the output file.
!
!    Input, integer ( kind = 4 ) Q_INIT, the order of the field.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) SIEVE_MAX, the size of the sieve.  
!
!    Array MONPOL holds monic polynomials.
!
!    Array SIEVE says whether the polynomial is still OK.
!
!    Local, integer ( kind = 4 ) DEG_MAX, the highest degree of polynomial
!    to be handled.  
!
!    Local, integer ( kind = 4 ) Q_MAX, the order of the largest field to
!    be handled.
!
!    Local, integer ( kind = 4 ) NPOLS, the number of irreducible polynomials to
!    be calculated for a given field.
!
!  Global Parameters:
!
!    The common block /FIELD/ contains:
!
!    Common, integer ( kind = 4 ) Q, the order of the field.
!
!    Common, integer ( kind = 4 ) P, the characteristic of the field.
!
!    Common, integer ( kind = 4 ) ADD(0:Q_MAX-1,0:Q_MAX-1), the field 
!    addition table. 
!
!    Common, integer ( kind = 4 ) MUL(0:Q_MAX-1,0:Q_MAX-1), the field 
!    multiplication table. 
!
!    Common, integer ( kind = 4 ) SUB(0:Q_MAX-1,0:Q_MAX-1), the field 
!    subtraction table.
!
  implicit none

  integer ( kind = 4 ), parameter :: deg = -1
  integer ( kind = 4 ), parameter :: deg_max = 50
  integer ( kind = 4 ), parameter :: q_max = 50
! integer ( kind = 4 ), parameter :: npols = 25
  integer ( kind = 4 ), parameter :: npols = 50
  integer ( kind = 4 ), parameter :: sieve_max = 400

  integer ( kind = 4 ) add(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) find
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_characteristic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ounit
  integer ( kind = 4 ) monpol(sieve_max)
  integer ( kind = 4 ) mul(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pi(-1:deg_max)
  integer ( kind = 4 ) pj(-1:deg_max)
  integer ( kind = 4 ) pk(-1:deg_max)
  integer ( kind = 4 ) ptoi
  integer ( kind = 4 ) q
  integer ( kind = 4 ) q_init
  logical              sieve(sieve_max)
  integer ( kind = 4 ) sub(0:q_max-1,0:q_max-1)

  common /field/ p, q, add, mul, sub

  save /field/

  if ( q_init <= 1 .or. q_max < q_init ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IRRED - Fatal error!'
    write ( *, '(a,i8)' ) '  Bad value of Q = ', q_init
    stop
  end if

  p = i4_characteristic ( q_init )
!
!  If no field of order Q_INIT exists, there is nothing to do.
!
  if ( p == 0 ) then
    return
  end if

  write ( *, '(a,i8)' ) '  IRRED setting up case for Q = ', q_init
!
!  Set up the field arithmetic tables.
!
  call setfld ( q_init )
!
!  Set up the sieve containing only monic polynomials.
!
  i = 0
  j = 1
  k = q

  do n = 1, sieve_max

    i = i + 1

    if ( i == j ) then
      i = k
      j = 2 * k
      k = q * k
    end if

    monpol(n) = i
    sieve(n) = .true.

  end do
!
!  Write out the irreducible polynomials as they are found.
!
  n = 0
  write ( ounit, '(i3)' ) q_init

  do i = 1, sieve_max

    if ( sieve(i) ) then

      call itop ( monpol(i), pi )
      k = pi(deg)
      write ( ounit, '(20i3)' ) k, pi(0:k)
      n = n + 1

      if ( n == npols ) then
        return
      end if

      do j = i, sieve_max

        call itop ( monpol(j), pj )
        call plymul ( pi, pj, pk )
        k = find ( ptoi ( pk ), monpol, j, sieve_max )

        if ( k /= -1 ) then
          sieve(k) = .false.
        end if

      end do

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'IRRED - Warning!'
  write ( *, '(a)' ) '  The sieve size SIEVE_MAX is too small.'
  write ( *, '(a,i8)' ) '  Number of irreducible polynomials found: ', n
  write ( *, '(a,i8)' ) '  Number needed: ', npols

  return
end
subroutine itop ( in, poly )

!*****************************************************************************80
!
!! ITOP converts an integer to a polynomial.
!
!  Discussion:
!
!    The polynomial will have coefficients in the field of order Q.
!
!    The integer is regarded as equivalent to a polynomial of the form
!      I = AN*Q**N + ... + A0.
!
!    Polynomials stored as arrays have the coefficient of degree n in 
!    POLY(N), and the degree of the polynomial in POLY(-1).  The parameter
!    DEG is just to remind us of this last fact.  A polynomial which is
!    identically 0 is given degree -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2007
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
!    Input, integer ( kind = 4 ) IN, the integer to be decoded.
!
!    Output, integer ( kind = 4 ) POLY(-1:DEG_MAX), the polynomial in the field
!    of order Q, represented by IN.  POLY(-1) contains the degree
!    of the polynomial.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) DEG_MAX, the highest degree of polynomial 
!    to be handled. 
!
!    Local, integer ( kind = 4 ) Q_MAX, the order of the largest field to be handled.
!
!  Global Parameters:
!
!    The common block /FIELD/ contains:
!
!    Common, integer ( kind = 4 ) Q, the order of the field.
!
!    Common, integer ( kind = 4 ) P, the characteristic of the field.
!
!    Common, integer ( kind = 4 ) ADD(0:Q_MAX-1,0:Q_MAX-1), the field 
!    addition table. 
!
!    Common, integer ( kind = 4 ) MUL(0:Q_MAX-1,0:Q_MAX-1), the field 
!    multiplication table. 
!
!    Common, integer ( kind = 4 ) SUB(0:Q_MAX-1,0:Q_MAX-1), the field 
!    subtraction table.
!
  implicit none

  integer ( kind = 4 ), parameter :: deg = -1
  integer ( kind = 4 ), parameter :: deg_max = 50
  integer ( kind = 4 ), parameter :: q_max = 50

  integer ( kind = 4 ) add(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) in
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mul(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) poly(-1:deg_max)
  integer ( kind = 4 ) q
  integer ( kind = 4 ) sub(0:q_max-1,0:q_max-1)

  common /field/ p, q, add, mul, sub

  save /field/

  do j = -1, deg_max
    poly(j) = 0
  end do

  i = in
  j = -1

  do

    if ( i == 0 ) then
      exit
    end if

    j = j + 1

    if ( deg_max < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ITOP - Fatal error!'
      write ( *, '(a)' ) '  Polynomial degree exceeds DEG_MAX.'
      stop
    end if

    poly(j) = mod ( i, q )
    i = i / q

  end do

  poly(deg) = j

  return
end
subroutine plymul ( pa, pb, pc )

!*****************************************************************************80
!
!! PLYMUL multiplies two polynomials over a field of order Q.
!
!  Discussion:
!
!    The parameter Q_MAX gives the order of the largest field to
!    be handled.
!
!    Polynomials stored as arrays have the coefficient of degree n in 
!    POLY(N), and the degree of the polynomial in POLY(-1).  The parameter
!    DEG is just to remind us of this last fact.  A polynomial which is
!    identically 0 is given degree -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2007
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
!    Input, integer ( kind = 4 ) PA(-1:DEG_MAX), PB(-1:DEG_MAX), the polynomials
!    to be multiplied.  Their degree is stored in entry (-1).
!
!    Output, integer ( kind = 4 ) PC(-1:DEG_MAX), the product of PA and PB,
!    in the field of order Q.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) DEG_MAX, the highest degree of polynomial
!    to be handled. 
!
!  Global Parameters:
!
!    The common block /FIELD/ contains:
!
!    Common, integer ( kind = 4 ) Q, the order of the field.
!
!    Common, integer ( kind = 4 ) P, the characteristic of the field.
!
!    Common, integer ( kind = 4 ) ADD(0:Q_MAX-1,0:Q_MAX-1), the field addition table. 
!
!    Common, integer ( kind = 4 ) MUL(0:Q_MAX-1,0:Q_MAX-1), the field 
!    multiplication table. 
!
!    Common, integer ( kind = 4 ) SUB(0:Q_MAX-1,0:Q_MAX-1), the field 
!    subtraction table.
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
    write ( *, '(a)' ) 'PLYMUL - Fatal error.'
    write ( *, '(a)' ) '  The degree of the product exceeds DEG_MAX.'
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

  do i = degc + 1, deg_max
    pc(i) = 0
  end do

  return
end
function ptoi ( poly )

!*****************************************************************************80
!
!! PTOI converts a polynomial in the field of order Q to an integer.
!
!  Discussion:
!
!    A polynomial in the field of order Q defines an integer
!      I = AN*Q**N + ... + A0.
!
!    Polynomials stored as arrays have the coefficient of degree n in 
!    POLY(N), and the degree of the polynomial in POLY(-1).  The parameter
!    DEG is just to remind us of this last fact.  A polynomial which is
!    identically 0 is given degree -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2007
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
!    Input, integer ( kind = 4 ) POLY(-1:DEG_MAX), the polynomial.
!
!    Output, integer ( kind = 4 ) PTOI, the corresponding integer.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) DEG_MAX, the highest degree of polynomial 
!    to be handled.
!
!  Global Parameters:
!
!    The common block /FIELD/ contains:
!
!    Common, integer ( kind = 4 ) Q, the order of the field.
!
!    Common, integer ( kind = 4 ) P, the characteristic of the field.
!
!    Common, integer ( kind = 4 ) ADD(0:Q_MAX-1,0:Q_MAX-1), the field 
!    addition table. 
!
!    Common, integer ( kind = 4 ) MUL(0:Q_MAX-1,0:Q_MAX-1), the field 
!    multiplication table. 
!
!    Common, integer ( kind = 4 ) SUB(0:Q_MAX-1,0:Q_MAX-1), the field 
!    subtraction table.
!
  implicit none

  integer ( kind = 4 ), parameter :: deg = -1
  integer ( kind = 4 ), parameter :: deg_max = 50
  integer ( kind = 4 ), parameter :: q_max = 50

  integer ( kind = 4 ) add(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mul(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) poly(-1:deg_max)
  integer ( kind = 4 ) ptoi
  integer ( kind = 4 ) q
  integer ( kind = 4 ) sub(0:q_max-1,0:q_max-1)

  common /field/ p, q, add, mul, sub

  save /field/

  i = 0
  do j = poly(deg), 0, -1
    i = q * i + poly(j)
  end do

  ptoi = i

  return
end
subroutine setfld ( q_init )

!*****************************************************************************80
!
!! SETFLD sets up arithmetic tables for the finite field of order QIN.
!
!  Discussion:
!
!    This subroutine sets up addition, multiplication, and
!    subtraction tables for the finite field of order QIN.
!
!    If necessary, it reads precalculated tables from the file
!    "gfarit.txt", which are supposed to have been created by GFARIT.
!
!    Polynomials stored as arrays have the coefficient of degree n in 
!    POLY(N), and the degree of the polynomial in POLY(-1).  The parameter
!    DEG is just to remind us of this last fact.  A polynomial which is
!    identically 0 is given degree -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2007
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
!    Input, Q_INIT, the order of the field.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) Q_MAX, the order of the largest field to
!    be handled.
!
!    Local, integer ( kind = 4 ) DEG_MAX, the highest degree of polynomial 
!    to be handled. 
!
!  Global Parameters:
!
!    The common block /FIELD/ contains:
!
!    Common, integer ( kind = 4 ) Q, the order of the field.
!
!    Common, integer ( kind = 4 ) P, the characteristic of the field.
!
!    Common, integer ( kind = 4 ) ADD(0:Q_MAX-1,0:Q_MAX-1), the field addition table. 
!
!    Common, integer ( kind = 4 ) MUL(0:Q_MAX-1,0:Q_MAX-1), the field 
!    multiplication table. 
!
!    Common, integer ( kind = 4 ) SUB(0:Q_MAX-1,0:Q_MAX-1), the field 
!    subtraction table.
!
  implicit none

  integer ( kind = 4 ), parameter :: deg = -1
  integer ( kind = 4 ), parameter :: deg_max = 50
  integer ( kind = 4 ), parameter :: q_max = 50

  integer ( kind = 4 ) add(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_characteristic
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mul(0:q_max-1,0:q_max-1)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) q_init
  integer ( kind = 4 ) sub(0:q_max-1,0:q_max-1)

  common /field/ p, q, add, mul, sub

  save /field/

  if ( q_init <= 1 .or. q_max < q_init ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETFLD - Fatal error!'
    write ( *, '(a,i8)' ) '  Bad value of Q = ', q_init
    stop
  end if

  q = q_init
  p = i4_characteristic ( q )

  if ( p == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETFLD - Fatal error!'
    write ( *, '(a,i8)' ) '  There is no field of order Q = ', q
    stop
  end if
!
!  Handle a field of prime order:  calculate ADD and MUL.
!
  if ( p == q ) then

    do i = 0, q - 1
      do j = 0, q - 1
        add(i,j) = mod ( i + j, p )
        mul(i,j) = mod ( i * j, p )
      end do
    end do
!
!  Handle a field of prime-power order:  tables for
!  ADD and MUL are in the file "gfarit.txt".
!
  else

    call get_unit ( iunit )

    open ( unit = iunit, file = 'gfarit.txt', status = 'old' )

    do

      read ( iunit, '(20i3)', iostat = ios ) n

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETFLD - Fatal error!'
        write ( *, '(a,i8)' ) '  Could not find tables for Q =', q
        stop
      end if

      do i = 0, n - 1
        read ( iunit, '(20i3)' ) add(i,0:n-1)
      end do

      do i = 0, n - 1
        read ( iunit, '(20i3)' ) mul(i,0:n-1)
      end do

      if ( n == q ) then
        exit
      end if

    end do

    close ( unit = iunit )

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
  integer   ( kind = 4 ) d
  character ( len = 8 )  date
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  character ( len = 10 ) time
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y
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
