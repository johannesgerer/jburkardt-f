program main

!*****************************************************************************80
!
!! MAIN is the main program for GFARIT.
!
!  Discussion:
!
!    GFARIT writes the arithmetic tables called "gfarit.txt".
!
!    The program calculates addition and multiplication tables
!    for arithmetic in finite fields, and writes them out to
!    the file "gfarit.txt".  Tables are only calculated for fields
!    of prime-power order Q, the other cases being trivial.
!
!    For each value of Q, the file contains first Q, then the
!    addition table, and lastly the multiplication table.
!
!    After "gfarit.txt" has been set up, run GFPLYS to set up 
!    the file "gfplys.txt".  That operation requires reading 
!    "gfarit.txt".  
!
!    The files "gfarit.txt" and "gfplys.txt" should be saved 
!    for future use.  
!
!    Thus, a user needs to run GFARIT and GFPLYS just once,
!    before running the set of programs associated with GENIN.  
!
!    The set of programs tailored for base 2, using the driver GENIN2, 
!    requires neither "gfarit.txt" nor "gfplys.txt", hence neither GFARIT
!    nor GFPLYS.
!
!  Modified:
!
!    06 September 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, pages 494-495, 1994.
!
  implicit none

  character ( len = 80 ) :: file_out_name = 'gfarit.txt'
  integer file_out_unit
  integer, parameter :: q_max = 50
  integer qin

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GFARIT'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A program which computes a set of arithmetic'
  write ( *, '(a)' ) '  tables, and writes them to a file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tables will be created for prime and prime-power fields'
  write ( *, '(a,i8)' ) '  with orders Q = 2 through Q = ', q_max
  write ( *, '(a)' ) ' '

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace' )
  
  do qin = 2, q_max
    call gftab ( file_out_unit, qin )
  end do

  close ( unit = file_out_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GFARIT'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
subroutine gftab ( iunit, qin )

!*****************************************************************************80
!
!! GFTAB computes and writes data for a particular field size QIN.
!
!  Discussion:
!
!    A polynomial with coefficients A(*) in the field of order Q
!    can also be stored in an integer I, with
!
!      I = AN*Q**N + ... + A0.
!
!    The COMMON block, used by many subroutines,
!    gives the order Q of a field, its characteristic P, and its
!    addition, multiplication, and subtraction tables.
!    The parameter q_max gives the order of the largest field to
!    be handled.
!
!    Polynomials stored as arrays have the
!    coefficient of degree n in POLY(N), and the degree of the
!    polynomial in POLY(-1).  The parameter DEG is just to remind
!    us of this last fact.  A polynomial which is identically 0
!    is given degree -1.
!
!    The common /FIELD/ will be set up to work modulo P, the
!    characteristic of field QIN.  We construct the tables for
!    the field of order QIN in GFADD and GFMUL.
!
!    IRRPLY holds irreducible polynomials for constructing
!    prime-power fields.  IRRPLY(-2,I) says which field this
!    row is used for, and then the rest of the row is a
!    polynomial (with the degree in IRRPLY(-1,I) as usual).
!    The chosen irreducible poly is copied into MODPLY for use.
!
!  Modified:
!
!    06 September 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
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
!    Input, integer IUNIT, the unit number to which the data
!    should be written.
!
!    Input, integer Q, the order of the field for which the
!    addition and multiplication tables are needed.
!
  implicit none

  integer, parameter :: deg_max = 50
  integer, parameter :: q_max = 50

  integer add(0:q_max-1,0:q_max-1)
  integer, parameter :: deg = -1
  integer gfadd(0:q_max-1,0:q_max-1)
  integer gfmul(0:q_max-1,0:q_max-1)
  integer i
  integer i4_characteristic
  integer irrply(8,-2:deg_max)
  integer iunit
  integer j
  integer modply(-1:deg_max)
  integer mul(0:q_max-1,0:q_max-1)
  integer p
  integer pi(-1:deg_max)
  integer pj(-1:deg_max)
  integer pk(-1:deg_max)
  integer ptoi
  integer q
  integer qin
  integer sub(0:q_max-1,0:q_max-1)

  common /field/ p, q, add, mul, sub

  save /field/
  save irrply

  data irrply(1,-2:2) /4, 2, 1, 1, 1/
  data irrply(2,-2:3) /8, 3, 1, 1, 0, 1/
  data irrply(3,-2:2) /9, 2, 1, 0, 1/
  data irrply(4,-2:4) /16, 4, 1, 1, 0, 0, 1/
  data irrply(5,-2:2) /25, 2, 2, 0, 1/
  data irrply(6,-2:3) /27, 3, 1, 2, 0, 1/
  data irrply(7,-2:5) /32, 5, 1, 0, 1, 0, 0, 1/
  data irrply(8,-2:2) /49, 2, 1, 0, 1/

  if ( qin <= 1 .or. q_max < qin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GFTAB - Fatal error!'
    write ( *, '(a)' ) '  Bad value of Q.'
    return
  end if

  p = i4_characteristic ( qin )
!
!  If QIN is not a prime-power, we are not interested.
!
  if ( p == 0 .or. p == qin ) then
    return
  end if

  write ( *, '(a,i8,a,i8)' ) '  GFTAB computing table for Q = ', qin, &
    '  characteristic P = ', p
!
!  Otherwise, we set up the elements of the common /FIELD/
!  ready to do arithmetic mod P, the characteristic of QIN.
!
  call setfld ( qin )
!
!  Next find a suitable irreducible polynomial and copy it to array MODPLY.
!
  i = 1

  do while ( irrply(i,-2) /= qin )
    i = i + 1
  end do

  do j = -1, irrply(i,deg)
    modply(j) = irrply(i,j)
  end do

  do j = irrply(i,deg)+1, deg_max
    modply(j) = 0
  end do
!
!  Deal with the trivial cases ...
!
  do i = 0, qin-1
    gfadd(i,0) = i
    gfadd(0,i) = i
    gfmul(i,0) = 0
    gfmul(0,i) = 0
  end do

  do i = 1, qin-1
    gfmul(i,1) = i
    gfmul(1,i) = i
  end do
!
!  Now deal with the rest.  Each integer from 1 to QIN-1
!  is treated as a polynomial with coefficients handled mod P.
!  Multiplication of polynomials is mod MODPLY.
!
  do i = 1, qin-1

    call itop ( i, p, pi )

    do j = 1, i

      call itop ( j, p, pj )
      call plyadd ( pi, pj, pk )
      gfadd(i,j) = ptoi ( pk, p )
      gfadd(j,i) = gfadd(i,j)

      if ( 1 < i .and. 1 < j ) then
        call plymul ( pi, pj, pk )
        call plydiv ( pk, modply, pj, pk )
        gfmul(i,j) = ptoi ( pk, p )
        gfmul(j,i) = gfmul(i,j)
      end if

    end do

  end do
!
!  Write out the tables.
!
  write ( iunit, '(20i3)' ) qin

  do i = 0, qin-1
    write ( iunit, '(20i3)' ) gfadd(i,0:qin-1)
  end do

  do i = 0, qin-1
    write ( iunit, '(20i3)' ) gfmul(i,0:qin-1)
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
!  Modified:
!
!    05 September 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
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
subroutine itop ( in, p, poly )

!*****************************************************************************80
!
!! ITOP converts an integer to a polynomial in the field of order P.
!
!  Discussion:
!
!    A nonnegative integer IN can be decomposed into a polynomial in
!    powers of P, with coefficients between 0 and P-1, by setting:
!
!      J = 0
!      do while ( 0 < IN )
!        POLY(J) = mod ( IN, P )
!        J = J + 1
!        IN = IN / P
!      end do
!
!  Modified:
!
!    05 September 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
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
!    Input, integer IN, the (nonnegative) integer containing the 
!    polynomial information.
!
!    Input, integer P, the order of the field.
!
!    Output, integer POLY(-1:deg_max), the polynomial information.
!    POLY(-1) contains the degree of the polynomial.  POLY(I) contains
!    the coefficient of degree I.  Each coefficient is an element of
!    the field of order P; in other words, each coefficient is
!    between 0 and P-1.
!
  implicit none

  integer, parameter :: deg = -1
  integer, parameter :: deg_max = 50

  integer i
  integer in
  integer j
  integer p
  integer poly(-1:deg_max)

  do j = -1, deg_max
    poly(j) = 0
  end do

  i = in
  j = -1

  do while ( 0 < i )

    j = j + 1

    if ( deg_max < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ITOP - Fatal error!'
      write ( *, '(a)' ) '  The polynomial degree exceeds deg_max.'
      stop
    end if

    poly(j) = mod ( i, p )
    i = i / p

  end do

  poly(deg) = j

  return
end
subroutine plyadd ( pa, pb, pc )

!*****************************************************************************80
!
!! PLYADD adds two polynomials.
!
!  Discussion:
!
!    POLY(-1) contains the degree of the polynomial.  POLY(I) contains
!    the coefficient of degree I.  Each coefficient is an element of
!    the field of order Q; in other words, each coefficient is
!    between 0 and Q-1.
!
!  Modified:
!
!    05 September 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
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
!    Input, integer PA(-1:deg_max), the first polynomial.
!
!    Input, integer PB(-1:deg_max), the second polynomial.
!
!    Output, integer PC(-1:deg_max), the sum polynomial.
!
  implicit none

  integer, parameter :: deg = -1
  integer, parameter :: deg_max = 50
  integer, parameter :: q_max = 50

  integer add(0:q_max-1,0:q_max-1)
  integer degc
  integer i
  integer maxab
  integer mul(0:q_max-1,0:q_max-1)
  integer p
  integer pa(-1:deg_max)
  integer pb(-1:deg_max)
  integer pc(-1:deg_max)
  integer q
  integer sub(0:q_max-1,0:q_max-1)

  common /field/ p, q, add, mul, sub

  save /field/

  maxab = max ( pa(deg), pb(deg) )

  degc = -1

  do i = 0, maxab
    pc(i) = add ( pa(i), pb(i) )
    if ( pc(i) /= 0 ) then
      degc = i
    end if
  end do

  pc(deg) = degc

  do i = maxab+1, deg_max
    pc(i) = 0
  end do

  return
end
subroutine plydiv ( pa, pb, pq, pr )

!*****************************************************************************80
!
!! PLYDIV divides one polynomial by another.
!
!  Discussion:
!
!    Polynomial coefficients are elements of the field of order Q.
!
!    The COMMON block, used by many subroutines,
!    gives the order Q of a field, its characteristic P, and its
!    addition, multiplication, and subtraction tables.
!    The parameter q_max gives the order of the largest field to
!    be handled.
!
!  Modified:
!
!    05 September 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
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
!    Input, integer PA(-1:deg_max), the first polynomial.
!
!    Input, integer PB(-1:deg_max), the second polynomial.
!
!    Output, integer PQ(-1:deg_max), the quotient polynomial.
!
!    Output, integer PR(-1:deg_max), the remainder polynomial.
!
  implicit none

  integer, parameter :: deg = -1
  integer, parameter :: deg_max = 50
  integer, parameter :: q_max = 50

  integer add(0:q_max-1,0:q_max-1)
  integer binv
  integer d
  integer degb
  integer degq
  integer degr
  integer i
  integer j
  integer m
  integer mul(0:q_max-1,0:q_max-1)
  integer p
  integer pa(-1:deg_max)
  integer pb(-1:deg_max)
  integer pq(-1:deg_max)
  integer pr(-1:deg_max)
  integer ptq(-1:deg_max)
  integer ptr(-1:deg_max)
  integer q
  integer sub(0:q_max-1,0:q_max-1)

  common /field/ p, q, add, mul, sub

  save /field/

  if ( pb(deg) == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLYDIV -  Fatal error!'
    write ( *, '(a)' ) '  Division by zero polynomial.'
    stop
  end if

  do i = -1, deg_max
    ptq(i) = 0
    ptr(i) = pa(i)
  end do

  degr = pa(deg)
  degb = pb(deg)
  degq = degr - degb

  if ( degq < 0 ) then
    degq = -1
  end if
!
!  Find the inverse of the leading coefficient of PB.
!
  j = pb(degb)
  do i = 1, p-1
    if ( mul(i,j) == 1 ) then
      binv = i
    end if
  end do

  do d = degq, 0, -1
    m = mul ( ptr(degr), binv )
    do i = degb, 0, -1
      ptr(degr+i-degb) = sub ( ptr(degr+i-degb), mul ( m, pb(i) ) )
    end do
    degr = degr - 1
    ptq(d) = m
  end do

  pq(0:deg_max) = ptq(0:deg_max)
  pr(0:deg_max) = ptr(0:deg_max)

  pq(deg) = degq

  do while ( pr(degr) == 0 .and. 0 <= degr ) 
    degr = degr - 1
  end do

  pr(deg) = degr

  return
end
subroutine plymul ( pa, pb, pc )

!*****************************************************************************80
!
!! PLYMUL multiplies one polynomial by another.
!
!  Discussion:
!
!    Polynomial coefficients are elements of the field of order Q.
!
!    The COMMON block, used by many subroutines,
!    gives the order Q of a field, its characteristic P, and its
!    addition, multiplication, and subtraction tables.
!    The parameter q_max gives the order of the largest field to
!    be handled.
!
!  Modified:
!
!    05 September 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
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
!    Input, integer PA(-1:deg_max), the first polynomial.
!
!    Input, integer PB(-1:deg_max), the second polynomial.
!
!    Output, integer PC(-1:deg_max), the product polynomial.
!
  implicit none

  integer, parameter :: deg = -1
  integer, parameter :: deg_max = 50
  integer, parameter :: q_max = 50

  integer add(0:q_max-1,0:q_max-1)
  integer dega
  integer degb
  integer degc
  integer i
  integer j
  integer mul(0:q_max-1,0:q_max-1)
  integer p
  integer pa(-1:deg_max)
  integer pb(-1:deg_max)
  integer pc(-1:deg_max)
  integer pt(-1:deg_max)
  integer q
  integer sub(0:q_max-1,0:q_max-1)
  integer term

  common /field/ p, q, add, mul, sub

  save /field/

  dega = pa(deg)
  degb = pb(deg)

  if ( dega == -1 .or. degb == -1) then
    degc = -1
  else
    degc = dega + degb
  end if

  if ( deg_max < degc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLYMUL - Fatal error!'
    write ( *, '(a)' ) '  The degree of the product exceeds deg_max.'
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

  do i = degc+1, deg_max
    pc(i) = 0
  end do

  return
end
function ptoi ( poly, q )

!*****************************************************************************80
!
!! PTOI converts a polynomial in the field of order Q to an integer.
!
!  Discussion:
!
!    A polynomial with coefficients A(*) in the field of order Q
!    can also be stored in an integer I, with
!
!      I = AN*Q**N + ... + A0.
!
!  Modified:
!
!    05 September 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
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
!    Input, integer POLY(-1:deg_max), the polynomial information.
!    POLY(-1) contains the degree of the polynomial.  POLY(I) contains
!    the coefficient of degree I.  Each coefficient is an element of
!    the field of order P; in other words, each coefficient is
!    between 0 and P-1.
!
!    Input, integer P, the order of the field.
!
!    Output, integer PTOI, the (nonnegative) integer containing the 
!    polynomial information.
!
  implicit none

  integer, parameter :: deg_max = 50

  integer, parameter :: deg = -1
  integer i
  integer j
  integer poly(-1:deg_max)
  integer ptoi
  integer q

  i = 0

  do j = poly(deg), 0, -1
    i = q * i + poly(j)
  end do

  ptoi = i

  return
end
subroutine setfld ( qin )

!*****************************************************************************80
!
!! SETFLD sets up the arithmetic tables for a finite field.
!
!  Discussion:
!
!    This subroutine sets up addition, multiplication, and
!    subtraction tables for the finite field of order QIN.
!
!    If necessary, it reads precalculated tables from the file
!    "gfarit.txt".  These precalculated tables
!    are supposed to have been put there by GFARIT.
!
!    A polynomial with coefficients A(*) in the field of order Q
!    can also be stored in an integer I, with
!
!      I = AN*Q**N + ... + A0.
!
!    The COMMON block, used by many subroutines,
!    gives the order Q of a field, its characteristic P, and its
!    addition, multiplication, and subtraction tables.
!    The parameter q_max gives the order of the largest field to
!    be handled.
!
!  Modified:
!
!    05 September 2007
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter.
!    FORTRAN90 version by John Burkardt.
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
  implicit none

  integer, parameter :: deg = -1
  integer, parameter :: deg_max = 50
  integer, parameter :: q_max = 50

  integer add(0:q_max-1,0:q_max-1)
  integer i
  integer i4_characteristic
  integer ios
  integer iunit
  integer j
  integer mul(0:q_max-1,0:q_max-1)
  integer n
  integer p
  integer q
  integer qin
  integer sub(0:q_max-1,0:q_max-1)

  common /field/ p, q, add, mul, sub

  save /field/

  if ( qin <= 1 .or. q_max < qin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETFLD - Fatal error!'
    write ( *, '(a,i8)' ) '  Bad value of Q = ', q
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
!  Set up to handle a field of prime or prime-power order.
!  Calculate the addition and multiplication tables.
!
  do i = 0, p-1
    do j = 0, p-1
      add(i,j) = mod ( i + j, p )
    end do
  end do

  do i = 0, p-1
    do j = 0, p-1
      mul(i,j) = mod ( i * j, p )
    end do
  end do

  do i = 0, p-1
    do j = 0, p-1
      sub ( add(i,j), i ) = j
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
