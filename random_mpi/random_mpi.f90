program main

!*****************************************************************************80
!
!! MAIN is the main program for RANDOM_MPI.
!
!  Discussion:
!
!    This program demonstrates how P processors can generate the same
!    sequence of random numbers as 1 processor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
   use mpi

  integer ( kind = 4 ) a
  integer ( kind = 4 ) an
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bn
  integer ( kind = 4 ) c
  integer ( kind = 4 ) error
  integer ( kind = 4 ) id
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p
  integer ( kind = 4 ) u
  integer ( kind = 4 ) v
!
!  Initialize MPI.
!
  call MPI_Init ( error )
!
!  Get the number of processes.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, p, error )
!
!  Get the rank of this processor.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
!
!  Print a message.
!
  if ( id == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_MPI - Master process:'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a,i8)' ) '  The number of processors is P = ', p
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This program shows how a stream of random numbers'
    write ( *, '(a)' ) '  can be computed "in parallel" in an MPI program.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  We assume we are using a linear congruential'
    write ( *, '(a)' ) '  random number generator or "LCRG", which takes'
    write ( *, '(a)' ) '  an integer input and returns a new integer output:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    U = ( A * V + B ) mod C'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  We assume that we want the MPI program to produce'
    write ( *, '(a)' ) '  the same sequence of random values as a sequential'
    write ( *, '(a)' ) '  program would - but we want each processor to compute'
    write ( *, '(a)' ) '  one part of that sequence.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  We do this by computing a new LCRG which can compute'
    write ( *, '(a)' ) '  every P''th entry of the original one.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Our LCRG works with integers, but it is easy to'
    write ( *, '(a)' ) '  turn each integer into a real number between [0,1].'
  end if
!
!  A, B and C define the linear congruential random number generator.
!
  a = 16807
  b = 0
  c = 2147483647

  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  LCRG parameters:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  A  = ', a
    write ( *, '(a,i12)' ) '  B  = ', b
    write ( *, '(a,i12)' ) '  C  = ', c
  end if

  k_hi = p * 10
!
!  Processor 0 generates 10 * P random values.
!
  if ( id == 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Let processor 0 generate the entire random number sequence.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     K    ID         Input        Output'
    write ( *, '(a)' ) ' '

    k = 0
    v = 12345
    write ( *, '(2x,i4,2x,i4,2x,12x,2x,i12)' ) k, id, v

    do k = 1, k_hi
      u = v
      call lcrg_evaluate ( a, b, c, u, v )
      write ( *, '(2x,i4,2x,i4,2x,i12,2x,i12)' ) k, id, u, v
    end do

  end if
!
!  Processor P now participates by computing the P-th part of the sequence.
!
  call lcrg_anbn ( a, b, c, p, an, bn )

  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  LCRG parameters for P processors:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  AN = ', an
    write ( *, '(a,i12)' ) '  BN = ', bn
    write ( *, '(a,i12)' ) '  C  = ', c
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Have ALL the processors participate in computing'
    write ( *, '(a)' ) '  the same random number sequence.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     K    ID         Input        Output'
    write ( *, '(a)' ) ' '
  end if
!
!  Use the basis LCRG to get the ID-th value in the sequence.
!
  v = 12345
  do j = 1, id
    u = v
    call lcrg_evaluate ( a, b, c, u, v )
  end do
  k = id

  write ( *, '(2x,i4,2x,i4,2x,12x,2x,i12)' ) k, id, v
!
!  Now use the "skipping" LCRG to compute the values with indices
!  ID, ID+P, ID+2P, ...,
!
  do k = id + p, k_hi, p
    u = v
    call lcrg_evaluate ( an, bn, c, u, v )
    write ( *, '(2x,i4,2x,i4,2x,i12,2x,i12)' ) k, id, u, v
  end do
!
!  Terminate MPI.
!
  call MPI_Finalize ( error )
!
!  Terminate.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_MPI - Master process:'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

  stop
end
subroutine congruence ( a, b, c, ierror, x )

!*****************************************************************************80
!
!! CONGRUENCE solves a congruence of the form ( A * X = C ) mod B.
!
!  Discussion:
!
!    A, B and C are given integers.  The equation is solvable if and only
!    if the greatest common divisor of A and B also divides C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, C, the coefficients of the
!    Diophantine equation.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error, X was computed.
!    1, A = B = 0, C is nonzero.
!    2, A = 0, B and C nonzero, but C is not a multiple of B.
!    3, A nonzero, B zero, C nonzero, but C is not a multiple of A.
!    4, A, B, C nonzero, but GCD of A and B does not divide C.
!    5, algorithm ran out of internal space.
!
!    Output, integer ( kind = 4 ) X, the solution of the Diophantine equation.
!    X will be between 0 and B-1.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 100

  integer ( kind = 4 ) a
  integer ( kind = 4 ) a_copy
  integer ( kind = 4 ) a_mag
  integer ( kind = 4 ) a_sign
  integer ( kind = 4 ) b
  integer ( kind = 4 ) b_copy
  integer ( kind = 4 ) b_mag
  integer ( kind = 4 ) b_sign
  integer ( kind = 4 ) c
  integer ( kind = 4 ) c_copy
  integer ( kind = 4 ) g
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: one = 1
  integer ( kind = 4 ) q(nmax)
  logical swap
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
  integer ( kind = 4 ) z
!
!  Defaults for output parameters.
!
  ierror = 0
  x = 0
  y = 0
!
!  Special cases.
!
  if ( a == 0 .and. b == 0 .and. c == 0 ) then
    x = 0
    return
  else if ( a == 0 .and. b == 0 .and. c /= 0 ) then
    ierror = 1
    x = 0
    return
  else if ( a == 0 .and. b /= 0 .and. c == 0 ) then
    x = 0
    return
  else if ( a == 0 .and. b /= 0 .and. c /= 0 ) then
    x = 0
    if ( mod ( c, b ) /= 0 ) then
      ierror = 2
    end if
    return
  else if ( a /= 0 .and. b == 0 .and. c == 0 ) then
    x = 0
    return
  else if ( a /= 0 .and. b == 0 .and. c /= 0 ) then
    x = c / a
    if ( mod ( c, a ) /= 0 ) then
      ierror = 3
    end if
    return
  else if ( a /= 0 .and. b /= 0 .and. c == 0 ) then
!   g = i4_gcd ( a, b )
!   x = b / g
    x = 0
    return
  end if
!
!  Handle the "general" case: A, B and C are nonzero.
!
!  Step 1: Compute the GCD of A and B, which must also divide C.
!
  g = i4_gcd ( a, b )

  if ( mod ( c, g ) /= 0 ) then
    ierror = 4
    return
  end if

  a_copy = a / g
  b_copy = b / g
  c_copy = c / g
!
!  Step 2: Split A and B into sign and magnitude.
!
  a_mag = abs ( a_copy )
  a_sign = sign ( one, a_copy )
  b_mag = abs ( b_copy )
  b_sign = sign ( one, b_copy )
!
!  Another special case, A_MAG = 1 or B_MAG = 1.
!
  if ( a_mag == 1 ) then
    x = a_sign * c_copy
    return
  else if ( b_mag == 1 ) then
    x = 0
    return
  end if
!
!  Step 3: Produce the Euclidean remainder sequence.
!
  if ( b_mag <= a_mag ) then

    swap = .false.
    q(1) = a_mag
    q(2) = b_mag

  else

    swap = .true.
    q(1) = b_mag
    q(2) = a_mag

  end if

  n = 3

  do

    q(n) = mod ( q(n-2), q(n-1) )

    if ( q(n) == 1 ) then
      exit
    end if

    n = n + 1

    if ( nmax < n ) then
      ierror = 5
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CONGRUENCE - Fatal error!'
      write ( *, '(a)' ) '  Exceeded number of iterations.'
      stop
    end if

  end do
!
!  Step 4: Go backwards to solve X * A_MAG + Y * B_MAG = 1.
!
  y = 0
  do k = n, 2, -1
    x = y
    y = ( 1 - x * q(k-1) ) / q(k)
  end do
!
!  Step 5: Undo the swapping.
!
  if ( swap ) then
    z = x
    x = y
    y = z
  end if
!
!  Step 6: Apply signs to X and Y so that X * A + Y * B = 1.
!
  x = x * a_sign
!
!  Step 7: Multiply by C, so that X * A + Y * B = C.
!
  x = x * c_copy
!
!  Step 8: Force 0 <= X < B.
!
  x = mod ( x, b )
!
!  Step 9: Force positivity.
!
  if ( x < 0 ) then
    x = x + b
  end if

  return
end
function i4_gcd ( i, j )

!*****************************************************************************80
!
!! I4_GCD finds the greatest common divisor of I and J.
!
!  Discussion:
!
!    Only the absolute values of I and J are
!    considered, so that the result is always nonnegative.
!
!    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
!
!    If I and J have no common factor, I4_GCD is returned as 1.
!
!    Otherwise, using the Euclidean algorithm, I4_GCD is the
!    largest common factor of I and J.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer  ( kind = 4 )I, J, two numbers whose greatest common
!    divisor is desired.
!
!    Output, integer  ( kind = 4 )I4_GCD, the greatest common divisor
!    of I and J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j

  i4_gcd = 1
!
!  Return immediately if either I or J is zero.
!
  if ( i == 0 ) then
    i4_gcd = max ( 1, abs ( j ) )
    return
  else if ( j == 0 ) then
    i4_gcd = max ( 1, abs ( i ) )
    return
  end if
!
!  Set IP to the larger of I and J, IQ to the smaller.
!  This way, we can alter IP and IQ as we go.
!
  ip = max ( abs ( i ), abs ( j ) )
  iq = min ( abs ( i ), abs ( j ) )
!
!  Carry out the Euclidean algorithm.
!
  do

    ir = mod ( ip, iq )

    if ( ir == 0 ) then
      exit
    end if

    ip = iq
    iq = ir

  end do

  i4_gcd = iq

  return
end
subroutine lcrg_anbn ( a, b, c, n, an, bn )

!*****************************************************************************80
!
!! LCRG_ANBN computes the "N-th power" of a linear congruential generator.
!
!  Discussion:
!
!    We are considering a linear congruential random number generator.
!    The LCRG takes as input an integer value called SEED, and returns
!    an updated value of SEED,
!
!      SEED(out) = ( a * SEED(in) + b ) mod c.
!
!    and an associated pseudorandom real value
!
!      U = SEED(out) / c.
!
!    In most cases, a user is content to call the LCRG repeatedly, with
!    the updating of SEED being taken care of automatically.
!
!    The purpose of this routine is to determine the values of AN and BN
!    that describe the LCRG that is equivalent to N applications of the
!    original LCRG.
!
!    One use for such a facility would be to do random number computations
!    in parallel.  If each of N processors is to compute many random values,
!    you can guarantee that they work with distinct random values
!    by starting with a single value of SEED, using the original LCRG to generate
!    the first N-1 "iterates" of SEED, so that you now have N "seed" values,
!    and from now on, applying the N-th power of the LCRG to the seeds.
!
!    If the K-th processor starts from the K-th seed, it will essentially
!    be computing every N-th entry of the original random number sequence,
!    offset by K.  Thus the individual processors will be using a random
!    number stream as good as the original one, and without repeating, and
!    without having to communicate.
!
!    To evaluate the N-th value of SEED directly, we start by ignoring
!    the modular arithmetic, and working out the sequence of calculations
!    as follows:
!
!      SEED(0)   =     SEED.
!      SEED(1)   = a * SEED      + b
!      SEED(2)   = a * SEED(1)   + b = a^2 * SEED           + a * b + b
!      SEED(3)   = a * SEED(2)   + b = a^3 * SEED + a^2 * b + a * b + b
!      ...
!      SEED(N-1) = a * SEED(N-2) + b
!
!      SEED(N) = a * SEED(N-1) + b = a^N * SEED
!                                    + ( a^(n-1) + a^(n-2) + ... + a + 1 ) * b
!
!    or, using the geometric series,
!
!      SEED(N) = a^N * SEED + ( a^N - 1) / ( a - 1 ) * b
!              = AN * SEED + BN
!
!    Thus, from any SEED, we can determine the result of N applications of the
!    original LCRG directly if we can solve
!
!      ( a - 1 ) * BN = ( a^N - 1 ) * b in modular arithmetic,
!
!    and evaluate:
!
!      AN = a^N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Barry Wilkinson, Michael Allen,
!    Parallel Programming:
!    Techniques and Applications Using Networked Workstations and Parallel Computers,
!    Prentice Hall,
!    ISBN: 0-13-140563-2,
!    LC: QA76.642.W54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the multiplier for the LCRG.
!
!    Input, integer ( kind = 4 ) B, the added value for the LCRG.
!
!    Input, integer ( kind = 4 ) C, the base for the modular arithmetic.
!    For 32 bit arithmetic, this is often 2^31 - 1, or 2147483647.  It is
!    required that 0 < C.
!
!    Input, integer ( kind = 4 ) N, the "index", or number of times that the
!    LCRG is to be applied.  It is required that 0 <= N.
!
!    Output, integer ( kind = 4 ) AN, BN, the multiplier and added value for
!    the LCRG that represent N applications of the original LCRG.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) am1
  integer ( kind = 4 ) an
  integer ( kind = 4 ) anm1tb
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bn
  integer ( kind = 4 ) c
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCRG_ANBN - Fatal error!'
    write ( *, '(a,i12)' ) '  Illegal input value of N = ', n
    stop
  end if

  if ( c <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LCRG_ANBN - Fatal error!'
    write ( *, '(a,i12)' ) '  Illegal input value of C = ', c
    stop
  end if

  if ( n == 0 ) then
    an = 1
    bn = 0
  else if ( n == 1 ) then
    an = a
    bn = b
  else
!
!  Compute A^N.
!
    call power_mod ( a, n, c, an )
!
!  Solve
!    ( a - 1 ) * BN = ( a^N - 1 ) mod B
!  for BN.
!
    am1 = a - 1
    anm1tb = ( an - 1 ) * b

    call congruence ( am1, c, anm1tb, ierror, bn )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LCRG_ANBN - Fatal error!'
      write ( *, '(a)' ) '  An error occurred in the CONGRUENCE routine.'
      write ( *, '(a,i8)' ) '  The error code was IERROR = ', ierror
      stop
    end if

  end if

  return
end
subroutine lcrg_evaluate ( a, b, c, x, y )

!*****************************************************************************80
!
!! LCRG_EVALUATE evaluates an LCRG, y = ( A * x + B ) mod C.
!
!  Discussion:
!
!    This routine cannot be recommended for production use.  Because we want
!    to do modular arithmetic, but the base is not a power of 2, we need to
!    use "double precision" integers to keep accuracy.
!
!    If we knew the base C, we could try to avoid overflow while not changing
!    precision.
!
!    If the base C was a power of 2, we could rely on the usual properties of
!    integer arithmetic on computers, in which overflow bits, which are always
!    ignored, don't actually matter.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the multiplier for the LCRG.
!
!    Input, integer ( kind = 4 ) B, the added value for the LCRG.
!
!    Input, integer ( kind = 4 ) C, the base for the modular arithmetic.
!    For 32 bit arithmetic, this is often 2^31 - 1, or 2147483647.  It is
!    required that 0 < C.
!
!    Input, integer ( kind = 4 ) X, the value to be processed.
!
!    Output, integer ( kind = 4 ) Y, the processed value.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 8 ) a8
  integer ( kind = 4 ) b
  integer ( kind = 8 ) b8
  integer ( kind = 4 ) c
  integer ( kind = 8 ) c8
  integer ( kind = 4 ) x
  integer ( kind = 8 ) x8
  integer ( kind = 4 ) y
  integer ( kind = 8 ) y8
!
!  To avoid roundoff issues, we need to go to "double precision" integers.
!  (Not available on all planets.)
!
  a8 = a
  b8 = b
  c8 = c
  x8 = x

  y8 = mod ( a8 * x8 + b8, c8 )

  y = int ( y8, kind = 4 )

  if ( y < 0 ) then
    y = y + c
  end if

  return
end
subroutine power_mod ( a, n, m, x )

!*****************************************************************************80
!
!! POWER_MOD computes ( A^N ) mod M.
!
!  Discussion:
!
!    Some programming tricks are used to speed up the computation, and to
!    allow computations in which the value A**N is much too large to
!    store in an integer word.
!
!    First, for efficiency, the power A**N is computed by determining
!    the binary expansion of N, then computing A, A**2, A**4, and so on
!    by repeated squaring, and multiplying only those factors that
!    contribute to A**N.
!
!    Secondly, the intermediate products are immediately "mod'ed", which
!    keeps them small.
!
!    For instance, to compute ( A^13 ) mod 11, we essentially compute
!
!       13 = 1 + 4 + 8
!
!       A^13 = A * A^4 * A^8
!
!       A^13 ( mod 11 ) = A ( mod 11 ) * A^4 ( mod 11 ) * A^8 ( mod 11 ).
!
!    Fermat's little theorem says that if P is prime, and A is not divisible
!    by P, then ( A^(P-1) - 1 ) is divisible by P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the base of the expression to be tested.
!    0 <= A is required.
!
!    Input, integer ( kind = 4 ) N, the power to which the base is raised.
!    0 <= N is required.
!
!    Input, integer ( kind = 4 ) M, the divisor against which the expression
!    is tested.  0 < M is required.
!
!    Output, integer ( kind = 4 ) X, the remainder when A^N is divided by M.
!    If any input quantity is unacceptable, then the nonsensical value
!    X = -1 is returned.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 8 ) a_square2
  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 8 ) m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ), parameter :: two = 2
  integer ( kind = 4 ) x
  integer ( kind = 8 ) x2

  if ( a < 0 ) then
    x = -1
    return
  end if

  if ( m <= 0 ) then
    x = -1
    return
  end if

  if ( n < 0 ) then
    x = -1
    return
  end if
!
!  A_SQUARE contains the successive squares of A.
!
  a_square2 =  int ( a, kind = 8 )
  x2 = int ( 1, kind = 8 )
  m2 = int ( m, kind = 8 )

  ncopy = n

  do while ( 0 < ncopy )

    d = mod ( ncopy, two )

    if ( d == 1 ) then
      x2 = mod ( x2 * a_square2, m2 )
    end if

    a_square2 = mod ( a_square2 * a_square2, m2 )
    ncopy = ( ncopy - d ) / 2

  end do
!
!  Fix up X so that it is nonnegative.
!
  do while ( x2 < 0 )
    x2 = x2 + m2
  end do

  x = int ( x2 )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
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
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

