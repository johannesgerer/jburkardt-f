subroutine backtrack ( l, iarray, indx, k, nstack, stack, maxstack )

!*****************************************************************************80
!
!! BACKTRACK supervises a backtrack search.
!
!  Discussion:
!
!    The routine builds a vector, one element at a time, which is
!    required to satisfy some condition.
!
!    At any time, the partial vector may be discovered to be
!    unsatisfactory, but the routine records information about where the
!    last arbitrary choice was made, so that the search can be
!    carried out efficiently, rather than starting out all over again.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 1999
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) L, the length of the completed
!    candidate vector.
!
!    Input/output, integer ( kind = 4 ) IARRAY(L), the candidate vector.
!
!    Input/output, integer ( kind = 4 ) INDX.
!    On input, set INDX = 0 to start a search.
!    On output:
!    1, a complete output vector has been determined.
!    2, candidates are needed.
!    3, no more possible vectors exist.
!
!    Input/output, integer ( kind = 4 ) K, the current length of the candidate
!    vector.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!
!    Input/output, integer ( kind = 4 ) STACK(MAXSTACK), a list of candidates
!    for positions 1 through K.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum length of the stack.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) iarray(l)
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) stack(maxstack)
!
!  If this is the first call, request a candidate for position 1.
!
  if ( indx == 0 ) then
    k = 1
    nstack = 0
    indx = 2
    return
  end if
!
!  Examine the stack.
!
  do

    nstack = nstack - 1
!
!  If there are candidates for position K, take the first available
!  one off the stack, and increment K.
!
!  This may cause K to reach the desired value of L, in which case
!  we need to signal the user that a complete set of candidates
!  is being returned.
!
    if ( stack(nstack+1) /= 0 ) then

      iarray(k) = stack(nstack)
      stack(nstack) = stack(nstack+1) - 1

      if ( k /= l ) then
        k = k + 1
        indx = 2
      else
        indx = 1
      end if

      exit
!
!  If there are no candidates for position K, then decrement K.
!  If K is still positive, repeat the examination of the stack.
!
    else

      k = k - 1

      if ( k <= 0 ) then
        indx = 3
        exit
      end if

    end if

  end do

  return
end
subroutine bal_seq_check ( n, t, ierror )

!*****************************************************************************80
!
!! BAL_SEQ_CHECK checks a balanced sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of 0's (and 1's) in the sequence.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(2*N), a balanced sequence.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    -1, N is not positive.
!    I, the I-th entry of T is illegal.
!    2*N+1, there are not the same number of 1's as 0's.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) one_count
  integer ( kind = 4 ) t(2*n)
  integer ( kind = 4 ) zero_count

  ierror = 0

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BAL_SEQ_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  one_count = 0
  zero_count = 0

  do i = 1, 2 * n

    if ( t(i) == 0 ) then
      zero_count = zero_count + 1
    else if ( t(i) == 1 ) then
      one_count = one_count + 1
    else
      ierror = i
      return
    end if

    if ( zero_count < one_count ) then
      ierror = 1
      return
    end if

  end do

  if ( one_count /= zero_count ) then
    ierror = 2 * n + 1
  end if

  return
end
subroutine bal_seq_enum ( n, nseq )

!*****************************************************************************80
!
!! BAL_SEQ_ENUM enumerates the balanced sequences.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of 0's (and 1's) in the sequence.
!    N must be nonnegative.
!
!    Output, integer ( kind = 4 ) NSEQ, the number of balanced sequences.
!
  implicit none

  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nseq

  nseq = binomial ( 2 * n, n ) / ( n + 1 )

  return
end
subroutine bal_seq_rank ( n, t, rank )

!*****************************************************************************80
!
!! BAL_SEQ_RANK ranks a balanced sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of 0's (and 1's) in the sequence.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(2*N), a balanced sequence.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the balanced sequence.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) mxy
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(2*n)
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
!
!  Check.
!
  call bal_seq_check ( n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BAL_SEQ_RANK - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  y = 0
  rank = 0

  do x = 1, 2 * n - 1

    if ( t(x) == 0 ) then
      y = y + 1
    else
      call mountain ( n, x, y + 1, mxy )
      rank = rank + mxy
      y = y - 1
    end if

  end do

  return
end
subroutine bal_seq_successor ( n, t, rank )

!*****************************************************************************80
!
!! BAL_SEQ_SUCCESSOR computes the lexical balanced sequence successor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of 0's (and 1's) in the sequence.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) T(2*N), on input, a balanced sequence,
!    and on output, its lexical successor.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) open
  integer ( kind = 4 ) open_index
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) slot
  integer ( kind = 4 ) slot_index
  integer ( kind = 4 ) slot_ones
  integer ( kind = 4 ) t(2*n)
!
!  Return the first element.
!
  if ( rank == -1 ) then
    t(1:n) = 0
    t(n+1:2*n) = 1
    rank = 0
    return
  end if
!
!  Check.
!
  call bal_seq_check ( n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BAL_SEQ_SUCCESSOR - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  After the I-th 0 there is a 'slot' with the capacity to
!  hold between 0 and I ones.
!
!  The first element of the sequence has all the 1's cowering
!  behind the N-th 0.
!
!  We seek to move a 1 to the left, and to do it lexically,
!  we will move a 1 to the rightmost slot that is under capacity.
!
!  Find the slot.
!
  slot = 0
  slot_index = 0
  slot_ones = 0

  open = 0
  open_index = 0

  do i = 1, 2 * n

    if ( t(i) == 0 ) then

      if ( 0 < slot ) then
        if ( slot_ones < slot ) then
          open = slot
          open_index = slot_index
        end if
      end if

      slot = slot + 1
      slot_index = i

!     slot_ones = 0

    else
      slot_ones = slot_ones + 1
    end if

  end do
!
!  If OPEN is not 0, then preserve the string up to the OPEN-th 0,
!  preserve the 1's that follow, but then write a 1, then
!  all the remaining 0's and all the remaining 1's.
!
  if ( open /= 0 ) then

    j = open_index + 1

    do while ( t(j) == 1 )
      j = j + 1
    end do

    t(j) = 1

    do i = open + 1, n
      j = j + 1
      t(j) = 0
    end do

    t(j+1:2*n) = 1
!
!  If OPEN is 0, the last element was input.
!  Return the first one.
!
  else

    t(1:n) = 0
    t(n+1:2*n) = 1
    rank = 0
    return

  end if

  rank = rank + 1

  return
end
subroutine bal_seq_to_tableau ( n, t, tab )

!*****************************************************************************80
!
!! BAL_SEQ_TO_TABLEAU converts a balanced sequence to a 2 by N tableau.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of 0's (and 1's) in the sequence.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(2*N), a balanced sequence.
!
!    Output, integer ( kind = 4 ) TAB(2,N), a 2 by N tableau.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) r
  integer ( kind = 4 ) t(2*n)
  integer ( kind = 4 ) tab(2,n)
!
!  Check.
!
  call bal_seq_check ( n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BAL_SEQ_TO_TABLEAU - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  c(1) = 0
  c(2) = 0
  do i = 1, 2 * n
    r = t(i) + 1
    c(r) = c(r) + 1
    tab(r,c(r)) = i
  end do

  return
end
subroutine bal_seq_unrank ( rank, n, t )

!*****************************************************************************80
!
!! BAL_SEQ_UNRANK unranks a balanced sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the balanced sequence.
!
!    Input, integer ( kind = 4 ) N, the number of 0's (and 1's) in the sequence.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) T(2*N), a balanced sequence.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) low
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nseq
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(2*n)
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BAL_SEQ_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call bal_seq_enum ( n, nseq )

  if ( rank < 0 .or. nseq < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BAL_SEQ_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  y = 0
  low = 0

  do x = 1, 2 * n

    call mountain ( n, x, y + 1, m )

    if ( rank <= low + m - 1 ) then
      y = y + 1
      t(x) = 0
    else
      low = low + m
      y = y - 1
      t(x) = 1
    end if

  end do

  return
end
subroutine bell_numbers ( m, b )

!*****************************************************************************80
!
!! BELL_NUMBERS computes the Bell numbers.
!
!  Discussion:
!
!    There are B(M) restricted growth functions of length M.
!
!    There are B(M) partitions of a set of M objects.
!
!    B(M) is the sum of the Stirling numbers of the second kind,
!    S(M,N), for N = 0 to M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, indicates how many Bell numbers are to
!    compute.  M must be nonnegative.
!
!    Output, integer ( kind = 4 ) B(0:M), the first M+1 Bell numbers.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) b(0:m)
  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  b(0) = 1
  do j = 1, m
    b(j) = 0
    do i = 0, j - 1
      b(j) = b(j) + binomial ( j - 1, i ) * b(i)
    end do
  end do

  return
end
subroutine bell_values ( n_data, n, c )

!*****************************************************************************80
!
!! BELL_VALUES returns some values of the Bell numbers.
!
!  Discussion:
!
!    The Bell number B(N) is the number of restricted growth functions on N.
!
!    Note that the Stirling numbers of the second kind, S^m_n, count the
!    number of partitions of N objects into M classes, and so it is
!    true that
!
!      B(N) = S^1_N + S^2_N + ... + S^N_N.
!
!    The Bell numbers were named for Eric Temple Bell.
!
!    In Mathematica, the function can be evaluated by
!
!      Sum[StirlingS2[n,m],{m,1,n}]
!
!  Definition:
!
!    The Bell number B(N) is defined as the number of partitions (of
!    any size) of a set of N distinguishable objects.
!
!    A partition of a set is a division of the objects of the set into
!    subsets.
!
!  Examples:
!
!    There are 15 partitions of a set of 4 objects:
!
!      (1234),
!      (123) (4),
!      (124) (3),
!      (12) (34),
!      (12) (3) (4),
!      (134) (2),
!      (13) (24),
!      (13) (2) (4),
!      (14) (23),
!      (1) (234),
!      (1) (23) (4),
!      (14) (2) (3),
!      (1) (24) (3),
!      (1) (2) (34),
!      (1) (2) (3) (4).
!
!    and so B(4) = 15.
!
!  First values:
!
!     N         B(N)
!     0           1
!     1           1
!     2           2
!     3           5
!     4          15
!     5          52
!     6         203
!     7         877
!     8        4140
!     9       21147
!    10      115975
!
!  Recursion:
!
!    B(I) = sum ( 1 <= J <=I ) Binomial ( I-1, J-1 ) * B(I-J)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the Bell number.
!
!    Output, integer ( kind = 4 ) C, the value of the Bell number.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( n_max ) :: c_vec = (/ &
    1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
function binomial ( n, k )

!*****************************************************************************80
!
!! BINOMIAL computes the binomial coefficient C(N,K).
!
!  Discussion:
!
!    BINOMIAL(N,K) = C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, integer ( kind = 4 ) BINOMIAL, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icnk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    icnk = 0

  else if ( mn == 0 ) then

    icnk = 1

  else

    mx = max ( k, n - k )
    icnk = mx + 1

    do i = 2, mn
      icnk = ( icnk * ( mx + i ) ) / i
    end do

  end if

  binomial = icnk

  return
end
subroutine combin ( n, k, cnk )

!*****************************************************************************80
!
!! COMBIN computes the combinatorial coefficient C(N,K).
!
!  Discussion:
!
!    Real arithmetic is used, and C(N,K) is computed directly, via
!    Gamma functions, rather than recursively.
!
!    C(N,K) is the number of distinct combinations of K objects
!    chosen from a set of N distinct objects.  A combination is
!    like a set, in that order does not matter.
!
!    C(N,K) = N! / ( (N-K)! * K! )
!
!  Example:
!
!    The number of combinations of 2 things chosen from 5 is 10.
!
!    C(5,2) = ( 5 * 4 * 3 * 2 * 1 ) / ( ( 3 * 2 * 1 ) * ( 2 * 1 ) ) = 10.
!
!    The actual combinations may be represented as:
!
!      (1,2), (1,3), (1,4), (1,5), (2,3),
!      (2,4), (2,5), (3,4), (3,5), (4,5).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the value of N.
!
!    Input, integer ( kind = 4 ) K, the value of K.
!
!    Output, real ( kind = 8 ) CNK, the value of C(N,K)
!
  implicit none

  real ( kind = 8 ) arg
  real ( kind = 8 ) cnk
  real ( kind = 8 ) fack
  real ( kind = 8 ) facn
  real ( kind = 8 ) facnmk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_gamma_log

  if ( n < 0 ) then

    cnk = 0.0D+00

  else if ( k == 0 ) then

    cnk = 1.0D+00

  else if ( k == 1 ) then

    cnk = real ( n, kind = 8 )

  else if ( 1 < k .and. k < n - 1 ) then

    arg = real ( n + 1, kind = 8 )
    facn = r8_gamma_log ( arg )

    arg = real ( k + 1, kind = 8 )
    fack = r8_gamma_log ( arg )

    arg = real ( n - k + 1, kind = 8 )
    facnmk = r8_gamma_log ( arg )

    cnk = real ( nint ( exp ( facn - fack - facnmk ) ),  kind = 8 )

  else if ( k == n - 1 ) then

    cnk = real ( n, kind = 8 )

  else if ( k == n ) then

    cnk = 1.0D+00

  else

    cnk = 0.0D+00

  end if

  return
end
subroutine cycle_check ( n, ncycle, t, index, ierror )

!*****************************************************************************80
!
!! CYCLE_CHECK checks a permutation in cycle form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items permuted.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NCYCLE, the number of cycles.
!    1 <= NCYCLE <= N.
!
!    Input, integer ( kind = 4 ) T(N), INDEX(NCYCLE), describes the permutation
!    as a collection of NCYCLE cycles.  The first cycle is
!    T(1) -> T(2) -> ... -> T(INDEX(1)) -> T(1).
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    -1, N is less than 1.
!    -2, NCYCLE is less than 1 or greater than N.
!    -3, an entry of INDEX is illegal.
!    -4, the entries of INDEX do not add up to N.
!    I, entry I of T is illegal.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncycle

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifind
  integer ( kind = 4 ) index(ncycle)
  integer ( kind = 4 ) iseek
  integer ( kind = 4 ) t(n)

  ierror = 0
!
!  N must be at least 1.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CYCLE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if
!
!  1 <= NCYCLE <= N.
!
  if ( ncycle < 1 .or. n < ncycle ) then
    ierror = -2
    return
  end if
!
!  1 <= INDEX(I) <= N.
!
  do i = 1, ncycle
    if ( index(i) < 1 .or. n < index(i) ) then
      ierror = -3
      return
    end if
  end do
!
!  The INDEX(I)'s sum to N.
!
  if ( sum ( index(1:ncycle) ) /= n ) then
    ierror = -4
    return
  end if
!
!  1 <= T(I) <= N.
!
  do i = 1, n
    if ( t(i) < 1 .or. n < t(i) ) then
      ierror = i
      return
    end if
  end do
!
!  Verify that every value from 1 to N occurs in T.
!
  do iseek = 1, n

    ifind = 0

    do i = 1, n
      if ( t(i) == iseek ) then
        ifind = i
        exit
      end if
    end do

    if ( ifind == 0 ) then
      ierror = iseek
      return
    end if

  end do

  return
end
subroutine cycle_to_perm ( n, ncycle, t, index, p )

!*****************************************************************************80
!
!! CYCLE_TO_PERM converts a permutation from cycle to array form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items permuted.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NCYCLE, the number of cycles.
!    1 <= NCYCLE <= N.
!
!    Input, integer ( kind = 4 ) T(N), INDEX(NCYCLE), describes the permutation
!    as a collection of NCYCLE cycles.  The first cycle is
!    T(1) -> T(2) -> ... -> T(INDEX(1)) -> T(1).
!
!    Output, integer ( kind = 4 ) P(N), describes the permutation using a
!    single array.  For each index I, I -> P(I).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncycle

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) index(ncycle)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) t(n)
!
!  Check.
!
  call cycle_check ( n, ncycle, t, index, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CYCLE_TO_PERM - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  jhi = 0

  do i = 1, ncycle

    jlo = jhi + 1
    jhi = jhi + index(i)

    do j = jlo, jhi

      if ( j < jhi ) then
        p(t(j)) = t(j+1)
      else
        p(t(j)) = t(jlo)
      end if

    end do

  end do

  return
end
subroutine dist_enum ( k, m, num_dist )

!*****************************************************************************80
!
!! DIST_ENUM returns the number of distributions of indistinguishable objects.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the number of distinguishable "slots".
!
!    Input, integer ( kind = 4 ) M, the number of indistinguishable objects.
!
!    Output, integer ( kind = 4 ) NUM_DIST, the number of distributions of M
!    indistinguishable objects about K distinguishable slots.
!
  implicit none

  real ( kind = 8 ) cnk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) num_dist

  call combin ( m + k - 1, m, cnk )

  num_dist = nint ( cnk )

  return
end
subroutine dist_next ( k, m, q, more )

!*****************************************************************************80
!
!! DIST_NEXT returns the next distribution of indistinguishable objects.
!
!  Discussion:
!
!    A distribution of M objects into K parts is an ordered sequence
!    of K nonnegative integers which sum to M.  This is similar to
!    a partition of a set into K subsets, except that here the order
!    matters.  That is, (1,1,2) and (1,2,1) are considered to be
!    different distributions.
!
!    On the first call to this routine, the user should set MORE = FALSE,
!    to signal that this is a startup for the given computation.  The routine
!    will return the first distribution, and set MORE = TRUE.
!
!    If the user calls again, with MORE = TRUE, the next distribution
!    is being requested.  If the routine returns with MORE = TRUE, then
!    that distribution was found and returned.  However, if the routine
!    returns with MORE = FALSE, then no more distributions were found;
!    the enumeration of distributions has terminated.
!
!    A "distribution of M indistinguishable objects into K slots" is
!    sometimes called a "composition of the integer M into K parts".
!
!  Example:
!
!    K = 3, M = 5
!
!    0           0           5
!    0           1           4
!    0           2           3
!    0           3           2
!    0           4           1
!    0           5           0
!    1           0           4
!    1           1           3
!    1           2           2
!    1           3           1
!    1           4           0
!    2           0           3
!    2           1           2
!    2           2           1
!    2           3           0
!    3           0           2
!    3           1           1
!    3           2           0
!    4           0           1
!    4           1           0
!    5           0           0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Fenichel,
!    Algorithm 329:
!    Distribution of Indistinguishable Objects into
!    Distinguishable Slots,
!    Communications of the ACM,
!    Volume 11, Number 6, June 1968, page 430.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the number of distinguishable "slots".
!
!    Input, integer ( kind = 4 ) M, the number of indistinguishable objects.
!
!    Input/output, integer ( kind = 4 ) Q(K), the number of objects in each
!    slot.
!
!    Input/output, logical MORE, used by the user to start the computation,
!    and by the routine to stop the computation.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ), save :: leftmost = 1
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) q(k)
!
!  The startup call.
!
  if ( .not. more ) then

    more = .true.

    q(1:k-1) = 0
    q(k) = m

    leftmost = k + 1
!
!  There are no more distributions.
!  Reset Q to the first distribution in the sequence.
!
  else if ( q(1) == m ) then

    more = .false.

    q(1:k-1) = 0
    q(k) = m

    leftmost = k + 1

  else if ( leftmost < k + 1 ) then

    leftmost = leftmost - 1
    q(k) = q(leftmost) - 1
    q(leftmost) = 0
    q(leftmost-1) = q(leftmost-1) + 1
    if ( q(k) /= 0 ) then
      leftmost = k + 1
    end if

  else

    if ( q(k) == 1 ) then
      leftmost = k
    end if

    q(k) = q(k) - 1
    q(k-1) = q(k-1) + 1

  end if

  return
end
subroutine edge_check ( n_node, n_edge, t, ierror )

!*****************************************************************************80
!
!! EDGE_CHECK checks a graph stored by edges.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_NODE, the number of nodes in the graph.
!    N_NODE must be positive.
!
!    Input, integer ( kind = 4 ) N_EDGE, the number of edges in the graph.
!    N_EDGE must be positive.
!
!    Input, integer ( kind = 4 ) T(2,N_EDGE), describes the edges of the tree
!    as pairs of nodes.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    -1, N_NODE is not positive.
!    -2, N_EDGE is not positive.
!    0, no error.
!    I, edge T(1,I), T(2,I) is illegal.
!
  implicit none

  integer ( kind = 4 ) n_edge
  integer ( kind = 4 ) n_node

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) t(2,n_edge)

  ierror = 0

  if ( n_node < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EDGE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N_NODE < 1.'
    stop
  end if

  if ( n_edge < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EDGE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N_EDGE < 1.'
    stop
  end if
!
!  Every edge must join two legal nodes.
!
  do i = 1, 2
    do j = 1, n_edge
      if ( t(i,j) < 1 .or. n_node < t(i,j) ) then
        ierror = i
        return
      end if
    end do
  end do
!
!  Every edge must join distinct nodes.
!
  do j = 1, n_edge
    if ( t(1,j) == t(2,j) ) then
      ierror = i
      return
    end if
  end do
!
!  Every edge must be distinct.
!
  do j = 1, n_edge - 1
    do j2 = j + 1, n_edge
      if ( t(1,j) == t(1,j2) .and. t(2,j) == t(2,j2) ) then
        ierror = j2
        return
      else if ( t(1,j) == t(2,j2) .and. t(2,j) == t(1,j2) ) then
        ierror = j2
        return
      end if
    end do
  end do

  return
end
subroutine edge_degree ( n_node, n_edge, t, d )

!*****************************************************************************80
!
!! EDGE_DEGREE returns the degree of the nodes of a graph stored by edges.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_NODE, the number of nodes in the graph.
!    N_NODE must be positive.
!
!    Input, integer ( kind = 4 ) N_EDGE, the number of edges in the graph.
!    N_EDGE must be positive.
!
!    Input, integer ( kind = 4 ) T(2,N_EDGE), describes the edges of the tree
!    as pairs of nodes.
!
!    Output, integer ( kind = 4 ) D(N_NODE), the degree of each node.
!
  implicit none

  integer ( kind = 4 ) n_edge
  integer ( kind = 4 ) n_node

  integer ( kind = 4 ) d(n_node)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) t(2,n_edge)
!
!  Check.
!
  call edge_check ( n_node, n_edge, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EDGE_DEGREE - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  Compute the degree of each node.
!
  d(1:n_node) = 0

  do j = 1, n_edge
    d(t(1,j)) = d(t(1,j)) + 1
    d(t(2,j)) = d(t(2,j)) + 1
  end do

  return
end
subroutine edge_enum ( n_node, nedge )

!*****************************************************************************80
!
!! EDGE_ENUM enumerates the maximum number of edges in a graph on N_NODE nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_NODE, the number of nodes in the graph.
!    N_NODE must be positive.
!
!    Output, integer ( kind = 4 ) NEDGE, the maximum number of edges in a graph
!    on N_NODE nodes.
!
  implicit none

  integer ( kind = 4 ) n_node
  integer ( kind = 4 ) nedge

  nedge = ( n_node * ( n_node - 1 ) ) / 2

  return
end
function fall ( x, n )

!*****************************************************************************80
!
!! FALL computes the falling factorial function [X]_N.
!
!  Discussion:
!
!    The number of "injections" or 1-to-1 mappings from
!    a set of N elements to a set of M elements is [M]_N.
!
!    The number of permutations of N objects out of M is [M}_N.
!
!    The Stirling numbers of the first kind can be used
!    to convert a falling factorial into a polynomial, as follows:
!
!      [X]_N = S^0_N + S^1_N * X + S^2_N * X^2 + ... + S^N_N X^N.
!
!  Formula:
!
!    [X]_N = X * ( X - 1 ) * ( X - 2 ) * ... * ( X - N + 1 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the falling factorial
!    function.
!
!    Input, integer ( kind = 4 ) N, the order of the falling factorial function.
!    If N = 0, FALL = 1, if N = 1, FALL = X.  Note that if N is
!    negative, a "rising" factorial will be computed.
!
!    Output, integer ( kind = 4 ) FALL, the falling factorial function.
!
  implicit none

  integer ( kind = 4 ) arg
  integer ( kind = 4 ) fall
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x

  fall = 1

  arg = x

  if ( 0 < n ) then

    do i = 1, n
      fall = fall * arg
      arg = arg - 1
    end do

  else if ( n < 0 ) then

    do i = n, -1
      fall = fall * arg
      arg = arg + 1
    end do

  end if

  return
end
subroutine gray_code_check ( n, t, ierror )

!*****************************************************************************80
!
!! GRAY_CODE_CHECK checks a Gray code element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of digits in each element.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(N), an element of the Gray code.
!    Each entry T(I) is either 0 or 1.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error, T represents a Gray code element.
!    -1, N is not positive.
!    I, error, T(I) is an illegal value for a Gray code element.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) t(n)

  ierror = 0

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAY_CODE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  do i = 1, n

    if ( t(i) /= 0 .and. t(i) /= 1 ) then
      ierror = i
      return
    end if

  end do

  return
end
subroutine gray_code_enum ( n, ngray )

!*****************************************************************************80
!
!! GRAY_CODE_ENUM enumerates the Gray codes on N digits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of digits in each element.
!    N must be nonnegative.
!
!    Output, integer ( kind = 4 ) NGRAY, the number of distinct elements.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ngray

  ngray = 2**n

  return
end
subroutine gray_code_rank ( n, t, rank )

!*****************************************************************************80
!
!! GRAY_CODE_RANK computes the rank of a Gray code element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of digits in each element.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(N), an element of the Gray code.
!    Each entry T(I) is either 0 or 1.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the element.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(n)
!
!  Check.
!
  call gray_code_check ( n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAY_CODE_RANK - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  rank = 0
  b = 0

  do i = n - 1, 0, -1

    if ( t(n-i) /= 0 ) then
      b = 1 - b
    end if

    if ( b == 1 ) then
      rank = rank + 2 ** i
    end if

  end do

  return
end
subroutine gray_code_successor ( n, t, rank )

!*****************************************************************************80
!
!! GRAY_CODE_SUCCESSOR computes the binary reflected Gray code successor.
!
!  Example:
!
!    000, 001, 011, 010, 110, 111, 101, 100,
!    after which the sequence repeats.
!
!  Discussion:
!
!    In the original code, the successor of the element that has an
!    initial 1 followed by N-1 zeroes is undefined.  In this version,
!    the successor is the element with N zeroes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of digits in each element.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) T(N).
!    On input, T contains an element of the Gray code, that is,
!    each entry T(I) is either 0 or 1.
!    On output, T contains the successor to the input value; this
!    is an element of the Gray code, which differs from the input
!    value in a single position.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(n)
  integer ( kind = 4 ) weight
!
!  Return the first element.
!
  if ( rank == -1 ) then

    t(1:n) = 0
    rank = 0
    return

  end if
!
!  Check.
!
  call gray_code_check ( n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAY_CODE_SUCCESSOR - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  weight = sum ( t(1:n) )

  if ( mod ( weight, 2 ) == 0 ) then

    if ( t(n) == 0 ) then
      t(n) = 1
    else
      t(n) = 0
    end if

    rank = rank + 1
    return

  else

    do i = n, 2, -1
      if ( t(i) == 1 ) then
        if ( t(i-1) == 0 ) then
          t(i-1) = 1
        else
          t(i-1) = 0
        end if
        rank = rank + 1
        return
      end if
    end do
!
!  The final element was input.
!  Return the first element.
!
    t(1:n) = 0
    rank = 0

  end if

  return
end
subroutine gray_code_unrank ( rank, n, t )

!*****************************************************************************80
!
!! GRAY_CODE_UNRANK computes the Gray code element of given rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the element.
!    0 <= RANK <= 2**N.
!
!    Input, integer ( kind = 4 ) N, the number of digits in each element.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) T(N), the element of the Gray code which has
!    the given rank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) b
  integer ( kind = 4 ) bprime
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ngray
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
  integer ( kind = 4 ) t(n)
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAY_CODE_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call gray_code_enum ( n, ngray )

  if ( rank < 0 .or. ngray < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAY_CODE_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  rank_copy = rank
  t(1:n) = 0
  bprime = 0

  do i = n - 1, 0, -1

    b = rank_copy / 2**i

    if ( b /= bprime ) then
      t(n-i) = 1
    end if

    bprime = b
    rank_copy = rank_copy - b * 2 ** i

  end do

  return
end
function i4_factorial ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL computes the factorial of an I4.
!
!  Formula:
!
!    Factorial ( N ) = Product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, I4_FACTORIAL is returned as 1.
!
!    Output, integer ( kind = 4 ) I4_FACTORIAL, the factorial of N.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_factorial
  integer ( kind = 4 ) n

  i4_factorial = 1

  do i = 1, n
    i4_factorial = i4_factorial * i
  end do

  return
end
subroutine i4_factorial_values ( n, x, fx )

!*****************************************************************************80
!
!! I4_FACTORIAL_VALUES returns values of the factorial function for testing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.
!    On input, if N is 0, the first test data is returned, and N is set
!    to the index of the test data.  On each subsequent call, N is
!    incremented and that test data is returned.  When there is no more
!    test data, N is set to 0.
!
!    Output, integer ( kind = 4 ) X, the argument of the function.
!
!    Output, integer ( kind = 4 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 13

  integer ( kind = 4 ), save, dimension ( nmax ) :: fxvec = (/ &
    1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, &
    3628800, 39916800, 479001600 /)
  integer ( kind = 4 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x
  integer ( kind = 4 ), save, dimension ( nmax ) :: xvec = (/ &
     0,  1,  2, 3, 4, 5, 6, 7, 8, 9, &
    10, 11, 12 /)

  if ( n < 0 ) then
    n = 0
  end if

  n = n + 1

  if ( nmax < n ) then
    n = 0
    x = 0
    fx = 0
    return
  end if

  x = xvec(n)
  fx = fxvec(n)

  return
end
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
!    use I4_HUGE() and HUGE interchangeably.
!
!    However, when using the G95, the values returned by HUGE were
!    not equal to 2147483647, apparently, and were causing severe
!    and obscure errors in my random number generator, which needs to
!    add I4_HUGE to the seed whenever the seed is negative.  So I
!    am backing away from invoking HUGE, whereas I4_HUGE is under
!    my control.
!
!    Explanation: because under G95 the default integer type is 64 bits!
!    So HUGE ( 1 ) = a very very huge integer indeed, whereas
!    I4_HUGE ( ) = the same old 32 bit big value.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Pierre LEcuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value,
!    which should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) I4_UNIFORM, a number between
!    A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine i4vec_backtrack ( n, x, indx, k, nstack, stack, maxstack, ncan )

!*****************************************************************************80
!
!! I4VEC_BACKTRACK supervises a backtrack search for an I4VEC.
!
!  Discussion:
!
!    The routine tries to construct an integer vector one index at a time,
!    using possible candidates as supplied by the user.
!
!    At any time, the partially constructed vector may be discovered to be
!    unsatisfactory, but the routine records information about where the
!    last arbitrary choice was made, so that the search can be
!    carried out efficiently, rather than starting out all over again.
!
!    First, call the routine with INDX = 0 so it can initialize itself.
!
!    Now, on each return from the routine, if INDX is:
!      1, you've just been handed a complete candidate vector;
!         Admire it, analyze it, do what you like.
!      2, please determine suitable candidates for position X(K).
!         Return the number of candidates in NCAN(K), adding each
!         candidate to the end of STACK, and increasing NSTACK.
!      3, you're done.  Stop calling the routine;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of positions to be filled in
!    the vector.
!
!    Input/output, integer ( kind = 4 ) X(N), the partial or complete candidate
!    vector.
!
!    Input/output, integer ( kind = 4 ) INDX, a communication flag.
!    On input,
!      0 to start a search.
!    On output:
!      1, a complete output vector has been determined and returned in X(1:N);
!      2, candidates are needed for position X(K);
!      3, no more possible vectors exist.
!
!    Output, integer ( kind = 4 ) K, if INDX=2, the current vector index
!    being considered.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!
!    Input, integer ( kind = 4 ) STACK(MAXSTACK), a list of all current
!    candidates for all positions 1 through K.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum length of the stack.
!
!    Input/output, integer ( kind = 4 ) NCAN(N), lists the current number of
!    candidates for positions 1 through K.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) indx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncan(n)
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) stack(maxstack)
  integer ( kind = 4 ) x(n)
!
!  If this is the first call, request a candidate for position 1.
!
  if ( indx == 0 ) then
    k = 1
    nstack = 0
    indx = 2
    return
  end if
!
!  Examine the stack.
!
  do
!
!  If there are candidates for position K, take the first available
!  one off the stack, and increment K.
!
!  This may cause K to reach the desired value of N, in which case
!  we need to signal the user that a complete set of candidates
!  is being returned.
!
    if ( 0 < ncan(k) ) then

      x(k) = stack(nstack)
      nstack = nstack - 1

      ncan(k) = ncan(k) - 1

      if ( k /= n ) then
        k = k + 1
        indx = 2
      else
        indx = 1
      end if

      exit
!
!  If there are no candidates for position K, then decrement K.
!  If K is still positive, repeat the examination of the stack.
!
    else

      k = k - 1

      if ( k <= 0 ) then
        indx = 3
        exit
      end if

    end if

  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_part1 ( n, npart, x )

!*****************************************************************************80
!
!! I4VEC_PART1 partitions an integer N into NPART parts.
!
!  Example:
!
!    Input:
!
!      N = 17, NPART = 5
!
!    Output:
!
!      X = ( 13, 1, 1, 1, 1 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.  N
!    may be positive, zero, or negative.
!
!    Input, integer ( kind = 4 ) NPART, the number of entries in the array.
!    1 <= NPART <= N.
!
!    Output, integer ( kind = 4 ) X(NPART), the partition of N.  The entries of
!    X add up to N.  X(1) = N + 1 - NPART, and all other entries
!    are equal to 1.
!
  implicit none

  integer ( kind = 4 ) npart

  integer ( kind = 4 ) n
  integer ( kind = 4 ) x(npart)

  if ( npart < 1 .or. n < npart ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_PART1 - Fatal error!'
    write ( *, '(a)' ) '  The input value of NPART is illegal.'
    stop
  end if

  x(1) = n + 1 - npart
  x(2:npart) = 1

  return
end
subroutine i4vec_part2 ( n, npart, x )

!*****************************************************************************80
!
!! I4VEC_PART2 partitions an integer N into NPART nearly equal parts.
!
!  Example:
!
!    Input:
!
!      N = 17, NPART = 5
!
!    Output:
!
!      X = ( 4, 4, 3, 3, 3 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.  N
!    may be positive, zero, or negative.
!
!    Input, integer ( kind = 4 ) NPART, the number of entries in the array.
!    1 <= NPART
!
!    Output, integer ( kind = 4 ) X(NPART), the partition of N.  The entries of
!    X add up to N.  The entries of X are either all equal, or
!    differ by at most 1.  The entries of X all have the same sign
!    as N, and the "largest" entries occur first.
!
  implicit none

  integer ( kind = 4 ) npart

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x(npart)

  if ( npart < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_PART2 - Fatal error!'
    write ( *, '(a)' ) '  The input value of NPART is illegal.'
    stop
  end if

  x(1:npart) = 0

  if ( 0 < n ) then

    j = 1
    do i = 1, n
      x(j) = x(j) + 1
      j = j + 1
      if ( npart < j ) then
        j = 1
      end if
    end do

  else if ( n < 0 ) then

    j = 1
    do i = n, - 1
      x(j) = x(j) - 1
      j = j + 1
      if ( npart < j ) then
        j = 1
      end if
    end do

  end if

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
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
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(i8,i10)' ) i, a(i)
  end do

  return
end
subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! I4VEC_REVERSE reverses the elements of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    In FORTRAN90, call I4VEC_REVERSE is equivalent to:
!
!      A(1:N) = A(N:1:-1)
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)

  a(1:n) = a(n:1:-1)

  return
end
subroutine i4vec_search_binary_a ( n, a, b, indx )

!*****************************************************************************80
!
!! I4VEC_SEARCH_BINARY_A searches the ascending sorted I4VEC for a value.
!
!  Discussion:
!
!    Binary search is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the array to be searched.  A must
!    be sorted in ascending order.
!
!    Input, integer ( kind = 4 ) B, the value to be searched for.
!
!    Output, integer ( kind = 4 ) INDX, the result of the search.
!    0, B does not occur in A.
!    I, A(I) = B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) high
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid

  indx = 0

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == b ) then
      indx = mid
      exit
    else if ( a(mid) < b ) then
      low = mid + 1
    else if ( b < a(mid) ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine i4vec_search_binary_d ( n, a, b, indx )

!*****************************************************************************80
!
!! I4VEC_SEARCH_BINARY_D searches a descending sorted I4VEC for a value.
!
!  Discussion:
!
!    Binary search is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the array to be searched.  A must
!    be sorted in descending order.
!
!    Input, integer ( kind = 4 ) B, the value to be searched for.
!
!    Output, integer ( kind = 4 ) INDX, the result of the search.
!    0, B does not occur in A.
!    I, A(I) = B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) high
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid

  indx = 0

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == b ) then
      indx = mid
      exit
    else if ( b < a(mid) ) then
      low = mid + 1
    else if ( a(mid) < b ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine i4vec_sort_insert_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_INSERT_A uses an ascending insertion sort on an I4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) A(N).
!
!    On input, A contains data to be sorted.
!    On output, the entries of A have been sorted in ascending order.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( a(j) <= x ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine i4vec_sort_insert_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_INSERT_D uses a descending insertion sort on an I4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) A(N).
!
!    On input, A contains data to be sorted.
!    On output, the entries of A have been sorted in descending order.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( x <= a(j) ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine knapsack_01 ( n, mass_limit, p, w, x, mass, profit )

!*****************************************************************************80
!
!! KNAPSACK_01 solves the 0/1 knapsack problem.
!
!  Discussion:
!
!    The 0/1 knapsack problem is as follows:
!
!      Given:
!        a set of N objects,
!        a profit P(I) and weight W(I) associated with each object,
!        and a weight limit MASS_LIMIT,
!      Determine:
!        a set of choices X(I) which are 0 or 1, that maximizes the profit
!          P = Sum ( 1 <= I <= N ) P(I) * X(I)
!        subject to the constraint
!          Sum ( 1 <= I <= N ) W(I) * X(I) <= MASS_LIMIT.
!
!    This routine assumes that the objects have already been sorted
!    in order of decreasing "profit density", P(I)/W(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, real ( kind = 8 ) MASS_LIMIT, the weight limit of the
!    chosen objects.
!
!    Input/output, real ( kind = 8 ) P(N), the "profit" or value of each object.
!    P is assumed to be nonnegative.
!
!    Input/output, real ( kind = 8 ) W(N), the "weight" or cost of each object.
!    W is assumed to be  nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the choice function for the objects.
!    0, the object was not taken.
!    1, the object was taken.
!
!    Output, real ( kind = 8 ) MASS, the total mass of the objects taken.
!
!    Output, real ( kind = 8 ) PROFIT, the total profit of the objects taken.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxstack = 100

  integer ( kind = 4 ) n

  integer ( kind = 4 ) indx
  integer ( kind = 4 ) k
  real ( kind = 8 ) mass
  real ( kind = 8 ) mass_1
  real ( kind = 8 ) mass_2
  real ( kind = 8 ) mass_best
  real ( kind = 8 ) mass_limit
  real ( kind = 8 ) mass_remaining
  integer ( kind = 4 ) ncan(n)
  integer ( kind = 4 ) nstack
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) profit
  real ( kind = 8 ) profit_1
  real ( kind = 8 ) profit_2
  real ( kind = 8 ) profit_best
  real ( kind = 8 ) stack(maxstack)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_best(n)

  nstack = 0
!
!  Initialize the "best so far" data.
!
  x_best(1:n) = 0.0D+00
  profit_best = 0.0D+00
  mass_best = 0
!
!  Begin the backtracking solution.
!
  indx = 0

  do

    call r8vec_backtrack ( n, x, indx, k, nstack, stack, maxstack, ncan )
!
!  Got a new candidate.  Compare it to the best so far.
!
    if ( indx == 1 ) then

      profit = dot_product ( p, x )
      mass = dot_product ( w, x )

      if ( profit_best < profit .or. &
         ( profit == profit_best .and. mass < mass_best ) ) then
        profit_best = profit
        mass_best = mass
        x_best(1:n) = x(1:n)
      end if
!
!  Need candidates for X(K).
!
!  X(K) = 1 is possible if:
!
!    * adding W(K) to our mass doesn't put us over our mass limit;
!    * and adding P(K) to our current profit, and taking the best we
!      could get using rational X for the remainder would put us over
!      our current best.
!
!  X(K) = 0 is always possible.
!
    else if ( indx == 2 ) then

      ncan(k) = 0

      mass_1 = dot_product ( w(1:k-1), x(1:k-1) ) + w(k)

      if ( mass_1 <= mass_limit ) then

        mass_remaining = mass_limit - mass_1

        profit_1 = dot_product ( p(1:k-1), x(1:k-1) ) + p(k)

        if ( k < n ) then
          call knapsack_rational ( n - k, mass_remaining, p(k+1), w(k+1), &
            x(k+1), mass_2, profit_2 )
        else
          profit_2 = 0.0D+00
        end if

        if ( profit_best < profit_1 + profit_2 ) then
          if ( maxstack <= nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'KNAPSACK_01 - Fatal error!'
            write ( *, '(a)' ) '  Exceeded stack space.'
            return
          end if
          ncan(k) = ncan(k) + 1
          nstack = nstack + 1
          stack(nstack) = 1.0D+00
        end if

      end if

      if ( maxstack <= nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KNAPSACK_01 - Fatal error!'
        write ( *, '(a)' ) '  Exceeded stack space.'
        return
      end if

      ncan(k) = ncan(k) + 1
      nstack = nstack + 1
      stack(nstack) = 0.0D+00
!
!  Done.  Return the best solution.
!
    else

      profit = profit_best
      mass = mass_best
      x(1:n) = x_best(1:n)
      exit

    end if

  end do

  return
end
subroutine knapsack_rational ( n, mass_limit, p, w, x, mass, profit )

!*****************************************************************************80
!
!! KNAPSACK_RATIONAL solves the rational knapsack problem.
!
!  Discussion:
!
!    The rational knapsack problem is a generalization of the 0/1 knapsack
!    problem.  It is mainly used to derive a bounding function for the
!    0/1 knapsack problem.
!
!    The 0/1 knapsack problem is as follows:
!
!      Given:
!        a set of N objects,
!        a profit P(I) and weight W(I) associated with each object,
!        and a weight limit MASS_LIMIT,
!      Determine:
!        a set of choices X(I) which are 0 or 1, that maximizes the profit
!          P = Sum ( 1 <= I <= N ) P(I) * X(I)
!        subject to the constraint
!          Sum ( 1 <= I <= N ) W(I) * X(I) <= MASS_LIMIT.
!
!    By contrast, the rational knapsack problem allows the values X(I)
!    to be any value between 0 and 1.  A solution for the rational knapsack
!    problem is known.  Arrange the objects in order of their "profit density"
!    ratios P(I)/W(I), and then take in order as many of these as you can.
!    If you still have "room" in the weight constraint, then you should
!    take the maximal fraction of the very next object, which will complete
!    your weight limit, and maximize your profit.
!
!    If should be obvious that, given the same data, a solution for
!    the rational knapsack problem will always have a profit that is
!    at least as high as for the 0/1 problem.  Since the rational knapsack
!    maximum profit is easily computed, this makes it a useful bounding
!    function.
!
!    Note that this routine assumes that the objects have already been
!    arranged in order of the "profit density".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, real ( kind = 8 ) MASS_LIMIT, the weight limit of the
!    chosen objects.
!
!    Input, real ( kind = 8 ) P(N), the "profit" or value of each object.
!    The entries of P are assumed to be nonnegative.
!
!    Input, real ( kind = 8 ) W(N), the "weight" or cost of each object.
!    The entries of W are assumed to be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the choice function for the objects.
!    0.0, the object was not taken.
!    1.0, the object was taken.
!    R, where 0 < R < 1, a fractional amount of the object was taken.
!
!    Output, real ( kind = 8 ) MASS, the total mass of the objects taken.
!
!    Output, real ( kind = 8 ) PROFIT, the total profit of the objects taken.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) mass
  real ( kind = 8 ) mass_limit
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) profit
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  mass = 0.0D+00
  profit = 0.0D+00

  do i = 1, n

    if ( mass_limit <= mass ) then
      x(i) = 0.0D+00
    else if ( mass + w(i) <= mass_limit ) then
      x(i) = 1.0D+00
      mass = mass + w(i)
      profit = profit + p(i)
    else
      x(i) = ( mass_limit - mass ) / w(i)
      mass = mass_limit
      profit = profit + p(i) * x(i)
    end if

  end do

  return
end
subroutine knapsack_reorder ( n, p, w )

!*****************************************************************************80
!
!! KNAPSACK_REORDER reorders the knapsack data by "profit density".
!
!  Discussion:
!
!    This routine must be called to rearrange the data before calling
!    routines that handle a knapsack problem.
!
!    The "profit density" for object I is defined as P(I)/W(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, real ( kind = 8 ) P(N), the "profit" or value of each object.
!
!    Input/output, real ( kind = 8 ) W(N), the "weight" or cost of each object.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) t
  real ( kind = 8 ) w(n)
!
!  Rearrange the objects in order of "profit density".
!
  do i = 1, n
    do j = i + 1, n
      if ( p(i) * w(j) < p(j) * w(i) ) then

        t    = p(i)
        p(i) = p(j)
        p(j) = t

        t    = w(i)
        w(i) = w(j)
        w(j) = t

      end if
    end do
  end do

  return
end
subroutine ksubset_colex_check ( k, n, t, ierror )

!*****************************************************************************80
!
!! KSUBSET_COLEX_CHECK checks a K subset in colex form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the number of elements each K subset must
!    have. 1 <= K <= N.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(K), describes a K subset.  T(I) is the I-th
!    element of the K subset.  The elements must be listed in
!    DESCENDING order.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    -1, N is not positive.
!    -2, K is not positive.
!    I, entry I is illegal.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t(k)
  integer ( kind = 4 ) tmax

  ierror = 0

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_COLEX_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( k < 1 .or. n < k ) then
    ierror = -2
    return
  end if

  tmax = n + 1

  do i = 1, k

    if ( t(i) <= 0 .or. tmax <= t(i) ) then
      ierror = i
      return
    end if

    tmax = t(i)

  end do

  return
end
subroutine ksubset_colex_rank ( k, n, t, rank )

!*****************************************************************************80
!
!! KSUBSET_COLEX_RANK computes the colex rank of a K subset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the number of elements each K subset must
!    have.  1 <= K <= N.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(K), describes a K subset.  T(I) is the I-th
!    element of the K subset.  The elements must be listed in DESCENDING order.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the subset.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(k)
!
!  Check.
!
  call ksubset_colex_check ( k, n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_COLEX_CHECK - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  rank = 0

  do i = 1, k
    rank = rank + binomial ( t(i) - 1, k + 1 - i )
  end do

  return
end
subroutine ksubset_colex_successor ( k, n, t, rank )

!*****************************************************************************80
!
!! KSUBSET_COLEX_SUCCESSOR computes the K subset colex successor.
!
!  Discussion:
!
!    In the original code, there is a last element with no successor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the number of elements each K subset must
!    have.  1 <= K <= N.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) T(K), describes a K subset.  T(I) is the
!    I-th element.  The elements must be listed in DESCENDING order.
!    On input, T describes a K subset.
!    On output, T describes the next K subset in the ordering.
!    If the input T was the last in the ordering, then the output T
!    will be the first.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(k)
!
!  Return the first element.
!
  if ( rank == -1 ) then
    do i = 1, k
      t(i) = k + 1 - i
    end do
    rank = 0
    return
  end if
!
!  Check.
!
  call ksubset_colex_check ( k, n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_COLEX_SUCCESSOR - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  do i = k - 1, 1, -1
    if ( t(k+1-i) + 1 < t(k-i) ) then
      t(k+1-i) = t(k+1-i) + 1
      rank = rank + 1
      return
    end if
  end do

  if ( t(1) < n ) then
    t(1) = t(1) + 1
    do i = 1, k - 1
      t(k+1-i) = i
    end do
    rank = rank + 1
    return
  end if
!
!  The last K subset was input.
!  Return the first one.
!
  do i = 1, k
    t(i) = k + 1 - i
  end do

  rank = 0

  return
end
subroutine ksubset_colex_unrank ( rank, k, n, t )

!*****************************************************************************80
!
!! KSUBSET_COLEX_UNRANK computes the K subset of given colex rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the K subset.
!
!    Input, integer ( kind = 4 ) K, the number of elements each K subset must
!    have.  1 <= K <= N.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) T(K), describes the K subset of the given
!    rank.  T(I) is the I-th element.  The elements must be listed in
!    DESCENDING order.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nksub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
  integer ( kind = 4 ) t(k)
  integer ( kind = 4 ) x
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_COLEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  if ( k < 1 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_COLEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input K is illegal.'
    stop
  end if

  call ksubset_enum ( k, n, nksub )

  if ( rank < 0 .or. nksub < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_COLEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if
!
  rank_copy = rank

  x = n

  do i = 1, k

    do while ( rank_copy < binomial ( x, k + 1 - i ) )
      x = x - 1
    end do

    t(i) = x + 1
    rank_copy = rank_copy - binomial ( x, k + 1 - i )

  end do

  return
end
subroutine ksubset_enum ( k, n, nksub )

!*****************************************************************************80
!
!! KSUBSET_ENUM enumerates the K element subsets of an N set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the number of elements each K subset must
!    have. 0 <= K <= N.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    0 <= N.
!
!    Output, integer ( kind = 4 ) NKSUB, the number of distinct elements.
!
  implicit none

  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nksub

  nksub = binomial ( n, k )

  return
end
subroutine ksubset_lex_check ( k, n, t, ierror )

!*****************************************************************************80
!
!! KSUBSET_LEX_CHECK checks a K subset in lex form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the number of elements each K subset must
!    have. 1 <= K <= N.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(K), describes a K subset.  T(I) is the I-th
!    element of the K subset.  The elements must be listed in
!    DESCENDING order.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    -1, N is illegal.
!    -2, K is illegal.
!    I, entry I is illegal.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t(k)
  integer ( kind = 4 ) tmin

  ierror = 0

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_LEX_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( k < 1 .or. n < k ) then
    ierror = -2
    return
  end if

  tmin = 0

  do i = 1, k

    if ( t(i) <= tmin .or. n < t(i) ) then
      ierror = i
      return
    end if

    tmin = t(i)

  end do

  return
end
subroutine ksubset_lex_rank ( k, n, t, rank )

!*****************************************************************************80
!
!! KSUBSET_LEX_RANK computes the lexicographic rank of a K subset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the number of elements each K subset must
!    have.  1 <= K <= N.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(K), describes a K subset.  T(I) is the I-th
!    element.  The elements must be listed in ascending order.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the K subset.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(k)
  integer ( kind = 4 ) tim1
!
!  Check.
!
  call ksubset_lex_check ( k, n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_LEX_RANK - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    stop
  end if

  rank = 0

  do i = 1, k

    if ( i == 1 ) then
      tim1 = 0
    else
      tim1 = t(i-1)
    end if

    if ( tim1 + 1 <= t(i) - 1 ) then
      do j = tim1 + 1, t(i) - 1
        rank = rank + binomial ( n - j, k - i )
      end do
    end if

  end do

  return
end
subroutine ksubset_lex_successor ( k, n, t, rank )

!*****************************************************************************80
!
!! KSUBSET_LEX_SUCCESSOR computes the K subset lexicographic successor.
!
!  Discussion:
!
!    In the original code, there is a last element with no successor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the number of elements each K subset must
!    have. 1 <= K <= N.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) T(K), describes a K subset.  T(I) is
!    the I-th element.  The elements must be listed in ascending order.
!    On input, T describes a K subset.
!    On output, T describes the next K subset in the ordering.
!    If the input T was the last in the ordering, then the output T
!    will be the first.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isave
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(k)
!
!  Return the first element.
!
  if ( rank == -1 ) then
    call i4vec_indicator ( k, t )
    rank = 0
    return
  end if
!
!  Check.
!
  call ksubset_lex_check ( k, n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_LEX_SUCCESSOR - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  isave = 0

  do i = k, 1, -1
    if ( t(i) /= n - k + i ) then
      isave = i
      exit
    end if
  end do
!
!  The last K subset was input.
!  Return the first one.
!
  if ( isave == 0 ) then
    call i4vec_indicator ( k, t )
    rank = 0
  else

    do j = k, isave, -1
      t(j) = t(isave) + 1 + j - isave
    end do

    rank = rank + 1

  end if

  return
end
subroutine ksubset_lex_unrank ( rank, k, n, t )

!*****************************************************************************80
!
!! KSUBSET_LEX_UNRANK computes the K subset of given lexicographic rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the K subset.
!
!    Input, integer ( kind = 4 ) K, the number of elements each K subset must
!    have.  1 <= K <= N.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) T(K), describes the K subset of the given
!    rank.  T(I) is the I-th element.  The elements must be listed in
!    ascending order.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nksub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
  integer ( kind = 4 ) t(k)
  integer ( kind = 4 ) x
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  if ( k < 1 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input K is illegal.'
    stop
  end if

  call ksubset_enum ( k, n, nksub )

  if ( rank < 0 .or. nksub < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input rank is illegal.'
    stop
  end if

  rank_copy = rank

  x = 1

  do i = 1, k

    do while ( binomial ( n - x, k - i ) <= rank_copy )
      rank_copy = rank_copy - binomial ( n - x, k - i )
      x = x + 1
    end do

    t(i) = x
    x = x + 1

  end do

  return
end
subroutine ksubset_revdoor_rank ( k, n, t, rank )

!*****************************************************************************80
!
!! KSUBSET_REVDOOR_RANK computes the revolving door rank of a K subset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the number of elements each K subset must
!    have.  1 <= K <= N.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(K), describes a K subset.  T(I) is the I-th
!    element.  The elements must be listed in ascending order.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the K subset.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) s
  integer ( kind = 4 ) t(k)
!
!  Check.
!
  call ksubset_lex_check ( k, n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_REVDOOR_RANK - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  if ( mod ( k, 2 ) == 0 ) then

    rank = 0

  else

    rank = - 1

  end if

  s = 1

  do i = k, 1, -1
    rank = rank + s * binomial ( t(i), i )
    s = - s
  end do

  return
end
subroutine ksubset_revdoor_successor ( k, n, t, rank )

!*****************************************************************************80
!
!! KSUBSET_REVDOOR_SUCCESSOR computes the K subset revolving door successor.
!
!  Discussion:
!
!    After numerous attempts to implement the algorithm published in
!    Kreher and Stinson, the Nijenhuis and Wilf version was implemented
!    instead.  The K and S algorithm is supposedly based on the N and W one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the number of elements each K subset must
!    have.  1 <= K <= N.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) T(K), describes a K subset.  T(I) is the
!    I-th element.  The elements must be listed in ascending order.
!    On input, T describes a K subset.
!    On output, T describes the next K subset in the ordering.
!    If the input T was the last in the ordering, then the output T
!    will be the first.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(k)
!
!  Return the first element.
!
  if ( rank == - 1 ) then
    call i4vec_indicator ( k, t )
    rank = 0
    return
  end if
!
!  Check.
!
  call ksubset_lex_check ( k, n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_REVDOOR_SUCCESSOR - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  j = 0

  do

    if ( 0 < j .or. mod ( k, 2 ) == 0 ) then

      j = j + 1

      if ( k < j ) then
        t(k) = k
        rank = 0
        return
      end if

      if ( t(j) /= j ) then

        t(j) = t(j) - 1

        if ( j /= 1 ) then
          t(j-1) = j - 1
        end if

        rank = rank + 1
        return

      end if

    end if

    j = j + 1

    if ( j < k ) then
      if ( t(j) /= t(j+1) - 1 ) then
        exit
      end if
    else
      if ( t(j) /= n ) then
        exit
      end if
    end if

  end do

  t(j) = t(j) + 1

  if ( j /= 1 ) then
    t(j-1) = t(j) - 1
  end if

  rank = rank + 1

  return
end
subroutine ksubset_revdoor_unrank ( rank, k, n, t )

!*****************************************************************************80
!
!! KSUBSET_REVDOOR_UNRANK computes the K subset of given revolving door rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the K subset.
!
!    Input, integer ( kind = 4 ) K, the number of elements each K subset must
!    have.  1 <= K <= N.
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) T(K), describes the K subset of the given
!    rank.  T(I) is the I-th element.  The elements must be listed in
!    ascending order.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nksub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
  integer ( kind = 4 ) t(k)
  integer ( kind = 4 ) x
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_REVDOOR_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  if ( k < 1 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_REVDOOR_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input K is illegal.'
    stop
  end if

  call ksubset_enum ( k, n, nksub )

  if ( rank < 0 .or. nksub < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_REVDOOR_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  rank_copy = rank

  x = n

  do i = k, 1, -1

    do while ( rank_copy < binomial ( x, i ) )
      x = x - 1
    end do

    t(i) = x + 1
    rank_copy = binomial ( x + 1, i ) - rank_copy - 1

  end do

  return
end
subroutine marriage ( n, prefer, rank, fiancee, next )

!*****************************************************************************80
!
!! MARRIAGE finds a stable set of marriages for given preferences.
!
!  Discussion:
!
!    Given a set of N men and N women who must be married in pairs,
!    and information defining the relative rankings that each person
!    assigns to the candidates of the opposite sex, this routine finds
!    a stable set of marriages for them.
!
!    A stable set of marriages is a pairing of the men and women with
!    the stability property: if M1 marries W1 and M2 marries W2, then
!    it is never the case that M1 and W2 would both prefer to be married
!    to each other.
!
!    An important application of stable marriage algorithms occurs in
!    the annual matching of medical residents to hospitals.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Sedgewick,
!    Algorithms in C,
!    Addison-Wesley, 1990,
!    ISBN: 0-201-51425-7,
!    LC: QA76.73.C15S43.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of pairs of men and women.
!
!    Input, integer ( kind = 4 ) PREFER(N,N); for man I, the value of
!    PREFER(I,J) represents his J-th preference for a wife.
!
!    Input, integer ( kind = 4 ) RANK(N,N); for woman I, the value of RANK(I,J)
!    represents her ranking of man number J.  A value of 1 for RANK(I,J)
!    means woman I ranks man J most preferable, while a value of N
!    would mean she ranked him least preferable.
!
!    Output, integer ( kind = 4 ) FIANCEE(N); for woman I, FIANCEE(I) is the
!    man to whom she is now engaged.
!
!    Output, integer ( kind = 4 ) NEXT(N); for man I, NEXT(I) is his preference
!    ranking for the woman to whom he is now engaged.  A value of 1 represents
!    his first choice, a value of N his last.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) fiancee(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) next(n)
  integer ( kind = 4 ) prefer(n,n)
  integer ( kind = 4 ) rank(n,n)
  integer ( kind = 4 ) temp
  integer ( kind = 4 ) w
!
!  For man I, NEXT(I) is the woman I has most recently proposed to,
!  and hence NEXT(I)+1 is the next one to try.
!
  next(1:n) = 0
!
!  For woman I, FIANCEE(I) is the man she has agree to marry,
!  or 0 if she has not agreed to any man yet.
!
  fiancee(1:n) = 0
!
!  Start with an unengaged man, and end with an engaged woman.
!
  do i = 1, n

    m = i

    do

      next(m) = next(m) + 1

      w = prefer(m,next(m))

      if ( fiancee(w) == 0 ) then
        fiancee(w) = m
        exit
      end if

      if ( rank (w,m) < rank(w,fiancee(w)) ) then
        temp       = fiancee(w)
        fiancee(w) = m
        m          = temp
      end if

    end do

  end do

  return
end
subroutine mountain ( n, x, y, mxy )

!*****************************************************************************80
!
!! MOUNTAIN enumerates the mountains.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, ...
!    N must be positive.
!
!    Input, integer ( kind = 4 ) X, Y, ...
!    0 <= X <= 2 * N,
!
!    Output, integer ( kind = 4 ) MXY, the value of the "mountain function"
!    M ( N, X, Y ), which is the number of all mountain ranges from
!    (X,Y) to (2*N,0) which do not drop below sea level.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) c
  integer ( kind = 4 ) mxy
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
!
!  Check.
!
  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MOUNTAIN - Fatal error!'
    write ( *, '(a)' ) '  N <= 0.'
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  else if ( x < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MOUNTAIN - Fatal error!'
    write ( *, '(a)' ) '  X < 0.'
    write ( *, '(a,i8)' ) '  X = ', x
    stop
  else if ( 2 * n < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MOUNTAIN - Fatal error!'
    write ( *, '(a)' ) '  2 * N < X.'
    write ( *, '(a,i8)' ) '  X = ', x
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  end if
!
!  Special cases.
!
  if ( y < 0 ) then
    mxy = 0
    return
  end if

  if ( 2 * n < x + y ) then
    mxy = 0
    return
  end if

  if ( mod ( x + y, 2 ) == 1 ) then
    mxy = 0
    return
  end if

  a = 2 * n - x
  b = n - ( x + y ) / 2
  c = n - 1 - ( x + y ) / 2

  mxy = binomial ( a, b ) - binomial ( a, c )

  return
end
subroutine npart_enum ( n, npart, npartitions )

!*****************************************************************************80
!
!! NPART_ENUM enumerates the number of partitions of N with NPART parts.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    Normally N must be positive, but for this routine any
!    N is allowed.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    Normally, 1 <= NPART <= N is required,
!    but for this routine any value of NPART is allowed.
!
!    Output, integer ( kind = 4 ) NPARTITIONS is the number of partitions of N
!    with NPART parts.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) npart
  integer ( kind = 4 ) npartitions
  integer ( kind = 4 ) p(0:n,0:n)

  if ( n <= 0 ) then

    npartitions = 0

  else if ( npart <= 0 .or. n < npart ) then

    npartitions = 0

  else

    call npart_table ( n, npart, n, p )

    npartitions = p(n,npart)

  end if

  return
end
subroutine npart_rsf_lex_random ( n, npart, seed, a )

!*****************************************************************************80
!
!! NPART_RSF_LEX_RANDOM returns a random RSF NPART partition.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) A(NPART), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) npart

  integer ( kind = 4 ) a(npart)
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) npartitions
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) seed

  call npart_enum ( n, npart, npartitions )

  rank = i4_uniform ( 1, npartitions, seed )

  call npart_rsf_lex_unrank ( rank, n, npart, a )

  return
end
subroutine npart_rsf_lex_rank ( n, npart, a, rank )

!*****************************************************************************80
!
!! NPART_RSF_LEX_RANK computes the lex rank of an RSF NPART partition.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Input, integer ( kind = 4 ) A(NPART), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the partition.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) npart

  integer ( kind = 4 ) a(npart)
  integer ( kind = 4 ) b(npart)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) npartcopy
  integer ( kind = 4 ) p(0:n,0:npart)
  integer ( kind = 4 ) rank
!
!  Check.
!
  call part_rsf_check ( n, npart, a, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NPART_RSF_LEX_RANK - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  Get the table of partitions of N with NPART parts.
!
  call npart_table ( n, npart, n, p )
!
!  Copy the partition "backwards".
!
  do i = 1, npart
    b(i) = a(npart+1-i)
  end do

  rank = 0
  ncopy = n
  npartcopy = npart

  do while ( 0 < ncopy .and. 0 < npartcopy )

    if ( b(npartcopy) == 1 ) then

      ncopy = ncopy - 1
      npartcopy = npartcopy - 1

    else

      do i = 1, npartcopy
        b(i) = b(i) - 1
      end do
      rank = rank + p(ncopy-1,npartcopy-1)
      ncopy = ncopy - npartcopy

    end if

  end do

  return
end
subroutine npart_rsf_lex_successor ( n, npart, a, rank )

!*****************************************************************************80
!
!! NPART_RSF_LEX_SUCCESSOR computes the RSF lex successor for NPART partitions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be at least 1.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Input/output, integer ( kind = 4 ) A(NPART), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) npart

  integer ( kind = 4 ) a(npart)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rank
!
!  Return the first element.
!
  if ( rank == -1 ) then

    if ( npart < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NPART_RSF_LEX_SUCCESSOR - Fatal error!'
      write ( *, '(a)' ) '  NPART < 1.'
      stop
    end if

    a(1:npart-1) = 1
    a(npart) = n - ( npart - 1 )

    rank = 0
    return

  end if
!
!  Check.
!
  call part_rsf_check ( n, npart, a, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NPART_RSF_LEX_SUCCESSOR - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  Find the first index I for which A(NPART+1-I) + 1 < A(NPART).
!
  i = 2

  do

    if ( npart < i ) then
      exit
    end if

    if ( a(npart+1-i) + 1 < a(npart) ) then
      exit
    end if

    i = i + 1

  end do
!
!  If no such index, we've reached the end of the line.
!
  if ( i == npart + 1 ) then

    a(1:npart-1) = 1
    a(npart) = n - ( npart - 1 )

    rank = 0
    return
!
!  Otherwise, increment A(NPART+1-I), and adjust other entries.
!
  else

    a(npart+1-i) = a(npart+1-i) + 1
    d = - 1

    do j = i - 1, 2, -1
      d = d + a(npart+1-j) - a(npart+1-i)
      a(npart+1-j) = a(npart+1-i)
    end do

    a(npart) = a(npart) + d

  end if

  rank = rank + 1

  return
end
subroutine npart_rsf_lex_unrank ( rank, n, npart, a )

!*****************************************************************************80
!
!! NPART_RSF_LEX_UNRANK unranks an RSF NPART partition in the lex ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the partition.
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Output, integer ( kind = 4 ) A(NPART), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) npart

  integer ( kind = 4 ) a(npart)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) npartcopy
  integer ( kind = 4 ) npartitions
  integer ( kind = 4 ) p(0:n,0:npart)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
!
!  Check.
!
  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NPART_RSF_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input N is illegal.'
    stop
  end if

  if ( npart < 1 .or. n < npart ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NPART_RSF_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input NPART is illegal.'
    stop
  end if

  call npart_enum ( n, npart, npartitions )

  if ( rank < 0 .or. npartitions < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NPART_RSF_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if
!
!  Get the table of partitions of N with NPART parts.
!
  call npart_table ( n, npart, n, p )

  a(1:npart) = 0

  rank_copy = rank
  ncopy = n
  npartcopy = npart

  do while ( 0 < ncopy )

    if ( rank_copy < p(ncopy-1,npartcopy-1) ) then
      a(npart+1-npartcopy) = a(npart+1-npartcopy) + 1
      ncopy = ncopy - 1
      npartcopy = npartcopy - 1
    else
      do i = 1, npartcopy
        a(npart+1-i) = a(npart+1-i) + 1
      end do
      rank_copy = rank_copy - p(ncopy-1,npartcopy-1)
      ncopy = ncopy - npartcopy
    end if

  end do

  return
end
subroutine npart_sf_lex_successor ( n, npart, a, rank )

!*****************************************************************************80
!
!! NPART_SF_LEX_SUCCESSOR computes SF NPART partition.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Input/output, integer ( kind = 4 ) A(NPART), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.  The values in A must be in DESCENDING order.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) npart

  integer ( kind = 4 ) a(npart)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) temp
!
!  Return the first element.
!
  if ( rank == -1 ) then
    call i4vec_part2 ( n, npart, a )
    rank = 0
    return
  end if
!
!  Check.
!
  call part_sf_check ( n, npart, a, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NPART_SF_LEX_SUCCESSOR - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  Find the last entry that is 2 or more.
!
  do i = npart, 1, - 1
    if ( 1 < a(i) ) then
      indx = i
      exit
    end if
  end do
!
!  As long as the last nonunit occurs after the first position,
!  have it donate 1 to the left.
!
  if ( 1 < indx ) then

    a(indx) = a(indx) - 1
    a(indx-1) = a(indx-1) + 1
    indx = indx - 1

    do

      if ( indx <= 1 ) then
        exit
      end if

      if ( a(indx) <= a(indx-1) ) then
        exit
      end if

      temp      = a(indx)
      a(indx)   = a(indx-1)
      a(indx-1) = temp

      indx = indx - 1

    end do
!
!  Sum the tail.
!
    temp = sum ( a(indx+1:npart) )
!
!  Partition the tail sum equally over the tail.
!
    call i4vec_part2 ( temp, npart - indx, a(indx+1) )

    rank = rank + 1
!
!  If A(2) through A(NPART) are 1, then this is the last element.
!  Return the first one.
!
  else

    call i4vec_part2 ( n, npart, a )
    rank = 0

  end if

  return
end
subroutine npart_table ( n, npart, nmax, p )

!*****************************************************************************80
!
!! NPART_TABLE tabulates the number of partitions of N having NPART parts.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Input, integer ( kind = 4 ) NMAX, the leading dimension of P.
!
!    Output, integer ( kind = 4 ) P(0:NMAX,0:NPART), P(I,J) is the number of
!    partitions of I having J parts.
!
  implicit none

  integer ( kind = 4 ) nmax
  integer ( kind = 4 ) npart

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p(0:nmax,0:npart)

  p(0,0) = 1
  p(1:n,0) = 0

  do i = 1, n
    do j = 1, npart
      if ( i < j ) then
        p(i,j) = 0
      else if ( i < 2 * j ) then
        p(i,j) = p(i-1,j-1)
      else
        p(i,j) = p(i-1,j-1) + p(i-j,j)
      end if
    end do
  end do

  return
end
subroutine part_enum ( n, npartitions )

!*****************************************************************************80
!
!! PART_ENUM enumerates the number of partitions of N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    Normally N must be positive, but for this routine any
!    N is allowed.
!
!    Output, integer ( kind = 4 ) NPARTITIONS is the number of partitions of N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) npartitions
  integer ( kind = 4 ) p(0:n)

  if ( n < 0 ) then

    npartitions = 0

  else

    call part_table ( n, p )

    npartitions = p(n)

  end if

  return
end
subroutine part_rsf_check ( n, npart, a, ierror )

!*****************************************************************************80
!
!! PART_RSF_CHECK checks a reverse standard form partition of an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Input, integer ( kind = 4 ) A(NPART), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.  The entries must be in ASCENDING order.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    -1, N is illegal.
!    -2, NPART is illegal.
!    -3, the entries do not add up to N.
!    I, the I-th entry of A is illegal.
!
  implicit none

  integer ( kind = 4 ) npart

  integer ( kind = 4 ) a(npart)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n

  ierror = 0

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PART_RSF_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( npart < 1 .or. n < npart ) then
    ierror = -2
    return
  end if
!
!  Every entry must lie between 1 and N.
!
  do i = 1, npart
    if ( a(i) < 1 .or. n < a(i) ) then
      ierror = i
      return
    end if
  end do
!
!  The entries must be in ascending order.
!
  do i = 2, npart
    if ( a(i) < a(i-1) ) then
      ierror = i
      return
    end if
  end do
!
!  The entries must add up to N.
!
  if ( sum ( a(1:npart) ) /= n ) then
    ierror = -3
  end if

  return
end
subroutine part_sf_check ( n, npart, a, ierror )

!*****************************************************************************80
!
!! PART_SF_CHECK checks a standard form partition of an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Input, integer ( kind = 4 ) A(NPART), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.  The entries must be in DESCENDING order.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    -1, N is illegal.
!    -2, NPART is illegal.
!    -3, the entries do not add up to N.
!    I, the I-th entry of A is illegal.
!
  implicit none

  integer ( kind = 4 ) npart

  integer ( kind = 4 ) a(npart)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n

  ierror = 0

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PART_SF_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( npart < 1 .or. n < npart ) then
    ierror = -2
    return
  end if
!
!  Every entry must lie between 1 and N.
!
  do i = 1, npart
    if ( a(i) < 1 .or. n < a(i) ) then
      ierror = i
      return
    end if
  end do
!
!  The entries must be in descending order.
!
  do i = 2, npart
    if ( a(i-1) < a(i) ) then
      ierror = i
      return
    end if
  end do
!
!  The entries must add up to N.
!
  if ( sum ( a(1:npart) ) /= n ) then
    ierror = -3
  end if

  return
end
subroutine part_sf_conjugate ( n, npart, a, npart2, b )

!*****************************************************************************80
!
!! PART_SF_CONJUGATE computes the conjugate of a partition.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Input, integer ( kind = 4 ) A(N), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.
!
!    Output, integer ( kind = 4 ) NPART2, the number of parts of the conjugate
!    partition.
!
!    Output, integer ( kind = 4 ) B(N), contains the conjugate partition.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) npart2
!
!  Check.
!
  call part_sf_check ( n, npart, a, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PART_SF_CONJUGATE - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  npart2 = a(1)
  b(1:npart2) = 0

  do i = 1, npart
    do j = 1, a(i)
      b(j) = b(j) + 1
    end do
  end do

  return
end
subroutine part_sf_majorize ( n, nparta, a, npartb, b, result )

!*****************************************************************************80
!
!! PART_SF_MAJORIZE determines if partition A majorizes partition B.
!
!  Discussion:
!
!    The partitions must be in standard form.
!
!    If A, with NPARTA parts, and B, with NPARTB parts, are both partitions
!    of the same positive integer N, then we say that A majorizes B if,
!    for every index K from 1 to N, it is true that
!
!      sum ( 1 <= I <= K ) B(I) <= sum ( 1 <= I <= K ) A(I)
!
!    where entries of A beyond index NPARTA, and of B beyond BPARTB
!    are assumed to be 0.  We say that A strictly majorizes B if
!    A majorizes B, and for at least one index K the inequality is strict.
!
!    For any two partitions of N, it is possible that A majorizes B,
!    B majorizes A, both partitions majorize each other (in which case
!    they are equal), or that neither majorizes the other.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack vanLint, Richard Wilson,
!    A Course in Combinatorics,
!    Cambridge, 1992,
!    ISBN: 0-521-42260-4,
!    LC: QA164.L56.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NPARTA, the number of parts in partition A.
!    1 <= NPARTA <= N.
!
!    Input, integer ( kind = 4 ) A(NPARTA), contains partition A in standard
!    form.  A(1) through A(NPARTA) contain nonzero integers which sum to N.
!
!    Input, integer ( kind = 4 ) NPARTB, the number of parts in partition B.
!    1 <= NPARTB <= N.
!
!    Input, integer ( kind = 4 ) B(NPARTB), contains partition B in standard
!    form.  B(1) through B(NPARTB) contain nonzero integers which sum to N.
!
!    Output, integer ( kind = 4 ) RESULT, the result of the comparison.
!    -2, A and B are incomparable, but would have been +1.
!    -1, A < B, (A is strictly majorized by B),
!     0, A = B, (A and B are identical),
!    +1, A > B, (A strictly majorizes B),
!    +2, A and B are incomparable, but would have been +1.
!
  implicit none

  integer ( kind = 4 ) nparta
  integer ( kind = 4 ) npartb

  integer ( kind = 4 ) a(nparta)
  integer ( kind = 4 ) b(npartb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n
  integer ( kind = 4 ) result
  integer ( kind = 4 ) suma
  integer ( kind = 4 ) sumb
!
!  Check.
!
  call part_sf_check ( n, nparta, a, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PART_SF_MAJORIZE - Fatal error!'
    write ( *, '(a)' ) '  The input array A is illegal.'
    stop
  end if

  call part_sf_check ( n, npartb, b, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PART_SF_MAJORIZE - Fatal error!'
    write ( *, '(a)' ) '  The input array B is illegal.'
    stop
  end if

  result = 0
  suma = 0
  sumb = 0

  do i = 1, min ( nparta, npartb )

    if ( i <= nparta ) then
      suma = suma + a(i)
    end if

    if ( i <= npartb ) then
      sumb = sumb + b(i)
    end if

    if ( result == -1 ) then

      if ( sumb < suma ) then
        result = -2
        return
      end if

    else if ( result == 0 ) then

      if ( suma < sumb ) then
        result = -1
      else if ( sumb < suma ) then
        result = +1
      end if

    else if ( result == + 1 ) then

      if ( suma < sumb ) then
        result = +2
        return
      end if

    end if

  end do

  return
end
subroutine part_successor ( n, npart, a, rank )

!*****************************************************************************80
!
!! PART_SUCCESSOR computes the lexicographic partition successor.
!
!  Discussion:
!
!    PART_SUCCESSOR is "inspired by" the GenPartitions algorithm,
!    but instead of relying on recursion, generates the partitions
!    one at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) NPART, the number of parts of the
!    partition.  1 <= NPART <= N.
!
!    Input/output, integer ( kind = 4 ) A(N), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) asum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) rank
!
!  Return the first element.
!
  if ( rank == -1 ) then
    a(1:n) = 1
    npart = n
    rank = 0
    return
  end if
!
!  Check.
!
  call part_sf_check ( n, npart, a, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PART_SUCCESSOR - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    stop
  end if
!
!  If possible, increment the first intermediate position that
!  is less than its left hand neighbor, and has at least one
!  right hand neighbor.
!
  ihi = npart - 1

  do i = ihi, 2, -1

    if ( a(i) < a(i-1) ) then
      asum = sum ( a(i+1:npart) ) - 1
      a(i) = a(i) + 1
      a(i+1:npart) = 0
      npart = i + asum
      a(i+1:npart) = 1
      rank = rank + 1
      return
    end if

  end do
!
!  A) there are two or more parts
!  Increment the first, replace the rest by 1's.
!
  if ( 2 <= npart ) then
    a(1) = a(1) + 1
    a(2:npart) = 0
    npart = n - a(1) + 1
    a(2:npart) = 1
    rank = rank + 1
!
!  B) there's only one part.
!  We've reached the last item.
!  Return the first one.
!
  else if ( npart == 1 ) then
    a(1:n) = 1
    npart = n
    rank = 0
  end if

  return
end
subroutine part_table ( n, p )

!*****************************************************************************80
!
!! PART_TABLE tabulates the number of partitions of N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) P(0:N), P(I) is the number of partitions of I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(0:n)
  integer ( kind = 4 ) psum
  integer ( kind = 4 ) sign
  integer ( kind = 4 ) w
  integer ( kind = 4 ) wprime

  p(0) = 1
  p(1) = 1

  do i = 2, n

    sign = 1
    psum = 0
    w = 1
    j = 1
    wprime = w + j

    do while ( w < n )

      if ( 0 <= i - w ) then
        if ( sign == 1 ) then
          psum = psum + p(i-w)
        else
          psum = psum - p(i-w)
        end if
      end if

      if ( wprime <= i ) then

        if ( sign == 1 ) then
          psum = psum + p(i-wprime)
        else
          psum = psum - p(i-wprime)
        end if

      end if

      w = w + 3 * j + 1
      j = j + 1
      wprime = w + j
      sign = - sign

    end do

    p(i) = psum

  end do

  return
end
subroutine partition_greedy ( n, a, indx )

!*****************************************************************************80
!
!! PARTITION_GREEDY attacks the partition problem with a greedy algorithm.
!
!  Discussion:
!
!    Given a collection of N not necessarily distinct positive integers A(I),
!    it is desired to partition the values into two groups, whose sums are
!    as close as possible.
!
!  Algorithm:
!
!    Begin with sets 1 and 2 empty.
!
!    Process the data in descending order of magnitude.
!
!    The next item A(I) is added to set 1 or set 2, whichever has the
!    smallest current sum.
!
!    Stop as soon as all items have been allocated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Brian Hayes,
!    The Easiest Hard Problem,
!    American Scientist,
!    Volume 90, Number 2, March-April 2002, pages 113-117.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.  N must be positive.
!
!    Input/output, integer ( kind = 4 ) A(N), a collection of positive values.
!    On output, A has been sorted into descending order.
!
!    Output, integer ( kind = 4 ) INDX(N); INDX(I) is 1 if A(I) is part of
!    set 1, and 2 if it is assigned to set 2.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) sums(2)

  sums(1:2) = 0

  call i4vec_sort_insert_d ( n, a )

  do i = 1, n

    if ( sums(1) < sums(2) ) then
      j = 1
    else
      j = 2
    end if

    indx(i) = j
    sums(j) = sums(j) + a(i)

  end do

  return
end
subroutine partn_enum ( n, nmax, npartitions )

!*****************************************************************************80
!
!! PARTN_ENUM enumerates the partitions of N with maximum element NMAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    Normally N must be positive, but for this routine any
!    N is allowed.
!
!    Input, integer ( kind = 4 ) NMAX, the maximum element in the partition.
!    Normally, 1 <= NMAX <= N is required,
!    but for this routine any value of NMAX is allowed.
!
!    Output, integer ( kind = 4 ) NPARTITIONS is the number of partitions of N
!    with maximum element NMAX.
!
  implicit none

  integer ( kind = 4 ), parameter :: nbig = 25

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nmax
  integer ( kind = 4 ) npartitions
  integer ( kind = 4 ) p(0:nbig,0:nbig)

  if ( n <= 0 ) then

    npartitions = 0

  else if ( nmax <= 0 .or. n < nmax ) then

    npartitions = 0

  else

    call npart_table ( n, nmax, nbig, p )

    npartitions = p(n,nmax)

  end if

  return
end
subroutine partn_sf_check ( n, nmax, npart, a, ierror )

!*****************************************************************************80
!
!! PARTN_SF_CHECK checks an SF partition of an integer with largest entry NMAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NMAX, the value of the largest entry.
!    1 <= NMAX <= N.
!
!    Input, integer ( kind = 4 ) NPART, the number of parts of the partition.
!    1 <= NPART <= N.
!
!    Input, integer ( kind = 4 ) A(NPART), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.  The entries must be in DESCENDING order.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    -1, N is illegal.
!    -2, NMAX is illegal.
!    -3, NPART is illegal.
!    -3, the entries do not add up to N.
!    I, the I-th entry of A is illegal.
!
  implicit none

  integer ( kind = 4 ) npart

  integer ( kind = 4 ) a(npart)
  integer ( kind = 4 ) asum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nmax

  ierror = 0

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARTN_SF_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( nmax < 1 .or. n < nmax ) then
    ierror = -2
    return
  end if

  if ( npart < 1 .or. n < npart ) then
    ierror = -3
    return
  end if
!
!  Entry 1 must be NMAX.
!
  if ( a(1) /= nmax ) then
    ierror = 1
    return
  end if
!
!  Every entry must lie between 1 and N.
!
  do i = 1, npart
    if ( a(i) < 1 .or. n < a(i) ) then
      ierror = i
      return
    end if
  end do
!
!  The entries must be in descending order.
!
  do i = 2, npart
    if ( a(i-1) < a(i) ) then
      ierror = i
      return
    end if
  end do
!
!  The entries must add up to N.
!
  asum = 0
  do i = 1, npart
    asum = asum + a(i)
    if ( n < asum ) then
      ierror = i
      return
    end if
  end do

  if ( asum /= n ) then
    ierror = -3
  end if

  return
end
subroutine partn_successor ( n, nmax, npart, a, rank )

!*****************************************************************************80
!
!! PARTN_SUCCESSOR computes partitions whose largest part is NMAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NMAX, the maximum size of any part of the
!    partition.  1 <= NMAX <= N.
!
!    Input/output, integer ( kind = 4 ) NPART, the number of parts of the
!    partition.  1 <= NPART <= N.
!
!    Input/output, integer ( kind = 4 ) A(N), contains the partition.
!    A(1) through A(NPART) contain the nonzero integers which
!    sum to N.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) index
  integer ( kind = 4 ) nmax
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) temp
!
!  Return the first element.
!
  if ( rank == -1 ) then
    a(1) = nmax
    npart = n + 1 - nmax
    a(2:npart) = 1
    rank = 0
    return
  end if
!
!  Check.
!
  call partn_sf_check ( n, nmax, npart, a, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARTN_SUCCESSOR - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  If there are at least two parts, and the next to last is not NMAX,
!  then rob the last part and pay the next to the last part.
!  Then, if the next to last part is too big, swap it leftwards.
!
  if ( 1 < npart ) then

    if ( a(npart-1) < nmax ) then

      a(npart) = a(npart) - 1
      a(npart-1) = a(npart-1) + 1
      index = npart - 1

      do

        if ( index <= 1 ) then
          exit
        end if

        if ( a(index) <= a(index-1) ) then
          exit
        end if

        temp       = a(index-1)
        a(index-1) = a(index)
        a(index)   = temp

        index = index - 1

      end do
!
!  Sum the tail.
!
      temp = sum ( a(index+1:npart) )
!
!  Spread the sum as 1's.
!
      npart = index + temp
      a(index+1:npart) = 1
      rank = rank + 1
      return

    end if
!
!  Otherwise, we've reached the last item.
!  Return the first one.
!
  else

    npart = n + 1 - nmax
    a(1) = nmax
    a(2:npart) = 1
    rank = 0
    return

  end if

  return
end
subroutine perm_check ( n, p )

!*****************************************************************************80
!
!! PERM_CHECK checks a representation of a permutation.
!
!  Discussion:
!
!    The routine is given N and P, a vector of length N.
!    P is a legal represention of a permutation of the integers from
!    1 to N if and only if every integer from 1 to N occurs
!    as a value of P(I) for some I between 1 and N.
!
!    If an error is observed in the input, this routine prints a message
!    and stops.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values being permuted.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) P(N), the array to check.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifind
  integer ( kind = 4 ) iseek
  integer ( kind = 4 ) p(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i12)' ) '  N = ', n
    stop
  end if

  do i = 1, n
    if ( p(i) < 1 .or. n < p(i) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  P(I) < 1 or N < P(I).'
      write ( *, '(a,i12)' ) '  N = ', n
      write ( *, '(a,i12)' ) '  I = ', i
      write ( *, '(a,i12)' ) '  P(I) = ', p(i)
      stop
    end if
  end do

  do iseek = 1, n

    ifind = -1

    do i = 1, n
      if ( p(i) == iseek ) then
        ifind = i
        exit
      end if
    end do

    if ( ifind == -1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  Every I from 1 to N must occur in P.'
      write ( *, '(a)' ) '  Could not find occurrence of ', iseek
      stop
    end if

  end do

  return
end
subroutine perm_enum ( n, nperm )

!*****************************************************************************80
!
!! PERM_ENUM enumerates the permutations on N digits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values being permuted.
!    N must be nonnegative.
!
!    Output, integer ( kind = 4 ) NPERM, the number of distinct elements.
!
  implicit none

  integer ( kind = 4 ) i4_factorial
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nperm

  nperm = i4_factorial ( n )

  return
end
subroutine perm_inv ( n, p, pinv  )

!*****************************************************************************80
!
!! PERM_INV computes the inverse of a permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values being permuted.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) P(N), describes the permutation.
!    P(I) is the item which is permuted into the I-th place
!    by the permutation.
!
!    Output, integer ( kind = 4 ) PINV(N), the inverse permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) pinv(n)
!
!  Check.
!
  call perm_check ( n, p )

  do i = 1, n
    pinv(p(i)) = i
  end do

  return
end
subroutine perm_lex_rank ( n, p, rank )

!*****************************************************************************80
!
!! PERM_LEX_RANK computes the lexicographic rank of a permutation.
!
!  Discussion:
!
!    The original code altered the input permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values being permuted.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) P(N), describes the permutation.
!    P(I) is the item which is permuted into the I-th place
!    by the permutation.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_factorial
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) pcopy(n)
  integer ( kind = 4 ) rank
!
!  Check.
!
  call perm_check ( n, p )

  rank = 0
  pcopy(1:n) = p(1:n)

  do j = 1, n

    rank = rank + ( pcopy(j) - 1 ) * i4_factorial ( n - j )

    do i = j + 1, n
      if ( pcopy(j) < pcopy(i) ) then
        pcopy(i) = pcopy(i) - 1
      end if
    end do

  end do

  return
end
subroutine perm_lex_successor ( n, p, rank )

!*****************************************************************************80
!
!! PERM_LEX_SUCCESSOR computes the lexicographic permutation successor.
!
!  Example:
!
!    RANK  Permutation
!
!       0  1 2 3 4
!       1  1 2 4 3
!       2  1 3 2 4
!       3  1 3 4 2
!       4  1 4 2 3
!       5  1 4 3 2
!       6  2 1 3 4
!       ...
!      23  4 3 2 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values being permuted.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) P(N), describes the permutation.
!    P(I) is the item which is permuted into the I-th place
!    by the permutation.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) temp
!
!  Return the first element.
!
  if ( rank == -1 ) then
    call i4vec_indicator ( n, p )
    rank = 0
    return
  end if
!
!  Check.
!
  call perm_check ( n, p )
!
!  Seek I, the highest index for which the next element is bigger.
!
  i = n - 1

  do

    if ( i <= 0 ) then
      exit
    end if

    if ( p(i) <= p(i+1) ) then
      exit
    end if

    i = i - 1

  end do
!
!  If no I could be found, then we have reach the final permutation,
!  N, N-1, ..., 2, 1.  Time to start over again.
!
  if ( i == 0 ) then
    call i4vec_indicator ( n, p )
    rank = 0
  else
!
!  Otherwise, look for the the highest index after I whose element
!  is bigger than I's.  We know that I+1 is one such value, so the
!  loop will never fail.
!
    j = n
    do while ( p(j) < p(i) )
      j = j - 1
    end do
!
!  Interchange elements I and J.
!
    temp = p(i)
    p(i) = p(j)
    p(j) = temp
!
!  Reverse the elements from I+1 to N.
!
    call i4vec_reverse ( n - i, p(i+1:n) )
    rank = rank + 1

  end if

  return
end
subroutine perm_lex_unrank ( rank, n, p )

!*****************************************************************************80
!
!! PERM_LEX_UNRANK computes the permutation of given lexicographic rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the permutation.
!
!    Input, integer ( kind = 4 ) N, the number of values being permuted.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) P(N), describes the permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_factorial
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nperm
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call perm_enum ( n, nperm )

  if ( rank < 0 .or. nperm < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  rank_copy = rank

  p(n) = 1

  do j = 1, n - 1

    d = mod ( rank_copy, i4_factorial ( j + 1 ) ) / i4_factorial ( j )
    rank_copy = rank_copy - d * i4_factorial ( j )
    p(n-j) = d + 1

    do i = n - j + 1, n

      if ( d < p(i) ) then
        p(i) = p(i) + 1
      end if

    end do

  end do

  return
end
subroutine perm_mul ( n, p, q, r  )

!*****************************************************************************80
!
!! PERM_MUL computes the product of two permutations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,inson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values being permuted.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) P(N), Q(N), describes the permutation factors.
!
!    Output, integer ( kind = 4 ) R(N), the product permutation P * Q.
!    R(I) = P(Q(I)).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) q(n)
  integer ( kind = 4 ) r(n)
  integer ( kind = 4 ) s(n)
!
!  Check.
!
  call perm_check ( n, p )

  call perm_check ( n, q )
!
!  Use a temporary vector for the result, to avoid problems if
!  some arguments are actually identified.
!
  s(1:n) = p(q(1:n))

  r(1:n) = s(1:n)

  return
end
subroutine perm_parity ( n, p, parity )

!*****************************************************************************80
!
!! PERM_PARITY computes the parity of a permutation.
!
!  Discussion:
!
!    The routine requires the use of a temporary array.
!
!    A permutation is called "even" or "odd", depending on whether
!    it is equivalent to an even or odd number of pairwise
!    transpositions.  This is known as the "parity" of the
!    permutation.
!
!    The "sign" of a permutation is +1 if it has even parity,
!    and -1 if it has odd parity.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values being permuted.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) P(N), describes the permutation.
!    P(I) is the item which is permuted into the I-th place
!    by the permutation.
!
!    Output, integer ( kind = 4 ) PARITY, the parity of the permutation.
!    0, the permutation has even parity.
!    1, the permutation has odd parity.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) parity
  integer ( kind = 4 ) p(n)
!
!  Check.
!
  call perm_check ( n, p )

  a(1:n) = 0

  c = 0

  do j = 1, n

    if ( a(j) == 0 ) then

      c = c + 1
      a(j) = 1
      i = j

      do while ( p(i) /= j )
        i = p(i)
        a(i) = 1
      end do

    end if

  end do

  parity = mod ( n - c, 2 )

  return
end
subroutine perm_print ( n, p, title )

!*****************************************************************************80
!
!! PERM_PRINT prints a permutation.
!
!  Example:
!
!    Input:
!
!      P = 7 2 4 1 5 3 6
!
!    Printed output:
!
!      "This is the permutation:"
!
!      1 2 3 4 5 6 7
!      7 2 4 1 5 3 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation, in standard index form.
!
!    Input, character ( len = * ) TITLE, a title.
!    If no title is supplied, then only the permutation is printed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), parameter :: inc = 20
  integer ( kind = 4 ) p(n)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )

    do ilo = 1, n, inc
      ihi = min ( n, ilo + inc - 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(20i4)' ) ( i, i = ilo, ihi )
      write ( *, '(20i4)' ) p(ilo:ihi)
    end do

  else

    do ilo = 1, n, inc
      ihi = min ( n, ilo + inc - 1 )
      write ( *, '(20i4)' ) p(ilo:ihi)
    end do

  end if

  return
end
subroutine perm_tj_rank ( n, p, rank )

!*****************************************************************************80
!
!! PERM_TJ_RANK computes the Trotter-Johnson rank of a permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values being permuted.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) P(N), describes the permutation.
!    P(I) is the item which is permuted into the I-th place
!    by the permutation.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) rank
!
!  Check.
!
  call perm_check ( n, p )

  rank = 0

  do j = 2, n

    k = 1
    i = 1

    do while ( p(i) /= j )
      if ( p(i) < j ) then
        k = k + 1
      end if
      i = i + 1
    end do

    if ( mod ( rank, 2 ) == 0 ) then
      rank = j * rank + j - k
    else
      rank = j * rank + k - 1
    end if

  end do

  return
end
subroutine perm_tj_successor ( n, p, rank )

!*****************************************************************************80
!
!! PERM_TJ_SUCCESSOR computes the Trotter-Johnson permutation successor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values being permuted.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) P(N), describes the permutation.
!    P(I) is the item which is permuted into the I-th place
!    by the permutation.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) d
  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) par
  integer ( kind = 4 ) q(n)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) st
  integer ( kind = 4 ) temp
!
!  Return the first element.
!
  if ( rank == -1 ) then
    call i4vec_indicator ( n, p )
    rank = 0
    return
  end if
!
!  Check.
!
  call perm_check ( n, p )

  st = 0
  q(1:n) = p(1:n)
  done = .false.
  m = n

  do while ( 1 < m .and. .not. done )

    d = 1
    do while ( q(d) /= m )
      d = d + 1
    end do

    do i = d, m - 1
      q(i) = q(i+1)
    end do

    call perm_parity ( m - 1, q, par )

    if ( par == 1 ) then

      if ( d == m ) then
        m = m - 1
      else
        temp      = p(st+d)
        p(st+d)   = p(st+d+1)
        p(st+d+1) = temp
        done = .true.
      end if

    else

      if ( d == 1 ) then
        m = m - 1
        st = st + 1
      else
        temp      = p(st+d)
        p(st+d)   = p(st+d-1)
        p(st+d-1) = temp
        done = .true.
      end if

    end if

  end do
!
!  Last element was input.  Return first one.
!
  if ( m == 1 ) then
    call i4vec_indicator ( n, p )
    rank = 0
    return
  end if

  rank = rank + 1

  return
end
subroutine perm_tj_unrank ( rank, n, p )

!*****************************************************************************80
!
!! PERM_TJ_UNRANK computes the permutation of given Trotter-Johnson rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the permutation.
!
!    Input, integer ( kind = 4 ) N, the number of values being permuted.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) P(N), describes the permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_factorial
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) nperm
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) r1
  integer ( kind = 4 ) r2
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_TJ_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call perm_enum ( n, nperm )

  if ( rank < 0 .or. nperm < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_TJ_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  p(1) = 1
  r2 = 0

  do j = 2, n
!
!  Replace this ratio of factorials!
!
    r1 = ( rank * i4_factorial ( j ) ) / i4_factorial ( n )
    k = r1 - j * r2

    if ( mod ( r2, 2 ) == 0 ) then
      jhi = j - k
    else
      jhi = k + 1
    end if

    do i = j - 1, jhi, -1
      p(i+1) = p(i)
    end do
    p(jhi) = j

    r2 = r1

  end do

  return
end
subroutine perm_to_cycle ( n, p, ncycle, t, index )

!*****************************************************************************80
!
!! PERM_TO_CYCLE converts a permutation from array to cycle form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values being permuted.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) P(N), describes the permutation using a
!    single array.  For each index I, I -> P(I).
!
!    Output, integer ( kind = 4 ) NCYCLE, the number of cycles.
!    1 <= NCYCLE <= N.
!
!    Output, integer ( kind = 4 ) T(N), INDEX(N), describes the permutation
!    as a collection of NCYCLE cycles.  The first cycle is
!    T(1) -> T(2) -> ... -> T(INDEX(1)) -> T(1).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncycle

  integer ( kind = 4 ) i
  integer ( kind = 4 ) index(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nset
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) t(n)
!
!  Check.
!
  call perm_check ( n, p )
!
!  Initialize.
!
  ncycle = 0
  index(1:n) = 0
  t(1:n) = 0
  nset = 0
!
!  Find the next unused entry.
!
  do i = 1, n

    if ( 0 < p(i) ) then

      ncycle = ncycle + 1
      index(ncycle) = 1

      nset = nset + 1
      t(nset) = p(i)
      p(i) = - p(i)

      do

        j = t(nset)

        if ( p(j) < 0 ) then
          exit
        end if

        index(ncycle) = index(ncycle) + 1

        nset = nset + 1
        t(nset) = p(j)
        p(j) = - p(j)

      end do

    end if

  end do
!
!  If no unused entries remain, we are done.
!  Restore the sign of the permutation and return.
!
  p(1:n) = - p(1:n)

  return
end
subroutine pruefer_check ( n, p, ierror )

!*****************************************************************************80
!
!! PRUEFER_CHECK checks a Pruefer code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the tree.
!    N must be at least 3.
!
!    Input, integer ( kind = 4 ) P(N-2), the Pruefer code for the tree.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    -1, N is less than 3.
!    J, the element P(J) is illegal.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) p(n-2)

  ierror = 0

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRUEFER_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop
  end if

  do i = 1, n - 2
    if ( p(i) < 1 .or. n < p(i) ) then
      ierror = i
      return
    end if
  end do

  return
end
subroutine pruefer_enum ( n, ncode )

!*****************************************************************************80
!
!! PRUEFER_ENUM enumerates the Pruefer codes on N-2 digits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of digits in the code, plus 2.
!    N must be at least 3.
!
!    Output, integer ( kind = 4 ) NCODE, the number of distinct elements.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncode

  ncode = n ** ( n - 2 )

  return
end
subroutine pruefer_rank ( n, p, rank )

!*****************************************************************************80
!
!! PRUEFER_RANK ranks a Pruefer code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the tree.
!    N must be at least 3.
!
!    Input, integer ( kind = 4 ) P(N-2), the Pruefer code for the tree.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the Pruefer code.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(n-2)
  integer ( kind = 4 ) rank
!
!  Check.
!
  call pruefer_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRUEFER_RANK - Fatal error!'
    write ( *, '(a)' ) '  Input array is illegal.'
    write ( *, '(a,i8)' ) '  Error code = ', ierror
    stop
  end if

  rank = 0
  k = 1
  do i = n - 2, 1, -1
    rank = rank + k * ( p(i) - 1 )
    k = k * n
  end do

  return
end
subroutine pruefer_successor ( n, p, rank )

!*****************************************************************************80
!
!! PRUEFER_SUCCESSOR computes the lexical Pruefer sequence successor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the tree.
!    N must be at least 3.
!
!    Input/output, integer ( kind = 4 ) P(N-2), on input, the Pruefer code
!    for a tree, and on output, its lexical successor.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) p(n-2)
!
!  Return the first element.
!
  if ( rank == -1 ) then
    p(1:n-2) = 1
    rank = 0
    return
  end if
!
!  Check.
!
  call pruefer_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRUEFER_SUCCESSOR - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if

  j = n - 2

  do

    if ( p(j) /= n ) then
      exit
    end if

    j = j - 1

    if ( j <= 0 ) then
      exit
    end if

  end do

  if ( j /= 0 ) then
    p(j) = p(j) + 1
    p(j+1:n-2) = 1
    rank = rank + 1
  else
    p(1:n-2) = 1
    rank = 0
  end if

  return
end
subroutine pruefer_to_tree ( n, p, t )

!*****************************************************************************80
!
!! PRUEFER_TO_TREE converts a Pruefer code to a tree.
!
!  Discussion:
!
!    The original code attempts to tack on an extra entry to P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the tree.
!    N must be at least 3.
!
!    Input, integer ( kind = 4 ) P(N-2), the Pruefer code for the tree.
!
!    Output, integer ( kind = 4 ) T(2,N-1), describes the edges of the tree
!    as pairs of nodes.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) p(n-2)
  integer ( kind = 4 ) t(2,n-1)
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
!
!  Check.
!
  call pruefer_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRUEFER_TO_TREE - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal!'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  Initialize the tree to 0.
!
  t(1:2,1:n-1) = 0

  d(1:n) = 1

  do i = 1, n - 2
    d(p(i)) = d(p(i)) + 1
  end do

  do i = 1, n - 1

    x = n
    do while ( d(x) /= 1 )
      x = x - 1
    end do

    if ( i == n - 1 ) then
      y = 1
    else
      y = p(i)
    end if

    d(x) = d(x) - 1
    d(y) = d(y) - 1

    t(1,i) = x
    t(2,i) = y

  end do

  return
end
subroutine pruefer_unrank ( rank, n, p )

!*****************************************************************************80
!
!! PRUEFER_UNRANK unranks a Pruefer code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the Pruefer code.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the tree.
!    N must be at least 3.
!
!    Output, integer ( kind = 4 ) P(N-2), the Pruefer code for the tree.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ncode
  integer ( kind = 4 ) p(n-2)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRUEFER_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call pruefer_enum ( n, ncode )

  if ( rank < 0 .or. ncode < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRUEFER_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  rank_copy = rank

  do i = n - 2, 1, -1
    p(i) = mod ( rank_copy, n ) + 1
    rank_copy = ( rank_copy - p(i) + 1 ) / n
  end do

  return
end
subroutine queens ( n, iarray, k, nstack, istack, maxstack )

!*****************************************************************************80
!
!! QUEENS finds possible positions for the K-th nonattacking queen.
!
!  Discussion:
!
!    The chessboard is N by N, and is being filled one column at a time,
!    with a tentative solution to the nonattacking queen problem.  So
!    far, K-1 rows have been chosen, and we now need to provide a list
!    of all possible rows that might be used in column K.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the total number of queens to place, and
!    the length of a side of the chessboard.
!
!    Input, integer ( kind = 4 ) IARRAY(N).  The first K-1 entries of IARRAY
!    record the rows into which queens have already been placed.
!
!    Input, integer ( kind = 4 ) K, the column for which we need possible
!    row positions for the next queen.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of stack.
!    On output, this has been updated.
!
!    Input/output, integer ( kind = 4 ) ISTACK(MAXSTACK).  On output, we
!    have added the candidates, and the number of candidates, to the end
!    of the stack.
!
!    Input, integer ( kind = 4 ) MAXSTACK, maximum dimension of ISTACK.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) maxstack

  logical diag
  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) istack(maxstack)
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncan
  integer ( kind = 4 ) nstack
  logical row

  ncan = 0

  do irow = 1, n
!
!  If row IROW has already been used, that's it.
!
    row = .false.

    do jcol = 1, k - 1
      if ( iarray(jcol) == irow ) then
        row = .true.
      end if
    end do

    if ( .not. row ) then

      diag = .false.

      do jcol = 1, k - 1

        if ( irow == iarray(jcol) + k - jcol .or. &
             irow == iarray(jcol) - ( k - jcol ) ) then

          diag = .true.

        end if

      end do

      if ( .not. diag ) then
        ncan = ncan + 1
        nstack = nstack + 1
        istack(nstack) = irow
      end if

    end if

  end do

  nstack = nstack + 1
  istack(nstack) = ncan

  return
end
function r4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R4_UNIFORM returns a scaled pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) ( kind = 4 ) SEED, the "seed" value,
!    which should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ) k
  real ( kind = 4 ) r4_uniform
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r4_uniform = a + ( b - a ) * real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
function r8_gamma_log ( x )

!*****************************************************************************80
!
!! R8_GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!  Discussion:
!
!    The program uses rational functions that theoretically approximate
!    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
!    approximation for 12 < X is from reference Hart et al, while approximations
!    for X < 12.0D+00 are similar to those in reference Cody and Hillstrom,
!    but are unpublished.  The accuracy achieved depends on the arithmetic
!    system, the compiler, intrinsic functions, and proper selection of the
!    machine-dependent constants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the
!    Gamma Function,
!    Mathematics of Computation,
!    Volume 21, Number 98, April 1967, pages 198-203.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thacher,
!    Christoph Witzgall,<br>
!    Computer Approximations,<br>
!    Wiley, 1968,<br>
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function.
!    X must be positive.
!
!    Output, real ( kind = 8 ) R8_GAMMA_LOG, the logarithm of the Gamma
!    function of X.
!    If X <= 0.0, or if overflow would occur, the program returns the
!    largest representable floating point number.
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
!  FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!  Approximate values for some important machines are:
!
!                            BETA      MAXEXP         XBIG
!
!  CRAY-1        (S.P.)        2        8191       9.62D+2461
!  Cyber 180/855
!    under NOS   (S.P.)        2        1070       1.72D+319
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)        2         128       4.08D+36
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!  IBM 3033      (D.P.)       16          63       4.29D+73
!  VAX D-Format  (D.P.)        2         127       2.05D+36
!  VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                           FRTBIG
!
!  CRAY-1        (S.P.)   3.13D+615
!  Cyber 180/855
!    under NOS   (S.P.)   6.44D+79
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   1.42D+9
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)   2.25D+76
!  IBM 3033      (D.P.)   2.56D+18
!  VAX D-Format  (D.P.)   1.20D+9
!  VAX G-Format  (D.P.)   1.89D+76
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = 8 ) corr
  real ( kind = 8 ), parameter :: d1 = - 5.772156649015328605195174D-01
  real ( kind = 8 ), parameter :: d2 =   4.227843350984671393993777D-01
  real ( kind = 8 ), parameter :: d4 =   1.791759469228055000094023D+00
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: frtbig = 1.42D+09
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real ( kind = 8 ), parameter :: pnt68 = 0.6796875D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 4.08D+36
  real ( kind = 8 ) xden
  real ( kind = 8 ) xm1
  real ( kind = 8 ) xm2
  real ( kind = 8 ) xm4
  real ( kind = 8 ) xnum
  real ( kind = 8 ) xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0D+00 .or. xbig < x ) then
    r8_gamma_log = huge ( r8_gamma_log )
    return
  end if

  if ( x <= epsilon ( x ) ) then

    res = - log ( x )

  else if ( x <= 1.5D+00 ) then

    if ( x < pnt68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0D+00
      xm1 = ( x - 0.5D+00 ) - 0.5D+00
    end if

    if ( x <= 0.5D+00 .or. pnt68 <= x ) then

      xden = 1.0D+00
      xnum = 0.0D+00

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5D+00 ) - 0.5D+00
      xden = 1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0D+00 ) then

    xm2 = x - 2.0D+00
    xden = 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0D+00 ) then

    xm4 = x - 4.0D+00
    xden = - 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0D+00

    if ( x <= frtbig ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + sqrtpi - 0.5D+00 * corr
    res = res + x * ( corr - 1.0D+00 )

  end if

  r8_gamma_log = res

  return
end
subroutine r8vec_backtrack ( n, x, indx, k, nstack, stack, maxstack, ncan )

!*****************************************************************************80
!
!! R8VEC_BACKTRACK supervises a backtrack search for an R8VEC.
!
!  Discussion:
!
!    The routine tries to construct a real vector one index at a time,
!    using possible candidates as supplied by the user.
!
!    At any time, the partially constructed vector may be discovered to be
!    unsatisfactory, but the routine records information about where the
!    last arbitrary choice was made, so that the search can be
!    carried out efficiently, rather than starting out all over again.
!
!    First, call the routine with INDX = 0 so it can initialize itself.
!
!    Now, on each return from the routine, if INDX is:
!      1, you've just been handed a complete candidate vector;
!         Admire it, analyze it, do what you like.
!      2, please determine suitable candidates for position X(K).
!         Return the number of candidates in NCAN(K), adding each
!         candidate to the end of STACK, and increasing NSTACK.
!      3, you're done.  Stop calling the routine;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of positions to be filled in
!    the vector.
!
!    Input/output, real ( kind = 8 ) X(N), the partial or complete
!    candidate vector.
!
!    Input/output, integer ( kind = 4 ) INDX, a communication flag.
!    On input,
!      0 to start a search.
!    On output:
!      1, a complete output vector has been determined and returned in X(1:N);
!      2, candidates are needed for position X(K);
!      3, no more possible vectors exist.
!
!    Input/output, integer ( kind = 4 ) K, if INDX=2, the current vector index
!    being considered.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!
!    Input/output, real ( kind = 8 ) STACK(MAXSTACK), a list of all current
!    candidates for all positions 1 through K.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum length of the stack.
!
!    Input/output, integer ( kind = 4 ) NCAN(N), lists the current number
!    of candidates for positions 1 through K.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) indx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncan(n)
  integer ( kind = 4 ) nstack
  real ( kind = 8 ) stack(maxstack)
  real ( kind = 8 ) x(n)
!
!  If this is the first call, request a candidate for position 1.
!
  if ( indx == 0 ) then
    k = 1
    nstack = 0
    indx = 2
    return
  end if
!
!  Examine the stack.
!
  do
!
!  If there are candidates for position K, take the first available
!  one off the stack, and increment K.
!
!  This may cause K to reach the desired value of N, in which case
!  we need to signal the user that a complete set of candidates
!  is being returned.
!
    if ( 0 < ncan(k) ) then

      x(k) = stack(nstack)
      nstack = nstack - 1

      ncan(k) = ncan(k) - 1

      if ( k /= n ) then
        k = k + 1
        indx = 2
      else
        indx = 1
      end if

      exit
!
!  If there are no candidates for position K, then decrement K.
!  If K is still positive, repeat the examination of the stack.
!
    else

      k = k - 1

      if ( k <= 0 ) then
        indx = 3
        exit
      end if

    end if

  end do

  return
end
subroutine rgf_check ( m, f, ierror )

!*****************************************************************************80
!
!! RGF_CHECK checks a restricted growth function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the domain of the RGF is the integers
!    from 1 to M.  M must be positive.
!
!    Input, integer ( kind = 4 ) F(M), the restricted growth function.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    -1, M is illegal.
!    I, entry I of the restricted growth function is illegal.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) f(m)
  integer ( kind = 4 ) fmax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror

  ierror = 0

  if ( m < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RGF_CHECK - Fatal error!'
    write ( *, '(a)' ) '  M < 1.'
    stop
  end if

  fmax = 0
  do i = 1, m
    if ( f(i) <= 0 .or. fmax + 1 < f(i) ) then
      ierror = i
      return
    end if
    fmax = max ( fmax, f(i) )
  end do

  return
end
subroutine rgf_enum ( m, nrgf )

!*****************************************************************************80
!
!! RGF_ENUM enumerates the restricted growth functions on M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the domain of the RGF is the integers
!    from 1 to M.  M must be positive.  However, for the enumeration routine
!    only, it is legal to call with any value of M.
!
!    Output, integer ( kind = 4 ) NRGF, the number of restricted growth
!    functions.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) b(0:m)
  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nrgf

  if ( m < 0 ) then

    nrgf = 0

  else if ( m == 0 ) then

    nrgf = 1

  else

    b(0) = 1
    do j = 1, m
      b(j) = 0
      do i = 0, j - 1
        b(j) = b(j) + binomial ( j - 1, i ) * b(i)
      end do
    end do

    nrgf = b(m)

  end if

  return
end
subroutine rgf_g_table ( m, mmax, d )

!*****************************************************************************80
!
!! RGF_G_TABLE tabulates the generalized restricted growth functions.
!
!  Example:
!
!    M = 6
!
!    D =  1    1    1    1    1    1    1
!         1    2    3    4    5    6    0
!         2    5   10   17   26    0    0
!         5   15   37   77    0    0    0
!        15   52  151    0    0    0    0
!        52  203    0    0    0    0    0
!       203    0    0    0    0    0    0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, indicates how many rows and columns are to
!    be computed.  M must be nonnegative.
!
!    Input, integer ( kind = 4 ) MMAX, the value used to allocate space for the
!    D array.  MMAX must be at least M.
!
!    Output, integer ( kind = 4 ) D(0:MMAX,0:MMAX), the first M+1 rows and
!    M+1 columns of the table of the number of generalized restricted growth
!    functions.  D(I,J) is the number of GRGF's of length I with restriction
!    parameter J.
!
  implicit none

  integer ( kind = 4 ) mmax

  integer ( kind = 4 ) d(0:mmax,0:mmax)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m

  d(0,0:m) = 1

  do i = 1, m
    do j = 0, m
      if ( j <= m - i ) then
        d(i,j) = j * d(i-1,j) + d(i-1,j+1)
      else
        d(i,j) = 0
      end if
    end do
  end do

  return
end
subroutine rgf_rank ( m, f, rank )

!*****************************************************************************80
!
!! RGF_RANK ranks a restricted growth function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the domain of the RGF is the integers
!    from 1 to M.  M must be positive.
!
!    Input, integer ( kind = 4 ) F(M), the restricted growth function.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the restricted growth
!    function.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) d(0:m,0:m)
  integer ( kind = 4 ) f(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) rank
!
!  Check.
!
  call rgf_check ( m, f, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RGF_RANK - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal!'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  Get the generalized restricted growth function table.
!
  call rgf_g_table ( m, m, d )

  rank = 0
  j = 1
  do i = 2, m
    rank = rank + ( f(i) - 1 ) * d(m-i,j)
    j = max ( j, f(i) )
  end do

  return
end
subroutine rgf_successor ( m, f, rank )

!*****************************************************************************80
!
!! RGF_SUCCESSOR generates the next restricted growth function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the domain of the RGF is the integers
!    from 1 to M.  M must be positive.
!
!    Input/output, integer ( kind = 4 ) F(M), the restricted growth function.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) f(m)
  integer ( kind = 4 ) fmax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) rank
!
!  Return the first element.
!
  if ( rank == -1 ) then
    f(1:m) = 1
    rank = 0
    return
  end if
!
!  Check.
!
  call rgf_check ( m, f, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RGF_SUCCESSOR - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal!'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  Find the first position from the right which can be incremented.
!
  do i = m, 2, -1

    fmax = 1
    do j = 2, i - 1
      fmax = max ( fmax, f(j) )
    end do
!
!  Increment the function at this position, and set later entries to 1.
!
    if ( f(i) /= fmax + 1 ) then
      f(i) = f(i) + 1
      f(i+1:m) = 1
      rank = rank + 1
      return
    end if

  end do
!
!  The final element was input.
!  Return the first element.
!
  f(1:m) = 1
  rank = 0

  return
end
subroutine rgf_to_setpart ( m, f, nsub, s, index )

!*****************************************************************************80
!
!! RGF_TO_SETPART converts a restricted growth function to a set partition.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the domain of the RGF is the integers
!    from 1 to M.  M must be positive.
!
!    Input, integer ( kind = 4 ) F(M), the restricted growth function.
!
!    Output, integer ( kind = 4 ) NSUB, the number of nonempty subsets into
!    which the set is partitioned.
!
!    Output, integer ( kind = 4 ) S(M), describes the partition of a set of
!    M objects into NSUB nonempty subsets.  If element I of the
!    superset belongs to subset J, then S(I) = J.
!
!    Output, integer ( kind = 4 ) INDEX(M), lists the location in S of the last
!    element of each subset.  Thus, the elements of subset 1
!    are S(1) through S(INDEX(1)), the elements of subset 2
!    are S(INDEX(1)+1) through S(INDEX(2)) and so on.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nsub

  integer ( kind = 4 ) f(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) index(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) s(m)
!
!  Check.
!
  call rgf_check ( m, f, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RGF_TO_SETPART - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal!'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  Determine the number of subsets.
!
  nsub = maxval ( f(1:m) )
!
!  Initialize.
!
  s(1:m) = 0
  index(1:nsub) = 0
!
!  For each subset I, collect the indices of F which have value I.
!  These are the elements of the I-th subset.
!
  k = 0
  do i = 1, nsub
    do j = 1, m
      if ( f(j) == i ) then
        k = k + 1
        s(k) = j
      end if
    end do
    index(i) = k
  end do

  return
end
subroutine rgf_unrank ( rank, m, f )

!*****************************************************************************80
!
!! RGF_UNRANK returns the restricted growth function of a given rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the restricted growth
!    function.
!
!    Input, integer ( kind = 4 ) M, the domain of the RGF is the integers
!    from 1 to M.  M must be positive.
!
!    Output, integer ( kind = 4 ) F(M), the restricted growth function.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) d(0:m,0:m)
  integer ( kind = 4 ) f(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nrgf
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
!
!  Check.
!
  if ( m < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RGF_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input M is illegal.'
    stop
  end if

  call rgf_enum ( m, nrgf )

  if ( rank < 0 .or. nrgf < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RGF_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if
!
!  Get the generalized restricted growth function table.
!
  call rgf_g_table ( m, m, d )

  rank_copy = rank
  j = 1
  f(1) = 1

  do i = 2, m

    if ( j * d(m-i,j) <= rank_copy ) then
      f(i) = j + 1
      rank_copy = rank_copy - j * d(m-i,j)
      j = j + 1
    else
      f(i) = 1 + ( rank_copy / d(m-i,j) )
      rank_copy = mod ( rank_copy, d(m-i,j) )
    end if

  end do

  return
end
subroutine setpart_check ( m, nsub, s, index, ierror )

!*****************************************************************************80
!
!! SETPART_CHECK checks a set partition.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of elements of the set.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) NSUB, the number of nonempty subsets into
!    which the set is partitioned.  1 <= NSUB <= M.
!
!    Input, integer ( kind = 4 ) INDEX(NSUB), lists the location in S of the
!    last element of each subset.  Thus, the elements of subset 1
!    are S(1) through S(INDEX(1)), the elements of subset 2
!    are S(INDEX(1)+1) through S(INDEX(2)) and so on.
!
!    Input, integer ( kind = 4 ) S(M), contains the integers from 1 to M,
!    grouped into subsets as described by INDEX.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    -I, the I-th element of INDEX is illegal.
!    +I, the I-th element of S is illegal.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nsub

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) index(nsub)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s(m)

  ierror = 0
!
!  Check M.
!
  if ( m < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETPART_CHECK - Fatal error!'
    write ( *, '(a)' ) '  M < 1.'
    stop
  end if
!
!  Check INDEX.
!
  imin = 0
  do i = 1, nsub
    if ( index(i) <= imin .or. m < index(i) ) then
      ierror = -i
      return
    end if
    imin = index(i)
  end do
!
!  Check the elements of S.
!
  do i = 1, nsub

    if ( s(i) <= 0 .or. m < s(i) ) then
      ierror = i
      return
    end if

    do j = 1, i - 1
      if ( s(j) == s(i) ) then
        ierror = i
        return
      end if
    end do

  end do

  return
end
subroutine setpart_enum ( m, npart )

!*****************************************************************************80
!
!! SETPART_ENUM enumerates the partitions of a set of M elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of elements in the set.
!    M must be positive.  However, for the enumeration routine only,
!    it is legal to call with any value of M.
!
!    Output, integer ( kind = 4 ) NPART, the number of partitions of the set.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) b(0:m)
  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) npart

  if ( m < 0 ) then

    npart = 0

  else if ( m == 0 ) then

    npart = 1

  else

    b(0) = 1
    do j = 1, m
      b(j) = 0
      do i = 0, j - 1
        b(j) = b(j) + binomial ( j - 1, i ) * b(i)
      end do
    end do

    npart = b(m)

  end if

  return
end
subroutine setpart_to_rgf ( m, nsub, s, index, f )

!*****************************************************************************80
!
!! SETPART_TO_RGF converts a set partition to a restricted growth function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of elements of the set.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) NSUB, the number of nonempty subsets into
!    which the set is partitioned.  1 <= NSUB <= M.
!
!    Input, integer ( kind = 4 ) INDEX(NSUB), lists the location in S of the
!    last element of each subset.  Thus, the elements of subset 1
!    are S(1) through S(INDEX(1)), the elements of subset 2
!    are S(INDEX(1)+1) through S(INDEX(2)) and so on.
!
!    Input, integer ( kind = 4 ) S(M), contains the integers from 1 to M,
!    grouped into subsets as described by INDEX.
!
!    Output, integer ( kind = 4 ) F(M), the restricted growth function from
!    M to NSUB.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nsub

  integer ( kind = 4 ) f(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) index(nsub)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) s(m)
!
!  Check.
!
  call setpart_check ( m, nsub, s, index, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETPART_TO_RGF - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    stop
  end if

  khi = 0
  do i = 1, nsub
    klo = khi + 1
    khi = index(i)
    do k = klo, khi
      f(s(k)) = i
    end do
  end do

  return
end
subroutine stirling_numbers1 ( m, n, s )

!*****************************************************************************80
!
!! STIRLING_NUMBERS1 computes Stirling numbers of the first kind.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the maximum row to compute.
!    M must be nonnegative.
!
!    Input, integer ( kind = 4 ) N, the maximum column to compute.
!    N must be nonnegative.
!
!    Output, integer ( kind = 4 ) S(0:M,0:N), the first M+1 rows and N+1 columns
!    of the table of Stirling numbers of the first kind.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s(0:m,0:n)

  s(0:m,0:n) = 0

  s(0,0) = 1

  do i = 1, m
    do j = 1, n
      if ( j <= i ) then
        s(i,j) = s(i-1,j-1) - ( i - 1 ) * s(i-1,j)
      end if
    end do
  end do

  return
end
subroutine stirling_numbers2 ( m, n, s )

!*****************************************************************************80
!
!! STIRLING_NUMBERS2 computes Stirling numbers of the second kind.
!
!  Discussion:
!
!    The reference has a typographical error, referring to
!    S(I-J,J-1) instead of S(I-1,J-1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the maximum row to compute.
!    M must be nonnegative.
!
!    Input, integer ( kind = 4 ) N, the maximum column to compute.
!    N must be nonnegative.
!
!    Output, integer ( kind = 4 ) S(0:M,0:N), the first M+1 rows and N+1 columns
!    of the table of Stirling numbers of the second kind.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s(0:m,0:n)

  s(0:m,0:n) = 0

  s(0,0) = 1

  do i = 1, m
    do j = 1, n
      if ( j <= i ) then
        s(i,j) = j * s(i-1,j) + s(i-1,j-1)
      end if
    end do
  end do

  return
end
subroutine subset_check ( n, t )

!*****************************************************************************80
!
!! SUBSET_CHECK checks a subset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(N), the subset.  If T(I) = 0, item I is
!    not in the subset; if T(I) = 1, item I is in the subset.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) t(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUBSET_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  do i = 1, n

    if ( t(i) /= 0 .and. t(i) /= 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SUBSET_CHECK - Fatal error!'
      write ( *, '(a)' ) '  T(I) is not 0 and not 1.'
      stop
    end if

  end do

  return
end
subroutine subset_colex_rank ( n, t, rank )

!*****************************************************************************80
!
!! SUBSET_COLEX_RANK computes the colexicographic rank of a subset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the master set.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(N), the subset.  If T(I) = 0, item I is
!    not in the subset; if T(I) = 1, item I is in the subset.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the subset.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(n)
!
!  Check.
!
  call subset_check ( n, t )

  rank = 0

  do i = 1, n

    if ( t(i) == 1 ) then
      rank = rank + 2 ** ( i - 1 )
    end if

  end do

  return
end
subroutine subset_colex_successor ( n, t, rank )

!*****************************************************************************80
!
!! SUBSET_COLEX_SUCCESSOR computes the subset colexicographic successor.
!
!  Discussion:
!
!    In the original code, there is a last element with no successor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) T(N), describes a subset.  T(I) is 0 if
!    the I-th element of the master set is not in the subset, and is
!    1 if the I-th element is part of the subset.
!    On input, T describes a subset.
!    On output, T describes the next subset in the ordering.
!    If the input T was the last in the ordering, then the output T
!    will be the first.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(n)
!
!  Return the first element.
!
  if ( rank == -1 ) then
    t(1:n) = 0
    rank = 0
    return
  end if
!
!  Check.
!
  call subset_check ( n, t )

  do i = 1, n

    if ( t(i) == 0 ) then
      t(i) = 1
      rank = rank + 1
      return
    else
      t(i) = 0
    end if

  end do

  rank = 0

  return
end
subroutine subset_colex_unrank ( rank, n, t )

!*****************************************************************************80
!
!! SUBSET_COLEX_UNRANK computes the subset of given colexicographic rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the subset.
!
!    Input, integer ( kind = 4 ) N, the number of items in the master set.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) T(N), the subsetof the given rank.
!    If T(I) = 0, item I is not in the subset; if T(I) = 1, item I is
!    in the subset.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
  integer ( kind = 4 ) t(n)
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUBSET_COLEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call subset_enum ( n, nsub )

  if ( rank < 0 .or. nsub < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUBSET_COLEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  rank_copy = rank

  do i = 1, n
    if ( mod ( rank_copy, 2 ) == 1 ) then
      t(i) = 1
    else
      t(i) = 0
    end if

    rank_copy = rank_copy / 2

  end do

  return
end
subroutine subset_complement ( n, a, b )

!*****************************************************************************80
!
!! SUBSET_COMPLEMENT computes the complement of a set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set, of which A is
!    a subset.  N must be positive.
!
!    Input, integer ( kind = 4 ) A(N), a subset of the master set.
!    A(I) = 0 if the I-th element is in the subset A, and is
!    1 otherwise.
!
!    Output, integer ( kind = 4 ) B(N), the complement of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)
!
!  Check.
!
  call subset_check ( n, a )

  b(1:n) = 1 - a(1:n)

  return
end
subroutine subset_distance ( n, t1, t2, dist )

!*****************************************************************************80
!
!! SUBSET_DISTANCE computes the Hamming distance between two sets.
!
!  Discussion:
!
!    The sets T1 and T2 are assumed to be subsets of a set of N elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set, of which T1 and
!    T2 are subsets.  N must be positive.
!
!    Input, integer ( kind = 4 ) T1(N), T2(N), two subsets of the master set.
!    T1(I) = 0 if the I-th element is in the subset T1, and is
!    1 otherwise; T2 is defined similarly.
!
!    Output, integer ( kind = 4 ) DIST, the Hamming distance between T1 and T2,
!    defined as the number of elements of the master set which are
!    in either T1 or T2 but not both.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) dist
  integer ( kind = 4 ) i
  integer ( kind = 4 ) t1(n)
  integer ( kind = 4 ) t2(n)
!
!  Check.
!
  call subset_check ( n, t1 )

  call subset_check ( n, t2 )

  dist = 0

  do i = 1, n

    if ( ( t1(i) == 0 .and. t2(i) /= 0 ) .or. &
         ( t1(i) /= 0 .and. t2(i) == 0 ) ) then
      dist = dist + 1
    end if

  end do

  return
end
subroutine subset_enum ( n, nsub )

!*****************************************************************************80
!
!! SUBSET_ENUM enumerates the subsets of a set with N elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the set.
!    N must be at least 0.
!
!    Output, integer ( kind = 4 ) NSUB, the number of distinct elements.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nsub

  nsub = 2**n

  return
end
subroutine subset_intersect ( n, a, b, c )

!*****************************************************************************80
!
!! SUBSET_INTERSECT computes the intersection of two sets.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set, of which A and
!    B are subsets.  N must be positive.
!
!    Input, integer ( kind = 4 ) A(N), B(N), two subsets of the master set.
!    A(I) = 0 if the I-th element is in the subset A, and is
!    1 otherwise; B is defined similarly.
!
!    Output, integer ( kind = 4 ) C(N), the intersection of A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
!
!  Check.
!
  call subset_check ( n, a )

  call subset_check ( n, b )

  do i = 1, n
    c(i) = min ( a(i), b(i) )
  end do

  return
end
subroutine subset_lex_rank ( n, t, rank )

!*****************************************************************************80
!
!! SUBSET_LEX_RANK computes the lexicographic rank of a subset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the master set.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(N), the subset.  If T(I) = 0, item I is
!    not in the subset; if T(I) = 1, item I is in the subset.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the subset.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(n)
!
!  Check.
!
  call subset_check ( n, t )

  rank = 0

  do i = 1, n

    if ( t(i) == 1 ) then
      rank = rank + 2**( n - i )
    end if

  end do

  return
end
subroutine subset_lex_successor ( n, t, rank )

!*****************************************************************************80
!
!! SUBSET_LEX_SUCCESSOR computes the subset lexicographic successor.
!
!  Discussion:
!
!    In the original code, there is a last element with no successor.
!
!  Example:
!
!    On initial call, N is 5 and the input value of RANK is -1.
!    Then here are the successive outputs from the program:
!
!   Rank   T1   T2   T3   T4   T5
!   ----   --   --   --   --   --
!      0    0    0    0    0    0
!      1    0    0    0    0    1
!      2    0    0    0    1    0
!      3    0    0    0    1    1
!     ..   ..   ..   ..   ..   ..
!     30    1    1    1    1    0
!     31    1    1    1    1    1
!      0    0    0    0    0    0  <-- Cycle restarts with first element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) T(N), describes a subset.  T(I) is 0 if
!    the I-th element of the master set is not in the subset, and is
!    1 if the I-th element is part of the subset.
!    On input, T describes a subset.
!    On output, T describes the next subset in the ordering.
!    If the input T was the last in the ordering, then the output T
!    will be the first.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(n)
!
!  Return the first element.
!
  if ( rank == -1 ) then
    t(1:n) = 0
    rank = 0
    return
  end if
!
!  Check.
!
  call subset_check ( n, t )

  do i = n, 1, -1

    if ( t(i) == 0 ) then
      t(i) = 1
      rank = rank + 1
      return
    else
      t(i) = 0
    end if

  end do

  rank = 0

  return
end
subroutine subset_lex_unrank ( rank, n, t )

!*****************************************************************************80
!
!! SUBSET_LEX_UNRANK computes the subset of given lexicographic rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the subset.
!
!    Input, integer ( kind = 4 ) N, the number of items in the master set.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) T(N), the subset of the given rank.
!    If T(I) = 0, item I is not in the subset; if T(I) = 1, item I is in
!    the subset.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
  integer ( kind = 4 ) t(n)
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUBSET_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call subset_enum ( n, nsub )

  if ( rank < 0 .or. nsub < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUBSET_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  rank_copy = rank

  do i = n, 1, -1

    if ( mod ( rank_copy, 2 ) == 1 ) then
      t(i) = 1
    else
      t(i) = 0
    end if

    rank_copy = rank_copy / 2

  end do

  return
end
subroutine subset_union ( n, a, b, c )

!*****************************************************************************80
!
!! SUBSET_UNION computes the union of two sets.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set, of which A and
!    B are subsets.  N must be positive.
!
!    Input, integer ( kind = 4 ) A(N), B(N), two subsets of the master set.
!    A(I) = 0 if the I-th element is in the subset A, and is
!    1 otherwise; B is defined similarly.
!
!    Output, integer ( kind = 4 ) C(N), the union of A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
!
!  Check.
!
  call subset_check ( n, a )

  call subset_check ( n, b )

  do i = 1, n
    c(i) = max ( a(i), b(i) )
  end do

  return
end
subroutine subset_weight ( n, t, weight )

!*****************************************************************************80
!
!! SUBSET_WEIGHT computes the Hamming weight of a set.
!
!  Discussion:
!
!    The Hamming weight is simply the number of elements in the set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set, of which T
!    is a subset.  N must be positive.
!
!    Input, integer ( kind = 4 ) T(N), defines the subset T.
!    T(I) is 1 if I is an element of T, and 0 otherwise.
!
!    Output, integer ( kind = 4 ) WEIGHT, the Hamming weight of the subset T.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) t(n)
  integer ( kind = 4 ) weight
!
!  Check.
!
  call subset_check ( n, t )

  weight = sum ( t(1:n) )

  return
end
subroutine subset_xor ( n, a, b, c )

!*****************************************************************************80
!
!! SUBSET_XOR computes the symmetric difference of two sets.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set, of which A and
!    B are subsets.  N must be positive.
!
!    Input, integer ( kind = 4 ) A(N), B(N), two subsets of the master set.
!    A(I) = 0 if the I-th element is in the subset A, and is
!    1 otherwise; B is defined similarly.
!
!    Output, integer ( kind = 4 ) C(N), the symmetric difference of A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
!
!  Check.
!
  call subset_check ( n, a )

  call subset_check ( n, b )

  do i = 1, n
    c(i) = max ( a(i), b(i) ) - min ( a(i), b(i) )
  end do

  return
end
subroutine subsetsum_swap ( n, a, sum_desired, index, sum_achieved )

!*****************************************************************************80
!
!! SUBSETSUM_SWAP seeks a solution of the subset sum problem by swapping.
!
!  Discussion:
!
!    Given a collection of N not necessarily distinct positive integers A(I),
!    and a positive integer SUM_DESIRED, select a subset of the values so that
!    their sum is as close as possible to SUM_DESIRED without exceeding it.
!
!  Algorithm:
!
!    Start with no values selected, and SUM_ACHIEVED = 0.
!
!    Consider each element A(I):
!
!      If A(I) is not selected and SUM_ACHIEVED + A(I) <= SUM_DESIRED,
!        select A(I).
!
!      If A(I) is still not selected, and there is a selected A(J)
!      such that SUM_GOT < SUM_ACHIEVED + A(I) - A(J),
!        select A(I) and deselect A(J).
!
!      If no items were selected on this sweep,
!        exit.
!      Otherwise,
!        repeat the search.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.  N must be positive.
!
!    Input/output, integer ( kind = 4 ) A(N), a collection of positive values.
!    On output, A has been sorted into descending order.
!
!    Input, integer ( kind = 4 ) SUM_DESIRED, the desired sum.
!
!    Output, integer ( kind = 4 ) INDEX(N); INDEX(I) is 1 if A(I) is part of the
!    sum, and 0 otherwise.
!
!    Output, integer ( kind = 4 ) SUM_ACHIEVED, the sum of the selected
!    elements.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) index(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nmove
  integer ( kind = 4 ) sum_achieved
  integer ( kind = 4 ) sum_desired
!
!  Initialize.
!
  sum_achieved = 0
  index(1:n) = 0
!
!  Sort into descending order.
!
  call i4vec_sort_insert_d ( n, a )

  do

    nmove = 0

    do i = 1, n

      if ( index(i) == 0 ) then

        if ( sum_achieved + a(i) <= sum_desired ) then
          index(i) = 1
          sum_achieved = sum_achieved + a(i)
          nmove = nmove + 1
          cycle
        end if

      end if

      if ( index(i) == 0 ) then

        do j = 1, n

          if ( index(j) == 1 ) then

            if ( sum_achieved < sum_achieved + a(i) - a(j) .and. &
              sum_achieved + a(i) - a(j) <= sum_desired ) then
              index(j) = 0
              index(i) = 1
              nmove = nmove + 2
              sum_achieved = sum_achieved + a(i) - a(j)
              exit
            end if

          end if

        end do

      end if

    end do

    if ( nmove <= 0 ) then
      exit
    end if

  end do

  return
end
subroutine tableau_check ( n, tab )

!*****************************************************************************80
!
!! TABLEAU_CHECK checks a 2 by N tableau.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of columns in the tableau.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) TAB(2,N), a 2 by N tableau.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) tab(2,n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TABLEAU_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
    return
  end if
!
!  The entries must be between 0 and 2*N.
!
  do i = 1, 2
    do j = 1, n
      if ( tab(i,j) < 1 .or. 2 * n < tab(i,j) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TABLEAU_CHECK - Fatal error!'
        write ( *, '(a)' ) '  TAB(I,J) < 1 or N < TAB(I,J).'
        stop
      end if
    end do
  end do
!
!  The entries must be increasing to the right.
!
  do i = 1, 2
    do j = 2, n
      if ( tab(i,j) <= tab(i,j-1) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TABLEAU_CHECK - Fatal error!'
        write ( *, '(a)' ) '  TAB(I,J) <= TAB(I,J-1).'
        stop
      end if
    end do
  end do
!
!  The entries must be increasing down.
!
  i = 2
  do j = 1, n
    if ( tab(i,j) <= tab(i-1,j) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TABLEAU_CHECK - Fatal error!'
      write ( *, '(a)' ) '  TAB(I,J) <= TAB(I-1,J).'
      stop
    end if
  end do

  return
end
subroutine tableau_enum ( n, ntab )

!*****************************************************************************80
!
!! TABLEAU_ENUM enumerates the 2 by N standard tableaus.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of columns in the tableau.
!    N must be nonnegative.
!
!    Output, integer ( kind = 4 ) NTAB, the number of 2 by N standard tableaus.
!
  implicit none

  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ntab

  ntab = binomial ( 2 * n, n ) / ( n + 1 )

  return
end
subroutine tableau_to_bal_seq ( n, tab, t )

!*****************************************************************************80
!
!! TABLEAU_TO_BAL_SEQ converts a 2 by N tableau to a balanced sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of 0's (and 1's) in the sequence.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) TAB(2,N), a 2 by N tableau.
!
!    Output, integer ( kind = 4 ) T(2*N), a balanced sequence.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) t(2*n)
  integer ( kind = 4 ) tab(2,n)
!
!  Check.
!
  call tableau_check ( n, tab )

  do i = 1, 2
    do j = 1, n
      t(tab(i,j)) = i - 1
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
subroutine tree_check ( n, t, ierror )

!*****************************************************************************80
!
!! TREE_CHECK checks a tree.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the tree.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(2,N-1), describes the edges of the tree
!    as pairs of nodes.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    -1, N was illegal.
!    J, the edge T(1,J) to T(2,J) is illegal.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) t(2,n-1)
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y

  ierror = 0

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  do i = 1, 2
    do j = 1, n - 1
      if ( t(i,j) < 1 .or. n < t(i,j) ) then
        ierror = j
        return
      end if
    end do
  end do
!
!  Compute the degree of each node.
!
  call edge_degree ( n, n - 1, t, d )
!
!  Delete a node of degree 1, N-1 times.
!
  do k = 1, n - 1

    x = 1

    do while ( d(x) /= 1 )
      x = x + 1
      if ( n < x ) then
        ierror = -1
        return
      end if
    end do
!
!  Find its neighbor.
!
    j = 1

    do

      if ( t(1,j) == x ) then
        y = t(2,j)
        exit
      end if

      if ( t(2,j) == x ) then
        y = t(1,j)
        exit
      end if

      j = j + 1

      if ( n < j ) then
        ierror = -1
        return
      end if

    end do
!
!  Delete the edge.
!
    t(1,j) = - t(1,j)
    t(2,j) = - t(2,j)

    d(x) = d(x) - 1
    d(y) = d(y) - 1

  end do

  t(1:2,1:n-1) = - t(1:2,1:n-1)

  return
end
subroutine tree_enum ( n, ntree )

!*****************************************************************************80
!
!! TREE_ENUM enumerates the trees on N nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in each tree.
!    N must normally be at least 3, but for this routine,
!    any value of N is allowed.
!
!    Output, integer ( kind = 4 ) NTREE, the number of distinct elements.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ntree

  if ( n < 1 ) then
    ntree = 0
  else if ( n == 1 ) then
    ntree = 1
  else if ( n == 2 ) then
    ntree = 1
  else
    ntree = n**( n - 2 )
  end if

  return
end
subroutine tree_rank ( n, t, rank )

!*****************************************************************************80
!
!! TREE_RANK ranks a tree.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the tree.
!    N must be at least 3.
!
!    Input, integer ( kind = 4 ) T(2,N-1), describes the edges of the tree
!    as pairs of nodes.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the tree.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(n-2)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(2,n-1)
!
!  Check the tree.
!
  call tree_check ( n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_RANK - Fatal error!'
    write ( *, '(a)' ) '  Input tree is illegal.'
    write ( *, '(a,i8)' ) '  Error code = ', ierror
    stop
  end if
!
!  Convert the tree to a Pruefer code.
!
  call tree_to_pruefer ( n, t, p )
!
!  Find the rank of the Pruefer code.
!
  call pruefer_rank ( n, p, rank )

  return
end
subroutine tree_successor ( n, t, rank )

!*****************************************************************************80
!
!! TREE_SUCCESSOR returns the successor of a tree.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the tree.
!    N must be at least 3.
!
!    Input/output, integer ( kind = 4 ) T(2,N-1), describes the edges of the
!    tree as pairs of nodes.  On output, the input tree has been replaced
!    by its successor.
!
!    Input/output, integer RANK, the rank of the tree.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(n-2)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(2,n-1)
!
!  Return the first element.
!
  if ( rank == -1 ) then
    p(1:n-2) = 1
    call pruefer_to_tree ( n, p, t )
    rank = 0
    return
  end if
!
!  Check the tree.
!
  call tree_check ( n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_SUCCESSOR - Fatal error!'
    write ( *, '(a)' ) '  Input tree is illegal.'
    write ( *, '(a,i8)' ) '  Error code = ', ierror
    stop
  end if
!
!  Convert the tree to a Pruefer code.
!
  call tree_to_pruefer ( n, t, p )
!
!  Find the successor of the Pruefer code.
!
  call pruefer_successor ( n, p, rank )
!
!  Convert the Pruefer code to the tree.
!
  call pruefer_to_tree ( n, p, t )

  return
end
subroutine tree_to_pruefer ( n, t, p )

!*****************************************************************************80
!
!! TREE_TO_PRUEFER converts a tree to a Pruefer code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the tree.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) T(2,N-1), describes the edges of the tree
!    as pairs of nodes.
!
!    Output, integer ( kind = 4 ) P(N-2), the Pruefer code for the tree.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) d(n)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(n-2)
  integer ( kind = 4 ) t(2,n-1)
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
!
!  Check.
!
  call tree_check ( n, t, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_TO_PRUEFER - Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal!'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  Compute the degree of each node.
!
  call edge_degree ( n, n - 1, t, d )

  do j = 1, n - 2
!
!  Find a node of degree 1.
!
    x = n
    do while ( d(x) /= 1 )
      x = x - 1
    end do
!
!  Find its neighbor.
!
    k = 1

    do

      if ( t(1,k) == x ) then
        y = t(2,k)
        exit
      end if

      if ( t(2,k) == x ) then
        y = t(1,k)
        exit
      end if

      k = k + 1

    end do
!
!  Store the neighbor.
!
    p(j) = y
!
!  Delete the edge from the tree.
!
    d(x) = d(x) - 1
    d(y) = d(y) - 1

    t(1,k) = - t(1,k)
    t(2,k) = - t(2,k)

  end do
!
!  Remove the negative signs from the first N-2 columns of the tree.
!
  t(1:2,1:n-2) = - t(1:2,1:n-2)

  return
end
subroutine tree_unrank ( rank, n, t )

!*****************************************************************************80
!
!! TREE_UNRANK unranks a tree.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the tree.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the tree.
!    N must be at least 3.
!
!    Output, integer ( kind = 4 ) T(2,N-1), describes the edges of the tree
!    as pairs of nodes.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n-2)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
  integer ( kind = 4 ) t(2,n-1)
  integer ( kind = 4 ) tree_num
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call tree_enum ( n, tree_num )

  if ( rank < 0 .or. tree_num < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if
!
!  Unrank the Pruefer code.
!
  call pruefer_unrank ( rank, n, p )
!
!  Convert the Pruefer code to a tree.
!
  call pruefer_to_tree ( n, p, t )

  return
end
