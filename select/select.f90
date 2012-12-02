function bnew ( n_max, k_max, n, k, family, mgo, b )

!*****************************************************************************80
!
!! BNEW returns either of two values needed to compute B(N,K).
!
!  Discussion:
!
!    The abstract form of the recurrence for B(N,K) is:
!
!      B(N,K) = PHI(N,K) * B(N1,K1) + PSI(N,K) * B(N2,K2)
!
!    where, typically, but not always:
!
!      N1 = N - 1
!      N2 = N - 1
!      K1 = K
!      K2 = K - 1
!
!    This routine is given the values N and K, and returns one
!    of the two values B(N1,K1) or B(N2,K2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer N_MAX, the maximum value of N to be considered.
!
!    Input, integer K_MAX, the maximum value of K to be considered.
!
!    Input, integer N, the size of the family.
!
!    Input, integer K, the size of the object.
!
!    Input, integer FAMILY, specifies the combinatorial family.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional space
!       over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!    Input, integer MGO, is
!    1 to request B(N1,K1), or
!    2 to request B(N2,K2).
!
!    Input, integer B(N_MAX,K_MAX), the precomputed matrix enumerating
!    the combinatorial family.
!
!    Output, integer BNEW, the value of B(N1,K1) or B(N2,K2).
!
  implicit none

  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) b(n_max,k_max)
  integer ( kind = 4 ) bnew
  integer ( kind = 4 ) family
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lx
  integer ( kind = 4 ) ly
  integer ( kind = 4 ) mgo
  integer ( kind = 4 ) n
  integer ( kind = 4 ) xnew

  lx = xnew ( n, k, family, mgo )

  ly = k - mgo + 1

  if ( 0 < lx * ly ) then
    bnew = b(lx,ly)
  else
    bnew = 1
  end if

  return
end
subroutine code_to_composition ( n, k, mu, nu, edge, m, composition )

!*****************************************************************************80
!
!! CODE_TO_COMPOSITION decodes a composition.
!
!  Discussion:
!
!    NCOMP = 0
!    Composition = ()
!
!    For J = M+1 downto 2,
!
!      if NU(J-1) /= NU(J)
!        NCOMP = NCOMP + 1
!        COMPOSITION(2:min(NCOMP,K)) := COMPOSITION(1:min(NCOMP,K)-1)
!        COMPOSITION(1) = 0
!      else
!        COMPOSITION(1) := COMPOSITION(1) + 1
!      end if
!
!    end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer N, the size of the family.
!
!    Input, integer K, the size of the object.
!
!    Input, integer MU(M+1), part of the code for the object.
!
!    Input, integer NU(M+1), part of the code for the object.
!
!    Input, integer EDGE(M), part of the code for the object.
!
!    Input, integer M, the size of the EDGE, MU and NU arrays.
!
!    Output, integer COMPOSITION(K), records the composition.
!    It should be the case that the entries of COMPOSITION
!    are nonnegative, and sum to N.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) composition(k)
  integer ( kind = 4 ) edge(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) mu(m+1)
  integer ( kind = 4 ) nu(m+1)

  composition(1:k) = 0

  ncomp = 1

  do j = m+1, 2, -1

    if ( nu(j-1) /= nu(j) ) then

      ncomp = ncomp + 1

      if ( ncomp <= k ) then
        composition(2:ncomp) = composition(1:ncomp-1)
      else
        composition(2:k) = composition(1:k-1)
      end if

      composition(1) = 0

    else

      composition(1) = composition(1) + 1

    end if

  end do

  if ( ncomp < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CODE_TO_COMPOSITION - Fatal error!'
    write ( *, '(a)' ) '  Fewer than K parts in the composition.'
    write ( *, '(a,i8)' ) '  NCOMP = ', ncomp
    write ( *, '(a,i8)' ) '  K =     ', k
    stop
  end if

  return
end
subroutine code_to_partition_num ( n, k, mu, nu, edge, m, partition_num )

!*****************************************************************************80
!
!! CODE_TO_PARTITION_NUM decodes a numeric partition.
!
!  Discussion:
!
!    Partition = partition whose largest part is 0.
!
!    For J = M+1 downto 2,
!
!      if NU(J-1) /= NU(J)
!        add 1 to the largest part of the partition.
!      else
!        make an extra copy of the largest part of the partition.
!      end if
!
!    end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer N, the size of the family.
!
!    Input, integer K, the size of the object.
!
!    Input, integer MU(M+1), part of the code for the object.
!
!    Input, integer NU(M+1), part of the code for the object.
!
!    Input, integer EDGE(M), part of the code for the object.
!
!    Input, integer M, the size of the EDGE, MU and NU arrays.
!
!    Output, integer PARTITION_NUM(M), records the integers
!    to be added together to get N.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) edge(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mu(m+1)
  integer ( kind = 4 ) nu(m+1)
  integer ( kind = 4 ) np
  integer ( kind = 4 ) partition_num(n)

  partition_num(1:m) = 0

  np = 1

  do j = m+1, 2, -1

    if ( nu(j-1) /= nu(j) ) then

      if ( k <= partition_num(np) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CODE_TO_PARTITION_NUM - Warning!'
        write ( *, '(a)' ) '  An entry of the partition is greater than K.'
      end if

      partition_num(np) = partition_num(np) + 1

    else

      np = np + 1

      partition_num(np) = partition_num(np-1)

    end if

  end do

  if ( partition_num(np) /= k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CODE_TO_PARTITION_NUM - Warning!'
    write ( *, '(a)' ) '  The largest entry of the partition is not K.'
  end if

  if ( sum ( partition_num(1:np) ) /= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CODE_TO_PARTITION_NUM - Warning!'
    write ( *, '(a)' ) '  Partition does not sum to N.'
    stop
  end if

  return
end
subroutine code_to_partition_set ( n, k, mu, nu, edge, m, partition_set )

!*****************************************************************************80
!
!! CODE_TO_PARTITION_SET decodes a set partition.
!
!  Discussion:
!
!    Partition = empty partition of NP = 0 parts.
!
!    For J = M+1 downto 2,
!
!      if NU(J-1) /= NU(J)
!        NP = NP + 1.
!        insert MU(J-1) as a singleton into part NP.
!      else
!        add MU(J-1) to part EDGE(J-1).
!      end if
!
!    end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2003
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer N, the size of the family.
!
!    Input, integer K, the size of the object.
!
!    Input, integer MU(M+1), part of the code for the object.
!
!    Input, integer NU(M+1), part of the code for the object.
!
!    Input, integer EDGE(M), part of the code for the object.
!
!    Input, integer M, the size of the EDGE, MU and NU arrays.
!
!    Output, integer PARTITION_SET(N), records, for each item, the
!    part into which it was placed.  Parts are numbered from 0
!    to K-1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) edge(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mu(m+1)
  integer ( kind = 4 ) nu(m+1)
  integer ( kind = 4 ) np
  integer ( kind = 4 ) partition_set(n)

  partition_set(1:n) = 0

  np = -1

  do j = m+1, 2, -1

    if ( nu(j-1) /= nu(j) ) then

      np = np + 1

      if ( k-1 < np ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CODE_TO_PARTITION_SET - Fatal error!'
        write ( *, '(a)' ) '  Adding more than K parts to partition.'
        stop
      end if

      partition_set(j-1) = np

    else

      partition_set(j-1) = edge(j-1)

    end if

  end do

  if ( np < k-1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CODE_TO_PARTITION_SET - Fatal error!'
    write ( *, '(a)' ) '  Fewer than K parts in the partition.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  K =  ', k
    write ( *, '(a,i8)' ) '  NP = ', np
    write ( *, '(a,i8)' ) '  M =  ', m
    write ( *, '(a)' ) '  MU:'
    write ( *, '(10i8)' ) mu(1:m)
    write ( *, '(a)' ) '  NU:'
    write ( *, '(10i8)' ) nu(1:m)
    write ( *, '(a)' ) '  EDGE:'
    write ( *, '(10i8)' ) edge(1:n)
    stop
  end if

  return
end
subroutine code_to_perm_cycle ( n, k, mu, nu, edge, m, perm_cycle )

!*****************************************************************************80
!
!! CODE_TO_PERM_CYCLE decodes a permutation of N objects with K cycles.
!
!  Discussion:
!
!    Permutation = empty permutation of no cycles.
!
!    For J = M+1 downto 2,
!
!      if NU(J-1) /= NU(J)
!        increase number of cycles by 1;
!        adjoin MU(J-1) as a singleton cycle;
!      else
!        insert MU(J-1) into the the EDGE(J-1)-th "position".
!      end if
!
!    end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer N, the size of the family.
!
!    Input, integer K, the size of the object.
!
!    Input, integer MU(M+1), part of the code for the object.
!
!    Input, integer NU(M+1), part of the code for the object.
!
!    Input, integer EDGE(M), part of the code for the object.
!
!    Input, integer M, the size of the EDGE, MU and NU arrays.
!
!    Output, integer PERM_CYCLE(N), records, for each item,
!    its image under the permutation.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) cycles
  integer ( kind = 4 ) edge(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mu(m+1)
  integer ( kind = 4 ) nu(m+1)
  integer ( kind = 4 ) perm_cycle(n)

  perm_cycle(1:n) = 0

  cycles = 0

  do j = n+1, 2, -1

    if ( nu(j-1) /= nu(j) ) then

      cycles = cycles + 1

      if ( k < cycles ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CODE_TO_PERM_CYCLE - Fatal error!'
        write ( *, '(a)' ) '  Adding more than K cycles to the permutation.'
        stop
      end if

      perm_cycle(n+2-j) = -mu(j-1)

    else

      i = edge(j-1) + 2

      perm_cycle(i+1:n+2-j) = perm_cycle(i:n+1-j)
      perm_cycle(i) = mu(j-1)

    end if

  end do

  if ( cycles < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CODE_TO_PERM_CYCLE - Fatal error!'
    write ( *, '(a)' ) '  Fewer than K cycles in the permutation.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  K =  ', k
    write ( *, '(a,i8)' ) '  CYCLES = ', cycles
    write ( *, '(a,i8)' ) '  M =  ', m
    write ( *, '(a)' ) '  MU:'
    write ( *, '(10i8)' ) mu(1:m)
    write ( *, '(a)' ) '  NU:'
    write ( *, '(10i8)' ) nu(1:m)
    write ( *, '(a)' ) '  EDGE:'
    write ( *, '(10i8)' ) edge(1:n)
    stop
  end if

  return
end
subroutine code_to_perm_runs ( n, k, mu, nu, edge, m, perm )

!*****************************************************************************80
!
!! CODE_TO_PERM_RUNS decodes a permutation of N objects with K runs.
!
!  Discussion:
!
!    Permutation = empty permutation.
!
!    For J = M+1 downto 2,
!
!      if NU(J-1) == NU(J)
!        append MU(J-1) to the end of the EDGE(J-1)th run.
!      else
!        Insert MU(J-1) into the EDGE(J-1)th space of an existing run.
!      end if
!
!    end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer N, the size of the family.
!
!    Input, integer K, the size of the object.
!
!    Input, integer MU(M+1), part of the code for the object.
!
!    Input, integer NU(M+1), part of the code for the object.
!
!    Input, integer EDGE(M), part of the code for the object.
!
!    Input, integer M, the size of the EDGE, MU and NU arrays.
!
!    Output, integer PERM(N); PERM(I) is the image of item I
!    under the permutation.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) count
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) edge(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) len
  integer ( kind = 4 ) mu(m+1)
  integer ( kind = 4 ) nu(m+1)
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) runs

  perm(1:n) = 0

  if ( debug ) then
    write ( *, '(a,1x,8i2)' ) 'EDGE:', edge(1:m)
    write ( *, '(a,   8i2)' ) 'MU:  ', mu(1:m+1)
    write ( *, '(a,   8i2)' ) 'NU:  ', nu(1:m+1)
  end if

  runs = 1

  do j = n+1, 2, -1

    len = n + 1 - j
!
!  Append MU(J-1) at the end of the EDGE(J-1)th run.
!
    if ( nu(j-1) == nu(j) ) then

      i = 0
      count = 0

      do

        i = i + 1

        if ( len < i ) then
          perm(i) = mu(j-1)
          if ( debug ) then
            write ( *,* ) 'E1: ', mu(j-1), ' to position ', i
          end if
          exit
        end if

        if ( 1 < i ) then
          if ( perm(i) < perm(i-1) ) then
            count = count + 1
          end if
        end if

        if ( count == edge(j-1)+1 ) then
          perm(i+1:len+1) = perm(i:len)
          perm(i) = mu(j-1)
          if ( debug ) then
            write ( *,* ) 'E2: ', mu(j-1), ' to position ', i
          end if
          exit
        end if

      end do
!
!  Insert MU(J-1) into the EDGE(J-1)th "interstitial" spot.
!
    else

      count = -1

      do i = 1, len+1

        if ( i == 1 ) then

          if ( edge(j-1) == 0 ) then
            count = count + 1
            perm(2:len+1) = perm(1:len)
            perm(1) = mu(j-1)
            if ( debug ) then
              write ( *,* ) 'N1: ', mu(j-1), ' to position ', i
            end if
            exit
          end if

        else if ( i <= len ) then

          if ( perm(i) < perm(i+1) ) then

            count = count + 1
            if ( edge(j-1) == count ) then
              perm(i+2:len+1) = perm(i+1:len)
              perm(i+1) = mu(j-1)
              if ( debug ) then
                write ( *,* ) 'N2: ', mu(j-1), ' to position ', i+1
              end if
              exit
            end if
          end if

        else

          perm(i) = mu(j-1)
          if ( debug ) then
            write ( *,* ) 'N3: ', mu(j-1), ' to position ', i+1
          end if
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine code_to_subset ( n, k, mu, nu, edge, m, subset )

!*****************************************************************************80
!
!! CODE_TO_SUBSET decodes a subset.
!
!  Discussion:
!
!    Initialize subset to the empty set.
!
!    For J = M+1 downto 2,
!
!      if NU(J-1) /= NU(J)
!        add MU(J-1) to the subset.
!      end if
!
!    end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer N, the size of the family.
!
!    Input, integer K, the size of the object.
!
!    Input, integer MU(M+1), part of the code for the object.
!
!    Input, integer NU(M+1), part of the code for the object.
!
!    Input, integer EDGE(M), is not used for this routine.
!
!    Input, integer M, determines the size of the EDGE, MU and NU arrays.
!
!    Output, integer SUBSET(K), the indices of the items in the subset.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) edge(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) mu(m+1)
  integer ( kind = 4 ) nu(m+1)
  integer ( kind = 4 ) subset(k)

  subset(1:k) = 0

  kk = 0

  do j = m+1, 2, -1

    if ( nu(j-1) /= nu(j) ) then

      kk = kk + 1

      if ( k < kk ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CODE_TO_SUBSET - Fatal error!'
        write ( *, '(a)' ) '  Adding more than K items to subset.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  K =  ', k
        write ( *, '(a,i8)' ) '  KK = ', kk
        write ( *, '(a,i8)' ) '  M =  ', m
        call i4vec_print ( m+1, mu, '  MU' )
        call i4vec_print ( m+1, nu, '  NU' )
        stop
      end if

      subset(kk) = mu(j-1)

    end if
  end do

  if ( kk < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CODE_TO_SUBSET - Fatal error!'
    write ( *, '(a)' ) '  Added fewer than K items to subset.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  K =  ', k
    write ( *, '(a,i8)' ) '  KK = ', kk
    write ( *, '(a,i8)' ) '  M =  ', m
    call i4vec_print ( m+1, mu, '  MU' )
    call i4vec_print ( m+1, nu, '  NU' )
    stop
  end if

  return
end
subroutine code_to_subspace ( n, k, mu, nu, edge, m, subspace )

!*****************************************************************************80
!
!! CODE_TO_SUBSPACE decodes a subspace.
!
!  Discussion:
!
!    NR = 0.
!    NC = 0.
!    Initialize subspace to the empty NR by NC matrix.
!
!    For J = M+1 downto 2,
!
!      if NU(J-1) /= NU(J)
!        NR = NR + 1;
!        NC = NC + 1;
!        append row NR and column NC, all zero except the (NR,NC)
!        entry, which is 1.
!      else
!        NC = NC + 1;
!        the entries of column NC are the digits of EDGE(J-1) base Q.
!      end if
!
!    end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer N, the size of the family.
!
!    Input, integer K, the size of the object.
!
!    Input, integer MU(M+1), part of the code for the object.
!
!    Input, integer NU(M+1), part of the code for the object.
!
!    Input, integer EDGE(M), part of the code for the object.
!
!    Input, integer M, the size of the EDGE, MU and NU arrays.
!
!    Output, integer SUBSPACE(K,N), contains K rows which
!    span a subspace of the vector space of dimension N,
!    over the field Q.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) e
  integer ( kind = 4 ) edge(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mu(m+1)
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nu(m+1)
  integer ( kind = 4 ), parameter :: q = 2
  integer ( kind = 4 ) subspace(k,n)

  subspace(1:k,1:n) = 0

  nr = 0
  nc = 0

  do j = m+1, 2, -1

    if ( nu(j-1) /= nu(j) ) then

      nr = nr + 1
      nc = nc + 1

      if ( k < nr ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CODE_TO_SUBSPACE - Fatal error!'
        write ( *, '(a)' ) '  Adding more than K rows to subspace matrix.'
        stop
      end if

      if ( n < nc ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CODE_TO_SUBSPACE - Fatal error!'
        write ( *, '(a)' ) '  Adding more than N columns to subspace matrix.'
        stop
      end if

      subspace(nr,nc) = 1

    else

      nc = nc + 1
      e = edge(j-1)

      do i = 1, nr
        subspace(i,nc) = mod ( e, q )
        e = e / q
      end do

    end if

  end do

  if ( nr < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CODE_TO_SUBSPACE - Fatal error!'
    write ( *, '(a)' ) '  Added fewer than K rows to subspace matrix.'
    stop
  end if

  if ( nc < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CODE_TO_SUBSPACE - Fatal error!'
    write ( *, '(a)' ) '  Added fewer than N columns to subspace matrix.'
    stop
  end if

  return
end
subroutine composition_to_string ( n, k, composition, s )

!*****************************************************************************80
!
!! COMPOSITION_TO_STRING writes a character representation of a composition.
!
!  Example:
!
!    N = 5
!    K = 3
!    COMPOSITION = ( 2, 0, 3 )
!
!    S = '5 = 2|0|3'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the set.
!
!    Input, integer K, the number of parts.
!
!    Input, integer COMPOSITION(K), the parts of the composition.
!
!    Output, character ( len = * ) S, a representation of the composition.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) composition(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  character ( len = * ) s

  s = ' '

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMPOSITION_TO_STRING - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( k < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMPOSITION_TO_STRING - Fatal error!'
    write ( *, '(a)' ) '  K < 1.'
    stop
  end if

  i = 0

  i = i + 1
  call i4_to_s_left ( n, s(i:) )
  i = len_trim ( s )

  i = i + 1
  s(i:i) = '='

  do j = 1, k

    if ( 1 < j ) then
      i = i + 1
      s(i:i) = '|'
    end if

    i = i + 1
    call i4_to_s_left ( composition(j), s(i:) )
    i = len_trim ( s )

  end do

  return
end
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
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
!    Input, integer DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine family_enumerate ( family, n_max, k_max, b )

!*****************************************************************************80
!
!! FAMILY_ENUMERATE produces an enumeration table for a given family.
!
!  Discussion:
!
!    The table B(1:N_MAX,1:K_MAX) is computed.  Entry B(N,K)
!    contains the number of objects of size K in the family of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer FAMILY, specifies the combinatorial family.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional space
!       over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!    Input, integer N_MAX, the maximum size of N to be considered.
!
!    Input, integer K_MAX, the maximum size of K to be considered.
!
!    Output, integer B(N_MAX,K_MAX), the table.
!
  implicit none

  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) b(n_max,k_max)
  integer ( kind = 4 ) bnew
  integer ( kind = 4 ) family
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) phi
  integer ( kind = 4 ) psi

  do n = 1, n_max
    do k = 1, k_max

      b(n,k) = phi ( n, k, family ) &
             * bnew ( n_max, k_max, n, k, family, 1, b ) &
             + psi ( n, k, family ) &
             * bnew ( n_max, k_max, n, k, family, 2, b )

    end do
  end do

  return
end
subroutine family_next ( family, n_max, k_max, b, n, k, newone, rank, m, mu, &
  nu, edge )

!*****************************************************************************80
!
!! FAMILY_NEXT carries out a task for a combinatorial family of order N, K.
!
!  Discussion:
!
!    Typically, this routine is called to return one object at a time from
!    the family.  The objects are returned in rank order.
!
!    Before the first call, the user should set N_MAX, K_MAX, N, and K,
!    set up the array B by calling FAMILY ENUMERATE, and set NEWONE to FALSE.
!
!    Upon calling the routine, the value of NEWONE will be reset to TRUE,
!    indicating that information about the first object is returned in
!    RANK, M, MU, NU, EDGE.
!
!    To get the next object, call again, passing in as input arguments
!    the same values of RANK, M, MU, NU and EDGE.  The routine needs these
!    previous values to compute the next ones.
!
!    Each subsequent call returns the next object, until NEWONE is returned with
!    the value FALSE, in which case there are no more objects (and so the final
!    object was returned on the PREVIOUS call.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer FAMILY, specifies the combinatorial family.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional space
!       over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!    Input, integer N_MAX, the maximum size of N to be considered.
!
!    Input, integer K_MAX, the maximum size of K to be considered.
!
!    Input, integer B(N_MAX,K_MAX), the enumeration table for the family.
!    The user must evaluate this table by calling FAMILY_ENUMERATE.
!
!    Input, integer N, the size of the family.
!    1 <= N <= N_MAX.
!
!    Input, integer K, the size of the object.
!    1 <= K <= K_MAX.
!
!    Input/output, logical NEWONE.
!    An input value of FALSE indicates startup, and an input value
!    of TRUE indicates a followup call.  An output value of TRUE indicates
!    a next object was returned, but an output value of FALSE indicates
!    there are no more objects.
!
!    Input/output, integer RANK, the rank of an object.
!
!    Input/output, integer M, the length of a code.
!
!    Input/output, integer MU(M+1), part of a code.
!
!    Input/output, integer NU(M+1), part of a code.
!
!    Input/output, integer EDGE(M), part of a code.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) b(n_max,k_max)
  integer ( kind = 4 ) b1
  integer ( kind = 4 ) b2
  integer ( kind = 4 ) bnew
  integer ( kind = 4 ) edge(n+k-1)
  integer ( kind = 4 ) family
  integer, save :: family_save = 0
  integer ( kind = 4 ) j
  integer, save :: k_max_save = 0
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mu(n+k)
  logical newone
  integer ( kind = 4 ) nu(n+k)
  integer, save :: n_max_save = 0
  integer ( kind = 4 ) phi
  integer ( kind = 4 ) psi
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rankp
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2
  integer ( kind = 4 ) xnew

  if ( .not. newone ) then
!
!  Check N.
!
    if ( n < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAMILY_NEXT - Fatal error!'
      write ( *, '(a)' ) '  N < 1.'
      write ( *, '(a,i12)' ) '  N = ', n
      stop
    end if

    if ( n_max < n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAMILY_NEXT - Fatal error!'
      write ( *, '(a)' ) '  N_MAX < N.'
      write ( *, '(a,i12)' ) '  N = ', n
      write ( *, '(a,i12)' ) '  N_MAX = ', n_max
      stop
    end if
!
!  Check K.
!
    if ( k < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAMILY_NEXT - Fatal error!'
      write ( *, '(a)' ) '  K < 1.'
      write ( *, '(a,i12)' ) '  K = ', k
      stop
    end if

    if ( n < k ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAMILY_NEXT - Fatal error!'
      write ( *, '(a)' ) '  N < K.'
      write ( *, '(a,i12)' ) '  N = ', n
      write ( *, '(a,i12)' ) '  K = ', k
      stop
    end if

    if ( k_max < k ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAMILY_NEXT - Fatal error!'
      write ( *, '(a)' ) '  K_MAX < K.'
      write ( *, '(a,i12)' ) '  K = ', k
      write ( *, '(a,i12)' ) '  K_MAX = ', k_max
      stop
    end if

    newone = .true.
    rank = 0
    rankp = rank
    j = 1
    m = j
    mu(j) = n
    nu(j) = k

  else

    rank = rank + 1
    j = m
!
!  Starting at J = M, seek EDGE(J) that is movable.
!  But if we reach J = 0, there are no more, and we
!  have exhausted the family list.
!
    do

      t1 = phi ( mu(j), nu(j), family )
      t2 = psi ( mu(j), nu(j), family )

      if ( edge(j) + 1 < t1 + t2 ) then
        exit
      end if

      j = j - 1

      if ( j == 0 ) then
        newone = .false.
        return
      end if

    end do

    edge(j) = edge(j) + 1
    l = j + 1

    if ( edge(j) == t1 ) then
      nu(l) = nu(l) - 1
      mu(l) = xnew ( mu(l-1), nu(l-1), family, 2 )
    end if

    rankp = 0
    m = l

  end if

  do

    t1 = phi ( mu(m), nu(m), family )

    t2 = psi ( mu(m), nu(m), family )

    if ( t1 + t2 == 0 ) then
      m = m - 1
      exit
    end if

    b1 = bnew ( n_max, k_max, mu(m), nu(m), family, 1, b )

    if ( rankp < t1 * b1 ) then

      edge(m) = rankp / b1
      rankp = rankp - edge(m) * b1
      mu(m+1) = xnew ( mu(m), nu(m), family, 1 )
      nu(m+1) = nu(m)

    else

      rankp = rankp - b1 * t1
      b2 = bnew ( n_max, k_max, mu(m), nu(m), family, 2, b )
      edge(m) = t1 + rankp / b2
      rankp = rankp - ( edge(m) - t1 ) * b2
      mu(m+1) = xnew ( mu(m), nu(m), family, 2 )
      nu(m+1) = nu(m) - 1

    end if

    m = m + 1

  end do

  return
end
subroutine family_rank ( family, n_max, k_max, b, n, k, m, mu, nu, edge, &
  rank )

!*****************************************************************************80
!
!! FAMILY_RANK returns the rank of an object of a given family.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer FAMILY, specifies the combinatorial family.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional space
!       over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!    Input, integer N_MAX, the maximum size of N to be considered.
!
!    Input, integer K_MAX, the maximum size of K to be considered.
!
!    Input, integer B(N_MAX,K_MAX), the enumeration table for the family.
!    The user must evaluate this table by calling FAMILY_ENUMERATE.
!
!    Input, integer N, the size of the family.
!    1 <= N <= N_MAX.
!
!    Input, integer K, the size of the object.
!    1 <= K <= K_MAX.
!
!    Input, integer M, the length of the code.  For
!    families 1 through 5, M is always equal to N.
!    For family 6, M varies, but is never greater than N.
!    For family 7, M = N+K-1.
!
!    Input, integer MU(M+1), part of the code for the object.
!
!    Input, integer NU(M+1), part of the code for the object.
!
!    Input, integer EDGE(M), part of the code for the object.
!
!    Output, integer RANK, the rank of the object.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) b(n_max,k_max)
  integer ( kind = 4 ) bnew
  integer ( kind = 4 ) edge(n+k-1)
  integer ( kind = 4 ) family
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mu(n+k)
  integer ( kind = 4 ) nu(n+k)
  integer ( kind = 4 ) phi
  integer ( kind = 4 ) rank
!
!  Check N.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_RANK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i12)' ) '  N = ', n
    stop
  end if

  if ( n_max < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_RANK - Fatal error!'
    write ( *, '(a)' ) '  N_MAX < N.'
    write ( *, '(a,i12)' ) '  N = ', n
    write ( *, '(a,i12)' ) '  N_MAX = ', n_max
    stop
  end if
!
!  Check K.
!
  if ( k < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_RANK - Fatal error!'
    write ( *, '(a)' ) '  K < 1.'
    write ( *, '(a,i12)' ) '  K = ', k
    stop
  end if

  if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_RANK - Fatal error!'
    write ( *, '(a)' ) '  N < K.'
    write ( *, '(a,i12)' ) '  N = ', n
    write ( *, '(a,i12)' ) '  K = ', k
    stop
  end if

  if ( k_max < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_RANK - Fatal error!'
    write ( *, '(a)' ) '  K_MAX < K.'
    write ( *, '(a,i12)' ) '  K = ', k
    write ( *, '(a,i12)' ) '  K_MAX = ', k_max
    stop
  end if
!
!  Should check M, MU, NU and EDGE!
!
  rank = 0

  do j = 1, m - 1

    if ( nu(j+1) == nu(j) ) then
      rank = rank + edge(j) &
        * bnew ( n_max, k_max, mu(j), nu(j), family, 1, b )
    else
      rank = rank &
        + phi ( mu(j), nu(j), family ) &
        * bnew ( n_max, k_max, mu(j), nu(j), family, 1, b ) &
        + ( edge(j) - phi ( mu(j), nu(j), family) ) &
        * bnew ( n_max, k_max, mu(j), nu(j), family, 2, b )
    end if

  end do

  return
end
subroutine family_sample ( family, n_max, k_max, b, n, k, seed, rank, m, mu, &
  nu, edge )

!*****************************************************************************80
!
!! FAMILY_SAMPLE produces a sample object of a given family.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer FAMILY, specifies the combinatorial family.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional space
!       over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!    Input, integer N_MAX, the maximum size of N to be considered.
!
!    Input, integer K_MAX, the maximum size of K to be considered.
!
!    Input, integer B(N_MAX,K_MAX), the enumeration table for the family.
!    The user must evaluate this table by calling FAMILY_ENUMERATE.
!
!    Input, integer N, the size of the family.
!    1 <= N <= N_MAX.
!
!    Input, integer K, the size of the object.
!    1 <= K <= K_MAX.
!
!    Input/output, integer SEED, a seed for the random number generator.
!
!    Output, integer RANK, the rank of the object.
!
!    Output, integer M, the length of the code.  For
!    families 1 through 5, M is always equal to N.
!    For family 6, M varies, but is never greater than N.
!    For family 7, M = N+K-1.
!
!    Output, integer MU(M+1), part of the code for the object.
!
!    Output, integer NU(M+1), part of the code for the object.
!
!    Output, integer EDGE(M), part of the code for the object.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) b(n_max,k_max)
  integer ( kind = 4 ) b1
  integer ( kind = 4 ) b2
  integer ( kind = 4 ) bnew
  real ( kind = 8 ) d
  integer ( kind = 4 ) edge(n+k-1)
  integer ( kind = 4 ) family
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mu(n+k)
  integer ( kind = 4 ) nu(n+k)
  integer ( kind = 4 ) phi
  integer ( kind = 4 ) psi
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rankp
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2
  integer ( kind = 4 ) xnew
!
!  Check N.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i12)' ) '  N = ', n
    stop
  end if

  if ( n_max < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  N_MAX < N.'
    write ( *, '(a,i12)' ) '  N = ', n
    write ( *, '(a,i12)' ) '  N_MAX = ', n_max
    stop
  end if
!
!  Check K.
!
  if ( k < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  K < 1.'
    write ( *, '(a,i12)' ) '  K = ', k
    stop
  end if

  if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  N < K.'
    write ( *, '(a,i12)' ) '  N = ', n
    write ( *, '(a,i12)' ) '  K = ', k
    stop
  end if

  if ( k_max < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  K_MAX < K.'
    write ( *, '(a,i12)' ) '  K = ', k
    write ( *, '(a,i12)' ) '  K_MAX = ', k_max
    stop
  end if
!
!  Choose RANK between 0 and B(N,K)-1.
!
  d = r8_uniform_01 ( seed )

  rank = int ( d * real ( b(n,k), kind = 8 ) )
!
!  Construct the object.
!
  rankp = rank
  j = 1
  m = j
  mu(j) = n
  nu(j) = k

  do

    t1 = phi ( mu(m), nu(m), family )

    t2 = psi ( mu(m), nu(m), family )

    if ( t1 + t2 == 0 ) then
      m = m - 1
      exit
    end if

    b1 = bnew ( n_max, k_max, mu(m), nu(m), family, 1, b )

    if ( rankp < t1 * b1 ) then

      edge(m) = rankp / b1
      rankp = rankp - edge(m) * b1
      mu(m+1) = xnew ( mu(m), nu(m), family, 1 )
      nu(m+1) = nu(m)

    else

      rankp = rankp - b1 * t1
      b2 = bnew ( n_max, k_max, mu(m), nu(m), family, 2, b )
      edge(m) = t1 + rankp / b2
      rankp = rankp - ( edge(m) - t1 ) * b2
      mu(m+1) = xnew ( mu(m), nu(m), family, 2 )
      nu(m+1) = nu(m) - 1

    end if

    m = m + 1

  end do

  return
end
subroutine family_unrank ( family, n_max, k_max, b, n, k, rank, m, mu, nu, &
  edge )

!*****************************************************************************80
!
!! FAMILY_UNRANK produces an object of given rank of a given family.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer FAMILY, specifies the combinatorial family.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional space
!       over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!    Input, integer N_MAX, the maximum size of N to be considered.
!
!    Input, integer K_MAX, the maximum size of K to be considered.
!
!    Input, integer B(N_MAX,K_MAX), the enumeration table for the family.
!    The user must evaluate this table by calling FAMILY_ENUMERATE.
!
!    Input, integer N, the size of the family.
!    1 <= N <= N_MAX.
!
!    Input, integer K, the size of the object.
!    1 <= K <= K_MAX.
!
!    Input, integer RANK, the rank of the object.
!
!    Output, integer M, the length of the code.  For
!    families 1 through 5, M is always equal to N.
!    For family 6, M varies, but is never greater than N.
!    For family 7, M = N+K-1.
!
!    Output, integer MU(M+1), part of the code for the object.
!
!    Output, integer NU(M+1), part of the code for the object.
!
!    Output, integer EDGE(M), part of the code for the object.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) b(n_max,k_max)
  integer ( kind = 4 ) b1
  integer ( kind = 4 ) b2
  integer ( kind = 4 ) bnew
  integer ( kind = 4 ) edge(n+k-1)
  integer ( kind = 4 ) family
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mu(n+k)
  integer ( kind = 4 ) nu(n+k)
  integer ( kind = 4 ) phi
  integer ( kind = 4 ) psi
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rankp
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2
  integer ( kind = 4 ) xnew
!
!  Check N.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i12)' ) '  N = ', n
    stop
  end if

  if ( n_max < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  N_MAX < N.'
    write ( *, '(a,i12)' ) '  N = ', n
    write ( *, '(a,i12)' ) '  N_MAX = ', n_max
    stop
  end if
!
!  Check K.
!
  if ( k < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  K < 1.'
    write ( *, '(a,i12)' ) '  K = ', k
    stop
  end if

  if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  N < K.'
    write ( *, '(a,i12)' ) '  N = ', n
    write ( *, '(a,i12)' ) '  K = ', k
    stop
  end if

  if ( k_max < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  K_MAX < K.'
    write ( *, '(a,i12)' ) '  K = ', k
    write ( *, '(a,i12)' ) '  K_MAX = ', k_max
    stop
  end if
!
!  Check RANK.
!
  if ( rank < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  RANK < 0.'
    write ( *, '(a,i12)' ) '  RANK = ', rank
    stop
  end if

  if ( b(n,k) <= rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FAMILY_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  B(N,K) < RANK.'
    write ( *, '(a,i12)' ) '  B(N,K) = ', b(n,k)
    write ( *, '(a,i12)' ) '  RANK = ', rank
    stop
  end if
!
!  Construct the object.
!
  rankp = rank
  j = 1
  m = j
  mu(j) = n
  nu(j) = k

  do

    t1 = phi ( mu(m), nu(m), family )

    t2 = psi ( mu(m), nu(m), family )

    if ( t1 + t2 == 0 ) then
      m = m - 1
      exit
    end if

    b1 = bnew ( n_max, k_max, mu(m), nu(m), family, 1, b )

    if ( rankp < t1 * b1 ) then

      edge(m) = rankp / b1
      rankp = rankp - edge(m) * b1
      mu(m+1) = xnew ( mu(m), nu(m), family, 1 )
      nu(m+1) = nu(m)

    else

      rankp = rankp - b1 * t1
      b2 = bnew ( n_max, k_max, mu(m), nu(m), family, 2, b )
      edge(m) = t1 + rankp / b2
      rankp = rankp - ( edge(m) - t1 ) * b2
      mu(m+1) = xnew ( mu(m), nu(m), family, 2 )
      nu(m+1) = nu(m) - 1

    end if

    m = m + 1

  end do

  return
end
subroutine i4_to_s_left ( intval, s )

!*****************************************************************************80
!
!! I4_TO_S_LEFT converts an integer to a left-justified string.
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  1
!        -1  -1
!         0  0
!      1952  1952
!    123456  123456
!   1234567  ******  <-- Not enough room!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be left-justified.  If there is not enough space,
!    the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * ) s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Find the last digit of IVAL, strip it off, and stick it into the string.
!
  do

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

    if ( ival == 0 ) then
      exit
    end if

  end do
!
!  Shift the string to the left.
!
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '

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
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) ( kind = 4 ) a
  integer ( kind = 4 ) ( kind = 4 ) b
  integer ( kind = 4 ) ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) ( kind = 4 ) seed
  integer ( kind = 4 ) ( kind = 4 ) value

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
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an integer matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, integer A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

  return
end
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an integer matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 7 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7)') j
    end do

    write ( *, '(''  Col '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a7)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an integer vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) big
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(2x,i8,2x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(2x,i8,2x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(2x,i8,2x,i12)' ) i, a(i)
    end do
  end if

  return
end
subroutine partition_num_to_string ( n, k, m, partition_num, s )

!*****************************************************************************80
!
!! PARTITION_NUM_TO_STRING represents a numeric partition.
!
!  Example:
!
!    N = 5
!    K = 3
!    M = 3
!    PARTITION_NUM = ( 1, 1, 3 )
!
!    S = '5=1+1+3'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the set.
!
!    Input, integer K, the size of the largest part.
!
!    Input, integer M, the number of parts.
!
!    Input, integer PARTITION_NUM(M) contains the (positive) summands.
!
!    Output, character ( len = * ) S, a representation of the partition.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) partition_num(m)
  character ( len = * ) s

  s = ' '

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARTITION_NUM_TO_STRING - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( k < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARTITION_NUM_TO_STRING - Fatal error!'
    write ( *, '(a)' ) '  K < 1.'
    stop
  end if

  i = 0

  i = i + 1
  call i4_to_s_left ( n, s(i:) )
  i = len_trim ( s )

  i = i + 1
  s(i:i) = '='

  do j = 1, m

    if ( partition_num(j) == 0 ) then
      if ( n /= 0 .or. 1 < j ) then
        cycle
      end if
    end if

    if ( 1 < j ) then
      i = i + 1
      s(i:i) = '+'
    end if

    i = i + 1
    call i4_to_s_left ( partition_num(j), s(i:) )
    i = len_trim ( s )

  end do

  return
end
subroutine partition_set_to_string ( n, k, partition_set, s )

!*****************************************************************************80
!
!! PARTITION_SET_TO_STRING writes a character representation of a set partition.
!
!  Example:
!
!    N = 5
!    K = 3
!    PARTITION_SET = ( 0, 1, 1, 2, 0 )
!
!    S = '(1,5)(2,3)(4)'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the set.
!
!    Input, integer K, the number of parts.
!
!    Input, integer PARTITION_SET(N), records, for each partition,
!    the index of the part to which it is assigned.  Parts are numbered
!    from 0 to K-1.
!
!    Output, character ( len = * ) S, a representation of the partition.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) part
  integer ( kind = 4 ) partition_set(n)
  character ( len = * ) s
  integer ( kind = 4 ) size

  s = ' '

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARTITION_SET_TO_STRING - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( k < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARTITION_SET_TO_STRING - Fatal error!'
    write ( *, '(a)' ) '  K < 1.'
    stop
  end if

  i = 0

  do part = 0, k-1

    i = i + 1
    s(i:i) = '('

    size = 0

    do j = 1, n

      if ( partition_set(j) == part ) then

        size = size + 1

        if ( 1 < size ) then
          i = i + 1
          s(i:i) = ','
        end if

        i = i + 1
        call i4_to_s_left ( j, s(i:) )
        i = len_trim ( s )

      end if

    end do

    i = i + 1
    s(i:i) = ')'

  end do

  return
end
subroutine perm_cycle_to_string ( n, k, perm_cycle, s )

!*****************************************************************************80
!
!! PERM_CYCLE_TO_STRING makes a string out of a permutation in cycle form.
!
!  Example:
!
!    N = 5
!    K = 3
!
!    PERM_CYCLE = ( -1, 5, 4, -2, -3 )
!
!    S = '(1,5,4)(2)(3)'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the set.
!
!    Input, integer K, the number of parts.
!
!    Input, integer PERM_CYCLE(N), lists the elements being permuted,
!    with the first element of a cycle negative, and the successive
!    images following immediately.
!
!    Output, character ( len = * ) S, a representation of the partition.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) perm_cycle(n)
  character ( len = * ) s

  s = ' '

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_CYCLE_TO_STRING - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( k < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_CYCLE_TO_STRING - Fatal error!'
    write ( *, '(a)' ) '  K < 1.'
    stop
  end if

  i = 0

  do j = 1, n

    if ( perm_cycle(j) < 0 ) then

      if ( 1 < j ) then
        i = i + 1
        s(i:i) = ')'
      end if

      i = i + 1
      s(i:i) = '('
      i = i + 1
      call i4_to_s_left ( -perm_cycle(j), s(i:) )

    else

      i = i + 1
      s(i:i) = ','
      i = i + 1
      call i4_to_s_left ( perm_cycle(j), s(i:) )

    end if

    i = len_trim ( s )

  end do

  i = i + 1
  s(i:i) = ')'

  return
end
function phi ( n, k, family )

!*****************************************************************************80
!
!! PHI returns the coefficient PHI(N,K) in the recurrence.
!
!  Discussion:
!
!    The abstract form of the recurrence is:
!
!      B(N,K) = PHI(N,K) * B(N1,K1) + PSI(N,K) * B(N2,K2)
!
!    For FAMILY = 1, the recurrence is
!
!      C(N,K) = C(N-1,K) + C(N-1,K-1)
!
!      so PHI(N,K) is always 1, unless K < 0 or N-1 < K.
!
!    For FAMILY = 2, the recurrence is
!
!      S1(N,K) = K * S1(N-1,K) + S1(N-1,K-1)
!
!    For FAMILY = 3, the recurrence is
!
!      S2(N,K) = (N-1) * S2(N-1,K) + S2(N-1,K-1)
!
!    For FAMILY = 4, the recurrence is
!
!      E(N,K,Q) = Q**K * E(N-1,K,Q) + E(N-1,K-1,Q)
!
!    For FAMILY = 5, the recurrence is
!
!      <N,K> = K * <N-1,K> + (N-K+1) * <N-1,K-1>
!
!    For FAMILY = 6, the recurrence is
!
!      P(N,K) = P(N-K,K) + P(N-1,K-1)
!
!    For FAMILY = 7, the recurrence is
!
!      B(N,K) = B(N-1,K) + B(N,K-1)
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
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer N, the size of the family.
!
!    Input, integer K, the size of the object.
!
!    Input, integer FAMILY, specifies the combinatorial family.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional space
!       over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!    Output, integer PHI, the value of the PHI coefficient
!    in the recurrence.
!
  implicit none

  integer ( kind = 4 ) family
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) phi
  integer ( kind = 4 ), parameter :: q = 2

  if ( family == 1 ) then

    if ( 0 <= k .and. k <= n-1 ) then
      phi = 1
    else
      phi = 0
    end if

  else if ( family == 2 ) then

    if ( 1 <= k .and. k <= n-1 ) then
      phi = k
    else
      phi = 0
    end if

  else if ( family == 3 ) then

    if ( 0 <= k .and. k <= n-1 ) then
      phi = n - 1
    else
      phi = 0
    end if

  else if ( family == 4 ) then

    if ( 0 <= k .and. k <= n-1 ) then
      phi = q**k
    else
      phi = 0
    end if

  else if ( family == 5 ) then

    if ( 1 <= k .and. k+1 <= n ) then
      phi = k
    else if ( k == 1 .and. n == 1 ) then
      phi = 1
    else
      phi = 0
    end if

  else if ( family == 6 ) then

    if ( 1 <= k .and. 2 * k <= n ) then
      phi = 1
    else
      phi = 0
    end if

  else if ( family == 7 ) then

    if ( 1 <= k .and. 1 <= n ) then
      phi = 1
    else
      phi = 0
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PHI - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of FAMILY.'
    stop

  end if

  return
end
function psi ( n, k, family )

!*****************************************************************************80
!
!! PSI returns the coefficient PSI(N,K) in the recurrence.
!
!  Discussion:
!
!    The abstract form of the recurrence is:
!
!      B(N,K) = PHI(N,K) * B(N1,K1) + PSI(N,K) * B(N2,K2)
!
!    For FAMILY = 1, the recurrence is
!
!      C(N,K) = C(N-1,K) + C(N-1,K-1)
!
!      so PSI(N,K) is always 1, unless K-1 < 0 or N-1 < K-1.
!
!    For FAMILY = 2, the recurrence is
!
!      S1(N,K) = K * S1(N-1,K) + S1(N-1,K-1)
!
!    For FAMILY = 3, the recurrence is
!
!      S2(N,K) = (N-1) * S2(N-1,K) + S2(N-1,K-1)
!
!    For FAMILY = 4, the recurrence is
!
!      E(N,K,Q) = Q**K * E(N-1,K,Q) + E(N-1,K-1,Q)
!
!    For FAMILY = 5, the recurrence is
!
!      <N,K> = K * <N-1,K> + (N-K+1) * <N-1,K-1>
!
!    For FAMILY = 6, the recurrence is
!
!      P(N,K) = P(N-K,K) + P(N-1,K-1)
!
!    For FAMILY = 7, the recurrence is
!
!      B(N,K) = B(N-1,K) + B(N,K-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2003
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer N, the size of the family.
!
!    Input, integer K, the size of the object.
!
!    Input, integer FAMILY, specifies the combinatorial family.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional space
!       over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!    Output, integer PSI, the value of the PSI coefficient
!    in the recurrence.
!
  implicit none

  integer ( kind = 4 ) family
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) psi

  if ( family == 1 )then

    if ( 0 <= k-1 .and. k-1 <= n-1 ) then
      psi = 1
    else
      psi = 0
    end if

  else if ( family == 2 ) then

    if ( k == 1 .and. n == 1 ) then
      psi = 1
    else if ( 2 <= k .and. k <= n ) then
      psi = 1
    else
      psi = 0
    end if

  else if ( family == 3 ) then

    if ( 2 <= k .and. k <= n ) then
      psi = 1
    else if ( k == 1 .and. n == 1 ) then
      psi = 1
    else
      psi = 0
    end if

  else if ( family == 4 ) then

    if ( 1 <= k .and. k <= n ) then
      psi = 1
    else
      psi = 0
    end if

  else if ( family == 5 ) then

    if ( 2 <= k .and. k <= n ) then
      psi = n - k + 1
    else
      psi = 0
    end if

  else if ( family == 6 ) then

    if ( k == 1 .and. n == 1 ) then
      psi = 1
    else if ( 2 <= k .and. k <= n ) then
      psi = 1
    else
      psi = 0
    end if

  else if ( family == 7 ) then

    if ( 2 <= k .and. 0 <= n ) then
      psi = 1
    else
      psi = 0
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PSI - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of FAMILY.'
    stop

  end if

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
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
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine select ( family, task, n_max, k_max, n, k, mu, nu, edge, m, newone, &
  rank, seed, b )

!*****************************************************************************80
!
!! SELECT carries out a task for a combinatorial family of order N, K.
!
!  Discussion:
!
!    Note that the algorithm, as published, specifies a dimension of
!    N for MU and NU.  This is incorrect.  MU and NU must be dimensioned
!    at least N+1, and in some cases, larger, otherwise illegal memory
!    references will result.
!
!    On the first call for a given set of values of FAMILY, N_MAX and K_MAX,
!    the array B is computed.  To save time, it is assumed that on subsequent
!    calls, this array is still available, and has not been changed by the user.
!    If any of FAMILY, N_MAX or K_MAX change on a call, then the routine
!    recomputes B.
!
!    For task #1, listing all the items in a family:
!
!      Set N_MAX, K_MAX, N and K, and set NEWONE = .false.
!
!      SELECT will return MU, NU, EDGE, M, NEWONE, RANK, B.
!
!      If NEWONE is false on return, then no object was returned.
!      Otherwise, MU, NU, EDGE, M and RANK describe an object, which
!      you may print or examine.  Pass these values back to SELECT
!      for another call, and you will get the next object.
!
!    For task #2, ranking a given object:
!
!      Set N_MAX, K_MAX, N, K, MU, NU, EDGE, M
!
!      SELECT will return RANK, B.
!
!    For task #3, producing an object of given rank:
!
!      Set N_MAX, K_MAX, N, K, RANK
!
!      SELECT will return MU, NU, EDGE, M, B.
!
!    For task #4, producing a random object:
!
!      Set N_MAX, K_MAX, N, K, SEED.
!
!      SELECT will return MU, NU, EDGE, M, RANK, B.
!
!    For task #5, enumerating the objects:
!
!      Set N_MAX, K_MAX.
!
!      SELECT will return the matrix B containing the desired table.
!
!  Errata:
!
!    In the reference, in the printed results for Family 1, the
!    values of NU seem to be wrong.  The correct values are:
!
!      Rank   Edge  Mu      Nu      Subset    M
!
!         0  00000  543210  333210  {1,2,3}   5
!         1  01000  543210  332210  {1,2,4}   5
!         2  01100  543210  332110  {1,3,4}   5
!         3  01110  543210  332100  {2,3,4}   5
!         4  10000  543210  322210  {1,2,5}   5
!         5  10100  543210  322110  {1,3,5}   5
!         6  10110  543210  322100  {2,3,5}   5
!         7  11000  543210  321110  {1,4,5}   5
!         8  11010  543210  321100  {2,4,5}   5
!         9  11100  543210  321000  {3,4,5}   5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer FAMILY, specifies the combinatorial family.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional space
!       over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!    Input, integer TASK, specifies the task.
!    1, Present each object of the family, one at a time.
!    2, Rank a given object of the family.
!    3, Produce an object of given rank.
!    4, Select an object at random.
!    5, Enumerate the objects.
!
!    Input, integer N_MAX, the maximum size of N to be considered.
!
!    Input, integer K_MAX, the maximum size of K to be considered.
!
!    Input, integer N, the size of the family.
!    Not required for TASK = 5.
!    1 <= N <= N_MAX otherwise.
!
!    Input, integer K, the size of the object.
!    Not required for TASK = 5.
!    1 <= K <= K_MAX otherwise.
!
!    Input/output, integer MU(M+1), part of the code for the object.
!    The use and meaning of MU depend on TASK and FAMILY.
!
!    Input/output, integer NU(M+1), part of the code for the object.
!    The use and meaning of NU depend on TASK and FAMILY.
!
!    Input/output, integer EDGE(M), part of the code for the object.
!    The use and meaning of EDGE depend on TASK and FAMILY.
!
!    Input/output, integer M, the length of the code.  For
!    families 1 through 5, M is always equal to N.
!    For family 6, M varies, but is never greater than N.
!    For family 7, M = N+K-1.
!
!    Input/output, logical NEWONE, is only needed for TASK 1.
!    It is set to FALSE by the user when asking for the first object.
!    The routine returns the first object, and returns NEWONE as TRUE.
!    If the user now calls again (leaving NEWONE TRUE), the next object
!    will be returned, and so on.  When there are no more objects to
!    return, the routine sets NEWONE to FALSE.
!
!    Input/output, integer RANK, the rank of an object.
!    For task 1, 2, and 4, RANK is an output quantity.
!    For task 3, RANK is an input quantity, and must be
!    specified by the user.  For TASK 5, RANK is not needed.
!
!    Input/output, integer SEED, a seed for the random number generator.
!    SEED is only needed for TASK 4.
!
!    Output, integer B(N_MAX,K_MAX), is computed by the program as needed.
!    B(I,J) is the number of objects of "size" J in the family of size I.
!    The user may examine B, but should not change its entries or dimension.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) b(n_max,k_max)
  integer ( kind = 4 ) b1
  integer ( kind = 4 ) b2
  integer ( kind = 4 ) bnew
  real ( kind = 8 ) d
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) edge(n+k-1)
  integer ( kind = 4 ) family
  integer, save :: family_save = 0
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer, save :: k_max_save = 0
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mu(n+k)
  integer ( kind = 4 ) n1
  logical newone
  integer ( kind = 4 ) nu(n+k)
  integer, save :: n_max_save = 0
  integer ( kind = 4 ) phi
  integer ( kind = 4 ) psi
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rankp
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2
  integer ( kind = 4 ) task
  integer ( kind = 4 ) xnew
!
!  Recompute the B matrix if the family or maximum dimensions have changed.
!
  if ( family_save /= family .or. &
       k_max_save < k_max .or. &
       n_max_save < n_max ) then

    n_max_save = n_max
    k_max_save = k_max
    family_save = family

    do n1 = 1, n_max
      do k1 = 1, k_max

        b(n1,k1) = phi ( n1, k1, family ) &
                 * bnew ( n_max, k_max, n1, k1, family, 1, b ) &
                 + psi ( n1, k1, family ) &
                 * bnew ( n_max, k_max, n1, k1, family, 2, b )

      end do
    end do

  end if

  if ( task == 1 ) then

    if ( .not. newone ) then

      newone = .true.
      rank = 0
      rankp = rank
      j = 1
      m = j
      mu(j) = n
      nu(j) = k

    else

      rank = rank + 1
      j = m
!
!  Starting at J = M, seek EDGE(J) that is movable.
!  But if we reach J = 0, there are no more, and we
!  have exhausted the family list.
!
      do

        t1 = phi ( mu(j), nu(j), family )
        t2 = psi ( mu(j), nu(j), family )

        if ( edge(j) + 1 < t1 + t2 ) then
          exit
        end if

        j = j - 1

        if ( j == 0 ) then
          newone = .false.
          return
        end if

      end do

      edge(j) = edge(j) + 1
      l = j + 1

      if ( edge(j) == t1 ) then
        nu(l) = nu(l) - 1
        mu(l) = xnew ( mu(l-1), nu(l-1), family, 2 )
      end if

      rankp = 0
      m = l

    end if

  else if ( task == 2 ) then

    rank = 0

    do j = 1, m - 1

      if ( nu(j+1) == nu(j) ) then
        rank = rank + edge(j) &
          * bnew ( n_max, k_max, mu(j), nu(j), family, 1, b )
      else
        rank = rank &
          + phi ( mu(j), nu(j), family ) &
          * bnew ( n_max, k_max, mu(j), nu(j), family, 1, b ) &
          + ( edge(j) - phi ( mu(j), nu(j), family) ) &
          * bnew ( n_max, k_max, mu(j), nu(j), family, 2, b )
      end if

    end do

  else if ( task == 3 ) then

    rankp = rank
    j = 1
    m = j
    mu(j) = n
    nu(j) = k

  else if ( task == 4 ) then

    d = r8_uniform_01 ( seed )

    rank = int ( d * real ( b(n,k), kind = 8 ) )
    rankp = rank
    j = 1
    m = j
    mu(j) = n
    nu(j) = k
!
!  For task 5, we simply return the current matrix B.
!
  else if ( task == 5 ) then

    return

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SELECT - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of TASK = ', task
    write ( *, '(a)' ) '  Legal values satisfy 1 <= TASK <= 5.'
    stop

  end if

  do

    t1 = phi ( mu(m), nu(m), family )

    t2 = psi ( mu(m), nu(m), family )

    if ( t1 + t2 == 0 ) then
      m = m - 1
      exit
    end if

    b1 = bnew ( n_max, k_max, mu(m), nu(m), family, 1, b )

    if ( rankp < t1 * b1 ) then

      edge(m) = rankp / b1
      rankp = rankp - edge(m) * b1
      mu(m+1) = xnew ( mu(m), nu(m), family, 1 )
      nu(m+1) = nu(m)

    else

      rankp = rankp - b1 * t1
      b2 = bnew ( n_max, k_max, mu(m), nu(m), family, 2, b )
      edge(m) = t1 + rankp / b2
      rankp = rankp - ( edge(m) - t1 ) * b2
      mu(m+1) = xnew ( mu(m), nu(m), family, 2 )
      nu(m+1) = nu(m) - 1

    end if

    m = m + 1

  end do

  return
end
subroutine subset_to_string ( k, subset, s )

!*****************************************************************************80
!
!! SUBSET_TO_STRING writes a character representation of a subset.
!
!  Example:
!
!    K = 3
!    SUBSET = (/ 3, 4, 12 /)
!
!    S = '{3,4,12}'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer K, the size of the object.  K = 0 is allowed.
!
!    Input, integer SUBSET(K), the indices of the items in the subset.
!
!    Output, character ( len = * ) S, a representation of the subset.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character ( len = * ) s
  integer ( kind = 4 ) subset(k)

  s = ' '

  if ( k < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUBSET_TO_STRING - Fatal error!'
    write ( *, '(a)' ) '  K < 0.'
    stop
  end if

  i = 0

  i = i + 1
  s(i:i) = '{'

  do j = 1, k

    if ( 1 < j ) then
      i = i + 1
      s(i:i)=','
    end if

    i = i + 1
    call i4_to_s_left ( subset(j), s(i:) )

    i = len_trim ( s )

  end do

  i = i + 1
  s(i:i) = '}'

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
function xnew ( n, k, family, mgo )

!*****************************************************************************80
!
!! XNEW returns an index for the recursive enumeration of a family.
!
!  Discussion:
!
!    The recursive enumerations have the form:
!
!      E(N,K) = E(N1,K1) + E(N2,K2)
!
!    where, typically, but not always:
!
!      N1 = N - 1
!      N2 = N - 1
!      K1 = K
!      K2 = K - 1
!
!    This routine returns the correct value of K1 or K2 for the
!    given family.
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
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    pages 110-113.
!
!  Parameters:
!
!    Input, integer N, the size of the family.
!
!    Input, integer K, the size of the object.
!
!    Input, integer FAMILY, specifies the combinatorial family.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional space
!       over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!    Input, integer MGO, is 1 to request index K1, or 2 to request
!    index K2.
!
!    Output, integer XNEW, the value of the requested index.
!
  implicit none

  integer ( kind = 4 ) family
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mgo
  integer ( kind = 4 ) n
  integer ( kind = 4 ) xnew

  if ( family == 6 .and. mgo == 1 ) then
    xnew = n - k
  else if ( family == 7 .and. mgo == 2 ) then
    xnew = n
  else
    xnew = n - 1
  end if

  return
end
