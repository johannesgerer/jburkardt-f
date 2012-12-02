program main

!*****************************************************************************80
!
!! MAIN is the main program for COMBO_PRB.
!
!  Discussion:
!
!    COMBO_PRB calls the COMBO tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMBO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the COMBO library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )

  call test20 ( )
  call test21 ( )
  call test22 ( )
  call test23 ( )
  call test24 ( )
  call test25 ( )
  call test26 ( )
  call test27 ( )
  call test28 ( )
  call test29 ( )

  call test30 ( )
  call test31 ( )
  call test32 ( )
  call test33 ( )
  call test34 ( )
  call test35 ( )
  call test36 ( )
  call test37 ( )
  call test38 ( )
  call test39 ( )

  call test40 ( )
  call test41 ( )
  call test42 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMBO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' )  ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BAL_SEQ_ENUM, BAL_SEQ_RANK, BAL_SEQ_SUCCESSOR, BAL_SEQ_UNRANK.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) nseq
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(2*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Balanced sequences:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BAL_SEQ_ENUM enumerates,'
  write ( *, '(a)' ) '  BAL_SEQ_RANK ranks,'
  write ( *, '(a)' ) '  BAL_SEQ_SUCCESSOR lists,'
  write ( *, '(a)' ) '  BAL_SEQ_UNRANK unranks.'
!
!  Enumerate.
!
  call bal_seq_enum ( n, nseq )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  the number of balanced sequences is ', nseq
  write ( *, '(a)' ) ' '
!
!  List.
!
  rank = -1

  do

    rank_old = rank

    call bal_seq_successor ( n, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(4x,i3,2x,10i2)' ) rank, t(1:2*n)

  end do
!
!  Unrank.
!
  rank = nseq / 2

  call bal_seq_unrank ( rank, n, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The element of rank ', rank
  write ( *, '(a)' ) ' '
  write ( *, '(4x,10i2)' ) t(1:2*n)
!
!  Rank.
!
  call bal_seq_rank ( n, t, rank )

  call i4vec_print ( 2*n, t, '  Element to be ranked:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Computed rank: ', rank

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests BAL_SEQ_TO_TABLEAU, TABLEAU_TO_BAL_SEQ.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(2*n)
  integer ( kind = 4 ) tab(2,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  BAL_SEQ_TO_TABLEAU converts a balanced'
  write ( *, '(a)' ) '  sequence to a tableau;'
  write ( *, '(a)' ) '  TABLEAU_TO_BAL_SEQ converts a tableau'
  write ( *, '(a)' ) '  to a balanced sequence.'
!
!  Pick a "random" balanced sequence.
!
  rank = 7

  call bal_seq_unrank ( rank, n, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "Random" balanced sequence:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,8i2)' ) t(1:2*n)
!
!  Convert to a tableau.
!
  call bal_seq_to_tableau ( n, t, tab )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Corresponding tableau'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,4i2)' ) tab(1,1:n)
  write ( *, '(4x,4i2)' ) tab(2,1:n)
!
!  Convert to a balanced sequence.
!
  call tableau_to_bal_seq ( n, tab, t )

  call i4vec_print ( 2*n, t, '  Corresponding balanced sequence:' )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests BELL_NUMBERS.
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
  implicit none

  integer ( kind = 4 ), allocatable :: b(:)
  integer ( kind = 4 ) bn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  BELL_NUMBERS computes Bell numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        BELL(N)    BELL_NUMBERS(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bell_values ( n_data, n, bn )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( b(0:n) )

    call bell_numbers ( n, b )

    write ( *, '(2x,i8,2x,i12,2x,i12)' ) n, bn, b(n)

    deallocate ( b )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests BINOMIAL.
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
  implicit none

  integer ( kind = 4 ) binomial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  BINOMIAL computes binomial coefficients.'

  do i = -1, 5
    do j = -1, 5
      write ( *, '(3i8)' ) i, j, binomial ( i, j )
    end do
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests CYCLE_TO_PERM, PERM_TO_CYCLE.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ) i
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) index(n)
  integer ( kind = 4 ) ncycle
  integer ( kind = 4 ) nperm
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  CYCLE_TO_PERM converts a permutation from'
  write ( *, '(a)' ) '  cycle to array form;'
  write ( *, '(a)' ) '  PERM_TO_CYCLE converts a permutation from'
  write ( *, '(a)' ) '  array to cycle form.'
!
!  Enumerate.
!
  call perm_enum ( n, nperm )
!
!  Choose a "random" permutation.
!
  rank = nperm / 2

  call perm_lex_unrank ( rank, n, p )

  call perm_print ( n, p, '  "Random" permutation:' )
!
!  Convert the permutation to cycle form.
!
  call perm_to_cycle ( n, p, ncycle, t, index )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Corresponding cycle form:'
  write ( *, '(a,i8)' ) '  Number of cycles is ', ncycle
  write ( *, '(a)' ) ' '
  jlo = 0
  do i = 1, ncycle
    write ( *, '(4x,20i4)' ) t(jlo+1:jlo+index(i))
    jlo = jlo + index(i)
  end do
!
!  Convert the set partition back to an RGF.
!
  call cycle_to_perm ( n, ncycle, t, index, p )

  call perm_print ( n, p, '  Corresponding permutation:' )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests DIST_ENUM and DIST_NEXT.
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
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) idist
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) num_dist
  integer ( kind = 4 ) q(k)

  m = 5
  more = .false.

  call dist_enum ( k, m, num_dist )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  For a distribution of M indistinguishable'
  write ( *, '(a)' ) '  objects among K distinguishable slots:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DIST_ENUM enumerates them;'
  write ( *, '(a)' ) '  DIST_NEXT produces the "next" one.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Number of:'
  write ( *, '(a,i8)' ) '    indistinguishable objects = ', m
  write ( *, '(a,i8)' ) '    distinguishable slots =     ', k
  write ( *, '(a,i8)' ) '    distributions is            ', num_dist
  write ( *, '(a)' ) ' '

  idist = 0

  do

    call dist_next ( k, m, q, more )

    if ( .not. more ) then
      exit
    end if

    idist = idist + 1
    write ( *, '(4x,6i5)' ) idist, q(1:k)

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests I4_FACTORIAL.
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
  implicit none

  integer ( kind = 4 ) i4_factorial
  integer ( kind = 4 ) fx
  integer ( kind = 4 ) fx2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  I4_FACTORIAL evaluates the factorial function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       FACTORIAL(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call i4_factorial_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x <= 0.0D+00 ) then
      cycle
    end if

    fx2 = i4_factorial ( x )

    write ( *, '(i4,2i12)' ) x, fx, fx2

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests GRAY_CODE_*.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) ngray
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Gray codes:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GRAY_CODE_ENUM enumerates,'
  write ( *, '(a)' ) '  GRAY_CODE_RANK ranks,'
  write ( *, '(a)' ) '  GRAY_CODE_SUCCESSOR lists,'
  write ( *, '(a)' ) '  GRAY_CODE_UNRANK unranks.'
!
!  Enumerate.
!
  call gray_code_enum ( n, ngray )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  the number of Gray code elements is ', ngray
  write ( *, '(a)' ) ' '
!
!  List
!
  rank = -1

  do

    rank_old = rank

    call gray_code_successor ( n, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(4x,6i5)' ) rank, t(1:n)

  end do
!
!  Unrank.
!
  rank = ngray / 2

  call gray_code_unrank ( rank, n, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The element of rank ', rank
  write ( *, '(a)' ) ' '
  write ( *, '(4x,6i5)' ) t(1:n)
!
!  Rank.
!
  call gray_code_rank ( n, t, rank )

  call i4vec_print ( n, t, '  Element to be ranked:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Computed rank: ', rank

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests I4VEC_SEARCH_BINARY_A and I4VEC_SORT_INSERT_A.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) index

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  Integer vectors:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I4VEC_SORT_INSERT_A ascending sorts;'
  write ( *, '(a)' ) '  I4VEC_SEARCH_BINARY_A searches a ascending sorted vector.'

  a(1:n) = (/ 6, 7, 1, 0, 4, 3, 2, 1, 5, 8 /)

  call i4vec_print ( n, a, '  Before ascending sort:' )

  call i4vec_sort_insert_a ( n, a )

  call i4vec_print ( n, a, '  After ascending sort:' )

  b = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Now search for an instance of the value ', b

  call i4vec_search_binary_a ( n, a, b, index )

  write ( *, '(a)' ) ' '
  if ( index == 0 ) then
    write ( *, '(a)' ) '  The value does not occur.'
  else
    write ( *, '(a,i8)' ) '  The value occurs at index = ', index
  end if

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests I4VEC_SEARCH_BINARY_D and I4VEC_SORT_INSERT_D.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) index

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  Integer vectors:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I4VEC_SORT_INSERT_D descending sorts;'
  write ( *, '(a)' ) '  I4VEC_SEARCH_BINARY_D searches a descending '
  write ( *, '(a)' ) ' sorted vector.'

  a(1:n) = (/ 6, 7, 1, 0, 4, 3, 2, 1, 5, 8 /)

  call i4vec_print ( n, a, '  Before descending sort:' )

  call i4vec_sort_insert_d ( n, a )

  call i4vec_print ( n, a, '  After descending sort:' )

  b = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Now search for an instance of the value ', b

  call i4vec_search_binary_d ( n, a, b, index )

  write ( *, '(a)' ) ' '
  if ( index == 0 ) then
    write ( *, '(a)' ) '  The value does not occur.'
  else
    write ( *, '(a,i8)' ) '  The value occurs at index = ', index
  end if

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests KNAPSACK_REORDER and KNAPSACK_01.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  real    ( kind = 8 ) :: mass
  real    ( kind = 8 ) :: mass_limit = 26.0
  real    ( kind = 8 ), dimension ( n ) :: p = (/ &
    24.0, 13.0, 23.0, 15.0, 16.0 /)
  real    ( kind = 8 ) :: profit
  real    ( kind = 8 ), dimension ( n ) :: w = (/ &
    12.0,  7.0, 11.0,  8.0,  9.0 /)
  real    ( kind = 8 ), dimension ( n ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  KNAPSACK_REORDER reorders the knapsack data.'
  write ( *, '(a)' ) '  KNAPSACK_01 solves the 0/1 knapsack problem.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Object, Profit, Mass, "Profit Density"'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f7.3,2x,f7.3,2x,f7.3)' ) i, p(i), w(i), p(i)/w(i)
  end do

  call knapsack_reorder ( n, p, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After reordering by Profit Density:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Object, Profit, Mass, "Profit Density"'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,f7.3,2x,f7.3,2x,f7.3)' ) i, p(i), w(i), p(i) / w(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,f7.3)' ) '  Total mass restriction is ', mass_limit

  call knapsack_01 ( n, mass_limit, p, w, x, mass, profit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Object, Density, Choice, Profit, Mass'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,f7.3,f7.3,2f7.3)' ) i, p(i)/w(i), x(i), &
      x(i) * p(i), x(i) * w(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,2f7.3)' ) '  Total:            ', profit, mass

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests KNAPSACK_REORDER and KNAPSACK_RATIONAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  real    ( kind = 8 ) :: mass
  real    ( kind = 8 ) :: mass_limit = 26.0
  real    ( kind = 8 ), dimension ( n ) :: p = (/ &
    24.0, 13.0, 23.0, 15.0, 16.0 /)
  real    ( kind = 8 ) :: profit
  real    ( kind = 8 ), dimension ( n ) :: w = (/ &
    12.0,  7.0, 11.0,  8.0,  9.0 /)
  real    ( kind = 8 ), dimension ( n ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  KNAPSACK_REORDER reorders the knapsack data.'
  write ( *, '(a)' ) '  KNAPSACK_RATIONAL solves the rational knapsack problem.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Object, Profit, Mass, "Profit Density"'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,3f7.3)' ) i, p(i), w(i), p(i) / w(i)
  end do

  call knapsack_reorder ( n, p, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After reordering by Profit Density:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Object, Profit, Mass, "Profit Density"'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,3f7.3)' ) i, p(i), w(i), p(i) / w(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,f7.3)' ) '  Total mass restriction is ', mass_limit

  call knapsack_rational ( n, mass_limit, p, w, x, mass, profit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Object, Density, Choice, Profit, Mass'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,f7.3,f7.3,2f7.3)' ) i, p(i) / w(i), x(i), &
      x(i) * p(i), x(i) * w(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,2f7.3)' ) '  Total:            ', profit, mass

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests KSUBSET_COLEX_RANK, _SUCCESSOR, _UNRANK, _ENUM.
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
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nksub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(k)

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  K-subsets of an N set,'
  write ( *, '(a)' ) '  using the colexicographic ordering:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KSUBSET_COLEX_RANK ranks,'
  write ( *, '(a)' ) '  KSUBSET_COLEX_SUCCESSOR lists,'
  write ( *, '(a)' ) '  KSUBSET_COLEX_UNRANK unranks.'
  write ( *, '(a)' ) '  KSUBSET_ENUM enumerates,'
!
!  Enumerate.
!
  call ksubset_enum ( k, n, nksub )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  the number of K subsets is ', nksub
  write ( *, '(a)' ) ' '
!
!  List
!
  rank = -1

  do

    rank_old = rank

    call ksubset_colex_successor ( k, n, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(4x,6i5)' ) rank, t(1:k)

  end do
!
!  Unrank.
!
  rank = nksub / 2

  call ksubset_colex_unrank ( rank, k, n, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The element of rank ', rank
  write ( *, '(a)' ) ' '
  write ( *, '(4x,6i5)' ) t(1:k)
!
!  Rank.
!
  call ksubset_colex_rank ( k, n, t, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The rank of the element:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,6i5)' ) t(1:k)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  is computed as ', rank

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests KSUBSET_ENUM, _LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK.
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
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nksub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(k)

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  K-subsets of an N set,'
  write ( *, '(a)' ) '  using the lexicographic ordering:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KSUBSET_ENUM enumerates,'
  write ( *, '(a)' ) '  KSUBSET_LEX_RANK ranks,'
  write ( *, '(a)' ) '  KSUBSET_LEX_SUCCESSOR lists,'
  write ( *, '(a)' ) '  KSUBSET_LEX_UNRANK unranks.'
!
!  Enumerate.
!
  call ksubset_enum ( k, n, nksub )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  the number of K subsets is ', nksub
  write ( *, '(a)' ) ' '
!
!  List
!
  rank = -1

  do

    rank_old = rank

    call ksubset_lex_successor ( k, n, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(6i5)' ) rank, t(1:k)

  end do
!
!  Unrank.
!
  rank = nksub / 2

  call ksubset_lex_unrank ( rank, k, n, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The element of rank ', rank
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) t(1:k)
!
!  Rank.
!
  call ksubset_lex_rank ( k, n, t, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The rank of the element:'
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) t(1:k)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  is computed as ', rank

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests KSUBSET_ENUM, _REVDOOR_RANK, _REVDOOR_SUCCESSOR, _REVDOOR_UNRANK.
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
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nksub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(k)

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  K-subsets of an N set,'
  write ( *, '(a)' ) '  using the revolving door ordering:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KSUBSET_ENUM enumerates,'
  write ( *, '(a)' ) '  KSUBSET_REVDOOR_RANK ranks,'
  write ( *, '(a)' ) '  KSUBSET_REVDOOR_SUCCESSOR lists,'
  write ( *, '(a)' ) '  KSUBSET_REVDOOR_UNRANK unranks.'
!
!  Enumerate.
!
  call ksubset_enum ( k, n, nksub )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  the number of K subsets is ', nksub
  write ( *, '(a)' ) ' '
!
!  List
!
  rank = -1

  do

    rank_old = rank

    call ksubset_revdoor_successor ( k, n, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(6i5)' ) rank, t(1:k)

  end do
!
!  Unrank.
!
  rank = nksub / 2

  call ksubset_revdoor_unrank ( rank, k, n, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The element of rank ', rank
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) t(1:k)
!
!  Rank.
!
  call ksubset_revdoor_rank ( k, n, t, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The rank of the element:'
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) t(1:k)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  is computed as ', rank

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests MARRIAGE.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) fiancee(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) next(n)
  integer ( kind = 4 ) prefer(n,n)
  integer ( kind = 4 ) rank(n,n)
!
!  PREFER(M,W) is the index of women W on man M's list.
!
  prefer(1:5,1:5) = reshape ( (/ &
    2, 1, 2, 1, 5, &
    5, 2, 3, 3, 3, &
    1, 3, 5, 2, 2, &
    3, 4, 4, 4, 1, &
    4, 5, 1, 5, 4  &
    /), (/ 5, 5 /) )
!
!  RANK(W,M) is the index of man M on woman W's list.
!
  rank(1:5,1:5) = reshape ( (/ &
    2, 4, 1, 4, 5, &
    4, 3, 3, 2, 2, &
    5, 5, 4, 1, 3, &
    3, 1, 2, 3, 1, &
    1, 2, 5, 5, 4  &
   /), (/ 5, 5 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  MARRIAGE arranges a set of stable marriages'
  write ( *, '(a)' ) '  given a set of preferences.'

  call marriage ( n, prefer, rank, fiancee, next )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Man, Wife''s rank, Wife'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(3i8)' ) i, next(i), prefer(i,next(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Woman, Husband''s rank, Husband'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(3i8)' ) i, rank(i,fiancee(i)), fiancee(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Correct result:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  M:W 1  2  3  4  5'
  write ( *, '(a)' ) '   1  +  .  .  .  .'
  write ( *, '(a)' ) '   2  .  .  .  +  .'
  write ( *, '(a)' ) '   3  .  .  .  .  +'
  write ( *, '(a)' ) '   4  .  .  +  .  .'
  write ( *, '(a)' ) '   5  .  +  .  .  .'

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests MOUNTAIN.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) mxy
  integer ( kind = 4 ) row(0:2*n)
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  MOUNTAIN computes mountain numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Y    MXY'
  write ( *, '(a)' ) ' '

  do y = 0, n
    do x = 0, 2*n
      call mountain ( n, x, y, mxy )
      row(x) = mxy
    end do
    write ( *, '(2x,i2,3x,11i4)' ) y, row(0:2*n)
  end do

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests NPART_ENUM, _RSF_LEX_RANK, _RSF_LEX_SUCCESSOR, _RSF_LEX_UNRANK.
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
  implicit none

  integer ( kind = 4 ), parameter :: npart = 3

  integer ( kind = 4 ) n
  integer ( kind = 4 ) npartitions
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(npart)

  n = 12

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  Partitions of N with NPART parts'
  write ( *, '(a)' ) '  in reverse standard form:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NPART_ENUM enumerates,'
  write ( *, '(a)' ) '  NPART_RSF_LEX_RANK ranks,'
  write ( *, '(a)' ) '  NPART_RSF_LEX_SUCCESSOR lists;'
  write ( *, '(a)' ) '  NPART_RSF_LEX_UNRANK unranks.'
!
!  Enumerate.
!
  call npart_enum ( n, npart, npartitions )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  and NPART = ', npart
  write ( *, '(a,i8)' ) '  the number of partitions is ', npartitions
  write ( *, '(a)' ) ' '
!
!  List.
!
  rank = -1

  do

    rank_old = rank

    call npart_rsf_lex_successor ( n, npart, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(6i5)' ) rank, t(1:npart)

  end do
!
!  Unrank.
!
  rank = npartitions / 3

  call npart_rsf_lex_unrank ( rank, n, npart, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The element of rank ', rank
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) t(1:npart)
!
!  Rank.
!
  call npart_rsf_lex_rank ( n, npart, t, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The rank of the element:'
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) t(1:npart)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  is computed as ', rank

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests NPART_RSF_LEX_RANDOM;
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
  implicit none

  integer ( kind = 4 ), parameter :: npart = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) t(npart)

  n = 12
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  Partitions of N with NPART parts'
  write ( *, '(a)' ) '  in reverse standard form:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NPART_RSF_LEX_RANDOM produces random examples.'

  do i = 1, 10

    call npart_rsf_lex_random ( n, npart, seed, t )

    write ( *, '(6i5)' ) t(1:npart)

  end do

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests NPART_ENUM and NPART_SF_SUCCESSOR;
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
  implicit none

  integer ( kind = 4 ), parameter :: npart = 3

  integer ( kind = 4 ) n
  integer ( kind = 4 ) npartitions
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(npart)

  n = 12

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  Partitions of N with NPART parts'
  write ( *, '(a)' ) '  in standard form:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NPART_ENUM enumerates,'
  write ( *, '(a)' ) '  NPART_SF_LEX_SUCCESSOR lists.'
!
!  Enumerate.
!
  call npart_enum ( n, npart, npartitions )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  and NPART = ', npart
  write ( *, '(a,i8)' ) '  the number of partitions is ', npartitions
  write ( *, '(a)' ) ' '
!
!  List.
!
  rank = -1

  do

    rank_old = rank

    call npart_sf_lex_successor ( n, npart, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(6i5)' ) rank, t(1:npart)

  end do

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests NPART_TABLE and PART_TABLE.
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
  implicit none

  integer ( kind = 4 ), parameter :: maxn = 10
  integer ( kind = 4 ), parameter :: maxpart = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(0:maxn,0:maxpart)
  integer ( kind = 4 ) p2(0:maxn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  NPART_TABLE tabulates partitions'
  write ( *, '(a)' ) '  of N with NPART parts;'
  write ( *, '(a)' ) '  PART_TABLE tabulates partitions of N.'

  call npart_table ( maxn, maxpart, maxn, p )

  call part_table ( maxn, p2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I P(I)  P(I,0) P(I,1) P(I,2) P(I,3) P(I,4) P(I,5)'
  write ( *, '(a)' ) ' '

  do i = 0, maxn
    write ( *, '(11i5)' ) i, p2(i), p(i,0:maxpart)
  end do

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests PART_ENUM and PART_SUCCESSOR.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 8

  integer ( kind = 4 ) npart
  integer ( kind = 4 ) npartitions
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  PART_SUCCESSOR produces partitions of N,'
  write ( *, '(a)' ) '  PART_ENUM enumerates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Partitions of N = ', n
!
!  Enumerate.
!
  call part_enum ( n, npartitions )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  the number of partitions is ', npartitions
  write ( *, '(a)' ) ' '
!
!  List.
!
  rank = -1

  do

    rank_old = rank

    call part_successor ( n, npart, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(2x,i2,3x,10i3)' ) rank, t(1:npart)

  end do

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests PART_SUCCESSOR and PART_SF_CONJUGATE.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 8

  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) npartb
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  PART_SUCCESSOR produces partitions of N,'
  write ( *, '(a)' ) '  PART_SF_CONJUGATE produces the conjugate of a partition.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Partitions of N = ', n
!
!  List.
!
  rank = -1

  do

    rank_old = rank

    call part_successor ( n, npart, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(2x,i2,4x,10i3)' ) rank, t(1:npart)
    call part_sf_conjugate ( n, npart, t, npartb, b )
    write ( *, '(2x,a4,2x,10i3)' ) 'Con:', b(1:npartb)

  end do

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests PART_SF_MAJORIZE.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 8

  integer ( kind = 4 ), parameter, dimension ( n ) :: a = (/ 2, 2, 2, 1, 1, 0, 0, 0 /)
  integer ( kind = 4 ), parameter, dimension ( n ) :: b = (/ 3, 1, 1, 1, 1, 1, 0, 0 /)
  integer ( kind = 4 ), parameter, dimension ( n ) :: c = (/ 2, 2, 1, 1, 1, 1, 0, 0 /)
  integer ( kind = 4 ), parameter :: nparta = 5
  integer ( kind = 4 ), parameter :: npartb = 6
  integer ( kind = 4 ), parameter :: npartc = 6
  integer ( kind = 4 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  PART_SF_MAJORIZE determines if one partition'
  write ( *, '(a)' ) '  majorizes another.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Partitions of N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a,2x,10i3)' ) 'A:', a(1:nparta)
  write ( *, '(2x,a,2x,10i3)' ) 'B:', b(1:npartb)
  write ( *, '(2x,a,2x,10i3)' ) 'C:', c(1:npartc)

  call part_sf_majorize ( n, nparta, a, npartb, b, result )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  A compare B: ', result
  call part_sf_majorize ( n, npartb, b, npartc, c, result )
  write ( *, '(a,i8)' ) '  B compare C: ', result
  call part_sf_majorize ( n, npartc, c, nparta, a, result )
  write ( *, '(a,i8)' ) '  C compare A: ', result
  call part_sf_majorize ( n, npartc, c, npartc, c, result )
  write ( *, '(a,i8)' ) '  C compare C: ', result

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests PARTITION_GREEDY.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) sums(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  PARTITION_GREEDY partitions an integer vector into'
  write ( *, '(a)' ) '  two subsets with nearly equal sum.'

  a = (/ 2, 10, 3, 8, 5, 7, 9, 5, 3, 2 /)

  call partition_greedy ( n, a, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Data set #1 partitioned:'
  write ( *, '(a)' ) ' '

  sums(1:2) = 0

  do i = 1, n

    if ( indx(i) == 1 ) then
      sums(1) = sums(1) + a(i)
      write ( *, '(2x,i6)' ) a(i)
    else
      write ( *, '(2x,6x,i6)' ) a(i)
      sums(2) = sums(2) + a(i)
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sums:'
  write ( *, '(a)' ) ' '
  write ( *, '(2i6)' ) sums(1), sums(2)

  a = (/ 771, 121, 281, 854, 885, 734, 486, 1003, 83, 62 /)

  call partition_greedy ( n, a, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set #2 partitioned:'
  write ( *, '(a)' ) ' '

  sums(1:2) = 0

  do i = 1, n

    if ( indx(i) == 1 ) then
      sums(1) = sums(1) + a(i)
      write ( *, '(i6)' ) a(i)
    else
      write ( *, '(6x,i6)' ) a(i)
      sums(2) = sums(2) + a(i)
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sums:'
  write ( *, '(a)' ) ' '
  write ( *, '(2i6)' ) sums(1), sums(2)

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26 tests PARTN_ENUM, PARTN_SUCCESSOR and PART_SF_CONJUGATE.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) nmax
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) npart2
  integer ( kind = 4 ) npartitions
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  Partitions of N with maximum element NMAX:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PARTN_SUCCESSOR lists;'
  write ( *, '(a)' ) '  PARTN_ENUM enumerates.'

  nmax = 4
!
!  Enumerate.
!
  call partn_enum ( n, nmax, npartitions )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  and NMAX = ', nmax
  write ( *, '(a,i8)' ) '  the number of partitions is ', npartitions
  write ( *, '(a)' ) ' '
!
!  List.
!
  rank = -1

  do

    rank_old = rank

    call partn_successor ( n, nmax, npart, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(2x,i2,3x,15i3)' ) rank, t(1:npart)

  end do
!
!  List conjugates.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat, but list RSF conjugated partitions.'
  write ( *, '(a)' ) ' '
  rank = -1

  do

    rank_old = rank

    call partn_successor ( n, nmax, npart, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    call part_sf_conjugate ( n, npart, t, npart2, b )
    call i4vec_reverse ( npart2, b )
    write ( *, '(2x,i2,3x,15i3)' ) rank, b(1:npart2)

  end do

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests PERM_INV and PERM_MUL.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) nperm
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) q(n)
  integer ( kind = 4 ) r(n)
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  Permutations of the integers:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PERM_INV computes an inverse permutation,'
  write ( *, '(a)' ) '  PERM_MUL multiplies two permutations.'
!
!  Enumerate.
!
  call perm_enum ( n, nperm )
!
!  Unrank.
!
  rank = nperm / 2

  call perm_lex_unrank ( rank, n, p )

  call perm_print ( n, p, '  The permutation P is ' )
!
!  Invert.
!
  call perm_inv ( n, p, q )

  call perm_print ( n, q, '  The inverse permutation Q is ' )
!
!  Multiply.
!
  call perm_mul ( n, p, q, r )

  call perm_print ( n, r, '  The product R = P * Q is ' )

  return
end
subroutine test28 ( )

!*****************************************************************************80
!
!! TEST28 tests PERM_ENUM, _LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) nperm
  integer ( kind = 4 ) pi(n)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  Permutations of the integers,'
  write ( *, '(a)' ) '  using the lexicographic ordering:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PERM_ENUM enumerates,'
  write ( *, '(a)' ) '  PERM_LEX_RANK ranks,'
  write ( *, '(a)' ) '  PERM_LEX_SUCCESSOR lists,'
  write ( *, '(a)' ) '  PERM_LEX_UNRANK unranks.'
!
!  Enumerate.
!
  call perm_enum ( n, nperm )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  the number of permutations is ', nperm
  write ( *, '(a)' ) ' '
!
!  List
!
  rank = -1

  do

    rank_old = rank

    call perm_lex_successor ( n, pi, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(6i5)' ) rank, pi(1:n)

  end do
!
!  Unrank.
!
  rank = nperm / 2

  call perm_lex_unrank ( rank, n, pi )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The element of rank ', rank

  call perm_print ( n, pi, ' ' )
!
!  Rank.
!
  call perm_lex_rank ( n, pi, rank )

  call perm_print ( n, pi, '  The rank of the element:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) 'is computed as ', rank

  return
end
subroutine test29 ( )

!*****************************************************************************80
!
!! TEST29 tests PERM_TJ_ENUM, _TJ_RANK, _TJ_SUCCESSOR, _TJ_UNRANK.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) nperm
  integer ( kind = 4 ) pi(n)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  Permutations of the integers'
  write ( *, '(a)' ) '  using the Trotter-Johnson ordering:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PERM_ENUM enumerates,'
  write ( *, '(a)' ) '  PERM_TJ_RANK ranks,'
  write ( *, '(a)' ) '  PERM_TJ_SUCCESSOR lists,'
  write ( *, '(a)' ) '  PERM_TJ_UNRANK unranks.'
!
!  Enumerate.
!
  call perm_enum ( n, nperm )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  the number of permutations is ', nperm
  write ( *, '(a)' ) ' '
!
!  List
!
  rank = -1

  do

    rank_old = rank

    call perm_tj_successor ( n, pi, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(6i5)' ) rank, pi(1:n)

  end do
!
!  Unrank.
!
  rank = nperm / 2

  call perm_tj_unrank ( rank, n, pi )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The element of rank ', rank

  call perm_print ( n, pi, ' ' )
!
!  Rank.
!
  call perm_tj_rank ( n, pi, rank )

  call perm_print ( n, pi, '  The rank of the element:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  is computed as ', rank

  return
end
subroutine test30 ( )

!*****************************************************************************80
!
!! TEST30 tests PRUEFER_ENUM, PRUEFER_RANK, PRUEFER_SUCCESSOR, PRUEFER_UNRANK.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) ncode
  integer ( kind = 4 ) p(n-2)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  Pruefer codes:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PRUEFER_ENUM enumerates,'
  write ( *, '(a)' ) '  PRUEFER_RANK ranks,'
  write ( *, '(a)' ) '  PRUEFER_SUCCESSOR lists,'
  write ( *, '(a)' ) '  PRUEFER_UNRANK unranks.'
!
!  Enumerate.
!
  call pruefer_enum ( n, ncode )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  the number of Pruefer codes is ', ncode
  write ( *, '(a)' ) ' '
!
!  List
!
  rank = -1

  do

    rank_old = rank

    call pruefer_successor ( n, p, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(6i5)' ) rank, p(1:n-2)

  end do
!
!  Unrank.
!
  rank = ncode / 2

  call pruefer_unrank ( rank, n, p )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The element of rank ', rank
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) p(1:n-2)
!
!  Rank.
!
  call pruefer_rank ( n, p, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The rank of the element:'
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) p(1:n-2)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  is computed as ', rank

  return
end
subroutine test31 ( )

!*****************************************************************************80
!
!! TEST31 tests PRUEFER_TO_TREE and TREE_TO_PRUEFER.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i4_hi
  integer ( kind = 4 ) i4_lo
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n-2)
  integer ( kind = 4 ) pruefer_num
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) t(2,n-1)
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  PRUEFER_TO_TREE converts a Pruefer code to'
  write ( *, '(a)' ) '  a tree;'
  write ( *, '(a)' ) '  TREE_TO_PRUEFER converts a tree to a Pruefer'
  write ( *, '(a)' ) '  code.'

  call pruefer_enum ( n, pruefer_num )

  i4_lo = 0
  i4_hi = pruefer_num - 1

  do test = 1, test_num
!
!  Pick a "random" Pruefer code.
!
    rank = i4_uniform ( i4_lo, i4_hi, seed )

    call pruefer_unrank ( rank, n, p )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Random Pruefer code of rank ', rank
    write ( *, '(6i5)' ) p(1:n-2)
!
!  Convert the Pruefer code to a tree.
!
    call pruefer_to_tree ( n, p, t )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Edge list for the corresponding tree:'
    write ( *, '(a)' ) ' '
    do j = 1, n - 1
      write ( *, '(6i5)' ) j, t(1:2,j)
    end do
!
!  Convert the tree to a Pruefer code.
!
    call tree_to_pruefer ( n, t, p )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Corresponding Pruefer code:'
    write ( *, '(6i5)' ) p(1:n-2)

  end do

  return
end
subroutine test32 ( )

!*****************************************************************************80
!
!! TEST32 tests QUEENS and BACKTRACK.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: maxstack = n * n

  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) istack(maxstack)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nstack

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32'
  write ( *, '(a)' ) '  QUEENS produces nonattacking queens'
  write ( *, '(a)' ) '  on a chessboard.'
  write ( *, '(a)' ) '  BACKTRACK supervises a backtrack search.'
  write ( *, '(a)' ) ' '

  indx = 0

  do

    call backtrack ( n, iarray, indx, k, nstack, istack, maxstack )

    if ( indx == 1 ) then

      write ( *, '(19i4)' ) iarray(1:n)

    else if ( indx == 2 ) then

      call queens ( n, iarray, k, nstack, istack, maxstack )

    else

      exit

    end if

  end do

  return
end
subroutine test33 ( )

!*****************************************************************************80
!
!! TEST33 tests RGF_G_TABLE.
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
  implicit none

  integer ( kind = 4 ), parameter :: MMAX = 8

  integer ( kind = 4 ) d(0:MMAX,0:MMAX)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m

  m = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33'
  write ( *, '(a)' ) '  RGF_G_TABLE tabulates generalized restricted'
  write ( *, '(a)' ) '  growth functions.'
  write ( *, '(a)' ) ' '

  call rgf_g_table ( m, MMAX, d )

  do i = 0, m
    write ( *, '(7i6)' ) d(i,0:m-i)
  end do

  return
end
subroutine test34 ( )

!*****************************************************************************80
!
!! TEST34 tests RGF_ENUM, RGF_RANK, RGF_SUCCESSOR, RGF_UNRANK.
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
  implicit none

  integer ( kind = 4 ), parameter :: m = 4

  integer ( kind = 4 ) f(m)
  integer ( kind = 4 ) nrgf
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST34'
  write ( *, '(a)' ) '  Restricted growth functions:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RGF_ENUM enumerates,'
  write ( *, '(a)' ) '  RGF_RANK ranks,'
  write ( *, '(a)' ) '  RGF_SUCCESSOR lists;'
  write ( *, '(a)' ) '  RGF_UNRANK unranks.'
!
!  Enumerate.
!
  call rgf_enum ( m, nrgf )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For M = ', m
  write ( *, '(a,i8)' ) '  the number of RGF''s is ', nrgf
  write ( *, '(a)' ) ' '
!
!  List.
!
  rank = -1

  do

    rank_old = rank

    call rgf_successor ( m, f, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(6i5)' ) rank, f(1:m)

  end do
!
!  Unrank.
!
  rank = nrgf / 2

  call rgf_unrank ( rank, m, f )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The element of rank ', rank
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) f(1:m)
!
!  Rank.
!
  call rgf_rank ( m, f, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The rank of the element:'
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) f(1:m)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  is computed as ', rank
  return
end
subroutine test35 ( )

!*****************************************************************************80
!
!! TEST35 tests RGF_TO_SETPART and SETPART_TO_RGF.
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
  implicit none

  integer ( kind = 4 ), parameter :: m = 8

  integer ( kind = 4 ) i
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) f(m)
  integer ( kind = 4 ) index(m)
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) s(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST35'
  write ( *, '(a)' ) '  RGF_TO_SETPART converts a balanced'
  write ( *, '(a)' ) '  sequence to a restricted growth function;'
  write ( *, '(a)' ) '  SETPART_TO_RGF converts a restricted growth'
  write ( *, '(a)' ) '  function to a balanced sequence.'
!
!  Choose a "random" RGF.
!
  rank = 7
  call rgf_unrank ( rank, m, f )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "Random" restricted growth function:'
  write ( *, '(a)' ) ' '
  write ( *, '(8i2)' ) f(1:m)
!
!  Convert the RGF to a set partition.
!
  call rgf_to_setpart ( m, f, nsub, s, index )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Corresponding set partition'
  write ( *, '(a)' ) ' '
  jlo = 1
  do i = 1, nsub
    write ( *, '(8i4)' ) s(jlo:index(i))
    jlo = index(i) + 1
  end do
!
!  Convert the set partition back to an RGF.
!
  call setpart_to_rgf ( m, nsub, s, index, f )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Corresponding restricted growth function:'
  write ( *, '(a)' ) ' '
  write ( *, '(8i2)' ) f(1:m)

  return
end
subroutine test36 ( )

!*****************************************************************************80
!
!! TEST36 tests SETPART_ENUM.
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
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) npart

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST36'
  write ( *, '(a)' ) '  Set partitions:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SETPART_ENUM enumerates,'
  write ( *, '(a)' ) ' '
!
!  Enumerate.
!
  do n = 1, 6
    call setpart_enum ( n, npart )
    write ( *, '(i6,i6)' ) n, npart
  end do

  return
end
subroutine test37 ( )

!*****************************************************************************80
!
!! TEST37 tests STIRLING_NUMBERS1.
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
  implicit none

  integer ( kind = 4 ), parameter :: maxm = 6
  integer ( kind = 4 ), parameter :: maxn = 6

  integer ( kind = 4 ) i
  integer ( kind = 4 ) s(0:maxm,0:maxn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST37'
  write ( *, '(a)' ) '  STIRLING_NUMBERS1 computes a table of Stirling'
  write ( *, '(a)' ) '  numbers of the first kind.'

  call stirling_numbers1 ( maxm, maxn, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I S(I,0) S(I,1) S(I,2) S(I,3) S(I,4) S(I,5)'
  write ( *, '(a)' ) ' '

  do i = 0, maxm
    write ( *, '(11i5)' ) i, s(i,0:maxn)
  end do

  return
end
subroutine test38 ( )

!*****************************************************************************80
!
!! TEST38 tests STIRLING_NUMBERS2.
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
  implicit none

  integer ( kind = 4 ), parameter :: maxm = 6
  integer ( kind = 4 ), parameter :: maxn = 6

  integer ( kind = 4 ) i
  integer ( kind = 4 ) s(0:maxm,0:maxn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST38'
  write ( *, '(a)' ) '  STIRLING_NUMBERS2 computes a table of Stirling'
  write ( *, '(a)' ) '  numbers of the second kind.'

  call stirling_numbers2 ( maxm, maxn, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I S(I,0) S(I,1) S(I,2) S(I,3) S(I,4) S(I,5)'
  write ( *, '(a)' ) ' '

  do i = 0, maxm
    write ( *, '(11i5)' ) i, s(i,0:maxn)
  end do

  return
end
subroutine test39 ( )

!*****************************************************************************80
!
!! TEST39 tests SUBSET_COLEX_RANK, _COLEX_SUCCESSOR, _COLEX_UNRANK, _ENUM.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST39'
  write ( *, '(a)' ) '  All subsets of a set,'
  write ( *, '(a)' ) '  using the colexicographic ordering:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SUBSET_COLEX_RANK ranks,'
  write ( *, '(a)' ) '  SUBSET_COLEX_SUCCESSOR lists,'
  write ( *, '(a)' ) '  SUBSET_COLEX_UNRANK unranks.'
  write ( *, '(a)' ) '  SUBSET_ENUM enumerates.'
!
!  Enumerate.
!
  call subset_enum ( n, nsub )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  the number of subsets is ', nsub
  write ( *, '(a)' ) ' '
!
!  List
!
  rank = -1

  do

    rank_old = rank

    call subset_colex_successor ( n, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(6i5)' ) rank, t(1:n)

  end do
!
!  Unrank.
!
  rank = nsub / 3

  call subset_colex_unrank ( rank, n, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The element of rank ', rank
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) t(1:n)
!
!  Rank.
!
  call subset_colex_rank ( n, t, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The rank of the element:'
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) t(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  is computed as ', rank

  return
end
subroutine test40 ( )

!*****************************************************************************80
!
!! TEST40 tests SUBSET_ENUM, _LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST40'
  write ( *, '(a)' ) '  All subsets of a set,'
  write ( *, '(a)' ) '  using the lexicographic ordering:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SUBSET_ENUM enumerates,'
  write ( *, '(a)' ) '  SUBSET_LEX_RANK ranks,'
  write ( *, '(a)' ) '  SUBSET_LEX_SUCCESSOR lists,'
  write ( *, '(a)' ) '  SUBSET_LEX_UNRANK unranks.'
!
!  Enumerate.
!
  call subset_enum ( n, nsub )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  the number of subsets is ', nsub
  write ( *, '(a)' ) ' '
!
!  List
!
  rank = -1

  do

    rank_old = rank

    call subset_lex_successor ( n, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(6i5)' ) rank, t(1:n)

  end do
!
!  Unrank.
!
  rank = nsub / 3

  call subset_lex_unrank ( rank, n, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The element of rank ', rank
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) t(1:n)
!
!  Rank.
!
  call subset_lex_rank ( n, t, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The rank of the element:'
  write ( *, '(a)' ) ' '
  write ( *, '(6i5)' ) t(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  is computed as ', rank

  return
end
subroutine test41 ( )

!*****************************************************************************80
!
!! TEST41 tests SUBSETSUM_SWAP.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) index(n)
  integer ( kind = 4 ) sum_achieved
  integer ( kind = 4 ) sum_desired

  sum_desired = 17

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST41'
  write ( *, '(a)' ) '  SUBSETSUM_SWAP seeks a solution of the subset'
  write ( *, '(a)' ) '  sum problem using pair swapping.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The desired sum is ', sum_desired

  a(1:7) = (/ 12, 8, 11, 30, 8, 3, 7 /)

  call subsetsum_swap ( n, a, sum_desired, index, sum_achieved )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    A(I), INDEX(I)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,2i5)' ) a(i), index(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The achieved sum is ', sum_achieved

  return
end
subroutine test42 ( )

!*****************************************************************************80
!
!! TEST42 tests TREE_ENUM, TREE_RANK, TREE_SUCCESSOR, TREE_UNRANK.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_old
  integer ( kind = 4 ) t(2,n-1)
  integer ( kind = 4 ) tree_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST42'
  write ( *, '(a)' ) '  Trees:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TREE_ENUM enumerates,'
  write ( *, '(a)' ) '  TREE_RANK ranks,'
  write ( *, '(a)' ) '  TREE_SUCCESSOR lists,'
  write ( *, '(a)' ) '  TREE_UNRANK unranks.'
!
!  Enumerate.
!
  call tree_enum ( n, tree_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  the number of trees is ', tree_num
  write ( *, '(a)' ) ' '
!
!  List
!
  rank = -1

  do

    rank_old = rank

    call tree_successor ( n, t, rank )

    if ( rank <= rank_old ) then
      exit
    end if

    write ( *, '(i5,2x,5i5)' ) rank, t(1,1:n-1)
    write ( *, '(5x,2x,5i5)' )       t(2,1:n-1)

  end do
!
!  Unrank.
!
  rank = tree_num / 2

  call tree_unrank ( rank, n, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The element of rank ', rank
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5i5)' ) t(1,1:n-1)
  write ( *, '(2x,5i5)' ) t(2,1:n-1)
!
!  Rank.
!
  call tree_rank ( n, t, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The rank of the element:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5i5)' ) t(1,1:n-1)
  write ( *, '(2x,5i5)' ) t(2,1:n-1)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  is computed as ', rank

  return
end
