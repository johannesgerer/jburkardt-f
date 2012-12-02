program main

!*****************************************************************************80
!
!! MAIN is the main program for SELECT_PRB.
!
!  Discussion:
!
!    SELECT_PRB carries out tests of the SELECT algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) family

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SELECT_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SELECT library.'
!
!  Try each family.
!
  do family = 1, 7

    write ( *, '(a)' ) ' '
    write ( *, '(a,i1)' ) '  FAMILY = ', family

    if ( family == 1 ) then
      write ( *, '(a)' ) '    K subsets of an N set.'
    else if ( family == 2 ) then
      write ( *, '(a)' ) '    Partitions of N objects into K classes.'
    else if ( family == 3 ) then
      write ( *, '(a)' ) '    Permutations of N objects with K cycles.'
    else if ( family == 4 ) then
      write ( *, '(a)' ) '    Vector subspaces of dimension K over'
      write ( *, '(a)' ) '    N-space over the field of order Q (where Q = 2 ).'
    else if ( family == 5 ) then
      write ( *, '(a)' ) '    Permutations of N letters with K increasing runs.'
    else if ( family == 6 ) then
      write ( *, '(a)' ) '    Partitions of N whose largest part is K.'
    else if ( family == 7 ) then
      write ( *, '(a)' ) '    Compositions of N into K parts.'
    end if

    call test01 ( family )

  end do
!
!  Termimate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SELECT_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( family )

!*****************************************************************************80
!
!! TEST01 does the tests for a particular FAMILY.
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
!  Parameters:
!
!    Input, integer FAMILY, indicates the family to be considered.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional
!       space over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K increasing runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!  Local parameters:
!
!    Local, integer K_MAX, the largest value of K, the object size,
!    to consider.
!
!    Local, integer N_MAX, the largest value of N, the family size,
!    to consider.
!
  implicit none

  integer ( kind = 4 ), parameter :: k_max = 10
  integer ( kind = 4 ), parameter :: n_max = 10

  integer ( kind = 4 ) b(n_max,k_max)
  integer ( kind = 4 ) family
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) task
!
!  Try each task.
!
  do task = 1, 5

    if ( task == 1 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  TASK = 1'
      write ( *, '(a)' ) '    Present each object of the family.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Maximum family size N_MAX = ', n_max
      write ( *, '(a,i12)' ) '  Maximum object size K_MAX = ', k_max

      call family_enumerate ( family, n_max, k_max, b )

      do n = 3, n_max-1, 3
        do k = 1, n, (n/3)
          call test03 ( family, task, n_max, k_max, b, n, k )
        end do
      end do

    else if ( task == 2 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  TASK = 2'
      write ( *, '(a)' ) '    Rank a given object of the family.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Maximum family size N_MAX = ', n_max
      write ( *, '(a,i12)' ) '  Maximum object size K_MAX = ', k_max

      call family_enumerate ( family, n_max, k_max, b )

      do n = 3, n_max-1, 3
        do k = 1, n, (n/3)
          call test03 ( family, task, n_max, k_max, b, n, k )
        end do
      end do

    else if ( task == 3 ) then

      do n = 3, n_max-1, 3
        do k = 1, n, (n/3)
          call test_family_unrank ( family, n_max, k_max, b, n, k )
        end do
      end do

    else if ( task == 4 ) then

      seed = 123456789

      do n = 3, n_max-1, 3
        do k = 1, n, (n/3)
          call test_family_sample ( family, n_max, k_max, n, k, seed )
        end do
      end do

    else if ( task == 5 ) then

      call test_family_enumerate ( family, n_max, k_max )

    end if

  end do

  return
end
subroutine test_family_enumerate ( family, n_max, k_max )

!*****************************************************************************80
!
!! TEST_FAMILY_ENUMERATE tests FAMILY_ENUMERATE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FAMILY, indicates the family to be considered.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional
!       space over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K increasing runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!    Input, integer N_MAX, the largest value of N, the family size,
!    to consider.
!
!    Input, integer K_MAX, the largest value of K, the object size,
!    to consider.
!
  implicit none

  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) b(n_max,k_max)
  integer ( kind = 4 ) family

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_FAMILY_ENUMERATE'
  write ( *, '(a)' ) '  FAMILY_ENUMERATE enumerates the objects of the family.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Maximum family size N_MAX = ', n_max
  write ( *, '(a,i12)' ) '  Maximum object size K_MAX = ', k_max

  call family_enumerate ( family, n_max, k_max, b )

  call i4mat_print ( n_max, k_max, b, '  Enumeration Matrix' )

  return
end
subroutine test_family_sample ( family, n_max, k_max, n, k, seed  )

!*****************************************************************************80
!
!! TEST_FAMILY_SAMPLE tests FAMILY_SAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FAMILY, indicates the family to be considered.
!    1, K subsets of an N set.
!    2, Partitions of N objects into K classes.
!    3, Permutations of N objects with K cycles.
!    4, Vector subspaces of dimension K over N-dimensional
!       space over the field of order Q (where Q = 2 ).
!    5, Permutations of N letters with K increasing runs.
!    6, Partitions of N whose largest part is K.
!    7, Compositions of N into K parts.
!
!    Input, integer N_MAX, the largest value of N, the family size,
!    to consider.
!
!    Input, integer K_MAX, the largest value of K, the object size,
!    to consider.
!
!    Input, integer N, the family size.
!
!    Input, integer K, the object size.
!
!    Input/output, integer SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) b(n_max,k_max)
  integer ( kind = 4 ) composition(k_max)
  integer ( kind = 4 ) edge(n_max+k_max)
  integer ( kind = 4 ) family
  character ( len = 80 ) form
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mu(n_max+k_max+1)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nu(n_max+k_max+1)
  integer ( kind = 4 ) partition_num(n_max)
  integer ( kind = 4 ) partition_set(n_max)
  integer ( kind = 4 ) perm_cycle(n_max)
  integer ( kind = 4 ) perm_runs(n_max)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) seed
  character ( len = 80 ) string
  integer ( kind = 4 ) subset(k_max)
  integer ( kind = 4 ) subspace(k_max,n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_FAMILY_SAMPLE'
  write ( *, '(a)' ) '  FAMILY_SAMPLE randomly samples the objects of the family.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Maximum family size N_MAX = ', n_max
  write ( *, '(a,i12)' ) '  Maximum object size K_MAX = ', k_max
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Family size N = ', n
  write ( *, '(a,i12)' ) '  Object size K = ', k
  write ( *, '(a,i12)' ) '  Initial SEED =  ', seed

  call family_enumerate ( family, n_max, k_max, b )

  write ( *, '(a)' ) ' '

  if ( family == 1 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Subset      M'
  else if ( family == 2 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      SetPart     M'
  else if ( family == 3 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Perm_Cycle  M'
  else if ( family == 4 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Subspace    M'
  else if ( family == 5 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Perm_Runs   M'
  else if ( family == 6 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      NumPart     M'
  else if ( family == 7 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Composition M'
  end if

  write ( *, '(a)' ) ' '

  do j = 1, 5

    call family_sample ( family, n_max, k_max, b, n, k, seed, rank, m, mu, &
      nu, edge )

    if ( family == 1 ) then

      call code_to_subset ( n, k, mu, nu, edge, m, subset )

      call subset_to_string ( k, subset, string )

      write ( form, '(a,i2,a,i2,a,i2,a)' ) &
        '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

      write ( *, form ) &
        rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

    else if ( family == 2 ) then

      call code_to_partition_set ( n, k, mu, nu, edge, m, partition_set )

      call partition_set_to_string ( n, k, partition_set, string )

      write ( form, '(a,i2,a,i2,a,i2,a)' ) &
        '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

      write ( *, form ) &
        rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

    else if ( family == 3 ) then

      call code_to_perm_cycle ( n, k, mu, nu, edge, m, perm_cycle )

      call perm_cycle_to_string ( n, k, perm_cycle, string )

      write ( form, '(a,i2,a,i2,a,i2,a)' ) &
        '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

      write ( *, form ) &
        rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

    else if ( family == 4 ) then

      call code_to_subspace ( n, k, mu, nu, edge, m, subspace )

      write ( *, '(2x,i6,2x,5i1,2x,6i1,2x,6i1,2x,5i2,2x,i2)' ) &
        rank, edge(1:n), mu(1:n+1), nu(1:n+1), subspace(1,1:n), m
      do i = 2, k
        write ( *, '(31x,5i2)' ) subspace(i,1:n)
      end do
      write ( *, '(a)' ) ' '

    else if ( family == 5 ) then

      call code_to_perm_runs ( n, k, mu, nu, edge, m, perm_runs )

      write ( *, '(2x,i6,2x,5i1,2x,6i1,2x,6i1,2x,5i2,2x,i2)' ) &
        rank, edge(1:n), mu(1:n+1), nu(1:n+1), perm_runs(1:n), m

    else if ( family == 6 ) then

      call code_to_partition_num ( n, k, mu, nu, edge, m, partition_num )

      call partition_num_to_string ( n, k, m, partition_num, string )

      write ( form, '(a,i2,a,i2,a,i2,a)' ) &
        '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

      write ( *, form ) &
        rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

    else if ( family == 7 ) then

      call code_to_composition ( n, k, mu, nu, edge, m, composition )

      call composition_to_string ( n, k, composition, string )

      write ( *, '(2x,i6,2x,7i1,2x,8i1,2x,8i1,2x,a,2x,i2)' ) &
        rank, edge(1:7), mu(1:7+1), nu(1:7+1), trim ( string ), m

    end if

  end do

  return
end
subroutine test_family_unrank ( family, n_max, k_max, b, n, k )

!*****************************************************************************80
!
!! TEST_FAMILY_UNRANK unranks an object from FAMILY with N <= N_MAX, K <= K_MAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n_max
  integer ( kind = 4 ) k_max

  integer ( kind = 4 ) b(n_max,k_max)
  integer ( kind = 4 ) composition(k_max)
  integer ( kind = 4 ) edge(n_max+k_max)
  integer ( kind = 4 ) family
  character ( len = 80 ) form
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mu(n_max+k_max+1)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nu(n_max+k_max+1)
  integer ( kind = 4 ) partition_num(n_max)
  integer ( kind = 4 ) partition_set(n_max)
  integer ( kind = 4 ) perm_cycle(n_max)
  integer ( kind = 4 ) perm_runs(n_max)
  integer ( kind = 4 ) rank
  integer, save :: seed = 123456789
  character ( len = 80 ) string
  integer ( kind = 4 ) subset(k_max)
  integer ( kind = 4 ) subspace(k_max,n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_FAMILY_UNRANK'
  write ( *, '(a)' ) '  FAMILY_UNRANK returns the object of given rank.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Maximum family size N_MAX = ', n_max
  write ( *, '(a,i12)' ) '  Maximum object size K_MAX = ', k_max
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Family size N =             ', n
  write ( *, '(a,i12)' ) '  Object size K =             ', k

  call family_enumerate ( family, n_max, k_max, b )

  write ( *, '(a)' ) ' '
  if ( family == 1 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Subset      M'
  else if ( family == 2 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      SetPart     M'
  else if ( family == 3 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Perm_Cycle  M'
  else if ( family == 4 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Subspace    M'
  else if ( family == 5 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Perm_Runs   M'
  else if ( family == 6 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      NumPart     M'
  else if ( family == 7 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Composition M'
  end if

  write ( *, '(a)' ) ' '

  do i = 1, 5

    rank = i4_uniform ( 0, b(n,k)-1, seed )

    call family_unrank ( family, n_max, k_max, b, n, k, rank, m, mu, nu, &
      edge )

    if ( family == 1 ) then

      call code_to_subset ( n, k, mu, nu, edge, m, subset )

      call subset_to_string ( k, subset, string )

      write ( form, '(a,i2,a,i2,a,i2,a)' ) &
        '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

      write ( *, form ) &
        rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

    else if ( family == 2 ) then

      call code_to_partition_set ( n, k, mu, nu, edge, m, partition_set )

      call partition_set_to_string ( n, k, partition_set, string )

      write ( form, '(a,i2,a,i2,a,i2,a)' ) &
        '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

      write ( *, form ) &
        rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

    else if ( family == 3 ) then

      call code_to_perm_cycle ( n, k, mu, nu, edge, m, perm_cycle )

      call perm_cycle_to_string ( n, k, perm_cycle, string )

      write ( form, '(a,i2,a,i2,a,i2,a)' ) &
        '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

      write ( *, form ) &
        rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

    else if ( family == 4 ) then

      call code_to_subspace ( n, k, mu, nu, edge, m, subspace )

      write ( *, '(2x,i6,2x,5i1,2x,6i1,2x,6i1,2x,5i2,2x,i2)' ) &
        rank, edge(1:n), mu(1:n+1), nu(1:n+1), subspace(1,1:n), m
      do ii = 2, k
        write ( *, '(31x,5i2)' ) subspace(ii,1:n)
      end do
      write ( *, '(a)' ) ' '

    else if ( family == 5 ) then

      call code_to_perm_runs ( n, k, mu, nu, edge, m, perm_runs )

      write ( *, '(2x,i6,2x,5i1,2x,6i1,2x,6i1,2x,5i2,2x,i2)' ) &
        rank, edge(1:n), mu(1:n+1), nu(1:n+1), perm_runs(1:n), m

    else if ( family == 6 ) then

      call code_to_partition_num ( n, k, mu, nu, edge, m, partition_num )

      call partition_num_to_string ( n, k, m, partition_num, string )

      write ( form, '(a,i2,a,i2,a,i2,a)' ) &
        '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

      write ( *, form ) &
        rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

    else if ( family == 7 ) then

      call code_to_composition ( n, k, mu, nu, edge, m, composition )

      call composition_to_string ( n, k, composition, string )

      write ( *, '(2x,i6,2x,7i1,2x,8i1,2x,8i1,2x,a,2x,i2)' ) &
        rank, edge(1:7), mu(1:7+1), nu(1:7+1), trim ( string ), m

    end if

  end do

  return
end
subroutine test03 ( family, task, n_max, k_max, b, n, k )

!*****************************************************************************80
!
!! TEST03 carries out TASK for FAMILY with N <= N_MAX, K <= K_MAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n_max
  integer ( kind = 4 ) k_max

  integer ( kind = 4 ) b(n_max,k_max)
  integer ( kind = 4 ) composition(k_max)
  integer ( kind = 4 ) edge(n_max+k_max)
  integer ( kind = 4 ) family
  character ( len = 80 ) form
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mu(n_max+k_max+1)
  integer ( kind = 4 ) n
  logical newone
  integer ( kind = 4 ) nu(n_max+k_max+1)
  integer ( kind = 4 ) partition_num(n_max)
  integer ( kind = 4 ) partition_set(n_max)
  integer ( kind = 4 ) perm_cycle(n_max)
  integer ( kind = 4 ) perm_runs(n_max)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank2
  integer, save :: seed = 123456789
  character ( len = 80 ) string
  integer ( kind = 4 ) subset(k_max)
  integer ( kind = 4 ) subspace(k_max,n_max)
  integer ( kind = 4 ) task

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Family size N =             ', n
  write ( *, '(a,i12)' ) '  Object size K =             ', k
  if ( task == 4 ) then
    write ( *, '(a,i12)' ) '  Random number SEED = ', seed
  end if
  write ( *, '(a)' ) ' '

  if ( family == 1 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Subset      M'
  else if ( family == 2 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      SetPart     M'
  else if ( family == 3 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Perm_Cycle  M'
  else if ( family == 4 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Subspace    M'
  else if ( family == 5 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Perm_Runs   M'
  else if ( family == 6 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      NumPart     M'
  else if ( family == 7 ) then
    write ( *, '(a)' ) '    Rank  Edge   Mu      Nu      Composition M'
  end if

  write ( *, '(a)' ) ' '

  if ( task == 1 ) then

    newone = .false.

    do

      call family_next ( family, n_max, k_max, b, n, k, newone, rank, m, mu, &
        nu, edge )

      if ( .not. newone ) then
        write ( *, '(a)' ) '  (Last Element)'
        exit
      end if

      if ( 10 <= rank ) then
        write ( *, '(a)' ) '  (That''s enough!)'
        exit
      end if

      if ( family == 1 ) then

        call code_to_subset ( n, k, mu, nu, edge, m, subset )

        call subset_to_string ( k, subset, string )

        write ( form, '(a,i2,a,i2,a,i2,a)' ) &
          '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

        write ( *, form ) &
          rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

      else if ( family == 2 ) then

        call code_to_partition_set ( n, k, mu, nu, edge, m, partition_set )

        call partition_set_to_string ( n, k, partition_set, string )

        write ( form, '(a,i2,a,i2,a,i2,a)' ) &
          '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

        write ( *, form ) &
          rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

      else if ( family == 3 ) then

        call code_to_perm_cycle ( n, k, mu, nu, edge, m, perm_cycle )

        call perm_cycle_to_string ( n, k, perm_cycle, string )

        write ( form, '(a,i2,a,i2,a,i2,a)' ) &
          '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

        write ( *, form ) &
          rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

      else if ( family == 4 ) then

        call code_to_subspace ( n, k, mu, nu, edge, m, subspace )

        write ( *, '(2x,i6,2x,5i1,2x,6i1,2x,6i1,2x,5i2,2x,i2)' ) &
          rank, edge(1:n), mu(1:n+1), nu(1:n+1), subspace(1,1:n), m

        do i = 2, k
          write ( *, '(31x,5i2)' ) subspace(i,1:n)
        end do
        write ( *, '(a)' ) ' '

        if ( 10 <= rank ) then
          exit
        end if

      else if ( family == 5 ) then

        call code_to_perm_runs ( n, k, mu, nu, edge, m, perm_runs )

        write ( *, '(2x,i6,2x,5i1,2x,6i1,2x,6i1,2x,5i2,2x,i2)' ) &
          rank, edge(1:n), mu(1:n+1), nu(1:n+1), perm_runs(1:n), m

      else if ( family == 6 ) then

        call code_to_partition_num ( n, k, mu, nu, edge, m, partition_num )

        call partition_num_to_string ( n, k, m, partition_num, string )

        write ( form, '(a,i2,a,i2,a,i2,a)' ) &
          '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

        write ( *, form ) &
          rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

      else if ( family == 7 ) then

        call code_to_composition ( n, k, mu, nu, edge, m, composition )

        call composition_to_string ( n, k, composition, string )

        write ( *, '(2x,i6,2x,7i1,2x,8i1,2x,8i1,2x,a,2x,i2)' ) &
          rank, edge(1:7), mu(1:7+1), nu(1:7+1), trim ( string ), m

      end if

    end do
!
!  Rank a particular object.
!  Here, we first select one at random.
!
  else if ( task == 2 ) then

    call family_sample ( family, n_max, k_max, b, n, k, seed, rank, m, mu, &
      nu, edge )

    rank2 = rank

    if ( family == 1 ) then

      call code_to_subset ( n, k, mu, nu, edge, m, subset )

      call subset_to_string ( k, subset, string )

      write ( form, '(a,i2,a,i2,a,i2,a)' ) &
        '(2x,a,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

      write ( *, form ) &
        '??????', edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

    else if ( family == 2 ) then

      call code_to_partition_set ( n, k, mu, nu, edge, m, partition_set )

      call partition_set_to_string ( n, k, partition_set, string )

        write ( form, '(a,i2,a,i2,a,i2,a)' ) &
          '(2x,a,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

        write ( *, form ) &
          '??????', edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

    else if ( family == 3 ) then

      call code_to_perm_cycle ( n, k, mu, nu, edge, m, perm_cycle )

      call perm_cycle_to_string ( n, k, perm_cycle, string )

      write ( form, '(a,i2,a,i2,a,i2,a)' ) &
        '(2x,a,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

      write ( *, form ) &
        '??????', edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

    else if ( family == 4 ) then

      call code_to_subspace ( n, k, mu, nu, edge, m, subspace )

      write ( *, '(2x,''???????'',2x,5i1,2x,6i1,2x,6i1,2x,5i2,2x,i2)' ) &
        edge(1:n), mu(1:n+1), nu(1:n+1), subspace(1,1:n), m

      do i = 2, k
        write ( *, '(31x,5i2)' ) subspace(i,1:n)
      end do

    else if ( family == 5 ) then

      call code_to_perm_runs ( n, k, mu, nu, edge, m, perm_runs )

      write ( *, '(2x,''??????'',2x,5i1,2x,6i1,2x,6i1,2x,5i2,2x,i2)' ) &
        edge(1:n), mu(1:n+1), nu(1:n+1), perm_runs(1:n), m

    else if ( family == 6 ) then

      call code_to_partition_num ( n, k, mu, nu, edge, m, partition_num )

      call partition_num_to_string ( n, k, m, partition_num, string )

      write ( form, '(a,i2,a,i2,a,i2,a)' ) &
        '(2x,a,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

      write ( *, form ) &
        '??????', edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

    else if ( family == 7 ) then

      call code_to_composition ( n, k, mu, nu, edge, m, composition )

      call composition_to_string ( n, k, composition, string )

      write ( *, '(2x,a4,2x,7i1,2x,8i1,2x,8i1,2x,a,2x,i2)' ) &
        '??????', edge(1:7), mu(1:7+1), nu(1:7+1), trim ( string ), m

    end if

    call family_rank ( family, n_max, k_max, b, n, k, m, mu, nu, edge, &
      rank )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  Computed rank is ', rank
    write ( *, '(a,i12)' ) '  Correct rank is  ', rank2

  else if ( task == 3 ) then

    do i = 1, 5

      rank = i4_uniform ( 0, b(n,k)-1, seed )

      call family_unrank ( family, n_max, k_max, b, n, k, rank, m, mu, nu, &
        edge )

      if ( family == 1 ) then

        call code_to_subset ( n, k, mu, nu, edge, m, subset )

        call subset_to_string ( k, subset, string )

        write ( form, '(a,i2,a,i2,a,i2,a)' ) &
          '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

        write ( *, form ) &
          rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

      else if ( family == 2 ) then

        call code_to_partition_set ( n, k, mu, nu, edge, m, partition_set )

        call partition_set_to_string ( n, k, partition_set, string )

        write ( form, '(a,i2,a,i2,a,i2,a)' ) &
          '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

        write ( *, form ) &
          rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

      else if ( family == 3 ) then

        call code_to_perm_cycle ( n, k, mu, nu, edge, m, perm_cycle )

        call perm_cycle_to_string ( n, k, perm_cycle, string )

        write ( form, '(a,i2,a,i2,a,i2,a)' ) &
          '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

        write ( *, form ) &
          rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

      else if ( family == 4 ) then

        call code_to_subspace ( n, k, mu, nu, edge, m, subspace )

        write ( *, '(2x,i6,2x,5i1,2x,6i1,2x,6i1,2x,5i2,2x,i2)' ) &
          rank, edge(1:n), mu(1:n+1), nu(1:n+1), subspace(1,1:n), m
        do ii = 2, k
          write ( *, '(31x,5i2)' ) subspace(ii,1:n)
        end do
        write ( *, '(a)' ) ' '

      else if ( family == 5 ) then

        call code_to_perm_runs ( n, k, mu, nu, edge, m, perm_runs )

        write ( *, '(2x,i6,2x,5i1,2x,6i1,2x,6i1,2x,5i2,2x,i2)' ) &
          rank, edge(1:n), mu(1:n+1), nu(1:n+1), perm_runs(1:n), m

      else if ( family == 6 ) then

        call code_to_partition_num ( n, k, mu, nu, edge, m, partition_num )

        call partition_num_to_string ( n, k, m, partition_num, string )

        write ( form, '(a,i2,a,i2,a,i2,a)' ) &
          '(2x,i6,2x,', n, 'i1,2x,', n+1, 'i1,2x,', n+1, 'i1,2x,a,2x,i2)'

        write ( *, form ) &
          rank, edge(1:n), mu(1:n+1), nu(1:n+1), trim ( string ), m

      else if ( family == 7 ) then

        call code_to_composition ( n, k, mu, nu, edge, m, composition )

        call composition_to_string ( n, k, composition, string )

        write ( *, '(2x,i6,2x,7i1,2x,8i1,2x,8i1,2x,a,2x,i2)' ) &
          rank, edge(1:7), mu(1:7+1), nu(1:7+1), trim ( string ), m

      end if

    end do

  end if

  return
end
