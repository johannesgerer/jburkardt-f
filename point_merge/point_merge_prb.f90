program main

!*****************************************************************************80
!
!! MAIN is the main program for POINT_MERGE_PRB.
!
!  Discussion:
!
!    Compare correctness of the codes.
!
!    Compare speed of the codes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_unique
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) tol

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POINT_MERGE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the POINT_MERGE library.'
!
!  TEST01 gives me some confidence that, at least for zero-tolerance,
!  the radial approach is accurate, as compared to the "Unique count"
!  (which cannot be extended to a tolerance version in multidimensions)
!  and the "Tol Unique Count", which is an O(N^2) algorithm.
!
  m = 3
  n = 10
  n_unique = 7
  seed = 123456789
  call test01 ( m, n, n_unique, seed )

  m = 4
  n = 20
  n_unique = 11
  seed = 987654321
  call test01 ( m, n, n_unique, seed )
!
!  In TEST02, I want to compute the same data, but with "blurred"
!  duplicates, and a tolerance version of the radial approach,
!  compared to "Tol Unique Count".
!
  m = 3
  n = 10
  n_unique = 7
  tol = 0.00001D+00
  seed = 123456789
  call test02 ( m, n, n_unique, tol, seed )

  m = 4
  n = 20
  n_unique = 11
  tol = 0.00001D+00
  seed = 987654321
  call test02 ( m, n, n_unique, tol, seed )
!
!  In TEST03, I want to measure the time required for a sequence
!  of increasingly hard problems.
!
  m = 3
  n = 100
  n_unique = n / 2
  tol = 0.00001D+00
  seed = 123456789
  call test03 ( m, n, n_unique, tol, seed )

  m = 3
  n = 1000
  n_unique = n / 2
  tol = 0.00001D+00
  seed = 123456789
  call test03 ( m, n, n_unique, tol, seed )

  m = 3
  n = 10000
  n_unique = n / 2
  tol = 0.00001D+00
  seed = 123456789
  call test03 ( m, n, n_unique, tol, seed )

  if ( .false. ) then
    m = 3
    n = 100000
    n_unique = n / 2
    tol = 0.00001D+00
    seed = 123456789
    call test03 ( m, n, n_unique, tol, seed )
  end if
!
!  In TEST04, repeat TEST02, but now compute the index vector.
!
  m = 3
  n = 10
  n_unique = 7
  tol = 0.00001D+00
  seed = 123456789
  call test04 ( m, n, n_unique, tol, seed )

  m = 4
  n = 20
  n_unique = 11
  tol = 0.00001D+00
  seed = 987654321
  call test04 ( m, n, n_unique, tol, seed )
!
!  In TEST05, I want to measure the time required for a sequence
!  of increasingly hard problems.
!
  m = 3
  n = 100
  n_unique = n / 2
  tol = 0.00001D+00
  seed = 123456789
  call test05 ( m, n, n_unique, tol, seed )

  m = 3
  n = 1000
  n_unique = n / 2
  tol = 0.00001D+00
  seed = 123456789
  call test05 ( m, n, n_unique, tol, seed )

  m = 3
  n = 10000
  n_unique = n / 2
  tol = 0.00001D+00
  seed = 123456789
  call test05 ( m, n, n_unique, tol, seed )

  if ( .false. ) then
    m = 3
    n = 100000
    n_unique = n / 2
    tol = 0.00001D+00
    seed = 123456789
    call test05 ( m, n, n_unique, tol, seed )
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POINT_MERGE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( m, n, n_unique, seed )

!*****************************************************************************80
!
!! TEST01 tests uniqueness counting with no tolerance.
!
!  Discussion:
!
!    POINT_UNIQUE_COUNT uses an O(N) algorithm.
!    POINT_RADIAL_UNIQUE_COUNT uses an algorithm that should be,
!      in general, O(N);
!    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
!
!    For this test, we just want to make sure the algorithms agree
!    in the counting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) n_unique
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) unique_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  To count the unique columns in an R8COL, we call'
  write ( *, '(a)' ) '  POINT_UNIQUE_COUNT,'
  write ( *, '(a)' ) '  POINT_RADIAL_UNIQUE_COUNT, (with random center)'
  write ( *, '(a)' ) '  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M =     ', m
  write ( *, '(a,i8)' ) '  N =     ', n
  write ( *, '(a,i12)' ) '  SEED =  ', seed

  call r8col_duplicates ( m, n, n_unique, seed, a )

  call r8mat_transpose_print ( m, n, a, '  Matrix with N_UNIQUE unique columns:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N_UNIQUE =                  ', n_unique

  call point_unique_count ( m, n, a, unique_num )
  write ( *, '(a,i8)' ) '  POINT_UNIQUE_COUNT =        ', unique_num

  call point_radial_unique_count ( m, n, a, seed, unique_num )
  write ( *, '(a,i8)' ) '  POINT_RADIAL_UNIQUE_COUNT = ', unique_num

  tol = 0.0D+00
  call point_tol_unique_count ( m, n, a, tol, unique_num )
  write ( *, '(a,i8)' ) '  POINT_TOL_UNIQUE_COUNT =    ', unique_num

  return
end
subroutine test02 ( m, n, n_unique, tol, seed )

!*****************************************************************************80
!
!! TEST02 tests uniqueness counting with a tolerance.
!
!  Discussion:
!
!    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
!      in general, O(N);
!    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
!
!    For this test, we just want to make sure the algorithms agree
!    in the counting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n_unique
  real    ( kind = 8 ) r(m)
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) unique_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  To count the unique columns in an R8COL, we call'
  write ( *, '(a)' ) '  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)'
  write ( *, '(a)' ) '  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M =     ', m
  write ( *, '(a,i8)' ) '  N =     ', n
  write ( *, '(a,g14.6)' ) '  TOL =  ', tol
  write ( *, '(a,i12)' ) '  SEED =  ', seed

  call r8col_duplicates ( m, n, n_unique, seed, a )

  call r8mat_transpose_print ( m, n, a, '  Matrix with N_UNIQUE unique columns:' )
!
!  The form of the tolerance test means that if two vectors are initially
!  equal, they remain "tolerably equal" after the addition of random
!  perturbation vectors whose 2-norm is no greater than TOL/2.
!
  do j = 1, n
    call r8vec_uniform_01 ( m, seed, r )
    r(1:m) = r(1:m) / sqrt ( sum ( r(1:m)**2 ) )
    a(1:m,j) = a(1:m,j) + 0.5D+00 * tol * r(1:m)
  end do

  call r8mat_transpose_print ( m, n, a, '  Blurred matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N_UNIQUE =                      ', n_unique

  call point_radial_tol_unique_count ( m, n, a, tol, seed, unique_num )
  write ( *, '(a,i8)' ) '  POINT_RADIAL_TOL_UNIQUE_COUNT = ', unique_num

  call point_tol_unique_count ( m, n, a, tol, unique_num )
  write ( *, '(a,i8)' ) '  POINT_TOL_UNIQUE_COUNT =        ', unique_num

  return
end
subroutine test03 ( m, n, n_unique, tol, seed )

!*****************************************************************************80
!
!! TEST03 compares timings for two uniqueness counters.
!
!  Discussion:
!
!    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
!      in general, O(N);
!    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) ctime1
  real    ( kind = 8 ) ctime2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n_unique
  real    ( kind = 8 ) r(m)
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) unique_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  To count the unique columns in an R8COL, we call'
  write ( *, '(a)' ) '  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)'
  write ( *, '(a)' ) '  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M =     ', m
  write ( *, '(a,i8)' ) '  N =     ', n
  write ( *, '(a,g14.6)' ) '  TOL =  ', tol
  write ( *, '(a,i12)' ) '  SEED =  ', seed

  call r8col_duplicates ( m, n, n_unique, seed, a )
!
!  The form of the tolerance test means that if two vectors are initially
!  equal, they remain "tolerably equal" after the addition of random
!  perturbation vectors whose 2-norm is no greater than TOL/2.
!
  do j = 1, n
    call r8vec_uniform_01 ( m, seed, r )
    r(1:m) = r(1:m) / sqrt ( sum ( r(1:m)**2 ) )
    a(1:m,j) = a(1:m,j) + 0.5D+00 * tol * r(1:m)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N_UNIQUE =                      ', n_unique

  call cpu_time ( ctime1 )
  call point_radial_tol_unique_count ( m, n, a, tol, seed, unique_num )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  POINT_RADIAL_TOL_UNIQUE_COUNT = ', unique_num
  write ( *, '(a,g14.6)' ) '  CPU_TIME = ', ctime2 - ctime1

  call cpu_time ( ctime1 )
  call point_tol_unique_count ( m, n, a, tol, unique_num )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  POINT_TOL_UNIQUE_COUNT =        ', unique_num
  write ( *, '(a,g14.6)' ) '  CPU_TIME = ', ctime2 - ctime1

  return
end
subroutine test04 ( m, n, n_unique, tol, seed )

!*****************************************************************************80
!
!! TEST04 tests uniqueness indexing with a tolerance.
!
!  Discussion:
!
!    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
!      in general, O(N);
!    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
!
!    For this test, we just want to make sure the algorithms agree
!    in the counting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) dist
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n_unique
  real    ( kind = 8 ) r(m)
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) undx(n)
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) xdnu(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  To index the unique columns in an R8COL, we call'
  write ( *, '(a)' ) '  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)'
  write ( *, '(a)' ) '  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M =     ', m
  write ( *, '(a,i8)' ) '  N =     ', n
  write ( *, '(a,g14.6)' ) '  TOL =  ', tol
  write ( *, '(a,i12)' ) '  SEED =  ', seed

  call r8col_duplicates ( m, n, n_unique, seed, a )

  call r8mat_transpose_print ( m, n, a, '  Matrix with N_UNIQUE unique columns:' )
!
!  The form of the tolerance test means that if two vectors are initially
!  equal, they remain "tolerably equal" after the addition of random
!  perturbation vectors whose 2-norm is no greater than TOL/2.
!
  do j = 1, n
    call r8vec_uniform_01 ( m, seed, r )
    r(1:m) = r(1:m) / sqrt ( sum ( r(1:m)**2 ) )
    a(1:m,j) = a(1:m,j) + 0.5D+00 * tol * r(1:m)
  end do

  call r8mat_transpose_print ( m, n, a, '  Blurred matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N_UNIQUE =                      ', n_unique

  call point_radial_tol_unique_index ( m, n, a, tol, seed, unique_num, undx, &
    xdnu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POINT_RADIAL_TOL_UNIQUE_INDEX'
  write ( *, '(a,i8)' ) '  Unique_num = ', unique_num

  call i4vec_print ( unique_num, undx, '  UNDX:' )

  call i4vec_print ( n, xdnu, '  XDNU:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  List of nonunique points P(J), represented by'
  write ( *, '(a)' ) '  point with index I(J).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  J, P(J)'
  write ( *, '(a)' ) '  I(J), P(I(J))'
  write ( *, '(a)' ) '  || P(J) - P(I(J)) || (should be <= TOL.)'
  write ( *, '(a)' ) ' '
  do j = 1, n
    k = undx(xdnu(j))
    if ( j .ne. k ) then
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i4,4(2x,g14.6))' ) j, a(1:m,j)
      write ( *, '(2x,i4,4(2x,g14.6))' ) k, a(1:m,k)
      dist = sqrt ( sum ( ( a(1:m,j) - a(1:m,k) )**2 ) )
      write ( *, '(10x,g10.4)' ) dist
    end if
  end do
!
!  The interpretation of XDNU is simpler for POINT_TOL_UNIQUE_INDEX.
!
  call point_tol_unique_index ( m, n, a, tol, unique_num, xdnu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POINT_TOL_UNIQUE_INDEX'
  write ( *, '(a,i8)' ) '  Unique_num = ', unique_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  List of nonunique points P(J), represented by'
  write ( *, '(a)' ) '  point with index I(J).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  J, P(J)'
  write ( *, '(a)' ) '  I(J), P(I(J))'
  write ( *, '(a)' ) '  || P(J) - P(I(J)) || (should be <= TOL.)'
  write ( *, '(a)' ) ' '
  do j = 1, n
    k = xdnu(j)
    if ( j .ne. k ) then
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i4,4(2x,g14.6))' ) j, a(1:m,j)
      write ( *, '(2x,i4,4(2x,g14.6))' ) k, a(1:m,k)
      dist = sqrt ( sum ( ( a(1:m,j) - a(1:m,k) )**2 ) )
      write ( *, '(10x,g10.4)' ) dist
    end if
  end do

  return
end
subroutine test05 ( m, n, n_unique, tol, seed )

!*****************************************************************************80
!
!! TEST05 times uniqueness indexing with a tolerance.
!
!  Discussion:
!
!    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
!      in general, O(N);
!    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
!
!    For this test, we just want to make sure the algorithms agree
!    in the counting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) ctime1
  real    ( kind = 8 ) ctime2
  real    ( kind = 8 ) dist
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n_unique
  real    ( kind = 8 ) r(m)
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) undx(n)
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) xdnu(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  We time the computations in TEST04, calling'
  write ( *, '(a)' ) '  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)'
  write ( *, '(a)' ) '  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M =     ', m
  write ( *, '(a,i8)' ) '  N =     ', n
  write ( *, '(a,g14.6)' ) '  TOL =  ', tol
  write ( *, '(a,i12)' ) '  SEED =  ', seed

  call r8col_duplicates ( m, n, n_unique, seed, a )
!
!  The form of the tolerance test means that if two vectors are initially
!  equal, they remain "tolerably equal" after the addition of random
!  perturbation vectors whose 2-norm is no greater than TOL/2.
!
  do j = 1, n
    call r8vec_uniform_01 ( m, seed, r )
    r(1:m) = r(1:m) / sqrt ( sum ( r(1:m)**2 ) )
    a(1:m,j) = a(1:m,j) + 0.5D+00 * tol * r(1:m)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N_UNIQUE =                      ', n_unique

  call cpu_time ( ctime1 )
  call point_radial_tol_unique_index ( m, n, a, tol, seed, unique_num, undx, &
    xdnu )
  call cpu_time ( ctime2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POINT_RADIAL_TOL_UNIQUE_INDEX'
  write ( *, '(a,i8)' ) '  Unique_num = ', unique_num
  write ( *, '(a,g14.6)' ) '  Time = ', ctime2 - ctime1

  call cpu_time ( ctime1 )
  call point_tol_unique_index ( m, n, a, tol, unique_num, xdnu )
  call cpu_time ( ctime2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  POINT_TOL_UNIQUE_INDEX'
  write ( *, '(a,i8)' ) '  Unique_num = ', unique_num
  write ( *, '(a,g14.6)' ) '  Time = ', ctime2 - ctime1

  return
end
