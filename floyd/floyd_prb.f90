program main

!*****************************************************************************80
!
!! MAIN is the main program for FLOYD_PRB.
!
!  Discussion:
!
!    FLOYD_PRB calls a set of problems for FLOYD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) wtime

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FLOYD_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FLOYD library.'

  call test01 ( )
  call test02 ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FLOYD_TEST03'
  write ( *, '(a)' ) '  Test I4MAT_FLOYD on the MOD(I,J) matrix.'
  write ( *, '(a)' ) '  The work is roughly N^3.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N   Time (seconds)  Time/N^3'
  write ( *, '(a)' ) ' '

  n = 1
  do while ( n <= 1024 )
    call test03 ( n, wtime )
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) &
      n, wtime, 1000000.0D+00 * wtime / real ( n**3, kind = 8 )
    n = n * 2
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FLOYD_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests I4MAT_FLOYD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ), dimension ( n, n ) :: a = reshape ( (/ &
     0, -1, -1, -1, -1, -1, &
     2,  0, -1, -1, -1,  5, &
     5,  7,  0, -1,  2, -1, &
    -1,  1,  4,  0, -1,  2, &
    -1, -1, -1,  3,  0,  4, &
    -1,  8, -1, -1,  3,  0  &
    /), (/ n, n /) )
  integer ( kind = 4 ) huge
  integer ( kind = 4 ) i4_huge

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  I4MAT_FLOYO uses Floyd''s algorithm to find the'
  write ( *, '(a)' ) '  shortest distance between all pairs of nodes'
  write ( *, '(a)' ) '  in a directed graph, starting from the initial array'
  write ( *, '(a)' ) '  of direct node-to-node distances.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the initial direct distance array, if'
  write ( *, '(a)' ) '    A(I,J) = -1,'
  write ( *, '(a)' ) '  this indicates there is NO directed link from'
  write ( *, '(a)' ) '  node I to node J.  In that case, the value of'
  write ( *, '(a)' ) '  of A(I,J) is essentially "infinity".'

  call i4mat_print ( n, n, a, '  Initial direct distance array:' )

  huge = i4_huge ( ) / 2

  where ( a(1:n,1:n) == - 1 )
    a = huge
  end where

  call i4mat_floyd ( n, a )

  where ( a(1:n,1:n) == huge )
    a = - 1
  end where

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the final shortest distance array, if'
  write ( *, '(a)' ) '    A(I,J) = -1,'
  write ( *, '(a)' ) '  this indicates there is NO directed path from'
  write ( *, '(a)' ) '  node I to node J.'

  call i4mat_print ( n, n, a, '  Final shortest distance array:' )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests R8MAT_FLOYD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ), dimension ( n, n ) :: a = reshape ( (/ &
     0.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, &
     2.0D+00,  0.0D+00, -1.0D+00, -1.0D+00, -1.0D+00,  5.0D+00, &
     5.0D+00,  7.0D+00,  0.0D+00, -1.0D+00,  2.0D+00, -1.0D+00, &
    -1.0D+00,  1.0D+00,  4.0D+00,  0.0D+00, -1.0D+00,  2.0D+00, &
    -1.0D+00, -1.0D+00, -1.0D+00,  3.0D+00,  0.0D+00,  4.0D+00, &
    -1.0D+00,  8.0D+00, -1.0D+00, -1.0D+00,  3.0D+00,  0.0D+00  &
    /), (/ n, n /) )
  real ( kind = 8 ) r8_huge

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  R8MAT_FLOYO uses Floyd''s algorithm to find the'
  write ( *, '(a)' ) '  shortest distance between all pairs of nodes'
  write ( *, '(a)' ) '  in a directed graph, starting from the initial array'
  write ( *, '(a)' ) '  of direct node-to-node distances.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the initial direct distance array, if'
  write ( *, '(a)' ) '    A(I,J) = -1,'
  write ( *, '(a)' ) '  this indicates there is NO directed link from'
  write ( *, '(a)' ) '  node I to node J.  In that case, the value of'
  write ( *, '(a)' ) '  of A(I,J) is essentially "infinity".'

  call r8mat_print ( n, n, a, '  Initial direct distance array:' )

  where ( a(1:n,1:n) == - 1.0D+00 )
    a = r8_huge ( )
  end where

  call r8mat_floyd ( n, a )

  where ( a(1:n,1:n) == r8_huge ( ) )
    a = - 1.0D+00
  end where

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the final shortest distance array, if'
  write ( *, '(a)' ) '    A(I,J) = -1,'
  write ( *, '(a)' ) '  this indicates there is NO directed path from'
  write ( *, '(a)' ) '  node I to node J.'

  call r8mat_print ( n, n, a, '  Final shortest distance array:' )

  return
end
subroutine test03 ( n, wtime )

!*****************************************************************************80
!
!! TEST03 tests I4MAT_FLOYD.
!
!  Discussion:
!
!    The matrix size is input by the user.
!
!    The matrix A has the property that
!
!      A(I,J) = 1 if I is divisible by J.
!
!  Example:
!
!    N = 6
!
!    1 0 0 0 0 0
!    1 1 0 0 0 0
!    1 0 1 0 0 0
!    1 1 0 1 0 0
!    1 0 0 0 1 0
!    1 1 1 0 0 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the matrix.
!
!    Output, real ( kind = 8 ) WTIME, the CPU  time required by I4MAT_FLOYD.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) huge
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) j
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) wtime

  huge = i4_huge ( ) / 2

  do j = 1, n
    do i = 1, n
      if ( mod ( i, j ) == 0 ) then
        a(i,j) = 1
      else
        a(i,j) = huge
      end if
    end do
  end do

  call cpu_time ( time1 )

  call i4mat_floyd ( n, a )

  call cpu_time ( time2 )

  wtime = time2 - time1

  return
end
