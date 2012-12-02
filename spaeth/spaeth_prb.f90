program main

!*****************************************************************************80
!
!! MAIN is the main program for SPAETH_PRB.
!
!  Discussion:
!
!    SPAETH_PRB tests routines from the SPAETH library.
!
!  Modified:
!
!    10 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPAETH_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPAETH library.'

  call test01 ( 'spaeth_01.txt' )
  call test01 ( 'spaeth_02.txt' )
  call test01 ( 'spaeth_03.txt' )
  call test01 ( 'spaeth_04.txt' )
  call test01 ( 'spaeth_05.txt' )
  call test01 ( 'spaeth_06.txt' )
  call test01 ( 'spaeth_07.txt' )
  call test01 ( 'spaeth_08.txt' )
  call test02 ( 'spaeth_09.txt' )
  call test03 ( 'spaeth_01.txt' )
  call test04
  call test05 ( 'spaeth_04.txt' )
  call test06 ( 'spaeth_04.txt' )
  call test07 ( 'spaeth_04.txt' )
  call test08 ( 'spaeth_05.txt' )
  call test09 ( 'spaeth_05.txt' )
  call test10 ( 'spaeth_09.txt' )
  call test11 ( 'spaeth_09.txt' )
  call test12 ( 'spaeth_09.txt' )
  call test13 ( 'spaeth_10.txt' )
  call test14 ( 'spaeth_11.txt', 1 )
  call test14 ( 'spaeth_11.txt', 2 )
  call test14 ( 'spaeth_11.txt', 3 )
  call test15 ( 'spaeth_12.txt' )
  call test16 ( 'spaeth_13.txt' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPAETH_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( file_name )

!*****************************************************************************80
!
!! TEST01 tests DATA_SIZE, DATA_R_READ, DATA_R_SHOW
!
  implicit none

  integer ( kind = 4 ) columns
  character ( len = * ) :: file_name
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rows
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  DATA_SIZE reports the size of a data set.'
  write ( *, '(a)' ) '  DATA_R_READ reads a real data set.'
  write ( *, '(a)' ) '  DATA_R_SHOW makes a plot of them.'

  call data_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a,i6)' ) '  Dimension of data items is ', n
  write ( *, '(a)' ) ' '

  allocate ( x(1:m,1:n) )

  call data_d_read ( file_name, m, n, x )

  j1 = 1
  j2 = 2
  rows = 25
  columns = 50

  call data_d_show ( m, n, x, j1, j2, rows, columns )

  deallocate ( x )

  return
end
subroutine test02 ( file_name )

!*****************************************************************************80
!
!! TEST02 tests DATA_SIZE, DATA_I_READ
!
  implicit none

  character ( len = * ) :: file_name
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  DATA_SIZE reports the size of a data set.'
  write ( *, '(a)' ) '  DATA_I_READ reads an integer data set.'

  call data_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a,i6)' ) '  Dimension of data items is ', n
  write ( *, '(a)' ) ' '

  allocate ( x(1:m,1:n) )

  call data_i_read ( file_name, m, n, x )

  deallocate ( x )

  return
end
subroutine test03 ( file_name )

!*****************************************************************************80
!
!! TEST03 tests TRWEXM.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 149.
!
!  Local Parameters:
!
!    N1, N2, the minimum and maximum number of clusters.
!
!    NRMAX, the number of different initial configurations to try.
!
!    S, the dimension of the data, which is 2 for this example.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) columns
  real ( kind = 8 ) d
  real ( kind = 8 ) e(n)
  character ( len = * ) :: file_name
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: m0 = 1
  integer ( kind = 4 ) mj(n)
  integer ( kind = 4 ) rows
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xbar
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  TRWEXM clusters the data.'
  write ( *, '(a)' ) '  CLUSTER_SHOW shows the clusters.'

  call data_size ( file_name, m, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a,i6)' ) '  Dimension of data items is ', s
  write ( *, '(a,i6)' ) '  Number of clusters is ', n
  write ( *, '(a)' ) ' '

  allocate ( x(1:m,1:s) )
  allocate ( xbar(1:n,1:s) )
  allocate ( z(1:m) )

  call data_d_read ( file_name, m, s, x )
!
!  Consider a number of clusters N = 3.
!
  seed = 37519
!
!  Initially partition the data randomly.
!
  call randp ( m, m0, n, z, seed )
!
!  Carry out the exchange algorithm for the variance criterion.
!
  call trwexm ( x, m, s, z, m0, mj, xbar, n, e, d, it, iflag )

  if ( iflag /= 0 ) then
    return
  end if

  j1 = 1
  j2 = 2
  rows = 25
  columns = 50

  call cluster_d_show ( m, s, x, j1, j2, z, rows, columns )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  M = ', m
  write ( *, '(a,i6)' ) '  N = ', n
  write ( *, '(a,i6)' ) '  IT = ', it
  write ( *, '(a,g14.6)' ) '  D = ', d
  write ( *, '(2x,i3,i6,g14.6)' ) ( j, mj(j), e(j), j = 1, n )

  deallocate ( x )
  deallocate ( xbar )
  deallocate ( z )

  return
end
subroutine test04

!*****************************************************************************80
!
!! TEST04 tests LDLT.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) d(n,n)
  real ( kind = 8 ), parameter :: eps = 0.00001D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  real ( kind = 8 ) l(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  LDLT computes the Cholesky decomposition of'
  write ( *, '(a)' ) '  a general storage positive definite symmetric matrix:'
  write ( *, '(a)' ) '    A = L * D * L'''
  write ( *, '(a)' ) '  where L is an unit lower triangular matrix'
  write ( *, '(a)' ) '  and D is a diagonal matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The matrix order is N = ', n
!
!  Choose a desired factor L:
!
  call random_number ( harvest = l(1:n,1:n) )

  d(1:n,1:n) = 0.0D+00

  do i = 1, n
    d(i,i) = l(i,i) + 0.25D+00
    l(i,i) = 1.0D+00
    do j = i+1, n
      l(i,j) = 0.0D+00
    end do
  end do

  call r8mat_print ( n, n, l, '  The desired Cholesky factor L:' )

  call r8mat_print ( n, n, d, '  The desired Cholesky factor D:' )
!
!  Construct the matrix A.
!
  a(1:n,1:n) = matmul ( &
    matmul ( l(1:n,1:n), d(1:n,1:n) ), transpose ( l(1:n,1:n) ) )

  call r8mat_print ( n, n, a, '  The matrix A:' )
!
!  Call LDLT to factor the matrix A.
!
  call ldlt ( n, a, eps, iflag )

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  LDLT returned IFLAG = ', iflag
    write ( *, '(a)' ) '  This means the matrix is not positive definite.'
    return
  end if
!
!  Now we print the factorization matrix U.
!
  l(1:n,1:n) = 0.0D+00

  do i = 1, n
    l(i,i) = 1.0D+00
  end do

  do i = 1, n
    do j = 1, i-1
      l(i,j) = a(i,j)
    end do
  end do

  d(1:n,1:n) = 0.0D+00
  do i = 1, n
    d(i,i) = a(i,i)
  end do

  call r8mat_print ( n, n, l, '  The computed Cholesky factor L:' )

  call r8mat_print ( n, n, d, '  The computed Cholesky factor D:' )
!
!  Compute the Cholesky product.
!
  a(1:n,1:n) = matmul ( &
    matmul ( l(1:n,1:n), d(1:n,1:n) ), transpose ( l(1:n,1:n) ) )

  call r8mat_print ( n, n, a, '  The computed product L * D * L'':' )

  return
end
subroutine test05 ( file_name )

!*****************************************************************************80
!
!! TEST05 tests TRWEXM.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 149.
!
!  Local Parameters:
!
!    N1, N2, the minimum and maximum number of clusters.
!
!    NRMAX, the number of different initial configurations to try.
!
!    S, the dimension of the data, which is 2 for this example.
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 2
  integer ( kind = 4 ), parameter :: n2 = 6
  integer ( kind = 4 ), parameter :: nrmax = 20

  real ( kind = 8 ) d
  real ( kind = 8 ) d_min
  real ( kind = 8 ) dnnr(n1:n2,nrmax)
  real ( kind = 8 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kit(n1:n2,nrmax)
  integer ( kind = 4 ), parameter :: kprint = 0
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: m0 = 1
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mj
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr
  real ( kind = 8 ), parameter :: r = 0.999D+00
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xbar
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  TRWEXM clusters the data.'

  call data_size ( file_name, m, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a,i6)' ) '  Dimension of data items is ', s
  write ( *, '(a)' ) ' '

  allocate ( x(1:m,1:s) )
  allocate ( z(1:m) )

  call data_d_read ( file_name, m, s, x )
!
!  Consider a number of clusters N.
!
  do n = n1, n2

    allocate ( xbar(1:n,1:s) )
    allocate ( e(1:n) )
    allocate ( mj(1:n) )

    d_min = huge ( d_min )
    seed = 37519
!
!  Try NRMAX different starting configurations.
!
    do nr = 1, nrmax
!
!  Initially partition the data randomly.
!
      call randp ( m, m0, n, z, seed )
!
!  Carry out the exchange algorithm for the variance criterion.
!
      call trwexm ( x, m, s, z, m0, mj, xbar, n, e, d, it, iflag )

      if ( iflag /= 0 ) then
        dnnr(n,nr) = 0.0D+00
        kit(n,nr) = -iflag
        cycle
      end if

      dnnr(n,nr) = d
      kit(n,nr) = it

      if ( d < d_min * r ) then

        if ( kprint /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) '  M = ', m
          write ( *, '(a,i6)' ) '  N = ', n
          write ( *, '(a,i6)' ) '  NR = ', nr
          write ( *, '(a,i6)' ) '  IT = ', it
          write ( *, '(a,g14.6)' ) '  D = ', d
          write ( *, '(2x,i3,i6,g14.6)' ) ( j, mj(j), e(j), j = 1, n )
        end if

        d_min = d

      end if

    end do

    deallocate ( e )
    deallocate ( mj )
    deallocate ( xbar )

  end do

  do j = n1, n2

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'Cluster size = ', j
    write ( *, '(a)' ) ' '

    do nr = 1, nrmax
      write ( *, '(2x,i6,i6,g14.6)' ) nr, kit(j,nr), dnnr(j,nr)
    end do

  end do

  deallocate ( x )
  deallocate ( z )

  return
end
subroutine test06 ( file_name )

!*****************************************************************************80
!
!! TEST06 tests TRWMDM.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 149.
!
!  Local Parameters:
!
!    N1, N2, the minimum and maximum number of clusters.
!
!    NRMAX, the number of different initial configurations to try.
!
!    S, the dimension of the data, which is 2 for this example.
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 2
  integer ( kind = 4 ), parameter :: n2 = 6
  integer ( kind = 4 ), parameter :: nrmax = 20

  real ( kind = 8 ) d
  real ( kind = 8 ) d_min
  real ( kind = 8 ) dnnr(n1:n2,nrmax)
  real ( kind = 8 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kit(n1:n2,nrmax)
  integer ( kind = 4 ), parameter :: kprint = 0
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: m0 = 1
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mj
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr
  real ( kind = 8 ), parameter :: r = 0.999D+00
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xbar
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  TRWMDM clusters the data.'

  call data_size ( file_name, m, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a,i6)' ) '  Dimension of data items is ', s
  write ( *, '(a)' ) ' '

  allocate ( x(1:m,1:s) )
  allocate ( z(1:m) )

  call data_d_read ( file_name, m, s, x )
!
!  Consider a number of clusters N.
!
  do n = n1, n2

    allocate ( xbar(1:n,1:s) )
    allocate ( e(1:n) )
    allocate ( mj(1:n) )

    d_min = huge ( d_min )
    seed = 37519
!
!  Try NRMAX different starting configurations.
!
    do nr = 1, nrmax
!
!  Initially partition the data randomly.
!
      call randp ( m, m0, n, z, seed )
!
!  Carry out the exchange algorithm for the variance criterion.
!
      call trwmdm ( x, m, s, z, m0, mj, xbar, n, e, d, it, iflag )

      if ( iflag /= 0 ) then
        dnnr(n,nr) = 0.0D+00
        kit(n,nr) = -iflag
        cycle
      end if

      dnnr(n,nr) = d
      kit(n,nr) = it

      if ( d < d_min * r ) then

        if ( kprint /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) '  M = ', m
          write ( *, '(a,i6)' ) '  N = ', n
          write ( *, '(a,i6)' ) '  NR = ', nr
          write ( *, '(a,i6)' ) '  IT = ', it
          write ( *, '(a,g14.6)' ) '  D = ', d
          write ( *, '(2x,i3,i6,g14.6)' ) ( j, mj(j), e(j), j = 1, n )
        end if

        d_min = d

      end if

    end do

    deallocate ( e )
    deallocate ( mj )
    deallocate ( xbar )

  end do

  do j = n1, n2

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Cluster size = ', j
    write ( *, '(a)' ) ' '

    do nr = 1, nrmax
      write ( *, '(2x,i6,i6,g14.6)' ) nr, kit(j,nr), dnnr(j,nr)
    end do

  end do

  deallocate ( x )
  deallocate ( z )

  return
end
subroutine test07 ( file_name )

!*****************************************************************************80
!
!! TEST07 tests TRWMDM, TRWEXM.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 149.
!
!  Local Parameters:
!
!    N1, N2, the minimum and maximum number of clusters.
!
!    NRMAX, the number of different initial configurations to try.
!
!    S, the dimension of the data, which is 2 for this example.
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 2
  integer ( kind = 4 ), parameter :: n2 = 6
  integer ( kind = 4 ), parameter :: nrmax = 20

  real ( kind = 8 ) d
  real ( kind = 8 ) d_min
  real ( kind = 8 ) dnnr(n1:n2,nrmax)
  real ( kind = 8 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kit(n1:n2,nrmax)
  integer ( kind = 4 ), parameter :: kprint = 0
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: m0 = 1
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mj
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr
  real ( kind = 8 ), parameter :: r = 0.999D+00
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xbar
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  TRWMDM clusters the data.'
  write ( *, '(a)' ) '  TRWEXM clusters the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The TRWMDM clusters are improved by TRWEXM.'

  call data_size ( file_name, m, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a,i6)' ) '  Dimension of data items is ', s
  write ( *, '(a)' ) ' '

  allocate ( x(1:m,1:s) )
  allocate ( z(1:m) )

  call data_d_read ( file_name, m, s, x )
!
!  Consider a number of clusters N.
!
  do n = n1, n2

    allocate ( xbar(1:n,1:s) )
    allocate ( e(1:n) )
    allocate ( mj(1:n) )

    d_min = huge ( d_min )
    seed = 37519
!
!  Try NRMAX different starting configurations.
!
    do nr = 1, nrmax
!
!  Initially partition the data randomly.
!
      call randp ( m, m0, n, z, seed )
!
!  Carry out the minimal distance algorithm for the variance criterion.
!
      call trwmdm ( x, m, s, z, m0, mj, xbar, n, e, d, it, iflag )

      if ( iflag /= 0 ) then
        dnnr(n,nr) = 0.0
        kit(n,nr) = -iflag
        cycle
      end if
!
!  The output clusters from the minimal distance algorithm are
!  fed into the exchange algorithm.
!
      call trwmdm ( x, m, s, z, m0, mj, xbar, n, e, d, it, iflag )

      if ( iflag /= 0 ) then
        dnnr(n,nr) = 0.0D+00
        kit(n,nr) = -iflag
        cycle
      end if

      dnnr(n,nr) = d
      kit(n,nr) = it

      if ( d < d_min * r ) then

        if ( kprint /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) '  M = ', m
          write ( *, '(a,i6)' ) '  N = ', n
          write ( *, '(a,i6)' ) '  NR = ', nr
          write ( *, '(a,i6)' ) '  IT = ', it
          write ( *, '(a,g14.6)' ) '  D = ', d
          write ( *, '(2x,i3,i6,g14.6)' ) ( j, mj(j), e(j), j = 1, n )
        end if

        d_min = d

      end if

    end do

    deallocate ( e )
    deallocate ( mj )
    deallocate ( xbar )

  end do

  do j = n1, n2

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Cluster size = ', j
    write ( *, '(a)' ) ' '

    do nr = 1, nrmax
      write ( *, '(2x,i6,i6,g14.6)' ) nr, kit(j,nr), dnnr(j,nr)
    end do

  end do

  deallocate ( x )
  deallocate ( z )

  return
end
subroutine test08 ( file_name )

!*****************************************************************************80
!
!! TEST08 tests DETEXM.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 164.
!
!  Local Parameters:
!
!    N1, N2, the minimum and maximum number of clusters.
!
!    NRMAX, the number of different initial configurations to try.
!
!    S, the dimension of the data, which is 2 for this example.
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 2
  integer ( kind = 4 ), parameter :: n2 = 6
  integer ( kind = 4 ), parameter :: nrmax = 20

  real ( kind = 8 ) detw_min
  real ( kind = 8 ) detw
  real ( kind = 8 ) dnnr(n1:n2,nrmax)
  character ( len = * ) :: file_name
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kit(n1:n2,nrmax)
  integer ( kind = 4 ), parameter :: kprint = 0
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: m0 = 1
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mj
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr
  real ( kind = 8 ), parameter :: r = 0.999D+00
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xbar
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  DETEXM clusters the data.'

  call data_size ( file_name, m, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a,i6)' ) '  Dimension of data items is ', s
  write ( *, '(a)' ) ' '

  allocate ( x(1:m,1:s) )
  allocate ( z(1:m) )

  call data_d_read ( file_name, m, s, x )
!
!  Normalize the data.
!
  call trafor ( x, m, s )
!
!  Consider a number of clusters N.
!
  do n = n1, n2

    allocate ( xbar(1:n,1:s) )
    allocate ( mj(1:n) )

    detw_min = huge ( detw_min )
    seed = 37519
!
!  Try NRMAX different starting configurations.
!
    do nr = 1, nrmax
!
!  Initially partition the data randomly.
!
      call randp ( m, m0, n, z, seed )
!
!  Carry out the exchange algorithm for the determinant criterion.
!
      call detexm ( x, m, s, z, m0, mj, xbar, n, detw, it, iflag )

      if ( iflag /= 0 ) then
        dnnr(n,nr) = 0.0D+00
        kit(n,nr) = -iflag
        cycle
      end if

      dnnr(n,nr) = detw
      kit(n,nr) = it

      if ( detw < detw_min * r ) then

        if ( kprint /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) '  M = ', m
          write ( *, '(a,i6)' ) '  N = ', n
          write ( *, '(a,i6)' ) '  NR = ', nr
          write ( *, '(a,i6)' ) '  IT = ', it
          write ( *, '(a,g14.6)' ) '  DETW = ', detw
          write ( *, '(2x,i3,i6)' ) ( j, mj(j), j = 1, n )
        end if

        detw_min = detw

      end if

    end do

    deallocate ( mj )
    deallocate ( xbar )

  end do

  do j = n1, n2

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Cluster size = ', j
    write ( *, '(a)' ) ' '

    do nr = 1, nrmax
      write ( *, '(2x,i6,i6,g14.6)' ) nr, kit(j,nr), dnnr(j,nr)
    end do

  end do

  deallocate ( x )
  deallocate ( z )

  return
end
subroutine test09 ( file_name )

!*****************************************************************************80
!
!! TEST09 tests DWBEXM.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 170.
!
!  Local Parameters:
!
!    N1, N2, the minimum and maximum number of clusters.
!
!    NRMAX, the number of different initial configurations to try.
!
!    S, the dimension of the data, which is 2 for this example.
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 2
  integer ( kind = 4 ), parameter :: n2 = 6
  integer ( kind = 4 ), parameter :: nrmax = 20

  real ( kind = 8 ) beta
  real ( kind = 8 ) d
  real ( kind = 8 ) d_min
  real ( kind = 8 ) dnnr(n1:n2,nrmax)
  real ( kind = 8 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kit(n1:n2,nrmax)
  integer ( kind = 4 ), parameter :: kprint = 0
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m0
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mj
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr
  real ( kind = 8 ), parameter :: r = 0.999D+00
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: xbar
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  DWBEXM clusters the data.'

  beta = 0.5D+00
  m0 = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Minimal cluster population M0 = ', m0
  write ( *, '(a,g14.6)' ) '  Determinant exponent BETA = ', beta

  call data_size ( file_name, m, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a,i6)' ) '  Dimension of data items is ', s
  write ( *, '(a)' ) ' '

  allocate ( x(1:m,1:s) )
  allocate ( z(1:m) )

  call data_d_read ( file_name, m, s, x )
!
!  Normalize the data.
!
  call trafor ( x, m, s )
!
!  Consider a number of clusters N.
!
  do n = n1, n2

    allocate ( xbar(1:n,1:s) )
    allocate ( e(1:n) )
    allocate ( mj(1:n) )

    d_min = huge ( d_min )
    seed = 37519
!
!  Try NRMAX different starting configurations.
!
    do nr = 1, nrmax
!
!  Initially partition the data randomly.
!
      call randp ( m, m0, n, z, seed )
!
!  Carry out the exchange algorithm for the determinant criterion.
!
      call dwbexm ( x, m, s, z, m0, mj, xbar, n, beta, e, d, it, iflag )

      if ( iflag /= 0 ) then
        dnnr(n,nr) = 0.0D+00
        kit(n,nr) = -iflag
        cycle
      end if

      dnnr(n,nr) = d
      kit(n,nr) = it

      if ( d < d_min * r ) then

        if ( kprint /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) '  M = ', m
          write ( *, '(a,i6)' ) '  N = ', n
          write ( *, '(a,i6)' ) '  NR = ', nr
          write ( *, '(a,i6)' ) '  IT = ', it
          write ( *, '(a,g14.6)' ) '  D = ', d
          write ( *, '(2x,i3,i6,g14.6)' ) ( j, mj(j), e(j), j = 1, n )
        end if

        d_min = d

      end if

    end do

    deallocate ( e )
    deallocate ( mj )
    deallocate ( xbar )

  end do

  do j = n1, n2

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Cluster size = ', j
    write ( *, '(a)' ) ' '

    do nr = 1, nrmax
      write ( *, '(2x,i6,i6,g14.6)' ) nr, kit(j,nr), dnnr(j,nr)
    end do

  end do

  deallocate ( x )
  deallocate ( z )

  return
end
subroutine test10 ( file_name )

!*****************************************************************************80
!
!! TEST10 tests OVSEXM.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 186.
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 4
  integer ( kind = 4 ), parameter :: n2 = 4
  integer ( kind = 4 ), parameter :: nrmax = 20

  integer ( kind = 4 ) d
  logical, parameter :: debug = .false.
  integer ( kind = 4 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kd
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: m0 = 1
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mj
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: x
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  OVSEXM clusters integer ordinal data.'

  call data_size ( file_name, m, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a,i6)' ) '  Dimension of data items is ', s
  write ( *, '(a)' ) ' '

  allocate ( p(1:m) )
  allocate ( x(1:m,1:s) )
  allocate ( z(1:m) )

  call data_i_read ( file_name, m, s, x )

  if ( debug ) then
    call data_i_print ( m, s, x, '  The data matrix:' )
  end if
!
!  Consider a number of clusters N.
!
  do n = n1, n2

    allocate ( e(1:n) )
    allocate ( mj(1:n) )

    seed = 37519
    kd = huge ( kd )
!
!  Try NRMAX different starting configurations.
!
    write ( *, '(a)' ) &
      '  M N NR IT      D     J MJ      E  J MJ      E  J MJ      E  J MJ      E'
    write ( *, '(a)' ) ' '

    do nr = 1, nrmax
!
!  Initially partition the data randomly.
!
      call randp ( m, m0, n, z, seed )
!
!  Carry out the exchange algorithm for the determinant criterion.
!
      call ovsexm ( x, m, s, z, mj, n, e, d, it, iflag )

      if ( iflag /= 0 ) then
        write ( *, '(a,i6)' ) '  OVSEXM returned IFLAG = ', iflag
        cycle
      end if

      write ( *, '(2x,i3,i2,i3,i3,i7,3x,9(i3,i4,i6))' ) &
        m, n, nr, it, d, ( j, mj(j), e(j), j = 1, n )

      if ( d < kd ) then
        kd = d
        p(1:m) = z(1:m)
      end if

    end do

    deallocate ( e )
    deallocate ( mj )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a,i6)' ) '  NRMAX = ', nrmax
    write ( *, '(a)' ) ' '
    do j = 1, n
      write ( *, '(2x,i6)' ) j
      do i = 1, m
        if ( j == p(i) ) then
          write ( *, '(2x,i4,3x,67i1)' ) i, x(i,1:s)
        end if
      end do
    end do

  end do

  deallocate ( p )
  deallocate ( x )
  deallocate ( z )

  return
end
subroutine test11 ( file_name )

!*****************************************************************************80
!
!! TEST11 tests OVREXM.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 189.
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 4
  integer ( kind = 4 ), parameter :: n2 = 4
  integer ( kind = 4 ), parameter :: nrmax = 20

  integer ( kind = 4 ) d
  logical, parameter :: debug = .false.
  integer ( kind = 4 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kd
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: m0 = 1
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mj
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: x
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  OVREXM clusters integer ordinal data.'

  call data_size ( file_name, m, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a,i6)' ) '  Dimension of data items is ', s
  write ( *, '(a)' ) ' '

  allocate ( p(1:m) )
  allocate ( x(1:m,1:s) )
  allocate ( z(1:m) )

  call data_i_read ( file_name, m, s, x )

  if ( debug ) then
    call data_i_print ( m, s, x, '  The data matrix:' )
  end if
!
!  Consider a number of clusters N.
!
  do n = n1, n2

    allocate ( e(1:n) )
    allocate ( mj(1:n) )

    seed = 37519
    kd = huge ( kd )
!
!  Try NRMAX different starting configurations.
!
    write ( *, '(a)' ) &
      '  M N NR IT      D     J MJ      E  J MJ      E  J MJ      E  J MJ      E'
    write ( *, '(a)' ) ' '

    do nr = 1, nrmax
!
!  Initially partition the data randomly.
!
      call randp ( m, m0, n, z, seed )
!
!  Carry out the exchange algorithm for the determinant criterion.
!
      call ovrexm ( x, m, s, z, mj, n, e, d, it, iflag )

      if ( iflag /= 0 ) then
        write ( *, '(a,i6)' ) '  OVREXM returns IFLAG = ', iflag
        cycle
      end if

      write ( *, '(2x,i3,i2,i3,i3,i7,3x,9(i3,i4,i6))' ) &
        m, n, nr, it, d, ( j, mj(j), e(j), j = 1, n )

      if ( d < kd ) then
        kd = d
        p(1:m) = z(1:m)
      end if

    end do

    deallocate ( e )
    deallocate ( mj )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a,i6)' ) '  NRMAX = ', nrmax
    write ( *, '(a)' ) ' '
    do j = 1, n
      write ( *, '(2x,i6)' ) j
      do i = 1, m
        if ( j == p(i) ) then
          write ( *, '(2x,i4,3x,67i1)' ) i, x(i,1:s)
        end if
      end do
    end do

  end do

  deallocate ( p )
  deallocate ( x )
  deallocate ( z )

  return
end
subroutine test12 ( file_name )

!*****************************************************************************80
!
!! TEST12 tests OVPEXM.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 189.
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 4
  integer ( kind = 4 ), parameter :: n2 = 4
  integer ( kind = 4 ), parameter :: nrmax = 20

  integer ( kind = 4 ) d
  logical, parameter :: debug = .false.
  integer ( kind = 4 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kd
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: m0 = 1
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mj
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: t = 4
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: x
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  OVPEXM clusters integer ordinal data.'

  call data_size ( file_name, m, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a,i6)' ) '  Dimension of data items is ', s
  write ( *, '(a,i6)' ) '  The ordinal range is 1 to ', t
  write ( *, '(a)' ) ' '

  allocate ( p(1:m) )
  allocate ( x(1:m,1:s) )
  allocate ( z(1:m) )

  call data_i_read ( file_name, m, s, x )

  if ( debug ) then
    call data_i_print ( m, s, x, '  The data matrix:' )
  end if
!
!  Consider a number of clusters N.
!
  do n = n1, n2

    allocate ( e(1:n) )
    allocate ( mj(1:n) )

    seed = 37519
    kd = huge ( kd )
!
!  Try NRMAX different starting configurations.
!
    write ( *, '(a)' ) &
      '  M N NR IFLAG IT      D     J MJ      E  J MJ      E  J MJ      E  J MJ      E'
    write ( *, '(a)' ) ' '

    do nr = 1, nrmax
!
!  Initially partition the data randomly.
!
      call randp ( m, m0, n, z, seed )
!
!  Carry out the exchange algorithm for the determinant criterion.
!
      call ovpexm ( x, m, s, z, mj, n, e, d, it, t, iflag )

      if ( iflag /= 0 ) then
        cycle
      end if

      write ( *, '(2x,i3,i2,i3,i6,i3,i7,3x,9(i3,i4,i6))' ) &
        m, n, nr, iflag, it, d, ( j, mj(j), e(j), j = 1, n )

      if ( d < kd ) then
        kd = d
        p(1:m) = z(1:m)
      end if

    end do

    deallocate ( e )
    deallocate ( mj )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a,i6)' ) '  NRMAX = ', nrmax
    write ( *, '(a)' ) ' '
    do j = 1, n
      write ( *, '(2x,i6)' ) j
      do i = 1, m
        if ( j == p(i) ) then
          write ( *, '(2x,i4,3x,67i1)' ) i, x(i,1:s)
        end if
      end do
    end do

  end do

  deallocate ( p )
  deallocate ( x )
  deallocate ( z )

  return
end
subroutine test13 ( file_name )

!*****************************************************************************80
!
!! TEST13 tests BVPEXM.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 192.
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 2
  integer ( kind = 4 ), parameter :: n2 = 6
  integer ( kind = 4 ), parameter :: nrmax = 20

  integer ( kind = 4 ) d
  logical, parameter :: debug = .false.
  integer ( kind = 4 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kd
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: m0 = 1
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mj
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: x
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  BVPEXM clusters binary (0,1) data.'

  call data_size ( file_name, m, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a,i6)' ) '  Dimension of data items is ', s
  write ( *, '(a)' ) ' '

  allocate ( p(1:m) )
  allocate ( x(1:m,1:s) )
  allocate ( z(1:m) )

  call data_i_read ( file_name, m, s, x )

  if ( debug ) then
    call data_i_print ( m, s, x, '  The data matrix:' )
  end if
!
!  Consider a number of clusters N.
!
  do n = n1, n2

    allocate ( e(1:n) )
    allocate ( mj(1:n) )

    seed = 37519
    kd = huge ( kd )
!
!  Try NRMAX different starting configurations.
!
    write ( *, '(a)' ) &
      '  M N NR IFLAG IT      D     J MJ      E  J MJ      E  J MJ      E  J MJ      E'
    write ( *, '(a)' ) ' '

    do nr = 1, nrmax
!
!  Initially partition the data randomly.
!
      call randp ( m, m0, n, z, seed )
!
!  Carry out the exchange algorithm.
!
      call bvpexm ( x, m, s, z, mj, n, e, d, it, iflag )

      if ( iflag /= 0 ) then
        cycle
      end if

      write ( *, '(2x,i3,i2,i3,i6,i3,i7,3x,9(i3,i4,i6))' ) &
        m, n, nr, iflag, it, d, ( j, mj(j), e(j), j = 1, n )

      if ( d < kd ) then
        kd = d
        p(1:m) = z(1:m)
      end if

    end do

    deallocate ( e )
    deallocate ( mj )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a,i6)' ) '  NRMAX = ', nrmax
    write ( *, '(a)' ) ' '
    do j = 1, n
      write ( *, '(2x,i6)' ) j
      do i = 1, m
        if ( j == p(i) ) then
          write ( *, '(2x,i4,3x,67i1)' ) i, x(i,1:s)
        end if
      end do
    end do

  end do

  deallocate ( p )
  deallocate ( x )
  deallocate ( z )

  return
end
subroutine test14 ( file_name, method )

!*****************************************************************************80
!
!! TEST14 tests TIHEXM.
!
!  Modified:
!
!    26 April 2002
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 197-198.
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 2
  integer ( kind = 4 ), parameter :: n2 = 6
  integer ( kind = 4 ), parameter :: nrmax = 20

  real ( kind = 8 ) d
  logical, parameter :: debug = .false.
  real ( kind = 8 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kd
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: m0 = 1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) method
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mj
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: t
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  TIHEXM clusters data using criteria that do not'
  write ( *, '(a)' ) '  rely on the computation of centers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Three different methods are available.'
  write ( *, '(a,i6)' ) '  In this run, we use method ', method
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We assume access to a pair-wise distance matrix.'
!
!  Read the distance matrix.
!
  call data_size ( file_name, m, m2 )

  if ( m /= m2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST14'
    write ( *, '(a)' ) '  Distance matrix file does not have same number'
    write ( *, '(a)' ) '  of rows and columns.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data items is ', m
  write ( *, '(a)' ) ' '

  allocate ( p(1:m) )
  allocate ( t(1:m,1:m) )
  allocate ( z(1:m) )

  call data_d_read ( file_name, m, m, t )

  if ( debug ) then
    call data_d_print ( m, m, t, '  The distance matrix:' )
  end if
!
!  Consider a number of clusters N.
!
  do n = n1, n2

    allocate ( e(1:n) )
    allocate ( mj(1:n) )

    seed = 37519
    kd = huge ( kd )
!
!  Try NRMAX different starting configurations.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  M N NR IFLAG IT      D     J MJ      E  J MJ      E  J MJ      E  J MJ      E'
    write ( *, '(a)' ) ' '

    do nr = 1, nrmax
!
!  Initially partition the data randomly.
!
      call randp ( m, m0, n, z, seed )
!
!  Carry out the exchange algorithm.
!
      call tihexm ( m, t, n, mj, method, z, e, d, it, iflag )

      if ( iflag /= 0 ) then
        cycle
      end if

      write ( *, '(2x,i3,i2,i3,i6,i3,f12.1,3x,9(i3,i4,f12.1))' ) &
        m, n, nr, iflag, it, d, ( j, mj(j), e(j), j = 1, n )

      if ( d < kd ) then
        kd = d
        p(1:m) = z(1:m)
      end if

    end do

    deallocate ( e )
    deallocate ( mj )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Clustering for N = ', n
    write ( *, '(a)' ) ' '

    write ( *, '(2x,20i3)' ) p(1:m)

  end do

  deallocate ( p )
  deallocate ( t )
  deallocate ( z )

  return
end
subroutine test15 ( file_name )

!*****************************************************************************80
!
!! TEST15 tests DATA_D_WRITE.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 201.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 42

  real ( kind = 8 ) f
  character ( len = * ) :: file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(m,m)
  real ( kind = 8 ) urand

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  DATA_D_WRITE writes real data to a file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we create a distance matrix'
  write ( *, '(a)' ) '  used in one of the Spaeth examples and write'
  write ( *, '(a)' ) '  it to a file.'
  write ( *, '(a)' ) ' '

  seed = 37519281

  do i = 1, m
    do j = 1, m

      if ( j < i ) then

        f = urand ( seed )

        if ( f <= 0.3D+00 ) then
          f = f * 15.0D+00
        else
          f = f * 99.9D+00
        end if

        t(i,j) = real ( int ( f ), kind = 8 )
        t(j,i) = t(i,j)

      else if ( j == i ) then

        t(i,j) = 0.0D+00

      else

      end if

    end do
  end do

  call data_d_write ( file_name, m, m, t )

  return
end
subroutine test16 ( file_name )

!*****************************************************************************80
!
!! TEST16 tests CLREXM.
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Dissection and Analysis,
!    Theory, FORTRAN Programs, Examples,
!    Ellis Horwood, 1985, page 206-207.
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 2
  integer ( kind = 4 ), parameter :: n2 = 5
  integer ( kind = 4 ), parameter :: nrmax = 20

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ), allocatable, dimension ( : ) :: b
  real ( kind = 8 ) d
  logical, parameter :: debug = .false.
  real ( kind = 8 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kd
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m0
  integer ( kind = 4 ), allocatable, dimension ( : ) :: mj
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nr
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) sp1
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: temp
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: y
  integer ( kind = 4 ), allocatable, dimension ( : ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  CLREXM clusters data using criteria that do not'
  write ( *, '(a)' ) '  rely on the computation of centers.'
!
!  Read the data matrix.
!
  call data_size ( file_name, m, sp1 )

  s = sp1 - 1
  m0 = s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  Number of data rows is ', m
  write ( *, '(a,i6)' ) '  Number of data columns is ', sp1
  write ( *, '(a)' ) ' '

  allocate ( a(1:m,1:s) )
  allocate ( b(1:m) )
  allocate ( p(1:m) )
  allocate ( temp(1:m,1:sp1) )
  allocate ( z(1:m) )

  call data_d_read ( file_name, m, sp1, temp )

  if ( debug ) then
    call data_d_print ( m, s, temp, '  The regression data:' )
  end if

  b(1:m) = temp(1:m,1)
  a(1:m,1:s) = temp(1:m,2:sp1)

  deallocate ( temp )
!
!  Consider a number of clusters N.
!
  do n = n1, n2

    allocate ( e(1:n) )
    allocate ( mj(1:n) )
    allocate ( x(1:n,1:s) )
    allocate ( y(1:n,1:s) )

    seed = 37519
    kd = huge ( kd )
!
!  Try NRMAX different starting configurations.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  M S N M0 NR IT           D     J MJ      E  J MJ      E  J MJ      E  J MJ      E'
    write ( *, '(a)' ) ' '

    do nr = 1, nrmax
!
!  Initially partition the data randomly.
!
      call randp ( m, m0, n, z, seed )
!
!  Carry out the exchange algorithm.
!
      call clrexm ( a, m, s, b, n, z, m0, mj, x, e, d, it, iflag )

      if ( iflag /= 0 ) then
        cycle
      end if

      write ( *, '(2x,i3,i2,i2,i3,i3,i3,g14.6,3x,9(i3,i4,g14.6))' ) &
        m, s, n, m0, nr, it, d, ( j, mj(j), e(j), j = 1, n )

      if ( d < kd ) then
        kd = d
        p(1:m) = z(1:m)
        y(1:n,1:s) = x(1:n,1:s)
      end if

    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Clustering for N = ', n
    write ( *, '(a)' ) ' '
    write ( *, '(2x,20i3)' ) p(1:m)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Cluster coefficients:'
    write ( *, '(a)' ) ' '
    do j = 1, n
      write ( *, '(2x,i3,5g14.6)' ) j, y(j,1:s)
    end do

    deallocate ( e )
    deallocate ( mj )
    deallocate ( x )
    deallocate ( y )

  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( p )
  deallocate ( z )

  return
end
