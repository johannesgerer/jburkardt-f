program main

!*****************************************************************************80
!
!! MAIN is the main program for SPAETH2_PRB.
!
!  Discussion:
!
!    SPAETH2_PRB tests routines from the SPAETH2 library.
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
  write ( *, '(a)' ) 'SPAETH2_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPAETH2 library.'

  call test01 ( 'spaeth2_03.txt' )
  call test01 ( 'spaeth2_04.txt' )
  call test01 ( 'spaeth2_05.txt' )

  call test02 ( 'spaeth2_03.txt' )
  call test03 ( 'spaeth2_03.txt' )
  call test04 ( 'spaeth2_03.txt' )
  call test05 ( 'spaeth2_03.txt' )
  call test06 ( 'spaeth2_03.txt' )
  call test07 ( 'spaeth2_03.txt' )
  call test08 ( )
  call test09 ( 'spaeth2_03.txt' )
  call test10 ( 'spaeth2_02.txt' )
  call test11 ( 'spaeth2_03.txt' )
  call test12 ( )
  call test13 ( )
  call test14 ( 'spaeth2_03.txt' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPAETH2_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( file_name )

!*****************************************************************************80
!
!! TEST01 tests DATA_SIZE, DATA_R_READ, DATA_R_SHOW.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
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
  write ( *, '(a,i6)' ) '  The number of data items is M =        ', m
  write ( *, '(a,i6)' ) '  The dimension of the data items is N = ', n
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
!! TEST02 tests LEADER.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = * ) :: file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  real ( kind = 8 ) rho
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  LEADER uses a simple clustering algorithm.'

  call data_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  The number of data items is M =        ', m
  write ( *, '(a,i6)' ) '  The dimension of the data items is N = ', n
  write ( *, '(a)' ) ' '

  allocate ( p(1:m) )
  allocate ( x(1:m,1:n) )

  call data_d_read ( file_name, m, n, x )

  do i = 1, 4

    rho = 10.0E+00 * real ( i, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Using a cluster parameter of RHO = ', rho

    call leader ( m, n, x, rho, p, nc )

    write ( *, '(a,i6)' ) '  Number of clusters created was ', nc

  end do

  deallocate ( p )
  deallocate ( x )

  return
end
subroutine test03 ( file_name )

!*****************************************************************************80
!
!! TEST03 tests JOINER.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = * ) :: file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  real ( kind = 8 ) rho
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  LEADER uses a simple clustering algorithm.'

  call data_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  The number of data items is M =        ', m
  write ( *, '(a,i6)' ) '  The dimension of the data items is N = ', n
  write ( *, '(a)' ) ' '

  allocate ( p(1:m) )
  allocate ( x(1:m,1:n) )

  call data_d_read ( file_name, m, n, x )

  do i = 1, 4

    rho = 10.0E+00 * real ( i, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Using a cluster parameter of RHO = ', rho

    call joiner ( m, n, x, rho, p, nc )

    write ( *, '(a,i6)' ) '  Number of clusters created was ', nc

  end do

  deallocate ( p )
  deallocate ( x )

  return
end
subroutine test04 ( file_name )

!*****************************************************************************80
!
!! TEST04 tests ZWEIGO.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = * ) :: file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: nc = 2
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  ZWEIGO groups data into TWO clusters.'

  call data_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  The number of data items is M =        ', m
  write ( *, '(a,i6)' ) '  The dimension of the data items is N = ', n
  write ( *, '(a)' ) ' '

  allocate ( p(1:m) )
  allocate ( x(1:m,1:n) )

  call data_d_read ( file_name, m, n, x )

  call zweigo ( m, n, x, p )

  write ( *, '(a,i6)' ) '  Number of clusters created was ', nc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I Cluster'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(2x,2i6)' ) i, p(i)
  end do

  deallocate ( p )
  deallocate ( x )

  return
end
subroutine test05 ( file_name )

!*****************************************************************************80
!
!! TEST05 tests HMEANS.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  HMEANS uses a variance diminishing procedure.'

  call data_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  The number of data items is M =        ', m
  write ( *, '(a,i6)' ) '  The dimension of the data items is N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For each cluster size from 1 to 10,'
  write ( *, '(a)' ) '  we carry out the experiment 3 times.'

  allocate ( p(1:m) )
  allocate ( x(1:m,1:n) )

  call data_d_read ( file_name, m, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NC   Initial       Final'
  write ( *, '(a)' ) '       Variance      Variance'
  write ( *, '(a)' ) ' '

  do nc = 1, 10

    allocate ( e(1:nc) )

    m0 = 1
    seed = 123456789

    write ( *, '(a)' ) ' '

    do k = 1, 3

      call randp ( m, m0, nc, p, seed )

      call cluster_variance ( m, n, x, nc, p, e )

      d1 = sum ( e(1:nc) )

      call hmeans ( m, n, x, nc, p )

      call cluster_variance ( m, n, x, nc, p, e )

      d2 = sum ( e(1:nc) )

      write ( *, '(2x,i2,2g14.6)' ) nc, d1, d2

    end do

    deallocate ( e )

  end do

  deallocate ( p )
  deallocate ( x )

  return
end
subroutine test06 ( file_name )

!*****************************************************************************80
!
!! TEST06 tests KMEANS.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) d1
  real ( kind = 8 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  KMEANS uses a variance diminishing procedure.'

  call data_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  The number of data items is M =        ', m
  write ( *, '(a,i6)' ) '  The dimension of the data items is N = ', n

  allocate ( p(1:m) )
  allocate ( x(1:m,1:n) )

  call data_d_read ( file_name, m, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NC   Initial       Final'
  write ( *, '(a)' ) '       Variance      Variance'
  write ( *, '(a)' ) ' '

  do nc = 1, 10

    allocate ( e(1:nc) )

    m0 = 1
    seed = 123456789

    write ( *, '(a)' ) ' '

    do k = 1, 3

      call randp ( m, m0, nc, p, seed )

      call cluster_variance ( m, n, x, nc, p, e )

      d1 = sum ( e(1:nc) )

      call kmeans ( m, n, x, nc, p, e, d )

      write ( *, '(2x,i2,2g14.6)' ) nc, d1, d

    end do

    deallocate ( e )

  end do

  deallocate ( p )
  deallocate ( x )

  return
end
subroutine test07 ( file_name )

!*****************************************************************************80
!
!! TEST07 tests CLUDIA.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  integer ( kind = 4 ) pi
  integer ( kind = 4 ) pj
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: t
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  CLUDIA carries out the KMEANS procedure for a'
  write ( *, '(a)' ) '    generalized distance function that is stored'
  write ( *, '(a)' ) '    in a matrix.  For this example, we simply wish'
  write ( *, '(a)' ) '    to repeat the particular case of standard distance'
  write ( *, '(a)' ) '    in the L2 norm.  But CLUDIA can handle more'
  write ( *, '(a)' ) '    interesting cases as well.'

  call data_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  The number of data items is M =        ', m
  write ( *, '(a,i6)' ) '  The dimension of the data items is N = ', n
  write ( *, '(a)' ) ' '

  allocate ( p(1:m) )
  allocate ( t(1:m,1:m) )
  allocate ( x(1:m,1:n) )

  call data_d_read ( file_name, m, n, x )
!
!  Determine the distance matrix.
!
  do i = 1, m
    do j = 1, m
      t(i,j) = sqrt ( sum ( ( x(i,1:n) - x(j,1:n) )**2 ) )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NC   Initial       Final'
  write ( *, '(a)' ) '       Variance      Variance'
  write ( *, '(a)' ) ' '

  do nc = 1, 10

    allocate ( e(1:nc) )

    m0 = 1
    seed = 123456789

    write ( *, '(a)' ) ' '

    do k = 1, 3

      call randp ( m, m0, nc, p, seed )

      e(1:nc) = 0.0E+00
      do i = 1, m
        pi = p(i)
        do j = i+1, m
          pj = p(j)
          if ( pj == pi ) then
            e(pj) = e(pj) + t(i,j)
          end if
        end do
      end do

      d1 = sum ( e(1:nc) )

      call cludia ( m, t, p, nc, e, d2 )

      write ( *, '(2x,i2,2g14.6)' ) nc, d1, d2

    end do

    deallocate ( e )

  end do

  deallocate ( p )
  deallocate ( t )
  deallocate ( x )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests R8MAT_DET
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) det

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  DMAT_DET computes the determinant of a'
  write ( *, '(a)' ) '  symmetric matrix.'

  call dif_inverse ( n, a )

  call r8mat_print ( n, n, a, '  The matrix to be analyzed:' )

  call r8mat_det ( n, a, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The computed determinant is ', det

  det = 1.0E+00 / real ( n + 1, kind = 8 )
  write ( *, '(a,g14.6)' ) '  The correct determinant is  ', det

  return
end
subroutine test09 ( file_name )

!*****************************************************************************80
!
!! TEST09 tests WMEANS.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) det
  character ( len = * ) :: file_name
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  WMEANS uses a determinant diminishing procedure.'

  call data_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  The number of data items is M =        ', m
  write ( *, '(a,i6)' ) '  The dimension of the data items is N = ', n
  write ( *, '(a)' ) ' '

  allocate ( p(1:m) )
  allocate ( x(1:m,1:n) )

  call data_d_read ( file_name, m, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NC   Determinant'
  write ( *, '(a)' ) ' '

  do nc = 1, 10

    m0 = 1
    seed = 123456789

    write ( *, '(a)' ) ' '

    do k = 1, 3

      call randp ( m, m0, nc, p, seed )

      call wmeans ( m, n, x, nc, p, det )

      write ( *, '(2x,i2,g14.6)' ) nc, det

    end do

  end do

  deallocate ( p )
  deallocate ( x )

  return
end
subroutine test10 ( file_name )

!*****************************************************************************80
!
!! TEST10 tests ORDERED.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nc = 5

  character ( len = * ) :: file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ), dimension ( nc, nc ) :: q
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: s
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  ORDERED clusters an ordered set of 1D data.'

  call data_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )

  if ( n /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST10 - Fatal error!'
    write ( *, '(a)' ) '  Input dataset is not 1 dimensional.'
    stop
  end if

  allocate ( s(1:m,1:nc) )
  allocate ( x(1:m) )

  call data_d_read ( file_name, m, n, x )

  call r8vec_sort_bubble_a ( m, x )

  call data_d_print ( m, n, x, '  The sorted 1D dataset:' )

  call ordered ( m, x, nc, q, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cluster  First'
  write ( *, '(a)' ) '           Value'
  write ( *, '(a)' ) ' '

  do i = 1, nc
    write ( *, '(2x,i6,i6)' ) i, q(nc,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Point  Variance'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(2x,i6,g14.6)' ) i, s(i,nc)
  end do

  deallocate ( s )
  deallocate ( x )

  return
end
subroutine test11 ( file_name )

!*****************************************************************************80
!
!! TEST11 tests CLUSTA.
!
!  Discussion:
!
!    This test is failing.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) d1
  real ( kind = 8 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: s
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  CLUSTA uses a multiple location allocation approach.'

  call data_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  The number of data items is M =        ', m
  write ( *, '(a,i6)' ) '  The dimension of the data items is N = ', n
  write ( *, '(a)' ) ' '

  allocate ( p(1:m) )
  allocate ( w(1:m) )
  allocate ( x(1:m,1:n) )

  call data_d_read ( file_name, m, n, x )
!
!  Use uniform weights.
!
  w(1:m) = 1.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NC   Initial       Final'
  write ( *, '(a)' ) '       Variance      Variance'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '

  do nc = 1, 10

    allocate ( e(1:nc) )
    allocate ( s(1:nc,1:n) )

    m0 = 1
    seed = 123456789

    write ( *, '(a)' ) ' '

    do k = 1, 3

      call randp ( m, m0, nc, p, seed )

      call cluster_variance ( m, n, x, nc, p, e )

      d1 = sum ( e(1:nc) )

      call clusta ( m, n, x, w, p, nc, s, e, d, iflag )

      write ( *, '(2x,i2,2g14.6)' ) nc, d1, d

    end do

    deallocate ( e )
    deallocate ( s )

  end do

  deallocate ( p )
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests I4VEC_PERML.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  logical first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) q(n)
  integer ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  I4VEC_PERML returns permutations of a vector.'

  do i = 1, n
    x(i) = i
  end do

  k = 0
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i2,2x,5i2)' ) k, x(1:n)

  first = .true.

  do

    call i4vec_perml ( n, x, q, first )

    k = k + 1
    write ( *, '(2x,i2,2x,5i2)' ) k, x(1:n)

    if ( first ) then
      exit
    end if

  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests I4VEC_PERMS.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  logical first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  I4VEC_PERMS returns the first half of the set of all'
  write ( *, '(a)' ) '  permutations of a vector.'

  do i = 1, n
    x(i) = i
  end do

  k = 0
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i2,2x,5i2)' ) k, x(1:n)

  first = .true.

  do

    call i4vec_perms ( n, x, first )

    k = k + 1
    write ( *, '(2x,i2,2x,5i2)' ) k, x(1:n)

    if ( first ) then
      exit
    end if

  end do

  return
end
subroutine test14 ( file_name )

!*****************************************************************************80
!
!! TEST14 tests KMEANS.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ), allocatable, dimension ( : ) :: e
  character ( len = * ) :: file_name
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nc
  integer ( kind = 4 ), allocatable, dimension ( : ) :: p
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  EMEANS uses a sum of absolute distances to the median.'

  call data_size ( file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data set is in file ' // trim ( file_name )
  write ( *, '(a,i6)' ) '  The number of data items is M =        ', m
  write ( *, '(a,i6)' ) '  The dimension of the data items is N = ', n

  allocate ( p(1:m) )
  allocate ( x(1:m,1:n) )

  call data_d_read ( file_name, m, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NC   Initial       Final'
  write ( *, '(a)' ) '       Median        Median'
  write ( *, '(a)' ) '       Distance      Distance'
  write ( *, '(a)' ) ' '

  do nc = 1, 10

    allocate ( e(1:nc) )

    m0 = 1
    seed = 123456789

    write ( *, '(a)' ) ' '

    do k = 1, 3

      call randp ( m, m0, nc, p, seed )

      call cluster_median_distance ( m, n, x, nc, p, d1 )

      call emeans ( m, n, x, nc, p )

      call cluster_median_distance ( m, n, x, nc, p, d2 )

      write ( *, '(2x,i2,2x,g14.6,2x,g14.6)' ) nc, d1, d2

    end do

    deallocate ( e )

  end do

  deallocate ( p )
  deallocate ( x )

  return
end
