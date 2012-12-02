program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA058_PRB.
!
!  Discussion:
!
!    ASA058_PRB tests the ASA058 clustering algorithm.
!
!  Modified:
!
!    23 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA058_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA058 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA058_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tries out the ASA058 routine.
!
!  Modified:
!
!    04 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 5
  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n = 100

  integer ( kind = 4 ) b(n)
  real ( kind = 8 ) d(k,m)
  real ( kind = 8 ) dev(k)
  real ( kind = 8 ) dev_sum
  integer ( kind = 4 ) e(k)
  integer ( kind = 4 ) e_sum
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) nz
  real ( kind = 8 ) x(n,m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test the CLUSTR algorithm.'
  write ( *, '(a)' ) '  Applied Statistics Algorithm 58'
!
!  Read the data.
!
  write (  *, '(a)' ) ' '
  write (  *, '(a)' ) '  Reading the data.'

  open ( unit = 1, file = 'points_100.txt', status = 'old' )

  do i = 1, n
    read ( 1, * ) ( x(i,j), j = 1, m )
  end do

  close ( unit = 1 )
!
!  Print a few data values.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 5 data values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, ( x(i,j), j = 1, m )
  end do
!
!  Initialize the cluster centers arbitrarily.
!
  do i = 1, k
    do j = 1, m
      d(i,j) = x(i,j)
    end do
  end do
!
!  Compute the clusters.
!
  nz = 1
  k2 = k

  call clustr ( x, d, dev, b, f, e, n, m, k, nz, k2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cluster  Population  Energy'
  write ( *, '(a)' ) ' '

  do i = 1, k
    write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) i, e(i), dev(i)
  end do

  e_sum = 0
  dev_sum = 0.0D+00

  do i = 1, k
    e_sum = e_sum + e(i)
    dev_sum = dev_sum + dev(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(2x,a8,2x,i8,2x,g14.6)' ) '   Total', e_sum, dev_sum

  return
end
