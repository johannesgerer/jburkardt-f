program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA136_PRB.
!
!  Discussion:
!
!    ASA136_PRB tests the ASA136 clustering algorithm.
!
!  Modified:
!
!    14 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA136_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA136 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA136_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ()

!*****************************************************************************80
!
!! TEST01 tries out the ASA136 routine.
!
!  Modified:
!
!    14 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 5
  integer ( kind = 4 ), parameter :: m = 100
  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) an1(k)
  real ( kind = 8 ) an2(k)
  real ( kind = 8 ) c(k,n)
  real ( kind = 8 ) d(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic1(m)
  integer ( kind = 4 ) ic2(m)
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) itran(k)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) live(k)
  integer ( kind = 4 ) nc(k)
  integer ( kind = 4 ) nc_sum
  integer ( kind = 4 ) ncp(k)
  real ( kind = 8 ) wss(k)
  real ( kind = 8 ) wss_sum

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test the KMNS algorithm,'
  write ( *, '(a)' ) '  Applied Statistics Algorithm #136.'
!
!  Read the data.
!
  open ( unit = 1, file = 'points_100.txt', status = 'old' )

  do i = 1, m
    read ( 1, * ) ( a(i,j), j = 1, n )
  end do

  close ( unit = 1 )
!
!  Print a few data values.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 5 data values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, ( a(i,j), j = 1, n )
  end do
!
!  Initialize the cluster centers.
!  Here, we arbitrarily make the first K data points cluster centers.
!
  do i = 1, k
    do j = 1, n
      c(i,j) = a(i,j)
    end do
  end do

  iter = 50
!
!  Compute the clusters.
!
  call kmns ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, &
    itran, live, iter, wss, ifault )

  if ( ifault /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a,i8)' ) '  KMNS returned IFAULT = ', ifault
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cluster  Population  Energy'
  write ( *, '(a)' ) ' '

  nc_sum = 0
  wss_sum = 0.0D+00

  do i = 1, k
    write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) i, nc(i), wss(i)
    nc_sum = nc_sum + nc(i)
    wss_sum = wss_sum + wss(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(2x,a8,2x,i8,2x,g14.6)' ) '   Total', nc_sum, wss_sum

  return
end
