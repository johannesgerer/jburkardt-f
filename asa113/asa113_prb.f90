program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA113_PRB.
!
!  Discussion:
!
!    ASA113_PRB tests the ASA113 clustering algorithm.
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA113_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA113 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA113_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tries out the ASA113 routine.
!
!  Modified:
!
!   16 February 2008
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
  integer ( kind = 4 ) c(m)
  real ( kind = 8 ) c_center(k,n)
  integer ( kind = 4 ) c_size(k)
  integer ( kind = 4 ) ci
  real ( kind = 8 ) critvl
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ntrans1
  integer ( kind = 4 ) ntrans2
  real ( kind = 8 ) wss(k)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test the ASA113 classification algorithm.'
!
!  Read the data
!
  open ( unit = 1, file = 'points_100.txt', status = 'old' )

  do i = 1, m
    read ( 1, * ) ( a(i,j), j = 1, n )
  end do

  close ( unit = 1 )
!
!  Print first five points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First five points:'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, ( a(i,j), j = 1, n )
  end do
!
!  Assign points randomly to classes.
!
  do i = 1, m
    c(i) = mod ( i, k ) + 1
  end do
!
!  Define the critical value as the sum of the squares of the distances
!  of the points to their cluster center.
!
  do i = 1, k
    c_size(i) = 0
    do j = 1, n
      c_center(i,j) = 0.0D+00
    end do
  end do

  do i = 1, m
    ci = c(i)
    c_size(ci) = c_size(ci) + 1
    do j = 1, n
      c_center(ci,j) = c_center(ci,j) + a(i,j)
    end do
  end do

  do i = 1, k
    do j = 1, n
      c_center(i,j) = c_center(i,j) / real ( c_size(i), kind = 8 )
    end do
  end do

  do i = 1, k
    wss(i) = 0.0D+00
  end do

  do i = 1, m
    ci = c(i)
    do j = 1, n
      wss(ci) = wss(ci) + ( a(i,j) - c_center(ci,j) )**2
    end do
  end do

  critvl = 0.0D+00
  do i = 1, k
    critvl = critvl + wss(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '        Initial CRITVL = ', critvl
!
!  Compute the clusters.
!
  ntrans1 = -1
  ntrans2 = -1

  do

    call trnsfr ( a, c, c_size, m, k, n, critvl, ntrans1, ifault )

    if ( ntrans1 == 0 .and. ntrans2 == 0 ) then
      exit
    end if

    write ( *, '(a,g14.6)' ) '  After TRNSFR, CRITVL = ', critvl

    call swap ( a, c, c_size, m, k, n, critvl, ntrans2, ifault )

    if ( ntrans1 == 0 .and. ntrans2 == 0 ) then
      exit
    end if

    write ( *, '(a,g14.6)' ) '    After SWAP, CRITVL = ', critvl

  end do
!
!  Define the critical value as the sum of the squares of the distances
!  of the points to their cluster center.
!
  do i = 1, k
    do j = 1, n
      c_center(i,j) = 0.0D+00
    end do
  end do

  do i = 1, m
    ci = c(i)
    do j = 1, n
      c_center(ci,j) = c_center(ci,j) + a(i,j)
    end do
  end do

  do i = 1, k
    do j = 1, n
      c_center(i,j) = c_center(i,j) / real ( c_size(i), kind = 8 )
    end do
  end do

  do i = 1, k
    wss(i) = 0.0D+00
  end do

  do i = 1, m
    ci = c(i)
    do j = 1, n
      wss(ci) = wss(ci) + ( a(i,j) - c_center(ci,j) )**2
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cluster  Population  Energy'
  write ( *, '(a)' ) ' '

  do i = 1, k
    write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) i, c_size(i), wss(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a8,2x,i8,2x,g14.6)' ) '   Total', m, critvl

  return
end
subroutine crswap ( a, c, c_size, m, k, n, critvl, i1, i2, c1, c2, iswitch, &
  inc )

!*****************************************************************************80
!
!! CRSWAP determines the effect of swapping two objects.
!
!  Discussion:
!
!    This computation is very inefficient.  It is only set up so that we
!    can compare algorithm ASA 113 to the K-means algorithms ASA 058 and
!    ASA 136.
!
!  Modified:
!
!    15 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Colin Banfield, LC Bassill,
!    Algorithm AS 113:
!    A transfer for non-hierarchichal classification,
!    Applied Statistics,
!    Volume 26, Number 2, 1977, pages 206-210.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(M,N), the data values.  There are M objects,
!    each having spatial dimension N.
!
!    Input, integer ( kind = 4 ) C(M), the classification of each object.
!
!    Input, integer ( kind = 4 ) C_SIZE(K), the number of objects in each class.
!
!    Input, integer ( kind = 4 ) M, the number of objects.
!
!    Input, integer ( kind = 4 ) K, the number of classes.
!
!    Input, integer ( kind = 4 ) N, the number of spatial dimensions, or variates,
!    of the objects.
!
!    Input, real ( kind = 8 ) CRITVL, the current value of the criterion.
!
!    Input, integer ( kind = 4 ) I1, I2, the objects to be swapped.
!
!    Input, integer ( kind = 4 ) C1, C2, the current classes of objects I1 and I2.
!
!    Input, integer ( kind = 4 ) ISWITCH:
!    1, indicates that I1 and I2 should be temporarily swapped, the
!       change in CRITVL should be computed, and then I1 and I2 restored.
!    2, indicates that I1 and I2 will be swapped.
!
!    Output, real ( kind = 8 ) INC, the change to CRITVL that would occur if I1 and
!    I2 were swapped.  This is only computed for ISWITCH = 1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) c(m)
  real ( kind = 8 ) c_center(k,n)
  integer ( kind = 4 ) c_size(k)
  integer ( kind = 4 ) c1
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) ci
  real ( kind = 8 ) critvl
  real ( kind = 8 ) critvl_new
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real ( kind = 8 ) inc
  integer ( kind = 4 ) iswitch
  integer ( kind = 4 ) j
  real ( kind = 8 ) wss(k)

  if ( iswitch == 2 ) then
    return
  end if
!
!  Move object I1 from class C1 to class C2.
!  Move object I2 from class C2 to class C1.
!
  c(i1) = c2
  c(i2) = c1
!
!  Define the critical value as the sum of the squares of the distances
!  of the points to their cluster center.
!
  do i = 1, k
    c_size(i) = 0
    do j = 1, n
      c_center(i,j) = 0.0D+00
    end do
  end do

  do i = 1, m
    ci = c(i)
    c_size(ci) = c_size(ci) + 1
    do j = 1, n
      c_center(ci,j) = c_center(ci,j) + a(i,j)
    end do
  end do

  do i = 1, k
    do j = 1, n
      c_center(i,j) = c_center(i,j) / real ( c_size(i), kind = 8 )
    end do
  end do

  do i = 1, k
    wss(i) = 0.0D+00
  end do

  do i = 1, m
    ci = c(i)
    do j = 1, n
      wss(ci) = wss(ci) + ( a(i,j) - c_center(ci,j) )**2
    end do
  end do

  critvl_new = 0.0D+00
  do i = 1, k
    critvl_new = critvl_new + wss(i)
  end do

  inc = critvl_new - critvl
!
!  Move object I1 from class C2 to class C1.
!  Move object I2 from class C1 to class C2.
!
  c(i1) = c1
  c(i2) = c2

  return
end
subroutine crtran ( a, c, c_size, m, k, n, critvl, i1, c1, c2, iswitch, inc )

!*****************************************************************************80
!
!! CRTRAN determines the effect of moving an object to another class.
!
!  Discussion:
!
!    This computation is very inefficient.  It is only set up so that we
!    can compare algorithm ASA 113 to the K-means algorithms ASA 058 and
!    ASA 136.
!
!  Modified:
!
!    15 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Colin Banfield, LC Bassill,
!    Algorithm AS 113:
!    A transfer for non-hierarchichal classification,
!    Applied Statistics,
!    Volume 26, Number 2, 1977, pages 206-210.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(M,N), the data values.  There are M objects,
!    each having spatial dimension N.
!
!    Input, integer ( kind = 4 ) C(M), the classification of each object.
!
!    Input, integer ( kind = 4 ) C_SIZE(K), the number of objects in each class.
!
!    Input, integer ( kind = 4 ) M, the number of objects.
!
!    Input, integer ( kind = 4 ) K, the number of classes.
!
!    Input, integer ( kind = 4 ) N, the number of spatial dimensions, or variates,
!    of the objects.
!
!    Input, real ( kind = 8 ) CRITVL, the current value of the criterion.
!
!    Input, integer ( kind = 4 ) I1, the object to be transferred.
!
!    Input, integer ( kind = 4 ) C1, C2, the current class of object I1, and the
!    class to which it may be transferred.
!
!    Input, integer ( kind = 4 ) ISWITCH:
!    1, indicates that I1 should be temporarily transferred, the
!       change in CRITVL should be computed, and then I1 restored.
!    2, indicates that I1 will be permanently transferred.
!
!    Output, real ( kind = 8 ) INC, the change to CRITVL that would occur if I1 were
!    transferred from class C1 to C2.  This is only computed for ISWITCH = 1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) c(m)
  real ( kind = 8 ) c_center(k,n)
  integer ( kind = 4 ) c_size(k)
  integer ( kind = 4 ) c1
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) ci
  real ( kind = 8 ) critvl
  real ( kind = 8 ) critvl_new
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  real ( kind = 8 ) inc
  integer ( kind = 4 ) iswitch
  integer ( kind = 4 ) j
  real ( kind = 8 ) wss(k)

  if ( iswitch == 2 ) then
    return
  end if
!
!  Move object I from class C1 to class C2.
!
  c(i1) = c2
  c_size(c1) = c_size(c1) - 1
  c_size(c2) = c_size(c2) + 1
!
!  Define the critical value as the sum of the squares of the distances
!  of the points to their cluster center.
!
  do i = 1, k
    c_size(i) = 0
    do j = 1, n
      c_center(i,j) = 0.0D+00
    end do
  end do

  do i = 1, m
    ci = c(i)
    c_size(ci) = c_size(ci) + 1
    do j = 1, n
      c_center(ci,j) = c_center(ci,j) + a(i,j)
    end do
  end do

  do i = 1, k
    do j = 1, n
      c_center(i,j) = c_center(i,j) / real ( c_size(i), kind = 8 )
    end do
  end do

  do i = 1, k
    wss(i) = 0.0D+00
  end do

  do i = 1, m
    ci = c(i)
    do j = 1, n
      wss(ci) = wss(ci) + ( a(i,j) - c_center(ci,j) )**2
    end do
  end do

  critvl_new = 0.0D+00
  do i = 1, k
    critvl_new = critvl_new + wss(i)
  end do

  inc = critvl_new - critvl
!
!  Move object I1 from class C2 to class C1.
!
  c(i1) = c1
  c_size(c1) = c_size(c1) + 1
  c_size(c2) = c_size(c2) - 1

  return
end
