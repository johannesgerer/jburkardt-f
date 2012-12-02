program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS526_PRB.
!
!  Discussion:
!
!    TOMS526_PRB runs tests on the TOMS526 bivariate scattered data routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS526_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS526 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS526_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests IDBVIP, interpolation at scattered output points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ncp = 4
  integer ( kind = 4 ), parameter :: ndp = 20
  integer ( kind = 4 ), parameter :: nip = 5

  real ( kind = 8 ) error
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwk(max(31,27+ncp)*ndp+nip)
  integer ( kind = 4 ) md
  real ( kind = 8 ) temp
  real ( kind = 8 ) wk(8*ndp)
  real ( kind = 8 ), dimension ( ndp ) :: xd = (/ &
     0.632000D+00, 0.295109D+00, 0.578769D+00, 0.140772D+00, 0.985549D+00, &
     0.209717D+00, 0.893022D+00, 0.928088D+00, 0.766803D+00, 0.819351D+00, &
     0.959191D+00, 0.516312D+00, 0.052801D+00, 0.626174D+00, 0.275320D+00, &
     0.659645D+00, 0.396438D+00, 0.204824D+00, 0.906341D+00, 0.880356D+00 /)
  real ( kind = 8 ) xi(nip)
  real ( kind = 8 ), dimension ( ndp ) :: yd = (/ &
     0.918675D+00, 0.821677D+00, 0.962992D+00, 0.477774D+00, 0.056257D+00, &
     0.654656D+00, 0.377358D+00, 0.263626D+00, 0.706323D+00, 0.584658D+00, & 
     0.154619D+00, 0.986449D+00, 0.197552D+00, 0.924616D+00, 0.788099D+00, &
     0.886828D+00, 0.945136D+00, 0.643341D+00, 0.335303D+00, 0.416050D+00 /)
  real ( kind = 8 ) yi(nip)
  real ( kind = 8 ) zd(ndp)
  real ( kind = 8 ) zi(nip)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  IDBVIP accepts scattered Z(X,Y) data, and interpolates '
  write ( *, '(a)' ) '  values of Z at new points.'
!
!  Set up the Z data.
!
  zd(1:ndp) = exp ( xd(1:ndp) ) * sin ( yd(1:ndp) )
!
!  Set the points at which we want interpolated values.
!
  do i = 1, nip
    xi(i) = real ( i - 1, kind = 8 ) / real ( nip - 1, kind = 8 )
  end do

  yi(1:nip) = 1.0D+00 - xi(1:nip)
!
!  Call IDBVIP to estimate the value of Z at the points (XI,YI).
!
  md = 1

  call idbvip ( md, ncp, ndp, xd, yd, zd, nip, xi, yi, zi, iwk, wk )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   X            Y             Z             Exact Z       Error'
  write ( *, '(a)' ) ' '

  do i = 1, nip
    temp = exp ( xi(i) ) * sin ( yi(i) )
    error = zi(i) - temp
    write ( *, '(5g14.6)' ) xi(i), yi(i), zi(i), temp, error
  end do
 
  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests IDSFFT, interpolation on a regular grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ncp = 4
  integer ( kind = 4 ), parameter :: ndp = 20
  integer ( kind = 4 ), parameter :: nxi = 5
  integer ( kind = 4 ), parameter :: nyi = 4

  real ( kind = 8 ) error
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwk(max(31,27+ncp)*ndp+nxi*nyi)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) md
  real ( kind = 8 ) temp
  real ( kind = 8 ) wk(5*ndp)
  real ( kind = 8 ), dimension ( ndp ) :: xd = (/ &
    0.632000D+00, 0.295109D+00, 0.578769D+00, 0.140772D+00, 0.985549D+00, &
    0.209717D+00, 0.893022D+00, 0.928088D+00, 0.766803D+00, 0.819351D+00, &
    0.959191D+00, 0.516312D+00, 0.052801D+00, 0.626174D+00, 0.275320D+00, &
    0.659645D+00, 0.396438D+00, 0.204824D+00, 0.906341D+00, 0.880356D+00 /)
  real ( kind = 8 ) xi(nxi)
  real ( kind = 8 ), dimension ( ndp ) :: yd = (/ &
      0.918675D+00, 0.821677D+00, 0.962992D+00, 0.477774D+00, 0.056257D+00, &
      0.654656D+00, 0.377358D+00, 0.263626D+00, 0.706323D+00, 0.584658D+00, &
      0.154619D+00, 0.986449D+00, 0.197552D+00, 0.924616D+00, 0.788099D+00, &
      0.886828D+00, 0.945136D+00, 0.643341D+00, 0.335303D+00, 0.416050D+00 /)
  real ( kind = 8 ) yi(nyi)
  real ( kind = 8 ) zd(ndp)
  real ( kind = 8 ) zi(nxi,nyi)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  IDSFFT accepts scattered Z(X,Y) data, and interpolates'
  write ( *, '(a)' ) '  values of Z on a regular grid of points.'
!
!  Set up the Z data
!
  zd(1:ndp) = exp ( xd(1:ndp) ) * sin ( yd(1:ndp) )
!
!  Set up the points at which we want interpolated values.
!
  do i = 1, nxi
    xi(i) = real ( i - 1, kind = 8 ) / real ( nxi - 1, kind = 8 )
  end do
 
  do j = 1, nyi
    yi(j) = real ( j - 1, kind = 8 ) / real ( nyi - 1, kind = 8 )
  end do
!
!  Call IDSFFT to estimate the value of Z.
!
  md = 1
 
  call idsfft ( md, ncp, ndp, xd, yd, zd, nxi, nyi, xi, yi, zi, iwk, wk )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   X            Y             Z             Exact Z       Error'
  write ( *, '(a)' ) ' '

  do i = 1, nxi
    do j = 1, nyi
      temp = exp ( xi(i) ) * sin ( yi(j) )
      error = zi(i,j) - temp
      write ( *, '(5g14.6)' ) xi(i), yi(j), zi(i,j), temp, error
    end do
    write ( *, '(a)' ) ' '
  end do
 
  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 is the test code supplied with the original algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ncp = 4
  integer ( kind = 4 ), parameter :: ndp = 30
  integer ( kind = 4 ), parameter :: nip = 1
  integer ( kind = 4 ), parameter :: nxi = 6
  integer ( kind = 4 ), parameter :: nyi = 5

  real ( kind = 8 ), dimension ( nxi, nyi ) :: dzi1
  real ( kind = 8 ), dimension ( nxi, nyi ) :: dzi2
  integer ( kind = 4 ) idp
  integer ( kind = 4 ) iwk1(max(31,27+ncp)*ndp+nip)
  integer ( kind = 4 ) iwk2(max(31,27+ncp)*ndp+nxi*nyi)
  integer ( kind = 4 ) ixi
  integer ( kind = 4 ) iyi
  integer ( kind = 4 ) md
  real ( kind = 8 ) wk(8*ndp)
  real ( kind = 8 ), dimension ( ndp ) :: xd = (/ &
    11.16D+00, 24.20D+00, 19.85D+00, 10.35D+00, 19.72D+00, &
     0.00D+00, 20.87D+00, 19.99D+00, 10.28D+00,  4.51D+00, &
     0.00D+00, 16.70D+00,  6.08D+00, 25.00D+00, 14.90D+00, &
     0.00D+00,  9.66D+00,  5.22D+00, 11.77D+00, 15.10D+00, &
    25.00D+00, 25.00D+00, 14.59D+00, 15.20D+00,  5.23D+00, &
     2.14D+00,  0.51D+00, 25.00D+00, 21.67D+00,  3.31D+00 /)
  real ( kind = 8 ), dimension ( nxi ) :: xi = (/ &
    0.00D+00,  5.00D+00, 10.00D+00, 15.00D+00, 20.00D+00, 25.00D+00 /)
  real ( kind = 8 ), dimension ( ndp ) :: yd = (/ &
     1.24D+00, 16.23D+00, 10.72D+00,  4.11D+00,  1.39D+00, &
    20.00D+00, 20.00D+00,  4.62D+00, 15.16D+00, 20.00D+00, &
     4.48D+00, 19.65D+00,  4.58D+00, 11.87D+00,  3.12D+00, &
     0.00D+00, 20.00D+00, 14.66D+00, 10.47D+00, 17.19D+00, &
     3.87D+00,  0.00D+00,  8.71D+00,  0.00D+00, 10.72D+00, &
    15.03D+00,  8.37D+00, 20.00D+00, 14.36D+00,  0.13D+00 /)
  real ( kind = 8 ), dimension ( nyi ) :: yi = (/ &
    0.00D+00,  5.00D+00, 10.00D+00, 15.00D+00, 20.00D+00 /)
  real ( kind = 8 ), dimension ( ndp ) :: zd = (/ &
    22.15D+00,  2.83D+00,  7.97D+00, 22.33D+00, 16.83D+00, &
    34.60D+00,  5.74D+00, 14.72D+00, 21.59D+00, 15.61D+00, &
    61.77D+00,  6.31D+00, 35.74D+00,  4.40D+00, 21.70D+00, &
    58.20D+00,  4.73D+00, 40.36D+00, 13.62D+00, 12.57D+00, &
     8.74D+00, 12.00D+00, 14.81D+00, 21.60D+00, 26.50D+00, &
    53.10D+00, 49.43D+00,  0.60D+00,  5.52D+00, 44.08D+00 /)
  real ( kind = 8 ), dimension ( nxi, nyi ) :: zi = reshape ( (/ &
    58.20D+00, 39.55D+00, 26.90D+00, 21.71D+00, 17.68D+00, 12.00D+00, &
    61.58D+00, 39.39D+00, 22.04D+00, 21.29D+00, 14.36D+00,  8.04D+00, &
    59.18D+00, 27.39D+00, 16.78D+00, 13.25D+00,  8.59D+00,  5.36D+00, &
    52.82D+00, 40.27D+00, 22.76D+00, 16.61D+00,  7.40D+00,  2.88D+00, &
    34.60D+00, 14.05D+00,  4.12D+00,  3.17D+00,  6.31D+00,  0.60D+00 /), &
    (/ nxi, nyi /) )
  real ( kind = 8 ), dimension ( nxi, nyi ) :: zi1
  real ( kind = 8 ), dimension ( nxi, nyi ) :: zi2  

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Original test program for ACM TOMS 526'
!
!  Print the input data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i3)' ) '   Input data        NDP =', ndp
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      I      XD     YD     ZD'
  write ( *, '(a)' ) ' '

  do idp = 1, ndp
    if ( mod ( idp, 5 ) == 1 ) then
      write ( *, '(a)' ) ' '
    end if
    write ( *, '(5x,i2,2x,3f7.2)' ) idp, xd(idp), yd(idp), zd(idp)
  end do
!
!  IDBVIP calculation.
!
  md = 1

  do iyi = 1, nyi
    do ixi = 1, nxi
      
      call idbvip ( md, ncp, ndp, xd, yd, zd, nip, xi(ixi), yi(iyi), zi1(ixi,iyi), iwk1, wk )

      md = 2
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Using IDBVIP for interpolation on arbitrary points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                          ZI1(XI,YI)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       XI    YI ='
  write ( *, '(12x,5f7.2)' ) yi(1:nyi)
  write ( *, '(a)' ) ' '
  do ixi = 1, nxi
    write ( *, '(1x,f9.2,2x,5f7.2)' ) xi(ixi), zi1(ixi,1:nyi)
  end do

  dzi1(1:nxi,1:nyi) = abs ( zi1(1:nxi,1:nyi) - zi(1:nxi,1:nyi) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Difference'
  write ( *, '(a)' ) '                         DZI1(XI,YI)'
  write ( *, '(a)' ) '       XI    YI ='
  write ( *, '(12x,5f7.2)' ) yi(1:nyi)
  write ( *, '(a)' ) ' '
  do ixi = 1, nxi
    write ( *, '(1x,f9.2,2x,5f7.2)' ) xi(ixi), dzi1(ixi,1:nyi)
  end do
!
!  IDSFFT calculation.
!
  md = 1

  call idsfft ( md, ncp, ndp, xd, yd, zd, nxi, nyi, xi, yi, zi2, iwk2, wk )

  dzi2(1:nxi,1:nyi) = abs ( zi2(1:nxi,1:nyi) - zi(1:nxi,1:nyi) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Using IDSFFT for interpolation on a regular grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                          ZI2(XI,YI)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       XI    YI ='
  write ( *, '(12x,5f7.2)' ) yi(1:nyi)
  write ( *, '(a)' ) ' '
  do ixi = 1, nxi
    write ( *, '(1x,f9.2,2x,5f7.2)' ) xi(ixi), zi2(ixi,1:nyi)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Difference'
  write ( *, '(a)' ) '                         DZI2(XI,YI)'
  write ( *, '(a)' ) '       XI    YI ='
  write ( *, '(12x,5f7.2)' ) yi(1:nyi)
  write ( *, '(a)' ) ' '
  do ixi = 1, nxi
    write ( *, '(1x,f9.2,2x,5f7.2)' ) xi(ixi), dzi2(ixi,1:nyi)
  end do
    
  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2005
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
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
