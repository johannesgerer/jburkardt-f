subroutine idbvip ( md, ncp, ndp, xd, yd, zd, nip, xi, yi, zi, iwk, wk )

!*****************************************************************************80
!
!! IDBVIP performs bivariate interpolation of irregular X, Y data.
!
!  Discussion:
!
!    This subroutine performs bivariate interpolation when the projections
!    of the data points in the XY plane are irregularly
!    distributed in the plane.
!
!    The very first call to this subroutine and the call with a new
!    NCP value, a new NDP value, or new contents of the XD and
!    YD arrays must be made with MD = 1.  The call with MD = 2 must be
!    preceded by another call with the same NCP and NDP values and
!    with the same contents of the XD and YD arrays.  The call with
!    MD = 3 must be preceded by another call with the same NCP, NDP,
!    and NIP values and with the same contents of the XD, YD, XI,
!    and YI arrays.  Between the call with MD = 2 or MD = 3 and its
!    preceding call, the IWK and WK arrays must not be disturbed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    Original FORTRAN77 version by Hiroshi Akima.
!    FORTRAN90 version by John Burkardt.
!
!  References:
!
!    Hiroshi Akima,
!    Algorithm 526:
!    A Method of Bivariate Interpolation and Smooth Surface Fitting
!    for Values Given at Irregularly Distributed Points,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, pages 160-164.
!
!    Hiroshi Akima,
!    On Estimating Partial Derivatives for Bivariate Interpolation
!    of Scattered Data,
!    Rocky Mountain Journal of Mathematics,
!    Volume 14, Number 1, Winter 1984, pages 41-51.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MD, mode of computation.  MD must be 1, 2 or 3.
!    1 for new NCP and/or new XD, YD;
!    2 for old NCP, old XD, YD, new XI, YI;
!    3 for old NCP, old XD, YD, old XI, YI.
!
!    Input, integer ( kind = 4 ) NCP, the number of additional data points used
!    for estimating partial derivatives at each data point.
!    2 < =  NCP < NDP.  Reasonable values are 3 < =  NCP < =  5.
!
!    Input, integer ( kind = 4 ) NDP, the number of data points.
!    4 < =  NDP.
!
!    Input, real ( kind = 8 ) XD(NDP), Y(NDP), the X and Y coordinates
!    of the data points.
!
!    Input, real ( kind = 8 ) ZD(NDP), the data values at the data points.
!
!    Input, integer ( kind = 4 ) NIP, the number of output points at which
!    interpolation is to be performed.
!    1 < =  NIP.
!
!    Input, real ( kind = 8 ) XI(NIP), YI(NIP), the coordinates of the
!    points at which interpolation is to be performed.
!
!    Output, real ( kind = 8 ) ZI(NIP), the interpolated data values.
!
!    Workspace, integer IWK(max(31,27+NCP)*NDP+NIP).
!
!    Workspace, real ( kind = 8 ) WK(8*NDP).
!
  implicit none

  integer ( kind = 4 ) ncp
  integer ( kind = 4 ) ndp
  integer ( kind = 4 ) nip

  integer ( kind = 4 ) iip
  integer ( kind = 4 ) itpv
  integer ( kind = 4 ) iwk(max(31,27+ncp)*ndp+nip)
  integer ( kind = 4 ) jwipc
  integer ( kind = 4 ) jwipl
  integer ( kind = 4 ) jwipt
  integer ( kind = 4 ) jwit
  integer ( kind = 4 ) jwit0
  integer ( kind = 4 ) jwiwk
  integer ( kind = 4 ) jwiwl
  integer ( kind = 4 ) jwiwp
  integer ( kind = 4 ) md
  integer ( kind = 4 ) md0
  integer ( kind = 4 ) ncp0
  integer ( kind = 4 ) ncppv
  integer ( kind = 4 ) ndp0
  integer ( kind = 4 ) ndppv
  integer ( kind = 4 ) nip0
  integer ( kind = 4 ) nippv
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nt
  real ( kind = 8 ) wk(8*ndp)
  real ( kind = 8 ) xd(ndp)
  real ( kind = 8 ) xi(nip)
  real ( kind = 8 ) yd(ndp)
  real ( kind = 8 ) yi(nip)
  real ( kind = 8 ) zd(ndp)
  real ( kind = 8 ) zi(nip)

  save /idlc/
  save /idpi/

  common /idlc/ nit
  common /idpi/ itpv
!
!  Setting of some input parameters to local variables,
!  for MD = 1, 2, 3.
!
  md0 = md
  ncp0 = ncp
  ndp0 = ndp
  nip0 = nip
!
!  Error check for MD = 1, 2, 3.
!
  if ( md0 < 1 .or. 3 < md0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IDBVIP - Fatal error!'
    write ( *, '(a)' ) '  Improper input parameters.'
    write ( *, '(a)' ) '  MD < 1 or 3 < MD.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  MD  = ', md0
    stop
  end if

  if ( ncp0 < 2 .or. ndp0 <= ncp0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IDBVIP - Fatal error!'
    write ( *, '(a)' ) '  Improper input parameters.'
    write ( *, '(a)' ) '  NCP < 2 or NDP <= NCP.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  NCP = ', ncp0
    write ( *, '(a,i8)' ) '  NDP = ', ndp0
    stop
  end if

  if ( ndp0 < 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IDBVIP - Fatal error!'
    write ( *, '(a)' ) '  Improper input parameters.'
    write ( *, '(a)' ) '  NDP < 4.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  NDP = ', ndp0
    stop
  end if

  if ( nip0 < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IDBVIP - Fatal error!'
    write ( *, '(a)' ) '  NIP < 1.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  NIP = ', nip0
    stop
  end if

  if ( md0 < 2 ) then

    iwk(1) = ncp0
    iwk(2) = ndp0

  else

    ncppv = iwk(1)
    ndppv = iwk(2)

    if ( ncp0 /= ncppv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'IDBVIP - Fatal error!'
      write ( *, '(a)' ) '  Improper input parameters.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  NCP   = ', ncp0
      write ( *, '(a,i8)' ) '  NCPPV = ', ncppv
      stop
    end if

    if ( ndp0 /= ndppv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'IDBVIP - Fatal error!'
      write ( *, '(a)' ) '  Improper input parameters.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  NDP   = ', ndp0
      write ( *, '(a,i8)' ) '  NDPPV = ', ndppv
      stop
    end if

  end if

  if ( md0 < 3 ) then

    iwk(3) = nip

  else

    nippv = iwk(3)

    if ( nip0 /= nippv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'IDBVIP - Fatal error!'
      write ( *, '(a)' ) '  Improper input parameters.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  NIP   = ', nip0
      write ( *, '(a,i8)' ) '  NIPPV = ', nippv
      stop
    end if

  end if
!
!  Allocation of storage areas in the IWK array, for MD = 1, 2, 3.
!
  jwipt = 16
  jwiwl = 6 * ndp0 + 1
  jwiwk = jwiwl
  jwipl = 24 * ndp0 + 1
  jwiwp = 30 * ndp0 + 1
  jwipc = 27 * ndp0 + 1
  jwit0 = max ( 31, 27 + ncp0 ) * ndp0
!
!  Triangulate the XY plane, for MD = 1.
!
  if ( md0 == 1 ) then

    call idtang ( ndp0, xd, yd, nt, iwk(jwipt), nl, iwk(jwipl), &
      iwk(jwiwl), iwk(jwiwp), wk )

    iwk(5) = nt
    iwk(6) = nl

    if ( nt == 0 ) then
      return
    end if

  end if
!
!  Determine NCP points closest to each data point, for MD = 1.
!
  if ( md0 <= 1 ) then

    call idcldp ( ndp0, xd, yd, ncp0, iwk(jwipc) )

    if ( iwk(jwipc) == 0 ) then
      return
    end if

  end if
!
!  Locate all points at which interpolation is to be performed,
!  for MD = 1, 2.
!
  if ( md0 /= 3 ) then

    nit = 0
    jwit = jwit0

    do iip = 1, nip0
      jwit = jwit + 1
      call idlctn ( ndp0, xd, yd, nt, iwk(jwipt), nl, iwk(jwipl), &
        xi(iip), yi(iip), iwk(jwit), iwk(jwiwk), wk )
    end do

  end if
!
!  Estimate partial derivatives at all data points,
!  for MD = 1, 2, 3.
!
  call idpdrv ( ndp0, xd, yd, zd, ncp0, iwk(jwipc), wk )
!
!  Interpolate the ZI values, for MD = 1,2,3;
!
  itpv = 0
  jwit = jwit0

  do iip = 1, nip0
    jwit = jwit + 1
    call idptip ( xd, yd, zd, nt, iwk(jwipt), nl, iwk(jwipl), wk, &
      iwk(jwit), xi(iip), yi(iip), zi(iip) )
  end do

  return
end
subroutine idcldp ( ndp, xd, yd, ncp, ipc )

!*****************************************************************************80
!
!! IDCLDP selects several neighbors of a given data point.
!
!  Discussion:
!
!    This subroutine selects several data points that are closest
!    to each of the data point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    Original FORTRAN77 version by Hiroshi Akima.
!    FORTRAN90 version by John Burkardt.
!
!  References:
!
!    Hiroshi Akima,
!    Algorithm 526:
!    A Method of Bivariate Interpolation and Smooth Surface Fitting
!    for Values Given at Irregularly Distributed Points,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, pages 160-164.
!
!    Hiroshi Akima,
!    On Estimating Partial Derivatives for Bivariate Interpolation
!    of Scattered Data,
!    Rocky Mountain Journal of Mathematics,
!    Volume 14, Number 1, Winter 1984, pages 41-51.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDP, the number of data points.
!    4 < =  NDP.
!
!    Input, real ( kind = 8 ) XD(NDP), Y(NDP), the X and Y coordinates
!    of the data points.
!
!    Input, integer ( kind = 4 ) NCP, the number of additional data points 
!    used for estimating partial derivatives at each data point.
!
!    Output, integer ( kind = 4 ) IPC(NCP*NDP), the indices of NCP data points 
!    closest to each of the NDP data points.
!
  implicit none

  integer ( kind = 4 ) ncp
  integer ( kind = 4 ) ndp

  real ( kind = 8 ) dsq0(ncp)
  real ( kind = 8 ) dsqf
  real ( kind = 8 ) dsqi
  real ( kind = 8 ) dsqmn
  real ( kind = 8 ) dsqmx
  real ( kind = 8 ) dx12
  real ( kind = 8 ) dx13
  real ( kind = 8 ) dy12
  real ( kind = 8 ) dy13
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip2
  integer ( kind = 4 ) ip2mn
  integer ( kind = 4 ) ip3
  integer ( kind = 4 ) ip3mn
  integer ( kind = 4 ) ipc(ncp*ndp)
  integer ( kind = 4 ) ipc0(ncp)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) j4
  integer ( kind = 4 ) jmx
  integer ( kind = 4 ) nclpt
  integer ( kind = 4 ) ncp0
  integer ( kind = 4 ) ndp0
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) x1
  real ( kind = 8 ) xd(ndp)
  real ( kind = 8 ) y1
  real ( kind = 8 ) yd(ndp)
!
!  Statement function.
!
  dsqf ( u1, v1, u2, v2 ) = ( u2 - u1 )**2 + ( v2 - v1 )**2
!
!  Preliminary processing.
!
  ndp0 = ndp
  ncp0 = ncp
  if(ndp0<2)  go to 90
  if(ncp0<1.or.ncp0>=ndp0)    go to 90
!
!  Calculation.
!
      do ip1 = 1, ndp0
!
!  Select NCP points.
!
        x1 = xd(ip1)
        y1 = yd(ip1)
        j1 = 0
        dsqmx = 0.0D+00

        do ip2 = 1,ndp0

          if(ip2==ip1) then
            cycle
          end if

          dsqi = dsqf(x1,y1,xd(ip2),yd(ip2))
          j1 = j1+1
          dsq0(j1) = dsqi
          ipc0(j1) = ip2

          if ( dsqmx < dsqi ) then
            dsqmx = dsqi
            jmx = j1
          end if

          if(j1>=ncp0) then
            exit
          end if

        end do

        ip2mn = ip2+1

        do ip2 = ip2mn,ndp0

          if(ip2==ip1) then
            cycle
          end if

          dsqi = dsqf(x1,y1,xd(ip2),yd(ip2))

          if(dsqi>=dsqmx) then
            cycle
          end if

          dsq0(jmx) = dsqi
          ipc0(jmx) = ip2
          dsqmx = 0.0D+00

          do j1 = 1,ncp0
            if ( dsqmx < dsq0(j1) ) then
              dsqmx = dsq0(j1)
              jmx = j1
            end if
          end do

        end do
!
!  Check if all the NCP+1 points are collinear.
!
        ip2 = ipc0(1)
        dx12 = xd(ip2)-x1
        dy12 = yd(ip2)-y1

        do j3 = 2,ncp0
          ip3 = ipc0(j3)
          dx13 = xd(ip3)-x1
          dy13 = yd(ip3)-y1
          if((dy13*dx12-dx13*dy12)/=0.0D+00 )    go to 50
        end do
!
!  Search for the closest noncollinear point.
!
   40   nclpt = 0

        do 43  ip3 = 1,ndp0

          if ( ip3 == ip1 ) then
            cycle
          end if

          do j4 = 1,ncp0
            if(ip3==ipc0(j4))     go to 43
          end do

          dx13 = xd(ip3)-x1
          dy13 = yd(ip3)-y1
          if((dy13*dx12-dx13*dy12)==0.0D+00 )    go to 43
          dsqi = dsqf(x1,y1,xd(ip3),yd(ip3))
          if(nclpt==0)       go to 42
          if(dsqi>=dsqmn)    go to 43
   42     nclpt = 1
          dsqmn = dsqi
          ip3mn = ip3

   43   continue

        if ( nclpt == 0 ) then
          go to 91
        end if

        dsqmx = dsqmn
        ipc0(jmx) = ip3mn
!
!  Replace the local array for the output array.
!
   50   continue

        j1 = (ip1-1)*ncp0

        do j2 = 1,ncp0
          j1 = j1+1
          ipc(j1) = ipc0(j2)
        end do

      end do

      return
!
!  error exit
!
   90 write ( *, 2090 )
      go to 92
   91 write ( *,2091)
   92 write ( *,2092)  ndp0,ncp0
      ipc(1) = 0
      return
 2090 format(1x/' ***   improper input parameter value(s).')
 2091 format(1x/' ***   all collinear data points.')
 2092 format('   ndp  = ',i5,5x,'ncp  = ',i5/ &
         ' error detected in routine   idcldp'/)
      end
      subroutine idgrid ( xd, yd, nt, ipt, nl, ipl, nxi, nyi, xi, yi, ngp, igp )

!*****************************************************************************80
!
!! IDGRID organizes grid points for surface fitting.
!
!  Discussion:
!
!    This subroutine organizes grid points for surface fitting by
!    sorting them in ascending order of triangle numbers and of the
!    border line segment number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2007
!
!  Author:
!
!    Original FORTRAN77 version by Hiroshi Akima.
!    FORTRAN90 version by John Burkardt.
!
!  References:
!
!    Hiroshi Akima,
!    Algorithm 526:
!    A Method of Bivariate Interpolation and Smooth Surface Fitting
!    for Values Given at Irregularly Distributed Points,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, pages 160-164.
!
!    Hiroshi Akima,
!    On Estimating Partial Derivatives for Bivariate Interpolation
!    of Scattered Data,
!    Rocky Mountain Journal of Mathematics,
!    Volume 14, Number 1, Winter 1984, pages 41-51.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XD(NDP), Y(NDP), the X and Y coordinates
!    of the data points.
!
!    Input, integer ( kind = 4 ) NT, the number of triangles.
!
!    Input, integer ( kind = 4 ) IPT(3*NT), the indices of the triangle vertexes.
!
!    Input, integer ( kind = 4 ) NL, the number of border line segments.
!
!    Input, integer ( kind = 4 ) IPL(3*NL), containing the point numbers of the
!    end points of the border line segments and their respective triangle
!    numbers.
!
!    Input, integer ( kind = 4 ) NXI, NYI, the number of grid points in the X and Y
!    coordinates.
!
!    Input, real ( kind = 8 ) XI(NXI), YI(NYI), the coordinates of the
!    grid points.
!
!    Output, integer ( kind = 4 ) NGP(2*(NT+2*NL)) where the
!    number of grid points that belong to each of the
!    triangles or of the border line segments are to be stored.
!
!    Output, integer ( kind = 4 ) IGP(NXI*NYI), where the grid point numbers
!    are to be stored in ascending order of the triangle number and the border
!    line segment number.
!
  implicit none

  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nxi
  integer ( kind = 4 ) nyi

  integer ( kind = 4 ) igp(nxi*nyi)
  integer ( kind = 4 ) il0
  integer ( kind = 4 ) il0t3
  integer ( kind = 4 ) ilp1
  integer ( kind = 4 ) ilp1t3
  integer ( kind = 4 ) insd
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip2
  integer ( kind = 4 ) ip3
  integer ( kind = 4 ) ipl(3*nl)
  integer ( kind = 4 ) ipt(3*nt)
  integer ( kind = 4 ) it0
  integer ( kind = 4 ) it0t3
  integer ( kind = 4 ) ixi
  integer ( kind = 4 ) iximx
  integer ( kind = 4 ) iximn
  integer ( kind = 4 ) iyi
  integer ( kind = 4 ) izi
  integer ( kind = 4 ) jigp0
  integer ( kind = 4 ) jigp1
  integer ( kind = 4 ) jigp1i
  integer ( kind = 4 ) jngp0
  integer ( kind = 4 ) jngp1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ngp(2*nt+4*nl)
  integer ( kind = 4 ) ngp0
  integer ( kind = 4 ) ngp1
  integer ( kind = 4 ) nl0
  integer ( kind = 4 ) nt0
  integer ( kind = 4 ) nxi0
  integer ( kind = 4 ) nxinyi
  integer ( kind = 4 ) nyi0
  real ( kind = 8 ) side
  real ( kind = 8 ) spdt
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) u3
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) v3
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xd(*)
  real ( kind = 8 ) xi(nxi)
  real ( kind = 8 ) xii
  real ( kind = 8 ) ximn
  real ( kind = 8 ) ximx
  real ( kind = 8 ) xmn
  real ( kind = 8 ) xmx
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) yd(*)
  real ( kind = 8 ) yi(nyi)
  real ( kind = 8 ) yii
  real ( kind = 8 ) yimn
  real ( kind = 8 ) yimx
  real ( kind = 8 ) ymn
  real ( kind = 8 ) ymx
!
!  Statement functions.
!
      side(u1,v1,u2,v2,u3,v3) = (u1-u3)*(v2-v3) - (v1-v3)*(u2-u3)
      spdt(u1,v1,u2,v2,u3,v3) = (u1-u2)*(u3-u2) + (v1-v2)*(v3-v2)
!
!  Preliminary processing.
!
      nt0 = nt
      nl0 = nl
      nxi0 = nxi
      nyi0 = nyi
      nxinyi = nxi0*nyi0
      ximn = min (xi(1),xi(nxi0))
      ximx = max (xi(1),xi(nxi0))
      yimn = min (yi(1),yi(nyi0))
      yimx = max (yi(1),yi(nyi0))
!
!  Determine grid points inside the data area.
!
      jngp0 = 0
      jngp1 = 2*(nt0+2*nl0) + 1
      jigp0 = 0
      jigp1 = nxinyi + 1

      do 160 it0 = 1,nt0

        ngp0 = 0
        ngp1 = 0
        it0t3 = it0*3
        ip1 = ipt(it0t3-2)
        ip2 = ipt(it0t3-1)
        ip3 = ipt(it0t3)
        x1 = xd(ip1)
        y1 = yd(ip1)
        x2 = xd(ip2)
        y2 = yd(ip2)
        x3 = xd(ip3)
        y3 = yd(ip3)
        xmn = min (x1,x2,x3)
        xmx = max (x1,x2,x3)
        ymn = min (y1,y2,y3)
        ymx = max (y1,y2,y3)
        insd = 0

        do ixi = 1,nxi0
          if (xi(ixi)>=xmn .and. xi(ixi)<=xmx) go to 10
          if (insd==0) go to 20
          iximx = ixi - 1
          go to 30

   10     continue

          if (insd==1) go to 20
          insd = 1
          iximn = ixi
   20     continue
        end do

        if (insd==0) go to 150
        iximx = nxi0

   30   do 140 iyi = 1,nyi0

          yii = yi(iyi)
          if (yii<ymn .or. yii>ymx) go to 140

          do 130 ixi = iximn,iximx

            xii = xi(ixi)
            l = 0
            if (side(x1,y1,x2,y2,xii,yii)) 130, 40, 50
   40       l = 1
   50       if (side(x2,y2,x3,y3,xii,yii)) 130, 60, 70
   60       l = 1
   70       if (side(x3,y3,x1,y1,xii,yii)) 130, 80, 90
   80       l = 1
   90       izi = nxi0*(iyi-1) + ixi
            if (l==1) go to 100
            ngp0 = ngp0 + 1
            jigp0 = jigp0 + 1
            igp(jigp0) = izi
            go to 130

  100       continue

            do jigp1i = jigp1,nxinyi
              if (izi==igp(jigp1i)) go to 130
            end do

  120       ngp1 = ngp1 + 1
            jigp1 = jigp1 - 1
            igp(jigp1) = izi
  130     continue
  140   continue
  150   jngp0 = jngp0 + 1
        ngp(jngp0) = ngp0
        jngp1 = jngp1 - 1
        ngp(jngp1) = ngp1

  160 continue
!
!  Determine grid points outside the data area.
!  in semi-infinite rectangular area.
!
      do 450 il0 = 1,nl0

        ngp0 = 0
        ngp1 = 0
        il0t3 = il0*3
        ip1 = ipl(il0t3-2)
        ip2 = ipl(il0t3-1)
        x1 = xd(ip1)
        y1 = yd(ip1)
        x2 = xd(ip2)
        y2 = yd(ip2)
        xmn = ximn
        xmx = ximx
        ymn = yimn
        ymx = yimx

        if (y2>=y1) xmn = min (x1,x2)
        if (y2<=y1) xmx = max (x1,x2)
        if (x2<=x1) ymn = min (y1,y2)
        if (x2>=x1) ymx = max (y1,y2)
        insd = 0

        do 180 ixi = 1,nxi0

          if (xi(ixi)>=xmn .and. xi(ixi)<=xmx) go to 170
          if (insd==0) go to 180
          iximx = ixi - 1
          go to 190
  170     continue

          if ( insd /= 1 ) then
            insd = 1
            iximn = ixi
          end if

  180   continue

        if (insd==0) go to 310
        iximx = nxi0

  190   do 300 iyi = 1,nyi0

          yii = yi(iyi)
          if (yii<ymn .or. yii>ymx) go to 300

          do 290 ixi = iximn,iximx

            xii = xi(ixi)
            l = 0
            if (side(x1,y1,x2,y2,xii,yii)) 210, 200, 290
  200       l = 1
  210       if (spdt(x2,y2,x1,y1,xii,yii)) 290, 220, 230
  220       l = 1
  230       if (spdt(x1,y1,x2,y2,xii,yii)) 290, 240, 250
  240       l = 1
  250       izi = nxi0*(iyi-1) + ixi
            if (l==1) go to 260
            ngp0 = ngp0 + 1
            jigp0 = jigp0 + 1
            igp(jigp0) = izi
            go to 290
  260       continue

            do jigp1i = jigp1,nxinyi
              if (izi==igp(jigp1i)) go to 290
            end do

            ngp1 = ngp1 + 1
            jigp1 = jigp1 - 1
            igp(jigp1) = izi

  290     continue

  300   continue

  310   jngp0 = jngp0 + 1
        ngp(jngp0) = ngp0
        jngp1 = jngp1 - 1
        ngp(jngp1) = ngp1
!
!  In semi-infinite triangular area.
!
        ngp0 = 0
        ngp1 = 0
        ilp1 = mod(il0,nl0) + 1
        ilp1t3 = ilp1*3
        ip3 = ipl(ilp1t3-1)
        x3 = xd(ip3)
        y3 = yd(ip3)
        xmn = ximn
        xmx = ximx
        ymn = yimn
        ymx = yimx
        if (y3>=y2 .and. y2>=y1) xmn = x2
        if (y3<=y2 .and. y2<=y1) xmx = x2
        if (x3<=x2 .and. x2<=x1) ymn = y2
        if (x3>=x2 .and. x2>=x1) ymx = y2
        insd = 0

        do 330 ixi = 1,nxi0

          if (xi(ixi)>=xmn .and. xi(ixi)<=xmx) go to 320
          if (insd==0) go to 330
          iximx = ixi - 1
          go to 340
  320     continue

          if ( insd /= 1 ) then
            insd = 1
            iximn = ixi
          end if

  330   continue

        if (insd==0) go to 440

        iximx = nxi0
  340   do 430 iyi = 1,nyi0
          yii = yi(iyi)
          if (yii<ymn .or. yii>ymx) go to 430
          do 420 ixi = iximn,iximx
            xii = xi(ixi)
            l = 0
            if (spdt(x1,y1,x2,y2,xii,yii)) 360, 350, 420
  350       l = 1
  360       if (spdt(x3,y3,x2,y2,xii,yii)) 380, 370, 420
  370       l = 1
  380       izi = nxi0*(iyi-1) + ixi
            if (l==1) go to 390
            ngp0 = ngp0 + 1
            jigp0 = jigp0 + 1
            igp(jigp0) = izi
            go to 420
  390       continue
            do jigp1i = jigp1,nxinyi
              if (izi==igp(jigp1i)) go to 420
            end do
            ngp1 = ngp1 + 1
            jigp1 = jigp1 - 1
            igp(jigp1) = izi
  420     continue

  430   continue

  440   jngp0 = jngp0 + 1
        ngp(jngp0) = ngp0
        jngp1 = jngp1 - 1
        ngp(jngp1) = ngp1

  450 continue

      return
      end
      subroutine idlctn ( ndp, xd, yd, nt, ipt, nl, ipl, xii, yii, iti, &
        iwk, wk)

!*****************************************************************************80
!
!! IDLCTN finds the triangle that contains a point.
!
!  Discussion:
!
!    This subroutine locates a point, that is, determines to what triangle
!    a given point (XII,YII) belongs.  When the given point
!    does not lie inside the data area, this subroutine determines
!    the border line segment when the point lies in an outside
!    rectangular area, and two border line segments when the point
!    lies in an outside triangular area.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2007
!
!  Author:
!
!    Original FORTRAN77 version by Hiroshi Akima.
!    FORTRAN90 version by John Burkardt.
!
!  References:
!
!    Hiroshi Akima,
!    Algorithm 526:
!    A Method of Bivariate Interpolation and Smooth Surface Fitting
!    for Values Given at Irregularly Distributed Points,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, pages 160-164.
!
!    Hiroshi Akima,
!    On Estimating Partial Derivatives for Bivariate Interpolation
!    of Scattered Data,
!    Rocky Mountain Journal of Mathematics,
!    Volume 14, Number 1, Winter 1984, pages 41-51.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDP, the number of data points.
!    4 < =  NDP.
!
!    Input, real ( kind = 8 ) XD(NDP), Y(NDP), the X and Y coordinates
!    of the data points.
!
!    Input, integer ( kind = 4 ) NT, the number of triangles.
!
!    Input, integer ( kind = 4 ) IPT(3*NT), the point numbers of the vertexes of
!    the triangles,
!
!    Input, integer ( kind = 4 ) NL, the number of border line segments.
!
!    Input, integer ( kind = 4 ) IPL(3*NL), the point numbers of the end points of
!    the border line segments and their respective triangle numbers.
!
!    Input, real ( kind = 8 ) XII, YII, the coordinates of the point
!    to be located.
!
!    Output, integer ( kind = 4 ) ITI, the triangle number, when the point is
!    inside the data area, or two border line segment numbers, il1 and il2,
!    coded to il1*(nt+nl)+il2, when the point is outside the data area.
!
!    Input/output, integer ( kind = 4 ) IWK(18*NDP), real ( kind = 8 ) WK(8*NDP),
!    arrays that contain information about the triangulation of the data.
!    These arrays are set up on the first call to IDLCTN for a given set of
!    data; the data should not be altered by the user between calls for the
!    same set of data, as it will be reused rather than being recomputed.
!
  implicit none

  integer ( kind = 4 ) ndp
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nl

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) idp
  integer ( kind = 4 ) idsc(9)
  integer ( kind = 4 ) il1
  integer ( kind = 4 ) il1t3
  integer ( kind = 4 ) il2
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip2
  integer ( kind = 4 ) ip3
  integer ( kind = 4 ) ipl(3*nl)
  integer ( kind = 4 ) ipt(3*nt)
  integer ( kind = 4 ) isc
  integer ( kind = 4 ) it0
  integer ( kind = 4 ) it0t3
  integer ( kind = 4 ) iti
  integer ( kind = 4 ) itipv
  integer ( kind = 4 ) itsc
  integer ( kind = 4 ) iwk(18*ndp)
  integer ( kind = 4 ) jiwk
  integer ( kind = 4 ) jwk
  integer ( kind = 4 ) ndp0
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nl0
  integer ( kind = 4 ) nt0
  integer ( kind = 4 ) ntl
  integer ( kind = 4 ) ntsc(9)
  integer ( kind = 4 ) ntsci
  real ( kind = 8 ) side
  real ( kind = 8 ) spdt
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) u3
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) v3
  real ( kind = 8 ) wk(8*ndp)
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xd(ndp)
  real ( kind = 8 ) xi
  real ( kind = 8 ) xii
  real ( kind = 8 ) xmn
  real ( kind = 8 ) xmx
  real ( kind = 8 ) xs1
  real ( kind = 8 ) xs2
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) yd(ndp)
  real ( kind = 8 ) yi
  real ( kind = 8 ) yii
  real ( kind = 8 ) ymn
  real ( kind = 8 ) ymx
  real ( kind = 8 ) ys1
  real ( kind = 8 ) ys2

  save /idlc/

  common /idlc/ nit
!
!  Statement functions.
!
      side(u1,v1,u2,v2,u3,v3) = (u1-u3)*(v2-v3) - (v1-v3)*(u2-u3)

      spdt(u1,v1,u2,v2,u3,v3) = (u1-u2)*(u3-u2) + (v1-v2)*(v3-v2)
!
!  Preliminary processing.
!
      ndp0 = ndp
      nt0 = nt
      nl0 = nl
      ntl = nt0 + nl0
      x0 = xii
      y0 = yii
!
!  Processing for a new set of data points.
!
      if ( nit/=0) go to 80

      nit = 1
!
!  Divide the XY plane into nine rectangular sections.
!
      xmn = xd(1)
      xmx = xmn
      ymn = yd(1)
      ymx = ymn

      do idp = 2,ndp0
        xi = xd(idp)
        yi = yd(idp)
        xmn = min (xi,xmn)
        xmx = max (xi,xmx)
        ymn = min (yi,ymn)
        ymx = max (yi,ymx)
      end do

      xs1 = ( xmn + xmn + xmx ) / 3.0D+00
      xs2 = ( xmn + xmx + xmx ) / 3.0D+00
      ys1 = ( ymn + ymn + ymx ) / 3.0D+00
      ys2 = ( ymn + ymx + ymx ) / 3.0D+00
!
!  Determine and store in the IWK array triangle numbers of
!  the triangles associated with each of the nine sections.
!
      ntsc(1:9) = 0
      idsc(1:9) = 0

      it0t3 = 0
      jwk = 0

      do 70 it0 = 1,nt0

        it0t3 = it0t3 + 3
        i1 = ipt(it0t3-2)
        i2 = ipt(it0t3-1)
        i3 = ipt(it0t3)
        xmn = min (xd(i1),xd(i2),xd(i3))
        xmx = max (xd(i1),xd(i2),xd(i3))
        ymn = min (yd(i1),yd(i2),yd(i3))
        ymx = max (yd(i1),yd(i2),yd(i3))
        if (ymn>ys1) go to 30
        if (xmn<=xs1) idsc(1) = 1
        if (xmx>=xs1 .and. xmn<=xs2) idsc(2) = 1
        if (xmx>=xs2) idsc(3) = 1
   30   if (ymx<ys1 .or. ymn>ys2) go to 40
        if (xmn<=xs1) idsc(4) = 1
        if (xmx>=xs1 .and. xmn<=xs2) idsc(5) = 1
        if (xmx>=xs2) idsc(6) = 1
   40   if (ymx<ys2) go to 50
        if (xmn<=xs1) idsc(7) = 1
        if (xmx>=xs1 .and. xmn<=xs2) idsc(8) = 1
        if (xmx>=xs2) idsc(9) = 1
   50   continue

        do isc = 1,9
          if ( idsc(isc) /= 0 ) then
            jiwk = 9*ntsc(isc) + isc
            iwk(jiwk) = it0
            ntsc(isc) = ntsc(isc) + 1
            idsc(isc) = 0
          end if
        end do
!
!  Store in the WK array the minimum and maximum of the X and
!  Y coordinate values for each of the triangle.
!
        jwk = jwk + 4
        wk(jwk-3) = xmn
        wk(jwk-2) = xmx
        wk(jwk-1) = ymn
        wk(jwk) = ymx

   70 continue

      go to 110
!
!  Check if in the same triangle as previous.
!
   80 it0 = itipv
      if (it0>nt0) go to 90
      it0t3 = it0*3
      ip1 = ipt(it0t3-2)
      x1 = xd(ip1)
      y1 = yd(ip1)
      ip2 = ipt(it0t3-1)
      x2 = xd(ip2)
      y2 = yd(ip2)
      if (side(x1,y1,x2,y2,x0,y0)<0.0D+00 ) go to 110
      ip3 = ipt(it0t3)
      x3 = xd(ip3)
      y3 = yd(ip3)
      if (side(x2,y2,x3,y3,x0,y0)<0.0D+00 ) go to 110
      if (side(x3,y3,x1,y1,x0,y0)<0.0D+00 ) go to 110
      go to 170
!
!  Check if on the same border line segment.
!
   90 il1 = it0/ntl
      il2 = it0 - il1*ntl
      il1t3 = il1*3
      ip1 = ipl(il1t3-2)
      x1 = xd(ip1)
      y1 = yd(ip1)
      ip2 = ipl(il1t3-1)
      x2 = xd(ip2)
      y2 = yd(ip2)
      if (il2/=il1) go to 100
      if (spdt(x1,y1,x2,y2,x0,y0)<0.0D+00 ) go to 110
      if (spdt(x2,y2,x1,y1,x0,y0)<0.0D+00 ) go to 110
      if (side(x1,y1,x2,y2,x0,y0)>0.0D+00 ) go to 110
      go to 170
!
!  Check if between the same two border line segments.
!
  100 if (spdt(x1,y1,x2,y2,x0,y0)>0.0D+00 ) go to 110
      ip3 = ipl(3*il2-1)
      x3 = xd(ip3)
      y3 = yd(ip3)
      if (spdt(x3,y3,x2,y2,x0,y0)<=0.0D+00 ) go to 170
!
!  Locate inside the data area.
!  Determine the section in which the point in question lies.
!
  110 isc = 1
      if (x0>=xs1) isc = isc + 1
      if (x0>=xs2) isc = isc + 1
      if (y0>=ys1) isc = isc + 3
      if (y0>=ys2) isc = isc + 3
!
!  Search through the triangles associated with the section.
!
      ntsci = ntsc(isc)
      if (ntsci<=0) go to 130
      jiwk = -9 + isc

      do 120 itsc = 1,ntsci
        jiwk = jiwk + 9
        it0 = iwk(jiwk)
        jwk = it0*4
        if (x0<wk(jwk-3)) go to 120
        if (x0>wk(jwk-2)) go to 120
        if (y0<wk(jwk-1)) go to 120
        if (y0>wk(jwk)) go to 120
        it0t3 = it0*3
        ip1 = ipt(it0t3-2)
        x1 = xd(ip1)
        y1 = yd(ip1)
        ip2 = ipt(it0t3-1)
        x2 = xd(ip2)
        y2 = yd(ip2)
        if (side(x1,y1,x2,y2,x0,y0)<0.0D+00 ) go to 120
        ip3 = ipt(it0t3)
        x3 = xd(ip3)
        y3 = yd(ip3)
        if (side(x2,y2,x3,y3,x0,y0)<0.0D+00 ) go to 120
        if (side(x3,y3,x1,y1,x0,y0)<0.0D+00 ) go to 120
        go to 170
  120 continue
!
!  Locate outside the data area.
!
  130 do 150 il1 = 1,nl0
        il1t3 = il1*3
        ip1 = ipl(il1t3-2)
        x1 = xd(ip1)
        y1 = yd(ip1)
        ip2 = ipl(il1t3-1)
        x2 = xd(ip2)
        y2 = yd(ip2)
        if (spdt(x2,y2,x1,y1,x0,y0)<0.0D+00 ) go to 150
        if (spdt(x1,y1,x2,y2,x0,y0)<0.0D+00 ) go to 140
        if (side(x1,y1,x2,y2,x0,y0)>0.0D+00 ) go to 150
        il2 = il1
        go to 160
  140   il2 = mod(il1,nl0) + 1
        ip3 = ipl(3*il2-1)
        x3 = xd(ip3)
        y3 = yd(ip3)
        if (spdt(x3,y3,x2,y2,x0,y0)<=0.0D+00) go to 160
  150 continue

      it0 = 1
      go to 170
  160 it0 = il1*ntl + il2
!
!  normal exit
!
  170 iti = it0
      itipv = it0

      return
      end
      subroutine idpdrv ( ndp, xd, yd, zd, ncp, ipc, pd )

!*****************************************************************************80
!
!! IDPDRV estimates first and second partial derivatives at data points.
!
!  Discussion:
!
!    This subroutine estimates partial derivatives of the first and
!    second order at the data points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2007
!
!  Author:
!
!    Original FORTRAN77 version by Hiroshi Akima.
!    FORTRAN90 version by John Burkardt.
!
!  References:
!
!    Hiroshi Akima,
!    Algorithm 526:
!    A Method of Bivariate Interpolation and Smooth Surface Fitting
!    for Values Given at Irregularly Distributed Points,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, pages 160-164.
!
!    Hiroshi Akima,
!    On Estimating Partial Derivatives for Bivariate Interpolation
!    of Scattered Data,
!    Rocky Mountain Journal of Mathematics,
!    Volume 14, Number 1, Winter 1984, pages 41-51.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDP, the number of data points.
!    4 < =  NDP.
!
!    Input, real ( kind = 8 ) XD(NDP), Y(NDP), Z(NDP), the X, Y and
!    Z coordinates of the data points.
!
!    Input, integer ( kind = 4 ) NCP, the number of additional data points used for
!    estimating partial derivatives at each data point,
!
!    Input, integer ( kind = 4 ) IPC(NCP*NDP), the indices of NCP data points
!    closest to each of the NDP data points.
!
!    Output, real ( kind = 8 ) PD(5*NDP), the estimated ZX, ZY, ZXX, ZXY,
!    and ZYY values at the Ith data point are to be stored as the
!    (5*I-4)th, (5*I-3)rd, (5*I-2)nd, (5*I-1)st and (5*I)th elements,
!    respectively, where I = 1, 2, ..., NDP.
!
  implicit none

  integer ( kind = 4 ) ncp
  integer ( kind = 4 ) ndp

  real ( kind = 8 ) dnmx
  real ( kind = 8 ) dnmxx
  real ( kind = 8 ) dnmxy
  real ( kind = 8 ) dnmy
  real ( kind = 8 ) dnmyx
  real ( kind = 8 ) dnmyy
  real ( kind = 8 ) dnmz
  real ( kind = 8 ) dx1
  real ( kind = 8 ) dx2
  real ( kind = 8 ) dy1
  real ( kind = 8 ) dy2
  real ( kind = 8 ) dz1
  real ( kind = 8 ) dz2
  real ( kind = 8 ) dzx1
  real ( kind = 8 ) dzx2
  real ( kind = 8 ) dzy1
  real ( kind = 8 ) dzy2
  integer ( kind = 4 ) ic1
  integer ( kind = 4 ) ic2
  integer ( kind = 4 ) ic2mn
  integer ( kind = 4 ) idp
  integer ( kind = 4 ) ip0
  integer ( kind = 4 ) ipc(ncp*ndp)
  integer ( kind = 4 ) ipi
  integer ( kind = 4 ) jipc
  integer ( kind = 4 ) jipc0
  integer ( kind = 4 ) jpd
  integer ( kind = 4 ) jpd0
  integer ( kind = 4 ) ncp0
  integer ( kind = 4 ) ncpm1
  integer ( kind = 4 ) ndp0
  real ( kind = 8 ) nmx
  real ( kind = 8 ) nmxx
  real ( kind = 8 ) nmxy
  real ( kind = 8 ) nmy
  real ( kind = 8 ) nmyy
  real ( kind = 8 ) nmyx
  real ( kind = 8 ) nmz
  real ( kind = 8 ) pd(5*ndp)
  real ( kind = 8 ) x0
  real ( kind = 8 ) xd(ndp)
  real ( kind = 8 ) y0
  real ( kind = 8 ) yd(ndp)
  real ( kind = 8 ) z0
  real ( kind = 8 ) zd(ndp)
  real ( kind = 8 ) zx0
  real ( kind = 8 ) zy0
!
!  Preliminary processing.
!
  ndp0 = ndp
  ncp0 = ncp
  ncpm1 = ncp0-1
!
!  Estimation of ZX and ZY.
!
      do 24  ip0 = 1,ndp0

        x0 = xd(ip0)
        y0 = yd(ip0)
        z0 = zd(ip0)
        nmx = 0.0D+00
        nmy = 0.0D+00
        nmz = 0.0D+00
        jipc0 = ncp0*(ip0-1)

        do 23  ic1 = 1,ncpm1

          jipc = jipc0+ic1
          ipi = ipc(jipc)
          dx1 = xd(ipi)-x0
          dy1 = yd(ipi)-y0
          dz1 = zd(ipi)-z0
          ic2mn = ic1+1

          do ic2 = ic2mn,ncp0

            jipc = jipc0+ic2
            ipi = ipc(jipc)
            dx2 = xd(ipi)-x0
            dy2 = yd(ipi)-y0
            dnmz = dx1*dy2-dy1*dx2

            if ( dnmz == 0.0D+00 ) then
              cycle
            end if

            dz2 = zd(ipi)-z0
            dnmx = dy1*dz2-dz1*dy2
            dnmy = dz1*dx2-dx1*dz2

            if ( dnmz < 0.0D+00 ) then
              dnmx = -dnmx
              dnmy = -dnmy
              dnmz = -dnmz
            end if

            nmx = nmx+dnmx
            nmy = nmy+dnmy
            nmz = nmz+dnmz

          end do

   23   continue

        jpd0 = 5*ip0
        pd(jpd0-4) = -nmx/nmz
        pd(jpd0-3) = -nmy/nmz

   24 continue
!
!  Estimation of ZXX, ZXY, and ZYY.
!
   30 continue

      do ip0 = 1,ndp0

        jpd0 = jpd0+5
        x0 = xd(ip0)
        jpd0 = 5*ip0
        y0 = yd(ip0)
        zx0 = pd(jpd0-4)
        zy0 = pd(jpd0-3)

        nmxx = 0.0D+00
        nmxy = 0.0D+00
        nmyx = 0.0D+00
        nmyy = 0.0D+00
        nmz  = 0.0D+00

        jipc0 = ncp0*(ip0-1)

        do ic1 = 1,ncpm1

          jipc = jipc0+ic1
          ipi = ipc(jipc)
          dx1 = xd(ipi)-x0
          dy1 = yd(ipi)-y0
          jpd = 5*ipi
          dzx1 = pd(jpd-4)-zx0
          dzy1 = pd(jpd-3)-zy0
          ic2mn = ic1+1

          do ic2 = ic2mn,ncp0

            jipc = jipc0+ic2
            ipi = ipc(jipc)
            dx2 = xd(ipi)-x0
            dy2 = yd(ipi)-y0
            dnmz  = dx1*dy2 -dy1*dx2

            if ( dnmz == 0.0D+00 ) then
              cycle
            end if

            jpd = 5 * ipi
            dzx2 = pd(jpd-4) - zx0
            dzy2 = pd(jpd-3) - zy0
            dnmxx = dy1 * dzx2 - dzx1 * dy2
            dnmxy = dzx1 * dx2 - dx1 * dzx2
            dnmyx = dy1 * dzy2 - dzy1 * dy2
            dnmyy = dzy1 * dx2 - dx1 * dzy2

            if ( dnmz < 0.0D+00 ) then
              dnmxx = -dnmxx
              dnmxy = -dnmxy
              dnmyx = -dnmyx
              dnmyy = -dnmyy
              dnmz  = -dnmz
            end if

            nmxx = nmxx+dnmxx
            nmxy = nmxy+dnmxy
            nmyx = nmyx+dnmyx
            nmyy = nmyy+dnmyy
            nmz  = nmz +dnmz

          end do

        end do

        pd(jpd0-2) = - nmxx / nmz
        pd(jpd0-1) = - ( nmxy + nmyx ) / ( 2.0D+00 * nmz )
        pd(jpd0)   = - nmyy / nmz

      end do

      return
      end
      subroutine idptip ( xd, yd, zd, nt, ipt, nl, ipl, pdd, iti, xii, yii, &
        zii )

!*****************************************************************************80
!
!! IDPTIP performs interpolation, determining a value of Z given X and Y.
!
!  Discussion:
!
!    This subroutine performs punctual interpolation or extrapolation
!    that is, determines the Z value at a point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2007
!
!  Author:
!
!    Original FORTRAN77 version by Hiroshi Akima.
!    FORTRAN90 version by John Burkardt.
!
!  References:
!
!    Hiroshi Akima,
!    Algorithm 526:
!    A Method of Bivariate Interpolation and Smooth Surface Fitting
!    for Values Given at Irregularly Distributed Points,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, pages 160-164.
!
!    Hiroshi Akima,
!    On Estimating Partial Derivatives for Bivariate Interpolation
!    of Scattered Data,
!    Rocky Mountain Journal of Mathematics,
!    Volume 14, Number 1, Winter 1984, pages 41-51.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XD(NDP), Y(NDP), Z(NDP), the X, Y and Z
!    coordinates of the data points.
!
!    Input, integer ( kind = 4 ) NT, the number of triangles.
!
!    Input, integer ( kind = 4 ) IPT(3*NT), the point numbers of the vertexes of
!    the triangles.
!
!    Input, integer ( kind = 4 ) NL, the number of border line segments.
!
!    Input, integer ( kind = 4 ) IPL(3*NL), the point numbers of the end points
!    of the border line segments and their respective triangle numbers.
!
!    Input, real ( kind = 8 ) PDD(5*NDP). the partial derivatives at
!    the data points,
!
!    Input, integer ( kind = 4 ) ITI, triangle number of the triangle in which lies
!    the point for which interpolation is to be performed,
!
!    Input, real ( kind = 8 ) XII, YII, the X and Y coordinates of the
!    point for which interpolation is to be performed.
!
!    Output, real ( kind = 8 ) ZII, the interpolated Z value.
!
  implicit none

  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nt

  real ( kind = 8 ) a
  real ( kind = 8 ) aa
  real ( kind = 8 ) ab
  real ( kind = 8 ) ac
  real ( kind = 8 ) act2
  real ( kind = 8 ) ad
  real ( kind = 8 ) adbc
  real ( kind = 8 ) ap
  real ( kind = 8 ) b
  real ( kind = 8 ) bb
  real ( kind = 8 ) bc
  real ( kind = 8 ) bdt2
  real ( kind = 8 ) bp
  real ( kind = 8 ) c
  real ( kind = 8 ) cc
  real ( kind = 8 ) cd
  real ( kind = 8 ) cp
  real ( kind = 8 ) csuv
  real ( kind = 8 ) d
  real ( kind = 8 ) dd
  real ( kind = 8 ) dlt
  real ( kind = 8 ) dp
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  real ( kind = 8 ) h3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idp
  integer ( kind = 4 ) il1
  integer ( kind = 4 ) il2
  integer ( kind = 4 ) ipl(3*nl)
  integer ( kind = 4 ) ipt(3*nt)
  integer ( kind = 4 ) it0
  integer ( kind = 4 ) iti
  integer ( kind = 4 ) itpv
  integer ( kind = 4 ) jipl
  integer ( kind = 4 ) jipt
  integer ( kind = 4 ) jpd
  integer ( kind = 4 ) jpdd
  integer ( kind = 4 ) kpd
  real ( kind = 8 ) lu
  real ( kind = 8 ) lv
  integer ( kind = 4 ) ntl
  real ( kind = 8 ) p0
  real ( kind = 8 ) p00
  real ( kind = 8 ) p01
  real ( kind = 8 ) p02
  real ( kind = 8 ) p03
  real ( kind = 8 ) p04
  real ( kind = 8 ) p05
  real ( kind = 8 ) p1
  real ( kind = 8 ) p10
  real ( kind = 8 ) p11
  real ( kind = 8 ) p12
  real ( kind = 8 ) p13
  real ( kind = 8 ) p14
  real ( kind = 8 ) p2
  real ( kind = 8 ) p20
  real ( kind = 8 ) p21
  real ( kind = 8 ) p22
  real ( kind = 8 ) p23
  real ( kind = 8 ) p3
  real ( kind = 8 ) p30
  real ( kind = 8 ) p31
  real ( kind = 8 ) p32
  real ( kind = 8 ) p4
  real ( kind = 8 ) p40
  real ( kind = 8 ) p41
  real ( kind = 8 ) p5
  real ( kind = 8 ) p50
  real ( kind = 8 ) pd(15)
  real ( kind = 8 ) pdd(*)
  real ( kind = 8 ) thsv
  real ( kind = 8 ) thus
  real ( kind = 8 ) thuv
  real ( kind = 8 ) thxu
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x(3)
  real ( kind = 8 ) x0
  real ( kind = 8 ) xd(*)
  real ( kind = 8 ) xii
  real ( kind = 8 ) y(3)
  real ( kind = 8 ) y0
  real ( kind = 8 ) yd(*)
  real ( kind = 8 ) yii
  real ( kind = 8 ) z(3)
  real ( kind = 8 ) zd(*)
  real ( kind = 8 ) zii
  real ( kind = 8 ) zu(3)
  real ( kind = 8 ) zuu(3)
  real ( kind = 8 ) zuv(3)
  real ( kind = 8 ) zv(3)
  real ( kind = 8 ) zvv(3)

  save /idpi/

  common/idpi/itpv

  equivalence (p5,p50)
!
!  Preliminary processing.
!
      it0 = iti
      ntl = nt+nl
      if(it0<=ntl)      go to 20
      il1 = it0/ntl
      il2 = it0-il1*ntl
      if(il1==il2)      go to 40
      go to 60
!
!  Calculation of ZII by interpolation.
!  Check if the necessary coefficients have been calculated.
!
   20 if(it0==itpv)     go to 30
!
!  Load coordinate and partial derivative values at the vertices.
!
   21 continue

      jipt = 3*(it0-1)
      jpd = 0

      do i = 1,3

        jipt = jipt+1
        idp = ipt(jipt)
        x(i) = xd(idp)
        y(i) = yd(idp)
        z(i) = zd(idp)
        jpdd = 5*(idp-1)

        do kpd = 1,5
          jpd = jpd+1
          jpdd = jpdd+1
          pd(jpd) = pdd(jpdd)
        end do

      end do
!
!  Determine the coefficients for the coordinate system
!  transformation from the XY system to the UV system
!  and vice versa.
!
   24 continue

      x0 = x(1)
      y0 = y(1)
      a = x(2)-x0
      b = x(3)-x0
      c = y(2)-y0
      d = y(3)-y0
      ad = a*d
      bc = b*c
      dlt = ad-bc
      ap =  d/dlt
      bp = -b/dlt
      cp = -c/dlt
      dp =  a/dlt
!
!  Convert the partial derivatives at the vertexes of the
!  triangle for the UV coordinate system.
!
   25 continue

      aa = a*a
      act2 = 2.0D+00 *a*c
      cc = c*c
      ab = a*b
      adbc = ad+bc
      cd = c*d
      bb = b*b
      bdt2 = 2.0D+00 *b*d
      dd = d*d

      do i = 1,3
        jpd = 5*i
        zu(i) = a*pd(jpd-4)+c*pd(jpd-3)
        zv(i) = b*pd(jpd-4)+d*pd(jpd-3)
        zuu(i) = aa*pd(jpd-2)+act2*pd(jpd-1)+cc*pd(jpd)
        zuv(i) = ab*pd(jpd-2)+adbc*pd(jpd-1)+cd*pd(jpd)
        zvv(i) = bb*pd(jpd-2)+bdt2*pd(jpd-1)+dd*pd(jpd)
      end do
!
!  Calculate the coefficients of the polynomial.
!
   27 continue

      p00 = z(1)
      p10 = zu(1)
      p01 = zv(1)
      p20 = 0.5D+00*zuu(1)
      p11 = zuv(1)
      p02 = 0.5D+00*zvv(1)
      h1 = z(2)-p00-p10-p20
      h2 = zu(2)-p10-zuu(1)
      h3 = zuu(2)-zuu(1)
      p30 =  10.0D+00*h1-4.0D+00*h2+0.5D+00*h3
      p40 = -15.0D+00*h1+7.0D+00*h2        -h3
      p50 =   6.0D+00*h1-3.0D+00*h2+0.5D+00*h3
      h1 = z(3)-p00-p01-p02
      h2 = zv(3)-p01-zvv(1)
      h3 = zvv(3)-zvv(1)
      p03 =  10.0D+00*h1-4.0D+00*h2+0.5D+00*h3
      p04 = -15.0D+00*h1+7.0D+00*h2        -h3
      p05 =   6.0D+00*h1-3.0D+00*h2+0.5D+00*h3
      lu = sqrt(aa+cc)
      lv = sqrt(bb+dd)
      thxu = atan2(c,a)
      thuv = atan2(d,b)-thxu
      csuv = cos(thuv)
      p41 = 5.0D+00*lv*csuv/lu*p50
      p14 = 5.0D+00*lu*csuv/lv*p05
      h1 = zv(2)-p01-p11-p41
      h2 = zuv(2)-p11-4.0D+00*p41
      p21 =  3.0D+00*h1-h2
      p31 = -2.0D+00*h1+h2
      h1 = zu(3)-p10-p11-p14
      h2 = zuv(3)-p11-4.0D+00*p14
      p12 =  3.0D+00*h1-h2
      p13 = -2.0D+00*h1+h2
      thus = atan2(d-c,b-a)-thxu
      thsv = thuv-thus
      aa =  sin(thsv)/lu
      bb = -cos(thsv)/lu
      cc =  sin(thus)/lv
      dd =  cos(thus)/lv
      ac = aa*cc
      ad = aa*dd
      bc = bb*cc
      g1 = aa*ac*(3.0D+00*bc+2.0D+00*ad)
      g2 = cc*ac*(3.0D+00*ad+2.0D+00*bc)
      h1 = -aa*aa*aa*(5.0D+00*aa*bb*p50+(4.0D+00*bc+ad)*p41) &
           -cc*cc*cc*(5.0D+00*cc*dd*p05+(4.0D+00*ad+bc)*p14)
      h2 = 0.5D+00*zvv(2)-p02-p12
      h3 = 0.5D+00*zuu(3)-p20-p21
      p22 = (g1*h2+g2*h3-h1)/(g1+g2)
      p32 = h2-p22
      p23 = h3-p22
      itpv = it0
!
!  Convert XII and YII to UV system.
!
   30 continue

      dx = xii-x0
      dy = yii-y0
      u = ap*dx+bp*dy
      v = cp*dx+dp*dy
!
!  Evaluate the polynomial.
!
   31 continue

      p0 = p00+v*(p01+v*(p02+v*(p03+v*(p04+v*p05))))
      p1 = p10+v*(p11+v*(p12+v*(p13+v*p14)))
      p2 = p20+v*(p21+v*(p22+v*p23))
      p3 = p30+v*(p31+v*p32)
      p4 = p40+v*p41
      zii = p0+u*(p1+u*(p2+u*(p3+u*(p4+u*p5))))
      return
!
!  Calculation of zii by extrapolation in the rectangle.
!  Check if the necessary coefficients have been calculated.
!
   40 continue

      if(it0==itpv)     go to 50
!
!  Load coordinate and partial derivative values at the end
!  points of the border line segment.
!
   41 continue

      jipl = 3*(il1-1)
      jpd = 0

      do i = 1,2

        jipl = jipl+1
        idp = ipl(jipl)
        x(i) = xd(idp)
        y(i) = yd(idp)
        z(i) = zd(idp)
        jpdd = 5*(idp-1)

        do kpd = 1,5
          jpd = jpd+1
          jpdd = jpdd+1
          pd(jpd) = pdd(jpdd)
        end do

      end do
!
!  Determine the coefficients for the coordinate system
!  transformation from the XY system to the UV system
!  and vice versa.
!
   44 continue

      x0 = x(1)
      y0 = y(1)
      a = y(2)-y(1)
      b = x(2)-x(1)
      c = -b
      d = a
      ad = a*d
      bc = b*c
      dlt = ad-bc
      ap =  d/dlt
      bp = -b/dlt
      cp = -bp
      dp =  ap
!
!  Convert the partial derivatives at the end points of the
!  border line segment for the UV coordinate system.
!
   45 continue

      aa = a*a
      act2 = 2.0D+00*a*c
      cc = c*c
      ab = a*b
      adbc = ad+bc
      cd = c*d
      bb = b*b
      bdt2 = 2.0D+00*b*d
      dd = d*d

      do i = 1,2
        jpd = 5*i
        zu(i) = a*pd(jpd-4)+c*pd(jpd-3)
        zv(i) = b*pd(jpd-4)+d*pd(jpd-3)
        zuu(i) = aa*pd(jpd-2)+act2*pd(jpd-1)+cc*pd(jpd)
        zuv(i) = ab*pd(jpd-2)+adbc*pd(jpd-1)+cd*pd(jpd)
        zvv(i) = bb*pd(jpd-2)+bdt2*pd(jpd-1)+dd*pd(jpd)
      end do
!
!  Calculate the coefficients of the polynomial.
!
   47 continue

      p00 = z(1)
      p10 = zu(1)
      p01 = zv(1)
      p20 = 0.5D+00*zuu(1)
      p11 = zuv(1)
      p02 = 0.5D+00*zvv(1)
      h1 = z(2)-p00-p01-p02
      h2 = zv(2)-p01-zvv(1)
      h3 = zvv(2)-zvv(1)
      p03 =  10.0D+00*h1-4.0D+00*h2+0.5D+00*h3
      p04 = -15.0D+00*h1+7.0D+00*h2        -h3
      p05 =   6.0D+00*h1-3.0D+00*h2+0.5D+00*h3
      h1 = zu(2)-p10-p11
      h2 = zuv(2)-p11
      p12 =  3.0D+00*h1-h2
      p13 = -2.0D+00*h1+h2
      p21 = 0.0D+00
      p23 = -zuu(2)+zuu(1)
      p22 = -1.5D+00*p23
      itpv = it0
!
!  Convert XII and YII to UV system.
!
   50 continue

      dx = xii-x0
      dy = yii-y0
      u = ap*dx+bp*dy
      v = cp*dx+dp*dy
!
!  Evaluate the polynomial.
!
   51 continue

      p0 = p00+v*(p01+v*(p02+v*(p03+v*(p04+v*p05))))
      p1 = p10+v*(p11+v*(p12+v*p13))
      p2 = p20+v*(p21+v*(p22+v*p23))
      zii = p0+u*(p1+u*p2)
      return
!
!  Calculation of ZII by extrapolation in the triangle.
!  Check if the necessary coefficients have been calculated.
!
   60 continue

      if ( it0 == itpv ) then
        go to 70
      end if
!
!  Load coordinate and partial derivative values at the vertex
!  of the triangle.
!
   61 jipl = 3*il2-2
      idp = ipl(jipl)
      x(1) = xd(idp)
      y(1) = yd(idp)
      z(1) = zd(idp)
      jpdd = 5*(idp-1)

      do kpd = 1,5
        jpdd = jpdd+1
        pd(kpd) = pdd(jpdd)
      end do
!
!  Calculate the coefficients of the polynomial.
!
      p00 = z(1)
      p10 = pd(1)
      p01 = pd(2)
      p20 = 0.5D+00*pd(3)
      p11 = pd(4)
      p02 = 0.5D+00*pd(5)
      itpv = it0
!
!  Convert XII and YII to UV system.
!
   70 continue

      u = xii-x(1)
      v = yii-y(1)
!
!  Evaluate the polynomial.
!
      p0 = p00+v*(p01+v*p02)
      p1 = p10+v*p11
      zii = p0+u*(p1+u*p20)

      return
      end
      subroutine idsfft ( md, ncp, ndp, xd, yd, zd, nxi, nyi, xi, yi, zi, &
        iwk, wk )

!*****************************************************************************80
!
!! IDSFFT fits a smooth surface Z(X,Y) given irregular (X,Y,Z) data.
!
!  Discussion:
!
!    This subroutine performs smooth surface fitting when the projections
!    of the data points in the XY plane are irregularly
!    distributed in the plane.
!
!    The very first call to this subroutine and the call with a new
!    NCP value, a new NDP value, or new contents of the XD and
!    YD arrays must be made with MD = 1.  The call with MD = 2 must be
!    preceded by another call with the same NCP and NDP values and
!    with the same contents of the XD and YD arrays.  The call with
!    MD = 3 must be preceded by another call with the same NCP, NDP,
!    NXI, and NYI values and with the same contents of the XD, YD,
!    XI, and YI arrays.  Between the call with MD = 2 or MD = 3 and its
!    preceding call, the IWK and WK arrays must not be disturbed.
!    Use of a value between 3 and 5 (inclusive) for NCP is recommended
!    unless there is evidence that dictate otherwise.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2007
!
!  Author:
!
!    Original FORTRAN77 version by Hiroshi Akima.
!    FORTRAN90 version by John Burkardt.
!
!  References:
!
!    Hiroshi Akima,
!    Algorithm 526:
!    A Method of Bivariate Interpolation and Smooth Surface Fitting
!    for Values Given at Irregularly Distributed Points,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, pages 160-164.
!
!    Hiroshi Akima,
!    On Estimating Partial Derivatives for Bivariate Interpolation
!    of Scattered Data,
!    Rocky Mountain Journal of Mathematics,
!    Volume 14, Number 1, Winter 1984, pages 41-51.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MD, mode of computation (must be 1, 2, or 3,
!    else an error return will occur).
!    * 1, if this is the first call to this routine, or if the value of
!    NDP has been changed from the previous call, or if the contents of
!    the XD or YD arrays have been changed from the previous call.
!    * 2, if the values of NDP and the XD, YD arrays are unchanged from
!    the previous call, but new values for XI, YI are being used.  If
!    MD = 2 and NDP has been changed since the previous call to IDSFFT,
!    an error return occurs.
!    * 3, if the values of NDP, NXI, NYI, XD, YD, XI, YI are unchanged
!    from the previous call, i.e. if the only change on input to idsfft
!    is in the ZD array.  If MD = 3 and NDP, nxi or nyi has been changed
!    since the previous call to IDSFFT, an error return occurs.
!    Between the call with MD = 2 or MD = 3 and the preceding call, the
!    IWK and WK work arrays should not be disturbed.
!
!    Input, integer ( kind = 4 ) NCP, the number of data points used to
!    estimate partial derivatives.  2 <= NCP < NDP.  A value between 3 and
!    5 is recommended.
!
!    Input, integer ( kind = 4 ) NDP, the number of data points.  NDP must be
!    at least 4.
!
!    Input, real ( kind = 8 ) XD(NDP), Y(NDP), Z(NDP), the X, Y and Z
!    coordinates of the data points.
!
!    Input, integer ( kind = 4 ) NXI, NYI, the number of output grid points in the
!    X and Y directions.  NXI and NYI must each be at least 1.
!
!    Input, real ( kind = 8 ) XI(NXI), YI(NYI), the X and Y coordinates
!    of the grid points.
!
!    Output, real ( kind = 8 ) ZI(NXI,NYI), contains the interpolated Z
!    values at the grid points.
!
!    Workspace, integer IWK(max(31,27+NCP)*NDP+NXI*NYI).
!
!    Workspace, real ( kind = 8 ) WK(5*NDP).
!
  implicit none

  integer ( kind = 4 ) ncp
  integer ( kind = 4 ) ndp
  integer ( kind = 4 ) nxi
  integer ( kind = 4 ) nyi

  integer ( kind = 4 ) il1
  integer ( kind = 4 ) il2
  integer ( kind = 4 ) iti
  integer ( kind = 4 ) itpv
  integer ( kind = 4 ) iwk(max(31,27+ncp)*ndp+nxi*nyi)
  integer ( kind = 4 ) ixi
  integer ( kind = 4 ) iyi
  integer ( kind = 4 ) izi
  integer ( kind = 4 ) jig0mx
  integer ( kind = 4 ) jig0mn
  integer ( kind = 4 ) jig1mx
  integer ( kind = 4 ) jig1mn
  integer ( kind = 4 ) jigp
  integer ( kind = 4 ) jngp
  integer ( kind = 4 ) jwigp
  integer ( kind = 4 ) jwigp0
  integer ( kind = 4 ) jwipc
  integer ( kind = 4 ) jwipl
  integer ( kind = 4 ) jwipt
  integer ( kind = 4 ) jwiwl
  integer ( kind = 4 ) jwiwp
  integer ( kind = 4 ) jwngp
  integer ( kind = 4 ) jwngp0
  integer ( kind = 4 ) md
  integer ( kind = 4 ) md0
  integer ( kind = 4 ) ncp0
  integer ( kind = 4 ) ncppv
  integer ( kind = 4 ) ndp0
  integer ( kind = 4 ) ndppv
  integer ( kind = 4 ) ngp0
  integer ( kind = 4 ) ngp1
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nngp
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nxi0
  integer ( kind = 4 ) nxipv
  integer ( kind = 4 ) nyi0
  integer ( kind = 4 ) nyipv
  real ( kind = 8 ) wk(5*ndp)
  real ( kind = 8 ) xd(ndp)
  real ( kind = 8 ) xi(nxi)
  real ( kind = 8 ) yd(ndp)
  real ( kind = 8 ) yi(nyi)
  real ( kind = 8 ) zd(ndp)
  real ( kind = 8 ) zi(nxi*nyi)

  save /idpi/

  common /idpi/ itpv
!
!  Setting of some input parameters to local variables, for MD = 1, 2, 3.
!
  md0 = md
  ncp0 = ncp
  ndp0 = ndp
  nxi0 = nxi
  nyi0 = nyi
!
!  Error check, for MD = 1, 2, 3.
!
      if(md0<1.or.md0>3)           go to 90
      if(ncp0<2.or.ncp0>=ndp0)      go to 90
      if(ndp0<4)                      go to 90
      if(nxi0<1.or.nyi0<1)         go to 90

      if(md0>=2)        go to 21
      iwk(1) = ncp0
      iwk(2) = ndp0
      go to 22

   21 ncppv = iwk(1)
      ndppv = iwk(2)
      if(ncp0/=ncppv)   go to 90
      if(ndp0/=ndppv)   go to 90

   22 continue

      if(md0>=3)        go to 23
      iwk(3) = nxi0
      iwk(4) = nyi0
      go to 30
   23 nxipv = iwk(3)
      nyipv = iwk(4)
      if(nxi0/=nxipv)   go to 90
      if(nyi0/=nyipv)   go to 90
!
!  Allocation of storage areas in the IWK array, for MD = 1,2,3.
!
   30 jwipt = 16
      jwiwl = 6*ndp0+1
      jwngp0 = jwiwl-1
      jwipl = 24*ndp0+1
      jwiwp = 30*ndp0+1
      jwipc = 27*ndp0+1
      jwigp0 = max ( 31, 27 + ncp0 ) * ndp0
!
!  Triangulate the XY plane, for MD = 1.
!
   40 continue

      if ( md0 == 1 ) then

        call idtang(ndp0,xd,yd,nt,iwk(jwipt),nl,iwk(jwipl), &
          iwk(jwiwl),iwk(jwiwp),wk)

        iwk(5) = nt
        iwk(6) = nl

        if ( nt == 0 ) then
          return
        end if

      end if
!
!  Determine NCP points closest to each data point, for MD = 1.
!
      if ( md0 == 1 ) then

        call idcldp(ndp0,xd,yd,ncp0,iwk(jwipc))

        if ( iwk(jwipc) == 0 ) then
          return
        end if

      end if
!
!  Sort output grid points in ascending order of the triangle
!  number and the border line segment number, for md = 1,2.
!
      if ( md0 /= 3 ) then
        call idgrid(xd,yd,nt,iwk(jwipt),nl,iwk(jwipl),nxi0,nyi0, &
          xi,yi,iwk(jwngp0+1),iwk(jwigp0+1))
      end if
!
!  Estimate partial derivatives at all data points,
!  for md = 1,2,3.
!
      call idpdrv(ndp0,xd,yd,zd,ncp0,iwk(jwipc),wk)
!
!  Interpolate the ZI values, for md = 1,2,3.
!
   80 itpv = 0
      jig0mx = 0
      jig1mn = nxi0*nyi0+1
      nngp = nt+2*nl

      do jngp = 1,nngp

        iti = jngp

        if ( nt < jngp ) then
          il1 = (jngp-nt+1)/2
          il2 = (jngp-nt+2)/2
          if(il2>nl)     il2 = 1
          iti = il1*(nt+nl)+il2
        end if

        jwngp = jwngp0+jngp
        ngp0 = iwk(jwngp)

        if ( 0 < ngp0 ) then
          jig0mn = jig0mx+1
          jig0mx = jig0mx+ngp0

          do jigp = jig0mn,jig0mx

            jwigp = jwigp0+jigp
            izi = iwk(jwigp)
            iyi = (izi-1)/nxi0+1
            ixi = izi-nxi0*(iyi-1)

            call idptip(xd,yd,zd,nt,iwk(jwipt),nl,iwk(jwipl),wk, &
              iti,xi(ixi),yi(iyi),zi(izi))

          end do

        end if

        jwngp = jwngp0+2*nngp+1-jngp
        ngp1 = iwk(jwngp)

        if ( 0 < ngp1 ) then
          jig1mx = jig1mn-1
          jig1mn = jig1mn-ngp1

          do jigp = jig1mn,jig1mx

            jwigp = jwigp0+jigp
            izi = iwk(jwigp)
            iyi = (izi-1)/nxi0+1
            ixi = izi-nxi0*(iyi-1)
            call idptip(xd,yd,zd,nt,iwk(jwipt),nl,iwk(jwipl),wk, &
              iti,xi(ixi),yi(iyi),zi(izi))
          end do

        end if

      end do

      return

   90 write ( *, 2090 ) md0,ncp0,ndp0,nxi0,nyi0
      return
 2090 format(1x/' ***   improper input parameter value(s).'/ &
         '   md  = ',i4,10x,'ncp  = ',i6,10x,'ndp  = ',i6, &
         10x,'nxi  = ',i6,10x,'nyi  = ',i6/ &
         ' error detected in routine   idsfft'/)
      end
      subroutine idtang ( ndp, xd, yd, nt, ipt, nl, ipl, iwl, iwp, wk )

!*****************************************************************************80
!
!! IDTANG performs triangulation.
!
!  Discussion:
!
!    This subroutine performs triangulation.  It divides the XY
!    plane into a number of triangles according to given data
!    points in the plane, determines line segments that form the
!    border of data area, and determines the triangle numbers
!    corresponding to the border line segments.
!
!    At completion, point numbers of the vertexes of each triangle
!    are listed counter-clockwise.  Point numbers of the end points
!    of each border line segment are listed counter-clockwise,
!    listing order of the line segments being counter-clockwise.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2007
!
!  Author:
!
!    Original FORTRAN77 version by Hiroshi Akima.
!    FORTRAN90 version by John Burkardt.
!
!  References:
!
!    Hiroshi Akima,
!    Algorithm 526:
!    A Method of Bivariate Interpolation and Smooth Surface Fitting
!    for Values Given at Irregularly Distributed Points,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, pages 160-164.
!
!    Hiroshi Akima,
!    On Estimating Partial Derivatives for Bivariate Interpolation
!    of Scattered Data,
!    Rocky Mountain Journal of Mathematics,
!    Volume 14, Number 1, Winter 1984, pages 41-51.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDP, the number of data points.
!
!    Input, real ( kind = 8 ) XD(NDP), Y(NDP), the X and Y coordinates
!    of the data points.
!
!    Output, integer ( kind = 4 ) NT, the number of triangles,
!
!    Output, integer ( kind = 4 ) IPT(6*NDP-15), where the point numbers of the
!    vertexes of the IT-th triangle are to be stored as entries
!    3*IT-2, 3*IT-1, and 3*IT, for IT = 1 to NT.
!
!    Output, integer ( kind = 4 ) NL, the number of border line segments.
!
!    Output, integer ( kind = 4 ) IPL(6*NDP), where the point numbers of the end
!    points of the (il)th border line segment and its respective triangle
!    number are to be stored as the (3*il-2)nd, (3*il-1)st, and (3*il)th
!    elements, il = 1,2,..., nl.
!
!    Workspace, integer IWL(18*NDP),
!
!    Workspace, integer IWP(NDP),
!
!    Workspace, real ( kind = 8 ) WK(NDP).
!
  implicit none

  integer ( kind = 4 ) ndp

  real ( kind = 8 ) ar
  real ( kind = 8 ) armn
  real ( kind = 8 ) armx
  real ( kind = 8 ) dsq12
  real ( kind = 8 ) dsqf
  real ( kind = 8 ) dsqi
  real ( kind = 8 ) dsqmn
  real ( kind = 8 ) dsqmx
  real ( kind = 8 ) dx
  real ( kind = 8 ) dx21
  real ( kind = 8 ) dxmn
  real ( kind = 8 ) dxmx
  real ( kind = 8 ) dy
  real ( kind = 8 ) dy21
  real ( kind = 8 ) dymn
  real ( kind = 8 ) dymx
  integer ( kind = 4 ) idxchg
  integer ( kind = 4 ) ilf
  integer ( kind = 4 ) ilft2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip1p1
  integer ( kind = 4 ) ip2
  integer ( kind = 4 ) ip3
  integer ( kind = 4 ) ipl(6*ndp)
  integer ( kind = 4 ) ipl1
  integer ( kind = 4 ) ipl2
  integer ( kind = 4 ) iplj1
  integer ( kind = 4 ) iplj2
  integer ( kind = 4 ) ipmn1
  integer ( kind = 4 ) ipmn2
  integer ( kind = 4 ) ipt(6*ndp-15)
  integer ( kind = 4 ) ipt1
  integer ( kind = 4 ) ipt2
  integer ( kind = 4 ) ipt3
  integer ( kind = 4 ) ipti
  integer ( kind = 4 ) ipti1
  integer ( kind = 4 ) ipti2
  integer ( kind = 4 ) irep
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it1t3
  integer ( kind = 4 ) it2t3
  integer ( kind = 4 ) itf(2)
  integer ( kind = 4 ) its
  integer ( kind = 4 ) itt3
  integer ( kind = 4 ) itt3r
  integer ( kind = 4 ) iwl(18*ndp)
  integer ( kind = 4 ) iwp(ndp)
  integer ( kind = 4 ) jlt3
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) jp2
  integer ( kind = 4 ) jp2t3
  integer ( kind = 4 ) jp3t3
  integer ( kind = 4 ) jpc
  integer ( kind = 4 ) jpmn
  integer ( kind = 4 ) jpmx
  integer ( kind = 4 ) jwl
  integer ( kind = 4 ) jwl1
  integer ( kind = 4 ) jwl1mn
  integer ( kind = 4 ) ndp0
  integer ( kind = 4 ) ndpm1
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nl0
  integer ( kind = 4 ) nlf
  integer ( kind = 4 ) nlfc
  integer ( kind = 4 ) nlft2
  integer ( kind = 4 ) nln
  integer ( kind = 4 ) nlnt3
  integer ( kind = 4 ) nlt3
  integer ( kind = 4 ), parameter :: nrep = 100
  integer ( kind = 4 ) nsh
  integer ( kind = 4 ) nsht3
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nt0
  integer ( kind = 4 ) ntf
  integer ( kind = 4 ) ntt3
  integer ( kind = 4 ) ntt3p3
  real ( kind = 8 ), parameter :: ratio = 1.0D-06
  real ( kind = 8 ) side
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) u3
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) v3
  real ( kind = 8 ) wk(ndp)
  real ( kind = 8 ) x1
  real ( kind = 8 ) xd(ndp)
  real ( kind = 8 ) xdmp
  real ( kind = 8 ) y1
  real ( kind = 8 ) yd(ndp)
  real ( kind = 8 ) ydmp
!
!  Statement functions
!
      dsqf(u1,v1,u2,v2) = (u2-u1)**2+(v2-v1)**2
      side(u1,v1,u2,v2,u3,v3) = (v3-v1)*(u2-u1)-(u3-u1)*(v2-v1)
!
!  Preliminary processing.
!
   10 ndp0 = ndp
      ndpm1 = ndp0-1
      if(ndp0<4)       go to 90
!
!  Determine the closest pair of data points and their midpoint.
!
   20 dsqmn = dsqf(xd(1),yd(1),xd(2),yd(2))
      ipmn1 = 1
      ipmn2 = 2

      do ip1 = 1, ndpm1

        x1 = xd(ip1)
        y1 = yd(ip1)
        ip1p1 = ip1+1

        do ip2 = ip1p1,ndp0

          dsqi = dsqf(x1,y1,xd(ip2),yd(ip2))
          if ( dsqi == 0.0D+00 )      go to 91

          if ( dsqi < dsqmn ) then
            dsqmn = dsqi
            ipmn1 = ip1
            ipmn2 = ip2
          end if

        end do

      end do

      dsq12 = dsqmn
      xdmp = (xd(ipmn1)+xd(ipmn2))/2.0D+00
      ydmp = (yd(ipmn1)+yd(ipmn2))/2.0D+00
!
!  Sort the other (NDP-2) data points in ascending order of
!  distance from the midpoint and store the sorted data point
!  numbers in the IWP array.
!
   30 jp1 = 2

      do ip1 = 1,ndp0
        if ( ip1 /= ipmn1 .and. ip1 /= ipmn2 ) then
          jp1 = jp1+1
          iwp(jp1) = ip1
          wk(jp1) = dsqf(xdmp,ydmp,xd(ip1),yd(ip1))
        end if
      end do

      do jp1 = 3,ndpm1

        dsqmn = wk(jp1)
        jpmn = jp1

        do jp2 = jp1,ndp0
          if ( wk(jp2) < dsqmn ) then
            dsqmn = wk(jp2)
            jpmn = jp2
          end if
        end do

        its = iwp(jp1)
        iwp(jp1) = iwp(jpmn)
        iwp(jpmn) = its
        wk(jpmn) = wk(jp1)

      end do
!
!  If necessary, modify the ordering in such a way that the
!  first three data points are not collinear.
!
   35 ar = dsq12*ratio
      x1 = xd(ipmn1)
      y1 = yd(ipmn1)
      dx21 = xd(ipmn2)-x1
      dy21 = yd(ipmn2)-y1

      do jp = 3,ndp0
        ip = iwp(jp)
        if(abs((yd(ip)-y1)*dx21-(xd(ip)-x1)*dy21)>ar) &
                     go to 37
      end do

      go to 92
   37 continue

      if ( jp /= 3 ) then

        jpmx = jp
        jp = jpmx+1

        do jpc = 4,jpmx
          jp = jp-1
          iwp(jp) = iwp(jp-1)
        end do

        iwp(3) = ip

      end if
!
!  Form the first triangle.  Store point numbers of the vertices
!  of the triangle in the IPT array, and store point numbers
!  of the border line segments and the triangle number in
!  the IPL array.
!
      ip1 = ipmn1
      ip2 = ipmn2
      ip3 = iwp(3)

      if(side(xd(ip1),yd(ip1),xd(ip2),yd(ip2),xd(ip3),yd(ip3)) &
           >=0.0D+00)       go to 41

      ip1 = ipmn2
      ip2 = ipmn1

   41 continue

      nt0 = 1
      ntt3 = 3
      ipt(1) = ip1
      ipt(2) = ip2
      ipt(3) = ip3
      nl0 = 3
      nlt3 = 9
      ipl(1) = ip1
      ipl(2) = ip2
      ipl(3) = 1
      ipl(4) = ip2
      ipl(5) = ip3
      ipl(6) = 1
      ipl(7) = ip3
      ipl(8) = ip1
      ipl(9) = 1
!
!  Add the remaining (NDP-3) data points, one by one.
!
   50 continue

      do 79  jp1 = 4,ndp0

        ip1 = iwp(jp1)
        x1 = xd(ip1)
        y1 = yd(ip1)
!
!  Determine the visible border line segments.
!
        ip2 = ipl(1)
        jpmn = 1
        dxmn = xd(ip2)-x1
        dymn = yd(ip2)-y1
        dsqmn = dxmn**2+dymn**2
        armn = dsqmn*ratio
        jpmx = 1
        dxmx = dxmn
        dymx = dymn
        dsqmx = dsqmn
        armx = armn

        do 52  jp2 = 2,nl0
          ip2 = ipl(3*jp2-2)
          dx = xd(ip2)-x1
          dy = yd(ip2)-y1
          ar = dy*dxmn-dx*dymn
          if(ar>armn)       go to 51
          dsqi = dx**2+dy**2
          if(ar>=(-armn).and.dsqi>=dsqmn)      go to 51
          jpmn = jp2
          dxmn = dx
          dymn = dy
          dsqmn = dsqi
          armn = dsqmn*ratio
   51     ar = dy*dxmx-dx*dymx
          if(ar<(-armx))    go to 52
          dsqi = dx**2+dy**2
          if(ar<=armx.and.dsqi>=dsqmx)    go to 52
          jpmx = jp2
          dxmx = dx
          dymx = dy
          dsqmx = dsqi
          armx = dsqmx*ratio

   52   continue

        if(jpmx<jpmn) then
          jpmx = jpmx+nl0
        end if

        nsh = jpmn-1
        if(nsh<=0)      go to 60
!
!  Shift (rotate) the IPL array to have the invisible border
!  line segments contained in the first part of the IPL array.
!
        nsht3 = nsh*3

        do jp2t3 = 3,nsht3,3
          jp3t3 = jp2t3+nlt3
          ipl(jp3t3-2) = ipl(jp2t3-2)
          ipl(jp3t3-1) = ipl(jp2t3-1)
          ipl(jp3t3)   = ipl(jp2t3)
        end do

        do jp2t3 = 3,nlt3,3
          jp3t3 = jp2t3+nsht3
          ipl(jp2t3-2) = ipl(jp3t3-2)
          ipl(jp2t3-1) = ipl(jp3t3-1)
          ipl(jp2t3)   = ipl(jp3t3)
        end do

        jpmx = jpmx-nsh
!
!  Add triangles to the IPT array, updates border line
!  segments in the IPL array, and sets flags for the border
!  line segments to be reexamined in the IWL array.
!
   60   jwl = 0

        do 64  jp2 = jpmx,nl0

          jp2t3 = jp2*3
          ipl1 = ipl(jp2t3-2)
          ipl2 = ipl(jp2t3-1)
          it   = ipl(jp2t3)
!
!  Add a triangle to the IPT array.
!
          nt0 = nt0+1
          ntt3 = ntt3+3
          ipt(ntt3-2) = ipl2
          ipt(ntt3-1) = ipl1
          ipt(ntt3)   = ip1
!
!  Update border line segments in the IPL array.
!
          if(jp2/=jpmx)      go to 61
          ipl(jp2t3-1) = ip1
          ipl(jp2t3)   = nt0
   61     if(jp2/=nl0)       go to 62
          nln = jpmx+1
          nlnt3 = nln*3
          ipl(nlnt3-2) = ip1
          ipl(nlnt3-1) = ipl(1)
          ipl(nlnt3)   = nt0
!
!  Determine the vertex that does not lie on the border
!  line segments.
!
   62     itt3 = it*3
          ipti = ipt(itt3-2)
          if(ipti/=ipl1.and.ipti/=ipl2)   go to 63
          ipti = ipt(itt3-1)
          if(ipti/=ipl1.and.ipti/=ipl2)   go to 63
          ipti = ipt(itt3)
!
!  Check if the exchange is necessary.
!
   63     if(idxchg(xd,yd,ip1,ipti,ipl1,ipl2)==0)     go to 64
!
!  Modify the IPT array when necessary.
!
          ipt(itt3-2) = ipti
          ipt(itt3-1) = ipl1
          ipt(itt3)   = ip1
          ipt(ntt3-1) = ipti
          if(jp2==jpmx)      ipl(jp2t3) = it
          if(jp2==nl0.and.ipl(3)==it)     ipl(3) = nt0
!
!  Set flags in the IWL array.
!
          jwl = jwl+4
          iwl(jwl-3) = ipl1
          iwl(jwl-2) = ipti
          iwl(jwl-1) = ipti
          iwl(jwl)   = ipl2

   64   continue

        nl0 = nln
        nlt3 = nlnt3
        nlf = jwl/2
        if(nlf==0)      go to 79
!
!  Improve triangulation.
!
   70   ntt3p3 = ntt3+3

        do 78  irep = 1,nrep

          do 76  ilf = 1,nlf

            ilft2 = ilf*2
            ipl1 = iwl(ilft2-1)
            ipl2 = iwl(ilft2)
!
!  Locate in the IPT array two triangles on both sides of
!  the flagged line segment.
!
            ntf = 0

            do 71  itt3r = 3,ntt3,3
              itt3 = ntt3p3-itt3r
              ipt1 = ipt(itt3-2)
              ipt2 = ipt(itt3-1)
              ipt3 = ipt(itt3)
              if(ipl1/=ipt1.and.ipl1/=ipt2.and. &
                 ipl1/=ipt3)      go to 71
              if(ipl2/=ipt1.and.ipl2/=ipt2.and. &
                 ipl2/=ipt3)      go to 71
              ntf = ntf+1
              itf(ntf) = itt3/3
              if(ntf==2)     go to 72
   71       continue

            if(ntf<2)       go to 76
!
!  Determine the vertexes of the triangles that do not lie
!  on the line segment.
!
   72       it1t3 = itf(1)*3
            ipti1 = ipt(it1t3-2)
            if(ipti1/=ipl1.and.ipti1/=ipl2)    go to 73
            ipti1 = ipt(it1t3-1)
            if(ipti1/=ipl1.and.ipti1/=ipl2)    go to 73
            ipti1 = ipt(it1t3)
   73       it2t3 = itf(2)*3
            ipti2 = ipt(it2t3-2)
            if(ipti2/=ipl1.and.ipti2/=ipl2)    go to 74
            ipti2 = ipt(it2t3-1)
            if(ipti2/=ipl1.and.ipti2/=ipl2)    go to 74
            ipti2 = ipt(it2t3)
!
!  Check if the exchange is necessary.
!
   74       if(idxchg(xd,yd,ipti1,ipti2,ipl1,ipl2)==0) then
              go to 76
            end if
!
!  Modify the IPT array when necessary.
!
            ipt(it1t3-2) = ipti1
            ipt(it1t3-1) = ipti2
            ipt(it1t3)   = ipl1
            ipt(it2t3-2) = ipti2
            ipt(it2t3-1) = ipti1
            ipt(it2t3)   = ipl2
!
!  Set new flags.
!
            jwl = jwl+8
            iwl(jwl-7) = ipl1
            iwl(jwl-6) = ipti1
            iwl(jwl-5) = ipti1
            iwl(jwl-4) = ipl2
            iwl(jwl-3) = ipl2
            iwl(jwl-2) = ipti2
            iwl(jwl-1) = ipti2
            iwl(jwl)   = ipl1

            do jlt3 = 3,nlt3,3
              iplj1 = ipl(jlt3-2)
              iplj2 = ipl(jlt3-1)
              if((iplj1==ipl1.and.iplj2==ipti2).or. &
                 (iplj2==ipl1.and.iplj1==ipti2)) &
                               ipl(jlt3) = itf(1)
              if((iplj1==ipl2.and.iplj2==ipti1).or. &
                 (iplj2==ipl2.and.iplj1==ipti1)) &
                               ipl(jlt3) = itf(2)
            end do

   76     continue
          nlfc = nlf
          nlf = jwl/2
          if(nlf==nlfc)      go to 79
!
!  Reset the IWL array for the next round.
!
          jwl = 0
          jwl1mn = (nlfc+1)*2
          nlft2 = nlf*2

          do jwl1 = jwl1mn,nlft2,2
            jwl = jwl+2
            iwl(jwl-1) = iwl(jwl1-1)
            iwl(jwl)   = iwl(jwl1)
          end do

          nlf = jwl/2
   78   continue
   79 continue
!
!  Rearrange the IPT array so that the vertexes of each triangle
!  are listed counter-clockwise.
!
   80 do itt3 = 3,ntt3,3

        ip1 = ipt(itt3-2)
        ip2 = ipt(itt3-1)
        ip3 = ipt(itt3)

        if(side(xd(ip1),yd(ip1),xd(ip2),yd(ip2),xd(ip3),yd(ip3)) &
              < 0.0D+00 ) then
          ipt(itt3-2) = ip2
          ipt(itt3-1) = ip1
        end if

      end do

      nt = nt0
      nl = nl0
      return
!
!  Error exit.
!
   90 write (*,2090)  ndp0
      go to 93
   91 write ( *,2091)  ndp0,ip1,ip2,x1,y1
      go to 93
   92 write ( *,2092)  ndp0
   93 write ( *,2093)
      nt = 0
      return

 2090 format(1x/' ***   ndp less than 4.'/'   ndp  = ',i5)
 2091 format(1x/' ***   identical data points.'/ &
         '   ndp  = ',i5,5x,'ip1  = ',i5,5x,'ip2  = ',i5, &
         5x,'xd  = ',e12.4,5x,'yd  = ',e12.4)
 2092 format(1x/' ***   all collinear data points.'/ &
         '   ndp  = ',i5)
 2093 format(' error detected in routine   idtang'/)
      end
      function idxchg ( x, y, i1, i2, i3, i4 )

!*****************************************************************************80
!
!! IDXCHG determines whether two triangles should be exchanged.
!
!  Discussion:
!
!    This function determines whether or not the exchange of two
!    triangles is necessary on the basis of max-min-angle criterion
!    by Charles Lawson.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2007
!
!  Author:
!
!    Original FORTRAN77 version by Hiroshi Akima.
!    FORTRAN90 version by John Burkardt.
!
!  References:
!
!    Hiroshi Akima,
!    Algorithm 526:
!    A Method of Bivariate Interpolation and Smooth Surface Fitting
!    for Values Given at Irregularly Distributed Points,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, pages 160-164.
!
!    Hiroshi Akima,
!    On Estimating Partial Derivatives for Bivariate Interpolation
!    of Scattered Data,
!    Rocky Mountain Journal of Mathematics,
!    Volume 14, Number 1, Winter 1984, pages 41-51.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(*), Y(*), the coordinates of the data.
!
!    Input, integer ( kind = 4 ) I1, I2, I3, I4, indices of four points
!    P1, P2, P3, and P4 that form a quadrilateral with P3 and P4 connected
!    diagonally.
!
!    Output, integer ( kind = 4 ) IDXCHG, is 1 if an exchange is recommended,
!    and 0 otherwise.
!
  implicit none

  real ( kind = 8 ) a1sq
  real ( kind = 8 ) a2sq
  real ( kind = 8 ) a3sq
  real ( kind = 8 ) a4sq
  real ( kind = 8 ) b1sq
  real ( kind = 8 ) b2sq
  real ( kind = 8 ) b3sq
  real ( kind = 8 ) b4sq
  real ( kind = 8 ) c1sq
  real ( kind = 8 ) c2sq
  real ( kind = 8 ) c3sq
  real ( kind = 8 ) c4sq
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) idx
  integer ( kind = 4 ) idxchg
  real ( kind = 8 ) s1sq
  real ( kind = 8 ) s2sq
  real ( kind = 8 ) s3sq
  real ( kind = 8 ) s4sq
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) u3
  real ( kind = 8 ) u4
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
!
!  Preliminary processing.
!
  x1 = x(i1)
  y1 = y(i1)
  x2 = x(i2)
  y2 = y(i2)
  x3 = x(i3)
  y3 = y(i3)
  x4 = x(i4)
  y4 = y(i4)
!
!  Calculation.
!
  idx = 0
  u3 = ( y2 - y3 ) * ( x1 - x3 ) - ( x2 - x3 ) * ( y1 - y3 )
  u4 = ( y1 - y4 ) * ( x2 - x4 ) - ( x1 - x4 ) * ( y2 - y4 )

  if ( 0.0D+00 < u3 * u4 ) then

    u1 = (y3-y1)*(x4-x1)-(x3-x1)*(y4-y1)
    u2 = (y4-y2)*(x3-x2)-(x4-x2)*(y3-y2)

    a1sq = (x1-x3)**2+(y1-y3)**2
    b3sq = a1sq
    b1sq = (x4-x1)**2+(y4-y1)**2
    a4sq = b1sq
    c1sq = (x3-x4)**2+(y3-y4)**2
    c2sq = c1sq
    a2sq = (x2-x4)**2+(y2-y4)**2
    b4sq = a2sq
    b2sq = (x3-x2)**2+(y3-y2)**2
    a3sq = b2sq
    c3sq = (x2-x1)**2+(y2-y1)**2
    c4sq = c3sq

    s1sq = u1 * u1 / ( c1sq * max ( a1sq, b1sq ) )
    s2sq = u2 * u2 / ( c2sq * max ( a2sq, b2sq ) )
    s3sq = u3 * u3 / ( c3sq * max ( a3sq, b3sq ) )
    s4sq = u4 * u4 / ( c4sq * max ( a4sq, b4sq ) )

    if ( min ( s1sq, s2sq ) < min ( s3sq, s4sq ) ) then
      idx = 1
    end if

  end if

  idxchg = idx

  return
end
