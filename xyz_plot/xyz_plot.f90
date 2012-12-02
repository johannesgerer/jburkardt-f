program main

!*****************************************************************************80
!
!! MAIN is the main program for XYZ_PLOT.
!
!  Discussion:
!
!    XYZPLOT is an interactive graphics program.
!
!  Issues:
!
!    Convert Least Squares calculations to ORTPOL/ORTVAL if possible.
!
!    Clean up scatterplot triangulation.
!
!    CON doesn't seem to be set up.
!
!    "Dev = ps" fails, whereas "dev=ps" works.
!
!    Notes: There is something wrong with the CLIP routine.  I noticed this
!    because perfectly vertical lines which were entirely within the range
!    were getting clipped.  In exasperation, I simply iced out the clipping
!    for now.
!
!    3D axis, and labelling, are messed up.
!
!    Add CON routine.
!
!  Note:
!
!    MAXFIX specifies the maximum length of a formula that a user
!    can type in,in characters.
!
!    MAXRPN specifies the maximum length of a "compiled" RPN formula,
!    in symbols.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: BIG = 60000

  integer, parameter :: maxdat = BIG
  integer, parameter :: maxiwork = 1000
  integer, parameter :: maxrwork = 110
  integer, parameter :: maxtri = 1000
  integer, parameter :: val_max = BIG
  integer, parameter :: maxvl3 = BIG

  character ( len = 80 ) carray
  character ( len = 10 ) coord
  real dat(maxdat,9)
  character ( len = 10 ) dev
  character ( len = 80 ) filexyz
  integer i
  integer icom
  integer icurve
  integer idat(maxdat)
  integer ierror
  integer iplot
  integer iplt1
  integer iplt2
  character ( len = 10 ) isay
  integer ival(val_max)
  integer iwork(maxiwork)
  integer ixplt1
  integer ixplt2
  integer iyplt1
  integer iyplt2
  integer j
  integer k
  logical l3d
  logical ldat
  logical lframe
  integer m
  integer marray
  integer mdat
  integer ncon
  integer ndat
  integer nfine
  integer ntri
  integer nval
  integer nxgrid
  integer nygrid
  real plhite
  real plwide
  logical s_eqi
  real tdat(maxdat)
  real theta
  integer nodtri(3,maxtri)
  character ( len = 80 ) title
  integer tnbr(3,maxtri)
  real udat(maxdat)
  real umax
  real umin
  real uval
  real vdat(maxdat)
  real vscale
  real rwork(maxrwork)
  real x1
  real x2
  real x3
  real x3val(maxvl3)
  real x3dat(maxdat)
  real xdat(maxdat)
  real xf
  real xmax
  real xmaxw
  real xmin
  real xminw
  real xplt1
  real xplt2
  real xval(val_max)
  real y1
  real y2
  real y3
  real y3val(maxvl3)
  real y3dat(maxdat)
  real ydat(maxdat)
  real yf
  real ymax
  real ymaxw
  real ymin
  real yminw
  real yplt1
  real yplt2
  real yval(val_max)
  real z1
  real z2
  real z3
  real z3val(maxvl3)
  real z3dat(maxdat)
  real zdat(maxdat)
  real zf
!
!  anyplt stuff
!
  common /anycom/ iplt1,iplt2,ixplt1,ixplt2,iyplt1, &
                  iyplt2,marray,xplt1,xplt2,yplt1,yplt2
  common /anychr/ carray

  call timestamp ( )
!
!  Initialize.
!
  dev = '?'
  idat(1:maxdat) = 0
  ierror = 0
  iplot = 0
  ldat = .true.
  lframe = .false.
  mdat = 0
  ncon = 11
  nfine = 101
  nxgrid = 2
  nygrid = 2
  plhite = 1.0E+00
  plwide = 1.0E+00
  theta = 0.0E+00
  title = ' '
  vscale = 1.0E+00
  x3dat(1:maxdat) = 0.0E+00
  xdat(1:maxdat) = 0.0E+00
  xmax = 1.0E+00
  xmaxw = 1.0E+00
  xmin = 0.0E+00
  xminw = 0.0E+00
  x3val(1:val_max) = 0.0E+00
  xval(1:val_max) = 0.0E+00
  y3dat(1:maxdat) = 0.0E+00
  ydat(1:maxdat) = 0.0E+00
  ymax = 1.0E+00
  ymaxw = 1.0E+00
  ymin = 0.0E+00
  yminw = 0.0E+00
  y3val(1:val_max) = 0.0E+00
  yval(1:val_max) = 0.0E+00
  z3dat(1:maxdat) = 0.0E+00
  zdat(1:maxdat) = 0.0E+00
  z3val(1:val_max) = 0.0E+00
!
!  Say hello.
!
  call hello
!
!  Choose a device.
!
  call device_choose ( dev )
!
!  Set options.
!
10    continue

  call option ( dev, iplot, ldat, lframe, ncon, nfine, &
    nxgrid, nygrid, theta, title, vscale )

  if ( iplot == 0 ) then
    icom = 0
    carray = dev
    call anyplt ( icom )

    icom = 13
    call anyplt ( icom )
    write ( *, '(a)' ) trim ( carray )
  end if

  if ( iplot > 0 ) then
    go to 60
  end if
!
!  Initialize for this plot
!
20    continue

  iplot = iplot + 1
  nval = 0
  ndat = 0
  icurve = 0
!
!  Start here if a curve is to be added to the plot.
!
30    continue

  icurve = icurve + 1
!
!  Get the coordinate system option.
!
40    continue

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Coordinate choices:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CAU, CON, FX, FXY, LXY, LTXY, PWP, Q, RT, SCAT, '
  write ( *, '(a)' ) 'STAR, SXY, TXY, TXYZ, UVXY, XYS, XYZ, or ?: '
  write ( *, '(a)' ) 'Enter coordinate choice.'
  read ( *, '(a)' ) coord

  if ( s_eqi(coord,'CAU') ) then

    l3d = .false.

    call cau ( ival, val_max, nval, xval, yval )

  else if ( s_eqi(coord,'CON') ) then

    l3d = .false.
    mdat = 3

    call getdat ( dat, idat, ierror, maxdat, mdat, 'xyu', ndat )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) 'XYZ_PLOT - Fatal error!'
      go to 10
    end if

    xdat(1:ndat) = dat(1:ndat,1)
    ydat(1:ndat) = dat(1:ndat,2)
    udat(1:ndat) = dat(1:ndat,3)
!
!  Now triangulate the data.
!
    if ( 2*ndat > maxiwork ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Integer workspace needs exceeded.'
      go to 10
    end if

    call ivec_identity ( ndat, iwork )

    if ( 2*ndat > maxrwork ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Real workspace needs exceeded.'
      go to 10
    end if

    j = 0
    do i = 1, ndat
      j = j+1
      rwork(j) = xdat(i)
      j = j+1
      rwork(j) = ydat(i)
    end do

    call rtris2 ( ndat, ndat, rwork, iwork(1), ntri, nodtri, tnbr, &
      iwork(ndat+1), ierror )
!
!  Now call TRICON to get the levels.
!
    m = ncon + 2
    umin = minval ( udat(1:ndat) )
    umax = maxval ( udat(1:ndat) )

    do i = 2, m-1

      call rvec_even_select ( umin, umax, m, i, uval )

      call tricon ( ival, val_max, ndat, nodtri, ntri, nval, udat, &
        uval, xdat, xval, ydat, yval )

    end do

  else if ( s_eqi ( coord, 'FX' ) ) then

    l3d = .false.
    mdat = 2

    call getdat ( dat, idat, ierror, maxdat, mdat, 'xy', ndat )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) 'XYZ_PLOT - Fatal error!'
      go to 10
    end if

    nval = ndat

    xval(1:nval) = dat(1:nval,1)
    yval(1:nval) = dat(1:nval,2)

    ival(1:nval-1) = 1
    ival(nval) = 0

  else if ( s_eqi(coord,'FXY') ) then

    l3d = .false.

    call ffxy ( ierror, ival, val_max, nval, xval, yval )

  else if ( s_eqi(coord,'LXY') ) then

    l3d = .false.
    mdat = 2

    call getdat ( dat, idat, ierror, maxdat, mdat, 'xy', ndat )

    xdat(1:ndat) = dat(1:ndat,1)
    ydat(1:ndat) = dat(1:ndat,2)

    call lsqxy ( ierror, ival, maxrwork, val_max, ndat, nfine, &
      nval, rwork, xdat, xval, ydat, yval )

    if ( ierror /= 0 ) then
      go to 10
    end if

  else if ( s_eqi(coord,'LTXY') ) then

    l3d = .false.
    mdat = 2

    call getdat ( dat, idat, ierror, maxdat, mdat, 'xy', ndat )

    xdat(1:ndat) = dat(1:ndat,1)
    ydat(1:ndat) = dat(1:ndat,2)

    do i = 1, ndat
      tdat(i) = real ( i - 1 ) / real ( ndat - 1 )
    end do

    call lsqtxy ( ierror, ival, maxrwork, val_max, ndat, nfine, &
      nval, rwork, tdat, xdat, xval, ydat, yval )

    if ( ierror /= 0 ) then
      go to 10
    end if

  else if ( s_eqi(coord,'PWP') ) then

    l3d = .false.
    mdat = 2

    call getdat ( dat, idat, ierror, maxdat, mdat, 'xy', ndat )

    xdat(1:ndat) = dat(1:ndat,1)
    ydat(1:ndat) = dat(1:ndat,2)

    call pwp ( ierror, ival, maxrwork, val_max, ndat, nfine, &
      nval, rwork, xdat, xval, ydat, yval )

  else if ( s_eqi ( coord, 'RT' ) ) then

    l3d = .false.
    mdat = 2

    call getdat ( dat, idat, ierror, maxdat, mdat, 'rt', ndat )

    xdat(1:ndat) = dat(1:ndat,1) * cos ( dat(1:ndat,2) )
    ydat(1:ndat) = dat(1:ndat,1) * sin ( dat(1:ndat,2) )

    idat(1:ndat-1) = 1
    idat(ndat) = 0

    nval = ndat
    do i = 1, nval
      xval(i) = dat(i,1) * cos ( dat(i,2) )
      yval(i) = dat(i,1) * sin ( dat(i,2) )
    end do

    ival(1:nval-1) = 1
    ival(nval) = 0

  else if ( s_eqi(coord,'SCAT') ) then

    l3d = .false.
    mdat = 2

    call getdat ( dat, idat, ierror, maxdat, mdat, 'xy', ndat )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) 'XYZ_PLOT - Fatal error!'
      go to 10
    end if

    xdat(1:ndat) = dat(1:ndat,1)
    ydat(1:ndat) = dat(1:ndat,2)
!
!  Offer to triangulate the points.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Want to triangulate the points?'
    read ( *, '(a)') isay

    if ( s_eqi ( isay(1:1), 'Y' ) ) then

      if ( 2*ndat > maxiwork ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Integer workspace needs exceeded.'
        go to 10
      end if

      call ivec_identity ( ndat, iwork )

      if ( 2*ndat > maxrwork ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Real workspace needs exceeded.'
        go to 10
      end if

      j = 0
      do i = 1, ndat
        j = j+1
        rwork(j) = xdat(i)
        j = j+1
        rwork(j) = ydat(i)
      end do

      call rtris2 ( ndat, ndat, rwork, iwork(1), ntri, nodtri, tnbr, &
        iwork(ndat+1), ierror )

      if ( ierror /= 0 ) then
        go to 10
      end if

      do i = 1, ntri
        do j = 1, 4
          if ( j < 4 ) then
            k = nodtri(j,i)
          else
            k = nodtri(1,i)
          end if

          nval = nval + 1
          xval(nval) = rwork(2*k-1)
          yval(nval) = rwork(2*k)
          ival(nval) = 1
        end do
        ival(nval) = 0
      end do

    end if

  else if ( s_eqi(coord,'STAR') ) then

    l3d = .false.

    call star ( ierror, ival, val_max, nval, xval, yval )

  else if ( s_eqi(coord,'SXY') ) then

    l3d = .false.
    mdat = 2

    call getdat ( dat, idat, ierror, maxdat, mdat, 'xy', ndat )

    xdat(1:ndat) = dat(1:ndat,1)
    ydat(1:ndat) = dat(1:ndat,2)

    call sxy ( ival, val_max, ndat, nfine, nval, xdat, xval, ydat, yval )

  else if ( s_eqi(coord,'TXY') ) then

    l3d = .false.
    mdat = 3

    call getdat ( dat, idat, ierror, maxdat, mdat, 'txy', ndat )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) 'XYZ_PLOT - Fatal error!'
      go to 10
    end if

    xdat(1:ndat) = dat(1:ndat,2)
    ydat(1:ndat) = dat(1:ndat,3)

    nval = ndat

    xval(1:nval) = dat(1:nval,2)
    yval(1:nval) = dat(1:nval,3)
    ival(1:nval-1) = 1
    ival(nval) = 0

  else if ( s_eqi(coord,'TXYZ') ) then

    l3d = .true.
    mdat = 4

    call getdat ( dat, idat, ierror, maxdat, mdat, 'txyz', ndat )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) 'XYZ_PLOT - Fatal error!'
      go to 10
    end if

    x3dat(1:ndat) = dat(1:ndat,2)
    y3dat(1:ndat) = dat(1:ndat,3)
    z3dat(1:ndat) = dat(1:ndat,4)

    nval = ndat

    do i = 1, nval
      x3val(i) = dat(i,2)
      y3val(i) = dat(i,3)
      z3val(i) = dat(i,4)
    end do

    ival(1:nval-1) = 1
    ival(nval) = 0

  else if ( s_eqi(coord,'UVXY') ) then

    l3d = .false.
    mdat = 4

    call getdat ( dat, idat, ierror, maxdat, mdat, 'xyuv', ndat )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) 'XYZ_PLOT - Fatal error!'
      go to 10
    end if

    xdat(1:ndat) = dat(1:ndat,1)
    ydat(1:ndat) = dat(1:ndat,2)
    udat(1:ndat) = dat(1:ndat,3)
    vdat(1:ndat) = dat(1:ndat,4)

    do i = 1, ndat

      nval = nval + 1
      xval(nval) = xdat(i)
      yval(nval) = ydat(i)
      ival(nval) = 2

      nval = nval+1
      xval(nval) = xdat(i) + udat(i)
      yval(nval) = ydat(i) + vdat(i)
      ival(nval) = 0

    end do

  else if ( s_eqi(coord,'XYS') ) then

    l3d = .false.
    mdat = 2

    do

      call getdat ( dat, idat, ierror, maxdat, mdat, 'xy', ndat )

      if ( ndat > 3 ) then
        exit
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'You must define at least 4 data values!'

    end do

    xdat(1:ndat) = dat(1:ndat,1)
    ydat(1:ndat) = dat(1:ndat,2)

    call xys ( ival, maxdat, val_max, ndat, nfine, nval, tdat, xdat, &
      xval, ydat, yval )

  else if ( s_eqi(coord,'XYZ') ) then

    l3d = .true.
    mdat = 3

    call getdat ( dat, idat, ierror, maxdat, mdat, 'xyz', ndat )

    if ( ierror /= 0 ) then
      go to 10
    end if

    x3dat(1:ndat) = dat(1:ndat,1)
    y3dat(1:ndat) = dat(1:ndat,2)
    z3dat(1:ndat) = dat(1:ndat,3)

    nval = ndat
    do i = 1, nval
      x3val(i) = dat(i,1)
      y3val(i) = dat(i,2)
      z3val(i) = dat(i,3)
    end do
    ival(1:nval) = idat(1:nval)

  else if ( s_eqi(coord,'?') ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CAU     Caustic plot'
    write ( *, '(a)' ) 'CON     U(X,Y) contour plot.'
    write ( *, '(a)' ) 'FX      (X,Y) data.'
    write ( *, '(a)' ) 'FXY     (X,Y) defined by F(X,Y) = 0.'
    write ( *, '(a)' ) 'LXY     Least squares curve, (X,Y) data'
    write ( *, '(a)' ) 'LTXY    Least squares curve, (X(T), Y(T) ).'
    write ( *, '(a)' ) 'PWP     Piecewise polynomial, (X,Y) data'
    write ( *, '(a)' ) 'Q       Cancel this plot.'
    write ( *, '(a)' ) 'RT      ( R, THETA) polar plot.'
    write ( *, '(a)' ) 'SCAT    Scatter plot, (X,Y) data.'
    write ( *, '(a)' ) 'STAR    Plot a star.'
    write ( *, '(a)' ) 'SXY     Cubic spline, X, Y data.'
    write ( *, '(a)' ) 'TXY     ( X(T), Y(T) ).'
    write ( *, '(a)' ) 'TXYZ    ( X(T), Y(T), Z(T) ).'
    write ( *, '(a)' ) 'UVXY    2D vector plot, (U(X,Y), V(X,Y) ).'
    write ( *, '(a)' ) 'XYS     Plot ???'
    write ( *, '(a)' ) 'XYZ     (X, Y, Z) data.'
    go to 40

  else if ( s_eqi ( coord(1:1), 'Q' ) ) then
    go to 10
  else

    write ( *, '(a)' ) 'That was not a legal choice!'
    go to 10

  end if

  if ( ierror  /=  0 ) then
    write ( *, '(a)' ) 'Some error condition has occurred!'
    go to 10
  end if

!
!  Project 3D data.
!
50    continue

  if ( l3d ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Choose 3D -> 2D data projection method:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '-X     drop X coordinate, display YZ;'
    write ( *, '(a)' ) '-Y     drop Y coordinate, display XZ;'
    write ( *, '(a)' ) '-Z     drop Z coordinate, display XY;'
    write ( *, '(a)' ) 'OPLANE orthographic projection into plane.'
    write ( *, '(a)' ) 'PPLANE perspective projection into plane.'
    write ( *, '(a)' ) 'PX     project X into YZ using THETA;'
    write ( *, '(a)' ) 'PY     project Y into XZ using THETA;'
    write ( *, '(a)' ) 'PZ     project Z into XY using THETA;'

    read ( *, '(a)' ) isay

    if ( s_eqi ( isay(1:2), '-X' ) ) then

      xval(1:nval) = y3val(1:nval)
      yval(1:nval) = z3val(1:nval)

      xdat(1:ndat) = y3dat(1:ndat)
      ydat(1:ndat) = z3dat(1:ndat)

    else if ( s_eqi(isay,'-Y') ) then

      xval(1:nval) = x3val(1:nval)
      yval(1:nval) = z3val(1:nval)

      xdat(1:ndat) = x3dat(1:ndat)
      ydat(1:ndat) = z3dat(1:ndat)

    else if ( s_eqi(isay,'-Z') ) then

      xval(1:nval) = x3val(1:nval)
      yval(1:nval) = y3val(1:nval)

      xdat(1:ndat) = x3dat(1:ndat)
      ydat(1:ndat) = y3dat(1:ndat)

    else if ( s_eqi(isay,'OPLANE') ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Enter 3 points on the plane (X,Y,Z):'
      read ( *, * ) x1, y1, z1
      read ( *, * ) x2, y2, z2
      read ( *, * ) x3, y3, z3

      call proplane3 ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
        nval, x3val, y3val, z3val, x3val, y3val, z3val )

      call proplane2 ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
        nval, x3val, y3val, z3val, xval, yval )

      call proplane3 ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
        ndat, x3dat, y3dat, z3dat, x3dat, y3dat, z3dat )

      call proplane2 ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
        ndat, x3dat, y3dat, z3dat, xdat, ydat )

    else if ( s_eqi(isay,'PPLANE') ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Enter 3 points on the plane (X,Y,Z):'
      read ( *, * ) x1, y1, z1
      read ( *, * ) x2, y2, z2
      read ( *, * ) x3, y3, z3
      write ( *, '(a)') 'Enter focus point (X,Y,Z):'
      read ( *, * ) xf, yf, zf

      call plane_exp_project_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
        xf, yf, zf, nval, x3val, y3val, z3val, x3val, y3val, z3val, iwork )

      call proplane2 ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
        nval, x3val, y3val, z3val, xval, yval )

      call plane_exp_project_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
        xf, yf, zf, ndat, x3dat, y3dat, z3dat, x3dat, y3dat, z3dat, iwork )

      call proplane2 ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
        ndat, x3dat, y3dat, z3dat, xdat, ydat )

    else if ( s_eqi(isay,'PX') ) then

      call conv3d ( nval, theta, xval, yval, x3val, y3val, z3val )

      call conv3d ( ndat, theta, xdat, ydat, x3dat, y3dat, z3dat )

    else if ( s_eqi(isay,'PY') ) then

      call conv3d ( nval, theta, xval, yval, y3val, x3val, z3val )

      call conv3d ( ndat, theta, xdat, ydat, y3dat, x3dat, z3dat )

    else if ( s_eqi(isay,'PZ') ) then

      call conv3d ( nval, theta, xval, yval, z3val, x3val, y3val )

      call conv3d ( ndat, theta, xdat, ydat, z3dat, x3dat, y3dat )

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Unrecognized projection option!'
      go to 50

    end if

  end if
!
!  Draw the plot.
!
  call drwplt ( ierror, ival, l3d, ldat, lframe, ndat, nval, nxgrid, &
    nygrid, plhite, plwide, theta, title, vscale, xdat, xmax, xmaxw, xmin, &
    xminw, xval, ydat, ymax, ymaxw, ymin, yminw, yval )
!
!  Finish up
!
  icom = 9
  iplt1 = 1
  call anyplt ( icom )
!
!  Now that that plot's done, what next?
!
60    continue

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C  another curve on this plot'
  write ( *, '(a)' ) 'CHKDAT Print the data.'
  write ( *, '(a)' ) 'CHKVAL Print the values.'
  write ( *, '(a)' ) 'F  write plot data to file.'
  write ( *, '(a)' ) 'P  another plot'
  write ( *, '(a)' ) 'R  redraw'
  write ( *, '(a)' ) 'O  change option'
  write ( *, '(a)' ) 'Q  quit'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter option'

  read ( *, '(a)' ) isay

  if ( s_eqi ( isay, 'C' ) ) then

    go to 30

  else if ( s_eqi(isay,'F') ) then

    if ( l3d ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'What should the output file be called?'
      read ( *, '(a)' ) filexyz

      call xyz_write ( filexyz, ival, nval, x3val, y3val, z3val )

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'This option not ready for 2D plots.'

    end if

    go to 60

  else if ( s_eqi ( isay(1:1), 'O' ) ) then

    call option ( dev, iplot, ldat, lframe, ncon, nfine, &
      nxgrid, nygrid, theta, title, vscale )

    go to 60

  else if ( s_eqi ( isay(1:1), 'P' ) ) then

    go to 20

  else if ( s_eqi ( isay(1:1),'R' ) ) then

    go to 50

  else if ( s_eqi ( isay(1:1), 'Q' ) ) then

    icom = 1
    call anyplt ( icom )

    if ( .not. s_eqi(dev,'cgmb') ) then
      call file_delete ( 'cgmout' )
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYZ_PLOT:'
    write ( *, '(a)' ) '  Normal end of execution.'

    stop

  else if ( s_eqi(isay,'CHKVAL') ) then

    if ( .not. l3d ) then
      call chkval ( ival, l3d, nval, xval, yval, yval )
    else
      call chkval ( ival, l3d, nval, x3val, y3val, z3val )
    end if
    go to 60

  else if ( s_eqi(isay,'CHKDAT') ) then

    if ( .not. l3d ) then
      call chkdat ( idat, mdat, ndat, xdat, ydat, zdat )
    else
      call chkdat ( idat, mdat, ndat, x3dat, y3dat, z3dat )
    end if
    go to 60

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYZ_PLOT did not recognize that command!'
    go to 10

  end if

  stop
end
function angle_rad_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )

!*******************************************************************************
!
!! ANGLE_RAD_3D returns the angle in radians between two rays in 3D.
!
!  Discussion:
!
!    The routine always computes the SMALLER of the two angles between
!    two rays.  Thus, if the rays make an (exterior) angle of
!    1.5 radians, the (interior) angle of 0.5 radians will be reported.
!
!  Formula:
!
!    X dot Y = Norm(X) * Norm(Y) * Cos ( Angle(X,Y) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are three points
!    which define the rays.  The rays are:
!    ( X1-X2, Y1-Y2, Z1-Z2 ) and ( X3-X2, Y3-Y2, Z3-Z2 ).
!
!    Output, real ANGLE_RAD_3D, the angle between the two rays, in radians.
!    This value will always be between 0 and PI.  If either ray has
!    zero length, then the angle is returned as zero.
!
  implicit none

  real angle_rad_3d
  real dot
  real dot0_3d
  real enorm0_3d
  real v1norm
  real v2norm
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3

  dot = dot0_3d ( x2, y2, y2, x1, y1, z1, x3, y3, z3 )
  v1norm = enorm0_3d ( x1, y1, z1, x2, y2, z2 )
  v2norm = enorm0_3d ( x3, y3, z3, x2, y2, z2 )

  if ( v1norm == 0.0E+00 .or. v2norm == 0.0E+00 ) then
    angle_rad_3d = 0.0E+00
  else
    angle_rad_3d = acos ( dot / ( v1norm * v2norm ) )
  end if

  return
end
subroutine arrow ( xstart, ystart, xtip, ytip, ndraw, xdraw, ydraw )

!***********************************************************************
!
!! ARROW returns points that specify an arrow from one point to another.
!
!  Discussion:
!
!    The arrow will stretch between two user specified points.
!
!    The "head" of the arrow may be fatter or thinner than expected
!    if the X and Y scales of the graph are not in the same
!    proportions.
!
!
!                       left(3)
!                        |\
!                        | \
!                        |  \
!    start(1)*-----base(2,6) * tip(4)
!                        |  /
!                        | /
!                        |/
!                       rite(5)
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XSTART, YSTART, the starting point for the arrow.
!
!    Input, real XTIP, YTIP, the end point for the arrow.
!
!    Output, integer NDRAW, the number of points to draw, which will be 6.
!
!    Output, real XDRAW(NDRAW), YDRAW(NDRAW), the X and Y coordinates
!    of the points to connect to draw the arrow.
!
  implicit none

  real alpha
  real del
  real dist
  integer ndraw
  real, parameter :: pi = 3.14159265358979E+00
  real theta
  real xbase
  real xdraw(*)
  real xleft
  real xrite
  real xstart
  real xtip
  real ybase
  real ydraw(*)
  real yleft
  real yrite
  real ystart
  real ytip

  theta = 0.5E+00 * pi - atan2 ( 2.0E+00, 1.0E+00 )
  dist = sqrt ( ( xtip - xstart )**2 + ( ytip - ystart )**2 )
  alpha = atan2 ( ytip - ystart, xtip - xstart )
  del = sqrt ( 5.0E+00 ) / 3.0E+00

  xbase = ( xstart + 2.0E+00 * xtip ) / 3.0E+00
  ybase = ( ystart + 2.0E+00 * ytip ) / 3.0E+00

  xleft = xstart + del * dist * cos ( alpha - theta )
  yleft = ystart + del * dist * sin ( alpha - theta )

  xrite = xstart + del * dist * cos ( alpha + theta )
  yrite = ystart + del * dist * sin ( alpha + theta )

  ndraw = 6
  xdraw(1:6) = (/ xstart, xbase, xleft, xtip, xrite, xbase /)
  ydraw(1:6) = (/ ystart, ybase, yleft, ytip, yrite, ybase /)

  return
end
subroutine ch_cap ( c )

!*******************************************************************************
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*******************************************************************************
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
function ch_is_alpha ( c )

!*******************************************************************************
!
!! CH_IS_ALPHA returns TRUE if C is an alphabetic character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, a character to check.
!
!    Output, logical CH_IS_ALPHA is TRUE if C is an alphabetic character.
!
  implicit none

  character c
  logical ch_is_alpha

  if ( ( lle ( 'a', c ) .and. lle ( c, 'z' ) ) .or. &
       ( lle ( 'A', c ) .and. lle ( c, 'Z' ) ) ) then
    ch_is_alpha = .true.
  else
    ch_is_alpha = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*******************************************************************************
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine cau ( ival, val_max, nval, xval, yval )

!***********************************************************************
!
!! CAU draws a caustic plot.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer IVAL(VAL_MAX).  IVAL(I) is 1 if point I
!    is to be connected to point I+1.
!
!    Input, integer VAL_MAX, the maximum number of points that can be 
!    used for plots.
!
!    Input/output, integer NVAL, the current number of points.
!
!    Input/output, real XVAL(VAL_MAX), YVAL(VAL_MAX), the coordinates
!    of points to be used in graphs.
!
  implicit none

  integer val_max

  real angle
  integer i
  integer ios
  integer ip
  integer iq
  integer ival(val_max)
  integer j
  integer nval
  real, parameter :: pi = 3.14159265358979E+00
  real xval(val_max)
  real yval(val_max)

  do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the number of points on the circle.'
    read ( *, *, iostat = ios ) iq

    if ( ios /= 0 ) then
      return
    end if

    if ( iq < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Please choose a value that is more than 3!'
    else if ( nval+2*iq > val_max ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'Choose a value no more than ', val_max
    else
      exit
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter point to which point 1 is connected.'
  read ( *, *, iostat = ios ) ip

  if ( ios /= 0 ) then
    return
  end if
!
!  Draw lines between points:
!    1 and   P
!    2 and 2*P
!    ...
!    Q and Q*P
!
  do i = 1, iq

    nval = nval + 1
    angle = ( 2.0E+00 * pi * i ) / real ( iq )
    xval(nval) = cos ( angle )
    yval(nval) = sin ( angle )
    ival(nval) = 1

    j = mod ( ip * i, iq )

    nval = nval + 1
    angle = ( 2.0E+00 * pi * j ) / real ( iq )
    xval(nval) = cos ( angle )
    yval(nval) = sin ( angle )
    ival(nval) = 0

  end do

  return
end
subroutine chkdat ( idat, mdat, ndat, xdat, ydat, zdat )

!***********************************************************************
!
!! CHKDAT prints out the value of the data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ndat

  integer i
  integer idat(ndat)
  integer mdat
  real xdat(ndat)
  real ydat(ndat)
  real zdat(ndat)

  if ( ndat  <=  0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ChkDat: No data is defined!'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ChkDat: User data values:'
  write ( *, '(a)' ) ' '

  do i = 1, ndat

    if ( mdat == 2 ) then
      write ( *, '(i5,2g15.6)' ) i, xdat(i), ydat(i)
    else
      write ( *, '(i5,3g15.6)' ) i, xdat(i), ydat(i), zdat(i)
    end if

    if ( idat(i) == 0 ) then
      write ( *, '(a)' ) ' '
    end if

  end do

  return
end
subroutine chkval ( ival, l3d, nval, xval, yval, zval )

!***********************************************************************
!
!! CHKVAL prints out the value of the plot points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL(VAL_MAX).  IVAL(I) is 1 if point I
!    is to be connected to point I+1.
!
!    Input/output, integer NVAL, the current number of plot points.
!
  implicit none

  integer nval

  integer i
  integer ival(nval)
  logical l3d
  real xval(nval)
  real yval(nval)
  real zval(nval)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ChkVal'
  write ( *, '(a)' ) '  Plot points:'
  write ( *, '(a)' ) ' '

  do i = 1, nval
    if ( .not. l3d ) then
      write(*,'(1x,i5,2g15.6)') i, xval(i), yval(i)
    else
      write(*,'(1x,i5,3g15.6)') i, xval(i), yval(i), zval(i)
    end if
    if ( ival(i) == 0 ) then
      write ( *, '(a)' ) ' '
    end if
  end do

  return
end
subroutine comrpn ( com, ierror, ifrm, infix, irpn, maxrpn, namvar, value )

!*******************************************************************************
!
!! COMRPN can translate formulas you type in and evaluate them.
!
!  This means that you can design interactive FORTRAN programs
!  which can input their formulas at run time, rather than
!  recompiling every time you want a new formula.
!
!  If you declare vectors or matrices, you may reference entries
!  such as X(3) or A(3,2) or even X(2+3).
!
!
!  Formulas are made from operators, constants, punctuation,
!  variables, and  functions.
!
!
!  The list of legal operators includes:
!
!    +  -  *  /  **  ^  = \
!
!  * means multiplication, and is standard matrix multiplication if
!  both arguments are (conformable) matrices.  If one argument is a
!  scalar, it multiplies all the entries of the other argument.
!  If both arguments are row or column vectors, the dot product
!  is taken.
!
!  Thus,
!
!    x is a scalar equal to 2,
!    y is a scalar equal to 3,
!    u is a vector equal to (1, 2, 3),
!    v is a vector equal to (1, 0, 2),
!
!  then
!
!    2*3 results in 6.
!
!    x*y results in 6.
!
!    u*v results in 7.
!
!    x*v would result in the vector (2, 0, 4).
!
!
!  / means division, as in A/B, but for matrices the only allowed
!  form requires that B be a scalar, in which case each element of
!  A is divided by B.
!
!
!  + and - are allowed for pairs of scalars, vectors or matrices of
!  the same order, and also for a square matrix and a scalar,
!  in which case A+2 is interpreted as A+2*I.
!
!
!  ** or ^ means exponentiation.  In the scalar case, any
!  nonnegative scalar can be taken to any power, positive, zero, or
!  negative.  A negative scalar may only be taken to an integer
!  power.
!
!  ** is also supported for square matrices, which can be taken
!  to any nonnegative integer power.  Nonsingular matrices may
!  also be taken to negative integer powers.
!
!
!  = is used to assign values.  For vectors or matrices, this may
!  also be used to assign a single entry, as in A(3,2)=7.
!
!
!  / is standard scalar division.  It is used in a formula like
!
!    X/Y
!
!  or
!
!    (A*B)/(X+1)
!
!  Normally, X would be a scalar quantity.  However, as NONSTANDARD
!  shorthand, you may 'divide' a matrix by a scalar, in which case
!  each entry of the matrix is divided by the scalar, e.g.
!
!    A/2 = (1/2)*A.
!
!  MATALG will NOT allow you to use the "/" operator to represent
!  multiplication by an inverse matrix.  Thus, if A is a matrix,
!  the statement
!
!    B/A
!
!  is illegal.  However, the statement
!
!    INV(A)*B
!
!  will compute the inverse of A and multiply by B, and the
!  statement
!
!    A \ B
!
!  will compute the LU factorization of A, and use that to
!  compute the inverse of A times B.
!
!
!  Constants:
!
!
!  Real and integer constants may appear in your formulas, as well
!  as the symbol 'PI' and the machine unit roundoff number 'EPS'
!  which is the power of two such that 1+EPS>1 but 1+(EPS/2)=1.
!
!
!  Punctuation:
!
!
!    Blanks may be used anywhere.
!
!    Commas are used to separate arguments in a function, such as
!      MAX(X,Y).
!
!    Parentheses are used:
!
!      to group quantities:  (X+Y)*Z
!
!      to reference an entry of a vector or matrix: A(2,2)
!
!      to mark the argument of a function: SIN(X), MAX(X,Y)
!
!
!  Functions and operators:
!
!
!  S, S1, S2: Arguments which may only be scalar.
!  V        : Arguments may only be a scalar or vector.
!  M        : Arguments may only be scalar or square matrix.
!  MV       : Arguments may only be scalar, vector, or square matrix.
!  *        : Arguments which may be scalar, vector, or matrix.
!  I        : Arguments which may only be positive integers.
!
!  ABS(*)         Absolute value of *.
!
!  ACOS(S)        The arc cosine of S.
!                 -1 <= S <=1.
!
!  ASIN(S)        The arc sine of S.
!                 -1 <= S <=1.
!
!  ATAN(S)        The arc tangent of S.
!
!  ATAN2(S1,S2)   Arc tangent of (S1/S2)
!                 Correctly computes ATAN2(1.0,0.0)=PI/2.
!
!  COS(S)         The cosine of S, with S measured in radians.
!
!  COSD(S)        The cosine of S, with S measured in degrees.
!
!  COSH(S)        Hyperbolic cosine of S.
!
!  DET(M)         Determinant of square matrix M.
!
!  DIAG(*)        The diagonal entries of *, stored in a vector.
!
!  E              The base of the natural logarithm system.
!                 E=2.71828...
!                 You may not change the value of E.
!
!  EPS            The machine epsilon, or unit roundoff number.
!                 EPS is the power of 2 such that
!
!                   1 < 1+EPS and 1=1+(EPS/2).
!
!                 You may not change the value of EPS.
!
!  EVAL(M)        Real and imaginary parts of eigenvalues of matrix
!                 M, stored in a matrix of N rows and 2 columns,
!                 with the real parts in column 1, and the
!                 imaginary parts in column 2.
!
!  EVEC(M)        Eigenvectors of square matrix M, stored in a
!                 square matrix of same size as M.
!
!                 (Eigenvectors corresponding to a complex pair:
!                 If eigenvalues j and j+1 are a complex pair,
!                 then the eigenvector for eigenvalue j is
!                 column j + i*column j+1, and the eigenvector
!                 for eigenvalue j+1 is column j - i*column j+1.)
!
!  EXP(S)         Exponential of S.
!
!  GCF(I,J)       Greatest common factor of two integers.
!
!  HILBERT(I)     The Hilbert matrix of order I.
!
!  HILBINV(I)     The inverse of the Hilbert matrix of order I.
!
!  HOUSE(V)       The Householder elementary reflector matrix H,
!                 defined as
!
!                   H = I - 2*(V * TRANSPOSE(V)) / (TRANS(V)*V)
!
!  ID(I)          The matrix identity of order I.
!                 ID(3) is the 3 by 3 identity, for instance.
!
!  INT(*)         Truncates real values to their integer part.
!                 INT(1.0) = INT(1.1) = INT(1.9) = 1.
!
!  INV(M)         The inverse matrix of the square matrix M.
!                 If M is singular, there will be no inverse.
!
!  IVAL(M)        Imaginary parts of eigenvalues of square matrix M
!                 stored in a column vector.
!
!  LCM(I,J)       Least common multiple of two integers.
!
!  LN(S)          The natural logarithm of S.
!                 S must be greater than 0.
!
!  LOG(S)         The natural logarithm of S.
!                 S must be greater than 0.
!
!  LOG10(S)       The logarithm base 10 of S.
!                 S must be greater than 0.
!
!  LOG2(S)        The logarithm base 2 of S.
!                 S must be greater than 0.
!
!  MAX(S1,S2)     The maximum of S1 and S2.
!
!  MIN(S1,S2)     Minimum of S1 and S2.
!
!  NEG(*)         Changes sign of *.
!
!  NINT(*)        Nearest integer value to entries of *.
!
!  NORM0(*)       Maximum or infinity norm of a vector or matrix.
!                 NORM0(S) = ABS(S)
!                 NORM0(V) returns the maximum of the absolute
!                   values of the entries of V.
!                 NORM0(M) returns the maximum of the sum of the
!                   absolute values of the entries in each row.
!
!  NORM1(*)       The L1-norm of a vector or matrix.
!                 NORM1(S) = ABS(S).
!                 NORM1(V) returns the sum of the absolute values
!                   of the entries of V.
!                 NORM0(M) returns the maximum of the sum of the
!                   absolute values of the entries in each column.
!
!  NORM2(MV)      L2-norm, Euclidean norm or root-mean-square norm
!                 of a vector or a square matrix M.  NORM2(V) returns
!                 the square root of the sum of the squares of the
!                 entries of V.  NORM2(M) returns the maximum magnitude
!                 of the eigenvalues of Transpose(M)*M.
!
!  NORMF(*)       The Frobenius norm.  NORMF(M) returns the square
!                 root of the sum of the squares of all the entries
!                 of a matrix M.  NORMF may also be applied to a
!                 vector, giving the same results as NORM2(V).
!
!  NORMS(MV)      The spectral "norm".  NORMS(M) returns the
!                 value of the maximum of the absolute values of
!                 all the eigenvalues of the square matrix M.
!
!                 Note that the spectral norm is NOT a true
!                 norm, although it is frequently used as though
!                 it were.
!
!  PI             3.14159265358979...
!                 You may not change the value of PI.
!
!  POLY(V,M)      Polynomial evaluation.  V contains the
!                 coefficients of the polynomial, with V(1) the
!                 coefficient of M**(N-1) and V(N) the constant
!                 term.
!
!                 M is the scalar or square matrix argument of the
!                 polynomial.
!
!  RAN(*)         Fills * with random numbers between 0 and 1.
!
!  ROUND(*)       Replaces all "small" entries of * by 0.
!                 Here "small" means less than RTOL in absolute
!                 value.
!
!  RTOL           The tolerance for the ROUND command, initially set
!                 to EPS.  To change the value of RTOL, simply
!                 reassign it with a statement like "RTOL=0.1E-6".
!
!  RVAL(M)        Real parts of eigenvalues of square matrix M,
!                 stored in a column vector.
!
!  SIN(S)         The sine of S, with S measured in radians.
!
!  SIND(S)        The sine of S, with S measured in degrees.
!
!  SINH(S)        Hyperbolic sine of S.
!
!  SQRT(S)        Square root of S.  S must be non-negative.
!
!  STEP(S)        Heavyside step function.
!
!                   STEP=0 if S < 0,
!                   STEP=1 if 0 <= S.
!
!  TAN(S)         Tangent of S.
!
!  TAND(S)        Tangent of S measured in degrees.
!
!  TANH(S)        Hyperbolic tangent of S.
!
!  TRACE(*)       Sum of diagonal elements of *.
!
!  TRANS(*)       Transpose of matrix or vector.
!
!                 A=TRANS(A) is a legal formula if A is square.
!
!  ZEROS(I)       Matrix zero of order I.
!                 ZEROS(4) is the 4 by 4 zero matrix.
!
!  I!             Factorial of I.  I should be an integer from 0 to
!                 25.  0!=1, 1!=1, 2!=1*2=2, 3!=1*2*3=6, and so on.
!                 For large I, the exact value is not returned.
!
!*******************************************************************************
!
!  To add a new function to COMRPN:
!
!    Update the information in the IOPSYM, IPRSYM, and SYMBOL
!      arrays.
!
!    Insert text into FUNVAL or FUNSCL which evaluates the
!      function, and determines its dimensions.
!
!*******************************************************************************
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Parameters:
!
!  COM    Input, character COM.
!
!         The command to be executed, which determines what other
!         input must be supplied.
!
!
!         COM='A'  Add (=declare) a variable:
!
!         Causes the program to add the variable named in NAMVAR
!         to its list of symbols.  From now on, it may be used in
!         formulas, and values may be assigned to it.  If you
!         declare a name which already exists, then if that name
!         was already declared by the program, your command is
!         ignored.  Otherwise, your new definition overrides your
!         old one.
!
!         Input for the "A" command includes NAMVAR, NROW,
!         and NCOL.
!
!
!         COM='D'  delete a variable:
!
!         This command can be used to free up memory.  You
!         must also give the name of the variable to be deleted.
!         The variable must be one you have already defined.
!
!         Input for the "D" command includes NAMVAR.
!
!
!         COM='E'  evaluate a formula:
!
!         Using current values of variables, evaluate the formula.
!         The formula number is given by IFRM.  The value of the
!         formula is returned in VALUE, and its dimensions in NROW
!         and NCOL.
!
!         Input to the "E" command includes IFRM.
!         Output from the "E" command includes VALUE.
!
!
!         COM='F', Formula is to be compiled:
!
!         The user assigns an index, IFRM, to the formula, so that
!         it can be referred to later.  The formula is stored in
!         INFIX on input.  The formula may refer to variables which
!         have not yet been declared.
!
!         Output from the "F" command includes IFRM.
!
!
!         COM='G', Formula is to be compiled, all variables
!         in the formula have already been declared:
!
!         The same as 'F' except that the formula may not refer to
!         any variables which have not been declared.
!
!         Output from the "G" command includes IFRM.
!
!
!         COM='I'  Initialize:
!
!         This must be the first call, to initialize COMRPN's internal
!         data.  Also, if old data is to be flushed out, this command
!         can be used.
!
!
!         COM='R'  Read value of symbol:
!
!         The program returns the current value of the variable
!         NAMVAR.  VALUE will contain the value, and NROW and NCOL
!         the dimensions.
!
!
!         COM='V'  assign value to symbol:
!
!         This command assigns the value of the numbers in VALUE
!         to the variable NAMVAR, with NROW and NCOL specifying
!         the dimension of the variable.
!
!  IERROR Output, integer IERROR.
!
!         IERROR is an error flag.
!
!         0 Means no error occurred.
!
!         1 means some kind of error.  These errors are usually
!         signaled by a printed message.  They can include the
!         following problems:
!
!         On adding a variable, if the name was already in use,
!         either by you or the program, the program returns
!         with this error warning.
!
!         On setting a variable, if the name supplied does
!         not represent any recognized variable, the program
!         ignores the command and returns.
!
!         On compiling a formula, the program will return with
!         an error warning if any of the following occur:
!
!         The formula does not compile.  An unknown variable,
!         missing parentheses, garbled information or mistyping
!         may be responsible.  For that matter, the compiler
!         may fail to compile perfectly good, but complicated
!         formulas.
!
!         There is not enough space in IRPN to store the RPN code.
!         In this case, you must either re-initialize, or increase
!         the storage available in IRPN and signal this with
!         an increased value of MAXRPN.
!
!         On evaluating a formula, an error can occur if there
!         is no formula corresponding to the given index IFRM.
!
!         Certain arithmetic errors, such as division
!         by zero, may be caught and flagged by the compiler.
!
!  IFRM   Input/output, integer IFRM.
!
!         COMRPN can store more than one formula internally, at
!         one time.  COMRPN refers to the formulas by an index
!         number.  When a formula is entered, with the "F" or "G"
!         command, it is assigned an index number, whose value is
!         returned to the user.  Then, whenever the formula is to
!         be evaluated, the user must input the value of IFRM, so
!         that COMRPN knows which formula to evaluate.  In many
!         cases, only one formula is being used, so that IFRM is
!         usually just 1.
!
!         For the "E" command, IFRM is the index of the formula to
!         evaluate.
!
!         For the "F" or "G" commands, the index assigned to the
!         formula.
!
!  INFIX  character ( len = * ) INFIX, for the F or G commands,
!         contains the user's formula to be compiled and evaluated.
!
!  IRPN   integer IRPN(MAXRPN), for the F, G or E commands,
!         a vector used to store the compiled versions of the input
!         formulas.
!
!  MAXRPN Input, integer MAXRPN.
!
!         The maximum length of IRPN.  80 should be enough.
!
!  NAMVAR Input, character ( len = * ) NAMVAR, used for the A, D, L, R and
!         V commands, containing the name of the variable to be
!         added, deleted, listed, read or valued.
!
!         NAMVAR should not be a blank, nor should it be the name
!         'VOID'.
!
!  VALUE  real VALUE, used to pass values in and out
!         of the program. 
!
!         For the V command, VALUE is the input value to be
!         assigned to a variable.
!
!         For the E command, VALUE is the output value of the formula.
!
!         For the R command, VALUE is the output value of the variable.
!
!*******************************************************************************
!
!  Internal variables:
!
!  IOPSYM Internal, integer IOPSYM(MAXSYM), contains, for each
!         symbol, the number of operands.  In particular, if
!         a symbol represents a constant, IOPSYM(I) is 0.
!         If a symbol represents a unary operator such as ABS,
!         IOPSYM(I) is 1.  IOPSYM may be as large as 3.
!
!  IPRSYM Internal, integer IPRSYM(MAXSYM), contains, for each
!         symbol, the operator "priority" of that symbol.  This
!         is simply to help COMRPN deal with ambiguous formulas like
!         A*B+C.
!
!  IPVAL  Internal, integer IPVAL(MAXSYM), contains, for each
!         symbol, the address in VALSYM where the value of the
!         symbol is stored, if it is a scalar.  If the symbol
!         represents a vector or matrix, IPVAL(IARG) points to
!         the location of the first entry.
!
!  MAXSYM Internal, integer MAXSYM, controls the maximum number
!         of symbols allowed.  This includes permanent symbols
!         defined by the program, as well as symbols the user
!         declares while running the program.
!
!*******************************************************************************
!
  implicit none

  integer, parameter :: maxsym = 150
  integer, parameter :: maxval = 2500
  integer, parameter :: maxlen = 20

  integer maxrpn

  character com
  integer ierror
  integer, save :: ifinis
  integer, save :: ifree
  integer ifreesv
  integer ifrm
  integer ihi
  integer ilo
  integer imin
  integer implic
  integer, save :: indx1
  integer, save :: indx2
  integer, save :: ineg
  character ( len = * ) infix
  integer, dimension ( 80 ), save :: intsym
  integer, dimension ( maxsym ), save :: iopsym
  integer, dimension ( maxsym ), save :: iprsym
  integer, dimension ( maxsym ), save :: ipval
  integer, save :: irad
  integer irpn(maxrpn)
  integer, save :: irtol
  integer, dimension ( maxsym ), save :: istack
  integer maxfix
  integer maxfrm
  character ( len = * ) namvar
  character ( len = maxlen ) namvr
  integer, save :: nints
  integer nrpn
  integer, save :: nsym
  integer, save :: nsymp
  integer, save :: nsyms
  integer nuvar
  logical s_eqi
  logical s_paren_check
  character ( len = maxlen ), dimension ( maxsym ), save :: symbol
  real value
  real, dimension ( maxval ), save :: valsym
!
  ierror = 0
  implic = 0
  namvr = namvar
  maxfix = len ( infix )
  if ( maxfix > 0 ) then
    maxfrm = ( maxrpn / maxfix )
  else
    maxfrm = 0
  end if
  call s_blank_delete(namvr)

  if ( s_eqi ( com, 'E' ) .or. s_eqi ( com, 'G' ) ) then

    if ( ifrm <= 0 .or. ifrm > maxfrm ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMRPN - Error!'
      write ( *, '(a,i6)' ) '  Illegal formula index= ', ifrm
      ierror = 1
      return
    end if

  end if
!
!  COM=I initialize.
!
  if ( s_eqi ( com, 'I' ) ) then

    call inicom ( ifinis, ifree, indx1, indx2, ineg, infix, intsym, iopsym, &
      iprsym, ipval, irad, irpn, irtol, istack, maxrpn, maxsym, &
      maxval, nints, nsym, nsymp, nsyms, symbol, valsym, value )

    return
!
!  COM=A, Add variable.
!
  else if (s_eqi(com,'a') ) then

    call symadd ( ierror, ifree, iopsym, iprsym, ipval, &
      maxsym, maxval, namvr, nsym, nsymp, symbol, valsym )

    if ( ierror == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMRPN - Note:'
      write ( *, '(a)' ) '  Added the variable ' // trim ( namvr )
    end if

    return
!
!  COM=V, Set variable value.
!  COM=R, Get variable value.
!
  else if ( s_eqi ( com, 'r' ) .or. s_eqi ( com, 'v' ) ) then

    call symbol_value ( com, ierror, ifree, iopsym, iprsym, ipval, maxsym, &
      maxval, namvr, nsym, nsymp, symbol, valsym, value )

    return
!
!  COM=F, compile formula.
!
  else if ( s_eqi ( com, 'f' ) ) then

    nuvar=1
!
!  Remove all blanks from the formula.
!
    call s_blank_delete ( infix )
!
!  Reject the formula if it has no information in it.
!
    if ( len_trim ( infix ) <= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMRPN - Fatal error!'
      write ( *, '(a)' ) '  The formula has zero length.'
      return
    end if
!
!  Do a simple check on parentheses.
!
    if ( .not. s_paren_check ( infix ) ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMRPN - Fatal error!'
      write ( *, '(a)' ) '  The formula did not pass the parentheses checks!'
      return
    end if
!
!  Convert INFIX, the string of characters, into INTSYM, an array of
!  integers which are indices of tokens.
!
    call tokens ( ierror, ifinis, ifree, implic, indx1, indx2, ineg, infix, &
      intsym, iopsym, iprsym, ipval, maxsym, maxval, nints, nsym, nsymp, &
      nuvar, symbol, valsym )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMRPN - Fatal error!'
      write ( *, '(a)' ) '  Could not convert formula into tokens.'
      if ( implic /= 0 ) then
        write ( *, '(a)' ) '  Cancelling the variable ' // symbol(implic)
        symbol(implic) = 'VOID'
      end if
      return
    end if
!
!  Convert INTSYM, which is an infix formula, into IRPN, which
!  is an RPN formula.
!
    imin = (ifrm-1) * 80 + 1

    call rpnset ( ierror, ifinis, imin, intsym, iopsym, iprsym, irpn, &
      istack, maxrpn, maxsym, nints, nrpn, symbol )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMRPN - Fatal error!'
      write ( *, '(a)' ) '  Error!  Could not compile formula.'
      if ( implic /= 0 ) then
        write ( *, '(a)' ) '  Cancelling the variable ' // symbol(implic)
        symbol(implic) = 'VOID'
      end if
      return
    end if
!
!  Check that operators and arguments are balanced.
!
    ihi = imin + nrpn
    call rpnchk ( ierror, ihi, ilo, iopsym, irpn, maxrpn, maxsym )

    if ( ierror /= 0 .or. ilo /= imin + 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMRPN - Error!'
      write ( *, '(a)' ) '  Parentheses mismatch.'
      write ( *, '(a,i6)' ) '  ILO = ', ilo
      write ( *, '(a,i6)' ) '  IMIN = ', imin
      ierror = 1
      if ( implic /= 0 ) then
        write ( *, '(a)' ) '  Cancelling the variable ' // symbol(implic)
        symbol(implic) = 'VOID'
      end if
      return
    end if
!
!  For implicit definition via equality, evaluate formula to
!  get dimensions.
!
    if ( implic /= 0 ) then
      go to 10
    end if

  else if ( s_eqi ( com, 'E' ) ) then

    go to 10

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMRPN - Fatal error!'
    write ( *, '(a)' ) 'Error!  Unknown command = ' // trim ( com )
    ierror = 1

  end if

  return
!
!  COM=E, Evaluate formula.
!
10    continue

  imin = (ifrm-1) * 80 + 1
  nrpn = 80
  ifreesv = ifree

  call rpnval ( ierror, ifree, imin, iopsym, iprsym, ipval, irad, irpn, &
    irtol, istack, maxrpn, maxsym, maxval, nrpn, nsym, &
    nsyms, symbol, valsym, value )

  ifree = ifreesv

  if ( ierror /= 0 ) then
    value = 0.0E+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMRPN - Fatal error!'
    write ( *, '(a)' ) '  The formula could not be evaluated!'
    if ( implic /= 0 )  then
      write ( *, '(a)' ) '  Cancelling the variable ' // symbol(implic)
      symbol(implic) = 'VOID'
    end if
    return
  end if

  return
end
subroutine conv3d ( nval, theta, xval, yval, x3val, y3val, z3val )

!***********************************************************************
!
!! CONV3D converts 3D data to a 2D projection.
!
!  Discussion:
!
!    A "presentation angle" THETA is used to project the 3D point
!    (X3D, Y3D, Z3D) to the 2D projection (XVAL,YVAL).
!
!    The formula used is
!
!      X2D = - sin ( theta ) * X3D + Y3D
!      Y2D = - sin ( theta ) * X3D + Z3D
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NVAL, the number of 3D values to be converted.
!
!    Input, real THETA, the presentation angle.
!
!    Output, real XVAL(NVAL), YVAL(NVAL), the 2D projections of the
!    input data.
!
!    Input, real X3VAL(NVAL), Y3VAL(NVAL), Z3VAL(NVAL), the 3D points
!    to be projected.
!
  implicit none

  integer nval

  real angle
  real, parameter :: pi = 3.14159265358979E+00
  real s
  real theta
  real x3val(nval)
  real xval(nval)
  real y3val(nval)
  real yval(nval)
  real z3val(nval)

  angle = pi * theta / 180.0E+00
  s = sin ( angle )

  xval(1:nval) = - s * x3val(1:nval) + y3val(1:nval)
  yval(1:nval) = - s * x3val(1:nval) + z3val(1:nval)

  return
end
subroutine data_to_dif ( diftab, ntab, xtab, ytab )
!
!*******************************************************************************
!
!! DATA_TO_DIF sets up a divided difference table from raw data.
!
!  Discussion:
!
!    Space can be saved by using a single array for both the DIFTAB and
!    YTAB dummy parameters.  In that case, the divided difference table will
!    overwrite the Y data without interfering with the computation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real DIFTAB(NTAB), the divided difference coefficients
!    corresponding to the input (XTAB,YTAB).
!
!    Input, integer NTAB, the number of pairs of points
!    (XTAB(I),YTAB(I)) which are to be used as data.  The
!    number of entries to be used in DIFTAB, XTAB and YTAB.
!
!    Input, real XTAB(NTAB), the X values at which data was taken.
!    These values must be distinct.
!
!    Input, real YTAB(NTAB), the corresponding Y values.
!
  implicit none

  integer ntab

  real diftab(ntab)
  integer i
  integer j
  logical rvec_distinct
  real xtab(ntab)
  real ytab(ntab)

  if ( .not. rvec_distinct ( ntab, xtab ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_TO_DIF - Fatal error!'
    write ( *, '(a)' ) '  Two entries of XTAB are equal!'
    return
  end if
!
!  Copy the data values into DIFTAB.
!
  diftab(1:ntab) = ytab(1:ntab)
!
!  Compute the divided differences.
!
  do i = 2, ntab
    do j = ntab, i, -1

      diftab(j) = ( diftab(j) - diftab(j-1) ) / ( xtab(j) - xtab(j+1-i) )

    end do
  end do

  return
end
subroutine device_choose ( dev )

!***********************************************************************
!
!! DEVICE_CHOOSE has the user choose the output device.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  character ( len = * ) dev
  logical s_eqi

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Please choose a graphics device:'
  write ( *, '(a)' ) '  Options are CGM, PS or XWS.'
  read ( *, '(a)' ) dev

  call s_blank_delete ( dev )

  if ( s_eqi ( dev(1:3), 'cgm' ) ) then
    dev = 'cgmb'
    write ( *, '(a)' ) 'Output to a CGM binary file "anyplt.cgm".'
  else if ( s_eqi ( dev, 'ps' ) ) then
    write ( *, '(a)' ) 'Output to a PostScript file "anyplt.ps".'
  else if ( s_eqi ( dev, 'xws' ) ) then
    write ( *, '(a)' ) 'Output to an X window screen.'
  else
    write ( *, '(a)' ) 'Your device ' // trim ( dev ) // ' was not recognized!'
    dev = 'xws'
  end if

  return
end
function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )
!
!***********************************************************************
!
!! DIAEDG chooses a diagonal edge.
!
!  Discussion:
!
!    The routine determines whether 0--2 or 1--3 is the diagonal edge 
!    that should be chosen, based on the circumcircle criterion, where 
!    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple 
!    quadrilateral in counterclockwise order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!    coordinates of the vertices of a quadrilateral, given in
!    counter clockwise order.
!
!    Output, integer DIAEDG, chooses a diagonal:
!    +1, if diagonal edge 02 is chosen;
!    -1, if diagonal edge 13 is chosen;
!     0, if the four vertices are cocircular.
!
  implicit none

  real ca
  real cb
  integer diaedg
  real dx10
  real dx12
  real dx30
  real dx32
  real dy10
  real dy12
  real dy30
  real dy32
  real s
  real tol
  real tola
  real tolb
  real x0
  real x1
  real x2
  real x3
  real y0
  real y1
  real y2
  real y3

  tol = 100.0E+00 * epsilon ( tol )

  dx10 = x1 - x0
  dy10 = y1 - y0
  dx12 = x1 - x2
  dy12 = y1 - y2
  dx30 = x3 - x0
  dy30 = y3 - y0
  dx32 = x3 - x2
  dy32 = y3 - y2

  tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
  tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

  ca = dx10 * dx30 + dy10 * dy30
  cb = dx12 * dx32 + dy12 * dy32

  if ( ca > tola .and. cb > tolb ) then

    diaedg = -1

  else if ( ca < -tola .and. cb < -tolb ) then

    diaedg = 1

  else

    tola = max ( tola, tolb )
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( s > tola ) then
      diaedg = -1
    else if ( s < -tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

  return
end
subroutine dif_val ( diftab, ntab, xtab, xval, yval )

!*******************************************************************************
!
!! DIF_VAL evaluates a divided difference polynomial at a point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real DIFTAB(NTAB), the divided difference polynomial coefficients.
!
!    Input, integer NTAB, the number of divided difference
!    coefficients in DIFTAB, and the number of points XTAB.
!
!    Input, real XTAB(NTAB), the X values upon which the
!    divided difference polynomial is based.
!
!    Input, real XVAL, a value of X at which the polynomial
!    is to be evaluated.
!
!    Output, real YVAL, the value of the polynomial at XVAL.
!
  implicit none

  integer ntab

  real diftab(ntab)
  integer i
  real xtab(ntab)
  real xval
  real yval

  yval = diftab(ntab)
  do i = 1, ntab-1
    yval = diftab(ntab-i) + ( xval - xtab(ntab-i) ) * yval
  end do

  return
end
subroutine digit_to_ch ( digit, c )

!*******************************************************************************
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
function dot0_3d ( x0, y0, z0, x1, y1, z1, x2, y2, z2 )

!*******************************************************************************
!
!! DOT0_3D computes the dot product of (P1-P0) and (P2-P0) in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, Z0, the coordinates of the point P0.
!
!    Input, real X1, Y1, Z1, the coordinates of the point P1.
!
!    Input, real X2, Y2, Z2, the coordinates of the point P2.
!
!    Output, real DOT0_3D, the dot product of (P1-P0) and (P2-P0).
!
  implicit none

  real dot0_3d
  real x0
  real x1
  real x2
  real y0
  real y1
  real y2
  real z0
  real z1
  real z2

  dot0_3d = ( x1 - x0 ) * ( x2 - x0 ) + ( y1 - y0 ) * ( y2 - y0 ) + &
            ( z1 - z0 ) * ( z2 - z0 )

  return
end
subroutine drwplt ( ierror, ival, l3d, ldat, lframe, ndat, nval, &
  nxgrid, nygrid, plhite, plwide, theta, title, vscale, xdat, xmax, &
  xmaxw, xmin, xminw, xval, ydat, ymax, ymaxw, ymin, yminw, yval )

!***********************************************************************
!
!! DRWPLT displays the graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
!    Input, integer IVAL(VAL_MAX).  IVAL(I) is 1 if point I
!    is to be connected to point I+1.
!
!    Input, integer NVAL, the number of plot points.
!
  implicit none

  real, parameter :: hletr = 0.1E+00

  integer ndat
  integer nval

  integer i
  integer icom
  integer ierror
  integer iplt1
  integer iplt2
  integer ipoint
  character isay
  integer ival(nval)
  integer ival1
  integer ival2
  integer ixplt1
  integer ixplt2
  integer iyplt1
  integer iyplt2
  logical l3d
  logical ldat
  logical s_eqi
  logical lframe
  integer marray
  integer ndraw
  integer nxgrid
  integer nygrid
  real r_pi
  real plhite
  real plwide
  real space
  real start
  real stheta
  real theta
  character ( len = 80 ) title
  real vscale
  real x1
  real x1max
  real x1min
  real x1rang
  real x2
  real x2max
  real x2min
  real x2rang
  real xdat(ndat)
  real xdif
  real xdraw(6)
  real, save :: xmarg = 0.05E+00
  real xmax
  real xmaxw
  real xmin
  real xminw
  real xplt1
  real xplt2
  real xrange
  real xskal
  real xval(nval)
  real xvpos1
  real xvpos2
  real y1
  real y1max
  real y1min
  real y1rang
  real y2
  real y2max
  real y2min
  real y2rang
  real ydat(ndat)
  real ydif
  real ydraw(6)
  real, save :: ymargb = 0.05E+00
  real, save :: ymargt = 0.05E+00
  real ymax
  real ymaxw
  real ymin
  real yminw
  real yplt1
  real yplt2
  real yrange
  real yskal
  real yval(nval)
  real yvpos1
  real yvpos2
!
!  ANYPLT stuff.
!
  character ( len = 80 ) carray
  common /anycom/ iplt1,iplt2,ixplt1,ixplt2,iyplt1, &
                     iyplt2,marray,xplt1,xplt2,yplt1,yplt2
  common /anychr/ carray
!
  stheta = sin ( r_pi ( ) * theta / 180.0E+00 )
!
!  Set the data range.
!
  xmin = min ( minval ( xval(1:nval) ), minval ( xdat(1:ndat) ) )
  xmax = max ( maxval ( xval(1:nval) ), maxval ( xdat(1:ndat) ) )
  ymin = min ( minval ( yval(1:nval) ), minval ( ydat(1:ndat) ) )
  ymax = max ( maxval ( yval(1:nval) ), maxval ( ydat(1:ndat) ) )
!
!  Set the window range.
!
  xmaxw = xmax + 0.1E+00 * ( xmax - xmin )
  xminw = xmin - 0.1E+00 * ( xmax - xmin )
  ymaxw = ymax + 0.1E+00 * ( ymax - ymin )
  yminw = ymin - 0.1E+00 * ( ymax - ymin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter Y" if you want the same scale on both axes.'

  read ( *, '(a)' ) isay

  if ( ierror /= 0 ) then
    return
  end if

  if ( s_eqi ( isay(1:1),'Y' ) ) then

    xdif = xmaxw - xminw
    ydif = ymaxw - yminw

    if ( xdif > ydif ) then
      ymaxw = ymaxw + 0.5E+00 * ( xdif - ydif )
      yminw = yminw - 0.5E+00 * ( xdif - ydif )
    else
      xmaxw = xmaxw + 0.5E+00 * ( ydif - xdif )
      xminw = xminw - 0.5E+00 * ( ydif - xdif )
    end if

  end if
!
!  Set the X and Y scaling factors.
!
  xrange = xmaxw-xminw
  if ( xrange  <=  0.0E+00 ) then
    xrange = plwide-2.0E+00 *xmarg
  end if
  xskal = (plwide-2.0E+00 *xmarg)/xrange

  yrange = ymaxw-yminw
  if ( yrange  <=  0.0E+00 ) then
    yrange = plhite-ymargb-ymargt
  end if
  yskal = (plhite-ymargb-ymargt)/yrange
!
!  Set minimum, range, and maximum for X and Y physical dimensions.
!
  x1min = xskal*xminw-xmarg
  x1rang = plwide
  x1max = x1min+x1rang
  y1min = yskal*yminw-ymargb
  y1rang = plhite
  y1max = y1min+y1rang
  x2min = xskal*xminw
  x2rang = plwide-2.0E+00 *xmarg
  x2max = x2min+x2rang
  y2min = yskal*yminw
  y2rang = plhite-ymargb-ymargt
  y2max = y2min+y2rang
  icom = 2
  call anyplt ( icom )
!
!  Set the physical range of the plot.
!  Set the origin at minimum X, minimum Y.
!
  icom = 3
  xplt1 = x1min
  xplt2 = x1rang
  yplt1 = y1min
  yplt2 = y1rang
  call anyplt ( icom )
!
!  Draw a frame around the plot.
!
  if ( lframe ) then

    write ( *, '(a)' ) 'DEBUG: Framing the plot.'

    icom = 4
    xplt1 = x1min
    yplt1 = y1min
    call anyplt ( icom )

    icom = 5
    yplt1 = y1max
    call anyplt ( icom )

    xplt1 = x1max
    call anyplt ( icom )

    yplt1 = y1min
    call anyplt ( icom )

    xplt1 = x1min
    call anyplt ( icom )

  end if
!
!  Draw the grid lines.
!
  if ( nxgrid > 0 ) then

    do i = 1, nxgrid
      call rvec_even_select ( x2min, x2max, nxgrid, i, xplt1 )
      yplt1 = y2max
      icom = 4
      call anyplt ( icom )
      yplt1 = y2min
      icom = 5
      call anyplt ( icom )
    end do

  end if

  if ( nygrid > 0 ) then

    do i = 1, nygrid
      xplt1 = x2min
      call rvec_even_select ( y2min, y2max, nygrid, i, yplt1 )
      icom = 4
      call anyplt ( icom )
      xplt1 = x2max
      icom = 5
      call anyplt ( icom )
    end do

  end if
!
!  Label the plot
!
  if ( title /= ' ' ) then
    carray = title
    call s_blanks_delete(carray)
    marray = len_trim ( carray )
    xplt2 = hletr
    xplt1 = 0.5E+00 *(x1min+x1max)
    yplt1 = y1max-2.0E+00 *hletr
    yplt2 = 0.0E+00
    icom = 7
    call anyplt ( icom )
  end if

  if ( .not. l3d ) then
    write ( carray, '(a,f12.4,a,f12.4)' ) 'XMin=', xminw, ', XMax=', xmaxw
  else
    carray = ' '
  end if

  call s_blanks_delete(carray)
  marray = len_trim ( carray )
  xplt2 = hletr
  xplt1 = 0.5E+00 *(x1min+x1max)
  yplt1 = y1max-3.0E+00 *hletr
  icom = 7
!     call anyplt ( icom )

  if ( .not. l3d ) then
    write ( carray, '(a,f12.4,a,f12.4)' ) 'YMin=', yminw, ', YMax=', ymaxw
  else
    carray = ' '
  end if

  call s_blanks_delete(carray)
  marray = len_trim ( carray )
  xplt1 = 0.5E+00 *(x1min+x1max)
  yplt1 = y1max-4.0E+00 *hletr
  icom = 7
!     call anyplt ( icom )

  if ( l3d ) then
    carray = ' '
    call s_blanks_delete(carray)
    marray = len_trim ( carray )
    xplt1 = 0.5E+00 *(x1min+x1max)
    yplt1 = y1max-5.0E+00 *hletr
    icom = 7
!       call anyplt ( icom )
    write ( carray, '(a,f12.4)' ) '3D presentation angle=', theta
    call s_blanks_delete(carray)
    marray = len_trim ( carray )
    xplt1 = 0.5E+00 *(x1min+x1max)
    yplt1 = y1max-6.0E+00 *hletr
    icom = 7
!       call anyplt ( icom )
  end if
!
!  Draw the 3D axes.
!
  if ( l3d ) then

    xplt1 = x1min+0.1E+00
    yplt1 = y1min+0.1E+00
    icom = 4
    call anyplt ( icom )

    xplt1 = x1min+0.1E+00-0.05E+00*stheta
    yplt1 = y1min+0.1E+00-0.05E+00*stheta
    icom = 5
    call anyplt ( icom )

    carray = 'x'
    marray = 1
    xplt1 = x1min+0.1E+00-0.06E+00*stheta
    yplt1 = y1min+0.1E+00-0.06E+00*stheta
    icom = 7
    call anyplt ( icom )

    xplt1 = x1min+0.1E+00
    yplt1 = y1min+0.1E+00
    icom = 4
    call anyplt ( icom )

    xplt1 = x1min+0.15E+00
    yplt1 = y1min+0.1E+00
    icom = 5
    call anyplt ( icom )

    carray = 'y'
    marray = 1
    xplt1 = x1min+0.16E+00
    yplt1 = y1min+0.1E+00
    icom = 7
    call anyplt ( icom )

    xplt1 = x1min+0.1E+00
    yplt1 = y1min+0.1E+00
    icom = 4
    call anyplt ( icom )

    xplt1 = x1min+0.1E+00
    yplt1 = y1min+0.15E+00
    icom = 5
    call anyplt ( icom )

    carray = 'z'
    marray = 1
    xplt1 = x1min+0.1E+00
    yplt1 = y1min+0.16E+00
    icom = 7
    call anyplt ( icom )
  end if
!
!  Begin drawing
!
  ipoint = 1
  xvpos2 = xskal*xval(ipoint)
  yvpos2 = yskal*yval(ipoint)
  ival2 = ival(ipoint)
  x1 = xvpos2
  y1 = yvpos2
  ival1 = ival2
  x2 = xvpos2
  y2 = yvpos2

  icom = 4
  xplt1 = x2
  yplt1 = y2
  call anyplt ( icom )

  do while ( ipoint < nval )

    ipoint = ipoint+1

    ival1 = ival2
    ival2 = ival(ipoint)
    xvpos1 = xvpos2
    yvpos1 = yvpos2
    xvpos2 = xskal*xval(ipoint)
    yvpos2 = yskal*yval(ipoint)
    x1 = xvpos1
    y1 = yvpos1
    x2 = xvpos2
    y2 = yvpos2

    if ( ival1 == 1 ) then

      icom = 4
      xplt1 = x1
      yplt1 = y1
      call anyplt ( icom )
      icom = 5
      xplt1 = x2
      yplt1 = y2
      call anyplt ( icom )
!
!  Draw an arrow representing a vector quantity.
!
    else if ( ival1 == 2 ) then

      x2 = x1 + vscale * ( x2 - x1 )
      y2 = y1 + vscale * ( y2 - y1 )

      call arrow ( x1, y1, x2, y2, ndraw, xdraw, ydraw )

      if ( ndraw > 0 ) then
        icom = 4
        xplt1 = xdraw(1)
        yplt1 = ydraw(1)
        call anyplt ( icom )
        do i = 2, ndraw
          icom = 5
          xplt1 = xdraw(i)
          yplt1 = ydraw(i)
          call anyplt ( icom )
        end do

      end if
!
!         uu = vscale*(x2-x1)
!         vv = vscale*(y2-y1)
!         tnorm = sqrt(uu**2+vv**2)
!
!         if ( tnorm > 0.0E+00 ) then
!
!           thet = 0.5 * r_pi ( ) - atan2(2.0,1.0)
!           alph = atan2(vv,uu)
!           del = sqrt(5.0)*tnorm/3.0E+00
!
!           u1 = x1+del*cos(alph-thet)
!           v1 = y1+del*sin(alph-thet)
!
!           u2 = x1+del*cos(alph+thet)
!           v2 = y1+del*sin(alph+thet)
!
!           icom = 5
!           xplt1 = x1+vscale*(x2-x1)
!           yplt1 = y1+vscale*(y2-y1)
!           call anyplt ( icom )
!
!           icom = 4
!           xplt1 = u1
!           yplt1 = v1
!           call anyplt ( icom )
!
!           icom = 5
!           xplt1 = x1+vscale*(x2-x1)
!           yplt1 = y1+vscale*(y2-y1)
!           call anyplt ( icom )
!
!           icom = 5
!           xplt1 = u2
!           yplt1 = v2
!           call anyplt ( icom )
!
!         end if
!
    end if

  end do
!
!  Mark the data points.
!
  if ( ldat .and. ndat <= 100 ) then

    do i = 1, ndat

      x1 = xskal*xdat(i)
      y1 = yskal*ydat(i)

      if ( x2min  <=  x1 .and. x1 <= x2max .and. &
           y2min  <=  y1 .and. y1 <= y2max ) then
        xplt1 = x1
        yplt1 = y1
        icom = 11
        call anyplt ( icom )
      end if

    end do

  end if

  return
end
subroutine rsftdw ( l, u, k, lda, a, map )

!***********************************************************************
!
!! RSFTDW shifts A(*,MAP(L)) down a heap of size U.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer L, U, the lower and upper indices of part of the heap.
!
!    Input, integer K, the spatial dimension of the points.
!
!    Input, integer LDA, the leading dimension of A in the calling routine.
!
!    Input, real A(LDA,N); A(I,J) contains the I-th coordinate
!    of point J.
!
!    Input/output, integer MAP(N).
!    On input, the points of A with indices MAP(1), MAP(2), ...,
!    MAP(N) are to be sorted
!
!    On output, MAP contains a permutation of its input values, with the
!    property that, lexicographically,
!      A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N))
!
  implicit none

  integer lda

  real a(lda,*)
  integer i
  integer j
  integer k
  integer l
  integer map(*)
  logical rless
  integer t
  integer u

  i = l
  j = 2 * i
  t = map(i)

  do while ( j <= u )

    if ( j < u ) then
      if ( rless ( k, a(1,map(j)), a(1,map(j+1)) ) ) then
        j = j + 1
      end if
    end if

    if ( rless ( k, a(1,map(j)), a(1,t) ) ) then
      exit
    end if

    map(i) = map(j)
    i = j
    j = 2 * i

  end do

  map(i) = t

  return
end
subroutine rtris2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierr )

!***********************************************************************
!
!! RTRIS2 constructs a Delaunay triangulation of 2D vertices.
!
!  Discussion:
!
!    The routine constructs the Delaunay triangulation of a set of 2D vertices
!    using an incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (X,Y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer NPT, the number of vertices.
!
!    Input, integer MAXST, the maximum size available for the STACK array; 
!    should be about NPT to be safe, but MAX(10,2*LOG2(NPT)) is usually enough.
!
!    Input, real VCL(2,NPT), the coordinates of the vertices.
!
!    Input/output, integer IND(NPT), the indices in VCL of the vertices 
!    to be triangulated.  On output, IND has been permuted by the sort.
!
!    Output, integer NTRI, the number of triangles in the triangulation; 
!    NTRI is equal to 2*NPT - NB - 2, where NB is the number of boundary 
!    vertices.
!
!    Output, integer TIL(3,NTRI), the triangle incidence list.
!    The elements are indices of VCL.  The vertices of the triangles are 
!    in counter clockwise order
!
!    Output, integer TNBR(3,NTRI), the triangle neighbor list.
!    Positive elements are indices of TIL; negative elements are used for links
!    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
!    where I, J = triangle, edge index; TNBR(J,I) refers to
!    the neighbor along edge from vertex J to J+1 (mod 3)
!
!    Workspace, integer STACK(MAXST), used for a stack of triangles for which
!    circumcircle test must be made.
!
!    Output, integer IERR, an error flag, nonzero if an error occurred.
!
  implicit none

  integer maxst
  integer npt

  real cmax
  integer e
  integer i
  integer ierr
  integer ind(npt)
  integer j
  integer k
  integer l
  integer ledg
  integer lr
  integer lrline
  integer ltri
  integer m
  integer m1
  integer m2
  integer n
  integer ntri
  integer redg
  integer rtri
  integer stack(maxst)
  integer t
  real temp
  integer til(3,npt*2)
  integer tnbr(3,npt*2)
  real tol
  integer top
  real vcl(2,*)
!
  ierr = 0
  tol = 100.0E+00 * epsilon ( tol )
!
!  Sort the vertices by increasing (x,y) and obtain initial triangle(s).
!
  call rhpsrt ( 2, npt, 2, vcl, ind )

  m1 = ind(1)

  do i = 2, npt

    m = m1
    m1 = ind(i)

    k = 0

    do j = 1, 2

      cmax = max ( abs ( vcl(j,m) ), abs ( vcl(j,m1) ) )

      if ( abs ( vcl(j,m) - vcl(j,m1) ) > tol * cmax .and. cmax > tol ) then
        k = j
        exit
      end if

    end do

    if ( k == 0 ) then
      ierr = 224
      return
    end if

  end do

  m1 = ind(1)
  m2 = ind(2)
  j = 3

  do

    if ( j > npt ) then
      ierr = 225
      return
    end if

    m = ind(j)
    lr = lrline ( vcl(1,m), vcl(2,m), vcl(1,m1), vcl(2,m1), vcl(1,m2), &
      vcl(2,m2), 0.0E+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1
  
  end do

  ntri = j - 2

  if ( lr == -1 ) then

    til(1,1) = m1
    til(2,1) = m2
    til(3,1) = m
    tnbr(3,1) = -3

    do i = 2, ntri

      m1 = m2
      m2 = ind(i+1)
      til(1,i) = m1
      til(2,i) = m2
      til(3,i) = m
      tnbr(1,i-1) = -3 * i
      tnbr(2,i-1) = i
      tnbr(3,i) = i - 1

    end do

    tnbr(1,ntri) = -3 * ntri - 1
    tnbr(2,ntri) = -5
    ledg = 2
    ltri = ntri

  else

    til(1,1) = m2
    til(2,1) = m1
    til(3,1) = m
    tnbr(1,1) = -4

    do i = 2, ntri
      m1 = m2
      m2 = ind(i+1)
      til(1,i) = m2
      til(2,i) = m1
      til(3,i) = m
      tnbr(3,i-1) = i
      tnbr(1,i) = -3 * i - 3
      tnbr(2,i) = i - 1
    end do

    tnbr(3,ntri) = -3 * ntri
    tnbr(2,1) = -3 * ntri - 2
    ledg = 2
    ltri = 1

  end if
!
!  Insert vertices one at a time from outside convex hull, determine
!  visible boundary edges, and apply diagonal edge swaps until
!  Delaunay triangulation of vertices (so far) is obtained.
!
  top = 0

  do i = j+1, npt

    m = ind(i)
    m1 = til(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = til(ledg+1,ltri)
    else
      m2 = til(1,ltri)
    end if

    lr = lrline ( vcl(1,m), vcl(2,m), vcl(1,m1), vcl(2,m1), &
      vcl(1,m2), vcl(2,m2), 0.0E+00 )

    if ( lr > 0 ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -tnbr(ledg,ltri)
      rtri = l / 3
      redg = mod(l,3) + 1
    end if

    call vbedg ( vcl(1,m), vcl(2,m), vcl, til, tnbr, ltri, ledg, rtri, redg )

    n = ntri + 1
    l = -tnbr(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -tnbr(e,t)
      m2 = til(e,t)

      if ( e <= 2 ) then
        m1 = til(e+1,t)
      else
        m1 = til(1,t)
      end if

      ntri = ntri + 1
      tnbr(e,t) = ntri
      til(1,ntri) = m1
      til(2,ntri) = m2
      til(3,ntri) = m
      tnbr(1,ntri) = t
      tnbr(2,ntri) = ntri - 1
      tnbr(3,ntri) = ntri + 1
      top = top + 1

      if ( top > maxst ) then
        ierr = 8
        return
      end if

      stack(top) = ntri

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    tnbr(ledg,ltri) = -3 * n - 1
    tnbr(2,n) = -3 * ntri - 2
    tnbr(3,ntri) = -l
    ltri = n
    ledg = 2

    call swapec ( m, top, maxst, ltri, ledg, vcl, til, tnbr, stack, ierr )

    if ( ierr /= 0 ) then
      return
    end if

  end do

  return
end
function enorm0_3d ( x0, y0, z0, x1, y1, z1 )

!*******************************************************************************
!
!! ENORM0_3D computes the Euclidean norm of (P1-P0) in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, Z0, X1, Y1, Z1, the coordinates of the points
!    P0 and P1.
!
!    Output, real ENORM0_3D, the Euclidean norm of (P1-P0).
!
  implicit none

  real enorm0_3d
  real x0
  real x1
  real y0
  real y1
  real z0
  real z1

  enorm0_3d = sqrt ( ( x1 - x0 )**2 + ( y1 - y0 )**2 + ( z1 - z0 )**2 )

  return
end
function enorm_nd ( n, x )

!*******************************************************************************
!
!! ENORM_ND computes the Euclidean norm of a vector in ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Input, real X(N), the coordinates of the vector.
!
!    Output, real ENORM_ND, the Euclidean norm of the vector.
!
  implicit none

  integer n

  real enorm_nd
  integer i
  real x(n)

  enorm_nd = sqrt ( sum ( x(1:n)**2 ) )

  return
end
subroutine evalf ( dat, ierror, idat, maxdat, mdat, namdat, ndat )

!***********************************************************************
!
!! EVALF gets a formula from the user and evaluates it.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real DAT(MAXDAT,MDAT).
!
!    On input, DAT(I,J) contains the value of variable J at point I,
!    for all the independent variables needed.  However, column
!    IDAT is probably not defined.
!
!    On output, DAT(I,IDAT) contains the value of variable IDAT
!    at point I, for I = 1 to NDAT, as determined by a formula
!    entered by the user.  The formula may involve the independent
!    variables as stored in DAT.
!
!    Output, integer IERROR, error flag.  If IERROR is nonzero, then
!    an error occurred in reading, compiling, or evaluating the
!    formula.  The entries of column IDAT of DAT have not been
!    properly set.
!
!    Input, integer MAXDAT, the first dimension of DAT, which is
!    at least equal to NDAT.
!
!    Input, integer MDAT, the number of columns of DAT, and the number
!    of dependent and independent variables.
!
!    Input, character ( len = MDAT ) NAMDAT, a string containing the
!    one character names for the variables.  The order of the names
!    must correspond to the order of the corresponding columns in DAT.
!
!    Input, integer NDAT, the number of points at which the values of
!    the independent variables were given, and at which the value
!    of the dependent variable is desired.
!
  implicit none

  integer maxdat
  integer, parameter :: maxrpn = 100
  integer ndat

  real dat(maxdat,9)
  integer i
  integer, save :: icall = 0
  integer idat
  integer ierror
  integer ifrm
  character ( len = 80 ) infix
  integer irpn(maxrpn)
  integer j
  integer lenc
  integer lenchr
  character ( len = 80 ) line
  integer mdat
  character ( len = 6 ) namdat
  real value

  icall = icall + 1

  if ( icall == 1 ) then
!
!  Initialize the compiler.
!
    ifrm = 0
    call comrpn ( 'I', ierror, ifrm, infix, irpn, maxrpn, ' ', value )
!
!  Declare variables to the compiler.
!
    line = 'rstuvwxyz'
    lenc = len_trim ( line )
    do i = 1, lenc
      call comrpn ( 'A', ierror, ifrm, infix, irpn, maxrpn, line(i:i), value )
    end do

  end if
!
!  Get the formula.
!
  do

    write ( *, '(a)' ) 'Enter formula for ' // namdat(idat:idat)
    read ( *, '(a)' ) infix
!
!  Compile the formula.
!
    ifrm = 1
    call comrpn ( 'F', ierror, ifrm, infix, irpn, maxrpn, ' ', value )

    if ( ierror == 0 ) then
      exit
    end if

    write ( *, '(a)' ) 'Your formula could not be compiled.'

  end do
!
!  Evaluate the formula at every data point.
!
  do i = 1, ndat
!
!  Tell the compiler the current value of all variables.
!
    do j = 1, mdat

      call comrpn ( 'V', ierror, ifrm, infix, irpn, maxrpn, namdat(j:j), dat(i,j) )

    end do
!
!  Now evaluate the formula.
!
    ifrm = 1
    call comrpn ( 'E', ierror, ifrm, infix, irpn, maxrpn, ' ', dat(i,idat) )

  end do

  return
end
subroutine ffxy ( ierror, ival, val_max, nval, xval, yval )

!***********************************************************************
!
!! FFXY sets up data for F(X,Y)=0 plots.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
!    Input/output, integer IVAL(VAL_MAX).  IVAL(I) is 1 if point I
!    is to be connected to point I+1.
!
!    Input, integer VAL_MAX, the maximum number of points that can be 
!    used for plots.
!
!    Input/output, integer NVAL, the current number of plot points.
!
  implicit none

  integer, parameter :: maxfix = 80
  integer, parameter :: maxrpn = 300
  integer val_max

  real dir
  integer ierror
  integer ifrm
  character ( len = 80 ) infix
  integer irpn(maxrpn)
  integer ival(val_max)
  integer ndraw
  integer nval
  real value
  real xstart
  real xval(val_max)
  real ystart
  real yval(val_max)
!
!  Initialize the compiler.
!
  ifrm = 0
  call comrpn ( 'I', ierror, ifrm, infix, irpn, maxrpn, ' ', value )
!
!  Declare X and Y to the compiler.
!
  call comrpn ( 'A', ierror, ifrm, infix, irpn, maxrpn, 'X', value )

  call comrpn ( 'A', ierror, ifrm, infix, irpn, maxrpn, 'Y', value )
!
  write ( *, '(a)' ) 'Enter starting point (X, Y)'
  read ( *, * ) xstart, ystart

  write ( *, '(a)' ) 'Enter direction of +1 or -1'
  read ( *, * ) dir

  if ( dir /= -1.0E+00 ) then
    dir = 1.0E+00
  end if

  do

    write ( *, '(a)' ) 'Enter formula for F(X,Y)'
    read ( *, '(a)' ) infix

    ifrm = 1
    call comrpn ( 'F', ierror, ifrm, infix, irpn, maxrpn, ' ', value )

    if ( ierror == 0 ) then
      exit
    end if

  end do

  do

    write ( *, '(a)' ) 'Enter formula for dF/dX'
    read ( *, '(a)' ) infix

    ifrm = 2
    call comrpn ( 'F', ierror, ifrm, infix, irpn, maxrpn, ' ', value )
    if ( ierror == 0 ) then
      exit
    end if

  end do

  do

    write ( *, '(a)' ) 'Enter formula for dF/dY'
    read ( *, '(a)' ) infix

    ifrm = 3
    call comrpn ( 'F', ierror, ifrm, infix, irpn, maxrpn, ' ', value )

    if ( ierror == 0 ) then
      exit
    end if

  end do

  write ( *, '(a)' ) 'Enter number of points to draw.'
  read ( *, * ) ndraw

  call fxy ( dir, ierror, infix, irpn, ival, maxfix, &
    maxrpn, val_max, ndraw, nval, xstart, xval, ystart, yval )

  return
end
subroutine file_delete ( file_name )
!
!*******************************************************************************
!
!! FILE_DELETE deletes a named file if it exists.
!
!  Discussion:
!
!    You might want to call this routine to get rid of any old copy
!    of a file, before trying to open a new copy with the OPEN argument:
!      status = 'new'.
!
!    It's not always safe to open a file with " STATUS = 'UNKNOWN' ".
!    For instance, on the SGI, the most recent version of the FORTRAN
!    compiler seems to go crazy when I open an unformatted direct
!    access file this way.  It creates an enormous file (of somewhat
!    random size).  The problem goes away if I delete any old copy
!    using this routine, and then open a fresh copy with
!    " STATUS = 'NEW' ".  It's a scary world.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to be deleted.
!
  implicit none

  character ( len = * ) file_name
  integer ios
  integer iunit
  logical lfile
!
!  Does the file exist?
!
  inquire ( file = file_name, exist = lfile )

  if ( .not. lfile ) then
    return
  end if
!
!  Get a free unit number.
!
  call get_unit ( iunit )

  if ( iunit == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE: Warning!'
    write ( *, '(a)' ) '  A free FORTRAN unit could not be found.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_DELETE: deleting old version of ' // &
    trim ( file_name )

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE: Warning!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( file_name )
    return
  end if

  close ( unit = iunit, status = 'delete' )

  return
end
subroutine funscl ( angle_to_radian, arg1, arg2, ierror, result, rtol, sym )
!
!*******************************************************************************
!
!! FUNSCL evaluates a scalar function of one or more scalar arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ARG2, ARG2, the values of the arguments of the function.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Output, real RESULT, the value of the scalar function.
!
!    Input, real RTOL, the rounding tolerance.
!
!    Input, character ( len = * ) SYM, the symbolic name of the
!    function to be evaluated.
!
  implicit none

  real, parameter :: degrad = 3.14159265358979E+00 / 180.0E+00

  real angle_to_radian
  real arg1
  real arg2
  integer i_lcm
  integer iarg
  integer iarg1
  integer iarg2
  integer ierror
  integer i_gcd
  real result
  real rtol
  logical s_eqi
  real sarg
  character ( len = * ) sym
  real temp
!
  ierror = 0
  result = 0.0E+00
!
!  Addition.
!
  if ( sym == '+' ) then

    result = arg1 + arg2
!
!  Subtraction.
!
  else if ( sym == '-' ) then

    result = arg1 - arg2
!
!  Division.
!
  else if ( sym == '/' ) then

    if ( arg2 == 0.0E+00 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FUNSCL - Error!'
      write ( *, '(a)' ) '  Attempt to divide by 0!'
      return
    end if

    result = arg1 / arg2
!
!  Multiplication.
!
  else if ( sym == '*' ) then
    result = arg1 * arg2
!
!  Exponentiation.
!
  else if ( sym == '**' .or. sym == '^' ) then

    if ( arg1 == 0.0E+00 .and. arg2 == 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FUNSCL - Error!'
      write ( *, '(a)' ) '  Attempt to compute 0**0 !'
      ierror = 1
    else if ( arg1 < 0.0E+00 .and. real(int(arg2)) /= arg2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FUNSCL - Error!'
      write ( *, '(a)' ) '  Illegal exponentiation:'
      write ( *, '(a,g14.6,a,g14.6)' ) '  ', arg1, ' ** ', arg2
      ierror = 1
    else

      sarg = 1.0E+00

      if ( arg1 < 0.0E+00 ) then
        iarg = int ( arg2 )
        if ( 2 * ( iarg / 2 ) /= iarg ) then
          sarg = - 1.0E+00
        end if
      end if

      result = sarg * abs ( arg1 )**arg2

    end if
!
!  Assignment.
!
  else if ( sym == '=' ) then
    result = arg2
!
!  Absolute value.
!
  else if ( s_eqi ( sym, 'ABS' ) ) then
    result = abs ( arg1 )
!
!  Arc Cosine.
!
  else if ( s_eqi ( sym, 'ACOS' ) ) then

    if ( arg1 < -1.0E+00 .or. arg1 > 1.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FUNSCL - Error!'
      write ( *, '(a,g14.6)' ) '  Illegal inverse cosine of ', arg1
      ierror = 1
      return
    end if

    result = acos ( arg1 ) / angle_to_radian
!
!  Arc Sine.
!
  else if ( s_eqi ( sym, 'ASIN' ) ) then

    if ( arg1 < - 1.0E+00 .or. arg1 > 1.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FUNSCL - Error!'
      write ( *, '(a,g14.6)' ) '  Illegal inverse sine of ', arg1
      ierror = 1
    end if

    result = asin ( arg1 ) / angle_to_radian
!
!  Arc Tangent.
!
  else if ( s_eqi ( sym, 'ATAN' ) ) then
    result = atan ( arg1 ) / angle_to_radian
!
!  Arc Tangent of ratio.
!
  else if ( s_eqi ( sym, 'ATAN2' ) ) then

    if ( arg1 == 0.0E+00 .and. arg2 == 0.0E+00 ) then
      result = 0.0E+00
    else
      result = atan2 ( arg1, arg2 ) / angle_to_radian
    end if
!
!  Cosine.
!
  else if ( s_eqi ( sym, 'COS' ) ) then
    result = cos ( angle_to_radian * arg1 )
!
!  Cotangent of angle.
!
  else if ( s_eqi ( sym, 'COT' ) ) then
    result = cos ( angle_to_radian * arg1 ) / sin ( angle_to_radian * arg1 )
!
!  Hyperbolic cosine.
!
  else if ( s_eqi ( sym, 'COSH' ) ) then
    result = cosh ( angle_to_radian * arg1 )
!
!  Cosecant of angle.
!
  else if ( s_eqi ( sym, 'CSC' ) ) then
    result = 1.0E+00 / sin ( angle_to_radian * arg1 )
!
!  Determinant (of a scalar).
!
  else if ( s_eqi ( sym, 'DET' ) ) then
    result = arg1
!
!  Diagonal (of a scalar).
!
  else if ( s_eqi ( sym, 'DIAG' ) ) then
    result = arg1
!
!  Exponential.
!
  else if ( s_eqi ( sym, 'EXP' ) ) then
    result = exp ( arg1 )
!
!  Greatest common factor.
!
  else if ( s_eqi ( sym, 'GCD' ) ) then
    iarg1 = nint ( arg1 )
    iarg2 = nint ( arg2 )
    result = real ( i_gcd ( iarg1, iarg2 ) )
!
!  Hilbert.
!
  else if ( s_eqi ( sym, 'HILBERT' ) ) then
    result = 1.0E+00
!
!  Hilbert inverse.
!
  else if ( s_eqi ( sym, 'HILBINV' ) ) then
    result = 1.0E+00
!
!  Identity.
!
  else if ( s_eqi ( sym, 'ID' ) ) then
    result = 1.0E+00
!
!  Integer value.
!
  else if ( s_eqi ( sym, 'INT' ) ) then
    result = real ( int ( arg1 ) )
!
!  Inverse (of a scalar).
!
  else if ( s_eqi ( sym, 'INV' ) ) then

    if ( arg1 == 0.0E+00 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FUNSCL - Error!'
      write ( *, '(a)' ) '  Attempt to compute 1/0!'
      return
    end if

    result = 1.0E+00 / arg1
!
!  Least common multiple.
!
  else if ( s_eqi ( sym, 'LCM' ) ) then
    iarg1 = nint ( arg1 )
    iarg2 = nint ( arg2 )
    result = real ( i_lcm ( iarg1, iarg2 ) )
!
!  Natural logarithm.
!
  else if ( s_eqi ( sym, 'ALOG' ) .or. s_eqi ( sym, 'LN' ) .or. &
            s_eqi ( sym, 'LOG' ) ) then

    if ( arg1 <= 0.0E+00 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FUNSCL - Error!'
      write ( *, '(a,g14.6)' ) '  Illegal LOG of ', arg1
      return
    end if

    result = log ( arg1 )
!
!  Logarithm base 10.
!
  else if ( s_eqi ( sym, 'ALOG10' ) .or. s_eqi ( sym, 'LOG10' ) ) then

    if ( arg1 <= 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FUNSCL - Error!'
      write ( *, '(a,g14.6)' ) '  Illegal LOG10 of ', arg1
      ierror = 1
      return
    end if

    result = log10 ( arg1 )
!
!  Logarithm base 2.
!
  else if ( s_eqi ( sym, 'LOG2' ) ) then

    if ( arg1 <= 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FUNSCL - Error!'
      write ( *, '(a,g14.6)' ) '  Illegal LOG2 of ', arg1
      ierror = 1
      return
    end if

    result = log ( arg1 ) / log ( 2.0E+00 )
!
!  Maximum of two values.
!
  else if ( s_eqi ( sym, 'MAX' ) ) then

    if ( arg1 > arg2 ) then
      result = arg1
    else
      result = arg2
    end if
!
!  Minimum of two values.
!
  else if ( s_eqi ( sym, 'MIN' ) ) then

    if ( arg1 < arg2 ) then
      result = arg1
    else
      result = arg2
    end if
!
!  Negation
!
  else if ( s_eqi ( sym, 'NEG' ) ) then
    result = - arg1
!
!  Nearest integer value
!
  else if ( s_eqi ( sym, 'NINT' ) ) then
    result = anint ( arg1 )
!
!  NORM0, NORM1, NORM2, NORME or NORMF of a scalar.
!
  else if ( s_eqi ( sym, 'NORM0' ) .or. s_eqi ( sym, 'NORM1' ) .or. &
            s_eqi ( sym, 'NORM2' ) .or. s_eqi ( sym, 'NORME' ) .or. &
            s_eqi ( sym, 'NORMF' ) ) then
    result = abs ( arg1 )
!
!  POLY (scalar coefficient array, means constant polynomial).
!
  else if ( s_eqi ( sym, 'POLY' ) ) then
    result = arg1
!
!  Random value.
!
  else if ( s_eqi ( sym, 'RAN' ) ) then

    call random_number ( result )
!
!  Rounding.
!
  else if ( s_eqi ( sym, 'ROUND' ) ) then

    if ( abs ( arg1 ) < rtol ) then
      result = 0.0E+00
    else
      result = arg1
    end if
!
!  Secant.
!
  else if ( s_eqi ( sym, 'SEC' ) ) then
    result = 1.0E+00 / cos ( angle_to_radian * arg1 )
!
!  Sine.
!
  else if ( s_eqi ( sym, 'SIN' ) ) then
    result = sin ( angle_to_radian * arg1 )
!
!  Hyperbolic sine.
!
  else if ( s_eqi ( sym, 'SINH' ) ) then
    result = sinh ( angle_to_radian * arg1 )
!
!  Square root.
!
  else if ( s_eqi ( sym, 'SQRT' ) ) then

    if ( arg1 < 0.0E+00 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FUNSCL - Error!'
      write ( *, '(a,g14.6)' ) '  Illegal SQRT of ', arg1
      return
    end if

    result = sqrt ( arg1 )
!
!  Step function.
!
  else if ( s_eqi ( sym, 'STEP' ) ) then

    if ( arg1 < 0.0E+00 ) then
      result = 0.0E+00
    else
      result = 1.0E+00
    end if
!
!  Tangent.
!
  else if ( s_eqi ( sym, 'TAN' ) ) then
    result = tan ( angle_to_radian * arg1 )
!
!  Hyperbolic tangent.
!
  else if ( s_eqi ( sym, 'TANH' ) ) then
    result = tanh ( angle_to_radian * arg1 )
!
!  Trace (of a scalar)
!
  else if ( s_eqi ( sym, 'TRACE' ) ) then
    result = arg1
!
!  Transpose (of a scalar)
!
  else if ( s_eqi ( sym, 'TRANS' ) ) then
    result = arg1
!
!  Zero
!
  else if ( s_eqi ( sym, 'ZEROS' ) ) then
    result = 0.0E+00
!
!  Factorial function.
!
  else if ( sym == '!' ) then

    if ( arg1 < 0.0E+00 .or. arg1 > 25.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FUNSCL - Error!'
      write ( *, '(a)' ) '  Argument out of range for factorial!'
      ierror = 1
      return
    end if

    temp = arg1
    result = 1.0E+00

    do while ( temp > 1.0E+00 )
      result = result * temp
      temp = temp - 1.0E+00
    end do
!
!  Unknown function.
!
  else
    result = 0.0E+00
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FUNSCL - Error!'
    write ( *, '(a)' ) '  Unrecognized function name:' // trim ( sym )
  end if

  return
end
subroutine funval ( iarg1, iarg2, ierror, ifree, iopsym, iprsym, &
  ipset, ipval, irad, irtol, itemp, maxsym, maxval, nsyms, &
  sym, symbol, valsym )

!*******************************************************************************
!
!! FUNVAL evaluates a scalar, vector or matrix valued function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IARG1, IARG2, are the indices
!    of the arguments to the current function.
!
!    Output, integer IERROR.
!
!    IERROR is an error flag.  If it is zero on return,
!    then no error was detected, and the formula was
!    evaluated successfully.  If it is nonzero on return,
!    then an error was found, and the evaluation of the
!    formula could not be carried out.
!
!    Input, integer IFREE, the address of the next free
!    memory location in VALSYM.
!
!    Input, integer IOPSYM(MAXSYM), contains, for each
!    symbol, the number of operands.  In particular, if
!    a symbol represents a constant, IOPSYM(I) is 0.
!    If a symbol represents a unary operator such as ABS,
!    IOPSYM(I) is 1.  IOPSYM may be as large as 3.
!
!    Input, integer IPRSYM(MAXSYM), the priorities of the functions.
!
!    Output, integer IPSET.  If the function was an 'INDX1' or
!    'INDX2', then IPSET is the relative offset of the value.
!
!    Input, integer IPVAL(MAXSYM), contains, for each symbol,
!    the address in VALSYM where the value of the symbol
!    is stored, if it is a scalar.  If the symbol represents
!    a vector or matrix, IPVAL(IARG) points to location of
!    the first entry.
!
!    Input, integer ITEMP, the index of the constant.
!
!    Input, integer MAXIW, the amount of space in IWORK, which
!    should be at least MAXROW.
!
!    Input, integer MAXSYM, the maximum number of symbols allowed.
!
!    Input, integer MAXVAL, the maximum number of values allowed
!    in VALSYM.
!
!    Input, integer NSYMS, the number of declared symbols.
!
!    Input, integer NUMDIM(2,MAXSYM).
!    For each symbol I, NUMDIM(1,I) is the number of rows in its
!    value, and NUMDIM(2,I) is the number of columns.
!
!    Input, character ( len = MAXLEN ) SYM, the symbolic name of the function.
!
!    Input, character ( len = MAXLEN ) SYMBOL(MAXSYM), the symbolic names
!    of the internally defined constants and functions.
!
!    Input, real VALSYM(MAXVAL), contains the values of all the
!    symbolic variables.
!
  implicit none

  integer, parameter :: maxlen = 20
  integer maxsym
  integer maxval

  real angle_to_radian
  real arg1
  real arg2
  character ( len = 3 ) ctemp
  integer i
  integer iarg1
  integer iarg2
  integer ierror
  integer ifree
  integer index1
  integer index4
  integer indx
  integer info
  integer iopsym(maxsym)
  integer iprsym(maxsym)
  integer ipset
  integer ipval(maxsym)
  integer irad
  integer irtol
  integer itemp
  integer jseed
  integer nsyms
  real result
  real rtol
  logical s_eqi
  character ( len = maxlen ) sym
  character ( len = maxlen ) symbol(maxsym)
  real valsym(maxval)

  ipset = 0
  ierror = 0
!
!  If the operator is assignment, then the user may be trying to
!  reset the values of certain reserved variables, including
!  PI and EPS.
!
  if ( sym == '=' ) then

    if ( s_eqi ( symbol(iarg1), 'E' ) .or. &
         s_eqi ( symbol(iarg1), 'EPS' ) .or. &
         s_eqi ( symbol(iarg1), 'PI' )  ) then

      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FUNVAL - Error!'
      write ( *, '(a)' ) '  You may not change the value of ' // &
        trim ( symbol(iarg1) )

      return
    end if

  end if

  index1 = ipval(iarg1)

  if ( nsyms >= maxsym ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FUNVAL - Error!'
    write ( *, '(a)' ) '  Not enough free memory left!'
    write ( *, '(a)' ) '  The KILL or INIT commands may help!'
    return
  end if

  nsyms = nsyms + 1
  iopsym(nsyms) = 0
  iprsym(nsyms) = 10
  itemp = itemp + 1
  call i_to_s_zero ( itemp, ctemp )
  symbol(nsyms) = 'STK000'
  symbol(nsyms)(4:6) = ctemp
  ipval(nsyms) = ifree
  index4 = ipval(nsyms)

  arg1 = valsym(index1)
  if ( iarg2 == 0 ) then
    arg2 = 0.0E+00
  else
    arg2 = valsym(ipval(iarg2))
  end if

  angle_to_radian = valsym(irad)
  rtol = valsym(irtol)

  call funscl ( angle_to_radian, arg1, arg2, ierror, result, rtol, sym )

  valsym(index4) = result

  if ( sym == '=' ) then
    arg1 = arg2
    valsym(index1) = arg2
  end if

  return
end
subroutine fxy ( dir, ierror, infix, irpn, ival, maxfix, &
  maxrpn, val_max, ndraw, nval, xstart, xval, ystart, yval )
!
!***********************************************************************
!
!! FXY computes points on a implicit curve defined by F(X,Y) = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer IVAL(VAL_MAX).  IVAL(I) is 1 if point I
!    is to be connected to point I+1.
!
!    Input, integer VAL_MAX, the maximum number of points that can be 
!    used for plots.
!
!    Input/output, integer NVAL, the current number of plot points.
!
  implicit none

  integer maxfix
  integer maxrpn
  integer val_max

  real dfdx
  real dfdy
  real dir
  real h
  real hnew
  integer ierror
  integer ifrm
  character ( len = 80 ) infix
  integer ipc
  integer irpn(maxrpn)
  integer ival(val_max)
  integer jval
  integer ndraw
  integer nval
  real tdot
  real tnorm
  real tx
  real txnew
  real ty
  real tynew
  real x
  real xnew
  real xstart
  real xval(val_max)
  real y
  real ynew
  real ystart
  real yval(val_max)
!
  jval = 0
  h = 0.05
!
!  First find the starting point.
!
  ipc = 1

  do

    xnew = xstart
    ynew = ystart

    call newton ( ierror, infix, ipc, irpn, maxfix, maxrpn, xnew, ynew )

    if ( ierror == 0 ) then
      exit 
    end if

    if ( ipc == 2 ) then
      ierror = 1
      return
    end if

    ipc = 2

  end do

  if ( ipc == 1 ) then
    txnew = dir
    tynew = 0.0E+00
  else
    txnew = 0.0E+00
    tynew = dir
  end if
!
!  main loop
!
  do

    jval = jval+1
    nval = nval+1
    xval(nval) = xnew
    yval(nval) = ynew
    ival(nval) = 1

    if ( jval >= ndraw ) then
      ival(nval) = 0
      exit
    end if

    y = ynew
    x = xnew
    tx = txnew
    ty = tynew
!
!  Compute the tangent vector
!
    ifrm = 1
    call comrpn ( 'V', ierror, ifrm, infix, irpn, maxrpn, 'X', x )

    call comrpn ( 'V', ierror, ifrm, infix, irpn, maxrpn, 'Y', y )

    ifrm = 2
    call comrpn ( 'E', ierror, ifrm, infix, irpn, maxrpn, ' ', dfdx )

    ifrm = 3
    call comrpn ( 'E', ierror, ifrm, infix, irpn, maxrpn, ' ', dfdy )

    if ( abs ( dfdx ) > abs ( dfdy ) ) then
      tynew = 1.0E+00
      txnew = - dfdy / dfdx
      ipc = 2
    else
      txnew = 1.0E+00
      tynew = - dfdx / dfdy
      ipc = 1
    end if

    tnorm = sqrt ( txnew**2 + tynew**2 )
    txnew = txnew / tnorm
    tynew = tynew / tnorm
    tdot = txnew * tx + tynew * ty

    if ( tdot < 0.0E+00 ) then
      txnew = -txnew
      tynew = -tynew
    end if
!
!  Set the stepsize, compute the starting point, and call NEWTON.
!
    hnew = h

    do

      xnew = x + hnew * txnew
      ynew = y + hnew * tynew

      call newton ( ierror, infix, ipc, irpn, maxfix, maxrpn, xnew, ynew )

      if ( ierror == 0 ) then
        exit
      end if

      if ( hnew <= 0.125E+00 * h ) then
        ierror = 2
        return
      end if

      hnew = 0.5E+00 * hnew
 
    end do

  end do

  return
end
subroutine get_unit ( iunit )

!*******************************************************************************
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine getdat ( dat, idat, ierror, maxdat, mdat, namdat, ndat )

!***********************************************************************
!
!! GETDAT gets information from the user defining the data vectors.
!
!  Discussion:
!
!    The user is allowed to specify the values of the variables in
!    a variety of ways, including:
!
!    * equally spaced between given limits;
!    * Chebyshev points between given limits;
!    * to be typed in by the user;
!    * to be read in as columns of a disk file;
!    * to be evaluated by a formula based on data already entered.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real DAT(MAXDAT,NDAT), contains in DAT(I,J) the I-th
!    value of variable J.
!
!    Output, integer IERROR, an error flag.  If IERROR is nonzero, then
!    DAT could not be set up properly because of some error.
!
!    Input, integer MAXDAT, the maximum number of points at which a
!    variable value is defined.
!
!    Input, integer MDAT, the number of variables to be defined.
!
!    Input, character ( len = MDAT ) NAMDAT, a string containing the one
!    character names of the variables, in an order that corresponds to the
!    order of the columns in DAT to which their values will be
!    assigned.
!
!    Output, integer NDAT, the number of points at which the value
!    of each variable was defined.
!
  implicit none

  integer maxdat

  character ( len = 9 ) code
  real dat(maxdat,9)
  logical done
  character ( len = 80 ) filnam
  integer i
  integer ierror
  integer indx(9)
  integer idat(maxdat)
  integer j
  integer k
  integer l
  logical s_eqi
  character ( len = 80 ) line
  character ( len = 9 ) list
  integer mdat
  integer n
  character ( len = * ) namdat
  integer ndat
  integer ndat2
  integer ngrid
  real rval
  character ( len = 80 ) string
  real xmax
  real xmin
!
  if ( mdat  <=  0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GETDAT - Fatal error!'
    write ( *, '(a,i6)' ) '  Nonpositive input value of MDAT = ', mdat
    ierror = 1
    return
  end if
!
!  Say hello, and get data codes.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Data is needed for the following:'
  write ( *, '(1x,9(a1,1x))' ) ( namdat(i:i), i = 1, mdat )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' )'There is a code for each way of defining data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '*: Cartesian product;'
  write ( *, '(a)' ) 'E: Equally spaced data, user limits;'
  write ( *, '(a)' ) 'C: Chebyshev spaced data, user limits;'
  write ( *, '(a)' ) 'F: Formula in terms of other data;'
  write ( *, '(a)' ) 'D: Data is a column in a file;'
  write ( *, '(a)' ) 'T: Data will be typed in.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter a code for each data item:'
  write ( *, '(a)' ) '  (Example: "ETF means '
  write ( *, '(a)' ) '  X is equally spaced,'
  write ( *, '(a)' ) '  Y is typed in, and'
  write ( *, '(a)' ) '  Z is a formula in X and Y.)'

  read ( *, '(a)' ) code

  do i = 1, mdat

    if ( code(i:i) /= '*' .and. &
      .not. s_eqi ( code(i:i), 'E' ) .and. &
      .not. s_eqi ( code(i:i), 'C' ) .and. &
      .not. s_eqi ( code(i:i), 'F' ) .and. &
      .not. s_eqi ( code(i:i), 'D' ) .and. &
      .not. s_eqi ( code(i:i), 'T' ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GETDAT - Fatal error!'
      write ( *, '(a,i6,a)' ) '  Code letter ', i, ' is illegal.'
      ierror = 1
      return
    end if

  end do
!
!  The user may specify the number of data points implicitly or explicitly.
!  We have to be prepared to catch the specification when it comes,
!  and enforce it later on other data items.
!
  if ( index ( code(1:mdat), '*' ) == 0 ) then
    ndat = 0
  else
    ndat = 1
  end if
!
!  0: Handle any "*" codes.
!
  do i = 1, mdat

    if ( code(i:i) == '*' ) then

      string = 'Enter the number of grid points for ' // namdat(i:i)
      write ( *, '(a)' ) trim ( string )

      read ( *, * ) ngrid

      string = 'Enter lower and upper limits for ' // namdat(i:i)
      write ( *, '(a)' ) trim ( string )

      read ( *, * ) xmin, xmax

    else

      ngrid = 1
      xmin = 0.0E+00
      xmax = 0.0E+00

    end if

    do j = 1, ngrid

      if ( j == 1 ) then

        dat(1:ndat,i) = xmin

      else

        do l = 1, i-1
          do k = 1, ndat
            dat(k+(j-1)*ndat,l) = dat(k,l)
          end do
        end do

        call rvec_even_select ( xmin, xmax, ngrid, j, rval )

        do k = 1, ndat
          dat(k+(j-1)*ndat,i) = rval
        end do

      end if

    end do

    ndat = ndat * ngrid

  end do
!
!  1: Handle any "E" or "C" codes.
!
  do i = 1, mdat

    if ( s_eqi (code(i:i), 'E' ) .or. s_eqi (code(i:i), 'C' ) ) then

      if ( ndat == 0 ) then
        write ( *, '(a)' ) 'Enter the number of data points.'
        read ( *, * ) ndat
      end if

      string = 'Enter lower and upper limits for ' // namdat(i:i)
      write ( *, '(a)' ) trim ( string )

      read ( *, * ) xmin, xmax

      if ( s_eqi ( code(i:i), 'E' ) ) then
        call rvec_even ( xmin, xmax, ndat, dat(1,i) )
      else
        call r2_cheby ( ndat, xmin, xmax, dat(1,i) )
      end if

    end if

  end do
!
!  2: Handle any "T" codes:
!
  n = 0
  list = ' '
  do i = 1, mdat
    if ( s_eqi(code(i:i),'T') ) then
      n = n + 1
      list(2*n-1:2*n-1) = namdat(i:i)
      indx(n) = i
    end if
  end do

  if ( n > 0 ) then

    ndat2 = 0

    write ( *, '(a)' ) ' '
    string = 'Type in values of ' // list(1:2*n-1)
    write ( *, '(a)' ) trim ( string )

    if ( ndat == 0 ) then
      write ( *, '(a)' ) 'Terminate with a blank line.'
    else
      write ( *, '(a,i6)' )  '  The number of values needed is ', ndat
    end if

    j = 1

10      continue

    read ( *, '(a)' ) line

    if ( ndat == 0 .and. line == ' ' ) then
      ndat = ndat2
      write ( *, '(a,i6)' ) '  Number of values entered was ', ndat
      go to 30
    end if

    done = .true.

20      continue

    call r_next ( line, rval, done )

    if ( .not. done ) then

      dat ( ndat2+1, indx(j) ) = rval

      j = j + 1
      if ( j > n ) then
        j = 1
        ndat2 = ndat2 + 1
      end if

      if ( ndat /= 0 .and. ndat2 == ndat ) then
        go to 30
      end if
      go to 20
    end if

    go to 10

30      continue

  end if
!
!  3:  Handle disk file data.
!  If NDAT is still 0, read til end of file.
!  Otherwise, only read NDAT lines of file.
!
  n = 0
  list = ' '

  do i = 1, mdat
    if ( s_eqi ( code(i:i), 'D') ) then
      n = n + 1
      list(2*n-1:2*n-1) = namdat(i:i)
      indx(n) = i
    end if
  end do

  if ( n > 0 ) then

    write ( *, '(a)' ) ' '
    string = 'Enter file name containing data for ' // list(1:2*n-1)
    write ( *, '(a)' ) trim ( string )
    read ( *,'(a)' ) filnam
    open ( unit = 1, file = filnam, status = 'old' )

    ndat2 = 0


40      continue

    read ( 1, '(a)', end = 60 ) line
!
!  A line beginning with '#' is a comment.
!
    if ( line(1:1) == '#' ) then
      go to 40
    end if
!
!  A blank line means that the previously read point terminates a line.
!
    if ( line == ' ' ) then
      if ( ndat2 > 0 ) then
        idat(ndat2) = 0
      end if
      go to 40
    end if
!
!  Prepare to extract real values from the line of data.
!
    done = .true.
    j = 0

50      continue

    call r_next ( line, rval, done )

    if ( done ) then
      go to 40
    end if

    j = j+1
!
!  Store the value in DAT.
!
    dat(ndat2+1,indx(j)) = rval

    if ( j >= n ) then
      j = 0
      ndat2 = ndat2 + 1
      idat(ndat2) = 1
      if ( ndat /= 0 .and. ndat2 >= ndat ) then
        go to 70
      end if
    end if

    go to 50
!
!  On end of file, make sure we read enough data.
!
60      continue

    if ( ndat /= 0 ) then
      if ( ndat2 < ndat ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GETDAT2 - Fatal error!'
        write ( *, '(a)' ) '  Not enough data was in the file.'
        ierror = 1
        return
      end if
    else
      ndat = ndat2
      write ( *, '(a,i6)' ) '  Number of sets of data read was ', ndat
    end if
!
!  Jump here if we've successfully read NDAT sets of data.
!
70      continue

    close ( unit = 1 )

  else

    do i = 1, ndat
      idat(i) = 1
    end do

    idat(ndat) = 0

  end if
!
!  4: Ready to handle the "F" option?
!
  if ( ndat == 0 ) then
    ierror = 1
    return
  end if

  do i = 1, mdat

    if ( s_eqi ( code(i:i), 'F' ) ) then

      call evalf ( dat, ierror, i, maxdat, mdat, namdat, ndat )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) 'GETDAT - Fatal error!'
        return
      end if
    end if

  end do

  return
end
subroutine hello
!
!***********************************************************************
!
!! HELLO reports the program name, and other information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XYZ_PLOT'
  write ( *, '(a)' ) 'The interactive plotter'
  write ( *, '(a)' ) 'Version 1.06'
  write ( *, '(a)' ) 'Last modified on 25 April 1998'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '18 November 2000'
  write ( *, '(a)' ) '  Resurrected yet again, starting version 1.06.'
  write ( *, '(a)' ) '25 April 1998:'
  write ( *, '(a)' ) '  Added arbitrary orthographic plane projection.'
  write ( *, '(a)' ) '17 April 1998:'
  write ( *, '(a)' ) '  Finished plane projection setup.'
  write ( *, '(a)' ) '16 April 1998:'
  write ( *, '(a)' ) '  Using ARROW to draw vectors (untested)'
  write ( *, '(a)' ) '14 April 1998:'
  write ( *, '(a)' ) '  Added XYZ output file option.'
  write ( *, '(a)' ) '13 April 1998: X, Y, Z 3D projection options.'
  write ( *, '(a)' ) '31 March 1998: Sketched Cartesian products.'
  write ( *, '(a)' ) '31 March 1998: Began CON implementation.'
  write ( *, '(a)' ) '27 March 1998: Parameterized least squares.'
  write ( *, '(a)' ) '27 March 1998: Preliminary SCAT triangulation.'
  write ( *, '(a)' ) '26 March 1998: UVXY uses GETDAT.'
  write ( *, '(a)' ) '25 March 1998: Added scatterplot.'
  write ( *, '(a)' ) '24 March 1998: RT, TXY plots use GETDAT.'
  write ( *, '(a)' ) '19 March 1998: XYZ files allow # comments.'
  write ( *, '(a)' ) '19 March 1998: Purged CHRWRT, etc.'
  write ( *, '(a)' ) '18 March 1998: Revised XYZ call.'
  write ( *, '(a)' ) '17 March 1998: Updated PWP, XYS calls.'
  write ( *, '(a)' ) '16 March 1998: Updated FX, LSQ calls.'
  write ( *, '(a)' ) '12 March 1998: Renamed Z(X,Y) option.'

  return
end
function i_gcd ( i, j )
!
!*******************************************************************************
!
!! I_GCD finds the greatest common divisor of I and J.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, two numbers whose greatest common divisor
!    is desired.
!
!    Output, integer I_GCD, the greatest common divisor of I and J.
!
!    Note that only the absolute values of I and J are
!    considered, so that the result is always nonnegative.
!
!    If I or J is 0, I_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
!
!    If I and J have no common factor, I_GCD is returned as 1.
!
!    Otherwise, using the Euclidean algorithm, I_GCD is the
!    largest common factor of I and J.
!
  implicit none

  integer i
  integer i_gcd
  integer ip
  integer iq
  integer ir
  integer j

  i_gcd = 1
!
!  Return immediately if either I or J is zero.
!
  if ( i == 0 ) then
    i_gcd = max ( 1, abs ( j ) )
    return
  else if ( j == 0 ) then
    i_gcd = max ( 1, abs ( i ) )
    return
  end if
!
!  Set IP to the larger of I and J, IQ to the smaller.
!  This way, we can alter IP and IQ as we go.
!
  ip = max ( abs ( i ), abs ( j ) )
  iq = min ( abs ( i ), abs ( j ) )
!
!  Carry out the Euclidean algorithm.
!
  do

    ir = mod ( ip, iq )

    if ( ir == 0 ) then
      exit
    end if

    ip = iq
    iq = ir

  end do

  i_gcd = iq

  return
end
function i_lcm ( i, j )

!*******************************************************************************
!
!! I_LCM computes the least common multiple of two integers.
!
!  Discussion:
!
!    The least common multiple may be defined as
!
!      LCM(I,J) = ABS( I * J ) / GCF(I,J)
!
!    where GCF(I,J) is the greatest common factor of I and J.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the integers whose I_LCM is desired.
!
!    Output, integer I_LCM, the least common multiple of I and J.
!    I_LCM is never negative.  I_LCM is 0 if either I or J is zero.
!
  implicit none

  integer i
  integer i_gcd
  integer j
  integer i_lcm

  i_lcm = abs ( i * ( j / i_gcd ( i, j ) ) )

  return
end
subroutine i_swap ( i, j )

!*******************************************************************************
!
!! I_SWAP swaps two integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer i
  integer j
  integer k

  k = i
  i = j
  j = k

  return
end
subroutine i_to_s_zero ( intval, s )

!*******************************************************************************
!
!! I_TO_S_ZERO converts an integer to a string, with zero padding.
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  000001
!        -1  -00001
!         0  000000
!      1952  001952
!    123456  123456
!   1234567  ******  <-- Not enough room!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be right justified, and zero padded.
!    If there is not enough space, the string will be filled with stars.
!
  implicit none

  character c
  integer i
  integer idig
  integer ihi
  integer ilo
  integer intval
  integer ipos
  integer ival
  character ( len = * ) s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  Working from right to left, strip off the digits of the integer
!  and place them into S(ILO:IHI).
!
  ipos = ihi

  do while ( ival /= 0 .or. ipos == ihi )

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

  end do
!
!  Fill the empties with zeroes.
!
  do i = ilo, ipos
    s(i:i) = '0'
  end do

  return
end
subroutine inicom ( ifinis, ifree, indx1, indx2, ineg, infix, intsym, iopsym, &
  iprsym, ipval, irad, irpn, irtol, istack, maxrpn, maxsym, &
  maxval, nints, nsym, nsymp, nsyms, symbol, valsym, value )

!*******************************************************************************
!
!! INICOM initializes data for COMRPN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IFINIS, the index in SYMBOL of the symbol
!    '$', meaning the end of the formula.
!
!    Output, integer IFREE, the index of the first free address in VALSYM.
!
!    Output, integer INDX1, the index of the symbol 'INDEX1'.
!
!    Output, integer INDX2, the index of the symbol 'INDEX2'.
!
!    Output, integer INEG, the index of the symbol 'NEG'.
!
!    Output, character ( len = * ) INFIX, space for an infix formula.
!
!    Output, integer INTSYM(80), a set of integers which
!    are the indices of symbols, representing an infix formula.
!
!    Output, integer IPVAL(MAXSYM), contains, for each symbol,
!    the address in VALSYM where the value of the symbol
!    is stored, if it is a scalar.  If the symbol represents
!    a vector or matrix, IPVAL(IARG) points to location of
!    the first entry.
!
!    Output, integer IRAD, the index of the symbol 'ANGLE_TO_RADIAN'.
!
!    Output, integer IRPN(MAXRPN), used to store the compiled
!    versions of user formulas.
!
!    Workspace, integer ISTACK(MAXSYM), workspace for interpreting
!    the formula.
!
!    Input, integer MAXRPN, specifies the length of IRPN.
!
!    Input, integer MAXSYM, the maximum number of symbols allowed.
!
!    Input, integer MAXVAL, the maximum number of values allowed in VALSYM.
!
!    Output, integer NINTS, the number of useful entries in INTSYM.
!
!    Output, integer NSYM, the total number of symbols, including temporaries.
!
!    Output, integer NSYMP, the number of permanent symbols.
!
!    Output, integer NYMS, the number of declared symbols.
!
!    Output, integer NUMDIM(2,MAXSYM).
!    For each symbol I, NUMDIM(1,I) is the number of rows in its
!    value, and NUMDIM(2,I) is the number of columns.
!
!    Output, character ( len = MAXLEN ) SYMBOL(MAXSYM), the symbolic names
!    of the internally defined constants and functions.
!
!    Output, real VALSYM(MAXVAL), contains the values of all the
!    symbolic variables.
!
!    Workspace, real VALUE(MAXROW,MAXROW), space used to hold
!    the value of the variable to be saved.
!
  implicit none

  integer, parameter :: maxlen = 20
  integer maxrpn
  integer maxsym
  integer maxval

  integer i
  integer ifinis
  integer ifree
  integer indx1
  integer indx2
  integer ineg
  character ( len = * ) infix
  integer intsym(80)
  integer iopsym(maxsym)
  integer iprsym(maxsym)
  integer ipval(maxsym)
  integer irad
  integer irpn(maxrpn)
  integer irtol
  integer istack(maxsym)
  integer j
  integer nints
  integer nsym
  integer nsymp
  integer nsyms
  logical s_eqi
  character ( len = maxlen ) symbol(maxsym)
  real temp
  real valsym(maxval)
  real value

  nsym = 0

  nsym = nsym + 1
  symbol(nsym) = 'ABS'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'ACOS'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'ALOG'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'ALOG10'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'ANGLE_TO_RADIAN'
  iopsym(nsym) = 0
  iprsym(nsym) = 10
  valsym(nsym) = 1.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'ASIN'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'ATAN'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'ATAN2'
  iopsym(nsym) = 2
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'COS'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'COSH'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'COT'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'CSC'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'DET'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'DIAG'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'E'
  iopsym(nsym) = 0
  iprsym(nsym) = 10
  valsym(nsym) = exp ( 1.0E+00 )

  nsym = nsym + 1
  symbol(nsym) = 'EPS'
  iopsym(nsym) = 0
  iprsym(nsym) = 10
  valsym(nsym) = epsilon ( 1.0E+00 )

  nsym = nsym + 1
  symbol(nsym) = 'EVAL'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'EVEC'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'EXP'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'GCD'
  iopsym(nsym) = 2
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'HILBERT'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'HILBINV'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'HOUSE'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'ID'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'INT'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'INV'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'IVAL'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'KRYLOV'
  iopsym(nsym) = 3
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'LCM'
  iopsym(nsym) = 2
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'LN'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'LOG'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'LOG10'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'LOG2'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'MAX'
  iopsym(nsym) = 2
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'MIN'
  iopsym(nsym) = 2
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'NEG'
  iopsym(nsym) = 1
  iprsym(nsym) = 5
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'NINT'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'NORM0'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'NORM1'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'NORM2'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'NORME'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'NORMF'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'NORMS'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'PI'
  iopsym(nsym) = 0
  iprsym(nsym) = 10
  valsym(nsym) = 4.0E+00 * atan2 ( 1.0, 1.0E+00 )

  nsym = nsym + 1
  symbol(nsym) = 'POLY'
  iopsym(nsym) = 2
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'RAN'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'RCOND'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'ROUND'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'RTOL'
  iopsym(nsym) = 0
  iprsym(nsym) = 10
  valsym(nsym) = epsilon ( 1.0E+00 )

  nsym = nsym + 1
  symbol(nsym) = 'RVAL'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'SEC'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'SEED'
  iopsym(nsym) = 0
  iprsym(nsym) = 10
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'SIN'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'SINH'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'SQRT'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'STEP'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'TAN'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'TANH'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'TRACE'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'TRANS'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'ZEROS'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = '**'
  iopsym(nsym) = 2
  iprsym(nsym) = 8
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = '^'
  iopsym(nsym) = 2
  iprsym(nsym) = 8
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = '+'
  iopsym(nsym) = 2
  iprsym(nsym) = 5
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = '-'
  iopsym(nsym) = 2
  iprsym(nsym) = 5
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = '/'
  iopsym(nsym) = 2
  iprsym(nsym) = 7
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = '*'
  iopsym(nsym) = 2
  iprsym(nsym) = 6
  valsym(nsym) = 0.0E+00
!
!  Take care of fact that ninny UNIX FORTRAN compilers won't let us
!  explicitly use a backslash character!
!
  nsym = nsym + 1
  symbol(nsym) = char ( 92 )
  iopsym(nsym) = 2
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = '('
  iopsym(nsym) = -1
  iprsym(nsym) = 1
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = ')'
  iopsym(nsym) = -1
  iprsym(nsym) = 2
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = '='
  iopsym(nsym) = 2
  iprsym(nsym) = 3
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = ','
  iopsym(nsym) = -1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = '!'
  iopsym(nsym) = 1
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'INDEX1'
  iopsym(nsym) = 2
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = 'INDEX2'
  iopsym(nsym) = 3
  iprsym(nsym) = 9
  valsym(nsym) = 0.0E+00

  nsym = nsym + 1
  symbol(nsym) = '$'
  iopsym(nsym) = -1
  iprsym(nsym) = 0
  valsym(nsym) = 0.0E+00

  nsymp = nsym
  nsyms = nsymp

  ifree = nsymp + 1

  infix = ' '

  intsym(1:80) = 0
  irpn(1:maxrpn) = 0
  istack(1:maxsym) = 0

  nints = 0

  call ivec_identity ( nsymp, ipval )

  do i = 1, nsymp

    if ( s_eqi ( symbol(i), 'ANGLE_TO_RADIAN' ) ) then
      irad = i
    else if ( s_eqi ( symbol(i), 'NEG' ) ) then
      ineg = i
    else if ( symbol(i) == '$' ) then
      ifinis = i
    else if ( s_eqi ( symbol(i), 'INDEX1' ) ) then
      indx1 = i
    else if ( s_eqi ( symbol(i), 'INDEX2' ) ) then
      indx2 = i
    else if ( s_eqi ( symbol(i), 'RTOL' ) ) then
      irtol = i
    end if

  end do

  call random_seed ( )
 
  value = 0.0E+00

  return
end
subroutine ivec_identity ( n, a )
!
!*******************************************************************************
!
!! IVEC_IDENTITY sets an integer vector to the identity vector A(I)=I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, integer A(N), the array to be initialized.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer i
!
  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine least_set ( ntab, xtab, ytab, ndeg, ptab, array, eps, ierror )
!
!*******************************************************************************
!
!! LEAST_SET constructs the least squares polynomial approximation to data.
!
!  Discussion:
!
!    The routine LEAST_EVAL must be used to evaluate the approximation at a
!    point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Parameters:
!
!    Input, integer NTAB, the number of data points.
!
!    Input, real XTAB(NTAB), the X data.  The values in XTAB
!    should be distinct, and in increasing order.
!
!    Input, real YTAB(NTAB), the Y data values corresponding
!    to the X data in XTAB.
!
!    Input, integer NDEG, the degree of the polynomial which the
!    program is to use.  NDEG must be at least 1, and less than or
!    equal to NTAB-1.
!
!    Output, real PTAB(NTAB).  PTAB(I) is the value of the
!    least squares polynomial at the point XTAB(I).
!
!    Output, real ARRAY(2*NTAB+3*NDEG), an array containing data about
!    the polynomial.
!
!    Output, real EPS, the root-mean-square discrepancy of the
!    polynomial fit.
!
!    Output, integer IERROR, error flag.
!    zero, no error occurred;
!    nonzero, an error occurred, and the polynomial could not be computed.
!
  implicit none
!
  integer ndeg
  integer ntab
!
  real array(2*ntab+3*ndeg)
  real eps
  real error
  integer i
  integer i0l1
  integer i1l1
  integer ierror
  integer it
  integer k
  integer mdeg
  real ptab(ntab)
  real rn0
  real rn1
  real s
  real sum2
  real xtab(ntab)
  real y_sum
  real ytab(ntab)
!
  ierror = 0
!
!  Check NDEG.
!
  if ( ndeg < 1 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEAST_SET - Fatal error!'
    write ( *, '(a)' ) '  NDEG < 1.'
    return
  else if ( ndeg >= ntab ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEAST_SET - Fatal error!'
    write ( *, '(a)' ) '  NDEG >= NTAB.'
    return
  end if
!
!  Check that the abscissas are strictly increasing.
!
  do i = 1, ntab-1
    if ( xtab(i) >= xtab(i+1) ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEAST_SET - Fatal error!'
      write ( *, '(a)' ) '  XTAB must be strictly increasing, but'
      write ( *, '(a,i6,a,g14.6)' ) '  XTAB(', i, ') = ', xtab(i)
      write ( *, '(a,i6,a,g14.6)' ) '  XTAB(', i+1, ') = ', xtab(i+1)
      return
    end if
  end do

  i0l1 = 3 * ndeg
  i1l1 = 3 * ndeg + ntab

  y_sum = sum ( ytab )
  rn0 = ntab
  array(2*ndeg) = y_sum / real ( ntab )

  ptab(1:ntab) = y_sum / real ( ntab )

  error = 0.0E+00
  do i = 1, ntab
    error = error + ( y_sum / real ( ntab ) - ytab(i) )**2
  end do

  if ( ndeg == 0 ) then
    eps = sqrt ( error / real ( ntab ) )
    return
  end if

  array(1) = sum ( xtab ) / real ( ntab )

  s = 0.0E+00
  sum2 = 0.0E+00
  do i = 1, ntab
    array(i1l1+i) = xtab(i) - array(1)
    s = s + array(i1l1+i)**2
    sum2 = sum2 + array(i1l1+i) * ( ytab(i) - ptab(i) )
  end do

  rn1 = s
  array(2*ndeg+1) = sum2 / s

  do i = 1, ntab
    ptab(i) = ptab(i) + sum2 * array(i1l1+i) / s
  end do

  error = 0.0E+00
  do i = 1, ntab
    error = error + ( ptab(i) - ytab(i) )**2
  end do

  if ( ndeg == 1 ) then
    eps = sqrt ( error / real ( ntab ) )
    return
  end if

  do i = 1, ntab
    array(3*ndeg+i) = 1.0E+00
  end do

  mdeg = 2
  k = 2

  do

    array(ndeg-1+k) = rn1 / rn0

    sum2 = 0.0E+00
    do i = 1, ntab
      sum2 = sum2 + xtab(i) * array(i1l1+i)**2
    end do

    array(k) = sum2 / rn1

    s = 0.0E+00
    sum2 = 0.0E+00
    do i = 1, ntab
      array(i0l1+i) = ( xtab(i) - array(k) ) * array(i1l1+i) &
        - array(ndeg-1+k) * array(i0l1+i)
      s = s + array(i0l1+i)**2
      sum2 = sum2 + array(i0l1+i) * ( ytab(i) - ptab(i) )
    end do

    rn0 = rn1
    rn1 = s
    it = i0l1
    i0l1 = i1l1
    i1l1 = it
    array(2*ndeg+k) = sum2 / rn1

    do i = 1, ntab
      ptab(i) = ptab(i) + sum2 * array(i1l1+i) / rn1
    end do

    error = 0.0E+00
    do i = 1, ntab
      error = error + ( ptab(i) - ytab(i) )**2
    end do

    if ( mdeg >= ndeg ) then
      exit
    end if

    mdeg = mdeg + 1
    k = k + 1

  end do

  eps = sqrt ( error / real ( ntab ) )

  return
end
subroutine least_val ( x, ndeg, array, value )
!
!*******************************************************************************
!
!! LEAST_VAL evaluates a least squares polynomial defined by LEAST_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Parameters:
!
!    Input, real X, the point at which the polynomial is to be evaluated.
!
!    Input, integer NDEG, the degree of the polynomial fit used.
!    This is the value of NDEG as returned from LEAST_SET.
!
!    Input, real ARRAY(*), an array of a certain dimension.
!    See LEAST_SET for details on the size of ARRAY.
!
!    ARRAY contains information about the polynomial, as set up by LEAST_SET.
!
!    Output, real VALUE, the value of the polynomial at X.
!
  implicit none
!
  real array(*)
  real dk
  real dkp1
  real dkp2
  integer k
  integer l
  integer ndeg
  real value
  real x
!
  if ( ndeg <= 0 ) then

    value = array(2*ndeg)

  else if ( ndeg == 1 ) then

    value = array(2*ndeg) + array(2*ndeg+1) * ( x - array(1) )

  else

    dkp2 = array(3*ndeg)
    dkp1 = array(3*ndeg-1) + ( x - array(ndeg) ) * dkp2

    do l = 1, ndeg-2

      k = ndeg - 1 - l

      dk = array(2*ndeg+k) + ( x - array(k+1) ) * dkp1 - array(ndeg+1+k) * dkp2

      dkp2 = dkp1

      dkp1 = dk

    end do

    value = array(2*ndeg) + ( x - array(1) ) * dkp1 - array(ndeg+1) * dkp2

  end if

  return
end
subroutine line_seg_contains_point_1d ( x1, x2, x3, u )
!
!*******************************************************************************
!
!! LINE_SEG_CONTAINS_POINT_1D reports if a line segment contains a point in 1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, X2, two points defining a line segment.
!    The line segment has origin at X1, and unit at X2.
!
!    Input, real X3, a point to be tested.
!
!    Output, real U, the coordinate of X3 in units of (X2-X1).
!    The point X3 is contained in the line segment if 0 <= U <= 1.
!
  implicit none
!
  real u
  real unit
  real x1
  real x2
  real x3
!
  unit = x2 - x1

  if ( unit == 0.0E+00 ) then

    if ( x3 == x1 ) then
      u = 0.5E+00
    else if ( x3 < x1 ) then
      u = - huge ( u )
    else if ( x3 > x1 ) then
      u = huge ( u )
    end if

  else

    u = ( x3 - x1 ) / unit

  end if

  return
end
function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )
!
!***********************************************************************
!
!! LRLINE determines where a point lies in relation to a directed line.
!
!  Discussion:
!
!    LRLINE determines whether a point is to the left of, right of,
!    or on a directed line parallel to a line through given points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Parameters:
!
!    Input, XU,YU, XV1,YV1, XV2,YV2 - vertex coordinates; the directed
!              line is parallel to and at signed distance DV to the
!              left of the directed line from (XV1,YV1) to (XV2,YV2);
!              (XU,YU) is the vertex for which the position
!              relative to the directed line is to be determined
!
!    Input, DV - signed distance (positive for left)
!
!    Output, LRLINE - +1, 0, or -1 depending on whether (XU,YU) is
!              to the right of, on, or left of the directed line
!              (0 if line degenerates to a point)
!
  implicit none
!
  real, parameter :: tol = 0.0000001
!
  real dv
  real dx
  real dxu
  real dy
  real dyu
  integer lrline
  real t
  real temp
  real tolabs
  real xu
  real xv1
  real xv2
  real yu
  real yv1
  real yv2
!
  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( abs(dx), abs(dy), abs(dxu), abs(dyu), abs(dv) )

  t = dy * dxu - dx * dyu

  if ( dv /= 0.0E+00 ) then
    t = t + dv * sqrt ( dx**2 + dy**2 )
  end if

  temp = 1.0E+00
  lrline = int ( sign(temp,t) )

  if ( abs(t) <= tolabs ) then
    lrline = 0
  end if

  return
end
subroutine lsqtxy ( ierror, ival, maxrwork, val_max, ndat, nfine, &
  nval, rwork, tdat, xdat, xval, ydat, yval )
!
!***********************************************************************
!
!! LSQTXY sets up a parameterized least squares polynomial plot.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
!    Input/output, integer IVAL(VAL_MAX).  IVAL(I) is 1 if point I
!    is to be connected to point I+1.
!
!    Input, integer VAL_MAX, the maximum number of points that can be 
!    used for plots.
!
!    Input/output, integer NVAL, the current number of plot points.
!
  implicit none
!
  integer val_max
  integer ndat
  integer maxrwork
!
  real eps
  integer i
  integer ierror
  integer ival(val_max)
  integer iw1
  logical s_eqi
  integer ndeg
  integer nfine
  integer nval
  real rwork(maxrwork)
  real tdat(ndat)
  real thi
  real tlo
  real tval
  real xdat(ndat)
  real xval(val_max)
  real ydat(ndat)
  real yval(val_max)
!
  ierror = 0
  nval = nfine
!
!  Get NDEG, the degree of the polynomial.
!
  write ( *, '(a)' ) 'Enter the degree of the polynomial:'
  read ( *, * ) ndeg

  if ( ndeg < 0 ) then
    write ( *, '(a)' ) 'LSQTXY - Fatal error!'
    write ( *, '(a)' ) '  The polynomial degree must be nonnegative.'
    ierror = 1
    return
  end if
!
!  Make sure we have enough work space.
!
  if ( 3*ndat+3*ndeg > maxrwork ) then
    write ( *, '(a)' ) 'LSQTXY - Fatal error!'
    write ( *, '(a)' ) '  There is not enough work space to compute'
    write ( *, '(a)' ) '  the least squares polynomial.'
    ierror = 1
    return
  end if
!
!  Set up the least squares polynomial for X.
!
  iw1 = 2*ndat+3*ndeg+1
  call least_set ( ndat, tdat, xdat, ndeg, rwork(iw1), rwork, eps, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Fill in the graph values, XVAL.
!
  tlo = minval ( tdat(1:ndat) )
  thi = maxval ( tdat(1:ndat) )

  do i = 1, nval
    call rvec_even_select ( tlo, thi, nval, i, tval )
    call least_val ( tval, ndeg, rwork, xval(i) )
    ival(i) = 1
  end do

  call least_set ( ndat, tdat, ydat, ndeg, rwork(iw1), rwork, eps, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Fill in the graph values, YVAL.
!
  do i = 1, nval
    call rvec_even_select ( tlo, thi, nval, i, tval )
    call least_val ( tval, ndeg, rwork, yval(i) )
    ival(i) = 1
  end do

  ival(nval) = 0

  return
end
subroutine lsqxy ( ierror, ival, maxrwork, val_max, ndat, nfine, &
  nval, rwork, xdat, xval, ydat, yval )
!
!***********************************************************************
!
!! LSQXY sets up the data defining a least squares polynomial plot.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
!    Input/output, integer IVAL(VAL_MAX).  IVAL(I) is 1 if point I
!    is to be connected to point I+1.
!
!    Input, integer VAL_MAX, the maximum number of points that can be 
!    used for plots.
!
!    Input/output, integer NVAL, the current number of plot points.
!
  implicit none
!
  integer val_max
  integer ndat
  integer maxrwork
!
  real eps
  integer i
  integer ierror
  integer ival(val_max)
  integer iw1
  logical s_eqi
  integer ndeg
  integer nfine
  integer nval
  real rwork(maxrwork)
  real xdat(ndat)
  real xhi
  real xlo
  real xval(val_max)
  real ydat(ndat)
  real yval(val_max)
!
!  Get NDEG, the degree of the polynomial.
!
  write ( *, '(a)' ) 'Enter the degree of the polynomial:'
  read ( *, * ) ndeg

  if ( ndeg < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LSQ - Fatal error!'
    write ( *, '(a)' ) '  The polynomial degree must be nonnegative.'
    ierror = 1
    return
  end if
!
!  Make sure we have enough work space.
!
  if ( 3 * ndat + 3 * ndeg > maxrwork ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LSQ - Fatal error!'
    write ( *, '(a)' ) '  There is not enough work space to compute'
    write ( *, '(a)' ) '  the least squares polynomial.'
    ierror = 1
    return
  end if
!
!  Set up the least squares polynomial.
!
  iw1 = 2 * ndat + 3 * ndeg + 1
  call least_set ( ndat, xdat, ydat, ndeg, rwork(iw1), rwork, eps, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Fill in the graph values, XVAL, YVAL.
!
  xlo = minval ( xdat(1:ndat) )
  xhi = maxval ( xdat(1:ndat) )

  do i = 1, nfine

    nval = nval+1
    call rvec_even_select ( xlo, xhi, nfine, i, xval(nval) )
    call least_val ( xval(nval), ndeg, rwork, yval(nval) )
    ival(nval) = 1

  end do

  ival(nval) = 0

  return
end
subroutine newton ( ierror, infix, ipc, irpn, maxfix, maxrpn, x, y )
!
!***********************************************************************
!
!! NEWTON applies Newton's method to solve F(X,Y)=0, with X or Y fixed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none
!
  integer maxfix
  integer maxrpn
!
  real dfdx
  real dfdy
  real fxy
  integer ierror
  integer ifrm
  character ( len = 80 ) infix
  integer ipc
  integer irpn(maxrpn)
  integer istep
  integer, parameter :: maxstp = 20
  real x
  real y
!
  ierror = 0
  istep = 0

  do
!
!  Evaluate F(X,Y).
!
    ifrm = 1
    call comrpn ( 'V', ierror, ifrm, infix, irpn, maxrpn, 'X', x )

    call comrpn ( 'V', ierror, ifrm, infix, irpn, maxrpn, 'Y', y )

    ifrm = 1
    call comrpn ( 'E', ierror, ifrm, infix, irpn, maxrpn, ' ', fxy )

    if ( abs ( fxy ) < 0.00001 ) then
      exit
    end if
!
!  Increment the step counter.
!
    istep = istep + 1

    if ( istep > maxstp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NEWTON - Fatal error!'
      write ( *, '(a)' ) '  Too many Newton steps!'
      ierror = 1
      exit
    end if
!
!  Apply the Newton correction.
!
    if ( ipc == 2 ) then

      ifrm = 2
      call comrpn ( 'E', ierror, ifrm, infix, irpn, maxrpn, ' ', dfdx )

      x = x - fxy/dfdx

    else

      ifrm = 3
      call comrpn ( 'E', ierror, ifrm, infix, irpn, maxrpn, ' ', dfdy )

      y = y - fxy/dfdy

    end if

  end do

  return
end
subroutine option ( dev, iplot, ldat, lframe, ncon, &
  nfine, nxgrid, nygrid, theta, title, vscale )
!
!***********************************************************************
!
!! OPTION allows the user to set various options.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none
!
  character ( len = * ) dev
  integer icom
  integer ierror
  integer ios
  integer iplot
  character ( len = 80 ) isay
  integer itemp
  logical ldat
  logical lframe
  integer ncon
  integer nfine
  integer nxgrid
  integer nygrid
  logical s_eqi
  real theta
  character ( len = 80 ) title
  real vscale
!
  do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter option to change, or H for help, or RETURN if done:'

    read ( *, '(a)', iostat = ios ) isay

    if ( isay == ' ' ) then
      exit
    end if
!
!  C means cancel the current plot.
!
    if ( s_eqi ( isay(1:1), 'C' ) ) then

      exit
!
!  'DATA/NODATA' means show or don't show the data points.
!
    else if ( s_eqi ( isay(1:4), 'data' ) ) then

      ldat = .true.
      write ( *, '(a)' ) 'Data points WILL be shown.'

    else if ( s_eqi ( isay(1:6), 'nodata' ) ) then

      ldat = .false.
      write ( *, '(a)' ) 'Data points will NOT be shown.'
!
!  FRAME/NOFRAME
!
    else if ( s_eqi ( isay(1:1), 'F' ) ) then

      lframe = .true.
      write ( *, '(a)' ) 'A plot frame will be drawn.'

    else if ( s_eqi ( isay(1:3), 'NOF' ) ) then

      lframe = .false.
      write ( *, '(a)' ) 'A plot frame will NOT be drawn.'
!
!  G: Draw the graph.
!
    else if ( s_eqi ( isay(1:1), 'G') ) then

      exit
!
!  Help: Print the list of commands.
!
    else if ( s_eqi ( isay(1:1), 'H' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Your options are:'
      write ( *, '(a)' ) ' '

      if ( iplot > 0 ) then
        write ( *, '(a)' ) 'C  cancel this plot'
      end if

      if ( ldat ) then
        write ( *, '(a)' ) 'NODATA   Do not show the data points.'
      else
        write ( *, '(a)' ) 'DATA     Show the data points.'
      end if

      write ( *, '(a)' ) 'DEV=     Set the graphics output device.'
 
      if ( lframe ) then
        write ( *, '(a)' ) 'NOFRAME  Do NOT frame the plot.'
      else
        write ( *, '(a)' ) 'FRAME    Frame the plot.'
      end if

      write ( *, '(a)' ) 'G        Ready to plot'
      write ( *, '(a)' ) 'H        Help (list these options).'
      write ( *, '(a)' ) 'NCON=    Set number of contour levels.'
      write ( *, '(a)' ) 'NFINE=   Set number of points on curves.'
      write ( *, '(a)' ) 'NXGRID=  Set number of X grid lines.'
      write ( *, '(a)' ) 'NYGRID=  Set number of Y grid lines.'
      write ( *, '(a)' ) 'Q        Quit, stop this program, now!'
      write ( *, '(a)' ) 'THETA=   Set the 3d presentation angle.'
      write ( *, '(a)' ) 'TITLE=   Set the title.'
      write ( *, '(a)' ) 'VSCALE=  Set vector scale size.'
!
!  NCON =
!
    else if ( s_eqi(isay(1:5),'ncon=' ) ) then

      read ( isay(6:), * ) itemp
      ncon = itemp
      write ( *, '(a,i6)' ) 'Number of contour levels set to ', ncon
!
!  NFINE =
!
    else if ( s_eqi(isay(1:6),'nfine=') ) then

      read ( isay(7:), * ) itemp
      nfine = itemp
      write ( *, '(a,i6)' ) 'Number of points on curves set to ', nfine
!
!  NXGRID =
!
    else if ( s_eqi(isay(1:7),'nxgrid=') ) then

      read ( isay(8:), * ) nxgrid
      write ( *, '(a,i6)' ) 'Number of X grid lines set to ', nxgrid
!
!  NYGRID =
!
    else if ( s_eqi(isay(1:7),'nygrid=') ) then
 
      read ( isay(8:), * ) nygrid
      write ( *, '(a,i6)' ) 'Number of Y grid lines set to ', nygrid
!
!  Q: Quit
!
    else if ( s_eqi(isay(1:1),'Q' ) ) then

      write ( *, '(a)' ) 'Enter "Y" so I know you really mean it!'
      read ( *, '(a1)' ) isay

      if ( s_eqi(isay,'y') ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'XYZ_PLOT is quitting now.'

        icom = 1
        call anyplt ( icom )

        if ( .not.s_eqi(dev,'cgmb') ) then
          call file_delete ( 'cgmout' )
        end if

        stop

      end if
!
!  THETA =
!
    else if ( s_eqi(isay(1:6),'theta=') ) then

      read ( isay(7:), * ) theta
!
!  TITLE =
!
    else if ( s_eqi(isay(1:6),'title=') ) then

      title = isay(7:)
!
!  VSCALE =
!
    else if ( s_eqi(isay(1:7),'vscale=') ) then

      read ( isay(8:), * ) vscale
!
!  Blank.
!
    else if ( s_eqi(isay,' ') ) then

    else

      write ( *, '(a)' ) 'That option did not make sense!'

    end if

  end do

  return
end
function r_pi ( )
!
!*******************************************************************************
!
!! R_PI returns the value of pi.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real R_PI, the value of pi.
!
  implicit none
!
  real r_pi
!
  r_pi = 3.14159265358979323846264338327950288419716939937510E+00

  return
end
subroutine plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
!
!*******************************************************************************
!
!! PLANE_EXP2IMP_3D converts an explicit plane to implicit form in 3D.
!
!  Definition:
!
!    The explicit form of a plane in 3D is
!
!      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
!
!    The implicit form of a plane in 3D is
!
!      A * X + B * Y + C * Z + D = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, X2, X3, Y3, Z3, are three points
!    on the plane, which must be distinct, and not collinear.
!
!    Output, real A, B, C, D, coefficients which describe the plane.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3

  a = ( y2 - y1 ) * ( z3 - z1 ) - ( z2 - z1 ) * ( y3 - y1 )
  b = ( z2 - z1 ) * ( x3 - x1 ) - ( x2 - x1 ) * ( z3 - z1 )
  c = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 )
  d = - x2 * a - y2 * b - z2 * c

  return
end
subroutine plane_exp_project_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  xf, yf, zf, npnt, xo, yo, zo, xp, yp, zp, ivis )
!
!*******************************************************************************
!
!! PLANE_EXP_PROJECT_3D projects points through a point onto a plane in 3D.
!
!  Discussion:
!
!    The explicit form of a plane in 3D is
!
!      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the
!    coordinates of three points on the plane.
!
!    Input, real XF, YF, ZF, are the coordinates of the focus point.
!
!    Input, integer NPNT, the number of points to project.
!
!    Input, real XO(NPNT), YO(NPNT), ZO(NPNT), are the coordinates of
!    the object points.
!
!    Output, real XP(NPNT), YP(NPNT), ZP(NPNT), are the coordinates of the
!    projections of the object points through the focus point onto
!    the plane.  XP, YP, and ZP may share the same memory as XO, YO,
!    and ZO, in which case the projections will overwrite the original data.
!
!    Output, integer IVIS(NPNT), visibility indicator:
!    3, the object was behind the plane;
!    2, the object was already on the plane;
!    1, the object was between the focus and the plane;
!    0, the line from the object to the focus is parallel to the plane,
!    so the object is "invisible".
!    -1, the focus is between the object and the plane.  The object
!    might be considered invisible.
!
  implicit none
!
  integer npnt
!
  real a
  real alpha
  real angle_rad_3d
  real b
  real beta
  real c
  real d
  real disfo
  real disfn
  integer i
  integer ivis(npnt)
  real x1
  real x2
  real x3
  real xf
  real xn
  real xo(npnt)
  real xp(npnt)
  real y1
  real y2
  real y3
  real yf
  real yn
  real yo(npnt)
  real yp(npnt)
  real z1
  real z2
  real z3
  real zf
  real zn
  real zo(npnt)
  real zp(npnt)
!
!  Put the plane into ABCD form.
!
  call plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
!
!  Get the nearest point on the plane to the focus.
!
  call plane_imp_point_near_3d ( a, b, c, d, xf, yf, zf, xn, yn, zn )
!
!  Get the distance from the focus to the plane.
!
  call points_dist_3d ( xf, yf, zf, xn, yn, zn, disfn )
!
!  If the focus lies in the plane, this is bad.  We could still
!  project points that actually lie in the plane, but we'll
!  just bail out.
!
  if ( disfn == 0.0E+00 ) then
    ivis(1:npnt) = 0
    xp(1:npnt) = xf
    yp(1:npnt) = yf
    zp(1:npnt) = zf
    return
  end if
!
!  Process the points.
!
  do i = 1, npnt
!
!  Get the distance from the focus to the object.
!
    call points_dist_3d ( xf, yf, zf, xo(i), yo(i), zo(i), disfo )

    if ( disfo == 0.0E+00 ) then

      ivis(i) = 0
      xp(i) = xn
      yp(i) = yn
      zp(i) = zn

    else
!
!  Compute ALPHA, the angle between (OBJECT-FOCUS) and (NEAREST-FOCUS).
!
      alpha = angle_rad_3d ( xo(i), yo(i), zo(i), xf, yf, zf, xn, yn, zn )

      if ( cos ( alpha ) == 0.0E+00 ) then

        ivis(i) = 0
        xp(i) = xn
        yp(i) = yn
        zp(i) = zn

      else
!
!  BETA is Dist(NEAREST-FOCUS) / ( Cos(ALPHA)*Dist(OBJECT-FOCUS) )
!
        beta = disfn / ( cos ( alpha ) * disfo )

        if ( beta > 1.0E+00 ) then
          ivis(i) = 1
        else if ( beta == 1.0E+00 ) then
          ivis(i) = 2
        else if ( beta > 0.0E+00 ) then
          ivis(i) = 3
        else
          ivis(i) = -1
        end if
!
!  Set the projected point.
!
        xp(i) = xf + beta * ( xo(i) - xf )
        yp(i) = yf + beta * ( yo(i) - yf )
        zp(i) = zf + beta * ( zo(i) - zf )

      end if

    end if

  end do

  return
end
subroutine plane_imp_point_near_3d ( a, b, c, d, x, y, z, xn, yn, zn )
!
!*******************************************************************************
!
!! PLANE_IMP_POINT_NEAR_3D: nearest point on a implicit plane to a point in 3D.
!
!  Definition:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, coefficients that define the plane as
!    the set of points for which A*X+B*Y+C*Z+D = 0.
!
!    Input, real X, Y, Z, the coordinates of the point.
!
!    Output, real XN, YN, ZN, the coordinates of the nearest point on
!    the plane.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real t
  real x
  real xn
  real y
  real yn
  real z
  real zn
!
  if ( a == 0.0E+00 .and. b == 0.0E+00 .and. c == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_POINT_NEAR_3D - Fatal error!'
    write ( *, '(a)' ) '  A = B = C = 0.'
    stop
  end if
!
!  The normal N to the plane is (A,B,C).
!
!  The line defined by (XN-X)/A = (YN-Y)/B = (ZN-Z)/C = T
!  goes through (X,Y,Z) and is parallel to N.
!
!  Solving for the point (XN,YN,ZN) we get
!
!    XN = A*T+X
!    YN = B*T+Y
!    ZN = C*T+Z
!
!  Now place these values in the equation for the plane:
!
!    A*(A*T+X) + B*(B*T+Y) + C*(C*T+Z) + D = 0
!
!  and solve for T:
!
!    T = (-A*X-B*Y-C*Z-D) / (A * A + B * B + C * C )
!
  t = - ( a * x + b * y + c * z + d ) / ( a * a + b * b + c * c )

  xn = x + a * t
  yn = y + b * t
  zn = z + c * t

  return
end
subroutine points_dist_3d ( x1, y1, z1, x2, y2, z2, dist )
!
!*******************************************************************************
!
!! POINTS_DIST_3D finds the distance between two points in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, determines the pair of points
!    (X1,Y1,Z1) and (X2,Y2,Z2) whose distance apart is be determined.
!
!    Output, real DIST, the distance between the points.
!
  implicit none
!
  real dist
  real x1
  real x2
  real y1
  real y2
  real z1
  real z2
!
  dist = sqrt ( ( x1 - x2 )**2 + ( y1 - y2 )**2 + ( z1 - z2 )**2 )

  return
end
subroutine proplane2 ( x1, y1, z1, x2, y2, z2, x3, y3, z3, npnt, xp, yp, zp, &
  alpha, beta )
!
!*******************************************************************************
!
!! PROPLANE2 produces 2D coordinates of points that lie in a plane, in 3D.
!
!  Discussion:
!
!    The plane is specified by three non-colinear points, which we will
!    call P1, P2 and P3.
!
!    The first thing to do is to compute two orthonormal vectors V1 and
!    V2, so that any point P that lies in the plane may be written as
!
!      P = P1 + alpha * V1 + beta * V2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the
!    coordinates of three points on the plane.  These three points
!    should not lie in a straight line, but this condition is not
!    checked.
!
!    Input, integer NPNT, the number of points to project.
!
!    Input, real XP(NPNT), YP(NPNT), ZP(NPNT), are the Cartesian
!    coordinates of points which lie on the plane spanned by the
!    three points.  These points are not checked to ensure that
!    they lie on the plane.
!
!    Output, real ALPHA(NPNT), BETA(NPNT), the "in-plane" coordinates of
!    the points.
!
  implicit none
!
  integer npnt
!
  real alpha(npnt)
  real beta(npnt)
  real dot
  integer i
  real v1(3)
  real v2(3)
  real x1
  real x2
  real x3
  real xp(npnt)
  real y1
  real y2
  real y3
  real yp(npnt)
  real z1
  real z2
  real z3
  real zp(npnt)
!
!  Compute the two basis vectors for the affine plane.
!
  v1(1) = x2 - x1
  v1(2) = y2 - y1
  v1(3) = z2 - z1

  call vector_unit_nd ( 3, v1 )

  v2(1) = x3 - x1
  v2(2) = y3 - y1
  v2(3) = z3 - z1

  dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)

  v2(1:3) = v2(1:3) - dot * v1(1:3)

  call vector_unit_nd ( 3, v2 )
!
!  Now decompose each (X,Y,Z).
!
  do i = 1, npnt

    alpha(i) = ( xp(i) - x1 ) * v1(1) +  ( yp(i) - y1 ) * v1(2) + &
               ( zp(i) - z1 ) * v1(3)

    beta(i) =  ( xp(i) - x1 ) * v2(1) + ( yp(i) - y1 ) * v2(2) + &
               ( zp(i) - z1 ) * v2(3)

  end do

  return
end
subroutine proplane3 ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  npnt, xo, yo, zo, xp, yp, zp )
!
!*******************************************************************************
!
!! PROPLANE3 projects points orthographically onto a plane, in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the
!    coordinates of three points on the plane.
!
!    Input, integer NPNT, the number of points to project.
!
!    Input, real XO(NPNT), YO(NPNT), ZO(NPNT), are the coordinates of
!    the object points.
!
!    Output, real XP(NPNT), YP(NPNT), ZP(NPNT), are the coordinates of
!    the projections of the object points through the focus point onto
!    the plane.
!
!    XP, YP, and ZP may share the same memory as XO, YO, and ZO, in
!    which case the projections will overwrite the original data.
!
  implicit none
!
  integer npnt
!
  real a
  real b
  real c
  real d
  integer i
  real x1
  real x2
  real x3
  real xo(npnt)
  real xp(npnt)
  real y1
  real y2
  real y3
  real yo(npnt)
  real yp(npnt)
  real z1
  real z2
  real z3
  real zo(npnt)
  real zp(npnt)
!
!  Put the plane into ABCD form.
!
  call plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
!
!  For each point, its image in the plane is the nearest point
!  in the plane.
!
  do i = 1, npnt

    call plane_imp_point_near_3d ( a, b, c, d, xo(i), yo(i), &
      zo(i), xp(i), yp(i), zp(i) )

  end do

  return
end
subroutine pwp ( ierror, ival, maxrwork, val_max, ndat, nfine, &
  nval, rwork, xdat, xval, ydat, yval )
!
!***********************************************************************
!
!! PWP sets up the data defining a piecewise polynomial plot.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
!    Input/output, integer IVAL(VAL_MAX).  IVAL(I) is 1 if point I
!    is to be connected to point I+1.
!
!    Input, integer VAL_MAX, the maximum number of points that can be 
!    used for plots.
!
!    Input/output, integer NVAL, the current number of plot points.
!
  implicit none
!
  integer val_max
  integer ndat
  integer maxrwork
!
  integer i
  integer ierror
  integer ihi
  integer ilo
  integer ival(val_max)
  integer j
  integer jhi
  logical s_eqi
  integer ndeg
  integer nfine
  integer nint
  integer nval
  real rwork(maxrwork)
  real xdat(ndat)
  real xval(val_max)
  real ydat(ndat)
  real yval(val_max)
!
!  Get NDEG, the polynomial degree of the pieces.
!
  ierror = 0

  write ( *, '(a)' ) 'Enter the polynomial degree of the pieces.'
  read ( *, * ) ndeg

  if ( ndeg < 0 ) then
    write ( *, '(a)' ) 'PWP - Fatal error!'
    write ( *, '(a)' ) '  The polynomial degree must be nonnegative.'
    ierror = 1
    return
  end if

  nint = (ndat-1) / ndeg
!
!  Make sure we have enough room in WORK.
!
  if ( ndat > maxrwork ) then
    write ( *, '(a)' ) 'PWP - Fatal error!'
    write ( *, '(a)' ) '  There is not enough work space to set up'
    write ( *, '(a)' ) '  the piecewise polynomial.'
    ierror = 1
    return
  end if
!
!  Fill in the graph values, XVAL, YVAL, NVAL.
!
!  The code for NDEG = 0 ignores the fact that we need NINT+1 X values
!  and NINT Y values.
!
  if ( ndeg == 0 ) then

    do i = 1, nint

      nval = nval+1
      xval(nval) = xdat(i)
      yval(nval) = ydat(i)

      nval = nval+1
      xval(nval) = xdat(i+1)
      yval(nval) = ydat(i)

    end do

  else if ( ndeg == 1 ) then

    do i = 1, ndat
      nval = nval+1
      xval(nval) = xdat(i)
      yval(nval) = ydat(i)
    end do

  else

    nval = nval+1
    xval(nval) = xdat(1)
    yval(nval) = ydat(1)

    do i = 1, nint

      ilo = 1 + (i-1) * ndeg
      ihi = 1 + i * ndeg
      call data_to_dif ( rwork, ndeg+1, xdat(ilo), ydat(ilo) )
      jhi = (nfine/nint) + 1

      do j = 1, jhi
        nval = nval+1
        call rvec_even_select ( xdat(ilo), xdat(ihi), jhi+1, j, xval(nval) )
        call dif_val ( rwork, ndeg+1, xdat(ilo), xval(nval), yval(nval) )
        ival(nval) = 1
      end do

    end do

  end if

  ival(nval) = 0

  return
end
subroutine r2_cheby ( n, alo, ahi, a )
!
!*******************************************************************************
!
!! R2_CHEBY sets up the Chebyshev abscissas in a real interval.
!
!  Discussion:
!
!    The routine sets up a vector of X values spaced between the values
!    XLO and XHI in a similar way to the spacing of the Chebyshev
!    points of the same order in the interval [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points to compute.
!
!    Input, real ALO, AHI, the range.
!
!    Output, real A(N), the computed X values.
!
  implicit none
!
  integer n
!
  real a(n)
  real ahi
  real alo
  real arg
  integer i
  real, parameter :: pi = 3.14159265358979E+00

  if ( n == 1 ) then

    a(1) = 0.5E+00 * ( alo + ahi )

  else if ( n > 1 ) then

    do i = 1, n

      arg = real ( 2 * i - 1 ) * pi / real ( 2 * n )

      a(i) = 0.5E+00 * ( ( 1.0E+00 + cos ( arg ) ) * alo &
        + ( 1.0E+00 - cos ( arg ) ) * ahi )

    end do

  end if

  return
end
subroutine r_next ( s, r, done )
!
!*******************************************************************************
!
!! R_NEXT "reads" real numbers from a string, one at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string, presumably containing real
!    numbers.  These may be separated by spaces or commas.
!
!    Output, real R.  If DONE is FALSE, then R contains the
!    "next" real value read from the string.  If DONE is TRUE, then
!    R is zero.
!
!    Input/output, logical DONE.
!    On input with a fresh string, the user should set DONE to TRUE.
!    On output, the routine sets DONE to FALSE if another real
!    value was read, or TRUE if no more reals could be read.
!
  implicit none
!
  logical done
  integer ierror
  integer lchar
  integer, save :: next = 1
  real r
  character ( len = * ) s
!
  r = 0.0E+00

  if ( done ) then
    next = 1
    done = .false.
  end if

  if ( next > len ( s ) ) then
    done = .true.
    return
  end if

  call s_to_r ( s(next:), r, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    done = .true.
    next = 1
  else
    done = .false.
    next = next + lchar
  end if

  return
end
subroutine r_to_s_left ( r, s )
!
!*******************************************************************************
!
!! R_TO_S_LEFT writes a real into a left justified character string.
!
!  Discussion:
!
!    A 'G14.6' format is used with a WRITE statement.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the real number to be written into the string.
!
!    Output, character ( len = * ) S, the string into which
!    the real number is to be written.  If the string is less than 14
!    characters long, it will will be returned as a series of
!    asterisks.
!
  implicit none
!
  integer i
  integer nchar
  real r
  character ( len = * ) s
  character ( len = 14 ) s2
!
  nchar = len ( s )

  if ( nchar < 14 ) then

    do i = 1, nchar
      s(i:i) = '*'
    end do

  else if ( r == 0.0E+00 ) then
    s(1:14) = '     0.0E+00      '
  else
    write ( s2, '(g14.6)' ) r
    s(1:14) = s2
  end if
!
!  Shift the string left.
!
  s = adjustl ( s )

  return
end
subroutine rhpsrt ( k, n, lda, a, map )
!
!***********************************************************************
!
!! RHPSRT sorts points into lexicographic order using heap sort.
!
!  Discussion:
!
!    This routine uses heapsort to obtain the permutation of N K-dimensional
!    points so that the points are in lexicographic increasing order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Phone: (403) 492-5757
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer K, the dimension of the points (for instance, 2
!    for points in the plane).
!
!    Input, integer N, the number of points.
!
!    Input, integer LDA, the leading dimension of array A in the calling
!    routine; LDA should be at least K.
!
!    Input, real A(LDA,N); A(I,J) contains the I-th coordinate
!    of point J.
!
!    Input/output, integer MAP(N).
!    On input, the points of A with indices MAP(1), MAP(2), ...,
!    MAP(N) are to be sorted
!
!    On output, MAP contains a permutation of its input values, with the
!    property that, lexicographically,
!      A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N))
!
  implicit none
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer k
  integer map(n)
!
  do i = n/2, 1, -1
    call rsftdw ( i, n, k, lda, a, map )
  end do

  do i = n, 2, -1
    call i_swap ( map(1), map(i) )
    call rsftdw ( 1, i-1, k, lda, a, map )
  end do

  return
end
function rless ( k, p, q )
!
!***********************************************************************
!
!! RLESS determines whether P is lexicographically less than Q.
!
!  Discussion:
!
!    P and Q are K-dimensional points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Phone: (403) 492-5757
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer K, the spatial dimension of the points.
!
!    Input, real P(K), Q(K), the points to be compared.
!
!    Output, logical RLESS, is TRUE if P < Q, FALSE otherwise.
!
  implicit none
!
  integer k
!
  real cmax
  integer i
  real p(k)
  real q(k)
  logical rless
  real tol
!
  tol = 100.0E+00 * epsilon ( tol )

  do i = 1, k

    cmax = max ( abs ( p(i) ), abs ( q(i) ) )

    if ( abs ( p(i) - q(i) ) > tol * cmax .and. cmax > tol ) then

      if ( p(i) < q(i) ) then
        rless = .true.
      else
        rless = .false.
      end if

      return

    end if

  end do

  rless = .false.

  return
end
subroutine rpnchk ( ierror, ihi, ilo, iopsym, irpn, maxrpn, maxsym )
!
!*******************************************************************************
!
!! RPNCHK examines an RPN formula, looking for a complete RPN expression.
!
!  Discussion:
!
!    The routine starts at location IHI, and finds the position ILO
!    such that IRPN(ILO)...IRPN(IHI) represents a single argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, 0, no error; 1 an error.
!
!    Input, integer IHI, the location in IRPN where the search begins.
!
!    Output, integer ILO, the location in IRPN such that IRPN(ILO) through
!    IRPN(IHI) represents a single argument, or IHI+1 if no such location
!    was found.
!
!    Input, integer IOPSYM(MAXSYM), contains, for each symbol, the number
!    of operands.  In particular, if a symbol represents a constant,
!    IOPSYM(I) is 0.  If a symbol represents a unary operator such as ABS,
!    IOPSYM(I) is 1.  IOPSYM may be as large as 3.
!
!    Workspace, integer IRPN(MAXRPN), used to store the compiled
!    versions of user formulas.
!
!    Input, integer MAXRPN, specifies the length of IRPN.
!
!    Input, integer MAXSYM, the maximum number of symbols allowed.
!
  implicit none
!
  integer maxrpn
  integer maxsym
!
  integer ierror
  integer ihi
  integer ilo
  integer iopsym(maxsym)
  integer irpn(maxrpn)
  integer isum
!
  isum = 0
  ierror = 0
  ilo = ihi + 1

  do

    ilo = ilo - 1

    if ( ilo <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RPNCHK - Error!'
      write ( *, '(a)' ) '  Reached beginning of formula without matchup.'
      ierror = 1
      return
    end if

    if ( iopsym(irpn(ilo)) < 0 .and. isum == 0 ) then
      cycle
    end if

    isum = isum + 1 - iopsym(irpn(ilo))

    if ( isum == 1 ) then
      exit
    end if

  end do

  return
end
subroutine rpnset ( ierror, ifinis, imin, intsym, iopsym, iprsym, irpn, &
  istack, maxrpn, maxsym, nints, nrpn, symbol )
!
!*******************************************************************************
!
!! RPNSET converts the infix formula into an RPN formula.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer IFINIS, the index in SYMBOL of the symbol
!    '$', meaning the end of the formula.
!
!    Input, integer IMIN, an offset for accessing IRPN.
!
!    Input, integer INTSYM(80), a set of integers which
!    are the indices of symbols, representing an infix formula.
!
!    Input, integer IOPSYM(MAXSYM), contains, for each
!    symbol, the number of operands.  In particular, if
!    a symbol represents a constant, IOPSYM(I) is 0.
!    If a symbol represents a unary operator such as ABS,
!    IOPSYM(I) is 1.  IOPSYM may be as large as 3.
!
!    Input, integer IPRSYM(MAXSYM), the priority of each symbol.
!
!    Output, integer IRPN(MAXRPN), the RPN version of the infix formula.
!
!    Workspace, integer ISTACK(MAXSYM).
!
!    Input, integer MAXRPN, the maximum number of RPN symbols allowed.
!
!    Input, integer MAXSYM, the maximum number of symbols allowed.
!
!    Input, integer NINTS, the number of useful entries in INTSYM.
!
!    Output, integer NRPN, the number of useful entries in IRPN.
!
!    Input, character ( len = MAXLEN ) SYMBOL(MAXSYM), the symbolic names
!    of the internally defined constants and functions.
!
  implicit none
!
  integer, parameter :: maxlen = 20
  integer maxrpn
  integer maxsym
!
  integer ierror
  integer ifinis
  integer imin
  integer intsym(80)
  integer iopsym(maxsym)
  integer iprsym(maxsym)
  integer iread
  integer irpn(maxrpn)
  integer istack(maxsym)
  integer istak
  integer isym
  integer jsym
  integer nints
  integer nrpn
  character ( len = maxlen ) sym
  character ( len = maxlen ) sym2
  character ( len = maxlen ) symbol(maxsym)
!
  ierror = 0

  nrpn = 0
!
!  An error if the formula seems to have nothing in it.
!
  if ( nints <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RPNSET - Error!'
    write ( *, '(a)' ) '  This formula seems to be blank!'
    ierror = 1
    return
  end if

  iread = 0
  istak = 0
!
!  Read symbol number IREAD from INTSYM.
!
  do

    iread = iread + 1

    if ( iread > nints ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RPNSET - Error!'
      write ( *, '(a)' ) '  This formula does not make sense!'
      return
    end if

    isym = intsym(iread)
    if ( isym == 0 ) then
      cycle
    end if

    sym = symbol(isym)
    if ( sym == ',' ) then
      cycle
    end if
!
!  If the symbol is "$", then it's time to pop the stack.
!
    if ( sym == '$' ) then

      do

        if ( istak <= 0 ) then

          nrpn = nrpn + 1
          irpn(nrpn+imin) = ifinis

          if ( iread < nints ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'RPNSET - Error!'
            write ( *, '(a)' ) '  Some of the formula is left over!'
            ierror = 1
          end if

          return
        end if

        nrpn = nrpn + 1
        irpn(nrpn+imin) = istack(istak)
        istak = istak - 1

      end do

    end if
!
!  A left parenthesis goes onto the stack.
!
    if ( sym == '(' ) then
      istak = istak + 1
      istack(istak) = isym
      cycle
    end if
!
!  A variable or constant goes immediately into IRPN.
!
    if ( iopsym(isym) == 0 ) then
      nrpn = nrpn + 1
      irpn(nrpn+imin) = isym
      cycle
    end if
!
!  Put operators onto the stack.
!
    do

      if ( istak <= 0 ) then
        istak = istak + 1
        istack(istak) = isym
        exit
      end if

      jsym = istack(istak)

      if ( iprsym(jsym) <= iprsym(isym) ) then
        jsym = istack(istak)
        sym2 = symbol(jsym)

        if ( sym == ')' .and. sym2 == '(' ) then
          istak = istak - 1
        else
          istak = istak + 1
          istack(istak) = isym
        end if
        exit
      end if

      nrpn = nrpn + 1
      irpn(nrpn+imin) = istack(istak)
      istak = istak - 1

    end do

  end do

  return
end
subroutine rpnval ( ierror, ifree, imin, iopsym, iprsym, ipval, &
  irad, irpn, irtol, istack, maxrpn, maxsym, maxval, &
  nrpn, nsym, nsyms, symbol, valsym, value )
!
!*******************************************************************************
!
!! RPNVAL evaluates the symbolic functions in an RPN formula.
!
!  Discussion:
!
!    RPNVAL determines the number of arguments, gathering their values,
!    and calling FUNVAL to evaluate the given function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real EPS, the value of the machine precision.
!
!    Output, integer IERROR.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer IFREE, the index of the next free location in VALSYM.
!
!    Input, integer IMIN, an offset for indexing IRPN.
!
!    Input, integer IOPSYM(MAXSYM), contains, for each
!    symbol, the number of operands.  In particular, if
!    a symbol represents a constant, IOPSYM(I) is 0.
!    If a symbol represents a unary operator such as ABS,
!    IOPSYM(I) is 1.  IOPSYM may be as large as 3.
!
!    Input, integer IPRSYM(MAXSYM), the relative priorities
!    of the functions.
!
!    Input, integer IPVAL(MAXSYM), contains, for each symbol,
!    the address in VALSYM where the value of the symbol
!    is stored, if it is a scalar.  If the symbol represents
!    a vector or matrix, IPVAL(IARG) points to location of
!    the first entry.
!
!    Workspace, integer IRPN(MAXRPN), used to store the compiled
!    versions of user formulas.
!
!    Workspace, integer ISTACK(MAXSYM), workspace for interpreting
!    the formula.
!
!    Input, integer MAXRPN, specifies the length of IRPN.
!
!    Input, integer MAXSYM, the maximum number of symbols allowed.
!
!    Input, integer MAXVAL, the maximum number of values allowed in VALSYM.
!
!    Output, integer NCOL, the number of columns in the result.
!
!    Output, integer NROW, the number of rows in the result.
!
!    Input, integer NRPN, the number of useful entries in IRPN.
!
!    Input, integer NSYM, the total number of symbols, including temporaries.
!
!    Input, integer NSYMS, the number of declared symbols.
!
!    Input, integer NUMDIM(2,MAXSYM).
!    For each symbol I, NUMDIM(1,I) is the number of rows in its
!    value, and NUMDIM(2,I) is the number of columns.
!
!    Input, character ( len = * ) SYMBOL(MAXSYM), the symbolic names
!    of the internally defined constants and functions.
!
!    Input, real VALSYM(MAXVAL), contains the values of all the
!    symbolic variables.
!
!    Workspace, real VALUE, space used to hold
!    the value of the variable to be saved.
!
  implicit none
!
  integer, parameter :: maxlen = 20
  integer maxrpn
  integer maxsym
!
  character ( len = 3 ) ctemp
  integer iarg1
  integer iarg2
  integer ierror
  integer ifree
  integer imin
  integer index1
  integer index4
  integer iopsym(maxsym)
  integer iprsym(maxsym)
  integer ipset
  integer ipval(maxsym)
  integer irad
  integer iread
  integer irpn(maxrpn)
  integer irtol
  integer istack(maxsym)
  integer istak
  integer isym
  integer isymo
  integer itemp
  integer maxval
  integer nrpn
  integer nsym
  integer nsyms
  character ( len = maxlen ) sym
  character ( len = maxlen ) symbol(maxsym)
  real valsym(maxval)
  real value
!
  ierror = 0
  value = 0.0E+00

  if ( nrpn <= 1 ) then
    return
  end if

  istak = 0
  isymo = 0
  nsyms = nsym
  itemp = 0
  iread = 0

  do

    iread = iread + 1

    if ( iread > nrpn ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RPNVAL - Error!'
      write ( *, '(a)' ) '  The formula ended unexpectedly!'
      exit
    end if

    isym = irpn(iread+imin)

    if ( 0 < isym .and. isym <= maxsym ) then
      sym = symbol(isym)
    else
      sym = '$'
    end if

    if ( sym == '$' .or. iread > nrpn ) then
      index4 = ipval(isymo)
      value = valsym(index4)
      exit
    end if
!
!  Constants and variables go into stack.
!
    if ( iopsym(isym) == 0 ) then

      if ( isym <= nsym ) then
        istak = istak+1
        istack(istak) = isym
        isymo = isym
      else

        if ( nsyms >= maxsym ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'RPNVAL - Error!'
          write ( *, '(a)' ) '  There is no memory left for more symbols!'
          write ( *, '(a)' ) '  Perhaps the KILL or INIT command would help.'
          ierror = 1
          exit
        end if

        nsyms = nsyms+1
        isymo = nsyms
        istak = istak+1
        istack(istak) = nsyms
        iprsym(nsyms) = 10
        iopsym(nsyms) = 0
        itemp = itemp+1
        ctemp = ' '
        call i_to_s_zero ( itemp, ctemp )

        symbol(nsyms) = 'STK000'
        symbol(nsyms)(4:6) = ctemp
        ipval(nsyms) = ifree
        ifree = ifree + 1
        index1 = ipval(isym)
        index4 = ipval(nsyms)
        valsym(index4) = valsym(index1)
      end if

      cycle

    end if
!
!  Pull off arguments.
!
    iarg1 = istack(istak)
    iarg2 = 0

    if ( iopsym(isym) == 2 ) then
      iarg2 = istack(istak)
      istak = istak-1
      iarg1 = istack(istak)
    else if (iopsym(isym) == 3 ) then
      istak = istak-1
      iarg2 = istack(istak)
      istak = istak-1
      iarg1 = istack(istak)
    end if

    sym = symbol(isym)

    call funval ( iarg1, iarg2, ierror, ifree, iopsym, iprsym, &
      ipset, ipval, irad, irtol, itemp, maxsym, maxval, nsyms, &
      sym, symbol, valsym )

    isymo = nsyms

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RPNVAL - Warning!'
      write ( *, '(a)' ) '  Evaluation of this formula is abandoned.'
      value = 0.0E+00
      exit
    end if

    if ( ipset == 0 ) then
      ipval(nsyms) = ifree
      ifree = ifree + 1
    else
      ipval(nsyms) = ipset
      ipset = 0
    end if

    istack(istak) = nsyms

  end do

  return
end
subroutine rvec_bracket ( n, x, xval, left, right )
!
!*******************************************************************************
!
!! RVEC_BRACKET searches a sorted array for successive brackets of a value.
!
!  Discussion:
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of input array.
!
!    Input, real X(N), an array sorted into ascending order.
!
!    Input, real XVAL, a value to be bracketed.
!
!    Output, integer LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      XVAL > X(N), when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none
!
  integer n
!
  integer i
  integer left
  integer right
  real x(n)
  real xval
!
  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end
function rvec_distinct ( n, a )
!
!*******************************************************************************
!
!! RVEC_DISTINCT is true if the entries in a real vector are distinct.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real A(N), the vector to be checked.
!
!    Output, logical RVEC_DISTINCT is .TRUE. if all N elements of A
!    are distinct.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  integer j
  logical rvec_distinct
!
  rvec_distinct = .false.

  do i = 2, n
    do j = 1, i - 1
      if ( a(i) == a(j) ) then
        return
      end if
    end do
  end do

  rvec_distinct = .true.

  return
end
subroutine rvec_even ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! RVEC_EVEN returns N real values, evenly spaced between ALO and AHI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the low and high values.
!
!    Input, integer N, the number of values.
!
!    Output, real A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none
!
  integer n
!
  real a(n)
  real ahi
  real alo
  integer i
!
  if ( n == 1 ) then

    a(1) = 0.5E+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i ) * alo + real ( i - 1 ) * ahi ) / real ( n - 1 )
    end do

  end if

  return
end
subroutine rvec_even_select ( xlo, xhi, n, ival, xval )
!
!*******************************************************************************
!
!! RVEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
!
!  Discussion:
!
!    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / real ( N - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XLO, XHI, the low and high values.
!
!    Input, integer N, the number of values.
!
!    Input, integer IVAL, the index of the desired point.
!    IVAL is normally between 1 and N, but may be any
!    integer value.
!
!    Output, real XVAL, the IVAL-th of N evenly spaced values
!    between XLO and XHI.
!
!    Unless N = 1, X(1) = XLO and X(N) = XHI.
!
!    If N = 1, then X(1) = 0.5*(XLO+XHI).
!
  implicit none
!
  integer n
!
  integer ival
  real xhi
  real xlo
  real xval
!
  if ( n == 1 ) then

    xval = 0.5E+00 * ( xlo + xhi )

  else

    xval = ( real ( n - ival ) * xlo + real ( ival - 1 ) * xhi ) &
      / real ( n - 1 )

  end if

  return
end
subroutine rvec_unit_euclidean ( n, a )
!
!*******************************************************************************
!
!! RVEC_UNIT_EUCLIDEAN Euclidean normalizes a vector in ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vector.
!
!    Input/output, A(N), the vector to be normalized.
!
  implicit none
!
  integer n
!
  real a(n)
  real norm
!
  norm = sqrt ( sum ( a(1:n)**2 ) )

  if ( norm /= 0.0E+00 ) then
    a(1:n) = a(1:n) / norm
  end if

  return
end
subroutine s3_fs ( a1, a2, a3, n, b, x )
!
!*******************************************************************************
!
!! S3_FS factors and solves a tridiagonal linear system.
!
!  Discussion:
!
!    This algorithm requires that each diagonal entry be nonzero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real A1(2:N), A2(1:N), A3(1:N-1).
!    On input, the nonzero diagonals of the linear system.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, integer N, the order of the linear system.
!
!    Input/output, real B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B has been overwritten by factorization information.
!
!    Output, real X(N), the solution of the linear system.
!
  implicit none
!
  integer n
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real b(n)
  integer i
  real x(n)
  real xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a2(i) == 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'S3_FS - Fatal error!'
      write ( *, '(a,i6,a)' ) '  A2(', i, ') = 0.'
      return
    end if
  end do

  do i = 2, n-1

    xmult = a1(i) / a2(i-1)
    a2(i) = a2(i) - xmult * a3(i-1)

    b(i) = b(i) - xmult * b(i-1)

  end do

  xmult = a1(n) / a2(n-1)
  a2(n) = a2(n) - xmult * a3(n-1)

  x(n) = ( b(n) - xmult * b(n-1) ) / a2(n)
  do i = n-1, 1, -1
    x(i) = ( b(i) - a3(i) * x(i+1) ) / a2(i)
  end do

  return
end
subroutine s_blank_delete ( s )
!
!*******************************************************************************
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discssion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none
!
  character c
  integer iget
  integer iput
  integer nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
!
  iput = 0
  nchar = len_trim ( s )

  do iget = 1, nchar

    c = s(iget:iget)

    if ( c /= ' ' .and. c /= TAB ) then
      iput = iput + 1
      s(iput:iput) = c
    end if

  end do

  s(iput+1:nchar) = ' '

  return
end
subroutine s_blanks_delete ( s )
!
!*******************************************************************************
!
!! S_BLANKS_DELETE replaces consecutive blanks by one blank.
!
!  Discussion:
!
!    The remaining characters are left justified and right padded with blanks.
!    TAB characters are converted to spaces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none
!
  integer i
  integer j
  character newchr
  character oldchr
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
!
  j = 0
  newchr = ' '

  do i = 1, len ( s )

    oldchr = newchr
    newchr = s(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    s(i:i) = ' '

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      s(j:j) = newchr
    end if

  end do

  return
end
subroutine s_cap ( s )
!
!*******************************************************************************
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none
!
  character c
  integer i
  integer nchar
  character ( len = * ) s
!
  nchar = len_trim ( s )

  do i = 1, nchar

    c = s(i:i)
    call ch_cap ( c )
    s(i:i) = c

  end do

  return
end
function s_eqi ( s1, s2 )
!
!*******************************************************************************
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none
!
  character c1
  character c2
  integer i
  integer len1
  integer len2
  integer lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2
!
  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
function s_is_alpha ( s )
!
!*******************************************************************************
!
!! S_IS_ALPHA returns .TRUE. if the string contains only alphabetic characters.
!
!  Discussion:
!
!    Here, alphabetic characters are 'A' through 'Z' and 'a' through 'z'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be checked.
!
!    Output, logical S_IS_ALPHA, .TRUE. if the string contains only
!    alphabetic characters, .FALSE. otherwise.
!
  implicit none
!
  logical ch_is_alpha
  integer i
  character ( len = * ) s
  logical s_is_alpha
!
  s_is_alpha = .false.

  do i = 1, len ( s )

    if ( .not. ch_is_alpha ( s(i:i) ) ) then
      return
    end if

  end do

  s_is_alpha = .true.

  return
end
function s_paren_check ( s )
!
!*******************************************************************************
!
!! S_PAREN_CHECK checks the parentheses in a string.
!
!  Discussion:
!
!    Blanks are removed from the string, and then the following checks
!    are made:
!
!    1) as we read the string left to right, there must never be more
!       right parentheses than left ones;
!    2) there must be an equal number of left and right parentheses;
!    3) there must be no occurrences of the abutting packages '...)(...'.
!    4) there must be no occurrences of the empty package '()'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to check.
!
!    Output, logical S_PAREN_CHECK is TRUE if the string passed the checks.
!
  implicit none
!
  integer i
  integer isum
  character ( len = * ) s
  character ( len = 256 ) s_copy
  integer s_len
  logical s_paren_check
!
  s_copy = s
  call s_blank_delete ( s_copy)

  s_len = len_trim ( s_copy )
!
!  1) Letting '(' = +1 and ')' = -1, check that the running parentheses sum
!  is always nonnegative.
!
  isum = 0
  do i = 1, s_len

    if ( s_copy(i:i) == '(' ) then
      isum = isum + 1
    end if

    if ( s_copy(i:i) == ')' ) then

      isum = isum - 1

      if ( isum < 0 ) then
        s_paren_check = .false.
        return
      end if

    end if

  end do
!
!  2) Check that the final parentheses sum is zero.
!
  if ( isum /= 0 ) then
    s_paren_check = .false.
    return
  end if
!
!  3) Check that there are no ")(" pairs.
!
  do i = 2, s_len
    if ( s_copy(i-1:i) == ')(' ) then
      s_paren_check = .false.
      return
    end if
  end do
!
!  4) Check that there are no "()" pairs.
!
  do i = 2, s_len

    if ( s_copy(i-1:i) == '()' ) then
      s_paren_check = .false.
      return
    end if

  end do
!
!  The checks were passed.
!
  s_paren_check = .true.

  return
end
subroutine s_to_r ( s, r, ierror, lchar )
!
!*******************************************************************************
!
!! S_TO_R reads a real number from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0E+00
!    '     1   '       1.0E+00
!    '1A'              1.0E+00
!    '12,34,56'        12.0E+00
!    '  34 7'          34.0E+00
!    '-1E2ABCD'        -100.0E+00
!    '-1X2ABCD'        -1.0E+00
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0E+00
!    '17d2'            1700.0E+00
!    '-14e-2'         -0.14
!    'e2'              100.0E+00
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real R, the real value that was read from the string.
!
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none
!
  logical ch_eqi
  character c
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer lchar
  integer nchar
  integer ndig
  real r
  real rbot
  real rexp
  real rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
!
  nchar = len_trim ( s )
  ierror = 0
  r = 0.0E+00
  lchar = - 1
  isgn = 1
  rtop = 0.0E+00
  rbot = 1.0E+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( ihave > 1 ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( ihave >= 6 .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0E+00 * rtop + real ( ndig )
      else if ( ihave == 5 ) then
        rtop = 10.0E+00 * rtop + real ( ndig )
        rbot = 10.0E+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 .or. lchar+1 >= nchar ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0E+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0E+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0E+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )
!
!*******************************************************************************
!
!! SPLINE_CUBIC_SET computes the second derivatives of a cubic spline.
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to
!    determine the second derivative data, passing in the data to be
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) )
!             + C(IVAL) * ( T - T(IVAL) )**2
!             + D(IVAL) * ( T - T(IVAL) )**3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      =
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL)
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!
!    Boundary conditions must be applied at the first and last knots.
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points; N must be at least 2.
!
!    Input, real T(N), the points where data is specified.
!    The values should be distinct, and increasing.
!
!    Input, real Y(N), the data values to be interpolated.
!
!    Input, integer IBCBEG, the left boundary condition flag:
!
!      0: the spline should be a quadratic over the first interval;
!      1: the first derivative at the left endpoint should be YBCBEG;
!      2: the second derivative at the left endpoint should be YBCBEG.
!
!    Input, real YBCBEG, the left boundary value, if needed.
!
!    Input, integer IBCEND, the right boundary condition flag:
!
!      0: the spline should be a quadratic over the last interval;
!      1: the first derivative at the right endpoint should be YBCEND;
!      2: the second derivative at the right endpoint should be YBCEND.
!
!    Input, real YBCEND, the right boundary value, if needed.
!
!    Output, real YPP(N), the second derivatives of the cubic spline.
!
  implicit none
!
  integer n
!
  real diag(n)
  integer i
  integer ibcbeg
  integer ibcend
  real sub(2:n)
  real sup(1:n-1)
  real t(n)
  real y(n)
  real ybcbeg
  real ybcend
  real ypp(n)
!
!  Check.
!
  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The number of knots must be at least 2.'
    write ( *, '(a,i6)' ) '  The input value of N = ', n
    stop
  end if

  do i = 1, n-1
    if ( t(i) >= t(i+1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
      write ( *, '(a)') '  The knots must be strictly increasing, but'
      write ( *, '(a,i6,a,g14.6)' ) '  T(',  i,') = ', t(i)
      write ( *, '(a,i6,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
      stop
    end if
  end do
!
!  Set the first equation.
!
  if ( ibcbeg == 0 ) then
    ypp(1) = 0.0E+00
    diag(1) = 1.0E+00
    sup(1) = -1.0E+00
  else if ( ibcbeg == 1 ) then
    ypp(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
    diag(1) = ( t(2) - t(1) ) / 3.0E+00
    sup(1) = ( t(2) - t(1) ) / 6.0E+00
  else if ( ibcbeg == 2 ) then
    ypp(1) = ybcbeg
    diag(1) = 1.0E+00
    sup(1) = 0.0E+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCBEG must be 0, 1 or 2.'
    write ( *, '(a,i6)' ) '  The input value is IBCBEG = ', ibcbeg
    stop
  end if
!
!  Set the intermediate equations.
!
  do i = 2, n-1
    ypp(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
           - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
    sub(i) = ( t(i) - t(i-1) ) / 6.0E+00
    diag(i) = ( t(i+1) - t(i-1) ) / 3.0E+00
    sup(i) = ( t(i+1) - t(i) ) / 6.0E+00
  end do
!
!  Set the last equation.
!
  if ( ibcend == 0 ) then
    ypp(n) = 0.0E+00
    sub(n) = -1.0E+00
    diag(n) = 1.0E+00
  else if ( ibcend == 1 ) then
    ypp(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
    sub(n) = ( t(n) - t(n-1) ) / 6.0E+00
    diag(n) = ( t(n) - t(n-1) ) / 3.0E+00
  else if ( ibcend == 2 ) then
    ypp(n) = ybcend
    sub(n) = 0.0E+00
    diag(n) = 1.0E+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCEND must be 0, 1 or 2.'
    write ( *, '(a,i6)' ) '  The input value is IBCEND = ', ibcend
    stop
  end if
!
!  Special case:
!    N = 2, IBCBEG = IBCEND = 0.
!
  if ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then

    ypp(1) = 0.0E+00
    ypp(2) = 0.0E+00
!
!  Solve the linear system.
!
  else

    call s3_fs ( sub, diag, sup, n, ypp, ypp )

  end if

  return
end
subroutine spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )
!
!*******************************************************************************
!
!! SPLINE_CUBIC_VAL evaluates a cubic spline at a specific point.
!
!  Discussion:
!
!    SPLINE_CUBIC_SET must have already been called to define the
!    values of YPP.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A
!             + B * ( T - T(IVAL) )
!             + C * ( T - T(IVAL) )**2
!             + D * ( T - T(IVAL) )**3
!
!    Here:
!      A = Y(IVAL)
!      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C = YPP(IVAL) / 2
!      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real T(N), the knot values.
!
!    Input, real Y(N), the data values at the knots.
!
!    Input, real YPP(N), the second derivatives of the spline at
!    the knots.
!
!    Input, real TVAL, a point, typically between T(1) and T(N), at
!    which the spline is to be evalulated.  If TVAL lies outside
!    this range, extrapolation is used.
!
!    Output, real YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
  implicit none
!
  integer n
!
  real dt
  real h
  integer left
  integer right
  real t(n)
  real tval
  real y(n)
  real ypp(n)
  real yppval
  real ypval
  real yval
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!
  call rvec_bracket ( n, t, tval, left, right )
!
!  Evaluate the polynomial.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  yval = y(left) &
       + dt * ( ( y(right) - y(left) ) / h &
              - ( ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00 ) * h &
       + dt * ( 0.5E+00 * ypp(left) &
       + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0E+00 * h ) ) ) )

  ypval = ( y(right) - y(left) ) / h &
       - ( ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00 ) * h &
       + dt * ( ypp(left) &
       + dt * ( 0.5E+00 * ( ypp(right) - ypp(left) ) / h ) )

  yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h

  return
end
subroutine star ( ierror, ival, val_max, nval, xval, yval )
!
!***********************************************************************
!
!! STAR draws a star.
!
!  Discussion:
!
!    Point I on the circle is connected to point I+ISTAR, where
!    ISTAR is specified by the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
!    Input/output, integer IVAL(VAL_MAX); IVAL(I) is 1 if point I is to
!    be connected to point I+1.
!
!    Input, integer VAL_MAX, the maximum number of plot points.
!
!    Input/output, integer NVAL, the current number of plot points.
!
!    Input/output, real XVAL(VAL_MAX), YVAL(VAL_MAX), the plot points.
!
  implicit none
!
  integer val_max
!
  real angle
  integer i
  integer ierror
  integer ios
  integer ip
  integer iq
  integer ival(val_max)
  integer j
  integer nval
  real, parameter :: pi = 3.14159265358979E+00
  real xval(val_max)
  real yval(val_max)
!
  ierror = 0

  do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STAR'
    write ( *, '(a)' ) '  The star is constructed by connecting points'
    write ( *, '(a)' ) '  that lie on a circle.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  To specify the star you want, please enter'
    write ( *, '(a)' ) '  the number of points on the circle.'

    read ( *, *, iostat = ios ) iq

    if ( ios /= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STAR - Fatal error!'
      write ( *, '(a)' ) '  I/O error.'
      return
    else if ( ierror > 5 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STAR - Fatal error!'
      write ( *, '(a)' ) '  Too many errors in a row!'
      return
    else if ( iq  <=  2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Please choose a value that is more than 2!'
      ierror = ierror + 1
    else if ( nval + 2 * iq > val_max ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STAR'
      write ( *, '(a)' ) '  There is not enough memory left.'
      write ( *, '(a,i6)' ) '  Choose a value no more than ', &
        ( val_max - nval ) / 2
      ierror = ierror + 1
    else
      exit
    end if

  end do

  ierror = 0

  do 

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the point to which point 1 is connected.'
    write ( *, '(a,i6)' ) 'This must be a value between 2 and ', iq

    read ( *, *, iostat = ios ) ip

    if ( ios /= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STAR - Fatal error!'
      write ( *, '(a)' ) '  I/O error.'
      return
    else if ( ierror > 5 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STAR - Fatal error!'
      write ( *, '(a)' ) '  Too many errors in a row!'
      return
    else if ( ip < 2 ) then
      ierror = ierror + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Please choose a value that is at least 2!'
    else if ( ip > iq ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'Please choose a value that is no more than ', iq
      ierror = ierror + 1
    else
      exit
    end if

  end do
!
!  Connect points
!    1 and P
!    2 and P+1
!      ...
!    Q and P+Q
!
  do i = 1, iq

    nval = nval+1
    angle = ( 2.0E+00 * pi * real ( i ) ) / real ( iq )
    xval(nval) = cos ( angle )
    yval(nval) = sin ( angle )
    ival(nval) = 1

    j = mod ( i + ( ip - 1 ), iq )

    nval = nval+1
    angle = ( 2.0E+00 * pi * real ( j ) ) / real ( iq )
    xval(nval) = cos ( angle )
    yval(nval) = sin ( angle )
    ival(nval) = 0

  end do

  return
end
subroutine swapec ( i, top, maxst, btri, bedg, vcl, til, tnbr, stack, ierr )
!
!***********************************************************************
!
!! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!
!  Discussion: 
!
!    The routine swaps diagonal edges in a 2-D triangulation, based on
!    the empty circumcircle criterion, until all triangles are Delaunay, 
!    given that I is the index of the new vertex added to triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Phone: (403) 492-5757
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer I, the index in VCL of the new vertex.
!
!    Input/output, integer TOP, the index of the top of the stack.
!    On output, TOP is zero.
!
!    Input, integer MAXST, the maximum size available for the STACK array.
!
!    Input/output, integer BTRI, BEDG; on input, if positive, are the 
!    triangle and edge indices of a boundary edge whose updated indices 
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, real VCL(2,*), the coordinates of the vertices.
!
!    Input/output, integer TIL(3,*), the triangle incidence list.  May be
!    updated on output because of swaps.
!
!    Input/output, integer TNBR(3,*), the triangle neighbor list; negative
!    values are used for links of the counter-clockwise linked list of 
!    boundary edges; May be updated on output because of swaps.
!
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer STACK(1:MAXST); on input, entries 1 through TOP
!    contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer IERR is set to 8 for abnormal return.
!
  implicit none
!
  integer maxst
!
  integer a
  integer b
  integer bedg
  integer btri
  integer c
  integer diaedg
  integer e
  integer ee
  integer em1
  integer ep1
  integer f
  integer fm1
  integer fp1
  integer i
  integer ierr
  integer l
  integer r
  integer s
  integer stack(maxst)
  integer swap
  integer t
  integer til(3,*)
  integer tnbr(3,*)
  integer top
  integer tt
  integer u
  real vcl(2,*)
  real x
  real y
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  ierr = 0
  x = vcl(1,i)
  y = vcl(2,i)

  do

    if ( top <= 0 ) then
      exit
    end if

    t = stack(top)
    top = top - 1

    if ( til(1,t) == i ) then
      e = 2
      b = til(3,t)
    else if ( til(2,t) == i ) then
      e = 3
      b = til(1,t)
    else
      e = 1
      b = til(2,t)
    end if

    a = til(e,t)
    u = tnbr(e,t)

    if ( tnbr(1,u) == t ) then
      f = 1
      c = til(3,u)
    else if ( tnbr(2,u) == t ) then
      f = 2
      c = til(1,u)
    else
      f = 3
      c = til(2,u)
    end if

    swap = diaedg ( x, y, vcl(1,a), vcl(2,a), vcl(1,c), vcl(2,c), &
      vcl(1,b), vcl(2,b) )

    if ( swap == 1 ) then

      em1 = e - 1
      if ( em1 == 0 ) then
        em1 = 3
      end if

      ep1 = e + 1
      if ( ep1 == 4 ) then
        ep1 = 1
      end if

      fm1 = f - 1
      if ( fm1 == 0 ) then
        fm1 = 3
      end if

      fp1 = f + 1
      if (fp1 == 4) then
        fp1 = 1
      end if

      til(ep1,t) = c
      til(fp1,u) = i
      r = tnbr(ep1,t)
      s = tnbr(fp1,u)
      tnbr(ep1,t) = u
      tnbr(fp1,u) = t
      tnbr(e,t) = s
      tnbr(f,u) = r

      if ( tnbr(fm1,u) > 0 ) then
        top = top + 1
        stack(top) = u
      end if

      if ( s > 0 ) then

        if ( tnbr(1,s) == u ) then
          tnbr(1,s) = t
        else if ( tnbr(2,s) == u ) then
          tnbr(2,s) = t
        else
          tnbr(3,s) = t
        end if

        top = top + 1

        if ( top > maxst ) then
          ierr = 8
          return
        end if

        stack(top) = t

      else

        if ( u == btri .and. fp1 == bedg ) then
          btri = t
          bedg = e
        end if

        l = -(3*t + e-1)
        tt = t
        ee = em1

        do while ( tnbr(ee,tt) > 0 )

          tt = tnbr(ee,tt)

          if ( til(1,tt) == a ) then
            ee = 3
          else if ( til(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tnbr(ee,tt) = l

      end if

      if ( r > 0 ) then

        if ( tnbr(1,r) == t ) then
          tnbr(1,r) = u
        else if ( tnbr(2,r) == t ) then
          tnbr(2,r) = u
        else
          tnbr(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( tnbr(ee,tt) > 0 )

          tt = tnbr(ee,tt)

          if ( til(1,tt) == b ) then
            ee = 3
          else if ( til(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tnbr(ee,tt) = l

      end if

    end if

  end do

  return
end
subroutine sxy ( ival, val_max, ndat, nfine, nval, xdat, xval, ydat, yval )
!
!***********************************************************************
!
!! SXY sets up a cubic spline ( X, S(X) ), interpolating given data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer IVAL(VAL_MAX).  IVAL(I) is 1 if point I
!    is to be connected to point I+1.
!
!    Input, integer VAL_MAX, the maximum number of points that can be 
!    used for plots.
!
!    Input/output, integer NVAL, the current number of plot points.
!
  implicit none
!
  integer val_max
  integer ndat
!
  integer i
  integer ibcbeg
  integer ibcend
  integer iderv
  integer ival(val_max)
  integer nfine
  integer nval
  logical s_eqi
  real xdat(ndat)
  real xhi
  real xlo
  real xval(val_max)
  real y0
  real y1
  real y2
  real ybcbeg
  real ybcend
  real ydat(ndat)
  real ypp(ndat)
  real yval(val_max)
!
!  Get the boundary conditions
!
  ibcbeg = 0
  ybcbeg = 0.0E+00

  do

    write ( *, '(a)' ) 'Enter left boundary condition code.'
    write ( *, '(a)' ) '0=not-a-knot, 1=first derivative, 2=second.'
    read ( *, * ) ibcbeg

    if ( 0 <= ibcbeg .and. ibcbeg <= 2 ) then
      exit
    end if

    write ( *, '(a)' ) 'That was not acceptable.'

  end do

  if ( ibcbeg == 1 ) then
    write ( *, '(a)' ) 'Enter first derivative value:'
    read ( *, * ) ybcbeg
  else if ( ibcbeg == 2 ) then
    write ( *, '(a)' ) 'Enter second derivative value:'
    read ( *, * ) ybcbeg
  end if

  ibcend = 0
  ybcend = 0.0E+00

  do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'specify right boundary.'
    write ( *, '(a)' ) '0=not-a-knot, 1=first derivative, 2=second.'
    read ( *, * ) ibcend

    if ( 0 <= ibcend  .and. ibcend <= 2 ) then
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'That was not acceptable.'

  end do

  if ( ibcend == 1 ) then
    write ( *, '(a)' ) 'Enter first derivative value:'
    read ( *, * ) ybcend
  else if ( ibcend == 2 ) then
    write ( *, '(a)' ) 'Enter second derivative value:'
    read ( *, * ) ybcend
  end if
!
!  Set up the cubic spline.
!
  call spline_cubic_set ( ndat, xdat, ydat, ibcbeg, ybcbeg, ibcend, ybcend, ypp )
!
!  Fill in the graph values, NVAL, XVAL, YVAL.
!
  write ( *, '(a)' ) 'Enter 0 to graph spline, 1 for derivative, etc.'
  read ( *, * ) iderv

  xlo = xdat(1)
  xhi = xdat(ndat)

  do i = 1, nfine

    nval = nval+1
    call rvec_even_select ( xlo, xhi, nfine, i, xval(nval) )
    call spline_cubic_val ( ndat, xdat, ydat, ypp, xval(i), y0, y1, y2 )

    if ( iderv == 0 ) then
      yval(i) = y0
    else if ( iderv == 1 ) then
      yval(i) = y1
    else if ( iderv == 2 ) then
      yval(i) = y2
    else
      yval(i) = 0.0E+00
    end if

    ival(nval) = 1

  end do

  ival(nval) = 0

  return
end
subroutine symadd ( ierror, ifree, iopsym, iprsym, ipval, &
  maxsym, maxval, namvr, nsym, nsymp, symbol, valsym )
!
!***********************************************************************
!
!! SYMADD adds a symbol name to the list of symbolic names.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred, variable was added to list.
!    1, error occurred, variable was not added to list.
!
!    Input/output, integer IFREE, the index of the next free
!    memory location in VALSYM.
!
!    Input, integer IOPSYM(MAXSYM), contains, for each
!    symbol, the number of operands.  In particular, if
!    a symbol represents a constant, IOPSYM(I) is 0.
!    If a symbol represents a unary operator such as ABS,
!    IOPSYM(I) is 1.  IOPSYM may be as large as 3.
!
!    Input, integer IPRSYM(MAXSYM), the relative priorities
!    of the functions.
!
!    Input, integer IPVAL(MAXSYM), contains, for each symbol,
!    the address in VALSYM where the value of the symbol
!    is stored, if it is a scalar.  If the symbol represents
!    a vector or matrix, IPVAL(IARG) points to location of
!    the first entry.
!
!    Input, integer MAXSYM, the maximum number of symbols allowed.
!
!    Input, integer MAXVAL, the maximum number of values allowed
!    in VALSYM.
!
!    Input, character ( len = MAXLEN ) NAMVR, the name of the variable.
!
!    Input, integer NSYM, the total number of symbols, including temporaries.
!
!    Input, integer NSYMP, the number of permanent symbols.
!
!    Input, integer NUMDIM(2,MAXSYM).
!    For each symbol I, NUMDIM(1,I) is the number of rows in its
!    value, and NUMDIM(2,I) is the number of columns.
!
!    Input, character ( len = MAXLEN ) SYMBOL(MAXSYM), the symbolic names
!    of the internally defined constants and functions.
!
!    Input, real VALSYM(MAXVAL), contains the values of all the
!    symbolic variables.
!
  implicit none
!
  integer, parameter :: maxlen = 20
  integer maxsym
  integer maxval
!
  integer i
  integer ierror
  integer ifree
  integer indx
  integer iopsym(maxsym)
  integer iprsym(maxsym)
  integer ipval(maxsym)
  integer lennam
  character ( len = maxlen ) namvr
  integer nsym
  integer nsymp
  logical s_eqi
  character ( len = maxlen ) symbol(maxsym)
  real valsym(maxval)
!
!  Get the length of the variable name.
!
  ierror = 0
  lennam = len_trim ( namvr )

  if ( lennam <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SYMADD - Error!'
    write ( *, '(a)' ) '  This variable name has zero length!'
    ierror = 1
    return
  end if
!
!  Check if the name is already in use.
!
  do i = 1, nsymp

    if ( s_eqi ( namvr, symbol(i) ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SYMADD - Error!'
      write ( *, '(a)' ) '  The name ' // trim ( namvr ) // &
        ' is already in use.'
      ierror = 1
      return
    end if

  end do
!
!  Check for user overriding earlier use of same name.
!
  do i = nsymp+1, nsym

    if ( s_eqi ( namvr, symbol(i) ) ) then

      symbol(i) = 'VOID'

    end if

  end do
!
!  Is there enough space in SYMBOL for another symbol name?
!
  if ( nsym >= maxsym ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SYMADD - Error!'
    write ( *, '(a)' ) '  There is not enough memory for ' // trim ( namvr )
    write ( *, '(a)' ) '  Perhaps the KILL or INIT commands would help.'
    ierror = 1
    return
  end if
!
!  Is there enough space in VALSYM for the symbol's value?
!
  if ( ifree > maxval ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SYMADD - Error!'
    write ( *, '(a)' ) '  There is not enough memory for ' // trim ( namvr )
    write ( *, '(a)' ) '  Perhaps the KILL or INIT commands would help.'
    return
  end if
!
!  Insert information about the symbol into various tables.
!
  nsym = nsym + 1
  symbol(nsym) = namvr
  iopsym(nsym) = 0
  iprsym(nsym) = 10
  ipval(nsym) = ifree
  ifree = ifree + 1
  indx = ipval(nsym)
  valsym(indx) = 0.0E+00

  return
end
subroutine symbol_value ( com, ierror, ifree, iopsym, iprsym, ipval, maxsym, &
  maxval, namvar, nsym, nsymp, symbol, valsym, value )
!
!*******************************************************************************
!
!! SYMBOL_VALUE evaluates, sets or deletes a variable.
!
!  Discussion:
!
!    SYMVAL accepts the name of a symbol in NAMVAR, and does one
!    of three things, depending on the value of COM:
!
!      'R' - "reads" the value of the variable, returning it in VALUE.
!      'V' - "sets" the value of the variable from the array VALUE.
!      'D' - "deletes" the variable from memory.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character COM, determines what task is to be done:
!    'R', "reads" the variable, returning its value and dimensions.
!    'V', "sets" the value and dimensions of the variable.
!    'D', "deletes" the variable from memory.
!
!    Output, integer ierror, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input/output, integer IFREE, the index of the next free
!    memory location in VALSYM.
!
!    Input, integer IOPSYM(MAXSYM), contains, for each
!    symbol, the number of operands.  In particular, if
!    a symbol represents a constant, IOPSYM(i) is 0.
!    If a symbol represents a unary operator such as abs,
!    IOPSYM(i) is 1.  IOPSYM may be as large as 3.
!
!    Input, integer IPRSYM(MAXSYM), the function priorities.
!
!    Input, integer IPVAL(MAXSYM), contains, for each symbol,
!    the address in VALSYM where the value of the symbol
!    is stored, if it is a scalar.  If the symbol represents
!    a vector or matrix, IPVAL(IARG) points to location of
!    the first entry.
!
!    Input, integer MAXSYM, the maximum number of symbols allowed.
!
!    Input, integer MAXVAL, the maximum number of values allowed.
!
!    Input, character ( len = maxlen ) namvar, the name of the variable.
!
!    Input, integer NSYM, the total number of symbols, including temporaries.
!
!    Input, integer NSYMP, the number of permanent symbols.
!
!    Input, character ( len = maxlen ) SYMBOL(MAXSYM), symbolic names.
!
!    Input, real VALSYM(MAXVAL), the values of the symbolic variables.
!
!    Input/output, real VALUE, holds the value of the variable.
!
  implicit none
!
  integer, parameter :: maxlen = 10
  integer maxsym
  integer maxval
!
  character com
  integer i
  integer ierror
  integer ifree
  integer ihi
  integer ilo
  integer indx
  integer iopsym(maxsym)
  integer iprsym(maxsym)
  integer ipval(maxsym)
  integer j
  integer lennam
  integer match
  character ( len = maxlen ) namvar
  integer ncol
  integer nrow
  integer nsym
  integer nsymp
  logical s_eqi
  character ( len = maxlen ) symbol(maxsym)
  real valsym(maxval)
  real value
!
  ierror = 0

  if ( len_trim ( namvar ) <= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SYMBOL_VALUE - Error!'
    write ( *, '(a)' ) '  Bad variable name ' // trim ( namvar )
    return
  end if
!
!  Search for a match between NAMVAR and any predefined symbol.
!
  match = 0

  do i = 1, nsym
    if ( s_eqi ( namvar, symbol(i) ) ) then
      match = i
      exit
    end if
  end do
!
!  If no such variable seems to exist, but we have been asked
!  to DELETE or READ such a variable, then exit.
!
  if ( match == 0 ) then

    if ( s_eqi ( com, 'D' ) .or. s_eqi ( com, 'R' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SYMBOL_VALUE - Error!'
      write ( *, '(a)' ) '  No variable named ' // trim ( namvar )

      if ( s_eqi ( com, 'R' ) ) then
        ierror = 1
        return
      end if

    end if
!
!  Add the name of the new variable to the internal list.
!
    call symadd ( ierror, ifree, iopsym, iprsym, ipval, &
      maxsym, maxval, namvar, nsym, nsymp, symbol, valsym )

    if ( ierror /= 0 ) then
      return
    end if

    match = nsym
!
!  The matched symbol has index MATCH.
!
  else
!
!  If the symbol is not the name of a variable, but
!  rather the name of an operator, then this is an error.
!
    if ( iopsym(match) /= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SYMBOL_VALUE - Error!'
      write ( *, '(a)' ) '  You are trying to assign a value to'
      write ( *, '(a)' ) '  the name of a function or operator:'
      write ( *, '(a)' ) symbol(match)
      return
    end if
!
!  Some variables may not be revalued.
!
    if ( s_eqi ( com, 'V' ) ) then
      if ( s_eqi ( symbol(match), 'E' ) ) ierror = 1
      if ( s_eqi ( symbol(match), 'EPS' ) ) ierror = 1
      if ( s_eqi ( symbol(match), 'PI' ) ) ierror = 1
      if ( symbol(match)(1:1) == '"' ) ierror = 1
    end if
!
!  Some variables may not be deleted.
!
    if ( s_eqi ( com, 'D' ) ) then
      if ( s_eqi ( symbol(match), 'E' ) ) ierror = 1
      if ( s_eqi ( symbol(match), 'EPS' ) ) ierror = 1
      if ( s_eqi ( symbol(match), 'PI' ) ) ierror = 1
      if ( s_eqi ( symbol(match), 'RTOL' ) ) ierror = 1
      if ( s_eqi ( symbol(match), 'SEED' ) ) ierror = 1
      if ( symbol(match)(1:1) == '"' ) ierror = 1
    end if

    if ( ierror /= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SYMBOL_VALUE - Error!'
      write ( *, '(a)' ) '  Try to change special variable:'
      write ( *, '(a)' ) symbol(match)
      return
    end if

    nrow = 1
    ncol = 1

  end if

  indx = ipval(match)

  if ( s_eqi ( com, 'V' ) ) then
    valsym(indx) = value
  else if ( s_eqi ( com, 'R' ) ) then
    value = valsym(indx)
  else if ( s_eqi ( com, 'D' ) ) then

    ilo = indx
    ihi = ifree-1
    do i = 1, ihi+1-ilo
      valsym(ilo-nrow*ncol+i-1) = valsym(ilo+i-1)
    end do
    ifree = ifree-nrow*ncol

    do i = match+1, nsym
      iopsym(i-1) = iopsym(i)
      iprsym(i-1) = iprsym(i)
      ipval(i-1) = ipval(i)-nrow*ncol
      symbol(i-1) = symbol(i)
    end do

    nsym = nsym-1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SYMBOL_VALUE:'
    write ( *, '(a)' ) '  Removed the variable ' // trim ( namvar )
    write ( *, '(a,i6)' ) '  Free memory now: ', maxval+1-ifree

  end if

  return
end
subroutine timestamp ( )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
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
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tokens ( ierror, ifinis, ifree, implic, indx1, indx2, ineg, infix, &
  intsym, iopsym, iprsym, ipval, maxsym, maxval, nints, nsym, nsymp, &
  nuvar, symbol, valsym )
!
!*******************************************************************************
!
!! TOKENS parses a character string into recognized symbols and constants.
!
!  Discussion:
!
!    The routine produces the output array of integers INTSYM, containing
!    whose entries stand for recognized symbols, variable names or
!    constants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer IFINIS, the index in SYMBOL of the symbol
!    '$', meaning the end of the formula.
!
!    Input, integer IFREE, the address of the next free
!    memory location in VALSYM.
!
!    Output, integer IMPLIC, implicit variable flag.
!    0, the assignment variable was not implicit.
!    nonzero, the assignment variable was implicit, and has been
!    added to the symbol table, and assigned the symbol index of IMPLIC.
!
!    Input, integer INDX1, the index of the symbol 'INDEX1'.
!
!    Input, integer INDX2, the index of the symbol 'INDEX2'.
!
!    Input, integer INEG, the index of the symbol 'NEG'.
!
!    Input, character ( len = * ) INFIX, a sequence of characters representing
!    an infix formula.
!
!    Output, integer INTSYM(80), a sequence of integers which are
!    the indices of the symbols used in INFIX.
!
!    Input, integer IOPSYM(MAXSYM), contains, for each
!    symbol, the number of operands.  In particular, if
!    a symbol represents a constant, IOPSYM(I) is 0.
!    If a symbol represents a unary operator such as ABS,
!    IOPSYM(I) is 1.  IOPSYM may be as large as 3.
!
!    Input, integer IPRSYM(MAXSYM), the relative priorities
!    of the functions.
!
!    Input, integer IPVAL(MAXSYM), contains, for each symbol,
!    the address in VALSYM where the value of the symbol
!    is stored, if it is a scalar.  If the symbol represents
!    a vector or matrix, IPVAL(IARG) points to location of
!    the first entry.
!
!    Input, integer MAXSYM, the maximum number of symbols allowed.
!
!    Input, integer MAXVAL, the maximum number of values allowed in VALSYM.
!
!    Output, integer NINTS, the number of useful entries in INTSYM.
!
!    Input, integer NSYM, the total number of symbols, including temporaries.
!
!    Input, integer NSYMP, the number of permanent symbols.
!
!    Input, integer NUVAR.
!    If NUVAR is 0, then the formula may refer to variables which
!    have not yet been declared.
!    Otherwise, all variables on the right hand side of the "="
!    sign must have been previously declared.
!
!    Input, character ( len = MAXLEN ) SYMBOL(MAXSYM), the symbolic names
!    of the internally defined constants and functions.
!
!    Input, real VALSYM(MAXVAL), contains the values of all the
!    symbolic variables.
!
  implicit none
!
  integer, parameter :: maxlen = 20
  integer maxsym
!
  character ctemp
  integer i
  integer iback
  integer icom1
  integer icom2
  integer ierror
  integer ifinis
  integer ifree
  integer ilpr
  integer implic
  integer indx1
  integer indx2
  integer ineg
  character ( len = * ) infix
  integer intsym(80)
  integer iopsym(maxsym)
  integer ipoint
  integer ipos
  integer iprsym(maxsym)
  integer ipval(maxsym)
  integer irpr
  integer j
  integer jback
  integer lchar
  integer lenfix
  integer lenmat
  integer lens
  integer loc
  integer loc1
  integer loc2
  integer match
  integer matcho
  integer maxval
  character ( len = maxlen ) namvr
  integer nints
  logical notnum
  integer nsym
  integer nsymp
  integer nuvar
  real rval
  logical s_eqi
  logical s_is_alpha
  character sym
  character ( len = maxlen ) sym1
  character ( len = maxlen ) sym2
  character ( len = maxlen ) sym3
  character ( len = maxlen ) symbol(maxsym)
  real valsym(maxval)
!
  implic = 0
  ierror = 0
  nints = 0
  match = 0
  ipos = 1
  lenfix = len_trim ( infix )
!
!  1: Does this formula begin with a variable name followed
!     by an equals sign?
!
  sym = '='
  loc1 = index ( infix, sym )
  if ( loc1 <= 1 ) then
    go to 20
  end if
!
!  The formula begins with a variable name followed by an equals sign.
!
!  Pick off the name of the variable, but first make sure
!  you stop before any left hand parenthesis, in case the
!  formula is something like "alfred(7)=5".
!
  sym = '('
  loc2 = index ( infix, sym )

  if ( loc2 <= 0 ) then
    loc = loc1
  else
    loc = min ( loc1, loc2 )
  end if

  if ( loc <=1 .or. loc > maxlen+1 ) then
    go to 20
  end if

  loc = loc - 1
  namvr = infix(1:loc)
!
!  Check to see whether the name of the assignment variable has
!  already been defined.
!
  do i = 1, nsym
    if ( s_eqi ( namvr, symbol(i) ) ) then
      go to 20
    end if
  end do
!
!  The assignment variable has not previously been defined.
!  Prepare to add it to the table of symbols.
!
  ipos = loc + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOKENS:'
  write ( *, '(a)' ) '  Implicit definition of ' // trim ( namvr )

  call symadd ( ierror, ifree, iopsym, iprsym, ipval, maxsym, maxval, namvr, &
    nsym, nsymp, symbol, valsym )

  if ( ierror /= 0 ) then
    return
  end if

  nints = nints + 1
  intsym(nints) = nsym
  match = nsym
  implic = nsym
!
!  Process the next chunk of the formula.
!
20    continue
!
!  If we have reached the end of the input string, we're done.
!  Append an "End-of-formula" marker "$" and return.
!
  if ( ipos > lenfix ) then
    nints = nints + 1
    intsym(nints) = ifinis
    return
  end if
!
!  We're not done yet.  Remember the last match in MATCHO, and
!  look for the next match.  Go through all the symbols, and
!  try to find the longest match possible.
!
  matcho = match
  match = 0
  lenmat = 0

  do i = 1, nsym

    lens = len_trim ( symbol(i) )

    if ( lens > lenmat .and. ipos + lens - 1 <= lenfix  ) then

      if ( s_eqi ( infix(ipos:ipos+lens-1), symbol(i) ) ) then
        lenmat = lens
        match = i
      end if

    end if

  end do

  if ( match == 0 ) then
    go to 110
  end if
!
!  We found a match.
!
!  But watch out!  The user might be implicitly defining a name,
!  but because the first part matches an earlier symbol, you don't
!  realize it!  For example, in a formula involving SINGLE, you
!  might think you see "SIN".  So make sure that the match you've
!  got is not followed by an alphabetic character.
!
  if ( iopsym(match) /= 0 ) then
    go to 40
  end if

  if ( ipos+lenmat-1 >= lenfix ) then
    go to 40
  end if

  ctemp = infix(ipos+lenmat:ipos+lenmat)

  if ( s_is_alpha ( ctemp ) ) then
    go to 140
  end if
!
!  We have matched an old symbol.
!
40    continue
  sym1 = symbol(match)
  sym2 = ' '
  if ( matcho /= 0 ) then
    sym2 = symbol(matcho)
  end if
!
!  Check for unary minus or plus.
!
  if ( sym1 == '-' ) then
    if ( sym2 == '**' .or. sym2 == ',' .or. &
      sym2 == '(' .or. sym2 == '=' .or. &
         ipos == 1 ) then

      sym = infix(ipos+1:ipos+1)

      if ( lge ( sym, '0' ) .and. lle ( sym, '9' ) ) then
        go to 110
      end if

      if ( sym == '.' ) then
        go to 110
      end if

    end if
  end if

  if ( sym1 == '-' ) then

    if ( ipos == 1 .or. sym2 == '(' .or. sym2 == ',' .or. sym2 == '=' ) then
      nints = nints + 1
      intsym(nints) = ineg
      ipos = ipos + 1
      go to 20
    end if

  end if

  if ( sym1 == '+' ) then

    if ( ipos == 1 .or. sym2 == '(' .or. sym2 == '=' ) then
      ipos = ipos + 1
      go to 20
    end if

  end if
!
!  If we've hit a right parenthesis, then we assume we're dealing with:
!
!    a vector index: A(3)
!    a matrix index: A(1,2)
!    a two place operator: MAX(A,B), or MIN or POLY
!    a three place operator: KRYLOV(A,X,N)
!
!  We need to find the matching left parenthesis, count the
!  arguments by assuming that commas separate them, and read
!  off the name preceding the left parenthesis.
!
  if ( sym1 == ')' ) then

    icom1 = 0
    icom2 = 0
    ilpr = 0
    irpr = 1

    do iback = 1, nints

      i = nints + 1 - iback
      if ( intsym(i) <= 0 ) then
        go to 90
      end if

      sym3 = symbol ( intsym(i) )

      if ( sym3 == ',' .and. irpr-ilpr == 1 ) then

        if ( icom1 == 0 ) then
          icom1 = i
        else if ( icom2 == 0 ) then
          icom2 = i
        else
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TOKENS - Error!'
          write ( *, '(a)' ) '  Too many commas in formula.'
          return
        end if

      else if ( sym3 == '(' ) then
        ilpr = ilpr + 1
      else if ( sym3 == ')' ) then
        irpr = irpr + 1
      end if
!
!  If ILPR=IRPR, we've reached the matching left parenthesis.
!
      if ( ilpr == irpr ) then

        if ( i <= 1 ) then
          go to 100
        end if

        if ( intsym(i-1) <= 0 ) then
          go to 100
        end if

        sym3 = symbol ( intsym(i-1) )

        if ( iopsym(intsym(i-1)) /= 0 .and. &
             .not. s_eqi ( sym3, 'MAX' ) .and. &
             .not. s_eqi ( sym3, 'MIN' ) .and. &
             .not. s_eqi ( sym3, 'ATAN2' ) .and. &
             .not. s_eqi ( sym3, 'POLY' ) .and. &
             .not. s_eqi ( sym3, 'KRYLOV' ) ) then
          go to 100
        end if
!
!  Copy ) into INTSYM.
!
        nints = nints + 1
        intsym(nints) = match
!
!  INDX1 is our name for a vector reference.
!
        if ( icom1 == 0 ) then

          match = indx1
!
!  INDX2 is our name for a matrix reference.
!
        else if ( icom2 == 0 ) then

          intsym(icom1) = match

          do jback = 1, nints-icom1
            j = nints + 1 - jback
            intsym(j+1) = intsym(j)
          end do

          nints = nints + 1
          intsym(icom1+1) = match - 1

          if ( iopsym(intsym(i-1)) /= 0 ) then
            ipos = ipos + 1
            go to 20
          end if

          match = indx2
!
!  Insert )( to divide a trio of arguments.
!
        else

          intsym(icom1) = match

          do jback = 1, nints-icom1
            j = nints + 1 - jback
            intsym(j+1) = intsym(j)
          end do

          intsym(icom1+1) = match - 1

          nints = nints + 1

          intsym(icom2) = match

          do jback = 1, nints-icom2
            j = nints + 1 - jback
            intsym(j+1) = intsym(j)
          end do

          intsym(icom2+1) = match - 1

          nints = nints + 1

          if ( iopsym(intsym(i-1)) /= 0 ) then
            ipos = ipos + 1
            go to 20
          end if

          match = indx2

        end if

        ipos = ipos + 1
        nints = nints + 1
        intsym(nints) = match
        go to 20
      end if

90    continue

    end do

  end if

100   continue

  lens = len_trim ( symbol(match) )
  nints = nints + 1
  intsym(nints) = match
  ipos = ipos + lens
  go to 20
!
!  Check for a constant.
!
110   continue

  call s_to_r ( infix(ipos:), rval, ierror, lchar )

  if ( lchar <= 0 ) then
    go to 140
  end if

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TOKENS - Error!'
    write ( *, '(a)' ) '  Illegal real number.'
    return
  end if
!
!  Look for a constant already defined with the same value.
!
  do i = nsymp+1, nsym

    if ( iprsym(i) == 10 ) then
      if ( symbol(i)(1:1) == '"' ) then
        ipoint = ipval(i)
        if ( valsym(ipoint) == rval ) then
          match = i
          go to 130
        end if
      end if
    end if

  end do
!
!  This is a new constant.  Make up a name for it.
!
  if ( real ( int ( rval ) ) == rval ) then
    write ( namvr, '(''"'',i19)' ) int ( rval )
  else
    write ( namvr, '(''"'', g19.10)' ) rval
  end if

  call s_blank_delete ( namvr )
!
!  Add the constant to the symbol list.
!
  call symadd ( ierror, ifree, iopsym, iprsym, ipval, maxsym, maxval, namvr, &
    nsym, nsymp, symbol, valsym )

  if ( ierror /= 0 ) then
    return
  end if

  ipoint = ipval(nsym)
  valsym(ipoint) = rval
  match = nsym

130   continue
!
!  Advance LCHAR characters, except that CHRCTR will read in a
!  comma as the end of the number, and we have to back up in
!  such a case.
!
  if ( infix(ipos+lchar-1:ipos+lchar-1) /= ',' ) then
    ipos = ipos + lchar
  else
    ipos = ipos + lchar - 1
  end if

  nints = nints + 1
  intsym(nints) = match
  go to 20
!
!  Consider implicit variable declaration on right hand side.
!
140   continue

  if ( nuvar == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TOKENS - Error!'
    write ( *, '(a)' ) '  Undeclared variable in formula.'
    return
  end if

  namvr = ' '

  do i = 1, maxlen
    sym = infix(ipos-1+i:ipos-1+i)
    notnum = llt ( sym, '0' ) .or. lgt ( sym, '0' )

    if ( ( .not. s_is_alpha(sym) ) .and. notnum ) then
      exit
    end if

    namvr(i:i) = sym

  end do

  lchar = len_trim ( namvr )

  if ( lchar <= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TOKENS - Error!'
    write ( *, '(a)' ) '  The string starting at location ', ipos, &
    ' is unreadable.'
    return
  end if

  ipos = ipos + lchar

  call symadd ( ierror, ifree, iopsym, iprsym, ipval, maxsym, maxval, namvr, &
    nsym, nsymp, symbol, valsym )


  nints = nints + 1
  intsym(nints) = nsym
  match = nsym

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOKENS:'
  write ( *, '(a)' ) '  The undeclared variable ' // trim ( namvr ) // &
    ' will be assumed to be a scalar.'

  go to 20
end
subroutine tricon ( ival, val_max, ndat, nodtri, ntri, nval, udat, &
  uval, xdat, xval, ydat, yval )
!
!***********************************************************************
!
!! TRICON draws a contour line of a scalar quantity.
!
!  Discussion:
!
!    The scalar quantity is assumed to be known at the nodes of a
!    triangulation of the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer IVAL(VAL_MAX).  IVAL(I) is 1 if point I
!    is to be connected to point I+1.
!
!    Input, integer VAL_MAX, the maximum number of points that can be 
!    used for plots.
!
!    Input, integer NODTRI(3,NTRI).  Triangle I is made up of nodes
!    NODTRI(1,I), NODTRI(2,I) and NODTRI(3,I).
!
!    Input, integer NTRI, the number of triangles.
!
!    Input/output, integer NVAL, the current number of plot points.
!
!    Input, real UDAT(NDAT).  UDAT(I) is the function value at
!    the point ( XDAT(I), YDAT(I) ).
!
!    Input, real UVAL, the value associated with the contour line.
!
!    Input, real XDAT(NDAT), the X coordinates of the nodes.
!
!    Input, real YDAT(NDAT), the Y coordinates of the nodes.
!
  implicit none
!
  integer val_max
  integer ndat
  integer ntri
!
  integer i
  integer icross
  integer ival(val_max)
  integer j
  integer jp1
  integer n1
  integer n2
  integer nodtri(3,ntri)
  integer nval
  real t
  real u1
  real u2
  real udat(ndat)
  real uval
  real x1
  real x2
  real xdat(ndat)
  real xval(val_max)
  real xx(3)
  real y1
  real y2
  real ydat(ndat)
  real yval(val_max)
  real yy(3)
!
!  Consider each triangle:
!
  do i = 1, ntri

    icross = 0
!
!  Consider each side of the triangle:
!
    do j = 1, 3

      n1 = nodtri(j,i)
      x1 = xdat(n1)
      y1 = ydat(n1)
      u1 = udat(n1)

      jp1 = j+1
      if ( jp1 > 3 ) then
        jp1 = 1
      end if

      n2 = nodtri(jp1,i)
      x2 = xdat(n2)
      y2 = ydat(n2)
      u2 = udat(n2)

      call line_seg_contains_point_1d ( u1, u2, uval, t )

      if ( 0.0E+00 <= t .and. t <= 1.0E+00 ) then
        icross = icross+1
        xx(icross) = x1 + (uval-u1)*(x2-x1)/(u2-u1)
        yy(icross) = y1 + (uval-u1)*(y2-y1)/(u2-u1)
      end if

    end do
!
!  If we found two crossings, draw the line.
!
    if ( icross == 2 ) then

      nval = nval +1
      xval(nval) = xx(1)
      yval(nval) = yy(1)
      ival(nval) = 1

      nval = nval + 1
      xval(nval) = xx(2)
      yval(nval) = yy(2)
      ival(nval) = 0

    end if

  end do

  return
end
subroutine vbedg ( x, y, vcl, til, tnbr, ltri, ledg, rtri, redg )
!
!***********************************************************************
!
!! VBEDG determines which boundary edges are visible to a point.
!
!  Discussion:
!
!    The point (X,Y) is assumed to be outside the convex hull of the
!    region covered by the 2D triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Parameters:
!
!        X,Y - 2-D point outside convex hull.
!
!        VCL(1:2,1:*) - coordinates of 2-D vertices
!        TIL(1:3,1:*) - triangle incidence list
!        TNBR(1:3,1:*) - triangle neighbor list; negative values are
!              used for links of CCW linked list of boundary edges;
!              LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input/output, LTRI,LEDG - if LTRI /= 0 then they are assumed to be
!    as defined below and are not changed, else they are updated.
!    On output, LTRI is the index of boundary triangle to left of leftmost
!    boundary triangle visible from (X,Y)
!    LEDG is the boundary edge of triangle LTRI to left of leftmost
!    boundary edge visible from (X,Y).  1 <= LEDG <= 3.
!
!    Input/output, RTRI.  On input, the index of boundary triangle to
!    begin search at.  On output, the index of the rightmost boundary
!    triangle visible from (X,Y)
!
!    Input/output, REDG.  The edge of triangle RTRI that is visible from (X,Y)
!    1 <= REDG <= 3.
!
  implicit none
!
  integer a
  integer b
  integer e
  integer l
  logical ldone
  integer ledg
  integer lr
  integer lrline
  integer ltri
  integer redg
  integer rtri
  integer t
  real temp
  integer til(3,*)
  integer tnbr(3,*)
  real vcl(2,*)
  real x
  real y
!
!  Find rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor information.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

   10 continue

  l = -tnbr(redg,rtri)
  t = l/3
  e = mod(l,3) + 1
  a = til(e,t)

  if ( e <= 2 ) then
    b = til(e+1,t)
  else
    b = til(1,t)
  end if

  temp = 0.0E+00
  lr = lrline ( x, y, vcl(1,a), vcl(2,a), vcl(1,b), vcl(2,b), temp )

  if ( lr > 0 ) then
    rtri = t
    redg = e
    go to 10
  end if

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

   20 continue

  b = til(e,t)

  if ( e >= 2 ) then
    e = e - 1
  else
    e = 3
  end if

30    continue

  if ( tnbr(e,t) > 0 ) then

    t = tnbr(e,t)

    if ( til(1,t) == b ) then
      e = 3
    else if ( til(2,t) == b ) then
      e = 1
    else
      e = 2
    end if

    go to 30

  end if

  a = til(e,t)
  temp = 0.0E+00
  lr = lrline ( x, y, vcl(1,a), vcl(2,a), vcl(1,b), vcl(2,b), temp )

  if ( lr > 0 ) then
    go to 20
  end if

  ltri = t
  ledg = e

  return
end
subroutine vector_unit_nd ( n, v )
!
!*******************************************************************************
!
!! VECTOR_UNIT_ND normalizes a vector in ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, real V(N), the vector to be normalized.  On output,
!    V should have unit Euclidean norm.  However, if the input vector
!    has zero Euclidean norm, it is not altered.
!
  implicit none
!
  integer n
!
  real enorm_nd
  integer i
  real temp
  real v(n)
!
  temp = enorm_nd ( n, v )

  if ( temp /= 0.0E+00 ) then
    v(1:n) = v(1:n) / temp
  end if

  return
end
subroutine xys ( ival, maxdat, val_max, ndat, nfine, nval, tdat, xdat, &
  xval, ydat, yval )
!
!***********************************************************************
!
!! XYS computes splines through X and Y data.
!
!  Discussion:
!
!    An independent parameter T is used to trace out the curves.
!    T varies from 0 to 1.
!
!    The "not-a-knot" boundary condition is used at the endpoints of
!    the splines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer IVAL(VAL_MAX).  IVAL(I) is 1 if point I
!    is to be connected to point I+1.
!
!    Input, integer VAL_MAX, the maximum number of line segments that
!    can be drawn.
!
!    Input/output, integer NVAL, the current number of plot points.
!
  implicit none
!
  integer maxdat
  integer val_max
!
  integer i
  integer ibcbeg
  integer ibcend
  character isay
  integer ival(val_max)
  integer maxrwork
  integer ndat
  integer nfine
  integer nval
  real rwork(4,ndat)
  logical s_eqi
  real thi
  real tlo
  real tdat(maxdat)
  real tval
  real xbcbeg
  real xbcend
  real xdat(maxdat)
  real xpp(maxdat)
  real xppval
  real xpval
  real xval(val_max)
  real ybcbeg
  real ybcend
  real ydat(maxdat)
  real ypp(maxdat)
  real yppval
  real ypval
  real yval(val_max)
!
  tlo = 0.0E+00
  thi = 1.0E+00
!
!  Assign the T values.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter "Y" to join curve ends.'

  read ( *, '(a)' ) isay

  if ( s_eqi ( isay, 'Y' ) ) then
    ndat = ndat + 1
    xdat(ndat) = xdat(1)
    ydat(ndat) = ydat(1)
  end if

  call rvec_even ( tlo, thi, ndat, tdat)
!
!  Set up the spline data
!
  nval = nfine

  ibcbeg = 0
  xbcbeg = 0.0E+00
  ibcend = 0
  xbcend = 0.0E+00

  call spline_cubic_set ( ndat, tdat, xdat, ibcbeg, xbcbeg, ibcend, xbcend, &
    xpp )

  ibcbeg = 0
  ybcbeg = 0.0E+00
  ibcend = 0
  ybcend = 0.0E+00

  call spline_cubic_set ( ndat, tdat, ydat, ibcbeg, ybcbeg, ibcend, ybcend, &
    ypp )

  do i = 1, nval
    call rvec_even_select ( tlo, thi, nval, i, tval )
    call spline_cubic_val ( ndat, tdat, xdat, xpp, tval, xval(i), xpval, &
      xppval )
    call spline_cubic_val ( ndat, tdat, ydat, ypp, tval, yval(i), ypval, &
      yppval )
  end do

  ival(1:nval-1) = 1
  ival(nval) = 0

  return
end
subroutine xyz_write ( filexyz, ival, npnt, x, y, z )
!
!***********************************************************************
!
!! XYZ_WRITE writes graphics data to an XYZ file.
!
!  Discussion:
!
!    The XYZ format is very simple:
!
!    Comment lines begin with a '#' in column 1.
!    All other lines are either XYZ coordinates, or blank.
!    Two successive XYZ coordinates represent a line segment.
!    A blank line terminates the current line segment.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 80 ) FILEXYZ, the name of the XYZ file being
!    created.
!
!    Input, integer IVAL(MAXPNT), describes lines.  If IVAL(I) is:
!    0, no line goes through point I;
!    1, point I is the last point in a line;
!    2, point I is the first point in a line;
!    3, point I is an interior point of a line.
!
!    Input, integer NPNT, the number of points used in the internal data
!    format.
!
!    Input, real X(NPNT), Y(NPNT), Z(NPNT), contains the X, Y, and Z
!    components of points used in the internal data format.
!
  implicit none

  integer npnt

  character ( len = * ) filexyz
  integer iunit
  integer ival(npnt)
  integer j
  integer lenc
  integer ntline
  real x(npnt)
  real y(npnt)
  real z(npnt)

  call get_unit ( iunit )

  open ( unit = iunit, file = filexyz, status = 'replace' )

  write ( iunit, '(a)' ) '# ' // trim ( filexyz ) // ' created by XYZ_WRITE.'
  write ( iunit, '(a)' ) '#'
  ntline = 3

  do j = 1, npnt

    if ( j > 1 ) then
      if ( ival(j) == 0 .or. ival(j) == 2 ) then
        write ( iunit, '(a)' ) ' '
        ntline = ntline + 1
      end if
    end if

    write ( iunit, '(3g14.6)' ) x(j), y(j), z(j)
    ntline = ntline + 1

  end do

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XYZ_WRITE:'
  write ( *, '(a,i6,a)' ) '  Output file contains ', ntline, ' lines.'

  return
end
