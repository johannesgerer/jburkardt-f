program main

!*****************************************************************************80
!
!! MAIN is the main program for CRYSTAL_PLOT.
!
!  Discussion:
!
!    CRYSTAL_PLOT creats plots from the output of the CRYSTAL program.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxbot = 20
  integer ( kind = 4 ), parameter :: maxl = 64
  integer ( kind = 4 ), parameter :: maxobj = 46

  integer ( kind = 4 ), parameter :: maxm = 2*maxl-1

  real b1jbl(maxm)
  real b2jbl(maxm)
  real b3jbl(maxl)
  character ( len = 40 ) command
  real cost
  real cost2
  character ( len = 8 ) date
  real delx
  real dely
  character ( len = 10 ) dev
  real dxcdp(maxl,maxm)
  real dxdp(maxl,maxm)
  real dycdp(maxl,maxm)
  real dydp(maxl,maxm)
  real e(maxl,maxm)
  logical echo
  character ( len = 80 ) fildat
  character ( len = 80 ) filgrf
  character ( len = 80 ) filinp
  real gamt(maxl,maxm)
  real grace
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icmax
  integer ( kind = 4 ) icmin
  integer ( kind = 4 ) icolor(maxobj)
  integer ( kind = 4 ) icrys1
  integer ( kind = 4 ) icrys2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) itable
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) jcmax
  integer ( kind = 4 ) jcmin
  integer ( kind = 4 ) jcrys1
  integer ( kind = 4 ) jcrys2
  integer ( kind = 4 ) jtemp
  integer ( kind = 4 ) kcrys(maxl,maxm)
  integer ( kind = 4 ) kmelt(maxl,maxm)
  integer ( kind = 4 ) kvoid(maxl,maxm)
  integer ( kind = 4 ) l
  logical lbar
  integer ( kind = 4 ) lenc
  logical lflag(maxl,maxm)
  logical lgopen
  logical lnei
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nbot
  logical ncflag(maxl,maxm)
  integer ( kind = 4 ) ncon
  logical nflag(maxl,maxm)
  logical npflag(maxl,maxm)
  integer ( kind = 4 ) nxskip
  integer ( kind = 4 ) nyskip
  character ( len = 25 ) object(maxobj)
  logical ovrlay
  real p(maxl,maxm)
  real pc(maxl,maxm)
  real psi(maxl,maxm)
  logical reflec
  real rueta(maxl,maxm)
  real ruksi(maxl,maxm)
  logical s_eqi
  real s1(maxl,maxm)
  real s2(maxl,maxm)
  real scalecv
  real scalenb
  real scalenc
  real scalenp
  real scalev
  logical show(maxobj)
  real srange
  real t(maxl,maxm)
  real te(maxl,maxm)
  real temp
  real tempa(maxl,maxm)
  character ( len = 10 ) time
  character ( len = 80 ) title
  character ( len = 80 ) title2
  real tk(maxl,maxm)
  real tnow
  real tnow2
  real u(maxl,maxm)
  real v(maxl,maxm)
  logical vpflag(maxl,maxm)
  real vmag(maxl,maxm)
  real vort(maxl,maxm)
  real w(maxl,maxm)
  real x1max
  real x1min
  real x2max
  real x2min
  real xbot(maxbot)
  real xc(maxl,maxm)
  real xmax
  real xmin
  real xp(maxl,maxm)
  real xsmax
  real xsmin
  real xtmax
  real xtmin
  real y1max
  real y1min
  real y2max
  real y2min
  real ybot(maxbot)
  real yc(maxl,maxm)
  real ymax
  real ymin
  real yp(maxl,maxm)
  real ysmax
  real ysmin
  real ytmax
  real ytmin

  call timestamp ( )

  call hello
!
!  Initialize data.
!
  call init(b1jbl,b2jbl,b3jbl,cost,cost2,delx,dely,dev,echo,fildat,filgrf, &
    filinp,grace,icmax,icmin,icolor,icrys1,icrys2,iplot,itable,jcmax, &
    jcmin,jcrys1,jcrys2,l,lbar,m,maxbot,maxl,maxm,maxobj,nbot,ncflag, &
    ncon,nflag,npflag,nxskip,nyskip,object,ovrlay,p,pc,psi,reflec,scalecv, &
    scalenb,scalenc,scalenp,scalev,show,title,title2,u,v,vpflag,vmag,vort, &
    x1max,x1min,x2max,x2min,xbot,xc,xmax,xmin,xp,xsmax,xsmin,xtmax,xtmin, &
    y1max,y1min,y2max,y2min,ybot,yc,ymax,ymin,yp,ysmax,ysmin,ytmax,ytmin)

  open ( unit = 17, file = filinp, status = 'replace' )
!
!  Get the next command.
!
10    continue

  write ( *, * ) ' '
  write ( *, * ) '? ("H" for help)'
  read(*,'(a)',end = 50,err=50)command
  call flushl(command)
  write(17,'(a)') trim ( command )
  if ( echo ) then
    write(*,'(a)') trim ( command )
  end if
  if ( command == ' ')go to 10

11 continue
!
!  AXis: show axis of symmetry.
!
  if ( s_eqi ( command(1:2),'ax') ) then

    show(28) = .not.show(28)
    if ( show(28) ) then
      write ( *, '(a)' ) 'The axis of symmetry will be shown.'
    else
      write ( *, '(a)' ) 'The axis of symmetry will NOT be shown.'
    end if
!
!  BOundary: show boundary.
!
  else if ( s_eqi ( command(1:2),'bo') ) then

    show(1) = .not.show(1)
    if ( show(1) ) then
      write ( *, '(a)' ) 'The boundary will be shown.'
    else
      write ( *, '(a)' ) 'The boundary will NOT be shown.'
    end if
!
!  BAr: Switch display of color bar.
!
  else if ( s_eqi ( command(1:2),'ba') ) then
 
    lbar = .not.lbar
    if ( lbar ) then
      write ( *, '(a)' ) 'The color bar will be shown.'
    else
      write ( *, '(a)' ) 'The color bar will NOT be shown.'
    end if
!
!  BH: bottom half.
!
  else if ( s_eqi ( command,'bh') ) then

    ysmax = ysmin+0.5*(ysmax-ysmin)

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  BL: bottom left quarter.
!
  else if ( s_eqi ( command,'bl') ) then

    xsmax = xsmin+0.5*(xsmax-xsmin)
    ysmax = ysmin+0.5*(ysmax-ysmin)

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  BC: bottom center quarter.
!
  else if ( s_eqi ( command,'bc') ) then

    temp = 0.25*(xsmax-xsmin)
    xsmin = xsmin+temp
    xsmax = xsmax-temp
    ysmax = ysmin+0.5*(ysmax-ysmin)

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  BR: bottom right quarter.
!
  else if ( s_eqi ( command,'br') ) then

    xsmin = xsmin+0.5*(xsmax-xsmin)
    ysmax = ysmin+0.5*(ysmax-ysmin)

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  C: choose colors.
!
  else if ( s_eqi ( command,'c') ) then

    write ( *, * ) ' '
    write ( *, * ) ' Number  Color  Name'
    write ( *, * ) ' '
    do i = 1,maxobj
      write(*,'(1x,i2,2x,i3,2x,a)')i,icolor(i),object(i)
    end do

    write ( *, * ) ' '
    write ( *, * ) 'Enter an object number, and a color number.'
    read(*,*,end = 50,err=50)itemp,jtemp
    write(17,*)itemp,jtemp
    if ( echo ) then
      write ( *, * ) itemp,jtemp
    end if

    if ( 1 <= itemp .and. itemp <= maxobj ) then
      icolor(itemp) = jtemp
    else
      write ( *, * ) 'Your object number was out of bounds.'
    end if
!
!  CC: choose color contour labels
!
!  For some strange reason, in order to make the color table
!  active, we have to call NEWFRM!
!
  else if ( s_eqi ( command(1:2),'cc') ) then

    if ( dev == ' ' ) then
      write ( *, * ) ' '
      write ( *, * ) 'CRYSTAL_PLOT - Warning!'
      write ( *, * ) '  Please use the "DEV" command to choose'
      write ( *, * ) '  a device, and THEN the "CC" command.'
      go to 10
    end if

    if ( s_eqi ( command(1:3),'cc = ') ) then
      read(command(4:),*,err = 50,end=50)itable
    else
      write ( *, * ) ' '
      write ( *, * ) 'Built in color tables include:'
      write ( *, * ) ' '
      write ( *, * ) '1  low black to high white.'
      write ( *, * ) '2  low blue to high yellow.'
      write ( *, * ) '3  low red, high blue, with bands between.'
      write ( *, * ) '4  low red, yellow, green, blue, high white.'
      write ( *, * ) '5  low white, blue, green, yellow, high red'
      write ( *, * ) '6  low blue to high red.'
      write ( *, * ) '7  linear scale between 2 user colors.'
      write ( *, * ) '8  linear scale between N user colors.'
      write ( *, * ) '9  low white to high black.'
      write ( *, * ) ' '
      write ( *, * ) 'Enter a color table index between 1 and 9,'
      write ( *, * ) 'or 0 to enter a color table from a file.'

      read(*,*,end = 50,err=50)itable
      write(17,*)itable
      if ( echo ) then
        write ( *, * ) itable
      end if
    end if

    call gettab(dev,echo,filgrf,grace,icmax,icmin,ierror, &
      iplot,itable,lgopen,ovrlay)

    if ( itable == 1 .or. itable == 9 ) then
      jcmax = 200
      jcmin = 32
    else
      jcmax = 255
      jcmin = 2
    end if
 
    write ( *, * ) ' '
    write ( *, * ) 'Lowest color used will be JCMIN =  ',jcmin
    write ( *, * ) 'Highest color used will be JCMAX = ',jcmax
!
!  CH: center half.
!
  else if ( s_eqi ( command,'ch') ) then

    temp = 0.25*(xsmax-xsmin)
    xsmin = xsmin+temp
    xsmax = xsmax-temp

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  CP: show corrected pressure contours.
!
  else if ( s_eqi ( command,'cp') ) then

    show(5) = .not.show(5)
    if ( show(5) ) then
      write ( *, * ) 'Corrected pressures will be shown.'
    else
      write ( *, * ) 'Corrected pressures will NOT be shown.'
    end if
!
!  CPC: show corrected pressure color plots.
!
  else if ( s_eqi ( command,'cpc') ) then

    show(12) = .not.show(12)
    if ( show(12) ) then
      write ( *, * ) 'Corrected pressure colors will be shown.'
    else
      write ( *, * ) 'Corrected pressure colors will NOT be shown.'
    end if
!
!  CRUcible: show crucible wall.
!
  else if ( s_eqi ( command(1:3),'cru') ) then

    show(19) = .not.show(19)
    if ( show(19) ) then
      write ( *, * ) 'The crucible wall will be shown.'
    else
      write ( *, * ) 'The crucible wall will NOT be shown.'
    end if
!
!  CRYB: show the crystal boundary.
!
  else if ( s_eqi ( command(1:4),'cryb') ) then

    show(20) = .not.show(20)
    if ( show(20) ) then
      write ( *, * ) 'The crystal boundary will be shown.'
    else
      write ( *, * ) 'The crystal boundary will NOT be shown.'
    end if
!
!  CRYS: display items in the crystal.
!
  else if ( s_eqi ( command,'crys') .or. s_eqi ( command,'crystal') ) then

    show(39) = .not.show(39)
    if ( show(39) ) then
      write ( *, * ) 'Items within the crystal will be shown.'
    else
      write ( *, * ) 'Items within the crystal will NOT be shown.'
    end if
!
!  CRYSC: display shaded crystal.
!
  else if ( s_eqi ( command,'crysc') ) then

    show(44) = .not.show(44)
    if ( show(44) ) then
      write ( *, * ) 'The crystal interior will be colored.'
    else
      write ( *, * ) 'The crystal interior will NOT be colored.'
    end if
!
!  CV: show control volumes.
!
  else if ( s_eqi ( command,'cv') ) then

    show(2) = .not.show(2)
    if ( show(2) ) then
      write ( *, * ) 'Control volumes will be shown.'
    else
      write ( *, * ) 'Control volumes will NOT be shown.'
    end if
!
!  CVN: show control volume numbers.
!
  else if ( s_eqi ( command,'cvn') ) then

    show(43) = .not.show(43)
    if ( show(43) ) then
      write ( *, * ) 'Control volume numbers will be shown.'
    else
      write ( *, * ) 'Control volume numbers will NOT be shown.'
    end if
!
!  DAT = 
!
  else if ( s_eqi ( command(1:3),'dat') ) then

20      continue

    if ( s_eqi ( command(1:4),'dat = ') ) then

      fildat = command(5:)

    else
   
      write ( *, * ) ' '
      write ( *, * ) 'Enter the name of the new input data file:'
      read(*,'(a)',err = 50,end=50)fildat
      write(17,'(a)') trim ( fildat )
      if ( echo ) then
        write(*,'(a)') trim ( fildat )
      end if

    end if
 
    open(unit = 10,file=fildat,status='old',err=60)

    call rsread(b1jbl,b2jbl,b3jbl,cost,e,gamt,icrys1,icrys2,jcrys1,jcrys2, &
      kcrys,kmelt,kvoid,l,m,maxbot,maxl,maxm,nbot,p,pc,psi,rueta,ruksi,t, &
      te,tk,tnow,u,v,vmag,vort,w,xbot,xc,xp,ybot,yc,yp)

    close(unit = 10)
 
    reflec = .false.

    write ( *, * ) ' '
    write ( *, * ) 'CRYSTAL_PLOT - Note:'
    write ( *, * ) ' '
    write ( *, * ) '  Region is ',l,' rows wide and ',m,' rows high.'
    write ( *, * ) '  Time is ',tnow
    write ( *, * ) '  Cost function is ',cost
    write ( *, * ) '  ICRYS1 = ',icrys1
    write ( *, * ) '  ICRYS2 = ',icrys2
    write ( *, * ) '  JCRYS1 = ',jcrys1
    write ( *, * ) '  JCRYS2 = ',jcrys2
 
    call rsize(delx,dely,grace,l,m,maxl,maxm,ncflag,npflag,srange,vpflag, &
      xc,x1max,x1min,x2max,x2min,xmax,xmin,xsmax,xsmin,yc,y1max, &
      y1min,y2max,y2min,ymax,ymin,ysmax,ysmin)

    xtmax = xsmax
    xtmin = xsmin
    ytmax = ysmax
    ytmin = ysmin
!
!  'DEV = ' Choose the graphics device.
!
  else if ( s_eqi ( command(1:3),'dev') ) then

    if ( dev /= ' ' ) then
      write ( *, * ) ' '
      write ( *, * ) 'CRYSTAL_PLOT - Error!'
      write ( *, * ) '  You have already chosen device '//dev
      write ( *, * ) '  You may not change your mind!'
      go to 10
    end if

    if ( s_eqi ( command(1:4),'dev = ') ) then
      dev = trim ( command(5:) )
    else
      write ( *, * ) ' '
      write ( *, * ) 'Enter the graphics device desired.'
      write ( *, * ) ' '
      write ( *, * ) 'Options include:'
      write ( *, * ) ' '
      write ( *, * ) 'CGMB  output to a CGM binary file.'
      write ( *, * ) 'PS    output to a PostScript file.'
      write ( *, * ) 'XWS   output to an X window screen.'
      read(*,'(a)',end = 50,err=50)dev
      write(17,'(a)') trim ( dev )
      if ( echo ) then
        write(*,'(a)') trim ( dev )
      end if
    end if

    if ( s_eqi ( dev(1:3),'cgm') ) then
      dev = 'cgmb'
      write ( *, * ) 'Output will be to a CGM binary file "cplot.cgm".'
    else if ( s_eqi ( dev,'ps') ) then
      write ( *, * ) 'Output will be to a PostScript file "cplot.ps".'
    else if ( s_eqi ( dev,'xws') ) then
      write ( *, * ) 'Output will be to an X window screen.'
    else
      write ( *, * ) 'Your device '//dev//' was not recognized!'
      dev = ' '
    end if
!
!  DIF
!
  else if ( s_eqi ( command(1:3),'dif') ) then
   
30      continue

    write ( *, * ) ' '
    write ( *, * ) 'Enter the name of the first input data file:'
    read(*,'(a)',err = 50,end=50)fildat
    write(17,'(a)') trim ( fildat )
    if ( echo ) then
      write(*,'(a)') trim ( fildat )
    end if
 
    open(unit = 10,file=fildat,status='old',err=70)

    call rsread(b1jbl,b2jbl,b3jbl,cost,e,gamt,icrys1,icrys2,jcrys1,jcrys2, &
      kcrys,kmelt,kvoid,l,m,maxbot,maxl,maxm,nbot,p,pc,psi,rueta,ruksi,t, &
      te,tk,tnow,u,v,vmag,vort,w,xbot,xc,xp,ybot,yc,yp)

    close(unit = 10)

    reflec = .false.
 
40      continue

    write ( *, * ) ' '
    write ( *, * ) 'Enter the name of the second input data file:'
    read(*,'(a)',err = 50,end=50)fildat
    write(17,'(a)') trim ( fildat )
    if ( echo ) then
      write(*,'(a)') trim ( fildat )
    end if
 
    open(unit = 10,file=fildat,status='old',err=80)

    call rsdiff(cost2,dxcdp,dxdp,dycdp,dydp,e,gamt,l,m,maxl,maxm,nbot,p,pc, &
      psi,rueta,ruksi,t,te,tempa,tk,tnow2,u,v,vmag,vort,w,xc,xp,yc,yp)

    close(unit = 10)
    write ( *, * ) 'Region is ',l,' rows wide and ',m,' rows high.'
    write ( *, * ) 'Time 1 is ',tnow
    write ( *, * ) 'Time 2 is ',tnow2
    write ( *, * ) 'Cost function 1 is ',cost
    write ( *, * ) 'Cost function 2 is ',cost2
    write ( *, * ) 'ICRYS1 = ',icrys1
    write ( *, * ) 'JCRYS2 = ',jcrys2
 
    call rsize(delx,dely,grace,l,m,maxl,maxm,ncflag,npflag,srange,vpflag, &
      xc,x1max,x1min,x2max,x2min,xmax,xmin,xsmax,xsmin,yc,y1max, &
      y1min,y2max,y2min,ymax,ymin,ysmax,ysmin)
!
!  DOUBLE
!
  else if ( s_eqi ( command(1:3),'dou') ) then

    reflec = .true.

    write ( *, * ) ' '
    write ( *, * ) 'The region will be reflected around the Y axis.'

    call double(e,gamt,jcrys1,jcrys2,kcrys,kmelt,kvoid,l,m, &
      maxbot,maxl,maxm,nbot,ncflag,npflag,p,pc,psi,rueta,ruksi, &
      t,te,tk,u,v,vmag,vort,vpflag,w,xbot,xc,xp,ybot,yc,yp)

    write ( *, * ) ' '
    write ( *, * ) 'CRYSTAL_PLOT - Note:'
    write ( *, * ) '  The doubled region is ',l,' rows wide and ', m,' rows high.'
!
!  Force recomputation of the plot window by setting XSMIN = XSMAX=0.
!
    xsmin = 0.0
    xsmax = 0.0

    call rsize(delx,dely,grace,l,m,maxl,maxm,ncflag,npflag, &
      srange,vpflag,xc,x1max,x1min,x2max,x2min,xmax,xmin,xsmax,xsmin,yc,y1max, &
      y1min,y2max,y2min,ymax,ymin,ysmax,ysmin)
!
!  DXCDP: show (dXCdP,dYCdP) vector plots.
!
  else if ( s_eqi ( command(1:3),'dxc') ) then

    show(27) = .not.show(27)
    if ( show(27) ) then
      write ( *, * ) '(dXCdP,dYCdP) vector plots will be shown.'
    else
      write ( *, * ) '(dXCdP,dYCdP) vector plots will NOT be shown.'
    end if
!
!  ECHO  User input will be echoed to the output file.
!
  else if ( s_eqi ( command,'echo') ) then
    echo = .not.echo
    if ( echo ) then
      write(17,'(a)')command
      if ( echo ) then
        write(*,'(a)')command
      end if
      write ( *, * ) 'User input will be echoed to the output file.'
    else
      write ( *, * ) 'User input will NOT be echoed to the output file.'
    end if
!
!  FILE = : set the name of the graphics output file.
!
  else if ( s_eqi ( command(1:4),'file') ) then

    if ( iplot > 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'Display - Warning:'
      write ( *, * ) '  It is too late to specify a plot file name.'
    else if ( s_eqi ( command(1:5),'file = ') ) then

      filgrf = trim ( command(6:) ) 
    else
      write ( *, * ) ' '
      write ( *, * ) 'Enter the plot file name.'
      read(*,'(a)',err = 50,end=50)filgrf 
      write(17,'(a)') trim ( filgrf )
      if ( echo ) then
        write(*,'(a)') trim ( filgrf )
      end if
    end if
!
!  FRAME: show the frame
!
  else if ( s_eqi ( command,'frame') ) then

    show(3) = .true.
    write ( *, * ) 'A frame will be shown around the picture.'
!
!  FS: show free surface.
!
  else if ( s_eqi ( command(1:2),'fs') ) then

    show(21) = .not.show(21)
    if ( show(21) ) then
      write ( *, * ) 'The free surface will be shown.'
    else
      write ( *, * ) 'The free surface will NOT be shown.'
    end if
!
!  FULL: show full picture.
!
  else if ( s_eqi ( command,'full') ) then
!
!  Force recomputation of the plot window by setting XSMIN = XSMAX=0.
!
    xsmin = 0.0
    xsmax = 0.0

    call rsize(delx,dely,grace,l,m,maxl,maxm,ncflag,npflag,srange,vpflag, &
      xc,x1max,x1min,x2max,x2min,xmax,xmin,xsmax,xsmin,yc,y1max, &
      y1min,y2max,y2min,ymax,ymin,ysmax,ysmin)

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  GRACE = : set the grace margin.
!
  else if ( s_eqi ( command(1:5),'grace') ) then

    if ( s_eqi ( command(1:6),'grace = ') ) then
      read(command(7:),*,err = 90,end=90)grace
    else
      write ( *, * ) 'Enter the grace margin:'
      read(*,*)grace
      write(17,*)grace
      if ( echo ) then
        write ( *, * ) grace
      end if
    end if

    write ( *, * ) ' '
    write ( *, * ) 'CRYSTAL_PLOT - Note:'
    write ( *, * ) '  The grace margin was set to GRACE  =  ',grace

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  GRAPH
!
  else if ( s_eqi ( command(1:1),'g') ) then

    if ( l < 1 .or. m.lt.1 ) then
      write ( *, * ) ' '
      write ( *, * ) 'Whoa!  Please read data with the DAT command!'
    else if ( dev == ' ' ) then
      write ( *, * ) ' '
      write ( *, * ) 'Please use the DEV command first!'
    else

      call graph(delx,dely,dev,dxcdp,dycdp,echo,filgrf,gamt,icmax, &
        icmin,icolor,icrys1,icrys2,iplot,itable,jcmax,jcmin,jcrys1, &
        jcrys2,kcrys,kmelt,kvoid,l,lbar,lflag,lgopen,m,maxbot,maxl, &
        maxm,maxobj,nbot,ncflag, &
        ncon,nflag,npflag,nxskip,nyskip,object,ovrlay,p,pc,psi,s1,s2, &
        scalecv,scalenb,scalenc,scalenp,scalev,show, &
        srange,t,te,title,title2,tk,u,v,vpflag,vmag,vort,w, &
        x1max,x1min,x2max,x2min,xbot,xc,xp,xsmax,xsmin,y1max,y1min, &
        y2max,y2min,ybot,yc,yp,ysmax,ysmin)

    end if
!
! HALF
!
  else if ( s_eqi ( command(1:4),'half') ) then

    reflec = .false.

    write ( *, * ) ' '
    write ( *, * ) 'The region will NOT be reflected around the Y axis.'

    call half(e,gamt,jcrys1,jcrys2,kcrys,kmelt,kvoid,l,m, &
      maxbot,maxl,maxm,nbot,ncflag,npflag,p,pc,psi,rueta,ruksi, &
      t,te,tk,u,v,vmag,vort,vpflag,w,xbot,xc,xp,ybot,yc,yp)

    write ( *, * ) ' '
    write ( *, * ) 'CRYSTAL_PLOT - Note:'
    write ( *, * ) '  The undoubled region is ',l,' rows wide and ', &
      m,' rows high.'
!
!  Force recomputation of the plot window by setting XSMIN = XSMAX=0.
!
    xsmin = 0.0
    xsmax = 0.0

    call rsize(delx,dely,grace,l,m,maxl,maxm,ncflag,nflag,srange,vpflag, &
      xc,x1max,x1min,x2max,x2min,xmax,xmin,xsmax,xsmin,yc,y1max, &
      y1min,y2max,y2min,ymax,ymin,ysmax,ysmin)
!
!  HELP
!
  else if ( s_eqi ( command(1:1),'h') ) then
    call help
!
!  ICMAX = : set the maximum available color index.
!
  else if ( s_eqi ( command(1:6),'icmax = ') ) then

    read(command(7:),*,err = 50,end=50)icmax
    if ( icmax > 255 ) then
      write ( *, * ) 'ICMAX must be no more than 255'
      icmax = 255
    end if
    write ( *, * ) 'Maximum available color set to ',icmax
!
!  ICMIN = : set the minimum available color index.
!
  else if ( s_eqi ( command(1:6),'icmin = ') ) then

    read(command(7:),*,err = 50,end=50)icmin
    if ( icmin < 2 ) then
      write ( *, * ) 'ICMIN must be no less than 2'
      icmin = 2
    end if
    write ( *, * ) 'Minimum color set to ',icmin
!
!  JCMAX = : set the maximum used color index.
!
  else if ( s_eqi ( command(1:6),'jcmax = ') ) then

    read(command(7:),*,err = 50,end=50)jcmax
    if ( jcmax > 255 ) then
      write ( *, * ) 'JCMAX must be no more than 255.'
      jcmax = 255
    end if
    write ( *, * ) 'Maximum used color set to ',jcmax
!
!  JCMIN = : set the minimum used color index.
!
  else if ( s_eqi ( command(1:6),'jcmin = ') ) then

    read(command(7:),*,err = 50,end=50)jcmin
    if ( jcmin < 2 ) then
      write ( *, * ) 'JCMIN must be no less than 2.'
      jcmin = 2
    end if
    write ( *, * ) 'Minimum used color set to ',jcmin
!
!  LH: left half.
!
  else if ( s_eqi ( command,'lh') ) then

    xsmax = xsmin+0.5*(xsmax-xsmin)
    temp = 0.25*(ysmax-ysmin)

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  MC: middle center quarter.
!
  else if ( s_eqi ( command,'mc') ) then

    temp = 0.25*(xsmax-xsmin)
    xsmin = xsmin+temp
    xsmax = xsmax-temp
    temp = 0.25*(ysmax-ysmin)
    ysmax = ysmax-temp
    ysmin = ysmin+temp

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  MELT: display items in the melt.
!
  else if ( s_eqi ( command,'melt') ) then

    show(40) = .not.show(40)
    if ( show(40) ) then
      write ( *, * ) 'Items within the melt will be shown.'
    else
      write ( *, * ) 'Items within the melt will NOT be shown.'
    end if
!
!  MELTC: display shaded melt.
!
  else if ( s_eqi ( command,'meltc') ) then

    show(45) = .not.show(45)
    if ( show(45) ) then
      write ( *, * ) 'The interior of the melt will be colored.'
    else
      write ( *, * ) 'The interior of the melt will NOT be colored.'
    end if
!
!  MH: middle half.
!
  else if ( s_eqi ( command,'mh') ) then

    temp = 0.25*(ysmax-ysmin)
    ysmax = ysmax-temp
    ysmin = ysmin+temp

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  ML: middle left quarter.
!
  else if ( s_eqi ( command,'ml') ) then

    xsmax = xsmin+0.5*(xsmax-xsmin)
    temp = 0.25*(ysmax-ysmin)
    ysmax = ysmax-temp
    ysmin = ysmin+temp

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)

!
!  MR: middle right quarter.
!
  else if ( s_eqi ( command,'mr') ) then

    xsmin = xsmin+0.5*(xsmax-xsmin)
    temp = 0.25*(ysmax-ysmin)
    ysmax = ysmax-temp
    ysmin = ysmin+temp

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  NB: show nodes on bottom that define crucible shape.
!
  else if ( s_eqi ( command,'nb') ) then

    show(34) = .not.show(34)
    if ( show(34) ) then
      write ( *, * ) 'Nodes on bottom will be shown.'
    else
      write ( *, * ) 'Nodes on bottom will NOT be shown.'
    end if
!
!  NC: show corner nodes
!
  else if ( s_eqi ( command,'nc') ) then

    show(29) = .not.show(29)
    if ( show(29) ) then
      write ( *, * ) 'Corner nodes will be shown.'
    else
      write ( *, * ) 'Corner nodes will NOT be shown.'
    end if
!
!  NCON = : Set the number of contour lines.
!
  else if ( s_eqi ( command(1:4),'ncon') ) then

    if ( s_eqi ( command(1:5),'ncon = ') ) then
      read(command(6:),*,err = 50,end=50)ncon
      write ( *, * ) 'Number of contour lines set to ',ncon
    else
      write ( *, * ) 'Enter number of contour lines.'
      read(*,*)ncon
      write(17,*)ncon
      if ( echo ) then
        write ( *, * ) ncon
      end if
    end if
!
!  NOFRAME: no frame, please
!
  else if ( s_eqi ( command,'noframe') ) then

    show(3) = .false.
    write ( *, * ) 'No frame will be shown around the picture.'
!
!  NP: show primary nodes
!
  else if ( s_eqi ( command,'np') ) then

    show(4) = .not.show(4)
    if ( show(4) ) then
      write ( *, * ) 'Primary nodes will be shown.'
    else
      write ( *, * ) 'Primary nodes will NOT be shown.'
    end if
!
!  NXSKIP = : skip value for column nodes.
!
  else if ( s_eqi ( command(1:7),'nxskip = ') ) then

    read(command(8:),*,err = 50,end=50)itemp
    if ( itemp <= 0 ) then
      write ( *, * ) 'Error!  NXSKIP must be 1 or more!'
    else
      nxskip = itemp
      write ( *, * ) 'Skip value for column nodes is ',nxskip
    end if
!
!  NYSKIP = : skip value for row nodes.
!
  else if ( s_eqi ( command(1:7),'nyskip = ') ) then

    read(command(8:),*,err = 50,end=50)itemp
    if ( itemp <= 0 ) then
      write ( *, * ) 'Error!  NYSKIP must be 1 or more!'
    else
      nyskip = itemp
      write ( *, * ) 'Skip value for row nodes is ',nyskip
    end if
!
!  OVERLAY: Switch the overlay value.
!
  else if ( s_eqi ( command,'overlay') ) then
    ovrlay = .not.ovrlay
    if ( ovrlay ) then
      write ( *, * ) 'Plots will be overlayed until next OVERLAY.'
    else
      write ( *, * ) 'This overlay plot is done.'
      call newfrm
    end if
!
!  P: show pressure contours.
!
  else if ( s_eqi ( command,'p') ) then

    show(5) = .not.show(5)
    if ( show(5) ) then
      write ( *, * ) 'Pressures will be shown.'
    else
      write ( *, * ) 'Pressures will NOT be shown.'
    end if
!
!  PC: show pressure color plots.
!
  else if ( s_eqi ( command,'pc') ) then

    show(12) = .not.show(12)
    if ( show(12) ) then
      write ( *, * ) 'Pressure colors will be shown.'
    else
      write ( *, * ) 'Pressure colors will NOT be shown.'
    end if
!
!  Q  =  Quit.
!  QY  =  QUIT NOW!
!
  else if ( s_eqi ( command(1:1),'q') ) then

    if ( s_eqi ( command(2:2),'y') ) then
      command = 'y'
    else
      write ( *, * ) 'Enter "y" to confirm you want to quit!'
      read(*,'(a)')command
      call flushl(command)
      write(17,'(a)') trim ( command )
      if ( echo ) then
        write(*,'(a)') trim ( command )
      end if
    end if

    if ( s_eqi ( command,'y') ) then

      if ( lgopen ) then
        call grfcls

        if ( lnei ( dev, 'cgmb' ) ) then
          call delete ( 'cgmout' )
        end if
      end if

      close(unit = 17)

      stop

    end if
!
!  RH: right half.
!
  else if ( s_eqi ( command,'rh') ) then

    temp = 0.50*(xsmax-xsmin)
    xsmin = xsmin+temp

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  S: show stream lines
!
  else if ( s_eqi ( command,'s') ) then

    show(6) = .not.show(6)
    if ( show(6) ) then
      write ( *, * ) 'Stream lines will be shown.'
    else
      write ( *, * ) 'Stream lines will NOT be shown.'
    end if
!
!  SC: show stream colors
!
  else if ( s_eqi ( command,'sc') ) then

    show(35) = .not.show(35)
    if ( show(35) ) then
      write ( *, * ) 'Stream colors will be shown.'
    else
      write ( *, * ) 'Stream colors will NOT be shown.'
    end if
!
!  SCALENB = : set the bottom node scale.
!
  else if ( s_eqi ( command(1:8),'scalenb = ') ) then

    read(command(9:),*,err = 50,end=50)scalenb
    write ( *, * ) 'Scale factor for bottom nodes set to ',scalenb
!
!  SCALENC = : set the corner node scale.
!
  else if ( s_eqi ( command(1:8),'scalenc = ') ) then

    read(command(9:),*,err = 50,end=50)scalenc
    write ( *, * ) 'Scale factor for corner nodes set to ',scalenc
!
!  SCALENP = : set the primary node scale.
!
  else if ( s_eqi ( command(1:8),'scalenp = ') ) then

    read(command(9:),*,err = 50,end=50)scalenp
    write ( *, * ) 'Scale factor for primary nodes set to ',scalenp
!
!  SCALEV = : set the velocity vector scale.
!
  else if ( s_eqi ( command(1:7),'scalev = ') ) then

    read(command(8:),*,err = 50,end=50)scalev
    write ( *, * ) 'Velocity vector scale factor set to ',scalev
!
!  TC: top center quarter.
!
  else if ( s_eqi ( command,'tc') ) then

    temp = 0.25*(xsmax-xsmin)
    xsmin = xsmin+temp
    xsmax = xsmax-temp
    ysmin = ysmin+0.5*(ysmax-ysmin)

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  TE: show turbulent epsilon contours.
!
  else if ( s_eqi ( command,'te') ) then

    show(22) = .not.show(22)
    if ( show(22) ) then
      write ( *, * ) 'Turbulent epsilon will be shown.'
    else
      write ( *, * ) 'Turbulent epsilon will NOT be shown.'
    end if
!
!  TEMP: show temperature contours.
!
  else if ( s_eqi ( command,'temp') ) then

    show(11) = .not.show(11)
    if ( show(11) ) then
      write ( *, * ) 'Temperatures will be shown.'
    else
      write ( *, * ) 'Temperatures will NOT be shown.'
    end if
!
!  TEMPC: show temperature colors.
!
  else if ( s_eqi ( command,'tempc') ) then

    show(13) = .not.show(13)
    if ( show(13) ) then
      write ( *, * ) 'Temperature colors will be shown.'
    else
      write ( *, * ) 'Temperature colors will NOT be shown.'
    end if
!
!  TH: top half
!
  else if ( s_eqi ( command,'th') ) then

    ysmin = ysmin+0.5*(ysmax-ysmin)

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  TITLE2 = : enter subtitle
!
  else if ( s_eqi ( command(1:7),'title2 = ') ) then
    title2 = trim ( command(8:) )
  else if ( s_eqi ( command(1:6),'title2') ) then
    write ( *, * ) 'Enter the subtitle:'
    read(*,'(a)',end = 50)title2
    write(17,'(a)') trim ( title2 )
    if ( echo ) then
      write(*,'(a)') trim ( title2 )
    end if
!
!  TITLE = : enter title
!
  else if ( s_eqi ( command(1:6),'title = ') ) then
    title = trim ( command(7:) )
  else if ( s_eqi ( command(1:5),'title') ) then
    write ( *, * ) 'Enter the plot title:'
    read(*,'(a)',end = 50)title
    write(17,'(a)') trim ( title )
    if ( echo ) then
      write(*,'(a)') trim ( title )
    end if
!
!  TK: show turbulent kinetic energy contours.
!
  else if ( s_eqi ( command,'tk') ) then

    show(23) = .not.show(23)
    if ( show(23) ) then
      write ( *, * ) 'Turbulent KE will be shown.'
    else
      write ( *, * ) 'Turbulent KE will NOT be shown.'
    end if
!
!  TL: top left quarter.
!
  else if ( s_eqi ( command,'tl') ) then

    xsmax = xsmin+0.5*(xsmax-xsmin)
    ysmin = ysmin+0.5*(ysmax-ysmin)

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  TR: top right quarter.
!
  else if ( s_eqi ( command,'tr') ) then

    xsmin = xsmin+0.5*(xsmax-xsmin)
    ysmin = ysmin+0.5*(ysmax-ysmin)

    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
      y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  U: show unit velocity direction vectors.
!
  else if ( s_eqi ( command,'u') ) then

    show(9) = .not.show(9)
    if ( show(9) ) then
      write ( *, * ) 'Unit velocity direction field will be shown.'
    else
      write ( *, * ) 'Unit velocity direction field will NOT be shown.'
    end if
!
!  V: show velocity vectors.
!
  else if ( s_eqi ( command,'v') ) then

    show(8) = .not.show(8)
    if ( show(8) ) then
      write ( *, * ) 'Velocities will be shown.'
    else
      write ( *, * ) 'Velocities will NOT be shown.'
    end if
!
!  VIS: show viscosity contours.
!
  else if ( s_eqi ( command,'vis') ) then

    show(24) = .not.show(24)
    if ( show(24) ) then
      write ( *, * ) 'Viscosity contours will be in next plot.'
    else
      write ( *, * ) 'Viscosity contours will NOT be in next plot.'
    end if
!
!  VM: show velocity magnitude contours.
!
  else if ( s_eqi ( command,'vm') ) then

    show(36) = .not.show(36)
    if ( show(36) ) then
      write ( *, * ) 'Velocity magnitudes will be shown.'
    else
      write ( *, * ) 'Velocity magnitudes will NOT be shown.'
    end if
!
!  VMC: show velocity magnitude colors.
!
  else if ( s_eqi ( command,'vmc') ) then

    show(37) = .not.show(37)
    if ( show(37) ) then
      write ( *, * ) 'Velocity magnitude colors will be shown.'
    else
      write ( *, * ) 'Velocity magnitude colors will NOT be shown.'
    end if
!
!  VOID: display items in the void.
!
  else if ( s_eqi ( command,'void') ) then

    show(38) = .not.show(38)
    if ( show(38) ) then
      write ( *, * ) 'Items in the void will be shown.'
    else
      write ( *, * ) 'Items in the void will NOT be shown.'
    end if
!
!  VOIDC: display shaded void.
!
  else if ( s_eqi ( command,'voidc') ) then

    show(46) = .not.show(46)
    if ( show(46) ) then
      write ( *, * ) 'The void will be colored.'
    else
      write ( *, * ) 'The void will NOT be colored.'
    end if
!
!  VORT: show vorticity contours.
!
  else if ( s_eqi ( command,'vort') ) then

    show(41) = .not.show(41)
    if ( show(41) ) then
      write ( *, * ) 'Vorticity contours will be shown.'
    else
      write ( *, * ) 'Vorticity contours will NOT be shown.'
    end if
!
!  VORTC: show vorticity colors.
!
  else if ( s_eqi ( command,'vortc') ) then

    show(42) = .not.show(42)
    if ( show(42) ) then
      write ( *, * ) 'Vorticity colors will be shown.'
    else
      write ( *, * ) 'Vorticity colors will NOT be shown.'
    end if
!
!  VPN: set visible primary nodes.
!
  else if ( s_eqi ( command,'vpn') ) then
 
    call vizpn(echo,l,m,maxl,maxm,vpflag,xp,yp)
!
!  W: show angular momentum contours.
!
  else if ( s_eqi ( command,'w') ) then

    show(10) = .not.show(10)
    if ( show(10) ) then
      write ( *, * ) 'Angular momentum contours will be in next plot.'
    else
      write ( *, * ) 'Angular momentum contours will NOT be shown.'
    end if
!
!  WC: show angular momentum color plots.
!
  else if ( s_eqi ( command,'wc') ) then

    show(14) = .not.show(14)
    if ( show(14) ) then
      write ( *, * ) 'Angular momentum colors will be shown.'
    else
      write ( *, * ) 'Angular momentum colors will NOT be shown.'
    end if
!
!  X: set the data window.
!
  else if ( s_eqi ( command,'x') ) then

    call getwin(echo,grace,srange,xmax,xmin,x1max,x1min,x2max, &
      x2min,xsmax,xsmin,ymax,ymin,y1max,y1min,y2max,y2min,ysmax,ysmin)
!
!  XC: show X coordinate contours.
!
  else if ( s_eqi ( command,'xc') ) then

    show(30) = .not.show(30)
    if ( show(30) ) then
      write ( *, * ) 'X coordinate contours will be shown.'
    else
      write ( *, * ) 'X coordinate contours will NOT be shown.'
    end if
!
!  XCC: show X coordinate colors.
!
  else if ( s_eqi ( command,'xcc') ) then

    show(31) = .not.show(31)
    if ( show(31) ) then
      write ( *, * ) 'X coordinate colors will be shown.'
    else
      write ( *, * ) 'X coordinate colors will NOT be shown.'
    end if
!
!  YC: show Y coordinate contours.
!
  else if ( s_eqi ( command,'yc') ) then

    show(32) = .not.show(32)
    if ( show(32) ) then
      write ( *, * ) 'Y coordinate contours will be shown.'
    else
      write ( *, * ) 'Y coordinate contours will NOT be shown.'
    end if
!
!  YCC: show Y coordinate colors.
!
  else if ( s_eqi ( command,'ycc') ) then

    show(33) = .not.show(33)
    if ( show(33) ) then
      write ( *, * ) 'Y coordinate colors will be shown.'
    else
      write ( *, * ) 'Y coordinate colors will NOT be shown.'
    end if
!
!  Unrecognized command.
!
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CRYSTAL_PLOT - Warning!'
    write ( *, '(a)' ) '  Unrecognized command:' // trim ( command )
  end if

  go to 10

50    continue
  write ( *, * ) ' '
  write ( *, * ) 'CRYSTAL_PLOT - Fatal error!'
  write ( *, * ) '  Input error or end of file.'
  write ( *, * ) ' '
  write ( *, * ) '  Warning!  The graphics file may not be'
  write ( *, * ) '  properly terminated.'
  stop

60    continue
  write ( *, * ) ' '
  write ( *, * ) 'CRYSTAL_PLOT - Warning!'
  write ( *, * ) '  The file you named could not be opened!'
  fildat = ' '
  command = 'dat'
  go to 11

70    continue
  write ( *, * ) ' '
  write ( *, * ) 'CRYSTAL_PLOT - Warning!'
  write ( *, * ) '  The file you named could not be opened!'
  fildat = ' '
  command = 'dif'
  go to 11

80    continue
  write ( *, * ) ' '
  write ( *, * ) 'CRYSTAL_PLOT - Warning!'
  write ( *, * ) '  The file you named could not be opened!'
  fildat = ' '
  command = 'dif'
  go to 11

90    continue
  write ( *, * ) ' '
  write ( *, * ) 'CRYSTAL_PLOT - Warning!'
  write ( *, * ) '  Your input was not understood.'
  go to 10

end
subroutine arrow ( xstart, ystart, xtip, ytip )

!*****************************************************************************80
!
!! ARROW can be used to draw an arrow at any point on a graph.
!
!  Discussion:
!
!    The arrow will stretch between two user specified points.
!
!    The "head" of the arrow may be fatter or thinner than expected
!    if the X and Y scales of the graph are not in the same
!    proportions.
!
!                     left
!                     |\
!                     | \
!    start ------- base  tip
!                     | /
!                     |/
!                     rite
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XSTART, YSTART, the starting point for the
!    arrow.
!
!    Input, real XTIP, YTIP, the end point for the arrow.
!
  implicit none

  real alpha
  real del
  real dist
  real, parameter :: pi = 3.1415926E+00
  real theta
  real xbase
  real xleft
  real xrite
  real xstart
  real xtip
  real ybase
  real yleft
  real yrite
  real ystart
  real ytip

  if ( xstart == xtip .and. ystart == ytip)return

  theta = 0.5E+00*pi-atan2(2.0E+00,1.0E+00)
  dist = sqrt((xtip-xstart)**2+(ytip-ystart)**2)
  alpha = atan2(ytip-ystart,xtip-xstart)
  del = sqrt(5.0E+00)*dist/3.0E+00

  call movcgm(xstart,ystart)

  xbase = xstart+dist*cos(alpha)*2.0E+00/3.0E+00
  ybase = ystart+dist*sin(alpha)*2.0E+00/3.0E+00
  call drwcgm(xbase,ybase)

  xleft = xstart+del*cos(alpha-theta)
  yleft = ystart+del*sin(alpha-theta)
  call drwcgm(xleft,yleft)

  call drwcgm(xtip,ytip)

  xrite = xstart+del*cos(alpha+theta)
  yrite = ystart+del*sin(alpha+theta)
  call drwcgm(xrite,yrite)

  call drwcgm(xbase,ybase)

  return
end
subroutine box(xmin,xmax,ymin,ymax)

!*****************************************************************************80
!
!! BOX draws a rectangle whose corners are specified by the user.
!
!  Discussion:
!
!    The rectangle drawn by box has the corners:
!
!      (XMIN,YMAX)   (XMAX,YMAX)
!
!      (XMIN,YMIN)   (XMAX,YMIN)
!
!    BOX can be used to place a rectangle anywhere in the picture.  
!    However, BOX may also be used to place a rectangle around the 
!    entire picture, producing a "frame".
!
!    The DRAWCGM routine PLYLIN is used to draw the box, and hence 
!    the color of the line may be set by calling the DRAWCGM routine 
!    LINCLR first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XMIN, XMAX, the minimum and maximum X 
!    coordinates of the box.
!
!    Input, real YMIN, YMAX, the minimum and maximum Y 
!    coordinates of the box.
!
  implicit none
!
  integer ( kind = 4 ), parameter :: npoints = 5

  real x(npoints)
  real xmax
  real xmin
  real y(npoints)
  real ymax
  real ymin

  x(1) = xmin
  y(1) = ymin

  x(2) = xmax
  y(2) = ymin

  x(3) = xmax
  y(3) = ymax

  x(4) = xmin
  y(4) = ymax

  x(5) = xmin
  y(5) = ymin

  call plylin(npoints,x,y)

  return
end
subroutine buzz(dev,x1max,x1min,y1max,y1min)

!*****************************************************************************80
!
!! BUZZ forces the graphical system to slow down a bit.
!
!  Discussion:
!
!    BUZZ is just "busy work" which seems to fix a problem that occurs 
!    when the XWS interface is used.  In that case, the last bit of the 
!    graph is not drawn.  So here, we just make the last bit of the graph 
!    something we don't care about.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*10 DEV, the graphics output device to be used. 
!    Current legal include:
! 
!      cgmb - CGM binary file.
!      ps   - PostScript file.
!      xws  - X window screen (interactive).
!
!    Output, real X1MAX, X1MIN, the maximum and minimum X 
!    coordinates of the plot, which includes a small grace margin.
!
!    Output, real Y1MAX, Y1MIN, the maximum and minimum Y 
!    coordinates of the plot, which includes a small grace margin.
!
  implicit none

  character ( len = 10 ) dev
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icolor
  real x1max
  real x1min
  real y1max
  real y1min

  if ( dev == 'xws' ) then

    icolor = 0
    call linclr(icolor)

    do i = 1,100

      call box(x1min,x1max,y1min,y1max)

    end do

  end if

  return
end
subroutine cbar(icolor,jcmax,jcmin,maxobj,ncon,smax,smin,srange,x1,x2,y1,y2)

!*****************************************************************************80
!
!! CBAR draws a color bar.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  ICOLOR Input, integer ( kind = 4 ) ICOLOR(MAXOBJ).
!         Contains the color indexes for each object.
!         However, in some cases, ICOLOR is actual a color table
!         index.
! 
!  JCMAX,
!  JCMIN  Input, integer ( kind = 4 ) JCMAX, JCMIN, the maximum and 
!         minimum color indices to use in the color bar.
!
!  MAXOBJ Input, integer ( kind = 4 ) MAXOBJ.
!         The number of graphical "objects".
! 
!  NCON   Input, integer ( kind = 4 ) NCON, the number of color contour
!         regions drawn, and hence, the number of colors
!         to be displayed in the color bar.
!
!  SMAX,
!  SMIN   Input, real SMAX, SMIN, the maximum and minimum
!         values of the quantity whose color contours are
!         being drawn.  These numbers will be printed along
!         with the color bar.
!
!  SRANGE Input, real SRANGE.
!         The maximum of XSMAX-XSMIN and YSMAX-YSMIN.
!         This gives the size of a square containing the data
!         window.
!
!  X1,
!  X2,
!  Y1,
!  Y2     Input, real X1, X2, Y1, Y2, specify the minimum and
!         maximum X and Y coordinates of the color bar.
!
  implicit none

  integer ( kind = 4 ) maxobj

  real angle
  character chrrel*14
  character chrtmp*14
  real cwide
  character ( len = 6 ) flush
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icolor(maxobj)
  integer ( kind = 4 ) jcmax
  integer ( kind = 4 ) jcmin
  integer ( kind = 4 ) jcolor
  integer ( kind = 4 ) lent
  integer ( kind = 4 ) ncon
  real pwide
  real smax
  real smin
  real srange
  real x
  real x1
  real x2
  real xcorn(5)
  real xl
  real xr
  real y
  real y1
  real y2
  real ycorn(5)

  ycorn(1) = y1
  ycorn(2) = y1
  ycorn(3) = y2
  ycorn(4) = y2
  ycorn(5) = y1

  call linclr(icolor(1))

  do i = 0,ncon

    xl = ((ncon+1-i)*x1+i*x2)/real(ncon+1)
    xr = ((ncon-i)*x1+(i+1)*x2)/real(ncon+1)

    xcorn(1) = xl
    xcorn(2) = xr
    xcorn(3) = xr
    xcorn(4) = xl
    xcorn(5) = xl

    jcolor = int(((ncon-i)*jcmin+i*jcmax)/real(ncon))
    call filclr(jcolor)

    call plygon(4,xcorn,ycorn)

    call plylin(5,xcorn,ycorn)

  end do
!
!  Print labels for the lowest and highest contours.
!
  cwide = 0.9*srange/40.0

  chrtmp = chrrel(smin)
  call s_blank_delete ( chrtmp )

  angle = 0.0
  x = x1
  y = y1-1.5*cwide
  pwide = srange
  flush = 'left'
  call s_plot(angle,cwide,pwide,chrtmp,x,y,flush)

  chrtmp = chrrel(smax)
  call s_blank_delete ( chrtmp )

  angle = 0.0
  x = x2
  y = y1-1.5*cwide
  pwide = srange
  flush = 'right'
  call s_plot(angle,cwide,pwide,chrtmp,x,y,flush)

  return
end
subroutine cbox(grace,x1max,x1min,y1max,y1min)

!*****************************************************************************80
!
!! CBOX draws a 16 by 16 color box.
!
!  Discussion:
!
!    You may want to call NEWFRM before or after calling CBOX, in
!    order to clear the frame!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  GRACE  Input, real GRACE.
!         The size of the "grace" margin on the plot.
!
!  X1MAX,
!  X1MIN  Output, real X1MAX, X1MIN, the maximum and minimum X 
!         coordinates of the plot, which includes a small grace margin.
!
!  Y1MAX,
!  Y1MIN  Output, real Y1MAX, Y1MIN, the maximum and minimum Y 
!         coordinates of the plot, which includes a small grace margin.
!
  implicit none

  integer ( kind = 4 ), parameter :: npoly = 5

  real grace
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icolor
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  real x(npoly)
  real x1max
  real x1min
  real x2max
  real x2min
  real y(npoly)
  real y1max
  real y1min
  real y2max
  real y2min
!
!  Set the coordinate system to be 0 < =  X <= 16.0, 
!  and 0.0 < =  Y <= 16.0
!
  x2min = 0.0E+00
  x2max = 16.0E+00
  y2min = 0.0E+00
  y2max = 16.0E+00

  x1min = x2min-grace*(x2max-x2min)
  x1max = x2max+grace*(x2max-x2min)
  y1min = y2min-grace*(y2max-y2min)
  y1max = y2max+grace*(y2max-y2min)

  call setwcd(x1min,y1min,x1max,y1max,ierror)
! 
!  Draw the color boxes.
!
  icolor = 0

  do i = 1,16

    y(1) = 16-i
    y(2) = 16-i
    y(3) = 17-i
    y(4) = 17-i

    do j = 1,16

      call filclr(icolor)
      icolor = icolor+1

      x(1) = j-1
      x(2) = j
      x(3) = j
      x(4) = j-1

      call plygon(4,x,y)

    end do
  end do
!
!  Draw black lines around the boxes.
!
  icolor = 1
  call linclr(icolor)

  do i = 1,16

    y(1) = 16-i
    y(2) = 17-i
    y(3) = 17-i
    y(4) = 16-i
    y(5) = y(1)

    do j = 1,16

      x(1) = j-1
      x(2) = j-1
      x(3) = j
      x(4) = j
      x(5) = x(1)

      call plylin(5,x,y)

    end do
  end do

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
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
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function chrrel(rval)

!*****************************************************************************80
!
!! CHRREL converts a real number to a right-justified string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  RVAL   Input, real RVAL, a real number.
!
!  CHRREL Output (through function value), character*14 CHRREL,
!         a right-justified character variable containing the
!         representation of RVAL, using a G14.6 format.
!
  implicit none

  character ( len = 14 ) chrrel
  character ( len = 14 ) chrtmp
  real rval
!
!  We can't seem to write directly into CHRREL because of compiler
!  quibbles.
!
  if ( real(int(rval)) == rval .and. abs(rval) < 1.0E+13 ) then

    write(chrtmp,'(i14)')int(rval)

  else

    write(chrtmp,'(g14.6)')rval

  end if

  chrrel = chrtmp
  return
end
subroutine colcon(a,aname,echo,icolor,jcmax,jcmin,l,lbar,m,maxl, &
  maxm,maxobj,ncon,nflag,srange,xp,x2max,x2min,yp,y1max,y2max)

!*****************************************************************************80
!
!! COLCON supervises the creation of a color contour plot.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  A      Input, real A(MAXL,MAXM), the quantity whose contour
!         is desired.
!
!  ANAME  Input, character*(*) ANAME, the name of the quantity
!         to be contoured.
!
!  ICOLOR Input, integer ( kind = 4 ) ICOLOR(MAXOBJ).
!         Contains the color indexes for each object.
!         However, in some cases, ICOLOR is actual a color table
!         index.
! 
!  JCMAX,
!  JCMIN  Input, integer ( kind = 4 ) JCMAX, JCMIN, the maximum and 
!         minimum color indices to use in the color bar.
!
!  L      Input, integer ( kind = 4 ) L, the number of rows of data.
!
!  LBAR   Input, logical LBAR, is .TRUE. if the color bar should
!         be shown.
!
!  M      Input, integer ( kind = 4 ) M, the number of columns of data.
!
!  MAXL,
!  MAXM   Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!         rows and columns of data.
!
!  MAXOBJ Input, integer ( kind = 4 ) MAXOBJ.
!         The number of graphical "objects".
! 
!  NCON   Input, integer ( kind = 4 ) NCON.
!         The number of contour lines to be drawn.
! 
!  NFLAG  Input, logical NFLAG(MAXL,MAXM).
!
!         NFLAG is used to "flag" which nodes are active,
!         that is, to be displayed, and which not, in the graph.
!
!  SRANGE Input, real SRANGE.
!         The maximum of XSMAX-XSMIN and YSMAX-YSMIN.
!         This gives the size of a square containing the data
!         window.
!
!  X      Input, real X(MAXL,MAXM).
!         The X coordinates of the nodes.
! 
!  X2MAX,
!  X2MIN  Input, real X2MAX, X2MIN, the maximum and minimum X 
!         coordinates that should be used for plotting.  No plotting 
!         commands should exceed these values.  This is where the 
!         "frame" might be drawn.
!
!  Y      Input, real Y(MAXL,MAXM).
!         The Y coordinates of the nodes.
! 
!  Y1MAX  Input, real Y1MAX, the maximum Y coordinates of the plot, 
!         which includes a small grace margin.
!
!  Y2MAX  Input, real Y2MAX, the maximum Y coordinates that should be 
!         used for plotting. 
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm
  integer ( kind = 4 ) maxobj

  real a(maxl,maxm)
  character ( len = * )  aname
  logical echo
  integer ( kind = 4 ) icolor(maxobj)
  integer ( kind = 4 ) jcmax
  integer ( kind = 4 ) jcmin
  integer ( kind = 4 ) l
  logical lbar
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ncon
  logical nflag(maxl,maxm)
  real smax
  real smin
  real srange
  real xp(maxl,maxm)
  real x1
  real x2
  real x2max
  real x2min
  real yp(maxl,maxm)
  real y1
  real y1max
  real y2
  real y2max
!
!  Get the minimum and maximum values of the data.
!
  call fsize(l,m,maxl,maxm,nflag,a,smax,smin)
!
!  If there is no variation in the data, don't try to draw a 
!  contour plot.
!
  if ( smax <= smin ) then
    write ( *, * ) ' '
    write ( *, * ) 'COLCON - Warning!'
    write ( *, * ) '  No color contours, all values equal.'
    return
  end if
!
!  Allow the user to adjust the range of the data.
!
  call setsiz(echo,aname,smax,smin)
!
!  Draw the color contour plot.
!
  call tricol(jcmax,jcmin,l,m,maxl,maxm,ncon,nflag,a,smax,smin,xp,yp)
!
!  Draw the color bar.
!
  if ( lbar ) then
    x1 = x2min
    x2 = x2max
    y1 = y2max
    y2 = y1max
    call cbar(icolor,jcmax,jcmin,maxobj,ncon,smax,smin,srange,x1,x2,y1,y2)
  end if

  return
end
subroutine cross(px,py,qx,qy,sl,sm,sh,sval,xl,xm,xh,yl,ym,yh)

!*****************************************************************************80
!
!! CROSS finds two places where a given value occurs on a triangle.  
!
!  Discussion:
!
!    The corners of the triangle are (XL,YL), (XM,YM) and
!    (XH,YH), and the associated S values are SL, SM and SH.  It 
!    must be the case that SL < =  SM <= SH.
!
!    CROSS returns two points:
!
!     (PX,PY), which occurs on one of the two sides that include 
!           (XM,YM), and 
!     (QX,QY), which occurs on the side between (XL,YL) and (XH,YH).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  PX,
!  PY     Output, real PX, PY, the X and Y coordinates of a point
!         at which S = SVAL, lying on a side of the triangle which 
!         ends at (XM,YM).
!
!  QX,
!  QY     Output, real QX, QY, the X and Y coordinates of a point
!         at which S = SVAL, lying on the side of the triangle which
!         lies between (XL,YL) and (XH,YH).
!
!  SL,
!  SM,
!  SH     Input, real SL, SM, SH, the low, medium, and high values
!         of S, associated with the three corners.
!
!  SVAL   Input, real SVAL, the value of S for which a contour line
!         is sought.
!
!  XL,
!  XM,
!  XH     Input, real XL, XM, XH, the X coordinates of the nodes
!         at which the low, medium and high values of S occur.
!
!  YL,
!  YM,
!  YH     Input, real YL, YM, YH, the Y coordinates of the nodes
!         at which the low, medium and high values of S occur.
!
  implicit none

  real px
  real py
  real qx
  real qy
  real sl
  real sm
  real sh
  real sval
  real xl
  real xm
  real xh
  real yl
  real ym
  real yh

  if ( sval < sl ) then

    px = xl
    py = yl
    qx = xl
    qy = yl

  else if ( sval >= sh ) then

    px = xh
    py = yh
    qx = xh
    qy = yh

  else

    if ( sval < sm ) then
      px = xl+(sval-sl)*(xm-xl)/(sm-sl)
      py = yl+(sval-sl)*(ym-yl)/(sm-sl)
    else
      px = xm+(sval-sm)*(xh-xm)/(sh-sm)
      py = ym+(sval-sm)*(yh-ym)/(sh-sm)
    end if

    qx = xl+(sval-sl)*(xh-xl)/(sh-sl)
    qy = yl+(sval-sl)*(yh-yl)/(sh-sl)

  end if

  return
end
subroutine delete(filnam)

!*****************************************************************************80
!
!! DELETE deletes a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*(*) FILNAM, the file to be deleted.
!
  implicit none

  character ( len = * )  filnam

  open(unit = 99,file=filnam,status='old',err=10)
  close(unit = 99,status='delete',err=10)
10    continue

  return
end
subroutine diamnd(xcentr,ycentr,radius,filled)

!*****************************************************************************80
!
!! DIAMND may be used to draw an open or filled diamond of a given size.
!
!  Discussion:
!
!    DIAMND calls PLYLIN to draw an open diamond, or PLYGON to draw a 
!    filled diamond.
!
!    To control the color of the diamond, simply call the DRAWCGM
!    routine LINCLR before drawing open diamonds or FILCOR before
!    drawing closed diamonds.
!
!    If the X and Y coordinate systems do not have the same scale,
!    the diamond drawn will be "flattened".  This can happen if the
!    routine SETWCD has been used to set the X and Y dimensions to
!    different extents.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XCENTR, YCENTR, the X and Y coordinates of
!    the center of the diamond.
!
!    Input, real RADIUS, the radius of the diamond.
!
!    Input, logical FILLED, .TRUE. if the diamond is to be filled.
!
  implicit none

  integer ( kind = 4 ), parameter :: npts = 5

  logical filled
  real radius
  real xcentr
  real xpoint(npts)
  real ycentr
  real ypoint(npts)

  xpoint(1) = xcentr+radius
  ypoint(1) = ycentr

  xpoint(2) = xcentr
  ypoint(2) = ycentr+radius

  xpoint(3) = xcentr-radius
  ypoint(3) = ycentr

  xpoint(4) = xcentr
  ypoint(4) = ycentr-radius

  xpoint(5) = xcentr+radius
  ypoint(5) = ycentr

  if ( filled ) then
    call plygon(npts-1,xpoint,ypoint)
  else
    call plylin(npts,xpoint,ypoint)
  end if

  return
end
subroutine double(e,gamt,jcrys1,jcrys2,kcrys,kmelt,kvoid,l,m, &
  maxbot,maxl,maxm,nbot,ncflag,npflag,p,pc,psi,rueta,ruksi, &
  t,te,tk,u,v,vmag,vort,vpflag,w,xbot,xc,xp,ybot,yc,yp)

!*****************************************************************************80
!
!! DOUBLE "reflects" the data around the Y axis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real E(MAXL,MAXM), the magnetic stream function.
!
!    Input/output, real GAMT(MAXL,MAXM), the diffusion coefficient.
!
!    Input/output, integer ( kind = 4 ) JCRYS1, the J coordinate of the column of 
!    corner nodes that define the left side of the crystal.
!
!    Input/output, integer ( kind = 4 ) JCRYS2, the J coordinate of the column of 
!    corner nodes that define the right side of the crystal.
!
!    Input/output, integer ( kind = 4 ) KCRYS(MAXL,MAXM), 
!    0, if control volume (I,J) is away from the crystal.
!    1, if control volume (I,J) is on the external crystal boundary.
!    2, if control volume (I,J) is on the internal crystal boundary.
!    3, if control volume (I,J) is in the crystal interior.
!
!    Input/output, integer ( kind = 4 ) KMELT(MAXL,MAXM),
!    0, if control volume (I,J) is away from the melt.
!    1, if control volume (I,J) is on the external melt boundary.
!    2, if control volume (I,J) is on the internal melt boundary.
!    3, if control volume (I,J) is in the melt interior.
! 
!    Input/output, integer ( kind = 4 ) KVOID(MAXL,MAXM), 
!    0, if control volume (I,J) is away from the void.
!    1, if control volume (I,J) is on the external void boundary.
!    2, if control volume (I,J) is on the internal void boundary.
!    3, if control volume (I,J) is in the void interior.
! 
!    Input, integer ( kind = 4 ) L, the number of rows of data.
!
!    Input/output, integer ( kind = 4 ) M, the number of columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
!    Input/output, real P(MAXL,MAXM), the pressure.
!
!    Input/output, real PC(MAXL,MAXM), the corrected pressure.
!
!    Input/output, real PSI(MAXL,MAXM), the stream function.
!
!    Input/output, real RUETA(MAXL,MAXM), RUKSI(MAXL,MAXM), the
!    the momentum in the ETA and KSI directions.
!
!    Input/output, real T(MAXL,MAXM), the temperature.
!
!    Input/output, real TE(MAXL,MAXM), the turbulent epsilon.
!
!    Input/output, real TK(MAXL,MAXM), the turbulent K.
!
!    Input/output, real U(MAXL,MAXM), the horizontal velocity.
!
!    Input/output, real V(MAXL,MAXM), the vertical velocity.
!
!    Input/output, real W(MAXL,MAXM), the axial velocity.
!
!    Input/output, real X(MAXL,MAXM), the X coordinates of the primary
!    nodes.
!
!    Input/output, real XC(MAXL,MAXM), the X coordinates of the corner
!    nodes.
!
!    Input/output, real Y(MAXL,MAXM), the Y coordinates of the primary
!    nodes.
!
!    Input/output, real YC(MAXL,MAXM), the Y coordinates of the corner
!    nodes.
!
  implicit none

  integer ( kind = 4 ) maxbot
  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  real e(maxl,maxm)
  real gamt(maxl,maxm)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys1
  integer ( kind = 4 ) jcrys2
  integer ( kind = 4 ) kcrys(maxl,maxm)
  integer ( kind = 4 ) kmelt(maxl,maxm)
  integer ( kind = 4 ) kvoid(maxl,maxm)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nbot
  logical ncflag(maxl,maxm)
  logical npflag(maxl,maxm)
  real p(maxl,maxm)
  real pc(maxl,maxm)
  real psi(maxl,maxm)
  real rueta(maxl,maxm)
  real ruksi(maxl,maxm)
  real t(maxl,maxm)
  real te(maxl,maxm)
  real tk(maxl,maxm)
  real u(maxl,maxm)
  real v(maxl,maxm)
  real vmag(maxl,maxm)
  real vort(maxl,maxm)
  logical vpflag(maxl,maxm)
  real w(maxl,maxm)
  real xbot(maxbot)
  real xc(maxl,maxm)
  real xp(maxl,maxm)
  real ybot(maxbot)
  real yc(maxl,maxm)
  real yp(maxl,maxm)

  if ( 2*m-1 > maxm ) then
    write ( *, * ) ' '
    write ( *, * ) 'DOUBLE - Fatal error!'
    write ( *, * ) '  The reflected data has M = ',2*m-1
    write ( *, * ) '  but CRYSTAL_PLOT has MAXM  =  ',maxm
    stop
  end if

  call rduble(e,l,m,maxl,maxm,'e')
  call rduble(gamt,l,m,maxl,maxm,'gamt')
  call lduble(ncflag,l,m,maxl,maxm)
  call lduble(npflag,l,m,maxl,maxm)
  call rduble(p,l,m,maxl,maxm,'p')
  call rduble(pc,l,m,maxl,maxm,'pc')
  call rduble(psi,l,m,maxl,maxm,'psi')
  call rduble(rueta,l,m,maxl,maxm,'rueta')
  call rduble(ruksi,l,m,maxl,maxm,'ruksi')
  call rduble(t,l,m,maxl,maxm,'t')
  call rduble(te,l,m,maxl,maxm,'te')
  call rduble(tk,l,m,maxl,maxm,'tk')
  call rduble(v,l,m,maxl,maxm,'v')
  call rduble(u,l,m,maxl,maxm,'u')
  call rduble(vort,l,m,maxl,maxm,'vort')
  call lduble(vpflag,l,m,maxl,maxm)
  call rduble(w,l,m,maxl,maxm,'w')
  call rduble(yp,l,m,maxl,maxm,'yp')
  call rduble(yc,l,m,maxl,maxm,'yc')
  call rduble(xp,l,m,maxl,maxm,'xp')
  call rduble(xc,l,m,maxl,maxm,'xc')
  call iduble(kcrys,l,m,maxl,maxm)
  call iduble(kmelt,l,m,maxl,maxm)
  call iduble(kvoid,l,m,maxl,maxm)
!
!  Take care of data defining the crystal extent.
!
  j = jcrys2
  jcrys1 = m-j+1
  jcrys2 = m+j-1
!
!  Update M to its new value.
!
  m = 2*m-1
!
!  Take care of data associated with the crucible shape.
!
  do i = 2*nbot-1,nbot,-1
    xbot(i) = xbot(i+1-nbot)
    ybot(i) = ybot(i+1-nbot)
  end do

  do i = 1,nbot-1
    xbot(i) = -xbot(2*nbot-i)
    ybot(i) = ybot(2*nbot-i)
  end do

  nbot = 2*nbot-1
!
!  Set the velocity magnitudes.
!
  do i = 1,l
    do j = 1,m
      vmag(i,j) = sqrt(u(i,j)**2+v(i,j)**2)
    end do
  end do

  return
end
subroutine dshlin(n,x,y,dshsiz)

!*****************************************************************************80
!
!! DSHLIN draws a dashed line connecting a series of points.
!
!  Discussion:
!
!    If the X and Y coordinates use different scale
!    factors, then dashes at different angles will seem to
!    have different lengths.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points to be connected.
!
!    Input, real X(N), Y(N), the X and Y coordinates of the
!    points.
!
!    Input, real DSHSIZ, the length, in the X, Y coordinate
!    system, of the dashed lines.  If it is negative or zero,
!    an error occurs.
!
  implicit none

  integer ( kind = 4 ) n

  real dist
  real dshsiz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ndash
  real x(n)
  real xnew
  real xold
  real xxnew
  real xxold
  real y(n)
  real ynew
  real yold
  real yynew
  real yyold
!
!  Make sure that DSHSIZ is positive.
!
  if ( dshsiz <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DshLin - Fatal error!'
    write ( *, * ) '  The parameter DSHSIZ must be positive.'
    write ( *, * ) '  but the input value is ',dshsiz
    return
  end if

  xnew = x(1)
  ynew = y(1)
  
  do i = 2,n
  
    xold = xnew
    yold = ynew
    xnew = x(i)
    ynew = y(i)
    dist = sqrt( (xnew-xold)**2 + (ynew-yold)**2 )
    
    if ( dist > 0.0 ) then
    
      ndash = int((dist/dshsiz)+1)
      
      if ( mod(ndash,4) /= 0 ) then
        ndash = ndash+(4-mod(ndash,4))
      end if
      
      if ( ndash <= 3)ndash = 4
!
!  Desired pattern is:
!
!  X0 - dash - blank - blank - dash - dash - blank - blank - dash - X1
!
      do j = 1,ndash
      
        if ( mod(j,4) == 0 .or. mod(j,4) == 1 ) then
          xxold = ( (ndash+1-j)*xold + (j-1)*xnew) / real(ndash)
          yyold = ( (ndash+1-j)*yold + (j-1)*ynew) / real(ndash)
          xxnew = ( (ndash-j)*xold + j*xnew) / real(ndash)
          yynew = ( (ndash-j)*yold + j*ynew) / real(ndash)
          call movcgm(xxold,yyold)
          call drwcgm(xxnew,yynew)
        end if
        
      end do
      
    end if
    
  end do

  return
end
subroutine flushl(string)

!*****************************************************************************80
!
!! FLUSHL flushes a string left.  
!
!  Example:
!
!    Input             Output
!
!    '     Hello'      'Hello     '
!    ' Hi there!  '    'Hi there!   '
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  STRING Input/output, character*(*) STRING.
!
!         On input, STRING is a string of characters.
!
!         On output, any initial blank characters in STRING
!         have been cut, and pasted back onto the end.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nonb
  character null
  character ( len = * )  string
!
!  Check the length of the string to the last nonblank.
!  If nonpositive, return.
!
  lchar = len_trim (string)
  if ( lchar <= 0)return
  null = char(0)
!
!  Find the occurrence of the first nonblank.
!
  do i = 1,lchar
    nonb = i
    if ( string(i:i) /= ' ' .and. string(i:i) /= null)go to 10
  end do

  return

10    continue
!
!  Shift the string left.
!
  do i = 1,lchar+1-nonb
    string(i:i) = string(i+nonb-1:i+nonb-1)
  end do
!
!  Blank out the end of the string.
!
  do i = lchar+2-nonb,lchar
    string(i:i) = ' '
  end do

  return
end
subroutine fsize(l,m,maxl,maxm,nflag,r,rmax,rmin)

!*****************************************************************************80
!
!! FSIZE computes the range of a real array defined at primary nodes.  
!
!  Discussion:
!
!    Only "flagged" values are considered.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
!    Input, logical NFLAG(MAXL,MAXM).
!    NFLAG is used to "flag" which nodes are active,
!    that is, to be displayed, and which not, in the graph.
!
!    Input, real R(MAXL,MAXM), the array to be checked.
!
!    Output, real RMAX, RMIN, the maximum and minimum values
!    of R.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical nflag(maxl,maxm)
  real r(maxl,maxm)
  real rmax
  real rmin
  logical started

  started = .false.
  rmax = 0.0
  rmin = 0.0
 
  do i = 1,l
    do j = 1,m
 
      if ( nflag(i,j) ) then
 
        if ( started ) then
          rmin = min(rmin,r(i,j))
          rmax = max(rmax,r(i,j))
        else
          started = .true.
          rmin = r(i,j)
          rmax = r(i,j)
        end if
 
      end if
 
    end do
  end do

  return
end
subroutine gettab(dev,echo,filgrf,grace,icmax,icmin,ierror,iplot,itable, &
  lgopen,ovrlay)

!*****************************************************************************80
!
!! GETTAB gets the color table choice from the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*10 DEV, the graphics output device to be used.  
!    Current legal values for DEV include:
! 
!      cgmb - CGM binary file.
!      ps   - PostScript file.
!      xws  - X window screen (interactive).
!
!    Input, character*80 FILGRF, the name of the output
!    graphics file.
!
!    Input, real GRACE.
!    The size of the "grace" margin on the plot.
!
!    Input, integer ( kind = 4 ) ICMAX, ICMIN, the maximum and minimum color 
!    indices to use in color contour graphics.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer ( kind = 4 ) IPLOT, the number of plots made so far.
! 
!    Output, integer ( kind = 4 ) ITABLE, the desired color table.
!    1: low black to high white
!    2: low blue to high yellow
!    3: low red, high blue, with bands between.
!    4: low red, yellow, green, blue, high white.
!    5: low white, blue, green, yellow, high red.
!    6: low blue to high red
!    7: linear table between 2 user colors.
!    8: linear table between N user colors.
!    9: low white to high black.
! 
!    Input, logical OVRLAY.
!    If OVRLAY is true, then the next time that a plot is
!    requested, a "new frame" command is suppressed, so that
!    the new plot is shown on top of the previous one.
!
  implicit none

  character ( len = 10 ) dev
  logical echo
  character ( len = 80 ) filcol
  character ( len = 80 ) filgrf
  real grace
  integer ( kind = 4 ) icmax
  integer ( kind = 4 ) icmin
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) itable
  logical lgopen
  logical ovrlay
  real x1max
  real x1min
  real y1max
  real y1min

  if ( itable == 0 ) then
    write ( *, * ) 'Enter name of color table file.'
    read(*,'(a)',err = 90,end=90)filcol
    write(17,'(a)')filcol
    if ( echo ) then
      write(*,'(a)')filcol
    end if
    call getctb(icmin,icmax,filcol,ierror)
  end if

  call preplt(dev,echo,filgrf,icmax,icmin,iplot,itable,lgopen,ovrlay)

  call settab(echo,icmax,icmin,itable)

  if ( dev == 'xws' ) then

    call cbox(grace,x1max,x1min,y1max,y1min)

    call buzz(dev,x1max,x1min,y1max,y1min)

  end if

  return

90    continue
  write ( *, * ) ' '
  write ( *, * ) 'GetTab - Error'
  ierror = 1
  return
end
subroutine getwin(echo,grace,srange,xmax,xmin,x1max,x1min,x2max, &
  x2min,xsmax,xsmin,ymax,ymin,y1max,y1min,y2max,y2min,ysmax,ysmin)

!*****************************************************************************80
!
!! GETWIN responds to the "W" command by telling the user
!  the current window, and getting the new value from the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real GRACE.
!    The size of the "grace" margin on the plot.
!
!  SRANGE Output, real SRANGE.
!         The maximum of XSMAX-XSMIN and YSMAX-YSMIN.
!         This gives the size of a square containing the data
!         window.
!
!  XMAX   Input, real XMAX.
!         The maximum X coordinate of all the nodes.
!         The maximum entry in the XC array.
! 
!  XMIN   Input, real XMIN.
!         The minimum X coordinate of all the nodes.
!         The minimum entry in the XC array.
!   
!  X1MAX,
!  X1MIN  Output, real X1MAX, X1MIN, the maximum and minimum X 
!         coordinates of the plot, which includes a small grace margin.
!
!  X2MAX,
!  X2MIN  Input/output, real X2MAX, X2MIN, the maximum and minimum X 
!         coordinates that should be used for plotting.  No plotting 
!         commands should exceed these values.  This is where the 
!         "frame" might be drawn.
!
!  XSMAX  Input/output, real XSMAX.
!         The maximum X coordinate of the data to be displayed.
!         XSMAX defaults to XMAX, but can be made smaller to
!         focus on a portion of the region.
! 
!  XSMIN  Input/output, real XSMIN.
!         The minimum X coordinate of the data to be displayed.
!         XSMIN defaults to XMIN, but can be made larger to
!         focus on a portion of the region.
! 
!  Y1MAX,
!  Y1MIN  Output, real Y1MAX, Y1MIN, the maximum and minimum Y 
!         coordinates of the plot, which includes a small grace margin.
!
!  Y2MAX,
!  Y2MIN  Input/output, real Y2MAX, Y2MIN, the maximum and minimum Y 
!         coordinates that should be used for plotting.  No plotting commands 
!         should  exceed these values.  This is where the "frame" might be 
!         drawn.
!
!  YMAX   Input, real YMAX.
!         The maximum Y coordinate of all the nodes.
!         The maximum value attained by the YC array.
! 
!  YMIN   Input, real YMIN.
!         The minimum Y coordinate of all the nodes.
!         The minimum value attained by the YC array.
!  
!  YSMAX  Input/output, real YSMAX.
!         The maximum Y coordinate of the data to be displayed.
!         YSMAX defaults to YMAX, but can be made smaller to
!         focus on a portion of the region.
! 
!  YSMIN  Input/output, real YSMIN.
!         The minimum Y coordinate of the data to be displayed.
!         YSMIN defaults to YMIN, but can be made larger to
!         focus on a portion of the region.
!
  implicit none

  logical echo
  real grace
  real srange
  real x1max
  real x1min
  real x2max
  real x2min
  real xmax
  real xmin
  real xsmax
  real xsmin
  real y1max
  real y1min
  real y2max
  real y2min
  real ymax
  real ymin
  real ysmax
  real ysmin
!
10    continue
  write ( *, * ) ' '
  write ( *, * ) 'GetWin:'
  write ( *, * ) ' '
  write ( *, * ) '  Data coordinates:'
  write ( *, * ) ' '
  write ( *, * ) xmin,'   =  XMIN  <= X <= XMAX =  ',xmax
  write ( *, * ) ymin,'   =  YMIN  <= Y <= YMAX =  ',ymax
  write ( *, * ) ' '
  write ( *, * ) '  Displayed data coordinates:'
  write ( *, * ) ' '
  write ( *, * ) xsmin,'  =  XSMIN <= X <= XSMAX = ',xsmax
  write ( *, * ) ysmin,'  =  YSMIN <= Y <= YSMAX = ',ysmax
  write ( *, * ) ' '
  write ( *, * ) '  Picture window coordinates:'
  write ( *, * ) ' '
  write ( *, * ) x2min,'  =  X2MIN <= X <= X2MAX = ',x2max
  write ( *, * ) y2min,'  =  Y2MIN <= Y <= Y2MAX = ',y2max
  write ( *, * ) ' '
  write ( *, * ) '  Enter new displayed data coordinates for X and Y:'
  write ( *, * ) '  Use the order XSMIN, XSMAX, YSMIN, YSMAX:'

  read(*,*,err = 20,end=20)xsmin,xsmax,ysmin,ysmax
  write(17,*)xsmin,xsmax,ysmin,ysmax
  if ( echo ) then
    write ( *, * ) xsmin,xsmax,ysmin,ysmax
  end if
!
!  Compute box containing data.
!
  call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax,xsmin, &
    y1max,y1min,y2max,y2min,ysmax,ysmin)

  return

20    continue
  write ( *, * ) ' '
  write ( *, * ) 'GetWin - Warning!'
  write ( *, * ) '  Your input was not acceptable!'
  write ( *, * ) '  Please try again!'
  write ( *, * ) ' '
  go to 10

end
subroutine graph(delx,dely,dev,dxcdp,dycdp,echo,filgrf,gamt,icmax,icmin, &
  icolor,icrys1,icrys2,iplot,itable,jcmax,jcmin,jcrys1,jcrys2,kcrys,kmelt, &
  kvoid,l,lbar,lflag,lgopen,m,maxbot,maxl,maxm,maxobj,nbot,ncflag,ncon, &
  nflag,npflag,nxskip,nyskip,object,ovrlay,p,pc,psi,s1,s2,scalecv,scalenb, &
  scalenc,scalenp,scalev,show,srange,t,te,title,title2,tk,u,v,vpflag,vmag, &
  vort,w,x1max,x1min,x2max,x2min,xbot,xc,xp,xsmax,xsmin,y1max,y1min, &
  y2max,y2min,ybot,yc,yp,ysmax,ysmin)

!*****************************************************************************80
!
!! GRAPH draws the graph, after the user has specified what is to be shown.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real DELX.
!    The X spacing between nodes.  In some cases,
!    this spacing is modified to create isoparametric elements.
! 
!  DELY   Input, real DELY.
!         The Y spacing between nodes.  In some cases,
!         this spacing is modified to create isoparametric elements.
! 
!  DEV    Input, character*10 DEV.
!         The graphics output device to be used.  Current legal
!         values for DEV include:
! 
!         cgmb - CGM binary file.
!         ps   - PostScript file.
!         xws  - X window screen (interactive).
!
!  DXCDP  Input, real DXCDP(MAXL,MAXM), the difference between
!         the second and first values of XC.
!
!  DYCDP  Input, real DYCDP(MAXL,MAXM), the difference between
!         the second and first values of YC.
!
!  FILGRF Input, character*80 FILGRF, the name of the output
!         graphics file.
!
!  GAMT   Output, real GAMT(MAXL,MAXM), the diffusion coefficient.
!
!  ICMAX,
!  ICMIN  Input, integer ( kind = 4 ) ICMAX, ICMIN, the maximum and 
!         minimum color indices to use in the color bar.
!
!  ICOLOR Input, integer ( kind = 4 ) ICOLOR(MAXOBJ).
!         Contains the color indexes for each object.
!         However, in some cases, ICOLOR is actual a color table
!         index.
!
!  ICRYS1 Input, integer ( kind = 4 ) ICRYS1, the I coordinate of the row of corner
!         nodes that define the bottom of the crystal.
!
!  ICRYS2 Input, integer ( kind = 4 ) ICRYS2, the I coordinate of the row of corner
!         nodes that define the top of the crystal.
!
!  IPLOT  Input, integer ( kind = 4 ) IPLOT.
!         The number of plots made so far.
! 
!  ITABLE Input, integer ( kind = 4 ) ITABLE, the desired color table.
!  
!         1: low black to high white
!         2: low blue to high yellow
!         3: low red, high blue, with bands between.
!         4: low red, yellow, green, blue, high white.
!         5: low white, blue, green, yellow, high red.
!         6: low blue to high red
!         7: linear table between 2 user colors.
!         8: linear table between N user colors.
!         9: low white to high black.
! 
!
!  JCMAX,
!  JCMIN  Input, integer ( kind = 4 ) JCMAX, JCMIN, the maximum and 
!         minimum color indices to use in the color bar.
!
!  JCRYS1 Input/output, integer ( kind = 4 ) JCRYS1, the J coordinate of the column of 
!         corner nodes that define the left side of the crystal.
!
!  JCRYS2 Input/output, integer ( kind = 4 ) JCRYS2, the J coordinate of the column of 
!         corner nodes that define the right side of the crystal.
!
!  KCRYS  Input, integer ( kind = 4 ) KCRYS(MAXL,MAXM), 
!         0, if control volume (I,J) is away from the crystal.
!         1, if control volume (I,J) is on the external crystal boundary.
!         2, if control volume (I,J) is on the internal crystal boundary.
!         3, if control volume (I,J) is in the crystal interior.
!
!  LFLAG  Workspace, logical LFLAG(MAXL,MAXM), .TRUE. if primary node
!         (I,J) lies within the visible data, .FALSE. otherwise.
!
!  KMELT  Input, integer ( kind = 4 ) KMELT(MAXL,MAXM),
!         0, if control volume (I,J) is away from the melt.
!         1, if control volume (I,J) is on the external melt boundary.
!         2, if control volume (I,J) is on the internal melt boundary.
!         3, if control volume (I,J) is in the melt interior.
! 
!  KVOID  Input, integer ( kind = 4 ) KVOID(MAXL,MAXM), 
!         0, if control volume (I,J) is away from the void.
!         1, if control volume (I,J) is on the external void boundary.
!         2, if control volume (I,J) is on the internal void boundary.
!         3, if control volume (I,J) is in the void interior.
! 
!  L      Input, integer ( kind = 4 ) L, the number of rows of data.
!
!  LBAR   Input, logical LBAR, is .TRUE. if the color bar should
!         be shown.
!
!  M      Input, integer ( kind = 4 ) M, the number of columns of data.
!
!  MAXBOT Input, integer ( kind = 4 ) MAXBOT, the maximum number of crucible nodes
!         allowed.
!
!  MAXL,
!  MAXM   Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!         rows and columns of data.
!
!  MAXOBJ Input, integer ( kind = 4 ) MAXOBJ.
!         The number of graphical "objects".
! 
!  NBOT   Input, integer ( kind = 4 ) NBOT, the number of crucible nodes.
!
!  NCON   Input, integer ( kind = 4 ) NCON.
!         The number of contour lines to be drawn.  This is
!         initialized to 12, but may be changed by the user.
! 
!  NFLAG  Input, logical NFLAG(MAXL,MAXM).
!
!         NFLAG is used to "flag" which nodes are active,
!         that is, to be displayed, and which not, in the graph.
!
!  NXSKIP Input, integer ( kind = 4 ) NXSKIP.
!         NXSKIP is used to "thin" out a vector plot.  
!
!         If NXSKIP = 1, then a standard vector plot is made.
!  
!         Otherwise, in the X direction, vectors are drawn only
!         in columns 1, 1+NXSKIP, 1+2*NXSKIP and so on.
!
!  NYSKIP Input, integer ( kind = 4 ) NYSKIP.
!         NYSKIP is used to "thin" out a vector plot.  
!
!         If NYSKIP = 1, then a standard vector plot is made.
!
!         Otherwise, in the Y direction, vectors are drawn only
!         in rows 1, 1+NYSKIP, 1+2*NYSKIP and so on.
!  
!  OBJECT Input, character*25 OBJECT(MAXOBJ), the names of the 
!         graphical objects.
!
!  OVRLAY Input, logical OVRLAY.
!         If OVRLAY is true, then the next time that a plot is
!         requested, a "new frame" command is suppressed, so that
!         the new plot is shown on top of the previous one.
! 
!  P      Input, real P(MAXL,MAXM), the pressure.
!
!  PC     Input, real PC(MAXL,MAXM), the corrected pressure.
!
!  PSI    Input, real PSI(MAXL,MAXM), the stream function.
!
!  S1     Workspace, real S1(MAXL,MAXM).
!
!  S2     Workspace, real S2(MAXL,MAXM).
!
!  SHOW   Input, logical SHOW(MAXOBJ).
!         Contains, for each object, a flag determining whether it
!         is to be shown or not.
! 
!  SRANGE Input, real SRANGE.
!         The maximum of XSMAX-XSMIN and YSMAX-YSMIN.
!         This gives the size of a square containing the data
!         window.
!
!  SCALECV Input, real SCALECV.
!         A scale factor for the control volume numbers, which has a
!         default value of 1.
!
!  SCALENB Input, real SCALENB.
!         A scale factor for the size of the nodes that define the
!         crucible shape, which has a default value of 1.
!
!  SCALENC Input, real SCALENC.
!         A scale factor for the size of the corner nodes, which has a
!         default value of 1.
!
!  SCALENP Input, real SCALENP.
!         A scale factor for the size of the primary nodes, which has a
!         default value of 1.
!
!  SCALEV Input, real SCALEV.
!         A scale factor for velocity vectors.  This starts out at 1.0.
! 
!  T      Input, real T(MAXL,MAXM), the temperature.
!
!  TE     Input, real TE(MAXL,MAXM), the turbulent epsilon.
!
!  TITLE  Input, character*40 TITLE.
!         A title for the plots.
! 
!  TITLE2 Input, character*40 TITLE2.
!         A subtitle used in the profile plots.
!
!  TK     Input, real TK(MAXL,MAXM), the turbulent K.
!
!  U      Input, real U(MAXL,MAXM), the horizontal velocity.
!
!  V      Input, real V(MAXL,MAXM), the vertical velocity.
!
!  VPFLAG Input, logical VPFLAG(MAXL,MAXM), a flag which determines
!         whether vector data at a given node should be displayed.
!
!  VMAG   Input, real VMAG(MAXL,MAXM), the velocity magnitude.
!
!  W      Input, real W(MAXL,MAXM), the axial velocity.
!
!  X      Output, real X(MAXL,MAXM), the X coordinates of the primary
!         nodes.
!
!  X1MAX,
!  X1MIN  Output, real X1MAX, X1MIN, the maximum and minimum X 
!         coordinates of the plot, which includes a small grace margin.
!
!  X2MAX,
!  X2MIN  Input, real X2MAX, X2MIN, the maximum and minimum X 
!         coordinates that should be used for plotting.  No plotting 
!         commands should exceed these values.  This is where the 
!         "frame" might be drawn.
!
!  XBOT   Input, real XBOT(NBOT), the X coordinates of the crucible
!         bottom nodes.
!
!  XC     Input, real XC(MAXL,MAXM), the X coordinates of the corner
!         nodes.
!
!  XSMAX  Input, real XSMAX.
!         The maximum X coordinate of the data to be displayed.
!         XSMAX defaults to XMAX, but can be made smaller to
!         focus on a portion of the region.
! 
!  XSMIN  Input, real XSMIN.
!         The minimum X coordinate of the data to be displayed.
!         XSMIN defaults to XMIN, but can be made larger to
!         focus on a portion of the region.
!
!  Y      Output, real Y(MAXL,MAXM), the Y coordinates of the primary
!         nodes.
!
!  Y1MAX,
!  Y1MIN  Output, real Y1MAX, Y1MIN, the maximum and minimum Y 
!         coordinates of the plot, which includes a small grace margin.
!
!  Y2MAX,
!  Y2MIN  Input, real Y2MAX, Y2MIN, the maximum and minimum Y 
!         coordinates that should be used for plotting.  No plotting 
!         commands should exceed these values.  This is where the 
!         "frame" might be drawn.
!
!  YBOT   Input, real YBOT(NBOT), the Y coordinates of the crucible
!         bottom nodes.
!
!  YC     Input, real YC(MAXL,MAXM), the Y coordinates of the corner
!         nodes.
!
!  YSMAX  Input, real YSMAX.
!         The maximum Y coordinate of the data to be displayed.
!         YSMAX defaults to YMAX, but can be made smaller to
!         focus on a portion of the region.
! 
!  YSMIN  Input, real YSMIN.
!         The minimum Y coordinate of the data to be displayed.
!         YSMIN defaults to YMIN, but can be made larger to
!         focus on a portion of the region.
!
  implicit none

  integer ( kind = 4 ) maxbot
  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm
  integer ( kind = 4 ) maxobj

  character ( len = 25 ) aname
  real angle
  character ( len = 6 ) chrtmp
  real csize
  real cwide
  real delx
  real dely
  character ( len = 10 ) dev
  real dshsiz
  real dxcdp(maxl,maxm)
  real dycdp(maxl,maxm)
  logical echo
  character ( len = 80 ) filgrf
  logical filled
  character flush
  real gamt(maxl,maxm)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icmax
  integer ( kind = 4 ) icmin
  integer ( kind = 4 ) icolor(maxobj)
  integer ( kind = 4 ) icrys1
  integer ( kind = 4 ) icrys2
  integer ( kind = 4 ) ido
  logical inside
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) itable
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcmax
  integer ( kind = 4 ) jcmin
  integer ( kind = 4 ) jcolor
  integer ( kind = 4 ) jcrys1
  integer ( kind = 4 ) jcrys2
  integer ( kind = 4 ) kcrys(maxl,maxm)
  integer ( kind = 4 ) kmelt(maxl,maxm)
  integer ( kind = 4 ) kvoid(maxl,maxm)
  integer ( kind = 4 ) l
  logical lbar
  logical lbar2
  integer ( kind = 4 ) lent
  logical lflag(maxl,maxm)
  logical lgopen
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nbot
  logical ncflag(maxl,maxm)
  integer ( kind = 4 ) ncon
  logical nflag(maxl,maxm)
  integer ( kind = 4 ) nonzer
  logical npflag(maxl,maxm)
  integer ( kind = 4 ) npts
  integer ( kind = 4 ) nval
  integer ( kind = 4 ) nxskip
  integer ( kind = 4 ) nyskip
  character ( len = 25 ) object(maxobj)
  logical ovrlay
  real p(maxl,maxm)
  real pc(maxl,maxm)
  real psi(maxl,maxm)
  real pwide
  real s1(maxl,maxm)
  real s2(maxl,maxm)
  real scalecv
  real scalenb
  real scalenc
  real scalenp
  real scalev
  logical show(maxobj)
  real srange
  real t(maxl,maxm)
  real te(maxl,maxm)
  character ( len = 40 ) title
  character ( len = 40 ) title2
  real tk(maxl,maxm)
  real u(maxl,maxm)
  real uu
  real uvmag
  real v(maxl,maxm)
  real velscl
  logical vpflag(maxl,maxm)
  real vmag(maxl,maxm)
  real vmax
  real vort(maxl,maxm)
  real vv
  real w(maxl,maxm)
  real x1max
  real x1min
  real x2max
  real x2min
  real xbot(maxbot)
  real xc(maxl,maxm)
  real xp(maxl,maxm)
  real xpoly(4)
  real xsmax
  real xsmin
  real xval
  real xvec(2)
  real xx(300)
  real y1max
  real y1min
  real y2max
  real y2min
  real ybot(maxbot)
  real yc(maxl,maxm)
  real yp(maxl,maxm)
  real ypoly(4)
  real ysmax
  real ysmin
  real yval
  real yvec(2)
  real yy(300)

  call preplt(dev,echo,filgrf,icmax,icmin,iplot,itable,lgopen,ovrlay)

  cwide = 0.9E+00 * srange / 40.0E+00
  pwide = srange
!
!  Set the scale of the picture.
!  We allow ourselves more or less a one percent margin.
!
  xvec(1) = x2min-0.01*(x2max-x2min)
  yvec(1) = y2min-0.01*(y2max-y2min)
  xvec(2) = x2max+0.01*(x2max-x2min)
  yvec(2) = y2max+0.01*(y2max-y2min)
 
  call setscl(xvec,yvec,2)
!
!  Mark those primary nodes which are visible.
!
!  This marker is based on visibility (inside the screen limits)
!  and on whether the user has asked to hide the crystal,
!  melt, or void regions.
!
  do i = 1,l
    do j = 1,m

      if ( (.not.inside(xsmin,xp(i,j),xsmax)) .or. &
           (.not.inside(ysmin,yp(i,j),ysmax)) ) then
        lflag(i,j) = .false.
      else if ( kcrys(i,j) == 3 .and. (.not.show(39)) ) then
        lflag(i,j) = .false.
      else if ( kmelt(i,j) == 3 .and. (.not.show(40)) ) then
        lflag(i,j) = .false.
      else if ( kvoid(i,j) == 3 .and. (.not.show(38)) ) then
        lflag(i,j) = .false.
      else
        lflag(i,j) = .true.
      end if

    end do
  end do
!
!  FRAME: Draw a box around the region.
!
  if ( show(3) ) then

    call linclr(icolor(3))
    call box(x2min,x2max,y2min,y2max)

  end if
!
!  TITLE/TITLE2: Draw the titles.
!
  if ( show(7) ) then

    call linclr(icolor(7))

    lent = len_trim ( title )

    if ( lent > 0 ) then
      angle = 0.0
      xval = 0.5*(xsmax+xsmin)
      yval = ysmax+3.0*cwide
      flush = 'center'
      call s_plot(angle,cwide,pwide,title,xval,yval,flush)
    end if

    lent = len_trim ( title2 )

    if ( lent > 0 ) then
      angle = 0.0
      xval = 0.5*(xsmax+xsmin)
      yval = ysmax+1.5*cwide
      flush = 'center'
      call s_plot(angle,cwide,pwide,title2,xval,yval,flush)
    end if

  end if
!
!  CPC: Draw corrected pressure color contours for nodes which are visible,
!  and which are in or on the melt region.
!
  if ( show(26) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    aname = object(26)
    call colcon(pc,aname,echo,icolor,jcmax,jcmin,l,lbar,m,maxl, &
      maxm,maxobj,ncon,nflag,srange,xp,x2max,x2min,yp,y1max,y2max)

  end if
!
!  CRYSC: color in the crystal.
!
  if ( show(44) ) then

    jcolor = icolor(44)
    call filclr(jcolor)

    do i = 1,l-1
      do j = 1,m-1

        if ( kcrys(i,j) >= 2 ) then

          npts = 4

          xpoly(1) = xc(i,j)
          ypoly(1) = yc(i,j)

          xpoly(2) = xc(i+1,j)
          ypoly(2) = yc(i+1,j)

          xpoly(3) = xc(i+1,j+1)
          ypoly(3) = yc(i+1,j+1)

          xpoly(4) = xc(i,j+1)
          ypoly(4) = yc(i,j+1)

          call plygon(npts,xpoly,ypoly)

        end if
      end do
    end do
  end if
!
!  MELTC: color in the melt.
!
  if ( show(45) ) then

    jcolor = icolor(45)
    call filclr(jcolor)

    do i = 1,l-1
      do j = 1,m-1

        if ( kmelt(i,j) >= 2 ) then

          npts = 4

          xpoly(1) = xc(i,j)
          ypoly(1) = yc(i,j)

          xpoly(2) = xc(i+1,j)
          ypoly(2) = yc(i+1,j)

          xpoly(3) = xc(i+1,j+1)
          ypoly(3) = yc(i+1,j+1)

          xpoly(4) = xc(i,j+1)
          ypoly(4) = yc(i,j+1)

          call plygon(npts,xpoly,ypoly)

        end if
      end do
    end do
  end if
!
!  PC: Draw pressure color contours for nodes which are visible,
!  and which are in or on the melt region.
!
  if ( show(12) ) then

    do i = 1,l
      do j = 1,m

        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if

      end do
    end do

    aname = object(12)
    call colcon(p,aname,echo,icolor,jcmax,jcmin,l,lbar,m,maxl, &
      maxm,maxobj,ncon,nflag,srange,xp,x2max,x2min,yp,y1max,y2max)

  end if
!
!  SC: Draw stream color contours in the melt region.
!
  if ( show(35) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    aname = object(35)
    call colcon(psi,aname,echo,icolor,jcmax,jcmin,l,lbar,m,maxl, &
      maxm,maxobj,ncon,nflag,srange,xp,x2max,x2min,yp,y1max,y2max)

  end if
!
!  TEMPC: Draw temperature color contours for nodes which are visible,
!  and which are in or on the melt or crystal regions.
!
  if ( show(13) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. (kmelt(i,j) >= 2 .or. kcrys(i,j) >= 2) ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    aname = object(13)
    call colcon(t,aname,echo,icolor,jcmax,jcmin,l,lbar,m,maxl, &
      maxm,maxobj,ncon,nflag,srange,xp,x2max,x2min,yp,y1max,y2max)

  end if
!
!  VMAGC: Draw velocity magnitude color contours for nodes which are visible,
!  and which are in or on the melt region.
!
  if ( show(37) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    aname = object(37)
    call colcon(vmag,aname,echo,icolor,jcmax,jcmin,l,lbar,m,maxl, &
      maxm,maxobj,ncon,nflag,srange,xp,x2max,x2min,yp,y1max,y2max)

  end if
!
!  VOIDC: color in the void.
!
  if ( show(46) ) then

    jcolor = icolor(46)
    call filclr(jcolor)

    do i = 1,l-1
      do j = 1,m-1

        if ( kvoid(i,j) >= 2 ) then

          npts = 4

          xpoly(1) = xc(i,j)
          ypoly(1) = yc(i,j)

          xpoly(2) = xc(i+1,j)
          ypoly(2) = yc(i+1,j)

          xpoly(3) = xc(i+1,j+1)
          ypoly(3) = yc(i+1,j+1)

          xpoly(4) = xc(i,j+1)
          ypoly(4) = yc(i,j+1)

          call plygon(npts,xpoly,ypoly)

        end if
      end do
    end do
  end if
!
!  VORTC: draw vorticity color contours in melt region.
!
  if ( show(42) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    aname = object(42)
    call colcon(vort,aname,echo,icolor,jcmax,jcmin,l,lbar,m,maxl, &
      maxm,maxobj,ncon,nflag,srange,xp,x2max,x2min,yp,y1max,y2max)

  end if
!
!  WC: Draw angular momentum color contours in melt region.
!
  if ( show(14) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    aname = object(14)
    call colcon(w,aname,echo,icolor,jcmax,jcmin,l,lbar,m,maxl, &
      maxm,maxobj,ncon,nflag,srange,xp,x2max,x2min,yp,y1max,y2max)

  end if
!
!  XCC: draw XC coordinate color contours.
!
  if ( show(31) ) then

    do i = 1,l
      do j = 1,m

        nflag(i,j) = ncflag(i,j)

        if ( .not.lflag(i,j) ) then
          nflag(i,j) = .false.
        end if

      end do
    end do

    aname = object(31)
    call colcon(xc,aname,echo,icolor,jcmax,jcmin,l,lbar,m,maxl, &
      maxm,maxobj,ncon,nflag,srange,xp,x2max,x2min,yp,y1max,y2max)

  end if
!
!  YCC: draw YC coordinate color contours.
!
  if ( show(33) ) then

    do i = 1,l
      do j = 1,m

        nflag(i,j) = ncflag(i,j)

        if ( .not.lflag(i,j) ) then
          nflag(i,j) = .false.
        end if

      end do
    end do

    aname = object(33)
    call colcon(yc,aname,echo,icolor,jcmax,jcmin,l,lbar,m,maxl, &
      maxm,maxobj,ncon,nflag,srange,xp,x2max,x2min,yp,y1max,y2max)

  end if
!
!  LINE DRAWINGS BEGIN HERE.
!
!
!  AX: draw the axis of symmetry.
!
  if ( show(28) ) then

    if ( xsmin <= 0.0 .and. 0.0 <= xsmax ) then

      call linclr(icolor(28))

      nval = 2
      xx(1) = 0.0
      yy(1) = ysmin
      xx(2) = 0.0
      yy(2) = ysmax

      dshsiz = 0.025
      call dshlin(nval,xx,yy,dshsiz)

    end if

  end if
!
!  BO: Draw the boundary.
!
  if ( show(1) ) then

    call linclr(icolor(1))

    dshsiz = 0.05
    nval = 2

    do i = 1,l-1

      if ( ncflag(i,1) .and. ncflag(i+1,1) ) then

        xx(1) = xc(i,1)
        yy(1) = yc(i,1)
        xx(2) = xc(i+1,1)
        yy(2) = yc(i+1,1)

        call dshlin(nval,xx,yy,dshsiz)

      end if

    end do

    do j = 1,m-1

      if ( ncflag(l,j) .and. ncflag(l,j+1) ) then

        xx(1) = xc(l,j)
        yy(1) = yc(l,j)
        xx(2) = xc(l,j+1)
        yy(2) = yc(l,j+1)

        call dshlin(nval,xx,yy,dshsiz)

      end if

    end do

    do i = 1,l-1

      if ( ncflag(i,m) .and. ncflag(i+1,m) ) then

        xx(1) = xc(i,m)
        yy(1) = yc(i,m)
        xx(2) = xc(i+1,m)
        yy(2) = yc(i+1,m)

        call dshlin(nval,xx,yy,dshsiz)

      end if

    end do

    do j = 1,m-1

      if ( ncflag(1,j) .and. ncflag(1,j+1) ) then

        xx(1) = xc(1,j)
        yy(1) = yc(1,j)
        xx(2) = xc(1,j+1)
        yy(2) = yc(1,j+1)

        call dshlin(nval,xx,yy,dshsiz)

      end if

    end do

  end if
!
!  CP: Draw corrected pressure contours in melt region.
!
  if ( show(25) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    if ( show(26) ) then
      ido = 0
      call linclr(icolor(25))
    else
      ido = 1
    end if

    aname = object(25)
    lbar2 = lbar .and. .not.show(26)

    call lincon(pc,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
      maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

  end if
!
!  CRU: Draw the crucible boundary.
!
  if ( show(19) ) then

    call linclr(icolor(19))

    i = 1
    do j = 1,m-1

      if ( xsmin <= xc(i,j) .and. xc(i,j) <= xsmax.and. &
           ysmin <= yc(i,j) .and. yc(i,j) <= ysmax.and. &
           xsmin <= xc(i,j+1) .and. xc(i,j+1) <= xsmax.and. &
           ysmin <= yc(i,j+1) .and. yc(i,j+1) <= ysmax ) then
        call movcgm(xc(i,j),yc(i,j))
        call drwcgm(xc(i,j+1),yc(i,j+1))
      end if

    end do
!
!  Do this loop only if reflected.
!
    if ( xc(1,1) == -xc(1,m) ) then

      j = 1
      do i = 1,icrys1-1

        if ( xsmin <= xc(i,j) .and. xc(i,j) <= xsmax.and. &
             ysmin <= yc(i,j) .and. yc(i,j) <= ysmax.and. &
             xsmin <= xc(i+1,j) .and. xc(i+1,j) <= xsmax.and. &
             ysmin <= yc(i+1,j) .and. yc(i+1,j) <= ysmax ) then
          call movcgm(xc(i,j),yc(i,j))
          call drwcgm(xc(i+1,j),yc(i+1,j))
        end if

      end do

    end if

    j = m
    do i = 1,icrys1-1

      if ( xsmin <= xc(i,j) .and. xc(i,j) <= xsmax.and. &
           ysmin <= yc(i,j) .and. yc(i,j) <= ysmax.and. &
           xsmin <= xc(i+1,j) .and. xc(i+1,j) <= xsmax.and. &
           ysmin <= yc(i+1,j) .and. yc(i+1,j) <= ysmax ) then
        call movcgm(xc(i,j),yc(i,j))
        call drwcgm(xc(i+1,j),yc(i+1,j))
      end if

    end do

  end if
!
!  CRYS: Draw the crystal.
!
  if ( show(20) ) then

    call linclr(icolor(20))

    call movcgm(xc(icrys1,jcrys1),yc(icrys1,jcrys1))

    do j = jcrys1+1,jcrys2
      if ( ncflag(icrys1,j-1) .and. ncflag(icrys1,j) ) then
        call drwcgm(xc(icrys1,j),yc(icrys1,j))
      end if
    end do

    do i = icrys1+1,icrys2
      if ( ncflag(i-1,jcrys2) .and. ncflag(i,jcrys2) ) then
        call drwcgm(xc(i,jcrys2),yc(i,jcrys2))
      end if
    end do

    do j = jcrys2-1,jcrys1,-1
      if ( ncflag(icrys2,j+1) .and. ncflag(icrys2,j) ) then
        call drwcgm(xc(icrys2,j),yc(icrys2,j))
      end if
    end do

    do i = icrys2,icrys1,-1
      if ( ncflag(i+1,jcrys1) .and. ncflag(i,jcrys1) ) then
        call drwcgm(xc(i,jcrys1),yc(i,jcrys1))
      end if
    end do

  end if
!
!  CV: Draw the control volumes.
!
  if ( show(2) ) then

    call linclr(icolor(2))

    do i = 1,l
      do j = 1,m

        if ( (.not.inside(xsmin,xc(i,j),xsmax)) .or. & 
             (.not.inside(ysmin,yc(i,j),ysmax)) ) then
          nflag(i,j) = .false.
        else
          nflag(i,j) = .true.
        end if

      end do
    end do

    do i = 1,l-1
      do j = 1,m-1

        if ( (kcrys(i,j) == 3 .and. (.not.show(39))) .or. &
             (kmelt(i,j) == 3 .and. (.not.show(40))) .or. &
             (kvoid(i,j) == 3 .and. (.not.show(38))) ) then
          nflag(i,j) = .false.
          nflag(i+1,j) = .false.
          nflag(i,j+1) = .false.
          nflag(i+1,j+1) = .false.
        end if

      end do
    end do

!       do i = 1,l
!         do j = 1,m
!           nflag(i,j) = lflag(i,j)
!         end do
!       end do

    call showcv(l,m,maxl,maxm,nflag,xc,yc)

  end if
!
!  CVN: Draw the control volume numbers.
!
  if ( show(43) ) then
    angle = 0.0
    cwide = 0.05*srange*scalecv
    pwide = srange
    flush = 'center'
    call linclr(icolor(43))
    do i = 1,l-1
      do j = 1,m-1

        if ( lflag(i,j) ) then

          write(chrtmp,'(i6)')i
          call flushl(chrtmp)

          call s_plot(angle,cwide,pwide,chrtmp,xp(i,j),yp(i,j),flush)

        end if

      end do
    end do

  end if
!
!  FS: Draw the free surface.
!
  if ( show(21) ) then

    call linclr(icolor(21))

    do i = 2,l-1
      do j = 1,m-1
        if (  kmelt(i,j) == 1 .and. kmelt(i-1,j) == 2) then

          if ( xsmin <= xc(i,j) .and. xc(i,j) <= xsmax .and. &
               ysmin <= yc(i,j) .and. yc(i,j) <= ysmax .and. &
               xsmin <= xc(i,j+1) .and. xc(i,j+1) <= xsmax .and. &
               ysmin <= yc(i,j+1) .and. yc(i,j+1) <= ysmax ) then
            call movcgm(xc(i,j),yc(i,j))
            call drwcgm(xc(i,j+1),yc(i,j+1))
          end if

        end if
      end do
    end do

  end if
!
!  NB: Draw the crucible nodes.
!
  if ( show(34) ) then

    call linclr(icolor(34))
    filled = .false.
    csize = 0.0125*srange*scalenb

    do i = 1,nbot
      if ( inside(xsmin,xbot(i),xsmax) .and. &
           inside(ysmin,ybot(i),ysmax) ) then
        call circle(xbot(i),ybot(i),csize,filled)
      end if
    end do

  end if
!
!  NC: Draw the corner nodes.
!
  if ( show(29) ) then

    call linclr(icolor(29))
    filled = .false.
    csize = 0.0125*srange*scalenc

    do i = 1,l
      do j = 1,m

        nflag(i,j) = ncflag(i,j)

        if ( mod(i-1,nxskip) /= 0 ) then
          nflag(i,j) = .false.
        end if

        if ( mod(j-1,nyskip) /= 0 ) then
          nflag(i,j) = .false.
        end if

        if ( nflag(i,j) ) then
          call diamnd(xc(i,j),yc(i,j),csize,filled)
        end if

      end do
    end do

  end if
!
!  NP: Draw the primary nodes.
!
  if ( show(4) ) then

    call linclr(icolor(4))
    filled = .false.
    csize = 0.0125*srange*scalenp

    do i = 1,l
      do j = 1,m

        nflag(i,j) = npflag(i,j)

        if ( mod(i-1,nxskip) /= 0 ) then
          nflag(i,j) = .false.
        end if

        if ( mod(j-1,nyskip) /= 0 ) then
          nflag(i,j) = .false.
        end if

        if ( nflag(i,j) ) then
          call circle(xp(i,j),yp(i,j),csize,filled)
        end if

      end do
    end do

  end if
!
!  P: Draw pressure contours in melt region.
!
  if ( show(5) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    if ( show(12) ) then
      ido = 0
      call linclr(icolor(5))
    else
      ido = 1
    end if

    aname = object(5)
    lbar2 = lbar .and. .not.show(12)

    call lincon(p,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
      maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

  end if
!
!  S: Draw stream lines in melt region.
!
  if ( show(6) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    if ( show(35) ) then
      ido = 0
      call linclr(icolor(6))
    else
      ido = 1
    end if
!
!  Temporary
!
    ido = 0

    aname = object(6)
    lbar2 = lbar

    call lincon(psi,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
      maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

  end if
!
!  TE: Draw turbulent epsilon contours in melt region.
!
  if ( show(22) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    ido = 1

    aname = object(22)
    lbar2 = lbar

    call lincon(te,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
      maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

  end if
!
!  TEMP: Draw temperature contours in melt and crystal regions.
!
  if ( show(11) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. (kmelt(i,j) >= 2 .or. kcrys(i,j) >= 2) ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    if ( show(13) ) then
      ido = 0
      call linclr(icolor(11))
    else
      ido = 1
    end if

    aname = object(11)
    lbar2 = lbar .and. .not.show(13)

    call lincon(t,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
      maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

  end if
!
!  TK: Draw turbulent KE contours in melt region.
!
  if ( show(23) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    ido = 1

    aname = object(23)
    lbar2 = lbar

    call lincon(tk,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
      maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

  end if
!
!  VIS: Draw viscosity contours in melt region.
!
  if ( show(24) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    ido = 1

    aname = object(24)
    lbar2 = lbar

    call lincon(gamt,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
      maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

  end if
!
!  VMAG: Draw velocity magnitude contours in melt region.
!
  if ( show(36) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    if ( show(37) ) then
      ido = 0
      call linclr(icolor(36))
    else
      ido = 1
    end if

    aname = object(36)
    lbar2 = lbar .and. .not.show(37)

    call lincon(vmag,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
      maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

  end if
!
!  VORT: Draw vorticity contours in melt.
!
  if ( show(41) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    if ( show(42) ) then
      ido = 0
      call linclr(icolor(41))
    else
      ido = 1
    end if

    aname = object(41)
    lbar2 = lbar .and. .not.show(42)

    call lincon(vort,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
  maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

  end if
!
!  W: Draw angular momentum contours in melt region.
!
  if ( show(10) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) .and. kmelt(i,j) >= 2 ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    if ( show(14) ) then
      ido = 0
      call linclr(icolor(10))
    else
      ido = 1
    end if

    aname = object(10)
    lbar2 = lbar .and. .not.show(14)

    call lincon(w,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
      maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

  end if

!
!  XC: Draw X coordinate contours.
!
  if ( show(30) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    if ( show(31) ) then
      ido = 0
      call linclr(icolor(30))
    else
      ido = 1
    end if

    aname = object(30)
    lbar2 = lbar .and. .not.show(31)

    call lincon(xp,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
      maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

  end if
!
!  YC: Draw Y coordinate contours.
!
  if ( show(32) ) then

    do i = 1,l
      do j = 1,m
        if ( lflag(i,j) ) then
          nflag(i,j) = .true.
        else
          nflag(i,j) = .false.
        end if
      end do
    end do

    if ( show(33) ) then
      ido = 0
      call linclr(icolor(32))
    else
      ido = 1
    end if

    aname = object(32)
    lbar2 = lbar .and. .not.show(33)

    call lincon(yp,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
      maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

  end if


!
!  V: Draw velocity vectors in melt region.
!
  if ( show(8) ) then

    call linclr(icolor(8))

    do i = 1,l
      do j = 1,m

        nflag(i,j) = npflag(i,j)

        if ( mod(i-1,nxskip) /= 0 ) then
          nflag(i,j) = .false.
        end if

        if ( mod(j-1,nyskip) /= 0 ) then
          nflag(i,j) = .false.
        end if

        if ( .not.vpflag(i,j) ) then
          nflag(i,j) = .false.
        end if

        if ( kmelt(i,j) <= 1 ) then
          nflag(i,j) = .false.
        end if

        if ( .not.lflag(i,j) ) then
          nflag(i,j) = .false.
        end if

      end do
    end do
!
!  Get the maximum velocity norm of the displayed data.
!
    call vsize(l,m,maxl,maxm,nflag,u,v,vmax)

    if ( vmax == 0.0 ) then

      write ( *, * ) ' '
      write ( *, * ) 'GRAPH - Warning!'
      write ( *, * ) '  No nonzero velocity vectors were visible.'

    else

      velscl = scalev*0.5*min(nxskip*delx,nyskip*dely)/vmax

      call vector(l,m,maxl,maxm,nflag,velscl,u,v,xp,yp)

    end if

  end if
!
!  U: Draw velocity direction vectors in melt region.
!
  if ( show(9) ) then

    call linclr(icolor(9))

    do i = 1,l
      do j = 1,m

        nflag(i,j) = npflag(i,j)

        if ( mod(i-1,nxskip) /= 0 ) then
          nflag(i,j) = .false.
        end if

        if ( mod(j-1,nyskip) /= 0 ) then
          nflag(i,j) = .false.
        end if

        if ( .not.vpflag(i,j) ) then
          nflag(i,j) = .false.
        end if

        if ( kmelt(i,j) <= 1 ) then
          nflag(i,j) = .false.
        end if

        if ( .not.lflag(i,j) ) then
          nflag(i,j) = .false.
        end if

      end do
    end do

    nonzer = 0

    do i = 1,l
      do j = 1,m

        uu = u(i,j)
        vv = v(i,j)
        uvmag = sqrt(uu*uu+vv*vv)

        if ( uvmag /= 0.0 ) then
          nonzer = nonzer+1
          s1(i,j) = uu/uvmag
          s2(i,j) = vv/uvmag
        else
          s1(i,j) = 0.0
          s2(i,j) = 0.0
        end if

      end do
    end do

    if ( nonzer > 0 ) then

      velscl = scalev*0.5*min(nyskip*delx,nxskip*dely)

      call vector(l,m,maxl,maxm,nflag,velscl,s1,s2,xp,yp)

    else

      write ( *, * ) ' '
      write ( *, * ) 'GRAPH - Warning!'
      write ( *, * ) '  No nonzero velocity vectors were visible.'

    end if

  end if
!
!  Draw (dXCdP, dYCdP) direction vectors in entire region.
!
  if ( show(27) ) then

    call linclr(icolor(27))

    do i = 1,l
      do j = 1,m

        nflag(i,j) = ncflag(i,j)

        if ( mod(i-1,nxskip) /= 0 ) then
          nflag(i,j) = .false.
        end if

        if ( mod(j-1,nyskip) /= 0 ) then
          nflag(i,j) = .false.
        end if

        if ( .not.lflag(i,j) ) then
          nflag(i,j) = .false.
        end if

      end do
    end do

    nonzer = 0

    do i = 1,l
      do j = 1,m

        uu = dxcdp(i,j)
        vv = dycdp(i,j)
        uvmag = sqrt(uu*uu+vv*vv)

        if ( uvmag /= 0.0 ) then
          nonzer = nonzer+1
          s1(i,j) = uu/uvmag
          s2(i,j) = vv/uvmag
        else
          s1(i,j) = 0.0
          s2(i,j) = 0.0
        end if

      end do
    end do

    if ( nonzer > 0 ) then

      velscl = scalev*0.5*min(nyskip*delx,nxskip*dely)

      call vector(l,m,maxl,maxm,nflag,velscl,s1,s2,xp,yp)

    else

      write ( *, * ) ' '
      write ( *, * ) 'GRAPH - Warning!'
      write ( *, * ) '  No (dXCdP,dYCdP) directions, all zero.'

    end if

  end if
!
!  Pause, if we are doing X-Windows.
!
  call buzz(dev,x1min,x1max,y1min,y1max)

  return
end
subroutine half(e,gamt,jcrys1,jcrys2,kcrys,kmelt,kvoid,l,m, &
  maxbot,maxl,maxm,nbot,ncflag,npflag,p,pc,psi,rueta,ruksi, &
  t,te,tk,u,v,vmag,vort,vpflag,w,xbot,xc,xp,ybot,yc,yp)

!*****************************************************************************80
!
!! HALF removes the reflected data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real E(MAXL,MAXM), the magnetic stream function.
!
!    Input/output, real GAMT(MAXL,MAXM), the diffusion coefficient.
!
!    Input/output, integer ( kind = 4 ) JCRYS1, the J coordinate of the column of 
!    corner nodes that define the left side of the crystal.
!
!    Input/output, integer ( kind = 4 ) JCRYS2, the J coordinate of the column of 
!    corner nodes that define the right side of the crystal.
!
!    Input/output, integer ( kind = 4 ) KCRYS(MAXL,MAXM), 
!    0, if control volume (I,J) is away from the crystal.
!    1, if control volume (I,J) is on the external crystal boundary.
!    2, if control volume (I,J) is on the internal crystal boundary.
!    3, if control volume (I,J) is in the crystal interior.
!
!    Input/output, integer ( kind = 4 ) KMELT(MAXL,MAXM),
!    0, if control volume (I,J) is away from the melt.
!    1, if control volume (I,J) is on the external melt boundary.
!    2, if control volume (I,J) is on the internal melt boundary.
!    3, if control volume (I,J) is in the melt interior.
! 
!    Input/output, integer ( kind = 4 ) KVOID(MAXL,MAXM), 
!    0, if control volume (I,J) is away from the void.
!    1, if control volume (I,J) is on the external void boundary.
!    2, if control volume (I,J) is on the internal void boundary.
!    3, if control volume (I,J) is in the void interior.
! 
!    Input, integer ( kind = 4 ) L, the number of rows of data.
!
!    Input/output, integer ( kind = 4 ) M, the number of columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
!    Input/output, real P(MAXL,MAXM), the pressure.
!
!    Input/output, real PC(MAXL,MAXM), the corrected pressure.
!
!    Input/output, real PSI(MAXL,MAXM), the stream function.
!
!    Input/output, real RUETA(MAXL,MAXM), RUKSI(MAXL,MAXM), the
!    the momentum in the ETA and KSI directions.
!
!    Input/output, real T(MAXL,MAXM), the temperature.
!
!    Input/output, real TE(MAXL,MAXM), the turbulent epsilon.
!
!    Input/output, real TK(MAXL,MAXM), the turbulent K.
!
!    Input/output, real U(MAXL,MAXM), the horizontal velocity.
!
!    Input/output, real V(MAXL,MAXM), the vertical velocity.
!
!    Input/output, real W(MAXL,MAXM), the axial velocity.
!
!    Input/output, real X(MAXL,MAXM), the X coordinates of the primary
!    nodes.
!
!    Input/output, real XC(MAXL,MAXM), the X coordinates of the corner
!    nodes.
!
!    Input/output, real Y(MAXL,MAXM), the Y coordinates of the primary
!    nodes.
!
!    Input/output, real YC(MAXL,MAXM), the Y coordinates of the corner
!    nodes.
!
  implicit none

  integer ( kind = 4 ) maxbot
  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  real e(maxl,maxm)
  real gamt(maxl,maxm)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys1
  integer ( kind = 4 ) jcrys2
  integer ( kind = 4 ) l
  integer ( kind = 4 ) kcrys(maxl,maxm)
  integer ( kind = 4 ) kmelt(maxl,maxm)
  integer ( kind = 4 ) kvoid(maxl,maxm)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nbot
  logical ncflag(maxl,maxm)
  logical npflag(maxl,maxm)
  real p(maxl,maxm)
  real pc(maxl,maxm)
  real psi(maxl,maxm)
  real rueta(maxl,maxm)
  real ruksi(maxl,maxm)
  real t(maxl,maxm)
  real te(maxl,maxm)
  real tk(maxl,maxm)
  real u(maxl,maxm)
  real v(maxl,maxm)
  real vmag(maxl,maxm)
  real vort(maxl,maxm)
  logical vpflag(maxl,maxm)
  real w(maxl,maxm)
  real xbot(maxbot)
  real xc(maxl,maxm)
  real xp(maxl,maxm)
  real ybot(maxbot)
  real yc(maxl,maxm)
  real yp(maxl,maxm)

  m = (m+1)/2

  call rhalf(e,l,m,maxl,maxm)
  call rhalf(gamt,l,m,maxl,maxm)
  call ihalf(kcrys,l,m,maxl,maxm)
  call ihalf(kmelt,l,m,maxl,maxm)
  call ihalf(kvoid,l,m,maxl,maxm)
  call lhalf(ncflag,l,m,maxl,maxm)
  call lhalf(npflag,l,m,maxl,maxm)
  call rhalf(p,l,m,maxl,maxm)
  call rhalf(pc,l,m,maxl,maxm)
  call rhalf(psi,l,m,maxl,maxm)
  call rhalf(rueta,l,m,maxl,maxm)
  call rhalf(ruksi,l,m,maxl,maxm)
  call rhalf(t,l,m,maxl,maxm)
  call rhalf(te,l,m,maxl,maxm)
  call rhalf(tk,l,m,maxl,maxm)
  call rhalf(v,l,m,maxl,maxm)
  call rhalf(u,l,m,maxl,maxm)
  call rhalf(vort,l,m,maxl,maxm)
  call lhalf(vpflag,l,m,maxl,maxm)
  call rhalf(w,l,m,maxl,maxm)
  call rhalf(xp,l,m,maxl,maxm)
  call rhalf(xc,l,m,maxl,maxm)
  call rhalf(yp,l,m,maxl,maxm)
  call rhalf(yc,l,m,maxl,maxm)
!
!  Take care of data defining the extent of the crystal.
!
  jcrys1 = 1
  jcrys2 = jcrys2-m+1
!
!  Take care of data defining the crucible shape.
!
  nbot = (nbot+1)/2

  do i = 1,nbot
    xbot(i) = xbot(i+nbot-1)
    ybot(i) = ybot(i+nbot-1)
  end do
!
!  Set the velocity magnitudes.
!
  do i = 1,l
    do j = 1,m
      vmag(i,j) = sqrt(u(i,j)**2+v(i,j)**2)
    end do
  end do

  return
end
subroutine hello

!*****************************************************************************80
!
!! HELLO prints out the program name, date, and purpose.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
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

  write ( *, * ) ' '
  write ( *, * ) 'HELLO!'
  write ( *, * ) ' '
  write ( *, * ) '  CRYSTAL_PLOT'
  write ( *, * ) '  Last revised on 12 June 1996'
  write ( *, * ) ' '
  write ( *, * ) '  CRYSTAL_PLOT produces X Window, PostScript or CGM'
  write ( *, * ) '  graphics output from the data of the 2D finite '
  write ( *, * ) '  volume program MASTRAPP2D.'
  write ( *, * ) ' '
  write ( *, * ) '  Scalar quantities can be displayed using contour'
  write ( *, * ) '  lines or contour colors.'
  write ( *, * ) ' '
  write ( *, * ) '  Vector quantities can be displayed as vectors,'
  write ( *, * ) '  vector direction fields, stream lines or stream'
  write ( *, * ) '  contours, and magnitude contours.'
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) '  CHECK OUT THE NEW CRYSC command!'
  write ( *, * ) ' '

  return
end
subroutine help

!*****************************************************************************80
!
!! HELP prints out a list of commands.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
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

  write ( *, * ) ' '
  write ( *, * ) 'Commands:'
  write ( *, * ) ' '
  write ( *, * ) 'AXis     Show the axis of symmetry.'
  write ( *, * ) 'BO       Show the boundary.'
  write ( *, * ) 'BAr      Show color bar.'
  write ( *, * ) 'C        Choose colors.'
  write ( *, * ) 'CC =       Choose color table.'
  write ( *, * ) 'CP       Corrected pressure contours.'
  write ( *, * ) 'CPC      Corrected pressure colors.'
  write ( *, * ) 'CRUcible Show the crucible wall.'
  write ( *, * ) 'CRYB     Show the crystal boundary.'
  write ( *, * ) 'CRYS     Display items in crystal.'
  write ( *, * ) 'CRYSC    Color in the crystal.'
  write ( *, * ) 'CV       Show control volumes.'
  write ( *, * ) 'CVN      Show control volume numbers.'
  write ( *, * ) 'DAT =      Specify the plot data file.'
  write ( *, * ) 'DEV =      Specify the graphics output device.'
  write ( *, * ) 'DIFF     Show the difference between two files.'
  write ( *, * ) 'DOUBLE   Reflect the region about the Y axis.'
  write ( *, * ) 'DXCDP    Plot (dXCdP,dYCdP) for a parameter.'
  write ( *, * ) 'ECHO     Echo user input to output file.'
  write ( *, * ) 'FILE =     Set the output file name.'
  write ( *, * ) 'FRAME    Show the frame.'
  write ( *, * ) 'FS       Show the free surface.'
  write ( *, * ) 'G        Show the current graph.'
  write ( *, * ) 'GRACE =    Set the grace margin.'
  write ( *, * ) 'H        Help (print this list).'
  write ( *, * ) 'HALF     Don''t reflect the region about the Y axis.'
  write ( *, * ) 'ICMAX    Set maximum available color.'
  write ( *, * ) 'ICMIN    Set minimum available color.'
  write ( *, * ) 'JCMAX    Set maximum used color.'
  write ( *, * ) 'JCMIN    Set minimum used color.'
  write ( *, * ) 'MELT     Display items in the melt.'
  write ( *, * ) 'MELTC    Color in the melt (NOT IMPLEMENTED YET).'
  write ( *, * ) 'NB       Show bottom nodes.'
  write ( *, * ) 'NC       Show corner nodes.'
  write ( *, * ) 'NP       Show primary nodes.'
  write ( *, * ) 'NCON =     Set number of contour lines.'
  write ( *, * ) 'NOFRAME  Do not show frame.'
  write ( *, * ) 'NXSKIP   Set X node skip.'
  write ( *, * ) 'NYSKIP   Set Y node skip.'
  write ( *, * ) 'OV       Overlay next graphs.'
  write ( *, * ) 'P        Pressure contours.'
  write ( *, * ) 'PC       Pressure colors.'
  write ( *, * ) 'Q        Quit.'
  write ( *, * ) 'S        Stream lines.'
  write ( *, * ) 'SC       Stream colors.'
  write ( *, * ) 'SCALECV =  Set scale factor for control volume numbers.'
  write ( *, * ) 'SCALENB =  Set scale factor for bottom nodes.'
  write ( *, * ) 'SCALENC =  Set scale factor for corner nodes.'
  write ( *, * ) 'SCALENP =  Set scale factor for primary nodes.'
  write ( *, * ) 'SCALEV =   Set velocity vector scale.'
  write ( *, * ) 'TE       Turbulent epsilon.'
  write ( *, * ) 'TEMP     Temperatures.'
  write ( *, * ) 'TEMPC    Temperature colors.'
  write ( *, * ) 'TITLE =    Set title.'
  write ( *, * ) 'TITLE2 =   Set subtitle.'
  write ( *, * ) 'TK       Turbulent K.'
  write ( *, * ) 'U        Unit velocities.'
  write ( *, * ) 'V        Velocity vectors.'
  write ( *, * ) 'VIS      Viscosity.'
  write ( *, * ) 'VM       Velocity magnitude contours.'
  write ( *, * ) 'VMC      Velocity magnitude colors.'
  write ( *, * ) 'VOID     Display items in the void.'
  write ( *, * ) 'VOIDC    Color in the void (NOT IMPLEMENTED YET).'
  write ( *, * ) 'VORT     Vorticity contours.'
  write ( *, * ) 'VORTC    Vorticity colors.'
  write ( *, * ) 'VPN      Set visibile primary nodes.'
  write ( *, * ) 'W        Angular momentum contours.'
  write ( *, * ) 'WC       Angular momentum colors.'
  write ( *, * ) 'XC       X coordinate contours.'
  write ( *, * ) 'XCC      X coordinate colors.'
  write ( *, * ) 'YC       Y coordinate contours.'
  write ( *, * ) 'YCC      Y coordinate colors.'
  write ( *, * ) ' '
  write ( *, * ) 'TH MH BH  top, middle, bottom halves.'
  write ( *, * ) 'LH CH RH  left, center, right halves.'
  write ( *, * ) ' '
  write ( *, * ) 'TL TC TR  top left, center, right quarters.'
  write ( *, * ) 'ML MC MR  middle left, center, right quarters.'
  write ( *, * ) 'BL BC BR  bottom left, center, right quarters.'
  write ( *, * ) 'X         choose arbitrary X, Y subwindow.'
  write ( *, * ) ' '
  write ( *, * ) 'FULL      return to full picture.'

  return
end
subroutine iduble(ia,l,m,maxl,maxm)

!*****************************************************************************80
!
!! IDUBLE shifts and copies integer data when the region is "doubled".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) IA(MAXL,MAXM), the array to
!    be adjusted.
!
!    Input, integer ( kind = 4 ) L, M, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  integer ( kind = 4 ) ia(maxl,maxm)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
!
!  Shift data to the right M-1 columns.
!
  do i = 1,l
    do jj = 2*m-1,m,-1
      ia(i,jj) = ia(i,jj+1-m)
    end do
  end do
!
!  Now "reflect" data into columns 1 to M-1.
!
  do i = 1,l
    do j = 1,m-1
      ia(i,j) = ia(i,2*m-j)
    end do
  end do

  return
end
subroutine ihalf(ia,l,m,maxl,maxm)

!*****************************************************************************80
!
!! IHALF adjusts an integer array when the region is cut in half.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) IA(MAXL,MAXM), the array to be shifted.
!
!    Input, integer ( kind = 4 ) L, M, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  integer ( kind = 4 ) ia(maxl,maxm)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
!
!  Shift data to the left M-1 columns.
!
  do i = 1,l
    do j = 1,m
      ia(i,j) = ia(i,j+m-1)
    end do
  end do

  return
end
subroutine init(b1jbl,b2jbl,b3jbl,cost,cost2,delx,dely,dev,echo,fildat, &
  filgrf,filinp,grace,icmax,icmin,icolor,icrys1,icrys2,iplot,itable,jcmax, &
  jcmin,jcrys1,jcrys2,l,lbar,m,maxbot,maxl,maxm,maxobj,nbot,ncflag,ncon, &
  nflag,npflag,nxskip,nyskip,object,ovrlay,p,pc,psi,reflec,scalecv,scalenb, &
  scalenc,scalenp,scalev,show,title,title2,u,v,vpflag,vmag,vort,x1max,x1min, &
  x2max,x2min,xbot,xc,xmax,xmin,xp,xsmax,xsmin,xtmax,xtmin,y1max,y1min,y2max, &
  y2min,ybot,yc,ymax,ymin,yp,ysmax,ysmin,ytmax,ytmin)

!*****************************************************************************80
!
!! INIT initializes the values of data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  DEV    Output, character*10 DEV.
!         The graphics output device to be used.  Current legal
!         values for DEV include:
! 
!         cgmb - CGM binary file.
!         ps   - PostScript file.
!         xws  - X window screen (interactive).
!
!  FILDAT Output, character*80 FILDAT.
!         The name of the data file to be read in, which contains
!         the information defining the mesh and the physical
!         parameters.
! 
!  FILGRF Output, character*80 FILGRF, the name of the output
!         graphics file.
!
!  FILINP Output, character*80 FILINP.
!         The name of a file in which will be placed a copy of the
!         input typed by the user while running Display.
!
!  GRACE  Output, real GRACE.
!         The size of the "grace" margin on the plot.
! 
!  ICMAX,
!  ICMIN  Output, integer ( kind = 4 ) ICMAX, ICMIN, the maximum and 
!         minimum color indices to use in the color bar.
!
!  ICOLOR Output, integer ( kind = 4 ) ICOLOR(MAXOBJ).
!         Contains the color indexes for each object.
!         However, in some cases, ICOLOR is actual a color table
!         index.
! 
!  ICRYS1 Output, integer ( kind = 4 ) ICRYS1, the I coordinate of the row of corner
!         nodes that define the bottom of the crystal.
!
!  ICRYS2 Output, integer ( kind = 4 ) ICRYS2, the I coordinate of the row of corner
!         nodes that define the top of the crystal.
!
!  IPLOT  Output, integer ( kind = 4 ) IPLOT.
!         The number of plots made so far.
! 
!  ITABLE Output, integer ( kind = 4 ) ITABLE, the desired color table.
!  
!         1: low black to high white
!         2: low blue to high yellow
!         3: low red, high blue, with bands between.
!         4: low red, yellow, green, blue, high white.
!         5: low white, blue, green, yellow, high red.
!         6: low blue to high red
!         7: linear table between 2 user colors.
!         8: linear table between N user colors.
!         9: low white to high black.
! 
!  JCMAX,
!  JCMIN  Output, integer ( kind = 4 ) JCMAX, JCMIN, the maximum and 
!         minimum color indices to use in the color bar.
!
!  LBAR   Output, logical LBAR, is .TRUE. if the color bar should
!         be shown.
!
!  MAXBOT Input, integer ( kind = 4 ) MAXBOT, the maximum storage for the nodes
!         that define the crucible bottom.
!
!  MAXL,
!  MAXM   Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!         rows and columns of data.
!
!  MAXOBJ Input, integer ( kind = 4 ) MAXOBJ.
!         The number of graphical "objects".
! 
!  NBOT   Output, integer ( kind = 4 ) NBOT, the number of crucible nodes.
!
!  NCON   Output, integer ( kind = 4 ) NCON.
!         The number of contour lines to be drawn.  This is
!         initialized to 12, but may be changed by the user.
! 
!  NFLAG  Output, logical NFLAG(MAXL,MAXM).
!
!         NFLAG is used to "flag" which nodes are active,
!         that is, to be displayed, and which not, in the graph.
!
!  NXSKIP Output, integer ( kind = 4 ) NXSKIP.
!         NXSKIP is used to "thin" out a vector plot.  
!
!         If NXSKIP = 1, then a standard vector plot is made.
!  
!         Otherwise, in the X direction, vectors are drawn only
!         in columns 1, 1+NXSKIP, 1+2*NXSKIP and so on.
!
!  NYSKIP Output, integer ( kind = 4 ) NYSKIP.
!         NYSKIP is used to "thin" out a vector plot.  
!
!         If NYSKIP = 1, then a standard vector plot is made.
!
!         Otherwise, in the Y direction, vectors are drawn only
!         in rows 1, 1+NYSKIP, 1+2*NYSKIP and so on.
!  
!  OBJECT Input, character*25 OBJECT(MAXOBJ), the names of the 
!         graphical objects.
!
!  OVRLAY Input, logical OVRLAY.
!         If OVRLAY is true, then the next time that a plot is
!         requested, a "new frame" command is suppressed, so that
!         the new plot is shown on top of the previous one.
! 
!  P      Output, real P(MAXL,MAXM), the pressure.
!
!  PC     Output, real PC(MAXL,MAXM), the corrected pressure.
!
!  PSI    Output, real PSI(MAXL,MAXM), the stream function.
!
!  REFLEC Output, logical REFLEC, .TRUE. if the region is
!         reflected, .FALSE. otherwise.
!
!  SCALECV Output, real SCALECV.
!         A scale factor for the control volume numbers, which has a
!         default value of 1.
!
!  SCALENB Output, real SCALENB.
!         A scale factor for the size of the nodes that define the
!         crucible shape, which has a default value of 1.
!
!  SCALENC Output, real SCALENC.
!         A scale factor for the size of the corner nodes, which has a
!         default value of 1.
!
!  SCALENP Output, real SCALENP.
!         A scale factor for the size of the primary nodes, which has a
!         default value of 1.
!
!  SCALEV Output, real SCALEV.
!         A scale factor for velocity vectors.  This starts out at 1.0.
! 
!  SHOW   Output, logical SHOW(MAXOBJ).
!         Contains, for each object, a flag determining whether it
!         is to be shown or not.
!
!  TITLE  Output, character*40 TITLE.
!         A title for the plots.
! 
!  TITLE2 Output, character*40 TITLE2.
!         A subtitle used in the profile plots.
!
!  U      Output, real U(MAXL,MAXM), the horizontal velocity.
!
!  V      Output, real V(MAXL,MAXM), the vertical velocity.
!
!  VPFLAG Output, logical VPFLAG(MAXL,MAXM), a flag which determines
!         whether vector data at a given node should be displayed.
!
!  VMAG   Output, real VMAG(MAXL,MAXM), the velocity magnitudes.
!
!  X      Output, real X(MAXL,MAXM), the X coordinates of the primary
!         nodes.
!
!  X1MAX,
!  X1MIN  Output, real X1MAX, X1MIN, the maximum and minimum X 
!         coordinates of the plot, which includes a small grace margin.
!
!  X2MAX,
!  X2MIN  Input, real X2MAX, X2MIN, the maximum and minimum X 
!         coordinates that should be used for plotting.  No plotting 
!         commands should exceed these values.  This is where the 
!         "frame" might be drawn.
!
!  XBOT   Output, real XBOT(NBOT), the X coordinates of the crucible
!         bottom nodes.
!
!  XSMAX  Input, real XSMAX.
!         The maximum X coordinate of the data to be displayed.
!         XSMAX defaults to XMAX, but can be made smaller to
!         focus on a portion of the region.
! 
!  XSMIN  Input, real XSMIN.
!         The minimum X coordinate of the data to be displayed.
!         XSMIN defaults to XMIN, but can be made larger to
!         focus on a portion of the region.
! 
!  Y      Output, real Y(MAXL,MAXM), the Y coordinates of the primary
!         nodes.
!
!  Y1MAX,
!  Y1MIN  Output, real Y1MAX, Y1MIN, the maximum and minimum Y 
!         coordinates of the plot, which includes a small grace margin.
!
!  Y2MAX,
!  Y2MIN  Input, real Y2MAX, Y2MIN, the maximum and minimum Y 
!         coordinates that should be used for plotting.  No plotting 
!         commands should exceed these values.  This is where the 
!         "frame" might be drawn.
!
!  XBOT   Output, real XBOT(NBOT), the X coordinates of the crucible
!         bottom nodes.
!
!  YSMAX  Input, real YSMAX.
!         The maximum Y coordinate of the data to be displayed.
!         YSMAX defaults to YMAX, but can be made smaller to
!         focus on a portion of the region.
! 
!  YSMIN  Input, real YSMIN.
!         The minimum Y coordinate of the data to be displayed.
!         YSMIN defaults to YMIN, but can be made larger to
!         focus on a portion of the region.
!
  implicit none

  integer ( kind = 4 ) maxbot
  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm
  integer ( kind = 4 ) maxobj

  real b1jbl(maxm)
  real b2jbl(maxm)
  real b3jbl(maxl)
  real cost
  real cost2
  real delx
  real dely
  character ( len = 10 ) dev
  logical echo
  character ( len = 80 ) fildat
  character ( len = 80 ) filgrf
  character ( len = 80 ) filinp
  real grace
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icmax
  integer ( kind = 4 ) icmin
  integer ( kind = 4 ) icolor(maxobj)
  integer ( kind = 4 ) icrys1
  integer ( kind = 4 ) icrys2
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) itable
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcmax
  integer ( kind = 4 ) jcmin
  integer ( kind = 4 ) jcrys1
  integer ( kind = 4 ) jcrys2
  integer ( kind = 4 ) l
  logical lbar
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nbot
  logical ncflag(maxl,maxm)
  integer ( kind = 4 ) ncon
  logical nflag(maxl,maxm)
  logical npflag(maxl,maxm)
  integer ( kind = 4 ) nxskip
  integer ( kind = 4 ) nyskip
  character ( len = 25 ) object(maxobj)
  logical ovrlay
  real p(maxl,maxm)
  real pc(maxl,maxm)
  real psi(maxl,maxm)
  logical reflec
  real scalecv
  real scalenb
  real scalenc
  real scalenp
  real scalev
  logical show(maxobj)
  character ( len = 40 ) title
  character ( len = 40 ) title2
  real u(maxl,maxm)
  real v(maxl,maxm)
  logical vpflag(maxl,maxm)
  real vmag(maxl,maxm)
  real vort(maxl,maxm)
  real x1max
  real x1min
  real x2max
  real x2min
  real xbot(maxbot)
  real xc(maxl,maxm)
  real xmax
  real xmin
  real xp(maxl,maxm)
  real xsmax
  real xsmin
  real xtmax
  real xtmin
  real y1max
  real y1min
  real y2max
  real y2min
  real ybot(maxbot)
  real yc(maxl,maxm)
  real ymax
  real ymin
  real yp(maxl,maxm)
  real ysmax
  real ysmin
  real ytmax
  real ytmin

  do i = 1,maxm
    b1jbl(i) = 0.0
    b2jbl(i) = 0.0
  end do

  do i = 1,maxl
    b3jbl(i) = 0.0
  end do

  cost = 0.0
  cost2 = 0.0
  delx = 0.0
  dely = 0.0
  dev = ' '
  echo = .false.
  fildat = 'cplot.dat'
  filgrf = ' '
  filinp = 'cplot.inp'
  grace = 0.01
  icmax = 255
  icmin = 2

  do i = 1,maxobj
    icolor(i) = 1
  end do
  icolor(2) = 2

  icrys1 = 1
  icrys2 = 1
  iplot = 0
  itable = 9
  jcmax = 200
  jcmin = 32
  jcrys1 = 1
  jcrys2 = 1
  l = 1
  lbar = .true.
  m = 1
  nbot = 0

  do i = 1,maxl
    do j = 1,maxm
      ncflag(i,j) = .false.
    end do
  end do

  ncon = 9

  do i = 1,maxl
    do j = 1,maxm
      nflag(i,j) = .false.
    end do
  end do

  do i = 1,maxl
    do j = 1,maxm
      npflag(i,j) = .false.
    end do
  end do

  nxskip = 1
  nyskip = 1

  object(1) = 'Boundary'
  object(2) = 'Element'
  object(3) = 'Frame'
  object(4) = 'Primary nodes'
  object(5) = 'Pressure'
  object(6) = 'Stream lines'
  object(7) = 'Title'
  object(8) = 'Velocities'
  object(9) = 'Directions'
  object(10) = 'Ang Mom'
  object(11) = 'Temperature'
  object(12) = 'Pressure colors'
  object(13) = 'Temperature colors'
  object(14) = 'Ang Mom colors'
  object(15) = 'X velocity'
  object(16) = 'Y velocity'
  object(17) = 'X velocity colors'
  object(18) = 'Y velocity colors'
  object(19) = 'Crucible'
  object(20) = 'Crystal side'
  object(21) = 'Free surface'
  object(22) = 'Turbulent Epsilon'
  object(23) = 'Turbulent KE'
  object(24) = 'Viscosity'
  object(25) = 'Corrected pressure'
  object(26) = 'Corrected pressure colors'
  object(27) = '(dXCdP, dYCdP) vectors'
  object(28) = 'Axis of symmetry'
  object(29) = 'Corner nodes'
  object(30) = 'X coordinate'
  object(31) = 'X coordinate colors'
  object(32) = 'Y coordinate'
  object(33) = 'Y coordinate colors'
  object(34) = 'Crucible nodes'
  object(35) = 'Stream colors'
  object(36) = 'Velocity magnitude'
  object(37) = 'Velocity magnitude colors'
  object(38) = 'Void interior'
  object(39) = 'Crystal interior'
  object(40) = 'Melt interior'
  object(41) = 'Vorticity contours'
  object(42) = 'Vorticity colors'
  object(43) = 'control volume numbers'
  object(44) = 'Crystal color'
  object(45) = 'Melt color'
  object(46) = 'Void color'

  ovrlay = .false.

  do i = 1,maxl
    do j = 1,maxm
      p(i,j) = 0.0
    end do
  end do

  do i = 1,maxl
    do j = 1,maxm
      pc(i,j) = 0.0
    end do
  end do

  do i = 1,maxl
    do j = 1,maxm
      psi(i,j) = 0.0
    end do
  end do

  reflec = .false.

  scalecv = 1.0
  scalenb = 1.0
  scalenc = 1.0
  scalenp = 1.0
  scalev = 1.0

  do i = 1,maxobj
    show(i) = .false.
  end do

  show(1) = .false.
  show(3) = .false.
  show(7) = .true.
  show(19) = .true.
  show(20) = .true.
  show(21) = .true.
  show(28) = .true.
  show(38) = .false.
  show(39) = .true.
  show(40) = .true.
 
  title = ' '
  title2 = ' '

  do i = 1,maxl
    do j = 1,maxm
      u(i,j) = 0.0
      v(i,j) = 0.0
    end do
  end do

  do i = 1,maxl
    do j = 1,maxm
      vpflag(i,j) = .true.
    end do
  end do

  do i = 1,maxl
    do j = 1,maxm
      vort(i,j) = 0.0
    end do
  end do

  x2max = 1.0
  x2min = 0.0

  do i = 1,maxbot
    xbot(i) = 0.0
  end do

  do i = 1,maxl
    do j = 1,maxm
      xc(i,j) = 0.0
    end do
  end do

  xmax = 0.0
  xmin = 0.0

  do i = 1,maxl
    do j = 1,maxm
      xp(i,j) = 0.0
    end do
  end do

  xsmax = 0.0
  xsmin = 0.0
  xtmax = 0.0
  xtmin = 0.0

  y2max = 1.0
  y2min = 0.0

  do i = 1,maxbot
    ybot(i) = 0.0
  end do

  do i = 1,maxl
    do j = 1,maxm
      yc(i,j) = 0.0
    end do
  end do

  ymax = 0.0
  ymin = 0.0
  yp(1:maxl,1:maxm) = 0.0
  ysmax = 0.0
  ysmin = 0.0
  ytmax = 0.0
  ytmin = 0.0
!
!  Set things that depend on other things.
!
  do i = 1,maxl
    do j = 1,maxm
      vmag(i,j) = sqrt(u(i,j)**2+v(i,j)**2)
    end do
  end do

  x1min = x2min-grace*(x2max+x2min)
  x1max = x2max+grace*(x2max+x2min)
  y1min = y2min-grace*(y2max+y2min)
  y1max = y2max+(grace+0.05)*(y2max+y2min)

  return
end
function inside(x1,xmid,x2)

!*****************************************************************************80
!
!! INSIDE reports whether XMID is, or is not, between X1 and X2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, XMID, X2, three values to be tested.
!         
!  INSIDE Output, logical INSIDE.
!
!         INSIDE will be TRUE if XMID is "inside" the interval
!         spanned by X1 and X2.  That is, if
!
!           X1 < =  XMID <= X2
!         or
!           X2 < =  XMID <= X1
!
!         Otherwise, INSIDE will be FALSE.
!
  implicit none

  logical inside
  real x1
  real x2
  real xmid

  if ( (x1 <= xmid .and. xmid <= x2) .or. (x1 >= xmid .and. xmid >= x2) ) then
    inside = .true.
  else
    inside = .false.
  end if

  return
end
subroutine lduble(a,l,m,maxl,maxm)

!*****************************************************************************80
!
!! LDUBLE shifts and copies logical data when the region is "doubled".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, logical A(MAXL,MAXM), the array to
!    be adjusted.
!
!    Input, integer ( kind = 4 ) L, M, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  logical a(maxl,maxm)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
!
!  Shift data to right M-1 columns.
!
  do i = 1,l
    do jj = 2*m-1,m,-1
      a(i,jj) = a(i,jj+1-m)
    end do
  end do
!
!  Now "reflect" data into columns 1 to M-1.
!
  do i = 1,l
    do j = 1,m-1
      a(i,j) = a(i,2*m-j)
    end do
  end do

  return
end
subroutine lhalf(a,l,m,maxl,maxm)

!*****************************************************************************80
!
!! LHALF adjusts a logical array when the region is cut in half.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, logical A(MAXL,MAXM), the array to be shifted.
!
!    Input, integer ( kind = 4 ) L, M, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  logical a(maxl,maxm)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
!
!  Shift data to left M-1 columns.
!
  do i = 1,l
    do j = 1,m
      a(i,j) = a(i,j+m-1)
    end do
  end do

  return
end
subroutine lincon(a,aname,echo,icolor,ido,jcmax,jcmin,l,lbar2,m, &
  maxl,maxm,maxobj,ncon,nflag,srange,x2max,x2min,xp,y1max,y2max,yp)

!*****************************************************************************80
!
!! LINCON supervises the creation of a line contour plot.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IDO.
!    0, draw the contours in black.
!    1, draw the contours in a sequence of colors.
!
!    Input, real A(MAXL,MAXM), the quantity to be contoured.
!
!    Input, character*(*) ANAME, the name of the quantity
!    to be contoured.
!
!    Input, integer ( kind = 4 ) ICOLOR(MAXOBJ).
!    Contains the color indexes for each object.
!    However, in some cases, ICOLOR is actual a color table
!    index.
! 
!    Input, integer ( kind = 4 ) JCMAX, JCMIN, the maximum and 
!    minimum color indices to use in the color bar.
!
!    Input, integer ( kind = 4 ) L, the number of rows of data.
!
!    Input, logical LBAR2, is .TRUE. if the color bar should
!    be shown.
!
!    Input, integer ( kind = 4 ) M, the number of columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
!    Input, integer ( kind = 4 ) MAXOBJ.
!    The number of graphical "objects".
! 
!    Input, integer ( kind = 4 ) NCON, the number of color contour
!    regions drawn, and hence, the number of colors
!    to be displayed in the color bar.
!
!    Input, logical NFLAG(MAXL,MAXM).
!    NFLAG is used to "flag" which nodes are active,
!    that is, to be displayed, and which not, in the graph.
!
!    Input, real SRANGE.
!    The maximum of XSMAX-XSMIN and YSMAX-YSMIN.
!    This gives the size of a square containing the data
!    window.
!
!    Input, real X(MAXL,MAXM), the X coordinates of the primary
!    nodes.
!
!    Input, real X2MAX, X2MIN, the maximum and minimum X 
!    coordinates that should be used for plotting.  No plotting 
!    commands should exceed these values.  This is where the 
!    "frame" might be drawn.
!
!    Input, real Y(MAXL,MAXM), the Y coordinates of the primary
!    nodes.
!
!    Input, real Y2MAX, Y2MIN, the maximum and minimum Y 
!    coordinates that should be used for plotting.  No plotting 
!    commands should exceed these values.  This is where the 
!    "frame" might be drawn.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm
  integer ( kind = 4 ) maxobj

  real a(maxl,maxm)
  character ( len = * )  aname
  logical echo
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icolor(maxobj)
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) jcmax
  integer ( kind = 4 ) jcmin
  integer ( kind = 4 ) jcolor
  integer ( kind = 4 ) l
  logical lbar2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ncon
  logical nflag(maxl,maxm)
  real smax
  real smin
  real srange
  real sval

  real x1
  real x2
  real x2max
  real x2min
  real xp(maxl,maxm)
  real y1
  real y1max
  real y2
  real y2max
  real yp(maxl,maxm)
!
!  Get the minimum and maximum values of the data.
!
  call fsize(l,m,maxl,maxm,nflag,a,smax,smin)
!
!  Don't draw the contours if the maximum is not greater than the minimum.
!
  if ( smax <= smin ) then
    write ( *, * ) ' '
    write ( *, * ) 'LINCON - Warning!'
    write(*,'(a)')aname
    write ( *, * ) '  No contours, all values equal.'
    return
  end if
!
!  Allow the user to adjust the range of the data.
!
  call setsiz(echo,aname,smax,smin)
!
!  Draw each contour line.
!
  do i = 1,ncon

    if ( ido == 1 ) then
      jcolor = int(((ncon-i)*jcmin+i*jcmax)/real(ncon))
      call linclr(jcolor)
    end if

    sval = ((ncon+1-i)*smin+i*smax)/real(ncon+1)

    call tricon(l,m,maxl,maxm,nflag,a,sval,xp,yp)

  end do

  if ( lbar2 ) then
    x1 = x2min
    x2 = x2max

    y1 = y2max
    y2 = y1max

    call cbar(icolor,jcmax,jcmin,maxobj,ncon,smax,smin,srange, x1,x2,y1,y2)

  end if

  return
end
function lnei ( s1, s2 )

!*****************************************************************************80
!
!! LNEI compares two strings for non-equality, ignoring case.
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
!    Output, logical LNEI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical lnei
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  lnei = .true.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc+1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc+1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  lnei = .false.

  return
end
subroutine pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax, &
  xsmin,y1max,y1min,y2max,y2min,ysmax,ysmin)

!*****************************************************************************80
!
!! PLTBOX computes a square box containing the data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  GRACE  Input, real GRACE.
!         The size of the "grace" margin on the plot.
! 
!  SRANGE Output, real SRANGE.
!         The maximum of XSMAX-XSMIN and YSMAX-YSMIN.
!         This gives the size of a square containing the data
!         window.
!
!  X1MAX,
!  X1MIN  Output, real X1MAX, X1MIN, the maximum and minimum X 
!         coordinates of the plot, which includes a small grace margin.
!
!  X2MAX,
!  X2MIN  Output, real X2MAX, X2MIN, the maximum and minimum X 
!         coordinates that should be used for plotting.  No plotting 
!         commands should exceed these values.  This is where the 
!         "frame" might be drawn.
!
!  XSMAX  Input, real XSMAX.
!         The maximum X coordinate of the data to be displayed.
!         XSMAX defaults to XMAX, but can be made smaller to
!         focus on a portion of the region.
! 
!  XSMIN  Input, real XSMIN.
!         The minimum X coordinate of the data to be displayed.
!         XSMIN defaults to XMIN, but can be made larger to
!         focus on a portion of the region.
! 
!  Y1MAX,
!  Y1MIN  Output, real Y1MAX, Y1MIN, the maximum and minimum Y 
!         coordinates of the plot, which includes a small grace margin.
!
!  Y2MAX,
!  Y2MIN  Output, real Y2MAX, Y2MIN, the maximum and minimum Y coordinates that
!         should be used for plotting.  No plotting commands should 
!         exceed these values.  This is where the "frame" might be drawn.
!  
!  YSMAX  Input, real YSMAX.
!         The maximum Y coordinate of the data to be displayed.
!         YSMAX defaults to YMAX, but can be made smaller to
!         focus on a portion of the region.
! 
!  YSMIN  Input, real YSMIN.
!         The minimum Y coordinate of the data to be displayed.
!         YSMIN defaults to YMIN, but can be made larger to
!         focus on a portion of the region.
!
  implicit none

  real grace
  real srange
  real x1max
  real x1min
  real x2max
  real x2min
  real xsmax
  real xsmin
  real y1max
  real y1min
  real y2max
  real y2min
  real ysmax
  real ysmin

  srange = max(xsmax-xsmin,ysmax-ysmin)
 
  if ( srange <= 0.0 ) then
    srange = 1.0
  end if
 
  x2min = 0.5*(xsmin+xsmax)-0.5*srange
  x2max = 0.5*(xsmin+xsmax)+0.5*srange
  y2min = 0.5*(ysmin+ysmax)-0.5*srange
  y2max = 0.5*(ysmin+ysmax)+0.5*srange

  x1min = x2min-grace*(x2max-x2min)
  x1max = x2max+grace*(x2max-x2min)
  y1min = y2min-grace*(y2max-y2min)
  y1max = y2max+(grace+0.05)*(y2max-y2min)

  write ( *, * ) ' '
  write ( *, * ) 'PLTBOX - Note:'
  write ( *, * ) ' '
  write ( *, * ) '  New total picture coordinates:'
  write ( *, * ) ' '
  write ( *, * ) x1min,' X1MIN < =  X <= X1MAX ',x1max
  write ( *, * ) y1min,' Y1MIN < =  Y <= Y1MAX ',y1max
  write ( *, * ) ' '
  write ( *, * ) '  New graphing area coordinates:'
  write ( *, * ) ' '
  write ( *, * ) x2min,' X2MIN < =  X <= X2MAX ',x2max
  write ( *, * ) y2min,' Y2MIN < =  Y <= Y2MAX ',y2max
  write ( *, * ) ' '
  write ( *, * ) '  New data window coordinates:'
  write ( *, * ) ' '
  write ( *, * ) xsmin,' XSMIN < =  X <= XSMAX ',xsmax
  write ( *, * ) ysmin,' YSMIN < =  Y <= YSMAX ',ysmax

  return
end
subroutine preplt ( dev, echo, filgrf, icmax, icmin, iplot, itable, &
  lgopen, ovrlay )

!*****************************************************************************80
!
!! PREPLT should be called before doing each plot.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  DEV    Input, character*10 DEV.
!         The graphics output device to be used.  Current legal
!         values for DEV include:
! 
!         cgmb - CGM binary file.
!         ps   - PostScript file.
!         xws  - X window screen (interactive).
! 
!  FILGRF Input, character*80 FILGRF, the name of the output
!         graphics file.
!
!  ICMAX,
!  ICMIN  Input, integer ( kind = 4 ) ICMAX, ICMIN.
!         The maximum and minimum color indices to be used from
!         the color table.  The color table contains 256 colors,
!         but color indices 1 and 2 are black and white, and for some
!         reason, the predefined DRAWCGM tables generally only
!         use 2-200 for sensible colors.
! 
!         Of course the entries in the color table are "off by one".
!         The first entry is for color 0, and the 256-th entry for
!         color 255.
!
!  IPLOT  Input/output, integer ( kind = 4 ) IPLOT.
!         The number of plots made so far.
!
!         If IPLOT is less than or equal to 1, it is reset to 1.
!
!  ITABLE Input, integer ( kind = 4 ) ITABLE, the desired color table.
!  
!         1: low black to high white
!         2: low blue to high yellow
!         3: low red, high blue, with bands between.
!         4: low red, yellow, green, blue, high white.
!         5: low white, blue, green, yellow, high red.
!         6: low blue to high red
!         7: linear table between 2 user colors.
!         8: linear table between N user colors.
!         9: low white to high black.
! 
!  OVRLAY Input, logical OVRLAY.
!         If OVRLAY is true, then the next time that a plot is
!         requested, a "new frame" command is suppressed, so that
!         the new plot is shown on top of the previous one.
!
  implicit none

  integer ( kind = 4 ), parameter :: nval = 2

  character ( len = * )  dev
  logical echo
  character ( len = 80 ) filgrf
  integer ( kind = 4 ) icmax
  integer ( kind = 4 ) icmin
  integer ( kind = 4 ) icolor
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) itable
  logical, save :: linit = .false.
  logical lgopen
  logical ovrlay
  real xval(nval)
  real yval(nval)
!
!  If it's the first picture, then
!
!    Choose an output device,
!    Give the output file a name,
!    Initialize the graphics package.
!
  if ( .not.linit ) then

    linit = .true.

    if ( dev == ' ' ) then
      write ( *, * ) ' '
      write ( *, * ) 'PrePlt - Warning!'
      write ( *, * ) '  No output device was specified.'
      write ( *, * ) '  PostScript output will be generated.'
      dev = 'ps'
    end if

    call device(dev)

    if ( dev == 'cgmb' ) then

      if ( filgrf == ' ' ) then
        filgrf = 'cplot.cgm'
      end if

      call outfil(filgrf)

    else if ( dev == 'ps' ) then

      if ( filgrf == ' ' ) then
        filgrf = 'cplot.ps'
      end if

      call outfil(filgrf)

    end if

    if ( itable == 0 ) then
      itable = 1
    end if

    call settab(echo,icmax,icmin,itable)

    call grfini

    lgopen = .true.

    xval(1) = 0.0
    xval(2) = 1.0
    yval(1) = 0.0
    yval(2) = 1.0
    call setscl(xval,yval,nval)

    icolor = 1
    call linclr(icolor)
    call filclr(icolor)
!
!  Else, if it's a later picture,
!
!    Issue a "new frame" command, unless an overlay was requested.
!
  else

    if ( .not.ovrlay ) then
      call newfrm
    end if

  end if

  iplot = iplot+1

  return
end
subroutine rduble(a,l,m,maxl,maxm,name)

!*****************************************************************************80
!
!! RDUBLE shifts and copies real data when the region is doubled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real A(MAXL,MAXM), the array of data to be
!    modified.
!
!    Input, integer ( kind = 4 ) L, the number of rows of data.
!
!    Input, integer ( kind = 4 ) M, the number of columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
!    Input, character*(*) NAME, the name of the array being
!    shifted.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  real a(maxl,maxm)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  character ( len = * )  name
!
!  Shift data to the right M-1 columns.
!
  do i = 1,l
    do jj = 2*m-1,m,-1
      a(i,jj) = a(i,jj+1-m)
    end do
  end do
!
!  Now "reflect" data into columns 1 to M-1.
!
  do i = 1,l
    do j = 1,m-1
      a(i,j) = a(i,2*m-j)
    end do
  end do
!
!  Special adjustments for certain arrays.
!  (I'm just absolutely guessing that RUETA is the horizontal-like
!  velocity here).
!
  if ( name == 'rueta' .or. name == 'u' .or. name == 'w' .or. &
       name == 'xp' .or. name == 'xc' ) then
    do i = 1,l
      do j = 1,m-1
        a(i,j) = -a(i,j)
      end do
    end do
  end if
    
  return
end
subroutine rhalf(a,l,m,maxl,maxm)

!*****************************************************************************80
!
!! RHALF shifts the information in a real array when the region is halved.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real A(MAXL,MAXM), the array to be shifted.
!
!    Input, integer ( kind = 4 ) L, the number of rows of data.
!
!    Input, integer ( kind = 4 ) M, the number of columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  real a(maxl,maxm)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
!
!  Shift data to left M-1 columns.
!
  do i = 1,l
    do j = 1,m
      a(i,j) = a(i,j+m-1)
    end do
  end do

  return
end
subroutine rsdiff(cost2,dxcdp,dxdp,dycdp,dydp,e,gamt,l,m,maxl, &
  maxm,nbot,p,pc,psi,rueta,ruksi,t,te,tempa,tk,tnow2, &
  u,v,vmag,vort,w,xc,xp,yc,yp)

!*****************************************************************************80
!
!! RSDIFF reads in restart information and computes differences.  
!
!  Discussion:
!
!    We assume one set of information was previously read in.  RSDIFF
!    reads in a second set, and replaces the two sets of information
!    by their differences.
!
!    It has to interchange X and Y data, because of some strange 
!    convention in CRYSTAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  COST2  Output, real COST2, the value of the cost functional 
!         associated with the second solution.
!
!  DXCDP  Output, real DXCDP(MAXL,MAXM), the difference between
!         the second and first values of XC.
!
!  DXDP   Output, real DXDP(MAXL,MAXM), the difference between
!         the second and first values of X.
!
!  DYCDP  Output, real DYCDP(MAXL,MAXM), the difference between
!         the second and first values of YC.
!
!  DYDP   Output, real DYDP(MAXL,MAXM), the difference between
!         the second and first values of Y.
!
!  E      Input/output, real E(MAXL,MAXM).
!
!         On input, the magnetic stream function.
!         On output, the difference between the second and first values.
!
!  GAMT   Input, real GAMT(MAXL,MAXM).
!
!         On input, the diffusion coefficient.
!         On output, the difference between the second and first values.
!
!  L      Input, integer ( kind = 4 ) L, the number of rows of data.
!
!  M      Input, integer ( kind = 4 ) M, the number of columns of data.
!
!  MAXL,
!  MAXM   Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!         rows and columns of data.
!
!  P      Input/output, real P(MAXL,MAXM).
!
!         On input, the pressure.
!         On output, the difference between the second and first values.
!
!  PC     Input/output, real PC(MAXL,MAXM).
!
!         On input, the corrected pressure.
!         On output, the difference between the second and first values.
!
!  PSI    Input/output, real PSI(MAXL,MAXM).
!
!         On input, the stream function.
!         On output, the difference between the second and first values.
!
!  RUETA,
!  RUKSI  Input/output, real RUETA(MAXL,MAXM), RUKSI(MAXL,MAXM).
!
!         On input, the momentum in the ETA and KSI directions.
!         On output, the difference between the second and first values.
!
!  T      Input/output, real T(MAXL,MAXM).
!
!         On input, the temperature.
!         On output, the difference between the second and first values.
!
!  TE     Input/output, real TE(MAXL,MAXM).
!
!         On input, the turbulent epsilon.
!         On output, the difference between the second and first values.
!
!  TEMP   Workspace, real TEMP(MAXL,MAXM).
!
!  TK     Input/output, real TK(MAXL,MAXM).
!
!         On input, the turbulent K.
!         On output, the difference between the second and first values.
!
!  TNOW2  Output, real TNOW2, the second time.
!
!  U      Input/output, real U(MAXL,MAXM).
!
!         On input, the horizontal velocity.
!         On output, the difference between the second and first values.
!
!  V      Input/output, real V(MAXL,MAXM).
!
!         On input, the vertical velocity.
!         On output, the difference between the second and first values.
!
!  W      Input/output, real W(MAXL,MAXM).
!
!         On input, the axial velocity.
!         On output, the difference between the second and first values.
!
!  X      Input, real X(MAXL,MAXM), the X coordinates of the primary
!         nodes.
!
!  XC     Input, real XC(MAXL,MAXM), the X coordinates of the corner
!         nodes.
!
!  Y      Input, real Y(MAXL,MAXM), the Y coordinates of the primary
!         nodes.
!
!  YC     Input, real YC(MAXL,MAXM), the Y coordinates of the corner
!         nodes.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  real cost2
  real dxcdp(maxl,maxm)
  real dxdp(maxl,maxm)
  real dycdp(maxl,maxm)
  real dydp(maxl,maxm)
  real e(maxl,maxm)
  real gamt(maxl,maxm)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) nbot
  real p(maxl,maxm)
  real pc(maxl,maxm)
  real psi(maxl,maxm)
  real rueta(maxl,maxm)
  real ruksi(maxl,maxm)
  real t(maxl,maxm)
  real te(maxl,maxm)
  real tempa(maxl,maxm)
  real tk(maxl,maxm)
  real tnow2
  real u(maxl,maxm)
  real v(maxl,maxm)
  real vmag(maxl,maxm)
  real vort(maxl,maxm)
  real w(maxl,maxm)
  real xc(maxl,maxm)
  real xp(maxl,maxm)
  real yc(maxl,maxm)
  real yp(maxl,maxm)

  read(10,*)cost2
  read(10,*)l2

  if ( l2 /= l ) then
    write ( *, * ) ' '
    write ( *, * ) 'RSDIFF - Fatal error!'
    write ( *, * ) '  Second input file has L2 = ',l2
    write ( *, * ) '  but first input file has L  =  ',l
    stop
  end if
!
!  The next two read statements are just "dummies".
!
  read(10,*)i
  read(10,*)i
  read(10,*)m2

  if ( m2 /= m ) then
    write ( *, * ) ' '
    write ( *, * ) 'RSDIFF - Fatal error!'
    write ( *, * ) '  Second input file has M2 = ',m2
    write ( *, * ) '  but first input file has M  =  ',m
    stop
  end if

  read(10,*)i
  read(10,*)tnow2

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      e(i,j) = tempa(i,j)-e(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      gamt(i,j) = tempa(i,j)-gamt(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      p(i,j) = tempa(i,j)-p(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      pc(i,j) = tempa(i,j)-pc(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      psi(i,j) = tempa(i,j)-psi(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      rueta(i,j) = tempa(i,j)-rueta(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      ruksi(i,j) = tempa(i,j)-ruksi(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      t(i,j) = tempa(i,j)-t(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      te(i,j) = tempa(i,j)-te(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      tk(i,j) = tempa(i,j)-tk(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      v(i,j) = tempa(i,j)-v(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      u(i,j) = tempa(i,j)-u(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      vort(i,j) = tempa(i,j)-vort(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      w(i,j) = tempa(i,j)-w(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      dydp(i,j) = tempa(i,j)-yp(i,j)
    end do
  end do

  read(10,'(6g14.5)')(tempa(i,1),i = 1,nbot)

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      dycdp(i,j) = tempa(i,j)-yc(i,j)
    end do
  end do

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      dxdp(i,j) = tempa(i,j)-xp(i,j)
    end do
  end do

  read(10,'(6g14.5)')(tempa(i,1),i = 1,nbot)

  read(10,'(6g14.5)')((tempa(i,j),i = 1,l),j=1,m)
  do i = 1,l
    do j = 1,m
      dxcdp(i,j) = tempa(i,j)-xc(i,j)
    end do
  end do
!
!  Set the velocity magnitudes.
!
  do i = 1,l
    do j = 1,m
      vmag(i,j) = sqrt(u(i,j)**2+v(i,j)**2)
    end do
  end do

  return
end
subroutine rsize(delx,dely,grace,l,m,maxl,maxm,ncflag,npflag, &
  srange,vpflag,xc,x1max,x1min,x2max,x2min,xmax,xmin,xsmax,xsmin, &
  yc,y1max,y1min,y2max,y2min,ymax,ymin,ysmax,ysmin)

!*****************************************************************************80
!
!! RSIZE computes the size of the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  DELX   Input, real DELX.
!         The X spacing between nodes.  In some cases,
!         this spacing is modified to create isoparametric elements.
! 
!  DELY   Input, real DELY.
!         The Y spacing between nodes.  In some cases,
!         this spacing is modified to create isoparametric elements.
! 
!  GRACE  Input, real GRACE.
!         The size of the "grace" margin on the plot.
! 
!  L      Input, integer ( kind = 4 ) L, the number of rows of data.
!
!  M      Input, integer ( kind = 4 ) M, the number of columns of data.
!
!  MAXL,
!  MAXM   Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!         rows and columns of data.
!
!  SRANGE Input, real SRANGE.
!         The maximum of XSMAX-XSMIN and YSMAX-YSMIN.
!         This gives the size of a square containing the data
!         window.
!
!  VPFLAG Input, logical VPFLAG(MAXL,MAXM), a flag which determines
!         whether vector data at a given node should be displayed.
!
!  X1MAX,
!  X1MIN  Output, real X1MAX, X1MIN, the maximum and minimum X 
!         coordinates of the plot, which includes a small grace margin.
!
!  X2MAX,
!  X2MIN  Output, real X2MAX, X2MIN, the maximum and minimum X 
!         coordinates that should be used for plotting.  No plotting
!         commands should  exceed these values.  This is where the 
!         "frame" might be drawn.
!
!  XC     Input, real XC(MAXL,MAXM).
!         The X coordinates of the nodes.
! 
!  XMAX   Output, real XMAX.
!         The maximum X coordinate of all the nodes.
!         The maximum entry in the XC array.
! 
!  XMIN   Output, real XMIN.
!         The minimum X coordinate of all the nodes.
!         The minimum entry in the XC array.
! 
!  XSMAX  Output, real XSMAX.
!         The maximum X coordinate of the data to be displayed.
!         XSMAX defaults to XMAX, but can be made smaller to
!         focus on a portion of the region.
! 
!  XSMIN  Output, real XSMIN.
!         The minimum X coordinate of the data to be displayed.
!         XSMIN defaults to XMIN, but can be made larger to
!         focus on a portion of the region.
! 
!  Y1MAX,
!  Y1MIN  Output, real Y1MAX, Y1MIN, the maximum and minimum Y 
!         coordinates of the plot, which includes a small grace margin.
!
!  Y2MAX,
!  Y2MIN  Output, real Y2MAX, Y2MIN, the maximum and minimum Y 
!         coordinates that should be used for plotting.  No plotting 
!         commands should exceed these values.  This is where the 
!         "frame" might be drawn.
!  
!  YC     Input, real YC(MAXL,MAXM).
!         The Y coordinates of the nodes.
! 
!  YMAX   Output, real YMAX.
!         The maximum Y coordinate of all the nodes.
!         The maximum value attained by the YC array.
! 
!  YMIN   Output, real YMIN.
!         The minimum Y coordinate of all the nodes.
!         The minimum value attained by the YC array.
!  
!  YSMAX  Output, real YSMAX.
!         The maximum Y coordinate of the data to be displayed.
!         YSMAX defaults to YMAX, but can be made smaller to
!         focus on a portion of the region.
! 
!  YSMIN  Output, real YSMIN.
!         The minimum Y coordinate of the data to be displayed.
!         YSMIN defaults to YMIN, but can be made larger to
!         focus on a portion of the region.
! 
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  real delx
  real dely
  real grace
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical ncflag(maxl,maxm)
  logical npflag(maxl,maxm)
  real srange
  logical vpflag(maxl,maxm)
  real xc(maxl,maxm)
  real x1max
  real x1min
  real x2max
  real x2min
  real xmax
  real xmin
  real xsmax
  real xsmin
  real yc(maxl,maxm)
  real y1max
  real y1min
  real y2max
  real y2min
  real ymax
  real ymin
  real ysmax
  real ysmin
!
!  Reset the nodal flags to TRUE.
!
  ncflag(1:l,1:m) = .true.

  do i = 1,l
    do j = 1,m
      npflag(i,j) = .true.
    end do
  end do

  do i = 1,l
    do j = 1,m
      vpflag(i,j) = .true.
    end do
  end do
!
!  Find the maximum and minimum values of X and Y at the corner nodes.
!
  call fsize(l,m,maxl,maxm,ncflag,xc,xmax,xmin)

  call fsize(l,m,maxl,maxm,ncflag,yc,ymax,ymin)
!
!  DELX and DELY are estimates for the X and Y dimensions of the
!  region.  For various reasons, we are going to behave as though
!  the region were a square.
!
  delx = (xmax-xmin)/real(l)
  dely = (ymax-ymin)/real(m)

  write ( *, * ) ' '
  write ( *, * ) 'RSize - Note:'
  write ( *, * ) '  Physical data coordinates:'
  write ( *, * ) ' '
  write ( *, * ) xmin,' =  XMIN <= XC <= XMAX = ',xmax
  write ( *, * ) ymin,' =  YMIN <= YC <= YMAX = ',ymax
!
!  In order to allow display of data from repeated steps with a
!  fixed window, we will only reset the window if it hasn't been
!  set already.
!
  if ( xsmin == xsmax .or. ysmin == ysmax .or. &
       xmin > xsmax .or. xmax < xsmin .or. &
       ymin > ysmax .or. ymax < ysmin ) then

    write ( *, * ) ' '
    write ( *, * ) 'RSIZE - Note:'
    write ( *, * ) '  Setting data window to physical window!'
!
!  Data window starts out the same as physical window.
!
    xsmax = xmax
    xsmin = xmin
    ysmax = ymax
    ysmin = ymin
!
!  Compute box containing data.
!
    call pltbox(grace,srange,x1max,x1min,x2max,x2min,xsmax, &
      xsmin,y1max,y1min,y2max,y2min,ysmax,ysmin)

  end if

  return
end
subroutine rsread(b1jbl,b2jbl,b3jbl,cost,e,gamt,icrys1,icrys2,jcrys1,jcrys2, &
  kcrys,kmelt,kvoid,l,m,maxbot,maxl,maxm,nbot,p,pc,psi,rueta,ruksi,t, &
  te,tk,tnow,u,v,vmag,vort,w,xbot,xc,xp,ybot,yc,yp)

!*****************************************************************************80
!
!! RSREAD reads in restart information.  
!
!  Discussion:
!
!    It has to interchange X and Y data, because of some 
!    strange convention in CRYSTAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real COST, the value of the cost functional 
!    associated with this solution.
!
!    Output, real E(MAXL,MAXM), the magnetic stream function.
!
!    Output, real GAMT(MAXL,MAXM), the diffusion coefficient.
!
!    Output, integer ( kind = 4 ) ICRYS1, the I coordinate of the row of corner
!    nodes that define the bottom of the crystal.
!
!    Output, integer ( kind = 4 ) ICRYS2, the I coordinate of the row of corner
!    nodes that define the top of the crystal.
!
!    Output, integer ( kind = 4 ) JCRYS1, the J coordinate of the column of corner
!    nodes that define the left side of the crystal.
!
!    Output, integer ( kind = 4 ) JCRYS2, the J coordinate of the column of corner
!    nodes that define the right side of the crystal.
!
!    Output, integer ( kind = 4 ) KCRYS(MAXL,MAXM), 
!    0, if control volume (I,J) is away from the crystal.
!    1, if control volume (I,J) is on the external crystal boundary.
!    2, if control volume (I,J) is on the internal crystal boundary.
!    3, if control volume (I,J) is in the crystal interior.
!
!    Output, integer ( kind = 4 ) KMELT(MAXL,MAXM),
!    0, if control volume (I,J) is away from the melt.
!    1, if control volume (I,J) is on the external melt boundary.
!    2, if control volume (I,J) is on the internal melt boundary.
!    3, if control volume (I,J) is in the melt interior.
! 
!    Output, integer ( kind = 4 ) KVOID(MAXL,MAXM), 
!    0, if control volume (I,J) is away from the void.
!    1, if control volume (I,J) is on the external void boundary.
!    2, if control volume (I,J) is on the internal void boundary.
!    3, if control volume (I,J) is in the void interior.
! 
!    Output, integer ( kind = 4 ) L, the number of rows of data.
!
!    Output, integer ( kind = 4 ) M, the number of columns of data.
!
!    Input, integer ( kind = 4 ) MAXBOT, the maximum storage for the nodes
!    that define the crucible bottom.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
!    Output, integer ( kind = 4 ) NBOT, the number of crucible nodes.
!
!    Output, real P(MAXL,MAXM), the pressure.
!
!    Output, real PC(MAXL,MAXM), the corrected pressure.
!
!    Output, real PSI(MAXL,MAXM), the stream function.
!
!    Output, real RUETA(MAXL,MAXM), RUKSI(MAXL,MAXM), the
!    the momentum in the ETA and KSI directions.
!
!    Output, real T(MAXL,MAXM), the temperature.
!
!    Output, real TE(MAXL,MAXM), the turbulent epsilon.
!
!    Output, real TK(MAXL,MAXM), the turbulent K.
!
!    Output, real TNOW, the current time.
!
!    Output, real U(MAXL,MAXM), the horizontal velocity.
!
!    Output, real V(MAXL,MAXM), the vertical velocity.
!
!    Output, real VMAG(MAXL,MAXM), the velocity magnitude.
!
!    Output, real VORT(MAXL,MAXM), the vorticity.
!
!    Output, real W(MAXL,MAXM), the axial velocity.
!
!    Output, real X(MAXL,MAXM), the X coordinates of the primary nodes.
!
!    Output, real XBOT(NBOT), the X coordinates of the crucible bottom nodes.
!
!    Output, real XC(MAXL,MAXM), the X coordinates of the corner nodes.
!
!    Output, real Y(MAXL,MAXM), the Y coordinates of the primary nodes.
!
!    Output, real YBOT(NBOT), the Y coordinates of the crucible bottom nodes.
!
!    Output, real YC(MAXL,MAXM), the Y coordinates of the corner nodes.
!
  implicit none

  integer ( kind = 4 ) maxbot
  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  real b1jbl(maxm)
  real b2jbl(maxm)
  real b3jbl(maxl)
  real cost
  real e(maxl,maxm)
  real gamt(maxl,maxm)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icrys1
  integer ( kind = 4 ) icrys2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys1
  integer ( kind = 4 ) jcrys2
  integer ( kind = 4 ) kcrys(maxl,maxm)
  integer ( kind = 4 ) kmelt(maxl,maxm)
  integer ( kind = 4 ) kvoid(maxl,maxm)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  character ( len = 6 ) name
  integer ( kind = 4 ) nbot
  real p(maxl,maxm)
  real pc(maxl,maxm)
  real psi(maxl,maxm)
  real rueta(maxl,maxm)
  real ruksi(maxl,maxm)
  real t(maxl,maxm)
  real te(maxl,maxm)
  real tk(maxl,maxm)
  real tnow
  real u(maxl,maxm)
  real v(maxl,maxm)
  real vmag(maxl,maxm)
  real vort(maxl,maxm)
  real w(maxl,maxm)
  real xbot(maxbot)
  real xc(maxl,maxm)
  real xp(maxl,maxm)
  real ybot(maxbot)
  real yc(maxl,maxm)
  real yp(maxl,maxm)

  name = 'cost'
  read(10,*,end = 10,err=20)cost

  name = 'l'
  read(10,*,end = 10,err=20)l

  if ( l > maxl ) then
    write ( *, * ) ' '
    write ( *, * ) 'RSREAD - Fatal error!'
    write ( *, * ) '  Input data file has L = ',l
    write ( *, * ) '  but CRYSTAL_PLOT has MAXL  =  ',maxl
    stop
  end if

  name = 'jcrys2'
  read(10,*,end = 10,err=20)jcrys2

  name = 'icrys1'
  read(10,*,end = 10,err=20)icrys1

  name = 'm'
  read(10,*,end = 10,err=20)m

  if ( m > maxm ) then
    write ( *, * ) ' '
    write ( *, * ) 'RSREAD - Fatal error!'
    write ( *, * ) '  Input data file has M = ',m
    write ( *, * ) '  but CRYSTAL_PLOT has MAXM  =  ',maxm
    stop
  end if

  name = 'nbot'
  read(10,*,end = 10,err=20)nbot

  name = 'tnow'
  read(10,*,end = 10,err=20)tnow

  name = 'b1jbl'
  read(10,'(6g14.5)',end = 10,err=20)(b1jbl(i),i=1,m)

  name = 'b2jbl'
  read(10,'(6g14.5)',end = 10,err=20)(b2jbl(i),i=1,m)

  name = 'b3jbl'
  read(10,'(6g14.5)',end = 10,err=20)(b3jbl(i),i=1,l)

  name = 'e'
  read(10,'(6g14.5)',end = 10,err=20)((e(i,j),i=1,l),j=1,m)

  name = 'gamt'
  read(10,'(6g14.5)',end = 10,err=20)((gamt(i,j),i=1,l),j=1,m)

  name = 'p'
  read(10,'(6g14.5)',end = 10,err=20)((p(i,j),i=1,l),j=1,m)

  name = 'pc'
  read(10,'(6g14.5)',end = 10,err=20)((pc(i,j),i=1,l),j=1,m)

  name = 'psi'
  read(10,'(6g14.5)',end = 10,err=20)((psi(i,j),i=1,l),j=1,m)

  name = 'rueta'
  read(10,'(6g14.5)',end = 10,err=20)((rueta(i,j),i=1,l),j=1,m)

  name = 'ruksi'
  read(10,'(6g14.5)',end = 10,err=20)((ruksi(i,j),i=1,l),j=1,m)

  name = 't'
  read(10,'(6g14.5)',end = 10,err=20)((t(i,j),i=1,l),j=1,m)

  name = 'te'
  read(10,'(6g14.5)',end = 10,err=20)((te(i,j),i=1,l),j=1,m)

  name = 'tk'
  read(10,'(6g14.5)',end = 10,err=20)((tk(i,j),i=1,l),j=1,m)

  name = 'v'
  read(10,'(6g14.5)',end = 10,err=20)((v(i,j),i=1,l),j=1,m)

  name = 'u'
  read(10,'(6g14.5)',end = 10,err=20)((u(i,j),i=1,l),j=1,m)

  name = 'vort'
  read(10,'(6g14.5)',end = 10,err=20)((vort(i,j),i=1,l),j=1,m)

  name = 'w'
  read(10,'(6g14.5)',end = 10,err=20)((w(i,j),i=1,l),j=1,m)

  name = 'y'
  read(10,'(6g14.5)',end = 10,err=20)((yp(i,j),i=1,l),j=1,m)

  name = 'ybot'
  read(10,'(6g14.5)',end = 10,err=20)(ybot(i),i=1,nbot)

  name = 'yc'
  read(10,'(6g14.5)',end = 10,err=20)((yc(i,j),i=1,l),j=1,m)

  name = 'x'
  read(10,'(6g14.5)',end = 10,err=20)((xp(i,j),i=1,l),j=1,m)

  name = 'xbot'
  read(10,'(6g14.5)',end = 10,err=20)(xbot(i),i=1,nbot)

  name = 'xc'
  read(10,'(6g14.5)',end = 10,err=20)((xc(i,j),i=1,l),j=1,m)
!
!  Set the velocity magnitudes.
!
  do i = 1,l
    do j = 1,m
      vmag(i,j) = sqrt(u(i,j)**2+v(i,j)**2)
    end do
  end do
!
!  Set the last crystal row (L), and the first crystal column (1).
!
  icrys2 = l
  jcrys1 = 1
!
!  Set a marker for where each node is.
!  Regions are:
!
!    C - crystal
!    M - melt
!    V - void
!
  do i = 1,l
    do j = 1,m

      if ( i < icrys1-1 .or. j > jcrys2+1 ) then
        kcrys(i,j) = 0
      else if ( i == icrys1-1 .and. j <= jcrys2+1 ) then
        kcrys(i,j) = 1
      else if ( i >= icrys1-1 .and. j == jcrys2+1 ) then
        kcrys(i,j) = 1
      else if ( i == icrys1 .and. j <= jcrys2 ) then
        kcrys(i,j) = 2
      else if ( i >= icrys1 .and. j == jcrys2 ) then
        kcrys(i,j) = 2
      else if ( i == l .and. j <= jcrys2 ) then
        kcrys(i,j) = 2
      else if ( i >= icrys1 .and. j == 1 ) then
        kcrys(i,j) = 2
      else
        kcrys(i,j) = 3
      end if

    end do
  end do

  do i = 1,l
    do j = 1,m

      if ( i > icrys1 ) then
        kmelt(i,j) = 0
      else if ( i == icrys1 ) then
        kmelt(i,j) = 1
      else if ( i == icrys1-1 ) then
        kmelt(i,j) = 2
      else if ( i <= icrys1-1 .and. j == m ) then
        kmelt(i,j) = 2
      else if ( i == 1 ) then
        kmelt(i,j) = 2
      else
        kmelt(i,j) = 3
      end if

    end do
  end do

  do i = 1,l
    do j = 1,m

      if ( i < icrys1-1 ) then
        kvoid(i,j) = 0
      else if ( j < jcrys2 ) then
        kvoid(i,j) = 0
      else if ( i == icrys1-1 .and. j >= jcrys2 ) then
        kvoid(i,j) = 1
      else if ( i >= icrys1-1 .and. j == jcrys2 ) then
        kvoid(i,j) = 1
      else if ( i == icrys1 .and. j >= jcrys2+1 ) then
        kvoid(i,j) = 2
      else if ( i >= icrys1 .and. j == jcrys2+1 ) then
        kvoid(i,j) = 2
      else
        kvoid(i,j) = 3
      end if
    end do
  end do
 
  return

10    continue
  write ( *, * ) ' '
  write ( *, * ) 'RSREAD - Fatal error!'
  write ( *, * ) '  End-of-file while reading data.'
  write ( *, * ) '  RSREAD was trying to read variable ' // name
  stop

20    continue
  write ( *, * ) ' '
  write ( *, * ) 'RSREAD - Fatal error!'
  write ( *, * ) '  READ error while reading data.'
  write ( *, * ) '  RSREAD was trying to read ' // name
  stop
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
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
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) nchar
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
subroutine s_cap ( s )

!*****************************************************************************80
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

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) s

  nchar = len_trim ( s )

  do i = 1, nchar

    c = s(i:i)
    call ch_cap ( c )
    s(i:i) = c

  end do

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
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

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

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
subroutine s_plot ( angle, cwide, pwide, s, x, y, flush )

!*****************************************************************************80
!
!! S_PLOT plots a character string onto a graphics image.
!
!  Discussion:
!
!    The string can be at any angle and at any size.
!
!    The plot is assumed to be of size PWIDE by PHITE, although PHITE
!    itself is not input.
!
!    This routine must be modified to work with a particular graphics package.
!    The current code calls two routines:
!      MOVCGM ( X, Y ) moves to a point (X,Y) in the plot;
!      DRWCGM ( X, Y ) draws a line from the current point to (X,Y).
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
!    Input, real ANGLE, the angle in degrees at which the
!    string is to be drawn.  0 is typical.  90 degrees would
!    cause the string to be written from top to bottom.
!
!    Input, real CWIDE, the width of the characters.  This
!    is measured in the same units as the plot width PWIDE.
!    For PWIDE = 1, a plot size of 0.025 would be reasonable,
!    since 40 characters would fit, but 2.0 would be nonsense.
!
!    Input, real PWIDE, the width of the plot, in the same
!    units as CWIDE.
!
!    Input, character ( len = * ) S, contains the text to be plotted.
!    Only characters with ASCII codes between 32 and 126 will actually
!    be plotted.  Any other characters are "unprintable", and will be
!    plotted as blanks.
!
!    Input, real X, Y, the coordinates of a point which
!    determines where the string is drawn.  The string will
!    be drawn starting at, centered or, or ending at (X,Y),
!    depending on the value of FLUSH.
!
!    Input, character ( len = * ) FLUSH, a string which specifies how to place
!    the string.  Only the first character of FLUSH is examined, and the case of
!    the character is not important.
!
!    'L' - the string will be drawn flush left.
!    'C' - the string will be centered.
!    'R' - the string will be drawn flush right.
!
  implicit none

  real, parameter :: PI = 3.1415926535
  real, parameter :: DEG_TO_RAD = PI / 180.0

  real angle
  real ca
  character c
  real cwide
  character ( len = * ) flush
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iascii
  integer ( kind = 4 ) icr
  integer ( kind = 4 ), save, dimension ( 1617 ) :: ifont
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipen
  integer ( kind = 4 ), save, dimension ( 95 ) :: ipoint
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) nmax
  integer ( kind = 4 ) nvec
  real pwide
  logical rotate
  character ( len = * ) s
  real sa
  real scl2
  real x
  real xb
  real xc
  real xcopy
  real xnew
  real xold
  real xrot
  real xt
  real y
  real yb
  real yc
  real ycopy
  real ynew
  real yold
  real yrot
  real yt
!
!  IPOINT is a pointer array into IFONT.
!
!  IPOINT(I) records where the "strokes" for character I begin
!  in the IFONT array.
!
  data ( ipoint(i), i = 1, 95 ) / &
       1,   3,  26,  45,  66, 102, 130, 156, 166, 186, 206, 222, 233, &
     249, 255, 267, 273, 293, 306, 328, 353, 363, 383, 411, 423, 457, &
     483, 506, 533, 541, 552, 560, 587, 625, 638, 665, 683, 699, 714, &
     727, 754, 770, 786, 805, 818, 826, 838, 848, 868, 884, 909, 930, &
     956, 967, 981, 989,1001,1012,1025,1035,1045,1051,1061,1069,1075, &
    1081,1108,1131,1149,1172,1194,1214,1243,1260,1284,1307,1323,1336, &
    1364,1381,1401,1424,1447,1464,1486,1499,1516,1524,1536,1547,1560, &
    1570,1586,1592,1608 /
!
!  IFONT contains the strokes defining the various symbols.
!
  data ( ifont(i), i = 1, 396 ) / &
     1, 0, 2,10,11, 9,22,10,23,11,22,10,11, 0, 9, 7, 9, 9,11, 9,11, 7, &
     9, 7, 0, 2, 8,17, 7,23, 9,23, 8,17, 0,14,17,13,23,15,23,14,17, 0, &
     4, 9,23, 7, 7, 0,13,23,11, 7, 0, 5,17,15,17, 0, 5,13,15,13, 0, 3, &
    15,19,13,21, 9,21, 7,19, 7,17, 9,15,13,15,15,13,15,11,13, 9, 9, 9, &
     7,11, 0, 9,23, 9, 7, 0,13,23,13, 7, 0, 3, 5,23, 9,23, 9,19, 5,19, &
     5,23, 0,17,23, 5, 7, 0,13, 7,13,11,17,11,17, 7,13, 7, 0, 1,17, 7, &
     7,17, 7,19, 9,21,13,21,15,19,15,17, 5,13, 5,11, 9, 7,13, 7,17,15, &
     0, 1,10,17, 9,23,11,23,10,17, 0, 1,12,23,11,21,10,19, 9,17, 9,15, &
     9,13,10,11,11, 9,12, 7, 0, 1,12,23,13,21,14,19,15,17,15,15,15,13, &
    14,11,13, 9,12, 7, 0, 3, 7,15,15,15, 0,13,19, 9,11, 0, 9,19,13,11, &
     0, 2, 7,15,15,15, 0,11,19,11,11, 0, 1,11, 7, 9, 7, 9, 9,11, 9,11, &
     7,11, 6,10, 4, 0, 1, 7,15,15,15, 0, 1, 9, 7, 9, 9,11, 9,11, 7, 9, &
     7, 0, 1,15,23, 7, 7, 0, 1, 9,23,13,23,15,19,15,11,13, 7, 9, 7, 7, &
    11, 7,19, 9,23, 0, 2, 7,21, 9,23, 9, 7, 0, 7, 7,11, 7, 0, 1, 5,21, &
     9,23,15,23,17,21,17,19,15,17, 7,13, 5,10, 5, 7,17, 7, 0, 2, 5,23, &
    17,23,15,17,13,15, 9,15, 0,13,15,17,13,17,10,14, 7, 8, 7, 5,10, 0, &
     1,13, 7,13,23, 5,13,17,13, 0, 1,17,23, 5,23, 5,17,13,17,17,15,17, &
    11,13, 7, 9, 7, 5,11, 0, 1,17,19,13,23, 9,23, 5,19, 5,13, 9,15,13 /

  data ( ifont(i), i =  397, 792 ) / &
    15,17,13,17,11,13, 7, 9, 7, 5,11, 5,13, 0, 1, 5,19, 5,23,17,23,11, &
    15,11, 7, 0, 1, 8,15, 6,17, 6,21, 8,23,14,23,16,21,16,17,14,15, 8, &
    15, 5,13, 5, 9, 8, 7,14, 7,17, 9,17,13,14,15, 0, 1,17,17,15,15, 7, &
    15, 5,17, 5,21, 7,23,15,23,17,21,17,11,15, 7, 7, 7, 5,11, 0, 2, 9, &
    13, 9,15,11,15,11,13, 9,13, 0, 9, 7, 9, 9,11, 9,11, 7, 9, 7, 0, 2, &
     9,13, 9,15,11,15,11,13, 9,13, 0,11, 7, 9, 7, 9, 9,11, 9,11, 7,11, &
     6,10, 4, 0, 1,17,21, 5,15,17, 9, 0, 2, 7,15,15,15, 0, 7, 9,15, 9, &
     0, 1, 5,21,17,15, 5, 9, 0, 2, 7,21, 9,23,13,23,15,21,15,19,11,15, &
    11,11, 0,10, 7,10, 9,12, 9,12, 7,10, 7, 0, 1,13, 7, 9, 7, 5,11, 5, &
    19, 9,23,13,23,17,19,17,11,15, 9,13,11,12,10,10,10, 9,11, 9,15,10, &
    16,12,16,13,15,13,11, 0, 2, 5, 7,11,23,17, 7, 0, 8,15,14,15, 0, 2, &
     5, 7, 5,23,15,23,17,21,17,17,15,15, 5,15, 0,15,15,17,13,17, 9,15, &
     7, 5, 7, 0, 1,17,19,13,23, 9,23, 5,19, 5,11, 9, 7,13, 7,17,11, 0, &
     1, 5, 7, 5,23,13,23,17,19,17,11,13, 7, 5, 7, 0, 2,17,23, 5,23, 5, &
     7,17, 7, 0, 5,15,12,15, 0, 2, 5, 7, 5,23,17,23, 0, 5,15,12,15, 0, &
     2,17,19,13,23, 9,23, 5,19, 5,11, 9, 7,13, 7,17,11,17,15,13,15, 0, &
    17,11,17, 7, 0, 3, 5, 7, 5,23, 0, 5,15,17,15, 0,17,23,17, 7, 0, 3, &
     9,23,13,23, 0,11,23,11, 7, 0, 9, 7,13, 7, 0, 2,15,23,15,11,12, 7 /

  data ( ifont(i), i =  793, 1188 ) / &
     8, 7, 5,11, 5,13, 0,13,23,17,23, 0, 2, 5, 7, 5,23, 0,17,23, 5,15, &
    17, 7, 0, 1, 5,23, 5, 7,17, 7, 0, 1, 5, 7, 5,23,11,11,17,23,17, 7, &
     0, 1, 5, 7, 5,23,17, 7,17,23, 0, 1,17,19,13,23, 9,23, 5,19, 5,11, &
     9, 7,13, 7,17,11,17,19, 0, 1, 5, 7, 5,23,13,23,17,21,17,17,13,15, &
     5,15, 0, 2,17,19,13,23, 9,23, 5,19, 5,11, 9, 7,13, 7,17,11,17,19, &
     0,13,11,17, 7, 0, 2, 5, 7, 5,23,13,23,17,21,17,17,13,15, 5,15, 0, &
    13,15,17, 7, 0, 1,17,19,13,23, 9,23, 5,20, 5,18, 9,15,13,15,17,12, &
    17,10,13, 7, 9, 7, 5,10, 0, 2, 5,23,17,23, 0,11,23,11, 7, 0, 1, 5, &
    23, 5,10, 8, 7,14, 7,17,10,17,23, 0, 1, 5,23,11, 7,17,23, 0, 1, 5, &
    23, 8, 7,11,17,14, 7,17,23, 0, 2, 5,23,17, 7, 0,17,23, 5, 7, 0, 2, &
     5,23,11,13,17,23, 0,11,13,11, 7, 0, 1, 5,23,17,23, 5, 7,17, 7, 0, &
     1,11,23, 7,23, 7, 7,11, 7, 0, 1, 7,23,15, 7, 0, 1, 7,23,11,23,11, &
     7, 7, 7, 0, 1, 7,21,11,23,15,21, 0, 1, 5, 3,17, 3, 0, 1, 9,23,13, &
    19, 0, 2, 7,14, 9,15,13,15,15,14,15, 7, 0,15,12, 9,12, 7,11, 7, 8, &
     9, 7,13, 7,15, 8, 0, 2, 7,23, 7, 7, 0, 7,13, 9,15,13,15,15,13,15, &
     9,13, 7, 9, 7, 7, 9, 0, 1,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7,13, &
     7,15, 9, 0, 2,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7,13, 7,15, 9, 0, &
    15,23,15, 7, 0, 1, 7,11,15,11,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7 /

  data ( ifont(i), i = 1189, 1584 ) / &
    13, 7,15, 9, 0, 3, 9, 7, 9,23,13,23,13,22, 0, 8,15,12,15, 0, 8, 7, &
    11, 7, 0, 2,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7,13, 7,15, 9, 0,15, &
    13,15, 3,13, 1, 9, 1, 7, 3, 0, 2, 7, 7, 7,23, 0, 7,14, 9,15,13,15, &
    15,14,15, 7, 0, 3, 9,15,11,15,11, 7, 0, 9, 7,13, 7, 0, 9,17, 9,19, &
    11,19,11,17, 9,17, 0, 2, 9,15,11,15,11, 1, 7, 1, 7, 3, 0, 9,17,11, &
    17,11,19, 9,19, 9,17, 0, 3, 7, 7, 7,23, 0,15,15, 7,10, 0, 9,11,15, &
     7, 0, 2, 9,23,11,23,11, 7, 0, 9, 7,13, 7, 0, 3, 7,15, 7, 7, 0, 7, &
    14, 8,15,10,15,11,14,11, 7, 0,11,14,12,15,14,15,15,14,15, 7, 0, 2, &
     7, 7, 7,15, 0, 7,14, 9,15,13,15,15,14,15, 7, 0, 1, 7,13, 9,15,13, &
    15,15,13,15, 9,13, 7, 9, 7, 7, 9, 7,13, 0, 2, 7,13, 9,15,13,15,15, &
    13,15, 9,13, 7, 9, 7, 7, 9, 0, 7,14, 7, 1, 0, 2,15,13,13,15, 9,15, &
     7,13, 7, 9, 9, 7,13, 7,15, 9, 0,15,14,15, 1, 0, 2, 7,15, 9,15, 9, &
     7, 0, 9,13,11,15,13,15,15,13, 0, 1,15,13,13,15, 9,15, 7,13, 9,11, &
    13,11,15, 9,13, 7, 9, 7, 7, 9, 0, 2, 9,23, 9, 7,11, 7, 0, 7,17,11, &
    17, 0, 2, 7,15, 7, 9, 9, 7,13, 7,15, 9, 0,15,15,15, 7, 0, 1, 7,15, &
    11, 7,15,15, 0, 1, 7,15, 9, 7,11,11,13, 7,15,15, 0, 2, 7,15,15, 7, &
     0, 7, 7,15,15, 0, 2, 7,15,11, 7, 0,15,15,10, 5, 7, 1, 0, 1, 7,15, &
    15,15, 7, 7,15, 7, 0, 1,11,23, 7,23, 9,17, 7,15, 9,13, 7, 7,11, 7 /

  data ( ifont(i), i = 1585, 1617 ) / &
     0, 1, 9,23, 9, 7, 0, 1, 7,23,11,23, 9,17,11,15, 9,13,11, 7, 7, 7, &
     0, 1, 5,21, 7,23,15,21,17,23, 0 /

  nchar = len_trim ( s )

  if ( pwide <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S_PLOT - Serious error!'
    write ( *, * ) '  The plot width PWIDE is negative!'
    write ( *, * ) '  PWIDE = ', pwide
    return
  end if
!
!  Chop titles that are too long.  To do this, we need to know the
!  width of the plot (PWIDE) in same units as CWIDE.
!
  nmax = int ( pwide / cwide )

  if ( nchar > nmax ) then
    nchar = nmax
  end if
!
!  Shift string if centering or right flush option used.
!
  if ( flush(1:1) == 'l' .or. flush(1:1) == 'L' ) then
    xcopy = x
    ycopy = y
  else if ( flush(1:1) == 'c' .or. flush(1:1) == 'C' ) then
    xcopy = x - 0.5 * nchar * cwide * cos ( angle * DEG_TO_RAD )
    ycopy = y - 0.5 * nchar * cwide * sin ( angle* DEG_TO_RAD )
  else if ( flush(1:1) == 'r' .or. flush(1:1) == 'R' ) then
    xcopy = x - nchar * cwide * cos ( angle * DEG_TO_RAD )
    ycopy = y - nchar * cwide * sin ( angle * DEG_TO_RAD )
  else
    xcopy = x
    ycopy = y
  end if
!
!  Note that screen coordinates are used.
!  Thus a width of 0.1 is intended to mean 1/10 of screen size.
!
!  Set the scale factor for character height.
!
  scl2 = cwide / 16.0
!
!  Set the starting point for the line of text, the lower left
!  corner of the first character.
!
!  Set the origin about which rotation is performed.
!
  xb = xcopy
  xrot = xcopy
  yb = ycopy
  yrot = ycopy
!
!  Get trig functions if rotation required, converting from
!  degrees to radians.
!
  if ( angle == 0.0 ) then
    rotate = .false.
  else
    ca = cos ( angle * DEG_TO_RAD )
    sa = sin ( angle * DEG_TO_RAD )
    rotate = .true.
  end if
!
!  Loop over all characters in the string.
!
  do icr = 1, nchar

    xold = x
    yold = y
    xnew = x
    ynew = y
!
!  Get the ASCII code for the character and shift by 31 so that
!  the first printable character becomes code 1.
!
    c = s(icr:icr)
    iascii = ichar ( c ) - 31
!
!  Replace any nonprintable characters with blanks.
!
    if ( iascii < 1 .or. iascii > 95 ) then
      iascii = 1
    end if
!
!  Get the pointer to this character in font table.
!
    ip = ipoint(iascii)
!
!  Get the number of "vectors" required to draw the character.
!  Here "vectors" means the number of times the pen is lowered, not
!  the number of pen strokes.
!
!  For blanks, this number is 1, due to the way the
!  algorithm is coded.
!
    nvec = ifont(ip)
!
!  Loop over all required pen movements.
!
    do iv = 1, nvec

      ipen = 3
      ip = ip + 1

10        continue

      if ( ifont(ip) == 0 ) then
        go to 20
      end if

      xc = xb + ( ifont(ip) - 1 ) * scl2
      yc = yb + ( ifont(ip+1) - 7 ) * scl2
!
!  Apply rotation if necessary.
!
      if ( rotate ) then
        xt = xc - xrot
        yt = yc - yrot
        xc = xrot + ca * xt - sa * yt
        yc = yrot + sa * xt + ca * yt
      end if
!
!  Plot the pen stroke.
!
      if ( ipen == 3 ) then
        xnew = xc
        ynew = yc
      else
        xold = xnew
        yold = ynew
        xnew = xc
        ynew = yc
!
!  Call the user supplied routine to draw a line from
!  (XOLD,YOLD) to (XNEW,YNEW).
!
        call movcgm ( xold, yold )
        call drwcgm ( xnew, ynew )

      end if

      ipen = 2
      ip = ip + 2
      go to 10

20    continue

    end do
!
!  Advance the base to compensate for character just drawn.
!
    xb = xb + cwide

  end do

  return
end
subroutine setsiz(echo,name,smax,smin)

!*****************************************************************************80
!
!! SETSIZ allows the user to adjust the range of a contour plot.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*(*), the name of the variable.
!
!  SMAX,
!  SMIN   Input/output, real SMAX, SMIN.  
!
!         On input, SMAX and SMIN are the maximum and minimum values
!         attained by the variable.
!
!         On output, SMAX and SMIN may have been changed by the
!         user to more suitable values.
!
  implicit none

  logical echo
  character isay
  character ( len = * )  name
  logical s_eqi
  real smax
  real smin
  real tmax
  real tmin

  write ( *, * ) ' '
  write(*,'(1x,a)')name
  write ( *, * ) ' '
  write ( *, * ) 'SetSiz - Current range is ',smin,' to ',smax
  write ( *, * ) ' '
  write ( *, * ) '  Do you want to change the range?'
  read(*,'(a)',err = 10,end=10)isay
  write(17,'(a)')isay
  if ( echo ) then
    write(*,'(a)')isay
  end if
  if ( .not.s_eqi ( isay,'y'))return

  write ( *, * ) ' '
  write ( *, * ) '  Enter new minimum and maximum values for contours.'
  read(*,*,end = 10,err=10)tmin,tmax
  write(17,*)tmin,tmax
  if ( echo ) then
    write ( *, * ) tmin,tmax
  end if
  smin = tmin
  smax = tmax
  return
 
10    continue

  write ( *, * ) ' '
  write ( *, * ) 'SetSiz - Warning!'
  write ( *, * ) '  There was trouble reading your input,'
  write ( *, * ) '  so it was ignored.'
 
  return
end
subroutine settab(echo,icmax,icmin,itable)

!*****************************************************************************80
!
!! SETTAB replaces SETCTB, the DRAWCGM routine for setting up
!  the color tables.
!
!  Discussion:
!
!    So SETTAB sets the colors between ICMIN and ICMAX, which
!    should typically be 2 and 255.
!
!    SETTAB will also set the values of color 0 to white, and
!    color 1 to black.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ICMAX, the maximum color index to be set.
!
!  ICMIN  Input, integer ( kind = 4 ) ICMIN, the minimum color index to be set.
!
!  ITABLE Input, integer ( kind = 4 ) ITABLE, the desired color table.
!
!         1: low black to high white
!         2: low blue to high yellow
!         3: low red, high blue, with bands between.
!         4: low red, yellow, green, blue, high white.
!         5: low white, blue, green, yellow, high red.
!         6: low blue to high red
!         7: linear table between 2 user colors.
!         8: linear table between n user colors.
!         9: low white to high black.
!
  implicit none

  real bhi
  real blo
  real bval
  logical echo
  real ghi
  real glo
  real gval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icmax
  integer ( kind = 4 ) icmin
  integer ( kind = 4 ) icol1
  integer ( kind = 4 ) icol2
  integer ( kind = 4 ) itable
  real, parameter :: pi = 3.1415926E+00
  real rhi
  real rlo
  real rval
  real theta
!
!  1: Low black to high white
!
  if ( itable == 1 ) then
    do i = icmin,icmax
      bval = real(i-icmin)/real(icmax-icmin)
      gval = real(i-icmin)/real(icmax-icmin)
      rval = real(i-icmin)/real(icmax-icmin)
      call setclr(i,bval,gval,rval)
    end do
!
!  2: Low blue to high yellow.
!
  else if ( itable == 2 ) then
    do i = icmin,icmax
      rval = real(i-icmin)/real(icmax-icmin)
      gval = real(i-icmin)/real(icmax-icmin)
      bval = (icmax-i)/real(icmax-icmin)
      call setclr(i,bval,gval,rval)
    end do
!
!  3: Low red, high blue, with bands between.
!
  else if ( itable == 3 ) then
    do i = icmin,icmax
      theta = 0.5*pi*real(i-icmin)/real(icmax-icmin)
      rval = cos(theta)**2
      bval = sin(theta)**2
      gval = 0.8*sin(10.0*theta)**6
      call setclr(i,bval,gval,rval)
    end do
!
!  4: Low red, yellow, green, blue, high white.
!
  else if ( itable == 4 ) then
    do i = icmin,icmax
      theta = 4.0*real(i-icmin)/real(icmax-icmin)
      rval = exp(-(theta-1.0)**2)+exp(-(theta-4.0)**2)
      gval = exp(-(theta-2.0)**2)+exp(-(theta-4.0)**2)
      bval = exp(-(theta-3.0)**2)+exp(-(theta-4.0)**2)
      if ( rval > 1.0)rval = 1.0
      if ( gval > 1.0)gval = 1.0
      if ( bval > 1.0)bval = 1.0
      call setclr(i,bval,gval,rval)
    end do
!
!  5: Low white, blue, green, yellow, high red.
!
  else if ( itable == 5 ) then
    do i = icmin,icmax
      theta = 4.0*real(icmax-i)/real(icmax-icmin)
      rval = exp(-(theta-1.0)**2)+exp(-(theta-4.0)**2)
      gval = exp(-(theta-2.0)**2)+exp(-(theta-4.0)**2)
      bval = exp(-(theta-3.0)**2)+exp(-(theta-4.0)**2)
      if ( rval > 1.0)rval = 1.0
      if ( gval > 1.0)gval = 1.0
      if ( bval > 1.0)bval = 1.0
      call setclr(i,bval,gval,rval)
    end do
!
!  6: Low blue to high red
!
  else if ( itable == 6 ) then
    do i = icmin,icmax
      rval = real(i-icmin)/real(icmax-icmin)
      gval = 0.0
      bval = (icmax-i)/real(icmax-icmin)
      call setclr(i,bval,gval,rval)
    end do
!
!  7: Interpolate between two values.
!
  else if ( itable == 7 ) then

1996    continue

    write ( *, * ) ' '
    write ( *, * ) 'Enter (Rlo, Glo, Blo), (Rhi, Ghi, Bhi)'
    write ( *, * ) 'Note: (0,0,0) is black, (1,1,1) is white!'
    write ( *, * ) ' '

    read(*,*,end = 1952,err=1964)rlo,glo,blo,rhi,ghi,bhi

    write(17,*)rlo,glo,blo,rhi,ghi,bhi
    if ( echo ) then
      write ( *, * ) rlo,glo,blo,rhi,ghi,bhi
    end if
    if ( rlo < 0.0)rlo = 0.0
    if ( rhi > 1.0)rhi = 1.0
    if ( glo < 0.0)glo = 0.0
    if ( ghi > 1.0)ghi = 1.0
    if ( blo < 0.0)blo = 0.0
    if ( bhi > 1.0)bhi = 1.0
    do i = icmin,icmax
      rval = (rlo*(icmax-i)+rhi*(i-icmin))/real(icmax-icmin)
      gval = (glo*(icmax-i)+ghi*(i-icmin))/real(icmax-icmin)
      bval = (blo*(icmax-i)+bhi*(i-icmin))/real(icmax-icmin)
      call setclr(i,bval,gval,rval)
    end do
!
!  8: Interpolate between several values.
!
  else if ( itable == 8 ) then

    icol1 = icmin
    write ( *, * ) 'Enter (R, G, B) for color index ',icol1
    write ( *, * ) '      (0, 0, 0) is black.'
    write ( *, * ) '      (1, 1, 1) is white.'
    read(*,*)rlo,glo,blo
    write(17,*)rlo,glo,blo
    if ( echo ) then
      write ( *, * ) rlo,glo,blo
    end if
    if ( rlo < 0.0)rlo = 0.0
    if ( glo < 0.0)glo = 0.0
    if ( blo < 0.0)blo = 0.0

10      continue

    write ( *, * ) 'Enter index of next color to set'
    write ( *, * ) 'between ',icol1+1,' and ',icmax
    read(*,*)icol2
    write(17,*)icol2
    if ( echo ) then
      write ( *, * ) icol2
    end if

    if ( icol2 <= icol1 .or. icol2 > icmax ) then
      write ( *, * ) 'SetTab - Warning!'
      write ( *, * ) '  Your color index was not accepted!'
      go to 10
    end if

    write ( *, * ) ' '
    write ( *, * ) 'Enter (R, G, B) for color index ',icol2
    read(*,*)rhi,ghi,bhi
    write(17,*)rhi,ghi,bhi
    if ( echo ) then
      write ( *, * ) rhi,ghi,bhi
    end if

    if ( rhi > 1.0)rhi = 1.0
    if ( ghi > 1.0)ghi = 1.0
    if ( bhi > 1.0)bhi = 1.0

    do i = icol1,icol2
      rval = (rlo*(icol2-i)+rhi*(i-icol1))/real(icol2-icol1)
      gval = (glo*(icol2-i)+ghi*(i-icol1))/real(icol2-icol1)
      bval = (blo*(icol2-i)+bhi*(i-icol1))/real(icol2-icol1)
      call setclr(i,bval,gval,rval)
    end do

    if ( icol2 < icmax ) then
      icol1 = icol2
      rlo = rhi
      glo = ghi
      blo = bhi
      go to 10
    end if
!
!  9: Low white to high black
!
  else if ( itable == 9 ) then
    do i = icmin,icmax
      bval = real(icmax-i)/real(icmax-icmin)
      gval = real(icmax-i)/real(icmax-icmin)
      rval = real(icmax-i)/real(icmax-icmin)
      call setclr(i,bval,gval,rval)
    end do
!
!  Unknown table.
!
  else
    write ( *, * ) ' '
    write ( *, * ) 'SetTab - Fatal error!'
    write ( *, * ) '  Legal color table indices are '
    write ( *, * ) '  between 1 and 8.  Your value was ',itable
  end if
!
!  Background color 0 is to be white.
!
  i = 0
  rval = 1.0
  gval = 1.0
  bval = 1.0
  call setclr(i,rval,gval,bval)
!
!  Foreground color 1 is to be black.
!
  i = 1
  rval = 0.0
  gval = 0.0
  bval = 0.0
  call setclr(i,rval,gval,bval)
  return

1952  continue
  write ( *, * ) ' '
  write ( *, * ) 'SETTAB - Fatal error!'
  write ( *, * ) '  Unexpected end of file!'
  stop

1964  continue
  write ( *, * ) ' '
  write ( *, * ) 'SETTAB - Warning!'
  write ( *, * ) '  Illegal format for input data!'
  write ( *, * ) '  Try again!'
  return
end
subroutine showcv ( l, m, maxl, maxm, nflag, xc, yc )

!*****************************************************************************80
!
!! SHOWCV shows the control volumes.
!
!  Discussion:
!
!    If ALL of the four corner nodes are visible, then this routine
!    will draw all four sides of the control volume.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
!    Input, logical NFLAG(MAXL,MAXM).
!    NFLAG is used to "flag" which nodes are active,
!    that is, to be displayed, and which not, in the graph.
!
!    Input, real XC(MAXL,MAXM).
!    The X coordinates of the corner nodes.
! 
!    Input, real YC(MAXL,MAXM).
!    The Y coordinates of the corner nodes.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical nflag(maxl,maxm)
  real xc(maxl,maxm)
  real yc(maxl,maxm)

  do i = 1,l-1
    do j = 1,m-1
      if ( nflag(i,j) .and. nflag(i,j+1) .and. nflag(i+1,j) .and. &
           nflag(i+1,j+1) ) then
        call movcgm(xc(i,j),yc(i,j))
        call drwcgm(xc(i,j+1),yc(i,j+1))
        call drwcgm(xc(i+1,j+1),yc(i+1,j+1))
        call drwcgm(xc(i+1,j),yc(i+1,j))
        call drwcgm(xc(i,j),yc(i,j))
      end if
    end do
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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

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
subroutine trichk(i1,j1,i2,j2,i3,j3,maxl,maxm,s,sval,xp,yp)

!*****************************************************************************80
!
!! TRICHK is given a triangle, formed by nodes (I1,J1), (I2,J2), and
!  (I3,J3), and the value of the quantity S at these nodes.
!
!  Discussion:
!
!    TRICHK must draw a contour line for S = SVAL through this triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1, J1, I2, J2, I3, J3, the I and J indices
!    of the three nodes which form the triangle.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
!    Input real S(MAXL,MAXM), the value of S at each node.
!
!    Input, real SVAL, the value of S for which the contour
!    line is desired.
!
!    Input, real XP(MAXL,MAXM).
!    The X coordinates of the nodes.
! 
!    Input, real YP(MAXL,MAXM).
!    The Y coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) icross
  logical inside
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  real s(maxl,maxm)
  real sval
  real xp(maxl,maxm)
  real xx(3)
  real yp(maxl,maxm)
  real yy(3)

  icross = 0

  if ( inside(s(i1,j1),sval,s(i2,j2)) ) then
    icross = icross+1
    xx(icross) = xp(i1,j1)+(sval-s(i1,j1))*(xp(i2,j2)-xp(i1,j1)) &
      /(s(i2,j2)-s(i1,j1))
    yy(icross) = yp(i1,j1)+(sval-s(i1,j1))*(yp(i2,j2)-yp(i1,j1)) &
      /(s(i2,j2)-s(i1,j1))
  end if

  if ( inside(s(i2,j2),sval,s(i3,j3)) ) then
    icross = icross+1
    xx(icross) = xp(i2,j2)+(sval-s(i2,j2))*(xp(i3,j3)-xp(i2,j2)) &
      /(s(i3,j3)-s(i2,j2))
    yy(icross) = yp(i2,j2)+(sval-s(i2,j2))*(yp(i3,j3)-yp(i2,j2)) &
      /(s(i3,j3)-s(i2,j2))
  end if

  if ( inside(s(i3,j3),sval,s(i1,j1)) ) then
    icross = icross+1
    xx(icross) = xp(i3,j3)+(sval-s(i3,j3))*(xp(i1,j1)-xp(i3,j3)) &
      /(s(i1,j1)-s(i3,j3))
    yy(icross) = yp(i3,j3)+(sval-s(i3,j3))*(yp(i1,j1)-yp(i3,j3)) &
     /(s(i1,j1)-s(i3,j3))
  end if

  if ( icross == 2 ) then
    call plylin(2,xx,yy)
  end if

  return
end
subroutine tricol(jcmax,jcmin,l,m,maxl,maxm,ncon,nflag,s,smax, &
  smin,xp,yp)

!*****************************************************************************80
!
!! TRICOL uses color to indicate all the points which have a
!  function value greater than a given value.
!
!  Discussion:
!
!    TRICOL is used for quantities associated with the three
!    corner nodes of a 6 node finite element, in particular,
!    pressure.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) JCMAX, JCMIN, the maximum and 
!    minimum color indices to use in the color bar.
!
!    Input, integer ( kind = 4 ) L, M, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
!    Input, integer ( kind = 4 ) NCON, the number of contour lines to draw.
!
!    Input, logical NFLAG(MAXL,MAXM).
!    NFLAG is used to "flag" which nodes are active,
!    that is, to be displayed, and which not, in the graph.
!
!    Input real S(MAXL,MAXM), the value of S at each node.
!
!    Input, real SMAX, SMIN, the maximum and minimum
!    values of the quantity whose color contours are
!    being drawn.  These numbers will be printed along
!    with the color bar.
!
!    Input, real XP(MAXL,MAXM).
!    The X coordinates of the nodes.
! 
!    Input, real YP(MAXL,MAXM).
!    The Y coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) jcmax
  integer ( kind = 4 ) jcmin
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ncon
  logical nflag(maxl,maxm)
  real s(maxl,maxm)
  real smax
  real smin
  real xp(maxl,maxm)
  real yp(maxl,maxm)
!
!  Draw the contour line by searching over each element.
!
  do i = 1,l-1
    do j = 1,m-1

      if ( nflag(i,j+1) .and. nflag(i+1,j).and.nflag(i,j) ) then
        i1 = i
        j1 = j
        i2 = i
        j2 = j+1
        i3 = i+1
        j3 = j
        call trilhk(i1,j1,i2,j2,i3,j3,jcmax,jcmin,maxl,maxm, &
          ncon,s,smax,smin,xp,yp)
      end if

      if ( nflag(i,j+1) .and. nflag(i+1,j).and.nflag(i+1,j+1) ) then
        i1 = i+1
        j1 = j+1
        i2 = i
        j2 = j+1
        i3 = i+1
        j3 = j
        call trilhk(i1,j1,i2,j2,i3,j3,jcmax,jcmin,maxl,maxm, &
          ncon,s,smax,smin,xp,yp)
      end if

    end do
  end do

  return
end
subroutine tricon(l,m,maxl,maxm,nflag,s,sval,xp,yp)

!*****************************************************************************80
!
!! TRICON draws a single contour line in a linear finite element.
!
!  Discussion:
!
!    The contour line is the set of points (X,Y) for which S(X,Y) = SVAL,
!    where S is some scalar quantity, generally pressure, associated with 
!    the three corner nodes of a triangular finite element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
!    Input, logical NFLAG(MAXL,MAXM).
!    NFLAG is used to "flag" which nodes are active,
!    that is, to be displayed, and which not, in the graph.
!
!    Input real S(MAXL,MAXM), the value of S at each node.
!
!    Input, real SVAL, the value of S for which the contour
!    line is desired.
!
!    Input, real XP(MAXL,MAXM).
!    The X coordinates of the nodes.
! 
!    Input, real YP(MAXL,MAXM).
!    The Y coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical nflag(maxl,maxm)
  real s(maxl,maxm)
  real sval
  real xp(maxl,maxm)
  real yp(maxl,maxm)
!
!  Draw the contour line by searching over each element.
!
  do i = 1,l-1
    do j = 1,m-1

      if ( nflag(i,j+1) .and. nflag(i+1,j) ) then

        if ( nflag(i,j) ) then
          call trichk(i,j,i,j+1,i+1,j,maxl,maxm,s,sval,xp,yp)
        end if

        if ( nflag(i+1,j+1) ) then
          call trichk(i+1,j+1,i,j+1,i+1,j,maxl,maxm,s,sval,xp,yp)
        end if

      end if

    end do
  end do

  return
end
subroutine trilhk(i1,j1,i2,j2,i3,j3,jcmax,jcmin,maxl,maxm, &
       ncon,s,smax,smin,xp,yp)

!*****************************************************************************80
!
!! TRILHK is given a triangle, formed by nodes (I1,J1), (I2,J2), and
!  (I3,J3), and the value of the quantity S at these nodes.
!
!  Discussion:
!
!    That portion of the triangle which is greater than contour value
!    value SVAL(I) but less than SVAL(I+1) is to be filled in with
!    color I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1, J1, I2, J2, I3, J3, the I and J indices
!    of the three nodes which form the triangle.
!
!  JCMAX,
!  JCMIN  Input, integer ( kind = 4 ) JCMAX, JCMIN, the maximum and 
!         minimum color indices to use in the color bar.
!
!  MAXL,
!  MAXM   Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!         rows and columns of data.
!
!  NCON   Input, integer ( kind = 4 ) NCON, the number of color contour
!         regions drawn, and hence, the number of colors
!         to be displayed in the color bar.
!
!  S      Input real S(MAXL,MAXM), the value of S at each node.
!
!  SMAX,
!  SMIN   Input, real SMAX, SMIN, the maximum and minimum
!         values of the quantity whose color contours are
!         being drawn.  These numbers will be printed along
!         with the color bar.
!
!  X      Input, real X(MAXL,MAXM).
!         The X coordinates of the nodes.
! 
!  Y      Input, real Y(MAXL,MAXM).
!         The Y coordinates of the nodes.
! 
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ih
  integer ( kind = 4 ) il
  integer ( kind = 4 ) im
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) jcmax
  integer ( kind = 4 ) jcmin
  integer ( kind = 4 ) jcolor
  integer ( kind = 4 ) jh
  integer ( kind = 4 ) jl
  integer ( kind = 4 ) jm
  integer ( kind = 4 ) ncon
  integer ( kind = 4 ) npts
  real px
  real pxold
  real py
  real pyold
  real qx
  real qxold
  real qy
  real qyold
  real s(maxl,maxm)
  real s1
  real s2
  real s3
  real sc1
  real sc2
  real sh
  real sl
  real sm
  real smax
  real smin
  real xp(maxl,maxm)
  real xh
  real xl
  real xm
  real xpoly(5)
  real yp(maxl,maxm)
  real yh
  real yl
  real ym
  real ypoly(5)

  s1 = s(i1,j1)
  s2 = s(i2,j2)
  s3 = s(i3,j3)

  if ( s1 <= s2 .and. s2 <= s3 ) then
    il = i1
    jl = j1
    im = i2
    jm = j2
    ih = i3
    jh = j3
  else if ( s1 <= s3 .and. s3 <= s2 ) then
    il = i1
    jl = j1
    im = i3
    jm = j3
    ih = i2
    jh = j2
  else if ( s2 <= s1 .and. s1 <= s3 ) then
    il = i2
    jl = j2
    im = i1
    jm = j1
    ih = i3
    jh = j3
  else if ( s2 <= s3 .and. s3 <= s1 ) then
    il = i2
    jl = j2
    im = i3
    jm = j3
    ih = i1
    jh = j1
  else if ( s3 <= s1 .and. s1 <= s2 ) then
    il = i3
    jl = j3
    im = i1
    jm = j1
    ih = i2
    jh = j2
  else if ( s3 <= s2 .and. s2 <= s1 ) then
    il = i3
    jl = j3
    im = i2
    jm = j2
    ih = i1
    jh = j1
  end if

  sl = s(il,jl)
  sm = s(im,jm)
  sh = s(ih,jh)

  xl = xp(il,jl)
  xm = xp(im,jm)
  xh = xp(ih,jh)

  yl = yp(il,jl)
  ym = yp(im,jm)
  yh = yp(ih,jh)

  do i = 0,ncon

    sc1 = ((ncon+1-i)*smin+i*smax)/real(ncon+1)
    sc2 = ((ncon-i)*smin+(i+1)*smax)/real(ncon+1)
!
!  Check that some data in the triangle lies in the range [SC1,SC2).
!
!  The following comparison was corrected from ".LT." which would
!  fail on the case where SL = SH.
!
    if ( max(sl,sc1) <= min(sh,sc2) ) then

      jcolor = int(((ncon-i)*jcmin+i*jcmax)/real(ncon))
      call filclr(jcolor)
!
!  Take care of possibility that entire triangle lies in the contour.
!
      if ( sc1 <= sl .and. sh < sc2 ) then

        npts = 3
        xpoly(1) = xl
        ypoly(1) = yl
        xpoly(2) = xm
        ypoly(2) = ym
        xpoly(3) = xh
        ypoly(3) = yh

        call plygon(npts,xpoly,ypoly)
!
!  Find (PXOLD,PYOLD) and (QXOLD,QYOLD), where the line S = SC1 crosses 
!  the triangle.
!
      else

        call cross(px,py,qx,qy,sl,sm,sh,sc1,xl,xm,xh,yl,ym,yh)

        pxold = px
        pyold = py
        qxold = qx
        qyold = qy
!
!  Find (PX,PY) and (QX,QY), where the line S = SC2 crosses the triangle.
!
        call cross(px,py,qx,qy,sl,sm,sh,sc2,xl,xm,xh,yl,ym,yh)
!
!  Now draw the polygon formed by these four points, plus possibly
!  the point (XM,YM).
!
        npts = 4
        xpoly(1) = pxold
        ypoly(1) = pyold
        xpoly(2) = qxold
        ypoly(2) = qyold
        xpoly(3) = qx
        ypoly(3) = qy
        xpoly(4) = px
        ypoly(4) = py

        if ( sc1 <= sm .and. sm <= sc2 ) then
          npts = 5
          xpoly(5) = xm
          ypoly(5) = ym
        end if

        call plygon(npts,xpoly,ypoly)

      end if

    end if

  end do
 
  return
end
subroutine vector(l,m,maxl,maxm,nflag,vecscl,u,v,xp,yp)

!*****************************************************************************80
!
!! VECTOR draws a vector field.
!
!  Discussion:
!
!    An arrow pointing in the direction (U(I), V(I)) is drawn at the 
!    point (X(I),Y(I)).  The arrow's length is scaled by VECSCL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, the number of rows and columns of data.
!
!  MAXL,
!  MAXM   Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!         rows and columns of data.
!
!  NFLAG  Input, logical NFLAG(MAXL,MAXM).
!         NFLAG is used to "flag" which nodes are to be visible.
!
!  VECSCL Input, real VECSCL.
!         A scale factor for velocity vectors.  This starts out at
!         1.0.
! 
!  U      Input, real U(MAXL,MAXM).
! 
!         U(I,J) is the horizontal fluid velocity at node (I,J).
!   
!  V      Input, real V(MAXL,MAXM).
! 
!         V(I,J) is the vertical fluid velocity at node (I,J).
! 
!  X      Input, real XP(MAXL,MAXM).
!         The X coordinates of the nodes.
! 
!  Y      Input, real YP(MAXL,MAXM).
!         The Y coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical nflag(maxl,maxm)
  real u(maxl,maxm)
  real v(maxl,maxm)
  real vecscl
  real xp(maxl,maxm)
  real xtip
  real yp(maxl,maxm)
  real ytip
!
  do i = 1,l
    do j = 1,m

      if ( nflag(i,j) ) then

        xtip = xp(i,j)+vecscl*u(i,j)
        ytip = yp(i,j)+vecscl*v(i,j)

        call arrow(xp(i,j),yp(i,j),xtip,ytip)

      end if

    end do
  end do

  return
end
subroutine vizpn(echo,l,m,maxl,maxm,vpflag,xp,yp)

!*****************************************************************************80
!
!! VIZPN sets the visibility of primary nodes.
!
!  Discussion:
!
!    The user sets a minimum distance.  The routine picks the nodes
!    to display based on the requirement that they be at least this
!    distance apart.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
!    Input, real XP(MAXL,MAXM), the X coordinates of the nodes.
! 
!    Input, real YP(MAXL,MAXM), the Y coordinates of the nodes.
! 
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  real dismin
  real disnod
  logical echo
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) noff
  real temp
  logical vpflag(maxl,maxm)
  real xp(maxl,maxm)
  real yp(maxl,maxm)

  write ( *, * ) ' '
  write ( *, * ) 'VIZPN - Input, please!'
  write ( *, * ) '  Give a minimum separation for visible primary nodes.'

  read(*,*,end = 20,err=20)dismin
  write(17,*)dismin
  if ( echo ) then
    write ( *, * ) dismin
  end if

  temp = dismin**2

  noff = 0
  do i = 1,l
    do j = 1,m

      vpflag(i,j) = .true.

      do ii = 1,l
        do jj = 1,j-1

          if ( vpflag(ii,jj) ) then
            disnod = (xp(i,j)-xp(ii,jj))**2+(yp(i,j)-yp(ii,jj))**2
            if ( disnod < temp ) then
              vpflag(i,j) = .false.
              noff = noff+1
              go to 10
            end if
          end if

        end do
      end do

10        continue

    end do
  end do

  write ( *, * ) ' '
  write ( *, * ) 'VIZPN - Note:'
  write ( *, * ) '  Out of ',l*m,' nodes,'
  write ( *, * ) '  you have turned off ',noff
  write ( *, * ) '  leaving ',l*m-noff

  return

20    continue

  write ( *, * ) ' '
  write ( *, * ) 'VIZPN - Warning'
  write ( *, * ) '  Could not read your value of DISMIN.'
  write ( *, * ) '  Aborting this command!'
  return
end
subroutine vsize ( l, m, maxl, maxm, nflag, u, v, vmax )

!*****************************************************************************80
!
!! VSIZE computes the maximum visible velocity magnitude.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) MAXL, MAXM, the maximum allowed number of
!    rows and columns of data.
!
!    Input, logical NFLAG(MAXL,MAXM).
!    NFLAG(I,J) is .TRUE. if the data at node (I,J) may be displayed,
!    and .FALSE. if it should not be displayed.
!
!    Input, real U(MAXL,MAXM), the horizontal fluid velocity
!    at primary node (I,J).
!   
!    Input, real V(MAXL,MAXM), the vertical fluid velocity 
!    at primary node (I,J).
! 
!    Output, real VMAX.
!    VMAX is the largest Euclidean norm of the displayed
!    velocity vectors.
!
  implicit none

  integer ( kind = 4 ) maxl
  integer ( kind = 4 ) maxm

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical nflag(maxl,maxm)
  real u(maxl,maxm)
  real v(maxl,maxm)
  real vmax

  vmax = 0.0E+00

  do i = 1,l
    do j = 1,m
      if ( nflag(i,j) ) then
        vmax = max ( vmax,sqrt(u(i,j)**2+v(i,j)**2) )
      end if
    end do
  end do

  return
end
