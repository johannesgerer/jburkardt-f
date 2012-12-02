subroutine anyplt ( icom )

!*****************************************************************************80
!
!! ANYPLT is an interface routine to a variety of graphics packages.
!
!  Discussion:
!
!    ANYPLT is a subroutine which provides a simple, standard interface
!    between FORTRAN programs and various output devices.  To run a
!    program which calls ANYPLT on a different machine, the program
!    is not modified in any way, but a different version of the ANYPLT
!    program is provided.  Currently, the following versions are available:
!
!    ANYATT - AT&T PC6300 graphics (640 by 400).  Requires ATTPLT.ASM.
!    ANYBUG - Simple debugging output to a file.
!    ANYCAL - CALCOMP file output.  Available on many mainframes.
!    ANYIBM - IBM PC hi resolution (640 by 200).  Requires IBMPLT.ASM.
!    ANYMAC - Macintosh graphics.  Requires auxilliary routine TOOLBX.SUB.
!    ANYNCR - NCAR graphics package.
!    ANYNUL - Does nothing.
!    ANYP10 - PLOT10 interactive graphics. (1024 by 768)
!    ANYTTY - Simple 'typewriter' graphics (80 by 24)
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
!    Input, integer ( kind = 4 ) ICOM, specifies the graphics request being made.
!    0, enable graphics.
!    1, disable graphics.
!    2, begin plot.
!    3, define plot size.
!    4, move to a point.
!    5, draw to a point.
!    6, clear screen.
!    7, write string at position.
!    8, use virtual cursor.
!    9, end plot.
!    10, ring bell.
!    11, mark data.
!    12, return screen data.
!    13, return version.
!    14, draw an arrow.
!
  implicit none

  integer ( kind = 4 )*4 lineto,moveto
  parameter (lineto = z'89109000')
  parameter (moveto = z'89309000')

  save ifont
  save ipoint
  save ixmn,ixmx,iymn,iymx
  save ixmin,ixmax,iymin,iymax
  save rotate
  save xmin,xmax,ymin,ymax
  save xold
  save xsmin,xsmax,ysmin,ysmax
  save xsmn,xsmx,ysmn,ysmx
  save xstart,ystart
  save yold

  character carray*80
  character cong*1
  integer ( kind = 4 )*4 grafptr
  integer ( kind = 4 )   ifont(1617)
  integer ( kind = 4 )   ipoint(95)
  logical   rotate

  common /anycom/ iplt1,iplt2,ixplt1,ixplt2,iyplt1, &
                  iyplt2,marray,xplt1,xplt2,yplt1,yplt2
  common /anychr/ carray
!
!  Pointer array into IFONT
!
  data (ipoint(i),i = 1,95) /
     &   1,   3,  26,  45,  66, 102, 130, 156, 166, 186, 206, 222, 233,
     & 249, 255, 267, 273, 293, 306, 328, 353, 363, 383, 411, 423, 457,
     & 483, 506, 533, 541, 552, 560, 587, 625, 638, 665, 683, 699, 714,
     & 727, 754, 770, 786, 805, 818, 826, 838, 848, 868, 884, 909, 930,
     & 956, 967, 981, 989,1001,1012,1025,1035,1045,1051,1061,1069,1075,
     &1081,1108,1131,1149,1172,1194,1214,1243,1260,1284,1307,1323,1336,
     &1364,1381,1401,1424,1447,1464,1486,1499,1516,1524,1536,1547,1560,
     &1570,1586,1592,1608/
!
!  IFONT contains the strokes defining the various symbols.
!
  data (ifont(i),i =    1, 396)/
     & 1, 0, 2,10,11, 9,22,10,23,11,22,10,11, 0, 9, 7, 9, 9,11, 9,11, 7,
     & 9, 7, 0, 2, 8,17, 7,23, 9,23, 8,17, 0,14,17,13,23,15,23,14,17, 0,
     & 4, 9,23, 7, 7, 0,13,23,11, 7, 0, 5,17,15,17, 0, 5,13,15,13, 0, 3,
     &15,19,13,21, 9,21, 7,19, 7,17, 9,15,13,15,15,13,15,11,13, 9, 9, 9,
     & 7,11, 0, 9,23, 9, 7, 0,13,23,13, 7, 0, 3, 5,23, 9,23, 9,19, 5,19,
     & 5,23, 0,17,23, 5, 7, 0,13, 7,13,11,17,11,17, 7,13, 7, 0, 1,17, 7,
     & 7,17, 7,19, 9,21,13,21,15,19,15,17, 5,13, 5,11, 9, 7,13, 7,17,15,
     & 0, 1,10,17, 9,23,11,23,10,17, 0, 1,12,23,11,21,10,19, 9,17, 9,15,
     & 9,13,10,11,11, 9,12, 7, 0, 1,12,23,13,21,14,19,15,17,15,15,15,13,
     &14,11,13, 9,12, 7, 0, 3, 7,15,15,15, 0,13,19, 9,11, 0, 9,19,13,11,
     & 0, 2, 7,15,15,15, 0,11,19,11,11, 0, 1,11, 7, 9, 7, 9, 9,11, 9,11,
     & 7,11, 6,10, 4, 0, 1, 7,15,15,15, 0, 1, 9, 7, 9, 9,11, 9,11, 7, 9,
     & 7, 0, 1,15,23, 7, 7, 0, 1, 9,23,13,23,15,19,15,11,13, 7, 9, 7, 7,
     &11, 7,19, 9,23, 0, 2, 7,21, 9,23, 9, 7, 0, 7, 7,11, 7, 0, 1, 5,21,
     & 9,23,15,23,17,21,17,19,15,17, 7,13, 5,10, 5, 7,17, 7, 0, 2, 5,23,
     &17,23,15,17,13,15, 9,15, 0,13,15,17,13,17,10,14, 7, 8, 7, 5,10, 0,
     & 1,13, 7,13,23, 5,13,17,13, 0, 1,17,23, 5,23, 5,17,13,17,17,15,17,
     &11,13, 7, 9, 7, 5,11, 0, 1,17,19,13,23, 9,23, 5,19, 5,13, 9,15,13/
  data (ifont(i),i =  397, 792)/
     &15,17,13,17,11,13, 7, 9, 7, 5,11, 5,13, 0, 1, 5,19, 5,23,17,23,11,
     &15,11, 7, 0, 1, 8,15, 6,17, 6,21, 8,23,14,23,16,21,16,17,14,15, 8,
     &15, 5,13, 5, 9, 8, 7,14, 7,17, 9,17,13,14,15, 0, 1,17,17,15,15, 7,
     &15, 5,17, 5,21, 7,23,15,23,17,21,17,11,15, 7, 7, 7, 5,11, 0, 2, 9,
     &13, 9,15,11,15,11,13, 9,13, 0, 9, 7, 9, 9,11, 9,11, 7, 9, 7, 0, 2,
     & 9,13, 9,15,11,15,11,13, 9,13, 0,11, 7, 9, 7, 9, 9,11, 9,11, 7,11,
     & 6,10, 4, 0, 1,17,21, 5,15,17, 9, 0, 2, 7,15,15,15, 0, 7, 9,15, 9,
     & 0, 1, 5,21,17,15, 5, 9, 0, 2, 7,21, 9,23,13,23,15,21,15,19,11,15,
     &11,11, 0,10, 7,10, 9,12, 9,12, 7,10, 7, 0, 1,13, 7, 9, 7, 5,11, 5,
     &19, 9,23,13,23,17,19,17,11,15, 9,13,11,12,10,10,10, 9,11, 9,15,10,
     &16,12,16,13,15,13,11, 0, 2, 5, 7,11,23,17, 7, 0, 8,15,14,15, 0, 2,
     & 5, 7, 5,23,15,23,17,21,17,17,15,15, 5,15, 0,15,15,17,13,17, 9,15,
     & 7, 5, 7, 0, 1,17,19,13,23, 9,23, 5,19, 5,11, 9, 7,13, 7,17,11, 0,
     & 1, 5, 7, 5,23,13,23,17,19,17,11,13, 7, 5, 7, 0, 2,17,23, 5,23, 5,
     & 7,17, 7, 0, 5,15,12,15, 0, 2, 5, 7, 5,23,17,23, 0, 5,15,12,15, 0,
     & 2,17,19,13,23, 9,23, 5,19, 5,11, 9, 7,13, 7,17,11,17,15,13,15, 0,
     &17,11,17, 7, 0, 3, 5, 7, 5,23, 0, 5,15,17,15, 0,17,23,17, 7, 0, 3,
     & 9,23,13,23, 0,11,23,11, 7, 0, 9, 7,13, 7, 0, 2,15,23,15,11,12, 7/
  data (ifont(i),i =  793,1188)/
     & 8, 7, 5,11, 5,13, 0,13,23,17,23, 0, 2, 5, 7, 5,23, 0,17,23, 5,15,
     &17, 7, 0, 1, 5,23, 5, 7,17, 7, 0, 1, 5, 7, 5,23,11,11,17,23,17, 7,
     & 0, 1, 5, 7, 5,23,17, 7,17,23, 0, 1,17,19,13,23, 9,23, 5,19, 5,11,
     & 9, 7,13, 7,17,11,17,19, 0, 1, 5, 7, 5,23,13,23,17,21,17,17,13,15,
     & 5,15, 0, 2,17,19,13,23, 9,23, 5,19, 5,11, 9, 7,13, 7,17,11,17,19,
     & 0,13,11,17, 7, 0, 2, 5, 7, 5,23,13,23,17,21,17,17,13,15, 5,15, 0,
     &13,15,17, 7, 0, 1,17,19,13,23, 9,23, 5,20, 5,18, 9,15,13,15,17,12,
     &17,10,13, 7, 9, 7, 5,10, 0, 2, 5,23,17,23, 0,11,23,11, 7, 0, 1, 5,
     &23, 5,10, 8, 7,14, 7,17,10,17,23, 0, 1, 5,23,11, 7,17,23, 0, 1, 5,
     &23, 8, 7,11,17,14, 7,17,23, 0, 2, 5,23,17, 7, 0,17,23, 5, 7, 0, 2,
     & 5,23,11,13,17,23, 0,11,13,11, 7, 0, 1, 5,23,17,23, 5, 7,17, 7, 0,
     & 1,11,23, 7,23, 7, 7,11, 7, 0, 1, 7,23,15, 7, 0, 1, 7,23,11,23,11,
     & 7, 7, 7, 0, 1, 7,21,11,23,15,21, 0, 1, 5, 3,17, 3, 0, 1, 9,23,13,
     &19, 0, 2, 7,14, 9,15,13,15,15,14,15, 7, 0,15,12, 9,12, 7,11, 7, 8,
     & 9, 7,13, 7,15, 8, 0, 2, 7,23, 7, 7, 0, 7,13, 9,15,13,15,15,13,15,
     & 9,13, 7, 9, 7, 7, 9, 0, 1,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7,13,
     & 7,15, 9, 0, 2,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7,13, 7,15, 9, 0,
     &15,23,15, 7, 0, 1, 7,11,15,11,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7/
  data (ifont(i),i = 1189,1584)/
     &13, 7,15, 9, 0, 3, 9, 7, 9,23,13,23,13,22, 0, 8,15,12,15, 0, 8, 7,
     &11, 7, 0, 2,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7,13, 7,15, 9, 0,15,
     &13,15, 3,13, 1, 9, 1, 7, 3, 0, 2, 7, 7, 7,23, 0, 7,14, 9,15,13,15,
     &15,14,15, 7, 0, 3, 9,15,11,15,11, 7, 0, 9, 7,13, 7, 0, 9,17, 9,19,
     &11,19,11,17, 9,17, 0, 2, 9,15,11,15,11, 1, 7, 1, 7, 3, 0, 9,17,11,
     &17,11,19, 9,19, 9,17, 0, 3, 7, 7, 7,23, 0,15,15, 7,10, 0, 9,11,15,
     & 7, 0, 2, 9,23,11,23,11, 7, 0, 9, 7,13, 7, 0, 3, 7,15, 7, 7, 0, 7,
     &14, 8,15,10,15,11,14,11, 7, 0,11,14,12,15,14,15,15,14,15, 7, 0, 2,
     & 7, 7, 7,15, 0, 7,14, 9,15,13,15,15,14,15, 7, 0, 1, 7,13, 9,15,13,
     &15,15,13,15, 9,13, 7, 9, 7, 7, 9, 7,13, 0, 2, 7,13, 9,15,13,15,15,
     &13,15, 9,13, 7, 9, 7, 7, 9, 0, 7,14, 7, 1, 0, 2,15,13,13,15, 9,15,
     & 7,13, 7, 9, 9, 7,13, 7,15, 9, 0,15,14,15, 1, 0, 2, 7,15, 9,15, 9,
     & 7, 0, 9,13,11,15,13,15,15,13, 0, 1,15,13,13,15, 9,15, 7,13, 9,11,
     &13,11,15, 9,13, 7, 9, 7, 7, 9, 0, 2, 9,23, 9, 7,11, 7, 0, 7,17,11,
     &17, 0, 2, 7,15, 7, 9, 9, 7,13, 7,15, 9, 0,15,15,15, 7, 0, 1, 7,15,
     &11, 7,15,15, 0, 1, 7,15, 9, 7,11,11,13, 7,15,15, 0, 2, 7,15,15, 7,
     & 0, 7, 7,15,15, 0, 2, 7,15,11, 7, 0,15,15,10, 5, 7, 1, 0, 1, 7,15,
     &15,15, 7, 7,15, 7, 0, 1,11,23, 7,23, 9,17, 7,15, 9,13, 7, 7,11, 7/
  data (ifont(i),i = 1585,1617)/
     & 0, 1, 9,23, 9, 7, 0, 1, 7,23,11,23, 9,17,11,15, 9,13,11, 7, 7, 7,
     & 0, 1, 5,21, 7,23,15,21,17,23, 0/
!
!  ICOM = 0  Enable graphics.  For interactive graphics, clear the screen.
!  For graphics which create a file, open the file.
!  Also, user will input the portion of the (0,1) by (0,1) output
!  screen that is to be used.  At this point, the effective screen
!  (IXMIN,IXMAX) by (IYMIN,IYMAX) must be computed.
!
!  Note that the whole Macintosh screen is NOT available to Microsoft FORTRAN.
!  Note also the inversion that causes YSMN to be greater than YSMX!
!
  if ( icom == 0 ) then
    xsmn = 116
    xsmx = 396
    ysmn = 311
    ysmx = 31
    ixmn = xsmn
    ixmx = xsmx
    iymn = ysmn
    iymx = ysmx
    if ( xplt1 == xplt2 ) then
      write(*,*)'anymac - warning!  icom = 0 command is asking for'
      write(*,*)'         zero screen width (xplt1 = xplt2).'
      xplt1 = 0.0
      xplt2 = 1.0
    end if
    xsmin = xsmn+xplt1*(xsmx-xsmn)
    xsmax = xsmn+xplt2*(xsmx-xsmn)
    if ( yplt1 == yplt2 ) then
      write(*,*)'anymac - warning!  icom = 0 command is asking for'
      write(*,*)'         zero screen height (yplt1 = yplt2).'
      yplt1 = 0.0
      yplt2 = 1.0
    end if
    ysmin = ysmn+yplt1*(ysmx-ysmn)
    ysmax = ysmn+yplt2*(ysmx-ysmn)
    ixmin = xsmin
    ixmax = xsmax
    iymin = ysmin
    iymax = ysmax
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
!
!  ICOM = 1  Disable graphics.  For interactive graphics, clear the screen.
!  For graphics which create a file, close the file.
!  (Doesn't really apply for Macintosh)
!
!  ICOM = 2  Begin plot
!
  else if ( icom == 2 ) then
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
    do i = 1,48
      write(*,'(1x)')
    end do
!
!  ICOM = 3  Define user coordinate system (XMIN,XMAX), (YMIN,YMAX).
!
  else if ( icom == 3 ) then
    xmin = xplt1
    xmax = xmin+xplt2
    ymin = yplt1
    ymax = ymin+yplt2
!
!  ICOM = 4  Move to point.  Only purpose is to begin a line.
!  Input is X, Y in user coordinates.
!
  else if ( icom == 4 ) then
    xold = xplt1
    yold = yplt1
!
!  ICOM = 5  Draw to point.  Only called after a previous move or draw.
!  Input is X, Y in user coordinates.
!
  else if ( icom == 5 ) then
    call drwlin(xold,yold,xplt1,yplt1,xmin,xmax,ymin,ymax,
     &  ixmin,ixmax,iymin,iymax)
    xold = xplt1
    yold = yplt1
!
!  ICOM = 6  Clear screen.  Must be a better way on the Macintosh, but I
!  don't know it yet.
!
  else if ( icom == 6 ) then
    do i = 1,48
      write(*,'(1x)')
    end do
!
!  ICOM = 7,  Write string at position.
!  Variable height and angle should be supported.
!  Note that for this command, screen coordinates are used.
!  Thus a width of 0.1 is intended to mean 1/10 of screen size.
!
  else if ( icom == 7 ) then
!
!  Set scale factor for character height
!
    csize = xplt2
    angle = yplt2
    scl2 = csize*.0625
!
!  Set starting point for line of text (lower left corner of first
!  character) and origin about which rotation is performed.
!
    xb = xplt1
    xrot = xplt1
    yb = yplt1
    yrot = yplt1
    rotate = .false.
!
!  Get trig functions if rotation required, converting from
!  degrees to radians.
!
    if ( angle /= 0.0 ) then
      ca = cos(angle*.017453)
      sa = sin(angle*.017453)
      rotate = .true.
    end if
!
!  Loop over all characters in string
!
    do 30 icr = 1,marray
!
!  Get ASCII code for character and shift by 31 so first printable
!  character becomes code 1
!
      iascii = ichar(carray(icr:icr))-31
!
!  Replace any non-printable characters with blanks
!
      if ( iascii < 1 .or. isscii > 95 )iascii = 1
!
!  Get pointer to this character in font table
!
      ip = ipoint(iascii)
!
!  Get number of "vectors" required to draw character.
!  Here "vectors" means number of times pen is lowered, not the
!  number of pen strokes. ( = 1 for blanks, due to way algorithm is coded)
!
      nvec = ifont(ip)
!
!  Loop over all required pen movements
!
      do 20 iv = 1,nvec

        ipen = 3
        ip = ip+1
15      continue
        if ( ifont(ip)==0)go to 20
        x = xb+(ifont(ip)-1)*scl2
        y = yb+(ifont(ip+1)-7)*scl2
!
!  Apply rotation if necessary
!
        if ( rotate ) then
          xt = x-xrot
          yt = y-yrot
          x = ca*xt-sa*yt+xrot
          y = sa*xt+ca*yt+yrot
          endif
!
!  Plot the pen stroke
!
       if ( ipen==3 ) then
          call xtoix(ix2,ixmax,ixmin,x,1.0,0.0)
          call xtoix(iy2,iymax,iymin,y,1.0,0.0)
          ix = ix2
          iy = iy2
          call toolbx ( moveto, ix, iy )
        else
          ix1 = ix2
          iy1 = iy2
          call xtoix(ix2,ixmax,ixmin,x,1.0,0.0)
          call xtoix(iy2,iymax,iymin,y,1.0,0.0)
          ix = ix2
          iy = iy2
          call toolbx ( lineto, ix, iy )
          endif
        ipen = 2
        ip = ip+2
        go to 15
20          continue
!
!  Advance base to compensate for character just drawn
!
      xb = xb+csize
30        continue
!
!  ICOM = 8  Use virtual cursor.  Not implemented.
!
  else if ( icom == 8 ) then
!
!  ICOM = 9  End plot
!
  else if ( icom == 9 ) then
    write(*,'(1x)')
    read(*,'(1x)')
    do i  =  1, 48
      write(*,'(1x)')
    end do
!
!  ICOM = 10  Ring bell
!
  else if ( icom == 10 ) then
    cong = char(7)
    write(*,'(1x,a1)')cong
!
!  ICOM = 11  Mark data with a set of strokes that are like
!  a *, a +, or an X.  If a '.' is requested, actually try to
!  draw a point (pixel) if possible.
!
  else if ( icom == 11 ) then
    if ( carray(1:1) == ' ' ) then
      return
    else if ( carray(1:1) == '.' ) then
      call xtoix(ix1,ixmax,ixmin,xplt1,xmax,xmin)
      call xtoix(iy1,iymax,iymin,yplt1,ymax,ymin)
      ix = ix1
      iy = iy1
      call toolbx ( moveto, ix, iy )
      call toolbx ( lineto, ix, iy )
    else
      delt = 0.5*xlen*14.0/real(ixlen)
      x1 = xplt1+delt
      x2 = xplt1-delt
      y1 = yplt1+delt
      y2 = yplt1-delt
      call xtoix(ix1,ixmax,ixmin,x1,xmax,xmin)
      call xtoix(ix2,ixmax,ixmin,x2,xmax,xmin)
      call xtoix(ix3,ixmax,ixmin,x3,xmax,xmin)
      call xtoix(iy1,iymax,iymin,y1,ymax,ymin)
      call xtoix(iy2,iymax,iymin,y2,ymax,ymin)
      call xtoix(iy3,iymax,iymin,y3,ymax,ymin)

      if ( carray(1:1) /= '+' ) then
        ix = ix1
        iy = iy2
        call toolbx ( moveto, ix, iy )
        ix = ix2
        iy = iy1
        call toolbx ( lineto, ix, iy )
        ix = ix1
        iy = iy1
        call toolbx ( moveto, ix, iy )
        ix = ix2
        iy = iy2
        call toolbx ( lineto, ix, iy )
      end if

      if ( carray(1:1) /=  'X' .and. carray(1:1) .ne. 'x' ) then
        ix = ix1
        iy = iy3
        call toolbx ( moveto, ix, iy )
        ix = ix2
        iy = iy3
        call toolbx ( lineto, ix, iy )
        ix = ix3
        iy = iy1
        call toolbx ( moveto, ix, iy )
        ix = ix3
        iy = iy2
        call toolbx ( lineto, ix, iy )
        ix2 = ix3
      end if

    end if
!
!  ICOM = 12, Return screen X, Y maximum and minimum in pixels
!  or other 'natural' coordinates.
!
  else if ( icom == 12 ) then
    xplt1 = xsmn
    xplt2 = xsmx
    yplt1 = ysmn
    yplt2 = ysmx
!
!  ICOM = 13  Return version number and date of code.
!
  else if ( icom == 13 ) then
    carray = 'anyplt - version 1.08  09 october 1990  macintosh'
!
!  ICOM = 14, draw an arrow of given screen length and screen angle.
!  Even if the X and Y axes are of vastly different scales, this
!  allows one to draw vectors of a predictable direction.
!
  else if ( icom == 14 ) then
    x1 = xplt1
    y1 = yplt1
    angle = xplt2
    alen = yplt2
    x2 = x1+(xmax-xmin)*alen*cos(angle)
    y2 = y1+(ymax-ymin)*alen*sin(angle)
    call drwlin(x1,y1,x2,y2,xmin,xmax,ymin,ymax,
     &  ixmin,ixmax,iymin,iymax)
!
!  Unknown command.
!
  else
    write ( *, *) ' '
    write ( *, * ) 'ANYMAC - Fatal error!'
    write ( *, * ) '  Unknown command  =  ', icom
    stop
  end if

  return
end
subroutine drwlin ( xold, yold, xplt1, yplt1, xmin, xmax, ymin, ymax, &
  ixmin, ixmax, iymin, iymax )

!*****************************************************************************80
!
!! DRWLIN
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 )*4 ix
  integer ( kind = 4 )*4 iy
  integer ( kind = 4 )*4 lineto,moveto

  parameter (lineto = z'89109000')
  parameter (moveto = z'89309000')

  call clip ( xold, yold, xplt1, yplt1, xc, yc, xd, yd, idraw, xmin, &
    xmax, ymin, ymax )

  if ( idraw == 1 ) then
    call xtoix ( ix1, ixmax, ixmin, xc, xmax, xmin )
    call xtoix ( iy1, iymax, iymin, yc, ymax, ymin )
    call xtoix ( ix2, ixmax, ixmin, xd, xmax, xmin )
    call xtoix ( iy2, iymax, iymin, yd, ymax, ymin )

    ix = ix1
    iy = iy1
    call toolbx ( moveto, ix, iy )
    ix = ix2
    iy = iy2
    call toolbx ( lineto, ix, iy )
  end if
 
  return
end
subroutine xtoix(ix,ixmax,ixmin,x,xmax,xmin)

!*****************************************************************************80
!
!! XTOIX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
  implicit none

  ix = int(real(ixmax-ixmin)*(x-xmin)/(xmax-xmin))+ixmin
  if ( ix < ixmin .and. ix < ixmax)ix = min(ixmin,ixmax)
  if ( ix > ixmax .and. ix > ixmax)ix = max(ixmin,ixmax)

  return
end
subroutine clip(xa,ya,xb,yb,xc,yc,xd,yd,idraw,x0,x1,y0,y1)

!*****************************************************************************80
!
!! CLIP
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
  implicit none

  real xval(2)
  real yval(2)

  call clip1 ( xa, ya, xb, yb, xval, yval, next, x0, x1, y0, y1 )

  if ( next < 2 ) then
    idraw = 0
  else
    idraw = 1
  end if

  xc = xval(1)
  yc = yval(1)
  xd = xval(2)
  yd = yval(2)

  return
end
subroutine clip1(xa,ya,xb,yb,xval,yval,next,x0,x1,y0,y1)

!*****************************************************************************80
!
!! CLIP1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
  implicit none

  dimension xval(2),yval(2)
!
!  Check to see if both points are interior to the box, and hence
!  nothing need be done.
!
  x00 = min(x0,x1)
  x11 = max(x0,x1)
  y00 = min(y0,y1)
  y11 = max(y0,y1)

  next = 0
  if ( (x00 <= xa .and. xa <= x11).and. &
       (y00 <= ya .and. ya <= y11) ) then
    next = next+1
    xval(next) = xa
    yval(next) = ya
  end if

  if ( (x00 <= xb .and. xb <= x11).and. &
       (y00 <= yb .and. yb <= y11) ) then
    next = next+1
    xval(next) = xb
    yval(next) = yb
  end if

  if ( next ==2 ) then
    return
  end if

  call clip2(xa,ya,xb,yb,x00,y00,y11,y,xval,yval,next)

  if ( next ==2 ) then
    return
  end if

  call clip2(xa,ya,xb,yb,x11,y00,y11,y,xval,yval,next)

  if ( next ==2 ) then
    return
  end if

  call clip2(ya,xa,yb,xb,y00,x00,x11,x,yval,xval,next)

  if ( next ==2 ) then
    return
  end if

  call clip2(ya,xa,yb,xb,y11,x00,x11,x,yval,xval,next)

  return
end
subroutine clip2 ( xa, ya, xb, yb, x00, y00, y11, y, xval, yval, next )

!*****************************************************************************80
!
!! CLIP2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
  implicit none

  integer ( kind = 4 ) next
  real t
  real x00
  real xa
  real xb
  real xval(*)
  real y
  real y00
  real y11
  real yval(*)

  if ( xb == xa ) then
    return
  end if

  t = ( x00 - xa ) / ( xb - xa )

  if ( 0.0 <= t .and. t <= 1.0 ) then

    y = ( 1.0 - t ) * ya + t * yb

    if ( y00 <= y .and. y <= y11 ) then
      next = next + 1
      xval(next) = x00
      yval(next) = y
    end if

  end if

  return
end
