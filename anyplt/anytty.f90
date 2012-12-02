subroutine anyplt ( icom )

!*****************************************************************************80
!
!! ANYPLT is a generic graphics interface routine.
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
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ICOM, the index of the graphics command.
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

  character ( len = 80 ) carray
  real, parameter :: degrees_to_radians = 3.14159265 / 180.0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) icom
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) iplt1
  integer ( kind = 4 ) iplt2
  character isay
  integer ( kind = 4 ) isig
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) ixend
  integer ( kind = 4 ), save :: ixmax
  integer ( kind = 4 ), save :: ixmin
  integer ( kind = 4 ), parameter :: ixmn = 1
  integer ( kind = 4 ), parameter :: ixmx = 80
  integer ( kind = 4 ) ixplt1
  integer ( kind = 4 ) ixplt2
  integer ( kind = 4 ) ixstr
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) iyend
  integer ( kind = 4 ), save :: iymax
  integer ( kind = 4 ), save :: iymin
  integer ( kind = 4 ), parameter :: iymn = 24
  integer ( kind = 4 ), parameter :: iymx =  1
  integer ( kind = 4 ) iyplt1
  integer ( kind = 4 ) iyplt2
  integer ( kind = 4 ) iystr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) marray
  character, save, dimension ( iymx:iymn, ixmn:ixmx ) :: screen
  real x
  real x1
  real x2
  real xend
  real, save :: xmax
  real, save :: xmin
  real xplt1
  real xplt2
  real, save :: xsmax
  real, save :: xsmin
  real, save :: xsmn
  real, save :: xsmx
  real, save :: xstart
  real y
  real y1
  real y2
  real yend
  real, save :: ymax
  real, save :: ymin
  real yplt1
  real yplt2
  real, save :: ysmax
  real, save :: ysmin
  real, save :: ysmn
  real, save :: ysmx
  real, save :: ystart

  common /anycom/ iplt1, iplt2, ixplt1, ixplt2, iyplt1, &
                  iyplt2, marray, xplt1, xplt2, yplt1, yplt2

  common /anychr/ carray

  xsmn = 1
  xsmx = 80
  ysmn = 24
  ysmx = 1
!
!  ICOM = 0  Enable graphics
!
  if ( icom == 0 ) then

    xsmin = xsmn + xplt1 * ( xsmx - xsmn )
    xsmax = xsmn + xplt2 * ( xsmx - xsmn )
    ysmin = ysmn + yplt1 * ( ysmx - ysmn )
    ysmax = ysmn + yplt2 * ( ysmx - ysmn )

    ixmin = int ( xsmin )
    ixmax = int ( xsmax )
    iymin = int ( ysmin )
    iymax = int ( ysmax )
!
!  ICOM = 1  Disable graphics
!
  else if ( icom == 1 ) then
!
!  ICOM = 2  Begin plot
!
  else if ( icom == 2 ) then
    screen(iymx:iymn,ixmn:ixmx) = ' '
!
!  ICOM = 3  Define the plot size
!
  else if ( icom == 3 ) then

    xmin = xplt1
    xmax = xmin + xplt2
    ymin = yplt1
    ymax = ymin + yplt2
!
!  ICOM = 4  Move to point
!
  else if ( icom == 4 ) then
    xstart = xplt1
    ystart = yplt1
!
!  ICOM = 5  Draw to point
!
  else if ( icom == 5 ) then

    xend = xplt1
    yend = yplt1

    call r4_to_i4 ( xend,   xmax, xmin, ixend, ixmax, ixmin )
    call r4_to_i4 ( xstart, xmax, xmin, ixstr, ixmax, ixmin )
    call r4_to_i4 ( yend,   ymax, ymin, iyend, iymax, iymin )
    call r4_to_i4 ( ystart, ymax, ymin, iystr, iymax, iymin )

    if ( abs ( ixend - ixstr) > abs ( iyend - iystr ) ) then

      ihi = abs ( ixend - ixstr ) + 1

      if ( ixend < ixstr ) then
        isig = -1
      else
        isig = +1
      end if

      do icount = 1, ihi
        j = ixstr + isig * ( icount - 1 )
        i = iystr + int ( real ( &
          ( iyend - iystr ) * ( j - ixstr ) / ( ixend - ixstr ) ) )
        if ( iymx <= i .and. i <= iymn .and. &
             ixmn <= j .and. j <= ixmx ) then
          if ( screen(i,j) == ' ') then
            screen(i,j) = '*'
          end if
        end if
      end do

    else if ( abs ( iyend - iystr ) > 0 ) then

      jhi = abs ( iyend - iystr ) + 1

      if ( iyend < iystr ) then
        isig = -1
      else
        isig = +1
      end if

      do icount = 1, jhi

        i = iystr + isig * ( icount - 1 )

        j = ixstr + int ( real ( &
          ( ixend - ixstr ) * ( i - iystr ) / ( iyend - iystr ) ) )

        if ( iymx <= i .and. i <= iymn .and. &
             ixmn <= j .and. j <= ixmx ) then
          if ( screen(i,j) == ' ' ) then
            screen(i,j) = '*'
          end if
        end if
      end do

    else

      i = iystr
      j = ixstr
      if ( iymx <= i .and. i <= iymn .and. &
           ixmn <= j .and. j <= ixmx ) then
        if ( screen(i,j) == ' ' ) then
          screen(i,j) = '*'
        end if
      end if

    end if

    xstart = xend
    ystart = yend
!
!  ICOM = 6  Clear screen
!
  else if ( icom == 6 ) then

    screen(iymx:iymn,ixmn:ixmx) = ' '

    do i = iymx, iymn
      write ( *, * ) ' '
    end do
!
!  ICOM = 7,  Write string at position
!
  else if ( icom == 7 ) then

    call r4_to_i4 ( xplt1, xmax, xmin, ix, ixmax, ixmin )
    call r4_to_i4 ( yplt1, ymax, ymin, iy, iymax, iymin )

    do i = 1, marray
      if ( (ix+i-1) <= ixmax ) then
        screen(iy,ix+i-1) = carray(i:i)
      end if
    end do
!
!  ICOM = 8  Use virtual cursor
!
  else if ( icom == 8 ) then
!
!  ICOM = 9  End plot
!
  else if ( icom == 9 ) then

    do i = iymax, iymin
      write ( *, '(80a)' ) screen(i,1:ixmax)
    end do

!   write ( *, * ) ' '
!   write ( *, * ) 'Type RETURN to continue.'
!   read ( *, * )

!
!  ICOM = 10  Ring bell
!
  else if ( icom == 10 ) then
    write ( *, '(a)' ) char ( 7 )
!
!  ICOM = 11  Mark data
!
  else if ( icom == 11 ) then
    call r4_to_i4 ( xplt1, xmax, xmin, j, ixmax, ixmin )
    call r4_to_i4 ( yplt1, ymax, ymin, i, iymax, iymin )
    if ( iymx <= i .and. i <= iymn .and. &
         ixmn <= j .and. j <= ixmx ) then
      screen(i,j) = carray(1:1)
    end if
!
!  ICOM = 12  Return screen data
!
  else if ( icom == 12 ) then
    xplt1 = xsmn
    xplt2 = xsmx
    yplt1 = ysmn
    yplt2 = ysmx
!
!  ICOM = 13  Return version
!
  else if ( icom == 13 ) then
    carray = 'ANYPLT - Version 1.04  10 April 2001  TTY Graphics'
!
!  ICOM = 14, Draw an arrow.
!
  else if ( icom == 14 ) then

    x1 = xplt1
    y1 = yplt1

    call r4_to_i4 ( x1, xmax, xmin, i1, ixmax, ixmin )
    call r4_to_i4 ( y1, ymax, ymin, j1, iymax, iymin )

    x2 = xplt1 + yplt2 * cos ( degrees_to_radians * xplt2 )
    y2 = yplt1 + yplt2 * sin ( degrees_to_radians * xplt2 )

    call r4_to_i4 ( x2, xmax, xmin, i2, ixmax, ixmin )
    call r4_to_i4 ( y2, ymax, ymin, j2, iymax, iymin )

    if ( i1 /= i2 ) then

      if ( i1 > i2 ) then
        call i4_swap ( i1, i2 )
        call i4_swap ( j1, j2 )
      end if

      do i = i1, i2
        y = real ( ( i2 - i ) * j1 + ( i - i1 ) * j2 ) / real ( i2 - i1 )
        j = nint ( y )
        if ( iymx <= i .and. i <= iymn .and. &
             ixmn <= j .and. j <= ixmx ) then
          screen(i,j) = '@'
        end if
      end do

    end if

    if ( j1 /= j2 ) then

      if ( j1 > j2 ) then
        call i4_swap ( i1, i2 )
        call i4_swap ( j1, j2 )
      end if

      do j = j1, j2

        x = real ( ( j2 - j      ) * i1   &
                 + (      j - j1 ) * i2 ) &
          / real   ( j2     - j1 )

        i = nint ( x )

        if ( iymx <= i .and. i <= iymn .and. &
             ixmn <= j .and. j <= ixmx ) then
          screen(i,j) = '@'
        end if

      end do

    end if

    if ( i1 == i2 .and. j1 == j2 ) then
      if ( iymx <= i1 .and. i1 <= iymn .and. &
           ixmn <= j1 .and. j1 <= ixmx ) then
        screen(i1,j1) = '@'
      end if
    end if
!
!  Unknown command.
!
  else
    write ( *, * ) 'ANYPLT - Fatal error!'
    write ( *, * ) '  Unknown value of ICOM = ', icom
    stop
  end if

  return
end
subroutine r4_to_i4 ( x, xmax, xmin, i, imax, imin )

!*****************************************************************************80
!
!! R4_TO_I4 maps real X in [XMIN, XMAX] to integer I in [IMIN, IMAX].
!
!  Discussion:
!
!    I := IMIN + ( IMAX - IMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
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
!    Input, real X, the real number to be converted.
!
!    Input, real XMAX, XMIN, the real range.
!
!    Output, integer ( kind = 4 ) I, the value in the range [IMIN,IMAX] that
!    corresponds to X.
!
!    Input, integer ( kind = 4 ) IMAX, IMIN, the allowed range of the output
!    variable.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real x
  real xmax
  real xmin

  if ( xmax == xmin ) then

    i = ( imax + imin ) / 2

  else

    i = nint ( ( ( xmax - x ) * real ( imin ) + ( x - xmin ) * real ( imax ) ) &
      / ( xmax - xmin ) )

  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
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
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
