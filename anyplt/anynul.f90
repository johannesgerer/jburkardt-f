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
!    21 November 2000
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
  integer ( kind = 4 ) icom
  integer ( kind = 4 ) iplt1
  integer ( kind = 4 ) iplt2
  integer ( kind = 4 ) ixplt1
  integer ( kind = 4 ) ixplt2
  integer ( kind = 4 ) iyplt1
  integer ( kind = 4 ) iyplt2
  integer ( kind = 4 ) marray
  real xplt1
  real xplt2
  real yplt1
  real yplt2

  common /anycom/ iplt1, iplt2, ixplt1, ixplt2, iyplt1, &
                  iyplt2, marray, xplt1, xplt2, yplt1, yplt2
  common /anychr/ carray

  if ( icom == 13 ) then
    carray = 'ANYNUL - Version 1.02  21 November 2000  Null Graphics'
  end if

  return
end
