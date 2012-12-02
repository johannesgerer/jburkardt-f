program main

!*****************************************************************************80
!
!! MAIN is the main program for HB_TO_ST.
!
!  Discussion:
!
!    HB_TO_ST converts a Harwell-Boeing file to sparse triplet format.
!
!    Read a sparse matrix in the Harwell/Boeing format and output a matrix
!    in zero-based sparse triplet format.  Only the lower triangular part of a
!    symmetric matrix is provided.  Does not handle skew-symmetric
!    matrices.
!
!  Usage:
!
!    hb_to_st hb_file st_file
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Timothy Davis,
!    Direct Methods for Sparse Linear Systems,
!    SIAM, Philadelphia, 2006.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 10000
  integer ( kind = 4 ), parameter :: nzmax = 30000
!
!  The values of NMAX and NZMAX specified here determine the maximum size
!  of the matrices that can be read by the program.  If the program fails
!  to read a matrix, and the error is caused by insufficient space, a message
!  will be printed suggesting the correct larger value of NMAX or NZMAX
!  that should be used here.
!
  integer   ( kind = 4 )  arg_num
  integer   ( kind = 4 )  col
  integer   ( kind = 4 )  iarg
  integer   ( kind = 4 )  iargc
  integer   ( kind = 4 )  indcrd
  integer   ( kind = 4 )  index(nzmax)
  character ( len = 16 )  indfmt
  character ( len = 255 ) input_filename
  integer   ( kind = 4 )  input_unit
  integer   ( kind = 4 )  ios
  character ( len = 30 )  key
  integer   ( kind = 4 )  n
  integer   ( kind = 4 )  ncol
  integer   ( kind = 4 )  nel
  integer   ( kind = 4 )  nrhs
  integer   ( kind = 4 )  nrow
  integer   ( kind = 4 )  nz
  integer   ( kind = 4 )  nzrhs
  character ( len = 255 ) output_filename
  integer   ( kind = 4 )  output_unit
  integer   ( kind = 4 )  p
  integer   ( kind = 4 )  ptr(nmax)
  integer   ( kind = 4 )  ptrcrd
  character ( len = 16 )  ptrfmt
  integer   ( kind = 4 )  rhscrd
  character ( len = 20 )  rhsfmt
  character ( len = 3 )   rhstyp
  integer   ( kind = 4 )  row
  real ( kind = 8 )  skew
  integer   ( kind = 4 )  stype
  logical sym
  character ( len = 72 )  title
  integer   ( kind = 4 )  totcrd
  character ( len = 3 )   type
  integer   ( kind = 4 )  valcrd
  character ( len = 20 )  valfmt
  real ( kind = 8 )  value(nzmax)

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HB_TO_ST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Convert sparse matrix from HB to ST format.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, input_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_TO_ST:'
    write ( *, '(a)' ) '  Please enter the name of the input file.'

    read ( *, '(a)' ) input_filename

  end if
!
!  Second command line argument is the output file name.
!
  if ( 2 <= arg_num ) then

    iarg = 2
    call getarg ( iarg, output_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_TO_ST:'
    write ( *, '(a)' ) '  Please enter the name of the output file.'

    read ( *, '(a)' ) output_filename

  end if
!
!  Open the input file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_TO_ST - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( input_filename ) // '".'
    stop
  end if
!
!  Read header information.
!
  read ( input_unit, '(a72,a8)', iostat = ios ) title, key

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_TO_ST - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading HB header record #1.'
    stop
  end if

  read ( input_unit, '(5i14)', iostat = ios ) totcrd, ptrcrd, indcrd, valcrd, rhscrd

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_TO_ST - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading HB header record #2.'
    stop
  end if

  read ( input_unit, '(a3,11x,4i14)', iostat = ios ) type, nrow, ncol, nz, nel

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_TO_ST - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading HB header record #3.'
    stop
  end if

  read ( input_unit, '(2a16,2a20)', iostat = ios ) ptrfmt, indfmt, valfmt, rhsfmt

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_TO_ST - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading HB header record #4.'
    stop
  end if

  if ( 0 < rhscrd ) then

    read ( input_unit, '(a3,11x,i14,i14)', iostat = ios ) rhstyp, nrhs, nzrhs

   if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HB_TO_ST - Fatal error!'
      write ( *, '(a)' ) '  I/O error reading HB header record #5.'
      stop
    end if

  end if

  if ( type(2:2) == 'Z' .or. type(2:2) == 'z' ) then
    skew = -1.0D+00
  else if ( type(2:2) == 'S' .or. type(2:2) == 's' ) then
    skew =  1.0D+00
  else
    skew = 0.0D+00
  end if

  sym = ( skew /= 0.0D+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,a)'   ) 'Title:  ', title
  write ( *, '(a,a)'   ) 'Key:    ', key
  write ( *, '(a,a)'   ) 'Type:   ', type
  write ( *, '(a,i14)' ) 'NROW:   ', nrow
  write ( *, '(a,i14)' ) 'NCOL:   ', ncol
  write ( *, '(a,i14)' ) 'NZ:     ', nz

  if ( 0 < rhscrd ) then
    write ( *, '(a,a)'   ) 'RHSTYP: ', rhstyp
    write ( *, '(a,i14)' ) 'NRHS:   ', nrhs
    write ( *, '(a,i14)' ) 'NZRHS:  ', nzrhs
  end if

  write ( *, '(a,l1)'  ) 'SYM:   ', sym
  write ( *, '(a,g14.6)' ) 'SKEW:  ', skew

  if ( skew == -1.0D+00 ) then
   write ( *, '(a)' ) ' '
   write ( *, '(a)' ) 'HB_TO_ST - Fatal error!'
   write ( *, '(a)' ) '  Cannot handle skew-symmetric matrices.'
   stop
  end if

  n = max ( nrow, ncol )

  if ( nmax <= ncol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_TO_ST - Fatal error!'
    write ( *, '(a)' ) '  The internal arrays are too small to'
    write ( *, '(a)' ) '  handle this matrix.  To proceed,'
    write ( *, '(a,i14)' ) '  increase NMAX to more than ', ncol
    stop
  end if

  if ( nzmax < nz ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_TO_ST - Fatal error!'
    write ( *, '(a)' ) '  The internal arrays are too small to'
    write ( *, '(a)' ) '  handle this matrix.  To proceed,'
    write ( *, '(a,i14)' ) '  increase NZMAX to at least ', nz
    stop
  end if
!
!  Read the pointers, indices, and possibly data values.
!
  read ( input_unit, ptrfmt, iostat = ios ) ptr(1:ncol+1)

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_TO_ST - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading HB pointer data.'
    stop
  end if

  read ( input_unit, indfmt, iostat = ios ) index(1:nz)

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_TO_ST - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading HB index data.'
    stop
  end if

  if ( 0 < valcrd ) then

    read ( input_unit, valfmt, iostat = ios ) value(1:nz)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HB_TO_ST - Fatal error!'
      write ( *, '(a)' ) '  I/O error reading HB value data.'
      stop
    end if

  end if

  close ( unit = input_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HB_TO_ST:'
  write ( *, '(a)' ) '  READ HBP file: "' // trim ( input_filename ) // '".'
!
!  Open the output file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace' )
!
!  Create the triplet form of the input matrix
!  STYPE = 0: unsymmetric
!  STYPE = -1: symmetric, lower triangular part present
!
  stype = - skew

  do col = 1, ncol

    do p = ptr(col), ptr(col+1) - 1

      row = index(p)

      if ( 0 < valcrd ) then
        write ( output_unit, '(i8,i8,e30.18)' ) row-1, col-1, value(p)
        if ( sym .and. row /= col ) then
  	      write ( output_unit, '(i8,i8,e30.18)' ) col-1, row-1, skew * value(p)
  	    end if
      else
        write ( output_unit, '(i8,i8,e30.18)' ) row-1, col-1, 1.0D+00
      end if

    end do

  end do

  close ( unit = output_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HB_TO_ST:'
  write ( *, '(a)' ) '  Created DSP file: "' // trim ( output_filename ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HB_TO_ST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

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

  character ( len = 8  ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
