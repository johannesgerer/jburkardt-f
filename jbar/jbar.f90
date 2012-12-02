program main

!*****************************************************************************80
!
!! MAIN reads in a data file and creates a bar graph in PPMA format.
!
!  Discussion:
!
!    A tiny little error in DBAR_TO_X, in which a "YEAR1" had
!    replaced a "YEAR2", had been messing up all my plots, but
!    it's fixed now.
!
!  Usage:
!
!    jbar FILE.DAT TITLE FILE.PPMA
!
!    where
!
!      FILE.DAT is the name of an input data file, formatted as described
!      in routine DATA_READ;  the extension "DAT" determines the type of data: 
!        'current', weekly data from today back;
!        8 characters long, 'YYYMMDD', weekly data from YYYY/MM/DD back;
!        6 characters long, 'YYYYMM', monthly data; 
!        4 characters long, 'YYYY', or other, yearly data is assumed.
!
!      TITLE is a brief title to display on the graph;
!
!      FILE.PPMA is the name of the PPMA (ASCII portable pixel map) file 
!      in which the created graph should be stored;  If this filename
!      is omitted, then it will normally be constructed by replacing
!      the extension of the input data file by ".PPMA".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: MAX_BAR = 370
  integer ( kind = 4 ), parameter :: MAX_CHUNK = 80
  integer ( kind = 4 ), parameter :: MAX_RGB = 100
  integer ( kind = 4 ), parameter :: MAX_COL = 900
  integer ( kind = 4 ), parameter :: MAX_ROW = 300
  integer ( kind = 4 ), parameter :: MAX_TITLE = 70

  integer ( kind = 4 ) b(MAX_ROW,MAX_COL)
  integer ( kind = 4 ) day
  integer ( kind = 4 ) day1
  integer ( kind = 4 ) day2
  logical, parameter :: debug = .false.
  character ( len = 14 ) dbar(MAX_BAR)
  character ( len = 100 ) filename_in
  character ( len = 100 ) filename_out
  integer ( kind = 4 ) g(MAX_ROW,MAX_COL)
  character ( len = 10 )hms
  integer ( kind = 4 ) hour
  integer ( kind = 4 ) hour1
  integer ( kind = 4 ) hour2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) inunit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ipxfargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilabel(MAX_BAR,MAX_CHUNK)
  integer ( kind = 4 ) ilen
  integer ( kind = 4 ) iperm(MAX_RGB)
  integer ( kind = 4 ) j
  character ( len = 30 ) label(MAX_BAR,MAX_CHUNK)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) minute
  integer ( kind = 4 ) minute1
  integer ( kind = 4 ) minute2
  integer ( kind = 4 ) month
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) month1
  integer ( kind = 4 ) month2
  integer ( kind = 4 ) nbar
  integer ( kind = 4 ) nchunk(MAX_BAR)
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrgb
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) numarg
  integer ( kind = 4 ) nuniq
  integer ( kind = 4 ) r(MAX_ROW,MAX_COL)
  integer ( kind = 4 ) rgb(3,MAX_RGB)
  logical s_eqi
  integer ( kind = 4 ) second
  integer ( kind = 4 ) second1
  integer ( kind = 4 ) second2
  character ( len = 10 )time
  character ( len = MAX_TITLE ) title
  character ( len = 8 ) today
  character ( len = 10 )type
  integer ( kind = 4 ) values(8)
  real ( kind = 4 ) xbar(MAX_BAR)
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ychunk(MAX_BAR,MAX_CHUNK)
  integer ( kind = 4 ) year
  integer ( kind = 4 ) year1
  integer ( kind = 4 ) year2
  real ( kind = 4 ) ymax
  real ( kind = 4 ) ymin
  character ( len = 5 ) zone
!
!  Initialize.
!
  ierror = 0

  call date_and_time ( today, time, zone, values )

  hms = time
  today = today
!
!  Initialize data.
!
  call init ( b, dbar, g, ilabel, label, MAX_BAR, MAX_CHUNK, &
    MAX_COL, MAX_RGB, MAX_ROW, MAX_TITLE, nbar, nchunk, ncol, &
    nrgb, nrow, r, rgb, title, xbar, ychunk )
!
!  Say hello.
!
  if ( debug ) then

    call timestamp ( )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JBAR'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) '  Simple PPMA graphs of Internet traffic times.'

  end if
!
!  Count the number of command-line arguments.
!
!  Old style:
!
  numarg = iargc ( )
!
!  New style:
!
! numarg = ipxfargc ( )
!
!  The input file name is the first optional command line argument.
!
  if ( numarg >= 1 ) then

    iarg = 1
!
!  Old style:
!
    call getarg ( iarg, filename_in )
!
!  New style:
!
!   call pxfgetarg ( iarg, filename_in, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'JBAR - Fatal error!'
!     write ( *, '(a)' ) '  Could not read command line argument.'
!     stop
!   end if

  else
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Please enter the input data file name.'
    read ( *, '(a)' ) filename_in

  end if
!
!  Extract the file name extension.
!
  call file_ext ( filename_in, i1, i2 )
!
!  Assume the time range based on the file extension.
!
  if ( s_eqi ( filename_in(i1:i2), 'CURRENT' ) ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  File extension is "current", so...'
    end if

    type = 'CURRENT'
    
    call s_to_ymd ( today, 'YYYYMMDD', year2, month2, day2 )
    hour2 = 24
    minute2 = 0
    second2 = 0

    year1 = year2
    month1 = month2
    day1 = day2 - 6
    hour1 = 0
    minute1 = 0
    second1 = 0

    call ymd_check_common ( year1, month1, day1, ierror )

  else if ( i2 + 1 - i1 == 4 ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  File extension is 4 characters long, so...'
    end if

    type = 'YEAR'

    call s_to_i ( filename_in(i1:i2), year1, ierror, last )

    month1 = 1
    day1 = 1
    hour1 = 0
    minute1 = 0
    second1 = 0

    year2 = year1
    month2 = 12
    day2 = 31
    hour2 = 24
    minute2 = 0
    second2 = 0

  else if ( i2 + 1 - i1 == 6 ) then

    type = 'MONTH'

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  File extension is 6 characters long, so...'
    end if
 
    call s_to_i ( filename_in(i1:i1+3), year1, ierror, last )

    call s_to_i ( filename_in(i1+4:i1+5), month1, ierror, last )

    day1 = 1
    hour1 = 0
    minute1 = 0
    second1 = 0

    year2 = year1
    month2 = month1
    day2 = month_length_common ( year1, month1 )
    hour2 = 24
    minute2 = 0
    second2 = 0

  else if ( i2 + 1 - i1 == 8 ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  File extension is 8 characters long, so...'
    end if

    type = 'WEEK'

    call s_to_ymd ( filename_in(i1:i2), 'YYYYMMDD', year2, month2, day2 )

    hour2 = 24
    minute2 = 0
    second2 = 0

    year1 = year2
    month1 = month2
    day1 = day2 - 6
    hour1 = 0
    minute1 = 0
    second1 = 0

    call ymd_check_common ( year1, month1, day1, ierror )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JBAR - Fatal error.'
    write ( *, '(a)' ) '  The file extension is unintelligible.'
    stop

  end if

  if ( debug ) then
    write ( *, '(a)' ) '...the time range is assumed to be ' // trim ( type )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'The data begins on '
    write ( *, '(i4,2x,i2,2x,i2,2x,i2,2x,i2,2x,i2)' ) &
      year1, month1, day1, hour1, minute1, second1
    write ( *, '(a)' ) 'and ends on '
    write ( *, '(i4,2x,i2,2x,i2,2x,i2,2x,i2,2x,i2)' ) &
      year2, month2, day2, hour2, minute2, second2
  end if
!
!  Get the plot title.
!  For now, we need to capitalize it entirely.
!  Also, if the title is DEBUG, that sets the DEBUG switch on.
!
  if ( numarg >= 2 ) then

    iarg = 2
!
!  Old style:
!
    call getarg ( iarg, title )
!
!  New style:
!
!   call pxfgetarg ( iarg, title, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'JBAR - Fatal error!'
!     write ( *, '(a)' ) '  Could not read command line argument.'
!     stop
!   end if

  else
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Please enter a title for the plot.'
    read ( *, '(a)' ) title

  end if

  call s_cap ( title )

  if ( title == 'DEBUG' ) then
  end if
!
!  Get the output file name.
!
  if ( numarg >= 3 ) then

    iarg = 3
!
!  Old style:
!
    call getarg ( iarg, filename_out )
!
!  New style:
!
!   call pxfgetarg ( iarg, filename_out, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'JBAR - Fatal error!'
!     write ( *, '(a)' ) '  Could not read the command line argument.'
!     stop
!   end if

  else

    if ( i1 > 1 ) then
      filename_out = filename_in(1:i2) // '.ppma'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'The output file name:'
      write ( *, '(a)' ) trim ( filename_out )

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Please enter the output PPMA graphics file name.'
      read ( *, '(a)' ) filename_out

    end if

  end if
!
!  Read the input data file.
!
  call get_unit ( inunit )

  open ( unit = inunit, file = filename_in, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JBAR - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    stop
  end if

  call data_read ( dbar, debug, ierror, inunit, label, MAX_BAR, &
    MAX_CHUNK, nbar, nchunk, ychunk )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JBAR - Fatal error!'
    write ( *, '(a)' ) '  DATA_READ error return = ', ierror
  end if

  close ( unit = inunit )
!
!  Possibly print the date data.
!
  if ( debug ) then
    call data_write ( dbar, label, MAX_BAR, MAX_CHUNK, nbar, nchunk, ychunk )
  end if
!
!  Convert the date coordinates to appropriate X values.
!
  call dbar_to_x ( year1, month1, day1, hour1, minute1, second1, dbar, nbar, &
    xbar )
!
!  Possibly print the X,Y data.
!
  if ( debug ) then
    call x_write ( dbar, nbar, xbar )
  end if
!
!  Sort the data by X value.
!
  call x_sort ( label, MAX_BAR, MAX_CHUNK, nbar, nchunk, xbar, ychunk )
!
!  For each X value, sort the data by Y value.
!
  call y_sort ( label, MAX_BAR, MAX_CHUNK, nbar, nchunk, ychunk )
!
!  Compute the X and Y ranges.
!
  xmax = xbar(1)
  xmin = xbar(1)
  do i = 2, nbar
    xmax = max ( xmax, xbar(i) )
    xmin = min ( xmin, xbar(i) )
  end do

  ymax = ychunk(1,1)
  ymin = ychunk(1,1)
  do i = 1, nbar
    do j = 1, nchunk(i)
      ymax = max ( ymax, ychunk(i,j) )
      ymin = min ( ymin, ychunk(i,j) )
    end do
  end do

  ymin = min ( ymin, 0.0E+00 )

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Data ranges:'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,2x,2g14.6)' ) 'X: ', xmin, xmax
    write ( *, '(2x,a,2x,2g14.6)' ) 'Y: ', ymin, ymax
  end if
!
!  Assign a unique index to each label.
!
  call chrlab ( MAX_BAR*MAX_CHUNK, nuniq, label, ilabel )

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'Number of unique label entries is ', nuniq
  end if

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Labels:'
    write ( *, '(a)' ) ' '

    do j = 1, MAX_CHUNK
      do i = 1, MAX_BAR
        if ( ilabel(i,j) > 0 ) then
          write ( *, '(i6,2x,a30)' ) ilabel(i,j), label(i,j)
        end if
      end do
    end do

  end if
!
!  Set colors for each label.
!
  if ( nuniq > MAX_RGB ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JBAR - Warning!'
    write ( *, '(a)' ) '  Not enough colors for all the items.'
  end if

  nrgb = min ( nuniq, MAX_RGB )

  call color_set ( MAX_RGB, nrgb, rgb )
!
!  Get a random color permutation.
!
  call perm_random2 ( nrgb, iperm )
!
!  Create the RGB data.
!
  ncol = MAX_COL
  nrow = MAX_ROW

  if ( type == 'WEEK' .or. type == 'CURRENT' ) then

    call week_plot ( b, day1, debug, g, ilabel, iperm, MAX_BAR, MAX_CHUNK, &
      MAX_RGB, month1, nbar, nchunk, ncol, nrgb, nrow, r, rgb, title, xbar, &
      xmax, xmin, ychunk, year1 )

  else if ( type == 'MONTH' ) then

    call month_plot ( b, debug, g, ilabel, iperm, MAX_BAR, MAX_CHUNK, MAX_RGB, &
      nbar, nchunk, ncol, nrgb, nrow, r, rgb, title, xbar, xmax, xmin, ychunk )

  else if ( type == 'YEAR' ) then

    call year_plot ( b, g, ilabel, iperm, MAX_BAR, MAX_CHUNK, MAX_RGB, &
      nbar, nchunk, ncol, nrgb, nrow, r, rgb, title, xbar, xmax, xmin, ychunk )

  end if
!
!  Write the RGB data to the PPMA output file.
!
  call ppma_write ( filename_out, ierror, nrow, ncol, r, g, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JBAR - Fatal error!'
    write ( *, '(a,i6)' ) '  PPMA_WRITE error return = ', ierror
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'JBAR'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine bitchr75 ( c, pattern )

!*****************************************************************************80
!
!! BITCHR75 returns a 7 by 5 bit pattern for a given character.
!
!  Discussion:
!
!    The data statements used here were generated by FONTDT.
!
!  Example:
!
!    C = 'A'
!
!    PATTERN =
!
!      0 0 1 0 0
!      0 1 0 1 0
!      1 1 0 1 1
!      1 0 0 0 1
!      1 1 1 1 1
!      1 0 0 0 1
!      1 0 0 0 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, a character whose bit pattern is desired.
!
!    Output, integer ( kind = 4 ) PATTERN(7,5), the bit pattern for the character,
!    which will be all 0's if the character is not available.
!
  implicit none

  integer ( kind = 4 ) bits(7,5,68)
  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) ipoint(0:255)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) pattern(7,5)

  data ( ipoint(i), i = 0, 255 ) / &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, &
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, &
    32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, &
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, &
    64,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 65, 66, 67, 68,  0, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 1), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0, &
      1,  1,  1,  1,  1,  0,  1, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 2), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      1,  1,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0, &
      1,  1,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 3), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  1,  0,  1,  0,  0, &
      0,  1,  1,  1,  1,  1,  0, &
      0,  0,  1,  0,  1,  0,  0, &
      0,  1,  1,  1,  1,  1,  0, &
      0,  0,  1,  0,  1,  0,  0 /

  data ((bits(i,j, 4), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  1,  0,  0,  0,  0, &
      0,  1,  0,  1,  0,  1,  0, &
      1,  1,  1,  1,  1,  1,  1, &
      0,  1,  0,  1,  0,  1,  0, &
      0,  0,  0,  0,  1,  0,  0 /

  data ((bits(i,j, 5), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  1,  0,  0,  1,  0, &
      0,  1,  1,  0,  1,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  1,  0,  1,  1,  0, &
      0,  1,  0,  0,  1,  1,  0 /

  data ((bits(i,j, 6), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  1,  0,  1,  1,  0, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  0,  1,  0,  1, &
      0,  1,  0,  0,  0,  1,  0, &
      0,  0,  0,  0,  1,  0,  1 /

  data ((bits(i,j, 7), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  1,  0,  0,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 8), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  1,  1,  1,  0,  0, &
      0,  1,  0,  0,  0,  1,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 9), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      0,  1,  0,  0,  0,  1,  0, &
      0,  0,  1,  1,  1,  0,  0 /

  data ((bits(i,j, 10), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  0,  0,  0,  0,  0, &
      0,  0,  1,  0,  0,  0,  0, &
      1,  1,  0,  0,  0,  0,  0, &
      0,  0,  1,  0,  0,  0,  0, &
      0,  1,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 11), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  1,  1,  1,  1,  1,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  0,  1,  0,  0,  0 /

  data ((bits(i,j, 12), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  1,  1, &
      0,  0,  0,  0,  0,  1,  0, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 13), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 14), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  1,  1, &
      0,  0,  0,  0,  0,  1,  1, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 15), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  1,  1, &
      0,  0,  0,  0,  1,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  1,  0,  0,  0,  0, &
      1,  1,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 16), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  1,  1,  1,  1,  0, &
      1,  0,  0,  0,  1,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  1,  0,  0,  0,  1, &
      0,  1,  1,  1,  1,  1,  0 /

  data ((bits(i,j, 17), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      0,  1,  0,  0,  0,  0,  0, &
      1,  1,  1,  1,  1,  1,  1, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 18), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  0,  1,  1, &
      1,  0,  0,  0,  1,  0,  1, &
      0,  1,  0,  1,  0,  0,  1, &
      0,  0,  1,  0,  0,  0,  1 /

  data ((bits(i,j, 19), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  0,  0,  0,  1,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      0,  1,  1,  1,  1,  1,  0 /

  data ((bits(i,j, 20), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  0,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      1,  1,  1,  1,  1,  1,  1 /

  data ((bits(i,j, 21), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  0,  0,  1,  0, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  0,  1,  1,  0 /

  data ((bits(i,j, 22), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  1,  1,  1,  1,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      0,  0,  0,  1,  0,  0,  1, &
      0,  0,  0,  1,  0,  0,  1, &
      0,  0,  0,  0,  1,  1,  0 /

  data ((bits(i,j, 23), i = 1, 7 ), j = 1, 5 ) / &
      1,  0,  0,  0,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  0, &
      1,  0,  0,  1,  1,  1,  1, &
      1,  0,  1,  0,  0,  0,  0, &
      1,  1,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 24), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  1,  0,  1,  1,  0, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      0,  1,  1,  0,  1,  1,  0 /

  data ((bits(i,j, 25), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  1,  0,  0,  1,  0, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      0,  1,  1,  1,  1,  1,  0 /

  data ((bits(i,j, 26), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  1,  0,  1,  0, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 27), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  1,  0,  1,  0, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 28), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  1,  0,  1,  0,  0, &
      0,  1,  0,  0,  0,  1,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 29), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  1,  0,  1,  0,  0, &
      0,  0,  1,  0,  1,  0,  0, &
      0,  0,  1,  0,  1,  0,  0, &
      0,  0,  1,  0,  1,  0,  0, &
      0,  0,  1,  0,  1,  0,  0 /

  data ((bits(i,j, 30), i = 1, 7 ), j = 1, 5 ) / &
      1,  0,  0,  0,  0,  0,  1, &
      0,  1,  0,  0,  0,  1,  0, &
      0,  0,  1,  0,  1,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 31), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  0,  0,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  0, &
      1,  0,  0,  1,  1,  0,  1, &
      1,  0,  1,  0,  0,  0,  0, &
      0,  1,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 32), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  1,  1,  1,  1,  0, &
      0,  1,  0,  1,  0,  0,  1, &
      0,  1,  1,  0,  1,  0,  1, &
      0,  1,  1,  0,  1,  0,  0, &
      0,  0,  1,  1,  1,  0,  0 /

  data ((bits(i,j, 33), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  1,  1,  1,  1,  1, &
      0,  1,  1,  0,  1,  0,  0, &
      1,  0,  0,  0,  1,  0,  0, &
      0,  1,  1,  0,  1,  0,  0, &
      0,  0,  1,  1,  1,  1,  1 /

  data ((bits(i,j, 34), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      0,  1,  1,  0,  1,  1,  0 /

  data ((bits(i,j, 35), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  1,  1,  1,  1,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      0,  1,  0,  0,  0,  1,  0 /

  data ((bits(i,j, 36), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      0,  1,  1,  1,  1,  1,  0 /

  data ((bits(i,j, 37), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  0,  0,  0,  1 /

  data ((bits(i,j, 38), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  1, &
      1,  0,  0,  1,  0,  0,  0, &
      1,  0,  0,  1,  0,  0,  0, &
      1,  0,  0,  1,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 39), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  1,  1,  1,  1,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  1,  0,  1, &
      0,  1,  0,  0,  1,  1,  0 /

  data ((bits(i,j, 40), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  1, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      1,  1,  1,  1,  1,  1,  1 /

  data ((bits(i,j, 41), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  1,  1,  1,  1,  1,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 42), i = 1, 7 ), j = 1, 5 ) / &
      1,  0,  0,  0,  0,  1,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  1,  1,  1,  1,  1,  0, &
      1,  0,  0,  0,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 43), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  1, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  1,  0,  1,  0,  0, &
      0,  1,  0,  0,  0,  1,  0, &
      1,  0,  0,  0,  0,  0,  1 /

  data ((bits(i,j, 44), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  1, &
      0,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  1 /

  data ((bits(i,j, 45), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  1, &
      1,  1,  0,  0,  0,  0,  0, &
      0,  0,  1,  1,  0,  0,  0, &
      1,  1,  0,  0,  0,  0,  0, &
      1,  1,  1,  1,  1,  1,  1 /

  data ((bits(i,j, 46), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  1, &
      1,  1,  0,  0,  0,  0,  0, &
      0,  0,  1,  1,  1,  0,  0, &
      0,  0,  0,  0,  0,  1,  1, &
      1,  1,  1,  1,  1,  1,  1 /

  data ((bits(i,j, 47), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  1,  1,  1,  1,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      0,  1,  1,  1,  1,  1,  0 /

  data ((bits(i,j, 48), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  1, &
      1,  0,  0,  1,  0,  0,  0, &
      1,  0,  0,  1,  0,  0,  0, &
      1,  0,  0,  1,  0,  0,  0, &
      0,  1,  1,  0,  0,  0,  0 /

  data ((bits(i,j, 49), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  1,  1,  1,  1,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  1,  0,  1, &
      1,  0,  0,  0,  0,  1,  0, &
      0,  1,  1,  1,  1,  0,  1 /

  data ((bits(i,j, 50), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  1, &
      1,  0,  0,  1,  0,  0,  0, &
      1,  0,  0,  1,  1,  0,  0, &
      1,  0,  0,  1,  0,  1,  0, &
      0,  1,  1,  0,  0,  0,  1 /

  data ((bits(i,j, 51), i = 1, 7 ), j = 1, 5 ) / &
      0,  1,  1,  0,  0,  1,  0, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      0,  1,  0,  0,  1,  1,  0 /

  data ((bits(i,j, 52), i = 1, 7 ), j = 1, 5 ) / &
      1,  0,  0,  0,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  0, &
      1,  1,  1,  1,  1,  1,  1, &
      1,  0,  0,  0,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 53), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  0, &
      0,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  1, &
      1,  1,  1,  1,  1,  1,  0 /

  data ((bits(i,j, 54), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  0,  0,  0,  0, &
      0,  0,  0,  1,  1,  1,  0, &
      0,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  1,  1,  1,  0, &
      1,  1,  1,  0,  0,  0,  0 /

  data ((bits(i,j, 55), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  0,  0, &
      0,  0,  0,  0,  1,  1,  1, &
      0,  0,  1,  1,  1,  1,  0, &
      0,  0,  0,  0,  1,  1,  1, &
      1,  1,  1,  1,  1,  0,  0 /

  data ((bits(i,j, 56), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  0,  0,  0,  1,  1, &
      0,  0,  1,  0,  1,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  1,  0,  1,  0,  0, &
      1,  1,  0,  0,  0,  1,  1 /

  data ((bits(i,j, 57), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  0,  0,  0,  0,  0, &
      0,  0,  1,  0,  0,  0,  0, &
      0,  0,  0,  1,  1,  1,  1, &
      0,  0,  1,  0,  0,  0,  0, &
      1,  1,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 58), i = 1, 7 ), j = 1, 5 ) / &
      1,  0,  0,  0,  0,  1,  1, &
      1,  0,  0,  0,  1,  0,  1, &
      1,  0,  0,  1,  0,  0,  1, &
      1,  0,  1,  0,  0,  0,  1, &
      1,  1,  0,  0,  0,  0,  1 /

  data ((bits(i,j, 59), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  1,  1,  1,  1,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 60), i = 1, 7 ), j = 1, 5 ) / &
      1,  1,  0,  0,  0,  0,  0, &
      0,  0,  1,  0,  0,  0,  0, &
      0,  0,  0,  1,  0,  0,  0, &
      0,  0,  0,  0,  1,  0,  0, &
      0,  0,  0,  0,  0,  1,  1 /

  data ((bits(i,j, 61), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  1,  1,  1,  1,  1,  1 /

  data ((bits(i,j, 62), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  1,  0,  0,  0,  0, &
      0,  1,  0,  0,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  0, &
      0,  1,  0,  0,  0,  0,  0, &
      0,  0,  1,  0,  0,  0,  0 /

  data ((bits(i,j, 63), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  1 /

  data ((bits(i,j, 64), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  0, &
      0,  1,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 65), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  1,  0,  0,  0, &
      0,  1,  1,  0,  1,  1,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 66), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0, &
      1,  1,  1,  0,  1,  1,  1, &
      0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  0,  0 /

  data ((bits(i,j, 67), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  0,  0,  0,  0,  0, &
      1,  0,  0,  0,  0,  0,  1, &
      1,  0,  0,  0,  0,  0,  1, &
      0,  1,  1,  0,  1,  1,  0, &
      0,  0,  0,  1,  0,  0,  0 /

  data ((bits(i,j, 68), i = 1, 7 ), j = 1, 5 ) / &
      0,  0,  1,  0,  0,  0,  0, &
      0,  1,  0,  0,  0,  0,  0, &
      0,  1,  0,  0,  0,  0,  0, &
      0,  0,  1,  0,  0,  0,  0, &
      0,  1,  0,  0,  0,  0,  0 /

  indx = ichar ( c )
  k = ipoint ( indx )

  if ( k == 0 ) then

    pattern(1:7,1:5) = 0

  else

    pattern(1:7,1:5) = bits(1:7,1:5,k)

  end if

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
function ch_eqi ( c1, c2 )

!*****************************************************************************80
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
!    14 August 1999
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
  character c2
  character cc1
  character cc2

  cc1 = c1
  cc2 = c2

  call ch_cap ( cc1 )
  call ch_cap ( cc2 )

  if ( cc1 == cc2 ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine chrlab ( n, nuniq, string, ident )

!*****************************************************************************80
!
!! CHRLAB makes an index array for an array of (repeated) strings.
!
!  Discussion:
!
!    CHRLAB is given an array of strings.  It assigns an integer 
!    to each unique string, and returns an equivalent array of
!    these values.  
!
!    Note that blank strings are treated specially.  Any blank
!    string gets an identifier of 0.  Blank strings are not
!    counted in the value of NUNIQ.
!
!  Example:
!
!    String    IDENT
!
!    ALPHA       1
!    ALPHA      -1
!    BETA        2
!    ALPHA      -1
!    BETA       -2
!    GAMMA       3
!    ALPHA      -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in STRING.  
!
!    Output, integer ( kind = 4 ) NUNIQ, the number of unique nonblank strings
!    in STRING.
!
!    Input, character ( len = * ) STRING(N), the list of strings to be
!    checked.
!
!    Output, integer ( kind = 4 ) IDENT(N), the identifiers assigned to the
!    strings.  If STRING(I) is blank, then IDENT(I) is 0.
!    Otherwise, if STRING(I) is the first occurrence of a 
!    given string, then it is assigned a positive identifier.
!    If STRING(I) is a later occurrence of a string, then
!    it is assigned a negative identifier, whose absolute
!    value is the identifier of the first occurrence.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ident(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nuniq
  character ( len = * ) string(n)

  nuniq = 0

  do i = 1, n

    if ( string(i) == ' ' ) then

      ident(i) = 0

    else

      do j = 1, i-1
        if ( ident(j) > 0 ) then
          if ( string(j) == string(i) ) then
            ident(i) = - ident(j)
            go to 10
          end if
        end if
      end do

      nuniq = nuniq + 1
      ident(i) = nuniq

10    continue

    end if

  end do

  return
end
subroutine color_set ( MAX_RGB, nrgb, rgb )

!*****************************************************************************80
!
!! COLOR_SET sets RGB colors for each unique item in the bar graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAX_RGB, the maximum number of colors we will compute.
!
!    Input, integer ( kind = 4 ) NRGB, the actual number of colors we will compute.
!
!    Output, real ( kind = 4 ) RGB(3,MAX_RGB), contains the RGB values, each between
!    0 and 1, for the colors we have computed.
!
  implicit none

  integer ( kind = 4 ) MAX_RGB

  real ( kind = 4 ) angle
  real ( kind = 4 ) b
  real ( kind = 4 ) g
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nrgb
  real ( kind = 4 ) r
  integer ( kind = 4 ) rgb(3,MAX_RGB)

  do j = 1, nrgb

    angle = 360.0 * real ( j - 1 ) / real ( nrgb )

    call hexcol ( angle, r, g, b )

    rgb(1,j) = int ( 255.0 * r )
    rgb(2,j) = int ( 255.0 * g )
    rgb(3,j) = int ( 255.0 * b ) 

  end do

  return
end
subroutine data_read ( dbar, debug, ierror, inunit, label, MAX_BAR, MAX_CHUNK, &
  nbar, nchunk, ychunk )

!*****************************************************************************80
!
!! DATA_READ reads in data for a bar plot.
!
!  Discussion:
!
!    The format of a single line of the data file is:
!
!      YVAL, LABEL, YMDHMS
!
!    where
!
!      YVAL is an integer Y coordinate (integer); 
!      LABEL is a label;
!      YMDHMS is a 14 character year/month/date/hour/minute/second string
!
!    Blank lines, and comments (beginning with '#') are allowed.
!
!    The X coordinate defines the location of the bar.
!
!    The bars should be drawn as follows:
!
!      Sort in ascending order all the Y coordinates associated with
!      a given X coordinate.  
!
!      Color the bar in horizontal stripes as you proceed from 0 to
!      the first Y coordinate, the second and so on.  The color used
!      on the first stripe is that associated with the first Y coordinate,
!      and so on.
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
!    Output, character ( len = 14 ) DBAR(MAX_BAR), the date coordinate values
!    associated with the bars.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    nonzero, an error occurred.
!
!    Input, integer ( kind = 4 ) INUNIT, the FORTRAN unit number from which the
!    data will be read.
!
!    Output, character ( len = 30 ) LABEL(MAX_BAR,MAX_CHUNK); LABEL(I,J) is the
!    label or tag associated with "chunk" J in bar I.
!
!    Input, integer ( kind = 4 ) MAX_BAR, the maximum number of bars.
!
!    Input, integer ( kind = 4 ) MAX_CHUNK, the maximum number of "chunks", that is,
!    Y values associated with a single bar.
!
!    Output, integer ( kind = 4 ) NBAR, the number of bars.
!
!    Output, integer ( kind = 4 ) NCHUNK(MAX_BAR).  NCHUNK(I) is the number of Y values
!    associated with bar I.
!
!    Output, real ( kind = 4 ) YCHUNK(MAX_BAR,MAX_CHUNK); YCHUNK(I,J) is a Y coordinate
!    value associate with bar I.
!
  implicit none

  integer ( kind = 4 ) MAX_BAR
  integer ( kind = 4 ) MAX_CHUNK

  character ( len = 14 ) dbar(MAX_BAR)
  logical debug
  character ( len = 14 ) dval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibar
  integer ( kind = 4 ) ichunk
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inunit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ival
  character ( len = 30 ) label(MAX_BAR,MAX_CHUNK)
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nbad
  integer ( kind = 4 ) nbadx
  integer ( kind = 4 ) nbady
  integer ( kind = 4 ) nbar
  integer ( kind = 4 ) nchunk(MAX_BAR)
  integer ( kind = 4 ) ntext
  character ( len = 80 ) string
  character ( len = 30 ) sval
  real ( kind = 4 ) ychunk(MAX_BAR,MAX_CHUNK)

  ibar = 0
  ichunk = 0
  ierror = 0
  ilo = 0
  nbad = 0
  nbadx = 0
  nbady = 0
  nbar = 0
  ntext = 0
!
!  Read a line of data.
!
  do

    read ( inunit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      exit
    end if

    ntext = ntext + 1
!
!  Ignore comments.
!
    if ( string /= ' ' ) then

      if ( string(1:1) /= '#' ) then

        ihi = 0
!
!  Read the Y coordinate, IVAL.
!
        call word_next ( string, ilo, ihi )

        if ( ihi == 0 ) then
          nbad = nbad + 1
        else

          call s_to_i ( string(ilo:ihi), ival, ierror, lchar )

          if ( ierror /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DATA_READ - Fatal error!'
            write ( *, '(a)' ) '  Could not interpret a Y coordinate.'
            return
          end if
!
!  Read the label, SVAL.
!
          call word_next ( string, ilo, ihi )

          if ( ihi == 0 ) then
            nbad = nbad + 1
          else

            sval = string(ilo:ihi)
!
!  Read the date coordinate, DVAL.
!
            call word_next ( string, ilo, ihi )

            if ( ihi == 0 ) then
              nbad = nbad + 1
            else

              dval = string(ilo:ihi)
!
!  If the X coordinate is new, add it to the XBAR list.
!  Set IBAR to the index in XBAR of the X value.
!
              ibar = 0

              do i = 1, nbar
                if ( dval == dbar(i) ) then
                  ibar = i
                end if
              end do

              if ( ibar == 0 ) then
                ibar = nbar + 1
                if ( ibar <= MAX_BAR ) then
                  nbar = ibar
                  dbar(ibar) = dval
                else
                  nbadx = nbadx + 1
                end if
              end if
!
!  Add the Y value to the list of Y values associated with the X value,
!  and the label to the label array.
!
              if ( ibar <= MAX_BAR ) then
                ichunk = nchunk(ibar) + 1
                if ( ichunk <= MAX_CHUNK ) then
                  nchunk(ibar) = nchunk(ibar) + 1
                  ychunk(ibar,ichunk) = real ( ival, kind = 4 )
                  label(ibar,ichunk) = sval
                else
                  nbady = nbady + 1
                end if
              end if

            end if

          end if

        end if

      end if

    end if

  end do

20    continue
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_READ - Report:'
    write ( *, '(a,i6)' ) '  Number of input text lines NTEXT =   ', ntext
    write ( *, '(a,i6)' ) '  Number of bad lines NBAD =           ', nbad
    write ( *, '(a,i6)' ) '  Number of bars NBAR =                ', nbar
    write ( *, '(a,i6)' ) '  Number of discarded bars =           ', nbadx
    write ( *, '(a,i6)' ) '  Number of discarded subbars =        ', nbady
  end if

  return
end
subroutine data_write ( dbar, label, MAX_BAR, MAX_CHUNK, nbar, nchunk, ychunk )

!*****************************************************************************80
!
!! DATA_WRITE prints out the data for a bar plot.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 14 ) DBAR(MAX_BAR), the date coordinate values
!    associated with the bars.
!
!    Input, character ( len = 30 ) LABEL(MAX_BAR,MAX_CHUNK); LABEL(I,J) is the
!    label or tag associated with "chunk" J in bar I.
!
!    Input, integer ( kind = 4 ) MAX_BAR, the maximum number of bars.
!
!    Input, integer ( kind = 4 ) MAX_CHUNK, the maximum number of "chunks", that is,
!    Y values associated with a single bar.
!
!    Input, integer ( kind = 4 ) NBAR, the number of bars.
!
!    Input, integer ( kind = 4 ) NCHUNK(MAX_BAR).  NCHUNK(I) is the number of Y values
!    associated with bar I.
!
!    Input, real ( kind = 4 ) YCHUNK(MAX_BAR,MAX_CHUNK); YCHUNK(I,J) is a Y coordinate
!    value associate with bar I.
!
  implicit none

  integer ( kind = 4 ) MAX_BAR
  integer ( kind = 4 ) MAX_CHUNK

  character ( len = 14 ) dbar(MAX_BAR)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character ( len = 30 ) label(MAX_BAR,MAX_CHUNK)
  integer ( kind = 4 ) nbar
  integer ( kind = 4 ) nchunk(MAX_BAR)
  real ( kind = 4 ) ychunk(MAX_BAR,MAX_CHUNK)

  do i = 1, nbar

    write ( *, '(a)' ) ' '
    write ( *, '(i3,2x,a14,2x,i3)' ) i, dbar(i), nchunk(i)

    do j = 1, nchunk(i)

      write ( *, '(2x,i3,2x,f8.3,2x,a30)' ) j, ychunk(i,j), label(i,j)

    end do

  end do

  return
end
subroutine day_borrow_common ( y, m, d )

!*****************************************************************************80
!
!! DAY_BORROW_COMMON borrows days from months in a common date.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, a year, month, and day
!    representing a date.  On input, D might be negative.  On output, 
!    M should have decreased by one month, and D gone up by the 
!    number of days in the month we "cashed in".  Y may be affected 
!    if the input value of M was 1.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) y

  do while ( d <= 0 )
 
    m = m - 1
    d = d + month_length_common ( y, m )
  
    do while ( m <= 0 ) 
      call month_borrow ( y, m )
    end do
 
  end do
 
  return
end
subroutine day_carry_common ( y, m, d )

!*****************************************************************************80
!
!! DAY_CARRY_COMMON carries days to months in a common date.
!
!  Algorithm:
!
!    While D > number of days in M:
!      decrease the day D by the number of days in the month M;
!      increase M by 1;
!      if necessary, adjust Y.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, the year, month and day.  
!    On output, D is between 1 and the number of days in M.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) days
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) y

  days = month_length_common ( y, m )

  do while ( d > days ) 
 
    d = d - days
    m = m + 1
    days = month_length_common ( y, m )
 
    do while ( m > 12 ) 
      call month_carry ( y, m )
    end do
 
  end do
 
  return
end
subroutine dbar_to_x ( year1, month1, day1, hour1, minute1, second1, dbar, &
  nbar, xbar )

!*****************************************************************************80
!
!! DBAR_TO_X converts a date into an appropriate X coordinate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) YEAR1, MONTH1, DAY1, HOUR1, MINUT1, SECOND1,
!    the YMDHMS coordinates of the first date in the plotting range.
!
!    Input, character ( len = 14 ) DBAR(NBAR), the dates associated with the
!    bars to be plotted.
!
!    Input, integer ( kind = 4 ) NBAR, the number of bars.
!
!    Output, real ( kind = 4 ) XBAR(NBAR), the X coordinates of the bars.
!
  implicit none

  integer ( kind = 4 ) nbar

  integer ( kind = 4 ) day_dif
  integer ( kind = 4 ) day1
  integer ( kind = 4 ) day2
  character ( len = 14 ) dbar(nbar)
  integer ( kind = 4 ) hour_dif
  integer ( kind = 4 ) hour1
  integer ( kind = 4 ) hour2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) minute_dif
  integer ( kind = 4 ) minute1
  integer ( kind = 4 ) minute2
  integer ( kind = 4 ) month1
  integer ( kind = 4 ) month2
  integer ( kind = 4 ) second_dif
  integer ( kind = 4 ) second1
  integer ( kind = 4 ) second2
  real ( kind = 4 ) xbar(nbar)
  integer ( kind = 4 ) year1
  integer ( kind = 4 ) year2

  do i = 1, nbar

    call s_to_ymdhms ( dbar(i), 'YYYYMMDDhhmmss', year2, month2, &
      day2, hour2, minute2, second2 )

    call ymdhms_dif_dhms ( year1, month1, day1, hour1, minute1, second1, &
                           year2, month2, day2, hour2, minute2, second2, &
                       day_dif, hour_dif, minute_dif, second_dif, ierror )

    xbar(i) = ( ( real ( second_dif, kind = 4 ) / 60.0 &
      + real ( minute_dif, kind = 4 ) ) / 60.0 &
      + real ( hour_dif, kind = 4 ) ) / 24.0 &
      + real ( day_dif, kind = 4 )

  end do

  return
end
subroutine file_ext ( filnam, i, j )

!*****************************************************************************80
!
!! FILE_EXT determines the "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!    Blanks are unusual in filenames.  This routine ignores all
!    trailing blanks, but will treat initial or internal blanks
!    as regular characters acceptable in a file name.
!
!  Example:
!
!    FILNAM     I  J
!
!    bob.for    5  7
!    N.B.C.D    7  7
!    Naomi.     0  0
!    Arthur     0  0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILNAM, a file name to be examined.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the first and last characters
!    in the file extension.
!
!    If at least one period occurs in the filename, and at least one
!    nonblank character follows that period, then I will be the index
!    of the first character after the period, and J the index of the
!    last nonblank character after the period.  The extension is
!    therefore equal to FILNAM(I:J).
!
!    Otherwise, I and J will be returned as 0, indicating that the file
!    has no extension.
!
  implicit none

  character ( len = * ) filnam
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s_index_last

  i = s_index_last ( filnam, '.' )

  if ( i /= 0 ) then

    j = len_trim ( filnam )

    if ( i == j ) then
      i = 0
      j = 0
    else
      i = i + 1
    end if

  else

    j = 0

  end if

  return
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
!    Output, integer ( kind = 4 ) IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
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
subroutine hour_borrow ( y, m, d, h )

!*****************************************************************************80
!
!! HOUR_BORROW "borrows" a day of hours.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, H, the year, month, day
!    and hour of the date.  The value of H is presumably negative, and
!    so hours will be "borrowed" to make H positive.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  do while ( h <= 0 ) 
 
    h = h + 24
    d = d - 1
 
    do while ( d <= 0 )
      call day_borrow_common ( y, m, d )
    end do
  
  end do
 
  return
end
subroutine hour_carry ( y, m, d, h )

!*****************************************************************************80
!
!! HOUR_CARRY is given a YMDH date, and carries hours to days.
!
!  Algorithm:
!
!    While H > 24:
!
!      decrease H by the number of hours in a day;
!      increase D by 1;
!      if necessary, adjust M and Y.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, H, the year, month, day
!    and hour of the date.  On input, H is presumably 24 or greater.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) y

  do while ( h > 24 )
 
    h = h - 24
    d = d + 1
 
    do while ( d > month_length_common ( y, m ) )
      call day_carry_common ( y, m, d )
    end do
 
  end do
 
  return
end
subroutine hexcol ( angle, r, g, b )

!*****************************************************************************80
!
!! HEXCOL returns a color on the perimeter of the color hexagon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) ANGLE, the angle in the color hexagon.  The sextants are
!    defined by the following points:
!        0 degrees, 1, 0, 0, red;
!       60 degrees, 1, 1, 0, yellow;
!      120 degrees, 0, 1, 0, green;
!      180 degrees, 0, 1, 1, cyan;
!      240 degrees, 0, 0, 1, blue;
!      300 degrees, 1, 0, 1, magenta.
!
!    Output, real R, G, B, RGB specifications for the color that lies
!    at the given angle, on the perimeter of the color hexagon.  One
!    value will be 1, and one value will be 0.
!
  implicit none

  real ( kind = 4 ) angle
  real ( kind = 4 ) b
  real ( kind = 4 ) g
  real ( kind = 4 ) r

  angle = mod ( angle, 360.0 )

  if ( angle < 0.0 ) then
    angle = angle + 360.0
  end if

  if ( angle <= 60.0 ) then
    r = 1.0
    g = angle / 60.0
    b = 0.0
  else if ( angle <= 120.0 ) then
    r = ( 120.0 - angle ) / 60.0
    g = 1.0
    b = 0.0
  else if ( angle <= 180.0 ) then
    r = 0.0
    g = 1.0
    b = ( angle - 120.0 ) / 60.0
  else if ( angle <= 240.0 ) then
    r = 0.0
    g = ( 240.0 - angle ) / 60.0
    b = 1.0
  else if ( angle <= 300.0 ) then
    r = ( angle - 240.0 ) / 60.0
    g = 0.0
    b = 1.0
  else if ( angle <= 360.0 ) then
    r = 1.0
    g = 0.0
    b = ( 360.0 - angle ) / 60.0
  end if

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
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
!    Input, integer ( kind = 4 ) ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) ( kind = 4 ) i
  integer ( kind = 4 ) ( kind = 4 ) i4_modp
  integer ( kind = 4 ) ( kind = 4 ) j
  integer ( kind = 4 ) ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
subroutine i_random ( ilo, ihi, i )

!*****************************************************************************80
!
!! I_RANDOM returns a random integer in a given range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ILO, IHI, the minimum and maximum acceptable values.
!
!    Output, integer ( kind = 4 ) I, the randomly chosen integer.
!
  implicit none

  logical, save :: seed = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  real ( kind = 4 ) r
  real ( kind = 4 ) rhi
  real ( kind = 4 ) rlo

  if ( .not. seed ) then
    call random_seed
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = r )
!
!  Set a real interval [RLO,RHI] which contains the integers [ILO,IHI],
!  each with a "neighborhood" of width 1.
!
  rlo = real ( ilo, kind = 4 ) - 0.5
  rhi = real ( ihi, kind = 4 ) + 0.5
!
!  Set I to the integer that is nearest the scaled value of R.
!
  i = nint ( ( 1.0 - r ) * rlo + r * rhi )
!
!  In case of oddball events at the boundary, enforce the limits.
!
  i = max ( i, ilo )
  i = min ( i, ihi )

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
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
!    Input/output, integer ( kind = 4 ) ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) ( kind = 4 ) i
  integer ( kind = 4 ) ( kind = 4 ) j
  integer ( kind = 4 ) ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine init ( b, dbar, g, ilabel, label, MAX_BAR, MAX_CHUNK, MAX_COL, &
  MAX_RGB, MAX_ROW, MAX_TITLE, nbar, nchunk, ncol, nrgb, nrow, r, rgb, title, &
  xbar, ychunk )

!*****************************************************************************80
!
!! INIT initializes the data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) B(MAX_ROW,MAX_COL), contains a value between 0 and 255
!    for the blue component of the pixel colors.
!
!    Output, character ( len = 14 ) DBAR(NBAR), the date coordinate values
!    associated with the bars.
!
!    Output, integer ( kind = 4 ) G(MAX_ROW,MAX_COL), contains a value between 0 and 255
!    for the green component of the pixel colors.
!
!    Output, integer ( kind = 4 ) ILABEL(MAX_BAR,MAX_CHUNK); ILABEL(I,J) contains the
!    integer ( kind = 4 ) tag for the corresponding string in the LABEL array.
!
!    Output, character ( len = 30 ) LABEL(MAX_BAR,MAX_CHUNK); LABEL(I,J) is the
!    label or tag associated with "chunk" J in bar I.
!
!    Input, integer ( kind = 4 ) MAX_BAR, the maximum number of bars.
!
!    Input, integer ( kind = 4 ) MAX_CHUNK, the maximum number of "chunks", that is,
!    Y values associated with a single bar.
!
!    Input, integer ( kind = 4 ) MAX_RGB, the maximum number of colors we will compute,
!    to distinguish "chunks" associated with different labels.
!
!    Input, integer ( kind = 4 ) MAX_TITLE, the maximum number of characters in TITLE.
!
!    Output, integer ( kind = 4 ) NBAR, the number of bars.
!
!    Output, integer ( kind = 4 ) NCHUNK(MAX_BAR).  NCHUNK(I) is the number of Y values
!    associated with bar I.
!
!    Output, integer ( kind = 4 ) NCOL, the actual number of columns of pixels to use
!    in the output PPMA data.
!
!    Output, integer ( kind = 4 ) NRGB, the actual number of colors we will compute
!    to distinguish "chunks" associated with different labels.
!
!    Output, integer ( kind = 4 ) NROW, the actual number of rows of pixels to use
!    in the output PPMA data.
!
!    Output, integer ( kind = 4 ) R(MAX_ROW,MAX_COL), contains a value between 0 and 255
!    for the red component of the pixel colors.
!
!    Output, integer ( kind = 4 ) RGB(3,MAX_RGB), contains the RGB values, each between
!    0 and 255, for the colors we have computed.
!
!    Output, real ( kind = 4 ) XBAR(MAX_BAR), the X coordinate values associated with
!    the bars.
!
!    Output, real ( kind = 4 ) YCHUNK(MAX_BAR,MAX_CHUNK); YCHUNK(I,J) is a Y coordinate
!    value associate with bar I.
!
  implicit none

  integer ( kind = 4 ) MAX_BAR
  integer ( kind = 4 ) MAX_CHUNK
  integer ( kind = 4 ) MAX_RGB
  integer ( kind = 4 ) MAX_COL
  integer ( kind = 4 ) MAX_ROW
  integer ( kind = 4 ) MAX_TITLE

  integer ( kind = 4 ) b(MAX_ROW,MAX_COL)
  character ( len = 14 ) dbar(MAX_BAR)
  integer ( kind = 4 ) g(MAX_ROW,MAX_COL)
  integer ( kind = 4 ) ilabel(MAX_BAR,MAX_CHUNK)
  character ( len = 30 ) label(MAX_BAR,MAX_CHUNK)
  integer ( kind = 4 ) nbar
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nchunk(MAX_BAR)
  integer ( kind = 4 ) nrgb
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) r(MAX_ROW,MAX_COL)
  integer ( kind = 4 ) rgb(3,MAX_RGB)
  character ( len = * ) title
  real ( kind = 4 ) xbar(MAX_BAR)
  real ( kind = 4 ) ychunk(MAX_BAR,MAX_CHUNK)

  b(1:MAX_ROW,1:MAX_COL) = 255
  dbar(1:MAX_BAR) = ' '
  g(1:MAX_ROW,1:MAX_COL) = 255
  ilabel(1:MAX_BAR,1:MAX_CHUNK) = 0
  label(1:MAX_BAR,1:MAX_CHUNK) = ' '
  nbar = 0
  nchunk(1:MAX_BAR) = 0
  ncol = 0
  nrgb = 0
  nrow = 0
  r(1:MAX_ROW,1:MAX_COL) = 255
  rgb(1:3,1:MAX_RGB) = 255
  title = ' '
  xbar(1:MAX_BAR) = 0.0
  ychunk(1:MAX_BAR,1:MAX_CHUNK) = 0.0

  return
end
subroutine minute_borrow ( y, m, d, h, n )

!*****************************************************************************80
!
!! MINUTE_BORROW "borrows" an hour of minutes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, H, N, the year,
!    month, day, hour and minute representing a date.  On input, N 
!    might be negative.
!    On output, H should have decreased by one, and N gone up by 60.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) y

  do while ( n < 0 ) 
 
    n = n + 60
    h = h - 1
 
    do while ( h < 0 )
      call hour_borrow ( y, m, d, h )
    end do
 
  end do
 
  return
end
subroutine minute_carry ( y, m, d, h, n )

!*****************************************************************************80
!
!! MINUTE_CARRY is given a YMDHMS date, and carries minutes to hours.
!
!  Algorithm:
!
!    While N >= 60:
!
!      decrease N by the number of minutes in an hour;
!      increase H by 1;
!      if necessary, adjust Y, M and D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, H, N, the date.  
!    On output, N is between 0 and 59.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) y

  do while ( n >= 60 ) 
 
    n = n - 60
    h = h + 1
 
    do while ( h > 24 )
      call hour_carry ( y, m, d, h )
    end do
 
  end do
 
  return
end
subroutine month_borrow ( y, m )

!*****************************************************************************80
!
!! MONTH_BORROW "borrows" a year of months while M is nonpositive.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the month and year representing a date.
!    On input, M might be negative.  On output, Y should have decreased by
!    one, and M gone up by the number of months in the year that we
!    "cashed in".  The routine knows there was no year 0.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  do while ( m <= 0 ) 
 
    m = m + 12
    y = y - 1
 
    if ( y == 0 ) then
      y = - 1
    end if
 
  end do
 
  return
end
subroutine month_carry ( y, m )

!*****************************************************************************80
!
!! MONTH_CARRY is given a YM date, and carries months to years.
!
!  Algorithm:
!
!    While M > 12:
!
!      decrease M by 12;
!      increase Y by 1;
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
!    Input/output, integer ( kind = 4 ) Y, M, the year and month.  
!    On output, M is no greater than 12.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  do while ( m > 12 ) 
 
    m = m - 12
    y = y + 1

    if ( y == 0 ) then
      y = 1
    end if
 
  end do
 
  return
end
function month_length_common ( y, m )

!*****************************************************************************80
!
!! MONTH_LENGTH_COMMON returns the number of days in a given common month.
!
!  Discussion:
!
!    The "common" calendar is meant to be the calendar which is Julian before
!    the transition date, and Gregorian afterwards, with the transition date
!    best specified as as JED = 2299160.
!
!    The routine knows that February has 28 days, except in leap years,
!    when it has 29.
!
!    In the common calendar, October 1582 had only 21 days
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year in which the month occurred.
!
!    Input, integer ( kind = 4 ) M, the number of the month for which information
!    is desired.  M should be between 1 (January) and 12 (December).
!
!    Output, integer ( kind = 4 ) MONTH_LENGTH_COMMON, the number of days in the month.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter, dimension(12) :: mdays = (/ &
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_common

  month_length_common = 0
!
!  Since this is a function, make a copy of the input arguments.
!
  m2 = m
  y2 = y
!
!  Check the input date.
!
  call ym_check_common ( y2, m2, ierror )

  if ( ierror /= 0 ) then
    month_length_common = 0
    return
  end if
!
!  Take care of the special case.
!
  if ( y2 == 1582 ) then
    if ( m2 == 10 ) then
      month_length_common = 21
      return
    end if
  end if
!
!  Get the number of days in the month.
!
  month_length_common = mdays ( m2 )
!
!  If necessary, add 1 day for February 29.
!
  if ( m2 == 2 .and. year_is_leap_common ( y2 ) ) then
    month_length_common = month_length_common + 1
  end if
 
  return
end
subroutine month_name_to_month ( month_name, month )

!*****************************************************************************80
!
!! MONTH_NAME_TO_MONTH returns the month number of a given month
!
!  Discussion:
!
!    Capitalization is ignored.  The month name has to match up to
!    the unique beginning of a month name, and the rest is ignored.
!    Here are the limits:
!
!      JAnuary
!      February
!      MARch
!      APril
!      MAY
!      JUNe
!      JULy
!      AUgust
!      September
!      October
!      November
!      December
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) MONTH_NAME, a string containing a month
!    name or abbreviation.
!
!    Output, integer ( kind = 4 ) MONTH, the number of the month, or -1 if the name
!    could not be recognized.
!
  implicit none

  integer ( kind = 4 ) month
  character ( len = * ) month_name
  character ( len = 3 ) string

  string = month_name
  call s_cap ( string )

       if ( string(1:2) == 'JA' ) then
    month = 1
  else if ( string(1:1) == 'F' ) then
    month = 2
  else if ( string(1:3) == 'MAR' ) then
    month = 3
  else if ( string(1:2) == 'AP' ) then
    month = 4
  else if ( string(1:3) == 'MAY' ) then
    month = 5
  else if ( string(1:3) == 'JUN' ) then
    month = 6
  else if ( string(1:3) == 'JUL' ) then
    month = 7
  else if ( string(1:2) == 'AU' ) then
    month = 8
  else if ( string(1:1) == 'S' ) then
    month = 9
  else if ( string(1:1) == 'O' ) then
    month = 10
  else if ( string(1:1) == 'N' ) then
    month = 11
  else if ( string(1:1) == 'D' ) then
    month = 12
  else
    month = - 1
  end if

  return
end
subroutine month_plot ( b, debug, g, ilabel, iperm, MAX_BAR, MAX_CHUNK, &
  MAX_RGB, nbar, nchunk, ncol, nrgb, nrow, r, rgb, title, xbar, xmax, xmin, &
  ychunk )

!*****************************************************************************80
!
!! MONTH_PLOT sets the RGB data for a bar graph of a month's data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) B(NROW,NOL), contains a value between 0 and 255
!    for the blue component of the pixel colors.
!
!    Output, integer ( kind = 4 ) G(NROW,NCOL), contains a value between 0 and 255
!    for the green component of the pixel colors.
!
!    Input, integer ( kind = 4 ) ILABEL(MAX_BAR,MAX_CHUNK); ILABEL(I,J) contains the
!    integer ( kind = 4 ) tag for the corresponding string in the LABEL array.
!   
!    Input, integer ( kind = 4 ) MAX_BAR, the maximum number of bars.
!
!    Input, integer ( kind = 4 ) MAX_CHUNK, the maximum number of "chunks", that is,
!    Y values associated with a single bar.
!
!    Input, integer ( kind = 4 ) MAX_RGB, the maximum number of colors we will compute,
!    to distinguish "chunks" associated with different labels.
!
!    Input, integer ( kind = 4 ) NBAR, the number of bars.
!
!    Input, integer ( kind = 4 ) NCHUNK(MAX_BAR).  NCHUNK(I) is the number of Y values
!    associated with bar I.
!
!    Input, integer ( kind = 4 ) NCOL, the actual number of columns of pixels to use
!    in the output PPMA data.
!
!    Input, integer ( kind = 4 ) NRGB, the actual number of colors we will compute
!    to distinguish "chunks" associated with different labels.
!
!    Input, integer ( kind = 4 ) NROW, the actual number of rows of pixels to use
!    in the output PPMA data.
!
!    Output, integer ( kind = 4 ) R(NROW,NCOL), contains a value between 0 and 255
!    for the red component of the pixel colors.
!
!    Input, integer ( kind = 4 ) RGB(3,MAX_RGB), contains the RGB values, each between
!    0 and 255, for the colors we have computed for labeling.
!
!    Input, character ( len = MAX_TITLE ) TITLE, a title for the plot.
!
!    Input, real ( kind = 4 ) XBAR(MAX_BAR), the X coordinate values associated with
!    the bars.
!
!    Input, real ( kind = 4 ) XMAX, XMIN, the maximum and minimum X coordinate
!    values associated with the bars.
!
!    Input, real ( kind = 4 ) YCHUNK(MAX_BAR,MAX_CHUNK); YCHUNK(I,J) is a Y coordinate
!    value associate with bar I.
!
  implicit none

  integer ( kind = 4 ) MAX_BAR
  integer ( kind = 4 ) MAX_CHUNK
  integer ( kind = 4 ) MAX_RGB
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) b(nrow,ncol)
  character c
  logical debug
  character ( len = 2 ) name2
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) icolor
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imid
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) ilabel(MAX_BAR,MAX_CHUNK)
  integer ( kind = 4 ) iperm(MAX_RGB)
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jleft
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmid
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) jrite
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) nbar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) nchunk(MAX_BAR)
  integer ( kind = 4 ) nrgb
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) pattern(7,5)
  integer ( kind = 4 ) r(nrow,ncol)
  integer ( kind = 4 ) rgb(3,MAX_RGB)
  character ( len = 80 ) string
  character ( len = * ) title
  real ( kind = 4 ) x
  real ( kind = 4 ) xbar(MAX_BAR)
  real ( kind = 4 ) xdif
  real ( kind = 4 ) xleft
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmaxg
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xming
  real ( kind = 4 ) xrite
  real ( kind = 4 ) xwide
  real ( kind = 4 ) xwide_min
  real ( kind = 4 ) y
  real ( kind = 4 ) ybot
  real ( kind = 4 ) ychunk(MAX_BAR,MAX_CHUNK)
  real ( kind = 4 ) ymaxg
  real ( kind = 4 ) yming
  real ( kind = 4 ) ytop
!
!  Draw the vertical axis.
!
  do j = 13, 14
    do i = 1, nrow
      r(i,j) = 0
      g(i,j) = 0
      b(i,j) = 0
    end do
  end do
!
!  Draw the horizontal axis.
!
  do i = nrow-15, nrow-14
    do j = 1, ncol
      r(i,j) = 0
      g(i,j) = 0
      b(i,j) = 0
    end do
  end do
!
!  Choose a width for the bars.
!
  xdif = xmax - xmin
  do i = 1, nbar
    do j = i+1, nbar
      xdif = min ( xdif, abs ( xbar(i) - xbar(j) ) )
    end do
  end do

  xwide = 0.8 * xdif
  xwide_min = 0.14

  xwide = max ( xwide, xwide_min )
!
!  The mapping is from 
!
!  ( XMIN-0.5*XWIDE-0.1, XMAX+0.5*XWIDE+0.1 ), ( YMIN-0.1, YMAX+0.1 ) to 
!  (                 31,            NCOL-5 ),  (  NROW-20,       21 )
!
  xmaxg = 31.0
  xming =  0.0
  ymaxg =  20.0
  yming =   0.0

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Graph ranges:'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,2x,2g14.6)' ) 'X: ', xming, xmaxg
    write ( *, '(2x,a,2x,2g14.6)' ) 'Y: ', yming, ymaxg
  end if

  imin = nrow - 20
  imax = 21

  jmin = 31
  jmax = ncol - 5
!
!  Make pale gray vertical lines between the days.
!
  do k = 0, 31

    x = real ( k, kind = 4 )
    call r_to_i ( j, jmax, jmin, x, xmaxg, xming )

    do i = imax, imin
      r(i,j) = 127
      g(i,j) = 127
      b(i,j) = 127
    end do

  end do
!
!  Put the day number below each month's data.
!
  i = nrow - 3

  do k = 1, 31

    x = real ( k, kind = 4 ) - 0.5 

    call r_to_i ( jmid, jmax, jmin, x, xmaxg, xming )

    write ( name2, '(i2)' ) k

    nchar = 2

    do kk = 1, nchar

      j = jmid - 0.5 * (5+1) * nchar + (kk-1) * (5+1)
      c = name2(kk:kk)
      call bitchr75 ( c, pattern )

      do i1 = 1, 7
        do j1 = 1, 5

          if ( pattern(i1,j1) == 1 ) then
            r(i+i1-7,j+j1-1) = 0
            g(i+i1-7,j+j1-1) = 0
            b(i+i1-7,j+j1-1) = 0
          end if

        end do
      end do

    end do

  end do
!
!  Put the title on top.
!  EFFECTIVELY THIS PUTS A LIMIT ON THE TITLE LENGTH, BUT YOU DON'T CHECK.
!
  i = 18
  x = 0.5 * ( xming + xmaxg )

  call r_to_i ( jmid, jmax, jmin, x, xmaxg, xming )

  call s_cap ( title )
  nchar = len_trim ( title )

  do kk = 1, nchar

    j = jmid - 0.5 * (5+1) * 2 * nchar + 2 * (kk-1) * (5+1)
    c = title(kk:kk)
    call bitchr75 ( c, pattern )

    do i1 = 1, 7
      do j1 = 1, 5

        if ( pattern(i1,j1) == 1 ) then

          r(i+2*i1-14,j+2*j1-2) = 0
          r(i+2*i1-14,j+2*j1-1) = 0
          r(i+2*i1-13,j+2*j1-2) = 0
          r(i+2*i1-13,j+2*j1-1) = 0

          g(i+2*i1-14,j+2*j1-2) = 0
          g(i+2*i1-14,j+2*j1-1) = 0
          g(i+2*i1-13,j+2*j1-2) = 0
          g(i+2*i1-13,j+2*j1-1) = 0

          b(i+2*i1-14,j+2*j1-2) = 0
          b(i+2*i1-14,j+2*j1-1) = 0
          b(i+2*i1-13,j+2*j1-2) = 0
          b(i+2*i1-13,j+2*j1-1) = 0

        end if

      end do
    end do

  end do
!
!  Put the words "HOP COUNT" along the vertical axis.
!
  j = 9

  y = 0.5 * ( ymaxg + yming )

  call r_to_i ( imid, imax, imin, y, ymaxg, yming )

  string = 'HOP COUNT'
  nchar = len_trim ( string )

  do kk = 1, nchar

    i = imid + 0.5*(5+1)*nchar - (kk-1) * (5+1)
    c = string(kk:kk)
    call bitchr75 ( c, pattern )

    do i1 = 1, 7
      do j1 = 1, 5

        if ( pattern(i1,j1) == 1 ) then
          r(i+1-j1,j+i1-7) = 0
          g(i+1-j1,j+i1-7) = 0
          b(i+1-j1,j+i1-7) = 0
        end if

      end do
    end do

  end do
!
!  Make horizontal grid lines.
!
  ny = 6

  do k = 1, ny

    y = ( ( ny - k ) * yming + ( k - 1 ) * ymaxg ) / real ( ny - 1, kind = 4 )

    call r_to_i ( i, imax, imin, y, ymaxg, yming )

    do j = jmin, jmax
      r(i,j) = 127
      g(i,j) = 127
      b(i,j) = 127
    end do

  end do
!
!  Put the hop count to the left of each horizontal grid line.
!
  jmid = 23

  do k = 2, ny

    y = ( ( ny - k ) * yming + ( k - 1 ) * ymaxg ) / real ( ny - 1, kind = 4 )

    call r_to_i ( i, imax, imin, y, ymaxg, yming )

    write ( name2, '(i2)' ) 5 * ( k - 1 )

    nchar = 2

    do kk = 1, nchar

      j = jmid - 0.5 * (5+1) * nchar + (kk-1) * (5+1)
      c = name2(kk:kk)
      call bitchr75 ( c, pattern )

      do i1 = 1, 7
        do j1 = 1, 5

          if ( pattern(i1,j1) == 1 ) then
            r(i+i1-7,j+j1-1) = 0
            g(i+i1-7,j+j1-1) = 0
            b(i+i1-7,j+j1-1) = 0
          end if

        end do
      end do

    end do

  end do
!
!  Draw the bars.
!
  do i = 1, nbar

    xleft = xbar(i) - 0.5 * xwide
    xrite = xbar(i) + 0.5 * xwide

    call r_to_i ( jleft, jmax, jmin, xleft, xmaxg, xming )
    call r_to_i ( jrite, jmax, jmin, xrite, xmaxg, xming )

    ytop = 0.0
    call r_to_i ( itop, imax, imin, ytop, ymaxg, yming )

    do j = 1, nchunk(i)

      ybot = ytop
      ibot = itop

      ytop = ychunk(i,j)

      icolor = abs ( ilabel(i,j) )

      if ( icolor > nrgb ) then
        icolor = mod ( icolor, nrgb ) + 1
      end if

      icolor = iperm(icolor)

      if ( ytop > ybot ) then

        call r_to_i ( itop, imax, imin, ytop, ymaxg, yming )

        if ( ibot <= itop ) then
          isgn = 1
        else
          isgn = - 1
        end if

        if ( jleft <= jrite ) then
          jsgn = 1
        else
          jsgn = - 1
        end if

        do ii = ibot, itop, isgn
          do jj = jleft, jrite, jsgn
            r(ii,jj) = rgb(1,icolor)
            g(ii,jj) = rgb(2,icolor)
            b(ii,jj) = rgb(3,icolor)
          end do
        end do

      end if

    end do
  end do

  return
end
subroutine month_to_month_name ( m, month_name )

!*****************************************************************************80
!
!! MONTH_TO_MONTH_NAME returns the name of a given month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of the month, which should
!    be between 1 and 12.
!
!    Output, character ( len = * ) MONTH_NAME, a string containing as much of
!    the month's name as will fit.  To get the typical 3-letter abbreviations
!    for the months, simply declare
!      character ( len = 3 ) MONTH_NAME
!    or pass in MONTH_NAME(1:3).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 9 ), parameter, dimension(12) :: name = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)

  if ( m < 1 .or. m > 12 ) then

    do i = 1, len ( month_name )
      month_name(i:i) = '?'
    end do

  else

    month_name = name(m)

  end if

  return
end
subroutine perm_random2 ( n, iarray )

!*****************************************************************************80
!
!! PERM_RANDOM2 selects a random permutation of N elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    James Filliben
!    National Bureau of Standards.  
!
!  Reference:
!
!    K L Hoffman and D R Shier,
!    Algorithm 564,
!    A Test Problem Generator for Discrete Linear L1 Approximation Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 615-617.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of the array.
!
!    Output, integer ( kind = 4 ) IARRAY(N), a random permutation of the
!    digits 1 though N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) iadd
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) j
  real ( kind = 4 ) u

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_RANDOM2 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal input value of N  = ', n
    write ( *, '(a)' ) '  N must be at least 1!'
    stop
  end if
 
  if ( n == 1 ) then
    iarray(1) = 1
    return
  end if
 
  do i = 1, n
    iarray(i) = i
  end do
 
  do i = 1, n

    call i_random ( 1, n, iadd )

    j = i + iadd

    if ( j > n ) then
      j = j - n
    end if
 
    if ( i /= j ) then
      call i4_swap ( iarray(i), iarray(j) )
    end if
 
  end do
 
  return
end
subroutine ppm_check_data ( r, g, b, ierror, maxcol, ncol, nrow )

!*****************************************************************************80
!
!! PPM_CHECK_DATA checks pixel data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) R(NROW,NCOL), G(NROW,NCOL), B(NROW,NCOL), contains the 
!    RGB pixel data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, the data is illegal.
!
!    Input, integer ( kind = 4 ) MAXCOL, the maximum value.
!
!    Input, integer ( kind = 4 ) NCOL, NROW, the number of rows and columns of data.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) b(nrow,ncol)
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) r(nrow,ncol)

  ierror = 0
!
!  Make sure no color is negative nor greater than MAXCOL.
!
  do i = 1, nrow
    do j = 1, ncol

      if ( r(i,j) < 0 .or. g(i,j) < 0 .or. b(i,j) < 0 ) then

        ierror = 1
        return

      end if

     if ( r(i,j) > maxcol .or. g(i,j) > maxcol .or. b(i,j) > maxcol ) then

        ierror = 1
        return

      end if

    end do
  end do

  return
end
subroutine ppma_write ( filename, ierror, nrow, ncol, r, g, b )

!*****************************************************************************80
!
!! PPMA_WRITE writes an ASCII portable pixel map file.
!
!  Discussion:
!
!    PPMA files can be viewed by XV.
!
!    Programs to convert files to this format include:
!
!      GIFTOPPM  - GIF file
!      PGMTOPPM  - Portable Gray Map file
!      PICTTOPPM - Macintosh PICT file
!      XPMTOPPM  - X11 pixmap file
!
!    Various programs can convert other formats to PPM format, including:
!
!      BMPTOPPM - Microsoft Windows BMP file.
!
!    A PPM file can also be converted to other formats, by programs:
!
!      PPMTOACAD - AutoCAD file
!      PPMTOGIF  - GIF file
!      PPMTOPGM  - Portable Gray Map file
!      PPMTOPICT - Macintosh PICT file
!      PPMTOPUZZ - X11 puzzle file
!      PPMTORGB3 - 3 Portable Gray Map files
!      PPMTOXPM  - X11 pixmap file
!      PPMTOYUV  - Abekas YUV file
!    
!  Example:
!
!    P3
!    # feep.ppm
!    4 4
!    15
!     0  0  0    0  0  0    0  0  0   15  0 15
!     0  0  0    0 15  7    0  0  0    0  0  0
!     0  0  0    0  0  0    0 15  7    0  0  0
!    15  0 15    0  0  0    0  0  0    0  0  0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file to which
!    the data should be written.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) R(NROW,NCOL), G(NROW,NCOL), B(NROW,NCOL), contain
!    the red, green and blue values of each pixel.  These should
!    be positive.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) b(nrow,ncol)
  logical debug
  character ( len = * ) filename
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = 2 ) magic
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) r(nrow,ncol)

  debug = .false.
  ierror = 0
!
!  Compute the maximum color value.
!
  maxcol = r(1,1)

  do i = 1, nrow
    do j = 1, ncol

      maxcol = max ( maxcol, r(i,j) )
      maxcol = max ( maxcol, g(i,j) )
      maxcol = max ( maxcol, b(i,j) )

    end do
  end do
!
!  Check the data.
!
  call ppm_check_data ( r, g, b, ierror, maxcol, ncol, nrow )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Bad data!'
    ierror = 1
    return
  end if
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    ierror = 2
    return
  end if
!
!  Write the data.
!
  magic = 'P3'

  write ( output_unit, '(a2)' ) magic
  write ( output_unit, '(a)' ) '# ASCII PPM file created by PPMA_WRITE.'
  write ( output_unit, '(i5,2x,i5)' ) ncol, nrow
  write ( output_unit, '(i5)' ) maxcol

  do i = 1, nrow
    do jlo = 1, ncol, 4
      jhi = min ( jlo + 3, ncol )
      write ( output_unit, '(12i5)' ) ( r(i,j), g(i,j), b(i,j), j = jlo, jhi )
    end do
  end do
!
!  Close the file.
!
  close ( unit = output_unit )
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE - Note:'
    write ( *, '(a)' ) '  The data was checked and written.'
    write ( *, '(a,i6)' ) '  Number of data rows NROW =    ', nrow
    write ( *, '(a,i6)' ) '  Number of data columns NCOL = ', ncol
    write ( *, '(a,i6)' ) '  Maximum color value MAXCOL =  ', maxcol
  end if

  return
end
subroutine r_to_i ( ix, ixmax, ixmin, x, xmax, xmin )

!*****************************************************************************80
!
!! R_TO_I maps real X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
!
!  Discussion:
!
!    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
!    IX := min ( IX, max ( IXMIN, IXMAX ) )
!    IX := max ( IX, min ( IXMIN, IXMAX ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IX, the value in the range [IXMIN,IXMAX] that
!    corresponds to X.
!
!    Input, integer ( kind = 4 ) IXMAX, IXMIN, the allowed range of the output
!    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
!    It is not necessary that IXMIN be less than IXMAX.
!
!    Input, real ( kind = 4 ) X, the real number to be converted.
!
!    Input, real ( kind = 4 ) XMAX, XMIN, the real range.  XMAX and XMIN must not be
!    equal.  It is not necessary that XMIN be less than XMAX.
!
  implicit none

  integer ( kind = 4 ) ix
  integer ( kind = 4 ) ixmax
  integer ( kind = 4 ) ixmin
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin

  if ( xmax == xmin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R_TO_I - Fatal error!'
    write ( *, '(a)' ) '  XMAX = XMIN, making a zero divisor!'
    write ( *, '(a,g14.6)' ) '  XMAX = ', xmax
    write ( *, '(a,g14.6)' ) '  XMIN = ', xmin
    stop
  end if

  ix = nint 
  ( 
    real 
    (
      ( ( xmax - x ) * ixmin + ( x - xmin ) * ixmax ) &
      / ( xmax - xmin ), kind = 4
    ) 
  )

  ix = min ( ix, max ( ixmin, ixmax ) )
  ix = max ( ix, min ( ixmin, ixmax ) )

  return
end
subroutine s_cap ( string )

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
!    16 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) STRING, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) string

  nchar = len ( string )

  do i = 1, nchar

    c = string(i:i)
    call ch_cap ( c )
    string(i:i) = c

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
function s_index_last ( s, sub )

!*****************************************************************************80
!
!! S_INDEX_LAST finds the LAST occurrence of a given substring.
!
!  Discussion:
!
!    It returns the location in the string at which the substring SUB is
!    first found, or 0 if the substring does not occur at all.
!
!    The routine is also trailing blank insensitive.  This is very
!    important for those cases where you have stored information in
!    larger variables.  If S is of length 80, and SUB is of
!    length 80, then if S = 'FRED' and SUB = 'RED', a match would
!    not be reported by the standard FORTRAN INDEX, because it treats
!    both variables as being 80 characters long!  This routine assumes that
!    trailing blanks represent garbage!
!
!    This means that this routine cannot be used to find, say, the last
!    occurrence of a substring 'A ', since it assumes the blank space
!    was not specified by the user, but is, rather, padding by the
!    system.  However, as a special case, this routine can properly handle
!    the case where either S or SUB is all blanks.
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
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character ( len = * ) SUB, the substring to search for.
!
!    Output, integer ( kind = 4 ) S_INDEX_LAST.  0 if SUB does not occur in
!    the string.  Otherwise S_INDEX_LAST = I, where S(I:I+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the last place
!    this happens.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) llen1
  integer ( kind = 4 ) llen2
  character ( len = * ) s
  integer ( kind = 4 ) s_index_last
  character ( len = * ) sub

  s_index_last = 0

  llen1 = len_trim ( s )
  llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN
!
  if ( llen1 == 0 ) then
    llen1 = len ( s )
  end if

  if ( llen2 == 0 ) then
    llen2 = len ( sub )
  end if

  if ( llen2 > llen1 ) then
    return
  end if

  do j = 1, llen1+1-llen2

    i = llen1 + 2 - llen2 - j

    if ( s(i:i+llen2-1) == sub ) then
      s_index_last = i
      return
    end if

  end do

  return
end
subroutine s_to_hms ( s, pat, hour, minute, second )

!*****************************************************************************80
!
!! S_TO_HMS converts a string into a H:M:S date.
!
!  Discussion:
!
!    The characters in PAT indicate where the data is stored.  A particular
!    letter, such as "H", indicates, an hour field, while the number of "H"
!    characters indicates the width of the field.
!
!    The codes are:
!
!    'H' or 'h' an hour field
!    'M' or 'm' a minute field
!    'S' or 's' a second field
!
!  Example:
!
!    S                    PAT
!    ------------         ------------
!    '230859'             'hhmmss'
!    '10:30'              'hh:mm'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string containing the data.
!
!    Input, character ( len = * ) PAT, describes how the data is stored.
!    PAT must be the same length as S.
!
!    Output, integer ( kind = 4 ) HOUR, MINUTE, SECOND, the hour, minute and second
!    represented by the string.  Any item not read from the string will
!    have a value of -1.
!
  implicit none

  integer ( kind = 4 ) hour
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  integer ( kind = 4 ) minute
  character ( len = * ) pat
  character ( len = * ) s
  integer ( kind = 4 ) second

  hour = - 1
  minute = - 1
  second = - 1

  length = min ( len ( s ), len ( pat ) )

  ihi = 0

10    continue

  ilo = ihi + 1
  ihi = ilo

20    continue

  if ( ihi + 1 <= length ) then
    if ( pat(ihi+1:ihi+1) == pat(ilo:ilo) ) then
      ihi = ihi + 1
      go to 20
    end if
  end if

       if ( pat(ilo:ilo) == 'H' .or. pat(ilo:ilo) == 'h' ) then
    call s_to_i ( s(ilo:ihi), hour, ierror, last )
  else if ( pat(ilo:ilo) == 'M' .or. pat(ilo:ilo) == 'm' ) then
    call s_to_i ( s(ilo:ihi), minute, ierror, last )
  else if ( pat(ilo:ilo) == 'S' .or. pat(ilo:ilo) == 's' ) then
    call s_to_i ( s(ilo:ihi), second, ierror, last )
  end if

  if ( ihi < length ) then
    go to 10
  end if

  return
end
subroutine s_to_i ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I reads an integer value from a string.
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
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If STRING is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine s_to_r ( s, r, ierror, lchar )

!*****************************************************************************80
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
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
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
!    Output, real ( kind = 4 ) R, the real value that was read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 4 ) r
  real ( kind = 4 ) rbot
  real ( kind = 4 ) rexp
  real ( kind = 4 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

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
        rtop = 10.0E+00 * rtop + real ( ndig, kind = 4 )
      else if ( ihave == 5 ) then
        rtop = 10.0E+00 * rtop + real ( ndig, kind = 4 )
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
subroutine s_to_ymd ( s, pat, year, month, day )

!*****************************************************************************80
!
!! S_TO_YMD converts a string into a YMD date.
!
!  Discussion:
!
!    The characters in PAT indicate where the day, month and year data
!    is stored.  A particular letter, such as "Y", indicates, a year
!    field, while the number of "Y" characters indicates the width of
!    the field.  
!
!    The codes are:
!
!    'Y' or 'y', a year field
!    'M' or 'm', a numeric month field
!    'N' or 'n', a literal month field
!    'D' or 'd', a day field
!
!  Example:
!
!    S              PAT
!    ------------   ------------
!    '19991031'     'YYYYMMDD'
!    '10-31-99'     'MM-DD-YY'
!    '10-31-99'     'MM/DD/YY'
!    'Oct 31 1999'  'NNN DD YYYY'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string containing the data.
!
!    Input, character ( len = * ) PAT, describes how the data is stored.
!    PAT must be the same length as S. 
!
!    Output, integer ( kind = 4 ) YEAR, MONTH, DAY, the YMD date
!    represented by the string.  Any item not read from the string will
!    have a value of -1.
!
  implicit none

  integer ( kind = 4 ) day
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  integer ( kind = 4 ) month
  character ( len = * ) pat
  character ( len = * ) s
  integer ( kind = 4 ) year

  day = - 1
  month = - 1
  year = - 1

  length = min ( len ( s ), len ( pat ) )

  ihi = 0

10    continue

  ilo = ihi + 1
  ihi = ilo

20    continue

  if ( ihi + 1 <= length ) then
    if ( pat(ihi+1:ihi+1) == pat(ilo:ilo) ) then
      ihi = ihi + 1
      go to 20
    end if
  end if

       if ( pat(ilo:ilo) == 'Y' .or. pat(ilo:ilo) == 'y' ) then
    call s_to_i ( s(ilo:ihi), year, ierror, last )
  else if ( pat(ilo:ilo) == 'M' .or. pat(ilo:ilo) == 'm' ) then
    call s_to_i ( s(ilo:ihi), month, ierror, last )
  else if ( pat(ilo:ilo) == 'N' .or. pat(ilo:ilo) == 'n' ) then
    call month_name_to_month ( s(ilo:ihi), month )
  else if ( pat(ilo:ilo) == 'D' .or. pat(ilo:ilo) == 'd' ) then
    call s_to_i ( s(ilo:ihi), day, ierror, last )
  end if

  if ( ihi < length ) then
    go to 10
  end if

  return
end
subroutine s_to_ymdhms ( s, pat, year, month, day, hour, minute, second )

!*****************************************************************************80
!
!! S_TO_YMDHMS converts a string into a YMD H:M:S date.
!
!  Discussion:
!
!    The characters in PAT indicate where the data is stored.  A particular 
!    letter, such as "Y", indicates, a year field, while the number of "Y" 
!    characters indicates the width of the field.  
!
!    The codes are:
!
!    'Y' a year field
!    'M' a numeric month field
!    'N' a literal month field
!    'D' a day field
!    'h' an hour field
!    'm' a minute field
!    's' a second field
!
!  Example:
!
!    S                    PAT
!    ------------         ------------
!    '19991031230859'     'YYYYMMDDhhmmss'
!    '10-31-99'           'MM-DD-YY'
!    '10-31-99'           'MM/DD/YY'
!    'Oct 31 1999'        'NNN DD YYYY'
!    '10:30'              'hh:mm'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string containing the data.
!
!    Input, character ( len = * ) PAT, describes how the data is stored.
!    PAT must be the same length as S.
!
!    Output, integer ( kind = 4 ) YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, the YMDHMS
!    date represented by the string.  Any item not 
!    read from the string will have a value of -1.
!
  implicit none

  integer ( kind = 4 ) day
  integer ( kind = 4 ) hour
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  integer ( kind = 4 ) minute
  integer ( kind = 4 ) month
  character ( len = * ) pat
  character ( len = * ) s
  integer ( kind = 4 ) second
  integer ( kind = 4 ) year

  day = - 1
  hour = - 1
  minute = - 1
  month = - 1
  second = - 1
  year = - 1

  length = min ( len ( s ), len ( pat ) )

  ihi = 0

10    continue

  ilo = ihi + 1
  ihi = ilo

20    continue

  if ( ihi + 1 <= length ) then
    if ( pat(ihi+1:ihi+1) == pat(ilo:ilo) ) then
      ihi = ihi + 1
      go to 20
    end if
  end if

       if ( pat(ilo:ilo) == 'Y' ) then
    call s_to_i ( s(ilo:ihi), year, ierror, last )
  else if ( pat(ilo:ilo) == 'M' ) then
    call s_to_i ( s(ilo:ihi), month, ierror, last )
  else if ( pat(ilo:ilo) == 'N' ) then
    call month_name_to_month ( s(ilo:ihi), month )
  else if ( pat(ilo:ilo) == 'D' ) then
    call s_to_i ( s(ilo:ihi), day, ierror, last )
  else if ( pat(ilo:ilo) == 'h' ) then
    call s_to_i ( s(ilo:ihi), hour, ierror, last )
  else if ( pat(ilo:ilo) == 'm' ) then
    call s_to_i ( s(ilo:ihi), minute, ierror, last )
  else if ( pat(ilo:ilo) == 's' ) then
    call s_to_i ( s(ilo:ihi), second, ierror, last )
  end if

  if ( ihi < length ) then
    go to 10
  end if

  return
end
subroutine second_borrow ( y, m, d, h, n, s )

!*****************************************************************************80
!
!! SECOND_BORROW "borrows" a minute of seconds.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, H, N, S, the YMDHMS date.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) y

  do while ( s < 0 )
 
    s = s + 60
    n = n - 1
 
    do while ( n < 0 ) 
      call minute_borrow ( y, m, d, h, n )
    end do
 
  end do
 
  return
end
subroutine second_carry ( y, m, d, h, n, s )

!*****************************************************************************80
!
!! SECOND_CARRY is given a YMDHMS date, and carries seconds to minutes.
!
!  Algorithm:
!
!    While S >= 60:
!
!      decrease S by 60;
!      increase N by 1;
!      if necessary, adjust H, D, M and Y.
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
!    Input/output, integer ( kind = 4 ) Y, M, D, H, N, S, 
!    the year, month, day, hours, minutes, seconds,
!    On output, S is between 0 and 59.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) y

  do while ( s >= 60 ) 
 
    s = s - 60
    n = n + 1
 
    do while ( n >= 60 )
      call minute_carry ( y, m, d, h, n )
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
subroutine week_plot ( b, day1, debug, g, ilabel, iperm, MAX_BAR, MAX_CHUNK, &
  MAX_RGB, month1, nbar, nchunk, ncol, nrgb, nrow, r, rgb, title, xbar, xmax, &
  xmin, ychunk, year1 )

!*****************************************************************************80
!
!! WEEK_PLOT sets the RGB data for a bar graph of a week's data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) B(NROW,NOL), contains a value between 0 and 255
!    for the blue component of the pixel colors.
!
!    Output, integer ( kind = 4 ) G(NROW,NCOL), contains a value between 0 and 255
!    for the green component of the pixel colors.
!
!    Input, integer ( kind = 4 ) ILABEL(MAX_BAR,MAX_CHUNK); ILABEL(I,J) contains the
!    integer ( kind = 4 ) tag for the corresponding string in the LABEL array.
!
!    Input, integer ( kind = 4 ) MAX_BAR, the maximum number of bars.
!
!    Input, integer ( kind = 4 ) MAX_CHUNK, the maximum number of "chunks", that is,
!    Y values associated with a single bar.
!
!    Input, integer ( kind = 4 ) MAX_RGB, the maximum number of colors we will compute,
!    to distinguish "chunks" associated with different labels.
!
!    Input, integer ( kind = 4 ) NBAR, the number of bars.
!
!    Input, integer ( kind = 4 ) NCHUNK(MAX_BAR).  NCHUNK(I) is the number of Y values
!    associated with bar I.
!
!    Input, integer ( kind = 4 ) NCOL, the actual number of columns of pixels to use
!    in the output PPMA data.
!
!    Input, integer ( kind = 4 ) NRGB, the actual number of colors we will compute
!    to distinguish "chunks" associated with different labels.
!
!    Input, integer ( kind = 4 ) NROW, the actual number of rows of pixels to use
!    in the output PPMA data.
!
!    Output, integer ( kind = 4 ) R(NROW,NCOL), contains a value between 0 and 255
!    for the red component of the pixel colors.
!
!    Input, integer ( kind = 4 ) RGB(3,MAX_RGB), contains the RGB values, each between
!    0 and 255, for the colors we have computed for labeling.
!
!    Input, character ( len = MAX_TITLE ) TITLE, a title for the plot.
!
!    Input, real ( kind = 4 ) XBAR(MAX_BAR), the X coordinate values associated with
!    the bars.
!
!    Input, real ( kind = 4 ) XMAX, XMIN, the maximum and minimum X coordinate
!    values associated with the bars.
!
!    Input, real ( kind = 4 ) YCHUNK(MAX_BAR,MAX_CHUNK); YCHUNK(I,J) is a Y coordinate
!    value associate with bar I.
!
  implicit none

  integer ( kind = 4 ) MAX_BAR
  integer ( kind = 4 ) MAX_CHUNK
  integer ( kind = 4 ) MAX_RGB
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) b(nrow,ncol)
  character c
  integer ( kind = 4 ) day
  integer ( kind = 4 ) day1
  logical debug
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) icolor
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imid
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) ilabel(MAX_BAR,MAX_CHUNK)
  integer ( kind = 4 ) iperm(MAX_RGB)
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jleft
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmid
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) jrite
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) month
  integer ( kind = 4 ) month1
  character ( len = 2 ) name2
  integer ( kind = 4 ) nbar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) nchunk(MAX_BAR)
  integer ( kind = 4 ) nrgb
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) pattern(7,5)
  integer ( kind = 4 ) r(nrow,ncol)
  integer ( kind = 4 ) rgb(3,MAX_RGB)
  character ( len = 80 ) string
  character ( len = * ) title
  real ( kind = 4 ) x
  real ( kind = 4 ) xbar(MAX_BAR)
  real ( kind = 4 ) xdif
  real ( kind = 4 ) xleft
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmaxg
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xming
  real ( kind = 4 ) xrite
  real ( kind = 4 ) xwide
  real ( kind = 4 ) xwide_min
  real ( kind = 4 ) y
  real ( kind = 4 ) ybot
  real ( kind = 4 ) ychunk(MAX_BAR,MAX_CHUNK)
  integer ( kind = 4 ) year
  integer ( kind = 4 ) year1
  real ( kind = 4 ) ymaxg
  real ( kind = 4 ) yming
  real ( kind = 4 ) ytop
!
!  Draw the vertical axis.
!
  do j = 13, 14
    do i = 1, nrow
      r(i,j) = 0
      g(i,j) = 0
      b(i,j) = 0
    end do
  end do
!
!  Draw the horizontal axis.
!
  do i = nrow-15, nrow-14
    do j = 1, ncol
      r(i,j) = 0
      g(i,j) = 0
      b(i,j) = 0
    end do
  end do
!
!  Round XMIN down, XMAX up, to whole numbers.
!
  if ( xmin < 0.0 ) then
    xmin = real ( int ( xmin - 0.99 ), kind = 4 )
  else
    xmin = real ( int ( xmin ), kind = 4 )
  end if

  if ( xmax + 0.99 < 0.0 ) then
    xmax = real ( int ( xmax ), kind = 4 )
  else
    xmax = real ( int ( xmax + 0.99 ), kind = 4 )
  end if

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Adjusted X data range:'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,2x,2g14.6)' ) 'X: ', xmin, xmax
  end if
!
!  Choose XWIDE, a "reasonable" width for the bars.
!
  xdif = xmax - xmin
  do i = 1, nbar
    do j = i+1, nbar
      xdif = min ( xdif, abs ( xbar(i) - xbar(j) ) )
    end do
  end do

  xwide = 0.8 * xdif
  xwide_min = 0.14

  xwide = max ( xwide, xwide_min )
!
!  The mapping is from
!
!  ( XMIN-0.5*XWIDE-0.1, XMAX+0.5*XWIDE+0.1 ), ( YMIN-0.1, YMAX+0.1 ) to
!
!                 JMIN               JMAX           IMIN      IMAX
!  (                31,            NCOL-5 ),  (  NROW-20,       21 )
!
  xming =   0.0
  xmaxg =   7.0

  yming =   0.0
  ymaxg =  20.0

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Graph ranges:'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,2x,2g14.6)' ) 'X: ', xming, xmaxg
    write ( *, '(2x,a,2x,2g14.6)' ) 'Y: ', yming, ymaxg
  end if

  imin = nrow - 20
  imax = 21

  jmin = 31
  jmax = ncol - 5
!
!  Make pale gray vertical lines between the days.
!
  do k = 0, 7

    x = real ( k, kind = 4 )
    call r_to_i ( j, jmax, jmin, x, xmaxg, xming )

    do i = imax, imin
      r(i,j) = 127
      g(i,j) = 127
      b(i,j) = 127
    end do

  end do
!
!  Put the day number below each day's data.
!
  i = nrow - 3

  do k = 1, 7

    x = real ( k, kind = 4 ) - 0.5

    call r_to_i ( jmid, jmax, jmin, x, xmaxg, xming )

    day = day1 + k - 1
    month = month1
    year = year1

    call ymd_check_common ( year, month, day, ierror )

    write ( name2, '(i2)' ) day

    nchar = 2

    do kk = 1, nchar

      j = jmid - 0.5 * (5+1) * nchar + (kk-1) * (5+1)
      c = name2(kk:kk)
      call bitchr75 ( c, pattern )

      do i1 = 1, 7
        do j1 = 1, 5

          if ( pattern(i1,j1) == 1 ) then
            r(i+i1-7,j+j1-1) = 0
            g(i+i1-7,j+j1-1) = 0
            b(i+i1-7,j+j1-1) = 0
          end if

        end do
      end do

    end do

  end do
!
!  Put the title on top.
!  EFFECTIVELY THIS PUTS A LIMIT ON THE TITLE LENGTH, BUT YOU DON'T CHECK.
!
  i = 18
  x = 0.5 * ( xming + xmaxg )

  call r_to_i ( jmid, jmax, jmin, x, xmaxg, xming )

  call s_cap ( title )
  nchar = len_trim ( title )

  do kk = 1, nchar

    j = jmid - 0.5 * (5+1) * 2 * nchar + 2 * (kk-1) * (5+1)
    c = title(kk:kk)
    call bitchr75 ( c, pattern )

    do i1 = 1, 7
      do j1 = 1, 5

        if ( pattern(i1,j1) == 1 ) then

          r(i+2*i1-14,j+2*j1-2) = 0
          r(i+2*i1-14,j+2*j1-1) = 0
          r(i+2*i1-13,j+2*j1-2) = 0
          r(i+2*i1-13,j+2*j1-1) = 0

          g(i+2*i1-14,j+2*j1-2) = 0
          g(i+2*i1-14,j+2*j1-1) = 0
          g(i+2*i1-13,j+2*j1-2) = 0
          g(i+2*i1-13,j+2*j1-1) = 0

          b(i+2*i1-14,j+2*j1-2) = 0
          b(i+2*i1-14,j+2*j1-1) = 0
          b(i+2*i1-13,j+2*j1-2) = 0
          b(i+2*i1-13,j+2*j1-1) = 0

        end if

      end do
    end do

  end do
!
!  Put the words "HOP COUNT" along the vertical axis.
!
  j = 9

  y = 0.5 * ( ymaxg + yming )

  call r_to_i ( imid, imax, imin, y, ymaxg, yming )

  string = 'HOP COUNT'
  nchar = len_trim ( string )

  do kk = 1, nchar

    i = imid + 0.5*(5+1)*nchar - (kk-1) * (5+1)
    c = string(kk:kk)
    call bitchr75 ( c, pattern )

    do i1 = 1, 7
      do j1 = 1, 5

        if ( pattern(i1,j1) == 1 ) then
          r(i+1-j1,j+i1-7) = 0
          g(i+1-j1,j+i1-7) = 0
          b(i+1-j1,j+i1-7) = 0
        end if

      end do
    end do

  end do
!
!  Make horizontal grid lines.
!
  ny = 6

  do k = 1, ny

    y = ( ( ny - k ) * yming + ( k - 1 ) * ymaxg ) / real ( ny - 1, kind = 4 )

    call r_to_i ( i, imax, imin, y, ymaxg, yming )

    do j = jmin, jmax
      r(i,j) = 127
      g(i,j) = 127
      b(i,j) = 127
    end do

  end do
!
!  Put the hop count to the left of each horizontal grid line.
!
  jmid = 23

  do k = 2, ny

    y = ( ( ny - k ) * yming + ( k - 1 ) * ymaxg ) / real ( ny - 1, kind = 4 )

    call r_to_i ( i, imax, imin, y, ymaxg, yming )

    write ( name2, '(i2)' ) 5 * ( k - 1 )

    nchar = 2

    do kk = 1, nchar

      j = jmid - 0.5 * (5+1) * nchar + (kk-1) * (5+1)
      c = name2(kk:kk)
      call bitchr75 ( c, pattern )

      do i1 = 1, 7
        do j1 = 1, 5

          if ( pattern(i1,j1) == 1 ) then
            r(i+i1-7,j+j1-1) = 0
            g(i+i1-7,j+j1-1) = 0
            b(i+i1-7,j+j1-1) = 0
          end if

        end do
      end do

    end do

  end do
!
!  Draw the bars.
!
  do i = 1, nbar

    xleft = xbar(i) - 0.5 * xwide - xmin
    xrite = xbar(i) + 0.5 * xwide - xmin

    call r_to_i ( jleft, jmax, jmin, xleft, xmaxg, xming )
    call r_to_i ( jrite, jmax, jmin, xrite, xmaxg, xming )

    ytop = 0.0
    call r_to_i ( itop, imax, imin, ytop, ymaxg, yming )

    do j = 1, nchunk(i)

      ybot = ytop
      ibot = itop

      ytop = ychunk(i,j)

      icolor = abs ( ilabel(i,j) )

      if ( icolor > nrgb ) then
        icolor = mod ( icolor, nrgb ) + 1
      end if

      icolor = iperm(icolor)

      if ( ytop > ybot ) then

        call r_to_i ( itop, imax, imin, ytop, ymaxg, yming )

        if ( ibot <= itop ) then
          isgn = 1
        else
          isgn = - 1
        end if

        if ( jleft <= jrite ) then
          jsgn = 1
        else
          jsgn = - 1
        end if

        do ii = ibot, itop, isgn
          do jj = jleft, jrite, jsgn
            r(ii,jj) = rgb(1,icolor)
            g(ii,jj) = rgb(2,icolor)
            b(ii,jj) = rgb(3,icolor)
          end do
        end do

      end if

    end do
  end do

  return
end
subroutine word_next ( s, ilo, ihi )

!*****************************************************************************80
!
!! WORD_NEXT finds the next (blank separated) word in a string.
!
!  Discussion:
!
!    This routine is usually used repetitively on a fixed string.  On each
!    call, it accepts IHI, the index of the last character of the
!    previous word extracted from the string.
!
!    It then computes ILO and IHI, the first and last characters of
!    the next word in the string.
!
!    It is assumed that words are separated by one or more spaces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string of words to be analyzed.
!
!    Output, integer ( kind = 4 ) ILO is the location of the first character of the
!    next word, or 0 if there was no next word.
!
!    Input/output, integer ( kind = 4 ) IHI.
!
!    On input, IHI is taken to be the LAST character of the
!    PREVIOUS word, or 0 if the first word is sought.
!
!    On output, IHI is the index of the last character of
!    the next word, or 0 if there was no next word.
!
  implicit none

  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) length
  character ( len = * ) s

  length = len_trim ( s )
!
!  Find ILO, the index of the first nonblank character after
!  (the old value of) IHI.
!
  if ( ihi < 0 ) then
    ilo = 0
  else
    ilo = ihi
  end if

  do

    ilo = ilo + 1

    if ( length < ilo ) then
      ilo = 0
      ihi = 0
      return
    end if

    if ( s(ilo:ilo) /= ' ') then
      exit
    end if

  end do
!
!  Find IHI, the index of the next blank character, or end of line.
!
  ihi = ilo

  do

    ihi = ihi + 1

    if ( length <= ihi ) then
      ihi = length
      return
    end if

    if ( s(ihi:ihi) == ' ' ) then
      exit
    end if

  end do
!
!  Decrement IHI to point to the previous, nonblank, character.
!
  ihi = ihi - 1

  return
end
subroutine x_sort ( label, MAX_BAR, MAX_CHUNK, nbar, nchunk, xbar, ychunk )

!*****************************************************************************80
!
!! X_SORT rearranges the data so that the X values are increasing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = 30 ) LABEL(MAX_BAR,MAX_CHUNK); LABEL(I,J) 
!    is the label or tag associated with "chunk" J in bar I.
!
!    Input, integer ( kind = 4 ) MAX_BAR, the maximum number of bars.
!
!    Input, integer ( kind = 4 ) MAX_CHUNK, the maximum number of "chunks", that is,
!    Y values associated with a single bar.
!
!    Input, integer ( kind = 4 ) NBAR, the number of bars.
!
!    Input/output, integer ( kind = 4 ) NCHUNK(MAX_BAR).  NCHUNK(I) is the number 
!    of Y values associated with bar I.
!
!    Input/output, real ( kind = 4 ) XBAR(MAX_BAR), the X coordinate values associated 
!    with the bars.
!
!    Input/output, real ( kind = 4 ) YCHUNK(MAX_BAR,MAX_CHUNK); YCHUNK(I,J) is a Y 
!    coordinate value associate with bar I.
!
  implicit none

  integer ( kind = 4 ) MAX_BAR
  integer ( kind = 4 ) MAX_CHUNK

  character ( len = 30 ) ctemp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 30 ) label(MAX_BAR,MAX_CHUNK)
  integer ( kind = 4 ) nbar
  integer ( kind = 4 ) nchunk(MAX_BAR)
  real ( kind = 4 ) temp
  real ( kind = 4 ) xbar(MAX_BAR)
  real ( kind = 4 ) ychunk(MAX_BAR,MAX_CHUNK)

  do i = 1, nbar
    do j = 1, i-1

      if ( xbar(j) > xbar(i) ) then

        temp = xbar(j)
        xbar(j) = xbar(i)
        xbar(i) = temp

        itemp = nchunk(j)
        nchunk(j) = nchunk(i)
        nchunk(i) = itemp

        do k = 1, MAX_CHUNK

          temp = ychunk(i,k)
          ychunk(i,k) = ychunk(j,k)
          ychunk(j,k) = temp

          ctemp = label(i,k)
          label(i,k) = label(j,k)
          label(j,k) = ctemp

        end do

      end if

    end do
  end do

  return
end
subroutine x_write ( dbar, nbar, xbar )

!*****************************************************************************80
!
!! X_WRITE prints out the X data for the plot.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 14 ) DBAR(NBAR), the date coordinate values
!    associated  with the bars.
!
!    Input, integer ( kind = 4 ) NBAR, the number of bars.
!
!    Input, real ( kind = 4 ) XBAR(NBAR), the X coordinate values associated 
!    with the bars.
!
  implicit none

  integer ( kind = 4 ) nbar

  character ( len = 14 ) dbar(nbar)
  integer ( kind = 4 ) i
  real ( kind = 4 ) xbar(nbar)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  YYYYMMDDhhmmss      X'
  write ( *, '(a)' ) ' '

  do i = 1, nbar

    write ( *, '(i3,2x,a14,2x,f8.3)' ) i, dbar(i), xbar(i)

  end do

  return
end
subroutine y_sort ( label, MAX_BAR, MAX_CHUNK, nbar, nchunk, ychunk )

!*****************************************************************************80
!
!! Y_SORT sorts the Y data in each bar so that it is increasing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = 30 ) LABEL(MAX_BAR,MAX_CHUNK); LABEL(I,J) is 
!    the label or tag associated with "chunk" J in bar I.
!
!    Input, integer ( kind = 4 ) MAX_BAR, the maximum number of bars.
!
!    Input, integer ( kind = 4 ) MAX_CHUNK, the maximum number of "chunks", that is,
!    Y values associated with a single bar.
!
!    Input, integer ( kind = 4 ) NBAR, the number of bars.
!
!    Input/output, integer ( kind = 4 ) NCHUNK(MAX_BAR).  NCHUNK(I) is the number 
!    of Y values associated with bar I.
!
!    Input/output, real ( kind = 4 ) YCHUNK(MAX_BAR,MAX_CHUNK); YCHUNK(I,J) is a 
!    Y coordinate value associate with bar I.
!
  implicit none

  integer ( kind = 4 ) MAX_BAR
  integer ( kind = 4 ) MAX_CHUNK

  character ( len = 30 ) ctemp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 30 ) label(MAX_BAR,MAX_CHUNK)
  integer ( kind = 4 ) nbar
  integer ( kind = 4 ) nchunk(MAX_BAR)
  real ( kind = 4 ) temp
  real ( kind = 4 ) ychunk(MAX_BAR,MAX_CHUNK)

  do i = 1, nbar
    do j = 1, nchunk(i)
      do k = 1, j-1

        if ( ychunk(i,j) < ychunk(i,k) ) then

          temp = ychunk(i,j)
          ychunk(i,j) = ychunk(i,k)
          ychunk(i,k) = temp

          ctemp = label(i,j)
          label(i,j) = label(i,k)
          label(i,k) = ctemp

        end if

      end do
    end do
  end do

  return
end
subroutine year_check ( year, ierror )

!*****************************************************************************80
!
!! YEAR_CHECK checks a year for reasonableness.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) YEAR, the year, which must not be 0.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if YEAR is legal, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) year

  if ( year /= 0 ) then
    ierror = 0
  else
    ierror = 1
  end if

  return
end
function year_is_leap_common ( y )

!*****************************************************************************80
!
!! YEAR_IS_LEAP_COMMON returns TRUE if the given common year was a leap year.
!
!  Discussion:
!
!    The "common" calendar is meant to be the calendar which is Julian before
!    the transition date, and Gregorian afterwards, with the transition date
!    best specified as as JED = 2299160.
!
!  Algorithm:
!
!    If ( the year is less than 0 ) then
!
!      if the year+1 is divisible by 4 then
!        the year is a leap year.
!
!    else if ( the year is 0 ) then
!
!      the year is not a leap year ( in fact, it's illegal )
!
!    else if ( the year is no greater than 1582 ) then 
!
!      if the year is divisible by 4 then
!        the year is a leap year.
!
!    else if (
!      the year is divisible by 4 and
!      ( the year is not divisible by 100
!      or
!      the year is divisible by 400 )
!      ) then
!        the year is a leap year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, logical YEAR_IS_LEAP_COMMON, TRUE if the year was a leap year,
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2
  logical year_is_leap_common

  if ( y == 0 ) then
    year_is_leap_common = .false.
    return
  end if
!
!  BC years have to have 1 added to them to make a proper leap year
!  evaluation.
!
  if ( y < 0 ) then
    y2 = y + 1
  else
    y2 = y
  end if
 
  if ( y2 <= 1582 ) then
 
    if ( i4_modp ( y2, 4 ) == 0 ) then
      year_is_leap_common = .true.
    else
      year_is_leap_common = .false.
    end if
 
  else
 
    if ( i4_modp ( y2, 400 ) == 0 ) then
      year_is_leap_common = .true.
    else if ( i4_modp ( y2, 100 ) == 0 ) then
      year_is_leap_common = .false.
    else if ( i4_modp ( y2, 4 ) == 0 ) then
      year_is_leap_common = .true.
    else
      year_is_leap_common = .false.
    end if
 
  end if
 
  return
end
function year_length_common ( y )

!*****************************************************************************80
!
!! YEAR_LENGTH_COMMON returns the number of days in a given common year.
!
!  Discussion:
!
!    The "common" calendar is meant to be the calendar which is Julian before
!    the transition date, and Gregorian afterwards, with the transition date
!    best specified as as JED = 2299160.
!
!    If Y is 0, then the routine returns 0, reflecting the fact that
!    there was officially no year 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y, the year to be checked.
!
!    Output, integer ( kind = 4 ) YEAR_LENGTH_COMMON, the number of days in the year.
!
  implicit none

  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_common
  logical year_is_leap_common

  if ( y == 0 ) then
    year_length_common = 0
  else if ( y == 1582 ) then
    year_length_common = 355
  else if ( year_is_leap_common ( y ) ) then
    year_length_common = 366
  else
    year_length_common = 365
  end if
 
  return
end
subroutine year_plot ( b, g, ilabel, iperm, MAX_BAR, MAX_CHUNK, &
  MAX_RGB, nbar, nchunk, ncol, nrgb, nrow, r, rgb, title, xbar, xmax, xmin, &
  ychunk )

!*****************************************************************************80
!
!! YEAR_PLOT sets the RGB data for the bar graph of a year's data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) B(NROW,NOL), contains a value between 0 and 255
!    for the blue component of the pixel colors.
!
!    Output, integer ( kind = 4 ) G(NROW,NCOL), contains a value between 0 and 255
!    for the green component of the pixel colors.
!
!    Input, integer ( kind = 4 ) ILABEL(MAX_BAR,MAX_CHUNK); ILABEL(I,J) contains the
!    integer ( kind = 4 ) tag for the corresponding string in the LABEL array.
!
!    Input, integer ( kind = 4 ) MAX_BAR, the maximum number of bars.
!
!    Input, integer ( kind = 4 ) MAX_CHUNK, the maximum number of "chunks", that is,
!    Y values associated with a single bar.
!
!    Input, integer ( kind = 4 ) MAX_RGB, the maximum number of colors we will compute,
!    to distinguish "chunks" associated with different labels.
!
!    Input, integer ( kind = 4 ) NBAR, the number of bars.
!
!    Input, integer ( kind = 4 ) NCHUNK(MAX_BAR).  NCHUNK(I) is the number of Y values
!    associated with bar I.
!
!    Input, integer ( kind = 4 ) NCOL, the actual number of columns of pixels to use
!    in the output PPMA data.
!
!    Input, integer ( kind = 4 ) NRGB, the actual number of colors we will compute
!    to distinguish "chunks" associated with different labels.
!
!    Input, integer ( kind = 4 ) NROW, the actual number of rows of pixels to use
!    in the output PPMA data.
!
!    Output, integer ( kind = 4 ) R(NROW,NCOL), contains a value between 0 and 255
!    for the red component of the pixel colors.
!
!    Input, integer ( kind = 4 ) RGB(3,MAX_RGB), contains the RGB values, each between
!    0 and 1, for the colors we have computed for labeling.
!
!    Input, character ( len = _TITLE ) TITLE, a title for the plot.
!
!    Input, real ( kind = 4 ) XBAR(MAX_BAR), the X coordinate values associated with
!    the bars.
!
!    Input, real ( kind = 4 ) XMAX, XMIN, the maximum and minimum X coordinate
!    values associated with the bars.
!
!    Input, real ( kind = 4 ) YCHUNK(MAX_BAR,MAX_CHUNK); YCHUNK(I,J) is a Y coordinate
!    value associate with bar I.
!
  implicit none

  integer ( kind = 4 ) MAX_BAR
  integer ( kind = 4 ) MAX_CHUNK
  integer ( kind = 4 ) MAX_RGB
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) b(nrow,ncol)
  character c
  integer ( kind = 4 ) daysum
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) icolor
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imid
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) ilabel(MAX_BAR,MAX_CHUNK)
  integer ( kind = 4 ) iperm(MAX_RGB)
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jleft
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmid
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) jrite
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ), parameter, dimension ( 12 ) :: monlen = (/ &
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  character ( len = 3 ) monname
  character ( len = 2 ) name2
  integer ( kind = 4 ) nbar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) nchunk(MAX_BAR)
  integer ( kind = 4 ) nrgb
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) pattern(7,5)
  integer ( kind = 4 ) r(nrow,ncol)
  integer ( kind = 4 ) rgb(3,MAX_RGB)
  character ( len = 80 ) string
  character ( len = * ) title
  real ( kind = 4 ) x
  real ( kind = 4 ) xbar(MAX_BAR)
  real ( kind = 4 ) xdif
  real ( kind = 4 ) xleft
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmaxg
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xming
  real ( kind = 4 ) xrite
  real ( kind = 4 ) xwide
  real ( kind = 4 ) xwide_min
  real ( kind = 4 ) y
  real ( kind = 4 ) ybot
  real ( kind = 4 ) ychunk(MAX_BAR,MAX_CHUNK)
  real ( kind = 4 ) ymaxg
  real ( kind = 4 ) yming
  real ( kind = 4 ) ytop
!
!  Draw the vertical axis.
!
  do j = 13, 14
    do i = 1, nrow
      r(i,j) = 0
      g(i,j) = 0
      b(i,j) = 0
    end do
  end do
!
!  Draw the horizontal axis.
!
  do i = nrow-15, nrow-14
    do j = 1, ncol
      r(i,j) = 0
      g(i,j) = 0
      b(i,j) = 0
    end do
  end do
!
!  Choose a width for the bars.
!
  xdif = xmax - xmin
  do i = 1, nbar
    do j = i+1, nbar
      xdif = min ( xdif, abs ( xbar(i) - xbar(j) ) )
    end do
  end do

  xwide = 0.8 * xdif
  xwide_min = 1.0
  xwide = max ( xwide, xwide_min )
!
!  The mapping is from
!  ( XMIN-0.5*XWIDE-0.1, XMAX+0.5*XWIDE+0.1 ), ( YMIN-0.1, YMAX+0.1 ) to
!  (                 31,            NCOL-5 ),  (  NROW-20,       21 )
!
  xmaxg = 365.0
  xming =  0.0
  ymaxg =  20.0
  yming =   0.0

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Graph ranges:'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,2x,2g14.6)' ) 'X: ', xming, xmaxg
    write ( *, '(2x,a,2x,2g14.6)' ) 'Y: ', yming, ymaxg
  end if

  imin = nrow - 20
  imax = 21

  jmin = 31
  jmax = ncol - 5
!
!  Make pale gray vertical lines between the months.
!
  daysum = 0

  do k = 1, 13

    x = real ( daysum, kind = 4 ) + 0.5
    call r_to_i ( j, jmax, jmin, x, xmaxg, xming )

    do i = imax, imin
      r(i,j) = 127
      g(i,j) = 127
      b(i,j) = 127
    end do

    if ( k <= 12 ) then
      daysum = daysum + monlen(k)
    end if

  end do
!
!  Put the month names below each month's data.
!
  i = nrow - 3
  daysum = 0

  do k = 1, 12

    x = real ( daysum, kind = 4 ) + 0.5 * monlen(k)
    daysum = daysum + monlen(k)

    call r_to_i ( jmid, jmax, jmin, x, xmaxg, xming )
    call month_to_month_name ( k, monname )

    nchar = 3

    do kk = 1, nchar

      j = jmid - 0.5 * (5+1) * nchar + (kk-1) * (5+1)
      c = monname(kk:kk)
      call bitchr75 ( c, pattern )

      do i1 = 1, 7
        do j1 = 1, 5

          if ( pattern(i1,j1) == 1 ) then
            r(i+i1-7,j+j1-1) = 0
            g(i+i1-7,j+j1-1) = 0
            b(i+i1-7,j+j1-1) = 0
          end if

        end do
      end do

    end do

  end do
!
!  Put the title on top.
!
  i = 18
  x = 0.5 * ( xming + xmaxg )

  call r_to_i ( jmid, jmax, jmin, x, xmaxg, xming )

  call s_cap ( title )
  nchar = len_trim ( title )

  do kk = 1, nchar

    j = jmid - 0.5 * (5+1) * 2 * nchar + 2 * (kk-1) * (5+1)
    c = title(kk:kk)
    call bitchr75 ( c, pattern )

    do i1 = 1, 7
      do j1 = 1, 5

        if ( pattern(i1,j1) == 1 ) then

          r(i+2*i1-14,j+2*j1-2) = 0
          r(i+2*i1-14,j+2*j1-1) = 0
          r(i+2*i1-13,j+2*j1-2) = 0
          r(i+2*i1-13,j+2*j1-1) = 0

          g(i+2*i1-14,j+2*j1-2) = 0
          g(i+2*i1-14,j+2*j1-1) = 0
          g(i+2*i1-13,j+2*j1-2) = 0
          g(i+2*i1-13,j+2*j1-1) = 0

          b(i+2*i1-14,j+2*j1-2) = 0
          b(i+2*i1-14,j+2*j1-1) = 0
          b(i+2*i1-13,j+2*j1-2) = 0
          b(i+2*i1-13,j+2*j1-1) = 0

        end if

      end do
    end do

  end do
!
!  Put the words "HOP COUNT" along the vertical axis.
!
  j = 9

  y = 0.5 * ( ymaxg + yming )

  call r_to_i ( imid, imax, imin, y, ymaxg, yming )

  string = 'HOP COUNT'
  nchar = len_trim ( string )

  do kk = 1, nchar

    i = imid + 0.5*(5+1)*nchar - (kk-1) * (5+1)
    c = string(kk:kk)
    call bitchr75 ( c, pattern )

    do i1 = 1, 7
      do j1 = 1, 5

        if ( pattern(i1,j1) == 1 ) then
          r(i+1-j1,j+i1-7) = 0
          g(i+1-j1,j+i1-7) = 0
          b(i+1-j1,j+i1-7) = 0
        end if

      end do
    end do

  end do
!
!  Make horizontal grid lines.
!
  ny = 6

  do k = 1, ny

    y = ( ( ny - k ) * yming + ( k - 1 ) * ymaxg ) / real ( ny - 1, kind = 4 )

    call r_to_i ( i, imax, imin, y, ymaxg, yming )

    do j = jmin, jmax
      r(i,j) = 127
      g(i,j) = 127
      b(i,j) = 127
    end do

  end do
!
!  Put the hop count to the left of each horizontal grid line.
!
  jmid = 23

  do k = 2, ny

    y = ( ( ny - k ) * yming + ( k - 1 ) * ymaxg ) / real ( ny - 1, kind = 4 )

    call r_to_i ( i, imax, imin, y, ymaxg, yming )

    write ( name2, '(i2)' ) 5 * ( k - 1 )

    nchar = 2

    do kk = 1, nchar

      j = jmid - 0.5 * (5+1) * nchar + (kk-1) * (5+1)
      c = name2(kk:kk)
      call bitchr75 ( c, pattern )

      do i1 = 1, 7
        do j1 = 1, 5

          if ( pattern(i1,j1) == 1 ) then
            r(i+i1-7,j+j1-1) = 0
            g(i+i1-7,j+j1-1) = 0
            b(i+i1-7,j+j1-1) = 0
          end if

        end do
      end do

    end do

  end do
!
!  Draw the bars.
!
  do i = 1, nbar

    xleft = xbar(i) - 0.5 * xwide
    xrite = xbar(i) + 0.5 * xwide

    call r_to_i ( jleft, jmax, jmin, xleft, xmaxg, xming )
    call r_to_i ( jrite, jmax, jmin, xrite, xmaxg, xming )

    ytop = 0.0
    call r_to_i ( itop, imax, imin, ytop, ymaxg, yming )

    if ( debug ) then
      write ( *, '(a)'       ) ' '
      write ( *, '(a,i6)'    ) 'Bar    ', i
      write ( *, '(a,i6)'    ) 'Chunks ', nchunk(i)
      write ( *, '(a,g14.6)' ) 'X =    ', xbar(i)
    end if

    do j = 1, nchunk(i)

      ybot = ytop
      ibot = itop

      ytop = ychunk(i,j)

      icolor = abs ( ilabel(i,j) )

      if ( nrgb < icolor ) then
        icolor = mod ( icolor, nrgb ) + 1
      end if

      icolor = iperm(icolor)

      if ( ybot < ytop ) then

        call r_to_i ( itop, imax, imin, ytop, ymaxg, yming )

        if ( ibot <= itop ) then
          isgn = 1
        else
          isgn = - 1
        end if

        if ( jleft <= jrite ) then
          jsgn = 1
        else
          jsgn = - 1
        end if

        do ii = ibot, itop, isgn
          do jj = jleft, jrite, jsgn
            r(ii,jj) = rgb(1,icolor)
            g(ii,jj) = rgb(2,icolor)
            b(ii,jj) = rgb(3,icolor)
          end do
        end do

      end if

    end do
  end do

  return
end
subroutine ym_check_common ( y, m, ierror )

!*****************************************************************************80
!
!! YM_CHECK_COMMON checks a YM date for reasonableness.
!
!  Discussion:
!
!    If the month is less than 1, then the month is incremented
!    by 12, and the year decremented by 1, repeatedly, until
!    the month is greater than or equal to 1.
!
!    If the month is greater than 12, then the month is decremented
!    by 12, and the year incremented by 1, repeatedly, until the
!    month is less than or equal to 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, the YM date.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was found in the date,
!    and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y
!
!  Check the year.
!
  call year_check ( y, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  The output month must be at least 1 and no more than 12.
!
  do while ( m < 1 ) 
    call month_borrow ( y, m )
  end do
 
  do while ( m > 12 ) 
    call month_carry ( y, m )
  end do
 
  return
end
subroutine ymd_check_common ( y, m, d, ierror )

!*****************************************************************************80
!
!! YMD_CHECK_COMMON checks a common YMD date for reasonableness.
!
!  Discussion:
!
!    Certain simple errors in dates will be corrected, such as
!      "31 September 1996"
!    which will become
!      "1 October 1996".
!
!    The routine also knows that in the common calendar, the dates 
!    15 October 1582 through 24 October 1582 are illegal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, the year of the date, which may
!    be altered when correcting M or D.  Y should not
!    be zero.
!
!    Input/output, integer ( kind = 4 ) M, the month of the date, which may
!    be altered.  M should be between 1 and 12, but
!    the routine will try to correct values that are too big or too small.
!
!    Input/output, integer ( kind = 4 ) D, the day of the date, which may be
!    altered.  D should be positive, and between 1 and
!    the number of days in the month M.  However, the routine will
!    try to correct values that are too big or small.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if the date is legal.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) y
!
!  Check the month and year.  After this call, M is
!  guaranteed to be between 1 and 12.
!
  call ym_check_common ( y, m, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Check if the day number is too small.
!
  do while ( d <= 0 ) 
    call day_borrow_common ( y, m, d )
  end do
 
  do while ( d > month_length_common ( y, m ) ) 
    call day_carry_common ( y, m, d )
  end do
!
!  Now make sure that the date does not fall in the
!  Julian-to-Gregorian calendar switchover limbo.
!
  if ( y == 1582 ) then
    if ( m == 10 ) then
      if ( 5 <= d .and. d <= 14 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'YMD_CHECK_COMMON - Warning!'
        write ( *, '(a)' ) '  Illegal date!'
      end if
    end if
  end if

  return
end
subroutine ymdhms_check ( y, m, d, h, n, s, ierror )

!*****************************************************************************80
!
!! YMDHMS_CHECK checks an YMDHMS date for reasonableness.
!
!  Discussion:
!
!    The routine will correct certain simple errors in dates, such as
!      "11:03:42 31 September 1996"
!    which will become
!      "11:03:42 1 October 1996".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y, M, D, H, N, S.
!    These items have the obvious meanings.
!    The routine may change any of these values to more reasonable
!    values.
!
!    Output, integer ( kind = 4 ) IERROR, is 0 if no error was detected in the
!    date, and 1 otherwise.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) y
!
!  Check that the second is between 0 and 59.
!  N may get bumped up or down.
!
  do while ( s < 0 ) 
    call second_borrow ( y, m, d, h, n, s )
  end do

  do while ( s >= 60 ) 
    call second_carry ( y, m, d, h, n, s )
  end do
!
!  Check that the minute is between 0 and 59.
!  H may get bumped up or down.
!
  do while ( n < 0 ) 
    call minute_borrow ( y, m, d, h, n )
  end do

  do while ( n >= 60 )
    call minute_carry ( y, m, d, h, n )
  end do
!
!  Check that the hour is between 1 and 24.
!  D may get bumped up or down.
!
  do while ( h < 1 )
    call hour_borrow ( y, m, d, h )
  end do

  do while ( h > 24 ) 
    call hour_carry ( y, m, d, h )
  end do
!
!  Now make adjustments to D, M, and Y.
!
  call ymd_check_common ( y, m, d, ierror )

  return
end
subroutine ymdhms_compare ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, n2, s2, &
  cmp, ierror )

!*****************************************************************************80
!
!! YMDHMS_COMPARE compares two YMDHMS dates.
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
!    Input, integer ( kind = 4 ) Y1, M1, D1, H1, N1, S1, the first date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2, H2, N2, S2, the second date.
!
!    Output, character CMP:
!    '<' if date 1 precedes date 2;
!    '=' if date 1 equals date 2;
!    '>' if date 1 follows date 2;
!    '?' if one of the dates was illegal.
!
!    Output, integer ( kind = 4 ) IERROR, is 1 if either date was illegal, and
!    0 otherwise.
!
  implicit none

  character cmp
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) h1
  integer ( kind = 4 ) h2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  cmp = '?'
!
!  Make local copies of the input, and then check them.
!  We need local copies because the checking routine can
!  change the input values.
!
  call ymdhms_check ( y1, m1, d1, h1, n1, s1, ierror )

  if ( ierror /= 0 ) then
    return
  end if

  call ymdhms_check ( y2, m2, d2, h2, n2, s2, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Compare years...
!
  if ( y1 < y2 ) then
    cmp = '<'
  else if ( y1 > y2 ) then
    cmp = '>'
  else
!
!  ...if necessary, compare months in equal years...
!
    if ( m1 < m2 ) then
      cmp = '<'
    else if ( m1 > m2 ) then
      cmp = '>'
    else
!
!  ...if necessary, compare days in equal months...
!
      if ( d1 < d2 ) then
        cmp = '<'
      else if ( d1 > d2 ) then
        cmp = '>'
      else
!
!  ...if necessary, compare hours in equal days...
!
        if ( h1 < h2 ) then
          cmp = '<'
        else if ( h1 > h2 ) then
          cmp = '>'
        else
!
!  ...if necessary, compare minutes in equal hours...
!
          if ( n1 < n2 ) then
            cmp = '<'
          else if ( n1 > n2 ) then
            cmp = '>'
          else
!
!  ...if necessary, compare seconds in equal minutes...
!
            if ( s1 < s2 ) then
              cmp = '<'
            else if ( s1 > s2 ) then
              cmp = '>'
            else
              cmp = '='
            end if
 
          end if
        end if
      end if
    end if
  end if
 
  return
end
subroutine ymdhms_dif_dhms ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, n2, s2, &
  nday, nhour, nminute, nsecond, ierror )

!*****************************************************************************80
!
!! YMDHMS_DIF_DHMS computes the DHMS difference between two YMDHMS dates.
!
!  Discussion:
!
!    The result is reported in days, minutes, hours and seconds.
!    The result is POSITIVE if the second date is later than the first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Y1, M1, D1, H1, N1, S1, the first date.
!
!    Input, integer ( kind = 4 ) Y2, M2, D2, H2, N2, S2, the second date.
!
!    Output, integer ( kind = 4 ) NDAY, NHOUR, NMINUTE, NSECOND, the number of 
!    days, hours, minutes, seconds between the two dates.
!
!    Output, integer ( kind = 4 ) IERROR, is 1 if either date is illegal, 
!    0 otherwise.
!
  implicit none

  character cmp
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) h1
  integer ( kind = 4 ) h2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) month_length_common
  integer ( kind = 4 ) nday
  integer ( kind = 4 ) nhour
  integer ( kind = 4 ) nminute
  integer ( kind = 4 ) nsecond
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) year_length_common
!
!  Compare the dates.
!
  call ymdhms_compare ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, n2, s2, &
    cmp, ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  We swap dates, if necessary, so that date 1 is never greater
!  than date 2.
!
  if ( cmp == '>' ) then
    write ( *, * ) 'Swapping dates'
    call ymdhms_swap ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, n2, s2 )
  end if
!
!  Set beginning to 00:00:00 1 Jan Y1.
!      end       to 00:00:00 1 Jan Y1.
!
  nhour   = h2 - h1
  nminute = n2 - n1
  nsecond = s2 - s1
  nday    = 0
!
!  Subtract seconds, minutes, hours
!
  do while ( nsecond < 0 ) 
    nsecond = nsecond + 60
    nminute = nminute - 1
  end do
 
  do while ( nminute < 0 )
    nminute = nminute + 60
    nhour = nhour - 1
  end do
 
  do while ( nhour < 0 ) 
    nhour = nhour + 24
    nday = nday - 1
  end do
!
!  Set beginning to H1:M1:S1 1 Jan Y1
!      end          H2:M2:S2 1 Jan Y2.
!
  do y = y1, y2-1
    nday = nday + year_length_common ( y )
  end do
!
!  Set beginning to H1:M1:S1 1 Jan Y1
!      end          H2:M2:S2 1 M2  Y2.
!
  do m = 1, m2-1
    nday = nday + month_length_common ( y2, m )
  end do
!
!  Set beginning to H1:M1:S1 1  Jan Y1
!      end          H2:M2:S2 D2 M2  Y2.
!
  nday = nday + ( d2 - 1 )
!
!  Set beginning to H1:M1:S1 1  M1 Y1
!      end          H2:M2:S2 D2 M2 Y2.
!
  do m = 1, m1-1
    nday = nday - month_length_common ( y1, m )
  end do
!
!  Set beginning to H1:M1:S1 D1 M1 Y1
!      end          H2:M2:S2 D2 M2 Y2.
!
  nday = nday - ( d1 - 1 )
!
!  If we swapped dates, then the differences are
!  correct, but the signs should be negated.
!
!  Set beginning to H2:M2:S2 D2 M2 Y2,
!    end       to H1:M1:S1 D1 M1 Y1.
!
  if ( cmp == '>' ) then

    nday = - nday
    nhour = - nhour
    nminute = - nminute
    nsecond = - nsecond

    call ymdhms_swap ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, n2, s2 )

  end if
 
  return
end
subroutine ymdhms_swap ( y1, m1, d1, h1, n1, s1, y2, m2, d2, h2, n2, s2 )

!*****************************************************************************80
!
!! YMDHMS_SWAP swaps the data defining two YMDHMS dates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) Y1, M1, D1, H1, M1, S1, Y2, M2, D2, H2, M2, S2;
!    the data representing the two dates.
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) h1
  integer ( kind = 4 ) h2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  call i4_swap ( y1, y2 )

  call i4_swap ( m1, m2 )
 
  call i4_swap ( d1, d2 )
 
  call i4_swap ( h1, h2 )
 
  call i4_swap ( n1, n2 )
 
  call i4_swap ( s1, s2 )
 
  return
end
