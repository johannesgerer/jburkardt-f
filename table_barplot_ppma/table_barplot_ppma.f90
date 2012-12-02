program main

!*****************************************************************************80
!
!! MAIN is the main program for TABLE_BARPLOT_PPMA.
!
!  Discussion:
!
!    TABLE_BARPLOT_PPMA creates a PPMA barplot from a table file.
!
!    BAR_PPMA does a sort of bar plot, given a data file.
!
!    The program has been tailored for my specific needs, so it
!    has some quirks.  
!
!    First, we assume that the data file contains N rows of data,
!    each with M items.
!
!    We are going to plot N tall skinny bars (one pixel wide!),
!    and the bar is divided into M colors, whose relative widths
!    reflect the relative sizes of the items.
!
!    The same scale is used for each bar, so some bars may not reach
!    as high as others.  A light gray background will show up, in that case.
!
!    Only the magnitude of the items is important.
!
!  Usage:
!
!    table_barplot_ppma file.txt
!
!    where
!
!      file.txt is the name of an input data file.
!
!      Output is written to file.ppma
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: b
  integer ( kind = 4 ) bar_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: g
  real ( kind = 8 ) height_max
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  character ( len = 255 ) input_file_name
  integer ( kind = 4 ) item_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) numarg
  character ( len = 255 ) output_file_name
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: r
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: values
!
!  Say hello.
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TABLE_BARPLOT_PPMA'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Simple PPMA bar graphs from input data'
  write ( *, '(a)' ) '  in the TABLE file format.'
!
!  Count the number of command-line arguments.
!
  numarg = iargc ( )
!
!  The input file name is the first optional command line argument.
!
  if ( 1 <= numarg ) then

    iarg = 1
    call getarg ( iarg, input_file_name )

  else
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Please enter the input data file name.'
    read ( *, '(a)' ) input_file_name

  end if
!
!  Count the number of columns in the file.
!  This will be the number of separate items in each bar.
!
  call file_column_count ( input_file_name, item_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The number of items in each bar will be ', item_num
!
!  Count the number of lines in the file.
!  This will be the number of columns or bars in the chart.
!
  call file_row_count ( input_file_name, bar_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The number of bars will be ', bar_num
!
!  Read the input data file.
!
  allocate ( values(1:item_num,1:bar_num) )

  call data_read ( input_file_name, item_num, bar_num, values )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TABLE_BARPLOT_PPMA - Fatal error!'
    write ( *, '(a,i8)' ) '  DATA_READ error return = ', ierror
  end if
!
!  The data is not scaled, and I want a common scale over all bars.
!
  height_max = 0.0D+00

  do j = 1, bar_num
    height_max = max ( height_max, sum ( values(1:item_num,j) ) )
  end do

  values(1:item_num,1:bar_num) = values(1:item_num,1:bar_num) / height_max

  nrow = 300
  ncol = bar_num

  allocate ( r(1:nrow,1:ncol) )
  allocate ( g(1:nrow,1:ncol) )
  allocate ( b(1:nrow,1:ncol) )

  call bar_data_to_rgb ( bar_num, item_num, transpose ( values ), &
    nrow, ncol, r, g, b )

  deallocate ( values )
!
!  Write the RGB data to the PPMA output file.
!
  call file_ext ( input_file_name, i1, i2 )

  if ( i1 /= 0 ) then
    output_file_name = input_file_name(1:i1-2) // '.ppma'
  else
    output_file_name = trim ( input_file_name ) // '.ppma'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output file name will be "' &
    // trim ( output_file_name ) // '".'

  call ppma_write ( output_file_name, ierror, nrow, ncol, r, g, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TABLE_BARPLOT_PPMA - Fatal error!'
    write ( *, '(a,i8)' ) '  PPMA_WRITE error return = ', ierror
    stop
  end if
!
!  Terminate.
!
  deallocate ( r )
  deallocate ( g )
  deallocate ( b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TABLE_BARPLOT_PPMA'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine angle_to_rgb ( angle, r, g, b )

!*****************************************************************************80
!
!! ANGLE_TO_RGB returns a color on the perimeter of the color hexagon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle in the color hexagon.
!    The sextants are defined by the following points:
!        0 degrees, 1, 0, 0, red;
!       60 degrees, 1, 1, 0, yellow;
!      120 degrees, 0, 1, 0, green;
!      180 degrees, 0, 1, 1, cyan;
!      240 degrees, 0, 0, 1, blue;
!      300 degrees, 1, 0, 1, magenta.
!
!    Output, real ( kind = 8 ) R, G, B, RGB specifications for the color 
!    that lies at the given angle, on the perimeter of the color hexagon. 
!    One value will be 1, and one value will be 0.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle2
  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ), parameter :: degrees_to_radians = &
    3.14159265D+00 / 180.0D+00
  real ( kind = 8 ) r

  angle = mod ( angle, 360.0D+00 )

  if ( angle < 0.0D+00 ) then
    angle = angle + 360.0D+00
  end if

  if ( angle <= 60.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = 1.0D+00
    g = tan ( angle2 )
    b = 0.0D+00

  else if ( angle <= 120.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = cos ( angle2 ) / sin ( angle2 )
    g = 1.0D+00
    b = 0.0D+00

  else if ( angle <= 180.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = 1.0D+00
    b = tan ( angle2 )

  else if ( angle <= 240.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = cos ( angle2 ) / sin ( angle2 )
    b = 1.0D+00

  else if ( angle <= 300.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = tan ( angle2 )
    g = 0.0D+00
    b = 1.0D+00

  else if ( angle <= 360.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = 1.0D+00
    g = 0.0D+00
    b = cos ( angle2 ) / sin ( angle2 )

  end if

  return
end
subroutine bar_data_example ( nx, ny, ytab )

!*****************************************************************************80
!
!! BAR_DATA_EXAMPLE returns some sample bar data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NX, the number of X sample points.
!
!    Input, integer ( kind = 4 ) NY, the number of Y divisions.
!
!    Output, real ( kind = 8 ) YTAB(NX,NY), a set of values with the 
!    property that, for each I (row index or "X"), the values YTAB(I,1:NY) 
!    are all nonnegative, and sum to no more than 1.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ytab(nx,ny)

  ymax = 0.0D+00

  do i = 1, nx

    x = 2.0D+00 * pi * real ( i - 1, kind = 8 ) &
    / real ( nx - 1, kind = 8 )

    do j = 1, ny
      ytab(i,j) = 1.0D+00 + sin ( real ( j, kind = 8 ) * x )
    end do

    ymax = max ( ymax, sum ( ytab(i,1:ny) ) )

  end do
  
  ytab(1:nx,1:ny) = ytab(1:nx,1:ny) / ymax

  return
end
subroutine bar_data_to_rgb ( nx, ny, ytab, nrow, ncol, r, g, b )

!*****************************************************************************80
!
!! BAR_DATA_TO_RGB makes RGB arrays from bar graph data.
!
!  Discussion:
!
!    This routine makes some fairly strict assumptions.  
!
!    First, no X data is supplied, because it is assumed that the
!    data is sampled at equal spacings, so particular X values are
!    irrelevant.
!
!    Secondly, it is assumed that the Y data is all nonnegative,
!    and scaled in such a way that the maximum sum for any particular
!    I index is 1.  Under this assumption, the program needs to do
!    no scaling of the Y data, and can determine the color of each
!    pixel simply by adding up just enough Y's to exceed the
!    imputed Y value of the pixel.
!
!    Once the RGB arrays have been set up, a plot can be made by
!    calling, for example, the routine PPMA_WRITE from PPMA_IO.
!
!  Example:
!
!    Try really hard, and you can see a bar graph in the "plot" below.
!    We have to imagine that the three symbols $, * and - represent
!    colors.  
!
!    *-----------
!    *-------*---
!    **-----***--
!    $***-***$***
!    $$*****$$$*$
!    $$$**$$$$$$$
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NX, the number of X sample points.
!
!    Input, integer ( kind = 4 ) NY, the number of Y divisions.
!
!    Input, real ( kind = 8 ) YTAB(NX,NY), a set of values with the 
!    property that, for each I (row index or "X"), the values YTAB(I,1:NY) 
!    are all nonnegative, and sum to no more than 1.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of pixels per row and
!    column of the image.
!
!    Output, integer ( kind = 4 ) R(NROW,NCOL), G(NROW,NCOL), B(NROW,NCOL), the
!    RGB arrays for the plot.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  integer ( kind = 4 ) b(nrow,ncol)
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) pixel_column
  integer ( kind = 4 ) pixel_row
  integer ( kind = 4 ) r(nrow,ncol)
  integer ( kind = 4 ) rb
  integer ( kind = 4 ) rg
  integer ( kind = 4 ) rr
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) ysum
  real ( kind = 8 ) ytab(nx,ny)
!
!  Proceed through the columns of the graph.  Each column
!  corresponds to a particular value of X, and a set of Y values.
!
  do pixel_column = 1, ncol

    x = real ( 2 * pixel_column - 1, kind = 8 ) &
      / real ( 2 * ( ncol - 1 ), kind = 8 )

    it = nint ( real ( nx * ( 2 * pixel_column - 1 ) + ncol, kind = 8 ) &
              / real ( 2 * ncol, kind = 8 ) )

    jt = 0
    ysum = 0.0D+00

    do pixel_row = 1, nrow

      y = real ( 2 * pixel_row - 1, kind = 8 ) &
        / real ( 2 * ( nrow - 1 ), kind = 8 )

      do while ( ysum < y )

        jt = jt + 1

        if ( ny < jt ) then
          exit 
        end if

        ysum = ysum + ytab(it,jt)

      end do

      if ( ny < jt ) then

        rr = 204
        rg = 204
        rb = 217

      else

        call i4_to_rgb ( jt, rr, rg, rb )

      end if

      r(nrow+1-pixel_row,pixel_column) = rr
      g(nrow+1-pixel_row,pixel_column) = rg
      b(nrow+1-pixel_row,pixel_column) = rb

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
subroutine data_read ( file_in_name, m, n, coord )

!*****************************************************************************80
!
!! DATA_READ reads generator coordinate data from a file.
!
!  Discussion:
!
!    The file is assumed to contain one record per line.
!
!    Records beginning with the '#' character are comments, and are ignored.
!    Blank lines are also ignored.
!
!    Each line that is not ignored is assumed to contain exactly (or at least)
!    M real numbers, representing the coordinates of a point.
!
!    There are assumed to be exactly (or at least) N such records.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) N, the number of points.  The program
!    will stop reading data once N values have been read.
!
!    Output, real ( kind = 8 ) COORD(M,N), the point coordinates.
!
  implicit none
!
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) coord(m,n)
  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  real ( kind = 8 ) x(m)

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( file_in_name )
    stop
  end if

  i = 0

  do while ( i < n )

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = i
      exit
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    i = i + 1

    coord(1:m,i) = x(1:m)

  end do

  close ( unit = file_in_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DATA_READ:'
  write ( *, '(a,i8)' ) '  Read coordinate data from file.'

  return
end
subroutine file_column_count ( file_name, ncolumn )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of NCOLUMN words, separated
!    by spaces.  There may also be some blank lines, and some comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NCOLUMN, the number of columns assumed to be in the file.
!
  implicit none

  character ( len = * ) file_name
  logical got_one
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 256 ) line
  integer ( kind = 4 ) ncolumn
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    ncolumn = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    got_one = .true.
    exit

  end do

  if ( .not. got_one ) then

    rewind ( iunit )

    do 

      read ( iunit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = iunit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    ncolumn = 0
    return
  end if

  call s_word_count ( line, ncolumn )

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
subroutine file_row_count ( file_name, nline )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of lines in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Blank lines and comment lines, which begin with '#', are not counted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NLINE, the number of lines found in the file.
!    If the file could not be opened, then NLINE is returned as -1.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 256 ) line
  integer ( kind = 4 ) nline
  logical, parameter :: verbose = .false.

  nline = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then

    nline = - 1

    if ( verbose ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
      write ( *, '(a)' ) '  Could not open the file:'
      write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    end if

    return

  end if
!
!  Count the lines.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    nline = nline + 1

  end do

  close ( unit = iunit )

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
!    18 September 2005
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
function i4_log_2 ( i )

!*****************************************************************************80
!
!! I4_LOG_2 returns the integer part of the logarithm base 2 of |I|.
!
!  Example:
!
!     I  I4_LOG_2
!
!     0  -1
!     1,  0
!     2,  1
!     3,  1
!     4,  2
!     5,  2
!     6,  2
!     7,  2
!     8,  3
!     9,  3
!    10,  3
!   127,  6
!   128,  7
!   129,  7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 2 is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_2, the integer part of the logarithm base 2 of
!    the absolute value of I.
!    For positive I4_LOG_2(I), it should be true that
!      2**I4_LOG_2(X) <= |I| < 2**(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
  implicit none

  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs

  if ( i == 0 ) then

    i4_log_2 = - huge ( i4_log_2 )

  else

    i4_log_2 = 0

    i_abs = abs ( i )

    do while ( 2 <= i_abs )
      i_abs = i_abs / 2
      i4_log_2 = i4_log_2 + 1
    end do

  end if

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the positive remainder when I is divided by J.
!
!  Discussion:
!
!    NREM = I4_MODP ( I, J )
!    NMULT = ( I - NREM ) / J
!
!    I = J * NMULT + NREM
!
!  Example:
!
!        I         J   NMULT  NREM    Factorization
!
!      107        50       2     7    107 =  2 *  50 + 7
!      107       -50      -2     7    107 = -2 * -50 + 7
!     -107        50      -3    43   -107 = -3 *  50 + 43
!     -107       -50       3    43   -107 =  3 * -50 + 43
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
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the positive remainder when I is divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) i4_modp

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two integers.
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
!    Input/output, integer ( kind = 4 ) I, J, the two integers to be swapped.
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
subroutine i4_to_angle ( i, angle )

!*****************************************************************************80
!
!! I4_TO_ANGLE maps integers to points on a circle.
!
!  Discussion:
!
!    The angles are intended to be used to select colors on a color
!    hexagon whose 6 vertices are red, yellow, green, cyan, blue,
!    magenta.
!
!  Example:
!
!     I   X      ANGLE
!
!     0   0/3      0
!     1   1/3    120
!     2   2/3    240
!
!     3   1/6     60
!     4   3/6    180
!     5   5/6    300
!
!     6   1/12    30
!     7   3/12    90
!     8   5/12   150
!     9   7/12   210
!    10   9/12   270
!    11  11/12   330
!
!    12   1/24    15
!    13   3/24    45
!    14   5/24    75
!    etc
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the desired color.
!
!    Output, real ( kind = 8 ) ANGLE, an angle, measured in degrees,
!    between 0 and 360.
!
  implicit none

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4

  if ( 0 <= abs ( i ) .and. abs ( i ) <= 2 ) then

    angle = 120.0D+00 * real ( abs ( i ), kind = 8 )

  else

    i1 = i4_log_2 ( abs ( i ) / 3 )
    i2 = abs ( i ) + 1 - 3 * 2**i1
    i3 = 2 * ( i2 - 1 ) + 1
    i4 = 3 * 2**( i1 + 1 )

    angle = 360.0D+00 * real ( i3, kind = 8 ) / real ( i4, kind = 8 )

  end if

  return
end
subroutine i4_to_rgb ( i, r, g, b )

!*****************************************************************************80
!
!! I4_TO_RGB maps integers to RGB colors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the desired color.
!
!    Output, integer ( kind = 4 ) R, G, B, the RGB specifications for a color.
!
  implicit none

  real ( kind = 8 ) angle
  integer ( kind = 4 ) b
  integer ( kind = 4 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) r
  real ( kind = 8 ) rb
  real ( kind = 8 ) rg
  real ( kind = 8 ) rr
!
!  Red
!
  if ( i == 1 ) then
    r = 255
    g = 0
    b = 0
!
!  Green
!
  else if ( i == 2 ) then
    r = 0
    g = 255
    b = 0
!
!  Blue
!
  else if ( i == 3 ) then
    r = 0
    g = 0
    b = 255
!
!  Cyan
!
  else if ( i == 4 ) then
    r = 0
    g = 255
    b = 255
!
!  Magenta
!
  else if ( i == 5 ) then
    r = 255
    g = 0
    b = 255
!
!  Yellow
!
  else if ( i == 6 ) then
    r = 255
    g = 255
    b = 0
!
!  Brown5
!
  else if ( i == 7 ) then
    r = 139
    g =  35
    b =  35
!
!  Orange
!
  else if ( i == 8 ) then
    r = 255
    g = 165
    b = 0
!
!  Goldenrod5
!
  else if ( i == 9 ) then
    r = 139
    g = 105
    b =  20
!
!  Medium Purple
!
  else if ( i == 10 ) then
    r = 147
    g = 112
    b = 219
!
!  Coral
!
   else if ( i == 11 ) then

    r = 255
    g = 127
    b =  80
!
!  Pink5
!
  else if ( i == 12 ) then
    r = 139
    g =  99
    b = 108
!
!  GreenYellow
!
  else if ( i == 13 ) then
    r = 173
    g = 255
    b =  47
!
!  Aquamarine
!
  else if ( i == 14 ) then
    r = 127
    g = 255
    b = 212
!
!  Pale Green3
!
  else if ( i == 15 ) then
    r = 124
    g = 205
    b = 124
!
!  Burlywood
!
  else if ( i == 16 ) then
    r = 222
    g = 184
    b = 135
!
!  Cornsilk3
!
  else if ( i == 17 ) then
    r = 205
    g = 200
    b = 177
!
!  Lemon_Chiffon3
!
  else if ( i == 18 ) then
    r = 205
    g = 201
    b = 165
!
!  Maroon
!
  else if ( i == 19 ) then
    r = 176
    g = 48
    b = 96
!
!  Slate_Blue2
!
  else if ( i == 20 ) then
    r = 131
    g = 111
    b = 255

  else

    call i4_to_angle ( i, angle )

    call angle_to_rgb ( angle, rr, rg, rb )

    r = min ( int ( rr * 255 ), 255 )
    g = min ( int ( rg * 255 ), 255 )
    b = min ( int ( rb * 255 ), 255 )

  end if

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
      write ( output_unit, '(12i6)' ) ( r(i,j), g(i,j), b(i,j), j = jlo, jhi )
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
    write ( *, '(a,i8)' ) '  Number of data rows NROW =    ', nrow
    write ( *, '(a,i8)' ) '  Number of data columns NCOL = ', ncol
    write ( *, '(a,i8)' ) '  Maximum color value MAXCOL =  ', maxcol
  end if

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

  if ( llen1 < llen2 ) then
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
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
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
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
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
!    '-4.2D+2'         -420.0
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
!    Output, real ( kind = 8 ) R, the real value that was read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!
!    0, no errors occurred.
!
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
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
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
        rtop = 10.0D+00 * rtop + real ( ndig )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig )
        rbot = 10.0D+00 * rbot
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
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
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
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) length
  real ( kind = 8 ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

  end do

  return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
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
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) nword
  character ( len = * ) s

  nword = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      nword = nword + 1
      blank = .false.
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
