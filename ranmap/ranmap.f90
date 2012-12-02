program main

!*****************************************************************************80
!
!! MAIN is the main program for RANMAP.
!
!  Discussion:
!
!    RANMAP carries out a weighted affine mapping iteration.
!
!    The map coefficients for the dragon are "slightly" wrong, and
!    those for the tree seem very wrong.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Scott Bailey, Theodore Kim, Robert Strichartz,
!    Inside the Levy Dragon,
!    American Mathematical Monthly,
!    Volume 109, Number 8, October 2002, pages 689-703.
!
!    Michael Barnsley, Alan Sloan,
!    A Better Way to Compress Images,
!    Byte Magazine,
!    Volume 13, Number 1, January 1988, pages 215-224.
!
!    Michael Barnsley,
!    Fractals Everywhere,
!    Academic Press, 1988,
!    ISBN: 0120790696,
!    LC: QA614.86.B37.
!
!    Michael Barnsley, Lyman Hurd,
!    Fractal Image Compression,
!    Peters, 1993,
!    ISBN: 1568810008,
!    LC: TA1632.B353
!
!    Alexander Dewdney,
!    Mathematical Recreations,
!    Scientific American,
!    Volume 262, Number 5, May 1990, pages 126-129.
!
!    Bernt Wahl, Peter VanRoy, Michael Larsen, Eric Kampman,
!    Exploring Fractals on the Mac,
!    Addison Wesley, 1995,
!    ISBN: 0201626306,
!    LC: QA614.86.W34.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable, dimension (:,:) :: coord
  real ( kind = 8 ), dimension ( 2 ) :: coord_max = (/ 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension ( 2 ) :: coord_min = (/ 0.0D+00, 0.0D+00 /)
  character ( len = 80 ) :: filedot_name = 'ranmap.eps'
  character ( len = 80 ) :: filepoint_name = 'ranmap.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) map
  real ( kind = 8 ), allocatable, dimension ( : ) :: map_bin
  real ( kind = 8 ), allocatable, dimension (:,:,:) :: map_matrix
  real ( kind = 8 ), allocatable, dimension ( : ) :: map_weight
  integer ( kind = 4 ) n
  integer ( kind = 4 ) map_num
  real ( kind = 8 ) r8_uniform_01
  logical r8vec_ascends
  real ( kind = 8 ) s
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) weight_sum
  real ( kind = 8 ) x_diff
  real ( kind = 8 ) y_diff

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANMAP'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the use of weighted affine'
  write ( *, '(a)' ) '  mappings to draw fractals.'
!
!  Get the number of points to plot.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Choose the number of points to plot'
  read ( *, * ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of points to plot is ', n
!
!  Get the map option.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Choose a map:'
  write ( *, '(a)' ) '  1 = the cross'
  write ( *, '(a)' ) '  2 = the dragon'
  write ( *, '(a)' ) '  3 = the fern'
  write ( *, '(a)' ) '  4 = the leaf'
  write ( *, '(a)' ) '  5 = the Levy dragon'
  write ( *, '(a)' ) '  6 = the tree'
  write ( *, '(a)' ) '  7 = the triangle'
  write ( *, '(a)' ) '  8 = a random map (not debugged)'
  write ( *, '(a)' ) '  9 = specify your own.'

  read ( *, * ) map

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The map choice is ', map
!
!  Set the number of maps.
!
  if ( map == 1 ) then
    map_num = 5
  else if ( map == 2 ) then
    map_num = 2
  else if ( map == 3 ) then
    map_num = 4
  else if ( map == 4 ) then
    map_num = 4
  else if ( map == 5 ) then
    map_num = 2
  else if ( map == 6 ) then
    map_num = 4
  else if ( map == 7 ) then
    map_num = 3
  else if ( map == 8 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Choosing the number of maps at random.'
    map_num = i4_uniform ( 2, 6, seed )
  else if ( map == 9 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'How many maps will be used?'
    read ( *, * ) map_num
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of maps to use is ', map_num
!
!  Set the map bins.
!
  allocate ( map_bin(1:map_num-1) )

  if ( map == 1 ) then
    map_bin(1:map_num-1) = (/ 0.20D+00, 0.40D+00, 0.60D+00, 0.80D+00 /)
  else if ( map == 2 ) then
    map_bin(1:map_num-1) = (/ 0.50D+00 /)
  else if ( map == 3 ) then
    map_bin(1:map_num-1) = (/ 0.75D+00, 0.85D+00, 0.95D+00 /)
  else if ( map == 4 ) then
    map_bin(1:map_num-1) = (/ 0.25D+00, 0.50D+00, 0.75D+00 /)
  else if ( map == 5 ) then
    map_bin(1:map_num-1) = (/ 0.50D+00 /)
  else if ( map == 6 ) then
    map_bin(1:map_num-1) = (/ 0.05D+00, 0.20D+00, 0.60D+00 /)
  else if ( map == 7 ) then
    map_bin(1:map_num-1) = (/ 0.33D+00, 0.66D+00 /)
  else if ( map == 8 ) then

    allocate ( map_weight(1:map_num) )

    do i = 1, map_num
      map_weight(i) = r8_uniform_01 ( seed )
    end do

    weight_sum = sum ( map_weight(1:map_num) )
    map_weight(1:map_num) = map_weight(1:map_num) / weight_sum

    s = 0.0D+00
    do i = 1, map_num-1
      s = s + map_weight(i)
      map_bin(i) = s
    end do

    deallocate ( map_weight )

  else if ( map == 9 ) then

    do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6,a)' ) '  Enter an increasing set of ', &
        map_num-1, ' values.'
      write ( *, '(a)' ) '  These, along with 0 and 1, define the intervals'
      write ( *, '(a)' ) '  associated with each map.'
      write ( *, '(a)' ) ' '
      do i = 1, map_num-1
        read ( *, * ) map_bin(i)
      end do

      if ( .not. r8vec_ascends ( map_num-1, map_bin ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Your bin values are not acceptable because'
        write ( *, '(a)' ) '  they are not ascending.'
        cycle
      end if

      if ( map_bin(1) < 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Your bin values are not acceptable because'
        write ( *, '(a)' ) '  MAP_BIN(1) < 0.'
        cycle
      end if

      if ( 1.0D+00 < map_bin(map_num-1) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Your bin values are not acceptable because'
        write ( *, '(a)' ) '  1 < MAP_BIN(MAP_NUM-1).'
        cycle
      end if

      exit

    end do

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The map bin values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Left    Right  Length'
  write ( *, '(a)' ) ' '

  do i = 1, map_num
    if ( i == 1 ) then
      a = 0.0D+00
    else
      a = map_bin(i-1)
    end if
    if ( i < map_num ) then
      b = map_bin(i)
    else
      b = 1.0D+00
    end if
    write ( *, '(3f8.4)' ) a, b, b - a
  end do
!
!  Set the map coefficients.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Setting the map coefficients.'
  write ( *, '(a)' ) ' '

  allocate ( map_matrix(2,3,map_num) )

  if ( map == 1 ) then
    call cross_map ( map_matrix )
  else if ( map == 2 ) then
    call dragon_map ( map_matrix )
  else if ( map == 3 ) then
    call fern_map ( map_matrix )
  else if ( map == 4 ) then
    call leaf_map ( map_matrix )
  else if ( map == 5 ) then
    call levy_dragon_map ( map_matrix )
  else if ( map == 6 ) then
    call tree_map ( map_matrix )
  else if ( map == 7 ) then
    call triangle_map ( map_matrix )
  else if ( map == 8 ) then

    do k = 1, map_num
      do j = 1, 3
        do i = 1, 2
          map_matrix(i,j,k) = r8_uniform_01 ( seed )
        end do
      end do
    end do

  else if ( map == 9 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  For each map, enter 2 rows of 3 values.'
    write ( *, '(a)' ) ' '
    do k = 1, map_num
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  For map K = ', k
      write ( *, '(a)' ) ' '
      do i = 1, 2
        read ( *, * ) map_matrix(i,1:3,k)
      end do
    end do

  end if

  do k = 1, map_num
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) ' Map # ', k
    write ( *, '(a)' ) ' '
    do i = 1, 2
      write ( *, '(3f8.4)' ) map_matrix(i,1:3,k)
    end do
  end do

  allocate ( coord(2,0:n) )
!
!  Open a file, carry out the mapping, and close the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = filepoint_name, status = 'replace' )

  call random_map ( iunit, n, map_num, map_bin, map_matrix, seed, coord )

  close ( unit = iunit )

  deallocate ( map_bin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Points 0 through 10:'
  write ( *, '(a)' ) ' '
  do j = 0, 10
    write ( *, '(2x,i6,2x,f8.4,2x,f8.4)' ) j, coord(1:2,j)
  end do
!
!  Make an EPS file of the results.
!
  coord_min(1) = min ( coord_min(1), minval ( coord(1,0:n) ) )
  coord_max(1) = max ( coord_max(1), maxval ( coord(1,0:n) ) )
  coord_min(2) = min ( coord_min(2), minval ( coord(2,0:n) ) )
  coord_max(2) = max ( coord_max(2), maxval ( coord(2,0:n) ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X coordinate range:'
  write ( *, '(2x,2g14.6)' ) coord_min(1), coord_max(1)
  write ( *, '(a)' )  ' '
  write ( *, '(a)' ) '  Y coordinate range:'
  write ( *, '(2x,2g14.6)' ) coord_min(2), coord_max(2)
!
!  Preserve the aspect ratio by extending one of the ranges.
!
  x_diff = coord_max(1) - coord_min(1)
  y_diff = coord_max(2) - coord_min(2)

  if ( y_diff < x_diff ) then
    coord_min(2) = ( coord_min(2) + coord_max(2) ) - 0.5D+00 * x_diff
    coord_max(2) = ( coord_min(2) + coord_max(2) ) + 0.5D+00 * x_diff
  else if ( x_diff < y_diff ) then
    coord_min(1) = ( coord_min(1) + coord_max(1) ) - 0.5D+00 * y_diff
    coord_max(1) = ( coord_min(1) + coord_max(1) ) + 0.5D+00 * y_diff
  end if

  call dot_plot ( filedot_name, n+1, coord, coord_min, coord_max )

  deallocate ( coord )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The point data has been saved in the file "' // &
    trim ( filepoint_name ) // '".'
  write ( *, '(a)' ) '  A point plot has been saved in the file "' // &
    trim ( filedot_name ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANMAP:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
subroutine cross_map ( map_matrix )

!*****************************************************************************80
!
!! CROSS_MAP sets the cross map matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MAP_MATRIX(2,3,MAP_NUM), the set of MAP_NUM
!    affine maps.  Each map has the form of a 2 by 3 matrix, so that
!      Xnew = A11 * X + A12 * Y + A13
!      Ynew = A21 * X + A22 * Y + A23
!
  implicit none

  real ( kind = 8 ) map_matrix(2,3,5)

  map_matrix(1:2,1:3,1:5) = reshape ( source = &
    (/ 0.333D+00, 0.000D+00, &
       0.000D+00, 0.333D+00, &
       0.333D+00, 0.000D+00, &
       0.333D+00, 0.000D+00, &
       0.000D+00, 0.333D+00, &
       0.000D+00, 0.333D+00, &
       0.333D+00, 0.000D+00, &
       0.000D+00, 0.333D+00, &
       0.333D+00, 0.333D+00, &
       0.333D+00, 0.000D+00, &
       0.000D+00, 0.333D+00, &
       0.666D+00, 0.333D+00, &
       0.333D+00, 0.000D+00, &
       0.000D+00, 0.333D+00, &
       0.333D+00, 0.666D+00 /), shape = (/ 2, 3, 5 /) )

  return
end
subroutine dot_plot ( filedot_name, num_pts, coord, coord_min, coord_max )

!*****************************************************************************80
!
!! DOT_PLOT plots the individual points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILEDOT_NAME, the name of the dot file.
!
!    Input, integer ( kind = 4 ) NUM_PTS, the number of points.
!
!    Input, real ( kind = 8 ) COORD(2,NUM_PTS), the point coordinates.
!
!    Input, real ( kind = 8 ) COORD_MIN(2), COORD_MAX(2), the data ranges.
!
  implicit none

  integer ( kind = 4 ) num_pts

  real ( kind = 8 ) coord(2,num_pts)
  real ( kind = 8 ) coord_max(2)
  real ( kind = 8 ) coord_min(2)
  character ( len = * ) filedot_name
  integer ( kind = 4 ) filedot_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ), parameter :: x_ps_max = 576
  integer ( kind = 4 ), parameter :: x_ps_min = 36
  integer ( kind = 4 ), parameter :: y_ps_max = 576
  integer ( kind = 4 ), parameter :: y_ps_min = 36

  call get_unit ( filedot_unit )

  call ps_file_open ( filedot_name, filedot_unit, ierror )

  call eps_file_head ( filedot_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )

  call ps_page_head ( coord_min(1), coord_min(2), coord_max(1), coord_max(2) )

  do i = 1, num_pts
    call ps_mark_point ( coord(1,i), coord(2,i) )
  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( filedot_unit )

  return
end
subroutine dragon_map ( map_matrix )

!*****************************************************************************80
!
!! DRAGON_MAP sets the dragon map matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MAP_MATRIX(2,3,MAP_NUM), the set of MAP_NUM
!    affine maps.  Each map has the form of a 2 by 3 matrix, so that
!      Xnew = A11 * X + A12 * Y + A13
!      Ynew = A21 * X + A22 * Y + A23
!
  implicit none

  real ( kind = 8 ) map_matrix(2,3,2)

  map_matrix(1:2,1:3,1:2) = reshape ( source = &
    (/ 0.500D+00, -0.500D+00, &
       0.500D+00,  0.500D+00, &
       0.125D+00,  0.675D+00, &
       0.500D+00, -0.500D+00, &
       0.500D+00,  0.500D+00, &
      -0.125D+00,  0.325D+00 /), shape = (/ 2, 3, 2 /) )

  return
end
subroutine eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
  y_ps_max )

!*****************************************************************************80
!
!! EPS_FILE_HEAD writes header information to an encapsulated PostScript file.
!
!  Discussion:
!
!    The file should contain the description of only one page, but this
!    is not currently checked.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) X_PS_MIN, Y_PS_MIN, X_PS_MAX, Y_PS_MAX, the minimum
!    and maximum X and Y values of the data, in PostScript units.  Any data
!    that lies outside this range will not show up properly.  A reasonable
!    set of values might be 0, 0, 612, 792, or, for a half inch margin,
!    36, 36, 576, 756.
!
  implicit none

  character ( len = 8 )  date
  character ( len = * )  file_name
  real      ( kind = 8 ) line_blue
  real      ( kind = 8 ) line_green
  real      ( kind = 8 ) line_red
  integer   ( kind = 4 ) state
  integer   ( kind = 4 ) unit
  integer   ( kind = 4 ) x_ps_max
  integer   ( kind = 4 ) x_ps_min
  integer   ( kind = 4 ) y_ps_max
  integer   ( kind = 4 ) y_ps_min
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_HEAD - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 1 is required.'
    return
  end if
!
!  Initialization
!
  call ps_default ( )
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call date_and_time ( date )
!
!  Write the prolog.
!
  write ( unit, '(a)' )     '%!PS-Adobe-3.0 EPSF-3.0'
  write ( unit, '(a)' )     '%%Creator: ps_write.f90'
  write ( unit, '(a)' )     '%%Title: ' // trim ( file_name )
  write ( unit, '(a)' )     '%%CreationDate: '// trim ( date )
  write ( unit, '(a)' )     '%%Pages: 1'
  write ( unit, '(a,4i6)' ) '%%BoundingBox:', &
    x_ps_min, y_ps_min, x_ps_max, y_ps_max
  write ( unit, '(a)' )     '%%Document-Fonts: Times-Roman'
  write ( unit, '(a)' )     '%%LanguageLevel: 1'
  write ( unit, '(a)' )     '%%EndComments'
  write ( unit, '(a)' )     '%%BeginProlog'
  write ( unit, '(a)' )     '/inch {72 mul} def'
  write ( unit, '(a)' )     '%%EndProlog'
!
!  Set the font.
!
  write ( unit, '(a)' ) '/Times-Roman findfont'
  write ( unit, '(a)' ) '1.00 inch scalefont'
  write ( unit, '(a)' ) 'setfont'
!
!  Set the line color.
!
  line_red = 0.0D+00
  line_green = 0.0D+00
  line_blue = 0.0D+00

  call ps_color_line ( 'SET', line_red, line_green, line_blue )
!
!  Reset the state.
!
  state = 2

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine eps_file_tail ( )

!*****************************************************************************80
!
!! EPS_FILE_TAIL writes trailer information to an encapsulated PostScript file.
!
!  Discussion:
!
!    Looks like that penultimate 'end' line is not wanted, so I commented
!    it out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) num_pages
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state == 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_TAIL - Warning!'
    write ( *, '(a)' ) '  A page was open.  It is being forced closed.'
    state = 2
    call ps_setting_int ( 'SET', 'STATE', state )
  end if

  if ( state /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_TAIL - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 is required.'
    return
  end if
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Retrieve the number of pages.
!
  call ps_setting_int ( 'GET', 'NUM_PAGES', num_pages )

  if ( 1 < num_pages ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_TAIL - Warning!'
    write ( *, '(a)' ) '  An encapsulated PostScript file describes ONE page.'
    write ( *, '(a,i9,a)' ) '  This file describes ', num_pages, ' pages.'
    write ( *, '(a)' ) '  It is not a legal EPS file.'
  end if
!
!  Write the epilog.
!
  write ( unit, '(a)' ) '%%Trailer'
! write ( unit, '(a)' ) 'end'
  write ( unit, '(a)' ) '%%EOF'
!
!  Zero out the number of pages.
!
  num_pages = 0

  call ps_setting_int ( 'SET', 'NUM_PAGES', num_pages )
!
!  Reset the state.
!
  state = 4

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine fern_map ( map_matrix )

!*****************************************************************************80
!
!! FERN_MAP sets the fern map matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MAP_MATRIX(2,3,MAP_NUM), the set of MAP_NUM
!    affine maps.  Each map has the form of a 2 by 3 matrix, so that
!      Xnew = A11 * X + A12 * Y + A13
!      Ynew = A21 * X + A22 * Y + A23
!
  implicit none

  real ( kind = 8 ) map_matrix(2,3,4)

  map_matrix(1:2,1:3,1:4) = reshape ( source = &
    (/ 0.8560D+00, -0.0205D+00, &
       0.0414D+00,  0.8580D+00, &
       0.0700D+00,  0.1470D+00, &
       0.2440D+00,  0.1760D+00, &
      -0.3850D+00,  0.2240D+00, &
       0.3930D+00,  0.1020D+00, &
      -0.1440D+00,  0.1810D+00, &
       0.3900D+00,  0.2590D+00, &
       0.5270D+00, -0.0140D+00, &
       0.0000D+00,  0.0310D+00, &
       0.0000D+00,  0.2160D+00, &
       0.4860D+00,  0.0500D+00 /)  , shape = (/ 2, 3, 4 /) )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer ( kind = 4 ) between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer ( kind = 4 ) between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
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
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer ( kind = 4 ) between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine leaf_map ( map_matrix )

!*****************************************************************************80
!
!! LEAF_MAP sets the leaf map matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MAP_MATRIX(2,3,MAP_NUM), the set of MAP_NUM
!    affine maps.  Each map has the form of a 2 by 3 matrix, so that
!      Xnew = A11 * X + A12 * Y + A13
!      Ynew = A21 * X + A22 * Y + A23
!
  implicit none

  real ( kind = 8 ) map_matrix(2,3,4)

  map_matrix(1:2,1:3,1:4) = reshape ( source = &
    (/ 0.800D+00,  0.000D+00, &
       0.000D+00,  0.800D+00, &
       0.100D+00,  0.040D+00, &
       0.500D+00,  0.000D+00, &
       0.000D+00,  0.500D+00, &
       0.250D+00,  0.400D+00, &
       0.355D+00,  0.355D+00, &
      -0.355D+00,  0.355D+00, &
       0.266D+00,  0.078D+00, &
       0.355D+00, -0.355D+00, &
       0.355D+00,  0.355D+00, &
       0.378D+00,  0.434D+00 /), shape = (/ 2, 3, 4 /) )

  return
end
subroutine levy_dragon_map ( map_matrix )

!*****************************************************************************80
!
!! LEVY_DRAGON_MAP sets the Levy dragon map matrix.
!
!  Discussion:
!
!    You can figure out the mapping based on the following idea:
!
!    Start with a triangle T with coordinates (-1,0), (1,0), (0,1).
!
!    There are two mappings, whose images of T are:
!
!    M1: ( 0,1), (1,0), ( 1,1)
!    M2: (-1,0), (0,1), (-1,1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Scott Bailey, Theodore Kim, Robert Strichartz,
!    Inside the Levy Dragon,
!    American Mathematical Monthly,
!    Volume 109, Number 8,
!    October 2002, pages 689-703.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MAP_MATRIX(2,3,MAP_NUM), the set of MAP_NUM
!    affine maps.  Each map has the form of a 2 by 3 matrix, so that
!      Xnew = A11 * X + A12 * Y + A13
!      Ynew = A21 * X + A22 * Y + A23
!
  implicit none

  real ( kind = 8 ) map_matrix(2,3,2)

  map_matrix(1:2,1:3,1:2) = reshape ( source = &
    (/ 0.500D+00, -0.500D+00, &
       0.500D+00,  0.500D+00, &
       0.500D+00,  0.500D+00, &
       0.500D+00,  0.500D+00, &
      -0.500D+00,  0.500D+00, &
      -0.500D+00,  0.500D+00 /), shape = (/ 2, 3, 2 /) )

  return
end
function point_inside_box_2d ( x1, y1, x2, y2, x, y )

!*****************************************************************************80
!
!! POINT_INSIDE_BOX_2D determines if a point is inside a box in 2D.
!
!  Definition:
!
!    A "box" is defined by its "left down" corner and its
!    "right up" corner, and all the points between.  It is
!    assumed that the sides of the box align with coordinate directions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the two corners of the box.
!
!    Input, real ( kind = 8 ) X, Y, the point to be checked.
!
!    Output, logical POINT_INSIDE_BOX_2D, is .TRUE. if (X,Y) is inside the
!    box, or on its boundary, and .FALSE. otherwise.
!
  implicit none

  logical point_inside_box_2d
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  if ( x1 <= x .and. x <= x2 .and. &
       y1 <= y .and. y <= y2 ) then
    point_inside_box_2d = .true.
  else
    point_inside_box_2d = .false.
  end if

  return
end
subroutine ps_color_line ( action, r, g, b )

!*****************************************************************************80
!
!! PS_COLOR_LINE handles the line color.
!
!  Discussion:
!
!    By calling this routine, you can temporarily set the line color,
!    draw some lines, and then restore it to whatever it was.
!
!    An earlier version of this routine did not use the SAVE command for
!    the stack arrrays, meaning the stored data was lost.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the desired action.
!    'SET', set the line color to RGB.
!    'GET', set RGB to the current line color.
!    'PUSH', push a value onto the RGB stack.
!    'POP', pop the RGB stack.
!
!    Input, real ( kind = 8 ) R, G, B, the RGB values for the new line color.
!
  implicit none

  integer ( kind = 4 ), parameter :: nstack = 10

  character ( len = * ) action
  real ( kind = 8 ) b
  real ( kind = 8 ) b_old
  real ( kind = 8 ), save, dimension ( nstack) :: b_stack
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ), save, dimension ( nstack) :: g_stack
  integer ( kind = 4 ), save :: istack = 0
  real ( kind = 8 ) r
  real ( kind = 8 ) r_old
  real ( kind = 8 ), save, dimension ( nstack) :: r_stack
  logical s_eqi

  if ( s_eqi ( action, 'SET' ) ) then

    call ps_color_line_set ( r, g, b )

  else if ( s_eqi ( action, 'GET' ) ) then

    call ps_setting_real ( 'GET', 'LINE_RED', r )
    call ps_setting_real ( 'GET', 'LINE_GREEN', g )
    call ps_setting_real ( 'GET', 'LINE_BLUE', b )

  else if ( s_eqi ( action, 'POP' ) ) then

    if ( 0 < istack ) then
      r = r_stack(istack)
      g = g_stack(istack)
      b = b_stack(istack)
      istack = istack - 1
    end if

    call ps_color_line_set ( r, g, b )

  else if ( s_eqi ( action, 'PUSH' ) ) then

    call ps_setting_real ( 'GET', 'LINE_RED', r_old )
    call ps_setting_real ( 'GET', 'LINE_GREEN', g_old )
    call ps_setting_real ( 'GET', 'LINE_BLUE', b_old )

    if ( istack <= nstack ) then
      istack = istack + 1
      r_stack(istack) = r_old
      g_stack(istack) = g_old
      b_stack(istack) = b_old
    end if

    call ps_color_line_set ( r, g, b )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_COLOR_LINE - Fatal error!'
    write ( *, '(a)' ) '  Unexpected ACTION.'
    stop

  end if

  return
end
subroutine ps_color_line_set ( r, g, b )

!*****************************************************************************80
!
!! PS_COLOR_LINE_SET sets the line color.
!
!  Discussion:
!
!    By calling this routine, you guarantee that a check will be made
!    of the current line color.  If the current and new line colors are
!    the same, then we skip the extraneous action of setting the color.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB values for the new line color.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b_old
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ) r
  real ( kind = 8 ) r_old
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Check the state.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_COLOR_LINE_SET - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  A PostScript state of at least 1 is required.'
    return
  end if
!
!  Get the current colors.
!
  call ps_setting_real ( 'GET', 'LINE_RED', r_old )
  call ps_setting_real ( 'GET', 'LINE_GREEN', g_old )
  call ps_setting_real ( 'GET', 'LINE_BLUE', b_old )
!
!  If any color has changed, we need to reset them.
!
  if ( r_old /= r .or. g_old /= g .or. b_old /= b ) then

    call ps_setting_int ( 'GET', 'UNIT', unit )

    call ps_comment ( 'Set RGB line color.' )

    write ( unit, '(3f7.4,a)' ) r, g, b, ' setrgbcolor'

    call ps_setting_real ( 'SET', 'LINE_RED', r )
    call ps_setting_real ( 'SET', 'LINE_GREEN', g )
    call ps_setting_real ( 'SET', 'LINE_BLUE', b )

  end if

  return
end
subroutine ps_comment ( string )

!*****************************************************************************80
!
!! PS_COMMENT inserts a comment into the PostScript file.
!
!  Discussion:
!
!    A comment begins with a percent sign in column 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the comment.
!
  implicit none

  character ( len = * ) string
  integer ( kind = 4 ) unit
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Write the comment.
!
  if ( len_trim ( string ) == 0 ) then
    write ( unit, '(a)' ) '%'
  else
    write ( unit, '(a)' ) '%'
    write ( unit, '(a2,a)' ) '% ', trim ( string )
    write ( unit, '(a)' ) '%'
  end if

  return
end
subroutine ps_default ( )

!*****************************************************************************80
!
!! PS_DEFAULT sets the internal settings to their default values
!
!  Discussion:
!
!    Certain variables are not reset, including the number of pages,
!    the unit number, the internal state, and variables relating to
!    the size and shape of the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  real ( kind = 8 ) fill_blue
  real ( kind = 8 ) fill_green
  real ( kind = 8 ) fill_red
  real ( kind = 8 ) font_size
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer ( kind = 4 ) line_width
  integer ( kind = 4 ) marker_size

  line_width = 1
  marker_size = 5

  call ps_setting_int ( 'SET', 'LINE_WIDTH', line_width )
  call ps_setting_int ( 'SET', 'MARKER_SIZE', marker_size )

  fill_blue = 0.7D+00
  fill_green = 0.7D+00
  fill_red = 0.7D+00
  font_size = 0.1D+00
  line_blue = 0.0D+00
  line_green = 0.0D+00
  line_red = 0.0D+00

  call ps_setting_real ( 'SET', 'FILL_BLUE', fill_blue )
  call ps_setting_real ( 'SET', 'FILL_GREEN', fill_green )
  call ps_setting_real ( 'SET', 'FILL_RED', fill_red )
  call ps_setting_real ( 'SET', 'FONT_SIZE', font_size )
  call ps_setting_real ( 'SET', 'LINE_BLUE', line_blue )
  call ps_setting_real ( 'SET', 'LINE_GREEN', line_green )
  call ps_setting_real ( 'SET', 'LINE_RED', line_red )

  return
end
subroutine ps_file_open ( file_name, unit, ierror )

!*****************************************************************************80
!
!! PS_FILE_OPEN opens a new version of a PostScript file with a given name.
!
!  Discussion:
!
!    If a file of the given name already exists, it is deleted.
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
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) UNIT, the FORTRAN unit to which output should
!    be written.
!
!    Input, character ( len = 80 ) FILE_NAME, the name of the output file.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    nonzero, the file could not be created.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_OPEN - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 0 is required.'
    write ( *, '(a)' ) '  Call PS_FILE_CLOSE first!'
    return
  end if

  ierror = 0
!
!  Now create a new empty file of the given name.
!
  open ( unit = unit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    ierror = ios
    return
  end if

  state = 1
  call ps_setting_int ( 'SET', 'STATE', state )

  call ps_setting_int ( 'SET', 'UNIT', unit )

  return
end
subroutine ps_file_close ( unit )

!*****************************************************************************80
!
!! PS_FILE_CLOSE closes a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) UNIT, the FORTRAN unit to which output was written.
!
  implicit none

  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state < 1 .or. 4 < state ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_CLOSE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 1, 2, 3 or 4 is required.'
    return
  end if

  close ( unit = unit )

  state = 0
  call ps_setting_int ( 'SET', 'STATE', state )

  unit = 0
  call ps_setting_int ( 'SET', 'UNIT', unit )

  return
end
subroutine ps_line_closed ( npoint, x, y )

!*****************************************************************************80
!
!! PS_LINE_CLOSED adds the graph of a closed line to a PostScript file.
!
!  Discussion:
!
!    A "closed" line is one in which the last point is connected back
!    to the first one.
!
!    The current point is set to the first (and logically last) point
!    in the list.
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
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOINT, the number of points in the line.
!
!    Input, real ( kind = 8 ) X(NPOINT), Y(NPOINT), the X and Y components
!    of the points.
!
  implicit none

  integer ( kind = 4 ) npoint

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x(npoint)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(npoint)
  real ( kind = 8 ) ymin
!
!  Refuse to handle fewer than 2 points.
!
  if ( npoint < 2 ) then
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LINE_CLOSED - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
!
!  Draw lines.
!
  write ( unit, '(a)' ) 'newpath'

  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  do i = 2, npoint
    px = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    py = plotymin2 + nint ( alpha * ( y(i) - ymin ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'
  end do
!
!  Add the final extra segment to the initial point.
!
  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto'
!
!  Draw the line.
!
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x(1) )
  call ps_setting_real ( 'SET', 'YCUR', y(1) )

  return
end
subroutine ps_mark_point ( x, y )

!*****************************************************************************80
!
!! PS_MARK_POINT marks a point with a tiny point.
!
!  Discussion:
!
!    The current point is set to the point.
!
!    The point is drawn with the current RGB line colors.
!
!    If the point is outside the region, the command is ignored.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of the point to mark.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) marker_size
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  logical point_inside_box_2d
  integer ( kind = 4 ) pxcen
  integer ( kind = 4 ) pycen
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MARK_POINT - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'MARKER_SIZE', marker_size )
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'XMAX', xmax )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
  call ps_setting_real ( 'GET', 'YMAX', ymax )
!
!  If the point is outside the plot box, don't draw it.
!
  if ( .not. point_inside_box_2d ( xmin, ymin, xmax, ymax, x, y ) ) then
    return
  end if

  call ps_comment ( 'Draw a point' )

  write ( unit, '(a)' ) 'newpath'

  pxcen = plotxmin2 + nint ( alpha * ( x - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y - ymin ) )

  write ( unit, '(2i6,a)' ) pxcen, pycen, ' moveto'
  write ( unit, '(2i6,a)' ) pxcen+1, pycen, ' lineto'
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_page_head ( xmin, ymin, xmax, ymax )

!*****************************************************************************80
!
!! PS_PAGE_HEAD writes header information on a new page.
!
!  Discussion:
!
!    I think an earlier version of this code, which wrote
!    "%% Page:" rather than "%%Page:" may have caused problems
!    for some interpreters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XMIN, YMIN, XMAX, YMAX, the minimum and maximum X
!    and Y values of the data to be drawn on this page.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) num_pages
  integer ( kind = 4 ) state
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer ( kind = 4 ) margin
  integer ( kind = 4 ) pagexmax
  integer ( kind = 4 ) pagexmin
  integer ( kind = 4 ) pageymax
  integer ( kind = 4 ) pageymin
  integer ( kind = 4 ) plotxmax
  integer ( kind = 4 ) plotxmin
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymax
  integer ( kind = 4 ) plotymin
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) unit
  real ( kind = 8 ) xcur
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmax2
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xmin2
  real ( kind = 8 ) xvec(4)
  real ( kind = 8 ) ycur
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymax2
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ymin2
  real ( kind = 8 ) yvec(4)
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state == 3 ) then
    state = 2
    call ps_setting_int ( 'SET', 'STATE', state )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_PAGE_HEAD - Warning!'
    write ( *, '(a)' ) '  The current open page is forced closed.'
  end if

  if ( state /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_PAGE_HEAD - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'NUM_PAGES', num_pages )

  num_pages = num_pages + 1

  call ps_setting_int ( 'SET', 'NUM_PAGES', num_pages )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a,i6,i6)' ) '%%Page: ', num_pages, num_pages
  write ( unit, '(a)' ) 'save'
!
!  Reset the state.
!
  state = 3

  call ps_setting_int ( 'SET', 'STATE', state )
!
!  Determine and store parameters.
!
  if ( xmax == xmin ) then
    xmax2 = xmax + 1.0D+00
    xmin2 = xmax - 1.0D+00
  else
    xmax2 = xmax
    xmin2 = xmin
  end if

  if ( ymax == ymin ) then
    ymax2 = ymax + 1.0D+00
    ymin2 = ymax - 1.0D+00
  else
    ymax2 = ymax
    ymin2 = ymin
  end if
!
!  Set the value of "current point".
!
  xcur = xmin
  ycur = ymin
!
!  Set the conversion factors.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  margin = 36

  plotxmax = pagexmax - margin
  plotxmin = pagexmin + margin
  plotymax = pageymax - margin
  plotymin = pageymin + margin

  alpha = min ( real ( plotxmax - plotxmin, kind = 8 ) / ( xmax2 - xmin2 ), &
                real ( plotymax - plotymin, kind = 8 ) / ( ymax2 - ymin2 ) )
!
!  Adjust PLOTXMIN and PLOTYMIN to center the image.
!
  plotxmin2 = nint ( 0.5D+00 * &
    ( real ( plotxmin + plotxmax, kind = 8 ) - alpha * ( xmax2 - xmin2 ) ) )

  plotymin2 = nint ( 0.5D+00 * &
    ( real ( plotymin + plotymax, kind = 8 ) - alpha * ( ymax2 - ymin2 ) ) )
!
!  Store data.
!
  call ps_setting_int ( 'SET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'SET', 'PYMIN', plotymin2 )

  call ps_setting_real ( 'SET', 'ALPHA', alpha )
  call ps_setting_real ( 'SET', 'XCUR', xcur )
  call ps_setting_real ( 'SET', 'XMIN', xmin )
  call ps_setting_real ( 'SET', 'XMAX', xmax )
  call ps_setting_real ( 'SET', 'YCUR', ycur )
  call ps_setting_real ( 'SET', 'YMIN', ymin )
  call ps_setting_real ( 'SET', 'YMAX', ymax )
!
!  Draw a gray border around the page.
!
  line_red = 0.9D+00
  line_green = 0.9D+00
  line_blue = 0.9D+00

  call ps_color_line ( 'PUSH', line_red, line_green, line_blue )

  call ps_comment ( 'Draw a gray border around the page.' )

  xvec(1:4) = (/ xmin, xmax, xmax, xmin /)
  yvec(1:4) = (/ ymin, ymin, ymax, ymax /)

  call ps_line_closed ( 4, xvec, yvec )

  call ps_color_line ( 'POP', line_red, line_green, line_blue )

  return
end
subroutine ps_page_tail ( )

!*****************************************************************************80
!
!! PS_PAGE_TAIL writes tail information at the end of a page.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_PAGE_TAIL - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a)' ) 'restore showpage'

  call ps_comment ( 'End of page' )
!
!  Reset the state.
!
  state = 2

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine ps_setting_int ( action, variable, value )

!*****************************************************************************80
!
!! PS_SETTING_INT sets, gets, or prints integer ( kind = 4 ) internal PS_WRITE parameters.
!
!  Discussion:
!
!    Normally, the user does not call this routine.  It is a utility
!    used by the package.
!
!    I'd like a more sophisticated pop and push.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the desired action:
!    'GET' to get the current value of VARIABLE, or
!    'POP' to return the current value and set a new value;
!    'SET' to set a new value of VARIABLE, or
!    'PUSH' to return the current value and set a new value;
!    'PRINT' to print the current value of VARIABLE.
!
!    Input, character ( len = * ) VARIABLE, the variable to get or set:
!    'LINE_WIDTH', the line width.
!      0 is the very thinnest line possible,
!      1 is more usual, 2 is thicker, and so on.
!    'MARKER_SIZE', the size of marker circles and disks, in PostScript points;
!    'NUM_PAGES', the number of pages begun or completed;
!    'PXMIN', the location of the left hand margin of the region
!       in PostScript points;
!    'PYMIN', the location of the lower margin of the region
!       in PostScript points;
!    'STATE', the current internal state,
!      0, file not open,
!      1, file open, no header written, no page open,
!      2, file open, header written, no page open,
!      3, file open, header written, page open.
!      4, file open, header written, trailer written.
!    'UNIT', the FORTRAN output unit associated with the PostScript file.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    If ACTION = 'GET', then VALUE is an output quantity, and is the
!    current internal value of the variable.
!
!    If ACTION = 'SET', then VALUE is an input quantity, and the
!    current internal value of the variable is set to this value.
!
!    If ACTION = 'PRINT', then VALUE is ignored.
!
  implicit none

  character ( len = * ) action
  integer ( kind = 4 ), save :: line_width = 1
  integer ( kind = 4 ), save :: marker_size = 0
  integer ( kind = 4 ), save :: num_pages = 0
  integer ( kind = 4 ), save :: pxmin = 0
  integer ( kind = 4 ), save :: pymin = 0
  integer ( kind = 4 ), save :: state = 0
  integer ( kind = 4 ), save :: unit = 0
  integer ( kind = 4 ) value
  character ( len = * ) variable

  if ( variable == 'LINE_WIDTH' ) then

    if ( action == 'GET' ) then
      value = line_width
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Line width, LINE_WIDTH = ', line_width
    else if ( action == 'SET' ) then
      line_width = value
    else if ( action == 'POP' ) then
      call i4_swap ( line_width, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( line_width, value )
    end if

  else if ( variable == 'MARKER_SIZE' ) then

    if ( action == 'GET' ) then
      value = marker_size
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Marker size, MARKER_SIZE = ', marker_size
    else if ( action == 'SET' ) then
      marker_size = value
    else if ( action == 'POP' ) then
      call i4_swap ( marker_size, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( marker_size, value )
    end if

  else if ( variable == 'NUM_PAGES' ) then

    if ( action == 'GET' ) then
      value = num_pages
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Number of pages, NUM_PAGES = ', num_pages
    else if ( action == 'SET' ) then
      num_pages = value
    end if

  else if ( variable == 'PXMIN' ) then

    if ( action == 'GET' ) then
      value = pxmin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'PostScript minimum X point, PXMIN = ', pxmin
    else if ( action == 'SET' ) then
      pxmin = value
    else if ( action == 'POP' ) then
      call i4_swap ( pxmin, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( pxmin, value )
    end if

  else if ( variable == 'PYMIN' ) then

    if ( action == 'GET' ) then
      value = pymin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'PostScript minimum Y point, PYMIN = ', pymin
    else if ( action == 'SET' ) then
      pymin = value
    else if ( action == 'POP' ) then
      call i4_swap ( pymin, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( pymin, value )
    end if

  else if ( variable == 'STATE' ) then

    if ( action == 'GET' ) then
      value = state
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Current internal state, STATE = ', state
    else if ( action == 'SET' ) then
      state = value
    else if ( action == 'POP' ) then
      call i4_swap ( state, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( state, value )
    end if

  else if ( variable == 'UNIT' ) then

    if ( action == 'GET' ) then
      value = unit
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Current FORTRAN unit, UNIT = ', unit
    else if ( action == 'SET' ) then
      unit = value
    else if ( action == 'POP' ) then
      call i4_swap ( unit, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( unit, value )
    end if

  end if

  return
end
subroutine ps_setting_real ( action, variable, value )

!*****************************************************************************80
!
!! PS_SETTING_REAL sets, gets, or prints real internal PS_WRITE parameters.
!
!  Discussion:
!
!    I'd like a more sophisticated pop and push.
!
!    This routine has been revised to print an error message and stop
!    if the ACTION or VARIABLE is unexpected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is either:
!    'GET' to get the current value, or
!    'POP' to return the current value and set a new one;
!    'PRINT' to print the current value, or
!    'SET' to set the current value or
!    'PUSH' to set a new value and return the current one.
!
!    Input, character ( len = * ) VARIABLE, the variable to get or set:
!    'ALPHA', the scale factor from XY user space to PostScript points;
!    'FILL_BLUE', the intensity of the blue fill color, between 0.0 and 1.0.
!    'FILL_GREEN', the intensity of the green fill color, between 0.0 and 1.0.
!    'FILL_RED', the intensity of the red fill color, between 0.0 and 1.0.
!    'FONT_SIZE', the font size, in inches.
!    'LINE_BLUE', the blue component of the line color, between 0.0 and 1.0.
!    'LINE_GREEN', the green component of the line color, between 0.0 and 1.0.
!    'LINE_RED', the red component of the line color, between 0.0 and 1.0.
!    'XCUR', the current X location.
!    'XMAX', maximum X value of the data.
!    'XMIN', minimum X value of the data.
!    'YCUR', the current Y location.
!    'YMAX', maximum Y value of the data.
!    'YMIN', minimum Y value of the data.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If ACTION = 'GET', then VALUE is an output quantity, and is the
!    current internal value of the variable.
!
!    If ACTION = 'SET', then VALUE is an input quantity, and the
!    current internal value of the variable is set to this value.
!
!    If ACTION = 'PRINT', then VALUE is ignored.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = 0.0D+00
  real ( kind = 8 ), save :: fill_blue = 0.7D+00
  real ( kind = 8 ), save :: fill_green = 0.7D+00
  real ( kind = 8 ), save :: fill_red = 0.7D+00
  real ( kind = 8 ), save :: font_size = 0.1D+00
  real ( kind = 8 ), save :: line_blue = 0.0D+00
  real ( kind = 8 ), save :: line_green = 0.0D+00
  real ( kind = 8 ), save :: line_red = 0.0D+00
  real ( kind = 8 ) value
  character ( len = * ) variable
  real ( kind = 8 ), save :: xcur = 0.0D+00
  real ( kind = 8 ), save :: xmax = 1.0D+00
  real ( kind = 8 ), save :: xmin = 0.0D+00
  real ( kind = 8 ), save :: ycur = 0.0D+00
  real ( kind = 8 ), save :: ymax = 0.0D+00
  real ( kind = 8 ), save :: ymin = 0.0D+00

  if ( variable == 'ALPHA' ) then

    if ( action == 'GET' ) then
      value = alpha
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Scale factor from user to PS, ALPHA = ', alpha
    else if ( action == 'SET' ) then
      alpha = value
    else if ( action == 'POP' ) then
      call r8_swap ( alpha, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( alpha, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'FILL_BLUE' ) then

    if ( action == 'GET' ) then
      value = fill_blue
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Blue fill RGB value, FILL_BLUE = ', fill_blue
    else if ( action == 'SET' ) then
      fill_blue = value
    else if ( action == 'POP' ) then
      call r8_swap ( fill_blue, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( fill_blue, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'FILL_GREEN' ) then

    if ( action == 'GET' ) then
      value = fill_green
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Green fill RGB value, FILL_GREEN = ', fill_green
    else if ( action == 'SET' ) then
      fill_green = value
    else if ( action == 'POP' ) then
      call r8_swap ( fill_green, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( fill_green, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'FILL_RED' ) then

    if ( action == 'GET' ) then
      value = fill_red
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'RED fill RGB value, FILL_RED = ', fill_red
    else if ( action == 'SET' ) then
      fill_red = value
    else if ( action == 'POP' ) then
      call r8_swap ( fill_red, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( fill_red, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'FONT_SIZE' ) then

    if ( action == 'GET' ) then
      value = font_size
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Font size, FONT_SIZE = ', font_size
    else if ( action == 'SET' ) then
      font_size = value
    else if ( action == 'POP' ) then
      call r8_swap ( font_size, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( font_size, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'LINE_BLUE' ) then

    if ( action == 'GET' ) then
      value = line_blue
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Blue line RGB value, LINE_BLUE = ', line_blue
    else if ( action == 'SET' ) then
      line_blue = value
    else if ( action == 'POP' ) then
      call r8_swap ( line_blue, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( line_blue, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'LINE_GREEN' ) then

    if ( action == 'GET' ) then
      value = line_green
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Green line RGB value, LINE_GREEN = ', line_green
    else if ( action == 'SET' ) then
      line_green = value
    else if ( action == 'POP' ) then
      call r8_swap ( line_green, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( line_green, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'LINE_RED' ) then

    if ( action == 'GET' ) then
      value = line_red
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Red line RGB value, LINE_RED = ', line_red
    else if ( action == 'SET' ) then
      line_red = value
    else if ( action == 'POP' ) then
      call r8_swap ( line_red, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( line_red, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'XCUR' ) then

    if ( action == 'GET' ) then
      value = xcur
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Current X location, XCUR = ', xcur
    else if ( action == 'SET' ) then
      xcur = value
    else if ( action == 'POP' ) then
      call r8_swap ( xcur, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( xcur, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'XMAX' ) then

    if ( action == 'GET' ) then
      value = xmax
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Maximum X value, XMAX = ', xmax
    else if ( action == 'SET' ) then
      xmax = value
    else if ( action == 'POP' ) then
      call r8_swap ( xmax, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( xmax, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'XMIN' ) then

    if ( action == 'GET' ) then
      value = xmin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Minimum X value, XMIN = ', xmin
    else if ( action == 'SET' ) then
      xmin = value
    else if ( action == 'POP' ) then
      call r8_swap ( xmin, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( xmin, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'YCUR' ) then

    if ( action == 'GET' ) then
      value = ycur
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Current Y location, YCUR = ', ycur
    else if ( action == 'SET' ) then
      ycur = value
    else if ( action == 'POP' ) then
      call r8_swap ( ycur, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( ycur, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'YMAX' ) then

    if ( action == 'GET' ) then
      value = ymax
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Maximum Y value, YMAX = ', ymax
    else if ( action == 'SET' ) then
      ymax = value
    else if ( action == 'POP' ) then
      call r8_swap ( ymax, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( ymax, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'YMIN' ) then

    if ( action == 'GET' ) then
      value = ymin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Minimum Y value, YMIN = ', ymin
    else if ( action == 'SET' ) then
      ymin = value
    else if ( action == 'POP' ) then
      call r8_swap ( ymin, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( ymin, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
    write ( *, '(a)' ) '  Unexpected variable!'
    stop

  end if

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r8_to_bin_uneven ( bin_num, xbin, x, bin )

!*****************************************************************************80
!
!! R8_TO_BIN_UNEVEN places X in one of several unevenly spaced bins.
!
!  Discussion:
!
!    The XBIN array is assumed to be sorted.
!
!  Example:
!
!    BIN_NUM = 5
!    XBIN(1:4) = (/ 0.0, 2.0, 8.0, 9.0 /)
!
!    so bins are
!
!    1  ( -Inf,   0 )
!    2  (    0,   2 )
!    3  (    2,   8 )
!    4  (    8,   9 )
!    5  (    9, Inf )
!
!    X   BIN
!
!   -7    1
!   -3    1
!    0    1
!    0.1  2
!    1    2
!    3    3
!    8    3
!    9.5  5
!   13    5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BIN_NUM, the number of bins.
!
!    Input, real ( kind = 8 ) XBIN(BIN_NUM-1), the dividing values for the bins.
!
!    Input, real ( kind = 8 ) X, a value to be placed in a bin.
!
!    Output, integer ( kind = 4 ) BIN, the index of the bin to which X is assigned.
!
  implicit none

  integer ( kind = 4 ) bin_num

  integer ( kind = 4 ) bin
  real ( kind = 8 ) x
  real ( kind = 8 ) xbin(bin_num-1)

  bin = 1

  do while ( bin < bin_num )

    if ( x <= xbin(bin) ) then
      return
    end if

    bin = bin + 1

  end do

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer ( kind = 4 ) arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer ( kind = 4 ),
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function r8vec_ascends ( n, x )

!*****************************************************************************80
!
!! R8VEC_ASCENDS determines if an R8VEC is (weakly) ascending.
!
!  Example:
!
!    X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
!
!    R8VEC_ASCENDS = TRUE
!
!    The sequence is not required to be strictly ascending.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, real ( kind = 8 ) X(N), the array to be examined.
!
!    Output, logical R8VEC_ASCENDS, is TRUE if the entries of X ascend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical r8vec_ascends
  real ( kind = 8 ) x(n)

  r8vec_ascends = .false.

  do i = 1, n-1
    if ( x(i+1) < x(i) ) then
      return
    end if
  end do

  r8vec_ascends = .true.

  return
end
subroutine random_map ( iunit, n, map_num, map_bin, map_matrix, seed, coord )

!*****************************************************************************80
!
!! RANDOM_MAP randomly applies one of several affine maps to a point.
!
!  Discussion:
!
!    The user specifies a set of affine maps on [0,1]x[0,1], and their
!    associated weights or selection probabilities.
!
!    This routine then repeatedly selects a map at random, and applies
!    it to the current point to get the next point.
!
!    The computed points are written to an open I/O unit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which the points should
!    be written.
!
!    Input, integer ( kind = 4 ) N, the number of times to apply the transformation.
!
!    Input, integer ( kind = 4 ) MAP_NUM, the number of affine maps.
!
!    Input, real ( kind = 8 ) MAP_BIN(MAP_NUM-1), the cumulative probabilities
!    of the maps.  The entries of MAP_BIN should be nonnegative, increasing,
!    and no greater than 1.
!    MAP_BIN(1) is the probability that map 1 will be chosen.
!    MAP_BIN(2) is the probability that map 1 or 2 will be chosen, and so on.
!    MAP_BIN(MAP_NUM) is implicitly 1.
!
!    Input, real ( kind = 8 ) MAP_MATRIX(2,3,MAP_NUM), the set of affine maps.
!    Each map has the form of a 2 by 3 matrix, so that
!      Xnew = A11 * X + A12 * Y + A13
!      Ynew = A21 * X + A22 * Y + A23
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) COORD(2,0:N), the coordinates of the points.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) map_num

  real ( kind = 8 ) coord(2,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) map
  real ( kind = 8 ) map_matrix(2,3,map_num)
  real ( kind = 8 ) map_bin(map_num-1)
  real ( kind = 8 ) p
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(3)
!
!  Pick a random starting point.
!
!  The third component is 1, to allow us to implement affine mapping
!  using matrix multiplication.
!
  x(1) = r8_uniform_01 ( seed )
  x(2) = r8_uniform_01 ( seed )
  x(3) = 1.0D+00

  coord(1:2,0) = x(1:2)
!
!  Choose a random affine map and apply it to the current point.
!
  do i = 1, n

    p = r8_uniform_01 ( seed )

    call r8_to_bin_uneven ( map_num, map_bin, p, map )

    x(1:2) = matmul ( map_matrix(1:2,1:3,map), x(1:3) )

    coord(1:2,i) = x(1:2)

    write ( iunit, '(2x,g14.6,2x,g14.6)' ) x(1:2)

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
subroutine tree_map ( map_matrix )

!*****************************************************************************80
!
!! TREE_MAP sets the tree map matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MAP_MATRIX(2,3,MAP_NUM), the set of MAP_NUM
!    affine maps.  Each map has the form of a 2 by 3 matrix, so that
!      Xnew = A11 * X + A12 * Y + A13
!      Ynew = A21 * X + A22 * Y + A23
!
  implicit none

  real ( kind = 8 ) map_matrix(2,3,4)

  map_matrix(1:2,1:3,1:4) = reshape ( source = &
    (/  0.000D+00,  0.000D+00, &
        0.000D+00,  0.500D+00, &
        0.500D+00,  0.000D+00, &
        0.100D+00,  0.000D+00, &
        0.000D+00,  0.100D+00, &
        0.450D+00,  0.150D+00, &
        0.420D+00,  0.420D+00, &
       -0.420D+00,  0.420D+00, &
        0.290D+00, -0.010D+00, &
        0.420D+00, -0.420D+00, &
        0.420D+00,  0.420D+00, &
        0.290D+00,  0.410D+00 /), shape = (/ 2, 3, 4 /) )

  return
end
subroutine triangle_map ( map_matrix )

!*****************************************************************************80
!
!! TRIANGLE_MAP sets the triangle map matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MAP_MATRIX(2,3,MAP_NUM), the set of MAP_NUM
!    affine maps.  Each map has the form of a 2 by 3 matrix, so that
!      Xnew = A11 * X + A12 * Y + A13
!      Ynew = A21 * X + A22 * Y + A23
!
  implicit none

  real ( kind = 8 ) map_matrix(2,3,3)

  map_matrix(1:2,1:3,1:3) = reshape ( source = &
    (/ 0.500D+00, 0.000D+00, &
       0.000D+00, 0.500D+00, &
       0.000D+00, 0.000D+00, &
       0.500D+00, 0.000D+00, &
       0.000D+00, 0.500D+00, &
       0.500D+00, 0.000D+00, &
       0.500D+00, 0.000D+00, &
       0.000D+00, 0.500D+00, &
       0.250D+00, 0.330D+00 /), shape = (/ 2, 3, 3 /) )

  return
end
