function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!
  implicit none

  character c
  logical ch_is_digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
!    If C was 'illegal', then DIGIT is -1.
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
subroutine digit_inc ( c )

!*****************************************************************************80
!
!! DIGIT_INC increments a decimal digit.
!
!  Example:
!
!    Input  Output
!    -----  ------
!    '0'    '1'
!    '1'    '2'
!    ...
!    '8'    '9'
!    '9'    '0'
!    'A'    'A'
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
!    Input/output, character C, a digit to be incremented.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  call ch_to_digit ( c, digit )

  if ( digit == -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, c )

  return
end
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
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
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine edge_list ( element_num, element_node, edge_num, edge_nodes )

!*****************************************************************************80
!
!! EDGE_LIST creates a list of the unique edges in a graph.
!
!  Discussion:
!
!    The routine extracts the successive pairs of vertices that
!    define each edge of a face.  It reorders each pair so that
!    the lesser element is listed first.  It sorts the entire list.
!    Then it counts the unique entries and tosses out the repeats.
!
!    The routine does not distinguish "directions" of an edge.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes
!    making faces.
!
!    Output, integer ( kind = 4 ) EDGE_NUM, the number of unique edges.
!
!    Output, integer ( kind = 4 ) EDGE_NODES(2,3*ELEMENT_NUM), contains in
!    the first EDGE_NUM columns the unique edges.
!
  implicit none

  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) edge_nodes(2,3*element_num)
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) edge_num_old
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(3,element_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) vert
  integer ( kind = 4 ) vert2
!
!  First count the number of edges with duplication.
!
  edge_num = 3 * element_num
!
!  Store the edges.
!
  edge_num = 0
  do element = 1, element_num
    do vert = 1, 3
      edge_num = edge_num + 1
      i = element_node(vert,element)
      vert2 = i4_wrap ( vert+1, 1, 3 )
      j = element_node(vert2,element)
      edge_nodes(1,edge_num) = min ( i, j )
      edge_nodes(2,edge_num) = max ( i, j )
    end do
  end do
!
!  Sort the edges.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( edge_num, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4_swap ( edge_nodes(1,i), edge_nodes(1,j) )
      call i4_swap ( edge_nodes(2,i), edge_nodes(2,j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      if ( edge_nodes(1,i) < edge_nodes(1,j) ) then
        isgn = -1
      else if ( edge_nodes(1,i) == edge_nodes(1,j) ) then
        if ( edge_nodes(2,i) < edge_nodes(2,j) ) then
          isgn = -1
        else if ( edge_nodes(2,i) == edge_nodes(2,j) ) then
          isgn = 0
        else
          isgn = 1
        end if
      else
        isgn = 1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do
!
!  Count the unique entries and squash the array.
!
  edge_num_old = edge_num

  edge_num = 0

  do i = 1, edge_num_old

    if ( i == 1 ) then

      edge_num = 1

    else

      if ( edge_nodes(1,edge_num) /= edge_nodes(1,i) .or. &
           edge_nodes(2,edge_num) /= edge_nodes(2,i) ) then

        edge_num = edge_num + 1
        edge_nodes(1,edge_num) = edge_nodes(1,i)
        edge_nodes(2,edge_num) = edge_nodes(2,i)

      end if

    end if

  end do

  return
end
subroutine element3_eps ( file_name, node_num, node_x, node_y, element_num, &
  element_mask, element_node, title )

!*****************************************************************************80
!
!! ELEMENT3_EPS creates an EPS file containing an image of the mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to create.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_X(NODE_NUM), NODE_Y(NODE_NUM),
!    the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, logical ELEMENT_MASK(ELEMENT_NUM), a mask for the elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) ave_x
  real ( kind = 8 ) ave_y
  integer ( kind = 4 ), parameter :: circle_size = 3
  real ( kind = 8 ) dif
  integer ( kind = 4 ) element
  logical element_mask(element_num)
  integer ( kind = 4 ) element_node(3,element_num)
  integer ( kind = 4 ) eps_unit
  integer ( kind = 4 ) eps_x
  integer ( kind = 4 ) eps_y
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  logical node_mask(node_num)
  real ( kind = 8 ) node_x(node_num)
  real ( kind = 8 ) node_x_max
  real ( kind = 8 ) node_x_min
  real ( kind = 8 ) node_y(node_num)
  real ( kind = 8 ) node_y_max
  real ( kind = 8 ) node_y_min
  real ( kind = 8 ) scale
  character ( len = 40 ) string
  character ( len = * ) title
!
!  Determine the range of the unmasked elements.
!
  node_x_min =  huge ( node_x_min )
  node_x_max = -huge ( node_x_max )
  node_y_min =  huge ( node_y_min )
  node_y_max = -huge ( node_y_max )

  node_mask(1:node_num) = .false.

  do element = 1, element_num
    if ( element_mask(element) ) then
      do j = 1, 3
        node = element_node(j,element)
        node_mask(node) = .true.
        node_x_min = min ( node_x_min, node_x(node) )
        node_x_max = max ( node_x_max, node_x(node) )
        node_y_min = min ( node_y_min, node_y(node) )
        node_y_max = max ( node_y_max, node_y(node) )
      end do
    end if
  end do

  if ( node_y_max - node_y_min < node_x_max - node_x_min ) then
    scale = node_x_max - node_x_min
    dif = ( node_x_max - node_x_min ) - ( node_y_max - node_y_min )
    node_y_max = node_y_max + 0.5D+00 * dif
    node_y_min = node_y_min - 0.5D+00 * dif
  else
    scale = node_y_max - node_y_min
    dif = ( node_y_max - node_y_min ) - ( node_x_max - node_x_min )
    node_x_max = node_x_max + 0.5D+00 * dif
    node_x_min = node_x_min - 0.5D+00 * dif
  end if

  call get_unit ( eps_unit )

  open ( unit = eps_unit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ELEMENT3_EPS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output EPS file.'
    stop
  end if

  write ( eps_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( eps_unit, '(a)' ) '%%Creator: element3_eps(test_mesh.f90)'
  write ( eps_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( eps_unit, '(a)' ) '%%Pages: 1'
  write ( eps_unit, '(a)' ) '%%BoundingBox:    36    36   576   756'
  write ( eps_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( eps_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( eps_unit, '(a)' ) '%%EndComments'
  write ( eps_unit, '(a)' ) '%%BeginProlog'
  write ( eps_unit, '(a)' ) '/inch {72 mul} def'
  write ( eps_unit, '(a)' ) '%%EndProlog'
  write ( eps_unit, '(a)' ) '%%Page:      1     1'
  write ( eps_unit, '(a)' ) 'save'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Set RGB line color.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.9000 0.9000 0.9000 setrgbcolor'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Draw a gray border around the page.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'newpath'
  write ( eps_unit, '(a)' ) '    36   126 moveto'
  write ( eps_unit, '(a)' ) '   576   126 lineto'
  write ( eps_unit, '(a)' ) '   576   666 lineto'
  write ( eps_unit, '(a)' ) '    36   666 lineto'
  write ( eps_unit, '(a)' ) '    36   126 lineto'
  write ( eps_unit, '(a)' ) 'stroke'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Set RGB line color.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Label the plot:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'
  write ( eps_unit, '(a)' ) '/Times-Roman findfont 0.50 inch scalefont setfont'
  write ( eps_unit, '(a)' ) '    36   666 moveto'
  write ( eps_unit, '(a)' ) '(' // trim ( title ) // ') show'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Define a clipping polygon'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '    36   126 moveto'
  write ( eps_unit, '(a)' ) '   576   126 lineto'
  write ( eps_unit, '(a)' ) '   576   666 lineto'
  write ( eps_unit, '(a)' ) '    36   666 lineto'
  write ( eps_unit, '(a)' ) '    36   126 lineto'
  write ( eps_unit, '(a)' ) 'clip newpath'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw filled dots at each node:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.9000 setrgbcolor'

  do node = 1, node_num

    if ( node_mask(node) ) then

      eps_x = int &
        ( ( node_x_max - node_x(node)              ) *  61.0D+00   &
        + (            + node_x(node) - node_x_min ) * 551.0D+00 ) &
        / scale

      eps_y = int &
        ( ( node_y_max - node_y(node)              ) * 151.0D+00   &
        + (              node_y(node) - node_y_min ) * 641.0D+00 ) &
        / scale

      write ( eps_unit, '(a,i4,2x,i4,2x,i4,a)' ) &
        'newpath  ', eps_x, eps_y, circle_size, ' 0 360 arc closepath fill'

    end if

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Label the nodes:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 1.0000 setrgbcolor'
  write ( eps_unit, '(a)' ) '/Times-Roman findfont 0.20 inch scalefont setfont'

  do node = 1, node_num

    if ( node_mask(node) ) then

      eps_x = int &
        ( ( node_x_max - node_x(node)              ) *  61.0D+00   &
        + (            + node_x(node) - node_x_min ) * 551.0D+00 ) &
        / scale

      eps_y = int &
        ( ( node_y_max - node_y(node)              ) * 151.0D+00   &
        + (              node_y(node) - node_y_min ) * 641.0D+00 ) &
        / scale

      write ( string, '(i4)' ) node
      string = adjustl ( string )

      write ( eps_unit, '(i4,2x,i4,a)' ) eps_x, eps_y+5, &
        ' moveto (' // trim ( string ) // ') show'

    end if

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw the element sides:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.9000 0.0000 0.0000 setrgbcolor'

  do element = 1, element_num

    if ( .not. element_mask(element) ) then
      cycle
    end if

    node = element_node(1,element)

    eps_x = int &
      ( ( node_x_max - node_x(node)              ) *  61.0D+00   &
      + (            + node_x(node) - node_x_min ) * 551.0D+00 ) &
      / scale

    eps_y = int &
      ( ( node_y_max - node_y(node)              ) * 151.0D+00   &
      + (              node_y(node) - node_y_min ) * 641.0D+00 ) &
      / scale

    write ( eps_unit, '(a,i4,2x,i4,a)' ) 'newpath ', eps_x, eps_y, ' moveto'

    do i = 1, 3

      ip1 = mod ( i, 3 ) + 1;
      node = element_node(ip1,element)

      eps_x = int &
        ( ( node_x_max - node_x(node)              ) *  61.0D+00   &
        + (            + node_x(node) - node_x_min ) * 551.0D+00 ) &
        / scale

      eps_y = int &
        ( ( node_y_max - node_y(node)              ) * 151.0D+00   &
        + (              node_y(node) - node_y_min ) * 641.0D+00 ) &
        / scale

      write ( eps_unit, '(i4,2x,i4,a)' ) eps_x, eps_y, ' lineto'

    end do

    write ( eps_unit, '(a)' ) 'stroke'

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Label the elements:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 1.0000 0.0000 0.0000 setrgbcolor'
  write ( eps_unit, '(a)' ) '/Times-Roman findfont 0.30 inch scalefont setfont'

  do element = 1, element_num

    if ( .not. element_mask(element) ) then
      cycle
    end if

    ave_x = 0.0D+00
    ave_y = 0.0D+00

    do i = 1, 3

      node = element_node(i,element)

      ave_x = ave_x + node_x(node)
      ave_y = ave_y + node_y(node)

    end do

    ave_x = ave_x / 3.0D+00
    ave_y = ave_y / 3.0D+00

    eps_x = int &
      ( ( node_x_max - ave_x              ) *  61.0D+00   &
      + (            + ave_x - node_x_min ) * 551.0D+00 ) &
      / scale

    eps_y = int &
      ( ( node_y_max - ave_y              ) * 151.0D+00   &
      + (              ave_y - node_y_min ) * 641.0D+00 ) &
      / scale

    write ( string, '(i4)' ) element
    string = adjustl ( string )

    write ( eps_unit, '(i4,2x,i4,a)' ) eps_x, eps_y, ' moveto (' &
      // trim ( string ) // ') show'

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'restore showpage'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% End of page'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%%Trailer'
  write ( eps_unit, '(a)' ) '%%EOF'

  close ( unit = eps_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELEMENT3_EPS:'
  write ( *, '(a)' ) '  An encapsulated PostScript file was created'
  write ( *, '(a)' ) '  containing an image of the nodes and elements.'
  write ( *, '(a)' ) '  The file is named "' // trim ( file_name ) // '".'

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC generates the next filename in a series.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Example:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  logical ch_is_digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( ch_is_digit ( c ) ) then

      call digit_inc ( c )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

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
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer division.
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
!  Example:
!
!        I     J     MOD  I4_MODP    Factorization
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
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j

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
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an integer to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) wide

  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
    i4_wrap = ilo
  else
    i4_wrap = ilo + i4_modp ( ival - ilo, wide )
  end if

  return
end
subroutine mesh_node_matrix ( element_num, node_num, element_node, node_adj )

!*****************************************************************************80
!
!! MESH_NODE_MATRIX returns the node adjacency matrix of a mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
!    Output, integer ( kind = 4 ) NODE_ADJ(NODE_NUM,NODE_NUM), the node 
!    adjacency matrix.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) element_node(3,element_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) node_adj(node_num,node_num)

  node_adj(1:node_num,1:node_num) = 0

  do j = 1, element_num
    do i = 1, 3

      ip1 = i4_wrap ( i+1, 1, 3 )
      n1 = element_node(i,j)
      n2 = element_node(ip1,j)

      node_adj(n1,n1) = node_adj(n1,n1) + 1
      node_adj(n1,n2) = node_adj(n1,n2) - 1
      node_adj(n2,n1) = node_adj(n2,n1) - 1
      node_adj(n2,n2) = node_adj(n2,n2) + 1

    end do
  end do

  return
end
subroutine mesh_poly ( file_name, node_num, x, y, element_num, &
  element_node, edge_num, edge_nodes, hole_num, hole_x, hole_y )

!*****************************************************************************80
!
!! MESH_POLY creates a POLY file of a mesh, for input to TRIANGLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to create.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the XY
!    coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Input, integer ( kind = 4 ) EDGE_NODES(2,EDGE_NUM), the nodes that
!    form each edge.
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of holes in the mesh.
!
!    Input, real ( kind = 8 ) HOLE_X(NODE_NUM), HOLE_Y(NODE_NUM),
!    the XY coordinates of a point in each hold.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) node_num

  character ( len = 8 ) date
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_nodes(2,edge_num)
  integer ( kind = 4 ) element_node(3,element_num)
  character ( len = * ) file_name
  integer ( kind = 4 ) hole
  real ( kind = 8 ) hole_x(hole_num)
  real ( kind = 8 ) hole_y(hole_num)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) node
  integer ( kind = 4 ) poly_unit
  integer ( kind = 4 ), parameter :: region_num = 0
  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) y(node_num)

  call get_unit ( poly_unit )

  open ( unit = poly_unit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_POLY - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output POLY file.'
    stop
  end if

  call date_and_time ( date )

  write ( poly_unit, '(a)' ) '#  ' // trim ( file_name )
  write ( poly_unit, '(a)' ) '#  Created by mesh_poly(testmesh.f90)'
  write ( poly_unit, '(a)' ) '#  Creation date:' // trim ( date )
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Vertex  Dimension  Attribute  Marker'
  write ( poly_unit, '(a)' ) '#  Count              Count      0/1'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(2x,i8,a)' ) node_num, '  2  0  0'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Vertex  X  Y  Attributes  Marker'
  write ( poly_unit, '(a)' ) '#  Index'
  write ( poly_unit, '(a)' ) '#'
  do node = 1, node_num
    write ( poly_unit, '(2x,i8,2x,f10.6,2x,f10.6)' ) node, x(node), y(node)
  end do
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Segment  Marker'
  write ( poly_unit, '(a)' ) '#  Count    0/1'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(2x,i8,a)' ) edge_num, '  0'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Segment_index  Node1  Node2  Marker'
  write ( poly_unit, '(a)' ) '#'
  do edge = 1, edge_num
    write ( poly_unit, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) &
      edge, edge_nodes(1,edge), edge_nodes(2,edge), 0
  end do
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Hole'
  write ( poly_unit, '(a)' ) '#  Count'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(2x,i8)' ) hole_num
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Hole_index  X  Y'
  write ( poly_unit, '(a)' ) '#'
  do hole = 1, hole_num
    write ( poly_unit, '(2x,i8,2x,f10.6,2x,f10.6)' ) &
      hole, hole_x(hole), hole_y(hole)
  end do
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(a)' ) '#  Region'
  write ( poly_unit, '(a)' ) '#  Count'
  write ( poly_unit, '(a)' ) '#'
  write ( poly_unit, '(i8)' ) region_num

  close ( unit = poly_unit )

  return
end
subroutine mesh00_element_node ( number, element_num, element_node )

!*****************************************************************************80
!
!! MESH00_ELEMENT_NODE returns the element->node data of any mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NUMBER, the number of the mesh.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
  implicit none

  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) element_node(3,element_num)
  integer ( kind = 4 ) number

  if ( number == 1 ) then
    call mesh01_element_node ( element_num, element_node )
  else if ( number == 2 ) then
    call mesh02_element_node ( element_num, element_node )
  else if ( number == 3 ) then
    call mesh03_element_node ( element_num, element_node )
  else if ( number == 4 ) then
    call mesh04_element_node ( element_num, element_node )
  else if ( number == 5 ) then
    call mesh05_element_node ( element_num, element_node )
  else if ( number == 6 ) then
    call mesh06_element_node ( element_num, element_node )
  else if ( number == 7 ) then
    call mesh07_element_node ( element_num, element_node )
  else if ( number == 8 ) then
    call mesh08_element_node ( element_num, element_node )
  else if ( number == 9 ) then
    call mesh09_element_node ( element_num, element_node )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH00_ELEMENT_NODE - Fatal error!'
    write ( *, '(a,i8)' ) '  Unknown mesh number = ', number
    stop
  end if

  return
end
subroutine mesh00_element3_eps ( number, file_name )

!*****************************************************************************80
!
!! MESH00_ELEMENT3_EPS creates an image of any mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NUMBER, the number of the mesh.
!
!    Input, character ( len = * ) FILE_NAME, the name for the file.
!
  implicit none

  logical, allocatable, dimension ( : ) :: element_mask
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  character ( len = * ) file_name
  character ( len = 80 ) name
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: node_x
  real ( kind = 8 ), allocatable, dimension ( : ) :: node_y
  integer ( kind = 4 ) number

  call mesh00_name ( number, name )

  call mesh00_node_num ( number, node_num )

  allocate ( node_x(1:node_num) )
  allocate ( node_y(1:node_num) )

  call mesh00_element_num ( number, element_num )

  allocate ( element_mask(1:element_num) )
  allocate ( element_node(1:3,1:element_num) )

  call mesh00_node_xy ( number, node_num, node_x, node_y )

  call mesh00_element_node ( number, element_num, element_node )

  element_mask(1:element_num) = .true.

  call element3_eps ( file_name, node_num, node_x, node_y, element_num, &
    element_mask, element_node, 'The Elements:' )

  deallocate ( element_mask)
  deallocate ( element_node )
  deallocate ( node_x )
  deallocate ( node_y )

  return
end
subroutine mesh00_element_num ( number, element_num )

!*****************************************************************************80
!
!! MESH00_ELEMENT_NUM returns the number of elements of any mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NUMBER, the number of the mesh.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) number

  if ( number == 1 ) then
    call mesh01_element_num ( element_num )
  else if ( number == 2 ) then
    call mesh02_element_num ( element_num )
  else if ( number == 3 ) then
    call mesh03_element_num ( element_num )
  else if ( number == 4 ) then
    call mesh04_element_num ( element_num )
  else if ( number == 5 ) then
    call mesh05_element_num ( element_num )
  else if ( number == 6 ) then
    call mesh06_element_num ( element_num )
  else if ( number == 7 ) then
    call mesh07_element_num ( element_num )
  else if ( number == 8 ) then
    call mesh08_element_num ( element_num )
  else if ( number == 9 ) then
    call mesh09_element_num ( element_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH00_ELEMENT_NUMBER - Fatal error!'
    write ( *, '(a,i8)' ) '  Unknown mesh number = ', number
    stop
  end if

  return
end
subroutine mesh00_hole_num ( number, hole_num )

!*****************************************************************************80
!
!! MESH00_HOLE_NUM returns the number of holes in a mesh given its number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NUMBER, the number of the mesh.
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes in the mesh.
!
  implicit none

  integer ( kind = 4 ) hole_num
  integer ( kind = 4 ) number

  if ( number == 1 ) then
    call mesh01_hole_num ( hole_num )
  else if ( number == 2 ) then
    call mesh02_hole_num ( hole_num )
  else if ( number == 3 ) then
    call mesh03_hole_num ( hole_num )
  else if ( number == 4 ) then
    call mesh04_hole_num ( hole_num )
  else if ( number == 5 ) then
    call mesh05_hole_num ( hole_num )
  else if ( number == 6 ) then
    call mesh06_hole_num ( hole_num )
  else if ( number == 7 ) then
    call mesh07_hole_num ( hole_num )
  else if ( number == 8 ) then
    call mesh08_hole_num ( hole_num )
  else if ( number == 9 ) then
    call mesh09_hole_num ( hole_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH00_HOLE_NUM - Fatal error!'
    write ( *, '(a,i8)' ) '  Unknown mesh number = ', number
    stop
  end if

  return
end
subroutine mesh00_hole_xy ( number, hole_num, hole_x, hole_y )

!*****************************************************************************80
!
!! MESH00_HOLE_XY returns hole coordinates of a mesh given its number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NUMBER, the number of the mesh.
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
!    Output, real ( kind = 8 ) HOLE_X(HOLE_NUM), HOLE_Y(HOLE_NUM),
!    the coordinates of a point in each hole.
!
  implicit none

  integer ( kind = 4 ) hole_num

  integer ( kind = 4 ) number
  real ( kind = 8 ) hole_x(hole_num)
  real ( kind = 8 ) hole_y(hole_num)

  if ( number == 1 ) then
    call mesh01_hole_xy ( hole_num, hole_x,  hole_y )
  else if ( number == 2 ) then
    call mesh02_hole_xy ( hole_num, hole_x,  hole_y )
  else if ( number == 3 ) then
    call mesh03_hole_xy ( hole_num, hole_x,  hole_y )
  else if ( number == 4 ) then
    call mesh04_hole_xy ( hole_num, hole_x,  hole_y )
  else if ( number == 5 ) then
    call mesh05_hole_xy ( hole_num, hole_x,  hole_y )
  else if ( number == 6 ) then
    call mesh06_hole_xy ( hole_num, hole_x,  hole_y )
  else if ( number == 7 ) then
    call mesh07_hole_xy ( hole_num, hole_x,  hole_y )
  else if ( number == 8 ) then
    call mesh08_hole_xy ( hole_num, hole_x,  hole_y )
  else if ( number == 9 ) then
    call mesh09_hole_xy ( hole_num, hole_x,  hole_y )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH00_HOLE_XY - Fatal error!'
    write ( *, '(a,i8)' ) '  Unknown mesh number = ', number
    stop
  end if

  return
end
subroutine mesh00_name ( number, name )

!*****************************************************************************80
!
!! MESH00_NAME returns the name of a mesh given its number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NUMBER, the number of the mesh.
!
!    Output, character ( len = * ) NAME, the name of the mesh.
!
  implicit none

  character ( len = * ) name
  integer ( kind = 4 ) number

  name = '?'

  if ( number == 1 ) then
    call mesh01_name ( name )
  else if ( number == 2 ) then
    call mesh02_name ( name )
  else if ( number == 3 ) then
    call mesh03_name ( name )
  else if ( number == 4 ) then
    call mesh04_name ( name )
  else if ( number == 5 ) then
    call mesh05_name ( name )
  else if ( number == 6 ) then
    call mesh06_name ( name )
  else if ( number == 7 ) then
    call mesh07_name ( name )
  else if ( number == 8 ) then
    call mesh08_name ( name )
  else if ( number == 9 ) then
    call mesh09_name ( name )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH00_NAME - Fatal error!'
    write ( *, '(a,i8)' ) '  Unknown mesh number = ', number
    stop
  end if

  return
end
subroutine mesh00_node_eps ( number, file_name )

!*****************************************************************************80
!
!! MESH00_NODE_EPS creates an image of the nodes of a mesh, given its number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NUMBER, the number of the mesh.
!
!    Input, character ( len = * ) FILE_NAME, the name for the file.
!
  implicit none

  character ( len = * ) file_name
  character ( len = 80 ) name
  logical, allocatable, dimension ( : ) :: node_mask
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: node_x
  real ( kind = 8 ), allocatable, dimension ( : ) :: node_y
  integer ( kind = 4 ) number

  call mesh00_name ( number, name )

  call mesh00_node_num ( number, node_num )

  allocate ( node_x(1:node_num) )
  allocate ( node_y(1:node_num) )
  allocate ( node_mask(1:node_num) )

  call mesh00_node_xy ( number, node_num, node_x, node_y )

  node_mask(1:node_num) = .true.

  call node_eps ( file_name, node_num, node_mask, node_x, node_y, 'The Nodes:' )

  deallocate ( node_mask)
  deallocate ( node_x )
  deallocate ( node_y )

  return
end
subroutine mesh00_node_num ( number, node_num )

!*****************************************************************************80
!
!! MESH00_NODE_NUM returns the number of nodes of a mesh given its number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NUMBER, the number of the mesh.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of elements in the mesh.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) number

  if ( number == 1 ) then
    call mesh01_node_num ( node_num )
  else if ( number == 2 ) then
    call mesh02_node_num ( node_num )
  else if ( number == 3 ) then
    call mesh03_node_num ( node_num )
  else if ( number == 4 ) then
    call mesh04_node_num ( node_num )
  else if ( number == 5 ) then
    call mesh05_node_num ( node_num )
  else if ( number == 6 ) then
    call mesh06_node_num ( node_num )
  else if ( number == 7 ) then
    call mesh07_node_num ( node_num )
  else if ( number == 8 ) then
    call mesh08_node_num ( node_num )
  else if ( number == 9 ) then
    call mesh09_node_num ( node_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH00_NODE_NUMBER - Fatal error!'
    write ( *, '(a,i8)' ) '  Unknown mesh number = ', number
    stop
  end if

  return
end
subroutine mesh00_node_xy ( number, node_num, x, y )

!*****************************************************************************80
!
!! MESH00_NODE_XY returns the node coordinates of a mesh given its number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NUMBER, the number of the mesh.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the coordinates
!    of the nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) number
  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) y(node_num)

  if ( number == 1 ) then
    call mesh01_node_xy ( node_num, x,  y )
  else if ( number == 2 ) then
    call mesh02_node_xy ( node_num, x,  y )
  else if ( number == 3 ) then
    call mesh03_node_xy ( node_num, x,  y )
  else if ( number == 4 ) then
    call mesh04_node_xy ( node_num, x,  y )
  else if ( number == 5 ) then
    call mesh05_node_xy ( node_num, x,  y )
  else if ( number == 6 ) then
    call mesh06_node_xy ( node_num, x,  y )
  else if ( number == 7 ) then
    call mesh07_node_xy ( node_num, x,  y )
  else if ( number == 8 ) then
    call mesh08_node_xy ( node_num, x,  y )
  else if ( number == 9 ) then
    call mesh09_node_xy ( node_num, x,  y )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH00_NODE_XY - Fatal error!'
    write ( *, '(a,i8)' ) '  Unknown mesh number = ', number
    stop
  end if

  return
end
subroutine mesh00_num ( mesh_num )

!*****************************************************************************80
!
!! MESH00_NUM returns the number of meshes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) MESH_NUM, the number of meshes.
!
  implicit none

  integer ( kind = 4 ) mesh_num

  mesh_num = 9

  return
end
subroutine mesh00_poly ( number, file_name )

!*****************************************************************************80
!
!! MESH00_POLY creates a POLY file of a mesh, given its number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NUMBER, the number of the mesh.
!
!    Input, character ( len = * ) FILE_NAME, the name for the file.
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: edge_nodes
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  character ( len = * ) file_name
  integer ( kind = 4 ) hole_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: hole_x
  real ( kind = 8 ), allocatable, dimension ( : ) :: hole_y
  character ( len = 80 ) name
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) number
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  real ( kind = 8 ), allocatable, dimension ( : ) :: y

  call mesh00_name ( number, name )

  call mesh00_node_num ( number, node_num )

  call mesh00_element_num ( number, element_num )

  allocate ( edge_nodes(2,3*element_num) )
  allocate ( x(1:node_num) )
  allocate ( y(1:node_num) )

  call mesh00_node_xy ( number, node_num, x, y )

  allocate ( element_node(1:3,1:element_num) )

  call mesh00_element_node ( number, element_num, element_node )

  call edge_list ( element_num, element_node, edge_num, edge_nodes )

  call mesh00_hole_num ( number, hole_num )

  allocate ( hole_x(hole_num) )
  allocate ( hole_y(hole_num) )

  call mesh00_hole_xy ( number, hole_num, hole_x, hole_y )

  call mesh_poly ( file_name, node_num, x, y, element_num, element_node, &
    edge_num, edge_nodes, hole_num, hole_x, hole_y )

  deallocate ( edge_nodes )
  deallocate ( element_node )
  deallocate ( hole_x )
  deallocate ( hole_y )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine mesh01_element_node ( element_num, element_node )

!*****************************************************************************80
!
!! MESH01_ELEMENT_NODE returns the element->node data of mesh #01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
  implicit none

  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) element_node(3,2)
  integer ( kind = 4 ), dimension(3,2) :: element_node_save &
    = reshape ( (/ &
    1, 2, 3, &
    4, 3, 2 /), (/ 3, 2 /) )

  element_node(1:3,1:2) = element_node_save(1:3,1:2)

  return
end
subroutine mesh01_element_num ( element_num )

!*****************************************************************************80
!
!! MESH01_ELEMENT_NUM returns the number of elements of mesh #01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
  implicit none

  integer ( kind = 4 ) element_num

  element_num = 2

  return
end
subroutine mesh01_hole_num ( hole_num )

!*****************************************************************************80
!
!! MESH01_HOLE_NUM returns the number of holes in mesh #01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 0

  return
end
subroutine mesh01_hole_xy ( hole_num, hole_x, hole_y )

!*****************************************************************************80
!
!! MESH01_HOLE_XY returns hole coordinates of mesh #01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) HOLE_X(HOLE_NUM), HOLE_Y(HOLE_NUM),
!    the XY coordinates a point in each hole.
!
  implicit none

  integer ( kind = 4 ) hole_num

  real ( kind = 8 ) hole_x(hole_num)
  real ( kind = 8 ) hole_y(hole_num)

  return
end
subroutine mesh01_name ( name )

!*****************************************************************************80
!
!! MESH01_NAME returns the name of mesh #01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the mesh.
!
  implicit none

  character ( len = * ) name

  name = 'The Square'

  return
end
subroutine mesh01_node_element_num ( node_element_num )

!*****************************************************************************80
!
!! MESH01_NODE_ELEMENT_NUM returns the number of nodes-element items of mesh #01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_ELEMENT_NUM, the number of nodes-element
!    data items.
!
  implicit none

  integer ( kind = 4 ) node_element_num

  node_element_num = 6

  return
end
subroutine mesh01_node_num ( node_num )

!*****************************************************************************80
!
!! MESH01_NODE_NUM returns the number of nodes of mesh #01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  node_num = 4

  return
end
subroutine mesh01_node_element ( node_num, node_element_num, &
  node_element_index, node_element )

!*****************************************************************************80
!
!! MESH01_NODE_ELEMENT returns the node-element data of mesh #01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) NODE_ELEMENT_NUM, the number of node-element
!    data items.
!
!    Output, integer ( kind = 4 ) NODE_ELEMENT_INDEX(NODE_NUM+1), for each 
!    node I, entry I points to the first entry in NODE_ELEMENT, and entry I+1
!    to the entry after the least entry in NODE_ELEMENT.
!
!    Output, integer ( kind = 4 ) NODE_ELEMENT(NODE_ELEMENT_NUM), for node I,
!    a list of the elements to which it belongs, in entries
!    NODE_ELEMENT_INDEX(I) through NODE_ELEMENT_INDEX(I+1)-1.
!
  implicit none

  integer ( kind = 4 ) node_element_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) node_element(node_element_num)
  integer ( kind = 4 ) node_element_index(node_num+1)

  node_element_index(1:node_num+1) = (/ 1, 2, 4, 6, 7 /)

  node_element = (/ 1, 1, 2, 1, 2, 2 /)

  return
end
subroutine mesh01_node_xy ( node_num, x, y )

!*****************************************************************************80
!
!! MESH01_NODE_XY returns the nodes of mesh #01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the XY
!    coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) y(node_num)

  x(1:node_num) = (/ 0.0D+00, 1.0D+00, 0.0D+00, 1.0D+00 /)
  y(1:node_num) = (/ 0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00 /)

  return
end
subroutine mesh02_element_node ( element_num, element_node )

!*****************************************************************************80
!
!! MESH02_ELEMENT_NODE returns the element->node data of mesh #02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
  implicit none

  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) element_node(3,30)
  integer ( kind = 4 ), dimension(3,30) :: element_node_save = reshape ( (/ &
     4,  1,  2, &
     6,  2,  3, &
     4, 10,  1, &
     5,  4,  2, &
     5,  2,  6, &
     6,  3, 16, &
     7,  4,  5, &
     9,  5,  6, &
    11, 10,  4, &
     7, 11,  4, &
     8,  7,  5, &
     8,  5,  9, &
     9,  6, 15, &
    15,  6, 16, &
     7,  8, 12, &
     9, 14,  8, &
    12, 11,  7, &
    13, 12,  8, &
    13,  8, 14, &
    14,  9, 15, &
    11, 19, 10, &
    12, 17, 11, &
    13, 18, 12, &
    19, 21, 10, &
    17, 19, 11, &
    17, 12, 18, &
    18, 20, 17, &
    17, 20, 19, &
    20, 22, 19, &
    19, 22, 21 /), (/ 3, 30 /) )

  element_node(1:3,1:30) = element_node_save(1:3,1:30)

  return
end
subroutine mesh02_element_num ( element_num )

!*****************************************************************************80
!
!! MESH02_ELEMENT_NUM returns the number of elements of mesh #02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
  implicit none

  integer ( kind = 4 ) element_num

  element_num = 30

  return
end
subroutine mesh02_hole_num ( hole_num )

!*****************************************************************************80
!
!! MESH02_HOLE_NUM returns the number of holes in mesh #02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 0

  return
end
subroutine mesh02_hole_xy ( hole_num, hole_x, hole_y )

!*****************************************************************************80
!
!! MESH02_HOLE_XY returns hole coordinates of mesh #02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) HOLE_X(HOLE_NUM), HOLE_Y(HOLE_NUM),
!    the XY coordinates a point in each hole.
!
  implicit none

  integer ( kind = 4 ) hole_num

  real ( kind = 8 ) hole_x(hole_num)
  real ( kind = 8 ) hole_y(hole_num)

  return
end
subroutine mesh02_name ( name )

!*****************************************************************************80
!
!! MESH02_NAME returns the name of mesh #02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the mesh.
!
  implicit none

  character ( len = * ) name

  name = 'The Graded L'

  return
end
subroutine mesh02_node_element_num ( node_element_num )

!*****************************************************************************80
!
!! MESH02_NODE_ELEMENT_NUM returns the number of nodes-element items of mesh #02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_ELEMENT_NUM, the number of nodes-element
!    data items.
!
  implicit none

  integer ( kind = 4 ) node_element_num

  node_element_num = 90

  return
end
subroutine mesh02_node_element ( node_num, node_element_num, &
  node_element_index, node_element )

!*****************************************************************************80
!
!! MESH02_NODE_ELEMENT returns the node-element data of mesh #02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) NODE_ELEMENT_NUM, the number of node-element
!    data items.
!
!    Output, integer ( kind = 4 ) NODE_ELEMENT_INDEX(NODE_NUM+1), for each 
!    node I, entry I points to the first entry in NODE_ELEMENT, and entry I+1
!    to the entry after the least entry in NODE_ELEMENT.
!
!    Output, integer ( kind = 4 ) NODE_ELEMENT(NODE_ELEMENT_NUM), for node I,
!    a list of the elements to which it belongs, in entries
!    NODE_ELEMENT_INDEX(I) through NODE_ELEMENT_INDEX(I+1)-1.
!
  implicit none

  integer ( kind = 4 ) node_element_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) node_element(node_element_num)
  integer ( kind = 4 ) node_element_index(node_num+1)

  node_element_index(1:node_num+1) = (/ &
    1,  3,  7,  9, 15, 21, 27, 32, 38, 43, &
   47, 53, 59, 62, 65, 68, 70, 75, 78, 84, &
   87, 89, 91 /)

  node_element = (/ &
     1,  3,  1,  2,  4,  5,  2,  6,  1,  3, &
     4,  7,  9, 10,  4,  5,  7,  8, 11, 12, &
     2,  5,  6,  8, 13, 14,  7, 10, 11, 15, &
    17, 11, 12, 15, 16, 18, 19,  8, 12, 13, &
    16, 20,  3,  9, 21, 24,  9, 10, 17, 21, &
    22, 25, 15, 17, 18, 22, 23, 26, 18, 19, &
    23, 16, 19, 20, 13, 14, 20,  6, 14, 22, &
    25, 26, 27, 28, 23, 26, 27, 21, 24, 25, &
    28, 29, 30, 27, 28, 29, 24, 30, 29, 30 /)

  return
end
subroutine mesh02_node_num ( node_num )

!*****************************************************************************80
!
!! MESH02_NODE_NUM returns the number of nodes of mesh #02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  node_num = 22

  return
end
subroutine mesh02_node_xy ( node_num, x, y )

!*****************************************************************************80
!
!! MESH02_NODE_XY returns the nodes of mesh #02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the XY
!    coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) y(node_num)

  x(1:node_num) = (/ &
    0.0D+00, 4.0D+00, 8.0D+00, 2.0D+00, 4.0D+00, &
    6.0D+00, 3.0D+00, 4.0D+00, 5.0D+00, 0.0D+00, &
    2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00, 6.0D+00, &
    8.0D+00, 3.0D+00, 4.0D+00, 2.0D+00, 4.0D+00, &
    0.0D+00, 4.0D+00 /) / 8.0D+00

  y(1:node_num) = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 2.0D+00, 2.0D+00, &
    2.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, 4.0D+00, &
    4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, &
    4.0D+00, 5.0D+00, 5.0D+00, 6.0D+00, 6.0D+00, &
    8.0D+00, 8.0D+00 /) / 8.0D+00

  return
end
subroutine mesh03_element_node ( element_num, element_node )

!*****************************************************************************80
!
!! MESH03_ELEMENT_NODE returns the element->node data of mesh #03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
  implicit none

  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) element_node(3,26)
  integer ( kind = 4 ), dimension ( 3, 26 ) :: element_node_save = reshape ( (/ &
     3,  1,  4, &
     2,  4,  1, &
     5,  3,  6, &
     4,  6,  3, &
    10,  5, 11, &
     6, 11,  5, &
    15,  7, 16, &
     8, 16,  7, &
    16,  8, 17, &
     9, 17,  8, &
    17,  9, 18, &
    10, 18,  9, &
    18, 10, 19, &
    11, 19, 10, &
    19, 11, 20, &
    12, 20, 11, &
    20, 12, 21, &
    13, 21, 12, &
    21, 13, 22, &
    14, 22, 13, &
    23, 18, 24, &
    19, 24, 18, &
    25, 23, 26, &
    24, 26, 23, &
    27, 25, 28, &
    26, 28, 25 /), (/ 3, 26 /) )

  element_node(1:3,1:26) = element_node_save(1:3,1:26)

  return
end
subroutine mesh03_element_num ( element_num )

!*****************************************************************************80
!
!! MESH03_ELEMENT_NUM returns the number of elements of mesh #03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
  implicit none

  integer ( kind = 4 ) element_num

  element_num = 26

  return
end
subroutine mesh03_hole_num ( hole_num )

!*****************************************************************************80
!
!! MESH03_HOLE_NUM returns the number of holes in mesh #03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 0

  return
end
subroutine mesh03_hole_xy ( hole_num, hole_x, hole_y )

!*****************************************************************************80
!
!! MESH03_HOLE_XY returns hole coordinates of mesh #03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) HOLE_X(HOLE_NUM), HOLE_Y(HOLE_NUM),
!    the XY coordinates a point in each hole.
!
  implicit none

  integer ( kind = 4 ) hole_num

  real ( kind = 8 ) hole_x(hole_num)
  real ( kind = 8 ) hole_y(hole_num)

  return
end
subroutine mesh03_name ( name )

!*****************************************************************************80
!
!! MESH03_NAME returns the name of mesh #03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the mesh.
!
  implicit none

  character ( len = * ) name

  name = 'The Plus Shaped Domain'

  return
end
subroutine mesh03_node_num ( node_num )

!*****************************************************************************80
!
!! MESH03_NODE_NUM returns the number of nodes of mesh #03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  node_num = 28

  return
end
subroutine mesh03_node_xy ( node_num, x, y )

!*****************************************************************************80
!
!! MESH03_NODE_XY returns the nodes of mesh #03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the XY
!    coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) y(node_num)

  x(1:node_num) = (/ &
    3.0D+00, 4.0D+00, 3.0D+00, 4.0D+00, 3.0D+00, &
    4.0D+00, 0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00, &
    4.0D+00, 5.0D+00, 6.0D+00, 7.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00, &
    6.0D+00, 7.0D+00, 3.0D+00, 4.0D+00, 3.0D+00, &
    4.0D+00, 3.0D+00, 4.0D+00 /) / 7.0D+00

  y(1:node_num) = (/ &
    0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00, 2.0D+00, &
    2.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, &
    3.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, 4.0D+00, &
    4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, &
    4.0D+00, 4.0D+00, 5.0D+00, 5.0D+00, 6.0D+00, &
    6.0D+00, 7.0D+00, 7.0D+00 /) / 7.0D+00

  return
end
subroutine mesh04_element_node ( element_num, element_node )

!*****************************************************************************80
!
!! MESH04_ELEMENT_NODE returns the element->node data of mesh #04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
  implicit none

  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) element_node(3,element_num)

  element_node(1:3, 1) = (/  1,  2,  5 /)
  element_node(1:3, 2) = (/  6,  5,  2 /)
  element_node(1:3, 3) = (/  3,  4,  7 /)
  element_node(1:3, 4) = (/  8,  7,  4 /)
  element_node(1:3, 5) = (/  5,  6,  9 /)
  element_node(1:3, 6) = (/ 10,  9,  6 /)
  element_node(1:3, 7) = (/  7,  8, 11 /)
  element_node(1:3, 8) = (/ 12, 11,  8 /)
  element_node(1:3, 9) = (/  9, 10, 13 /)
  element_node(1:3,10) = (/ 14, 13, 10 /)
  element_node(1:3,11) = (/ 11, 12, 19 /)
  element_node(1:3,12) = (/ 20, 19, 12 /)
  element_node(1:3,13) = (/ 13, 14, 21 /)
  element_node(1:3,14) = (/ 22, 21, 14 /)
  element_node(1:3,15) = (/ 14, 15, 22 /)
  element_node(1:3,16) = (/ 23, 22, 15 /)
  element_node(1:3,17) = (/ 15, 16, 23 /)
  element_node(1:3,18) = (/ 24, 23, 16 /)
  element_node(1:3,19) = (/ 16, 17, 24 /)
  element_node(1:3,20) = (/ 25, 24, 17 /)
  element_node(1:3,21) = (/ 17, 18, 25 /)
  element_node(1:3,22) = (/ 26, 25, 18 /)
  element_node(1:3,23) = (/ 18, 19, 26 /)
  element_node(1:3,24) = (/ 27, 26, 19 /)
  element_node(1:3,25) = (/ 19, 20, 27 /)
  element_node(1:3,26) = (/ 28, 27, 20 /)
  element_node(1:3,27) = (/ 21, 22, 29 /)
  element_node(1:3,28) = (/ 30, 29, 22 /)
  element_node(1:3,29) = (/ 27, 28, 31 /)
  element_node(1:3,30) = (/ 32, 31, 28 /)
  element_node(1:3,31) = (/ 29, 30, 33 /)
  element_node(1:3,32) = (/ 34, 33, 30 /)
  element_node(1:3,33) = (/ 31, 32, 35 /)
  element_node(1:3,34) = (/ 36, 35, 32 /)
  element_node(1:3,35) = (/ 33, 34, 37 /)
  element_node(1:3,36) = (/ 38, 37, 34 /)
  element_node(1:3,37) = (/ 35, 36, 39 /)
  element_node(1:3,38) = (/ 40, 39, 36 /)

  return
end
subroutine mesh04_element_num ( element_num )

!*****************************************************************************80
!
!! MESH04_ELEMENT_NUM returns the number of elements of mesh #04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
  implicit none

  integer ( kind = 4 ) element_num

  element_num = 38

  return
end
subroutine mesh04_hole_num ( hole_num )

!*****************************************************************************80
!
!! MESH04_HOLE_NUM returns the number of holes in mesh #04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 0

  return
end
subroutine mesh04_hole_xy ( hole_num, hole_x, hole_y )

!*****************************************************************************80
!
!! MESH04_HOLE_XY returns hole coordinates of mesh #04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) HOLE_X(HOLE_NUM), HOLE_Y(HOLE_NUM),
!    the XY coordinates a point in each hole.
!
  implicit none

  integer ( kind = 4 ) hole_num

  real ( kind = 8 ) hole_x(hole_num)
  real ( kind = 8 ) hole_y(hole_num)

  return
end
subroutine mesh04_name ( name )

!*****************************************************************************80
!
!! MESH04_NAME returns the name of mesh #04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the mesh.
!
  implicit none

  character ( len = * ) name

  name = 'The H Shaped Domain'

  return
end
subroutine mesh04_node_num ( node_num )

!*****************************************************************************80
!
!! MESH04_NODE_NUM returns the number of nodes of mesh #04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  node_num = 40

  return
end
subroutine mesh04_node_xy ( node_num, x, y )

!*****************************************************************************80
!
!! MESH04_NODE_XY returns the nodes of mesh #04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the XY
!    coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) y(node_num)

  x(1:node_num) = (/ &
    0.0D+00, 1.0D+00, 6.0D+00, 7.0D+00, 0.0D+00, &
    1.0D+00, 6.0D+00, 7.0D+00, 0.0D+00, 1.0D+00, &
    6.0D+00, 7.0D+00, 0.0D+00, 1.0D+00, 2.0D+00, &
    3.0D+00, 4.0D+00, 5.0D+00, 6.0D+00, 7.0D+00, &
    0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, &
    5.0D+00, 6.0D+00, 7.0D+00, 0.0D+00, 1.0D+00, &
    6.0D+00, 7.0D+00, 0.0D+00, 1.0D+00, 6.0D+00, &
    7.0D+00, 0.0D+00, 1.0D+00, 6.0D+00, 7.0D+00 /) / 7.0D+00

  y(1:node_num) = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 2.0D+00, 2.0D+00, &
    2.0D+00, 2.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, &
    3.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, &
    4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, &
    4.0D+00, 4.0D+00, 4.0D+00, 5.0D+00, 5.0D+00, &
    5.0D+00, 5.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
    6.0D+00, 7.0D+00, 7.0D+00, 7.0D+00, 7.0D+00 /) / 7.0D+00

  return
end
subroutine mesh05_element_node ( element_num, element_node )

!*****************************************************************************80
!
!! MESH05_ELEMENT_NODE returns the element->node data of mesh #05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
  implicit none

  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) element_node(3,element_num)

  element_node(1:3, 1) = (/  4,  1,  2 /)
  element_node(1:3, 2) = (/  2,  5,  4 /)
  element_node(1:3, 3) = (/  5,  2,  3 /)
  element_node(1:3, 4) = (/  5,  3,  7 /)
  element_node(1:3, 5) = (/  7,  9,  5 /)
  element_node(1:3, 6) = (/  9,  7, 12 /)
  element_node(1:3, 7) = (/  9, 12, 11 /)
  element_node(1:3, 8) = (/ 11,  8,  9 /)
  element_node(1:3, 9) = (/  8, 11, 10 /)
  element_node(1:3,10) = (/  8, 10,  6 /)
  element_node(1:3,11) = (/  6,  4,  8 /)
  element_node(1:3,12) = (/  4,  6,  1 /)

  return
end
subroutine mesh05_element_num ( element_num )

!*****************************************************************************80
!
!! MESH05_ELEMENT_NUM returns the number of elements of mesh #05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
  implicit none

  integer ( kind = 4 ) element_num

  element_num = 12

  return
end
subroutine mesh05_hole_num ( hole_num )

!*****************************************************************************80
!
!! MESH05_HOLE_NUM returns the number of holes in mesh #05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 1

  return
end
subroutine mesh05_hole_xy ( hole_num, hole_x, hole_y )

!*****************************************************************************80
!
!! MESH05_HOLE_XY returns hole coordinates of mesh #05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) HOLE_X(HOLE_NUM), HOLE_Y(HOLE_NUM),
!    the XY coordinates a point in each hole.
!
  implicit none

  integer ( kind = 4 ) hole_num

  real ( kind = 8 ) hole_x(hole_num)
  real ( kind = 8 ) hole_y(hole_num)

  hole_x(1:hole_num) = (/ 0.5D+00 /)
  hole_y(1:hole_num) = (/ 0.5D+00 /)

  return
end
subroutine mesh05_name ( name )

!*****************************************************************************80
!
!! MESH05_NAME returns the name of mesh #05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the mesh.
!
  implicit none

  character ( len = * ) name

  name = 'The Hollow Square (Small Hole)'

  return
end
subroutine mesh05_node_num ( node_num )

!*****************************************************************************80
!
!! MESH05_NODE_NUM returns the number of nodes of mesh #05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  node_num = 12

  return
end
subroutine mesh05_node_xy ( node_num, x, y )

!*****************************************************************************80
!
!! MESH05_NODE_XY returns the nodes of mesh #05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the XY
!    coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) y(node_num)

  x(1:node_num) = (/ &
    0.0D+00, 2.0D+00, 4.0D+00, 1.0D+00, 3.0D+00, &
    0.0D+00, 4.0D+00, 1.0D+00, 3.0D+00, 0.0D+00, &
    2.0D+00, 4.0D+00 /) / 4.0D+00
  y(1:node_num) = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00, &
    2.0D+00, 2.0D+00, 3.0D+00, 3.0D+00, 4.0D+00, &
    4.0D+00, 4.0D+00 /) / 4.0D+00

  return
end
subroutine mesh06_element_node ( element_num, element_node )

!*****************************************************************************80
!
!! MESH06_ELEMENT_NODE returns the element->node data of mesh #06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
  implicit none

  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) element_node(3,element_num)

  element_node(1:3, 1) = (/  1,  2,  7 /)
  element_node(1:3, 2) = (/  8,  7,  2 /)
  element_node(1:3, 3) = (/  2,  3,  8 /)
  element_node(1:3, 4) = (/  9,  8,  3 /)
  element_node(1:3, 5) = (/  3,  4,  9 /)
  element_node(1:3, 6) = (/ 10,  9,  4 /)
  element_node(1:3, 7) = (/  4,  5, 10 /)
  element_node(1:3, 8) = (/ 11, 10,  5 /)
  element_node(1:3, 9) = (/  5,  6, 11 /)
  element_node(1:3,10) = (/ 12, 11,  6 /)
  element_node(1:3,11) = (/ 11, 12, 15 /)
  element_node(1:3,12) = (/ 16, 15, 12 /)
  element_node(1:3,13) = (/ 15, 16, 19 /)
  element_node(1:3,14) = (/ 20, 19, 16 /)
  element_node(1:3,15) = (/ 19, 20, 25 /)
  element_node(1:3,16) = (/ 26, 25, 20 /)
  element_node(1:3,17) = (/ 25, 26, 31 /)
  element_node(1:3,18) = (/ 32, 31, 26 /)
  element_node(1:3,19) = (/ 31, 30, 25 /)
  element_node(1:3,20) = (/ 24, 25, 30 /)
  element_node(1:3,21) = (/ 30, 29, 24 /)
  element_node(1:3,22) = (/ 23, 24, 29 /)
  element_node(1:3,23) = (/ 29, 28, 23 /)
  element_node(1:3,24) = (/ 22, 23, 28 /)
  element_node(1:3,25) = (/ 28, 27, 22 /)
  element_node(1:3,26) = (/ 21, 22, 27 /)
  element_node(1:3,27) = (/ 22, 21, 18 /)
  element_node(1:3,28) = (/ 17, 18, 21 /)
  element_node(1:3,29) = (/ 18, 17, 14 /)
  element_node(1:3,30) = (/ 13, 14, 17 /)
  element_node(1:3,31) = (/ 14, 13,  8 /)
  element_node(1:3,32) = (/  7,  8, 13 /)

  return
end
subroutine mesh06_element_num ( element_num )

!*****************************************************************************80
!
!! MESH06_ELEMENT_NUM returns the number of elements of mesh #06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
  implicit none

  integer ( kind = 4 ) element_num

  element_num = 32

  return
end
subroutine mesh06_hole_num ( hole_num )

!*****************************************************************************80
!
!! MESH06_HOLE_NUM returns the number of holes in mesh #06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 1

  return
end
subroutine mesh06_hole_xy ( hole_num, hole_x, hole_y )

!*****************************************************************************80
!
!! MESH06_HOLE_XY returns hole coordinates of mesh #06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) HOLE_X(HOLE_NUM), HOLE_Y(HOLE_NUM),
!    the XY coordinates a point in each hole.
!
  implicit none

  integer ( kind = 4 ) hole_num

  real ( kind = 8 ) hole_x(hole_num)
  real ( kind = 8 ) hole_y(hole_num)

  hole_x(1:hole_num) = (/ 0.5D+00 /)
  hole_y(1:hole_num) = (/ 0.5D+00 /)

  return
end
subroutine mesh06_name ( name )

!*****************************************************************************80
!
!! MESH06_NAME returns the name of mesh #06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the mesh.
!
  implicit none

  character ( len = * ) name

  name = 'The Hollow Square (Large Hole)'

  return
end
subroutine mesh06_node_num ( node_num )

!*****************************************************************************80
!
!! MESH06_NODE_NUM returns the number of nodes of mesh #06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  node_num = 32

  return
end
subroutine mesh06_node_xy ( node_num, x, y )

!*****************************************************************************80
!
!! MESH06_NODE_XY returns the nodes of mesh #06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the XY
!    coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) y(node_num)

  x(1:node_num) = (/ &
    0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, &
    5.0D+00, 0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00, &
    4.0D+00, 5.0D+00, 0.0D+00, 1.0D+00, 4.0D+00, &
    5.0D+00, 0.0D+00, 1.0D+00, 4.0D+00, 5.0D+00, &
    0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, &
    5.0D+00, 0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00, &
    4.0D+00, 5.0D+00 /) / 5.0D+00
  y(1:node_num) = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 2.0D+00, 2.0D+00, 2.0D+00, &
    2.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, &
    4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, &
    4.0D+00, 5.0D+00, 5.0D+00, 5.0D+00, 5.0D+00, &
    5.0D+00, 5.0D+00 /) / 5.0D+00

  return
end
subroutine mesh07_element_node ( element_num, element_node )

!*****************************************************************************80
!
!! MESH07_ELEMENT_NODE returns the element->node data of mesh #07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
  implicit none

  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) element_node(3,element_num)

  element_node(1:3, 1) = (/  8,  7,  1 /)
  element_node(1:3, 2) = (/  1,  2,  8 /)
  element_node(1:3, 3) = (/  2,  9,  8 /)
  element_node(1:3, 4) = (/  9,  2,  3 /)
  element_node(1:3, 5) = (/ 10,  9,  3 /)
  element_node(1:3, 6) = (/ 10,  3, 11 /)
  element_node(1:3, 7) = (/ 11,  3,  4 /)
  element_node(1:3, 8) = (/ 12, 11,  4 /)
  element_node(1:3, 9) = (/ 12,  4, 13 /)
  element_node(1:3,10) = (/ 13,  4,  5 /)
  element_node(1:3,11) = (/ 13,  5, 21 /)
  element_node(1:3,12) = (/ 14, 21,  5 /)
  element_node(1:3,13) = (/ 16, 15,  6 /)
  element_node(1:3,14) = (/  6,  7, 16 /)
  element_node(1:3,15) = (/ 17, 16,  7 /)
  element_node(1:3,16) = (/ 17,  7, 18 /)
  element_node(1:3,17) = (/  8, 18,  7 /)
  element_node(1:3,18) = (/ 10, 11, 19 /)
  element_node(1:3,19) = (/ 20, 19, 11 /)
  element_node(1:3,20) = (/ 11, 12, 20 /)
  element_node(1:3,21) = (/ 21, 26, 13 /)
  element_node(1:3,22) = (/ 21, 14, 27 /)
  element_node(1:3,23) = (/ 22, 28, 15 /)
  element_node(1:3,24) = (/ 16, 22, 15 /)
  element_node(1:3,25) = (/ 17, 18, 23 /)
  element_node(1:3,26) = (/ 23, 18, 30 /)
  element_node(1:3,27) = (/ 24, 30, 18 /)
  element_node(1:3,28) = (/ 19, 20, 25 /)
  element_node(1:3,29) = (/ 25, 20, 31 /)
  element_node(1:3,30) = (/ 26, 21, 33 /)
  element_node(1:3,31) = (/ 27, 33, 21 /)
  element_node(1:3,32) = (/ 29, 28, 22 /)
  element_node(1:3,33) = (/ 30, 35, 23 /)
  element_node(1:3,34) = (/ 30, 24, 36 /)
  element_node(1:3,35) = (/ 31, 37, 25 /)
  element_node(1:3,36) = (/ 33, 27, 41 /)
  element_node(1:3,37) = (/ 29, 34, 28 /)
  element_node(1:3,38) = (/ 34, 29, 35 /)
  element_node(1:3,39) = (/ 35, 30, 43 /)
  element_node(1:3,40) = (/ 36, 43, 30 /)
  element_node(1:3,41) = (/ 37, 31, 46 /)
  element_node(1:3,42) = (/ 38, 46, 31 /)
  element_node(1:3,43) = (/ 38, 31, 32 /)
  element_node(1:3,44) = (/ 39, 38, 32 /)
  element_node(1:3,45) = (/ 40, 33, 41 /)
  element_node(1:3,46) = (/ 35, 42, 34 /)
  element_node(1:3,47) = (/ 42, 35, 43 /)
  element_node(1:3,48) = (/ 44, 43, 36 /)
  element_node(1:3,49) = (/ 44, 36, 45 /)
  element_node(1:3,50) = (/ 45, 36, 37 /)
  element_node(1:3,51) = (/ 37, 46, 45 /)
  element_node(1:3,52) = (/ 38, 47, 46 /)
  element_node(1:3,53) = (/ 39, 47, 38 /)
  element_node(1:3,54) = (/ 47, 39, 48 /)
  element_node(1:3,55) = (/ 40, 48, 39 /)
  element_node(1:3,56) = (/ 40, 41, 48 /)

  return
end
subroutine mesh07_element_num ( element_num )

!*****************************************************************************80
!
!! MESH07_ELEMENT_NUM returns the number of elements of mesh #07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
  implicit none

  integer ( kind = 4 ) element_num

  element_num = 56

  return
end
subroutine mesh07_hole_num ( hole_num )

!*****************************************************************************80
!
!! MESH07_HOLE_NUM returns the number of holes in mesh #07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 3

  return
end
subroutine mesh07_hole_xy ( hole_num, hole_x, hole_y )

!*****************************************************************************80
!
!! MESH07_HOLE_XY returns hole coordinates of mesh #07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) HOLE_X(HOLE_NUM), HOLE_Y(HOLE_NUM),
!    the XY coordinates a point in each hole.
!
  implicit none

  integer ( kind = 4 ) hole_num

  real ( kind = 8 ) hole_x(hole_num)
  real ( kind = 8 ) hole_y(hole_num)

  hole_x(1:hole_num) = (/ 0.135D+00, 0.454D+00, 0.727D+00 /)
  hole_y(1:hole_num) = (/ 0.583D+00, 0.500D+00, 0.583D+00 /)

  return
end
subroutine mesh07_name ( name )

!*****************************************************************************80
!
!! MESH07_NAME returns the name of mesh #07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the mesh.
!
  implicit none

  character ( len = * ) name

  name = 'The 3-Hole Problem'

  return
end
subroutine mesh07_node_num ( node_num )

!*****************************************************************************80
!
!! MESH07_NODE_NUM returns the number of nodes of mesh #07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  node_num = 48

  return
end
subroutine mesh07_node_xy ( node_num, x, y )

!*****************************************************************************80
!
!! MESH07_NODE_XY returns the nodes of mesh #07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the XY
!    coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) y(node_num)

  x(1:node_num) = (/ &
    3.0D+00, 4.0D+00, 6.0D+00, 8.0D+00,10.0D+00, &
    1.0D+00, 2.0D+00, 3.0D+00, 5.0D+00, 6.0D+00, &
    7.0D+00, 8.0D+00, 9.0D+00,11.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00, 3.0D+00, 6.0D+00, 7.0D+00, &
   10.0D+00, 1.0D+00, 2.0D+00, 4.0D+00, 6.0D+00, &
    9.0D+00,11.0D+00, 0.0D+00, 1.0D+00, 3.0D+00, &
    7.0D+00, 9.0D+00,10.0D+00, 1.0D+00, 2.0D+00, &
    4.0D+00, 6.0D+00, 8.0D+00, 9.0D+00,10.0D+00, &
   11.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00, &
    7.0D+00, 9.0D+00,10.0D+00 /) / 11.0D+00

  y(1:node_num) = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 2.0D+00, &
    2.0D+00, 2.0D+00, 2.0D+00, 2.0D+00, 2.0D+00, &
    2.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, &
    3.0D+00, 3.0D+00, 4.0D+00, 4.0D+00, 4.0D+00, &
    4.0D+00, 4.0D+00, 4.0D+00, 5.0D+00, 5.0D+00, &
    5.0D+00, 5.0D+00, 5.0D+00, 5.0D+00, 5.0D+00, &
    5.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
    6.0D+00, 6.0D+00, 6.0D+00 /) / 6.0D+00

  return
end
subroutine mesh08_element_node ( element_num, element_node )

!*****************************************************************************80
!
!! MESH08_ELEMENT_NODE returns the element->node data of mesh #08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
  implicit none

  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) element_node(3,element_num)

  element_node(1:3, 1) = (/  1,  7,  6 /)
  element_node(1:3, 2) = (/  7,  1,  2 /)
  element_node(1:3, 3) = (/  2,  8,  7 /)
  element_node(1:3, 4) = (/  8,  2,  3 /)
  element_node(1:3, 5) = (/  3,  9,  8 /)
  element_node(1:3, 6) = (/  9,  3,  4 /)
  element_node(1:3, 7) = (/  4, 10,  9 /)
  element_node(1:3, 8) = (/ 10,  4,  5 /)
  element_node(1:3, 9) = (/  5, 11, 10 /)
  element_node(1:3,10) = (/ 12,  6, 13 /)
  element_node(1:3,11) = (/  7, 13,  6 /)
  element_node(1:3,12) = (/  8,  9, 14 /)
  element_node(1:3,13) = (/ 15, 14,  9 /)
  element_node(1:3,14) = (/ 16, 10, 17 /)
  element_node(1:3,15) = (/ 11, 17, 10 /)
  element_node(1:3,16) = (/ 18, 12, 13 /)
  element_node(1:3,17) = (/ 13, 19, 18 /)
  element_node(1:3,18) = (/ 19, 13, 14 /)
  element_node(1:3,19) = (/ 14, 20, 19 /)
  element_node(1:3,20) = (/ 20, 14, 15 /)
  element_node(1:3,21) = (/ 15, 21, 20 /)
  element_node(1:3,22) = (/ 21, 15, 16 /)
  element_node(1:3,23) = (/ 16, 22, 21 /)
  element_node(1:3,24) = (/ 22, 16, 17 /)
  element_node(1:3,25) = (/ 17, 23, 22 /)
  element_node(1:3,26) = (/ 24, 18, 25 /)
  element_node(1:3,27) = (/ 19, 25, 18 /)
  element_node(1:3,28) = (/ 26, 20, 27 /)
  element_node(1:3,29) = (/ 21, 27, 20 /)
  element_node(1:3,30) = (/ 22, 23, 28 /)
  element_node(1:3,31) = (/ 29, 28, 23 /)
  element_node(1:3,32) = (/ 24, 31, 30 /)
  element_node(1:3,33) = (/ 31, 24, 25 /)
  element_node(1:3,34) = (/ 25, 32, 31 /)
  element_node(1:3,35) = (/ 32, 25, 26 /)
  element_node(1:3,36) = (/ 26, 33, 32 /)
  element_node(1:3,37) = (/ 33, 26, 27 /)
  element_node(1:3,38) = (/ 27, 34, 33 /)
  element_node(1:3,39) = (/ 34, 27, 28 /)
  element_node(1:3,40) = (/ 28, 35, 34 /)
  element_node(1:3,41) = (/ 35, 28, 29 /)
  element_node(1:3,42) = (/ 36, 30, 37 /)
  element_node(1:3,43) = (/ 31, 37, 30 /)
  element_node(1:3,44) = (/ 32, 33, 38 /)
  element_node(1:3,45) = (/ 39, 38, 33 /)
  element_node(1:3,46) = (/ 40, 34, 41 /)
  element_node(1:3,47) = (/ 35, 41, 34 /)
  element_node(1:3,48) = (/ 42, 36, 37 /)
  element_node(1:3,49) = (/ 37, 43, 42 /)
  element_node(1:3,50) = (/ 43, 37, 38 /)
  element_node(1:3,51) = (/ 38, 44, 43 /)
  element_node(1:3,52) = (/ 44, 38, 39 /)
  element_node(1:3,53) = (/ 39, 45, 44 /)
  element_node(1:3,54) = (/ 45, 39, 40 /)
  element_node(1:3,55) = (/ 40, 46, 45 /)
  element_node(1:3,56) = (/ 46, 40, 41 /)

  return
end
subroutine mesh08_element_num ( element_num )

!*****************************************************************************80
!
!! MESH08_ELEMENT_NUM returns the number of elements of mesh #08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
  implicit none

  integer ( kind = 4 ) element_num

  element_num = 56

  return
end
subroutine mesh08_hole_num ( hole_num )

!*****************************************************************************80
!
!! MESH08_HOLE_NUM returns the number of holes in mesh #08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 6

  return
end
subroutine mesh08_hole_xy ( hole_num, hole_x, hole_y )

!*****************************************************************************80
!
!! MESH08_HOLE_XY returns hole coordinates of mesh #08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) HOLE_X(HOLE_NUM), HOLE_Y(HOLE_NUM),
!    the XY coordinates a point in each hole.
!
  implicit none

  integer ( kind = 4 ) hole_num

  real ( kind = 8 ) hole_x(hole_num)
  real ( kind = 8 ) hole_y(hole_num)

  hole_x(1:hole_num) = (/ &
    0.272D+00, 0.636D+00, 0.363D+00, 0.727D+00, 0.272D+00, 0.636D+00 /)

  hole_y(1:hole_num) = (/ &
    0.200D+00, 0.200D+00, 0.500D+00, 0.500D+00, 0.800D+00, 0.800D+00 /)

  return
end
subroutine mesh08_name ( name )

!*****************************************************************************80
!
!! MESH08_NAME returns the name of mesh #08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the mesh.
!
  implicit none

  character ( len = * ) name

  name = 'The 6-Hole Problem'

  return
end
subroutine mesh08_node_num ( node_num )

!*****************************************************************************80
!
!! MESH08_NODE_NUM returns the number of nodes of mesh #08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  node_num = 46

  return
end
subroutine mesh08_node_xy ( node_num, x, y )

!*****************************************************************************80
!
!! MESH08_NODE_XY returns the nodes of mesh #08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the XY
!    coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) y(node_num)

  x(1:node_num) = (/ &
     1.0D+00,  3.0D+00,  5.0D+00,  7.0D+00,  9.0D+00, &
     0.0D+00,  2.0D+00,  4.0D+00,  6.0D+00,  8.0D+00, &
    10.0D+00,  0.0D+00,  2.0D+00,  4.0D+00,  6.0D+00, &
     8.0D+00, 10.0D+00,  1.0D+00,  3.0D+00,  5.0D+00, &
     7.0D+00,  9.0D+00, 11.0D+00,  1.0D+00,  3.0D+00, &
     5.0D+00,  7.0D+00,  9.0D+00, 11.0D+00,  0.0D+00, &
     2.0D+00,  4.0D+00,  6.0D+00,  8.0D+00, 10.0D+00, &
     0.0D+00,  2.0D+00,  4.0D+00,  6.0D+00,  8.0D+00, &
    10.0D+00,  1.0D+00,  3.0D+00,  5.0D+00,  7.0D+00, &
     9.0D+00 /) / 11.0D+00

  y(1:node_num) = (/ &
     0.0D+00,  0.0D+00,  0.0D+00,  0.0D+00,  0.0D+00, &
     1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00, &
     1.0D+00,  3.0D+00,  3.0D+00,  3.0D+00,  3.0D+00, &
     3.0D+00,  3.0D+00,  4.0D+00,  4.0D+00,  4.0D+00, &
     4.0D+00,  4.0D+00,  4.0D+00,  6.0D+00,  6.0D+00, &
     6.0D+00,  6.0D+00,  6.0D+00,  6.0D+00,  7.0D+00, &
     7.0D+00,  7.0D+00,  7.0D+00,  7.0D+00,  7.0D+00, &
     9.0D+00,  9.0D+00,  9.0D+00,  9.0D+00,  9.0D+00, &
     9.0D+00, 10.0D+00, 10.0D+00, 10.0D+00, 10.0D+00,&
    10.0D+00 /) / 10.0D+00

  return
end
subroutine mesh09_element_node ( element_num, element_node )

!*****************************************************************************80
!
!! MESH09_ELEMENT_NODE returns the element->node data of mesh #09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    element->node data.
!
  implicit none

  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) element_node(3,element_num)

  element_node(1:3, 1) = (/  4,  2,  1 /)
  element_node(1:3, 2) = (/  5,  1,  3 /)
  element_node(1:3, 3) = (/  4,  6,  2 /)
  element_node(1:3, 4) = (/  5,  3,  8 /)
  element_node(1:3, 5) = (/  4,  7,  6 /)
  element_node(1:3, 6) = (/  7,  4,  5 /)
  element_node(1:3, 7) = (/  5,  8,  7 /)

  return
end
subroutine mesh09_element_num ( element_num )

!*****************************************************************************80
!
!! MESH09_ELEMENT_NUM returns the number of elements of mesh #09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
  implicit none

  integer ( kind = 4 ) element_num

  element_num = 7

  return
end
subroutine mesh09_hole_num ( hole_num )

!*****************************************************************************80
!
!! MESH09_HOLE_NUM returns the number of holes in mesh #09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) HOLE_NUM, the number of holes.
!
  implicit none

  integer ( kind = 4 ) hole_num

  hole_num = 1

  return
end
subroutine mesh09_hole_xy ( hole_num, hole_x, hole_y )

!*****************************************************************************80
!
!! MESH09_HOLE_XY returns hole coordinates of mesh #09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HOLE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) HOLE_X(HOLE_NUM), HOLE_Y(HOLE_NUM),
!    the XY coordinates a point in each hole.
!
  implicit none

  integer ( kind = 4 ) hole_num

  real ( kind = 8 ) hole_x(hole_num)
  real ( kind = 8 ) hole_y(hole_num)

  hole_x(1:hole_num) = (/ 0.500D+00 /)
  hole_y(1:hole_num) = (/ 0.444D+00 /)

  return
end
subroutine mesh09_name ( name )

!*****************************************************************************80
!
!! MESH09_NAME returns the name of mesh #09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the mesh.
!
  implicit none

  character ( len = * ) name

  name = 'The Pinched Hole Problem'

  return
end
subroutine mesh09_node_num ( node_num )

!*****************************************************************************80
!
!! MESH09_NODE_NUM returns the number of nodes of mesh #09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  node_num = 8

  return
end
subroutine mesh09_node_xy ( node_num, x, y )

!*****************************************************************************80
!
!! MESH09_NODE_XY returns the nodes of mesh #09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981,
!    ISBN: 0131652745,
!    LC: QA188.G46.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the XY
!    coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) y(node_num)

  x(1:node_num) = (/ &
    2.0D+00, 0.0D+00, 4.0D+00, 1.0D+00, 3.0D+00, &
    0.0D+00, 2.0D+00, 4.0D+00 /) / 4.0D+00

  y(1:node_num) = (/ &
    0.0D+00, 1.0D+00, 1.0D+00, 2.0D+00, 2.0D+00, &
    3.0D+00, 3.0D+00, 3.0D+00 /) / 3.0D+00

  return
end
subroutine node_eps ( file_name, node_num, node_mask, node_x, node_y, title )

!*****************************************************************************80
!
!! NODE_EPS creates an EPS file containing an image of the nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to create.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, logical NODE_MASK(NODE_NUM), is TRUE for those nodes to be plotted.
!
!    Input, real ( kind = 8 ) NODE_X(NODE_NUM), NODE_Y(NODE_NUM),
!    the coordinates of the nodes.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
  implicit none

  integer ( kind = 4 ) node_num

  integer ( kind = 4 ), parameter :: circle_size = 3
  real ( kind = 8 ) dif
  integer ( kind = 4 ) eps_unit
  integer ( kind = 4 ) eps_x
  integer ( kind = 4 ) eps_y
  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) node
  logical node_mask(node_num)
  real ( kind = 8 ) node_x(node_num)
  real ( kind = 8 ) node_x_max
  real ( kind = 8 ) node_x_min
  real ( kind = 8 ) node_y(node_num)
  real ( kind = 8 ) node_y_max
  real ( kind = 8 ) node_y_min
  real ( kind = 8 ) scale
  character ( len = 40 ) string
  character ( len = * ) title
!
!  Determine the range of the unmasked nodes.
!
  node_x_min =  huge ( node_x_min )
  node_x_max = -huge ( node_x_max )
  node_y_min =  huge ( node_y_min )
  node_y_max = -huge ( node_y_max )

  do node = 1, node_num
    if ( node_mask(node) ) then
      node_x_min = min ( node_x_min, node_x(node) )
      node_x_max = max ( node_x_max, node_x(node) )
      node_y_min = min ( node_y_min, node_y(node) )
      node_y_max = max ( node_y_max, node_y(node) )
    end if
  end do

  if ( node_y_max - node_y_min < node_x_max - node_x_min ) then
    scale = node_x_max - node_x_min
    dif = ( node_x_max - node_x_min ) - ( node_y_max - node_y_min )
    node_y_max = node_y_max + 0.5 * dif
    node_y_min = node_y_min - 0.5 * dif
  else
    scale = node_y_max - node_y_min
    dif = ( node_y_max - node_y_min ) - ( node_x_max - node_x_min )
    node_x_max = node_x_max + 0.5 * dif
    node_x_min = node_x_min - 0.5 * dif
  end if

  call get_unit ( eps_unit )

  open ( unit = eps_unit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NODE_EPS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output EPS file.'
    stop
  end if

  write ( eps_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( eps_unit, '(a)' ) '%%Creator: node_eps(aitch.f90)'
  write ( eps_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( eps_unit, '(a)' ) '%%Pages: 1'
  write ( eps_unit, '(a)' ) '%%BoundingBox:    36    36   576   756'
  write ( eps_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( eps_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( eps_unit, '(a)' ) '%%EndComments'
  write ( eps_unit, '(a)' ) '%%BeginProlog'
  write ( eps_unit, '(a)' ) '/inch {72 mul} def'
  write ( eps_unit, '(a)' ) '%%EndProlog'
  write ( eps_unit, '(a)' ) '%%Page:      1     1'
  write ( eps_unit, '(a)' ) 'save'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Set RGB line color.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.9000 0.9000 0.9000 setrgbcolor'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Draw a gray border around the page.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'newpath'
  write ( eps_unit, '(a)' ) '    36   126 moveto'
  write ( eps_unit, '(a)' ) '   576   126 lineto'
  write ( eps_unit, '(a)' ) '   576   666 lineto'
  write ( eps_unit, '(a)' ) '    36   666 lineto'
  write ( eps_unit, '(a)' ) '    36   126 lineto'
  write ( eps_unit, '(a)' ) 'stroke'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Set RGB line color.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Label the plot:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'
  write ( eps_unit, '(a)' ) '/Times-Roman findfont 0.50 inch scalefont setfont'
  write ( eps_unit, '(a)' ) '    36   666 moveto'
  write ( eps_unit, '(a)' ) '(' // trim ( title ) // ') show'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Define a clipping polygon'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '    36   126 moveto'
  write ( eps_unit, '(a)' ) '   576   126 lineto'
  write ( eps_unit, '(a)' ) '   576   666 lineto'
  write ( eps_unit, '(a)' ) '    36   666 lineto'
  write ( eps_unit, '(a)' ) '    36   126 lineto'
  write ( eps_unit, '(a)' ) 'clip newpath'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw filled dots at each node:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 1.0000 setrgbcolor'

  do node = 1, node_num

    if ( node_mask(node) ) then

      eps_x = int &
        ( ( node_x_max - node_x(node)              ) *  61.0   &
        + (            + node_x(node) - node_x_min ) * 551.0 ) &
        / scale

      eps_y = int &
        ( ( node_y_max - node_y(node)              ) * 151.0   &
        + (              node_y(node) - node_y_min ) * 641.0 ) &
        / scale

      write ( eps_unit, '(a,i4,2x,i4,2x,i4,a)' ) &
        'newpath  ', eps_x, eps_y, circle_size, ' 0 360 arc closepath fill'

    end if

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Label the nodes:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'
  write ( eps_unit, '(a)' ) '/Times-Roman findfont 0.20 inch scalefont setfont'

  do node = 1, node_num

    if ( node_mask(node) ) then

      eps_x = int &
        ( ( node_x_max - node_x(node)              ) *  61.0   &
        + (            + node_x(node) - node_x_min ) * 551.0 ) &
        / scale

      eps_y = int &
        ( ( node_y_max - node_y(node)              ) * 151.0   &
        + (              node_y(node) - node_y_min ) * 641.0 ) &
        / scale

      write ( string, '(i4)' ) node
      string = adjustl ( string )

      write ( eps_unit, '(i4,2x,i4,a)' ) eps_x, eps_y+5, &
        ' moveto (' // trim ( string ) // ') show'

    end if

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'restore showpage'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% End of page'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%%Trailer'
  write ( eps_unit, '(a)' ) '%%EOF'

  close ( unit = eps_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NODE_EPS:'
  write ( *, '(a)' ) '  An encapsulated PostScript file was created'
  write ( *, '(a)' ) '  containing an image of the nodes.'
  write ( *, '(a)' ) '  The file is named "' // trim ( file_name ) // '".'

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )
!
!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I precedes J, ISGN = +1 if J precedes I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I precedes J;
!    ISGN => 0 means J precedes I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0

  if ( n <= 1 ) then
    indx = 0
    return
  end if
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    n1 = n
    k = n / 2
    k1 = k
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i = i + 1
      end if

      j = k1
      k1 = i
      indx = - 1
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        indx = 0
      else
        i = n1
        n1 = n1 - 1
        j = 1
        indx = 1
      end if

      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i = 2 * k1

    if ( i == n1 ) then
      j = k1
      k1 = i
      indx = - 1
      return
    else if ( i <= n1 ) then
      j = i + 1
      indx = - 2
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    indx = 0
  else
    i = n1
    n1 = n1 - 1
    j = 1
    indx = 1
  end if

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
