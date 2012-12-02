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
!    CH_EQI ( 'A', 'a' ) is TRUE.
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

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

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
function ch_is_control ( c )

!*****************************************************************************80
!
!! CH_IS_CONTROL is TRUE if C is a control character.
!
!  Discussion:
!
!    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
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
!    Input, character C, the character to be tested.
!
!    Output, logical CH_IS_CONTROL, TRUE if C is a control character, and
!    FALSE otherwise.
!
  implicit none

  character c
  logical ch_is_control

  if ( ichar ( c ) <= 31 .or. 127 <= ichar ( c ) ) then
    ch_is_control = .true.
  else
    ch_is_control = .false.
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

  if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

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
subroutine obj_face_node_print ( face_num, order_max, face_order, face_node )

!*****************************************************************************80
!
!! OBJ_FACE_NODE_PRINT prints the node indices for each face.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices
!    per face.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Input, integer ( kind = 4 ) FACE_NODE(ORDER_MAX,FACE_NUM), the nodes that
!    make up each face.
!
  implicit none

  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_node(order_max,face_num)
  integer ( kind = 4 ) face_order(face_num)
  integer ( kind = 4 ) order
  character ( len = 26 ) string

  write ( string, '(a,i2,a)' ) '(2x,i6,2x,i6,2x,', order_max, '(2x,i6))'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Face   Order      Nodes'
  write ( *, '(a)' ) ' '

  do face = 1, face_num

    order = face_order(face)

    write ( *, string ) face, order, face_node(1:order,face)

  end do

  return
end
subroutine obj_normal_vector_print ( normal_num, normal_vector )

!*****************************************************************************80
!
!! OBJ_NORMAL_VECTOR_PRINT prints the normal vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NORMAL_NUM, the number of normal vectors.
!
!    Input, real ( kind = 8 ) NORMAL_VECTOR(3,NORMAL_NUM), the normal vectors.
!
  implicit none

  integer ( kind = 4 ) normal_num

  integer ( kind = 4 ) face
  integer ( kind = 4 ) normal
  real ( kind = 8 ) normal_vector(3,normal_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal Vectors:'
  write ( *, '(a)' ) ' '

  do normal = 1, normal_num

    write ( *, '(2x,i6,3(2x,g14.6))' ) normal, normal_vector(1:3,normal)

  end do

  return
end
subroutine obj_node_xyz_print ( node_num, node_xyz )

!*****************************************************************************80
!
!! OBJ_NODE_XYZ_PRINT prints the node coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the coordinates
!    of the nodes.
!
  implicit none

  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xyz(3,node_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Node         Coordinates'
  write ( *, '(a)' ) ' '

  do node = 1, node_num

    write ( *, '(2x,i6,3(2x,g14.6))' ) node, node_xyz(1:3,node)

  end do

  return
end
subroutine obj_read ( input_filename, node_num, face_num, normal_num, &
  order_max, node_xyz, face_order, face_node, normal_vector, vertex_normal )

!*****************************************************************************80
!
!! OBJ_READ reads graphics information from a Wavefront OBJ file.
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.
!
!    This is controlled by whether the input values have been zeroed
!    out or not.  This routine simply tacks on the information it
!    finds to the current graphics object.
!
!  Example:
!
!    #  magnolia.obj
!
!    v -3.269770 -39.572201 0.876128
!    v -3.263720 -39.507999 2.160890
!    ...
!    v 0.000000 -9.988540 0.000000
!    vn 1.0 0.0 0.0
!    ...
!    vn 0.0 1.0 0.0
!
!    f 8 9 11 10
!    f 12 13 15 14
!    ...
!    f 788 806 774
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) NORMAL_NUM, the number of normal vectors.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices
!    per face.
!
!    Output, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the coordinates of points.
!
!    Output, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Output, integer ( kind = 4 ) FACE_NODE(ORDER_MAX,FACE_NUM), the nodes
!    making faces.
!
!    Output, real ( kind = 8 ) NORMAL_VECTOR(3,NORMAL_NUM), normal vectors.
!
!    Output, integer ( kind = 4 ) VERTEX_NORMAL(ORDER_MAX,FACE_NUM), the indices
!    of normal vectors per vertex.
!
  implicit none

  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) normal_num
  integer ( kind = 4 ) order_max

  logical done
  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_node(order_max,face_num)
  integer ( kind = 4 ) face_order(face_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_file_unit
  integer ( kind = 4 ) input_file_status
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) lchar
  character ( len = 255 ) line
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xyz(3,node_num)
  integer ( kind = 4 ) normal
  real ( kind = 8 ) normal_vector(3,normal_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  logical s_eqi
  integer ( kind = 4 ) text_num
  real ( kind = 8 ) temp
  integer ( kind = 4 ) vertex
  integer ( kind = 4 ) vertex_normal(order_max,face_num)
  character ( len = 255 ) word
  integer ( kind = 4 ) word_index
  character ( len = 255 ) word_one

  ierror = 0

  face = 0
  node = 0
  normal = 0
  text_num = 0

  face_node(1:order_max,1:face_num) = 0
  face_order(1:face_num) = 0
  node_xyz(1:3,1:node_num) = 0.0D+00
  normal_vector(1:3,1:normal_num) = 0.0D+00
  vertex_normal(1:order_max,1:face_num) = 0

  call get_unit ( input_file_unit )

  open ( unit = input_file_unit, file = input_filename, status = 'old', &
    iostat = input_file_status )

  if ( input_file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OBJ_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' &
     // trim ( input_filename ) // '".'
    stop
  end if

  word = ' '
!
!  Read a line of text from the file.
!
  do

    read ( input_file_unit, '(a)', iostat = input_file_status ) line

    if ( input_file_status /= 0 ) then
      exit
    end if

    text_num = text_num + 1
!
!  Replace any control characters (in particular, TAB's) by blanks.
!
    call s_control_blank ( line )

    done = .true.
    word_index = 0
!
!  Read a word from the line.
!
    call word_next_read ( line, word, done )
!
!  If no more words in this line, read a new line.
!
    if ( done ) then
      cycle
    end if
!
!  If this word begins with '#' or '$', then it's a comment.  Read a new line.
!
    if ( word(1:1) == '#' .or. word(1:1) == '$' ) then
      cycle
    end if

    word_index = word_index + 1

    if ( word_index == 1 ) then
      word_one = word
    end if
!
!  BEVEL
!  Bevel interpolation.
!
    if ( s_eqi ( word_one, 'BEVEL' ) ) then
!
!  BMAT
!  Basis matrix.
!
    else if ( s_eqi ( word_one, 'BMAT' ) ) then
!
!  C_INTERP
!  Color interpolation.
!
    else if ( s_eqi ( word_one, 'C_INTERP' ) ) then
!
!  CON
!  Connectivity between free form surfaces.
!
    else if ( s_eqi ( word_one, 'CON' ) ) then
!
!  CSTYPE
!  Curve or surface type.
!
    else if ( s_eqi ( word_one, 'CSTYPE' ) ) then
!
!  CTECH
!  Curve approximation technique.
!
    else if ( s_eqi ( word_one, 'CTECH' ) ) then
!
!  CURV
!  Curve.
!
    else if ( s_eqi ( word_one, 'CURV' ) ) then
!
!  CURV2
!  2D curve.
!
    else if ( s_eqi ( word_one, 'CURV2' ) ) then
!
!  D_INTERP
!  Dissolve interpolation.
!
    else if ( s_eqi ( word_one, 'D_INTERP' ) ) then
!
!  DEG
!  Degree.
!
    else if ( s_eqi ( word_one, 'DEG' ) ) then
!
!  END
!  End statement.
!
    else if ( s_eqi ( word_one, 'END' ) ) then
!
!  F V1 V2 V3 ...
!    or
!  F V1/VT1/VN1 V2/VT2/VN2 ...
!    or
!  F V1//VN1 V2//VN2 ...
!
!  Face.
!  A face is defined by the vertices.
!  Optionally, slashes may be used to include the texture vertex
!  and vertex normal indices.
!
    else if ( s_eqi ( word_one, 'F' ) ) then

      face = face + 1

      vertex = 0

      do

        call word_next_read ( line, word, done )

        if ( done ) then
          exit
        end if

        vertex = vertex + 1
!
!  Locate the slash characters in the word, if any.
!
        i1 = index ( word, '/' )
        if ( 0 < i1 ) then
          i2 = index ( word(i1+1:), '/' ) + i1
        else
          i2 = 0
        end if
!
!  Read the vertex index.
!
        call s_to_i4 ( word, itemp, ierror, lchar )

        if ( ierror /= 0 ) then
          itemp = -1
          ierror = 0
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'OBJ_READ - Error!'
          write ( *, '(a)' ) '  Bad FACE field.'
          write ( *, '(a)' ) trim ( word )
        end if

        face_node(vertex,face) = itemp
        face_order(face) = face_order(face) + 1
!
!  If there are two slashes, then read the data following the second one.
!
        if ( 0 < i2 ) then

          call s_to_i4 ( word(i2+1:), itemp, ierror, lchar )

          vertex_normal(vertex,face) = itemp

        end if

      end do
!
!  G
!  Group name.
!
    else if ( s_eqi ( word_one, 'G' ) ) then
!
!  HOLE
!  Inner trimming loop.
!
    else if ( s_eqi ( word_one, 'HOLE' ) ) then
!
!  L
!  A line, described by a sequence of vertex indices.
!  Are the vertex indices 0 based or 1 based?
!
    else if ( s_eqi ( word_one, 'L' ) ) then
!
!  LOD
!  Level of detail.
!
    else if ( s_eqi ( word_one, 'LOD' ) ) then
!
!  MG
!  Merging group.
!
    else if ( s_eqi ( word_one, 'MG' ) ) then
!
!  MTLLIB
!  Material library.
!
    else if ( s_eqi ( word_one, 'MTLLIB' ) ) then
!
!  O
!  Object name.
!
    else if ( s_eqi ( word_one, 'O' ) ) then
!
!  P
!  Point.
!
    else if ( s_eqi ( word_one, 'P' ) ) then
!
!  PARM
!  Parameter values.
!
    else if ( s_eqi ( word_one, 'PARM' ) ) then
!
!  S
!  Smoothing group.
!
    else if ( s_eqi ( word_one, 'S' ) ) then
!
!  SCRV
!  Special curve.
!
    else if ( s_eqi ( word_one, 'SCRV' ) ) then
!
!  SHADOW_OBJ
!  Shadow casting.
!
    else if ( s_eqi ( word_one, 'SHADOW_OBJ' ) ) then
!
!  SP
!  Special point.
!
    else if ( s_eqi ( word_one, 'SP' ) ) then
!
!  STECH
!  Surface approximation technique.
!
    else if ( s_eqi ( word_one, 'STECH' ) ) then
!
!  STEP
!  Stepsize.
!
    else if ( s_eqi ( word_one, 'STEP' ) ) then
!
!  SURF
!  Surface.
!
    else if ( s_eqi ( word_one, 'SURF' ) ) then
!
!  TRACE_OBJ
!  Ray tracing.
!
    else if ( s_eqi ( word_one, 'TRACE_OBJ' ) ) then
!
!  TRIM
!  Outer trimming loop.
!
    else if ( s_eqi ( word_one, 'TRIM' ) ) then
!
!  USEMTL
!  Material name.
!
    else if ( s_eqi ( word_one, 'USEMTL' ) ) then
!
!  V X Y Z
!  Geometric vertex.
!
    else if ( s_eqi ( word_one, 'V' ) ) then

      node = node + 1

      do i = 1, 3
        call word_next_read ( line, word, done )
        call s_to_r8 ( word, temp, ierror, lchar )
        node_xyz(i,node) = temp
      end do
!
!  VN
!  Vertex normals.
!
    else if ( s_eqi ( word_one, 'VN' ) ) then

      normal = normal + 1

      do i = 1, 3
        call word_next_read ( line, word, done )
        call s_to_r8 ( word, temp, ierror, lchar )
        normal_vector(i,normal) = temp
      end do
!
!  VT
!  Vertex texture.
!
    else if ( s_eqi ( word_one, 'VT' ) ) then
!
!  VP
!  Parameter space vertices.
!
    else if ( s_eqi ( word_one, 'VP' ) ) then
!
!  Unrecognized keyword.
!
    else

    end if

  end do

  close ( unit = input_file_unit )

  return
end
subroutine obj_size ( input_filename, node_num, face_num, normal_num, &
  order_max )

!*****************************************************************************80
!
!! OBJ_SIZE determines sizes of graphics objects in an Alias OBJ file.
!
!  Discussion:
!
!    The only items of interest to this routine are vertices,
!    faces, and normal vectors.
!
!  Example:
!
!    #  magnolia.obj
!
!    v -3.269770 -39.572201 0.876128
!    v -3.263720 -39.507999 2.160890
!    ...
!    v 0.000000 -9.988540 0.000000
!
!    vn 1.0 0.0 0.0
!    ...
!    vn 0.0 1.0 0.0
!
!    f 8 9 11 10
!    f 12 13 15 14
!    ...
!    f 788 806 774
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 December 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the input file name.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) NORMAL_NUM, the number of normal vectors.
!
!    Output, integer ( kind = 4 ) ORDER_MAX, the maximum face order.
!
  implicit none

  logical done
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_file_status
  integer ( kind = 4 ) input_file_unit
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) lchar
  character ( len = 255 ) line
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) normal_num
  integer ( kind = 4 ) order_max
  logical s_eqi
  real ( kind = 8 ) temp
  integer ( kind = 4 ) text_num
  integer ( kind = 4 ) vertex
  character ( len = 255 ) word
  integer ( kind = 4 ) word_index
  character ( len = 255 ) word_one

  ierror = 0

  face_num = 0
  node_num = 0
  normal_num = 0
  order_max = 0
  text_num = 0

  call get_unit ( input_file_unit )

  open ( unit = input_file_unit, file = input_filename, status = 'old', &
    iostat = input_file_status )

  if ( input_file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OBJ_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' &
     // trim ( input_filename ) // '".'
    stop
  end if

  word = ' '
!
!  Read a line of text from the file.
!
  do

    read ( input_file_unit, '(a)', iostat = input_file_status ) line

    if ( input_file_status /= 0 ) then
      exit
    end if

    text_num = text_num + 1
!
!  Replace any control characters (in particular, TAB's) by blanks.
!
    call s_control_blank ( line )

    done = .true.
    word_index = 0
!
!  Read a word from the line.
!
    call word_next_read ( line, word, done )
!
!  If no more words in this line, read a new line.
!
    if ( done ) then
      cycle
    end if
!
!  If this word begins with '#' or '$', then it's a comment.  Read a new line.
!
    if ( word(1:1) == '#' .or. word(1:1) == '$' ) then
      cycle
    end if

    word_index = word_index + 1

    if ( word_index == 1 ) then
      word_one = word
    end if
!
!  F V1 V2 V3 ...
!    or
!  F V1/VT1/VN1 V2/VT2/VN2 ...
!    or
!  F V1//VN1 V2//VN2 ...
!
!  Face.
!  A face is defined by the vertices.
!  Optionally, slashes may be used to include the texture vertex
!  and vertex normal indices.
!
    if ( s_eqi ( word_one, 'F' ) ) then

      face_num = face_num + 1

      vertex = 0

      do

        call word_next_read ( line, word, done )

        if ( done ) then
          exit
        end if

        vertex = vertex + 1
        order_max = max ( order_max, vertex )
!
!  Locate the slash characters in the word, if any.
!
        i1 = index ( word, '/' )
        if ( 0 < i1 ) then
          i2 = index ( word(i1+1:), '/' ) + i1
        else
          i2 = 0
        end if
!
!  Read the vertex index.
!
        call s_to_i4 ( word, itemp, ierror, lchar )

        if ( ierror /= 0 ) then
          itemp = -1
          ierror = 0
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'OBJ_SIZE - Error!'
          write ( *, '(a)' ) '  Bad FACE field.'
          write ( *, '(a)' ) trim ( word )
        end if
!
!  If there are two slashes, then read the data following the second one.
!
        if ( 0 < i2 ) then
          call s_to_i4 ( word(i2+1:), itemp, ierror, lchar )
        end if

      end do
!
!  V X Y Z W
!  Geometric vertex.
!
    else if ( s_eqi ( word_one, 'V' ) ) then

      node_num = node_num + 1
      cycle
!
!  VN
!  Vertex normals.
!
    else if ( s_eqi ( word_one, 'VN' ) ) then

      normal_num = normal_num + 1
      cycle

    end if

  end do

  close ( unit = input_file_unit )

  return
end
subroutine obj_size_print ( filename, node_num, face_num, normal_num, &
  order_max )

!*****************************************************************************80
!
!! OBJ_SIZE_PRINT prints sizes associated with an OBJ file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of vertices defined.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces defined.
!
!    Input, integer ( kind = 4 ) NORMAL_NUM, the number of normal
!    vectors defined.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices
!    per face.
!
  implicit none

  integer ( kind = 4 ) face_num
  character ( len = * ) filename
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) normal_num
  integer ( kind = 4 ) order_max

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Object sizes for OBJ file "' // &
    trim ( filename ) // '":'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Nodes              = ', node_num
  write ( *, '(a,i6)' ) '  Faces              = ', face_num
  write ( *, '(a,i6)' ) '  Maximum face order = ', order_max
  write ( *, '(a,i6)' ) '  Normal vectors     = ', normal_num

  return
end
subroutine obj_vertex_normal_print ( order_max, face_num, face_order, &
  vertex_normal )

!*****************************************************************************80
!
!! OBJ_VERTEX_NORMAL_PRINT prints the normal vectors indices per vertex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices
!    per face.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Input, integer ( kind = 4 ) VERTEX_NORMAL(ORDER_MAX,FACE_NUM), the
!    indices of normal vectors per vertex.
!
  implicit none

  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_order(face_num)
  integer ( kind = 4 ) order
  character ( len = 26 ) string
  integer ( kind = 4 ) vertex_normal(order_max,face_num)

  write ( string, '(a,i2,a)' ) '(2x,i6,2x,i6,2x,', order_max, '(2x,i6))'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal Vector Indices:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Face   Order'
  write ( *, '(a)' ) ' '

  do face = 1, face_num

    order = face_order(face)

    write ( *, string ) face, order, vertex_normal(1:order,face)

  end do

  return
end
subroutine obj_write ( output_filename, node_num, face_num, normal_num, &
  order_max, node_xyz, face_order, face_node, normal_vector, vertex_normal )

!*****************************************************************************80
!
!! OBJ_WRITE writes graphics information to an Alias OBJ file.
!
!  Discussion:
!
!    If no normal vectors are supplied (NORMAL_NUM <= 0) then
!    a simple format is used for the "F" records.  Otherwise,
!    the "v//vn" format is used.
!
!  Example:
!
!    #  no_normals.obj
!
!    g Group002
!
!    v -3.269770 -39.572201 0.876128
!    v -3.263720 -39.507999 2.160890
!    ...
!    v 0.000000 -9.988540 0.000000
!
!    f 8 9 11 10
!    f 12 13 15 14
!    ...
!    f 788 806 774
!
!    #  normals_supplied.obj
!
!    g Group001
!
!    v -3.269770 -39.572201 0.876128
!    v -3.263720 -39.507999 2.160890
!    ...
!    v 0.000000 -9.988540 0.000000
!
!    vn 0.0 1.0 0.0
!    vn 1.0 0.0 0.0
!    ...
!    vn 0.0 0.0 1.0
!
!    f 8//1 9//2 11//3 10//4
!    f 12//5 13//6 15//7 14//8
!    ...
!    f 788//800 806//803 774//807
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) NORMAL_NUM, the number of normal vectors.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices
!    per face.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the coordinates of points.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Input, integer ( kind = 4 ) FACE_NODE(ORDER_MAX,FACE_NUM), the nodes
!    making faces.
!
!    Input, real ( kind = 8 ) NORMAL_VECTOR(3,NORMAL_NUM), normal vectors.
!
!    Input, integer ( kind = 4 ) VERTEX_NORMAL(ORDER_MAX,FACE_NUM), the
!    indices of normal vectors per vertex.
!
  implicit none

  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) normal_num
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_node(order_max,face_num)
  integer ( kind = 4 ) face_order(face_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xyz(3,node_num)
  integer ( kind = 4 ) normal
  real ( kind = 8 ) normal_vector(3,normal_num)
  integer ( kind = 4 ) output_file_unit
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_file_status
  character ( len = 255 ) text
  integer ( kind = 4 ) text_num
  character ( len = 255 ) text2
  integer ( kind = 4 ) vertex
  integer ( kind = 4 ) vertex_normal(order_max,face_num)
  real ( kind = 8 ) w

  call get_unit ( output_file_unit )

  open ( unit = output_file_unit, file = output_filename, &
    status = 'replace', iostat = output_file_status )

  if ( output_file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OBJ_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' &
      // trim ( output_filename ) // '".'
    return
  end if

  text_num = 0

  write ( output_file_unit, '(a)' ) '# ' // trim ( output_filename )
  write ( output_file_unit, '(a)' ) '# ' // 'created by OBJ_WRITE.'
  write ( output_file_unit, '(a)' ) ' '
  write ( output_file_unit, '(a)' ) 'g Group001'

  text_num = text_num + 4
!
!  V: vertex coordinates.
!  For some reason, a fourth "coordinate" may be recommended.
!  What is its meaning?
!
  if ( 0 < node_num ) then
    write ( output_file_unit, '(a)' ) ' '
    text_num = text_num + 1
  end if

  w = 1.0D+00
  do node = 1, node_num
    write ( text, '(a1,2x,4g14.6)' ) 'v', node_xyz(1:3,node), w
    call s_blanks_delete ( text )
    write ( output_file_unit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do
!
!  VN: normal vectors.
!
  if ( 0 < normal_num ) then

    write ( output_file_unit, '(a)' ) ' '
    text_num = text_num + 1

    do normal = 1, normal_num

      write ( text, '(a2,2x,3f7.3)' ) 'vn', normal_vector(1:3,normal)
      call s_blanks_delete ( text )
      write ( output_file_unit, '(a)' ) trim ( text )
      text_num = text_num + 1

    end do

  end if
!
!  F: Faces, specified as a list of triples, one triple for each vertex:
!    vertex index/vertex texture index/vertex normal index
!
  if ( 0 < face_num ) then
    write ( output_file_unit, '(a)' ) ' '
    text_num = text_num + 1
  end if

  do face = 1, face_num

    text = 'f'

    if ( normal_num <= 0 ) then

      do vertex = 1, face_order(face)
        text2 = ' '
        write ( text2(2:), '(i8)' ) face_node(vertex,face)
        call s_blank_delete ( text2(2:) )
        call s_cat ( text, text2, text )
      end do

    else

      do vertex = 1, face_order(face)
        text2 = ' '
        write ( text2(2:), '(i8, ''//'', i8 )' ) face_node(vertex,face), &
          vertex_normal(vertex,face)
        call s_blank_delete ( text2(2:) )
        call s_cat ( text, text2, text )
      end do

    end if

    write ( output_file_unit, '(a)' ) trim ( text )
    text_num = text_num + 1

  end do

  close ( unit = output_file_unit )
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'OBJ_WRITE:'
  write ( *, '(a,i6,a)' ) '  Wrote ', text_num, ' text lines to "' &
    // trim ( output_filename ) // '".'

  return
end
subroutine r8vec_cross_product_3d ( v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_CROSS_PRODUCT_3D computes the cross product of two vectors in 3D.
!
!  Discussion:
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), the two vectors.
!
!    Output, real ( kind = 8 ) V3(3), the cross product vector.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v2(dim_num)
  real ( kind = 8 ) v3(dim_num)

  v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

  return
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

  character c
  integer ( kind = 4 ) get
  integer ( kind = 4 ) put
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  put = 0
  nchar = len_trim ( s )

  do get = 1, nchar

    c = s(get:get)

    if ( c /= ' ' .and. c /= TAB ) then
      put = put + 1
      s(put:put) = c
    end if

  end do

  s(put+1:nchar) = ' '

  return
end
subroutine s_blanks_delete ( s )

!*****************************************************************************80
!
!! S_BLANKS_DELETE replaces consecutive blanks by one blank.
!
!  Discussion:
!
!    Thanks to Bill Richmond for pointing out a programming flaw which
!    meant that, as characters were slid to the left through multiple
!    blanks, their original images were not blanked out.  This problem
!    is easiest resolved by using a copy of the string.
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
!    30 September 2004
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

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character newchr
  character oldchr
  character ( len = * ) s
  integer ( kind = 4 ) s_len
  character, parameter :: TAB = char ( 9 )
  character ( len = len ( s ) ) s_copy

  s_len = len ( s )

  j = 0
  s_copy(1:s_len) = s(1:s_len)
  s(1:s_len) = ' '

  newchr = ' '

  do i = 1, s_len

    oldchr = newchr
    newchr = s_copy(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      s(j:j) = newchr
    end if

  end do

  return
end
subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT concatenates two strings to make a third string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  if ( s1 == ' ' .and. s2 == ' ' ) then
    s3 = ' '
  else if ( s1 == ' ' ) then
    s3 = s2
  else if ( s2 == ' ' ) then
    s3 = s1
  else
    s3 = trim ( s1 ) // trim ( s2 )
  end if

  return
end
subroutine s_control_blank ( s )

!*****************************************************************************80
!
!! S_CONTROL_BLANK replaces control characters with blanks.
!
!  Discussion:
!
!    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
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
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  logical ch_is_control
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) s

  nchar = len_trim ( s )

  do i = 1, nchar
    if ( ch_is_control ( s(i:i) ) ) then
      s(i:i) = ' '
    end if
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
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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
subroutine s_to_i4 ( s, ival, ierror, length )

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
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S
!    used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) length
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
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
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
!    07 September 2004
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
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
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
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
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
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

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
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
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
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length+1 == nchar ) then
    length = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
subroutine word_next_read ( s, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_READ "reads" words from a string, one at a time.
!
!  Discussion:
!
!    The following characters are considered to be a single word,
!    whether surrounded by spaces or not:
!
!      " ( ) { } [ ]
!
!    Also, if there is a trailing comma on the word, it is stripped off.
!    This is to facilitate the reading of lists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string, presumably containing words
!    separated by spaces.
!
!    Output, character ( len = * ) WORD.
!    If DONE is FALSE, then WORD contains the "next" word read.
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!    On input with a fresh string, set DONE to TRUE.
!    On output, the routine sets DONE:
!      FALSE if another word was read,
!      TRUE if no more words could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), save :: lenc = 0
  integer ( kind = 4 ), save :: next = 1
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
  character ( len = * ) word
!
!  We "remember" LENC and NEXT from the previous call.
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
  if ( done ) then

    next = 1
    done = .false.
    lenc = len_trim ( s )

    if ( lenc <= 0 ) then
      done = .true.
      word = ' '
      return
    end if

  end if
!
!  Beginning at index NEXT, search the string for the next nonblank,
!  which signals the beginning of a word.
!
  ilo = next
!
!  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  do

    if ( lenc < ilo ) then
      word = ' '
      done = .true.
      next = lenc + 1
      return
    end if
!
!  If the current character is blank, skip to the next one.
!
    if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
  if ( s(ilo:ilo) == '"' .or. &
       s(ilo:ilo) == '(' .or. &
       s(ilo:ilo) == ')' .or. &
       s(ilo:ilo) == '{' .or. &
       s(ilo:ilo) == '}' .or. &
       s(ilo:ilo) == '[' .or. &
       s(ilo:ilo) == ']' ) then

    word = s(ilo:ilo)
    next = ilo + 1
    return

  end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
  next = ilo + 1

  do while ( next <= lenc )

    if ( s(next:next) == ' ' ) then
      exit
    else if ( s(next:next) == TAB ) then
      exit
    else if ( s(next:next) == '"' ) then
      exit
    else if ( s(next:next) == '(' ) then
      exit
    else if ( s(next:next) == ')' ) then
      exit
    else if ( s(next:next) == '{' ) then
      exit
    else if ( s(next:next) == '}' ) then
      exit
    else if ( s(next:next) == '[' ) then
      exit
    else if ( s(next:next) == ']' ) then
      exit
    end if

    next = next + 1

  end do

  if ( s(next-1:next-1) == ',' ) then
    word = s(ilo:next-2)
  else
    word = s(ilo:next-1)
  end if

  return
end
