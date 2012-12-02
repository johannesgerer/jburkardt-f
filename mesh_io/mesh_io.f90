subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
!    which guarantee the ASCII collating sequence.
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
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  character              ch
  integer   ( kind = 4 ) itemp

  itemp = iachar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Discussion:
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
  logical   ch_eqi

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
subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!     CH  DIGIT
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
!    Input, character CH, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
!    If CH was 'illegal', then DIGIT is -1.
!
  implicit none

  character              ch
  integer   ( kind = 4 ) digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then

    digit = iachar ( ch ) - 48

  else if ( ch == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine cyl248_data ( dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, &
  edge_vertex, edge_label, triangle_vertex, triangle_label, &
  quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
  tetrahedron_label, hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! CYL248_DATA defines the data for a 3D tetrahedral mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pascal Frey,
!    MEDIT: An interactive mesh visualization software,
!    Technical Report RT-0253,
!    Institut National de Recherche en Informatique et en Automatique,
!    03 December 2001.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 2 or 3.
!
!    Input, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Input, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Input, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
!    Output, real ( kind = 8 ) VERTEX_COORDINATE(DIM,VERTICES), the
!    coordinates of each vertex.
!
!    Output, integer ( kind = 4 ) VERTEX_LABEL(VERTICES), a label for
!    each vertex.
!
!    Output, integer ( kind = 4 ) EDGE_VERTEX(2,EDGES), the vertices that form
!    each edge.
!
!    Output, integer ( kind = 4 ) EDGE_LABEL(EDGES), a label for each edge.
!
!    Output, integer ( kind = 4 ) TRIANGLE_VERTEX(3,TRIANGLES), the vertices
!    that form each triangle.
!
!    Output, integer ( kind = 4 ) TRIANGLE_LABEL(TRIANGLES), a label for
!    each triangle.
!
!    Output, integer ( kind = 4 ) QUADRILATERAL_VERTEX(4,QUADRILATERALS), the
!    vertices that form each quadrilateral.
!
!    Output, integer ( kind = 4 ) QUADRILATERAL_LABEL(QUADRILATERALS), a label
!    for each quadrilateral.
!
!    Output, integer ( kind = 4 ) TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
!    vertices that form each tetrahedron.
!
!    Output, integer ( kind = 4 ) TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
!    each tetrahedron.
!
!    Output, integer ( kind = 4 ) HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the
!    vertices that form each hexahedron.
!
!    Output, integer ( kind = 4 ) HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
!    each hexahedron.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ), parameter :: dim_save = 3
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ), parameter :: tetrahedrons_save = 248
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ), parameter :: triangles_save = 154
  integer ( kind = 4 ) vertices
  integer ( kind = 4 ), parameter :: vertices_save = 92

  integer ( kind = 4 ) edge_label(edges)
  integer ( kind = 4 ) edge_vertex(2,edges)
  integer ( kind = 4 ) hexahedron_label(hexahedrons)
  integer ( kind = 4 ) hexahedron_vertex(8,hexahedrons)
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals)
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  integer ( kind = 4 ), save :: tetrahedron_vertex_save(4,tetrahedrons_save) &
    = reshape ( (/ &
    23, 1, 9, 8, &
    27, 9, 23, 1, &
    26, 8, 23, 9, &
    26, 9, 7, 8, &
    2, 9, 27, 1, &
    26, 9, 10, 7, &
    26, 28, 7, 10, &
    11, 29, 3, 2, &
    7, 6, 10, 28, &
    10, 6, 31, 28, &
    11, 29, 30, 3, &
    11, 30, 4, 3, &
    11, 30, 32, 4, &
    10, 6, 5, 31, &
    11, 5, 4, 32, &
    19, 33, 34, 20, &
    39, 22, 40, 16, &
    39, 17, 36, 22, &
    39, 22, 16, 17, &
    40, 22, 15, 16, &
    12, 19, 20, 33, &
    19, 20, 34, 18, &
    12, 33, 20, 35, &
    38, 37, 14, 21, &
    36, 22, 17, 18, &
    38, 14, 15, 21, &
    13, 14, 37, 21, &
    12, 20, 13, 35, &
    80, 32, 11, 30, &
    80, 28, 10, 31, &
    80, 31, 59, 28, &
    80, 58, 57, 26, &
    80, 28, 58, 26, &
    80, 59, 58, 28, &
    80, 28, 26, 10, &
    80, 10, 26, 9, &
    80, 9, 11, 10, &
    80, 9, 26, 23, &
    80, 23, 26, 57, &
    80, 23, 27, 9, &
    80, 23, 56, 27, &
    80, 30, 11, 29, &
    80, 5, 10, 11, &
    80, 5, 11, 32, &
    80, 5, 32, 31, &
    80, 31, 10, 5, &
    80, 2, 11, 9, &
    80, 29, 11, 2, &
    80, 2, 9, 27, &
    80, 27, 29, 2, &
    81, 40, 39, 22, &
    81, 22, 39, 36, &
    81, 18, 36, 34, &
    81, 34, 20, 18, &
    81, 22, 36, 18, &
    81, 20, 22, 18, &
    81, 37, 38, 21, &
    81, 20, 33, 35, &
    81, 13, 21, 20, &
    81, 13, 20, 35, &
    81, 13, 37, 21, &
    81, 35, 37, 13, &
    81, 20, 21, 22, &
    81, 34, 33, 20, &
    81, 21, 38, 15, &
    81, 38, 40, 15, &
    81, 22, 21, 15, &
    81, 15, 40, 22, &
    82, 60, 74, 59, &
    82, 74, 25, 59, &
    82, 73, 72, 58, &
    82, 25, 73, 58, &
    82, 59, 25, 58, &
    82, 58, 72, 57, &
    82, 57, 80, 58, &
    82, 58, 80, 59, &
    83, 71, 79, 70, &
    83, 70, 76, 78, &
    83, 79, 76, 70, &
    83, 79, 60, 76, &
    83, 82, 60, 74, &
    84, 54, 64, 55, &
    84, 64, 65, 55, &
    84, 65, 63, 55, &
    84, 65, 71, 63, &
    85, 29, 62, 30, &
    85, 80, 29, 30, &
    85, 29, 61, 62, &
    85, 78, 83, 76, &
    85, 78, 76, 30, &
    85, 62, 78, 30, &
    85, 76, 83, 60, &
    85, 76, 32, 30, &
    85, 32, 80, 30, &
    85, 32, 76, 60, &
    85, 27, 61, 29, &
    85, 80, 27, 29, &
    85, 83, 82, 60, &
    85, 77, 78, 62, &
    85, 60, 82, 59, &
    85, 59, 82, 80, &
    85, 32, 60, 31, &
    85, 80, 32, 31, &
    85, 60, 59, 31, &
    85, 59, 80, 31, &
    86, 51, 68, 52, &
    86, 69, 68, 51, &
    86, 68, 67, 52, &
    86, 52, 67, 53, &
    86, 67, 66, 53, &
    86, 53, 66, 54, &
    87, 50, 70, 49, &
    87, 71, 70, 50, &
    87, 63, 71, 50, &
    87, 63, 84, 71, &
    87, 70, 69, 49, &
    87, 71, 83, 70, &
    87, 49, 69, 51, &
    87, 69, 86, 51, &
    88, 64, 66, 73, &
    88, 72, 73, 66, &
    88, 72, 82, 73, &
    88, 24, 72, 66, &
    88, 64, 73, 25, &
    88, 73, 82, 25, &
    88, 66, 64, 54, &
    88, 84, 54, 64, &
    88, 87, 86, 84, &
    88, 67, 24, 66, &
    88, 66, 86, 67, &
    88, 64, 25, 65, &
    88, 65, 84, 64, &
    88, 25, 74, 65, &
    88, 25, 82, 74, &
    88, 83, 87, 71, &
    88, 71, 87, 84, &
    88, 82, 83, 74, &
    88, 74, 83, 71, &
    88, 65, 74, 71, &
    88, 71, 84, 65, &
    89, 86, 87, 84, &
    89, 39, 48, 44, &
    89, 44, 49, 43, &
    89, 44, 43, 36, &
    89, 44, 48, 50, &
    89, 48, 63, 50, &
    89, 86, 84, 54, &
    89, 51, 87, 86, &
    89, 44, 50, 49, &
    89, 50, 87, 49, &
    89, 43, 49, 51, &
    89, 49, 87, 51, &
    89, 39, 44, 36, &
    89, 36, 81, 39, &
    89, 63, 48, 47, &
    89, 47, 48, 40, &
    89, 46, 55, 47, &
    89, 38, 46, 47, &
    89, 55, 63, 47, &
    89, 55, 84, 63, &
    89, 43, 42, 34, &
    89, 43, 51, 42, &
    89, 45, 53, 54, &
    89, 53, 86, 54, &
    89, 45, 54, 46, &
    89, 42, 52, 41, &
    89, 41, 52, 53, &
    89, 52, 86, 53, &
    89, 42, 51, 52, &
    89, 51, 86, 52, &
    89, 46, 54, 55, &
    89, 54, 84, 55, &
    90, 56, 75, 61, &
    90, 24, 75, 56, &
    90, 27, 56, 61, &
    90, 61, 85, 27, &
    90, 75, 77, 61, &
    90, 80, 82, 57, &
    90, 85, 82, 80, &
    90, 57, 24, 56, &
    90, 72, 24, 57, &
    90, 57, 82, 72, &
    90, 80, 56, 27, &
    90, 85, 80, 27, &
    91, 85, 90, 77, &
    91, 86, 87, 69, &
    91, 78, 77, 69, &
    91, 83, 88, 82, &
    91, 90, 82, 88, &
    91, 67, 88, 86, &
    91, 88, 87, 86, &
    91, 87, 88, 83, &
    91, 83, 85, 78, &
    91, 78, 85, 77, &
    91, 77, 75, 68, &
    91, 77, 90, 75, &
    91, 69, 77, 68, &
    91, 68, 86, 69, &
    91, 68, 75, 67, &
    91, 67, 86, 68, &
    91, 24, 88, 67, &
    91, 90, 88, 24, &
    91, 69, 87, 70, &
    91, 87, 83, 70, &
    91, 75, 24, 67, &
    91, 75, 90, 24, &
    92, 89, 46, 45, &
    92, 41, 53, 45, &
    92, 89, 45, 53, &
    92, 89, 53, 41, &
    92, 89, 41, 42, &
    92, 35, 41, 45, &
    92, 33, 41, 35, &
    92, 35, 81, 33, &
    92, 35, 45, 37, &
    92, 81, 35, 37, &
    92, 34, 89, 42, &
    92, 81, 89, 34, &
    92, 33, 42, 41, &
    92, 37, 45, 46, &
    92, 37, 46, 38, &
    92, 81, 37, 38, &
    92, 33, 34, 42, &
    92, 33, 81, 34, &
    83, 74, 60, 71, &
    83, 60, 79, 71, &
    89, 39, 40, 48, &
    89, 39, 81, 40, &
    89, 36, 43, 34, &
    89, 34, 81, 36, &
    89, 63, 87, 50, &
    89, 84, 87, 63, &
    54, 88, 66, 86, &
    54, 88, 86, 84, &
    90, 72, 88, 24, &
    90, 82, 88, 72, &
    38, 47, 89, 40, &
    38, 89, 81, 40, &
    92, 46, 89, 38, &
    92, 89, 81, 38, &
    80, 23, 57, 56, &
    80, 57, 90, 56, &
    61, 85, 62, 77, &
    61, 90, 85, 77, &
    82, 85, 91, 83, &
    82, 90, 91, 85, &
    70, 91, 78, 83, &
    70, 78, 91, 69 /), (/ 4, tetrahedrons_save /) )
  integer ( kind = 4 ) triangle_label(triangles)
  integer ( kind = 4 ), save :: triangle_label_save(triangles_save) &
    = (/ &
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
    4, 4, 2, 2, 2, 2, 2, 2, 2, 2, &
    2, 2, 2, 2, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3 /)
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  integer ( kind = 4 ), save &
    :: triangle_vertex_save(3,triangles_save) = reshape ( (/ &
    12, 20, 19, &
    12, 13, 20, &
    19, 20, 18, &
    20, 22, 18, &
    22, 17, 18, &
    22, 16, 17, &
    13, 21, 20, &
    13, 14, 21, &
    14, 15, 21, &
    22, 15, 16, &
    22, 21, 15, &
    20, 21, 22, &
    1, 9, 8, &
    2, 9, 1, &
    9, 7, 8, &
    2, 11, 9, &
    11, 2, 3, &
    11, 3, 4, &
    9, 10, 7, &
    7, 10, 6, &
    10, 5, 6, &
    11, 4, 5, &
    5, 10, 11, &
    9, 11, 10, &
    23, 1, 8, &
    26, 23, 8, &
    26, 8, 7, &
    27, 1, 23, &
    2, 1, 27, &
    26, 7, 28, &
    7, 6, 28, &
    27, 29, 2, &
    29, 3, 2, &
    29, 30, 3, &
    30, 4, 3, &
    6, 31, 28, &
    6, 5, 31, &
    5, 32, 31, &
    5, 4, 32, &
    12, 19, 33, &
    19, 34, 33, &
    19, 18, 34, &
    12, 33, 35, &
    12, 35, 13, &
    18, 36, 34, &
    36, 18, 17, &
    35, 37, 13, &
    13, 37, 14, &
    38, 14, 37, &
    38, 15, 14, &
    39, 36, 17, &
    39, 17, 16, &
    38, 40, 15, &
    40, 16, 15, &
    39, 16, 40, &
    33, 41, 35, &
    33, 42, 41, &
    33, 34, 42, &
    36, 43, 34, &
    43, 42, 34, &
    39, 44, 36, &
    44, 43, 36, &
    35, 45, 37, &
    35, 41, 45, &
    37, 46, 38, &
    37, 45, 46, &
    38, 47, 40, &
    38, 46, 47, &
    39, 48, 44, &
    39, 40, 48, &
    47, 48, 40, &
    44, 49, 43, &
    44, 50, 49, &
    44, 48, 50, &
    43, 51, 42, &
    43, 49, 51, &
    42, 52, 41, &
    42, 51, 52, &
    41, 53, 45, &
    41, 52, 53, &
    45, 54, 46, &
    45, 53, 54, &
    46, 55, 47, &
    46, 54, 55, &
    30, 32, 4, &
    23, 56, 27, &
    23, 57, 56, &
    23, 26, 57, &
    28, 58, 26, &
    58, 57, 26, &
    31, 59, 28, &
    59, 58, 28, &
    32, 60, 31, &
    60, 59, 31, &
    27, 61, 29, &
    27, 56, 61, &
    29, 62, 30, &
    29, 61, 62, &
    55, 63, 47, &
    63, 48, 47, &
    48, 63, 50, &
    54, 64, 55, &
    64, 65, 55, &
    65, 63, 55, &
    53, 66, 54, &
    66, 64, 54, &
    52, 67, 53, &
    67, 66, 53, &
    51, 68, 52, &
    68, 67, 52, &
    49, 69, 51, &
    69, 68, 51, &
    50, 70, 49, &
    70, 69, 49, &
    63, 71, 50, &
    71, 70, 50, &
    65, 71, 63, &
    64, 25, 65, &
    64, 73, 25, &
    64, 66, 73, &
    67, 24, 66, &
    24, 72, 66, &
    72, 73, 66, &
    68, 75, 67, &
    75, 24, 67, &
    69, 77, 68, &
    77, 75, 68, &
    70, 78, 69, &
    78, 77, 69, &
    62, 78, 30, &
    78, 76, 30, &
    76, 32, 30, &
    32, 76, 60, &
    61, 77, 62, &
    77, 78, 62, &
    56, 75, 61, &
    75, 77, 61, &
    57, 24, 56, &
    24, 75, 56, &
    58, 72, 57, &
    72, 24, 57, &
    59, 25, 58, &
    25, 73, 58, &
    73, 72, 58, &
    60, 74, 59, &
    74, 25, 59, &
    25, 74, 65, &
    65, 74, 71, &
    70, 76, 78, &
    71, 79, 70, &
    79, 76, 70, &
    79, 60, 76, &
    74, 60, 71, &
    60, 79, 71 /), (/ 3, triangles_save /) )
  real ( kind = 8 ) vertex_coordinate(dim,vertices)
  real ( kind = 8 ), save :: vertex_coordinate_save(dim_save,vertices_save) &
    = reshape ( (/ &
  1.0,       0.2,         0.0, &
  1.0,       0.141421,    0.141421, &
  1.0,       0.0,         0.2, &
  1.0,      -0.141421,    0.141421, &
  1.0,      -0.2,         0.0, &
  1.0,      -0.141421,   -0.141421, &
  1.0,       0.0,        -0.2, &
  1.0,       0.141421,   -0.141421, &
  1.0,       0.066163,   -0.0302872, &
  1.0,      -0.0615154,  -0.0610739, &
  1.0,      -0.0306985,   0.0668017, &
  0.0,       0.2,         0.0, &
  0.0,       0.141421,   -0.141421, &
  0.0,       0.0,        -0.2, &
  0.0,      -0.141421,   -0.141421, &
  0.0,      -0.2,         0.0, &
  0.0,      -0.141421,    0.141421, &
  0.0,       0.0,         0.2, &
  0.0,       0.141421,    0.141421, &
  0.0,       0.0686748,   0.0255359, &
  0.0,       0.0,        -0.0865993, &
  0.0,      -0.0686749,   0.0255359, &
  0.8816,    0.185522,   -0.0747102, &
  0.642415,  0.187806,   -0.0687668, &
  0.627606, -0.0696445,  -0.187482, &
  0.876431,  0.0811908,  -0.182779, &
  0.881613,  0.186118,    0.0732131, &
  0.872048, -0.0699008,  -0.187387, &
  0.878318,  0.0844232,   0.181308, &
  0.845861, -0.0716063,   0.186742, &
  0.866503, -0.182493,   -0.0818307, &
  0.859402, -0.186751,    0.0715813, &
  0.131355,  0.18477,     0.0765501, &
  0.13317,   0.077694,    0.184292, &
  0.130862,  0.185301,   -0.0752567, &
  0.135181, -0.0749468,   0.185426, &
  0.130839,  0.0781729,  -0.18409, &
  0.131856, -0.0754694,  -0.185214, &
  0.135683, -0.184121,    0.0780993, &
  0.134207, -0.184959,   -0.0760928, &
  0.261923,  0.199982,    0.00264585, &
  0.263928,  0.144161,    0.138627, &
  0.268645,  0.00535339,  0.199928, &
  0.272346, -0.137646,    0.145098, &
  0.26108,   0.144683,   -0.138082, &
  0.260772,  0.00498797, -0.199938, &
  0.264253, -0.139152,   -0.143655, &
  0.270288, -0.199962,    0.00389323, &
  0.408181, -0.0730357,   0.186187, &
  0.411818, -0.184374,    0.0774991, &
  0.397539,  0.080738,    0.182979, &
  0.39192,   0.185619,    0.0744699, &
  0.392192,  0.184438,   -0.0773479, &
  0.389194,  0.0770141,  -0.184577, &
  0.38786,  -0.0747817,  -0.185493, &
  0.762413,  0.199986,   -0.0023425, &
  0.762987,  0.151152,   -0.13097, &
  0.741526,  0.0187858,  -0.199116, &
  0.746899, -0.128364,   -0.153371, &
  0.720076, -0.19917,    -0.0182053, &
  0.7628,    0.152219,    0.129728, &
  0.763882,  0.0434475,   0.195224, &
  0.399903, -0.1841,     -0.0781489, &
  0.506331, -0.00579066, -0.199916, &
  0.514514, -0.133894,   -0.148568, &
  0.526121,  0.135152,   -0.147424, &
  0.517967,  0.199953,   -0.0043215, &
  0.520585,  0.147847,    0.13469, &
  0.533956,  0.0124181,   0.199614, &
  0.558316, -0.136902,    0.145801, &
  0.549126, -0.199624,   -0.0122659, &
  0.657307,  0.117735,   -0.161674, &
  0.611189,  0.041829,   -0.195577, &
  0.631917, -0.164669,   -0.113508, &
  0.641444,  0.187001,    0.0709267, &
  0.720251, -0.155557,    0.125706, &
  0.647345,  0.0932963,   0.176906, &
  0.677484, -0.0430068,   0.195321, &
  0.635293, -0.188734,    0.0661777, &
  0.888023, -0.00868364, -0.00818647, &
  0.112146,  0.0,        -0.0118425, &
  0.676228,  0.0124197,  -0.0856487, &
  0.638436, -0.0639898,   0.0525795, &
  0.452586, -0.0410297,  -0.0704842, &
  0.762004, -0.0188614,   0.0693717, &
  0.463368,  0.0649048,   0.0262133, &
  0.473921, -0.0356443,   0.0388516, &
  0.557002,  0.0123705,  -0.0932599, &
  0.290986, -0.0200898,   0.00857934, &
  0.7038,    0.0856777,   0.0182744, &
  0.576134,  0.0436218,   0.0828782, &
  0.215187,  0.080855,   -0.0314946 /), (/ 3, vertices_save /) )
  integer ( kind = 4 ) vertex_label(vertices)
  integer ( kind = 4 ), save :: vertex_label_save(vertices_save) = (/ &
    3, 3, 3, 3, 3, 3, 3, 3, 2, 2, &
    2, 3, 3, 3, 3, 3, 3, 3, 3, 4, &
    4, 4, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0 /)

  call r8mat_copy ( dim, vertices, vertex_coordinate_save, vertex_coordinate )
  call i4vec_copy ( vertices, vertex_label_save, vertex_label )
  call i4mat_copy ( 3, triangles, triangle_vertex_save, &
    triangle_vertex )
  call i4vec_copy ( triangles, triangle_label_save, &
    triangle_label )
  call i4mat_copy ( 4, tetrahedrons, tetrahedron_vertex_save, &
    tetrahedron_vertex )

  tetrahedron_label(1:tetrahedrons) = 1

  return
end
subroutine cyl248_size ( dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons )

!*****************************************************************************80
!
!! CYL248_SIZE defines the sizes for a 3D tetrahedral mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pascal Frey,
!    MEDIT: An interactive mesh visualization software,
!    Technical Report RT-0253,
!    Institut National de Recherche en Informatique et en Automatique,
!    03 December 2001.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM, the spatial dimension, which should be 2 or 3.
!
!    Output, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Output, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Output, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Output, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Output, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Output, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  dim = 3
  vertices = 92
  edges = 0
  triangles = 154
  quadrilaterals = 0
  tetrahedrons = 248
  hexahedrons = 0

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
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
  logical              lopen

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
subroutine hexahexa_2x2x2_data ( dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, &
  edge_vertex, edge_label, triangle_vertex, triangle_label, &
  quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
  tetrahedron_label, hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! HEXAHEXA_2X2X2_DATA defines the data for a 3D hexahedral mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pascal Frey,
!    MEDIT: An interactive mesh visualization software,
!    Technical Report RT-0253,
!    Institut National de Recherche en Informatique et en Automatique,
!    03 December 2001.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 2 or 3.
!
!    Input, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Input, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Input, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
!    Output, real ( kind = 8 ) VERTEX_COORDINATE(DIM,VERTICES), the
!    coordinates of each vertex.
!
!    Output, integer ( kind = 4 ) VERTEX_LABEL(VERTICES), a label for
!    each vertex.
!
!    Output, integer ( kind = 4 ) EDGE_VERTEX(2,EDGES), the vertices that form
!    each edge.
!
!    Output, integer ( kind = 4 ) EDGE_LABEL(EDGES), a label for each edge.
!
!    Output, integer ( kind = 4 ) TRIANGLE_VERTEX(3,TRIANGLES), the vertices
!    that form each triangle.
!
!    Output, integer ( kind = 4 ) TRIANGLE_LABEL(TRIANGLES), a label for
!    each triangle.
!
!    Output, integer ( kind = 4 ) QUADRILATERAL_VERTEX(4,QUADRILATERALS), the
!    vertices that form each quadrilateral.
!
!    Output, integer ( kind = 4 ) QUADRILATERAL_LABEL(QUADRILATERALS), a label
!    for each quadrilateral.
!
!    Output, integer ( kind = 4 ) TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
!    vertices that form each tetrahedron.
!
!    Output, integer ( kind = 4 ) TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
!    each tetrahedron.
!
!    Output, integer ( kind = 4 ) HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the
!    vertices that form each hexahedron.
!
!    Output, integer ( kind = 4 ) HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
!    each hexahedron.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ), parameter :: dim_save = 3
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ), parameter :: hexahedrons_save = 8
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ), parameter :: quadrilaterals_save = 24
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices
  integer ( kind = 4 ), parameter :: vertices_save = 27

  integer ( kind = 4 ) edge_label(edges)
  integer ( kind = 4 ) edge_vertex(2,edges)
  integer ( kind = 4 ) hexahedron_label(hexahedrons)
  integer ( kind = 4 ), save :: hexahedron_label_save(hexahedrons_save) = (/ &
    1, 1, 1, 1, 1, 1, 1, 1 /)
  integer ( kind = 4 ) hexahedron_vertex(8,hexahedrons)
  integer ( kind = 4 ), save :: hexahedron_vertex_save(8,hexahedrons_save) &
    = reshape ( (/ &
      1,  2,  5,  4, 10, 11, 14, 13, &
      2,  3,  6,  5, 11, 12, 15, 14, &
      4,  5,  8,  7, 13, 14, 17, 16, &
      5,  6,  9,  8, 14, 15, 18, 17, &
     10, 11, 14, 13, 19, 20, 23, 22, &
     11, 12, 15, 14, 20, 21, 24, 23, &
     13, 14, 17, 16, 22, 23, 26, 25, &
     14, 15, 18, 17, 23, 24, 27, 26 /), (/ 8, hexahedrons_save /) )
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ), save :: quadrilateral_label_save(quadrilaterals_save) &
    = (/ &
     1, 1, 1, 1, 2, 2, 2, 2, 3, 3, &
     3, 3, 4, 4, 4, 4, 5, 5, 5, 5, &
     6, 6, 6, 6 /)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals)
  integer ( kind = 4 ), save &
    :: quadrilateral_vertex_save(4,quadrilaterals_save) = reshape ( (/ &
     1,  4,  5,  2, &
     2,  5,  6,  3, &
     4,  7,  8,  5, &
     5,  8,  9,  6, &
     1,  2, 11, 10, &
     2,  3, 12, 11, &
    10, 11, 20, 19, &
    11, 12, 21, 20, &
     3,  6, 15, 12, &
     6,  9, 18, 15, &
    12, 15, 24, 21, &
    15, 18, 27, 24, &
     7, 16, 17,  8, &
     8, 17, 18,  9, &
    16, 25, 26, 17, &
    17, 26, 27, 18, &
     1, 10, 13,  4, &
     4, 13, 16,  7, &
    10, 19, 22, 13, &
    13, 22, 25, 16, &
    19, 20, 23, 22, &
    20, 21, 24, 23, &
    22, 23, 26, 25, &
    23, 24, 27, 26 /), (/ 4, quadrilaterals_save /) )
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  integer ( kind = 4 ) triangle_label(triangles)
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  real ( kind = 8 ) vertex_coordinate(dim,vertices)
  real ( kind = 8 ), save :: vertex_coordinate_save(dim_save,vertices_save) &
    = reshape ( (/ &
    0.0, 0.0, 0.0, &
    0.5, 0.0, 0.0, &
    1.0, 0.0, 0.0, &
    0.0, 0.5, 0.0, &
    0.5, 0.5, 0.0, &
    1.0, 0.5, 0.0, &
    0.0, 1.0, 0.0, &
    0.5, 1.0, 0.0, &
    1.0, 1.0, 0.0, &
    0.0, 0.0, 0.5, &
    0.5, 0.0, 0.5, &
    1.0, 0.0, 0.5, &
    0.0, 0.5, 0.5, &
    0.5, 0.5, 0.5, &
    1.0, 0.5, 0.5, &
    0.0, 1.0, 0.5, &
    0.5, 1.0, 0.5, &
    1.0, 1.0, 0.5, &
    0.0, 0.0, 1.0, &
    0.5, 0.0, 1.0, &
    1.0, 0.0, 1.0, &
    0.0, 0.5, 1.0, &
    0.5, 0.5, 1.0, &
    1.0, 0.5, 1.0, &
    0.0, 1.0, 1.0, &
    0.5, 1.0, 1.0, &
    1.0, 1.0, 1.0 /), (/ 3, vertices_save /) )
  integer ( kind = 4 ) vertex_label(vertices)
  integer ( kind = 4 ), save :: vertex_label_save(vertices_save) = (/ &
    5, 2, 3, 5, 1, 3, 5, 4, 4, 5, &
    2, 3, 5, 0, 3, 5, 4, 4, 6, 6, &
    6, 6, 6, 6, 6, 6, 6 /)

  call r8mat_copy ( dim, vertices, vertex_coordinate_save, vertex_coordinate )
  call i4vec_copy ( vertices, vertex_label_save, vertex_label )
  call i4mat_copy ( 4, quadrilaterals, quadrilateral_vertex_save, &
    quadrilateral_vertex )
  call i4vec_copy ( quadrilaterals, quadrilateral_label_save, &
    quadrilateral_label )
  call i4mat_copy ( 8, hexahedrons, hexahedron_vertex_save, hexahedron_vertex )
  call i4vec_copy ( hexahedrons, hexahedron_label_save, hexahedron_label )

  return
end
subroutine hexahexa_2x2x2_size ( dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons )

!*****************************************************************************80
!
!! HEXAHEXA_2X2X2_SIZE defines the sizes for a 3D hexahedral mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pascal Frey,
!    MEDIT: An interactive mesh visualization software,
!    Technical Report RT-0253,
!    Institut National de Recherche en Informatique et en Automatique,
!    03 December 2001.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM, the spatial dimension, which should be 2 or 3.
!
!    Output, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Output, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Output, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Output, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Output, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Output, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  dim = 3
  vertices = 27
  edges = 0
  triangles = 0
  quadrilaterals = 24
  tetrahedrons = 0
  hexahedrons = 8

  return
end
subroutine i4mat_copy ( m, n, a1, a2 )

!*****************************************************************************80
!
!! I4MAT_COPY copies an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A1(M,N), the matrix to be copied.
!
!    Output, integer ( kind = 4 ) A2(M,N), the copied matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(m,n)
  integer ( kind = 4 ) a2(m,n)

  a2(1:m,1:n) = a1(1:m,1:n)

  return
end
subroutine i4vec_copy ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC_COPY copies an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) A1(N), the vector to be copied.
!
!    Output, integer ( kind = 4 ) A2(N), a copy of A1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)

  a2(1:n) = a1(1:n)

  return
end
subroutine mesh_data_print ( dim, vertices, edges, triangles, quadrilaterals, &
  tetrahedrons, hexahedrons, vertex_coordinate, vertex_label,  edge_vertex, &
  edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, &
  quadrilateral_label, tetrahedron_vertex, tetrahedron_label, &
  hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! MESH_DATA_PRINT prints mesh data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pascal Frey,
!    MEDIT: An interactive mesh visualization software,
!    Technical Report RT-0253,
!    Institut National de Recherche en Informatique et en Automatique,
!    03 December 2001.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 2 or 3.
!
!    Input, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Input, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Input, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
!    Input, real VERTEX_COORDINATE(DIM,VERTICES), the coordinates
!    of each vertex.
!
!    Input, integer ( kind = 4 ) VERTEX_LABEL(VERTICES), a label for
!    each vertex.
!
!    Input, integer ( kind = 4 ) EDGE_VERTEX(2,EDGES), the vertices that form
!    each edge.
!
!    Input, integer ( kind = 4 ) EDGE_LABEL(EDGES), a label for each edge.
!
!    Input, integer ( kind = 4 ) TRIANGLE_VERTEX(3,TRIANGLES), the vertices
!    that form each triangle.
!
!    Input, integer ( kind = 4 ) TRIANGLE_LABEL(TRIANGLES), a label for each
!    triangle.
!
!    Input, integer ( kind = 4 ) QUADRILATERAL_VERTEX(4,QUADRILATERALS), the
!    vertices that form each quadrilateral.
!
!    Input, integer ( kind = 4 ) QUADRILATERAL_LABEL(QUADRILATERALS), a label
!    for each quadrilateral.
!
!    Input, integer ( kind = 4 ) TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
!    vertices that form each tetrahedron.
!
!    Input, integer ( kind = 4 ) TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
!    each tetrahedron.
!
!    Input, integer ( kind = 4 ) HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the vertices
!    that form each hexahedron.
!
!    Input, integer ( kind = 4 ) HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
!    each hexahedron.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  integer ( kind = 4 ) edge_label(edges)
  integer ( kind = 4 ) edge_vertex(2,edges)
  integer ( kind = 4 ) hexahedron_label(hexahedrons)
  integer ( kind = 4 ) hexahedron_vertex(8,hexahedrons)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals)
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  integer ( kind = 4 ) triangle_label(triangles)
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  real    ( kind = 8 ) vertex_coordinate(dim,vertices)
  integer ( kind = 4 ) vertex_label(vertices)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertices:'
  write ( *, '(a)' ) ' '
  if ( dim == 2 ) then
    do j = 1, vertices
      write ( *, '(2(2x,f10.4),2x,''('',i4,'')'')' ) &
        vertex_coordinate(1:dim,j), vertex_label(j)
    end do
  else
    do j = 1, vertices
      write ( *, '(3(2x,f10.4),2x,''('',i4,'')'')' ) &
        vertex_coordinate(1:dim,j), vertex_label(j)
    end do
  end if

  if ( 0 < edges ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Edges:'
    write ( *, '(a)' ) ' '
    do j = 1, edges
      write ( *, '(2(2x,i8),2x,''('',i4,'')'')' ) &
        edge_vertex(1:2,j), edge_label(j)
    end do
  end if

  if ( 0 < triangles ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Triangles:'
    write ( *, '(a)' ) ' '
    do j = 1, triangles
      write ( *, '(3(2x,i8),2x,''('',i4,'')'')' ) &
        triangle_vertex(1:3,j), triangle_label(j)
    end do
  end if

  if ( 0 < quadrilaterals ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Quadrilaterals:'
    write ( *, '(a)' ) ' '
    do j = 1, quadrilaterals
      write ( *, '(4(2x,i8),2x,''('',i4,'')'')' ) &
        quadrilateral_vertex(1:4,j), quadrilateral_label(j)
    end do
  end if

  if ( 0 < tetrahedrons ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Tetrahedrons:'
    write ( *, '(a)' ) ' '
    do j = 1, tetrahedrons
      write ( *, '(4(2x,i8),2x,''('',i4,'')'')' ) &
        tetrahedron_vertex(1:4,j), tetrahedron_label(j)
    end do
  end if

  if ( 0 < hexahedrons ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Hexahedrons:'
    write ( *, '(a)' ) ' '
    do j = 1, hexahedrons
      write ( *, '(8(2x,i8),2x,''('',i4,'')'')' ) &
        hexahedron_vertex(1:8,j), hexahedron_label(j)
    end do
  end if

  return
end
subroutine mesh_data_read ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
  vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
  quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
  tetrahedron_label, hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! MESH_DATA_READ reads data from a MESH file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pascal Frey,
!    MEDIT: An interactive mesh visualization software,
!    Technical Report RT-0253,
!    Institut National de Recherche en Informatique et en Automatique,
!    03 December 2001.
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the MESH file.
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should
!    be 2 or 3.
!
!    Input, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Input, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Input, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
!    Output, real VERTEX_COORDINATE(DIM,VERTICES), the coordinates
!    of each vertex.
!
!    Output, integer ( kind = 4 ) VERTEX_LABEL(VERTICES), a label for
!    each vertex.
!
!    Output, integer ( kind = 4 ) EDGE_VERTEX(2,EDGES), the vertices that form
!    each edge.
!
!    Output, integer ( kind = 4 ) EDGE_LABEL(EDGES), a label for each edge.
!
!    Output, integer ( kind = 4 ) TRIANGLE_VERTEX(3,TRIANGLES), the vertices
!    that form each triangle.
!
!    Output, integer ( kind = 4 ) TRIANGLE_LABEL(TRIANGLES), a label for each
!    triangle.
!
!    Output, integer ( kind = 4 ) QUADRILATERAL_VERTEX(4,QUADRILATERALS), the
!    vertices that form each quadrilateral.
!
!    Output, integer ( kind = 4 ) QUADRILATERAL_LABEL(QUADRILATERALS), a label
!    for each quadrilateral.
!
!    Output, integer ( kind = 4 ) TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
!    vertices that form each tetrahedron.
!
!    Output, integer ( kind = 4 ) TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
!    each tetrahedron.
!
!    Output, integer ( kind = 4 ) HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the vertices
!    that form each hexahedron.
!
!    Output, integer ( kind = 4 ) HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
!    each hexahedron.
!
  implicit none

  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_label(edges)
  integer ( kind = 4 ) edge_vertex(2,edges)
  character ( len = * ) filename
  integer ( kind = 4 ) fileunit
  integer ( kind = 4 ) hexahedron
  integer ( kind = 4 ) hexahedron_label(hexahedrons)
  integer ( kind = 4 ) hexahedron_vertex(8,hexahedrons)
  integer ( kind = 4 ) i4vec(9)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  character ( len = 80 ) keyword
  integer ( kind = 4 ) length
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) quadrilateral
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals)
  real    ( kind = 8 ) r8vec(9)
  logical s_begin
  logical s_eqi
  integer ( kind = 4 ) tetrahedron
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  character ( len = 255 ) text
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle_label(triangles)
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  integer ( kind = 4 ) vertex
  real    ( kind = 8 ) vertex_coordinate(dim,vertices)
  integer ( kind = 4 ) vertex_label(vertices)
!
!  Initialize everything to nothing.
!
  vertex_coordinate(1:dim,1:vertices) = 0.0D+00
  vertex_label(1:vertices) = 0
  edge_vertex(1:2,1:edges) = 0
  edge_label(1:edges) = 0
  triangle_vertex(1:3,1:triangles) = 0
  triangle_label(1:triangles) = 0
  quadrilateral_vertex(1:4,1:quadrilaterals) = 0
  quadrilateral_label(1:quadrilaterals) = 0
  tetrahedron_vertex(1:4,1:tetrahedrons) = 0
  tetrahedron_label(1:tetrahedrons) = 0
  hexahedron_vertex(1:8,1:hexahedrons) = 0
  hexahedron_label(1:hexahedrons) = 0
!
!  Open the file.
!
  call get_unit ( fileunit )

  open ( unit = fileunit, file = filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open file.'
    stop
  end if
!
!  Read lines til you get alphanumerics and determine a "mode"
!
  line_num = 0
  keyword = 'NONE'

  do

    read ( fileunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

    if ( len_trim ( text ) == 0 ) then
      keyword = 'NONE'
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if
!
!  Remove initial blanks.
!
    text = adjustl ( text )
!
!  Expecting a keyword.
!
        if ( s_eqi ( text, 'CORNERS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'DIMENSION' ) ) then

      keyword = 'DIMENSION'

    else if ( s_eqi ( text, 'EDGES' ) ) then

      keyword = 'EDGES'

    else if ( s_eqi ( text, 'END' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  END statement encountered.'
      exit

    else if ( s_eqi ( text, 'HEXAHEDRA' ) .or. &
              s_eqi ( text, 'HEXAHEDRONS' ) ) then

      keyword = 'HEXAHEDRONS'

    else if ( s_begin ( text, 'MESHVERSIONFORMATTED' ) ) then

    else if ( s_eqi ( text, 'NORMALATQUADRILATERALVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALATTRIANGLEVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALATVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'QUADRILATERALS' ) ) then

      keyword = 'QUADRILATERALS'

    else if ( s_eqi ( text, 'REQUIREDEDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'REQUIREDVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'RIDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TANGENTATEDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TANGENTS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TETRAHEDRA' ) .or. &
              s_eqi ( text, 'TETRAHEDRONS' ) ) then

      keyword = 'TETRAHEDRONS'

    else if ( s_eqi ( text, 'TRIANGLES' ) ) then

      keyword = 'TRIANGLES'

    else if ( s_eqi ( text, 'VERTICES' ) ) then

      keyword = 'VERTICES'
!
!  Presumably, numeric data to be processed by keyword.
!
    else if ( s_eqi ( keyword, 'DIMENSION' ) ) then

      call s_to_i4 ( text, dim, ierror, length )

      keyword = 'NONE'

    else if ( s_eqi ( keyword, 'EDGES' ) ) then

      call s_to_i4 ( text, edges, ierror, length )

      keyword = 'EDGE_VERTEX'
      edge = 0

    else if ( s_eqi ( keyword, 'EDGE_VERTEX' ) ) then

      call s_to_i4vec ( text, 3, i4vec, ierror )
      edge = edge + 1
      edge_vertex(1:2,edge) = i4vec(1:2)
      edge_label(edge) = i4vec(3)

    else if ( s_eqi ( keyword, 'HEXAHEDRONS' ) ) then

      call s_to_i4 ( text, hexahedrons, ierror, length )

      keyword = 'HEXAHEDRON_VERTEX'
      hexahedron = 0

    else if ( s_eqi ( keyword, 'HEXAHEDRON_VERTEX' ) ) then

      call s_to_i4vec ( text, 9, i4vec, ierror )
      hexahedron = hexahedron + 1
      hexahedron_vertex(1:8,hexahedron) = i4vec(1:8)
      hexahedron_label(hexahedron) = i4vec(9)

    else if ( s_eqi ( keyword, 'QUADRILATERALS' ) ) then

      call s_to_i4 ( text, quadrilaterals, ierror, length )

      keyword = 'QUADRILATERAL_VERTEX'
      quadrilateral = 0

    else if ( s_eqi ( keyword, 'QUADRILATERAL_VERTEX' ) ) then

      call s_to_i4vec ( text, 5, i4vec, ierror )
      quadrilateral = quadrilateral + 1
      quadrilateral_vertex(1:4,quadrilateral) = i4vec(1:4)
      quadrilateral_label(quadrilateral) = i4vec(5)

    else if ( s_eqi ( keyword, 'TETRAHEDRONS' ) ) then

      call s_to_i4 ( text, tetrahedrons, ierror, length )

      keyword = 'TETRAHEDRON_VERTEX'
      tetrahedron = 0

    else if ( s_eqi ( keyword, 'TETRAHEDRON_VERTEX' ) ) then

      call s_to_i4vec ( text, 5, i4vec, ierror )
      tetrahedron = tetrahedron + 1
      tetrahedron_vertex(1:4,tetrahedron) = i4vec(1:4)
      tetrahedron_label(tetrahedron) = i4vec(5)

    else if ( s_eqi ( keyword, 'TRIANGLES' ) ) then

      call s_to_i4 ( text, triangles, ierror, length )

      keyword = 'TRIANGLE_VERTEX'
      triangle = 0

    else if ( s_eqi ( keyword, 'TRIANGLE_VERTEX' ) ) then

      call s_to_i4vec ( text, 4, i4vec, ierror )
      triangle = triangle + 1
      triangle_vertex(1:3,triangle) = i4vec(1:3)
      triangle_label(triangle) = i4vec(4)

    else if ( s_eqi ( keyword, 'VERTICES' ) ) then

      call s_to_i4 ( text, vertices, ierror, length )

      keyword = 'VERTEX_COORDINATE'
      vertex = 0

    else if ( s_eqi ( keyword, 'VERTEX_COORDINATE' ) ) then

      call s_to_r8vec ( text, dim + 1, r8vec, ierror )
      vertex = vertex + 1
      vertex_coordinate(1:dim,vertex) = r8vec(1:dim)
      vertex_label(vertex) = int ( r8vec(dim+1) )

    else if ( s_eqi ( keyword, 'SKIP' ) ) then

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MESH_DATA_READ - Fatal error!'
      write ( *, '(a,i8)' ) &
        '  Could not find keyword while reading line ', line_num
      write ( *, '(a)' ) '"' // trim ( text ) // '".'
      stop

    end if

  end do
!
!  Close the file.
!
  close ( unit = fileunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Read ', line_num, &
    ' lines from "' // trim ( filename ) // '".'

  return
end
subroutine mesh_size_print ( dim, vertices, edges, triangles, quadrilaterals, &
  tetrahedrons, hexahedrons )

!*****************************************************************************80
!
!! MESH_SIZE_PRINT prints mesh sizes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pascal Frey,
!    MEDIT: An interactive mesh visualization software,
!    Technical Report RT-0253,
!    Institut National de Recherche en Informatique et en Automatique,
!    03 December 2001.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 2 or 3.
!
!    Input, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Input, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Input, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of dimensions = ', dim
  write ( *, '(a,i8)' ) '  Number of vertices = ', vertices
  write ( *, '(a,i8)' ) '  Number of edges = ', edges
  write ( *, '(a,i8)' ) '  Number of triangles = ', triangles
  write ( *, '(a,i8)' ) '  Number of quadrilaterals = ', quadrilaterals
  write ( *, '(a,i8)' ) '  Number of tetrahedrons = ', tetrahedrons
  write ( *, '(a,i8)' ) '  Number of hexahedrons = ', hexahedrons

  return
end
subroutine mesh_size_read ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons )

!*****************************************************************************80
!
!! MESH_SIZE_READ reads sizes from a MESH file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pascal Frey,
!    MEDIT: An interactive mesh visualization software,
!    Technical Report RT-0253,
!    Institut National de Recherche en Informatique et en Automatique,
!    03 December 2001.
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the MESH file.
!
!    Output, integer ( kind = 4 ) DIM, the spatial dimension, which should be 2 or 3.
!
!    Output, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Output, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Output, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Output, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Output, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Output, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  character ( len = * ) filename
  integer ( kind = 4 ) fileunit
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 80 ) keyword
  integer ( kind = 4 ) length
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) quadrilaterals
  logical s_begin
  logical s_eqi
  integer ( kind = 4 ) tetrahedrons
  character ( len = 255 ) text
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices
!
!  Initialize everything to nothing.
!
  dim = 0
  vertices = 0
  edges = 0
  triangles = 0
  quadrilaterals = 0
  tetrahedrons = 0
  hexahedrons = 0
!
!  Open the file.
!
  call get_unit ( fileunit )

  open ( unit = fileunit, file = filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_SIZE_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open file.'
    stop
  end if
!
!  Read lines til you get alphanumerics and determine a "mode"
!
  line_num = 0
  keyword = 'NONE'

  do

    read ( fileunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

    if ( len_trim ( text ) == 0 ) then
      keyword = 'NONE'
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if
!
!  Remove initial blanks.
!
    text = adjustl ( text )
!
!  Expecting a keyword.
!
        if ( s_eqi ( text, 'CORNERS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'DIMENSION' ) ) then

      keyword = 'DIMENSION'

    else if ( s_eqi ( text, 'EDGES' ) ) then

      keyword = 'EDGES'

    else if ( s_eqi ( text, 'END' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  END statement encountered.'
      exit

    else if ( s_eqi ( text, 'HEXAHEDRA' ) .or. &
              s_eqi ( text, 'HEXAHEDRONS' ) ) then

      keyword = 'HEXAHEDRONS'

    else if ( s_begin ( text, 'MESHVERSIONFORMATTED' ) ) then

    else if ( s_eqi ( text, 'NORMALATQUADRILATERALVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALATTRIANGLEVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALATVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'QUADRILATERALS' ) ) then

      keyword = 'QUADRILATERALS'

    else if ( s_eqi ( text, 'REQUIREDEDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'REQUIREDVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'RIDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TANGENTATEDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TANGENTS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TETRAHEDRA' ) .or. &
              s_eqi ( text, 'TETRAHEDRONS' ) ) then

      keyword = 'TETRAHEDRONS'

    else if ( s_eqi ( text, 'TRIANGLES' ) ) then

      keyword = 'TRIANGLES'

    else if ( s_eqi ( text, 'VERTICES' ) ) then

      keyword = 'VERTICES'
!
!  Presumably, numeric data to be processed by keyword.
!
    else if ( s_eqi ( keyword, 'DIMENSION' ) ) then

      call s_to_i4 ( text, dim, ierror, length )

      keyword = 'NONE'

    else if ( s_eqi ( keyword, 'EDGES' ) ) then

      call s_to_i4 ( text, edges, ierror, length )

      keyword = 'EDGE_VERTEX'

    else if ( s_eqi ( keyword, 'EDGE_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'HEXAHEDRONS' ) ) then

      call s_to_i4 ( text, hexahedrons, ierror, length )

      keyword = 'HEXAHEDRON_VERTEX'

    else if ( s_eqi ( keyword, 'HEXAHEDRON_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'QUADRILATERALS' ) ) then

      call s_to_i4 ( text, quadrilaterals, ierror, length )

      keyword = 'QUADRILATERAL_VERTEX'

    else if ( s_eqi ( keyword, 'QUADRILATERAL_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'TETRAHEDRONS' ) ) then

      call s_to_i4 ( text, tetrahedrons, ierror, length )

      keyword = 'TETRAHEDRON_VERTEX'

    else if ( s_eqi ( keyword, 'TETRAHEDRON_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'TRIANGLES' ) ) then

      call s_to_i4 ( text, triangles, ierror, length )

      keyword = 'TRIANGLE_VERTEX'

    else if ( s_eqi ( keyword, 'TRIANGLE_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'VERTICES' ) ) then

      call s_to_i4 ( text, vertices, ierror, length )

      keyword = 'VERTEX_COORDINATE'

    else if ( s_eqi ( keyword, 'VERTEX_COORDINATE' ) ) then

    else if ( s_eqi ( keyword, 'SKIP' ) ) then

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MESH_SIZE_READ - Fatal error!'
      write ( *, '(a,i8)' ) &
        '  Could not find keyword while reading line ', line_num
      write ( *, '(a)' ) '"' // trim ( text ) // '".'
      stop

    end if

  end do
!
!  Close the file.
!
  close ( unit = fileunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Read ', line_num, &
    ' lines from "' // trim ( filename ) // '".'

  return
end
subroutine mesh_write ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
  vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
  quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
  tetrahedron_label, hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! MESH_WRITE writes sizes and data to a MESH file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pascal Frey,
!    MEDIT: An interactive mesh visualization software,
!    Technical Report RT-0253,
!    Institut National de Recherche en Informatique et en Automatique,
!    03 December 2001.
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file to be created.
!    Ordinarily, the name should include the extension ".mesh".
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should
!    be 2 or 3.
!
!    Input, integer ( kind = 4 ) VERTICES, the number of vertices.
!
!    Input, real ( kind = 8 ) VERTEX_COORDINATE(DIM,VERTICES), the coordinates
!    of each vertex.
!
!    Input, integer ( kind = 4 ) VERTEX_LABEL(VERTICES), a label for
!    each vertex.
!
!    Input, integer ( kind = 4 ) EDGES, the number of edges (may be 0).
!
!    Input, integer ( kind = 4 ) EDGE_VERTEX(2,EDGES), the vertices that form
!    each edge.
!
!    Input, integer ( kind = 4 ) EDGE_LABEL(EDGES), a label for each edge.
!
!    Input, integer ( kind = 4 ) TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer ( kind = 4 ) TRIANGLE_VERTEX(3,TRIANGLES), the vertices
!    that form each triangle.
!
!    Input, integer ( kind = 4 ) TRIANGLE_LABEL(TRIANGLES), a label for each
!    triangle.
!
!    Input, integer ( kind = 4 ) QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer ( kind = 4 ) QUADRILATERAL_VERTEX(4,QUADRILATERALS), the
!    vertices that form each quadrilateral.
!
!    Input, integer ( kind = 4 ) QUADRILATERAL_LABEL(QUADRILATERALS), a label
!    for each quadrilateral.
!
!    Input, integer ( kind = 4 ) TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
!    vertices that form each tetrahedron.
!
!    Input, integer ( kind = 4 ) TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
!    each tetrahedron.
!
!    Input, integer ( kind = 4 ) HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
!    Input, integer ( kind = 4 ) HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the vertices
!    that form each hexahedron.
!
!    Input, integer ( kind = 4 ) HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
!    each hexahedron.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  integer ( kind = 4 ) edge_label(edges)
  integer ( kind = 4 ) edge_vertex(2,edges)
  character ( len = * ) filename
  integer ( kind = 4 ) fileunit
  integer ( kind = 4 ) hexahedron_label(hexahedrons)
  integer ( kind = 4 ) hexahedron_vertex(8,hexahedrons)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals)
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  integer ( kind = 4 ) triangle_label(triangles)
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  real ( kind = 8 ) vertex_coordinate(dim,vertices)
  integer ( kind = 4 ) vertex_label(vertices)
!
!  Open the file.
!
  call get_unit ( fileunit )

  open ( unit = fileunit, file = filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open file.'
    stop
  end if

  write ( fileunit, '(a)' ) 'MeshVersionFormatted 1'
  write ( fileunit, '(a)' ) '#  Created by mesh_write.f90'
!
!  Dimension information.
!
  write ( fileunit, '(a)' ) ' '
  write ( fileunit, '(a)' ) 'Dimension'
  write ( fileunit, '(i8)' ) dim
!
!  Vertices.
!
  write ( fileunit, '(a)' ) ' '
  write ( fileunit, '(a)' ) 'Vertices'
  write ( fileunit, '(i8)' ) vertices
  if ( dim == 2 ) then
    do j = 1, vertices
      write ( fileunit, '(2(2x,f10.6),2x,i8)' ) &
        vertex_coordinate(1:dim,j), vertex_label(j)
    end do
  else if ( dim == 3 ) then
    do j = 1, vertices
      write ( fileunit, '(3(2x,f10.6),2x,i8)' ) &
        vertex_coordinate(1:dim,j), vertex_label(j)
    end do
  end if
!
!  Edges.
!
  if ( 0 < edges ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Edges'
    write ( fileunit, '(i8)' ) edges
    do j = 1, edges
      write ( fileunit, '(2(2x,i8),2x,i8)' ) &
        edge_vertex(1:2,j), edge_label(j)
    end do
  end if
!
!  Triangles.
!
  if ( 0 < triangles ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Triangles'
    write ( fileunit, '(i8)' ) triangles
    do j = 1, triangles
      write ( fileunit, '(3(2x,i8),2x,i8)' ) &
        triangle_vertex(1:3,j), triangle_label(j)
    end do
  end if
!
!  Quadrilaterals.
!
  if ( 0 < quadrilaterals ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Quadrilaterals'
    write ( fileunit, '(i8)' ) quadrilaterals
    do j = 1, quadrilaterals
      write ( fileunit, '(4(2x,i8),2x,i8)' ) &
        quadrilateral_vertex(1:4,j), quadrilateral_label(j)
    end do
  end if
!
!  Tetrahedron.
!
  if ( 0 < tetrahedrons ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Tetrahedra'
    write ( fileunit, '(i8)' ) tetrahedrons
    do j = 1, tetrahedrons
      write ( fileunit, '(4(2x,i8),2x,i8)' ) &
        tetrahedron_vertex(1:4,j), tetrahedron_label(j)
    end do
  end if
!
!  Hexahedron.
!
  if ( 0 < hexahedrons ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Hexahedra'
    write ( fileunit, '(i8)' ) hexahedrons
    do j = 1, hexahedrons
      write ( fileunit, '(8(2x,i8),2x,i8)' ) &
        hexahedron_vertex(1:8,j), hexahedron_label(j)
    end do
  end if
!
!  End.
!
  write ( fileunit, '(a)' ) ' '
  write ( fileunit, '(a)' ) 'End'

  close ( unit = fileunit )

  return
end
subroutine r8mat_copy ( m, n, a, b )

!*****************************************************************************80
!
!! R8MAT_COPY copies an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix to be copied.
!
!    Output, real ( kind = 8 ) B(M,N), a copy of the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) b(m,n)

  b(1:m,1:n) = a(1:m,1:n)

  return
end
subroutine r8vec_copy ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_COPY copies an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), the vector to be copied.
!
!    Output, real ( kind = 8 ) A2(N), a copy of A1.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)

  a2(1:n) = a1(1:n)

  return
end
function s_begin ( s1, s2 )

!*****************************************************************************80
!
!! S_BEGIN is TRUE if one string matches the beginning of the other.
!
!  Discussion:
!
!    The strings are compared, ignoring blanks, spaces and capitalization.
!
!  Example:
!
!     S1              S2      S_BEGIN
!
!    'Bob'          'BOB'     TRUE
!    '  B  o b '    ' bo b'   TRUE
!    'Bob'          'Bobby'   TRUE
!    'Bobo'         'Bobb'    FALSE
!    ' '            'Bob'     FALSE    (Do not allow a blank to match
!                                       anything but another blank string.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to be compared.
!
!    Output, logical S_BEGIN, is TRUE if the strings match up to
!    the end of the shorter string, ignoring case.
!
  implicit none

  logical                ch_eqi
  integer   ( kind = 4 ) i1
  integer   ( kind = 4 ) i2
  logical                s_begin
  character ( len = * )  s1
  integer   ( kind = 4 ) s1_length
  character ( len = * )  s2
  integer   ( kind = 4 ) s2_length

  s1_length = len_trim ( s1 )
  s2_length = len_trim ( s2 )
!
!  If either string is blank, then both must be blank to match.
!  Otherwise, a blank string matches anything, which is not
!  what most people want.
!
  if ( s1_length == 0 .or. s2_length == 0 ) then

    if ( s1_length == 0 .and. s2_length == 0 ) then
      s_begin = .true.
    else
      s_begin = .false.
    end if

    return

  end if

  i1 = 0
  i2 = 0
!
!  Find the next nonblank in S1.
!
  do

    do

      i1 = i1 + 1

      if ( s1_length < i1 ) then
        s_begin = .true.
        return
      end if

      if ( s1(i1:i1) /= ' ' ) then
        exit
      end if

    end do
!
!  Find the next nonblank in S2.
!
    do

      i2 = i2 + 1

      if ( s2_length < i2 ) then
        s_begin = .true.
        return
      end if

      if ( s2(i2:i2) /= ' ' ) then
        exit
      end if

    end do
!
!  If the characters match, get the next pair.
!
    if ( .not. ch_eqi ( s1(i1:i1), s2(i2:i2) ) ) then
      exit
    end if

  end do

  s_begin = .false.

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Discussion:
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

  character              c1
  character              c2
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lenc
  logical                s_eqi
  character ( len = *  ) s1
  integer   ( kind = 4 ) s1_length
  character ( len = *  ) s2
  integer   ( kind = 4 ) s2_length

  s1_length = len ( s1 )
  s2_length = len ( s2 )
  lenc = min ( s1_length, s2_length )

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

  do i = lenc + 1, s1_length
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, s2_length
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_to_i4 ( s, value, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) VALUE, the integer value read from the string.
!    If the string is blank, then VALUE will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters
!    of S used to make the integer.
!
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) length
  character ( len = * )  s
  integer   ( kind = 4 ) state
  character              :: TAB = achar ( 9 )
  integer   ( kind = 4 ) value

  value = 0
  ierror = 0
  length = 0

  state = 0
  isgn = 1

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  STATE = 0, haven't read anything.
!
    if ( state == 0 ) then

      if ( c == ' ' .or. c == TAB ) then

      else if ( c == '-' ) then
        state = 1
        isgn = -1
      else if ( c == '+' ) then
        state = 1
        isgn = +1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' .or. c == TAB ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 2, have read at least one digit, expecting more.
!
    else if ( state == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + iachar ( c ) - iachar ( '0' )

      else

        value = isgn * value
        ierror = 0
        length = i - 1
        return

      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( state == 2 ) then

    value = isgn * value
    ierror = 0
    length = len_trim ( s )

  else

    value = 0
    ierror = 1
    length = 0

  end if

  return
end
subroutine s_to_i4vec ( s, n, i4vec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an integer vector from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2003
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
!    Output, integer ( kind = 4 ) I4VEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) i4vec(n)
  integer   ( kind = 4 ) length
  character ( len = * )  s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_i4 ( s(ilo:), i4vec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

  end do

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 value from a string.
!
!  Discussion:
!
!    An "R8" value is simply a real number to be stored as a
!    variable of type "real ( kind = 8 )".
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
!    12 January 2009
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

  character              c
  logical                ch_eqi
  real      ( kind = 8 ) dval
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ihave
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) iterm
  integer   ( kind = 4 ) jbot
  integer   ( kind = 4 ) jsgn
  integer   ( kind = 4 ) jtop
  integer   ( kind = 4 ) length
  integer   ( kind = 4 ) ndig
  real      ( kind = 8 ) rbot
  real      ( kind = 8 ) rexp
  real      ( kind = 8 ) rtop
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  character           :: TAB = achar ( 9 )

  s_length = len_trim ( s )

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

    if ( s_length < length + 1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' .or. c == TAB ) then

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
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
  if ( iterm /= 1 .and. length + 1 == s_length ) then
    length = s_length
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
subroutine s_to_r8vec ( s, n, r8vec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Discussion:
!
!    An R8VEC is a vector of real values, of type "real ( kind = 8 )".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2005
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
!    Output, real ( kind = 8 ) R8VEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) lchar
  real      ( kind = 8 ) r8vec(n)
  character ( len = * )  s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), r8vec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

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
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
