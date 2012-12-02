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
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 3.
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
!    Output, double VERTEX_COORDINATE(3,VERTICES), the XYZ coordinates
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
  integer ( kind = 4 ), save :: triangle_label_save(triangles_save) = (/ &
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
  integer ( kind = 4 ), save :: triangle_vertex_save(3,triangles_save) &
    = reshape ( (/ &
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
  real ( kind = 8 ), save :: vertex_coordinate_save(3,VERTICES_SAVE) &
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

  call r8mat_copy ( 3, vertices, vertex_coordinate_save, vertex_coordinate )

  call i4vec_copy ( vertices, vertex_label_save, vertex_label )

  call i4mat_copy ( 3, triangles, triangle_vertex_save, triangle_vertex )

  call i4vec_copy ( triangles, triangle_label_save, triangle_label )

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
!    Output, integer DIM, the spatial dimension, which should be 3.
!
!    Output, integer VERTICES, the number of vertices.
!
!    Output, integer EDGES, the number of edges (may be 0).
!
!    Output, integer TRIANGLES, the number of triangles (may be 0).
!
!    Output, integer QUADRILATERALS, the number of quadrilaterals (may be 0).
!
!    Output, integer TETRAHEDRONS, the number of tetrahedrons (may be 0).
!
!    Output, integer HEXAHEDRONS, the number of hexahedrons (may be 0).
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
subroutine data_print ( dim, vertices, edges, triangles, quadrilaterals, &
  tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex, &
  edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, &
  quadrilateral_label, tetrahedron_vertex, tetrahedron_label, &
  hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! DATA_PRINT prints the data of an ICE grid dataset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 November 2010
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
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 3.
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
!    Input, double VERTEX_COORDINATE(3,VERTICES), the XYZ coordinates
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
!    Input, integer ( kind = 4 ) TRIANGLE_LABEL(TRIANGLES), a label for 
!    each triangle.
!
!    Input, integer ( kind = 4 ) QUADRILATERAL_VERTEX(4,QUADRILATERALS), 
!    the vertices that form each quadrilateral.
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
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals) 
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  integer ( kind = 4 ) triangle_label(triangles) 
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  real ( kind = 8 ) vertex_coordinate(3,vertices)
  integer ( kind = 4 ) vertex_label(vertices)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertices:'
  write ( *, '(a)' ) ' '
  do j = 1, vertices
    write ( *, '(3(2x,f10.4),2x,''('',i4,'')'')' ) &
      vertex_coordinate(1:3,j), vertex_label(j)
  end do

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
subroutine data_read ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
  vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
  quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
  tetrahedron_label, hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! DATA_READ reads ICE data from a NETCDF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 November 2010
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
!    Ordinarily, the name should include the extension '.nc'.
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 3.
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
!    Output, real VERTEX_COORDINATE(3,VERTICES), the XYZ coordinates
!    of each vertex.
!
!    Output, integer ( kind = 4 ) VERTEX_LABEL(VERTICES), a label for 
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
!    Input, integer ( kind = 4 ) TRIANGLE_LABEL(TRIANGLES), a label for 
!    each triangle.
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
  use netcdf

  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  integer ( kind = 4 ) dim_dimension
  integer ( kind = 4 ) dim_edges
  integer ( kind = 4 ) dim_eight
  integer ( kind = 4 ) dim_four
  integer ( kind = 4 ) dim_hexahedrons
  integer ( kind = 4 ) dim_quadrilaterals
  integer ( kind = 4 ) dim_tetrahedrons
  integer ( kind = 4 ) dim_three
  integer ( kind = 4 ) dim_triangles
  integer ( kind = 4 ) dim_two
  integer ( kind = 4 ) dim_vertices
  integer ( kind = 4 ) dimid
  integer ( kind = 4 ) dimids(2)
  integer ( kind = 4 ) edge_label(edges)
  integer ( kind = 4 ) edge_vertex(2,edges) 
  character ( len = * ) filename
  integer ( kind = 4 ) hexahedron_label(hexahedrons)
  integer ( kind = 4 ) hexahedron_vertex(8,hexahedrons)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) ncid
  integer ( kind = 4 ) ndims
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals)
  integer ( kind = 4 ) status
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  integer ( kind = 4 ) triangle_label(triangles) 
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  integer ( kind = 4 ) var_edge_label
  integer ( kind = 4 ) var_edge_vertex
  integer ( kind = 4 ) var_hexahedron_label
  integer ( kind = 4 ) var_hexahedron_vertex
  integer ( kind = 4 ) var_quadrilateral_label
  integer ( kind = 4 ) var_quadrilateral_vertex
  integer ( kind = 4 ) var_tetrahedron_label
  integer ( kind = 4 ) var_tetrahedron_vertex
  integer ( kind = 4 ) var_triangle_label
  integer ( kind = 4 ) var_triangle_vertex
  integer ( kind = 4 ) var_vertex_coordinate
  integer ( kind = 4 ) var_vertex_label
  integer ( kind = 4 ) varid
  real ( kind = 8 ) vertex_coordinate(3,vertices)
  integer ( kind = 4 ) vertex_label(vertices)
  integer ( kind = 4 ) xtype
!
!  Open the file.
!
  mode = NF90_NOCLOBBER
  status = nf90_open ( filename, mode, ncid )
!
!  Vertices.
!
  status = nf90_inq_varid ( ncid, 'Vertex_Coordinate', varid )
  status = nf90_get_var ( ncid, varid, vertex_coordinate )

  status = nf90_inq_varid ( ncid, 'Vertex_Label', varid )
  status = nf90_get_var ( ncid, varid, vertex_label )
!
!  Edges.
!
  if ( 0 < edges ) then
    status = nf90_inq_varid ( ncid, 'Edge_Vertex', varid )
    status = nf90_get_var ( ncid, varid, edge_vertex )

    status = nf90_inq_varid ( ncid, 'Edge_Label', varid )
    status = nf90_get_var ( ncid, varid, edge_label )
  end if
!
!  Triangles.
!
  if ( 0 < triangles ) then
    status = nf90_inq_varid ( ncid, 'Triangle_Vertex', varid )
    status = nf90_get_var ( ncid, varid, triangle_vertex )

    status = nf90_inq_varid ( ncid, 'Triangle_Label', varid )
    status = nf90_get_var ( ncid, varid, triangle_label )
  end if
!
!  Quadrilaterals.
!
  if ( 0 < quadrilaterals ) then
    status = nf90_inq_varid ( ncid, 'Quadrilateral_Vertex', varid )
    status = nf90_get_var ( ncid, varid, quadrilateral_vertex )

    status = nf90_inq_varid ( ncid, 'Quadrilateral_Label', varid )
    status = nf90_get_var ( ncid, varid, quadrilateral_label )
  end if
!
!  Tetrahedrons.
!
  if ( 0 < tetrahedrons ) then
    status = nf90_inq_varid ( ncid, 'Tetrahedron_Vertex', varid )
    status = nf90_get_var ( ncid, varid, tetrahedron_vertex )

    status = nf90_inq_varid ( ncid, 'Tetrahedron_Label', varid )
    status = nf90_get_var ( ncid, varid, tetrahedron_label )
  end if
!
!  Hexahedrons.
!
  if ( 0 < hexahedrons ) then
    status = nf90_inq_varid ( ncid, 'Hexahedron_Vertex', varid )
    status = nf90_get_var ( ncid, varid, hexahedron_vertex )

    status = nf90_inq_varid ( ncid, 'Hexahedron_Label', varid )
    status = nf90_get_var ( ncid, varid, hexahedron_label )
  end if
!
!  Close the file.
!
  status = nf90_close ( ncid )

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
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 3.
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
!    Output, double VERTEX_COORDINATE(3,VERTICES), the XYZ coordinates
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

  call r8vec_copy ( 3 * vertices, vertex_coordinate_save, vertex_coordinate )
  call i4vec_copy ( vertices, vertex_label_save, vertex_label )
  call i4vec_copy ( 4 * quadrilaterals, quadrilateral_vertex_save, &
    quadrilateral_vertex )
  call i4vec_copy ( quadrilaterals, quadrilateral_label_save, &
    quadrilateral_label )
  call i4vec_copy ( 8 * hexahedrons, hexahedron_vertex_save, hexahedron_vertex )
  call i4vec_copy ( hexahedrons, hexahedron_label_save, hexahedron_label )

  return
end
subroutine hexahexa_2x2x2_size ( dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons )

!*****************************************************************************80
!
!! HEXAHEXA_2x2x2_SIZE defines the sizes for a 3D hexahedral mesh.
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
!    Output, integer DIM, the spatial dimension, which should be 3.
!
!    Output, integer VERTICES, the number of vertices.
!
!    Output, integer EDGES, the number of edges (may be 0).
!
!    Output, integer TRIANGLES, the number of triangles (may be 0).
!
!    Output, integer QUADRILATERALS, the number of quadrilaterals (may be 0).
!
!    Output, integer TETRAHEDRONS, the number of tetrahedrons (may be 0).
!
!    Output, integer HEXAHEDRONS, the number of hexahedrons (may be 0).
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
!    An I4MAT is an array of I4's.
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
!    Input, integer ( kind = 4 ) A1(M,N), the vector to be copied.
!
!    Output, integer ( kind = 4 ) A2(M,N), a copy of A1.
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
subroutine ice_write ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, &
  edge_vertex, edge_label, triangle_vertex, triangle_label, &
  quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
  tetrahedron_label, hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! ICE_WRITE writes 3D ICE sizes and data to a NETCDF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 November 2010
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
!    Ordinarily, the name should include the extension '.nc'.
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 3.
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
!    Input, real VERTEX_COORDINATE(3,VERTICES), the XYZ coordinates
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
!    Input, integer ( kind = 4 ) TRIANGLE_LABEL(TRIANGLES), a label for 
!    each triangle.
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
  use netcdf

  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) edges
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) tetrahedrons
  integer ( kind = 4 ) triangles
  integer ( kind = 4 ) vertices

  integer ( kind = 4 ) dim_dimension
  integer ( kind = 4 ) dim_edges
  integer ( kind = 4 ) dim_eight
  integer ( kind = 4 ) dim_four
  integer ( kind = 4 ) dim_hexahedrons
  integer ( kind = 4 ) dim_quadrilaterals
  integer ( kind = 4 ) dim_tetrahedrons
  integer ( kind = 4 ) dim_three
  integer ( kind = 4 ) dim_triangles
  integer ( kind = 4 ) dim_two
  integer ( kind = 4 ) dim_vertices
  integer ( kind = 4 ) dimids(2)
  integer ( kind = 4 ) edge_label(edges)
  integer ( kind = 4 ) edge_vertex(2,edges) 
  character ( len = * ) filename
  integer ( kind = 4 ) hexahedron_label(hexahedrons)
  integer ( kind = 4 ) hexahedron_vertex(8,hexahedrons)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) ncid
  integer ( kind = 4 ) ndims
  integer ( kind = 4 ) quadrilateral_label(quadrilaterals)
  integer ( kind = 4 ) quadrilateral_vertex(4,quadrilaterals)
  integer ( kind = 4 ) status
  integer ( kind = 4 ) tetrahedron_label(tetrahedrons)
  integer ( kind = 4 ) tetrahedron_vertex(4,tetrahedrons)
  integer ( kind = 4 ) triangle_label(triangles) 
  integer ( kind = 4 ) triangle_vertex(3,triangles)
  integer ( kind = 4 ) var_edge_label
  integer ( kind = 4 ) var_edge_vertex
  integer ( kind = 4 ) var_hexahedron_label
  integer ( kind = 4 ) var_hexahedron_vertex
  integer ( kind = 4 ) var_quadrilateral_label
  integer ( kind = 4 ) var_quadrilateral_vertex
  integer ( kind = 4 ) var_tetrahedron_label
  integer ( kind = 4 ) var_tetrahedron_vertex
  integer ( kind = 4 ) var_triangle_label
  integer ( kind = 4 ) var_triangle_vertex
  integer ( kind = 4 ) var_vertex_coordinate
  integer ( kind = 4 ) var_vertex_label
  real ( kind = 8 ) vertex_coordinate(3,vertices)
  integer ( kind = 4 ) vertex_label(vertices)
  integer ( kind = 4 ) xtype
!
!  Create the file.  This automatically 'opens' it as well.
!
  mode = NF90_CLOBBER
  status = nf90_create ( filename, mode, ncid )
!
!  Put NETCDF into 'define' mode.
!
  status = nf90_redef ( ncid )
!
!  Dimension information.
!
!  If a dimension has length 0, it seems to be taken to be the unlimited
!  dimension (not what you want) and then if you have two such dimensions,
!  you get a complaint that you have tried to define the unlimited dimension
!  twice.  The fix requires the programmer not to write anything whose
!  dimension is zero.
!
  status = nf90_def_dim ( ncid, 'Dimension', dim, dim_dimension )

  status = nf90_def_dim ( ncid, 'Vertices', vertices, dim_vertices )

  if ( 0 < edges ) then
    status = nf90_def_dim ( ncid, 'Edges', edges, dim_edges )
  end if

  if ( 0 < triangles ) then
    status = nf90_def_dim ( ncid, 'Triangles', triangles, dim_triangles )
  end if

  if ( 0 < quadrilaterals ) then
    status = nf90_def_dim ( ncid, 'Quadrilaterals', quadrilaterals, &
      dim_quadrilaterals )
  end if

  if ( 0 < tetrahedrons ) then
    status = nf90_def_dim ( ncid, 'Tetrahedrons', tetrahedrons, &
      dim_tetrahedrons )
  end if

  if ( 0 < hexahedrons ) then
    status = nf90_def_dim ( ncid, 'Hexahedrons', hexahedrons, dim_hexahedrons )
  end if

  status = nf90_def_dim ( ncid, 'Two', 2, dim_two )
  status = nf90_def_dim ( ncid, 'Three', 3, dim_three )
  status = nf90_def_dim ( ncid, 'Four', 4, dim_four )
  status = nf90_def_dim ( ncid, 'Eight', 8, dim_eight )
!
!  Define variables.
!
  ndims = 2
  dimids(1) = dim_three
  dimids(2) = dim_vertices
  status = nf90_def_var ( ncid, 'Vertex_Coordinate', NF90_DOUBLE, &
    dimids, var_vertex_coordinate )

  ndims = 1
  dimids(1) = dim_vertices
  status = nf90_def_var ( ncid, 'Vertex_Label', NF90_INT, dimids, &
    var_vertex_label )

  if ( 0 < edges ) then
    ndims = 2
    dimids(1) = dim_two
    dimids(2) = dim_edges
    status = nf90_def_var ( ncid, 'Edge_Vertex', NF90_INT, dimids, &
      var_edge_vertex )

    ndims = 1
    dimids(1) = dim_edges
    status = nf90_def_var ( ncid, 'Edge_Label', NF90_INT, dimids, &
      var_edge_label )
  end if

  if ( 0 < triangles ) then
    ndims = 2
    dimids(1) = dim_three
    dimids(2) = dim_triangles
    status = nf90_def_var ( ncid, 'Triangle_Vertex', NF90_INT, dimids, &
      var_triangle_vertex )

    ndims = 1
    dimids(1) = dim_triangles
    status = nf90_def_var ( ncid, 'Triangle_Label', NF90_INT, dimids, &
      var_triangle_label )
  end if

  if ( 0 < quadrilaterals ) then
    ndims = 2
    dimids(1) = dim_four
    dimids(2) = dim_quadrilaterals
    status = nf90_def_var ( ncid, 'Quadrilateral_Vertex', NF90_INT, &
      dimids, var_quadrilateral_vertex )

    ndims = 1
    dimids(1) = dim_quadrilaterals
    status = nf90_def_var ( ncid, 'Quadrilateral_Label', NF90_INT, &
      dimids, var_quadrilateral_label )
  end if

  if ( 0 < tetrahedrons ) then
    ndims = 2
    dimids(1) = dim_four
    dimids(2) = dim_tetrahedrons
    status = nf90_def_var ( ncid, 'Tetrahedron_Vertex', NF90_INT, &
      dimids, var_tetrahedron_vertex )

    ndims = 1
    dimids(1) = dim_tetrahedrons
    status = nf90_def_var ( ncid, 'Tetrahedron_Label', NF90_INT, &
      dimids, var_tetrahedron_label )
  end if

  if ( 0 < hexahedrons ) then
    ndims = 2
    dimids(1) = dim_eight
    dimids(2) = dim_hexahedrons
    status = nf90_def_var ( ncid, 'Hexahedron_Vertex', NF90_INT, &
      dimids, var_hexahedron_vertex )

    ndims = 1
    dimids(1) = dim_hexahedrons
    status = nf90_def_var ( ncid, 'Hexahedron_Label', NF90_INT, dimids, &
      var_hexahedron_label )
  end if
!
!  Terminate the definition phase.
!
  status = nf90_enddef ( ncid )
!
!  Write the data.
!
  status = nf90_put_var ( ncid, var_vertex_coordinate, vertex_coordinate )
  status = nf90_put_var ( ncid, var_vertex_label, vertex_label )

  if ( 0 < edges ) then
    status = nf90_put_var ( ncid, var_edge_vertex, edge_vertex )
    status = nf90_put_var ( ncid, var_edge_label, edge_label )
  end if

  if ( 0 < triangles ) then
    status = nf90_put_var ( ncid, var_triangle_vertex, triangle_vertex )
    status = nf90_put_var ( ncid, var_triangle_label, triangle_label )
  end if

  if ( 0 < quadrilaterals ) then
    status = nf90_put_var ( ncid, var_quadrilateral_vertex, &
      quadrilateral_vertex )
    status = nf90_put_var ( ncid, var_quadrilateral_label, quadrilateral_label )
  end if

  if ( 0 < tetrahedrons ) then
    status = nf90_put_var ( ncid, var_tetrahedron_vertex, tetrahedron_vertex )
    status = nf90_put_var ( ncid, var_tetrahedron_label, tetrahedron_label )
  end if

  if ( 0 < hexahedrons ) then
    status = nf90_put_var ( ncid, var_hexahedron_vertex, hexahedron_vertex )
    status = nf90_put_var ( ncid, var_hexahedron_label, hexahedron_label )
  end if
!
!  Close the file.
!
  status = nf90_close ( ncid )

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
subroutine size_print ( dim, vertices, edges, triangles, quadrilaterals, &
  tetrahedrons, hexahedrons )

!*****************************************************************************80
!
!! SIZE_PRINT prints the sizes of an ICE grid dataset.
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
!    Input, integer ( kind = 4 ) DIM, the spatial dimension, which should be 3.
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
subroutine size_read ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons )

!*****************************************************************************80
!
!! SIZE_READ reads ICE sizes from a NETCDF file.
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
!   John Burkardt
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
!    Input, string FILENAME, the name of the file to be read.
!    Ordinarily, the name should include the extension '.nc'.
!
!    Output, integer ( kind = 4 ) DIM, the spatial dimension, which should be 3.
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
  use netcdf

  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dimid
  integer ( kind = 4 ) edges
  character ( len = * ) filename
  integer ( kind = 4 ) hexahedrons
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) ncid
  integer ( kind = 4 ) quadrilaterals
  integer ( kind = 4 ) status
  integer ( kind = 4 ) tetrahedrons
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
  mode = NF90_NOWRITE
  status = nf90_open ( filename, mode, ncid )
!
!  Get the dimension information.
!
!  In an act of perplexing effect, the F90 NETCDF does not
!  include a "NF90_INQ_DIMLEN" function to allow you to determine
!  the length of a dimension from its ID.  Instead, NF90_INQUIRE_DIMENSION
!  returns both name and length.  Sadly, then, I must use the name
!  to get the ID, and then the ID to get the dimension...and the name again.
!  Oh, I see, this lets them use OPTIONAL arguments.  Whoopee!
!
  if ( nf90_inq_dimid ( ncid, 'Dimension', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = dim )
  end if

  if ( nf90_inq_dimid ( ncid, 'Vertices', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = vertices )
  end if

  if ( nf90_inq_dimid ( ncid, 'Edges', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = edges )
  end if

  if ( nf90_inq_dimid ( ncid, 'Triangles', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = triangles )
  end if

  if ( nf90_inq_dimid ( ncid, 'Quadrilaterals', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = quadrilaterals )
  end if

  if ( nf90_inq_dimid ( ncid, 'Tetrahedra', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = tetrahedrons )
  end if

  if ( nf90_inq_dimid ( ncid, 'Tetrahedrons', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = tetrahedrons )
  end if

  if ( nf90_inq_dimid ( ncid, 'Hexahedra', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = hexahedrons )
  end if

  if ( nf90_inq_dimid ( ncid, 'Hexahedrons', dimid ) == NF90_NOERR ) then
    status = nf90_inquire_dimension ( ncid, dimid, len = hexahedrons )
  end if
!
!  Close the file.
!
  status = nf90_close ( ncid )

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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
