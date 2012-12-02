# graph_mst_tree.plot  04 January 2011
#
file graph_mst_tree.ps

  space -0.5 -0.5 5.5 5.5

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 5.0 5.0 6 6
#
#  Red things
#
    fill_rgb 1.0 0.0 0.0
    line_width 8

    line   1.0  4.4  0.4  2.6
#   line   1.0  4.4  3.0  3.4
#   line   0.4  2.6  1.0  1.0
    line   0.4  2.6  2.0  2.0
    line   1.0  1.0  2.0  2.0
#   line   1.0  1.0  3.0  0.8
#   line   2.0  2.0  4.0  2.0
    line   2.0  2.0  3.0  0.8
    line   4.0  2.0  3.0  3.4
    line   4.0  2.0  4.6  4.0
    line   4.0  2.0  4.8  1.6
    line   4.0  2.0  3.0  0.8
#   line   3.0  3.4  4.6  4.0
#   line   4.6  4.0  4.8  1.6
#   line   4.8  1.6  3.0  0.8
    line   4.8  1.6  4.6  0.2
#   line   3.0  0.8  4.6  0.2

    fill_rgb 0.5 0.5 0.5
    line_width 2

#   line   1.0  4.4  0.4  2.6
    line   1.0  4.4  3.0  3.4
    line   0.4  2.6  1.0  1.0
#   line   0.4  2.6  2.0  2.0
#   line   1.0  1.0  2.0  2.0
    line   1.0  1.0  3.0  0.8
    line   2.0  2.0  4.0  2.0
#   line   2.0  2.0  3.0  0.8
#   line   4.0  2.0  3.0  3.4
#   line   4.0  2.0  4.6  4.0
#   line   4.0  2.0  4.8  1.6
#   line   4.0  2.0  3.0  0.8
    line   3.0  3.4  4.6  4.0
    line   4.6  4.0  4.8  1.6
    line   4.8  1.6  3.0  0.8
#   line   4.8  1.6  4.6  0.2
    line   3.0  0.8  4.6  0.2
#
#  Blue things
#
    fill_rgb 0.0 0.0 1.0

    circle_fill 1.0  4.4  0.25
    circle_fill 0.4  2.6  0.25
    circle_fill 1.0  1.0  0.25
    circle_fill 2.0  2.0  0.25
    circle_fill 4.0  2.0  0.25
    circle_fill 3.0  3.4  0.25
    circle_fill 4.6  4.0  0.25
    circle_fill 4.8  1.6  0.25
    circle_fill 3.0  0.8  0.25
    circle_fill 4.6  0.2  0.25
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0
#
#  Node labels
#
    font_size 0.40

    moveto  1.0  4.70
    label A
    moveto  0.3  2.90
    label B
    moveto  1.0  1.25
    label C
    moveto  2.0  2.25
    label D
    moveto  3.9  2.35
    label E
    moveto  3.0  3.65
    label F
    moveto  4.6  4.25
    label G
    moveto  4.8  1.85
    label H
    moveto  3.0  1.10
    label I
    moveto  4.4  0.45
    label J
#
#  Line labels
#
    font_size 0.30

    moveto 0.5 3.5
    label  3
    moveto 2.0 4.0
    label  2
    moveto 0.7  1.8
    label  17
    moveto 1.2 2.4
    label  16
    moveto 1.4  1.6
    label  8
    moveto 2.0 1.0
    label  18
    moveto 3.0 2.1
    label  11
    moveto 2.5 1.4
    label  4
    moveto 3.5 2.7
    label  1
    moveto 4.1 3.0
    label  6
    moveto 4.4 1.9
    label  5
    moveto 3.3  1.4
    label  10
    moveto 3.8  3.8
    label  7
    moveto 4.7  2.8
    label  15
    moveto 3.9  1.3
    label  12
    moveto 4.7 0.9
    label  13
    moveto 3.8 0.6
    label  9

  endpage

endfile
