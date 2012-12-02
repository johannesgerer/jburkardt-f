# graph_tsp.plot  30 December 2010
#
file graph_tsp.ps

  space -0.5 -0.5 6.5 4.5

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 6.0 4.0 7 5
#
#  Green things (lines)
#
    fill_rgb 0.0 1.0 0.0
    line_width 3

    line   0.0  2.0  6.0  1.0
    arrow  0.0  2.0  3.0  1.5
    line   1.0  0.0  4.0  4.0
    arrow  1.0  0.0  2.5  2.0
#
#  Red things (lines)
#
    fill_rgb 1.0 0.0 0.0
    line_width 3

    line   0.0  2.0  1.0  0.0
    arrow  0.0  2.0  0.5  1.0
    line   6.0  1.0  4.0  4.0
    arrow  6.0  1.0  5.0  2.5
#
#  Blue things (nodes)
#
    fill_rgb 0.0 0.0 1.0

    circle_fill  0.0  2.0  0.25
    circle_fill  6.0  1.0  0.25
    circle_fill  1.0  0.0  0.25
    circle_fill  4.0  4.0  0.25
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0
#
#  Node labels
#
    font_size 0.40

    moveto  0.0  2.3
    label A
    moveto  6.0  1.3
    label B
    moveto  1.0  0.3
    label C
    moveto  4.0  4.3
    label D

  endpage

endfile
