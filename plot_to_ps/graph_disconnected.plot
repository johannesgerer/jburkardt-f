# graph_disconnected.plot  30 December 2010
#
file graph_disconnected.ps

  space -0.5 -0.5 4.5 2.5

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 4.0 2.0 5 3
#
#  Red things (lines)
#
    fill_rgb 0.9 0.0 0.1
    line_width 3

    line  0.0  2.0  1.0  1.0
    line  0.0  2.0  0.0  0.0
    line  2.0  2.0  4.0  2.0
    line  4.0  2.0  4.0  0.0
    line  1.0  1.0  3.0  1.0
    line  1.0  1.0  0.0  0.0
#
#  Blue things (nodes)
#
    fill_rgb 0.0 0.1 0.9

    circle_fill  0.0  2.0  0.20
    circle_fill  2.0  2.0  0.20
    circle_fill  4.0  2.0  0.20
    circle_fill  1.0  1.0  0.20
    circle_fill  3.0  1.0  0.20
    circle_fill  0.0  0.0  0.20
    circle_fill  2.0  0.0  0.20
    circle_fill  4.0  0.0  0.20
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0
#
#  Node labels
#
    font_size 0.40

    moveto  0.0  2.25
    label A
    moveto  2.0  2.25
    label B
    moveto  4.0  2.25
    label C
    moveto  1.0  1.25
    label D
    moveto  3.0  1.25
    label E
    moveto  0.0  0.25
    label F
    moveto  2.0  0.25
    label G
    moveto  4.0  0.25
    label H

  endpage

endfile
