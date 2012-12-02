# graph_dijkstra.plot  06 January 2011
#
file graph_dijkstra.ps

  space -0.5 -0.5 5.5 3.5

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 5.0 3.0 6 4
#
#  Red things (lines)
#
    fill_rgb 0.9 0.0 0.1
    line_width 3

    line  0.0  1.0  2.0  1.0
    line  0.0  1.0  1.0  3.0
    line  2.0  1.0  1.0  3.0
    line  1.0  3.0  4.0  2.0
    line  2.0  1.0  4.0  2.0
    line  2.0  1.0  5.0  0.0
    line  2.0  1.0  1.0  0.0
    line  5.0  0.0  1.0  0.0
#
#  Blue things (nodes)
#
    fill_rgb 0.0 0.1 0.9

    circle_fill  0.0  1.0  0.20
    circle_fill  2.0  1.0  0.20
    circle_fill  1.0  3.0  0.20
    circle_fill  4.0  2.0  0.20
    circle_fill  5.0  0.0  0.20
    circle_fill  1.0  0.0  0.20
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0
#
#  Node labels
#
    font_size 0.40

    moveto  -0.25  1.25
    label A
    moveto  2.0  1.25
    label B
    moveto  1.0  3.25
    label C
    moveto  4.0  2.25
    label D
    moveto  5.0  0.25
    label E
    moveto  1.0  0.25
    label F

    moveto  1.0  1.1
    label 40
    moveto 0.0 2.0
    label 15
    moveto  1.5  2.0
    label 20
    moveto  2.5  2.5
    label 100
    moveto  2.8  1.6
    label 10
    moveto  3.5  0.6
    label 25
    moveto 1.3 0.5
    label 6
    moveto 3.0  0.1
    label 8

  endpage

endfile
