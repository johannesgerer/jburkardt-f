# graph_tsp.plot  30 December 2010
#
file graph_tsp.ps

  space -0.5 -0.5 6.5 6.5

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 6.0 6.0 7 7
#
#  Red things (lines)
#
    fill_rgb 1.0 0.0 0.0
    line_width 3

    line  3.0  6.0  6.0  5.0
    line  3.0  6.0  5.0  0.0
    line  3.0  6.0  1.0  1.0
    line  3.0  6.0  0.0  4.0
    line  6.0  5.0  5.0  0.0
    line  6.0  5.0  1.0  1.0
    line  6.0  5.0  0.0  4.0
    line  5.0  0.0  1.0  1.0
    line  5.0  0.0  0.0  4.0
    line  1.0  1.0  0.0  4.0
#
#  Blue things (nodes)
#
    fill_rgb 0.0 0.0 1.0

    circle_fill  3.0  6.0  0.25
    circle_fill  6.0  5.0  0.25
    circle_fill  5.0  0.0  0.25
    circle_fill  1.0  1.0  0.25
    circle_fill  0.0  4.0  0.25
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0
#
#  Node labels
#
    font_size 0.40

    moveto  3.3  6.2
    label A
    moveto  6.0  5.3
    label B
    moveto  5.3  0.0
    label C
    moveto  1.0  0.4
    label D
    moveto  0.0  4.3
    label E
#
#  Line labels
#
    font_size 0.30

    moveto  4.5  5.6
    label 3
    moveto  4.2  2.8
    label 4
    moveto  2.2  3.5
    label 2
    moveto  1.5  5.2
    label 7
    moveto  5.6  2.5
    label 4
    moveto  3.2  2.9
    label 6
    moveto  3.0  4.6
    label 3
    moveto  3.0  0.6
    label 5
    moveto  2.9  1.7
    label 8
    moveto  0.6  2.5
    label 6

  endpage

endfile
