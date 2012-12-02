# graph_simple.plot  23 December 2010
#
file graph_simple.ps

  space -0.5 -0.5 3.5 3.5

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 3.0 3.0 4 4
#
#  Red things
#
    fill_rgb 1.0 0.0 0.0
    line_width 3

    line   2.0  1.0  1.0  0.0
    line   2.0  1.0  0.0  2.0
    line   1.0  0.0  0.0  2.0
    line   0.0  2.0  3.0  3.0
#
#  Blue things
#
    fill_rgb 0.0 0.0 1.0

    circle_fill  2.0  1.0 0.25
    circle_fill  1.0  0.0 0.25
    circle_fill  0.0  2.0 0.25
    circle_fill  3.0  3.0 0.25
    circle_fill  3.0  0.0 0.25
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0

    font_size 0.40

    moveto 2.0 1.35
    label A
    moveto 1.0 0.35
    label B
    moveto 0.0 2.35
    label C
    moveto 3.0 3.35
    label D
    moveto 3.0 0.35
    label E


  endpage

endfile
