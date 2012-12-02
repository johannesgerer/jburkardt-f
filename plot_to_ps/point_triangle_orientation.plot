# point_triangle_orientation.plot  03 February 2007
#
file point_triangle_orientation.ps

  space -1.5 -1.5 5.5 6.5

  page

    line_rgb 0.5 0.5 0.5

    grid -1.0 -1.0 5.0 6.0 7 8

    line_thick 3
#
#  Blue things are neutral.
#
    line_rgb 0.0 0.0 1.0

    moveto  4.0  0.0
    drawto  3.0  4.0
    drawto  0.0  1.0
    drawto  4.0  0.0

    circle       4.0  0.0 0.125
    circle       3.0  4.0 0.125
    circle       0.0  1.0 0.125

    circle_fill  2.0  3.0 0.125

    circle_fill  3.0  4.0 0.125
#
#  Red things are negative
#
    line_rgb 1.0 0.0 0.0

    circle_fill  1.0  4.0 0.125
    circle_fill  3.0  5.0 0.125
    circle_fill  4.0  1.0 0.125
    circle_fill  4.0  5.0 0.125
#
#  Green things are positive.
#
    line_rgb 0.0 1.0 0.0


    circle_fill  3.0  2.0 0.125
    circle_fill  2.0  1.0 0.125
#
#  Black things.
#
    line_rgb 0.0 0.0 0.0

    moveto  0.0 -1.0
    drawto  0.0  6.0
    moveto -1.0  0.0
    drawto  5.0  0.0

    font_size 0.35
    moveto 1.0 -0.5
    label T={{4,0},{3,4},{0,1}}

  endpage

endfile
