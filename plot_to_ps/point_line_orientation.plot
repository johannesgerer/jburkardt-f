# point_line_orientation.plot  03 February 2007
#
file point_line_orientation.ps

  space -1.0 -1.0 6.0 12.0

  page

    grid 0.0 0.0 5.0 11.0 6 12

    font_size 0.35
    line_thick 3
#
#  Blue things are "neutral"
#
    line_rgb 0.0 0.0 1.0

    circle_fill 1.0  3.0 0.125
    moveto 2.0 3.0
    circle_fill 4.0  9.0 0.125
    moveto 0.5 9.0

    moveto  0.0  1.0
    drawto  5.0 11.0

    circle_fill  3.0  7.0 0.125
#
#  Red things are negative
#
    line_rgb 1.0 0.0 0.0

    circle_fill  3.0  3.0 0.125
    circle_fill  3.0  4.0 0.125
    circle_fill  3.0  5.0 0.125
    circle_fill  3.0  6.0 0.125

    moveto 1.5 2.0
    label Negative distance
#
#  Green things are positive.
#
    line_rgb 0.3 1.0 0.3

    circle_fill  3.0  8.0 0.125
    circle_fill  3.0  9.0 0.125
    circle_fill  3.0 10.0 0.125

    moveto 0.5 8.5
    label Positive distance
#
#  Black things are titles.
#
    line_rgb 0.0 0.0 0.0
    moveto 0.5 0.5
    label L1={{1,3}{4,9}}

  endpage

endfile
