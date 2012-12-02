# triangle_distance.plot  21 February 2011
#
file triangle_distance.ps

  space -2.5 -1.5 6.5 7.5

  page

    line_rgb 0.5 0.5 0.5
    grid -2.0 -1.0 6.0 7.0 9 9

    line_thick 4
    line_rgb 0.0 0.0 0.0
    moveto -2.0 0.0
    drawto  6.0 0.0
    moveto  0.0 -1.0
    drawto  0.0 7.0
    circle_fill 0.0 0.0 0.125

    line_thick 3

    fill_rgb 0.9 0.95 0.35
    triangle_fill 4.0 1.0 3.0 5.0 0.0 2.0

    fill_rgb 0.5 0.7 0.5

    circle_fill  4.0   1.0  0.125
    circle_fill  3.0   5.0  0.125
    circle_fill  0.0   2.0  0.125
#
#  Grid lines on mapped triangle.
#
    line_thick 3

    fill_rgb 0.2 0.2 0.9
    moveto -1.42  3.42
    drawto  1.58  6.42
    moveto  0.0   2.0
    drawto -1.42  3.42
    moveto  3.0   5.0
    drawto  1.58  6.42
    moveto  0.0   2.0
    drawto -0.48  0.06
    moveto  4.0   1.0
    drawto  3.52 -0.94
    moveto -0.48  0.06
    drawto  3.52 -0.94
    moveto  3.0  5.0
    drawto  4.94 5.48
    moveto  4.0  1.0
    drawto  5.94 1.48
    moveto  4.94 5.48
    drawto  5.94 1.48
    arc     0.0  2.0  2.0 135.0 -104.0
    arc     3.0  5.0  2.0 14.0 135.0
    arc     4.0  1.0  2.0 -104.0 14.0

    fill_rgb 0.0 0.7 0.0
    moveto -0.71  2.71
    drawto  2.29  5.71
    moveto  0.0   2.0
    drawto -0.71  2.71
    moveto  3.0   5.0
    drawto  2.29  5.71
    moveto  0.0   2.0
    drawto -0.24  1.03
    moveto  4.0   1.0
    drawto  3.76  0.03
    moveto -0.24  1.03
    drawto  3.76  0.03
    moveto  3.0  5.0
    drawto  3.97 5.24
    moveto  4.0  1.0
    drawto  4.97 1.24
    moveto  3.97 5.24
    drawto  4.97 1.24
    arc     0.0  2.0  1.0 135.0 -104.0
    arc     3.0  5.0  1.0 14.0 135.0
    arc     4.0  1.0  1.0 -104.0 14.0
#
#  Draw triangle boundaries.
#
    fill_rgb 0.0 0.0 0.0

    moveto   4.0  1.0
    drawto   3.0  5.0
    drawto   0.0  2.0
    drawto   4.0  1.0
#
#  Draw test points.
#
    fill_rgb 1.0 0.1 0.1
    circle_fill -2.0   1.0  0.125
    circle_fill  0.0   2.0  0.125
    circle_fill  2.0   3.0  0.125
    circle_fill  4.0   4.0  0.125
    circle_fill  6.0   5.0  0.125

  endpage

endfile
