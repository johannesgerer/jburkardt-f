# line_dist.plot  30 January 2011
#
file line_dist.ps

  space -4.0 -3.0 7.0 8.0

  page

    line_rgb 0.8 0.8 0.8
    grid -4.0 -3.0 7.0 8.0 12 12
#
#  Black things.
#
    line_width 3

    fill_rgb 0.0 0.0 0.0
    moveto -4.0  8.0
    drawto  7.0 -3.0

    moveto  0.0  4.0
    drawto  4.0  8.0

#   moveto  7.0  2.0
#   drawto  2.5  6.5

    font_size 0.40

    moveto 6.5 2.5
    label Q

    moveto -1.0 4.0
    label P1

    moveto 3.5 1.0
    label P2

    moveto 3.0 6.5
    label Normal vector

    moveto 2.5 6.8
    label t

    moveto 4.5 0.0
    label s
#
#  Blue things
#
    fill_rgb 0.0 0.0 1.0

    circle_fill  0.0  4.0  0.40
    circle_fill  3.0  1.0  0.40
#
#  Green things
#
    fill_rgb 0.0 1.0 0.0

    circle_fill  2.5  6.5  0.25
    circle_fill  4.5 -0.5  0.25
#
#  Gray things (axis)
#
    fill_rgb 0.5 0.5 0.5

    line_width 3

    moveto  0.0 -3.0
    drawto  0.0  8.0
    moveto -4.0  0.0
    drawto  7.0  0.0
#
#  Red things.
#
    fill_rgb 1.0 0.0 0.0

#   moveto  0.0  4.0
#   drawto  2.5  6.5

    moveto  0.0  4.0
    drawto  4.5 -0.5

    moveto  0.0  4.0
    drawto  7.0  2.0

    moveto  7.0  2.0
    drawto  4.5 -0.5

    circle_fill 7.0 2.0 0.25



  endpage

endfile
