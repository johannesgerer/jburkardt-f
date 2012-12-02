# line_perp.plot  26 January 2011
#
file line_perp.ps

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

    circle_fill -3.0  7.0  0.25
    circle_fill  1.0  3.0  0.25
    circle_fill  1.5  2.5  0.25
    circle_fill  2.0  2.0  0.25
    circle_fill  4.5 -0.5  0.25
    circle_fill  6.0 -2.0  0.25
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0

    line_width 3

    moveto  0.0 -3.0
    drawto  0.0  8.0
    moveto -4.0  0.0
    drawto  7.0  0.0



    font_size 0.40

    moveto -3.0 7.5
    label -1.0
    moveto 0.0 4.5
    label P1 = 0.0
    moveto 1.0 3.5
    label 1/3
    moveto 1.5 3.0
    label 1/2
    moveto 2.0 2.5
    label 2/3
    moveto 3.0 1.5
    label P2 = 1.0
    moveto 4.5 0.0
    label 1.5
    moveto 6.0 -1.5
    label 2.0
#
#  Another red thing.
#
    fill_rgb 1.0 0.0 0.0
    moveto  7.0  2.0
    drawto  2.0 -3.0
    circle_fill 7.0 2.0 0.25
    moveto 5.5 2.3
    label Q=1.5


  endpage

endfile
