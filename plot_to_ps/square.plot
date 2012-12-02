# square.plot  11 April 2005
#
file square.ps

  space -0.05 -0.05 1.05 1.05

  page

    moveto  0.0  0.0
    drawto  1.0  0.0
    drawto  1.0  1.0
    drawto  0.0  1.0
    drawto  0.0  0.0

    fill_rgb 0.8 0.8 1.0
    square_fill 0.5 0.5 0.5

    line_thick 5

    line_rgb 1.0 0.0 0.0
    line 0.0 0.0 1.0 0.0
    line 1.0 0.0 1.0 1.0
    line 1.0 1.0 0.0 1.0
    line 0.0 1.0 0.0 0.0

    fill_rgb 0.5 0.7 0.5
    circle_fill 0.0 0.0 0.05
    circle_fill 1.0 0.0 0.05
    circle_fill 0.0 1.0 0.05
    circle_fill 1.0 1.0 0.05

  endpage

endfile
