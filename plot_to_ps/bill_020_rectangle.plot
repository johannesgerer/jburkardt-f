# bill_020_rectangle.plot  08 August 2006
#
file bill_020_rectangle.ps

  space -0.5 -0.5 2.5 3.5

  page

    line_thick 2

    line_rgb 0.9 0.0 0.0

    moveto  0.00   0.00
    drawto  2.00   0.00
    moveto  0.00   1.00
    drawto  2.00   1.00
    moveto  0.00   2.00
    drawto  2.00   2.00

    moveto  0.00   0.00
    drawto  0.00   2.00
    moveto  1.00   0.00
    drawto  1.00   2.00
    moveto  2.00   0.00
    drawto  2.00   2.00

    moveto  0.00   1.00
    drawto  1.00   2.00
    moveto  0.00   0.00
    drawto  2.00   2.00
    moveto  1.00   0.00
    drawto  2.00   1.00

    line_rgb 0.0 0.0 1.0

    circle_fill  0.00   0.00 0.03
    circle_fill  1.00   0.00 0.03
    circle       2.00   0.00 0.03
    circle_fill  0.00   1.00 0.03
    circle_fill  1.00   1.00 0.03
    circle_fill  2.00   1.00 0.03
    circle       0.00   2.00 0.03
    circle_fill  1.00   2.00 0.03
    circle_fill  2.00   2.00 0.03

    line_rgb 0.0 1.0 0.0

    moveto  0.0   0.0
    drawto  1.0   2.8
    moveto  0.0   1.0
    drawto  1.0   2.8
    moveto  1.0   0.0
    drawto  1.0   2.8
    moveto  1.0   1.0
    drawto  1.0   2.8
    moveto  1.0   2.0
    drawto  1.0   2.8
    moveto  2.0   1.0
    drawto  1.0   2.8
    moveto  2.0   2.0
    drawto  1.0   2.8

    line_rgb 0.0 1.0 0.0

    circle_fill  1.00   2.80 0.03


  endpage

endfile
