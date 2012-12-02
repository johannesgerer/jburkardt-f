# bill_060.plot  05 January 2006
#
file bill_060.ps

  space 0.0 0.0 2.0 1.0

  page

    line_thick 2

    line_rgb 0.9 0.0 0.0

    moveto  0.2  0.15
    drawto  0.8  0.15
    drawto  0.5  0.75
    drawto  0.2  0.15

    line_rgb 0.0 0.0 1.0

    circle_fill 0.2 0.15 0.02
    circle_fill 0.8 0.15 0.02
    circle_fill 0.5 0.75 0.02
    circle_fill 0.5 0.35 0.02

    moveto 0.2 0.5
    line_rgb 0.0 0.0 0.0
    font_size 0.25
    label X

    moveto 0.25 0.55
    line_rgb 0.0 0.0 0.0
    font_size 0.25
    label h

    line_rgb 0.9 0.0 0.0

    moveto  1.2  0.15
    drawto  1.8  0.15
    drawto  1.5  0.75
    drawto  1.2  0.15

    line_rgb 0.0 0.0 1.0

    circle_fill 1.2 0.15 0.02
    circle_fill 1.8 0.15 0.02
    circle_fill 1.5 0.75 0.02

    moveto 1.2 0.5
    line_rgb 0.0 0.0 0.0
    font_size 0.25
    label Q

    moveto 1.25 0.55
    line_rgb 0.0 0.0 0.0
    font_size 0.25
    label h

  endpage

endfile
