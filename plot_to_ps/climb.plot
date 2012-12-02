# climb.plot  28 September 2001
#
#  Illustrates the problem of a train climbing a hill.
#
file climb.ps

  space 0.0 0.0 100.0 100.0

  page

    line_width 3

    line_rgb 0.0 0.5 0.5

    moveto 10.0 10.0
    drawto 90.0 10.0
    drawto 90.0 50.0
    drawto 10.0 10.0

    line_rgb 0.5 0.0 0.0

    arrow 15.0 32.5 85.0 67.5

    font_size 0.35

    line_rgb 0.0 0.0 0.0

    moveto 95.0 30.0
    label ?

    moveto 10.0 90.0
    label Distance is measured along the slope!

  endpage

endfile
