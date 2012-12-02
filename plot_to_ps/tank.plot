# tank.plot  17 February 2001
#
#  Illustrates the problem of determining the water volume
#  in a tank with circular cross section.
#
file tank.ps

  space 0.0 0.0 100.0 100.0

  page

    circle 50.0 50.0 50.0

    point 50.0 50.0

    fill_rgb 0.00 0.00 1.0
    sector_fill 50.0 50.0 50.0 240.0 300.0

    fill_rgb 0.95 0.95 1.0
    triangle_fill 25.0  6.7 50.0 50.0 75.0  6.7

    line_rgb 0.0 0.8 0.0

    moveto 25.0  6.7
    drawto 50.0 50.0

    moveto 50.0  6.7
    drawto 50.0 50.0

    moveto 75.0  6.7
    drawto 50.0 50.0

    line_rgb 0.0 0.0 1.0
    moveto 25.0  6.7
    drawto 75.0  6.7

    fill_rgb 0.0 0.0 0.0

    font_size 0.25

    moveto 36 30
    label R

    moveto 47 30
    label R - H

    moveto 49 45
    label A

    font_size 0.50

    moveto 17 75
    label The water tank problem

  endpage

endfile
