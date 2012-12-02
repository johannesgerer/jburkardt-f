# circle.plot  02 August 1999
#
file circle.ps

  space -10.0 -10.0 10.0 10.0

  page

    moveto -9.0 -9.0
    drawto  9.0 -9.0
    drawto  9.0  9.0
    drawto -9.0  9.0
    drawto -9.0 -9.0

    circle 2.0 -5.0 4.0

    ellipse 0.0 -6.0 5.0 1.0 20

    grid -5.0 -5.0 5.0 5.0 11 11

    point 0.0 0.0

    line_rgb 0.0 1.0 0.0

    arc 3.0 3.0 2.0 0.0 90.0

    line_rgb 0.5 0.5 0.0
    star 2.0 5.0 1.0

    line_rgb 0.0 0.5 0.5
    square -3.0 7.0 1.0

    fill_rgb 1.0 0.8 0.8
    square_fill -4.0 6.0 1.0

    fill_rgb 0.0 0.6 0.2
    ellipse_fill -5.0 -6.0 2.0 4.0 20

    moveto -5.0 -5.0
    line_rgb 0.0 0.0 0.0
    font_size 0.25
    label This is "circle.plot".

    moveto -5.0 -2.0
    font_size 0.125
    label This is small print.

  endpage

endfile
