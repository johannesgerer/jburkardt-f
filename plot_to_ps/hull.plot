# hull.plot  12 September 2000
#
file hull.eps

  space -3.0 -1.0 12.0 10.0

  page

    line_rgb 0.2 0.2 0.2
    grid -3.0 -1.0 12.0 10.0 16 12

    line_rgb 0.5 0.5 0.5

    moveto  0.0  0.0
    drawto  2.0  2.0
    drawto -1.0  3.0
    drawto -2.0  2.0
    drawto  0.0  0.0

    moveto  8.0  2.0
    drawto  9.0  5.0
    drawto  7.0  4.0
    drawto  8.0  2.0

    line_rgb 1.0 0.0 0.0

    moveto  5.0  6.0
    drawto  6.0  7.0
    drawto  8.0  8.0
    drawto 11.0  7.0
    drawto 10.0  4.0
    drawto  8.0  2.0
    drawto  6.0  4.0
    drawto  5.0  6.0

    moveto  0.0  8.0
    font_size 0.50
    label Minkowski sum in red.

  endpage

endfile
