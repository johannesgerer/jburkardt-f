# triangle.plot  12 April 2005
#
file triangle.ps

  space -0.05 -0.05 1.05 0.90

  page

    line_thick 5

    fill_rgb 0.60 0.95 0.90

    polygon_fill
      0.0  0.0
      1.0  0.2
      0.3  0.8
    endpolygon

    line_rgb 0.0 0.0 0.0

    moveto  0.0  0.0
    drawto  1.0  0.2
    drawto  0.3  0.8
    drawto  0.0  0.0

    line_rgb 1.0 0.0 0.0

    moveto  0.00  0.00
    drawto  0.65  0.50
    moveto  1.00  0.20
    drawto  0.15  0.40
    moveto  0.30  0.80
    drawto  0.50  0.10

    fill_rgb 0.0 0.0 0.0
    circle_fill 0.0 0.0 0.025
    circle_fill 1.0 0.2 0.025
    circle_fill 0.3 0.8 0.025

    fill_rgb 1.0 0.0 0.0
    circle_fill 0.65 0.50 0.025
    circle_fill 0.50 0.10 0.025
    circle_fill 0.15 0.40 0.025

    fill_rgb 0.0 0.0 1.0
    circle_fill 0.43 0.33 0.025

    font_size 0.35
    fill_rgb 0.0 0.0 0.0

    moveto  0.50  0.33
    label Centroid
    moveto 0.68 0.50
    label Midpoint
    moveto 0.53 0.05
    label Midpoint
    moveto 0.17 0.40
    label Midpoint

  endpage

endfile
