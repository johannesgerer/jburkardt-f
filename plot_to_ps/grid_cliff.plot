# grid_cliff.plot  26 August 2009
#
file grid_cliff.ps

  space -0.1 -0.1 2.2 2.2

  page
#
#  Shade the region.
#
    line_thick 2

    fill_rgb 0.95 0.95 0.6

    polygon_fill
      0.000 0.000
      2.000 0.000
      2.000 2.000
      0.000 2.000
    endpolygon
#
#  Draw the grid lines.
#
    line_rgb 0.9 0.0 0.0
#
#  3 Diagonals.
#
    moveto 0.0 1.0
    drawto 1.0 2.0

    moveto 0.0 0.0
    drawto 2.0 2.0

    moveto 1.0 0.0
    drawto 2.0 1.0
#
#  3 Horizontals.
#
    moveto 0.0 0.0
    drawto 2.0 0.0

    moveto 0.0 1.0
    drawto 2.0 1.0

    moveto 0.0 2.0
    drawto 2.0 2.0
#
#  3 Verticals.
#
    moveto 0.0 0.0
    drawto 0.0 2.0

    moveto 1.0 0.0
    drawto 1.0 2.0

    moveto 2.0 0.0
    drawto 2.0 2.0
#
#  Green nodes.
#
    line_rgb 0.0 1.0 0.0

    circle_fill 0.0 0.0 0.04
    circle_fill 0.5 0.0 0.04
    circle_fill 1.0 0.0 0.04
    circle_fill 1.5 0.0 0.04
    circle_fill 2.0 0.0 0.04
    circle_fill 0.0 0.5 0.04
    circle_fill 0.5 0.5 0.04
#
#  Blue nodes
#
    line_rgb 0.0 0.0 1.0

    circle_fill 1.0 0.5 0.04
    circle_fill 1.5 0.5 0.04
    circle_fill 2.0 0.5 0.04
    circle_fill 0.0 1.0 0.04
    circle_fill 0.5 1.0 0.04
    circle_fill 1.0 1.0 0.04
#
#  Red nodes.
#
    line_rgb 1.0 0.0 0.0

    circle_fill 1.5 1.0 0.04
    circle_fill 2.0 1.0 0.04
    circle_fill 0.0 1.5 0.04
    circle_fill 0.5 1.5 0.04
    circle_fill 1.0 1.5 0.04
    circle_fill 1.5 1.5 0.04
#
#  Black Nodes.
#
    line_rgb 0.0 0.0 0.0

    circle_fill 2.0 1.5 0.04
    circle_fill 0.0 2.0 0.04
    circle_fill 0.5 2.0 0.04
    circle_fill 1.0 2.0 0.04
    circle_fill 1.5 2.0 0.04
    circle_fill 2.0 2.0 0.04
#
#  Text
#
    font_size 0.25
    line_rgb 0.0 0.0 0.0

    moveto 0.05 0.05
    label 1
    moveto 0.55 0.05
    label 2
    moveto 1.05 0.05
    label 3
    moveto 1.55 0.05
    label 4
    moveto 2.05 0.05
    label 5

    moveto 0.05 0.55
    label 6
    moveto 0.55 0.55
    label 7
    moveto 1.05 0.55
    label 8
    moveto 1.55 0.55
    label 9
    moveto 2.05 0.55
    label 10

    moveto 0.05 1.05
    label 11
    moveto 0.55 1.05
    label 12
    moveto 1.05 1.05
    label 13
    moveto 1.55 1.05
    label 14
    moveto 2.05 1.05
    label 15

    moveto 0.05 1.55
    label 16
    moveto 0.55 1.55
    label 17
    moveto 1.05 1.55
    label 18
    moveto 1.55 1.55
    label 19
    moveto 2.05 1.55
    label 20

    moveto 0.05 2.05
    label 21
    moveto 0.55 2.05
    label 22
    moveto 1.05 2.05
    label 23
    moveto 1.55 2.05
    label 24
    moveto 2.05 2.05
    label 25

    font_size 0.55

    moveto 0.33 0.66
    label 1
    moveto 0.66 0.33
    label 2
    moveto 1.33 0.66
    label 3
    moveto 1.66 0.33
    label 4
    moveto 0.33 1.66
    label 5
    moveto 0.66 1.33
    label 6
    moveto 1.33 1.66
    label 7
    moveto 1.66 1.33
    label 8
    
  endpage

endfile
