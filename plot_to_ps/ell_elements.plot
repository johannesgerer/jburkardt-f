# triangulation_04.plot  04 June 2009
#
file triangulation_04.ps

  space -0.1 -0.1 2.1 2.1

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 2.0 2.0 9 9

    line_thick 2

    fill_rgb 0.7 0.9 0.7

    polygon_fill
      0.000 0.000
      2.000 0.000
      2.000 1.000
      1.000 1.000
      1.000 2.000
      0.000 2.000
    endpolygon

    line_rgb 0.9 0.0 0.0

    moveto 0.0 0.0
    drawto 1.0 0.0

    moveto 1.0 0.0
    drawto 2.0 0.0

    moveto 0.0 1.0
    drawto 1.0 1.0

    moveto 1.0 1.0
    drawto 2.0 1.0

    moveto 0.0 2.0
    drawto 1.0 2.0


    moveto 0.0 0.0
    drawto 0.0 1.0

    moveto 0.0 1.0
    drawto 0.0 2.0

    moveto 1.0 0.0
    drawto 1.0 1.0

    moveto 1.0 1.0
    drawto 1.0 2.0


    moveto 2.0 0.0
    drawto 2.0 1.0


    moveto 0.0 1.0
    drawto 1.0 0.0

    moveto 1.0 1.0
    drawto 0.0 2.0

    moveto 2.0 0.0
    drawto 1.0 1.0

    moveto 0.0 0.5
    drawto 2.0 0.5

    moveto 0.0 1.5
    drawto 1.0 1.5

    moveto 0.5 0.0
    drawto 0.5 2.0

    moveto 1.5 0.0
    drawto 1.5 1.0

    moveto 0.0 0.5
    drawto 0.5 0.0

    moveto 0.0 1.5
    drawto 1.5 0.0

    moveto 1.5 1.0
    drawto 2.0 0.5

    moveto 0.5 2.0
    drawto 1.0 1.5

    line_rgb 0.9 0.0 0.0

    moveto 0.0 0.0
    drawto 2.0 0.0
    drawto 2.0 1.0
    drawto 1.0 1.0
    drawto 1.0 2.0
    drawto 0.0 2.0
    drawto 0.0 0.0

    line_rgb 0.0 0.0 1.0

    circle_fill 0.0 0.0 0.02
    circle_fill 0.5 0.0 0.02
    circle_fill 1.0 0.0 0.02
    circle_fill 1.5 0.0 0.02
    circle_fill 2.0 0.0 0.02

    circle_fill 0.0 0.5 0.02
    circle_fill 0.5 0.5 0.02
    circle_fill 1.0 0.5 0.02
    circle_fill 1.5 0.5 0.02
    circle_fill 2.0 0.5 0.02

    circle_fill 0.0 1.0 0.02
    circle_fill 0.5 1.0 0.02
    circle_fill 1.0 1.0 0.02
    circle_fill 1.5 1.0 0.02
    circle_fill 2.0 1.0 0.02

    circle_fill 0.0 1.5 0.02
    circle_fill 0.5 1.5 0.02
    circle_fill 1.0 1.5 0.02

    circle_fill 0.0 2.0 0.02
    circle_fill 0.5 2.0 0.02
    circle_fill 1.0 2.0 0.02

  endpage

endfile
