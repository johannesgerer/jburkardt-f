# triangulation_02.plot  04 June 2009
#
file triangulation_02.ps

  space 0.0 0.0 2.0 2.0

  page

    line_thick 2

    fill_rgb 0.7 0.9 0.7

    polygon_fill
      0.000 0.000
      1.000 0.000
      0.500 0.866
    endpolygon

    polygon_fill
      1.000 1.732
      1.500 0.866
      0.500 0.866
    endpolygon

    polygon_fill
      1.000 0.000
      2.000 0.000
      1.500 0.866
    endpolygon

    line_rgb 0.9 0.0 0.0

    moveto 0.000 0.000
    drawto 1.000 0.000

    moveto 0.000 0.000
    drawto 0.500 0.866

    moveto 0.500 0.866
    drawto 1.000 1.732

    moveto 0.500 0.866
    drawto 1.500 0.866

    moveto 1.500 0.866
    drawto 1.000 1.732

    moveto 1.000 0.000
    drawto 0.500 0.866

    moveto 1.000 0.000
    drawto 1.500 0.866

    moveto 1.000 0.000
    drawto 2.000 0.000

    moveto 2.000 0.000
    drawto 1.500 0.866

    line_rgb 0.0 0.0 1.0

    circle_fill 0.000 0.000 0.02
    circle_fill 1.000 0.000 0.02
    circle_fill 2.000 0.000 0.02
    circle_fill 0.500 0.866 0.02
    circle_fill 1.500 0.866 0.02
    circle_fill 1.000 1.732 0.02


  endpage

endfile
