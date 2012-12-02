# snake_polygon.plot  05 February 2011

file snake_polygon.ps

  space -5.0 -5.0 55.0 25.0

  page

    line_rgb 0.0 0.0 0.0
    line_width 4

    polygon
      10   0
      20  10
      30   0
      40  10
      50   0
      50  10
      40  20
      30  10
      20  20
      10  10
       0  20
       0  10
    endpolygon

    line_rgb 0.9 0.5 0.9

    polygon_fill
      10   0
      20  10
      30   0
      40  10
      50   0
      50  10
      40  20
      30  10
      20  20
      10  10
       0  20
       0  10
    endpolygon

  endpage

endfile
