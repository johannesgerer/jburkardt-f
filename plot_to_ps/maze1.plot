# maze1.plot  21 January 2011
#
file maze.ps

  space -0.50 -0.50 10.50 11.50

  page

    line_rgb 0.7 0.7 0.7

    grid 0.0 0.0 10.0 11.0 11 12

    line_thick 3

    fill_rgb 0.90 0.90 0.40

    polygon_fill
      0  0
     10  0
     10  3
      5  3
      5  6
      8  6
      8  7
      4  7
      4  3
      2  3
      2  2
      9  2
      9  1
      1  1
      1  4
      3  4
      3  5
      1  5
      1  9
      0  9
    end

    polygon_fill
      0 11
      0 10
      2 10
      2  6
      3  6
      3  8
      5  8
      5  9
      3  9
      3 10
      9 10
      9  9
      6  9
      6  8
      9  8
      9  5
      6  5
      6  4
     10  4
     10 11
    end
#
#  Label the plot.
#
    font_size 0.35
    line_rgb 0.0 0.0 0.0

    moveto 0.0 9.25
    label Entrance
    moveto 9.0 3.25
    label Exit

  endpage

endfile
