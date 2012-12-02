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
#  Red things
#
    line_rgb 1.0 0.0 0.0

    moveto 0.5 9.5
    drawto 1.5 9.5
    drawto 1.5 5.5
    drawto 3.5 5.5

    moveto 3.5 5.5
    drawto 3.5 3.5
    drawto 1.5 3.5
    drawto 1.5 1.5
    drawto 8.5 1.5

    moveto 3.5 5.5
    drawto 3.5 7.5
    drawto 5.5 7.5

    moveto 5.5 7.5
    drawto 5.5 9.5

    moveto 5.5 9.5
    drawto 8.5 9.5

    moveto 5.5 9.5
    drawto 3.5 9.5

    moveto 5.5 7.5
    drawto 8.5 7.5
    drawto 8.5 5.5
    drawto 5.5 5.5
    drawto 5.5 3.5
    drawto 9.5 3.5
#
#  Blue things
#
    fill_rgb 0.0 0.0 1.0
    circle_fill 0.5 9.5 0.25
    circle_fill 9.5 3.5 0.25
    circle_fill 5.5 7.5 0.25
    circle_fill 3.5 5.5 0.25
    circle_fill 5.5 9.5 0.25
    circle_fill 8.5 1.5 0.25
    circle_fill 8.5 9.5 0.25
    circle_fill 3.5 9.5 0.25
#
#  Label the plot.
#
    font_size 0.35
    line_rgb 0.0 0.0 0.0

    moveto 0.5 10.0
    label 1
    moveto 9.5 4.0
    label 2
    moveto 6.0 8.0
    label 3
    moveto 4.0 5.5
    label 4
    moveto 5.5 10.0
    label 5
    moveto 9.0 1.5
    label 6
    moveto 9.0 10.0
    label 7
    moveto 3.5 10.0
    label 8
  endpage

endfile
