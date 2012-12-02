# museum2.plot  24 January 2011
#
file museum2.ps

  space -0.50 -0.50 12.50 6.50

  page

    line_rgb 0.5 0.5 0.5

    grid 0.0 0.0 12.0 6.0 7 4

    fill_rgb 0.94 0.94 0.96
#
#  Fill room with light gray.
#
    polygon_fill
      0.0 0.0
     12.0 0.0
     12.0 6.0
      0.0 6.0
    end
#
#  Draw room centers
#
    line_thick 3

    line_rgb 0.0 0.0 1.0
    circle  1.0 1.0 0.50
    circle  3.0 1.0 0.50
    circle  5.0 1.0 0.50
    circle  7.0 1.0 0.50
    circle  9.0 1.0 0.50
    circle 11.0 1.0 0.50

    circle  1.0 3.0 0.50
    circle  3.0 3.0 0.50
    circle  5.0 3.0 0.50
    circle  7.0 3.0 0.50
    circle  9.0 3.0 0.50
    circle 11.0 3.0 0.50

    circle  1.0 5.0 0.50
    circle  3.0 5.0 0.50
    circle  5.0 5.0 0.50
    circle  7.0 5.0 0.50
    circle  9.0 5.0 0.50
    circle 11.0 5.0 0.50
#
#  Edges.
#
    line_thick 4
    line_rgb 1.0 0.0 0.0

    moveto 1.0  1.5
    drawto 1.0  2.5
    moveto 1.0  3.5
    drawto 1.0  4.5

    moveto 3.0  1.5
    drawto 3.0  2.5
    moveto 3.0  3.5
    drawto 3.0  4.5

    moveto 5.0  1.5
    drawto 5.0  2.5


    moveto 9.0 -0.5
    drawto 9.0  0.5
    moveto 9.0  1.5
    drawto 9.0  2.5
    moveto 9.0  3.5
    drawto 9.0  4.5

    moveto 11.0  1.5
    drawto 11.0  2.5

    moveto  3.5 5.0
    drawto  4.5 5.0
    moveto  7.5 5.0
    drawto  8.5 5.0
    moveto  9.5 5.0
    drawto 10.5 5.0

    moveto  3.5 3.0
    drawto  4.5 3.0
    moveto  5.5 3.0
    drawto  6.5 3.0
    moveto  9.5 3.0
    drawto 10.5 3.0

    moveto  1.5 1.0
    drawto  2.5 1.0
    moveto  5.5 1.0
    drawto  6.5 1.0
    moveto  7.5 1.0
    drawto  8.5 1.0
#
#  Draw room boundaries.
#
    line_thick 4

    line_rgb 0.1 0.6 0.1

    moveto  8.0 0.0
    drawto  0.0 0.0
    drawto  0.0 6.0
    drawto 12.0 6.0
    drawto 12.0 0.0
    drawto 10.0 0.0
    drawto 10.0 2.0

    moveto 10.0 4.0
    drawto 12.0 4.0

    moveto 2.0 2.0
    drawto 2.0 6.0

    moveto 4.0 0.0
    drawto 4.0 2.0

    moveto 4.0 4.0
    drawto 6.0 4.0
    drawto 6.0 6.0

    moveto 6.0 4.0
    drawto 8.0 4.0
    drawto 8.0 2.0
    drawto 6.0 2.0
#
#  Label the plot.
#
    font_size 0.40
    line_rgb 0.0 0.0 0.0

    moveto  0.8 0.8
    label   M
    moveto  2.8 0.8
    label   N
    moveto  4.8 0.8
    label   O
    moveto  6.8 0.8
    label   P
    moveto  8.8 0.8
    label   Q
    moveto 10.8 0.8
    label   R

    moveto  0.8 2.8
    label   G
    moveto  2.8 2.8
    label   H
    moveto  4.8 2.8
    label   I
    moveto  6.8 2.8
    label   J
    moveto  8.8 2.8
    label   K
    moveto 10.8 2.8
    label   L

    moveto  0.8 4.8
    label   A
    moveto  2.8 4.8
    label   B
    moveto  4.8 4.8
    label   C
    moveto  6.8 4.8
    label   D
    moveto  8.8 4.8
    label   E
    moveto 10.8 4.8
    label   F

  endpage

endfile
