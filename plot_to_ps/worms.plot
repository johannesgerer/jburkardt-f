# worms.plot  02 August 1999
#
file worms.ps

  space -0.5 -0.5 10.5 6.5

  page

    line_rgb 0.9 0.9 0.9

    grid  0.0  0.0 10.0 6.0 11 7

    line_thick 3

    line_rgb 0.5 0.0 0.3

    moveto  0.5  0.5
    drawto  1.5  0.5
    drawto  2.5  0.5
    drawto  3.5  0.5
    drawto  4.5  0.5
    drawto  5.5  0.5
    drawto  6.5  0.5
    drawto  7.5  0.5
    drawto  8.5  0.5
    drawto  9.5  0.5

    fill_rgb 0.9 0.75 0.75

    square_fill  0.5  0.5 0.50
    square_fill  1.5  0.5 0.4
    square_fill  2.5  0.5 0.4
    square_fill  3.5  0.5 0.4
    square_fill  4.5  0.5 0.4
    square_fill  5.5  0.5 0.4
    square_fill  6.5  0.5 0.4
    square_fill  7.5  0.5 0.4
    square_fill  8.5  0.5 0.4
    square_fill  9.5  0.5 0.3

    moveto  0.5  3.5
    drawto  0.5  4.5
    drawto  1.5  4.5
    drawto  1.5  3.5
    drawto  2.5  3.5
    drawto  2.5  2.5
    drawto  3.5  2.5
    drawto  3.5  3.5
    drawto  4.5  3.5
    drawto  4.5  4.5

    fill_rgb 0.95 0.75 0.75

    square_fill  0.5  3.5 0.50
    square_fill  0.5  4.5 0.4
    square_fill  1.5  4.5 0.4
    square_fill  1.5  3.5 0.4
    square_fill  2.5  3.5 0.4
    square_fill  2.5  2.5 0.4 
    square_fill  3.5  2.5 0.4
    square_fill  3.5  3.5 0.4
    square_fill  4.5  3.5 0.4
    square_fill  4.5  4.5 0.3

    moveto  7.5  4.5
    drawto  7.5  3.5
    drawto  6.5  3.5
    drawto  6.5  4.5
    drawto  6.5  5.5
    drawto  7.5  5.5
    drawto  8.5  5.5
    drawto  8.5  4.5
    drawto  8.5  3.5
    drawto  8.5  2.5

    fill_rgb 0.95 0.75 0.8

    square_fill  7.5  4.5 0.50
    square_fill  7.5  3.5 0.4
    square_fill  6.5  3.5 0.4
    square_fill  6.5  4.5 0.4
    square_fill  6.5  5.5 0.4
    square_fill  7.5  5.5 0.4
    square_fill  8.5  5.5 0.4
    square_fill  8.5  4.5 0.4
    square_fill  8.5  3.5 0.4
    square_fill  8.5  2.5 0.3

    moveto 0.0 5.5
    line_rgb 0.0 0.0 0.0
    line_thick 1
    font_size 0.25
    label 3 model worm houses with 10 rooms.

  endpage

endfile
