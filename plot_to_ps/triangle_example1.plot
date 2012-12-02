# triangle_example1.plot  07 June 2009
#
file triangle.ps

  space -1.0 -1.0 11.0 11.0

  page

    line_rgb 0.3 0.3 0.3

    grid 0.0 0.0 10.0 10.0 11 11

    fill_rgb 0.35 0.95 0.9

    triangle_fill 4.0 1.0 8.0 3.0 0.0 9.0

    line_thick 3

    moveto   4.0  1.0
    drawto  8.0  3.0
    drawto   0.0  9.0
    drawto   4.0  1.0

    line_thick 5

    fill_rgb 0.5 0.7 0.5

    circle_fill  4.0 1.0 0.25
    circle_fill 8.0 3.0 0.25
    circle_fill 0.0  9.0 0.25

    font_size 0.35
    fill_rgb 0.0 0.0 0.0

    moveto  4.5  0.5
    label A = {4,1}
    moveto  6.5.0  4.0
    label B = {8,3}
    moveto  0.0  9.5
    label C = {0,9}

  endpage

endfile
