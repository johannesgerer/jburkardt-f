# prototein.plot  02 August 1999
#
file prototein.ps

  space -1.0 -1.0 5.0 5.0

  page

    line_rgb 0.9 0.9 0.9

    grid -0.0 -0.0 5.0 5.0 6 6

    line_thick 3

    line_rgb 0.0 1.0 0.0

    moveto  1.0  2.0
    drawto  2.0  2.0
    drawto  2.0  1.0
    drawto  2.0  0.0
    drawto  1.0  0.0
    drawto  1.0  1.0
    drawto  0.0  1.0
    drawto  0.0  2.0
    drawto  0.0  3.0
    drawto  1.0  3.0
    drawto  2.0  3.0
    drawto  2.0  4.0
    drawto  3.0  4.0
    drawto  3.0  3.0
    drawto  4.0  3.0
    drawto  4.0  2.0
    drawto  3.0  2.0
    drawto  3.0  1.0
    drawto  4.0  1.0
    drawto  4.0  0.0
    drawto  3.0  0.0
 
    fill_rgb 0.0 0.0 1.0

    circle_fill  1.0  2.0 0.1
    circle_fill  2.0  2.0 0.1
    circle_fill  2.0  1.0 0.1
    circle_fill  2.0  0.0 0.1
    circle_fill  1.0  0.0 0.1
    circle_fill  1.0  1.0 0.1
    circle_fill  0.0  1.0 0.1
    circle_fill  0.0  2.0 0.1
    circle_fill  0.0  3.0 0.1
    circle_fill  1.0  3.0 0.1
    circle_fill  2.0  3.0 0.1
    circle_fill  2.0  4.0 0.1
    circle_fill  3.0  4.0 0.1
    circle_fill  3.0  3.0 0.1
    circle_fill  4.0  3.0 0.1
    circle_fill  4.0  2.0 0.1
    circle_fill  3.0  2.0 0.1
    circle_fill  3.0  1.0 0.1
    circle_fill  4.0  1.0 0.1
    circle_fill  4.0  0.0 0.1
    circle_fill  3.0  0.0 0.1
 
    moveto 0.0 4.5
    line_rgb 0.0 0.0 0.0
    line_thick 1
    font_size 0.25
    label A prototein with 21 atoms.

  endpage

endfile
