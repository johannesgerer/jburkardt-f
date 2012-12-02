# spiral.plot
#
file spiral.eps

  space -1.0 3.0 6.0 10.0

  page

    line_rgb 0.2 0.2 0.2
    grid -1.0 3.0 6.0 10.0 8 8

    line_rgb 1.0 0.0 0.0

    font_size 0.25

    moveto  0.25  8.5
    label 10

    moveto  1.25  8.5
    label 11

    moveto  2.25  8.5
    label 12

    moveto  3.25  8.5
    label 13

    moveto  0.5  7.5
    label 9

    moveto  1.5  7.5
    label 2

    moveto  2.5  7.5
    label 3

    moveto  3.25  7.5
    label 14

    moveto  0.5  6.5
    label 8

    moveto  1.5  6.5
    label 1

    point 1.25 6.25

    moveto  2.5  6.5
    label 4

    moveto  3.25  6.5
    label 15

    moveto  0.5  5.5
    label 7

    moveto  1.5  5.5
    label 6

    moveto  2.5  5.5
    label 5

    line_rgb 0.7 0.7 0.7

    moveto  1.4  6.5
    drawto  1.4  7.4
    drawto  2.4  7.4
    drawto  2.4  6.5
    drawto  2.4  5.4
    drawto  1.5  5.4
    drawto  0.4  5.4
    drawto  0.4  6.5
    drawto  0.4  7.5
    drawto  0.4   8.4
    drawto  1.25  8.4
    drawto  2.25  8.4
    drawto  3.2   8.4
    drawto  3.2   7.5
    drawto  3.2   6.5

    line_rgb 0.0 0.0 0.0

    point 1.25 6.25
    point 3.75 6.10
    point  2.4  3.2
    point -0.3  7.1
    point  4.1  8.3
    point -0.2 5.3
    point  5.4 4.3

  endpage

endfile
