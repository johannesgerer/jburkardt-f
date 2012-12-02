# tasks3.plot  29 June 2008
#
file tasks3.ps

  space -1.0 -1.0 21.0 11.0

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 20.0 10.0 21 11
#
#  Gray things
#
    fill_rgb 0.5 0.5 0.5

    circle_fill  8.0 10.0 0.25

    circle_fill 10.0 10.0 0.25
    circle_fill 10.0  8.0 0.25

    circle_fill 12.0 10.0 0.25
    circle_fill 12.0  8.0 0.25
    circle_fill 12.0  6.0 0.25

    circle_fill 14.0 10.0 0.25
    circle_fill 14.0  8.0 0.25
    circle_fill 14.0  6.0 0.25
    circle_fill 14.0  4.0 0.25

    circle_fill 16.0 10.0 0.25
    circle_fill 16.0  8.0 0.25
    circle_fill 16.0  6.0 0.25
    circle_fill 16.0  4.0 0.25
    circle_fill 16.0  2.0 0.25
#
#  Red things
#
    fill_rgb 1.0 0.0 0.0

    circle_fill  6.0 10.0 0.25

    circle_fill  8.0  8.0 0.25

    circle_fill 10.0  6.0 0.25

    circle_fill 12.0  4.0 0.25

    circle_fill 14.0  2.0 0.25

    circle_fill 16.0  0.0 0.25
#
#  Blue things
#
    fill_rgb 0.0 0.0 1.0

    circle_fill  6.0  8.0 0.25
    circle_fill  6.0  6.0 0.25
    circle_fill  6.0  4.0 0.25
    circle_fill  6.0  2.0 0.25
    circle_fill  6.0  0.0 0.25

    circle_fill  8.0  6.0 0.25
    circle_fill  8.0  4.0 0.25
    circle_fill  8.0  2.0 0.25
    circle_fill  8.0  0.0 0.25

    circle_fill 10.0  4.0 0.25
    circle_fill 10.0  2.0 0.25
    circle_fill 10.0  0.0 0.25

    circle_fill 12.0  2.0 0.25
    circle_fill 12.0  0.0 0.25

    circle_fill 14.0  0.0 0.25
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0

    font_size 0.20

    moveto 5.9 10.5
    label 1

    moveto  5.9  8.5
    label 2
    moveto  5.9  6.5
    label 2
    moveto  5.9  4.5
    label 2
    moveto  5.9  2.5
    label 2
    moveto  5.9  0.5
    label 2

    moveto 7.9 8.5
    label 3

    moveto  7.9  6.5
    label 4
    moveto  7.9  4.5
    label 4
    moveto  7.9  2.5
    label 4
    moveto  7.9  0.5
    label 4

    moveto 9.9 6.5
    label 5

    moveto  9.9  4.5
    label 6
    moveto  9.9  2.5
    label 6
    moveto  9.9  0.5
    label 6

    moveto 11.9 4.5
    label 7

    moveto  11.9  2.5
    label 8
    moveto  11.9  0.5
    label 8

    moveto 13.9 2.5
    label 9

    moveto 13.8 0.5
    label 10

    moveto 15.8 0.5
    label 11

  endpage

endfile
