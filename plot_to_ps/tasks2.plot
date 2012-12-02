# tasks2.plot  28 June 2008
#
file tasks2.ps

  space -1.0 -1.0 21.0 11.0

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 20.0 10.0 21 11
#
#  Blue things
#
    fill_rgb 0.0 0.0 1.0

    circle_fill   0.0  5.0 0.25

    circle_fill  10.0 10.0 0.25
    circle_fill  10.0  8.0 0.25
    circle_fill  10.0  6.0 0.25

    circle_fill  10.0  0.0 0.25

    circle_fill  20.0  5.0 0.25
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0

    font_size 0.20

    moveto -0.10 5.5
    label Start

    moveto 9.9 10.5
    label 2
    moveto 9.9 8.5
    label 3
    moveto 9.9 6.5
    label 4

    circle_fill  10.0  4.0 0.25
    circle_fill  10.0  3.0 0.25
    circle_fill  10.0  2.0 0.25

    moveto 9.9 0.5
    label 100

    moveto 19.9 5.5
    label End
#
#  Red things
#
    fill_rgb 1.0 0.0 0.0

    arrow  0.5 5.0  9.5 10.0
    arrow  0.5 5.0  9.5  8.0
    arrow  0.5 5.0  9.5  6.0

    arrow  0.5 5.0  9.5  0.0

    arrow  10.5 10.0  19.75 5.0
    arrow  10.5  8.0  19.75 5.0
    arrow  10.5  6.0  19.75 5.0

    arrow  10.5  0.0  19.75 5.0

  endpage

endfile
