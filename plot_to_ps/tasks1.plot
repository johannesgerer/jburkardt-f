# tasks1.plot  28 June 2008
#
file tasks1.ps

  space -1.0 -1.0 21.0 11.0

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 20.0 10.0 21 11
#
#  Blue things
#
    fill_rgb 0.0 0.0 1.0
    circle_fill  0.0 5.0 0.25
    circle_fill  2.0 5.0 0.25
    circle_fill  4.0 5.0 0.25
    circle_fill  6.0 5.0 0.25
    circle_fill  8.0 5.0 0.25
    circle_fill 10.0 5.0 0.25
    circle_fill 12.0 5.0 0.25
    circle_fill 14.0 5.0 0.25
    circle_fill 16.0 5.0 0.25
    circle_fill 18.0 5.0 0.25
    circle_fill 20.0 5.0 0.25
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0

    font_size 0.20

    moveto -0.10 6.0
    label 1
    moveto 1.9 6.0
    label 2
    moveto 3.9 6.0
    label 3
    moveto 5.9 6.0
    label 4
    moveto 7.9 6.0
    label 5
    moveto 9.9 6.0
    label 6
    moveto 11.9 6.0
    label 7
    moveto 13.9 6.0
    label 8
    moveto 15.9 6.0
    label 9
    moveto 17.8 6.0
    label 10
    moveto 19.8 6.0
    label 11
#
#  Red things
#
    fill_rgb 1.0 0.0 0.0

    arrow  0.25 5.0  1.75 5.0
    arrow  2.25 5.0  3.75 5.0
    arrow  4.25 5.0  5.75 5.0
    arrow  6.25 5.0  7.75 5.0
    arrow  8.25 5.0  9.75 5.0
    arrow 10.25 5.0 11.75 5.0
    arrow 12.25 5.0 13.75 5.0
    arrow 14.25 5.0 15.75 5.0
    arrow 16.25 5.0 17.75 5.0
    arrow 18.25 5.0 19.75 5.0
  endpage

endfile
