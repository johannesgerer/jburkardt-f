# centroid.plot  28 June 2008
#
file centroid.ps

  space -1.0 -1.0 21.0 11.0

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 20.0 10.0 21 11
#
#  Blue things
#
    fill_rgb 0.0 0.0 1.0

    circle_fill  1.0 10.0 0.25
    circle_fill  5.0  4.0 0.25
    circle_fill  3.0  1.0 0.25
    circle_fill  7.0  7.0 0.25
    circle_fill 12.0  5.0 0.25
    circle_fill 20.0  3.0 0.25
#
#  Green things.
#
    fill_rgb 0.0 1.0 0.0
    circle_fill  8.0  5.0 0.60
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0

    font_size 0.20

    moveto 1.0 10.5
    label 1
    moveto 5.0 4.5
    label 1
    moveto 3.0 1.5
    label 1
    moveto 7.0 7.5
    label 1
    moveto 12.0 5.5
    label 1
    moveto 20.0 3.5
    label 1
    moveto 8.0 6.0
    label C = ( 8.0, 5.0 )

    font_size 0.30

    moveto 9 10
    label Centroid
#
#  Red things
#
    fill_rgb 1.0 0.0 0.0

    line   1.0 10.0 8.0 5.0
    line   5.0  4.0 8.0 5.0
    line   3.0  1.0 8.0 5.0
    line   7.0  7.0 8.0 5.0
    line  12.0  5.0 8.0 5.0
    line  20.0  3.0 8.0 5.0

  endpage

endfile
