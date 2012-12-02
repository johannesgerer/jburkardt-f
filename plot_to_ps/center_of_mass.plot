# center_of_mass.plot  23 September 2009
#
file center_of_mass.ps

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
    circle_fill  2.8  7.1 0.60
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0

    font_size 0.20

    moveto 1.0 10.5
    label 50
    moveto 5.0 4.5
    label 2
    moveto 3.0 1.5
    label 20
    moveto 7.0 7.5
    label 1
    moveto 12.0 5.5
    label 5
    moveto 20.0 3.5
    label 2
    moveto 2.8 7.6
    label C = ( 2.8, 7.1)
    font_size 0.30
    moveto 8.5 10.0
    label Center of Mass
#
#  Red things
#
    fill_rgb 1.0 0.0 0.0

    line   1.0 10.0 2.8  7.1
    line   5.0  4.0 2.8  7.1
    line   3.0  1.0 2.8  7.1
    line   7.0  7.0 2.8  7.1
    line  12.0  5.0 2.8  7.1
    line  20.0  3.0 2.8  7.1

  endpage

endfile
