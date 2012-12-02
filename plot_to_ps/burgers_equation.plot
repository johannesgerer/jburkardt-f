# burgers_equation.plot  09 April 2012
#
#  11 nodes in X direction, 6 in Y direction.
#  Bottom row is solid red circles.
#  Left and right columns are solid green circles.
#  Remaining are light gray.
#
file heat_equation.ps

  space -1.0 -1.0 12.0 7.0

  page

    fill_rgb 0.0 0.0 0.0

    grid 0.0 0.0 11.0 6.0 12 7
#
#  Red Circles
#
    fill_rgb 1.0 0.0 0.0

    circle_fill  0.0 0.0 0.25
    circle_fill  1.0 0.0 0.25
    circle_fill  2.0 0.0 0.25
    circle_fill  3.0 0.0 0.25
    circle_fill  4.0 0.0 0.25
    circle_fill  5.0 0.0 0.25
    circle_fill  6.0 0.0 0.25
    circle_fill  7.0 0.0 0.25
    circle_fill  8.0 0.0 0.25
    circle_fill  9.0 0.0 0.25
    circle_fill 10.0 0.0 0.25
    circle_fill 11.0 0.0 0.25
#
#  Green circles.
#
    fill_rgb 0.0 1.0 0.0

    circle_fill 11.0 1.0 0.25
    circle_fill 11.0 2.0 0.25
    circle_fill 11.0 3.0 0.25
    circle_fill 11.0 4.0 0.25
    circle_fill 11.0 5.0 0.25
    circle_fill 11.0 6.0 0.25
#
#  Gray circles.
#
    fill_rgb 0.8 0.8 0.8

    circle_fill  0.0 1.0 0.25
    circle_fill  0.0 2.0 0.25
    circle_fill  0.0 3.0 0.25
    circle_fill  0.0 4.0 0.25
    circle_fill  0.0 5.0 0.25
    circle_fill  0.0 6.0 0.25

    circle_fill  1.0 1.0 0.25
    circle_fill  2.0 1.0 0.25
    circle_fill  3.0 1.0 0.25
    circle_fill  4.0 1.0 0.25
    circle_fill  5.0 1.0 0.25
    circle_fill  6.0 1.0 0.25
    circle_fill  7.0 1.0 0.25
    circle_fill  8.0 1.0 0.25
    circle_fill  9.0 1.0 0.25
    circle_fill 10.0 1.0 0.25

    circle_fill  1.0 2.0 0.25
    circle_fill  2.0 2.0 0.25
    circle_fill  3.0 2.0 0.25
    circle_fill  4.0 2.0 0.25
    circle_fill  5.0 2.0 0.25
    circle_fill  6.0 2.0 0.25
    circle_fill  7.0 2.0 0.25
    circle_fill  8.0 2.0 0.25
    circle_fill  9.0 2.0 0.25
    circle_fill 10.0 2.0 0.25

    circle_fill  1.0 3.0 0.25
    circle_fill  2.0 3.0 0.25
    circle_fill  3.0 3.0 0.25
    circle_fill  4.0 3.0 0.25
    circle_fill  5.0 3.0 0.25
    circle_fill  6.0 3.0 0.25
    circle_fill  7.0 3.0 0.25
    circle_fill  8.0 3.0 0.25
    circle_fill  9.0 3.0 0.25
    circle_fill 10.0 3.0 0.25

    circle_fill  1.0 4.0 0.25
    circle_fill  2.0 4.0 0.25
    circle_fill  3.0 4.0 0.25
    circle_fill  4.0 4.0 0.25
    circle_fill  5.0 4.0 0.25
    circle_fill  6.0 4.0 0.25
    circle_fill  7.0 4.0 0.25
    circle_fill  8.0 4.0 0.25
    circle_fill  9.0 4.0 0.25
    circle_fill 10.0 4.0 0.25

    circle_fill  1.0 5.0 0.25
    circle_fill  2.0 5.0 0.25
    circle_fill  3.0 5.0 0.25
    circle_fill  4.0 5.0 0.25
    circle_fill  5.0 5.0 0.25
    circle_fill  6.0 5.0 0.25
    circle_fill  7.0 5.0 0.25
    circle_fill  8.0 5.0 0.25
    circle_fill  9.0 5.0 0.25
    circle_fill 10.0 5.0 0.25

    circle_fill  1.0 6.0 0.25
    circle_fill  2.0 6.0 0.25
    circle_fill  3.0 6.0 0.25
    circle_fill  4.0 6.0 0.25
    circle_fill  5.0 6.0 0.25
    circle_fill  6.0 6.0 0.25
    circle_fill  7.0 6.0 0.25
    circle_fill  8.0 6.0 0.25
    circle_fill  9.0 6.0 0.25
    circle_fill 10.0 6.0 0.25
#
#  Pure green is not that readable!
#
    fill_rgb 0.0 0.5 0.4
    font_size 0.25
    moveto 11.75 1.0
    label_slant 90 Periodic boundary condition

    fill_rgb 1.0 0.0 0.0
    font_size 0.25
    moveto 4.0 -0.75
    label_slant 0 Initial conditions


  endpage

endfile
