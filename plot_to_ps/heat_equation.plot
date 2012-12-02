# heat_equation.plot  16 April 2008
#
file heat_equation.ps

  space -1.0 -1.0 21.0 11.0

  page

    grid 0.0 0.0 20.0 10.0 21 11
#
#  The Red region
#
    fill_rgb 1.0 0.0 0.0

    circle_fill  0.0 0.0 0.25
    circle_fill  1.0 0.0 0.25
    circle_fill  2.0 0.0 0.25
    circle_fill  3.0 0.0 0.25
    circle_fill  4.0 0.0 0.25
    circle_fill  5.0 0.0 0.25
    circle_fill  6.0 0.0 0.25

    circle_fill  0.0 1.0 0.25
    circle_fill  1.0 1.0 0.25
    circle_fill  2.0 1.0 0.25
    circle_fill  3.0 1.0 0.25
    circle_fill  4.0 1.0 0.25
    circle_fill  5.0 1.0 0.25
    circle_fill  6.0 1.0 0.25

    circle_fill  0.0 2.0 0.25
    circle_fill  1.0 2.0 0.25
    circle_fill  2.0 2.0 0.25
    circle_fill  3.0 2.0 0.25
    circle_fill  4.0 2.0 0.25
    circle_fill  5.0 2.0 0.25
    circle_fill  6.0 2.0 0.25

    circle       0.0 3.0 0.25
    circle       1.0 3.0 0.25
    circle       2.0 3.0 0.25
    circle       3.0 3.0 0.25
    circle       4.0 3.0 0.25
    circle       5.0 3.0 0.25
    circle       6.0 3.0 0.25
#
#  A typical stencil
#
    circle_fill  2.0 6.0 0.25
    circle_fill  3.0 6.0 0.25
    circle_fill  4.0 6.0 0.25

    circle       3.0 7.0 0.25
#
#  The Green region
#
    fill_rgb 0.0 1.0 0.0

    circle_fill  7.0 0.0 0.25
    circle_fill  8.0 0.0 0.25
    circle_fill  9.0 0.0 0.25
    circle_fill 10.0 0.0 0.25
    circle_fill 11.0 0.0 0.25
    circle_fill 12.0 0.0 0.25
    circle_fill 13.0 0.0 0.25

    circle_fill  7.0 1.0 0.25
    circle_fill  8.0 1.0 0.25
    circle_fill  9.0 1.0 0.25
    circle_fill 10.0 1.0 0.25
    circle_fill 11.0 1.0 0.25
    circle_fill 12.0 1.0 0.25
    circle_fill 13.0 1.0 0.25

    circle_fill  7.0 2.0 0.25
    circle_fill  8.0 2.0 0.25
    circle_fill  9.0 2.0 0.25
    circle_fill 10.0 2.0 0.25
    circle_fill 11.0 2.0 0.25
    circle_fill 12.0 2.0 0.25
    circle_fill 13.0 2.0 0.25

    circle       7.0 3.0 0.25
    circle       8.0 3.0 0.25
    circle       9.0 3.0 0.25
    circle      10.0 3.0 0.25
    circle      11.0 3.0 0.25
    circle      12.0 3.0 0.25
    circle      13.0 3.0 0.25
#
#  A sample stencil that needs one value from the right.
#
    circle_fill 12.0 6.0 0.25
    circle_fill 13.0 6.0 0.25

    circle      13.0 7.0 0.25
#
#  The Blue Region
#
    fill_rgb 0.0 0.0 1.0

    circle_fill 14.0 0.0 0.25
    circle_fill 15.0 0.0 0.25
    circle_fill 16.0 0.0 0.25
    circle_fill 17.0 0.0 0.25
    circle_fill 18.0 0.0 0.25
    circle_fill 19.0 0.0 0.25
    circle_fill 20.0 0.0 0.25

    circle_fill 14.0 1.0 0.25
    circle_fill 15.0 1.0 0.25
    circle_fill 16.0 1.0 0.25
    circle_fill 17.0 1.0 0.25
    circle_fill 18.0 1.0 0.25
    circle_fill 19.0 1.0 0.25
    circle_fill 20.0 1.0 0.25

    circle_fill 14.0 2.0 0.25
    circle_fill 15.0 2.0 0.25
    circle_fill 16.0 2.0 0.25
    circle_fill 17.0 2.0 0.25
    circle_fill 18.0 2.0 0.25
    circle_fill 19.0 2.0 0.25
    circle_fill 20.0 2.0 0.25

    circle      14.0 3.0 0.25
    circle      15.0 3.0 0.25
    circle      16.0 3.0 0.25
    circle      17.0 3.0 0.25
    circle      18.0 3.0 0.25
    circle      19.0 3.0 0.25
    circle      20.0 3.0 0.25
#
#  One circle that interacts with the interior set.
#
    circle_fill 14.0 6.0 0.25
#
#  A sample stencil at the boundary.
#
    circle_fill 18.0 6.0 0.25
    circle_fill 19.0 6.0 0.25
    circle_fill 20.0 6.0 0.25

    circle      19.0 7.0 0.25

  endpage

endfile
