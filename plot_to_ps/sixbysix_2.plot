# sixbysix_2.plot  17 February 2010
#
file sixbysix_2.ps

  space -1.0 -1.0 6.0 6.0

  page

    line_rgb 0.6 0.6 0.6
    grid 0.0 0.0 5.0 5.0 6 6

    line_rgb 0.0 0.0 0.0
    linewidth 2
    grid 0.5 0.0 4.5 5.0 5 0
    grid 0.0 0.5 5.0 4.5 0 5

    linewidth 1
#
#  The interior points.
#
    fill_rgb 0.0 0.0 1.0

    circle_fill  1.0 1.0 0.15
    circle_fill  2.0 1.0 0.15
    circle_fill  3.0 1.0 0.15
    circle_fill  4.0 1.0 0.15

    circle_fill  1.0 2.0 0.15
    circle_fill  2.0 2.0 0.15
    circle_fill  3.0 2.0 0.15
    circle_fill  4.0 2.0 0.15

    circle_fill  1.0 3.0 0.15
    circle_fill  2.0 3.0 0.15
    circle_fill  3.0 3.0 0.15
    circle_fill  4.0 3.0 0.15

    circle_fill  1.0 4.0 0.15
    circle_fill  2.0 4.0 0.15
    circle_fill  3.0 4.0 0.15
    circle_fill  4.0 4.0 0.15
#
#  The boundary points.
#
    fill_rgb 0.6 0.6 1.0

    circle_fill  0.0 0.0 0.15
    circle_fill  1.0 0.0 0.15
    circle_fill  2.0 0.0 0.15
    circle_fill  3.0 0.0 0.15
    circle_fill  4.0 0.0 0.15
    circle_fill  5.0 0.0 0.15

    circle_fill  0.0 1.0 0.15
    circle_fill  5.0 1.0 0.15

    circle_fill  0.0 2.0 0.15
    circle_fill  5.0 2.0 0.15

    circle_fill  0.0 3.0 0.15
    circle_fill  5.0 3.0 0.15

    circle_fill  0.0 4.0 0.15
    circle_fill  5.0 4.0 0.15

    circle_fill  0.0 5.0 0.15
    circle_fill  1.0 5.0 0.15
    circle_fill  2.0 5.0 0.15
    circle_fill  3.0 5.0 0.15
    circle_fill  4.0 5.0 0.15
    circle_fill  5.0 5.0 0.15

  endpage

endfile
