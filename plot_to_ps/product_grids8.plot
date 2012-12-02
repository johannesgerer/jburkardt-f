# product_grids4.plot  21 April 2010
#
file product_grids4.eps

  space -0.33 -0.33 6.33 6.33

  page

#
#  Color the 6x6 boxes.
#
    fill_rgb 1.000  0.750  0.750
    fill_rgb 1.000  0.875  0.375

    polygon_fill
       0.0  0.0
       6.0  0.0
       6.0  6.0
       0.0  6.0
       0.0  0.0
    endpolygon
#
#  Lay down a pale grid.
#
    line_rgb 0.6 0.6 0.6

    grid  1.0  1.0   5.0   5.0  5 5
#
#  Draw the points.
#
    fill_rgb 0.0 0.0 1.0

    circle_fill  1.0  1.0 0.2
    circle_fill  1.0  2.0 0.2
    circle_fill  1.0  3.0 0.2
    circle_fill  1.0  4.0 0.2
    circle_fill  1.0  5.0 0.2
    circle_fill  2.0  1.0 0.2
    circle_fill  2.0  2.0 0.2
    circle_fill  2.0  3.0 0.2
    circle_fill  2.0  4.0 0.2
    circle_fill  2.0  5.0 0.2
    circle_fill  3.0  1.0 0.2
    circle_fill  3.0  2.0 0.2
    circle_fill  3.0  3.0 0.2
    circle_fill  3.0  4.0 0.2
    circle_fill  3.0  5.0 0.2
    circle_fill  4.0  1.0 0.2
    circle_fill  4.0  2.0 0.2
    circle_fill  4.0  3.0 0.2
    circle_fill  4.0  4.0 0.2
    circle_fill  4.0  5.0 0.2
    circle_fill  5.0  1.0 0.2
    circle_fill  5.0  2.0 0.2
    circle_fill  5.0  3.0 0.2
    circle_fill  5.0  4.0 0.2
    circle_fill  5.0  5.0 0.2

  endpage

endfile
