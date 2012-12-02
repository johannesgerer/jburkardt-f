# product_grids2.plot  21 April 2010
#
file product_grids.eps

  space -1.0 -1.0 19.0 19.0

  page

#
#  Color the 6x6 boxes.
#
    fill_rgb 0.900  0.900  0.900

    polygon_fill
       0.0  0.0
       6.0  0.0
       6.0  6.0
       0.0  6.0
       0.0  0.0
    endpolygon

    fill_rgb 1.000  1.000  0.000

    polygon_fill
       6.0  0.0
      12.0  0.0
      12.0  6.0
       6.0  6.0
       6.0  0.0
    endpolygon

    polygon_fill
       0.0  6.0
       6.0  6.0
       6.0 12.0
       0.0 12.0
       0.0  6.0
    endpolygon

    fill_rgb 1.000  0.750  0.750

    polygon_fill
      12.0  0.0
      18.0  0.0
      18.0  6.0
      12.0  6.0
      12.0  0.0
    endpolygon

    polygon_fill
       6.0  6.0
      12.0  6.0
      12.0 12.0
       6.0 12.0
       6.0  6.0
    endpolygon

    polygon_fill
       0.0 12.0
       6.0 12.0
       6.0 18.0
       0.0 18.0
       0.0 12.0
    endpolygon

    fill_rgb 0.900  0.900  0.900

    polygon_fill
      12.0  6.0
      18.0  6.0
      18.0 12.0
      12.0 12.0
      12.0  6.0
    endpolygon

    polygon_fill
       6.0 12.0
      12.0 12.0
      12.0 18.0
       6.0 18.0
       6.0 12.0
    endpolygon

    fill_rgb 0.900  0.900  0.900

    polygon_fill
      12.0 12.0
      18.0 12.0
      18.0 18.0
      12.0 18.0
      12.0 12.0
    endpolygon
#
#  Lay down a pale grid.
#
    line_rgb 0.6 0.6 0.6

    grid  1.0  1.0   5.0   5.0  5 5
    grid  7.0  1.0  11.0   5.0  5 5
    grid 13.0  1.0  17.0   5.0  5 5

    grid  1.0  7.0   5.0  11.0  5 5
    grid  7.0  7.0  11.0  11.0  5 5
    grid 13.0  7.0  17.0  11.0  5 5

    grid  1.0 13.0   5.0  17.0  5 5
    grid  7.0 13.0  11.0  17.0  5 5
    grid 13.0 13.0  17.0  17.0  5 5
#
#  Draw the points.
#
    fill_rgb 0.5 0.5 0.5

    circle_fill  3.0  3.0 0.2

    fill_rgb 0.0 0.0 0.0

    circle_fill  7.0  3.0 0.2
    circle_fill  9.0  3.0 0.2
    circle_fill 11.0  3.0 0.2
    circle_fill 13.0  3.0 0.2
    circle_fill 14.0  3.0 0.2
    circle_fill 15.0  3.0 0.2
    circle_fill 16.0  3.0 0.2
    circle_fill 17.0  3.0 0.2

    circle_fill  3.0  7.0 0.2
    circle_fill  7.0  7.0 0.2
    circle_fill  9.0  7.0 0.2
    circle_fill 11.0  7.0 0.2

    fill_rgb 0.5 0.5 0.5

    circle_fill 13.0  7.0 0.2
    circle_fill 14.0  7.0 0.2
    circle_fill 15.0  7.0 0.2
    circle_fill 16.0  7.0 0.2
    circle_fill 17.0  7.0 0.2

    fill_rgb 0.0 0.0 0.0

    circle_fill  3.0  9.0 0.2
    circle_fill  7.0  9.0 0.2
    circle_fill  9.0  9.0 0.2
    circle_fill 11.0  9.0 0.2

    fill_rgb 0.5 0.5 0.5

    circle_fill 13.0  9.0 0.2
    circle_fill 14.0  9.0 0.2
    circle_fill 15.0  9.0 0.2
    circle_fill 16.0  9.0 0.2
    circle_fill 17.0  9.0 0.2


    fill_rgb 0.0 0.0 0.0

    circle_fill  3.0 11.0 0.2
    circle_fill  7.0 11.0 0.2
    circle_fill  9.0 11.0 0.2
    circle_fill 11.0 11.0 0.2

    fill_rgb 0.5 0.5 0.5

    circle_fill 13.0 11.0 0.2
    circle_fill 14.0 11.0 0.2
    circle_fill 15.0 11.0 0.2
    circle_fill 16.0 11.0 0.2
    circle_fill 17.0 11.0 0.2

    fill_rgb 0.0 0.0 0.0
    circle_fill  3.0 13.0 0.2
    fill_rgb 0.5 0.5 0.5
    circle_fill  7.0 13.0 0.2
    circle_fill  9.0 13.0 0.2
    circle_fill 11.0 13.0 0.2
    circle_fill 13.0 13.0 0.2
    circle_fill 14.0 13.0 0.2
    circle_fill 15.0 13.0 0.2
    circle_fill 16.0 13.0 0.2
    circle_fill 17.0 13.0 0.2

    fill_rgb 0.0 0.0 0.0
    circle_fill  3.0 14.0 0.2
    fill_rgb 0.5 0.5 0.5
    circle_fill  7.0 14.0 0.2
    circle_fill  9.0 14.0 0.2
    circle_fill 11.0 14.0 0.2
    circle_fill 13.0 14.0 0.2
    circle_fill 14.0 14.0 0.2
    circle_fill 15.0 14.0 0.2
    circle_fill 16.0 14.0 0.2
    circle_fill 17.0 14.0 0.2

    fill_rgb 0.0 0.0 0.0
    circle_fill  3.0 15.0 0.2
    fill_rgb 0.5 0.5 0.5
    circle_fill  7.0 15.0 0.2
    circle_fill  9.0 15.0 0.2
    circle_fill 11.0 15.0 0.2
    circle_fill 13.0 15.0 0.2
    circle_fill 14.0 15.0 0.2
    circle_fill 15.0 15.0 0.2
    circle_fill 16.0 15.0 0.2
    circle_fill 17.0 15.0 0.2

    fill_rgb 0.0 0.0 0.0
    circle_fill  3.0 16.0 0.2
    fill_rgb 0.5 0.5 0.5
    circle_fill  7.0 16.0 0.2
    circle_fill  9.0 16.0 0.2
    circle_fill 11.0 16.0 0.2
    circle_fill 13.0 16.0 0.2
    circle_fill 14.0 16.0 0.2
    circle_fill 15.0 16.0 0.2
    circle_fill 16.0 16.0 0.2
    circle_fill 17.0 16.0 0.2

    fill_rgb 0.0 0.0 0.0
    circle_fill  3.0 17.0 0.2
    fill_rgb 0.5 0.5 0.5
    circle_fill  7.0 17.0 0.2
    circle_fill  9.0 17.0 0.2
    circle_fill 11.0 17.0 0.2
    circle_fill 13.0 17.0 0.2
    circle_fill 14.0 17.0 0.2
    circle_fill 15.0 17.0 0.2
    circle_fill 16.0 17.0 0.2
    circle_fill 17.0 17.0 0.2

  endpage

endfile
