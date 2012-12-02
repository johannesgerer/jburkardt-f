# chebyshev.plot  08 August 2012
#
file chebyshev.ps

  space -1.1  -0.1  1.1  1.1

  page

    line_thick 3
    line_rgb 0.7 0.7 0.7
    grid -1.0 -0.1 1.0 1.0 21 11

    line_thick 4
    line_rgb 0.0 0.0 1.0
    arc 0.0 0.0 1.0 0.0 180.0

    line_thick 3
    line_rgb 0.5 0.5 1.0

    radius 0.0 0.0 1.0   0

    radius 0.0 0.0 1.0  20
    line 0.94 0.34 0.94 0.0

    radius 0.0 0.0 1.0  40
    line 0.77 0.64 0.77 0.0

    radius 0.0 0.0 1.0  60
    line 0.50 0.87 0.50 0.0

    radius 0.0 0.0 1.0  80
    line 0.17 0.98 0.17 0.0

    radius 0.0 0.0 1.0 100
    line -0.17 0.98 -0.17 0.0

    radius 0.0 0.0 1.0 120
    line -0.50 0.87 -0.50 0.0

    radius 0.0 0.0 1.0 140
    line -0.77 0.64 -0.77 0.0

    radius 0.0 0.0 1.0 160
    line -0.94 0.34 -0.94 0.0

    radius 0.0 0.0 1.0 180

    line_thick 4
    line_rgb 0.0 0.0 1.0
    line -1.0  0.0  +1.0  0.0

    line_rgb 1.0 0.0 0.0
    circle_fill  1.00  0.0  0.04
    circle_fill  0.94  0.0  0.04
    circle_fill  0.77  0.0  0.04
    circle_fill  0.50  0.0  0.04
    circle_fill  0.17  0.0  0.04
    circle_fill -0.17  0.0  0.04
    circle_fill -0.50  0.0  0.04
    circle_fill -0.77  0.0  0.04
    circle_fill -0.94  0.0  0.04
    circle_fill -1.00  0.0  0.04


  endpage

endfile
