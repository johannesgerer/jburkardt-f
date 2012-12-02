# sixbysix_4.plot  27 February 2010
#
file sixbysix_4.ps

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

    fill_rgb 1.0 0.0 0.0

    circle_fill  1.0 3.0 0.15
    circle_fill  2.0 3.0 0.15
    circle_fill  3.0 3.0 0.15
    circle_fill  4.0 3.0 0.15

    fill_rgb 0.9 0.9 0.9

    circle_fill  1.0 4.0 0.15
    circle_fill  2.0 4.0 0.15
    circle_fill  3.0 4.0 0.15
    circle_fill  4.0 4.0 0.15

    circle_fill  1.0 5.0 0.15
    circle_fill  2.0 5.0 0.15
    circle_fill  3.0 5.0 0.15
    circle_fill  4.0 5.0 0.15
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
    circle_fill  5.0 5.0 0.15
#
#  Intermediate X values.
#
    fill_rgb 0.0 0.7 0.3

    circle_fill  0.5 0.0 0.10
    circle_fill  1.5 0.0 0.10
    circle_fill  2.5 0.0 0.10
    circle_fill  3.5 0.0 0.10
    circle_fill  4.5 0.0 0.10

    circle_fill  0.5 1.0 0.10
    circle_fill  1.5 1.0 0.10
    circle_fill  2.5 1.0 0.10
    circle_fill  3.5 1.0 0.10
    circle_fill  4.5 1.0 0.10

    circle_fill  0.5 2.0 0.10
    circle_fill  1.5 2.0 0.10
    circle_fill  2.5 2.0 0.10
    circle_fill  3.5 2.0 0.10
    circle_fill  4.5 2.0 0.10

    fill_rgb 0.9 0.9 0.9

    circle_fill  0.5 3.0 0.10
    circle_fill  1.5 3.0 0.10
    circle_fill  2.5 3.0 0.10
    circle_fill  3.5 3.0 0.10
    circle_fill  4.5 3.0 0.10

    circle_fill  0.5 4.0 0.10
    circle_fill  1.5 4.0 0.10
    circle_fill  2.5 4.0 0.10
    circle_fill  3.5 4.0 0.10
    circle_fill  4.5 4.0 0.10
#
#  Arrow
#
    linewidth 3

    fill_rgb 1.0 0.0 0.0
    arrow 1.0 2.0 1.0 2.75
    fill_rgb 1.0 0.0 0.0
    arrow 2.0 2.0 2.0 2.75
    fill_rgb 1.0 0.0 0.0
    arrow 3.0 2.0 3.0 2.75
    fill_rgb 1.0 0.0 0.0
    arrow 4.0 2.0 4.0 2.75
#
#  Labels
#
    fill_rgb 0.0 0.0 0.0
    font_size 0.30
    moveto 5.2 -.15
    label T = 0
    moveto 5.2 0.85
    label T = 1
    moveto 5.2 1.85
    label T = 2
    moveto 5.2 2.85
    label T = 3
    moveto 5.2 3.85
    label T = 4
    moveto 5.2 4.85
    label T = 5
    moveto 1.5 -0.50
    label Initial Conditions
    moveto 1.0  5.50
    label Marching Forward in Time
 
  endpage

endfile
