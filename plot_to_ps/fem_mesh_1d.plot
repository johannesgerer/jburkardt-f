# fem_mesh_1d.plot  18 October 2011
#
file fem_mesh_1d.ps

  space -1.0 -1.0 21.0 11.0

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 20.0 10.0 21 11

    linewidth 3
#
#  Blue things
#
    fill_rgb 0.0 0.0 1.0

    circle_fill  0.0 9.0 0.25
    circle_fill 20.0 9.0 0.25

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

    circle_fill  0.0 1.0 0.25
    circle_fill  2.0 1.0 0.25
    circle_fill  4.0 1.0 0.25
    circle_fill  6.0 1.0 0.25
    circle_fill  8.0 1.0 0.25
    circle_fill 10.0 1.0 0.25
    circle_fill 12.0 1.0 0.25
    circle_fill 14.0 1.0 0.25
    circle_fill 16.0 1.0 0.25
    circle_fill 18.0 1.0 0.25
    circle_fill 20.0 1.0 0.25
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0

    font_size 0.20

    moveto -0.10 10.0
    label A
    moveto 19.8 10.0
    label B

    moveto -0.10 6.0
    label N1
    moveto 1.9 6.0
    label N2
    moveto 3.9 6.0
    label N3
    moveto 5.9 6.0
    label N4
    moveto 7.9 6.0
    label N5
    moveto 9.9 6.0
    label N6
    moveto 11.9 6.0
    label N7
    moveto 13.9 6.0
    label N8
    moveto 15.9 6.0
    label N9
    moveto 17.8 6.0
    label N10
    moveto 19.8 6.0
    label N11

    moveto 0.8 2.0
    label E1
    moveto 2.8 2.0
    label E2
    moveto 4.8 2.0
    label E3
    moveto 6.8 2.0
    label E4
    moveto 8.8 2.0
    label E5
    moveto 10.8 2.0
    label E6
    moveto 12.8 2.0
    label E7
    moveto 14.8 2.0
    label E8
    moveto 16.8 2.0
    label E9
    moveto 18.8 2.0
    label E10
#
#  Red things
#
    fill_rgb 1.0 0.0 0.0

    line 0.25 9.0  19.75 9.0

    line  0.25 1.0  1.75 1.0
    line  2.25 1.0  3.75 1.0
    line  4.25 1.0  5.75 1.0
    line  6.25 1.0  7.75 1.0
    line  8.25 1.0  9.75 1.0
    line 10.25 1.0 11.75 1.0
    line 12.25 1.0 13.75 1.0
    line 14.25 1.0 15.75 1.0
    line 16.25 1.0 17.75 1.0
    line 18.25 1.0 19.75 1.0



  endpage

endfile
