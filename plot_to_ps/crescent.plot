# crescent.plot  02 August 1999
#
file crescent.ps

  space -6.0 -6.0 6.0 6.0

  page

    fill_rgb 0.0 0.6 0.9
    circle_fill 0.0 0.0 4.0

    fill_rgb 1.0 1.0 1.0
    circle_fill 2.0 0.0 3.0


    line_thick 15
    line_rgb 1.0 0.0 0.0
    arc 0.0 0.0 4.02 0.0 360.0
    line_rgb 0.0 1.0 0.0
    arc 2.0 0.0 2.95 0.0 360.0

    fill_rgb 0.0 0.0 0.0
    line_thick 1
    grid -5.0 -5.0 5.0 5.0 11 11

    fill_rgb 1.0 0.0 0.0
    circle_fill 0.0 0.0 0.25
    fill_rgb 0.0 1.0 0.0
    circle_fill 2.0 0.0 0.25

  endpage

endfile
