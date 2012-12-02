# horizon.plot  19 September 2002
#
#  Illustrates the problem of determining the distance from an observer
#  of a given height to the horizon.
#
file horizon.ps

  space 0.0 0.0 100.0 100.0

  page

    line_thick 3
 
    line_rgb 0.0 0.0 1.0

    arc   50.0 30.0 40.0 10.0 120.0

    point 50.0 30.0
    point 50.0 70.0
    point 50.0 80.0
    point 74.0 62.0

    line_rgb 1.0 0.0 0.0

    moveto 50.0  30.0
    drawto 50.0  80.0
    drawto 74.0  62.0
    drawto 50.0  30.0

    font_size 0.25
    moveto 40.0 60.0
    label Side C

    moveto 65.0 70.0
    label Side A

    moveto 65.0 46.0
    label Side B

    font_size 0.50
    moveto 17.0 10.0
    label The horizon problem

  endpage

endfile
