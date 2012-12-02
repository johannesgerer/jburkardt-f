# egyptian.plot  02 August 1999
#
file egyptian.ps

  space -10.0 -10.0 10.0 10.0

  page

    line_rgb 1.0 0.0 0.0

    circle 0.0 0.0 9.0

    line_rgb 0.5 0.5 0.5

    grid -9.0 -9.0 9.0 9.0 19 19

    line_rgb 0.0 1.0 0.0

    square 0.0 0.0 8.0

    line_rgb 0.0 0.0 0.0
    font_size 0.25

    moveto -8.0 10.0
    label The Egyptians approximated the area of a circle
    moveto -8.0 9.5
    label by the area of a square whose side was 8/9 the
    moveto -8.0 9.0
    label diameter of the circle.


  endpage

endfile
