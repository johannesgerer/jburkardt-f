# polygon_fill.plot  02 February 2006
#
file polygon_fill.ps

  space -10.0 -10.0 10.0 10.0

  page

    fill_gray 0.9

    polygon_fill
      -9.0 -3.0
      -9.0  3.0
      -3.0  3.0
      -3.0  9.0
       3.0  9.0
       3.0  3.0
       9.0  3.0
       9.0 -3.0
       3.0 -3.0
       3.0 -9.0
      -3.0 -9.0
      -3.0 -3.0
    endpolygon

    font_size 0.2

    fill_gray 0.0
    moveto -9.0 -10.0
    label This is "polygon_fill.plot".

    moveto 4.0 4.0
    label_slant 45.0 This is slanted text.

  endpage

endfile
