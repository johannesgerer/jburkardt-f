# triangle_orientation.plot  21 January 2011
#
file triangle_orientation.ps

  space -0.50 -0.50 10.50 7.50

  page

    line_rgb 0.5 0.5 0.5

    grid 0.0 0.0 10.0 7.0 11 8

    line_thick 3

    fill_rgb 0.30 0.95 0.95
    triangle_fill 2.0 2.0 8.0 2.0 5.0 5.0

    fill_rgb 0.90 0.90 0.40

    polygon_fill
      0.0 2.0
      2.0 2.0
      5.0 5.0
      3.0 7.0
      0.0 7.0
    end

    polygon_fill
     10.0 2.0
      8.0 2.0
      5.0 5.0
      7.0 7.0
     10.0 7.0
    end

    polygon_fill
      0.0 0.0
      2.0 2.0
      8.0 2.0
     10.0 0.0
    end

    fill_rgb 0.95 0.50 0.50
    triangle_fill  0.0 0.0  2.0 2.0 0.0 2.0
    triangle_fill  0.0 0.0  2.0 2.0 0.0 2.0
    triangle_fill  5.0 5.0  7.0 7.0 3.0 7.0
    triangle_fill 10.0 0.0 10.0 2.0 8.0 2.0
#
#  Draw triangle boundaries.
#
    fill_rgb 0.8 0.8 0.8

    moveto   0.0  0.0
    drawto   7.0  7.0
    moveto  10.0  0.0
    drawto   3.0  7.0
    moveto   0.0  2.0
    drawto  10.0  2.0
#
#  Label the plot.
#
    font_size 0.35
    line_rgb 0.0 0.0 0.0

    moveto 1.0 3.0
    label ++-
    moveto 4.6 3.0
    label +++
    moveto 8.6 3.0
    label +-+

    moveto 4.6 1.0
    label -++
    moveto 4.6 6.0
    label +--

    moveto 0.3 1.2
    label -+-
    moveto 9.4 1.2
    label --+

  endpage

endfile
