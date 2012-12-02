# points.plot  12 February 2001
#
debug
file points.eps

  space -3.0 -1.0 12.0 10.0

  page

    line_rgb 0.2 0.2 0.2
    grid -3.0 -1.0 12.0 10.0 16 12

    line_rgb 0.5 0.5 0.5

    point_list
      0.0  0.0
      2.0  2.0
     -1.0  3.0
     -2.0  2.0
      0.0  0.0
      8.0  2.0
      9.0  5.0
      7.0  4.0
      8.0  2.0
      5.0  6.0
      6.0  7.0
      8.0  8.0
     11.0  7.0
     10.0  4.0
      8.0  2.0
      6.0  4.0
      5.0  6.0
    end

    moveto  0.0  8.0    
    font_size 0.50
    label Point List plot

  endpage

endfile
