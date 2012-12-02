# histogram.plot  28 August 1999
#
file histogram.ps

  space -2.0 -10.0 12.0 10.0

  page

    line_gray 0.7

    grid -1.0 -10.0 12.0 10.0 14 11

    fill_rgb 0.8 0.3 0.1
    line_rgb 0.1 0.8 0.2

    histo_fat = 0.40

    histogram
       0.0 -3.0
       1.0  3.0
       2.0  4.0
       3.0  8.0
       4.0  9.0
       5.0  7.0
       6.0  3.0
       7.0 -3.0
       8.0 -2.0
       9.0 -9.0
      10.0 -5.0
      11.0 -3.0
    endhistogram

    font_size 0.2

    fill_gray 0.0
    moveto 0.0 -10.0
    label This is "histogram.plot".

  endpage

endfile
