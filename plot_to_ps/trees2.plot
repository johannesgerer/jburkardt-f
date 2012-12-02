# trees2.plot  13 October 2000
#
file trees2.ps

  space -10.0 -10.0 10.0 10.0

  page

    line_thick 2

    fill_rgb 0.0 0.8 0.2
    line_rgb 0.2 0.2 0.8

#   circle_fill 0.0 0.0 0.2

    star      0.0 0.0  10.00
    star_disk 0.0 0.0  10.00
#
#  The smaller star is shrunk by a factor of PHI**2.
#
    star      0.0 0.0  -3.82
    star_disk 0.0 0.0  -3.82

  endpage

endfile
