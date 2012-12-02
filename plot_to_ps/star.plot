# star.plot  19 October 2007
#
#  Draw a star with its vertices marked by disks.
#
file star.ps

  space 0.0 0.0 100.0 100.0

  page

    line_width 3
    star 50.0 50.0 45.0
    star_disk 50.0 50.0 45.0

  endpage

endfile
