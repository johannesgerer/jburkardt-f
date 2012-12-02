# manhole.plot  29 March 2001
#
file manhole.ps

  space -1.0 -1.0 1.0 1.0

  page

    line_rgb 1.0 0.0 0.0

    moveto  0.5  0.0
    drawto  0.0  0.866 
    drawto -0.5  0.0
    drawto  0.5  0.0

    line_rgb 0.0 1.0 0.0

    arc  0.5 0.0   1.0 120 180
    arc  0.0 0.866 1.0 240 300
    arc -0.5 0.0   1.0   0  60

  endpage

endfile
