# double_circle.plot  14 November 2005
#
file double_circle.eps

  space -2.0 -2.0 2.0 2.0

  page

    moveto 0.0 0.707
    font_size 0.125
    label Center #1
    point  0.0 0.707
    circle 0.0 0.707 1.0

    moveto 0.0 -0.707
    font_size 0.125
    label Center #2
    point  0.0 -0.707
    circle 0.0 -0.707 1.0

    moveto -0.707 0.0
    font_size 0.125
    label P1
    point -0.707 0.0

    moveto 0.707 0.0
    font_size 0.125
    label P2
    point  0.707 0.0

    moveto -1.0 2.0
    font_size 0.50
    label P1 and P2 lie on two circles.

  endpage

endfile
