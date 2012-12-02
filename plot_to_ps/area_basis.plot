# area_basis.plot  01 February 2007
#
file area_basis.ps

  space -0.05 -0.05 1.05 1.05

  page

    line_rgb 1.0 0.0 0.0
#
#  Draw the triangle edges.
#
    moveto  0.0  0.0
    drawto  1.0  0.2
    drawto  0.3  0.9
    drawto  0.0  0.0
#
#  Draw lines from the vertices to P.
#
    moveto  0.0  0.0
    drawto  0.5  0.5
    moveto  1.0  0.2
    drawto  0.5  0.5
    moveto  0.3  0.9
    drawto  0.5  0.5
#
#  Mark the vertices and P.
#
    line_thick 5

    fill_rgb 0.5 0.7 0.5
    circle_fill 0.0 0.0 0.025
    circle_fill 1.0 0.2 0.025
    circle_fill 0.3 0.9 0.025
    circle_fill 0.5 0.5 0.025
#
#  Label the vertices and P.
#
    font_size 0.35
    fill_rgb 0.0 0.0 0.0

    moveto  0.0  0.05
    label A
    moveto  1.0  0.25
    label B
    moveto  0.3  0.95
    label C
    moveto  0.5  0.55
    label P

    moveto  0.60  0.60
    label {\xi_A}
    moveto  0.30  0.50
    label {\xi_B}
    moveto  0.5  0.28
    label \{xi_C}

  endpage

endfile
