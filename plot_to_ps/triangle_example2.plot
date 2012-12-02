# triangle_example2.plot  09 June 2009
#
file triangle.ps

  space -1.0 -1.0 11.0 11.0

  page

    line_rgb 0.3 0.3 0.3

    grid -1.0 -1.0 11.0 11.0 13 13

    fill_rgb 0.35 0.95 0.9

    triangle_fill 0.0 0.0 10.0 0.0 5.0 8.66

    line_thick 3

    moveto   0.0  0.0
    drawto  10.0  0.0
    drawto   5.0  8.66
    drawto   0.0  0.0

    fill_rgb 0.5 0.7 0.5

    circle_fill  0.0  0.0  0.25
    circle_fill 10.0  0.0  0.25
    circle_fill  5.0  8.66 0.25

    font_size 0.35
    line_rgb 0.0 0.0 0.0

    moveto  0.5  -0.5
    label Tex2.Va = {0, 0}
    moveto  6.5 1.0
    label Tex2.Vb = {10, 0}
    moveto  5.0  9.0
    label Tex2.Vc = {5, 8.66}

  endpage

endfile
