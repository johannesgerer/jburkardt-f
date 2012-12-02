# graph_paths.plot  30 December 2010
#
file graph_disconnected.ps

  space -0.5 -0.5 4.5 3.5

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 4.0 3.0 5 4
#
#  Red things (lines)
#
    fill_rgb 0.9 0.0 0.1
    line_width 3

    line  0.0  3.0  1.0  3.0
    line  3.0  2.0  4.0  2.0
    line  2.0  1.0  3.0  1.0
    line  0.0  0.0  1.0  0.0
    line  2.0  0.0  3.0  0.0

    line  0.0  2.0  0.0  3.0
    line  0.0  1.0  0.0  2.0
    line  0.0  0.0  0.0  1.0
    line  1.0  2.0  1.0  3.0
    line  2.0  2.0  2.0  3.0
    line  3.0  0.0  3.0  1.0
    line  3.0  2.0  3.0  3.0
    line  4.0  0.0  4.0  1.0
    line  4.0  1.0  4.0  2.0

    line  0.0  2.0  1.0  3.0
    line  0.0  1.0  1.0  0.0
    line  1.0  0.0  2.0  1.0
    line  1.0  1.0  2.0  2.0
    line  1.0  1.0  2.0  3.0
    line  3.0  1.0  4.0  0.0
    line  3.0  1.0  4.0  2.0
    line  3.0  2.0  4.0  3.0
#
#  Blue things (nodes)
#
    fill_rgb 0.0 0.1 0.9

    circle_fill  0.0  0.0  0.15
    circle_fill  1.0  0.0  0.15
    circle_fill  2.0  0.0  0.15
    circle_fill  3.0  0.0  0.15
    circle_fill  4.0  0.0  0.15

    circle_fill  0.0  1.0  0.15
    circle_fill  1.0  1.0  0.15
    circle_fill  2.0  1.0  0.15
    circle_fill  3.0  1.0  0.15
    circle_fill  4.0  1.0  0.15

    circle_fill  0.0  2.0  0.15
    circle_fill  1.0  2.0  0.15
    circle_fill  2.0  2.0  0.15
    circle_fill  3.0  2.0  0.15
    circle_fill  4.0  2.0  0.15

    circle_fill  0.0  3.0  0.15
    circle_fill  1.0  3.0  0.15
    circle_fill  2.0  3.0  0.15
    circle_fill  3.0  3.0  0.15
    circle_fill  4.0  3.0  0.15
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0
#
#  Node labels
#
    font_size 0.40

    moveto  0.1  3.25
    label A
    moveto  1.1  3.25
    label B
    moveto  2.1  3.25
    label C
    moveto  3.1  3.25
    label D
    moveto  4.1  3.25
    label E
    moveto  0.1  2.25
    label F
    moveto  1.1  2.25
    label G
    moveto  2.1  2.25
    label H
    moveto  3.1  2.25
    label I
    moveto  4.1  2.25
    label J
    moveto  0.1  1.25
    label K
    moveto  1.1  1.25
    label L
    moveto  2.1  1.25
    label M
    moveto  3.1  1.25
    label N
    moveto  4.1  1.25
    label O
    moveto  0.1  0.25
    label P
    moveto  1.1  0.25
    label Q
    moveto  2.1  0.25
    label R
    moveto  3.1  0.25
    label S
    moveto  4.1  0.25
    label T

  endpage

endfile
