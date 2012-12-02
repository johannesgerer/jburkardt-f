# graph_mst2.plot  19 January 2011
#
file graph_mst2.ps

  space -0.5 -0.5 9.5 6.5

  page

    line_rgb 0.8 0.8 0.8
    grid 0.0 0.0 9.0 6.0 10 7
#
#  Red things
#
    fill_rgb 1.0 0.0 0.0
    line_width 3

    line 0.0  6.0  1.0  4.0
    line 0.0  6.0  0.0  0.0
    line 0.0  6.0  5.0  4.0
    line 1.0  4.0  3.0  4.0
    line 1.0  4.0  1.0  2.0
    line 1.0  4.0  3.0  2.0
    line 3.0  4.0  3.0  2.0
    line 1.0  2.0  3.0  2.0
    line 1.0  2.0  0.0  0.0
    line 3.0  2.0  5.0  4.0
    line 3.0  2.0  7.0  0.0
    line 0.0  0.0  7.0  0.0
    line 5.0  4.0  7.0  5.0
    line 5.0  4.0  7.0  2.0
    line 5.0  4.0  7.0  0.0
    line 7.0  5.0  9.0  5.0
    line 9.0  5.0  9.0  2.0
    line 7.0  2.0  9.0  2.0
    line 7.0  2.0  7.0  0.0
    line 7.0  2.0  9.0  0.0
    line 7.0  0.0  9.0  0.0

#
#  Blue things
#
    fill_rgb 0.0 0.0 1.0

    circle_fill 0.0  6.0  0.25
    circle_fill 1.0  4.0  0.25
    circle_fill 3.0  4.0  0.25
    circle_fill 1.0  2.0  0.25
    circle_fill 3.0  2.0  0.25
    circle_fill 0.0  0.0  0.25
    circle_fill 5.0  4.0  0.25
    circle_fill 7.0  5.0  0.25
    circle_fill 9.0  5.0  0.25
    circle_fill 7.0  2.0  0.25
    circle_fill 9.0  2.0  0.25
    circle_fill 7.0  0.0  0.25
    circle_fill 9.0  0.0  0.25
#
#  Black things
#
    fill_rgb 0.0 0.0 0.0
#
#  Node labels
#
    font_size 0.40

    moveto  0.3  6.0
    label A
    moveto  1.2  4.2
    label B
    moveto  3.2  4.2
    label C
    moveto  1.2  2.2
    label D
    moveto  3.2  2.2
    label E
    moveto  0.2  0.2
    label F
    moveto  5.2  4.2
    label G
    moveto  7.2  5.2
    label H
    moveto  9.2  5.2
    label I
    moveto  7.2  2.2
    label J
    moveto  9.1  2.3
    label K
    moveto  7.2  0.2
    label L
    moveto  9.1  0.3
    label M
#
#  Line labels
#
    font_size 0.30

    moveto 0.5 5.0
    label  23
    moveto 0.0  3.0
    label 60
    moveto 2.5  5.1
    label 54
    moveto 2.0  4.1
    label 19
    moveto 1.1  3.0
    label 20
    moveto 2.0  3.0
    label 26
    moveto 3.1  3.0
    label 20
    moveto 2.0  2.1
    label 21
    moveto 0.6  1.0
    label 22
    moveto 4.3  3.0
    label 29
    moveto 5.0  1.0
    label 46
    moveto 3.5  0.1
    label 70
    moveto 5.7  4.7
    label 22
    moveto 6.1  3.0
    label 28
    moveto 6.0  2.0
    label 44
    moveto 8.1  5.1
    label 20
    moveto 8.4  3.5
    label 30
    moveto 8.1  2.1
    label 20
    moveto 7.0  1.0
    label 10
    moveto 8.1  1.0
    label 27
    moveto 8.0  0.1
    label 18

  endpage

endfile
