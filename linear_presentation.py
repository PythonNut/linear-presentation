import numpy as np
from gauss_codes import gknot
import os
import copy

def sign(i):
    if i > 0:
        return 1
    elif i < 0:
        return -1
    return 0

def get_cnum(c):
    """
    Crossings in our gauss code are formatted like

    (+/-) (crossing number).(0/0.5)

    where + indicates we're on an overstrand, - indicates we're on an
    understrand, 0 indicates the crossing is positive, and .5 indicates the
    crossing is negative.

    We just want to extract the crossing number here
    """

    if c < 0:
        if c % 1 != 0:
            c += 0.5
    else:
        if c % 1 != 0:
            c -= 0.5

    return int(abs(c))

def get_y_in(cross_seen, c2):
    if not cross_seen:
        return 0

    if c2 < 0:
        if c2 % 1 != 0:
            # vertical understrand, negative crossing
            #
            #      y
            #      |
            #
            # ----------->
            #
            #      |
            #      v
            #
            c2_y = 1
        else:
            # vertical understrand, positive crossing
            #
            #      ^
            #      |
            #
            # ----------->
            #
            #      |
            #      y
            c2_y = -1
    else:
        if c2 % 1 != 0:
            # vertical overstrand, negative crossing
            #
            #      ^
            #      |
            #      |
            # ---- | ---->
            #      |
            #      |
            #      y
            #
            c2_y = -1
        else:
            # vertical overstrand, positive crossing
            #
            #      y
            #      |
            #      |
            # ---- | ---->
            #      |
            #      |
            #      v
            #
            c2_y = 1

    return c2_y

def add_c_shape(path, x1, x2, y_sign):
    y = y_sign * abs(x2 - x1)
    path.append((x1, y))
    path.append((x2, y))
    path.append((x2, 0))

def get_c_shapes(path):
    intervals = []
    for p1, p2 in zip(path, path[1:]):
       (p1x, p1y), (p2x, p2y) = p1, p2
       if p1y == p2y:
           intervals.append((min(p1x, p2x), max(p1x, p2x), sign(p1y)))

    return intervals

def get_envelope(c_shapes, x, y):
    y = sign(y)
    assert y != 0
    min_size, min_c = float('inf'), None
    for (cx1, cx2, cy) in c_shapes:
        if cy == y and cx1 < x and x < cx2:
            size = cx2 - cx1
            assert size != min_size
            if size < min_size:
                min_size = size
                min_c = (cx1, cx2)

    return min_c

def get_path(path, cross_x, gaps, c1i, c2i, x1, x2, y1, y2):
    c_shapes = get_c_shapes(path)
    delta = abs(c2i - c1i)
    if y1 == y2:
        if y1 == 0:
            # easy
            path.append((x2, y2))

        else:
            e1 = get_envelope(c_shapes, x1, y1)
            e2 = get_envelope(c_shapes, x2, y2)
            if e1 == e2:
                # base case
                add_c_shape(path, x1, x2, y1)
            else:
                print('oh no1')
            # y = y1 * (delta//2 + 1)
            # path.append((x1, y))
            # path.append((x2, y))
            # path.append((x2, 0))

    else:
        if y1 == 0:
            path.append((x1 + 1, 0))
            e1 = get_envelope(c_shapes, x1, y2)
            e2 = get_envelope(c_shapes, x2, y2)
            if e1 == e2:
                # base case
                add_c_shape(path, x1+1, x2, y2)
            else:
                print('oh no2')

            # y = y2 * (delta//2 + 1)
            # add_c_shape(path, x1 + 1, x2, y2)
            # path.append((x1 + 1, y))
            # path.append((x2, y))
            # path.append((x2, 0))

        elif y2 == 0:
            # path.append((x1, y1))
            # y = y1 * (delta//2 + 1)
            # path.append((x1, y))

            # # We want to go into the second crossing horizontally
            e1 = get_envelope(c_shapes, x1, y1)
            e2 = get_envelope(c_shapes, x2 - 1, y1)
            if e1 == e2:
                # base case
                add_c_shape(path, x1, x2 - 1, y1)
                path.append((x2, y2))

            else:
                print('oh no3')

            # path.append((x, y))
            # path.append((x,0))
            # path.append((x+1,0))

            # gaps.append(c2i)
            # path.append((x2, y))
            # path.append((x2, 0))

        else:
            print("oh no4")
            # path.append((x1, y1))
            # y = y1 * (delta//2 + 1)
            # path.append((x1, y))

            # print(y1, y2)
            # if x1 < x2:
            #     possible_gaps = list(filter(lambda x: x1 <= x and x <= x2, gaps))
            #     cross_gap = max(possible_gaps)

            #     x_gap = cross_x[cross_gap] - delta
            # else:
            #     possible_gaps = list(filter(lambda x: x2 <= x and x <= x1, gaps))
            #     cross_gap = min(possible_gaps)
            #     x_gap = cross_x[cross_gap] + delta

            # path.append((x_gap, y))
            # path.append((x_gap, -1 * y))
            # path.append((x2, y))
            # path.append((x2, 0))



def build_stupid_graph(gcode):
    # The total number of crossings (n) is going to be 1/2 the total length of
    # the gauss code. We use integer divide to implicitly cast to an int so that
    # we can put it into `range(n)` later.
    n = len(gcode)//2

    # Initialize the drawing path, which is currently a list of tuples
    # representing coordinates along the path that we are going to draw.
    #
    # The current form for these coordinates is (<x>, <y>)
    path = [(0, 0)]

    # A list giving us the x coordinates at which each crossing will be drawn.
    # Initially, we will make the drawing grid _way_ larger than it needs to be,
    # to account for some edge cases where we need a lot of horizontal spacing.
    #
    # This occurs when we have performed some "backtracking" in the knot, i.e.
    # when we encounter a crossing we've seen already before we've encountered
    # all crossings in the diagram. E.g., in the knot (6,4).
    cross_x = [i*n for i in range(n)]

    # Keep track of crossings we've seen already.
    cross_seen = [False for i in range(n)]

    # Keep track of the places where we left gaps for the backtracking process
    gaps = []

    gcode = copy.deepcopy(gcode)
    gcode.append(gcode[0])

    for i in range(len(gcode) - 1):
        # Get information from the gauss code
        c1, c2 = gcode[i], gcode[i+1]
        c1i = get_cnum(c1) - 1
        c2i = get_cnum(c2) - 1
        x1, x2 = cross_x[c1i], cross_x[c2i]
        y1 = get_y_in(cross_seen[c1i], c1)
        y2 = -1 * get_y_in(cross_seen[c2i], c2)
        # print(cross_seen)
        cross_seen[c1i] = True

        print(f"c: {c1} → {c2}, x: {x1} → {x2}, y: {y1} → {y2}")

        get_path(path, cross_x, gaps, c1i, c2i, x1, x2, y1, y2)

    # fix the termination FIXME!
    path[-2] = (-1, path[-2][1])
    path[-1] = (-1, path[-1][1])
    path.append((0, 0))

    signs = []
    for cnum in range(1,n):
        for c in gcode:
            if get_cnum(c) == cnum:
                signs.append((c/abs(c)))

    return path, cross_x, signs


def draw_presentation(draw_paths, x_vals, signs, fname="test_gauss"):
    """
    draw

    draw_paths a list of lists where each sublist is a list of tuples specifying
    a path to draw in TikZ

    x_vals is
    """
    break_width = 0.15

    # Get the LaTeX preamble
    preamble = "\\documentclass[border=10pt]{standalone}\n"
    preamble += "\\usepackage{tikz}\n"
    preamble += "\\begin{document}\n"
    preamble += "\\begin{tikzpicture}\n"

    drawing = ""
    for path in draw_paths:
        drawing += "\\draw[-latex]"
        for vert in path:
            drawing += str(vert) + " -- "
        # Remove the ` -- ` from the last coordinate, terminate the path, and
        # break the line
        drawing = drawing[:-4] + ";\n"
    for cnum, x in enumerate(x_vals):
        x = x_vals[cnum]
        sign = signs[cnum]
        drawing += f"\\fill[white] ({x-break_width}, -{break_width}) rectangle ({x + break_width}, {break_width});\n"
        if sign < 0:
            drawing += f"\\draw ({x}, {break_width}) -- ({x}, -{break_width});\n"
        else:
            drawing += f"\\draw ({x-break_width}, 0) -- ({x+break_width}, 0);\n"

        drawing += f"\\node[circle, fill=black, draw=black, inner sep=1pt] () at ({x}, 0) " + "{};\n"
        drawing += f"\\node[above left] () at ({x}, 0) " + "{\\small $" + str(int(cnum+1)) + "$};\n"

    out_str = preamble + drawing + "\\end{tikzpicture}\n\\end{document}"

    # Need to chdir else pdflatex will pollute the parent directory with aux
    # files and stuff
    os.chdir("tests")
    with open(f"{fname}.tex", "w") as f:
        f.write(out_str)
    os.system(f"pdflatex -halt-on-error {fname}.tex | grep '^!.*' -A200 --color=always")
    print(fname)
    os.system(f"zathura {fname}.pdf")
    os.chdir("..")
    return out_str

if __name__ == '__main__':
    # knot = gknot[(6,5)]
    # knot = gknot[(8, 4)]
    # knot = gknot[(3, 1)]
    # knot = gknot[(6,4)]
    # knot = gknot[(11,2)]
    # knot = gknot[(5,2)]
    pathological_test = [
        -1.5, 2, -3, 4.5, -5.5, 3, -6, 1.5, 7, 5.5, -4.5, 6, -2, -7
    ]
    path, cross_x, signs = build_stupid_graph(pathological_test)
    draw_presentation([path], cross_x, signs)
