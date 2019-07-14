import numpy as np
from gauss_codes import gknot


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

def get_c2_y(cross_seen, c2):
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
            #
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

def get_path(path, cross_seen, x1, x2, y1, y2):
    pass


def build_stupid_graph(gcode):
    n = len(gcode)//2
    path = [(-1, 0, None)]
    cross_x = [i*n for i in range(n)]
    cross_seen = [False for i in range(n)]

    # if gcode[0] < 0:
    #     # going under
    #     path.append((0, 0, -1))

    # else:
    #     # going over
    #     path.append((0, 0, +1))

    # path.append((1, 0, None))

    for i in range(len(gcode) - 1):
        c1, c2 = gcode[i], gcode[i+1]
        c1i = get_cnum(c1) - 1
        c2i = get_cnum(c2) - 1
        x1, x2 = cross_x[c1i], cross_x[c2i]
        y1 = get_c2_y(cross_seen[c1i], c1)
        y2 = get_c2_y(cross_seen[c2i], c2)
        cross_seen[c1i-1] = True

        print(f"c: {c1} → {c2}, x: {x1} → {x2}, y: {y1} → {y2}")

        get_path(path, cross_seen, x1, x2, y1, y2)



def draw_presentation(draw_paths, x_vals, fname="test_gauss"):
    """
    draw

    draw_paths a list of lists where each sublist is a list of tuples specifying
    a path to draw in TikZ

    x_vals is
    """
    # Get the LaTeX preamble
    preamble = "\\documentclass[border=10pt]{standalone}\n"
    preamble += "\\usepackage{tikz}\n"
    preamble += "\\begin{document}\n"
    preamble += "\\begin{tikzpicture}\n"



    drawing = ""
    for path in draw_paths:
        drawing += "\\draw[-latex]"
        # drawing += "\\draw"
        for vert in path:
            drawing += str(vert) + " -- "
        # Remove the ` -- ` from the last coordinate, terminate the path, and
        # break the line
        drawing = drawing[:-4] + ";\n"
    for cnum, x in x_vals:
        x = x_vals[cnum]
        drawing += f"\\node[circle, draw=black, inner sep=1pt] () at ({x}, 0) " + "{};\n"
        drawing += f"\\node[above left] () at ({x}, 0) " + "{\\small $" + str(int(cnum)) + "$};\n"
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
    knot = gknot[(8, 4)]
    build_stupid_graph(knot)
