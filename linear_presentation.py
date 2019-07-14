import numpy as np



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
            return abs(c + 0.5)
        else:
            return abs(c)
    else:
        if c % 1 != 0:
            return abs(c - 0.5)
        else:
            return c

def get_c2_y(c2):
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
