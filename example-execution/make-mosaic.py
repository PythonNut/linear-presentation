from math import ceil, sqrt
from os import listdir, chdir
import subprocess


def get_prefix(fname):
    """
    Given something like
        `3-1.pdf`
    Extract and return the tuple
        `(3,1)`
    """
    # Crossing number and index of knot in rolfsen table.
    # OK: the `.split(".")[0]` gets the portion before the
    # `.pdf` in the filename.
    #
    # The .split("-") separates the two pieces of data (e.g.,
    # "8-18" gets split to ["8", "18"])
    #
    # The list comprehension is just to convert both to ints.
    return tuple([int(glyph) for glyph in fname.split(".")[0].split("-")])


# def sort_knots(knotlist):
#     """
#     Sort the knots by their indices in the rolfsen table.

#     Use mergesort.

#     I think there's a way to do this with the built-in `sorted()`
#     function but I'm too lazy so I'll just write a merge sort.
#     """
#     # Empty list and 1-element lists are sorted
#     if len(knotlist) <= 1:
#         return knotlist
#     else:
#         # knotlist length
#         kl_len = len(knotlist)

#         # Halfway index
#         hi = kl_len // 2

#         left = sort_knots[knotlist[:hi]]
#         right = sort_knots[knotlist[hi:]]

#         while (left and right):
#             lpref = get_prefix(left[0])
#             rpref = get_prefix(right[0])

#             # Python supports tuple comparison out-of-the-box
#             if lpref < rpref:


#     return


# Filter for pdfs with a `-` in the name
knots = [fname for fname in listdir("finals") if (".pdf" in fname) and ("-" in fname)]


knots = sorted(knots, key=get_prefix)

l = ceil(sqrt(len(knots)))
spacing = 0.9 / l


out_str = r"""\documentclass{article}
\usepackage{fkmath}
\usepackage{float}
\usepackage{subcaption}
\usepackage{graphicx}

\allowdisplaybreaks
\pagestyle{empty}
"""

# Comment this out for standard colors:
out_str += r"""\usepackage{xcolor} % For customizing page color and such
\definecolor{bcol}{HTML}{1D252C}
\definecolor{tcol}{HTML}{D6D7D9}
\pagecolor{bcol}
\color{tcol}
"""
out_str += (
    "\\usepackage[margin=.5in, landscape, papersize={"
    + str(1 + 2 * l)
    + "in, "
    + str(1 + 2 * l)
    + "in}]{geometry}"
)

out_str += r"\begin{document}"

# Lay out the knots in a grid table
for row in range(l):
    out_str += "\\begin{figure}[H]\n\\centering\n"
    for col in range(l):
        index = row * l + col
        if index == len(knots):
            out_str += r"\hfill"
        elif index > len(knots):
            continue
        else:
            fname = knots[index]
            label = get_prefix(fname)

            out_str += r"\begin{subfigure}[t]{" + str(spacing) + r"\textwidth}" + "\n"
            out_str += r"\centering" + "\n"
            out_str += r"\includegraphics[width = .95in]{" + fname + "}\n"
            out_str += r"\caption*{$" + f"{label}" + "$}" + "\n"
            out_str += r"\end{subfigure}" + "\n"

    out_str += "\\end{figure}\n\\vfill\n"
out_str += r"\end{document}"

with open("finals/mosaic.tex", "w") as f:
    f.write(out_str)

chdir("finals")
subprocess.run(["pdflatex", "mosaic.tex"])
chdir("..")
