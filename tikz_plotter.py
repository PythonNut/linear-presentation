import os
import itertools as it
import subprocess
from math import sqrt, ceil


def crossing_tex(cnum, orient, x, xcurr, add_labels=False):
    """
    `cnum`: an int representing which crossing we're drawing. cnum < 0
            corresponds to the horizontal strand being an understrand,
            cnum > 0 is the same but an overstrand.
    `orient`: either +1 or -1, corresponding to positive, negative
              crossings respectively.
    `x`: the x value where we're drawing our crossing.
    `xcurr`: the x value for the vertical dashed line
    """
    out_str = ""
    gap = 0.2
    if cnum < 0:
        out_str += f"    \\draw ({x - 1}, 0) -- ({x - gap}, 0);\n"
        out_str += f"    \\draw[-latex] ({x + gap}, 0) -- ({x + .5}, 0);\n"
        out_str += f"    \\draw ({x + .5}, 0) -- ({x + 1}, 0);\n"
        if orient == 1:
            out_str += f"    \\draw[-latex] ({x}, .5) -- ({x}, -.5);\n"
        elif orient == -1:
            out_str += f"    \\draw[-latex] ({x}, -.5) -- ({x}, .5);\n"
        else:
            assert False
    else:
        out_str += f"    \\draw ({x -1}, 0) -- ({x + 1}, 0);\n"
        out_str += f"    \\draw[-latex] ({x + gap}, 0) -- ({x + .5}, 0);\n"
        if orient == 1:
            out_str += f"    \\draw ({x}, -.5) -- ({x}, -{gap});\n"
            out_str += f"    \\draw[-latex] ({x}, {gap}) -- ({x}, .5);\n"
        elif orient == -1:
            out_str += f"    \\draw ({x}, .5) -- ({x}, {gap});\n"
            out_str += f"    \\draw[-latex] ({x}, -{gap}) -- ({x}, -.5);\n"
        else:
            assert False

    if add_labels:

        lab = f"$\\small {abs(cnum)}"
        if cnum < 0:
            lab += "_u"
        else:
            lab += "_o"

        if orient == 1:
            lab += "^+"
        else:
            lab += "^-"

        lab += "$"

        if orient == 1:
            out_str += (
                f"    \\node[circle, fill=black, draw=black, inner sep=1pt] ({abs(cnum)}circ) at ({x}, 0) "
                + "{};\n"
            )

            out_str += (
                f"    \\node[above left] ({abs(cnum)}lab) at ({x}, 0) "
                + "{"
                + lab
                + "};\n"
            )

        else:
            out_str += (
                f"    \\node[circle, fill=white, draw=black, inner sep=1pt] ({abs(cnum)}circ) at ({x}, 0) "
                + "{};\n"
            )

            out_str += (
                f"    \\node[above left] ({abs(cnum)}lab) at ({x}, 0) "
                + "{"
                + lab
                + "};\n"
            )

    return out_str


def get_thalf(x0, x1, xcurr):
    """
    Given x0 <= xcurr <= x1, compute which theta in [0,180] yields
    xcurr when we take ((x1 - x0)/2)*cos(theta)
    """
    xavg = (x0 + x1) / 2
    arg = 2 * (xcurr - xavg) / (x1 - x0)
    return gp(180 * arccos(arg) / pi)  # convert to degrees and cast to float


def arc_tex(
    x0, x1, xcurr, y, s, x0_through=False, x1_through=False, label_semiarcs=False
):
    """

    """
    out_str = ""
    r = abs(x0 - x1) / 2
    xavg = (x1 + x0) / 2

    if y > 0:
        # if xcurr < x1:
        #     nstyle = ", opacity=.2"
        # else:
        nstyle = ""

        if label_semiarcs:

            if r == 0.5:
                out_str += (
                    f"    \\node[above{nstyle}] ({s}lab) at ({xavg}, 0) "
                    + "{$s_{"
                    + str(s)
                    + "}$};\n"
                )
            else:
                out_str += (
                    f"    \\node[above{nstyle}] ({s}lab) at ({(x1+x0)/2}, {r+.5}) "
                    + "{$s_{"
                    + str(s)
                    + "}$};\n"
                )

    # Handle the arcs that connect things directly along the horizontal
    if r == 0.5:
        # if xcurr < x0:
        #     out_str += f"    \\draw[opacity=.2] ({x0-.5}, 0) -- ({x1+.5}, 0);\n"
        # elif x1 < xcurr:
        out_str += f"    \\draw ({x0 - .5}, 0) -- ({x1 + .5}, 0);\n"
        # else:
        #     out_str += f"    \\draw ({x0 - .5}, 0) -- ({xcurr}, 0);\n"
        #     out_str += f"    \\draw[opacity=.2] ({xcurr}, 0) -- ({x1 + .5}, 0);\n"
        return out_str

    if x0 < x1:
        t0 = 180
        t1 = 0
    else:
        t0 = 0
        t1 = 180

    # if xcurr < x0:
    #     out_str += f"    \\draw[opacity=.2] ({x0}, {y}) arc ({t0}:{t1}:{r});\n"
    #     # print(x0, x1, xcurr, y, s)
    #     # assert False

    # elif xcurr < x1:
    #     # Find the angle to draw until such that our arc meets the
    #     # dotted line
    #     thalf =
    #     out_str += f"    \\draw[opacity=.6] ({x0},  {y}) arc ({t0}:{thalf}:{r});\n"
    #     out_str += f"    \\draw[opacity=.2] ({x1},  {y}) arc ({t1}:{thalf}:{r});\n"

    # else:
    out_str += f"    \\draw ({x0}, {y}) arc ({t0}:{t1}:{r});\n"

    if x0_through:
        if xcurr < x0:
            out_str += f"\\draw ({x0}, {y}) -- ({x0}, 0);\n"
        else:
            out_str += f"\\draw ({x0}, {y}) -- ({x0}, 0);\n"

    if x1_through:
        if xcurr < x1:
            out_str += f"\\draw ({x1}, {y}) -- ({x1}, 0);\n"
        else:
            out_str += f"\\draw ({x1}, {y}) -- ({x1}, 0);\n"

    return out_str


def stacks_tex(upper, lower, nmax):
    """
    Generate the code for drawing the TikZ representation of our
    stacks (upper and lower)
    """
    nmax += 2
    out_str = ""

    # Upper stack
    out_str += f"    \\draw[thick] (-2, .5) -- (-2, {nmax/2+.5}) -- (-1, {nmax/2+.5}) -- (-1, .5) ;\n"
    # Lower stack
    out_str += f"    \\draw[thick] (-2, -.5) -- (-2, -{nmax/2+.5}) -- (-1, -{nmax/2+.5}) -- (-1, -.5) ;\n"

    # Draw the horizontal lines demarking parts of the stack
    out_str += "    \\draw "
    for h in range(1, nmax):
        h /= 2
        h += 0.5
        out_str += f"(-2, {h}) -- (-1, {h}) (-2, -{h}) -- (-1, -{h})"
    out_str += ";\n"

    for (h, s) in enumerate(upper[::-1]):
        h /= 2
        h += 0.75
        out_str += (
            f"    \\node () at (-1.5, {h}) " + "{\\small $s_{" + str(s) + "}$};\n"
        )

    for (h, s) in enumerate(lower[::-1]):
        h /= 2
        h += 0.75
        out_str += (
            f"    \\node () at (-1.5, -{h}) " + "{\\small $s_{" + str(s) + "}$};\n"
        )

    return out_str


def plot_one_step(
    gc,
    orient,
    state,
    i,
    upper_cs_f,
    lower_cs_f,
    crossings_f,
    dirname,
    draw_stacks=True,
    draw_line=True,
    use_fk_colors=False,
):
    """
    Plot a single step in the execution of the algorithm.

    `gc`: Together with `orient`, comprises the sage-format gauss code
    `orient`: see above
    `state`: a single entry in the `states` list defined in `route()`.
    `i`: the current step of the computation. Used to generate the
         filename.
    `<whatever>_f`: the version of <whatever> that gets spat out at
                    the end of `route()`. `_f` is for "final"

    """
    fname = f"{dirname}/{i:03}.tex"

    nmax = len(gc) // 2
    # print(gc, len(gc), nmax)

    out_str = r"""\documentclass[border=1pt]{standalone}
\usepackage{tikz}"""

    if use_fk_colors:
        out_str += r"""\usepackage{xcolor} % For customizing page color and such
\definecolor{bcol}{HTML}{1D252C}
\definecolor{tcol}{HTML}{D6D7D9}
\pagecolor{bcol}
\color{tcol}"""

    out_str += r"""\begin{document}
  \begin{tikzpicture}[every draw/.style={line width=20pt}]
    """
    x, upper, lower, upper_cs, lower_cs, crossings = state
    # print(state)
    if draw_line:
        out_str += f"    \\draw[dotted] ({x}, -{2*nmax}) -- ({x}, {2*nmax});\n"

    num_cs_so_far = len(crossings)

    seen_cs = set()
    i = 0
    # Only draw the crossings we've seen so far in black
    while len(seen_cs) < num_cs_so_far:
        c = gc[i]
        if abs(c) not in seen_cs:
            out_str += crossing_tex(c, orient[len(seen_cs)], crossings[len(seen_cs)], x)
            seen_cs.add(abs(c))
        i += 1

    # Draw the stacks
    if draw_stacks:
        out_str += stacks_tex(upper, lower, nmax)

    # Ok now we draw the remaining crossings in gray
    out_str += "    \\begin{scope}[every path/.style={opacity=.2}, every node/.style={opacity=.2}]\n"
    while i < len(gc):
        c = gc[i]
        if abs(c) not in seen_cs:
            out_str += crossing_tex(
                c, orient[len(seen_cs)], crossings_f[len(seen_cs)], x
            )
            seen_cs.add(abs(c))
        i += 1

    out_str += "    \\end{scope}"

    # Now we draw the semicircles
    for s in upper_cs_f.keys():
        # Special case for if we get the first arc
        if upper_cs_f[s] == []:
            pairs = lower_cs_f[s]
        else:
            pairs = upper_cs_f[s]

        for (x0, x1) in pairs:
            if x0 > x1:
                (x0, x1) = (x1, x0)

            # Check wheter we need to draw the arc down through the horizontal
            if x0 in crossings_f:
                x0_t = False
            else:
                x0_t = True

            if x1 in crossings_f:
                x1_t = False
            else:
                x1_t = True

            # If not in the special case, do the regular call
            if upper_cs_f[s] != []:
                out_str += arc_tex(x0, x1, x, 0.5, s, x0_through=x0_t, x1_through=x1_t)

    out_str += "\\begin{scope}[yscale=-1]"
    for s in lower_cs_f.keys():
        # Special case for if we get the first arc
        if lower_cs_f[s] == []:
            pairs = upper_cs_f[s]
        else:
            pairs = lower_cs_f[s]

        for (x0, x1) in pairs:
            if x0 > x1:
                (x0, x1) = (x1, x0)

            # Check wheter we need to draw the arc down through the horizontal
            if x0 in crossings_f:
                x0_t = False
            else:
                x0_t = True

            if x1 in crossings_f:
                x1_t = False
            else:
                x1_t = True

            # If not in the special case, do the regular call
            if lower_cs_f[s] != []:
                out_str += arc_tex(x0, x1, x, 0.5, s, x0_through=x0_t, x1_through=x1_t)

    out_str += "\\end{scope}"
    out_str += "  \\end{tikzpicture}\n\\end{document}"
    with open(fname, "w") as f:
        f.write(out_str)


def circle_int_point(xl1, xr1, xl2, xr2):
    # print(xl1, xr1)
    # print(xl2, xr2)
    c1 = (xr1 + xl1) / 2
    c2 = (xr2 + xl2) / 2

    r1 = abs(xr1 - c1)
    r2 = abs(xr2 - c2)

    d = abs(c2 - c1)

    try:
        assert d != 0
    except AssertionError as e:
        print(f"xl1: {xl1}, xr1: {xr1}\nxl2: {xl2}, xr2: {xr2}")
        raise (e)

    xint = (d ** 2 - r2 ** 2 + r1 ** 2) / (2 * d) + c1  # - 0.5
    yint = (
        sqrt(abs(4 * (d ** 2) * (r1 ** 2) - (d ** 2 - r2 ** 2 + r1 ** 2) ** 2))
        / (2 * d)
        + 0.5
    )

    # print(f"yint: {yint}")

    return (xint, yint)


def compute_intersections(upper_cs, lower_cs):
    # Get all of the upper / lower semiarcs
    all_upper = []

    int_pts = []

    for key in upper_cs:
        all_upper.extend(upper_cs[key])

    for (xl1, xr1), (xl2, xr2) in it.combinations(all_upper, 2):
        if (xl1 <= xl2) and (xr2 <= xr1):
            continue
        elif (xl2 < xl1) and (xr1 < xr2):
            continue
        elif (xr1 < xl2) or (xr2 < xl1):
            continue
        else:
            int_pts += [circle_int_point(xl1, xr1, xl2, xr2)]

    all_lower = []
    for key in lower_cs:
        all_lower.extend(lower_cs[key])

    for (xl1, xr1), (xl2, xr2) in it.combinations(all_lower, 2):
        if (xl1 <= xl2) and (xr2 <= xr1):
            continue
        elif (xl2 < xl1) and (xr1 < xr2):
            continue
        elif (xr1 < xl2) or (xr2 < xl1):
            continue
        else:
            xint, yint = circle_int_point(xl1, xr1, xl2, xr2)
            int_pts += [(xint, -yint)]

    return int_pts


def plot_virtual(
    gc,
    orient,
    upper_cs,
    lower_cs,
    crossings,
    dirname,
    dir_prefix="example-execution/virtuals/",
    use_fk_colors=False,
    warn_classical=True,
):
    """
    Plot a single step in the execution of the algorithm.

    `gc`: Together with `orient`, comprises the sage-format gauss code
    `orient`: see above
    `state`: a single entry in the `states` list defined in `route()`.
    `i`: the current step of the computation. Used to generate the
         filename.
    """
    expanded_dirname = f"{dir_prefix}{dirname}/"
    if not os.path.exists(expanded_dirname):
        os.mkdir(expanded_dirname)

    os.chdir(expanded_dirname)

    fname = f"{dirname}.tex"

    nmax = len(gc) // 2

    out_str = r"""\documentclass[border=1pt]{standalone}
\usepackage{tikz}"""

    if use_fk_colors:
        out_str += r"""\usepackage{xcolor} % For customizing page color and such
\definecolor{bcol}{HTML}{1D252C}
\definecolor{tcol}{HTML}{D6D7D9}
\pagecolor{bcol}
\color{tcol}"""

    out_str += r"""\begin{document}
  \begin{tikzpicture}[every draw/.style={line width=20pt}]
    """
    x_vals = set()
    for key in upper_cs:
        for xtup in upper_cs[key]:
            x_vals |= set(xtup)
    for key in lower_cs:
        for xtup in lower_cs[key]:
            x_vals |= set(xtup)

    # Set the max x value because I didn't make the code modular and
    # it's expecting to draw the vertical line
    x = max(x_vals)

    int_pts = compute_intersections(upper_cs, lower_cs)
    if warn_classical and (not int_pts):
        print(f"knot {dirname} had no virtual crossings!")
    for (x, y) in int_pts:
        out_str += f"\\draw ({x}, {y}) circle (.25);\n"

    seen_cs = set()
    i = 0
    # Only draw the crossings we've seen so far in black
    while len(seen_cs) < nmax:
        c = gc[i]
        ind = abs(c)
        if ind % 1 != 0:
            ind -= 0.5

        if abs(c) not in seen_cs:
            out_str += crossing_tex(c, orient[len(seen_cs)], crossings[len(seen_cs)], x)
            seen_cs.add(abs(c))
        i += 1

    # Now we draw the semicircles
    for s in upper_cs.keys():
        # Special case for if we get the first arc
        if upper_cs[s] == []:
            pairs = lower_cs[s]
        else:
            pairs = upper_cs[s]

        for (x0, x1) in pairs:
            if x0 > x1:
                (x0, x1) = (x1, x0)

            # Check wheter we need to draw the arc down through the horizontal
            if x0 in crossings:
                x0_t = False
            else:
                x0_t = True

            if x1 in crossings:
                x1_t = False
            else:
                x1_t = True

            # If not in the special case, do the regular call
            if upper_cs[s] != []:
                out_str += arc_tex(x0, x1, x, 0.5, s, x0_through=x0_t, x1_through=x1_t)

    out_str += "\\begin{scope}[yscale=-1]\n"
    for s in lower_cs.keys():
        # Special case for if we get the first arc
        if lower_cs[s] == []:
            pairs = upper_cs[s]
        else:
            pairs = lower_cs[s]

        for (x0, x1) in pairs:
            if x0 > x1:
                (x0, x1) = (x1, x0)

            # Check wheter we need to draw the arc down through the horizontal
            if x0 in crossings:
                x0_t = False
            else:
                x0_t = True

            if x1 in crossings:
                x1_t = False
            else:
                x1_t = True

            # If not in the special case, do the regular call
            if lower_cs[s] != []:
                out_str += arc_tex(x0, x1, x, 0.5, s, x0_through=x0_t, x1_through=x1_t)

    out_str += "\\end{scope}"
    out_str += "  \\end{tikzpicture}\n\\end{document}"
    with open(fname, "w") as f:
        f.write(out_str)

    # pdflatex seems to be faster than lua for these files
    subprocess.run(["pdflatex", fname], capture_output=True)
    os.chdir("../../..")


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


def make_mosaic(flavor="virtuals"):
    os.chdir(f"example-execution/{flavor}")
    # Filter for pdfs with a `-` in the name
    knots = [fname for fname in os.listdir() if "-" in fname]

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
    # out_str += r"""\usepackage{xcolor} % For customizing page color and such
    # \definecolor{bcol}{HTML}{1D252C}
    # \definecolor{tcol}{HTML}{D6D7D9}
    # \pagecolor{bcol}
    # \color{tcol}
    # """
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

                out_str += (
                    r"\begin{subfigure}[t]{" + str(spacing) + r"\textwidth}" + "\n"
                )
                out_str += r"\centering" + "\n"
                out_str += (
                    r"\includegraphics[width = .95in]{"
                    + fname
                    + "/"
                    + fname
                    + ".pdf"
                    + "}\n"
                )
                out_str += r"\caption*{$" + f"{label}" + "$}" + "\n"
                out_str += r"\end{subfigure}" + "\n"

        out_str += "\\end{figure}\n\\vfill\n"
    out_str += r"\end{document}"

    with open(f"{flavor}_mosaic.tex", "w") as f:
        f.write(out_str)

    # Use lualatex in case the mosaic exceeds pdflatex memory limit
    subprocess.run(["lualatex", f"{flavor}_mosaic.tex"], capture_output=True)
    os.chdir("../..")

    return
