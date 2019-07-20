import numpy as np
from gauss_codes import gknot
import os
import copy
from collections import defaultdict
import itertools as it
import random

def sign(i):
    if i > 0:
        return 1
    elif i < 0:
        return -1
    return 0

def get_cnum(c):
    """
    Extract the crossing number corresponding to an entry in a gauss code.

    Args:

    Crossings in our gauss code are formatted like

    (+/-) (crossing number).(0/0.5)

    where (+) indicates we're on an overstrand, (-) indicates we're on an
    understrand, a trailing (.0) indicates the crossing is positive, and a
    trailing (.5) indicates the crossing is negative.
    """

    if c < 0:
        if c % 1 != 0:
            c += 0.5
    else:
        if c % 1 != 0:
            c -= 0.5

    return int(abs(c))

def get_y_in(cross_seen, c):
    """
    Get the incident y value corresponding to the term
    """
    if not cross_seen:
        return 0

    if c < 0:
        if c % 1 != 0:
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
            c_y = 1
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
            c_y = -1
    else:
        if c % 1 != 0:
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
            c_y = -1
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
            c_y = 1

    return c_y


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


# def get_gap_x(c1, c2, r_gap_c, n):
#     l_gap_c = r_gap_c - 1

#     # In this case, we're entering _in_ through the gap, so we want to hug a
#     # different side.
#     if l_gap_c < c1:
#         x_shift = l_gap_c - c2 + 1
#         return (l_gap_c * (n + 3)) + x_shift
#     # Exiting gap
#     else:
#         x_shift = c2 - r_gap_c + 1
#         return r_gap_c * (n + 3) - x_shift

def get_gap_x(c1, c2, n):
    assert c1 != c2
    if c1 < c2:
        # exiting
        return c2 * 3 - 2
    else:
        # Suspect this never actually occurs
        return c2 * (n + 3) + 1 + n - c1

def find_gap(cross_x, envelope):
    print(envelope)
    assert envelope != None
    _, ex2 = envelope

    # stupid algorithm
    for i, x in enumerate(cross_x):
        if x > ex2:
            return i

    return len(cross_x)

def check_gap(c_shapes, cross_x, gapi):
    gx1 = cross_x[gapi-1]

    # TODO: this is sketch
    if gapi < len(cross_x):
        gx2 = cross_x[gapi]
    else:
        gx2 = cross_x[-1] + cross_x[1]

    # checking for a line segment
    for (cx1, cx2, cy) in c_shapes:
        if cy == 0 and cx1 <= gx1 and gx2 <= cx2:
            return False

    return True

def escape(c_shapes, cross_x, x, y):
    n = len(cross_x)
    e = get_envelope(c_shapes, x, y)
    if e is None:
        return [None]

    gapi = find_gap(cross_x, e)

    if check_gap(c_shapes, cross_x, gapi):
        # calculate c1 from x and c2 from gap
        c1 = x//3
        c2 = gapi
        x2 = get_gap_x(c1, c2, n)
        return [(*e, x2)] + escape(c_shapes, cross_x, x2, -y)

    else:
        return [e]

def get_path(path, cross_x, c1i, c2i, x1, x2, y1, y2):
    n = len(cross_x)
    c_shapes = get_c_shapes(path)
    delta = abs(c2i - c1i)
    if y1 == y2:
        if y1 == 0:
            # CASE 0
            # easy
            path.append((x2, 0))

        else:
            # CASE 1
            e1 = get_envelope(c_shapes, x1, y1)
            e2 = get_envelope(c_shapes, x2, y2)
            if e1 == e2:
                # base case, they share the same envelope and are on
                # the same side, so we know that we can just connect
                # them directly.
                add_c_shape(path, x1, x2, y1)
            else:
                # The exit and approach are on opposite sides and at least one has an containing envelope.

                esc1 = escape(c_shapes, cross_x, x1, y1)
                esc2 = escape(c_shapes, cross_x, x2, y2)

                print(esc1, esc2)
                # We already know that we aren't in a matching
                # envelope, so we have two cases to consider.
                # TODO: Validate that this is correct
                # TODO: This is the same as case 4
                if len(esc1) >= len(esc2):
                    # we break symmetry arbitrarily
                    x_psuedo = esc1[0][2]
                    print('recurse!')
                    add_c_shape(path, x1, x_psuedo, y1)
                    get_path(path, cross_x, find_gap(cross_x, (0, x_psuedo)), c2i, x_psuedo, x2, -y1, y2)

                else:
                    x_psuedo = esc2[0][2]
                    # Ok, so this time we go in reverse
                    print('recurse!')
                    get_path(path, cross_x, c1i, find_gap(cross_x, (0, x_psuedo)), x1, x_psuedo, y1, -y2)
                    add_c_shape(path, x_psuedo, x2, y2)

                print('oh no1')

    else:
        if y1 == 0:
            # CASE 2
            # We are coming out of the first crossing horizontally,
            # but into the second vertically, so we must be
            # backtracking.

            path.append((x1 + 1, 0))
            e1 = get_envelope(c_shapes, x1, y2)
            e2 = get_envelope(c_shapes, x2, y2)
            if e1 == e2:
                # base case
                add_c_shape(path, x1+1, x2, y2)

            else:
                # We're backtracking, so the source must be "open to
                # the air". In fact, it must be open on both sides of
                # the main line.
                assert escape(c_shapes, cross_x, x1 + 1, -1) == [None]
                assert escape(c_shapes, cross_x, x1 + 1, +1) == [None]

                esc2 = escape(c_shapes, cross_x, x2, 1)
                print(esc2)

                x_psuedo = esc2[0][2]
                print('recurse!')
                get_path(path, cross_x, c1i, find_gap(cross_x, (0, x_psuedo)), x1, x_psuedo, y1, -y2)
                add_c_shape(path, x_psuedo, x2, y2)
                print('oh no2')

        elif y2 == 0:
            # CASE 3
            # We want to go into the second crossing horizontally

            # We go just short of the path, then we go in.
            # Importantly, there is no restriction on which direction
            # we can approach the main line from.
            e1 = get_envelope(c_shapes, x1, y1)
            e2 = get_envelope(c_shapes, x2 - 1, y1)
            if e1 == e2:
                # base case
                add_c_shape(path, x1, x2 - 1, y1)

            else:
                # Two cases to consider here. Either the destination
                # is a new vertex, or it's 1, but in either case the
                # destination is "open to the air", although we don't
                # know which direction is open. Therefore, all of the
                # escaping is done on the source end.
                esc2u = escape(c_shapes, cross_x, x2-1, +1)
                esc2d = escape(c_shapes, cross_x, x2-1, -1)
                assert esc2u == [None] or esc2d == [None]

                esc1 = escape(c_shapes, cross_x, x1, y1)

                # Thus, we escape from the source
                x_psuedo = esc1[0][2]
                add_c_shape(path, x1, x_psuedo, y1)
                print("recurse!")
                get_path(path, cross_x, find_gap(cross_x, (0, x_psuedo)), c2i, x_psuedo, x2-1, -y1, y2)

                print('oh no3')

            path.append((x2, 0))

        else:
            # CASE 4
            # We must exit and approach from different sides, so we
            # must cross the main line at some point. The envelopes
            # cannot be equal unless they are both the None (aka.
            # largest) envelope.

            e1 = get_envelope(c_shapes, x1, y1)
            e2 = get_envelope(c_shapes, x2, y2)

            if e1 == e2 == None:
                # If both are "open to the air", we can go to the
                # right or left, we choose to go to the right for
                # aesthetic reasons. There can be no more than min(n -
                # c1i, n - c2i) other paths that must go around the
                # right, and we add 1 to make room for the horizontal
                # exit from the last crossing.

                x_psuedo = cross_x[-1] + 1 + min(n - c1i, n - c2i)
                add_c_shape(path, x1, x_psuedo, y1)
                add_c_shape(path, x_psuedo, x2, y2)
            else:
                # The exit and approach are on opposite sides and at least one has an containing envelope.

                esc1 = escape(c_shapes, cross_x, x1, y1)
                esc2 = escape(c_shapes, cross_x, x2, y2)

                print(esc1, esc2)
                # We already know that we aren't in a matching
                # envelope, so we have two cases to consider.
                # TODO: Validate that this is correct
                if len(esc1) >= len(esc2):
                    # we break symmetry arbitrarily
                    x_psuedo = esc1[0][2]
                    print('recurse!')
                    add_c_shape(path, x1, x_psuedo, y1)
                    get_path(path, cross_x, find_gap(cross_x, (0, x_psuedo)), c2i, x_psuedo, x2, -y1, y2)

                else:
                    x_psuedo = esc2[0][2]
                    # Ok, so this time we go in reverse
                    print('recurse!')
                    get_path(path, cross_x, c1i, find_gap(cross_x, (0, x_psuedo)), x1, x_psuedo, y1, -y2)
                    add_c_shape(path, x_psuedo, x2, y2)

def find_x_conflicts(bad_path, n):
    # First, get rid of duplicate points in the path
    path = [bad_path[0]]
    for p in bad_path:
        if path[-1] != p:
            path.append(p)

    x_lookup = defaultdict(lambda: [])
    for p0, p1, p2, p3, p4 in zip(path, path[1:], path[2:], path[3:], path[4:]):
        assert p0 != p1
        (p1x, p1y), (p2x, p2y), (p3x, p3y) = p1, p2, p3
        if not (p1x == p2x == p3x): continue
        # We're going across the main line
        # we know the segment is vertical

        # TODO: prove we can't get a conflict at a crossing
        if p1x%3 == 0: continue
        (p0x, p0y), (p4x, p4y) = p0, p4

        assert p0x != p1x
        assert p4x != p3x
        assert p0y == p1y
        assert p4y == p3y

        # we know the end segments are horizontal
        d0 = -1 if p0x < p1x else 1
        d4 = -1 if p4x < p3x else 1

        x_lookup[p1x].append(
            (p1y, p3y, d0, d4)
            # (p0, p1, p2, p3, p4)
        )

    import pprint
    pprint.pprint(dict(x_lookup.items()))

    nudge_map = {}
    for x, shapes in x_lookup.items():
        for p in it.permutations(shapes):
            fail = False
            for p1, p2 in zip(p, p[1:]):
                p1y1, p1y2, d1x1, d1x2 = p1
                p2y1, p2y2, d2x1, d2x2 = p2

                # Orient the flips
                if p1y1 > p1y2:
                    p1y1, p1y2 = p1y2, p1y1
                    d1x1, d1x2 = d1x2, d1x1

                if p2y1 > p2y2:
                    p2y1, p2y2 = p2y2, p2y1
                    d2x1, d2x2 = d2x2, d2x1

                # Detect crosses
                if p1y1 > p2y1 and d1x1 == 1: # p1y1 above p2y1 crosses it forward
                    break
                if p1y1 < p2y1 and d2x1 == -1: # p2y1 above p1y1 cross it backward
                    break

                if p1y2 < p2y2 and d1x2 == 1: # p1y2 below p2y2 crosses it forward
                    break
                if p1y2 > p2y2 and d2x2 == -1: # p2y2 below p1y2 cross it backward
                    break
            else:
                # this ordering is valid
                break
        else:
            assert False

        for i, shape in enumerate(p):
            nudge_map[shape[:2]] = (i+1)/(len(p) + 1)

    print(nudge_map)

    for i in range(len(path) - 4):
        p0, p1, p2, p3, p4 = path[i:i+5]
        (p1x, p1y), (p2x, p2y), (p3x, p3y) = p1, p2, p3
        if not (p1x == p2x == p3x): continue
        if not p2y == 0: continue
        # rounding is a horrible hack
        if (round(p1y), round(p3y)) in nudge_map:
            dx = nudge_map[(round(p1y), round(p3y))]
            (p0x, p0y), (p4x, p4y) = p0, p4
            # TODO: Idk if this is quite right
            dy1 = dx/2 if (p0x < p1x) == (p1y < 0) else -dx/2
            dy3 = dx/2 if (p4x < p3x) == (p1y < 0) else -dx/2
            p1, p2, p3 = (p1x+dx, p1y+dy1), (p2x+dx, p2y), (p3x+dx, p3y + dy3)
            p0 = (p0x, p0y + dy1)
            p4 = (p4x, p4y + dy3)
            path[i:i+5] = p0, p1, p2, p3, p4
            print(p0, p1, p2, p3, p4)

    return path


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
    cross_x = [i*3 for i in range(n)]

    # Keep track of crossings we've seen already.
    cross_seen = [False for i in range(n)]

    for i in range(len(gcode) - 1):
        # Get information from the gauss code
        c1, c2 = gcode[i], gcode[i+1]
        c1i = get_cnum(c1) - 1
        c2i = get_cnum(c2) - 1
        x1, x2 = cross_x[c1i], cross_x[c2i]

        # y_out is -1 * y_in
        y1 = -1 * get_y_in(cross_seen[c1i], c1)
        y2 = get_y_in(cross_seen[c2i], c2)
        cross_seen[c1i] = True

        print(f"c: {c1} → {c2}, x: {x1} → {x2}, y: {y1} → {y2}")

        get_path(path, cross_x, c1i, c2i, x1, x2, y1, y2)

    # fix the termination FIXME!
    x1, _ = path[-1]
    y1 = -1 * get_y_in(True, gcode[-1])
    print(f"c: {gcode[-1]} → {gcode[0]}, x: {x1} → {0}, y: {y1} → {0}")
    get_path(path, cross_x, 0, 0, x1, 0, y1, 0)

    signs = []
    for cnum in range(1,n+1):
        # This is gross
        for c in gcode:
            stop = False
            if get_cnum(c) == cnum:
                signs.append((c/abs(c)))
                stop = True
            if stop:
                break

    path = find_x_conflicts(path, n)

    return path, cross_x, signs


def draw_presentation(draw_paths, x_vals, signs, fname="test_gauss"):
    """
    draw

    draw_paths a list of lists where each sublist is a list of tuples specifying
    a path to draw in TikZ

    x_vals is
    """
    break_width = 0.4

    # Get the LaTeX preamble
    preamble = "\\documentclass[border=10pt]{standalone}\n"
    preamble += "\\usepackage{tikz}\n"
    preamble += "\\begin{document}\n"
    preamble += "\\begin{tikzpicture}\n"

    drawing = ""
    for path in draw_paths:
        drawing += "\\draw[-latex, line width=1mm]"
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
            drawing += f"\\draw[line width=1mm] ({x}, {break_width}) -- ({x}, -{break_width});\n"
        else:
            drawing += f"\\draw[line width=1mm] ({x-break_width}, 0) -- ({x+break_width}, 0);\n"

        drawing += f"\\node[circle, fill=black, draw=black, inner sep=3pt] () at ({x}, 0) " + "{};\n"
        drawing += f"\\node[above left] () at ({x}, 0) " + "{\\large $" + str(int(cnum+1)) + "$};\n"

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
    knot = gknot[(8, 4)]
    # knot = gknot[(3, 1)]
    # knot = gknot[(6,4)]
    knot = gknot[(8,16)]
    # knot = gknot[(5,2)]
    # knot = gknot[(4,1)]
    path, cross_x, signs = build_stupid_graph(knot)

    # pathological_test = [
    #     -1.5, 2, -3, 4.5, -5.5, 3, -6, 1.5, 7, 5.5, -4.5, 6, -2, -7
    # ]
    # path, cross_x, signs = build_stupid_graph(pathological_test)

    draw_presentation([path], cross_x, signs)
