# ============================================================================ #
#                                                                              #
#                              Linear Presentation                             #
#                                                                              #
# ============================================================================ #
#
#
# ---------------------------------------------------------------------------- #
#                                                                              #
#                     ┌───────────┐                                            #
#                     │           │                                            #
#        ┌───┐        │   ┌───┐   │                                            #
#        │   │        │   │   │   │                                            #
#    ┌───────┘    ┌───│───────│───┘                                            #
#    │   │        │   │   │   │                                                #
#    └───┘        │   └───┘   │                                                #
#                 │           │                                                #
#                 └───────────┘                                                #
#                                                                              #
# ---------------------------------------------------------------------------- #
#
# If you want to build your own, try these
# ┌ └ ┐ ┘ ├ ┤ ┬ ┴ ┼ │ ─

import itertools as it
from gauss_codes import gknot, conn_sum
import networkx as nx
import os


################################################################################
#                                                                              #
#                               Helper Functions                               #
#                                                                              #
################################################################################


def rotate_right(l, n):
    return l[-n:] + l[:-n]


def irotate_right(l, n):
    return it.chain(l[-n:], l[:-n])


def sign(i):
    """
    FIXME why is there a return 0?
    """
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


def decompose_c(c):
    """
    In our initial path construction, routes can take many shapes:

         ┌──┐           ┌───────┐
         │  │           │       │
    x    │  y           x       y
    │    │
    └────┘

    Note, we can decompose all of these paths into c shapes

    """
    over_under = 1 if c > 0 else -1
    c = abs(c)
    pos_neg = 1 if c % 1 == 0 else -1
    cnum = int(abs(c))
    return cnum, over_under, pos_neg


def recompose_c(cnum, over_under, pos_neg):
    if pos_neg == -1:
        cnum += 0.5
    if over_under == -1:
        cnum *= -1
    return cnum


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


def build_embedding_graph(gcode):
    n = len(gcode) // 2

    G = nx.Graph()

    crossings = []
    knot_order = []
    # First, build crossing gadgets
    for i in range(n):
        nodes = range(5 * i, 5 * i + 5)
        G.add_nodes_from(nodes)

        t, l, c, r, b = nodes
        G.add_edge(l, t)
        G.add_edge(r, t)
        G.add_edge(b, l)
        G.add_edge(b, r)

        crossings.append(nodes)

    for i in range(len(gcode)):
        # Get information from the gauss code
        cross1, o1, p1 = decompose_c(gcode[i])
        cross2, o2, p2 = decompose_c(gcode[(i + 1) % len(gcode)])
        print(cross1, cross2, o1, o2, p1, p2)

        t1, l1, c1, r1, b1 = crossings[cross1 - 1]
        t2, l2, c2, r2, b2 = crossings[cross2 - 1]

        if o1 == 1:
            src = t1
        elif p1 == 1:
            src = l1
        else:
            src = r1

        if o2 == 1:
            dst = b2
        elif p2 == 1:
            dst = r2
        else:
            dst = l2

        G.add_edge(c1, src)
        G.add_edge(src, dst)
        G.add_edge(dst, c2)
        knot_order.append([c1, src, dst])

    assert nx.check_planarity(G)[0]

    print(crossings)

    return G, knot_order


def walk_face(embed, src, dst):
    face = [src, dst]

    while True:
        last, new = embed.next_face_half_edge(*face[-2:])
        if [last, new] == face[:2]:
            break
        face.append(new)

    return face[:-1]


def get_all_faces(G, embed):
    faces = []
    half_edges = set(it.chain(G.edges, map(lambda e: e[::-1], G.edges)))
    while half_edges:
        e = next(iter(half_edges))
        face = walk_face(embed, *e)
        faces.append(face)
        for fe in zip(face, irotate_right(face, -1)):
            half_edges.remove(fe)

    return faces


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
    min_size, min_c = float("inf"), None
    for (cx1, cx2, cy) in c_shapes:
        if cy == y and cx1 < x and x < cx2:
            size = cx2 - cx1
            assert size != min_size
            if size < min_size:
                min_size = size
                min_c = (cx1, cx2)

    return min_c


def find_gap(cross_x, path, envelope):
    print("finding gap for", envelope)
    assert envelope != None
    _, ex2 = envelope
    c_shapes = get_c_shapes(path)

    # For sanity checking, we ensure that the crossing is valid by
    # checking for horizontal lines between the two nearest crossings.

    # TODO: This is sketch, we are advancing past the last crossing
    # and back from the first by cross_x[1], which is supposed to be
    # the delta.
    next_x = min((x for x in cross_x if x > ex2), default=cross_x[-1] + cross_x[1])
    prev_x = max((x for x in cross_x if x < ex2), default=-cross_x[1])

    # Check for the connecting line segment on the main line
    for (cx1, cx2, cy) in c_shapes:
        if cy == 0 and prev_x + 1 <= cx1 <= cx2 <= next_x - 1:
            print(f"blocking segment {cx1}, {cx2} found inside {prev_x}, {next_x}")
            return [False, -1, -1, -1]

    # Ok the gap exists, so where is it? ex2 is the right side of the
    # envelope, so it's the right bound. We need to search for the
    # left bound.
    right_bound = ex2
    left_bound = max((x for (x, y) in path if x < right_bound), default=-cross_x[1])
    midpoint = (right_bound + left_bound) / 2

    print(f"gap from {left_bound} to {right_bound}, midpoint {midpoint}")
    return True, left_bound, midpoint, right_bound


def escape(cross_x, path, x, y):
    print(f"Escaping from {x} at y = {y}")
    c_shapes = get_c_shapes(path)
    e = get_envelope(c_shapes, x, y)
    if e is None:
        return [None]

    gap_found, _, midpoint, _ = find_gap(cross_x, path, e)

    if gap_found:
        return [(*e, midpoint)] + escape(cross_x, path, midpoint, -y)

    return [e]


def normalize_gauss_order(gcode):
    crossings_map = {}
    top_cross = 0
    new_gcode = []
    for c in gcode:
        cnum, over_under, pos_neg = decompose_c(c)
        if cnum not in crossings_map:
            top_cross += 1
            crossings_map[cnum] = top_cross

        new_gcode.append(recompose_c(crossings_map[cnum], over_under, pos_neg))

    return new_gcode


def get_path(path, cross_x, x1, x2, y1, y2):
    print(f"Getting path from {x1} to {x2}")
    n = len(cross_x)
    c_shapes = get_c_shapes(path)
    if y1 == y2:
        if y1 == 0:
            # CASE 0
            # easy
            path.append((x2, y2))
            # TODO: This case may not be as easy as we though
            assert get_envelope(c_shapes, x1 + 1, 1) == None
            assert get_envelope(c_shapes, x2 - 1, 1) == None
            assert get_envelope(c_shapes, x1 + 1, -1) == None
            assert get_envelope(c_shapes, x2 - 1, -1) == None

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

                esc1 = escape(cross_x, path, x1, y1)
                esc2 = escape(cross_x, path, x2, y2)

                print(esc1, esc2)
                # We already know that we aren't in a matching
                # envelope, so we have two cases to consider.
                # TODO: Validate that this is correct
                # TODO: This is the same as case 4
                if len(esc1) >= len(esc2):
                    # we break symmetry arbitrarily
                    x_psuedo = esc1[0][2]
                    print("recurse! 1.1")
                    add_c_shape(path, x1, x_psuedo, y1)
                    get_path(path, cross_x, x_psuedo, x2, -y1, y2)

                else:
                    x_psuedo = esc2[0][2]
                    # Ok, so this time we go in reverse
                    print("recurse! 1.2")
                    get_path(path, cross_x, x1, x_psuedo, y1, -y2)
                    add_c_shape(path, x_psuedo, x2, y2)

                print("oh no1")

    else:
        if y1 == 0:
            # CASE 2
            # We are coming out of the first crossing horizontally,
            # but into the second vertically, so we must be
            # backtracking.

            # The question is, which way do we go, up or down?

            # First, we need to compute the orientation of the gadgets
            # edge_order = [0, 3, 4, 1]
            # src_rotation_map = {(+1, +1): 1, (+1, -1): 1, (-1, +1): 2, (-1, -1): 0}
            # src_edge_order = rotate_right(edge_order, src_rotation_map[decompose_c(c)])

            ## Old implementation here
            path.append((x1 + 1, 0))
            e1 = get_envelope(c_shapes, x1, y2)
            e2 = get_envelope(c_shapes, x2, y2)
            if e1 == e2:
                # base case
                add_c_shape(path, x1 + 1, x2, y2)

            else:
                # We're backtracking, so the source must be "open to
                # the air". In fact, it must be open on both sides of
                # the main line.
                assert escape(cross_x, path, x1 + 1, -1) == [None]
                assert escape(cross_x, path, x1 + 1, +1) == [None]

                esc2 = escape(cross_x, path, x2, y2)
                print(esc2)

                x_psuedo = esc2[0][2]
                print(f"recurse! 2 {x_psuedo}")
                get_path(path, cross_x, x1, x_psuedo, y1, -y2)
                add_c_shape(path, x_psuedo, x2, y2)
                print("oh no2")

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
                esc2u = escape(cross_x, path, x2 - 1, +1)
                esc2d = escape(cross_x, path, x2 - 1, -1)
                assert esc2u == [None] or esc2d == [None]

                esc1 = escape(cross_x, path, x1, y1)

                # Thus, we escape from the source
                x_psuedo = esc1[0][2]
                add_c_shape(path, x1, x_psuedo, y1)
                print("recurse! 3")
                get_path(path, cross_x, x_psuedo, x2 - 1, -y1, y2)

                print("oh no3")

            path.append((x2, y2))

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

                # x_psuedo = cross_x[-1] + 1 + min(n - c1i, n - c2i)

                x_max = cross_x[-1] + cross_x[1]
                x_min = max(x for (x, y) in path)
                x_psuedo = (x_max + x_min) / 2

                add_c_shape(path, x1, x_psuedo, y1)
                add_c_shape(path, x_psuedo, x2, y2)
            else:
                # The exit and approach are on opposite sides and at least one has an containing envelope.

                esc1 = escape(cross_x, path, x1, y1)
                esc2 = escape(cross_x, path, x2, y2)

                print(esc1, esc2)
                # We already know that we aren't in a matching
                # envelope, so we have two cases to consider.
                # TODO: Validate that this is correct
                if len(esc1) >= len(esc2):
                    # we break symmetry arbitrarily
                    x_psuedo = esc1[0][2]
                    print("recurse! 4.1")
                    add_c_shape(path, x1, x_psuedo, y1)
                    get_path(path, cross_x, x_psuedo, x2, -y1, y2)

                else:
                    x_psuedo = esc2[0][2]
                    # Ok, so this time we go in reverse
                    print("recurse! 4.2")
                    get_path(path, cross_x, x1, x_psuedo, y1, -y2)
                    add_c_shape(path, x_psuedo, x2, y2)


def build_stupid_graph(gcode):
    norm_gcode = normalize_gauss_order(gcode)
    if gcode != norm_gcode:
        print(f"Reordered\n{gcode} to\n{norm_gcode}.")
        gcode = norm_gcode

    # The total number of crossings (n) is going to be 1/2 the total length of
    # the gauss code. We use integer divide to implicitly cast to an int so that
    # we can put it into `range(n)` later.
    n = len(gcode) // 2

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
    cross_x = [i * 3 for i in range(n)]

    # Keep track of crossings we've seen already.
    cross_seen = [False for i in range(n)]

    G, _ = build_embedding_graph(gcode)
    embed = nx.check_planarity(G)[1]
    faces = get_all_faces(G, embed)

    for i in range(len(gcode) - 1):
        # Get information from the gauss code
        c1, c2 = gcode[i], gcode[i + 1]
        c1i = get_cnum(c1) - 1
        c2i = get_cnum(c2) - 1
        x1, x2 = cross_x[c1i], cross_x[c2i]

        # y_out is -1 * y_in
        y1 = -1 * get_y_in(cross_seen[c1i], c1)
        y2 = get_y_in(cross_seen[c2i], c2)
        cross_seen[c1i] = True

        print(f"c: {c1} → {c2}, x: {x1} → {x2}, y: {y1} → {y2}")
        get_path(path, cross_x, x1, x2, y1, y2)
    else:
        # fix the termination FIXME!
        x1, _ = path[-1]
        y1 = -1 * get_y_in(True, gcode[-1])
        print(f"c: {gcode[-1]} → {gcode[0]}, x: {x1} → {0}, y: {y1} → {0}")
        get_path(path, cross_x, x1, 0, y1, 0)

    signs = []
    for cnum in range(1, n + 1):
        # This is gross
        for c in gcode:
            stop = False
            if get_cnum(c) == cnum:
                signs.append((c / abs(c)))
                stop = True
            if stop:
                break

    return path, cross_x, signs


def space_xy(draw_path, cross_x):
    """
    gets it good enough to fix by hand
    """
    x_vals = set()
    y_vals = set()
    for x, y in draw_path:
        x_vals.add((x))
        y_vals.add((y))

    x_vals = sorted(x_vals)
    y_vals = sorted(y_vals)

    x0 = x_vals.index(0)
    y0 = y_vals.index(0)

    x_map = {x: i - x0 for i, x in enumerate(x_vals)}
    y_map = {y: i - y0 for i, y in enumerate(y_vals)}

    new_path = []
    for (x, y) in draw_path:
        new_path += [(x_map[x], y_map[y])]

    new_cross_x = [x_map[x] for x in cross_x]

    return new_path, new_cross_x


def draw_presentation(draw_paths, x_vals, signs, fname="test_gauss", display=True):
    """
    draw

    draw_paths a list of lists where each sublist is a list of tuples specifying
    a path to draw in TikZ

    x_vals is
    """
    n = len(x_vals)
    break_width = 0.2

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
            drawing += f"\\draw[] ({x}, {break_width}) -- ({x}, -{break_width});\n"
        else:
            drawing += f"\\draw[] ({x-break_width}, 0) -- ({x+break_width}, 0);\n"

        drawing += (
            f"\\node[circle, fill=black, draw=black, inner sep=1pt] () at ({x}, 0) "
            + "{};\n"
        )
        drawing += (
            f"\\node[above left] () at ({x}, 0) " + "{$" + str(int(cnum + 1)) + "$};\n"
        )

    out_str = preamble + drawing + "\\end{tikzpicture}\n\\end{document}"

    # Need to chdir else pdflatex will pollute the parent directory with aux
    # files and stuff
    os.chdir("tests")
    with open(f"{fname}.tex", "w") as f:
        f.write(out_str)
    os.system(f"pdflatex -halt-on-error {fname}.tex | grep '^!.*' -A200 --color=always")
    print(fname)
    if display:
        os.system(f"zathura {fname}.pdf")
    os.chdir("..")
    return out_str


def aggressive_tightener(path, cross_x):

    return


def get_shift_dict(gcode):
    """
    Create a dict where <val> indicates how many unique crossings we'll
    encounter before backtracking if we perform <key> cyclic permutations of the
    input Gauss code
    """
    n = len(gcode) // 2
    len_dict = {i: 0 for i in range(n)}

    for i in len_dict:
        i_gcode = gcode[i:] + gcode[:i]
        i_seen_cs = set()
        for c in i_gcode:
            cnum = get_cnum(c)
            if cnum in i_seen_cs:
                break
            else:
                len_dict[i] += 1
                i_seen_cs.add(cnum)

    sorted_dict = sorted(len_dict.items(), key=lambda kv: kv[1])
    print("sorted_dict is", sorted_dict)
    return sorted_dict


def construct_all():
    for gc in gknot.keys():
        knot = gknot[gc]
        sorted_dict = get_shift_dict(knot)

        # We want to get all cyclic shift permutations as well
        for shift_num in range(len(knot) // 2):
            print(knot)
            path, cross_x, signs = build_stupid_graph(knot)
            path, cross_x = space_xy(path, cross_x)
            try:
                c, i = gc
                draw_presentation([path], cross_x, signs, fname=f"{c}_{i}")
            except TypeError:
                draw_presentation([path], cross_x, signs, fname="0")

            knot = knot[1:] + [knot[0]]


if __name__ == "__main__":
    # construct_all()

    # Bad apples currently:
    # knot = gknot[(6,2)]
    # knot = gknot[(7,2)]
    # knot = gknot[(6,6)]
    # knot = gknot[(7,3)]
    # knot = gknot[(8,10)]

    # knot_inds_to_sum = [(8,10), (8,10), (8,10), (8,10)]
    # inds = [4, 1, 2, 7]
    # knot = gknot[knot_inds_to_sum.pop()]
    # for i, ind in zip(inds, knot_inds_to_sum):
    #     knot = conn_sum(knot, gknot[ind], ind=i)

    # path, cross_x, signs = build_stupid_graph(knot)
    # path, cross_x = space_xy(path, cross_x)
    # print(path)
    # draw_presentation([path], cross_x, signs, fname="test_gauss")

    knot = gknot[(11, 2)]
    # knot = gknot[(11, 42)]
    # knot = gknot[(11,2)]

    path, cross_x, signs = build_stupid_graph(knot)
    # path, cross_x = space_xy(path, cross_x)
    draw_presentation([path], cross_x, signs, fname="0")
    # ...and probably more, but that's where we get wrecked rn.

    # pathological_test = [
    #     -1.5, 2, -3, 4.5, -5.5, 3, -6, 1.5, 7, 5.5, -4.5, 6, -2, -7
    # ]
    # path, cross_x, signs = build_stupid_graph(pathological_test)
    # path, cross_x = space_xy(path, cross_x)
    # draw_presentation([path], cross_x, signs, fname="0")
