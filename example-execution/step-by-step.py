from sage.all import *
from enum import Enum
import itertools as it
import matplotlib.pyplot as plt
import matplotlib as mpl

# For the step-by-step plotting garbage
from copy import deepcopy
from os import mkdir, listdir, chdir
from os.path import isdir

import gauss_codes

# For the parallel compilation
from multiprocessing import Pool

# For executing the compile command
from subprocess import run


class Dir(Enum):
    UP = 0
    RIGHT = 1
    DOWN = 2
    LEFT = 3

    def flip(self):
        """
        Perform one of
              UP <-> DOWN
            LEFT <-> RIGHT
        depending on the current value.
        """
        return Dir((self.value + 2) % 4)


def fix_gc_order(gc, orient):
    """
    If there exist i,j with i < j such that crossing j is encountered before
    crossing i, we say the gauss code is improperly ordered.

    This function relabels the gauss code to fix that.
    """
    new_gc = []
    mapping = {}  # Initialize the mapping from "old" -> "new" labels

    for elt in gc:
        c, s = abs(elt), sign(elt)  # "sign" gives over/under info
        if c not in mapping:
            mapping[c] = len(mapping) + 1  # relabel by order of encounter
        new_c = mapping[c]
        new_gc.append(new_c * s)

    rev_mapping = {v: k for k, v in mapping.items()}
    assert len(rev_mapping) == len(mapping)
    new_orient = [orient[rev_mapping[c] - 1] for c in range(1, len(mapping) + 1)]

    return new_gc, new_orient


def linear_layout(gc, orient):
    """
    Given the sage-format gauss code (gc, orient) as input, construct and return

      1. `crossings`, a dictionary indexed by crossing number with values an
         associated "connection information" dictionary. Basically, the
         dictionary for crossing `i` has

           a) keys: the [RIGHT, UP, LEFT, DOWN] strands involved in crossing `i`

           b) values: the other endpoint of said strand, represented as a tuple
              (j, DIR), where j is the crossing number of the endpoint, and DIR
              describes which of [RIGHT, UP, LEFT, DOWN] this strand is for
              crossing j.

      2. `semiarcs`, a list of tuples of the form (i, DIR_i, j, DIR_j),
         representing all of the oriented semiarcs appearing in our diagram.
    """
    n = max(gc)
    crossings = {i: {} for i in range(1, n + 1)}
    semiarcs = []

    # Add crossing 1 and the outgoing direction manually.
    seen = set([1])
    exit_dir = Dir.RIGHT

    # Iterate through all subsequent pairs of characters in the gauss code,
    # which essentially gets us all semiarcs directly.
    for a, b in zip(gc, gc[1:]):
        # The way the algorithm works, whenever we see a new crossing for the
        # first time, we always orient the incoming strand to lie on the
        # horizontal line. This corresponds to coming in from the LEFT.
        if abs(b) not in seen:
            entry_dir = Dir.LEFT

        # Otherwise, we check to see whether
        #  - `b` is positive and we're entering on an overstrand, or
        #  - `b` is negative and we're entering on an understarnd.
        # In both these cases, the entering strand comes in from the top.
        elif sign(b) == orient[abs(b) - 1]:
            entry_dir = Dir.UP

        # Only remaining case is entering from the bottom
        else:
            entry_dir = Dir.DOWN

        crossings[abs(a)][exit_dir] = (abs(b), entry_dir)
        crossings[abs(b)][entry_dir] = (abs(a), exit_dir)
        semiarcs.append((abs(a), exit_dir, abs(b), entry_dir))
        seen.add(abs(a))
        exit_dir = entry_dir.flip()  # Exit opposite of the way we came in

    # Manually add the incoming connection for crossing 1
    crossings[abs(b)][exit_dir] = (1, Dir.LEFT)
    crossings[1][Dir.LEFT] = (abs(b), exit_dir)
    semiarcs.append((abs(b), exit_dir, abs(1), Dir.LEFT))

    return crossings, semiarcs


def unpack_paths(paths, n, bound, m):
    result = []
    for i in range(n):
        path = [paths[i, j] for j in range(m) if paths[i, j] > 0]
        result.append(path)
    return result


def compact_paths(paths, crossings=[]):
    all_points = []
    for path in paths:
        all_points.extend(path)

    all_points = sorted(set(all_points))
    result = []
    for path in paths:
        new_semiarc = [all_points.index(x) for x in path]
        result.append(new_semiarc)

    new_crossings = [all_points.index(x) for x in crossings]
    return result, new_crossings


def plot(upper_cs, lower_cs, crossings=[], straight=0):
    def c_height(cs):
        all_arcs = []
        for arcs in cs.values():
            all_arcs.extend(arcs)

        dominate = {}

        for a, b in all_arcs:
            for c, d in all_arcs:
                if a < c < d < b:
                    dominate.setdefault((a, b), []).append((c, d))

        height = {}
        for a, b in sorted(all_arcs, key=lambda t: len(dominate.get(t, []))):
            height[a, b] = 0
            for dominated in dominate.get((a, b), []):
                height[a, b] = max(height[a, b], height[dominated])

            if b - a > 1:
                height[a, b] += 1

        return height

    upper_heights = c_height(upper_cs)
    lower_heights = c_height(lower_cs)

    def add_arc(a, b, y=0, down=False):
        h = abs(a - b)
        if straight >= 3:
            if down:
                h = lower_heights[a, b]
            else:
                h = upper_heights[a, b]

        if straight >= 1 and abs(a - b) <= 1:
            plt.plot([a, b], [0, 0], "k")

        elif straight >= 2:
            if down:
                plt.plot([a, a, b, b], [0, -h, -h, 0], "k")
            else:
                plt.plot([a, a, b, b], [0, h, h, 0], "k")

        else:
            c = (a + b) / 2
            if down:
                plt.gca().add_patch(mpl.patches.Arc((c, y), h, h, theta1=180))
            else:
                plt.gca().add_patch(mpl.patches.Arc((c, y), h, h, theta2=180))

    for arcs in upper_cs.values():
        for (a, b) in arcs:
            add_arc(a, b, 0, False)

    for arcs in lower_cs.values():
        for (a, b) in arcs:
            add_arc(a, b, 0, True)

    # for i, path in enumerate(paths):
    #     par = parity[i]
    #     for j, (a, b) in enumerate(zip(path, path[1:])):
    #         add_arc(a, b, 0, (par + j) % 2 == 0)

    for x in crossings:
        plt.plot([x - 1, x + 1], [0, 0], "k")

    plt.gcf().set_size_inches(10.5, 18.5)
    ax = plt.gcf().gca()
    ax.set_aspect("equal", "box")
    ax.set_axis_off()
    plt.show()


def nelson_gc_to_sage_gc(gc):
    """
    Convert the gauss code format given in `gauss_codes.py` (Prof. Nelson's
    gauss codes) to the format used internally by `sage`.

    Prof. Nelson's gauss codes use `sign` to indicate under/over, and a trailing
    `.5` to indicate negative crossings.

    Sage's internal Gauss codes separate the crossing sign into a second list.
    """
    new_gc = [int(x) for x in gc]  # Remove all trailing .5's
    assert len(new_gc) % 2 == 0
    n = len(new_gc) // 2  # Number of crossings
    new_orient = [1 if i in gc else -1 for i in range(1, n + 1)]
    return [[new_gc], new_orient]


def knot_to_layout(K):
    """
    Fix the gauss code ordering if needed, and call to `linear_layout()` to get
    the connection information for the crossings, and the list of semiarcs.
    """
    [gc], orient = K.oriented_gauss_code()
    gc, orient = fix_gc_order(gc, orient)
    crossings, semiarcs = linear_layout(gc, orient)
    return crossings, semiarcs


def plain_semiarcs(semiarcs):
    return [(a, b.value, c, d.value) for a, b, c, d in semiarcs]


def route(semiarcs):
    """
    Perform the actual algorithm itself.
    """

    # Initialize the stacks representing the upper / lower half-planes. We'll be
    # pushing / popping semiarcs IDs (see: semiarc_map) to these stacks as we
    # perform the algorithm.
    upper, lower = [], []
    assert len(semiarcs) % 2 == 0

    # Create a dictionary for associating (crossing, dir) pairs with
    # an abstract ID representing the semiarc
    semiarc_map = {}

    # `a` = crossing a, `da` = exit  direction for `a`
    # `b` = crossing b, `db` = enter direction for `b`
    for i, (a, da, b, db) in enumerate(semiarcs):
        semiarc_map[a, da] = semiarc_map[b, db] = i

    # Check if the upper element is v
    def peek(l, v):
        return len(l) >= 1 and l[-1] == v

    # Initialize the dictionaries giving our semicircle positioning information.
    #
    # Keys: semiarc id's
    #
    # Values: a list of numbers. A protypical element x in this list encodes
    #         information as follows:
    #         - abs(x) gives the x coordinate of one of the endpoints of a
    #           semicircle that's part of our arc of interest
    #         - sign(x) describes whether we have an arc that's going out from
    #           the horizontal line or returning to it. x < 0 --> returning,
    #           x > 0 --> leaving.
    upper_cs = {}
    lower_cs = {}

    x = 1
    global states
    states = []

    # In the below, in general, pushing a semiarc to stack corresponds to
    # recording an x coordinate where it leaves the spine. Popping corresponds
    # to recording an x coordinate where it returns to the spine. Note, x will
    # be a variable in the outer scope when these functions are called.
    def push_upper(s):
        """
        Pushes s (a semiarc ID, as obtained through semiarc_map) to the stack
        representing the upper half-plane, and adds the current x value to the
        list of associated x values for the semiarc. See the initialization
        comment for `upper_cs` for more.
        """
        global states
        upper_cs.setdefault(s, []).append(x)
        upper.append(s)
        states += [
            (
                x,
                deepcopy(upper),
                deepcopy(lower),
                deepcopy(upper_cs),
                deepcopy(lower_cs),
                deepcopy(crossings),
            )
        ]

    def push_lower(s):
        global states
        lower_cs.setdefault(s, []).append(x)
        lower.append(s)
        # states += [(x, deepcopy(upper), deepcopy(lower), upper_cs, lower_cs, crossings)]
        states += [
            (
                x,
                deepcopy(upper),
                deepcopy(lower),
                deepcopy(upper_cs),
                deepcopy(lower_cs),
                deepcopy(crossings),
            )
        ]

    def pop_upper(expect=None):
        global states
        s = upper.pop()
        if expect is not None:
            assert expect == s
        upper_cs.setdefault(s, []).append(-x)
        # states += [(x, deepcopy(upper), deepcopy(lower), upper_cs, lower_cs, crossings)]
        states += [
            (
                x,
                deepcopy(upper),
                deepcopy(lower),
                deepcopy(upper_cs),
                deepcopy(lower_cs),
                deepcopy(crossings),
            )
        ]
        return s

    def pop_lower(expect=None):
        global states
        s = lower.pop()
        if expect is not None:
            assert expect == s
        lower_cs.setdefault(s, []).append(-x)
        # states += [(x, deepcopy(upper), deepcopy(lower), upper_cs, lower_cs, crossings)]
        states += [
            (
                x,
                deepcopy(upper),
                deepcopy(lower),
                deepcopy(upper_cs),
                deepcopy(lower_cs),
                deepcopy(crossings),
            )
        ]
        return s

    def push_or_pop_upper(s):
        """
        If s is the top element of upper, pop it off.

        Else, push it on.
        """
        if peek(upper, s):
            pop_upper(s)
        else:
            push_upper(s)

    def push_or_pop_lower(s):
        if peek(lower, s):
            pop_lower(s)
        else:
            push_lower(s)

    crossings = []

    # Each crossing appears in 4 semiarcs: in 2 as the start point, and in 2 as
    # the endpoint. Hence, len(semiarcs) // 2 + 1 gives `max_crossing_number` +
    # 1. So range(1, len(semiarcs) // 2 + 1) gives us an iterator for all the
    # crossings.
    #
    # i = crossing number
    for i in range(1, len(semiarcs) // 2 + 1):
        # l is the ID of the semiarc entering crossing i from the left.
        l = semiarc_map[i, Dir.LEFT]

        print(x, 1, lower, upper[::-1], l)

        # If this assertion fails, we would need to push the semiarc.
        # (like for the first arc special case), but how do we know
        # which side to generate it on?
        assert l in upper or l in lower or not (upper or lower)

        # If this assertion fails, how do we determine which direction
        # to shift in and how much to shift?
        assert (upper + lower).count(l) <= 1

        # If `l` has already been pushed to `upper` in the past, pop elements
        # from `upper` and push them to `lower` until we reach `l`, and then pop
        # `l`.
        if l in upper:
            while not peek(upper, l):
                push_or_pop_lower(pop_upper())
                x += 1

            pop_upper(l)

        # The same for lower
        elif l in lower:
            while not peek(lower, l):
                push_or_pop_upper(pop_lower())
                x += 1

            pop_lower(l)

        # In the event that we're treating the first crossing, both stacks are
        # empty. In this case we push to _both_, and clean up the issue at the
        # end.
        elif i == 1:
            push_lower(l)
            push_upper(l)

        # We're now done processing the left strand for crossing i, so move to
        # addressing the up/down strands. These are located one x value to the
        # right of the left strand's endpoint.
        x += 1

        # Get the semiarc IDs for the UP / DOWN strands
        u = semiarc_map[i, Dir.UP]
        d = semiarc_map[i, Dir.DOWN]
        print(x, 2, lower, upper[::-1], d, u)

        push_or_pop_upper(u)
        push_or_pop_lower(d)
        crossings.append(x)  # up / down share same x as the crossing

        x += 1

        # Process the rightward (exiting) strand for the crossing.
        r = semiarc_map[i, Dir.RIGHT]
        print(x, 3, lower, upper[::-1], r)

        # If this assertion fails, how will we know which side to pop
        # from? If we just choose an arbitrary side we won't fail to
        # be planar due to shifting, but we cannot guarantee
        # optimality.
        assert not (peek(upper, r) and peek(lower, r))

        if peek(upper, r):
            pop_upper(r)
        elif peek(lower, r):
            pop_lower(r)
        else:
            # If this assertion fails, then it is because we are on
            # the last crossing, but we failed to find the final
            # semiarc.
            assert (i + 1, Dir.LEFT) in semiarc_map
            l2 = semiarc_map[i + 1, Dir.LEFT]

            # If this assertion fails, how do we know which side to
            # push the semiarc on? If we just choose an arbitrary side
            # we won't fail to be planar due to shifting, but we
            # cannot guarantee optimality.
            assert l2 in upper + lower + [r]

            if l2 == r:
                # Arbitrary decision, either choice is optimal.
                push_upper(r)
            elif l2 in upper:
                push_lower(r)
            elif l2 in lower:
                push_upper(r)

        x += 1

    # The only things remaining in the stacks at this point should be
    #   (a) the copy of the first crossing's left strand (recall we pushed it
    #       once to each stack), and
    #   (b) strands that go around the right endpoint of the spine.
    # Hence, abs(len(upper) - len(lower)) <= 1.
    assert abs(len(upper) - len(lower)) <= 1
    while upper and lower:
        print(lower, upper[::-1])
        assert pop_upper() == pop_lower()
        x += 1

    print(lower, upper[::-1])

    # Cleanup the "leftover" return c
    if peek(upper, semiarc_map[1, Dir.LEFT]):
        upper.pop()
    elif peek(lower, semiarc_map[1, Dir.LEFT]):
        lower.pop()

    print(lower, upper[::-1])

    # Ensure that there's nothing else remaining to pop.
    assert upper == lower == []

    # Honestly this probably shouldn't be an inner function
    def deparenthesize(seq):
        """
        Pair up endpoints for departing / returning to the spine
        """
        result = []
        stack = []
        for x in seq:
            if x > 0:
                stack.append(x)
            else:
                y = stack.pop()
                result.append((y, abs(x)))
        return result

    upper_cs = {k: deparenthesize(v) for k, v in upper_cs.items()}
    lower_cs = {k: deparenthesize(v) for k, v in lower_cs.items()}

    return upper_cs, lower_cs, crossings, states


def crossing_tex(cnum, orient, x, xcurr):
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

    if cnum < 0:
        out_str += f"    \\draw ({x - 1}, 0) -- ({x - .1}, 0);\n"
        out_str += f"    \\draw[-latex] ({x + .1}, 0) -- ({x + .5}, 0);\n"
        out_str += f"    \\draw ({x + .5}, 0) -- ({x + 1}, 0);\n"
        if orient == 1:
            out_str += f"    \\draw[-latex] ({x}, .5) -- ({x}, -.5);\n"
        elif orient == -1:
            out_str += f"    \\draw[-latex] ({x}, -.5) -- ({x}, .5);\n"
        else:
            assert False
    else:
        out_str += f"    \\draw[-latex] ({x -1}, 0) -- ({x + 1}, 0);\n"
        if orient == 1:
            out_str += f"    \\draw ({x}, -.5) -- ({x}, -.1);\n"
            out_str += f"    \\draw[-latex] ({x}, .1) -- ({x}, .5);\n"
        elif orient == -1:
            out_str += f"    \\draw ({x}, .5) -- ({x}, .1);\n"
            out_str += f"    \\draw[-latex] ({x}, -.1) -- ({x}, -.5);\n"
        else:
            assert False

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
            f"    \\node[above left] ({abs(cnum)}lab) at ({abs(cnum)}circ) "
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
            f"    \\node[above left] ({abs(cnum)}lab) at ({abs(cnum)}circ) "
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


def arc_tex(x0, x1, xcurr, y, s, x0_through=False, x1_through=False):
    """

    """
    out_str = ""
    r = abs(x0 - x1) / 2
    xavg = (x1 + x0) / 2

    if y > 0:
        if xcurr < x1:
            nstyle = ", opacity=.2"
        else:
            nstyle = ""

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
        if xcurr < x0:
            out_str += f"    \\draw[opacity=.2] ({x0-.5}, 0) -- ({x1+.5}, 0);\n"
        elif x1 < xcurr:
            out_str += f"    \\draw ({x0 - .5}, 0) -- ({x1 + .5}, 0);\n"
        else:
            out_str += f"    \\draw ({x0 - .5}, 0) -- ({xcurr}, 0);\n"
            out_str += f"    \\draw[opacity=.2] ({xcurr}, 0) -- ({x1 + .5}, 0);\n"
        return out_str

    if x0 < x1:
        t0 = 180
        t1 = 0
    else:
        t0 = 0
        t1 = 180

    if xcurr < x0:
        out_str += f"    \\draw[opacity=.2] ({x0}, {y}) arc ({t0}:{t1}:{r});\n"
        # print(x0, x1, xcurr, y, s)
        # assert False

    elif xcurr < x1:
        # Find the angle to draw until such that our arc meets the
        # dotted line
        thalf = get_thalf(x0, x1, xcurr)
        out_str += f"    \\draw[opacity=.6] ({x0},  {y}) arc ({t0}:{thalf}:{r});\n"
        out_str += f"    \\draw[opacity=.2] ({x1},  {y}) arc ({t1}:{thalf}:{r});\n"

    else:
        out_str += f"    \\draw ({x0}, {y}) arc ({t0}:{t1}:{r});\n"

    if x0_through:
        if xcurr < x0:
            out_str += f"\\draw[opacity=.2] ({x0}, {y}) -- ({x0}, 0);\n"
        else:
            out_str += f"\\draw ({x0}, {y}) -- ({x0}, 0);\n"

    if x1_through:
        if xcurr < x1:
            out_str += f"\\draw[opacity=.2] ({x1}, {y}) -- ({x1}, 0);\n"
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


def plot_one_step(gc, orient, state, i, upper_cs_f, lower_cs_f, crossings_f, dirname):
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
\usepackage{tikz}
\begin{document}
  \begin{tikzpicture}
    """
    x, upper, lower, upper_cs, lower_cs, crossings = state
    print(state)
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


def build_frames(gkey):
    """
    gkey: an index into gauss_codes.gknot
    """
    cn, ind = gkey
    dirname = f"{cn}-{ind}"
    if not isdir(dirname):
        mkdir(dirname)

    [gc], orient = nelson_gc_to_sage_gc(gauss_codes.gknot[gkey])

    K = Knot([[gc], orient])
    crossings, semiarcs = knot_to_layout(K)

    # Awful. Keep a separate copy of the final state so that we can
    # draw it all in gray.
    upper_cs_f, lower_cs_f, crossings_f, states = route(semiarcs)
    for i, state in enumerate(states):
        plot_one_step(
            gc, orient, state, i, upper_cs_f, lower_cs_f, crossings_f, dirname
        )


def compile_single_file(fname):

    # run(["pdflatex", fname])
    run(
        [
            # "true",
            # "|",
            # ":",
            # "|",
            "pdflatex",
            "-halt-on-error",
            fname,
            "|",
            "grep",
            "'^!.*'",
            "-A200",
            "--color=always",
        ]
    )


def convert_single_file(fname):
    pdf = fname
    png = pdf[:-4] + ".png"
    run(
        ["convert", "-density", "300", pdf, "-flatten", "-quality", "96", png,]
    )


def compile_all():
    tex_files = [fname for fname in listdir() if ".tex" in fname]
    with Pool(5) as p:
        p.map(compile_single_file, tex_files)


def amalgamate(dirname, re_png=False, gif=False, mp4=True):

    if re_png:
        pdf_files = [fname for fname in listdir() if ".pdf" in fname]
        with Pool(5) as p:
            p.map(convert_single_file, pdf_files)
    if gif:
        gif_name = f"{dirname}.gif"
        run(
            [
                "convert",
                "-layers",
                "OptimizePlus",
                "-delay",
                "30",
                "*.png",
                "-loop",
                "0",
                gif_name,
            ]
        )
    if mp4:
        mov_name = f"{dirname}.mp4"
        run(
            [
                "ffmpeg",
                "-r",
                "2",
                "-i",
                "%03d.png",
                "-vcodec",
                "libx264",
                "-y",
                "-an",
                mov_name,
            ]
        )


if __name__ == "__main__":
    gkey = (7, 2)
    cn, ind = gkey
    dirname = f"{cn}-{ind}"
    build_frames(gkey)
    chdir(dirname)
    compile_all()
    amalgamate(dirname, re_png=True)
    chdir("..")
