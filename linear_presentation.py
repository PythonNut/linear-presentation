from sage.all import *
from enum import Enum
import itertools as it
import matplotlib.pyplot as plt
import matplotlib as mpl

from pprint import PrettyPrinter as PrettyPrinter


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


def sgn(x):
    "Return the sign of x"
    return 0 if x == 0 else -1 if x < 0 else 1


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
        elif sgn(b) == orient[abs(b) - 1]:
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


def virtual_route(semiarcs):
    """
    Perform the algorithm for virtual knots
    """
    upper, lower = [], []
    assert len(semiarcs) % 2 == 0

    # Create a dictionary for associating (crossing, dir) pairs with
    # an abstract ID representing the semiarc
    semiarc_map = {}

    # `a` = crossing a, `da` = exit  direction for `a`
    # `b` = crossing b, `db` = enter direction for `b`
    for i, (a, da, b, db) in enumerate(semiarcs):
        semiarc_map[a, da] = semiarc_map[b, db] = i

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

    # In the below, in general, pushing a semiarc to stack corresponds to
    # recording an x coordinate where it leaves the spine. Popping corresponds
    # to recording an x coordinate where it returns to the spine. Note, x will
    # be a variable in the outer scope when these functions are called.
    def peek(l, v):
        """
        Check if top element of the stack is v
        """
        return len(l) >= 1 and l[-1] == v

    def push_upper(s):
        """
        Pushes s (a semiarc ID, as obtained through semiarc_map) to the stack
        representing the upper half-plane, and adds the current x value to the
        list of associated x values for the semiarc. See the initialization
        comment for `upper_cs` for more.
        """
        upper_cs.setdefault(s, []).append(x)
        upper.append(s)

    def push_lower(s):
        lower_cs.setdefault(s, []).append(x)
        lower.append(s)

    def pop_upper(expect=None):
        s = upper.pop()
        if expect is not None:
            assert expect == s
        upper_cs.setdefault(s, []).append(-x)
        return s

    def pop_lower(expect=None):
        s = lower.pop()
        if expect is not None:
            assert expect == s
        lower_cs.setdefault(s, []).append(-x)
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

    def test_lower(s):
        if peek(lower, s):
            print("hi")

    def delpop_lower(s):
        """
        Remove the first occurence of s from lower, acting as if we
        were popping
        """
        # assert s in lower
        if peek(lower, s):
            pop_lower(s)
        else:
            for (i, el) in enumerate(lower):
                if el == s:
                    lower.pop(i)
                    lower_cs.setdefault(s, []).append(-x)
                    break

    def delpop_upper(s):
        """
        Remove the first occurence of s from lower, acting as if we
        were popping
        """
        # assert s in upper
        if peek(upper, s):
            pop_upper(s)
        else:
            for (i, el) in enumerate(upper):
                if el == s:
                    upper.pop(i)
                    upper_cs.setdefault(s, []).append(-x)
                    break

    crossings = []

    # Each crossing appears in 4 semiarcs: in 2 as the start point, and in 2 as
    # the endpoint. Hence, len(semiarcs) // 2 + 1 gives `max_crossing_number` +
    # 1. So range(1, len(semiarcs) // 2 + 1) gives us an iterator for all the
    # crossings.
    n = len(semiarcs) // 2 + 1
    #
    # i = crossing number
    #
    # We treat i = 1 as a special edge case.
    #
    # For the first crossing, guess whether we place the arc in
    # the upper or lower plane. We do this here instead of later
    # to avoid triggering the "if L in upper" cases, etc.
    i = 1
    L, U, D, R = (
        semiarc_map[(i, Dir.LEFT)],
        semiarc_map[(i, Dir.UP)],
        semiarc_map[(i, Dir.DOWN)],
        semiarc_map[(i, Dir.RIGHT)],
    )
    push_upper(L)
    push_lower(L)
    x += 1
    crossings.append(x)
    push_upper(U)
    push_lower(D)
    x += 1
    push_upper(R)
    x += 1

    for i in range(2, n):
        L, U, D, R = (
            semiarc_map[(i, Dir.LEFT)],
            semiarc_map[(i, Dir.UP)],
            semiarc_map[(i, Dir.DOWN)],
            semiarc_map[(i, Dir.RIGHT)],
        )

        print(i, "R  ", lower, upper[::-1], R)

        if L in upper:
            Lind = upper.index(L)
            # print("L in upper!")
            # print(L, U, D)
            # print(lower, upper[::-1])

            if U not in upper and U in lower:
                delpop_lower(U)

                # We want to insert U in the stack before L
                upper_cs.setdefault(U, []).append(x)
                upper.insert(Lind - 1, U)
                x += 1
                # print("U not in upper, U in lower")
                # print(lower, upper[::-1])
            elif U not in lower and U in upper and upper.index(U) > Lind:
                upper.remove(U)
                upper.insert(Lind - 1, U)
                # x += 1
                # print("elif case")
                # print(lower, upper[::-1])

            Lind = upper.index(L)
            if D not in lower and D in upper and upper.index(D) < Lind:
                # print("D not in lower blah blah")
                upper.remove(D)
                upper.insert(Lind, D)
                # x += 1
                # print(lower, upper[::-1])

            # elif D in lower and lower.index(D) > Lind:
            #     lower.remove(D)
            #     lower.insert(Lind, D)
            #     x += 1

            # Move stuff from upper to lower until we get to the
            # desired element
            while not peek(upper, L):
                push_or_pop_lower(pop_upper())
                x += 1

            print(L, D)

            # print(L, U)
            pop_upper(L)

        elif L in lower:
            Lind = lower.index(L)
            if D not in lower and D in upper:
                delpop_upper(D)
                lower_cs.setdefault(D, []).append(x)
                lower.insert(Lind - 1, D)
                x += 1

            elif D not in upper and D in lower and lower.index(D) > Lind:
                lower.remove(D)
                lower.insert(Lind - 1, D)

            Lind = lower.index(L)
            if U not in upper and U in lower and lower.index(U) < Lind:
                lower.remove(U)
                lower.insert(Lind, U)
                # x += 1

            # elif U in lower and lower.index(U) > Lind:
            #     lower.remove(U)
            #     lower.insert(Lind, U)
            #     x += 1

            print(L, D)

            while not peek(lower, L):
                push_or_pop_upper(pop_lower())
                x += 1

            # pop_upper(L)
            pop_lower(L)

        # If we're coming in from the left and have to unnest a bunch
        # of arcs on the top, then do it
        # if L in upper:
        #     # Move stuff from upper to lower until we get to the
        #     # desired element
        #     while not peek(upper, L):
        #         push_or_pop_lower(pop_upper())
        #         x += 1

        #     pop_upper(L)

        # elif L in lower:
        #     while not peek(lower, L):
        #         push_or_pop_upper(pop_lower())
        #         x += 1

        #     pop_lower(L)

        print(i, "U,D", lower, upper[::-1], U, D)

        x += 1

        # Handle the vertical strands
        if U in upper:
            delpop_upper(U)
        else:
            push_upper(U)

        if D in lower:
            delpop_lower(D)
        else:
            push_lower(D)

        crossings.append(x)

        x += 1

        print(i, "R  ", lower, upper[::-1], R)

        # This should only happen for the final strand
        if R in upper:
            delpop_upper(R)
        elif R in lower:
            delpop_lower(R)
        else:
            # if peek(upper, R):
            #     pop_upper(R)
            # elif peek(lower, R):
            #     pop_lower(R)
            # else:

            # ((i + 1) % n) + 1 gets us the next crossing, rolling
            # over to crossing 1 if we're on the final crossing.
            # print((i % n + 1))
            assert ((i % n) + 1, Dir.LEFT) in semiarc_map
            L2 = semiarc_map[((i % n) + 1), Dir.LEFT]

            # If this assertion fails, how do we know which side to
            # push the semiarc on? If we just choose an arbitrary side
            # we won't fail to be planar due to shifting, but we
            # cannot guarantee optimality.
            assert L2 in upper + lower + [R]

            if L2 == R:
                # Arbitrary decision, either choice is optimal.
                push_upper(R)
            elif L2 in upper:
                push_lower(R)
            elif L2 in lower:
                push_upper(R)

        x += 1

    # The only things remaining in the stacks at this point should be
    #   (a) the copy of the first crossing's left strand (recall we pushed it
    #       once to each stack), and
    #   (b) strands that go around the right endpoint of the spine.
    # Hence, abs(len(upper) - len(lower)) <= 1.
    # assert abs(len(upper) - len(lower)) <= 1

    print(lower, upper[::-1])
    for s in sorted(set(upper)):
        print(lower, upper[::-1])
        upper = [p for p in upper if p != s]
        lower = [p for p in lower if p != s]
        x += 1

    print(lower, upper[::-1])

    # Cleanup the "leftover" return c
    if peek(upper, semiarc_map[1, Dir.LEFT]):
        upper.pop()
    elif peek(lower, semiarc_map[1, Dir.LEFT]):
        lower.pop()

    print(lower, upper[::-1])

    # Ensure that there's nothing else remaining to pop.
    # assert upper == lower == []

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

    return upper_cs, lower_cs, crossings


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
        upper_cs.setdefault(s, []).append(x)
        upper.append(s)

    def push_lower(s):
        lower_cs.setdefault(s, []).append(x)
        lower.append(s)

    def pop_upper(expect=None):
        s = upper.pop()
        if expect is not None:
            assert expect == s
        upper_cs.setdefault(s, []).append(-x)
        return s

    def pop_lower(expect=None):
        s = lower.pop()
        if expect is not None:
            assert expect == s
        lower_cs.setdefault(s, []).append(-x)
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

    return upper_cs, lower_cs, crossings


if __name__ == "__main__":
    B = BraidGroup(4)
    K = Knot(B([1, 1, 1]))
    # K = Knot(B([1,1,1,2,-1,2,-3,2,-3]))

    # K = Knot([[3,1,2,4], [8,9,1,7], [5,6,7,3], [4,18,6,5],
    #           [17,19,8,18], [9,10,11,14], [10,12,13,11],
    #           [12,19,15,13], [20,16,14,15], [16,20,17,2]])

    K = Link(
        [
            [
                [
                    1,
                    -2,
                    -3,
                    -8,
                    -12,
                    13,
                    -14,
                    15,
                    -7,
                    -1,
                    2,
                    -4,
                    10,
                    11,
                    -13,
                    12,
                    -11,
                    -16,
                    4,
                    3,
                    -5,
                    6,
                    -9,
                    7,
                    -15,
                    14,
                    16,
                    -10,
                    8,
                    9,
                    -6,
                    5,
                ]
            ],
            [-1, -1, 1, 1, 1, 1, -1, 1, 1, -1, 1, -1, -1, -1, -1, -1],
        ]
    )

    # K = Knot([[[-1,2,-3,4,-5,6,-7,1,-4,8,-6,5,-8,3,-2,7]], [-1,-1,-1,-1,1,1,-1,-1]])
    # import gauss_codes
    # K = Knot(nelson_gc_to_sage_gc(gauss_codes.gknot[11,2]))
    # K = Knot(
    #     [
    #         [1, 5, 2, 4],
    #         [3, 8, 4, 9],
    #         [5, 11, 6, 10],
    #         [14, 7, 15, 8],
    #         [9, 2, 10, 3],
    #         [18, 12, 19, 11],
    #         [6, 13, 7, 14],
    #         [22, 15, 23, 16],
    #         [20, 18, 21, 17],
    #         [12, 20, 13, 19],
    #         [24, 21, 1, 22],
    #         [16, 23, 17, 24],
    #     ]
    # )
    # K = Knot(
    #     [
    #         [4, 2, 5, 1],
    #         [8, 4, 9, 3],
    #         [12, 9, 1, 10],
    #         [10, 5, 11, 6],
    #         [6, 11, 7, 12],
    #         [2, 8, 3, 7],
    #     ]
    # )

    # crossings, semiarcs = knot_to_layout(K)
    import gauss_codes

    for key in (
        # [(4, i) for i in range(82, 91)] + [(4, i) for i in range(96, 99)] + [(4, 107)]
        # [(4, 84 + i * 0.25) for i in range(4)]
        # (4, 26),
        gauss_codes.vknot.keys()
    ):
        gc = gauss_codes.vknot[(key)]
        # gc = gauss_codes.conn_sum(gauss_codes.gknot[10, 132], gauss_codes.gknot[8, 19])
        # gc = gauss_codes.conn_sum(gc, gauss_codes.gknot[6, 2], ind=10)
        # gc = gauss_codes.conn_sum(gc, gauss_codes.gknot[8, 13], ind=5)
        # K = Knot(nelson_gc_to_sage_gc(gc))
        K = Knot(nelson_gc_to_sage_gc(gc))

        print(key)
        crossings, semiarcs = knot_to_layout(K)
        # pp = PrettyPrinter(width=70, compact=True)
        # pp.pprint(crossings)
        # pp.pprint(semiarcs)

        try:
            upper_cs, lower_cs, crossings = virtual_route(semiarcs)
        except AssertionError as e:
            raise (e)
        plot(upper_cs, lower_cs, crossings, straight=0)

    # for name, n_gc in gauss_codes.gknot.items():
    #     print("="*10 + " " + str(name))
    #     K = Knot(nelson_gc_to_sage_gc(n_gc))
    #     crossings, semiarcs = knot_to_layout(K)
    #     plot(*route(semiarcs), straight=0)
