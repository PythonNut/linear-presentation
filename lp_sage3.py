from sage.all import *
from enum import Enum
import itertools as it
import matplotlib.pyplot as plt
import matplotlib as mpl

class Dir(Enum):
    UP = 0
    RIGHT = 1
    DOWN = 2
    LEFT = 3

    def flip(self):
        return Dir((self.value + 2)%4)

def fix_gc_order(gc, orient):
    new_gc = []
    mapping = {}

    for elt in gc:
        c, s = abs(elt), sign(elt)
        if c not in mapping:
            mapping[c] = len(mapping) + 1
        new_c = mapping[c]
        new_gc.append(new_c * s)

    rev_mapping = {v:k for k, v in mapping.items()}
    assert len(rev_mapping) == len(mapping)
    new_orient = [orient[rev_mapping[c]-1] for c in range(1, len(mapping)+1)]

    return new_gc, new_orient

def linear_layout(gc, orient):
    n = max(gc)
    crossings = {i:{} for i in range(1,n+1)}
    semiarcs = []
    # gc = gc + [gc[0]]
    seen = set([1])
    exit_dir = Dir.RIGHT
    for a, b in zip(gc, gc[1:]):
        if abs(b) not in seen:
            entry_dir = Dir.LEFT
        elif sign(b) == orient[abs(b)-1]:
            entry_dir = Dir.UP
        else:
            entry_dir = Dir.DOWN

        crossings[abs(a)][exit_dir] = (abs(b), entry_dir)
        crossings[abs(b)][entry_dir] = (abs(a), exit_dir)
        semiarcs.append((abs(a), exit_dir, abs(b), entry_dir))
        seen.add(abs(a))
        exit_dir = entry_dir.flip()

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

def plot(top_demiarcs, bot_demiarcs, crossings=[], straight=0):

    def demiarc_height(demiarcs):
        all_arcs = []
        for arcs in demiarcs.values():
            all_arcs.extend(arcs)

        dominate = {}
        for a, b in all_arcs:
            for c, d in all_arcs:
                if a < c < d < b:
                    dominate.setdefault((a, b), []).append((c, d))

        height = {}
        for a, b in sorted(all_arcs, key=lambda t:len(dominate.get(t, []))):
            height[a, b] = 0
            for dominated in dominate.get((a, b), []):
                height[a, b] = max(height[a, b], height[dominated])

            if b - a > 1:
                height[a, b] += 1

        return height

    top_heights = demiarc_height(top_demiarcs)
    bot_heights = demiarc_height(bot_demiarcs)

    def add_arc(a, b, y=0, down=False):
        h = abs(a-b)
        if straight >= 3:
            if down:
                h = bot_heights[a, b]
            else:
                h = top_heights[a, b]

        if straight >= 1 and abs(a-b) <= 1:
            plt.plot([a, b], [0, 0], 'k')

        elif straight >= 2:
            if down:
                plt.plot([a, a, b, b], [0, -h, -h, 0], 'k')
            else:
                plt.plot([a, a, b, b], [0, h, h, 0], 'k')

        else:
            c = (a + b)/2
            if down:
                plt.gca().add_patch(mpl.patches.Arc((c, y), h, h, theta1=180))
            else:
                plt.gca().add_patch(mpl.patches.Arc((c, y), h, h, theta2=180))

    for arcs in top_demiarcs.values():
        for (a, b) in arcs:
            add_arc(a, b, 0, False)

    for arcs in bot_demiarcs.values():
        for (a, b) in arcs:
            add_arc(a, b, 0, True)

    # for i, path in enumerate(paths):
    #     par = parity[i]
    #     for j, (a, b) in enumerate(zip(path, path[1:])):
    #         add_arc(a, b, 0, (par + j) % 2 == 0)

    for x in crossings:
        plt.plot([x-1, x+1], [0, 0], 'k')

    plt.gcf().set_size_inches(10.5, 18.5)
    plt.axis('off')
    plt.axes().set_aspect('equal', 'datalim')
    plt.show()

def nelson_gc_to_sage_gc(gc):
    new_gc = [int(x) for x in gc]
    assert len(new_gc)%2 == 0
    n = len(new_gc)//2
    new_orient = [1 if i in gc else -1 for i in range(1, n+1)]
    return [[new_gc], new_orient]

def knot_to_layout(K):
    [gc], orient = K.oriented_gauss_code()
    gc, orient = fix_gc_order(gc, orient)
    crossings, semiarcs = linear_layout(gc, orient)
    return crossings, semiarcs

def plain_semiarcs(semiarcs):
    return [(a, b.value, c, d.value) for a,b,c,d in semiarcs]

def route(semiarcs):
    top, bot = [], []
    assert len(semiarcs)%2 == 0

    semiarc_map = {}
    for i, (a, da, b, db) in enumerate(semiarcs):
        semiarc_map[a, da] = semiarc_map[b, db] = i

    def peek(l, v):
        return len(l) >= 1 and l[-1] == v

    top_demiarcs = {}
    bot_demiarcs = {}

    x = 1
    def push_top(v):
        top_demiarcs.setdefault(v, []).append(x)
        top.append(v)

    def push_bot(v):
        bot_demiarcs.setdefault(v, []).append(x)
        bot.append(v)

    def pop_top():
        v = top.pop()
        top_demiarcs.setdefault(v, []).append(-x)
        return v

    def pop_bot():
        v = bot.pop()
        bot_demiarcs.setdefault(v, []).append(-x)
        return v

    def push_or_pop_top(v):
        if peek(top, v):
            pop_top()
        else:
            push_top(v)

    def push_or_pop_bot(v):
        if peek(bot, v):
            pop_bot()
        else:
            push_bot(v)

    crossings = []

    for i in range(1, len(semiarcs)//2 + 1):
        l = semiarc_map[i, Dir.LEFT]
        print(x, 1, bot, top[::-1], l)
        # If this assertion fails, we would need to push the semiarc.
        # (like for the first arc special case), but how do we know
        # which side to generate it on?
        assert l in top or l in bot or not (top or bot)

        # If this assertion fails, how do we determine which direction
        # to shift in?
        assert not (l in top and l in bot)
        if l in top:
            while not peek(top, l):
                push_or_pop_bot(pop_top())
                x += 1

            pop_top()

        elif l in bot:
            while not peek(bot, l):
                push_or_pop_top(pop_bot())
                x += 1

            pop_bot()

        elif not (top or bot):
            # Special case: we don't know anything so we hedge our
            # bets and cleanup at the end
            push_bot(l)
            push_top(l)

        x += 1

        u = semiarc_map[i, Dir.UP]
        d = semiarc_map[i, Dir.DOWN]
        print(x, 2, bot, top[::-1], d, u)
        push_or_pop_top(u)
        push_or_pop_bot(d)
        crossings.append(x)

        x += 1

        r = semiarc_map[i, Dir.RIGHT]
        print(x, 3, bot, top[::-1], r)

        # If this assertion fails, how will we know which side to pop
        # from? If we just choose an arbitrary side we won't fail to
        # be planar due to shifting, but we cannot guarantee
        # optimality.
        assert not (peek(top, r) and peek(bot, r))

        if peek(top, r):
            pop_top()
        elif peek(bot, r):
            pop_bot()
        else:
            # If this assertion fails, then it is because we are on
            # the last crossing, but we failed to find the final
            # semiarc.
            assert (i+1, Dir.LEFT) in semiarc_map
            l2 = semiarc_map[i+1, Dir.LEFT]

            # If this assertion fails, how do we know which side to
            # push the semiarc on? If we just choose an arbitrary side
            # we won't fail to be planar due to shifting, but we
            # cannot guarantee optimality.
            assert l2 in top + bot + [r]

            if l2 == r:
                # Arbitrary decision, either choice is optimal.
                push_top(r)
            elif l2 in top:
                push_bot(r)
            elif l2 in bot:
                push_top(r)

        x += 1

    assert -1 <= len(top) - len(bot) <= 1
    while top and bot:
        print(bot, top[::-1])
        assert pop_top() == pop_bot()
        x += 1

    print(bot, top[::-1])

    # Cleanup the "leftover" return demiarc
    if peek(top, semiarc_map[1, Dir.LEFT]):
        top.pop()
    elif peek(bot, semiarc_map[1, Dir.LEFT]):
        bot.pop()

    print(bot, top[::-1])
    assert top == bot == []

    # Honestly this probably shouldn't be an inner function
    def deparenthesize(seq):
        result = []
        stack = []
        for x in seq:
            if x > 0:
                stack.append(x)
            else:
                y = stack.pop()
                result.append((y, abs(x)))

        return result

    top_demiarcs = {k:deparenthesize(v) for k, v in top_demiarcs.items()}
    bot_demiarcs = {k:deparenthesize(v) for k, v in bot_demiarcs.items()}

    return top_demiarcs, bot_demiarcs, crossings


if __name__ == "__main__":
    B = BraidGroup(4)
    K = Knot(B([1,1,1]))
    # K = Knot(B([1,1,1,2,-1,2,-3,2,-3]))

    # K = Knot([[3,1,2,4], [8,9,1,7], [5,6,7,3], [4,18,6,5],
    #           [17,19,8,18], [9,10,11,14], [10,12,13,11],
    #           [12,19,15,13], [20,16,14,15], [16,20,17,2]])

    K = Link([[[1,-2,-3,-8,-12,13,-14,15,-7,-1,2,-4,10,11,-13,12,
                -11,-16,4,3,-5,6,-9,7,-15,14,16,-10,8,9,-6,5]],
              [-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1]])

    # K = Knot([[[-1,2,-3,4,-5,6,-7,1,-4,8,-6,5,-8,3,-2,7]], [-1,-1,-1,-1,1,1,-1,-1]])
    # import gauss_codes
    # K = Knot(nelson_gc_to_sage_gc(gauss_codes.gknot[11,2]))
    # K = Knot([[1,5,2,4],[3,8,4,9],[5,11,6,10],[14,7,15,8],[9,2,10,3],[18,12,19,11],[6,13,7,14],[22,15,23,16],[20,18,21,17],[12,20,13,19],[24,21,1,22],[16,23,17,24]])
    # K = Knot([[4,2,5,1],[8,4,9,3],[12,9,1,10],[10,5,11,6],[6,11,7,12],[2,8,3,7]])

    # crossings, semiarcs = knot_to_layout(K)
    import gauss_codes
    K = Knot(nelson_gc_to_sage_gc(gauss_codes.gknot[10, 132]))
    crossings, semiarcs = knot_to_layout(K)
    plot(*route(semiarcs), straight=0)

    # for name, n_gc in gauss_codes.gknot.items():
    #     print("="*10 + " " + str(name))
    #     K = Knot(nelson_gc_to_sage_gc(n_gc))
    #     crossings, semiarcs = knot_to_layout(K)
    #     plot(*route(semiarcs), straight=0)
