from sage.all import *
from enum import Enum
import itertools as it
from sage_numerical_backends_gurobi.gurobi_backend import GurobiBackend
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

def build_paths(MLP, n, bound, m):
    # Just create a set of n sequences of length m with one-hot values
    # in the interval [1 ... bound] along with a bitarray to track
    # parity values

    # shape n x m x bound
    paths = MLP.new_variable(integer=True, nonnegative=True)
    paths.set_max(bound)

    inds = MLP.new_variable(binary=True)

    ne_inds = MLP.new_variable(binary=True)

    for i in range(n):
        for j in range(m):
            MLP.add_constraint(inds[i, j] * bound >= paths[i, j])
            MLP.add_constraint(paths[i, j] >= inds[i, j])

        for j in range(m-1):
            MLP.add_constraint(inds[i, j] >= inds[i, j+1])

            hot = -2*bound*(2 - inds[i, j] - inds[i, j+1])
            MLP.add_constraint(paths[i, j] - paths[i, j+1] + 2*bound*ne_inds[i, j] >= 1 + hot)
            MLP.add_constraint(paths[i, j+1] - paths[i, j] + 2*bound*(1-ne_inds[i, j]) >= 1 + hot)


    # Track the parity of the demiarcs: A true value means the first
    # demiarc was in the top half.

    # shape n x (m-1)
    parity = MLP.new_variable(binary=True)

    return paths, parity, inds

def require_planarity(MLP, paths, parity, inds, n, bound, m):
    pln_inds = MLP.new_variable(binary=True)
    maxs = MLP.new_variable(integer=True, nonnegative=True)
    maxs.set_max(bound)
    max_inds = MLP.new_variable(binary=True)

    for i in range(n):
        for j in range(m-1):
            MLP.add_constraint(maxs[i, j] >= paths[i, j])
            MLP.add_constraint(maxs[i, j] >= paths[i, j+1])

            MLP.add_constraint(paths[i, j] + bound*max_inds[i, j] >= maxs[i, j])
            MLP.add_constraint(paths[i, j+1] + bound*(1-max_inds[i, j]) >= maxs[i, j])

    for i, k in it.combinations(range(n), 2):
        for j, l in it.product(range(m-1), range(m-1)):
            b = maxs[i, j]
            a = paths[i, j] + paths[i, j+1] - b
            d = maxs[k, l]
            c = paths[k, l] + paths[k, l+1] - d
            M = 2*bound + 1

            flipi = parity[i] if j %2 == 0 else (1 - parity[i])
            flipk = parity[k] if l %2 == 0 else (1 - parity[k])
            hot = M*(4 - inds[i, j] - inds[i, j+1] - inds[k, l] - inds[k, l+1])
            both_bot = M * (flipi + flipk)
            both_top = M * (2 - flipi - flipk)

            inside = pln_inds[i,j,k,l,0]
            branch = pln_inds[i,j,k,l,1]

            # a < b < c < d
            MLP.add_constraint(b <= c-1 + M*(1-inside) + M*branch + hot + both_top)
            MLP.add_constraint(b <= c-1 + M*(1-inside) + M*branch + hot + both_bot)

            # c < d < a < b
            MLP.add_constraint(d <= a-1 + M*(1-inside) + M*(1-branch) + hot + both_top)
            MLP.add_constraint(d <= a-1 + M*(1-inside) + M*(1-branch) + hot + both_bot)

            # a < c < d < b
            MLP.add_constraint(a <= c-1 + M*inside + M*branch + hot + both_top)
            MLP.add_constraint(d <= b-1 + M*inside + M*branch + hot + both_top)
            MLP.add_constraint(a <= c-1 + M*inside + M*branch + hot + both_bot)
            MLP.add_constraint(d <= b-1 + M*inside + M*branch + hot + both_bot)

            # c < a < b < d
            MLP.add_constraint(c <= a-1 + M*inside + M*(1-branch) + hot + both_top)
            MLP.add_constraint(b <= d-1 + M*inside + M*(1-branch) + hot + both_top)
            MLP.add_constraint(c <= a-1 + M*inside + M*(1-branch) + hot + both_bot)
            MLP.add_constraint(b <= d-1 + M*inside + M*(1-branch) + hot + both_bot)

    return pln_inds

def require_connections(MLP, semiarcs, paths, parity, inds, step_size, shift, m, bound):
    for i, (a, da, b, db) in enumerate(semiarcs):
        MLP.add_constraint(paths[i, m-1] == 0)

        # Taking care of the "a" conditions is easy
        a_pos = step_size * a + shift
        if da == Dir.UP:
            MLP.add_constraint(parity[i] == 1)
        elif da == Dir.DOWN:
            MLP.add_constraint(parity[i] == 0)
        elif da == Dir.RIGHT:
            a_pos += 1
            if db == Dir.LEFT:
                MLP.add_constraint(parity[i] == 1)

        MLP.add_constraint(paths[i, 0] == a_pos)

        # Now, for the "b" conditions
        b_pos = step_size * b + shift
        if db == Dir.LEFT: b_pos -= 1

        for j in range(0, m-1):
            this_slice = inds[i, j]
            next_slice = inds[i, j+1]
            MLP.add_constraint(2*bound * ((1 - this_slice) + next_slice) >= paths[i, j] - b_pos)
            MLP.add_constraint(2*bound * ((1 - this_slice) + next_slice) >= b_pos - paths[i, j])

            flip = parity[i] if j %2 == 0 else (1 - parity[i])
            if db == Dir.UP:
                MLP.add_constraint(1-flip >= this_slice - next_slice)
            if db == Dir.DOWN:
                MLP.add_constraint(flip >= this_slice - next_slice)

def unpack_paths(paths, n, bound, m):
    result = []
    for i in range(n):
        path = [paths[i, j] for j in range(m) if paths[i, j] > 0]
        result.append(path)
    return result


def plot(paths, parity, crossings=[]):
    def add_arc(a, b, y=0, down=False):
        h = abs(a-b)
        c = (a + b)/2
        if down:
            plt.gca().add_patch(mpl.patches.Arc((c, y), h, h, theta1=180))
        else:
            plt.gca().add_patch(mpl.patches.Arc((c, y), h, h, theta2=180))

    for i, path in enumerate(paths):
        par = parity[i]
        for j, (a, b) in enumerate(zip(path, path[1:])):
            # print(j, (par + j - 1) % 2 == 0)
            add_arc(a, b, 0, (par + j) % 2 == 0)

    for x in crossings:
        plt.plot([x-1, x+1], [0, 0], 'k')

    plt.gcf().set_size_inches(10.5, 18.5)
    plt.axis('scaled')
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

if __name__ == "__main__":
    B = BraidGroup(4)
    K = Knot(B([1,1,1]))
    # K = Knot(B([1,1,1,2,-1,2,-3,2,-3]))

    # K = Knot([[3,1,2,4], [8,9,1,7], [5,6,7,3], [4,18,6,5],
    #           [17,19,8,18], [9,10,11,14], [10,12,13,11],
    #           [12,19,15,13], [20,16,14,15], [16,20,17,2]])

    # K = Link([[[1,-2,-3,-8,-12,13,-14,15,-7,-1,2,-4,10,11,-13,12,
    #             -11,-16,4,3,-5,6,-9,7,-15,14,16,-10,8,9,-6,5]],
    #           [-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1]])

    # K = Knot([[[-1,2,-3,4,-5,6,-7,1,-4,8,-6,5,-8,3,-2,7]], [-1,-1,-1,-1,1,1,-1,-1]])
    K = Knot([[1,5,2,4],[3,8,4,9],[5,11,6,10],[14,7,15,8],[9,2,10,3],[18,12,19,11],[6,13,7,14],[22,15,23,16],[20,18,21,17],[12,20,13,19],[24,21,1,22],[16,23,17,24]])
    # K = Knot([[4,2,5,1],[8,4,9,3],[12,9,1,10],[10,5,11,6],[6,11,7,12],[2,8,3,7]])

    crossings, semiarcs = knot_to_layout(K)

    n = len(crossings)
    m = 8

    MLP = MixedIntegerLinearProgram(maximization=False, solver=GurobiBackend)
    paths, parity, inds = build_paths(MLP, 2*n, n**2, m)
    pln_inds = require_planarity(MLP, paths, parity, inds, 2*n, n**2 + 2, m)
    require_connections(MLP, semiarcs, paths, parity, inds, n, -n+2, m, n**2 + 2)
    MLP.set_objective(MLP.sum(inds.values()))

    MLP.solve(log=1)

    got_paths = unpack_paths(MLP.get_values(paths), len(semiarcs), n ** 2, m)
    got_parity = MLP.get_values(parity)
