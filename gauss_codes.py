def gc_to_tex_str(gc):
    """
    Convert the sage gauss code to a string for formatting in LaTeX
    """
    out_str = ""

    for c in gc[0][0]:
        csign = gc[1][abs(c) - 1]

        if c < 0:
            uo = "u"
        else:
            uo = "o"

        if csign < 0:
            sgn = "-"
        else:
            sgn = "+"

        out_str += f"{abs(c)}_{uo}^{sgn} "

    return out_str


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


def OU_to_nums(gc):
    """
    convert an input gauss code of the form
    (str):  "O1-O2-U1-U2-..."
    to the nelson gauss code format
    (list): [1.5, 2.5, -1.5, -2.5, ...]
    """
    new_gc = []
    segs = [gc[3 * i : 3 * i + 3] for i in range(len(mystr) // 3)]
    for seg in segs:
        n = int(seg[1])
        if seg[2] == "-":
            n += 0.5
        if seg[0] == "U":
            n *= -1
        new_gc += n
    return new_gc


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


def conn_sum(gcode0, gcode1, ind=0):
    c0nums = {get_cnum(c) for c in gcode0}
    cmax = max(c0nums)
    print(cmax)
    gcode1 = [c + cmax if c > 0 else c - cmax for c in gcode1]
    return gcode0[:ind] + gcode1 + gcode0[ind:]


def multi_sum_from_gknot(idxs, inds=None):
    inds = inds or [0] * (len(idxs) - 1)
    knot = gknot[idxs[0]]
    for i, idx in zip(inds, idxs[1:]):
        knot = conn_sum(knot, gknot[idx], ind=i)
    return knot


def dt2gauss(dt):
    x, y = [], []
    for i, dti in enumerate(dt):
        x.extend((2 * i + 1, dti))

    for i in range(len(dt)):
        y.extend("ou" if x[2 * i + 1] > 0 else "uo")

    gauss = []
    for i in range(1, len(x) + 1):
        j = x.index(i)
        gauss.append((j // 2 + 1) * (-1 if y[j] == "u" else 1))

    return gauss


def gauss2signed_gauss(gcode):
    import networkx as nx

    n = len(gcode) // 2

    G = nx.Graph()
    crossings = []
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
        j = (i + 1) % len(gcode)
        cross1, o1 = abs(gcode[i]), gcode[i] > 0
        cross2, o2 = abs(gcode[j]), gcode[j] > 0

        t1, l1, c1, r1, b1 = crossings[cross1 - 1]
        t2, l2, c2, r2, b2 = crossings[cross2 - 1]

        src = b1 if o1 else r1
        dst = t2 if o2 else l2

        G.add_edge(c1, src)
        G.add_edge(src, dst)
        G.add_edge(dst, c2)

    embed = nx.check_planarity(G)[1]
    orient = []

    for t, l, c, r, b in crossings:
        o = list(embed.neighbors_cw_order(c))
        if o[(o.index(t) + 1) % 4] == r:
            orient.append(+1)
        else:
            orient.append(-1)

    return [gcode], orient


def get_cknots(max_cnum=10):
    """
    Get all classical prime knots up to 10 crossings
    """
    from data.classical_gauss_codes import gknot

    return dict(
        filter(lambda key_val_pair: key_val_pair[0][0] <= max_cnum, gknot.items())
    )


def get_vknots(max_cnum=5):
    from data.virtual_gauss_codes import vknot

    return dict(
        filter(lambda key_val_pair: key_val_pair[0][0] <= max_cnum, vknot.items())
    )
