from copy import deepcopy
import math
import itertools as it

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl

from gauss_codes import gknot, conn_sum
import linear_presentation as lp

def rotate_right(l, n):
    return l[-n:] + l[:-n]

def irotate_right(l, n):
    return it.chain(l[-n:], l[:-n])

def build_embedding_graph(gcode):
    n = len(gcode) // 2

    G = nx.Graph()

    crossings = []
    knot_order = []
    # First, build crossing gadgets
    for i in range(n):
        nodes = range(5*i, 5*i+5)
        G.add_nodes_from(nodes)

        t, l, c, r, b = nodes
        G.add_edge(l, t)
        G.add_edge(r, t)
        G.add_edge(b, l)
        G.add_edge(b, r)

        crossings.append(nodes)

    for i in range(len(gcode)):
        # Get information from the gauss code
        cross1, o1, p1 = lp.decompose_c(gcode[i])
        cross2, o2, p2 = lp.decompose_c(gcode[(i + 1) % len(gcode)])
        print(cross1, cross2, o1, o2, p1, p2)

        t1, l1, c1, r1, b1 = crossings[cross1-1]
        t2, l2, c2, r2, b2 = crossings[cross2-1]

        if o1 == 1: src = t1
        elif p1 == 1: src = l1
        else: src = r1

        if o2 == 1: dst = b2
        elif p2 == 1: dst = r2
        else: dst = l2

        G.add_edge(c1, src)
        G.add_edge(src, dst)
        G.add_edge(dst, c2)
        knot_order.append([c1, src, dst])

    assert nx.check_planarity(G)[0]

    print(crossings)

    return G, knot_order

def clean_graph(G):
    G = deepcopy(G)
    n = len(G.nodes)//5
    for i in range(n):
        t, l, c, r, b = range(5*i, 5*i+5)
        G.remove_edge(l, t)
        G.remove_edge(r, t)
        G.remove_edge(b, l)
        G.remove_edge(b, r)
    return G

def planar_plot(G, pos=None):
    if pos is None:
        embed = nx.check_planarity(G)[1]
        pos = nx.combinatorial_embedding_to_pos(embed)
        return planar_plot(G, pos)

    print(pos)
    # G = clean_graph(G)
    # trace_graph(G)
    nx.draw(G, pos, node_color = [['red', 'blue', 'green', 'yellow', 'pink'][i%5] for i in G.nodes], with_labels=True)
    plt.show()

def walk_face(embed, src, dst):
    face = [src, dst]

    while True:
        last, new = embed.next_face_half_edge(*face[-2:])
        if [last, new] == face[:2]: break
        face.append(new)

    return face[:-1]

def get_outer_face(G, embed=None, pos=None):
    if pos is None or embed is None:
        embed = nx.check_planarity(G)[1]
        pos = nx.combinatorial_embedding_to_pos(embed)
        return get_outer_face(G, embed, pos)

    leftmost = min(G.nodes, key=lambda v: pos[v][0])
    topmost_neighbor = max(G[leftmost], key=lambda v: pos[v][1])
    return walk_face(embed, topmost_neighbor, leftmost)

def get_all_faces(G, embed):
    faces = []
    half_edges = set(it.chain(G.edges, map(lambda e:e[::-1], G.edges)))
    while half_edges:
        e = next(iter(half_edges))
        face = walk_face(embed, *e)
        faces.append(face)
        for fe in zip(face, irotate_right(face, -1)):
            half_edges.remove(fe)

    return faces

def merge_paths(h1, h2):
    g1 = [h//5 for h in h1]
    g2 = [h//5 for h in h2]
    g_common = next(iter(set(g1) & set(g2)))
    i1 = [(i, h - g_common*5) for i, h in enumerate(h1) if h//5 == g_common]
    i2 = [(i, h - g_common*5) for i, h in enumerate(h2) if h//5 == g_common]

    # Enforce convenient constraints that apply WLOG
    if len(i1) == 2:
        print("swapping paths")
        return merge_paths(h2, h1)

    if i1[0][0] != 0:
        print("rotate path 1")
        return merge_paths(rotate_right(h1, -i1[0][0]), h2)

    if i2[0][0] != 0:
        print("rotate path 2")
        return merge_paths(h1, rotate_right(h2, -i2[0][0]))

    edge_order = [0, 1, 4, 3]
    pre_flip = i1[2][1]
    post_flip = edge_order[(edge_order.index(i1[2][1]) + 2)%len(edge_order)]
    # print(i1, i2)

    if i2[0][1] == post_flip:
        print("inverting path 2")
        insert_seq = rotate_right(h2[::-1], 1)
    else:
        insert_seq = rotate_right(h2, -1)

    print(pre_flip, post_flip, insert_seq)
    assert insert_seq[0] == post_flip + g_common*5
    assert insert_seq[-1] not in (
        pre_flip + g_common*5,
        post_flip + g_common*5,
        i1[0][1] + g_common*5
    )

    return [*h1[:2], *insert_seq, *h1[2:]]

def extract_cycle(G, start=None):
    n = G.number_of_nodes()
    G_tmp = deepcopy(G)
    v = start or next(iter(G_tmp.nodes))
    assert v%5 != 2
    end = set(G_tmp[v])
    path = []
    while True:
        if not G_tmp[v]:
            break

        print(f"At {v}")
        # planar_plot(G_tmp)
        path.append(v)
        G_tmp.remove_node(v)
        gadget, offset = divmod(v, 5)
        g_i = gadget * 5
        edge_order = [0, 1, 4, 3]

        r = g_i + edge_order[(offset + 1)%4]
        l = g_i + edge_order[(offset - 1)%4]

        out = r if r in G_tmp else l
        print(f"Out {out}")
        c = g_i + 2
        if c in G_tmp:
            print(f'Through center {c}')
            path.append(c)
            G_tmp.remove_node(c)
            # planar_plot(G_tmp)

        saved_neighbors = G_tmp[out]

        path.append(out)
        G_tmp.remove_node(out)

        # Find the edge that actually leaves this gadget
        for after in saved_neighbors:
            if after < g_i or after >= g_i + 5:
                v = after
                break
        else:
            # TODO: Do I still need this clause?
            break

    return path, G_tmp

def check_hamiltonian(G, h):
    assert set(h) == set(G.nodes)

    for v1, v2 in zip(h, h[1:]):
        assert (v1, v2) in G.edges

    assert (h[-1], h[0]) in G.edges

def hamiltonian_cycle(G):
    G_tmp = deepcopy(G)
    components = []

    # TODO: argument for why decomposing a graph into cycles yields
    # cycles that may always be fused by swapping gadget parities
    while True:
        path, G_tmp = extract_cycle(G_tmp)
        components.append(path)
        if G_tmp.number_of_nodes() == 0:
            break

    # TODO: can we do this in linear time?
    while len(components) > 1:
        for (i1, c1), (i2, c2) in it.combinations(enumerate(components), 2):
            if set(h//5 for h in c1) & set(h//5 for h in c2):
                # careful, pop the last one first
                components.pop(i2)
                components.pop(i1)
                components.append(merge_paths(c1, c2))
                break
        else:
            assert False

    [ham] = components
    check_hamiltonian(G, ham)
    return ham

def right_outwards_hamiltonian(G, ham):
    # Note that this is only because the outer face is right outwards
    # oriented, due to the embedding being ccw ordered.
    outer_face = get_outer_face(G)
    ham_dedges = set(zip(ham, rotate_right(ham, -1)))
    outer_face_dedges = set(zip(outer_face, rotate_right(outer_face, -1)))

    # TODO: topological argument why at least one edge of the
    # hamiltonian cycle must fall on the outside face
    if ham_dedges & outer_face_dedges:
        return ham

    return ham[::-1]

def pos_angle(theta):
    if theta < 0:
        return theta + 2*math.pi
    return theta

def check_closed_under_flip(edges):
    for u, v in edges:
        assert (v, u) in edges

def deduplicate_edges(edges):
    result = set()
    for u, v in edges:
        result.add((min(u, v), max(u, v)))
    return result

def classify_edges(G, pos, ham):
    ham =  right_outwards_hamiltonian(G, ham)
    cycle = []
    inside = []
    outside = []
    for h1, h2, h3 in set(zip(ham, rotate_right(ham, -1), rotate_right(ham, -2))):
        (h1x, h1y), (h2x, h2y), (h3x, h3y) = pos[h1], pos[h2], pos[h3]
        prev_slope = math.atan2(h1y - h2y, h1x - h2x)
        next_slope = math.atan2(h3y - h2y, h3x - h2x)
        # h1 <---- h2 ----> h3
        for n in G[h2]:
            if n in (h1, h3):
                cycle.append((h2, n))
                continue

            x, y = pos[n]
            slope = math.atan2(y - h2y, x - h2x)
            r1 = pos_angle(slope - next_slope)
            r2 = pos_angle(prev_slope - next_slope)
            if r1 > r2:
                inside.append((h2, n))
            else:
                outside.append((h2, n))

    check_closed_under_flip(cycle)
    check_closed_under_flip(inside)
    check_closed_under_flip(outside)

    return (
        deduplicate_edges(cycle),
        deduplicate_edges(inside),
        deduplicate_edges(outside)
    )

def zerotone_book_embedding_arcs(G):
    embed = nx.check_planarity(G)[1]
    pos = nx.combinatorial_embedding_to_pos(embed)
    ham = hamiltonian_cycle(G)

    arcs = []

    cycle, inside, outside = classify_edges(G, pos, ham)
    n = G.number_of_nodes()
    x_map = {x:i for i, x in enumerate(ham)}

    arcs.append((n-1, 0))
    for i in range(n-1):
        arcs.append((i+1, i))

    for a, b in inside:
        x1, x2 = x_map[a], x_map[b]
        arcs.append((min(x1, x2), max(x1, x2)))

    for a, b in outside:
        x1, x2 = x_map[a], x_map[b]
        arcs.append((max(x1, x2), min(x1, x2)))

    return arcs

def draw_zerotone_book_embedding(arcs, ham):
    # embed = nx.check_planarity(G)[1]
    # pos = nx.combinatorial_embedding_to_pos(embed)
    # ham = hamiltonian_cycle(G)

    # cycle, inside, outside = classify_edges(G, pos, ham)
    # arcs = zerotone_book_embedding_arcs(G)

    for i, v in enumerate(ham):
        plt.plot([i], [0], 'ro-')
        plt.annotate(str(v), (i, 0))

    def add_arc(a, b, down=False):
        h = abs(a-b)
        c = (a + b)/2
        if down:
            plt.gca().add_patch(mpl.patches.Arc((c, 0), h, h, theta1=180))
        else:
            plt.gca().add_patch(mpl.patches.Arc((c, 0), h, h, theta2=180))

    for a, b in arcs:
        add_arc(a, b, b<a)

    # x_map = {x:i for i, x in enumerate(ham)}

    # n = G.number_of_nodes()
    # add_arc(0, n-1, True)
    # for i in range(n-1):
    #     add_arc(i, i+1, True)

    # for a, b in inside:
    #     add_arc(x_map[a], x_map[b])

    # for a, b in outside:
    #     add_arc(x_map[a], x_map[b], True)

    plt.gcf().set_size_inches(15, 15)
    plt.axis('scaled')
    plt.show()

def ham_plot(G, cycle, inside, outside):
    embed = nx.check_planarity(G)[1]
    pos = nx.combinatorial_embedding_to_pos(embed)
    edge_colors = []
    for edge in G.edges:
        if edge in cycle:
            edge_colors.append('black')
        elif edge in inside:
            edge_colors.append('blue')
        elif edge in outside:
            edge_colors.append('red')
        else:
            assert False

    print(edge_colors)
    nx.draw(G, pos, node_color = [['red', 'blue', 'green', 'yellow', 'pink'][i%5] for i in G.nodes], with_labels=True, edge_color=edge_colors)
    plt.show()

def linear_order_draw_plot(ham, upper, inter, lower, target_ordering, hide_upper_spine=False):
    n = len(ham)
    if not hide_upper_spine:
        for i, v in enumerate(ham):
            plt.plot([i], [0], 'ro-')
            plt.annotate(str(v), (i, 0))

    for i, v in enumerate(target_ordering):
        plt.plot([i], [-2*n], 'ro-')
        plt.annotate(str(v), (i, -2*n))

    def add_arc(a, b, y=0, down=False):
        h = abs(a-b)
        c = (a + b)/2
        if down:
            plt.gca().add_patch(mpl.patches.Arc((c, y), h, h, theta1=180))
        else:
            plt.gca().add_patch(mpl.patches.Arc((c, y), h, h, theta2=180))

    for a, b in upper:
        add_arc(a, b, 0, a > b)

    for a, b in lower:
        add_arc(a, b, -2*n, a > b)

    for u, l in inter:
        plt.plot([u, l], [0, -2*n], 'k-')

    plt.gcf().set_size_inches(10.5, 18.5)
    plt.axis('scaled')
    plt.show()

def linear_order_draw(G, target_ordering, plot=False):
    n = G.number_of_nodes()
    embed = nx.check_planarity(G)[1]
    pos = nx.combinatorial_embedding_to_pos(embed)
    ham = hamiltonian_cycle(G)

    cycle, inside, outside = classify_edges(G, pos, ham)

    def add_arc(a, b, down=False):
        h = abs(a-b)
        c = (a + b)/2
        if down:
            plt.gca().add_patch(mpl.patches.Arc((c, 0), h, h, theta1=180))
        else:
            plt.gca().add_patch(mpl.patches.Arc((c, 0), h, h, theta2=180))

    upper_seq = set()
    upper_x_map = {x:i for i, x in enumerate(ham)}
    lower_x_map = {x:i for i, x in enumerate(target_ordering)}
    inter_seg = set()
    lower_seq = set()

    # add_arc(0, n-1)
    upper_seq.add((n-1, 0))
    for i in range(n-1):
        upper_seq.add((i+1, i))

    for a, b in inside:
        l = min(upper_x_map[a], upper_x_map[b])
        r = max(upper_x_map[a], upper_x_map[b])
        upper_seq.add((l, r))

    for a, b in outside:
        l = min(upper_x_map[a], upper_x_map[b])
        r = max(upper_x_map[a], upper_x_map[b])
        upper_seq.add((r, l))

    # Let's goooo
    for i in range(n):

        v = ham[i]
        dest = lower_x_map[v]
        print(f"Moving vertex {ham[i]} from {i} to {dest}")
        if plot:
            linear_order_draw_plot(ham, upper_seq, inter_seg, lower_seq, target_ordering)

        # Step 1, move the trajectory edges
        print('is', inter_seg)
        upper_right_map = {a: b for (a, b) in upper_seq if a < b}
        upper_left_map = {b: a for (a, b) in upper_seq if a < b}
        Il = {(a, b) for a, b in inter_seg if a > i and b < dest and (a in range(n) or upper_left_map[a] != i)}
        Ir = {(a, b) for a, b in inter_seg if a < i and b > dest and upper_right_map[a] != i}
        Oa = {(a, b) for a, b in lower_seq if a < dest < b}
        assert not (Il and Ir) # otherwise not planar

        intermediates = len(Il) + len(Ir) + len(Oa)
        all_lower = [p for e in lower_seq for p in e] + list(range(n))
        t1prime = max(x for x in all_lower if x < dest) if dest != 0 else -1
        t2prime = min(x for x in all_lower if x > dest) if dest != n-1 else n
        print('Il', Il)
        print('Ir', Ir)
        print("t'", t1prime, "t''", t2prime)

        p = [dest - (dest-t1prime)/(intermediates + 1) * (k+1) for k in range(intermediates)]
        q = [dest + (t2prime-dest)/(intermediates + 1) * (k+1) for k in range(intermediates)]
        print('p', p, 'q', q)

        lower_seq -= Oa
        for j, (a, b) in enumerate(sorted(Oa, key=lambda ab:ab[1]-ab[0])):
            pa = p[-j-1]
            qb = q[-j-1]
            lower_seq.add((a, pa))
            lower_seq.add((qb, pa))
            lower_seq.add((qb, b))

        lower_left_map = {a: b for (a, b) in lower_seq if a > b}
        lower_right_map = {b: a for (a, b) in lower_seq if a > b}
        if Il:
            if not Oa and False:
                # haha
                inter_seg -= Il

                for j, (a, b) in enumerate(sorted(Il, key=lambda ab:ab[0]+ab[1])):
                    pb = p[j]
                    qa = q[j]
                    if b in lower_left_map:
                        c = lower_left_map[b]
                        lower_seq.remove((b, c))
                        lower_seq.add((qa, c))
                        inter_seg.add((a, qa))
                    else:
                        lower_seq.add((qa, b))
                        # lower_seq.add((b, pb))
                        # lower_seq.add((qa, pb))
                        inter_seg.add((a, qa))
            else:
                inter_seg -= Il
                for j, (a, b) in enumerate(sorted(Il, key=lambda ab:ab[0]+ab[1])):
                    pb = p[j]
                    qa = q[j]
                    lower_seq.add((b, pb))
                    lower_seq.add((qa, pb))
                    inter_seg.add((a, qa))

        elif Ir:
            if not Oa and False:
                inter_seg -= Ir
                for j, (a, b) in enumerate(sorted(Ir, key=lambda ab:-ab[0]-ab[1])):
                    pa = p[j]
                    qb = q[j]
                    if b in lower_right_map:
                        c = lower_right_map[b]
                        lower_seq.remove((c, b))
                        inter_seg.add((a, pa))
                        lower_seq.add((c, pa))
                    else:
                        inter_seg.add((a, pa))
                        lower_seq.add((b, pa))

            else:
                inter_seg -= Ir
                for j, (a, b) in enumerate(sorted(Ir, key=lambda ab:-ab[0]-ab[1])):
                    pa = p[j]
                    qb = q[j]
                    inter_seg.add((a, pa))
                    lower_seq.add((qb, pa))
                    lower_seq.add((qb, b))

        if plot:
            linear_order_draw_plot(ham, upper_seq, inter_seg, lower_seq, target_ordering)
        # Step 2, move the incident edges
        all_upper = [p for e in upper_seq for p in e]
        Ebru = {(a, b) for (a, b) in upper_seq if a > b == i}
        upper_seq -= Ebru
        for a, _ in Ebru:
            inter_seg.add((a, dest))

        # Also covers Ebl
        Ei = {(a, b) for (a, b) in inter_seg if a == i}
        print('Ei', Ei)
        inter_seg -= Ei
        for _, b in Ei:
            if b < dest:
                lower_seq.add((b, dest))
            else:
                lower_seq.add((dest, b))

        Etr = {(a, b) for (a, b) in upper_seq if i == a < b}
        upper_seq -= Etr
        sprime = min(x for x in all_upper if x > i) if i != n-1 else n
        print('sprime', sprime)
        print('Etr', Etr)
        p = [i + (sprime-i)/(len(Etr) + 1) * (k+1) for k in range(len(Etr))]
        print('p', p)
        for j, (_, b) in enumerate(sorted(Etr, key=lambda ab:-ab[1])):
            upper_seq.add((p[j], b))
            inter_seg.add((p[j], dest))

        Etl = {(a, b) for (a, b) in upper_seq if a < b == i}
        print("Etl", Etl)
        upper_seq -= Etl
        for a, _ in Etl:
            # find the real b
            inter_map = {ai: bi for ai, bi in inter_seg}
            b = inter_map[a]
            inter_seg.remove((a, b))
            if b < dest:
                lower_seq.add((b, dest))
            else:
                lower_seq.add((dest, b))

    if plot:
        linear_order_draw_plot(ham, upper_seq, inter_seg, lower_seq, target_ordering)
    return lower_seq

def gauss_target_order(gcode):
    n = len(gcode) // 2
    target_order = []
    cross_seen = set()
    g_i = 0
    for c in gcode:
        cross, o, p = lp.decompose_c(c)
        if cross in cross_seen:
            continue

        cross_seen.add(cross)

        # l u c d r
        if o == 1:
            gadget_order = [4, 3, 2, 1, 0]
            # if p == 1:
            #     gadget_order = [4, 3, 2, 1, 0]
            # else:
            #     gadget_order = [4, 1, 2, 3, 0]
        else:
            if p == 1:
                # comes out on l
                gadget_order = [3, 4, 2, 0, 1]
            else:
                gadget_order = [1, 0, 2, 4, 3]

        target_order.extend([g_i*5 + o for o in gadget_order])
        g_i += 1

    return target_order

def combine_arcs(arcs, target_ordering):
    edges = [(*arc, [arc]) for arc in arcs]
    # TODO: make this fast by being less dumb
    while True:
        for (i1, e1), (i2, e2) in it.combinations(enumerate(edges), 2):
            a1, b1, arcs1 = e1
            a2, b2, arcs2 = e2
            common = set((a1, b1)) & set((a2, b2))
            if common and len(common) == 1 and next(iter(common)) not in target_ordering:
                a3, b3 = set((a1, b1)) ^ set((a2, b2))
                arcs3 = arcs1 + arcs2
                # careful, pop the last one first
                edges.pop(i2)
                edges.pop(i1)
                edges.append((a3, b3, arcs3))
                break
        else:
            break

    return edges

def clean_edges(edges):
    resulting_edges = []
    resulting_edges = edges
    # for (a, b, arcs) in edges:
    #     if a//5 != b//5:
    #         resulting_edges.append((a, b, arcs))

    "pop"
    while True:
        found = False
        for (e1, e2, arcs) in resulting_edges:
            down_right = {b: (i, a) for i, (a, b) in enumerate(arcs) if a > b}
            down_left = {a: (i, b) for i, (a, b) in enumerate(arcs) if a > b}
            for i, (a, b) in enumerate(arcs):
                if a%1 != 0 and b%1 != 0 and a < b and b-a < 1 and a in down_left and b in down_right:
                    j, c = down_left[a]
                    k, d = down_right[b]
                    # if c%1 == 0 or d%1 == 0: continue
                    for l in sorted([i, j, k], reverse=True):
                        arcs.pop(l)
                    arcs.append((d, c))
                    print(f"Popped {c} {a} {b} {d} to {c} {d}")
                    found = True
                    break

                if a%1 == 0 and b%1 != 0 and a < b and b-a < 1 and b in down_right:
                    j, c = down_right[b]
                    for l in sorted([i, j], reverse=True):
                        arcs.pop(l)
                    arcs.append((c, a))
                    found=True
                    break

                if b%1 == 0 and a%1 != 0 and a < b and b-a < 1 and a in down_left:
                    j, c = down_left[a]
                    for l in sorted([i, j], reverse=True):
                        arcs.pop(l)
                    arcs.append((b, c))
                    found=True
                    break

        if not found: break

    return resulting_edges

def flatten_edges(edges):
    resulting_arcs = []
    for (a, b, arcs) in edges:
        resulting_arcs.extend(arcs)
    return resulting_arcs

def flip_vertex(edges, v, target_order):
    x = target_order.index(v)
    print(f"Flipping vertex {v} at position {x}")

    for (e1, e2, arcs) in edges:
        up_right = {a: (i, b) for i, (a, b) in enumerate(arcs) if a < b}
        up_left = {b: (i, a) for i, (a, b) in enumerate(arcs) if a < b}
        for i, (a, b) in enumerate(arcs):
            if x in (e1, e2) and len(arcs) > 1:
                print(f"Incident edge {e1} {e2} will be replaced")
                print(arcs)
                arcs.clear()
                arcs.append((max(e1, e2), min(e1, e2)))
                print(arcs)
                break
            if b < x < a:
                if a - b < 1:
                    print(f"Underskip {a} {b} will be fused")
                    j, c = up_left[b]
                    k, d = up_right[a]
                    for l in sorted([i, j, k], reverse=True):
                        arcs.pop(l)
                    arcs.append((c, d))
                    break

                elif a - b <= 2:
                    print(f"Underpass {a} {b} will be flipped")

                else:
                    assert False

def zerotone_advance_vertex(ham, arcs, src, dst):
    v = ham[src]
    print(f"Moving vertex {v} from {src} to {dst}")
    assert dst < src
    incident = []
    arcs = set(arcs)
    for a, b in arcs:
        if a == src:
            incident.append(b)
        elif b == src:
            incident.append(a)

    for i in incident:
        arcs.discard((i, src))
        arcs.discard((src, i))
    print(f"Incident edges: {incident}")

    # shift
    new_arcs = set()

    def interval_inc(i):
        if dst <= i < src:
            return i+1
        return i

    for a, b in arcs:
        new_arcs.add((interval_inc(a), interval_inc(b)))

    new_ham = deepcopy(ham)
    assert new_ham.pop(src) == v
    print(v)
    new_ham.insert(dst, v)
    for i in incident:
        if i >= dst:
            new_arcs.add((dst, interval_inc(i)))
        else:
            new_arcs.add((i, dst))

    return new_ham, new_arcs

if __name__ == "__main__":
    # for t, knot in gknot.items():
    #     print('=' * 10 + str(t) + '=' * 10)
    #     # knot = gknot[(11, 42)]
    #     G, knot_order = build_embedding_graph(lp.normalize_gauss_order(knot))
    #     print(hamiltonian_cycle(G))
    #     # assert nx.is_biconnected(G)
    #     # assert nx.is_k_edge_connected(G, 4)

    #     # embed = nx.check_planarity(G)[1]
    #     # pos = nx.combinatorial_embedding_to_pos(embed)

    knot = gknot[(8, 10)]
    G, knot_order = build_embedding_graph(lp.normalize_gauss_order(knot))
    print(knot_order)
    embed = nx.check_planarity(G)[1]
    pos = nx.combinatorial_embedding_to_pos(embed)
    ham = hamiltonian_cycle(G)
    print(ham)
    # cycle, inside, outside = classify_edges(G, pos, ham)
    # print(cycle, inside, outside)
    # ham_plot(G, cycle, inside, outside)

    # draw_zerotone_book_embedding(G, ham)
    zarcs = zerotone_book_embedding_arcs(G)

    # ham2, arcs2 = deepcopy(ham), deepcopy(zarcs)
    # while True:
    #     draw_zerotone_book_embedding(arcs2, ham2)
    #     inp = input(">>> ")
    #     if not inp: break
    #     src, dst = map(int, inp.strip().split(','))
    #     src = ham2.index(src)
    #     dst = ham2.index(dst)
    #     ham2, arcs2 = zerotone_advance_vertex(ham2, arcs2, src, dst)

    target_ordering = gauss_target_order(knot)
    # target_ordering = [0, 2, 1, 9, 7, 6, 29, 27, 26, 35, 37, 36, 30, 32, 31, 25, 28, 34, 33, 14, 12, 11, 19, 17, 16, 39, 38, 24, 22, 21, 4, 3, 20, 23, 15, 18, 10, 13, 5, 8][::-1]
    # target_ordering = [3,0,2,4,1,9,8,7,6,5,13,10,12,14,11]
    # target_ordering = [0,2,1,4,3,9,7,6,5,8,14,12,11,10,13]
    # print(target_ordering)
    # import random
    # target_ordering = list(G.nodes)
    # random.shuffle(target_ordering)
    arcs = linear_order_draw(G, target_ordering, plot=True)
    edges = combine_arcs(arcs, target_ordering)
    cleaned_edges = clean_edges(edges)
    cleaned_arcs = flatten_edges(cleaned_edges)
    linear_order_draw_plot(ham, [], [], cleaned_arcs, target_ordering, hide_upper_spine=True)
    # flip_vertex(cleaned_edges, 8, target_ordering)
    # cleaned_arcs = flatten_edges(cleaned_edges)
    # linear_order_draw_plot(ham, [], [], cleaned_arcs, target_ordering, hide_upper_spine=True)
    # linear_order_draw(G, [14, 13, 1, 3, 5, 4, 10, 8, 7, 11, 6, 9, 2, 0, 12])

    # V, E = semicircle_embed(Gm)
    # print(V, E)
    # assert nx.is_biconnected(G)
    # planar_plot(G)
    # plt.show()
