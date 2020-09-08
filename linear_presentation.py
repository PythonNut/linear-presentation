from tqdm import tqdm
from multiprocessing import Pool

# Import all the routing functions
from routing import *

# Ok we really gotta clean these things up at some point
from gauss_codes import get_cknots, get_vknots, nelson_gc_to_sage_gc

c_dict = get_cknots(max_cnum=11)
v_dict = get_vknots(max_cnum=5)

# Drawing
from tikz_plotter import plot_virtual, make_mosaic

######################################################################
#                                                                    #
#                 Bulk processing w/ multithreading                  #
#                                                                    #
######################################################################
def proc_single_cknot(key):
    # Get the Gauss code for the knot
    kcode = list(c_dict[key])

    # The knot dict already has in sage gc format
    K = Knot(kcode)
    crossings, semiarcs = knot_to_layout(K)

    # Experimental: Use virtual routing for classical knots as well
    try:
        upper_cs, lower_cs, crossings = route(semiarcs)
    except AssertionError as e:
        print(key)
        raise (e)

    [gc], orient = kcode
    cn, kind = key
    dirname = f"{cn}-{kind}"
    try:
        plot_virtual(
            gc,
            orient,
            upper_cs,
            lower_cs,
            crossings,
            dirname,
            dir_prefix="example-execution/classicals/",
            use_fk_colors=False,
            warn_classical=False,
        )
    except ZeroDivisionError as e:
        print(key)
        raise (e)


def proc_single_vknot(key):
    # Get the Gauss code for the knot
    gc = v_dict[key]
    K = Knot(nelson_gc_to_sage_gc(gc))
    crossings, semiarcs = knot_to_layout(K)

    try:
        upper_cs, lower_cs, crossings = virtual_route(semiarcs)
    except AssertionError as e:
        print(key)
        raise (e)

    [gc], orient = nelson_gc_to_sage_gc(gc)
    cn, kind = key
    dirname = f"{cn}-{kind}"
    try:
        plot_virtual(
            gc, orient, upper_cs, lower_cs, crossings, dirname, use_fk_colors=False,
        )
    except ZeroDivisionError as e:
        print(key)
        raise (e)


def proc_all(flavor="classicals"):
    if flavor == "classicals":
        knot_dict = get_cknots(max_cnum=11)
        proc_func = proc_single_cknot
    elif flavor == "virtuals":
        knot_dict = get_vknots(max_cnum=5)
        proc_func = proc_single_vknot

    _n = len(knot_dict)
    with Pool(7) as p:
        with tqdm(total=_n) as pbar:
            for i, _ in enumerate(p.imap_unordered(proc_func, knot_dict)):
                pbar.update()

    make_mosaic(flavor)


# if __name__ == "__main__":
#    from tikz_plotter import plot_virtual, make_mosaic
#    proc_all(flavor="virtuals")
