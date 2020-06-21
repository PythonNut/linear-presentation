import os
import re
from urllib import request
from tqdm import tqdm


def OU_to_nums(gc):
    """
    convert an input gauss code of the form
    (str):  "O1-O2-U1-U2-..."
    to the nelson gauss code format
    (list): [1.5, 2.5, -1.5, -2.5, ...]
    """
    new_gc = []
    segs = [gc[3 * i : 3 * i + 3] for i in range(len(gc) // 3)]
    for seg in segs:
        n = int(seg[1])
        if seg[2] == "-":
            n += 0.5
        if seg[0] == "U":
            n *= -1
        new_gc += [n]
    return new_gc


def url(key):
    (n, cnum) = key
    return f"https://www.math.toronto.edu/drorbn/Students/GreenJ/{n}.{cnum}.html"


def get_from_file(fname, write_name="data/vknot_full.py"):
    """
    Results are given at
    https://www.math.toronto.edu/drorbn/Students/GreenJ/results.html
    in the form of a gzipped archive. Download it, use `gunzip` to
    extract, and run this function on the result to convert the OU
    codes to nelson codes.
    """
    with open(fname, "r") as f:
        codes = f.read()
    codes = [code for code in codes.split("\n") if ("O" in code) or ("U" in code)]

    vknot = {}
    k_index = 0
    cnum = 0
    for code in tqdm(codes):
        nelson_code = OU_to_nums(code)
        new_cnum = len(nelson_code) // 2

        if new_cnum != cnum:
            k_index = 1
            cnum = new_cnum
        else:
            k_index += 1

        vknot[(new_cnum, k_index)] = nelson_code

    with open(write_name, "w") as f:
        f.write(f"vknot_full = {vknot}")

    return vknot


def pull_from_website():
    vknot = {}
    keys = [(2, 1)] + [(3, i) for i in range(1, 8)] + [(4, i) for i in range(1, 109)]

    for key in tqdm(keys):
        # Get the url
        k_url = url(key)

        # Fetch the html contents as a string
        k_page = request.urlopen(k_url).read().decode("utf-8")

        gcstr = re.findall(r"Gauss code</a></dt><dd>[OU\d\+\-]*</dd>", k_page)[0]
        gc = gcstr[23:-5]
        gc = OU_to_nums(gc)
        vknot[key] = gc

    with open("vcodes.py", "w") as f:
        f.write("vknot = " + str(vknot))
