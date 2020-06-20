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
