def convert_gc(gc):
    """
    Convert from the gauss code standard
    """

    return

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

# def convert_

gknot = {
  (0):[-1,1],
(3,1):[-1,2,-3,1,-2,3],
(4,1):[-1.5,2.5,-3,4,-2.5,1.5,-4,3],
(5,1):[-1,2,-3,4,-5,1,-2,3,-4,5],
(5,2):[-1,2,-3,4,-5,3,-2,1,-4,5],
(6,1):[-1,2,-3.5,4.5,-5.5,6.5,-2,1,-6.5,5.5,-4.5,3.5],
(6,2):[-1,2.5,-3.5,4.5,-5.5,1,-6,3.5,-4.5,5.5,-2.5,6],
(6,3):[1,-2.5,3.5,-4.5,2.5,-5,6,-1,5,-3.5,4.5,-6],
(6,4):[-1,2,-3,1,-2,3,-4,5,-6,4,-5,6],
(6,5):[-1,2,-3,1,-2,3,4.5,-5.5,6.5,-4.5,5.5,-6.5],
(7,1):[-1,2,-3,4,-5,6,-7,1,-2,3,-4,5,-6,7],
(7,2):[-1.5,6.5,-7.5,1.5,-2.5,3.5,-4.5,5.5,-6.5,7.5,-5.5,4.5,-3.5,2.5],
(7,3):[-1,2,-3,4,-7,6,-5,1,-2,3,-4,5,-6,7],
(7,4):[-1,2,-3,4,-5,6,-7,1,-4,3,-2,7,-6,5],
(7,5):[-1.5,2.5,-3.5,4.5,-5.5,1.5,-2.5,3.5,-6.5,7.5,-4.5,5.5,-7.5,6.5],
(7,6):[-1.5,2.5,-3,4,-5.5,6.5,-4,3,-7.5,1.5,-6.5,5.5,-2.5,7.5],
(7,7):[-1.5,2.5,-3,4,-2.5,5.5,-6,7,-5.5,1.5,-4,3,-7,6],
(8,1):[-1.5,2.5,-3,4,-2.5,1.5,-5.5,6.5,-7.5,8.5,-4,3,-8.5,7.5,-6.5,5.5],
(8,2):[1,-2,3.5,-4.5,5.5,-6.5,7.5,-8.5,2,-1,8.5,-3.5,4.5,-5.5,6.5,-7.5],
(8,3):[-1,2,-3,4,-5.5,6.5,-7.5,8.5,-4,3,-2,1,-8.5,7.5,-6.5,5.5],
(8,4):[-1.5,2,-3,4,-5,1.5,-6.5,7.5,-8.5,5,-4,3,-2,6.5,-7.5,8.5],
(8,5):[-1.5,2.5,-3,4,-5,6,-7,8,-2.5,1.5,-6,7,-8,3,-4,5],
(8,6):[1,-2,3.5,-4.5,5.5,-6.5,7.5,-8.5,2,-1,8.5,-7.5,6.5,-3.5,4.5,-5.5],
(8,7):[-1.5,2.5,-3.5,1.5,-4,5,-6,7,-8,4,-2.5,3.5,-5,6,-7,8],
(8,8):[1,-2,3,-4.5,5.5,-6.5,4.5,-7,8,-5.5,6.5,-3,2,-1,7,-8],
(8,9):[-1.5,2.5,-3.5,4.5,-5,6,-7,8,-2.5,3.5,-4.5,1.5,-8,5,-6,7],
(8,10):[-1,2,-3,4,-5,1,-2,6.5,-7.5,3,-4,8.5,-6.5,7.5,-8.5,5],
(8,11):[1,-2,3.5,-4.5,5.5,-6.5,7.5,-8.5,2,-1,8.5,-5.5,4.5,-3.5,6.5,-7.5],
(8,12):[-1,2,-3,4,-5.5,6.5,-4,3,-7.5,8.5,-2,1,-8.5,7.5,-6.5,5.5],
(8,13):[1.5,-2.5,3,-4,5,-6,7,-3,8.5,-1.5,2.5,-8.5,4,-7,6,-5],
(8,14):[-1,2,-3.5,4.5,-5.5,6.5,-7.5,5.5,-4.5,8.5,-2,1,-6.5,7.5,-8.5,3.5],
(8,15):[1.5,-2.5,3.5,-4.5,5.5,-3.5,6.5,-7.5,8.5,-6.5,2.5,-1.5,7.5,-8.5,4.5,-5.5],
(8,16):[1.5,-2,3,-4.5,5.5,-6,2,-7.5,4.5,-5.5,8.5,-1.5,7.5,-3,6,-8.5],
(8,17):[-1.5,2.5,-3,4,-5,6,-2.5,7.5,-4,5,-8.5,1.5,-6,3,-7.5,8.5],
(8,18):[-1,2.5,-3.5,4,-5,6.5,-2.5,7,-4,8.5,-6.5,1,-7,3.5,-8.5,5],
(8,19):[-1,2,-3,-4,5,1,-2,-6,4,7,-8,-5,6,3,-7,8],
(8,20):[-1,2.5,3,-4,-5.5,1,6.5,-3,4,-7.5,8.5,5.5,-2.5,-6.5,7.5,-8.5],
(8,21):[-1,2.5,-3.5,-4.5,5.5,-6,-7.5,1,4.5,-5.5,8.5,7.5,-2.5,3.5,6,-8.5],
(9,2):[-1.5,2.5,-3.5,4.5,-5.5,6.5,-7.5,8.5,-9.5,1.5,-2.5,9.5,-8.5,7.5,-6.5,5.5,-4.5,3.5],
(9,24):[-1.5,2.5,-3,4.5,-5.5,6,-7,8,-2.5,1.5,-6,7,-8,9.5,-4.5,5.5,-9.5,3],
(9,32):[-1.5,2.5,-3,4,-5,6,-7,8,-2.5,9.5,-4,7,-8,3,-9.5,1.5,-6,5],
(10,132):[1.5,-2.5,3.5,-1.5,-4,5.5,6,-7.5,-8,4,9.5,-6,-10.5,8,2.5,-3.5,7.5,10.5,-5.5,-9.5],
(11,1):[-1,2.5,-3.5,-4,5,6.5,-2.5,7,-8,3.5,-6.5,1,-7,8,9.5,-10.5,11.5,-5,4,-9.5,10.5,-11.5],
(11,2):[-1,2.5,-3.5,-4.5,5.5,-6.5,7,-8,4.5,-5.5,6.5,9.5,-2.5,10,-11,3.5,-9.5,1,-10,11,8,-7],
(11,42):[-1,2,3.5,-4,5,-6,7,-3.5,8.5,-5,6,9.5,-10.5,11.5,-9.5,1,-2,10.5,-11.5,-7,4,-8.5]
}

# Values represent where to slice gauss code to get a nice presentation of the knot
cycle_rules = {
    (5,2):4,
    (6,2):1,
    (6,4):3,
    (6,5):3,
    (7,0):0,
    (7,2):2,
    (7,3):2,
    (7,4):4,
    (7,5):5,
    # (7,6):
    (8,1):2,
    (8,4):5,
    (8,6):2,
    (8,8):6
}

# def apply_cycle_rules(gcode_dict):
#     new_dict = {}
#     for key in gcode_dict:
#         knot = gcode_dict[key]
#         i
#         new_dict[key] =
