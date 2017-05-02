# -*- coding: utf8 -*-
from testrun import cfiles
import os

TRAJECTORY = os.path.join(os.path.dirname(__file__), "data", "nt.xyz")
OUTPUT = "tmp.dat"


def read_data(path):
    data = []
    with open(path) as fd:
        for line in fd:
            if line.startswith("#"):
                continue
            dist, value = map(float, line.split())
            data.append((dist, value))
    return data


def check_density(data):
    for (radius, value) in data:
        # Check the inside of the nanotube is empty
        if (radius < 5.0):
            assert(value == 0)
        # Check density is zero between Oint and Si
        if (radius > 7.1 and radius < 7.7):
            assert(value == 0)

def check_max(data):
    # Check the maximal value
    max_value = max(data, key=lambda u: u[1])
    assert(max_value[0] > 9.5)
    assert(max_value[0] < 9.9)

def density(selection):
    out, err = cfiles(
        "density",
        "-c", "24:24:25.458:90:90:120",
        "--axis=Z",
        "--max=20",
        "--points=200",
        "--profile=radial",
        "-s", selection,
        TRAJECTORY, "-o", OUTPUT
    )
    assert(out == "")
    assert(err == "")

    data = read_data(OUTPUT)
    check_density(data)
    return data

if __name__ == '__main__':
    tot = density("atoms: all")
    al = density("atoms: type Al")
    si = density("atoms: type Si")
    o = density("atoms: type O")
    h = density("atoms: type H")
    # Check tot is the sum of all the elements
    for radius in range(200):
        assert(tot[radius][1] - (al[radius][1] + si[radius][1] + o[radius][1] + h[radius][1]) < 0.1)
    os.unlink(OUTPUT)
