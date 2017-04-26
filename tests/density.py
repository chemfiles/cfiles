# -*- coding: utf8 -*-
from testrun import cfiles
import os

TRAJECTORY = os.path.join(os.path.dirname(__file__), "data", "water.xyz")
OUTPUT = "tmp.dat"


def read_data(path):
    data = []
    with open(path) as fd:
        for line in fd:
            if line.startswith("#"):
                continue
            phi, value = map(float, line.split())
            data.append((phi, value))
    return data


def check_density(data):
    # Check the maximal value
    max_value = max(data, key=lambda u: u[1])
    assert(max_value[0] < 7.5)
    assert(max_value[1] == 1)

def density(selection):
    out, err = cfiles(
        "density",
        "-c", "15",
        "--axis=Z",
        "-s", selection,
        TRAJECTORY, "-o", OUTPUT
    )
    assert(out == "")
    assert(err == "")

    data = read_data(OUTPUT)
    check_density(data)


if __name__ == '__main__':
    density("atoms: all")
    density("atoms: name(#1) O")
    os.unlink(OUTPUT)
