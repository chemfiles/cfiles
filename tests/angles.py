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


def check_angles(data):
    # Check the maximal value
    max_value = max(data, key=lambda u: u[1])
    assert(max_value[0] == 104.85)
    assert(max_value[1] == 1)

    # Check that other values are roughlty zero
    for (theta, value) in data:
        if (theta != 104.85):
            assert(abs(value) < 0.02)


def angles(selection):
    out, err = cfiles(
        "angles",
        "--guess-bonds",
        "-c", "15",
        "-s", selection,
        TRAJECTORY, "-o", OUTPUT
    )
    assert(out == "")
    assert(err == "")

    data = read_data(OUTPUT)
    check_angles(data)


if __name__ == '__main__':
    angles("angles: all")
    angles("angles: name(#1) H and name(#2) O and name(#3) H")
    os.unlink(OUTPUT)
