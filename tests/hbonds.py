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
            dist, angle = map(float, line.split()[-2:])
            data.append((dist, angle))
    return data


def check_angles(data):
    # Check the maximal value
    max_value = max(data, key=lambda u: u[0])
    assert(max_value[0] <= 3)
    max_value = max(data, key=lambda u: u[1])
    assert(max_value[1] <= 30)
    assert(max_value[1] >= 0)


def hbonds():
    out, err = cfiles(
        "hbonds",
        "--guess-bonds",
        "-c", "15",
        TRAJECTORY, "-o", OUTPUT
    )
    assert(out == "")
    assert(err == "")

    data = read_data(OUTPUT)
    check_angles(data)


if __name__ == '__main__':
    hbonds()
    os.unlink(OUTPUT)
