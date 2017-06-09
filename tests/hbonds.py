# -*- coding: utf8 -*-
from testrun import cfiles
import os

TRAJECTORY = os.path.join(os.path.dirname(__file__), "data", "water.xyz")
OUTPUT = "tmp.dat"
EXPECTED = [
    (39, 0, 1), (270, 0, 2), (69, 3, 4), (75, 3, 5), (228, 6, 8), (66, 9, 10),
    (3, 9, 11), (24, 18, 19), (27, 18, 20), (42, 21, 22)
]


def read_data(path):
    data = []
    indexes = []
    frame = True
    with open(path) as fd:
        for line in fd:
            if line.startswith("#"):
                if line.startswith("# Frame: 1"):
                    frame = False
                continue
            dist, angle = map(float, line.split()[-2:])
            data.append((dist, angle))
            if (frame):
                donor = int(line.split()[1])
                acceptor = int(line.split()[3])
                hydrogen = int(line.split()[5])
                indexes.append((donor, acceptor, hydrogen))
    return data, indexes


def check_angles(data):
    # Check the maximal value
    max_value = max(data, key=lambda u: u[0])
    assert(max_value[0] <= 3)
    max_value = max(data, key=lambda u: u[1])
    assert(max_value[1] <= 30)
    assert(max_value[1] >= 0)


def check_hbonds(indexes):
    # Check the first hbonds
    for bond in EXPECTED:
        assert(bond in indexes)


def hbonds():
    out, err = cfiles(
        "hbonds",
        "--guess-bonds",
        "-c", "15",
        TRAJECTORY, "-o", OUTPUT
    )
    assert(out == "")
    assert(err == "")

    data, indexes = read_data(OUTPUT)
    check_angles(data)
    check_hbonds(indexes)


if __name__ == '__main__':
    hbonds()
    os.unlink(OUTPUT)
