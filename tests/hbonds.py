# -*- coding: utf8 -*-
from testrun import cfiles
import tempfile
import os

TRAJECTORY = os.path.join(os.path.dirname(__file__), "data", "water.xyz")
EXPECTED = [
    (0, 1, 39), (0, 2, 270), (3, 4, 69), (3, 5, 75), (6, 8, 228), (9, 10, 66),
    (9, 11, 3), (18, 19, 24), (18, 20, 27), (21, 22, 42)
]


def read_data(path):
    indexes = []
    with open(path) as fd:
        for line in fd:
            if line.startswith("#"):
                continue
            splitted = map(int, line.split())
            if len(splitted) == 2:
                continue
            donor, acceptor, hydrogen = splitted
            indexes.append((donor, acceptor, hydrogen))
    return indexes


def check_hbonds(indexes):
    print(indexes)
    # Check the first hbonds
    for bond in EXPECTED:
        assert(bond in indexes)


def hbonds(output):
    out, err = cfiles(
        "hbonds",
        "--guess-bonds",
        "-c", "15",
        TRAJECTORY, "-o", output
    )
    assert(out == "")
    assert(err == "")

    indexes = read_data(output)
    check_hbonds(indexes)


if __name__ == '__main__':
    with tempfile.NamedTemporaryFile() as file:
        hbonds(file.name)
