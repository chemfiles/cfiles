# -*- coding: utf8 -*-
from testrun import cfiles
import tempfile
import os

TRAJECTORY = os.path.join(os.path.dirname(__file__), "data", "water.xyz")


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
    # Check that the integral is 1
    dr = 0.9
    sum = 0
    for _, value in data:
        sum += dr * value
    assert abs(sum - 1) < 1e-5

    # Check that other values are roughlty zero
    for (theta, value) in data:
        if theta != 104.85:
            assert abs(value) < 0.02


def angles(selection, output):
    out, err = cfiles(
        "angles", "--guess-bonds", "-c", "15", "-s", selection, TRAJECTORY, "-o", output
    )
    assert out == ""
    assert err == ""

    data = read_data(output)
    check_angles(data)


if __name__ == "__main__":
    with tempfile.NamedTemporaryFile() as file:
        angles("angles: all", file.name)
        angles("angles: name(#1) H and name(#2) O and name(#3) H", file.name)
