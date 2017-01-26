# -*- coding: utf8 -*-
from testrun import cfiles
import os

TRAJECTORY = os.path.join(os.path.dirname(__file__), "data", "water.xyz")
OUTPUT = "tmp.rdf"


def read_rdf(path):
    data = []
    with open(path) as fd:
        for line in fd:
            if line.startswith("#"):
                continue
            r, g_r, n_r = map(float, line.split())
            data.append((r, g_r, n_r))
    return data


def check_oxygen_rdf(data):
    # Check the maximal value of the rdf
    max_value = max(data, key=lambda u: u[1])
    assert(max_value[0] == 2.75)
    assert(max_value[1] > 3)
    # Check that the coordination number in the first sphere is around 2
    assert(abs(max_value[2] - 2) < 1)

    # Check the last zero value
    last_zero = [u for u in data if u[1] == 0 and u[0] < max_value[0]][-1]
    assert(last_zero[0] == 2.45)

    # Check that the rdf converges to 1
    end = [g[1] for g in data[len(data)/2:]]
    ave = sum(end) / len(end)
    assert(abs(ave - 1) < 1e-1)

    # Check that the coordination number goes to N ⪞ 99
    n_max = data[-1][2]
    assert(abs(n_max - 99) < 3)


def check_oh_rdf(data):
    # Check the maximal value of the rdf
    max_value = max(data, key=lambda u: u[1])
    assert(max_value[0] == 0.95)
    assert(max_value[1] > 25)
    # Check that the coordination number in the first sphere is 3
    assert(abs(max_value[2] - 3) < 0.1)

    # Check the last zero value
    last_zero = [u for u in data if u[1] == 0 and u[0] < max_value[0]][-1]
    assert(last_zero[0] == 0.9)

    # Check that the rdf converges to 1
    end = [g[1] for g in data[len(data)/2:]]
    ave = sum(end) / len(end)
    assert(abs(ave - 1) < 1e-1)

    # Check that the coordination number goes to N ⪞ 148.5
    n_max = data[-1][2]
    assert(abs(n_max - 148.5) < 7)


def oxygen_rdf_all():
    '''Oxygen rdf for the whole trajectory'''
    out, err = cfiles(
        "rdf",
        "-c", "15",          # Set cell
        "-p", "150",         # Use 150 points in the histogram
        "-s", "name O",      # Compute rdf between O
        TRAJECTORY, "-o", OUTPUT
    )
    assert(out == "")
    assert(err == "")

    data = read_rdf(OUTPUT)
    check_oxygen_rdf(data)


def oxygen_rdf_partial():
    '''Oxygen rdf for the second half of the trajectory'''
    out, err = cfiles(
        "rdf",
        "--start", "50",
        "-c", "15",          # Set cell
        "-p", "150",         # Use 150 points in the histogram
        "-s", "name O",      # Compute rdf between O
        TRAJECTORY, "-o", OUTPUT
    )
    assert(out == "")
    assert(err == "")

    data = read_rdf(OUTPUT)
    check_oxygen_rdf(data)


def OH_rdf_all():
    '''Oxygen-Hydrogen rdf for the whole trajectory'''
    out, err = cfiles(
        "rdf",
        "-c", "15",
        "-p", "150",
        "-s", "pairs: name(#1) O and name(#2) H",
        TRAJECTORY, "-o", OUTPUT
    )
    assert(out == "")
    assert(err == "")

    data = read_rdf(OUTPUT)
    check_oh_rdf(data)


def OH_rdf_partial():
    '''Oxygen-Hydrogen rdf for the second half of the trajectory'''
    out, err = cfiles(
        "rdf",
        "--start", "50",
        "-c", "15",
        "-p", "150",
        "-s", "pairs: name(#1) O and name(#2) H",
        TRAJECTORY, "-o", OUTPUT
    )
    assert(out == "")
    assert(err == "")

    data = read_rdf(OUTPUT)
    check_oh_rdf(data)

    out, err = cfiles(
        "rdf",
        "--start", "50",
        "-c", "15",
        "-p", "150",
        "-s", "pairs: name(#1) H and name(#2) O",
        TRAJECTORY, "-o", OUTPUT
    )
    assert(out == "")
    assert(err == "")

    data = read_rdf(OUTPUT)
    check_oh_rdf(data)


if __name__ == '__main__':
    oxygen_rdf_all()
    oxygen_rdf_partial()
    OH_rdf_all()
    OH_rdf_partial()
    os.unlink(OUTPUT)
