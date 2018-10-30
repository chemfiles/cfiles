# -*- coding: utf8 -*-
from testrun import cfiles
import tempfile
import os

TRAJECTORY = os.path.join(os.path.dirname(__file__), "data", "water.xyz")


def read_rdf(path):
    data = []
    with open(path) as fd:
        for line in fd:
            if line.startswith("#"):
                continue
            r, g_r, n_ij, n_ji = map(float, line.split())
            data.append((r, g_r, n_ij, n_ji))
    return data


def check_oxygen_rdf(data):
    # Check the maximal value of the rdf
    max_index = max(enumerate(data), key=lambda u: u[1][1])[0]
    max_value = data[max_index]
    assert(max_value[0] == 2.775)
    assert(max_value[1] > 3)
    # Check that the coordination number in the first sphere is around 2.6
    assert(abs(data[max_index + 3][2] - 2.6) < 0.2)
    assert(abs(data[max_index + 3][3] - 2.6) < 0.2)

    # Check the last zero value
    last_zero = [u for u in data if u[1] == 0 and u[0] < max_value[0]][-1]
    assert(last_zero[0] == 2.475)

    # Check that the rdf converges to 1
    end = [g[1] for g in data[len(data)/2:]]
    ave = sum(end) / len(end)
    assert(abs(ave - 1) < 1e-1)


def check_oh_rdf(data):
    # Check the maximal value of the rdf
    max_index = max(enumerate(data), key=lambda u: u[1][1])[0]
    max_value = data[max_index]
    assert(max_value[0] == 0.975)
    assert(max_value[1] > 25)
    # Check the coordination number in the first sphere
    assert(abs(data[max_index + 3][2] - 2) < 0.01)
    assert(abs(data[max_index + 3][3] - 1) < 0.01)

    # Check the last zero value
    last_zero = [u for u in data if u[1] == 0 and u[0] < max_value[0]][-1]
    assert(last_zero[0] == 0.925)

    # Check that the rdf converges to 1
    end = [g[1] for g in data[len(data)/2:]]
    ave = sum(end) / len(end)
    assert(abs(ave - 1) < 1e-1)


def check_ho_rdf(data):
    # Check the maximal value of the rdf
    max_index = max(enumerate(data), key=lambda u: u[1][1])[0]
    max_value = data[max_index]
    assert(max_value[0] == 0.975)
    assert(max_value[1] > 25)
    # Check the coordination number in the first sphere
    assert(abs(data[max_index + 3][2] - 1) < 0.01)
    assert(abs(data[max_index + 3][3] - 2) < 0.01)

    # Check the last zero value
    last_zero = [u for u in data if u[1] == 0 and u[0] < max_value[0]][-1]
    assert(last_zero[0] == 0.925)

    # Check that the rdf converges to 1
    end = [g[1] for g in data[len(data)/2:]]
    ave = sum(end) / len(end)
    assert(abs(ave - 1) < 1e-1)


def oxygen_rdf_all(output):
    '''Oxygen rdf for the whole trajectory'''
    out, err = cfiles(
        "rdf",
        "-c", "15",          # Set cell
        "-p", "150",         # Use 150 points in the histogram
        "-s", "name O",      # Compute rdf between O
        TRAJECTORY, "-o", output
    )
    assert(out == "")
    assert(err == "")

    data = read_rdf(output)
    check_oxygen_rdf(data)


def oxygen_rdf_partial(output):
    '''Oxygen rdf for the second half of the trajectory'''
    out, err = cfiles(
        "rdf",
        "--steps", "50:",
        "-c", "15",          # Set cell
        "-p", "150",         # Use 150 points in the histogram
        "-s", "name O",      # Compute rdf between O
        TRAJECTORY, "-o", output
    )
    assert(out == "")
    assert(err == "")

    data = read_rdf(output)
    check_oxygen_rdf(data)


def OH_rdf_all(output):
    '''Oxygen-Hydrogen rdf for the whole trajectory'''
    out, err = cfiles(
        "rdf",
        "-c", "15",
        "-p", "150",
        "-s", "pairs: name(#1) O and name(#2) H",
        TRAJECTORY, "-o", output
    )
    assert(out == "")
    assert(err == "")

    data = read_rdf(output)
    check_oh_rdf(data)

    out, err = cfiles(
        "rdf",
        "-c", "15",
        "-p", "150",
        "-s", "pairs: name(#1) H and name(#2) O",
        TRAJECTORY, "-o", output
    )
    assert(out == "")
    assert(err == "")

    data = read_rdf(output)
    check_ho_rdf(data)


def OH_rdf_partial(output):
    '''Oxygen-Hydrogen rdf for half of the trajectory'''
    out, err = cfiles(
        "rdf",
        "--steps", "::2",
        "-c", "15",
        "-p", "150",
        "-s", "pairs: name(#1) O and name(#2) H",
        TRAJECTORY, "-o", output
    )
    assert(out == "")
    assert(err == "")

    data = read_rdf(output)
    check_oh_rdf(data)

    out, err = cfiles(
        "rdf",
        "--steps", ":50",
        "-c", "15",
        "-p", "150",
        "-s", "pairs: name(#1) H and name(#2) O",
        TRAJECTORY, "-o", output
    )
    assert(out == "")
    assert(err == "")

    data = read_rdf(output)
    check_ho_rdf(data)


if __name__ == '__main__':
    with tempfile.NamedTemporaryFile() as file:
        oxygen_rdf_all(file.name)
        oxygen_rdf_partial(file.name)
        OH_rdf_all(file.name)
        OH_rdf_partial(file.name)
