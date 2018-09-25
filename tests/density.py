# -*- coding: utf8 -*-
from testrun import cfiles
import tempfile
import os

TRAJECTORY = os.path.join(os.path.dirname(__file__), "data", "nt.xyz")


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


def density(selection, output):
    out, err = cfiles(
        "density",
        "-c", "24:24:25.458:90:90:120",
        "--max=20",
        "--points=200",
        "--radial=Z",
        "-s", selection,
        TRAJECTORY, "-o", output
    )
    assert(out == "")
    assert(err == "")

    data = read_data(output)
    check_density(data)
    return data


if __name__ == '__main__':
    with tempfile.NamedTemporaryFile() as file:
        tot = density("atoms: all", file.name)
        al = density("atoms: type Al", file.name)
        si = density("atoms: type Si", file.name)
        o = density("atoms: type O", file.name)
        h = density("atoms: type H", file.name)

        # Check tot is the sum of all the elements
        for radius in range(200):
            sum = al[radius][1] + si[radius][1] + o[radius][1] + h[radius][1]
            assert(abs(tot[radius][1] - sum) < 1e-3)
