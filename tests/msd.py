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
            t, msd = map(float, line.split())
            data.append((t, msd))
    return data


def msd(output):
    out, err = cfiles(
        "msd",
        "-c", "15",
        "--unwrap",
        "--selection", "name O",
        TRAJECTORY, "-o", output
    )
    assert(out == "")
    assert(err == "")

    data = read_data(output)
    check_msd(data)


def msd_no_cell(output):
    out, err = cfiles(
        "msd",
        "--selection", "name O",
        TRAJECTORY, "-o", output
    )
    assert(out == "")
    assert(err == "")


def check_msd(data):
    # This is only a regression test, checking that the right output is
    # generated.
    expected = read_data(os.path.join(
        os.path.dirname(__file__), "data", "water.msd.dat"
    ))

    for ((r, msd), (exp_r, exp_msd)) in zip(data, expected):
        assert(r == exp_r)
        assert(abs((msd - exp_msd) / msd) < 1e-3)


if __name__ == '__main__':
    with tempfile.NamedTemporaryFile() as file:
        msd(file.name)

    with tempfile.NamedTemporaryFile() as file:
        msd_no_cell(file.name)
