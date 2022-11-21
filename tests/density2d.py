import os
import tempfile

from testrun import cfiles

TRAJECTORY = os.path.join(os.path.dirname(__file__), "data", "nt.xyz")


def read_data(path):
    data = []
    with open(path) as fd:
        for line in fd:
            if line.startswith("#"):
                continue
            z, r, value = map(float, line.split())
            data.append((z, r, value))
    return data


def check_density(data):
    for (z, r, value) in data:
        # Check the inside of the nanotube is empty
        if r < 5.0:
            assert value == 0
        # Check density is zero between Oint and Si
        if r > 7.1 and r < 7.7:
            assert value == 0
        # Check density is zero between Si
        if z > 1 and z < 3 and r < 8 and r > 7:
            assert value == 0


def check_max(data):
    # Check the maximal value
    max_value = max(data, key=lambda u: u[1])
    assert max_value[1] > 0


def density(selection, output):
    out, err = cfiles(
        "density",
        "-c",
        "24:24:25.458:90:90:120",
        "--axis=Z",
        "--radial=Z",
        "--max=20:25",
        "--min=-20:0",
        "--points=200",
        "--origin",
        "0:0:0",
        "-s",
        selection,
        TRAJECTORY,
        "-o",
        output,
    )
    assert out == ""
    assert err == ""

    data = read_data(output)
    check_density(data)
    check_max(data)
    return data


if __name__ == "__main__":
    with tempfile.NamedTemporaryFile() as file:
        tot = density("atoms: all", file.name)
        al = density("atoms: type Al", file.name)
        si = density("atoms: type Si", file.name)
        o = density("atoms: type O", file.name)
        h = density("atoms: type H", file.name)

        # Check tot is the sum of all the elements
        for point in range(200):
            sum = al[point][2] + si[point][2] + o[point][2] + h[point][2]
            assert abs(tot[point][2] - sum) < 1e-3
