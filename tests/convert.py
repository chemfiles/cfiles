# -*- coding: utf8 -*-
from testrun import cfiles
import os

PDB_CONTENT = """CRYST1   15.000   15.000   15.000  90.00  90.00  90.00
ATOM      1  O       X   1       0.417   8.303  11.737  0.00  0.00
ATOM      2  H       X   1       1.320   8.480  12.003  0.00  0.00
ATOM      3  H       X   1       0.332   8.726  10.882  0.00  0.00
ATOM      4  O       X   2       9.733   8.950   5.501  0.00  0.00
ATOM      5  H       X   2      10.058   9.331   6.317  0.00  0.00
ATOM      6  H       X   2       8.922   9.426   5.320  0.00  0.00
ATOM      7  O       X   3       1.486   6.471   0.946  0.00  0.00
ATOM      8  H       X   3       1.497   5.859   1.682  0.00  0.00
ATOM      9  H       X   3       0.560   6.554   0.717  0.00  0.00
END
"""

XYZ_CONTENT = """9
Written by the chemfiles library
O 0.417 8.303 11.737
H 1.32 8.48 12.003
H 0.332 8.726 10.882
O 9.733 8.95 5.501
H 10.058 9.331 6.317
H 8.922 9.426 5.32
O 1.486 6.471 0.946
H 1.497 5.859 1.682
H 0.56 6.554 0.717
"""

XYZ_SEL_CONTENT = """3
Written by the chemfiles library
O 0.417 8.303 11.737
O 9.733 8.95 5.501
O 1.486 6.471 0.946
"""


class isolate_files(object):
    '''
    Mock the filesystem for convertion testing, by creating and removing files.
    This class must be used as a context manager, as in

        with isolate_files("name"):
            ...
    '''
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        try:
            os.unlink(self.name + ".pdb")
        except OSError:
            pass
        with open(self.name + ".pdb", "w") as fd:
            fd.write(PDB_CONTENT)

    def __exit__(self, type, value, traceback):
        try:
            os.unlink(self.name + ".pdb")
        except OSError:
            pass
        try:
            os.unlink(self.name + ".xyz")
        except OSError:
            pass


if __name__ == '__main__':
    filename = "convert"
    with isolate_files(filename):
        out, err = cfiles("convert", filename + ".pdb", filename + ".xyz")
        assert(out == "")
        assert(err == "")
        with open(filename + ".xyz") as fd:
            assert(fd.read() == XYZ_CONTENT)
    with isolate_files(filename):
        out, err = cfiles("convert", filename + ".pdb", filename + ".xyz","--selection", 'atoms: type O')
        assert(out == "")
        assert(err == "")
        with open(filename + ".xyz") as fd:
            assert(fd.read() == XYZ_SEL_CONTENT)
