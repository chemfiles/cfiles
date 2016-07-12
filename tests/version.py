# -*- coding: utf8 -*-
from testrun import cfiles


if __name__ == '__main__':
    stdout, stderr = cfiles("--version")

    assert(stderr == "")
    assert("version" in stdout)
    assert("using chemfiles" in stdout)

    stdout, stderr = cfiles("-V")

    assert(stderr == "")
    assert("version" in stdout)
    assert("using chemfiles" in stdout)
