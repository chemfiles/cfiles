# cfiles, multi-formats analysis tools for chemistry

[![Build Status](https://travis-ci.org/chemfiles/cfiles.svg?branch=master)](https://travis-ci.org/chemfiles/cfiles)
[![codecov](https://codecov.io/gh/chemfiles/cfiles/branch/master/graph/badge.svg)](https://codecov.io/gh/chemfiles/cfiles)

`cfiles` is an analysis frontend for the
[Chemfiles](https://github.com/chemfiles/chemfiles/) library, implementing
analysis algorithms for post-processing data from molecular or quantum
simulations.

This is still in alpha mode, largely untested, created for dog-feeding purpose
and as such it will surely eat your laundry and do nasty things. Feel free to
use it, without any guarantee about the result.

## Functionalities

* Radial distribution functions;
* Angles and dihedral angles distributions;
* Convertion from one file format to another.

## Get it, build it

To build it you can run

```bash
git clone --recursive https://github.com/chemfiles/cfiles
cd cfiles
mkdir build
cd build
cmake ..
make
```

You can then copy the `cfiles` executable to some place in your path

**Note:** When using the C++ standard library from GCC 4.8, `cfiles` rely on
the `boost_regex` library, which can be installed using your favorite package
manager.

## Contributing

You are welcome to contribute to cfiles, whatever your skill level. You can help
with new analysis algorithms, improving existing ones, adding or improving
documentation, *etc.* If you plan adding new code, please open an [issue] about
it for discussion.

Don't forget to run the tests when changing the code; and to add new tests when
adding new algorithms. The tests are Python scripts in the `tests` directory,
checking the output of the algorithms.

```bash
cd cfiles/build
ctest
```

Here is a short check list to contribute to cfiles. If there is anything you
don't understand, or if you have any question, please ask! You can reach me on
github issues, by email, or in the [gitter] chat.

- [ ] Fork cfiles;
- [ ] Create a local branch;
- [ ] Add code / correct typos / ...;
    - [ ] Add new tests with your code;
    - [ ] Add documentation for your code;
- [ ] Check that the tests still pass;
- [ ] Push to Github;
- [ ] Create a Pull-Request;
- [ ] Discuss your changes with the reviewers;
- [ ] Have your code merged;
- [ ] Celebrate! :tada: :cake: :tada:


[gitter]: https://gitter.im/chemfiles/chemfiles
[issue]: https://github.com/chemfiles/cfiles/issues
