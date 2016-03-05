# cfiles, multi-formats analysis tools for chemistry

`cfiles` is an analysis frontend for the
[Chemfiles](https://github.com/chemfiles/chemfiles/) library, implementing analysis
algorithms for post-processing data from molecular or quantum simulations.

This is still in alpha mode, largely untested, created for dog-feeding purpose and as such
it will surely eat your laundry and do nasty things. Feel free to use it, without any
guarantee about the result.

## Functionalities

* Radial distribution functions;
* Convertion from one file format to another.

## Get it, build it

To build it you can run

```bash
git clone --recursive https://github.com/chemfiles/cfiles
cd cfiles
mkdir build
cd build
cmake ..
make install
```

> Note: on CentOS 7 with GCC 4.8, `cfiles` rely on the `boost_regex` library,
> which can be installed using `yum install boost-regex`. The cmake invocation
> must then be `cmake -DBoost_REGEX_LIBRARY_RELEASE=/usr/lib64 ..`.
