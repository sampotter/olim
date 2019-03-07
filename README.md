# libolim

## Downloading

To get the library, the simplest thing to do is clone this repository:

``` shell
# git clone http://github.com/sampotter/olim
# cd olim
# git submodule update --init
```

The second two lines are unnecessary if you don't want to run the
tests.

## Building

To build the C++ library, run the following from the repository's root
directory:

``` shell
# mkdir -p build/Release
# cd build/Release
# cmake -DCMAKE_BUILD_TYPE=Release ../..
# make
```

To build the Python library, after running the above, run (from the
root directory of the project):

``` shell
# python setup.py build_ext -L build/Release -t build/Release
```

### Tests

To build and run the tests:

``` shell
# mkdir -p build/Debug
# cd build/Debug
# cmake -DCMAKE_BUILD_TYPE=Debug ../..
# make
# ctest
```

The final line runs the tests. If the tests take a long time to run,
you can run the above with `RelWithDebInfo` in place of `Debug`. The
downside of this is that it if you try to step through a test's
execution using a debugger like gdb or lldb, it may be hard to
understand what's going on.

## Installing

To install the Python library, from this repository's root directory,
run:

``` shell
# python setup.py install
```

There isn't currently a good way to install the C++ library
system-wide (this will be fixed at some point). Instead, we recommend
using CMake and `add_subdirectory` to consume this library from
another CMake project.

