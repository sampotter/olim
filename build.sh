#!/usr/bin/env bash

if [ `uname -s` == "Linux" ]; then
	rm -rf build
	for build_type in Debug RelWithDebInfo Release
	do
		mkdir -p build/$build_type
		cd build/$build_type
		CXX=clang++ cmake -DCMAKE_CXX_FLAGS="-stdlib=libc++" \
		   -DCMAKE_BUILD_TYPE=$build_type -GNinja ../..
		ninja -j `nproc`
		cd -
	done
elif [ `uname -s` == "Darwin" ]; then
	rm -rf build
	for build_type in Debug RelWithDebInfo Release
	do
		mkdir -p build/$build_type
		cd build/$build_type
		cmake -DCMAKE_BUILD_TYPE=$build_type -GNinja ../..
		ninja -j `sysctl -n hw.ncpu`
		cd -
	done
else
	echo "Unsupported platform" 1>&2
	exit 1
fi

# Local Variables:
# indent-tabs-mode: nil
# End:
