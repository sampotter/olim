#!/usr/bin/env bash

for build_type in Debug RelWithDebInfo Release
do
	if [ -d "build/$build_type" ]; then
		cd build/$build_type
		ninja -j `nproc`
		cd -
	fi
done
