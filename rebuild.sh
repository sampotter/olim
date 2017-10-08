#!/usr/bin/env bash

function error_exit {
    echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2
    exit 1
}

case `uname` in
	'Darwin')
		NPROC=`sysctl -n hw.ncpu`
		;;
	'Linux')
		NPROC=`nproc`
		;;
	*)
		error_exit "unsupported platform: `uname`"
		;;
esac

for build_type in Debug RelWithDebInfo Release
do
	if [ -d "build/$build_type" ]; then
		cd build/$build_type
		ninja -j $NPROC
		cd -
	fi
done
