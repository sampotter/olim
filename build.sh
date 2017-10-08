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

function build_target {
    CMAKE="cmake -DCMAKE_BUILD_TYPE=${1} -GNinja ../.."
    case `uname` in
        'Linux')
            CMAKE='CXX=clang++ $CMAKE_CMD -DCMAKE_CXX_FLAGS="-stdlib=libc++"'
            ;;
        *)
            ;;
    esac
    mkdir -p build/${1}
    cd build/${1}
    $CMAKE
    ninja -j $NPROC
    cd -
}

for BUILD_TYPE in Debug RelWithDebInfo Release
do
    build_target $BUILD_TYPE
done

# Local Variables:
# indent-tabs-mode: nil
# End:
