BUILD_DIR=./build
DEBUG_DIR=$(BUILD_DIR)/Debug
RELEASE_DIR=$(BUILD_DIR)/Release

all: release debug test

debug:
	@mkdir -p $(DEBUG_DIR)
	cd $(DEBUG_DIR); cmake -DCMAKE_BUILD_TYPE=Debug ../..; make

release:
	@mkdir -p $(RELEASE_DIR)
	cd $(RELEASE_DIR); cmake -DCMAKE_BUILD_TYPE=Release ../..; make

test:
	cd $(DEBUG_DIR); ctest
	cd $(RELEASE_DIR); ctest
