#!/usr/bin/env bash

# description
# Run all the operation to find code coverag and generate an html report

#set -x #echo all commands
make clean
make

./Tests.exe

gcov --relative-only *.cpp

lcov --capture --directory . --output-file gtest_coverage.info

genhtml gtest_coverage.info --output-directory CODE_COVERAGE

open CODE_COVERAGE/index.html