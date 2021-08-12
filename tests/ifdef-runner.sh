#!/usr/bin/env bash

# description

#set -x #echo all commands

export IFDEF_DEFINES=-DFUNC_1

make clean
make

# Run everything except the testfunc2 test
./Tests.exe --gtest_filter=-ifdef.*Two

export IFDEF_DEFINES=-DFUNC_2

make clean
make

# Run only the testfunc2 test
./Tests.exe --gtest_filter=ifdef.*Two