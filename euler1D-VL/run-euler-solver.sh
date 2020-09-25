#!/usr/bin/env bash

# description
# Builds, runs, and generates plots for a solver in this directory
#
# Dependencies:
#   it2dl
#
# Inputs:
#   1: (optional) The framerate of the animation in fps

# set -x #echo all commands

cd "${0%/*}"  # cd to hydro-sandbox directory

REPO_ROOT=$(git rev-parse --show-toplevel)

cd "${REPO_ROOT}/euler1D-VL"

echo -e "\nCompiling..."
make

echo -e "\nRunning..."
"./euler-solver.exe"

cd "${REPO_ROOT}/visualization"

echo -e "\n\nGenerating Plot..."
./euler-solver-plotter.py ${1}

echo -e "\nDownloading Plot..."
"${HOME}/.iterm2/it2dl" output-euler.mp4
echo -e "Download Complete."
