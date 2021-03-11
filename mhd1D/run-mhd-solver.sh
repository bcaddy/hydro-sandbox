#!/usr/bin/env bash

# description
# Builds, runs, and generates plots for a solver in this directory
#
# Dependencies:
#   it2dl
#
# Inputs:
#   1: -f (optional) The framerate of the animation in fps
#   2: -t (optional) Test the riemann solver using the riemannTester.cpp program

# set -x #echo all commands

EXE="./mhd-solver.exe"

# Get options
while getopts "f:t" opt; do
    case $opt in
        f)  # Set the FPS
            fps=${OPTARG}
            ;;
        t)  # Run the Riemann Tester
            riemannTester="riemannTester"
            EXE="./riemannTester.exe"
            ;;
        \?)
            echo "Invalid option: -${OPTARG}" >&2
            exit 1
            ;;
        :)
            echo "Option -${OPTARG} requires an argument." >&2
            exit 1
            ;;
    esac
done

cd "${0%/*}"  # cd to hydro-sandbox directory

REPO_ROOT=$(git rev-parse --show-toplevel)

cd "${REPO_ROOT}/mhd1D"

echo -e "\nCompiling..."
make ${riemannTester}

echo -e "\nRunning..."
${EXE}

cd "${REPO_ROOT}/visualization"

echo -e "\n\nGenerating Plot..."
./mhd-plotter.py ${fps}

echo -e "\nDownloading Plot..."
"${HOME}/.iterm2/it2dl" output-mhd.mp4
echo -e "Download Complete."

# Get attention
~/.iterm2/it2attention fireworks