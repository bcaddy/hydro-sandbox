#!/usr/bin/env bash

# description
# Builds, runs, and generates plots for a solver in this repo
#
# Dependencies:
#   it2dl

# set -x #echo all commands

if [ -z "${1}" ]; then
    echo "No solver chosen"
    exit 1
fi

cd "${0%/*}"  # cd to hydro-sandbox directory

REPO_ROOT=$(git rev-parse --show-toplevel)

cd "${REPO_ROOT}/basic-solvers"

echo -e "\nCompiling..."
make $1

echo -e "\nRunning..."
"./${1}.exe"

cd "${REPO_ROOT}/visualization"

echo -e "\n\nGenerating Plot..."
./basic-solvers-plotting.py

echo -e "\nDownloading Plot..."
"${HOME}/.iterm2/it2dl" output.mp4
echo -e "Download Complete."
