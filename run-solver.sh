#!/usr/bin/env bash

# description
# Builds, runs, and generates plots for a solver in this repo
#
# Dependencies:
#   it2dl

set -x #echo all commands

if [ -z "${1}"]; then
    echo "No solver chosen"
    exit 1
fi

cd "${0%/*}"  # cd to hydro-sandbox directory

REPO_ROOT=$(git rev-parse --show-toplevel) 

cd "${REPO_ROOT}/solver"

make $1

"./${1}.exe"

cd "${REPO_ROOT}/visualization"

./plotter.py

"${HOME}/.iterm2/it2dl" output.mp4