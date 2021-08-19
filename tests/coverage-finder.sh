#!/usr/bin/env bash

# description
# Run all the operation to find code coverag and generate an html report

# Get options
while getopts "hb" opt; do
    case $opt in
        h)  # Choose whether or not to generate html results
            HTML=true
            ;;
        b)  # Allow building and running
            BUILD=true
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

# Define the root of the repo
repo_root=$(git rev-parse --show-toplevel)

if [[ $BUILD == true ]]; then
    #set -x #echo all commands
    make clean
    make

    # make euler1D-VL and mhd1D files
    cd ${repo_root}/euler1D-VL
    rm -rf *.o *.gcno
    g++ -std=c++17 -Wall -Wextra -Wpedantic -fasynchronous-unwind-tables -fexceptions -D_GLIBCXX_ASSERTIONS -g --coverage -lgtest -lgtest_main -I/usr/local/Cellar/googletest/1.11.0/include -L/usr/local/Cellar/googletest/1.11.0/lib -lpthread -c *.cpp
    make

    cd ${repo_root}/mhd1D
    rm -rf *.o *.gcno
    g++ -std=c++17 -Wall -Wextra -Wpedantic -fasynchronous-unwind-tables -fexceptions -D_GLIBCXX_ASSERTIONS -g --coverage -lgtest -lgtest_main -I/usr/local/Cellar/googletest/1.11.0/include -L/usr/local/Cellar/googletest/1.11.0/lib -lpthread -c *.cpp

    cd ${repo_root}/tests

    ./Tests.exe
fi

# generate initial report
capture_directories=(--directory ${repo_root}/tests
                     --directory ${repo_root}/euler1D-VL
                     --directory ${repo_root}/mhd1D)
lcov --capture --initial ${capture_directories[@]} --output-file coverage_base.info

# Generate test results
# --capture = get the data
# --directory = where the program is? can be done more than once
# --output-file = results output file. Should end in ".info"
lcov --capture ${capture_directories[@]} --output-file coverage_test.info

# Combine base and test results
lcov --add-tracefile coverage_base.info --add-tracefile coverage_test.info --output-file coverage_all.info

# Extract data from only the files within Repo_root. This should exclude any
# system or external libraries
lcov --extract coverage_all.info "${repo_root}/*" --output-file coverage_all.info

exclude_patterns=('*-tests.cpp' # Remove traces of the tests themselves
                  '*-test.cpp') # Remove traces of the tests themselves
# --remove TRACEFILE PATTERN = remove all things associated with PATTERN in TRACEFILE
lcov --remove coverage_all.info "${exclude_patterns[@]}" --output-file coverage_all.info




lcov --list coverage_all.info

if [[ $HTML == true ]]; then
    echo -e "\n\n===== Generate HTML ====================================================="
    genhtml coverage_all.info --output-directory Code-coverage-html
    open Code-coverage-html/index.html
fi
