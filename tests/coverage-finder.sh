#!/usr/bin/env bash

# description
# Run all the operation to find code coverag and generate an html report

#set -x #echo all commands
make clean
make


# make euler1D-VL and mhd1D files
repo_root=$(git rev-parse --show-toplevel)
echo $repo_root
cd ${repo_root}/euler1D-VL
rm -rf *.o *.gcno
g++ -std=c++17 -Wall -Wextra -Wpedantic -fasynchronous-unwind-tables -fexceptions -D_GLIBCXX_ASSERTIONS -g --coverage -lgtest -lgtest_main -I/usr/local/Cellar/googletest/1.11.0/include -L/usr/local/Cellar/googletest/1.11.0/lib -lpthread -c *.cpp

cd ${repo_root}/mhd1D
rm -rf *.o *.gcno
g++ -std=c++17 -Wall -Wextra -Wpedantic -fasynchronous-unwind-tables -fexceptions -D_GLIBCXX_ASSERTIONS -g --coverage -lgtest -lgtest_main -I/usr/local/Cellar/googletest/1.11.0/include -L/usr/local/Cellar/googletest/1.11.0/lib -lpthread -c *.cpp

cd ${repo_root}/tests

# generate initial report
capture_directories=(--directory ${repo_root}/tests
                     --directory ${repo_root}/euler1D-VL
                     --directory ${repo_root}/mhd1D)
lcov --capture --initial ${capture_directories[@]} --output-file coverage_base.info

./Tests.exe

# Generate test results
# --capture = get the data
# --directory = where the program is? can be done more than once
# --output-file = results output file. Should end in ".info"
lcov --capture ${capture_directories[@]} --output-file coverage_test.info

# Combine base and test results
lcov --add-tracefile coverage_base.info --add-tracefile coverage_test.info --output-file coverage_all.info

exclude_patterns=('/usr/*'      # Remove everything from /usr/
                  '/Library/*'  # Remove everything from /Library/
                  '*-tests.cpp' # Remove traces of the tests themselves
                  '*-test.cpp') # Remove traces of the tests themselves
# --remove TRACEFILE PATTERN = remove all things associated with PATTERN in TRACEFILE
lcov --remove coverage_all.info "${exclude_patterns[@]}" --output-file coverage_all.info  # Remove traces of the tests themselves

lcov --list coverage_all.info

if [[ $1 == "html" ]]; then
    echo -e "\n\n===== Generate HTML ====================================================="
    genhtml coverage_all.info --output-directory Code-coverage-html
    open Code-coverage-html/index.html
fi