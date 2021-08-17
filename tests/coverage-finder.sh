#!/usr/bin/env bash

# description
# Run all the operation to find code coverag and generate an html report

#set -x #echo all commands
make clean
make

./Tests.exe
# --capture = get the data
# --directory = where the program is? can be done more than once
# --output-file = results output file. Should end in ".info"
lcov --capture --directory . --output-file coverage.info

exclude_patterns=('/usr/*'      # Remove everything from /usr/
                  '/Library/*'  # Remove everything from /Library/
                  '*-tests.cpp' # Remove traces of the tests themselves
                  '*-test.cpp') # Remove traces of the tests themselves
# --remove TRACEFILE PATTERN = remove all things associated with PATTERN in TRACEFILE
lcov --remove coverage.info "${exclude_patterns[@]}" --output-file coverage.info  # Remove traces of the tests themselves

lcov --list coverage.info

if [[ $1 == "html" ]]; then
    echo -e "\n\n===== Generate HTML ====================================================="
    genhtml coverage.info --output-directory Code-coverage-html
    open Code-coverage-html/index.html
fi