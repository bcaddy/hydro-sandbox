#.DEFAULT_GOAL := generate #choose the default goal instead of it being the first one
.PHONY: clean #tells make that these goals are not files but some other thing

# compiler options, debug and performance
cppDebug = g++ \
           -std=c++17 \
		   -Wall \
		   -Wextra \
		   -Wpedantic \
		   -fasynchronous-unwind-tables \
		   -fexceptions \
		   -fprofile-arcs \
		   -ftest-coverage \
		   -D_GLIBCXX_ASSERTIONS \
		   -g
cppPerf = g++ -std=c++17 -O3
cppCompiler = ${cppDebug} #the version used

# Linking in google test. Linking `gtest_main` implicitly includes GTests
# default main function. On some builds/distributions `gtest_main` includes the
# contents of `gtest` as well but that isn't something you can count on so you
# should include/link both
GTEST_ROOT = /usr/local/Cellar/googletest/1.11.0
INCLUDES = -L$(GTEST_ROOT)/lib \
           -lgtest \
		   -lgmock \
		   -lgtest_main \
		   -I$(GTEST_ROOT)/include/gtest \
		   -I$(GTEST_ROOT)/include/gmock

# list of other files to clean as well
CLEAN_LIST =

SRC = $(wildcard *.cpp)
OBJ = $(SRC:%.cpp=%.o)

# To switch between compiling the whole program every time and only compiling
# what has changed simply replace NAME.cpp with NAME.o and $(SRC) with
# $(OBJ). This will generate a lot of object files though
Tests.exe: $(OBJ)
	$(cppCompiler) $? -o $@ $(INCLUDES)

%.o: %.cpp
	$(cppCompiler) -c $<

clean:
	@echo "Cleaning up..."
	@rm -f $(CLEAN_LIST)
	@rm -f *.mod
	@rm -f *.o
	@rm -rf *.gch
	@rm -f *.exe
	@rm -rf *.dSYM
	@rm -f *.gcda
	@rm -f *.gcno
	@rm -rf CODE_COVERAGE
	@rm -f *.info
	@rm -f *.gcov
	@rm -f a.out
	@rm -rf doxygen-html
	@echo "Done"
