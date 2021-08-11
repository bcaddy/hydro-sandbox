#.DEFAULT_GOAL := generate #choose the default goal instead of it being the first one
.PHONY: clean #tells make that these goals are not files but some other thing

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
	cppCompiler = /usr/local/bin/g++-11
else
	cppCompiler = g++
endif

# compiler options, debug and performance
cppDebug = -std=c++17 \
		   -Wall \
		   -Wextra \
		   -Wpedantic \
		   -fasynchronous-unwind-tables \
		   -fexceptions \
		   -D_GLIBCXX_ASSERTIONS \
		   -g
# Add these back later when it's code coverage time
#    -fprofile-arcs \
#    -ftest-coverage \
cppPerf = g++ -std=c++17 -O3
cppFlags = ${cppDebug} #the version used

# Linking in google test. Linking `gtest_main` implicitly includes GTests
# default main function. On some builds/distributions `gtest_main` includes the
# contents of `gtest` as well but that isn't something you can count on so you
# should include/link both
# GTEST_ROOT = /usr/local/Cellar/googletest/1.11.0
GTEST_ROOT = ${GOOGLETEST_ROOT}
# GTEST_ROOT = /ihome/crc/install/googletest/1.11.0
# GTEST_ROOT = OLCF_GOOGLETEST_ROOT
INCLUDES = -lpthread \
		   -L$(GTEST_ROOT)/lib64 \
		   -lgtest \
		   -lgtest_main \
		   -I$(GTEST_ROOT)/include

# Defines
# DEFINES = -DFUNC_1
# DEFINES = -DFUNC_2

# list of other files to clean as well
CLEAN_LIST =

SRC = $(wildcard *.cpp)
OBJ = $(SRC:%.cpp=%.o)

# To switch between compiling the whole program every time and only compiling
# what has changed simply replace NAME.cpp with NAME.o and $(SRC) with
# $(OBJ). This will generate a lot of object files though
Tests.exe: $(OBJ)
	$(cppCompiler) $(cppFlags) $? -o $@ $(INCLUDES) $(DEFINES)


%.o: %.cpp
	$(cppCompiler) $(cppFlags)    -c $< $(INCLUDES) $(DEFINES)

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
