#.DEFAULT_GOAL := generate #choose the default goal instead of it being the first one
.PHONY: all clean #tells make that these goals are not files but some other thing

# compiler options, debug and performance
cppDebug = g++ -std=c++17 -Wall -Wextra -Wpedantic
cppPerf = g++ -std=c++17 -O3
cppCompiler = ${cppDebug} #the version used

# list of binaries to clean as well
binaries = advection.exe burgers.exe

Helper_SRC = helper.cpp initial-conditions.cpp interface.cpp limiters.cpp
Helper_OBJ = $(Helper_SRC:%.cpp=%.o)

all:
	@echo "This does nothing. See Makefile"

# To switch between compiling the whole program every time and only compiling
# what has changed simply replace NAME.cpp with NAME.o and $(Helper_SRC) with
# $(Helper_OBJ). This will generate a lot of object files though
advection: advection-solver.cpp $(Helper_SRC)
	$(cppCompiler) $? -o $@.exe

burgers: burgers-solver.cpp $(Helper_SRC)
	$(cppCompiler) $? -o $@.exe

%.o: %.cpp
	$(cppCompiler) -c $<

clean:
	@echo "Cleaning up..."
	@rm -f $(binaries)
	@rm -f *.mod
	@rm -f *.o
	@rm -f *.exe
	@rm -rf *.dSYM
	@rm -f a.out
	@echo "Done"
