#.DEFAULT_GOAL := generate #choose the default goal instead of it being the first one
.PHONY: advection.exe all clean #tells make that these goals are not files but some other thing

# compiler options, debug and performance
cppDebug = g++ -std=c++17 -Wall -Wextra -Wpedantic
cppPerf = g++ -std=c++17 -O3
cppCompiler = ${cppDebug} #the version used

# list of binaries to clean as well
binaries = tally.exe

Advect_SRC = 2nd-Order-linear-advection.cpp helper.cpp
Advect_OJB = $(SRC:.cpp = .o)


all:
	@echo "This does nothing. See Makefile"

advection: $(Advect_OBJ)
	$(cppCompiler) $(Advect_SRC) -o advection.exe $(Advect_OBJ)

clean:
	@echo "Cleaning up..."
	@rm -f $(binaries)
	@rm -f *.mod
	@rm -f *o
	@rm -f *.exe
	@rm -rf *.dSYM
	@rm -f a.out
	@echo "Done"
