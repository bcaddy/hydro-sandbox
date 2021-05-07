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
		   -D_GLIBCXX_ASSERTIONS \
		   -g \
		   -O0
cppPerf = g++ -std=c++17 -O3
cppCompiler = ${cppDebug} #the version used

# list of binaries to clean as well
binaries = euler-solver.exe test-Riemann.exe

Helper_SRC = Grid1D.cpp MhdSimulation1D.cpp HlldRiemannSolver.cpp
Helper_OBJ = $(Helper_SRC:%.cpp=%.o)

# To switch between compiling the whole program every time and only compiling
# what has changed simply replace NAME.cpp with NAME.o and $(Helper_SRC) with
# $(Helper_OBJ). This will generate a lot of object files though
mhd-solver: mhd-main.cpp $(Helper_SRC)
	$(cppCompiler) $? -o $@.exe


test-Helper_SRC = HlldRiemannSolver.cpp Grid1D.cpp
riemannTester: riemannTester.cpp $(test-Helper_SRC)
	$(cppCompiler) $? -o $@.exe

%.o: %.cpp
	$(cppCompiler) -c $<

clean:
	@echo "Cleaning up..."
	@rm -f $(binaries)
	@rm -f *.mod
	@rm -f *.o
	@rm -rf *.gch
	@rm -f *.exe
	@rm -rf *.dSYM
	@rm -f a.out
	@rm -rf doxygen-html
	@echo "Done"
