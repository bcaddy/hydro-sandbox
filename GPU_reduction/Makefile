#.DEFAULT_GOAL := generate #choose the default goal instead of it being the first one
.PHONY: clean #tells make that these goals are not files but some other thing

# compiler options, debug and performance
cudaCompiler = nvcc -g -O3 -std=c++11 -arch sm_70 -fmad=false -I/usr/local/cuda-11.4/include -DCUDA_BUILD
hipCompiler = hipcc --offload-arch=gfx90a -g -O3 -std=c++11 -Wall -ferror-limit=1 -DHIP_BUILD

chollaFlags = -DCUDA -DMPI_CHOLLA -DPRECISION=2 -DPPMP -DHLLC -DVL -DDENSITY_FLOOR -DTEMPERATURE_FLOOR -DDE -DCPU_TIME -DOUTPUT -DHDF5 -DMPI_GPU -DPARALLEL_OMP -DN_OMP_THREADS=8 -DO_HIP

# list of binaries to clean as well
binaries = gpu_reductions.exe

# To switch between compiling the whole program every time and only compiling
# what has changed simply replace NAME.cpp with NAME.o and $(Helper_SRC) with
# $(Helper_OBJ). This will generate a lot of object files though
cuda: gpu_reductions.cu
	$(cudaCompiler) $? -o gpu_reductions.exe

hip: gpu_reductions.cu
	$(hipCompiler) $? -o gpu_reductions.exe

# %.o: %.cu
# 	$(cudaCompiler) -c $<

clean:
	@echo "Cleaning up..."
	@rm -f $(binaries)
	@rm -f *.mod
	@rm -f *.o
	@rm -f *.gcno
	@rm -rf *.gch
	@rm -f *.exe
	@rm -rf *.dSYM
	@rm -f a.out
	@rm -rf doxygen-html
	@echo "Done"
