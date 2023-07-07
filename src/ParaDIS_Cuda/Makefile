#Makefile 
#define variables
objects= ParallelDataRecovery.o DTW.o imputeValue.o profiler.o ReadData.o Z_norm.o
NVCC= /usr/local/cuda-11.7/bin/nvcc               #cuda c compiler
opt= -use_fast_math -O2 -g -G -Xcompiler -fopenmp          #optimization flag

LIBS=  
execname= ParallelDataRecovery
L1CACHE = -Xptxas -dlcm=cg


#compile
$(execname): $(objects)
	#$(NVCC) $(opt) $(L1CACHE) -o $(execname) $(objects) $(LIBS)
	$(NVCC) $(opt) -o $(execname) $(objects) $(LIBS)

ParallelDataRecovery.o: ParallelDataRecovery.cu
	$(NVCC) $(opt) $(ARCH) -c ParallelDataRecovery.cu
DTW.o: DTW.cpp
	$(NVCC) $(opt) $(ARCH) -c DTW.cpp 
imputeValue.o: imputeValue.cu
	$(NVCC) $(opt) $(ARCH) -c imputeValue.cu 
profiler.o: profiler.cpp
	$(NVCC) $(opt) $(ARCH) -c profiler.cpp 
ReadData.o: ReadData.cpp
	$(NVCC) $(opt) $(ARCH) -c ReadData.cpp 
Z_norm.o: Z_norm.cpp
	$(NVCC) $(opt) $(ARCH) -c Z_norm.cpp 


#clean Makefile
clean:
	rm -rf *.o ./ParallelDataRecovery
#end of Makefile



