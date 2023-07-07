CFLAGS += -O3 -w -fopenmp

all: clean ParallelDataRecovery

ParallelDataRecovery:
	g++ $(CFLAGS) -o ParallelDataRecovery ParallelDataRecovery.cpp DTW.cpp imputeValue.cpp profiler.cpp lowerBounds.cpp ReadData.cpp Z_norm.cpp

clean:
	rm -rf *.o ./ParallelDataRecovery

