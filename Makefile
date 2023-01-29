CXXFLAGS=-g -Wall -std=c++11 -O2 -I/usr/local/include -L/usr/local/lib -lhdf5
CC=g++
CC_M=mpic++
CXX_MFLAGS=-fopenmp

# Note that those variables will be (usually) overwritten by one of the python scripts
EPS?=5000.0
MIN_PTS?=5
FILE?=input/50Kpts.ds
INDEX?=0

DEPTH?=1
STRATEGY?=0
OUT_DIR?=inputSplit
NO_MACHINES?=2
ORIG_FILE?=50Kpts.ds
VERBOSE?=0
OPTIMIZE_COMM?=0

# Generating some .o files ; this is useful to avoid recompiling everything (apparently)
UtilGen.o: UtilGen.cpp
	$(CC) $(CXXFLAGS) -c UtilGen.cpp

inputGenerator.o: inputGenerator.cpp
	$(CC) $(CXXFLAGS) $(CXX_MFLAGS) -c inputGenerator.cpp

dbScanDSU.o: dbScanDSU.cpp
	$(CC) $(CXXFLAGS) -c dbScanDSU.cpp

bruteDBScan.o: bruteDBScan.cpp
	$(CC) $(CXXFLAGS) -c bruteDBScan.cpp

dbScanMultiNode.o: dbScanMultiNode.cpp
	$(CC_M) $(CXXFLAGS) $(CXX_MFLAGS) -c dbScanMultiNode.cpp

# Compile targets, yielding the binary executables. Something tells me this is not how a Makefile
# should be used, but meh...
compile: bruteDBScan.o
	$(CC) $(CXXFLAGS) -o bruteDBScan bruteDBScan.cpp

compile_dsu: dbScanDSU.o
	$(CC) $(CXXFLAGS) dbScanDSU.o -o dbScanDSU

compile_gen: inputGenerator.o UtilGen.o
	$(CC) $(CXXFLAGS) $(CXX_MFLAGS) inputGenerator.o UtilGen.o -o inputGenerator

compile_multinode: dbScanMultiNode.o
	$(CC_M) $(CXXFLAGS) $(CXX_MFLAGS) dbScanMultiNode.o -o dbScanMulti


# The executable targets.
exec:
	./bruteDBScan -EPS $(EPS) -MIN_PTS $(MIN_PTS) -FILE $(FILE) -INDEX $(INDEX) -VERBOSE $(VERBOSE)

exec_dsu:
	./dbScanDSU -EPS $(EPS) -MIN_PTS $(MIN_PTS) -FILE $(FILE) -INDEX $(INDEX) -VERBOSE $(VERBOSE)

exec_gen:
	./inputGenerator -EPS $(EPS) -DEPTH $(DEPTH) -FILE $(FILE) -OUT_DIR $(OUT_DIR) \
	 -STRATEGY $(STRATEGY) -VERBOSE $(VERBOSE)

exec_multinode:
	mpiexec -machinefile hosts -np ${NO_MACHINES} ./dbScanMulti -INPUT_DIR $(OUT_DIR) \
	-ORIG_FILE $(ORIG_FILE) -EPS $(EPS) -MIN_PTS $(MIN_PTS) -VERBOSE $(VERBOSE) -OPTIMIZE_COMM $(OPTIMIZE_COMM)

clean:
	rm -rf bruteDBScan inputGenerator dbScanDSU dbScanMulti *.o
