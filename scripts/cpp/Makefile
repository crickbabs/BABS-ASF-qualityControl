#
# nourdine.bah@crick.ac.uk
#

# -static-libstdc++
# because conda cannot find the GLIBCXX_3.4.21 in the shared library
#
CXX = g++
CFLAG=-I /camp/stp/babs/working/bahn/code/cpp/seqan/include
LIB=-lpthread -lrt -lz -static-libstdc++
OPT=-std=c++14 -O3 -W -Wall -pedantic -fopenmp
DEF=-DNDEBUG -DSEQAN_HAS_ZLIB=1 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0

TRANS=transcript_directionality
QC=qc_directionality

.PHONY: .FORCE
.FORCE:

all: $(TRANS) $(QC)

$(TRANS): .FORCE
	$(CXX) $(CFLAG) $(LIB) $(OPT) $(DEF) -o $(TRANS) $(TRANS).cpp

$(QC): .FORCE
	$(CXX) $(CFLAG) $(LIB) $(OPT) $(DEF) -o $(QC) $(QC).cpp

