.PHONY: clean dependencies dependencies_test
all: libhodlr

OBJS = src/hodlr.o src/HODLR_Matrix.o src/HODLR_Node.o src/HODLR_Tree.o src/HODLR_Tree_NonSPD.o src/HODLR_Tree_SPD.o src/KDTree.o

CXXFLAGS = -std=c++11 -Wall -g -fPIC -DUSE_DOUBLE -I include/ -I src/

# Define the NDEBUG macro when compiling: this disables some run-time checks
# which slightly speeds up your program. As the speedup is negligible, it is
# recommended to keep run-time checks enabled.
# CXXFLAGS += -DNDEBUG

# ignore omp pragmas
CXXFLAGS += -Wno-unknown-pragmas

# The HODLR library works with -ffast-math but the performance gain is
# negligible.
#CXXFLAGS += -ffast-math

# optimize
OPT = -O3

# This file contains header dependencies. dependencies and dependencies_test will update this file.
DEPEND_FILE = .depend

dependencies: $(OBJS:.o=.cpp)
	rm -f $(DEPEND_FILE)
	$(CXX) $(OPT) $(CXXFLAGS) -MM $^ > $(DEPEND_FILE)
dependencies_test: $(TESTS:.o=.c) dependencies
	$(CXX) $(OPT) $(CXXFLAGS) -MM $(TESTS:.o=.c) | sed 's/.*/tests\/&/g' >> $(DEPEND_FILE)
-include $(DEPEND_FILE)

%.o : %.cpp
	$(CXX) $(OPT) $(CXXFLAGS) -o $@ -c $<

libhodlr: $(OBJS) dependencies
	$(CXX) -o libhodlr.so $(OBJS) -shared

clean:
	rm -f libhodlr.so src/*.o $(DEPEND_FILE)
