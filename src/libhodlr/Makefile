.PHONY: clean dependencies dependencies_test
all: libhodlr

OBJS = \
    src/hodlr.o \
    src/HODLR_Matrix.o \
    src/HODLR_Node.o \
    src/HODLR_Tree.o

CXXFLAGS = -Wall -Wno-unknown-pragmas -g -fPIC -I include/ -I src/
CXXFLAGS += -ffast-math

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
	rm -f libhodlr.so src/HODLR_Matrix.o src/HODLR_Node.o src/hodlr.o src/HODLR_Tree.o $(DEPEND_FILE)
