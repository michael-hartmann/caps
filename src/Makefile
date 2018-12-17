.PHONY: clean doc casimir_logdetD casimir casimir_tests dependencies dependencies_test
all: casimir casimir_logdetD

# get head and branch from git, and machine name
GIT := $(shell git)
ifdef GIT
    GIT_HEAD=$(shell git rev-parse HEAD)
    GIT_BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
endif

MACHINE=$(shell uname -a)

ifeq ($(CC),icc)
    MPICC ?= mpiicc
else
    MPICC ?= mpicc
endif

OBJS = \
    cquadpack/src/dqagi.o \
    cquadpack/src/dqk15.o \
    cquadpack/src/dqk61.o \
    cquadpack/src/dqk51.o \
    cquadpack/src/dqk41.o \
    cquadpack/src/dqk31.o \
    cquadpack/src/dqage.o \
    cquadpack/src/dqags.o \
    cquadpack/src/dqext.o \
    cquadpack/src/dqk15i.o \
    cquadpack/src/dqk21.o \
    cquadpack/src/dqsort.o \
    bessel.o \
    fcqs.o \
    plm.o \
    logfac.o \
    integration.o \
    libcasimir.o \
    material.o \
    matrix.o \
    misc.o \
    cache.o \
    psd.o \
    utils.o

TESTS = \
    tests/test_lfac.o \
    tests/test_logi.o \
    tests/test_Lambda.o \
    tests/tests.o \
    tests/unittest.o \
    tests/test_bessels.o \
    tests/test_fresnel.o \
    tests/test_mie.o \
    tests/test_Plm.o \
    tests/test_logdetD.o \
    tests/test_mie_drude.o

# use C99 standard
CFLAGS = -std=c99 -fPIC
CFLAGS += -I include/ -I cquadpack/include -I libhodlr/include/

ifdef GIT_HEAD
    CFLAGS += -DGIT_HEAD='"$(GIT_HEAD)"'
endif
ifdef GIT_BRANCH
    CFLAGS += -DGIT_BRANCH='"$(GIT_BRANCH)"'
endif
ifdef MACHINE
    CFLAGS += -DMACHINE='"$(MACHINE)"'
endif

# enable warnings
CFLAGS += -g -Wall -Wextra -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wno-unused-parameter -fstrict-aliasing

# By default, icc violates strict IEEE floating point behaviour - similar to
# GCC's --fast-math option. The code, however, relies on strict IEEE floating
# point behaviour. The option "-fp-model precise" sets the correct floating
# point behaviour.
# Also, when using the math kernel library (MKL), use the non-threaded library.
ifeq ($(CC),icc)
    CFLAGS += -fp-model precise -mkl=sequential
endif

# Activate for profiling information using gprof.
#CFLAGS += -pg

CXXFLAGS = -std=c++11 -Wall -Wextra -I include/

# ldflags
LDFLAGS += -lm -L . -L libhodlr/ -lhodlr
LDFLAGS += -lgfortran -lblas -llapack

# clang supports various runtime checks. This will significantly slow down the
# program, but the options are useful for debugging
# https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html
# CFLAGS   += -fsanitize=undefined -fno-omit-frame-pointer
# LDFLAGS  += -lubsan

# optimize
OPT = -O3

# If you are running your code on the same machine where you compile the code,
# or the target machine supports the same instruction set, this option may
# increase performance. When there are problems, the program will crash with a
# message similar to "illegal instruction code".
# On my Intel Core i7 machine, this option gives an performance boost of ~5%.
#OPT += -march=native

# Link time optimization. This might further improve performance. On my
# machine, it has next to no impact on performance.
#OPT += -flto

# This file contains header dependencies. dependencies and dependencies_test will update this file.
DEPEND_FILE = .depend

dependencies: $(OBJS:.o=.c)
	rm -f $(DEPEND_FILE)
	$(CC) $(CFLAGS) -MM $^ > $(DEPEND_FILE)
dependencies_test: $(TESTS:.o=.c) dependencies
	$(CC) $(CFLAGS) -MM $(TESTS:.o=.c) | sed 's/.*/tests\/&/g' >> $(DEPEND_FILE)
-include $(DEPEND_FILE)

%.o : %.c
	$(CC) $(CFLAGS) $(OPT) -o $@ -c $<

libcasimir.so: libhodlr.so $(OBJS) dependencies
	$(CC) -o libcasimir.so $(OBJS) -shared

libhodlr.so:
	cd libhodlr; make

casimir_logdetD: libcasimir.so
	$(CC) $(CFLAGS) $(OPT) casimir_logdetD.c -lcasimir $(LDFLAGS) -o casimir_logdetD

casimir_tests: $(TESTS) libcasimir.so
	$(CC) $(CFLAGS) $(OPT) $(TESTS) -lcasimir $(LDFLAGS) -o casimir_tests
	./casimir_tests

casimir: libcasimir.so
	$(MPICC) $(CFLAGS) $(OPT) casimir.c -lcasimir $(LDFLAGS) -o casimir

cylinder: libcasimir.so cylinder.o argparse.o
	$(CXX) $(OPT) argparse.o cylinder.o -o cylinder -lcasimir $(LDFLAGS)


clean:
	rm -f cylinder casimir_tests casimir_logdetD casimir *.o *.so cquadpack/src/*.o tests/*.o $(DEPEND_FILE)
	rm -rf doc
	cd libhodlr; make clean

doc:
	doxygen doxygen.conf
