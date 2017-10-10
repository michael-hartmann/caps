.PHONY: clean doc casimir_logdetD casimir casimir_tests casimir_scalar matrix_elems dependencies dependencies_test
all: casimir_logdetD

HODLRDIR            = ../../hodlr_wrapper
HODLRINCLUDEDIR     = ../../hodlr_wrapper/include
CQUADPACKINCLUDEDIR = cquadpack/include

SUPPORT_LAPACK = FALSE

GIT_HEAD=$(shell git rev-parse HEAD)
GIT_BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
MACHINE=$(shell uname -a)

# Choose a compiler.
# At the moment gcc, clang and icc are supported.
ifeq ($(CC),)
    CC = gcc
endif

ifeq ($(CC),icc)
    MPICC = mpiicc
else
    MPICC = mpicc
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
    plm.o \
    logfac.o \
    integration.o \
    libcasimir.o \
    material.o \
    matrix.o \
    misc.o \
    cache.o \
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
    tests/test_mie_drude.o

# use C99 standard
CFLAGS = -std=c99 -fPIC

# use traditional GNU semantics for inline
CFLAGS += -I . -I ${CQUADPACKINCLUDEDIR} -I ${HODLRINCLUDEDIR}
ifdef GIT_HEAD
    CFLAGS += -DGIT_HEAD='"$(GIT_HEAD)"'
endif
ifdef GIT_BRANCH
    CFLAGS += -DGIT_BRANCH='"$(GIT_BRANCH)"'
endif
ifdef GIT_HEAD
    CFLAGS += -DMACHINE='"$(MACHINE)"'
endif

# enable warnings
CFLAGS += -g -Wall -Wextra -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wno-unused-parameter -fstrict-aliasing
ifeq ($(CC),icc)
    CFLAGS += -fp-model precise -mkl=sequential
endif

# profiling
#CFLAGS += -pg

# ldflags
LDFLAGS += -lm -L ${HODLRDIR} -lhodlr

ifeq ($(SUPPORT_LAPACK),TRUE)
    LDFLAGS += -lgfortran -lblas -llapack
    CFLAGS  += -DSUPPORT_LAPACK
endif

# clang supports various runtime checks. This will slow down the program, but
# it's a good and easy way to find bugs.
# https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html
#CFLAGS += -fsanitize=undefined,memory -fno-omit-frame-pointer

# optimize
OPT = -O3 -msse2

# link time optimization
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

library: $(OBJS) dependencies
	$(CC) -o libcasimir.so $(OBJS) -shared

casimir_logdetD: library
	$(CC) $(CFLAGS) $(OPT) casimir_logdetD.c -lcasimir $(LDFLAGS) -o casimir_logdetD

casimir_scalar: library
	$(CC) -fopenmp $(CFLAGS) $(OPT) casimir_scalar.c -lcasimir $(LDFLAGS) -o casimir_scalar

matrix_elems: library
	$(CC) $(CFLAGS) $(OPT) matrix_elems.c -lcasimir $(LDFLAGS) -o matrix_elems

integrate_K: library
	$(CC) $(CFLAGS) $(OPT) integrate_K.c -lcasimir $(LDFLAGS) -o integrate_K

casimir_tests: $(TESTS) library
	$(CC) $(CFLAGS) $(OPT) $(TESTS) -lcasimir $(LDFLAGS) -o casimir_tests
	./casimir_tests

casimir: library
	$(MPICC) $(CFLAGS) $(OPT) casimir.c -lcasimir $(LDFLAGS) -o casimir

clean:
	rm -f casimir_tests casimir_logdetD casimir_scalar casimir integrate_K matrix_elems *.o *.so cquadpack/src/*.o tests/*.o $(DEPEND_FILE)
	rm -rf doc

doc:
	doxygen doxygen.conf
