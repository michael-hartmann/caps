.PHONY: clean doc casimir_logdetD casimir casimir_tests matrix_elems dependencies dependencies_test
all: casimir_logdetD

HODLRDIR            = ../../hodlr_wrapper
HODLRINCLUDEDIR     = ../../hodlr_wrapper/include
CQUADPACKINCLUDEDIR = cquadpack/include

SUPPORT_LAPACK = FALSE

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
    hash-table.o \
    logfac.o \
    integration.o \
    libcasimir.o \
    material.o \
    matrix.o \
    misc.o \
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
#   tests/test_integration_perf.o \
#   tests/test_logdet_HT.o \
#   tests/test_logdet.o \
# tests/test_integration_drude.o

# use C99 standard
CFLAGS = -std=c99

# use traditional GNU semantics for inline
#CFLAGS += -fgnu89-inline
CFLAGS += -I . -I ${CQUADPACKINCLUDEDIR} -I ${HODLRINCLUDEDIR}

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

casimir.o:
	$(MPICC) $(CFLAGS) $(OPT) casimir.c -o $@ -c $<

casimir_logdetD: casimir_logdetD.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) casimir_logdetD.o $(OBJS) $(LDFLAGS) -o casimir_logdetD

test: test.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) test.o $(OBJS) $(LDFLAGS) -o test

matrix_elems: matrix_elems.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) matrix_elems.o $(OBJS) $(LDFLAGS) -o matrix_elems

integrate_K: integrate_K.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) integrate_K.o $(OBJS) $(LDFLAGS) -o integrate_K

casimir_tests: $(TESTS) $(OBJS) dependencies_test
	$(CC) $(CFLAGS) $(OPT) $(TESTS) $(OBJS) $(LDFLAGS) -o casimir_tests
	./casimir_tests

casimir: casimir.o $(OBJS) dependencies
	$(MPICC) $(CFLAGS) $(OPT) casimir.o $(OBJS) $(LDFLAGS) -o casimir

clean:
	rm -f casimir_tests casimir_hiT casimir_logdetD casimir matrix_elems *.o tests/*.o $(DEPEND_FILE)
	rm -rf doc

doc:
	doxygen doxygen.conf
