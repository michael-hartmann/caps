.PHONY: clean doc casimir casimir_logdetD casimir_hiT casimir_T0 casimir_tests dependencies dependencies_test
all: casimir casimir_logdetD

# Choose a compiler.
# At the moment gcc, clang and icc are supported.
CC = gcc # icc, clang
MPICC = mpicc # mpiicc

OBJS = \
    cquadpack/src/dqagi.c \
    cquadpack/src/dqk15.c \
    cquadpack/src/dqk61.c \
    cquadpack/src/dqk51.c \
    cquadpack/src/dqk41.c \
    cquadpack/src/dqk31.c \
    cquadpack/src/dqage.c \
    cquadpack/src/dqags.c \
    cquadpack/src/dqext.c \
    cquadpack/src/dqk15i.c \
    cquadpack/src/dqk21.c \
    cquadpack/src/dqsort.c \
    hash-table.c \
    lookup.c \
    integration.o \
    libcasimir.o \
    matrix.o \
    sfunc.o \
    gausslaguerre.c \
    utils.o

TESTS = \
    tests/test_Lambda.o \
    tests/tests.o \
    tests/unittest.o \
    tests/test_bessels.o \
    tests/test_doublefact.o \
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
CFLAGS += -I . -I cquadpack/include

# enable warnings
CFLAGS += -g -Wall -Wextra -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wno-unused-parameter -fstrict-aliasing
# profiling
#CFLAGS += -pg

# ldflags
LDFLAGS += -lm -pthread -llapack -lblas -lgfortran

# clang supports various runtime checks. This will slow down the program, but
# it's a good and easy way to find bugs.
#CFLAGS += -fsanitize=undefined,memory

TARGET = "../bin/"

# optimize for native CPU
OPT = -march=native -O3 -mmmx -msse -msse2 -msse3 -msse4

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

casimir_T0.o:
	$(MPICC) $(CFLAGS) $(OPT) casimir_T0.c -o $@ -c $<

casimir: casimir.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) casimir.o $(OBJS) $(LDFLAGS) -o casimir

casimir_logdetD: casimir_logdetD.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) casimir_logdetD.o $(OBJS) $(LDFLAGS) -o casimir_logdetD

integrate_K: integrate_K.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) integrate_K.o $(OBJS) $(LDFLAGS) -o integrate_K

casimir_hiT: casimir_hiT.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) casimir_hiT.o $(OBJS) $(LDFLAGS) -o casimir_hiT

casimir_tests: $(TESTS) $(OBJS) dependencies_test
	$(CC) $(CFLAGS) $(OPT) $(TESTS) $(OBJS) $(LDFLAGS) -o casimir_tests
	./casimir_tests

casimir_T0: casimir_T0.o $(OBJS) dependencies
	$(MPICC) $(CFLAGS) $(OPT) casimir_T0.o $(OBJS) $(LDFLAGS) -o casimir_T0

install:
	cp casimir casimir_tests casimir_hiT casimir_logdetD casimir_T0 ${TARGET}

clean:
	rm -f casimir casimir_tests casimir_hiT casimir_logdetD casimir_T0 *.o tests/*.o $(DEPEND_FILE)
	rm -rf doc

doc:
	doxygen doxygen.conf
