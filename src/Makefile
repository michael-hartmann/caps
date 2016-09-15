.PHONY: clean doc casimir casimir_logdetD casimir_hiT casimir_Fn ld_drude casimir_T0 casimir_tests dependencies dependencies_test
all: casimir casimir_logdetD casimir_tests

# Choose a compiler.
#
# At the moment gcc, clang and icc are supported.
CC = gcc
#CC = icc
MPICC = mpicc
#MPICC = mpiicc

OBJS = \
    gausslaguerre.o \
    integration_drude.o \
    integration_perf.o \
    libcasimir.o \
    matrix.o \
    sfunc.o \
    utils.o

TESTS = \
    tests/test_bessels.o \
    tests/test_doublefact.o \
    tests/test_fresnel.o \
    tests/test_integration_perf.o \
    tests/test_lbinom.o \
    tests/test_mie.o \
    tests/test_Plm.o \
    tests/test_Xi.o \
    tests/test_epsilon.o \
    tests/test_integration_drude.o \
    tests/test_Lambda.o \
    tests/test_logdet.o \
    tests/test_logdet_HT.o \
    tests/test_mie_drude.o \
    tests/tests.o \
    tests/unittest.o

# use C99 standard
CFLAGS  = -std=c99

# use traditional GNU semantics for inline
CFLAGS += -fgnu89-inline
CFLAGS += -I .

# enable warnings and make them fatal
CFLAGS += -Wall -Werror -Wextra -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wno-unused-parameter

# ldflags
LDFLAGS += -lm -pthread

# support float128
#CFLAGS += -DFLOAT128
#LDFLAGS += -lquadmath

# clang supports various runtime checks. This will slow down the program, but
# it's a good and easy way to find bugs.
#CFLAGS += -fsanitize=undefined,memory

TARGET = "../bin/"

# optimize for native CPU
OPT = -march=native -O3
OPT +=-flto

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
	$(MPICC) $(CFLAGS) $(OPT) -o $@ -c $<

casimir: casimir.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) casimir.o $(OBJS) $(LDFLAGS) -o casimir

casimir_logdetD: casimir_logdetD.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) casimir_logdetD.o $(OBJS) $(LDFLAGS) -o casimir_logdetD

casimir_hiT: casimir_hiT.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) casimir_hiT.o $(OBJS) $(LDFLAGS) -o casimir_hiT

casimir_Fn: casimir_Fn.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) casimir_Fn.o $(OBJS) $(LDFLAGS) -o casimir_Fn

ld_drude: ld_drude.o $(OBJS) dependencies
	$(CC) $(CFLAGS) $(OPT) ld_drude.o $(OBJS) $(LDFLAGS) -o ld_drude

casimir_tests: $(TESTS) $(OBJS) dependencies_test
	$(CC) $(CFLAGS) $(OPT) $(TESTS) $(OBJS) $(LDFLAGS) -o casimir_tests
	./casimir_tests

casimir_T0: casimir_T0.o $(OBJS) dependencies
	$(MPICC) $(CFLAGS) $(OPT) casimir_T0.o $(OBJS) $(LDFLAGS) -o casimir_T0

install:
	cp casimir casimir_tests casimir_hiT casimir_logdetD casimir_Fn ld_drude casimir_T0 ${TARGET}

clean:
	rm -f casimir casimir_tests casimir_hiT casimir_logdetD casimir_Fn ld_drude casimir_T0 *.o tests/*.o $(DEPEND_FILE)
	rm -rf doc

doc:
	doxygen doxygen.conf
