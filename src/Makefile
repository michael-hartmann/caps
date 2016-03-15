.PHONY: clean doc casimir casimir_logdetD casimir_hiT casimir_Fn ld_drude casimir_T0 casimir_tests
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
    tests/test_determinant.o \
    tests/test_epsilon.o \
    tests/test_integration_drude.o \
    tests/test_Lambda.o \
    tests/test_logdet.o \
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


%.o : %.c
	$(CC) $(CFLAGS) $(OPT) -o $@ -c $<

casimir_T0.o:
	$(MPICC) $(CFLAGS) $(OPT) -o $@ -c $<

casimir: casimir.o $(OBJS)
	$(CC) $(CFLAGS) $(OPT) $^ $(LDFLAGS) -o casimir

casimir_logdetD: casimir_logdetD.o $(OBJS)
	$(CC) $(CFLAGS) $(OPT) $^ $(LDFLAGS) -o casimir_logdetD

casimir_hiT: casimir_hiT.o $(OBJS)
	$(CC) $(CFLAGS) $(OPT) $^ $(LDFLAGS) -o casimir_hiT

casimir_Fn: casimir_Fn.o $(OBJS)
	$(CC) $(CFLAGS) $(OPT) $^ $(LDFLAGS) -o casimir_Fn

ld_drude: ld_drude.o $(OBJS)
	$(CC) $(CFLAGS) $(OPT) $^ $(LDFLAGS) -o ld_drude

casimir_tests: $(TESTS) $(OBJS)
	$(CC) $(CFLAGS) $(OPT) $^ $(LDFLAGS) -o casimir_tests
	./casimir_tests

casimir_T0: casimir_T0.o $(OBJS)
	$(MPICC) $(CFLAGS) $(OPT) $^ $(LDFLAGS) -o casimir_T0

install:
	cp casimir casimir_tests casimir_hiT casimir_logdetD casimir_Fn ld_drude casimir_T0 ${TARGET}

clean:
	rm -f casimir casimir_tests casimir_hiT casimir_logdetD casimir_Fn ld_drude casimir_T0 *.o tests/*.o
	rm -rf doc

doc:
	doxygen doxygen.conf
