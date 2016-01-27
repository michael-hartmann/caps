.PHONY: clean doc casimir casimir_logdetD casimir_hiT casimir_Fn ld_drude casimir_T0 casimir_tests
all: casimir
#casimir_hiT casimir_logdetD casimir_Fn ld_drude casimir_T0 casimir_tests

# Choose a compiler.
#
# At the moment gcc, clang and icc are supported.
CC = gcc
CXX = g++
#CC = icc
MPICC = mpicc
#MPICC = mpiicc

# use C99 standard
CFLAGS  = -std=c99
# use traditional GNU semantics for inline
CFLAGS += -fgnu89-inline
CFLAGS += -I .
# enable warnings and make them fatal
CFLAGS += -Werror -Wall -Wextra -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wno-unused-parameter

CXXFLAGS = 
LDFLAGS += -lm -pthread -lqd

TARGET = "../bin/"

# optimize for native CPU
OPT = -march=native -flto -O3

OBJS = \
    gausslaguerre.o \
    integration_drude.o \
    integration_perf.o \
    libcasimir.o \
    matrix.o \
    sfunc.o \
    utils.o \
    dd.o

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

# support float128
CFLAGS += -DFLOAT128
LDFLAGS += -lquadmath

QD_INC = qd-2.3.17/include
QD_LIB = qd-2.3.17/src/.libs/

# Disable additional checks and assertions.
#
# This option will cause the program to perform additional checks. If an error
# occures, a error message will be printed and the program will be terminated.
# Without this option, the program will proceed calculating and return
# incorrect results, most likely NANs.
#
# Disabling this option will slightly increase performance.
#
# If you don't know what to say here, comment out the line.
# CFLAGS += -DNCHECKS

CFLAGS += -D CASIMIR_DETALG="\"DEBUG\""

%.o : %.c
	$(CC) $(CFLAGS) $(OPT) -o $@ -c $<

casimir_T0.o:
	$(MPICC) $(CFLAGS) $(OPT) -o $@ -c $<

dd.o :
	$(CXX) $(CXXFLAGS) $(OPT) -I $(QD_INC) -o dd.o -c dd.cpp

casimir: casimir.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) -L $(QD_LIB) $^ $(LDFLAGS) -o casimir

casimir_logdetD: casimir_logdetD.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) -L $(QD_LIB) $^ $(LDFLAGS) -o casimir_logdetD

casimir_hiT: casimir_hiT.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) -L $(QD_LIB) $^ $(LDFLAGS) -o casimir_hiT

casimir_Fn: casimir_Fn.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) -L $(QD_LIB) $^ $(LDFLAGS) -o casimir_Fn

ld_drude: ld_drude.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) -L $(QD_LIB) $^ $(LDFLAGS) -o ld_drude

casimir_tests: $(TESTS) $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) -L $(QD_LIB) $^ $(LDFLAGS) -o casimir_tests
	time ./casimir_tests

casimir_T0: casimir_T0.o $(OBJS)
	$(MPICC) $(CXXFLAGS) $(OPT) -L $(QD_LIB) $^ $(LDFLAGS) -o casimir_T0

install:
	cp casimir casimir_tests casimir_hiT casimir_logdetD casimir_Fn ld_drude casimir_T0 ${TARGET}

clean:
	rm -f casimir casimir_tests casimir_hiT casimir_logdetD casimir_Fn ld_drude casimir_T0 *.o tests/*.o

doc:
	doxygen doxygen.conf
