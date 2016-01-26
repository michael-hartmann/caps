.PHONY: clean doc casimir casimir_tests casimir_hiT casimir_logdetD casimir_Fn ld_drude casimir_T0
all: casimir casimir_hiT casimir_logdetD casimir_Fn ld_drude casimir_T0 casimir_tests

SOURCES  = sfunc.c integration_drude.c integration_perf.c libcasimir.c matrix.c utils.c gausslaguerre.c
# use C99 standard
CFLAGS  += -std=c99
# optimize for native CPU
CFLAGS  += -march=native
# use traditional GNU semantics for inline
CFLAGS  += -fgnu89-inline
# enable warnings and make them fatal
CFLAGS  += -Werror -Wall -Wextra -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wno-unused-parameter

LDFLAGS += -lm -pthread -lquadmath

TARGET = "../bin/"

# Choose a compiler.
#
# At the moment gcc, clang and icc are supported.
CC = gcc
#CC = icc
MPICC = mpicc
#MPICC = mpiicc


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
# CFLAGS += -DNDEBUG

# use optimization
ifeq ($(CC),gcc)
    OPTIMIZATION = -O3
endif
ifeq ($(CC),icc)
    OPTIMIZATION = -fast
endif


casimir:
	$(CC) $(OPTIMIZATION) $(CFLAGS) $(SOURCES) casimir.c $(LDFLAGS) -o casimir

casimir_logdetD:
	$(CC) $(OPTIMIZATION) $(CFLAGS) $(SOURCES) casimir_logdetD.c $(LDFLAGS) -o casimir_logdetD

casimir_hiT:
	$(CC) $(OPTIMIZATION) $(CFLAGS) $(SOURCES) casimir_hiT.c $(LDFLAGS) -o casimir_hiT

casimir_Fn:
	$(CC) $(OPTIMIZATION) $(CFLAGS) $(SOURCES) casimir_Fn.c $(LDFLAGS) -o casimir_Fn

ld_drude:
	$(CC) $(OPTIMIZATION) $(CFLAGS) $(SOURCES) ld_drude.c $(LDFLAGS) -o ld_drude

casimir_T0:
	$(MPICC) $(OPTIMIZATION) $(CFLAGS) $(SOURCES) casimir_T0.c $(LDFLAGS) -o casimir_T0

casimir_tests:
	$(CC) $(OPTIMIZATION) -I . $(CFLAGS) $(SOURCES) tests/unittest.c tests/test*.c $(LDFLAGS) -o casimir_tests
	./casimir_tests

install:
	cp casimir casimir_tests casimir_hiT casimir_logdetD casimir_Fn ld_drude casimir_T0 ${TARGET}

clean:
	rm -f casimir casimir_tests casimir_hiT casimir_logdetD *.o casimir_Fn ld_drude casimir_T0

doc:
	doxygen doxygen.conf
