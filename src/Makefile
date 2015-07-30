.PHONY: clean doc casimir casimir_tests casimir_hiT casimir_logdetD casimir_Fn ld_drude casimir_psd casimir_T0
all: casimir casimir_hiT casimir_logdetD casimir_tests

SOURCES  = sfunc.c integration_drude.c integration_perf.c libcasimir.c matrix.c utils.c gausslaguerre.c
CFLAGS  += -Werror -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -fno-common
LDFLAGS += -lm -pthread


# Choose a compiler.
#
# At the moment gcc and icc are supported.
CC = gcc
#CC = icc
MPICC = mpicc


# Use extended doubles for calculations.
#
# gcc and icc implement long double as 80-bit extended precision type on x86
# hardware. Extended doubles are less precise than quadrupule precision, but
# are much faster (approximately by a factor of 40).
#
# To use quadruple precision, comment out the line.
#
# If you don't know what to do here, set EXTENDED_DOUBLE to 1.
USE_EXTENDED_DOUBLE = 1


# Use LAPACK to calculate the determinant.
#
# If you define USE_LAPACK, the LAPACK routine dgetrf will be used to perform
# LU decompositions. The use of LAPACK is slightly faster, but LAPACK only
# operates on doubles. This will increase the precision. For very small values
# of the determinant, this will cause wrong results.
#
# Example: (no LAPACK support)
# For T=10, x=0.0075, lmax=800, n=131 and m=64:
#   log(det D) = 0                                      using LAPACK
#   log(det D) = -5.45180220412611e-15          without using LAPACK
#
# If you don't know what to say here, comment out the line.
# USE_LAPACK = 1


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


ifdef USE_LAPACK
    CFLAGS += -DUSE_LAPACK
    LDFLAGS += -llapack -lblas
endif

# use optimization
ifeq ($(CC),gcc)
    OPTIMIZATION = -O3
endif
ifeq ($(CC),icc)
    OPTIMIZATION = -fast
endif

# if we need quadrupole precision
ifdef USE_EXTENDED_DOUBLE
    CFLAGS += -DEXTENDED_DOUBLE
else
    ifeq ($(CC),gcc)
        LDFLAGS += -lquadmath
    endif
    ifeq ($(CC),icc)
        LDFLAGS += -Qoption,c,--extended_float_type
    endif
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

casimir_psd:
	$(CC) $(OPTIMIZATION) $(CFLAGS) $(SOURCES) casimir_psd.c $(LDFLAGS) -o casimir_psd

casimir_T0:
	$(MPICC) $(OPTIMIZATION) $(CFLAGS) $(SOURCES) casimir_T0.c $(LDFLAGS) -o casimir_T0

casimir_tests:
	$(CC) $(OPTIMIZATION) -I . $(CFLAGS) $(SOURCES) tests/unittest.c tests/test*.c $(LDFLAGS) -o casimir_tests
	./casimir_tests


clean:
	rm -f casimir casimir_tests casimir_hiT casimir_logdetD *.o casimir_Fn ld_drude casimir_psd

doc:
	doxygen doxygen.conf
