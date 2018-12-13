.PHONY: clean doc casimir_logdetD casimir casimir_tests dependencies dependencies_test
all: casimir casimir_logdetD

# get head and branch from git, and machine name
GIT_HEAD=$(shell git rev-parse HEAD)
GIT_BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
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
    argparse.o \
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
ifdef GIT_HEAD
    CFLAGS += -DMACHINE='"$(MACHINE)"'
endif

# enable warnings
CFLAGS += -g -Wall -Wextra -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wno-unused-parameter -fstrict-aliasing

# this is necessary to get correct results when using icc
ifeq ($(CC),icc)
    CFLAGS += -fp-model precise -mkl=sequential
endif

# profiling
#CFLAGS += -pg

CXXFLAGS = -Wall -Wextra -I include/

# ldflags
LDFLAGS += -lm -L . -L libhodlr/ -lhodlr
LDFLAGS += -lgfortran -lblas -llapack

# clang supports various runtime checks. This will slow down the program, but
# it's a good and easy way to find bugs.
# https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html
ifdef SANITIZE
    CFLAGS   += -fsanitize=undefined -fno-omit-frame-pointer
    LDFLAGS  += -lubsan
endif

# optimize
OPT = -O3

# link time optimization
ifdef LTO
    OPT += -flto
endif

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

cylinder: libcasimir.so
	$(CXX) $(CXXFLAGS) $(OPT) $(LDFLAGS) -lcasimir cylinder.cpp -o cylinder_bla


clean:
	rm -f casimir_cylinder casimir_tests casimir_logdetD casimir *.o *.so cquadpack/src/*.o tests/*.o $(DEPEND_FILE)
	rm -rf doc
	cd libhodlr; make clean

doc:
	doxygen doxygen.conf
