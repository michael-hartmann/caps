.PHONY: test

OPT = -O3
CXXFLAGS = -DUSE_DOUBLE -g -w -L .. -lhodlr -I ../include/
LDFLAGS = -L .. -lhodlr

# undefined behaviour sanitizer supported by clang++
#CXXFLAGS += -fsanitize=undefined,unsigned-integer-overflow -fno-omit-frame-pointer

all:
	$(CXX) $(CXXFLAGS) ${OPT} test.c ${LDFLAGS} -o test

test:
	./test

clean:
	rm -f test
