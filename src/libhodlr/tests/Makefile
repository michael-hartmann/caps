.PHONY: test

OPT = -O3
CXXFLAGS = -DUSE_DOUBLE -g -w -L .. -lhodlr -I ../include/

# undefined behaviour sanitizer supported by clang++
#CXXFLAGS += -fsanitize=undefined,unsigned-integer-overflow -fno-omit-frame-pointer

all:
	$(CXX) $(CXXFLAGS) ${OPT} test.c -o test

test:
	./test

clean:
	rm -f test
