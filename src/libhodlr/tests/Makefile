OPT = -O3

all:
	clang++ -fsanitize=undefined -fno-omit-frame-pointer -w ${OPT} test.c -I ../include/ ../src/hodlr.cpp -I ../src/header/ -o test

clean:
	rm -f test
