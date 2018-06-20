CFLAGS = -std=c99 -g -Wall -Wextra -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wno-unused-parameter -fstrict-aliasing
OPT = -O3

all:
	$(CC) $(OPT) $(CFLAGS) -fPIC -shared -o integrand.so integrand.c -lm

clean:
	rm -f integrand.so
