CFLAGS = -Wall -I ../src/include/
LDFLAGS = -L ../src/ -L ../src/libhodlr/ -lm -llapack -lblas -lgfortran -lcasimir -lhodlr
OPT = -O3

AuAl:
	$(CC) $(OPT) $(CFLAGS) AuAl.c -c
	$(CC) AuAl.o $(LDFLAGS) -o AuAl

clean:
	rm -f *o AuAl
