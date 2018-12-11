CXXFLAGS = -Wall -Wextra -I.. -L.. -L../libhodlr
LDFLAGS = -lcasimir -lhodlr -lgfortran -lblas -llapack
OPT = -O3

all:
	$(CXX) $(CXXFLAGS) $(OPT) $(LDFLAGS) cylinder.cpp -o cylinder
