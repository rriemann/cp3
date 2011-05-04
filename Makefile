CXXFLAGS=`root-config --cflags` -Wall -g
# -fopenmp
LDLIBS=`root-config --libs`
#  -lgomp

all:	main

main:	geom_pbc.o

clean:
	rm -f main main.o geom_pbc.o

.PHONY:	clean