CXXFLAGS=`root-config --cflags` -Wall -g
# -fopenmp
LDLIBS=`root-config --libs`
#  -lgomp

all:	main

clean:
	rm -f main main.o

.PHONY:	clean