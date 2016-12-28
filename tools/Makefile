BSHELL = /bin/sh

CC = $(CXX)
CFLAGS =  -Wall -Wextra -O2 -std=c++14 

all: BWT

BWT:
	$(CC) -I ../parallel-divsufsort/include $(CFLAGS) -o BWT bwt.cpp -L ../parallel-divsufsort/lib -ldivsufsort -lcilkrts -fcilkplus -llibprange
clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f *.o BWT
