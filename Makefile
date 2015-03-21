PREFIX=$(HOME)
CC=g++
AR=ar
CCFLAGS= -Wall -std=gnu++0x -Ofast 
LDFLAGS= -lntl -lgmp 

SRCS=$(wildcard *.cc)
OBJS=$(SRCS:.cc=.o)

.PHONY: all clean mrproper

all: IBE

IBE: $(OBJS)
	$(CC) $(CCFLAGS) -o IBE $(OBJS) $(LDFLAGS)

%.o: %.cc params.h
	$(CC) $(CCFLAGS) -c $< 

clean:
	rm -f *.o

mrproper:
	rm -f IBE
