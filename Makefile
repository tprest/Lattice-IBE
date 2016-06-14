PREFIX=$(HOME)
CC=g++
AR=ar

# Case 1: These are the standard compilation flags CCFLAGS and linker flags LDFLAGS.
CCFLAGS= -Wall -std=gnu++0x -Ofast 
LDFLAGS= -lntl -lgmp 

# Case 2: If NTL is installed in a specific location, say /path/to/ntl, you must specify it by using the CCFLAGS and LDFLAGS below instead.
# CCFLAGS= -Wall -I/path/to/ntl/include/ -std=gnu++0x -Ofast 
# LDFLAGS= -L/path/to/ntl/lib/ -lntl -lgmp 

# Case 3: In some cases, Unix doesn't find NTL even if it is installed in the standard location /usr/local. Then, uncomment the following lines.
# CCFLAGS= -Wall -I/usr/local/include/ -std=gnu++0x -Ofast 
# LDFLAGS= -L/usr/local/lib/ -lntl -lgmp 

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
