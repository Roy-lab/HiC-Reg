# the compiler: gcc for C program, define as g++ for C++
CC = g++

# compiler flags:
#-gadds debugging information to the executable file
#-Wall turns on most, but not all, compiler warnings
CFLAGS= -g  -L/usr/lib/ -lm #-Wall 
GSLFLAGS=-lgsl -lgslcblas 

all: aggregate_Signal

aggregate_Signal: Framework.C
	$(CC) $(CFLAGS) -o aggregateSignal *.C $(GSLFLAGS) 

clean:
	$(RM) aggregate_Signal
