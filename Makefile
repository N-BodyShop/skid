#
# Makefile for SKID (Spline Kernel Interpolative Denmax).
#
# on an SGI R4400 add -mips2 to CFLAGS
# on an SGI R8000 (TFP) add -mips4 -O3 to CFLAGS
#
CC=cc
CFLAGS	=	-O2
LIBS	=   -lm

default:	skid totipnat
	@echo To try the demo, type demo!

clean:
	rm -f *.o

skid: main.o kd.o smooth1.o grav.o
	$(CC) $(CFLAGS) -o skid main.o kd.o smooth1.o grav.o $(LIBS)

main.o: main.c kd.h smooth1.h
	$(CC) $(CFLAGS) -c main.c

kd.o: kd.c kd.h tipsydefs.h
	$(CC) $(CFLAGS) -c kd.c

smooth1.o: smooth1.c kd.h smooth1.h
	$(CC) $(CFLAGS) -c smooth1.c

grav.o: grav.c grav.h kd.h
	$(CC) $(CFLAGS) -c grav.c

#
#	May need to specify -lrpc on some systems
#
totipnat: totipnat.c tipsydefs.h
	$(CC) $(CFLAGS) -o totipnat totipnat.c $(LIBS)


