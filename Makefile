#
# Makefile for SKID (Spline Kernel Interpolative Denmax).
#
# on an SGI R4400 add -mips2 to CFLAGS
# on an SGI R8000 (TFP) add -mips4 -O3 to CFLAGS
#

#
# Defines for pgcc (Portland Group Compiler)
# Warning: generates incorrect code for kd.c with optimize, in tree build!
# Compile kd.c with no optimize and things are okay.
#
#CC = pgcc
#CFLAGS = -fast
#CFLAGSKD = -O0

#
# Defines for gcc, no optimization problems here.
#
CFLAGS = -O3 -funroll-loops
CFLAGSKD = $(CFLAGS)


LIBS	=   -lm

default:	skid totipnat
	@echo To try the demo, type demo!

clean:
	rm -f *.o

skid: main.o kd.o smooth1.o grav.o cosmo.o romberg.o
	$(CC) $(CFLAGS) -o skid main.o kd.o smooth1.o grav.o cosmo.o romberg.o $(LIBS)

main.o: main.c kd.h smooth1.h cosmo.h

kd.o: kd.c kd.h tipsydefs.h cosmo.h
	$(CC) $(CFLAGSKD) -c kd.c

smooth1.o: smooth1.c kd.h smooth1.h cosmo.h

grav.o: grav.c grav.h kd.h cosmo.h

#
#	May need to specify -lrpc on some systems
#
totipnat: totipnat.c tipsydefs.h
	$(CC) $(CFLAGS) -o totipnat totipnat.c $(LIBS)


