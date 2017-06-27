# Choose compiler
CC = mpicc
# setenv GSL_RNG_TYPE r250 -- to set environment variable
# Optimization flags: cc is generic, gcc and icc are for PIII
CFLAGS_gcc = -g -Wall -DHAVE_INLINE #Usually gcc is needed for debugging
CFLAGS_icc = -O3 -xP  -DHAVE_INLINE
CFLAGS_mpicc = -O3 -xP -DHAVE_INLINE

# Set linker flags: i.e. mpe libraries or options like -static

#LFLAGS_gcc = -static
#LFLAGS_icc = -static

LFLAGS_gcc = 
LFLAGS_icc = 
LFLAGS_mpiicc = 

# Compilation and linking flags 
CFLAGS = $(CFLAGS_$(CC)) -c
LFLAGS = $(LFLAGS_$(CC)) -lm -lgsl -lgslcblas


C_SRC     = main.c init.c polymerize.c misc.c interface.c rod.c forces.c propagator.c rod_misc.c eulerstep.c
OBJ_flex  = main.o init.o polymerize.o misc.o interface.o rod.o forces.o propagator.o rod_misc.o eulerstep.o

.c.o:
	$(CC)  $(CFLAGS) $*.c

flex: $(OBJ_flex)
	$(CC)  $(LFLAGS) -o mt $(OBJ_flex)

clean:
	rm -f *.o mt
