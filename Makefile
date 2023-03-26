# Arquivo makefile para rodar mais facilmente o que necessito
FC = gfortran
FFLAGS = -march=native -mtune=native -O3
FDEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,underflow,denormal
#IDIR = -I/usr/local/include/fgsl/
#LIBS = -lgsl -lfgsl -lopenblas -lm

all: difusao1dpermanente.exe

difusao1dpermanente.exe: difusao1dpermanente.f90 parameters.o
	$(FC) $(FFLAGS) -o $@ $(IDIR) $^ $(LIBS)

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $(IDIR) $<

debug: difusao1dpermanente.f90 parameters.f90
	$(FC) $(FDEBUGFLAGS) -o debug_difusao1dpermanente.exe $(IDIR) $^ $(LIBS)

clean:
	rm -f *.o *.mod *.exe *.txt
