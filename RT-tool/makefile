#--------RT tool makefile-------
#--------Thomas G. Bisbas-------

FC = gfortran
FFLAGS = -g -fbacktrace -fcheck=all
DUST = 0

ifeq ($(DUST),1)
   OPT += -DDUST
endif

SRC = modules.F90 readcoolants.F90 read_input.F90 readparams.F90 \
      readpdrfiles.F90 initializations.F90 writeoutputs.F90 excitation.F90 \
      readspopfile.F90 opticaldepth.F90 radiativetransfer.F90 RTtool.F90

EXEC = RTtool

all: $(EXEC)
	mv $(EXEC) ../

$(EXEC): $(SRC)
	$(FC) $(FFLAGS) $(OPT) -o $(EXEC) $(SRC)

clean:
	rm -f $(EXEC) *.o *.mod

