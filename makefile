#----------------------------------------------------------------------
# MIK: Building cascade of executable files --------------------------
#----------------------------------------------------------------------
# Updates: 11/26/2018- added dlsode.f to OBJS, S.M. Levine
#----------------------------------------------------------------------

OBJS =  local.f95 dlsode.f

DOFL  = ifort    -O3 -fp-model precise -fp-model source
DOFL2 = gfortran -O3 -fbounds-check
DOFL3 = gfortran -O0 -fbounds-check -mcmodel=medium -g -Wall -fbacktrace

all: mik clean

mik: $(OBJS)
	$(DOFL2) $(OBJS) -o mik.exe

clean:
	rm -rf *.mod

### end of makefile
