#----------------------------------------------------------------------
# MIK: Building cascade of executable files --------------------------
#----------------------------------------------------------------------
# 12/04/2018 by S.M. Levine
#----------------------------------------------------------------------

OBJS =  local.f95 opkdmain.f opkda1.f opkda2.f

DOFL  = ifort    -O3 -fp-model precise -fp-model source
DOFL2 = gfortran -O3 -fbounds-check
DOFL3 = gfortran -O0 -fbounds-check -mcmodel=medium -g -Wall -fbacktrace

all: mik clean

mik: $(OBJS)
	$(DOFL2) $(OBJS) -o mik.exe

clean:
	rm -rf *.mod

### end of makefile
