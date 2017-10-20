#CPP = -DSELF_START -DSTEPOUT
#CPP = -DUNIX -DVECTOR -DFAST_SQRT
CPP = -DUNIX 
MSP_OPTIONS = -DMSP -DFACOM -DLARGEN -DFOOLISH_OS
SX2_OPTIONS = -DSX2 -DLARGEN -DFOOLISH_OS
#FOPTIONS= -g -u  $(CPP)
FOPTIONS=  -u  $(CPP) -O2
MISCS = Makefile cdel.awk run.csh index.html testin
HEADERS = pcinc.h pcheader.h
SUBSOURCES = intgrt.F timestep.F terror.F intgrt2.F hermite.F intgrt3.F \
		startup.F mkmodel.F output.F
SUBOBJS = intgrt.o timestep.o terror.o intgrt2.o hermite.o intgrt3.o \
		startup.o mkmodel.o output.o

MAINSOURCES = pccode.F testbed.F testbed3.F
F77 = gfortran
pccode : pccode.o $(SUBOBJS)
	$(F77) $(FOPTIONS) pccode.o $(SUBOBJS) -o pccode
pccodeb : pccode.o $(SUBOBJS)
	$(F77) $(FOPTIONS) pccode.o $(SUBOBJS) -o pccodeb
testbed3 : testbed3.F $(SUBOBJS)
	$(F77) $(FOPTIONS) testbed3.F $(SUBOBJS) -o testbed3
testbed3b : testbed3.F intgrtb.F $(SUBOBJS)
	$(F77) $(FOPTIONS) -DROUND_OFF_COMP -DFACTOR2DT \
	 testbed3.F intgrtb.F $(SUBOBJS) -o testbed3b
testbed : testbed.o $(SUBOBJS)
	$(F77) $(FOPTIONS) testbed.o $(SUBOBJS) -o testbed


pccode.tgz : $(MAINSOURCES)  $(SUBSOURCES) $(HEADERS) $(MISCS)
	tar cvzf $@   $(MAINSOURCES)  $(SUBSOURCES) $(HEADERS) $(MISCS)

pccode.shar : $(MAINSOURCES)  $(SUBSOURCES) $(HEADERS) $(MISCS)
	bundle $(MAINSOURCES)  $(SUBSOURCES) $(HEADERS) $(MISCS)\
	 >pccode.shar

msp : pccode.MSP
	ls -l pccode.MSP
testbed.MSP : testbed.F $(SUBSOURCES) $(HEADERS)
	cat testbed.F $(SUBSOURCES)  | /lib/cpp -P $(MSP_OPTIONS) >temp
	dd conv=ucase <temp | awk -f cdel.awk >testbed.MSP
	rm temp
pccode.SX2 : pccode.F $(SUBSOURCES) $(HEADERS)
	cat pccode.F $(SUBSOURCES)  | /lib/cpp -P $(SX2_OPTIONS) >temp
	dd conv=ucase <temp | awk -f cdel.awk >pccode.SX2
	rm temp

hermite :hermite.F $(HEADERS) $(SUBOBJS)
	$(F77)  -DTEST $(CPP)  hermite.F $(SUBOBJS) -o hermite

clean :
	rm *.o *.for

intgrt.o :intgrt.F pcheader.h
	$(F77) -c $(FOPTIONS)  intgrt.F
intgrt2.o :intgrt2.F pcheader.h pcinc.h
	$(F77) -c $(FOPTIONS)  intgrt2.F
intgrt3.o :intgrt3.F pcheader.h
	$(F77) -c $(FOPTIONS)  intgrt3.F
timestep.o :timestep.F pcheader.h
	$(F77) -c $(FOPTIONS)  timestep.F
testbed.o :testbed.F pcheader.h
	$(F77) -c $(FOPTIONS)  testbed.F
terror.o :terror.F pcheader.h
	$(F77) -c $(FOPTIONS)  terror.F
hermite.o :hermite.F pcheader.h
	$(F77) -c $(FOPTIONS)  hermite.F
pccode.o :pccode.F pcheader.h pcinc.h
	$(F77) -c $(FOPTIONS)  pccode.F
startup.o :startup.F pcheader.h 
	$(F77) -c $(FOPTIONS)  startup.F
mkmodel.o :mkmodel.F pcheader.h 
	$(F77) -c $(FOPTIONS)  mkmodel.F
output.o :output.F pcheader.h pcinc.h
	$(F77) -c $(FOPTIONS)  output.F

#---------------------------------------------------------------------
# local part
INSTALLDIR=/usr2/makino/WWW/softwares/pccode
INSTALLFILES= index.html pccode.tgz
install: $(INSTALLFILES)
	rsync -avprog $(INSTALLFILES) $(INSTALLDIR)
