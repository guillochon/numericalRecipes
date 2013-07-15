#Makefile for NR

flashUtilities += nr.o nrutil.o nrtype.o newt.o annewt.o fmin.o fdjac.o eosjac.o \
                  ludcmp.o lubksb.o lnsrch.o zbrac.o zbrent.o newt_wrappers.o

nr.o : nrtype.o
nrutil.o : nrtype.o
newt.o : nr.o nrutil.o fdjac.o
annewt.o : nr.o nrutil.o eosjac.o
fmin.o : nrutil.o
fdjac.o : nrutil.o
eosjac.o : nrutil.o Eos_data.o
ludcmp.o : nrutil.o
lubksb.o : nrutil.o
lnsrch.o : nrutil.o
zbrac.o : nrutil.o
zbrent.o : nrutil.o
newt_wrappers.o : nrtype.o newt.o annewt.o
