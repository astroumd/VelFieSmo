#

MKGALCUBE = $(shell readlink -f ../scripts/mkgalcube)

.PHONY: model0 test1

#        model0 is special, we always make it
#        it optionally runs Bbarolo, if it exists
model0:
	rm -rf model0
	$(MKGALCUBE)
	-(cd model0; Bbarolo -p ../../data/model0_galmod.par)


# test1:  see if bsc=begeman can work (so far, it cannot, or too much numerical pixel noise)

# parameters for test1
INC   = 60
M     = 4
R0    = 1
RE    = 1000
BEAM1 = 1
BEAM2 = 2
NBODY = 1000000

test1:
	rm -rf test1
	mkdir test1
	(cd test1; $(MKGALCUBE) run=m0 beam=0        inc=$(INC) m=$(M) r0=$(R0) re=$(RE) nbody=$(NBODY))
	(cd test1; $(MKGALCUBE) run=m1 beam=$(BEAM1) inc=$(INC) m=$(M) r0=$(R0) re=$(RE) nbody=$(NBODY))
	(cd test1; $(MKGALCUBE) run=m2 beam=$(BEAM2) inc=$(INC) m=$(M) r0=$(R0) re=$(RE) nbody=$(NBODY))

#   1:  V
#   2:  V - v
#   3:  V + delta         (should be close to v)
#   4:  V + delta - v     (should be close to 0)
test1a:
	-testbsc test1/m0/m0.33d test1/m0/m0.33v test1/m1/m1.33d test1/m1/m1.33v test1/m1/bsc1 diff=f beeg=f
	-testbsc test1/m0/m0.33d test1/m0/m0.33v test1/m1/m1.33d test1/m1/m1.33v test1/m1/bsc2 diff=t beeg=f
	-testbsc test1/m0/m0.33d test1/m0/m0.33v test1/m1/m1.33d test1/m1/m1.33v test1/m1/bsc3 diff=f beeg=t 
	-testbsc test1/m0/m0.33d test1/m0/m0.33v test1/m1/m1.33d test1/m1/m1.33v test1/m1/bsc4 diff=t beeg=t
	nds9 -1 test1/m1/bsc1
	nds9 -2 test1/m1/bsc2
	nds9 -3 test1/m1/bsc3
	nds9 -4 test1/m1/bsc4

test1b:
	-testbsc test1/m0/m0.33d test1/m0/m0.33v test1/m2/m2.33d test1/m2/m2.33v test1/m2/bsc1 diff=f beeg=f
	-testbsc test1/m0/m0.33d test1/m0/m0.33v test1/m2/m2.33d test1/m2/m2.33v test1/m2/bsc2 diff=t beeg=f
	-testbsc test1/m0/m0.33d test1/m0/m0.33v test1/m2/m2.33d test1/m2/m2.33v test1/m2/bsc3 diff=f beeg=t 
	-testbsc test1/m0/m0.33d test1/m0/m0.33v test1/m2/m2.33d test1/m2/m2.33v test1/m2/bsc4 diff=t beeg=t 
	nds9 -1 test1/m2/bsc1
	nds9 -2 test1/m2/bsc2
	nds9 -3 test1/m2/bsc3
	nds9 -4 test1/m2/bsc4

ngc6503:
	$(MKGALCUBE) run=ngc6503 step=10 pa=-60 v0=115 r0=0.15 m=2 beam=0.23 vsys=28 sigmav=8

fig1:
	rm -rf fig1; $(MKGALCUBE) run=fig1 step=10 beam=2 m=2


test2:
	for m in 1 2 3 4; do \
	  for beam in 0 0.5 1 2 ; do \
	     echo M=$$m BEAM=$$beam ;\
	     $(MKGALCUBE) run=run_$$m:$$beam step=10 m=$$m beam=$$beam > run.log 2>&1 ;\
	     mv run.log run_$$m:$$beam ;\
          done ;\
	done
