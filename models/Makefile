#

MKGALCUBE = $(shell readlink -f ../scripts/mkgalcube)

.PHONY: model0 test1

#        model0 is special, we always make it
#        it optionally runs Bbarolo, if it exists
model0:
	rm -rf model0
	$(MKGALCUBE)
	-(cd model0; Bbarolo -p ../../data/model0_galmod.par > b1.log)
	-(cd model0; Bbarolo -p ../../data/model0_2dfit.par  > b2.log)
	-(cd model0; Bbarolo -p ../../data/model0_3dfit.par  > b3.log)


# test1:  see if bsc=begeman can work (so far, it cannot, or too much numerical pixel noise)

# parameters for test1 (yes, upper case)
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

# attempt to get a pretty nice model for N6503, exactly how it's on the sky etc.etc.
ngc6503:
	$(MKGALCUBE) run=ngc6503 step=10 pa=-60 v0=115 r0=0.15 m=2 beam=0.23 vsys=28 sigmav=8

# rms=0.0007  robust=0.00057
ngc6503fit:
	fitsccd ../data/ngc6503vrad.fits ngc6503.ccd error=1
	ccdmom ngc6503.ccd ngc6503den.ccd mom=0 clip=0.001 error=1
	ccdmom ngc6503.ccd ngc6503vel.ccd mom=1 rngmsk=t clip=0.003 error=1
	rotcurshape ngc6503vel.ccd 0,600 -60 74 28 162,122 rotcur1=core,115,50,2,1,1,1  fix=vsys,xpos,ypos,pa,inc
# exp  -> 117.2 (0.1)  27.0 (0.3)               rms=4.0
# core -> 117.2 (0.1)  48.8 (0.6)  4.1 (0.2)    rms=3.9
# density fit shows exp scalelength about 358 \pm 5, i.e. RE/R1=7.3
# B/R1 = 0.287 -> m0 ~ 4.5 ->  R1/R0 ~ 1.1 -> R0 ~ 44
# (44,4.5) --[14]-> (49,4.1)



fig1:
	rm -rf fig1; $(MKGALCUBE) run=fig1 step=10 beam=2 m=2


# a series of runs for different 'm' and 'beam' for fixed re
test2:
	rm -rf run_*
	for m in 1 2 3 4; do \
	  for beam in 0 0.5 1 2 ; do \
	     echo M=$$m BEAM=$$beam ;\
	     $(MKGALCUBE) run=run_$$m:$$beam m=$$m beam=$$beam > run.log 2>&1 ;\
	     mv run.log run_$$m:$$beam ;\
          done ;\
	done

# test if scaling produces the same (scaled) results [now it does]
test3:
	rm -rf model0a; $(MKGALCUBE) run=model0a r0=2 re=2 beam=2 v0=200 vrange=300 step=20 rmax=30 range=32 sigmav=10
	rm -rf model0b; $(MKGALCUBE) run=model0b r0=1 re=1 beam=1 v0=100 vrange=150 step=10 rmax=15 range=16 sigmav=5

# Example to get the large table for plotting
# grep FITTING run4_*/run.log  | grep -v beam | awk -F: '{print $5}' | awk '{if (NF==12) print $0}' | grep -v "100 120"

#       a disk at inc=75 (the default in mkgalcube actually).
run4:
	rm -rf run4_*
	for re in 0.5 1 2 3 4 5 6 10 100 1000; do \
	 for m in 1 1.5 2 2.5 3 3.5 4 4.5 5; do \
	  for beam in 0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.5 3 3.5 4 5 6 8 ; do \
	     echo M=$$m RE=$$re BEAM=$$beam ;\
	     $(MKGALCUBE) run=run4_$$m:$$re:$$beam m=$$m re=$$re beam=$$beam inc=75 > run.log 2>&1 ;\
	     mv run.log run4_$$m:$$re:$$beam ;\
          done ;\
	 done ;\
	done

# M51, from the CASAGUIDES, but with different inclination and rotation curve 
SEED = 123
BEAM = 1

m51:
	rm -rf m51a m51b
	$(MKGALCUBE) run=m51a cmap=../data/M51ha.ccd inc=45 beam=$(BEAM) seed=$(SEED) nbody=$(NBODY)
	$(MKGALCUBE) run=m51b                        inc=45 beam=$(BEAM) seed=$(SEED) nbody=$(NBODY)
	ccdmath m51a/m51a.33v,m51b/m51b.33v m51b/diff.vel %2-%1
	nds9 -1 m51a/m51a.33d
	nds9 -2 m51a/m51a.33v
	nds9 -3 m51b/m51b.33d
	nds9 -4 m51b/m51b.33v
	# nds9 -5 m51b/diff.vel
