#

MKGALCUBE = $(shell readlink -f ../scripts/mkgalcube)


#        model0 is special, we always make it
model0: .model0

.model0:
	rm -rf model0
	$(MKGALCUBE)


# test1:  see if bsc=begeman can work (so far, it cannot, or too much numerical pixel noise)

# parameters for test1
INC   = 60
M     = 4
R0    = 1
BEAM1 = 1
BEAM2 = 2

test1:
	mkdir test1
	(cd test1; $(MKGALCUBE) run=m0 beam=0        inc=$(INC) m=$(M) r0=$(R0))
	(cd test1; $(MKGALCUBE) run=m1 beam=$(BEAM1) inc=$(INC) m=$(M) r0=$(R0))
	(cd test1; $(MKGALCUBE) run=m2 beam=$(BEAM2) inc=$(INC) m=$(M) r0=$(R0))

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
