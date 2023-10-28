# Open FUSION Toolkit (OFT) testing makefile

.PHONY: tests

#---Compile test cases
tests:
	cd grid; $(MAKE)
	cd lin_alg; $(MAKE)
	cd fem; $(MAKE)
	cd physics; $(MAKE)

clean:
	cd grid; $(MAKE) clean
	cd lin_alg; $(MAKE) clean
	cd fem; $(MAKE) clean
	cd physics; $(MAKE) clean