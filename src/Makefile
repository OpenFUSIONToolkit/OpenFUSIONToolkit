# OpenFUSIONToolkit makefile
# All platform specific customization is done
# in associated make_libs.mk files for each system

include make_cmd.mk
include make_inc.mk

.IGNORE:
.PHONY: oftphysics oftfem oftla oftgrid oftbase programs tests

# Make alias for eclipse
all: programs

# Compile executables
programs: oftlib
	cd bin; $(MAKE)

# Make libraries
oftlib: oftphysics oftfem oftla oftgrid oftbase oftinclude oftext

# Compile external libraries
oftext:
	cd ext_libs; $(MAKE)

# Create includes
oftinclude:
	cd include; $(MAKE)

# Compile base runtime
oftbase:
	cd base; $(MAKE)

# Compile gridding
oftgrid: oftbase
	cd grid; $(MAKE)

# Compile linear algebra
oftla: oftgrid
	cd lin_alg; $(MAKE)

# Compile fem
oftfem: oftla oftgrid
	cd fem; $(MAKE)

# Compile physics
oftphysics: oftfem oftla oftgrid
	cd physics; $(MAKE)

# Compile Tests
test: $(TW_DRIVERS) $(GS_DRIVERS) oftlib
	cd tests; $(MAKE)

# Compile Examples
example: oftlib
	cd examples; $(MAKE)

# Remove built elements
clean:
	cd ext_libs; $(MAKE) clean
	cd base; $(MAKE) clean
	cd grid; $(MAKE) clean
	cd lin_alg; $(MAKE) clean
	cd fem; $(MAKE) clean
	cd physics; $(MAKE) clean
	cd bin; $(MAKE) clean
	cd examples; $(MAKE) clean
	cd tests; $(MAKE) clean

# Print help information
help:
	@echo
	@echo 'Make targets for OpenFUSIONToolkit (OFT):'
	@echo '  all             Build all available driver programs'
	@echo '    taylor          Build Taylor driver programs'
	@echo '    grad_shaf       Build Grad-Shafranov driver programs'
	@echo '    thin_wall       Build Thin-Wall driver programs'
	@echo
	@echo '  oftlib          Build OFT libraries'
	@echo '  test            Compile OFT test suite'
	@echo '  example         Compile OFT examples'
	@echo '  clean           Clean OFT object, module and executable files'
	@echo
	@echo 'New executables can be built using a generic build method by using'
	@echo 'the filename as the target without its extension. For example to'
	@echo 'compile a program contained in the file "example.F90" use the'
	@echo 'following command, which will create an executable named "example".'
	@echo
	@echo '  ~$$ make example '
	@echo

FORCE:
