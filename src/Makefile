# Makefile for UMAT Builder modular library
#
# Usage:
#   make              Build all object files for linking with ABAQUS
#   make clean        Remove build artifacts
#
# To generate a material law:
#   python ../generate.py config.json
#   python ../generate.py --example neo_hooke

FC       = gfortran
FFLAGS   = -O2 -ffree-form -fPIC
BUILDDIR = build

# Module compilation order (dependency chain)
MODULES = \
	$(BUILDDIR)/mod_constants.o \
	$(BUILDDIR)/mod_tensor.o \
	$(BUILDDIR)/mod_kinematics.o \
	$(BUILDDIR)/mod_continuum.o \
	$(BUILDDIR)/mod_hyperelastic.o \
	$(BUILDDIR)/mod_anisotropic.o \
	$(BUILDDIR)/mod_icosahedron.o \
	$(BUILDDIR)/mod_network.o \
	$(BUILDDIR)/mod_damage.o \
	$(BUILDDIR)/mod_viscosity.o

UMAT_OBJ = $(BUILDDIR)/umat_builder.o $(BUILDDIR)/uexternaldb.o

.PHONY: all clean

all: $(BUILDDIR) $(MODULES) $(UMAT_OBJ)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Module dependencies (order matters)
$(BUILDDIR)/mod_constants.o: mod_constants.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILDDIR)

$(BUILDDIR)/mod_tensor.o: mod_tensor.f90 $(BUILDDIR)/mod_constants.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_kinematics.o: mod_kinematics.f90 $(BUILDDIR)/mod_tensor.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_continuum.o: mod_continuum.f90 $(BUILDDIR)/mod_tensor.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_hyperelastic.o: mod_hyperelastic.f90 $(BUILDDIR)/mod_tensor.o $(BUILDDIR)/mod_kinematics.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_anisotropic.o: mod_anisotropic.f90 $(BUILDDIR)/mod_constants.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_icosahedron.o: mod_icosahedron.f90 $(BUILDDIR)/mod_constants.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_network.o: mod_network.f90 $(BUILDDIR)/mod_constants.o $(BUILDDIR)/mod_icosahedron.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_damage.o: mod_damage.f90 $(BUILDDIR)/mod_constants.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_viscosity.o: mod_viscosity.f90 $(BUILDDIR)/mod_constants.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/umat_builder.o: umat_builder.f90 $(MODULES)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/uexternaldb.o: uexternaldb.f90
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

clean:
	rm -rf $(BUILDDIR)
